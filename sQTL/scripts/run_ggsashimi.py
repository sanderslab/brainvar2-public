#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Sashimi plot pipeline using ggsashimi (Docker).
Usage: python run_ggsashimi.py <TSV>

The TSV filename must contain the sQTL region in the format:
    ...rsXXXXXX__CHROM_START_END...

The TSV must have columns: sample_id, genotype, stage
"""

import pandas as pd
import os
import subprocess
import re
import sys

# =========================
# CONFIG
# =========================

SRC_BAM_DIR   = "/mnt/s3/FXN/SoowheeKim/AlignedBam1/"
LOCAL_BAM_DIR = "/mnt/s3/eqtl/Bulk_RNA/sqtl_wasp/"
GTF           = "/home/ubuntu/gencode.v41.annotation.gtf"
STRAND        = "-"   # strand of the sQTL gene; change to "+" if needed

MIN_SAMPLES   = 2

# ggsashimi Docker image
GGSASHIMI_PY = "./ggsashimi.py"   # path to ggsashimi.py

# Palette: genotype 0/1/2 -> colour (used in the BAM manifest)
GENOTYPE_COLORS = {
    "0": "#2166ac",   # homozygous ref  -> blue
    "1": "#f4a582",   # heterozygous    -> orange
    "2": "#d6604d",   # homozygous alt  -> red
}

# ggsashimi plot options (passed as CLI flags)
GGSASHIMI_OPTS = [
    "--min-coverage", "1",
    "--alpha",        "0.5",
    "--base-size",    "12",
    "--ann-height",   "1.5",
    "--height",       "3",
    "--width",        "10",
]

# ── sanity checks ──────────────────────────────────────────────────────────────

if not os.path.exists(LOCAL_BAM_DIR):
    os.makedirs(LOCAL_BAM_DIR)
    print("[CONFIG] Created LOCAL_BAM_DIR: %s" % LOCAL_BAM_DIR)

for lbl, val in [
    ("SRC_BAM_DIR",   SRC_BAM_DIR),
    ("LOCAL_BAM_DIR", LOCAL_BAM_DIR),
    ("GTF",           GTF),
    ("STRAND",        STRAND),
    ("MIN_SAMPLES",   MIN_SAMPLES),
]:
    print("[CONFIG] %-14s: %s" % (lbl, val))
print("[CONFIG] GTF exists   : %s" % os.path.exists(GTF))

# =========================
# HELPERS
# =========================

def parse_region(tsv):
    """Extract rsid, chrom, start, end from the TSV filename."""
    print("\n[parse_region] %s" % tsv)
    base = os.path.basename(tsv)
    m = re.search(r"(rs\d+)__(\d+)_(\d+)_(\d+)", base)
    if not m:
        raise ValueError("[parse_region] Filename format incorrect: %s" % base)
    rsid, chrom, start, end = m.groups()
    print("[parse_region] rsid=%s | chrom=%s | start=%s | end=%s" % (rsid, chrom, start, end))
    return rsid, chrom, start, end


def get_gene_coords(chrom, start, end):
    """
    Find the gene in GFF3 that overlaps the leafcutter sQTL intron cluster.
    Filters by STRAND. Returns (gene_start, gene_end) as strings.
    """
    print("[get_gene_coords] Overlapping chr%s:%s-%s strand=%s" % (chrom, start, end, STRAND))
    chrom_str   = "chr%s" % chrom
    query_start = int(start)
    query_end   = int(end)
    best = None

    with open(GTF) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            if parts[0] != chrom_str or parts[2] != "gene" or parts[6] != STRAND:
                continue
            g_start, g_end = int(parts[3]), int(parts[4])
            if g_start <= query_end and g_end >= query_start:
                print("[get_gene_coords] Found: %s:%d-%d:%s" % (chrom_str, g_start, g_end, STRAND))
                best = (str(g_start), str(g_end))
                break

    if best is None:
        print("[get_gene_coords] WARNING: no gene found, falling back to sQTL coords")
        best = (start, str(int(end) + 1))

    print("[get_gene_coords] Using: chr%s:%s-%s" % (chrom, best[0], best[1]))
    return best


def get_bam(sample):
    """Copy BAM from SRC to LOCAL, always re-index to ensure fresh .bai. Returns local BAM path."""
    bam_name = "%s.STARAligned.sortedByCoord.WASPfiltered.out.bam" % sample
    src = os.path.join(SRC_BAM_DIR, bam_name)
    dst = os.path.join(LOCAL_BAM_DIR, bam_name)

    if not os.path.exists(dst):
        print("  [get_bam] Copying %s ..." % bam_name)
        subprocess.check_call(["cp", src, dst])
    else:
        print("  [get_bam] BAM already local: %s" % bam_name)

    # Validate index with samtools quickcheck; only re-index if missing or stale.
    bai = dst + ".bai"
    needs_index = False
    if not os.path.exists(bai):
        print("  [get_bam] No index found, will index.")
        needs_index = True
    else:
        ret = subprocess.call(["samtools", "quickcheck", dst], stderr=open(os.devnull, "w"))
        if ret != 0:
            print("  [get_bam] quickcheck failed, will re-index.")
            needs_index = True
        else:
            print("  [get_bam] Index OK, skipping.")
    if needs_index:
        subprocess.check_call(["samtools", "index", "-@", "4", dst])
        print("  [get_bam] Index written: %s.bai" % dst)

    return dst


def prepare_groups(df_stage):
    """
    Build genotype groups {genotype_str: [bam_path, ...]}.
    Only includes groups with >= MIN_SAMPLES samples.
    """
    print("\n[prepare_groups] Genotype counts: %s" % df_stage["genotype"].value_counts().to_dict())
    groups = {}
    for g in [0, 1, 2]:
        samples = df_stage[df_stage["genotype"] == g]["sample_id"].tolist()
        print("[prepare_groups] Genotype %d -> %d sample(s)" % (g, len(samples)))
        if len(samples) >= MIN_SAMPLES:
            groups[str(g)] = [get_bam(s) for s in samples]
        else:
            print("[prepare_groups] Genotype %d SKIPPED (< %d samples)" % (g, MIN_SAMPLES))
    return groups


def write_bam_manifest(groups, outdir):
    """
    Write a ggsashimi BAM manifest TSV:
        SampleName <TAB> /path/to/bam <TAB> colour
    Samples are named GT{g}_{i} (e.g. GT0_1, GT0_2, GT1_1 ...).
    Returns the manifest file path.
    """
    manifest_path = os.path.join(outdir, "bam_manifest.tsv")
    with open(manifest_path, "w") as fh:
        for g, bams in sorted(groups.items()):
            color = GENOTYPE_COLORS.get(g, "#999999")
            for i, bam in enumerate(bams, 1):
                sample_label = "GT%s_%d" % (g, i)
                fh.write("%s\t%s\t%s\n" % (sample_label, bam, color))
    print("[write_bam_manifest] Written: %s" % manifest_path)
    return manifest_path


def run_ggsashimi(groups, chrom, gene_start, gene_end, outdir, label):
    """
    Run ggsashimi.py directly (no Docker).
    """
    print("\n[run_ggsashimi] label=%s | chr%s:%s-%s" % (label, chrom, gene_start, gene_end))

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    manifest = write_bam_manifest(groups, outdir)

    # ggsashimi coordinate format: chrN:start-end (1-based, closed)
    coord = "chr%s:%s-%s" % (chrom, gene_start, gene_end)

    abs_manifest = os.path.abspath(manifest)
    abs_gtf      = os.path.abspath(GTF)
    out_prefix   = os.path.join(os.path.abspath(outdir), "sashimi_%s" % label)

    cmd = (
        [
            "python", GGSASHIMI_PY,
            "-b", abs_manifest,
            "-c", coord,
            "-g", abs_gtf,
            "-o", out_prefix,
            "-S", STRAND,
        ]
        + GGSASHIMI_OPTS
    )

    print("[run_ggsashimi] Command:\n  %s" % " ".join(cmd))

    try:
        subprocess.check_call(cmd)
        print("[run_ggsashimi] SUCCESS -> %s" % os.path.abspath(outdir))
    except subprocess.CalledProcessError as e:
        print("[run_ggsashimi] ERROR: return code %d" % e.returncode)
        print("[run_ggsashimi] Failed command:\n  %s" % " ".join(cmd))


# =========================
# MAIN
# =========================

def main():
    print("\n[main] Starting ggsashimi pipeline ...")

    if len(sys.argv) != 2:
        print("[main] Usage: python run_ggsashimi.py <TSV>")
        sys.exit(1)

    tsv = sys.argv[1]
    if not os.path.exists(tsv):
        print("[main] ERROR: TSV not found: %s" % tsv)
        sys.exit(1)

    rsid, chrom, start, end = parse_region(tsv)

    df = pd.read_csv(tsv, sep="\t")
    print("[main] TSV: %d rows | stages: %s" % (len(df), df["stage"].unique().tolist()))

    for stage in ["Prenatal", "Postnatal"]:
        print("\n[main] ========== Stage: %s ==========" % stage)
        df_stage = df[df["stage"] == stage]

        if df_stage.empty:
            print("[main] No samples, skipping.")
            continue

        groups = prepare_groups(df_stage)

        if len(groups) < 2:
            print("[main] Only %d genotype group(s) - need >=2, skipping %s." % (len(groups), stage))
            continue

        gene_start, gene_end = get_gene_coords(chrom, start, end)

        # Option A: one combined plot for the stage (all genotypes together)
        outdir = "sashimi_%s_%s_all_genotypes" % (rsid, stage)
        label  = "%s_%s_GT%s" % (rsid, stage, "_".join(sorted(groups.keys())))
        run_ggsashimi(groups, chrom, gene_start, gene_end, outdir, label)

        # Option B: pairwise plots (0v1, 0v2, 1v2)
        pairs = [("0", "1"), ("0", "2"), ("1", "2")]
        for g1, g2 in pairs:
            if g1 not in groups or g2 not in groups:
                print("[main] Skipping pair %s vs %s (group missing)" % (g1, g2))
                continue
            pair_groups = {g1: groups[g1], g2: groups[g2]}
            outdir_pair = "sashimi_%s_%s_%s_vs_%s" % (rsid, stage, g1, g2)
            label_pair  = "%s_%s_GT%s_vs_GT%s" % (rsid, stage, g1, g2)
            run_ggsashimi(pair_groups, chrom, gene_start, gene_end, outdir_pair, label_pair)

    print("\n[main] Done!")


if __name__ == "__main__":
    main()
