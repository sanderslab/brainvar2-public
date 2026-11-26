"""Standalone NMF analysis matching the major-cell-type RNA/ATAC Methods.

Run: python nmf_major_celltypes_only.py --workbook '/path/to/Supplementary Tables.xlsx'
Install: pip install numpy pandas scipy scikit-learn matplotlib seaborn requests openpyxl

Only Table S5 is used: no subtype panels, no LDSC.
Final k=5 after exploration.
"""
import argparse
import json
import re
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import requests
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram
from sklearn.decomposition import NMF
from sklearn.metrics.pairwise import cosine_similarity

WORKBOOK_PATH = 'Supplementary Tables.xlsx'  # overridden by --workbook
SHEET_S5 = 'Table S5'
S5_TRAIT_COL, S5_CELLTYPE_COL, S5_STAGE_COL = 'Trait', 'Celltype', 'Stage'
S5_NES_COL, S5_PADJ_COL, S5_MODALITY_COL = 'NES', 'padj', 'Modality'
S5_LEADINGEDGE_COL = 'LeadingEdge'
RANDOM_STATE = 0
NMF_K_FIXED = 5
NMF_STABILITY_RUNS = 100
GO_API_URL = 'https://biit.cs.ut.ee/gprofiler/api/gost/profile/'
GO_ORGANISM = 'hsapiens'
GO_FDR_THRESHOLD = 0.05
GO_TOP_TERMS = 15
GO_TOP_GENES = 100
GO_GENE_CUTOFFS = (25, 50, 100)
ALLOWED_MODALITIES = {'RNA', 'ATAC'}
MAJOR_CELLTYPES = {'Neuron', 'Micro', 'OPC', 'ExN', 'InN', 'Glia', 'NNC'}
MAJOR_SPECIFICITY_PREFIXES = ('1.', '2.')
FILL_MISSING_ENRICHMENT_WITH_ZERO = True
MISSING_NES_VALUE = 0.0
MISSING_PADJ_VALUE = 1.0
INCLUDE_MHC_GENES = True
RUN_MHC_SENSITIVITY = True
MHC_GENE_SYMBOLS = {'KIFC1', 'TAP2', 'POU5F1', 'HLA-DQA2', 'TRIM10', 'BTNL2', 'HLA-T', 'CDSN', 'TCF19', 'HLA-W', 'LTB', 'TNXB', 'HLA-DRB6', 'HLA-C', 'AGPAT1', 'HLA-A', 'RGL2', 'PPP1R10', 'HLA-DRB5', 'HLA-L', 'HLA-B', 'GPX5', 'HLA-DOB', 'PPT2', 'LTA', 'PSORS1C1', 'HLA-DOA', 'HLA-P', 'TAP1', 'RNF5', 'VARS', 'CFB', 'C6ORF47', 'HCG27', 'HLA-DRA', 'DDAH2', 'PRRT1', 'PSMB9', 'VARS1', 'TRIM31', 'C2', 'NFKBIL1', 'TRIM26', 'C6ORF10', 'HLA-DQA1', 'NCR3', 'SFTA2', 'HLA-K', 'GNL1', 'PPT2-EGFL8', 'HLA-DQB2', 'HLA-V', 'MRPS18B', 'EHMT2', 'DHX16', 'TRIM40', 'PPP1R11', 'STK19', 'ATP6V1G2', 'C4A', 'HLA-H', 'HLA-DQB1', 'C6ORF48', 'HLA-F', 'ZBTB12', 'TNF', 'HLA-G', 'HLA-DPB2', 'TRIM15', 'HLA-E', 'HLA-S', 'CSNK2B', 'HLA-DPB1', 'CYP21A2', 'MICA', 'TRIM39', 'BAG6', 'HLA-DMB', 'NOTCH4', 'C4B', 'MICB', 'MSH5', 'SKIV2L', 'HLA-DPA1', 'PRRC2A', 'PSMB8', 'TAPBP', 'EGFL8', 'HLA-J', 'FKBPL', 'PSMB8-AS1', 'VWA7', 'GPSM3', 'HLA-N', 'ABHD16A', 'HLA-DRB1', 'PBX2', 'MSH5-SAPCD1', 'AGER', 'HLA-DRB9', 'HLA-DQB3', 'LST1', 'AIF1', 'APOM', 'HLA-U', 'HLA-DMA'}
MHC_GENE_PREFIXES = ("HLA-",)



def is_mhc_gene(gene):
    gene = str(gene).strip().upper()
    if gene in MHC_GENE_SYMBOLS:
        return True
    return gene.startswith(MHC_GENE_PREFIXES)

def load_s5_wide(path=WORKBOOK_PATH, sheet=SHEET_S5, allowed_celltypes=None, specificity_prefixes=None):
    allowed_celltypes = allowed_celltypes or MAJOR_CELLTYPES
    specificity_prefixes = specificity_prefixes or MAJOR_SPECIFICITY_PREFIXES

    raw = pd.read_excel(path, sheet_name=sheet, header=None)

    header_row_idx = None
    for i in range(min(10, len(raw))):
        row_vals = raw.iloc[i].astype(str).str.strip()
        if (row_vals == "Celltype").any() and (row_vals == "Trait").any():
            header_row_idx = i
            break
    if header_row_idx is None:
        raise ValueError("Could not find S5 header row.")

    header_row = raw.iloc[header_row_idx].astype(str).str.strip()

    title_row_idx = None
    for i in range(header_row_idx - 1, -1, -1):
        if raw.iloc[i].notna().any():
            title_row_idx = i
            break
    titles = raw.iloc[title_row_idx].ffill() if title_row_idx is not None else None

    panel_starts = [c for c in range(raw.shape[1]) if header_row.iloc[c] == "Celltype"]
    panel_starts.append(raw.shape[1])

    panels = []
    for p in range(len(panel_starts) - 1):
        c0, c1 = panel_starts[p], panel_starts[p + 1]
        block_headers = header_row.iloc[c0:c1].tolist()

        seen = {}
        clean_headers = []
        blank_cols = []
        for j, h in enumerate(block_headers):
            h = "" if h == "nan" else h
            if h == "":
                name = f"__blank_{j}"
                blank_cols.append(name)
            else:
                seen[h] = seen.get(h, 0) + 1
                name = h if seen[h] == 1 else f"{h}__dup{seen[h]}"
            clean_headers.append(name)

        if S5_CELLTYPE_COL not in clean_headers:
            continue

        block = raw.iloc[header_row_idx + 1:, c0:c1].copy()
        block.columns = clean_headers
        block = block.drop(columns=blank_cols, errors="ignore")
        block = block.dropna(how="all")
        block = block[block[S5_CELLTYPE_COL].notna()]
        if block.empty:
            continue

        label = str(titles.iloc[c0]) if titles is not None else f"panel_{p+1}"
        block["specificity"] = label
        panels.append(block)

    df = pd.concat(panels, ignore_index=True)

    df["specificity"] = df["specificity"].astype(str).str.strip()
    raw_celltypes = sorted(df[S5_CELLTYPE_COL].astype(str).str.strip().unique())
    raw_panels = sorted(df["specificity"].unique())
    n_before = len(df)

    df = df[df["specificity"].str.startswith(specificity_prefixes)]
    df = df[df[S5_CELLTYPE_COL].astype(str).str.strip().isin(allowed_celltypes)]
    df[S5_MODALITY_COL] = df[S5_MODALITY_COL].astype(str).str.strip().str.upper()
    df = df[df[S5_MODALITY_COL].isin(ALLOWED_MODALITIES)]

    kept_celltypes = sorted(df[S5_CELLTYPE_COL].astype(str).str.strip().unique())
    missing_allowed = set(allowed_celltypes) - set(kept_celltypes)

    print(f"[load_s5_wide] celltype/panel filter: kept {len(df)}/{n_before} rows.")
    print(f"  All celltypes found in sheet before filtering: {raw_celltypes}")
    print(f"  All specificity panels found in sheet before filtering: {raw_panels}")
    print(f"  Celltypes retained after filtering: {kept_celltypes}")
    if missing_allowed:
        print(f"  [warning] these allowed_celltypes had ZERO matching rows - check spelling against the raw celltype list above: {missing_allowed}")

    if df.empty:
        raise ValueError(
            "No rows left after celltype/specificity filtering - allowed_celltypes or "
            "ALLOWED_SPECIFICITY_PREFIXES likely don't match your sheet's actual labels."
        )

    for col in (S5_NES_COL, S5_PADJ_COL, "pval"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    # ------------------------------------------------------------
    # Missing-enrichment rule: NA NES means no detected enrichment.
    # Fill NES=0 and padj=1.0 BEFORE calculating signed_score.
    # ------------------------------------------------------------
    row_level_missing = df[[S5_NES_COL, S5_PADJ_COL]].isna().any(axis=1)
    if FILL_MISSING_ENRICHMENT_WITH_ZERO:
        df[S5_NES_COL] = df[S5_NES_COL].fillna(MISSING_NES_VALUE)
        df[S5_PADJ_COL] = df[S5_PADJ_COL].fillna(MISSING_PADJ_VALUE)
    else:
        row_level_missing = row_level_missing.astype(bool)

    df["feature"] = (
        df[S5_CELLTYPE_COL].astype(str) + "_" +
        df[S5_STAGE_COL].astype(str) + "_" +
        df[S5_MODALITY_COL].astype(str) + "_" +
        df["specificity"].astype(str)
    )

    # signed_score = sign(NES) * -log10(padj).
    # With NES=0 and padj=1.0, the score is exactly 0.
    df["signed_score"] = np.sign(df[S5_NES_COL]) * -np.log10(
        df[S5_PADJ_COL].clip(lower=1e-300)
    )
    df["signed_score"] = df["signed_score"].clip(-20, 20)
    df["was_row_level_missing_enrichment"] = row_level_missing

    wide_nes_raw = df.pivot_table(
        index=S5_TRAIT_COL, columns="feature", values=S5_NES_COL, aggfunc="mean"
    )
    wide_signed_raw = df.pivot_table(
        index=S5_TRAIT_COL, columns="feature", values="signed_score", aggfunc="mean"
    )

    # ------------------------------------------------------------
    # Some trait x feature combinations may have no row at all in the
    # source table. Under the user's assumption, those also mean
    # "no enrichment detected": NES=0, padj=1.0, signed_score=0.
    # Save an audit log before filling.
    # ------------------------------------------------------------
    complete_features = sorted(set(wide_nes_raw.columns) | set(wide_signed_raw.columns))
    wide_nes_raw = wide_nes_raw.reindex(columns=complete_features)
    wide_signed_raw = wide_signed_raw.reindex(columns=complete_features)

    na_mask = wide_nes_raw.isna() | wide_signed_raw.isna()
    na_rows = []
    for trait in wide_nes_raw.index:
        for feature in complete_features:
            if bool(na_mask.loc[trait, feature]):
                parts = str(feature).split("_")
                na_rows.append({
                    "Trait": trait,
                    "feature": feature,
                    "original_NES": np.nan,
                    "replacement_NES": MISSING_NES_VALUE,
                    "replacement_padj": MISSING_PADJ_VALUE,
                    "replacement_signed_score": 0.0,
                    "reason": "no enrichment row / NA enrichment result",
                })

    na_audit = pd.DataFrame(na_rows)

    if FILL_MISSING_ENRICHMENT_WITH_ZERO:
        wide_nes = wide_nes_raw.fillna(MISSING_NES_VALUE)
        wide_signed = wide_signed_raw.fillna(0.0)
    else:
        wide_nes = wide_nes_raw.copy()
        wide_signed = wide_signed_raw.copy()

    print(
        f"[load_s5_wide] missing-enrichment audit: {len(na_audit)} trait x feature cells "
        f"encoded as NES={MISSING_NES_VALUE}, padj={MISSING_PADJ_VALUE}, signed_score=0."
        if FILL_MISSING_ENRICHMENT_WITH_ZERO
        else f"[load_s5_wide] missing-enrichment audit: {len(na_audit)} cells retained as NA."
    )

    # Store the audit table as an attribute for optional saving without
    # changing the historical three-object return signature.
    wide_nes.attrs["missing_enrichment_audit"] = na_audit
    wide_signed.attrs["missing_enrichment_audit"] = na_audit

    return wide_nes, wide_signed, df

def leadingedge_to_gene_sets(
    df,
    trait_col=S5_TRAIT_COL,
    le_col=S5_LEADINGEDGE_COL,
    padj_col=S5_PADJ_COL,
    sig_threshold=0.05,
    min_genes=10,
    drop_mhc=True,
):
    sig_df = df[df[padj_col] < sig_threshold]

    gene_sets = {}
    for trait, sub in sig_df.groupby(trait_col):
        genes = set()
        for entry in sub[le_col].dropna():
            for g in re.split(r"[,;]\s*", str(entry)):
                g = g.strip()
                if not g:
                    continue
                if drop_mhc and is_mhc_gene(g):
                    continue
                genes.add(g)

        if len(genes) >= min_genes:
            gene_sets[trait] = genes

    return gene_sets

def _nmf_binary_matrix(gene_sets):
    all_genes = sorted(set().union(*gene_sets.values()))
    traits = list(gene_sets.keys())
    M = pd.DataFrame(0, index=traits, columns=all_genes, dtype=float)
    for t in traits:
        M.loc[t, list(gene_sets[t])] = 1.0
    return M

def consensus_nmf(gene_sets, n_programs, n_runs=50, random_state=RANDOM_STATE):
    M = _nmf_binary_matrix(gene_sets)
    traits = M.index.tolist()
    n = len(traits)

    co_occurrence = np.zeros((n, n))
    recon_errors = []
    all_W = []

    rng = np.random.RandomState(random_state)
    for run in range(n_runs):
        seed = rng.randint(0, 1_000_000)
        model = NMF(n_components=n_programs, init="random", random_state=seed, max_iter=500)
        W = model.fit_transform(M.values)
        recon_errors.append(model.reconstruction_err_)
        all_W.append(W)

        dominant = np.argmax(W, axis=1)
        same = (dominant[:, None] == dominant[None, :]).astype(float)
        co_occurrence += same

    consensus = co_occurrence / n_runs
    consensus_df = pd.DataFrame(consensus, index=traits, columns=traits)
    return consensus_df, float(np.mean(recon_errors)), all_W

def nmf_gene_programs_stable(gene_sets, n_programs, n_runs=50, random_state=RANDOM_STATE, max_iter=1000):
    """
    Run multiple random NMF restarts and KEEP THE BEST RUN (lowest
    reconstruction error). Returns the exact W and H from that same run.

    This is important for downstream interpretation: top genes, GO
    enrichment, and disease loadings all come from the same NMF solution,
    so program IDs cannot be accidentally permuted by a second NMF fit.
    """
    M = _nmf_binary_matrix(gene_sets)
    all_genes = M.columns.tolist()
    traits = M.index.tolist()

    best_err = np.inf
    best_W, best_H = None, None
    rng = np.random.RandomState(random_state)

    for _ in range(n_runs):
        seed = int(rng.randint(0, 1_000_000_000))
        model = NMF(
            n_components=n_programs,
            init="random",
            random_state=seed,
            max_iter=max_iter,
        )
        W = model.fit_transform(M.values)
        H = model.components_
        if model.reconstruction_err_ < best_err:
            best_err = model.reconstruction_err_
            best_W = W.copy()
            best_H = H.copy()

    if best_W is None or best_H is None:
        raise RuntimeError("NMF failed to produce a valid solution.")

    program_names = [f"program_{i+1}" for i in range(n_programs)]
    W_df = pd.DataFrame(best_W, index=traits, columns=program_names)

    top_genes_per_program = {}
    for i, prog in enumerate(program_names):
        top_idx = np.argsort(best_H[i])[::-1][:GO_TOP_GENES]
        top_genes_per_program[prog] = [all_genes[j] for j in top_idx]

    consensus_df, mean_err, _ = consensus_nmf(
        gene_sets,
        n_programs,
        n_runs=n_runs,
        random_state=random_state,
    )

    return W_df, best_H, top_genes_per_program, consensus_df, float(best_err), float(mean_err)

def _align_components_to_reference(H_ref, H):
    sim = cosine_similarity(H_ref, H)
    ref_idx, cur_idx = linear_sum_assignment(-sim)
    perm = np.empty(H.shape[0], dtype=int)
    perm[ref_idx] = cur_idx
    matched_sims = sim[ref_idx, cur_idx]
    return perm, matched_sims

def nmf_restart_stability(gene_sets, k, n_runs=100, random_state=0, max_iter=1000):
    M = _nmf_binary_matrix(gene_sets)
    trait_names = M.index.tolist()
    program_names = [f"program_{i+1}" for i in range(k)]

    rng = np.random.default_rng(random_state)
    runs = []

    for _ in range(n_runs):
        seed = int(rng.integers(0, 1_000_000_000))
        model = NMF(n_components=k, init="random", random_state=seed, max_iter=max_iter)
        W = model.fit_transform(M.values)
        H = model.components_
        runs.append({"W": W, "H": H, "err": model.reconstruction_err_})

    ref = min(runs, key=lambda d: d["err"])["H"]

    aligned_W = []
    aligned_H = []
    matched_sims = []

    for r in runs:
        perm, sims = _align_components_to_reference(ref, r["H"])
        aligned_W.append(r["W"][:, perm])
        aligned_H.append(r["H"][perm, :])
        matched_sims.append(sims)

    aligned_W = np.stack(aligned_W, axis=0)
    aligned_H = np.stack(aligned_H, axis=0)
    matched_sims = np.stack(matched_sims, axis=0)

    aligned_W_norm = aligned_W / (aligned_W.sum(axis=2, keepdims=True) + 1e-12)

    W_mean_norm = aligned_W_norm.mean(axis=0)
    W_sd_norm = aligned_W_norm.std(axis=0)

    dominant = aligned_W_norm.argmax(axis=2)
    dominant_mode = np.array([
        np.bincount(dominant[:, i], minlength=k).argmax()
        for i in range(dominant.shape[1])
    ])
    dominant_stability = np.array([
        (dominant[:, i] == dominant_mode[i]).mean()
        for i in range(dominant.shape[1])
    ])

    component_stability = pd.DataFrame(
        {
            "mean_matched_cosine_to_ref": matched_sims.mean(axis=0),
            "sd_matched_cosine_to_ref": matched_sims.std(axis=0),
        },
        index=program_names,
    )

    return {
        "W_mean_norm": pd.DataFrame(W_mean_norm, index=trait_names, columns=program_names),
        "W_sd_norm": pd.DataFrame(W_sd_norm, index=trait_names, columns=program_names),
        "dominant_program_stability": pd.Series(dominant_stability, index=trait_names, name="dominant_program_stability"),
        "component_stability": component_stability,
        "mean_reconstruction_error": float(np.mean([r["err"] for r in runs])),
        "all_runs": runs,
    }

def clean_gene_symbols(genes):
    out = []
    for g in genes:
        if g is None:
            continue
        g = str(g).strip()
        if not g or g.upper() in {"NAN", "NONE", "NA"}:
            continue
        out.append(g)
    return sorted(set(out))

def extract_top_genes_from_H(H, gene_names, top_n=100):
    """Top genes for each NMF program from the EXACT final H matrix."""
    gene_names = list(gene_names)
    top = {}
    for i in range(H.shape[0]):
        idx = np.argsort(H[i])[::-1][:top_n]
        top[f"program_{i+1}"] = [gene_names[j] for j in idx]
    return top

def run_go_enrichment(query_genes, background_genes,
                      fdr=GO_FDR_THRESHOLD, top_terms=GO_TOP_TERMS):
    """GO Biological Process enrichment via g:Profiler with a custom background."""
    query_genes = clean_gene_symbols(query_genes)
    background_genes = clean_gene_symbols(background_genes)

    bg_set = set(background_genes)
    query_genes = [g for g in query_genes if g in bg_set]

    print(f"GO input: {len(query_genes)} query genes / {len(background_genes)} background genes")

    if len(query_genes) < 3 or len(background_genes) < 10:
        print("GO enrichment skipped: insufficient query/background genes.")
        return pd.DataFrame()

    payload = {
        "organism": GO_ORGANISM,
        "query": query_genes,
        "sources": ["GO:BP"],
        "domain_scope": "custom",
        "background": background_genes,
        "significance_threshold_method": "fdr",
        "user_threshold": float(fdr),
        "all_results": False,
        "no_evidences": True,
        "output": "json",
    }

    try:
        response = requests.post(
            GO_API_URL,
            json=payload,
            headers={"User-Agent": "disease-NMF-GO-analysis/2.0"},
            timeout=120,
        )
    except Exception as exc:
        print("GO enrichment request failed:")
        print(repr(exc))
        return pd.DataFrame()

    if response.status_code != 200:
        print(f"GO API returned HTTP {response.status_code}")
        print("Server response:")
        print(response.text[:3000])
        return pd.DataFrame()

    try:
        result = response.json()
    except Exception as exc:
        print("Could not parse GO API response as JSON:", repr(exc))
        return pd.DataFrame()

    if "result" not in result:
        print("GO API response has no 'result' field.")
        print("Response keys:", list(result.keys()))
        return pd.DataFrame()

    df = pd.DataFrame(result["result"])
    if df.empty:
        print("No significant GO Biological Process terms found.")
        return df

    keep = [
        c for c in [
            "source", "native", "name", "description", "p_value",
            "term_size", "intersection_size", "effective_domain_size",
            "intersections"
        ] if c in df.columns
    ]
    df = df[keep].copy()

    if "p_value" in df.columns:
        df = df.sort_values("p_value", ascending=True)

    return df.head(top_terms).reset_index(drop=True)

def _theme_scores(go_df):
    if go_df is None or go_df.empty or "name" not in go_df.columns:
        return {}
    text = " ".join(go_df["name"].astype(str).str.lower().tolist())
    groups = {
        "Immune / antigen-presentation signaling": [
            "immune", "antigen", "leukocyte", "myeloid", "microgl",
            "interferon", "cytokine", "inflamm", "complement",
        ],
        "Synaptic / neurotransmission": [
            "synapse", "synaptic", "neurotransmitter", "presynaptic",
            "postsynaptic", "vesicle",
        ],
        "Neuronal development / connectivity": [
            "neuron development", "neurogenesis", "axon", "neuron projection",
            "cell adhesion", "neuronal", "neuron differentiation",
            "nervous system development",
        ],
        "Chromatin / histone biology": [
            "chromatin", "histone", "nucleosome", "dna packaging",
            "chromosome organization", "epigenetic",
        ],
        "RNA processing / translation": [
            "rna processing", "mrna", "splicing", "ribosome", "translation",
        ],
        "Cell cycle / proliferation": [
            "cell cycle", "mitotic", "cell division", "chromosome segregation",
            "dna replication",
        ],
    }
    return {label: sum(text.count(term) for term in terms) for label, terms in groups.items()}

def propose_go_label(go_by_cutoff):
    """
    Choose a label only when the same broad biological theme is supported
    across multiple gene cutoffs. This avoids overinterpreting a single GO term.
    """
    labels = []
    for cutoff, go_df in go_by_cutoff.items():
        scores = _theme_scores(go_df)
        if scores:
            best = max(scores, key=scores.get)
            if scores[best] >= 2:
                labels.append(best)

    if len(labels) >= 2:
        counts = pd.Series(labels).value_counts()
        best_label = counts.index[0]
        if counts.iloc[0] >= 2:
            return best_label, "robust_across_cutoffs"

    # If no broad theme is stable, use the strongest GO term at 100 genes.
    go100 = go_by_cutoff.get(100, pd.DataFrame())
    if go100 is not None and not go100.empty and "name" in go100.columns:
        return str(go100.iloc[0]["name"]), "top_GO_term_only"

    return "Unresolved / inspect top genes", "unresolved"

def go_label_programs_exact_H(
    H,
    gene_sets,
    out_dir,
    label_prefix,
    cutoffs=GO_GENE_CUTOFFS,
):
    """
    IMPORTANT: GO is run directly on the EXACT H matrix from the final NMF run.
    No second NMF fit is performed, so program IDs cannot be permuted between
    NMF and GO.
    """
    M = _nmf_binary_matrix(gene_sets)
    all_genes = M.columns.tolist()
    background = clean_gene_symbols(all_genes)

    all_top = extract_top_genes_from_H(H, all_genes, top_n=max(cutoffs))

    pd.DataFrame({
        "program": list(all_top.keys()),
        f"top_{max(cutoffs)}_genes": [", ".join(v) for v in all_top.values()],
    }).to_csv(
        f"{out_dir}/{label_prefix}_top_{max(cutoffs)}_genes_from_exact_H.csv",
        index=False,
    )

    labels = []
    combined = []

    for program, genes100 in all_top.items():
        print("\n" + "-" * 100)
        print(f"GO ENRICHMENT: {label_prefix} / {program}")
        print("-" * 100)

        go_by_cutoff = {}
        for cutoff in cutoffs:
            genes = genes100[:cutoff]
            print(f"\n--- Top {cutoff} genes ---")
            go_df = run_go_enrichment(genes, background)
            go_by_cutoff[cutoff] = go_df

            go_df.to_csv(
                f"{out_dir}/{label_prefix}_{program}_GO_BP_top{cutoff}.csv",
                index=False,
            )

            if not go_df.empty and "name" in go_df.columns:
                print(go_df[[c for c in ["name", "p_value", "intersection_size"] if c in go_df.columns]].head(10).to_string(index=False))
            else:
                print("No significant GO Biological Process terms.")

            if not go_df.empty:
                tmp = go_df.copy()
                tmp.insert(0, "gene_cutoff", cutoff)
                tmp.insert(0, "program", program)
                tmp.insert(0, "analysis", label_prefix)
                combined.append(tmp)

        proposed, evidence = propose_go_label(go_by_cutoff)
        print(f"\nPROPOSED ROBUST LABEL: {proposed}")
        print(f"Label evidence: {evidence}")

        labels.append({
            "analysis": label_prefix,
            "program": program,
            "proposed_label": proposed,
            "label_evidence": evidence,
            "n_genes_top100": len(genes100),
            "background_genes": len(background),
            "n_significant_terms_top25": len(go_by_cutoff.get(25, pd.DataFrame())),
            "n_significant_terms_top50": len(go_by_cutoff.get(50, pd.DataFrame())),
            "n_significant_terms_top100": len(go_by_cutoff.get(100, pd.DataFrame())),
        })

    combined_df = pd.concat(combined, ignore_index=True) if combined else pd.DataFrame()
    combined_df.to_csv(f"{out_dir}/{label_prefix}_all_GO_BP_by_cutoff.csv", index=False)

    labels_df = pd.DataFrame(labels)
    labels_df.to_csv(f"{out_dir}/{label_prefix}_GO_based_program_labels.csv", index=False)

    return all_top, combined_df, labels_df

def align_programs_across_mhc_runs(H_in, genes_in, H_out, genes_out):
    """Align MHC-included vs MHC-excluded programs using common genes only."""
    common = sorted(set(genes_in) & set(genes_out))
    if len(common) < 10:
        return pd.DataFrame(), common

    idx_in = [genes_in.index(g) for g in common]
    idx_out = [genes_out.index(g) for g in common]

    A = H_in[:, idx_in]
    B = H_out[:, idx_out]
    sim = cosine_similarity(A, B)
    ref_idx, cur_idx = linear_sum_assignment(-sim)

    rows = []
    for r, c in zip(ref_idx, cur_idx):
        rows.append({
            "MHC_included_program": f"program_{r+1}",
            "MHC_excluded_program": f"program_{c+1}",
            "cosine_similarity_common_genes": float(sim[r, c]),
        })
    return pd.DataFrame(rows), common

def make_unique_labels(labels):
    counts = {}
    out = []
    for lab in labels:
        counts[lab] = counts.get(lab, 0) + 1
        out.append(lab if counts[lab] == 1 else f"{lab} ({counts[lab]})")
    return out

def plot_program_loadings_stacked(W_df, outpath):
    normalized = W_df.div(W_df.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)

    dominant = normalized.idxmax(axis=1)
    dom_strength = normalized.max(axis=1)

    order = (
        pd.DataFrame({"_dom": dominant, "_dom_strength": dom_strength}, index=normalized.index)
        .sort_values(by=["_dom", "_dom_strength"], ascending=[True, False])
        .index
    )
    normalized = normalized.loc[order]

    fig, ax = plt.subplots(figsize=(8, max(6, len(normalized) * 0.3)))
    left = np.zeros(len(normalized))
    colors = sns.color_palette("Set2", normalized.shape[1])

    for i, col in enumerate(normalized.columns):
        ax.barh(normalized.index, normalized[col], left=left, color=colors[i], label=col)
        left += normalized[col].values

    ax.set_xlabel("Fraction of loading")
    ax.set_title("Disease composition across regulatory programs")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left")
    ax.set_xlim(0, 1)
    plt.tight_layout()
    plt.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close()



def run_one(s5, output_dir, include_mhc):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    label = 'MHC_included' if include_mhc else 'MHC_excluded'
    gene_sets = leadingedge_to_gene_sets(
        s5, sig_threshold=0.05, min_genes=10, drop_mhc=not include_mhc,
    )
    M = _nmf_binary_matrix(gene_sets)
    print(f'\n{label}: {M.shape[0]} traits x {M.shape[1]} genes')
    M.to_csv(output_dir / 'binary_trait_gene_matrix.csv')
    W, H, top, consensus, best_err, mean_err = nmf_gene_programs_stable(
        gene_sets, n_programs=NMF_K_FIXED, n_runs=NMF_STABILITY_RUNS,
        random_state=RANDOM_STATE, max_iter=1000,
    )
    H_df = pd.DataFrame(H, index=W.columns, columns=M.columns)
    W.to_csv(output_dir / 'W_trait_program_loadings.csv')
    H_df.to_csv(output_dir / 'H_program_gene_loadings.csv')
    consensus.to_csv(output_dir / 'NMF_consensus_matrix.csv')
    W_norm = W.div(W.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
    W_norm.to_csv(output_dir / 'W_normalized_trait_program_loadings.csv')
    plot_program_loadings_stacked(W, str(output_dir / 'normalized_program_loadings.pdf'))
    pd.DataFrame({p:pd.Series(g) for p,g in top.items()}).to_csv(output_dir / 'top100_genes_by_program.csv', index=False)

    stability = nmf_restart_stability(
        gene_sets, k=NMF_K_FIXED, n_runs=NMF_STABILITY_RUNS,
        random_state=RANDOM_STATE, max_iter=1000,
    )
    stability['component_stability'].to_csv(output_dir / 'component_stability.csv')
    stability['dominant_program_stability'].to_csv(output_dir / 'trait_assignment_stability.csv')
    stability['W_mean_norm'].to_csv(output_dir / 'W_mean_normalized_across_restarts.csv')
    stability['W_sd_norm'].to_csv(output_dir / 'W_SD_normalized_across_restarts.csv')
    pd.DataFrame({'reconstruction_error':[r['err'] for r in stability['all_runs']]}).to_csv(
        output_dir / 'restart_reconstruction_errors.csv', index=False)

    _, go_df, labels = go_label_programs_exact_H(
        H, gene_sets, str(output_dir), label, cutoffs=GO_GENE_CUTOFFS)
    summary = {
        'analysis':label, 'major_celltypes_only':True, 'modalities':['RNA','ATAC'],
        'specificity_panels':['1.','2.'], 'FGSEA_FDR_cutoff':0.05,
        'minimum_leading_edge_genes':10, 'n_traits':int(M.shape[0]),
        'n_genes':int(M.shape[1]), 'nmf_rank_fixed':NMF_K_FIXED,
        'random_state':RANDOM_STATE, 'restarts':NMF_STABILITY_RUNS,
        'maximum_iterations':1000, 'initialization':'random',
        'selection':'minimum reconstruction error across restarts',
        'best_reconstruction_error':best_err,
        'mean_reconstruction_error_best_solution_search':mean_err,
        'component_cosine_mean':stability['component_stability']['mean_matched_cosine_to_ref'].to_dict(),
        'GO_background':'full NMF gene universe', 'GO_FDR_threshold':0.05,
        'GO_top_gene_cutoffs':list(GO_GENE_CUTOFFS),
    }
    (output_dir / 'parameters_and_summary.json').write_text(json.dumps(summary, indent=2))
    print('\nComponent stability:\n', stability['component_stability'])
    print('\nTrait assignment stability:\n', stability['dominant_program_stability'].sort_values())
    return {'H':H,'genes':M.columns.tolist(), 'summary':summary}


def main():
    parser=argparse.ArgumentParser(description='Major cell type RNA/ATAC NMF only')
    parser.add_argument('--workbook', required=True, help='Excel file containing Table S5')
    parser.add_argument('--output', default='./NMF_major_celltypes_only')
    args=parser.parse_args()
    s5=load_s5_wide(args.workbook, allowed_celltypes=MAJOR_CELLTYPES,
                    specificity_prefixes=MAJOR_SPECIFICITY_PREFIXES)[2]
    output=Path(args.output)
    included=run_one(s5, output / 'MHC_included', include_mhc=True)
    if RUN_MHC_SENSITIVITY:
        excluded=run_one(s5, output / 'MHC_excluded', include_mhc=False)
        comparisons, common=align_programs_across_mhc_runs(
            included['H'], included['genes'], excluded['H'], excluded['genes'])
        comparisons.to_csv(output / 'MHC_program_alignment.csv', index=False)
        print('\nMHC sensitivity comparison:\n',comparisons)
    print(f'\nSaved NMF-only outputs to: {output}')

if __name__ == '__main__':
    main()
