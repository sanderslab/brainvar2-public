#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: run_mbv.sh --vcf FILE --outdir DIR [options] BAM [BAM ...]

Run QTLtools MBV once per indexed RNA BAM against one indexed cohort VCF.

Required:
  --vcf FILE              Whole-genome genotype VCF/BCF
  --outdir DIR            Directory for per-BAM reports and logs
  BAM [BAM ...]           RNA alignments; shell globs are supported

Options:
  --mapping-quality INT   Minimum read mapping quality (default: 15)
  --qtltools COMMAND      QTLtools executable (default: QTLtools)
  --dry-run               Print commands without running QTLtools
  -h, --help              Show this help
EOF
}

vcf=''
outdir=''
mapping_quality='15'
qtltools='QTLtools'
dry_run=false
bams=()

while (($#)); do
    case "$1" in
        --vcf) vcf=${2:?"ERROR: --vcf requires a value"}; shift 2 ;;
        --outdir) outdir=${2:?"ERROR: --outdir requires a value"}; shift 2 ;;
        --mapping-quality) mapping_quality=${2:?"ERROR: --mapping-quality requires a value"}; shift 2 ;;
        --qtltools) qtltools=${2:?"ERROR: --qtltools requires a value"}; shift 2 ;;
        --dry-run) dry_run=true; shift ;;
        -h|--help) usage; exit 0 ;;
        -*) echo "ERROR: unknown option: $1" >&2; usage >&2; exit 2 ;;
        *) bams+=("$1"); shift ;;
    esac
done

[[ -n "$vcf" && -n "$outdir" && ${#bams[@]} -gt 0 ]] || {
    echo 'ERROR: --vcf, --outdir, and at least one BAM are required' >&2
    usage >&2
    exit 2
}
[[ -f "$vcf" ]] || { echo "ERROR: VCF/BCF not found: $vcf" >&2; exit 2; }
[[ -f "$vcf.tbi" || -f "$vcf.csi" ]] || {
    echo "ERROR: index not found beside $vcf" >&2
    exit 2
}
[[ "$dry_run" == true ]] || command -v "$qtltools" >/dev/null || {
    echo "ERROR: QTLtools executable not found: $qtltools" >&2
    exit 2
}

mkdir -p "$outdir"

for bam in "${bams[@]}"; do
    [[ -f "$bam" ]] || { echo "ERROR: BAM not found: $bam" >&2; exit 2; }
    [[ "$bam" == *.bam ]] || { echo "ERROR: expected a .bam file: $bam" >&2; exit 2; }
    [[ -f "$bam.bai" || -f "${bam%.bam}.bai" ]] || {
        echo "ERROR: index not found beside $bam" >&2
        exit 2
    }

    name=$(basename "$bam")
    report="$outdir/$name.mbv.txt"
    log="$outdir/$name.mbv.log"
    [[ ! -e "$report" && ! -e "$log" ]] || {
        echo "ERROR: output already exists for $name" >&2
        exit 2
    }

    command=(
        "$qtltools" mbv
        --bam "$bam"
        --vcf "$vcf"
        --out "$report"
        --log "$log"
        --filter-mapping-quality "$mapping_quality"
    )
    printf '%q ' "${command[@]}"
    printf '\n'
    [[ "$dry_run" == true ]] || "${command[@]}"
done
