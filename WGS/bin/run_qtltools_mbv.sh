#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: run_qtltools_mbv.sh --bam FILE --vcf FILE --out FILE [options]

Required:
  --bam FILE              Indexed RNA-seq BAM/SAM input
  --vcf FILE              Indexed genotype VCF/BCF input
  --out FILE              QTLtools MBV report path

Optional:
  --reg REGION            Region accepted by QTLtools --reg
  --include-samples FILE  Genotype sample inclusion list
  --exclude-samples FILE  Genotype sample exclusion list
  --log FILE              QTLtools log path
  --seed INTEGER          Random seed (QTLtools default: 15112011)
  --qtltools COMMAND      QTLtools executable (default: QTLtools)
  --dry-run               Print the shell-escaped command without running it
  -h, --help              Show this help
EOF
}

bam=''
vcf=''
out=''
reg=''
include_samples=''
exclude_samples=''
log=''
seed=''
qtltools='QTLtools'
dry_run=false

while (($#)); do
    case "$1" in
        --bam|--vcf|--out|--reg|--include-samples|--exclude-samples|--log|--seed|--qtltools)
            (($# >= 2)) || { echo "ERROR: $1 requires a value" >&2; exit 2; }
            key=${1#--}
            key=${key//-/_}
            printf -v "$key" '%s' "$2"
            shift 2
            ;;
        --dry-run)
            dry_run=true
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "ERROR: unknown option: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

[[ -n "$bam" && -n "$vcf" && -n "$out" ]] || {
    echo 'ERROR: --bam, --vcf, and --out are required' >&2
    usage >&2
    exit 2
}

command=("$qtltools" mbv --bam "$bam" --vcf "$vcf" --out "$out")
[[ -z "$reg" ]] || command+=(--reg "$reg")
[[ -z "$include_samples" ]] || command+=(--include-samples "$include_samples")
[[ -z "$exclude_samples" ]] || command+=(--exclude-samples "$exclude_samples")
[[ -z "$log" ]] || command+=(--log "$log")
[[ -z "$seed" ]] || command+=(--seed "$seed")

if [[ "$dry_run" == true ]]; then
    printf '%q ' "${command[@]}"
    printf '\n'
    exit 0
fi

command -v "$qtltools" >/dev/null || {
    echo "ERROR: QTLtools executable not found: $qtltools" >&2
    exit 1
}
"${command[@]}"
