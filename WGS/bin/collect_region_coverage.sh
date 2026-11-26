#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: collect_region_coverage.sh --alignment FILE --region REGION --out FILE [options]

Required:
  --alignment FILE    Indexed BAM/CRAM input
  --region REGION     Region accepted by samtools coverage
  --out FILE          Tabular coverage output

Optional:
  --reference FILE    Reference FASTA for CRAM decoding
  --min-mq INTEGER    Minimum mapping quality (default: 0)
  --min-bq INTEGER    Minimum base quality (default: 0)
  --samtools COMMAND  samtools executable (default: samtools)
  --dry-run           Print the shell-escaped command without running it
  -h, --help          Show this help
EOF
}

alignment=''
region=''
out=''
reference=''
min_mq=0
min_bq=0
samtools='samtools'
dry_run=false

while (($#)); do
    case "$1" in
        --alignment|--region|--out|--reference|--min-mq|--min-bq|--samtools)
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

[[ -n "$alignment" && -n "$region" && -n "$out" ]] || {
    echo 'ERROR: --alignment, --region, and --out are required' >&2
    usage >&2
    exit 2
}
[[ "$min_mq" =~ ^[0-9]+$ && "$min_bq" =~ ^[0-9]+$ ]] || {
    echo 'ERROR: --min-mq and --min-bq must be non-negative integers' >&2
    exit 2
}

command=("$samtools" coverage --region "$region" --min-MQ "$min_mq" --min-BQ "$min_bq")
[[ -z "$reference" ]] || command+=(--reference "$reference")
command+=(--output "$out" "$alignment")

if [[ "$dry_run" == true ]]; then
    printf '%q ' "${command[@]}"
    printf '\n'
    exit 0
fi

command -v "$samtools" >/dev/null || {
    echo "ERROR: samtools executable not found: $samtools" >&2
    exit 1
}
"${command[@]}"
