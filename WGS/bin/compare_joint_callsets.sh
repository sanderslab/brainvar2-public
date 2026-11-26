#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'USAGE'
Usage:
  compare_joint_callsets.sh \
    --candidate CANDIDATE.vcf.gz \
    --reference REFERENCE.vcf.gz \
    [--output comparison.tsv]

Compare sample order, variant keys, genotype calls, and complete per-sample
FORMAT payloads without writing sample identifiers or genotype rows to disk.
The command exits nonzero when any reported metric differs.
USAGE
}

candidate=''
reference=''
output=''

while (($#)); do
    case "$1" in
        --candidate)
            candidate=${2:?missing value for --candidate}
            shift 2
            ;;
        --reference)
            reference=${2:?missing value for --reference}
            shift 2
            ;;
        --output)
            output=${2:?missing value for --output}
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "ERROR: unknown argument: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

[[ -n "$candidate" ]] || {
    echo 'ERROR: --candidate is required' >&2
    exit 2
}
[[ -n "$reference" ]] || {
    echo 'ERROR: --reference is required' >&2
    exit 2
}

for command_name in bcftools sha256sum awk; do
    command -v "$command_name" >/dev/null || {
        echo "ERROR: required command not found: $command_name" >&2
        exit 2
    }
done

for vcf in "$candidate" "$reference"; do
    [[ -r "$vcf" ]] || {
        echo "ERROR: VCF is not readable: $vcf" >&2
        exit 2
    }
done

export LC_ALL=C

sample_count() {
    bcftools query -l "$1" | awk 'END { print NR }'
}

record_count() {
    bcftools query -f '%CHROM\n' "$1" | awk 'END { print NR }'
}

stream_sha256() {
    sha256sum | awk '{ print $1 }'
}

sample_order_sha256() {
    bcftools query -l "$1" | stream_sha256
}

record_genotype_sha256() {
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' "$1" |
        stream_sha256
}

record_format_sha256() {
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FORMAT\n' "$1" |
        stream_sha256
}

declare -A candidate_metrics reference_metrics
metrics=(
    sample_count
    record_count
    sample_order_sha256
    record_genotype_sha256
    record_format_sha256
)

candidate_metrics[sample_count]=$(sample_count "$candidate")
reference_metrics[sample_count]=$(sample_count "$reference")
candidate_metrics[record_count]=$(record_count "$candidate")
reference_metrics[record_count]=$(record_count "$reference")
candidate_metrics[sample_order_sha256]=$(sample_order_sha256 "$candidate")
reference_metrics[sample_order_sha256]=$(sample_order_sha256 "$reference")
candidate_metrics[record_genotype_sha256]=$(record_genotype_sha256 "$candidate")
reference_metrics[record_genotype_sha256]=$(record_genotype_sha256 "$reference")
candidate_metrics[record_format_sha256]=$(record_format_sha256 "$candidate")
reference_metrics[record_format_sha256]=$(record_format_sha256 "$reference")

if [[ -n "$output" ]]; then
    output_dir=$(dirname "$output")
    [[ -d "$output_dir" ]] || {
        echo "ERROR: output directory does not exist: $output_dir" >&2
        exit 2
    }
    exec >"$output"
fi

printf 'metric\tcandidate\treference\tmatch\n'
all_match=true
for metric in "${metrics[@]}"; do
    match=false
    if [[ "${candidate_metrics[$metric]}" == "${reference_metrics[$metric]}" ]]; then
        match=true
    else
        all_match=false
    fi
    printf '%s\t%s\t%s\t%s\n' \
        "$metric" \
        "${candidate_metrics[$metric]}" \
        "${reference_metrics[$metric]}" \
        "$match"
done

[[ "$all_match" == true ]]
