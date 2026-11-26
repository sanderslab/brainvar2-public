#!/usr/bin/env bash
set -euo pipefail

project_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$project_dir"

tmpdir=$(mktemp -d "${TMPDIR:-/tmp}/brainvar2-wgs.XXXXXX")
trap 'rm -rf "$tmpdir"' EXIT

python3 -m json.tool nextflow_schema.json >/dev/null
rg -q '"additionalProperties": false' nextflow_schema.json
bash -n scripts/utility/download_resources.sh
bash scripts/utility/download_resources.sh --help >/dev/null
bash -n \
    bin/run_qtltools_mbv.sh \
    bin/collect_region_coverage.sh \
    bin/compare_joint_callsets.sh
bash bin/run_qtltools_mbv.sh --help >/dev/null
bash bin/collect_region_coverage.sh --help >/dev/null
bash bin/compare_joint_callsets.sh --help >/dev/null
python3 tests/test_reconcile_cohort_counts.py
python3 bin/reconcile_cohort_counts.py \
    --manifest tests/fixtures/cohort/joint_manifest.csv \
    --library-metadata tests/fixtures/cohort/library_metadata_v1.tsv \
    --subject-metadata tests/fixtures/cohort/subject_metadata.tsv \
    --draft-library-metadata tests/fixtures/cohort/library_metadata_v2.tsv \
    --output "$tmpdir/cohort_counts.json"
python3 -m json.tool "$tmpdir/cohort_counts.json" >/dev/null

mbv_command=$(bash bin/run_qtltools_mbv.sh \
    --bam reviewer.bam --vcf reviewer.vcf.gz --out reviewer.mbv.txt --dry-run)
rg -q -- 'QTLtools mbv --bam reviewer\.bam --vcf reviewer\.vcf\.gz --out reviewer\.mbv\.txt' \
    <<<"$mbv_command"

coverage_command=$(bash bin/collect_region_coverage.sh \
    --alignment reviewer.cram --region chrY --out reviewer.coverage.tsv --dry-run)
rg -q -- 'samtools coverage .*--region chrY.*reviewer\.cram' <<<"$coverage_command"

bash bin/compare_joint_callsets.sh \
    --candidate tests/fixtures/validation/joint_callset.vcf \
    --reference tests/fixtures/validation/joint_callset.vcf \
    --output "$tmpdir/joint_call_comparison.tsv"
test "$(awk -F '\t' 'NR > 1 && $4 == "true" { count++ } END { print count + 0 }' \
    "$tmpdir/joint_call_comparison.tsv")" -eq 5
if bash bin/compare_joint_callsets.sh \
    --candidate tests/fixtures/validation/joint_callset.vcf \
    --reference tests/fixtures/validation/joint_callset_dp_mismatch.vcf \
    --output "$tmpdir/joint_call_mismatch.tsv"; then
    echo 'ERROR: FORMAT-only mismatch was not detected' >&2
    exit 1
fi
awk -F '\t' '$1 == "record_genotype_sha256" { exit $4 == "true" ? 0 : 1 }' \
    "$tmpdir/joint_call_mismatch.tsv"
awk -F '\t' '$1 == "record_format_sha256" { exit $4 == "false" ? 0 : 1 }' \
    "$tmpdir/joint_call_mismatch.tsv"
python3 bin/validate_samplesheet.py \
    parabricks tests/fixtures/parabricks_samplesheet.csv
python3 bin/validate_samplesheet.py \
    joint tests/fixtures/joint_samplesheet.csv --contigs chr1

if python3 bin/validate_samplesheet.py \
    parabricks tests/fixtures/invalid_parabricks_samplesheet.csv >/dev/null 2>&1; then
    echo 'ERROR: duplicate Parabricks sample-FASTQ record unexpectedly passed' >&2
    exit 1
fi

if python3 bin/validate_samplesheet.py \
    joint tests/fixtures/invalid_joint_samplesheet.csv --contigs chr1,chr2 >/dev/null 2>&1; then
    echo 'ERROR: incomplete joint samplesheet unexpectedly passed' >&2
    exit 1
fi

NXF_OFFLINE=true nextflow run . \
    -stub-run \
    -c conf/test_parabricks.config \
    --outdir "$tmpdir/results-parabricks" \
    -work-dir "$tmpdir/work-parabricks" \
    -ansi-log false

if NXF_OFFLINE=true nextflow run . \
    -stub-run \
    -c conf/test.config \
    --unknown_parameter reject_me \
    --outdir "$tmpdir/results-unknown-param" \
    -work-dir "$tmpdir/work-unknown-param" \
    -ansi-log false >"$tmpdir/unknown-param.log" 2>&1; then
    echo 'ERROR: unknown workflow parameter unexpectedly passed' >&2
    exit 1
fi
rg -q 'Unknown workflow parameter.*--unknown_parameter' "$tmpdir/unknown-param.log"

NXF_OFFLINE=true nextflow run . \
    -stub-run \
    -c conf/test.config \
    --outdir "$tmpdir/results-joint" \
    -work-dir "$tmpdir/work-joint" \
    -ansi-log false

test -e "$tmpdir/results-parabricks/parabricks/gvcf/reviewer_sample_1.grch38.gatk.haplotypecaller.g.vcf.gz"
test -e "$tmpdir/results-joint/joint/final_callset/reviewer_cohort.jointcalled.recalibrated.GRCh38.vcf.gz"
test -e "$tmpdir/results-joint/pipeline_info/software_versions.yml"

rg -q -- '--markdups-picard-version-2182' modules/local/parabricks/fq2bam/main.nf
rg -q -- '--optical-duplicate-pixel-distance 2500' modules/local/parabricks/fq2bam/main.nf
rg -q -- 'task.memory.giga / 2' modules/local/parabricks/fq2bam/main.nf
rg -q -- 'accelerator = 1' nextflow.config
test "$(rg -l -- 'def numGpus = task\.ext\.num_gpus \?: 1' modules/local/parabricks/*/main.nf | wc -l)" -eq 2
rg -q -- 'BWA-MEM 0\.7\.15' README.md
rg -q -- 'GATK 4\.2\.0\.0' README.md
rg -q -- '-K 10000000' README.md
rg -q -- '-G StandardHCAnnotation' modules/local/parabricks/haplotypecaller/main.nf
rg -q -- '--batch-size 50' nextflow.config
rg -q -- "Channel.value\(99.8\)" workflows/brainvar2_wgs.nf
rg -q -- "^cat > sample_map.tsv <<'END_SAMPLE_MAP'$" modules/local/gatk4/genomicsdbimport/main.nf
test "$(rg -c -- '^END_SAMPLE_MAP$' modules/local/gatk4/genomicsdbimport/main.nf)" -eq 1
test "$(rg -c -- "^cat > versions.yml <<'END_VERSIONS'$" modules/local/gatk4/genomicsdbimport/main.nf)" -eq 2
test "$(rg -c -- '^END_VERSIONS$' modules/local/gatk4/genomicsdbimport/main.nf)" -eq 2
test "$(rg -l 'parabricks:4\.0\.1-1@sha256:1aed52b6730a213b8bf7a1aa6998668ad26f421fb49230c04e4340c6f68c4d9e' modules/local/parabricks/*/main.nf | wc -l)" -eq 2
test "$(rg -l 'gatk4:4\.2\.6\.1--hdfd78af_0@sha256:2036bc919ce30eda27195ecd2639d6d556d4af03c04a66fcd5628211b2cd13db' modules/local/gatk4/*/main.nf | wc -l)" -eq 5

git diff --check
git diff --cached --check

required_files=(
    .github/workflows/ci.yml
    README.md LICENSE CHANGELOG.md
    THIRD_PARTY_NOTICES.md nextflow.config nextflow_schema.json main.nf
    docs/methods.md docs/provenance.md docs/cohort_accounting.md
    docs/identity_qc.md docs/verification.md docs/execution_environment.md
    bin/compare_joint_callsets.sh validation/chr21_joint_call_equivalence.tsv
)
for required_file in "${required_files[@]}"; do
    [[ -e "$required_file" ]] || {
        echo "ERROR: required release file is missing: $required_file" >&2
        exit 1
    }
done

mapfile -d '' public_files < <(
    git ls-files -z --cached --others --exclude-standard
)

existing_public_files=()
for public_file in "${public_files[@]}"; do
    [[ -f "$public_file" ]] && existing_public_files+=("$public_file")
done
public_files=("${existing_public_files[@]}")

if ((${#public_files[@]} == 0)); then
    echo 'ERROR: no public release files were found' >&2
    exit 1
fi

scan_files=()
for public_file in "${public_files[@]}"; do
    [[ "$public_file" == scripts/check_release.sh ]] || scan_files+=("$public_file")
done

protected_fastq_found=false
while IFS= read -r tracked_path; do
    case "${tracked_path}" in
        tests/fixtures/fastq/reviewer_sample_1_R1.fastq.gz|\
        tests/fixtures/fastq/reviewer_sample_1_R2.fastq.gz|\
        tests/fixtures/fastq/reviewer_sample_1_chunk2_R1.fastq.gz|\
        tests/fixtures/fastq/reviewer_sample_1_chunk2_R2.fastq.gz)
            ;;
        *)
            echo "ERROR: protected FASTQ is tracked: ${tracked_path}" >&2
            protected_fastq_found=true
            ;;
    esac
done < <(git ls-files | rg '\.(fastq|fq)(\.gz)?$' || true)

if [[ "${protected_fastq_found}" == true ]]; then
    exit 1
fi

if git ls-files | rg -q '^resources/GRCh38/|\.part$'; then
    echo 'ERROR: downloaded reference content or a partial download is tracked' >&2
    exit 1
fi

protected_genomic_found=false
while IFS= read -r tracked_path; do
    case "$tracked_path" in
        tests/fixtures/*)
            ;;
        *)
            echo "ERROR: individual-level genomic artifact is tracked: $tracked_path" >&2
            protected_genomic_found=true
            ;;
    esac
done < <(
    printf '%s\n' "${public_files[@]}" | \
        rg '\.(bam|bai|cram|crai|g\.vcf|vcf|bcf)(\.gz|\.tbi|\.csi)?$' || true
)
if [[ "$protected_genomic_found" == true ]]; then
    exit 1
fi

git check-ignore -q --no-index private/example.fastq.gz || {
    echo 'ERROR: protected FASTQ ignore rule is missing' >&2
    exit 1
}

if git check-ignore -q --no-index tests/fixtures/fastq/reviewer_sample_1_R1.fastq.gz; then
    echo 'ERROR: synthetic reviewer FASTQ fixtures are unexpectedly ignored' >&2
    exit 1
fi

if rg -n '/mnt/(ssd|data)|/staging/lalli|nf_stage/' "${scan_files[@]}"; then
    echo 'ERROR: machine-local path found in public release content' >&2
    exit 1
fi

if rg -ni \
    'AKIA[0-9A-Z]{16}|aws_secret_access_key|BEGIN (RSA |OPENSSH |EC )?PRIVATE KEY|password[[:space:]]*[:=]' \
    "${scan_files[@]}"; then
    echo 'ERROR: possible credential or private key found in public release content' >&2
    exit 1
fi

if printf '%s\n' "${public_files[@]}" | rg -q 'consolidate_mbv_reports\.py$'; then
    echo 'ERROR: legacy MBV consolidation script must not be release content' >&2
    exit 1
fi

if rg -ni 'GRCh37|hg19|CHM13|T2T' "${scan_files[@]}"; then
    echo 'ERROR: out-of-scope reference-build text found in public release content' >&2
    exit 1
fi

echo 'Repository checks passed.'
