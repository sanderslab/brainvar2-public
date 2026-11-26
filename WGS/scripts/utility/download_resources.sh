#!/usr/bin/env bash
# Download and verify the public GRCh38 resources used by the BrainVar2
# production workflow.
#
# Usage:
#   bash scripts/utility/download_resources.sh
#   bash scripts/utility/download_resources.sh --verify-only
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
resource_root="${repo_root}/resources/GRCh38"
source_root='https://ngi-igenomes.s3.amazonaws.com/igenomes/Homo_sapiens/GATK/GRCh38'
verify_only=false

usage() {
    sed -n '2,/^set -euo pipefail/{ /^#/s/^# \?//p }' "$0"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --verify-only)
            verify_only=true
            shift
            ;;
        --help|-h)
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

for command_name in curl sha256sum wc tr readlink cmp; do
    if ! command -v "${command_name}" >/dev/null 2>&1; then
        echo "ERROR: required command not found: ${command_name}" >&2
        exit 1
    fi
done

verify_file() {
    local path="$1"
    local expected_size="$2"
    local expected_sha256="$3"
    local actual_size
    local actual_sha256

    if [[ ! -f "${path}" ]]; then
        echo "ERROR: missing resource: ${path}" >&2
        return 1
    fi

    actual_size=$(wc -c < "${path}" | tr -d '[:space:]')
    if [[ "${actual_size}" != "${expected_size}" ]]; then
        echo "ERROR: size mismatch for ${path}" >&2
        echo "       expected ${expected_size}; observed ${actual_size}" >&2
        return 1
    fi

    actual_sha256=$(sha256sum "${path}")
    actual_sha256=${actual_sha256%% *}
    if [[ "${actual_sha256}" != "${expected_sha256}" ]]; then
        echo "ERROR: SHA-256 mismatch for ${path}" >&2
        echo "       expected ${expected_sha256}; observed ${actual_sha256}" >&2
        return 1
    fi
}

fetch_file() {
    local relative_path="$1"
    local expected_size="$2"
    local expected_sha256="$3"
    local destination="${resource_root}/${relative_path}"
    local partial="${destination}.part"
    local url="${source_root}/${relative_path}"

    if [[ -f "${destination}" ]]; then
        if verify_file "${destination}" "${expected_size}" "${expected_sha256}"; then
            echo "Verified: ${relative_path}"
            return 0
        fi
        echo "ERROR: refusing to replace an existing invalid resource." >&2
        echo "       Remove ${destination} explicitly, then rerun." >&2
        return 1
    fi

    mkdir -p "$(dirname "${destination}")"
    echo "Downloading: ${relative_path}"
    if ! curl \
        --fail \
        --location \
        --retry 5 \
        --retry-delay 2 \
        --show-error \
        --silent \
        --continue-at - \
        --output "${partial}" \
        "${url}"; then
        echo "ERROR: download failed: ${url}" >&2
        return 1
    fi

    if ! verify_file "${partial}" "${expected_size}" "${expected_sha256}"; then
        rm -f "${partial}"
        echo "ERROR: discarded invalid partial download for ${relative_path}" >&2
        return 1
    fi

    mv "${partial}" "${destination}"
    echo "Verified: ${relative_path}"
}

verify_manifest_entry() {
    local relative_path="$1"
    local expected_size="$2"
    local expected_sha256="$3"
    local destination="${resource_root}/${relative_path}"

    if ! verify_file "${destination}" "${expected_size}" "${expected_sha256}"; then
        return 1
    fi
    echo "Verified: ${relative_path}"
}

failures=0
while IFS=$'\t' read -r relative_path expected_size expected_sha256; do
    [[ -n "${relative_path}" ]] || continue
    if [[ "${verify_only}" == true ]]; then
        verify_manifest_entry "${relative_path}" "${expected_size}" "${expected_sha256}" ||
            failures=$((failures + 1))
    else
        fetch_file "${relative_path}" "${expected_size}" "${expected_sha256}" ||
            failures=$((failures + 1))
    fi
done <<'MANIFEST'
Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta	3249912778	93157a161863464c9435062fd67c173fdaf99cb8b32f1455018361387ffa5564
Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta.fai	160928	edefd93c489dc1baefad312f40388089f8db5cf6dcc3ba0955669ead274e8b6b
Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.dict	581712	07523b4a0afd1127ae909f92f8b3911aff096461aed7eff3f1acf5bdfc1dc103
Sequence/BWAIndex/Homo_sapiens_assembly38.fasta.64.alt	487553	d9254da07b8030e26129dc29d9d02b9c30a360b233a367ce041691c00407d510
Sequence/BWAIndex/Homo_sapiens_assembly38.fasta.64.amb	20199	4199db9e913e372b2dae4ac17723def657f696daf71024b66f6b695947918c77
Sequence/BWAIndex/Homo_sapiens_assembly38.fasta.64.ann	455474	49231cd8462826afe424aeb20b6b2b9ecff32ff1a90f88f5d5f3c1a0cfb9ba99
Sequence/BWAIndex/Homo_sapiens_assembly38.fasta.64.bwt	3217347004	a01ecc6bfae7abcf3542f68314da843d5594f8b8cf3eb7799371ac3e99231754
Sequence/BWAIndex/Homo_sapiens_assembly38.fasta.64.pac	804336731	5dd6629f95e9765a19821fd1ffc3e883e5af177527a220a6b72a4e5e492114dc
Sequence/BWAIndex/Homo_sapiens_assembly38.fasta.64.sa	1608673512	b87d29d937b58dfde9943d8300ab0f9490bae5ba4ddc54726da023aeebd48da4
Annotation/GATKBundle/dbsnp_138.hg38.vcf.gz	1560889937	04d74d585c2adbcc601e159d6c610d479598779749d7e67a26e35af6cf58a97c
Annotation/GATKBundle/dbsnp_138.hg38.vcf.gz.tbi	2321422	12c86c71249bd7a29059b213807e3681aadb3237c0ef60188fc144eabad7514f
Annotation/GATKBundle/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz	3053999	8bca5d7af481f8bb0cf657d0665a6a84be7f0b24c254bf6a912c9733e39f4217
Annotation/GATKBundle/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi	421689	593660b27c1edca1201f2b779dec189837dee9dab993baf2e7b3200f87fc038c
Annotation/GATKBundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz	20685880	e11195ef119927518f01965bfe63591ec4c761b727ed26e72f17fed23a6739e5
Annotation/GATKBundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi	1500013	6856f5ad5690e1f095960ecee450e26a67f9486a8336aa08be4f599d958c1bac
Annotation/GATKBundle/hapmap_3.3.hg38.vcf.gz	62043448	6c0402a3b77e97e1067f52fe1b13bbaf39eb4d8b846d24e76b3ea8fe82237188
Annotation/GATKBundle/hapmap_3.3.hg38.vcf.gz.tbi	1552123	2aed605c4b96aebd480d9ef515a59777387e17b6988423688b1993aadf2ea155
Annotation/GATKBundle/1000G_omni2.5.hg38.vcf.gz	53238342	9c4aa8e66972591500cc46c2b6792a74d81ae856cc467e22f5f58169dc987399
Annotation/GATKBundle/1000G_omni2.5.hg38.vcf.gz.tbi	1544108	60aaea8e7b62831671c4c129dbb1deb364e991afe185bf3a0a220916fb7e9cc3
Annotation/GATKBundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz	1888262073	734d6634a98094ada4eb365a3a6aa4babcbecf3617a95fb5a23ed3ff98a595d1
Annotation/GATKBundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi	2128536	7abef1e31a3d513f393c0bccbfb80df12ab7eb141d175707e3dee87981a2f4bc
MANIFEST

if [[ "${failures}" -gt 0 ]]; then
    echo "ERROR: ${failures} resource file(s) failed validation." >&2
    exit 1
fi

whole_genome_dir="${resource_root}/Sequence/WholeGenomeFasta"
bwa_dir="${resource_root}/Sequence/BWAIndex"

ensure_alignment_link() {
    local filename="$1"
    local source_path="${whole_genome_dir}/${filename}"
    local destination="${bwa_dir}/${filename}"
    local relative_target="../WholeGenomeFasta/${filename}"

    if [[ -L "${destination}" ]]; then
        if [[ "$(readlink "${destination}")" != "${relative_target}" ]]; then
            echo "ERROR: unexpected link target: ${destination}" >&2
            return 1
        fi
        return 0
    fi

    if [[ -e "${destination}" ]]; then
        if ! cmp -s "${source_path}" "${destination}"; then
            echo "ERROR: ${destination} differs from ${source_path}" >&2
            return 1
        fi
        return 0
    fi

    if [[ "${verify_only}" == true ]]; then
        echo "ERROR: missing Parabricks reference link: ${destination}" >&2
        return 1
    fi

    ln -s "${relative_target}" "${destination}"
}

for filename in \
    Homo_sapiens_assembly38.fasta \
    Homo_sapiens_assembly38.fasta.fai \
    Homo_sapiens_assembly38.dict; do
    ensure_alignment_link "${filename}" || failures=$((failures + 1))
done

if [[ "${failures}" -gt 0 ]]; then
    echo "ERROR: ${failures} Parabricks reference link(s) failed validation." >&2
    exit 1
fi

echo "GRCh38 resource bundle verified at ${resource_root}"
