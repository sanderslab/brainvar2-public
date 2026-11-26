#!/usr/bin/env bash
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1

project_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$project_dir"

tmpdir=$(mktemp -d "${TMPDIR:-/tmp}/brainvar2-mbv.XXXXXX")
trap 'rm -rf "$tmpdir"' EXIT

PYTHONPYCACHEPREFIX="$tmpdir/pycache" python3 -m py_compile \
    bin/collect_mbv_reports.py notebooks/plot_mbv_results.py
python3 -m unittest discover -s tests -p 'test_*.py'
bash -n bin/run_mbv.sh
bash bin/run_mbv.sh --help >/dev/null

python3 bin/collect_mbv_reports.py \
    --reports tests/fixtures/reports/RNA_A.mbv.txt tests/fixtures/reports/RNA_B.mbv.txt \
    --output "$tmpdir/mbv_all.tsv.gz"
test "$(gzip -cd "$tmpdir/mbv_all.tsv.gz" | wc -l)" -eq 7

python3 - <<'PY'
import json
import hashlib
from pathlib import Path

notebook = json.loads(Path("notebooks/plot_mbv_results.ipynb").read_text())
for number, cell in enumerate(notebook["cells"], start=1):
    if cell["cell_type"] == "code":
        compile("".join(cell["source"]), f"notebook-cell-{number}", "exec")

decision_pdf = Path(
    "docs/decision_records/Joseph_022124 Brainvar2 Sample Matching - documenting.pdf"
)
pdf_bytes = decision_pdf.read_bytes()
assert pdf_bytes.startswith(b"%PDF-"), "decision record is not a PDF"
assert len(pdf_bytes) == 10051768, "decision PDF size differs from supplied artifact"
assert hashlib.sha256(pdf_bytes).hexdigest() == (
    "da0324aa3535d801beeb353c2d39fcb9f2d04c08eddcf4a6410f1679bdfe674a"
), "decision PDF checksum differs from supplied artifact"
PY

mapfile -t public_files < <(
    find . -path './.git' -prune -o -path '*/__pycache__' -prune -o \
        -path './results*' -prune -o -path './_internal' -prune -o \
        -type f ! -path './scripts/check_release.sh' -print | sort
)

if printf '%s\n' "${public_files[@]}" | grep -q 'consolidate_mbv_reports.py$'; then
    echo 'ERROR: legacy consolidator found in release content' >&2
    exit 1
fi

if grep -En '/mnt/(ssd|data)|/staging/lalli|AKIA[0-9A-Z]{16}|aws_secret_access_key|BEGIN (RSA |OPENSSH |EC )?PRIVATE KEY' "${public_files[@]}"; then
    echo 'ERROR: local path, credential, or private key found in release content' >&2
    exit 1
fi

if find . -path './.git' -prune -o -path './_internal' -prune -o -type f \
    \( -name '*.bam' -o -name '*.bai' -o -name '*.cram' -o -name '*.crai' \
       -o -name '*.vcf' -o -name '*.vcf.gz' -o -name '*.bcf' \) -print -quit | grep -q .; then
    echo 'ERROR: protected genomic file found in release content' >&2
    exit 1
fi

echo 'Release checks passed.'
