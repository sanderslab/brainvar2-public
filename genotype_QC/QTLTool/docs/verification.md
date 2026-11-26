# Verification

`scripts/check_release.sh` checks Bash and Python syntax, executes focused unit
tests, validates the notebook JSON, runs the collector on synthetic reports, and
scans release files for protected data, credentials, and the legacy consolidator.

The collector tests verify that:

- every native candidate row is retained;
- native numeric strings are not recalculated;
- report and candidate identifiers are represented explicitly;
- candidate-set mismatches are rejected by default;
- compressed output is deterministic.

The runner test verifies historical report naming and the recovered mapping
quality setting without requiring genomic software.

## Controlled archive validation

On 2026-07-15, the collector successfully read all 261 retained bulk reports.
Each report contained 271 candidate rows, yielding exactly 70,731 rows. No
source report was empty or omitted.

An independent read-only comparison established that the retained historical
primary summary agrees with the native maxima for all 255 consensus matches and
all 6 unresolved cases. The historical second- and third-candidate fields do not
represent the second- and third-highest concordance values and are outside the
new collector’s scope.

The exact 2023 QTLtools executable hash was not recovered from a production log.
QTLtools 1.3.1 is supported by manuscript wording and tracked project code; the
published repository pins that version for portable reruns.
