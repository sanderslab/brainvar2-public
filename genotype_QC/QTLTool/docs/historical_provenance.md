# Historical provenance

## Recovered execution evidence

The controlled project archive retains 261 non-empty native bulk MBV reports:
46 BrainVar2 reports and 215 BrainVar reports. Every native report contains one
header and 271 candidate genotype rows. A December 2023 consolidated summary
contains the same 261 RNA queries. These artifacts demonstrate successful MBV
execution; the presence of source code alone is not the basis for that claim.

The exact executable hash used for those 2023 jobs was not recorded in a
retained run log. A tracked 2022 BrainVar RNA Nextflow module specified the
`naotokubota/qtltools:1.3.1` container and mapping-quality threshold 15, and the
manuscript records QTLtools 1.3.1. Those facts support the selected reproducible
rerun environment but do not prove that every historical report came from that
tracked module. Historical filenames differ from that module's output suffix.

## Controlled equivalence audit

On 2026-07-15, an independent read-only audit in the controlled workspace
compared all 261 native reports with the retained historical summary:

| Check | Result |
| --- | ---: |
| Native reports linked to historical rows | 261/261 |
| Unique heterozygous/homozygous consensus candidates | 255 |
| Historical consensus matches reproduced | 255/255 |
| Discordant metric winners | 6 |
| Historical unresolved statuses reproduced | 6/6 |
| Overall primary status equivalence | 261/261 |
| Legacy “second closest” equal to true second-highest | 0/261 |
| Legacy “third closest” equal to true third-highest | 0/261 |

This establishes equivalence for the primary consortium-era match/unresolved
result while independently demonstrating that the legacy secondary-ranking
fields are unsuitable for reuse.

## Legacy code boundary

The February 2024 `consolidate_mbv_reports.py` remains unchanged in controlled
storage and is not copied into this repository. It selected second and third
indices from the ascending end of an argsort, among other limitations. The new
collector replaces only its data-loading role and does not calculate rankings.
The controlled audit is provenance evidence, not release code; this repository
neither edits the legacy file nor rewrites historical summaries.

Historical native reports, summaries, notebooks, and complete sample-identifier
inventories remain in the approved controlled environment. Their machine-local
locations are recorded in the ignored internal audit, not in public Git history.

## Published decision evidence

The approved 212-page master deck under `docs/decision_records/` preserves the
MBV methods explanation, non-synthetic comparison plots, cohort-level findings,
and individual edge-case discussions used during manual review. It is direct
evidence of the adjudication process, while the native reports and candidate-level
tables remain controlled.
