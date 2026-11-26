# Cohort accounting

BrainVar2 WGS counts refer to different observational units and should not be
used interchangeably.

| Cohort frame | Count | Definition |
| --- | ---: | --- |
| Joint-calling WGS records | 271 | Unique WGS `LibraryID` values; all are present in the final joint manifest |
| WGS records linked to a donor | 263 | WGS library rows with a non-empty `SubjectID` |
| Directly linked unique donors | 254 | Unique non-empty `SubjectID` values among WGS rows |
| Unassigned WGS records | 8 | WGS rows without a `SubjectID` |
| Additional resequencing records | 9 | Records beyond one WGS library per linked donor |
| Putative donor-level WGS entities | 262 | 254 linked donors plus 8 unassigned records, only if each unassigned record represents a distinct donor |
| Usable WGS library records | 234 | WGS rows with library-level `Usable=Yes` |
| Usable unique donors in production metadata v1.3 | 232 | Subject rows with `Usable_DNA=Yes` in the analysis metadata release |
| Usable unique donors in corrected metadata v1.3.1 | 233 | One donor-level usability flag changed after v1.3 |

The complete library accounting is:

```text
271 WGS records = 254 linked unique donors
                + 8 unassigned WGS records
                + 9 additional records from resequenced donors
```

Seven linked donors have two WGS library records and one linked donor has
three. The metadata do not assign `SubjectID` values to the eight unassigned
records, so the 262 value is strongly supported as a working donor-level count
but is not a directly counted unique-ID total.

The 271-sample joint manifest has 6,775 rows: 271 unique WGS library IDs by 25
contigs. Its sample set is identical to the 271 WGS `LibraryID` values in both
metadata v1.3 and corrected metadata v1.3.1.

## Metadata versions

Metadata v1.3 is the production analysis freeze. The v1.3.1 correction
contains the same 842 library IDs and the same 271 WGS library IDs, but one WGS
row's non-empty
`SubjectID` was reassigned. Aggregate WGS library, linked-record, unique-donor,
unassigned-record, and resequencing counts are unchanged. The v1.3.1
subject table also changes one `Usable_DNA` flag from `No` to `Yes`, producing
233 usable donors instead of the 232 in v1.3.

The unreleased draft library metadata v1.4 has 841 rows rather than 842 because
one bulk-RNA row was removed. Its 271 WGS rows are identical to v1.3; it
does not carry the v1.3.1 WGS `SubjectID` reassignment. Consequently,
v1.4 cannot explain the 271-versus-262 WGS discrepancy. Version differences
must be considered for bulk-RNA and eQTL sample selection separately.

## Reproducing the aggregate counts

Run the controlled metadata and manifest through:

```bash
python3 bin/reconcile_cohort_counts.py \
  --manifest PRIVATE_JOINT_MANIFEST.csv \
  --library-metadata PRIVATE_TEAM_LIBRARY_METADATA_v1.3.tsv \
  --subject-metadata PRIVATE_TEAM_SUBJECT_METADATA_v1.3.tsv \
  --draft-library-metadata PRIVATE_DRAFT_LIBRARY_METADATA_v1.4.tsv \
  --output cohort_counts.json
```

The JSON includes exact input paths, sizes, modification times, aggregate
counts, matrix-completeness checks, and metadata/manifest set-overlap counts.
It does not emit library or donor identifiers. Input metadata remain protected
and must not be committed.
