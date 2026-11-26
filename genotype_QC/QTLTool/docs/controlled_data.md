# Controlled-data boundary

Do not commit or publicly upload:

- RNA BAM/SAM files and indexes;
- WGS VCF/BCF files and indexes;
- BAM lists containing project paths;
- native MBV reports and logs;
- collected candidate-level tables;
- working PowerPoint presentations and unapproved reviewer notes;
- generated figures outside the curated public decision-record directory.

The tracked test reports contain synthetic identifiers and values. Public code
contains no production paths, object-storage credentials, or raw genomic data.

The project has approved public release of non-synthetic rendered figures used
for pairing decisions and the supplied master decision PDF. Approved artifacts
must be placed under `docs/decision_records/`, where their provenance and scope
are documented. This approval does not extend to the source MBV reports,
candidate-level tables, sample metadata, BAMs, VCFs, or unrestricted review
working files.
