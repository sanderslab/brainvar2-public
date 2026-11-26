# Manual review

The notebook generates evidence for review; it does not make a pairing decision.
For the historical analysis, the per-alignment figures were manually organized
in a PowerPoint presentation and the corresponding WGS donor was identified by
visual inspection.

For a reproducible rerun:

1. Run every notebook cell against the collected candidate-level table.
2. Place the PNG or SVG for each RNA alignment into the controlled review deck.
3. Record the selected WGS donor, or mark the alignment unresolved, in the deck.
4. Have a second team member review the figure-to-decision transcription.
5. Preserve the completed deck with the run outputs in controlled storage.
6. Copy only team-approved rendered figures or decision records into
   `docs/decision_records/` for public release.

No historical numerical decision threshold was documented. The public repository
therefore must not present a software-generated best candidate as the consortium
pairing. Raw inputs and unrestricted review files remain controlled, but approved
rendered evidence and decision documentation may be committed under
`docs/decision_records/`.

The current approved master record is the 212-page
[BrainVar2 sample-matching deck](decision_records/Joseph_022124%20Brainvar2%20Sample%20Matching%20-%20documenting.pdf).
It contains the methods overview, real comparison plots, cohort-level results,
and individual-case narratives used by the working group.
