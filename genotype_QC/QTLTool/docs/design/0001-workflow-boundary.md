# Decision 0001: preserve the historical analysis shape

- Status: Accepted
- Date: 2026-07-15

## Decision

The public repository implements the same four-stage analysis used for the
project: Bash MBV loop, Python concatenation, notebook plotting, and manual
PowerPoint review. It does not incorporate MBV into the WGS pipeline and does
not automate final pairing decisions.

The legacy consolidator remains unchanged in controlled storage. A new collector
replaces its data-loading role only: it concatenates every native candidate row
and computes no nearest-neighbor or separation fields.
