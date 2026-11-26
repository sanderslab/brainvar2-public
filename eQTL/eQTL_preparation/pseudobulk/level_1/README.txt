Level 1 RNA pseudobulk utilities
================================

The Python script in this directory reads a period-specific level 1 H5AD,
joins cell-type and library annotations, and uses `decoupler.get_pseudobulk`
to aggregate expression by library and cell type. The current script sets
`mode="mean"`; to produce sum-based data, use an explicitly reviewed version
with `mode="sum"`. `combine_pseudobulk.R` combines per-period mean files into a
matrix with gene rows and sample/cell-type/period columns.

## Inputs and outputs

The scripts currently contain fixed paths for the H5AD input, cell metadata,
and output directory. Update those paths for the execution environment. The
Python script assumes an H5AD with `Barcode` values in `obs` and a metadata CSV
with `cell_id`, `celltype`, and `library` columns. The combine script expects
files named like `BrainVar2_period_<number>_level_1_pseudobulk_mean.txt` in its
configured working directory.

## Environment and reproducibility

The historical environment notes record Python 3.9.16, anndata 0.10.3,
decoupler 1.5.0, pandas 2.1.3, scanpy 1.9.6, and NumPy 1.26.2. These versions
are not pinned by an environment file in this directory. The historical run
notes also report a local modification to decoupler's pseudobulk utility to
retain zero-count samples/features. That modification is not included here;
the exact installed package source must therefore be documented to reproduce
that behavior. Do not edit an installed package without preserving its source
and recording the change.

H5AD data, metadata, and generated pseudobulk matrices are not included.
Confirm aggregation mode, zero-count handling, sample metadata joins, and
output dimensions before downstream analysis.
