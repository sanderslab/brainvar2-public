# Gene functional categories

This folder contains the GO annotation/clustering and expression-plotting
scripts for the functional categories.

## Workflow

1. `01_GO_annotation_and_clusters.R` selects the configured significant gene
   set, maps gene identifiers, queries GO annotations, summarizes terms, and
   writes the GO-cluster workbooks.
2. `02_classify_GO_function_categories.R` assigns genes to the Methods
   categories using the cellular-component GO terms: plasma membrane
   (GO:0005886), synapse (GO:0045202), and nucleus (GO:0005634). Neuronal
   Communication requires plasma membrane or synapse annotation and no
   nucleus annotation. Transcriptional Regulation requires nucleus annotation
   and neither plasma membrane nor synapse. Remaining genes, including genes
   without returned GO annotations, are assigned `None/Other`.
3. `03_plot_expression_groups.R` plots expression for the neuronal
   communication and transcriptional regulation groups across the configured
   developmental periods and cell types.
