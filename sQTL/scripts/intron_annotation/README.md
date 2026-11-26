# Intron Annotation Utilities

This folder contains scripts for annotating splicing junctions using reference annotations (e.g., GENCODE GTF). The goal is to determine whether detected introns are known, cryptic, or novel, and to assign gene-level information for downstream interpretation and visualization.

---

## 🧬 Overview of Annotation Steps

### Step 1: **Prepare Reference-Based Annotations**

* **Script**: `prepare_leafviz_annot.sh`
* **Input**:

  * GENCODE GTF file
  * Path to LeafCutter directory
* **Output**:

  * `<prefix>_all_introns.bed.gz`
  * `<prefix>_threeprime.bed.gz`
  * `<prefix>_fiveprime.bed.gz`
  * `<prefix>_all_exons.txt.gz`

---

### Step 2: **Annotate Introns**

* **Scripts**: `annotate_intron.sh` (wrapper), `annotate_intron.R` (core)
* **Input**:

  * Annotation prefix (from step 1)
  * LeafCutter intron count file (`leafcutter_perind.counts.gz`)
* **Output**:

  * `leafviz.RData` — R object containing cluster-level verdict, genes, transcripts, and constitutive scores
* **Details**:

  * Labels each intron with a "verdict":

    * `annotated`: matches a known intron
    * `novel annotated pair`: both splice sites are annotated but not as a pair
    * `cryptic_fiveprime` / `cryptic_threeprime`: one known splice site
    * `cryptic_unanchored`: both sites unknown
    * `unknown_strand` / `cryptic`: partial info
  * Assigns gene name, Ensembl ID, transcript overlap
  * Computes a "constitutive score" for each intron

---

### Step 3: **Summarize Annotation Table**

* **Script**: `summarise_annot_v1.0.R`
* **Input**:

  * `leafviz.RData` (from Step 2)
* **Output**:

  * `all.introns.RData`: `data.frame` with merged verdict, gene, transcript, constitutive score, coordinates

---

## 📥 Input Requirements

* GENCODE-compatible GTF file
* `leafcutter_perind.counts.gz` file from junction clustering
* LeafCutter scripts (for Perl helper)
* R packages: `data.table`, `dplyr`, `argparser`, `parallel`, `stringr`, `magrittr`

