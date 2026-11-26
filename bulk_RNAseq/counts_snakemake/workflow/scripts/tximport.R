log <- file(snakemake@log[[1]], open="wt")
sink(log ,type = "output")
sink(log, type = "message")

library(txdbmaker)
library(tximport)
library(data.table)

run_tximport <- function(data_path, out_path, myparam) {
    print(data_path)
    print(out_path)
    txdb = txdbmaker::makeTxDbFromGFF("/data/GRCh38-references/gencode.v41.annotation.gtf")
    k <- keys(txdb, keytype = "TXNAME")
    tx2gene <- select(txdb, k, "GENEID", "TXNAME")
    
    files <- file.path(data_path)
    txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

    mf <- data.frame(txi$counts)
    setDT(mf, keep.rownames = TRUE)[]
    fwrite(mf, out_path)

}

run_tximport(snakemake@input[[1]], snakemake@output[[1]], snakemake@config[["myparam"]])