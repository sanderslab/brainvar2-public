# mkdir bulk_rna/

# goofys --region us-east-1 brainvar2:bulk_RNA_counts/ bulk_rna/

# R

local_libs = "/home/ec2-user/R_libraries"
.libPaths(local_libs)

library(data.table)

bv1 <- fread("~/bulk_rna/brainvar1/merged_count_matrix.txt", data.table=F)
bv2 <- fread("~/bulk_rna/brainvar2/merged_count_matrix.txt", data.table=F)

all <- merge(bv1, bv2, by="gene")

write.table(all, "~/all_bams_matrix.txt", col.names=T, row.names=F, quote=F, sep="\t")
