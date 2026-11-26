library(data.table)

broad_main <- fread("~/Documents/brainvar/data/broad_b38_contigs_no_alt_hla_decoy.txt", data.table=F, fill=TRUE)
all_gencode <- fread("~/Documents/brainvar/data/GRCh38.p13.genome.all_contigs.txt", data.table=F)

broad_main$contig_simp <- gsub("AC:", "", broad_main$V2)

gencode_overlap <- all_gencode[which(all_gencode$V2 %in% broad_main$contig_simp),]
main_gencode_contigs <- all_gencode[1:25,]

final_gencode <- rbind(main_gencode_contigs, gencode_overlap)
final_gencode$V1 <- gsub(">", "", final_gencode$V1)


write.table(as.data.frame(final_gencode$V1), "~/Documents/brainvar/data/contigs_in_broad_b38_NO_alt_hla_decoy_contigs.txt", col.names=F, row.names=F, quote=F, sep="\t")


# This output was copied to this file:
# GRCh38.p13.Broad.no_alt_no_decoy.contigs.txt


awk '/^>/ {print substr($1, 2)}' /home/ec2-user/GRCh38.p13.GTEx.no-alt-hla-decoy.genome.fa | wc -l
