library(data.table)

broad_all <- fread("~/Documents/brainvar/data/broad_b38_all_contigs.txt", data.table=F, fill=TRUE)
broad_main <- fread("~/Documents/brainvar/data/broad_b38_contigs_no_alt_hla_decoy.txt", data.table=F, fill=TRUE)
#all_gencode <- fread("~/Documents/brainvar/data/all_contigs_GRCh38.p14", data.table=F)
all_gencode <- fread("~/Documents/brainvar/data/GRCh38.p13.genome.all_contigs.txt", data.table=F)

# 639 contigs in p.13 and 706 in p.14 

# Extract ALT, HLA and decoy contigs from Broad assembly
broad_rm <- broad_all[-which(broad_all$V1 %in% broad_main$V1),]

broad_rm$contig_simp <- gsub("AC:", "", broad_rm$V2)
broad_all$contig_simp <- gsub("AC:", "", broad_all$V2)
broad_main$contig_simp <- gsub("AC:", "", broad_main$V2)

sum(all_gencode$V2 %in% broad_all$contig_simp) # 169 + 25 canonical contigs = 194 (this is true for both p.13 and p.14)
sum(all_gencode$V2 %in% broad_rm$contig_simp) # 0
sum(all_gencode$V2 %in% broad_main$contig_simp) # 169 + 25 canonical contigs = 194 (this is true for both p.13 and p.14) (this is redundant)

# Look at contigs in filtered Broad assembly not in gencode
broad_main[-which(broad_main$contig_simp %in% all_gencode$V2),]

# Look at contigs in GENCODE seemingly not in Broad (N=512 in p.14, 445 in p.13 - ignoring 25 canonical that clearly do match)
# could be a string matching thing rather than truly not being in there
gencode_unique <- all_gencode[-which(all_gencode$V2 %in% broad_all$contig_simp),]
gencode_unique$contig_simp <- gsub(">", "", gencode_unique$V1)

sum(gencode_unique$contig_simp %in% broad_all$contig_simp) # 261 recovered in broad_all (true for both p13 and p14)
table(broad_all[which(broad_all$contig_simp %in% gencode_unique$contig_simp),1]) # all labeled as alt (true for both p13 and p14)

gencode_unique <- gencode_unique[-which(gencode_unique$contig_simp %in% broad_all$contig_simp),] # grep some of the ones from gencode in broad, (184 in p13, 251 in p14) just to be sure

write.table(broad_rm, "~/Documents/brainvar/data/broad_b38_alt_hla_decoy_contigs.txt", col.names=F, row.names=F, quote=F, sep="\t")