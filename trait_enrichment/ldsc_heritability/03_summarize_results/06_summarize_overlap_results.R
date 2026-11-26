setwd("/data3/Bert/BrainVar-NarjesRohani-2025-11-11")
rm(list=ls()); gc()
options(stringsAsFactors = F)

require(data.table)
require(dendextend)


trts=list.files("SumstatData",pattern="sumstat",recursive = T,full.names = F)
trts=unique(unlist(lapply(strsplit(trts,"/",fixed=T),`[[`,1)))

# get the categories
cats=fread(paste0("ldscData_2026-02-11/",trts[1],"/",trts[1],".annot.gz"),header=F,data.table=F,sep="\t",nrows=1)[,-c(1:4)]
cats=gsub("_all","",cats)

proportion.h2=NULL
enrichment=NULL
coef.Z=NULL
rank.Z=NULL
RESULTS=list()
for(trt in trts){
  tmp=fread(paste0("ldscData_2026-02-11/",trt,"/",trt,"_overlap_h2.results"),header=T,data.table=F)
  proportion.h2=rbind.data.frame(proportion.h2,data.frame(t(round(tmp$Prop._h2,4))))
  enrichment=rbind.data.frame(enrichment,data.frame(t(round(tmp$Enrichment,4))))
  coef.Z=rbind.data.frame(coef.Z,data.frame(t(round(tmp$`Coefficient_z-score`,4))))
  rank.Z=rbind.data.frame(rank.Z,data.frame(t(rank(-tmp$`Coefficient_z-score`))))
  tmp$Category=gsub("_allL2_0","",tmp$Category)
  tmp=tmp[1:12,]
  RESULTS[[trt]]=tmp
}
# fix the OCD2025 issue
RESULTS=RESULTS[names(RESULTS) != "OCD"]
names(RESULTS)=gsub("OCD2025","OCD",names(RESULTS))
WriteXLS::WriteXLS(RESULTS,paste0("Results_2026-02-11/BrainVar-LDSC-H2-Output-",Sys.Date(),".xlsx"))

colnames(proportion.h2)=cats; rownames(proportion.h2)=trts
colnames(enrichment)=cats; rownames(enrichment)=trts
colnames(coef.Z)=cats; rownames(coef.Z)=trts
colnames(rank.Z)=cats; rownames(rank.Z)=trts

# fix OCD2025
proportion.h2=proportion.h2[rownames(proportion.h2) != "OCD",]; rownames(proportion.h2)=gsub("OCD2025","OCD",rownames(proportion.h2))
enrichment=enrichment[rownames(enrichment) != "OCD",]; rownames(enrichment)=gsub("OCD2025","OCD",rownames(enrichment))
coef.Z=coef.Z[rownames(coef.Z) != "OCD",]; rownames(coef.Z)=gsub("OCD2025","OCD",rownames(coef.Z))
rank.Z=rank.Z[rownames(rank.Z) != "OCD",]; rownames(rank.Z)=gsub("OCD2025","OCD",rownames(rank.Z))

h2=as.data.frame(read_excel("Results_2026-02-11/BrainVar-LDSC-H2-Results-2026-02-16.xlsx",sheet=1))
rownames(h2)=h2$trt; h2=h2[,-1]


WriteXLS::WriteXLS(list(proportion.h2=data.frame(round(h2[,c("h2","Z")],4),proportion.h2[,1:(ncol(proportion.h2)-1)]),
                        enrichment=data.frame(round(h2[,c("h2","Z")],4),enrichment[,1:(ncol(proportion.h2)-1)]),
                        Z.regr.coef=data.frame(h2[,c("h2","Z")],coef.Z[,1:(ncol(proportion.h2)-1)])),
                        paste0("Results_2026-02-11/BrainVar-LDSC-H2-Selected-Results-",Sys.Date(),".xlsx"),row.names = T)


df=data.frame(round(h2[,c("h2","Z")],4),proportion.h2[,1:12])
df=df[df$Z > 6,]
df=df[,3:ncol(df)]*df$h2
H=hclust(dist(df))
pdf(paste0("Results_2026-02-11/BrainVar-Heritability-per-cell-type-plot-",nrow(df),"-selected-phenotypes-",Sys.Date(),".pdf"),height=8,width=8)
plot(H, main = "Heritability",xlab="Selected Phenotypes",sub="")
dev.off()

df=data.frame(h2[,c("h2","Z")],coef.Z[,1:(ncol(proportion.h2)-1)])
df=df[df$Z > 6,]
#df=df[,3:ncol(df)]*df$h2
H=hclust(dist(df),method="average")
pdf(paste0("Results_2026-02-11/BrainVar-Z-score-per-cell-type-plot-",nrow(df),"-selected-phenotypes-clust-method-AVERAGE-",Sys.Date(),".pdf"),height=8,width=8)
plot(H, main = "Z(regr.coef.)",xlab="Selected Phenotypes",sub="")
dev.off()
