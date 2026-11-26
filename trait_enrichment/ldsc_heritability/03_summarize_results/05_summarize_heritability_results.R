setwd("/data3/Bert/BrainVar-NarjesRohani-2025-11-11")
rm(list=ls()); gc()
options(stringsAsFactors = F)

require(data.table)


trts=list.files("SumstatData",pattern="sumstat",recursive = T,full.names = F)
trts=unique(unlist(lapply(strsplit(trts,"/",fixed=T),`[[`,1)))

# get the categories
cats=fread(paste0("ldscData_2026-02-11/",trts[1],"/",trts[1],".annot.gz"),header=F,data.table=F,sep="\t",nrows=1)[,-c(1:4)]
cats=gsub("_all","",cats)

# get the counts for each classification
M=NULL
for(trt in trts){
  tmp=fread(paste0("ldscData_2026-02-11/",trt,"/",trt,"_ref_ld.l2.M_5_50"),header=F,data.table=F)
  M=rbind.data.frame(M,data.frame(tmp))
}
rownames(M)=trts
colnames(M)=cats


#### nor read the output
re=NULL
h2.est=NULL
se.est=NULL
prop.h2.est=NULL
enrich.est=NULL
coef.est=NULL
coef.se.est=NULL

for(trt in trts){
  in.df=fread(paste0("ldscData_2026-02-11/",trt,"/",trt,"_h2.log"),header=F,data.table=F,sep="\t")
  tmp=in.df[grep("scale h2:",in.df[,1])[1],]
  scale=tolower(unlist(strsplit(tmp," ",fixed=T))[2])
  h2=as.numeric(unlist(strsplit(unlist(strsplit(tmp,": "))[2]," (",fixed=T))[1])
  se=as.numeric(gsub(")","",unlist(strsplit(tmp," (",fixed=T))[2],fixed=T))
  # categories
  tmp=in.df[grep("Categories:",in.df[,1])[1],]
  cat.list=unlist(strsplit(unlist(strsplit(tmp,": "))[2]," "))
  cat.list=gsub("_allL2_0","",cat.list)
  # heritabilities
  tmp=in.df[grep("scale h2:",in.df[,1])[2],]
  h2.list=unlist(strsplit(unlist(strsplit(tmp,": "))[2]," "))
  h2.list=as.numeric(h2.list[h2.list != ""])
  names(h2.list)=paste(cat.list,"h2",sep=".")
  # standard error
  tmp=in.df[grep("scale h2 SE:",in.df[,1])[1],]
  se.list=unlist(strsplit(unlist(strsplit(tmp,": "))[2]," "))
  se.list=as.numeric(se.list[se.list != ""])
  names(se.list)=paste(cat.list,"se",sep=".")

  # proportion
  tmp=in.df[grep("Proportion",in.df[,1])[1],]
  prop.h2.list=unlist(strsplit(unlist(strsplit(tmp,": "))[2]," "))
  prop.h2.list=as.numeric(prop.h2.list[prop.h2.list != ""])
  names(prop.h2.list)=paste(cat.list,"ph2",sep=".")
  
  
  
  # enrichment
  tmp=in.df[grep("Enrichment:",in.df[,1])[1],]
  enrich.list=unlist(strsplit(unlist(strsplit(tmp,": "))[2]," "))
  enrich.list=as.numeric(enrich.list[enrich.list != ""])
  names(enrich.list)=paste(cat.list,"enr",sep=".")
  
  # coefficient
  tmp=in.df[grep("Coefficients:",in.df[,1])[1],]
  coef.list=unlist(strsplit(unlist(strsplit(tmp,": "))[2]," "))
  coef.list=as.numeric(enrich.list[coef.list != ""])
  names(coef.list)=paste(cat.list,"coef",sep=".")
  
  # standard error
  tmp=in.df[grep("Coefficient SE",in.df[,1])[1],]
  coef.se.list=unlist(strsplit(unlist(strsplit(tmp,": "))[2]," "))
  coef.se.list=as.numeric(coef.se.list[se.list != ""])
  names(coef.se.list)=paste(cat.list,"se",sep=".")
  
  
  re=rbind.data.frame(re,data.frame(trt=trt,scale=scale,h2=h2,se=se)) #,t(h2.list),t(se.list)))
  h2.est=rbind.data.frame(h2.est,data.frame(t(h2.list)))
  se.est=rbind.data.frame(se.est,data.frame(t(se.list)))
  prop.h2.est=rbind.data.frame(prop.h2.est,data.frame(t(prop.h2.list)))
  enrich.est=rbind.data.frame(enrich.est,data.frame(t(enrich.list)))
  coef.est=rbind.data.frame(coef.est,data.frame(t(coef.list)))
  coef.se.est=rbind.data.frame(coef.se.est,data.frame(t(coef.se.list)))

}
rownames(h2.est)=trts; colnames(h2.est)=cats
rownames(se.est)=trts; colnames(se.est)=cats
rownames(prop.h2.est)=trts; colnames(prop.h2.est)=cats
rownames(enrich.est)=trts; colnames(enrich.est)=cats
rownames(coef.est)=trts; colnames(coef.est)=cats
rownames(coef.se.est)=trts; colnames(coef.se.est)=cats

# fix the OCD issue
total.h2=cbind.data.frame(re,data.frame(Z=round(re$h2/re$se,4),N=M$Int))
rownames(total.h2)=total.h2$trt
total.h2=total.h2[rownames(total.h2) != "OCD",]; rownames(total.h2)=gsub("OCD2025","OCD",rownames(total.h2))
total.h2$trt=rownames(total.h2)

h2.est=h2.est[rownames(h2.est) != "OCD",]; rownames(h2.est)=gsub("OCD2025","OCD",rownames(h2.est))
se.est=se.est[rownames(se.est) != "OCD",]; rownames(se.est)=gsub("OCD2025","OCD",rownames(se.est))
prop.h2.est=prop.h2.est[rownames(prop.h2.est) != "OCD",]; rownames(prop.h2.est)=gsub("OCD2025","OCD",rownames(prop.h2.est))
enrich.est=enrich.est[rownames(enrich.est) != "OCD",]; rownames(enrich.est)=gsub("OCD2025","OCD",rownames(enrich.est))
coef.est=coef.est[rownames(coef.est) != "OCD",]; rownames(coef.est)=gsub("OCD2025","OCD",rownames(coef.est))
coef.se.est=coef.se.est[rownames(coef.se.est) != "OCD",]; rownames(coef.est)=gsub("OCD2025","OCD",rownames(coef.est))

RE=list(total.h2=total.h2,partial.h2=round(h2.est,5),proportion.h2=round(prop.h2.est,4),enrichment=round(enrich.est,4))
WriteXLS::WriteXLS(RE,paste0("Results_2026-02-11/BrainVar-LDSC-H2-Results-",Sys.Date(),".xlsx"),row.names = T)


rownames(re)=re$trt
re$Z=re$h2/re$se
re$N=rowSums(M[rownames(re),])
h2.est=h2.est[rownames(re),]
se.est=se.est[rownames(re),]

# additional matrices
prop.snp=M[rownames(re),]/re$N
prop.h2=h2.est/rowSums(h2.est)
prop.snp=prop.snp[rownames(re),]
prop.h2=prop.h2[rownames(re),]

tmp=prop.h2
tmp[tmp < 0]=0
rowSums(tmp)

tmp=prop.h2
tmp[tmp > 0]=0
rowSums(tmp)

enrich=prop.h2/prop.snp

coeff=h2.est/M

fread("ldscData/ADHD/ADHD_h2.log",header=F,data.table=F,sep="\t",skip=26,nrows=1)
grep("Categories",input[,1])
read.table(input[28,])
