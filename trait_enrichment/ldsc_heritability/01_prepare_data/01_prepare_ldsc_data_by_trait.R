setwd("/data3/Bert/BrainVar-NarjesRohani-2025-11-11")
rm(list=ls()); gc()
options(stringsAsFactors=F)

###### added QC 2026-01-16
# remove SNP with maf < 0.05 (using plink) and Chi2(GWAS) > max(80,.001max(N))
# remove MHC (chr6:28,510,120-33,480,577)

require(data.table)
require(doParallel)
require(viridis)

# function
selectData <- function(chr,trt,work.dir){
  # sumstats
  sumstats=fread(paste0("SumstatData/",trt,"/",trt,".qc.sumstats.gz"),header=T,data.table=F)
  sumstats=sumstats[sumstats$Z^2 <= max(80,0.001*max(sumstats$N)),]
  sumstats=sumstats[!duplicated(sumstats$SNP),]

  rownames(sumstats)=sumstats$SNP
  
  # select the snps for this chromosome
  snps=rownames(sumstats)[gsub(":","",substr(rownames(sumstats),1,2)) == chr]
  if(chr == 6){
    snps=snps[as.numeric(unlist(lapply(strsplit(snps,":"),`[[`,2))) < 28510120 | as.numeric(unlist(lapply(strsplit(snps,":"),`[[`,2))) > 33480577]
  }
  fwrite(data.frame(snps=snps),paste0(work.dir,"/",trt,"-chr",chr,".snp"),row.names=F,col.names=F,quote=F,sep="\t")
  
  # select the sumstats 
  sumstats=sumstats[snps,]
  sumstats=sumstats[!is.na(sumstats$SNP),]
  fwrite(sumstats,paste0(work.dir,"/",trt,"-chr",chr,".sumstats"),row.names=F,col.names=T,quote=F,sep="\t")
  
  # select the annotations
  annot=fread(paste0("AnnotationData_2026-02-11/BrainVar-annotations-chr",chr,".annot.gz"),header=T,data.table=F)
  rownames(annot)=annot$SNP
  annot=annot[snps,]
  annot=annot[!is.na(annot$SNP),]
  fwrite(annot,paste0(work.dir,"/",trt,"-chr",chr,".annot"),row.names=F,col.names=T,quote=F,sep="\t")
  
  # select the genotypes
  plink.command=paste("/data3/Software/Plink2/plink2",
                      "--bfile",paste0("GenotypeData/BrainVar-HG38-imputed-cM-chr",chr),
                      "--extract",paste0(work.dir,"/",trt,"-chr",chr,".snp"),
                      "--maf","0.05",
                      "--make-bed",
                      "--out",paste0(work.dir,"/",trt,"-chr",chr),
                      sep=" ")
  system(plink.command)
  bim=fread(paste0(work.dir,"/",trt,"-chr",chr,".bim"),header=F,data.table=F)
  return(data.frame(chr=chr,n.snp=length(snps),n.stat=nrow(sumstats),n.annot=nrow(annot),n.gnt=nrow(bim)))
}

### 
trts=list.files("SumstatData",pattern="sumstat",recursive = T,full.names = F)
trts=unique(unlist(lapply(strsplit(trts,"/",fixed=T),`[[`,1)))
chrs=1:22

RE=NULL
for(trt in trts){
  # create the directories to use
  print(trt)
  trt.dir=paste0("ldscData_2026-02-11/",trt)
  dir.create(trt.dir,recursive = T,showWarnings = F)
  work.dir=paste0("ldscWork_2026-02-11/",trt)
  dir.create(work.dir,recursive = T,showWarnings = F)
  
  # select the data for each chromosome
  re=rbindlist(mclapply(chrs,selectData,trt=trt,work.dir=work.dir,mc.cores=8,mc.preschedule=T))
  
  # combine the chromosome specific values
  #### combine the plink files
  files=paste0(work.dir,"/",trt,"-chr",chrs)
  fwrite(data.frame(files=files[-1]),paste0(work.dir,"/",trt,".gntfiles"),row.names=F,col.names=F,quote=F,sep="\t")
  plink.command=paste("/data3/Software/Plink/plink",
                      "--bfile",files[1],
                      "--merge-list",paste0(work.dir,"/",trt,".gntfiles"),
                      "--make-bed",
                      "--out",paste0(trt.dir,"/",trt),
                      sep=" ")
  system(plink.command)
  # read the bim file for reference
  bim=fread(paste0(trt.dir,"/",trt,".bim"),header=F,data.table=F)
  
  ### sumstats
  sumstats=NULL
  for(chr in chrs){
    sumstats=rbind.data.frame(sumstats,fread(paste0(work.dir,"/",trt,"-chr",chr,".sumstats"),header=T,data.table=F))
  }
  rownames(sumstats)=sumstats$SNP
  sumstats=sumstats[bim$V2,]
  fwrite(sumstats,paste0(trt.dir,"/",trt,".sumstats.gz"),row.names=F,col.names=T,quote=F,sep="\t")
  
  ### annotations
  annot=NULL
  for(chr in chrs){
    annot=rbind.data.frame(annot,fread(paste0(work.dir,"/",trt,"-chr",chr,".annot"),header=T,data.table=F))
  }
  rownames(annot)=annot$SNP
  annot=annot[bim$V2,]
  # add an intercept to the annotations
  annot$Int=1
  
  fwrite(annot,paste0(trt.dir,"/",trt,".annot.gz"),row.names=F,col.names=T,quote=F,sep="\t")
  
  # calculate the minor allele frequency
  plink.command=paste("/data3/Software/Plink/plink",
                      "--bfile",paste0(trt.dir,"/",trt),
                      "--freq",
                      "--out",paste0(trt.dir,"/",trt),
                      sep=" ")
  system(plink.command)
  
  
  RE=rbind.data.frame(RE,data.frame(trt,n.snp=nrow(bim)))
  

}
RE
WriteXLS::WriteXLS(RE,paste0("ldscData_2026-02-11","/","snp_counts.xlsx"))
