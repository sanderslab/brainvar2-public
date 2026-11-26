setwd("/data3/Bert/BrainVar-NarjesRohani-2025-11-11")
rm(list=ls()); gc()
options(stringsAsFactors=F)

require(data.table)
require(doParallel)
require(viridis)

### 
trts=list.files("SumstatData",pattern="sumstat",recursive = T,full.names = F)
trts=unique(unlist(lapply(strsplit(trts,"/",fixed=T),`[[`,1)))
chrs=1:22

master.scr.file="ldscData_2026-02-11/h2.master.scr"
write.table("#!/bin/bash",file=master.scr.file,row.names=F,col.names=F,quote=F,sep="\t",append=F)
write.table("",file=master.scr.file,row.names=F,col.names=F,quote=F,sep="\t",append=T)

# read the prevalence data
KP=as.data.frame(read_excel("/data3/DownLoadedData/BrainVar-NarjesRohani-2025-11-07/traits-2026-01-01/LDSC.parameters.xlsx",na="NA"))
rownames(KP)=KP$Pheno


for(trt in trts){
  trt.dir=paste0("ldscData_2026-02-11/",trt)
  work.dir=paste0("ldscWork_2026-02-11/",trt)
  
  
  scr.file=paste0(trt.dir,"/",trt,"-h2.scr")
  write.table("#!/bin/bash",file=scr.file,row.names=F,col.names=F,quote=F,sep="\t",append=F)
  write.table("",file=scr.file,row.names=F,col.names=F,quote=F,sep="\t",append=T)
  
  # and the heritability
  if(is.na(KP[trt,"K"])){
    # continuous trait
    h2.command=paste("/data3/Software/ldsc/ldsc.py",
                     "--h2",paste0(trt.dir,"/",trt,".sumstats.gz"),
                     "--ref-ld",paste0(trt.dir,"/",trt,"_ref_ld"),
                     "--w-ld",paste0(trt.dir,"/",trt,"_w_ld"),
                     "--out",paste0(trt.dir,"/",trt,"_h2"),
                     sep=" ")
  }else{
    # binary trait
    h2.command=paste("/data3/Software/ldsc/ldsc.py",
                     "--h2",paste0(trt.dir,"/",trt,".sumstats.gz"),
                     "--ref-ld",paste0(trt.dir,"/",trt,"_ref_ld"),
                     "--w-ld",paste0(trt.dir,"/",trt,"_w_ld"),
                     "--samp-prev",KP[trt,"P",],
                     "--pop-prev",KP[trt,"K",],
                     "--out",paste0(trt.dir,"/",trt,"_h2"),
                     sep=" ")
  }
  h2.command
  write.table("# h2",file=scr.file,row.names=F,col.names=F,quote=F,sep="\t",append=T)
  write.table(h2.command,file=scr.file,row.names=F,col.names=F,quote=F,sep="\t",append=T)
  write.table("",file=scr.file,row.names=F,col.names=F,quote=F,sep="\t",append=T)
  system(paste("chmod +x",scr.file))
  
  # add to the master script file
  write.table(paste(scr.file),file=master.scr.file,row.names=F,col.names=F,quote=F,sep="\t",append=T)
  write.table("",file=master.scr.file,row.names=F,col.names=F,quote=F,sep="\t",append=T)
  
}

system(paste("chmod +x",master.scr.file))




