setwd("/data3/Bert/BrainVar-NarjesRohani-2025-11-11")
rm(list=ls()); gc()
options(stringsAsFactors=F)

### run "source activate ldsc" if not yet started


require(data.table)
require(doParallel)
require(viridis)

### 
trts=list.files("SumstatData",pattern="sumstat",recursive = T,full.names = F)
trts=unique(unlist(lapply(strsplit(trts,"/",fixed=T),`[[`,1)))
chrs=1:22

master.scr.file="ldscData_2026-02-11/ldsc.master.scr"
write.table("#!/bin/bash",file=master.scr.file,row.names=F,col.names=F,quote=F,sep="\t",append=F)
write.table("",file=master.scr.file,row.names=F,col.names=F,quote=F,sep="\t",append=T)


# check the sumstat file for abs(Z) > 7
re=NULL
for(trt in trts){
  print(trt)
  tmp=fread(paste0("ldscData_2026-02-11/",trt,"/",trt,".sumstats.gz"),header=T,data.table=F)
  re=rbind.data.frame(re,data.frame(trt,n.snp=nrow(tmp),n.gt.7=sum(abs(tmp$Z)> 7,n.sign=sum(abs(tmp$Z) > abs(qnorm(5e-8/2))))))
}



for(trt in trts){
  trt.dir=paste0("ldscData_2026-02-11/",trt)
  work.dir=paste0("ldscWork_2026-02-11/",trt)
  
  
  scr.file=paste0(trt.dir,"/",trt,"-ldsc.scr")
  write.table("#!/bin/bash",file=scr.file,row.names=F,col.names=F,quote=F,sep="\t",append=F)
  write.table("",file=scr.file,row.names=F,col.names=F,quote=F,sep="\t",append=T)
  
  # w_ld
  ldsc.command=paste("/data3/Software/ldsc/ldsc.py",
                     "--l2",
                     "--bfile",paste0(trt.dir,"/",trt),
                     "--ld-wind-cm","1",
                     "--out",paste0(trt.dir,"/",trt,"_w_ld"),
                     sep=" ")
  #ldsc.command
  write.table("# w_ld",file=scr.file,row.names=F,col.names=F,quote=F,sep="\t",append=T)
  write.table(ldsc.command,file=scr.file,row.names=F,col.names=F,quote=F,sep="\t",append=T)
  write.table("",file=scr.file,row.names=F,col.names=F,quote=F,sep="\t",append=T)
  
  # partioned
  ldsc.command=paste("/data3/Software/ldsc/ldsc.py",
                     "--l2",
                     "--bfile",paste0(trt.dir,"/",trt),
                     "--ld-wind-cm","1",
                     "--annot",paste0(trt.dir,"/",trt,".annot.gz"),
                     "--out",paste0(trt.dir,"/",trt,"_ref_ld"),
                     sep=" ")
  ldsc.command
  write.table("# ref_ld",file=scr.file,row.names=F,col.names=F,quote=F,sep="\t",append=T)
  write.table(ldsc.command,file=scr.file,row.names=F,col.names=F,quote=F,sep="\t",append=T)
  write.table("",file=scr.file,row.names=F,col.names=F,quote=F,sep="\t",append=T)
  
  system(paste("chmod +x",scr.file))
  
  # add to the master script file
  write.table(paste(scr.file),file=master.scr.file,row.names=F,col.names=F,quote=F,sep="\t",append=T)
  write.table("",file=master.scr.file,row.names=F,col.names=F,quote=F,sep="\t",append=T)
  
}

system(paste("chmod +x",master.scr.file))


# run the ldsc script on the console


