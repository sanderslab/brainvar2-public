ml tabix/0.2.6

for f in *.bed
do
bgzip $f && tabix -p bed $f.gz
done

for f in PC*
do
cd $f
gzip *
ln -s ../*.bed.gz .
ln -s ../*.bed.gz.tbi .
ln -s ../../brainvar2.022924.joint_called.recalibrated.varQC.AF.vcf.gz .
ln -s ../../brainvar2.022924.joint_called.recalibrated.varQC.AF.vcf.gz.tbi .
cd ..
done

source ~/.bashrc

ml Rmath/3.6.3
ml eigen/3.3.7
ml boost
ml fastqtl

for f in PC*
do
cd $f
for bed in *.bed.gz
do
cell=`echo $bed | sed 's/.bed.gz//'`
echo $cell
for chr in {1..22}
do
fastQTL --vcf ../../brainvar2.022924.joint_called.recalibrated.varQC.AF.vcf.gz --bed $cell.bed.gz --permute 1000 --out $cell.quantile.txt.gz --normal --cov $cell.cov.txt.gz --commands 22 $cell.commands.22.tmp
sed -n ${chr}p $cell.commands.22.tmp  > $cell.commands.chr$chr.lsf
rm $cell.commands.22.tmp
done
cat $cell.commands.chr*.lsf > $cell.commands.22.lsf
rm $cell.commands.chr*.lsf
done
echo $f
cd ..
done


for f in PC*
do
cd $f
for bed in *.bed.gz
do
cell=`echo $bed | sed 's/.bed.gz//'`
echo $cell
fastQTL --vcf ../../brainvar2.022924.joint_called.recalibrated.varQC.AF.vcf.gz --bed $bed --permute 1000 --out $cell.quantile.txt.gz --normal --cov $cell.cov.txt.gz --commands 22 $cell.commands.22.tmp
sed 's/^ //; s/$/ /' $cell.commands.22.tmp > $cell.commands.22.lsf
rm $cell.commands.22.tmp
done
cd ..
done

cp ../*.pre.command.tmp ./

for f in PC*
do
cd $f
for bed in *.bed.gz
do
cell=`echo $bed | sed 's/.bed.gz//'`
echo $cell
cat ../$cell.pre.command.tmp $cell.commands.22.lsf > $cell.commands.lsf
bsub < $cell.commands.lsf
done
cd ..
done


ml R/4.2.0
R

library(tidyverse)
library(qvalue)
library(parallel)

fastqtl_out_files = list.files('./', pattern = '.gz.chr[0-9]', full.names = T, recursive = TRUE)
d <- tibble(files=fastqtl_out_files,
            cell_type=basename(files) %>% 
              gsub('.quantile.txt.gz..+','',.) %>% 
              gsub('\\.',' ',.) %>% 
              gsub('OPCs   COPs','OPCs / COPs',.)) %>% 
  mutate(PCs=dirname(files) %>% basename(.)) %>% 
  mutate(file_content=map(files,read.table,header=F)) %>% 
  unnest(file_content) %>% 
  filter(!is.na(V11)) %>% 
  select(-files) %>% 
  setNames(c("cell_type","PCs","pid", "nvar", "shape1", "shape2", "dummy", "sid", "dist", "npval", "slope","ppval", "bpval")) %>% 
  group_by(cell_type,PCs) %>% 
  mutate(adj_p=qvalue(bpval)$qvalues) %>% 
  ungroup() %>% 
  arrange(cell_type,bpval)


dir.create('output/',showWarnings = FALSE,recursive = TRUE)
write_tsv(d,'output/eqtl.allPCs.txt')

selected_pcs <- d %>% 
  filter(adj_p<0.05) %>% 
  dplyr::count(PCs) %>% filter(n==max(n))

d %>% 
  filter(PCs==selected_pcs$PCs) %>% 
  dplyr::select(-PCs) %>% 
  write_tsv(.,paste0('output/eqtl.',selected_pcs$PCs,'.txt'))

q()
n


# Clean up
# Remove tmp files
for d in PC*
do
cd $d
rm *quantile*
cd ..
done

####
# Nominal Pass
source ~/.bashrc
cd $eqtl_path


ml Rmath/3.6.3
ml eigen/3.3.7
ml boost
ml fastqtl

mkdir -p 'PC10_nominal'
cd PC10_nominal
ln -s ../PC10/*.bed.gz .
ln -s ../PC10/*.bed.gz.tbi .
ln -s ../PC10/*.cov.txt.gz .
ln -s ../PC10/brainvar2.022924.joint_called.recalibrated.varQC.AF.vcf.gz .
ln -s ../PC10/brainvar2.022924.joint_called.recalibrated.varQC.AF.vcf.gz.tbi .

# bulk 
mkdir -p 'PC50_nominal'
cd PC50_nominal
ln -s ../PC50/*.bed.gz .
ln -s ../PC50/*.bed.gz.tbi .
ln -s ../PC50/*.cov.txt.gz .
ln -s ../PC50/brainvar2.022924.joint_called.recalibrated.varQC.AF.vcf.gz .
ln -s ../PC50/brainvar2.022924.joint_called.recalibrated.varQC.AF.vcf.gz.tbi .


for bed in *.bed.gz
do
cell=`echo $bed | sed 's/.bed.gz//'`
fastQTL --vcf brainvar2.022924.joint_called.recalibrated.varQC.AF.vcf.gz --bed $bed --out $cell.quantile.txt.gz --normal --cov $cell.cov.txt.gz --commands 22 $cell.commands.22.tmp
sed 's/^ //; s/$/ /' $cell.commands.22.tmp > $cell.commands.22.lsf
rm $cell.commands.22.tmp
done


for bed in *.bed.gz
do
cell=`echo $bed | sed 's/.bed.gz//'`
echo $cell
cat ../$cell.pre.command.tmp $cell.commands.22.lsf > $cell.commands.lsf
bsub < $cell.commands.lsf
done

# Clean data

cd $eqtl_path'/PC10_nominals'
rm *quantile.txt.gz
rm *.bed.gz
rm *.bed.gz.tbi  
rm *.cov.txt.gz
rm *.commands.lsf
rm *.vcf.gz
rm *.vcf.gz.tbi

for f in *.quantile.txt.gz.*
do
gzip $f
x=`echo $f | sed 's/:.\+$/.gz/' | sed 's/.quantile.txt.gz//'`
mv $f.gz $x
done

cd /PC5_indep
ml tabix/0.2.6

for f in *.qtltools.bed
do
bgzip $f && tabix -p bed $f.gz
done


ln -s ../../new.vcf.gz .
ln -s ../../new.vcf.gz.tbi .
ln -s ../PC5/*.cov.txt.gz  .
