# test - 3 samples
nohup parallel --joblog bam2junc_joblog.txt -j3 '
  SAMPLE=$(basename {})
  OUTDIR=/home/ec2-user/outputs/junctions
  OUTFILE=${OUTDIR}/${SAMPLE}.junc
  bash /home/ec2-user/leafcutter/scripts/bam2junc.sh {} ${OUTFILE}
' :::: /home/ec2-user/outputs/bam_list_junc.txt > parallel_bam2junc.log 2>&1 &


# all samples
ls -1 ~/data/AlignedBam1/*.bam > bam_list_all_samples.txt

nohup bash -c '
export PATH=/home/ec2-user/leafcutter/scripts:$PATH
parallel --joblog bam2junc_joblog.txt -j8 "
  SAMPLE=\$(basename {})
  OUTDIR=/home/ec2-user/outputs/junctions
  OUTFILE=\${OUTDIR}/\${SAMPLE}.junc
  bash /home/ec2-user/leafcutter/scripts/bam2junc.sh {} \${OUTFILE}
" :::: /home/ec2-user/outputs/bam_list_all_samples.txt
echo "[INFO] All junction files written to juncfiles.txt"
' > parallel_bam2junc.log 2>&1 &