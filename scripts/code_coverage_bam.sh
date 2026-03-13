REGION="Pf3D7_13_v3:1725000-1728000"

while read sample
do
    samtools depth -r $REGION ${sample}.bqsr.bam | \
    awk -v s=$sample '{sum+=$3; n++} END {if(n>0) print s, sum/n; else print s, "NA"}'
done < /mnt/storage13/ahri/training_FA5/admixture/samples_er.txt \
> /mnt/storage13/ahri/training_FA5/coverage_comaprision/depth_k13_pipeline.txt