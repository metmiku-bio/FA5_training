bam_dir="/mnt/storage13/ahri/hrp2_analysis/bam_eritrea"
output="/mnt/storage13/ahri/training_FA5/bam"
ref="/mnt/storage13/ahri/plasmodium_falciparum/Pfalciparum.genome.fasta"

# cd "$bam_dir"

# ls *.bam | parallel -j 8 '
# sample=$(basename {} .bam)
# samtools view -b {} Pf3D7_13_v3 > '"$output"'/${sample}.bam
# echo "{} extraction finished"
# '

# it is also interesting to check by lowering the used quality to 20 but for now i am using q 30
# cd "$output"

# ls *.bam | parallel -j 10 '
# sample=$(basename {} .bam)

# samtools view -b -q 30 -f 2 {} \
# | samtools sort -o '"$output"'/${sample}.sorted.mapq30.properpaired.bam

# samtools index '"$output"'/${sample}.sorted.mapq30.properpaired.bam

# echo "{} all filtering finished"
# '

#mkdir -p variant
variant="/mnt/storage13/ahri/training_FA5/variant"
# cd $output



# parallel -j 6 '
# sample=$(basename {} .sorted.mapq30.properpaired.bam)

# gatk HaplotypeCaller \
# -R '"$ref"' \
# -I {} \
# -O '"$variant"'/${sample}.g.vcf.gz \
# -ERC GVCF \
# --native-pair-hmm-threads 6

# echo "{} variant calling is finished"
# ' ::: *.sorted.mapq30.properpaired.bam

# gatk CombineGVCFs -R $ref $(for f in ${variant}/*g.vcf.gz ; do echo --variant $f; done) -O ${variant}/combined.old_way.g.vcf.gz
# echo "combined gvcf finish"
# gatk GenotypeGVCFs -R $ref -V ${variant}/combined.old_way.g.vcf.gz -O combined.variants.genotype.gz
# echo "genotype finished "
# bcftools view -m2 -M2 -i 'type="SNP"' /mnt/storage13/ahri/training_FA5/combined.variants.genotype.gz -Oz -o ${variant}/combined.variants.snp.genotype.gz
# echo "snp filtering finish"

# bcftools index ${variant}/combined.variants.snp.genotype.gz
#gatk IndexFeatureFile \
# -I /mnt/storage13/ahri/training_FA5/variant/combined.variants.snp.genotype.gz
# gatk VariantFiltration -V ${variant}/combined.variants.snp.genotype.gz \
#     -filter "QD < 2.0" --filter-name "QD2"  \
#     -filter "QUAL < 200.0" --filter-name "QUAL200" \
#     -filter "SOR > 3.0" --filter-name "SOR3"  \
#     -filter "FS > 60.0" --filter-name "FS60"  \
#     -filter "MQ < 40.0" --filter-name "MQ40"  \
#     -O ${variant}/combined.variants.hard.filtered.genotype.vcf.gz

# echo "variant filteration finish"

#activate this environment since having high jfava
#conda activate gatk_gcnv
 java -Xmx8g -jar /mnt/storage13/ahri/snpEff/snpEff.jar Plasmodium_falciparum ${variant}/combined.variants.hard.filtered.genotype.vcf.gz > ${variant}/combined.variants.hard.filtered.genotype.ann.vcf

# for making the vcfphylip 

python vcf2phylip.py -i ../variant/combined.variants.hard.filtered.genotype.vcf.gz --fasta
# iqtree input file 
iqtree -s combined.variants.hard.filtered.genotype.min4.fasta -m MFP -B 1000


# admixture

plink2 --vcf 012026_pf_mgen_ug_india_eritrea_ethiopia_popgen.filt.csq.bi.GT.miss0.2.vqslod.filt.snps.erit.vcf.gz --make-bed --out output_er
plink2 --bfile output_er --indep-pairwise 50 5 0.5 --out output.QC.pruned --set-all-var-ids \@:# --bad-ld 
plink2 --bfile output_er --freq --out output_er
plink2 --bfile output_er --read-freq output_er.afreq --pca 5 --out output_er.pca --set-all-var-ids \@:#

awk '{$1="0" ;print $0}' output_er.bim >  output_er.bim.tmp
mv output_er.bim output_er.bim.source
mv output_er.bim.tmp output.filtered.bim

BED_FILE="/mnt/storage13/ahri/training_FA5/admixture/output_er.bed"
OUTDIR="/mnt/storage13/ahri/training_FA5/admixture"

cd $OUTDIR

for k in {3..12}; do
    admixture --cv -j8 $BED_FILE $k | tee log${k}.out
done