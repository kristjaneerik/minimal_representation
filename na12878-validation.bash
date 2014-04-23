cd /humgen/atgu1/fs03/eminikel/045atgu/minrep

na12878_alone_all=/humgen/gsa-hpprojects/NA12878Collection/NIST/v2.17_09112013/NISTIntegratedCalls_13datasets_130719_allcall_UGHapMerge_HetHomVarPASS_VQSRv2.17_all_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs.minimal.vcf
na12878_joint_all=/seq/dax/macarthur_joint_calling/v1/macarthur_joint_calling.site_only_plus_NA12878.vcf.gz

# "joint" has many more multi-allelic rows
cat $na12878_alone | \grep -v ^# | \grep ',' | wc -l
# 9886
zcat $na12878_joint | \grep -v ^# | \grep ',' | wc -l
# 1005151

# find a couple of examples
cat $na12878_alone | \grep -v ^# | \grep ',' | head -1
# 1	1827835	.	CTT	CT,C	9402	PASS	set=Intersection	GT	2|1
# multi-allelic rows in "alone" are those where NA12878 is het alt, like above
zcat $na12878_joint | \grep -v ^# | \grep ',' | head -10
# here's a good example:
# 1	861223	.	AG	A,CG	2044.43	VQSRTrancheINDEL94.00to95.00	AC=2,4;AF=1.177e-05,2.354e-05;AN=169918

# need to subset both to just the Broad exome interval list.
exome=/humgen/gsa-hpprojects/GATK/bundle/current/b37/Broad.human.exome.b37.interval_list
# 1       30366   30503   +       target_1
# 1       69089   70010   +       target_2
# 1       367657  368599  +       target_3
# 1       621094  622036  +       target_4
# 1       861320  861395  +       target_5
# 1       865533  865718  +       target_6

java -Xmx8g -jar /humgen/gsa-hpprojects/GATK/bin/current/GenomeAnalysisTK.jar \
   -R /humgen/gsa-hpprojects/GATK/bundle/current/b37/human_g1k_v37.fasta \
   -T SelectVariants \
   --variant $na12878_alone_all \
   -L $exome \
   -o na12878_alone_exome.vcf

java -Xmx8g -jar /humgen/gsa-hpprojects/GATK/bin/current/GenomeAnalysisTK.jar \
   -R /humgen/gsa-hpprojects/GATK/bundle/current/b37/human_g1k_v37.fasta \
   -T SelectVariants \
   --variant $na12878_joint_all \
   -L $exome \
   -o na12878_joint_exome.vcf

na12878_alone=na12878_alone_exome.vcf
na12878_joint=na12878_joint_exome.vcf

# how many sites?
cat $na12878_joint | \grep -v ^# | wc -l
# 6386149
cat $na12878_alone | \grep -v ^# | wc -l
# 14980

# how many multi-allelic rows?
cat $na12878_joint | \grep -v ^# | awk '{print $1,$2,$4,$5}' | awk '$4 ~ /,/ {print $0}' | wc -l
# 566440
cat $na12878_alone | \grep -v ^# | awk '{print $1,$2,$4,$5}' | awk '$4 ~ /,/ {print $0}' | wc -l
# 5

# decompose each allele into its own row and print only chr, pos, ref, alt to two files
cat $na12878_joint | \grep -v ^# | awk '{print $1,$2,$4,$5}' | awk -F'[ ,]' '{for(i=4;i<=NF;i++) {print $1,$2,$3,$i} }' > alleles_joint.txt
cat $na12878_alone | \grep -v ^# | awk '{print $1,$2,$4,$5}' | awk -F'[ ,]' '{for(i=4;i<=NF;i++) {print $1,$2,$3,$i} }' > alleles_alone.txt

wc -l alleles_joint.txt
# 6996872 alleles_joint.txt
wc -l alleles_alone.txt
# 14985 alleles_alone.txt

# 15k exonic variants probably reasonable given the strictness of this callset - see readme
# ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/variant_calls/NIST/README.NIST.v2.17.txt