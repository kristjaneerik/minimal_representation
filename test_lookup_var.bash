# example script to test lookup_var.py

na12878_alone_all=/humgen/atgu1/fs03/eminikel/045atgu/minrep/NISTIntegratedCalls_13datasets_130719_allcall_UGHapMerge_HetHomVarPASS_VQSRv2.17_all_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs.minimal.vcf.gz
na12878_joint_all=/seq/dax/macarthur_joint_calling/v1/macarthur_joint_calling.site_only_plus_NA12878.vcf.gz
full_86k_vcf=#***real file location deliberately obscured***
reftable_86k=/humgen/atgu1/fs03/eminikel/045atgu/minrep/86k_gvcfs_and_bams.txt

# example 1: not using reftable. 
lookup_var.py 1 1827835 CTTGG CTGG $na12878_alone_all -i

# further examples using reftable against 86k dataset

# three ways to find the same variant
lookup_var.py 1 861223 AG CG $full_86k_vcf -t $reftable_86k -i # as found in reference VCF
lookup_var.py 1 861222 TAG TCG $full_86k_vcf -t $reftable_86k -i # a different non-minimal representation of same variant
lookup_var.py 1 861223 A C $full_86k_vcf -t $reftable_86k -i # searching for the minimal representation
# search for the other allele at the same site
lookup_var.py 1 861223 AG A $full_86k_vcf -t $reftable_86k -i

# three ways to find another variant
lookup_var.py 1 878677 CCT TCT $full_86k_vcf -t $reftable_86k -i # as found in reference VCF
lookup_var.py 1 878676 ACCT ATCT $full_86k_vcf -t $reftable_86k -i
lookup_var.py 1 878677 C T $full_86k_vcf -t $reftable_86k -i # searching for minimal representation

# and another
lookup_var.py 1 879058 T TC $full_86k_vcf -i

# also works on bi-allelic sites
lookup_var.py 1 878680 G T $full_86k_vcf -i

# test multi-variant mode
echo -e "1 861223 AG CG\n1 878677 CCT TCT" > variant_list.txt
lookup_var.py $full_86k_vcf -l variant_list.txt -t $reftable_86k -i
# or without returning individual-level variants
lookup_var.py $full_86k_vcf -l variant_list.txt -t $reftable_86k