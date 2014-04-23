# example script to test lookup_var.py

na12878_alone_all=/humgen/atgu1/fs03/eminikel/045atgu/minrep/NISTIntegratedCalls_13datasets_130719_allcall_UGHapMerge_HetHomVarPASS_VQSRv2.17_all_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs.minimal.vcf.gz
na12878_joint_all=/seq/dax/macarthur_joint_calling/v1/macarthur_joint_calling.site_only_plus_NA12878.vcf.gz
full_86k_vcf=#***real file location deliberately obscured***
reftable_86k=/humgen/atgu1/fs03/eminikel/045atgu/minrep/86k_gvcfs_and_bams.txt

# example 1: not using reftable. note use of NA where reftable param would go
lookup_var.py $na12878_alone_all NA 1 1827835 CTTGG CTGG
# output:
# # You searched for:  1827835 CT C
# # Relevant line from VCF:  1 1827835 . CT CT,C 9402 PASS set=Intersection GT
# #SAMPLE	CALL
# NA12878	2|1

# further examples using reftable against 86k dataset

# three ways to find the same variant
lookup_var.py $full_86k_vcf $reftable_86k 1 861223 AG CG # as found in reference VCF
lookup_var.py $full_86k_vcf $reftable_86k 1 861222 TAG TCG # a different non-minimal representation of same variant
lookup_var.py $full_86k_vcf $reftable_86k 1 861223 A C # searching for the minimal representation

# three ways to find another variant
lookup_var.py $full_86k_vcf $reftable_86k 1 878677 CCT TCT # as found in reference VCF
lookup_var.py $full_86k_vcf $reftable_86k 1 878676 ACCT ATCT
lookup_var.py $full_86k_vcf $reftable_86k 1 878677 C T # searching for minimal representation

