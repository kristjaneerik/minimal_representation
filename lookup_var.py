#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/python_2.7.1-sqlite3-rtrees/bin/python

# look up a variant of interest in a reference VCF
# usage: python lookup_var.py reference.vcf reftable.txt chr pos ref alt

# prerequisites in path: tabix, awk, grep

import sys
import re
import os.path
import subprocess
import minimal_representation as mr
import gzip
import pysam

# indices of relevant columns in reftable, in 1-based numbering
# (1-based because these will go to awk, not python
sample_name_colno = 5 # NAME_IN_VCF column
project_name_colno = 3 # PROJECT_OR_COHORT_NAME column

# courtesy of Konrad
def get_vcf_colnames(refvcf):
	with gzip.open(refvcf) as f:
		for line in f:
			if line.startswith('#CHROM'):
				column_names = line.strip().split()
				break
	return column_names

# gets lines in a VCF that are within +- pos_buffer of chr:pos
# returns a list of strings (lines) i.e. [line1, line2, ...]
def get_vcf_lines(refvcf,pos_buffer,chr,pos):
	startpos = int(pos) - pos_buffer
	endpos = int(pos) + pos_buffer
	tabixfile = pysam.Tabixfile(refvcf)
	vcfline_generator = tabixfile.fetch(chr,startpos,endpos)
	lines = list(vcfline_generator)
	return lines

def find_var_indivs(refvcf,reftable,chr,pos,ref,alt):
	# dictionary to hold info on people with the variant allele
	variant_indivs = {}
	# convert input variants to minimal representation
	pos, ref, alt = mr.get_minimal_representation(pos,ref,alt)
	print "# Minimal representation of your search: ",
	print pos, ref, alt
	# use tabix to grab 100 bp on either side of putative variant
	lines = get_vcf_lines(refvcf,100,chr,pos)
	# get the #CHROM line from the gzipped VCF
	column_names = get_vcf_colnames(refvcf)
	# now search the lines for the variant of interest
	match_found = False # default is you haven't found a matching variant
	for line in lines:
		cols = line.split("\t")
		if len(cols) <= 9: # skip any extra non-VCF lines that appear in output
			continue
		vchr, vpos, vid, vref, valt, vqual, vfilter, vinfo, vformat = cols[:9]
		vpos = int(vpos) # must cast to into to match incoming pos variable
		valt_alleles = valt.split(",")
		for valt_allele in valt_alleles:
			vpos_mr, vref_mr, valt_allele_mr = mr.get_minimal_representation(vpos, vref, valt_allele)
			# check if we've found a match
			if vchr == chr and vpos_mr == pos and vref_mr == ref and valt_allele_mr == alt:
				match_found = True
				# output the variant info as called in the reference VCF
				print "# Relevant line from VCF: ",
				print '\t'.join(cols[:9])
				allele_no = valt_alleles.index(valt_allele)
				format_fields = vformat.split(":")
				gt_idx = format_fields.index("GT") # in what order does genotype appear
				for column_no in range(9,len(cols)):
					call = cols[column_no]
					call_fields = call.split(":")
					genotype = call_fields[gt_idx]
					alleles = re.split("/|\|",genotype) # split on / or | to support UG or HC calls.
					# check if this individual has the allele in question
					if str(allele_no) in alleles:
						sample_name = column_names[column_no]
						call_info = call
						# store this person and their call in the dict
						variant_indivs[sample_name] = call_info
				break # stop looking for more matching alleles
		if match_found:
			break # stop looking for more matching sites
	if not match_found:
		print "# No matches found."
	else:
		# at this point we have printed the variant call, and stored the info
		# on each individual with the variant.
		# now if possible we also want to look up what study they're from.
		# the reftable parameter is optional, so we'll check if the table exists
		if os.path.isfile(reftable):
			print "#SAMPLE\tPROJECT\tCALL"
			for sample_name, call_info in variant_indivs.iteritems():
				awk_command = "awk '$"+str(sample_name_colno)+" == \"" + \
					sample_name + "\" {printf $" + str(project_name_colno) + \
					"}' " + reftable
				# print awk_command
				awk_output = subprocess.check_output(awk_command,
					stderr=subprocess.STDOUT,shell=True)
				project_name = awk_output
				print sample_name+"\t"+project_name+"\t"+call_info
		else:
			print "#SAMPLE\tCALL"
			for sample_name, call_info in variant_indivs.iteritems():
				print sample_name+"\t"+call_info

if __name__ == '__main__':
	if len(sys.argv) < 7:
		print "usage: python lookup_var.py reference.vcf.gz reftable.txt chr pos ref alt"
		exit()
	else:
		refvcf, reftable, chr, pos, ref, alt = sys.argv[1:7]
		pos = int(pos) 
		refidx = refvcf + ".tbi"
		if not os.path.isfile(refidx):
			print "looks like your vcf isn't tabix'ed. couldn't find: "+refidx
			exit()
		find_var_indivs(refvcf,reftable,chr,pos,ref,alt)

