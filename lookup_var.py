#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/python_2.7.1-sqlite3-rtrees/bin/python

# look up a variant of interest in a reference VCF
# single-variant usage: python lookup_var.py chr pos ref alt reference.vcf.gz [-t reftable.txt]
# multi-variant usage:  python lookup_var.py reference.vcf.gz -l variant_list.txt [-t reftable.txt]

import sys
import re
import os.path
import minimal_representation as mr
import gzip
import pysam
import argparse

# indices of relevant columns in reftable, in 0-based numbering
sample_name_colno = 4 # NAME_IN_VCF column
project_name_colno = 2 # PROJECT_OR_COHORT_NAME column

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
    startpos = int(pos) - int(pos_buffer)
    endpos = int(pos) + int(pos_buffer)
    tabixfile = pysam.Tabixfile(refvcf)
    vcfline_generator = tabixfile.fetch(chr,startpos,endpos)
    lines = list(vcfline_generator)
    return lines

def get_project_name(reftable,sample_name):
    project_name = None
    open_function = gzip.open if reftable.endswith('.gz') else open
    with open_function(reftable) as f:
        for line in f.readlines():
            cols = line.strip().split("\t")
            if cols[sample_name_colno] == sample_name:
                project_name = cols[project_name_colno]
                break
    return project_name

def find_var_indivs(refvcf,reftable,chr,pos,ref,alt,find_indivs):
    # dictionary to hold info on people with the variant allele
    variant_indivs = {}
    # convert input variants to minimal representation
    pos, ref, alt = mr.get_minimal_representation(pos,ref,alt)
    print "##Minimal representation of your search: ",
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
                print "##Relevant line from VCF: ",
                print '\t'.join(cols[:9])
                if find_indivs:
                    allele_no = valt_alleles.index(valt_allele) + 1 # first alt allele is 1 (ref is 0)
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
        print "##No matches found."
    else:
        # at this point we have printed the variant call, and stored the info
        # on each individual with the variant.
        # now if possible we also want to look up what study they're from.
        # the reftable parameter is optional, so we'll check if the table exists
        if find_indivs and reftable is not None and os.path.isfile(reftable):
            print "#SAMPLE\tPROJECT\tCALL"
            for sample_name, call_info in variant_indivs.iteritems():
                project_name = get_project_name(reftable,sample_name)
                if project_name is None:
                    project_name = ""
                print "%s\t%s\t%s" % (sample_name, project_name, call_info)
        elif find_indivs:
            print "#SAMPLE\tCALL"
            for sample_name, call_info in variant_indivs.iteritems():
                print "%s\t%s" % (sample_name, call_info)
        # else do nothing

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Look up a variant in a reference VCF.')
    parser.add_argument('chr', metavar='chr', type=str,  nargs='?',
                       help='chromosome of variant of interest')
    parser.add_argument('pos', metavar='pos', type=int,  nargs='?',
                       help='position of variant of interest')
    parser.add_argument('ref', metavar='ref', type=str,  nargs='?',
                       help='ref allele of variant of interest')
    parser.add_argument('alt', metavar='alt', type=str,  nargs='?',
                       help='alt allele of variant of interest')
    parser.add_argument('refvcf', metavar='reference.vcf.gz', type=str, 
                       help='path to gzipped, tabixed reference VCF')
    parser.add_argument('-l','--varlist', metavar='variant_list.txt',
                       type=str, nargs='?',
                       help='path to whitespace-delimited txt file of variants\n \
                             with 4 columns: chr pos ref alt')
    parser.add_argument('-t','--reftable', metavar='reftable.txt', 
                       dest='reftable', type=str, nargs='?',
                       help='path to table about samples in VCF')
    parser.add_argument('-i','--indivs', action='store_true',
                       help='Return list of individuals with the variant')
    args = parser.parse_args()
    # check that the reference VCF is tabixed
    refidx = args.refvcf + ".tbi"
    if not os.path.isfile(refidx):
        print "looks like your vcf isn't tabix'ed. couldn't find: "+refidx
        exit()
    # check whether to run in single-variant mode or multi-variant mode
    if args.varlist is not None:
        if os.path.isfile(args.varlist):
            open_function = gzip.open if args.varlist.endswith('.gz') else open
            with open_function(args.varlist) as f:
                for line in f.readlines():
                    if len(line.strip()) < 1:
                        continue # skip empty lines
                    chr, pos, ref, alt = line.strip().split()
                    pos = int(pos)
                    find_var_indivs(args.refvcf,args.reftable,chr,pos,ref,alt,args.indivs)
        else:
            print "##Couldn't find file you specified: %s" % (args.varlist)
            exit()
    elif (not any (x is None for x in [args.chr,args.pos,args.ref,args.alt])):
        find_var_indivs(args.refvcf,args.reftable,args.chr,args.pos,args.ref,args.alt,args.indivs)
    else:
        print "##You specified neither a single variant nor a list of variants."
        exit()

