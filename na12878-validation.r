require(sqldf)
options(stringsAsFactors=FALSE)

alone = read.table("alleles_alone.txt")
joint = read.table("alleles_joint.txt")

colnames(alone) = c("chr","pos","ref","alt")
colnames(joint) = c("chr","pos","ref","alt")

# how many "alone" variants cannot be found in the "joint" dataset?
sqldf("
select   count(*)
from     alone a
where    not exists (
	     select   null
	     from     joint j
	     where    j.chr = a.chr
	     and      j.pos = a.pos
	     and      j.ref = a.ref
	     and      j.alt = a.alt
	     )
;")
#   count(*)
# 1      135

# get the minimal representation for one allele
minrep_onerow = function(pos, ref, alt) {
	ref = unlist(strsplit(ref,split="")) # string to vector
	alt = unlist(strsplit(alt,split="")) # string to vector
	if (length(ref) == 1 && length(alt) == 1) {
		# do nothing
	} else {
		# strip off identical suffixes
		while(alt[length(alt)] == ref[length(ref)] && min(length(alt),length(ref)) > 1) {
			alt = alt[1:length(alt)-1]
			ref = ref[1:length(ref)-1]
		}
		# strip off identical prefixes and increment position
		while(alt[1] == ref[1] && min(length(alt),length(ref)) > 1) {
			alt = alt[2:length(alt)]
			ref = ref[2:length(ref)]
			pos = pos + 1
		}
	}
	ref = paste(ref,collapse="")
	alt = paste(alt,collapse="")
	return (list(pos, ref, alt))
}

# vectorized version of above
minrep_vectorized = function(pos, ref, alt) {
	stopifnot(length(pos) == length(ref) && length(ref) == length(alt))
	for (i in 1:length(pos)) {
		list_of_results = minrep_onerow(pos[i],ref[i],alt[i])
		pos[i] = list_of_results[[1]]
		ref[i] = list_of_results[[2]]
		alt[i] = list_of_results[[3]]
		# if (i %% 1000 == 0) {
		# 	print(i)
		# }
	}
	return (data.frame(pos, ref, alt))
}

# convert "alone" and "joint" allele sets to minimal representation
alone_mr = alone
alone_mr[,c("pos","ref","alt")] = minrep_vectorized(alone_mr$pos,alone_mr$ref,alone_mr$alt)

joint_mr = joint
joint_mr[,c("pos","ref","alt")] = minrep_vectorized(joint_mr$pos,joint_mr$ref,joint_mr$alt)

# how many "alone" variants are missing from "joint" once both are in minrep?
sqldf("
select   count(*)
from     alone_mr a
where    not exists (
	     select   null
	     from     joint_mr j
	     where    j.chr = a.chr
	     and      j.pos = a.pos
	     and      j.ref = a.ref
	     and      j.alt = a.alt
	     )
;")
#   count(*)
# 1        4

# only 4 - which are they?
sqldf("
select   *
from     alone_mr a
where    not exists (
	     select   null
	     from     joint_mr j
	     where    j.chr = a.chr
	     and      j.pos = a.pos
	     and      j.ref = a.ref
	     and      j.alt = a.alt
	     )
;")
#   chr       pos ref                                         alt
# 1   1 152681680   C                         CAGCTCTGGGGGCTGCTGT
# 2  12  53207583   C CCACCAAAGCCACCAGTGCCGAAACCAGCTCCGAAGCCGCCGG
# 3  17  26708547  CT                                           C
# 4  19  18392236   C                                           A

# manually looked these up and they are simply not called.
# to check without tabixing, use grep -P "1\t1526816" etc.
# Two are called slightly differently:
# 1       152681684       .       T       TCTGGGGGCTGCTGTAGCC
# 12      53207585        .       A       ACCAAAGCCACCAGTGCCGAAACCAGCTCCGAAGCCGCCGGCG
# And two are not called at all - for proof, here are the nearby variants in the joint VCF:
# 17      26708544        .       C       A       
# 17      26708547        .       C       A       
# 17      26708550        rs11658194      T       C
# 19      18392112        .       C       T       
# 19      18392255        .       G       C       
# 19      18392283        .       G       C
