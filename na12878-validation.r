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

source("minimal_representation.r")

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
