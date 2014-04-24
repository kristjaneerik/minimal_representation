import sqlite3
import pandas as pd
import pandas.io.sql as pd_sql
from minimal_representation import *

# NA12878 called alone
alone = pd.read_table("alleles_alone.txt",sep=" ",
    names=['chr','pos','ref','alt'],
    dtype={'chr':str,'pos':int,'ref':str,'alt':str})

# NA12878 joint-called with 86K samples
joint = pd.read_table("alleles_joint.txt",sep=" ",
    names=['chr','pos','ref','alt'],
    dtype={'chr':str,'pos':int,'ref':str,'alt':str})

alone_mr = alone
for row in range(0,alone_mr.shape[0]):
    chr, pos, ref, alt = alone_mr.loc[row]
    pos, ref, alt = get_minimal_representation(pos,ref,alt)
    alone_mr.loc[row] = chr, pos, ref, alt
    if row > 100:
        break

joint_mr = joint
for row in range(0,joint_mr.shape[0]):
    chr, pos, ref, alt = joint_mr.loc[row]
    pos, ref, alt = get_minimal_representation(pos,ref,alt)
    joint_mr.loc[row] = chr, pos, ref, alt
    if row > 100:
        break

# above was too slow, try using apply instead of .loc
def mr_by_row(row):
    chr, pos, ref, alt = row
    pos, ref, alt = mr.get_minimal_representation(pos,ref,alt)
    return chr, pos, ref, alt
alone_mr = alone.apply(mr_by_row,axis=1)
# much faster but resulting df has 1 column which is a 4-tuple,
# instead of 4 columns. Googled and couldn't figure out how to
# unpack tuples in apply


cnx = sqlite3.connect(':memory:')

# can write the tables directly to sqlite3
pd_sql.write_frame(alone, name='alone', con=cnx)
pd_sql.read_sql("select * from alone limit 10;", cnx)

pd_sql.write_frame(joint, name='joint', con=cnx)
pd_sql.read_sql("select * from joint limit 10;", cnx)

# but what we want to do is convert the site-based 
# representation to allele-based representation

# create new tables to hold alleles from each dataset
pd_sql.execute("create table alone_al \
    (chr text, pos int, ref text, alt text);", cnx)
pd_sql.execute("create table joint_al \
    (chr text, pos int, ref text, alt text);", cnx)

# and another pair to hold minimal representation alleles
pd_sql.execute("create table alone_mr \
    (chr text, pos int, ref text, alt text);", cnx)
pd_sql.execute("create table joint_mr \
    (chr text, pos int, ref text, alt text);", cnx)

# for performance, let's do all the inserts and
# then commit only once when we're done.
# this requires having a cursor
c = cnx.cursor()

# walk through rows (genomic sites) splitting alleles
for site in range(0,alone.shape[0]):
    chr, pos, ref, alt = alone.loc[site]
    if (',' in alone['alt'][site]): # multi-allelic sites
        alt_alleles = alt.split(',')
        for alt_allele in alt_alleles: # loop over each allele
            # "al" table gets the allele as-is
            sil = c.execute("insert into alone_al (chr,pos,ref,alt) \
                values('%s',%s,'%s','%s')" % (chr,pos,ref,alt_allele))
            # "mr" table gets the minimal representation of the allele
            pos_mr, ref_mr, alt_mr = get_minimal_representation(pos,ref,alt_allele)
            sil = c.execute("insert into alone_al (chr,pos,ref,alt) \
                values('%s',%s,'%s','%s')" % (chr,pos_mr,ref_mr,alt_mr))
    else: # bi-allelic sites need no special treatment
        sil = c.execute("insert into alone_al (chr,pos,ref,alt) \
            values('%s',%s,'%s','%s')" % (chr,pos,ref,alt))
        sil = c.execute("insert into alone_mr (chr,pos,ref,alt) \
            values('%s',%s,'%s','%s')" % (chr,pos,ref,alt))
    if (site > 10000):
        break
    if (site % 1000 == 0):
        print site

cnx.commit()


