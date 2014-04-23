import sqlite3
import pandas as pd
import pandas.io.sql as pd_sql
import minimal_representation as mr

# NA12878 called alone
alone = pd.read_table("alleles_alone.txt",sep=" ",
	names=['chr','pos','ref','alt'],
	dtype={'chr':str,'pos':int,'ref':str,'alt':str})

# NA12878 joint-called with 86K samples
joint = pd.read_table("alleles_joint.txt",sep=" ",
	names=['chr','pos','ref','alt'],
	dtype={'chr':str,'pos':int,'ref':str,'alt':str})

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
