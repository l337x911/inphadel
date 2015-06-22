"Merges bed files based on chr, start, end" 

import pandas as pd
import sys

def ucsc_fmt(s):
	return "{0}:{1}-{2}".format(s.ix[0],s.ix[1],s.ix[2])

def merge(beda,chain_ab,bedb):
	
	beda_df = pd.read_csv(beda, sep='\t', index_col=None, header=None, names=['chrA','startA','endA','mid'])
	chain_df = pd.read_csv(chain_ab, sep='\t', index_col=None, header=None, names=['chrB','startB','endB','mid'])
	bedb_df = pd.read_csv(bedb, sep='\t', index_col=None, header=None, names=['chrB','startB','endB','label'])
	
	df = pd.merge(chain_df, bedb_df, on=['chrB','startB','endB'], how='outer')
	df = pd.merge(beda_df, df, on=['mid'], how='inner')	
	
	df = df.apply(lambda s:pd.Series([ucsc_fmt(s[['chrA','startA','endA']]),ucsc_fmt(s[['chrB','startB','endB']]),s['label'],s['label']], index=['posA','posB','label','edit']), axis=1)
	
	df.to_csv(sys.stdout, sep='\t',index=False, header=True) 

if __name__=='__main__':
	import argparse
	
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('hg18_bed', help='Bed')
	parser.add_argument('coord_to_mid_chain', help='Coord Bed')
	parser.add_argument('hg19_bed_anno', help='Annotated bed')
	
	args = parser.parse_args()
	merge(args.hg18_bed, args.coord_to_mid_chain, args.hg19_bed_anno)
