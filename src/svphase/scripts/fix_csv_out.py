import pandas as pd
import numpy as np

def fix(in_fpath, out_fpath):
	d = pd.read_csv(in_fpath, sep='\t', index_col=False, header=0)
	d['contig'] = map(lambda x:eval(x)[0], d.ix[:,'index'])
	d['start'] = map(lambda x:eval(x)[1], d.ix[:,'index'])
	d['end'] = map(lambda x:eval(x)[2], d.ix[:,'index'])
	d = d[list(d.columns[-3:])+list(d.columns[:-3])]
	d = d.drop('index', axis=1)

	d.to_csv(out_fpath, sep='\t', index=False, float_format='%0.5f', header=True)
if __name__ == '__main__':
	import argparse

	parser = argparse.ArgumentParser(description="Fixing pandas sorting issue")
	parser.add_argument('csv_fpath', help='Input pandas csv to fix')
	parser.add_argument('out_csv_fpath', help='Output pandas csv file with "contig, start, end" prefix')

	args = parser.parse_args()

	fix(args.csv_fpath, args.out_csv_fpath)  
