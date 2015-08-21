import pandas as pd
import numpy as np
import os
from collections import defaultdict
from operator import itemgetter
from svphase.utils.common import sort_chr
from svphase.utils.config import PRECISION
from svphase.learn.evaluation import ClassLabel

class ConfusionMatrix():
	def __init__(self):
		pass

	def latex_cmat(self, truth, pred):
		pred = pd.read_csv(pred, sep='\t', index_col=None, header=0).rename(columns={'best':'label1'}).set_index(['contig','start','end'])
		truth = pd.read_csv(truth, sep='\t', index_col=None, header=None, names=['contig','start','end','label3']).set_index(['contig','start','end'])
		
		cols = ClassLabel().classes + ['no_data',]
	
		df = pd.concat([pred,truth],axis=1).reset_index()
		g = df.groupby(['label1','label3'])
		cmat_s = g.apply(lambda df:df.index.size)
		total =float( cmat_s.sum()) 
		cmat_df = cmat_s.reset_index()
		print cmat_s, cmat_df
		total -= cmat_df.loc[cmat_df['label1']=='no_data',0].sum()

		cmat = pd.DataFrame(np.zeros((len(cols), len(cols))), index=pd.Index(cols, name='label3'), columns=pd.Index(cols, name='label1'))

		for (c1,c3), cell_count in cmat_s.iteritems():
			cmat.loc[c3,c1] = cell_count
		print cmat, total

		print "\\hline \\hline"
		print "{0} & {1} \\\\".format(cmat.columns.name," & ".join(cmat.columns))
		print "\\hline"
		print "{0} & {1} \\\\".format(cmat.index.name," &"*(len(cmat.columns)-1))
		for idx, s in cmat.iterrows():
			print "{0} & {1} \\\\".format(idx, " & ".join(["{0:d} ({1:.2f})".format(int(s[c]),s[c]/total) for c in cols]))
		 
if __name__ == '__main__':
	import sys
	import argparse
	
	parser = argparse.ArgumentParser(description='Creates a confusion matrix')
	parser.add_argument('truth', help='Truth bed file with annotations')
	parser.add_argument('pred', help='Prediction csv file with best column')
	
	args = parser.parse_args()
	
	w = ConfusionMatrix()
	w.latex_cmat(args.truth, args.pred)

