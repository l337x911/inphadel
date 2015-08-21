import pandas as pd
import numpy as np
import os
from collections import defaultdict
from operator import itemgetter
from svphase.utils.common import sort_chr
from svphase.utils.config import PRECISION


class AssignParentToScaffold():
	def __init__(self,df):
		self.df = df
		self.contig_assignments = {}
	def _assign_parent_to_scaffold(self, df):
		# Expects pA,pB in column label1, and P,M in column label2
		counts = df.groupby(['label1','label2'])

		collapse = { ('pA','P','pB','M') : 0, ('pA','M','pB','P'):0}

		for k in collapse.keys():
			if tuple(k[:2]) in counts.groups:
				collapse[k] += len(counts.groups[tuple(k[:2])])
			if tuple(k[-2:]) in counts.groups:
				collapse[k] += len(counts.groups[tuple(k[-2:])])
		return collapse	

	def set_assignment(self):
		#self.df = self.df[np.logical_and(self.df['label1']!='hom', self.df['label2']!='H')]
		self.df['label3'] = 'no_data'
		by_contig = self.df.groupby('contig1')
		for contig, cdf in by_contig:
			collapse = self._assign_parent_to_scaffold(cdf)
			collapse_list = sorted(collapse.iteritems(),key=itemgetter(1),reverse=True)

			assignment = collapse_list[0][0]
			assignment_dict = dict(((i,j) for i,j in zip(assignment[1::2],assignment[:-1:2])))
			assignment_dict['H'] = 'hom'
			self.contig_assignments[contig] = assignment_dict
			self.df.loc[cdf.index,'label3'] = cdf['label2'].apply(assignment_dict.get)

	def assign_bed_fpath(self, bed_fpath):
		assert len(self.contig_assignments)>0
		df = pd.read_csv(bed_fpath, sep='\t', index_col=None, header=None, names=['contig','start','end','parent'])
		df['truth'] = '-'
		for contig, cdf in df.groupby('contig'):
			df.loc[cdf.index,'truth'] = cdf['parent'].apply(self.contig_assignments[contig].get)
		return df

	def latex_cmat(self, original):
		cols = ['pA','pB','hom']
		b1 = pd.read_csv(original[0], sep='\t', index_col=None, header=None, names=['contig1','start1','end1','label1'])
		b2 = pd.read_csv(original[1], sep='\t', index_col=None, header=None, names=['contig2','start2','end2','label3'])
		b1_counts = b1.groupby('label1').apply(lambda df:df.index.size)[cols]	
		b2_counts = b2.groupby('label3').apply(lambda df:df.index.size)[cols]	

		
		g = self.df.groupby(['label1','label3'])
		cmat_s = g.apply(lambda df:df.index.size)
		total =float( cmat_s.sum())
		cmat_df = cmat_s.reset_index()

		cmat = pd.DataFrame(np.zeros((len(cols)+1, len(cols)+1)), index=pd.Index(cols+['disjoint',], name='label3'), columns=pd.Index(cols+['disjoint',], name='label1'))
		cmat.loc['disjoint',cols] = b1_counts-cmat_df.groupby('label1').sum()[0]
		cmat.loc[cols,'disjoint'] = b2_counts-cmat_df.groupby('label3').sum()[0]

		for (c1,c3), cell_count in cmat_s.iteritems():
			cmat.loc[c3,c1] = cell_count
		print cmat, total

		print "\\hline \\hline"
		print "{0} & {1} \\\\".format(cmat.columns.name," & ".join(cmat.columns))
		print "\\hline"
		print "{0} & {1} \\\\".format(cmat.index.name," &"*(len(cmat.columns)-1))
		for idx, s in cmat.iterrows():
			if idx=='disjoint':
				print "{0} & {1} & \\\\".format(idx, " & ".join(map("{0:d}".format, s[cols].astype(int))))
			else:
				print "{0} & {1} & {2:d} \\\\".format(idx, " & ".join(["{0:d} ({1:.2f})".format(int(s[c]),s[c]/total) for c in cols]), s['disjoint'].astype(int))
		
if __name__ == '__main__':
	import sys
	import argparse
	
	parser = argparse.ArgumentParser(description='Creates a confusion matrix')
	parser.add_argument('bed', help='Bed intersect file, (output from intersectBed -wao)')
	parser.add_argument('--assign-bed', dest='assign', default=None, nargs=1, help='Assign bed to a new label')
	parser.add_argument('--origin-bed', dest='original', default=None, nargs=2, help='Outputs confusion matrix')
	
	args = parser.parse_args()
	
	w = AssignParentToScaffold(pd.read_csv(args.bed, sep='\t', header=None, index_col=None, names=['contig1','start1','end1','label1','contig2','start2','end2','label2','overlap']))
	w.set_assignment()

	if args.assign is None:
		w.df['match'] = w.df['label1']==w.df['label3']
		if args.original is None:
			w.df[['contig1','start1','start2','label1','label2','label3','match']].to_csv(sys.stdout, sep='\t', index=None, header=None)
		else:
			w.latex_cmat(args.original)

	else:
		df = w.assign_bed_fpath(args.assign[0])
		df.to_csv(sys.stdout, sep='\t', index=None, header=None)
