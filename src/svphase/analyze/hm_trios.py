""" Assign chromosomal haplotypes A,B to parent with maximum shared alleles
"""

import sys
import numpy as na
from operator import itemgetter
import pandas as pd
from itertools import chain

pd.set_option("display.max_rows",10000)
pd.set_option("display.max_columns",10000)

class SNPTransmit(object):
	def __init__(self, pat_id, mat_id, hapmap_fpath_fmt, vcf_fpath_fmt):
		self.current_chrom = None

		#self.hapmap_fpath = '/media/T02/data/hapmap/phased/hapmap3_r2_b36_fwd.consensus.qc.poly.%s_ceu.phased'%chrom
		#self.vcf_fpath = '/media/T02/data/hic/hic_phasing_vcf/%s_step2_inference'%chrom

		self.hapmap_fpath_fmt = hapmap_fpath_fmt
		self.vcf_fpath_fmt = vcf_fpath_fmt

		self.hapmap_fpath = None
		self.vcf_fpath = None
		
		self.pat_id = pat_id
		self.mat_id = mat_id	
		self.transmit_t = None
		self.transmit_s = None

	def get_hets(self):
		df = pd.read_csv(self.hapmap_fpath, sep=' ')
		df['chr'] = self.current_chrom
		df = df.set_index(['chr','phys_position'])
		return df["{0}_A|{0}_B|{1}_A|{1}_B".format(self.pat_id,self.mat_id).split('|')],	df['rsID'] 
	
	def vcf_phased(self):
		df = pd.read_csv(self.vcf_fpath, sep='\t', names=['chr', 'phys_position', 'rsID', 'ref', 'alt', 'qual', 'filter', 'info', 'format', 'sample'])
		df = df.set_index(['chr', 'phys_position'])
		#print df['sample'] 
		phased_idx =	df['sample'].map(lambda x: x[1]=='|')
		df = df[phased_idx]
		ref = df['ref'].values
		alt = df['alt'].values
	
		a = df['sample'].map(lambda x:x[0]).astype(int)
		b = df['sample'].map(lambda x:x[2]).astype(int)
	
		df['A'] = '.'
		df['B'] = '.'
		df['A'][a.values==0] = ref[a.values==0]
		df['A'][a.values==1] = alt[a.values==1]
		df['B'][b.values==0] = ref[b.values==0]
		df['B'][b.values==1] = alt[b.values==1]
	
		return df[['A','B']], df['rsID'] 
	
	def snp_transmit(self, chrom='chr20', debug_flag=False):
		if chrom==self.current_chrom:
			return self.transmit_t

		self.current_chrom = chrom
		self.hapmap_fpath = self.hapmap_fpath_fmt.format(chrom=chrom)
		self.vcf_fpath = self.vcf_fpath_fmt.format(chrom=chrom)
		#print self.hapmap_fpath, self.vcf_fpath

		hm_df, hm_id_df= self.get_hets()
		sel_df,sel_pos_df = self.vcf_phased()
		#print hm_df
		#print sel_df
		subset_id = hm_df.index.intersection(sel_df.index)
		if debug_flag:
			print subset_id, hm_df.shape	
			print "SidPhased: {0}/{1}".format(len(subset_id), len(sel_df.index))

		hm_fdf=hm_df.ix[subset_id,:]
		sel_fdf=sel_df.ix[subset_id,:]
		#print hm_fdf
		#print sel_fdf
		t = sel_fdf[['A','B']]
		t['A_pat'] = hm_fdf['NA12891_A']==t['A']
		t['A_mat'] = hm_fdf['NA12892_A']==t['A']
		t['B_pat'] = hm_fdf['NA12891_A']==t['B']
		t['B_mat'] = hm_fdf['NA12892_A']==t['B']
	
		# remove SNPs where the alleles do not match
		# the parents 
		t = t.ix[na.logical_and((t['A_pat'].astype(int)+t['A_mat'].astype(int))==1,
														(t['B_pat'].astype(int)+t['B_mat'].astype(int))==1),:]
		self.transmit_t = t 
		#print t 
		self.transmit_s = t[['A_pat','A_mat','B_pat','B_mat']].astype(int).sum()
		# asserts one parent didn't donate both homologs
		a = na.argmax(self.transmit_s[['A_pat','A_mat']])
		b = na.argmax(self.transmit_s[['B_pat','B_mat']])
		assert a!=b
		
		return t
	def _map_parent_to_homolog(self, parent):
		# assumes transmit is stored
		return 'p'+na.argmax(self.transmit_s.loc[['A_%s'%parent, 'B_%s'%parent]])[0]

	def assign_dels_to_homolog(self, df, parent):		
		g = df.groupby('chr')
		df['class'] = '.'
		for chrom, g_df in g:
			if not str(chrom).startswith('chr'):
				chrom = "chr{0}".format(chrom)
			try: 
				self.snp_transmit(chrom)
				c = self._map_parent_to_homolog(parent)
				df.ix[g_df.index,'class']=c
			except IOError:
				pass
		return df

def cnv_hets(hapmap_fpath, c_id, pat_id, mat_id):
	df = pd.read_csv(hapmap_fpath, sep='\t', index_col=0)
	#v = df[[c_id,pat_id,mat_id]]
	v = df[[c_id,pat_id,mat_id]].dropna(axis=0).astype(int)
	return v, df[['chr','start','end']].ix[v.index,:]

def cnv_transmit(any_dels=False, inc_size=50, quake_flag=False):
	c_id = 'NA12878'
	pat_id = 'NA12891'
	mat_id = 'NA12892'

	with open('/media/T02/data/hapmap/fan_whole_haplotype_ids.txt', 'rb') as f:
		quake_hm_ids = [cnp.strip() for cnp in f]
	hm_df, hm_pos_df = cnv_hets('/media/T02/data/hapmap/hm3_cnv_submission.txt', c_id, pat_id, mat_id) 
	#print hm_df, hm_pos_df
	valid= pd.concat((hm_df, hm_pos_df), axis=1) 

	g_c = hm_df.groupby([c_id,])
	g = hm_df.groupby([c_id, pat_id, mat_id])
	#print "HetCases", g_c.size()
	#print "Cases", g.size()

	hom = g_c.get_group(0)
	hom_set = set([(0,0),(0,1),(1,0),(1,1)])
	#hom_set = set([(0,0),(0,1),(1,0)])

	hom_idx = hom.apply(lambda x: tuple(x[[pat_id, mat_id]]) in hom_set, axis=1)
	#print hom_idx.sum(), hom_idx.size
	hom_v = hm_pos_df.ix[hom_idx.index,:]
	hom_v['class'] = 'homA'
	#print hom_v	
	#print valid.ix[hom_v.index,:]

	inc = g.get_group((2,2,2))
	
	inc_iidx = na.arange(inc.index.size)
	na.random.shuffle(inc_iidx)
	inc_iidx = na.sort(inc_iidx[:inc_size])
	inc_idx = inc.ix[inc.index[inc_iidx],:]
	inc_v = hm_pos_df.ix[inc_idx.index,:]
	inc_v['class'] = 'inc'

	het = g_c.get_group(1)

	if any_dels:
		final = pd.concat([hm_pos_df.ix[hom_v.index,:], hm_pos_df.ix[het.index,:]]).sort(['chr','start','end'])
		final['chr'] = final['chr'].map(lambda x:"chr{0}".format(x))
		final.to_csv(sys.stdout, header=False, index=False,sep='\t')

		return

	snp = SNPTransmit('NA12891','NA12892')
	pat_set = set([(0,1),(0,2),(1,2)])
	pat_idx = het.apply(lambda x: tuple(x[[pat_id, mat_id]]) in pat_set, axis=1)
	pat_v = hm_pos_df.ix[pat_idx.index[pat_idx],:]
	#print pat_idx
	#print pat_v
	pat_v = snp.assign_dels_to_homolog(pat_v, 'pat')

	mat_set = set([(1,0),(2,0),(2,1)])
	mat_idx = het.apply(lambda x: tuple(x[[pat_id, mat_id]]) in mat_set, axis=1)
	mat_v = hm_pos_df.ix[mat_idx.index[mat_idx],:]
	#print mat_idx
	mat_v = snp.assign_dels_to_homolog(mat_v, 'mat')
	valid['pat'] = pat_idx
	valid['mat'] = mat_idx
	valid['hom'] = hom_idx
	valid['inc'] = inc_idx.apply(lambda x:(x==2).all(), axis=1)
	#print na.logical_and(pat_idx==mat_idx, pat_idx)
	#print "Intersection", mat_v.index.intersection(pat_v.index) 
	#print valid.dropna(how='all', subset=['pat','mat','hom','inc'])
	#final = pd.concat([hom_v,]).sort(['chr', 'start', 'end'])
	final = pd.concat([inc_v,hom_v,pat_v,mat_v]).sort(['chr', 'start', 'end'])
	final['chr'] = final['chr'].map(lambda x:"chr{0}".format(x))
	if quake_flag:
		final = final.ix[pd.Index(quake_hm_ids).intersection(final.index),:]
	final.to_csv(sys.stdout, header=False, index=False, sep='\t')

def mills_transmit(somatic_flag=False):
	mills_mid_inherit = pd.read_csv('/home/anand/Projects/assembly/data/mills_2010-08-table4.CEU1kb.mid_to_NA12878inherit.map', sep='\t', header=None, names=['MERGED_ID','INHERIT']).set_index('MERGED_ID')
	mills_mid = pd.read_csv('/home/anand/Projects/assembly/data/trio.2010_06.deletions.sites.1kb.bed', sep='\t', header=None, names=['chr','start','end','MERGED_ID']).set_index('MERGED_ID')
	#print mills_mid_inherit.ix[:10,:]
	#print mills_mid.ix[:10,:]
	hom_idx = mills_mid_inherit.index[mills_mid_inherit['INHERIT']=='homA']
	hom_v = mills_mid.ix[hom_idx,:]
	hom_v['class'] = 'homA'

	snp = SNPTransmit('NA12891','NA12892')
	pat_idx = mills_mid_inherit.index[mills_mid_inherit['INHERIT']=='pat']
	pat_v = mills_mid.ix[pat_idx,:]
	pat_v = snp.assign_dels_to_homolog(pat_v, 'pat')

	mat_idx = mills_mid_inherit.index[mills_mid_inherit['INHERIT']=='mat']
	mat_v = mills_mid.ix[mat_idx,:]
	mat_v = snp.assign_dels_to_homolog(mat_v, 'mat')

	final = pd.concat([hom_v,pat_v,mat_v]).sort(['chr', 'start', 'end'])
	final = final.ix[final['class']!='.',:]

	chrom_name_map = dict(zip(map("chr{0:d}".format, range(1,23)), range(1,23)))
	chrom_name_map['chrX']=23 
	chrom_name_map['chrY']=24 

	final['new_chr'] = final['chr'].map(chrom_name_map.get)
	final = final.sort(['new_chr','start','end']).drop(['new_chr',],axis=1)

	if somatic_flag:
		final = mills_mid.ix[mills_mid.index - final.index,:]		

	final.to_csv(sys.stdout, header=False, index=False, sep='\t')

def snp_transmit():
	w = SNPTransmit('NA12891','NA12892')
	t= w.snp_transmit('chr20', debug_flag=True) 
	s = t[['A_pat','A_mat','B_pat','B_mat']].astype(int).sum()
	
	print s
	parents = ['pat','mat']
	a = parents[na.argmax(s[['A_pat','A_mat']])]
	b = parents[na.argmax(s[['B_pat','B_mat']])]
	assert a!=b
	
if __name__ == '__main__':


	#snp_transmit()
	cnv_transmit(any_dels=True, quake_flag=False)
	#mills_transmit(somatic_flag=True)
