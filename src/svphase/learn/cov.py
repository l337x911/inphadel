""" Counts reads from a SAM or simplified binary files.
"""
import pysam
from svphase.utils.config import REFERENCE_FASTA, HIND3_STR
from svphase.parser import ReadsParserSAM, ReadsParserDat
from collections import defaultdict

from itertools import izip, chain, repeat, count
import bisect
import re
import os
import pandas as pd
import numpy as na
# from matplotlib import pyplot as plt

'''
 Module computes feats for a chunk of data.
 A "data" object calculates the read counts from various separated
 mapped read files.
'''
def line_format(vlist):
	return '\t'.join(["% 10d" % i for i in vlist])

def float_format(vlist):
	return ' '.join(["% 11d" % vlist[0], ] + ["% 1.5f" % i for i in vlist[1:]])

orient_order = (('+', '-'), ('-', '+'), ('-', '-'), ('+', '+'),)
parser_bam = ReadsParserSAM()
parser_dat = ReadsParserDat()

class NearestCutSite(object):
	def __init__(self, ref_obj, contig):
		seq = ref_obj.fetch(reference=contig).upper()

		hind3 = re.compile(HIND3_STR)
		self.pos = sorted([m.start() + 1 for m in re.finditer(hind3, seq)])
	def _find_below(self, a, lo=0):
		idx = bisect.bisect_right(self.pos, a, lo=lo)
		if idx == 0:
			return None
		return idx, self.pos[idx - 1]
	def _find_above(self, a, hi=None):
		if hi is None: hi = len(self.pos)
		idx = bisect.bisect_left(self.pos, a, hi=hi)
		if idx == len(self.pos) - 1:
			return None
		return idx, self.pos[idx]
	def find_below(self, a):
		val = self._find_below(a)
		if val is None:
			return na.NaN
		return val[1]
	def find_above(self, a):
		val = self._find_above(a)
		if val is None:
			return na.NaN
		return val[1]
	def find_closest_to_b_in(self, a, b):
		val1 = self._find_above(a)
		if val1 is None:
			return na.NaN
		val2 = self._find_below(b, lo=val1[0])
		if val2 is None or val1[0] == val2[0]:
			return na.NaN
		return val2[1]
	def find_closest_to_a_in(self, a, b):
		val1 = self._find_below(b)
		if val1 is None:
			return na.NaN
		val2 = self._find_above(a, hi=val1[0])
		if val2 is None or val1[0] == val2[0]:
			return na.NaN
		return val2[1]			

def test_ncutsites():
	pos = [100000, 105000, 200000, 500000]
	
	finder = NearestCutSite(pysam.FastaFile(REFERENCE_FASTA), 'chr20')

	print "Test Above and Below"
	for p in pos:
		print finder.find_above(p), p, finder.find_below(p)
	
	print "Test Ranges"
	for p1, p2 in zip(pos[::2], pos[1::2]):
		print p1, finder.find_closest_to_a_in(p1, p2), finder.find_closest_to_b_in(p1, p2), p2

def get_parser(reads_fpath):
	if reads_fpath.endswith('.bam'):
		return parser_bam
	elif reads_fpath.endswith('.dat'):
		return parser_dat
	else:
		raise

class PairedReads(object):
	def __init__(self):
		self.counts = pd.DataFrame(na.zeros((len(orient_order), 2), dtype=na.uint16), index=orient_order, columns=(1, 2))
	def load(self, bam_fpath, c, a, b):
		links = defaultdict(list)
		parser = get_parser(bam_fpath)
		for row in parser.get_reads(bam_fpath, c, a, b):
			links[parser.get_id(row)].append(row)

		for pe_read in links.itervalues():
			strd = parser.get_strandedness(pe_read[0])
			if len(pe_read) > 1:
				self.counts[2][strd] += 1
			self.counts[1][strd] += 1

class PairedReadsAcross(object):
	def __init__(self):
		self.counts = pd.Series(na.zeros(len(orient_order), dtype=na.uint16), index=orient_order)
	def load(self, bam_fpath, c1, a1, b1, c2, a2, b2):
		parser = get_parser(bam_fpath)

		leftlink = dict([(parser.get_id(row), row) for row in parser.get_reads(bam_fpath, c1, a1, b1)])
		rightlink = dict([(parser.get_id(row), row) for row in parser.get_reads(bam_fpath, c2, a2, b2)])
		for row in map(leftlink.__getitem__, set(leftlink.viewkeys()).intersection(set(rightlink.viewkeys()))):
			self.counts[parser.get_strandedness(row)] += 1

class ReadCollectionFactory(object):
	def __init__(self):
		pass
	def manufacture(self, obj_class, args):
		obj = obj_class()
		obj.load(*args)
		return obj

factory = ReadCollectionFactory()
cutsiteFinder = {}
def get_cutsite_finder(ref_fpath, contig):
	if not cutsiteFinder.has_key((ref_fpath, contig)):
		cutsiteFinder[(ref_fpath, contig)] = NearestCutSite(pysam.FastaFile(ref_fpath), contig)
	return cutsiteFinder[(ref_fpath, contig)]

class OneBreak(object):
	def __init__(self, window):
		self.w = window
		self.break_info = ["c1", "a1", "b1"]
		self.signal_label = ["cov", "w.around_w.around", "w.around_w.inside", "w.inside_w.around"]
		self.signal = None
		orient1, orient2 = zip(*orient_order)
		self.signal_header = list(izip(self.break_info, repeat(''), repeat(''))) + list(chain(*[izip(repeat(i), orient1, orient2) for i in self.signal_label]))
		self.signal_header = pd.DataFrame(self.signal_header, columns=('wtype', 'orient1', 'orient2',))
		self.header_depth = 3
	def _norm(self, a, b):
		return 1000. / (b - a)
	def load(self, bam_fpath, c, a, b):
		cov = factory.manufacture(PairedReads, (bam_fpath, c, a, b))
		norm = self._norm(a, b)
	
		a_a = factory.manufacture(PairedReadsAcross, (bam_fpath, c, a - self.w, a, c, b, b + self.w))
		a_i = factory.manufacture(PairedReadsAcross, (bam_fpath, c, a - self.w, a, c, a, a + self.w))
		i_a = factory.manufacture(PairedReadsAcross, (bam_fpath, c, b - self.w, b, c, b, b + self.w))

		self.signal = self.signal_header.copy()
		self.signal[0] = [c, a, b, ] + list(cov.counts[1] * norm) + list(a_a.counts * norm) + list(a_i.counts * norm) + list(i_a.counts * norm)

class OneBreakCutSite(OneBreak):
	def __init__(self, ref_dir, contig, window):
		OneBreak.__init__(self, window)
		self.ref_fpath = ref_dir
		self.contig = contig
		self.w /= 2
		self.dist_label = ["dist_cut1", "dist_cut2"]
		
		self.signal_header['dist'] = ''
		for wtype in self.signal_label[1:]:
			td1 = pd.DataFrame([[wtype, '', '', self.dist_label[0], ]], columns=['wtype', 'orient1', 'orient2', 'dist'])
			td2 = pd.DataFrame([[wtype, '', '', self.dist_label[1], ]], columns=['wtype', 'orient1', 'orient2', 'dist'])

			self.signal_header = pd.concat([self.signal_header, td1, td2], axis=0, ignore_index=True)
		
		self.header_depth += 1
	def load(self, bam_fpath, c, a, b):
		cov = factory.manufacture(PairedReads, (bam_fpath, c, a, b))
		norm = self._norm(a, b)
		
		nan_series = pd.Series([na.nan] * 4, index=orient_order)
		cf = get_cutsite_finder(self.ref_fpath, self.contig)
		cut_a1 = cf.find_below(a)
		cut_a2 = cf.find_above(b)
		cut_i1 = cf.find_closest_to_a_in(a, b)
		cut_i2 = cf.find_closest_to_b_in(a, b)
	
		if cut_a1 is na.NaN or cut_a2 is na.NaN:
			a_a = nan_series
		else:
			a_a = factory.manufacture(PairedReadsAcross, (bam_fpath, c, cut_a1 - self.w, cut_a1 + self.w, c, cut_a2 - self.w, cut_a2 + self.w)).counts * norm

		if cut_a1 is na.NaN or cut_i1 is na.NaN:
			a_i = nan_series
		else:
			a_i = factory.manufacture(PairedReadsAcross, (bam_fpath, c, cut_a1 - self.w, cut_a1 + self.w, c, cut_i1 - self.w, cut_i1 + self.w)).counts * norm
	
		if cut_i2 is na.NaN or cut_a2 is na.NaN:
			i_a = nan_series
		else:
			i_a = factory.manufacture(PairedReadsAcross, (bam_fpath, c, cut_i2 - self.w, cut_i2 + self.w, c, cut_a2 - self.w, cut_a2 + self.w)).counts * norm

		self.signal = self.signal_header.copy()

		self.signal[0] = [c, a, b, ] + list(cov.counts[1] * norm) + list(a_a * norm) + list(a_i * norm) + list(i_a * norm) + [a - cut_a1, cut_a2 - b, a - cut_a1, cut_i1 - a, b - cut_i2, cut_a2 - b] 

class TwoBreak(object):
	def __init__(self):
		pass

def perform_one_break_allele_cut(ref_dir, loci, bama_fpath, bamb_fpath, out_flag=True):
	signal_df = None
	for idx, (c, a, b) in enumerate(loci):
		objA = OneBreakCutSite(ref_dir, c, 1000)
		objA.load(bama_fpath, c, a, b)
		objB = OneBreakCutSite(ref_dir, c, 1000)
		objB.load(bamb_fpath, c, a, b)
		objA.signal['allele'] = 'A'
		objB.signal['allele'] = 'B'
		if signal_df is None:
			signal_df = pd.concat([objA.signal, objB.signal[3:]], axis=0, ignore_index=True)
			continue
		signal_df[idx] = pd.concat([objA.signal, objB.signal[3:]], axis=0, ignore_index=True)[0]

	cols = list(signal_df.columns)
	cols[1:6] = cols[0:5]
	cols[0] = 'allele'
	signal_df = signal_df[cols]
	if out_flag:
		print signal_df
	return signal_df

def perform_one_break_allele(loci, bama_fpath, bamb_fpath, out_flag=True):
	signal_df = None
	for idx, (c, a, b) in enumerate(loci):
		objA = OneBreak(1000)
		objA.load(bama_fpath, c, a, b)
		objB = OneBreak(1000)
		objB.load(bamb_fpath, c, a, b)
		objA.signal['allele'] = 'A'
		objB.signal['allele'] = 'B'

		if signal_df is None:
			signal_df = pd.concat([objA.signal, objB.signal[3:]], axis=0, ignore_index=True)
			continue

		signal_df[idx] = pd.concat([objA.signal, objB.signal[3:]], axis=0, ignore_index=True)[0]
	signal_df = signal_df.rename(index=dict(izip(signal_df.index, na.arange(signal_df.index.size))))	
	cols = list(signal_df.columns)
	cols[1:5] = cols[0:4]
	cols[0] = 'allele'
	signal_df = signal_df[cols]
	if out_flag:
		# print signal_df
		orients = signal_df.groupby(['orient1', 'orient2'])
		print signal_df.ix[sorted(orients.groups[('', '')] + orients.groups[('+', '-')])]
	return signal_df

def perform_one_break(loci, bam_fpath, out_flag=True):
	signal_df = None
	for idx, (c, a, b) in enumerate(loci):
		objA = OneBreak(1000)
		objA.load(bam_fpath, c, a, b)
		if signal_df is None:
			signal_df = objA.signal.copy()
			continue
		signal_df[idx] = objA.signal[0]

	if out_flag:
		print signal_df
	return signal_df

def orient_wrap_hic(df):
	col_start = int(na.where(df.columns == 0)[0])
	new_df = df.set_index(df.columns.tolist()[:col_start])
	
	info = new_df.xs(('', ''), level=('orient1', 'orient2')) 
	info['orient'] = ''
	oppo = new_df.xs(('+', '-'), level=('orient1', 'orient2')) + new_df.xs(('-', '+'), level=('orient1', 'orient2'))
	oppo['orient'] = 'opposite'
	same = new_df.xs(('-', '-'), level=('orient1', 'orient2')) + new_df.xs(('+', '+'), level=('orient1', 'orient2'))
	same['orient'] = 'same'

	new_df = pd.concat([info, oppo, same, ]).sortlevel(0)
	new_df = new_df.reset_index()
	cols = set(df.columns.tolist()[:col_start])

	cols.discard('orient2')
	cols.discard('orient1')
	new_cols = [i for i in df.columns[:col_start] if i in cols] + ['orient', ] 

	new_df.set_index(new_cols, inplace=True)
	# print df
	return new_df.reset_index()

class DataPreloadedException(Exception):
	pass

class Data(object):
	def __init__(self, ref_fpath=None, contig=None):
		self.loci = []
		self.wgs = None
		self.wgs_allele = None
		self.hic_allele = None
		self.hic_allele_ncut = None
		self.break_info = ['c1', 'a1', 'b1']
		self.row_index_header = ['data_type', 'allele', 'wtype', 'orient', 'orient1', 'orient2', 'dist']
		self.prep_once_flag = False
		if ref_fpath is None:
			self.ref_fpath = REFERENCE_FASTA
		else:
			self.ref_fpath = ref_fpath
		self.contig = contig

	def print_dataframes(self, show=7):
		print "# WGS"
		print self.wgs.ix[:, :show]
		print "# WGS Allele"
		print self.wgs_allele.ix[:, :show]
		print "# HiC Allele"
		print self.hic_allele.ix[:, :show]
		print "# HiC Allele Neighbor Cuts"
		print self.hic_allele_ncut.ix[:, :show]

	def fill(self, adaptor, save_flag=False):
		if adaptor is None:
			adaptor = DataAdaptor()
		try:
			self.wgs = perform_one_break(self.loci, self.wgs_fpath, out_flag=False)
			self.wgs_allele = perform_one_break_allele(self.loci, self.wgsa_fpath, self.wgsb_fpath, out_flag=False)
			self.hic_allele = orient_wrap_hic(perform_one_break_allele(self.loci, self.hica_fpath, self.hicb_fpath, out_flag=False))
			self.hic_allele_ncut = orient_wrap_hic(perform_one_break_allele_cut(self.ref_fpath, self.loci, self.hica_fpath, self.hicb_fpath, out_flag=False))
			adaptor.adapt(self)
		except:
			print "Errored data fill:", self.loci
			raise		

		if save_flag and not self.check_on_disk():
			self.save()

	def _set_loci(self, sv_fpath):
		with open(sv_fpath, 'rb') as f:
			for line in f:
				if line.startswith('track'): continue
				tokens = line.strip().split('\t')
				self.loci.append((tokens[0], int(tokens[1]), int(tokens[2])))
		self.loci = filter(lambda x:x[0] == self.contig, self.loci)

	def _rm_break_info(self, df):
		#if self.prep_once_flag: raise DataPreloadedException()
		break_info_idx = [i for i, v in enumerate(df['wtype']) if v in self.break_info]
		return df.drop(df.index[break_info_idx], axis=0)

	def _rm_wrong_orients(self, df=None):
		if self.prep_once_flag: raise DataPreloadedException()
		if not df is None:
			disc = set([('-', '-'), ('+', '+'), ('-', '+')])
			idx = [i for i, v in df[['orient1', 'orient2']].iterrows() if tuple(v) in disc]
			return df.drop(idx, axis=0)
		else:
			self.wgs = self._rm_wrong_orients(self.wgs)
			self.wgs_allele = self._rm_wrong_orients(self.wgs_allele)

	def _rm_hic_arounds(self, df=None):
		if self.prep_once_flag: raise DataPreloadedException()
		if not df is None:
			# print "Rm arounds"
			# print df
			idx = [i for i, v in enumerate(df['wtype'] == 'w.around_w.around') if v]
			return df.drop(df.index[idx], axis=0)
		else:
			self.hic_allele = self._rm_hic_arounds(self.hic_allele)
			self.hic_allele_ncut = self._rm_hic_arounds(self.hic_allele_ncut)

	def _rm_hic_ncut_cov(self, df=None):
		if self.prep_once_flag: raise DataPreloadedException()
		if not df is None:
			# print "Rm cov"
			# print df
			idx = [i for i, v in enumerate(df['wtype'] == 'cov') if v]
			return df.drop(df.index[idx], axis=0)
		else:
			self.hic_allele_ncut = self._rm_hic_ncut_cov(self.hic_allele_ncut)

	def _rm_dist(self, df=None):
		if self.prep_once_flag: raise DataPreloadedException()
		if not df is None:
			s = set(['dist_cut1', 'dist_cut2'])
			idx = [i for i, v in enumerate(df['dist']) if v in s]
			df = df.drop(idx, axis=0)
			return df.drop(['dist', ], 1)
		else:
			self.hic_allele_ncut = self._rm_dist(self.hic_allele_ncut)

	def _add_data_type(self):
		if self.prep_once_flag: raise DataPreloadedException()
		self.wgs['data_type'] = 'wgs'
		self.wgs_allele['data_type'] = 'wgs_a'
		self.hic_allele['data_type'] = 'hic_a'
		self.hic_allele_ncut['data_type'] = 'hic_ac'

	def _apply_func_for_allele_pair(self, df, func):
		df = self._rm_break_info(df)
		
		a = (df.ix[df['allele'] == 'A', :]).drop(['allele', ], 1)
		b = (df.ix[df['allele'] == 'B', :]).drop(['allele', ], 1)
		b = b.rename(index=dict(zip(b.index, a.index)))		

		assert (a.columns == b.columns).all()
		assert (a.index == b.index).all()
		col_start = int(na.where(a.columns == 0)[0])
		new_idx_cols = list(a.columns[:col_start])
		a = a.set_index(new_idx_cols).astype(float)
		b = b.set_index(new_idx_cols).astype(float)
		return func(a, b)

	def _push_last_col_to_first(self, df):
		cols = list(df.columns)
		cols = cols[-1:] + cols[:-1]
		return df[cols]


	def get_sq_diff_allele_features(self, combines=None):
		self._add_data_type()
		self._rm_wrong_orients()
		self._rm_dist()
		self._rm_hic_arounds()
		self._rm_hic_ncut_cov()
	 
		nonallele_df = self._push_last_col_to_first(self._rm_break_info(self.wgs))
		allele_norm_dfs = []
		sq_diff_func = lambda x, y:na.power((x - y), 2)

		if combines is None:
			c = [self.wgs_allele, self.hic_allele, self.hic_allele_ncut]
		elif combines == 'hic_only':
			c = [self.hic_allele, self.hic_allele_ncut]
		elif combines == 'wgs_only':
			c = [self.wgs_allele]
		elif combines == 'only':
			c = []
		
		for df in map(self._rm_break_info, c):
			df = self._push_last_col_to_first(df)
			allele_norm_dfs.append(self._apply_func_for_allele_pair(df, sq_diff_func).reset_index())
		
		features = pd.concat([self._rm_break_info(self.wgs), ] + allele_norm_dfs, axis=0, ignore_index=True)

		h = [i for i in self.row_index_header if i in features.columns]
		cols = h + range(0, len(self.loci))
		features = features[cols]
		return features.rename(index=dict(zip(features.index, count(0))))
		# return feats 

	def index_features(self, combines=None, feature_typeset=None):
		if feature_typeset == 'simple_sum':
			features = self.get_simple_sum_features(combines)
		elif feature_typeset == 'binary_allele':
			features = self.get_sq_diff_allele_features(combines)
		else:
			features = self.get_features(combines)

		cols_set = set(list(features.columns))
		h = [i for i in self.row_index_header if i in cols_set]
		return features.set_index(h)

class FileStructureData(Data):
	def __init__(self, wgs_prefix_path, hic_prefix_path, sv_fpath, contig, file_fmt='bam', ref_fpath=None):

		fmt_str = os.path.join("{prefix}","{contig}.{type}.{ext}")

		Data.__init__(self, ref_fpath=ref_fpath, contig=contig)
		self.wgs_fpath = fmt_str.format(prefix=wgs_prefix_path, contig=contig, ext=file_fmt, type='all')
		self.wgsa_fpath = fmt_str.format(prefix=wgs_prefix_path, contig=contig, ext=file_fmt, type='A')
		self.wgsb_fpath = fmt_str.format(prefix=wgs_prefix_path, contig=contig, ext=file_fmt, type='B')

		self.hica_fpath = fmt_str.format(prefix=hic_prefix_path, contig=contig, ext=file_fmt, type='A')
		self.hicb_fpath = fmt_str.format(prefix=hic_prefix_path, contig=contig, ext=file_fmt, type='B')
		
		self._set_loci(sv_fpath)


class FileStructureDataWithTruth(Data):
	def __init__(self, wgs_prefix_path, hic_prefix_path, sv_fpath, contig, file_fmt='bam', ref_fpath=None):

		fmt_str = os.path.join("{prefix}","{contig}.{type}.{ext}")

		Data.__init__(self, ref_fpath=ref_fpath, contig=contig)
		self.wgs_fpath = fmt_str.format(prefix=wgs_prefix_path, contig=contig, ext=file_fmt, type='all')
		self.wgsa_fpath = fmt_str.format(prefix=wgs_prefix_path, contig=contig, ext=file_fmt, type='A')
		self.wgsb_fpath = fmt_str.format(prefix=wgs_prefix_path, contig=contig, ext=file_fmt, type='B')

		self.hica_fpath = fmt_str.format(prefix=hic_prefix_path, contig=contig, ext=file_fmt, type='A')
		self.hicb_fpath = fmt_str.format(prefix=hic_prefix_path, contig=contig, ext=file_fmt, type='B')
		
		self.truth = []
		self._set_loci_with_truth(sv_fpath)

	def _set_loci_with_truth(self, sv_fpath):
		with open(sv_fpath, 'rb') as f:
			for line in f:
				if line.startswith('track'): continue
				tokens = line.strip().split('\t')
				if tokens[0]!=self.contig: continue
				self.loci.append((tokens[0], int(tokens[1]), int(tokens[2])))
				self.truth.append(tokens[3])

class DataAdaptor(object):
	def __init__(self):
		pass
	def adapt(self, data_obj):
		pass

class RPKMAdaptor(DataAdaptor):
	def __init__(self, num_wgs_reads, num_hic_reads):
		DataAdaptor.__init__(self)
		self.num_wgs_reads = num_wgs_reads
		self.num_hic_reads = num_hic_reads

		self.norm_wgs = 1000000. / num_wgs_reads
		self.norm_hic = 1000000. / num_hic_reads	

	def _apply_to_df(self, data_obj, df, norm):	
		col_start = int(na.where(df.columns == 0)[0])
		# print df.index
		# print df[df.columns[col_start:]]
		break_info_idx = max([i for i, v in enumerate(df['wtype']) if v in data_obj.break_info]) + 1
		i, c = df.index[break_info_idx:], df.columns[col_start:]
		df.ix[i, c] = df.ix[i, c].astype(float) * norm

	def adapt(self, data_obj):
		self._apply_to_df(data_obj, data_obj.wgs, self.norm_wgs) 
		self._apply_to_df(data_obj, data_obj.wgs_allele, self.norm_wgs) 
		self._apply_to_df(data_obj, data_obj.hic_allele, self.norm_hic) 
		self._apply_to_df(data_obj, data_obj.hic_allele_ncut, self.norm_hic) 
		return data_obj


if __name__ == '__main__':
	pass
