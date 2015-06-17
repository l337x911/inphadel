from svphase.utils.config import PRECISION
from itertools import izip, count
import pandas as pd
import numpy as np

class FeatureSubset(object):
	def __init__(self):
		pass
	def _nonzero(self, df):
		return ((df > PRECISION)).apply(any, axis=0)
	def _nonzero_per_data_subset(self, d):
		cols = range(0, len(d.loci))
		wgs = self._nonzero(d._rm_break_info(d.wgs).loc[:, cols].astype(float))
		wgs_allele = self._nonzero(d._rm_break_info(d.wgs_allele).loc[:, cols].astype(float))
		hic_allele = self._nonzero(d._rm_break_info(d.hic_allele).loc[:, cols].astype(float))
		hic_allele_ncut = self._nonzero(d._rm_break_info(d.hic_allele_ncut).loc[:, cols].fillna(0).astype(float))
		return wgs, wgs_allele, hic_allele,hic_allele_ncut
	def data_nonzero(self, d):
		wgs, wgs_allele, hic_allele, hic_allele_ncut = self._nonzero_per_data_subset(d)
		return np.logical_or(np.logical_or(hic_allele, hic_allele_ncut), wgs_allele)
		
	def data_subset(self, d):
		""" get combination of data subets 
		Args:
			d - data object """
		return [d.wgs, d.wgs_allele, d.hic_allele, d.hic_allele_ncut]
	
class HicOnlySubset(FeatureSubset):
	def __init__(self):
		super(FeatureSubset)
	def data_nonzero(self, d):
		wgs, wgs_allele, hic_allele, hic_allele_ncut = self._nonzero_per_data_subset(d)
		return np.logical_or(hic_allele, hic_allele_ncut)
	def data_subset(self, d):
		return [d.wgs, d.hic_allele, d.hic_allele_ncut]

class WgsOnlySubset(FeatureSubset):
	def __init__(self):
		super(FeatureSubset)
	def data_nonzero(self, d):
		wgs, wgs_allele, hic_allele, hic_allele_ncut = self._nonzero_per_data_subset(d)
		return wgs_allele.values
	def data_subset(self, d):
		return [d.wgs, d.wgs_allele]

class AlleleLessSubset(FeatureSubset):
	def __init__(self):
		super(FeatureSubset)
	def data_nonzero(self, d):
		wgs, wgs_allele, hic_allele, hic_allele_ncut = self._nonzero_per_data_subset(d)
		return wgs.values
	def data_subset(self, d):
		return [d.wgs,]

class Features(object):
	""" Contains features from filled data features """
	def __init__(self, filled_data_objs, feature_subset=None):
		"""
		
		Args:
			filled_data_objs - Data Objects with filled/adapted members
		
		"""
		self._data = filled_data_objs
		self._nonzero_loci = None
		self._features = None # Samples are on axis 0, and Features are on axis 1
		if feature_subset is None:
			self._feature_subset = FeatureSubset()
			
	def _prep_data(self, d):
		d._add_data_type()
		d._rm_wrong_orients()
		d._rm_dist()
		d._rm_hic_arounds()
		d._rm_hic_ncut_cov()
	def get_nonzero(self):
		if not self._nonzero_loci is None:
			return self._nonzero_loci
		
		self._nonzero_loci = np.concatenate(map(self._feature_subset.data_nonzero, self._data))
		return self._nonzero_loci
	
	def get_features(self):
		if not self._features is None:
			return self._features
		
		feats = []
		
		for d in self._data:
			self._prep_data(d)
			c = self._feature_subset.data_subset(d)
			f = pd.concat(map(d._rm_break_info,c), axis=0, ignore_index=True)
			cols_set = set(list(f.columns))
			h = [i for i in d.row_index_header if i in cols_set]
		
			f = f.set_index(h)
			d.prep_once_flag = True
			feats.append(f)
					
		for i,j in izip(feats, feats[1:]):
			assert (i.index == j.index).all()
		features = pd.concat(feats, axis=1, ignore_index=True)
		features = features.fillna(0)
		features = features.astype(float)
		 
		self._features = features.T
		return self._features
	def get_nonzero_features(self):
		return self.get_features().loc[self.get_nonzero(),:]

class SimpleSumFeatures(Features):
	def __init__(self, filled_data_objs, feature_subset=None):
		Features.__init__(self, filled_data_objs, feature_subset)
	def get_features(self):
		if not self._features is None:
			return self._features

		feats = []
		
		for d in self._data:
			self._prep_data()
			c = self._feature_subset.data_subset(d)
			f = pd.concat(map(d._rm_break_info,c), axis=0, ignore_index=True)
			cols_set = set(list(f.columns))
			h = [i for i in d.row_index_header if i in cols_set]

			f = f.set_index(h)
			d.prep_once_flag = True
			feats.append(f)
		
			
		raise NotImplementedError()

class SqDiffAlleleFeatures(Features):
	def __init__(self, filled_data_objs, feature_subset=None):
		Features.__init__(self, filled_data_objs, feature_subset)
	def get_features(self):
		raise NotImplementedError()


