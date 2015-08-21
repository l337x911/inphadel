import sys
import pandas as pd
from itertools import chain
import numpy as np
import os
from svphase.utils.common import logger
from svphase.learn.features import Features, PreloadFeatures, SimpleSumFeatures

class ClassLabel(object):
	def __init__(self):
		self.classes = ['pA','pB','hom', 'inc']
		self.int_to_str_dict = dict(zip(range(len(self.classes)), self.classes))
		self.str_to_int_dict = dict(zip(self.classes, range(len(self.classes))))
	def is_deletion_on_a(self, c):
		if c=='pA' or c=='hom':
			return True
		else: 
			return False
	def is_deletion_on_b(self, c):
		if c=='pB' or c=='hom':
			return True
		else: 
			return False
	def is_class(self, c):
		return c in self.classes

class Evaluation(object):
	def __init__(self, feature_subset=None):
		self.data = []
		self.feats = None
		self.labels = None
		self.feature_subset = feature_subset
		self.label_obj = ClassLabel()
	def add_data(self, data):
		self.data.append(data)
	def set_features(self, save_prefix=None, preload_prefix=None, simple_sum_flag=False):
		if not self.feats is None:
			logger.warning('Reloading features')

		if (not preload_prefix is None) and os.path.isfile(preload_prefix + ".pkl"):
			logger.info("Using preloaded features.")
			self.feats = PreloadFeatures(preload_prefix) 
		elif simple_sum_flag:
			self.feats = SimpleSumFeatures(self.data, self.feature_subset)	
		else:
			self.feats = Features(self.data, self.feature_subset)

		feats = self.feats.get_features()
		logger.info("Loaded samples: %d, features: %d", feats.shape[0], feats.shape[1])
		feats = self.feats.get_nonzero_features()
		logger.info("Non-zero samples: %d, features: %d", feats.shape[0], feats.shape[1])
		if not save_prefix is None:
			self.feats.save(save_prefix)

	def set_labels(self, save_prefix=None, preload_prefix=None):
		if not self.labels is None:
			logger.warning('Reloading labels')
		if self.feats is None:
			logger.error('Truth labels need features to be loaded')
			sys.exit(1)

		if (not preload_prefix is None) and  os.path.isfile(preload_prefix + '.txt'):
			with open(preload_prefix + '.txt', 'rb') as f:
				self.labels = [line.strip() for line in f]
		else:
			self.labels = list(chain(*(d.truth for d in self.data)))
			if not save_prefix is None:
				with open(save_prefix+'.txt', 'wb') as f:
					f.write('\n'.join(self.labels) + '\n')

		logger.info("Loaded %d truth labels", len(self.labels))

	def manual_check(self):
		pd.set_option('display.max_rows', 999)
		pd.set_option('display.max_columns', 999)
		pd.set_option('display.height', 999)
		pd.set_option('display.width', 999)
		data = self.feats.get_nonzero_features()
		data = data.fillna(0).astype(np.float64)
		labels_str = np.array(self.labels)[self.feats.get_nonzero()]
		labels_int = map(self.label_obj.str_to_int_dict.get, labels_str)
		label_str_idx = pd.MultiIndex.from_tuples([tuple(['truth_str',] + [np.nan,]*(len(data.columns.names)-1)),], names=data.columns.names)
		label_int_idx = pd.MultiIndex.from_tuples([tuple(['truth_int',] + [np.nan,]*(len(data.columns.names)-1)),], names=data.columns.names)

		inspect = pd.concat([pd.DataFrame(labels_str, index=data.index, columns=label_str_idx), pd.DataFrame(labels_int, index=data.index, columns=label_int_idx), data], names=data.columns.names, axis=1)
		print inspect
		#logger.info('\n%s', inspect)

if __name__=='__main__':
	import argparse 
	from pkg_resourcs import Requirement, resource_filename

	parser = argparse.ArgumentParser(description='Compare predictions with a truth set.')
