"""InPhaDel: Genotypes and phase deletions on a single chromosome using a specific classification model

Trains models for phasing deletions using underlying WGS+HiC data 

"""
import sys
import os
import pickle
import pandas as pd
import numpy as np
import warnings
from itertools import chain

from sklearn import svm
from sklearn import ensemble
from sklearn.cross_validation import StratifiedKFold
from sklearn.grid_search import GridSearchCV
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import accuracy_score


from svphase.utils.config import RANDOM_STATE,DATA_PREFIX
from svphase.learn.cov import FileStructureDataWithTruth, RPKMAdaptor
from svphase.learn.features import Features, PreloadFeatures
from svphase.utils.common import logger
from svphase.inphadel import default_arguments, ClassLabel

class Model(object):
	def __init__(self, model_pkl, random_state):
		self.pkl = model_pkl
		self.random_state = random_state
		self.clf = None
		self.params = {}

class SVMModel(Model):
	def __init__(self, model_pkl, random_state):
		Model.__init__(self, model_pkl, random_state)
		self.clf = svm.SVC(kernel='linear', probability=True, random_state=self.random_state)
		self.params = {'C':[1,10,100]}

class Trainer(object):
	def __init__(self, k=5, feature_subset=None): 
		self.data = []
		self.feats = None
		self.labels = None
		self.k = k
		self.feature_subset = feature_subset		
		self.inner_models = []
		self.inner_accuracies = []
		self.label_obj = ClassLabel()		

	def add_data(self, data):
		self.data.append(data)
	def set_features(self, save_prefix=None, preload_prefix=None):
		if not self.feats is None:
			logger.warning('Reloading features')

		if (not preload_prefix is None) and os.path.isfile(preload_prefix + ".pkl"):
			logger.info("Using preloaded features.")
			self.feats = PreloadFeatures(preload_prefix) 
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
	
	def check_stable(self, scores, params_df):
		thresh = 0.10
		deviation = max(scores)-min(scores)

		if deviation>thresh:
			logger.warning('Model test accuracies deviate more than {thresh:0.1f}%, Deviation {dev:0.4f}'.format(thresh=thresh*100, dev=deviation))

		common_params = params_df.groupedby(list(params_df.columns))
		if len(common_params.groups)>1:
			logger.warning('Model had unstable parameters\n{len:%s}'.format(common_params.apply(len)))

	def train(self, model):
		data = self.feats.get_nonzero_features()
		labels_str = np.array(self.labels)[self.feats.get_nonzero()]
		labels_int = map(self.label_obj.str_to_int_dict.get, labels_str)
		logger.debug(labels_str)
		logger.debug(labels_int)
		labels = np.array(labels_str, dtype=int)

		# Applies Nested Cross Validation on data + labels
		outer_skf = StratifiedKFold(labels, self.k, shuffle=True, random_state=model.random_state)

		scores = []
		best_params = []

		for fold, (train_index, test_index) in enumerate(outer_skf):
			# Perform a new inner cross validation
			outer_skf = StratifiedKFold(labels[train_index], self.k, shuffle=True, random_state=model.random_state)
			grid_clf = GridSearchCV(estimator=model.clf, param_grid=model.params, scoring=accuracy_score, n_jobs=-1, cv=outer_skf, refit=True, verbose=True)
			grid_clf.fit(data[train_index], labels[train_index])
			
			test_accuracy = accuracy_score(labels[test_index], grid_clf.best_estimator_.predict(data[test_index]))
			scores.append(test_accuracy)
			best_params.append(pd.Series(grid_clf.best_params_))

			logger.debug('OUTER: %0.4f test accuracy with params %s on fold %d', test_accuracy, grid_clf.best_params_, fold)
			for pt in grid_clf.grid_scores_:
				logger.debug(' INNER: %.4f avg accuracy with params %s and scores %s', pt.mean_validation_score, pt.parameters, ','.join(map('{0:0.4f}'.format,pt.cv_validation_scores)))
			
		self.check_stable(scores, pd.concat(best_params, axis=1).T)	

class PrepTrainingFileStructure(object):
	def __init__(self, idir, ftype):
		assert ftype=='dat' or ftype=='bam'
		self.ftype = ftype
		if ftype == 'bam':
			self.idx_ftype = 'bai'
		elif ftype == 'dat':
			self.idx_ftype = 'npz'
		else:
			logger.error('No valid read filetype: %s', ftype)
			sys.exit(1)

		self.idir = idir
		self.sv_fpath = None
		self.wgs_read_count = None	
		self.hic_read_count = None	
		self.loci = None

	def load(self):
		""" Validate File Structure """
		if not os.path.isdir(os.path.join(self.idir, 'wgs')):
			logger.error('WGS directory not found: %s/wgs', self.idir)
			sys.exit(1)
		if not os.path.isdir(os.path.join(self.idir, 'hic')):
			logger.error('HiC directory not found: %s/hic', self.idir )
			sys.exit(1)
	
		self.sv_fpath = os.path.join(self.idir, 'truth','truth.bed')
		if not os.path.isfile(self.sv_fpath):
			logger.error('Truth bed file not found: %s', self.sv_fpath)
			sys.exit(1)

		self._set_read_counts()		
		self._set_loci()

		for contig in self.loci.contig.unique():
			if not os.path.isfile(os.path.join(self.idir, 'vcf', contig+'.vcf')):
				logger.error('Did not find VCF file for contig: %s/vcf/%s.vcf', self.idir, contig)
				sys.exit(1)
			self._check_for_contig('wgs', contig)	
			self._check_for_contig('hic', contig)	
	
	def _check_for_contig(self, data_src, contig):
		contig_all = os.path.join(self.idir, data_src, contig+'.all.{0}'.format(self.ftype))
		if not os.path.isfile(contig_all):
			logger.error('Missing %s contig : %s', data_src, contig_all)
			sys.exit(1)
		if not os.path.isfile('{0}.{1}'.format(contig_all,self.idx_ftype)):
			logger.error('Missing %s contig index: %s.%s)', data_src, contig_all, self.idx_ftype)
			sys.exit(1)
	
		contig_allele = os.path.join(self.idir, data_src, contig+'.{allele}.'+self.ftype)
		if not (os.path.isfile(contig_allele.format(allele='A')) and os.path.isfile(contig_allele.format(allele='B'))):
			logger.error('Missing allele split for %s on %s...',data_src, contig)

	def _set_read_counts(self):
		""" Check for idxstats file """
		wgs_stat_fpath = os.path.join(self.idir, 'wgs','all.stat')
		if not os.path.isfile(wgs_stat_fpath):
			logger.error('Samtools idxstats file not found: %s/wgs/all.stat', self.idir)
			sys.exit(1)
	
		self.wgs_read_count = pd.read_csv(wgs_stat_fpath, sep='\t', header=None, index_col=0).astype(int)
	
		hic_stat_fpath = os.path.join(self.idir, 'hic','all.stat')
		if not os.path.isfile(hic_stat_fpath):
			logger.error('Samtools idxstats file not found: %s/hic/all.stat', self.idir)
			sys.exit(1)
		self.hic_read_count = pd.read_csv(hic_stat_fpath, sep='\t', header=None, index_col=0).astype(int)
	
	def _set_loci(self):
		with open(self.sv_fpath, 'rb') as f:
			skiprows=0
			for line in f:
				if line.startswith('track'):
					skiprows+=1
				break	
	
		# Loads truth values into dataframe
		self.loci = pd.read_csv(self.sv_fpath, sep='\t', header=None, skiprows=skiprows, index_col=None).rename(columns={0:'contig',1:'start',2:'end',3:'truth'})
	
def main():
	import argparse
	import logging
	import time
	from pkg_resources import Requirement, resource_filename

	parser = argparse.ArgumentParser(description=__doc__)
	
	default_arguments(parser)
	parser.add_argument('version', help='version of training model generated')
	parser.add_argument('ftype', help='file format for reads', choices=['bam','dat'])
	# Sets the input file directories, model, reference, and debug
	parser.add_argument('input_dirs', nargs='+', help='directories containing input hic bam, wgs bam, and idxstat files.')

	parser.add_argument('--seed', type=int, default=None, help='Random initial state for some training procedures')
	parser.add_argument('-P', '--preload-features', dest='preload_feats', action='store_true', default=False, help='Preload features from disk. Fpath is derived from model/version')
	parser.add_argument('-S', '--save-features', dest='save_feats', action='store_true', default=False, help='Save features to disk  for speed. Fpath is derived from model/version')
	parser.add_argument('-k', default=5, help='# of folds in nested cross validation')	

	args = parser.parse_args()

	if args.debug:
		logger.setLevel(logging.DEBUG)
	else:
		logger.setLevel(logging.INFO)

	random_state=int(time.time())%1000000

	model_pkl = resource_filename(Requirement.parse('InPhaDel'), 'models/{m}.{v}.pkl'.format(m=args.model.upper(), v=args.version))

	ref_fpath=args.reference_fasta
	feats_to_disk = resource_filename(Requirement.parse('InPhaDel'), 'models/feats.{v}'.format(v=args.version))

	save_prefix = feats_to_disk if args.save_feats else None
	preload_prefix = feats_to_disk if args.preload_feats else None

	# Process for each contig
	trainer = Trainer(k=args.k)
	#print loci
	for idir in args.input_dirs:	
		fs = PrepTrainingFileStructure(idir, args.ftype)
		fs.load()
		
		for contig in sorted(fs.loci.contig.unique()):
			d = FileStructureDataWithTruth(os.path.join(idir, 'wgs'), os.path.join(idir, 'hic'), fs.sv_fpath, contig, file_fmt=args.ftype, ref_fpath=ref_fpath)

			adaptor = RPKMAdaptor(fs.wgs_read_count.loc[contig,1], fs.hic_read_count.loc[contig,1])
			if not args.preload_feats: 	
				d.fill(adaptor=adaptor)
			trainer.add_data(d)
			logger.info('Added data from contig %s', contig)

	#with warnings.catch_warnings():
	#	warnings.simplefilter("ignore")
	trainer.set_features(save_prefix=save_prefix, preload_prefix=preload_prefix)
	trainer.set_labels(save_prefix=save_prefix, preload_prefix=preload_prefix)

	if args.model=='svm':
			m = SVMModel(model_pkl, args.seed)

	acc = trainer.train(m)

if __name__ == '__main__':
	main()
