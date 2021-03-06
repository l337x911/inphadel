"""InPhaDel: Genotypes and phase deletions on a single chromosome using a specific classification model

Trains models for phasing deletions using underlying WGS+HiC data 

"""
import sys
import os
import pickle
import pandas as pd
import numpy as np
import warnings
from itertools import izip

from sklearn import svm
from sklearn import ensemble
from sklearn.cross_validation import StratifiedKFold
from sklearn.grid_search import GridSearchCV
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import accuracy_score

from svphase.utils.common import logger
from svphase.utils.config import RANDOM_STATE,DATA_PREFIX
from svphase.learn.cov import FileStructureDataWithTruth, RPKMAdaptor
from svphase.learn.evaluation import Evaluation
from svphase.learn.features import HicOnlySubset, WgsOnlySubset
from svphase.inphadel import default_arguments

class Model(object):
	def __init__(self, model_pkl, random_state):
		self.pkl = model_pkl
		self.random_state = random_state
		self.clf = None
		self.params = {}
	def clf_stats(self):
		pass

class SVMModel(Model):
	def __init__(self, model_pkl, random_state):
		Model.__init__(self, model_pkl, random_state)
		self.clf = svm.SVC(kernel='linear', probability=True, random_state=self.random_state)
		self.params = {'C':[.1,1,10,100]}

class RFModel(Model):
	def __init__(self, model_pkl, random_state):
		Model.__init__(self, model_pkl, random_state)
		self.clf = ensemble.RandomForestClassifier(oob_score=True, random_state=self.random_state)
		self.params = {'n_estimators':[10,20,50,100], 'max_depth':[2,5,10,20]}
	def clf_stats(self):
		logger.info('RFModel: OOB_Score {0:0.4f}'.format(self.clf.oob_score_))
		

class KNNModel(Model):
	def __init__(self, model_pkl, random_state):
		Model.__init__(self, model_pkl, random_state)
		self.clf = KNeighborsClassifier(n_neighbors=5, weights='uniform', algorithm='brute')
		self.params = {'n_neighbors':[2,4,8,16,32]}


class Trainer(Evaluation):
	def __init__(self, k=5, feature_subset=None): 
		Evaluation.__init__(self, feature_subset=feature_subset)
		self.k = k
		self.inner_models = []
		self.inner_accuracies = []
		self._cpred_index =map(lambda x:'cpred:'+x,self.label_obj.classes)
		self._tpred_index =map(lambda x:'tpred:'+x,self.label_obj.classes)


	def check_stable(self, scores, common_params_df):
		thresh = 0.10
		deviation = max(scores)-min(scores)

		if deviation>thresh:
			logger.warning('Model test accuracies deviate more than {thresh:0.1f}%, Deviation {dev:0.4f}'.format(thresh=thresh*100, dev=deviation))

		if len(common_params_df.index)>1:
			logger.warning('Model had unstable parameters\n{len:s}'.format(len=common_params_df))

	def _to_series(self, outer_fold, inner_fold, test_accuracy, correct_preds, total_preds, params):

		if correct_preds is None:
			c = pd.Series([0,]*len(self.label_obj.classes), index=self._cpred_index, dtype=int)
		else:
			c = pd.Series(correct_preds, index=self._cpred_index, dtype=int)
		if total_preds is None:
			t = pd.Series([0,]*len(self.label_obj.classes), index=self._tpred_index, dtype=int)
		else:
			t = pd.Series(total_preds, index=self._tpred_index, dtype=int)

		return pd.concat([pd.Series([outer_fold, inner_fold, test_accuracy], index=['outer_fold','inner_fold', 'test_accuracy']), c,t,pd.Series(params, dtype=int)])

	def _correct_preds_per_class(self, labels, preds):
		# Expects labels, and preds to be in ints
		cpreds = [0,]*len(self.label_obj.classes)
		tpreds = [0,]*len(self.label_obj.classes)
		
		for l,p in izip(labels, preds):
			tpreds[l] += 1
			cpreds[l] += int(p==l)

		return cpreds, tpreds

	def train(self, model):
		data = self.feats.get_nonzero_features()
		data = data.fillna(0).astype(np.float64).values
		assert np.isfinite(data).all()
		labels_str = np.array(self.labels)[self.feats.get_nonzero()]
		labels_int = map(self.label_obj.str_to_int_dict.get, labels_str)
		logger.debug(labels_str)
		logger.debug(labels_int)
		labels = np.array(labels_int, dtype=int)
		logger.debug(labels)
		# Applies Nested Cross Validation on data + labels
		#outer_skf = StratifiedKFold(labels, n_folds=self.k)
		outer_skf = StratifiedKFold(labels, n_folds=self.k, shuffle=True, random_state=model.random_state)

		scores = []
	
		skf_stats = []	
		logger.debug('Data shape %s, label shape %s', data.shape, labels.shape)
		for fold, (train_index, test_index) in enumerate(outer_skf):
			# Perform a new inner cross validation
			logger.debug('Train shape %s, test shape %s', train_index.shape, test_index.shape)

			train_data, test_data = data[train_index,:], data[test_index,:]
			train_label, test_label = labels[train_index], labels[test_index]
			logger.debug('Train Data shape %s, Train Label shape %s', train_data.shape, train_label.shape)
			inner_skf = StratifiedKFold(train_label, n_folds=self.k, shuffle=True, random_state=model.random_state)
			grid_clf = GridSearchCV(model.clf, param_grid=model.params, scoring='accuracy', n_jobs=-1, cv=inner_skf, refit=True, verbose=1)
			grid_clf.fit(train_data, train_label)
		
			outer_predict = grid_clf.best_estimator_.predict(test_data)
			test_accuracy = accuracy_score(test_label, outer_predict )
			scores.append(test_accuracy)
			cpreds, tpreds = self._correct_preds_per_class(test_label, outer_predict)			
			#logger.debug('OUTER: %0.4f test accuracy on fold %d', test_accuracy, fold)
			logger.info('OUTER: %0.4f test accuracy with params %s on fold %d', test_accuracy, grid_clf.best_params_, fold)
			skf_stats.append(self._to_series(fold, None, test_accuracy, cpreds, tpreds, grid_clf.best_params_))
			for pt in grid_clf.grid_scores_:
				logger.info(' INNER: %.4f avg accuracy with params %s and scores %s', pt.mean_validation_score, pt.parameters, ','.join(map('{0:0.4f}'.format,pt.cv_validation_scores)))
				for inner_fold, cv_score in enumerate(pt.cv_validation_scores):
					skf_stats.append(self._to_series(fold, inner_fold, cv_score, None, None, pt.parameters))

		skf_stats_df =  pd.concat(skf_stats, axis=1).T
		param_cols = skf_stats_df.columns[3+2*len(self.label_obj.classes):]
		pred_cols = skf_stats_df.columns[3:3+2*len(self.label_obj.classes)]
		skf_stats_df[pred_cols] = skf_stats_df[pred_cols].astype(int)
		skf_stats_df[param_cols] = skf_stats_df[param_cols].astype(int)
		#skf_stats_df[['outer_fold','inner_fold']] = skf_stats_df[['outer_fold','inner_fold']].astype(int)

		skf_stats_df.to_csv(model.pkl + '.stat', sep='\t', index=False, header=True, float_format="%0.4f") 
		# Assumes skf_stats_df columns are outer_fold, inner_fold, accuracy, cpreds:,tpreds:, param1, param2, ...
		outer_params_df = skf_stats_df.ix[skf_stats_df['inner_fold'].isnull(),param_cols]
		common_params = outer_params_df.groupby(list(outer_params_df.columns)).apply(len)
		self.check_stable(scores, common_params)
		argmax = common_params.argmax()
		try:
			iter(argmax)
		except TypeError:
			argmax = [argmax,]		
	
		outer_best_params = dict(zip(common_params.index.names, argmax))
		logger.info('Final Params %s', outer_best_params)
		# Final Fit!
		model.clf.set_params(**outer_best_params)
		model.clf.fit(data, labels)
		model.clf_stats()

		final_accuracy = accuracy_score(labels, model.clf.predict(data))
		logger.info('Final accuracy: %0.4f with params %s', final_accuracy, outer_best_params)
		with open(model.pkl, 'wb') as f:
			pickle.dump(model.clf, f)

class PrepTrainingFileStructure(object):
	def __init__(self, idir):
		self.idir = idir
		self.sv_fpath = None
		self.wgs_read_count = None	
		self.hic_read_count = None	
		self.loci = None
		self.ref_fpath = None
		self.ftype = None
		self.idx_ftype = None

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

		self.ref_fpath = os.path.join(self.idir, 'reference.fa')
		if not os.path.isfile(self.ref_fpath):
			logger.error('Reference file not found: %s', self.ref_fpath)
			sys.exit(1)
		if not os.path.isfile(self.ref_fpath + '.fai'):
			logger.error('Reference Index file not found: %s', self.ref_fpath)
			sys.exit(1)
		
		if os.path.isfile(os.path.join(self.idir, 'dat')):
			self.ftype = 'dat'
		elif os.path.isfile(os.path.join(self.idir, 'bam')):
			self.ftype = 'bam'
		else:
			logger.error('No read filetype file found. idir should contain a bam or dat file')
			sys.exit(1)

		if self.ftype == 'bam':
			self.idx_ftype = 'bai'
		elif self.ftype == 'dat':
			self.idx_ftype = 'npz'
		else:
			logger.error('No valid read filetype extensions: %s', self.ftype)
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
	#parser.add_argument('ftype', help='file format for reads', choices=['bam','dat'])
	# Sets the input file directories, model, reference, and debug
	parser.add_argument('input_dirs', nargs='+', help='directories containing input hic bam, wgs bam, and idxstat files.')

	parser.add_argument('--seed', type=int, default=None, help='Random initial state for some training procedures')
	parser.add_argument('-k', type=int, default=5, help='# of folds in nested cross validation')	
	parser.add_argument('-C', '--check-features', dest='check_feats', action='store_true', default=False, help='Outputs features for training data')	

	args = parser.parse_args()

	if args.debug:
		logger.setLevel(logging.DEBUG)
	else:
		logger.setLevel(logging.INFO)

	random_state=int(time.time())%1000000

	model_pkl = resource_filename(Requirement.parse('InPhaDel'), 'models/{m}.{v}.pkl'.format(m=args.model, v=args.version))

	#ref_fpath=args.reference_fasta
	feats_to_disk = resource_filename(Requirement.parse('InPhaDel'), 'models/feats.{v}'.format(v=args.version))

	save_prefix = feats_to_disk if args.save_feats else None
	preload_prefix = feats_to_disk if args.preload_feats else None

	if args.subset is None:
		feature_subset = None
	elif args.subset=='hic-only':
		feature_subset = HicOnlySubset()
	elif args.subset=='wgs-only':
		feature_subset = WgsOnlySubset()

	# Process for each contig
	trainer = Trainer(k=args.k, feature_subset=feature_subset)
	#print loci
	for idir in args.input_dirs:	
		fs = PrepTrainingFileStructure(idir)
		fs.load()
		
		for contig in sorted(fs.loci.contig.unique()):
			d = FileStructureDataWithTruth(os.path.join(idir, 'wgs'), os.path.join(idir, 'hic'), fs.sv_fpath, contig, file_fmt=fs.ftype, ref_fpath=fs.ref_fpath)

			adaptor = RPKMAdaptor(fs.wgs_read_count.loc[contig,1], fs.hic_read_count.loc[contig,1])
			if not args.preload_feats: 	
				d.fill(adaptor=adaptor)
			trainer.add_data(d)
			logger.info('Added data from contig %s', contig)

	#with warnings.catch_warnings():
	#	warnings.simplefilter("ignore")
	trainer.set_features(save_prefix=save_prefix, preload_prefix=preload_prefix, simple_sum_flag=args.simple_sum)
	trainer.set_labels(save_prefix=save_prefix, preload_prefix=preload_prefix)

	if args.model=='svm':
		M = SVMModel
	elif args.model=='rf':
		M = RFModel
	elif args.model=='knn':
		M = KNNModel	
	if args.check_feats:
		trainer.manual_check()
	else:
		acc = trainer.train(M(model_pkl, args.seed))

if __name__ == '__main__':
	main()
