from sklearn import svm, ensemble
from sklearn.cross_validation import StratifiedKFold, cross_val_score
from sklearn.grid_search import GridSearchCV
from sklearn import metrics

from svphase.vis import rect_comp
from svphase.learn.classify import realistic_classes_simple, realistic_classes_wgs2x, realistic_classes_hic2x, realistic_classes_sim_hm_mix, PRECISION
from svphase.learn.cov import DataAdaptor, RPKMAdaptor, SimulationDelData, RealDelData, Chr20RealDelData
from svphase.inphadel import Predictor

from operator import itemgetter
import pickle
from itertools import chain, compress
from collections import defaultdict
import numpy as na
import pandas as pd
import time
import os

#PLOT_DIR = '/home/anand/Dropbox/SV_Detection/Results/working'
#PLOT_DIR = '/home/anand/Dropbox/SV_Detection/Results/working'
#PLOT_DIR = '/home/anand/Dropbox/SVPhasePaper/fig'
PLOT_DIR = '/home/anand/Projects/InPhaDelPaper/fig/adjust/'
#PLOT_DIR = None

#DATA_PREFIX = '/media/T02/data/hic/'
DATA_PREFIX = '/home/anand/T02/data/hic/'
FONT = {'family': 'Liberation Sans', 'size': 14}

pd.set_option('display.max_columns',500)
pd.set_option('display.max_rows',1000)
#pd.set_printoptions(max_columns=500, max_rows=1000)

#def FPATH_TO_NAME(prefix_fpath, feature_typeset):
#	if feature_typeset is None:
#		feature_typeset = ''
#	return '{0}.{1}'.format(feature_typeset, os.path.basename(prefix_fpath).replace('sim.del-chr19-20.5','')[:-4])

def FPATH_TO_NAME(prefix_fpath, feature_typeset=None, eval_name=''):
	if feature_typeset is None:
		feature_typeset = ''
	name = os.path.basename(prefix_fpath).split('.')
	if name[-1]=='pkl':
		name = name[:-1]

	path = '' if PLOT_DIR is None else PLOT_DIR
	
	f_s = None
	ml = None
	if len(name)==3:
		model,f = name[1:]
	elif len(name)==5:
		model,f,nfolds,ml = name[1:]
	elif len(name)==6:
		model,f,f_s,nfolds,ml = name[1:]
	else:
		print name
		raise

	path = os.path.join(path, '_'.join([model,f]), eval_name)
	try:
		os.makedirs(path)
	except OSError:
		pass

	if ml is None:
		return path
	if f_s is None:
		return os.path.join(path, '-'.join([nfolds,ml]))
	else:
		return os.path.join(path, '-'.join([nfolds,f_s,ml]))
		
def plot_predict_multi(clf_prefix_fpath, classifier, evaluator=None):
	import matplotlib as mpl
	from matplotlib import pyplot as plt
	from matplotlib import gridspec as gridspec
	mpl.rc('font', **FONT)

	if evaluator is None:
		evaluator = EvalSim

	if '.' in classifier:
		cv,classifier = classifier.split('.')
	else:
		cv = None

	if '-' in classifier:
		feature_typeset, classifier = classifier.split('-')
	else:
		feature_typeset = None

	assert classifier in set(['SVM', 'RandomForest', 'GNB', 'KN'])
	#categories = ['', '.only', '.wgs_only', '.hic_only']
	#titles = ['All', 'WGS', "WGS + WGS Allele Filter", "WGS + HiC Allele Filter"]

	categories = ['', '.wgs_only', '.hic_only']
	titles = ['All', "WGS + WGS Allele Filter", "WGS + HiC Allele Filter"]

	fig = plt.figure(figsize=(8,6.8),facecolor='white', edgecolor='none')
	gs = gridspec.GridSpec(len(categories),1)
	gs.update(left=0.08, right=0.92, hspace=0.25)
	for i, cat in enumerate(categories):
		ax = fig.add_subplot(gs[i,0])

		clf_fpath = "{0}{1}.".format(clf_prefix_fpath, cat)
		if not cv is None:
			clf_fpath += "%s."%cv
		if not feature_typeset is None:
			clf_fpath += "%s-"%feature_typeset
		clf_fpath += "%s.pkl"%classifier
		
		c = cat[1:] if cat!='' else None
		t = evaluator(clf_fpath, combines=c, feature_typeset=feature_typeset)
		t.plot(ax=ax, min_accuracy=0.4)
		if i ==0:
			ax.set_xlabel('% correct classifications')
			ax.set_title("%s - %s"%(classifier, titles[i]))
		else:
			ax.set_title(titles[i])
			ax.set_ylabel('Accuracy')
	gs.tight_layout(fig)
	if PLOT_DIR is None:
		plt.show()
	else:
		if cv is None:
			plt.savefig("{0}/{1}-error_all.pdf".format(FPATH_TO_NAME(clf_prefix_fpath, eval_name=t.name), classifier))
		else:
			plt.savefig("{0}/{2}-{1}-error_all.pdf".format(FPATH_TO_NAME(clf_prefix_fpath, eval_name=t.name), classifier, cv))

def plot_predict_comp(clf_hic_only_fpath, clf_wgs_only_fpath, evaluator=None):
	import matplotlib as mpl	
	from matplotlib import pyplot as plt
	if evaluator is None:
		evaluator = EvalSim
	mpl.rc('font', **FONT)
	
	hic_ev = evaluator(clf_hic_only_fpath, combines='hic_only')
	hic_df, idx_to_class = hic_ev.plot(ax='')
	wgs_ev = evaluator(clf_wgs_only_fpath, combines='wgs_only')
	wgs_df, idx_to_class = wgs_ev.plot(ax='')
	# note truth now has nan_idx being class 4
	nan_idx = [k for k,v in idx_to_class.items() if v=='nan'][0]

	truth = hic_df['truth']
	wgs_additional_idx = na.logical_and(hic_df['truth']!=wgs_df['truth'], truth==nan_idx)
	truth.ix[wgs_additional_idx] = wgs_df.ix[wgs_additional_idx, 'truth']

	#quick_check = pd.DataFrame(truth)
	#quick_check['hic'] = hic_df['truth']
	#quick_check['wgs'] = wgs_df['truth']
	#print quick_check.T.ix[:,::2]
	#assert (hic_df['truth']==wgs_df['truth']).all()

	df = pd.DataFrame([truth, hic_df['correct'], wgs_df['correct']], index=['truth', 'hic_only', 'wgs_only']).T
	df['both'] = ((df['hic_only']+df['wgs_only'])==2).astype(int)
	grouped = df.groupby('truth')
	agg_count = grouped.aggregate(na.sum)
	both_count = agg_count['both'].rename(idx_to_class)
	hic_only_count = (agg_count['hic_only']-agg_count['both']).rename(idx_to_class)
	wgs_only_count = (agg_count['wgs_only']-agg_count['both']).rename(idx_to_class)
	grouped_size = grouped.size().rename(idx_to_class)
	print agg_count
	df_count = pd.DataFrame([hic_only_count, both_count, wgs_only_count], index=['a', 'both', 'b'], columns=both_count.index).T
	print df_count
	print "Columns:\n", df_count.sum(axis=0)/float(df_count.sum().sum())
	print "Classes:\n", df_count.sum(axis=1)
	print "Total:\n", df_count.sum().sum()
	#phased_combined = (df_count.T[['pA',]].values+df_count.T[['pB',]].values)
	#print phased_combined
	#print phased_combined[0]/float(na.sum(phased_combined))
	#print phased_combined[2]/float(na.sum(phased_combined))

	fig = plt.figure(figsize=(7,4), facecolor='white', edgecolor='none')
	ax = fig.add_subplot(111)
	ax.set_title(clf_hic_only_fpath.split('.')[-2])
	rect_comp.model_diff_ax(ax, df_count,grouped_size, text_ylabels=['HiC-Allele\nonly', 'Both', 'WGS-Allele\nonly'])	
	#plt.tight_layout()
	fig.subplots_adjust(left=0.10, bottom=0.08, right=0.85, top=0.87)
	if PLOT_DIR is None:
		plt.show()
	else:
		hic_only = FPATH_TO_NAME(clf_hic_only_fpath, eval_name=hic_ev.name)
		wgs_only = FPATH_TO_NAME(clf_wgs_only_fpath, eval_name=hic_ev.name)
		plt.savefig("{0}-comp-{1}.pdf".format(hic_only, os.path.basename(wgs_only)))

class EvalPrediction(object):
	def __init__(self, clf_fpath, combines, feature_typeset=None,	score_flag=False, name=''):
		self.pclf = realistic_classes_simple(False)
		self.predr = Predictor(self.pclf, clf_fpath, feature_typeset=feature_typeset)
		self.name = name
		self.predr_name = FPATH_TO_NAME(clf_fpath, feature_typeset, self.name)
		print self.predr_name
		#self.sv_fpath = None
		self.pred = None
		self.lpred = None
		self.truth = None
		self.combines = combines
		self._predict()
		if not score_flag:
			self._truth()
	def _predict(self):
		pass
	def _truth(self):
		pass

	def _get_inputs(self):
		 return pd.DataFrame([self.predr.class_int[i] for i in self.truth], columns=['truth',])
	def _save_preds(self):
		idx_to_class = dict(self.predr.class_int.items())
		pl_df = self._get_inputs()
		#pl_df = self._get_inputs().drop_duplicates(['chr','start','end'])
		if self.pred is None:
			assert not self.lpred is None
			d = pd.DataFrame(self.lpred, columns=map(itemgetter(1), sorted(self.predr.class_int.items())))
			d = d.fillna(-100)
			pl_df['pred'] = d.idxmax(axis=1)
			pl_df.ix[na.logical_not(self.nonzero_idx),'pred'] = 'nan'
			pl_df = pd.concat([pl_df,d],axis=1)		 
		else:
			pl_df['pred'] = map(idx_to_class.get,self.pred)
			pl_df.ix[na.logical_not(self.nonzero_idx),'pred'] = 'nan'
		pl_df.to_csv("{0}.txt".format(self.predr_name), sep='\t', header=True, index=False)

	def _set_nonzero(self, pclf):
		self.nonzero_idx = pclf.get_nonzero(self.combines)
	def write_log_proba(self, csv_fpath=None):
		if csv_fpath is None:
			print self.pred_df.T.ix[:,::2]
		else:
			self.pred_df.to_csv(csv_fpath)
	def get_accuracy(self, classes, with_nan_flag=False):
		assert not self.truth is None
		assert not self.pred is None
		assert self.pred.size == self.truth.size
		nan_idx = int(na.max(self.truth)+1)
		
		idx_to_class = dict(self.predr.class_int.items())
		idx_to_class[nan_idx]='nan'

		df_compare = pd.DataFrame([self.truth,self.pred], index=['truth', 'pred']).T
		df_compare.ix[na.logical_not(self.nonzero_idx),'pred'] = nan_idx+1
		df_compare.ix[na.logical_not(self.nonzero_idx),'truth'] = nan_idx

		#df_compare = pd.concat([df_compare, nan_df], ignore_index=True)
		df_compare['correct'] = (df_compare['truth']==df_compare['pred']).astype(int)
		grouped = df_compare.groupby('truth')
		correct_count = grouped.aggregate(na.sum)['correct'].rename(idx_to_class)
		grouped_size = grouped.size().rename(idx_to_class)
		
		#print correct_count
		#print grouped_size
		#print "Fraction", grouped_size[['nan',]]/float(na.sum(grouped_size))
		#print "# of used examples ", (na.sum(grouped_size)-grouped_size[['nan',]]).values
		gclasses = list(classes)
		if with_nan_flag:
			gclasses.append('nan')
		acc =	na.sum(correct_count[classes])/float(na.sum(grouped_size[gclasses]))
		return acc 

	def plot(self, ax=None, min_accuracy=0.1):
		try:
			assert not self.pred is None
		except AssertionError:
			print >>sys.stderr, "Plot requires non log probabilities"
			return
		assert not self.truth is None
		assert self.pred.size == self.truth.size
		nan_idx = int(na.max(self.truth)+1)
		
		idx_to_class = dict(self.predr.class_int.items())
		idx_to_class[nan_idx]='nan'
		print idx_to_class
		het_hack_nan = na.logical_or(self.truth==1, self.truth==3)
		print 'HackNAN', na.sum(na.logical_and(na.logical_not(self.nonzero_idx), het_hack_nan))
		df_compare = pd.DataFrame([self.truth,self.pred], index=['truth', 'pred']).T

		df_compare.ix[na.logical_not(self.nonzero_idx),'pred'] = nan_idx+1
		df_compare.ix[na.logical_not(self.nonzero_idx),'truth'] = nan_idx

		#df_compare = pd.concat([df_compare, nan_df], ignore_index=True)
		df_compare['correct'] = (df_compare['truth']==df_compare['pred']).astype(int)
		grouped = df_compare.groupby('truth')
		correct_count = grouped.aggregate(na.sum)['correct'].rename(idx_to_class)
		grouped_size = grouped.size().rename(idx_to_class)
		
		#print correct_count
		print "Fraction", grouped_size[['nan',]]/float(na.sum(grouped_size))
		print "# of used examples ", (na.sum(grouped_size)-grouped_size[['nan',]]).values
		print "Classification Accuracy", na.sum(correct_count[['pA','pB','inc','homA']])/float(na.sum(grouped_size[['pA','pB','inc','homA']]))
	 
		print correct_count.astype(float)/grouped_size.astype(float)
		print correct_count.astype(float)
		print grouped_size
		#return df_compare, idx_to_class
		if ax is None:
			import matplotlib as mpl
			from matplotlib import pyplot as plt
			fig = plt.figure(figsize=(6,3),facecolor='white', edgecolor='none')
			mpl.rc('font', **FONT)
			ax = fig.add_subplot(111)
			#ax.set_title(self.predr_name)
			#ax.set_title('.'.join(clf_fpath.split('.')[1:-1]))
			rect_comp.error_ax(ax, correct_count, grouped_size, min_accuracy=min_accuracy, background_flag=False)
			ax.set_ylabel('Accuracy')
			fig.tight_layout()
			if PLOT_DIR is None:
				plt.show()
			else:
				plt.savefig("{0}-error.pdf".format(self.predr_name))

		elif ax!='':
			rect_comp.error_ax(ax, correct_count, grouped_size) 
		#print correct_count.rename(idx_to_class)
		#print grouped_size.rename(idx_to_class)
		return df_compare, idx_to_class

class EvalReal(EvalPrediction):
	def __init__(self, clf_fpath, combines, 
							 #sv_fpath='/home/anand/Projects/assembly/data/gm12878/sid-del-truth-man.bed',
							 #sv_fpath='/home/anand/Projects/assembly/data/gm12878/hapmap3_r2_b36_fwd-quake-truth.bed',
							 #sv_fpath='/home/anand/Projects/assembly/data/mills_2010-08-table4.CEU1kb.NA12878.bed',
							 sv_fpath='/home/anand/Projects/assembly/data/mills_2010-08-table4.CEU1kb.NA12878.anyp.bed',
							 #sv_fpath='/home/anand/Projects/assembly/data/mills_2010-08-table4.CEU1kb.NA12878.any.bed',
							 #sv_fpath='/home/anand/Projects/assembly/data/mills_na12878-del-gs.NA12878.any.bed',
							 #sv_fpath='/home/anand/Projects/assembly/data/mills_na12878-del-gs.NA12878.anyp.bed',
							 #sv_fpath='/home/anand/Projects/assembly/data/mills_na12878-del-gs.NA12878.max.bed',
							 #sv_fpath='/home/anand/Projects/assembly/data/mills_2010-08-table4.CEU1kb.NA12878.max.bed',
							 #sv_fpath='/home/anand/Projects/assembly/data/gm12878/hapmap3_r2_b36_fwd-truth.bed',
							 #sv_fpath='/home/anand/Projects/assembly/data/gm12878/hapmap3_r2_b36_fwd.dels_only.bed',
							 score_flag=False, feature_typeset=None, name=None):
		if name is None:
			name = os.path.basename(sv_fpath)[:5]
		self.sv_fpath = sv_fpath 
		EvalPrediction.__init__(self, clf_fpath, combines, score_flag=score_flag, feature_typeset=feature_typeset, name=name)
		self._set_nonzero(self.predr)
		self._save_preds()

	def _get_inputs(self):
		t = pd.read_csv(self.sv_fpath, sep='\t', header=None, names=['chr','start','end','truth'])
		print t.shape
		return t

	def _get_mapped_read_count(self, fpath):
		t = {}
		with open(fpath, 'rb') as f:
			for line in f:
				tokens = line.strip().split('\t')
				t[tokens[0]] = int(tokens[2])
		return t
	def _predict(self, log_proba_flag=False):
		from svphase.utils.common import chrom_name_map 

		prefix_dir = DATA_PREFIX + 'splitfiles'
		wgs_fpath = DATA_PREFIX + 'unsplit_files/wgs_bamfiles_unsplit/NA12878.hiseq.wgs.bwa.recal.bam'
		hic_stat_fpath = DATA_PREFIX + 'unsplit_files/hic_bamfiles_unsplit/chr.all.stat'
		wgs_stat_fpath = DATA_PREFIX + 'unsplit_files/wgs_bamfiles_unsplit/NA12878.hiseq.wgs.bwa.recal.stat'
		wgs_read_count = self._get_mapped_read_count(wgs_stat_fpath)
		hic_read_count = self._get_mapped_read_count(hic_stat_fpath)
 
		avail_chr = [c for c in os.listdir(prefix_dir) if c in chrom_name_map.viewkeys()]
		
		with open(self.sv_fpath, 'rb') as f:
			contigs_with_loci = set([t.split('\t')[0] for t in f if t.startswith('chr')])	 
		#print wgs_read_count, hic_read_count 
		#adaptor = DataAdaptor()
		adaptor = RPKMAdaptor(wgs_read_count['chr19']+wgs_read_count['chr20'], hic_read_count['chr19']+hic_read_count['chr20'])

		for contig in sorted(avail_chr, key=chrom_name_map.get):
			if not contig in contigs_with_loci: continue
			#if not contig=='chr20': continue
			#print wgs_read_count[contig], hic_read_count[contig]
			rdata = RealDelData(wgs_fpath, prefix_dir, self.sv_fpath, contig)
			#rdata.fill(adaptor=RPKMAdaptor(wgs_read_count[contig], hic_read_count[contig]), load_flag=True)
			rdata.fill(adaptor=adaptor, load_flag=True)
			self.predr.add_data(rdata)
		if log_proba_flag:
			self.lpred = self.predr.predict_log_proba(combines=self.combines, nonzero_flag=False)
			#self.pred_df = self.predr.get_pred_df(self.pred)
		else:
			self.pred = self.predr.predict(combines=self.combines, format_flag=False, nonzero_flag=False)
			#self.pred_s = self.predr.get_pred_series(self.pred)

 
	def _truth(self):

		class_str_dict = defaultdict(lambda:-1)
		for k,v in self.pclf.class_str.iteritems():
			class_str_dict[k] = v

		self.truth = na.array(map(class_str_dict.get, self._get_inputs()['truth']))

class EvalScoreReal(EvalReal):
	def __init__(self, clf_fpath, combines, feature_typeset=None, sv_fpath=None):
		if sv_fpath is None:
			#sv_fpath = '/home/anand/Projects/assembly/data/trios/gm12878.del1kb.svforest.bed'
			#sv_fpath='/home/anand/Projects/assembly/data/gm12878/sid-del-truth-man.bed'
			sv_fpath='/home/anand/Projects/assembly/data/gm12878/hapmap3_r2_b36_fwd-truth.bed'
			#sv_fpath='/home/anand/Projects/assembly/data/trio.2010_06.deletions.sites.1kb.bed'
			#sv_fpath='/home/anand/Projects/assembly/data/mills_2010-08-table4.CEU1kb.NA12878.anyp.bed'
			#sv_fpath='/home/anand/Projects/assembly/data/mills_na12878-del-gs.NA12878.anyp.bed'
			#sv_fpath='/home/anand/Projects/assembly/data/gm12878/hapmap3_r2_b36_fwd.dels_only.bed'
			#sv_fpath='/home/anand/Projects/assembly/data/mills_2010-08-table4.CEU1kb.NA12878.any.bed'
			#sv_fpath='/home/anand/Projects/assembly/data/mills_2010-08-table4.CEU1kb.hets_only.bed'
			#sv_fpath='/home/anand/Projects/assembly/data/mills_na12878-del-gs.NA12878.any.bed'
		name = os.path.basename(sv_fpath)[:5] + "-score-h"
		EvalReal.__init__(self, clf_fpath, combines, feature_typeset=feature_typeset, score_flag=False, sv_fpath=sv_fpath, name=name)
		#EvalReal.__init__(self, clf_fpath, combines, score_flag=True)
	def _predict(self):
		EvalReal._predict(self, log_proba_flag=True)
		#print self.pred_df

class EvalSim(EvalPrediction):
	def __init__(self, clf_fpath, combines, feature_typeset=None, score_flag=False, name='simul'):
		EvalPrediction.__init__(self, clf_fpath, combines, feature_typeset=feature_typeset, score_flag=score_flag, name=name)
		self._save_preds()

	def _predict(self, log_proba_flag=False):		 
		for data in self.pclf.data:
			self.predr.add_data(data)
		if log_proba_flag:
			self.lpred = self.predr.predict_log_proba(combines=self.combines, nonzero_flag=False)
		else:
			self.pred = self.predr.predict(combines=self.combines, format_flag=False, nonzero_flag=False)
		 
	def _truth(self):
		pclf2 = realistic_classes_simple(False)
		pclf2.get_features(combines=self.combines, nonzero_flag=False)
		self._set_nonzero(pclf2)
		self.truth = pclf2.labels

class EvalScoreSim(EvalSim):
	def __init__(self, clf_fpath, combines, feature_typeset=None):
		EvalSim.__init__(self, clf_fpath, combines, feature_typeset=feature_typeset, score_flag=False, name='simul-score')
	def _predict(self):
		EvalSim._predict(self, log_proba_flag=True)
		#print self.pred_df

class EvalSimWGS2X(EvalSim):
	def __init__(self, clf_fpath, combines, feature_typeset=None):
		EvalSim.__init__(self,clf_fpath, combines, feature_typeset=feature_typeset, name='sim-wgs2x')
	def _predict(self):
		for data in realistic_classes_wgs2x(False).data:
			self.predr.add_data(data)
		self.pred = self.predr.predict(combines=self.combines, format_flag=False, nonzero_flag=False)
	def _truth(self):
		pclf2 = realistic_classes_wgs2x(False)
		pclf2.get_features(combines=self.combines, nonzero_flag=False)
		self._set_nonzero(pclf2)
		self.truth = pclf2.labels

class EvalSimHiC2X(EvalSim):
	def __init__(self, clf_fpath, combines, feature_typeset=None):
		EvalSim.__init__(self,clf_fpath, combines, feature_typeset=feature_typeset, name='sim-hic2x')
	def _predict(self):
		for data in realistic_classes_hic2x(False).data:
			self.predr.add_data(data)
		self.pred = self.predr.predict(combines=self.combines, format_flag=False, nonzero_flag=False)
	def _truth(self):
		pclf2 = realistic_classes_hic2x(False)
		pclf2.get_features(combines=self.combines, nonzero_flag=False)
		self._set_nonzero(pclf2)
		self.truth = pclf2.labels

class EvalMix(EvalPrediction):
	def __init__(self, clf_fpath, combines, feature_typeset=None, name='sim-hm'):
		EvalPrediction.__init__(self, clf_fpath, combines, feature_typeset=feature_typeset, name=name)
	def _predict(self,log_proba_flag=False):		 
		for data in realistic_classes_sim_hm_mix(False).data:
			self.predr.add_data(data)
		if log_proba_flag:
			self.lpred = self.predr.predict_log_proba(combines=self.combines, nonzero_flag=False)
		else:
			self.pred = self.predr.predict(combines=self.combines, format_flag=False, nonzero_flag=False)
	def _truth(self):
		pclf2 = realistic_classes_sim_hm_mix(False)
		pclf2.get_features(combines=self.combines, nonzero_flag=False)
		self._set_nonzero(pclf2)
		self.truth = pclf2.labels

class EvalScoreMix(EvalMix):
	def __init__(self, clf_fpath, combines, feature_typeset=None):
		EvalMix.__init__(self, clf_fpath, combines, feature_typeset=feature_typeset, name='sim-hm-score')
		self._save_preds()
	def _predict(self):
		EvalMix._predict(self, log_proba_flag=True)

def _feature_typeset_from_clf_fpath(fpath):
	tokens = clf_fpath.split('.')[-2].split('-')
	if len(tokens)==1:
		return None
	else:
		assert tokens[0] in ['simple_sum', 'binary_allele']
		return tokens[0]

if __name__=='__main__':
	import sys, os
	evaluator = EvalMix
	#evaluator = EvalScoreMix
	#evaluator = EvalSim
	#evaluator = EvalScoreSim
	#evaluator = EvalReal
	#evaluator = EvalScoreReal

	#evaluator = EvalSimWGS2X
	#evaluator = EvalSimHiC2X

	clf_fpath = sys.argv[1]
	if clf_fpath.endswith('.pkl'):
		try:
			assert os.path.isfile(clf_fpath)
		except AssertionError:
			raise BaseException("Classifier fpath is not correct")	
		if len(sys.argv)==2:
			out_fpath = '{0}.output.csv'.format(os.path.basename(clf_fpath[:-4]))
			print out_fpath

			w = evaluator(clf_fpath, clf_fpath.split('.')[-4], feature_typeset=_feature_typeset_from_clf_fpath(clf_fpath))
			w.plot(ax=None)
			#w.write_log_proba(csv_fpath=out_fpath)
			#plot_simulation_predict(clf_fpath, ax='')
		else:
			clf_wgs_only_fpath = sys.argv[2]
			plot_predict_comp(clf_fpath, clf_wgs_only_fpath, evaluator)
	else:
		clf_prefix_fpath, classifier = clf_fpath, sys.argv[2]
		plot_predict_multi(clf_prefix_fpath, classifier, evaluator)
