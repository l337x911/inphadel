from sklearn import svm, ensemble
from sklearn.cross_validation import StratifiedKFold, cross_val_score
from sklearn.grid_search import GridSearchCV
from sklearn.naive_bayes import GaussianNB
from sklearn.neighbors import KNeighborsClassifier
from sklearn.decomposition import PCA


from svphase.utils.config import RANDOM_STATE,DATA_PREFIX
from svphase.learn.cov import SimulationDelData, RealDelData, RPKMAdaptor

import pickle
from itertools import chain
import numpy as na
import pandas as pd
import time


class SVPhaseClassification(object):
	def __init__(self, feature_typeset=None, nfolds=5):
		self.nfolds = nfolds
		self.data= []
		self.class_str = {}
		self.feature_typeset = feature_typeset
		self.name = ''
		self.clf = None
		self.feats = None
		self.nonzero_idx = None
		self.labels = None
	def get_nonzero(self, combines=None):
			return na.concatenate([d.get_nonzero(combines) for d in self.data])
	def get_features(self, combines=None, nonzero_flag=True):
		assert len(self.data)>0
		featured_data = [d.index_features(combines, self.feature_typeset).T for d in self.data]
		for i,j in zip(featured_data[:-1], featured_data[1:]):
			assert (i.columns==j.columns).all()

		labels = list(chain(*[d.get_truth_labels() for d in self.data]))
		features = pd.concat(featured_data, axis=0, ignore_index=True)		
		assert features.shape[0]==len(labels)
		features = features.fillna(0)
		features = features.astype(float)
		#featuresw = feats.T.reset_index()
		#featuresw.to_csv('t.csv')
		labels = na.array([self.class_str[l] for l in labels])
		#print feats
		if nonzero_flag:
			nonzero_idx = self.get_nonzero(combines)

			#print nonzero_idx
			#print feats.shape, nonzero_idx.shape
			features,labels = features.ix[nonzero_idx,:], labels[nonzero_idx]
			self.nonzero_idx = nonzero_idx

		self.feats = features
		self.labels = labels
		return features, labels
	def preload_class_labels(self, class_labels):
		assert len(self.class_str)==0
		for l in class_labels:
			self.class_str[l] = len(self.class_str)

	def add(self, data):
		self.data.append(data)
		for l in data.get_truth_labels():
			assert l in self.class_str
#		if not self.class_str.has_key(data.class_label):
#			self.class_str[data.class_label] = len(self.class_str)
	def _remap_classes(self):
		v = sorted(set(self.class_str.values()))
		map_old_to_new = dict([(j,i) for i,j in enumerate(v)])
		for k in self.class_str.keys():
			self.class_str[k] = map_old_to_new[self.class_str[k]]
 
	def merge_hom(self):
		assert 'homA' in self.class_str.viewkeys()
		assert 'homB' in self.class_str.viewkeys()
		self.class_str['homB'] = self.class_str['homA'] 
		self._remap_classes()
	def merge_phase(self):
		assert 'pA' in self.class_str.viewkeys()
		assert 'pB' in self.class_str.viewkeys()
		self.class_str['p'] = self.class_str['pA']		
		self.class_str['pB'] = self.class_str['pB']
		self._remap_classes()

	#def _nonzero_features(self, feature):
	#	# assumes feats are positive floats
	#	nonzeros_idx = feature.apply(na.sum, axis=1)>PRECISION
	#	return nonzeros_idx
	#	#return feats.ix[nonzeros_idx,:], labels[nonzeros_idx]

	def check_features(self, features, labels):
		assert len(self.data)>0
		assert len(self.class_str)>0
		assert features.values.shape[0]==len(labels)		
		print "# of Features:", features.shape[1]
		print "# of Samples:", features.shape[0]
		print "# of Distinct Labels:", len(set(list(labels)))
		print "# Bin Counts, class", na.bincount(labels, minlength=5), self.class_str
		print "# Random States", RANDOM_STATE
	def _get_score_freq(self, features, labels, bins=None):
		if bins is None:
			bins = list(na.linspace(0,4,65))+[100,]
		p = self.clf.predict_log_proba(features)
		pred = na.argmax(p,axis=1)
		max_val = na.max(p,axis=1)
		d = na.copy(p)

		d[na.arange(d.shape[0]),pred] = -na.inf
		max_vald = na.max(d,axis=1)
		delta_val = max_val-max_vald

		true_pred = labels==pred
		false_pred = na.logical_not(true_pred)

		vt,b = na.histogram(-max_val[true_pred], bins=bins)
		vf,b = na.histogram(-max_val[false_pred], bins=bins)
		dt,b = na.histogram(delta_val[true_pred], bins=bins)
		df,b = na.histogram(delta_val[false_pred], bins=bins)
		#print p.shape
		#print self.class_str
		#print p[na.logical_and(delta_val<0.5, true_pred),:]
		#print "ZeroFeature TP", na.sum(feats.ix[true_pred,:].apply(na.sum, axis=1)==0)
		#print "ZeroFeature FP", na.sum(feats.ix[false_pred,:].apply(na.sum, axis=1)==0)
		return pd.DataFrame([vt,vf,dt,df], index=['NegLogProbAT', 'NegLogProbAF', 'DeltaLogProbAT', 'DeltaLogProbAF'], columns=bins[:-1], dtype=na.int), bins
	
	def write_scores_to_csv(self, features, labels, csv_prefix_fpath=None, class_label=''):
		#feats, labels = self._rm_zero_features(feats, labels)
		if class_label== '':
			frame, bins = self._get_score_freq(features, labels)
		else: 
			label_idx = labels==self.class_str[class_label[1:]]
			frame, bins = self._get_score_freq(features.ix[label_idx], labels[label_idx])
		if not csv_prefix_fpath is None:
			frame.to_csv('%s.%s%s.csv'%(csv_prefix_fpath, self.name,class_label))

		#print frame.T
		#print frame.T.describe()

	def _pickle_clf(self, pickle_prefix_fpath):
		if pickle_prefix_fpath is None or self.clf is None:
			return
		with open("%s.%s.pkl"%(pickle_prefix_fpath, self.name), 'wb') as f:
			pickle.dump(self.clf, f)
	def plot_pca(self, features, labels, combines):
		from matplotlib import pyplot as plt
		p = PCA(n_components=2)
		color_v = {'pA':'r', 'pB':'b', 'homA':'0.2', 'homB':'g', 'inc':'0.4'}

		show = ['pA', 'pB']
		label_idx = labels==self.class_str['pA']
		for lbl in show:
			label_idx = na.logical_or(label_idx, labels==self.class_str[lbl])
		
		pca_features = p.fit_transform(features.values[label_idx,:])		
		print "PCA Features", p.explained_variance_ratio_
		print pca_features.shape
		
		fig = plt.figure(facecolor='white', edgecolor=None)
		
		for lbl in show:
			pf = p.transform(features.values[labels==self.class_str[lbl],:])

			plt.scatter(pf[:,0], pf[:,1], c=color_v[lbl],
								linewidths=0, label=lbl)
		plt.legend()
		
		if combines is None:
			name = "%s all"%(self.name)
		else:
			name = "%s %s"%(self.name, combines)
		
		plt.title(name)
		plt.show()

	def perform(self, pickle_prefix_fpath=None, csv_prefix_fpath=None):
		raise

class SVMClassification(SVPhaseClassification):
	def __init__(self, feature_typeset=None, nfolds=5):
		SVPhaseClassification.__init__(self, feature_typeset, nfolds)
		self.name = '{0:02d}.SVM'.format(nfolds)

	def perform(self, pickle_prefix_fpath=None, csv_prefix_fpath=None, combines=None):
		features, labels = self.get_features(combines)
		#self.plot_pca(feats, labels, combines)
				
		#shuffle_idx = na.arange(len(labels))
		#na.random.shuffle(shuffle_idx)

		self.check_features(features, labels)
		clf = svm.SVC(kernel='linear', probability=True, random_state=RANDOM_STATE )
		parameters = {'C':[1,10,100]}

		cv = StratifiedKFold(labels, n_folds=self.nfolds)
		#scoring = metrics.make_scorer(metrics.recall_score, average='micro')		 
		scoring = 'accuracy'

		grid_clf = GridSearchCV(estimator=clf, param_grid=parameters, scoring=scoring, n_jobs=-1, cv=cv, refit=True, verbose=1)
		grid_clf.fit(features, labels)
		print "SVM Best Estimator", grid_clf.best_score_, grid_clf.best_params_
		scores = cross_val_score(grid_clf.best_estimator_, features, labels, scoring=scoring, cv=cv,
															n_jobs=-1, verbose=True)
		print scores
		self.clf = grid_clf.best_estimator_

		if not csv_prefix_fpath is None:
			if not combines is None:
				csv_prefix_fpath = "%s.%s"%(csv_prefix_fpath, combines)
			self.write_scores_to_csv(features, labels, csv_prefix_fpath=csv_prefix_fpath, class_label='.pA')
			self.write_scores_to_csv(features, labels, csv_prefix_fpath=csv_prefix_fpath, class_label='.pB')
			self.write_scores_to_csv(features, labels, csv_prefix_fpath=csv_prefix_fpath)
		if not combines is None:		
			pickle_prefix_fpath = "%s.%s"%(pickle_prefix_fpath, combines)
		self._pickle_clf(pickle_prefix_fpath)	

		return scores, grid_clf.best_estimator_

class KNeighborsClassification(SVPhaseClassification):
	def __init__(self, feature_typeset=None, nfolds=5):
		SVPhaseClassification.__init__(self, feature_typeset, nfolds)
		if feature_typeset is None:
			self.name = '{0:02d}.KN'.format(nfolds)
		else:
			self.name = '{0:02d}.{1}-KN'.format(nfolds,feature_typeset)

	def perform(self, pickle_prefix_fpath=None, csv_prefix_fpath=None, combines=None):
		features, labels = self.get_features(combines)
				
		self.check_features(features, labels)
		clf = KNeighborsClassifier(n_neighbors=5, weights='uniform', algorithm='brute')

		parameters = {'n_neighbors':[2,4,8,16,32]}

		cv = StratifiedKFold(labels, n_folds=self.nfolds)
		#scoring = metrics.make_scorer(metrics.recall_score, average='micro')		 
		scoring = 'accuracy'

		grid_clf = GridSearchCV(estimator=clf, param_grid=parameters, scoring=scoring, n_jobs=-1, cv=cv, refit=True, verbose=1)
		grid_clf.fit(features, labels)
		print "KN Best Estimator", grid_clf.best_score_, grid_clf.best_params_
		scores = cross_val_score(grid_clf.best_estimator_, features, labels, scoring=scoring, cv=cv,
															n_jobs=-1, verbose=True)
		print scores
		self.clf = grid_clf.best_estimator_

		if not csv_prefix_fpath is None:
			if not combines is None:
				csv_prefix_fpath = "%s.%s"%(csv_prefix_fpath, combines)
			#self.write_scores_to_csv(feats, labels, csv_prefix_fpath=csv_prefix_fpath, class_label='.pA')
			#self.write_scores_to_csv(feats, labels, csv_prefix_fpath=csv_prefix_fpath, class_label='.pB')
			#self.write_scores_to_csv(feats, labels, csv_prefix_fpath=csv_prefix_fpath)
		if not combines is None:		
			pickle_prefix_fpath = "%s.%s"%(pickle_prefix_fpath, combines)
		self._pickle_clf(pickle_prefix_fpath)	

		return scores, grid_clf.best_estimator_


class NaiveBayesClassification(SVPhaseClassification):
	def __init__(self, feature_typeset=None, nfolds=5):
		SVPhaseClassification.__init__(self, feature_typeset, nfolds)
		if feature_typeset is None:
			self.name = '{0:02d}.GNB'.format(nfolds)
		else:
			self.name = '{0:02d}.{1}-GNB'.format(nfolds,feature_typeset)

	def perform(self, pickle_prefix_fpath=None, csv_prefix_fpath=None, combines=None):
		features, labels = self.get_features(combines)
		self.check_features(features, labels)
		clf = GaussianNB()
 
		cv = StratifiedKFold(labels, n_folds=self.nfolds)
		#scoring = metrics.make_scorer(metrics.recall_score, average='micro')		 
		scoring = 'accuracy'

		scores = cross_val_score(clf, features, labels, scoring=scoring, cv=cv, n_jobs=-1, verbose=True)
		print scores
		self.clf = clf
		self.clf.fit(features, labels)

		if not csv_prefix_fpath is None:
			if not combines is None:
				csv_prefix_fpath = "%s.%s"%(csv_prefix_fpath, combines)
			self.write_scores_to_csv(features, labels, csv_prefix_fpath=csv_prefix_fpath, class_label='.pA')
			self.write_scores_to_csv(features, labels, csv_prefix_fpath=csv_prefix_fpath, class_label='.pB')
			self.write_scores_to_csv(features, labels, csv_prefix_fpath=csv_prefix_fpath)
		if not combines is None:		
			pickle_prefix_fpath = "%s.%s"%(pickle_prefix_fpath, combines)
		self._pickle_clf(pickle_prefix_fpath)	

		return scores, clf

class RandomForestClassification(SVPhaseClassification):
	def __init__(self, feature_typeset=None, nfolds=5):
		SVPhaseClassification.__init__(self, feature_typeset, nfolds)
		self.name = '{0:02d}.RandomForest'.format(nfolds)
	def perform(self, pickle_prefix_fpath=None, csv_prefix_fpath=None, combines=None):
		features, labels = self.get_features(combines)
		#self.plot_pca(feats, labels, combines)
		self.check_features(features, labels)

		#for feat in feats.columns:
		#	print feat

		clf = ensemble.RandomForestClassifier(oob_score=True, random_state=RANDOM_STATE)
		parameters = {'n_estimators':[10,20,50,100], 'max_depth':[2,5,10,20]}
		#parameters = {'n_estimators':[10,100,1000,10000]}

		cv = StratifiedKFold(labels, n_folds=self.nfolds)

		#scoring = metrics.make_scorer(metrics.recall_score, average='micro')		 
		scoring = 'accuracy'

		grid_clf = GridSearchCV(estimator=clf, param_grid=parameters, scoring=scoring, n_jobs=-1, cv=cv, refit=True, verbose=2)
		grid_clf.fit(features, labels)
		print "Rand Forest Best Estimator", grid_clf.best_score_, grid_clf.best_params_
		importances = pd.Series(grid_clf.best_estimator_.feature_importances_, index=features.columns)
		importances.sort()
		#print importances[::-1]
		
		scores = cross_val_score(grid_clf.best_estimator_, features, labels, scoring=scoring, cv=cv,
															n_jobs=-1, verbose=True)
		print scores
		self.clf = grid_clf.best_estimator_
		print "OOB_Score", self.clf.oob_score_
		if not csv_prefix_fpath is None:
			if not combines is None:
				csv_prefix_fpath = "%s.%s"%(csv_prefix_fpath, combines)
			self.write_scores_to_csv(features, labels, csv_prefix_fpath=csv_prefix_fpath, class_label='.pA')
			self.write_scores_to_csv(features, labels, csv_prefix_fpath=csv_prefix_fpath, class_label='.pB')
			self.write_scores_to_csv(features, labels, csv_prefix_fpath=csv_prefix_fpath)
		if not combines is None:		
			pickle_prefix_fpath = "%s.%s"%(pickle_prefix_fpath, combines)
		self._pickle_clf(pickle_prefix_fpath)
		return scores, grid_clf.best_estimator_

def get_classification_obj(classifier,nfolds, feature_typeset=None):
	assert classifier in set([None, 'SVM', 'RandomForest', 'GNB', 'KN'])
	if classifier is None or classifier == 'SVM':
		pclf = SVMClassification(feature_typeset=feature_typeset, nfolds=nfolds)
	elif classifier == 'RandomForest':
		pclf = RandomForestClassification(feature_typeset=feature_typeset, nfolds=nfolds)
	elif classifier == 'GNB':
		feature_typeset = 'simple_sum'
		pclf = NaiveBayesClassification(feature_typeset=feature_typeset, nfolds=nfolds)
	elif classifier == 'KN':
		#feature_typeset = 'simple_sum'
		pclf = KNeighborsClassification(feature_typeset=feature_typeset, nfolds=nfolds)
	else:
		raise
	return pclf

def realistic_classes_simple(perform_flag=True, classifier=None, evidence=None, nfolds=7):
	#assert evidence in set([None, 'only', 'hic_only', 'wgs_only'])
	#assert classifier in set([None, 'SVM', 'RandomForest'])
	# Only do pA,pB,hom and inc (four class) n_fold=6with SVM classifier.	
	pclf = get_classification_obj(classifier, nfolds)
	sv_fpath = '/home/anand/Projects/assembly/data/simulations/sim.chr20.%d.bed'
	wgs_prefix_fpath = '/media/T02/data/gm12878/sim.wgs/sim.chr20'
	#wgs_test_prefix_fpath = '/media/T01/data/gm12878/hic/sim.wgs/repA_hypoth/sim.chr20'
	hic_prefix_fpath = '/media/T02/data/gm12878/hic/dat2/sim.chr20'
	wgs2_prefix_fpath = '/media/T02/data/gm12878/sim.wgs/sim.chr19'
	hic2_prefix_fpath = '/media/T02/data/gm12878/hic/dat2/sim.chr19'
	sim_order = [4,2,0,1,3]
	
	#out_prefix_fpath = '/home/anand/Projects/assembly/data/simulations/models/sim.del-chr19-20.5norm'
	#out_prefix_fpath = '/home/anand/Projects/assembly/data/simulations/models/sim.del-chr19-20.5-2c'
	out_prefix_fpath = '/home/anand/Projects/assembly/data/simulations/models/sim.del-chr19-20.10-2c'
	#out_prefix_fpath = '/home/anand/Projects/assembly/data/simulations/models/sim.del-chr19-20.5'
	load_flag = True
	
	#adaptor19 = RPKMAdaptor(46197064,6192316)
	#adaptor20 = RPKMAdaptor(51100213,8521856)
	#adaptor19 = adaptor20 = DataAdaptor()	
	#adaptor = DataAdaptor()	
	adaptor = RPKMAdaptor(46197064+51100213,6192316+8521856)

	pclf.preload_class_labels(['inc','pA','homA', 'pB', 'homB']) 

	for i in range(5):
		mismatch_sim_num = sim_order[i]
	
		inc_data = SimulationDelData(wgs_prefix_fpath, hic_prefix_fpath, sv_fpath%i, mismatch_sim_num, 'inc', wgs_sample='.down4', hic_sample='.down2')
		inc_data.fill(adaptor=adaptor,load_flag=load_flag)
		pclf.add(inc_data)		

		for allele in ['A','B']:	
			allele_sv_fpath = '/home/anand/Projects/assembly/data/simulations/sim.chr20.%d.%s.bed'%(i, allele)
		
			p_data = SimulationDelData(wgs_prefix_fpath, hic_prefix_fpath, allele_sv_fpath, i, "p%s"%allele, wgs_sample='.down4', hic_sample='.down2')
			p_data.fill(adaptor=adaptor, load_flag=load_flag)
			pclf.add(p_data)
			hom_data = SimulationDelData(wgs_prefix_fpath, hic_prefix_fpath, allele_sv_fpath, i, "hom%s"%allele, wgs_sample='.down2', hic_sample='')			
			hom_data.fill(adaptor=adaptor, load_flag=load_flag)
			pclf.add(hom_data)

			allele2_sv_fpath = '/home/anand/Projects/assembly/data/simulations/sim.chr19.%d.%s.bed'%(i, allele)
			
			p_data = SimulationDelData(wgs2_prefix_fpath, hic2_prefix_fpath, allele2_sv_fpath, i, "p%s"%allele, wgs_sample='.down4', hic_sample='.down2')

			p_data.fill(adaptor=adaptor, load_flag=load_flag)
			pclf.add(p_data)
		
	pclf.merge_hom()

	if perform_flag:
		pclf.perform(pickle_prefix_fpath=out_prefix_fpath, csv_prefix_fpath=out_prefix_fpath, combines=evidence)	
		#pclf.perform(csv_prefix_fpath=out_prefix_fpath, combines=evidence)
		#pclf.perform(combines=evidence)
	return pclf

def _get_mapped_read_count(fpath):
	t = {}
	with open(fpath, 'rb') as f:
		for line in f:
			tokens = line.strip().split('\t')
			t[tokens[0]] = int(tokens[2])
	return t

def realistic_classes_hm(perform_flag=True, classifier=None, evidence=None, nfolds=7):
	import os
	#assert evidence in set([None, 'only', 'hic_only', 'wgs_only'])
	#assert classifier in set([None, 'SVM', 'RandomForest', 'GNB'])
	# Only do pA,pB,hom and inc (four class) n_fold=6with SVM classifier.	
	pclf = get_classification_obj(classifier, nfolds)
	
	sv_fpath='/home/anand/Projects/assembly/data/gm12878/hapmap3_r2_b36_fwd-truth.bed'
	prefix_dir = DATA_PREFIX + 'splitfiles'
	wgs_fpath = DATA_PREFIX + 'unsplit_files/wgs_bamfiles_unsplit/NA12878.hiseq.wgs.bwa.recal.bam'
	hic_stat_fpath = DATA_PREFIX + 'unsplit_files/hic_bamfiles_unsplit/chr.all.stat'
	wgs_stat_fpath = DATA_PREFIX + 'unsplit_files/wgs_bamfiles_unsplit/NA12878.hiseq.wgs.bwa.recal.stat'
	wgs_read_count = _get_mapped_read_count(wgs_stat_fpath)
	hic_read_count = _get_mapped_read_count(hic_stat_fpath)

	#out_prefix_fpath = '/home/anand/Projects/assembly/data/simulations/models/sim.del-chr19-20.5norm'
	out_prefix_fpath = '/home/anand/Projects/assembly/data/simulations/models/real.hm.5-2c'
	#out_prefix_fpath = '/home/anand/Projects/assembly/data/simulations/models/sim.del-chr19-20.5'
	load_flag = True

	chrom_name_map = dict(zip(map("chr{0:d}".format, range(1,23)), range(1,23)))
	chrom_name_map['chrX']=23 
	chrom_name_map['chrY']=24 

	avail_chr = [c for c in os.listdir(prefix_dir) if c in chrom_name_map.viewkeys()]
		
	with open(sv_fpath, 'rb') as f:
		contigs_with_loci = set([t.split('\t')[0] for t in f if t.startswith('chr')])	 

	#print wgs_read_count, hic_read_count 
	#adaptor = DataAdaptor()
	adaptor = RPKMAdaptor(wgs_read_count['chr19']+wgs_read_count['chr20'], hic_read_count['chr19']+hic_read_count['chr20'])

	pclf.preload_class_labels(['inc','pA','homA', 'pB', 'homB']) 
	for contig in sorted(avail_chr, key=chrom_name_map.get):
		if not contig in contigs_with_loci: continue
		#if not contig=='chr20': continue
		#print wgs_read_count[contig], hic_read_count[contig]
		rdata = RealDelData(wgs_fpath, prefix_dir, sv_fpath, contig, truth_flag=True)
		#rdata.fill(adaptor=RPKMAdaptor(wgs_read_count[contig], hic_read_count[contig]), load_flag=True)
		rdata.fill(adaptor=adaptor, load_flag=load_flag)
		pclf.add(rdata) 

	
	pclf.merge_hom()

	if perform_flag:
		pclf.perform(pickle_prefix_fpath=out_prefix_fpath, csv_prefix_fpath=out_prefix_fpath, combines=evidence)	
		#pclf.perform(csv_prefix_fpath=out_prefix_fpath, combines=evidence)
		#pclf.perform(combines=evidence)
	return pclf

def realistic_classes_sim_hm_mix(perform_flag=True, classifier=None, evidence=None, nfolds=7):
	#assert evidence in set([None, 'only', 'hic_only', 'wgs_only'])

	pclf = get_classification_obj(classifier, nfolds)
	
	out_prefix_fpath = '/home/anand/Projects/assembly/data/simulations/models/sim.del-chr19-20-hm.10-2c'

	sim_pclf = realistic_classes_simple(perform_flag=False, classifier=classifier, evidence=evidence, nfolds=nfolds)
	hm_pclf = realistic_classes_hm(perform_flag=False, classifier=classifier, evidence=evidence, nfolds=nfolds)
	
	pclf.preload_class_labels(['inc','pA','homA', 'pB', 'homB']) 
	
	for data in chain(sim_pclf.data, hm_pclf.data):
		pclf.add(data)

	pclf.merge_hom()

	if perform_flag:
		pclf.perform(pickle_prefix_fpath=out_prefix_fpath, csv_prefix_fpath=out_prefix_fpath, combines=evidence)	
		#pclf.perform(csv_prefix_fpath=out_prefix_fpath, combines=evidence)
		#pclf.perform(combines=evidence)
	return pclf

def realistic_classes_wgs2x(perform_flag=False, classifier=None, evidence=None, nfolds=7):
	#assert evidence in set([None, 'only', 'hic_only', 'wgs_only'])
	#assert classifier in set([None, 'SVM', 'RandomForest', 'GNB'])
	# Only do pA,pB,hom and inc (four class) n_fold=6with SVM classifier.	
	pclf = get_classification_obj(classifier, nfolds)

	wgs_prefix_fpath = '/media/T02/data/gm12878/sim.wgs/sim.chr20'
	hic_prefix_fpath = '/media/T02/data/gm12878/hic/dat2/sim.chr20'
	wgs2_prefix_fpath = '/media/T02/data/gm12878/sim.wgs/sim.chr19'
	hic2_prefix_fpath = '/media/T02/data/gm12878/hic/dat2/sim.chr19'
	sim_order = [4,2,0,1,3]
	
	out_prefix_fpath = '/home/anand/Projects/assembly/data/simulations/models/sim.del-chr19-20.wgs2xnorm'
	load_flag = True
	
	#adaptor19 = RPKMAdaptor(2*46197064,6192316)
	#adaptor20 = RPKMAdaptor(2*51100213,8521856)
	#adaptor19 = adaptor20 = DataAdaptor()	
	adaptor = RPKMAdaptor(2*(46197064+51100213),6192316+8521856)

	pclf.preload_class_labels(['inc','pA','homA', 'pB', 'homB']) 

	for i in range(5):

		for allele in ['A','B']:	
			allele_sv_fpath = '/home/anand/Projects/assembly/data/simulations/sim.chr20.%d.%s.bed'%(i, allele)
		
			p_data = SimulationDelData(wgs_prefix_fpath, hic_prefix_fpath, allele_sv_fpath, i, "p%s"%allele, wgs_sample='.down2', hic_sample='.down2')
			p_data.fill(adaptor=adaptor, load_flag=load_flag)
			pclf.add(p_data)

			allele2_sv_fpath = '/home/anand/Projects/assembly/data/simulations/sim.chr19.%d.%s.bed'%(i, allele)
			
			p_data = SimulationDelData(wgs2_prefix_fpath, hic2_prefix_fpath, allele2_sv_fpath, i, "p%s"%allele, wgs_sample='.down2', hic_sample='.down2')

			p_data.fill(adaptor=adaptor, load_flag=load_flag)
			pclf.add(p_data)
		
	pclf.merge_hom()

	if perform_flag:
		pclf.perform(pickle_prefix_fpath=out_prefix_fpath, csv_prefix_fpath=out_prefix_fpath, combines=evidence)	
		#pclf.perform(csv_prefix_fpath=out_prefix_fpath, combines=evidence)
		#pclf.perform(combines=evidence)
	return pclf

def realistic_classes_hic2x(perform_flag=False, classifier=None, evidence=None, nfolds=7):
	#assert evidence in set([None, 'only', 'hic_only', 'wgs_only'])
	#assert classifier in set([None, 'SVM', 'RandomForest', 'GNB'])
	# Only do pA,pB,hom and inc (four class) n_fold=6with SVM classifier.	
	pclf = get_classification_obj(classifier, nfolds)

	sv_fpath = '/home/anand/Projects/assembly/data/simulations/sim.chr20.%d.bed'
	wgs_prefix_fpath = '/media/T02/data/gm12878/sim.wgs/sim.chr20'
	hic_prefix_fpath = '/media/T02/data/gm12878/hic/dat2/sim.chr20'
	wgs2_prefix_fpath = '/media/T02/data/gm12878/sim.wgs/sim.chr19'
	hic2_prefix_fpath = '/media/T02/data/gm12878/hic/dat2/sim.chr19'
	sim_order = [4,2,0,1,3]
	
	out_prefix_fpath = '/home/anand/Projects/assembly/data/simulations/models/sim.del-chr19-20.hic2xnorm'
	load_flag = True
	
	#adaptor19 = RPKMAdaptor(46197064,2*6192316)
	#adaptor20 = RPKMAdaptor(51100213,2*8521856)
	#adaptor19 = adaptor20 = DataAdaptor()	
	adaptor = RPKMAdaptor(46197064+51100213,2*(6192316+8521856))

	pclf.preload_class_labels(['inc','pA','homA', 'pB', 'homB']) 

	for i in range(5):
		mismatch_sim_num = sim_order[i]

		for allele in ['A','B']:	
			allele_sv_fpath = '/home/anand/Projects/assembly/data/simulations/sim.chr20.%d.%s.bed'%(i, allele)
		 
			p_data = SimulationDelData(wgs_prefix_fpath, hic_prefix_fpath, allele_sv_fpath, i, "p%s"%allele, wgs_sample='.down4', hic_sample='')
			p_data.fill(adaptor=adaptor, load_flag=load_flag)
			pclf.add(p_data)

			allele2_sv_fpath = '/home/anand/Projects/assembly/data/simulations/sim.chr19.%d.%s.bed'%(i, allele)
			p_data = SimulationDelData(wgs2_prefix_fpath, hic2_prefix_fpath, allele2_sv_fpath, i, "p%s"%allele, wgs_sample='.down4', hic_sample='')

			p_data.fill(adaptor=adaptor, load_flag=load_flag)
			pclf.add(p_data)

	pclf.merge_hom()

	if perform_flag:
		pclf.perform(pickle_prefix_fpath=out_prefix_fpath, csv_prefix_fpath=out_prefix_fpath, combines=evidence)	
		#pclf.perform(csv_prefix_fpath=out_prefix_fpath, combines=evidence)
		#pclf.perform(combines=evidence)
	return pclf

def phased_classes_break_error(start_size_range, end_size_range, right_break_error, left_break_error, classifier=None, evidence=None, nfolds=10):
	assert evidence in set([None, 'only', 'hic_only', 'wgs_only'])
	assert classifier in set([None, 'SVM', 'RandomForest', 'GNB'])
	# Only do pA,pB,hom and inc (four class) n_fold=6with SVM classifier.	
	pclf = get_classification_obj(classifier, nfolds)
	
	wgs_prefix_fpath = '/media/T02/data/gm12878/sim.wgs/sim.chr20'
	hic_prefix_fpath = '/media/T02/data/gm12878/hic/dat2/sim.chr20'
	
	load_flag = True
	
	#adaptor20 = RPKMAdaptor(51100213,8521856)
	#adaptor19 = adaptor20 = DataAdaptor()	
	adaptor20 = RPKMAdaptor(46197064+51100213,6192316+8521856)

	pclf.preload_class_labels(['inc','pA','homA', 'pB', 'homB']) 

	for i in range(5):
		for allele in ['A','B']:	
			allele_sv_fpath = '/home/anand/Projects/assembly/data/simulations/break_error/sim.chr20.{0:d}.{1:s}.{2:03d}.{3:03d}.{4:04d}.{5:04d}.bed'.format(i, allele, start_size_range, end_size_range, right_break_error, left_break_error)
		
			p_data = SimulationDelData(wgs_prefix_fpath, hic_prefix_fpath, allele_sv_fpath, i, "p%s"%allele, wgs_sample='.down4', hic_sample='.down2')
			p_data.fill(adaptor=adaptor20, load_flag=load_flag)
			pclf.add(p_data)

	pclf.merge_hom()

	return pclf

if __name__=='__main__':
	import sys
	RANDOM_STATE=int(time.time())%1000000
	#realistic_classes_simple(classifier='RandomForest', evidence=None,perform_flag=True, nfolds=10)
	realistic_classes_sim_hm_mix(classifier='KN',evidence=None,nfolds=10)
	sys.exit(1)	
	
	#for c in ['SVM','RandomForest']:
		#for c in ['GNB',]:
	for c in ['KN','SVM','RandomForest']:
		for e in [None, 'hic_only', 'wgs_only', 'only']:
			#for e in ['wgs_only',]:
			#for nfold in [5,7,10]:
			for nfold in [10,]:
				RANDOM_STATE += 1
				#realistic_classes_simple(classifier=c,evidence=e,nfolds=nfold)
				#realistic_classes_hm(classifier=c,evidence=e,nfolds=nfold)
				realistic_classes_sim_hm_mix(classifier=c,evidence=e,nfolds=nfold)
 
