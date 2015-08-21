import pandas as pd
from matplotlib import pyplot as plt, gridspec
from matplotlib.patches import Wedge
from itertools import izip
import numpy as np
import os
from scipy.stats import norm, stats
from svphase.learn.evaluation import ClassLabel

LABELS = ClassLabel()

#pd.set_printoptions(max_rows=10000,max_columns=10000)

def trap_area(x1,x2,y1,y2):
	assert x2>=x1
	assert y2>=y1
	
	w = x2-x1
	a = w*y1
	return a + 0.5 * w * (y2-y1)

def get_mann_whitney_stat(df):
	df['acc'] = df['pred']==df['truth']
	df['nacc'] = df['pred']!=df['truth']
	df = df.sort(['scores', 'nacc'], ascending=False)
	#print df[['scores','acc']]
	#print df['acc']
	rank_value = np.cumsum(np.logical_not(df['acc'].values))
	#print rank_value
	rank_sum = np.sum(rank_value[df['acc']])
	print (rank_value[df['acc']]), rank_sum
	n1 = np.sum(df['acc'])
	n2 = np.sum(np.logical_not(df['acc'].values))
	return pd.Series([rank_sum, n1, n2], index=['u','n1','n2'])

def standard_normal_mw(s):
	n1,n2 = s['n1'], s['n2']
	mean = n1*n2 * 0.5
	var = (n1*n2*(n1+n2+1)/12.0)**0.5
	return (s['u']-mean)/var

def compute_mw_cutoff(x_s,y_s, sig_flag=True):
	#print x_s
	#print y_s
	x_z = standard_normal_mw(x_s)
	y_z = standard_normal_mw(y_s)
	d = x_z - y_z
	if sig_flag:
		return norm.cdf(d,0,2)
	else:
		return d


class Scoring(object):
	def __init__(self, truth_fpath, bed_fpath, reset_min=None):
		self.truth_fpath = truth_fpath
		self.bed_fpath = bed_fpath
		self.truth = None
		self.all_scores = None
		self.reset_min = reset_min
		self.df = None
	def get_mw(self, df=None):
		if df is None:
			df = self.find_scores()
		return get_mann_whitney_stat(df)

	def find_scores(self, scores='pred', exclude_no_data_flag=True):
		assert scores=='pred' or scores=='truth'
		truth_df = pd.read_csv(self.truth_fpath, sep='\t', index_col=None, header=None, names=['contig','start','end', 'truth'])
		pred_df = pd.read_csv(self.bed_fpath, sep='\t', index_col=None, header=0).rename(columns={'best':'pred'})
		df = pd.concat([truth_df.set_index(['contig','start','end']), pred_df.set_index(['contig','start','end'])], axis=1)	
		df = df.loc[df['truth'].apply(LABELS.is_class),:]
		
		df = df.dropna(subset=['pred',])
		if exclude_no_data_flag:
			df = df.loc[df['pred']!='no_data',:]

		self.all_scores = df[LABELS.classes]
		self.truth = pd.DataFrame(df[['truth','pred']])
		df['scores'] = [self.all_scores.ix[i,s] for i,s in self.truth[scores].iteritems()]
		if not self.reset_min is None:
			df.ix[df['scores']<self.reset_min,'scores'] = self.reset_min

		self.df = df
		#print df
		return df

class ROC(Scoring):
	def __init__(self, bed_fpath):
		Scoring.__init__(self, bed_fpath)
		self.auc = {}
	def _set_tp_fp_rates(self,method='pred'):
		ndf = self.truth.ix[np.logical_or(self.truth['truth']=='pA',self.truth['truth']=='pB'),:]
		ndf['binary_p'] = 0.0
		ndf['binary_n'] = 0.0
		ndf.ix[ndf['truth']==ndf['pred'],'binary_p'] = 1.0
		ndf.ix[ndf['truth']!=ndf['pred'],'binary_n'] = 1.0
		
		if method == 'pred':
			ndf['scores'] = [self.all_scores.ix[i,s] for i,s in ndf['pred'].iteritems()]
		elif method == 'diff':
			t = [sorted(s) for i,s in self.all_scores.iterrows()]
			ndf['scores'] = map(lambda x:x[-1]-x[-2], t)
		else:
			raise "No method selected"
		
		ndf = ndf.sort(column=['scores'],ascending=False)
		ndf['tpr'] = np.cumsum(ndf['binary_p'])/ndf['binary_p'].sum()
		ndf['fpr'] = np.cumsum(ndf['binary_n'])/ndf['binary_n'].sum()
		#ndf['fpr'] = 1-(np.cumsum(ndf['binary_fp'][::-1])[::-1])/n
		return ndf

	def _auc(self, df):
		t = df[['fpr','tpr']]
		auc = 0
		for (i,s1),(j,s2) in izip(t.ix[:-1,:].iterrows(), t.ix[1:,:].iterrows()):
			auc += trap_area(s1['fpr'],s2['fpr'],s1['tpr'],s2['tpr'])
		return auc

	def find_scores(self):
		Scoring.find_scores(self)
		pred_roc = self._set_tp_fp_rates(method='pred')
		diff_roc = self._set_tp_fp_rates(method='diff')

		print self._auc(pred_roc)
		print self._auc(diff_roc)
		
		return {'pred':pred_roc, 'diff':diff_roc}
	
	def generate_plot(self, ax=None):
		rocs = self.find_scores()
		print "="*20, 'pred', "="*20
		print rocs['pred']
		print "="*20, 'diff', "="*20
		print rocs['diff']
	
		fig = None
		if ax is None:
			fig = plt.figure()
			ax = fig.add_subplot(111)
		

		ax.plot(rocs['pred']['fpr'],rocs['pred']['tpr'],c='r',lw=5,label='max prediction score\nAUC={0:0.2f}'.format(self._auc(rocs['pred'])))
		ax.plot(rocs['diff']['fpr'],rocs['diff']['tpr'],c='b',lw=5,label='difference in top two scores\nAUC={0:0.2f}'.format(self._auc(rocs['diff'])))
		
		ax.set_xlabel('False positive rate')
		ax.set_ylabel('True positive rate')
		ax.set_title("{0} ROC (P={1:d},N={2:d})".format(os.path.basename(self.bed_fpath).split('-')[1], int(sum(rocs['pred']['binary_p'])), int(sum(rocs['pred']['binary_n']))))
		plt.legend(loc='lower right')
		#ax.set_xlim(-0.01,0.2)
		#ax.set_ylim(0.80,1.01)
	
		return fig

class ClassReferenceROC(ROC):
	def __init__(self, bed_fpath):
		ROC.__init__(self, bed_fpath)
	def _set_tp_fp_rates(self, c):
		ndf = self.truth.ix[np.logical_or(self.truth['truth']=='pA',self.truth['truth']=='pB'),:]
		ndf['binary_p'] = 0.0
		ndf['binary_n'] = 0.0
		ndf.ix[ndf['truth']==c,'binary_p'] = 1.0
		ndf.ix[ndf['truth']!=c,'binary_n'] = 1.0
		ndf['scores'] = self.all_scores[c]
		ndf = ndf.sort(column=['scores'],ascending=False)
		ndf['tpr'] = np.cumsum(ndf['binary_p'])/ndf['binary_p'].sum()
		ndf['fpr'] = np.cumsum(ndf['binary_n'])/ndf['binary_n'].sum()
		return ndf

	def find_scores(self):
		Scoring.find_scores(self)
		pA = self._set_tp_fp_rates('pA')
		pB = self._set_tp_fp_rates('pB')

		print self._auc(pA)
		print self._auc(pB)
		
		return pA, pB
	
	def generate_plot(self, ax=None):
		pA,pB = self.find_scores()
		print "="*20, 'pA', "="*20
		print pA
		print "="*20, 'pB', "="*20
		print pB  
	
		fig = None
		if ax is None:
			fig = plt.figure()
			ax = fig.add_subplot(111)
		ax.plot(pA['fpr'],pA['tpr'],c='b',lw=5,label='pA AUC={0:0.2f}\n(P={1:d} N={2:d})'.format(self._auc(pA), int(sum(pA['binary_p'])), int(sum(pA['binary_n']))))
		ax.plot(pB['fpr'],pB['tpr'],c='r',lw=5,label='pB AUC={0:0.2f}\n(P={1:d} N={2:d})'.format(self._auc(pB), int(sum(pB['binary_p'])), int(sum(pB['binary_n']))))
		
		ax.set_xlabel('False positive rate')
		ax.set_ylabel('True positive rate')
		ax.set_title("{0} ROC".format(os.path.basename(self.bed_fpath).split('-')[1]))
		plt.legend(loc='lower right')
		#ax.set_xlim(-0.01,0.2)
		#ax.set_ylim(0.80,1.01)
	
		return fig

def dual_half_circle(ax, center, radius, angle=0, colors=('w','k'), **kwargs):
	theta1, theta2 = angle, angle+180
	w1 = Wedge(center, radius, theta1, theta2, fc=colors[0], **kwargs)
	w2 = Wedge(center, radius, theta2, theta1, fc=colors[1], **kwargs)
	for wedge in [w1, w2]:
		ax.add_artist(wedge)  
	return [w1, w2]

class ScoringCompare(object):
	def __init__(self, bed1_fpath, bed2_fpath, bed1_name, bed2_name, reset_min=None):
		self.xlabel = bed1_name
		self.ylabel = bed2_name 
		self.reset_min = reset_min
		self.x = Scoring(bed1_fpath, reset_min=self.reset_min)
		self.y = Scoring(bed2_fpath, reset_min=self.reset_min)
	def generate_plot(self, ax=None, min_log_prob=-1.2, radius=0.05):
		if not self.reset_min is None:
			min_log_prob = self.reset_min

		fig = None
		if ax is None:
			fig = plt.figure()
			ax = fig.add_subplot(111)
		x_df = self.x.find_scores()
		y_df = self.y.find_scores()
		assert (x_df.index==y_df.index).all()

		colors = {True:'b', False:'r'}
		x_s = self.x.get_mw(x_df)
		y_s = self.y.get_mw(y_df)
		#ax.scatter(x_df['scores'],y_df['scores'],lw=0)
		print "X {0} Y {1}".format(self.xlabel, self.ylabel)
		print 'X Truth', np.sum(x_df['pred']==x_df['truth'])
		print "X", x_df['scores'].describe()
		print "X U and n", x_s['u']
		print 'Y Truth', np.sum(y_df['pred']==y_df['truth'])
		print "Y", y_df['scores'].describe()
		print "Y U and n", y_s['u']
		compute_mw_cutoff(x_s,y_s)

		print "pearson r", stats.pearsonr(x_df['scores'],y_df['scores'])

		for c1,c2,s1,s2 in izip(x_df['pred']==x_df['truth'], y_df['pred']==y_df['truth'], x_df['scores'],y_df['scores']):
			dual_half_circle(ax, (s1,s2), radius=radius, angle=225, colors=(colors[c1], colors[c2]), lw=0, alpha=0.6)

		ax.set_xlabel(self.xlabel)
		ax.set_ylabel(self.ylabel)
		#min_log_prob = -1.2
		#min_log_prob = -9
		ax.set_xlim(min_log_prob,0)
		ax.set_ylim(min_log_prob,0)
		return fig

class ScoringCompareWithTruth(object):
	def __init__(self, truth_fpath, csv_fpath, xname, reset_min=None):
		self.xlabel = xname
		self.reset_min = reset_min
		self.x = Scoring(truth_fpath,csv_fpath, reset_min=self.reset_min)
	def generate_plot(self, clabel, ax=None, num_of_bins=25, min_log_prob=-1.2, xaxis_off=False):
		if not self.reset_min is None:
			min_log_prob = self.reset_min

		fig = None
		if ax is None:
			fig = plt.figure()
			ax = fig.add_subplot(111)
		x_df = self.x.find_scores('pred')
		#y_df = self.x.find_scores('truth')
		#assert (x_df.index==y_df.index).all()

		colors = {True:'b', False:'r'}
		#ax.scatter(x_df['scores'],y_df['scores'],lw=0)
		#print "X {0} Y {1}".format(self.xlabel, self.ylabel)
		print 'X Truth', np.sum(x_df['pred']==x_df['truth'])
		print "X", x_df['scores'].describe()

		#print "pearson r", stats.pearsonr(x_df['scores'],y_df['scores'])

		#for c1,c2,s1,s2 in izip(x_df['pred']==x_df['truth'], y_df['pred']==y_df['truth'], x_df['scores'],y_df['scores']):
		#	dual_half_circle(ax, (s1,s2), radius=radius, angle=225, colors=(colors[c1], colors[c2]), lw=0, alpha=0.6)

		x_df = x_df.loc[x_df['truth']==clabel,:]
		x_df_p = x_df.loc[x_df['pred']==x_df['truth'],:]
		x_df_n = x_df.loc[x_df['pred']!=x_df['truth'],:]

		bins = np.linspace(min_log_prob, 0, num_of_bins)
		ax.hist(x_df_p['scores'].values, bins=bins, color=colors[True], lw=0, alpha=0.6)
		ax.hist(x_df_n['scores'].values, bins=bins, color=colors[False], lw=0, alpha=0.6)

		ax.set_ylabel("{0} deletion\ncount".format(clabel))
		#min_log_prob = -1.2
		#min_log_prob = -9
		if xaxis_off:
			plt.setp(ax.get_xticklabels(), visible=False)
		else:
			ax.set_xlabel(self.xlabel)
			ax.set_xlim(min_log_prob,0)
		ax.set_ylim(0,13)
		return fig

class ScoringMWMatrix(object):
	def __init__(self, bed_fpaths=[], bed_names=[], reset_min=None):
		self.bed_fpaths = bed_fpaths
		self.bed_names = bed_names 
		self.reset_min = reset_min
		self.s = [Scoring(b, reset_min=self.reset_min) for b in self.bed_fpaths]
	def generate_matrix(self):
		u = pd.DataFrame([x.get_mw() for x in self.s], index=self.bed_names, columns=['u','n1','n2']).T
		print "MW Statistic"
		print u

		n = len(self.bed_names)
		significance = pd.DataFrame(np.zeros((n,n)), index=u.columns, columns=u.columns)
		d = pd.DataFrame(np.zeros((n,n)), index=u.columns, columns=u.columns)
		for i in xrange(n):
			i_name = u.columns[i]
			for j in xrange(i+1,n):
				j_name = u.columns[j]
				print i_name, j_name, compute_mw_cutoff(u[i_name],u[j_name])

				significance.ix[i_name,j_name] = compute_mw_cutoff(u[i_name],u[j_name])
				significance.ix[j_name,i_name] = compute_mw_cutoff(u[j_name],u[i_name])
				d.ix[i_name,j_name] = compute_mw_cutoff(u[i_name],u[j_name], sig_flag=False)

		print significance
		print d


def create_supplemental_roc_figure(rf_fpath, svm_fpath, fig_fpath=None):
	fig = plt.figure(figsize=(11,6), facecolor=None, edgecolor=None)
	gs = gridspec.GridSpec(1,2)

	if class_ref_flag:
		svm_obj = ClassReferenceROC(svm_fpath)
		rf_obj = ClassReferenceROC(rf_fpath) 
	else:
		svm_obj = ROC(svm_fpath)
		rf_obj = ROC(rf_fpath) 

	svm_ax = fig.add_subplot(gs[0,0])
	#print svm_obj.find_scores()
	svm_obj.generate_plot(ax=svm_ax)
	
	rf_ax = fig.add_subplot(gs[0,1])
	#print rf_obj.find_scores()
	rf_obj.generate_plot(ax=rf_ax)
	
	gs.tight_layout(fig)
	if fig_fpath is None:
		plt.show()
	else:
		plt.savefig(fig_fpath)

def create_scat_figure(rf_fpath, svm_fpath, rf_label, svm_label, title, fig_fpath=None):
	fig = plt.figure(figsize=(11,6), facecolor=None, edgecolor=None)
	gs = gridspec.GridSpec(1,1)
	
	c = ScoringCompare(rf_fpath, svm_fpath, rf_label, svm_label)
	
	ax = fig.add_subplot(gs[0,0], aspect='equal')
	c.generate_plot(ax=ax)
	t1 = ax.set_title(title)
	t1.set_style('italic')

	gs.tight_layout(fig)
	if fig_fpath is None:
		plt.show()
	else:
		plt.savefig(fig_fpath) 

def create_scat_class_figure(truth_fpath, csv_fpath, label, title, fig_fpath=None):
	fig = plt.figure(figsize=(11,6), facecolor=None, edgecolor=None)
	gs = gridspec.GridSpec(3,1)
		

	c = ScoringCompareWithTruth(truth_fpath, csv_fpath, label)
	
	ax1 = fig.add_subplot(gs[0,0])
	ax2 = fig.add_subplot(gs[1,0], sharex=ax1)
	ax3 = fig.add_subplot(gs[2,0], sharex=ax1)
	c.generate_plot('pA',ax1, xaxis_off=True)
	c.generate_plot('pB',ax2, xaxis_off=True)
	c.generate_plot('hom',ax3)
	t1 = ax1.set_title(title)
	#t1.set_style('italic')

	plt.legend(loc='upper left')
	gs.tight_layout(fig)
	if fig_fpath is None:
		plt.show()
	else:
		plt.savefig(fig_fpath) 
 

"""
def create_mw_mat():
	rf_fpath, svm_fpath = get_input_fpaths(confident_flag=False)
	sskn_fpath, kn_fpath = get_knn_input_fpaths()
	gs_df = get_gs_df()

	w = ScoringMWMatrix(bed_fpaths=[sskn_fpath, kn_fpath, rf_fpath, svm_fpath], bed_names=['ssKNN', 'KN', 'RandomForest', 'SVM'])
	w.generate_matrix()

	print "Mills Gold Standard"

	w = ScoringMWMatrix(bed_fpaths=[sskn_fpath, kn_fpath, rf_fpath, svm_fpath], bed_names=['ssKNN', 'KN', 'RandomForest', 'SVM'], gs_df=gs_df)
	w.generate_matrix()

def get_null_mw_ex_df(n=100):
	df = pd.DataFrame(['A',]*n, columns=['pred'])
	df['truth'] = df['pred']
	df['scores'] = np.arange(n)
	return df
"""

def mw_ex():
	print "Examples"
	print "50, 50, 100, missclassification rank 1"
	ex1 = get_null_mw_ex_df()
	diff = 10
	ex1.ix[diff:49+diff,'truth'] = 'B'
	ex1_u = get_mann_whitney_stat(ex1)
	print "Ex1", ex1_u
	ex2 = get_null_mw_ex_df()
	ex2.ix[:49,'truth'] = 'B'
	ex2_u = get_mann_whitney_stat(ex2)
	print "Ex2", ex2_u
	print compute_mw_cutoff(ex1_u,ex2_u)

def mw_ex_simple():
	print "Examples"
	print "1, 99, 100, missclassification rank 1"
	ex1 = get_null_mw_ex_df()
	ex1.ix[99,'truth'] = 'B'
	ex1_u = get_mann_whitney_stat(ex1)
	print "Ex1", ex1_u
	ex2 = get_null_mw_ex_df()
	ex2.ix[0,'truth'] = 'B'
	ex2_u = get_mann_whitney_stat(ex2)
	print "Ex2", ex2_u
	print compute_mw_cutoff(ex1_u,ex2_u)


if __name__ == '__main__':
	import argparse

	parser = argparse.ArgumentParser(description="Plot scatter plot with scores")

	parser.add_argument('truth_fpath', help='Truth fpath')
	parser.add_argument('csv_fpath', help='Prediction file with log probabilities for each class')
	parser.add_argument('xlabel', help='Label on X axis (first csv file)')
	parser.add_argument('title', help='Title of figure')
	parser.add_argument('-o', dest='save', default=None, help='Save figure')

	args = parser.parse_args()
	
	create_scat_class_figure(args.truth_fpath, args.csv_fpath, args.xlabel, args.title, fig_fpath=args.save)
	
	
	
