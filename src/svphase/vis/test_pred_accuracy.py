import os
import numpy as np
import pandas as pd
import matplotlib as mpl
from matplotlib import pyplot as plt
from svphase.utils.config import FONT, COLORS
from svphase.learn.evaluation import ClassLabel

#def _get_model_and_version(model_stat):
#	fname = os.path.basename(model_stat).split('.')
#	return fname[0], '.'.join(fname[1:3])

def plot(truth_bed, pred_csv, out_fpath=None, low_ylim=0.40, high_ylim=1.00):
	label_obj = ClassLabel()
	truth_df = pd.read_csv(truth_bed, sep='\t', index_col=None, header=None, names=['contig', 'start', 'end','truth'])
	pred_df = pd.read_csv(pred_csv, sep='\t', index_col=None, header=0)

	df = pd.concat([truth_df.set_index(['contig','start','end']), pred_df.set_index(['contig','start','end'])], axis=1)	

	# Only include truth labeled items
	df = df.loc[df['truth'].apply(label_obj.is_class),:]

	# Only include items in predictions
	df = df.loc[np.logical_and(pd.notnull(df['truth']),pd.notnull(df['best'])),:]
	
	# count excluded deletions
	excluded_count = (df['best']=='no_data').sum()
	print "Excluded", excluded_count
	# Only include truth labeled items
	df = df.loc[df['best']!='no_data',:]

	df['correct'] = df.apply(lambda s:s['best']==s['truth'], axis=1)
	#for idx, s in df.iterrows():
	#	print s[['best','truth','correct']]
	
	truth_classes_g = df.groupby('truth')
	print "Number of errors", sum(df['correct']==False)
	acc_df = pd.concat([truth_classes_g.apply(lambda d:d['correct'].sum()),df.groupby('truth').apply(lambda d:len(d.index))], axis=1).rename(columns={0:'correct',1:'total'})
	
	acc_df = acc_df.T
	acc_df['excl.'] = 0
	acc_df.loc['total','excl.'] = excluded_count 

	acc_df = acc_df.T
	print acc_df, acc_df.columns, acc_df.index
	acc_df['acc'] = acc_df['correct'].astype(float)/acc_df['total']
	print acc_df, acc_df.columns, acc_df.index

	classes_count = [c for c in acc_df.index if c in label_obj.classes]
	print 'Overall Accuracy: %0.4f'%(acc_df.loc[classes_count,'correct'].sum()/float(acc_df.loc[classes_count,'total'].sum()))
	print 'Number of Deletions %d, Number of deletions with data %d'%(acc_df.loc[:,'total'].sum(),acc_df.loc[classes_count,'total'].sum())

	fig = plt.figure(figsize=(8,3), facecolor='white', edgecolor='none')
	mpl.rc('font', **FONT)
	ax = fig.add_subplot(111)

	order = [c for c in label_obj.classes]
	order = order[:-1] + ['excl.'] + order[-1:]
	# filter order for classes with data (except excl.)
	order = [c for c in order if c in acc_df.index]
		
	rename_labels = dict(zip(order, order))
	rename_labels['hom']='hom.'
	rename_labels['inc']='incorrect'
	
	class_labels = map(rename_labels.get, order)
	white_space = 0.01

	width = acc_df.loc[order,'total']
	cumsum = width.cumsum()
	width /= float(cumsum[-1])
	print width

	left = np.zeros(len(width), dtype=float)
	left[1:] = width.cumsum()[:-1]
	left += white_space
	width = width-2*white_space

	ax.bar(left, np.ones(len(left),dtype=float), width=width, color=COLORS['bg'], linewidth=0) 
	ax.bar(left, acc_df.loc[order,'acc'].values, width=width, color=COLORS['cor'], linewidth=0) 

	for px, w, c,c_w,c_c in zip(left, width, class_labels, acc_df.loc[order,'total'], acc_df.loc[order,'acc']):
		ax.text(px+ w/2, low_ylim, c, ha='center', va='bottom', size='large')
		ax.text(px+ w/2, high_ylim, "%d"%c_w, ha='center', va='bottom', size='small')
		ax.text(px+ w/2, c_c-0.01, "%0.2f"%c_c, ha='center', va='top', size='small')
	
	ax.set_ylabel('Accuracy')
	ax.set_xticks([])
	ax.set_ylim(low_ylim, high_ylim)

	ax.set_frame_on(False)
	ax.get_yaxis().tick_left()
	ax.yaxis.set_major_formatter(plt.FormatStrFormatter('%0.2f'))
	ax.get_xaxis().tick_bottom()
	#ax.axes.get_xaxis().set_visible(False)
	fig.tight_layout()
	if out_fpath is None:
		plt.show()
	else:
		plt.savefig(out_fpath)

if __name__=='__main__':
	import argparse
	
	parser = argparse.ArgumentParser(description='Plot accuracies for each class from a truth file and a set of predictions')
	
	parser.add_argument('-o', dest='save', default=None, help="Path to save figure to")
	parser.add_argument('truth', help="truth.bed file, with loci+truth labels")
	parser.add_argument('pred', help="prediction csv file from InPhadel, (likely used truth.bed to run InPhadel)")

	args = parser.parse_args()
	plot(args.truth, args.pred, out_fpath=args.save)
