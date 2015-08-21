import os
import numpy as np
import pandas as pd
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import gridspec

from svphase.learn.features import PreloadFeatures
from svphase.analyze.compare_model_stats import get_model_and_version, main_per_class as compare_model_stats
from svphase.learn.evaluation import ClassLabel
from svphase.utils.config import FONT, COLORS

def plot(models_dir, model, model_versions, rlabels, title, low_ylim=0.4, high_ylim=1.0):
	assert len(rlabels)==len(model_versions)

	print models_dir, model, model_versions
	n = len(rlabels)
	model_stats = map("{0}/{1}.{2}.pkl.stat".format, [models_dir,]*n, [model,]*n, model_versions)

	excluded_counts = {}
	for v in model_versions:
		feat = PreloadFeatures("{0}/feats.{1}".format(models_dir, v))
		excluded_counts[v] = 1
		excluded_counts[v] = sum(feat.get_nonzero()==False)

	print "Comparing Models", model_stats
	statsdf = compare_model_stats(model_stats)
	original_order = pd.DataFrame(np.arange(len(model_stats)), index=pd.MultiIndex.from_tuples(map(get_model_and_version, model_stats)), columns=['original_idx'])
	#print original_order
	statsdf = pd.concat([statsdf, original_order], axis=1).reset_index().set_index('original_idx').sort_index()
	statsdf = statsdf.set_index('version')
	
	excl_c_k, excl_c_v = zip(*excluded_counts.items())
	print excl_c_k, excl_c_v 
	excl_c = pd.DataFrame([excl_c_v,], index=['tpred:excl.',], columns=excl_c_k).T
	statsdf = pd.concat([statsdf, excl_c], axis=1)

	print statsdf
	fig = plt.figure(figsize=(8,6.8),facecolor='white', edgecolor='none')
	mpl.rc('font', **FONT)

	outer_grid = gridspec.GridSpec(1,1,wspace=0.0, hspace=0.0)
	outer_grid.update(left=0.08, right=0.92, hspace=0.25)

	tax = fig.add_subplot(outer_grid[0])
	tax.set_title(title, size='x-large', y=1.05)
	tax.set_ylabel('Avg. test accuracy', size='x-large', labelpad=40)
	tax.set_frame_on(False)
	tax.set_xticks([])
	tax.set_yticks([])	
	
	gs = gridspec.GridSpecFromSubplotSpec(len(rlabels),1, subplot_spec=outer_grid[0],hspace=0.35)

	label_obj = ClassLabel()	

	order = ['tpred:'+c for c in label_obj.classes]
	order = order[:-1] + ['tpred:excl.',] + order[-1:]
	class_labels = [c[6:] for c in order]
	
	rename_labels = dict(zip(class_labels, class_labels))
	rename_labels['hom']='hom.'
	rename_labels['inc']='incorrect'
	
	class_labels = map(rename_labels.get, class_labels)

	white_space = 0.01
	for i, (v,stats), rlabel in zip(xrange(len(rlabels)), statsdf.iterrows(), rlabels):
		ax = fig.add_subplot(gs[i,0])
		print rlabel
		width = stats[order]
		cumsum = width.cumsum()
		width /= float(cumsum[-1])

		left = np.zeros(len(width), dtype=float)
		left[1:] = width.cumsum()[:-1]
		left += white_space
		width = width-2*white_space

		avg_c = map(lambda c:'avg'+c[5:], order)	
		std_c = map(lambda c:'std'+c[5:], order)	
		print v, stats, left, stats[avg_c].values, stats[std_c].values

		ax.bar(left, np.ones(len(left),dtype=float), width=width, color=COLORS['bg'], linewidth=0) 
		ax.bar(left, stats[avg_c].values, width=width, yerr=stats[std_c].values, color=COLORS['cor'], linewidth=0) 

		#ylim_for_count = min(stats[avg_c][stats[avg_c].values>0])-0.05

		for px, w, c,c_w in zip(left, width, class_labels, stats[order]):
			ax.text(px+ w/2, low_ylim, c, ha='center', va='bottom', size='large')
			ax.text(px+ w/2, high_ylim, "%d"%c_w, ha='center', va='bottom', size='small') 	

		ax.set_xticks([])
		#ax.set_xticklabels(tuple(labels), rotation=80, multialignment='center')

		ax.set_ylim(low_ylim, high_ylim)
		ax.set_frame_on(False)
		ax.get_yaxis().tick_left()
		ax.yaxis.set_major_formatter(plt.FormatStrFormatter('%0.2f'))
		ax.get_xaxis().tick_bottom()
		
		ax.set_title(rlabel, loc='left', y=1.09)


	outer_grid.tight_layout(fig)
	plt.show()

if __name__=='__main__':
	import argparse
	
	parser = argparse.ArgumentParser(description='Plot test accuracies for 3 subsets')
	
	parser.add_argument('--excl', type=str, default=None, help='Exclusion count per column(separated by commas)') 
	parser.add_argument('title', help="title")
	parser.add_argument('rlabels', help="row labels for each model. (separated by commas)")
	parser.add_argument('model_dir', help="model directory")
	parser.add_argument('model', help="model directory")
	parser.add_argument('model_versions', nargs='+', help="stat files generated by each model")

	args = parser.parse_args()
	
	plot(args.model_dir, args.model, args.model_versions, args.rlabels.split(','), args.title)
