import os
import numpy as np
import pandas as pd
import matplotlib as mpl
from matplotlib import pyplot as plt
from svphase.utils.config import FONT, COLORS
from svphase.learn.evaluation import ClassLabel

def _get_model_and_version(model_stat):
	fname = os.path.basename(model_stat).split('.')
	return fname[0], '.'.join(fname[1:3])

def plot(labels, truth_bed, csv_pred1, csv_pred2, fig_fpath=None):
	label_obj = ClassLabel()
	assert len(labels) == 2
	index_cols = ['contig','start','end']
	labels = [l.replace('\\n','\n') for l in labels]
	print labels

	truth_df = pd.read_csv(truth_bed, sep='\t', index_col=None, header=None, names=['contig', 'start', 'end','truth'])
	pred_df1 = pd.read_csv(csv_pred1, sep='\t', index_col=None, header=0).set_index(index_cols).rename(columns={'best':'best1'})
	pred_df2 = pd.read_csv(csv_pred2, sep='\t', index_col=None, header=0).set_index(index_cols).rename(columns={'best':'best2'})

	df = pd.concat([truth_df.set_index(index_cols), pred_df1.loc[:,['best1',]], pred_df2.loc[:,['best2',]]], axis=1)
	df = df.loc[df['truth'].apply(label_obj.is_class),:]

	# Only include items in predictions
	df = df.loc[np.logical_and(pd.notnull(df['best1']), pd.notnull(df['best2'])),:]

	df['a'] = df['truth']==df['best1']
	df['b'] = df['truth']==df['best2']
	df['both'] = np.logical_and(df['a'],df['b'])	

	df['a'] = np.logical_and(df['a'], np.logical_not(df['both']))
	df['b'] = np.logical_and(df['b'], np.logical_not(df['both']))

	#for s in df.iterrows():
	#	print s
	gdf = df.groupby('truth').sum()
	gdf['total'] = gdf.sum(axis=1)
	print gdf
	# normalize by class
	ndf = gdf.loc[:,:]
	ndf = ndf.apply(lambda s:s[['a','b','both']]/s['total'], axis=1)

	print gdf[['a','b','both']].sum(axis=0)/float(gdf['total'].sum())
	
	ndf['total'] = gdf['total']/gdf['total'].sum().astype(float)
	print ndf	

	fig = plt.figure(figsize=(6.5,3), facecolor='white', edgecolor='none')
	mpl.rc('font', **FONT)
	ax = fig.add_subplot(111)

	order = [c for c in label_obj.classes if c in gdf.index]

	rename_labels = dict(zip(order, order))
	rename_labels['hom']='hom.'
	rename_labels['inc']='incorrect'
	
	class_labels = map(rename_labels.get, order)
	white_space = 0.001

	width = gdf.loc[order,'total']
	cumsum = width.cumsum()
	width /= float(cumsum[-1])
	print width

	left = np.zeros(len(width), dtype=float)
	left[1:] = width.cumsum()[:-1]
	left += white_space
	width = width-2*white_space

	ax.bar(left, ndf.loc[order,'b'].values, bottom=0, width=width, color=COLORS['b'], linewidth=0) 
	ax.bar(left, ndf.loc[order,'both'].values, bottom=ndf.loc[order,'b'].values, width=width, color=COLORS['both'], linewidth=0) 
	ax.bar(left, ndf.loc[order,'a'].values, bottom=ndf.loc[order,'b'].values+ndf.loc[order,'both'].values, width=width, color=COLORS['a'], linewidth=0) 

	for px, w, c,c_w in zip(left, width, class_labels, gdf.loc[order,'total']):
		ax.text(px+ w/2, -.02, c, ha='center', va='top', size='large')
		#ax.text(px+ w/2, 1.0, "%d"%c_w, ha='center', va='bottom', size='small')
	
	pad = 0.01
	ha = 'center'
	tcolor = 'w'
	tsize = 'medium'
	tfmt = "{0:0.2f}"
	for px, w, c in zip(left, width, order):
		ax.text(px+ w/2, 1-ndf.loc[c,'a']+pad, tfmt.format(ndf.loc[c,'a']), ha=ha, va='bottom', size=tsize, color=tcolor)
		ax.text(px+ w/2, ndf.loc[c,'b']+ndf.loc[c,'both']/2, tfmt.format(ndf.loc[c,'both']), ha=ha, va='center', size=tsize, color=tcolor)
		ax.text(px+ w/2, ndf.loc[c,'b']-pad, tfmt.format(ndf.loc[c,'b']), ha=ha, va='top', size=tsize, color=tcolor)
		

	py = np.linspace(0.2,0.8, num=3)
	px = 1.01
	ax.text(px, py[0], labels[1], ha='left', va='center',color=COLORS['b'], size='medium')
	ax.text(px, py[1], 'Both', ha='left', va='center',color=COLORS['both'], size='medium')
	ax.text(px, py[2], labels[0], ha='left', va='center',color=COLORS['a'], size='medium')
	
	ax.set_ylabel('Fraction of correct predictions')
	ax.set_xticks([])
	ax.set_ylim(0, 1.0)
	print ax.get_xlim()
	ax.set_frame_on(False)
	ax.get_yaxis().tick_left()
	ax.yaxis.set_major_formatter(plt.FormatStrFormatter('%0.2f'))
	#ax.get_xaxis().tick_bottom()
	ax.axes.get_xaxis().set_visible(False)
	fig.tight_layout()
	fig.subplots_adjust(right=0.82, bottom=0.12)
	
	if fig_fpath is None:
		plt.show()
	else:
		plt.savefig(fig_fpath)

if __name__=='__main__':
	import argparse
	
	parser = argparse.ArgumentParser(description='Plot number of correct predictions per class for two classes')
	
	parser.add_argument('-o', dest='save', default=None, help='Save figure')
	parser.add_argument('labels', help="Labels for each model separated by commas")
	parser.add_argument('truth', help="Truth bed file")
	parser.add_argument('pred1', help="Prediction csv file model 1")
	parser.add_argument('pred2', help="Prediction csv file model 2")

	args = parser.parse_args()
	plot(args.labels.split(','), args.truth, args.pred1, args.pred2, fig_fpath=args.save)
