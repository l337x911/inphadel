import os
import numpy as np
import pandas as pd
import matplotlib as mpl
from  matplotlib import pyplot as plt, gridspec
from operator import itemgetter
from collections import defaultdict

break_error_ranges = [(1,5),(5,10),(10,300)]

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100, N=100):
	new_cmap = mpl.colors.LinearSegmentedColormap.from_list('trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval), cmap(np.linspace(minval, maxval, n)),N)
	return new_cmap

def _get_acc_dict_for_truth_fpath(truth_fpath, n):
	shifts = np.arange(0,100*n,100)
	truth_dir = os.path.dirname(truth_fpath)	
	truth = pd.read_csv(truth_fpath, sep='\t', names=['contig','start','end','truth'])

	fmt = "{d}/a{a:03d}.b{b:03d}.bed.csv"

	acc_dict = {}
	for i,a in enumerate(shifts):
		for j,b in enumerate(shifts):
			shift_df = pd.read_csv(fmt.format(d=truth_dir, a=a, b=b), sep='\t', index_col=None, header=0).drop('contig', axis=1).rename(columns={'start':'start_err','end':'end_err'})

			acc_df = pd.concat([truth, shift_df], axis=1)
			assert ((acc_df['start']-acc_df['start_err'])==a).all()
			assert ((acc_df['end_err']-acc_df['end'])==b).all()

			acc_df = acc_df.reset_index()
			acc_df['len'] = acc_df['end']-acc_df['start']
			acc_df['correct'] = acc_df['truth'] == acc_df['best']
			acc_dict[(i,j)]	= acc_df
	#print acc_dict
	return acc_dict

def compute_accuracy(truth_fpaths, n=3):
	shifts = np.arange(0,100*n,100)
	acc_grid = {}

	assert len(truth_fpaths)>0
	
	acc_dict = _get_acc_dict_for_truth_fpath(truth_fpaths[0], n)

	if len(truth_fpaths)>1:
		for truth_fpath in truth_fpaths[1:]:
			t_acc_dict = _get_acc_dict_for_truth_fpath(truth_fpath, n)
			for k in t_acc_dict.iterkeys():
				acc_dict[k] = pd.concat([acc_dict[k],t_acc_dict[k]], ignore_index=True)

	for r1,r2, in break_error_ranges:
		acc_grid_range = np.zeros((n,n),dtype=np.float)
		for i,a in enumerate(shifts):
			for j,b in enumerate(shifts):
				acc_df = acc_dict[(i,j)]
				fidx = np.logical_and(acc_df['len']>=(r1*1000),acc_df['len']<(r2*1000))
				acc_df = acc_df.loc[fidx,:]

				acc_grid_range[i,j] = acc_df['correct'].sum()/float(len(acc_df.index))
		acc_grid[(r1,r2)] = acc_grid_range		
		print "Accuracy for Deletions in {0}kb-{1}kb".format(r1,r2)
		print acc_grid_range
	return acc_grid

def plot_heatmap(npz, h_r=12):
	#npz = np.load('/home/anand/Projects/assembly/src/break_error_acc.npz')

	fig = plt.figure(figsize=(11,6),facecolor=None, edgecolor=None)
	
	gs = gridspec.GridSpec(2,3, height_ratios=[h_r,1])
	#gs.update(left=0.05, bottom=0.05, top=0.95, right=0.95, wspace=0.01, hspace=0.01)

	z_min, z_max = 0,1
	#cmap = truncate_colormap(plt.cm.YlGnBu,0.2,0.95,N=50)
	cmap = plt.cm.YlGnBu
	 
	for idx, (r, a) in enumerate(zip(break_error_ranges, map(itemgetter(1), sorted(npz.items())))):
		ax = plt.subplot(gs[0,idx], aspect='equal')
		offsets = map('{0:d}'.format, np.arange(0,a.shape[0]*100,100))
		hm = ax.pcolor(a, cmap=cmap, vmin=z_min, vmax=z_max)

		# put the major ticks at the middle of each cell
		ax.set_xticks(np.arange(a.shape[0])+0.5, minor=False)
		ax.set_yticks(np.arange(a.shape[1])+0.5, minor=False)
		
		ax.set_xticklabels(offsets, minor=False)
		ax.set_yticklabels(offsets, minor=False)
		ax.set_xlabel('Deletions with size\n{0}kb-{1}kb'.format(*r))
		ax.set_ylabel('Offset from true break a')

		ax.invert_yaxis()
		ax.xaxis.tick_top()    
	
	cax =  plt.subplot(gs[1,:])
	cb = fig.colorbar(hm, cax=cax, orientation='horizontal', format='%.2f')  
	cax.set_xlabel('accuracy')
	fig.tight_layout()
	plt.show()

if __name__ == '__main__':
	import argparse 
	
	parser = argparse.ArgumentParser(description="Plots heatmap of accuracy for aXXX.bXXX.bed.csv error calls (see learn/break_error.py)")

	parser.add_argument('truths', nargs='+', help='bed(s) file containing actual calls')
	parser.add_argument('n', type=int, help='number of 100 deviations for a and b')
	
	args = parser.parse_args()
	assert args.n<10 and args.n>0
	break_error_grid_acc = compute_accuracy(args.truths, args.n)
	plot_heatmap(break_error_grid_acc)

