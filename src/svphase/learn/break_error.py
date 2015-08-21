import sys
import os
import pandas as pd
import numpy as np

def create_grid(truth_fpath, n=3):
	shifts = np.arange(0,100*n,100)
	truth_dir = os.path.dirname(truth_fpath)	
	truth = pd.read_csv(truth_fpath, sep='\t', names=['contig','start','end','truth'])

	fmt = "{d}/a{a:03d}.b{b:03d}.bed"
	for i,a in enumerate(shifts):
		for j,b in enumerate(shifts):
			d = truth[['contig','start','end']]
			d.loc[:,'start'] = d['start'].astype(int)-a
			d.loc[:,'end'] = d['end'].astype(int)+b
			d.to_csv(fmt.format(d=truth_dir, a=a, b=b), sep='\t', index=False, header=False)


"""
def split_accuracy_by_size():
  size_bins = [(1,5),(5,10),(10,0)]
  acc_grids = []
  for r1,r2 in size_bins :
    acc_grid = na.zeros((n,n),dtype=na.float)
    for i,v1 in enumerate(shifts):
      for j,v2 in enumerate(shifts):
        evaluator = EvalSimPhaseError(clf_fpath, combines, r1,r2,v1,v2)
        acc_grid[i,j] = evaluator.get_accuracy(['pA','pB'], with_nan_flag=True)
    acc_grids.append(acc_grid)

  for (r1,r2),acc_grid in zip(size_bins, acc_grids):
    print "Accuracy for Deletions in {0}kb-{1}kb".format(r1,r2)
    print acc_grid  
  na.savez("break_error_acc.npz", *acc_grids)
"""

if __name__ == '__main__':
	import argparse
	
	parser = argparse.ArgumentParser(description="Produce aXXX.bXXX.bed error calls")

	parser.add_argument('truth', help='bed file containing real calls')
	parser.add_argument('n', type=int, help='number of 100 deviations for a and b')
	
	args = parser.parse_args()
	assert args.n<10 and args.n>0
	create_grid(args.truth, args.n)
  
