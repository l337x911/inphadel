import sys
import numpy as na
from svphase.learn.predict import EvalSim
from svphase.learn.classify import phased_classes_break_error

class EvalSimPhaseError(EvalSim):
  def __init__(self, clf_fpath, combines, start_size_range, end_size_range, right_break_error, left_break_error):
    self.start_size_range = start_size_range
    self.end_size_range = end_size_range
    self.right_break_error = right_break_error
    self.left_break_error = left_break_error
    EvalSim.__init__(self,clf_fpath, combines)
  def _get_data_collection(self):
    return phased_classes_break_error(self.start_size_range, self.end_size_range, self.right_break_error, self.left_break_error)
  def _predict(self):
    for data in self._get_data_collection().data:
      self.predr.add_data(data)
    self.pred = self.predr.predict(combines=self.combines, format_flag=False, nonzero_flag=False)
  def _truth(self):
    pclf2 = self._get_data_collection()
    pclf2.get_features(combines=self.combines, nonzero_flag=False)
    self._set_nonzero(pclf2)
    self.truth = pclf2.labels

def test(clf_fpath, combines):
  if combines == 'None':
    combines = None
  evaluator = EvalSimPhaseError(clf_fpath, combines, 1,5,0,0)
  print evaluator.get_accuracy(['pA','pB'], with_nan_flag=True)

def grid(clf_fpath, combines, n=3):
  if combines == 'None':
    combines = None
  assert n<=10 and n>0
  shifts = na.arange(0,100*n,100)

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

if __name__ == '__main__':
  clf_fpath = sys.argv[1]
  feature_subset = sys.argv[2]
  n = int(sys.argv[3])

  grid(clf_fpath, feature_subset, n)
  #test(clf_fpath, feature_subset)
  
