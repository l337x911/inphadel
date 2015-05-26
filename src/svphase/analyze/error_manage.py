import pandas as pd
from svphase.analyze.confusion_matrix import ConfusionMatrix, is_not_same
import os
import numpy as na

def find_scores(bed_fpath):
  df = pd.read_csv(bed_fpath, sep='\t', header=0)
  df = df.dropna(subset=['pred',])
  w = ConfusionMatrix(df)
  cmat = w.get_cmat(score_flag=True) 
  #print cmat
  #for k,v in w.get_score_distr_grouped(method='diff'):
  #  print k
  #  print v
  diff = w.get_score_distr(method='diff')
  return diff


def ri(df):
  rn = dict([('level_%s'%i, v) for i,v in enumerate(['chr', 'start', 'end'])])
  return df.reset_index().rename(columns=rn)

def separate_errors(gs_flag=False):
  gs = '-gs' if gs_flag else ''

  #skeleton = '/home/anand/Dropbox/SV_Detection/Results/working/del-chr19-20-hm_5-2c/mills-any-score/10-{0}{1}-p2error.txt' 
  #error_prefix = '/home/anand/Dropbox/SV_Detection/Results/working/del-chr19-20-hm_5-2c/mills-any-score/error' 
  skeleton =     '/home/anand/Dropbox/SVPhasePaper/fig/del-chr19-20-hm_10-2c/mills-score-h/10-{0}{1}.txt' 
  error_prefix = '/home/anand/Dropbox/SVPhasePaper/fig/del-chr19-20-hm_10-2c/mills-score-h/error' 
  #skeleton =     '/home/anand/Dropbox/SVPhasePaper/fig/del-chr19-20-hm_10-2c/mills-anyp-score-h/10-{0}{1}.txt' 
  #error_prefix = '/home/anand/Dropbox/SVPhasePaper/fig/del-chr19-20-hm_10-2c/mills-anyp-score-h/error' 

  rf = find_scores(skeleton.format('RandomForest', gs)).set_index(['chr','start','end'])
  svm = find_scores(skeleton.format('SVM', gs)).set_index(['chr','start','end'])
  #print rf.shape
  #print svm.shape
  intersect_idx = rf.index.intersection(svm.index)
  #print intersect_idx.shape

  svm_rename = dict([(i,'SVM_%s'%i) for i in ('err_label','score')])
  rf_rename = dict([(i,'RF_%s'%i) for i in ('err_label','score')])
  rfi = rf[['err_label','score']].ix[intersect_idx,:].rename(columns=rf_rename)
  svmi = svm[['err_label','score']].ix[intersect_idx,:].rename(columns=svm_rename)
  cb = pd.concat([rfi,svmi], axis=1)
  cb = cb.ix[na.logical_or(map(is_not_same, cb['RF_err_label']),map(is_not_same, cb['SVM_err_label'])),:]
  cbe = cb.ix[na.logical_and(map(is_not_same, cb['RF_err_label']),map(is_not_same, cb['SVM_err_label'])),:]
  ri(svm.ix[svm.index.drop(cbe.index),:]).to_csv('{0}/{1}-confident.txt'.format(error_prefix,os.path.basename(skeleton[:-4].format('SVM',gs))), sep='\t', index=False)
  ri(rf.ix[rf.index.drop(cbe.index),:]).to_csv('{0}/{1}-confident.txt'.format(error_prefix,os.path.basename(skeleton[:-4].format('RandomForest',gs))), sep='\t', index=False)


  same_err = cbe.ix[cbe['RF_err_label']==cbe['SVM_err_label'],:]
  diff_err = cb.ix[cb['RF_err_label']!=cb['SVM_err_label'],:]
 
  hom_err = same_err.ix[map(lambda x:x[3:]=='homA', same_err['RF_err_label']),:]
  inc_err = same_err.ix[map(lambda x:x[3:]=='inc', same_err['RF_err_label']),:]
  print same_err
  print diff_err
  #p2p_err = same_err.ix[map(lambda x:x[3, same_err['RF_err_label']),:]
 
  hom_err.reset_index().to_csv("{0}/diff{1}-hom.txt".format(error_prefix, gs), sep='\t', index=False) 
  inc_err.reset_index().to_csv("{0}/diff{1}-inc.txt".format(error_prefix, gs), sep='\t', index=False) 
  #print cb

  diff_err = diff_err.ix[na.logical_and(cb['RF_score']<-1,cb['SVM_score']<-1),:]
  diff_err.reset_index().to_csv("{0}/diff{1}-low.txt".format(error_prefix, gs), sep='\t', index=False) 

def rm_dup(df):
  g = df.groupby(['chr','start','end'])
  #return g.first().set_index(['chr','start','end']) 
  return g.agg(lambda x:x.ix[0])

def get_confident_df(df, tc_idx):
  confident = df
  truth_s = confident['truth']
  truth_s[tc_idx] = confident.ix[tc_idx,'pred']
  confident['truth'] = truth_s
  return confident

def report_consistent_calls(txt1_fpath, txt2_fpath):
  dfa = rm_dup(pd.read_csv(txt1_fpath, sep='\t', header=0))
  dfb = rm_dup(pd.read_csv(txt2_fpath, sep='\t', header=0))
  print dfa.shape
  calls = pd.notnull(dfa['pred'])
  calls_count = na.sum(calls)
  assert (calls==pd.notnull(dfb['pred'])).all() 
 
  print "Calls in A: {0}/{1}, calls {2:0.4f}, nans {3}".format(calls_count, dfa.index.size, float(calls_count)/dfa.index.size, dfa.index.size-calls_count)
  print "Histogram of A"
  print dfa.ix[calls,'pred'].value_counts()

  consistent_idx = dfa.ix[calls, 'pred']==dfb.ix[calls,'pred']
  truth_idx = dfa.ix[calls, 'pred']==dfa.ix[calls,'truth']
  tc_idx = na.logical_or(consistent_idx, truth_idx)
  print "Consistent: {0}/{1}, {2:0.4f}".format(na.sum(consistent_idx), consistent_idx.index.size, float(na.sum(consistent_idx))/consistent_idx.index.size)
  print "Consistent between A and B:"
  print dfa.ix[calls,'pred'][consistent_idx].value_counts()
  print "Truth or Consistent: {0}/{1}, {2:0.4f}".format(na.sum(tc_idx), tc_idx.index.size, float(na.sum(tc_idx))/tc_idx.index.size)
  print "Matching A and Truth or Consistent between A and B :"
  print dfa.ix[calls,'pred'][tc_idx].value_counts()

  confident_a = get_confident_df(dfa.ix[calls,:], tc_idx)
  confident_a.to_csv("{0}-confident.txt".format(txt1_fpath[:-4]), sep='\t', index=True)

  confident_b = get_confident_df(dfb.ix[calls,:], na.logical_or(consistent_idx, dfb.ix[calls,'pred']==dfb.ix[calls,'truth']))
  confident_b.to_csv("{0}-confident.txt".format(txt2_fpath[:-4]), sep='\t', index=True)

 
def get_rank(s):
  #return s[s['pred']]-max(s[s.index[2:]].drop([s['pred'],]))
  return s[s['pred']]-s[s['truth']]

 
def report_inconsistent_calls(txt_fpath):
  dfa = rm_dup(pd.read_csv(txt_fpath, sep='\t', header=0))
  calls = pd.notnull(dfa['pred'])
  dfa = dfa.ix[calls,:]
  dclasses = dfa.columns[2:]

  #dfa['rank'] = dfa.apply(lambda s: s[s['pred']]-max(s[dclasses].drop([s['pred'],])))
  dfa['rank'] = dfa.apply(get_rank, axis=1)

  dfa = dfa.ix[dfa['truth']!=dfa['pred'],:]
  for idx, s in dfa.sort(['rank',], ascending=False).iterrows():
    idx_str ="{0}:{1}-{2}".format(*idx)
    print  "{0} & {1}-{2} & {3} \\\\".format(idx_str, s['truth'],s['pred'], s['rank'])
  

if __name__ == '__main__':
  import sys
  #report_inconsistent_calls(sys.argv[1])
  #report_consistent_calls(sys.argv[1], sys.argv[2])
  
  separate_errors(gs_flag=False)
