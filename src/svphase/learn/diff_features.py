from svphase.learn.classify import realistic_classes_simple, realistic_classes_wgs2x, realistic_classes_hic2x
from svphase.learn.predict import EvalReal, EvalScoreReal, EvalSim, EvalMix
from svphase.utils.config import PRECISION
import sys
from operator import itemgetter
import pandas as pd
import numpy as na
from collections import defaultdict
PRECISION = PRECISION * 0.000000001

pd.set_printoptions(max_columns=1500, max_rows=1500)

def count_nans(df):
  return na.sum((df < PRECISION).apply(all,axis=0))  

def get_nans(df):
  return ((df < PRECISION).apply(all,axis=0))  

def nans_per_feat_set(df):
  g = df.T.groupby(level='data_type')
  nan_cats = g.aggregate(count_nans).ix[:,0]

  if 'hic_a' in g.groups and 'hic_ac' in g.groups:
    nan_cats = nan_cats.append(pd.Series(count_nans(pd.concat((g.get_group('hic_a'),g.get_group('hic_ac')))), index=['hic',]))
  
  alls = [t for t in g.groups if t != 'wgs']
  nan_cats = nan_cats.append(pd.Series(count_nans(pd.concat([g.get_group(t) for t in alls])), index=['all',]))
  #print nan_cats
  #sys.exit(1)
  return nan_cats

def to_latex(df, cols=None, title=''):
  if cols is None: cols = df.columns
  print "{0} & {1} \\\\".format(title," & ".join(cols))
  print "\\hline"
  for idx, s in df.iterrows():
    print "{0} & {1} \\\\".format(idx, " & ".join(["{0:d}".format(int(s[c])) for c in cols]))

def nans_in_set(feats):
  g = feats.T.groupby(level='data_type')
  f = g.aggregate(get_nans)
  return f

def naive_classifier(feats):
  t = feats.T
  #t['default'] = na.arange(len(t.index))
  #cols_names = ['default',]+t.index.names
  cols_names = t.index.names
  t = t.reset_index()
  t = t.drop(['orient', 'orient1', 'orient2'], axis=1)
  #t = t.T
  #print "Columns:", t.columns, "Index:", t.index
  g = t.groupby(['data_type','wtype','allele'])

  wgs_cov = []
  wgs_around = []
  new_feats = pd.DataFrame([na.zeros(len(feats.index))], index=['wgs-cov',]).T
  allele_support = defaultdict(list)
  for (dtype,wtype,allele),values in g:
    values = values.drop(['data_type','wtype','allele'], axis=1)
    if dtype=='wgs':
      if wtype=='cov':
        new_feats['wgs-cov'] = values.values.ravel()
        wgs_cov.append(values)
      if wtype=='w.around_w.around':
        wgs_around.append(values)
        new_feats['wgs-around'] = values.values.ravel()
    else:
      if wtype=='w.around_w.around':
        allele_support[(allele,-1)].append(values)
      else:
        allele_support[(allele,1)].append(values)
 
  for (allele,support),v in allele_support.iteritems():
    new_feats['%d*%s'%(support, allele)]= pd.concat(v).sum(axis=0).values.ravel()
  new_feats['A-B'] = new_feats['1*A']-new_feats['1*B']
  new_feats['npred'] = 'inc'
  new_feats.ix[new_feats['A-B']>0.001,'npred']='pB'
  new_feats.ix[new_feats['A-B']<-0.001,'npred']='pA'
  new_feats.ix[new_feats['wgs-cov']<1, 'npred']='homA'
  #print feats.ix[:,na.logical_and(feats.T['data_type']=='wgs',feats.T['wtype']!='w.around_w.around')].T
  return new_feats

def get_break_info(p):
  #print p.data[0].wgs.set_index(['data_type', 'wtype', 'orient1', 'orient2'])

  break_info_data = []
  for idx, d in enumerate(p.data):
    nd = d.wgs.set_index(['data_type', 'wtype', 'orient1', 'orient2']).T
    nd['class_label'] = d.class_label
    nd['data_idx'] = idx/4
    break_info_data.append(nd)
  #print break_info_data
  #return ""
  return pd.concat(break_info_data,axis=0, ignore_index=True).T

def diff_feature_set(df1,df2,data_type):
  #print df1.index
  #print df2.index
  g1 = df1.T.groupby(level='data_type')
  g2 = df2.T.groupby(level='data_type')
  f1 = g1.aggregate(get_nans)
  f2 = g2.aggregate(get_nans)
  diff_idx = (f2.T[[data_type,]]!=f1.T[[data_type,]]).T.values.ravel()
  print df1.index[diff_idx]
  #print g1.get_group(data_type).shape
  print g1.get_group(data_type).ix[:,df1.index[diff_idx]]
  print g2.get_group(data_type).ix[:,df2.index[diff_idx]]
  
  return g2.get_group(data_type)-g1.get_group(data_type)

def sim_test(evidence, classifier, nfolds):
  #clf_fpath = '/home/anand/Projects/assembly/data/simulations/models/sim.del-chr19-20.5norm.%02d.%s.pkl'%(nfolds, classifier)
  #clf_only_fpath = '/home/anand/Projects/assembly/data/simulations/models/sim.del-chr19-20-hm.10-2c.%s.%02d.%s.pkl'
  #r1 = EvalReal(clf_fpath, evidence)
  p1 = realistic_classes_simple(False, classifier=classifier, evidence=evidence, nfolds=nfolds)
  p2 = realistic_classes_wgs2x(False, classifier=classifier, evidence=evidence, nfolds=nfolds)
  p3 = realistic_classes_hic2x(False, classifier=classifier, evidence=evidence, nfolds=nfolds)
  
  #r1_f = r1.predr.feats
  #print r1_f.T.ix[:,::2]

  p1_f, p1_l = p1.get_features(evidence, nonzero_flag=False) 
  p2_f, p2_l = p2.get_features(evidence, nonzero_flag=False) 
  p3_f, p3_l = p3.get_features(evidence, nonzero_flag=False) 
  #assert (p1_l==p2_l).all()
  p1_l = pd.DataFrame(p1_l)
  p2_l = pd.DataFrame(p2_l)  
  #print p1_l.groupby(0).apply(len)
  #print p1_l.T.ix[:,::2]
  #print p1_f.T.ix[:,::2]
  #print p2_l.groupby(0).apply(len)
  #print p2_l.T.ix[:,::2]
  #print p2_f.T.xs('wgs').ix[:,:10]
  #print (p2_f-p1_f).T.ix[:,:90]
  hets = na.logical_or(p1_l.values==1, p1_l.values==3).reshape((p1_l.values.size,))
  p1_f_h = p1_f.ix[hets,:]
  p1_f_h = p1_f_h.rename(dict(zip(p1_f_h.index, p2_f.index)))
  nan_cols = nans_per_feat_set(p1_f_h).index
  print get_break_info(p2)
  print p1_f_h.shape
  nan_counts = pd.DataFrame([nans_per_feat_set(p1_f_h), nans_per_feat_set(p2_f),  nans_per_feat_set(p3_f)], index=['normal','wgs2x','hic2x'], columns=nan_cols)
  
  print to_latex(500-nan_counts, cols=['all','wgs_a','hic'])
  #print p2_f.mean(axis=0).astype(float), p1_f_h.mean(axis=0).astype(float)
  #print (p2_f.mean(axis=0).astype(float)/p1_f_h.mean(axis=0).astype(float))
  #print diff_feature_set(p1_f_h,p2_f,'wgs_a')
  feats = p1.feats.copy()
  nans_df =  nans_in_set(feats)
  nans_idx = na.logical_and(nans_df.T['wgs_a'],nans_df.T['hic_a'])

  idx_to_class = dict([(j,i) for i,j in p1.class_str.items()])
  print feats.shape, len(p1.labels)
  feats['class'] = [idx_to_class[i] for i in p1.labels]
  print "Number of Nans:", na.sum(nans_idx)
  nans_c =  pd.Series(feats.ix[nans_idx,'class']).value_counts()
  groups_c =  pd.Series(feats.ix[:,'class']).value_counts()

  print "Nans;\n", nans_c
  print "Used;\n", groups_c-nans_c
  print "Used Sum:", (groups_c-nans_c).sum()
  print "Groups;\n", groups_c


def real_comp_test(classifier, nfolds):
  clf_only_fpath = '/home/anand/Projects/assembly/data/simulations/models/sim.del-chr19-20-hm.10-2c.%s.%02d.%s.pkl'
  clf_fpath = '/home/anand/Projects/assembly/data/simulations/models/sim.del-chr19-20-hm.10-2c.%02d.%s.pkl'
  sv_fpath='/home/anand/Projects/assembly/data/mills_na12878-del-gs.NA12878.anyp.bed'
  
  p = EvalReal(clf_fpath%(nfolds, classifier), None, sv_fpath=sv_fpath)
  feats = p.predr.feats.copy()
  idx_to_class = dict(p.predr.class_int.items())
  nans_df =  nans_in_set(feats)
  feats['class'] = [idx_to_class[i] for i in p.truth]
  feats['pred'] = [idx_to_class[i] for i in p.pred]
  
  hic = EvalScoreReal(clf_only_fpath%('hic_only', nfolds, classifier), 'hic_only', sv_fpath=sv_fpath)
  feats_hic = hic.predr.feats.copy()
  d = pd.DataFrame(hic.lpred, columns=map(itemgetter(1), sorted(hic.predr.class_int.items())))
  d = d.fillna(-100)
  #print d
  feats['hic_pred'] = d.idxmax(axis=1)
  feats['hic_scores'] = d.max(axis=1)
  wgs = EvalReal(clf_only_fpath%('wgs_only', nfolds, classifier), 'wgs_only', sv_fpath=sv_fpath)
  feats_wgs = hic.predr.feats.copy()
  feats['wgs_pred'] = [idx_to_class[i] for i in wgs.pred]
  
  print nans_df.sum(axis=1)
  print feats.shape, nans_df.shape  

  #print feats['wgs_pred']
  hic_only_correct = na.logical_and(na.logical_and(feats['class']==feats['hic_pred'], na.logical_not(nans_df.T['hic_a'])), nans_df.T['wgs_a'])
  #print feats.ix[hic_only_correct,:]
  #print feats.shape, feats_hic.shape
 
  #print feats.columns
  
  #print get_break_info(p.predr).ix[:3,:]
  assert (get_break_info(p.predr).ix[:3,:].values == get_break_info(hic.predr).ix[:3,:].values).all()
  t =  get_break_info(p.predr).ix[:3,:].T
  for c in ['class', 'hic_pred', 'hic_scores', 'wgs_pred']:
    t[c] = feats[c]

  print t.shape
  hic_only_correct_bed = t.ix[hic_only_correct,:].sort(column='hic_scores', ascending=False)
  print hic_only_correct_bed
  hic_only_correct_bed.to_csv('/home/anand/Projects/assembly/data/mills_na12878-del-gs.NA12878.anyp.validate.bed', index=False, header=['contig', 'start', 'end', 'truth', 'hic_pred', 'hic_scores', 'wgs_pred'], sep='\t')

  #pf, pl = p.get_features(None, nonzero_flag=False)
  #pl = pd.DataFrame(pl)
  #print nans_per_feat_set(feats)
  #hic_p = EvalReal(clf_fpath%('hic_only', nfolds, classifier), 'hic_only', sv_fpath=sv_fpath)
  #wgs_p = EvalReal(clf_fpath%('wgs_only', nfolds, classifier), 'wgs_only', sv_fpath=sv_fpath)

  #hic_pf, hic_pl = hic_p.get_features('hic_only', nonzero_flag=False)
  #wgs_pf, wgs_pl = wgs_p.get_features('wgs_only', nonzero_flag=False)

  #hic_pl = pd.DataFrame(hic_pl)
  #wgs_pl = pd.DataFrame(wgs_pl)

def real_test(evaluator, evidence, classifier, nfolds):
  if not evidence is None:
    clf_fpath = '/home/anand/Projects/assembly/data/simulations/models/sim.del-chr19-20-hm.10-2c.%s.%02d.%s.pkl'%(evidence, nfolds, classifier)
  else:
    clf_fpath = '/home/anand/Projects/assembly/data/simulations/models/sim.del-chr19-20-hm.10-2c.%02d.%s.pkl'%(nfolds, classifier)
  
  sv_fpath='/home/anand/Projects/assembly/data/gm12878/hapmap3_r2_b36_fwd-truth.bed'
  #sv_fpath='/home/anand/Projects/assembly/data/mills_2010-08-table4.CEU1kb.NA12878.any.bed'
  #sv_fpath='/home/anand/Projects/assembly/data/mills_na12878-del-gs.NA12878.anyp.bed'
  r = evaluator(clf_fpath, clf_fpath.split('.')[-4], sv_fpath=sv_fpath)
  #r = evaluator(clf_fpath, clf_fpath.split('.')[-4])
  feats = r.predr.feats.copy()
  idx_to_class = dict(r.predr.class_int.items())

  naive = naive_classifier(feats)

  feats['naive'] = naive['A-B']
  feats['class'] = [idx_to_class[i] for i in r.truth]
  feats['pred'] = [idx_to_class[i] for i in r.pred]
  feats['npred'] = naive['npred']

  feats['nans'] = 'n'
  nans_df =  nans_in_set(feats)
  nans_idx = na.logical_and(nans_df.T['wgs_a'],nans_df.T['hic_a'])
  feats.ix[nans_idx,'nans'] = 'y'
  feats.ix[nans_idx,'pred'] = 'nan'
  feats.ix[nans_idx,'npred'] = 'nan'
  feats['sv_size'] = [b-a for c,a,b in r.predr.get_loci()]

  break_idx = get_break_info(r.predr)
  null_idx = lambda x:('wgs',x,'+','-')
  coi_idx = break_idx.ix[null_idx('c1'),:]=='chr13'
  #soi_idx = break_idx.ix[null_idx('a1'),:]>20000000
  #eoi_idx = break_idx.ix[null_idx('b1'),:]<20700000
  #aoi = na.logical_and(coi_idx,na.logical_and(soi_idx,eoi_idx))
  aoi = coi_idx
  #print break_idx.ix[:,aoi]
  #print feats.ix[399:404,:].T
  #print "No data counts"
  #print  nans_per_feat_set(feats)
  print "Number of Nans:", na.sum(nans_idx)
  nans_c =  pd.Series(feats.ix[nans_idx,'class']).value_counts()
  groups_c =  pd.Series(feats.ix[:,'class']).value_counts()

  print "Nans;\n", nans_c
  print "Used;\n", groups_c-nans_c
  print "Used Sum:", (groups_c-nans_c).sum()
  print "Groups;\n", groups_c
  
  print "Number of calls uniquely correct to naive"
  print na.sum(na.logical_and(na.logical_and(feats['npred']==feats['class'],feats['pred']!=feats['class']), na.logical_or(feats['npred']=='pA',feats['npred']=='pB')))

if __name__ == '__main__':
  evidence = sys.argv[1]
  classifier = 'RandomForest'
  #classifier = 'SVM'
  nfolds = 10
  if evidence == 'None':
    evidence = None 
  #real_comp_test(classifier, nfolds)
  sim_test(evidence, classifier, nfolds)
  #real_test(EvalReal, evidence, classifier, nfolds)
  #real_test(EvalMix, evidence, classifier, nfolds)

