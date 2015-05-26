import pandas as pd
import numpy as na
import os
from collections import defaultdict
from svphase.utils.common import sort_chr
from svphase.utils.config import PRECISION

def is_not_same(v):
  t = (len(v)-1)/2
  return v[:t]!=v[-t:]

class ConfusionMatrix():
  def __init__(self,df):
    self.df = df
    self.cmat = None
    self.cmat_map = defaultdict(list)
  def get_cmat(self, score_flag=False):
    if not self.cmat is None:
      return self.cmat
    truth_v = set(list(self.df['truth']))
    pred_v = set(list(self.df['pred']))
    cmat = pd.DataFrame(na.zeros((len(truth_v), len(pred_v))), index=sorted(list(truth_v)), columns=sorted(list(pred_v)), dtype=int)
    for idx,s in self.df[['truth','pred']].iterrows():
      cmat.ix[s['truth'],s['pred']] += 1
      self.cmat_map[(s['truth'],s['pred'])].append(idx)

    cmat['Total'] = cmat.sum(axis=1)
    cmat = cmat.T
    cmat['Total'] = cmat.sum(axis=1)
    cmat = cmat.T
    for c in cmat.columns:
      if c=='Total': continue
      cmat['pr-%s'%c] = cmat[c].astype(float)/cmat['Total'].astype(float)
    self.cmat = cmat
    return cmat
  def get_score_distr_grouped(self, method='truth'):
    assert method in set(['truth','pred','diff'])
    assert not self.cmat_map is None
    l = []
    for k,idx in sorted(self.cmat_map.iteritems()):
      t = self.df.ix[idx,k[0]]
      p = self.df.ix[idx,k[1]]
      if method == 'truth':
        l.append((k,t.values))
      elif method == 'pred':
        l.append((k,p.values))
      elif method == 'diff':
        l.append((k,t.values-p.values))
    return l
  def get_score_distr(self, func=None, method='truth'):
    assert method in set(['truth','pred','diff'])
    assert not self.cmat_map is None
    l = []
    for k,idx in sorted(self.cmat_map.iteritems()):
      t = self.df.ix[idx,k[0]]
      p = self.df.ix[idx,k[1]]

      if method == 'truth':
        v = t.values
      elif method == 'pred':
        v = p.values
      elif method == 'diff':
        v = t.values-p.values

      v = pd.Series(v, index=idx)
      if func is None:
        nidx = idx
      else:
        nidx = v.index[func(v)]   
      
      if len(nidx)>0:
        l.append(self.df.ix[nidx,:])
        l[-1]['err_label'] = "{0}-{1}".format(*k)
        l[-1]['score'] = v[nidx]
    m = pd.concat(l)
    if 'chr' in m.columns:
      return sort_chr(m)
    return m

  def plot_score_distr(self, method='truth', bin_count=50, savefig=None):
    from matplotlib import pyplot as plt, gridspec
    
    df = self.get_score_distr(method=method)
    df = df.ix[abs(df['score'])>PRECISION,:]
    truth_g = df.groupby('truth')
    fig = plt.figure(facecolor=None, edgecolor=None, linewidth=0, frameon=False)
    gs = gridspec.GridSpec(len(truth_g),1)
    axes = []
    
    cmap = plt.cm.jet
    pred_c = ['pA','nan','pB', 'homA','inc']
    cdict = dict(zip(pred_c, map(cmap, na.linspace(0,1,len(pred_c)))))
    if method=='diff' or method=='truth': 
      bins = na.linspace(-13, 0, bin_count)
    elif method=='pred':
      bins = na.linspace(-3, 0, bin_count)
   
    #bins = na.linspace(df['score'].min(), df['score'].max(), bin_count)

    lbl_fix = dict(zip(pred_c,pred_c))
    lbl_fix['homA'] = 'hom.'
    lbl_fix['inc'] = 'incorrect' 

    for i, (k,sdf) in enumerate(sorted(truth_g)):
      if len(axes)==0:
        ax = plt.subplot(gs[i,0])
      else:
        ax = plt.subplot(gs[i,0], sharex=axes[0])
      ax.set_title(k)
      err_g = sdf.groupby('pred')
      max_count = 0
      for err, score_distr in err_g:
       counts, bins2, patches = ax.hist(score_distr['score'], bins=bins, color=cdict[err], label=lbl_fix[err], lw=0, alpha=0.4)
       ax.plot(bins[:-1]+(bins[1:]-bins[:-1])/2,counts, color=cdict[err], lw=2, ls='--')
       max_count = max(max_count, max(counts))

      ax.set_ylabel('Count')
      if len(axes)==0:
        plt.setp(ax.get_xticklabels(), visible=False)
      ax.legend(loc='upper left')
      ax.yaxis.set_ticks(na.arange(0,max_count+2))
      ax.yaxis.tick_right()
      ax.yaxis.set_label_position("right")
      axes.append(ax)


    ax.set_xlabel('LogProb: {0}'.format(method))
    #plt.legend(loc='upper left')
    fig.tight_layout()
    if savefig is None:
      plt.show()
    else:
      plt.savefig(savefig)

  def latex_cmat(self, title):
    assert not self.cmat is None
    cols = ['homA', 'inc', 'pA', 'pB']
    for i in cols:
      if i in self.cmat.columns: continue
      self.cmat[i] = 0
      self.cmat['pr-%s'%i] = 0

    print "\\hline \\hline"
    print "{0} & {1} \\\\".format(title," & ".join(cols))
    print "\\hline"
    for idx, s in self.cmat.iterrows():
      print "{0} & {1} \\\\".format(idx, " & ".join(["{0:d} ({1:.2f})".format(int(s[c]),s['pr-%s'%c]) for c in cols]))
       
 
def cmat_from_bed_to_latex(bed_fpath):
  #df = pd.read_csv(bed_fpath, sep='\t', names=['chr','start','end','truth','pred'])
  df = pd.read_csv(bed_fpath, sep='\t', header=0)
  df = df.dropna(subset=['pred',])
  w = ConfusionMatrix(df)
  cmat = w.get_cmat(score_flag=True) 
  fpath_tokens = bed_fpath.split(os.sep)
  data_set = fpath_tokens[-2]
  #if fpath_tokens[2]=='gs':
  #  data_set = "Mills-gs"
  #else:
  #  data_set = "Mills"

  w.latex_cmat(title="{0} {1}".format(fpath_tokens[-1].split('-')[1], data_set))

def cmat_from_txt_to_latex(bed_fpath):
  #df = pd.read_csv(bed_fpath, sep='\t', names=['chr','start','end','truth','pred'])
  df = pd.read_csv(bed_fpath, sep='\t', header=0)
  df = df.dropna(subset=['pred',])
  w = ConfusionMatrix(df)
  cmat = w.get_cmat(score_flag=True) 
  #print cmat
  fpath_tokens = bed_fpath.split(os.sep)
  #fpath_tokens = os.path.basename(bed_fpath).split('-')
  ml = fpath_tokens[-1][:-4].split('-')[-1]
  train_set = fpath_tokens[-3]
  test_set = fpath_tokens[-2]

  w.latex_cmat(title="{0} {1}|{2}".format(ml, train_set, test_set))


def cmat_from_bed(bed_fpath):
  df = pd.read_csv(bed_fpath, sep='\t', header=0)
  df = df.dropna(subset=['pred',])
  w = ConfusionMatrix(df)
  cmat = w.get_cmat(score_flag=True) 
  #print cmat
  #for k,v in w.get_score_distr_grouped(method='diff'):
  #  print k
  #  print v
  diff = w.get_score_distr(method='diff')

  diff = diff.ix[diff['err_label'].map(is_not_same),:]
  diff.to_csv("{0}-diff-err.txt".format(bed_fpath[:-4]), sep='\t', header=True, index=False)
  #print w.get_score_distr(method='diff')[['chr','start','end', 'err_label', 'score']]
  #print w.get_score_distr(method='diff')[['err_label', 'score']]
  #w.plot_score_distr(method='diff', savefig="{0}-diff.png".format(bed_fpath[:-4]))
  #w.plot_score_distr(method='pred', savefig="{0}-pred.png".format(bed_fpath[:-4]))
  #w.plot_score_distr(method='truth', savefig="{0}-truth.png".format(bed_fpath[:-4]))

def cmat_per_approach(bed_fpath, approach_bed):
  #pred = pd.read_csv(bed_fpath, sep='\t', names=['chr','start','end','truth','pred']).set_index(['chr','start','end'])
  pred = pd.read_csv(bed_fpath, sep='\t', header=0).set_index(['chr','start','end'])
  approach_cmp = pd.read_csv(approach_bed, sep='\t', header=0).set_index(['chr','start','end'])
  pred_approach = pd.concat([pred,approach_cmp],axis=1)
  for approach in approach_cmp.columns:
    pred_approach.ix[pred_approach[approach]=='.',approach] = na.nan
    df = pred_approach[['pred',approach]].rename(columns={approach:'truth'}).dropna()
    cmat_obj = ConfusionMatrix(df)
    cmat = cmat_obj.get_cmat()
    if cmat.ix['Total','Total']>80:
      print approach
      print cmat

if __name__ == '__main__':
  import sys
  if len(sys.argv)==2:
    #cmat_from_bed(sys.argv[1])
    cmat_from_bed_to_latex(sys.argv[1])
    #cmat_from_txt_to_latex(sys.argv[1])
  elif len(sys.argv)==3:
    cmat_per_approach(sys.argv[1], sys.argv[2])
  else:
    print >>sys.stderr, "Nothing to do..."
