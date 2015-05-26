import pandas as pd
from itertools import chain
import numpy as na
from collections import defaultdict
from operator import itemgetter
from svphase.analyze.hm_trios import SNPTransmit
from svphase.utils.common import sort_chr
import sys

def get_mids_with_max_sample_labels(cs,samples):
  # max of computational_approaches
  g_id = cs.groupby(['MERGED_ID']+samples)
  approach_id_mid = g_id.apply(lambda x:len(x['TYPE_OF_COMPUTATIONAL_APPROACH'].value_counts()))
  mid_to_inherit = approach_id_mid.reset_index().groupby('MERGED_ID').apply(lambda s:s[samples+[0,]].set_index(samples).idxmax())
  for idx, s in enumerate(samples):
    mid_to_inherit[s] = mid_to_inherit[0].map(itemgetter(idx))
  return mid_to_inherit.drop([0,], axis=1)

def get_mids_with_any_sample_labels(cs,samples):
  # max of computational_approaches
  g_id = cs.groupby(['MERGED_ID'])
  approach_id_mid = g_id.apply(lambda x:(x[samples].sum()>0).astype(int))
  return approach_id_mid
  #mid_to_inherit = approach_id_mid.reset_index().groupby('MERGED_ID').apply(lambda s:s[samples+[0,]].set_index(samples).idxmax())
  #return mid_to_inherit

def algorithm_to_dels(approach, all_dels=False, hets_only=False, somatic_flag=False):
  supp_csv_fpath = '/home/anand/Projects/assembly/data/mills_2010-08-table4.CEU.csv'
  pm_to_m_map_fpath = '/home/anand/Projects/assembly/data/mills_2010-08-table4.CEU1kb.map'
  out_bed_fpath = '/home/anand/Projects/assembly/data/mills_2010-08-table4.CEU1kb.NA12878.{0}.bed'.format(approach)
  m_sites = pd.read_csv( '/home/anand/Projects/assembly/data/trio.2010_06.deletions.sites.1kb.bed',sep='\t', names=['chr','start','end','MERGED_ID'], index_col=None).set_index('MERGED_ID')

  alg = pd.read_csv(supp_csv_fpath, sep='\t', header=0, index_col=None).set_index('ID')
  
  chr_to_str = dict([(float(i), 'chr{0:d}'.format(i)) for i in xrange(0,23)])
  chr_to_str['X'] = 'chrX'
  chr_to_str['Y'] = 'chrY'

  alg['CHR'] = alg['CHR'].map(chr_to_str.get)

  alg = alg.ix[alg['EVENT_SIZE']>=1000,:]

  pm_to_m = pd.read_csv(pm_to_m_map_fpath, sep='\t', header=None, index_col=0, names=['ID','MERGED_ID'])
  alg_m = pd.concat((alg, pm_to_m), axis=1).dropna(axis=0, how='any', subset=['MERGED_ID','SAMPLES'])

  samples = sorted(list(set(chain(*[k.split(',') for k in alg_m.groupby('SAMPLES').groups.keys()]))))
   
  #x,y = samples[1:]
  #samples[1:] = y,x

  for s in samples:
    alg_m[s]=alg_m['SAMPLES'].map(lambda x:x.count(s))
  
  # Check if each center and approach made calls for each sample
  g = alg_m.groupby(['CENTER','TYPE_OF_COMPUTATIONAL_APPROACH'])
  trio_approaches = g.apply(lambda x:(x[samples].sum()>0).all())
 
  trio_index = None
  for trio_approach in trio_approaches.index[trio_approaches]:
    if trio_index == None:
      trio_index = g.get_group(trio_approach).index
    else:
      trio_index = trio_index.union(g.get_group(trio_approach).index)
  #print trio_approaches 
  #print alg_m.ix[trio_index,:].groupby(['CENTER','TYPE_OF_COMPUTATIONAL_APPROACH']).groups.keys()
  alg_m = alg_m.ix[trio_index,:]

  #print alg_m.ix[:10,['SAMPLES',]+samples]
  #print trio_approaches
  off,pat,mat = samples  
  if approach=='max':
    mids_s = get_mids_with_max_sample_labels(alg_m, samples)
  elif approach=='any':
    mids_s = get_mids_with_any_sample_labels(alg_m, samples)
  #print mids_s.ix[:30,:]
  #print mids_s.ix[(mids_s[off]+mids_s[pat]+mids_s[mat])==3,:]
  #print m_sites.ix[:10,:]
  #print mids_with_samples.shape
  #print mids_s.apply(lambda x:'{0}{1}{2}'.format(*x), axis=1).value_counts()
  if approach=='sep':
    g = alg_m.groupby(['CENTER','TYPE_OF_COMPUTATIONAL_APPROACH'])
    for p, alg_m_s in g:
      key = '{0}-{1}'.format(*p)
      mids_s = get_mids_with_any_sample_labels(alg_m_s, samples)
      final = infer_transmit(m_sites, mids_s, off,pat,mat,somatic_flag=False, sort_flag=False) 
      #m_sites[key]='.'
      new_idx = final.index.intersection(m_sites.index)
      #final['class'] = final['class'].astype(str)
      f_c =  pd.DataFrame(final.ix[new_idx,'class'].values, index=new_idx, columns=[key,])
      m_sites = pd.concat([m_sites, f_c], axis=1)
      m_sites[key] = m_sites[key].fillna('.')
    m_sites = sort_chr(m_sites)
    m_sites.to_csv(out_bed_fpath, header=True, index=False, sep='\t')
  elif hets_only:
    #print mids_s
    #print "="*10
    sort_chr(m_sites.ix[(mids_s[off]>=1).index,:]).to_csv(sys.stdout, header=False, index=False, sep='\t')
  elif all_dels:
    #sort_chr(m_sites.ix[(mids_s[off]>=1).index,:]).to_csv(sys.stdout, header=False, index=False, sep='\t')
    pd.Series((mids_s[off]>=1).index).to_csv(sys.stdout, header=False, index=False, sep='\t')
  else:
    final = infer_transmit(m_sites, mids_s, off,pat,mat,somatic_flag=somatic_flag) 
    final.to_csv(out_bed_fpath, header=False, index=False, sep='\t')

def infer_transmit(m_sites, mids_s, off,pat,mat, somatic_flag=False, sort_flag=True): 
  hom_idx = mids_s.index[(mids_s[off]+mids_s[pat]+mids_s[mat])==3]
  hom_v = m_sites.ix[hom_idx,:]
  hom_v['class'] = 'homA'

  snp = SNPTransmit(pat,mat)
  pat_idx = mids_s.index[na.logical_and((mids_s[off]+mids_s[pat])==2, mids_s[mat]==0)]
  pat_v = m_sites.ix[pat_idx,:]
  pat_v = snp.assign_dels_to_homolog(pat_v, 'pat')

  mat_idx = mids_s.index[na.logical_and((mids_s[off]+mids_s[mat])==2, mids_s[pat]==0)]
  mat_v = m_sites.ix[mat_idx,:]
  mat_v = snp.assign_dels_to_homolog(mat_v, 'mat')

  final = pd.concat([hom_v,pat_v,mat_v])
  final = final.ix[final['class']!='.',:]
  if sort_flag:
    final = sort_chr(final)
  if somatic_flag:
    final = m_sites.ix[na.logical_and(mids_s[off]==1,(mids_s[mat]+mids_s[pat])==0),:]    
  return final

if __name__ == '__main__':
  #algorithm_to_dels(approach='max')
  algorithm_to_dels(all_dels=True, approach='any')
  #algorithm_to_dels(approach='sep')
