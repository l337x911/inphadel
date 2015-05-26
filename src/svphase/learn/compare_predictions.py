import sys
from svphase.utils.config import BED_BIN 
import numpy as na
import pandas as pd
from subprocess import Popen, PIPE
from itertools import izip

pd.set_printoptions(max_rows=1000, max_columns=1000)

class VariantDB(object):
  def __init__(self, csv_fpath):
    self.fpath = csv_fpath
    self.df = pd.read_csv(csv_fpath, sep='\t', header=0)
  def get_deletions(self):
    pass


def reindex_bed_multi(df):
  m = pd.DataFrame(map(eval, df.index), columns=['contig','a','b'])
  df2 = pd.concat((m,df.reset_index(drop=True)), axis=1).set_index(['contig','a','b'])
  return df2

def get_bed_multi_from_str(bed_str):
  entries = [(i[0], int(i[1]), int(i[2])) for i in [t.split('\t') for t in bed_str.split('\n')[:-1]]]
  r = pd.MultiIndex.from_tuples(entries, names=['contig','a','b']) 
  return r

def trios_stat(score_df, pbed,mbed, overlap_fracs=[0.5,]):
  obed_num = "\n".join(["{0}\t{1}\t{2}\t{3}".format(i[0],i[1],i[2],c) for c,i in enumerate(score_df.index)])
  print '{0}/mergeBed -i'.format(BED_BIN)
  p = Popen('{0}/mergeBed -nms -i'.format(BED_BIN), shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
  merged_entries, merged_entries_err = p.communicate(obed_num)
  if merged_entries_err != '':
    print merged_entries_err
    raise

  entries = [t.split('\t')[-1].split(';') for t in merged_entries.split('\n')[:-1]]
  entries_dict = {}
  for c,vs in enumerate(entries):
    for v in vs:
      entries_dict[int(v)] = c
  score_idf = score_df.reset_index()
  overlapping_g = score_idf[['inc','pA','pB','homA']].max(axis=1).groupby(entries_dict)
  agg_idx = sorted([overlapping_g.groups[g_i][i] for g_i,i in overlapping_g.aggregate(na.argmax).iteritems()])
  
  merged_obed_str = '\n'.join(score_idf.ix[agg_idx, ['contig','a','b']].apply(lambda x:"{0}\t{1}\t{2}".format(*x), axis=1)) 
  merged_score_df = score_idf.ix[agg_idx,:]
  merged_score_df = merged_score_df.set_index(['contig','a','b'])
  merged_class_df = pd.DataFrame(merged_score_df.idxmax(axis=1), columns=['class',])

  for frac in overlap_fracs:
    m_str, p_str = '{0:0.2f}m'.format(frac), '{0:0.2f}p'.format(frac)
 
    merged_class_df[m_str] = 0
    pm = Popen('{0}/intersectBed -wa -f {2:0.4f} -a stdin -b {1}'.format(BED_BIN, mbed, frac), shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    pm_bed_str, pm_bed_str_err = pm.communicate(merged_obed_str)
    if pm_bed_str_err != '':
      print pm_bed_str_err
      raise
    elif pm_bed_str != '':
      merged_class_df.ix[get_bed_multi_from_str(pm_bed_str), m_str] = 1

    merged_class_df[p_str] = 0
    pp = Popen('{0}/intersectBed -wa -f {2:0.4f} -a stdin -b {1}'.format(BED_BIN, pbed, frac), shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    pp_bed_str, pp_bed_str_err = pp.communicate(merged_obed_str)

    if pp_bed_str_err != '':
      print pp_bed_str_err
      raise
    elif pp_bed_str != '':
      merged_class_df.ix[get_bed_multi_from_str(pp_bed_str), p_str] = 1    
 
  print merged_class_df


if __name__ == '__main__':
  score_df = pd.read_csv(sys.argv[1], header=0, index_col=0)
  score_df = reindex_bed_multi(score_df)
  paternal_bed, maternal_bed = sys.argv[2:4]
  trios_stat(score_df, paternal_bed, maternal_bed, overlap_fracs=[0.4, 0.5, 0.6, 0.7, 0.8, 0.9])

