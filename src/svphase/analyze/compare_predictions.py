import pandas as pd
import sys, os

def print_prediction_counts(df):
  print df.name
  nan_count =  pd.isnull(df['pred']).sum()
  print "NanCount", nan_count 
  freq = df['pred'].value_counts()
  call_count = freq.sum()
  freq = freq.append(pd.Series(freq['pA']+freq['pB'], index=['het',]))
  print freq
  print "TotalCount:", nan_count+call_count, call_count/float(nan_count+call_count)

def remove_nans(s):
  return s[pd.notnull(s)]

def compare(bed1_fpath, bed2_fpath):
  df1 = pd.read_csv(bed1_fpath, sep='\t', header=0).set_index(['chr','start','end'])
  df2 = pd.read_csv(bed2_fpath, sep='\t', header=0).set_index(['chr','start','end'])
  df1.name=os.path.basename(bed1_fpath)
  df2.name=os.path.basename(bed2_fpath)

  print_prediction_counts(df1)
  print_prediction_counts(df2)

  print "ComparingPredictions"
  # remove nans first
  df1_p = remove_nans(df1['pred'])
  df2_p = remove_nans(df2['pred'])
  assert all(df1_p.index.values==df2_p.index.values)

  print (df1_p==df2_p).sum()/float(df1_p.size)

if __name__=='__main__':
  compare(sys.argv[1], sys.argv[2])
