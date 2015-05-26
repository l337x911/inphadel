import pandas as pd
import sys

def filter_bed(mid_fpath, pred_fpath):
  with open(mid_fpath, 'rb') as f:
    mids = set([line.strip() for line in f])
  
  df = pd.read_csv(pred_fpath, sep='\t', header=0)
  #print df['truth']
  df[df['truth'].isin(mids)].to_csv(sys.stdout, sep='\t', header=True, index=False)


if __name__ == '__main__':
  filter_bed(sys.argv[1], sys.argv[2])
