import sys
import random

def shuffle_bed_split(bed_fpath):
  with open(bed_fpath, 'rb') as f:
    bed = [line.strip() for line in f]
      
  idx = range(len(bed))
  random.shuffle(idx)
  n = len(bed)/2

  beda_fpath,bedb_fpath = "%s.A.bed"%bed_fpath[:-4], "%s.B.bed"%bed_fpath[:-4]
  beda = open(beda_fpath, 'wb')
  for i in sorted(idx[:n]):
    print >>beda, bed[i]
  beda.close()

  bedb = open(bedb_fpath, 'wb')
  for i in sorted(idx[n:]):
    print >>bedb, bed[i]
  bedb.close()

if __name__=='__main__':
  shuffle_bed_split(sys.argv[1])

