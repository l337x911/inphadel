import sys
from svphase.parser import ReadsParserDat

def number_of_reads(fpath):
  r = ReadsParserDat()
  return r.get_single_read_count(fpath)   

def number_of_reads_at_dist(fpath, x):
  r = ReadsParserDat()
  n = 0
  for row in r.get_reads(fpath, ''):
    if abs(int(row[1])-int(row[3]))>x:
      n += 1
  return n

if __name__ == '__main__':
  assert len(sys.argv)>1
  dat_fpaths = sys.argv[1:]

  for dat_fpath in dat_fpaths:
    print dat_fpath, number_of_reads(dat_fpath)
    #print dat_fpath, number_of_reads_at_dist(dat_fpath, 25000)
