import numpy as na
from svphase.parser import ReadsParserDat
import time


dist = na.array([0,1000,3000,5000,10000,25000,50000,1000000,10000000,60000000,100000000])
resolute = 1
contig = 'chr20'

def get_read_counts(dat_fpath, bed_fpath):
  parser = ReadsParserDat()

  intervals = []
  with open(bed_fpath,'rb') as f:
    for line in f:
      tokens = line.strip().split('\t')
      intervals.append((int(tokens[1])*resolute, int(tokens[2])*resolute+resolute-1))

  read_counts = na.zeros((len(intervals), dist.size), dtype=na.float)
  
  prev_time = time.clock()
  for i_idx, (start,end) in enumerate(intervals):
    #if i_idx>2500: break
    inserts = [abs(row[3]-row[1]) for row in parser.get_reads(dat_fpath, contig, start, end)]

    #bbins=na.bincount(na.searchsorted(dist,inserts))
    #read_counts[i_idx, :len(bbins)] += bbins

    for i in na.searchsorted(dist,inserts):
      read_counts[i_idx,i] += 1000./(end-start)

    if (i_idx+1)%100 == 0:
      curr_time = time.clock()
      print "%06d/%06d - %0.3f"%(i_idx+1, len(intervals), curr_time-prev_time)
      prev_time = curr_time
  return read_counts

def print_read_counts(c):

  format_dstr = "{: >12d}"
  format_str = "{: >12.2f}"
  print ' '.join([format_dstr.format(i) for i in dist])
  print "Mean:"
  print ' '.join([format_str.format(i) for i in na.mean(c, axis=0)])
  print "Variance"
  print ' '.join([format_str.format(i) for i in na.var(c, axis=0)])

def check_dat(dat_chr_fpath, bed_chr_fpath):
  assert bed_chr_fpath.endswith('.bed')
  load_fpath =  "%s_comp_stat.npz"%(dat_chr_fpath[:-4]) 
  force_new_flag = True

  if not force_new_flag and os.path.isfile(load_fpath):
    npzfile = na.load(load_fpath)
    original = npzfile['read_count']
  else:
    original = get_read_counts(dat_chr_fpath, bed_chr_fpath)  
    na.savez(load_fpath, read_count=original)
  #print original 
  print_read_counts(original)

def read_from_file(npzfpath):
  npzfile = na.load(npzfpath)
  print_read_counts(npzfile['read_count'])
 
if __name__=='__main__':
  import os, sys
  if len(sys.argv)>2:
    check_dat(sys.argv[1], sys.argv[2])
  else:
    read_from_file(sys.argv[1])

