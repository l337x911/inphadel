'''
Created on Mar 24, 2011

@author: anand
'''

from collections import namedtuple
import sys

SAM_FLAGS_H = ['multi', 'proper_align',
      'query_unmapped', 'mate_unmapped',
      'rev_strand_of_query', 'rev_strand_of_next_mate',
      'first_frag', 'last_frag',
      'secondary_alignment', 'fail_quality_check'
      'duplication']

def get_flag(i,x):
  return bool((2**i) & x)
# pulls out flags from FLAG INT
SAM_FLAGS = dict([(flag, lambda j,x:get_flag(j,x)) for i, flag in enumerate(SAM_FLAGS_H)])

SAM_FLAGS_H = dict([(flag, i) for i, flag in enumerate(SAM_FLAGS_H)])

SAM_COLS = 'qname flags rname pos mapq cigar rnext pnext tlen seq qual'.split(' ')
SAMFragment = namedtuple('SAMFragment', SAM_COLS)

def filter_seqN(*fragments):
  for fragment in fragments:
    if fragment.seq.count('N')>0.4*len(fragment.seq):
      return True
  return False

def filter_unmapped(*fragments):
  for fragment in fragments:
    if (int(fragment.flags) & 0x4 ) == 0x4:
      return True
  return False

class SAMParser(object):
  def __init__(self):
    self.file = None
  def __iter__(self):
    if self.file is None:
      raise StopIteration
    return self
  def next(self):
    row = self.file.next().strip()
    while row[0] == '@':
      row = self.file.next().strip()
    
    return SAMFragment._make(row.split('\t')[:11])
    
  def parse(self, sam_fpath):
    if sam_fpath=="-":
      self.file = sys.stdin
    else:
      self.file = open(sam_fpath, 'rb')
    return iter(self)

  def parse_file(self, sam_file):
    self.file = sam_file
    return iter(self)  

  def parse_string(self, sam_read_str):    
    self.file = iter(sam_read_str.split('\n')[:-1])
    return iter(self)    

class SAMFiltering(object):
  def __init__(self, sam_fpath, filter_fpath):
    self.sam_fpath = sam_fpath
    self.filter_fpath = filter_fpath
  def filter(self, out=None, required=True):
    if out is None:
      out = sys.stdout
    
    with open(self.filter_fpath, 'rb') as f:
      header = f.next().strip()[1:].split('\t')
      cols_idx = []
      ignore_cols = []
      for idx, h in enumerate(header):
        try:
          cols_idx.append(SAM_COLS.index(h))
        except:
          ignore_cols.append(idx)
      
      frags_of_interest = set()
      for line in f:
        tokens = line.strip().split('\t')
        for idx in ignore_cols[::-1]:
          tokens.pop(idx)
        frags_of_interest.add(tuple(tokens))
        
    parser = SAMParser()
    if required is True:
      for frag in parser.parse(self.sam_fpath):
        key = tuple([frag[i] for i in cols_idx])
        if key in frags_of_interest:
          frags_of_interest.remove(key)
          print >>out, '\t'.join(frag)
    else:
      for frag in parser.parse(self.sam_fpath):
        key = tuple([frag[i] for i in cols_idx])
        if not key in frags_of_interest:
          #frags_of_interest.remove(key)
          print >>out, ('\t'.join(frag)).strip()
    
if __name__ == '__main__':
  #sam_fpath, filter_fpath = sys.argv[1], sys.argv[2]
  sam_fpath = '/home/anand/data/ambre/A01_1/aligned_blasr_h5_contigs_kmodes_min50_mod_filtered.sam'
  filter_fpath = '/home/anand/data/ambre/A01_1/aligned_blasr_h5_contigs_kmodes_min50_mod_filtered_q500.txt'
  
  filterer = SAMFiltering(sam_fpath, filter_fpath)
  
  with open('/home/anand/data/ambre/A01_1/aligned_blasr_h5_contigs_kmodes_min50_mod_filtered_q500.sam', 'wb') as f:
  #with open("%s_filtered.sam"%(sam_fpath[:-4]), 'wb') as f:
    filterer.filter(f, required=True)
  
