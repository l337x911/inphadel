import numpy as na
from svphase.parser import ReadsParserDat, ReadsWriterDat
from svphase.chromosome import Chromosome, DelModifier
from svphase.utils.config import READLENGTH

def filter_dat(dat_fpath, vcf_fpath):
  p = ReadsParserDat()
  
  het_pos = []
  with open(vcf_fpath, 'rb') as f:
    for line in f:
      tokens = line.strip().split('\t')
      try:
        pos = int(tokens[1])
      except:
        print line
        raise
      gt = tokens[9]
      if gt[1]=='|':
        het_pos.append(pos)

  het_pos.append(10000000000)
  het_pos = na.sort(het_pos)

  q = ReadsWriterDat('%s.filt.dat'%dat_fpath[:-4])
  for row in p.get_reads(dat_fpath, ''):
    idx,posa,rev_ca,posb,rev_cb = row
    
    snpa, snpb = het_pos[na.searchsorted(het_pos, (posa, posb), side='right')]
    snpa -= posa
    snpb -= posb
    if (snpa>0 and snpa<=READLENGTH) or (snpb>0 and snpb<=READLENGTH):
      q.write(posa,p.revc_to_strd[rev_ca],posb,p.revc_to_strd[rev_cb])
  q.close()


if __name__=='__main__':
  import sys

  dat_fpath = sys.argv[1]
  vcf_fpath = sys.argv[2]

  filter_dat(dat_fpath, vcf_fpath)


