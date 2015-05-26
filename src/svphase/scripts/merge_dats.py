from svphase.parser import ReadsParserDat, ReadsWriterDat
from itertools import chain

def merge(dat1_fpath, dat2_fpath, out_dat_fpath):
  parser = ReadsParserDat()
  writer = ReadsWriterDat(out_dat_fpath)
  
  for row in chain(parser.get_reads(dat1_fpath, contig=''), parser.get_reads(dat2_fpath, contig='')):
    idx,posa,revca,posb,revcb = row
    strda, strdb = parser.revc_to_strd[revca], parser.revc_to_strd[revcb]
    writer.write(posa,strda,posb,strdb)

  writer.close() 


if __name__=='__main__':
  import sys
  
  merge(sys.argv[1], sys.argv[2], sys.argv[3])
  pass
