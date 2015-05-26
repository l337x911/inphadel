from svphase.parser import ReadsParserDat, ReadsWriterDat
from svphase.simulate.hic import rand_gen

def split(dat_fpath, out1_dat_fpath, out2_dat_fpath):
  parser = ReadsParserDat()
  writers = (ReadsWriterDat(out1_dat_fpath), ReadsWriterDat(out2_dat_fpath))
  binary_gen = rand_gen(0,2)
  
  for row in parser.get_reads(dat_fpath, contig=''):  
    idx,posa,revca,posb,revcb = row
    strda, strdb = parser.revc_to_strd[revca], parser.revc_to_strd[revcb]
    writers[binary_gen.next()].write(posa,strda,posb,strdb)

  writers[0].close() 
  writers[1].close()

if __name__=='__main__':
  import sys
  
  split(sys.argv[1], sys.argv[2], sys.argv[3])
