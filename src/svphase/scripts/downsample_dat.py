from svphase.parser import ReadsParserDat, ReadsWriterDat
from svphase.simulate.hic import rand_gen
import sys

def downsample(dat_fpath, out_dat_fpath, ds=2):
  assert ds>1
  parser = ReadsParserDat()
  writer = ReadsWriterDat(out_dat_fpath)

  binary_rand = rand_gen(0,ds)
  for row in parser.get_reads(dat_fpath, contig=''):
    idx,posa,revca,posb,revcb = row
    if binary_rand.next()!=1: continue

    writer.write(posa, parser.revc_to_strd[revca], posb, parser.revc_to_strd[revcb])

  writer.close()

if __name__=='__main__':
  dat_fpath = sys.argv[1]
  ds = int(sys.argv[2])

  downsample(dat_fpath, "%s.down%d.dat"%(dat_fpath[:-4], ds), ds)
