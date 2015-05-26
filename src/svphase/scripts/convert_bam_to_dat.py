import sys
from svphase.parser import ReadsParserSAM, write_hic_to_dat

if __name__ == '__main__':
  bam_fpath = sys.argv[1]
  dat_fpath = sys.argv[2]
  contig = sys.argv[3]
  parser = ReadsParserSAM()
  write_hic_to_dat(parser, bam_fpath, dat_fpath, contig)
