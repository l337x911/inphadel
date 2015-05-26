from svphase.parser import ReadsParserDat
import sys

def print_reads(dat_fpath, a,b):
  parser = ReadsParserDat()
  for i in parser.get_reads(dat_fpath, '', a,b):
    print i 

if __name__ == '__main__':
  dat_fpath, a,b = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
  print_reads(dat_fpath, a,b)


