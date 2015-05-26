from svphase import parser
from itertools import count
from svphase.utils.config import READLENGTH

def convert(fq_fpath):
  writer = parser.ReadsWriterDat("%s.dat"%(fq_fpath[:-5]))
  
  with open(fq_fpath, 'rb') as f:
    pecount = count()
    for line in f:
      if not line.startswith('@'): continue
      # E.g. header "@del-del-chr20_14470409_14470948_2:0:0_2:0:0_b/1" 
      tokens = line[1:].strip().split(':')

      # guessing that is num_of_error
      contig, tstart, tend, num_of_errors = tokens[0].split('_')
      tstart,tend = int(tstart), int(tend)
      assert tstart<tend
      writer.write(tstart,'+',tend-READLENGTH,'-')
      
      c = pecount.next()
      if (c+1)%1000000==0:
        print "Finished %09d/%09d"%(c,0)

      #if c==1000000: break
  writer.close()

if __name__=='__main__':
  import sys
  convert(sys.argv[1])
  
