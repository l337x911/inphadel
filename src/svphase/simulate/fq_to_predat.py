import sys
from svphase import parser
from itertools import count
from svphase.utils.config import READLENGTH
import logging

logger = logging.getLogger('inphadel')

logger.setLevel(logging.DEBUG)
# Log to Console
ch = logging.StreamHandler()
#ch.setLevel(logging.INFO)
ch.setLevel(logging.DEBUG)

# Specifies format of log
ch.setFormatter(logging.Formatter('%(asctime)s [%(levelname) 9s] - %(message)s'))
logger.addHandler(ch)

def convert(fq, index_start):
	fmt = "{0:d}\t{1:d}\t{2:d}\t{3:d}\t{4:d}\n"	
	pecount = count(index_start)
	for line in fq:
		if not line.startswith('@'): continue
		#E.g. header "@del-del-chr20_14470409_14470948_2:0:0_2:0:0_b/1" 
		tokens = line[1:].strip().split(':')

		# guessing that is num_of_error
		contig, tstart, tend, num_of_errors = tokens[0].split('_')
		tstart,tend = int(tstart), int(tend)
		assert tstart<tend
		c = pecount.next()
		sys.stdout.write(fmt.format(c,tstart,0,tend-READLENGTH,1))
		sys.stdout.write(fmt.format(c,tend-READLENGTH,1,tstart,0))
			
		if ((c+1)%1000000)==0:
			logger.info("Finished %09d/%09d"%(c,0))

if __name__=='__main__':
	import argparse
	logger.info('Converting FQ to Predat')
	
	parser = argparse.ArgumentParser(description='Convert FQ to Predat format')
	parser.add_argument('index_start', type=int, help='Index start')

	args = parser.parse_args()
	convert(sys.stdin, args.index_start)
	
