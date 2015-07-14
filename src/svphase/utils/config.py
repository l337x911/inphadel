'''
Created on Jan 29, 2014

@author: anand
'''
import numpy as np
READLENGTH=100
RAND_BUFFER_SIZE = 50000
CONTIG_POS_TYPE = np.uint32
CONTIG_TYPE_MAX = np.iinfo(CONTIG_POS_TYPE).max
ZONE_TYPE = np.uint8
REFERENCE_FASTA="/home/anand/data/hg/18/hg18.fa"
REFERENCE_HG38="/home/anand/data/hg/38/hg38.fa"
#REFERENCE_DIR="/home/adp002/data/hg/hg18/%s.fa"
BED_BIN = "/home/anand/software/BEDTools-Version-2.12.0/bin"
THREADS = 7
PRECISION=0.000001
HIND3_STR = 'AAGCTT'
RANDOM_STATE=0
DATA_PREFIX = '/home/anand/T02/data/hic/'
SV_DATA = '/home/anand/Projects/inphadel/data'
TEST_DATA_PREFIX = '/home/anand/Projects/inphadel/abbrev_test'

#def get_open_fds():
#  '''
#    return the number of open file descriptors for current process
#    .. warning: will only work on UNIX-like os-es.
#  '''
#  import subprocess
#  import os

#  pid = os.getpid()
#  procs = subprocess.check_output([ "lsof", '-w', '-Ff', "-p", str( pid ) ] )

#  nprocs = len(filter(lambda s: s and s[ 0 ] == 'f' and s[1: ].isdigit(), procs.split( '\n' ) ))
#  return nprocs
