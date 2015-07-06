#!/bin/python
import struct, sys
import shlex
import os
from subprocess import Popen, PIPE

s = struct.Struct('<LL?L?')

def main(out):
	env = dict(os.environ)
	env['LC_ALL'] = 'C'
	sort_cmd = 'sort -n -S 10G -k 2,2'
	p = Popen(shlex.split(sort_cmd), stdout=PIPE, stdin=sys.stdin, env=env, bufsize=1) 

	with p.stdout:
		for line in iter(p.stdout.readline, b''):
			line = line.strip().split('\t')
			out.write(s.pack(int(line[0]),int(line[1]), line[2]=='1',int(line[3]),line[4]=='1'))
	p.wait()

if __name__=='__main__':
	import argparse
	parser = argparse.ArgumentParser(description="Converts an unsorted PREDAT file to a DAT (sorted binary)")
	
	parser.add_argument('-o,--out', dest='out', default=None, help='path to output file')
	args = parser.parse_args()
	
	if args.out is None:
		out = sys.stdout
	else:
		out = args.out
	main(out)
	if not args.out is None:
		out.close()
	
