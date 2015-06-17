#!/bin/python

import struct, sys

s = struct.Struct('<LL?L?')

out = sys.stdout
if len(sys.argv)>1:
	print "outputting to %s"%sys.argv[1]
	out = open(sys.argv[1], 'wb')
	
for line in sys.stdin:
	line = line.strip().split('\t')
	out.write(s.pack(int(line[0]),int(line[1]), line[2]=='1',int(line[3]),line[4]=='1'))

if len(sys.argv)>1:
	out.close()
