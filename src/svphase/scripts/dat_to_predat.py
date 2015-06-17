#!/bin/python

import struct, sys

s = struct.Struct('<LL?L?')

out = sys.stdout
if len(sys.argv)>1:
	print "outputting to %s"%sys.argv[1]
	out = open(sys.argv[1], 'wb')
	
while True:
	line = sys.stdin.read(s.size)
	if len(line)==0:
		break
	elif len(line) != s.size:
		raise Exception('Line is not the right size, %d. Try again', len(line))
		
	out.write("%d\t%d\t%d\t%d\t%d\n"%s.unpack(line))

if len(sys.argv)>1:
	out.close()
