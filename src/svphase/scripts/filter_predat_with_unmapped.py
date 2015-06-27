"""Mapper of liftover files to original predat files
Filters predat fmt files (see to_dat.sh) with an ordered liftover bed and liftover unmapped file list.

"""


import sys
from svphase.utils.common import iter_paired_lines


def predat_filter_and_convert(predat, liftover, contig_unmapped):
	unmapped = set()
	with open(contig_unmapped, 'rb') as f:
		for line in f:
			if line.startswith('#'): continue
			pos = line.strip().split('\t')[1]
			unmapped.add(pos)

	with open(predat, 'rb') as predat_f, open(liftover, 'rb') as liftover_f:
		
		for p1,p2 in iter_paired_lines(predat_f):
			row1 = p1.split('\t')
			row2 = p2.split('\t')
			assert row1[0]==row2[0]
			u1 = row1[1] in unmapped
			u2 = row2[1] in unmapped

			# Do not lift over reads mapping to position 0
			n1 = row1[1] == '0'
			n2 = row2[1] == '0'

			write_flag = False
			if (not u1) and (not u2):
				q1 = liftover_f.next().split('\t')
				q2 = liftover_f.next().split('\t')
				if not n1:
					row1[1] = q1[1]
					row2[3] = q1[1]
				if not n2:
					row1[3] = q2[1]
					row2[1] = q2[1]
				write_flag = True
			elif u1 and u2:
				# skip both bad reads
				pass
			else:
				# maintain pairing of read if only one coord was properly lifted
				q = liftover_f.next().split('\t')

				# In case a 0 gets unmapped
				if u1 and n1:
					row1[3] = q[1]
					row2[1] = q[1] 
					write_flag = True
				elif u2 and n2:
					row1[1] = q[1]
					row2[3] = q[1]
					write_flag = True
			if write_flag:	
				sys.stdout.write('\t'.join(row1))
				sys.stdout.write('\t'.join(row2))
if __name__ =='__main__':
	import argparse
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('predat', default=str, help='.predat file')
	parser.add_argument('liftover', default=str, help='bed generated from predat file')
	parser.add_argument('unmapped', default=str, help='unmapped file generated from predat file')

	args = parser.parse_args()
	
	assert args.predat.endswith('.predat')
	predat_filter_and_convert(args.predat, args.liftover, args.unmapped)

