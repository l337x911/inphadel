""" Only Filters VCF file to row-wise match generated liftover files """

import sys

# Read in VCF from stdin
contig_unmapped = sys.argv[1]

unmapped = set()
with open(contig_unmapped, 'rb') as f:
	for line in f:
		if line.startswith('#'): continue
		unmapped.add(tuple(line.strip().split('\t')[:2]))

for line in sys.stdin:
	if not tuple(line.split('\t')[:2]) in unmapped:
		sys.stdout.write(line)

