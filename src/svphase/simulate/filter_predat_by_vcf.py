import numpy as np
import sys
from svphase.chromosome import Chromosome, DelModifier
from svphase.utils.config import READLENGTH
from svphase.utils.common import iter_paired_lines

def filter_predat(vcf_fpath):
	p = ReadsParserDat()
	
	het_pos = []
	with open(vcf_fpath, 'rb') as f:
		for line in f:
			tokens = line.strip().split('\t')
			try:
				pos = int(tokens[1])
			except:
				print line
				raise
			gt = tokens[9]
			if gt[1] == '|':
				het_pos.append(pos)

	het_pos.append(10000000000)
	het_pos = np.sort(het_pos)

	for rowa,rowb in iter_paired_lines(sys.stdin):
		idx, posa, rev_ca, posb, rev_cb = rowa.strip().split('\t')
		
		snpa, snpb = het_pos[np.searchsorted(het_pos, (posa, posb), side='right')]
		snpa -= posa
		snpb -= posb
		if (snpa > 0 and snpa <= READLENGTH) or (snpb > 0 and snpb <= READLENGTH):
			sys.stdout.write(rowa)
			sys.stdout.write(rowb)


if __name__ == '__main__':
	import argparse

	parser = argparse.ArgumentParser(description='Filters reads from an unsorted predat file (stdin) and filters by overlapping het SNPs in vcf. (No sanity of chrom)')
	parser.add_argument('vcf', help='VCF 4.2 containing scaffold of phased variants (allele A | allele B)')
	args = parser.parse_args()

	filter_predat(args.vcf)

