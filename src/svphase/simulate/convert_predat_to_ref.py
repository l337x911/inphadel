from svphase.chromosome import Chromosome, DelModifierRefMap 
from svphase.utils import reference
import sys
from svphase.utils.common import logger

def convert_to_ref(contig, ref_fpath, sv_fpath):
	ref_obj = reference.Reference(ref_fpath)	
	
	chrom = Chromosome(ref_obj, contig)
	chrom.get_seq()
	derchrom = DelModifierRefMap(ref_obj, contig, chrom, sv_fpath)

	for row in sys.stdin:
		idx, posa,rev_ca,posb,rev_cb = row.split('\t')
		posa,posb = int(posa), int(posb)
		ref_posa, ref_posb = derchrom.derive_reference_position(posa), derchrom.derive_reference_position(posb)
		#print ref_posa, ref_posb
		if ref_posa is None or ref_posb is None: continue
		sys.stdout.write("{0}\t{1:d}\t{2}\t{3:d}\t{4}".format(idx,ref_posa, rev_ca, ref_posb, rev_cb))

if __name__ == '__main__':
	import argparse

	parser = argparse.ArgumentParser(description='Converts read coordinates from a delmodified object to reference coordinates')
	parser.add_argument('contig', help='chromosome to delete position')	
	parser.add_argument('bed', help='Simple bed file containing deletions (coordinates on reference)')
	parser.add_argument('reference_fasta',  help='path to reference fasta (must match contig)')
	args = parser.parse_args()


	#ref_fpath = '/home/adp002/data/hg/grch38/chroms/{0}.fa'
	convert_to_ref(args.contig, args.reference_fasta, args.bed) 
