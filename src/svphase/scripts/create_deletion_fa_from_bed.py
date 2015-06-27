""" Takes a reference fasta sequence and removes bed positions """

from svphase.chromosome import Chromosome, DelModifier
from svphase.inphadel import ClassLabel
from svphase.utils import reference
from svphase.utils.common import logger
import sys
import os
import tempfile 

def create_del_chromosomes(contig, ref_fpath, sv_fpath, out_dir):
	ref_obj = reference.Reference(ref_fpath)  

	classes = ClassLabel()

	chrom = Chromosome(ref_obj, contig)
	chrom.get_seq()
	with tempfile.NamedTemporaryFile() as contigA_f, open(sv_fpath, 'rb') as f:
		for line in f:
			if line.startswith('track name'): continue
			tokens = line.strip().split('\t')
			if tokens[0]==contig and classes.is_deletion_on_a(tokens[3]):
				contigA_f.write(line)
		contigA_f.flush()
		
		derchromA = DelModifier(ref_obj, contig, chrom, contigA_f.name, no_split_flag=False)
		derchromA.out_fasta("{out}/{contig}.chromA.fa".format(out=out_dir, contig=contig))
	derchromA = None

	with tempfile.NamedTemporaryFile() as contigB_f, open(sv_fpath, 'rb') as f:
		for line in f:
			if line.startswith('track name'): continue
			tokens = line.strip().split('\t')
			if tokens[0]==contig and classes.is_deletion_on_b(tokens[3]):
				contigB_f.write(line)
		contigB_f.flush()
		
		derchromB = DelModifier(ref_obj, contig, chrom, contigB_f.name, no_split_flag=False)
		derchromB.out_fasta("{out}/{contig}.chromB.fa".format(out=out_dir, contig=contig))

if __name__ == '__main__':
	import argparse
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('contig', help='chromosome to delete position')	
	parser.add_argument('bed', help='Simple bed file containing deletions (coordinates on reference)')
	parser.add_argument('reference_fasta',  help='path to reference fasta (must match contig)')
	parser.add_argument('outdir', help='directory to output contig.chromA.fa and contig.chromB.fa')

	args = parser.parse_args()
	
	create_del_chromosomes(args.contig, args.reference_fasta, args.bed, args.outdir) 
