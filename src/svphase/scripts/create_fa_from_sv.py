from svphase.chromosome import Chromosome, DelModifier, _get_contig_from_fpath
from svphase.utils import reference
import sys

def fa_del_chromosome(ref_fpath, sv_fpath):
  contig = _get_contig_from_fpath(sv_fpath)
  ref_obj = reference.Reference(ref_fpath.format(contig))  

  chrom = Chromosome(ref_obj, contig)
  chrom.get_seq()
  derchrom = DelModifier(ref_obj, contig, chrom, sv_fpath)
  derchrom.out_fasta("%s.fa"%(sv_fpath[:-4]))

if __name__ == '__main__':
  ref_fpath = "/media/T02/data/hg/grch38/chroms/{0}.fa"
  svfpath = sys.argv[1]
  fa_del_chromosome(ref_fpath, svfpath) 
