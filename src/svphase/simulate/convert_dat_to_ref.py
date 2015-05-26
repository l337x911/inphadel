from svphase.chromosome import _get_contig_from_fpath, Chromosome, DelModifier 
from svphase.utils import reference
from svphase.parser import ReadsParserDat, ReadsWriterDat
import sys

def convert_to_ref(ref_fpath, sv_fpath, dat_fpath):
  
  contig = _get_contig_from_fpath(sv_fpath)
  #ref_obj = reference.Reference('/media/T01/data/hg/hg18/chr20.fa')  
  ref_obj = reference.Reference(ref_fpath.format(contig))  

  chrom = Chromosome(ref_obj, contig)
  chrom.get_seq()
  derchrom = DelModifier(ref_obj, contig, chrom, sv_fpath)
 
  dat_ref_fpath = "%s.ref.dat"%(dat_fpath[:-4])
  w = ReadsWriterDat(dat_ref_fpath)
  parser = ReadsParserDat()
  for row in parser.get_reads(dat_fpath, "del-%s"%contig):
    idx, posa,rev_ca,posb,rev_cb = row
    ref_posa, ref_posb = derchrom.derive_reference_position(posa), derchrom.derive_reference_position(posb)
    #print ref_posa, ref_posb
    if ref_posa is None or ref_posb is None: continue
    w.write(ref_posa, parser.revc_to_strd[rev_ca], ref_posb, parser.revc_to_strd[rev_cb])
  w.close() 

if __name__ == '__main__':
  #ref_fpath = '/home/adp002/data/hg/grch38/chroms/{0}.fa'
  ref_fpath = '/media/T02/data/hg/grch38/chroms/{0}.fa'
  svfpath = sys.argv[1]
  dat_fpath = sys.argv[2]
  convert_to_ref(ref_fpath, svfpath, dat_fpath) 
