'''
Created on Jan 29, 2014

@author: anand
'''
import numpy as na
import re

from svphase.utils.config import CONTIG_POS_TYPE, CONTIG_TYPE_MAX, READLENGTH, REFERENCE_HG18, HIND3_STR

class Chromosome(object):
	def __init__(self, ref_obj, contig):
		self.contig = contig
		# 0-indexed
		self.cutsites = []
		self.cutsites_w = []
		self.het_snp_pos = []
		self.het_snp_info = []
		self.length = 0
		self.ref_obj = ref_obj
		self.seq = None
		self.re_pat = None

	def get_mappable(self):
		if self.seq is None:
			self.get_seq()
		s = na.fromstring(self.seq, dtype='S1')
		return s!='N'

	def get_seq(self):
		if self.seq is None:
			self.seq = self.ref_obj.get_contig(self.contig).upper()
			self.length = len(self.seq)
		return self.seq
	def modify(self, modifier, het_flag=False):
		raise
		#self.load_cutsites(modifier.origin_chr.re_pat)

		if not het_flag: return
		# Possibly slow and useless
		self.het_snp_pos = [na.where(modifier.map_to_origin==i)[0] for i in modifier.origin_chr.het_snp_pos]
		self.het_snp_pos = [i[0] for i in self.het_snp_pos if len(i)>0]
	def load_cutsites_from_bed(self, bed_fpath):
		if self.seq is None:
			self.get_seq()

		with open(bed_fpath, 'rb') as f:
			for line in f:
				tokens = line.strip().split('\t')
				self.cutsites_w.append((int(tokens[1]), int(tokens[2])))
		#self.cutsites_w = na.array(self.cutsites_w, dtype=CONTIG_POS_TYPE)
	def load_cutsites(self, re_pat):
		self.re_pat = re_pat

		if self.seq is None:
			self.get_seq()

		self.cutsites = na.array([m.start()+1 for m in re.finditer(re_pat, self.seq)], dtype=CONTIG_POS_TYPE)
	def load_snps(self, vcf_fpath):
		with open(vcf_fpath, 'rb') as f:
			for line in f:
				tokens = line.strip().split('\t')
				if tokens[0] != self.contig: continue
				self.het_snp_pos.append(int(tokens[1])-1)
		self.het_snp_pos = na.array(self.het_snp_pos, dtype=CONTIG_POS_TYPE)

	def out_fasta(self, fpath):
		with open(fpath, 'wb') as f:
			print >> f, ">%s-%s\n%s" % (self.origin_chr.contig, self.seq)

class SVModifier(Chromosome):
	def __init__(self, ref_obj, contig, chr_obj, sv_fpath):
		Chromosome.__init__(self, ref_obj, contig)
		self.sv_name = None
		# Note CONTIG_MAX_TYPE-1 are positions overlapping breakpoints
		self.map_origin_to_new = None
		self.map_new_to_origin = None
		self.map_new_to_origin_no_split = None
		self.origin_chr = chr_obj
		self.sv = []
		self._load(sv_fpath)
	#def get_mappable(self):
	#	assert not self.seq is None
	#	return Chromosome.get_mappable()
	def get_seq(self):
		assert not self.seq is None
		return self.seq
	def _load(self, sv_fpath):
		pass

class DelModifier(SVModifier):
	def __init__(self, ref_obj, contig, chr_obj, sv_fpath):
		SVModifier.__init__(self, ref_obj, contig, chr_obj, sv_fpath)
		self.sv_name = "del"
		self.contig = self.contig
	def _load(self, sv_fpath):
		self.map_origin_to_new = na.arange(self.origin_chr.length, dtype=CONTIG_POS_TYPE)

		with open(sv_fpath, 'rb') as f:
			for line in f:
				if line.startswith('track name'): continue

				tokens = line.strip().split('\t')
				p1, p2 = int(tokens[1]), int(tokens[2])
				self.sv.append((p1, p2))
				self.map_origin_to_new[p1:p2] = CONTIG_TYPE_MAX
		
		sv_idx = na.nonzero(self.map_origin_to_new < CONTIG_TYPE_MAX)

		self.map_new_to_origin = self.map_origin_to_new[sv_idx]
		
		self.map_origin_to_new[sv_idx] = na.arange(self.map_new_to_origin.size, dtype=CONTIG_POS_TYPE)

	def set_seq(self):
		origin_seq = na.fromstring(self.origin_chr.get_seq(), dtype='S1')
		self.seq = ''.join(origin_seq[self.map_new_to_origin])
		self.length = len(self.seq)

	def out_fasta(self, fpath):
		with open(fpath, 'wb') as f:
			print >> f, ">%s\n%s" % (self.contig, self.seq)

class DelModifierRefMap(DelModifier):
	def __init__(self, ref_obj, contig, chr_obj, sv_fpath):
		DelModifier.__init__(self, ref_obj, contig, chr_obj, sv_fpath)
		self.contig = self.contig
	def _load(self, sv_fpath):
		DelModifier._load(self,sv_fpath)

		# Since it's useless.
		self.seq = None
		self.origin_chr.seq = None

		self.map_new_to_origin_no_split = na.copy(self.map_new_to_origin)		
		# Remove split reads
		for i,j in self.sv:
			self.map_new_to_origin_no_split[i - READLENGTH-1:i] = CONTIG_TYPE_MAX

	def derive_reference_position(self, pos):
		if pos<0 or pos>=self.map_new_to_origin_no_split.size:
			return None

		p = self.map_new_to_origin_no_split[pos]
		if p==CONTIG_TYPE_MAX:
			return None
		return p
	
def print_cutsites():
	from svphase.utils import reference
	chrom = Chromosome(reference.Reference('/media/ET01/data/hg/hg18/chr19.fa'),'chr19')
	
	hind3 = re.compile('AAGCTT')
	chrom.load_cutsites(hind3)
	print "\n".join(["%s\t%d\t%d"%(chrom.contig,i,i+1) for i in chrom.cutsites])

def _get_contig_from_fpath(sv_fpath):
	#import os
	tokens = sv_fpath.split('.')
	contig = tokens[1]
	print contig
	return contig

def print_sv_cutsites(sv_fpath):
	from svphase.utils import reference
	contig = _get_contig_from_fpath(sv_fpath)
	ref_obj = reference.Reference(REFERENCE_HG18%(contig))	
	#sv_fpath = '/home/anand/Projects/assembly/data/gm12878/del.chr20.bed'

	chrom = Chromosome(ref_obj, contig)
	chrom.get_seq()
	derchrom = DelModifier(ref_obj, contig, chrom, sv_fpath)

	hind3 = re.compile(HIND3_STR)
	derchrom.load_cutsites(hind3)
	with open(sv_fpath[:-4]+'.hg18.hind3.bed', 'wb') as f:
		print >>f, "\n".join(["%s\t%d\t%d"%(derchrom.contig,i,i+1) for i in derchrom.cutsites])

if __name__ == '__main__':
	print_cutsites()
