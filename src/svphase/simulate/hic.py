import numpy as na
from itertools import izip, chain
import time
from collections import defaultdict

import heapq
import re

from svphase.utils import reference
from svphase.utils.common import logger
from svphase.utils.config import THREADS, RAND_BUFFER_SIZE, CONTIG_POS_TYPE, ZONE_TYPE, HIND3_STR
from svphase.parser import ReadsParserDat, ReadsWriterPredat
from svphase.chromosome import Chromosome, DelModifier

HIND3 = re.compile(HIND3_STR)

def rand_gen(low=0, high=100):
	while True:
		for i in na.random.randint(low, high=high, size=RAND_BUFFER_SIZE):
			yield i
def rand_gen_float():
	while True:
		for i in na.random.random_sample(size=RAND_BUFFER_SIZE):
			yield i

class Zoning(object):
	def __init__(self, n, k, cutsites):
		self.n = n
		self.k = k
		self.cutsites = cutsites

		# initialize diag is 2 for within k cutsites
		# and 0 otherwise
		self.init = na.zeros(n, dtype=ZONE_TYPE)
		for c in self.cutsites:
			self.init[c - k:c + k] = 1

	def get_zone(self, x, y):
		# # assumes sorted_coord is sorted by ascending y-x, where y>x
		# for (x,y) in sorted_coord:
		# 	yield x,y, self.init[x]+self.init[y]
		return self.init[x] + self.init[y]

class ZoningRand(Zoning):
	def __init__(self, n, k, cutsites, mappable):
		Zoning.__init__(self, n, k, cutsites)
		self.mappable = mappable
		self.offset = None
		self.diag = None
		self.rand_count = None
		self.rand_dict = None
		self.mappable_diag = None
		self.zone_find_failures = defaultdict(set)
		
	def _load_offset(self, offset):
		self.diag = na.zeros(self.n - offset, dtype=ZONE_TYPE)
		idx = na.arange(self.n - offset)
		self.diag += self.init[:self.n - offset]
		self.diag += self.init[offset:self.n]
		
		# create an unmappable diag
		self.mappable_diag = na.logical_and(self.mappable[:self.n - offset], self.mappable[offset:self.n])

		self.rand_count = {}
		self.rand_dict = {}
		self.rand_iters = [None] * 3
		for z in range(3):
			self.rand_dict[z] = idx[na.logical_and(self.diag == z, self.mappable_diag)]
			self.rand_count[z] = self.rand_dict[z].size
		self.offset = offset

	def get_rand_pos(self, offset, z):
		if offset >= self.n:
			self.zone_find_failures['over'].add(offset)
			return None
		elif offset != self.offset:
			self._load_offset(offset)
		try:
			assert self.rand_count[z] > 0
		except AssertionError:
			self.zone_find_failures[z].add(offset)
			return None

		if self.rand_iters[z] is None:
			self.rand_iters[z] = rand_gen(low=0, high=self.rand_count[z])
		j = self.rand_dict[z][self.rand_iters[z].next()]
		return (j, j + offset)

	def test(self):
		for offset in xrange(self.n - 1, 0, -1):
			print offset
			self._load_offset(offset)
			print "".join(["%d" % i for i in self.diag])
	 
class HiCReadShuffler(object):
	def __init__(self, chr_obj, new_chr_obj, cutsite_d=500, resolute=100, zoning_same_flag=False, sample_rate=1.0):
		self.chr_obj = chr_obj
		self.rand_chr_obj = new_chr_obj
		self.resolute = resolute
		self.rand_add = rand_gen(0, high=self.resolute)
		self.rand_float = rand_gen_float()
		self.sample_rate = sample_rate
		self.cutsite_k = cutsite_d / self.resolute
		
		self.chr_diag = None
		self.rand_chr_diag = None
		
		self._load(zoning_same_flag=zoning_same_flag)

	def _load(self, zoning_same_flag=False):
		self.load_rand()
		if zoning_same_flag:
			self.chr_diag = self.rand_chr_diag
		else:
			n = self.chr_obj.length / self.resolute
			cutsites = na.array(self.chr_obj.cutsites) / self.resolute
			self.chr_diag = Zoning(n, self.cutsite_k, cutsites)

	def load_rand(self, rand_chr_obj=None):
		if not rand_chr_obj is None:
			self.rand_chr_obj = rand_chr_obj
		n = self.rand_chr_obj.length / self.resolute
		mappable = self.rand_chr_obj.get_mappable()
		
		mappable_res = mappable[:n * self.resolute].reshape(n, self.resolute).sum(axis=1) == self.resolute
		assert len(mappable_res) == n 

		cutsites = na.array(self.rand_chr_obj.cutsites) / self.resolute
		self.rand_chr_diag = ZoningRand(n, self.cutsite_k, cutsites, mappable_res)

	def stat(self):
		# print "Read Shuffle Error Report:"
		if len(self.rand_chr_diag.zone_find_failures) > 0:
			for z, v in self.rand_chr_diag.zone_find_failures.iteritems():
				logger.debug("Z:{0}|{1:09d}".format(z, len(v)))

	def simple_shuffle_to(self, dat_fpath, out_predat_fpath):
		reader = ReadsParserDat()
	 
		reads = []
		for row in reader.get_reads(dat_fpath, ''):
			if self.rand_float.next() > self.sample_rate: continue

			posa = row[1]  
			posb = row[3]
			if posa > posb:
				posa, posb = posb, posa

			# se end reads have a pair that is 0 
			if posa==0 or posb==0:
				pos = max(posa,posb)
				posa = pos
				posb = pos
			posa = posa / self.resolute
			posb = posb / self.resolute

			reva, revb = reader.get_strandedness(row)
			heapq.heappush(reads, (posb - posa, posa, posb, reva, revb))

		logger.info('Starting to shuffle reads...')
		with ReadsWriterPredat(out_predat_fpath) as f:
			while True:
				try:
					dist, posa, posb, reva, revb = heapq.heappop(reads)
				except IndexError:
					break
				
				z = self.chr_diag.get_zone(posa, posb)
				r_pos = self.rand_chr_diag.get_rand_pos(posb - posa, z)

				if r_pos is None: continue

				# include single end (se) reads
				if posa==posb:
					assert r_pos[0]==r_pos[1]			
					f.write(r_pos[0] * self.resolute + self.rand_add.next(), reva, 0, revb)
				else:
					f.write(r_pos[0] * self.resolute + self.rand_add.next(), reva, r_pos[1] * self.resolute + self.rand_add.next(), revb)
				if len(reads)%1000000==0:
					logger.debug("Processed {0:010.0f} / {1:010.0f}".format((reader.n/2)*self.sample_rate-len(reads), self.sample_rate*reader.n / 2))

		if logger.getEffectiveLevel()==logging.DEBUG:
			logging.debug("Error Report;")
			self.stat()
		 
def test_zoning_large(res=10):
	import pandas as pd
	n = 100000000 / res
	k = 500 / res
	read_count = 10000
	cutsites_df = pd.read_csv('/home/anand/Projects/assembly/data/chr19.hg18.hind3.bed', sep='\t', header=None, names=['chr', 'start', 'end'])
	cutsites = cutsites_df['start'].astype(int) / res

	xs = na.random.randint(low=0, high=n, size=read_count)
	ys = na.array([na.random.randint(low=x + 100 / res, high=x + 500 / res) for x in xs])
	
	coords = sorted(izip(ys - xs, xs, ys))
	w = Zoning(n, k, cutsites)
	wr = ZoningRand(n, k, cutsites)

	offsets = set()
	prev = time.clock() 
	start_time = prev
	for idx, (offset, x, y) in enumerate(coords):
		if (len(offsets)) % 50 == 1 and not offset in offsets: 
			# if (idx+1)%500==1: 
			curr = time.clock()
			print "Processed: {0:d} reads @ {2:d} offsets in {1:03.2f} sec".format(idx, curr - prev, len(offsets))
			prev = curr
		try:
			z = w.get_zone(x, y)
		except:
			print "Get Zone error", x, y, offset, idx
			raise
		try:
			nx, ny = wr.get_rand_pos(offset, z)
		except:
			print "Random position error"
			raise
		offsets.add(offset)
	print "Total time: {0:03.2f} seconds".format(time.clock() - start_time)

def test_zoning():
	n = 50
	k = 3
	cutsites = [4, 10, 25, 30]
	w = Zoning(n, k, cutsites)
	a, b = na.meshgrid(na.arange(n), na.arange(n))
	coords = sorted([(x, y) for x, y in izip(a.flatten(), b.flatten()) if y >= x])
	d = {}
	zone_size = defaultdict(int)
	for x, y in coords:
		z = w.get_zone(x, y)
		d[(x, y)] = z
		zone_size[z] += 1
	# print d
	for y in xrange(n - 1, 0, -1):
		print "".join(["%d" % d[(x, y)] for x in xrange(0, y + 1)])
	wr = ZoningRand(n, k, cutsites)
	# wr.test()
	for offset in xrange(n - 1, 0, -1): 
		# print offset
		wr._load_offset(offset)
		zone_rand = "".join(["%d" % i for i in wr.diag])
		zone_reg = "".join(["%d" % d[(j, j + offset)] for j in xrange(0, n - offset)])
		assert zone_rand == zone_reg

	# randomize test
	freq_d = defaultdict(int)
	zone_bias = rand_gen(low=-1, high=4)
	for x, y in sorted(chain(*[coords, ] * 50), key=lambda old_coord:old_coord[1] - old_coord[0]):
		z = w.get_zone(x, y)
		if z < zone_bias.next(): continue

		nx, ny = wr.get_rand_pos(y - x, z)
		freq_d[(nx, ny)] += 1

	for y in xrange(n - 1, 0, -1):
		print " ".join(["%03d" % freq_d[(x, y)] for x in xrange(0, y + 1)])
	# zone_frequencies 
	zone_freq = defaultdict(int)
	for (x, y), v in freq_d.iteritems():
		zone_freq[w.get_zone(x, y)] += v
	
	print "Zone Frequencies"
	for z in [0, 1, 2]:
		print "{3}: {0:3d} {1:3d} {2:.2f}".format(zone_freq[z], zone_size[z], float(zone_freq[z]) / zone_size[z], z)

if __name__ == '__main__':
	import argparse
	import logging
	
	parser = argparse.ArgumentParser(description="Shuffles reads from a HiC file, and assuming Hind3 bias")
	
	parser.add_argument('contig', help='chromosome to delete position')	
	parser.add_argument('bed', help='Simple bed file containing deletions (coordinates on reference)')
	parser.add_argument('reference_fasta',  help='path to original reference fasta (must match contig)')
	parser.add_argument('reference_der_fasta',  help='path to derived reference fasta (must match contig)')
	parser.add_argument('dat', help='HiC reads in DAT format')
	parser.add_argument('out_predat', help='Output HiC reads in PREDAT format')
	parser.add_argument('-r,--resolute', dest='resolution', type=int, default=10,  help='Resolution fo read shuffling')
	parser.add_argument('-f,--frac-reads', dest='frac_reads', type=float, default=1.0,  help='Fraction of reads to shuffle')
	parser.add_argument('--debug', action='store_true', default=False,  help='Show debug messages')
	args = parser.parse_args()

	if args.debug:
		logger.setLevel(logging.DEBUG)
	else:
		logger.setLevel(logging.INFO)

	origin_ref_obj = reference.Reference(args.reference_fasta)
	new_ref_obj = reference.Reference(args.reference_der_fasta)

	chrom = Chromosome(origin_ref_obj, args.contig)
	chrom.load_cutsites(HIND3)

	derchrom = DelModifier(new_ref_obj, args.contig, chrom, args.bed)
	derchrom.set_seq()
	derchrom.load_cutsites(HIND3)

	shuffler = HiCReadShuffler(chrom, derchrom, resolute=args.resolution, sample_rate=args.frac_reads)
	logger.info('Loaded HiC structures')
	shuffler.simple_shuffle_to(args.dat, args.out_predat) 

	# test_zoning()
	# simulate_chromosome()
	#simulate_del_chromosome()
	# contig, hap, outnum = sys.argv[1], sys.argv[2], sys.argv[3]
	# simulate_del_chromosome_debruijn(contig, hap, int(outnum))
	# test_zoning_large()
