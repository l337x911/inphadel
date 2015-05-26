import numpy as na
from itertools import izip, chain
import time
from collections import defaultdict

import heapq
import re

from svphase.utils import reference
from svphase.utils.config import THREADS, RAND_BUFFER_SIZE, CONTIG_POS_TYPE, ZONE_TYPE, HIND3_STR
from svphase.parser import ReadsParserDat, ReadsWriterDat
from svphase.chromosome import Chromosome, DelModifier

HIND3 = re.compile(HIND3_STR)

def rand_gen(low=0,high=100):
  while True:
    for i in na.random.randint(low,high=high, size=RAND_BUFFER_SIZE):
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
      self.init[c-k:c+k] = 1

  def get_zone(self, x,y):
    ## assumes sorted_coord is sorted by ascending y-x, where y>x
    #for (x,y) in sorted_coord:
    #  yield x,y, self.init[x]+self.init[y]
    return self.init[x]+self.init[y]

class ZoningRand(Zoning):
  def __init__(self, n, k, cutsites, mappable):
    Zoning.__init__(self,n,k,cutsites)
    self.mappable = mappable
    self.offset = None
    self.diag = None
    self.rand_count = None
    self.rand_dict = None
    self.mappable_diag = None
    self.zone_find_failures = defaultdict(set)
    
  def _load_offset(self, offset):
    self.diag = na.zeros(self.n-offset, dtype=ZONE_TYPE)
    idx = na.arange(self.n-offset)
    self.diag += self.init[:self.n-offset]
    self.diag += self.init[offset:self.n]
    
    # create an unmappable diag
    self.mappable_diag = na.logical_and(self.mappable[:self.n-offset],self.mappable[offset:self.n])

    self.rand_count = {}
    self.rand_dict = {}
    self.rand_iters = [None]*3
    for z in range(3):
      self.rand_dict[z] = idx[na.logical_and(self.diag==z, self.mappable_diag)]
      self.rand_count[z] = self.rand_dict[z].size
    self.offset = offset

  def get_rand_pos(self, offset, z):
    if offset>=self.n:
      self.zone_find_failures['over'].add(offset)
      return None
    elif offset!=self.offset:
      self._load_offset(offset)
    try:
      assert self.rand_count[z]>0
    except AssertionError:
      self.zone_find_failures[z].add(offset)
      return None

    if self.rand_iters[z] is None:
      self.rand_iters[z] = rand_gen(low=0, high=self.rand_count[z])
    j = self.rand_dict[z][self.rand_iters[z].next()]
    return (j,j+offset)

  def test(self):
    for offset in xrange(self.n-1, 0, -1):
      print offset
      self._load_offset(offset)
      print "".join(["%d"%i for i in self.diag])
   
class HiCReadShuffler(object):
  def __init__(self, chr_obj, new_chr_obj, cutsite_d=500, resolute=100, zoning_same_flag=False, debug_flag=False):
    self.chr_obj = chr_obj
    self.rand_chr_obj = new_chr_obj
    self.resolute = resolute
    self.rand_add = rand_gen(0,high=self.resolute)
    self.rand_binary = rand_gen(0,high=2)
    self.DEBUG_FLAG = debug_flag
    self.cutsite_k = cutsite_d/self.resolute
    
    self.chr_diag = None
    self.rand_chr_diag = None
    
    self._load(zoning_same_flag=zoning_same_flag)

  def _load(self, zoning_same_flag=False):
    self.load_rand()
    if zoning_same_flag:
      self.chr_diag = self.rand_chr_diag
    else:
      n = self.chr_obj.length/self.resolute
      cutsites = na.array(self.chr_obj.cutsites)/self.resolute
      self.chr_diag = Zoning(n, self.cutsite_k, cutsites)

  def load_rand(self, rand_chr_obj=None):
    if not rand_chr_obj is None:
      self.rand_chr_obj = rand_chr_obj
    n = self.rand_chr_obj.length/self.resolute
    mappable = self.rand_chr_obj.get_mappable()
    
    mappable_res = mappable[:n*self.resolute].reshape(n,self.resolute).sum(axis=1)==self.resolute
    assert len(mappable_res)==n 

    cutsites = na.array(self.rand_chr_obj.cutsites)/self.resolute
    self.rand_chr_diag = ZoningRand(n, self.cutsite_k, cutsites, mappable_res)

  def stat(self):
    #print "Read Shuffle Error Report:"
    if len(self.rand_chr_diag.zone_find_failures)>0:
      for z,v in self.rand_chr_diag.zone_find_failures.iteritems():
        print "Z:{0}|{1:09d}".format(z,len(v))

  def simple_shuffle_to(self, dat_fpath, out_dat_fpath, het_flag=False):
    writer = ReadsWriterDat(out_dat_fpath)
    reader = ReadsParserDat()
   
    reads = []
    for row in reader.get_reads(dat_fpath, ''):
      if het_flag and self.rand_binary.next()==1: continue

      posa, posb = row[1]/self.resolute, row[3]/self.resolute
      if posa>posb:
        posa,posb = posb,posa
      reva, revb = reader.get_strandedness(row)
      heapq.heappush(reads, (posb-posa, posa,posb, reva,revb))

    start_time = time.clock()
    prev_time = start_time
    while True:
      try:
        dist, posa, posb, reva, revb = heapq.heappop(reads)
      except IndexError:
        break

      z = self.chr_diag.get_zone(posa,posb)
      r_pos = self.rand_chr_diag.get_rand_pos(posb-posa,z)

      if r_pos is None: continue
      writer.write(r_pos[0]*self.resolute+self.rand_add.next(), reva, r_pos[1]*self.resolute+self.rand_add.next(), revb)
      if self.DEBUG_FLAG and (writer.idx+1) %100000 == 0:
        curr_time = time.clock()
        print "Wrote {0:010d} / {1:010d}: {2:04.2f} min".format(writer.idx+1,reader.n/2, (curr_time-prev_time)/60.0)
        #break
        prev_time = curr_time
    writer.close()
    if self.DEBUG_FLAG:
      print "Error Report;"
      self.stat()
      print "Total Time: {0:04.2f} min".format((time.clock()-start_time) / 60.0) 
     
def test_zoning_large(res=10):
  import pandas as pd
  n = 100000000/res
  k = 500/res
  read_count = 10000
  cutsites_df = pd.read_csv('/home/anand/Projects/assembly/data/chr19.hg18.hind3.bed', sep='\t', header=None, names=['chr','start','end'])
  cutsites = cutsites_df['start'].astype(int)/res

  xs = na.random.randint(low=0,high=n, size=read_count)
  ys = na.array([na.random.randint(low=x+100/res,high=x+500/res) for x in xs])
  
  coords = sorted(izip(ys-xs, xs, ys))
  w = Zoning(n,k,cutsites)
  wr = ZoningRand(n,k,cutsites)

  offsets = set()
  prev = time.clock() 
  start_time = prev
  for idx, (offset, x,y) in enumerate(coords):
    if (len(offsets))%50==1 and not offset in offsets: 
      #if (idx+1)%500==1: 
      curr = time.clock()
      print "Processed: {0:d} reads @ {2:d} offsets in {1:03.2f} sec".format(idx, curr-prev, len(offsets))
      prev = curr
    try:
      z = w.get_zone(x,y)
    except:
      print "Get Zone error", x,y, offset, idx
      raise
    try:
      nx,ny = wr.get_rand_pos(offset,z)
    except:
      print "Random position error"
      raise
    offsets.add(offset)
  print "Total time: {0:03.2f} seconds".format(time.clock()-start_time)

def test_zoning():
  n = 50
  k = 3
  cutsites = [4,10,25,30]
  w = Zoning(n,k,cutsites)
  a,b = na.meshgrid(na.arange(n), na.arange(n))
  coords = sorted([(x,y) for x,y in izip(a.flatten(), b.flatten()) if y>=x])
  d = {}
  zone_size = defaultdict(int)
  for x,y in coords:
    z = w.get_zone(x,y)
    d[(x,y)] = z
    zone_size[z] += 1
  #print d
  for y in xrange(n-1,0,-1):
    print "".join(["%d"%d[(x,y)] for x in xrange(0,y+1)])
  wr = ZoningRand(n,k,cutsites)
  #wr.test()
  for offset in xrange(n-1, 0,-1): 
    #print offset
    wr._load_offset(offset)
    zone_rand = "".join(["%d"%i for i in wr.diag])
    zone_reg = "".join(["%d"%d[(j,j+offset)] for j in xrange(0,n-offset)])
    assert zone_rand==zone_reg

  # randomize test
  freq_d = defaultdict(int)
  zone_bias = rand_gen(low=-1,high=4)
  for x,y in sorted(chain(*[coords,]*50), key=lambda old_coord:old_coord[1]-old_coord[0]):
    z = w.get_zone(x,y)
    if z<zone_bias.next(): continue

    nx,ny = wr.get_rand_pos(y-x,z)
    freq_d[(nx,ny)] += 1

  for y in xrange(n-1,0,-1):
    print " ".join(["%03d"%freq_d[(x,y)] for x in xrange(0,y+1)])
  # zone_frequencies 
  zone_freq = defaultdict(int)
  for (x,y), v in freq_d.iteritems():
    zone_freq[w.get_zone(x,y)] += v
  
  print "Zone Frequencies"
  for z in [0,1,2]:
    print "{3}: {0:3d} {1:3d} {2:.2f}".format(zone_freq[z],zone_size[z], float(zone_freq[z])/zone_size[z], z)

def simulate_chromosome():
  contig = 'chr20'
  t_n = 11
  hic_dat_fpath = '/media/T02/data/gm12878/hic/dat/%s.dat'%contig
  new_hic_dat_fpath = '/media/T02/data/gm12878/hic/dat/%s_%03d.dat'%(contig,t_n)
  #bed_fpath = '/home/anand/Projects/assembly/data/chr20.hg18.hind3.bed'
  chrom = Chromosome(reference.Reference('/media/ET01/data/hg/hg18/%s.fa'%contig),contig)
  chrom.get_seq()
  # HindIII regular expression
  chrom.load_cutsites(HIND3)
  #chrom.load_cutsites_from_bed(bed_fpath)

  shuffler = HiCReadShuffler(chrom, chrom, debug_flag=True)
  print "Loaded HiC Struct"
  shuffler.simple_shuffle_to(hic_dat_fpath, new_hic_dat_fpath) 

def simulate_del_chromosome():
  contig = 'chr20'
  origin_ref_obj = reference.Reference('/media/ET01/data/hg/hg18/%s.fa'%contig)  
  new_ref_obj = reference.Reference('/media/T02/data/hg/grch38/chroms/%s.fa'%contig)  
  sv_fpath = '/home/anand/Projects/assembly/data/simulations/sim.chr20.0.bed'
  t_n = 11
  resolute = 10
  hic_dat_fpath = '/media/T02/data/gm12878/hic/dat/%s.dat'%contig
  new_hic_dat_fpath = '/media/T02/data/gm12878/hic/dat/del-%s_%03d.dat'%(contig,t_n)
  
  chrom = Chromosome(origin_ref_obj, contig)
  chrom.load_cutsites(HIND3)
  derchrom = DelModifier(new_ref_obj, contig, chrom, sv_fpath)
  derchrom.load_cutsites(HIND3)

  shuffler = HiCReadShuffler(chrom, derchrom, resolute=resolute, debug_flag=True)
  print "Loaded HiC Struct"
  shuffler.simple_shuffle_to(hic_dat_fpath, new_hic_dat_fpath) 

def simulate_chromosome_debruijn():
  contig = 'chr20'
  t_n = 7
  hic_dat_fpath = '/scratchfast2/adp002/gm12878/hic/dat/%s.dat'%contig
  new_hic_dat_fpath = '/scratchfast2/adp002/gm12878/hic/dat/%s_%03d.dat'%(contig,t_n)
  #bed_fpath = '/home/adp002/Projects/assembly/data/chr20.hg18.hind3.bed'
  chrom = Chromosome(reference.Reference('/home/adp002/data/hg/hg18/%s.fa'%contig),contig)
  chrom.get_seq()
  chrom.load_cutsites(HIND3)
  #chrom.load_cutsites_from_bed(bed_fpath)

  shuffler = HiCReadShuffler(chrom, chrom, debug_flag=True)
  print "Loaded HiC Struct"
  shuffler.simple_shuffle_to(hic_dat_fpath, new_hic_dat_fpath) 
  
def simulate_del_chromosome_debruijn(contig, letter, outnum):
  origin_ref_obj = reference.Reference('/home/adp002/data/hg/hg18/%s.fa'%contig)  
  new_ref_obj = reference.Reference('/home/adp002/data/hg/grch38/chroms/%s.fa'%contig)  
  sv_fpath = '/home/adp002/Projects/assembly/data/simulations/sim.{0}.{1}.{2}.bed'.format(contig, outnum, letter)
  resolute = 100
  hic_dat_fpath = '/scratchfast2/adp002/gm12878/hic/dat/%s.dat'%contig
  new_hic_dat_fpath = '/scratchfast2/adp002/gm12878/hic/dat2/sim.{0}.{1:03d}.{2}.hic.dat'.format(contig,outnum, letter)
  
  chrom = Chromosome(origin_ref_obj, contig)
  chrom.load_cutsites(HIND3)
  derchrom = DelModifier(new_ref_obj, contig, chrom, sv_fpath)
  derchrom.load_cutsites(HIND3)

  shuffler = HiCReadShuffler(chrom, derchrom, resolute=resolute, debug_flag=True)
  print "Loaded HiC Struct"
  shuffler.simple_shuffle_to(hic_dat_fpath, new_hic_dat_fpath) 

# Deprecated.
def simulate_many_del_chromosome_debruijn(het_flag=False):
  from itertools import product, repeat
  contig = 'chr20'
  ref_obj = reference.Reference('/home/adp002/data/hg/hg18/%s.fa'%contig)  
  t_n = 1
  resolute = 10
  if het_flag:
    iter_num = product(['.A','.B'], range(0,5))
  else:
    iter_num = zip(repeat(''), range(0,5))
  # assume .A.bed and .B.bed terminology if het
  
  hic_dat_fpath = '/scratchfast2/adp002/gm12878/hic/dat/%s.dat'%contig
 
  chrom = Chromosome(ref_obj, contig)
  chrom.load_cutsites(HIND3)
  shuffler = None

  for het,t_i in iter_num:
    sv_fpath = '/home/adp002/Projects/assembly/data/simulations/sim.%s.%d%s.bed'%(contig, t_i,het)
    new_hic_dat_fpath = '/scratchfast2/adp002/gm12878/hic/dat/sim.del-%s_%03d%s.dat'%(contig,t_i,het)

    derchrom = DelModifier(ref_obj, contig, chrom, sv_fpath)
    chrom.load_cutsites(HIND3)
   
    if shuffler is None:
      shuffler = HiCReadShuffler(chrom, derchrom, resolute=resolute, debug_flag=False)
    else:
      shuffler.load_rand(derchrom, new_load_prefix_fpath)
    print "Loaded HiC Struct,", t_i
    shuffler.simple_shuffle_to(hic_dat_fpath, new_hic_dat_fpath) 

if __name__ == '__main__':
  import sys
  #test_zoning()
  #simulate_chromosome()
  simulate_del_chromosome()
  #contig, hap, outnum = sys.argv[1], sys.argv[2], sys.argv[3]
  #simulate_del_chromosome_debruijn(contig, hap, int(outnum))
  #test_zoning_large()
