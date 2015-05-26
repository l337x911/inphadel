import SAM
from collections import defaultdict
from subprocess import Popen, PIPE
from itertools import chain
import numpy as na
#from matplotlib import pyplot as plt

bam_fpath = '/media/T01/data/gm12878/hic/NA12878.chr20.hiseq.wgs.bwa.recal.srtd.rmdup.eq.bam'
#bam_fpath = '/media/T01/data/gm12878/hic/chr20o.all.srtd.rmdup.bam'
#sv_fpath = '/home/anand/Projects/assembly/data/gm12878/nondel1kb.rand1kb.chr20.bed'
#sv_fpath = '/home/anand/Projects/assembly/data/gm12878/nondel1kb.rand1kb.chr20.srtd.bed'
sv_fpath = '/home/anand/Projects/assembly/data/gm12878/del1kb.1kb-around.chr20.bed'
#sv_fpath = '/home/anand/Projects/assembly/data/gm12878/del1kb.500b-arounds.chr20.bed'
#sv_fpath = '/home/anand/Projects/assembly/data/gm12878/del1kb.500b-around.chr20.bed'
#sv_fpath = '/home/anand/Projects/assembly/data/gm12878/del1kb.200b-around.chr20.bed'
#sv_fpath = '/home/anand/Projects/assembly/data/gm12878/del1kb.chr20.bed'

orients_pe = {'++':[], '+-':[],'-+':[],'--':[]}

def get_reads(fpath, contig, start=None, end=None):
  p = SAM.SAMParser()

  region = contig
  if not start is None:
    region = "%s:%d"%(region,start)
  if not end is None:
    region = "%s-%d"%(region,end)

  proc = Popen(['samtools', 'view', fpath, region], stdout=PIPE)
  output = proc.communicate()[0]
  
  uniqify = set()
  for row in p.parse_string(output):
    #if int(row.tlen)<=0: continue
    if not row.qname in uniqify:
      yield row
    uniqify.add(row.qname)


loci = []
with open(sv_fpath, 'rb') as f:
  for line in f:
    if line.startswith('track'): continue
    tokens = line.strip().split('\t')
    loci.append((tokens[0], int(tokens[1]), int(tokens[2])))

def line_format(vlist):
  return '\t'.join(["% 10d"%i for i in vlist])

def float_format(vlist):
  return ' '.join(["% 11d"%vlist[0],]+["% 1.5f"%i for i in vlist[1:]])

def get_strandedness(row):
  flags = int(row.flags)
  stranda = '-' if SAM.get_flag(SAM.SAM_FLAGS_H['rev_strand_of_query'], flags) else '+'
  strandb = '-' if SAM.get_flag(SAM.SAM_FLAGS_H['rev_strand_of_next_mate'], flags) else '+'

  if int(row.tlen)<0:
    strandb,stranda = stranda, strandb
  return (stranda, strandb)

#bins = na.concatenate((na.linspace(1,1000,num=5,endpoint=False),na.linspace(1000,10000,num=5,endpoint=False),na.logspace(4,8,num=5)))
#bins = na.logspace(0,5,num=6)
bins = na.concatenate((na.linspace(1,1500,num=40), [10**9]))
#bins = na.concatenate((-1*bins[::-1], bins))

orients = dict(zip(['+-','-+','--','++'], range(4)))

loci = zip(loci[::2], loci[1::2])
counts_arr = na.zeros((4,len(loci), len(bins)-1), dtype=na.float)

#for idx, (c1,a1,b1) in enumerate(loci):
for idx, ((c1,a1,b1),(c2,a2,b2)) in enumerate(loci):
  counts = defaultdict(list)

  leftlink = dict([(row.qname, row) for row in get_reads(bam_fpath, c1, a1, b1)])
  rightlink = dict([(row.qname, row) for row in get_reads(bam_fpath, c2, a2, b2)])
  #for row in leftlink.itervalues():
  for row in map(leftlink.__getitem__, set(leftlink.viewkeys()).intersection(set(rightlink.viewkeys()))):
    #print row
    #counts['%s%s'%get_strandedness(row)].append(int(row.tlen))
    if int(row.tlen)<0: continue
    counts['%s%s'%get_strandedness(row)].append(int(row.tlen)-(a2-b1))
  
  for v in counts.itervalues():
    v.sort()
  print c1, a1,b1, a2,b2, counts
  counts_hist = {}
  #print line_format(bins)
  for k,v in counts.items():
    counts_hist[k] = na.histogram(v,bins=bins)[0]
    counts_arr[orients[k],idx,:] = counts_hist[k]

#print counts_arr.shape

orients_srtd = sorted(orients.items(), key=lambda x:x[1])
orients_str = " ".join(zip(*orients_srtd)[0])
print orients_str
read_counts = na.sum(na.sum(counts_arr, axis=2), axis=0)

counts_arr2 = counts_arr[:,read_counts>0,:]
reado_counts = read_counts[read_counts>0].reshape((1,na.sum(read_counts>0)))

orients_arr = na.sum(counts_arr2, axis=2)
orients_narr = na.concatenate([reado_counts]*4,axis=0)
print na.mean(orients_arr/orients_narr, axis=1)
print na.std(orients_arr/orients_narr, axis=1)
print read_counts
out_table = na.zeros((len(bins)-1, 9))

#print line_format(bins)
#print out_table.shape, bins.shape
out_table[:,0] = bins[1:]

for idx, (k,i) in enumerate(orients_srtd):
  readd_counts = na.sum(counts_arr[i,:,:], axis=1)
  counts_arr3 = counts_arr[i,readd_counts>0,:]
  readd_counts = readd_counts[readd_counts>0]
  counts_narr3 = na.concatenate([readd_counts.reshape(readd_counts.size,1)]*(len(bins)-1), axis=1)
  #print counts_arr3.shape, counts_narr3.shape, readd_counts.shape
  #print k
  #print float_format(na.mean(counts_arr3/counts_narr3, axis=0))
  #print float_format(na.std(counts_arr3/counts_narr3, axis=0))
  out_table[:,idx+1] = na.mean(counts_arr3/counts_narr3, axis=0)
  out_table[:,idx+5] = na.std(counts_arr3/counts_narr3, axis=0)
 
print "Bins: mean %s std %s"%(orients_str,orients_str)
for idx in xrange(len(bins)-1):
  print float_format(out_table[idx,:])


