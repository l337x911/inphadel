import numpy as na
from matplotlib import pyplot as plt

hic_sv_fpath = '/home/anand/Projects/assembly/data/gm12878/del1kb.chr20o.hic.bed'
wgs_sv_fpath = '/home/anand/Projects/assembly/data/gm12878/del1kb.chr20.wgs.bed'
hic_nonsv_fpath = '/home/anand/Projects/assembly/data/gm12878/nondel1kb.rand.chr20o.hic.bed'
wgs_nonsv_fpath = '/home/anand/Projects/assembly/data/gm12878/nondel1kb.rand.chr20.wgs.bed'

class GenomeDataset(object):
  def __init__(self, x, c):
    self.x = x
    self.c = c

def get_cov_per_base(fpath):
  ratios = []
  with open(fpath, 'rb') as f:
    for line in f:
      tokens = line.strip().split('\t')
      ratios.append(float(tokens[3])/float(tokens[5]))
  return ratios

fig = plt.figure()
ax1 = fig.add_subplot(111)

hic = GenomeDataset(+0.25, 'r')
wgs = GenomeDataset(-0.25,  'b')

bg = 0.2
fg = 0.8

markersize = 60
# background 
hic_ratios = get_cov_per_base(hic_nonsv_fpath)
l = len(hic_ratios)
ax1.scatter(na.ones(l)+hic.x*na.random.random(l), hic_ratios, s=markersize, alpha=bg, linewidths=0, c=hic.c)
wgs_ratios = get_cov_per_base(wgs_nonsv_fpath)
l = len(wgs_ratios)
ax1.scatter(na.ones(l)+wgs.x*na.random.random(l), wgs_ratios, s=markersize, alpha=bg, linewidths=0, c=wgs.c)

hic_ratios = get_cov_per_base(hic_sv_fpath)
l = len(hic_ratios)
ax1.scatter(na.ones(l)+hic.x*na.random.random(l)+1, hic_ratios, s=markersize, alpha=fg, linewidths=0, c=hic.c, label='HiC')
wgs_ratios = get_cov_per_base(wgs_sv_fpath)
l = len(wgs_ratios)
ax1.scatter(na.ones(l)+wgs.x*na.random.random(l)+1, wgs_ratios, s=markersize, alpha=fg, linewidths=0, c=wgs.c, label='WGS')

ax1.set_ylabel('Reads per bp\n (Background gives expected distribution)', horizontalalignment='center')

plt.title('chr20 deletions called by Mills+Conrad')
plt.legend()
plt.xticks((1,2), ('Background', 'SV'))
plt.show()

