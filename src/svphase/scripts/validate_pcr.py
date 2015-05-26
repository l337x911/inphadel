import pandas as pd
import sys

class PCRTesting(object):
  def __init__(self, contig, start,end,snp_pos, window):
    self.contig = contig
    self.start = start
    self.end = end
    self.snp_pos = snp_pos
    self.window = window
    assert window<(end-start)
    self.regions = []
    self.labels = {}
  def set_regions(self):
    min_del_f = self.start-self.window
    max_del_r = self.end+self.window
    
    if self.snp_pos < self.start:
      self.regions.append((self.snp_pos-self.window/2,self.window/2,'forward','snp'))
      self.regions.append((self.start,self.window,'reverse','inner'))
      min_del_f = max(min_del_f, self.snp_pos)
    elif self.snp_pos > self.end:
      self.regions.append((self.end-self.window,self.window,'forward','inner'))
      self.regions.append((self.snp_pos,self.window/2,'reverse','snp'))
      max_del_r = min(max_del_r, self.snp_pos)
    else:
      raise Exception('SNP position falls within deletion')
    
    self.regions.append((min_del_f,self.start-min_del_f,'forward','outerF'))
    self.regions.append((self.end,max_del_r-self.end,'reverse','outerR'))
    
    s,e,o,name = zip(*self.regions)
    self.labels = dict([(j,i) for i,j in enumerate(name)])
  def out_regions(self, fpath=None):
    out = sys.stdout if fpath is None else open(fpath, 'wb')

    for s,e,o,l in self.regions:
      print >>out, "{0}\t{1}\t{2}\t{3}".format(self.contig, s,e,o)

    if not fpath is None:
      out.close()
  def out_region_links(self, fpath=None):
    out = sys.stdout if fpath is None else open(fpath, 'wb')
    
    print >>out, "{0}\t{1}\t{2}".format(self.labels['snp'], self.labels['inner'], 'snp-inner')
    c = 'outerR' if self.snp_pos<self.start else 'outerF'
    print >>out, "{0}\t{1}\t{2}".format(self.labels['snp'], self.labels[c], 'snp-outer')
    print >>out, "{0}\t{1}\t{2}".format(self.labels['outerF'], self.labels['outerR'], 'outerF-outerR')

    if not fpath is None:
      out.close()

def pcrable_regions(bed_fpath, w):
  ambre_prefix = '/home/anand/Projects/ambre/na12878p/hic_het_{0:02d}{1}'
  relabel_cols = {'X.1':'contig', 'X.2':'start', 'X.3':'end', 'X.8':'snp_pos', 'X.17':'dist'}
  dels = pd.read_csv(bed_fpath, sep='\t', header=None, index_col=None).rename(columns=relabel_cols)
  for idx, row in dels.iterrows():
    tester = PCRTesting(row['contig'], row['start'], row['end'], row['snp_pos'], w)
    print idx
    tester.set_regions()
    tester.out_regions(ambre_prefix.format(idx+1,'.txt'))
    tester.out_region_links(ambre_prefix.format(idx+1, '_links.txt'))
    #break

if __name__ == '__main__':
  pcrable_regions(sys.argv[1], w=1000)
