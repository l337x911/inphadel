import sys
import pandas as pd
from svphase.scripts.hm_trios import SNPTransmit

def dels_to_ase():
  # reads sv intervals from stdin
  out_bed_fpath = '/home/anand/Projects/assembly/data/ase_sites.bed'
  ase_genes = pd.read_csv('/home/anand/Projects/assembly/data/ASE_genes.txt', sep='\t', names=['gene','inherit']).set_index('gene').dropna()
  gb_genes = pd.read_csv('/home/anand/Projects/assembly/data/knownGene.hg18.txt', sep='\t', header=0).set_index('name')
  ase_genes.index.name = 'name'
  #print len(ase_genes.index.intersection(gb_genes.index))
  gb_genes = pd.concat((gb_genes, ase_genes), axis=1).dropna(axis=0, subset=['inherit']).rename(columns={'chrom':'chr'})
  snp = SNPTransmit('NA12891','NA12892')

  pat_v = gb_genes.ix[gb_genes['inherit']=='pat',:]
  snp.assign_dels_to_homolog(pat_v,'pat')

  mat_v = gb_genes.ix[gb_genes['inherit']=='mat',:]
  snp.assign_dels_to_homolog(mat_v,'mat')
  
  final = pd.concat([pat_v, mat_v])

  chrom_name_map = dict(zip(map("chr{0:d}".format, range(1,23)), range(1,23)))
  chrom_name_map['chrX']=23 
  chrom_name_map['chrY']=24 

  final['new_chr'] = final['chr'].map(chrom_name_map.get)
  final = final.sort(['new_chr','txStart','txEnd']).drop(['new_chr',],axis=1)

  #gb_genes = gb_genes.ix[ase_genes.index,:]
  #gb_genes['ASE'] = ''
  #print ase_genes.ix[:20,:]
  #gb_genes.ix[ase_genes.index,'ASE'] = ase_genes

  #print gb_genes.ix[:10,['inherit','chrom','txStart','txEnd']]

  final[['chr','txStart','txEnd', 'class']].to_csv(out_bed_fpath, header=False, index=False, sep='\t')

if __name__ == '__main__':
  dels_to_ase()
