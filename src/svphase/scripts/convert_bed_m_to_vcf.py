import pandas as pd
import sys

def convert(vcf_fpath, bed_fpath):

  vcf_header = ['contig','snp_id', 'ref', 'don', 'score', 'valid', 'info', 'format', 'sample']
  bed_header = ['contig', 'new_pos', 'new_posPlus']

  vcf = pd.read_csv(vcf_fpath, sep='\t', header=None, index_col=1, names=vcf_header)
  #print vcf.ix[:50000,:]

  bed = pd.read_csv(bed_fpath, sep='\t', header=None, index_col=3, names=bed_header)

  contig = vcf['contig'].values[0]
  #print contig

  #print len(vcf.index.intersection(bed.index))
  #print len(bed.index)

  #vcf['new_pos'] = 0
  #print vcf.ix[bed.index, :] 
  #vcf.ix[bed.index, ['new_pos',]] = bed['new_pos'].astype(int)

  #print vcf.columns 
  new_vcf = pd.concat([bed[['contig','new_pos']], vcf[vcf_header[1:]]], axis=1)
  new_vcf['new_pos'] = new_vcf['new_pos'].astype(int)

  new_header = ['contig', 'new_pos'] + vcf_header[1:]
  new_vcf.ix[bed.index,new_header].to_csv(sys.stdout, sep='\t', index=False, header=False)


if __name__ == '__main__':
  vcf_fpath = sys.argv[1]
  bed_fpath = sys.argv[2]
  convert(vcf_fpath, bed_fpath)
