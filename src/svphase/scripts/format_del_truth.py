import sys

fpath = sys.argv[1]

out_svfpath = '/home/anand/Projects/assembly/data/gm12878/sid-del.bed'

def default_class(samples_str, allele):
  if allele=='NA':
    return 'inc'
  elif allele=='B':
    return 'pB'
  elif allele=='A':
    return 'pA'
  elif allele=='Both':
    return 'homA'

with open(fpath, 'rb') as f:
  for line in f:
    tokens = line.strip().split('\t')
    chrom, start, end, samples, allele = tokens
    print "chr{0}\t{1}\t{2}\t{3}".format(chrom, start, end, default_class(samples, allele))

