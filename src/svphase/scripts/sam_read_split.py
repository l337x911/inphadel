import pysam
import pandas as pd
import os.path

def get_next_pileup_of_reads(reads, v):
	for p in reads.pileup(v.chrom, v.pos-1, v.pos, truncate=True):
		if p.pileups is None: 
			return iter([])
		return p.pileups
	return iter([])

def split_reads_by_allele(vcf_path, bam_path, out_bamA_path, out_bamB_path):
	""" Splits reads in BAM file that overlap phased variants in vcf file to bamA and bamB """

	#  Loops of vcf files
	comment_count = 0 
	with open(vcf_path,'rb') as f:
		for line in f:
			if not line.startswith('#'):
				break
			comment_count += 1
	
	vcf = pd.read_csv(vcf_path, sep='\t', skiprows=comment_count, header=None, index_col=False)
	vcf = vcf.rename(columns=dict(zip(range(len(vcf.columns)),['chrom', 'pos', 'id','ref','alt','qual','filter','info','format']+range(0,len(vcf.columns)-9))))
	vcf.loc[:,'pos'] = vcf.loc[:,'pos'].astype(int)
	vcf.loc[:,'ref'] = vcf.loc[:,'ref'].apply(str.upper)
	vcf.loc[:,'alt'] = vcf.loc[:,'alt'].apply(str.upper)

	#print vcf.loc[:,0].apply(lambda x:x[1]=='|')
	phased_vcf = vcf.loc[vcf.loc[:,0].apply(lambda x:x[1]=='|'),:]

	with pysam.AlignmentFile(bam_path, 'rb') as reads, pysam.AlignmentFile(out_bamA_path, 'wb', template=reads) as out_a, pysam.Samfile(out_bamB_path, 'wb', template=reads) as out_b:
		 for idx, v in phased_vcf.iterrows():
			#print v
			is_alt_variant_on_a = v[0]=='1|0'
			if not ( is_alt_variant_on_a or v[0]=='0|1'):
				raise Exception('Variant %s,%d does not have a valid genotype %s'%(v.chrom,v.pos,v[0]))

			for p_read in get_next_pileup_of_reads(reads,v):
				if p_read.qpos is None or p_read.alignment.cigar is None: continue
				qpos = int(p_read.qpos)
				cigar_type, cigar_count = p_read.alignment.cigar[0]
				if cigar_type==4:
					qpos -= cigar_count

				try:	
					allele = p_read.alignment.query[qpos]
				except IndexError:
					#continue
					raise Exception('%s Alignment %s (%d) does not have a query position %d'%(p_read.alignment.qname, p_read.alignment.query, p_read.alignment.alen, p_read.qpos))
				is_ref_allele = allele==v.ref
				is_alt_allele = allele==v.alt
				if not (is_ref_allele or is_alt_allele): continue

				if is_alt_allele:
					if is_alt_variant_on_a:
						out_a.write(p_read.alignment)
					else:
						out_b.write(p_read.alignment)
				else:
					# read is a reference match
					if is_alt_variant_on_a:
						out_b.write(p_read.alignment)
					else:
						out_a.write(p_read.alignment)
						

def file_check(parser, arg):
	if not os.path.exists(arg):
		parser.error('The file %s does not exist!'%arg)
	else:
		return arg

if __name__ == '__main__':
	import argparse
	out_a = '_Aallele.bam'
	out_b = '_Ballele.bam'

	parser = argparse.ArgumentParser(description='Partitions reads from a BAM file into a bam file with reads A and B')
	parser.add_argument('vcf', help='VCF 4.2 containing scaffold of phased variants (allele A | allele B)')
	parser.add_argument('bam', help='BAM file to split reads into Aallele and Ballele')
	parser.add_argument('out_prefix', default='prefix', help='prefix to path to output bam files (default: prefix_{a} and prefix_{b})'.format(a=out_a, b=out_b))
	
	args = parser.parse_args()
	split_reads_by_allele(args.vcf, args.bam, args.out_prefix+out_a, args.out_prefix+out_b) 
	
