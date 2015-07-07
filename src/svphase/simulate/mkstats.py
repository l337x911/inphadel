import os
from svphase.parser import ReadsParserDat


def compute_stats(sim_dir, sample_dir, datadir):
	sim_dir = os.path.join(sim_dir, datadir)
	sample_all = os.path.join(sample_dir, datadir, 'all.stat')
	sim_all = os.path.join(sim_dir, 'all.stat')

	with open(sample_all, 'rb') as sf:
		sample_d = {}
		for line in sf:
			tokens = line.strip().split('\t')
			sample_d[tokens[0]] = tokens
	
	dat_fmt = "{dir}/{c}.all.dat"
	r = ReadsParserDat()
	for contig in sample_d.keys():
		dat = dat_fmt.format(dir=sim_dir, c=contig)
		assert os.path.isfile(dat)
		n = r.get_single_read_count(dat)
		sample_d[contig][2] = "{0:d}".format(n)

	with open(sim_all, 'wb') as sf:
		for k,tokens in sorted(sample_d.iteritems()):
			sf.write("\t".join(tokens)+"\n")

if __name__ == '__main__':
	import argparse
	
	parser = argparse.ArgumentParser(description="Script to make all.stat files for simulation")
	
	parser.add_argument('sample_dir', help='Directory to Sample all.stat files')	
	parser.add_argument('sim_dir', help='Directory to Simulate all.stat files')

	args = parser.parse_args()
	assert os.path.isdir(args.sample_dir)
	assert os.path.isdir(args.sim_dir)

	compute_stats(args.sim_dir, args.sample_dir, 'wgs')	
	compute_stats(args.sim_dir, args.sample_dir, 'hic')	
