import os
from svphase.parser import ReadsParserDat


def compute_stats_from_sample_dir(sim_dir, sample_dir, datadir):
	sample_all = os.path.join(sample_dir, datadir, 'all.stat')

	with open(sample_all, 'rb') as sf:
		sample_d = {}
		for line in sf:
			tokens = line.strip().split('\t')
			sample_d[tokens[0]] = tokens
	count_and_write(sample_d, sim_dir, datadir)

def compute_stats_from_faidx(sim_dir, fai_fpath, datadir):
	
	with open(fai_fpath, 'rb') as sf:
		sample_d = {}
		for line in sf:
			tokens = line.strip().split('\t')
			sample_d[tokens[0]] = tokens[:2] + ['0','0']
	count_and_write(sample_d, sim_dir, datadir)


def count_and_write(sample_d, sim_dir, datadir):
	sim_dir = os.path.join(sim_dir, datadir)
	sim_all = os.path.join(sim_dir, 'all.stat')

	dat_fmt = "{dir}/{c}.all.dat"
	r = ReadsParserDat()
	for contig in sample_d.keys():
		dat = dat_fmt.format(dir=sim_dir, c=contig)
		try:
			assert os.path.isfile(dat)
		except:
			print "Did not find contig for %s"%contig
			continue
		n = r.get_single_read_count(dat)
		sample_d[contig][2] = "{0:d}".format(n)

	with open(sim_all, 'wb') as sf:
		for k,tokens in sorted(sample_d.iteritems()):
			sf.write("\t".join(tokens)+"\n")



if __name__ == '__main__':
	import argparse
	
	parser = argparse.ArgumentParser(description="Script to make all.stat files for simulation")
	
	parser.add_argument('-f', '--faidx', dest='faidx', default=None, help='Fai file containing contig lengths')	
	parser.add_argument('-s', '--sample-dir', dest='sample_dir', default=None, help='Directory to Sample all.stat files')	
	parser.add_argument('sim_dir', help='Directory to Simulate all.stat files')

	args = parser.parse_args()
	assert os.path.isdir(args.sim_dir)

	if not (args.sample_dir is None):
		assert os.path.isdir(args.sample_dir)
		compute_stats_from_sample_dir(args.sim_dir, args.sample_dir, 'wgs')	
		compute_stats_from_sample_dir(args.sim_dir, args.sample_dir, 'hic')
	elif not (args.faidx is None):
		assert os.path.isfile(args.faidx)
		compute_stats_from_faidx(args.sim_dir, args.faidx, 'wgs')
		compute_stats_from_faidx(args.sim_dir, args.faidx, 'hic')
	else:
		raise Exception("No template file specified, requires fai or sample dir")
