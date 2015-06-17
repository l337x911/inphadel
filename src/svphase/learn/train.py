"""InPhaDel: Genotypes and phase deletions on a single chromosome using a specific classification model

Trains models for phasing deletions using underlying WGS+HiC data 

"""

import pickle
import pandas as pd
import numpy as np
import warnings

from sklearn import svm, ensemble
from sklearn.cross_validation import StratifiedKFold, cross_val_score
from sklearn.grid_search import GridSearchCV
from sklearn.naive_bayes import GaussianNB
from sklearn.neighbors import KNeighborsClassifier
from sklearn.decomposition import PCA


from svphase.utils.config import RANDOM_STATE,DATA_PREFIX
from svphase.learn.cov import FileStructureData, RPKMAdaptor
from svphase.learn.features import Features
from svphase.utils.common import logger
from svphase.inphadel import default_arguments

class Trainer(object):
	def __init__(self, model_type='', model_pkl='', k=5, feature_subset=None): 
		self.data = []
		self.feats = None
		
	def add_data(self, data):
		self.data.append(data)
	def set_features(self):
		if self.feats is None:
			self.feats = Features(self.data, self.feature_subset)

def main():
	import argparse
	import logging
	import time
	from pkg_resources import Requirement, resource_filename

	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('--seed', type=int, default=None, help='Random initial state for some training procedures')
	
	default_arguments(parser)
	parser.add_argument('version', help='Version of training model generated.')
	# Sets the input file directories, model, reference, and debug
	parser.add_argument('input_dirs', nargs='+', help='directories containing input hic bam, wgs bam, and idxstat files.')

	parser.add_argument('-k', default=5, help='# of folds in nested cross validation')	

	args = parser.parse_args()

	if args.debug:
		logger.setLevel(logging.DEBUG)
	else:
		logger.setLevel(logging.INFO)

	random_state=int(time.time())%1000000


	if args.model=='rf':
		model = resource_filename(Requirement.parse('InPhaDel'), 'models/RF.{0}.pkl'.format(version))
	elif args.model=='svm':
		model = resource_filename(Requirement.parse('InPhaDel'), 'models/SVM.{0}.pkl'.format(version))
	elif args.model=='knn':
		model = resource_filename(Requirement.parse('InPhaDel'), 'models/KN.{0}.pkl'.format(version))
	else:
		logger.error('No training model chosen')
		sys.exit(1)

	ref_fpath=args.reference_fasta

	input_dir = args.input_dir
	# Validate File Structure
	if not os.path.isdir(os.path.join(input_dir, 'wgs')):
		logger.error('WGS directory not found: %s/wgs', input_dir)
		sys.exit(1)
	if not os.path.isdir(os.path.join(input_dir, 'hic')):
		logger.error('HiC directory not found: %s/hic', input_dir )
		sys.exit(1)
	
	with open(sv_fpath, 'rb') as f:
		skiprows=0
		for line in f:
			if line.startswith('track'):
				skiprows+=1
			break	

	# Loads truth values into dataframe
	loci = pd.read_csv(sv_fpath, sep='\t', header=None, skiprows=skiprows, index_col=None).rename(columns={0:'contig',1:'start',2:'end',3:'truth'})

	# Check for idxstats file
	wgs_stat_fpath = os.path.join(input_dir, 'wgs','all.stat')
	if not os.path.isfile(wgs_stat_fpath):
		logger.error('Samtools idxstats file not found: %s/wgs/all.stat', input_dir)
		sys.exit(1)

	wgs_read_count = pd.read_csv(wgs_stat_fpath, sep='\t', header=None, index_col=0).astype(int)

	hic_stat_fpath = os.path.join(input_dir, 'hic','all.stat')
	if not os.path.isfile(hic_stat_fpath):
		logger.error('Samtools idxstats file not found: %s/hic/all.stat', input_dir)
		sys.exit(1)
	hic_read_count = pd.read_csv(hic_stat_fpath, sep='\t', header=None, index_col=0).astype(int)

	
	sv_fpath = args.bed

	# Process for each contig
	trainer = Trainer(model_type=args.model, model_pkl=model, k=args.k, random_state=random_state)
	#print loci
	for contig in loci.contig.unique():
		if not os.path.isfile(os.path.join(input_dir, 'vcf', contig+'.vcf')):
			logger.error('Did not find VCF file for contig: %s/vcf/%s.vcf', input_dir, contig)
			sys.exit(1)
		check_for_contig_bam_and_split(input_dir, 'wgs', contig)	
		check_for_contig_bam_and_split(input_dir, 'hic', contig)	
		
		d = FileStructureData(os.path.join(input_dir, 'wgs'), os.path.join(input_dir, 'hic'), sv_fpath, contig, ref_fpath=ref_fpath)

		adaptor = RPKMAdaptor(wgs_read_count.loc[contig,1], hic_read_count.loc[contig,1])
		
		d.fill(adaptor=adaptor)
		trainer.add_data(d)

	with warnings.catch_warnings():
		warnings.simplefilter("ignore")
    
		acc = trainer.train()


if __name__ == '__main__':
	main()
