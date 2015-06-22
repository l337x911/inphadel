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

class PrepTrainingFileStructure(object):
	def __init__(self, idir, ftype):
		assert ftype=='dat' or ftype=='bam'
		self.ftype = ftype
		if ftype == 'bam':
			self.idx_ftype = 'bai'
		elif ftype == 'dat':
			self.idx_ftype = 'npz'
		else:
			logger.error('No valid read filetype: %s', ftype)
			sys.exit(1)

		self.idir = idir
		self.sv_fpath = None
		self.wgs_read_count = None	
		self.hic_read_count = None	
		self.loci = None

	def load(self):
		""" Validate File Structure """
		if not os.path.isdir(os.path.join(self.idir, 'wgs')):
			logger.error('WGS directory not found: %s/wgs', self.idir)
			sys.exit(1)
		if not os.path.isdir(os.path.join(self.idir, 'hic')):
			logger.error('HiC directory not found: %s/hic', self.idir )
			sys.exit(1)
	
		self.sv_fpath = os.path.join(self.idir, 'truth','truth.bed')
		if not os.path.isfile(sv_fpath):
			logger.error('Truth bed file not found: %s', self.sv_fpath)
			sys.exit(1)

		self._set_read_counts()		
			
		self._set_loci()

		for contig in self.loci.contig.unique():
			if not os.path.isfile(os.path.join(self.idir, 'vcf', contig+'.vcf')):
				logger.error('Did not find VCF file for contig: %s/vcf/%s.vcf', self.idir, contig)
				sys.exit(1)
			self._check_for_contig('wgs', contig)	
			self._check_for_contig('hic', contig)	
	
	def _check_for_contig(self, data_src, contig):
		contig_all = os.path.join(self.idir, data_src, contig+'.all.{0}'.format(self.ftype))
		if not os.path.isfile(contig_all):
			logger.error('Missing %s contig : %s', data_src, contig_all)
			sys.exit(1)
		if not os.path.isfile('{0}.{1}'.format(contig_all_bam,self.idx_ftype)):
			logger.error('Missing %s contig index: %s.%s)', data_src, contig_all_bam, self.idx_ftype)
			sys.exit(1)
	
		contig_allele = os.path.join(self.idir, data_src, contig+'.{allele}.'+self.ftype)
		if not (os.path.isfile(contig_allele.format(allele='A')) and os.path.isfile(contig_allele.format(allele='B'))):
			logger.error('Missing allele split for %s on %s...',data_src, contig)

	def _set_read_counts(self):
		""" Check for idxstats file """
		wgs_stat_fpath = os.path.join(self.idir, 'wgs','all.stat')
		if not os.path.isfile(wgs_stat_fpath):
			logger.error('Samtools idxstats file not found: %s/wgs/all.stat', self.idir)
			sys.exit(1)
	
		self.wgs_read_count = pd.read_csv(wgs_stat_fpath, sep='\t', header=None, index_col=0).astype(int)
	
		hic_stat_fpath = os.path.join(self.idir, 'hic','all.stat')
		if not os.path.isfile(hic_stat_fpath):
			logger.error('Samtools idxstats file not found: %s/hic/all.stat', self.idir)
			sys.exit(1)
		self.hic_read_count = pd.read_csv(hic_stat_fpath, sep='\t', header=None, index_col=0).astype(int)
	
	def _set_loci(self):
		with open(self.sv_fpath, 'rb') as f:
			skiprows=0
			for line in f:
				if line.startswith('track'):
					skiprows+=1
				break	
	
		# Loads truth values into dataframe
		self.loci = pd.read_csv(self.sv_fpath, sep='\t', header=None, skiprows=skiprows, index_col=None).rename(columns={0:'contig',1:'start',2:'end',3:'truth'})
	
def main():
	import argparse
	import logging
	import time
	from pkg_resources import Requirement, resource_filename

	parser = argparse.ArgumentParser(description=__doc__)
	
	default_arguments(parser)
	parser.add_argument('version', help='version of training model generated')
	parser.add_argument('ftype', help='file format for reads', choices=['bam','dat'])
	# Sets the input file directories, model, reference, and debug
	parser.add_argument('input_dirs', nargs='+', help='directories containing input hic bam, wgs bam, and idxstat files.')

	parser.add_argument('--seed', type=int, default=None, help='Random initial state for some training procedures')
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

	# Process for each contig
	trainer = Trainer(model_type=args.model, model_pkl=model, k=args.k, random_state=random_state)
	#print loci
	for idir in args.input_dirs:	
		fs = PrepTrainingFileStructure(idir)
		fs.load()
		
		for contig in fs.loci.contig.unique():
			d = FileStructureData(os.path.join(idir, 'wgs'), os.path.join(idir, 'hic'), fs.sv_fpath, contig, ref_fpath=ref_fpath)

			adaptor = RPKMAdaptor(fs.wgs_read_count.loc[contig,1], fs.hic_read_count.loc[contig,1])	
			d.fill(adaptor=adaptor)
			trainer.add_data(d)

	with warnings.catch_warnings():
		warnings.simplefilter("ignore")
    
		acc = trainer.train()


if __name__ == '__main__':
	main()
