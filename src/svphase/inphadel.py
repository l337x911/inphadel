"""InPhaDel: Genotypes and phase deletions on a single chromosome using a specific classification model"""

import os
import sys
import pandas as pd
import numpy as np
import pickle
import argparse
import warnings
from itertools import chain, compress
from svphase.learn.cov import FileStructureData, RPKMAdaptor
from svphase.learn.features import Features
from svphase.scripts.sam_read_split import split_reads_by_allele
from svphase.utils.common import logger

__author__ = 'Anand Patel'
__company__ = "University of California, San Diego"
__email__ = "adp002@ucsd.edu"
__version__ = '1.0'

class ClassLabel(object):
	def __init__(self):
		self.classes = ['pA','pB','hom', 'inc']
		self.int_to_str_dict = dict(zip(range(len(self.classes)), self.classes))
		self.str_to_int_dict = dict(zip(self.classes, range(len(self.classes))))
	def is_deletion_on_a(self, c):
		if c=='pA' or c=='hom':
			return True
		else: 
			return False
	def is_deletion_on_b(self, c):
		if c=='pB' or c=='hom':
			return True
		else: 
			return False

class Predictor(object):
	def __init__(self, pickle_fpath, feature_subset=None):
		self.feature_subset = feature_subset
		with open(pickle_fpath, 'rb') as f:
			self.clf = pickle.load(f)
		print "Clf Params:"
		print self.clf.get_params()		
		self.data = []
		self.feats = None # Features Object

		self.labels = ClassLabel()
		
	def add_data(self, data):
		self.data.append(data)
	def get_loci(self):
		return list(chain(*[d.loci for d in self.data]))

	def check_features(self):
		assert len(self.data)>0
		features = self.feats.get_features()
		print "# of Features:", features.shape[1]
		print "# of Samples:", features.shape[0]
		print "# of Zero Samples:", sum(self.feats.get_nonzero()==False)
		
	def set_features(self):
		if self.feats is None:
			self.feats = Features(self.data, self.feature_subset)

	def get_pred_series(self, pred):
		loci = self.get_loci()

		zero_loci = list(compress(loci, np.logical_not(self.feats.get_nonzero())))
		s = pd.Series([self.labels.int_to_str_dict[i] for i in pred], index=compress(loci, self.feats.get_nonzero()))
		if len(zero_loci) > 0:
			s = s.append(pd.Series(['no_data',]*len(zero_loci), index=zero_loci))
		return s

	def get_pred_dataframe(self, pred):
		loci = self.get_loci()
		zero_loci = list(compress(loci, np.logical_not(self.feats.get_nonzero())))

		df = pd.DataFrame(pred, index=compress(loci,self.feats.get_nonzero())).rename(columns=self.labels.int_to_str_dict)
		nan_df = pd.DataFrame([[-np.inf]*len(df.columns),]*len(zero_loci), index=zero_loci).rename(columns=self.labels.int_to_str_dict)
		df['best'] = df.apply(lambda s:s.argmax(), axis=1)
		nan_df['best'] = 'no_data'

		df = pd.concat((df,nan_df))
		return df

	def predict(self):
		self.set_features()
		self.check_features()
		
		pred = self.clf.predict(self.feats.get_nonzero_features())
		return self.get_pred_series(pred) 
			
	def predict_log_proba(self):
		self.set_features()
		self.check_features()
		
		features = self.feats.get_nonzero_features()
		try:	
			pred = self.clf.predict_log_proba(features)
		except AttributeError:
			pred = np.log(self.clf.predict_proba(features))
		return self.get_pred_dataframe(pred)

def check_for_contig_bam_and_split(input_dir, data_src, contig):
	contig_all_bam = os.path.join(input_dir, data_src, contig+'.all.bam')
	if not os.path.isfile(contig_all_bam):
		logger.error('Missing %s contig bam: %s', data_src, contig_all_bam)
		sys.exit(1)
	if not os.path.isfile(contig_all_bam+'.bai'):
		logger.error('Missing %s contig index bam: %s.bai)', data_src, contig_all_bam)
		sys.exit(1)

	contig_allele_bam = os.path.join(input_dir, data_src, contig+'.{allele}.bam')
	if not (os.path.isfile(contig_allele_bam.format(allele='A')) and os.path.isfile(contig_allele_bam.format(allele='B'))):
		logger.warning('Generating allele split for %s on %s...',data_src, contig)
		split_reads_by_allele(os.path.join(input_dir, 'vcf', contig+'.vcf'), contig_all_bam, contig_allele_bam.format(allele='A'), contig_allele_bam.format(allele='B'))

	
def default_arguments(parser):
	parser.add_argument('model', default='rf', help='Model used to classify deletions (rf - RandomForest, svm - Support Vector Machine, nn - Nearest Neighbors)', choices=['rf','svm','knn'])
	#parser.add_argument('--test', dest='test', action='store_true', default=False, help='Run unit tests to ensure proper installation')
	parser.add_argument('-r,--reference', dest='reference_fasta', default=None, help='path to reference fasta (defaults to reference in config.py)')
	parser.add_argument('--debug', dest='debug', action='store_true', default=False, help=argparse.SUPPRESS)


def main():
	import logging
	from pkg_resources import Requirement, resource_filename

	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('bed', help='Simple bed file containing deletions (coordinates on reference)')
	parser.add_argument('input_dir', help='directory containing input hic bam, wgs bam, and idxstat files.')
	default_arguments(parser)
	parser.add_argument('-o,--out-csv', dest='out_csv', default='out.txt', help='tab-deliminted file to output phasing predictions (default: out.txt)')

	args = parser.parse_args()

	if args.debug:
		logger.setLevel(logging.DEBUG)
	else:
		logger.setLevel(logging.INFO)
	input_dir = args.input_dir
	sv_fpath = args.bed

	if args.model=='rf':
		model = resource_filename(Requirement.parse('InPhaDel'), 'models/sim.del-chr19-20-hm.10-2c.10.RandomForest.pkl')
	elif args.model=='svm':
		model = resource_filename(Requirement.parse('InPhaDel'), 'models/sim.del-chr19-20-hm.10-2c.10.SVM.pkl')
	elif args.model=='knn':
		model = resource_filename(Requirement.parse('InPhaDel'), 'models/sim.del-chr19-20-hm.10-2c.10.KN.pkl')
	else:
		logger.error('No classifying model found')
		sys.exit(1)


	out_csv = args.out_csv
	ref_fpath=args.reference_fasta

	# Validate File Structure
	if not os.path.isdir(os.path.join(input_dir, 'wgs')):
		logger.error('WGS directory not found: %s/wgs', input_dir)
		sys.exit(1)
	if not os.path.isdir(os.path.join(input_dir, 'hic')):
		logger.error('HiC directory not found: %s/hic', input_dir )
		sys.exit(1)
	if not (os.path.isfile(model) and model.endswith('.pkl')):
		logger.error('Classifier model not found: %s', model)
		sys.exit(1) 
	
	with open(sv_fpath, 'rb') as f:
		skiprows=0
		for line in f:
			if line.startswith('track'):
				skiprows+=1
			break	

	loci = pd.read_csv(sv_fpath, sep='\t', header=None, skiprows=skiprows, index_col=None).rename(columns={0:'contig',1:'start',2:'end'})

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

	# Process for each contig
	predr = Predictor(model)
	#print loci
	for contig in loci.contig.unique():
		if not os.path.isfile(os.path.join(input_dir, 'vcf', contig+'.vcf')):
			logger.error('Did not find VCF file for contig: %s/vcf/%s.vcf', input_dir, contig)
			sys.exit(1)
		check_for_contig_bam_and_split(input_dir, 'wgs', contig)	
		check_for_contig_bam_and_split(input_dir, 'hic', contig)	
		
		d = FileStructureData(os.path.join(input_dir, 'wgs'), os.path.join(input_dir, 'hic'), sv_fpath, contig, ref_fpath=ref_fpath)

		# Change in the future to per "contig" when models are trained that way
		adaptor = RPKMAdaptor(wgs_read_count.loc['chr19',1]+wgs_read_count.loc['chr20',1], hic_read_count.loc['chr19',1]+hic_read_count.loc['chr20',1])
		
		d.fill(adaptor=adaptor)
		predr.add_data(d)

	with warnings.catch_warnings():
		warnings.simplefilter("ignore")
    
		preds = predr.predict_log_proba()
		preds = preds.sort_index().reset_index()

		preds.to_csv(out_csv, sep='\t', float_format='%0.5f', index=False)
		

if __name__ == '__main__':
	main()
	
