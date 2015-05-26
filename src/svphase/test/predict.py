'''
Created on May 19, 2015

@author: root
'''
import unittest
from pkg_resources import Requirement, resource_filename

from svphase.learn.cov import FileStructureData, RPKMAdaptor
from svphase.learn.features import Features
#from svphase.utils.config import TEST_DATA_PREFIX
from svphase.inphadel import Predictor
import pandas as pd

TEST_DATA_PREFIX=resource_filename(Requirement.parse('InPhaDel'),'abbrev_test')
MODEL=resource_filename(Requirement.parse('InPhaDel'),'models/sim.del-chr19-20-hm.10-2c.10.RandomForest.pkl')
#TEST_DATA_PREFIX='abbrev_test'
class PredictionTest(unittest.TestCase):

	def testHasFeatures(self):
		#sv_fpath = SV_DATA+'/mills_2010-08-table4.CEU1kb.NA12878som.bed'
		sv_fpath = TEST_DATA_PREFIX + '/test_mills_chr21.bed'
		d = FileStructureData(TEST_DATA_PREFIX+'/wgs', TEST_DATA_PREFIX+'/hic', sv_fpath, 'chr21')
		d2 = FileStructureData(TEST_DATA_PREFIX+'/wgs', TEST_DATA_PREFIX+'/hic', sv_fpath, 'chr21')
		self.assertGreater(len(d.loci), 0, "No deletions to compute")
		
		d.fill(adaptor=None)
		d2.fill(adaptor=None)
		feats = Features([d,d2])
		#print feats.get_features()
		#print feats.get_nonzero()
		#print feats.get_nonzero_features().T
		self.assertGreater(len(feats.get_nonzero_features()), 0, "Has 0 nonzero features")
		
	def testPredict(self):
		#sv_fpath = SV_DATA+'/mills_2010-08-table4.CEU1kb.NA12878som.bed'
		sv_fpath = TEST_DATA_PREFIX + '/test_mills_chr21.bed'
		d = FileStructureData(TEST_DATA_PREFIX+'/wgs', TEST_DATA_PREFIX+'/hic', sv_fpath, 'chr21')
		d2 = FileStructureData(TEST_DATA_PREFIX+'/wgs', TEST_DATA_PREFIX+'/hic', sv_fpath, 'chr21')
		self.assertGreater(len(d.loci), 0, "No deletions to compute")
		
		adaptor = None
		wgs_read_count = pd.read_csv(TEST_DATA_PREFIX+'/wgs/all.stat', sep='\t', header=None, index_col=0).astype(int)
		hic_read_count = pd.read_csv(TEST_DATA_PREFIX+'/hic/all.stat', sep='\t', header=None, index_col=0).astype(int)
		adaptor = RPKMAdaptor(wgs_read_count.loc['chr19',1]+wgs_read_count.loc['chr20',1], hic_read_count.loc['chr19',1]+hic_read_count.loc['chr20',1])
		
		d.fill(adaptor=adaptor)
		d2.fill(adaptor=adaptor)
		predr = Predictor(MODEL)
		predr.add_data(d)
		predr.add_data(d2)
		#print feats.get_features()
		#print feats.get_nonzero()
		#print feats.get_nonzero_features().T
		preds = predr.predict_log_proba()
		#preds = predr.predict()
		preds = preds.sort_index()
		check_same_pred = {}
		for k,v in preds.iteritems():
			if check_same_pred.has_key(k):
				self.assertEqual(v, check_same_pred[k], "Predictions %s and %s differ for loci %s"%(v,check_same_pred[k],k))
		
		#print preds
def main():
	import argparse

		
	
if __name__ == "__main__":
	#import sys;sys.argv = ['', 'Test.testName']
	unittest.main()
