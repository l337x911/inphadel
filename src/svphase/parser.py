'''
Created on Jan 29, 2014

@author: anand
'''
import struct
import os
import numpy as np
from itertools import imap
from operator import itemgetter
from subprocess import Popen, PIPE

from svphase.utils.config import CONTIG_POS_TYPE, READLENGTH
from svphase.utils import SAM

class ReadsParserSAM(object):
	def __init__(self):
		self.p = SAM.SAMParser()

	def get_reads(self, fpath, contig, start=None, end=None):
		region = contig
		if not start is None:
			region = "%s:%d" % (region, start)
		if not end is None:
			region = "%s-%d" % (region, end)
	
		proc = Popen(['samtools', 'view', fpath, region], stdout=PIPE)
		# output = proc.communicate()[0]
		# print contig, start, end, 
		uniqify = set()
		for row in self.p.parse_file(proc.stdout):
		# for row in self.p.parse_string(output):
			# if int(row.tlen)<=0: continue
			if not row.qname in uniqify:
				yield row
			uniqify.add(row.qname)
		# print len(uniqify)
		proc.communicate()

	def get_strandedness(self, row):
		flags = int(row.flags)
		stranda = '-' if SAM.get_flag(SAM.SAM_FLAGS_H['rev_strand_of_query'], flags) else '+'
		strandb = '-' if SAM.get_flag(SAM.SAM_FLAGS_H['rev_strand_of_next_mate'], flags) else '+'
	
		if int(row.tlen) < 0:
			strandb, stranda = stranda, strandb
		return (stranda, strandb)
	def get_id(self, row):
		return row.qname

# Expect only one contig
class ReadsParserDat(object):
	def __init__(self):
		self.struct = struct.Struct('<LL?L?')
		self.current_fpath = None
		self.index = None
		self.n = None
		self.revc_to_strd = {True:'-', False:'+'}
		self.f = None
		
	def _index(self, fpath):
		self.current_fpath = fpath
		self.index = []
		self.n = os.path.getsize(fpath) / self.struct.size
		if os.path.isfile("%s.npz" % fpath):
			npzfile = np.load("%s.npz" % fpath)
			self.index = npzfile['index']
			assert self.n == len(self.index)
			npzfile.close()
			return
				
		with open(fpath, 'rb') as f:
			while True:
				try:
					self.index.append(self.struct.unpack(f.read(self.struct.size))[1])
				except struct.error:
					break
 
		self.index = np.array(self.index, dtype=CONTIG_POS_TYPE)
		# print self.index, self.n, self.struct.size
		#print self.index
		#print np.where(self.index[1:]-self.index[:-1]<0)
		assert all(self.index[i] <= self.index[i + 1] for i in xrange(len(self.index) - 1))
		np.savez("%s.npz" % fpath, index=self.index)

	def get_max_index_count(self, fpath):

		if self.f is None or self.f.closed:
			self.f = open(fpath, 'rb')
		self.f.seek(0)

		max_idx = 0
		row_size = self.struct.size
		while True:
			
			row_str = self.f.read(self.struct.size)
			if row_str == '':
				break

			row = self.struct.unpack(row_str)
			max_idx = max(row[0],max_idx)
		return max_idx
		
	
	def get_single_read_count(self, fpath):
		if fpath != self.current_fpath:
			self._index(fpath)
			if not self.f is None:
				self.f.close()

		# count single reads which did not map to the chromosome
		unmapped_count = np.searchsorted(self.index, 0, side='right')

		return self.n - unmapped_count
	
	def get_reads(self, fpath, contig, start=None, end=None):
		if fpath != self.current_fpath:
			self._index(fpath)
			if not self.f is None:
				self.f.close()

		if self.f is None or self.f.closed:
			self.f = open(fpath, 'rb')

		i, j = 0, self.n
		if not start is None:
			i = np.searchsorted(self.index, start - READLENGTH + 1, side='left')
		if not end is None:
			j = np.searchsorted(self.index, end + READLENGTH - 1, side='right')
		
		uniqify = set()
		# print i,j
		for pos in xrange(i, j):
			# print i,j, i*self.struct.size
			self.f.seek(pos * self.struct.size)
			row = self.struct.unpack(self.f.read(self.struct.size))

			if pos == j - 1: 
				self.f.close()

			# if int(row.tlen)<=0: continue
			if not row[0] in uniqify:
				yield row
			uniqify.add(row[0])
		self.f.close()

	def get_strandedness(self, row):
		idx, posa, revca, posb, revcb = row
		if posa > posb:
			posa, revca, posb, revcb = posb, revcb, posa, revca
		return self.revc_to_strd[revca], self.revc_to_strd[revcb]
	def get_id(self, row):
		return row[0]

class ReadsWriterDat(object):
	def __init__(self, fpath):
		self.fpath = fpath
		self.struct = struct.Struct('<LL?L?')
		self.reads = []
		self.idx = 0
	def write(self, posa, strda, posb, strdb):
		
		rev_ca, rev_cb = False, False
		if strda == '-':
			rev_ca = True
		if strdb == '-': 
			rev_cb = True

		self.reads.append((posa, self.struct.pack(self.idx, posa, rev_ca, posb, rev_cb)))
		self.reads.append((posb, self.struct.pack(self.idx, posb, rev_cb, posa, rev_ca)))
		self.idx += 1

	def close(self):
		self.reads.sort()
		# print [x for x in imap(itemgetter(1), self.reads)]
		with open(self.fpath, 'wb') as f:
			f.write(''.join(imap(itemgetter(1), self.reads)))

class ReadsWriterPredat(object):
	def __init__(self, fpath):
		self.fpath = fpath
		self.fmt = "{0:d}\t{1:d}\t{2:d}\t{3:d}\t{4:d}\n"
	def __enter__(self):
		self.f = open(self.fpath, 'wb')
		self.idx = 0
		return self
	def write(self, posa, strda, posb, strdb):
		rev_ca, rev_cb = False, False
		if strda == '-':
			rev_ca = True
		if strdb == '-': 
			rev_cb = True

		self.f.write(self.fmt.format(self.idx, posa, rev_ca, posb, rev_cb))
		self.f.write(self.fmt.format(self.idx, posb, rev_cb, posa, rev_ca))
		self.idx += 1

	def __exit__(self, tp, value, tb):
		self.f.close()

def test_dat_write():
	temp_out = 'te.dat'
	a = ReadsWriterDat(temp_out)

	a.write(50, '+', 100, '-')
	a.write(55, '-', 105, '-')
	a.write(60, '+', 300, '-')
	a.write(1000, '+', 100, '-')
	a.write(1002, '-', 999, '-')
	a.write(1, '+', 3000, '-')
	a.close()

	b = ReadsParserDat()
	print "100 to 901"
	for row in b.get_reads(temp_out, 'chr20', 100, 901):
		print row
		print b.get_strandedness(row)
	print b.index 
	
	print "101 to 500"
	for row in b.get_reads(temp_out, 'chr20', 101, 500):
		print row
		print b.get_strandedness(row)

	print "501 to 5000"
	for row in b.get_reads(temp_out, 'chr20', 501, 5000):
		print row
		print b.get_strandedness(row)

def write_hic_to_dat(read_parser, bam_fpath, dat_fpath, contig):
	write_parser = ReadsWriterDat(dat_fpath)

	for row in read_parser.get_reads(bam_fpath, contig):
		p, t = int(row.pos), int(row.tlen)
		if t <= 0: continue
	
		pa, pb = p, p + t - READLENGTH
		strda, strdb = read_parser.get_strandedness(row)
		write_parser.write(pa, strda, pb, strdb)				
	write_parser.close()

if __name__ == '__main__':
	import sys
	# Converts bam to dat
	if len(sys.argv) > 1:
		bam_fpath = sys.argv[1]
		dat_fpath = sys.argv[2]
		parser = ReadsParserSAM()
		write_hic_to_dat(parser, bam_fpath, dat_fpath, 'chr19')
	else:
		test_dat_write()
