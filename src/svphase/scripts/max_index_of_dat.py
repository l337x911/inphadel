import sys
from svphase.parser import ReadsParserDat 


def main(dat):
	r = ReadsParserDat()
	return r.get_max_index_count(dat)

if __name__ == '__main__':
	import argparse

	parser = argparse.ArgumentParser(description="Gets maximum index from a dat file")
	parser.add_argument('dat', help='path to DAT file')
	args = parser.parse_args()

	print main(args.dat)
		 
