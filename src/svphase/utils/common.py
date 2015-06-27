import logging 
import sys

chrom_name_map = dict(zip(map("chr{0:d}".format, range(1,23)), range(1,23)))
chrom_name_map['chrX']=23 
chrom_name_map['chrY']=24 

def sort_chr(df):
  df['new_chr'] = df['chr'].map(chrom_name_map.get)
  df = df.sort(['new_chr','start','end']).drop(['new_chr',],axis=1)
  return df

logger = logging.getLogger('inphadel')

# Define Logging parameters"
logger.setLevel(logging.DEBUG)

# Log to Console
ch = logging.StreamHandler(sys.stdout)
#ch.setLevel(logging.INFO)
ch.setLevel(logging.DEBUG)

# Specifies format of log
ch.setFormatter(logging.Formatter('%(asctime)s [%(levelname) 9s] - %(message)s'))
logger.addHandler(ch)

def iter_paired_lines(stream):
	prev = None
	for idx, line in enumerate(stream):
		if idx%2==1:
			yield prev, line
		else:
			prev = line

	if idx%2==0:
		logger.error('Iter over paired lines ended on an odd line #')

