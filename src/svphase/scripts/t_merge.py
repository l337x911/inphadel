import pandas as pd
import numpy as np

intersect_fpath = "/var/opt/data/NA12878_HG19_Bashir/truth/truth.1kb.intersect_NA12878.txt"
pred_fpath = "/var/opt/data/NA12878_HG19_Bashir/pred/rf.2.3_5.0.csv"

idf = pd.read_csv(intersect_fpath, sep='\t', index_col=None, header=None, names=['contig1','start1','end1','label1','contig','start','end','label2','overlap']).set_index(['contig','start','end'])
pdf = pd.read_csv(pred_fpath, sep='\t', index_col=None, header=0).set_index(['contig','start','end'])

df = pd.concat([idf,pdf],axis=1)

df = df[np.logical_and(pd.notnull(df['best']),pd.notnull(df['label1']))]
df = df[df['best']!='no_data']
print df[['label1','label2','best']]

print "{0}/{1}".format(sum(df['label1']==df['best']), len(df['best']))

df = df[df['label2']=='H']
print "{0}/{1}".format(sum(df['label1']==df['best']), len(df['best']))
