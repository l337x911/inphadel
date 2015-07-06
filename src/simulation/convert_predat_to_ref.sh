#!/bin/sh
c=$1;
# should be "wgs" or "hic"
data=$2;
# should be chromA or chromB
hap=$3;
del=$4;

gzip -c -d /root/data/sim_deletions_0_${del}/${data}/${c}.${hap}.predat.gz | python -m svphase.simulate.convert_predat_to_ref ${c} /root/data/sim_deletions_0_${del}/truth/truth.${hap:5:1}.bed /home/anand/data/hg/19/hg19.fa | gzip -c >/root/data/sim_deletions_0_${del}/${data}/${c}.${hap}.ref.predat.gz;
