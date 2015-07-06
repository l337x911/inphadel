#!/bin/sh
c=$1;
# should be "wgs" or "hic"
data=$2;
# should be chromA or chromB
hap=$3;
del=$4;

gzip -c -d /root/data/sim_deletions_0_${del}/${data}/${c}.${hap}.ref.predat.gz | python -m svphase.simulate.filter_predat_by_vcf /root/data/sim_deletions_0_${del}/vcf/${c}.vcf | python -m svphase.scripts.predat_to_dat >/root/data/sim_deletions_0_${del}/${data}/${c}.${hap:5:1}.dat
