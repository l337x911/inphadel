#!/bin/sh
c=$1;
# should be "wgs" or "hic"
data=$2;
# should be chromA or chromB
del=$3;

cat <(gzip -c -d /root/data/sim_deletions_0_${del}/${data}/${c}.chromA.ref.predat.gz) <(gzip -c -d /root/data/sim_deletions_0_${del}/${data}/${c}.chromB.ref.predat.gz) | python -m svphase.scripts.predat_to_dat  >/root/data/sim_deletions_0_${del}/${data}/${c}.all.dat
