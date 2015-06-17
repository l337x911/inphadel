#!/bin/sh -u
# Expect's sam file
fname=$1;
newdir=$2;
predatdir=$3;
required="0x0042"

tmpfname=$(basename ${fname});
newfname=${tmpfname:0:${#tmpfname}-4};

echo ${newfname}
#samtools view -f ${required} ${fname} | cut -f 2,4,7,8 | grep "=" | awk '{print NR-1 "\t" $2 "\t" and($1,0x0010)/0x0010 "\t" $4 "\t" and($1,0x0020)/0x0020; print NR-1 "\t" $4 "\t" and($1,0x0020)/0x0020 "\t" $2 "\t" and($1,0x0010)/0x0010;}' >${newdir}/$(basename ${fname}).predat;
#sort -n -k 2,2b $(basename ${fname}).predat | python -m svphase.scripts.predat_to_dat >${newdir}/$(basename ${fname}).0_0.dat 

samtools view -f ${required} ${fname} | cut -f 2,4,7,8 | grep "=" | awk 'BEGIN {x=0;} {print x "\t" $2 "\t" and($1,0x0010)/0x0010 "\t" $4 "\t" and($1,0x0020)/0x0020; print x "\t" $4 "\t" and($1,0x0020)/0x0020 "\t" $2 "\t" and($1,0x0010)/0x0010; x+=1;}' | tee ${predatdir}/${newfname}.0_0.predat | sort -n -S 5G -k 2,2b | python -m svphase.scripts.predat_to_dat >${newdir}/${newfname}.0_0.dat 
