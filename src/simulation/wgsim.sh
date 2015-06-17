
prefix=$1;
letter=$2;
i=$3;

#chr19   63811651        46197064        976095
#chr20   62435964        51100213        695564

sleep 5;
python -m svphase.scripts.create_fa_from_sv ${prefix}.${i}.${letter}.bed;
#/home/anand/software/samtools-0.1.18/misc/wgsim ${prefix}.${i}.${letter}.fa -N 46200000 -1 100 -2 100 -h ${prefix}.${i}.${letter}.wgs_1.fq ${prefix}.${i}.${letter}.wgs_2.fq >${prefix}.${i}.${letter}.wgs.out ;
/home/anand/software/samtools-0.1.18/misc/wgsim ${prefix}.${i}.${letter}.fa -N 51100000 -1 100 -2 100 -h ${prefix}.${i}.${letter}.wgs_1.fq ${prefix}.${i}.${letter}.wgs_2.fq >${prefix}.${i}.${letter}.wgs.out ;
gzip ${prefix}.${i}.${letter}.wgs_2.fq;
mv ${prefix}.${i}.${letter}.wgs_1.fq /home/anand/dscratchfast2/gm12878/sim.wgs/; 
mv ${prefix}.${i}.${letter}.wgs_2.fq.gz /media/T02/data/gm12878/sim.wgs/;
echo ${prefix}.${i}.${letter}.fa;
sleep 2; 
