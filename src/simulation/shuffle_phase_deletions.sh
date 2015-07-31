contig=$1;
datadir=$2;
outpath=$3;

seed=$(date +%N | sed -e 's/000$//' -e 's/^0//');

num=75

#/home/anand/software/BEDTools-Version-2.12.0/bin/shuffleBed -i ${datadir}/trio.2010_06.deletions.001kb_300kb.hg19.bed -seed ${seed} -excl ${datadir}/gaps.hg19.bed -chrom -g ${datadir}/hg19.genome | /home/anand/software/BEDTools-Version-2.12.0/bin/mergeBed | shuf -n 50 | /home/anand/software/BEDTools-Version-2.12.0/bin/sortBed -i stdin >${outpath}.1;
grep "${contig}"$'\t' ${datadir}/trio.2010_06.deletions.001kb_300kb.hg19.bed | /home/anand/software/BEDTools-Version-2.12.0/bin/shuffleBed -i stdin -seed ${seed} -excl ${datadir}/gaps.hg19.bed -chrom -g ${datadir}/hg19.genome | /home/anand/software/BEDTools-Version-2.12.0/bin/mergeBed | shuf -n ${num} | /home/anand/software/BEDTools-Version-2.12.0/bin/sortBed -i stdin >${outpath}.1;


shuf -n ${num} <(for i in {1..64}; do echo -e "pA\npB"; done;) >${outpath}.2;
paste ${outpath}.1 ${outpath}.2 >${outpath};
rm ${outpath}.1;
rm ${outpath}.2;

#sed -n 'p;n' data/simulations/sim.${contig}.${outnum}.bed >data/simulations/sim.${contig}.${outnum}.A.bed
#sed -n 'n;p' data/simulations/sim.${contig}.${outnum}.bed >data/simulations/sim.${contig}.${outnum}.B.bed
