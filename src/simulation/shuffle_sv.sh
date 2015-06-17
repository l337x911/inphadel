contig=$1;
outnum=$2;

seed=$(date +%N | sed -e 's/000$//' -e 's/^0//');

/home/anand/software/BEDTools-Version-2.12.0/bin/shuffleBed -i data/trio.2010_06.delsizes.1kb-300kb.${contig}.bed -seed ${seed} -excl data/gaps.hg38.bed -chrom -g data/hg38.genome | /home/anand/software/BEDTools-Version-2.12.0/bin/mergeBed | shuf -n 50 | /home/anand/software/BEDTools-Version-2.12.0/bin/sortBed -i stdin >data/simulations/sim.${contig}.${outnum}.bed;

sed -n 'p;n' data/simulations/sim.${contig}.${outnum}.bed >data/simulations/sim.${contig}.${outnum}.A.bed
sed -n 'n;p' data/simulations/sim.${contig}.${outnum}.bed >data/simulations/sim.${contig}.${outnum}.B.bed
