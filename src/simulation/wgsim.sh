all_stat=$1;
fa=$2;
outdir=$3;
prev_index=$4;

#chr19   63811651        46197064        976095
#chr20   62435964        51100213        695564

chrom=$(basename ${fa:0:${#fa}-3});
contig=${chrom:0:4};
index_start=$((prev_index+1));

# Simulate half the number of paired end reads.
template_count=$(grep -e "${contig}"$'\t' ${all_stat} | awk '{print int($3/4);}');
echo $chrom $contig $template_count
echo $index_start

rm ${outdir}/${chrom}_1.fq
mkfifo ${outdir}/${chrom}_1.fq
sed -n 1~4p ${outdir}/${chrom}_1.fq | python -m svphase.scripts.fq_to_predat ${index_start} | gzip 1>${outdir}/${chrom}.predat.gz &

/home/anand/software/samtools-0.1.18/misc/wgsim -N ${template_count} -1 100 -2 100 ${fa} ${outdir}/${chrom}_1.fq /dev/null >${outdir}/${chrom}.out;
sleep 5;


#python -m svphase.scripts.create_fa_from_sv ${prefix}.${i}.${letter}.bed;
##/home/anand/software/samtools-0.1.18/misc/wgsim ${prefix}.${i}.${letter}.fa -N 46200000 -1 100 -2 100 -h ${prefix}.${i}.${letter}.wgs_1.fq ${prefix}.${i}.${letter}.wgs_2.fq >${prefix}.${i}.${letter}.wgs.out ;
#/home/anand/software/samtools-0.1.18/misc/wgsim ${prefix}.${i}.${letter}.fa -N 51100000 -1 100 -2 100 -h ${prefix}.${i}.${letter}.wgs_1.fq ${prefix}.${i}.${letter}.wgs_2.fq >${prefix}.${i}.${letter}.wgs.out ;
#gzip ${prefix}.${i}.${letter}.wgs_2.fq;
#mv ${prefix}.${i}.${letter}.wgs_1.fq /home/anand/dscratchfast2/gm12878/sim.wgs/; 
#mv ${prefix}.${i}.${letter}.wgs_2.fq.gz /media/T02/data/gm12878/sim.wgs/;
#echo ${prefix}.${i}.${letter}.fa;
#sleep 2; 
