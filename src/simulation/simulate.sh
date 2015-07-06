count=$1; # Should be 00?
threads=$2;
# expected truth, and vcf directories to have been created.
# expects vcf files and original NA12878 dataset

mkdir -p /root/data/sim_deletions_0_${count}/truth;
mkdir -p /root/data/sim_deletions_0_${count}/wgs;
mkdir -p /root/data/sim_deletions_0_${count}/hic;
ln -s /root/data/NA12878_HG19_GS_0_1/vcf /root/data/sim_deletions_0_${count}/vcf;

cd /root/data/sim_deletions_0_${count}/truth;
echo "Working on truth...";
for c in chr2 chr3 chr4; 
do 
if [ ! -f truth.${c}.bed ];
then 
	bash ~/Projects/inphadel-min/src/simulation/shuffle_deletions.sh ${c} ~/Projects/inphadel-min/data truth.${c}.bed; 
fi
done;

if [ ! -f truth.bed ];
then 
	cat truth.chr?.bed >truth.bed;
	grep -e "pA" -e "hom" truth.bed >truth.A.bed;
	grep -e "pB" -e "hom" truth.bed >truth.B.bed;
fi

#Verify bedsize distribution
bed=truth.bed && awk '{print $3-$2;}' ${bed} | sort -n -k 1,1 | awk -v wc="$(wc -l ${bed})" '{print $0, NR/wc;}' >${bed}.cumsum
awk 'BEGIN {x=0;} {if($3<0){print x, $0;}}' truth.bed >truth.bed.overlap

echo "Working on creating FASTA...";

cd /home/anand/Projects/inphadel-min/src
# Create chrom A fasta sequence with hom and pA deletions 
# Create chrom B fasta sequence with hom and pB deletions
for c in chr2 chr3 chr4; 
do 
	if [ ! -f /root/data/sim_deletions_0_${count}/truth/${c}.chromA.fa ];
	then
		python -m svphase.scripts.create_deletion_fa_from_bed ${c} /root/data/sim_deletions_0_${count}/truth/truth.bed ~/data/hg/19/hg19.fa /root/data/sim_deletions_0_${count}/truth; 
	fi
done;

# WGS simulate half the # of pe reads found in the original wgs/all.stat for chromA and chromB
for c in chr2 chr3 chr4;
do
	echo "Simulating WGS reads for ${c}...";
	if [ ! -f /root/data/sim_deletions_0_${count}/wgs/${c}.chromA.predat.gz ];
	then 
	bash simulation/wgsim.sh /root/data/NA12878_HG19_GS_0_1/wgs/all.stat /root/data/sim_deletions_0_${count}/truth/${c}.chromA.fa /root/data/sim_deletions_0_${count}/wgs -1;
	fi

	if [ ! -f /root/data/sim_deletions_0_${count}/wgs/${c}.chromB.predat.gz ];
	then 
	bash simulation/wgsim.sh /root/data/NA12878_HG19_GS_0_1/wgs/all.stat /root/data/sim_deletions_0_${count}/truth/${c}.chromB.fa /root/data/sim_deletions_0_${count}/wgs $(gzip -c -d /root/data/sim_deletions_0_${count}/wgs/${c}.chromA.predat.gz | tail -1 | cut -f 1);
	fi
done;

echo "Converting WGS reads to Ref...";
# CONVERT TO ORIGINAL REFERENCE COORDINATES FIRST!
if [ ! -f /root/data/sim_deletions_0_${count}/wgs/${c}.chromA.ref.predat.gz ]
then
parallel -j ${threads}	bash simulation/convert_predat_to_ref.sh {1} wgs {2} ${count} ::: chr2 chr3 chr4 ::: chromA chromB;
fi

echo "Converting WGS reads to all.dat...";
# merge chromA and chromB --> chrom.all.dat
parallel -j ${threads} bash simulation/merge_to_all_dat.sh {1} wgs ${count} ::: chr2 chr3 chr4;

echo "Filtering WGS reads with overlap to VCF...";
# filter chromA with vcf --> chrom.A.dat
# filter chromB with vcf --> chrom.B.dat
parallel -j ${threads} bash simulation/filter_predat_by_vcf.sh {1} wgs {2} ${count} ::: chr2 chr3 chr4 ::: chromA chromB;

parallel -j ${threads} rm ::: $(find /root/data/sim_deletions_0_${count}/wgs -name "*.predat.gz")
parallel -j ${threads} rm ::: $(find /root/data/sim_deletions_0_${count}/wgs -name "*.out")

# HiC shuffle half the # of chrom mapped se reads found in the original hic/*.dat count_dat --se for chromA and chromB
# HiC shuffle half the # of same chrom mapped pe reads found in the original hic/*.dat count_dat --pe for chromA and chromB
# reassign se.chromA, pe.chromA, se.chromB,  pe.chromB, index 
for c in chr2 chr3 chr4;
do
	tmpdir=$(mktemp -d -t inphadel.XXXX --tmpdir=/tmp/);
	echo "Simulating HiC reads for ${c}...";

	if [ ! -f /root/data/sim_deletions_0_${count}/hic/${c}.chromA.predat.gz ];
	then 
	mkfifo ${tmpdir}/fifo1;

	cat ${tmpdir}/fifo1 | gzip -c >/root/data/sim_deletions_0_${count}/hic/${c}.chromA.predat.gz &
	python -m svphase.simulate.hic -r 100 -f 0.5 --debug ${c} /root/data/sim_deletions_0_${count}/truth/truth.A.bed /home/anand/data/hg/19/hg19.fa /root/data/sim_deletions_0_${count}/truth/${c}.chromA.fa /root/data/NA12878_HG19_GS_0_1/hic/${c}.all.dat ${tmpdir}/fifo1;
	wait
	rm ${tmpdir}/fifo1;
	fi


	if [ ! -f /root/data/sim_deletions_0_${count}/hic/${c}.chromB.predat.gz ];
	then 
	mkfifo ${tmpdir}/fifo2;
	awk -v offset=$(gzip -c -d /root/data/sim_deletions_0_${count}/hic/${c}.chromA.predat.gz | tail -1 | cut -f 1) '{print $1+1+offset "\t" $2 "\t" $3 "\t" $4 "\t" $5;}' ${tmpdir}/fifo2 | gzip -c >/root/data/sim_deletions_0_${count}/hic/${c}.chromB.predat.gz &
	python -m svphase.simulate.hic -r 100 -f 0.5 --debug ${c} /root/data/sim_deletions_0_${count}/truth/truth.B.bed /home/anand/data/hg/19/hg19.fa /root/data/sim_deletions_0_${count}/truth/${c}.chromB.fa /root/data/NA12878_HG19_GS_0_1/hic/${c}.all.dat ${tmpdir}/fifo2;
	wait
	rm ${tmpdir}/fifo2;
	fi

rmdir ${tmpdir};
done;

echo "Converting WGS reads to Ref...";
# CONVERT TO ORIGINAL REFERENCE COORDINATES FIRST!
if [ ! -f /root/data/sim_deletions_0_${count}/hic/${c}.chromA.ref.predat.gz ]
then
parallel -j ${threads}	bash simulation/convert_predat_to_ref.sh {1} hic {2} ${count} ::: chr2 chr3 chr4 ::: chromA chromB;
fi

echo "Converting WGS reads to all.dat...";
# merge chromA and chromB --> chrom.all.dat
parallel -j ${threads} bash simulation/merge_to_all_dat.sh {1} hic ${count} ::: chr2 chr3 chr4;

echo "Filtering WGS reads with overlap to VCF...";
# filter chromA with vcf --> chrom.A.dat
# filter chromB with vcf --> chrom.B.dat
parallel -j ${threads} bash simulation/filter_predat_by_vcf.sh {1} hic {2} ${count} ::: chr2 chr3 chr4 ::: chromA chromB;

parallel -j ${threads} rm ::: $(find /root/data/sim_deletions_0_${count}/hic -name "*.predat.gz")
parallel -j ${threads} rm ::: $(find /root/data/sim_deletions_0_${count}/truth -name "*.fa")
parallel -j ${threads} rm ::: $(find /root/data/sim_deletions_0_${count}/truth -name "*.fa.fai")
