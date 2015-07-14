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
	if [ ! -f /root/data/sim_deletions_0_${count}/wgs/${c}.A.dat ];
	then 
	bash simulation/wgsim.sh /root/data/NA12878_HG19_GS_0_1/wgs/all.stat /root/data/sim_deletions_0_${count}/truth/${c}.chromA.fa /root/data/sim_deletions_0_${count}/wgs -1;
	fi

	if [ ! -f /root/data/sim_deletions_0_${count}/wgs/${c}.B.dat ];
	then 
	bash simulation/wgsim.sh /root/data/NA12878_HG19_GS_0_1/wgs/all.stat /root/data/sim_deletions_0_${count}/truth/${c}.chromB.fa /root/data/sim_deletions_0_${count}/wgs $(gzip -c -d /root/data/sim_deletions_0_${count}/wgs/${c}.chromA.predat.gz | tail -1 | cut -f 1);
	fi
done;

echo "Converting WGS reads to Ref...";
# CONVERT TO ORIGINAL REFERENCE COORDINATES FIRST!
parallel -j ${threads}	bash simulation/convert_predat_to_ref.sh {1} wgs {2} ${count} ::: chr2 chr3 chr4 ::: chromA chromB;

echo "Converting WGS reads to all.dat...";
# merge chromA and chromB --> chrom.all.dat

parallel -j ${threads} bash simulation/merge_to_all_dat.sh {1} wgs ${count} ::: chr2 chr3 chr4;

echo "Filtering WGS reads with overlap to VCF...";
# filter chromA with vcf --> chrom.A.dat
# filter chromB with vcf --> chrom.B.dat
parallel -j ${threads} bash simulation/filter_predat_by_vcf.sh {1} wgs {2} ${count} ::: chr2 chr3 chr4 ::: chromA chromB;

#parallel -j ${threads} rm ::: $(find /root/data/sim_deletions_0_${count}/wgs -name "*.predat.gz")
#parallel -j ${threads} rm ::: $(find /root/data/sim_deletions_0_${count}/wgs -name "*.out")

# HiC shuffle half the # of chrom mapped se reads found in the original hic/*.dat count_dat --se for chromA and chromB
# HiC shuffle half the # of same chrom mapped pe reads found in the original hic/*.dat count_dat --pe for chromA and chromB
# reassign se.chromA, pe.chromA, se.chromB,  pe.chromB, index 
bash simulation/hic.sh ${count} chr2 
bash simulation/hic.sh ${count} chr3 
bash simulation/hic.sh ${count} chr4

echo "Converting HiC reads to Ref...";
# CONVERT TO ORIGINAL REFERENCE COORDINATES FIRST!
parallel -j ${threads} bash simulation/convert_predat_to_ref.sh {1} hic {2} ${count} ::: chr2 chr3 chr4 ::: chromA chromB;

echo "Converting HiC reads to all.dat...";
# merge chromA and chromB --> chrom.all.dat
parallel -j ${threads} bash simulation/merge_to_all_dat.sh {1} hic ${count} ::: chr2 chr3 chr4;

echo "Filtering HiC reads with overlap to VCF...";
# filter chromA with vcf --> chrom.A.dat
# filter chromB with vcf --> chrom.B.dat
parallel -j ${threads} bash simulation/filter_predat_by_vcf.sh {1} hic {2} ${count} ::: chr2 chr3 chr4 ::: chromA chromB;

python -m svphase.simulate.mkstats /root/data/NA12878_HG19_GS_0_1/ /root/data/sim_deletions_0_${count}

parallel -j ${threads} rm ::: $(find /root/data/sim_deletions_0_${count}/hic -name "*.predat.gz")
parallel -j ${threads} rm ::: $(find /root/data/sim_deletions_0_${count}/truth -name "*.fa")
parallel -j ${threads} rm ::: $(find /root/data/sim_deletions_0_${count}/truth -name "*.fa.fai")
