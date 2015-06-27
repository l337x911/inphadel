
cd /root/data/sim_deletions_0_000/truth;
for c in chr2 chr3 chr4; 
do 
bash ~/Projects/inphadel-min/src/simulation/shuffle_deletions.sh ${c} ~/Projects/inphadel-min/data truth.${c}.bed; 
done;

cat truth.chr?.bed >truth.bed;

#Verify bedsize distribution
bed=truth.bed && awk '{print $3-$2;}' ${bed} | sort -n -k 1,1 | awk -v wc="$(wc -l ${bed})" '{print $0, NR/wc;}' >${bed}.cumsum
awk 'BEGIN {x=0;} {if($3<0){print x, $0;}}' truth.bed >truth.bed.overlap

cd /home/anand/Projects/inphadel-min/src
# Create chrom A fasta sequence with hom and pA deletions 
# Create chrom B fasta sequence with hom and pB deletions
for c in chr2 chr3 chr4; do python -m svphase.scripts.create_deletion_fa_from_bed ${c} /root/data/sim_deletions_0_000/truth/truth.bed ~/data/hg/19/hg19.fa /root/data/sim_deletions_0_000/truth; done;

# WGS simulate half the # of pe reads found in the original wgs/all.stat for chromA and chromB
for c in chr2 chr3 chr4;
do
	simulation/wgsim.sh /root/data/NA12878_HG19_GS_0_1/wgs/all.stat /root/data/sim_deletions_0_000/truth/${c}.chromA.fa /root/data/sim_deletions_0_000/wgs -1;
	simulation/wgsim.sh /root/data/NA12878_HG19_GS_0_1/wgs/all.stat /root/data/sim_deletions_0_000/truth/${c}.chrom${h}.fa /root/data/sim_deletions_0_000/wgs $(gzip -c -d /root/data/sim_deletions_0_000/wgs/${c}.chromA.predat.gz | tail -1 | cut -f 1);
done;

# CONVERT TO ORIGINAL REFERENCE COORDINATES FIRST!
python -m svphase.simulate.convert----


# merge chromA and chromB --> chrom.all.dat
for c in chr2 chr3 chr4;
do
cat <(gzip -c -d /root/data/sim_deletions_0_000/wgs/${c}.chromA.predat.gz) <(gzip -c -d /root/data/sim_deletions_0_000/wgs/${c}.chromB.predat.gz) | tee /root/data/sim_deletions_0_000/wgs/${c}.all.predat | LC_ALL=C sort -n -S 10G -k 2,2 | python -m svphase.scripts.predat_to_dat  >/root/data/sim_deletions_0_000/wgs/${c}.all.dat
done;

# filter chromA with vcf --> chrom.A.dat
# filter chromB with vcf --> chrom.B.dat
for c in chr2 chr3 chr4;
do
gzip -c -d /root/data/sim_deletions_0_000/wgs/${c}.chromA.predat.gz | python -m svphase.simulate.filter_predat_by_vcf /root/data/sim_deletions_0_000/vcf/${c}.vcf | LC_ALL=C sort -n -S 10G -k 2,2 | python -m svphase.scripts.predat_to_dat >/root/data/sim_deletions_0_000/wgs/${c}.A.dat
gzip -c -d /root/data/sim_deletions_0_000/wgs/${c}.chromB.predat.gz | python -m svphase.simulate.filter_predat_by_vcf /root/data/sim_deletions_0_000/vcf/${c}.vcf | LC_ALL=C sort -n -S 10G -k 2,2 | python -m svphase.scripts.predat_to_dat >/root/data/sim_deletions_0_000/wgs/${c}.B.dat
done;


# HiC shuffle half the # of same chrom mapped pe reads found in the original hic/*.dat count_dat --pe for chromA and chromB
# HiC shuffle half the # of chrom mapped se reads found in the original hic/*.dat count_dat --se for chromA and chromB
# reassign pe.chromB, se.chromA, se.chromB index 

# merge chromA and chromB --> chrom.all.dat
# filter chromA with vcf --> chrom.A.dat
# filter chromB with vcf --> chrom.B.dat

