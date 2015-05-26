
prefix=$1;
start=$2;
end=$3;
#letter=$4;

for i in $(seq ${start} 1 ${end});
do
  python -m svphase.scripts.shuffle_bed_split ${prefix}.${i}.bed;
done;


for letter in A B;
do
  for i in $(seq ${start} 1 ${end});
  do
    echo ${prefix}.${i}.${letter}.bed;
    python -m svphase.chromosome ${prefix}.${i}.${letter}.bed >${prefix}.${i}.${letter}.hg18.hind3.bed;
    bash resolute10_cutsites.sh ${prefix}.${i}.${letter}.hg18.hind3.bed 
    bash resolute100_cutsites.sh ${prefix}.${i}.${letter}.hg18.hind3.bed 
  done;
done;
