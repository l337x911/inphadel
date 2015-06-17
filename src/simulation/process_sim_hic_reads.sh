
svfpath=$1;
dat=${2};
vcffpath=${3};
letter=${4};

sleep $(date +%N | sed -e 's/000$//' -e 's/^0//' | tail -c 3)

echo ${svfpath};
echo ${dat};
echo ${vcffpath};
py=python;
#py=~/bin/python;

#$py -m svphase.simulate.convert_dat_to_ref ${svfpath}.${letter}.bed ${dat}.${letter}.hic.dat;
$py -m svphase.scripts.downsample_dat ${dat}.${letter}.hic.ref.dat 2;
$py -m svphase.scripts.filter_dat ${dat}.${letter}.hic.ref.down2.dat ${vcffpath};

# for homozygous callsets
$py -m svphase.scripts.het_split_dat ${dat}.${letter}.hic.ref.dat ${dat}.${letter}1.hic.ref.dat ${dat}.${letter}2.hic.ref.dat
$py -m svphase.scripts.filter_dat ${dat}.${letter}1.hic.ref.dat ${vcffpath};
$py -m svphase.scripts.filter_dat ${dat}.${letter}2.hic.ref.dat ${vcffpath};

