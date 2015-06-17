
svfpath=$1;
#fq=${2:0:${#2}-11};
fq=${2};
vcffpath=${3};
letter=${4};

echo ${svfpath};
echo ${fq};
echo ${vcffpath};
py=python;
#py=~/bin/python;

#$py -m svphase.simulate.fq_to_dat ${fq}.${letter}.wgs_1.fq;
#$py -m svphase.simulate.convert_dat_to_ref ${svfpath}.${letter}.bed ${fq}.${letter}.wgs.dat;
$py -m svphase.scripts.downsample_dat ${fq}.${letter}.wgs.ref.dat 4;
$py -m svphase.scripts.filter_dat ${fq}.${letter}.wgs.ref.down4.dat ${vcffpath};
#  gzip ${fq}.${letter}.wgs_1.fq &

# for homozygous callsets
$py -m svphase.scripts.downsample_dat ${fq}.${letter}.wgs.ref.dat 2;
$py -m svphase.scripts.het_split_dat ${fq}.${letter}.wgs.ref.down2.dat ${fq}.${letter}1.wgs.ref.down2.dat ${fq}.${letter}2.wgs.ref.down2.dat
$py -m svphase.scripts.filter_dat ${fq}.${letter}1.wgs.ref.down2.dat ${vcffpath};
$py -m svphase.scripts.filter_dat ${fq}.${letter}2.wgs.ref.down2.dat ${vcffpath};
#done;


$py -m svphase.scripts.merge_dats ${fq}.A.wgs.ref.down4.dat ${fq}.B.wgs.ref.down4.dat ${fq}.wgs.ref.down4.dat;
$py -m svphase.scripts.merge_dats ${fq}.A.wgs.ref.down2.dat ${fq}.B.wgs.ref.down2.dat ${fq}.wgs.ref.down2.dat;

