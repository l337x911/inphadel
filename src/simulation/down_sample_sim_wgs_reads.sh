
svfpath=$1;
fq=${2};
vcffpath=${3};
down=8;
down2=4;

echo ${svfpath};
echo ${fq};
echo ${vcffpath};
#py=python;
py=~/bin/python;

for letter in A B;
do
  #$py -m svphase.simulate.fq_to_dat ${fq}.${letter}.wgs_1.fq;
  #$py -m svphase.simulate.convert_dat_to_ref ${svfpath} ${fq}.${letter}.wgs.dat;
  $py -m svphase.scripts.downsample_dat ${fq}.${letter}.wgs.ref.dat ${down};
  $py -m svphase.scripts.filter_dat ${fq}.${letter}.wgs.ref.down${down}.dat ${vcffpath};

  # for homozygous callsets
  $py -m svphase.scripts.downsample_dat ${fq}.${letter}.wgs.ref.dat ${down2};
  $py -m svphase.scripts.het_split_dat ${fq}.${letter}.wgs.ref.down${down2}.dat ${fq}.${letter}1.wgs.ref.down${down2}.dat ${fq}.${letter}2.wgs.ref.down${down2}.dat
  $py -m svphase.scripts.filter_dat ${fq}.${letter}1.wgs.ref.down${down2}.dat ${vcffpath};
  $py -m svphase.scripts.filter_dat ${fq}.${letter}2.wgs.ref.down${down2}.dat ${vcffpath};
  sleep 2;
done;

# for whole unphased datasets
$py -m svphase.scripts.merge_dat ${fq}.A.wgs.ref.down${down}.dat ${fq}.B.wgs.ref.down${down}.dat ${fq}.wgs.ref.down${down}.dat;
$py -m svphase.scripts.merge_dat ${fq}.A.wgs.ref.down${down2}.dat ${fq}.B.wgs.ref.down${down2}.dat ${fq}.wgs.ref.down${down2}.dat;

