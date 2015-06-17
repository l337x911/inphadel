
svfpath=$1;
dat=${2};
vcffpath=${3};
letter=${4};
down=2;
#down2=2;

echo ${svfpath};
echo ${dat};
echo ${vcffpath};
py=python;
#py=~/bin/python;


#for letter in A B;
#do
  #$py -m svphase.simulate.convert_dat_to_ref ${svfpath} ${dat}.${letter}.dat;
#$py -m svphase.scripts.downsample_dat ${dat}.${letter}.ref.dat ${down};
$py -m svphase.scripts.filter_dat ${dat}.${letter}.hic.ref.dat ${vcffpath};

# for homozygous callsets
#$py -m svphase.scripts.downsample_dat ${dat}.${letter}.ref.dat ${down2};
#$py -m svphase.scripts.het_split_dat ${dat}.${letter}.ref.down${down2}.dat ${dat}.${letter}1.ref.down${down2}.dat ${dat}.${letter}2.ref.down${down2}.dat
#$py -m svphase.scripts.filter_dat ${dat}.${letter}1.ref.down${down2}.dat ${vcffpath};
#$py -m svphase.scripts.filter_dat ${dat}.${letter}2.ref.down${down2}.dat ${vcffpath};
#sleep 2;
#done;


