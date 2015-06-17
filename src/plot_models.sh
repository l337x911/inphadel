#prefixfpath=/home/anand/Projects/assembly/data/simulations/models/sim.del-chr18-20.5
#prefixfpath=/home/anand/Projects/assembly/data/simulations/models/sim.del-chr19-20.5norm
#prefixfpath=/home/anand/Projects/assembly/data/simulations/models/sim.del-chr19-20.5-2c
#prefixfpath=/home/anand/Projects/assembly/data/simulations/models/sim.del-chr19-20.10-2c
#prefixfpath=/home/anand/Projects/assembly/data/simulations/models/sim.del-chr19-20-hm.5-2c
prefixfpath=/home/anand/Projects/assembly/data/simulations/models/sim.del-chr19-20-hm.10-2c
#prefixfpath=/home/anand/Projects/assembly/data/simulations/models/real.hm.5-2c

#for cv in 05 07 10;
for cv in 10;
do
  #python -m svphase.learn.predict ${prefixfpath}.${cv}.simple_sum-KN.pkl;
  #python -m svphase.learn.predict ${prefixfpath}.${cv}.KN.pkl;
  #python -m svphase.learn.predict ${prefixfpath}.${cv}.RandomForest.pkl;
  #python -m svphase.learn.predict ${prefixfpath}.${cv}.SVM.pkl;

  #python -m svphase.learn.predict ${prefixfpath} ${cv}.simple_sum-KN;
  #python -m svphase.learn.predict ${prefixfpath} ${cv}.KN;
  #python -m svphase.learn.predict ${prefixfpath} ${cv}.simple_sum-GNB;
  python -m svphase.learn.predict ${prefixfpath} ${cv}.RandomForest;
  python -m svphase.learn.predict ${prefixfpath} ${cv}.SVM;

  #python -m svphase.learn.predict ${prefixfpath}.hic_only.${cv}.RandomForest.pkl ${prefixfpath}.wgs_only.${cv}.RandomForest.pkl;
  #python -m svphase.learn.predict ${prefixfpath}.hic_only.${cv}.SVM.pkl ${prefixfpath}.wgs_only.${cv}.SVM.pkl;
  #python -m svphase.learn.predict ${prefixfpath}.hic_only.${cv}.KN.pkl ${prefixfpath}.wgs_only.${cv}.KN.pkl;
done;
