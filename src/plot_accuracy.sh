#prefixfpath=/home/anand/Projects/assembly/data/simulations/sim.del-chr19-20.5
#prefixfpath=/home/anand/Projects/assembly/data/simulations/models/sim.del-chr19-20.5norm
#prefixfpath=/home/anand/Projects/assembly/data/simulations/models/sim.del-chr19-20-hm.5-2c
prefixfpath=/home/anand/Projects/assembly/data/simulations/models/sim.del-chr19-20-hm.10-2c

python -m svphase.learn.break_error ${prefixfpath}.10.RandomForest.pkl None 7;

