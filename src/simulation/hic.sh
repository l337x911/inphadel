count=$1
c=$2

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

