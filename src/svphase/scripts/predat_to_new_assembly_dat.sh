predat=$1
parentdir=$(dirname ${predat})
chrom=$(basename ${predat} | sed 's/\..*//')
assembly_chain="/home/anand/data/hg/hg18ToHg19.over.chain"

echo "Working on ${predat}"
~/software/liftOver <(awk -v chrom=${chrom} '{print chrom "\t" $2 "\t" $2+1;}' ${predat} ) ${assembly_chain} ${predat}.liftOver ${predat}.unmapped

python -m svphase.scripts.filter_predat_with_unmapped ${predat} ${predat}.liftOver ${predat}.unmapped | sort -n -S 5G -k 2,2b | python -m svphase.scripts.predat_to_dat >${predat:0:${#predat}-7}.dat


