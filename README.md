
=========== InPhaDel: (In)tegrating whole genome and proximity ligation sequencing to (Pha)se (Del)etion variants. ===========

Developed by Anand D. Patel in Vineet Bafna's lab at University of California,
San Diego patel [dot] anandd [at] gmail [dot] com

InPhaDel takes as input, a scaffold of phased single nucleotide variants (vcf),
whole genome sequencing data (sorted bam), proximity ligation data (sorted
bam), and deletions (bed) and a classification model ('rf|svm|knn').

InPhaDel returns a prediction for the deletion being homozygous, scaffold A,
scaffold B, or unlikely to be a deletion (bed).  For details of the method, see
Patel, Selvaraj, Bafna 2015.

Installing 
=========

Software Requirements 
-------------

InPhaDel requires Python 2.7.x with packages numpy 1.9, pandas 1.13,
scikit-learn 2.7.6, matplotlib [optional].  This tool may work on other version
of pandas, scikit-learn, and numpy, but has not been tested.

The tool also requires samtools 

Installing InPhaDel 
------------- 

InPhaDel can be run from source or installed.

Download the tar or zip ball and install using the setup.py script.
	
	python setup.py install

To verify the installation/configuration is correct and run a small test case
for InPhaDel (takes <5min) ::

	python setup.py test 
	
If all is well, continue to usage.

Alternatively, to run from source directory use ::

	python -m svphase.test.predict

Usage 
===========

InPhaDel requires inputs to follow a specific naming convention as described
below,

 a) Deletions must be in a simple bed file. Chromosome names must match
reference, and below input bam filenames.

 b) Reference in fasta file format (with fai index). Use `samtools faidx` to
generate index.

 c) Phased scaffolds in VCF file format, WGS, and HiC bam files split by
chromosome. Some are required, and optional files will be generated. 

  - INPUT_DIR/wgs/all.stat - generated by `samtools idxstat` of composite wgs
    data [required]
  - INPUT_DIR/hic/all.stat - generated by `samtools idxstat` of composite hic
    data [required]
  - INPUT_DIR/wgs/CHROMOSOME.all.bam[.bai] - WGS data where .bai is generated
    by `samtools index` [required]
  - INPUT_DIR/hic/CHROMOSOME.all.bam[.bai] - HiC data where .bai is generated
    by `samtools index` [required]
  - INPUT_DIR/vcf/CHROMOSOME.vcf [optional if INPUT_DIR/wgs/CHROMOSOME.[AB].bam
    and INPUT_DIR/hic/CHROMOSOME.[AB].bam] have been generated.
  
INPUT_DIR is an argument to `inphadel` and CHROMOSOME is the name of a
chromosome in the reference fasta file. There must be CHROMOSOME bam files for
each CHROMOSOME in the bed file.

To see command line arguments, use ::

	inphadel -h

or from source ::

	python -m svphase.inphadel -h


