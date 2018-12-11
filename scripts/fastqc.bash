#!/bin/bash
#$ -S /bin/bash
#$ -cwd


module load fastqc/0.11.4

# change directory into the target folder
#cd $1

#for i in *.fastq.gz
#do
#mkdir ../qc/${i}
#fastqc -o ../qc/${i} ./${i}
#done

fastqc ./ JL01_S1_R1_001.fastq.gz