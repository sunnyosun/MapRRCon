#!/bin/bash                                                                                                                     

#SBATCH --job-name=MapRRCon
##SBATCH --nodes=1
##SBATCH --cpus-per-task=1
#SBATCH --mem=100GB
##SBATCH --gres=gpu:1
##SBATCH --partition=gpu4_medium
#SBATCH --partition=cpu_short
#SBATCH --error=/gpfs/data/proteomics/projects/Sunny/chipseq/20181203_FCHN3FWBGX7/err_out/%x_%j.err
#SBATCH --output=/gpfs/data/proteomics/projects/Sunny/chipseq/20181203_FCHN3FWBGX7/err_out/%x_%j.out
##SBATCH --dependency=afterany:job_id

####################################################################
####################################################################
########################### Modules ################################
####################################################################
####################################################################

module load bwa/0.7.17
module load samtools/1.3
module load bedtools/2.26.0


####################################################################
########################### Inputs #################################
####################################################################

# needs 2 inputs from commandline                                                                                             
# 1. folder that contains the fastq.gz files (full path)                            
# 2. L1 consensus (e.g. L1HS)
# 3. number of mismatches

dir=$1
cd ${dir}
filename=(`ls *.fastq.gz`)
str=${filename[0]}
IFS=’_’ read -ra NAMES <<< "$str"
name=${NAMES[0]}
# make a maprrcon directory
mkdir maprrcon_$2
bam=${name}.bam
sortbam=${name}.sort
coverage=${name}_coverage.bedgraph
l1hsbam=${sortbam}.bam.L1HsOnly.bam

# directories of reference genome and L1
genome=/gpfs/data/proteomics/projects/Sunny/genome/human/hg38.fa
l1hs=/gpfs/data/proteomics/projects/Sunny/genome/L1/L1HS.fa
dir2=/gpfs/data/proteomics/projects/Sunny/te_extraction
pathtoreg=/gpfs/data/proteomics/projects/Sunny/RepeatMaskerL1HsInfo
reg=ucsc.rmsk.hg38.L1HS.bed
threshmis=$3

# get the sample names and determine if paired end
alen=${#filename[@]}
if (($alen == 2));then
    echo "pair"
    fq_r1=${filename[0]}
    fq_r2=${filename[1]}
    echo "fq_r1=${fq_r1}"
    echo "fq_r2=${fq_r2}"
    cd ${dir}
    bwa mem -t 4 ${genome} ${dir}/${fq_r1} ${dir}/${fq_r2} | samtools view -b -S - > ${dir}/${bam}
else
    echo "single"
    fq_r1=${filename[0]}
    echo "fq_r1=${fq_r1}"
    cd ${dir}
    bwa mem -t 4 ${genome} ${dir}/${fq_r1} | samtools view -b -S - > ${dir}/${bam}
fi

# alignment to the reference genome
samtools sort ${dir}/${bam} > ${dir}/${sortbam}.bam
samtools index ${dir}/${sortbam}.bam
rm ${dir}/${bam}
bedtools genomecov -ibam ${dir}/${sortbam}.bam -bga > ${dir}/${coverage}

# L1HS feature extraction
cd ${dir2}
java -classpath sam-1.89.jar:. ExtractSubsetOfBAMBasedOnBedFile ${dir} ${pathtoreg} ${sortbam}.bam ${reg} ${threshmis}

cd ${dir}
# convert bam to fastq
samtools fastq ${dir}/${l1hsbam} > ${dir}/${l1hsbam}.fastq

# remap to L1HS
bwa mem -t 1 ${l1hs} ${dir}/${l1hsbam}.fastq | samtools view -b -S - > ${dir}/${l1hsbam}.fastq.bam
samtools sort ${dir}/${l1hsbam}.fastq.bam > ${dir}/${l1hsbam}.fastq.bam.sort.bam
rm ${dir}/${l1hsbam}.fastq.bam
samtools index ${dir}/${l1hsbam}.fastq.bam.sort.bam
bedtools genomecov -ibam ${dir}/${l1hsbam}.fastq.bam.sort.bam -d > ${dir}/${name}_L1HS_coverage.bedgraph
