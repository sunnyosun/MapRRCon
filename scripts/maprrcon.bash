#!/bin/bash                                                                                                                     
#$ -S /bin/bash                                                                                                                 
#$ -cwd                                                                                                                         

module load bwa/0.7.7
module load samtools/1.3
module load bedtools/2.26.0
module load java/1.8
module load jre/1.8

# needs 4 inputs in the commandline                                                                                             
# 1. read.fastq.gz                                                                                                              
# 2. destination folder (subfolder of Sunny)                                                                                    
# 3. number of mismatches


fq_r1=$1_R1_001.fastq.gz
fq_r2=$1_R2_001.fastq.gz
bam=$1.bam
sortbam=$1.sort
coverage=$1_coverage.bedgraph
l1hsbam=$sortbam.bam.L1HsOnly.bam

genome=/ifs/data/proteomics/projects/Sunny/genome/hg38.fa
l1hs=/ifs/data/proteomics/projects/Sunny/genome/L1HS.fa
dir=/ifs/data/proteomics/projects/Sunny/$2
dir2=/ifs/data/proteomics/projects/Sunny/te_extraction
pathtoreg=/ifs/data/proteomics/projects/Sunny/RepeatMaskerL1HsInfo
reg=ucsc.rmsk.hg38.L1HS.bed
threshmis=$3

echo "fq_r1=${fq_r1}"
echo "sortbam=$sortbam"

cd $dir
bwa mem -t 4 $genome ${dir}/${fq_r1} $dir/${fq_r2} | samtools view -b -S - > ${dir}/${bam}
samtools sort $dir/${bam} > $dir/${sortbam}.bam
samtools index $dir/${sortbam}.bam
rm ${dir}/${bam}
bedtools genomecov -ibam $dir/${sortbam}.bam -bga > $dir/${coverage}

#L1HS feature extraction
cd $dir2
java -classpath sam-1.89.jar:. ExtractSubsetOfBAMBasedOnBedFile $dir $pathtoreg $sortbam.bam $reg $threshmis

cd $dir
# convert bam to fastq
samtools fastq $dir/$l1hsbam > $dir/$l1hsbam.fastq

# map to L1HS
bwa mem -t 1 ${l1hs} ${dir}/${l1hsbam}.fastq | samtools view -b -S - > ${dir}/${l1hsbam}.fastq.bam
samtools sort ${dir}/${l1hsbam}.fastq.bam > ${dir}/${l1hsbam}.fastq.bam.sort.bam
rm ${dir}/${l1hsbam}.fastq.bam
samtools index ${dir}/${l1hsbam}.fastq.bam.sort.bam
bedtools genomecov -ibam ${dir}/${l1hsbam}.fastq.bam.sort.bam -d > ${dir}/$1_L1HS_coverage.bedgraph
