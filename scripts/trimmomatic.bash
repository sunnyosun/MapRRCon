#!/bin/bash

#SBATCH --job-name=trimmomatic
##SBATCH --nodes=1
##SBATCH --cpus-per-task=1
#SBATCH --mem=100GB
##SBATCH --gres=gpu:1
##SBATCH --partition=gpu4_medium
#SBATCH --partition=cpu_long
#SBATCH --error=/gpfs/data/proteomics/projects/Sunny//err_out/%x_%j.err
#SBATCH --output=/gpfs/data/proteomics/projects/Sunny/err_out/%x_%j.out
##SBATCH --dependency=afterany:job_id

trimmomaticpath=/gpfs/data/proteomics/home/xs338/Trimmomatic-0.36

infq1=$1
infq2=$2
adapterfa=/gpfs/data/proteomics/home/xs338/Trimmomatic-0.36/adapters/TruSeq-adapters.fa
outfq1=$infq1.cleaned.fastq
outfq2=$infq2.cleaned.fastq
outrmfq1=$infq1.removed.fastq
outrmfq2=$infq2.removed.fastq
logfile=$infq1.log
echo "infq1=$infq1"
echo "infq2=$infq2"
echo "adapterfa=$adapterfa"

java -jar ${trimmomaticpath}/trimmomatic-0.36.jar PE -phred33 -trimlog $logfile $infq1 $infq2 $outfq1 $outrmfq1 $outfq2 $outrmfq2 ILLUMINACLIP:${adapterfa}:3:30:7:1:TRUE LEADING:2 TRAILING:2 SLIDINGWINDOW:4:10 MINLEN:36

# if single end, using the following
# java -jar ${trimmomaticpath}/trimmomatic-0.36.jar SE -phred33 -trimlog $logfile $infq1 $outfq1 ILLUMINACLIP:${adapterfa}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

gzip $outfq1
gzip $outfq2
gzip $outrmfq1
gzip $outrmfq2
