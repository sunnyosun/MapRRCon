#!/bin/bash

#$ -N trimmomatic
#$ -S /bin/bash
#$ -cwd
#$ -l h_rt=72:00:00
#$ -o /ifs/data/proteomics/projects/Sunny/jingchuan/2018-07-24_syni/fastq/err_out/$JOB_NAME_$JOB_ID.out
#$ -e /ifs/data/proteomics/projects/Sunny/jingchuan/2018-07-24_syni/fastq/err_out/$JOB_NAME_$JOB_ID.err

trimmomaticpath=/ifs/home/xs338/Trimmomatic-0.36

infq1=$1
#infq2=$2
adapterfa=/ifs/home/xs338/Trimmomatic-0.36/adapters/TruSeq-adapters.fa
outfq1=$infq1.cleaned.fastq
#outfq2=$infq2.cleaned.fastq
outrmfq1=$infq1.removed.fastq
#outrmfq2=$infq2.removed.fastq
logfile=$infq1.log
echo "infq1=$infq1"
#echo "infq2=$infq2"
echo "adapterfa=$adapterfa"

#java -jar ${trimmomaticpath}/trimmomatic-0.36.jar PE -threads $NSLOTS -phred33 -trimlog $logfile $infq1 $infq2 $outfq1 $outrmfq1 $outfq2 $outrmfq2 ILLUMINACLIP:${adapterfa}:3:30:7:1:TRUE LEADING:2 TRAILING:2 SLIDINGWINDOW:4:10 MINLEN:36

java -jar ${trimmomaticpath}/trimmomatic-0.36.jar SE -phred33 -trimlog $logfile $infq1 $outfq1 ILLUMINACLIP:${adapterfa}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

