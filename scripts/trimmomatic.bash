#!/bin/bash
#$ -S /bin/bash
#$ -cwd


trimmomaticpath=/ifs/home/xs338/Trimmomatic-0.36

infq1=$1
infq2=$2
adapterfa=$3
outfq1=$infq1.cleaned.fastq
outfq2=$infq2.cleaned.fastq
outrmfq1=$infq1.removed.fastq
outrmfq2=$infq2.removed.fastq
logfile=$infq1.log
echo "infq1=$infq1"
echo "infq2=$infq2"
echo "adapterfa=$adapterfa"

java -jar ${trimmomaticpath}/trimmomatic-0.36.jar PE -threads $NSLOTS -phred33 -trimlog $logfile $infq1 $infq2 $outfq1 $outrmfq1 $outfq2 $outrmfq2 ILLUMINACLIP:${adapterfa}:3:30:7:1:TRUE LEADING:2 TRAILING:2 SLIDINGWINDOW:4:10 MINLEN:36