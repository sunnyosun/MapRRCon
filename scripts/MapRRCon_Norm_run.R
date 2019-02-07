# Script to run normalization
source('MapRRCon_Norm.R')
args = commandArgs(trailingOnly=TRUE)
args[1]->readme_file
args[2]->stats_file
args[3]->l1family
args[4]->inputdir
args[5]->outputdir
readme=read.table(readme_file,sep='\t',header=T,stringsAsFactors=F)
stats=read.table(stats_file,sep='\t',header=F,stringsAsFactors=F)

setwd(inputdir)

for (i in 1:nrow(readme)){
    chip_file=paste(readme[i,2],'_',l1family,'_coverage.bedgraph',sep='')
    ctr_file=paste(readme[i,3],'_',l1family,'_coverage.bedgraph',sep='')
    mapped_chip=stats[which(stats[,1]==readme[i,2]),2]
    mapped_control=stats[which(stats[,1]==readme[i,3]),2]
    normalize_chip(chip_file,ctr_file,mapped_chip,mapped_control,inputdir,outputdir)
}

