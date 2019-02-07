# MapRRCon step 5 - Normalization
# To normalize ChIP file against Input
# need the mapped reads file

# This function takes 6 varibles:
# 1) chip_file: name of the ChIP file
# 2) ctr_file: name of the INPUT/control file
# 3) mapped_chip: # of mapped reads in the ChIP file
# 4) mapped_control: # of mapped reads in the control file
# 5) input directory
# 6) output directory

####################################################################################################################
#### to run this function: normalize_chip(chip_file,ctr_file,mapped_chip,mapped_control,inputdir,outputdir)
####################################################################################################################
# function
normalize_chip = function(chip_file,ctr_file,mapped_chip,mapped_control,input_dir,output_dir)
{
  # reading input files
  setwd(input_dir)
  ctr_data=read.table(ctr_file,stringsAsFactors = F,sep='\t')
  chip_data=read.table(chip_file,stringsAsFactors = F,sep='\t')
  
  # filtering out positions with low coverage
  ctr_data=ctr_data[which(ctr_data$V3>=10),]
  chip_data=chip_data[which(chip_data$V3>=10),]

  # normalizing library sizes
  ctr_data[,3]=(ctr_data[,3]/mapped_control)*1000000
  chip_data[,3]=(chip_data[,3]/mapped_chip)*1000000

  # normalizing to median
  m1=median(ctr_data[,3])
  m2=median(chip_data[,3])
  ctr_data[,4]=ctr_data[,3]/m1
  chip_data[,4]=chip_data[,3]/m2
  data=merge(chip_data,ctr_data,by='V2')
  data[,8]=data[,4]-data[,7]

  # outputs the normalized signals
  setwd(paste(output_dir,'/norm',sep=''))
  write.table(cbind(data[,1],data[,8]),paste(chip_file,'.txt',sep=''),quote=F,sep='\t',row.names=F,col.names=F)
  
  # plotting
  setwd(paste(output_dir,'/plots',sep=''))
  setEPS()
  postscript(paste(chip_file,'.eps',sep=''),height=6,width=4,paper='special')
  par(mfrow=c(3,1),lab=c(x=5,y=5,len=1),las=1,tcl=-0.25,mgp=c(1.75,0.75,0.3),
      bty='n',mar=c(2,4,1,1),cex.lab=1.5,cex.axis=2,font.lab=1,cex=0.5,
      omi=c(0,0,0,0))
  plot(chip_data[,2],chip_data[,3],col='dodgerblue',lwd=2,cex=2,type='l',main=chip_file,xlab='',ylab='')
  plot(ctr_data[,2],ctr_data[,3],col='dodgerblue',lwd=2,cex=2,type='l',main=ctr_file,xlab='',ylab='')
  plot(data[,1],data[,8],col='dodgerblue',lwd=2,cex=2,type='l',main=chip_file,xlab='',ylab='')
  dev.off()
}

