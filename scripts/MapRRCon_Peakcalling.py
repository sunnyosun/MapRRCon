"""
MapRRCon Step 6 - Peak calling

This Script contains the algorithm of peak calling on a consensus sequence
The input is a directory containing normalized coverage file with the following two columns:
1) position
2) normalized signal at each position

The output files are:
1) profile plots with marked peaks that are called
2) a list of peaks containing all the features

This program can be run on the commandline as follows: (user_directory specifies the input directory)

python MapRRCon_Peakcalling.py user_directory

"""

##############################################################################
# import modules
import sys,os
from pandas import Series,DataFrame
import pandas as pd
import numpy as np
import re
from datetime import datetime
import matplotlib
matplotlib.use('agg') 
import matplotlib.pyplot as plt
import openpyxl
from matplotlib.pyplot import cm
from scipy import stats
pd.options.mode.chained_assignment = None  # default='warn'


##############################################################################
plotlist=[15] #[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15] #
dir='Negative_Cases' #'Positive_Cases' #'data' #
# the first input in the commandline
if len(sys.argv)>1:
	dir=sys.argv[1]
outdir='results-'+dir+'-rmsd'
if not os.path.isdir(outdir):
	os.mkdir(outdir)
dfs={}
peaklist={}
rmsd_th=1.3
rmsd_average=1
height_th=1.0
for filename in os.listdir(dir):
	if filename.endswith('.txt') and not filename.startswith('FOXP2'):
		print filename
		count=1
		dfs[filename]=pd.read_table(dir+'/'+filename,header=None)
		dfs[filename].columns=['pos','intensity']
		if count in plotlist:
			print count
			fig,ax = plt.subplots(1,1,figsize=(16,6), sharex=True, sharey=False)
			ax.set_title(filename[:-4],fontsize=20, fontweight='bold')
			ax.set_xlabel('position',fontsize=20, fontweight='bold')
			ax.set_ylabel('intensity',fontsize=20, fontweight='bold')
			ax.plot(dfs[filename]['pos'],dfs[filename]['intensity'],color='darkred',lw=2,alpha=1)
			ax.set_xlim([0,np.max(dfs[filename]['pos'])])
			ax.set_ylim([np.min(dfs[filename]['intensity']),1.02*np.max(dfs[filename]['intensity'])])
			fig.savefig(outdir+'/'+filename[:-4]+'.'+str(count)+'.'+'.png',dpi=72,bbox_inches='tight')
			plt.close(fig)
		count+=1
		for width in [80]: #[40,80,160]:
			count=2
			window='flat'
			w=np.ones(2*width+1,'d')
			dfs[filename]['avg']=np.convolve(w/w.sum(),dfs[filename]['intensity'],'same')	
			dfs[filename]['sd']=(dfs[filename]['intensity']-dfs[filename]['avg'])**2
			dfs[filename]['msd']=np.convolve(w/w.sum(),dfs[filename]['sd'],'same')	
			dfs[filename]['rmsd']=dfs[filename]['msd']**0.5
			if count in plotlist:
				print count
				fig,ax = plt.subplots(1,1,figsize=(16,6), sharex=True, sharey=False)
				ax.set_title(filename[:-4],fontsize=20, fontweight='bold')
				ax.set_xlabel('position',fontsize=20, fontweight='bold')
				ax.set_ylabel('rmsd',fontsize=20, fontweight='bold')
				ax.plot(dfs[filename]['pos'],dfs[filename]['rmsd'],color='black',lw=2,alpha=1)
				ax.set_xlim([0,np.max(dfs[filename]['pos'])])
				ax.set_ylim([0,1.02*np.max(dfs[filename]['rmsd'])])
				fig.savefig(outdir+'/'+filename[:-4]+'.'+str(count)+'.'+'-'+str(width)+'-rmsd.png',dpi=72,bbox_inches='tight')
				plt.close(fig)
			count+=1
						
			fig,ax = plt.subplots(1,1,figsize=(6,6), sharex=True, sharey=False)
			ax.set_title(filename[:-4],fontsize=20, fontweight='bold')
			ax.set_xlabel('Intensity',fontsize=20, fontweight='bold')
			n, bins, patches = ax.hist(dfs[filename]['intensity'],color='black',bins=100,lw=2,alpha=1,histtype='step')
			#ax.set_xlim([0,np.max(bins)])
			if count in plotlist:
				print count
				fig.savefig(outdir+'/'+filename[:-4]+'.'+str(count)+'.'+'-'+str(width)+'-int-hist.png',dpi=72,bbox_inches='tight')
			count+=1
			max_val=np.max(n)
			th=0.2
			th_=0.1
			th_i=0
			max_i=len(n)
			done=0
			for i in range(len(n)):
				if done==0:
					if max_val==n[i]:
						max_i=i
					if i>max_i:
						if max_val*th_>n[i]:
							done=1
					if max_val*th<n[i]:
						th_i=i+1
			ax.plot([bins[th_i],bins[th_i]],[0,max_val],color='darkred',lw=1,alpha=1)
			if count in plotlist:
				print count
				fig.savefig(outdir+'/'+filename[:-4]+'.'+str(count)+'.'+'-'+str(width)+'-int-hist_.png',dpi=72,bbox_inches='tight')
			count+=1
			plt.close(fig)
			intensity_level=bins[th_i]
						
			fig,ax = plt.subplots(1,1,figsize=(6,6), sharex=True, sharey=False)
			ax.set_title(filename[:-4],fontsize=20, fontweight='bold')
			ax.set_xlabel('rmsd',fontsize=20, fontweight='bold')
			n, bins, patches = ax.hist(dfs[filename]['rmsd'],color='black',bins=100,lw=2,alpha=1,histtype='step')
			ax.set_xlim([0,np.max(bins)])
			if count in plotlist:
				print count
				fig.savefig(outdir+'/'+filename[:-4]+'.'+str(count)+'.'+'-'+str(width)+'-rmsd-hist.png',dpi=72,bbox_inches='tight')
			count+=1
			max_val=np.max(n)
			th=0.2
			th_=0.1
			th_i=0
			max_i=len(n)
			done=0
			for i in range(len(n)):
				if done==0:
					if max_val==n[i]:
						max_i=i
					if i>max_i:
						if max_val*th_>n[i]:
							done=1
					if max_val*th<n[i]:
						th_i=i+1
			rmsd_level=bins[th_i]
			if rmsd_average==1:
				print 'rmsd avg'
				rmsd_mean=np.mean(dfs[filename]['rmsd'][dfs[filename]['rmsd']<=rmsd_level])
				rmsd_std=np.std(dfs[filename]['rmsd'][dfs[filename]['rmsd']<=rmsd_level])
				rmsd_level=rmsd_mean+rmsd_std*(-2*np.log(th))**0.5
			ax.plot([rmsd_level,rmsd_level],[0,max_val],color='darkred',lw=1,alpha=1)
			if count in plotlist:
				print count
				fig.savefig(outdir+'/'+filename[:-4]+'.'+str(count)+'.'+'-'+str(width)+'-rmsd-hist_.png',dpi=72,bbox_inches='tight')
			count+=1
			ax.plot([rmsd_level*rmsd_th,rmsd_level*rmsd_th],[0,max_val],color='darkred',lw=1,ls='--',alpha=1)
			if count in plotlist:
				print count
				fig.savefig(outdir+'/'+filename[:-4]+'.'+str(count)+'.'+'-'+str(width)+'-rmsd-hist__.png',dpi=72,bbox_inches='tight')
			count+=1
			plt.close(fig)
			
			if count in plotlist:
				print count
				fig,ax = plt.subplots(1,1,figsize=(16,6), sharex=True, sharey=False)
				ax.set_title(filename[:-4],fontsize=20, fontweight='bold')
				ax.set_xlabel('position',fontsize=20, fontweight='bold')
				ax.set_ylabel('intensity',fontsize=20, fontweight='bold')
				ax.plot(dfs[filename]['pos'],dfs[filename]['intensity'],color='black',lw=2,alpha=1)
				for idx in dfs[filename]['pos'].index:
					if dfs[filename]['rmsd'].ix[idx]>rmsd_level*rmsd_th:
						ax.scatter(dfs[filename]['pos'].ix[idx],1.01*np.max(dfs[filename]['intensity']),color='darkred',s=10,alpha=1)
				ax.set_xlim([0,np.max(dfs[filename]['pos'])])
				ax.set_ylim([np.min(dfs[filename]['intensity']),1.02*np.max(dfs[filename]['intensity'])])
				fig.savefig(outdir+'/'+filename[:-4]+'.'+str(count)+'.'+'-'+str(width)+'_.png',dpi=72,bbox_inches='tight')
				plt.close(fig)
			count+=1
								
			width2=width
			window2='hanning'
			intensity_max=0
			peak_width_threshold=0.25
			if window2 in ['hanning', 'hamming', 'bartlett', 'blackman']:
				w2=eval('np.'+window2+'(2*width2+1)')
			else:
				window2='flat'
				w2=np.ones(2*width2+1,'d')
			dfs[filename]['int-smooth']=np.convolve(w2/w2.sum(),dfs[filename]['intensity'],'same')
			if count in plotlist:
				print count
				fig,ax = plt.subplots(1,1,figsize=(16,6), sharex=True, sharey=False)
				ax.set_title(filename[:-4],fontsize=20, fontweight='bold')
				ax.set_xlabel('position',fontsize=20, fontweight='bold')
				ax.set_ylabel('intensity',fontsize=20, fontweight='bold')
				ax.plot(dfs[filename]['pos'],dfs[filename]['int-smooth'],color='black',lw=2,alpha=1)
				ax.set_xlim([0,np.max(dfs[filename]['pos'])])
				ax.set_ylim([np.min(dfs[filename]['intensity']),1.02*np.max(dfs[filename]['intensity'])])
				fig.savefig(outdir+'/'+filename[:-4]+'.'+str(count)+'.'+'-'+window+str(width)+'-int-smooth.png',dpi=72,bbox_inches='tight')
				plt.close(fig)
			count+=1
			dfs[filename]['int-diff']=dfs[filename]['int-smooth'].diff()
			if count in plotlist:
				print count
				fig,ax = plt.subplots(1,1,figsize=(16,6), sharex=True, sharey=False)
				ax.set_title(filename[:-4],fontsize=20, fontweight='bold')
				ax.set_xlabel('position',fontsize=20, fontweight='bold')
				ax.set_ylabel('intensity',fontsize=20, fontweight='bold')
				ax.plot(dfs[filename]['pos'],dfs[filename]['int-diff'],color='black',lw=2,alpha=1)
				ax.set_xlim([0,np.max(dfs[filename]['pos'])])
				fig.savefig(outdir+'/'+filename[:-4]+'.'+str(count)+'.'+'-'+window+str(width)+'-int-diff.png',dpi=72,bbox_inches='tight')
				plt.close(fig)
			count+=1
			dfs[filename]['int-diff-bin']=0
			dfs[filename]['int-diff-bin'][dfs[filename]['int-diff']>0]=1
			dfs[filename]['int-diff-bin-diff']=dfs[filename]['int-diff-bin'].diff()
			dfs[filename]['int-diff-bin-diff-max']=dfs[filename]['int-diff-bin-diff']
			dfs[filename]['int-diff-bin-diff-min']=dfs[filename]['int-diff-bin-diff']
			dfs[filename]['int-diff-bin-diff-max'][dfs[filename]['int-diff-bin-diff']>0]=0
			dfs[filename]['int-diff-bin-diff-min'][dfs[filename]['int-diff-bin-diff']<0]=0
			if count in plotlist:
				print count
				fig,ax = plt.subplots(1,1,figsize=(16,6), sharex=True, sharey=False)
				ax.set_title(filename[:-4],fontsize=20, fontweight='bold')
				ax.set_xlabel('position',fontsize=20, fontweight='bold')
				ax.set_ylabel('intensity',fontsize=20, fontweight='bold')
				ax.plot(dfs[filename]['pos'],dfs[filename]['intensity'],color='darkred',lw=2,alpha=1)
				ax.plot(dfs[filename]['pos'],(-dfs[filename]['int-diff-bin-diff-max'])*1.02*(np.max(dfs[filename]['intensity'])-np.min(dfs[filename]['intensity']))+np.min(dfs[filename]['intensity']),color='black',lw=1,alpha=1)
				ax.set_xlim([0,np.max(dfs[filename]['pos'])])
				ax.set_ylim([np.min(dfs[filename]['intensity']),1.02*np.max(dfs[filename]['intensity'])])
				fig.savefig(outdir+'/'+filename[:-4]+'.'+str(count)+'.'+'-'+window+str(width)+'-int-diff-bin-diff-max.png',dpi=72,bbox_inches='tight')
				plt.close(fig)
			count+=1
			if count in plotlist:
				print count
				fig,ax = plt.subplots(1,1,figsize=(16,6), sharex=True, sharey=False)
				ax.set_title(filename[:-4],fontsize=20, fontweight='bold')
				ax.set_xlabel('position',fontsize=20, fontweight='bold')
				ax.set_ylabel('intensity',fontsize=20, fontweight='bold')
				ax.plot(dfs[filename]['pos'],dfs[filename]['intensity'],color='darkred',lw=2,alpha=1)
				ax.plot(dfs[filename]['pos'],(dfs[filename]['int-diff-bin-diff-min'])*1.02*(np.max(dfs[filename]['intensity'])-np.min(dfs[filename]['intensity']))+np.min(dfs[filename]['intensity']),color='black',lw=1,alpha=1)
				ax.set_xlim([0,np.max(dfs[filename]['pos'])])
				ax.set_ylim([np.min(dfs[filename]['intensity']),1.02*np.max(dfs[filename]['intensity'])])
				fig.savefig(outdir+'/'+filename[:-4]+'.'+str(count)+'.'+'-'+window+str(width)+'-int-diff-bin-diff-min.png',dpi=72,bbox_inches='tight')
				plt.close(fig)
			count+=1
								
			peaklist[filename]=[]
			for idx in dfs[filename].index:
				idx_=idx
				if dfs[filename]['int-diff-bin-diff'].ix[idx]==-1:
					peak_pos=dfs[filename]['pos'].ix[idx]
					peak_height=dfs[filename]['intensity'].ix[idx]
					signal_to_noise=0
					done=0
					for i in range(width2):
						if idx-i in dfs[filename].index:
							if dfs[filename]['intensity'].ix[idx-i]<peak_height/2.0:
								done=1
							if dfs[filename]['int-diff-bin-diff'].ix[idx-i]==1:
								done=1
							if done==0:
								if peak_height<dfs[filename]['intensity'].ix[idx-i]:
									peak_height=dfs[filename]['intensity'].ix[idx-i]
									peak_pos=dfs[filename]['pos'].ix[idx-i]
									idx_=idx-i
					done=0
					for i in range(width2):
						if idx+i in dfs[filename].index:
							if dfs[filename]['intensity'].ix[idx+i]<peak_height/2.0:
								done=1
							if dfs[filename]['int-diff-bin-diff'].ix[idx+i]==1:
								done=1
							if done==0:
								if peak_height<dfs[filename]['intensity'].ix[idx+i]:
									peak_height=dfs[filename]['intensity'].ix[idx+i]
									peak_pos=dfs[filename]['pos'].ix[idx+i]
									idx_=idx+i
					
					i_min_minus=0
					done=0
					for i in range(width2*5):
						if idx_-i in dfs[filename].index:
							if dfs[filename]['int-diff-bin-diff'].ix[idx_-i]==1:
								done=1
							if done==0:
								i_min_minus=i
					i_min_plus=0
					done=0
					for i in range(width2*5):
						if idx_+i in dfs[filename].index:
							if dfs[filename]['int-diff-bin-diff'].ix[idx_+i]==1:
								done=1
							if done==0:
								i_min_plus=i
					local_bgr=(dfs[filename]['intensity'].ix[idx_-i_min_minus]+dfs[filename]['intensity'].ix[idx_+i_min_plus])/2.0
					peak_minus=-1	
					i_peak_minus=0
					done=0
					for i in range(i_min_minus):
						if idx_-i in dfs[filename].index:	
							if (dfs[filename]['intensity'].ix[idx_-i]-local_bgr)<(peak_height-local_bgr)*peak_width_threshold:
								done=1
							if done==0:
								peak_minus=dfs[filename]['pos'].ix[idx_-i]
								i_peak_minus=i
					peak_plus=-1
					i_peak_plus=0
					done=0
					for i in range(i_min_plus):
						if idx_+i in dfs[filename].index:
							if (dfs[filename]['intensity'].ix[idx_+i]-local_bgr)<(peak_height-local_bgr)*peak_width_threshold:
								done=1
							if done==0:
								peak_plus=dfs[filename]['pos'].ix[idx_+i]
								i_peak_plus=i
					avg=dfs[filename]['intensity'].ix[idx_]
					avg_count=1
					for i in range(i_peak_minus):
						if idx_-1-i in dfs[filename].index:
							avg+=dfs[filename]['intensity'].ix[idx_-1-i]
							avg_count+=1
					for i in range(i_peak_plus):
						if idx_+1+i in dfs[filename].index:
							avg+=dfs[filename]['intensity'].ix[idx_+1+i]
							avg_count+=1
					avg/=avg_count
					rmsd=(dfs[filename]['intensity'].ix[idx_]-avg)**2
					rmsd_count=1
					for i in range(i_peak_minus):
						if idx_-1-i in dfs[filename].index:
							rmsd+=(dfs[filename]['intensity'].ix[idx_-1-i]-avg)**2
							rmsd_count+=1
					for i in range(i_peak_plus):
						if idx_+1+i in dfs[filename].index:
							rmsd+=(dfs[filename]['intensity'].ix[idx_+1+i]-avg)**2
							rmsd_count+=1
					rmsd/=rmsd_count
					rmsd=rmsd**0.5
					signal_to_noise=rmsd/rmsd_level
					if 30<peak_pos and peak_pos<6020:
						peaklist[filename].append((peak_pos,peak_minus,peak_plus,peak_plus-peak_minus,peak_height,peak_height-intensity_level,avg,local_bgr,signal_to_noise))
			f=open(outdir+'/'+filename[:-4]+'-'+window+str(width)+'-peaklist.text','w')
			f.write('pos\twhm-\twhm+\tdelta\theight\theight_norm\tavg\tbgr\tsig_to_noise\n')
			for (pos,minus,plus,delta,height,height_norm,avg,bgr,sig_to_noise) in peaklist[filename]:
				f.write(str(pos)+'\t'+str(minus)+'\t'+str(plus)+'\t'+str(delta)+'\t'+str(height)+'\t'+str(height_norm)+'\t'+str(avg)+'\t'+str(bgr)+'\t'+str(sig_to_noise)+'\n')
			f.close()
			
			fig,ax = plt.subplots(1,1,figsize=(16,6), sharex=True, sharey=False)
			ax.set_title(filename[:-4],fontsize=20, fontweight='bold')
			ax.set_xlabel('position',fontsize=20, fontweight='bold')
			ax.set_ylabel('intensity',fontsize=20, fontweight='bold')
			ax.plot(dfs[filename]['pos'],dfs[filename]['intensity'],color='darkred',lw=2,alpha=1)
			for (pos,minus,plus,delta,height,height_norm,avg,bgr,sig_to_noise) in peaklist[filename]:
				if sig_to_noise>=rmsd_th and height_norm>height_th:
					ax.scatter([pos],[1.01*height],facecolor='blue',edgecolor='darkblue',s=60,lw=2,alpha=1)
			ax.set_xlim([0,np.max(dfs[filename]['pos'])])
			ax.set_ylim([np.min(dfs[filename]['intensity']),1.02*np.max(dfs[filename]['intensity'])])
			if count in plotlist:
				print count
				fig.savefig(outdir+'/'+filename[:-4]+'.'+str(count)+'.'+'-'+window+str(width)+'-peaks.png',dpi=72,bbox_inches='tight')
			count+=1
			for (pos,minus,plus,delta,height,height_norm,avg,bgr,sig_to_noise) in peaklist[filename]:
				if sig_to_noise>=rmsd_th and height_norm>height_th:
					ax.plot([minus,plus],[peak_width_threshold*(height-bgr)+bgr,peak_width_threshold*(height-bgr)+bgr],color='darkblue',lw=2,alpha=1)
			if count in plotlist:
				print count
				fig.savefig(outdir+'/'+filename[:-4]+'.'+str(count)+'.'+'-'+window+str(width)+'-peaks-widths.png',dpi=72,bbox_inches='tight')
			count+=1
			for (pos,minus,plus,delta,height,height_norm,avg,bgr,sig_to_noise) in peaklist[filename]:
				if sig_to_noise>=rmsd_th and height_norm>height_th:
					ax.plot([minus,plus],[bgr,bgr],color='gray',lw=1,alpha=1)
			if count in plotlist:
				print count
				fig.savefig(outdir+'/'+filename[:-4]+'.'+str(count)+'.'+'-'+window+str(width)+'-peaks-widths-bgr.png',dpi=72,bbox_inches='tight')
				plt.close(fig)
			count+=1
			for (pos,minus,plus,delta,height,height_norm,avg,bgr,sig_to_noise) in peaklist[filename]:
				if sig_to_noise>=rmsd_th and height_norm>height_th:
					ax.plot([minus,plus],[avg,avg],color='black',lw=1,alpha=1)
			if count in plotlist:
				print count
				fig.savefig(outdir+'/'+filename[:-4]+'.'+str(count)+'.'+'-'+window+str(width)+'-peaks-widths-bgr-avg.png',dpi=72,bbox_inches='tight')
				plt.close(fig)
			count+=1
		
