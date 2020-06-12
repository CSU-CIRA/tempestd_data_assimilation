#!/usr/bin/python
# -*- coding: utf-8 -*-

##################################################################
# Python Program Name: compute_anomaly_correlation.py
# Programmer: Ting-Chi Wu @ CIRA/CSU
# Date: 03/04/2020
# Program Log: compute anomaly correlation coefficient for the 
#              500 mb Geopotential Height and TPW using the following
#              equation:
# 
#              ACC = mean((f-c)(a-c))/sqrt(mean((f-c)^2)*mean((a-c)^2))
#   
#              with the following three datasets:
#              1) GFS Production Analysis (a)
#              2) Forecasts from FV3GFS cycled experiments (f)
#              3) NCEP-NCAR Reanalysis 1 Climatology (c)
#             
# Reference: ECMWF_ACC_definition.pdf (owner: Erik Andersson) 
# URL: http://www.atmos.albany.edu/daes/atmclasses/atm401/spring_2016/ppts_pdfs/ECMWF_ACC_definition.pdf
##################################################################

import numpy as np
import os
from os import path
from pylab import *
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.basemap import Basemap
import datetime as dt

gas_const = 287.0 # J/K/kg
one_rad = 180/np.pi

pwd = '/mnt/ssdatenas/ting-chi/scripts'
#caseperiod = 'dec2018'
caseperiod = 'may2019'
#cdump = 'gfs'
cdump = 'gdas'
#varname, varmax = '500HGT', 0.86
varname, varmax = 'TPW', 0.68
nfcst = 41


expname1, title1 = 'control_{}'.format(caseperiod), 'Control'
expname2, title2 = 'mhs_{}'.format(caseperiod), 'AddMHS'
##expname3, title3 = 'tempestd_{}_calidiff'.format(caseperiod), 'AddTEMPESTD(TEMP)'
#expname3, title3 = 'tempestd_[]'.format(caseperiod), 'AddTEMPEST-D (MixQC+NoCali)'
expname3, title3 = 'tempestd_{}_newcalidiff_cyc1'.format(caseperiod), 'AddTEMPESTD'
expname4, title4 = 'tempestd_{}_newcalidiff_cyc1noch2'.format(caseperiod), 'AddTEMPESTD2'
##expname3, title3 = 'tempestd_dec2018_newcalidiff_retocean', 'AddTEMPESTD(TEMPQC+NewCali)'
##expname4, title4 = 'tempestd_dec2018_newcalidiff_mixqc', 'AddTEMPESTD(TEMPQC+MHSQC+NewCali)'
##expname4, title4 = 'tempestd_dec2018_mhsqc', 'AddTEMPESTD(MHSQC)'
##expname5, title5 = 'tempestd_dec2018_retocean', 'AddTEMPESTD(TEMPQC)'
##expname6, title6 = 'tempestd_dec2018_calidiff', 'AddTEMPESTD(TEMPQC+MHSQC)'


#----------------------------------------------------#
#   Keep as is
#----------------------------------------------------#
#
#  ACC files 
#
lfname = 'list_of_acc_files.txt'
os.system('ls -1 acc.{}_analysis/acc_{}.{}*.txt > list_of_acc_files.txt'.format(cdump, varname, expname1))
#os.system('ls -1 acc_{}.{}*00.txt > list_of_acc_files.txt'.format(varname, expname1))
#os.system('ls -1 acc_{}.{}*12.txt > list_of_acc_files.txt'.format(varname, expname1))
#os.system('ls -1 acc_{}.{}*.txt > list_of_acc_files.txt'.format('500HGT', expname))
f1names = [line.rstrip('\n') for line in open('{}'.format(lfname))]
os.remove('{}'.format(lfname))

os.system('ls -1 acc.{}_analysis/acc_{}.{}*.txt > list_of_acc_files.txt'.format(cdump, varname, expname2))
#os.system('ls -1 acc_{}.{}*00.txt > list_of_acc_files.txt'.format(varname, expname2))
#os.system('ls -1 acc_{}.{}*12.txt > list_of_acc_files.txt'.format(varname, expname2))
#os.system('ls -1 acc_{}.{}*.txt > list_of_acc_files.txt'.format('TPW', expname))
f2names = [line.rstrip('\n') for line in open('{}'.format(lfname))]
os.remove('{}'.format(lfname))

os.system('ls -1 acc.{}_analysis/acc_{}.{}*.txt > list_of_acc_files.txt'.format(cdump, varname, expname3))
#os.system('ls -1 acc_{}.{}*00.txt > list_of_acc_files.txt'.format(varname, expname3))
#os.system('ls -1 acc_{}.{}*12.txt > list_of_acc_files.txt'.format(varname, expname3))
f3names = [line.rstrip('\n') for line in open('{}'.format(lfname))]
os.remove('{}'.format(lfname))

os.system('ls -1 acc.{}_analysis/acc_{}.{}*.txt > list_of_acc_files.txt'.format(cdump, varname, expname4))
#os.system('ls -1 acc_{}.{}*00.txt > list_of_acc_files.txt'.format(varname, expname4))
#os.system('ls -1 acc_{}.{}*12.txt > list_of_acc_files.txt'.format(varname, expname4))
f4names = [line.rstrip('\n') for line in open('{}'.format(lfname))]
os.remove('{}'.format(lfname))

#os.system('ls -1 acc.{}_analysis/acc_{}.{}*.txt > list_of_acc_files.txt'.format(cdump, varname, expname5))
##os.system('ls -1 acc_{}.{}*00.txt > list_of_acc_files.txt'.format(varname, expname5))
##os.system('ls -1 acc_{}.{}*12.txt > list_of_acc_files.txt'.format(varname, expname5))
#f5names = [line.rstrip('\n') for line in open('{}'.format(lfname))]
#os.remove('{}'.format(lfname))

#os.system('ls -1 acc.{}_analysis/acc_{}.{}*.txt > list_of_acc_files.txt'.format(cdump, varname, expname6))
##os.system('ls -1 acc_{}.{}*00.txt > list_of_acc_files.txt'.format(varname, expname6))
##os.system('ls -1 acc_{}.{}*12.txt > list_of_acc_files.txt'.format(varname, expname6))
#f6names = [line.rstrip('\n') for line in open('{}'.format(lfname))]
#os.remove('{}'.format(lfname))

nfnames = len(f1names)
#nfnames = 9
print('ACC files (# = {}) :'.format(nfnames))
for jj in range(nfnames):
    print('{}, {}, and {}'.format(f1names[jj], f2names[jj], f3names[jj]))

#----------------------------------------------------#
#   Define Functions to be Used Later
#----------------------------------------------------#

def read_acc(filename):

#    with open(filename) as f:
#         dat = [line.rstrip() for line in f]
#    f.close()

    fid = open(filename,'r')
    print('Open to read {}'.format(filename))
    ncol = 4
    dat=[[] for i in range(ncol)]
    count=0
    for line in fid.readlines(): # this reads one line at a time
#        print line
        spl_line = line.split(',') # split into a list of strings by whitespace
        for idat in range(ncol):
            dat[idat].append(float(spl_line[idat])) # convert from string to float and append
        count+=1 # same as count=count+1
    dat = np.array(dat)

    return(dat)

#----------------------------------------------------#
#   Read and Manipulate Data
#----------------------------------------------------#

#accs = [read_acc(ffname) for ffname in fnames]
#accs = np.zeros([nfnames,2,nfcst])
#for jj in range(nfnames):
#    accs[jj,0,:] = read_acc(f1names[jj])
#    accs[jj,1,:] = read_acc(f2names[jj])

#accs = np.zeros([nfnames,3,4,nfcst])
accs = np.zeros([nfnames,4,4,nfcst])
#accs = np.zeros([nfnames,5,4,nfcst])
#accs = np.zeros([nfnames,6,4,nfcst])
for jj in range(nfnames):
    accs[jj,0,:,:] = read_acc(f1names[jj])
    accs[jj,1,:,:] = read_acc(f2names[jj])
    accs[jj,2,:,:] = read_acc(f3names[jj])
    accs[jj,3,:,:] = read_acc(f4names[jj])
#    accs[jj,4,:,:] = read_acc(f5names[jj])
#    accs[jj,5,:,:] = read_acc(f6names[jj])

accm = np.mean(accs,axis=0)
regions = ['Global', 'North Hemisphere', 'South Hemisphere', 'Tropics']
regsn = ['gl', 'nh', 'sh', 'tr']

#----------------------------------------------------#
#   Make a plot of ACC as a function of time 
#----------------------------------------------------#

def plot_acc_all_regions():
    '''
    Make a 4-panel plot of ACC as a function of time for Global, NH, SH, and Tropics
    '''
    fig = plt.figure(num=None, figsize=(24.0, 16.0), dpi=80, facecolor='w', edgecolor='k')
    fig.suptitle('Anomaly Correlation Coefficient: {} (averaged over all {} cycles)'.format(varname, caseperiod),fontsize=15,fontweight='bold')
    #fig.suptitle('Anomaly Correlation Coefficient: {} (averaged over all 00 UTC cycles)'.format(varname),fontsize=15,fontweight='bold')
    #fig.suptitle('Anomaly Correlation Coefficient: {} (averaged over all 12 UTC cycles)'.format(varname),fontsize=15,fontweight='bold')

    #xx = np.linspace(0,240,nfcst)
    xx = np.linspace(0,120,21)
    xx_str = [np.str(np.int(x)) for x in xx]

    for iplot in range(4):
       plotnum = iplot + 1
       ax = fig.add_subplot(2,2,plotnum)
#       plt.plot(xx,accm[0,iplot,:],'-k')
#       plt.plot(xx,accm[1,iplot,:],'-b')
#       plt.plot(xx,accm[2,iplot,:],'-g')
       plt.plot(xx,accm[0,iplot,:21],'-k')
       plt.plot(xx,accm[1,iplot,:21],'-b')
       plt.plot(xx,accm[2,iplot,:21],'-g')
       plt.plot(xx,accm[3,iplot,:21],'-r')
#       plt.plot(xx,accm[4,iplot,:21],color='orange')
#       plt.plot(xx,accm[5,iplot,:21],'-m')
#       plt.legend(['{}'.format(title1),'{}'.format(title2)],loc='best',fontsize=16)
#       plt.legend(['{}'.format(title1),'{}'.format(title2), '{}'.format(title3)],loc='best',fontsize=16)
       plt.legend(['{}'.format(title1),'{}'.format(title2), '{}'.format(title3), '{}'.format(title4)],loc='best',fontsize=16)
#       plt.legend(['{}'.format(title1),'{}'.format(title2), '{}'.format(title3), '{}'.format(title4), '{}'.format(title5)],loc='best',fontsize=16)
#       plt.legend(['{}'.format(title1),'{}'.format(title2), '{}'.format(title3), '{}'.format(title4), '{}'.format(title5), '{}'.format(title6)],loc='best',fontsize=16)
#       plt.legend(['{}'.format(title1),'{}'.format(title2), '{}'.format(title3), '{}'.format(title5), '{}'.format(title6)],loc='best',fontsize=16)
       plt.xlabel('Forecast Lead Time (hour)',fontsize=12)
       plt.ylabel('ACC',fontsize=12)
#       plt.title('Anomaly Correlation Coefficient from {} (averaged over all cycles)'.format(expname),fontsize=20)
       plt.title('{}'.format(regions[iplot]),fontsize=20)
#       ax.set_xlim(0,240)
#       ax.set_ylim(0.3,1.0)
#       ax.set_ylim(0.9,1.0)
       ax.set_ylim(0.7,1.0)
       plt.xticks(xx,xx_str,fontsize=10,rotation=45)
#       plt.tight_layout()
       plt.grid()
#       plt.savefig('acc_{}_5days.averaged.png'.format(varname))
#       plt.savefig('acc_{}_00UTC.averaged.png'.format(varname))
#       plt.savefig('acc_{}_12UTC.averaged.png'.format(varname))

def plot_acc_one_region(ireg):
    '''
    Make a plot of ACC as a function of time for selected region (ireg)
    (ireg: 0. Global, 1. NH, 2. SH, and 3. Tropics)
    '''
    fig = plt.figure(num=None, figsize=(15.0, 12.0), dpi=80, facecolor='w', edgecolor='k')
#    fig.suptitle('Anomaly Correlation Coefficient: {} (averaged over all {} cycles)'.format(varname, caseperiod),fontsize=15,fontweight='bold')
    xx = np.linspace(0,120,21)
    xx_str = [np.str(np.int(x)) for x in xx]
    ax = fig.add_subplot(2,1,1)
    ax.plot(xx,accm[0,ireg,:21],'-b', label='{}'.format(title1))
    ax.plot(xx,accm[1,ireg,:21],'-g', label='{}'.format(title2))
    ax.plot(xx,accm[2,ireg,:21],color='orange', label='{}'.format(title3))
#    ax.plot(xx,accm[2,ireg,:21],'-r', label='{}'.format(title3))
#    plt.plot(xx,accm[3,ireg,:21],color='orange')
#    plt.legend(['{}'.format(title1),'{}'.format(title2), '{}'.format(title3), '{}'.format(title4)],loc='best',fontsize=16)
#    plt.legend(['{}'.format(title1),'{}'.format(title2), '{}'.format(title3)],loc='best',fontsize=16)
    plt.xlabel('Forecast Lead Time (hour)',fontsize=16)
#    plt.ylabel('ACC',fontsize=12)
    ax.set_ylabel('ACC',fontsize=16)
    plt.title('{}: {} (averaged over all {} cycles)'.format(regions[ireg], varname, caseperiod),fontsize=20)
#    ax.set_ylim(0.9,1.0)
    if ireg == 0:
       ax.set_ylim(varmax,1.0)
    elif ireg == 1:
       ax.set_ylim(varmax,1.0)
    elif ireg == 2:
       ax.set_ylim(varmax,1.0)
    elif ireg == 3:
       ax.set_ylim(varmax,1.0)
    plt.xticks(xx,xx_str,fontsize=16,rotation=45)
    ax.tick_params(axis='y',labelsize=16)
    ax.grid(linewidth=2)
    plt.legend(bbox_to_anchor=(0.5, -0.2),loc='upper center',ncol=3,fontsize=16)
    
    ax2 = ax.twinx()
    ax2.set_ylabel('ACC Difference Relative to Control',fontsize=16)
    yy1 = accm[1,ireg,:21]-accm[0,ireg,:21]
    yy2 = accm[2,ireg,:21]-accm[0,ireg,:21]
    ylimmax = max(yy1.max(),yy2.max())
    ylimmin = min(yy1.min(),yy2.min())
    ax2.plot(xx,yy1,'--g', label='{} - Control'.format(title2))
    ax2.plot(xx,yy2,linestyle='--',color='orange', label='{} - Control'.format(title3))
#    ax2.set_ylim(ylimmin,ylimmax)
    ax2.set_ylim(-0.005,0.015)
    ax2.tick_params(axis='y',labelsize=16)
    ax2.axhline(y=0,c='black')
    ax2.grid(linestyle='--')

    plt.legend(bbox_to_anchor=(0.5, -0.3),loc='upper center',ncol=3,fontsize=16)
#    ax = fig.add_subplot(2,1,2)
#    plt.plot(xx,accm[1,ireg,:21]-accm[0,ireg,:21],'-g')
#    plt.plot(xx,accm[2,ireg,:21]-accm[0,ireg,:21],'-r')
#    plt.plot(xx,accm[3,ireg,:21]-accm[0,ireg,:21],color='orange')
#    plt.legend(['{}'.format(title2), '{}'.format(title3), '{}'.format(title4)],loc='best',fontsize=16)
#    plt.xlabel('Forecast Lead Time (hour)',fontsize=12)
#    plt.ylabel('ACC Difference Relative to Control',fontsize=12)
#    plt.title('{}'.format(regions[ireg]),fontsize=20)
##    ax.set_ylim(0.9,1.0)
##    ax.set_ylim(0.7,1.0)
#    plt.xticks(xx,xx_str,fontsize=10,rotation=45)
#    plt.grid()
    plt.savefig('acc_{}_{}_5days.averaged_{}.png'.format(varname, regsn[ireg], caseperiod))


if __name__ == '__main__':
    plt.ion()
#    plot_acc_all_regions()
    plot_acc_one_region(0)
    plot_acc_one_region(1)
    plot_acc_one_region(2)
    plot_acc_one_region(3)
    plt.show()


