#!/usr/bin/python
# -*- coding: utf-8 -*-

##################################################################
# Python Program Name: analysis_bias_comparison.py
# Programmer: Ting-Chi Wu @ CIRA/CSU
# Date: 03/31/2020
#
# Program Log: compute bias against GFS Production Analysis
#              with the following two datasets:
#              1) GFS/GDAS Production Analysis (a)
#              2) Analysis from FV3GFS cycled experiments (f)
#
# Usage (in ipython): %run analysis_bias_comparison.py caseperiod cyctime varname varmax
# Example: %run analysis_bias_comparison.py dec2018 2018121000 SPH 0.45 
##################################################################

import numpy as np
import os
import sys
from os import path
from pylab import *
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.basemap import Basemap
import datetime as dt

#cdump='gfs'
cdump='gdas'
#apath = '/mnt/ssdatenas/ting-chi/gfs_production/keep'
apath = '/mnt/ssdatenas/ting-chi/{}_production'.format(cdump)
fpath = '/mnt/ssdatenas/ting-chi/fv3gfs_experiments'

caseperiod = sys.argv[1]
cyctime = sys.argv[2]
varname = sys.argv[3]
varmax = float(sys.argv[4])

#expnames = ['control_dec2018', 'mhs_dec2018', 'tempestd_dec2018']
#expnames = ['control_dec2018', 'mhs_dec2018', 'tempestd_dec2018_calidiff']
#expnames = ['control_dec2018', 'mhs_dec2018', 'tempestd_dec2018_newcalidiff_retocean']
#expnames = ['control_dec2018', 'mhs_dec2018', 'tempestd_dec2018_newcalidiff_mixqc']
expnames = ['control_{}'.format(caseperiod), 
            'mhs_{}'.format(caseperiod), 
            'tempestd_{}_newcalidiff_cyc1'.format(caseperiod), 
            'tempestd_{}_newcalidiff_cyc1noch2'.format(caseperiod)]

titles = ['Control',
          'AddMHS',
          'AddTEMPESTD',
          'AddTEMPESTD2']

nexp = len(expnames)

if caseperiod == 'dec2018':
   cyctimes = ['2018120812', '2018120900', '2018120912',
               '2018121000', '2018121012', '2018121100',
               '2018121112', '2018121200', '2018121212']
elif caseperiod == 'may2019':
   cyctimes = ['2019051212', '2019051300', '2019051312',
               '2019051400', '2019051412', '2019051500',
               '2019051512', '2019051600', '2019051612',
               '2019051700', '2019051712', '2019051800',
               '2019051812', '2019051900', '2019051912',
               '2019052000', '2019052012', '2019052100',
               '2019052112', '2019052200', '2019052212']

ncyc = len(cyctimes)

#----------------------------------------------------#
nx = 1536
ny = 768
if (varname == 'TPW' or varname == '500HGT' or varname == 'PSFC'):
  nz = 1
else:
  nz = 64

#
#     1. GFS Production Analysis Files (mapped from C768 to C384 grid)
#
laname = 'list_of_{}_files.production.txt'.format(varname)
for i in range(ncyc):
    yyyymmdd = cyctimes[i][:8]
    hh = cyctimes[i][-2:]
    if ( i == 0 ):
       os.system('ls -1 {}/{}.{}/{}/{}000000 > {}'.format(apath, cdump, yyyymmdd, hh, varname, laname))
    else:
       os.system('ls -1 {}/{}.{}/{}/{}000000 >> {}'.format(apath, cdump, yyyymmdd, hh, varname, laname))
afnames = [line.rstrip('\n') for line in open('{}'.format(laname))]
os.remove('{}'.format(laname))
nafnames = len(afnames)
print('Production Analysis files (# = {}) :'.format(nafnames))
for jj in range(nafnames):
    print('{}'.format(afnames[jj]))

#
#     2. FV3GFS Experiment Analysis Files (C384 grid)
#
for j in range(nexp):
   lfname = 'list_of_{}_files.{}.txt'.format(varname,expnames[j])
   for i in range(ncyc):
      if ( i == 0 ):
         os.system('ls -1 {}/{}/gfs.analysis/gfs.{}/{}000000 > {}'.format(fpath, expnames[j], cyctimes[i], varname, lfname))
      else:
         os.system('ls -1 {}/{}/gfs.analysis/gfs.{}/{}000000 >> {}'.format(fpath, expnames[j], cyctimes[i], varname, lfname))
   ffnames = [line.rstrip('\n') for line in open('{}'.format(lfname))]
   os.remove('{}'.format(lfname))
   if j == 0:
      f1fnames = ffnames
   elif j == 1:
      f2fnames = ffnames
   elif j == 2:
      f3fnames = ffnames
   elif j == 3:
      f4fnames = ffnames

nffnames = len(ffnames)
print('FV3GFS Experiment Analysis Files (# = {}) :'.format(nffnames))
for jj in range(nffnames):
    print('{}'.format(f1fnames[jj]))
    print('{}'.format(f2fnames[jj]))
    print('{}'.format(f3fnames[jj]))
    print('{}'.format(f4fnames[jj]))
    print('--------------------------------')
#-----------------------------------------------------------------------------------------#

latfname = '{}/control_{}/gfs.forecast_derived/gfs.{}/xlat.dat'.format(fpath,caseperiod,cyctimes[0])
lonfname = '{}/control_{}/gfs.forecast_derived/gfs.{}/xlon.dat'.format(fpath,caseperiod,cyctimes[0])
pressfname = '{}/control_{}/gfs.analysis/gfs.{}/PRESS000000'.format(fpath,caseperiod,cyctimes[0])
print('Forecast Latitude file: {}'.format(latfname))
print('Forecast Longitude file: {}'.format(latfname))
print('Forecast Pressure file: {}'.format(pressfname))

#----------------------------------------------------#
#   Define Functions to be Used Later
#----------------------------------------------------#

def readbin(filename,nlev,ncol,nrow):

   dat_h = open(filename,'rb')
   if (nlev == 1 ) :
     shape = (nrow,ncol)
     dat = np.fromfile(dat_h,dtype='f4',count=ncol*nrow).reshape(shape)
   else:
     shape = (nlev,nrow,ncol)
     dat = np.fromfile(dat_h,dtype='f4',count=nlev*ncol*nrow).reshape(shape)

     dat_h.close()

   return(dat)

def find_idx_closest_value(array1d,value):
    
    idx = min(range(len(array1d)), key=lambda i: abs(array1d[i]-value))

    return(idx)
    

#----------------------------------------------------#
#   Read Data
#----------------------------------------------------#

# gfs/gdas production analysis files: 
anal = [readbin(afname, nz, nx, ny) for afname in afnames]

# fv3gfs experiment analysis files:
fcst1 = [readbin(f1fname, nz, nx, ny) for f1fname in f1fnames]
fcst2 = [readbin(f2fname, nz, nx, ny) for f2fname in f2fnames]
fcst3 = [readbin(f3fname, nz, nx, ny) for f3fname in f3fnames]
fcst4 = [readbin(f4fname, nz, nx, ny) for f4fname in f4fnames]


# lat/lon files
lat = readbin(latfname, 1, nx, ny)
lon = readbin(lonfname, 1, nx, ny)

# presssure files
press = readbin(pressfname, 64, nx, ny)

#----------------------------------------------------#
#   Compute Anomaly Correlation Coefficients
#----------------------------------------------------#

f1mas = [f1 - a for f1, a in zip(fcst1, anal) ]
f2mas = [f2 - a for f2, a in zip(fcst2, anal) ]
f3mas = [f3 - a for f3, a in zip(fcst3, anal) ]
f4mas = [f4 - a for f4, a in zip(fcst4, anal) ]


if ( cyctime == 'all' ):
  mf1mas = np.mean(f1mas,axis=0)
  mf2mas = np.mean(f2mas,axis=0)
  mf3mas = np.mean(f3mas,axis=0)
  mf4mas = np.mean(f4mas,axis=0)
else:
  for i in range(ncyc):
      if cyctimes[i]==cyctime:
        icyc=i
  mf1mas = f1mas[icyc]
  mf2mas = f2mas[icyc]
  mf3mas = f3mas[icyc]
  mf4mas = f4mas[icyc]

if ( varname == 'SPH' ):
   mf1mas = mf1mas * 1000.0 # convert from kg/kg to g/kg
   mf2mas = mf2mas * 1000.0 # convert from kg/kg to g/kg
   mf3mas = mf3mas * 1000.0 # convert from kg/kg to g/kg
   mf4mas = mf4mas * 1000.0 # convert from kg/kg to g/kg

if (nz == 1): 
   zmf1mas = np.mean(mf1mas,axis=1)
   zmf2mas = np.mean(mf2mas,axis=1)
   zmf3mas = np.mean(mf3mas,axis=1)
   zmf4mas = np.mean(mf4mas,axis=1)
else:
   zmf1mas = np.mean(mf1mas,axis=2)
   zmf2mas = np.mean(mf2mas,axis=2)
   zmf3mas = np.mean(mf3mas,axis=2)
   zmf4mas = np.mean(mf4mas,axis=2)

#plotmax = np.int(max(np.abs(mamcs.min()),np.abs(mamcs.max())))
#plotmin = -plotmax

#----------------------------------------------------#
#   Make a plot of ACC as a function of time 
#----------------------------------------------------#

def plot_mean_anomaly(plotlev):
  
   print('plot mean bias map')

   fig = figure(num=None, figsize=(24.0, 18.0), dpi=80, facecolor='w', edgecolor='k')
   fig.suptitle('Analysis Bias: {} (averaged over all {} cycles)'.format(varname, caseperiod),fontsize=20,fontweight='bold')
   for iplot in range(4):
      plotnum = iplot + 1
 
      ax = fig.add_subplot(2,2,plotnum)
      m = Basemap(projection='mill',llcrnrlat=-90.,urcrnrlat=90.,
                  llcrnrlon=0.,urcrnrlon=360.0,resolution='l')

      m.drawcoastlines()
      m.drawcountries()
      m.drawstates(linewidth=1,linestyle='solid',color='k')
      # draw parallels.
      parallels = np.arange(-90.,90,10.)
      m.drawparallels(parallels,labels=[1,0,0,0],fontsize=15)
      # draw meridians
      meridians = np.arange(-180.,180.,30.)
      m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=15)

      xx, yy = m(lon,lat)
      plotmax=2.5
      plotmin=-plotmax
#      plotcmap = plt.cm.seismic
      plotcmap = plt.cm.BrBG
      if plotnum == 1:
         plotvar = mf1mas[plotlev,:,:]
         plottitle = '{} - Production Analysis'.format(titles[0])

      elif plotnum == 2:
         plotvar = mf2mas[plotlev,:,:]
         plottitle = '{} - Production Analysis'.format(titles[1])

      elif plotnum == 3:
         plotvar = mf3mas[plotlev,:,:]
         plottitle = '{} - Production Analysis'.format(titles[2])

      elif plotnum == 4:
         plotvar = mf4mas[plotlev,:,:]
         plottitle = '{} - Production Analysis'.format(titles[3])

#      m.pcolor(xx,yy,plotvar,vmin=plotmin,vmax=plotmax,cmap=plotcmap)
      m.contourf(xx,yy,plotvar,20,vmin=plotmin,vmax=plotmax,cmap=plotcmap)
#      m.pcolor(xx,yy,plotvar)
      m.colorbar()
      ax.set_title(plottitle,fontsize=18,fontweight='bold')
      fig.subplots_adjust(wspace=.2,hspace = .1)
   fig.savefig('analysis_mean_bias_{}_{}.{}.png'.format(varname, caseperiod, cyctime))

def plot_zmean_bias():

   print('plot zonal mean')

   if ( nz == 1 ): 
     fig = figure(num=None, figsize=(8.0, 12.0), dpi=80, facecolor='w', edgecolor='k')
     ax = fig.add_subplot(1,1,1)
     plt.plot(zmf1mas,lat[:,0],'-b')
     plt.plot(zmf2mas,lat[:,0],'-g')
     plt.plot(zmf3mas,lat[:,0],'-r')
     plt.plot(zmf4mas,lat[:,0],color='orange')
#     ax.set_xlim(-2.0,0.5)
#     ax.set_xlim(-3.0,1.0)
     plt.grid()
#     plt.legend(['Analysis','{}'.format(expnames[0]),'{}'.format(expnames[1]),'{}'.format(expnames[2])],loc='best')
#     plt.legend(['{}'.format(expnames[0]),'{}'.format(expnames[1]),'{}'.format(expnames[2])],loc='best')
     plt.legend(['{}'.format(expnames[0]),'{}'.format(expnames[1]),'{}'.format(expnames[2]), '{}'.format(expnames[3])],loc='best')
     plt.title('Zonal Mean {} Bias (f-a averaged over all {} cycles): Analysis'.format(varname, caseperiod))
   else:
#     varmax= 0.5
#     varmax= 0.45
     xx = np.zeros([nz,ny])
     yy = np.zeros([nz,ny])
#     zmpress = np.mean(press, axis=2)
     zmpress = press[:,394,700]

     for jj in range(nz):
       xx[jj,:] = lat[:,0]
     for ii in range(ny):
       yy[:,ii] = zmpress

#     fig = figure(num=None, figsize=(16.0, 14.0), dpi=80, facecolor='w', edgecolor='k')
     fig = figure(num=None, figsize=(10.0, 21.0), dpi=80, facecolor='w', edgecolor='k')
     fig.suptitle('Zonal Mean Analysis Bias: {} (averaged over all {} cycles)'.format(varname, caseperiod),fontsize=20,fontweight='bold')
#     for iplot in range(4):
     for iplot in range(3):
        plotnum = iplot + 1
#        ax = fig.add_subplot(2,2,plotnum)
        ax = fig.add_subplot(3,1,plotnum)
        
        if plotnum == 1:
           plotvar = zmf1mas 
           plottitle = '{} - Production'.format(titles[0])
        elif plotnum == 2:
           plotvar = zmf2mas 
           plottitle = '{} - Production'.format(titles[1])
        elif plotnum == 3:
           plotvar = zmf3mas 
           plottitle = '{} - Production'.format(titles[2])
#        elif plotnum == 4:
#           plotvar = zmf4mas 
#           plottitle = '{} - Production'.format(titles[3])

#        plt.pcolor(xx,yy,plotvar,vmin=-varmax,vmax=varmax,cmap=plt.cm.seismic)
#        plt.pcolor(xx[:50,:],yy[:50,:],plotvar[:50,:],vmin=-varmax,vmax=varmax,cmap=plt.cm.seismic)
#        plt.contourf(xx,yy,plotvar,10,vmin=-varmax,vmax=varmax,cmap=plt.cm.seismic)
#        plt.pcolor(xx,yy,plotvar,vmin=-varmax,vmax=varmax,cmap=plt.cm.BrBG)
#        plt.contourf(xx[5:50,:],yy[5:50,:],plotvar[5:50,:],vmin=-varmax,vmax=varmax,cmap=plt.cm.BrBG)
        plt.contourf(xx,yy,plotvar,vmin=-varmax,vmax=varmax,cmap=plt.cm.BrBG)
#        plt.contourf(xx,yy,plotvar,50,vmin=-varmax,vmax=varmax,cmap=plt.cm.BrBG)
#        plt.pcolor(xx,yy,plotvar)
#        ax.set_ylim(100,1000)
#        ax.set_ylim(500,1000)
        ax.set_ylim(300,1000)
        plt.gca().invert_yaxis()
        cbar = plt.colorbar()
        cbar.set_label(label='[g/kg]',size=16)
        cbar.ax.tick_params(labelsize=16) 
        plt.grid()
        ax.set_xlabel('Latitude (degree)',fontsize=16)
        ax.set_ylabel('Pressure (hPa)',fontsize=16)
        ax.set_title(plottitle,fontsize=16,fontweight='bold')
        ax.tick_params(axis="x",labelsize=16)
        ax.tick_params(axis="y",labelsize=16)
#        fig.subplots_adjust(wspace=.4,hspace = .2)
        fig.subplots_adjust(wspace=.4,hspace = .4)

#   fig.savefig('analysis_zmean_bias_{}_{}.{}.png'.format(varname, caseperiod, cyctime))

if __name__ == '__main__':
    plt.ion()
    plot_mean_anomaly(0)
#    plot_zmean_bias()
#    plt.show()
