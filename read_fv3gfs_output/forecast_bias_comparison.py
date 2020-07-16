#!/usr/bin/python
# -*- coding: utf-8 -*-

##################################################################
# Python Program Name: forecast_bias_comparison.py
# Programmer: Ting-Chi Wu @ CIRA/CSU
# Date: 03/24/2020
#
# Program Log: compute bias against GFS Production Analysis
#              with the following three datasets:
#              1) GFS Production Analysis (a)
#              2) Forecasts from FV3GFS cycled experiments (f)
#              3) NCEP-NCAR Reanalysis 1 Climatology (c)
#
# Usage (in ipython): %run forecast_bias_comparison.py varname fhour
# Example: %run forecast_bias_comparison.py TPW 000
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
cpath = '/mnt/ssdatenas/ting-chi/ncep_ncar_r1'
#fhour = '000'

varname = sys.argv[1]
fhour = sys.argv[2]

expnames = ['control_dec2018', 'mhs_dec2018', 'tempestd_dec2018']
#expnames = ['control_dec2018', 'mhs_dec2018', 'tempestd_dec2018_calidiff']
#expnames = ['control_dec2018', 'mhs_dec2018', 'tempestd_dec2018_newcalidiff_retocean']
#expnames = ['control_dec2018', 'mhs_dec2018', 'tempestd_dec2018_newcalidiff_mixqc']
#expnames = ['control_dec2018', 'mhs_dec2018', 'tempestd_dec2018_newcalidiff_cyc1']
#expnames = ['control_dec2018', 'mhs_dec2018', 'tempestd_dec2018_newcalidiff_cyc1noch2']
nexp = len(expnames)

cyctimes = ['2018120812', '2018120900', '2018120912',
            '2018121000', '2018121012', '2018121100',
            '2018121112', '2018121200', '2018121212']
ncyc = len(cyctimes)

#----------------------------------------------------#
nx = 1536
ny = 768

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
#     2. FV3GFS Experiment Forecast Files (C384 grid)
#
for j in range(nexp):
   lfname = 'list_of_{}_files.{}.txt'.format(varname,expnames[j])
   for i in range(ncyc):
      if ( i == 0 ):
         os.system('ls -1 {}/{}/gfs.forecast_derived/gfs.{}/{}.f{}.* > {}'.format(fpath, expnames[j], cyctimes[i], varname, fhour, lfname))
      else:
         os.system('ls -1 {}/{}/gfs.forecast_derived/gfs.{}/{}.f{}.* >> {}'.format(fpath, expnames[j], cyctimes[i], varname, fhour, lfname))
   ffnames = [line.rstrip('\n') for line in open('{}'.format(lfname))]
   os.remove('{}'.format(lfname))
   if j == 0:
      f1fnames = ffnames
   elif j == 1:
      f2fnames = ffnames
   elif j == 2:
      f3fnames = ffnames

nffnames = len(ffnames)
print('Forecast Files (# = {}) :'.format(nffnames))
for jj in range(nffnames):
    print('{}'.format(f1fnames[jj]))
    print('{}'.format(f2fnames[jj]))
    print('{}'.format(f3fnames[jj]))
    print('--------------------------------')
#-----------------------------------------------------------------------------------------#

latfname = '{}/control_dec2018/gfs.forecast_derived/gfs.{}/xlat.dat'.format(fpath,cyctimes[0])
lonfname = '{}/control_dec2018/gfs.forecast_derived/gfs.{}/xlon.dat'.format(fpath,cyctimes[0])
print('Forecast Latitude file: {}'.format(latfname))
print('Forecast Longitude file: {}'.format(latfname))

#
#     3. NCEP-NCAR Reanalysis 1 (climatology) File 
#        (mapped from 2.5 x 2.5 deg to C384 gird)
#
cfname = '{}/{}_December_Climatology.dat'.format(cpath, varname)
print('NCEP-NCAR Re-Analysis File : {}'.format(cfname))

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

# analysis files: 
anal = [readbin(afname, 1, nx, ny) for afname in afnames]

# forecast files:
#fcst = [readbin(ffname, 1, nx, ny) for ffname in ffnames]
fcst1 = [readbin(f1fname, 1, nx, ny) for f1fname in f1fnames]
fcst2 = [readbin(f2fname, 1, nx, ny) for f2fname in f2fnames]
fcst3 = [readbin(f3fname, 1, nx, ny) for f3fname in f3fnames]

# climate files
clim = readbin(cfname, 1, nx, ny)

# lat/lon files
lat = readbin(latfname, 1, nx, ny)
lon = readbin(lonfname, 1, nx, ny)

#----------------------------------------------------#
#   Compute Anomaly Correlation Coefficients
#----------------------------------------------------#

#fmcs = [f - clim for f in fcst]
#f1mcs = [f1 - clim for f1 in fcst1]
#f2mcs = [f2 - clim for f2 in fcst2]
#f3mcs = [f3 - clim for f3 in fcst3]
#amcs = [a - clim for a in anal]

f1mcs = [f1 - a for f1, a in zip(fcst1, anal) ]
f2mcs = [f2 - a for f2, a in zip(fcst2, anal) ]
f3mcs = [f3 - a for f3, a in zip(fcst3, anal) ]

#mfmcs = np.mean(fmcs,axis=0)
mf1mcs = np.mean(f1mcs,axis=0)
mf2mcs = np.mean(f2mcs,axis=0)
mf3mcs = np.mean(f3mcs,axis=0)
#mamcs = np.mean(amcs,axis=0)

zmf1mcs = np.mean(mf1mcs,axis=1)
zmf2mcs = np.mean(mf2mcs,axis=1)
zmf3mcs = np.mean(mf3mcs,axis=1)
#zmamcs = np.mean(mamcs,axis=1)

#plotmax = np.int(max(np.abs(mamcs.min()),np.abs(mamcs.max())))
#plotmin = -plotmax

#----------------------------------------------------#
#   Make a plot of ACC as a function of time 
#----------------------------------------------------#

def plot_mean_anomaly():
  
   print('plot mean bias map')

   fig = figure(num=None, figsize=(16.0, 12.0), dpi=80, facecolor='w', edgecolor='k')
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
      plotcmap = plt.cm.seismic
      if plotnum == 1:
         plotvar = mamcs
         plottitle = 'Analysis - Climatology'

      elif plotnum == 2:
         plotvar = mf1mcs
         plottitle = 'Forecast ({}) - Climatology'.format(expnames[0])

      elif plotnum == 3:
         plotvar = mf2mcs
         plottitle = 'Forecast ({}) - Climatology'.format(expnames[1])

      elif plotnum == 4:
         plotvar = mf3mcs
         plottitle = 'Forecast ({}) - Climatology'.format(expnames[2])

      m.pcolor(xx,yy,plotvar,vmin=plotmin,vmax=plotmax,cmap=plotcmap)
      m.colorbar()
      ax.set_title(plottitle,fontsize=18,fontweight='bold')
      fig.subplots_adjust(wspace=.2,hspace = .1)

def plot_zmean_bias():

   print('plot zonal mean')

   fig = figure(num=None, figsize=(8.0, 12.0), dpi=80, facecolor='w', edgecolor='k')
   ax = fig.add_subplot(1,1,1)
#   plt.plot(zmamcs,lat[:,0],'-k')
   plt.plot(zmf1mcs,lat[:,0],'-b')
   plt.plot(zmf2mcs,lat[:,0],'-g')
   plt.plot(zmf3mcs,lat[:,0],'-r')
#   ax.set_xlim(-2.0,0.5)
#   ax.set_xlim(-3.0,1.0)
   ax.set_xlim(-1.0,1.0)
   plt.grid()
#   plt.legend(['Analysis','{}'.format(expnames[0]),'{}'.format(expnames[1]),'{}'.format(expnames[2])],loc='best')
   plt.legend(['{}'.format(expnames[0]),'{}'.format(expnames[1]),'{}'.format(expnames[2])],loc='best')
   plt.title('Zonal Mean {} Bias (f-a averaged over all cycles): {} h Forecast'.format(varname,fhour))

if __name__ == '__main__':
    plt.ion()
#    plot_mean_anomaly()
    plot_zmean_bias()
    plt.show()
