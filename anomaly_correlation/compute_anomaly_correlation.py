#!/usr/bin/python
# -*- coding: utf-8 -*-

##################################################################
# Python Program Name: compute_anomaly_correlation.py
# Programmer: Ting-Chi Wu @ CIRA/CSU
# Date: 03/04/2020
#
# Program Log: compute anomaly correlation coefficient for the 
#              500 mb Geopotential Height and TPW using the following
#              equation:
# 
#              ACC = mean((f-c)(a-c))/sqrt(mean((f-c)^2)*mean((a-c)^2))
#   
#              with the following three datasets:
#              1) GFS/GDAS Production Analysis (a)
#              2) Forecasts from FV3GFS cycled experiments (f)
#              3) NCEP-NCAR Reanalysis 1 Climatology (c)
#
# Usage (in ipython): %run compute_anomaly_correlation.py expname cyctime varname
# Example: %run compute_anomaly_correlation.py control_dec2018 2018120812 500HGT
#
# Reference: ECMWF_Products_UserGuide_2013.pdf Appendix A-2.6 (owner: Erik Andersson) 
# URL: http://www.atmos.albany.edu/daes/atmclasses/atm401/spring_2015/ECMWF_Products_UserGuide_2013.pdf
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

#cdump = 'gfs'
cdump = 'gdas'
#apath = '/mnt/ssdatenas/ting-chi/gfs_production/keep'
apath = '/mnt/ssdatenas/ting-chi/{}_production'.format(cdump)
fpath = '/mnt/ssdatenas/ting-chi/fv3gfs_experiments'
cpath = '/mnt/ssdatenas/ting-chi/ncep_ncar_r1'

print('icount = {}'.format(sys.argv[1]))
expname = sys.argv[1]
cyctime = sys.argv[2]
varname = sys.argv[3]

#expname = 'control_dec2018'
#expname = 'mhs_dec2018'
#cyctime = '2018120812'
#cyctime = '2018120900'
#cyctime = '2018120912'
#cyctime = '2018121000'
#cyctime = '2018121012'
#cyctime = '2018121100'
#cyctime = '2018121112'
#cyctime = '2018121200'
#cyctime = '2018121212'
#varname = '500HGT'
#varname = 'TPW'

#----------------------------------------------------#
#   Keep as is
#----------------------------------------------------#
nx = 1536
ny = 768

#
#     1. GFS/GDAS Production Analysis Files (mapped from C768 to C384 grid)
#
laname = 'list_of_{}_files.production.txt'.format(varname)
os.system('ls -1 {}/{}.*/*/{}* > {}'.format(apath,cdump,varname, laname))
afnames = [line.rstrip('\n') for line in open('{}'.format(laname))]
os.remove('{}'.format(laname))
nafnames = len(afnames)
print('Production Analysis files (# = {}) :'.format(nafnames))
for jj in range(nafnames):
    print('{}'.format(afnames[jj]))

#
#     2. FV3GFS Experiment Forecast Files (C384 grid)
#
ffolder = '{}/{}/gfs.forecast_derived/gfs.{}'.format(fpath,expname,cyctime)
lfname = 'list_of_{}_files.{}.txt'.format(varname,expname)
os.system('ls -1 {}/{}* > {}'.format(ffolder, varname, lfname))
ffnames = [line.rstrip('\n') for line in open('{}'.format(lfname))]
os.remove('{}'.format(lfname))
nffnames = len(ffnames)
print('Forecast Files (# = {}) :'.format(nffnames))
for jj in range(nffnames):
    print('{}'.format(ffnames[jj]))
#-----------------------------------------------------------------------------------------#

if (cyctime[:4]=='2018'):
   latfname = '{}/control_dec2018/gfs.forecast_derived/gfs.{}/xlat.dat'.format(fpath,cyctime)
   lonfname = '{}/control_dec2018/gfs.forecast_derived/gfs.{}/xlon.dat'.format(fpath,cyctime)
elif (cyctime[:4]=='2019'):
   latfname = '{}/control_may2019/gfs.forecast_derived/gfs.{}/xlat.dat'.format(fpath,cyctime)
   lonfname = '{}/control_may2019/gfs.forecast_derived/gfs.{}/xlon.dat'.format(fpath,cyctime)
print('Forecast Latitude file: {}'.format(latfname))
print('Forecast Longitude file: {}'.format(latfname))

#
#     3. NCEP-NCAR Reanalysis 1 (climatology) File 
#        (mapped from 2.5 x 2.5 deg to C384 gird)
#
if (cyctime[:4]=='2018'):
   cfname = '{}/{}_December_climatology.dat'.format(cpath, varname)
elif (cyctime[:4]=='2019'):
   cfname = '{}/{}_May_climatology.dat'.format(cpath, varname)

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
# First, find the begin and the end index from the
# analysis file list (use cyctime and forecast file length)

for idx in range(len(afnames)):
#    rtime = afnames[idx][-24:-16]+afnames[idx][-15:-13]
#    rtime = afnames[idx][43:51]+afnames[idx][52:54]
#    rtime = afnames[idx][48:56]+afnames[idx][57:59]
    rtime = afnames[idx][45:53]+afnames[idx][54:56]
#    print('idx = {}, rtime = {}'.format(idx, rtime))
    if ( rtime == cyctime ):
       idx_begin = idx
idx_end = idx_begin + nffnames

print('Production Analysis idx_begin = {}'.format(idx_begin))
print('Production Analysis idx_end = {}'.format(idx_end))
for jj in range(idx_begin,idx_end):
    print('{}'.format(afnames[jj]))

anal = [readbin(afname, 1, nx, ny) for afname in afnames[idx_begin:idx_end]]
#anal = readbin(afname, 1, nx, ny)

# forecast files:
fcst = [readbin(ffname, 1, nx, ny) for ffname in ffnames]
#fcst = readbin(ffnames[ifile], 1, nx, ny)

# climate files
clim = readbin(cfname, 1, nx, ny)

#----------------------------------------------------#
#   Index for Regions
#----------------------------------------------------#

lat = readbin(latfname, 1, nx, ny)
lon = readbin(lonfname, 1, nx, ny)

regions = ['Global', 'North Hemisphere', 'South Hemisphere', 'Tropics']
nregions = len(regions)
rdx_begin = np.zeros([nregions])
rdx_end = np.zeros([nregions])

# global 
rdx_begin[0] = np.int(0)
rdx_end[0] = np.int(ny)

# northern hemisphere
rdx_begin[1] = np.int(find_idx_closest_value(lat[:,0],20.0))
rdx_end[1] = np.int(find_idx_closest_value(lat[:,0],80.0))

# southern hemisphere
rdx_begin[2] = np.int(find_idx_closest_value(lat[:,0],-80.0))
rdx_end[2] = np.int(find_idx_closest_value(lat[:,0],-20.0))

# tropics
rdx_begin[3] = np.int(find_idx_closest_value(lat[:,0],-20.0) + 1)
rdx_end[3] = np.int(find_idx_closest_value(lat[:,0],20.0) - 1)
 
rdx_begin = rdx_begin.astype(int)
rdx_end = rdx_end.astype(int)
print('Index to identify regions:')
for jj in range(nregions):
    print('{}: begins at {}, ends at {}'.format(regions[jj],rdx_begin[jj],rdx_end[jj]))

#----------------------------------------------------#
#   Compute Anomaly Correlation Coefficients
#----------------------------------------------------#

fmcs = [f - clim for f in fcst]
amcs = [a - clim for a in anal]

#acc = np.zeros([nffnames])
#for ii in range(nffnames):
#    acc_top = np.mean((fmcs[ii] * amcs[ii]).flatten())
#    acc_bl = np.mean((fmcs[ii] * fmcs[ii]).flatten())
#    acc_br = np.mean((amcs[ii] * amcs[ii]).flatten())
#    acc_bottom = np.sqrt(acc_bl * acc_br)  
#    acc[ii] = acc_top / acc_bottom
#    print('acc = {}'.format(acc[ii]))

acc = np.zeros([nffnames,nregions])
for ii in range(nffnames):
    for rr in range(nregions):
        fmcxamc = fmcs[ii][rdx_begin[rr]:rdx_end[rr],:] * amcs[ii][rdx_begin[rr]:rdx_end[rr],:]
        fmcxfmc = fmcs[ii][rdx_begin[rr]:rdx_end[rr],:] * fmcs[ii][rdx_begin[rr]:rdx_end[rr],:]
        amcxamc = amcs[ii][rdx_begin[rr]:rdx_end[rr],:] * amcs[ii][rdx_begin[rr]:rdx_end[rr],:]
        acc_top = np.mean(fmcxamc.flatten())
        acc_bl = np.mean(fmcxfmc.flatten())
        acc_br = np.mean(amcxamc.flatten())
        acc_bottom = np.sqrt(acc_bl * acc_br)  
        acc[ii,rr] = acc_top / acc_bottom
    print('acc = {}, {}, {}, {}'.format(acc[ii,0],acc[ii,1],acc[ii,2],acc[ii,3]))

outfname = 'acc.{}_analysis/acc_{}.{}.{}.txt'.format(cdump, varname, expname,cyctime)
if (path.exists(outfname)):
   os.remove('{}'.format(outfname))
with open (outfname,'a') as text_file:
   for ii in range(len(acc)):
       text_file.write('{}, {}, {}, {}\n'.format(acc[ii,0],acc[ii,1],acc[ii,2],acc[ii,3]))
text_file.close()

#----------------------------------------------------#
#   Make a plot of ACC as a function of time 
#----------------------------------------------------#

def plot_acc_vs_time():

   fig1 = plt.figure(num=None, figsize=(16.0, 8.0), dpi=80, facecolor='w', edgecolor='k')
   xx = np.linspace(0,240,nffnames)
   xx_str = [np.str(np.int(x)) for x in xx]
   ax = fig1.add_subplot(1,1,1)
   plt.plot(xx,acc[:,0],'-k')
   plt.plot(xx,acc[:,1],'-b')
   plt.plot(xx,acc[:,2],'-g')
   plt.plot(xx,acc[:,3],'-r')
   plt.legend(['Global','Northern Hemisphere','Southern Hemisphere','The Tropics'],loc='best',fontsize=16)
   plt.xlabel('Forecast Lead Time (hour)',fontsize=12)
   plt.ylabel('ACC',fontsize=12)
   plt.title('Anomaly Correlation Coefficient from {} cycle {}: {}'.format(expname, cyctime, varname),fontsize=20)
   #ax.set_xlim(0,240)
   plt.xticks(xx,xx_str,fontsize=12)
   plt.tight_layout()
   plt.grid()
   #plt.savefig('figs/acc_{}.{}.{}.png'.format(varname,expname,cyctime))


#----------------------------------------------------#
#   Make a plot of Hovemuller diagram of fmc vs. amc 
#----------------------------------------------------#

def plot_hovmuller():

   fmcs_zonal = np.mean(fmcs[:],axis=2)
   amcs_zonal = np.mean(amcs[:],axis=2)

   xx_hov = np.zeros([nffnames,ny])
   yy_hov = np.zeros([nffnames,ny])

   for jj in range(ny):
       yy_hov[:,jj] = xx[:]
   for ii in range(nffnames):
       xx_hov[ii,:] = lat[:,0]

   if (varname == '500HGT'):
      varmax = 250.0
   else:
      varmax = 10.0
   varmin = -varmax

   fig2 = plt.figure(num=None, figsize=(16.0, 10.0), dpi=80, facecolor='w', edgecolor='k')
   fig2.suptitle('Hovmoller Diagram of Zonal-Averaged {} Anomaly'.format(varname),fontsize=15,fontweight='bold')
   ax = fig2.add_subplot(1,2,1)
   plt.pcolor(xx_hov,yy_hov,fmcs_zonal,vmin=varmin,vmax=varmax,cmap=plt.cm.seismic)
   plt.colorbar()
   plt.xlabel('Latitude',fontsize=12)
   plt.ylabel('Forecast Lead Time (hour)',fontsize=12)
   plt.title('Forecast ({}: {}) - Climatology'.format(expname,cyctime))
   plt.yticks(xx,xx_str,fontsize=12)
   ax = fig2.add_subplot(1,2,2)
   plt.pcolor(xx_hov,yy_hov,amcs_zonal,vmin=varmin,vmax=varmax,cmap=plt.cm.seismic)
   plt.colorbar()
   plt.xlabel('Latitude',fontsize=12)
   plt.ylabel('Forecast Lead Time (hour)',fontsize=12)
   plt.title('Analysis - Climatology')
   plt.yticks(xx,xx_str,fontsize=12)

#if __name__ == '__main__':
#    plt.ion()
#    plot_acc_vs_time()
#    plot_hovmuller()
