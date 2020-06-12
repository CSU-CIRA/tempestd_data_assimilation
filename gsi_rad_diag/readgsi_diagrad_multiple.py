#!/usr/bin/python
# -*- coding: utf-8 -*-

#####################################################
# Python Program Name: readgsi_diagrad.py
# Programmer: Ting-Chi Wu @ CIRA/CSU
# Date: 05/22/2019
# Program Log: read ascii diag_tempest_cubesat_ges
#                  and/or diag_tempest_cubesat_anl                      
#              from gsi (after read_diag_rad.exe)
#####################################################

#from __future__ import print_function
import numpy as np
#import netCDF4
import matplotlib as mpl
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import norm
import matplotlib.gridspec as gridspec
from mpl_toolkits.basemap import Basemap

#####################################################
# User Input (Time, Path, Files, etc)
#####################################################

#path='/data/users/twu/SSDAE/Hera/control_dec2018/2018120812'
parent_path = '/mnt/ssdatenas/ting-chi/fv3gfs_experiments'
#expname = 'control_dec2018'
#expname = 'control_may2019'
#expname = 'mhs_dec2018'
#expname = 'tempestd_dec2018'
#expname = 'tempestd_dec2018_calidiff'
#expname = 'tempestd_dec2018_retocean'
#expname = 'tempestd_dec2018_newcalidiff_retocean'
#expname = 'tempestd_dec2018_newcalidiff_mixqc'
expname = 'tempestd_dec2018_newcalidiff_cyc1'
#expname = 'tempestd_dec2018_newcalidiff_cyc1noch2'
#expname = 'tempestd_dec2018_mhsqc'

#expname = 'mhs_may2019'
#expname = 'tempestd_may2019_newcalidiff_cyc1'
#expname = 'tempestd_may2019_newcalidiff_cyc1noch2'

path = '{}/{}/gfs.diag_rad'.format(parent_path,expname)

#datetime = 2018120812
#datetime = 2018120900
#datetime = 2018120912
#datetime = 2018121000
#datetime = 2018121012
#datetime = 2018121100
#datetime = 2018121112
#datetime = 2018121200
#datetime = 2018121212

#datetime = 2019051212
#datetime = 2019051300
#datetime = 2019051312
#datetime = 2019051400
#datetime = 2019051412
#datetime = 2019051500
#datetime = 2019051512

datetime = 'all'
#datetime = '00z'
#datetime = '12z'


#sensats = ['mhs_metop-a','mhs_metop-b','mhs_n19', 'tempest_cubesat']
sensats = ['mhs_n19','mhs_metop-b','mhs_metop-a', 'tempest_cubesat']
nchanl = 5
tbdiff = 5.0

format_dict = {'path': path,  
               'fname1': sensats[0],
               'fname2': sensats[1],
               'fname3': sensats[2],
               'fname4': sensats[3], 
               'cyctime': datetime}

titles = ['MHS NOAA19', 'MHS NOAA19',
          'MHS MetOp-B', 'MHS MetOp-B',
          'MHS MetOp-A', 'MHS MetOp-A',
          'TEMPEST-D CubeSat', 'TEMPEST-D CubeSat']
#fnames = ['{path}/results_{fname1}_ges.{cyctime}',
#          '{path}/results_{fname1}_ges.{cyctime}',
#          '{path}/results_{fname2}_ges.{cyctime}',
#          '{path}/results_{fname2}_ges.{cyctime}',
#          '{path}/results_{fname3}_ges.{cyctime}',
#          '{path}/results_{fname3}_ges.{cyctime}',
#          '{path}/results_{fname4}_ges.{cyctime}',
#          '{path}/results_{fname4}_ges.{cyctime}']
fnames = ['{path}/results_{fname1}_ges.{cyctime}',
          '{path}/results_{fname1}_anl.{cyctime}',
          '{path}/results_{fname2}_ges.{cyctime}',
          '{path}/results_{fname2}_anl.{cyctime}',
          '{path}/results_{fname3}_ges.{cyctime}',
          '{path}/results_{fname3}_anl.{cyctime}',
          '{path}/results_{fname4}_ges.{cyctime}',
          '{path}/results_{fname4}_anl.{cyctime}']

fnames = [fname.format(**format_dict) for fname in fnames]
for title, fname in zip(titles, fnames):
    print('{}: {}'.format(title, fname))


#####################################################
# Define Constants
#####################################################

gas_const = 287.0 # J/K/kg
one_rad = 180/np.pi
freq_tempestd = ['87 GHz', '164 GHz', '174 GHz', '178 GHz', '181 GHz']
#tbmax = [280., 290., 290., 280., 290.]
#tbmin = [200., 240., 240., 240., 240.]
tbmax = [295., 295., 295., 295., 295.]
tbmin = [190., 190., 190., 190., 190.]
freq_mhs = ['89 GHz', '157 GHz', '183.31$\pm$1 GHz', '183.31$\pm$3 GHz', '190.31 GHz']

#####################################################
# Define Function that Read
#####################################################

def read_diag_rad(fname,nchanl):
    '''
    read the ascii file daig_tempest_cubesat_ges/anl
    line by line and output into an array allvar
    allvar[ichanl,idat]
    idat = 0: channel number
    idat = 1: latitude (degree)
    idat = 2: longitude (degree)
    idat = 3: time relative to analysis (hour)
    idat = 4: ocean/land indicator
    idat = 5: sensor scan position
    idat = 6: observed brightness temperature (K)
    idat = 7: o-b brightness temperature w/ bias correction (K)
    idat = 8: o-b brightness temperature w/o bias correction (K)
    idat = 9: inverse observation error (K)
    idat = 10: quality flag (see getchqc for more info)
    idat = 11: angle bias-correction term
    idat = 12: IWP retrieval (CSU 1DVAR)
    idat = 13: LWP retrieval (CSU 1DVAR)
    idat = 14: IWP guess (FV3GFS)
    idat = 15: LWP guess (FV3GFS)
    '''
    fid = open(fname,'r') 
    print('Open to read {}'.format(fname))
    # varlist 0 is the character "channel"
    # varlist 2 is the ":" sign
    varlist = [1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]
    nvar = len(varlist)
    var=[[] for i in range(nvar)]
    count=0
    for line in fid.readlines(): # this reads one line at a time
#        print line
        spl_line = line.split() # split into a list of strings by whitespace
        if count >1 and count < 2+nchanl:
           print spl_line
        elif count >= 2+nchanl:
           for ivar in range(nvar):
               var[ivar].append(float(spl_line[varlist[ivar]])) # convert from string to float and append
        count+=1 # same as count=count+1
    var = np.array(var)
    count = count - (2+nchanl) # take away count for headers
    print('count = {}'.format(count))
    nobs = count/nchanl
    print('number of obs (per channel) = {}'.format(nobs))
    allvar = np.zeros([nchanl, nvar, nobs])
    for j in range(nvar):
        for i in range(nchanl):
            allvar[i,j,:] = var[j,i::nchanl]
            loc = np.where(allvar[i,2,:]>180.)
            allvar[i,2,loc] = allvar[i,2,loc] - 360.
#    return(var, allvar)
    return(allvar)

def nqc(QCarray,ich,nobs,qcid):
#   if qcid >= 0 :
#     x=[i for i in range(nobs) if QCarray[ich-1,i] == qcid]
#   elif qcid < 0:
#     x=[i for i in range(nobs) if QCarray[ich-1,i] < 0]
   x=[i for i in range(nobs) if QCarray[ich-1,i] == qcid]
   return len(x)

def getchqc(QCarray,ich,nobs):
   chqc=np.empty([13])
   chqc[0]=nqc(QCarray,ich,nobs,0) # igood_qc = 0
#  chqc[1]=nqc(QCarray,ich,nobs,1)
   chqc[1]=nqc(QCarray,ich,nobs,1) # ifail_satinfo_qc = 1 and set iuse = -1
   chqc[2]=nqc(QCarray,ich,nobs,2) # ifail_crtm_qc = 2
   chqc[3]=nqc(QCarray,ich,nobs,3) # ifail_gross_qc = 3
   chqc[4]=nqc(QCarray,ich,nobs,4) # ifail_interchan_qc = 4
   chqc[5]=nqc(QCarray,ich,nobs,6) # ifail_gross_routine_qc = 6
   chqc[6]=nqc(QCarray,ich,nobs,7) # ifail_cloud_qc = 7
   chqc[7]=nqc(QCarray,ich,nobs,8) # ifail_emiss_qc = 8
   chqc[8]=nqc(QCarray,ich,nobs,9) # ifail_range_qc = 9
   chqc[9]=nqc(QCarray,ich,nobs,50) # ifail_fact1_qc = 50 for mhs (TPW index)
   chqc[10]=nqc(QCarray,ich,nobs,51) # ifail_factch4_qc = 51
   chqc[11]=nqc(QCarray,ich,nobs,61) # ifail_iland_det = 61
   chqc[12]=chqc[0]+chqc[1]+chqc[2]+chqc[3]+chqc[4]+chqc[5]+chqc[6]+chqc[7]+chqc[8]+chqc[9]+chqc[10]+chqc[11]
   return chqc

#####################################################
# Call Function and Manipulate Data
#####################################################


diags = [read_diag_rad(fname,nchanl) for fname in fnames]

diag1 = diags[0] # GES: MHS NOAA 19
diag2 = diags[1] # GES: MHS NOAA 19
diag3 = diags[2] # GES: MHS MetOp-B
diag4 = diags[3] # GES: MHS MetOp-B
diag5 = diags[4] # GES: MHS MetOp-A
diag6 = diags[5] # GES: MHS MetOP-A
diag7 = diags[6] # GES: TEMPEST-D
diag8 = diags[7] # GES: TEMPEST-D

qc1 = diag1[:,10,:]
qc2 = diag2[:,10,:]
qc3 = diag3[:,10,:]
qc4 = diag4[:,10,:]
qc5 = diag5[:,10,:]
qc6 = diag6[:,10,:]
qc7 = diag7[:,10,:]
qc8 = diag8[:,10,:]
for ich in range(0,nchanl):
    qc_out1 = getchqc(qc1,ich+1,qc1.shape[1])
    qc_out2 = getchqc(qc2,ich+1,qc2.shape[1])
    qc_out3 = getchqc(qc3,ich+1,qc3.shape[1])
    qc_out4 = getchqc(qc4,ich+1,qc4.shape[1])
    qc_out5 = getchqc(qc5,ich+1,qc5.shape[1])
    qc_out6 = getchqc(qc6,ich+1,qc6.shape[1])
    qc_out7 = getchqc(qc7,ich+1,qc7.shape[1])
    qc_out8 = getchqc(qc8,ich+1,qc8.shape[1])
    print('Channel {} QC: good, satinfo, crtm, gross, interch, gross_mhs, cloud, emiss, n/a, fact1(tpw), n/a, land, total'.format(ich+1))
    print('{} : {}'.format(titles[0],qc_out1))
    print('{} : {}'.format(titles[1],qc_out2))
    print('{} : {}'.format(titles[2],qc_out3))
    print('{} : {}'.format(titles[3],qc_out4))
    print('{} : {}'.format(titles[4],qc_out5))
    print('{} : {}'.format(titles[5],qc_out6))
    print('{} : {}'.format(titles[6],qc_out7))
    print('{} : {}'.format(titles[7],qc_out8))

#biases = np.zeros([5])
#gloc_use = np.where(diags[6][ichanl,10,:]==)
#aoc_use = np.where(diags[7][ichanl,10,:]==)
#diags[6][ichanl,7,gloc_use].shape[1]


#####################################################
# Make Plots
#####################################################

def plot_omb_vs_oma(ichanl):
    '''
    make a 4x1 sub-panel figure to show o-b vs o-a for all 4 instruments
    '''
    fig = plt.figure(num=None, figsize=(24.0, 16.0), dpi=80, facecolor='w', edgecolor='k')
    fig.suptitle('FV3 GSI Global Analysis from {} experiment Valid at {}'.format(expname, datetime),fontsize=15,fontweight='bold')
    for iplot in range(4):
        plotnum = iplot + 1
        ax = fig.add_subplot(4,1,plotnum)
#        print('ichanl = {}'.format(ichanl))
#        print('iplot = {}'.format(iplot))
        loc_use = np.where(diags[iplot*2][ichanl,10,:]==0)
        nnb = diags[iplot*2][ichanl,7,loc_use].shape[1]
        nna = diags[iplot*2+1][ichanl,7,loc_use].shape[1]
        xx = np.linspace(1,nnb,nnb)
        plt.plot(xx,diags[iplot*2][ichanl,7,loc_use][0,:], color='blue')
        plt.plot(xx,diags[iplot*2+1][ichanl,7,loc_use][0,:], color='red')
        plt.legend(['O-B','O-A'],loc='upper right',fontsize=12)
        if plotnum < 4:
           ax.set_title('{} Ch {}: {} (# of Obs = {},{})'.format(titles[iplot*2][4:],ichanl+1,freq_mhs[ichanl],nnb, nna),fontsize=20,fontweight='bold')
        else:
           ax.set_title('{} Ch {}: {} (# of Obs = {},{})'.format(titles[iplot*2][4:],ichanl+1,freq_tempestd[ichanl],nnb, nna),fontsize=20,fontweight='bold')
        fig.subplots_adjust(wspace=.25)
        fig.subplots_adjust(hspace=.5)
        plt.grid()
#    fig.savefig('{}.{}.omb_oma_serial.ch{}.png'.format(expname,datetime,ichanl+1))

def plot_omb_vs_oma_hist(ichanl):
    '''
    similar to plot_omb_vs_oma,
    but instead of a serial plot, it's a histogram for all 4 instruments
    '''
    num_bins = 100
    fig = plt.figure(num=None, figsize=(12.0, 24.0), dpi=80, facecolor='w', edgecolor='k')
    fig.suptitle('FV3 GSI Global Analysis from {} experiment Valid at {}'.format(expname, datetime),fontsize=15,fontweight='bold')
    for iplot in range(4):
        plotnum = iplot + 1
        ax = fig.add_subplot(4,1,plotnum)
        loc_use = np.where(diags[iplot*2][ichanl,10,:]==0)
        inv_err = diags[iplot*2][ichanl,9,loc_use][0,:]
#        ghist = diags[iplot*2][ichanl,7,loc_use][0,:]
        ghist = diags[iplot*2][ichanl,7,loc_use][0,:]*inv_err
#        n, bins, patches = plt.hist(ghist, num_bins, facecolor='blue',alpha=0.5)
        bins, edges = np.histogram(ghist, num_bins)
        left, right = edges[:-1], edges[1:]
        X = np.array([left,right]).T.flatten()
        Y = np.array([bins,bins]).T.flatten()
        plt.plot(X,Y,'-b')
        ghist_nbc = diags[iplot*2][ichanl,8,loc_use][0,:]*inv_err
        gbias_before = np.mean(ghist_nbc)
        gbias_after = np.mean(ghist)
        gstd_before = np.std(ghist_nbc)
        gstd_after = np.std(ghist)

        loc_use = np.where(diags[iplot*2+1][ichanl,10,:]==0)
        inv_err = diags[iplot*2+1][ichanl,9,loc_use][0,:]
#        ahist = diags[iplot*2+1][ichanl,7,loc_use][0,:]
        ahist = diags[iplot*2+1][ichanl,7,loc_use][0,:]*inv_err
#        n, bins, patches = plt.hist(ahist, num_bins, facecolor='red',alpha=0.5)
        bins, edges = np.histogram(ahist, num_bins)
        left, right = edges[:-1], edges[1:]
        X = np.array([left,right]).T.flatten()
        Y = np.array([bins,bins]).T.flatten()
        plt.plot(X,Y,'-r')
        ahist_nbc = diags[iplot*2+1][ichanl,8,loc_use][0,:]*inv_err
        abias_before = np.mean(ahist_nbc)
        abias_after = np.mean(ahist)
        astd_before = np.std(ahist_nbc)
        astd_after = np.std(ahist)

#        plt.legend(['O-B','O-A'],loc='best',fontsize=20)
        plt.legend(['O-B ($\mu={}, \sigma={}$)'.format(round(gbias_after,2),round(gstd_after,2)),
                    'O-A ($\mu={}, \sigma={}$)'.format(round(abias_after,2),round(astd_after,2))],
                    loc='best',fontsize=16)
        ax.tick_params(axis='x', labelsize=16)
        ax.tick_params(axis='y', labelsize=16)
        ax.set_xlabel('O-B/O-A (normalized by obs error [K])',fontsize=16)
        ax.set_ylabel('Obs Count',fontsize=16)
        ax.set_xlim(-3.0,3.0)
        if plotnum < 4:
           ax.set_title('{} Ch {}: {} (n = {}, {})'.format(titles[iplot*2][4:],ichanl+1,freq_mhs[ichanl],len(ghist),len(ahist)),fontsize=20,fontweight='bold')
           print('{} Ch {}: {} Ges: Bias before/after BC = {}/{}'.format(titles[iplot*2][4:],ichanl+1,freq_mhs[ichanl],gbias_before,gbias_after))
           print('{} Ch {}: {} Anl: Bias before/after BC = {}/{}'.format(titles[iplot*2][4:],ichanl+1,freq_mhs[ichanl],abias_before,abias_after))
        else:
           ax.set_title('{} Ch {}: {} (n = {}, {})'.format(titles[iplot*2][4:],ichanl+1,freq_tempestd[ichanl],len(ghist),len(ahist)),fontsize=20,fontweight='bold')
           print('{} Ch {}: {} Ges: Bias before/after BC = {}/{}'.format(titles[iplot*2][4:],ichanl+1,freq_tempestd[ichanl],gbias_before,gbias_after))
           print('{} Ch {}: {} Anl: Bias before/after BC = {}/{}'.format(titles[iplot*2][4:],ichanl+1,freq_tempestd[ichanl],abias_before,abias_after))
        fig.subplots_adjust(wspace=.25)
        fig.subplots_adjust(hspace=.5)
        fig.subplots_adjust(bottom=0.15)
        plt.grid()
#    fig.savefig('{}.{}.omb_oma_histogram.ch{}.png'.format(expname,datetime,ichanl+1))

def plot_omb_hist(ichanl):
    '''
    similar to plot_omb_vs_oma_hist,
    but focus o-b before (b) and after (a) bias corrections
    '''
    num_bins = 100
    fig = plt.figure(num=None, figsize=(12.0, 24.0), dpi=80, facecolor='w', edgecolor='k')
    fig.suptitle('FV3 GSI Global Analysis from {} experiment Valid at {}'.format(expname, datetime),fontsize=15,fontweight='bold')
    for iplot in range(4):
        plotnum = iplot + 1
        ax = fig.add_subplot(4,1,plotnum)
        loc_use = np.where(diags[iplot*2][ichanl,10,:]==0)
        inv_err = diags[iplot*2][ichanl,9,loc_use][0,:]
#        bhist = diags[iplot*2][ichanl,8,loc_use][0,:]
        bhist = diags[iplot*2][ichanl,8,loc_use][0,:]*inv_err
        bins, edges = np.histogram(bhist, num_bins)
        left, right = edges[:-1], edges[1:]
        X = np.array([left,right]).T.flatten()
        Y = np.array([bins,bins]).T.flatten()
        plt.plot(X,Y,'-b')
        bbias = np.mean(bhist)
        bstd = np.std(bhist)

#        ahist = diags[iplot*2][ichanl,7,loc_use][0,:]
        ahist = diags[iplot*2][ichanl,7,loc_use][0,:]*inv_err
        bins, edges = np.histogram(ahist, num_bins)
        left, right = edges[:-1], edges[1:]
        X = np.array([left,right]).T.flatten()
        Y = np.array([bins,bins]).T.flatten()
        plt.plot(X,Y,'-r')
        abias = np.mean(ahist)
        astd = np.std(ahist)

#        plt.legend(['O-B','O-A'],loc='best',fontsize=20)
        plt.legend(['Before BC ($\mu={}, \sigma={}$)'.format(round(bbias,2),round(bstd,2)),
                    'After BC ($\mu={}, \sigma={}$)'.format(round(abias,2),round(astd,2))],
                    loc='best',fontsize=16)
        ax.tick_params(axis='x', labelsize=16)
        ax.tick_params(axis='y', labelsize=16)
        ax.set_xlabel('O-B (normalized by obs error [K])',fontsize=16)
        ax.set_ylabel('Obs Count',fontsize=16)
        ax.set_xlim(-3.0,3.0)
        if plotnum < 4:
           ax.set_title('{} Ch {}: {} (n = {}, {})'.format(titles[iplot*2][4:],ichanl+1,freq_mhs[ichanl],len(bhist),len(ahist)),fontsize=20,fontweight='bold')
           print('{} Ch {}: {} Ges: Bias before BC = {}'.format(titles[iplot*2][4:],ichanl+1,freq_mhs[ichanl],bbias))
           print('{} Ch {}: {} Ges: Bias after BC = {}'.format(titles[iplot*2][4:],ichanl+1,freq_mhs[ichanl],abias))
        else:
           ax.set_title('{} Ch {}: {} (n = {}, {})'.format(titles[iplot*2][4:],ichanl+1,freq_tempestd[ichanl],len(bhist),len(ahist)),fontsize=20,fontweight='bold')
           print('{} Ch {}: {} Ges: Bias before BC = {}'.format(titles[iplot*2][4:],ichanl+1,freq_tempestd[ichanl],bbias))
           print('{} Ch {}: {} Ges: Bias after BC = {}'.format(titles[iplot*2][4:],ichanl+1,freq_tempestd[ichanl],abias))
        fig.subplots_adjust(wspace=.25)
        fig.subplots_adjust(hspace=.5)
        fig.subplots_adjust(bottom=0.15)
        plt.grid()
    fig.savefig('{}.{}.omb_bc_histogram.ch{}.png'.format(expname,datetime,ichanl+1))

def plot_oma_hist(ichanl):
    '''
    similar to plot_oma_hist,
    but instead of o-b, it's o-a
    '''
    num_bins = 100
    fig = plt.figure(num=None, figsize=(12.0, 24.0), dpi=80, facecolor='w', edgecolor='k')
    fig.suptitle('FV3 GSI Global Analysis from {} experiment Valid at {}'.format(expname, datetime),fontsize=15,fontweight='bold')
    for iplot in range(4):
        plotnum = iplot + 1
        ax = fig.add_subplot(4,1,plotnum)
        loc_use = np.where(diags[iplot*2+1][ichanl,10,:]==0)
        inv_err = diags[iplot*2+1][ichanl,9,loc_use][0,:]
#        bhist = diags[iplot*2+1][ichanl,8,loc_use][0,:]
        bhist = diags[iplot*2+1][ichanl,8,loc_use][0,:]*inv_err
        bins, edges = np.histogram(bhist, num_bins)
        left, right = edges[:-1], edges[1:]
        X = np.array([left,right]).T.flatten()
        Y = np.array([bins,bins]).T.flatten()
        plt.plot(X,Y,'-b')
        bbias = np.mean(bhist)
        bstd = np.std(bhist)

#        ahist = diags[iplot*2+1][ichanl,7,loc_use][0,:]
        ahist = diags[iplot*2+1][ichanl,7,loc_use][0,:]*inv_err
        bins, edges = np.histogram(ahist, num_bins)
        left, right = edges[:-1], edges[1:]
        X = np.array([left,right]).T.flatten()
        Y = np.array([bins,bins]).T.flatten()
        plt.plot(X,Y,'-r')
        abias = np.mean(ahist)
        astd = np.std(ahist)

        plt.legend(['Before BC ($\mu={}, \sigma={}$)'.format(round(bbias,2),round(bstd,2)),
                    'After BC ($\mu={}, \sigma={}$)'.format(round(abias,2),round(astd,2))],
                    loc='best',fontsize=16)
        ax.tick_params(axis='x', labelsize=16)
        ax.tick_params(axis='y', labelsize=16)
        ax.set_xlabel('O-A (normalized by obs error [K])',fontsize=16)
        ax.set_ylabel('Obs Count',fontsize=16)
        ax.set_xlim(-3.0,3.0)
        if plotnum < 4:
           ax.set_title('{} Ch {}: {} (n = {}, {})'.format(titles[iplot*2][4:],ichanl+1,freq_mhs[ichanl],len(bhist),len(ahist)),fontsize=20,fontweight='bold')
           print('{} Ch {}: {} Anl: Bias before BC = {}'.format(titles[iplot*2][4:],ichanl+1,freq_mhs[ichanl],bbias))
           print('{} Ch {}: {} Anl: Bias after BC = {}'.format(titles[iplot*2][4:],ichanl+1,freq_mhs[ichanl],abias))
        else:
           ax.set_title('{} Ch {}: {} (n = {}, {})'.format(titles[iplot*2][4:],ichanl+1,freq_tempestd[ichanl],len(bhist),len(ahist)),fontsize=20,fontweight='bold')
           print('{} Ch {}: {} Anl: Bias before BC = {}'.format(titles[iplot*2][4:],ichanl+1,freq_tempestd[ichanl],bbias))
           print('{} Ch {}: {} Anl: Bias after BC = {}'.format(titles[iplot*2][4:],ichanl+1,freq_tempestd[ichanl],abias))
        fig.subplots_adjust(wspace=.25)
        fig.subplots_adjust(hspace=.5)
        fig.subplots_adjust(bottom=0.15)
        plt.grid()
    fig.savefig('{}.{}.oma_bc_histogram.ch{}.png'.format(expname,datetime,ichanl+1))

def plot_orbits(ichanl):
    '''
    make a total of 8 figures to show obs versus analysis for all 4 instruments
    '''
    fig = plt.figure(num=None, figsize=(32.0, 16.0), dpi=80, facecolor='w', edgecolor='k')
    fig.suptitle('FV3 GSI Global Analysis from {} experiment: Valid at {}: Ch {}'.format(expname, datetime,ichanl+1),fontsize=15,fontweight='bold')

    for iplot in range(8):
#    for iplot in range(6,8):
        plotnum = iplot + 1
#        plotnum = iplot - 6  + 1
        ax = fig.add_subplot(4,2,plotnum)
#        ax = fig.add_subplot(1,2,plotnum)

#        fig = plt.figure(num=None, figsize=(16.0, 8.0), dpi=80, facecolor='w', edgecolor='k')
#        ax = fig.add_subplot(1,1,1)
        m = Basemap(projection='mill',llcrnrlat=-90.,urcrnrlat=90.,
                    llcrnrlon=-180.,urcrnrlon=180.0,resolution='l')
#        latdeg = 20.0
#        m = Basemap(projection='mill',llcrnrlat=-latdeg,urcrnrlat=latdeg,
#                    llcrnrlon=-180.,urcrnrlon=180.0,resolution='l')
#        m = Basemap(projection='mill',llcrnrlat=-25.,urcrnrlat=-5.,
#                    llcrnrlon=130.,urcrnrlon=170.0,resolution='l')

        m.drawcoastlines()
        m.drawcountries()
        m.drawstates(linewidth=1,linestyle='solid',color='k')
        # draw parallels.
        parallels = np.arange(-90.,90,10.)
        m.drawparallels(parallels,labels=[1,0,0,0],fontsize=16)
        # draw meridians
        meridians = np.arange(-180.,180.,30.)
        m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=16)

#        loc_use = np.where(diags[iplot][ichanl,10,:]==0)
#        loc_use = np.where((diags[iplot][ichanl,10,:]==0) & (diags[iplot][ichanl,1,:]<=latdeg) & (diags[iplot][ichanl,1,:]>=-latdeg))
#       good qc, over ocean, within N/S 20 degrees

#        if (np.mod(plotnum,2) == 0):
#           lonmax = 0.0
#           lonmin = -120.0
#        else:
#           lonmax = 180.0
#           lonmin = 80.0
#        loc_use = np.where((diags[iplot][ichanl,10,:]==0) & (diags[iplot][ichanl,4,:]==1.0) & (diags[iplot][ichanl,1,:]<=latdeg) & (diags[iplot][ichanl,1,:]>=-latdeg) \
#                           & (diags[iplot][ichanl,2]>=lonmin) & (diags[iplot][ichanl,2]<=lonmax))
#        loc_use = np.where((diags[iplot][ichanl,10,:]==0) & (diags[iplot][ichanl,4,:]==1.0) \
#                           & (diags[iplot][ichanl,2]>=lonmin) & (diags[iplot][ichanl,2]<=lonmax))
        loc_all = np.where(diags[iplot][ichanl,10,:]<999)
        loc_use = np.where(diags[iplot][ichanl,10,:]==0)
#        loc_cloudy = np.where(diags[iplot][ichanl,10,:]==50.0)
        loc_cloudy = np.where((diags[iplot][ichanl,10,:]==50.0) & (diags[iplot][ichanl,4,:]==1.0))
#        if plotnum == 1:
##           loc_use = np.where(diags[iplot][ichanl,12,:]<999)
##           plotvar = diags[iplot][ichanl,12,loc_use]
##           plottitle = 'IWP retrievals'
#           loc_use = np.where(diags[iplot][ichanl,12,:]<999)
#           plotvar = diags[iplot][ichanl,14,:]
#           plottitle = 'FV3GFS First-Guess: IWP'
#        else:
##           loc_use = np.where(diags[iplot][ichanl,13,:]<999)
##           plotvar = diags[iplot][ichanl,13,loc_use]
##           plottitle = 'LWP retrievals'
#           loc_use = np.where(diags[iplot][ichanl,13,:]<999)
#           plotvar = diags[iplot][ichanl,15,:]
#           plottitle = 'FV3GFS First-Guess: LWP'
        if (np.mod(plotnum,2) == 0):
           plotvar = diags[iplot][ichanl,6,loc_use]
           xx, yy = m(diags[iplot][ichanl,2,loc_use],diags[iplot][ichanl,1,loc_use])
#           plotvar = diags[iplot][ichanl,6,loc_cloudy]
#           xx, yy = m(diags[iplot][ichanl,2,loc_cloudy],diags[iplot][ichanl,1,loc_cloudy])
        else:
           plotvar = diags[iplot][ichanl,6,loc_all]
           xx, yy = m(diags[iplot][ichanl,2,loc_all],diags[iplot][ichanl,1,loc_all])
        varmean = np.mean(plotvar)
#        plotvar = diags[iplot][ichanl,3,loc_use]
#        xx, yy = m(diags[iplot][ichanl,2,:],diags[iplot][ichanl,1,:])
        nn = xx.shape[1]
#        nn = xx.shape[0]
        plottitle = '{} Obs'.format(titles[iplot])
        plotvmin, plotvmax = tbmin[ichanl], tbmax[ichanl] 
#        plotvmin, plotvmax = 0, 1.0 
#        plotvmin, plotvmax = -3.0, 3.0 
        cs = m.scatter(xx,yy,c=plotvar,s=20,edgecolor='none',vmin=plotvmin, vmax=plotvmax,cmap=plt.cm.jet)
#        cs = m.scatter(xx,yy,c=plotvar,s=20,edgecolor='none',vmin=plotvmin, vmax=plotvmax,cmap=plt.cm.rainbow)

        divider = make_axes_locatable(ax)
#        plt.annotate('Data from: {}'.format(fnames[iplot]), (0,0), (0,-0.2), xycoords='axes fraction')
#        plt.annotate('Data from: {}'.format(fnames[iplot]), (0,0), (0,0.05), xycoords='axes fraction')

        cbar = m.colorbar(cs)
        cbar.set_label('Tb [K]',fontsize=12)
#        cbar.set_label('[kg m-2]',fontsize=16)
        cbar.ax.tick_params(labelsize=12)

        ax.set_title('{}: ch {} # = {} (mean = {})'.format(plottitle,ichanl+1,nn,round(varmean,2)),fontsize=20,fontweight='bold')
        fig.subplots_adjust(wspace=.15)
        fig.subplots_adjust(hspace=.15)
#    fig.savefig('tempestd_iwp_lwp_retrievals.png')
#    fig.savefig('tempestd_iwp_lwp_firstguess.png')
#    fig.savefig('tb_ch{}_{}.png'.format(ichanl+1,datetime))
#    fig.savefig('obs_count_tb_ch{}_{}.png'.format(ichanl+1,datetime))



if __name__ == '__main__':
    plt.ion()
#    plot_omb_vs_oma()
#    plot_omb_hist(0)
#    plot_omb_hist(1)
#    plot_omb_hist(2)
#    plot_omb_hist(3)
#    plot_omb_hist(4)
#    plot_oma_hist(0)
#    plot_oma_hist(1)
#    plot_oma_hist(2)
#    plot_oma_hist(3)
#    plot_oma_hist(4)
#    plot_omb_vs_oma_hist(0)
#    plot_omb_vs_oma_hist(1)
#    plot_omb_vs_oma_hist(2)
#    plot_omb_vs_oma_hist(3)
#    plot_omb_vs_oma_hist(4)
    plot_orbits(0)
#    plot_orbits(1)
#    plot_orbits(2)
#    plot_orbits(3)
#    plot_orbits(4)
    plt.show()
