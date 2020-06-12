#!/usr/bin/python

import numpy as np
import netCDF4
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap


nx_fv3 = 1536
ny_fv3 = 768

path_fv3 = '/mnt/ssdatenas/ting-chi/fv3gfs_experiments/control_dec2018/gfs.forecast_derived/gfs.2018120812'
#path_fv3 = '/mnt/ssdatenas/ting-chi/fv3gfs_experiments/control_may2019/gfs.forecast_derived/gfs.2019051212'
path_clima = '/mnt/ssdatenas/ting-chi/ncep_ncar_r1'

latfname = '{}/xlat.dat'.format(path_fv3)
lonfname = '{}/xlon.dat'.format(path_fv3)
f1fname = '{}/500HGT.f000.2018120812'.format(path_fv3)
f2fname = '{}/TPW.f000.2018120812'.format(path_fv3)
#f1fname = '{}/500HGT.f000.2019051212'.format(path_fv3)
#f2fname = '{}/TPW.f000.2019051212'.format(path_fv3)
phifname = '{}/data_phi.nc'.format(path_clima)
pwatfname = '{}/data_pwat.nc'.format(path_clima)

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

def writebin(filename,var):

    ndims = len(var.shape)
    if ( ndims == 2 ):
      n1, n2 = var.shape[0], var.shape[1]
      print 'ndims of var = ', n1, n2, ' file size = ', 4*n1*n2
      dat = np.zeros([n1,n2],dtype='f4')
      dat[:,:] = var[:,:]
    elif ( ndims == 3 ):
      n1, n2, n3 = var.shape[0], var.shape[1], var.shape[2]
      print 'ndims of var = ', n1, n2, n3, ' file size = ', 4*n1*n2*n3
      dat = np.zeros([n1,n2,n3],dtype='f4')
      dat[:,:,:] = var[:,:,:]

# Python flips row automatically :
# need to make a plot with GrADs to make sure
# the y-axis is not up side down
#
#    dat=np.flipud(dat)
    outfile = open(filename,'wb')
    dat.tofile(outfile)
    outfile.close()

    return()


#----------------------------------------------------#
#   Manipulate Data
#----------------------------------------------------#

fv3lat = readbin(latfname, 1, nx_fv3, ny_fv3)
fv3lon = readbin(lonfname, 1, nx_fv3, ny_fv3)
fv3p500hgt = readbin(f1fname, 1, nx_fv3, ny_fv3)
fv3tpw = readbin(f2fname, 1, nx_fv3, ny_fv3)

fid = netCDF4.Dataset(phifname, 'r')
phi_clima = fid.variables['phiclim'][...]
fid = netCDF4.Dataset(pwatfname, 'r')
pwat_clima = fid.variables['PWAT'][...]

#phi_clima.shape
#Out[6]: (12, 3, 73, 144)
# 3: 925, 500, 250 mb

p500hgt_dec2018 = phi_clima[-1,1,:,:]
pwat_dec2018 = pwat_clima[707,:,:]

p500hgt_may2019 = phi_clima[4,1,:,:]
pwat_may2019 = pwat_clima[712,:,:]

lon = fid.variables['X'][...]
lat = fid.variables['Y'][...]
nx_clima = len(lon)
ny_clima = len(lat)
lons, lats = np.meshgrid(lon,lat)

#loncoord_clima = np.zeros(nx_clima)
#
#for iin in range(nx_clima):
#    for iout in range(nx_fv3):
#        if ( iout+1 >= nx_fv3 ):
#           n3 = fv3lon[0,0] + 360.0
#        else:
#           n3 = fv3lon[0,iout+1]
##        print('{}: n3, fv3lon[0,iout] = {}, {}'.format(iout,n3, fv3lon[0,iout]))
#        if ( lon[iin] < n3 and lon[iin] >= fv3lon[0,iout] ):
#           n1 = np.float(iout) 
#           n2 = (lon[iin] - fv3lon[0,iout])/(n3 - fv3lon[0,iout])
#           print('{}: n1, n2 = {}, {}'.format(iin,n1,n2))
#           loncoord_clima[iin] = n1 + n2
#
#latcoord_clima = np.zeros(ny_clima)
#for iin in range(1,ny_clima-1):
#    for iout in range(ny_fv3-1):
#        n3 = fv3lat[iout+1,0]
#        if ( lat[iin] < n3 and lat[iin] >= fv3lat[iout,0] ):
#           n1 = np.float(iout) 
#           n2 = (lat[iin] - fv3lat[iout,0])/(n3 - fv3lat[iout,0])
##           print('{}: n1, n2 = {}, {}'.format(iin,n1,n2))
#           latcoord_clima[iin] = n1 + n2
#latcoord_clima[0] = ny_fv3 - 1
#latcoord_clima[-1] = 0.

#----------------------------------------------------#
#  coordinates for mapping climate to fv3 grids
#----------------------------------------------------#

lon_coord_fv3inclima = np.zeros(nx_fv3)
for iin in range(nx_fv3):
    for iout in range(nx_clima):
        if (iout+1 >= nx_clima):
           n3 = lon[0] + 360.0
        else:
           n3 = lon[iout+1]
#        print('{}: n3, fv3lon[0,iin] = {}, {}'.format(iout,n3, fv3lon[0,iin]))
        if (fv3lon[0,iin] < n3 and fv3lon[0,iin] >= lon[iout]):
           n1 = np.float(iout)
           n2 = (fv3lon[0,iin] - lon[iout])/(n3 - lon[iout]) 
           lon_coord_fv3inclima[iin] = n1 + n2

lat_coord_fv3inclima = np.zeros(ny_fv3)
for iin in range(ny_fv3):
    for iout in range(ny_clima-1):
        n3 = lat[iout+1]
 #       print('{}: n3, fv3lat[iin,0] = {}, {}'.format(iin,n3, fv3lat[iin,0]))
        if (fv3lat[iin,0] >= n3 and fv3lat[iin,0] < lat[iout]):
           n1 = np.float(iout)
           n2 = (fv3lat[iin,0] - lat[iout])/(n3 - lat[iout]) 
#           print('{}: n1, n2 = {}, {}'.format(iin,n1,n2))
           lat_coord_fv3inclima[iin] = n1 + n2

#----------------------------------------------------#
#  mapping climate to fv3 grids
#----------------------------------------------------#

varin = p500hgt_dec2018
#varin = pwat_dec2018
#varin = p500hgt_may2019
#varin = pwat_may2019
varout = np.zeros([ny_fv3,nx_fv3])

#for iyyr in range(ny_fv3):
#    iyy = ny_fv3 - iyyr - 1
#    idy = np.int(lat_coord_fv3inclima[iyy])
for iyy in range(ny_fv3):
    idy = np.int(lat_coord_fv3inclima[iyy])
#    idy = ny_clima - np.int(lat_coord_fv3inclima[iyy]) - 1
    if ( idy+1 >= ny_clima ):
       idyp = ny_clima - 1
    else:
       idyp = idy + 1
    dely = lat_coord_fv3inclima[iyy]-idy
    for ixx in range(nx_fv3):
        idx = np.int(lon_coord_fv3inclima[ixx])
        delx = lon_coord_fv3inclima[ixx]-idx
#       (idx,idy) is the upperleft corner or the 4 points unit        
        if ( idx+1 >= nx_clima ):
           idxp = 0
        else:
           idxp = idx + 1

        varout[iyy,ixx] = (1.0-delx)*(1.0-dely)*varin[idy,idx] + (1.0-delx)*dely*varin[idyp,idx] + delx*dely*varin[idyp,idxp] + delx*(1.0-dely)*varin[idy,idxp]

writebin('500HGT_December_climatology.dat',varout)
varout_readin = readbin('500HGT_December_climatology.dat', 1, nx_fv3, ny_fv3)
#writebin('TPW_December_climatology.dat',varout)
#varout_readin = readbin('TPW_December_climatology.dat', 1, nx_fv3, ny_fv3)

#writebin('500HGT_May_climatology.dat',varout)
#varout_readin = readbin('500HGT_May_climatology.dat', 1, nx_fv3, ny_fv3)
#writebin('TPW_May_climatology.dat',varout)
#varout_readin = readbin('TPW_May_climatology.dat', 1, nx_fv3, ny_fv3)

plt.ion()
#plt.figure();plt.pcolor(varin);plt.colorbar();plt.show()
#plt.figure();plt.pcolor(varout);plt.colorbar();plt.show()
#plt.figure();plt.pcolor(np.flip(varout,axis=0));plt.colorbar();plt.show()

fig = plt.figure(num=None, figsize=(24.0, 16.0), dpi=80, facecolor='w', edgecolor='k')

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

    if plotnum == 1:
       xx, yy = m(lons,lats)
       plotvar = varin
       plotcmap = plt.cm.terrain
       plottitle = 'NCEP-NCAR R1: Monthly 500 mb GeopHeight'
#       plottitle = 'NCEP-NCAR R1: Monthly TPW'

    elif plotnum == 2:
       xx, yy = m(fv3lon,fv3lat)
       plotvar = varout
#       plotcmap = plt.cm.seismic
       plotcmap = plt.cm.terrain
       plottitle = 'NCEP-NCAR R1: Monthly 500 mb GeopHeight on C384 grid'
#       plottitle = 'NCEP-NCAR R1: Monthly TPW on C384 grid'

    elif plotnum == 3:
       xx, yy = m(fv3lon,fv3lat)
       plotvar = fv3p500hgt
#       plotvar = fv3tpw
#       plotcmap = plt.cm.seismic
       plotcmap = plt.cm.terrain
       plottitle = 'FV3GFS: 500 mb GeopHeight on C384 grid'
#       plottitle = 'FV3GFS: TPW on C384 grid'

    elif plotnum == 4:
       xx, yy = m(fv3lon,fv3lat)
       plotvar = varout_readin
#       plotcmap = plt.cm.seismic
       plotcmap = plt.cm.terrain
       plottitle = 'NCEP-NCAR R1 (read dat): 500 mb GeopHeight on C384 grid'
#       plottitle = 'NCEP-NCAR R1 (read dat): TPW on C384 grid'

    plotmax = 5800.0
    plotmin = 5100.0
#    plotmax = 5.0
#    plotmin = 65.0

    m.pcolor(xx,yy,plotvar,vmin=plotmin,vmax=plotmax,cmap=plotcmap)
    plt.colorbar()

#    plt.annotate('Data from: {}'.format(ffnames[ifile]), (0,0), (0,0.05), xycoords='axes fraction')
    ax.set_title(plottitle,fontsize=20,fontweight='bold')
    fig.subplots_adjust(wspace=.25,hspace = .3)

#m.pcolor(xx, yy, p500hgt_dec2018)
#m.contour(xx, yy, p500hgt_dec2018)
#m.pcolor(xx, yy, pwat_dec2018)
#m.pcolor(xx, yy, varin)
#m.contour(xx, yy, pwat_dec2018)
#plt.colorbar()

plt.savefig('p500hgt_December_clima_on_c384.png')
#plt.savefig('p500hgt_May_clima_on_c384.png')
#plt.savefig('tpw_December_clima_on_c384.png')
#plt.savefig('tpw_May_clima_on_c384.png')
plt.show()
