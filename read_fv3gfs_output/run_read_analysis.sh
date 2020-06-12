#!/bin/sh

# Author: Ting-Chi Wu
# Date: 2020/3/25
# Purpose: run read_fv3gfs_atm_nemsio_anl on analysis from cycled FV3GFS experiments 
#          to extract 3D fields of u-/v-winds, temperatures, pressure, specific humidity,
#          and a 2D field of surface pressure
#          then make GrADs control files to visualize them (with time loop capability)
#       
# Reference: adpated from run_read_production.sh and run_read_forecast.sh

#------------------------------------
# User Defined Variables
#------------------------------------
#export stime='2018120812'
export stime='2019051212'
#export expname='control_dec2018'
export expname='control_may2019'
#export expname='mhs_dec2018'
#export expname='tempestd_dec2018'
##export expname='tempestd_dec2018_calidiff'
export utildir=/scratch1/NESDIS/nesdis-rdo2/Ting-Chi.Wu/utils
#export comrot=/scratch1/NESDIS/nesdis-rdo2/Ting-Chi.Wu/comrot
export comrot=/scratch1/NESDIS/nesdis-rdo2/Ting-Chi.Wu/output_postproc
export cdump='gfs'
#export cdump='gdas'
export cycleint=12 # hour
#export ncycles=9 #(gfs analysis in 12-hour interval)
export ncycles=21 #(gfs analysis in 12-hour interval)
#export ncycles=1 #(gfs analysis in 12-hour interval)

cd $comrot/${expname}/

#------------------------------------
# Create folders to store data
#------------------------------------

if [ ! -d "${cdump}.analysis" ]; then
   mkdir "${cdump}.analysis"
fi

cd ${cdump}.analysis/
export curdir=`pwd`

#------------------------------------
# Manipulate with Dates
#------------------------------------
export sYYYY=`echo $stime | cut -c1-4`
export sMM=`echo $stime | cut -c5-6`
export sDD=`echo $stime | cut -c7-8`
export sHH=`echo $stime | cut -c9-10`
export sMonth=`date --date="${sYYYY}${sMM}${sDD} +${sHH} hours" +%b`
echo "month = " ${sMonth}

#------------------------------------
# Loop over forecast hours (0-240)
#------------------------------------

count=1

while [ $count -le $ncycles ]
do
#  cychour=`expr $count \* $cycleint`
  cychour=`echo "($count -1)*$cycleint" | bc`
  runtime=`date --date="${sYYYY}${sMM}${sDD} +${sHH} hours +${cychour}hour" +%Y%m%d%H`
  rMonth=`date --date="${sYYYY}${sMM}${sDD} +${sHH} hours +${cychour}hour" +%b`
  echo $cychour $runtime
  export rYYYY=`echo $runtime | cut -c1-4`
  export rMM=`echo $runtime | cut -c5-6`
  export rDD=`echo $runtime | cut -c7-8`
  export rHH=`echo $runtime | cut -c9-10`

  if [ ! -d "${cdump}.${runtime}" ]; then
     mkdir "${cdump}.${runtime}"
  fi
  cd ${cdump}.${runtime}/
 
  echo "currently at ${curdir}/${cdump}.${runtime}"


  rm *000000
  rm *.ctl
  rm *.dat
  rm out.read_anl

  echo "link ${comrot}/${expname}/${cdump}.${rYYYY}${rMM}${rDD}/${rHH}/${cdump}.t${rHH}z.atmanl.nemsio"
  ln -s ${comrot}/${expname}/${cdump}.${rYYYY}${rMM}${rDD}/${rHH}/${cdump}.t${rHH}z.atmanl.nemsio .

  echo "${utildir}/read_nemsio/read_fv3gfs_atm_nemsio_anl ${cdump}.t${rHH}z.atmanl.nemsio > out.read_anl 2>&1"
  ${utildir}/read_nemsio/read_fv3gfs_atm_nemsio_anl ${cdump}.t${rHH}z.atmanl.nemsio > out.read_anl 2>&1 

#------------------------------------
# Generate GrADs control files
#------------------------------------

  cat > press.ctl << EOF
DSET PRESS000000
TITLE FV3GFS run
UNDEF  -9999999.0
xdef 1536 linear   0.00000   0.234375
ydef  768 linear -89.91032   0.234375
ZDEF 64 LINEAR 1 1
TDEF 1 LINEAR ${rHH}:${rHH}Z${rDD}${rMonth}${rYYYY} 1hr
VARS 1
press 64 99 pressure (hPa)
ENDVARS
EOF

  cat > tempk.ctl << EOF
DSET TEMPK000000
TITLE FV3GFS run
UNDEF  -9999999.0
xdef 1536 linear   0.00000   0.234375
ydef  768 linear -89.91032   0.234375
ZDEF 64 LINEAR 1 1
TDEF 1 LINEAR ${rHH}:${rHH}Z${rDD}${rMonth}${rYYYY} 1hr
VARS 1
tempk 64 99 temperature (K)
ENDVARS
EOF

  cat > sph.ctl << EOF
DSET SPH000000
TITLE FV3GFS run
UNDEF  -9999999.0
xdef 1536 linear   0.00000   0.234375
ydef  768 linear -89.91032   0.234375
ZDEF 64 LINEAR 1 1
TDEF 1 LINEAR ${rHH}:${rHH}Z${rDD}${rMonth}${rYYYY} 1hr
VARS 1
sph 64 99 specific humidity (kg/kg)
ENDVARS
EOF

  cat > uwind.ctl << EOF
DSET U000000
TITLE FV3GFS run
UNDEF  -9999999.0
xdef 1536 linear   0.00000   0.234375
ydef  768 linear -89.91032   0.234375
ZDEF 64 LINEAR 1 1
TDEF 1 LINEAR ${rHH}:${rHH}Z${rDD}${rMonth}${rYYYY} 1hr
VARS 1
uwind 64 99 u-component wind (m/s)
ENDVARS
EOF

  cat > vwind.ctl << EOF
DSET V000000
TITLE FV3GFS run
UNDEF  -9999999.0
xdef 1536 linear   0.00000   0.234375
ydef  768 linear -89.91032   0.234375
ZDEF 64 LINEAR 1 1
TDEF 1 LINEAR ${rHH}:${rHH}Z${rDD}${rMonth}${rYYYY} 1hr
VARS 1
vwind 64 99 v-component wind (m/s)
ENDVARS
EOF

  cat > psfc.ctl << EOF
DSET PSFC000000
TITLE FV3GFS run
UNDEF  -9999999.0
xdef 1536 linear   0.00000   0.234375
ydef 768 linear -89.91032   0.234375
ZDEF 1 LEVELS 0
TDEF 1 LINEAR ${rHH}:${rHH}Z${rDD}${rMonth}${rYYYY} 1hr
VARS 1
psfc    1 99 surface pressure (hPa)
ENDVARS
EOF

  count=`expr $count + 1`
  cd $curdir/
done


exit
