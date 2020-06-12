#!/bin/sh

##---- Hera JOBCARD
##---- Submit as: sbatch $script
##SBATCH -J read_production
##SBATCH -A nesdis-rdo1
##SBATCH --partition=hera
##SBATCH -o log.read_prod
##SBATCH -e log.read_prod
##SBATCH --ntasks=1
##SBATCH -q batch
##SBATCH -t 8:00:00
#
# Author: Ting-Chi Wu
# Date: 2020/2/17
# Purpose: run read_fv3gfs_atm_nemsio_derived on GFS production (operational analysis) 
#          to extract 500 mb geopotential height and total precipitable water field
#       
# Revisions:
# - T.-C. Wu 2020/3/25 add call to read_fv3gfs_atm_nemsio_anl to extract
#                      3D fields of u-/v-winds, temperatures, pressure, specific humidity,
#                      and surface pressure
# Reference: TBD

#------------------------------------
# User Defined Variables
##------------------------------------
#export stime='2018120812'
#export stime='2018120906'
#export stime='2018122118'
#export stime='2018122400'
export stime='2019051212'
export utildir=/scratch1/NESDIS/nesdis-rdo2/Ting-Chi.Wu/utils
#export cdump='gfs'
export cdump='gdas'
export cycleint=6
#export ncycles=60
#export ncycles=64
export ncycles=84
#export ncycles=0
#export ncycles=1

export curdir=`pwd`

#------------------------------------
# Manipulate with Dates
#------------------------------------
export sYYYY=`echo $stime | cut -c1-4`
export sMM=`echo $stime | cut -c5-6`
export sDD=`echo $stime | cut -c7-8`
export sHH=`echo $stime | cut -c9-10`

count=0

while [ $count -le $ncycles ]
do
  cychour=`expr $count \* $cycleint`
  runtime=`date --date="${sYYYY}${sMM}${sDD} +${sHH} hours +${cychour}hour" +%Y%m%d%H`
  rMonth=`date --date="${sYYYY}${sMM}${sDD} +${sHH} hours +${cychour}hour" +%b`
  echo $cychour $runtime
  export rYYYY=`echo $runtime | cut -c1-4`
  export rMM=`echo $runtime | cut -c5-6`
  export rDD=`echo $runtime | cut -c7-8`
  export rHH=`echo $runtime | cut -c9-10`

  echo "month = " ${rMonth}

  cd ${utildir}/${cdump}_production/${cdump}.${rYYYY}${rMM}${rDD}/${rHH}/
  echo "we are at ${utildir}/${cdump}_production/${cdump}.${rYYYY}${rMM}${rDD}/${rHH}"

  rm *000000
  rm *.ctl
  rm *.dat
  rm out.read_derived
  rm out.read_anl

  echo "${utildir}/read_nemsio/read_fv3gfs_atm_nemsio_derived ${cdump}.t${rHH}z.atmanl.nemsio > out.read_derived 2>&1"
  ${utildir}/read_nemsio/read_fv3gfs_atm_nemsio_derived ${cdump}.t${rHH}z.atmanl.nemsio > out.read_derived 2>&1 

  cat > hsl.ctl << EOF
DSET 500HGT000000
TITLE FV3GFS run
UNDEF  -9999999.0
xdef 1536 linear   0.00000   0.234375
ydef 768 linear -89.91032   0.234375
ZDEF 1 LEVELS 0
TDEF 1 LINEAR ${rHH}:${rHH}Z${rDD}${rMonth}${rYYYY} 1hr
VARS 1
hsl     1 99 500 mb height (m)
ENDVARS
EOF

  cat > tpw.ctl << EOF
DSET TPW000000
TITLE FV3GFS run
UNDEF  -9999999.0
xdef 1536 linear   0.00000   0.234375
ydef 768 linear -89.91032   0.234375
ZDEF 1 LEVELS 0
TDEF 1 LINEAR ${rHH}:${rHH}Z${rDD}${rMonth}${rYYYY} 1hr
VARS 1
tpw     1 99 total precipitable water (mm)
ENDVARS
EOF

  if [ ${rHH} == "00" ] || [ ${rHH} == "12" ]; then

    echo "${utildir}/read_nemsio/read_fv3gfs_atm_nemsio_anl ${cdump}.t${rHH}z.atmanl.nemsio > out.read_anl 2>&1"
    ${utildir}/read_nemsio/read_fv3gfs_atm_nemsio_anl ${cdump}.t${rHH}z.atmanl.nemsio > out.read_anl 2>&1 

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

  fi

  count=`expr $count + 1`
  cd $curdir/
done

exit
