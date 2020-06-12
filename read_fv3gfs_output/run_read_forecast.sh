#!/bin/sh

# Author: Ting-Chi Wu
# Date: 2020/2/20
# Purpose: run read_fv3gfs_atm_nemsio_derived on forecast from cycled FV3GFS experiments 
#          to extract 500 mb geopotential height and total precipitable water field
#          and make GrADs control files to visualize them (with time loop capability)
#       
# Reference: adpated from run_read_production.sh

#------------------------------------
# User Defined Variables
#------------------------------------
#export stime='2018120812'
#export stime='2019051212'
export stime='2019052000'
#export expname='control_dec2018'
export expname='control_may2019'
#export expname='mhs_dec2018'
#export expname='tempestd_dec2018'
#export expname='tempestd_dec2018_calidiff'
export utildir=/scratch1/NESDIS/nesdis-rdo2/Ting-Chi.Wu/utils
#export comrot=/scratch1/NESDIS/nesdis-rdo2/Ting-Chi.Wu/comrot
export comrot=/scratch1/NESDIS/nesdis-rdo2/Ting-Chi.Wu/output_postproc
export cdump='gfs'
export cycleint=12 # hour
#export ncycles=9 #(gfs analysis in 12-hour interval)
#export ncycles=17 #(gfs analysis in 12-hour interval)
export ncycles=6 #(gfs analysis in 12-hour interval)
#export cdump='gdas'
export fcstint=6 # hour
export ntime=41 #(all 240-h forecast in 6-hourly interval)

cd $comrot/${expname}/

#------------------------------------
# Create folders to store data
#------------------------------------

if [ ! -d "${cdump}.forecast_derived" ]; then
   mkdir "${cdump}.forecast_derived"
fi

cd ${cdump}.forecast_derived/
export topdir=`pwd`

#------------------------------------
# Manipulate with Dates
#------------------------------------
export sYYYY=`echo $stime | cut -c1-4`
export sMM=`echo $stime | cut -c5-6`
export sDD=`echo $stime | cut -c7-8`
export sHH=`echo $stime | cut -c9-10`

#------------------------------------
# Loop over cycles
#------------------------------------

cycle=1

while [ $cycle -le $ncycles ]
do
  cychour=`echo "($cycle -1)*$cycleint" | bc`
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

  export curdir=`pwd`
  echo "currently at ${curdir}"

  rm *000*
  rm *.ctl
  rm *.dat
  rm out.read_derived*

#------------------------------------
# Start building GrADs control files
#------------------------------------

  cat > hsl.ctl << EOF
DSET ^${comrot}/${expname}/${cdump}.forecast_derived/${cdump}.${runtime}/500HGT.f%ch
EOF
  cat > tpw.ctl << EOF
DSET ^${comrot}/${expname}/${cdump}.forecast_derived/${cdump}.${runtime}/TPW.f%ch
EOF

#------------------------------------
# Loop over forecast hours (0-240)
#------------------------------------
  count=1

  while [ $count -le $ntime ]
  do
#    fcsthour=`expr $count - 1 \* $fcstint`
    fcsthour=`echo "($count -1)*$fcstint" | bc`
    fHHH=$(printf "%.3d" "$fcsthour")
    fvalidtime=`date --date="${rYYYY}${rMM}${rDD} +${rHH} hours +${fHHH}hour" +%Y%m%d%H`
    echo "cycle: $runtime forecast hour = $fHHH valid at $fvalidtime"

    rm ${cdump}.t${rHH}z.atmf*.nemsio
    rm xlat.dat
    rm xlon.dat

    echo "link ${comrot}/${expname}/${cdump}.${rYYYY}${rMM}${rDD}/${rHH}/${cdump}.t${rHH}z.atmf${fHHH}.nemsio"
    ln -s ${comrot}/${expname}/${cdump}.${rYYYY}${rMM}${rDD}/${rHH}/${cdump}.t${rHH}z.atmf${fHHH}.nemsio .

    echo "${utildir}/read_nemsio/read_fv3gfs_atm_nemsio_derived ${cdump}.t${rHH}z.atmf${fHHH}.nemsio > out.read_derived.f${fHHH} 2>&1"
    ${utildir}/read_nemsio/read_fv3gfs_atm_nemsio_derived ${cdump}.t${rHH}z.atmf${fHHH}.nemsio > out.read_derived.f${fHHH} 2>&1 

    mv 500HGT000000 500HGT.f${fHHH}.$fvalidtime
    mv TPW000000 TPW.f${fHHH}.$fvalidtime

# making the CHSUB list for time loop (set loopdim t in GrADs)
    cat >> hsl.ctl << EOF
CHSUB $count $count ${fHHH}.${fvalidtime} 
EOF
    cat >> tpw.ctl << EOF
CHSUB $count $count ${fHHH}.${fvalidtime} 
EOF

    count=`expr $count + 1`
  done

#------------------------------------
# Wrap up GrADs control files
#------------------------------------

  cat >> hsl.ctl << EOF
OPTIONS template
TITLE FV3GFS run
UNDEF  -9999999.0
xdef 1536 linear   0.00000   0.234375
ydef 768 linear -89.91032   0.234375
ZDEF 1 LEVELS 0
TDEF ${ntime} LINEAR ${rHH}:${rHH}Z${rDD}${rMonth}${rYYYY} ${fcstint}hr
VARS 1
hsl     1 99 500 mb height (m)
ENDVARS
EOF

  cat >> tpw.ctl << EOF
OPTIONS template
TITLE FV3GFS run
UNDEF  -9999999.0
xdef 1536 linear   0.00000   0.234375
ydef 768 linear -89.91032   0.234375
ZDEF 1 LEVELS 0
TDEF ${ntime} LINEAR ${rHH}:${rHH}Z${rDD}${rMonth}${rYYYY} ${fcstint}hr
VARS 1
tpw     1 99 total precipitable water (mm)
ENDVARS
EOF

  cycle=`expr $cycle + 1`
  cd $topdir/

done
 
exit
