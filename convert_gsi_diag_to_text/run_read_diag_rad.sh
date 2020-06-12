#!/bin/sh

# Author: Ting-Chi Wu
# Date: 2019/12/16
# Purpose: untar gdas.tHHz.radstat to extract diag files
#          and use the gsi.fd/util/Analysis_Utilities/read_diag tool 
#          to convert diag (binary) to ascii (text) file for post-processing
#       
# Reference: check the flow diagram of a typical FV3GFS experiment

#------------------------------------
# User Defined Variables
#------------------------------------
#export cdate8='20181212'
#export chour2='12'
export cdate8='20190522'
export chour2='12'
#export expname='control_dec2018'
export expname='control_may2019'
#export expname='mhs_dec2018'
#export expname='tempestd_dec2018'
#export expname='tempestd_dec2018_calidiff'
#export comrot=/scratch1/NESDIS/nesdis-rdo2/Ting-Chi.Wu/comrot
export comrot=/scratch1/NESDIS/nesdis-rdo2/Ting-Chi.Wu/output_postproc
export gsidir=/scratch1/NESDIS/nesdis-rdo2/Ting-Chi.Wu/fv3gfs/global-workflow/sorc/gsi.fd
export cdump='gfs'
#export cdump='gdas'

echo $comrot
echo $gsidir

cd $comrot/${expname}/

if [ ! -d "${cdump}.diag_rad" ]; then
   mkdir "${cdump}.diag_rad"
fi

#----------------------------------------------
# Untar radstat and extract diag files
#----------------------------------------------
cd ${comrot}/${expname}/${cdump}.${cdate8}/${chour2}
echo `pwd`
if [ -f "$cdump.t${chour2}z.radstat" ]; then
   mkdir diag_rad
   cd diag_rad/
   cp ../$cdump.t${chour2}z.radstat .
   tar -xf $cdump.t${chour2}z.radstat diag_mhs_metop-a_ges.${cdate8}${chour2}.gz diag_mhs_metop-a_anl.${cdate8}${chour2}.gz 
   tar -xf $cdump.t${chour2}z.radstat diag_mhs_metop-b_ges.${cdate8}${chour2}.gz diag_mhs_metop-b_anl.${cdate8}${chour2}.gz 
   tar -xf $cdump.t${chour2}z.radstat diag_mhs_n19_ges.${cdate8}${chour2}.gz diag_mhs_n19_anl.${cdate8}${chour2}.gz 
   tar -xf $cdump.t${chour2}z.radstat diag_tempest_cubesat_ges.${cdate8}${chour2}.gz diag_tempest_cubesat_anl.${cdate8}${chour2}.gz 
   gunzip diag_*.gz
   ls -1
fi

#------------------------------------
# Manipulate with Dates
#------------------------------------

cd ${gsidir}/util/Analysis_Utilities/read_diag
sensors=('mhs_metop-a' 'mhs_metop-b' 'mhs_n19' 'tempest_cubesat')
types=('ges' 'anl')

#echo ${sensors[0]} ${sensors[1]}

nsensor=4
ntype=2

dd=0
while [ $dd -lt $ntype ]
do
  ss=0
  while [ $ss -lt $nsensor ]
  do 
    rm namelist.rad
    cat << EOF > namelist.rad
    &iosetup
    infilename='./diag_${sensors[$ss]}_${types[$dd]}',
    outfilename='./results_${sensors[$ss]}_${types[$dd]}',
    /
EOF
    ln -s ${comrot}/${expname}/${cdump}.${cdate8}/${chour2}/diag_rad/diag_${sensors[$ss]}_${types[$dd]}.${cdate8}${chour2} diag_${sensors[$ss]}_${types[$dd]}
    ./read_diag_rad.exe
    rm diag_${sensors[$ss]}_${types[$dd]}
    mv results_${sensors[$ss]}_${types[$dd]} results_${sensors[$ss]}_${types[$dd]}.${cdate8}${chour2}
    ss=$(( ss+1 ))
  done
  dd=$(( dd+1 ))
done
ls -1 results_*
#mv results_* ${comrot}/${expname}/${cdump}.${cdate8}/${chour2}/diag_rad
cp results_* ${comrot}/${expname}/${cdump}.${cdate8}/${chour2}/diag_rad
mv results_* ${comrot}/${expname}/${cdump}.diag_rad
