#!/bin/bash
#
### Define some variables
#
export hour=0
export count=0
#export count=36
export count_max=240
export atm_prefix=gfs.t12z.atmf
export sfc_prefix=gfs.t12z.sfcf
export file_suffix=.nemsio
export begin_yyyymmdd=20181212
export begin_hour=12
#
### Define absolute paths (requires trailing slash!)
#
top_directory='/scratch1/NESDIS/nesdis-rdo2/Heather.Cronk/control_dec2018/'
current_directory=$top_directory'crtm_scripts/'
nemsio_directory=$top_directory'read_nemsio/'
gmi_orbit='/scratch1/NESDIS/nesdis-rdo2/Lewis.Grasso/gmi_orbit/'
crtm=$top_directory'crtm_interface/'
case='control_dec2018/'

experiment=$nemsio_directory$case
nemsio_source='/scratch2/NESDIS/nesdis-rdo1/Milija.Zupanski/output_postproc/'
gfdl_diam_directory=$top_directory'gfdl_microphysics/'
#
### Purge old stuff
#echo 'cd to '$crtm$case
#cd $crtm$case
#echo 'Removing all files in '$crtm$case
#rm *
#
### Let's get the ball rolling, shall we?
#
echo ' '
echo ' '
echo 'Beginning crtm_gmi.sh'
echo 'top directory = '$top_directory
echo 'current_directory = '$current_directory
echo 'nemsio_directory = '$nemsio_directory
echo 'experiment = '$experiment
echo 'nemsio_source = '$nemsio_source
echo ' '


while [ $count -le $count_max ] ; do
  
  echo 'cd to '$current_directory
  cd $current_directory
  echo ' '

  if [ $count -le 9 ] ; then
     export hour=00$count
  fi
  if [ 12 -le $count -a $count -le 99 ] ; then
     export hour=0$count
  fi
  if [ 100 -le $count -a $count -le 240 ] ; then
     export hour=$count
  fi
  echo ' '
  echo '*****************************Processing***************************'
  echo ' '
  atm_file=$atm_prefix$hour$file_suffix
  sfc_file=$sfc_prefix$hour$file_suffix
  obs_date=`date --date="$begin_yyyymmdd + $begin_hour hours + $hour hours" +%Y%m%d%H`
  echo 'Observations date '$obs_date
  echo 'atm filename '$atm_file
  echo 'sfc filename '$sfc_file
  echo 'count = '$count


  if [ 1 -eq 1 ] ; then

  echo ' '
  echo '                     **********   Step 1: Run read_fv3gfs_atm_nemsio.out   **********'
  echo '                     **********               read_fv3gfs_sfc_nemsio.out   **********'
  echo ' '
  echo 'cd to '$experiment
  cd $experiment
  echo 'removing all files in '$experiment
  rm *
  echo 'create symbolic link to model output files'
  echo 'ln -s '$nemsio_source$case'gfs.'$begin_yyyymmdd'/'$begin_hour'/'$atm_file
  ln -s $nemsio_source$case'gfs.'$begin_yyyymmdd'/'$begin_hour'/'$atm_file
  #cp $nemsio_source$case'gfs.'$begin_yyyymmdd'/'$begin_hour'/'$atm_file $experiment$atm_file
  echo 'ln -s '$nemsio_source$case'gfs.'$begin_yyyymmdd'/'$begin_hour'/'$sfc_file
  ln -s $nemsio_source$case'gfs.'$begin_yyyymmdd'/'$begin_hour'/'$sfc_file
  #cp $nemsio_source$case'gfs.'$begin_yyyymmdd'/'$begin_hour'/'$sfc_file $experiment$sfc_file
  echo 'Lets copy nemsio executables from '$nemsio_directory
  echo 'to here '$experiment
  cp $nemsio_directory'read_fv3gfs_atm_nemsio.out' $experiment'read_fv3gfs_atm_nemsio.out'
  cp $nemsio_directory'read_fv3gfs_sfc_nemsio.out' $experiment'read_fv3gfs_sfc_nemsio.out'
  echo 'Lets extract atm schtuff.'
  nohup $experiment'read_fv3gfs_atm_nemsio.out' $atm_file > $experiment'read_fv3gfs_atm_nemsio.list' 2>&1
  #$experiment'read_fv3gfs_atm_nemsio.out' $atm_file > $experiment'read_fv3gfs_atm_nemsio.list' 2>&1
  echo 'Lets extract sfc schtuff.'
  nohup $experiment'read_fv3gfs_sfc_nemsio.out' $sfc_file > $experiment'read_fv3gfs_sfc_nemsio.list' 2>&1
  #$experiment'read_fv3gfs_sfc_nemsio.out' $sfc_file > $experiment'read_fv3gfs_sfc_nemsio.list' 2>&1
  fi
 

  if [ 1 -eq 1 ] ; then

  echo ' '
  echo '                     **********   Step 2: Run gfdl_diams.out  **********'
  echo ' '
  echo '                     You must edit gfdl_diams.f90 to set correct path'
  echo '                               Do not forget to recompile '
  echo ' '
  echo 'cd to '$gfdl_diam_directory
  cd $gfdl_diam_directory
  echo 'Lets compute GFDL diameters.'
  nohup $gfdl_diam_directory'gfdl_diams.out' > $gfdl_diam_directory'gfdl_diams.list' 2>&1
  #$gfdl_diam_directory'gfdl_diams.out' > $gfdl_diam_directory'gfdl_diams.list' 2>&1
  fi

  if [ 1 -eq 1 ] ; then

  echo ' '
  echo '                    **********   Step 3: Run remap_gmi_to_fv3.out   **********'
  echo ' '
#  echo 'cd to '$crtm$case
#  cd $crtm$case
#  echo 'Removing all files in '$crtm$case
#  rm *
#  echo ' '
  echo 'cd to '$crtm
  cd $crtm
  rm $crtm'GMI_lat.dat'
  rm $crtm'GMI_lon.dat'
  echo 'Copying '$gmi_orbit'GMI-CRTM-LAT.'$obs_date'.dat' $crtm'GMI_lat.dat'
  echo 'Copying '$gmi_orbit'GMI-CRTM-LON.'$obs_date'.dat' $crtm'GMI_lon.dat'
  cp $gmi_orbit'GMI-CRTM-LAT.'$obs_date'.dat' $crtm'GMI_lat.dat'
  cp $gmi_orbit'GMI-CRTM-LON.'$obs_date'.dat' $crtm'GMI_lon.dat'
  nohup $crtm'remap_gmi_to_fv3.out' > $crtm'remap_gmi_to_fv3.list' 2>&1
  #$crtm'remap_gmi_to_fv3.out' > $crtm'remap_gmi_to_fv3.list' 2>&1
  #make sure output directory exists
  echo 'Making output directory '$crtm$case$begin_yyyymmdd$begin_hour
  mkdir -p $crtm$case$begin_yyyymmdd$begin_hour
  echo 'Copying '$crtm'gmi_orbit_window.dat' $crtm$case$begin_yyyymmdd$begin_hour'/gmi_orbit_window_'$obs_date'_'$hour
  cp $crtm'gmi_orbit_window.dat' $crtm$case$begin_yyyymmdd$begin_hour'/gmi_orbit_window_'$obs_date'_'$hour
  echo 'Moving '$crtm'gmi_orbit_window.dat' to $experiment
  mv $crtm'gmi_orbit_window.dat' $experiment'gmi_orbit_window.dat'
 
  fi

  if [ 1 -eq 1 ] ; then

  echo ' '
  echo '                    **********   Step 4: Run Get_CRTM   **********'
  echo ' '
  echo '                     You must edit Get_CRTM.f90 to set correct path'
  echo '                               Do not forget to recompile '
  echo ' '
  echo 'cd to '$crtm
  cd $crtm
  echo 'Running Get_CRTM'
  nohup $crtm'Get_CRTM' > $crtm'Get_CRTM.list' 2>&1
  #$crtm'Get_CRTM' > $crtm'Get_CRTM.list' 2>&1
  echo 'Moving '$crtm'tb_gmi_89VGhz' to $crtm$case$begin_yyyymmdd$begin_hour'/tb_gmi_89VGhz_'$obs_date'_'$hour
  mv $crtm'tb_gmi_89VGhz' $crtm$case$begin_yyyymmdd$begin_hour'/tb_gmi_89VGhz_'$obs_date'_'$hour

  fi
  
  if [ $count -ge 240 ] ; then
    echo ' '
    echo ' '
    echo 'Stopping crtm_gmi.sh; count = '$count
    exit
    echo ' '
  fi
  

  #exit

  let "count=$count + 3" ; export count

done



exit
