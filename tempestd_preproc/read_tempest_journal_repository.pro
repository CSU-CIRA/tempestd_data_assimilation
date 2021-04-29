pro read_tempest_journal_repository
;
;*** Time periods for GSI.
;
;*** 0300 utc < centered at 0600 utc <= 0900 utc --------> HDF FILE 1
;*** 0900 utc < centered at 1200 utc <= 1500 utc --------> HDF FILE 2
;*** 1500 utc < centered at 1800 utc <= 2100 utc --------> HDF FILE 3
;*** 2100 utc < centered at 0000 utc <= 0300 utc --------> HDF FILE 4
;
;
;*************************************************************
;*                                                           *
;*                        Part I                             *
;*                                                           *
;*            Define absolute paths and filenames            *
;*                                                           *
;*************************************************************
;
;                                           HHMMSS_begin    HHMMSS_end
;                                           ||||||          ||||||
tempest_data_in  = 'TEMPEST_L1_pub_20190521T235923_20190522T104732_v1.42.h5'
r_data_in        = 'tempest_20190522_v1.42.out.nc'
begin_file_date  =          20190522 ; YYYYMMDD from r_data_in filename
;
;tempest_data_in  = Name of file that contains TEMPEST-D data
;r_data_in        = Name of file that contains retrieved LWP, IWP, and CHI.
;begin_file_date  =          20190522 ; YYYYMMDD from r_data_in filename
;
;tempest_data_out = 'GSI_TEMPEST_L1_pub_20190522T0300_20190522T0900_v1.42.h5'
;tempest_data_out = 'GSI_TEMPEST_L1_pub_20190522T0900_20190522T1500_v1.42.h5'
tempest_data_out = 'GSI_TEMPEST_L1_pub_20190522T1500_20190522T2100_v1.42.h5'
;
;
path_in         = '/mnt/data2/tempest_may2019_in/'
path_out        = '/mnt/data2/tempest_may2019_out/'
;
;path_in         = Path to TEMPETST-D files to be read in
;path_out        = Path to TEMPEST-D files ready for GSI'
;
;
gsi_year        = 2019
gsi_month       = 05
gsi_day         = 22
;
;gsi_begin_hour  = 3.0     ; HHMM     (e.g. 0630 for 6:30am UTC, 2100 for 9pm)
;gsi_end_hour    = 9.0     ; HHMM     (e.g. 0630 for 6:30am UTC, 2100 for 9pm)
;gsi_begin_hour  = 9.0     ; HHMM     (e.g. 0630 for 6:30am UTC, 2100 for 9pm)
;gsi_end_hour    = 15.0    ; HHMM     (e.g. 0630 for 6:30am UTC, 2100 for 9pm)
gsi_begin_hour  = 15.0    ; HHMM     (e.g. 0630 for 6:30am UTC, 2100 for 9pm)
gsi_end_hour    = 21.0    ; HHMM     (e.g. 0630 for 6:30am UTC, 2100 for 9pm)
;
gsi_start_time  = gsi_begin_hour*3600.0 ; in seconds from 0000 UTC
gsi_end_time    = gsi_end_hour*3600.0   ; in seconds from 0000 UTC 
;
;*************************************************************
;*                                                           *
;*                        Part II                            *
;*                                                           *
;*               Read in and get array dimensions            *
;*                                                           *
;*************************************************************
;
print,' '
print,' '
print,' '
print,format='(a,1x,a)', 'Reading',tempest_data_in
print,format='(a,1x,a)', 'Located in',path_in
print,' '
;
a    = read_hdf5(path_in+tempest_data_in)
help,a,/struct	; This will list all variables in the hdf file.
time  = a.scan_utctime  ; Seconds since the beginning of year 2000
tb182_in      = reform(a.scan_tb[*,*,0])  ; 182 GHz Tbs (nominal) channel
tb180_in      = reform(a.scan_tb[*,*,1])  ; 180 GHz Tbs
tb176_in      = reform(a.scan_tb[*,*,2])  ; 176 GHz Tbs
tb165_in      = reform(a.scan_tb[*,*,3])  ; 165 GHz Tbs
tb89_in       = reform(a.scan_tb[*,*,4])  ;  89 GHz Tbs
zenith_in     = a.scan_binc               ; Earth Incidence Angle (degrees) AKA: Zenith
pxl_lat_in    = a.scan_blat               ; Latitude  of radiometer footprint
pxl_lon_in    = a.scan_blon               ; Longitude of radiometer footprint
mask_in       = a.scan_landmask           ; Landmask   [nscan,npix]
scan_angle_in = a.scan_scanang            ; Scan angle (degrees)
;
;*** Get dimensions
;
s     = size(a.scan_utctime)
nscan = long(s[1]) ; Number of scans
npix  = long(s[2]) ; Number of pixels per scan
;
print,' '
print,'          Number of scans:',nscan
print,'Number of pixels per scan:',npix
;
;*** Done reading in time from tempest-d
;
print,' '
print,format='(a,1x,a)','Done reading',tempest_data_in
print,format='(a,1x,a)','Located in',path_in
print,' '
;
;*** Read in lwp, iwp, chi
;
print,' '
print,format='(a)','Now get lwp, iwp, and chi data'
print,format='(a,1x,a)','Reading',r_data_in
print,format='(a,1x,a)','Located in',path_in
print,' '
;
;*** Open netcdf file
;
id = ncdf_open(path_in+r_data_in)
;
;*** Read in cloud liquid water path to variable lwp (units g/m^2)
;
tpid = ncdf_varid(id,'atmosphere_mass_content_of_cloud_liquid_water')
ncdf_varget,id,tpid,lwp
print,'Extraced lwp (g/m^2)'
;
;*** Read in cloud ice water path to variable iwp (units g/m^2)
;
tpid = ncdf_varid(id,'iwp')
ncdf_varget,id,tpid,iwp
print,'Extraced iwp (g/m^2)'
;
;*** Read in chi squared value to variable chi (unitless)
;
tpid = ncdf_varid(id,'chi_squared')
ncdf_varget,id,tpid,chi
print,'Extraced chi (unitless)'
;
;*** Read in observed Tbs (K)
;
tpid = ncdf_varid(id,'tbobs')
ncdf_varget,id,tpid,tbobs
print,'Extraced observed Tb (K)'
;
;*** Read in synthetic Tbs (K)
;
tpid = ncdf_varid(id,'tbsim')
ncdf_varget,id,tpid,tbsim
print,'Extraced synthetic Tb (K)'
;
;*** Close netcdf file
;
ncdf_close,id
;
print,format='(a)',' '
print,format='(a,1x,a)','Done reading',r_data_in
print,format='(a,1x,a)','Located in',path_in
print,format='(a)',' '
;
;*** Create 2-D array 'noconverge' - value of 1 means the retrieval did not converge 
;*** due to either land or irreconcilable TB mismatch (often deep convection).
;*** This might be helpful information for a data assimilation person
;
noconverge = uint(lwp) * 0 ; "uint" is unsigned integer: Type Code = 12
ind = where(tbobs(*,*,0) GT 0.0 AND tbsim(*,*,0) LT 0.0)
noconverge(ind) = 1
;
;*************************************************************
;*                                                           *
;*                        Part III                           *
;*                                                           *
;*              Declare and initialize output arrays         *
;*                                                           *
;*************************************************************
;
;*** Declare arrays
;
     year_array = make_array(nscan,npix,/integer)
    month_array = make_array(nscan,npix,/integer)
      day_array = make_array(nscan,npix,/integer)
     hour_array = make_array(nscan,npix,/integer)
   minute_array = make_array(nscan,npix,/integer)
   second_array = make_array(nscan,npix,/integer)
   converge_out = make_array(nscan,npix,/integer)
;
    pxl_lat_out = make_array(nscan,npix,/float)
    pxl_lon_out = make_array(nscan,npix,/float)
 scan_angle_out = make_array(nscan,npix,/float)
     zenith_out = make_array(nscan,npix,/float)
       tb89_out = make_array(nscan,npix,/float)
      tb165_out = make_array(nscan,npix,/float)
      tb176_out = make_array(nscan,npix,/float)
      tb180_out = make_array(nscan,npix,/float)
      tb182_out = make_array(nscan,npix,/float)
       mask_out = make_array(nscan,npix,/float)
        lwp_out = make_array(nscan,npix,/float)
        iwp_out = make_array(nscan,npix,/float)
        chi_out = make_array(nscan,npix,/float)
 ;
 ;*** Initialize integer arrays with -999
 ;*** Initialize float arrays with -999.0
 ;
 int_missing   = -999
 float_missing = -999.0
 ;
     year_array[*,*] = int_missing
    month_array[*,*] = int_missing
      day_array[*,*] = int_missing
     hour_array[*,*] = int_missing
   minute_array[*,*] = int_missing
   second_array[*,*] = int_missing
   converge_out[*,*] = int_missing
    ;
    pxl_lat_out[*,*] = float_missing
    pxl_lon_out[*,*] = float_missing
 scan_angle_out[*,*] = float_missing
     zenith_out[*,*] = float_missing
       tb89_out[*,*] = float_missing
      tb165_out[*,*] = float_missing
      tb176_out[*,*] = float_missing
      tb180_out[*,*] = float_missing
      tb182_out[*,*] = float_missing
       mask_out[*,*] = float_missing
        lwp_out[*,*] = float_missing
        iwp_out[*,*] = float_missing
        chi_out[*,*] = float_missing
;
;*************************************************************
;*                                                           *
;*                        Part IV                            *
;*                                                           *
;*        Test for GSI time bounds and fill up arrays        *
;*                                                           *
;*************************************************************
;
;*** Compute the number of send lines from start_time and end_time.
;
base_year  = 0     ; Year  of current day
base_month = 0     ; Month of current day
base_day   = 0     ; Day   of current day
base_t     = 0L    ; Time, in seconds, of 0000 UTC for this date from year 2000.
;
t = convert_date(begin_file_date,base_year=2000) ; extract time
;
base_year  = t.year    ; Year of current day
base_month = t.month   ; Month of current day
base_day   = t.day     ; Day of current day
base_t     = t.itime   ; UTC time, in seconds, at start of current day
;
year  = base_year
month = base_month
day   = base_day
;
reset_day = 1
i_begin   = -1
i_end     = -1
for i = 0, nscan-1 do begin
for j = 0, npix-1 do begin
 if ( finite(time[i,j]) eq 1 ) then begin
  ;
  image_time         = (time[i,j]-base_t)/3600.0       ; hours as real number
  hour               = long(image_time)                ; integer part of hour
  fractional_hour    = image_time-float(hour)          ; decimal part of hour
  minutes            = long( fractional_hour*60.0 )    ; integer part of minutes
  fractional_minutes = fractional_hour*60.0-float(minutes)
  seconds            = long( fractional_minutes*60.0 )
  ;
  if ( hour ge 24 ) then begin ; this only works for one day!
    hour = hour - 24
    if ( reset_day eq 1 ) then begin
     day = day + 1
     reset_day = 0
    endif
  endif
  ;
  if ( (year eq gsi_year) and (month eq gsi_month) and (day eq gsi_day) ) then begin
  satellite_time = float(hour)*3600.0 + float(minutes)*60.0 + float(seconds)
  if ((gsi_start_time lt satellite_time) and (satellite_time lt gsi_end_time)) then begin
    ;
    if ( (-55.0 lt scan_angle_in[i,j]) and (scan_angle_in[i,j] lt 55.0) ) then begin
      ;
      if ( i_begin eq -1 ) then begin
           i_begin = i
      endif
      i_end = i
      ;
          year_array[i,j] = year
         month_array[i,j] = month
           day_array[i,j] = day
          hour_array[i,j] = hour
        minute_array[i,j] = minutes
        second_array[i,j] = seconds

         pxl_lat_out[i,j] =    pxl_lat_in[i,j]
         pxl_lon_out[i,j] =    pxl_lon_in[i,j]
      scan_angle_out[i,j] = scan_angle_in[i,j]
          zenith_out[i,j] =     zenith_in[i,j]
            tb89_out[i,j] =       tb89_in[i,j]
           tb165_out[i,j] =      tb165_in[i,j]
           tb176_out[i,j] =      tb176_in[i,j]
           tb180_out[i,j] =      tb180_in[i,j]
           tb182_out[i,j] =      tb182_in[i,j]
            mask_out[i,j] =       mask_in[i,j]
             lwp_out[i,j] =           lwp[i,j]
             iwp_out[i,j] =           iwp[i,j]
             chi_out[i,j] =           chi[i,j]
        converge_out[i,j] =    noconverge[i,j]
      ;
      print,format='(a,1x,i5,1x,a,1x,i3,3x,a,1x,f7.2)','i',i,'j',j, $
            'scan angle',scan_angle_in[i,j]
      print,format='(a,1x,i4,1x,a,1x,i2,1x,a,1x,i2,1x,a,1x,i2,1x,a,1x,i2,1x,a,1x,i2)', $
            'year =',year,'month =',month,'day =',day,'hour =',hour,'minutes',minutes, $
            'seconds',seconds
    endif
    ;
  endif
  endif
  ;
 endif
endfor
endfor
print,format='(a,1x,i6,1x,a,1x,i6)','i_begin =',i_begin,'; i_end =',i_end
;
;*** No data.
;
if ( i_begin eq -1 ) then begin
  print,format='(a)',' '
  print,format='(a)',' '
  print,format='(a)',' '
  print,format='(a)','No data available for GSI'
  print,format='(a)',' '
  print,format='(a,1x,a,1x,a)','file',tempest_data_out,'not written'
  print,format='(a,1x,a)','in',path_out
endif
;
;*** Write out data in hdf5 format.
;
if ( i_begin ne -1 ) then begin ; i_begin can equal zero. 
  print,' '
  print,format='(a,1x,a)','Writing',tempest_data_out
  print,format='(a,1x,a)','Located in',path_out
  print,' '
  ;
  fid = h5f_create(path_out+tempest_data_out)
  ;
  datatype_id  = h5t_idl_create(pxl_lat_out)
  dataspace_id = h5s_create_simple(size(pxl_lat_out,/dimensions))
  dataset_id   = h5d_create(fid,'pixel latitude',datatype_id,dataspace_id)
  h5d_write,dataset_id,pxl_lat_out
  ;
  datatype_id  = h5t_idl_create(pxl_lon_out)
  dataspace_id = h5s_create_simple(size(pxl_lon_out,/dimensions))
  dataset_id   = h5d_create(fid,'pixel longitude',datatype_id,dataspace_id)
  h5d_write,dataset_id,pxl_lon_out
  ;
  datatype_id  = h5t_idl_create(scan_angle_out)
  dataspace_id = h5s_create_simple(size(scan_angle_out,/dimensions))
  dataset_id   = h5d_create(fid,'scan_angle',datatype_id,dataspace_id)
  h5d_write,dataset_id,scan_angle_out
  ;
  datatype_id  = h5t_idl_create(zenith_out)
  dataspace_id = h5s_create_simple(size(zenith_out,/dimensions))
  dataset_id   = h5d_create(fid,'zenith_angle',datatype_id,dataspace_id)
  h5d_write,dataset_id,zenith_out
  ;
  datatype_id  = h5t_idl_create(tb89_out)
  dataspace_id = h5s_create_simple(size(tb89_out,/dimensions))
  dataset_id   = h5d_create(fid,'Tb 89 GHz',datatype_id,dataspace_id)
  h5d_write,dataset_id,tb89_out
  ;
  datatype_id  = h5t_idl_create(tb165_out)
  dataspace_id = h5s_create_simple(size(tb165_out,/dimensions))
  dataset_id   = h5d_create(fid,'Tb 165 GHz',datatype_id,dataspace_id)
  h5d_write,dataset_id,tb165_out
  ;
  datatype_id  = h5t_idl_create(tb176_out)
  dataspace_id = h5s_create_simple(size(tb176_out,/dimensions))
  dataset_id   = h5d_create(fid,'Tb 176 GHz',datatype_id,dataspace_id)
  h5d_write,dataset_id,tb176_out
  ;
  datatype_id  = h5t_idl_create(tb180_out)
  dataspace_id = h5s_create_simple(size(tb180_out,/dimensions))
  dataset_id   = h5d_create(fid,'Tb 180 GHz',datatype_id,dataspace_id)
  h5d_write,dataset_id,tb180_out
  ;
  datatype_id  = h5t_idl_create(tb182_out)
  dataspace_id = h5s_create_simple(size(tb182_out,/dimensions))
  dataset_id   = h5d_create(fid,'Tb 182 GHz',datatype_id,dataspace_id)
  h5d_write,dataset_id,tb182_out
  ;
  datatype_id  = h5t_idl_create(mask_out)
  dataspace_id = h5s_create_simple(size(mask_out,/dimensions))
  dataset_id   = h5d_create(fid,'Land Mask',datatype_id,dataspace_id)
  h5d_write,dataset_id,mask_out
  ;
  datatype_id  = h5t_idl_create(year_array)
  dataspace_id = h5s_create_simple(size(year_array,/dimensions))
  dataset_id   = h5d_create(fid,'year',datatype_id,dataspace_id)
  h5d_write,dataset_id,year_array
  ;
  datatype_id  = h5t_idl_create(month_array)
  dataspace_id = h5s_create_simple(size(month_array,/dimensions))
  dataset_id   = h5d_create(fid,'month',datatype_id,dataspace_id)
  h5d_write,dataset_id,month_array
  ;
  datatype_id  = h5t_idl_create(day_array)
  dataspace_id = h5s_create_simple(size(day_array,/dimensions))
  dataset_id   = h5d_create(fid,'day',datatype_id,dataspace_id)
  h5d_write,dataset_id,day_array
  ;
  datatype_id  = h5t_idl_create(hour_array)
  dataspace_id = h5s_create_simple(size(hour_array,/dimensions))
  dataset_id   = h5d_create(fid,'hour',datatype_id,dataspace_id)
  h5d_write,dataset_id,hour_array
  ;
  datatype_id  = h5t_idl_create(minute_array)
  dataspace_id = h5s_create_simple(size(minute_array,/dimensions))
  dataset_id   = h5d_create(fid,'minute',datatype_id,dataspace_id)
  h5d_write,dataset_id,minute_array
  ;
  datatype_id  = h5t_idl_create(second_array)
  dataspace_id = h5s_create_simple(size(second_array,/dimensions))
  dataset_id   = h5d_create(fid,'second',datatype_id,dataspace_id)
  h5d_write,dataset_id,second_array
  ;
  datatype_id  = h5t_idl_create(lwp_out)
  dataspace_id = h5s_create_simple(size(lwp_out,/dimensions))
  dataset_id   = h5d_create(fid,'lwp',datatype_id,dataspace_id)
  h5d_write,dataset_id,lwp_out
  ;
  datatype_id  = h5t_idl_create(iwp_out)
  dataspace_id = h5s_create_simple(size(iwp_out,/dimensions))
  dataset_id   = h5d_create(fid,'iwp',datatype_id,dataspace_id)
  h5d_write,dataset_id,iwp_out
  ;
  datatype_id  = h5t_idl_create(chi_out)
  dataspace_id = h5s_create_simple(size(chi_out,/dimensions))
  dataset_id   = h5d_create(fid,'chi',datatype_id,dataspace_id)
  h5d_write,dataset_id,chi_out
  ;
  datatype_id  = h5t_idl_create(converge_out)
  dataspace_id = h5s_create_simple(size(converge_out,/dimensions))
  dataset_id   = h5d_create(fid,'converge',datatype_id,dataspace_id)
  h5d_write,dataset_id,converge_out
  ;
  H5D_CLOSE,dataset_id
  H5S_CLOSE,dataspace_id
  H5T_CLOSE,datatype_id
  H5F_CLOSE,fid
endif
;
;*** The end
;
print,' '
print,' '
print,' '
print,' '
print,'End of code'
print,' '
print,' '
print,' '
print,' '
stop
end
