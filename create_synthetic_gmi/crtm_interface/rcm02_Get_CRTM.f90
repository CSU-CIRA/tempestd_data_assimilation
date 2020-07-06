!============================================================================
! Program to provide an example of CRTM Forward Model function usage for 
! optical depth, TB, and Radiance computations.
!
! CREATION HISTORY:
!    Written by:     Paul van Delst, May 20, 2008
!    Have been modified by: Manajit Sengupta & Louie Grasso at CIRA
!    Modified for CRTM V.2.0.2 and WRF-ARW by: Louie Grasso & Y.Noh
!             June 2011~ (Noh@cira.colostate.edu)
!    Modified for CRTM V.2.0.5 and gfortran compiler (worked with v.4.6.2)
!             Feb 2012 (Noh@cira.colostate.edu)
!    Modified for CRTM V.2.1.3 (channel selection available) with all three
!                 compilers
!             Jan 2014 (Yoo-Jeong.Noh@colostate.edu)
!============================================================================


                      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      !                            !
                      !     November 01 2019       !
                      !                            !
                      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM Get_CRTM
!
!*** Module usage
!
USE CRTM_Module
USE CRTM_Forward_Module
! 
IMPLICIT NONE
!
!*** Parameters
!
CHARACTER(*), PARAMETER :: PROGRAM_NAME  = "Get_CRTM.f90"
CHARACTER(*), PARAMETER :: PATH = "/data2/grasso/small_sat_fv3/2019090312/f003/"
!
!*** Pick a version, then edit Makefile to set up the version you specified here; go on...
!
!CHARACTER(*), PARAMETER :: COEF_PATH = "/home/grasso/rtm/crtm_v2.2.3/Coefficient/"
CHARACTER(*), PARAMETER :: COEF_PATH = "/home/grasso/rtm/crtm_v2.1.3/Coefficients/"
!============================================================================
! 0. **** SOME SATELLITE PARAMETERS ****
!============================================================================
!
! Use the actual band number to select a specific channel in V.2.2.3
!
!                   ======================================
!                   =          CRTM  Variables           =
!                   ======================================
!
!
!
!*** Declare some CRTM parameters
!
INTEGER, PARAMETER          :: N_CLOUDS     = 5
INTEGER, PARAMETER          :: N_PROFILES   = 1  
INTEGER, PARAMETER          :: N_ABSORBERS  = 2
INTEGER, PARAMETER          :: N_AEROSOLS   = 0
INTEGER, PARAMETER          :: N_SENSORS    = 1
!
!CHARACTER(*), PARAMETER     :: SENSOR_TYPE = "geo" ! Geostationary
CHARACTER(*), PARAMETER     :: SENSOR_TYPE = "leo" ! Low Earth Orbiting
!
CHARACTER(*), PARAMETER     :: PATH_ANGLES = "/home/grasso/metimage/fv3_sat_angles/"
CHARACTER(*), PARAMETER     :: PATH_LATLON = "/home/grasso/metimage/fv3_latlon/"
!
!CHARACTER(*), PARAMETER     :: Sensor_Id(N_SENSORS) = (/'ahi_himawari8'/)
!CHARACTER(*), PARAMETER     :: Sensor_Id(N_SENSORS) = (/'v.ahi_himawari8'/)
!CHARACTER(*), PARAMETER     :: Sensor_Id(N_SENSORS) = (/'imgr_g12'/)
!CHARACTER(*), PARAMETER     :: Sensor_Id(N_SENSORS) = (/'avhrr3_n18'/) !noaa-18
!CHARACTER(*), PARAMETER     :: Sensor_Id(N_SENSORS) = (/'abi_gr'/)
!CHARACTER(*), PARAMETER     :: Sensor_Id(N_SENSORS) = (/'v.abi_gr'/)
!CHARACTER(*), PARAMETER     :: Sensor_Id(N_SENSORS) = (/'seviri_m10'/)
!CHARACTER(*), PARAMETER     :: Sensor_Id(N_SENSORS) = (/'v.seviri_m10'/)
!CHARACTER(*), PARAMETER     :: Sensor_Id(N_SENSORS) = (/'tempest-D_cubesat'/) !Tempest-D
CHARACTER(*), PARAMETER     :: Sensor_Id(N_SENSORS) = (/'modis_aqua'/) ! Modis Aqua 
!CHARACTER(*), PARAMETER     :: Sensor_Id(N_SENSORS) = (/'v.modis_aqua'/) ! Modis Aqua reflective 
!
!*** Declare some CRTM Structures
!
TYPE(CRTM_ChannelInfo_type) :: ChannelInfo(N_SENSORS)
TYPE(CRTM_Geometry_type)    :: geo(N_PROFILES)
TYPE(CRTM_Options_type)     :: opt(N_PROFILES)
TYPE(CRTM_Atmosphere_type)  :: Atm(N_PROFILES)
TYPE(CRTM_Surface_type)     :: Sfc(N_PROFILES)
TYPE(CRTM_RTSolution_type), DIMENSION(:,:), ALLOCATABLE :: RTSolution
!
!
!
!                   ======================================
!                   =          Local  Variables          =
!                   ======================================
!
!
!
!
!*** Declare local parameters
!
!*** Begining with v2.1.3, use the satellite band number.
!*** Thus for GOES-R ABI 10.7 is band=13. Previous crtm versions used band=7
!
INTEGER,  PARAMETER :: ICHA= 1            ! Number of channels to process
!
!INTEGER,  PARAMETER :: CHA= 1, CHB= 1    ! tempest-D 89  GHz Tbs (nominal) channel
!INTEGER,  PARAMETER :: CHA= 2, CHB= 2    ! tempest-D 165 GHz Tbs
!INTEGER,  PARAMETER :: CHA= 3, CHB= 3    ! tempest-D 176 GHz Tbs
!INTEGER,  PARAMETER :: CHA= 4, CHB= 4    ! tempest-D 180 GHz Tbs
!INTEGER,  PARAMETER :: CHA= 5, CHB= 5    ! tempest-D 182 GHz Tbs
!
!INTEGER,  PARAMETER :: CHA= 11, CHB= 11    ! HIMAWARI AHI  8.6 um
!INTEGER,  PARAMETER :: CHA= 13, CHB= 13    ! HIMAWARI AHI 10.4 um
!INTEGER,  PARAMETER :: CHA= 14, CHB= 14    ! HIMAWARI AHI 11.2 um
!INTEGER,  PARAMETER :: CHA= 15, CHB= 15    ! HIMAWARI AHI 12.3 um
!INTEGER,  PARAMETER :: CHA= 16, CHB= 16    ! HIMAWARI AHI 13.3 um
!
!INTEGER,  PARAMETER :: CHA= 4, CHB= 4      ! GOES-16 ABI 1.38 um
!INTEGER,  PARAMETER :: CHA= 7, CHB= 7      ! GOES-16 ABI 3.9 um
!INTEGER,  PARAMETER :: CHA= 9, CHB= 9      ! GOES-16 ABI 6.95 um
!INTEGER,  PARAMETER :: CHA= 13, CHB= 13    ! GOES-16 ABI 10.35 um
!
!INTEGER,  PARAMETER :: CHA= 9, CHB= 9      ! SEVIRI 10.8 um
!
!INTEGER,  PARAMETER :: CHA= 19, CHB= 19    ! MODIS aqua 0.940 um
!INTEGER,  PARAMETER :: CHA= 26, CHB= 26    ! MODIS aqua 1.375 um
!INTEGER,  PARAMETER :: CHA= 27, CHB= 27    ! MODIS aqua 6.715 um
INTEGER,  PARAMETER :: CHA= 28, CHB= 28    ! MODIS aqua 7.325 um
!INTEGER,  PARAMETER :: CHA= 31, CHB= 31    ! MODIS aqua 11.03 um
!
INTEGER,  PARAMETER :: SURFACE_EMISS = 0  ! 1 user specified emissivity, 0 not.
INTEGER,  PARAMETER :: SURFACE_ALBEDO = 0 ! 1 user specified albedo, 0 not.
INTEGER,  PARAMETER :: LAND_TYPE = 0      ! 1 user specified land type, 0 not.
INTEGER,  PARAMETER :: MICRO_FLAG = 1     ! 1 use model micro; 0 set all to 0.0
!
!***  "X", "Y", and "Z" points in domain
!
INTEGER,  PARAMETER :: NXP=3072, NYP=1536,N_LAYERS=64
INTEGER,  PARAMETER :: Year=2019, Month=9, Day=3
!
!***  Number of "Z" points in ozone array
!
INTEGER,  PARAMETER :: NKTOP= 45
!
REAL,     PARAMETER :: BTIME= 0.0, ETIME= 0.0, DELTA = 300.0
REAL,     PARAMETER :: DELTAZ = 100.0 ! lowest DELTAZ.
REAL,     PARAMETER :: R_d = 287.05   ! constants
REAL,     PARAMETER :: Cp = 1004.0    ! 
REAL,     PARAMETER :: g_const = 9.81 ! 
!
!*** Declare local integers, reals, and characters
!
INTEGER             :: IRECLEN, IREC
INTEGER             :: ITER, KTOP, NITER, NTFIL, NT1, NT2, NT3, NT4, NT5, NT6
INTEGER             :: Allocate_Status, Error_Status
INTEGER             :: i, ii, j, k, l, m, n_Channels
!
REAL                :: TIME1
!
CHARACTER(1)        :: TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6
CHARACTER(40)       :: TITLE7
!
!*** Declare local 2D arrays
!
REAL, DIMENSION(NXP,NYP)          :: VEG_TYPE, LAND
REAL, DIMENSION(NXP,NYP)          :: CANTMP, LAT, LON, EMISS, ALBEDO
REAL, DIMENSION(NXP,NYP)          :: RAD, TB
REAL, DIMENSION(NKTOP,2)          :: OZONE
REAL, DIMENSION(NXP,NYP)          :: SOL_AZI, SOL_ZEN
REAL, DIMENSION(NXP,NYP)          :: SAT_AZI, SAT_SCAN, SAT_ZEN
!
!*** Declare local 3D arrays
!
REAL(kind=4), DIMENSION(:,:,:), ALLOCATABLE :: PRESS,RVP,TEMP,OPDEPTH,RHO_AIR
!
!*** Cloud type(see NCLOUDS above)               1     2    4     3     5
!
REAL(kind=4), DIMENSION(:,:,:), ALLOCATABLE :: QCLOUD,QICE,QSNOW,QRAIN,QGRAUP
REAL(kind=4), DIMENSION(:,:,:), ALLOCATABLE :: DIAC,  DIAI,DIAS, DIAR, DIAG

integer(kind=4) :: allocation_status
!
REAL(kind=4)    :: Flux_Zenith_Angle
REAL(KIND=4)    :: GEO_LON, SENSOR_ZENITH, SENSOR_AZIMUTH, SENSOR_SCAN
REAL(KIND=4)    :: SOLAR_ZENITH, SOLAR_AZIMUTH, SOLARAZ, SCATANG
REAL(KIND=4)    :: Level_P_Ground, DPRESS
!
!*** Fill ozone array
!
OZONE(:,1) = &
     (/1050.0, 1025.0, 1000.0, 950.0, 920.0, 850.0, 780.0, 700.0, 670.0,  &
        620.0,  570.0,  500.0, 475.0, 430.0, 400.0, 350.0, 300.0, 250.0,  &
        200.0,  150.0,  135.0, 115.0, 100.0,  85.0,  70.0,  60.0,  50.0,  &
         30.0,   25.0,   20.0,  15.0,  10.0,   7.0,   5.0,   4.0,   3.0,  &
          2.0,    1.5,    1.0,   0.5,   0.2,   0.1,   0.05,  0.02,  0.01/)

OZONE(:,2) = &
     (/0.02579, 0.02633, 0.02689, 0.02805, 0.02878, 0.03070, 0.03249,  &
       0.03319, 0.03342, 0.03384, 0.03614, 0.03966, 0.04097, 0.04718,  &
       0.05193, 0.06369, 0.09861, 0.16277, 0.29116, 0.46009, 0.54848,  &
       0.72277, 0.93973, 1.28988, 1.80627, 2.29259, 2.86792, 4.59906,  &
       5.15298, 5.55513, 6.10266, 6.92006, 7.56126, 7.82119, 7.75004,  &
       7.35964, 6.11313, 5.26768, 3.82386, 2.13548, 1.04797, 0.65318,  &
       0.29723, 0.26103, 0.31140/)
!
!*** Allocate 3D arrays.
!
ALLOCATE(PRESS(N_LAYERS,NXP,NYP),     RVP(N_LAYERS,NXP,NYP),  &
          TEMP(N_LAYERS,NXP,NYP),  QCLOUD(N_LAYERS,NXP,NYP),  &
          QICE(N_LAYERS,NXP,NYP),   QSNOW(N_LAYERS,NXP,NYP),  &
         QRAIN(N_LAYERS,NXP,NYP),  QGRAUP(N_LAYERS,NXP,NYP),  &
          DIAC(N_LAYERS,NXP,NYP),    DIAI(N_LAYERS,NXP,NYP),  &
          DIAS(N_LAYERS,NXP,NYP),    DIAR(N_LAYERS,NXP,NYP),  &
          DIAG(N_LAYERS,NXP,NYP), OPDEPTH(N_LAYERS,NXP,NYP),  &
       RHO_AIR(N_LAYERS,NXP,NYP),                             &
       stat=allocation_status)

if ( allocation_status /= 0 ) then
  print("(a)"),"Allocation error for 3d arrays, stopping code at line 167."
  stop
endif
!############################################################################
!#                                                                          #
!#           -- BEGIN LOOP to READ OVER WRF output --                       #
!#                                                                          #
!############################################################################
!============================================================================
! Read WRF ouput files
!============================================================================
!
!*** Indicate CloudCoeff.bin version number and ADA_Module.f90 solar error.
!
!*** rcm02:/home/grasso/rtm/crtm_v2.1.3/REL-2.1.3.ifort/libsrc/ADA_Module.f90
!*** rcm02:/home/grasso/rtm/crtm_v2.1.3/REL-2.1.3.ifort/fix/CloudCoeff/Little_Endian
!
print("(a)")," "
print("(a)")," "
!print("(a)"),"Error: IF( RTV%Visible_Flag_true ) THEN; See ADA_Module.f90 line 393'sh"
print("(a)"),"Fixed: IF( RTV%Solar_Flag_true ) THEN; See ADA_Module.f90 line 393'sh"
print("(a)"),"Look below to see CloudCoeff RELEASE.VERSION: "
print("(a)"),"Look in Get_CRTM.f90 to see absolute paths for ADA and CloudCoeff: Line 217'sh."
print("(a)")," "
print("(a)")," "
!
!*** Calculate the number of files to read.
! 
NITER=INT((ETIME-BTIME)/(DELTA)) + 1
PRINT("(A,1X,I4)"),"Number of model files to process:",NITER
TIME1 = BTIME

MAIN: DO ITER = 1, NITER
 NTFIL=INT(TIME1)
 PRINT("(A,1X,I6)"),"Model time:",NTFIL

 NT6=NTFIL/100000
 NT5=(NTFIL-100000*NT6)/10000
 NT4=(NTFIL-100000*NT6-10000*NT5)/1000
 NT3=(NTFIL-100000*NT6-10000*NT5-1000*NT4)/100
 NT2=(NTFIL-100000*NT6-10000*NT5-1000*NT4-100*NT3)/10
 NT1=(NTFIL-100000*NT6-10000*NT5-1000*NT4-100*NT3-10*NT2)/1
 WRITE(TITLE1,FMT="(I1)")NT1
 WRITE(TITLE2,FMT="(I1)")NT2
 WRITE(TITLE3,FMT="(I1)")NT3
 WRITE(TITLE4,FMT="(I1)")NT4
 WRITE(TITLE5,FMT="(I1)")NT5
 WRITE(TITLE6,FMT="(I1)")NT6
 !!!!!!!!!!!!!!!!!!!!!!!!
 !                      !
 ! Read in WRF 3D data  !
 !                      !
 !!!!!!!!!!!!!!!!!!!!!!!!
 IRECLEN = 4*NXP*NYP
 !
 !*** Read in 3D Pressure (mb)
 !
 TITLE7="PRESS"//TITLE6//TITLE5//TITLE4//TITLE3//TITLE2//TITLE1
 CALL READ_INPUT_3D_DATA(TRIM(PATH),TITLE7,PRESS,N_LAYERS,NXP,NYP,IRECLEN)
 !
 !*** Read in 3D Temperatures (K)
 !
 TITLE7="TEMPK"//TITLE6//TITLE5//TITLE4//TITLE3//TITLE2//TITLE1
 CALL READ_INPUT_3D_DATA(TRIM(PATH),TITLE7,TEMP,N_LAYERS,NXP,NYP,IRECLEN)
 !
 !*** Read in 3D water vapor mixing ratio (g/kg)
 !
 TITLE7="RVP"//TITLE6//TITLE5//TITLE4//TITLE3//TITLE2//TITLE1
 CALL READ_INPUT_3D_DATA(TRIM(PATH),TITLE7,RVP,N_LAYERS,NXP,NYP,IRECLEN)
 !
 !*** Compute moist air density (kg/m**3).
 !
 DO K = 1, N_LAYERS
 DO I = 1, NXP
 DO J = 1, NYP
    RHO_AIR(K,I,J) = 100.0*PRESS(K,I,J)/                             &
                     (287.0*TEMP(K,I,J)*(1.0+0.61*RVP(K,I,J)*0.001))
 ENDDO
 ENDDO
 ENDDO
 !
 CLOUDS: IF ( N_CLOUDS /= 0 ) THEN

 if ( micro_flag == 1 ) then
 !
 !*** Read qcloud:WRF-ARW (g/kg)
 !
 TITLE7="RCP"//TITLE6//TITLE5//TITLE4//TITLE3//TITLE2//TITLE1
 CALL READ_INPUT_3D_DATA(TRIM(PATH),TITLE7,QCLOUD,N_LAYERS,NXP,NYP,IRECLEN)
 !
 !*** Read qice:WRF-ARW (g/kg)
 !
 TITLE7="RPP"//TITLE6//TITLE5//TITLE4//TITLE3//TITLE2//TITLE1
 CALL READ_INPUT_3D_DATA(TRIM(PATH),TITLE7,QICE,N_LAYERS,NXP,NYP,IRECLEN)
 !
 !*** Read qrain:WRF-ARW (g/kg) 
 !
 TITLE7="RRP"//TITLE6//TITLE5//TITLE4//TITLE3//TITLE2//TITLE1
 CALL READ_INPUT_3D_DATA(TRIM(PATH),TITLE7,QRAIN,N_LAYERS,NXP,NYP,IRECLEN)
 !
 !*** Read qsnow:WRF-ARW (g/kg) 
 !
 TITLE7="RAP"//TITLE6//TITLE5//TITLE4//TITLE3//TITLE2//TITLE1
 CALL READ_INPUT_3D_DATA(TRIM(PATH),TITLE7,QSNOW,N_LAYERS,NXP,NYP,IRECLEN)
 !
 !*** Read qgraupel:WRF-ARW (g/kg)
 !
 TITLE7="RGP"//TITLE6//TITLE5//TITLE4//TITLE3//TITLE2//TITLE1
 CALL READ_INPUT_3D_DATA(TRIM(PATH),TITLE7,QGRAUP,N_LAYERS,NXP,NYP,IRECLEN)
 !
 !*** Read in 3D Cloud water diameters (microns)
 !
 TITLE7="DIAC"//TITLE6//TITLE5//TITLE4//TITLE3//TITLE2//TITLE1
 CALL READ_INPUT_3D_DATA(TRIM(PATH),TITLE7,DIAC,N_LAYERS,NXP,NYP,IRECLEN)
 !
 !*** Read in 3D Grauple diameters (microns)
 !
 TITLE7="DIAG"//TITLE6//TITLE5//TITLE4//TITLE3//TITLE2//TITLE1
 CALL READ_INPUT_3D_DATA(TRIM(PATH),TITLE7,DIAG,N_LAYERS,NXP,NYP,IRECLEN)
 !
 !*** Read in 3D qice diameters (microns)
 !
 TITLE7="DIAI"//TITLE6//TITLE5//TITLE4//TITLE3//TITLE2//TITLE1
 CALL READ_INPUT_3D_DATA(TRIM(PATH),TITLE7,DIAI,N_LAYERS,NXP,NYP,IRECLEN)
 !
 !*** Read in 3D Rain water diameters (microns)
 !
 TITLE7="DIAR"//TITLE6//TITLE5//TITLE4//TITLE3//TITLE2//TITLE1
 CALL READ_INPUT_3D_DATA(TRIM(PATH),TITLE7,DIAR,N_LAYERS,NXP,NYP,IRECLEN)
 !
 !*** Read in 3D snow diameters (microns)
 !
 TITLE7="DIAS"//TITLE6//TITLE5//TITLE4//TITLE3//TITLE2//TITLE1
 CALL READ_INPUT_3D_DATA(TRIM(PATH),TITLE7,DIAS,N_LAYERS,NXP,NYP,IRECLEN)
 !
 endif
 if ( micro_flag == 0 ) then
  print("(a)")," "
  print("(a)"),"Setting all condensate stuff to zero."
  print("(a)")," "
  qcloud(:,:,:) = 0.0
  qrain(:,:,:)  = 0.0
  qgraup(:,:,:) = 0.0
  qice(:,:,:)   = 0.0
  qsnow(:,:,:)  = 0.0
  !
  diac(:,:,:) = 0.0
  diag(:,:,:) = 0.0
  diai(:,:,:) = 0.0
  diar(:,:,:) = 0.0
  dias(:,:,:) = 0.0
 endif
 ENDIF  CLOUDS
 !
 !!!!!!!!!!!!!!!!!!!!!!!!
 !                      !
 !   Read in 2D data    !
 !                      !
 !!!!!!!!!!!!!!!!!!!!!!!!
IF(SENSOR_TYPE == "leo")THEN
  SOL_AZI(:,:) = 90.0 ! Sun is to your east.
   PRINT("(A)")," "
   PRINT("(A)"),"Setting solar azimuth to 90 degrees."
   PRINT("(A)"),"That is, the sun is 90 degrees clockwise from north."
   PRINT("(A)")," "
   SOL_ZEN(:,:) = 30.0 ! Sun is 30 degrees from above you.
   PRINT("(A)")," "
   PRINT("(A)"),"Setting solar zenith to 30 degrees."
   PRINT("(A)")," "
   PRINT("(A)")," "
   PRINT("(A)"),"If Sun down you must change these values."
   PRINT("(A)")," "
   PRINT("(A)")," "
ENDIF
!
IF(SENSOR_TYPE == "geo")THEN
 !
 !*** Read in solar azimuth angle
 !
 TITLE7="solar_azimuth.dat"
 CALL READ_INPUT_2D_DATA(TRIM(PATH_ANGLES),TITLE7,SOL_AZI,NXP,NYP,IRECLEN)
 !
 !*** Read in solar zenith angle
 !
 TITLE7="solar_zenith.dat"
 CALL READ_INPUT_2D_DATA(TRIM(PATH_ANGLES),TITLE7,SOL_ZEN,NXP,NYP,IRECLEN)
ENDIF
 !
 !*** Read in satellite azimuth angle
 !
 TITLE7="azimuth.dat"
 CALL READ_INPUT_2D_DATA(TRIM(PATH_ANGLES),TITLE7,SAT_AZI,NXP,NYP,IRECLEN)
 !
 !*** Read in satellite scan angle
 !
 TITLE7="scan.dat"
 CALL READ_INPUT_2D_DATA(TRIM(PATH_ANGLES),TITLE7,SAT_SCAN,NXP,NYP,IRECLEN)
 !
 !*** Read in satellite zenith angle
 !
 TITLE7="zenith.dat"
 CALL READ_INPUT_2D_DATA(TRIM(PATH_ANGLES),TITLE7,SAT_ZEN,NXP,NYP,IRECLEN)
 !
 TITLE7="xlat.dat"
 CALL READ_INPUT_2D_DATA(TRIM(PATH_LATLON),TITLE7,LAT,NXP,NYP,IRECLEN)
 !
 TITLE7="xlon.dat"
 CALL READ_INPUT_2D_DATA(TRIM(PATH_LATLON),TITLE7,LON,NXP,NYP,IRECLEN)
 !
 !*** Read in 2D canopy temperature (K)
 !
 TITLE7="CANTMP"//TITLE6//TITLE5//TITLE4//TITLE3//TITLE2//TITLE1
 CALL READ_INPUT_2D_DATA(TRIM(PATH),TITLE7,CANTMP,NXP,NYP,IRECLEN)
 !
 !*** Read in 2D land type from FV3; liquid water=0, land=1, ice=2
 !
 TITLE7="LAND"//TITLE6//TITLE5//TITLE4//TITLE3//TITLE2//TITLE1
 CALL READ_INPUT_2D_DATA(TRIM(PATH),TITLE7,LAND,NXP,NYP,IRECLEN)

!PRINT("(A)")," "
!PRINT("(A)")," "
!PRINT("(A)"),"No canopy temp for RAMS 27june05."
!PRINT("(A)"),"We'll simply set cantmp(i,j)=temp(1,i,j)."
!PRINT("(A)")," "
!PRINT("(A)")," "

!DO I = 1, NXP
!  DO J = 1, NYP
!    CANTMP(I,J) = TEMP(1,I,J)
!  ENDDO
!ENDDO

 !
 !*** Read in 2D spectrial emissivity (from Eva Borbas at CIMSS)
 !
 IF ( SURFACE_EMISS == 1 )THEN
  !
  IF (Sensor_Id(N_SENSORS) == "ahi_himawari8" .or. Sensor_Id(N_SENSORS) == "abi_gr")THEN
   !
   !*** Begining with v2.1.3, use the satellite band number.
   !*** Thus for GOES-R ABI 10.7 is band=13. Previous crtm versions used band=7
   !
   IF( CHA == 7)  TITLE7="EMIS39.DAT"  
   IF( CHA == 11) TITLE7="EMIS85.DAT"  
   IF( CHA == 13) TITLE7="EMIS1035.DAT"
   IF( CHA == 14) TITLE7="EMIS112.DAT"
   IF( CHA == 15) TITLE7="EMIS123.DAT"
   IF( CHA == 16) TITLE7="EMIS133.DAT"
   !
   CALL READ_INPUT_2D_DATA(TRIM(PATH),TITLE7,EMISS,NXP,NYP,IRECLEN)
  ELSE
   PRINT("(A)")," "
   PRINT("(A,1X,A)"),"Will use CRTM emissivity data for",Sensor_Id(N_SENSORS)
   PRINT("(A)")," "
  ENDIF
  !
 ENDIF
 !
 !*** Read in 2D spectrial MODIS-16-DAY albedo (from Crystal Schaaf at B.U.)
 !
 IF ( SURFACE_ALBEDO == 1 )THEN
!  TITLE7="M2R_RAW_REFL_1.64.DAT" ! /home/grasso/fanyou_goesr/modis_1.64_albedo
   TITLE7="M2R_RAW_REFL_2.1.DAT" ! /home/grasso/fanyou_goesr/modis_1.64_albedo
   CALL READ_INPUT_2D_DATA(TRIM(PATH),TITLE7,ALBEDO,NXP,NYP,IRECLEN)
 ENDIF
 !
 !*** Read in 2D land types
 !
 IF ( LAND_TYPE == 1 )THEN
   TITLE7="caps_vegtype.dat" 
   CALL READ_INPUT_2D_DATA(TRIM(PATH),TITLE7,VEG_TYPE,NXP,NYP,IRECLEN)
 ENDIF
!
!##############################################################################
!
!============================================================================
! 2. **** INITIALIZE THE CRTM ****
!============================================================================
!
! 2a. This initializes the CRTM for the sensors predefined in the example
!     SENSOR_ID parameter.
!
! More data selection options for 'EmisCoeff_File' for V.2.0.2 
! Three Nalli et al (2008)'s data can be found in the CRTM Coefficients
! directory and two options for TauCoeff (ODPS or ODAS). Check the filenames
! and links of the CRTM_Coefficients directory shown as below
! --------------------------------------------------
WRITE( *,'(/5x,"Initializing the CRTM...")' )
Error_Status = CRTM_Init( Sensor_Id, ChannelInfo,                    &
               IRlandCoeff_File ='IGBP.IRland.EmisCoeff.bin',        &
               IRwaterCoeff_File='WuSmith.IRwater.EmisCoeff.bin',    &
               MWwaterCoeff_File='FASTEM4.MWwater.EmisCoeff.bin',    &
               VISlandCoeff_File='IGBP.VISland.EmisCoeff.bin',       &
               File_Path=COEF_PATH)

 IF ( Error_Status /= SUCCESS ) THEN 
   CALL Display_Message( PROGRAM_NAME,"Error initializing CRTM",Error_Status)  
   STOP
 END IF
!------------------------------------------
! 2b. Determine the total number of channels
!     for which the CRTM was initialized
!------------------------------------------
n_Channels = SUM(ChannelInfo%n_Channels)
PRINT("(A)")," "
PRINT("(A)")," "
PRINT("(A,1X,A)"),"Processing data for",Sensor_Id(N_SENSORS)
PRINT("(A,1X,I2)"),"n_channels:",n_Channels
PRINT*,"ChannelInfo(1)%Sensor_Channels:",ChannelInfo(1)%Sensor_Channel
PRINT("(A)")," "
PRINT("(A)")," "
! ----------------------------------
! 2c. Select the channel subset
!     NOTE: Channel list is sorted
!     into ascending order for processing.
! ----------------------------------
  Error_Status = CRTM_ChannelInfo_Subset( ChannelInfo(1), &
                                          Channel_Subset = (/CHA/) ) 
!                               Channel_Subset = (/(ICHA,ICHA=CHA,CHB)/) )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Error subsetting channels', FAILURE )
    STOP
  END IF
!============================================================================
! 3. **** ALLOCATE STRUCTURE ARRAYS ****
!============================================================================
!-----------------------
! 3a. Allocate the ARRAYS
!-----------------------
ALLOCATE( RTSolution( n_Channels, N_PROFILES ), STAT = Allocate_Status )
 IF ( Allocate_Status /= 0 ) THEN 
   CALL Display_Message( PROGRAM_NAME, &
              "Allocation error for RTSolution",Error_Status)  
   STOP
 END IF
!---------------------------
! 3b. Allocate the STRUCTURES
!---------------------------
!
!*** RTSolution structure
!
CALL CRTM_RTSolution_Create(RTSolution, N_LAYERS)
IF ( Any(.NOT. CRTM_RTSolution_Associated( RTSolution )) ) THEN
   CALL Display_Message( PROGRAM_NAME, &
                         "Error with CRTM_RTSolution_Create",Error_Status)  
 STOP
END IF
!
!*** Atmosphere structure
!
CALL CRTM_Atmosphere_Create(Atm,N_LAYERS,N_ABSORBERS,N_CLOUDS,N_AEROSOLS)
IF ( Any(.NOT. CRTM_Atmosphere_Associated( Atm )) ) THEN 
   CALL Display_Message( PROGRAM_NAME, &
                         "Error with CRTM_Atmosphere_Create",Error_Status)  
 STOP
END IF
!
!*** Options structure. Allows a user to use a
!*** surface emissivity dataset, for example.
!
CALL CRTM_Options_Create(opt, n_Channels)
IF (ANY( .NOT. CRTM_Options_Associated(opt)) ) THEN
   CALL Display_Message( PROGRAM_NAME, &
                         "Error with CRTM_Options_Create",Error_Status)  
 STOP
END IF
!============================================================================
! 4. **** ASSIGN INPUT DATA ****
!============================================================================
!---------------
! 4a.1 Profile #1 
!---------------
!
!   What is the CRTM structure/variable name for a user defined
!   surface emissivity (2D file) or surface albedo (2D file)?
!
! We specified emissivity so no need of land type info
!
!!NPOESS Classification Scheme
!!Surface Type Name Classification Index
!compacted soil 1, tilled soil 2, sand 3, rock 4
!irrigated low vegetation 5, meadow grass 6, scrub 7
!broadleaf forest 8, pine forest 9, tundra 10, grass soil 11
!broadleaf pine forest 12, grass scrub 13, soil grass scrub 14
!urban concrete 15, pine brush 16, broadleaf brush 17,
!wet soil 18, scrub soil 19, broadleaf70 pine30 20

IF ( SURFACE_EMISS == 0 .AND. SURFACE_ALBEDO == 0 )THEN
 Sfc(1)%Land_Type = 7
ENDIF
!
Atm(1)%Climatology    = US_STANDARD_ATMOSPHERE
Atm(1)%Absorber_Id    = (/H2O_ID, O3_ID/)
Atm(1)%Absorber_Units = (/MASS_MIXING_RATIO_UNITS,VOLUME_MIXING_RATIO_UNITS/)
!
!*** Specifiy date
!
geo(1)%Year   = Year
geo(1)%Month  = Month
geo(1)%Day    = Day
!CALL rdate_to_jdate(Year,Month,Day,JDAY)  ! jday = 141 of 21 May 2011
!---------------------------------------------------------------------
! Loop through the model horizontal grid points.
!---------------------------------------------------------------------
!
!*** Set the TB and RAD arrays to 0.0
!
TB(:,:)  = 0.0
RAD(:,:) = 0.0
!
!
YPOINTS: DO J = 1, NYP
XPOINTS: DO I = 1, NXP 
  !
!  IF ( LAND_TYPE == 1 )THEN
!    Sfc(1)%Land_Type = VEG_TYPE(I,J)
!  ENDIF
  Sfc(1)%Land_Temperature = CANTMP(I,J)
  !
  if ( land(i,j) == 0 ) then       ! liquid water
    Sfc(1)%Water_Coverage = 1.0_fp
    Sfc(1)%Land_Coverage  = 0.0_fp
  elseif ( land(i,j) == 1 ) then   ! land
    Sfc(1)%Water_Coverage = 0.0_fp
    Sfc(1)%Land_Coverage  = 1.0_fp
  elseif ( land(i,j) == 2 ) then   ! ice
    Sfc(1)%Water_Coverage = 0.0_fp
    Sfc(1)%Land_Coverage  = 1.0_fp
  endif
  !
  !*** Assign satellite and solar angles
  !
  IF(SENSOR_TYPE == "geo")THEN
    !
    !*** CRTM expects 0.0 <= longitude <= 360.0.
    !
    if (lon(i,j) < 0.0) then
        lon(i,j) = lon(i,j) + 360.0 
    endif
    geo(1)%Latitude  = LAT(I,J)
    geo(1)%Longitude = LON(I,J)
    !
    !*** SENSOR_ZENITH max is 80.0 in CRTM.
    !
    if(SAT_ZEN(I,J) > 80.0) THEN
       SAT_ZEN(I,J) = 80.0
    endif
    !
    geo(1)%Sensor_Scan_Angle    = SAT_SCAN(I,J)
    geo(1)%Sensor_Zenith_Angle  = SAT_ZEN(I,J)
    geo(1)%Sensor_Azimuth_Angle = SAT_AZI(I,J)
    geo(1)%Source_Zenith_Angle  = SOL_ZEN(I,J)
    geo(1)%Source_Azimuth_Angle = SOL_AZI(I,J)
    !
  ELSEIF(SENSOR_TYPE == "leo")THEN
    !
    !*** Absolute value of Satellite scan must be less than 80.0.
    !
    if ( sat_scan(i,j) < -75.0 ) then
         sat_scan(i,j) = -75.0            ! sat_scan can be negative for LEO !!!
    endif
    !
    !*** Absolute value of Satellite zenith must be less than 80.0.
    !
    if ( sat_zen(i,j) < -75.0 ) then
         sat_zen(i,j) = -75.0             ! sat_zen can be negative for LEO !!!
    endif
    !
    geo(1)%Sensor_Scan_Angle    = SAT_SCAN(I,J)
    geo(1)%Sensor_Zenith_Angle  = SAT_ZEN(I,J)
    geo(1)%Sensor_Azimuth_Angle = SAT_AZI(I,J)
    geo(1)%Source_Zenith_Angle  = SOL_ZEN(I,J)
    geo(1)%Source_Azimuth_Angle = SOL_AZI(I,J)
    !
  ENDIF
  !
  !*** To allow user defined surface emissivity OR albedo values to be used.
  !
  !*** Must use opt(1)%Emissivity(1) <---- the value of "1" in Emissivity(1) is
  !*** because in CRTM_Forward_Module.f90, line 512, "ln" is the loop variable
  !*** used as the array index for 
  !*** SfcOptics%Emissivity(1,1) = Options(m)%Emissivity(ln),
  !*** line 640. Since we force the CRTM to process just one channel,
  !***  CHA=CHB, thus ln=1 always.
  !
  IF( SURFACE_EMISS == 1 )THEN
     opt(1)%Emissivity(1)  = EMISS(I,J)
     opt(1)%Use_Emissivity = .TRUE.
  ENDIF
  IF( SURFACE_ALBEDO == 1 )THEN
     opt(1)%Use_Emissivity          = .TRUE.
     opt(1)%Use_Direct_Reflectivity = .TRUE.
     opt(1)%Direct_Reflectivity(1)  = ALBEDO(I,J)
     opt(1)%Emissivity(1)           = 1.0 - ALBEDO(I,J)
  ENDIF
  !
  !*** Set n_streams to 6, 8, 16, or 32. See CRTM_CloudScatter.f90
  !*** line 219ish.
  !
  opt(1)%Use_n_Streams = .FALSE.
!  opt(1)%Use_n_Streams = .TRUE.
  opt(1)%n_streams =8 
!  opt(1)%n_streams = 16
!  opt(1)%n_streams = 32
  !
  !---------------------------------------------------------------------
  ! Assign each WRF atmos & cloud ouput to CRTM structure variables
  !---------------------------------------------------------------------
  DO K = 1, N_LAYERS
    Atm(1)%Pressure(N_LAYERS-k+1)    = PRESS(K,I,J)
    Atm(1)%Temperature(N_LAYERS-k+1) = TEMP(K,I,J)
    Atm(1)%Absorber(N_LAYERS-k+1,1)  = RVP(K,I,J)
    !
    IF(RVP(K,I,J) < 0.0)THEN
      Atm(1)%Absorber(N_LAYERS-k+1,1) = 1E-8
    ENDIF
    !
  ENDDO
  !
  !*** Define ground level_pressure.
  !
  Level_P_Ground = 0.5*(PRESS(2,I,J)+PRESS(1,I,J)) *         &
         EXP( (g_const*DELTAZ)/(R_d*TEMP(1,I,J)*(1.0+0.61*RVP(1,I,J)*0.001)) )

  Atm(1)%Level_Pressure(N_LAYERS) = Level_P_Ground
  !
  !*** Define TOA level_pressure, at k = 0, arbitrarily as 0.5*press(n_layers)
  !
  Atm(1)%Level_Pressure(0) = 0.5*PRESS(N_LAYERS,I,J)
  !
  !*** Fill level_pressure between ground and TOA by interpolation
  !*** from layer pressures.
  !
  DO K = 1, N_LAYERS-1
    Atm(1)%Level_Pressure(N_LAYERS-K) = 0.5*(PRESS(K+1,I,J)+PRESS(K,I,J))
  ENDDO
  !
  !*** Fill Reff [microns] and LWP between ground and TOA (layer-defined)
  !
  !*** Multiply qcloud, qrain, ect by 0.001 to convert from g/kg to kg/kg.
  !
MICRO:  DO K = 1, N_LAYERS
  !
  !*** Pressure in Water_Content [kg/m^2] calculations must be in pascals.
  !*** Convert from millibars to pascals by multiplying by 100.0.
  !
  IF(K .EQ. 1) THEN
    DPRESS = 100.0*( Level_P_Ground - 0.5*(PRESS(K+1,I,J)+PRESS(K,I,J)) )
  ELSE IF(K .EQ. N_LAYERS) THEN
    DPRESS = 0.
  ELSE
    DPRESS= 100.0*0.5*( PRESS(K-1,I,J) - PRESS(K+1,I,J) )
  END IF
  !
  IF ( N_CLOUDS /= 0 ) THEN
    !
    Atm(1)%Cloud(1)%Type = WATER_CLOUD
    Atm(1)%Cloud(1)%Effective_Radius(N_LAYERS-K+1) = 0.5*DIAC(K,I,J) 
    Atm(1)%Cloud(1)%Water_Content(N_LAYERS-K+1)= 0.001*QCLOUD(K,I,J) &
                                               * DPRESS/g_const 
   
    Atm(1)%Cloud(2)%Type = ICE_CLOUD
    Atm(1)%Cloud(2)%Effective_Radius(N_LAYERS-K+1) = 0.5*DIAI(K,I,J) 
    Atm(1)%Cloud(2)%Water_Content(N_LAYERS-K+1)= 0.001*QICE(K,I,J)   &
                                               * DPRESS/g_const

    Atm(1)%Cloud(3)%Type = RAIN_CLOUD
    Atm(1)%Cloud(3)%Effective_Radius(N_LAYERS-K+1) = 0.5*DIAR(K,I,J) 
    Atm(1)%Cloud(3)%Water_Content(N_LAYERS-K+1)= 0.001*QRAIN(K,I,J)  &
                                               * DPRESS/g_const 

    Atm(1)%Cloud(4)%Type = SNOW_CLOUD
    Atm(1)%Cloud(4)%Effective_Radius(N_LAYERS-K+1) = 0.5*DIAS(K,I,J) 
    Atm(1)%Cloud(4)%Water_Content(N_LAYERS-K+1)= 0.001*QSNOW(K,I,J)  &
                                               * DPRESS/g_const  

    Atm(1)%Cloud(5)%Type = GRAUPEL_CLOUD
    Atm(1)%Cloud(5)%Effective_Radius(N_LAYERS-K+1) = 0.5*DIAG(K,I,J) 
    Atm(1)%Cloud(5)%Water_Content(N_LAYERS-K+1)= 0.001*QGRAUP(K,I,J) &
                                               * DPRESS/g_const 
    !
    !*** For M-Y or NSSL microphysics
    !
!    IF ( N_CLOUDS == 6 ) THEN
!      Atm(1)%Cloud(6)%Type = HAIL_CLOUD
!      Atm(1)%Cloud(6)%Effective_Radius(N_LAYERS-K+1) = 0.5*DIAH(K,I,J)
!      Atm(1)%Cloud(6)%Water_Content(N_LAYERS-K+1)= 0.001*QHAIL(K,I,J)  &
!                                                        * DPRESS/g_const
    ENDIF
    !
    !*** To make sure habit sizes exist only where habit mass exists.
    !
    IF(ABS(QCLOUD(K,I,J))*RHO_AIR(K,I,J).LE.1E-8) THEN
       Atm(1)%Cloud(1)%Effective_Radius(N_LAYERS-K+1) = 0.0
       Atm(1)%Cloud(1)%Water_Content(N_LAYERS-K+1) = 0.0
    ENDIF
    IF(ABS(QICE(K,I,J))*RHO_AIR(K,I,J).LE.1E-8) THEN
       Atm(1)%Cloud(2)%Effective_Radius(N_LAYERS-K+1) = 0.0
       Atm(1)%Cloud(2)%Water_Content(N_LAYERS-K+1) = 0.0
    ENDIF
    IF(ABS(QRAIN(K,I,J))*RHO_AIR(K,I,J).LE.1E-8) THEN
       Atm(1)%Cloud(3)%Effective_Radius(N_LAYERS-K+1) = 0.0
       Atm(1)%Cloud(3)%Water_Content(N_LAYERS-K+1) = 0.0
    ENDIF
    IF(ABS(QSNOW(K,I,J))*RHO_AIR(K,I,J).LE.1E-8) THEN
       Atm(1)%Cloud(4)%Effective_Radius(N_LAYERS-K+1) = 0.0
       Atm(1)%Cloud(4)%Water_Content(N_LAYERS-K+1) = 0.0
    ENDIF
    IF(ABS(QGRAUP(K,I,J))*RHO_AIR(K,I,J).LE.1E-8) THEN
       Atm(1)%Cloud(5)%Effective_Radius(N_LAYERS-K+1) = 0.0
       Atm(1)%Cloud(5)%Water_Content(N_LAYERS-K+1) = 0.0
    ENDIF
    !
    !*** For M-Y or NSSL microphysics
    !
!    IF ( N_CLOUDS == 6 ) THEN
!      IF(ABS(QHAIL(K,I,J)).LE.1E-8) THEN
!         Atm(1)%Cloud(6)%Effective_Radius(N_LAYERS-K+1) = 0.0
!         Atm(1)%Cloud(6)%Water_Content(N_LAYERS-K+1) = 0.0
!       ENDIF
!    ENDIF
    !
ENDDO MICRO
  ! 
  !*** Interpolate ozone values to model pressure levels.
  !*** Then invert vertical index: K=1 is ground in WRF; K=1 is TOA in CRTM.
  !
  DO K = 1, N_LAYERS
   KTOP = NKTOP-1
   DO II = 1, KTOP
    !
    IF(Atm(1)%Pressure(N_LAYERS-k+1) .LE. OZONE(II,1) .AND.        &
        Atm(1)%Pressure(N_LAYERS-k+1) .GT. OZONE(II+1,1))THEN
        !
        Atm(1)%Absorber(N_LAYERS-k+1,2) = OZONE(II,2) +           &
                       ( (OZONE(II,1)-Atm(1)%Pressure(N_LAYERS-k+1) ) /   &
                         (OZONE(II,1) - OZONE(II+1,1)) ) *                &
                         (OZONE(II+1,2) - OZONE(II,2))
        KTOP = II
    ENDIF
    !
   ENDDO
  ENDDO
  !############################################################################
  !#                                                                          #
  !#                  **** CALL THE CRTM FORWARD MODEL ****                   #
  !#                                                                          #
  !############################################################################
  !call crtm_atmosphere_Inspect(atm(1))
  !read(*,*)
  !
  !*** For leo platforms, call crtm if scan angles are less than 50.0 degrees.
  !
  if (sensor_type == "leo") then
   !
   if( (-50.0 < sat_scan(i,j)) .and. (sat_scan(i,j) < 50.0) )then
    Error_Status =   CRTM_Forward( Atm          , &  
                                   Sfc          , &                            
                                   geo          , &  
                                   ChannelInfo  , & 
                                   RTSolution   , &
                                   Options = opt  )  !For emissivity data
    m = 1  ! Only one profile
    l = 1  ! Only one channel
    !  DO m = 1, N_PROFILES
    DO K = 1, N_LAYERS
     OPDEPTH(N_LAYERS-K+1,I,J) = RTSolution(l,m)%Layer_Optical_Depth(K)
    ENDDO
    TB(I,J)  = RTSolution(l,m)%Brightness_Temperature
    RAD(I,J) = RTSolution(l,m)%Radiance
    !  ENDDO !m
    !
   endif
   !
  endif
  !
  !*** Geostationary
  !
  if (sensor_type == "geo") then
   !
   if( sat_scan(i,j) < 8.6 )THEN
    Error_Status =   CRTM_Forward( Atm          , &  
                                   Sfc          , &                            
                                   geo          , &  
                                   ChannelInfo  , & 
                                   RTSolution   , &
                                   Options = opt  )  !For emissivity data

    m = 1  ! Only one profile
    l = 1  ! Only one channel
    !  DO m = 1, N_PROFILES
    DO K = 1, N_LAYERS
     OPDEPTH(N_LAYERS-K+1,I,J) = RTSolution(l,m)%Layer_Optical_Depth(K)
    ENDDO
    TB(I,J)  = RTSolution(l,m)%Brightness_Temperature
    RAD(I,J) = RTSolution(l,m)%Radiance
    !  ENDDO !m
    !
   endif
   !
  endif
  !
  IF ( Error_Status /= SUCCESS ) THEN 
      CALL Display_Message( PROGRAM_NAME,  &
                           "Error in CRTM RTM_Forward",Error_Status)  
      STOP
  END IF
  !============================================================================
  ! Get OPDEPTH, TB, and RAD from RTSolution structure
  !============================================================================
  !
!  m = 1  ! Only one profile
!  l = 1  ! Only one channel
  !  DO m = 1, N_PROFILES
!  DO K = 1, N_LAYERS
!   OPDEPTH(N_LAYERS-K+1,I,J) = RTSolution(l,m)%Layer_Optical_Depth(K)
!  ENDDO
!  TB(I,J)  = RTSolution(l,m)%Brightness_Temperature
!  RAD(I,J) = RTSolution(l,m)%Radiance
!  !  ENDDO !m
  !
ENDDO XPOINTS
ENDDO YPOINTS
!============================================================================
!------------------------------------------
! DESTROY THE CRTM
!------------------------------------------
WRITE( *, '( /5x, "Destroying the CRTM..." )' )
Error_Status = CRTM_Destroy( ChannelInfo )
!
IF ( Error_Status /= SUCCESS ) THEN
  CALL Display_Message( PROGRAM_NAME,"Error destroying CRTM", Error_Status )
END IF
!-------------------------------------------------------------
! Deallocate the structures.
!-------------------------------------------------------------
CALL CRTM_Atmosphere_Destroy( Atm )
CALL CRTM_RTSolution_Destroy( RTSolution )
CALL CRTM_Options_Destroy( opt )

!-------------------------------------------------------------
! Deallocate the arrays
!-------------------------------------------------------------
DEALLOCATE(RTSolution, STAT = Allocate_Status)


!DEALLOCATE( PRESS,RVP,TEMP )
!DEALLOCATE( QCLOUD,QICE,QSNOW,QRAIN,QGRAUP )
!DEALLOCATE( DIAC,DIAI,DIAS,DIAR,DIAG )

!============================================================================
! Write the results!
!============================================================================
!TITLE7="OPDEPTH"//TITLE6//TITLE5//TITLE4//TITLE3//TITLE2//TITLE1
!PRINT*,"Writing file ",TITLE7
!OPEN(UNIT=34,FILE=TITLE7,FORM='UNFORMATTED',  &
!     ACCESS='DIRECT',STATUS="REPLACE",RECL=IRECLEN)
!IREC = 1
!DO K = 1, N_LAYERS
! WRITE(34,REC=IREC)((OPDEPTH(K,I,J),I=1,NXP),J=1,NYP)
! IREC = IREC + 1
!ENDDO
!CLOSE(34)

!DEALLOCATE( OPDEPTH  )
!
!TITLE7="tb_tempestd_89GHz_noclouds.dat"
!TITLE7="tb_tempestd_165GHz_clouds.dat" 
!TITLE7="tb_tempestd_176GHz_clouds.dat" 
!TITLE7="tb_tempestd_180GHz_clouds.dat" 
!TITLE7="tb_tempestd_182GHz_clouds.dat" 
!
!TITLE7="tb_abi_1035_clouds.dat" 
!
if ( 1 == 1 ) then
  !TITLE7="tb_aqua_6715_clouds.dat" 
  TITLE7="tb_aqua_7325_clouds.dat" 
  !TITLE7="tb_aqua_1103_clouds.dat" 
  !TITLE7="tb_abi_0695_clouds.dat" 
  !TITLE7="tb_seviri_1080_clouds.dat" 
  !TITLE7="tb_himawari_1040_clouds.dat" 
  !
  PRINT("(A)")," "
  PRINT("(A)")," "
  PRINT*,"Writing file ",trim(TITLE7)
  OPEN(UNIT=34,FILE=(TITLE7),FORM='UNFORMATTED',  &
       ACCESS='DIRECT',STATUS="REPLACE",RECL=IRECLEN)
  WRITE(34,REC=1)((TB(I,J),I=1,NXP),J=1,NYP)
  CLOSE(34)
endif

if ( 1 == 2 ) then
  !TITLE7="RAD"//TITLE6//TITLE5//TITLE4//TITLE3//TITLE2//TITLE1
!  TITLE7="rad_aqua_0940_clouds.dat"
  TITLE7="rad_aqua_1375_clouds.dat"
  PRINT*,"Writing file ",TITLE7
  OPEN(UNIT=34,FILE=TITLE7,FORM='UNFORMATTED',  &
      ACCESS='DIRECT',STATUS='REPLACE',RECL=IRECLEN)
  WRITE(34,REC=1)((RAD(I,J),I=1,NXP),J=1,NYP)
  CLOSE(34)
endif
!
!============================================================================
!*** Update time to get next file.
!
TIME1 = TIME1 + DELTA
!
ENDDO MAIN
!
END PROGRAM Get_CRTM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                     !
!                      yjnoh added                                    !
! Change day in date from 'normal' (month/day) to julian day          !
!                                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine rdate_to_jdate (year,month,date,jday)

!      implicit none

!      integer :: year,month,date
!      integer :: iy,jday

!      iy = year - 1900
!      jday= date                                              &
!          +min(1,max(0,month-1))*31                           &
!          +min(1,max(0,month-2))*(28+(1-min(1,mod(iy,4))))    &
!          +min(1,max(0,month-3))*31                           &
!          +min(1,max(0,month-4))*30                           &
!          +min(1,max(0,month-5))*31                           &
!          +min(1,max(0,month-6))*30                           &
!          +min(1,max(0,month-7))*31                           &
!          +min(1,max(0,month-8))*31                           &
!          +min(1,max(0,month-9))*30                           &
!          +min(1,max(0,month-10))*31                          &
!          +min(1,max(0,month-11))*30                          &
!          +min(1,max(0,month-12))*31

!return
!end subroutine rdate_to_jdate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                     !
!                     READ IN 3D INPUT DATA                           !
!                                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE READ_INPUT_3D_DATA(PATH,TITLE6,VAR,NZP,NXP,NYP,LRECLEN)
IMPLICIT NONE
!
!*** Input variables
!
INTEGER, INTENT(IN)  :: NZP, NXP, NYP, LRECLEN
REAL (KIND=4), DIMENSION(NZP,NXP,NYP), INTENT(OUT) :: VAR
CHARACTER (LEN=*), INTENT(IN) :: PATH, TITLE6
!
!*** Local variables
!
Integer             :: IREC, I, J, K
!
PRINT("(T1,A)"),TRIM(PATH)//TRIM(TITLE6)
IREC = 1
OPEN(UNIT=29,FILE=TRIM(PATH)//TRIM(TITLE6),FORM='UNFORMATTED',   &
     ACCESS='DIRECT',STATUS='OLD',RECL=LRECLEN)
DO K = 1, NZP
  READ(29,REC=IREC)((VAR(K,I,J),I=1,NXP),J=1,NYP)
  IREC = IREC + 1
ENDDO
CLOSE(29)
!
END SUBROUTINE READ_INPUT_3D_DATA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                     !
!                     READ IN 2D INPUT DATA                           !
!                                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE READ_INPUT_2D_DATA(PATH,TITLE6,VAR,NXP,NYP,LRECLEN)
IMPLICIT NONE
!
!*** Input variables
!
INTEGER, INTENT(IN)  :: NXP, NYP, LRECLEN
REAL (KIND=4), DIMENSION(NXP,NYP), INTENT(OUT) :: VAR
CHARACTER (LEN=*), INTENT(IN) :: PATH, TITLE6
!
!*** Local variables
!
Integer             :: IREC, I, J
!
PRINT("(T1,A)"),PATH//trim(TITLE6)
IREC = 1
OPEN(UNIT=29,FILE=PATH//TITLE6,FORM='UNFORMATTED',   &
     ACCESS='DIRECT',STATUS='OLD',RECL=LRECLEN)
  READ(29,REC=IREC)((VAR(I,J),I=1,NXP),J=1,NYP)
CLOSE(29)
!
END SUBROUTINE READ_INPUT_2D_DATA


