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
CHARACTER(*), PARAMETER :: PATH = path to data
!
!*** Path to crtm coefficient files
!
CHARACTER(*), PARAMETER :: COEF_PATH = "/scratch2/NCEPDEV/nwprod/NCEPLIBS/fix/crtm_v2.2.6/"
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
CHARACTER(*), PARAMETER     :: SENSOR_TYPE = "leo" ! Low Earth Orbiting
!
CHARACTER(*), PARAMETER     :: PATH_LATLON = PATH 
CHARACTER(*), PARAMETER     :: PATH_ANGLES = PATH
!
CHARACTER(*), PARAMETER     :: Sensor_Id(N_SENSORS) = (/'gmi_gpm'/) 
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
INTEGER,  PARAMETER :: ICHA= 1            ! Number of channels to process
!
INTEGER,  PARAMETER :: CHA= 8, CHB= 8      ! GMI 89 GHz Vertical Polarization
!
INTEGER,  PARAMETER :: SURFACE_EMISS = 0  ! 1 user specified emissivity, 0 not.
INTEGER,  PARAMETER :: SURFACE_ALBEDO = 0 ! 1 user specified albedo, 0 not.
INTEGER,  PARAMETER :: LAND_TYPE = 0      ! 1 user specified land type, 0 not.
INTEGER,  PARAMETER :: MICRO_FLAG = 1     ! 1 use model micro; 0 set all to 0.0
!
!***  "X", "Y", and "Z" points in domain
!
INTEGER,  PARAMETER :: NXP=1536, NYP=768,N_LAYERS=64
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
REAL, DIMENSION(NXP,NYP)          :: VEG_TYPE, LAND, ORBIT_WINDOW
REAL, DIMENSION(NXP,NYP)          :: SNOW_DEPTH,SNOW_COVER, ICE_COVER
REAL, DIMENSION(NXP,NYP)          :: CANTMP, LAT, LON, EMISS, ALBEDO
REAL, DIMENSION(NXP,NYP)          :: RAD, TB
REAL, DIMENSION(NKTOP,2)          :: OZONE
REAL, DIMENSION(NXP,NYP)          :: SOL_AZI, SOL_ZEN
REAL, DIMENSION(NXP,NYP)          :: SAT_AZI, SAT_SCAN, SAT_ZEN
!
!*** Declare local 3D arrays
!
REAL(kind=4), DIMENSION(:,:,:), ALLOCATABLE :: PRESS,RVP,TEMP,OPDEPTH,RHO_AIR
REAL(kind=4), DIMENSION(:,:,:), ALLOCATABLE :: SPH
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
       RHO_AIR(N_LAYERS,NXP,NYP),     SPH(N_LAYERS,NXP,NYP),  &
       stat=allocation_status)

if ( allocation_status /= 0 ) then
  print("(a)"),"Allocation error for 3d arrays, stopping code at line 181."
  stop
endif
!############################################################################
!#                                                                          #
!#           -- BEGIN LOOP to READ OVER output --                           #
!#                                                                          #
!############################################################################
!============================================================================
! Read ouput files
!============================================================================
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
 ! Read in 3D data      !
 !                      !
 !!!!!!!!!!!!!!!!!!!!!!!!
 IRECLEN = 4*NXP*NYP
 !
 !*** Read in 3D Pressure (mb)
 !
 TITLE7="PRESS"//TITLE6//TITLE5//TITLE4//TITLE3//TITLE2//TITLE1
 CALL READ_INPUT_3D_DATA(TRIM(PATH),TITLE7,PRESS,N_LAYERS,NXP,NYP,IRECLEN)
 DO I = 1, NXP
 DO J = 1, NYP
 DO K = 1, N_LAYERS
    IF ( PRESS(K,I,J) < 0.0 ) THEN
       PRINT*,"PRESSURE =",PRESS(K,I,J)
       PRINT*,"AT K",K," I",I," J",J
    ENDIF
 ENDDO
 ENDDO
 ENDDO
 !
 !*** Read in 3D Temperatures (K)
 !
 TITLE7="TEMPK"//TITLE6//TITLE5//TITLE4//TITLE3//TITLE2//TITLE1
 CALL READ_INPUT_3D_DATA(TRIM(PATH),TITLE7,TEMP,N_LAYERS,NXP,NYP,IRECLEN)
 !
 !*** Read in 3D specific humidity (g/g)
 !
 TITLE7="SPH"//TITLE6//TITLE5//TITLE4//TITLE3//TITLE2//TITLE1
 CALL READ_INPUT_3D_DATA(TRIM(PATH),TITLE7,SPH,N_LAYERS,NXP,NYP,IRECLEN)
 !
 !*** Compute moist air density (kg/m**3).
 !
 DO K = 1, N_LAYERS
 DO I = 1, NXP
 DO J = 1, NYP

        RVP(K,I,J) = 1000.0*SPH(K,I,J)/(1.0-SPH(K,I,J)) ! RVP in (g/kg)
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
 TITLE7="CLWMR"//TITLE6//TITLE5//TITLE4//TITLE3//TITLE2//TITLE1
 CALL READ_INPUT_3D_DATA(TRIM(PATH),TITLE7,QCLOUD,N_LAYERS,NXP,NYP,IRECLEN)
 !
 !*** Read qice:WRF-ARW (g/kg)
 !
 TITLE7="ICMR"//TITLE6//TITLE5//TITLE4//TITLE3//TITLE2//TITLE1
 CALL READ_INPUT_3D_DATA(TRIM(PATH),TITLE7,QICE,N_LAYERS,NXP,NYP,IRECLEN)
 !
 !*** Read qrain:WRF-ARW (g/kg) 
 !
 TITLE7="RWMR"//TITLE6//TITLE5//TITLE4//TITLE3//TITLE2//TITLE1
 CALL READ_INPUT_3D_DATA(TRIM(PATH),TITLE7,QRAIN,N_LAYERS,NXP,NYP,IRECLEN)
 !
 !*** Read qsnow:WRF-ARW (g/kg) 
 !
 TITLE7="SNMR"//TITLE6//TITLE5//TITLE4//TITLE3//TITLE2//TITLE1
 CALL READ_INPUT_3D_DATA(TRIM(PATH),TITLE7,QSNOW,N_LAYERS,NXP,NYP,IRECLEN)
 !
 !*** Read qgraupel:WRF-ARW (g/kg)
 !
 TITLE7="GRLE"//TITLE6//TITLE5//TITLE4//TITLE3//TITLE2//TITLE1
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
  print("(a)"),"Setting all condensate zero."
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
 !TITLE7="azimuth.dat"
 !CALL READ_INPUT_2D_DATA(TRIM(PATH_ANGLES),TITLE7,SAT_AZI,NXP,NYP,IRECLEN)
 !
 !*** Read in satellite scan angle
 !
 !TITLE7="scan.dat"
 !CALL READ_INPUT_2D_DATA(TRIM(PATH_ANGLES),TITLE7,SAT_SCAN,NXP,NYP,IRECLEN)
 !
 !*** Read in satellite zenith angle
 !
 !TITLE7="zenith.dat"
 !CALL READ_INPUT_2D_DATA(TRIM(PATH_ANGLES),TITLE7,SAT_ZEN,NXP,NYP,IRECLEN)
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
 !*** Read in 2D snow depth (mm).
 !
 TITLE7="SNOWD"//TITLE6//TITLE5//TITLE4//TITLE3//TITLE2//TITLE1
 CALL READ_INPUT_2D_DATA(TRIM(PATH),TITLE7,SNOW_DEPTH,NXP,NYP,IRECLEN)

 !*** Read in 2D snow coverage, values are 0.0 to 100.0 %
 !
 TITLE7="SNOWC"//TITLE6//TITLE5//TITLE4//TITLE3//TITLE2//TITLE1
 CALL READ_INPUT_2D_DATA(TRIM(PATH),TITLE7,SNOW_COVER,NXP,NYP,IRECLEN)
 !
 !*** Convert snow coverage range to 0.0 to 1.0
 !
 SNOW_COVER(:,:) = SNOW_COVER(:,:)/100.0
 !
 !*** Read in 2D ice coverage, values are 0.0 to 1.0
 !
 TITLE7="ICEC"//TITLE6//TITLE5//TITLE4//TITLE3//TITLE2//TITLE1
 CALL READ_INPUT_2D_DATA(TRIM(PATH),TITLE7,ICE_COVER,NXP,NYP,IRECLEN)
 !
 !*** Read in 2D land type from FV3; liquid water=0, land=1, ice=2
 !
 TITLE7="LAND"//TITLE6//TITLE5//TITLE4//TITLE3//TITLE2//TITLE1
 CALL READ_INPUT_2D_DATA(TRIM(PATH),TITLE7,LAND,NXP,NYP,IRECLEN)
 !
 !*** Read in 2D LEO orbit window: 0 no image, 1 image
 !
 TITLE7="gmi_orbit_window.dat"
 CALL READ_INPUT_2D_DATA(TRIM(PATH),TITLE7,ORBIT_WINDOW,NXP,NYP,IRECLEN)
 !
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
   TITLE7="M2R_RAW_REFL_2.1.DAT"
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
               MWwaterCoeff_File='FASTEM6.MWwater.EmisCoeff.bin',    &
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
  Sfc(1)%Land_Temperature = CANTMP(I,J)
  Sfc(1)%snow_Temperature = CANTMP(I,J)
  Sfc(1)%ice_Temperature  = CANTMP(I,J)
  !
  if ( land(i,j) == 0 ) then       ! liquid water
    Sfc(1)%Water_Coverage = 1.0_fp
    Sfc(1)%Land_Coverage  = 0.0_fp
    Sfc(1)%ice_Coverage   = 0.0_fp
    Sfc(1)%snow_Coverage  = 0.0_fp
  elseif ( land(i,j) == 1 ) then   ! land
    Sfc(1)%Water_Coverage = 0.0_fp
    Sfc(1)%Land_Coverage  = 1.0_fp-SNOW_COVER(I,J)
    Sfc(1)%ice_Coverage   = 0.0_fp
    Sfc(1)%snow_Coverage  = SNOW_COVER(I,J)
    Sfc(1)%snow_Depth     = SNOW_DEPTH(I,J)
  elseif ( land(i,j) == 2 ) then   ! ice
    Sfc(1)%Water_Coverage = 1.0_fp-ICE_COVER(I,J)
    Sfc(1)%Land_Coverage  = 0.0_fp
    Sfc(1)%ice_Coverage   = ICE_COVER(I,J)
    Sfc(1)%snow_Coverage  = 0.0_fp
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
         sat_scan(i,j) = -75.0            ! sat_scan can be negative for LEO
    endif
    !
    !*** Absolute value of Satellite zenith must be less than 80.0.
    !
    if ( sat_zen(i,j) < -75.0 ) then
         sat_zen(i,j) = -75.0             ! sat_zen can be negative for LEO
    endif
    !
    !*** GMI_GPM settings
    ! 
    geo(1)%Sensor_Scan_Angle    = 48.5           ! for conical scanning e.g. GMI
    geo(1)%Sensor_Zenith_Angle  = 53.0           ! for conical scanning e.g. GMI
    geo(1)%Sensor_Azimuth_Angle = 90.0           ! Test for GMI
    geo(1)%Source_Zenith_Angle  = SOL_ZEN(I,J)
    geo(1)%Source_Azimuth_Angle = SOL_AZI(I,J)
    !
  ENDIF
  !
  !*** To allow user defined surface emissivity OR albedo values to be used.
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
  opt(1)%Use_n_Streams = .FALSE.
  opt(1)%n_streams =8 
  !
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
ENDDO MICRO
  ! 
  !*** Interpolate ozone values to model pressure levels.
  !*** Then invert vertical index: K=1 is ground in model; K=1 is TOA in CRTM.
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
  !
  !*** For leo platforms, call crtm if scan angles are less than 50.0 degrees.
  !
  if (sensor_type == "leo") then
   !
   if( (-50.0 < sat_scan(i,j)) .and. (sat_scan(i,j) < 50.0) )then
   if( orbit_window(i,j) == 1.0 )then 
    Error_Status =   CRTM_Forward( Atm          , &  
                                   Sfc          , &                            
                                   geo          , &  
                                   ChannelInfo  , & 
                                   RTSolution   , &
                                   Options = opt  )  !For emissivity data
    m = 1
    l = 1
    DO K = 1, N_LAYERS
     OPDEPTH(N_LAYERS-K+1,I,J) = RTSolution(l,m)%Layer_Optical_Depth(K)
    ENDDO
    TB(I,J)  = RTSolution(l,m)%Brightness_Temperature
    RAD(I,J) = RTSolution(l,m)%Radiance
    !
   endif
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

    m = 1
    l = 1
    DO K = 1, N_LAYERS
     OPDEPTH(N_LAYERS-K+1,I,J) = RTSolution(l,m)%Layer_Optical_Depth(K)
    ENDDO
    TB(I,J)  = RTSolution(l,m)%Brightness_Temperature
    RAD(I,J) = RTSolution(l,m)%Radiance
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
!
!============================================================================
! Write results
!============================================================================
TITLE7="tb_gmi_89VGhz"
!
PRINT("(A)")," "
PRINT("(A)")," "
PRINT*,"Writing file ",trim(TITLE7)
OPEN(UNIT=34,FILE=(TITLE7),FORM='UNFORMATTED',  &
     ACCESS='DIRECT',STATUS="REPLACE",RECL=IRECLEN)
WRITE(34,REC=1)((TB(I,J),I=1,NXP),J=1,NYP)
CLOSE(34)
!
TITLE7="rad_aqua_1375_clouds.dat"
PRINT*,"Writing file ",TITLE7
OPEN(UNIT=34,FILE=TITLE7,FORM='UNFORMATTED',  &
    ACCESS='DIRECT',STATUS='REPLACE',RECL=IRECLEN)
WRITE(34,REC=1)((RAD(I,J),I=1,NXP),J=1,NYP)
CLOSE(34)
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
