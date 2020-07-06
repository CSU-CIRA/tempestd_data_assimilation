PROGRAM GFDL_DIAMS
IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                         !
! Code to compute particle size from GFDL microphysics scheme in FV**3.   !
! 15 January 2019                                                         !
!                                                                         !
! See: Lin, Farley, Orville, 1983, J. Climate Appl. Meteo., vol 22        !
!      page 1068, Eqs. 4-6.                                               !
!                                                                         !
!   N(D) = No*exp(-lamda*D)                                               !
!                                                                         !
!   All particles are assumed to be spherical.                            !
!                                                                         !
!   m = rho*V = 1/6*( pi*rho*D**3 )                                       !
!                                                                         !
!   q = int(o,inf)m*N(D)dD = 1/6*pi*rho*No*int(o,inf)D**3 exp(-lambda*D)  !
!                                                                         !
!   After many integration by parts   q = (1/lambda**4)*(pi*rho*No)       !
!                                                                         !
!   lambda = ( pi*rho*No/q )**0.25, let q = rho_air*(q*), then            !
!   lambda = (pi*rho*No/(rho_air*(q*)) )**0.25                            !
!                                                                         !
!   Deff = ( int(o,inf)D**3 *N(D)dD )/( int(o,inf)D**2 *N(D)dD )          !
!                                                                         !
!   int(o,inf)D**3 *N(D)dD = 6*No/lambda**4                               !
!   int(o,inf)D**2 *N(D)dD = 2*No/lambda**3                               !
!                                                                         !
!                        A general formula                                !
!                                                                         !
!   Deff=3/lambda, where lambda=( pi*rho_c*Nco/(rho_air*(q*_c)) )**0.25   !
!                                                                         !
!   Units:                                                                !
!          q*  .......... habit mass mixing ratio (kg_cond/kg_air)        !
!          Nco .......... intercept (m**-4)                               !
!          rho_air ...... density of MOIST air (kg_air/m**3_air)          !
!          rho_c ........ density of habit (kg/m**3)                      !
!                                                                         !
!                                                                         !
!                  / rho_r = 1000 kg/m**3                                 !
!          rho_c = | rho_s = 100  kg/m**3                                 !
!                  \ rho_h = 917  kg/m**3                                 !
!                                                                         !
!                  / Nro =  8.0e6 m**-4                                   !
!          Nco =   | Nso =  3.0e6 m**-4                                   !
!                  \ Nho =  3.0e4 m**-4                                   !
!                                                                         !
!  Note: N(D) = No*exp(-lamda*D) is for precipitation sized particles.    !
!        What about the effective size of cloud drops and cloud ice???    !
!        Easy Sneezy...We'll follow Ferrier and Aligo from cloud_efr.f90  !
!                                                                         !
!                                                                         !
!  Cloud Ice: Line 320'sh                                                 !
!                                                                         !
!               efr_qi =  75 um                                           !
!             Deff_ice = 150 um (seems large)                             !
!                                                                         !
!  Cloud droplets: Line 312'sh                                            !
!                                                                         !
!           tem4 = max(zero,(t0c-T1d)*r0_05)                              !
!           indexw = five + five*min(one,tem4)                            !
!           efr_ql = 1.5*indexw                                           !
!                                                                         !
!    Constants defined in "use constants", line 28 in cloud_efr.f90       !
!                                                                         !
!           zero = 0.0                                                    !
!           one  = 1.0                                                    !
!           five = 5.0                                                    !
!           t0c  = 273.15                                                 !
!           r0_05 = 0.05                                                  !
!                                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!*** Model dimensions
!
character(len=80), parameter :: data="/data2/grasso/small_sat_nssl_wrf/"
character(len=6),  parameter :: time="001800"
INTEGER(kind=4),   parameter :: NZP=34, NXP=1199, NYP=799
INTEGER(kind=4),   parameter :: LRECLEN = 4*NXP*NYP
!
!***  LOCAL VARIABLES
!
INTEGER(kind=4) :: i,j,k, IREC
INTEGER(kind=4) :: ALLOCATION_STATUS
!
!*** Define PI
!
REAL(kind=4),  PARAMETER :: PI=4*ATAN(1.0)
!
!*** Minimum habit mass mixing ratio(kg/kg)
!
REAL(kind=4),  PARAMETER :: qmin = 1.0E-11
!
!*** Minimum habit size (meters)
!
REAL(kind=4),  PARAMETER :: min_c=2.0e-6, min_r=20.0e-6
REAL(kind=4),  PARAMETER :: min_s=20.0e-6, min_h=20.0e-6
!
!*** Maximum habit size (meters)
!
REAL(kind=4),  PARAMETER :: max_c=1.0e-2, max_r=1.0e-2
REAL(kind=4),  PARAMETER :: max_s=2.0e-2,  max_h=5.0e-0
!
!*** Intercept values (meters**-4)
!
real(kind=4) :: Nro=8.0E6, Nso=3.0E6, Nho=3.0E4
!
!*** Habit density (kg/m^3)
!
real, parameter :: rhor=1000.0, rhos=100.0, rhoh=917.0
!
!*** From Ferrier and Aligo micro in cloud_efr.f90
!
real(kind=4) :: tem4, indexw, efr_ql
!
!*** Lambda from GFDL MICRO
!
real(kind=4) :: lambda
!
REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE :: MASS, DIAM
REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE :: RVP, SPH, PRESS, TEMPK, RHO_AIR
!
!*** Allocate MASS
!
ALLOCATE(MASS(NZP,NXP,NYP), STAT=ALLOCATION_STATUS)
IF( ALLOCATION_STATUS /= 0 )THEN
  PRINT("(A)"),"Allocation error; stopping code at line 119."
  STOP
ENDIF
!
!*** Allocate 3-D variables
!
ALLOCATE(RVP(NZP,NXP,NYP), SPH(NZP,NXP,NYP), PRESS(NZP,NXP,NYP), &
         TEMPK(NZP,NXP,NYP), RHO_AIR(NZP,NXP,NYP), STAT=ALLOCATION_STATUS)
IF( ALLOCATION_STATUS /= 0 )THEN
  PRINT("(A)"),"Allocation error; stopping code at line 65."
  STOP
ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                               !
!  Read in pressure (mb), temperature (K),      !
!  specific humidity (g/g).                     !
!                                               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!*** Read in FV**3 pressure (mb)
!
print("(A,1X,A)"),"Reading",trim(data)//"PRESS"//time
IREC = 1
OPEN(UNIT=29,FILE=trim(data)//"PRESS"//time,FORM="UNFORMATTED",  &
     ACCESS="DIRECT",STATUS="OLD",RECL=LRECLEN)
DO K = 1, NZP
 READ(29,REC=IREC)((PRESS(K,I,J),I=1,NXP),J=1,NYP)
 IREC = IREC + 1
ENDDO
CLOSE(29)
!
!*** Read in FV**3 temperature (K)
!
print("(A,1X,A)"),"Reading",trim(data)//"TEMPK"//time
IREC = 1
OPEN(UNIT=29,FILE=trim(data)//"TEMPK"//time,FORM="UNFORMATTED", &
     ACCESS="DIRECT",STATUS="OLD",RECL=LRECLEN)
DO K = 1, NZP
 READ(29,REC=IREC)((TEMPK(K,I,J),I=1,NXP),J=1,NYP)
 IREC = IREC + 1
ENDDO
CLOSE(29)
!
!*** Read in FV**3 specific humidity (g/g).
!
print("(A,1X,A)"),"Reading",trim(data)//"SPH"//time
IREC = 1
OPEN(UNIT=29,FILE=trim(data)//"SPH"//time,FORM="UNFORMATTED", &
     ACCESS="DIRECT",STATUS="OLD",RECL=LRECLEN)
DO K = 1, NZP
 READ(29,REC=IREC)((SPH(K,I,J),I=1,NXP),J=1,NYP)
 IREC = IREC + 1
ENDDO
CLOSE(29)
!
!*** Compute moist air density (kg/m^3).
!
do k = 1,nzp
do i = 1,nxp
do j = 1,nyp
 !
 !*** Convert specific humidity to water vapor mixing ratio (g/kg).
 !
 rvp(k,i,j) = 1000.0*sph(k,i,j)/(1.0-sph(k,i,j))
 !
 !*** Compute moist air density
 !
 rho_air(k,i,j) = 100.0*press(k,i,j)/ &
                  (287.0*tempk(k,i,j)*(1.0 + 0.61*rvp(k,i,j)*0.001))
 !
enddo
enddo
enddo
!
!*** Deallocate water vapor mixing ratio and pressure arrays.
!
DEALLOCATE(RVP, SPH, PRESS)
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                               !
!          Compute CLOUD diameters              !
!                                               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!*** Read in FV**3 cloud water mass mixing ratio (g/kg).
!
print("(A)")," "
print("(A,1X,A)"),"Reading",trim(data)//"RCP"//time
IREC = 1
OPEN(UNIT=29,FILE=trim(data)//"RCP"//time,FORM="UNFORMATTED", &
     ACCESS="DIRECT",STATUS="OLD",RECL=LRECLEN)
DO K = 1, NZP
 READ(29,REC=IREC)((MASS(K,I,J),I=1,NXP),J=1,NYP)
 IREC = IREC + 1
ENDDO
CLOSE(29)
!
diam(:,:,:) = 0.0
!
!*** CLOUD DIAMETER CALCULATIONS
!
do k = 1,nzp
do i = 1,nxp
do j = 1,nyp
 !
 !*** Convert mass(k,i,j) from g/kg to kg/kg.
 !
 mass(k,i,j) = mass(k,i,j) * 0.001
 ! 
 !*** tem4 = max(0.0, (273.15-tempk(k,i,j)*0.05)
 !*** indexw = 5.0 + 5.0*min(1.0,tem4)
 !*** efr_ql = 1.5*indexw
 !*** diam = 2.0 * efr_ql
 !
 IF( MASS(K,I,J) > qmin )THEN
   !
   tem4 = max(0.0, (273.15-tempk(k,i,j)*0.05) )
   indexw = 5.0 + 5.0*min(1.0,tem4)
   efr_ql = 1.5*indexw
   diam(k,i,j) = 2.0 * efr_ql
   !
   !*** Constrain diameter size range between min_c = 2 um and max_c = 1 cm.
   !
   IF ( DIAM(K,I,J) > max_c ) THEN
     DIAM(K,I,J) = max_c
   ENDIF
   IF ( DIAM(K,I,J) < min_c ) THEN
     DIAM(K,I,J) = min_c
   ENDIF
 ENDIF
 !
 !*** Convert diameter from meters to microns.
 !
 diam(k,i,j) = 1.0e6 * diam(k,i,j)
 !
enddo
enddo
enddo
!
!*** Write out cloud droplet diameters in microns.
!
print("(A,1X,A)"),"Writing",trim(data)//"DIAC"//time
IREC = 1
OPEN(UNIT=29,FILE=trim(data)//"DIAC"//time,FORM="UNFORMATTED", &
     ACCESS="DIRECT",STATUS="REPLACE",RECL=LRECLEN)
  DO K = 1, NZP
   WRITE(29,REC=IREC)((DIAM(K,I,J),I=1,NXP),J=1,NYP)
   IREC = IREC + 1
  ENDDO
CLOSE(29)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                               !
!          Compute RAIN diameters               !
!                                               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!*** Read in FV**3 qrain: rain water mass mixing ratio (g/kg).
!
print("(A)")," "
print("(A,1X,A)"),"Reading",trim(data)//"RRP"//time
IREC = 1
OPEN(UNIT=29,FILE=trim(data)//"RRP"//time,FORM="UNFORMATTED", &
     ACCESS="DIRECT",STATUS="OLD",RECL=LRECLEN)
DO K = 1, NZP
 READ(29,REC=IREC)((MASS(K,I,J),I=1,NXP),J=1,NYP)
 IREC = IREC + 1
ENDDO
CLOSE(29)
!
diam(:,:,:) = 0.0
!
!*** RAIN DIAMETER CALCULATIONS
!
do k = 1,nzp
do i = 1,nxp
do j = 1,nyp
 !
 !*** Convert mass(k,i,j) from g/kg to kg/kg.
 !
 mass(k,i,j) = mass(k,i,j) * 0.001
 !
 !*** mass(k,i,j) and qmin have units of kg/kg.
 !
 IF( MASS(K,I,J) > qmin )THEN
   !
   !*** Compute Lambda and diam for rain.
   !
   lambda = ( (pi* rhor * Nro)/( mass(k,i,j) * rho_air(k,i,j) ) )**0.25 
   diam(k,i,j) = 3.0 / lambda
   !
   !*** Constrain diameter size range between min_r = 20 um and max_r = 1 cm.
   !
   IF ( DIAM(K,I,J) > max_r ) THEN
     DIAM(K,I,J) = max_r
   ENDIF
   IF ( DIAM(K,I,J) < min_r ) THEN
     DIAM(K,I,J) = min_r
   ENDIF
 ENDIF
 !
 !*** Convert diameter from meters to microns.
 !
 diam(k,i,j) = 1.0e6 * diam(k,i,j)
 !
enddo
enddo
enddo
!
!*** Write out rain droplet diameters in microns.
!
print("(A,1X,A)"),"Writing",trim(data)//"DIAR"//time
IREC = 1
OPEN(UNIT=29,FILE=trim(data)//"DIAR"//time,FORM="UNFORMATTED", &
     ACCESS="DIRECT",STATUS="REPLACE",RECL=LRECLEN)
  DO K = 1, NZP
   WRITE(29,REC=IREC)((DIAM(K,I,J),I=1,NXP),J=1,NYP)
   IREC = IREC + 1
  ENDDO
CLOSE(29)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                               !
!           Compute ICE diameters               !
!                                               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!*** Read in FV**3 qice: ice mass mixing ratio (g/kg).
!
print("(A)")," "
print("(A,1X,A)"),"Reading",trim(data)//"RPP"//time
IREC = 1
OPEN(UNIT=29,FILE=trim(data)//"RPP"//time,FORM="UNFORMATTED", &
     ACCESS="DIRECT",STATUS="OLD",RECL=LRECLEN)
DO K = 1, NZP
 READ(29,REC=IREC)((MASS(K,I,J),I=1,NXP),J=1,NYP)
 IREC = IREC + 1
ENDDO
CLOSE(29)
!
diam(:,:,:) = 0.0
!
!*** ICE DIAMETER CALCULATIONS
!
do k = 1,nzp
do i = 1,nxp
do j = 1,nyp
 !
 !*** Convert mass(k,i,j) from g/kg to kg/kg.
 !
 mass(k,i,j) = mass(k,i,j) * 0.001
 !
 !*** mass(k,i,j) and qmin have units of kg/kg.
 !
 IF( MASS(K,I,J) > qmin )THEN
   !
   DIAM(K,I,J) = 150.0 
   !
 ENDIF
 !
enddo
enddo
enddo
!
!*** Write out ice diameters in microns.
!
print("(A,1X,A)"),"Writing",trim(data)//"DIAI"//time
IREC = 1
OPEN(UNIT=29,FILE=trim(data)//"DIAI"//time,FORM="UNFORMATTED", &
     ACCESS="DIRECT",STATUS="REPLACE",RECL=LRECLEN)
  DO K = 1, NZP
   WRITE(29,REC=IREC)((DIAM(K,I,J),I=1,NXP),J=1,NYP)
   IREC = IREC + 1
  ENDDO
CLOSE(29)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                               !
!           Compute SNOW diameters              !
!                                               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!*** Read in FV**3 qsnow: snow mass mixing ratio (g/kg).
!
print("(A)")," "
print("(A,1X,A)"),"Reading",trim(data)//"RAP"//time
IREC = 1
OPEN(UNIT=29,FILE=trim(data)//"RAP"//time,FORM="UNFORMATTED", &
     ACCESS="DIRECT",STATUS="OLD",RECL=LRECLEN)
DO K = 1, NZP
 READ(29,REC=IREC)((MASS(K,I,J),I=1,NXP),J=1,NYP)
 IREC = IREC + 1
ENDDO
CLOSE(29)
!
diam(:,:,:) = 0.0
!
!*** SNOW DIAMETER CALCULATIONS
!
do k = 1,nzp
do i = 1,nxp
do j = 1,nyp
 !
 !*** Convert mass(k,i,j) from g/kg to kg/kg.
 !
 mass(k,i,j) = mass(k,i,j) * 0.001
 !
 !*** mass(k,i,j) and qmin have units of kg/kg.
 !
 if( mass(k,i,j) > qmin )then
   !
   !*** Compute Lambda and diam for snow.
   !
   lambda = ( (pi* rhos * Nso)/( mass(k,i,j) * rho_air(k,i,j) ) )**0.25 
   diam(k,i,j) = 3.0 / lambda
   !
   !*** Constrain diameter size range between min_s = 20 um and max_s = 2 cm
   !
   IF ( DIAM(K,I,J) > max_s ) THEN
     DIAM(K,I,J) = max_s
   ENDIF
   IF ( DIAM(K,I,J) < min_s ) THEN
     DIAM(K,I,J) = min_s
   ENDIF
 ENDIF
 !
 !*** Convert diameter from meters to microns.
 !
 diam(k,i,j) = 1.0e6 * diam(k,i,j)
 !
enddo
enddo
enddo
!
!*** Write out snow diameters in microns.
!
print("(A,1X,A)"),"Writing",trim(data)//"DIAS"//time
IREC = 1
OPEN(UNIT=29,FILE=trim(data)//"DIAS"//time,FORM="UNFORMATTED", &
     ACCESS="DIRECT",STATUS="REPLACE",RECL=LRECLEN)
  DO K = 1, NZP
   WRITE(29,REC=IREC)((DIAM(K,I,J),I=1,NXP),J=1,NYP)
   IREC = IREC + 1
  ENDDO
CLOSE(29)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                               !
!        Compute GRAUPEL/HAIL diameters         ! 
!                                               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!*** Read in FV**3 qgraup/hail: graupel/hail mass mixing ratio (g/kg).
!*** There is only one habit here, no convention as to a name.
!
print("(A)")," "
print("(A,1X,A)"),"Reading",trim(data)//"RGP"//time
IREC = 1
OPEN(UNIT=29,FILE=trim(data)//"RGP"//time,FORM="UNFORMATTED", &
     ACCESS="DIRECT",STATUS="OLD",RECL=LRECLEN)
DO K = 1, NZP
 READ(29,REC=IREC)((MASS(K,I,J),I=1,NXP),J=1,NYP)
 IREC = IREC + 1
ENDDO
CLOSE(29)
!
diam(:,:,:) = 0.0
!
!*** GRAUPEL DIAMETER CALCULATIONS: Single moment
!
do k = 1,nzp
do i = 1,nxp
do j = 1,nyp
 !
 !*** Convert mass(k,i,j) from g/kg to kg/kg.
 !
 mass(k,i,j) = mass(k,i,j) * 0.001
 !
 !*** mass(k,i,j) and qmin have units of kg/kg.
 !
 if( mass(k,i,j) > qmin )then
   !
   !*** Compute Lambda and diam for graupel/hail.
   !
   lambda = ( (pi* rhoh * Nho)/( mass(k,i,j) * rho_air(k,i,j) ) )**0.25 
   diam(k,i,j) = 3.0 / lambda
   !
   !*** Constrain diameter size range between min_g = 20 um and max_h = 5 m ???
   !
   IF ( DIAM(K,I,J) > max_h ) THEN
     DIAM(K,I,J) = max_h
   ENDIF
   IF ( DIAM(K,I,J) < min_h ) THEN
     DIAM(K,I,J) = min_h
   ENDIF
 endif
 !
 !*** Convert diameter from meters to microns.
 !
 diam(k,i,j) = 1.0e6 * diam(k,i,j)
 !
enddo
enddo
enddo
!
!*** Write out graupel/hail diameters in microns.
!
print("(A,1X,A)"),"Writing",trim(data)//"DIAG"//time
IREC = 1
OPEN(UNIT=29,FILE=trim(data)//"DIAG"//time,FORM="UNFORMATTED", &
     ACCESS="DIRECT",STATUS="REPLACE",RECL=LRECLEN)
  DO K = 1, NZP
   WRITE(29,REC=IREC)((DIAM(K,I,J),I=1,NXP),J=1,NYP)
   IREC = IREC + 1
  ENDDO
CLOSE(29)
!
!*** All done.
!
print("(a)")," "
print("(a)")," "
print("(a)"),"Hurray, you made it to the end of the code."
print("(a)")," "
print("(a)"),"Hold on...that doesn't mean the code was successful. "
print("(a)")," "
print("(a)")," "
!
END PROGRAM GFDL_DIAMS
