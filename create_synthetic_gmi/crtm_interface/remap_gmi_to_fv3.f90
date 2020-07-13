program remap_gmi_to_fv3
implicit none
!
!***
!
integer(kind=4), parameter :: nxp=221, nyp=2963 ! GMI Orbit dimensions
real(kind=4), dimension(nxp,nyp) :: tb, lat, lon ! GMI lat and lon;  0.0 <lon<=360.0
real(kind=4), parameter :: rect_lat_begin = -89.91032, rect_lon_begin = 0.0
real(kind=4), parameter :: rect_delta_lat=0.234375, rect_delta_lon=0.234375
integer(kind=4), parameter :: rect_nxp=1536, rect_nyp=768 ! rectilinear grid dimensions
real(kind=4), dimension(rect_nxp,rect_nyp) :: rect_lat, rect_lon, rect_tb
!
integer(kind=4) :: reclen, i, j, rect_i, rect_j
!
!***
!
print*," "
print*," "
PRINT("(A,1X,F9.5,A,1X,F6.1)"),"RECT_LAT_BEGIN:",RECT_LAT_BEGIN,  &
                             ", RECT_LON_BEGIN:",RECT_LON_BEGIN
PRINT("(A,1X,F8.6,A,1X,F8.6)"),"RECT_DELTA_LAT:",RECT_DELTA_LAT, &
                             ", RECT_DELTA_LON:",RECT_DELTA_LON
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                     !
!           Read tbs or reflectances, lat, and lon values             !
!                                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
RECLEN=4*NXP*NYP
PRINT*," "
PRINT*, "Reading: GMI_lat.dat"
OPEN(UNIT=29,FILE="GMI_lat.dat",FORM="UNFORMATTED", &
     ACCESS="DIRECT",STATUS="OLD",RECL=RECLEN)
READ(29,REC=1)(( lat(i,j) ,i=1,nxp),j=1,nyp)
CLOSE(29)
PRINT*," "
PRINT*, "Reading: GMI_lon.dat"
OPEN(UNIT=29,FILE="GMI_lon.dat" ,FORM="UNFORMATTED", &
     ACCESS="DIRECT",STATUS="OLD",RECL=RECLEN)
READ(29,REC=1)(( lon(i,j) ,i=1,nxp),j=1,nyp)
CLOSE(29)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                     !
!                   Map Tbs to rect_tb                                !
!                                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
rect_tb(:,:) = 0.0
print*,"rect_lat_begin = ",rect_lat_begin
print*," rect_lat_end = ",rect_lat_begin + rect_delta_lat*(rect_nyp-1)
print*,"rect_lon_begin = ",rect_lon_begin
print*," rect_lon_end = ",rect_lon_begin + rect_delta_lon*(rect_nxp-1)

do i = 1, nxp
  do j = 1, nyp

   if( rect_lat_begin < lat(i,j) .and.          &
       lat(i,j) < (rect_lat_begin + rect_delta_lat*(rect_nyp-1)) )then
   if( rect_lon_begin < lon(i,j) .and.          &
       lon(i,j) < (rect_lon_begin + rect_delta_lon*(rect_nxp-1)) )then

    rect_i = nint( abs(rect_lon_begin - lon(i,j))/rect_delta_lon ) + 1
    rect_j = nint( abs(rect_lat_begin - lat(i,j))/rect_delta_lat ) + 1
    rect_tb(rect_i,rect_j) = 1.0

   endif
   endif

  enddo
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                     !
!                   Write rect_tb values                              !
!                                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
RECLEN=4*RECT_NXP*RECT_NYP
PRINT*," "
PRINT*, "Writing gmi_orbit_window.dat of size ",RECLEN
PRINT*, "RECT_NXP = ",RECT_NXP," RECT_NYP = ",RECT_NYP
OPEN(UNIT=29,FILE="gmi_orbit_window.dat",FORM="UNFORMATTED",ACCESS="DIRECT",  &
     STATUS="REPLACE",RECL=RECLEN)
WRITE(29,REC=1)(( rect_tb(i,j),i=1,rect_nxp),j=1,rect_nyp)
CLOSE(29)
!
end program remap_gmi_to_fv3
