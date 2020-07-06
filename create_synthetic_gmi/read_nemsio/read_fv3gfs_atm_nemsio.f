!- - - - -- - -- - -- - -- - - -- - --  -- - -- - -- - - -- - - - -- - --
! the program reads and list the contents in a nemsio file
! - - - - -- - -- - -- - -- - - -- - --  -- - -- - -- - - -- - - - -- -
!  Revision history
!   Sep 2008:  Jun Wang Initial code
!   Jan 2019:  Ting-Chi Wu: add output to binary for GrADs plotting
! --
  program main
!
  use nemsio_module
  implicit none
!
  integer, parameter:: double=selected_real_kind(p=13,r=200)
  type(nemsio_gfile) :: gfile,gfilem2,gfilem3,gfiled2
!
  real (kind=8) timef
  character(255) cin
  character(8) gdatatype,modelname
  character(2) level
  real,allocatable  :: tmp(:)
!---------------------------------------------------------------------------
!--- nemsio meta data
  real  isecond,stime,etime,dummy
  integer nrec,im,jm,lm,l,idate(7),version, im2,jm2, nframe, &
          ntrac,irealf,nrec1,version1,nmeta1,nfhour,nfminute,nfsecond, &
          nfsecondn,nfsecondd,nmeta,tlmeta
  integer nsoil,jcap,ncld,idsl,idvc,idvm,idrt,rlon_min,rlon_max, &
          rlat_min,rlat_max
  integer nmetavari,nmetavarr,nmetavarl,nmetavarc,nmetavarr8,    &
          nmetaaryi,nmetaaryr,nmetaaryl,nmetaaryc,nmetaaryr8
  integer ihrst,idat(3),mp_physics,sf_surface_physics,icycle,fieldsize
  logical global, run,extrameta
  character(16),allocatable :: recname(:),reclevtyp(:)
  integer,allocatable       :: reclev(:)
  real,allocatable          :: cpi(:),ri(:)

!################################################################3  
! Added by Ting-Chi Wu 2019/1/10
  real(4), allocatable :: vcoord(:,:,:), lat(:), lon(:)
  integer nvcoord
  integer, parameter:: i_missing=-9999
  integer k
  real, parameter::  zero_001=0.001
  real, allocatable :: ak5(:), bk5(:), ck5(:), tref5(:)
  real, parameter :: cp=1.0046e+3 ! GSI constants.f90 (global)
  real, parameter :: rd=2.8705e+2 ! GSI constants.f90 (global)
  real, parameter :: rv=4.6150e+2 ! GSI constants.f90 (global)
  real rd_over_cp, trk, kapr, kap1, fv

  real, allocatable :: lat2d(:,:), lon2d(:,:)
  real, allocatable :: ugrd(:,:,:), vgrd(:,:,:)
  real, allocatable :: temp(:,:,:), spfh(:,:,:), tv(:,:,:), ps(:,:)
  real, allocatable :: prsi(:,:,:), prsl(:,:,:)
  real, allocatable :: clwmr(:,:,:), rwmr(:,:,:), icmr(:,:,:)
  real, allocatable :: snmr(:,:,:), grle(:,:,:)
  integer iflag_tmp, iflag_spfh, iflag_pres, iflag_ugrd, iflag_vgrd
  integer iflag_clwmr, iflag_rwmr, iflag_icmr, iflag_snmr, iflag_grle
  integer tmp_jrec, spfh_jrec, pres_jrec, ugrd_jrec, vgrd_jrec
  integer clwmr_jrec, rwmr_jrec, icmr_jrec, snmr_jrec, grle_jrec
  integer(4) reclen
  integer nx, ny
!################################################################3  
!---------------------------------------------------------------------------
!--- local vars
  character(16) vname
  character(32) gtype
  character(16) vlevtyp
  integer i,ii,j,jj,jrec,krec,vlev,iret,lev,ista,iend,jsta,jend
!---------------------------------------------------------------------------
!
  character(16),allocatable :: variname(:),varrname(:),varlname(:),varcname(:),varr8name(:), &
                               aryiname(:),aryrname(:),arylname(:),arycname(:),aryr8name(:)
  integer,allocatable :: varival(:),aryilen(:),aryrlen(:),aryllen(:),aryclen(:),aryr8len(:)
  integer,allocatable :: aryival(:,:)
  real,allocatable :: varrval(:),aryrval(:,:)
  real(8),allocatable :: varr8val(:),aryr8val(:,:)
  logical,allocatable :: varlval(:),arylval(:,:)
  character(16),allocatable :: varcval(:),arycval(:,:)
!
!---------------------------------------------------------------------------
!
!-------------set up nemsio write--------------------------
  call nemsio_init(iret=iret)
  print *,'nemsio_init, iret=',iret
!
!+++++++++++++++++ read nemsil file with 2 meta data
!+++++++++++++++++++++++++++
!
!--- open gfile for reading
  print *,'3b:: start reading nemsio file '
!  cin='nemsio_2meta_big'
  call getarg(1,cin)
  call nemsio_open(gfile,trim(cin),'read',iret=iret)
  if(iret/=0) print *,'3b:: after open read, ',trim(cin), ' iret=',iret
!
!--- get dimension
  im=0;jm=0;lm=0;nframe=0;nrec=0
  call nemsio_getfilehead(gfile,dimx=im,dimy=jm,dimz=lm,nframe=nframe,nrec=nrec,&
       gdatatype=gdatatype,modelname=modelname,nmeta=nmeta,ntrac=ntrac,tlmeta=tlmeta,iret=iret)
  print *,'3b:: gfilem2,im=',im,'jm=',jm,'lm=',lm,'nframe=',nframe,'nrec=',nrec, &
       'gdatatype=',gdatatype,' modelname=',modelname,' nmeta=',nmeta,'ntrac=',ntrac, &
       'tlmeta=',tlmeta,'iret=',iret
!--- meta data info
  call nemsio_getfilehead(gfile,nfhour=nfhour,nfminute=nfminute,nsoil=nsoil,ncldt=ncld,&
       idsl=idsl,idvc=idvc,idvm=idvm,idrt=idrt,iret=iret)
  print *,'3b:: gfilem2,nfhour=',nfhour,'jcap=',jcap,'ncld=',ncld,'idvc=',idvc,'idrt=',idrt,'idsl=',idsl
!
! call nemsio_getheadvar(gfile,'nfhour',nfhour,iret=iret)
! print *,'nfhour=',nfhour
! call nemsio_getheadvar(gfile,'latf', latf,iret=iret)
! print *,'latf=',latf
!
  call nemsio_getfilehead(gfile,nmetavari=nmetavari,nmetavarr=nmetavarr,nmetavarl=nmetavarl, &
       nmetavarc=nmetavarc,nmetavarr8=nmetavarr8,nmetaaryi=nmetaaryi,nmetaaryr=nmetaaryr,    &
       nmetaaryr8=nmetaaryr8,nmetaaryl=nmetaaryl, nmetaaryc=nmetaaryc)
  print *,'nmetavari=',nmetavari,'nmetavarr=',nmetavarr,'nmetavarl=',nmetavarl,            &
          'nmetavarc=',nmetavarc,'nmetavarr8=',nmetavarr8
  print *,'nmetaaryi=',nmetaaryi,'nmetaaryr=',nmetaaryr,'nmetaaryl=',nmetaaryl,            &
          'nmetaaryc=',nmetaaryc,'nmetaaryr8=',nmetaaryr8
  if(nmetavari>0) then
    allocate(variname(nmetavari),varival(nmetavari))
    call nemsio_getfilehead(gfile,variname=variname,varival=varival)
    print *,'variname=',variname,'varival=',varival
  endif
  if(nmetavarr>0) then
    allocate(varrname(nmetavarr),varrval(nmetavarr))
    call nemsio_getfilehead(gfile,varrname=varrname,varrval=varrval)
    print *,'varrname=',varrname,'varrval=',varrval
  endif
  if(nmetavarr8>0) then
    allocate(varr8name(nmetavarr8),varr8val(nmetavarr8))
    call nemsio_getfilehead(gfile,varr8name=varr8name,varr8val=varr8val)
    print *,'varr8name=',varr8name,'varr8val=',varr8val
  endif
  if(nmetavarl>0) then
    allocate(varlname(nmetavarl),varlval(nmetavarl))
    call nemsio_getfilehead(gfile,varlname=varlname,varlval=varlval)
    print *,'varlname=',varlname,'varlval=',varlval
  endif
  if(nmetavarc>0) then
    allocate(varcname(nmetavarc),varcval(nmetavarc))
    call nemsio_getfilehead(gfile,varcname=varcname,varcval=varcval)
    print *,'varcname=',varcname,'varcval=',varcval
  endif

  if(nmetaaryi>0) then
    allocate(aryiname(nmetaaryi),aryilen(nmetaaryi))
    call nemsio_getfilehead(gfile,aryiname=aryiname,aryilen=aryilen)
    print *,'aryiname=',aryiname,'aryilen=',aryilen
    allocate(aryival(maxval(aryilen),nmetaaryi))
    call nemsio_getfilehead(gfile,aryival=aryival)
    do i=1,nmetaaryi
      print *,'aryiname=',aryiname(i),aryilen(i),aryival(1:aryilen(i),i)
    enddo
  endif
  if(nmetaaryr>0) then
    allocate(aryrname(nmetaaryr),aryrlen(nmetaaryr))
    call nemsio_getfilehead(gfile,aryrname=aryrname,aryrlen=aryrlen)
    print *,'aryrname=',aryrname,'aryrlen=',aryrlen
    allocate(aryrval(maxval(aryrlen),nmetaaryr))
    call nemsio_getfilehead(gfile,aryrval=aryrval)
    do i=1,nmetaaryr
      print *,'aryrname=',aryrname(i),aryrlen(i),aryrval(1:aryrlen(i),i)
    enddo
  endif
  if(nmetaaryr8>0) then
    allocate(aryr8name(nmetaaryr8),aryr8len(nmetaaryr8))
    call nemsio_getfilehead(gfile,aryr8name=aryr8name,aryr8len=aryr8len)
    print *,'aryr8name=',aryr8name,'aryr8len=',aryr8len
    allocate(aryr8val(maxval(aryr8len),nmetaaryr8))
    call nemsio_getfilehead(gfile,aryr8val=aryr8val)
    do i=1,nmetaaryr8
      print *,'aryr8name=',aryr8name(i),aryr8len(i),aryr8val(1:aryr8len(i),i)
    enddo
  endif
  if(nmetaaryl>0) then
    allocate(arylname(nmetaaryl),aryllen(nmetaaryl))
    call nemsio_getfilehead(gfile,arylname=arylname,aryllen=aryllen)
    print *,'arylname=',arylname,'aryllen=',aryllen
    allocate(arylval(maxval(aryllen),nmetaaryl))
    call nemsio_getfilehead(gfile,arylval=arylval)
    do i=1,nmetaaryl
      print *,'arylname=',arylname(i),aryllen(i),arylval(1:aryllen(i),i)
    enddo
  endif
  if(nmetaaryc>0) then
    allocate(arycname(nmetaaryc),aryclen(nmetaaryc))
    call nemsio_getfilehead(gfile,arycname=arycname,aryclen=aryclen)
    print *,'arycname=',arycname,'aryclen=',aryclen
    allocate(arycval(maxval(aryclen),nmetaaryc))
    call nemsio_getfilehead(gfile,arycval=arycval)
    do i=1,nmetaaryc
      print *,'arycname=',arycname(i),aryclen(i),arycval(1:aryclen(i),i)
    enddo
  endif


!
!---read fields
!
  fieldsize=(im+2*nframe)*(jm+2*nframe)
  allocate(tmp(fieldsize))

!################################################################3  
! Added by Ting-Chi Wu 2019/1/10

!--- lat  
  allocate(lat(fieldsize))
  call nemsio_getfilehead(gfile,lat=lat,iret=iret)
  print *, 'lat=', maxval(lat), minval(lat)
  allocate(lat2d(im+2*nframe,jm+2*nframe))
  lat2d=reshape(lat,(/im+2*nframe,jm+2*nframe/))
!--- lon  
  allocate(lon(fieldsize))
  call nemsio_getfilehead(gfile,lon=lon,iret=iret)
  print *, 'lon=', maxval(lon), minval(lon)
  allocate(lon2d(im+2*nframe,jm+2*nframe))
  lon2d=reshape(lon,(/im+2*nframe,jm+2*nframe/))

!-- write lat/lon to a binary file  
  reclen=4*fieldsize
  print *, 'Writing latitudes of size: ',reclen
  open(unit=29,file='xlat.dat',form='unformatted',access='direct', &
       status='replace',recl=reclen)
!  write(29,rec=1) ((lat2d(i,j),i=1,im+2*nframe),j=1,jm+2*nframe)     
  write(29,rec=1) ((lat2d(i,j),i=1,im+2*nframe),j=jm+2*nframe,1,-1)     
  close(29)
  print*, 'Writing longitudes of size: ',reclen
  open(unit=29,file='xlon.dat',form='unformatted',access='direct', &
       status='replace',recl=reclen)
  write(29,rec=1) ((lon2d(i,j),i=1,im+2*nframe),j=1,jm+2*nframe)     
  close(29)


!-- vcoord  
  allocate(vcoord(lm+1,3,2))
  call nemsio_getfilehead(gfile,vcoord=vcoord,iret=iret)
  if(iret==0) then
    print *,'levs=',lm,'vcoord(1:levs+1,1,1)=',vcoord(:,1,1)
    print *,'levs=',lm,'vcoord(1:levs+1,2,1)=',vcoord(:,2,1)
    print *,'levs=',lm,'vcoord(1:levs+1,3,1)=',vcoord(:,3,1)
    print *,'levs=',lm,'vcoord(1:levs+1,1,2)=',vcoord(:,1,2)
    print *,'levs=',lm,'vcoord(1:levs+1,2,2)=',vcoord(:,2,2)
    print *,'levs=',lm,'vcoord(1:levs+1,3,2)=',vcoord(:,3,2)
  endif

!-- Determine the type of vertical coordinate used by model (GSI: gesinfo.f90)
  nvcoord=3
  if(maxval(vcoord(:,3,1))==0. .and. &
     minval(vcoord(:,3,1))==0. ) then
     nvcoord=2
     if(maxval(vcoord(:,2,1))==0. .and. &
       minval(vcoord(:,2,1))==0. ) then
       nvcoord=1
     endif  
  endif   
  if(idsl==i_missing .or. idsl < 1) then
    idsl=1  
    if (nvcoord==3) idsl=2  
  endif  
  print *,'nvcoord=',nvcoord

  allocate(ak5(lm+1),bk5(lm+1),ck5(lm+1),tref5(lm))
  do k=1,lm+1
    ak5(k)=0.0
    bk5(k)=0.0
    ck5(k)=0.0
  enddo 
  do k=1,lm
     tref5(k)=300.0
  enddo   
  if (nvcoord==1) then
      do k=1,lm+1 
         bk5(k) = vcoord(k,1,1)
      enddo   
  elseif (nvcoord==2) then    
      do k=1,lm+1 
         ak5(k) = vcoord(k,1,1)*zero_001
         bk5(k) = vcoord(k,2,1)
      enddo   
  elseif (nvcoord==3) then    
      do k=1,lm+1 
         ak5(k) = vcoord(k,1,1)*zero_001
         bk5(k) = vcoord(k,2,1)
         ck5(k) = vcoord(k,3,1)
      enddo   
  else    
    print *,'***ERROR*** INVALID value for nvcoord=',nvcoord  
  endif  


!################################################################3  

  iflag_ugrd = 0
  iflag_vgrd = 0
  iflag_tmp = 0
  iflag_spfh = 0
  iflag_pres = 0
  iflag_clwmr = 0
  iflag_rwmr = 0
  iflag_icmr = 0
  iflag_snmr = 0
  iflag_grle = 0
  do jrec=1,nrec
    call nemsio_getrechead(gfile,jrec,vname,vlevtyp,vlev,iret)
    call nemsio_readrec(gfile,jrec,tmp,iret=iret)
!    print *,'3b:: read,jrec=',jrec,'iret=',iret,' vname=',trim(vname), &
!       ' vlevtyp=',trim(vlevtyp),' vlev=',vlev,'data=',maxval(tmp),minval(tmp)
    if (trim(vname)=='ugrd' .and. iflag_ugrd==0) then
       ugrd_jrec=jrec
       iflag_ugrd=1
    endif   
    if (trim(vname)=='vgrd' .and. iflag_vgrd==0) then
       vgrd_jrec=jrec
       iflag_vgrd=1
    endif   
    if (trim(vname)=='tmp' .and. iflag_tmp==0) then
       tmp_jrec=jrec
       iflag_tmp=1
    endif   
    if (trim(vname)=='spfh' .and. iflag_spfh==0) then
       spfh_jrec=jrec
       iflag_spfh=1
    endif   
    if (trim(vname)=='pres' .and. iflag_pres==0) then
       pres_jrec=jrec
       iflag_pres=1
    endif   
    if (trim(vname)=='clwmr' .and. iflag_clwmr==0) then
       clwmr_jrec=jrec
       iflag_clwmr=1
    endif   
    if (trim(vname)=='rwmr' .and. iflag_rwmr==0) then
       rwmr_jrec=jrec
       iflag_rwmr=1
    endif   
    if (trim(vname)=='icmr' .and. iflag_icmr==0) then
       icmr_jrec=jrec
       iflag_icmr=1
    endif   
    if (trim(vname)=='snmr' .and. iflag_snmr==0) then
       snmr_jrec=jrec
       iflag_snmr=1
    endif   
    if (trim(vname)=='grle' .and. iflag_grle==0) then
       grle_jrec=jrec
       iflag_grle=1
    endif   
  enddo

!################################################################3  
! Added by Ting-Chi Wu 2019/1/10

  print *,'jrec for beginning of ugrd=',ugrd_jrec
  print *,'jrec for beginning of vgrd=',vgrd_jrec
  print *,'jrec for beginning of tmp=',tmp_jrec
  print *,'jrec for beginning of spfh=',spfh_jrec
  print *,'jrec for beginning of pres=',pres_jrec
  print *,'jrec for beginning of clwmr=',clwmr_jrec
  print *,'jrec for beginning of rwmr=',rwmr_jrec
  print *,'jrec for beginning of icmr=',icmr_jrec
  print *,'jrec for beginning of snmr=',snmr_jrec
  print *,'jrec for beginning of grle=',grle_jrec

  allocate(temp(lm,im,jm))
  allocate(spfh(lm,im,jm))
  allocate(tv(lm,im,jm))
  allocate(ps(im,jm))
  allocate(prsi(lm+1, im,jm))
  allocate(prsl(lm,im,jm))
  allocate(ugrd(lm,im,jm))
  allocate(vgrd(lm,im,jm))

  allocate(clwmr(lm,im,jm))
  allocate(rwmr(lm,im,jm))
  allocate(icmr(lm,im,jm))
  allocate(snmr(lm,im,jm))
  allocate(grle(lm,im,jm))

  vname='ugrd'
  do jrec=ugrd_jrec,ugrd_jrec+lm-1
!    print *, 'reading u-wind (m/s) jrec =', jrec
    call nemsio_getrechead(gfile,jrec,vname,vlevtyp,vlev,iret)
    call nemsio_readrec(gfile,jrec,tmp,iret=iret)
    ugrd(jrec-ugrd_jrec+1,:,:)=reshape(tmp,(/im+2*nframe,jm+2*nframe/))
  enddo  

  vname='vgrd'
  do jrec=vgrd_jrec,vgrd_jrec+lm-1
!    print *, 'reading v-wind (m/s) jrec =', jrec
    call nemsio_getrechead(gfile,jrec,vname,vlevtyp,vlev,iret)
    call nemsio_readrec(gfile,jrec,tmp,iret=iret)
    vgrd(jrec-vgrd_jrec+1,:,:)=reshape(tmp,(/im+2*nframe,jm+2*nframe/))
  enddo  

  vname='tmp'
  do jrec=tmp_jrec,tmp_jrec+lm-1
!    print *, 'reading temperature (K) jrec =', jrec
    call nemsio_getrechead(gfile,jrec,vname,vlevtyp,vlev,iret)
    call nemsio_readrec(gfile,jrec,tmp,iret=iret)
    temp(jrec-tmp_jrec+1,:,:)=reshape(tmp,(/im+2*nframe,jm+2*nframe/))
  enddo  

  vname='spfh'
  do jrec=spfh_jrec,spfh_jrec+lm-1
!    print *, 'reading specific humidity (kg/kg) jrec =', jrec
    call nemsio_getrechead(gfile,jrec,vname,vlevtyp,vlev,iret)
    call nemsio_readrec(gfile,jrec,tmp,iret=iret)
    spfh(jrec-spfh_jrec+1,:,:)=reshape(tmp,(/im+2*nframe,jm+2*nframe/))
  enddo  

  vname='clwmr'
  do jrec=clwmr_jrec,clwmr_jrec+lm-1
!    print *, 'reading cloud water mixing ratio (kg/kg) jrec =', jrec
    call nemsio_getrechead(gfile,jrec,vname,vlevtyp,vlev,iret)
    call nemsio_readrec(gfile,jrec,tmp,iret=iret)
    clwmr(jrec-clwmr_jrec+1,:,:)=reshape(tmp,(/im+2*nframe,jm+2*nframe/))
  enddo  
  clwmr = clwmr * 1000.0

  vname='rwmr'
  do jrec=rwmr_jrec,rwmr_jrec+lm-1
!    print *, 'reading rain water mixing ratio (kg/kg) jrec =', jrec
    call nemsio_getrechead(gfile,jrec,vname,vlevtyp,vlev,iret)
    call nemsio_readrec(gfile,jrec,tmp,iret=iret)
    rwmr(jrec-rwmr_jrec+1,:,:)=reshape(tmp,(/im+2*nframe,jm+2*nframe/))
  enddo  
  rwmr = rwmr * 1000.0

  vname='icmr'
  do jrec=icmr_jrec,icmr_jrec+lm-1
!    print *, 'reading ice water mixing ratio (kg/kg) jrec =', jrec
    call nemsio_getrechead(gfile,jrec,vname,vlevtyp,vlev,iret)
    call nemsio_readrec(gfile,jrec,tmp,iret=iret)
    icmr(jrec-icmr_jrec+1,:,:)=reshape(tmp,(/im+2*nframe,jm+2*nframe/))
  enddo  
  icmr = icmr * 1000.0

  vname='snmr'
  do jrec=snmr_jrec,snmr_jrec+lm-1
!    print *, 'reading snow mixing ratio (kg/kg) jrec =', jrec
    call nemsio_getrechead(gfile,jrec,vname,vlevtyp,vlev,iret)
    call nemsio_readrec(gfile,jrec,tmp,iret=iret)
    snmr(jrec-snmr_jrec+1,:,:)=reshape(tmp,(/im+2*nframe,jm+2*nframe/))
  enddo  
  snmr = snmr * 1000.0

  vname='grle'
  do jrec=grle_jrec,grle_jrec+lm-1
!    print *, 'reading graupel mixing ratio (kg/kg) jrec =', jrec
    call nemsio_getrechead(gfile,jrec,vname,vlevtyp,vlev,iret)
    call nemsio_readrec(gfile,jrec,tmp,iret=iret)
    grle(jrec-grle_jrec+1,:,:)=reshape(tmp,(/im+2*nframe,jm+2*nframe/))
  enddo  
  grle = grle * 1000.0

!-- compute virtual temperature tv (K)
!-- to be used to construct 3-D pressure field

  fv=rv/rd-1.0
  tv=temp*(1.0+fv*spfh)

  vname='pres'
  jrec=pres_jrec
!  print *, 'reading surface pressure (Pascal) jrec =', jrec
  call nemsio_getrechead(gfile,jrec,vname,vlevtyp,vlev,iret)
  call nemsio_readrec(gfile,jrec,tmp,iret=iret)
  ps(:,:)=reshape(tmp,(/im+2*nframe,jm+2*nframe/))
  ps=0.01*ps ! convert Pascal to Hecto-Pascal


!-- construct 3-d pressure field at (GSI guess_grid.f90: global)

  rd_over_cp = rd/cp
  kapr = 1.0/rd_over_cp
  kap1 = rd_over_cp+1.0

  print *, 'idvc=',idvc
  do k=1,lm+1
    print *, 'ak5(k), bk5(k) = ', k, ak5(k), bk5(k)
  enddo  

  ps=0.1*ps ! convert Hecto-Pascal to Centi-Pascal
  do k=1,lm+1
    do j=1,jm
      do i=1,im   
      
        ! idvc (idvc5 in GSI gridmod.f90) is vertical coordinate identifier
        ! idsl (idsl5 in GSI gridmod.f90) is midlayer pressure definition (1: Philips; 2: average)
        if (idvc==1 .or. idvc==2) then
          prsi(k,i,j) = ak5(k)+(bk5(k)*ps(i,j)) 
        else if (idvc==3) then  
          print *, 'hi!'  
          if (k==1) then  
            prsi(k,i,j)=ps(i,j)  
          else if (k==lm+1) then 
            prsi(k,i,j)=0.0  
          else  
            trk=(0.5*(tv(k-1,i,j)+tv(k,i,j))/tref5(k))**kapr   
            prsi(k,i,j)=ak5(k)+(bk5(k)*ps(i,j))+(ck5(k)*trk)    
          endif
        endif  
        prsi(k,i,j)=max(prsi(k,i,j),0.0)
      enddo
    enddo
  enddo

  print *, 'idsl=',idsl
  if (idsl/=2) then
    do k=1,lm
      do j=1,jm
        do i=1,im   
          prsl(k,i,j)=((prsi(k,i,j)**kap1-prsi(k+1,i,j)**kap1)/ &
                      (kap1*(prsi(k,i,j)-prsi(k+1,i,j))))**kapr
        enddo
      enddo
    enddo
  else  
    do k=1,lm
      do j=1,jm
        do i=1,im   
          prsl(k,i,j)=0.5*(prsi(k,i,j)+prsi(k+1,i,j))
        enddo
      enddo
    enddo
  endif

!  stop

  ! convert ps, prsi, prsl from Centi-Pascal to Hecto-Pascal
  prsi=10.*prsi
  prsl=10.*prsl
  ps=10.*ps

!-- write field to a binary file  

  print*, 'Writing U-Wind (m/s)'
  open(unit=29,file='U000000',form='unformatted',access='direct', &
       status='replace',recl=reclen)
  do k=1,lm     
    write(29,rec=k) ((ugrd(k,i,j),i=1,im+2*nframe),j=jm+2*nframe,1,-1)     
  enddo  
  close(29)

  print*, 'Writing V-Wind (m/s)'
  open(unit=29,file='V000000',form='unformatted',access='direct', &
       status='replace',recl=reclen)
  do k=1,lm     
    write(29,rec=k) ((vgrd(k,i,j),i=1,im+2*nframe),j=jm+2*nframe,1,-1)     
  enddo  
  close(29)

  print*, 'Writing Temperature (K)'
  open(unit=29,file='TEMPK000000',form='unformatted',access='direct', &
       status='replace',recl=reclen)
  do k=1,lm     
    write(29,rec=k) ((temp(k,i,j),i=1,im+2*nframe),j=jm+2*nframe,1,-1)     
  enddo  
  close(29)

  print*, 'Writing Specific Humidity (kg/kg)'
  open(unit=29,file='SPH000000',form='unformatted',access='direct', &
       status='replace',recl=reclen)
  do k=1,lm     
    write(29,rec=k) ((spfh(k,i,j),i=1,im+2*nframe),j=jm+2*nframe,1,-1)     
  enddo  
  close(29)

  print*, 'Writing Cloud Water Mass Mixing Ratio (g/kg)'
  open(unit=29,file='CLWMR000000',form='unformatted',access='direct', &
       status='replace',recl=reclen)
  do k=1,lm     
    write(29,rec=k) ((clwmr(k,i,j),i=1,im+2*nframe),j=jm+2*nframe,1,-1)     
  enddo  
  close(29)
  
  print*, 'Writing Rain Water Mass Mixing Ratio (g/kg)'
  open(unit=29,file='RWMR000000',form='unformatted',access='direct', &
       status='replace',recl=reclen)
  do k=1,lm     
    write(29,rec=k) ((rwmr(k,i,j),i=1,im+2*nframe),j=jm+2*nframe,1,-1)     
  enddo  
  close(29)

  print*, 'Writing Ice Water Mass Mixing Ratio (g/kg)'
  open(unit=29,file='ICMR000000',form='unformatted',access='direct', &
       status='replace',recl=reclen)
  do k=1,lm     
    write(29,rec=k) ((icmr(k,i,j),i=1,im+2*nframe),j=jm+2*nframe,1,-1)     
  enddo  
  close(29)
  
  print*, 'Writing Snow Water Mass Mixing Ratio (g/kg)'
  open(unit=29,file='SNMR000000',form='unformatted',access='direct', &
       status='replace',recl=reclen)
  do k=1,lm     
    write(29,rec=k) ((snmr(k,i,j),i=1,im+2*nframe),j=jm+2*nframe,1,-1)     
  enddo  
  close(29)

  print*, 'Writing Graupel Mass Mixing Ratio (g/kg)'
  open(unit=29,file='GRLE000000',form='unformatted',access='direct', &
       status='replace',recl=reclen)
  do k=1,lm     
    write(29,rec=k) ((grle(k,i,j),i=1,im+2*nframe),j=jm+2*nframe,1,-1)     
  enddo  
  close(29)
!  print*, 'Writing Surface Pressure (hPa)'
!  open(unit=29,file='PSFC000000',form='unformatted',access='direct', &
!       status='replace',recl=reclen)
!!  write(29,rec=1) ((ps(i,j),i=1,im+2*nframe),j=1,jm+2*nframe)     
!  write(29,rec=1) ((ps(i,j),i=1,im+2*nframe),j=jm+2*nframe,1,-1)     
!  close(29)

!  print*, 'Writing 3-D Pressure (prsi) (hPa (mb))'
!  open(unit=29,file='PRESSI000000',form='unformatted',access='direct', &
!       status='replace',recl=reclen)
!  do k=1,lm+1     
!    write(29,rec=k) ((prsi(i,j,k),i=1,im+2*nframe),j=jm+2*nframe,1,-1)     
!  enddo  
!  close(29)

  print*, 'Writing 3-D Pressure (hPa (mb))'
  open(unit=29,file='PRESS000000',form='unformatted',access='direct', &
       status='replace',recl=reclen)
  do k=1,lm     
    write(29,rec=k) ((prsl(k,i,j),i=1,im+2*nframe),j=jm+2*nframe,1,-1)     
  enddo  
  close(29)

!--- close nemsio file
  call nemsio_close(gfile,iret=iret)
  if ( iret .ne.0) print *,'iret=',iret

!!---------------------------------------------------------------------------
  deallocate(tmp)
  deallocate(lat)
  deallocate(lat2d)
  deallocate(lon)
  deallocate(lon2d)
  deallocate(vcoord)
  deallocate(ak5)
  deallocate(bk5)
  deallocate(ck5)
  deallocate(tref5)
  deallocate(ugrd)
  deallocate(vgrd)
  deallocate(temp)
  deallocate(spfh)
  deallocate(tv)
  deallocate(ps)
  deallocate(prsi)
  deallocate(prsl)
  deallocate(clwmr)
  deallocate(rwmr)
  deallocate(icmr)
  deallocate(snmr)
  deallocate(grle)
!
!---------------------------------------------------------------------------
!
  call nemsio_finalize()
! - - - - -- - -- - -- - -- - - -- - --  -- - -- - -- - - -- - - - -- -
! --
  stop

 end program

