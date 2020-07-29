program postprot

USE netcdf

IMPLICIT NONE

INTEGER :: fidT, status, dimID_time_counter, dimID_x, dimID_y, dimID_axis_nbounds, dimID_deptht, &
&          mtime_counter, mx, my, maxis_nbounds, mdeptht, toce_ID, ssh_ID, soce_ID, fwfisf_ID,   &
&          sbt_ID, sbs_ID, deptht_bounds_ID, deptht_ID, fidM, fidSBC, cnt, nx_CRS, ny_CRS,       &
&          isfgammas_ID, isfgammat_ID, isfhalindr_ID, isfthermdr_ID, i, j, k, l, ipr, jpr,       &
&          toce_e3t_ID, soce_e3t_ID, e3t_ID

INTEGER :: salinityYZ_ID, temperatureYZ_ID, salinityXZ_ID, temperatureXZ_ID, bottomSalinity_ID,  &
&          bottomTemperature_ID, halineDriving_ID, thermalDriving_ID, frictionVelocity_ID,       &
&          meltRate_ID, bathymetry_ID, iceDraft_ID, meanSalinity_ID, meanTemperature_ID,         &
&          totalOceanVolume_ID, totalMeltFlux_ID, meanMeltRate_ID

INTEGER :: fidMSH, dimID_t, dimID_z, mt, mz, tmask_ID, misf_ID, mbathy_ID, x_ID, y_ID, z_ID,     &
&          e3t_1d_ID, gdepw_1d_ID, gdept_1d_ID, e3t_0_ID, ff_ID, e2t_ID, knew, mznew, ksrf,      &
&          e1t_ID, time_counter_ID, gdepw_0_ID, gdept_0_ID, isfdraft_ID

CHARACTER(LEN=120) :: file_in_T, file_in_SBC, file_MSH, file_out

INTEGER*1,ALLOCATABLE,DIMENSION(:,:,:,:) :: tmask

INTEGER*1,ALLOCATABLE,DIMENSION(:,:) :: maskisf

INTEGER*2,ALLOCATABLE,DIMENSION(:,:,:) :: misf, mbathy

REAL*4,ALLOCATABLE,DIMENSION(:) :: deptht, ISOMIP_OceVol, time

REAL*8,ALLOCATABLE,DIMENSION(:) :: znew, ISOMIP_Tmean, ISOMIP_Smean, ISOMIP_mean_melt, ISOMIP_total_melt

REAL*4,ALLOCATABLE,DIMENSION(:,:) :: ISOMIP_isf_draft, ISOMIP_bathy, CRS_ISOMIP_isf_draft, CRS_ISOMIP_bathy

REAL*8,ALLOCATABLE,DIMENSION(:,:) :: e3t_1d, gdepw_1d, gdept_1d

REAL*4,ALLOCATABLE,DIMENSION(:,:,:) :: isfdraft, e1t, e2t, ssh, sbt, sbs,                       &
&                                      isfgammas, isfgammat, isfhalindr, isfthermdr, fwfisf,    &
&                                      ISOMIP_melt_rate, ISOMIP_ustar, ISOMIP_therm_driv,       &
&                                      ISOMIP_halin_driv, ISOMIP_Tbot, ISOMIP_Sbot,             &
&                                      CRS_ISOMIP_melt_rate, CRS_ISOMIP_ustar, CRS_ISOMIP_therm_driv,&
&                                      CRS_ISOMIP_halin_driv, CRS_ISOMIP_Tbot, CRS_ISOMIP_Sbot

REAL*8,ALLOCATABLE,DIMENSION(:,:,:) :: ISOMIP_TXZ, ISOMIP_SXZ, ISOMIP_TYZ, ISOMIP_SYZ,               &
&                                      CRS_ISOMIP_TXZ, CRS_ISOMIP_SXZ, CRS_ISOMIP_TYZ, CRS_ISOMIP_SYZ

REAL*4,ALLOCATABLE,DIMENSION(:,:,:,:) :: toce, soce, gdepw_0, gdept_0, toce_e3t, soce_e3t, e3t

REAL*8,ALLOCATABLE,DIMENSION(:,:,:,:) :: e3t_0

INTEGER*4, DIMENSION(12) :: daym1

INTEGER*4 :: an

REAL*4 :: Gt

REAL*8 :: isf_area, tmp_XZ, tmp_YZ, eps, deltaz

!---------------------------------------

file_in_T   = 'MONTHLY/ISOMIP1km-CCCC_monthly_YYYY_grid_T.nc'
file_in_SBC = 'MONTHLY/ISOMIP1km-CCCC_monthly_YYYY_SBC.nc'
file_MSH    = '../../input/nemo_ISOMIP1km/mesh_mask_ISOMIP1km_EEEE.nc'
file_out    = 'MONTHLY/ISOMIP1km-CCCC_monthly_YYYY_tmpT.nc'

Gt = GGGG  !! GammaT in NEMO's namelist

eps=1.d-9

an=YYYY

!---------------------------------------
! Read grid_T file

status = NF90_OPEN(TRIM(file_in_T),0,fidT)
call erreur(status,.TRUE.,"read grid_T")

status = NF90_INQ_DIMID(fidT,"time_counter",dimID_time_counter)
call erreur(status,.TRUE.,"inq_dimID_time_counter")
status = NF90_INQ_DIMID(fidT,"x",dimID_x)
call erreur(status,.TRUE.,"inq_dimID_x")
status = NF90_INQ_DIMID(fidT,"y",dimID_y)
call erreur(status,.TRUE.,"inq_dimID_y")
status = NF90_INQ_DIMID(fidT,"deptht",dimID_deptht)
call erreur(status,.TRUE.,"inq_dimID_deptht")

status = NF90_INQUIRE_DIMENSION(fidT,dimID_time_counter,len=mtime_counter)
call erreur(status,.TRUE.,"inq_dim_time_counter")
status = NF90_INQUIRE_DIMENSION(fidT,dimID_x,len=mx)
call erreur(status,.TRUE.,"inq_dim_x")
status = NF90_INQUIRE_DIMENSION(fidT,dimID_y,len=my)
call erreur(status,.TRUE.,"inq_dim_y")
status = NF90_INQUIRE_DIMENSION(fidT,dimID_deptht,len=mdeptht)
call erreur(status,.TRUE.,"inq_dim_deptht")

ALLOCATE(  time(mtime_counter)  )
ALLOCATE(  toce(mx,my,mdeptht,mtime_counter)  )
ALLOCATE(  toce_e3t(mx,my,mdeptht,mtime_counter)  )
ALLOCATE(  e3t(mx,my,mdeptht,mtime_counter)  )
ALLOCATE(  ssh(mx,my,mtime_counter)  )
ALLOCATE(  soce(mx,my,mdeptht,mtime_counter)  )
ALLOCATE(  soce_e3t(mx,my,mdeptht,mtime_counter)  )
ALLOCATE(  sbt(mx,my,mtime_counter)  )
ALLOCATE(  sbs(mx,my,mtime_counter)  )

status = NF90_INQ_VARID(fidT,"toce",toce_ID)
call erreur(status,.TRUE.,"inq_toce_ID")
status = NF90_INQ_VARID(fidT,"toce_e3t",toce_e3t_ID)
call erreur(status,.TRUE.,"inq_toce_e3t_ID")
status = NF90_INQ_VARID(fidT,"e3t",e3t_ID)
call erreur(status,.TRUE.,"inq_e3t_ID")
status = NF90_INQ_VARID(fidT,"ssh",ssh_ID)
call erreur(status,.TRUE.,"inq_ssh_ID")
status = NF90_INQ_VARID(fidT,"soce",soce_ID)
call erreur(status,.TRUE.,"inq_soce_ID")
status = NF90_INQ_VARID(fidT,"soce_e3t",soce_e3t_ID)
call erreur(status,.TRUE.,"inq_soce_e3t_ID")
status = NF90_INQ_VARID(fidT,"sbt",sbt_ID)
call erreur(status,.TRUE.,"inq_sbt_ID")
status = NF90_INQ_VARID(fidT,"sbs",sbs_ID)
call erreur(status,.TRUE.,"inq_sbs_ID")

status = NF90_GET_VAR(fidT,toce_ID,toce)
call erreur(status,.TRUE.,"getvar_toce")
status = NF90_GET_VAR(fidT,toce_e3t_ID,toce_e3t)
call erreur(status,.TRUE.,"getvar_toce_e3t")
status = NF90_GET_VAR(fidT,e3t_ID,e3t)
call erreur(status,.TRUE.,"getvar_e3t")
status = NF90_GET_VAR(fidT,ssh_ID,ssh)
call erreur(status,.TRUE.,"getvar_ssh")
status = NF90_GET_VAR(fidT,soce_ID,soce)
call erreur(status,.TRUE.,"getvar_soce_e3t")
status = NF90_GET_VAR(fidT,soce_e3t_ID,soce_e3t)
call erreur(status,.TRUE.,"getvar_soce")
status = NF90_GET_VAR(fidT,sbt_ID,sbt)
call erreur(status,.TRUE.,"getvar_sbt")
status = NF90_GET_VAR(fidT,sbs_ID,sbs)
call erreur(status,.TRUE.,"getvar_sbs")

status = NF90_CLOSE(fidT)
call erreur(status,.TRUE.,"fin_lecture")

!---------------------------------------
! Read SBC file 

write(*,*) 'Reading ', TRIM(file_in_SBC)

status = NF90_OPEN(TRIM(file_in_SBC),0,fidSBC)
call erreur(status,.TRUE.,"read SBC")

ALLOCATE(  isfgammas(mx,my,mtime_counter)  )
ALLOCATE(  isfgammat(mx,my,mtime_counter)  )
ALLOCATE(  isfhalindr(mx,my,mtime_counter)  )
ALLOCATE(  isfthermdr(mx,my,mtime_counter)  )
ALLOCATE(  fwfisf(mx,my,mtime_counter)  )

status = NF90_INQ_VARID(fidSBC,"isfgammas",isfgammas_ID)
call erreur(status,.TRUE.,"inq_isfgammas_ID")
status = NF90_INQ_VARID(fidSBC,"isfgammat",isfgammat_ID)
call erreur(status,.TRUE.,"inq_isfgammat_ID")
status = NF90_INQ_VARID(fidSBC,"isfhalindr",isfhalindr_ID)
call erreur(status,.TRUE.,"inq_isfhalindr_ID")
status = NF90_INQ_VARID(fidSBC,"isfthermdr",isfthermdr_ID)
call erreur(status,.TRUE.,"inq_isfthermdr_ID")
status = NF90_INQ_VARID(fidSBC,"fwfisf",fwfisf_ID)
call erreur(status,.TRUE.,"inq_fwfisf_ID")

status = NF90_GET_VAR(fidSBC,isfgammas_ID,isfgammas)
call erreur(status,.TRUE.,"getvar_isfgammas")
status = NF90_GET_VAR(fidSBC,isfgammat_ID,isfgammat)
call erreur(status,.TRUE.,"getvar_isfgammat")
status = NF90_GET_VAR(fidSBC,isfhalindr_ID,isfhalindr)
call erreur(status,.TRUE.,"getvar_isfhalindr")
status = NF90_GET_VAR(fidSBC,isfthermdr_ID,isfthermdr)
call erreur(status,.TRUE.,"getvar_isfthermdr")
status = NF90_GET_VAR(fidSBC,fwfisf_ID,fwfisf)
call erreur(status,.TRUE.,"getvar_fwfisf")

status = NF90_CLOSE(fidSBC)
call erreur(status,.TRUE.,"fin_lecture")

!---------------------------------------                   
! Read mesh_mask file :

write(*,*) 'Reading ', TRIM(file_MSH)

status = NF90_OPEN(TRIM(file_MSH),0,fidMSH)          
call erreur(status,.TRUE.,"read mesh_mask file") 

status = NF90_INQ_DIMID(fidMSH,"t",dimID_t)
call erreur(status,.TRUE.,"inq_dimID_t")
status = NF90_INQ_DIMID(fidMSH,"z",dimID_z)
call erreur(status,.TRUE.,"inq_dimID_z")
status = NF90_INQ_DIMID(fidMSH,"y",dimID_y)
call erreur(status,.TRUE.,"inq_dimID_y")
status = NF90_INQ_DIMID(fidMSH,"x",dimID_x)
call erreur(status,.TRUE.,"inq_dimID_x")

status = NF90_INQUIRE_DIMENSION(fidMSH,dimID_t,len=mt)
call erreur(status,.TRUE.,"inq_dim_t")
status = NF90_INQUIRE_DIMENSION(fidMSH,dimID_z,len=mz)
call erreur(status,.TRUE.,"inq_dim_z")

ALLOCATE(  tmask(mx,my,mz,mt)  ) 
ALLOCATE(  misf(mx,my,mt)  ) 
ALLOCATE(  mbathy(mx,my,mt)  ) 
ALLOCATE(  e3t_1d(mz,mt)  ) 
ALLOCATE(  gdepw_1d(mz,mt)  ) 
ALLOCATE(  gdept_1d(mz,mt)  ) 
ALLOCATE(  e3t_0(mx,my,mz,mt)  ) 
ALLOCATE(  e2t(mx,my,mt)  ) 
ALLOCATE(  e1t(mx,my,mt)  ) 
ALLOCATE(  gdepw_0(mx,my,mz,mt)  ) 
ALLOCATE(  gdept_0(mx,my,mz,mt)  ) 
ALLOCATE(  isfdraft(mx,my,mt)  ) 

status = NF90_INQ_VARID(fidMSH,"tmask",tmask_ID)
call erreur(status,.TRUE.,"inq_tmask_ID")
status = NF90_INQ_VARID(fidMSH,"misf",misf_ID)
call erreur(status,.TRUE.,"inq_misf_ID")
status = NF90_INQ_VARID(fidMSH,"mbathy",mbathy_ID)
call erreur(status,.TRUE.,"inq_mbathy_ID")
status = NF90_INQ_VARID(fidMSH,"e3t_1d",e3t_1d_ID)
call erreur(status,.TRUE.,"inq_e3t_1d_ID")
status = NF90_INQ_VARID(fidMSH,"gdepw_1d",gdepw_1d_ID)
call erreur(status,.TRUE.,"inq_gdepw_1d_ID")
status = NF90_INQ_VARID(fidMSH,"gdept_1d",gdept_1d_ID)
call erreur(status,.TRUE.,"inq_gdept_1d_ID")
status = NF90_INQ_VARID(fidMSH,"e3t_0",e3t_0_ID)
call erreur(status,.TRUE.,"inq_e3t_0_ID")
status = NF90_INQ_VARID(fidMSH,"e2t",e2t_ID)
call erreur(status,.TRUE.,"inq_e2t_ID")
status = NF90_INQ_VARID(fidMSH,"e1t",e1t_ID)
call erreur(status,.TRUE.,"inq_e1t_ID")
status = NF90_INQ_VARID(fidMSH,"gdepw_0",gdepw_0_ID)
call erreur(status,.TRUE.,"inq_gdepw_0_ID")
status = NF90_INQ_VARID(fidMSH,"gdept_0",gdept_0_ID)
call erreur(status,.TRUE.,"inq_gdept_0_ID")
status = NF90_INQ_VARID(fidMSH,"isfdraft",isfdraft_ID)
call erreur(status,.TRUE.,"inq_isfdraft_ID")

status = NF90_GET_VAR(fidMSH,tmask_ID,tmask)
call erreur(status,.TRUE.,"getvar_tmask")
status = NF90_GET_VAR(fidMSH,misf_ID,misf)
call erreur(status,.TRUE.,"getvar_misf")
status = NF90_GET_VAR(fidMSH,mbathy_ID,mbathy)
call erreur(status,.TRUE.,"getvar_mbathy")
status = NF90_GET_VAR(fidMSH,e3t_1d_ID,e3t_1d)
call erreur(status,.TRUE.,"getvar_e3t_1d")
status = NF90_GET_VAR(fidMSH,gdepw_1d_ID,gdepw_1d)
call erreur(status,.TRUE.,"getvar_gdepw_1d")
status = NF90_GET_VAR(fidMSH,gdept_1d_ID,gdept_1d)
call erreur(status,.TRUE.,"getvar_gdept_1d")
status = NF90_GET_VAR(fidMSH,e3t_0_ID,e3t_0)
call erreur(status,.TRUE.,"getvar_e3t_0")
status = NF90_GET_VAR(fidMSH,e2t_ID,e2t)
call erreur(status,.TRUE.,"getvar_e2t")
status = NF90_GET_VAR(fidMSH,e1t_ID,e1t)
call erreur(status,.TRUE.,"getvar_e1t")
status = NF90_GET_VAR(fidMSH,gdepw_0_ID,gdepw_0)
call erreur(status,.TRUE.,"getvar_gdepw_0")
status = NF90_GET_VAR(fidMSH,gdept_0_ID,gdept_0)
call erreur(status,.TRUE.,"getvar_gdept_0")
status = NF90_GET_VAR(fidMSH,isfdraft_ID,isfdraft)
call erreur(status,.TRUE.,"getvar_isfdraft")

status = NF90_CLOSE(fidMSH)                      
call erreur(status,.TRUE.,"fin_lecture")     

!------------------------------------------------------
! refined output vertical grid :
mznew=144
deltaz=5.d0  !! vertical resolution (m)
ALLOCATE( znew(mznew) )
write(*,*) 'New vertical grid for the outputs :'
do knew=1,mznew
  znew(knew) = 2.5d0 - knew * deltaz
  if ( - znew(knew) .le. gdept_1d(1,1) ) ksrf=knew
enddo

!----------------------------------------------------------------
! Calculate and interpolate variables on ISOMIP+ std output grid

daym1 = (/ 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 /)
! time at the beginning of the period used for monthly average
do l=1,mtime_counter
  time(l) = ( an * 365 + daym1(l) ) * 86400.0
enddo

ALLOCATE( maskisf(mx,my) )
where( misf(:,:,1) .ge. 2 )
  maskisf(:,:) = 1
elsewhere
  maskisf(:,:) = 0
endwhere

ALLOCATE( ISOMIP_isf_draft (mx-2,my-2) )
ALLOCATE( ISOMIP_bathy     (mx-2,my-2) )

do i=2,mx-1
do j=2,my-1
  if ( misf(i,j,1) .ne. 0 ) then
    ISOMIP_isf_draft (i-1,j-1) = gdepw_0(i,j,misf(i,j,1)    ,1)
    ISOMIP_bathy     (i-1,j-1) = gdepw_0(i,j,mbathy(i,j,1)+1,1)
  else
    ISOMIP_isf_draft (i-1,j-1) = NF90_FILL_FLOAT
    ISOMIP_bathy     (i-1,j-1) = NF90_FILL_FLOAT
  endif
enddo
enddo

!===== 1D =====

ALLOCATE( ISOMIP_OceVol (mtime_counter) )
ALLOCATE( ISOMIP_Tmean  (mtime_counter) )
ALLOCATE( ISOMIP_Smean  (mtime_counter) )
ALLOCATE( ISOMIP_mean_melt (mtime_counter) )
ALLOCATE( ISOMIP_total_melt(mtime_counter) )

ISOMIP_OceVol (:) = 0.e0
ISOMIP_Tmean  (:) = 0.e0
ISOMIP_Smean  (:) = 0.e0 
ISOMIP_mean_melt (:) = 0.e0
ISOMIP_total_melt(:) = 0.e0

DO l=1,mtime_counter
  isf_area = 0.e0
  do i=1,mx
  do j=1,my
    do k=1,mz
      ISOMIP_OceVol (l) = ISOMIP_OceVol (l) +      e3t(i,j,k,l) * e1t(i,j,1) * e2t(i,j,1) * tmask(i,j,k,1)
      ISOMIP_Tmean  (l) = ISOMIP_Tmean  (l) + toce_e3t(i,j,k,l) * e1t(i,j,1) * e2t(i,j,1) * tmask(i,j,k,1)
      ISOMIP_Smean  (l) = ISOMIP_Smean  (l) + soce_e3t(i,j,k,l) * e1t(i,j,1) * e2t(i,j,1) * tmask(i,j,k,1)
    enddo
    isf_area = isf_area + e1t(i,j,1) * e2t(i,j,1) * maskisf(i,j)
    ISOMIP_total_melt(l) = ISOMIP_total_melt(l) - fwfisf(i,j,l) * e1t(i,j,1) * e2t(i,j,1) * maskisf(i,j)   !!  kg/s
  enddo
  enddo
  !-
  ISOMIP_Tmean(l) = ISOMIP_Tmean(l) / ISOMIP_OceVol(l)
  ISOMIP_Smean(l) = ISOMIP_Smean(l) / ISOMIP_OceVol(l)
  !-
  ISOMIP_mean_melt(l) = ISOMIP_total_melt(l) * 1.e-3 / isf_area   !! m.w.e / s
ENDDO

DEALLOCATE( e3t_0, maskisf )

!===== 2D =====

ALLOCATE( ISOMIP_melt_rate  (mx-2,my-2,mtime_counter) )
ALLOCATE( ISOMIP_ustar      (mx-2,my-2,mtime_counter) )
ALLOCATE( ISOMIP_therm_driv (mx-2,my-2,mtime_counter) )
ALLOCATE( ISOMIP_halin_driv (mx-2,my-2,mtime_counter) )
ALLOCATE( ISOMIP_Tbot       (mx-2,my-2,mtime_counter) )
ALLOCATE( ISOMIP_Sbot       (mx-2,my-2,mtime_counter) )

DO l=1,mtime_counter

  do i=2,mx-1
  do j=2,my-1
    !-
    if ( misf(i,j,1) .ge. 2 ) then
      ISOMIP_melt_rate (i-1,j-1,l) = - fwfisf(i,j,l) * 1.e-3   !!  m/s
      ISOMIP_ustar     (i-1,j-1,l) = isfgammat(i,j,l) / Gt     !!  m/s
      ISOMIP_therm_driv(i-1,j-1,l) = isfthermdr(i,j,l)         !!  degC
      ISOMIP_halin_driv(i-1,j-1,l) = - isfhalindr(i,j,l)       !!  psu
    else
      ISOMIP_melt_rate (i-1,j-1,l) = NF90_FILL_FLOAT
      ISOMIP_ustar     (i-1,j-1,l) = NF90_FILL_FLOAT
      ISOMIP_therm_driv(i-1,j-1,l) = NF90_FILL_FLOAT
      ISOMIP_halin_driv(i-1,j-1,l) = NF90_FILL_FLOAT
    endif
    !--
    if ( sum(tmask(i,j,:,1)) .gt. 0 ) then
      ISOMIP_Tbot      (i-1,j-1,l) = sbt(i,j,l)                !!  degC
      ISOMIP_Sbot      (i-1,j-1,l) = sbs(i,j,l)                !!  psu
    else
      ISOMIP_Tbot      (i-1,j-1,l) = NF90_FILL_FLOAT
      ISOMIP_Sbot      (i-1,j-1,l) = NF90_FILL_FLOAT
    endif
    !-
  enddo
  enddo

ENDDO


!===== SECTIONS =====

!if ( 0.5*(gphit(1,21,1)+gphit(1,22,1)) .ne. 40.0 ) then
!  write(*,*) '~!@#$%^* ADAPT SCRIPT TO HAVE SECTION AT y = 40km    >>>>> STOP !!'
!  stop
!endif
!
!if ( 0.5*(glamt(101,1,1)+glamt(102,1,1)) .ne. 520.0 ) then
!  write(*,*) '~!@#$%^* ADAPT SCRIPT TO HAVE SECTION AT x = 520km   >>>>> STOP !!'
!  stop
!endif
write(*,*) ' '
write(*,*) 'WARNING : ADAPT SCRIPT IF NOT COM ISOMIP+ GEOMETRY !!!!!!!!!!!!!!!'
write(*,*) ' '

ALLOCATE( ISOMIP_TXZ(mx-2,mznew,mtime_counter) )
ALLOCATE( ISOMIP_TYZ(my-2,mznew,mtime_counter) )
ALLOCATE( ISOMIP_SXZ(mx-2,mznew,mtime_counter) )
ALLOCATE( ISOMIP_SYZ(my-2,mznew,mtime_counter) )

ISOMIP_TXZ(:,:,:) = 0.d0
ISOMIP_TYZ(:,:,:) = 0.d0
ISOMIP_SXZ(:,:,:) = 0.d0
ISOMIP_SYZ(:,:,:) = 0.d0

DO l=1,mtime_counter

  !-- XZ section at y = 40 km
  do i=2,mx-1
    !-
    do knew=1,ksrf
      ISOMIP_TXZ(i-1,knew,l) = ( tmask(i,21,1,1) * toce(i,21,1,l) + tmask(i,22,1,1) * toce(i,22,1,l) ) / ( tmask(i,21,1,1) + tmask(i,22,1,1) + eps ) &
      &                        + (1 - MIN(1,tmask(i,21,1,1)+tmask(i,22,1,1))) * NF90_FILL_FLOAT
      ISOMIP_SXZ(i-1,knew,l) = ( tmask(i,21,1,1) * soce(i,21,1,l) + tmask(i,22,1,1) * soce(i,22,1,l) ) / ( tmask(i,21,1,1) + tmask(i,22,1,1) + eps ) &
      &                        + (1 - MIN(1,tmask(i,21,1,1)+tmask(i,22,1,1))) * NF90_FILL_FLOAT
    enddo
    !-
    do knew=ksrf+1,mznew
      tmp_XZ = 0.d0
      do jpr=41,42
        k=1
        !- target level fully embedded in an original level :
        if ( -znew(knew)-0.5*deltaz .ge. gdepw_0(i,jpr,k,1) .and. -znew(knew)+0.5*deltaz .le. gdepw_0(i,jpr,k+1,1) ) then
          ISOMIP_TXZ(i-1,knew,l) = ISOMIP_TXZ(i-1,knew,l) + tmask(i,jpr,k,1) * toce(i,jpr,k,l) * deltaz
          ISOMIP_SXZ(i-1,knew,l) = ISOMIP_SXZ(i-1,knew,l) + tmask(i,jpr,k,1) * soce(i,jpr,k,l) * deltaz
          tmp_XZ                 = tmp_XZ                 + tmask(i,jpr,k,1)                   * deltaz
        !- target level over 2 original levels :
        !- NB: here a target level can be over a maximum of three original level (if thinner target levels, adapt script)
        elseif ( -znew(knew)-0.5*deltaz .le. gdepw_0(i,jpr,k+1,1) .and. -znew(knew)+0.5*deltaz .ge. gdepw_0(i,jpr,k+1,1) .and. -znew(knew)-0.5*deltaz .gt. gdepw_0(i,jpr,k,1) ) then
          ISOMIP_TXZ(i-1,knew,l) = ISOMIP_TXZ(i-1,knew,l) + tmask(i,jpr,k  ,1) * toce(i,jpr,k  ,l) * ( gdepw_0(i,jpr,k+1,1) - (-znew(knew)-0.5*deltaz) ) &
          &                                               + tmask(i,jpr,k+1,1) * toce(i,jpr,k+1,l) * (  -znew(knew)+0.5*deltaz  - gdepw_0(i,jpr,k+1,1) )
          ISOMIP_SXZ(i-1,knew,l) = ISOMIP_SXZ(i-1,knew,l) + tmask(i,jpr,k  ,1) * soce(i,jpr,k  ,l) * ( gdepw_0(i,jpr,k+1,1) - (-znew(knew)-0.5*deltaz) ) &
          &                                               + tmask(i,jpr,k+1,1) * soce(i,jpr,k+1,l) * (  -znew(knew)+0.5*deltaz  - gdepw_0(i,jpr,k+1,1) )
          tmp_XZ                 = tmp_XZ                 + tmask(i,jpr,k  ,1)                     * ( gdepw_0(i,jpr,k+1,1) - (-znew(knew)-0.5*deltaz) ) &
          &                                               + tmask(i,jpr,k+1,1)                     * (  -znew(knew)+0.5*deltaz  - gdepw_0(i,jpr,k+1,1) )
        endif
        !!
        do k=2,mz-1
          !- target level fully embedded in an original level :
          if ( -znew(knew)-0.5*deltaz .ge. gdepw_0(i,jpr,k,1) .and. -znew(knew)+0.5*deltaz .le. gdepw_0(i,jpr,k+1,1) ) then
            ISOMIP_TXZ(i-1,knew,l) = ISOMIP_TXZ(i-1,knew,l) + tmask(i,jpr,k,1) * toce(i,jpr,k,l) * deltaz
            ISOMIP_SXZ(i-1,knew,l) = ISOMIP_SXZ(i-1,knew,l) + tmask(i,jpr,k,1) * soce(i,jpr,k,l) * deltaz
            tmp_XZ                 = tmp_XZ                 + tmask(i,jpr,k,1)                   * deltaz
          !- target level over 2 original levels :
          !- NB: here a target level can be over a maximum of three original level (if thinner target levels, adapt script)
          elseif ( -znew(knew)-0.5*deltaz .le. gdepw_0(i,jpr,k+1,1) .and. -znew(knew)+0.5*deltaz .ge. gdepw_0(i,jpr,k+1,1) .and. -znew(knew)-0.5*deltaz .gt. gdepw_0(i,jpr,k,1) ) then
            ISOMIP_TXZ(i-1,knew,l) = ISOMIP_TXZ(i-1,knew,l) + tmask(i,jpr,k  ,1) * toce(i,jpr,k  ,l) * ( gdepw_0(i,jpr,k+1,1) - (-znew(knew)-0.5*deltaz) ) &
            &                                               + tmask(i,jpr,k+1,1) * toce(i,jpr,k+1,l) * (  -znew(knew)+0.5*deltaz  - gdepw_0(i,jpr,k+1,1) )
            ISOMIP_SXZ(i-1,knew,l) = ISOMIP_SXZ(i-1,knew,l) + tmask(i,jpr,k  ,1) * soce(i,jpr,k  ,l) * ( gdepw_0(i,jpr,k+1,1) - (-znew(knew)-0.5*deltaz) ) &
            &                                               + tmask(i,jpr,k+1,1) * soce(i,jpr,k+1,l) * (  -znew(knew)+0.5*deltaz  - gdepw_0(i,jpr,k+1,1) )
            tmp_XZ                 = tmp_XZ                 + tmask(i,jpr,k  ,1)                     * ( gdepw_0(i,jpr,k+1,1) - (-znew(knew)-0.5*deltaz) ) &
            &                                               + tmask(i,jpr,k+1,1)                     * (  -znew(knew)+0.5*deltaz  - gdepw_0(i,jpr,k+1,1) )
          !- target level over 3 original levels :
          !- NB: here a target level can be over a maximum of three original level (if thinner target levels, adapt script)
          elseif ( -znew(knew)-0.5*deltaz .le. gdepw_0(i,jpr,k+1,1) .and. -znew(knew)+0.5*deltaz .ge. gdepw_0(i,jpr,k+1,1) .and. -znew(knew)-0.5*deltaz .le. gdepw_0(i,jpr,k,1) ) then
            ISOMIP_TXZ(i-1,knew,l) = ISOMIP_TXZ(i-1,knew,l) + tmask(i,jpr,k-1,1) * toce(i,jpr,k-1,l) * ( gdepw_0(i,jpr,k  ,1) - (-znew(knew)-0.5*deltaz) ) &
            &                                               + tmask(i,jpr,k  ,1) * toce(i,jpr,k  ,l) * deltaz                                              &
            &                                               + tmask(i,jpr,k+1,1) * toce(i,jpr,k+1,l) * (  -znew(knew)+0.5*deltaz  - gdepw_0(i,jpr,k+1,1) )
            ISOMIP_SXZ(i-1,knew,l) = ISOMIP_SXZ(i-1,knew,l) + tmask(i,jpr,k-1,1) * soce(i,jpr,k-1,l) * ( gdepw_0(i,jpr,k  ,1) - (-znew(knew)-0.5*deltaz) ) &
            &                                               + tmask(i,jpr,k  ,1) * soce(i,jpr,k  ,l) * deltaz                                              &
            &                                               + tmask(i,jpr,k+1,1) * soce(i,jpr,k+1,l) * (  -znew(knew)+0.5*deltaz  - gdepw_0(i,jpr,k+1,1) )
            tmp_XZ                 = tmp_XZ                 + tmask(i,jpr,k-1,1)                     * ( gdepw_0(i,jpr,k  ,1) - (-znew(knew)-0.5*deltaz) ) &
            &                                               + tmask(i,jpr,k  ,1)                     * deltaz                                              &
            &                                               + tmask(i,jpr,k+1,1)                     * (  -znew(knew)+0.5*deltaz  - gdepw_0(i,jpr,k+1,1) )
          endif
        enddo
      enddo  ! jpr
      if ( tmp_XZ .lt. deltaz ) then  !! we require at least a full cell to fill new cell (the max. being 2.0*deltaz)
        ISOMIP_TXZ(i-1,knew,l) = NF90_FILL_FLOAT
        ISOMIP_SXZ(i-1,knew,l) = NF90_FILL_FLOAT
      else
        ISOMIP_TXZ(i-1,knew,l) = ISOMIP_TXZ(i-1,knew,l) / tmp_XZ
        ISOMIP_SXZ(i-1,knew,l) = ISOMIP_SXZ(i-1,knew,l) / tmp_XZ
      endif
    enddo  ! knew
  enddo  ! i

  !-- YZ section at x = 520 km
  do j=2,my-1
    !-
    do knew=1,ksrf
      ISOMIP_TYZ(j-1,knew,l) = ( tmask(101,j,1,1) * toce(101,j,1,l) + tmask(102,j,1,1) * toce(102,j,1,l) ) / ( tmask(101,j,1,1) + tmask(102,j,1,1) + eps ) &
      &                        + (1 - MIN(1,tmask(101,j,1,1)+tmask(102,j,1,1))) * NF90_FILL_FLOAT
      ISOMIP_SYZ(j-1,knew,l) = ( tmask(101,j,1,1) * soce(101,j,1,l) + tmask(102,j,1,1) * soce(102,j,1,l) ) / ( tmask(101,j,1,1) + tmask(102,j,1,1) + eps ) &
      &                        + (1 - MIN(1,tmask(101,j,1,1)+tmask(102,j,1,1))) * NF90_FILL_FLOAT
    enddo
    !-
    do knew=ksrf+1,mznew
      tmp_YZ = 0.d0
      do ipr=201,202
        k=1
        !- target level fully embedded in an original level :
        if ( -znew(knew)-0.5*deltaz .ge. gdepw_0(ipr,j,k,1) .and. -znew(knew)+0.5*deltaz .le. gdepw_0(ipr,j,k+1,1) ) then
          ISOMIP_TYZ(j-1,knew,l) = ISOMIP_TYZ(j-1,knew,l) + tmask(ipr,j,k,1) * toce(ipr,j,k,l) * deltaz
          ISOMIP_SYZ(j-1,knew,l) = ISOMIP_SYZ(j-1,knew,l) + tmask(ipr,j,k,1) * soce(ipr,j,k,l) * deltaz
          tmp_YZ                 = tmp_YZ                 + tmask(ipr,j,k,1)                   * deltaz
        !- target level over 2 original levels :
        !- NB: here a target level can be over a maximum of three original level (if thinner target levels, adapt script)
        elseif ( -znew(knew)-0.5*deltaz .le. gdepw_0(ipr,j,k+1,1) .and. -znew(knew)+0.5*deltaz .ge. gdepw_0(ipr,j,k+1,1) .and. -znew(knew)-0.5*deltaz .gt. gdepw_0(ipr,j,k,1) ) then
          ISOMIP_TYZ(j-1,knew,l) = ISOMIP_TYZ(j-1,knew,l) + tmask(ipr,j,k  ,1) * toce(ipr,j,k  ,l) * ( gdepw_0(ipr,j,k+1,1) - (-znew(knew)-0.5*deltaz) ) &
          &                                               + tmask(ipr,j,k+1,1) * toce(ipr,j,k+1,l) * (  -znew(knew)+0.5*deltaz  - gdepw_0(ipr,j,k+1,1) )
          ISOMIP_SYZ(j-1,knew,l) = ISOMIP_SYZ(j-1,knew,l) + tmask(ipr,j,k  ,1) * soce(ipr,j,k  ,l) * ( gdepw_0(ipr,j,k+1,1) - (-znew(knew)-0.5*deltaz) ) &
          &                                               + tmask(ipr,j,k+1,1) * soce(ipr,j,k+1,l) * (  -znew(knew)+0.5*deltaz  - gdepw_0(ipr,j,k+1,1) )
          tmp_YZ                 = tmp_YZ                 + tmask(ipr,j,k  ,1)                     * ( gdepw_0(ipr,j,k+1,1) - (-znew(knew)-0.5*deltaz) ) &
          &                                               + tmask(ipr,j,k+1,1)                     * (  -znew(knew)+0.5*deltaz  - gdepw_0(ipr,j,k+1,1) )
        endif
        !!
        do k=2,mz-1
          !- target level fully embedded in an original level :
          if ( -znew(knew)-0.5*deltaz .ge. gdepw_0(ipr,j,k,1) .and. -znew(knew)+0.5*deltaz .le. gdepw_0(ipr,j,k+1,1) ) then
            ISOMIP_TYZ(j-1,knew,l) = ISOMIP_TYZ(j-1,knew,l) + tmask(ipr,j,k,1) * toce(ipr,j,k,l) * deltaz
            ISOMIP_SYZ(j-1,knew,l) = ISOMIP_SYZ(j-1,knew,l) + tmask(ipr,j,k,1) * soce(ipr,j,k,l) * deltaz
            tmp_YZ                 = tmp_YZ                 + tmask(ipr,j,k,1)                   * deltaz
          !- target level over 2 original levels :
          !- NB: here a target level can be over a maximum of three original level (if thinner target levels, adapt script)
          elseif ( -znew(knew)-0.5*deltaz .le. gdepw_0(ipr,j,k+1,1) .and. -znew(knew)+0.5*deltaz .ge. gdepw_0(ipr,j,k+1,1) .and. -znew(knew)-0.5*deltaz .gt. gdepw_0(ipr,j,k,1) ) then
            ISOMIP_TYZ(j-1,knew,l) = ISOMIP_TYZ(j-1,knew,l) + tmask(ipr,j,k  ,1) * toce(ipr,j,k  ,l) * ( gdepw_0(ipr,j,k+1,1) - (-znew(knew)-0.5*deltaz) ) &
            &                                               + tmask(ipr,j,k+1,1) * toce(ipr,j,k+1,l) * (  -znew(knew)+0.5*deltaz  - gdepw_0(ipr,j,k+1,1) )
            ISOMIP_SYZ(j-1,knew,l) = ISOMIP_SYZ(j-1,knew,l) + tmask(ipr,j,k  ,1) * soce(ipr,j,k  ,l) * ( gdepw_0(ipr,j,k+1,1) - (-znew(knew)-0.5*deltaz) ) &
            &                                               + tmask(ipr,j,k+1,1) * soce(ipr,j,k+1,l) * (  -znew(knew)+0.5*deltaz  - gdepw_0(ipr,j,k+1,1) )
            tmp_YZ                 = tmp_YZ                 + tmask(ipr,j,k  ,1)                     * ( gdepw_0(ipr,j,k+1,1) - (-znew(knew)-0.5*deltaz) ) &
            &                                               + tmask(ipr,j,k+1,1)                     * (  -znew(knew)+0.5*deltaz  - gdepw_0(ipr,j,k+1,1) )
          !- target level over 3 original levels :
          !- NB: here a target level can be over a maximum of three original level (if thinner target levels, adapt script)
          elseif ( -znew(knew)-0.5*deltaz .le. gdepw_0(ipr,j,k+1,1) .and. -znew(knew)+0.5*deltaz .ge. gdepw_0(ipr,j,k+1,1) .and. -znew(knew)-0.5*deltaz .le. gdepw_0(ipr,j,k,1) ) then
            ISOMIP_TYZ(j-1,knew,l) = ISOMIP_TYZ(j-1,knew,l) + tmask(ipr,j,k-1,1) * toce(ipr,j,k-1,l) * ( gdepw_0(ipr,j,k  ,1) - (-znew(knew)-0.5*deltaz) ) &
            &                                               + tmask(ipr,j,k  ,1) * toce(ipr,j,k  ,l) * deltaz                                              &
            &                                               + tmask(ipr,j,k+1,1) * toce(ipr,j,k+1,l) * (  -znew(knew)+0.5*deltaz  - gdepw_0(ipr,j,k+1,1) )
            ISOMIP_SYZ(j-1,knew,l) = ISOMIP_SYZ(j-1,knew,l) + tmask(ipr,j,k-1,1) * soce(ipr,j,k-1,l) * ( gdepw_0(ipr,j,k  ,1) - (-znew(knew)-0.5*deltaz) ) &
            &                                               + tmask(ipr,j,k  ,1) * soce(ipr,j,k  ,l) * deltaz                                              &
            &                                               + tmask(ipr,j,k+1,1) * soce(ipr,j,k+1,l) * (  -znew(knew)+0.5*deltaz  - gdepw_0(ipr,j,k+1,1) )
            tmp_YZ                 = tmp_YZ                 + tmask(ipr,j,k-1,1)                     * ( gdepw_0(ipr,j,k  ,1) - (-znew(knew)-0.5*deltaz) ) &
            &                                               + tmask(ipr,j,k  ,1)                     * deltaz                                              &
            &                                               + tmask(ipr,j,k+1,1)                     * (  -znew(knew)+0.5*deltaz  - gdepw_0(ipr,j,k+1,1) )
          endif
        enddo
      enddo ! ipr
      !-
      if ( tmp_YZ .lt. deltaz ) then !! we require at least a full cell to fill new cell (the max. being 2.0*deltaz)
        ISOMIP_TYZ(j-1,knew,l) = NF90_FILL_FLOAT
        ISOMIP_SYZ(j-1,knew,l) = NF90_FILL_FLOAT
      else
        ISOMIP_TYZ(j-1,knew,l) = ISOMIP_TYZ(j-1,knew,l) / tmp_YZ
        ISOMIP_SYZ(j-1,knew,l) = ISOMIP_SYZ(j-1,knew,l) / tmp_YZ
      endif
      !-
    enddo
  enddo

ENDDO

DEALLOCATE( toce, soce, gdept_0, gdepw_0 )

!---------------------------------------
! Average onto coarser ISOMIP+ grid :

nx_CRS=(mx-2)/2
ny_CRS=(my-2)/2

ALLOCATE( CRS_ISOMIP_isf_draft  (nx_CRS,ny_CRS) )
ALLOCATE( CRS_ISOMIP_bathy      (nx_CRS,ny_CRS) )
ALLOCATE( CRS_ISOMIP_melt_rate  (nx_CRS,ny_CRS,mtime_counter) )
ALLOCATE( CRS_ISOMIP_ustar      (nx_CRS,ny_CRS,mtime_counter) )
ALLOCATE( CRS_ISOMIP_therm_driv (nx_CRS,ny_CRS,mtime_counter) )
ALLOCATE( CRS_ISOMIP_halin_driv (nx_CRS,ny_CRS,mtime_counter) )
ALLOCATE( CRS_ISOMIP_Tbot       (nx_CRS,ny_CRS,mtime_counter) )
ALLOCATE( CRS_ISOMIP_Sbot       (nx_CRS,ny_CRS,mtime_counter) )
ALLOCATE( CRS_ISOMIP_TXZ(nx_CRS,mznew,mtime_counter) )
ALLOCATE( CRS_ISOMIP_TYZ(ny_CRS,mznew,mtime_counter) )
ALLOCATE( CRS_ISOMIP_SXZ(nx_CRS,mznew,mtime_counter) )
ALLOCATE( CRS_ISOMIP_SYZ(ny_CRS,mznew,mtime_counter) )

CRS_ISOMIP_isf_draft=0.e0
CRS_ISOMIP_bathy=0.e0
CRS_ISOMIP_melt_rate=0.e0
CRS_ISOMIP_ustar=0.e0
CRS_ISOMIP_therm_driv=0.e0
CRS_ISOMIP_halin_driv=0.e0
CRS_ISOMIP_Tbot=0.e0
CRS_ISOMIP_Sbot=0.e0

do i=1,nx_CRS
do j=1,ny_CRS

  cnt=0
  if ( ISOMIP_isf_draft(2*i-1,2*j-1) .ne. NF90_FILL_FLOAT ) then
    CRS_ISOMIP_isf_draft(i,j) = CRS_ISOMIP_isf_draft(i,j) + ISOMIP_isf_draft(2*i-1,2*j-1)
    cnt=cnt+1
  endif
  if ( ISOMIP_isf_draft(2*i-1,2*j) .ne. NF90_FILL_FLOAT ) then
    CRS_ISOMIP_isf_draft(i,j) = CRS_ISOMIP_isf_draft(i,j) + ISOMIP_isf_draft(2*i-1,2*j)
    cnt=cnt+1
  endif
  if ( ISOMIP_isf_draft(2*i,2*j-1) .ne. NF90_FILL_FLOAT ) then
    CRS_ISOMIP_isf_draft(i,j) = CRS_ISOMIP_isf_draft(i,j) + ISOMIP_isf_draft(2*i,2*j-1)
    cnt=cnt+1
  endif
  if ( ISOMIP_isf_draft(2*i,2*j) .ne. NF90_FILL_FLOAT ) then
    CRS_ISOMIP_isf_draft(i,j) = CRS_ISOMIP_isf_draft(i,j) + ISOMIP_isf_draft(2*i,2*j)
    cnt=cnt+1
  endif
  if ( cnt .ne. 0 ) then
    CRS_ISOMIP_isf_draft(i,j) = CRS_ISOMIP_isf_draft(i,j) / cnt
  else
    CRS_ISOMIP_isf_draft(i,j) = NF90_FILL_FLOAT
  endif

  cnt=0
  if ( ISOMIP_bathy(2*i-1,2*j-1) .ne. NF90_FILL_FLOAT ) then
    CRS_ISOMIP_bathy(i,j) = CRS_ISOMIP_bathy(i,j) + ISOMIP_bathy(2*i-1,2*j-1)
    cnt=cnt+1
  endif
  if ( ISOMIP_bathy(2*i-1,2*j) .ne. NF90_FILL_FLOAT ) then
    CRS_ISOMIP_bathy(i,j) = CRS_ISOMIP_bathy(i,j) + ISOMIP_bathy(2*i-1,2*j)
    cnt=cnt+1
  endif
  if ( ISOMIP_bathy(2*i,2*j-1) .ne. NF90_FILL_FLOAT ) then
    CRS_ISOMIP_bathy(i,j) = CRS_ISOMIP_bathy(i,j) + ISOMIP_bathy(2*i,2*j-1)
    cnt=cnt+1
  endif
  if ( ISOMIP_bathy(2*i,2*j) .ne. NF90_FILL_FLOAT ) then
    CRS_ISOMIP_bathy(i,j) = CRS_ISOMIP_bathy(i,j) + ISOMIP_bathy(2*i,2*j)
    cnt=cnt+1
  endif
  if ( cnt .ne. 0 ) then
    CRS_ISOMIP_bathy(i,j) = CRS_ISOMIP_bathy(i,j) / cnt
  else
    CRS_ISOMIP_bathy(i,j) = NF90_FILL_FLOAT
  endif

  do l=1,mtime_counter

    cnt=0
    if ( ISOMIP_melt_rate(2*i-1,2*j-1,l) .ne. NF90_FILL_FLOAT ) then
      CRS_ISOMIP_melt_rate(i,j,l) = CRS_ISOMIP_melt_rate(i,j,l) + ISOMIP_melt_rate(2*i-1,2*j-1,l)
      cnt=cnt+1
    endif
    if ( ISOMIP_melt_rate(2*i-1,2*j,l) .ne. NF90_FILL_FLOAT ) then
      CRS_ISOMIP_melt_rate(i,j,l) = CRS_ISOMIP_melt_rate(i,j,l) + ISOMIP_melt_rate(2*i-1,2*j,l)
      cnt=cnt+1
    endif
    if ( ISOMIP_melt_rate(2*i,2*j-1,l) .ne. NF90_FILL_FLOAT ) then
      CRS_ISOMIP_melt_rate(i,j,l) = CRS_ISOMIP_melt_rate(i,j,l) + ISOMIP_melt_rate(2*i,2*j-1,l)
      cnt=cnt+1
    endif
    if ( ISOMIP_melt_rate(2*i,2*j,l) .ne. NF90_FILL_FLOAT ) then
      CRS_ISOMIP_melt_rate(i,j,l) = CRS_ISOMIP_melt_rate(i,j,l) + ISOMIP_melt_rate(2*i,2*j,l)
      cnt=cnt+1
    endif
    if ( cnt .ne. 0 ) then
      CRS_ISOMIP_melt_rate(i,j,l) = CRS_ISOMIP_melt_rate(i,j,l) / cnt
    else
      CRS_ISOMIP_melt_rate(i,j,l) = NF90_FILL_FLOAT
    endif

    cnt=0
    if ( ISOMIP_ustar(2*i-1,2*j-1,l) .ne. NF90_FILL_FLOAT ) then
      CRS_ISOMIP_ustar(i,j,l) = CRS_ISOMIP_ustar(i,j,l) + ISOMIP_ustar(2*i-1,2*j-1,l)
      cnt=cnt+1
    endif
    if ( ISOMIP_ustar(2*i-1,2*j,l) .ne. NF90_FILL_FLOAT ) then
      CRS_ISOMIP_ustar(i,j,l) = CRS_ISOMIP_ustar(i,j,l) + ISOMIP_ustar(2*i-1,2*j,l)
      cnt=cnt+1
    endif
    if ( ISOMIP_ustar(2*i,2*j-1,l) .ne. NF90_FILL_FLOAT ) then
      CRS_ISOMIP_ustar(i,j,l) = CRS_ISOMIP_ustar(i,j,l) + ISOMIP_ustar(2*i,2*j-1,l)
      cnt=cnt+1
    endif
    if ( ISOMIP_ustar(2*i,2*j,l) .ne. NF90_FILL_FLOAT ) then
      CRS_ISOMIP_ustar(i,j,l) = CRS_ISOMIP_ustar(i,j,l) + ISOMIP_ustar(2*i,2*j,l)
      cnt=cnt+1
    endif
    if ( cnt .ne. 0 ) then
      CRS_ISOMIP_ustar(i,j,l) = CRS_ISOMIP_ustar(i,j,l) / cnt
    else
      CRS_ISOMIP_ustar(i,j,l) = NF90_FILL_FLOAT
    endif

    cnt=0
    if ( ISOMIP_therm_driv(2*i-1,2*j-1,l) .ne. NF90_FILL_FLOAT ) then
      CRS_ISOMIP_therm_driv(i,j,l) = CRS_ISOMIP_therm_driv(i,j,l) + ISOMIP_therm_driv(2*i-1,2*j-1,l)
      cnt=cnt+1
    endif
    if ( ISOMIP_therm_driv(2*i-1,2*j,l) .ne. NF90_FILL_FLOAT ) then
      CRS_ISOMIP_therm_driv(i,j,l) = CRS_ISOMIP_therm_driv(i,j,l) + ISOMIP_therm_driv(2*i-1,2*j,l)
      cnt=cnt+1
    endif
    if ( ISOMIP_therm_driv(2*i,2*j-1,l) .ne. NF90_FILL_FLOAT ) then
      CRS_ISOMIP_therm_driv(i,j,l) = CRS_ISOMIP_therm_driv(i,j,l) + ISOMIP_therm_driv(2*i,2*j-1,l)
      cnt=cnt+1
    endif
    if ( ISOMIP_therm_driv(2*i,2*j,l) .ne. NF90_FILL_FLOAT ) then
      CRS_ISOMIP_therm_driv(i,j,l) = CRS_ISOMIP_therm_driv(i,j,l) + ISOMIP_therm_driv(2*i,2*j,l)
      cnt=cnt+1
    endif
    if ( cnt .ne. 0 ) then
      CRS_ISOMIP_therm_driv(i,j,l) = CRS_ISOMIP_therm_driv(i,j,l) / cnt
    else
      CRS_ISOMIP_therm_driv(i,j,l) = NF90_FILL_FLOAT
    endif

    cnt=0
    if ( ISOMIP_halin_driv(2*i-1,2*j-1,l) .ne. NF90_FILL_FLOAT ) then
      CRS_ISOMIP_halin_driv(i,j,l) = CRS_ISOMIP_halin_driv(i,j,l) + ISOMIP_halin_driv(2*i-1,2*j-1,l)
      cnt=cnt+1
    endif
    if ( ISOMIP_halin_driv(2*i-1,2*j,l) .ne. NF90_FILL_FLOAT ) then
      CRS_ISOMIP_halin_driv(i,j,l) = CRS_ISOMIP_halin_driv(i,j,l) + ISOMIP_halin_driv(2*i-1,2*j,l)
      cnt=cnt+1
    endif
    if ( ISOMIP_halin_driv(2*i,2*j-1,l) .ne. NF90_FILL_FLOAT ) then
      CRS_ISOMIP_halin_driv(i,j,l) = CRS_ISOMIP_halin_driv(i,j,l) + ISOMIP_halin_driv(2*i,2*j-1,l)
      cnt=cnt+1
    endif
    if ( ISOMIP_halin_driv(2*i,2*j,l) .ne. NF90_FILL_FLOAT ) then
      CRS_ISOMIP_halin_driv(i,j,l) = CRS_ISOMIP_halin_driv(i,j,l) + ISOMIP_halin_driv(2*i,2*j,l)
      cnt=cnt+1
    endif
    if ( cnt .ne. 0 ) then
      CRS_ISOMIP_halin_driv(i,j,l) = CRS_ISOMIP_halin_driv(i,j,l) / cnt
    else
      CRS_ISOMIP_halin_driv(i,j,l) = NF90_FILL_FLOAT
    endif

    cnt=0
    if ( ISOMIP_Tbot(2*i-1,2*j-1,l) .ne. NF90_FILL_FLOAT ) then
      CRS_ISOMIP_Tbot(i,j,l) = CRS_ISOMIP_Tbot(i,j,l) + ISOMIP_Tbot(2*i-1,2*j-1,l)
      cnt=cnt+1
    endif
    if ( ISOMIP_Tbot(2*i-1,2*j,l) .ne. NF90_FILL_FLOAT ) then
      CRS_ISOMIP_Tbot(i,j,l) = CRS_ISOMIP_Tbot(i,j,l) + ISOMIP_Tbot(2*i-1,2*j,l)
      cnt=cnt+1
    endif
    if ( ISOMIP_Tbot(2*i,2*j-1,l) .ne. NF90_FILL_FLOAT ) then
      CRS_ISOMIP_Tbot(i,j,l) = CRS_ISOMIP_Tbot(i,j,l) + ISOMIP_Tbot(2*i,2*j-1,l)
      cnt=cnt+1
    endif
    if ( ISOMIP_Tbot(2*i,2*j,l) .ne. NF90_FILL_FLOAT ) then
      CRS_ISOMIP_Tbot(i,j,l) = CRS_ISOMIP_Tbot(i,j,l) + ISOMIP_Tbot(2*i,2*j,l)
      cnt=cnt+1
    endif
    if ( cnt .ne. 0 ) then
      CRS_ISOMIP_Tbot(i,j,l) = CRS_ISOMIP_Tbot(i,j,l) / cnt
    else
      CRS_ISOMIP_Tbot(i,j,l) = NF90_FILL_FLOAT
    endif

    cnt=0
    if ( ISOMIP_Sbot(2*i-1,2*j-1,l) .ne. NF90_FILL_FLOAT ) then
      CRS_ISOMIP_Sbot(i,j,l) = CRS_ISOMIP_Sbot(i,j,l) + ISOMIP_Sbot(2*i-1,2*j-1,l)
      cnt=cnt+1
    endif
    if ( ISOMIP_Sbot(2*i-1,2*j,l) .ne. NF90_FILL_FLOAT ) then
      CRS_ISOMIP_Sbot(i,j,l) = CRS_ISOMIP_Sbot(i,j,l) + ISOMIP_Sbot(2*i-1,2*j,l)
      cnt=cnt+1
    endif
    if ( ISOMIP_Sbot(2*i,2*j-1,l) .ne. NF90_FILL_FLOAT ) then
      CRS_ISOMIP_Sbot(i,j,l) = CRS_ISOMIP_Sbot(i,j,l) + ISOMIP_Sbot(2*i,2*j-1,l)
      cnt=cnt+1
    endif
    if ( ISOMIP_Sbot(2*i,2*j,l) .ne. NF90_FILL_FLOAT ) then
      CRS_ISOMIP_Sbot(i,j,l) = CRS_ISOMIP_Sbot(i,j,l) + ISOMIP_Sbot(2*i,2*j,l)
      cnt=cnt+1
    endif
    if ( cnt .ne. 0 ) then
      CRS_ISOMIP_Sbot(i,j,l) = CRS_ISOMIP_Sbot(i,j,l) / cnt
    else
      CRS_ISOMIP_Sbot(i,j,l) = NF90_FILL_FLOAT
    endif

  enddo

enddo
enddo

!-----

do i=1,nx_CRS

  do k=1,mznew
  do l=1,mtime_counter

    if ( ISOMIP_TXZ(2*i-1,k,l) .ne. NF90_FILL_FLOAT .and. ISOMIP_TXZ(2*i,k,l) .ne. NF90_FILL_FLOAT ) then
      CRS_ISOMIP_TXZ(i,k,l) = 0.50 * ( ISOMIP_TXZ(2*i-1,k,l) + ISOMIP_TXZ(2*i,k,l) )
    elseif ( ISOMIP_TXZ(2*i-1,k,l) .ne. NF90_FILL_FLOAT ) then
      CRS_ISOMIP_TXZ(i,k,l) = ISOMIP_TXZ(2*i-1,k,l)
    elseif ( ISOMIP_TXZ(2*i,k,l) .ne. NF90_FILL_FLOAT ) then
      CRS_ISOMIP_TXZ(i,k,l) = ISOMIP_TXZ(2*i,k,l)
    else
      CRS_ISOMIP_TXZ(i,k,l) = NF90_FILL_FLOAT
    endif

    if ( ISOMIP_SXZ(2*i-1,k,l) .ne. NF90_FILL_FLOAT .and. ISOMIP_SXZ(2*i,k,l) .ne. NF90_FILL_FLOAT ) then
      CRS_ISOMIP_SXZ(i,k,l) = 0.50 * ( ISOMIP_SXZ(2*i-1,k,l) + ISOMIP_SXZ(2*i,k,l) )
    elseif ( ISOMIP_SXZ(2*i-1,k,l) .ne. NF90_FILL_FLOAT ) then
      CRS_ISOMIP_SXZ(i,k,l) = ISOMIP_SXZ(2*i-1,k,l)
    elseif ( ISOMIP_SXZ(2*i,k,l) .ne. NF90_FILL_FLOAT ) then
      CRS_ISOMIP_SXZ(i,k,l) = ISOMIP_SXZ(2*i,k,l)
    else
      CRS_ISOMIP_SXZ(i,k,l) = NF90_FILL_FLOAT
    endif

  enddo
  enddo

enddo

!-----

do j=1,ny_CRS

  do k=1,mznew
  do l=1,mtime_counter

    if ( ISOMIP_TYZ(2*j-1,k,l) .ne. NF90_FILL_FLOAT .and. ISOMIP_TYZ(2*j,k,l) .ne. NF90_FILL_FLOAT ) then
      CRS_ISOMIP_TYZ(j,k,l) = 0.50 * ( ISOMIP_TYZ(2*j-1,k,l) + ISOMIP_TYZ(2*j,k,l) )
    elseif ( ISOMIP_TYZ(2*j-1,k,l) .ne. NF90_FILL_FLOAT ) then
      CRS_ISOMIP_TYZ(j,k,l) = ISOMIP_TYZ(2*j-1,k,l)
    elseif ( ISOMIP_TYZ(2*j,k,l) .ne. NF90_FILL_FLOAT ) then
      CRS_ISOMIP_TYZ(j,k,l) = ISOMIP_TYZ(2*j,k,l)
    else
      CRS_ISOMIP_TYZ(j,k,l) = NF90_FILL_FLOAT
    endif

    if ( ISOMIP_SYZ(2*j-1,k,l) .ne. NF90_FILL_FLOAT .and. ISOMIP_SYZ(2*j,k,l) .ne. NF90_FILL_FLOAT ) then
      CRS_ISOMIP_SYZ(j,k,l) = 0.50 * ( ISOMIP_SYZ(2*j-1,k,l) + ISOMIP_SYZ(2*j,k,l) )
    elseif ( ISOMIP_SYZ(2*j-1,k,l) .ne. NF90_FILL_FLOAT ) then
      CRS_ISOMIP_SYZ(j,k,l) = ISOMIP_SYZ(2*j-1,k,l)
    elseif ( ISOMIP_SYZ(2*j,k,l) .ne. NF90_FILL_FLOAT ) then
      CRS_ISOMIP_SYZ(j,k,l) = ISOMIP_SYZ(2*j,k,l)
    else
      CRS_ISOMIP_SYZ(j,k,l) = NF90_FILL_FLOAT
    endif

  enddo
  enddo

enddo


!---------------------------------------
! Writing output file :

write(*,*) 'Writing ', TRIM(file_out)

status = NF90_CREATE(TRIM(file_out),NF90_NOCLOBBER,fidM)
call erreur(status,.TRUE.,'create new file')

status = NF90_DEF_DIM(fidM,"nTime",NF90_UNLIMITED,dimID_time_counter)
call erreur(status,.TRUE.,"def_dimID_time_counter")
status = NF90_DEF_DIM(fidM,"nz",mznew,dimID_z)
call erreur(status,.TRUE.,"def_dimID_z")
status = NF90_DEF_DIM(fidM,"nx",nx_CRS,dimID_x)
call erreur(status,.TRUE.,"def_dimID_x")
status = NF90_DEF_DIM(fidM,"ny",ny_CRS,dimID_y)
call erreur(status,.TRUE.,"def_dimID_y")

status = NF90_DEF_VAR(fidM,"salinityYZ",NF90_FLOAT,(/dimID_y,dimID_z,dimID_time_counter/),salinityYZ_ID)
call erreur(status,.TRUE.,"def_var_salinityYZ_ID")
status = NF90_DEF_VAR(fidM,"temperatureYZ",NF90_FLOAT,(/dimID_y,dimID_z,dimID_time_counter/),temperatureYZ_ID)
call erreur(status,.TRUE.,"def_var_temperatureYZ_ID")
status = NF90_DEF_VAR(fidM,"salinityXZ",NF90_FLOAT,(/dimID_x,dimID_z,dimID_time_counter/),salinityXZ_ID)
call erreur(status,.TRUE.,"def_var_salinityXZ_ID")
status = NF90_DEF_VAR(fidM,"temperatureXZ",NF90_FLOAT,(/dimID_x,dimID_z,dimID_time_counter/),temperatureXZ_ID)
call erreur(status,.TRUE.,"def_var_temperatureXZ_ID")
status = NF90_DEF_VAR(fidM,"bottomSalinity",NF90_FLOAT,(/dimID_x,dimID_y,dimID_time_counter/),bottomSalinity_ID)
call erreur(status,.TRUE.,"def_var_bottomSalinity_ID")
status = NF90_DEF_VAR(fidM,"bottomTemperature",NF90_FLOAT,(/dimID_x,dimID_y,dimID_time_counter/),bottomTemperature_ID)
call erreur(status,.TRUE.,"def_var_bottomTemperature_ID")
status = NF90_DEF_VAR(fidM,"halineDriving",NF90_FLOAT,(/dimID_x,dimID_y,dimID_time_counter/),halineDriving_ID)
call erreur(status,.TRUE.,"def_var_halineDriving_ID")
status = NF90_DEF_VAR(fidM,"thermalDriving",NF90_FLOAT,(/dimID_x,dimID_y,dimID_time_counter/),thermalDriving_ID)
call erreur(status,.TRUE.,"def_var_thermalDriving_ID")
status = NF90_DEF_VAR(fidM,"frictionVelocity",NF90_FLOAT,(/dimID_x,dimID_y,dimID_time_counter/),frictionVelocity_ID)
call erreur(status,.TRUE.,"def_var_frictionVelocity_ID")
status = NF90_DEF_VAR(fidM,"meltRate",NF90_FLOAT,(/dimID_x,dimID_y,dimID_time_counter/),meltRate_ID)
call erreur(status,.TRUE.,"def_var_meltRate_ID")
status = NF90_DEF_VAR(fidM,"bathymetry",NF90_FLOAT,(/dimID_x,dimID_y/),bathymetry_ID)
call erreur(status,.TRUE.,"def_var_bathymetry_ID")
status = NF90_DEF_VAR(fidM,"iceDraft",NF90_FLOAT,(/dimID_x,dimID_y/),iceDraft_ID)
call erreur(status,.TRUE.,"def_var_iceDraft_ID")
status = NF90_DEF_VAR(fidM,"meanSalinity",NF90_FLOAT,(/dimID_time_counter/),meanSalinity_ID)
call erreur(status,.TRUE.,"def_var_meanSalinity_ID")
status = NF90_DEF_VAR(fidM,"meanTemperature",NF90_FLOAT,(/dimID_time_counter/),meanTemperature_ID)
call erreur(status,.TRUE.,"def_var_meanTemperature_ID")
status = NF90_DEF_VAR(fidM,"totalOceanVolume",NF90_FLOAT,(/dimID_time_counter/),totalOceanVolume_ID)
call erreur(status,.TRUE.,"def_var_totalOceanVolume_ID")
status = NF90_DEF_VAR(fidM,"totalMeltFlux",NF90_FLOAT,(/dimID_time_counter/),totalMeltFlux_ID)
call erreur(status,.TRUE.,"def_var_totalMeltFlux_ID")
status = NF90_DEF_VAR(fidM,"meanMeltRate",NF90_FLOAT,(/dimID_time_counter/),meanMeltRate_ID)
call erreur(status,.TRUE.,"def_var_meanMeltRate_ID")

status = NF90_PUT_ATT(fidM,salinityYZ_ID,"description","salinity slice in y-z plane through x = 500 km")
call erreur(status,.TRUE.,"put_att_salinityYZ_ID")
status = NF90_PUT_ATT(fidM,salinityYZ_ID,"units","PSU")
call erreur(status,.TRUE.,"put_att_salinityYZ_ID")
status = NF90_PUT_ATT(fidM,salinityYZ_ID,"_FillValue",NF90_FILL_FLOAT)
call erreur(status,.TRUE.,"put_att_salinityYZ_ID")
!-
status = NF90_PUT_ATT(fidM,temperatureYZ_ID,"description","temperature slice in y-z plane through x = 500 km")
call erreur(status,.TRUE.,"put_att_temperatureYZ_ID")
status = NF90_PUT_ATT(fidM,temperatureYZ_ID,"units","deg C")
call erreur(status,.TRUE.,"put_att_temperatureYZ_ID")
status = NF90_PUT_ATT(fidM,temperatureYZ_ID,"_FillValue",NF90_FILL_FLOAT)
call erreur(status,.TRUE.,"put_att_temperatureYZ_ID")
!-
status = NF90_PUT_ATT(fidM,salinityXZ_ID,"description","salinity slice in x-z plane through the center of the domain (y = 40 km)")
call erreur(status,.TRUE.,"put_att_salinityXZ_ID")
status = NF90_PUT_ATT(fidM,salinityXZ_ID,"units","PSU")
call erreur(status,.TRUE.,"put_att_salinityXZ_ID")
status = NF90_PUT_ATT(fidM,salinityXZ_ID,"_FillValue",NF90_FILL_FLOAT)
call erreur(status,.TRUE.,"put_att_salinityXZ_ID")
!-
status = NF90_PUT_ATT(fidM,temperatureXZ_ID,"description","temperature slice in x-z plane through the center of the domain (y = 40 km)")
call erreur(status,.TRUE.,"put_att_temperatureXZ_ID")
status = NF90_PUT_ATT(fidM,temperatureXZ_ID,"units","deg C")
call erreur(status,.TRUE.,"put_att_temperatureXZ_ID")
status = NF90_PUT_ATT(fidM,temperatureXZ_ID,"_FillValue",NF90_FILL_FLOAT)
call erreur(status,.TRUE.,"put_att_temperatureXZ_ID")
!-
status = NF90_PUT_ATT(fidM,bottomSalinity_ID,"description","salinity in the bottom grid cell of each ocean column")
call erreur(status,.TRUE.,"put_att_bottomSalinity_ID")
status = NF90_PUT_ATT(fidM,bottomSalinity_ID,"units","PSU")
call erreur(status,.TRUE.,"put_att_bottomSalinity_ID")
status = NF90_PUT_ATT(fidM,bottomSalinity_ID,"_FillValue",NF90_FILL_FLOAT)
call erreur(status,.TRUE.,"put_att_bottomSalinity_ID")
!-
status = NF90_PUT_ATT(fidM,bottomTemperature_ID,"description","temperature in the bottom grid cell of each ocean column")
call erreur(status,.TRUE.,"put_att_bottomTemperature_ID")
status = NF90_PUT_ATT(fidM,bottomTemperature_ID,"units","deg C")
call erreur(status,.TRUE.,"put_att_bottomTemperature_ID")
status = NF90_PUT_ATT(fidM,bottomTemperature_ID,"_FillValue",NF90_FILL_FLOAT)
call erreur(status,.TRUE.,"put_att_bottomTemperature_ID")
!-
status = NF90_PUT_ATT(fidM,halineDriving_ID,"description","haline driving used in the melt calculation")
call erreur(status,.TRUE.,"put_att_halineDriving_ID")
status = NF90_PUT_ATT(fidM,halineDriving_ID,"units","PSU")
call erreur(status,.TRUE.,"put_att_halineDriving_ID")
status = NF90_PUT_ATT(fidM,halineDriving_ID,"_FillValue",NF90_FILL_FLOAT)
call erreur(status,.TRUE.,"put_att_halineDriving_ID")
!-
status = NF90_PUT_ATT(fidM,thermalDriving_ID,"description","thermal driving used in the melt calculation")
call erreur(status,.TRUE.,"put_att_thermalDriving_ID")
status = NF90_PUT_ATT(fidM,thermalDriving_ID,"units","deg C")
call erreur(status,.TRUE.,"put_att_thermalDriving_ID")
status = NF90_PUT_ATT(fidM,thermalDriving_ID,"_FillValue",NF90_FILL_FLOAT)
call erreur(status,.TRUE.,"put_att_thermalDriving_ID")
!-
status = NF90_PUT_ATT(fidM,frictionVelocity_ID,"description","friction velocity u* used in melt calculations")
call erreur(status,.TRUE.,"put_att_frictionVelocity_ID")
status = NF90_PUT_ATT(fidM,frictionVelocity_ID,"units","m/s")
call erreur(status,.TRUE.,"put_att_frictionVelocity_ID")
status = NF90_PUT_ATT(fidM,frictionVelocity_ID,"_FillValue",NF90_FILL_FLOAT)
call erreur(status,.TRUE.,"put_att_frictionVelocity_ID")
!-
status = NF90_PUT_ATT(fidM,meltRate_ID,"description","melt rate, positive for melting")
call erreur(status,.TRUE.,"put_att_meltRate_ID")
status = NF90_PUT_ATT(fidM,meltRate_ID,"units","m/s")
call erreur(status,.TRUE.,"put_att_meltRate_ID")
status = NF90_PUT_ATT(fidM,meltRate_ID,"_FillValue",NF90_FILL_FLOAT)
call erreur(status,.TRUE.,"put_att_meltRate_ID")
!-
status = NF90_PUT_ATT(fidM,bathymetry_ID,"description","elevation of the bathymetry")
call erreur(status,.TRUE.,"put_att_bathymetry_ID")
status = NF90_PUT_ATT(fidM,bathymetry_ID,"units","m")
call erreur(status,.TRUE.,"put_att_bathymetry_ID")
status = NF90_PUT_ATT(fidM,bathymetry_ID,"_FillValue",NF90_FILL_FLOAT)
call erreur(status,.TRUE.,"put_att_bathymetry_ID")
!-
status = NF90_PUT_ATT(fidM,iceDraft_ID,"description","elevation of the ice-ocean interface")
call erreur(status,.TRUE.,"put_att_iceDraft_ID")
status = NF90_PUT_ATT(fidM,iceDraft_ID,"units","m")
call erreur(status,.TRUE.,"put_att_iceDraft_ID")
status = NF90_PUT_ATT(fidM,iceDraft_ID,"_FillValue",NF90_FILL_FLOAT)
call erreur(status,.TRUE.,"put_att_iceDraft_ID")
!-
status = NF90_PUT_ATT(fidM,meanSalinity_ID,"description","the salinity averaged over the ocean volume")
call erreur(status,.TRUE.,"put_att_meanSalinity_ID")
status = NF90_PUT_ATT(fidM,meanSalinity_ID,"units","PSU")
call erreur(status,.TRUE.,"put_att_meanSalinity_ID")
!-
status = NF90_PUT_ATT(fidM,meanTemperature_ID,"description","the potential temperature averaged over the ocean volume")
call erreur(status,.TRUE.,"put_att_meanTemperature_ID")
status = NF90_PUT_ATT(fidM,meanTemperature_ID,"units","deg C")
call erreur(status,.TRUE.,"put_att_meanTemperature_ID")
!-
status = NF90_PUT_ATT(fidM,totalOceanVolume_ID,"description","total volume of ocean")
call erreur(status,.TRUE.,"put_att_totalOceanVolume_ID")
status = NF90_PUT_ATT(fidM,totalOceanVolume_ID,"units","m^3")
call erreur(status,.TRUE.,"put_att_totalOceanVolume_ID")
!-
status = NF90_PUT_ATT(fidM,totalMeltFlux_ID,"description","total flux of melt water summed over area of floating ice, positive for melting")
call erreur(status,.TRUE.,"put_att_totalMeltFlux_ID")
status = NF90_PUT_ATT(fidM,totalMeltFlux_ID,"units","kg/s")
call erreur(status,.TRUE.,"put_att_totalMeltFlux_ID")
!-
status = NF90_PUT_ATT(fidM,meanMeltRate_ID,"description","mean melt rate averaged over area of floating ice, positive for melting")
call erreur(status,.TRUE.,"put_att_meanMeltRate_ID")
status = NF90_PUT_ATT(fidM,meanMeltRate_ID,"units","m/s")
call erreur(status,.TRUE.,"put_att_meanMeltRate_ID")

status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"experiment","ISOMIP+ TYPb_EEEE")
call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"contact","nicolas.jourdain@univ-grenoble-alpes.fr")
call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"institute","IGE-CNRS, Grenoble, France")
call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"model","NEMO3.6_r6402")
call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"repository","http://www.nemo-ocean.eu")
call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"vertical_coordinates","Z* (variable-volume levels with partial steps)")
call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"momentum_advection","Vector form")
call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"tracer_advection","2nd order Flux Corrected Transport (FCT) scheme")
call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"lateral_diffusion_for_momentum","Laplacian, along geopotential (horizontal)")
call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"lateral_diffusion_for_tracers","Laplacian, along iso-neutral")
call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"vertical_diffusion","TKE dependent scheme")
call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"convection","Enhanced vertical diffusivity and viscosity (100 m2/s)")
call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"top_mixed_layer","T,S,u are averaged over the top 10m")
call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"equation_of_state","linear")
call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"GammaT",GGGG)
call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"GammaS",SSSS)
call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"Cd",2.5e-3)
call erreur(status,.TRUE.,"put_att_GLOBAL_ID")

status = NF90_ENDDEF(fidM)
call erreur(status,.TRUE.,"fin_definition")

status = NF90_PUT_VAR(fidM,salinityYZ_ID,CRS_ISOMIP_SYZ)
call erreur(status,.TRUE.,"var_salinityYZ_ID")
status = NF90_PUT_VAR(fidM,temperatureYZ_ID,CRS_ISOMIP_TYZ)
call erreur(status,.TRUE.,"var_temperatureYZ_ID")
status = NF90_PUT_VAR(fidM,salinityXZ_ID,CRS_ISOMIP_SXZ)
call erreur(status,.TRUE.,"var_salinityXZ_ID")
status = NF90_PUT_VAR(fidM,temperatureXZ_ID,CRS_ISOMIP_TXZ)
call erreur(status,.TRUE.,"var_temperatureXZ_ID")
status = NF90_PUT_VAR(fidM,bottomSalinity_ID,CRS_ISOMIP_Sbot)
call erreur(status,.TRUE.,"var_bottomSalinity_ID")
status = NF90_PUT_VAR(fidM,bottomTemperature_ID,CRS_ISOMIP_Tbot)
call erreur(status,.TRUE.,"var_bottomTemperature_ID")
status = NF90_PUT_VAR(fidM,halineDriving_ID,CRS_ISOMIP_halin_driv)
call erreur(status,.TRUE.,"var_halineDriving_ID")
status = NF90_PUT_VAR(fidM,thermalDriving_ID,CRS_ISOMIP_therm_driv)
call erreur(status,.TRUE.,"var_thermalDriving_ID")
status = NF90_PUT_VAR(fidM,frictionVelocity_ID,CRS_ISOMIP_ustar)
call erreur(status,.TRUE.,"var_frictionVelocity_ID")
status = NF90_PUT_VAR(fidM,meltRate_ID,CRS_ISOMIP_melt_rate)
call erreur(status,.TRUE.,"var_meltRate_ID")
status = NF90_PUT_VAR(fidM,bathymetry_ID,CRS_ISOMIP_bathy)
call erreur(status,.TRUE.,"var_bathymetry_ID")
status = NF90_PUT_VAR(fidM,iceDraft_ID,CRS_ISOMIP_isf_draft)
call erreur(status,.TRUE.,"var_iceDraft_ID")
status = NF90_PUT_VAR(fidM,meanSalinity_ID,ISOMIP_Smean)
call erreur(status,.TRUE.,"var_meanSalinity_ID")
status = NF90_PUT_VAR(fidM,meanTemperature_ID,ISOMIP_Tmean)
call erreur(status,.TRUE.,"var_meanTemperature_ID")
status = NF90_PUT_VAR(fidM,totalOceanVolume_ID,ISOMIP_OceVol)
call erreur(status,.TRUE.,"var_totalOceanVolume_ID")
status = NF90_PUT_VAR(fidM,totalMeltFlux_ID,ISOMIP_total_melt)
call erreur(status,.TRUE.,"var_totalMeltFlux_ID")
status = NF90_PUT_VAR(fidM,meanMeltRate_ID,ISOMIP_mean_melt)
call erreur(status,.TRUE.,"var_meanMeltRate_ID")

status = NF90_CLOSE(fidM)
call erreur(status,.TRUE.,"final")

end program postprot


SUBROUTINE erreur(iret, lstop, chaine)
! pour les messages d'erreur
USE netcdf
INTEGER, INTENT(in)                     :: iret
LOGICAL, INTENT(in)                     :: lstop
CHARACTER(LEN=*), INTENT(in)            :: chaine
!
CHARACTER(LEN=80)                       :: message
!
IF ( iret .NE. 0 ) THEN
  WRITE(*,*) 'ROUTINE: ', TRIM(chaine)
  WRITE(*,*) 'ERREUR: ', iret
  message=NF90_STRERROR(iret)
  WRITE(*,*) 'CA VEUT DIRE:',TRIM(message)
  IF ( lstop ) STOP
ENDIF
!
END SUBROUTINE erreur
