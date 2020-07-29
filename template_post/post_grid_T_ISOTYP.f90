program postprot

USE netcdf

IMPLICIT NONE

INTEGER :: fidT, status, dimID_time_counter, dimID_x, dimID_y, dimID_axis_nbounds, dimID_deptht, &
&          mtime_counter, mx, my, maxis_nbounds, mdeptht, toce_ID, ssh_ID, soce_ID, fwfisf_ID,   &
&          sbt_ID, sbs_ID, deptht_bounds_ID, deptht_ID, fidM, fidSBC, iiso, jiso, mxiso, myiso,  &
&          isfgammas_ID, isfgammat_ID, isfhalindr_ID, isfthermdr_ID, i, j, k, l, ipr, jpr,       &
&          toce_e3t_ID, soce_e3t_ID, e3t_ID, isec, jsec

INTEGER :: salinityYZ_ID, temperatureYZ_ID, salinityXZ_ID, temperatureXZ_ID, bottomSalinity_ID,  &
&          bottomTemperature_ID, halineDriving_ID, thermalDriving_ID, frictionVelocity_ID,       &
&          meltRate_ID, bathymetry_ID, iceDraft_ID, meanSalinity_ID, meanTemperature_ID,         &
&          totalOceanVolume_ID, totalMeltFlux_ID, meanMeltRate_ID, glamt_ID, gphit_ID

INTEGER :: fidMSH, dimID_t, dimID_z, mt, mz, tmask_ID, misf_ID, mbathy_ID, x_ID, y_ID, z_ID,     &
&          e3t_1d_ID, gdepw_1d_ID, gdept_1d_ID, e3t_0_ID, ff_ID, e2t_ID, knew, mznew, ksrf,      &
&          e1t_ID, time_counter_ID, gdepw_0_ID, gdept_0_ID, isfdraft_ID

CHARACTER(LEN=120) :: file_in_T, file_in_SBC, file_MSH, file_out

INTEGER*1,ALLOCATABLE,DIMENSION(:,:,:,:) :: tmask

INTEGER*1,ALLOCATABLE,DIMENSION(:,:) :: isfmskORI, watmskORI

INTEGER*2,ALLOCATABLE,DIMENSION(:,:,:) :: misf, mbathy

INTEGER*2,ALLOCATABLE,DIMENSION(:) :: iinfT, isupT, jinfT, jsupT

REAL*4,ALLOCATABLE,DIMENSION(:) :: deptht, ISOTYP_OceVol, time, wxinfT, wxsupT, wyinfT, wysupT, xiso, yiso

REAL*8,ALLOCATABLE,DIMENSION(:) :: znew, ISOTYP_Tmean, ISOTYP_Smean, ISOTYP_mean_melt, ISOTYP_total_melt

REAL*4,ALLOCATABLE,DIMENSION(:,:) :: ISOTYP_isf_draft, ISOTYP_bathy

REAL*8,ALLOCATABLE,DIMENSION(:,:) :: e3t_1d, gdepw_1d, gdept_1d

REAL*4,ALLOCATABLE,DIMENSION(:,:,:) :: isfdraft, e1t, e2t, ssh, sbt, sbs,                       &
&                                      isfgammas, isfgammat, isfhalindr, isfthermdr, fwfisf,    &
&                                      ISOTYP_melt_rate, ISOTYP_ustar, ISOTYP_therm_driv,       &
&                                      ISOTYP_halin_driv, ISOTYP_Tbot, ISOTYP_Sbot, glamt, gphit

REAL*8,ALLOCATABLE,DIMENSION(:,:,:) :: ISOTYP_TXZ, ISOTYP_SXZ, ISOTYP_TYZ, ISOTYP_SYZ, tmp_TXZ, tmp_SXZ, tmp_TYZ, tmp_SYZ

REAL*4,ALLOCATABLE,DIMENSION(:,:,:,:) :: toce, soce, gdepw_0, gdept_0, toce_e3t, soce_e3t, e3t

REAL*8,ALLOCATABLE,DIMENSION(:,:,:,:) :: e3t_0

INTEGER*4, DIMENSION(12) :: daym1

INTEGER*4 :: an, month0

REAL*4 :: Gt

REAL*8 :: isf_area, tmp_XZ, tmp_YZ, eps, deltaz

!---------------------------------------

file_in_T   = 'MONTHLY/ISOTYP-CCCC_monthly_YYYY_grid_T.nc'
file_in_SBC = 'MONTHLY/ISOTYP-CCCC_monthly_YYYY_SBC.nc'
file_MSH    = '../../input/nemo_ISOTYP/mesh_mask_ISOTYP_EEEE.nc'
file_out    = 'MONTHLY/ISOTYP-CCCC_monthly_YYYY_tmpT.nc'

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
ALLOCATE(  glamt(mx,my,mt), gphit(mx,my,mt) )

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
status = NF90_INQ_VARID(fidMSH,"glamt",glamt_ID)
call erreur(status,.TRUE.,"inq_glamt_ID")
status = NF90_INQ_VARID(fidMSH,"gphit",gphit_ID)
call erreur(status,.TRUE.,"inq_gphit_ID")

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
status = NF90_GET_VAR(fidMSH,glamt_ID,glamt)
call erreur(status,.TRUE.,"getvar_glamt")
status = NF90_GET_VAR(fidMSH,gphit_ID,gphit)
call erreur(status,.TRUE.,"getvar_gphit")

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
! Calculate and interpolate variables on ISOTYP+ std output grid

daym1 = (/ 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 /)
! time at the beginning of the period used for monthly average
do l=1,mtime_counter
  time(l) = ( an * 365 + daym1(l) ) * 86400.0
enddo

! ice-shelf mask :
ALLOCATE( isfmskORI(mx,my) )
where( misf(:,:,1) .ge. 2 )
  isfmskORI(:,:) = 1
elsewhere
  isfmskORI(:,:) = 0
endwhere

! mask for grounded ice sheet :
ALLOCATE( watmskORI(mx,my) )
where( misf(:,:,1) .eq. 0 )
  watmskORI(:,:) = 0
elsewhere
  watmskORI(:,:) = 1
endwhere

!===== 1D =====

ALLOCATE( ISOTYP_OceVol (mtime_counter) )
ALLOCATE( ISOTYP_Tmean  (mtime_counter) )
ALLOCATE( ISOTYP_Smean  (mtime_counter) )
ALLOCATE( ISOTYP_mean_melt (mtime_counter) )
ALLOCATE( ISOTYP_total_melt(mtime_counter) )

ISOTYP_OceVol (:) = 0.e0
ISOTYP_Tmean  (:) = 0.e0
ISOTYP_Smean  (:) = 0.e0 
ISOTYP_mean_melt (:) = 0.e0
ISOTYP_total_melt(:) = 0.e0

DO l=1,mtime_counter
  isf_area = 0.e0
  do i=1,mx
  do j=1,my
    do k=1,mz
      ISOTYP_OceVol (l) = ISOTYP_OceVol (l) +      e3t(i,j,k,l) * e1t(i,j,1) * e2t(i,j,1) * tmask(i,j,k,1)
      ISOTYP_Tmean  (l) = ISOTYP_Tmean  (l) + toce_e3t(i,j,k,l) * e1t(i,j,1) * e2t(i,j,1) * tmask(i,j,k,1)
      ISOTYP_Smean  (l) = ISOTYP_Smean  (l) + soce_e3t(i,j,k,l) * e1t(i,j,1) * e2t(i,j,1) * tmask(i,j,k,1)
    enddo
    isf_area = isf_area + e1t(i,j,1) * e2t(i,j,1) * isfmskORI(i,j)
    ISOTYP_total_melt(l) = ISOTYP_total_melt(l) - fwfisf(i,j,l) * e1t(i,j,1) * e2t(i,j,1) * isfmskORI(i,j)   !!  kg/s
  enddo
  enddo
  !-
  ISOTYP_Tmean(l) = ISOTYP_Tmean(l) / ISOTYP_OceVol(l)
  ISOTYP_Smean(l) = ISOTYP_Smean(l) / ISOTYP_OceVol(l)
  !-
  ISOTYP_mean_melt(l) = ISOTYP_total_melt(l) * 1.e-3 / isf_area   !! m.w.e / s
ENDDO

DEALLOCATE( e3t_0 )

!===== 2D =====

!-- Parameters for interpolation :
mxiso=240
myiso=40
!--
ALLOCATE( ISOTYP_isf_draft  (mxiso,myiso) )
ALLOCATE( ISOTYP_bathy      (mxiso,myiso) )
ALLOCATE( ISOTYP_melt_rate  (mxiso,myiso,mtime_counter) )
ALLOCATE( ISOTYP_ustar      (mxiso,myiso,mtime_counter) )
ALLOCATE( ISOTYP_therm_driv (mxiso,myiso,mtime_counter) )
ALLOCATE( ISOTYP_halin_driv (mxiso,myiso,mtime_counter) )
ALLOCATE( ISOTYP_Tbot       (mxiso,myiso,mtime_counter) )
ALLOCATE( ISOTYP_Sbot       (mxiso,myiso,mtime_counter) )
ALLOCATE( xiso(mxiso), yiso(myiso) )
!--
ALLOCATE( iinfT(mxiso), jinfT(myiso), isupT(mxiso), jsupT(myiso) )
ALLOCATE( wxinfT(mxiso), wyinfT(myiso), wxsupT(mxiso), wysupT(myiso) )
!--
do iiso=1,mxiso
  xiso(iiso) = 319000.000000 + 2000.000000 * iiso
  do i=1,mx-1
    if ( xiso(iiso) .ge. glamt(i,1,1)*1.e3 .and. xiso(iiso) .lt. glamt(i+1,1,1)*1.e3 ) then
      iinfT(iiso)  = i
      isupT(iiso)  = i+1
      wxinfT(iiso) = ( glamt(i+1,1,1)*1.e3 - xiso(iiso) ) / ( glamt(i+1,1,1)*1.e3 - glamt(i,1,1)*1.e3 )
      wxsupT(iiso) = ( xiso(iiso) - glamt(i,1,1)*1.e3   ) / ( glamt(i+1,1,1)*1.e3 - glamt(i,1,1)*1.e3 )
    endif
  enddo
  if ( xiso(iiso) .lt. glamt(1,1,1)*1.e3 ) then
    iinfT(iiso)  = 1
    isupT(iiso)  = 1
    wxinfT(iiso) = 1.000000
    wxsupT(iiso) = 0.000000
  elseif ( xiso(iiso) .ge. glamt(mx,1,1)*1.e3 ) then
    iinfT(iiso)  = mx
    isupT(iiso)  = mx
    wxinfT(iiso) = 1.000000
    wxsupT(iiso) = 0.000000
  endif
enddo
do jiso=1,myiso
  yiso(jiso) = -1000.000000 + 2000.000000 * jiso
  do j=1,my-1
    if ( yiso(jiso) .ge. gphit(1,j,1)*1.e3 .and. yiso(jiso) .lt. gphit(1,j+1,1)*1.e3 ) then
      jinfT(jiso)  = j
      jsupT(jiso)  = j+1
      wyinfT(jiso) = ( gphit(1,j+1,1)*1.e3 - yiso(jiso) ) / ( gphit(1,j+1,1)*1.e3 - gphit(1,j,1)*1.e3 )
      wysupT(jiso) = ( yiso(jiso) - gphit(1,j,1)*1.e3   ) / ( gphit(1,j+1,1)*1.e3 - gphit(1,j,1)*1.e3 )
    endif
  enddo
  if ( yiso(jiso) .lt. gphit(1,1,1)*1.e3 ) then
    jinfT(jiso)  = 1
    jsupT(jiso)  = 1
    wyinfT(jiso) = 1.000000
    wysupT(jiso) = 0.000000
  elseif ( yiso(jiso) .ge. gphit(1,my,1)*1.e3 ) then
    jinfT(jiso)  = my
    jsupT(jiso)  = my
    wyinfT(jiso) = 1.000000
    wysupT(jiso) = 0.000000
  endif
enddo

do iiso=1,mxiso
do jiso=1,myiso
    if (   isfmskORI( iinfT(iiso), jinfT(jiso) )          &
       & + isfmskORI( iinfT(iiso), jsupT(jiso) )          &
       & + isfmskORI( isupT(iiso), jinfT(jiso) )          &
       & + isfmskORI( isupT(iiso), jsupT(jiso) )  .gt. 0  ) then
           ISOTYP_isf_draft(iiso,jiso) =  (   gdepw_0(iinfT(iiso),jinfT(jiso),misf(iinfT(iiso),jinfT(jiso),1),1) * isfmskORI(iinfT(iiso),jinfT(jiso)) * wxinfT(iiso) * wyinfT(jiso)    &
           &                                + gdepw_0(iinfT(iiso),jsupT(jiso),misf(iinfT(iiso),jsupT(jiso),1),1) * isfmskORI(iinfT(iiso),jsupT(jiso)) * wxinfT(iiso) * wysupT(jiso)    &
           &                                + gdepw_0(isupT(iiso),jinfT(jiso),misf(isupT(iiso),jinfT(jiso),1),1) * isfmskORI(isupT(iiso),jinfT(jiso)) * wxsupT(iiso) * wyinfT(jiso)    &
           &                                + gdepw_0(isupT(iiso),jsupT(jiso),misf(isupT(iiso),jsupT(jiso),1),1) * isfmskORI(isupT(iiso),jsupT(jiso)) * wxsupT(iiso) * wysupT(jiso) )  &
           &                            / (                                       isfmskORI(iinfT(iiso),jinfT(jiso)) * wxinfT(iiso) * wyinfT(jiso)    &
           &                                +                                     isfmskORI(iinfT(iiso),jsupT(jiso)) * wxinfT(iiso) * wysupT(jiso)    &
           &                                +                                     isfmskORI(isupT(iiso),jinfT(jiso)) * wxsupT(iiso) * wyinfT(jiso)    &
           &                                +                                     isfmskORI(isupT(iiso),jsupT(jiso)) * wxsupT(iiso) * wysupT(jiso) )
    else
           ISOTYP_isf_draft(iiso,jiso) = NF90_FILL_FLOAT
    endif
    !-
    if (   watmskORI( iinfT(iiso), jinfT(jiso) )          &
       & + watmskORI( iinfT(iiso), jsupT(jiso) )          &
       & + watmskORI( isupT(iiso), jinfT(jiso) )          &
       & + watmskORI( isupT(iiso), jsupT(jiso) )  .gt. 0  ) then
           ISOTYP_bathy(iiso,jiso) =  (   gdepw_0(iinfT(iiso),jinfT(jiso),mbathy(iinfT(iiso),jinfT(jiso),1)+1,1) * watmskORI(iinfT(iiso),jinfT(jiso)) * wxinfT(iiso) * wyinfT(jiso)    &
           &                            + gdepw_0(iinfT(iiso),jsupT(jiso),mbathy(iinfT(iiso),jsupT(jiso),1)+1,1) * watmskORI(iinfT(iiso),jsupT(jiso)) * wxinfT(iiso) * wysupT(jiso)    &
           &                            + gdepw_0(isupT(iiso),jinfT(jiso),mbathy(isupT(iiso),jinfT(jiso),1)+1,1) * watmskORI(isupT(iiso),jinfT(jiso)) * wxsupT(iiso) * wyinfT(jiso)    &
           &                            + gdepw_0(isupT(iiso),jsupT(jiso),mbathy(isupT(iiso),jsupT(jiso),1)+1,1) * watmskORI(isupT(iiso),jsupT(jiso)) * wxsupT(iiso) * wysupT(jiso) )  &
           &                        / (                                       watmskORI(iinfT(iiso),jinfT(jiso)) * wxinfT(iiso) * wyinfT(jiso)    &
           &                            +                                     watmskORI(iinfT(iiso),jsupT(jiso)) * wxinfT(iiso) * wysupT(jiso)    &
           &                            +                                     watmskORI(isupT(iiso),jinfT(jiso)) * wxsupT(iiso) * wyinfT(jiso)    &
           &                            +                                     watmskORI(isupT(iiso),jsupT(jiso)) * wxsupT(iiso) * wysupT(jiso) )
    else
           ISOTYP_bathy(iiso,jiso) = NF90_FILL_FLOAT
    endif
enddo
enddo

!=====
DO l=1,mtime_counter

  do iiso=1,mxiso
  do jiso=1,myiso

      if (   isfmskORI( iinfT(iiso), jinfT(jiso) )          &
         & + isfmskORI( iinfT(iiso), jsupT(jiso) )          &
         & + isfmskORI( isupT(iiso), jinfT(jiso) )          &
         & + isfmskORI( isupT(iiso), jsupT(jiso) )  .gt. 0  ) then
          ! 
          ISOTYP_melt_rate(iiso,jiso,l) =  (   fwfisf(iinfT(iiso),jinfT(jiso),l) * isfmskORI(iinfT(iiso),jinfT(jiso)) * wxinfT(iiso) * wyinfT(jiso)              &  !!  m/s
          &                                  + fwfisf(iinfT(iiso),jsupT(jiso),l) * isfmskORI(iinfT(iiso),jsupT(jiso)) * wxinfT(iiso) * wysupT(jiso)              &
          &                                  + fwfisf(isupT(iiso),jinfT(jiso),l) * isfmskORI(isupT(iiso),jinfT(jiso)) * wxsupT(iiso) * wyinfT(jiso)              &
          &                                  + fwfisf(isupT(iiso),jsupT(jiso),l) * isfmskORI(isupT(iiso),jsupT(jiso)) * wxsupT(iiso) * wysupT(jiso) ) * (-1.e-3) &
          &                              / (                                       isfmskORI(iinfT(iiso),jinfT(jiso)) * wxinfT(iiso) * wyinfT(jiso)              &
          &                                  +                                     isfmskORI(iinfT(iiso),jsupT(jiso)) * wxinfT(iiso) * wysupT(jiso)              &
          &                                  +                                     isfmskORI(isupT(iiso),jinfT(jiso)) * wxsupT(iiso) * wyinfT(jiso)              &
          &                                  +                                     isfmskORI(isupT(iiso),jsupT(jiso)) * wxsupT(iiso) * wysupT(jiso) )
          !
          ISOTYP_ustar(iiso,jiso,l)     =  (   isfgammat(iinfT(iiso),jinfT(jiso),l) * isfmskORI(iinfT(iiso),jinfT(jiso)) * wxinfT(iiso) * wyinfT(jiso)           &  !!  m/s
          &                                  + isfgammat(iinfT(iiso),jsupT(jiso),l) * isfmskORI(iinfT(iiso),jsupT(jiso)) * wxinfT(iiso) * wysupT(jiso)           &
          &                                  + isfgammat(isupT(iiso),jinfT(jiso),l) * isfmskORI(isupT(iiso),jinfT(jiso)) * wxsupT(iiso) * wyinfT(jiso)           &
          &                                  + isfgammat(isupT(iiso),jsupT(jiso),l) * isfmskORI(isupT(iiso),jsupT(jiso)) * wxsupT(iiso) * wysupT(jiso) ) / Gt    &
          &                              / (                                          isfmskORI(iinfT(iiso),jinfT(jiso)) * wxinfT(iiso) * wyinfT(jiso)           &
          &                                  +                                        isfmskORI(iinfT(iiso),jsupT(jiso)) * wxinfT(iiso) * wysupT(jiso)           &
          &                                  +                                        isfmskORI(isupT(iiso),jinfT(jiso)) * wxsupT(iiso) * wyinfT(jiso)           &
          &                                  +                                        isfmskORI(isupT(iiso),jsupT(jiso)) * wxsupT(iiso) * wysupT(jiso) )
          !
          ISOTYP_therm_driv(iiso,jiso,l) = (   isfthermdr(iinfT(iiso),jinfT(jiso),l) * isfmskORI(iinfT(iiso),jinfT(jiso)) * wxinfT(iiso) * wyinfT(jiso)          &  !!  degC
          &                                  + isfthermdr(iinfT(iiso),jsupT(jiso),l) * isfmskORI(iinfT(iiso),jsupT(jiso)) * wxinfT(iiso) * wysupT(jiso)          &
          &                                  + isfthermdr(isupT(iiso),jinfT(jiso),l) * isfmskORI(isupT(iiso),jinfT(jiso)) * wxsupT(iiso) * wyinfT(jiso)          &
          &                                  + isfthermdr(isupT(iiso),jsupT(jiso),l) * isfmskORI(isupT(iiso),jsupT(jiso)) * wxsupT(iiso) * wysupT(jiso) )        &
          &                              / (                                           isfmskORI(iinfT(iiso),jinfT(jiso)) * wxinfT(iiso) * wyinfT(jiso)          &
          &                                  +                                         isfmskORI(iinfT(iiso),jsupT(jiso)) * wxinfT(iiso) * wysupT(jiso)          &
          &                                  +                                         isfmskORI(isupT(iiso),jinfT(jiso)) * wxsupT(iiso) * wyinfT(jiso)          &
          &                                  +                                         isfmskORI(isupT(iiso),jsupT(jiso)) * wxsupT(iiso) * wysupT(jiso) )
          !
          ISOTYP_halin_driv(iiso,jiso,l) = (   isfhalindr(iinfT(iiso),jinfT(jiso),l) * isfmskORI(iinfT(iiso),jinfT(jiso)) * wxinfT(iiso) * wyinfT(jiso)          & !!  psu
          &                                  + isfhalindr(iinfT(iiso),jsupT(jiso),l) * isfmskORI(iinfT(iiso),jsupT(jiso)) * wxinfT(iiso) * wysupT(jiso)          &
          &                                  + isfhalindr(isupT(iiso),jinfT(jiso),l) * isfmskORI(isupT(iiso),jinfT(jiso)) * wxsupT(iiso) * wyinfT(jiso)          &
          &                                  + isfhalindr(isupT(iiso),jsupT(jiso),l) * isfmskORI(isupT(iiso),jsupT(jiso)) * wxsupT(iiso) * wysupT(jiso) )        &
          &                              / (                                           isfmskORI(iinfT(iiso),jinfT(jiso)) * wxinfT(iiso) * wyinfT(jiso)          &
          &                                  +                                         isfmskORI(iinfT(iiso),jsupT(jiso)) * wxinfT(iiso) * wysupT(jiso)          &
          &                                  +                                         isfmskORI(isupT(iiso),jinfT(jiso)) * wxsupT(iiso) * wyinfT(jiso)          &
          &                                  +                                         isfmskORI(isupT(iiso),jsupT(jiso)) * wxsupT(iiso) * wysupT(jiso) )
      else
          ISOTYP_melt_rate(iiso,jiso,l)  = NF90_FILL_FLOAT
          ISOTYP_ustar(iiso,jiso,l)      = NF90_FILL_FLOAT
          ISOTYP_therm_driv(iiso,jiso,l) = NF90_FILL_FLOAT
          ISOTYP_halin_driv(iiso,jiso,l) = NF90_FILL_FLOAT
      endif

      !--

      if (   watmskORI( iinfT(iiso), jinfT(jiso) )          &
         & + watmskORI( iinfT(iiso), jsupT(jiso) )          &
         & + watmskORI( isupT(iiso), jinfT(jiso) )          &
         & + watmskORI( isupT(iiso), jsupT(jiso) )  .gt. 0  ) then
          ! 
          ISOTYP_Tbot(iiso,jiso,l) =  (   sbt(iinfT(iiso),jinfT(jiso),l) * watmskORI(iinfT(iiso),jinfT(jiso)) * wxinfT(iiso) * wyinfT(jiso)    &  !!  degC
          &                             + sbt(iinfT(iiso),jsupT(jiso),l) * watmskORI(iinfT(iiso),jsupT(jiso)) * wxinfT(iiso) * wysupT(jiso)    &
          &                             + sbt(isupT(iiso),jinfT(jiso),l) * watmskORI(isupT(iiso),jinfT(jiso)) * wxsupT(iiso) * wyinfT(jiso)    &
          &                             + sbt(isupT(iiso),jsupT(jiso),l) * watmskORI(isupT(iiso),jsupT(jiso)) * wxsupT(iiso) * wysupT(jiso) )  &
          &                         / (                                    watmskORI(iinfT(iiso),jinfT(jiso)) * wxinfT(iiso) * wyinfT(jiso)    &
          &                             +                                  watmskORI(iinfT(iiso),jsupT(jiso)) * wxinfT(iiso) * wysupT(jiso)    &
          &                             +                                  watmskORI(isupT(iiso),jinfT(jiso)) * wxsupT(iiso) * wyinfT(jiso)    &
          &                             +                                  watmskORI(isupT(iiso),jsupT(jiso)) * wxsupT(iiso) * wysupT(jiso) )
          ! 
          ISOTYP_Sbot(iiso,jiso,l) =  (   sbs(iinfT(iiso),jinfT(jiso),l) * watmskORI(iinfT(iiso),jinfT(jiso)) * wxinfT(iiso) * wyinfT(jiso)    &  !!  psu
          &                             + sbs(iinfT(iiso),jsupT(jiso),l) * watmskORI(iinfT(iiso),jsupT(jiso)) * wxinfT(iiso) * wysupT(jiso)    &
          &                             + sbs(isupT(iiso),jinfT(jiso),l) * watmskORI(isupT(iiso),jinfT(jiso)) * wxsupT(iiso) * wyinfT(jiso)    &
          &                             + sbs(isupT(iiso),jsupT(jiso),l) * watmskORI(isupT(iiso),jsupT(jiso)) * wxsupT(iiso) * wysupT(jiso) )  &
          &                         / (                                    watmskORI(iinfT(iiso),jinfT(jiso)) * wxinfT(iiso) * wyinfT(jiso)    &
          &                             +                                  watmskORI(iinfT(iiso),jsupT(jiso)) * wxinfT(iiso) * wysupT(jiso)    &
          &                             +                                  watmskORI(isupT(iiso),jinfT(jiso)) * wxsupT(iiso) * wyinfT(jiso)    &
          &                             +                                  watmskORI(isupT(iiso),jsupT(jiso)) * wxsupT(iiso) * wysupT(jiso) )
          !
      else
          ISOTYP_Tbot(iiso,jiso,l) = NF90_FILL_FLOAT
          ISOTYP_Sbot(iiso,jiso,l) = NF90_FILL_FLOAT
      endif

      !--
  enddo
  enddo

ENDDO


!===== SECTIONS =====

do j=1,my-1
  if ( 40.0 .ge. gphit(1,j,1) .and. 40.0 .lt. gphit(1,j+1,1) ) jsec=j
enddo
do i=1,mx-1
  if ( 520.0 .ge. glamt(i,1,1) .and. 520.0 .lt. glamt(i+1,1,1) ) isec=i
enddo
write(*,*) 'Section XZ (y= 40km) defined for j =', jsec, jsec+1
write(*,*) 'Section YZ (x=520km) defined for i =', isec, isec+1

ALLOCATE( tmp_TXZ(mx,mznew,mtime_counter) )
ALLOCATE( tmp_TYZ(my,mznew,mtime_counter) )
ALLOCATE( tmp_SXZ(mx,mznew,mtime_counter) )
ALLOCATE( tmp_SYZ(my,mznew,mtime_counter) )

ALLOCATE( ISOTYP_TXZ(mxiso,mznew,mtime_counter) )
ALLOCATE( ISOTYP_TYZ(myiso,mznew,mtime_counter) )
ALLOCATE( ISOTYP_SXZ(mxiso,mznew,mtime_counter) )
ALLOCATE( ISOTYP_SYZ(myiso,mznew,mtime_counter) )

tmp_TXZ(:,:,:) = 0.d0
tmp_TYZ(:,:,:) = 0.d0
tmp_SXZ(:,:,:) = 0.d0
tmp_SYZ(:,:,:) = 0.d0

DO l=1,mtime_counter

  !-- XZ section at y = 40 km
  do i=1,mx
    !-
    do knew=1,ksrf
      tmp_TXZ(i,knew,l) = ( tmask(i,jsec,1,1) * toce(i,jsec,1,l) + tmask(i,jsec+1,1,1) * toce(i,jsec+1,1,l) ) / ( tmask(i,jsec,1,1) + tmask(i,jsec+1,1,1) + eps ) &
      &                    + (1 - MIN(1,tmask(i,jsec,1,1)+tmask(i,jsec+1,1,1))) * NF90_FILL_FLOAT
      tmp_SXZ(i,knew,l) = ( tmask(i,jsec,1,1) * soce(i,jsec,1,l) + tmask(i,jsec+1,1,1) * soce(i,jsec+1,1,l) ) / ( tmask(i,jsec,1,1) + tmask(i,jsec+1,1,1) + eps ) &
      &                    + (1 - MIN(1,tmask(i,jsec,1,1)+tmask(i,jsec+1,1,1))) * NF90_FILL_FLOAT
    enddo
    !-
    do knew=ksrf+1,mznew
      tmp_XZ = 0.d0
      do jpr=jsec,jsec+1
        k=1
        !- target level fully embedded in an original level :
        if ( -znew(knew)-0.5*deltaz .ge. gdepw_0(i,jpr,k,1) .and. -znew(knew)+0.5*deltaz .le. gdepw_0(i,jpr,k+1,1) ) then
          tmp_TXZ(i,knew,l) = tmp_TXZ(i,knew,l) + tmask(i,jpr,k,1) * toce(i,jpr,k,l) * deltaz
          tmp_SXZ(i,knew,l) = tmp_SXZ(i,knew,l) + tmask(i,jpr,k,1) * soce(i,jpr,k,l) * deltaz
          tmp_XZ            = tmp_XZ            + tmask(i,jpr,k,1)                   * deltaz
        !- target level over 2 original levels :
        !- NB: here a target level can be over a maximum of three original level (if thinner target levels, adapt script)
        elseif ( -znew(knew)-0.5*deltaz .le. gdepw_0(i,jpr,k+1,1) .and. -znew(knew)+0.5*deltaz .ge. gdepw_0(i,jpr,k+1,1) .and. -znew(knew)-0.5*deltaz .gt. gdepw_0(i,jpr,k,1) ) then
          tmp_TXZ(i,knew,l) = tmp_TXZ(i,knew,l) + tmask(i,jpr,k  ,1) * toce(i,jpr,k  ,l) * ( gdepw_0(i,jpr,k+1,1) - (-znew(knew)-0.5*deltaz) ) &
          &                                     + tmask(i,jpr,k+1,1) * toce(i,jpr,k+1,l) * (  -znew(knew)+0.5*deltaz  - gdepw_0(i,jpr,k+1,1) )
          tmp_SXZ(i,knew,l) = tmp_SXZ(i,knew,l) + tmask(i,jpr,k  ,1) * soce(i,jpr,k  ,l) * ( gdepw_0(i,jpr,k+1,1) - (-znew(knew)-0.5*deltaz) ) &
          &                                     + tmask(i,jpr,k+1,1) * soce(i,jpr,k+1,l) * (  -znew(knew)+0.5*deltaz  - gdepw_0(i,jpr,k+1,1) )
          tmp_XZ            = tmp_XZ            + tmask(i,jpr,k  ,1)                     * ( gdepw_0(i,jpr,k+1,1) - (-znew(knew)-0.5*deltaz) ) &
          &                                     + tmask(i,jpr,k+1,1)                     * (  -znew(knew)+0.5*deltaz  - gdepw_0(i,jpr,k+1,1) )
        endif
        !!
        do k=2,mz-1
          !- target level fully embedded in an original level :
          if ( -znew(knew)-0.5*deltaz .ge. gdepw_0(i,jpr,k,1) .and. -znew(knew)+0.5*deltaz .le. gdepw_0(i,jpr,k+1,1) ) then
            tmp_TXZ(i,knew,l) = tmp_TXZ(i,knew,l) + tmask(i,jpr,k,1) * toce(i,jpr,k,l) * deltaz
            tmp_SXZ(i,knew,l) = tmp_SXZ(i,knew,l) + tmask(i,jpr,k,1) * soce(i,jpr,k,l) * deltaz
            tmp_XZ            = tmp_XZ            + tmask(i,jpr,k,1)                   * deltaz
          !- target level over 2 original levels :
          !- NB: here a target level can be over a maximum of three original level (if thinner target levels, adapt script)
          elseif ( -znew(knew)-0.5*deltaz .le. gdepw_0(i,jpr,k+1,1) .and. -znew(knew)+0.5*deltaz .ge. gdepw_0(i,jpr,k+1,1) .and. -znew(knew)-0.5*deltaz .gt. gdepw_0(i,jpr,k,1) ) then
            tmp_TXZ(i,knew,l) = tmp_TXZ(i,knew,l) + tmask(i,jpr,k  ,1) * toce(i,jpr,k  ,l) * ( gdepw_0(i,jpr,k+1,1) - (-znew(knew)-0.5*deltaz) ) &
            &                                     + tmask(i,jpr,k+1,1) * toce(i,jpr,k+1,l) * (  -znew(knew)+0.5*deltaz  - gdepw_0(i,jpr,k+1,1) )
            tmp_SXZ(i,knew,l) = tmp_SXZ(i,knew,l) + tmask(i,jpr,k  ,1) * soce(i,jpr,k  ,l) * ( gdepw_0(i,jpr,k+1,1) - (-znew(knew)-0.5*deltaz) ) &
            &                                     + tmask(i,jpr,k+1,1) * soce(i,jpr,k+1,l) * (  -znew(knew)+0.5*deltaz  - gdepw_0(i,jpr,k+1,1) )
            tmp_XZ            = tmp_XZ            + tmask(i,jpr,k  ,1)                     * ( gdepw_0(i,jpr,k+1,1) - (-znew(knew)-0.5*deltaz) ) &
            &                                     + tmask(i,jpr,k+1,1)                     * (  -znew(knew)+0.5*deltaz  - gdepw_0(i,jpr,k+1,1) )
          !- target level over 3 original levels :
          !- NB: here a target level can be over a maximum of three original level (if thinner target levels, adapt script)
          elseif ( -znew(knew)-0.5*deltaz .le. gdepw_0(i,jpr,k+1,1) .and. -znew(knew)+0.5*deltaz .ge. gdepw_0(i,jpr,k+1,1) .and. -znew(knew)-0.5*deltaz .le. gdepw_0(i,jpr,k,1) ) then
            tmp_TXZ(i,knew,l) = tmp_TXZ(i,knew,l) + tmask(i,jpr,k-1,1) * toce(i,jpr,k-1,l) * ( gdepw_0(i,jpr,k  ,1) - (-znew(knew)-0.5*deltaz) ) &
            &                                     + tmask(i,jpr,k  ,1) * toce(i,jpr,k  ,l) * deltaz                                              &
            &                                     + tmask(i,jpr,k+1,1) * toce(i,jpr,k+1,l) * (  -znew(knew)+0.5*deltaz  - gdepw_0(i,jpr,k+1,1) )
            tmp_SXZ(i,knew,l) = tmp_SXZ(i,knew,l) + tmask(i,jpr,k-1,1) * soce(i,jpr,k-1,l) * ( gdepw_0(i,jpr,k  ,1) - (-znew(knew)-0.5*deltaz) ) &
            &                                     + tmask(i,jpr,k  ,1) * soce(i,jpr,k  ,l) * deltaz                                              &
            &                                     + tmask(i,jpr,k+1,1) * soce(i,jpr,k+1,l) * (  -znew(knew)+0.5*deltaz  - gdepw_0(i,jpr,k+1,1) )
            tmp_XZ            = tmp_XZ            + tmask(i,jpr,k-1,1)                     * ( gdepw_0(i,jpr,k  ,1) - (-znew(knew)-0.5*deltaz) ) &
            &                                     + tmask(i,jpr,k  ,1)                     * deltaz                                              &
            &                                     + tmask(i,jpr,k+1,1)                     * (  -znew(knew)+0.5*deltaz  - gdepw_0(i,jpr,k+1,1) )
          endif
        enddo
      enddo  ! jpr
      if ( tmp_XZ .lt. deltaz ) then  !! we require at least a full cell to fill new cell (the max. being 2.0*deltaz)
        tmp_TXZ(i,knew,l) = NF90_FILL_FLOAT
        tmp_SXZ(i,knew,l) = NF90_FILL_FLOAT
      else
        tmp_TXZ(i,knew,l) = tmp_TXZ(i,knew,l) / tmp_XZ
        tmp_SXZ(i,knew,l) = tmp_SXZ(i,knew,l) / tmp_XZ
      endif
    enddo  ! knew
  enddo  ! i

  !-- YZ section at x = 520 km
  do j=1,my
    !-
    do knew=1,ksrf
      ISOTYP_TYZ(j,knew,l) = ( tmask(isec,j,1,1) * toce(isec,j,1,l) + tmask(isec+1,j,1,1) * toce(isec+1,j,1,l) ) / ( tmask(isec,j,1,1) + tmask(isec+1,j,1,1) + eps ) &
      &                      + (1 - MIN(1,tmask(isec,j,1,1)+tmask(isec+1,j,1,1))) * NF90_FILL_FLOAT
      ISOTYP_SYZ(j,knew,l) = ( tmask(isec,j,1,1) * soce(isec,j,1,l) + tmask(isec+1,j,1,1) * soce(isec+1,j,1,l) ) / ( tmask(isec,j,1,1) + tmask(isec+1,j,1,1) + eps ) &
      &                      + (1 - MIN(1,tmask(isec,j,1,1)+tmask(isec+1,j,1,1))) * NF90_FILL_FLOAT
    enddo
    !-
    do knew=ksrf+1,mznew
      tmp_YZ = 0.d0
      do ipr=isec,isec+1
        k=1
        !- target level fully embedded in an original level :
        if ( -znew(knew)-0.5*deltaz .ge. gdepw_0(ipr,j,k,1) .and. -znew(knew)+0.5*deltaz .le. gdepw_0(ipr,j,k+1,1) ) then
          tmp_TYZ(j,knew,l) = tmp_TYZ(j,knew,l) + tmask(ipr,j,k,1) * toce(ipr,j,k,l) * deltaz
          tmp_SYZ(j,knew,l) = tmp_SYZ(j,knew,l) + tmask(ipr,j,k,1) * soce(ipr,j,k,l) * deltaz
          tmp_YZ            = tmp_YZ            + tmask(ipr,j,k,1)                   * deltaz
        !- target level over 2 original levels :
        !- NB: here a target level can be over a maximum of three original level (if thinner target levels, adapt script)
        elseif ( -znew(knew)-0.5*deltaz .le. gdepw_0(ipr,j,k+1,1) .and. -znew(knew)+0.5*deltaz .ge. gdepw_0(ipr,j,k+1,1) .and. -znew(knew)-0.5*deltaz .gt. gdepw_0(ipr,j,k,1) ) then
          tmp_TYZ(j,knew,l) = tmp_TYZ(j,knew,l) + tmask(ipr,j,k  ,1) * toce(ipr,j,k  ,l) * ( gdepw_0(ipr,j,k+1,1) - (-znew(knew)-0.5*deltaz) ) &
          &                                     + tmask(ipr,j,k+1,1) * toce(ipr,j,k+1,l) * (  -znew(knew)+0.5*deltaz  - gdepw_0(ipr,j,k+1,1) )
          tmp_SYZ(j,knew,l) = tmp_SYZ(j,knew,l) + tmask(ipr,j,k  ,1) * soce(ipr,j,k  ,l) * ( gdepw_0(ipr,j,k+1,1) - (-znew(knew)-0.5*deltaz) ) &
          &                                     + tmask(ipr,j,k+1,1) * soce(ipr,j,k+1,l) * (  -znew(knew)+0.5*deltaz  - gdepw_0(ipr,j,k+1,1) )
          tmp_YZ            = tmp_YZ            + tmask(ipr,j,k  ,1)                     * ( gdepw_0(ipr,j,k+1,1) - (-znew(knew)-0.5*deltaz) ) &
          &                                     + tmask(ipr,j,k+1,1)                     * (  -znew(knew)+0.5*deltaz  - gdepw_0(ipr,j,k+1,1) )
        endif
        !!
        do k=2,mz-1
          !- target level fully embedded in an original level :
          if ( -znew(knew)-0.5*deltaz .ge. gdepw_0(ipr,j,k,1) .and. -znew(knew)+0.5*deltaz .le. gdepw_0(ipr,j,k+1,1) ) then
            tmp_TYZ(j,knew,l) = tmp_TYZ(j,knew,l) + tmask(ipr,j,k,1) * toce(ipr,j,k,l) * deltaz
            tmp_SYZ(j,knew,l) = tmp_SYZ(j,knew,l) + tmask(ipr,j,k,1) * soce(ipr,j,k,l) * deltaz
            tmp_YZ            = tmp_YZ            + tmask(ipr,j,k,1)                   * deltaz
          !- target level over 2 original levels :
          !- NB: here a target level can be over a maximum of three original level (if thinner target levels, adapt script)
          elseif ( -znew(knew)-0.5*deltaz .le. gdepw_0(ipr,j,k+1,1) .and. -znew(knew)+0.5*deltaz .ge. gdepw_0(ipr,j,k+1,1) .and. -znew(knew)-0.5*deltaz .gt. gdepw_0(ipr,j,k,1) ) then
            tmp_TYZ(j,knew,l) = tmp_TYZ(j,knew,l) + tmask(ipr,j,k  ,1) * toce(ipr,j,k  ,l) * ( gdepw_0(ipr,j,k+1,1) - (-znew(knew)-0.5*deltaz) ) &
            &                                     + tmask(ipr,j,k+1,1) * toce(ipr,j,k+1,l) * (  -znew(knew)+0.5*deltaz  - gdepw_0(ipr,j,k+1,1) )
            tmp_SYZ(j,knew,l) = tmp_SYZ(j,knew,l) + tmask(ipr,j,k  ,1) * soce(ipr,j,k  ,l) * ( gdepw_0(ipr,j,k+1,1) - (-znew(knew)-0.5*deltaz) ) &
            &                                     + tmask(ipr,j,k+1,1) * soce(ipr,j,k+1,l) * (  -znew(knew)+0.5*deltaz  - gdepw_0(ipr,j,k+1,1) )
            tmp_YZ            = tmp_YZ            + tmask(ipr,j,k  ,1)                     * ( gdepw_0(ipr,j,k+1,1) - (-znew(knew)-0.5*deltaz) ) &
            &                                     + tmask(ipr,j,k+1,1)                     * (  -znew(knew)+0.5*deltaz  - gdepw_0(ipr,j,k+1,1) )
          !- target level over 3 original levels :
          !- NB: here a target level can be over a maximum of three original level (if thinner target levels, adapt script)
          elseif ( -znew(knew)-0.5*deltaz .le. gdepw_0(ipr,j,k+1,1) .and. -znew(knew)+0.5*deltaz .ge. gdepw_0(ipr,j,k+1,1) .and. -znew(knew)-0.5*deltaz .le. gdepw_0(ipr,j,k,1) ) then
            tmp_TYZ(j,knew,l) = tmp_TYZ(j,knew,l) + tmask(ipr,j,k-1,1) * toce(ipr,j,k-1,l) * ( gdepw_0(ipr,j,k  ,1) - (-znew(knew)-0.5*deltaz) ) &
            &                                     + tmask(ipr,j,k  ,1) * toce(ipr,j,k  ,l) * deltaz                                              &
            &                                     + tmask(ipr,j,k+1,1) * toce(ipr,j,k+1,l) * (  -znew(knew)+0.5*deltaz  - gdepw_0(ipr,j,k+1,1) )
            tmp_SYZ(j,knew,l) = tmp_SYZ(j,knew,l) + tmask(ipr,j,k-1,1) * soce(ipr,j,k-1,l) * ( gdepw_0(ipr,j,k  ,1) - (-znew(knew)-0.5*deltaz) ) &
            &                                     + tmask(ipr,j,k  ,1) * soce(ipr,j,k  ,l) * deltaz                                              &
            &                                     + tmask(ipr,j,k+1,1) * soce(ipr,j,k+1,l) * (  -znew(knew)+0.5*deltaz  - gdepw_0(ipr,j,k+1,1) )
            tmp_YZ            = tmp_YZ            + tmask(ipr,j,k-1,1)                     * ( gdepw_0(ipr,j,k  ,1) - (-znew(knew)-0.5*deltaz) ) &
            &                                     + tmask(ipr,j,k  ,1)                     * deltaz                                              &
            &                                     + tmask(ipr,j,k+1,1)                     * (  -znew(knew)+0.5*deltaz  - gdepw_0(ipr,j,k+1,1) )
          endif
        enddo
      enddo ! ipr
      !-
      if ( tmp_YZ .lt. deltaz ) then !! we require at least a full cell to fill new cell (the max. being 2.0*deltaz)
        tmp_TYZ(j,knew,l) = NF90_FILL_FLOAT
        tmp_SYZ(j,knew,l) = NF90_FILL_FLOAT
      else
        tmp_TYZ(j,knew,l) = tmp_TYZ(j,knew,l) / tmp_YZ
        tmp_SYZ(j,knew,l) = tmp_SYZ(j,knew,l) / tmp_YZ
      endif
      !-
    enddo

  enddo

  !---

  do iiso=1,mxiso
    do knew=1,mznew
      if ( tmp_TXZ(iinfT(iiso),knew,l) .eq. NF90_FILL_FLOAT .and. tmp_TXZ(isupT(iiso),knew,l) .eq. NF90_FILL_FLOAT ) then
        ISOTYP_TXZ(iiso,knew,l) = NF90_FILL_FLOAT
        ISOTYP_SXZ(iiso,knew,l) = NF90_FILL_FLOAT
      elseif ( tmp_TXZ(iinfT(iiso),knew,l) .eq. NF90_FILL_FLOAT ) then
        ISOTYP_TXZ(iiso,knew,l) = tmp_TXZ(isupT(iiso),knew,l)
        ISOTYP_SXZ(iiso,knew,l) = tmp_SXZ(isupT(iiso),knew,l)
      elseif ( tmp_TXZ(isupT(iiso),knew,l) .eq. NF90_FILL_FLOAT ) then
        ISOTYP_TXZ(iiso,knew,l) = tmp_TXZ(iinfT(iiso),knew,l)
        ISOTYP_SXZ(iiso,knew,l) = tmp_SXZ(iinfT(iiso),knew,l)
      else
        ISOTYP_TXZ(iiso,knew,l) = ( tmp_TXZ(iinfT(iiso),knew,l) * wxinfT(iiso) + tmp_TXZ(isupT(iiso),knew,l) * wxsupT(iiso) ) / ( wxinfT(iiso) + wxsupT(iiso) )
        ISOTYP_SXZ(iiso,knew,l) = ( tmp_SXZ(iinfT(iiso),knew,l) * wxinfT(iiso) + tmp_SXZ(isupT(iiso),knew,l) * wxsupT(iiso) ) / ( wxinfT(iiso) + wxsupT(iiso) )
      endif
    enddo
  enddo

  do jiso=1,myiso
    do knew=1,mznew
      if ( tmp_TYZ(jinfT(jiso),knew,l) .eq. NF90_FILL_FLOAT .and. tmp_TYZ(jsupT(jiso),knew,l) .eq. NF90_FILL_FLOAT ) then
        ISOTYP_TYZ(jiso,knew,l) = NF90_FILL_FLOAT
        ISOTYP_SYZ(jiso,knew,l) = NF90_FILL_FLOAT
      elseif ( tmp_TYZ(jinfT(jiso),knew,l) .eq. NF90_FILL_FLOAT ) then
        ISOTYP_TYZ(jiso,knew,l) = tmp_TYZ(jsupT(jiso),knew,l)
        ISOTYP_SYZ(jiso,knew,l) = tmp_SYZ(jsupT(jiso),knew,l)
      elseif ( tmp_TYZ(jsupT(jiso),knew,l) .eq. NF90_FILL_FLOAT ) then
        ISOTYP_TYZ(jiso,knew,l) = tmp_TYZ(jinfT(jiso),knew,l)
        ISOTYP_SYZ(jiso,knew,l) = tmp_SYZ(jinfT(jiso),knew,l)
      else
        ISOTYP_TYZ(jiso,knew,l) = ( tmp_TYZ(jinfT(jiso),knew,l) * wyinfT(jiso) + tmp_TYZ(jsupT(jiso),knew,l) * wysupT(jiso) ) / ( wyinfT(jiso) + wysupT(jiso) )
        ISOTYP_SYZ(jiso,knew,l) = ( tmp_SYZ(jinfT(jiso),knew,l) * wyinfT(jiso) + tmp_SYZ(jsupT(jiso),knew,l) * wysupT(jiso) ) / ( wyinfT(jiso) + wysupT(jiso) )
      endif
    enddo
  enddo

ENDDO

DEALLOCATE( toce, soce, gdept_0, gdepw_0 )

!---------------------------------------
! Writing output file :

write(*,*) 'Writing ', TRIM(file_out)

status = NF90_CREATE(TRIM(file_out),NF90_NOCLOBBER,fidM)
call erreur(status,.TRUE.,'create new file')

status = NF90_DEF_DIM(fidM,"nTime",NF90_UNLIMITED,dimID_time_counter)
call erreur(status,.TRUE.,"def_dimID_time_counter")
status = NF90_DEF_DIM(fidM,"nz",mznew,dimID_z)
call erreur(status,.TRUE.,"def_dimID_z")
status = NF90_DEF_DIM(fidM,"nx",mxiso,dimID_x)
call erreur(status,.TRUE.,"def_dimID_x")
status = NF90_DEF_DIM(fidM,"ny",myiso,dimID_y)
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

status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"experiment","ISOMIP+ TYP_EEEE")
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
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"lateral_diffusion_for_momentum","Bi-laplacian, along geopotential (horizontal)")
call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"lateral_diffusion_for_tracers","Laplacian, along iso-neutral")
call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"vertical_diffusion","TKE scheme")
call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"convection","Enhanced vertical diffusivity and viscosity (100m2/s)")
call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"top_mixed_layer","T,S,u are averaged over the top 30m")
call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"equation_of_state","EOS80")
call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"GammaT",Gt)
call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"GammaS",Gt/35.0)
call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"Cd",1.0e-3)
call erreur(status,.TRUE.,"put_att_GLOBAL_ID")

status = NF90_ENDDEF(fidM)
call erreur(status,.TRUE.,"fin_definition")

status = NF90_PUT_VAR(fidM,salinityYZ_ID,ISOTYP_SYZ)
call erreur(status,.TRUE.,"var_salinityYZ_ID")
status = NF90_PUT_VAR(fidM,temperatureYZ_ID,ISOTYP_TYZ)
call erreur(status,.TRUE.,"var_temperatureYZ_ID")
status = NF90_PUT_VAR(fidM,salinityXZ_ID,ISOTYP_SXZ)
call erreur(status,.TRUE.,"var_salinityXZ_ID")
status = NF90_PUT_VAR(fidM,temperatureXZ_ID,ISOTYP_TXZ)
call erreur(status,.TRUE.,"var_temperatureXZ_ID")
status = NF90_PUT_VAR(fidM,bottomSalinity_ID,ISOTYP_Sbot)
call erreur(status,.TRUE.,"var_bottomSalinity_ID")
status = NF90_PUT_VAR(fidM,bottomTemperature_ID,ISOTYP_Tbot)
call erreur(status,.TRUE.,"var_bottomTemperature_ID")
status = NF90_PUT_VAR(fidM,halineDriving_ID,ISOTYP_halin_driv)
call erreur(status,.TRUE.,"var_halineDriving_ID")
status = NF90_PUT_VAR(fidM,thermalDriving_ID,ISOTYP_therm_driv)
call erreur(status,.TRUE.,"var_thermalDriving_ID")
status = NF90_PUT_VAR(fidM,frictionVelocity_ID,ISOTYP_ustar)
call erreur(status,.TRUE.,"var_frictionVelocity_ID")
status = NF90_PUT_VAR(fidM,meltRate_ID,ISOTYP_melt_rate)
call erreur(status,.TRUE.,"var_meltRate_ID")
status = NF90_PUT_VAR(fidM,bathymetry_ID,ISOTYP_bathy)
call erreur(status,.TRUE.,"var_bathymetry_ID")
status = NF90_PUT_VAR(fidM,iceDraft_ID,ISOTYP_isf_draft)
call erreur(status,.TRUE.,"var_iceDraft_ID")
status = NF90_PUT_VAR(fidM,meanSalinity_ID,ISOTYP_Smean)
call erreur(status,.TRUE.,"var_meanSalinity_ID")
status = NF90_PUT_VAR(fidM,meanTemperature_ID,ISOTYP_Tmean)
call erreur(status,.TRUE.,"var_meanTemperature_ID")
status = NF90_PUT_VAR(fidM,totalOceanVolume_ID,ISOTYP_OceVol)
call erreur(status,.TRUE.,"var_totalOceanVolume_ID")
status = NF90_PUT_VAR(fidM,totalMeltFlux_ID,ISOTYP_total_melt)
call erreur(status,.TRUE.,"var_totalMeltFlux_ID")
status = NF90_PUT_VAR(fidM,meanMeltRate_ID,ISOTYP_mean_melt)
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
