!------------------------------------------------------------------
! N. JOURDAIN, LGGE-CNRS, FEB. 2016
!
! Used to calculate :
!    - barotropic streamfunction on grid_T
!    - zonal overturning  streamfunction on grid_T
!    - uTBL and vTBL on grid_T
!
!------------------------------------------------------------------
program modif                                         

USE netcdf                                            

IMPLICIT NONE                                         

INTEGER :: fidU, fidV, status, dimID_time_counter, dimID_axis_nbounds, dimID_depthu, dimID_depthv,          &
&          mtime_counter, maxis_nbounds, mdepthu, mdepthv, mx, my, time_counter_bounds_ID,                  &
&          time_ID, utbl_ID, uoce_e3u_ID, vtbl_ID, voce_e3v_ID, depthu_bounds_ID, &
&          depthu_ID, depthv_bounds_ID, depthv_ID, fidM, mznew

INTEGER :: fidMSH, dimID_t, dimID_z, dimID_y, dimID_x, mt, mz, vmask_ID, umask_ID, tmask_ID, misf_ID, mbathy_ID, i, j, k, l,      &
&          e3w_1d_ID, e3t_1d_ID, gdepw_1d_ID, gdept_1d_ID, e3w_0_ID, e3v_0_ID, e3u_0_ID, e3t_0_ID, ff_ID, e2v_ID, e2u_ID, e2t_ID, &
&          e1v_ID, e1u_ID, e1t_ID, time_counter_ID, gdepw_0_ID, gdepv_ID, gdepu_ID, gdept_0_ID, isfdraft_ID, gphiv_ID, gphiu_ID,  &
&          gphit_ID, glamv_ID, glamu_ID, glamt_ID, x_ID, y_ID, z_ID, BSF_ID, MOC_ID, knew

CHARACTER(LEN=120) :: file_MSH, file_in_U, file_in_V, file_out                     

INTEGER*1,ALLOCATABLE,DIMENSION(:,:,:,:) :: vmask, umask, tmask

INTEGER*2,ALLOCATABLE,DIMENSION(:,:,:) :: misf, mbathy

INTEGER*2,ALLOCATABLE,DIMENSION(:) :: kwsup, kwinf

REAL*4,ALLOCATABLE,DIMENSION(:) :: depthu, depthv, znew, time

REAL*4,ALLOCATABLE,DIMENSION(:,:) :: depthu_bounds, depthv_bounds, trpv, trpz, psif, psiz

REAL*8,ALLOCATABLE,DIMENSION(:,:) :: e3w_1d, e3t_1d, gdepw_1d, gdept_1d

REAL*4,ALLOCATABLE,DIMENSION(:,:,:) :: utbl, vtbl, isfdraft, gphiv, gphiu, gphit, glamv, glamu, glamt,   &
&                                      e1u, e2u, e1v, e2v, e1t, e2t, ISOMIP_BSF, ISOMIP_MOC, ISOMIP_uTBL,&
&                                      ISOMIP_vTBL

REAL*4,ALLOCATABLE,DIMENSION(:,:,:,:) :: uoce_e3u, voce_e3v, gdepw_0, gdepv, gdepu, gdept_0

REAL*8,ALLOCATABLE,DIMENSION(:,:,:,:) :: e3u_0, e3v_0  !, e3w_0, e3t_0

INTEGER*4, DIMENSION(12) :: daym1

INTEGER*4 :: an

REAL*4 :: eps

file_MSH   = '../../input/nemo_ISOMIP/mesh_mask_ISOMIP_EEEE.nc'
file_in_U  = 'MONTHLY/ISOMIP-CCCC_monthly_YYYY_grid_U.nc'
file_in_V  = 'MONTHLY/ISOMIP-CCCC_monthly_YYYY_grid_V.nc'
file_out   = 'MONTHLY/ISOMIP-CCCC_monthly_YYYY_tmpUV.nc'

eps=1.e-9

an=YYYY

!---------------------------------------                   
! Read grid_U input file :

write(*,*) 'Reading ', TRIM(file_in_U)
 
status = NF90_OPEN(TRIM(file_in_U),0,fidU)          
call erreur(status,.TRUE.,"read grid U") 
 
status = NF90_INQ_DIMID(fidU,"time_counter",dimID_time_counter)
call erreur(status,.TRUE.,"inq_dimID_time_counter")
status = NF90_INQ_DIMID(fidU,"axis_nbounds",dimID_axis_nbounds)
call erreur(status,.TRUE.,"inq_dimID_axis_nbounds")
status = NF90_INQ_DIMID(fidU,"depthu",dimID_depthu)
call erreur(status,.TRUE.,"inq_dimID_depthu")
status = NF90_INQ_DIMID(fidU,"x_grid_U",dimID_x)
call erreur(status,.TRUE.,"inq_dimID_x")
status = NF90_INQ_DIMID(fidU,"y_grid_U",dimID_y)
call erreur(status,.TRUE.,"inq_dimID_y")

status = NF90_INQUIRE_DIMENSION(fidU,dimID_time_counter,len=mtime_counter)
call erreur(status,.TRUE.,"inq_dim_time_counter")
status = NF90_INQUIRE_DIMENSION(fidU,dimID_axis_nbounds,len=maxis_nbounds)
call erreur(status,.TRUE.,"inq_dim_axis_nbounds")
status = NF90_INQUIRE_DIMENSION(fidU,dimID_depthu,len=mdepthu)
call erreur(status,.TRUE.,"inq_dim_depthu")
status = NF90_INQUIRE_DIMENSION(fidU,dimID_x,len=mx)
call erreur(status,.TRUE.,"inq_dim_x")
status = NF90_INQUIRE_DIMENSION(fidU,dimID_y,len=my)
call erreur(status,.TRUE.,"inq_dim_y")

ALLOCATE(  time( mtime_counter) )
ALLOCATE(  utbl(mx,my,mtime_counter)  ) 
ALLOCATE(  uoce_e3u(mx,my,mdepthu,mtime_counter)  ) 
ALLOCATE(  depthu_bounds(maxis_nbounds,mdepthu)  ) 
ALLOCATE(  depthu(mdepthu)  ) 
 
status = NF90_INQ_VARID(fidU,"utbl",utbl_ID)
call erreur(status,.TRUE.,"inq_utbl_ID")
status = NF90_INQ_VARID(fidU,"uoce_e3u",uoce_e3u_ID)
call erreur(status,.TRUE.,"inq_uoce_e3u_ID")
status = NF90_INQ_VARID(fidU,"depthu_bounds",depthu_bounds_ID)
call erreur(status,.TRUE.,"inq_depthu_bounds_ID")
status = NF90_INQ_VARID(fidU,"depthu",depthu_ID)
call erreur(status,.TRUE.,"inq_depthu_ID")

status = NF90_GET_VAR(fidU,utbl_ID,utbl)
call erreur(status,.TRUE.,"getvar_utbl")
status = NF90_GET_VAR(fidU,uoce_e3u_ID,uoce_e3u)
call erreur(status,.TRUE.,"getvar_uoce_e3u")
status = NF90_GET_VAR(fidU,depthu_bounds_ID,depthu_bounds)
call erreur(status,.TRUE.,"getvar_depthu_bounds")
status = NF90_GET_VAR(fidU,depthu_ID,depthu)
call erreur(status,.TRUE.,"getvar_depthu")
 
status = NF90_CLOSE(fidU)                      
call erreur(status,.TRUE.,"fin_lecture")     

!---------------------------------------                   
! Read grid_V input file :

write(*,*) 'Reading ', TRIM(file_in_V)
 
status = NF90_OPEN(TRIM(file_in_V),0,fidV)          
call erreur(status,.TRUE.,"read grid_V") 
 
status = NF90_INQ_DIMID(fidV,"depthv",dimID_depthv)
call erreur(status,.TRUE.,"inq_dimID_depthv")

status = NF90_INQUIRE_DIMENSION(fidV,dimID_depthv,len=mdepthv)
call erreur(status,.TRUE.,"inq_dim_depthv")
 
ALLOCATE(  vtbl(mx,my,mtime_counter)  ) 
ALLOCATE(  voce_e3v(mx,my,mdepthv,mtime_counter)  ) 
ALLOCATE(  depthv(mdepthv)  ) 
 
status = NF90_INQ_VARID(fidV,"vtbl",vtbl_ID)
call erreur(status,.TRUE.,"inq_vtbl_ID")
status = NF90_INQ_VARID(fidV,"voce_e3v",voce_e3v_ID)
call erreur(status,.TRUE.,"inq_voce_e3v_ID")
status = NF90_INQ_VARID(fidV,"depthv_bounds",depthv_bounds_ID)
call erreur(status,.TRUE.,"inq_depthv_bounds_ID")
status = NF90_INQ_VARID(fidV,"depthv",depthv_ID)
call erreur(status,.TRUE.,"inq_depthv_ID")

status = NF90_GET_VAR(fidV,vtbl_ID,vtbl)
call erreur(status,.TRUE.,"getvar_vtbl")
status = NF90_GET_VAR(fidV,voce_e3v_ID,voce_e3v)
call erreur(status,.TRUE.,"getvar_voce_e3v")
status = NF90_GET_VAR(fidV,depthv_bounds_ID,depthv_bounds)
call erreur(status,.TRUE.,"getvar_depthv_bounds")
status = NF90_GET_VAR(fidV,depthv_ID,depthv)
call erreur(status,.TRUE.,"getvar_depthv")
 
status = NF90_CLOSE(fidV)                      
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

ALLOCATE(  vmask(mx,my,mz,mt)  ) 
ALLOCATE(  umask(mx,my,mz,mt)  ) 
ALLOCATE(  tmask(mx,my,mz,mt)  ) 
ALLOCATE(  misf(mx,my,mt)  ) 
ALLOCATE(  mbathy(mx,my,mt)  ) 
ALLOCATE(  e3w_1d(mz,mt)  ) 
ALLOCATE(  e3t_1d(mz,mt)  ) 
ALLOCATE(  gdepw_1d(mz,mt)  ) 
ALLOCATE(  gdept_1d(mz,mt)  ) 
!ALLOCATE(  e3w_0(mx,my,mz,mt)  ) 
ALLOCATE(  e3v_0(mx,my,mz,mt)  ) 
ALLOCATE(  e3u_0(mx,my,mz,mt)  ) 
!ALLOCATE(  e3t_0(mx,my,mz,mt)  ) 
ALLOCATE(  e2v(mx,my,mt)  ) 
ALLOCATE(  e2u(mx,my,mt)  ) 
ALLOCATE(  e2t(mx,my,mt)  ) 
ALLOCATE(  e1v(mx,my,mt)  ) 
ALLOCATE(  e1u(mx,my,mt)  ) 
ALLOCATE(  e1t(mx,my,mt)  ) 
ALLOCATE(  gdepw_0(mx,my,mz,mt)  ) 
ALLOCATE(  gdepv(mx,my,mz,mt)  ) 
ALLOCATE(  gdepu(mx,my,mz,mt)  ) 
ALLOCATE(  gdept_0(mx,my,mz,mt)  ) 
ALLOCATE(  isfdraft(mx,my,mt)  ) 
ALLOCATE(  gphiv(mx,my,mt)  ) 
ALLOCATE(  gphiu(mx,my,mt)  ) 
ALLOCATE(  gphit(mx,my,mt)  ) 
ALLOCATE(  glamv(mx,my,mt)  ) 
ALLOCATE(  glamu(mx,my,mt)  ) 
ALLOCATE(  glamt(mx,my,mt)  ) 

status = NF90_INQ_VARID(fidMSH,"vmask",vmask_ID)
call erreur(status,.TRUE.,"inq_vmask_ID")
status = NF90_INQ_VARID(fidMSH,"umask",umask_ID)
call erreur(status,.TRUE.,"inq_umask_ID")
status = NF90_INQ_VARID(fidMSH,"tmask",tmask_ID)
call erreur(status,.TRUE.,"inq_tmask_ID")
status = NF90_INQ_VARID(fidMSH,"misf",misf_ID)
call erreur(status,.TRUE.,"inq_misf_ID")
status = NF90_INQ_VARID(fidMSH,"mbathy",mbathy_ID)
call erreur(status,.TRUE.,"inq_mbathy_ID")
status = NF90_INQ_VARID(fidMSH,"e3w_1d",e3w_1d_ID)
call erreur(status,.TRUE.,"inq_e3w_1d_ID")
status = NF90_INQ_VARID(fidMSH,"e3t_1d",e3t_1d_ID)
call erreur(status,.TRUE.,"inq_e3t_1d_ID")
status = NF90_INQ_VARID(fidMSH,"gdepw_1d",gdepw_1d_ID)
call erreur(status,.TRUE.,"inq_gdepw_1d_ID")
status = NF90_INQ_VARID(fidMSH,"gdept_1d",gdept_1d_ID)
call erreur(status,.TRUE.,"inq_gdept_1d_ID")
!status = NF90_INQ_VARID(fidMSH,"e3w_0",e3w_0_ID)
!call erreur(status,.TRUE.,"inq_e3w_0_ID")
status = NF90_INQ_VARID(fidMSH,"e3v_0",e3v_0_ID)
call erreur(status,.TRUE.,"inq_e3v_0_ID")
status = NF90_INQ_VARID(fidMSH,"e3u_0",e3u_0_ID)
call erreur(status,.TRUE.,"inq_e3u_0_ID")
!status = NF90_INQ_VARID(fidMSH,"e3t_0",e3t_0_ID)
!call erreur(status,.TRUE.,"inq_e3t_0_ID")
status = NF90_INQ_VARID(fidMSH,"e2v",e2v_ID)
call erreur(status,.TRUE.,"inq_e2v_ID")
status = NF90_INQ_VARID(fidMSH,"e2u",e2u_ID)
call erreur(status,.TRUE.,"inq_e2u_ID")
status = NF90_INQ_VARID(fidMSH,"e2t",e2t_ID)
call erreur(status,.TRUE.,"inq_e2t_ID")
status = NF90_INQ_VARID(fidMSH,"e1v",e1v_ID)
call erreur(status,.TRUE.,"inq_e1v_ID")
status = NF90_INQ_VARID(fidMSH,"e1u",e1u_ID)
call erreur(status,.TRUE.,"inq_e1u_ID")
status = NF90_INQ_VARID(fidMSH,"e1t",e1t_ID)
call erreur(status,.TRUE.,"inq_e1t_ID")
status = NF90_INQ_VARID(fidMSH,"gdepw_0",gdepw_0_ID)
call erreur(status,.TRUE.,"inq_gdepw_0_ID")
status = NF90_INQ_VARID(fidMSH,"gdepv",gdepv_ID)
call erreur(status,.TRUE.,"inq_gdepv_ID")
status = NF90_INQ_VARID(fidMSH,"gdepu",gdepu_ID)
call erreur(status,.TRUE.,"inq_gdepu_ID")
status = NF90_INQ_VARID(fidMSH,"gdept_0",gdept_0_ID)
call erreur(status,.TRUE.,"inq_gdept_0_ID")
status = NF90_INQ_VARID(fidMSH,"isfdraft",isfdraft_ID)
call erreur(status,.TRUE.,"inq_isfdraft_ID")
status = NF90_INQ_VARID(fidMSH,"gphiv",gphiv_ID)
call erreur(status,.TRUE.,"inq_gphiv_ID")
status = NF90_INQ_VARID(fidMSH,"gphiu",gphiu_ID)
call erreur(status,.TRUE.,"inq_gphiu_ID")
status = NF90_INQ_VARID(fidMSH,"gphit",gphit_ID)
call erreur(status,.TRUE.,"inq_gphit_ID")
status = NF90_INQ_VARID(fidMSH,"glamv",glamv_ID)
call erreur(status,.TRUE.,"inq_glamv_ID")
status = NF90_INQ_VARID(fidMSH,"glamu",glamu_ID)
call erreur(status,.TRUE.,"inq_glamu_ID")
status = NF90_INQ_VARID(fidMSH,"glamt",glamt_ID)
call erreur(status,.TRUE.,"inq_glamt_ID")

status = NF90_GET_VAR(fidMSH,vmask_ID,vmask)
call erreur(status,.TRUE.,"getvar_vmask")
status = NF90_GET_VAR(fidMSH,umask_ID,umask)
call erreur(status,.TRUE.,"getvar_umask")
status = NF90_GET_VAR(fidMSH,tmask_ID,tmask)
call erreur(status,.TRUE.,"getvar_tmask")
status = NF90_GET_VAR(fidMSH,misf_ID,misf)
call erreur(status,.TRUE.,"getvar_misf")
status = NF90_GET_VAR(fidMSH,mbathy_ID,mbathy)
call erreur(status,.TRUE.,"getvar_mbathy")
status = NF90_GET_VAR(fidMSH,e3w_1d_ID,e3w_1d)
call erreur(status,.TRUE.,"getvar_e3w_1d")
status = NF90_GET_VAR(fidMSH,e3t_1d_ID,e3t_1d)
call erreur(status,.TRUE.,"getvar_e3t_1d")
status = NF90_GET_VAR(fidMSH,gdepw_1d_ID,gdepw_1d)
call erreur(status,.TRUE.,"getvar_gdepw_1d")
status = NF90_GET_VAR(fidMSH,gdept_1d_ID,gdept_1d)
call erreur(status,.TRUE.,"getvar_gdept_1d")
!status = NF90_GET_VAR(fidMSH,e3w_0_ID,e3w_0)
!call erreur(status,.TRUE.,"getvar_e3w_0")
status = NF90_GET_VAR(fidMSH,e3v_0_ID,e3v_0)
call erreur(status,.TRUE.,"getvar_e3v_0")
status = NF90_GET_VAR(fidMSH,e3u_0_ID,e3u_0)
call erreur(status,.TRUE.,"getvar_e3u_0")
!status = NF90_GET_VAR(fidMSH,e3t_0_ID,e3t_0)
!call erreur(status,.TRUE.,"getvar_e3t_0")
status = NF90_GET_VAR(fidMSH,e2v_ID,e2v)
call erreur(status,.TRUE.,"getvar_e2v")
status = NF90_GET_VAR(fidMSH,e2u_ID,e2u)
call erreur(status,.TRUE.,"getvar_e2u")
status = NF90_GET_VAR(fidMSH,e2t_ID,e2t)
call erreur(status,.TRUE.,"getvar_e2t")
status = NF90_GET_VAR(fidMSH,e1v_ID,e1v)
call erreur(status,.TRUE.,"getvar_e1v")
status = NF90_GET_VAR(fidMSH,e1u_ID,e1u)
call erreur(status,.TRUE.,"getvar_e1u")
status = NF90_GET_VAR(fidMSH,e1t_ID,e1t)
call erreur(status,.TRUE.,"getvar_e1t")
status = NF90_GET_VAR(fidMSH,gdepw_0_ID,gdepw_0)
call erreur(status,.TRUE.,"getvar_gdepw_0")
status = NF90_GET_VAR(fidMSH,gdepv_ID,gdepv)
call erreur(status,.TRUE.,"getvar_gdepv")
status = NF90_GET_VAR(fidMSH,gdepu_ID,gdepu)
call erreur(status,.TRUE.,"getvar_gdepu")
status = NF90_GET_VAR(fidMSH,gdept_0_ID,gdept_0)
call erreur(status,.TRUE.,"getvar_gdept_0")
status = NF90_GET_VAR(fidMSH,isfdraft_ID,isfdraft)
call erreur(status,.TRUE.,"getvar_isfdraft")
status = NF90_GET_VAR(fidMSH,gphiv_ID,gphiv)
call erreur(status,.TRUE.,"getvar_gphiv")
status = NF90_GET_VAR(fidMSH,gphiu_ID,gphiu)
call erreur(status,.TRUE.,"getvar_gphiu")
status = NF90_GET_VAR(fidMSH,gphit_ID,gphit)
call erreur(status,.TRUE.,"getvar_gphit")
status = NF90_GET_VAR(fidMSH,glamv_ID,glamv)
call erreur(status,.TRUE.,"getvar_glamv")
status = NF90_GET_VAR(fidMSH,glamu_ID,glamu)
call erreur(status,.TRUE.,"getvar_glamu")
status = NF90_GET_VAR(fidMSH,glamt_ID,glamt)
call erreur(status,.TRUE.,"getvar_glamt")

status = NF90_CLOSE(fidMSH)                      
call erreur(status,.TRUE.,"fin_lecture")     

!------------------------------------------------------
! refined output vertical grid :
mznew=144
ALLOCATE( znew(mznew), kwsup(mznew), kwinf(mznew) )
write(*,*) 'New vertical grid for the outputs :'
do knew=1,mznew
  znew(knew) = 2.5 - knew * 5.0
  do k=1,mz-1
    if ( -znew(knew) .ge. gdepw_1d(k,1) .and. -znew(knew) .lt. gdepw_1d(k+1,1) ) then
      kwsup(knew) = k+1
      kwinf(knew) = k
    endif
  enddo
  write(*,109) -znew(knew), gdepw_1d(kwinf(knew),1), kwinf(knew), gdepw_1d(kwsup(knew),1), kwsup(knew)
  109 FORMAT('  new depth = ',f7.2,' m is interpolated between ',f7.2,' m (k=',i2,') and ',f7.2,' m (k=',i2,')')
enddo
write(*,*) ' '

!---------------------------------------                      
! Calculate requested outputs :

daym1 = (/ 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 /)
! time at the beginning of the period used for monthly average
do l=1,mtime_counter
  time(l) = ( an * 365 + daym1(l) ) * 86400.0 
enddo

ALLOCATE( trpv(mx,my), psif(mx,my) )
ALLOCATE( trpz(mx,mz), psiz(mx,mz) )
ALLOCATE( ISOMIP_BSF(mx-2,my-2   ,mtime_counter) )
ALLOCATE( ISOMIP_MOC(mx-2,  mznew,mtime_counter) )
ALLOCATE( ISOMIP_uTBL(mx-2,my-2,mtime_counter), ISOMIP_vTBL(mx-2,my-2,mtime_counter) )

DO l=1,mtime_counter

  trpv(:,:) = 0.e0
  trpz(:,:) = 0.e0
  psif(:,:) = 0.e0
  psiz(:,:) = 0.e0
  !- meridional barotropic transport :
  do k=1,mz
    trpv(:,:) = trpv(:,:) + voce_e3v(:,:,k,l) * e1v(:,:,1)
  enddo
  !- meridionally-averaged zonal transport :
  do j=1,my
    do i=1,mx
    do k=1,mz
      trpz(i,k) = trpz(i,k) + uoce_e3u(i,j,k,l) * e2u(i,j,1)
    enddo
    enddo
  enddo
  ! BSF calculated by integrating transport along X
  do i=2,mx
    psif(i,:) = psif(i-1,:) - trpv(i,:)  ! BSF at F-points
  enddo
  ! MOC calculated by integrating transport along Z
  do k=mz-1,1,-1
    psiz(:,k) = psiz(:,k+1) - trpz(:,k) ! MOC at W-points
  enddo
  !- interpolate BSF and MOC on ISOMIP's T-points 
  !  (NB: for regular grid in X,Y only):
  !  -> directly on ISOMIP's smaller grid ( mx-2 * my-2 )
  do i=2,mx-1
    ! BSF on T-points (NB: for regular grid in X,Y only) -> directly on ISOMIP's smaller grid ( mx-2 * my-2 ):
    do j=2,my-1
      if ( misf(i,j,1) .le. 1 ) then
        ISOMIP_uTBL(i-1,j-1,l) = NF90_FILL_FLOAT
        ISOMIP_vTBL(i-1,j-1,l) = NF90_FILL_FLOAT
        ISOMIP_BSF (i-1,j-1,l) = -0.25 * ( psif(i,j) + psif(i-1,j) + psif(i,j-1) + psif(i-1,j-1) ) ! ISOMIP convention
      else
        ISOMIP_uTBL(i-1,j-1,l) = 0.50 * ( utbl(i,j,l) + utbl(i-1,j,l) )
        ISOMIP_vTBL(i-1,j-1,l) = 0.50 * ( vtbl(i,j,l) + vtbl(i,j-1,l) )
        ISOMIP_BSF (i-1,j-1,l) = -0.25 * ( psif(i,j) + psif(i-1,j) + psif(i,j-1) + psif(i-1,j-1) ) ! ISOMIP convention
      endif
    enddo
    ! MOC (interpolation from W-point to ISOMIP 5-m vertical grid )
    do knew=1,mznew
      ISOMIP_MOC(i-1,knew,l) = (   psiz(i,kwinf(knew)) * ( gdepw_1d(kwsup(knew),1)+znew(knew))   &
      &                          + psiz(i,kwsup(knew)) * (-znew(knew)-gdepw_1d(kwinf(knew),1)) ) &
      &                        / ( gdepw_1d(kwsup(knew),1) - gdepw_1d(kwinf(knew),1) )
    enddo
  enddo

ENDDO

!---------------------------------------
! Writing new netcdf file :

write(*,*) 'Writing ', TRIM(file_out)

status = NF90_CREATE(TRIM(file_out),NF90_NOCLOBBER,fidM)
call erreur(status,.TRUE.,'create')
 
status = NF90_DEF_DIM(fidM,"nTime",NF90_UNLIMITED,dimID_time_counter)
call erreur(status,.TRUE.,"def_dimID_time_counter")
status = NF90_DEF_DIM(fidM,"nz",mznew,dimID_z)
call erreur(status,.TRUE.,"def_dimID_z")
status = NF90_DEF_DIM(fidM,"nx",mx-2,dimID_x)
call erreur(status,.TRUE.,"def_dimID_x")
status = NF90_DEF_DIM(fidM,"ny",my-2,dimID_y)
call erreur(status,.TRUE.,"def_dimID_y")

status = NF90_DEF_VAR(fidM,"time",NF90_FLOAT,(/dimID_time_counter/),time_ID)
call erreur(status,.TRUE.,"def_var_time_ID")
status = NF90_DEF_VAR(fidM,"uBoundaryLayer",NF90_FLOAT,(/dimID_x,dimID_y,dimID_time_counter/),utbl_ID)
call erreur(status,.TRUE.,"def_var_utbl_ID")
status = NF90_DEF_VAR(fidM,"vBoundaryLayer",NF90_FLOAT,(/dimID_x,dimID_y,dimID_time_counter/),vtbl_ID)
call erreur(status,.TRUE.,"def_var_vtbl_ID")
status = NF90_DEF_VAR(fidM,"barotropicStreamfunction",NF90_FLOAT,(/dimID_x,dimID_y,dimID_time_counter/),BSF_ID)
call erreur(status,.TRUE.,"def_var_BSF_ID")
status = NF90_DEF_VAR(fidM,"overturningStreamfunction",NF90_FLOAT,(/dimID_x,dimID_z,dimID_time_counter/),MOC_ID)
call erreur(status,.TRUE.,"def_var_MOC_ID")
status = NF90_DEF_VAR(fidM,"z",NF90_FLOAT,(/dimID_z/),z_ID)
call erreur(status,.TRUE.,"def_var_z_ID")
status = NF90_DEF_VAR(fidM,"x",NF90_FLOAT,(/dimID_x/),x_ID)
call erreur(status,.TRUE.,"def_var_x_ID")
status = NF90_DEF_VAR(fidM,"y",NF90_FLOAT,(/dimID_y/),y_ID)
call erreur(status,.TRUE.,"def_var_y_ID")

status = NF90_PUT_ATT(fidM,utbl_ID,"units","m/s")
call erreur(status,.TRUE.,"put_att_utbl_ID")
status = NF90_PUT_ATT(fidM,utbl_ID,"description","x-velocity in the boundary layer used to compute u*")
call erreur(status,.TRUE.,"put_att_utbl_ID")
status = NF90_PUT_ATT(fidM,utbl_ID,"_FillValue",NF90_FILL_FLOAT)
call erreur(status,.TRUE.,"put_att_utbl_ID")
!-
status = NF90_PUT_ATT(fidM,vtbl_ID,"units","m/s")
call erreur(status,.TRUE.,"put_att_vtbl_ID")
status = NF90_PUT_ATT(fidM,vtbl_ID,"description","y-velocity in the boundary layer used to compute u*")
call erreur(status,.TRUE.,"put_att_vtbl_ID")
status = NF90_PUT_ATT(fidM,vtbl_ID,"_FillValue",NF90_FILL_FLOAT)
call erreur(status,.TRUE.,"put_att_vtbl_ID")
!-
status = NF90_PUT_ATT(fidM,BSF_ID,"units","m^3/s")
call erreur(status,.TRUE.,"put_att_BSF_ID")
status = NF90_PUT_ATT(fidM,BSF_ID,"description","barotropic streamfunction")
call erreur(status,.TRUE.,"put_att_BSF_ID")
status = NF90_PUT_ATT(fidM,BSF_ID,"_FillValue",NF90_FILL_FLOAT)
call erreur(status,.TRUE.,"put_att_BSF_ID")
!-
status = NF90_PUT_ATT(fidM,MOC_ID,"units","m^3/s")
call erreur(status,.TRUE.,"put_att_MOC_ID")
status = NF90_PUT_ATT(fidM,MOC_ID,"description","overturning (meridional) streamfunction")
call erreur(status,.TRUE.,"put_att_MOC_ID")
status = NF90_PUT_ATT(fidM,MOC_ID,"_FillValue",NF90_FILL_FLOAT)
call erreur(status,.TRUE.,"put_att_MOC_ID")
!-
status = NF90_PUT_ATT(fidM,time_ID,"units","seconds since 0000-01-01 00:00:00")
call erreur(status,.TRUE.,"put_att_time_ID")
status = NF90_PUT_ATT(fidM,time_ID,"calendar","noleap")
call erreur(status,.TRUE.,"put_att_time_ID")
status = NF90_PUT_ATT(fidM,time_ID,"description","time at the beginning of the period over which monthly means are calculated")
call erreur(status,.TRUE.,"put_att_time_ID")
!-
status = NF90_PUT_ATT(fidM,z_ID,"units","m")
call erreur(status,.TRUE.,"put_att_z_ID")
status = NF90_PUT_ATT(fidM,z_ID,"description","Vertical axis")
call erreur(status,.TRUE.,"put_att_z_ID")
!-
status = NF90_PUT_ATT(fidM,x_ID,"units","m")
call erreur(status,.TRUE.,"put_att_x_ID")
status = NF90_PUT_ATT(fidM,x_ID,"description","X axis")
call erreur(status,.TRUE.,"put_att_x_ID")
!-
status = NF90_PUT_ATT(fidM,y_ID,"units","m")
call erreur(status,.TRUE.,"put_att_y_ID")
status = NF90_PUT_ATT(fidM,y_ID,"description","Y axis")
call erreur(status,.TRUE.,"put_att_y_ID")

status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"experiment","ISOMIP+ COM_EEEE")
call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"history","Run and post-processed by Nicolas Jourdain")
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
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"momentum_advection","flux form - 3rd order UBS")
call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"tracer_advection","3rd order UBS and 2nd order FCT on the vertical")
call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"lateral_diffusion_for_momentum","Laplacian, along geopotential (horizontal)")
call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"lateral_diffusion_for_tracers","Laplacian, along iso-neutral")
call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"vertical_diffusion","Constant eddy diffusivity")
call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"convection","enhanced vertical diffusivity and viscosity (0.1 m2/s)")
call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"top_mixed_layer","T,S,u are averaged over the top 20m")
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
 
status = NF90_PUT_VAR(fidM,time_ID,time)
call erreur(status,.TRUE.,"var_time_ID")
status = NF90_PUT_VAR(fidM,utbl_ID,ISOMIP_uTBL)
call erreur(status,.TRUE.,"var_utbl_ID")
status = NF90_PUT_VAR(fidM,vtbl_ID,ISOMIP_vTBL)
call erreur(status,.TRUE.,"var_vtbl_ID")
status = NF90_PUT_VAR(fidM,BSF_ID,ISOMIP_BSF)
call erreur(status,.TRUE.,"var_BSF_ID")
status = NF90_PUT_VAR(fidM,MOC_ID,ISOMIP_MOC)
call erreur(status,.TRUE.,"var_MOC_ID")
status = NF90_PUT_VAR(fidM,z_ID,znew)
call erreur(status,.TRUE.,"var_z_ID")
status = NF90_PUT_VAR(fidM,x_ID,glamt(2:mx-1,1,1)*1.e3)
call erreur(status,.TRUE.,"var_x_ID")
status = NF90_PUT_VAR(fidM,y_ID,gphit(1,2:my-1,1)*1.e3)
call erreur(status,.TRUE.,"var_y_ID")
 
status = NF90_CLOSE(fidM)                    
call erreur(status,.TRUE.,"final")         

end program modif



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
