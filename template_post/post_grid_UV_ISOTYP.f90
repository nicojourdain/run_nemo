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
&          time_ID, utbl_ID, uoce_e3u_ID, vtbl_ID, voce_e3v_ID, depthu_bounds_ID, iiso, jiso,               &
&          depthu_ID, depthv_bounds_ID, depthv_ID, fidM, mxiso, myiso, mziso

INTEGER :: fidMSH, dimID_t, dimID_z, dimID_y, dimID_x, mt, mz, vmask_ID, umask_ID, tmask_ID, misf_ID, mbathy_ID, i, j, k, l,      &
&          e3w_1d_ID, e3t_1d_ID, gdepw_1d_ID, gdept_1d_ID, e3w_0_ID, e3v_0_ID, e3u_0_ID, e3t_0_ID, ff_ID, e2v_ID, e2u_ID, e2t_ID, &
&          e1v_ID, e1u_ID, e1t_ID, time_counter_ID, gdepw_0_ID, gdepv_ID, gdepu_ID, gdept_0_ID, isfdraft_ID, gphiv_ID, gphiu_ID,  &
&          gphit_ID, glamv_ID, glamu_ID, glamt_ID, glamf_ID, gphif_ID, x_ID, y_ID, z_ID, BSF_ID, MOC_ID, knew

CHARACTER(LEN=120) :: file_MSH, file_in_U, file_in_V, file_out                     

INTEGER*1,ALLOCATABLE,DIMENSION(:,:,:,:) :: vmask, umask, tmask

INTEGER*2,ALLOCATABLE,DIMENSION(:,:,:) :: misf, mbathy

INTEGER*2,ALLOCATABLE,DIMENSION(:,:) :: isfmskORI

INTEGER*2,ALLOCATABLE,DIMENSION(:) :: kwsup, kwinf, iinfT, iinfU, iinfV, iinfF, isupT, isupU, isupV, isupF, &
&                                     jinfT, jinfU, jinfV, jinfF, jsupT, jsupU, jsupV, jsupF

REAL*4,ALLOCATABLE,DIMENSION(:) :: depthu, depthv, znew, time, wxinfT, wxinfU, wxinfV, wxinfF, wxsupT, wxsupU, wxsupV, wxsupF, &
&                                  xiso, yiso, wyinfT, wyinfU, wyinfV, wyinfF, wysupT, wysupU, wysupV, wysupF

REAL*4,ALLOCATABLE,DIMENSION(:,:) :: depthu_bounds, depthv_bounds, trpv, trpz, psif, psiz

REAL*8,ALLOCATABLE,DIMENSION(:,:) :: e3w_1d, e3t_1d, gdepw_1d, gdept_1d

REAL*4,ALLOCATABLE,DIMENSION(:,:,:) :: utbl, vtbl, isfdraft, gphiv, gphiu, gphit, glamv, glamu, glamt,   &
&                                      e1u, e2u, e1v, e2v, e1t, e2t, ISOTYP_BSF, ISOTYP_MOC, ISOTYP_uTBL,&
&                                      ISOTYP_vTBL, tmp_MOC, gphif, glamf

REAL*4,ALLOCATABLE,DIMENSION(:,:,:,:) :: uoce_e3u, voce_e3v, gdepw_0, gdepv, gdepu, gdept_0

REAL*8,ALLOCATABLE,DIMENSION(:,:,:,:) :: e3u_0, e3v_0  !, e3w_0, e3t_0

INTEGER*4, DIMENSION(12) :: daym1

INTEGER*4 :: an, month0

REAL*4 :: Gt, eps

file_MSH   = '../../input/nemo_ISOTYP/mesh_mask_ISOTYP_EEEE.nc'                                
file_in_U  = 'MONTHLY/ISOTYP-CCCC_monthly_YYYY_grid_U.nc'                                
file_in_V  = 'MONTHLY/ISOTYP-CCCC_monthly_YYYY_grid_V.nc'
file_out   = 'MONTHLY/ISOTYP-CCCC_monthly_YYYY_tmpUV.nc'

eps=1.e-9

Gt = GGGG

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
ALLOCATE(  gphif(mx,my,mt)  ) 
ALLOCATE(  glamv(mx,my,mt)  ) 
ALLOCATE(  glamu(mx,my,mt)  ) 
ALLOCATE(  glamt(mx,my,mt)  ) 
ALLOCATE(  glamf(mx,my,mt)  )

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
status = NF90_INQ_VARID(fidMSH,"gphif",gphif_ID)
call erreur(status,.TRUE.,"inq_gphif_ID")
status = NF90_INQ_VARID(fidMSH,"glamv",glamv_ID)
call erreur(status,.TRUE.,"inq_glamv_ID")
status = NF90_INQ_VARID(fidMSH,"glamu",glamu_ID)
call erreur(status,.TRUE.,"inq_glamu_ID")
status = NF90_INQ_VARID(fidMSH,"glamt",glamt_ID)
call erreur(status,.TRUE.,"inq_glamt_ID")
status = NF90_INQ_VARID(fidMSH,"glamf",glamf_ID)
call erreur(status,.TRUE.,"inq_glamf_ID")

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
status = NF90_GET_VAR(fidMSH,gphif_ID,gphif)
call erreur(status,.TRUE.,"getvar_gphif")
status = NF90_GET_VAR(fidMSH,glamv_ID,glamv)
call erreur(status,.TRUE.,"getvar_glamv")
status = NF90_GET_VAR(fidMSH,glamu_ID,glamu)
call erreur(status,.TRUE.,"getvar_glamu")
status = NF90_GET_VAR(fidMSH,glamt_ID,glamt)
call erreur(status,.TRUE.,"getvar_glamt")
status = NF90_GET_VAR(fidMSH,glamf_ID,glamf)
call erreur(status,.TRUE.,"getvar_glamf")

status = NF90_CLOSE(fidMSH)                      
call erreur(status,.TRUE.,"fin_lecture")     

!------------------------------------------------------
! refined output vertical grid :
mziso=144
ALLOCATE( znew(mziso), kwsup(mziso), kwinf(mziso) )
!write(*,*) 'New vertical grid for the outputs :'
do knew=1,mziso
  znew(knew) = 2.5 - knew * 5.0
  do k=1,mz-1
    if ( -znew(knew) .ge. gdepw_1d(k,1) .and. -znew(knew) .lt. gdepw_1d(k+1,1) ) then
      kwsup(knew) = k+1
      kwinf(knew) = k
    endif
  enddo
  !write(*,109) -znew(knew), gdepw_1d(kwinf(knew),1), kwinf(knew), gdepw_1d(kwsup(knew),1), kwsup(knew)
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
ALLOCATE( tmp_MOC(mx,mziso,mtime_counter) )

!-- Parameters for interpolation :
mxiso=240
myiso=40
!--
ALLOCATE( ISOTYP_BSF (mxiso,myiso  ,mtime_counter) )
ALLOCATE( ISOTYP_MOC (mxiso,  mziso,mtime_counter) )
ALLOCATE( ISOTYP_uTBL(mxiso,myiso  ,mtime_counter) )
ALLOCATE( ISOTYP_vTBL(mxiso,myiso  ,mtime_counter) )
ALLOCATE( xiso(mxiso), yiso(myiso) )
!--
ALLOCATE( iinfT(mxiso), jinfT(myiso), isupT(mxiso), jsupT(myiso) )
ALLOCATE( wxinfT(mxiso), wyinfT(myiso), wxsupT(mxiso), wysupT(myiso) )
ALLOCATE( iinfU(mxiso), jinfU(myiso), isupU(mxiso), jsupU(myiso) )
ALLOCATE( wxinfU(mxiso), wyinfU(myiso), wxsupU(mxiso), wysupU(myiso) )
ALLOCATE( iinfV(mxiso), jinfV(myiso), isupV(mxiso), jsupV(myiso) )
ALLOCATE( wxinfV(mxiso), wyinfV(myiso), wxsupV(mxiso), wysupV(myiso) )
ALLOCATE( iinfF(mxiso), jinfF(myiso), isupF(mxiso), jsupF(myiso) )
ALLOCATE( wxinfF(mxiso), wyinfF(myiso), wxsupF(mxiso), wysupF(myiso) )
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
    if ( xiso(iiso) .ge. glamu(i,1,1)*1.e3 .and. xiso(iiso) .lt. glamu(i+1,1,1)*1.e3 ) then
      iinfU(iiso)  = i
      isupU(iiso)  = i+1
      wxinfU(iiso) = ( glamu(i+1,1,1)*1.e3 - xiso(iiso) ) / ( glamu(i+1,1,1)*1.e3 - glamu(i,1,1)*1.e3 )
      wxsupU(iiso) = ( xiso(iiso) - glamu(i,1,1)*1.e3   ) / ( glamu(i+1,1,1)*1.e3 - glamu(i,1,1)*1.e3 )
    endif
    if ( xiso(iiso) .ge. glamv(i,1,1)*1.e3 .and. xiso(iiso) .lt. glamv(i+1,1,1)*1.e3 ) then
      iinfV(iiso)  = i
      isupV(iiso)  = i+1
      wxinfV(iiso) = ( glamv(i+1,1,1)*1.e3 - xiso(iiso) ) / ( glamv(i+1,1,1)*1.e3 - glamv(i,1,1)*1.e3 )
      wxsupV(iiso) = ( xiso(iiso) - glamv(i,1,1)*1.e3   ) / ( glamv(i+1,1,1)*1.e3 - glamv(i,1,1)*1.e3 )
    endif
    if ( xiso(iiso) .ge. glamf(i,1,1)*1.e3 .and. xiso(iiso) .lt. glamf(i+1,1,1)*1.e3 ) then
      iinfF(iiso)  = i
      isupF(iiso)  = i+1
      wxinfF(iiso) = ( glamf(i+1,1,1)*1.e3 - xiso(iiso) ) / ( glamf(i+1,1,1)*1.e3 - glamf(i,1,1)*1.e3 )
      wxsupF(iiso) = ( xiso(iiso) - glamf(i,1,1)*1.e3   ) / ( glamf(i+1,1,1)*1.e3 - glamf(i,1,1)*1.e3 )
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
  if ( xiso(iiso) .lt. glamu(1,1,1)*1.e3 ) then
    iinfU(iiso)  = 1
    isupU(iiso)  = 1
    wxinfU(iiso) = 1.000000
    wxsupU(iiso) = 0.000000
  elseif ( xiso(iiso) .ge. glamu(mx,1,1)*1.e3 ) then
    iinfU(iiso)  = mx
    isupU(iiso)  = mx
    wxinfU(iiso) = 1.000000
    wxsupU(iiso) = 0.000000
  endif
  if ( xiso(iiso) .lt. glamv(1,1,1)*1.e3 ) then
    iinfV(iiso)  = 1
    isupV(iiso)  = 1
    wxinfV(iiso) = 1.000000
    wxsupV(iiso) = 0.000000
  elseif ( xiso(iiso) .ge. glamv(mx,1,1)*1.e3 ) then
    iinfV(iiso)  = mx
    isupV(iiso)  = mx
    wxinfV(iiso) = 1.000000
    wxsupV(iiso) = 0.000000
  endif
  if ( xiso(iiso) .lt. glamf(1,1,1)*1.e3 ) then
    iinfF(iiso)  = 1
    isupF(iiso)  = 1
    wxinfF(iiso) = 1.000000
    wxsupF(iiso) = 0.000000
  elseif ( xiso(iiso) .ge. glamf(mx,1,1)*1.e3 ) then
    iinfF(iiso)  = mx
    isupF(iiso)  = mx
    wxinfF(iiso) = 1.000000
    wxsupF(iiso) = 0.000000
  endif
enddo
!--
do jiso=1,myiso
  yiso(jiso) = -1000.000000 + 2000.000000 * jiso
  do j=1,my-1
    if ( yiso(jiso) .ge. gphit(1,j,1)*1.e3 .and. yiso(jiso) .lt. gphit(1,j+1,1)*1.e3 ) then
      jinfT(jiso)  = j
      jsupT(jiso)  = j+1
      wyinfT(jiso) = ( gphit(1,j+1,1)*1.e3 - yiso(jiso) ) / ( gphit(1,j+1,1)*1.e3 - gphit(1,j,1)*1.e3 )
      wysupT(jiso) = ( yiso(jiso) - gphit(1,j,1)*1.e3   ) / ( gphit(1,j+1,1)*1.e3 - gphit(1,j,1)*1.e3 )
    endif
    if ( yiso(jiso) .ge. gphiu(1,j,1)*1.e3 .and. yiso(jiso) .lt. gphiu(1,j+1,1)*1.e3 ) then
      jinfU(jiso)  = j
      jsupU(jiso)  = j+1
      wyinfU(jiso) = ( gphiu(1,j+1,1)*1.e3 - yiso(jiso) ) / ( gphiu(1,j+1,1)*1.e3 - gphiu(1,j,1)*1.e3 )
      wysupU(jiso) = ( yiso(jiso) - gphiu(1,j,1)*1.e3   ) / ( gphiu(1,j+1,1)*1.e3 - gphiu(1,j,1)*1.e3 )
    endif
    if ( yiso(jiso) .ge. gphiv(1,j,1)*1.e3 .and. yiso(jiso) .lt. gphiv(1,j+1,1)*1.e3 ) then
      jinfV(jiso)  = j
      jsupV(jiso)  = j+1
      wyinfV(jiso) = ( gphiv(1,j+1,1)*1.e3 - yiso(jiso) ) / ( gphiv(1,j+1,1)*1.e3 - gphiv(1,j,1)*1.e3 )
      wysupV(jiso) = ( yiso(jiso) - gphiv(1,j,1)*1.e3   ) / ( gphiv(1,j+1,1)*1.e3 - gphiv(1,j,1)*1.e3 )
    endif
    if ( yiso(jiso) .ge. gphif(1,j,1)*1.e3 .and. yiso(jiso) .lt. gphif(1,j+1,1)*1.e3 ) then
      jinfF(jiso)  = j
      jsupF(jiso)  = j+1
      wyinfF(jiso) = ( gphif(1,j+1,1)*1.e3 - yiso(jiso) ) / ( gphif(1,j+1,1)*1.e3 - gphif(1,j,1)*1.e3 )
      wysupF(jiso) = ( yiso(jiso) - gphif(1,j,1)*1.e3   ) / ( gphif(1,j+1,1)*1.e3 - gphif(1,j,1)*1.e3 )
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
  if ( yiso(jiso) .lt. gphiu(1,1,1)*1.e3 ) then
    jinfU(jiso)  = 1
    jsupU(jiso)  = 1
    wyinfU(jiso) = 1.000000
    wysupU(jiso) = 0.000000
  elseif ( yiso(jiso) .ge. gphiu(1,my,1)*1.e3 ) then
    jinfU(jiso)  = my
    jsupU(jiso)  = my
    wyinfU(jiso) = 1.000000
    wysupU(jiso) = 0.000000
  endif
  if ( yiso(jiso) .lt. gphiv(1,1,1)*1.e3 ) then
    jinfV(jiso)  = 1
    jsupV(jiso)  = 1
    wyinfV(jiso) = 1.000000
    wysupV(jiso) = 0.000000
  elseif ( yiso(jiso) .ge. gphiv(1,my,1)*1.e3 ) then
    jinfV(jiso)  = my
    jsupV(jiso)  = my
    wyinfV(jiso) = 1.000000
    wysupV(jiso) = 0.000000
  endif
  if ( yiso(jiso) .lt. gphif(1,1,1)*1.e3 ) then
    jinfF(jiso)  = 1
    jsupF(jiso)  = 1
    wyinfF(jiso) = 1.000000
    wysupF(jiso) = 0.000000
  elseif ( yiso(jiso) .ge. gphif(1,my,1)*1.e3 ) then
    jinfF(jiso)  = my
    jsupF(jiso)  = my
    wyinfF(jiso) = 1.000000
    wysupF(jiso) = 0.000000
  endif
  !write(*,*) '## ', jiso, yiso(jiso), jinfU(jiso), jsupU(jiso)
  !write(*,*) '   ', gphiu(1,1,1)*1.e3, gphiu(1,my,1)*1.e3
enddo

!- ice-shelf mask :
ALLOCATE( isfmskORI(mx,my) )
where ( misf(:,:,1) .le. 1 )
  isfmskORI(:,:) = 0
elsewhere
  isfmskORI(:,:) = 1
endwhere

!======================================
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
  do k=1,mz
    write(*,*) 'l, k, trpz(1:2,k) = ', l, k, trpz(1:2,k)
  enddo
  ! BSF calculated by integrating transport along X
  do i=2,mx
    psif(i,:) = psif(i-1,:) - trpv(i,:)  ! BSF at F-points
  enddo
  ! MOC calculated by integrating transport along Z
  do k=mz-1,1,-1
    write(*,*) 'l, k, psiz(1,k+1), trpz(1,k) = ', l, k, psiz(1,k+1), trpz(1,k)
    psiz(:,k) = psiz(:,k+1) - trpz(:,k) ! MOC at W-points
    write(*,*) 'l, k, psiz(1,k) = ', l, k, psiz(1,k)
  enddo

  ! MOC (interpolation from W-point to ISOTYP 5-m vertical grid )
  do i=1,mx  
    do knew=1,mziso
      tmp_MOC(i,knew,l) = (   psiz(i,kwinf(knew)) * ( gdepw_1d(kwsup(knew),1)+znew(knew))   &
      &                     + psiz(i,kwsup(knew)) * (-znew(knew)-gdepw_1d(kwinf(knew),1)) ) &
      &                 / ( gdepw_1d(kwsup(knew),1) - gdepw_1d(kwinf(knew),1) )
    enddo
  enddo
  write(*,*) '   znew(30), tmp_MOC(1:2,30) = ', znew(30), tmp_MOC(1:2,30,l)

  !- interpolate BSF and MOC on ISOMIP's T-points 
  !  (NB: for regular grid in X,Y only):
  !  -> directly on ISOMIP's smaller grid ( mxiso * myiso )
  do iiso=1,mxiso

    ISOTYP_MOC(iiso,:,l) =   (   tmp_MOC(iinfT(iiso),:,l) * wxinfT(iiso)   &
    &                          + tmp_MOC(isupT(iiso),:,l) * wxsupT(iiso) ) &
    &                      / (                              wxinfT(iiso)   &
    &                          +                            wxsupT(iiso) )

    do jiso=1,myiso
  
      ! BSF on ISOMIP grid :
      ISOTYP_BSF(iiso,jiso,l) =   (   psif(iinfF(iiso),jinfF(jiso)) * wxinfF(iiso) * wyinfF(jiso)   &
      &                             + psif(iinfF(iiso),jsupF(jiso)) * wxinfF(iiso) * wysupF(jiso)   &
      &                             + psif(isupF(iiso),jinfF(jiso)) * wxsupF(iiso) * wyinfF(jiso)   &
      &                             + psif(isupF(iiso),jsupF(jiso)) * wxsupF(iiso) * wysupF(jiso) ) &
      &                         / (                                   wxinfF(iiso) * wyinfF(jiso)   &
      &                             +                                 wxinfF(iiso) * wysupF(jiso)   &
      &                             +                                 wxsupF(iiso) * wyinfF(jiso)   &
      &                             +                                 wxsupF(iiso) * wysupF(jiso) )
      !write(*,*)  iinfU(iiso),jinfU(jiso),isupU(iiso),jsupU(jiso)
      !write(*,*)  utbl(iinfU(iiso),jinfU(jiso),l) , isfmskORI(iinfU(iiso),jinfU(jiso)) , wxinfU(iiso) , wyinfU(jiso)
      !write(*,*)  utbl(iinfU(iiso),jsupU(jiso),l) , isfmskORI(iinfU(iiso),jsupU(jiso)) , wxinfU(iiso) , wysupU(jiso)
      !write(*,*)  utbl(isupU(iiso),jinfU(jiso),l) , isfmskORI(isupU(iiso),jinfU(jiso)) , wxsupU(iiso) , wyinfU(jiso)
      !write(*,*)  utbl(isupU(iiso),jsupU(jiso),l) , isfmskORI(isupU(iiso),jsupU(jiso)) , wxsupU(iiso) , wysupU(jiso)
  
      ! uTBL on ISOMIP grid :
      if (   isfmskORI( iinfU(iiso), jinfU(jiso) )            &
         & + isfmskORI( iinfU(iiso), jsupU(jiso) )            &
         & + isfmskORI( isupU(iiso), jinfU(jiso) )            &
         & + isfmskORI( isupU(iiso), jsupU(jiso) )  .gt. 0.1  ) then
          ! 
          ISOTYP_uTBL(iiso,jiso,l) =  (   utbl(iinfU(iiso),jinfU(jiso),l) * isfmskORI(iinfU(iiso),jinfU(jiso)) * wxinfU(iiso) * wyinfU(jiso)   &
          &                             + utbl(iinfU(iiso),jsupU(jiso),l) * isfmskORI(iinfU(iiso),jsupU(jiso)) * wxinfU(iiso) * wysupU(jiso)   &
          &                             + utbl(isupU(iiso),jinfU(jiso),l) * isfmskORI(isupU(iiso),jinfU(jiso)) * wxsupU(iiso) * wyinfU(jiso)   &
          &                             + utbl(isupU(iiso),jsupU(jiso),l) * isfmskORI(isupU(iiso),jsupU(jiso)) * wxsupU(iiso) * wysupU(jiso) ) &
          &                         / (                                     isfmskORI(iinfU(iiso),jinfU(jiso)) * wxinfU(iiso) * wyinfU(jiso)   &
          &                             +                                   isfmskORI(iinfU(iiso),jsupU(jiso)) * wxinfU(iiso) * wysupU(jiso)   &
          &                             +                                   isfmskORI(isupU(iiso),jinfU(jiso)) * wxsupU(iiso) * wyinfU(jiso)   &
          &                             +                                   isfmskORI(isupU(iiso),jsupU(jiso)) * wxsupU(iiso) * wysupU(jiso) )
          !
      else
          ISOTYP_uTBL(iiso,jiso,l) = NF90_FILL_FLOAT
      endif
     
      ! vTBL on ISOMIP grid :
      if (   isfmskORI( iinfV(iiso), jinfV(jiso) )            &
         & + isfmskORI( iinfV(iiso), jsupV(jiso) )            &
         & + isfmskORI( isupV(iiso), jinfV(jiso) )            &
         & + isfmskORI( isupV(iiso), jsupV(jiso) )  .gt. 0.1  ) then
          ! 
          ISOTYP_vTBL(iiso,jiso,l) =  (   vtbl(iinfV(iiso),jinfV(jiso),l) * isfmskORI(iinfV(iiso),jinfV(jiso)) * wxinfV(iiso) * wyinfV(jiso)   &
          &                             + vtbl(iinfV(iiso),jsupV(jiso),l) * isfmskORI(iinfV(iiso),jsupV(jiso)) * wxinfV(iiso) * wysupV(jiso)   &
          &                             + vtbl(isupV(iiso),jinfV(jiso),l) * isfmskORI(isupV(iiso),jinfV(jiso)) * wxsupV(iiso) * wyinfV(jiso)   &
          &                             + vtbl(isupV(iiso),jsupV(jiso),l) * isfmskORI(isupV(iiso),jsupV(jiso)) * wxsupV(iiso) * wysupV(jiso) ) &
          &                         / (                                     isfmskORI(iinfV(iiso),jinfV(jiso)) * wxinfV(iiso) * wyinfV(jiso)   &
          &                             +                                   isfmskORI(iinfV(iiso),jsupV(jiso)) * wxinfV(iiso) * wysupV(jiso)   &
          &                             +                                   isfmskORI(isupV(iiso),jinfV(jiso)) * wxsupV(iiso) * wyinfV(jiso)   &
          &                             +                                   isfmskORI(isupV(iiso),jsupV(jiso)) * wxsupV(iiso) * wysupV(jiso) )
          !
      else
          ISOTYP_vTBL(iiso,jiso,l) = NF90_FILL_FLOAT
      endif
   
    enddo

  enddo

  write(*,*) '  ISOTYP_MOC(1:2,30,l) = ', ISOTYP_MOC(1:2,30,l)

ENDDO

!---------------------------------------
! Writing new netcdf file :

write(*,*) 'Writing ', TRIM(file_out)

status = NF90_CREATE(TRIM(file_out),NF90_NOCLOBBER,fidM)
call erreur(status,.TRUE.,'create')
 
status = NF90_DEF_DIM(fidM,"nTime",NF90_UNLIMITED,dimID_time_counter)
call erreur(status,.TRUE.,"def_dimID_time_counter")
status = NF90_DEF_DIM(fidM,"nz",mziso,dimID_z)
call erreur(status,.TRUE.,"def_dimID_z")
status = NF90_DEF_DIM(fidM,"nx",mxiso,dimID_x)
call erreur(status,.TRUE.,"def_dimID_x")
status = NF90_DEF_DIM(fidM,"ny",myiso,dimID_y)
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
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"convection","Enhanced vertical diffusivity and viscosity (100 m2/s)")
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
 
status = NF90_PUT_VAR(fidM,time_ID,time)
call erreur(status,.TRUE.,"var_time_ID")
status = NF90_PUT_VAR(fidM,utbl_ID,ISOTYP_uTBL)
call erreur(status,.TRUE.,"var_utbl_ID")
status = NF90_PUT_VAR(fidM,vtbl_ID,ISOTYP_vTBL)
call erreur(status,.TRUE.,"var_vtbl_ID")
status = NF90_PUT_VAR(fidM,BSF_ID,ISOTYP_BSF)
call erreur(status,.TRUE.,"var_BSF_ID")
status = NF90_PUT_VAR(fidM,MOC_ID,ISOTYP_MOC)
call erreur(status,.TRUE.,"var_MOC_ID")
status = NF90_PUT_VAR(fidM,z_ID,znew)
call erreur(status,.TRUE.,"var_z_ID")
status = NF90_PUT_VAR(fidM,x_ID,xiso)
call erreur(status,.TRUE.,"var_x_ID")
status = NF90_PUT_VAR(fidM,y_ID,yiso)
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
