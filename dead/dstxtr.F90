! $Id$

! Purpose: External data needed by dust model

! Usage:
! use dstxtr ! [mdl] External data

! dst.h needed for DST_NBR
#include <dst.h> /* Dust preprocessor tokens */

module dstxtr ! [mdl] External data
  use shr_kind_mod,only:r8=>shr_kind_r8 ! [mdl] Precision r8, i8, ...
  implicit none

  public::xtr_dat_get ! [sbr] Read external forcing data
  public::xtr_dat_get_EuF ! [sbr] Read external forcing data from Eureka Flats
  public::xtr_dat_set ! [sbr] Set time-varying forcing data, real types
  private::xtr_dat_set_flt_1D ! [sbr] Set time-varying forcing data, real types, 1D arrays
  private::xtr_dat_set_flt_scl ! [sbr] Set time-varying forcing data, real types, scalars
  private::xtr_dat_set_int_1D ! [sbr] Set time-varying forcing data, integer types, 1D arrays
  
  ! Overloaded data set routines
  interface xtr_dat_set
     module procedure xtr_dat_set_flt_scl,xtr_dat_set_flt_1D,xtr_dat_set_int_1D
  end interface ! xtr_dat_set()

  integer,public,parameter::xtr_idx_asp_rat_lps=1 ! [idx] Index of asp_rat_lps
  integer,public,parameter::xtr_idx_doy=2 ! [idx] Index of doy
  integer,public,parameter::xtr_idx_hgt_mdp=3 ! [idx] Index of hgt_mdp
  integer,public,parameter::xtr_idx_oro=4 ! [idx] Index of oro
  integer,public,parameter::xtr_idx_prs_mdp=5 ! [idx] Index of prs_mdp
  integer,public,parameter::xtr_idx_prs_ntf=6 ! [idx] Index of prs_ntf
  integer,public,parameter::xtr_idx_q_H2O_vpr=7 ! [idx] Index of q_H2O_vpr
  integer,public,parameter::xtr_idx_sfc_typ=8 ! [idx] Index of sfc_typ
  integer,public,parameter::xtr_idx_tpt_gnd=9 ! [idx] Index of tpt_gnd
  integer,public,parameter::xtr_idx_tpt_ice=10 ! [idx] Index of tpt_ice
  integer,public,parameter::xtr_idx_tpt_mdp=11 ! [idx] Index of tpt_mdp
  integer,public,parameter::xtr_idx_tpt_soi=12 ! [idx] Index of tpt_soi
  integer,public,parameter::xtr_idx_tpt_sst=13 ! [idx] Index of tpt_sst
  integer,public,parameter::xtr_idx_vai_dst=14 ! [idx] Index of vai_dst
  integer,public,parameter::xtr_idx_vwc_sfc=15 ! [idx] Index of vwc_sfc
  integer,public,parameter::xtr_idx_wnd_mrd_mdp=16 ! [idx] Index of wnd_mrd_mdp
  integer,public,parameter::xtr_idx_wnd_znl_mdp=17 ! [idx] Index of wnd_znl_mdp
  integer,public,parameter::xtr_idx_end=17 ! [idx] Last xtr_dat_idx

contains
  
  subroutine xtr_dat_get( & ! [sbr] Read external forcing data
       fl_in, & ! I [sng] netCDF file with external forcing data
       time_nbr, & ! I [nbr] Number of timesteps to simulate
       plond, & ! I [nbr] Number of longitudes
       xtr_dat_nbr, & ! I [nbr] Number of external forcing data variables
       xtr_dat_flg, & ! O [flg] External forcing data exists for this variable
       xtr_dat_frc) ! O [frc] External forcing data (time_nbr,plond,xtr_dat_nbr)
    ! Purpose: Get time-varying zonal or meridional wind speed
    use shr_kind_mod,only:r8=>shr_kind_r8 ! [mdl] Precision r8, i8, ...
    use netcdf    ! [mdl] netCDF interface
    use nf90_utl  ! [mdl] netCDF utilities
    use utl_mdl   ! [mdl] Utility functions (date_time_get,mnt_chk...)
    implicit none
    ! Parameters
    character(len=*),parameter::sbr_nm="xtr_dat_get" ! [sng] Subroutine name
    ! Input
    character(len=*),intent(in)::fl_in ! [sng] netCDF file with external forcing data
    integer,intent(in)::time_nbr ! [nbr] Number of timesteps to simulate
    integer,intent(in)::plond ! [nbr] Number of longitudes
    integer,intent(in)::xtr_dat_nbr ! [nbr] Number of external forcing data variables
    ! Output
    logical,intent(out)::xtr_dat_flg(xtr_dat_nbr) ! O [flg] External forcing data exists for this variable
    real(r8),intent(out)::xtr_dat_frc(time_nbr,plond,xtr_dat_nbr) ! O [frc] External forcing data
    ! Local
    ! File metadata and dimension IDs
    integer nc_id             ! [id] File ID
    integer rcd               ! [rcd] Return success code
    integer time_dmn_id       ! [enm] Dimension ID for time
    integer time_id           ! [enm] Variable ID
    integer time_nbr_in       ! [nbr] Dimension size
    integer xtr_dat_frc_id(xtr_idx_end)
    
    integer idx
    real(r8) time(time_nbr)   ! [day] Time coordinate (day of year)
    
    idx=0 ! CEWI

    ! Initialize external forcing data to missing value
    xtr_dat_flg(:)=.false.
    xtr_dat_frc(:,:,:)=1.0e36

    ! Open netCDF file
    rcd=nf90_wrp_open(fl_in,nf90_nowrite,nc_id,sbr_nm=sbr_nm)
    ! Get dimension IDs
    rcd=nf90_wrp_inq_dimid(nc_id,"time",time_dmn_id)
    rcd=nf90_wrp(nf90_inquire_dimension(nc_id,time_dmn_id,len=time_nbr_in),sbr_nm//": inquire_dim time")
    if (time_nbr > time_nbr_in) error stop "time_nbr > time_nbr_in"
    if (xtr_dat_nbr /= xtr_idx_end) error stop "ERROR: xtr_dat_nbr /= xtr_idx_end"
    
    rcd=nf90_wrp_inq_varid(nc_id,"time",time_id,rcd_opt=NF90_ENOTVAR)
    ! Fill time with monotonic array if not present in external forcing file
    if (rcd == nf90_noerr) then 
       rcd=nf90_wrp(nf90_get_var(nc_id,time_id,time),sbr_nm//": gv time")
    else
       ! Implied do loop initialization
       time=real((/(idx,idx=1,time_nbr)/))
    endif ! endif 

    ! If variable ID is present, then set flag
    ! Get variable later
    ! fxm: Removed hard-wired indices
    rcd=nf90_wrp_inq_varid(nc_id,"asp_rat_lps",xtr_dat_frc_id(xtr_idx_asp_rat_lps),rcd_opt=NF90_ENOTVAR)
    if (rcd == nf90_noerr) xtr_dat_flg(xtr_idx_asp_rat_lps)=.true.

    rcd=nf90_wrp_inq_varid(nc_id,"doy",xtr_dat_frc_id(xtr_idx_doy),rcd_opt=NF90_ENOTVAR)
    if (rcd == nf90_noerr) xtr_dat_flg(xtr_idx_doy)=.true.

    rcd=nf90_wrp_inq_varid(nc_id,"hgt_mdp",xtr_dat_frc_id(xtr_idx_hgt_mdp),rcd_opt=NF90_ENOTVAR)
    if (rcd == nf90_noerr) xtr_dat_flg(xtr_idx_hgt_mdp)=.true.

    rcd=nf90_wrp_inq_varid(nc_id,"oro",xtr_dat_frc_id(xtr_idx_oro),rcd_opt=NF90_ENOTVAR)
    if (rcd == nf90_noerr) xtr_dat_flg(xtr_idx_oro)=.true.

    rcd=nf90_wrp_inq_varid(nc_id,"prs_mdp",xtr_dat_frc_id(xtr_idx_prs_mdp),rcd_opt=NF90_ENOTVAR)
    if (rcd == nf90_noerr) xtr_dat_flg(xtr_idx_prs_mdp)=.true.

    rcd=nf90_wrp_inq_varid(nc_id,"prs_ntf",xtr_dat_frc_id(xtr_idx_prs_ntf),rcd_opt=NF90_ENOTVAR)
    if (rcd == nf90_noerr) xtr_dat_flg(xtr_idx_prs_ntf)=.true.

    rcd=nf90_wrp_inq_varid(nc_id,"q_H2O_vpr",xtr_dat_frc_id(xtr_idx_q_H2O_vpr),rcd_opt=NF90_ENOTVAR)
    if (rcd == nf90_noerr) xtr_dat_flg(xtr_idx_q_H2O_vpr)=.true.

    rcd=nf90_wrp_inq_varid(nc_id,"sfc_typ",xtr_dat_frc_id(xtr_idx_sfc_typ),rcd_opt=NF90_ENOTVAR)
    if (rcd == nf90_noerr) xtr_dat_flg(xtr_idx_sfc_typ)=.true.

    rcd=nf90_wrp_inq_varid(nc_id,"tpt_gnd",xtr_dat_frc_id(xtr_idx_tpt_gnd),rcd_opt=NF90_ENOTVAR)
    if (rcd == nf90_noerr) xtr_dat_flg(xtr_idx_tpt_gnd)=.true.

    rcd=nf90_wrp_inq_varid(nc_id,"tpt_ice",xtr_dat_frc_id(xtr_idx_tpt_ice),rcd_opt=NF90_ENOTVAR)
    if (rcd == nf90_noerr) xtr_dat_flg(xtr_idx_tpt_ice)=.true.

    rcd=nf90_wrp_inq_varid(nc_id,"tpt_mdp",xtr_dat_frc_id(xtr_idx_tpt_mdp),rcd_opt=NF90_ENOTVAR)
    if (rcd == nf90_noerr) xtr_dat_flg(xtr_idx_tpt_mdp)=.true.

    rcd=nf90_wrp_inq_varid(nc_id,"tpt_soi",xtr_dat_frc_id(xtr_idx_tpt_soi),rcd_opt=NF90_ENOTVAR)
    if (rcd == nf90_noerr) xtr_dat_flg(xtr_idx_tpt_soi)=.true.

    rcd=nf90_wrp_inq_varid(nc_id,"tpt_sst",xtr_dat_frc_id(xtr_idx_tpt_sst),rcd_opt=NF90_ENOTVAR)
    if (rcd == nf90_noerr) xtr_dat_flg(xtr_idx_tpt_sst)=.true.

    rcd=nf90_wrp_inq_varid(nc_id,"vai_dst",xtr_dat_frc_id(xtr_idx_vai_dst),rcd_opt=NF90_ENOTVAR)
    if (rcd == nf90_noerr) xtr_dat_flg(xtr_idx_vai_dst)=.true.

    rcd=nf90_wrp_inq_varid(nc_id,"vwc_sfc",xtr_dat_frc_id(xtr_idx_vwc_sfc),rcd_opt=NF90_ENOTVAR)
    if (rcd == nf90_noerr) xtr_dat_flg(xtr_idx_vwc_sfc)=.true.

    rcd=nf90_wrp_inq_varid(nc_id,"wnd_mrd_mdp",xtr_dat_frc_id(xtr_idx_wnd_mrd_mdp),rcd_opt=NF90_ENOTVAR)
    if (rcd == nf90_noerr) xtr_dat_flg(xtr_idx_wnd_mrd_mdp)=.true.

    rcd=nf90_wrp_inq_varid(nc_id,"wnd_znl_mdp",xtr_dat_frc_id(xtr_idx_wnd_znl_mdp),rcd_opt=NF90_ENOTVAR)
    if (rcd == nf90_noerr) xtr_dat_flg(xtr_idx_wnd_znl_mdp)=.true.
       
    ! Get variables which are present
    do idx=1,xtr_dat_nbr
       if (xtr_dat_flg(idx)) then
          rcd=nf90_wrp(nf90_get_var(nc_id,xtr_dat_frc_id(idx),xtr_dat_frc(:,1,idx)),sbr_nm//": gv idx")
       endif ! end if variable present
    end do ! end loop over possible variables
    
    ! Verify time coordinate is monotonic increasing
    if (.not.mnt_ncr_chk(time,time_nbr)) then
       ! fxm: should improve input dataset
       !error stop "time coordinate not monotonic increasing"
       ! Implied do loop initialization
       time=real((/(idx,idx=1,time_nbr)/))
    endif ! endif
    if (rcd > 0) write(6,'(a)') "WARNING: rcd > 0 in xtr_dat_get"
  end subroutine xtr_dat_get
  
  subroutine xtr_dat_get_EuF( & ! [sbr] Read external forcing data from Eureka Flats
       fl_in, & ! [sng] netCDF file with external forcing data
       time_nbr, & ! [nbr] Number of timesteps to simulate
       plond, & ! [nbr] Number of longitudes
       xtr_dat_nbr, & ! [nbr] Number of external forcing data variables
       xtr_dat_frc) ! [frc] External forcing data (time_nbr,plond,xtr_dat_nbr)
    ! Purpose: Set time-varying zonal or meridional wind speed
    ! Routine was initially used for Eureka Flats data
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use shr_kind_mod,only:r8=>shr_kind_r8 ! [mdl] Precision r8, i8, ...
    use netcdf    ! [mdl] netCDF interface
    use nf90_utl  ! [mdl] netCDF utilities
    use utl_mdl   ! [mdl] Utility functions (date_time_get,mnt_chk...)
    implicit none
    character(len=*),parameter::sbr_nm='xtr_dat_get_EuF' ! [sng] Subroutine name
    character(len=*),intent(in)::fl_in ! [sng] netCDF file with external forcing data
    integer,intent(in)::time_nbr ! [nbr] Number of timesteps to simulate
    integer,intent(in)::plond ! [nbr] Number of longitudes
    integer,intent(in)::xtr_dat_nbr ! [nbr] Number of external forcing data variables
    real(r8),intent(out)::xtr_dat_frc(time_nbr,plond,xtr_dat_nbr) ! O [frc] External forcing data
    
    integer nc_id             ! [id] File ID
    integer rcd               ! [rcd] Return success code
    integer time_dmn_id       ! [enm] Dimension ID for time
    integer time_id           ! [enm] Variable ID
    integer time_nbr_in       ! [nbr] Dimension size
    integer wnd_mrd_mdp_hst_id
    integer wnd_znl_mdp_hst_id
    real(r8) time(time_nbr)   ! [day] Time coordinate (day of year)
    
    ! Open netCDF file
    rcd=nf90_wrp_open(fl_in,nf90_nowrite,nc_id,sbr_nm=sbr_nm)
    
    ! Get dimension IDs
    rcd=nf90_wrp_inq_dimid(nc_id,"time",time_dmn_id)
    rcd=nf90_wrp(nf90_inquire_dimension(nc_id,time_dmn_id,len=time_nbr_in),sbr_nm//": inquire_dim time")
    if (time_nbr > time_nbr_in) error stop "time_nbr > time_nbr_in"
    
    ! Get variable IDs
    rcd=nf90_wrp_inq_varid(nc_id,"time",time_id)
    rcd=nf90_wrp_inq_varid(nc_id,"wnd_mrd_mdp_hst",wnd_mrd_mdp_hst_id)
    rcd=nf90_wrp_inq_varid(nc_id,"wnd_znl_mdp_hst",wnd_znl_mdp_hst_id)
    
    ! Get time and wind data
    rcd=nf90_wrp(nf90_get_var(nc_id,time_id,time),sbr_nm//": get_var time")
    rcd=nf90_wrp(nf90_get_var(nc_id,wnd_mrd_mdp_hst_id,xtr_dat_frc(:,1,1)),sbr_nm//": get_var wnd_mrd_mdp_hst")
    rcd=nf90_wrp(nf90_get_var(nc_id,wnd_znl_mdp_hst_id,xtr_dat_frc(:,1,2)),sbr_nm//": get_var wnd_znl_mdp_hst")
    
    ! Verify time coordinate is monotonic increasing
    if(.not.mnt_ncr_chk(time,time_nbr)) error stop "time coordinate not monotonic increasing"
    write (6,'(a)') "Ingested "//fl_in
    if (dbg_lvl >= dbg_fl) then
       write (6,'(a,i3)') "INFO: "//sbr_nm//"() reports time_nbr = ",time_nbr
       write (6,'(a,f6.3)') "INFO: "//sbr_nm//"() reports max mrd wind = ",maxval(xtr_dat_frc(:,:,1))
       write (6,'(a,f6.3)') "INFO: "//sbr_nm//"() reports min mrd wind = ",minval(xtr_dat_frc(:,:,1))
       write (6,'(a,f6.3)') "INFO: "//sbr_nm//"() reports max znl wind = ",maxval(xtr_dat_frc(:,:,2))
       write (6,'(a,f6.3)') "INFO: "//sbr_nm//"() reports  min znl wind = ",minval(xtr_dat_frc(:,:,2))
    endif ! endif dbg
  end subroutine xtr_dat_get_EuF
  
  subroutine xtr_dat_set_flt_1D( & ! [sbr] Set time-varying forcing data, real types, 1D arrays
       time_nbr, & ! I [nbr] Number of timesteps to simulate
       plond, & ! I [nbr] Number of longitudes
       nstep, & ! I [idx] Timestep index
       xtr_dat_nbr, & ! I [nbr] Number of external forcing data variables
       xtr_dat_frc, & ! I [frc] External forcing data (time_nbr,plond,xtr_dat_nbr)
       xtr_dat_idx, & ! I [idx] External forcing data variable index
       xtr_dat_now) ! O [frc] Instantaneous forcing data
    ! Purpose: Set time-varying zonal or meridional wind speed
    use shr_kind_mod,only:r8=>shr_kind_r8 ! [mdl] Precision r8, i8, ...
    implicit none
    integer,intent(in)::time_nbr ! [nbr] Number of timesteps to simulate
    integer,intent(in)::plond ! [nbr] Number of longitudes
    integer,intent(in)::nstep ! [idx] Timestep index
    integer,intent(in)::xtr_dat_nbr ! [nbr] Number of external forcing data variables
    real(r8),intent(in)::xtr_dat_frc(time_nbr,plond,xtr_dat_nbr) ! I [frc] External forcing data
    integer,intent(in)::xtr_dat_idx ! [idx] External forcing data variable index
    real(r8),intent(out)::xtr_dat_now(plond) ! [frc] Instantaneous forcing data
    xtr_dat_now=xtr_dat_frc(nstep,:,xtr_dat_idx) ! [frc] Instantaneous forcing data
  end subroutine xtr_dat_set_flt_1D
  
  subroutine xtr_dat_set_flt_scl( & ! [sbr] Set time-varying forcing data, real types, scalars
       time_nbr, & ! I [nbr] Number of timesteps to simulate
       plond, & ! I [nbr] Number of longitudes
       nstep, & ! I [idx] Timestep index
       xtr_dat_nbr, & ! I [nbr] Number of external forcing data variables
       xtr_dat_frc, & ! I [frc] External forcing data (time_nbr,plond,xtr_dat_nbr)
       xtr_dat_idx, & ! I [idx] External forcing data variable index
       xtr_dat_now) ! O [frc] Instantaneous forcing data
    ! Purpose: Set time-varying zonal or meridional wind speed
    use shr_kind_mod,only:r8=>shr_kind_r8 ! [mdl] Precision r8, i8, ...
    implicit none
    integer,intent(in)::time_nbr ! [nbr] Number of timesteps to simulate
    integer,intent(in)::plond ! [nbr] Number of longitudes
    integer,intent(in)::nstep ! [idx] Timestep index
    integer,intent(in)::xtr_dat_nbr ! [nbr] Number of external forcing data variables
    real(r8),intent(in)::xtr_dat_frc(time_nbr,plond,xtr_dat_nbr) ! I [frc] External forcing data
    integer,intent(in)::xtr_dat_idx ! [idx] External forcing data variable index
    real(r8),intent(out)::xtr_dat_now ! [frc] Instantaneous forcing data
    xtr_dat_now=xtr_dat_frc(nstep,1,xtr_dat_idx) ! [frc] Instantaneous forcing data
  end subroutine xtr_dat_set_flt_scl
  
  subroutine xtr_dat_set_int_1D( & ! [sbr] Set time-varying forcing data, integer types, 1D arrays
       time_nbr, & ! I [nbr] Number of timesteps to simulate
       plond, & ! I [nbr] Number of longitudes
       nstep, & ! I [idx] Timestep index
       xtr_dat_nbr, & ! I [nbr] Number of external forcing data variables
       xtr_dat_frc, & ! I [frc] External forcing data (time_nbr,plond,xtr_dat_nbr)
       xtr_dat_idx, & ! I [idx] External forcing data variable index
       xtr_dat_now) ! O [frc] Instantaneous forcing data
    ! Purpose: set time-varying zonal or meridional wind speed
    use shr_kind_mod,only:r8=>shr_kind_r8 ! [mdl] Precision r8, i8, ...
    implicit none
    integer,intent(in)::time_nbr ! [nbr] Number of timesteps to simulate
    integer,intent(in)::plond ! [nbr] Number of longitudes
    integer,intent(in)::nstep ! [idx] Timestep index
    integer,intent(in)::xtr_dat_nbr ! [nbr] Number of external forcing data variables
    real(r8),intent(in)::xtr_dat_frc(time_nbr,plond,xtr_dat_nbr) ! I External forcing data
    integer,intent(in)::xtr_dat_idx ! [idx] External forcing data variable index
    integer,intent(out)::xtr_dat_now(plond) ! [frc] Instantaneous forcing data
    xtr_dat_now=nint(xtr_dat_frc(nstep,:,xtr_dat_idx)) ! [frc] Instantaneous forcing data
    end subroutine xtr_dat_set_int_1D
  
end module dstxtr
