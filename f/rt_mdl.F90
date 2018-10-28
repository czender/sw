! $Id$

! Purpose: RT utility routines used by swnb, nbm, ck, and lbl

! Copyright (C) 1994--2018 Charlie Zender
! This software is distributed under the terms of the GNU General Public License
! See http://www.gnu.org/copyleft/gpl.html for full license text

! Contains routines formerly in csz_f77.F, csz_mdl.F90

! Usage: 
! use rt_mdl ! [mdl] Radiative transfer utilities

module rt_mdl ! [mdl] Radiative transfer utilities
  implicit none
  public ! [stt] Symbols are public unless individually qualified as private
  
contains
  
  subroutine abs_xsx_get(fl_in,bnd_nbr, &
       abs_xsx,wvl_grd, &
       abs_xsx_dadT,tpt_std, &
       qnt_yld, &
       wvl_ctr,wvl_min,wvl_max)
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use netcdf ! [mdl] netCDF interface
    use nf90_utl ! [mdl] netCDF utilities
    use wvl_mdl ! [mdl] Wavelength grid parameters
    implicit none
    ! Parameters
    character(len=*),parameter::sbr_nm='abs_xsx_get' ! [sng] Subroutine name
    ! Input
    character(len=*),intent(in)::fl_in
    ! Output
    integer,intent(out)::bnd_nbr
    real,dimension(:),allocatable,intent(out)::abs_xsx
    real,dimension(:),allocatable,intent(out)::wvl_grd
    real,optional,intent(out)::tpt_std
    real,dimension(:),allocatable,optional,intent(out)::abs_xsx_dadT
    real,dimension(:),allocatable,optional,intent(out)::qnt_yld
    real,dimension(:),allocatable,optional,intent(out)::wvl_ctr
    real,dimension(:),allocatable,optional,intent(out)::wvl_min
    real,dimension(:),allocatable,optional,intent(out)::wvl_max

    ! Local
    integer::bnd_dmn_id
    integer::cnt_bnd(1)
    integer::cnt_bndp(1)
    integer::nc_id
    integer::rcd=nf90_noerr ! [rcd] Return success code
    integer::srt_one(1)

    ! Cross-section file input variables
    integer::abs_xsx_id
    integer::abs_xsx_dadT_id
    integer::qnt_yld_id
    integer::tpt_std_id
    integer::wvl_ctr_id
    integer::wvl_grd_id
    integer::wvl_min_id
    integer::wvl_max_id
    
    ! Main Code
    if (dbg_lvl >= dbg_sbr) write (6,'(a)') 'Entering '//sbr_nm
    
    ! Ingest fl_in
    rcd=nf90_wrp_open(fl_in,nf90_nowrite,nc_id)
    ! Get dimension IDs
    rcd=nf90_wrp_inq_dimid(nc_id,'bnd',bnd_dmn_id)
    
    ! Get dimension sizes
    rcd=nf90_wrp(nf90_inquire_dimension(nc_id,bnd_dmn_id,len=bnd_nbr),sbr_nm//': inquire_dim bnd')
    !if (bnd_nbr>bnd_nbr_max) stop 'bnd_nbr>bnd_nbr_max'
    cnt_bnd(1)=bnd_nbr
    cnt_bndp(1)=bnd_nbr+1
    
    allocate(abs_xsx(bnd_nbr),stat=rcd)
    if(rcd /= 0) stop 'allocate() failed for abs_xsx'
    allocate(wvl_grd(bnd_nbr+1),stat=rcd)
    if(rcd /= 0) stop 'allocate() failed for wvl_grd'

    ! Get variable IDs
    rcd=nf90_wrp_inq_varid(nc_id,'abs_xsx',abs_xsx_id)
    rcd=nf90_wrp_inq_varid(nc_id,'wvl_grd',wvl_grd_id)

    ! Get data
    rcd=nf90_wrp(nf90_get_var(nc_id,abs_xsx_id,abs_xsx,srt_one,cnt_bnd),'gv abs_xsx')
    rcd=nf90_wrp(nf90_get_var(nc_id,wvl_grd_id,wvl_grd,srt_one,cnt_bndp),'gv wvl_grd')

    if(present(abs_xsx_dadT)) then
       allocate(abs_xsx_dadT(bnd_nbr),stat=rcd)
       if(rcd /= 0) stop 'allocate() failed for abs_xsx_dadT'
       rcd=nf90_wrp_inq_varid(nc_id,'abs_xsx_dadT',abs_xsx_dadT_id,rcd_opt=nf90_enotvar) ! Tolerate variable absence
       if(rcd /= nf90_noerr) then
          rcd=nf90_wrp(nf90_get_var(nc_id,abs_xsx_dadT_id,abs_xsx_dadT,srt_one,cnt_bnd),'gv abs_xsx_dadT')
       else
          abs_xsx_dadT(:)=0.0
       endif ! !rcd
    endif ! !abs_xsx_dadT
    if(present(tpt_std)) then
       rcd=nf90_wrp_inq_varid(nc_id,'tpt_std',tpt_std_id)
       rcd=nf90_wrp(nf90_get_var(nc_id,tpt_std_id,tpt_std),'gv tpt_std')
    endif ! !tpt_std
    if(present(qnt_yld)) then
       allocate(qnt_yld(bnd_nbr),stat=rcd)
       if(rcd /= 0) stop 'allocate() failed for qnt_yld'
       rcd=nf90_wrp_inq_varid(nc_id,'qnt_yld',qnt_yld_id,rcd_opt=nf90_enotvar) ! Tolerate variable absence
       if(rcd /= nf90_noerr) then
          rcd=nf90_wrp(nf90_get_var(nc_id,qnt_yld_id,qnt_yld,srt_one,cnt_bnd),'gv qnt_yld')
       else
          qnt_yld(:)=0.0
       endif ! !rcd
    endif ! !qnt_yld
    if(present(wvl_ctr)) then
       allocate(wvl_ctr(bnd_nbr),stat=rcd)
       if(rcd /= 0) stop 'allocate() failed for wvl_ctr'
       rcd=nf90_wrp_inq_varid(nc_id,'wvl_ctr',wvl_ctr_id)
       rcd=nf90_wrp(nf90_get_var(nc_id,wvl_ctr_id,wvl_ctr,srt_one,cnt_bnd),'gv wvl_ctr')
    endif ! !wvl_ctr
    if(present(wvl_min)) then
       allocate(wvl_min(bnd_nbr),stat=rcd)
       if(rcd /= 0) stop 'allocate() failed for wvl_min'
       rcd=nf90_wrp_inq_varid(nc_id,'wvl_min',wvl_min_id)
       rcd=nf90_wrp(nf90_get_var(nc_id,wvl_min_id,wvl_min,srt_one,cnt_bnd),'gv wvl_min')
    endif ! !wvl_min
    if(present(wvl_max)) then
       allocate(wvl_max(bnd_nbr),stat=rcd)
       if(rcd /= 0) stop 'allocate() failed for wvl_max'
       rcd=nf90_wrp_inq_varid(nc_id,'wvl_max',wvl_max_id)
       rcd=nf90_wrp(nf90_get_var(nc_id,wvl_max_id,wvl_max,srt_one,cnt_bnd),'gv wvl_max')
    endif ! !wvl_max

    ! Close file
    rcd=nf90_wrp_close(nc_id,fl_in,'Ingested') ! [fnc] Close file
    
    if (dbg_lvl >= dbg_sbr) write (6,'(a)') 'Exiting '//sbr_nm
    return
  end subroutine abs_xsx_get                       ! end abs_xsx_get()

  subroutine mlk_bnd_prm_get(fl_in,bnd_nbr, &
       A_phi,A_psi,B_phi,B_psi,S_d_abs_cff_mss,S_p_abs_cff_mss, &
       wvl_ctr,wvl_min,wvl_max)
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use netcdf ! [mdl] netCDF interface
    use nf90_utl ! [mdl] netCDF utilities
    use wvl_mdl ! [mdl] Wavelength grid parameters
    implicit none
    ! Parameters
    character(len=*),parameter::sbr_nm='mlk_bnd_prm_get' ! [sng] Subroutine name
    ! Input
    character(len=*),intent(in)::fl_in
    ! Output
    integer,intent(out)::bnd_nbr
    real,dimension(:),allocatable,intent(out)::A_phi
    real,dimension(:),allocatable,intent(out)::A_psi
    real,dimension(:),allocatable,intent(out)::B_phi
    real,dimension(:),allocatable,intent(out)::B_psi
    real,dimension(:),allocatable,intent(out)::S_d_abs_cff_mss
    real,dimension(:),allocatable,intent(out)::S_p_abs_cff_mss
    real,dimension(:),allocatable,intent(out)::wvl_max
    real,dimension(:),allocatable,intent(out)::wvl_min
    real,dimension(:),allocatable,intent(out)::wvl_ctr
    ! Local
    integer::bnd_dmn_id
    integer::cnt_bnd(1)
    integer::nc_id
    integer::rcd=nf90_noerr ! [rcd] Return success code
    integer::srt_one(1)

    ! Narrow band CO2 input variables
    integer::A_phi_id
    integer::A_psi_id
    integer::B_phi_id
    integer::B_psi_id
    integer::S_d_abs_cff_mss_id
    integer::S_p_abs_cff_mss_id
    integer::wvl_max_id
    integer::wvl_min_id
    integer::wvl_ctr_id
    
    ! Main Code
    if (dbg_lvl >= dbg_sbr) write (6,'(a)') 'Entering '//sbr_nm
    
    ! Ingest fl_in
    rcd=nf90_wrp_open(fl_in,nf90_nowrite,nc_id)
    ! Get dimension IDs
    rcd=nf90_wrp_inq_dimid(nc_id,'bnd',bnd_dmn_id)
    
    ! Get dimension sizes
    rcd=nf90_wrp(nf90_inquire_dimension(nc_id,bnd_dmn_id,len=bnd_nbr),sbr_nm//": inquire_dim bnd")
    !if (bnd_nbr>bnd_nbr_max) stop 'bnd_nbr>bnd_nbr_max'
    cnt_bnd(1)=bnd_nbr
    
    allocate(A_phi(bnd_nbr),stat=rcd)
    if(rcd /= 0) stop "allocate() failed for A_phi"
    allocate(A_psi(bnd_nbr),stat=rcd)
    if(rcd /= 0) stop "allocate() failed for A_psi"
    allocate(B_phi(bnd_nbr),stat=rcd)
    if(rcd /= 0) stop "allocate() failed for B_phi"
    allocate(B_psi(bnd_nbr),stat=rcd)
    if(rcd /= 0) stop "allocate() failed for B_psi"
    allocate(S_d_abs_cff_mss(bnd_nbr),stat=rcd)
    if(rcd /= 0) stop "allocate() failed for S_d_abs_cff_mss"
    allocate(S_p_abs_cff_mss(bnd_nbr),stat=rcd)
    if(rcd /= 0) stop "allocate() failed for S_p_abs_cff_mss"
    allocate(wvl_max(bnd_nbr),stat=rcd)
    if(rcd /= 0) stop "allocate() failed for wvl_max"
    allocate(wvl_min(bnd_nbr),stat=rcd)
    if(rcd /= 0) stop "allocate() failed for wvl_min"
    allocate(wvl_ctr(bnd_nbr),stat=rcd)
    if(rcd /= 0) stop "allocate() failed for wvl_ctr"
    
    ! Get variable IDs
    rcd=nf90_wrp_inq_varid(nc_id,'A_phi',A_phi_id)
    rcd=nf90_wrp_inq_varid(nc_id,'A_psi',A_psi_id)
    rcd=nf90_wrp_inq_varid(nc_id,'B_phi',B_phi_id)
    rcd=nf90_wrp_inq_varid(nc_id,'B_psi',B_psi_id)
    rcd=nf90_wrp_inq_varid(nc_id,'S_d_abs_cff_mss',S_d_abs_cff_mss_id)
    rcd=nf90_wrp_inq_varid(nc_id,'S_p_abs_cff_mss',S_p_abs_cff_mss_id)
    rcd=nf90_wrp_inq_varid(nc_id,'wvl_max',wvl_max_id)
    rcd=nf90_wrp_inq_varid(nc_id,'wvl_min',wvl_min_id)
    rcd=nf90_wrp_inq_varid(nc_id,'wvl_ctr',wvl_ctr_id)
    
    ! Get data
    rcd=nf90_wrp(nf90_get_var(nc_id,A_phi_id,A_phi,srt_one,cnt_bnd),"gv A_phi")
    rcd=nf90_wrp(nf90_get_var(nc_id,A_psi_id,A_psi,srt_one,cnt_bnd),"gv A_psi")
    rcd=nf90_wrp(nf90_get_var(nc_id,B_phi_id,B_phi,srt_one,cnt_bnd),"gv B_phi")
    rcd=nf90_wrp(nf90_get_var(nc_id,B_psi_id,B_psi,srt_one,cnt_bnd),"gv B_psi")
    rcd=nf90_wrp(nf90_get_var(nc_id,S_d_abs_cff_mss_id,S_d_abs_cff_mss,srt_one,cnt_bnd),"gv S_d_abs_cff_mss")
    rcd=nf90_wrp(nf90_get_var(nc_id,S_p_abs_cff_mss_id,S_p_abs_cff_mss,srt_one,cnt_bnd),"gv S_p_abs_cff_mss")
    rcd=nf90_wrp(nf90_get_var(nc_id,wvl_max_id,wvl_max,srt_one,cnt_bnd),"gv wvl_max")
    rcd=nf90_wrp(nf90_get_var(nc_id,wvl_min_id,wvl_min,srt_one,cnt_bnd),"gv wvl_min")
    rcd=nf90_wrp(nf90_get_var(nc_id,wvl_ctr_id,wvl_ctr,srt_one,cnt_bnd),"gv wvl_ctr")
    ! Close file
    rcd=nf90_wrp_close(nc_id,fl_in,'Ingested') ! [fnc] Close file
    
    if (dbg_lvl >= dbg_sbr) write (6,'(a)') 'Exiting '//sbr_nm
    return
  end subroutine mlk_bnd_prm_get                       ! end mlk_bnd_prm_get()

  subroutine slr_spc_get_CCM(fl_slr,wvl_min,wvl_max,wvl_nbr,flx_slr_frc,xtr_typ_LHS,xtr_typ_RHS)
    ! Purpose: Wrapper routine for slr_spc_get which takes care of non-monotonicity of CCM wavelength grid
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use wvl_mdl ! [mdl] Wavelength grid parameters
    implicit none
    ! Input
    character(len=*),intent(in)::fl_slr
    integer,intent(in)::wvl_nbr
    integer,intent(in)::xtr_typ_LHS
    integer,intent(in)::xtr_typ_RHS
    real,intent(in)::wvl_min(wvl_nbr) 
    real,intent(in)::wvl_max(wvl_nbr)
    ! Output
    real,intent(out)::flx_slr_frc(wvl_nbr)
    ! Local
    integer::wvl_nbr_tmp
    
    ! Main Code
    if (dbg_lvl >= dbg_sbr) write (6,'(a)') 'Entering slr_spc_get_CCM()'
    
    ! Initialize some variables
    if (wvl_nbr > bnd_nbr_CCM_SW_max) stop 'wvl_nbr > bnd_nbr_CCM_SW_max in slr_spc_get_CCM()'
    if (wvl_nbr /= bnd_nbr_CCM_SW_max) stop 'wvl_nbr /= bnd_nbr_CCM_SW_max in slr_spc_get_CCM()'
    
    ! Get first 16 bands
    wvl_nbr_tmp=16
    call slr_spc_get(fl_slr,wvl_min,wvl_max,wvl_nbr_tmp,flx_slr_frc,xtr_typ_LHS,xtr_typ_RHS)
    ! Get band 17
    wvl_nbr_tmp=1
    call slr_spc_get(fl_slr,wvl_min(17),wvl_max(17),wvl_nbr_tmp,flx_slr_frc(17),xtr_typ_LHS,xtr_typ_RHS)
    ! Get bands 18 and 19
    wvl_nbr_tmp=2
    call slr_spc_get(fl_slr,wvl_min(18),wvl_max(18),wvl_nbr_tmp,flx_slr_frc(18),xtr_typ_LHS,xtr_typ_RHS)
    
    if (dbg_lvl >= dbg_sbr) write (6,'(a)') 'Exiting slr_spc_get_CCM()'
    return
  end subroutine slr_spc_get_CCM                       ! end slr_spc_get_CCM()
  
  subroutine rbn_vec_CCM( &
       in_nbr,grd_in,dat_in, &
       out_nbr,crd_out_min,crd_out_max,dat_out, &
       xtr_typ_LHS,xtr_typ_RHS)
    ! Purpose: Wrapper routine for rbn_vec which takes care of non-monotonicity of CCM wavelength grid
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use vec_mdl ! [mdl] Vector manipulation, interpolation, rebinning
    use wvl_mdl ! [mdl] Wavelength grid parameters
    implicit none
    ! Input
    integer,intent(in)::xtr_typ_LHS
    integer,intent(in)::xtr_typ_RHS
    integer,intent(in)::in_nbr
    integer,intent(in)::out_nbr
    real,intent(in)::crd_out_min(out_nbr) ! Output grid
    real,intent(in)::crd_out_max(out_nbr) ! Output grid
    real,intent(in)::grd_in(in_nbr+1)     ! Input grid
    real,intent(in)::dat_in(in_nbr)       ! Input data
    ! Output
    real,intent(out)::dat_out(out_nbr)     ! Input data rebinned to output grid
    ! Local
    integer::idx               ! [idx] Counting index
    integer::out_idx_srt       ! [idx] Starting index of current output block
    integer::out_nbr_tmp       ! [nbr] Number of output bins in current output block
    real grd_out(bnd_nbr_CCM_SW_max+1) ! Output grid
    
    ! Main Code
    if (dbg_lvl >= dbg_sbr) write (6,'(a)') 'Entering rbn_vec_CCM()'
    
    ! Initialize some variables
    if (out_nbr > bnd_nbr_CCM_SW_max) stop 'out_nbr > bnd_nbr_CCM_SW_max in rbn_vec_CCM()'
    if (out_nbr /= bnd_nbr_CCM_SW_max) stop 'out_nbr /= bnd_nbr_CCM_SW_max in rbn_vec_CCM()'
    
    ! Get first 16 bands
    out_idx_srt=1
    out_nbr_tmp=16
    do idx=out_idx_srt,out_idx_srt+out_nbr_tmp-1
       grd_out(idx)=crd_out_min(idx)
    enddo                     ! end loop over idx
    grd_out(out_idx_srt+out_nbr_tmp)=crd_out_max(out_idx_srt+out_nbr_tmp-1)
    call rbn_vec(in_nbr,grd_in,dat_in,out_nbr_tmp,grd_out(out_idx_srt),dat_out(out_idx_srt),xtr_typ_LHS,xtr_typ_RHS)
    
    ! Get band 17
    out_idx_srt=17
    out_nbr_tmp=1
    do idx=out_idx_srt,out_idx_srt+out_nbr_tmp-1
       grd_out(idx)=crd_out_min(idx)
    enddo                     ! end loop over idx
    grd_out(out_idx_srt+out_nbr_tmp)=crd_out_max(out_idx_srt+out_nbr_tmp-1)
    call rbn_vec(in_nbr,grd_in,dat_in,out_nbr_tmp,grd_out(out_idx_srt),dat_out(out_idx_srt),xtr_typ_LHS,xtr_typ_RHS)
    
    ! Get bands 18 and 19
    out_idx_srt=18
    out_nbr_tmp=2
    do idx=out_idx_srt,out_idx_srt+out_nbr_tmp-1
       grd_out(idx)=crd_out_min(idx)
    enddo                     ! end loop over idx
    grd_out(out_idx_srt+out_nbr_tmp)=crd_out_max(out_idx_srt+out_nbr_tmp-1)
    call rbn_vec(in_nbr,grd_in,dat_in,out_nbr_tmp,grd_out(out_idx_srt),dat_out(out_idx_srt),xtr_typ_LHS,xtr_typ_RHS)
    
    if (dbg_lvl >= dbg_sbr) write (6,'(a)') 'Exiting rbn_vec_CCM()'
    return
  end subroutine rbn_vec_CCM                       ! end rbn_vec_CCM()
  
  subroutine slr_spc_get(fl_slr,wvl_min_out,wvl_max_out,wvl_out_nbr,flx_frc_out,xtr_typ_LHS,xtr_typ_RHS)
    ! Purpose: Retrieve fractional solar spectrum from file
    ! Suffix _in refers to arrays stored in file
    ! Suffix _out refers to arrays passed into and out of this subroutine
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use netcdf ! [mdl] netCDF interface
    use nf90_utl ! [mdl] netCDF utilities
    use sng_mdl,only:ftn_strlen ! [mdl] String manipulation
    use utl_mdl ! [mdl] Utility functions (mnt_chk,...)
    use vec_mdl ! [mdl] Vector manipulation, interpolation, rebinning
    use xtr_mdl ! [mdl] Extrapolation/interpolation handling
    implicit none
    ! Parameters
    character(len=*),parameter::sbr_nm='slr_spc_get' ! [sng] Subroutine name
    integer,parameter::wvl_in_nbr_max=49934 ! Kur95 1 cm-1 spectrum
    ! Input
    character(len=*),intent(in)::fl_slr
    integer,intent(in)::wvl_out_nbr       ! [nbr] max(wvl_out_nbr)=205001 for NO2.F max(bnd_nbr_JPL,bnd_nbr_NCAR,bnd_nbr_NOAA)
    integer,intent(in)::xtr_typ_LHS
    integer,intent(in)::xtr_typ_RHS
    real,intent(in)::wvl_min_out(wvl_out_nbr) 
    real,intent(in)::wvl_max_out(wvl_out_nbr)
    ! Output
    real,intent(out)::flx_frc_out(wvl_out_nbr)
    ! Local
    integer::wvl_in_nbr        ! [nbr] Dimension size
    integer::rcd
    integer::nc_id
    integer::wvl_dim_id
    real wvl_min_out_mnt(wvl_out_nbr) ! [m] Minimum wavelength on montonically increasing grid
    real wvl_max_out_mnt(wvl_out_nbr) ! [m] Maximum wavelength on montonically increasing grid
    ! Solar spectrum input variables
    integer::flx_frc_blr_in_id
    integer::wvl_max_in_id
    integer::wvl_min_in_id
    ! Solar Spectrum input variables
    ! Allocatable variables
    real,dimension(:),allocatable::flx_frc_blr_in
    real,dimension(:),allocatable::wvl_max_in 
    real,dimension(:),allocatable::wvl_min_in 
    ! Local
    integer::out_idx
    logical mnt               ! Monotonicity flag
    logical REVERSE
    real flx_frc_blr_max_out(wvl_out_nbr)
    real flx_frc_blr_min_out(wvl_out_nbr)
    
    ! Main Code
    if (dbg_lvl >= dbg_sbr) write (6,'(a)') 'Entering slr_spc_get()'
    ! Initialize default values
    rcd=nf90_noerr ! nf90_noerr == 0
    REVERSE=.false.
    
    ! Check for unusual extrapolation types
    if (xtr_typ_LHS /= xtr_fll_ngh+xtr_vrb.and.xtr_typ_LHS /= xtr_fll_ngh)  &
         write (6,'(a)') 'WARNING: xtr_typ_LHS != xtr_fll_ngh in slr_spc_get()'
    if (xtr_typ_RHS /= xtr_fll_ngh+xtr_vrb.and.xtr_typ_LHS /= xtr_fll_ngh)  &
         write (6,'(a)') 'WARNING: xtr_typ_RHS != xtr_fll_ngh in slr_spc_get()'
    
    ! Convert to order of increasing wavelength if necessary
    wvl_min_out_mnt(:)=wvl_min_out(:) ! [m] Minimum wavelength on montonically increasing grid
    wvl_max_out_mnt(:)=wvl_max_out(:) ! [m] Maximum wavelength on montonically increasing grid
    if (wvl_min_out(1) > wvl_min_out(2)) then
       REVERSE=.true.
       call rvr_vec(wvl_min_out_mnt,wvl_out_nbr)
       call rvr_vec(wvl_max_out_mnt,wvl_out_nbr)
    endif                     ! endif reversing direction
    
    if (dbg_lvl >= dbg_io) then
       write (6,'(a,/,a,a,a,i4)') 'slr_spc_get():', &
            'fl_slr = ',fl_slr, &
            'wvl_out_nbr = ',wvl_out_nbr
       do out_idx=1,wvl_out_nbr
          write (6,'(2(a,i4,a,es10.3))')  &
               'wvl_min_out_mnt(',out_idx,') = ',wvl_min_out_mnt(out_idx), &
               ', wvl_max_out_mnt(',out_idx,') = ',wvl_max_out_mnt(out_idx)
       end do                 ! end loop over out_idx
    endif                     ! endif dbg
    
    ! Check for monotonicity
    mnt=mnt_chk(wvl_min_out_mnt,wvl_out_nbr)
    if (.not.mnt) stop 'wvl_min_out_mnt not monotonic in slr_spc_get()'
    
    ! Sanity check: even boundaries?
    do out_idx=1,wvl_out_nbr-1 ! Loop ends at wvl_out_nbr-1
       if (wvl_min_out_mnt(out_idx+1) /= wvl_max_out_mnt(out_idx)) then
          write (6,'(a,a)') prg_nm(1:ftn_strlen(prg_nm)), &
               ': ERROR slr_spc_get() reports uneven boundaries in output wavelength grid' 
          write (6,'(2(a,i4,a,es10.3))')   &
               'wvl_min_out_mnt(',out_idx+1,') = ',wvl_min_out_mnt(out_idx+1), &
               ' != wvl_max_out_mnt(',out_idx,') = ',wvl_max_out_mnt(out_idx)
          stop
       endif
    end do                    ! end loop over out_idx
    
    ! Get solar flux data
    rcd=nf90_wrp_open(fl_slr,nf90_nowrite,nc_id)
    ! Get dimension IDs
    rcd=nf90_wrp_inq_dimid(nc_id,'wvl',wvl_dim_id)
    ! Get dimension sizes
    rcd=nf90_wrp(nf90_inquire_dimension(nc_id,wvl_dim_id,len=wvl_in_nbr),sbr_nm//": inquire_dim wvl")
    if (wvl_in_nbr > wvl_in_nbr_max) stop 'wvl_in_nbr > wvl_in_nbr_max in slr_spc_get()'
    
    ! Allocate space for dynamic arrays
    allocate(flx_frc_blr_in(wvl_in_nbr),stat=rcd)
    allocate(wvl_max_in(wvl_in_nbr),stat=rcd) 
    allocate(wvl_min_in(wvl_in_nbr),stat=rcd) 
    ! Get variable IDs
    rcd=nf90_wrp_inq_varid(nc_id,'wvl_max',wvl_max_in_id)
    rcd=nf90_wrp_inq_varid(nc_id,'wvl_min',wvl_min_in_id)
    rcd=nf90_wrp_inq_varid(nc_id,'flx_frc_blr',flx_frc_blr_in_id)
    ! Get data
    rcd=nf90_wrp(nf90_get_var(nc_id,wvl_max_in_id,wvl_max_in),sbr_nm//": get_var wvl_max_in")
    rcd=nf90_wrp(nf90_get_var(nc_id,wvl_min_in_id,wvl_min_in),sbr_nm//": get_var wvl_min_in")
    rcd=nf90_wrp(nf90_get_var(nc_id,flx_frc_blr_in_id,flx_frc_blr_in),sbr_nm//": get_var flx_frc_blr_in")
    rcd=nf90_wrp_close(nc_id,fl_slr,'Ingested')
    
    ! NB: Recall slr_spc() defines flx_frc_blr(idx) to refer to fractional flux bluer than wvl_max(idx), 
    ! NOT wvl_ctr(idx) or wvl_min(idx)
    
    ! Check for monotonicity
    mnt=mnt_chk(wvl_min_in,wvl_in_nbr)
    if (.not.mnt) stop 'wvl_min_in not monotonic in slr_spc_get()'
    
    ! Interpolate flx_frc_blr_in to wvl_min_out_mnt
    call ntp_vec(wvl_in_nbr,wvl_max_in,flx_frc_blr_in,wvl_out_nbr,wvl_min_out_mnt,flx_frc_blr_min_out,xtr_typ_LHS,xtr_typ_RHS)
    ! Interpolate flx_frc_blr_in to wvl_max_out_mnt
    call ntp_vec(wvl_in_nbr,wvl_max_in,flx_frc_blr_in,wvl_out_nbr,wvl_max_out_mnt,flx_frc_blr_max_out,xtr_typ_LHS,xtr_typ_RHS)
    
    ! Reset extrapolation flags and strings 
    ! This call is not necessary, but included for safety and defensive programming
    call xtr_ini(xtr_typ_LHS,xtr_typ_RHS)
    
    ! Sanity check
    do out_idx=1,wvl_out_nbr-1
       if (flx_frc_blr_min_out(out_idx+1) /= flx_frc_blr_max_out(out_idx)) &
            stop 'flx_frc_blr_min_out(out_idx+1) /= flx_frc_blr_max_out(out_idx) in slr_spc_get()'
    end do                    ! end loop over out_idx
    
    ! Convert flx_frc_blr_[min,max]_out to flx_frc_out
    do out_idx=1,wvl_out_nbr
       flx_frc_out(out_idx)=flx_frc_blr_max_out(out_idx)-flx_frc_blr_min_out(out_idx)
    end do                    ! end loop over out_idx
    
    ! Revert to original order (decreasing wavelength) if necessary
    if (REVERSE) then
       call rvr_vec(flx_frc_out,wvl_out_nbr)
       !     call rvr_vec(wvl_min_out_mnt,wvl_out_nbr)
       !     call rvr_vec(wvl_max_out_mnt,wvl_out_nbr)
    endif                     ! endif reversing direction
    
    ! De-allocate dynamic variables
    if (allocated(flx_frc_blr_in)) deallocate(flx_frc_blr_in,stat=rcd)
    if (allocated(wvl_max_in)) deallocate(wvl_max_in,stat=rcd) 
    if (allocated(wvl_min_in)) deallocate(wvl_min_in,stat=rcd) 
    
    if (dbg_lvl >= dbg_sbr) write (6,'(a)') 'Exiting slr_spc_get()'
    return
  end subroutine slr_spc_get                       ! end slr_spc_get()
  
  subroutine wgt_get(fl_in,var_nm,wvl_grd_out,wvl_out_nbr,var_out)
    ! Purpose: Retrieve an array of spectral weights from file and remap to specified grid
    ! Subscript _in refers to arrays stored in file
    ! Subscript _out refers to arrays passed into and out of this subroutine
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use netcdf ! [mdl] netCDF interface
    use nf90_utl ! [mdl] netCDF utilities
    use vec_mdl ! [mdl] Vector manipulation, interpolation, rebinning
    use xtr_mdl ! [mdl] Extrapolation/interpolation handling
    implicit none
    ! Parameters
    character(len=*),parameter::sbr_nm='wgt_get' ! [sng] Subroutine name
    integer,parameter::wvl_in_nbr_max=1690 ! SWNB resolution
    ! Input
    character(len=*),intent(in)::fl_in
    character(len=*),intent(in)::var_nm
    integer,intent(in)::wvl_out_nbr       ! [nbr] max(wvl_out_nbr)=8192 for H2OH2O.F max(bnd_nbr_TGC98)
    real,intent(in)::wvl_grd_out(wvl_out_nbr+1) 
    ! Output
    real,intent(out)::var_out(wvl_out_nbr)
    ! Local
    integer::wvl_in_nbr        ! dimension size
    integer::grd_in_nbr        ! dimension size
    integer::rcd
    integer::nc_id
    integer::wvl_dim_id
    integer::grd_dim_id
    integer::xtr_typ_LHS
    integer::xtr_typ_RHS
    ! Transmission input variables
    integer::var_in_id
    integer::wvl_grd_in_id
    ! Transmission input variables
    ! Allocatable variables
    real,dimension(:),allocatable::var_in
    real,dimension(:),allocatable::wvl_grd_in 
    ! Local
    integer::out_idx
    
    ! Main Code
    if (dbg_lvl >= dbg_sbr) write (6,'(a)') 'Entering '//sbr_nm
    rcd=nf90_noerr              ! nf90_noerr == 0
    
    ! Initialize default values
    xtr_typ_LHS=xtr_err       !+xtr_vrb
    xtr_typ_RHS=xtr_prt_wgt+xtr_fll_nil !+xtr_vrb (TCG98 goes well past 5 microns)
    
    ! Get Transmission data
    rcd=nf90_wrp_open(fl_in,nf90_nowrite,nc_id)
    ! Get dimension IDs
    rcd=nf90_wrp_inq_dimid(nc_id,'bnd',wvl_dim_id)
    rcd=nf90_wrp_inq_dimid(nc_id,'grd',grd_dim_id)
    ! Get dimension sizes
    rcd=nf90_wrp(nf90_inquire_dimension(nc_id,wvl_dim_id,len=wvl_in_nbr),sbr_nm//": inquire_dim wvl_dim")
    if (wvl_in_nbr > wvl_in_nbr_max) stop 'wvl_in_nbr > wvl_in_nbr_max in wgt_get()'
    rcd=nf90_wrp(nf90_inquire_dimension(nc_id,grd_dim_id,len=grd_in_nbr),sbr_nm//": inquire_dim grd_dim")
    if (grd_in_nbr /= wvl_in_nbr+1) stop 'grd_in_nbr /= wvl_in_nbr+1 in wgt_get()'
    ! Allocate space for dynamic arrays
    allocate(var_in(wvl_in_nbr),stat=rcd)
    allocate(wvl_grd_in(grd_in_nbr),stat=rcd) 
    ! Get variable IDs
    rcd=nf90_wrp_inq_varid(nc_id,'wvl_grd',wvl_grd_in_id)
    rcd=nf90_wrp_inq_varid(nc_id,var_nm,var_in_id)
    ! Get data
    rcd=nf90_wrp(nf90_get_var(nc_id,wvl_grd_in_id,wvl_grd_in),sbr_nm//": get_var wvl_grd_in")
    rcd=nf90_wrp(nf90_get_var(nc_id,var_in_id,var_in),sbr_nm//": get_var var_in")
    ! Close file
    rcd=nf90_wrp_close(nc_id,fl_in,'Ingested')
    
    ! Rebin var_in to wvl_grd_out
    ! rbn_vec() performs any necessary array reversals and monotonicity checks
    call rbn_vec(wvl_in_nbr,wvl_grd_in,var_in, &
         wvl_out_nbr,wvl_grd_out,var_out, &
         xtr_typ_LHS,xtr_typ_RHS)
    
    ! Reset extrapolation flags and strings 
    ! This call is not necessary, but included for safety and defensive programming
    call xtr_ini(xtr_typ_LHS,xtr_typ_RHS)
    
    ! Sanity check
    do out_idx=1,wvl_out_nbr-1
       if (var_out(out_idx) < 0.0) stop 'var_out(out_idx) < 0.0 in wgt_get()'
    end do                    ! end loop over out_idx
    
    ! De-allocate dynamic variables
    if (allocated(var_in)) deallocate(var_in,stat=rcd)
    if (allocated(wvl_grd_in)) deallocate(wvl_grd_in,stat=rcd) 
    
    if (dbg_lvl >= dbg_sbr) write (6,'(a)') 'Exiting '//sbr_nm
    return
  end subroutine wgt_get                       ! end wgt_get()
  
  subroutine wvl_grd_CCM_LW_mk(bnd_nbr,wvl,wvl_min,wvl_max,wvl_ctr,wvl_dlt)
    ! Purpose: Return the CCM SW wavelength grid to the calling program
    ! NB: The CCM SW grid is not regular, so no interface wavelength array is returned
    use wvl_mdl ! [mdl] Wavelength grid parameters
    implicit none
    ! Parameters
    ! Output
    integer,intent(out)::bnd_nbr           ! Number of spectral intervals
    real,intent(out)::wvl(bnd_nbr_CCM_LW_max) ! [m] Nominal wavelength coordinate
    real,intent(out)::wvl_ctr(bnd_nbr_CCM_LW_max) ! [m] Center of band in wavelength space
    real,intent(out)::wvl_min(bnd_nbr_CCM_LW_max) ! [m] Minimum wavelength in band
    real,intent(out)::wvl_max(bnd_nbr_CCM_LW_max) ! [m] Maximum wavelength in band
    real,intent(out)::wvl_dlt(bnd_nbr_CCM_LW_max) ! [m] Bandwidth
    ! Local
    integer::idx               ! Counting index
    real wvn_min(bnd_nbr_CCM_LW_max) ! [cm-1]
    real wvn_max(bnd_nbr_CCM_LW_max) ! [cm-1]
    data wvn_min / 250.0,500.0,650.0,800.0,1000.0,1200.0 / ! [cm-1]
    data wvn_max / 500.0,650.0,800.0,1000.0,1200.0,2200.0 / ! [cm-1]
    ! Main Code
    ! Convert CCM wavelengths to SI
    bnd_nbr=bnd_nbr_CCM_LW_max
    do idx=1,bnd_nbr
       wvl_min(idx)=1.0/(100.0*wvn_max(idx)) ! [cm-1] -> [m]
       wvl_max(idx)=1.0/(100.0*wvn_min(idx)) ! [cm-1] -> [m]
    enddo                     ! end loop over CCM bnd
    do idx=1,bnd_nbr
       wvl_ctr(idx)=0.5*(wvl_min(idx)+wvl_max(idx)) ! [m] 
       wvl_dlt(idx)=wvl_max(idx)-wvl_min(idx) ! [m] 
       wvl(idx)=wvl_ctr(idx)  ! [m] 
    enddo                     ! end loop over CCM bnd
    return 
  end subroutine wvl_grd_CCM_LW_mk                       ! end wvl_grd_CCM_LW_mk()
  
  subroutine wvl_grd_CCM_SW_mk(bnd_nbr,wvl,wvl_min,wvl_max,wvl_ctr,wvl_dlt)
    ! Purpose: Return the CCM SW wavelength grid to the calling program
    ! NB: The CCM SW grid is not regular, so no interface wavelength array is returned
    ! wvl_grd_CCM_SW_mk() is called by nc_out_CCM_SW(), O3(), H2OH2O(), NO2(), and O2X()
    use wvl_mdl ! [mdl] Wavelength grid parameters
    implicit none
    ! Parameters
    ! Output
    integer,intent(out)::bnd_nbr           ! Number of spectral intervals
    real,intent(out)::wvl(bnd_nbr_CCM_SW_max) ! Nominal wavelength coordinate
    real,intent(out)::wvl_ctr(bnd_nbr_CCM_SW_max) ! Center of band in wavelength space
    real,intent(out)::wvl_min(bnd_nbr_CCM_SW_max) ! Minimum wavelength in band
    real,intent(out)::wvl_max(bnd_nbr_CCM_SW_max) ! Maximum wavelength in band
    real,intent(out)::wvl_dlt(bnd_nbr_CCM_SW_max) ! Bandwidth
    ! Local
    integer::idx               ! Counting index
    real(selected_real_kind(p=12)) wvl_min_mcr(bnd_nbr_CCM_SW_max)
    real(selected_real_kind(p=12)) wvl_max_mcr(bnd_nbr_CCM_SW_max)
    ! Setup CCM wavelength grid
    data wvl_min_mcr / &
         0.2000e+00, 0.2450e+00, 0.2650e+00, 0.2750e+00, 0.2850e+00, &
         0.2950e+00, 0.3050e+00, 0.3500e+00, 0.6400e+00, 0.7000e+00, &
         1.0890e+00, 1.4080e+00, 1.7490e+00, 2.1520e+00, 2.7140e+00, &
         3.4570e+00, 2.6300e+00, 4.1600e+00, 4.4750e+00 /
    data wvl_max_mcr / &
         0.2450e+00, 0.2650e+00, 0.2750e+00, 0.2850e+00, 0.2950e+00, &
         0.3050e+00, 0.3500e+00, 0.6400e+00, 0.7000e+00, 1.0890e+00, &
         1.4080e+00, 1.7490e+00, 2.1520e+00, 2.7140e+00, 3.4570e+00, &
         5.0000e+00, 2.8600e+00, 4.4750e+00, 4.5500e+00 /
    ! Main Code
    ! Convert CCM wavelengths to SI
    bnd_nbr=bnd_nbr_CCM_SW_max
    do idx=1,bnd_nbr
       wvl_min(idx)=wvl_min_mcr(idx)*1.0e-6
       wvl_max(idx)=wvl_max_mcr(idx)*1.0e-6
       wvl_ctr(idx)=0.5*(wvl_min(idx)+wvl_max(idx))
       wvl_dlt(idx)=wvl_max(idx)-wvl_min(idx)
       wvl(idx)=wvl_ctr(idx)
    enddo                     ! end loop over CCM bnd
    return 
  end subroutine wvl_grd_CCM_SW_mk                       ! end wvl_grd_CCM_SW_mk()
  
  subroutine nc_out_CCM_SW(fl_out,bnd_dim_id)
    ! Purpose: Write the CCM SW wavelength grid to a netCDF file
    use netcdf ! [mdl] netCDF interface
    use nf90_utl ! [mdl] netCDF utilities
    use wvl_mdl ! [mdl] Wavelength grid parameters
    implicit none
    ! Parameters
    character(len=*),parameter::sbr_nm='nc_out_CCM_SW' ! [sng] Subroutine name
    ! Input
    character(len=*),intent(in)::fl_out      ! Name of netCDF output file
    ! Output
    integer,intent(out)::bnd_dim_id        ! Dimension ID for bnd
    ! Local
    integer::bnd_nbr           ! Number of spectral intervals
    integer::idx               ! Counting index
    integer::nc_id             ! netCDF file ID of output file
    integer::rcd               ! Return success code
    integer::bnd_id            ! Coordinate ID for bnd
    integer::wvl_ctr_id        
    integer::wvl_max_id
    integer::wvl_min_id
    integer::wvl_dlt_id
    integer::wvn_id        
    integer::wvn_ctr_id        
    integer::wvn_max_id
    integer::wvn_min_id
    integer::wvn_dlt_id
    real wvl(bnd_nbr_CCM_SW_max) ! Nominal wavelength coordinate
    real wvl_ctr(bnd_nbr_CCM_SW_max) ! Center of band in wavelength space
    real wvl_min(bnd_nbr_CCM_SW_max) ! Minimum wavelength in band
    real wvl_max(bnd_nbr_CCM_SW_max) ! Maximum wavelength in band
    real wvl_dlt(bnd_nbr_CCM_SW_max) ! Bandwidth
    real wvn(bnd_nbr_CCM_SW_max) ! Nominal wavenumber coordinate
    real wvn_ctr(bnd_nbr_CCM_SW_max) ! Center of band in wavenumber space
    real wvn_min(bnd_nbr_CCM_SW_max) ! Minimum wavenumber in band
    real wvn_max(bnd_nbr_CCM_SW_max) ! Maximum wavenumber in band
    real wvn_dlt(bnd_nbr_CCM_SW_max) ! Bandwidth
    
    ! Main Code
    rcd=nf90_noerr ! nf90_noerr == 0
    bnd_nbr=bnd_nbr_CCM_SW_max
    
    call wvl_grd_CCM_SW_mk(bnd_nbr,wvl,wvl_min,wvl_max,wvl_ctr,wvl_dlt)
    
    do idx=1,bnd_nbr
       wvn_min(idx)=1.0/(100.0*wvl_max(idx)) ! cm-1
       wvn_max(idx)=1.0/(100.0*wvl_min(idx)) ! cm-1
       wvn_dlt(idx)=wvn_max(idx)-wvn_min(idx) ! cm-1
       wvn_ctr(idx)=0.5*(wvn_max(idx)+wvn_min(idx)) ! cm-1
       wvn(idx)=1.0/(100.0*wvl(idx)) ! cm-1
    enddo
    
    rcd=nf90_wrp_open(fl_out,nf90_write,nc_id)
    ! Put output file in define mode 
    rcd=nf90_wrp(nf90_redef(nc_id),sbr_nm//': nf90_redef')
    if (rcd /= nf90_noerr) call nf90_err_exit(rcd,fl_out)
    ! Define dimension IDs
    rcd=nf90_wrp(nf90_def_dim(nc_id,'bnd_CCM',bnd_nbr,bnd_dim_id),sbr_nm//': def_dim bnd_CCM')
    if (rcd /= nf90_noerr) call nf90_err_exit(rcd,fl_out)
    ! Variable definitions
    rcd=nf90_wrp(nf90_def_var(nc_id,'bnd_CCM',nf90_real,bnd_dim_id,bnd_id),sbr_nm//': def_var bnd')
    rcd=nf90_wrp(nf90_def_var(nc_id,'wvl_ctr_CCM',nf90_real,bnd_dim_id,wvl_ctr_id),sbr_nm//': def_var wvl_ctr')
    rcd=nf90_wrp(nf90_def_var(nc_id,'wvl_max_CCM',nf90_real,bnd_dim_id,wvl_max_id),sbr_nm//': def_var wvl_max')
    rcd=nf90_wrp(nf90_def_var(nc_id,'wvl_min_CCM',nf90_real,bnd_dim_id,wvl_min_id),sbr_nm//': def_var wvl_min')
    rcd=nf90_wrp(nf90_def_var(nc_id,'wvl_dlt_CCM',nf90_real,bnd_dim_id,wvl_dlt_id),sbr_nm//': def_var wvl_dlt')
    rcd=nf90_wrp(nf90_def_var(nc_id,'wvn_CCM',nf90_real,bnd_dim_id,wvn_id),sbr_nm//': def_var wvn')
    rcd=nf90_wrp(nf90_def_var(nc_id,'wvn_ctr_CCM',nf90_real,bnd_dim_id,wvn_ctr_id),sbr_nm//': def_var wvn_ctr')
    rcd=nf90_wrp(nf90_def_var(nc_id,'wvn_max_CCM',nf90_real,bnd_dim_id,wvn_max_id),sbr_nm//': def_var wvn_max')
    rcd=nf90_wrp(nf90_def_var(nc_id,'wvn_min_CCM',nf90_real,bnd_dim_id,wvn_min_id),sbr_nm//': def_var wvn_min')
    rcd=nf90_wrp(nf90_def_var(nc_id,'wvn_dlt_CCM',nf90_real,bnd_dim_id,wvn_dlt_id),sbr_nm//': def_var wvn_dlt')
    ! Add global attributes
    ! Add english text descriptions
    rcd=rcd+nf90_put_att(nc_id,bnd_id,'long_name','Band nominal wavelength')
    rcd=rcd+nf90_put_att(nc_id,wvl_ctr_id,'long_name','Band center wavelength')
    rcd=rcd+nf90_put_att(nc_id,wvl_max_id,'long_name','Band maximum wavelength')
    rcd=rcd+nf90_put_att(nc_id,wvl_min_id,'long_name','Band minimum wavelength')
    rcd=rcd+nf90_put_att(nc_id,wvl_dlt_id,'long_name','Bandwidth')
    rcd=rcd+nf90_put_att(nc_id,wvn_id,'long_name','Band nominal wavenumber')
    rcd=rcd+nf90_put_att(nc_id,wvn_ctr_id,'long_name','Band center wavenumber')
    rcd=rcd+nf90_put_att(nc_id,wvn_max_id,'long_name','Band maximum wavenumber')
    rcd=rcd+nf90_put_att(nc_id,wvn_min_id,'long_name','Band minimum wavenumber')
    rcd=rcd+nf90_put_att(nc_id,wvn_dlt_id,'long_name','Bandwidth')
    ! Add units
    rcd=rcd+nf90_put_att(nc_id,bnd_id,'units','meter')
    rcd=rcd+nf90_put_att(nc_id,wvl_ctr_id,'units','meter')
    rcd=rcd+nf90_put_att(nc_id,wvl_ctr_id,'units','meter')
    rcd=rcd+nf90_put_att(nc_id,wvl_max_id,'units','meter')
    rcd=rcd+nf90_put_att(nc_id,wvl_min_id,'units','meter')
    rcd=rcd+nf90_put_att(nc_id,wvl_dlt_id,'units','meter')
    rcd=rcd+nf90_put_att(nc_id,wvn_id,'units','centimeter-1')
    rcd=rcd+nf90_put_att(nc_id,wvn_ctr_id,'units','centimeter-1')
    rcd=rcd+nf90_put_att(nc_id,wvn_max_id,'units','centimeter-1')
    rcd=rcd+nf90_put_att(nc_id,wvn_min_id,'units','centimeter-1')
    rcd=rcd+nf90_put_att(nc_id,wvn_dlt_id,'units','centimeter-1')
    ! All dimensions, variables, and attributes have been defined
    rcd=nf90_wrp(nf90_enddef(nc_id),sbr_nm//': enddef') 
    if (rcd /= nf90_noerr) call nf90_err_exit(rcd,fl_out)
    ! Write out data
    rcd=rcd+nf90_put_var(nc_id,bnd_id,wvl)
    if (rcd /= nf90_noerr) call nf90_err_exit(rcd,'wvl in '//__FILE__)
    rcd=rcd+nf90_put_var(nc_id,wvl_min_id,wvl_min)
    if (rcd /= nf90_noerr) call nf90_err_exit(rcd,'wvl_min in '//__FILE__)
    rcd=rcd+nf90_put_var(nc_id,wvl_max_id,wvl_max)
    if (rcd /= nf90_noerr) call nf90_err_exit(rcd,'wvl_max in '//__FILE__)
    rcd=rcd+nf90_put_var(nc_id,wvl_ctr_id,wvl_ctr)
    if (rcd /= nf90_noerr) call nf90_err_exit(rcd,'wvl_ctr in '//__FILE__)
    rcd=rcd+nf90_put_var(nc_id,wvl_dlt_id,wvl_dlt)
    if (rcd /= nf90_noerr) call nf90_err_exit(rcd,'wvl_dlt in '//__FILE__)
    rcd=rcd+nf90_put_var(nc_id,wvn_id,wvn)
    if (rcd /= nf90_noerr) call nf90_err_exit(rcd,'wvn in '//__FILE__)
    rcd=rcd+nf90_put_var(nc_id,wvn_min_id,wvn_min)
    if (rcd /= nf90_noerr) call nf90_err_exit(rcd,'wvn_min in '//__FILE__)
    rcd=rcd+nf90_put_var(nc_id,wvn_max_id,wvn_max)
    if (rcd /= nf90_noerr) call nf90_err_exit(rcd,'wvn_max in '//__FILE__)
    rcd=rcd+nf90_put_var(nc_id,wvn_ctr_id,wvn_ctr)
    if (rcd /= nf90_noerr) call nf90_err_exit(rcd,'wvn_ctr in '//__FILE__)
    rcd=rcd+nf90_put_var(nc_id,wvn_dlt_id,wvn_dlt)
    if (rcd /= nf90_noerr) call nf90_err_exit(rcd,'wvn_dlt in '//__FILE__)
    ! Close output file
    rcd=nf90_wrp_close(nc_id,fl_out,'Wrote CCM grid to')
    return 
  end subroutine nc_out_CCM_SW                       ! end nc_out_CCM_SW()
  
  subroutine odal_H2OH2O_Chy97(wvn_ctr,lvl_nbr,tpt,mpl_H2O,RH_lqd,odal_H2OH2O,dbg_lvl)
    ! Purpose:
    ! Compute H2O dimer absorption cross section from Petr Chylek (personal communication, 1997) formulation.
    ! Routine computes optical depths for an entire column at a time
    
    ! Wed Jun 18 14:10:51 MDT 1997: Petr wrote to change the RW to NEM values, and correct q q1 value:
    ! "We have now the third set of data which we call "NEM" since they are 
    ! based on the nonlinear equation of motion (denoted by RW previously), that 
    ! include coupling to dimer rotational states. 
    ! When I sent it to you a few weeks ago they were without rotations. 
    ! I have corrected the data in the tables. 
    ! There was also one wrong number in the q=1 set. 
    ! As you can see now the q=1 set and the NEM set are very close to each other and we believe that they are correct.
    ! The main uncertainty is now in the dimer concentration. 
    ! Our equation is a parametric approximation to numerical data with the dimer binding energy of 6.5kcal/mole. 
    ! This value was used in the papers published in the J. Atmos. Chem. and in Appl. Optics. 
    ! If we take 5.5kcal/mole the concentration is by a factor of five lower.
    ! I don't see currently any way how to narrow the range, except from laboratory measurements 
    ! or atmospheric observations, which may take a long time."
    
    ! Wed Jun 18 09:18:53 MDT 1997: Petr wrote to say there were two mistakes in the cross sections:
    ! "1. The third set of cross sections (probably labeled by RW for Reimer and Watts) 
    ! should be multiplied by a factor of two. 
    ! You know that two and pi and similar numbers love to disappear from complicated derivations. 
    ! Now this third set of cross sections and sum rules results with q=1 are quite close to each other. 
    ! 2. As you ask before and I never checked, in the concentration expression 
    ! there should be a factor of RH instead of RH squared."
    
    ! Thu, 29 May 1997 16:11:48 -0300: Petr now has three sets of cross sections:
    ! "We have finished a different method of calculating the absorption cross sections.
    ! We have now three estimates. 
    ! The first one is the one I gave you before. 
    ! We reffer to them as "q=2", "q=1" and "RW".
    ! The bands are the same as before.
    ! q=2: 3.6; 3.6; 2.4; 1.2; 0.4; 0.16
    ! q=1: 6.5; 3.2; 1.1; 0.14; 0.054; 0.009
    ! RW:  3.5; 1.0; 0.4; 0.1; 0.006; 0.0000025
    ! The RW is probably underestimate since it does not include the rotational couplings. 
    ! We expect that the cross sections will be somewhere between q=1 and RW results. 
    ! All cross sections are in 10*(-21) cm*2."
    
    ! This subroutine implements these three estimates with suffices _q2, _q1, and _NEM, respectively. 
    ! Switch between the estimates by setting sigma_cm2 to desired sigma_??_cm2.
    ! The following text described the original formulation (q=2): 
    
    ! "The H2O dimer absorption optical depth is given by
    
    ! tau = 1x10*e19 (1+BT+CT*2) W sigma (RH)*2
    
    ! T is in degrees Celsius, 
    ! W is the total water vapor column in g/cm*2,
    ! RH is the relative humidity with respect to liquid,
    ! sigma is the dimer absorption cross section,
    ! Constants B = 5.62x10*e-2, C = 1.75x10*e-3.
    
    ! The average absorption cross section has been calculated in six broad spectral bands:
    
    ! 500 to 4000cm-1   sigma = 3.6x10*e-21 cm*2
    ! 4000 - 7500cm-1   sigma = 3.6x10*e-21
    ! 7500 -11000cm-1   sigma = 2.4x10*e-21
    ! 11000-14500cm-1   sigma = 1.2x10*e-21
    ! 14500-18000cm-1   sigma = 0.48x10*e-21
    ! 18000-21500cm-1   sigma = 0.16x10*e-21 cm*2
    
    ! The temperature dependence represents the temperature dependence of the
    ! dimer concentration. The cross section's temperature dependence is
    ! currently neglected and the cross sections are for 300K. The accuracy of
    ! the cross sections is estimated to be within a factor of two in infrared
    ! and little bit worse in visible.
    
    ! The expression for tau should be applicable in the for T > -15C. 
    ! For temperatures below -15C the dimer concentration is low and can
    ! be neglected or approximated by the value at -15C. 
    ! The difference between these two approximations should be negligible 
    ! (as far as the absorption is concerned)."
    
    use phys_cst_mdl,only:tpt_frz_pnt ! [mdl] Fundamental and derived physical constants
    implicit none
    ! Parameters
    integer,parameter::max_lvl_nbr=200
    real,parameter::a=1.e19 ! Scale of xsx formula
    real,parameter::b=5.62e-2 ! Linear temperature dependence coefficient in xsx formula
    real,parameter::c=1.75e-3 ! Quadratic temperature dependence coefficient in xsx formula
    ! Input
    integer,intent(in):: &
         dbg_lvl, &
         lvl_nbr              ! dimension size
    real,intent(in):: &
         RH_lqd(lvl_nbr),     & ! fraction
         mpl_H2O(lvl_nbr),    & ! [kg m-2]
         tpt(lvl_nbr),        & ! [K]
         wvn_ctr              ! cm-1
    ! Input/Output
    ! Output
    real,intent(out)::odal_H2OH2O(lvl_nbr) ! [kg m-2]
    ! Local
    integer::idx
    
    real &
         abs_xsx_H2OH2O_cm2, &
         mpl_H2O_gxcm2, &
         npl_H2OH2O(max_lvl_nbr), &
         npl_H2OH2O_cm2, &
         tpt_cls
    ! Main code
    
    do idx=1,lvl_nbr
       
       ! Anything colder than -15 C is treated as -15 C
       tpt_cls=tpt(idx)-tpt_frz_pnt
       if (tpt_cls<-15.0) then 
          tpt_cls=-15.0
       endif
       mpl_H2O_gxcm2=mpl_H2O(idx)*1000./10000. ! [kg m-2] -> g/cm2
       
       npl_H2OH2O(idx)=1.0e4*a*(1.0+b*tpt_cls+c*tpt_cls*tpt_cls)*mpl_H2O_gxcm2*RH_lqd(idx) ! 1.0e4 converts #/cm2 -> #/m2
       npl_H2OH2O_cm2=npl_H2OH2O(idx)*1.0e-4 ! #/m2 -> #/cm2
       abs_xsx_H2OH2O_cm2=1.0e4*abs_xsx_H2OH2O_Chy97(wvn_ctr) ! m2 -> cm2
       odal_H2OH2O(idx)=npl_H2OH2O_cm2*abs_xsx_H2OH2O_cm2
    enddo                     ! end loop over lvl
    
    if (dbg_lvl==5) then
       write (6,'(a,i3)') 'lvl_nbr = ',lvl_nbr
       write (6,'(a3,1x,5(a11,1x))') 'idx','t','mpl_H2O','RH_lqd','npl_H2OH2O','odal_H2OH2O'
       write (6,'(a3,1x,5(a11,1x))') '','K','[kg m-2]','%','#/m2',''
       do idx=1,lvl_nbr
          write (6,'(i3,1x,3(f11.7,1x),1(es11.4,1x),1(f11.7,1x))')  &
               idx,tpt(idx),mpl_H2O(idx),RH_lqd(idx)*100.,npl_H2OH2O(idx),odal_H2OH2O(idx)
       enddo                  ! end loop over lvl
    endif                     ! end if dbg
    
    return
  end subroutine odal_H2OH2O_Chy97                       ! end odal_H2OH2O_Chy97()
  
  real function abs_xsx_H2OH2O_Chy97(wvn_ctr)
    ! Purpose:
    ! Return H2O dimer absorption cross section from Petr Chylek (personal communication, 1997) formulation.
    
    use phys_cst_mdl,only:tpt_frz_pnt ! [mdl] Fundamental and derived physical constants
    implicit none
    ! Parameters
    integer,parameter::max_nbr_H2OH2O_bnd=6
    ! Input
    real,intent(in)::wvn_ctr              ! cm-1
    ! Input/Output
    ! Output
    ! Local
    integer  &
         idx, &
         idx_Chy97, &
         nbr_H2OH2O_bnd
    
    real &
         sigma_NEM_cm2(max_nbr_H2OH2O_bnd), &
         sigma_q1_cm2(max_nbr_H2OH2O_bnd), &
         sigma_q2_cm2(max_nbr_H2OH2O_bnd), &
         sigma_cm2(max_nbr_H2OH2O_bnd), &
         wvn_max_H2OH2O(max_nbr_H2OH2O_bnd), &
         wvn_min_H2OH2O(max_nbr_H2OH2O_bnd)
    
    data wvn_min_H2OH2O /500, & ! 20 microns  
         4000,                & ! 2.5 microns
         7500,                & ! 1.33 microns
         11000,               & ! .91 microns
         14500,               & ! .690 microns
         18000/               ! .555 microns
    data wvn_max_H2OH2O /4000, & ! 2.5 microns
         7500,                & ! 1.33 microns
         11000,               & ! .91 microns
         14500,               & ! .690 microns
         18000,               & ! .555 microns
         21500/               ! .465 microns
    data sigma_NEM_cm2 /7.1e-21, &
         2.1e-21, &
         1.0e-21, &
         0.38e-21, &
         0.08e-21, &
         0.008e-21/
    data sigma_q1_cm2 /6.5e-21, &
         3.2e-21, &
         1.1e-21, &
         0.27e-21, &
         0.05e-21, &
         0.009e-21/
    data sigma_q2_cm2 /3.6e-21, &
         3.6e-21, &
         2.4e-21, &
         1.2e-21, &
         0.48e-21, &
         0.16e-21/
    ! Main code
    
    nbr_H2OH2O_bnd=max_nbr_H2OH2O_bnd
    
    ! Choose among Petr's various experiments
    do idx=1,nbr_H2OH2O_bnd
       ! sigma_cm2(idx)=sigma_NEM_cm2(idx)
       sigma_cm2(idx)=0.5*(sigma_NEM_cm2(idx)+sigma_q1_cm2(idx))
    end do                    ! end loop over bnd
    
    idx_chy97=-1 ! CEWI
    if (wvn_ctr<wvn_min_H2OH2O(1)) then
       idx_Chy97=1
    else if (wvn_ctr>=wvn_max_H2OH2O(nbr_H2OH2O_bnd)) then
       idx_Chy97=nbr_H2OH2O_bnd
    else
       do idx=1,nbr_H2OH2O_bnd
          if (wvn_ctr>=wvn_min_H2OH2O(idx).and.wvn_ctr<wvn_max_H2OH2O(idx)) then 
             idx_Chy97=idx
             goto 300
          endif
       end do                 ! end loop over bnd
    endif
300 continue
    
    abs_xsx_H2OH2O_Chy97=1.0e-4*sigma_cm2(idx_Chy97) ! cm2 -> m2
    
    return
  end function abs_xsx_H2OH2O_Chy97                       ! end abs_xsx_H2OH2O_Chy97()
  
  complex function Voigt_Hum82(z)
    
    ! Purpose: Compute complex probability function using method of Hum82
    ! See Lio92 p. 30 (2.2.11) for derivation of arguments
    
    ! The Voigt function used in spectroscopy and astronomy is the real component of the complex error function
    ! See HAW78 for further details
    
    ! "Computes the complex probability function w(z)=exp(-z**2)*erfc(-i*z)
    ! in the upper half-plane z=x+iy (i.e., for y >= 0).
    ! Maximum relative error of both real and imaginary parts is < 1.0e-4"
    
    implicit none
    ! Parameters
    ! Input
    complex,intent(in)::z
    ! Input/Output
    ! Output
    ! Local
    complex t
    complex u
    real s
    real x
    real y
    ! Main code
    x=real(z)
    y=aimag(z)
    t=cmplx(y,-x)
    s=abs(x)+y
    
    if(s<15.0) goto 1
    
    ! Region 1
    Voigt_Hum82=t*0.5641896/(0.5+t*t)
    return
    
1   if(s<5.5) goto 2
    
    ! Region 2
    u=t*t
    Voigt_Hum82=t*(1.410474+u*0.5641896)/(0.75+u*(3.0+u))
    return
    
2   if(y<0.195*abs(x)-0.176) goto 3
    
    ! Region 3
    Voigt_Hum82=(16.4955+t*(20.20933+t*(11.96482+t*(3.778987+t*0.5642236))))/ &
         (16.4955+t*(38.82363+t*(39.27121+t*(21.69274+t*(6.699398+t)))))
    return
    
    ! Region 4
3   u=t*t
    Voigt_Hum82=cexp(u)-t*(36183.31-u*(3321.9905-u*(1540.787-u*(219.0313-u* &
         (35.76683-u*(1.320522-u*0.56419))))))/(32066.6-u*(24322.84-u* &
         (9022.228-u*(2186.181-u*(364.2191-u*(61.57037-u*(1.841439-u)))))))
    
    return
  end function Voigt_Hum82                       ! end voigt_Hum82()
  
  subroutine voigt_Kun97(x,y,nx,fac,prb)
    ! Purpose: Compute Voigt function using method of Kun97
    ! See Lio92 p. 30 (2.2.11) for derivation of arguments
    
    ! The Voigt function used in spectroscopy and astronomy is the real component of the complex error function
    ! See HAW78 for further details
    ! Kun97 accelerates Humlicek's method (Hum82) by only computing the real component of the complex error function
    ! This avoids computing the imaginary part, which is always thrown away in LBL applications
    ! Kun97 claims to be a factor of ~3 faster than Hum82 while keeping errors relative to Hum82 below 2.0e-6
    ! Following is the Kun97 subroutine as sent to me by Martin Kuntz on 19971008:
    
    ! Calculates the Voigt-Function times a user-definied factor 
    ! fac with a relative accuracy better than 2*10-6.
    
    ! If this subroutine is called several times with similar 
    ! values y, the numerically expensive coefficents a1..t8 
    ! are only calculated once. The coefficients are only recalculated
    ! if the relative change in y is grater than the internal parameter 
    ! rel, which is set to 1e-4. 
    
    ! x(*)   (in)   : Distance from line center in units of Doppler
    !               : halfwidths times sqrt(ln 2)       
    ! y      (in)   : Ratio of the Doppler halfwidth to the Lorentz
    !               : halfwidth times sqrt(ln2)	
    ! nx     (in)   : Length of vector x and prb
    ! fac    (in)   : A user defined factor
    ! prb(*) (out)  : The Voigt-function times the user defined factor
    
    ! author: Dr. Martin Kuntz,
    !         Institut fuer Meteorologie und Klimaforschung,
    !         Forschungszentrum Karlsruhe,
    !         Postfach 3640,
    !         76021 Karlsruhe, Germany
    !         kuntz@imk.fzk.de
    
    ! ref:    A new implementation of the Humlicek algorithm for the 
    !         calculation of the Voigt profile function, J. Quant. 
    !         Spectrosc. Radiat. Transfer, 57, 819-824, 1997 
    implicit none
    ! Parameters
    real rel
    parameter (rel = 1.0e-4)
    ! Input
    integer,intent(in)::nx
    real,intent(in)::fac
    real,intent(in)::x(nx)
    real,intent(in)::y
    ! Input/Output
    ! Output
    real,intent(out)::prb(nx)
    ! Local
    integer bmax
    integer bmin 
    integer lauf(1:4,4) 
    integer stack(20,4) 
    integer stackp 
    
    real a8,b8,c8,d8,e8,f8,g8,h8,o8,p8,q8,r8,s8,t8,a7,b7,c7,d7,e7,f7, &
         g7,h7,o7,p7,q7,r7,s7,t7,a6,b6,c6,d6,e6,a5,b5,c5,d5,e5,a4,b4, &
         c4,d4,a3,b3,c3,d3,a2,a1,b1,yps4,yps3,yps2,yps1
    save a8,b8,c8,d8,e8,f8,g8,h8,o8,p8,q8,r8,s8,t8,a7,b7,c7,d7,e7,f7, &
         g7,h7,o7,p7,q7,r7,s7,t7,a6,b6,c6,d6,e6,a5,b5,c5,d5,e5,a4,b4, &
         c4,d4,a3,b3,c3,d3,a2,a1,b1,yps4,yps3,yps2,yps1
    data yps1,yps2,yps3,yps4 /-1.0,-1.0,-1.0,-1.0/
    
    integer i1
    integer i2
    integer imax
    integer imin
    integer imitte
    real b2
    real c1
    real x2
    real y2
    real ym2
    ! Main code
    y2 = y*y
    if (y>=15.0.or.x(1)>=15.0.or.x(nx)<=-15.0) then
       lauf(1,1) = 1
       lauf(1,2) = nx
       lauf(1,3) = nx
       lauf(1,4) = 0
       goto 7
    endif
    
    do i2 = 1,4
       do i1 = 1,4
          lauf(i1,i2) = mod(i2,2)*(nx+1)
       end do ! i1
    end do ! i2
    
    stackp = 1
    stack(stackp,1) = 1
    stack(stackp,2) = nx
    stack(stackp,3) = bfun(y,x(1))
    stack(stackp,4) = bfun(y,x(nx))
    
2   imin = stack(stackp,1)
    imax = stack(stackp,2)
    bmin = stack(stackp,3)
    bmax = stack(stackp,4)
    if (bmin==bmax) then
       if (x(imax)<0.0) then
          lauf(bmin,1) = min(imin,lauf(bmin,1))
          lauf(bmax,2) = max(imax,lauf(bmax,2))
          stackp = stackp-1
          goto 3
       elseif (x(imin)>=0.0) then
          lauf(bmin,3) = min(imin,lauf(bmin,3))
          lauf(bmax,4) = max(imax,lauf(bmax,4))
          stackp = stackp-1
          goto 3
       endif
    endif
    imitte = (imax+imin)/2
    stack(stackp,1) = imitte+1
    stack(stackp,2) = imax
    stack(stackp,3) = bfun(y,x(imitte+1))
    stack(stackp,4) = bmax
    stackp = stackp+1
    stack(stackp,1) = imin
    stack(stackp,2) = imitte
    stack(stackp,3) = bmin
    stack(stackp,4) = bfun(y,x(imitte))
3   if (stackp>0) goto 2
    
    ! Region 4
    if (lauf(4,2)>=lauf(4,1).or.lauf(4,4)>=lauf(4,3)) then
       if (abs((y-yps4)/yps4)>rel) then
          yps4 = y
          a7 = y*(1.16028e9+y2*(-9.86604e8+y2*(4.56662e8+y2* &
               (-1.53575e8+y2*(4.08168e7+y2*(- 9.69463e6+y2*(1.6841e6+y2* &
               (-320772.0+y2*(40649.2+y2*(-5860.68+y2*(571.687+y2*(-72.9359 &
               +y2*(2.35944-y2*0.56419)))))))))))))
          b7 = y*(-5.60505e8+y2*(-9.85386e8+y2*(8.06985e8+y2* &
               (-2.91876e8+y2*(8.64829e7+y2*(-7.72359e6+y2*(3.59915e6+y2* &
               (-234417.0+y2*(45251.3+y2*(-2269.19+y2*(-234.143+y2* &
               (23.0312-y2*7.33447))))))))))))
          c7 = y*(-6.51523e8+y2*(2.47157e8+y2*(2.94262e8+y2* &
               (-2.04467e8+y2*(2.29302e7+y2*(-2.3818e7+y2*(576054.0+y2* &
               (98079.1+y2*(-25338.3+y2*(1097.77+y2* &
               (97.6203-y2*44.0068)))))))))))
          d7 = y*(-2.63894e8+y2*(2.70167e8+y2*(-9.96224e7+y2* &
               (-4.15013e7+y2*(3.83112e7+y2*(2.2404e6+y2*(-303569.0+y2* &
               (-66431.2+y2*(8381.97+y2*(228.563-y2*161.358))))))))))
          e7 = y*(-6.31771e7+y2*(1.40677e8+y2*(5.56965e6+y2* &
               (2.46201e7+y2*(468142.0+y2*(-1.003e6+y2*(-66212.1+y2* &
               (23507.6+y2*(296.38-y2*403.396)))))))))
          f7 = y*(-1.69846e7+y2*(4.07382e6+y2*(-3.32896e7+y2* &
               (-1.93114e6+y2*(-934717.0+y2*(8820.94+y2*(37544.8+y2* &
               (125.591-y2*726.113))))))))
          g7 = y*(-1.23165e6+y2*(7.52883e6+y2*(-900010.0+y2*(-186682.0+ &
               y2*(79902.5+y2*(37371.9+y2*(-260.198-y2*968.15)))))))
          h7 = y*(-610622.0+y2*(86407.6+y2*(153468.0+y2*(72520.9+y2* &
               (23137.1+y2*(-571.645-y2*968.15))))))
          o7 = y*(-23586.5+y2*(49883.8+y2*(26538.5+y2*(8073.15+y2* &
               (-575.164-y2*726.113)))))
          p7 = y*(-8009.1+y2*(2198.86+y2*(953.655+y2* &
               (-352.467-y2*403.396))))
          q7 = y*(-622.056+y2*(-271.202+y2*(-134.792-y2*161.358)))
          r7 = y*(-77.0535+y2*(-29.7896-y2*44.0068))
          s7 = y*(-2.92264-y2*7.33447)
          t7 = y*(-0.56419)
          a8 = 1.02827e9+y2*(-1.5599e9+y2*(1.17022e9+y2*(-5.79099e8+y2* &
               (2.11107e8+y2*(-6.11148e7+y2*(1.44647e7+y2*(-2.85721e6+y2* &
               (483737.0+y2*(-70946.1+y2*(9504.65+y2*(-955.194+y2*(126.532 &
               +y2*(-3.68288+y2)))))))))))))
          b8 = 1.5599e9+y2*(-2.28855e9+y2*(1.66421e9+y2*(-7.53828e8+y2* &
               (2.89676e8+y2*(-7.01358e7+y2*(1.39465e7+y2*(-2.84954e6+y2* &
               (498334.0+y2*(-55600.0+y2*(3058.26+y2*(533.254+y2* &
               (-40.5117+y2*14.0))))))))))))
          c8 = 1.17022e9+y2*(-1.66421e9+y2*(1.06002e9+y2*(-6.60078e8+y2* &
               (6.33496e7+y2*(-4.60396e7+y2*(1.4841e7+y2*(-1.06352e6+y2* &
               (-217801.0+y2*(48153.3+y2*(-1500.17+y2* &
               (-198.876+y2*91)))))))))))
          d8 = 5.79099e8+y2*(-7.53828e8+y2*(6.60078e8+y2*(5.40367e7+y2* &
               (1.99846e8+y2*(-6.87656e6+y2*(-6.89002e6+y2*(280428.0+y2* &
               (161461.0+y2*(-16493.7+y2*(-567.164+y2*364))))))))))
          e8 = 2.11107e8+y2*(-2.89676e8+y2*(6.33496e7+y2*(-1.99846e8+y2* &
               (-5.01017e7+y2*(-5.25722e6+y2*(1.9547e6+y2*(240373.0+y2* &
               (-55582.0+y2*(-1012.79+y2*1001)))))))))
          f8 = 6.11148e7+y2*(-7.01358e7+y2*(4.60396e7+y2*(-6.87656e6+y2* &
               (5.25722e6+y2*(3.04316e6+y2*(123052.0+y2*(-106663.0+y2* &
               (-1093.82+y2*2002))))))))
          g8 = 1.44647e7+y2*(-1.39465e7+y2*(1.4841e7+y2* &
               (6.89002e6+y2*(1.9547e6+y2*(-123052.0+y2*(-131337.0+y2* &
               (-486.14+y2*3003)))))))
          h8 = 2.85721e6+y2*(-2.84954e6+y2*(1.06352e6+y2*(280428.0+y2* &
               (-240373.0+y2*(-106663.0+y2*(486.14+y2*3432))))))
          o8 = 483737.0+y2*(-498334.0+y2*(-217801.0+y2*(-161461.0+y2* &
               (-55582.0+y2*(1093.82+y2*3003)))))
          p8 = 70946.1+y2*(-55600.0+y2*(-48153.3+y2*(-16493.7+y2* &
               (1012.79+y2*2002))))
          q8 = 9504.65+y2*(-3058.26+y2*(-1500.17+y2*(567.164+y2*1001.0)))
          r8 = 955.194+y2*(533.254+y2*(198.876+y2*364))
          s8 = 126.532+y2*(40.5117+y2*91.0)
          t8 = 3.68288+y2*14.
       endif
       ym2 = y*2
       do i2 = 1,3,2
          do i1 = lauf(4,i2),lauf(4,i2+1)
             x2 = x(i1)*x(i1)
             prb(i1) = fac*(exp(y2-x2)*cos(x(i1)*ym2) - &
                  (a7+x2*(b7+x2*(c7+x2*(d7+x2*(e7+x2*(f7+x2*(g7+x2*(h7+x2* &
                  (o7+x2*(p7+x2*(q7+x2*(r7+x2*(s7+x2*t7))))))))))))) / &
                  (a8+x2*(b8+x2*(c8+x2*(d8+x2*(e8+x2*(f8+x2*(g8+x2*(h8+x2* &
                  (o8+x2*(p8+x2*(q8+x2*(r8+x2*(s8+x2*(t8+x2)))))))))))))))
          end do ! i1
       end do ! i2
    endif
    
    ! Region 3
    if (lauf(3,2)>=lauf(3,1).or.lauf(3,4)>=lauf(3,3)) then
       if (abs((y-yps3)/yps3)>rel) then
          yps3 = y
          a5 = (272.102+y*(973.778+y*(1629.76+y*(1678.33+y*(1174.8+y* &
               (581.746+y*(204.501+y*(49.5213+y*(7.55895+y*0.564224)))))))))
          b5 = (-60.5644+y*(-2.34403+y*(220.843+y*(336.364+y*(247.198 &
               +y*(100.705+y*(22.6778+y*2.25689)))))))
          c5 = (4.58029+y*(18.546+y*(42.5683+y*(52.8454+y*(22.6798+y* &
               3.38534)))))
          d5 = (-0.128922+y*(1.66203+y*(7.56186+y*2.25689)))
          e5 = (0.000971457+y*0.564224)
          a6 = 272.102+y*(1280.83+y*(2802.87+y*(3764.97+y*(3447.63+y* &
               (2256.98+y*(1074.41+y*(369.199+y*(88.2674+ &
               y*(13.3988+y)))))))))
          b6 = 211.678+y*(902.306+y*(1758.34+y*(2037.31+y*(1549.68+y* &
               (793.427+y*(266.299+y*(53.5952+y*5.0)))))))
          c6 = 78.866+y*(308.186+y*(497.302+y*(479.258+y*(269.292+y* &
               (80.3928+y*10.0)))))
          d6 = 22.0353+y*(55.0293+y*(92.7568+y*(53.5952+y*10.0)))
          e6 = 1.49645+y*(13.3988+y*5.0)
       endif
       do i2 = 1,3,2
          do i1 = lauf(3,i2),lauf(3,i2+1)
             x2 = x(i1)*x(i1)
             prb(i1) = fac*(a5+x2*(b5+x2*(c5+x2*(d5+x2*e5))))/ &
                  (a6+x2*(b6+x2*(c6+x2*(d6+x2*(e6+x2)))))
          end do ! i1
       end do ! i2
    endif
    
    ! Region 2
    if (lauf(2,2)>=lauf(2,1).or.lauf(2,4)>=lauf(2,3)) then
       if (abs((y-yps2)/yps2)>rel) then
          yps2 = y
          a3 = y*(1.05786+y2*(4.65456+y2*(3.10304+y2*0.56419)))
          b3 = y*(2.962+y2*(0.56419+y2*1.69257))
          c3 = y*(1.69257*y2-2.53885)
          d3 = y*(0.56419)
          a4 = 0.5625+y2*(4.5+y2*(10.5+y2*(6.0+y2)))
          b4 = -4.5+y2*(9.0+y2*(6.0+y2*4.0))
          c4 = 10.5+y2*(-6.0+y2*6.0)
          d4 = -6.0+y2*4.0
       endif
       do i2 = 1,3,2
          do i1 = lauf(2,i2),lauf(2,i2+1)
             x2 = x(i1)*x(i1)
             prb(i1) = fac*(a3+x2*(b3+x2*(c3+x2*d3)))/ &
                  (a4+x2*(b4+x2*(c4+x2*(d4+x2))))
          end do ! i1
       end do ! i2
    endif
    
    ! Region 1
7   if (lauf(1,2)>=lauf(1,1).or.lauf(1,4)>=lauf(1,3)) then
       if (abs((y-yps1)/yps1)>rel) then
          yps1 = y
          a1 = 0.5641896*y
          b1 = 0.5+y2
          a2 = 4*y2
       endif
       c1 = fac*a1
       do i2 = 1,3,2
          do i1 = lauf(1,i2),lauf(1,i2+1)
             x2 = x(i1)*x(i1)
             b2 = b1-x2
             prb(i1) = c1*(b1+x2)/(b2*b2+a2*x2)
          end do ! i1
       end do ! i2
    endif
    
    return
  end subroutine voigt_Kun97                        ! end voigt_Kun97()
  
  integer function bfun(y,x)
    ! bfun() is supplied as a function utilized by voigt_Kun97()
    implicit none
    ! Input
    real,intent(in)::x
    real,intent(in)::y
    ! Local
    real s
    ! Main code
    s = abs(x)+y
    if (s>=15) then
       bfun = 1
    elseif (s>=5.5) then
       bfun = 2
    elseif (y>=(0.195*abs(x))-0.176) then
       bfun = 3
    else
       bfun = 4
    endif
    return
  end function bfun                       ! end bfun()
  
  character(10) function htrn_mlc_sng(mlc_nbr)
    ! Purpose: Return string describing molecule based on HITRAN96 ordering 
    implicit none
    ! Parameters
    integer,parameter::mlc_nbr_max_htrn=37 ! HITRAN96
    ! Input
    integer,intent(in)::mlc_nbr
    ! Input/Output
    ! Output
    ! Local
    character(10)::mlc_sng(mlc_nbr_max_htrn)
    ! Main code
    mlc_sng(1)='H2O'
    mlc_sng(2)='CO2'
    mlc_sng(3)='O3'
    mlc_sng(4)='N2O'
    mlc_sng(5)='CO'
    mlc_sng(6)='CH4'
    mlc_sng(7)='O2'
    mlc_sng(8)='NO'
    mlc_sng(9)='SO2'
    mlc_sng(10)='NO2'
    mlc_sng(11)='NH3'
    mlc_sng(12)='HNO3'
    mlc_sng(13)='OH'
    mlc_sng(14)='HF'
    mlc_sng(15)='HCl'
    mlc_sng(16)='HBr'
    mlc_sng(17)='HI'
    mlc_sng(18)='ClO'
    mlc_sng(19)='OCS'
    mlc_sng(20)='H2CO'
    mlc_sng(21)='HOCl'
    mlc_sng(22)='N2'
    mlc_sng(23)='HCN'
    mlc_sng(24)='CH3Cl'
    mlc_sng(25)='H2O2'
    mlc_sng(26)='C2H2'
    mlc_sng(27)='C2H6'
    mlc_sng(28)='PH3'
    mlc_sng(29)='COF2'
    mlc_sng(30)='SF6'
    mlc_sng(31)='H2S'
    mlc_sng(32)='HCOOH'
    mlc_sng(33)='HO2' 
    mlc_sng(34)='O' 
    mlc_sng(35)='ClONO2'
    mlc_sng(36)='NO+' 
    mlc_sng(37)='HOBr'
    htrn_mlc_sng=mlc_sng(mlc_nbr)
    return
  end function htrn_mlc_sng                       ! end htrn_mlc_sng()
  
  function lorentz(HWHM,ln_ctr,wvn)
    ! Purpose:
    ! Return normalized shape factor of Lorentzian at specified distance from line center
    implicit none
    ! Parameters
    ! Input
    real(selected_real_kind(p=12)),intent(in)::HWHM     ! 
    real(selected_real_kind(p=12)),intent(in)::ln_ctr   ! 
    real(selected_real_kind(p=12)),intent(in)::wvn      ! 
    ! Output
    ! Local
    real(selected_real_kind(p=12))::lorentz
    integer idx
    ! Main code
    idx=1
    ! fxm: function is unfinished
    lorentz=HWHM+0.0*wvn+0.0*ln_ctr+idx*0.0
    return
  end function lorentz                       ! end lorentz()
  
  subroutine ln_nbr_get(ln_ctr,ln_nbr,wng_sz_wvn,ln_idx_srt,ln_idx_end)
    ! Purpose:
    ! Return indices which bracket all lines in input array within a specified 
    ! distance from the given line center
    implicit none
    ! Parameters
    ! Input
    real(selected_real_kind(p=12)),intent(in)::ln_ctr
    real(selected_real_kind(p=12)),intent(in)::wng_sz_wvn
    integer,intent(in)::ln_nbr
    ! Input/Output
    ! Output
    integer,intent(out)::ln_idx_end
    integer,intent(out)::ln_idx_srt
    ! Local
    integer idx
    ! Main code
    idx=1
    ! fxm: function is unfinished
    ln_idx_srt=idx+0*ln_nbr+int(0.0*ln_ctr+0.0*wng_sz_wvn) ! CEWI
    ln_idx_end=idx
    return 
  end subroutine ln_nbr_get                       ! end ln_nbr_get()
  
  subroutine wvl_grd_mk(wvl_grd_typ, & ! Input
       wvl_grd_min,wvl_grd_max,wvl_grd_rsn, & ! Input/Output
       wvl_ctr,wvl_dlt,wvl_min,wvl_max,wvl_grd,wvl_nbr) ! Output
    ! Purpose: Create and return wavelength grid specified by wvl_grd_typ
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use sng_mdl,only:ftn_strlen   ! [mdl] String manipulation
    use utl_mdl ! [mdl] Utility functions (date_time_get,mnt_chk...)
    use wvl_mdl ! [mdl] Wavelength grid parameters
    implicit none
    ! Parameters
    integer,parameter::wvl_nbr_max=1690 ! [nbr] Maximum number of wavelengths
    ! Input
    integer,intent(in)::wvl_grd_typ       ! [enm] Wavelength grid type
    ! Input/Output
    real,intent(inout)::wvl_grd_min          ! [m] Minimum wavelength
    real,intent(inout)::wvl_grd_max          ! [m] Maximum wavelength
    real,intent(inout)::wvl_grd_rsn          ! [m] Resolution
    ! Output
    ! fxm: 20020520 wvl_nbr cannot be intent out _and_ used to allocate arrays
    ! so temporarily make it intent(inout) and recode this procedure before using it
    integer,intent(inout)::wvl_nbr            ! Number of spectral intervals
    real,intent(out)::wvl_grd(wvl_nbr+1)   ! [m] Wavelength grid
    real,intent(out)::wvl_ctr(wvl_nbr)     ! [m] Center of band in wavelength space
    real,intent(out)::wvl_min(wvl_nbr)     ! [m] Minimum wavelength in band
    real,intent(out)::wvl_max(wvl_nbr)     ! [m] Maximum wavelength in band
    real,intent(out)::wvl_dlt(wvl_nbr)     ! [m] Bandwidth
    ! Local
    integer wvl_idx           ! Counting index
    integer nbr_pure_O3_bnd   ! Size used in wvl_grd_NBM_SW
    logical mnt_ncr           ! Monotonic and increasing flag
    real wvn_max              ! [cm-1] Maximum wavenumber in bin
    ! Main Code
    ! Initialize scalars
    ! Define wavelength grid
    ! Procedure is to define wvl_nbr, wvl_min, wvl_grd_min, wvl_grd_max as primary variables
    ! All other variables are subsequently derived from these 
    if (wvl_grd_typ==wvl_grd_CCM_SW) then
       ! Define wavelength grid approximately equal to CCM shortwave grid
       call wvl_grd_CCM_SW_mk(wvl_nbr,wvl_ctr,wvl_min,wvl_max,wvl_ctr,wvl_dlt)
    else if (wvl_grd_typ==wvl_grd_CCM_LW) then
       ! Define wavelength grid approximately equal to CCM longwave grid
       call wvl_grd_CCM_LW_mk(wvl_nbr,wvl_ctr,wvl_min,wvl_max,wvl_ctr,wvl_dlt)
    else if (wvl_grd_typ==wvl_grd_NBM_SW) then
       ! Define wavelength grid approximately equal to original NBM grid
       ! 1590 10 cm-1 bands from 2000--17900 cm-1
       ! 100 bands from 176--558.65 nm
       wvl_grd_min=175.4e-09  ! [m]
       wvl_grd_max=5.0e-06 ! [m]
       wvl_nbr=1690
       nbr_pure_O3_bnd=100
       do wvl_idx=1,wvl_nbr-nbr_pure_O3_bnd
          wvn_max=2000.0+wvl_idx*10 ! [cm-1]
          wvl_min(wvl_idx)=1.0/(100.0*wvn_max) ! [m]
       enddo                  ! end loop over wvl
       wvl_grd_rsn=(wvl_min(wvl_nbr-nbr_pure_O3_bnd)-wvl_grd_min)/nbr_pure_O3_bnd
       do wvl_idx=wvl_nbr,wvl_nbr-(nbr_pure_O3_bnd-1),-1
          wvl_min(wvl_idx)=wvl_grd_min+wvl_grd_rsn*(wvl_nbr-wvl_idx) ! [m]
       enddo                  ! end loop over wvl
    else if (wvl_grd_typ==wvl_grd_NBM_LW) then
       ! Define wavelength grid appropriate for longwave
       wvl_grd_min=5.0e-06    ! [m]
       wvl_grd_max=100.0e-06   ! [m]
       wvl_grd_rsn=1.0e-06    ! [m]
       wvl_nbr=nint((wvl_grd_max-wvl_grd_min)/wvl_grd_rsn)
       do wvl_idx=1,wvl_nbr
          wvl_min(wvl_idx)=wvl_grd_min+(wvl_idx-1)*wvl_grd_rsn ! [m]
       enddo                  ! end loop over wvl
    else if (wvl_grd_typ==wvl_grd_dfl) then
       ! Define default wavelength grid
       wvl_nbr=nint((wvl_grd_max-wvl_grd_min)/wvl_grd_rsn)
       do wvl_idx=1,wvl_nbr
          wvl_min(wvl_idx)=wvl_grd_min+(wvl_idx-1)*wvl_grd_rsn ! [m]
       enddo                  ! end loop over wvl
    else         
       stop 'wvl_grd_mk() reports unknown wvl_grd_typ'
    endif                     ! endelse
    
    ! Following derived grid values depend upon direction of monotonicity
    mnt_ncr=mnt_ncr_chk(wvl_min,wvl_nbr)
    ! Compute wvl_max then wvl_grd
    if (mnt_ncr) then
       do wvl_idx=1,wvl_nbr-1 ! Loop ends at wvl_nbr-1
          wvl_max(wvl_idx)=wvl_min(wvl_idx+1) ! [m]
       enddo                  ! end loop over wvl
       wvl_max(wvl_nbr)=wvl_grd_max ! [m]
       do wvl_idx=1,wvl_nbr
          wvl_grd(wvl_idx)=wvl_min(wvl_idx) ! [m]
       enddo                  ! end loop over wvl
       wvl_grd(wvl_nbr+1)=wvl_max(wvl_nbr)
       wvl_grd_max=wvl_max(wvl_nbr) ! [m]
       wvl_grd_min=wvl_min(1) ! [m]
    else                      ! not mnt_ncr
       do wvl_idx=2,wvl_nbr   ! Loop starts at 2
          wvl_max(wvl_idx)=wvl_min(wvl_idx-1) ! [m]
       enddo                  ! end loop over wvl
       wvl_max(1)=wvl_grd_max ! [m]
       do wvl_idx=1,wvl_nbr
          wvl_grd(wvl_idx)=wvl_max(wvl_idx) ! [m]
       enddo                  ! end loop over wvl
       wvl_grd(wvl_nbr+1)=wvl_min(wvl_nbr)
       wvl_grd_max=wvl_max(1) ! [m]
       wvl_grd_min=wvl_min(wvl_nbr) ! [m]
    endif                     ! not mnt_ncr
    wvl_grd_rsn=(wvl_grd_max-wvl_grd_min)/wvl_nbr
    
    ! Final derived grid values 
    do wvl_idx=1,wvl_nbr
       wvl_ctr(wvl_idx)=0.5*(wvl_min(wvl_idx)+wvl_max(wvl_idx)) ! [m]
       wvl_dlt(wvl_idx)=wvl_max(wvl_idx)-wvl_min(wvl_idx) ! [m]
    enddo                     ! end loop over wvl
    
    ! Sanity check: even boundaries?
    if (mnt_ncr) then
       do wvl_idx=1,wvl_nbr-1 ! Loop ends at wvl_nbr-1
          if (wvl_min(wvl_idx+1) /= wvl_max(wvl_idx)) then
             write (6,'(a,a)') prg_nm(1:ftn_strlen(prg_nm)), &
                  ': ERROR wvl_grd_mk(): Uneven boundaries monotonically increasing array' 
             write (6,'(2(a,i4,a,es10.3))')   &
                  'wvl_min(',wvl_idx+1,') = ',wvl_min(wvl_idx+1), &
                  ' != wvl_max(',wvl_idx,') = ',wvl_max(wvl_idx)
             stop
          endif               ! endif
       end do                 ! end loop over wvl
    else                      ! not mnt_ncr
       do wvl_idx=1,wvl_nbr-1 ! Loop ends at wvl_nbr-1
          if (wvl_min(wvl_idx) /= wvl_max(wvl_idx+1)) then
             write (6,'(a,a)') prg_nm(1:ftn_strlen(prg_nm)), &
                  ': ERROR wvl_grd_mk(): Uneven boundaries monotonically decreasing array' 
             write (6,'(2(a,i4,a,es10.3))')   &
                  'wvl_min(',wvl_idx,') = ',wvl_min(wvl_idx), &
                  ' != wvl_max(',wvl_idx+1,') = ',wvl_max(wvl_idx+1)
             stop
          endif               ! endif
       end do                 ! end loop over wvl
    endif                     ! not mnt_ncr
    
    if (dbg_lvl>dbg_scl) then
       write (6,'(2a,i4)') prg_nm(1:ftn_strlen(prg_nm)),': wvl_grd_mk() set wvl_nbr = ',wvl_nbr
    endif                     ! endif dbg
    if (dbg_lvl>=dbg_io) then
       ! Print wavelength grid
       write (6,'(2a)') prg_nm(1:ftn_strlen(prg_nm)),': DEBUG Grid output by wvl_grd_mk():'
       write (6,'(6(a,1x))') 'wvl_idx','wvl_dlt','wvl_ctr','wvl_min','wvl_max','wvl_grd'
       do wvl_idx=1,wvl_nbr
          write (6,'(i4,1x,5(es15.8,1x))')  &
               wvl_idx,wvl_dlt(wvl_idx),wvl_ctr(wvl_idx),wvl_min(wvl_idx),wvl_max(wvl_idx),wvl_grd(wvl_idx)
       end do                 ! end loop over wvl
       write (6,'(i4,1x,5(es15.8,1x))') wvl_idx,1.0e36,1.0e36,1.0e36,1.0e36,wvl_grd(wvl_idx)
    endif                     ! endif dbg
    
    ! Sanity check: Everything else
    do wvl_idx=1,wvl_nbr
       ! For some reason this check triggers quite often with g77
       ! if (wvl_ctr(wvl_idx) /= (0.5*(wvl_min(wvl_idx)+wvl_max(wvl_idx)))) write (6,'(a,i4)') 'Bad wvl_ctr: ',wvl_idx
       if (wvl_dlt(wvl_idx) /= (wvl_max(wvl_idx)-wvl_min(wvl_idx))) write (6,'(a,i4)')  &
            'wvl_grd_mk() reports bad wvl_dlt: ',wvl_idx
       if (wvl_max(wvl_idx)<0.0) write (6,'(a,i4)')  &
            'wvl_grd_mk() reports bad wvl_max: ',wvl_idx
       if (wvl_min(wvl_idx)<0.0) write (6,'(a,i4)')  &
            'wvl_grd_mk() reports bad wvl_min: ',wvl_idx
       if (wvl_dlt(wvl_idx)<0.0) write (6,'(a,i4)')  &
            'wvl_grd_mk() reports bad wvl_dlt: ',wvl_idx
       if (wvl_ctr(wvl_idx)<0.0) write (6,'(a,i4)')  &
            'wvl_grd_mk() reports bad wvl_ctr: ',wvl_idx
       ! write (6,'(a,a,i4)') prg_nm(1:ftn_strlen(prg_nm)),': ERROR wvl_grd_mk(): Bad wavelength bin at index ',wvl_idx
       ! stop
    end do                    ! end loop over wvl
    
    return 
  end subroutine wvl_grd_mk                       ! end wvl_grd_mk()
  
  subroutine mlk_abs( &
       lev_nbr,             & ! Level dimension size
       levp_nbr,            & ! Interface dimension size
       prs_HITRAN,          & ! [Pa] Reference pressure
       slr_zen_ngl_cos,     & ! [frc] Solar zenith angle cosine
       !     
       A_phi,               & ! Linear temperature dependence of line strengths
       A_psi,               & ! Linear temperature dependence of Lorentzian HWHM
       B_phi,               & ! Quadratic temperature dependence of line strengths
       B_psi,               & ! Quadratic temperature dependence of Lorentzian HWHM
       !
       S_d_abs_cff_mss,     & ! [m2 kg-1] Band average line strength
       S_p_abs_cff_mss,     & ! [m2 kg-1] Band average line strength amount over line width
       !
       grv,                 & ! [m s-2] Gravity
       mmmr,                & ! [frc] Moist mass mixing ratio of gas
       mpl,                 & ! [kg m-2] Mass path layer of gas
       prs,                 & ! [Pa] Pressure at layer midpoint
       prs_ntf,             & ! [Pa] Pressure at layer interface
       tpt_dlt,             & ! [K] Layer temperature minus Malkmus fit temperature
       tpt_dlt_sqr,         & ! [K2] Square of tpt_dlt
       odal                 & ! [frc] Optical depth absorption layer of gas
       )
    ! Purpose:
    ! Evaluate the absorption optical depth in each layer based on the input 
    ! Malkmus band parameters and the thermodynamic profile
    
    ! Evaluate Malkmus narrow band absorptance for H2O, OH, O2, CO2, CH4, N2O, etc.
    ! The same physics applies to all gases, but there are some subtle differences:
    ! First, Malkmus parameters for all gases except H2O are computed at 5 cm-1 resolution
    ! H2O parameters are computed at 10 cm-1 resolution 
    ! KiR83 show 5 cm-1 bands are optimal for CO2
    ! Bri922 p. 11477 uses 5 cm-1 bands for CO2, O3, CH4, N2O; 10 cm-1 bands for H2O
    ! Kie97 p. 113 recommends 5 cm-1 for CO2, 10 cm-1 for H2O, 5--10 cm-1 for O3, 5 cm-1 for others
    
    ! Second, some weakly absorbing gases (e.g., OH) have bandstrengths which cause underflows in single precision
    ! Thus the htrn2nb code is patched to warn about this occurance
    ! Trace gases (e.g., OH) may also have a zero concentration in some or all levels
    ! This causes singularities in the computation of absorption optical depth in those levels
    ! The fix is to set the pressure weighted mass path to a non-zero value when the absorber path is zero
    ! This fix is currently implemented so that all gases can have zero concentration
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    implicit none
    ! Parameters
    integer lev_nbr_max
    integer levp_nbr_max
    parameter(lev_nbr_max=110)
    parameter(levp_nbr_max=lev_nbr_max+1)
    ! Input
    real(selected_real_kind(p=12))::slr_zen_ngl_cos ! [frc] Solar zenith angle cosine
    integer,intent(in)::lev_nbr           ! Level dimension size
    integer,intent(in)::levp_nbr          ! Interface dimension size
    real,intent(in)::A_phi                ! Linear temperature dependence of line strength
    real,intent(in)::A_psi                ! Linear temperature dependence of Lorentzian HWHM
    real,intent(in)::B_phi                ! Quadratic temperature dependence of line strengths
    real,intent(in)::B_psi                ! Quadratic temperature dependence of Lorentzian HWHM
    real,intent(in)::S_d_abs_cff_mss      ! [m2 kg-1] Band average line strength
    real,intent(in)::S_p_abs_cff_mss      ! [m2 kg-1] Band average line strength amount over line width
    real,intent(in)::grv(lev_nbr)         ! [m s-2] Gravity
    real,intent(in)::mmmr(lev_nbr)        ! [frc] Moist mass mixing ratio of gas
    real,intent(in)::mpl(lev_nbr)         ! [kg m-2] Mass path layer of gas
    real,intent(in)::prs(lev_nbr)         ! [Pa] Pressure at layer midpoint
    real,intent(in)::prs_HITRAN           ! [Pa] Reference pressure
    real,intent(in)::prs_ntf(levp_nbr)    ! [Pa] Pressure at layer interface
    real,intent(in)::tpt_dlt(lev_nbr)     ! [K] Layer temperature minus Malkmus fit temperature
    real,intent(in)::tpt_dlt_sqr(lev_nbr) ! [K2] Square of tpt_dlt
    ! Input/Output
    ! Output
    real,intent(out)::odal(lev_nbr)        ! [frc] Optical depth absorption layer of gas
    ! Local
    integer lev_idx           ! Counting index
    real float_foo            ! Factor in optical depth calculation
    real opt_dep_ITOD(levp_nbr_max) ! Interface transmission optical depth
    real phi_wgt(lev_nbr_max) ! Line strength temperature dependence factor
    real prs_bar(levp_nbr_max) ! Pressure weighted mass path
    real psi_wgt(lev_nbr_max) ! Lorentzian HWHM temperature dependence factor
    real u_bar(levp_nbr_max)  ! Pressure weighted mass path
    ! Main code
    if (dbg_lvl>=dbg_sbr) write (6,'(a)') 'Entering mlk_abs()'
    ! Enough memory? 
    if (lev_nbr>lev_nbr_max) stop 'lev_nbr>lev_nbr_max in mlk_abs()'
    do lev_idx=1,lev_nbr
       ! Phi and psi weights are mid-layer quantities
       phi_wgt(lev_idx)=exp(A_phi*tpt_dlt(lev_idx)+B_phi*tpt_dlt_sqr(lev_idx))
       psi_wgt(lev_idx)=exp(A_psi*tpt_dlt(lev_idx)+B_psi*tpt_dlt_sqr(lev_idx))
    enddo                     ! end loop over lev
    ! u_bar and prs_bar are layer-interface quantities
    ! Compute top interface
    u_bar(1)=mmmr(1)*phi_wgt(1)*prs_ntf(1)/grv(1)
    prs_bar(1)=mmmr(1)*psi_wgt(1)*prs_ntf(1)*prs_ntf(1)/grv(1)
    ! Sum integrands
    ! Recall that dp/g = mpl_mst_air ([kg m-2]) so that, e.g., q*dp/g = mpl [kg m-2]
    ! Using stored values saves lots of floating point operations
    do lev_idx=2,levp_nbr
       u_bar(lev_idx)=u_bar(lev_idx-1)+mpl(lev_idx-1)*phi_wgt(lev_idx-1)
       prs_bar(lev_idx)=prs_bar(lev_idx-1)+mpl(lev_idx-1)*psi_wgt(lev_idx-1)*prs(lev_idx-1)
    enddo                     ! end loop over lev
    ! Normalize pressure weighted mass path
    do lev_idx=1,levp_nbr
       ! Avoid overflows when concentration is zero (currently allowed only for OH)
       if (u_bar(lev_idx)>0.0) then
          prs_bar(lev_idx)=prs_bar(lev_idx)/(u_bar(lev_idx)*prs_HITRAN)
       else
          ! Avoid overflow in optical depth later
          prs_bar(lev_idx)=1.0
       endif                  ! endif
    enddo                     ! end loop over lev
    ! Compute absorption optical depth 
    ! "interface transmission optical depth" (ITOD) is the argument in the exponential that determines transmission 
    ! ITOD has a one over slr_zen_ngl_cos dependence built into it by the Malkmus formula, so ITOD represents full slant path optical depth, NOT zenith optical depth
    ! Perhaps ITOD should be given a diffusivity factor correction for high scattering optical depths because the diffuse radiation does not follow the full slant path
    ! ITOD is an interface quantity
    ! We could also have computed the transmission of a layer as a mid-layer quantity directly, but BPB advises against this for subtle reasons
    do lev_idx=1,levp_nbr
       float_foo=sqrt(1.0+4.0*S_p_abs_cff_mss*u_bar(lev_idx)/ &
            (prs_bar(lev_idx)*slr_zen_ngl_cos))
       opt_dep_ITOD(lev_idx)=(float_foo-1.0)* &
            0.5*S_d_abs_cff_mss*prs_bar(lev_idx)/ &
            S_p_abs_cff_mss
    enddo                     ! end loop over lev
    ! Compute layer "monochromatic" absorption optical depth 
    ! directly from the difference of the column "monochromatic"
    ! optical depths rather than from the log of the difference of
    ! the column transmissions in order to avoid numerical 
    ! difficulties with saturated lines when the transmission is 0.0.
    ! This log transmission method fails in single precision
    do lev_idx=1,lev_nbr
       odal(lev_idx)=odal(lev_idx)+ &
            slr_zen_ngl_cos*(opt_dep_ITOD(lev_idx+1)-opt_dep_ITOD(lev_idx))
    enddo                     ! end loop over lev
    if (dbg_lvl>=dbg_sbr) write (6,'(a)') 'Exiting mlk_abs()'
    return 
  end subroutine mlk_abs                       ! end mlk_abs()
  
end module rt_mdl ! [mdl] Radiative transfer utilities
