! $Id$

! Purpose: Utility routines for column (CLM) processing

! Copyright (C) 1994--2017 Charlie Zender
! This software is distributed under the terms of the GNU General Public License
! See http://www.gnu.org/copyleft/gpl.html for full license text

! These routines are heavily used by clm

! Usage:
! use clm_mdl ! [mdl] Column (CLM) processing

module clm_mdl ! [mdl] Column (CLM) processing
  implicit none
  public ! [stt] Symbols are public unless individually qualified as private
  
contains 
  
  subroutine rfm_read( &
       fl_in_unit, & ! I []
       levp_nbr, & ! I []
       rfm_val & ! O []
       )
    implicit none
    ! Parameters
    integer,parameter::rfm_clm_nbr=5 ! [nbr] Number of values/columns per row
    integer,parameter::levp_nbr_max=150 ! [nbr] Number of values/columns per row
    ! Commons
    ! Input Arguments
    integer,intent(in)::fl_in_unit,levp_nbr
    !    real,dimension(:),allocatable::rfm_raw
    !real,dimension(:),allocatable::rfm_raw
    !real,dimension(levp_nbr_max)::rfm_raw
    real::rfm_raw(levp_nbr_max)
    real,dimension(:),intent(out)::rfm_val
    ! Input/Output Arguments
    ! Output Arguments
    ! Local workspace
    character(len=128)::rfm_ln ! [sng] Buffer into which each RFM line is read
    integer rfm_row_nbr
    integer rfm_rph_nbr
    integer rfm_row_idx
    integer rfm_idx
    integer int_foo           ! [nbr] Integer
    integer lev_nbr
    integer levp_idx
    integer rcd               ! [rcd] Return success code
    ! Main code
    rcd=0

    !allocate(rfm_raw(levp_nbr),stat=rcd)

    lev_nbr=levp_nbr-1 ! [nbr] dimension size
    rfm_row_nbr=levp_nbr/rfm_clm_nbr
    rfm_rph_nbr=mod(levp_nbr,rfm_clm_nbr)
    if (rfm_rph_nbr /= 0) rfm_row_nbr=rfm_row_nbr+1
    write (6,*) 'rfm_row_nbr = ',rfm_row_nbr,', rfm_clm_nbr = ',rfm_clm_nbr,', rfm_rph_nbr = ',rfm_rph_nbr
    
    do rfm_row_idx=1,rfm_row_nbr
       read (fl_in_unit,'(a)') rfm_ln
       if (rfm_row_idx < rfm_row_nbr) then
          read (rfm_ln,*) (rfm_raw(rfm_idx),rfm_idx=(rfm_row_idx-1)*rfm_clm_nbr,rfm_row_idx*rfm_clm_nbr-1)
       else
          read (rfm_ln,*) (rfm_raw(rfm_idx),rfm_idx=(rfm_row_idx-1)*rfm_clm_nbr,(rfm_row_idx-1)*rfm_clm_nbr+rfm_rph_nbr)
       endif
    enddo ! rfm_row_idx
    
    do levp_idx=1,levp_nbr
       write (6,*) 'idx = ',levp_idx,', rfm_raw = ',rfm_raw(levp_idx)
    enddo

    do levp_idx=1,levp_nbr
       int_foo=levp_nbr-levp_idx+1
       rfm_val(int_foo)=rfm_raw(levp_idx)
    enddo

    if (rcd /= 0 ) write (6,*) 'ERROR in rfm_read()'
    !    if (allocated(rfm_raw)) deallocate(rfm_raw,stat=rcd)

    return
  end subroutine rfm_read ! end rfm_read()
  
  subroutine slr_crd_Bri92( &
       lat,                 & ! I [rdn] Latitude
       lcl_yr_day,          & ! I [day] Local year day
       slr_zen_ngl_cos,     & ! O [frc] Solar zenith angle cosine
       xnt_fac)             ! O [frc] Eccentricity factor
    ! Purpose: Compute solar geometry
    ! Reference: B. P. Briegleb's routine used in CCM2, CCM3
    
    implicit none
    ! Parameters
    real,parameter::days_per_year=365.0 ! days
    ! Commons
    ! Input Arguments
    real(selected_real_kind(p=12)),intent(in)::lat      ! [rdn] 
    real(selected_real_kind(p=12)),intent(in)::lcl_yr_day ! [day] Local year day
    ! Input/Output Arguments
    ! Output Arguments
    real(selected_real_kind(p=12)),intent(out)::slr_zen_ngl_cos ! [frc] Solar zenith angle cosine
    real,intent(out)::xnt_fac              ! [frc]
    ! Local workspace
    real(selected_real_kind(p=12))::&
         delta,               & ! solar declination (radian)
         phi,                 & ! local phase angle (0 is midnight) (radian)
         pi, & ! [frc] 3
         theta                ! solar polar coordinate angle (radian)
    ! Main code
    pi=4.0*atan(1.0)
    
    ! Compute eccentricity factor (Sun-Earth distance factor)
    theta=2.0*pi*lcl_yr_day/days_per_year
    xnt_fac= &
         1.000110+0.034221*cos(theta)+0.001280*sin(theta)+ &
         0.000719*cos(2.0*theta)+0.000077*sin(2.0*theta)
    
    ! Solar declination in radians
    delta = 0.006918-0.399912*cos(theta)+0.070257*sin(theta)- &
         0.006758*cos(2.0*theta)+0.000907*sin(2.0*theta)- &
         0.002697*cos(3.0*theta)+0.001480*sin(3.0*theta)
    
    phi=2.0*pi*lcl_yr_day
    slr_zen_ngl_cos= &
         sin(lat)*sin(delta)- &
         cos(lat)*cos(delta)*cos(phi)
    
    return
  end subroutine slr_crd_Bri92                       ! end slr_crd_Bri92()
  
  real function q_O3_ntp(prs_ntp)
    ! Purpose:
    ! Compute mixing ratio of O3 based by interpolating given
    ! pressure into ICRCCM O3 profile. 
    ! I added an extra level beneath the ICRCCM surface pressure 
    ! to assure that all reasonable surface pressures are always bracketed.
    ! Input pressure in Pascal 
    ! O3 mixing ratio returned in kg kg-1
    use xtr_mdl ! [mdl] Extrapolation/interpolation handling
    use vec_mdl ! [mdl] Vector manipulation, interpolation, rebinning
    implicit none
    ! Parameters
    integer,parameter::lev_nbr=94
    ! Input Arguments
    real,intent(in)::prs_ntp              ! Pa
    ! Input/Output Arguments
    ! Output Arguments
    ! Local workspace
    real prs(lev_nbr)         ! Pa
    real q_O3(lev_nbr)        ! kg kg-1
    data prs/ &
         0.,     009.29,     010.83,     012.63,     014.72, &
         017.16,     020.01,     023.33,     027.20,     031.72, &
         036.98,     043.11,     050.27,     058.61,     068.33, &
         079.67,     092.89,     108.30,     126.26,     147.21, &
         171.64,     200.11,     233.32,     272.03,     317.16, &
         369.78,     431.13,     502.66,     586.06,     683.30, &
         796.67,     928.85,     1083.00,     1262.60,     1472.10, &
         1716.40,     2001.10,     2333.20,     2720.30,     3171.60, &
         3697.80,     4311.30,     5026.60,     5860.60,     6833.00, &
         7966.70,     9288.50,     11000.00,     13000.00,     15000.00, &
         17000.00,     19000.00,     21000.00,     23000.00,     25000.00, &
         27000.00,     29000.00,     31000.00,     33000.00,     35000.00, &
         37000.00,     39000.00,     41000.00,     43000.00,     45000.00, &
         47000.00,     49000.00,     51000.00,     53000.00,     55000.00, &
         57000.00,     59000.00,     61000.00,     63000.00,     65000.00, &
         67000.00,     69000.00,     71000.00,     73000.00,     75000.00, &
         77000.00,     79000.00,     81000.00,     83000.00,     85000.00, &
         87000.00,     89000.00,     91000.00,     93000.00,     95000.00, &
         97000.00,     99000.00,     100650.00,     110000.00/ 
    
    data q_O3/ &
         0.,     1.909E-06,     2.158E-06,     2.381E-06,     2.581E-06, &
         2.761E-06,     2.925E-06,     3.074E-06,     3.210E-06,     3.335E-06, &
         3.451E-06,     3.559E-06,     3.659E-06,     3.753E-06,     3.836E-06, &
         3.890E-06,     4.306E-06,     4.973E-06,     5.534E-06,     6.004E-06, &
         6.785E-06,     7.809E-06,     8.656E-06,     9.356E-06,     9.875E-06, &
         1.016E-05,     1.038E-05,     1.055E-05,     1.072E-05,     1.132E-05, &
         1.194E-05,     1.244E-05,     1.285E-05,     1.276E-05,     1.158E-05, &
         1.052E-05,     9.614E-06,     8.838E-06,     7.913E-06,     6.592E-06, &
         5.156E-06,     3.822E-06,     2.735E-06,     1.829E-06,     1.145E-06, &
         6.656E-07,     4.283E-07,     2.632E-07,     2.106E-07,     1.824E-07, &
         1.629E-07,     1.473E-07,     1.320E-07,     1.193E-07,     1.085E-07, &
         9.940E-08,     9.201E-08,     8.711E-08,     8.284E-08,     7.901E-08, &
         7.557E-08,     7.325E-08,     7.169E-08,     7.026E-08,     6.886E-08, &
         6.756E-08,     6.636E-08,     6.517E-08,     6.406E-08,     6.302E-08, &
         6.201E-08,     6.103E-08,     6.012E-08,     5.930E-08,     5.891E-08, &
         5.865E-08,     5.840E-08,     5.814E-08,     5.763E-08,     5.703E-08, &
         5.646E-08,     5.592E-08,     5.535E-08,     5.463E-08,     5.396E-08, &
         5.331E-08,     5.270E-08,     5.200E-08,     5.108E-08,     5.020E-08, &
         4.934E-08,     4.852E-08,     4.786E-08,     4.000E-08/
    ! Main code
    ! Initialize default values
    ! Interpolate q_O3 to prs_ntp
    q_O3_ntp=ntp_vec_one(lev_nbr,prs,q_O3,prs_ntp)
    ! call ntp_vec(lev_nbr,prs,q_O3,1,prs_ntp,q_O3_out,xtr_typ_LHS,xtr_typ_RHS)
    ! q_O3_ntp=q_O3_out
    return
  end function q_O3_ntp                       ! end q_O3_ntp()
  
  real function q_NO2_ntp(prs)
    ! Purpose:
    ! Compute mixing ratio of NO2 based on constant relative concentration.
    ! Input volume mixing ratio as a fraction, i.e., 350 ppmv = 350/1.0e6 = 3.5e-4
    ! Input pressure in pascal. 
    ! NO2 mixing ratio returned in kg/kg.
    use phys_cst_mdl,only:mmw_dry_air,mmw_NO2 ! [mdl] Fundamental and derived physical constants
    implicit none
    ! Parameters
    ! Clean troposphere NO2_vmr is 0.1--0.5 ppbv, average is 0.3 ppbv (Sei86 p. 37)
    ! Polluted troposphere NO2_vmr is 50--250 ppbv, average is 150 ppbv (Sei86 p. 37)
    real,parameter::NO2_vmr=0.30e-09 
    !  real,parameter::NO2_vmr=150.0e-09
    real,parameter::SO2_vmr=5.0e-09
    real,parameter::CO_vmr=120.0e-09
    real,parameter::NO_vmr=0.03e-09
    real,parameter::HNO3_vmr=0.15e-09
    real,parameter::NH3_vmr=1.0e-09
    ! Input Arguments
    real,intent(in)::prs                  ! [Pa]
    ! Input/Output Arguments
    ! Output Arguments
    ! Local workspace
    ! Main code
    
    if (prs > 20000.0) then
       ! Set NO2 mixing ratio to clean tropospheric value everywhere beneath 200 mb
       q_NO2_ntp=NO2_vmr*(mmw_NO2/mmw_dry_air)
    else
       ! BrS84 p. 443 give mid-latitude equinox values above 10 km (269 mb)
       ! from equilibrium in a 1-D radiative convective photochemical model
       q_NO2_ntp=0.0
    endif
    
    return
  end function q_NO2_ntp                       ! end q_NO2_ntp()
  
  real function q_OH_ntp(prs)
    ! Purpose:
    ! Compute mixing ratio of OH based on constant relative concentration.
    ! Input volume mixing ratio as a fraction, i.e., 350 ppmv = 350/1.0e6 = 3.5e-4
    ! Input pressure in pascal 
    ! OH mixing ratio returned in kg/kg
    implicit none
    ! Parameters
    ! Clean troposphere OH_mmr is XXX ppbv, average is XXX ppbv
    ! Polluted troposphere OH_mmr is XXX ppbv, average is XXX ppbv
    real,parameter::OH_mmr=1.0e-13
    ! Input Arguments
    real,intent(in)::prs                  ! pascal
    ! Input/Output Arguments
    ! Output Arguments
    ! Local workspace
    ! Main code
    
    if (prs > 20000.0) then
       ! BrS84 p. 443 give mid-latitude equinox values
       ! from equilibrium in a 1-D radiative convective photochemical model
       q_OH_ntp=OH_mmr
    else
       ! Set OH mixing ratio to clean tropospheric value everywhere beneath 200 mb
       q_OH_ntp=0.0
    endif
    
    return
  end function q_OH_ntp                       ! end q_OH_ntp()
  
  subroutine aer_info_get(fl_aer,wvl_obs_aer,dns_aer,ext_cff_mss_aer_spc)
    ! Purpose: Retrieve microphysical optical properties from aerosol file
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use netcdf ! [mdl] netCDF interface
    use nf90_utl,only:nf90_wrp,nf90_wrp_close,nf90_wrp_open,nf90_wrp_inq_dimid,nf90_wrp_inq_varid ! [mdl] netCDF utilities
    use sng_mdl,only:ftn_strlen   ! [mdl] String manipulation
    use vec_mdl ! [mdl] Vector manipulation, interpolation, rebinning
    use xtr_mdl ! [mdl] Extrapolation/interpolation handling
    implicit none
    ! Parameters
    character(len=*),parameter::sbr_nm='aer_info_get' ! [sng] Subroutine name
    integer,parameter::aer_bnd_nbr_max=200 ! 0.05 micron resolution from 0.3--5.0 microns
    ! Input Arguments
    character(80),intent(in)::fl_aer
    real,intent(in)::wvl_obs_aer
    ! Input/Output Arguments
    ! Output Arguments
    real,intent(out)::dns_aer
    real,intent(out)::ext_cff_mss_aer_spc
    ! Local workspace
    integer aer_bnd_nbr       ! dimension size
    integer rcd
    integer nc_id
    integer bnd_dim_id
    ! Aerosol input variables
    integer dns_aer_id
    integer asm_prm_aer_id
    integer ext_cff_mss_aer_id
    integer abs_cff_mss_aer_id
    integer sca_cff_mss_aer_id
    integer wvl_ctr_aer_id
    integer wvl_max_aer_id
    integer wvl_min_aer_id
    ! Allocatable variables
    ! Aerosol input variables
    real,dimension(:),allocatable::abs_cff_mss_aer
    real,dimension(:),allocatable::asm_prm_aer
    real,dimension(:),allocatable::ext_cff_mss_aer
    real,dimension(:),allocatable::sca_cff_mss_aer
    real,dimension(:),allocatable::wvl_ctr_aer
    real,dimension(:),allocatable::wvl_max_aer
    real,dimension(:),allocatable::wvl_min_aer
    ! Main code
    
    ! Initialize default values
    rcd=nf90_noerr              ! nf90_noerr == 0
    
    ! Get aerosol data
    rcd=nf90_wrp_open(fl_aer,nf90_nowrite,nc_id)
    ! Get dimension IDs
    rcd=nf90_wrp_inq_dimid(nc_id,'wvl',bnd_dim_id)
    ! Get dimension sizes
    rcd=nf90_wrp(nf90_inquire_dimension(nc_id,bnd_dim_id,len=aer_bnd_nbr),sbr_nm//": inquire_dim bnd")
    if (aer_bnd_nbr > aer_bnd_nbr_max) stop 'aer_bnd_nbr > aer_bnd_nbr_max'
    ! Allocate space for dynamic arrays
    allocate(abs_cff_mss_aer(aer_bnd_nbr),stat=rcd)
    allocate(asm_prm_aer(aer_bnd_nbr),stat=rcd)
    allocate(ext_cff_mss_aer(aer_bnd_nbr),stat=rcd)
    allocate(sca_cff_mss_aer(aer_bnd_nbr),stat=rcd)
    allocate(wvl_ctr_aer(aer_bnd_nbr),stat=rcd)
    allocate(wvl_max_aer(aer_bnd_nbr),stat=rcd)
    allocate(wvl_min_aer(aer_bnd_nbr),stat=rcd)
    ! Get variable IDs
    rcd=nf90_wrp_inq_varid(nc_id,'ext_cff_mss',ext_cff_mss_aer_id)
    rcd=nf90_wrp_inq_varid(nc_id,'abs_cff_mss',abs_cff_mss_aer_id)
    rcd=nf90_wrp_inq_varid(nc_id,'sca_cff_mss',sca_cff_mss_aer_id)
    rcd=nf90_wrp_inq_varid(nc_id,'dns_aer',dns_aer_id)
    rcd=nf90_wrp_inq_varid(nc_id,'asm_prm',asm_prm_aer_id)
    rcd=nf90_wrp_inq_varid(nc_id,'wvl_max',wvl_max_aer_id)
    rcd=nf90_wrp_inq_varid(nc_id,'wvl_min',wvl_min_aer_id)
    rcd=nf90_wrp_inq_varid(nc_id,'wvl',wvl_ctr_aer_id)
    ! Get data
    rcd=nf90_wrp(nf90_get_var(nc_id,ext_cff_mss_aer_id,ext_cff_mss_aer),sbr_nm//": get_var ext_cff_mss_aer")
    rcd=nf90_wrp(nf90_get_var(nc_id,abs_cff_mss_aer_id,abs_cff_mss_aer),sbr_nm//": get_var abs_cff_mss_aer")
    rcd=nf90_wrp(nf90_get_var(nc_id,sca_cff_mss_aer_id,sca_cff_mss_aer),sbr_nm//": get_var sca_cff_mss_aer")
    rcd=nf90_wrp(nf90_get_var(nc_id,dns_aer_id,dns_aer),sbr_nm//": get_var dns_aer")
    rcd=nf90_wrp(nf90_get_var(nc_id,asm_prm_aer_id,asm_prm_aer),sbr_nm//": get_var asm_prm_aer")
    rcd=nf90_wrp(nf90_get_var(nc_id,wvl_max_aer_id,wvl_max_aer),sbr_nm//": get_var wvl_max_aer")
    rcd=nf90_wrp(nf90_get_var(nc_id,wvl_min_aer_id,wvl_min_aer),sbr_nm//": get_var wvl_min_aer")
    rcd=nf90_wrp(nf90_get_var(nc_id,wvl_ctr_aer_id,wvl_ctr_aer),sbr_nm//": get_var wvl_ctr_aer")
    ! Close file
    rcd=nf90_wrp_close(nc_id,fl_aer,'Ingested') ! [fnc] Close file
    
    ! Find visible extinction efficiency of aerosol
    ! call ntp_vec(aer_bnd_nbr,wvl_ctr_aer,ext_cff_mss_aer,1,wvl_obs_aer,ext_cff_mss_aer_spc,xtr_typ_LHS,xtr_typ_RHS)
    ext_cff_mss_aer_spc=ntp_vec_one(aer_bnd_nbr,wvl_ctr_aer,ext_cff_mss_aer,wvl_obs_aer)
    
    ! De-allocate dynamic variables
    if (allocated(abs_cff_mss_aer)) deallocate(abs_cff_mss_aer,stat=rcd)
    if (allocated(asm_prm_aer)) deallocate(asm_prm_aer,stat=rcd)
    if (allocated(ext_cff_mss_aer)) deallocate(ext_cff_mss_aer,stat=rcd)
    if (allocated(sca_cff_mss_aer)) deallocate(sca_cff_mss_aer,stat=rcd)
    if (allocated(wvl_ctr_aer)) deallocate(wvl_ctr_aer,stat=rcd)
    if (allocated(wvl_max_aer)) deallocate(wvl_max_aer,stat=rcd)
    if (allocated(wvl_min_aer)) deallocate(wvl_min_aer,stat=rcd)
    
    return
  end subroutine aer_info_get                       ! end aer_info_get()
  
  subroutine aer_odxc_get(wvl_obs_aer,time_mdl,odxc_obs_aer)
    ! Purpose: Retrieve measured optical depth and interpolate to given time
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use netcdf ! [mdl] netCDF interface
    use nf90_utl,only:nf90_wrp,nf90_wrp_close,nf90_wrp_open,nf90_wrp_inq_dimid,nf90_wrp_inq_varid ! [mdl] netCDF utilities
    use sng_mdl,only:ftn_strlen,ftn_strnul ! [mdl] String manipulation
    implicit none
    ! Parameters
    character(len=*),parameter::sbr_nm='aer_odxc_get' ! [sng] Subroutine name
    integer,parameter::chn_nbr_max=5 ! 5 MFRSR channels
    integer,parameter::time_nbr_max=312 ! Michalsky's ARESE obs have 312 records
    ! Input Arguments
    real(selected_real_kind(p=12)),intent(in)::time_mdl
    real,intent(in)::wvl_obs_aer
    ! Input/Output Arguments
    ! Output Arguments
    real,intent(out)::odxc_obs_aer
    ! Local workspace
    character(80) fl_aer
    integer chn_nbr           ! dimension size
    integer time_nbr          ! dimension size
    integer brk_lft_idx_chn
    integer brk_rgt_idx_chn
    integer brk_lft_idx_time
    integer brk_rgt_idx_time
    integer rcd
    integer chn_dim_id
    integer time_dim_id
    integer nc_id
    real odxc_lft_ntp   
    real odxc_rgt_ntp   
    ! Aerosol input variables
    integer chn_id
    integer odxc_aer_MFRSR_id
    integer time_id
    ! Allocatable variables
    ! Aerosol input variables
    real(selected_real_kind(p=12)),dimension(:),allocatable::time_obs
    real,dimension(:),allocatable::chn
    real,dimension(:,:),allocatable::odxc_aer_MFRSR
    ! Main code
    
    ! Initialize default values
    rcd=nf90_noerr              ! nf90_noerr == 0
    fl_aer='/data/zender/arese/mfrsr/arese_mfrsr.nc'//char(0)
    call ftn_strnul(fl_aer)
    
    ! Get aerosol data
    rcd=nf90_wrp_open(fl_aer,nf90_nowrite,nc_id)
    ! Get dimension IDs
    rcd=nf90_wrp_inq_dimid(nc_id,'chn',chn_dim_id)
    rcd=nf90_wrp_inq_dimid(nc_id,'time',time_dim_id)
    
    ! Get dimension sizes
    rcd=nf90_wrp(nf90_inquire_dimension(nc_id,chn_dim_id,len=chn_nbr),sbr_nm//": inquire_dim chn")
    if (chn_nbr > chn_nbr_max) stop 'chn_nbr > chn_nbr_max'
    rcd=nf90_wrp(nf90_inquire_dimension(nc_id,time_dim_id,len=time_nbr),sbr_nm//": inquire_dim time")
    if (time_nbr > time_nbr_max) stop 'time_nbr > time_nbr_max'
    
    ! Allocate space for dynamic arrays
    allocate(chn(chn_nbr),stat=rcd)
    allocate(odxc_aer_MFRSR(chn_nbr,time_nbr),stat=rcd)
    allocate(time_obs(time_nbr),stat=rcd)
    
    ! Get variable IDs
    rcd=nf90_wrp_inq_varid(nc_id,'chn',chn_id)
    rcd=nf90_wrp_inq_varid(nc_id,'odxc_aer_MFRSR',odxc_aer_MFRSR_id)
    rcd=nf90_wrp_inq_varid(nc_id,'time',time_id)
    
    ! Get data
    rcd=nf90_wrp(nf90_get_var(nc_id,chn_id,chn),sbr_nm//": get_var chn")
    rcd=nf90_wrp(nf90_get_var(nc_id,odxc_aer_MFRSR_id,odxc_aer_MFRSR),sbr_nm//": get_var odxc_aer_MFRSR")
    rcd=nf90_wrp(nf90_get_var(nc_id,time_id,time_obs),sbr_nm//": get_var time_obs")
    
    ! Close file
    rcd=nf90_wrp_close(nc_id,fl_aer,'Ingested') ! [fnc] Close file
    ! Find bracketing channels
    brk_lft_idx_chn=1
    do while((chn(brk_lft_idx_chn) < wvl_obs_aer).and.(brk_lft_idx_chn <= chn_nbr)) 
       brk_lft_idx_chn=brk_lft_idx_chn+1
    end do                    ! end loop over chn
    if (brk_lft_idx_chn > chn_nbr) then
       stop 'brk_lft_idx_chn > chn_nbr' 
    else 
       brk_lft_idx_chn=brk_lft_idx_chn-1
    endif
    if (brk_lft_idx_chn < chn_nbr) then
       brk_rgt_idx_chn=brk_lft_idx_chn+1
    else
       brk_rgt_idx_chn=brk_lft_idx_chn 
    endif
    
    ! Find bracketing times
    brk_lft_idx_time=1
    do while((time_obs(brk_lft_idx_time) < time_mdl).and.(brk_lft_idx_time <= time_nbr)) 
       brk_lft_idx_time=brk_lft_idx_time+1
    end do                    ! end loop over time
    if (brk_lft_idx_time > time_nbr) then
       stop 'brk_lft_idx_time > time_nbr' 
    else 
       brk_lft_idx_time=brk_lft_idx_time-1
    endif
    if (brk_lft_idx_time < time_nbr) then
       brk_rgt_idx_time=brk_lft_idx_time+1
    else
       brk_rgt_idx_time=brk_lft_idx_time 
    endif
    
    ! Find odxc at bracketing channels at exact time (i.e., adjust for partial derivative in time)
    if (brk_lft_idx_time /= brk_rgt_idx_time) then
       odxc_lft_ntp= &
            odxc_aer_MFRSR(brk_lft_idx_chn,brk_lft_idx_time)+ &
            (time_mdl-time_obs(brk_lft_idx_time))* &
            (odxc_aer_MFRSR(brk_lft_idx_chn,brk_rgt_idx_time)-odxc_aer_MFRSR(brk_lft_idx_chn,brk_lft_idx_time))/ &
            (time_obs(brk_rgt_idx_time)-time_obs(brk_lft_idx_time))
       odxc_rgt_ntp= &
            odxc_aer_MFRSR(brk_rgt_idx_chn,brk_lft_idx_time)+ &
            (time_mdl-time_obs(brk_lft_idx_time))* &
            (odxc_aer_MFRSR(brk_rgt_idx_chn,brk_rgt_idx_time)-odxc_aer_MFRSR(brk_rgt_idx_chn,brk_lft_idx_time))/ &
            (time_obs(brk_rgt_idx_time)-time_obs(brk_lft_idx_time))
    else
       odxc_lft_ntp=odxc_aer_MFRSR(brk_lft_idx_chn,brk_lft_idx_time)
       odxc_rgt_ntp=odxc_aer_MFRSR(brk_rgt_idx_chn,brk_rgt_idx_time)
    endif
    
    ! Linearly interpolate to find exact value at exact time (i.e., adjust for partial derivative in wavelength)
    if (brk_lft_idx_chn /= brk_rgt_idx_chn) then
       odxc_obs_aer= &
            odxc_lft_ntp+ &
            (wvl_obs_aer-chn(brk_lft_idx_chn))* &
            (odxc_rgt_ntp-odxc_lft_ntp)/ &
            (chn(brk_rgt_idx_chn)-chn(brk_lft_idx_chn))
    else
       odxc_obs_aer=odxc_lft_ntp
    endif
    
    if (dbg_lvl > 2) then
       write (0,*) 'time_mdl  = ',time_mdl
       write (0,*) 'time_obs(brk_lft_idx_time)  = ',time_obs(brk_lft_idx_time)
       write (0,*) 'time_obs(brk_rgt_idx_time)  = ',time_obs(brk_rgt_idx_time)
       write (0,*) 'wvl_obs_aer  = ',wvl_obs_aer
       write (0,*) 'chn(brk_lft_idx_chn)  = ',chn(brk_lft_idx_chn)
       write (0,*) 'chn(brk_rgt_idx_chn)  = ',chn(brk_rgt_idx_chn)
       write (0,*) 'odxc_lft_ntp  = ',odxc_lft_ntp
       write (0,*) 'odxc_rgt_ntp  = ',odxc_rgt_ntp
       write (0,*) 'odxc_obs_aer  = ',odxc_obs_aer
    endif                     ! endif dbg
    
    ! De-allocate dynamic variables
    if (allocated(chn)) deallocate(chn,stat=rcd)
    if (allocated(odxc_aer_MFRSR)) deallocate(odxc_aer_MFRSR,stat=rcd)
    if (allocated(time_obs)) deallocate(time_obs,stat=rcd)
    
    ! Sanity check
    if (odxc_obs_aer < 0.0) stop 'odxc_obs_aer < 0.0'
    return
  end subroutine aer_odxc_get                       ! end aer_odxc_get()
  
end module clm_mdl ! [mdl] Column (CLM) processing
