! $Id$

program htrn2nb
  
  ! Purpose: Compute Malkmus narrow band parameters from HITRAN line data
  ! If input is from BPB, just convert Malkmus band parameters to netCDF format
  
  ! 20181005: Update Fortran style
  ! 19981128: Default to run htrn2nb in double precision
  ! 19981128: Noticed that new narrow band parameters are ~10% larger than old
  ! Compared S_p_abs_cff_mss in mlk_H2O.nc vs. swnb_H2O.nc
  ! Could this be since the addition of rotation partition function exponents?
  ! fxm: Should try to track this down, verify parameters are correct
  
  ! Compilation:
  ! cd ${HOME}/aca; make -W htrn2nb.F90 OPTS=D precision=double htrn2nb; cd -
  ! cd ${HOME}/aca; make -W htrn2nb.F90 precision=double htrn2nb; cd -
  ! cd ${HOME}/aca; make OPTS=D precision=double htrn2nb; cd -
  
  ! KiR83 show 5 cm-1 bands are optimal for CO2
  ! Bri92 p. 11477 uses 5 cm-1 bands for CO2, O3, CH4, N2O; 10 cm-1 bands for H2O
  ! Kie97 p. 113 recommends 5 cm-1 for CO2, 10 cm-1 for H2O, 5--10 cm-1 for O3, 5 cm-1 for others
  
  ! NB: Compiling and running program in single precision may cause inexplicable crashes and floating point errors 
  ! This depends on gas species and wavelength interval because some line strengths are simply too small ( < 10^-36) to represent in single precision
  ! Even when line strength can be represented in single precision, the arithmetic below (performed in SI units) may cause IEEE underflow, etc
  ! Thus our strategy is to read and process the data in double precision, and then store results in single precision
  
  ! Usage:
  ! For H2O vapor (default) one can just use htrn2nb
  ! Other gases require explicitly setting # of bands and names and molecule #:
  
  ! Production:
  ! Shell script ${HOME}/aca/htrn2nb.sh processes most useful gases
  
  ! Debugging:
  ! htrn2nb -t 28 -l 1.0e-35 -h 1.0 -b 2 -i ${DATA}/hitran/foo.nc -o ${DATA}/hitran/foo_out.nc
  
  ! BPB-style input:
  ! htrn2nb -b 3180 -i ${DATA}/aca/CO2.dat -o ${DATA}/aca/mlk_CO2.nc
  
  ! Water Vapor:
  ! htrn2nb -l 2000.0 -h 27000.0 -b 2500 -i ${DATA}/hitran/H2O.nc -o ${DATA}/aca/mlk_H2O.nc
  ! htrn2nb -l 0.0 -h 25000.0 -b 2500 -i ${DATA}/hitran/H2O.nc -o ${DATA}/aca/mlk_H2O_all.nc
  ! htrn2nb -l 2000.0 -h 27000.0 -b 2500 -i ${DATA}/hitran/1H2_16O.nc -o ${DATA}/aca/mlk_1H2_16O.nc
  ! htrn2nb -l 2000.0 -h 27000.0 -b 2500 -i ${DATA}/hitran/1H2_17O.nc -o ${DATA}/aca/mlk_1H2_17O.nc
  ! htrn2nb -l 2000.0 -h 27000.0 -b 2500 -i ${DATA}/hitran/1H2_18O.nc -o ${DATA}/aca/mlk_1H2_18O.nc
  ! htrn2nb -l 2000.0 -h 27000.0 -b 2500 -i ${DATA}/hitran/1H_2H_16O.nc -o ${DATA}/aca/mlk_1H_2H_16O.nc
  
  ! Carbon dioxide:
  ! htrn2nb -l 2000.0 -h 27000.0 -b 5000 -i ${DATA}/hitran/CO2.nc -o ${DATA}/aca/mlk_CO2.nc

  ! Carbon monoxide:
  ! htrn2nb -l 2000.0 -h 27000.0 -b 5000 -i ${DATA}/hitran/CO.nc -o ${DATA}/aca/mlk_CO.nc
  
  use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
  use drv_cst_mdl ! [mdl] Derived physical constants
  use flx_slr_mdl ! [mdl] Solar spectral fluxes
  use fnd_cst_mdl ! [mdl] Fundamental physical constants
  use htrn_mdl ! [mdl] HITRAN constants, subroutines
  use netcdf ! [mdl] netCDF interface
  use nf90_utl ! [mdl] netCDF utilities
  use sng_mdl ! [mdl] String manipulation
  use utl_mdl ! [mdl] Utility functions (date_time_get,mnt_chk...)

  implicit none
  ! Parameters
  character(len=*),parameter::CVS_Date='$Date$' ! [sng] Date string
  character(len=*),parameter::CVS_Id='$Id$' ! [sng] CVS Identification
  character(len=*),parameter::CVS_Name='$HeadURL$' ! [sng] File name string
  character(len=*),parameter::CVS_Revision='$Revision$' ! [sng] File revision string
  character(len=*),parameter::nlc=char(0) ! [sng] NUL character = ASCII 0 = char(0)
  character(len=*),parameter::sbr_nm='htrn2nb' ! [sng] Subroutine name
  integer,parameter::r4=selected_real_kind(p=6) ! r4: 4B (C float) default, 8B (C double) possible
  integer,parameter::bnd_nbr_max=50000   ! 0.0--50000.0 cm-1 in 1 cm-1 bands
  integer,parameter::fl_in_unit=73
  integer,parameter::ln_nbr_max=600000 ! 20181001 HITRAN16 O3 from 0.2--100.0 um has 449570, CO2 has 559874
  integer,parameter::tpt_nbr_max=140     ! 180.0--320.0 K in 1 dgr increments
  integer,parameter::iso_sng_htrn_lng_max=20 ! [nbr] Maximum length of HITRAN isopotomer string
  integer,parameter::mlc_sng_htrn_lng_max=10 ! [nbr] Maximum length of HITRAN molecule string
  real(selected_real_kind(p=12)),parameter::ln_lo_min=1.0e-10   ! [cm-1] Minimum wavenumber
  real(selected_real_kind(p=12)),parameter::tpt_Malkmus_BPB=250.0 ! [K] Reference temperature for Malkmus parameters

  ! Locals with simple initialization and no command-line override
  ! integer::exit_status=0 ! [enm] Program exit status (non-standard Fortran)
  integer::rcd=nf90_noerr ! [rcd] Return success code
  
  character(2)::dsh_key ! [sng] command-line dash and switch
  character(200)::cmd_ln ! [sng] command-line
  character(200)::prg_ID ! [sng] Program ID
  character(26)::lcl_date_time ! [sng] Time formatted as Day Mth DD HH:MM:SS TZ YYYY
  character(80)::arg_val ! [sng] command-line argument value
  character(80)::opt_sng ! [sng] Option string

  ! Command-line parsing
  integer::arg_idx ! [idx] Counting index
  integer::arg_nbr ! [nbr] Number of command-line arguments
  integer::opt_lng ! [nbr] Length of option

  ! netCDF4 
  integer::dfl_lvl=0 ! [enm] Deflate level
  integer::flg_shf=1 ! [flg] Turn on netCDF4 shuffle filter
  integer::flg_dfl=1 ! [flg] Turn on netCDF4 deflate filter
  integer::fl_out_fmt=nco_format_undefined ! [enm] Output file format
  integer::nf90_create_mode=nf90_clobber ! [enm] Mode flag for nf90_create() call

  logical nc_flg
  
  integer bnd_dmn_id        ! Dimension ID for band
  integer bnd_dbg           ! [idx] Debugging band
  integer bnd_idx           ! Counting index
  integer bnd_nbr           ! Dimension size
  integer cnt_iso(1)
  integer cnt_ln(1)
  integer cnt_mlc(1)
  integer grd_dmn_id        ! Dimension ID for grid
  integer idx               ! Counting index
  integer iso_dmn_id        ! Dimension ID for iso
  integer iso_nbr           ! Dimension size
  integer ln_ctr_dmn_id     ! Dimension ID for ln
  integer ln_idx            ! Counting index
  integer ln_nbr            ! Dimension size
  integer mlc_dmn_id        ! Dimension ID for mlc
  integer mlc_nbr           ! Dimension size
  integer typ_out           ! Output floating point type
  integer nc_id             ! file handle
  integer srt_one(1)
  integer tpt_dmn_id        ! Dimension ID for tpt
  integer tpt_idx           ! Counting index
  integer tpt_nbr           ! Dimension size
  
  ! HITRAN netCDF input variables
  integer HWHM_air_id
  integer HWHM_tpt_dpn_xpn_id
  integer iso_id_id
  integer ln_ctr_id
  integer ln_iso_id
  integer ln_nrg_lwr_id
  integer ln_str_id
  integer mlc_id_id
  
  integer iso_id(iso_per_mlc_nbr_max_htrn)
  integer mlc_id
  
  integer A_phi_id
  integer A_psi_id
  integer B_phi_id
  integer B_psi_id
  integer S_d_abs_cff_mss_id
  integer S_d_id
  integer S_p_abs_cff_mss_id
  integer S_p_id
  integer alpha_id
  integer bnd_id            ! Coordinate ID
  integer bnd_ln_nbr_id
  integer bnd_dlt_id
  integer frc_slr_flx_LaN68_id
  integer frc_slr_flx_ThD71_id
  integer idx_rfr_air_STP_id
  integer oneD_foo_id
  integer phi_exact_id
  integer phi_fit_id
  integer psi_exact_id
  integer psi_fit_id
  integer scalar_foo_id
  integer tpt_id
  integer wvl_ctr_id
  integer wvl_grd_id
  integer wvl_max_id
  integer wvl_min_id
  integer wvl_dlt_id
  integer wvn_ctr_id
  integer wvn_grd_id
  integer wvn_max_id
  integer wvn_min_id
  integer wvn_dlt_id
  
  ! Computational-precision netCDF output variables
  real(selected_real_kind(p=12)) bnd_dlt
  real(selected_real_kind(p=12)) scalar_foo
  
  ! Allocatable variables
  ! Computational-precision HITRAN netCDF input variables
  real(selected_real_kind(p=12)),dimension(:),allocatable::HWHM_air
  real(selected_real_kind(p=12)),dimension(:),allocatable::HWHM_tpt_dpn_xpn
  real(selected_real_kind(p=12)),dimension(:),allocatable::ln_ctr
  real(selected_real_kind(p=12)),dimension(:),allocatable::ln_nrg_lwr
  real(selected_real_kind(p=12)),dimension(:),allocatable::ln_str
  ! Computational-precision netCDF output variables
  real(selected_real_kind(p=12)),dimension(:),allocatable::A_phi
  real(selected_real_kind(p=12)),dimension(:),allocatable::A_psi
  real(selected_real_kind(p=12)),dimension(:),allocatable::B_phi
  real(selected_real_kind(p=12)),dimension(:),allocatable::B_psi
  real(selected_real_kind(p=12)),dimension(:),allocatable::S_d
  real(selected_real_kind(p=12)),dimension(:),allocatable::S_d_abs_cff_mss
  real(selected_real_kind(p=12)),dimension(:),allocatable::S_p
  real(selected_real_kind(p=12)),dimension(:),allocatable::S_p_abs_cff_mss
  real(selected_real_kind(p=12)),dimension(:),allocatable::alpha
  real(selected_real_kind(p=12)),dimension(:),allocatable::bnd ! coordinate variable
  real(selected_real_kind(p=12)),dimension(:),allocatable::frc_slr_flx_LaN68
  real(selected_real_kind(p=12)),dimension(:),allocatable::frc_slr_flx_ThD71
  real(selected_real_kind(p=12)),dimension(:),allocatable::idx_rfr_air_STP
  real(selected_real_kind(p=12)),dimension(:),allocatable::oneD_foo
  real(selected_real_kind(p=12)),dimension(:),allocatable::phi_exact
  real(selected_real_kind(p=12)),dimension(:),allocatable::phi_fit
  real(selected_real_kind(p=12)),dimension(:),allocatable::psi_exact
  real(selected_real_kind(p=12)),dimension(:),allocatable::psi_fit
  real(selected_real_kind(p=12)),dimension(:),allocatable::tpt
  real(selected_real_kind(p=12)),dimension(:),allocatable::wvl_ctr
  real(selected_real_kind(p=12)),dimension(:),allocatable::wvl_dlt
  real(selected_real_kind(p=12)),dimension(:),allocatable::wvl_grd
  real(selected_real_kind(p=12)),dimension(:),allocatable::wvl_max
  real(selected_real_kind(p=12)),dimension(:),allocatable::wvl_min
  real(selected_real_kind(p=12)),dimension(:),allocatable::wvn_ctr
  real(selected_real_kind(p=12)),dimension(:),allocatable::wvn_dlt
  real(selected_real_kind(p=12)),dimension(:),allocatable::wvn_grd
  real(selected_real_kind(p=12)),dimension(:),allocatable::wvn_max
  real(selected_real_kind(p=12)),dimension(:),allocatable::wvn_min
  ! Single-precision binary file input variables
  real(selected_real_kind(p=6)),dimension(:),allocatable::A_phi_r4
  real(selected_real_kind(p=6)),dimension(:),allocatable::A_psi_r4
  real(selected_real_kind(p=6)),dimension(:),allocatable::B_phi_r4
  real(selected_real_kind(p=6)),dimension(:),allocatable::B_psi_r4
  real(selected_real_kind(p=6)),dimension(:),allocatable::S_d_r4
  real(selected_real_kind(p=6)),dimension(:),allocatable::S_p_r4
  real(selected_real_kind(p=6)),dimension(:),allocatable::alpha_r4
  real(selected_real_kind(p=6)),dimension(:),allocatable::wvn_ctr_r4
  ! Local variables
  real(selected_real_kind(p=12)),dimension(:),allocatable::Boltzmann_wgt_tpt_rcp
  real(selected_real_kind(p=12)),dimension(:),allocatable::bnd_sum_str_ln_lmt_rfr_rcp_2
  real(selected_real_kind(p=12)),dimension(:),allocatable::bnd_sum_wk_ln_lmt_rfr_rcp
  real(selected_real_kind(p=12)),dimension(:),allocatable::prt_fnc_tpt_scl
  real(selected_real_kind(p=12)),dimension(:),allocatable::tpt_HITRAN_tpt_rcp
  real(selected_real_kind(p=12)),dimension(:),allocatable::tpt_dlt_1
  real(selected_real_kind(p=12)),dimension(:),allocatable::tpt_dlt_2
  real(selected_real_kind(p=12)),dimension(:),allocatable::tpt_dlt_3
  real(selected_real_kind(p=12)),dimension(:),allocatable::tpt_dlt_4
  real(selected_real_kind(p=12)),dimension(:),allocatable::tpt_rcp
  integer,dimension(:),allocatable::bnd_ln_nbr
  integer,dimension(:),allocatable::ln_idx_max
  integer,dimension(:),allocatable::ln_idx_min

  logical::grd_LW_SW=.true.
  integer::bnd_nbr_LW=199
  integer::bnd_idx_LW_end
  integer::bnd_nbr_SW=2500

  real(selected_real_kind(p=12))::bnd_dlt_LW=10
  real(selected_real_kind(p=12))::bnd_dlt_SW=10
  real(selected_real_kind(p=12))::wvn_min_SW
  real(selected_real_kind(p=12))::wvn_max_LW
  real(selected_real_kind(p=12))::wvn_min_LW=10
  real(selected_real_kind(p=12))::wvn_max_SW=27000
  real(selected_real_kind(p=12))::wvn_bnd_LW_SW=2000
  
  real(selected_real_kind(p=12))::Boltzmann_rcp
  real(selected_real_kind(p=12))::Boltzmann_wgt
  real(selected_real_kind(p=12))::Boltzmann_wgt_tpt
  real(selected_real_kind(p=12))::HWHM_air_rfr
  real(selected_real_kind(p=12))::HWHM_air_tpt_crr
  real(selected_real_kind(p=12))::bnd_sum_str_ln_lmt_rfr
  real(selected_real_kind(p=12))::bnd_sum_str_ln_lmt_tpt_crr
  real(selected_real_kind(p=12))::bnd_sum_wk_ln_lmt_rfr
  real(selected_real_kind(p=12))::bnd_sum_wk_ln_lmt_tpt_crr
  real(selected_real_kind(p=12))::float_foo
  real(selected_real_kind(p=12))::ln_hi
  real(selected_real_kind(p=12))::ln_lo
  real(selected_real_kind(p=12))::ln_str_rfr
  real(selected_real_kind(p=12))::ln_str_tpt_crr
  real(selected_real_kind(p=12))::mtx_a1_phi
  real(selected_real_kind(p=12))::mtx_a1_psi
  real(selected_real_kind(p=12))::mtx_a2_phi
  real(selected_real_kind(p=12))::mtx_a2_psi
  real(selected_real_kind(p=12))::mtx_b1_phi
  real(selected_real_kind(p=12))::mtx_b1_psi
  real(selected_real_kind(p=12))::mtx_b2_phi
  real(selected_real_kind(p=12))::mtx_b2_psi
  real(selected_real_kind(p=12))::mtx_c1_phi
  real(selected_real_kind(p=12))::mtx_c1_psi
  real(selected_real_kind(p=12))::mtx_c2_phi
  real(selected_real_kind(p=12))::mtx_c2_psi
  real(selected_real_kind(p=12))::phi
  real(selected_real_kind(p=12))::psi
  real(selected_real_kind(p=12))::prt_fnc_tpt_scl_rfr
  real(selected_real_kind(p=12))::stm_msn_crc
  real(selected_real_kind(p=12))::stm_msn_crc_tpt
  real(selected_real_kind(p=12))::tpt_HITRAN_rcp
  real(selected_real_kind(p=12))::tpt_Malkmus_rfr
  real(selected_real_kind(p=12))::tpt_Malkmus_rfr_rcp
  real(selected_real_kind(p=12))::tpt_max
  real(selected_real_kind(p=12))::tpt_min
  real(selected_real_kind(p=12))::tpt_ncr
  
  character(iso_sng_htrn_lng_max) iso_sng_htrn(iso_nbr_max_htrn) ! Contains all HITRAN isotopomer names
  character(mlc_sng_htrn_lng_max) mlc_sng_htrn(mlc_nbr_max_htrn) ! Contains all HITRAN molecule names
  character(iso_sng_htrn_lng_max) sng_tmp ! Temporary string 
  character*120 iso_sng
  character*120 mlc_sng
  integer map(iso_per_mlc_nbr_max_htrn,mlc_nbr_max_htrn) ! map of mlc,iso indices to isotopomer index
  real mmw_iso(iso_nbr_max_htrn) ! [kg mol-1] mean molecular weight of all HITRAN isotopomers
  real mmw_mlc(mlc_nbr_max_htrn) ! [kg mol-1] mean molecular weight of all HITRAN molecules
  real xpn_mlc(mlc_nbr_max_htrn) ! Exponent defining temperature dependence of rotational partition function
  
  ! Set defaults for command-line options 
  character(80)::drc_in='/data/zender/hitran'//nlc ! [sng] Input directory
  character(80)::drc_out='/data/zender/aca'//nlc ! [sng] Output directory
  character(80)::fl_in='H2O.nc'//nlc ! [sng] Input file
  character(80)::fl_out='mlk_H2O.nc'//nlc ! [sng] Output file
  integer::int_foo=1 ! [nbr] Integer

  ! Main code
  
  ! Initialize default values
  nc_flg=.true.
  tpt_min=180.0             ! [K]
  tpt_max=320.0             ! [K]
  
  ! Initialize options that may be overridden by command line
  bnd_dbg=800               ! Option -B
  dbg_lvl=0                 ! Option -D
  bnd_nbr=2699              ! Option -b
  ln_lo=wvn_min_LW          ! [cm-1] Option -l
  ln_hi=wvn_max_SW          ! [cm-1] Option -h
  typ_out=nf90_float        ! [enm] Output floating point type
  tpt_Malkmus_rfr=tpt_Malkmus_BPB ! Option -T
  tpt_nbr=140               ! Option -t
  
  ! Retrieve command line arguments
  call date_time_get(lcl_date_time)
  call ftn_cmd_ln_sng(cmd_ln)
  call ftn_prg_ID_mk(CVS_Id,CVS_Revision,CVS_Date,prg_ID)
  write (6,'(a)') prg_ID(1:ftn_strlen(prg_ID))
  arg_nbr=command_argument_count() ! [nbr] Number of command-line arguments
  arg_idx=1 ! [idx] Counting index
  loop_while_options: do while (arg_idx <= arg_nbr)
     call ftn_getarg_wrp(arg_idx,arg_val) ! [sbr] Call getarg(), increment arg_idx
     dsh_key=arg_val(1:2) ! [sng] First two characters of option
     if (dsh_key == '--') then
        opt_lng=ftn_opt_lng_get(arg_val) ! [nbr] Length of option
        if (opt_lng <= 0) stop 'Long option has no name'
        opt_sng=arg_val(3:2+opt_lng) ! [sng] Option string
        if (dbg_lvl >= dbg_io) write (6,'(5a,i3)') prg_nm(1:ftn_strlen(prg_nm)), &
             ': DEBUG Double hyphen indicates multi-character option: ', &
             'opt_sng = ',opt_sng(1:ftn_strlen(opt_sng)),', opt_lng = ',opt_lng
        ! fxm: Change if else if construct to select case but how to handle fall-through cases elegantly?
        if (opt_sng == 'bnd_nbr' .or. opt_sng == 'wvn_nbr' ) then
           call ftn_arg_get(arg_idx,arg_val,bnd_nbr) ! [nbr] Number of bands
           grd_LW_SW=.false.
        else if (opt_sng == 'dbg' .or. opt_sng == 'dbg_lvl' ) then
           call ftn_arg_get(arg_idx,arg_val,dbg_lvl) ! [enm] Debugging level
        else if (opt_sng == 'dbl' .or. opt_sng == 'double' ) then
           typ_out=nf90_double ! [frc] Double
        else if (opt_sng == 'dfl' .or. opt_sng == 'deflate' ) then
           call ftn_arg_get(arg_idx,arg_val,dfl_lvl) ! [enm] Deflate level
        else if (opt_sng == 'drc_in') then
           call ftn_arg_get(arg_idx,arg_val,drc_in) ! [sng] Input directory
        else if (opt_sng == 'drc_out') then
           call ftn_arg_get(arg_idx,arg_val,drc_out) ! [sng] Output directory
        else if (opt_sng == 'fl_in') then
           call ftn_arg_get(arg_idx,arg_val,fl_in) ! [sng] Input file
        else if (opt_sng == 'fl_out') then
           call ftn_arg_get(arg_idx,arg_val,fl_out) ! [sng] Output file
        else if (opt_sng == 'flt' .or. opt_sng == 'flt_foo' ) then
           typ_out=nf90_float ! [frc] Double
        else if (opt_sng == 'ln_lo' .or. opt_sng == 'wvn_min' ) then
           call ftn_arg_get(arg_idx,arg_val,ln_lo) ! [cm-1] Minimum wavenumber
        else if (opt_sng == 'ln_hi' .or. opt_sng == 'wvn_max' ) then
           call ftn_arg_get(arg_idx,arg_val,ln_hi) ! [cm-1] Maximum wavenumber
        else if (opt_sng == 'wvn_bnd' .or. opt_sng == 'wvn_bnd_LW_SW' ) then
           call ftn_arg_get(arg_idx,arg_val,wvn_bnd_LW_SW) ! [cm-1] Boundary wavenumber between LW and SW grids
        else if (opt_sng == 'rsn_LW' .or. opt_sng == 'bnd_dlt_LW' ) then
           call ftn_arg_get(arg_idx,arg_val,bnd_dlt_LW) ! [cm-1] Resolution in LW portion of grid
        else if (opt_sng == 'rsn_SW' .or. opt_sng == 'bnd_dlt_SW' ) then
           call ftn_arg_get(arg_idx,arg_val,bnd_dlt_SW) ! [cm-1] Resolution in SW portion of grid
        else ! Option not recognized
           arg_idx=arg_idx-1 ! [idx] Counting index
           call ftn_getarg_err(arg_idx,arg_val) ! [sbr] Error handler for getarg()
        endif ! endif option is recognized
        ! Jump to top of while loop
        cycle loop_while_options ! C, F77, and F90 use "continue", "goto", and "cycle"
     else if (dsh_key(1:1) == '-') then ! endif long option
        ! Handle short options
        if (dsh_key == '-3') then
           fl_out_fmt=nf90_format_classic ! [enm] Output file format
        else if (dsh_key == '-4') then
           fl_out_fmt=nf90_format_netcdf4 ! [enm] Output file format
        else if (dsh_key == '-B') then
           call ftn_arg_get(arg_idx,arg_val,bnd_dbg)
        else if (dsh_key == '-b') then
           call ftn_arg_get(arg_idx,arg_val,bnd_nbr)
           grd_LW_SW=.false.
        else if (dsh_key == '-D') then
           call ftn_arg_get(arg_idx,arg_val,dbg_lvl)
        else if (dsh_key == '-f') then
           call ftn_arg_get(arg_idx,arg_val,float_foo)
        else if (dsh_key == '-h') then
           call ftn_arg_get(arg_idx,arg_val,ln_hi)
        else if (dsh_key == '-i') then
           call ftn_arg_get(arg_idx,arg_val,fl_in)
        else if (dsh_key == '-l') then
           call ftn_arg_get(arg_idx,arg_val,ln_lo)
        else if (dsh_key == '-n') then
           nc_flg=.not.nc_flg
        else if (dsh_key == '-o') then
           call ftn_arg_get(arg_idx,arg_val,fl_out)
        else if (dsh_key == '-t') then
           call ftn_arg_get(arg_idx,arg_val,tpt_nbr)
        else if (dsh_key == '-T') then
           call ftn_arg_get(arg_idx,arg_val,tpt_Malkmus_rfr)
        else if (dsh_key == '-v') then
           write (6,'(a)') CVS_Id
           goto 1000 ! Goto exit with error status
        else ! Option not recognized
           arg_idx=arg_idx-1 ! [idx] Counting index
           call ftn_getarg_err(arg_idx,arg_val) ! [sbr] Error handler for getarg()
        endif ! endif arg_val
     else ! endif short option
        ! Last argument(s) may be positional
        if (arg_idx == arg_nbr) then
           call ftn_strcpylsc(fl_in,arg_val) ! [sng] Input file
        else if (arg_idx == arg_nbr+1) then
           call ftn_strcpylsc(fl_out,arg_val) ! [sng] Output file
        else ! Option not recognized
           arg_idx=arg_idx-1 ! [idx] Counting index
           call ftn_getarg_err(arg_idx,arg_val) ! [sbr] Error handler for getarg()
        endif ! arg_idx != arg_nbr
     end if ! end position arguments
  end do loop_while_options ! end while (arg_idx <= arg_nbr)
  
  ! Compute quantities that depend on command line input
  call ftn_strnul(fl_in)
  call ftn_strnul(fl_out)

  ! Compute any quantities that might depend on command-line input
  ! Prepend user-specified path, if any, to input data file names
  if (ftn_strlen(drc_in) > 0) call ftn_drcpfx(drc_in,fl_in) ! [sng] Input file
  ! Prepend user-specified path, if any, to output data file names
  if (ftn_strlen(drc_out) > 0) call ftn_drcpfx(drc_out,fl_out) ! [sng] Output file

  ! Enough memory?
  if (tpt_nbr > tpt_nbr_max) stop 'tpt_nbr > tpt_nbr_max'
  if (bnd_nbr > bnd_nbr_max) stop 'bnd_nbr > bnd_nbr_max'
  
  if(grd_LW_SW) then
     wvn_min_LW=ln_lo
     wvn_max_SW=ln_hi
     wvn_min_SW=wvn_bnd_LW_SW
     wvn_max_LW=wvn_bnd_LW_SW
     bnd_nbr_LW=nint((wvn_max_LW-wvn_min_LW)/bnd_dlt_LW)
     bnd_nbr_SW=nint((wvn_max_SW-wvn_min_SW)/bnd_dlt_SW)
     bnd_nbr=bnd_nbr_LW+bnd_nbr_SW
     bnd_idx_LW_end=bnd_nbr_LW
  endif ! !grd_LW_SW

  ! Computational-precision netCDF output variables
  allocate(A_phi(bnd_nbr),stat=rcd)
  allocate(A_psi(bnd_nbr),stat=rcd)
  allocate(B_phi(bnd_nbr),stat=rcd)
  allocate(B_psi(bnd_nbr),stat=rcd)
  allocate(S_d(bnd_nbr),stat=rcd)
  allocate(S_d_abs_cff_mss(bnd_nbr),stat=rcd)
  allocate(S_p(bnd_nbr),stat=rcd)
  allocate(S_p_abs_cff_mss(bnd_nbr),stat=rcd)
  allocate(alpha(bnd_nbr),stat=rcd)
  allocate(bnd(bnd_nbr),stat=rcd) ! coordinate variable
  allocate(frc_slr_flx_LaN68(bnd_nbr),stat=rcd)
  allocate(frc_slr_flx_ThD71(bnd_nbr),stat=rcd)
  allocate(idx_rfr_air_STP(bnd_nbr),stat=rcd)
  
  allocate(oneD_foo(tpt_nbr),stat=rcd)
  allocate(phi_exact(tpt_nbr),stat=rcd)
  allocate(phi_fit(tpt_nbr),stat=rcd)
  allocate(psi_exact(tpt_nbr),stat=rcd)
  allocate(psi_fit(tpt_nbr),stat=rcd)
  allocate(tpt(tpt_nbr),stat=rcd)
  allocate(wvl_ctr(bnd_nbr),stat=rcd)
  allocate(wvl_dlt(bnd_nbr),stat=rcd)
  allocate(wvl_grd(bnd_nbr+1),stat=rcd)
  allocate(wvl_max(bnd_nbr),stat=rcd)
  allocate(wvl_min(bnd_nbr),stat=rcd)
  allocate(wvn_ctr(bnd_nbr),stat=rcd)
  allocate(wvn_dlt(bnd_nbr),stat=rcd)
  allocate(wvn_grd(bnd_nbr+1),stat=rcd)
  allocate(wvn_max(bnd_nbr),stat=rcd)
  allocate(wvn_min(bnd_nbr),stat=rcd)
  ! Local variables
  allocate(bnd_ln_nbr(bnd_nbr),stat=rcd)
  allocate(ln_idx_max(bnd_nbr),stat=rcd)
  allocate(ln_idx_min(bnd_nbr),stat=rcd)
  ! Local variables
  allocate(Boltzmann_wgt_tpt_rcp(tpt_nbr),stat=rcd)
  allocate(bnd_sum_str_ln_lmt_rfr_rcp_2(bnd_nbr),stat=rcd)
  allocate(bnd_sum_wk_ln_lmt_rfr_rcp(bnd_nbr),stat=rcd)
  allocate(prt_fnc_tpt_scl(tpt_nbr),stat=rcd)
  allocate(tpt_HITRAN_tpt_rcp(tpt_nbr),stat=rcd)
  allocate(tpt_dlt_1(tpt_nbr),stat=rcd)
  allocate(tpt_dlt_2(tpt_nbr),stat=rcd)
  allocate(tpt_dlt_3(tpt_nbr),stat=rcd)
  allocate(tpt_dlt_4(tpt_nbr),stat=rcd)
  allocate(tpt_rcp(tpt_nbr),stat=rcd)
  
  ! Compute any quantities that might depend on command line input
  call ftn_strnul(fl_in)
  call ftn_strnul(fl_out)
  
  if (nc_flg) then 
     write(6,'(2a)') prg_nm(1:ftn_strlen(prg_nm)), &
          ': INFO Input presumed to be netCDF format HITRAN line data'
  else 
     write(6,'(2a)') prg_nm(1:ftn_strlen(prg_nm)), &
          ': INFO Input presumed to be IEEE binary BPB-format Malkmus random band data' 
  endif                     ! endif
  
  if (nc_flg) then
     
     ! Open file
     rcd=nf90_wrp_open(fl_in,nf90_nowrite,nc_id,sbr_nm=sbr_nm)
     ! Get dimension IDs
     rcd=nf90_wrp_inq_dimid(nc_id,'iso',iso_dmn_id)
     rcd=nf90_wrp_inq_dimid(nc_id,'ln_ctr',ln_ctr_dmn_id)
     rcd=nf90_wrp_inq_dimid(nc_id,'mlc',mlc_dmn_id)
     ! Get dimension sizes
     rcd=rcd+nf90_inquire_dimension(nc_id,ln_ctr_dmn_id,len=ln_nbr)
     ! Enough memory?
     if (ln_nbr > ln_nbr_max) stop 'ln_nbr > ln_nbr_max'
     rcd=rcd+nf90_inquire_dimension(nc_id,mlc_dmn_id,len=mlc_nbr)
     if (mlc_nbr > mlc_nbr_max_htrn) stop 'mlc_nbr > mlc_nbr_max_htrn'
     rcd=rcd+nf90_inquire_dimension(nc_id,iso_dmn_id,len=iso_nbr)
     if (iso_nbr > iso_nbr_max_htrn) stop 'iso_nbr > iso_nbr_max_htrn'
     srt_one(1)=1
     cnt_ln(1)=ln_nbr
     cnt_mlc(1)=mlc_nbr
     cnt_iso(1)=iso_nbr
     
     ! Computational-precision HITRAN netCDF input variables
     allocate(HWHM_air(ln_nbr),stat=rcd)
     allocate(HWHM_tpt_dpn_xpn(ln_nbr),stat=rcd)
     allocate(ln_ctr(ln_nbr),stat=rcd)
     allocate(ln_nrg_lwr(ln_nbr),stat=rcd)
     allocate(ln_str(ln_nbr),stat=rcd)
     
     ! htrn2nb requires that only one molecule be present in an input file
     if (mlc_nbr /= 1) then
        write(6,'(a,i3,a)') 'htrn2nb: ERROR ',mlc_nbr,' molecules in input file'
        rcd=-1
        call exit(rcd)
     endif
     ! Get variable IDs
     rcd=nf90_wrp_inq_varid(nc_id,'HWHM_air',HWHM_air_id)
     rcd=nf90_wrp_inq_varid(nc_id,'HWHM_tpt_dpn_xpn',HWHM_tpt_dpn_xpn_id)
     rcd=nf90_wrp_inq_varid(nc_id,'ln_iso',ln_iso_id)
     rcd=nf90_wrp_inq_varid(nc_id,'ln_ctr',ln_ctr_id)
     rcd=nf90_wrp_inq_varid(nc_id,'ln_nrg_lwr',ln_nrg_lwr_id)
     rcd=nf90_wrp_inq_varid(nc_id,'ln_str',ln_str_id)
     rcd=nf90_wrp_inq_varid(nc_id,'mlc_id',mlc_id_id)
     rcd=nf90_wrp_inq_varid(nc_id,'iso_id',iso_id_id)
     ! Get data
     rcd=nf90_wrp(nf90_get_var(nc_id,mlc_id_id,mlc_id),'gv mlc_id')
     rcd=nf90_wrp(nf90_get_var(nc_id,iso_id_id,iso_id,srt_one,cnt_iso),'gv iso_id')
     rcd=nf90_wrp(nf90_get_var(nc_id,HWHM_air_id,HWHM_air,srt_one,cnt_ln),'gv HWHM_air')
     rcd=nf90_wrp(nf90_get_var(nc_id,HWHM_tpt_dpn_xpn_id,HWHM_tpt_dpn_xpn,srt_one,cnt_ln),'gv HWHM_tpt_dpn_xpn')
     rcd=nf90_wrp(nf90_get_var(nc_id,ln_ctr_id,ln_ctr,srt_one,cnt_ln),'gv ln_ctr')
     rcd=nf90_wrp(nf90_get_var(nc_id,ln_nrg_lwr_id,ln_nrg_lwr,srt_one,cnt_ln),'gv ln_nrg_lwr')
     rcd=nf90_wrp(nf90_get_var(nc_id,ln_str_id,ln_str,srt_one,cnt_ln),'gv ln_str')
     ! Get global attributes
     rcd=rcd+nf90_get_att(nc_id,nf90_global,'molecule',mlc_sng)
     rcd=rcd+nf90_get_att(nc_id,nf90_global,'isotope',iso_sng)
     call ftn_strnul(mlc_sng)
     call ftn_strnul(iso_sng)
     ! Close file
     rcd=nf90_wrp_close(nc_id,fl_in,'Ingested')
     
     ! Convert input data to SI units where necessary 
     ! Normally I never store anything in units other than SI in a netCDF file, but HITRAN input is an exception 
     ! Line data is exactly as it appears in the HITRAN database, which is basically CGS 
     ! ln_str is used in all sorts of SI algebra involving fundamental constants
     ! This algebra is too confusing to do in any units besides SI
     ! Thus ln_str should be scaled to SI, even though this can result in a loss of precision for weak lines in single precision
     ! Putting ln_str in SI ensures the band parameters (S_d, S_p, ...) computed from ln_str will be in SI
     do ln_idx=1,ln_nbr
        ln_str(ln_idx)=ln_str(ln_idx)*1.0e-4 ! [cm-1 cm2 mlc-1] -> [cm-1 m2 mlc-1]
        ln_str(ln_idx)=ln_str(ln_idx)*Avagadro ! [cm-1 m2 mlc-1] -> [cm-1 m2 mol-1]
     enddo                  ! end loop over lines
     
     ! Initialize arrays of molecular properties
     call rtl_fnc_tpt_xpn_get(xpn_mlc)
     call mmw_iso_get(mmw_iso)
     call mmw_mlc_get(mmw_mlc)
     call iso_idx_map_get(map)
     call mlc_sng_get(mlc_sng_htrn)
     call iso_sng_get(iso_sng_htrn)
     
     ! Print molecular diagnostics
     write(6,'(a10,a)') 'mlc_sng = ',mlc_sng(1:ftn_strlen(mlc_sng))
     write(6,'(a10,a)') 'iso_sng = ',iso_sng(1:ftn_strlen(iso_sng))
     write(6,'(a9,i2,a3,a,a8,es8.1,a21,f3.1)') 'Molecule ',mlc_id,' = ',mlc_sng(1:ftn_strlen(mlc_sng)), &
          ', mmw = ',mmw_mlc(mlc_id),' kg mol-1, xpn_mlc = ',xpn_mlc(mlc_id)
     do idx=1,iso_nbr
        sng_tmp=iso_sng_htrn(iso_id(idx))
        call ftn_strnul(sng_tmp)
        write(6,'(a11,i2,a3,a,a8,es8.1,a9)') 'Isotopomer',iso_id(idx),' = ',sng_tmp(1:ftn_strlen(sng_tmp)), &
             ', mmw = ',mmw_iso(iso_id(idx)),' kg mol-1'
     enddo                  ! end loop over iso
     
     ! Compute wavenumber coordinates
     if(grd_LW_SW) then
        ! Split grid with distinct resolutions in SW/LW regions
        do bnd_idx=1,bnd_nbr_LW
           wvn_grd(bnd_idx)=wvn_min_LW+(bnd_idx-1)*bnd_dlt_LW
        enddo                  ! end loop over bnd
        do bnd_idx=bnd_nbr_LW+1,bnd_nbr_LW+bnd_nbr_SW
           wvn_grd(bnd_idx)=wvn_min_SW+(bnd_idx-bnd_nbr_LW-1)*bnd_dlt_SW
        enddo                  ! end loop over bnd
        wvn_grd(bnd_nbr+1)=wvn_max_SW
     else ! !grd_LW_SW
        ! Normal grid with uniform resolution in wvn space
        bnd_dlt=(ln_hi-ln_lo)/bnd_nbr
        float_foo=bnd_dlt/2.0
        do bnd_idx=1,bnd_nbr
           wvn_grd(bnd_idx)=ln_lo+(bnd_idx-1)*bnd_dlt
        enddo                  ! end loop over bnd
        wvn_grd(bnd_nbr+1)=ln_hi
     endif ! !grd_LW_SW

     if (dbg_lvl >= dbg_crr) then
        do bnd_idx=1,bnd_nbr+1 
           write (6,*) 'idx = ',bnd_idx,'wvn_grd = ',wvn_grd(bnd_idx)
        enddo                  ! end loop over bnd
     endif ! endif dbg
     
     do bnd_idx=1,bnd_nbr
        wvn_min(bnd_idx)=wvn_grd(bnd_idx)
        wvn_max(bnd_idx)=wvn_grd(bnd_idx+1)
     enddo                  ! end loop over bnd
     do bnd_idx=1,bnd_nbr
        wvn_ctr(bnd_idx)=0.5*(wvn_min(bnd_idx)+wvn_max(bnd_idx))
     enddo                  ! end loop over bnd
     
     ! Initialize output parameters: 
     ! Empty bands (bands containing zero lines) will contain these initialized values on output 
     ! When there are no lines in a band, the S_d parameter should be 0.0
     ! In this case S_d = 0.0 causes transmission to become 1.0 as long as S_p != 0.0 (which would cause a divide-by-zero errors) 
     ! Thus, set S_p = 1.0 by default
     do bnd_idx=1,bnd_nbr
        S_d(bnd_idx)=0.0
        S_p(bnd_idx)=1.0
        A_phi(bnd_idx)=0.0
        B_phi(bnd_idx)=0.0
        A_psi(bnd_idx)=0.0
        B_psi(bnd_idx)=0.0
     enddo                  ! end loop over bnd
     
     do bnd_idx=1,bnd_nbr
        ln_idx_min(bnd_idx)=0
        ln_idx_max(bnd_idx)=-1
     enddo                  ! end loop over bnd
     
     if (dbg_lvl > dbg_fl) then
        write (6,'(2(a8,f9.2,a5,/),a10,i6,/,a10,i4)')  &
             'ln_lo = ',ln_lo,' cm-1', &
             'ln_hi = ',ln_hi,' cm-1', &
             'ln_nbr = ',ln_nbr, &
             'bnd_nbr = ',bnd_nbr
     endif                     ! endif dbg
     
     ! Find first line within specified spectral interval
     ! This index bootstraps next loop
     do ln_idx=1,ln_nbr
        if (ln_ctr(ln_idx) >= ln_lo.and.ln_ctr(ln_idx) <= ln_hi) goto 500
        if (ln_idx==ln_nbr) then
           write(6,'(a41,/,a27,g15.8,a5,/,a26,g15.8,a5,a8)')  &
                'No lines found within specified interval.', &
                'First line in input data = ',ln_ctr(1),' cm-1', &
                'Last line in input data = ',ln_ctr(ln_nbr),' cm-1', &
                'Exiting.'
           rcd=-1
           call exit(rcd)
        endif
     enddo                  ! end loop over lines
500  continue
     
     ! Find first and last lines in each band. We assume that on entry ln_idx points
     ! to beginning of first band, and that ln is monotonically increasing.
     do bnd_idx=1,bnd_nbr
        do while (ln_idx <= ln_nbr.and.ln_ctr(ln_idx) >= wvn_min(bnd_idx).and.ln_ctr(ln_idx) < wvn_max(bnd_idx))
           if (ln_idx_min(bnd_idx) == 0) ln_idx_min(bnd_idx)=ln_idx
           ln_idx_max(bnd_idx)=ln_idx
           ln_idx=ln_idx+1
        end do              ! end loop over lines within each band
     enddo                  ! end loop over bnd
     int_foo=0
     do bnd_idx=1,bnd_nbr
        bnd_ln_nbr(bnd_idx)=ln_idx_max(bnd_idx)-ln_idx_min(bnd_idx)+1
        int_foo=int_foo+bnd_ln_nbr(bnd_idx)
     enddo                  ! end loop over bnd
     
     ! Now we know which lines are in each band
     ! For each band compute temperature-dependent band parameters in preparation for least-squares fit
     tpt_ncr=(tpt_max-tpt_min)/tpt_nbr
     tpt(1)=tpt_min
     do tpt_idx=2,tpt_nbr
        tpt(tpt_idx)=tpt(tpt_idx-1)+tpt_ncr
     enddo                  ! end loop over temperatures
     
     if(grd_LW_SW) then 
        write (6,'(a51,f12.6,a9,f12.6,a10,i5,a9,f12.6,a17)')  &
             'Malkmus random band model parameters computed from ', &
             wvn_min_LW,' cm-1 to ',wvn_max_LW,' cm-1 for ',bnd_nbr_LW,' regular ',bnd_dlt_LW,' cm-1 bands in LW'
        write (6,'(a51,f12.6,a9,f12.6,a10,i5,a9,f12.6,a17)')  &
             'Malkmus random band model parameters computed from ', &
             wvn_min_SW,' cm-1 to ',wvn_max_SW,' cm-1 for ',bnd_nbr_SW,' regular ',bnd_dlt_SW,' cm-1 bands in SW'
     else
        write (6,'(a51,f12.6,a9,f12.6,a10,i5,a9,f12.6,a11)')  &
             'Malkmus random band model parameters computed from ', &
             ln_lo,' cm-1 to ',ln_hi,' cm-1 for ',bnd_nbr_LW,' regular ',bnd_dlt,' cm-1 bands'
     endif
     write (6,'(i6,a46)') int_foo,' lines fall within specified spectral interval'
     write (6,'(a63,f6.2,a2)')  &
          'Line strength parameters scaled to and saved at tpt(reference) = ',tpt_Malkmus_rfr,' K'
     write (6,'(a75,i3,a19,f6.2,a6,f6.2,a9,f5.2,a2)')  &
          'Temperature dependence accounted for by least-squares-fit to exact data for ', &
          tpt_nbr,' temperatures from ',tpt_min,' K to ',tpt_max,' K every ',tpt_ncr,' K'
     write (6,'(a76,i4,a3,f12.6,a8,f12.6,a3)')  &
          'Fitted and exact weak and strong line temperature dependence saved for band ',bnd_dbg, &
          ' = ',wvn_ctr(bnd_dbg),' cm-1 = ',wvl_ctr(bnd_dbg)*1.0e6,' um'
     
     ! Precompute some frequently used factors and/or their reciprocals
     Boltzmann_rcp=1.0/Boltzmann
     tpt_HITRAN_rcp=1.0/tpt_HITRAN
     tpt_Malkmus_rfr_rcp=1.0/tpt_Malkmus_rfr
     do tpt_idx=1,tpt_nbr
        tpt_rcp(tpt_idx)=1.0/tpt(tpt_idx)
        prt_fnc_tpt_scl(tpt_idx)=(tpt_HITRAN*tpt_rcp(tpt_idx))**xpn_mlc(mlc_id) ! Lio92 p. 33 (2.2.21)
        Boltzmann_wgt_tpt_rcp(tpt_idx)=tpt_HITRAN_rcp-tpt_rcp(tpt_idx)
        tpt_HITRAN_tpt_rcp(tpt_idx)=tpt_HITRAN*tpt_rcp(tpt_idx)
        tpt_dlt_1(tpt_idx)=tpt(tpt_idx)-tpt_Malkmus_rfr
        tpt_dlt_2(tpt_idx)=tpt_dlt_1(tpt_idx)*tpt_dlt_1(tpt_idx)
        tpt_dlt_3(tpt_idx)=tpt_dlt_2(tpt_idx)*tpt_dlt_1(tpt_idx)
        tpt_dlt_4(tpt_idx)=tpt_dlt_3(tpt_idx)*tpt_dlt_1(tpt_idx)
     enddo                  ! end loop over temperatures
     
     ! Least-squares fit matrix coefficients on LHS are purely temperature-dependent
     ! Linear system matrices look like
     ! a1*x + b1*y = c1 
     ! a2*x + b2*y = c2
     ! where x and y are least-squares solutions to: phi(t) = exp[x*(t-t0) + y*(t-t0)**2]
     
     mtx_a1_phi=0.0
     mtx_b1_phi=0.0
     mtx_a2_phi=0.0
     mtx_b2_phi=0.0
     mtx_a1_psi=0.0
     mtx_b1_psi=0.0
     mtx_a2_psi=0.0
     mtx_b2_psi=0.0
     do tpt_idx=1,tpt_nbr
        mtx_a1_phi=mtx_a1_phi+tpt_dlt_2(tpt_idx)
        mtx_b1_phi=mtx_b1_phi+tpt_dlt_3(tpt_idx)
        mtx_a2_phi=mtx_a2_phi+tpt_dlt_3(tpt_idx)
        mtx_b2_phi=mtx_b2_phi+tpt_dlt_4(tpt_idx)
        mtx_a1_psi=mtx_a1_psi+tpt_dlt_2(tpt_idx)
        mtx_b1_psi=mtx_b1_psi+tpt_dlt_3(tpt_idx)
        mtx_a2_psi=mtx_a2_psi+tpt_dlt_3(tpt_idx)
        mtx_b2_psi=mtx_b2_psi+tpt_dlt_4(tpt_idx)
     enddo                  ! end loop over temperatures
     
     ! Adjust all HITRAN temperature (296 K) line strengths to temperature
     ! nearer mass-weighted mean temperature of atmosphere (tpt_Malkmus_rfr) 
     ! for use in temperature fit to line strengths and for later storage 
     ! and reference. 
     ! This block computes reference band parameters BPB calls S_d and S_p in Bri92

     ! Suffix _rfr pertains to quantity evaulated at reference (tpt_Malkmus_rfr)
     ! temperature (250 K) at which band parameters and least-squares fit parameters
     ! saved to disk refer to.
     ! 
     ! Suffix _std pertains to quantity evaulated at standard HITRAN (tpt_HITRAN)
     ! temperature (296 K) of input line parameters. _std quantities are used
     ! to compute temperature-dependent band parameters in following loop.
     prt_fnc_tpt_scl_rfr=(tpt_HITRAN/tpt_Malkmus_rfr)**xpn_mlc(mlc_id)
     do bnd_idx=1,bnd_nbr
        if (bnd_ln_nbr(bnd_idx) > 0) then
           if (grd_LW_SW) then
              if (bnd_idx <= bnd_idx_LW_end) then
                 bnd_dlt=bnd_dlt_LW
              else
                 bnd_dlt=bnd_dlt_SW
              endif
           endif ! !grd_LW_SW

           bnd_sum_wk_ln_lmt_rfr=0.0
           bnd_sum_str_ln_lmt_rfr=0.0
           do ln_idx=ln_idx_min(bnd_idx),ln_idx_max(bnd_idx)
              
              ! NB: Negative sign in Boltzmann_wgt_tpt exponent has been taken
              ! care of by reversing temperatures in Boltzmann_wgt_tpt_rcp. 
              ! Factor of 100.0*speed_of_light converts from HITRAN wavenumbers to SI frequency in Hertz
              
              Boltzmann_wgt_tpt=100.0*Planck*speed_of_light*ln_nrg_lwr(ln_idx)*Boltzmann_rcp
              stm_msn_crc_tpt=-100.0*Planck*speed_of_light*ln_ctr(ln_idx)*Boltzmann_rcp
              Boltzmann_wgt=exp(Boltzmann_wgt_tpt*(tpt_HITRAN_rcp-tpt_Malkmus_rfr_rcp))
              stm_msn_crc= &
                   (1.0-exp(stm_msn_crc_tpt*tpt_Malkmus_rfr_rcp))/ &
                   (1.0-exp(stm_msn_crc_tpt*tpt_HITRAN_rcp))
              !              stm_msn_crc=
              ! $                 expm1(stm_msn_crc_tpt*tpt_Malkmus_rfr_rcp)/
              ! $                 expm1(stm_msn_crc_tpt*tpt_HITRAN_rcp)
              ln_str_rfr=ln_str(ln_idx)*prt_fnc_tpt_scl_rfr* &
                   Boltzmann_wgt*stm_msn_crc
              HWHM_air_rfr=HWHM_air(ln_idx)*(tpt_HITRAN*tpt_Malkmus_rfr_rcp)**HWHM_tpt_dpn_xpn(ln_idx)
              
              bnd_sum_wk_ln_lmt_rfr=bnd_sum_wk_ln_lmt_rfr+ln_str_rfr
              bnd_sum_str_ln_lmt_rfr=bnd_sum_str_ln_lmt_rfr+ &
                   sqrt(ln_str_rfr*HWHM_air_rfr)
              
           end do           ! end loop over lines within each band
           bnd_sum_wk_ln_lmt_rfr_rcp(bnd_idx)= &
                1.0/bnd_sum_wk_ln_lmt_rfr
           bnd_sum_str_ln_lmt_rfr_rcp_2(bnd_idx)= &
                1.0/(bnd_sum_str_ln_lmt_rfr*bnd_sum_str_ln_lmt_rfr)
           
           S_d(bnd_idx)=bnd_sum_wk_ln_lmt_rfr/bnd_dlt
           S_p(bnd_idx)=bnd_sum_wk_ln_lmt_rfr*bnd_sum_wk_ln_lmt_rfr/ &
                (4.0*bnd_sum_str_ln_lmt_rfr*bnd_sum_str_ln_lmt_rfr)
           
           ! Rescale data back to per molecule now that computations are complete
           S_d(bnd_idx)=S_d(bnd_idx)/Avagadro ! [m2 mol-1] -> [m2 mlc-1]
           S_p(bnd_idx)=S_p(bnd_idx)/Avagadro ! [m2 mol-1] -> [m2 mlc-1]
           
        endif               ! end if there are lines in this band
     enddo                  ! end loop over bnd
     
     ! We now have band mean line strength (S_d) and band pressure weighted
     ! strong line limit strength (S_p) computed for reference temperature.
     ! Compute analogous properties (phi and psi) over desired temperature range
     ! of quadratic-in-temperature parameterization.
     
     ! Suffix _tpt_crr pertains to quantities evaluated at/scaled to the current 
     ! temperature in the temperature range over which we are parameterizing band parameters
     do bnd_idx=1,bnd_nbr
        write (6,'(a1)',advance='no') '.'
        if (bnd_ln_nbr(bnd_idx) > 0) then
           mtx_c1_phi=0.0
           mtx_c2_phi=0.0
           mtx_c1_psi=0.0
           mtx_c2_psi=0.0
           do tpt_idx=1,tpt_nbr
              bnd_sum_wk_ln_lmt_tpt_crr=0.0
              bnd_sum_str_ln_lmt_tpt_crr=0.0
              do ln_idx=ln_idx_min(bnd_idx),ln_idx_max(bnd_idx)
                 
                 ! NB: Negative sign in Boltzmann_wgt_tpt exponent has been taken
                 ! care of by reversing temperatures in Boltzmann_wgt_tpt_rcp. 
                 ! Factor of 100.0*speed_of_light converts from HITRAN wavenumbers to SI frequency in Hertz
                 
                 Boltzmann_wgt_tpt=100.0*Planck*speed_of_light*ln_nrg_lwr(ln_idx)*Boltzmann_rcp
                 stm_msn_crc_tpt=-100.0*Planck*speed_of_light*ln_ctr(ln_idx)*Boltzmann_rcp
                 Boltzmann_wgt=exp(Boltzmann_wgt_tpt*Boltzmann_wgt_tpt_rcp(tpt_idx))
                 stm_msn_crc= &
                      (1.0-exp(stm_msn_crc_tpt*tpt_rcp(tpt_idx)))/ &
                      (1.0-exp(stm_msn_crc_tpt*tpt_HITRAN_rcp))
                 ln_str_tpt_crr=ln_str(ln_idx)*prt_fnc_tpt_scl(tpt_idx)*Boltzmann_wgt*stm_msn_crc
                 HWHM_air_tpt_crr=HWHM_air(ln_idx)*tpt_HITRAN_tpt_rcp(tpt_idx)**HWHM_tpt_dpn_xpn(ln_idx)
                 
                 bnd_sum_wk_ln_lmt_tpt_crr=bnd_sum_wk_ln_lmt_tpt_crr+ln_str_tpt_crr
                 bnd_sum_str_ln_lmt_tpt_crr=bnd_sum_str_ln_lmt_tpt_crr+ &
                      sqrt(ln_str_tpt_crr*HWHM_air_tpt_crr)
                 
              end do        ! end loop over lines within each band
              
              phi=bnd_sum_wk_ln_lmt_tpt_crr*bnd_sum_wk_ln_lmt_rfr_rcp(bnd_idx)
              psi=bnd_sum_str_ln_lmt_tpt_crr*bnd_sum_str_ln_lmt_tpt_crr* &
                   bnd_sum_str_ln_lmt_rfr_rcp_2(bnd_idx)
              
              mtx_c1_phi=mtx_c1_phi+ &
                   log(phi)*tpt_dlt_1(tpt_idx)
              mtx_c2_phi=mtx_c2_phi+ &
                   log(phi)*tpt_dlt_2(tpt_idx)
              mtx_c1_psi=mtx_c1_psi+ &
                   log(psi)*tpt_dlt_1(tpt_idx)
              mtx_c2_psi=mtx_c2_psi+ &
                   log(psi)*tpt_dlt_2(tpt_idx)
              
              if (bnd_idx==bnd_dbg) then
                 phi_exact(tpt_idx)=phi
                 psi_exact(tpt_idx)=psi
              endif         ! end if
              
           enddo            ! end loop over temperatures
           A_phi(bnd_idx)= &
                (mtx_b2_phi*mtx_c1_phi-mtx_b1_phi*mtx_c2_phi)/ &
                (mtx_a1_phi*mtx_b2_phi-mtx_a2_phi*mtx_b1_phi)
           B_phi(bnd_idx)= &
                (mtx_a2_phi*mtx_c1_phi-mtx_a1_phi*mtx_c2_phi)/ &
                (mtx_a2_phi*mtx_b1_phi-mtx_a1_phi*mtx_b2_phi)
           
           A_psi(bnd_idx)= &
                (mtx_b2_psi*mtx_c1_psi-mtx_b1_psi*mtx_c2_psi)/ &
                (mtx_a1_psi*mtx_b2_psi-mtx_a2_psi*mtx_b1_psi)
           B_psi(bnd_idx)= &
                (mtx_a2_psi*mtx_c1_psi-mtx_a1_psi*mtx_c2_psi)/ &
                (mtx_a2_psi*mtx_b1_psi-mtx_a1_psi*mtx_b2_psi)
           
        endif               ! end if there are lines in this band
        
     enddo                  ! end loop over bnd
     write (6,'(/)')
     
     ! De-allocate dynamic variables
     ! Computational-precision HITRAN netCDF input variables
     if (allocated(HWHM_air)) deallocate(HWHM_air,stat=rcd)
     if (allocated(HWHM_tpt_dpn_xpn)) deallocate(HWHM_tpt_dpn_xpn,stat=rcd)
     if (allocated(ln_ctr)) deallocate(ln_ctr,stat=rcd)
     if (allocated(ln_nrg_lwr)) deallocate(ln_nrg_lwr,stat=rcd)
     if (allocated(ln_str)) deallocate(ln_str,stat=rcd)
     
  else                      ! endif input file is netCDF format
     
     ! Single-precision binary file input variables
     allocate(A_phi_r4(bnd_nbr),stat=rcd)
     allocate(A_psi_r4(bnd_nbr),stat=rcd)
     allocate(B_phi_r4(bnd_nbr),stat=rcd)
     allocate(B_psi_r4(bnd_nbr),stat=rcd)
     allocate(S_d_r4(bnd_nbr),stat=rcd)
     allocate(S_p_r4(bnd_nbr),stat=rcd)
     allocate(alpha_r4(bnd_nbr),stat=rcd)
     allocate(wvn_ctr_r4(bnd_nbr),stat=rcd)
     
     ! Read input quantities
     open (fl_in_unit,file=fl_in,form='unformatted',status='old',iostat=rcd,access='sequential')
     
     read (fl_in_unit) ( &
          wvn_ctr_r4(bnd_idx), &
          S_d_r4(bnd_idx), &
          S_p_r4(bnd_idx), &
          A_phi_r4(bnd_idx), &
          B_phi_r4(bnd_idx), &
          A_psi_r4(bnd_idx), &
          B_psi_r4(bnd_idx), &
          alpha_r4(bnd_idx), &
          bnd_idx=1,bnd_nbr)
     
     close (fl_in_unit)
     write (6,'(a20,1x,a40)') 'Read input data from',fl_in
     
     ! Convert input data to SI units where necessary
     do bnd_idx=1,bnd_nbr
        S_d_r4(bnd_idx)=S_d_r4(bnd_idx)*1.0e-4_r4 ! [cm2 mlc-1] -> [m2 mlc-1]
        S_p_r4(bnd_idx)=S_p_r4(bnd_idx)*1.0e-4_r4 ! [cm2 mlc-1] -> [m2 mlc-1]
     enddo                  ! end loop over bnd
     
     ! Move temporary storage into permanant storage to facilitate
     ! auto-doubling this routine with compiler -r8 option.
     ! Without these shenanigans, autodoubling would fail on binary read statement
     do bnd_idx=1,bnd_nbr
        wvn_ctr(bnd_idx)=wvn_ctr_r4(bnd_idx)
        S_d(bnd_idx)=S_d_r4(bnd_idx)
        S_p(bnd_idx)=S_p_r4(bnd_idx)
        A_phi(bnd_idx)=A_phi_r4(bnd_idx)
        B_phi(bnd_idx)=B_phi_r4(bnd_idx)
        A_psi(bnd_idx)=A_psi_r4(bnd_idx)
        B_psi(bnd_idx)=B_psi_r4(bnd_idx)
        alpha(bnd_idx)=alpha_r4(bnd_idx)
     enddo                  ! end loop over bnd
     
     ! De-allocate dynamic variables
     ! Single-precision binary file input variables
     if (allocated(A_phi_r4)) deallocate(A_phi_r4,stat=rcd)
     if (allocated(A_psi_r4)) deallocate(A_psi_r4,stat=rcd)
     if (allocated(B_phi_r4)) deallocate(B_phi_r4,stat=rcd)
     if (allocated(B_psi_r4)) deallocate(B_psi_r4,stat=rcd)
     if (allocated(S_d_r4)) deallocate(S_d_r4,stat=rcd)
     if (allocated(S_p_r4)) deallocate(S_p_r4,stat=rcd)
     if (allocated(alpha_r4)) deallocate(alpha_r4,stat=rcd)
     if (allocated(wvn_ctr_r4)) deallocate(wvn_ctr_r4,stat=rcd)
     
     ! Compute wavenumber coordinates
     bnd_dlt=wvn_ctr(2)-wvn_ctr(1)
     float_foo=bnd_dlt/2.0
     do bnd_idx=1,bnd_nbr
        wvn_min(bnd_idx)=wvn_ctr(bnd_idx)-float_foo
     enddo                  ! end loop over bnd
     do bnd_idx=1,bnd_nbr-1
        wvn_max(bnd_idx)=wvn_min(bnd_idx+1)
     enddo                  ! end loop over bnd
     wvn_max(bnd_nbr)=wvn_ctr(bnd_nbr)+float_foo
     
  endif                     ! endif input file is in BPB binary format
  
  ! Set S_p to 1.0 when there are no lines in band (i.e., when S_d is 0.0) 
  ! This ensures S_p from my HITRAN calculations equals S_p from BPB's
  do bnd_idx=1,bnd_nbr
     if (S_d(bnd_idx) == 0.0) then
        S_p(bnd_idx)=1.0
     endif
  enddo                  ! end loop over bnd
  
  if (dbg_lvl==dbg_io) then
     bnd_idx=1
     write (6,100) &
          'band = ',bnd_idx, &
          'wvn_ctr = ',wvn_ctr(bnd_idx),' units, ', &
          'S_d = ',S_d(bnd_idx),' units, ', &
          'S_p = ',S_p(bnd_idx),' units, ', &
          'A_phi = ',A_phi(bnd_idx),' units, ', &
          'B_phi = ',B_phi(bnd_idx),' units, ', &
          'A_psi = ',A_psi(bnd_idx),' units, ', &
          'B_psi = ',B_psi(bnd_idx),' units, ', &
          'alpha = ',alpha(bnd_idx),' units, '
100  format( &
          a7,i4,/, &
          a10,es8.1,a8,/, &
          a10,es8.1,a8,/, &
          a10,es8.1,a8,/, &
          a10,es8.1,a8,/, &
          a10,es8.1,a8,/, &
          a10,es8.1,a8,/, &
          a10,es8.1,a8,/, &
          a10,es8.1,a8,/)
  endif                     ! endif dbg
  
  if (dbg_lvl==dbg_io) then
     write (6,'(a4,1x,8(a8,1x))') &
          'bnd_idx', &
          'wvn_ctr', &
          'S_d', &
          'S_p', &
          'A_phi', &
          'B_phi', &
          'A_psi', &
          'B_psi', &
          'alpha'
     do bnd_idx=1,bnd_nbr
        write (6,'(i4,1x,8(es8.1,1x))') &
             bnd_idx, &
             wvn_ctr(bnd_idx), &
             S_d(bnd_idx), &
             S_p(bnd_idx), &
             A_phi(bnd_idx), &
             B_phi(bnd_idx), &
             A_psi(bnd_idx), &
             B_psi(bnd_idx), &
             alpha(bnd_idx)
     enddo                  ! end loop over bnd
  endif                     ! endif dbg
  
  do bnd_idx=1,bnd_nbr
     if (wvn_min(bnd_idx) == 0.0) then
        ! Prevent ln_lo from crashing program
        write(6,'(2a,es8.1,a,es8.1,a)') prg_nm(1:ftn_strlen(prg_nm)), &
             ': WARNING ln_lo is 0.0 cm-1 so setting wvl_max to ', &
             1.0/(100.0*ln_lo_min),' m, i.e., effective ln_lo is ', &
             ln_lo_min,' cm-1 to prevent divide by zero problems'
        wvl_max(bnd_idx)=1.0/(100.0*ln_lo_min)
     else   
        wvl_max(bnd_idx)=1.0/(100.0*wvn_min(bnd_idx))
     endif                  ! endif
     wvl_min(bnd_idx)=1.0/(100.0*wvn_max(bnd_idx))
     wvl_ctr(bnd_idx)=0.5*(wvl_min(bnd_idx)+wvl_max(bnd_idx))
     wvl_dlt(bnd_idx)=wvl_max(bnd_idx)-wvl_min(bnd_idx)
     wvl_grd(bnd_idx)=wvl_max(bnd_idx)
     wvn_dlt(bnd_idx)=wvn_max(bnd_idx)-wvn_min(bnd_idx)
     bnd(bnd_idx)=wvn_ctr(bnd_idx)
     frc_slr_flx_ThD71(bnd_idx)= & ! NB: arguments in wavenumbers
          slfwtd(wvn_max(bnd_idx),wvn_min(bnd_idx))
     frc_slr_flx_LaN68(bnd_idx)= & ! NB: arguments in wavenumbers
          slfwln(wvn_max(bnd_idx),wvn_min(bnd_idx))
  enddo                  ! end loop over bnd
  wvl_grd(bnd_nbr+1)=wvl_min(bnd_nbr)
  
  ! Index of refraction through dry air at STP (Len93 p. 155)
  do bnd_idx=1,bnd_nbr
     idx_rfr_air_STP(bnd_idx)= &
          1.0+ &
          1.0e-6*(77.46+0.459/(1.0e12*wvl_ctr(bnd_idx)**2))* &
          prs_STP*0.01/tpt_STP
  enddo                  ! end loop over bnd
  
  ! Create some useful SI versions of line strength parameters
  !  write(6,'(a)') 'Creating S_d_abs_cff_mss and S_p_abs_cff_mss:'
  do bnd_idx=1,bnd_nbr
     S_d_abs_cff_mss(bnd_idx)=S_d(bnd_idx)*Avagadro/mmw_mlc(mlc_id) ! [cm-1 m2 mlc-1] -> [cm-1 m2 kg-1]
     S_p_abs_cff_mss(bnd_idx)=S_p(bnd_idx)*Avagadro/mmw_mlc(mlc_id) ! [cm-1 m2 mlc-1] -> [cm-1 m2 kg-1]
  enddo                  ! end loop over bnd
  
  ! Compute phi and psi from least-squares fit parameters
  do tpt_idx=1,tpt_nbr
     phi_fit(tpt_idx)=exp(A_phi(bnd_dbg)*tpt_dlt_1(tpt_idx)+ &
          B_phi(bnd_dbg)*tpt_dlt_2(tpt_idx))
     psi_fit(tpt_idx)=exp(A_psi(bnd_dbg)*tpt_dlt_1(tpt_idx)+ &
          B_psi(bnd_dbg)*tpt_dlt_2(tpt_idx))
  enddo                  ! end loop over t
  
  do bnd_idx=1,bnd_nbr
     if (S_p(bnd_idx) == 0.0) then
        write(6,*) 'ERROR: S_p(',bnd_idx,') = 0.0 at wvn =', &
             wvn_ctr(bnd_idx),' cm-1, wvl = ',wvl_ctr(bnd_idx)* 1.0e9,' nm'
     endif
     !if (S_d(bnd_idx) == 0.0) then
     !        write(6,*) 'WARNING: S_d(',bnd_idx,') = 0.0 at wvn =', &
     !             wvn_ctr(bnd_idx),' cm-1, wvl = ',wvl_ctr(bnd_idx)* 1.0e9,' nm'
     !endif
  enddo                  ! end loop over bnd

  ! Check that single precision bounds were not exceeded
  do bnd_idx=1,bnd_nbr
     if (S_d(bnd_idx) < 1.0e-36.and.S_d(bnd_idx) > 0.0) then 
        write(6,*) 'WARNING: Single precision underflow S_d(',bnd_idx,') = ',S_d(bnd_idx)
        S_d(bnd_idx)=0.0
     endif
     if (S_d_abs_cff_mss(bnd_idx) < 1.0e-36.and.S_d_abs_cff_mss(bnd_idx) > 0.0) then 
        write(6,*) 'WARNING: Single precision underflow S_d_abs_cff_mss(',bnd_idx,') = ',S_d_abs_cff_mss(bnd_idx)
        S_d_abs_cff_mss(bnd_idx)=0.0
     endif
     if (S_p_abs_cff_mss(bnd_idx) < 1.0e-36.and.S_p_abs_cff_mss(bnd_idx) > 0.0) then 
        write(6,*) 'WARNING: Single precision underflow S_p_abs_cff_mss(',bnd_idx,') = ',S_p_abs_cff_mss(bnd_idx)
        S_p(bnd_idx)=1.0
        S_p_abs_cff_mss(bnd_idx)=S_p(bnd_idx)*Avagadro/mmw_mlc(mlc_id) ! [m2 mlc-1] -> [m2 kg-1]
     endif
     if (S_d(bnd_idx) < 0.0.or.S_p(bnd_idx) < 0.0) then 
        write(6,'(/,a,i4,a4,f9.6,a2,f9.6,a6,f9.3,a2,f9.3,a5)')  &
             'WARNING: Impossible band parameters at bnd(',bnd_idx,') = ', &
             wvl_min(bnd_idx)*1.0e6,'--',wvl_max(bnd_idx)*1.0e6,' um = ', &
             wvn_min(bnd_idx),'--',wvn_max(bnd_idx),' cm-1'
        write (6,'(4(a,es10.3))')  &
             'S_d = ',S_d(bnd_idx),', S_d_abs_cff_mss = ',S_d_abs_cff_mss(bnd_idx), &
             'S_p = ',S_p(bnd_idx),', S_p_abs_cff_mss = ',S_p_abs_cff_mss(bnd_idx)
        stop
     endif                  ! endif
  end do                    ! end loop over bnd
  
#ifdef ENABLE_NETCDF4
  if (fl_out_fmt == nco_format_undefined) fl_out_fmt=nf90_format_classic ! [enm] Output file format
  if (fl_out_fmt == nf90_format_64bit) then
     nf90_create_mode=nf90_create_mode+nf90_64bit_offset
  else if (fl_out_fmt == nf90_format_netcdf4) then
     nf90_create_mode=nf90_create_mode+nf90_netcdf4
  else if (fl_out_fmt == nf90_format_netcdf4_classic) then
     nf90_create_mode=nf90_create_mode+(nf90_classic_model+nf90_netcdf4)
  end if ! end else fl_out_fmt
#else /* !ENABLE_NETCDF4 */
  if (fl_out_fmt == nco_format_undefined) fl_out_fmt=nf90_format_classic ! [enm] Output file format
  if(fl_out_fmt == nf90_format_classic) nf90_create_mode=nf90_create_mode+0 ! CEWI
#endif /* !ENABLE_NETCDF4 */
  dfl_lvl=dfl_lvl+0 ! CEWI
  flg_dfl=flg_dfl+0 ! CEWI
  flg_shf=flg_shf+0 ! CEWI
  rcd=nf90_wrp_create(fl_out,nf90_create_mode,nc_id,sbr_nm=sbr_nm)
  
  ! Define dimension IDs
  rcd=nf90_wrp(nf90_def_dim(nc_id,'bnd',bnd_nbr,bnd_dmn_id),sbr_nm//': def_dim bnd in '//__FILE__)
  rcd=nf90_wrp(nf90_def_dim(nc_id,'grd',bnd_nbr+1,grd_dmn_id),sbr_nm//': def_dim grd in '//__FILE__)
  rcd=nf90_wrp(nf90_def_dim(nc_id,'t',tpt_nbr,tpt_dmn_id),sbr_nm//': def_dim tpt in '//__FILE__)
  
  ! Variable definitions
  rcd=nf90_wrp(nf90_def_var(nc_id,'A_phi',typ_out,bnd_dmn_id,A_phi_id),sbr_nm//': dv A_phi')
  rcd=nf90_wrp(nf90_def_var(nc_id,'A_psi',typ_out,bnd_dmn_id,A_psi_id),sbr_nm//': dv A_psi')
  rcd=nf90_wrp(nf90_def_var(nc_id,'B_phi',typ_out,bnd_dmn_id,B_phi_id),sbr_nm//': dv B_phi')
  rcd=nf90_wrp(nf90_def_var(nc_id,'B_psi',typ_out,bnd_dmn_id,B_psi_id),sbr_nm//': dv B_psi')
  rcd=nf90_wrp(nf90_def_var(nc_id,'S_d',typ_out,bnd_dmn_id,S_d_id),sbr_nm//': dv S_d')
  rcd=nf90_wrp(nf90_def_var(nc_id,'S_d_abs_cff_mss',typ_out,bnd_dmn_id,S_d_abs_cff_mss_id),sbr_nm//': dv S_d_abs_cff_mss')
  rcd=nf90_wrp(nf90_def_var(nc_id,'S_p',typ_out,bnd_dmn_id,S_p_id),sbr_nm//': dv S_p')
  rcd=nf90_wrp(nf90_def_var(nc_id,'S_p_abs_cff_mss',typ_out,bnd_dmn_id,S_p_abs_cff_mss_id),sbr_nm//': dv S_p_abs_cff_mss')
  rcd=nf90_wrp(nf90_def_var(nc_id,'alpha',typ_out,bnd_dmn_id,alpha_id),sbr_nm//': dv alpha')
  rcd=nf90_wrp(nf90_def_var(nc_id,'bnd',typ_out,bnd_dmn_id,bnd_id),sbr_nm//': dv bnd')
  rcd=nf90_wrp(nf90_def_var(nc_id,'bnd_ln_nbr',nf90_int,bnd_dmn_id,bnd_ln_nbr_id),sbr_nm//': dv bnd_ln_nbr')
  rcd=nf90_wrp(nf90_def_var(nc_id,'bnd_dlt',typ_out,bnd_dlt_id),sbr_nm//': dv bnd_dlt')
  rcd=nf90_wrp(nf90_def_var(nc_id,'frc_slr_flx_LaN68',typ_out,bnd_dmn_id,frc_slr_flx_LaN68_id),sbr_nm//': dv frc_slr_flx_LaN68')
  rcd=nf90_wrp(nf90_def_var(nc_id,'frc_slr_flx_ThD71',typ_out,bnd_dmn_id,frc_slr_flx_ThD71_id),sbr_nm//': dv frc_slr_flx_ThD71')
  rcd=nf90_wrp(nf90_def_var(nc_id,'idx_rfr_air_STP',typ_out,bnd_dmn_id,idx_rfr_air_STP_id),sbr_nm//': dv idx_rfr_air_STP')
  rcd=nf90_wrp(nf90_def_var(nc_id,'oneD_foo',typ_out,tpt_dmn_id,oneD_foo_id),sbr_nm//': dv oneD_foo')
  rcd=nf90_wrp(nf90_def_var(nc_id,'phi_exact',typ_out,tpt_dmn_id,phi_exact_id),sbr_nm//': dv phi_exact')
  rcd=nf90_wrp(nf90_def_var(nc_id,'phi_fit',typ_out,tpt_dmn_id,phi_fit_id),sbr_nm//': dv phi_fit')
  rcd=nf90_wrp(nf90_def_var(nc_id,'psi_exact',typ_out,tpt_dmn_id,psi_exact_id),sbr_nm//': dv psi_exact')
  rcd=nf90_wrp(nf90_def_var(nc_id,'psi_fit',typ_out,tpt_dmn_id,psi_fit_id),sbr_nm//': dv psi_fit')
  rcd=nf90_wrp(nf90_def_var(nc_id,'scalar_foo',typ_out,scalar_foo_id),sbr_nm//': dv scalar_foo')
  rcd=nf90_wrp(nf90_def_var(nc_id,'tpt',typ_out,tpt_dmn_id,tpt_id),sbr_nm//': dv tpt')
  rcd=nf90_wrp(nf90_def_var(nc_id,'wvl_ctr',typ_out,bnd_dmn_id,wvl_ctr_id),sbr_nm//': dv wvl_ctr')
  rcd=nf90_wrp(nf90_def_var(nc_id,'wvl_grd',typ_out,grd_dmn_id,wvl_grd_id),sbr_nm//': dv wvl_grd')
  rcd=nf90_wrp(nf90_def_var(nc_id,'wvl_max',typ_out,bnd_dmn_id,wvl_max_id),sbr_nm//': dv wvl_max')
  rcd=nf90_wrp(nf90_def_var(nc_id,'wvl_min',typ_out,bnd_dmn_id,wvl_min_id),sbr_nm//': dv wvl_min')
  rcd=nf90_wrp(nf90_def_var(nc_id,'wvl_dlt',typ_out,bnd_dmn_id,wvl_dlt_id),sbr_nm//': dv wvl_dlt')
  rcd=nf90_wrp(nf90_def_var(nc_id,'wvn_ctr',typ_out,bnd_dmn_id,wvn_ctr_id),sbr_nm//': dv wvn_ctr')
  rcd=nf90_wrp(nf90_def_var(nc_id,'wvn_grd',typ_out,grd_dmn_id,wvn_grd_id),sbr_nm//': dv wvn_grd')
  rcd=nf90_wrp(nf90_def_var(nc_id,'wvn_max',typ_out,bnd_dmn_id,wvn_max_id),sbr_nm//': dv wvn_max')
  rcd=nf90_wrp(nf90_def_var(nc_id,'wvn_min',typ_out,bnd_dmn_id,wvn_min_id),sbr_nm//': dv wvn_min')
  rcd=nf90_wrp(nf90_def_var(nc_id,'wvn_dlt',typ_out,bnd_dmn_id,wvn_dlt_id),sbr_nm//': dv wvn_dlt')
  rcd=nf90_wrp(nf90_def_var(nc_id,'iso_id',nf90_int,iso_id_id),sbr_nm//': dv iso_id')
  rcd=nf90_wrp(nf90_def_var(nc_id,'mlc_id',nf90_int,mlc_id_id),sbr_nm//': dv mlc_id')
  
  ! Add global attributes
  rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'CVS_Id',CVS_Id),sbr_nm//': pa CVS_Id')
  rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'creation_date',lcl_date_time),sbr_nm//': pa creation_date')
  rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'prg_ID',prg_ID(1:ftn_strlen(prg_ID))),sbr_nm//': pa prg_ID')
  rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'cmd_ln',cmd_ln(1:ftn_strlen(cmd_ln))),sbr_nm//': pa cmd_ln')
  rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'molecule',mlc_sng(1:ftn_strlen(mlc_sng))),sbr_nm//': pa molecule')
  rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'isotope',iso_sng(1:ftn_strlen(iso_sng))),sbr_nm//': pa isotope')
  
  ! Add English text descriptions
  rcd=nf90_wrp(nf90_put_att(nc_id,A_phi_id,'long_name','Linear temperature dependence of line strengths'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,A_psi_id,'long_name','Linear temperature dependence of Lorentzian HWHM'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,B_phi_id,'long_name','Quadratic temperature dependence of line strengths'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,B_psi_id,'long_name','Quadratic temperature dependence of Lorentzian HWHM'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,S_d_abs_cff_mss_id,'long_name','Band average line strength'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,S_d_id,'long_name','Band average line strength'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,S_p_abs_cff_mss_id,'long_name','Band average line strength amount over line width'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,S_p_id,'long_name','Band average line strength amount over line width'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,alpha_id,'long_name','Description'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,bnd_id,'long_name','Band center'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,bnd_ln_nbr_id,'long_name','# of HITRAN lines in each band'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,bnd_dlt_id,'long_name','Uniform width of bands'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,frc_slr_flx_LaN68_id,'long_name','Fraction of solar flux: Labs & Neckel 1968'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,frc_slr_flx_ThD71_id,'long_name','Fraction of solar fluxL Thekeakara & Drummond 1971'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,idx_rfr_air_STP_id,'long_name','Index of refraction at band center at STP'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,iso_id_id,'long_name','HITRAN isotope number (1..9)'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,mlc_id_id,'long_name','HITRAN molecule number (1..37)'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,oneD_foo_id,'long_name','Description'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,phi_exact_id,'long_name','Phi exactly computed'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,phi_fit_id,'long_name','Phi from least squares fit'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,psi_exact_id,'long_name','Psi exactly computed'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,psi_fit_id,'long_name','Psi from least squares fit'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,scalar_foo_id,'long_name','Description'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,tpt_id,'long_name','Temperature'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,wvl_ctr_id,'long_name','Band center wavelength'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,wvl_grd_id,'long_name','Wavelength grid'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,wvl_max_id,'long_name','Band maximum wavelength'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,wvl_min_id,'long_name','Band minimum wavelength'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,wvl_dlt_id,'long_name','Bandwidth'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,wvn_ctr_id,'long_name','Band center wavenumber'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,wvn_grd_id,'long_name','Wavenumber grid'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,wvn_max_id,'long_name','Band maximum wavenumber'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,wvn_min_id,'long_name','Band minimum wavenumber'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,wvn_dlt_id,'long_name','Bandwidth'),sbr_nm)
  
  ! Add units
  rcd=nf90_wrp(nf90_put_att(nc_id,A_phi_id,'units','kelvin-1'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,A_psi_id,'units','kelvin-1'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,B_phi_id,'units','kelvin-2'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,B_psi_id,'units','kelvin-2'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,S_d_abs_cff_mss_id,'units','centimeter-1 meter2 kilogram-1'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,S_d_id,'units','meter2 molecule-1'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,S_p_abs_cff_mss_id,'units','centimeter-1 meter2 kilogram-1'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,S_p_id,'units','meter2 molecule-1'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,alpha_id,'units','unknown'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,bnd_id,'units','centimeter-1'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,bnd_ln_nbr_id,'units','cardinal'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,bnd_dlt_id,'units','centimeter-1'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,frc_slr_flx_LaN68_id,'units','fraction'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,frc_slr_flx_ThD71_id,'units','fraction'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,idx_rfr_air_STP_id,'units','fraction'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,oneD_foo_id,'units','unknown'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,phi_exact_id,'units','fraction'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,phi_fit_id,'units','fraction'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,psi_exact_id,'units','fraction'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,psi_fit_id,'units','fraction'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,scalar_foo_id,'units','unknown'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,tpt_id,'units','kelvin'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,wvl_ctr_id,'units','meter'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,wvl_grd_id,'units','meter'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,wvl_max_id,'units','meter'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,wvl_min_id,'units','meter'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,wvl_dlt_id,'units','meter'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,wvn_ctr_id,'units','centimeter-1'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,wvn_grd_id,'units','centimeter-1'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,wvn_max_id,'units','centimeter-1'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,wvn_min_id,'units','centimeter-1'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,wvn_dlt_id,'units','centimeter-1'),sbr_nm)
  
  ! All dimensions, variables, and attributes have been defined
  rcd=nf90_wrp(nf90_enddef(nc_id),sbr_nm//': enddef in '//__FILE__) 
  
  ! Write data
  rcd=nf90_wrp(nf90_put_var(nc_id,bnd_dlt_id,bnd_dlt),sbr_nm//': pv bnd_dlt in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,iso_id_id,iso_id),sbr_nm//': pv iso in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,mlc_id_id,mlc_id),sbr_nm//': pv mlc in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,scalar_foo_id,scalar_foo),sbr_nm//': pv scalar_foo in '//__FILE__)
  
  rcd=nf90_wrp(nf90_put_var(nc_id,A_phi_id,A_phi),sbr_nm//': pv A_phi in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,A_psi_id,A_psi),sbr_nm//': pv A_psi in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,B_phi_id,B_phi),sbr_nm//': pv B_phi in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,B_psi_id,B_psi),sbr_nm//': pv B_psi in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,S_d_abs_cff_mss_id,S_d_abs_cff_mss),sbr_nm//': pv S_d_abs_cff_mss in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,S_d_id,S_d),sbr_nm//': pv S_d in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,S_p_abs_cff_mss_id,S_p_abs_cff_mss),sbr_nm//': pv S_p_abs_cff_mss in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,S_p_id,S_p),sbr_nm//': pv S_p in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,alpha_id,alpha),sbr_nm//': pv alpha in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,bnd_id,bnd),sbr_nm//': pv bnd in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,bnd_ln_nbr_id,bnd_ln_nbr),sbr_nm//': pv bnd_ln_nbr in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,frc_slr_flx_LaN68_id,frc_slr_flx_LaN68),sbr_nm//': pv frc_slr_flx_LaN68 in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,frc_slr_flx_ThD71_id,frc_slr_flx_ThD71),sbr_nm//': pv frc_slr_flx_ThD71 in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,idx_rfr_air_STP_id,idx_rfr_air_STP),sbr_nm//': pv idx_rfr_air_STP in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,oneD_foo_id,oneD_foo),sbr_nm//': pv oneD_foo in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,phi_exact_id,phi_exact),sbr_nm//': pv phi_exact in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,phi_fit_id,phi_fit),sbr_nm//': pv phi_fit in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,psi_exact_id,psi_exact),sbr_nm//': pv psi_exact in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,psi_fit_id,psi_fit),sbr_nm//': pv psi_fit in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,tpt_id,tpt),sbr_nm//': pv tpt in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,wvl_ctr_id,wvl_ctr),sbr_nm//': pv wvl_ctr in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,wvl_grd_id,wvl_grd),sbr_nm//': pv wvl_grd in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,wvl_max_id,wvl_max),sbr_nm//': pv wvl_max in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,wvl_min_id,wvl_min),sbr_nm//': pv wvl_min in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,wvl_dlt_id,wvl_dlt),sbr_nm//': pv wvl_dlt in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,wvn_ctr_id,wvn_ctr),sbr_nm//': pv wvn_ctr in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,wvn_grd_id,wvn_grd),sbr_nm//': pv wvn_grd in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,wvn_max_id,wvn_max),sbr_nm//': pv wvn_max in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,wvn_min_id,wvn_min),sbr_nm//': pv wvn_min in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,wvn_dlt_id,wvn_dlt),sbr_nm//': pv wvn_dlt in '//__FILE__)

  rcd=nf90_wrp_close(nc_id,fl_out,'Wrote results to') ! [fnc] Close file
  
  ! De-allocate dynamic variables
  ! Local variables
  if (allocated(bnd_ln_nbr)) deallocate(bnd_ln_nbr,stat=rcd)
  if (allocated(ln_idx_max)) deallocate(ln_idx_max,stat=rcd)
  if (allocated(ln_idx_min)) deallocate(ln_idx_min,stat=rcd)
  if (allocated(Boltzmann_wgt_tpt_rcp)) deallocate(Boltzmann_wgt_tpt_rcp,stat=rcd)
  if (allocated(bnd_sum_str_ln_lmt_rfr_rcp_2)) deallocate(bnd_sum_str_ln_lmt_rfr_rcp_2,stat=rcd)
  if (allocated(bnd_sum_wk_ln_lmt_rfr_rcp)) deallocate(bnd_sum_wk_ln_lmt_rfr_rcp,stat=rcd)
  if (allocated(prt_fnc_tpt_scl)) deallocate(prt_fnc_tpt_scl,stat=rcd)
  if (allocated(tpt_HITRAN_tpt_rcp)) deallocate(tpt_HITRAN_tpt_rcp,stat=rcd)
  if (allocated(tpt_dlt_1)) deallocate(tpt_dlt_1,stat=rcd)
  if (allocated(tpt_dlt_2)) deallocate(tpt_dlt_2,stat=rcd)
  if (allocated(tpt_dlt_3)) deallocate(tpt_dlt_3,stat=rcd)
  if (allocated(tpt_dlt_4)) deallocate(tpt_dlt_4,stat=rcd)
  if (allocated(tpt_rcp)) deallocate(tpt_rcp,stat=rcd)
  ! Computational-precision netCDF output variables
  if (allocated(A_phi)) deallocate(A_phi,stat=rcd)
  if (allocated(A_psi)) deallocate(A_psi,stat=rcd)
  if (allocated(B_phi)) deallocate(B_phi,stat=rcd)
  if (allocated(B_psi)) deallocate(B_psi,stat=rcd)
  if (allocated(S_d)) deallocate(S_d,stat=rcd)
  if (allocated(S_d_abs_cff_mss)) deallocate(S_d_abs_cff_mss,stat=rcd)
  if (allocated(S_p)) deallocate(S_p,stat=rcd)
  if (allocated(S_p_abs_cff_mss)) deallocate(S_p_abs_cff_mss,stat=rcd)
  if (allocated(alpha)) deallocate(alpha,stat=rcd)
  if (allocated(bnd)) deallocate(bnd,stat=rcd) ! coordinate variable
  if (allocated(frc_slr_flx_LaN68)) deallocate(frc_slr_flx_LaN68,stat=rcd)
  if (allocated(frc_slr_flx_ThD71)) deallocate(frc_slr_flx_ThD71,stat=rcd)
  if (allocated(idx_rfr_air_STP)) deallocate(idx_rfr_air_STP,stat=rcd)
  if (allocated(oneD_foo)) deallocate(oneD_foo,stat=rcd)
  if (allocated(phi_exact)) deallocate(phi_exact,stat=rcd)
  if (allocated(phi_fit)) deallocate(phi_fit,stat=rcd)
  if (allocated(psi_exact)) deallocate(psi_exact,stat=rcd)
  if (allocated(psi_fit)) deallocate(psi_fit,stat=rcd)
  if (allocated(tpt)) deallocate(tpt,stat=rcd)
  if (allocated(wvl_ctr)) deallocate(wvl_ctr,stat=rcd)
  if (allocated(wvl_dlt)) deallocate(wvl_dlt,stat=rcd)
  if (allocated(wvl_grd)) deallocate(wvl_grd,stat=rcd)
  if (allocated(wvl_max)) deallocate(wvl_max,stat=rcd)
  if (allocated(wvl_min)) deallocate(wvl_min,stat=rcd)
  if (allocated(wvn_ctr)) deallocate(wvn_ctr,stat=rcd)
  if (allocated(wvn_dlt)) deallocate(wvn_dlt,stat=rcd)
  if (allocated(wvn_grd)) deallocate(wvn_grd,stat=rcd)
  if (allocated(wvn_max)) deallocate(wvn_max,stat=rcd)
  if (allocated(wvn_min)) deallocate(wvn_min,stat=rcd)
  
1000 continue
  
  rcd=0
  call exit(rcd)
end program htrn2nb
