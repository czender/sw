! $Id$

program O2X
  
  ! Purpose: Convert O2-O2 absorption cross section data to netCDF format
  ! NB: O2X should be compiled/run in double precision
  ! 20181209: Until today, abs_xsx_O2O2 was multiplied by k_O2_O2 to allow single precision cross-sections
  ! This clever and obscure manipulation allows entire swnb2 suite to run in single precision
  ! It is equivalent to archiving absorption cross sections of O4, an equilibrium product of O2+O2
  ! This was important in 1990s with g77 compiler on small machines
  ! As of 2018, we always run in double precision anyway because some IR gases (e.g., CO) have weak cross-sections
  
  ! CCMCG stands for CCM_CGS, i.e., quantity is in CGS units not SI
  ! Abbreviation is required to keep symbol names <= 31 characters in some Fortran compilers
  
  ! Compilation:
  ! cd ${HOME}/sw/aca; make -W O2X.F OPTS=D O2X; cd -
  ! cd ${HOME}/sw/aca; make -W O2X.F O2X; cd -
  ! cd ${HOME}/sw/aca; make OPTS=D O2X; cd -
  
  ! O2X --GOB90 ! Greenblatt et al. (1990) data
  ! O2X --HTR16 --drc_out=${DATA}/aca ! HITRAN data
  ! O2X -i ${DATA}/aca/abs_xsx_O2O2.txt -o ${DATA}/aca/abs_xsx_O2O2_GOB90.nc
  ! Exclude 1.26 micron band
  ! O2X -x -i ${DATA}/aca/abs_xsx_O2O2.txt -o ${DATA}/aca/abs_xsx_O2O2_xcl_1260nm.nc

  ! Usage:
  ! ncks -H -C -F -v wvl_ctr_CCM,flx_slr_frc_CCM,abs_xsx_O2O2_CCM ${DATA}/aca/abs_xsx_O2O2.nc
  ! ncks -H -C -F -d bnd,1.26e-6 -d bnd_CCM,1.26 -v wvl_ctr_CCM,flx_slr_frc_CCM,abs_xsx_O2O2_CCM,abs_xsx_O2O2 ${DATA}/aca/abs_xsx_O2O2.nc
  ! ncwa -a bnd -d bnd,0.295e-6,0.305e-6 ${DATA}/aca/abs_xsx_O2O2.nc ${DATA}/aca/O2O2.nc
  ! ncks -H -C -F -d bnd,0.3e-6 -v abs_xsx_O2O2 ${DATA}/aca/abs_xsx_O2O2.nc
  ! ncks -H -C -F -d bnd,1.26e-6 -v odal_O2O2_PUMC_O2_PUMP_O2,odal_O2N2_PUMC_O2_PUMP_N2 ${DATA}/aca/abs_xsx_O2O2.nc
  ! ncks -H -C -F -d bnd,1.26e-6 -v odal_O2O2_PUMC_O2_PUMP_O2,odal_O2O2_PUNC_O2_PUNP_O2,odal_O2O2_PUMC_O2_PUMP_O2_CCM ${DATA}/aca/abs_xsx_O2O2.nc
  ! ncks -H -C -F -v abs_xsx_O2O2 ${DATA}/aca/abs_xsx_O2O2.nc
  ! ncks -H -C -F -d bnd_CCM,6 -v abs_xsx_O2O2_CCM ${DATA}/aca/abs_xsx_O2O2.nc
  
  ! NB: data were measured and are recorded as absorption cross section excess
  ! (relative to the absorption cross section of O2) per unit concentration of O2.
  ! Thus typical numbers are O(1.0e-46), and must be handled in double precision
  
  ! Currently the code processes input ASCII data files that look like:
  
  !#Greenblatt, et al. O2-O2 absorption (ref: JGR, 95, 18577, 1990)
  !#3182 points (note: 335.2-663.4 nm & 1005.0-1136.8 nm)
  !#Column 1: wavelength in nm
  !#column 2: "O4" absorption in cm^5 mol^-2
  !335.2   2.079e-47
  !335.3   1.975e-47
  !...
  
  ! or, as of 201812, processes input .cia ASCII data files from HITRAN that look like:
  ! Vertical bar placed on last character of each field
  ! 1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
  ! O2-O2_2018.cia:
  ! First format line, from RGR12, does not fit 2018 data:
  ! ...mlc_frm_htrn....| wvn_min | wvn_max |bndnbr| tpt  | prs  | xsx_max | rsn|...cmt_htrn.........|rf|
  ! Second format line is from https://hitran.org/data/CIA/CIA_Readme.pdf and fits 2018 data. No prs_htrn, rsn_nst_htrn is length 6, cmt_htrn is length 27:
  ! ...mlc_frm_htrn....| wvn_min | wvn_max |bndnbr| tpt  | xsx_max | rsn |...cmt_htrn...............|rf|
  !               O2-O2  1150.000  1950.000   4002  193.4 8.767E-45 0.500               P up to 4atm  2
  ! O2-N2_2018.cia:
  ! ...mlc_frm_htrn....| wvn_min | wvn_max |bndnbr| tpt  | xsx_max | rsn |...cmt_htrn...............|rf|
  !               O2-N2  7450.1322 8487.4645   4237 293.00 9.349e-46 -.999           Fit to Mate 1999 10

  !  1150.0000 -6.497E-48
  !  1150.2000 -5.749E-48
  !  1150.4000 -7.841E-48
  ! ...
 
#if 0
  From http://hitran.org/data/CIA/CIA_Readme.pdf
  "More complete details of the CIA data sets are presented in C. Richard, I.E. Gordon, L.S. Rothman, M. Abel, L. Frommhold, M. Gustafsson, et al, JQSRT 113, 1276-1285 (2012).
  The data sets in the individual files are ordered by wavenumber spectral intervals, and secondly by increasing temperature."

  CSZ modified following text from .xsc description to apply to .cia files:
  "In the HITRAN FTP site, the data are presented as separate files for each individual molecule. Each portion of the file corresponding to a particular temperature-pressure pair begins with a header (see Table 1) that contains information on the wavenumber (cm−1) range, number of cross-section data in this set, temperature (K), and pressure (Torr). The maximum value of the absorption cross sections (fxm) and additional information containing the reference to that observation are also presented in each header. The cross sections have been cast into an equal wavenumber interval grid. It should be noted that the initial and final wavenumbers, νmin and νmax, respectively, of each temperature–pressure set for a given wavenumber region are not always identical. They have been taken from the analysis of the observations. The sampling intervals are also not necessarily identical for each temperature–pressure set. The wavenumber interval of the grid is obtained by taking the difference of the initial and final wavenumber and dividing this quantity by the number of points, N, minus one, i.e., Δν=(νmax−νmin)/(N−1).

  This value of N is provided so that a user’s personal program can read the table of cross- sections that follows the header. Note that the use of the features of HITRANonline makes much of this discussion transparent.

  The table below illustrates the format of each header record. 

  Quantity	Field length	Type	Comment
  Molecule	20	Character	Chemical formula (right-justified)
  Minimum wavenumber, νmin	10	Real	Start of range (cm−1) (-0.999 is missing value)
  Maximum wavenumber, νmax	10	Real	End of range (cm−1) (-0.999 is missing value)
  Number of points, N	7	Integer	Number of cross-sections in set
  Temperature, T	7	Real	Temperature (K) of set
  Pressure, P	        7	Real	Pressure of set in Torr (NB: 7 digits in .cia files not 6 digits as in .xsc files)
  Maximum cross-section value in set, σmax	10	Real	Useful for scaling plots (cm2/molecule fxm)
  Instrument resolution	5	Real	See note
  Comments	21	Character	Right-adjusted
  Reference	3	Integer	Index pointing to source of data
  Note: Most cross sections have been taken from Fourier transform spectrometer (FTS) measurements. In that case the resolution is given in cm−1. There are some cross-sections taken from grating spectrometer measurements in the UV. In those cases, the resolution is given in milli-Ångströms in the form xxx mÅ, where xxx are up to three digits.

  It is to be noted that the files may have many temperature–pressure sets for different spectral regions, as indicated by headers throughout the file. While the temperature–pressure (T,p) sets are reasonably complete for many species for an adequate simulation of atmospheric transmission in the spectral regions where those species are active, for other species an insufficiency of the (T,p) sets may become apparent. It is hoped that future measurements at extended sets of (T,p) combinations may help broaden the coverage in the database."

20181029: Many CIA cross-sections are negative!
#endif /* !0 */

  use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
  use netcdf ! [mdl] netCDF interface
  use nf90_utl ! [mdl] netCDF utilities
  use phys_cst_mdl ! [mdl] Fundamental and derived physical constants
  use prn_mdl,only:dat_prn_f77,dat_prn_f90 ! [mdl] Print formatting & utilities
  use rt_mdl ! [mdl] Radiative transfer utilities
  use sng_mdl ! [mdl] String manipulation
  use utl_mdl ! [mdl] Utility functions (date_time_get,mnt_chk...)
  use vec_mdl ! [mdl] Vector manipulation, interpolation, rebinning
  use wvl_mdl ! [mdl] Wavelength grid parameters
  use xtr_mdl ! [mdl] Extrapolation/interpolation handling
  
  implicit none
  ! Parameters
  character(*),parameter::CVS_Date='$Date$' ! [sng] Date string
  character(*),parameter::CVS_Revision='$Revision$' ! [sng] File revision string
  character(len=*),parameter::CVS_Id='$Id$' ! [sng] CVS Identification
  character(len=*),parameter::sbr_nm='O2X' ! [sng] Subroutine name
  character(*),parameter::fl_in_GOB90='abs_xsx_O2O2.txt'
  character(*),parameter::fl_in_HTR16_cold='O2-O2_2018.cia'
  character(*),parameter::fl_out_GOB90='abs_xsx_O2O2_GOB90.nc'
  character(*),parameter::fl_out_HTR16='abs_xsx_O2O2_HTR16.nc'
  character(*),parameter::fl_slr_dfl='spc_Kur95_01wvn.nc'
  character(*),parameter::fl_wgt_dfl='/data/zender/tmp/swnb_trp.nc'
  character(*),parameter::wgt_nm_dfl='flx_spc_dwn_sfc'
  character(*),parameter::nlc=char(0) ! [sng] NUL character = ASCII 0 = char(0)
  
  integer,parameter::bnd_nbr_GOB90=4086
  integer,parameter::bnd_nbr_HTR16=5817 ! [nbr] 5818 is # of interfaces, 5817 is # of bands
  integer,parameter::fl_in_unit=73
  integer,parameter::sng_lng_dfl_fl=80 ! [nbr] Default filename string length
  integer,parameter::sng_lng_dfl_stt=200 ! [nbr] Default statement string length

  real,parameter::mlt_fct=1.0e20 ! [frc]
  real,parameter::mss_val=nf90_fill_float ! Missing value = missing_value and/or _FillValue
  real,parameter::tpt_cold_GOB90=296.0
  real,parameter::tpt_warm_GOB90=296.0
  real,parameter::tpt_cold_HTR16=193.4 ! [K] 
  real,parameter::tpt_warm_HTR16=193.4 ! [K] 

  ! Locals with simple initialization and no command-line override
  ! integer::exit_status=0 ! [enm] Program exit status (non-standard Fortran)
  integer::rcd=nf90_noerr ! [rcd] Return success code

  ! Input Arguments
  ! Input/Output Arguments
  ! Output Arguments

  ! Command-line parsing
  character(2)::dsh_key ! [sng] command-line dash and switch
  character(200)::cmd_ln ! [sng] command-line
  character(200)::prg_ID ! [sng] Program ID
  character(26)::lcl_date_time ! [sng] Time formatted as Day Mth DD HH:MM:SS TZ YYYY
  character(80)::arg_val ! [sng] command-line argument value
  character(80)::opt_sng=nlc ! [sng] Option string
  integer::arg_idx ! [idx] Counting index
  integer::arg_nbr ! [nbr] Number of command-line arguments
  integer::opt_lng ! [nbr] Length of option

  ! Set defaults for command-line options 
  character(sng_lng_dfl_fl)::drc_in='/data/zender/aca'//nlc       ! [sng] Input directory
  character(sng_lng_dfl_fl)::drc_out=nlc ! [sng] Output directory
  character(80)::fl_in='in.nc'//nlc ! [sng] Input file
  character(80)::fl_out='foo.nc'//nlc ! [sng] Output file
  character(80)::fl_slr=fl_slr_dfl//nlc ! [sng] Solar spectra file
  character(80)::fl_wgt=fl_wgt_dfl//nlc ! [sng] Weight file
  character(80)::wgt_nm=wgt_nm_dfl//nlc ! [sng] Solar spectra file

  logical::flg_GOB90=.true.
  logical::flg_HTR16=.false.
  logical::WGT_TRN=.false.
  logical::cmd_ln_fl_in=.false.
  logical::cmd_ln_fl_out=.false.
  logical::flg_ncl_1260nm_bnd=.true.
  logical::flg_ncl_1530nm_bnd=.false.
  logical::std_tpt=.true.

  real::N2_cll_prt_fsh=0.2 ! [frc] Efficiency of N2 (relative to O2) as a collision partner for O2 in the 1.26 micron band
  real::tpt_cold=tpt_cold_GOB90
  real::tpt_warm=tpt_warm_GOB90
  real::tpt_std=296.0 ! Temperature at which generic O2O2 cross sections will be archived
  
  ! Local workspace
  character(sng_lng_dfl_fl)::src_fl_sng
  character(sng_lng_dfl_fl)::src_rfr_sng
  character(300)::src_wgt_sng
  character(sng_lng_dfl_fl)::lbl
  
  ! HITRAN CIA format
  real,parameter::mss_val_rsn_htrn=-0.999 ! [cm-1] Missing value for resolution (used when cross-sections are theoretical not measured)
  character(27)::cmt_htrn ! [sng] Comments (such as pressure range of data)
  character(20)::mlc_frm_htrn ! [sng] Molecule chemical formula (right-justified)
  integer::bnd_nbr_htrn ! [nbr] Number of wavenumber bins
  integer::rfr_nbr_htrn ! [idx] Reference number (# of bibliographic reference in HITRAN GRH17 reference list)
  real::rsn_nst_htrn ! [cm-1] Instrument resolution (-0.999 is missing value)
  real::tpt_htrn ! [K] Temperature
  real::wvn_max_htrn ! [cm-1] Wavenumber at end of range
  real::wvn_min_htrn ! [cm-1] Wavenumber at start of range
  real::xsx_max_htrn ! [cm2 mlc-1] Maximum cross-section
  real::wvn_rsn ! [cm-1] Wavenumber resolution

  integer::bnd_dmn_id        ! dimension ID for bands
  integer::grd_dmn_id        ! dimension ID for grid
  integer::bnd_idx           ! counting index
  integer::nc_id             ! file handle
  integer::bnd_nbr           ! dimension size
  integer::tpt_cold_id
  integer::tpt_warm_id
  integer::tpt_std_id
  
  integer::abs_cff_mss_O2O2_id
  integer::abs_xsx_O2N2_id
  integer::abs_xsx_O2O2_dadT_id
  integer::abs_xsx_O2O2_id
  integer::abs_xsx_O2O2_tpt_rfr_id
  integer::bnd_id            ! coordinate ID
  integer::flx_bnd_dwn_TOA_id
  integer::flx_bnd_pht_dwn_TOA_id
  integer::flx_slr_frc_id
  integer::flx_spc_dwn_TOA_id
  integer::idx_rfr_air_STP_id
  integer::nrg_pht_id
  integer::odal_O2N2_PUMC_O2_PUMP_N2_id
  integer::odal_O2N2_PUNC_O2_PUNP_N2_id
  integer::odal_O2O2_PUMC_O2_PUMP_O2_id
  integer::odal_O2O2_PUNC_O2_PUNP_O2_id
  integer::wvl_ctr_id
  integer::wvl_grd_id
  integer::wvl_max_id
  integer::wvl_min_id
  integer::wvl_dlt_id
  integer::wvn_ctr_id
  integer::wvn_grd_id
  integer::wvn_max_id
  integer::wvn_min_id
  integer::wvn_dlt_id
  
  real(selected_real_kind(p=12))::double_foo
  
  ! netCDF4 
  integer::dfl_lvl=0 ! [enm] Deflate level
  integer::flg_shf=1 ! [flg] Turn on netCDF4 shuffle filter
  integer::flg_dfl=1 ! [flg] Turn on netCDF4 deflate filter
  integer::fl_out_fmt=nco_format_undefined ! [enm] Output file format
  integer::nf90_create_mode=nf90_clobber ! [enm] Mode flag for nf90_create() call

  ! Allocatable variables
  real(selected_real_kind(p=12)),dimension(:),allocatable::odal_O2O2_PUNC_O2_PUNP_O2 ! [m5 mlc-2]
  real(selected_real_kind(p=12)),dimension(:),allocatable::odal_O2N2_PUNC_O2_PUNP_N2 ! [m5 mlc-2]
  real,dimension(:),allocatable::abs_cff_mss_O2O2
  real,dimension(:),allocatable::abs_xsx_O2N2 ! [m5 mlc-2]
  real(selected_real_kind(p=12)),dimension(:),allocatable::abs_xsx_O2O2
  real,dimension(:),allocatable::abs_xsx_O2O2_dadT ! [m5 mlc-2 K-1]
  real,dimension(:),allocatable::abs_xsx_O2O2_tpt_rfr ! [K]
  real,dimension(:),allocatable::bnd     ! coordinate variable
  real,dimension(:),allocatable::flx_bnd_dwn_TOA
  real,dimension(:),allocatable::flx_bnd_pht_dwn_TOA
  real,dimension(:),allocatable::flx_slr_frc
  real,dimension(:),allocatable::flx_spc_dwn_TOA
  real,dimension(:),allocatable::idx_rfr_air_STP
  real,dimension(:),allocatable::nrg_pht
  real,dimension(:),allocatable::odal_O2N2_PUMC_O2_PUMP_N2
  real,dimension(:),allocatable::odal_O2O2_PUMC_O2_PUMP_O2
  real,dimension(:),allocatable::wgt_spc
  real,dimension(:),allocatable::wvl     ! coordinate variable
  real,dimension(:),allocatable::wvl_ctr
  real,dimension(:),allocatable::wvl_dlt
  real,dimension(:),allocatable::wvl_grd
  real,dimension(:),allocatable::wvl_max
  real,dimension(:),allocatable::wvl_min
  real,dimension(:),allocatable::wvn     ! coordinate variable
  real,dimension(:),allocatable::wvn_ctr
  real,dimension(:),allocatable::wvn_dlt
  real,dimension(:),allocatable::wvn_grd
  real,dimension(:),allocatable::wvn_max
  real,dimension(:),allocatable::wvn_min
  real,dimension(:),allocatable::xsx_wgt_flx
  
  integer::bnd_nbr_CCM
  integer::nbr_dat_per_ln
  integer::slr_spc_xtr_typ
  integer::xtr_typ_LHS
  integer::xtr_typ_RHS
  
  real(selected_real_kind(p=12)) odal_O2O2_PUNC_O2_PUNP_O2_CCM(bnd_nbr_CCM_SW_max)
  real(selected_real_kind(p=12)) odal_O2N2_PUNC_O2_PUNP_N2_CCM(bnd_nbr_CCM_SW_max)
  
  real abs_cff_mss_O2N2_CCM(bnd_nbr_CCM_SW_max)
  real abs_cff_mss_O2O2_CCM(bnd_nbr_CCM_SW_max)
  real abs_xsx_O2N2_CCM(bnd_nbr_CCM_SW_max)
  real abs_xsx_O2O2_CCM(bnd_nbr_CCM_SW_max)
  real bnd_CCM(bnd_nbr_CCM_SW_max)
  real flx_slr_frc_CCM(bnd_nbr_CCM_SW_max)
  real flx_spc_dwn_TOA_CCM(bnd_nbr_CCM_SW_max)
  real foo_CCM(bnd_nbr_CCM_SW_max)
  real odal_O2N2_PUMC_O2_PUMP_N2_CCM(bnd_nbr_CCM_SW_max)
  real odal_O2N2_PUMC_O2_PUMP_N2_CCMCG(bnd_nbr_CCM_SW_max)
  real odal_O2O2_PUMC_O2_PUMP_O2_CCM(bnd_nbr_CCM_SW_max)
  real odal_O2O2_PUMC_O2_PUMP_O2_CCMCG(bnd_nbr_CCM_SW_max)
  real wgt_spc_CCM(bnd_nbr_CCM_SW_max)
  real wvl_ctr_CCM(bnd_nbr_CCM_SW_max)
  real wvl_max_CCM(bnd_nbr_CCM_SW_max)
  real wvl_min_CCM(bnd_nbr_CCM_SW_max)
  real wvl_dlt_CCM(bnd_nbr_CCM_SW_max)
  real xsx_wgt_flx_CCM(bnd_nbr_CCM_SW_max)

  ! Main code
  
  ! Initialize default values
  dbg_lvl=dbg_off

  ! Retrieve command line arguments
  call date_time_get(lcl_date_time)
  call ftn_cmd_ln_sng(cmd_ln)
  call ftn_prg_ID_mk(CVS_Id,CVS_Revision,CVS_Date,prg_ID)
  write (6,'(a)') prg_ID(1:ftn_strlen(prg_ID))
  arg_nbr=command_argument_count()           ! [nbr] Number of command line arguments
  arg_idx=1                 ! [idx] Counting index
  loop_while_options: do while (arg_idx <= arg_nbr)
     call ftn_getarg_wrp(arg_idx,arg_val) ! [sbr] Call getarg, increment arg_idx
     dsh_key=arg_val(1:2)   ! [sng] First two characters of option
     if_dbl_dsh: if (dsh_key == '--') then
        opt_lng=ftn_opt_lng_get(arg_val) ! [nbr] Length of option
        if (opt_lng <= 0) stop 'Long option has no name'
        opt_sng=arg_val(3:2+opt_lng) ! [sng] Option string
        if (opt_sng == 'dbg' .or. opt_sng == 'dbg_lvl' ) then
           call ftn_arg_get(arg_idx,arg_val,dbg_lvl) ! [enm] Debugging level
        else if (opt_sng == 'drc_in') then
           call ftn_arg_get(arg_idx,arg_val,drc_in) ! [sng] Input directory
        else if (opt_sng == 'drc_out') then
           call ftn_arg_get(arg_idx,arg_val,drc_out) ! [sng] Output directory
        else if (opt_sng == 'input' .or. opt_sng == 'fl_O2O2' .or. opt_sng == 'O2O2') then
           call ftn_arg_get(arg_idx,arg_val,fl_in) ! [sng] O2X file
           cmd_ln_fl_in=.true.
        else if (opt_sng == 'HTR16') then
           flg_HTR16=.true.
           flg_GOB90=.false.
        else if (opt_sng == 'GOB90') then
           flg_GOB90=.true.
           flg_HTR16=.false.
        else                ! Option not recognized
           arg_idx=arg_idx-1 ! [idx] Counting index
           call ftn_getarg_err(arg_idx,arg_val) ! [sbr] Error handler for getarg()
        endif               ! endif option is recognized
        ! Jump to top of while loop
        cycle loop_while_options ! C, F77, and F90 use "continue", "goto", and "cycle"
     endif if_dbl_dsh            ! endif long option
     ! Handle short options
     if_sgl_dsh: if (dsh_key == '-3') then
        fl_out_fmt=nf90_format_classic ! [enm] Output file format
     else if (dsh_key == '-4') then
        fl_out_fmt=nf90_format_netcdf4 ! [enm] Output file formaty
     else if (dsh_key == '-D') then
        call ftn_arg_get(arg_idx,arg_val,dbg_lvl)
     else if (dsh_key == '-f') then
        call ftn_arg_get(arg_idx,arg_val,wgt_nm)
     else if (dsh_key == '-G') then
        flg_GOB90=.true.
        flg_HTR16=.false.
     else if (dsh_key == '-i') then
        call ftn_arg_get(arg_idx,arg_val,fl_in)
        cmd_ln_fl_in=.true.
     else if (dsh_key == '-H') then
        flg_HTR16=.true.
        flg_GOB90=.false.
     else if (dsh_key == '-o') then
        call ftn_arg_get(arg_idx,arg_val,fl_out)
        cmd_ln_fl_out=.true.
     else if (dsh_key == '-S') then
        call ftn_arg_get(arg_idx,arg_val,fl_slr)
     else if (dsh_key == '-T') then
        std_tpt=.not.std_tpt
     else if (dsh_key == '-t') then
        call ftn_arg_get(arg_idx,arg_val,tpt_std)
     else if (dsh_key == '-W') then
        WGT_TRN=.true.
     else if (dsh_key == '-w') then
        call ftn_arg_get(arg_idx,arg_val,fl_wgt)
     else if (dsh_key == '-v') then
        write (6,'(a)') CVS_Id
        goto 1000
     else if (dsh_key == '-X') then
        flg_ncl_1530nm_bnd=.not.flg_ncl_1530nm_bnd
     else if (dsh_key == '-x') then
        flg_ncl_1260nm_bnd=.not.flg_ncl_1260nm_bnd
     else                   ! Option not recognized
        arg_idx=arg_idx-1   ! [idx] Counting index
        call ftn_getarg_err(arg_idx,arg_val) ! [sbr] Error handler for getarg()
     endif if_sgl_dsh       ! endif arg_val
  end do loop_while_options ! end while (arg_idx <= arg_nbr)

  ! Compute any quantities that might depend on command line input
  call ftn_strnul(fl_in)
  call ftn_strnul(fl_out)
  call ftn_strnul(fl_slr)
  call ftn_strnul(fl_wgt)
  call ftn_strnul(wgt_nm)
  call ftn_strcpy(src_fl_sng,'Original data file is ' // fl_in)
  if (flg_GOB90) then
     bnd_nbr=bnd_nbr_GOB90
     tpt_cold=tpt_cold_GOB90
     if (.not.cmd_ln_fl_in) fl_in=fl_in_GOB90//nlc
     if (.not.cmd_ln_fl_out) fl_out=fl_out_GOB90//nlc
     call ftn_strcpy(src_rfr_sng,'Data reference is Greenblatt et al. (1990) (GOB90)')
     if (flg_ncl_1260nm_bnd) then
        write (6,'(a)') 'Including 1.26 micron absorption band'
     else
        write (6,'(a)') 'Excluding 1.26 micron absorption band'
     endif ! flg_ncl_1260nm_bnd
     if (flg_ncl_1530nm_bnd) then
        write (6,'(a)') 'Including 1.53 micron absorption band seen by MCB98'
        stop 'ERROR: 1.53 micron absorption band not supported yet'
     else
        write (6,'(a)') 'Excluding 1.53 micron absorption band seen by MCB98'
     endif ! flg_ncl_1530nm_bnd
  endif ! flg_GOB90
  if (flg_HTR16) then
     bnd_nbr=bnd_nbr_HTR16
     tpt_cold=tpt_cold_HTR16
     if (.not.cmd_ln_fl_in) fl_in=fl_in_HTR16_cold//nlc
     if (.not.cmd_ln_fl_out) fl_out=fl_out_HTR16//nlc
     call ftn_strcpy(src_rfr_sng,'Data reference is HITRAN (2017) (HTR16)')
  endif ! flg_HTR16
  if (WGT_TRN) then
     call ftn_strcpy(src_wgt_sng,'CCM cross sections are averages of high resolution cross-sections ' &
             //'weighted by variable '//wgt_nm(1:ftn_strlen(wgt_nm))// &
             ' from model atmosphere in '//fl_wgt(1:ftn_strlen(fl_wgt))//char(0))
  else
     call ftn_strcpy(src_wgt_sng,'CCM cross sections are averages of high resolution cross-sections ' &
             //'weighted by TOA solar spectral flux'//char(0))
  endif ! WGT_TRN

  ! Compute quantities that may depend on command line input
  ! Prepend user-specified path, if any, to input data file names
  if (ftn_strlen(drc_in) > 0) then
     call ftn_drcpfx(drc_in,fl_in) ! [sng] Input file
     call ftn_drcpfx(drc_in,fl_slr) ! [sng] Solar spectrum file
     call ftn_drcpfx(drc_in,fl_wgt) ! [sng] Weight file
  endif                     ! endif drc_in
  ! Prepend user-specified path, if any, to output data file names
  if (ftn_strlen(drc_out) > 0) call ftn_drcpfx(drc_out,fl_out) ! [sng] Output file

  ! Open file to read header and obtain bnd_nbr if necessary
  open (fl_in_unit,file=fl_in,status='old',iostat=rcd)

  if (flg_HTR16) then            ! HTR16 data

     ! HITRAN data are in .cia format described above
     ! First, read-in bnd_nbr in order to allocate memory
     read (fl_in_unit,'(a20,f10.3,f10.3,i7,f7.3,e10.3,f6.0,a27,i3)') &
          mlc_frm_htrn,wvn_min_htrn,wvn_max_htrn,bnd_nbr_htrn,tpt_htrn, &
          xsx_max_htrn,rsn_nst_htrn,cmt_htrn,rfr_nbr_htrn

     ! Sanity check
     if (dbg_lvl >= dbg_fl) then
        write (6,'(a7,a20,f10.3,f10.3,i7,f7.3,e10.3,f6.0,a27,i3)') &
             'DEBUG: ',mlc_frm_htrn,wvn_min_htrn,wvn_max_htrn,bnd_nbr_htrn,tpt_htrn, &
             xsx_max_htrn,rsn_nst_htrn,cmt_htrn,rfr_nbr_htrn
     endif                     ! endif dbg

     tpt_cold=tpt_htrn
     bnd_nbr=bnd_nbr_htrn
     write (6,'(a16,1x,a)') 'Read header from',fl_in(1:ftn_strlen(fl_in))

  endif                     ! HTR16 data

  ! Allocate space for dynamic arrays
  allocate(abs_cff_mss_O2O2(bnd_nbr),stat=rcd)
  allocate(abs_xsx_O2N2(bnd_nbr),stat=rcd)
  allocate(abs_xsx_O2O2(bnd_nbr),stat=rcd)
  allocate(abs_xsx_O2O2_dadT(bnd_nbr),stat=rcd)
  allocate(abs_xsx_O2O2_tpt_rfr(bnd_nbr),stat=rcd)
  allocate(bnd(bnd_nbr),stat=rcd)     ! coordinate variable
  allocate(flx_bnd_dwn_TOA(bnd_nbr),stat=rcd)
  allocate(flx_bnd_pht_dwn_TOA(bnd_nbr),stat=rcd)
  allocate(flx_slr_frc(bnd_nbr),stat=rcd)
  allocate(flx_spc_dwn_TOA(bnd_nbr),stat=rcd)
  allocate(idx_rfr_air_STP(bnd_nbr),stat=rcd)
  allocate(nrg_pht(bnd_nbr),stat=rcd)
  allocate(odal_O2N2_PUMC_O2_PUMP_N2(bnd_nbr),stat=rcd)
  allocate(odal_O2N2_PUNC_O2_PUNP_N2(bnd_nbr),stat=rcd)
  allocate(odal_O2O2_PUMC_O2_PUMP_O2(bnd_nbr),stat=rcd)
  allocate(odal_O2O2_PUNC_O2_PUNP_O2(bnd_nbr),stat=rcd)
  allocate(wgt_spc(bnd_nbr),stat=rcd)
  allocate(wvl(bnd_nbr),stat=rcd)     ! coordinate variable
  allocate(wvl_ctr(bnd_nbr),stat=rcd)
  allocate(wvl_dlt(bnd_nbr),stat=rcd)
  allocate(wvl_grd(bnd_nbr+1),stat=rcd)
  allocate(wvl_max(bnd_nbr),stat=rcd)
  allocate(wvl_min(bnd_nbr),stat=rcd)
  allocate(wvn(bnd_nbr),stat=rcd)     ! coordinate variable
  allocate(wvn_ctr(bnd_nbr),stat=rcd)
  allocate(wvn_dlt(bnd_nbr),stat=rcd)
  allocate(wvn_grd(bnd_nbr+1),stat=rcd)
  allocate(wvn_max(bnd_nbr),stat=rcd)
  allocate(wvn_min(bnd_nbr),stat=rcd)
  allocate(xsx_wgt_flx(bnd_nbr),stat=rcd)
  
  if (flg_GOB90) then

     ! File is already open, read all data
     do bnd_idx=1,7
        read (fl_in_unit,'(a80)') lbl
     enddo ! bnd_idx
     lbl(1:1)=lbl(1:1) ! CEWI
     do bnd_idx=1,bnd_nbr
        read (fl_in_unit,*) &
                   wvl_ctr(bnd_idx), & ! [nm]
                   odal_O2O2_PUNC_O2_PUNP_O2(bnd_idx) ! [cm5 mlc-2]
     enddo ! bnd_idx
     
     ! Convert input grid to SI units where necessary
     do bnd_idx=1,bnd_nbr
        wvl_ctr(bnd_idx)=wvl_ctr(bnd_idx)*1.0e-9 ! nm -> m
     enddo ! bnd_idx

     ! NB: The measurement of O2-O2 cross sections leads to non-standard nomenclature
     
     ! Shardanand (1969) is convinced O2-O2 absorption is attributable to covalently bonded O4
     ! Shardanand (1977) is convinced O2-O2, O2-N2, O2-Ar absorption is attributable to van der Waals complexes
     ! GOB90 and Blake and McCoy (1987) believe O2-O2 absorption is attributable to collision-induced absorption
     ! Short lived collision pairs, not bound dimers, are conclusively shown to be the cause of O2-O2 and O2-N2 absorption in SPS98
     ! Thus concentration of O2-O2 the concentration of O2-O2 collision pairs
     ! Unfortunately, no one knows how to directly calculate the concentration of O2-O2
     ! The physical meaningfulness of such a concentration is not clear in any case
     ! Labs measure total monochromatic absorption cross section at given O2 pressure
     ! Collision-pair induced absorption is inferred from pressure variation of this absorption
     ! O2 concentration is increased by increasing O2 partial pressure
     ! Change in measured absorption cross section is then correlated with square of O2 concentration
     ! GOB90 show once absorption that varies linearly with O2 concentration is removed, residual absorption correlates nearly exactly with square of O2 concentration
     ! This residual absorption is due to O2-O2 collision complexes
     ! Unfortunately, only convolution of O2-O2 concentration and O2-O2 absorption is known with certainty
     ! O2-O2 concentration and O2-O2 absorption cross sections are not known individually
     ! How to implement this absorption in RT codes?
     
     ! abs_xsx_O2O2: the standard measurement of molecular cross-section
     ! Multiply abs_xsx_O2O2 by column path of O2-O2 complexes (# m-2) to obtain absorption optical depth of O2-O2
     ! abs_xsx_O2O2 has the advantage of physical simplicity; it measures O2-O2 absorption in terms of O2-O2 concentration
     ! This allows, e.g., O2-O2 absorption cross-sections to be directly compared to, say, O3 cross-sections
     ! Disadvantage is the uncertainty of O2-O2 concentration
     ! Shardanand (Sha77 p. 436, 527) quotes an interaction coefficient for O2 + O2 <-> O2-O2 of K = 2.6e-23 cm3 mlc-1
     ! GOB90 are skeptical of applying any equilibrium reaction coefficients to back out O2-O2 concentration from O2 concentration
     ! Since SPS98 and MCB98 agree with GOB90, I disavow using O2-O2 concentrations on physical principle
     ! I continue to use faux O2-O2 concentrations in swnb because changing clm and swnb to use either of the following more physical quantities will be time-consuming
     
     ! The first physically based method is based on O2 number concentration
     ! odal_O2O2_PUNC_O2_PUNP_O2: measured absorption optical depth of O2-O2 per unit number concentration of O2 per unit number path of O2
     ! This is the quantity measured by GOB90
     ! PUNC stands for "per unit number concentration"
     ! PUNP stands for "per unit number path"
     ! Multiply odal_O2O2_PUNC_O2_PUNP_O2 by number concentration of O2 molecules (# m-3) to obtain absorption cross-section of O2-O2 per molecule O2
     ! Multiply odal_O2O2_PUNC_O2_PUNP_O2 by O2 number concentration (# m-3) then by O2 number path (# m-2) to obtain absorption optical depth of O2-O2
     ! Advantage of this method is it does not rely on uncertain equilibrium rate coefficients for O2 + O2 <-> O2-O2
     ! There are two disadvantages to this method:
     ! 1. Computation of O2-O2 absorption optical depth in RT program is slightly non-standard because it requires knowing O2 concentration and path in the O2-O2 loop
     ! 2. odal_O2O2_PUNC_O2_PUNP_O2 is O(1.0e-46) and thus requires double precision (eight byte) representation,
     ! Remainder of RT code has (barely) managed to avoid double precision until now
     
     ! The second, and ultimately better, physically based method is based on O2 mass concentration
     ! odal_O2O2_PUMC_O2_PUMP_O2: measured absorption optical depth of O2-O2 per unit mass concentration of O2 per unit mass path of O2
     ! odal_O2O2_PUMC_O2_PUMP_O2 translates to "layer absorption optical depth of O2-O2 per unit mass concentration of O2 per unit mass path of O2"
     ! PUMC stands for "per unit mass concentration"
     ! PUMP stands for "per unit mass path"
     ! Multiply odal_O2O2_PUMC_O2_PUMP_O2 by mass concentration (kg m-3) of O2 (not O2-O2) to obtain absorption optical depth of O2-O2 per unit mass path of O2
     ! Multiply odal_O2O2_PUMC_O2_PUMP_O2 by mass concentration (kg m-3) of O2 (not O2-O2) then by mass path (kg m-2) of O2 (not O2-O2) to obtain absorption optical depth of O2-O2
     ! A main advantage of odal_O2O2_PUMC_O2_PUMP_O2 is that using it requires no assumptions about O2-O2 concentration
     ! All one needs to know to use odal_O2O2_PUMC_O2_PUMP_O2 is the abundance of O2, which is trivial to derive
     ! Another advantage is that odal_O2O2_PUMC_O2_PUMP_O2 is always representable in single precision (4 byte) storage

     ! Compute diagnostic grid variables
     wvl_grd(1)=wvl_ctr(1)-0.5*(wvl_ctr(2)-wvl_ctr(1))
     do bnd_idx=2,bnd_nbr
        wvl_grd(bnd_idx)=0.5*(wvl_ctr(bnd_idx-1)+wvl_ctr(bnd_idx))
     enddo ! bnd_idx
     wvl_grd(bnd_nbr+1)=wvl_ctr(bnd_nbr)+0.5*(wvl_ctr(bnd_nbr)-wvl_ctr(bnd_nbr-1))

     do bnd_idx=1,bnd_nbr
        wvl_min(bnd_idx)=wvl_grd(bnd_idx)
        wvl_max(bnd_idx)=wvl_grd(bnd_idx+1)
        wvl_ctr(bnd_idx)=0.5*(wvl_min(bnd_idx)+wvl_max(bnd_idx))
        wvl(bnd_idx)=wvl_ctr(bnd_idx)
        wvl_dlt(bnd_idx)=wvl_max(bnd_idx)-wvl_min(bnd_idx)
     enddo

     do bnd_idx=1,bnd_nbr
        wvn_min(bnd_idx)=1.0/(100.0*wvl_max(bnd_idx))
        wvn_max(bnd_idx)=1.0/(100.0*wvl_min(bnd_idx))
        wvn_dlt(bnd_idx)=wvn_max(bnd_idx)-wvn_min(bnd_idx)
     enddo
     
  else if (flg_HTR16) then ! HTR16 data

     do bnd_idx=1,bnd_nbr
        read (fl_in_unit,*) &
                   wvn_ctr(bnd_idx), & ! [cm-1]
                   odal_O2O2_PUNC_O2_PUNP_O2(bnd_idx) ! [cm5 mlc-2]
     enddo ! bnd_idx

     ! Sanity check
     if (dbg_lvl >= dbg_fl) then
        do bnd_idx=1,bnd_nbr
           write (6,'(a10,f10.3,a36,es10.3,a10)') &
                'wvn_ctr = ',wvn_ctr(bnd_idx),' cm-1, odal_O2O2_PUNC_O2_PUNP_O2 = ',odal_O2O2_PUNC_O2_PUNP_O2(bnd_idx),' cm5 mlc-2' 
        enddo
     endif                     ! endif dbg

     ! Convert input grid to SI units where necessary (not necessary for HTR16)

     ! O2-O2 from reference 2 have 0.2 cm-1 resolution except at 1385.9 cm-1 which has 0.1 cm-1 resolution
     ! Hence 4002 points cover 800 cm-1, wvn_rsn is not uniform throughout the measurements
     ! wvn_rsn=(wvn_max_htrn-wvn_min_htrn)/(bnd_nbr-1)
     wvn_grd(1)=wvn_ctr(1)-0.5*(wvn_ctr(2)-wvn_ctr(1))
     do bnd_idx=2,bnd_nbr
        wvn_grd(bnd_idx)=0.5*(wvn_ctr(bnd_idx-1)+wvn_ctr(bnd_idx))
     enddo
     wvn_grd(bnd_nbr+1)=wvn_ctr(bnd_nbr)+0.5*(wvn_ctr(bnd_nbr)-wvn_ctr(bnd_nbr-1))

     do bnd_idx=1,bnd_nbr
        wvn_min(bnd_idx)=wvn_grd(bnd_idx)
        wvn_max(bnd_idx)=wvn_grd(bnd_idx+1)
        wvn_ctr(bnd_idx)=0.5*(wvn_min(bnd_idx)+wvn_max(bnd_idx))
        wvn(bnd_idx)=wvn_ctr(bnd_idx)
        wvn_dlt(bnd_idx)=wvn_max(bnd_idx)-wvn_min(bnd_idx)
     enddo
     do bnd_idx=1,bnd_nbr
        wvl_min(bnd_idx)=1.0/(100.0*wvn_max(bnd_idx))
        wvl_max(bnd_idx)=1.0/(100.0*wvn_min(bnd_idx))
        wvl_dlt(bnd_idx)=wvl_max(bnd_idx)-wvl_min(bnd_idx)
        wvl_ctr(bnd_idx)=0.5*(wvl_min(bnd_idx)+wvl_max(bnd_idx))
     enddo

  endif                     ! HTR16 data

  close (fl_in_unit)
  write (6,'(a20,1x,a)') 'Read input data from',fl_in(1:ftn_strlen(fl_in))
  
  ! 20181210: Data processing for GOB90 and HTR16 data are now shared
  ! Convert input data to SI units where necessary
  do bnd_idx=1,bnd_nbr
     odal_O2O2_PUNC_O2_PUNP_O2(bnd_idx)=odal_O2O2_PUNC_O2_PUNP_O2(bnd_idx)*1.0e-10 ! [cm5 mlc-2] -> [m5 mlc-2]
  enddo ! bnd_idx

  ! Sanity check: HITRAN data has numerous negative values presumably due to noisy instrumentation
  do bnd_idx=1,bnd_nbr
     odal_O2O2_PUNC_O2_PUNP_O2(bnd_idx)=max(odal_O2O2_PUNC_O2_PUNP_O2(bnd_idx),0.0)
  enddo ! bnd_idx
  
  ! Zero the 1260 nm band only if requested
  if (.not.flg_ncl_1260nm_bnd) then
     do bnd_idx=1,bnd_nbr
        ! The 1.26 micron band is currently the only O2-O2 absorption beyond 1.137 microns
        if (wvl_ctr(bnd_idx) >= 1.137e-6.and.wvl_ctr(bnd_idx) <= 1.350e-6) then
           odal_O2O2_PUNC_O2_PUNP_O2(bnd_idx)=0.0 ! [m5 mlc-2]
        endif ! wvl_ctr
     enddo ! bnd_idx
  endif ! flg_ncl_1260nm_bnd
  
  do bnd_idx=1,bnd_nbr
     wvl(bnd_idx)=wvl_ctr(bnd_idx)
     bnd(bnd_idx)=wvl_ctr(bnd_idx)
     ! Here we used to normalize by an equilibrium constant to keep numbers in single precision
     !        double_foo=odal_O2O2_PUNC_O2_PUNP_O2(bnd_idx)/dble(k_O2_O2)
     double_foo=odal_O2O2_PUNC_O2_PUNP_O2(bnd_idx) ! defensive programming
     abs_xsx_O2O2(bnd_idx)=double_foo ! [m5 mlc-2]
     abs_xsx_O2O2_tpt_rfr(bnd_idx)=tpt_std ! [K] All GOB90 data taken at 296 K
     ! Temperature dependence for GOB90 data not currently implemented
     abs_xsx_O2O2_dadT(bnd_idx)=0.0 ! [m5 mlc-2 K-1] GOB90 showed very little temperature dependence 
     flx_bnd_pht_dwn_TOA(bnd_idx)=mss_val ! #/cm2/s -> #/m2/s
  enddo ! bnd_idx
  
  ! Set O2-N2 collision-induced absorption cross sections
  do bnd_idx=1,bnd_nbr
     if (wvl_ctr(bnd_idx) >= 1.137e-6.and.wvl_ctr(bnd_idx) <= 1.350e-6) then
        odal_O2N2_PUNC_O2_PUNP_N2(bnd_idx)=N2_cll_prt_fsh*odal_O2O2_PUNC_O2_PUNP_O2(bnd_idx) ! [m5 mlc-2]
        abs_xsx_O2N2(bnd_idx)=N2_cll_prt_fsh*abs_xsx_O2O2(bnd_idx) ! [m5 mlc-2]
     else                ! endif wvl
        odal_O2N2_PUNC_O2_PUNP_N2(bnd_idx)=0.0 ! [m5 mlc-2]
        abs_xsx_O2N2(bnd_idx)=0.0 ! [m5 mlc-2]
     endif               ! endelse wvl
  enddo ! bnd_idx
  
  ! Get TOA solar spectrum
  slr_spc_xtr_typ=xtr_fll_ngh !+xtr_vrb
  call slr_spc_get(fl_slr,wvl_min,wvl_max,bnd_nbr,flx_slr_frc,slr_spc_xtr_typ,slr_spc_xtr_typ)
  
  ! Compute diagnostic variables
  do bnd_idx=1,bnd_nbr
     abs_cff_mss_O2O2(bnd_idx)=abs_xsx_O2O2(bnd_idx)*(dble(Avagadro)/dble(mmw_O2))**2 ! [m5 kg-2]
     double_foo= &
          odal_O2O2_PUNC_O2_PUNP_O2(bnd_idx)*(dble(Avagadro)/dble(mmw_O2))**2 ! [m5 kg-2]
     odal_O2O2_PUMC_O2_PUMP_O2(bnd_idx)=double_foo ! [m5 kg-2] Defensive programming
     double_foo= &
          odal_O2N2_PUNC_O2_PUNP_N2(bnd_idx)*dble(Avagadro)*dble(Avagadro)/(dble(mmw_O2)*dble(mmw_N2))
     odal_O2N2_PUMC_O2_PUMP_N2(bnd_idx)=double_foo ! [m5 kg-2] Defensive programming
     nrg_pht(bnd_idx)=Planck*speed_of_light/wvl_ctr(bnd_idx)
     flx_bnd_dwn_TOA(bnd_idx)=flx_slr_frc(bnd_idx)*slr_cst_CCM
     flx_spc_dwn_TOA(bnd_idx)=flx_slr_frc(bnd_idx)*slr_cst_CCM/wvl_dlt(bnd_idx)
     flx_bnd_pht_dwn_TOA(bnd_idx)=flx_bnd_dwn_TOA(bnd_idx)/nrg_pht(bnd_idx)
  enddo ! bnd_idx
  
  ! Compute index of refraction through dry air at STP (Len93 p. 155)
  do bnd_idx=1,bnd_nbr
     idx_rfr_air_STP(bnd_idx)= &
          1.0+ &
          1.0e-6*(77.46+0.459/(1.0e12*wvl_ctr(bnd_idx)**2))* &
          prs_STP*0.01/tpt_STP
  enddo ! bnd_idx
  
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
  
  ! Variable definitions
  rcd=nf90_wrp(nf90_def_var(nc_id,'tpt_cold',nf90_float,tpt_cold_id),sbr_nm//': dv tpt_cold')
  rcd=nf90_wrp(nf90_def_var(nc_id,'tpt_warm',nf90_float,tpt_warm_id),sbr_nm//': dv tpt_warm')
  rcd=nf90_wrp(nf90_def_var(nc_id,'tpt_std',nf90_float,tpt_std_id),sbr_nm//': dv tpt_std')
  rcd=nf90_wrp(nf90_def_var(nc_id,'abs_cff_mss_O2O2',nf90_float,bnd_dmn_id,abs_cff_mss_O2O2_id),sbr_nm//': dv abs_cff_mss_O2O2')
  rcd=nf90_wrp(nf90_def_var(nc_id,'abs_xsx_O2N2',nf90_float,bnd_dmn_id,abs_xsx_O2N2_id),sbr_nm//': dv abs_xsx_O2N2')
  rcd=nf90_wrp(nf90_def_var(nc_id,'abs_xsx_O2O2',nf90_double,bnd_dmn_id,abs_xsx_O2O2_id),sbr_nm//': dv abs_xsx_O2O2')
  rcd=nf90_wrp(nf90_def_var(nc_id,'abs_xsx_O2O2_dadT',nf90_float,bnd_dmn_id,abs_xsx_O2O2_dadT_id),sbr_nm//': dv abs_xsx_O2O2_dadT')
  rcd=nf90_wrp(nf90_def_var(nc_id,'bnd',nf90_float,bnd_dmn_id,bnd_id),sbr_nm//': dv bnd')
  rcd=nf90_wrp(nf90_def_var(nc_id,'flx_bnd_dwn_TOA',nf90_float,bnd_dmn_id,flx_bnd_dwn_TOA_id),sbr_nm//': dv flx_bnd_dwn_TOA')
  rcd=nf90_wrp(nf90_def_var(nc_id,'flx_slr_frc',nf90_float,bnd_dmn_id,flx_slr_frc_id),sbr_nm//': dv flx_slr_frc')
  rcd=nf90_wrp(nf90_def_var(nc_id,'flx_spc_dwn_TOA',nf90_float,bnd_dmn_id,flx_spc_dwn_TOA_id),sbr_nm//': dv flx_spc_dwn_TOA')
  rcd=nf90_wrp(nf90_def_var(nc_id,'idx_rfr_air_STP',nf90_float,bnd_dmn_id,idx_rfr_air_STP_id),sbr_nm//': dv idx_rfr_air_STP')
  rcd=nf90_wrp(nf90_def_var(nc_id,'nrg_pht',nf90_float,bnd_dmn_id,nrg_pht_id),sbr_nm//': dv nrg_pht')
  rcd=nf90_wrp(nf90_def_var(nc_id,'wvl_ctr',nf90_float,bnd_dmn_id,wvl_ctr_id),sbr_nm//': dv wvl_ctr')
  rcd=nf90_wrp(nf90_def_var(nc_id,'wvl_grd',nf90_float,grd_dmn_id,wvl_grd_id),sbr_nm//': dv wvl_grd')
  rcd=nf90_wrp(nf90_def_var(nc_id,'wvl_max',nf90_float,bnd_dmn_id,wvl_max_id),sbr_nm//': dv wvl_max')
  rcd=nf90_wrp(nf90_def_var(nc_id,'wvl_min',nf90_float,bnd_dmn_id,wvl_min_id),sbr_nm//': dv wvl_min')
  rcd=nf90_wrp(nf90_def_var(nc_id,'wvl_dlt',nf90_float,bnd_dmn_id,wvl_dlt_id),sbr_nm//': dv wvl_dlt')
  rcd=nf90_wrp(nf90_def_var(nc_id,'wvn_ctr',nf90_float,bnd_dmn_id,wvn_ctr_id),sbr_nm//': dv wvn_ctr')
  rcd=nf90_wrp(nf90_def_var(nc_id,'wvn_grd',nf90_float,grd_dmn_id,wvn_grd_id),sbr_nm//': dv wvn_grd')
  rcd=nf90_wrp(nf90_def_var(nc_id,'wvn_max',nf90_float,bnd_dmn_id,wvn_max_id),sbr_nm//': dv wvn_max')
  rcd=nf90_wrp(nf90_def_var(nc_id,'wvn_min',nf90_float,bnd_dmn_id,wvn_min_id),sbr_nm//': dv wvn_min')
  rcd=nf90_wrp(nf90_def_var(nc_id,'wvn_dlt',nf90_float,bnd_dmn_id,wvn_dlt_id),sbr_nm//': dv wvn_dlt')
  ! Wrap
  rcd=nf90_wrp(nf90_def_var(nc_id,'flx_bnd_pht_dwn_TOA',nf90_float,bnd_dmn_id,flx_bnd_pht_dwn_TOA_id), &
       sbr_nm//': dv flx_bnd_pht_dwn_TOA')
  rcd=nf90_wrp(nf90_def_var(nc_id,'abs_xsx_O2O2_tpt_rfr',nf90_float,bnd_dmn_id,abs_xsx_O2O2_tpt_rfr_id), &
       sbr_nm//': dv abs_xsx_O2O2_tpt_rfr')
  rcd=nf90_wrp(nf90_def_var(nc_id,'odal_O2N2_PUMC_O2_PUMP_N2',nf90_float,bnd_dmn_id,odal_O2N2_PUMC_O2_PUMP_N2_id), &
       sbr_nm//': dv odal_O2N2_PUMC_O2_PUMP_N2')
  rcd=nf90_wrp(nf90_def_var(nc_id,'odal_O2O2_PUMC_O2_PUMP_O2',nf90_float,bnd_dmn_id,odal_O2O2_PUMC_O2_PUMP_O2_id), &
       sbr_nm//': dv odal_O2O2_PUMC_O2_PUMP_O2')
  rcd=nf90_wrp(nf90_def_var(nc_id,'odal_O2N2_PUNC_O2_PUNP_N2',nf90_double,bnd_dmn_id,odal_O2N2_PUNC_O2_PUNP_N2_id), &
       sbr_nm//': dv odal_O2N2_PUNC_O2_PUNP_N2')
  rcd=nf90_wrp(nf90_def_var(nc_id,'odal_O2O2_PUNC_O2_PUNP_O2',nf90_double,bnd_dmn_id,odal_O2O2_PUNC_O2_PUNP_O2_id), &
       sbr_nm//': dv odal_O2O2_PUNC_O2_PUNP_O2')
  
  ! Add global attributes
  rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'CVS_Id',CVS_Id),sbr_nm//': pa CVS_Id')
  rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'creation_date',lcl_date_time),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'prg_ID',prg_ID(1:ftn_strlen(prg_ID))),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'cmd_ln',cmd_ln(1:ftn_strlen(cmd_ln))),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'src_rfr_sng',src_rfr_sng(1:ftn_strlen(src_rfr_sng))),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'src_fl_sng',src_fl_sng(1:ftn_strlen(src_fl_sng))),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'src_wgt_sng',src_wgt_sng(1:ftn_strlen(src_wgt_sng))),sbr_nm)
  
  ! Add english text descriptions
  rcd=nf90_wrp(nf90_put_att(nc_id,tpt_cold_id,'long_name','Temperature of coldest CFC11 measurements employed'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,tpt_warm_id,'long_name','Temperature of warmest CFC11 measurements employed'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,tpt_std_id,'long_name','Temperature at which interpolated cross sections are archived'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,abs_cff_mss_O2O2_id,'long_name','O2-O2 mass absorption coefficient (per kg O2-O2)'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,abs_xsx_O2O2_tpt_rfr_id,'long_name','Valid temperature for absorption cross section'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,bnd_id,'long_name','Band nominal wavelength'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,flx_bnd_dwn_TOA_id,'long_name','Solar Energy flux in band'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,flx_bnd_pht_dwn_TOA_id,'long_name','Photon flux in band'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,flx_slr_frc_id,'long_name','Fraction of solar flux in band: ' // fl_slr),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,flx_spc_dwn_TOA_id,'long_name','Spectral solar insolation at TOA'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,idx_rfr_air_STP_id,'long_name','Index of refraction at band center at STP'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,nrg_pht_id,'long_name','Energy of photon at band center'),sbr_nm)
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
  ! Wrap
  rcd=nf90_wrp(nf90_put_att(nc_id,abs_xsx_O2O2_dadT_id,'long_name', &
       'Slope of absorption cross section temperature dependence'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,abs_xsx_O2O2_id,'long_name', &
       'O2-O2 collision-induced absorption cross section'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,abs_xsx_O2N2_id,'long_name', &
       'O2-N2 collision-induced absorption cross section'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,odal_O2O2_PUMC_O2_PUMP_O2_id,'long_name', &
       'O2-O2 layer absorption optical depth (per kg m-3 O2 per kg m-2 O2)'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,odal_O2N2_PUMC_O2_PUMP_N2_id,'long_name', &
       'O2-N2 layer absorption optical depth (per kg m-3 O2 per kg m-2 N2)'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,odal_O2O2_PUNC_O2_PUNP_O2_id,'long_name', &
       'O2-O2 layer absorption optical depth (per mlc m-3 O2 per mlc m-2 O2)'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,odal_O2N2_PUNC_O2_PUNP_N2_id,'long_name', &
       'O2-N2 layer absorption optical depth (per mlc m-3 O2 per mlc m-2 N2)'),sbr_nm)
  
  ! Add units
  rcd=nf90_wrp(nf90_put_att(nc_id,tpt_cold_id,'units','kelvin'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,tpt_warm_id,'units','kelvin'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,tpt_std_id,'units','kelvin'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,abs_cff_mss_O2O2_id,'units','meter5 kilogram-2'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,abs_xsx_O2N2_id,'units','meter5 molecule-2'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,abs_xsx_O2O2_dadT_id,'units','meter5 molecule-2 kelvin-1'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,abs_xsx_O2O2_id,'units','meter5 molecule-2'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,abs_xsx_O2O2_tpt_rfr_id,'units','kelvin'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,bnd_id,'units','meter'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,flx_bnd_dwn_TOA_id,'units','watt meter-2'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,flx_bnd_pht_dwn_TOA_id,'units','photon meter-2 second-1'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,flx_slr_frc_id,'units','fraction'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,flx_spc_dwn_TOA_id,'units','watt meter-2 meter-1'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,idx_rfr_air_STP_id,'units','fraction'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,nrg_pht_id,'units','joule photon-1'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,odal_O2N2_PUMC_O2_PUMP_N2_id,'units','meter5 kilogram-2'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,odal_O2N2_PUNC_O2_PUNP_N2_id,'units','meter5 molecule-2'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,odal_O2O2_PUMC_O2_PUMP_O2_id,'units','meter5 kilogram-2'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,odal_O2O2_PUNC_O2_PUNP_O2_id,'units','meter5 molecule-2'),sbr_nm)
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
  rcd=nf90_wrp(nf90_put_var(nc_id,tpt_cold_id,tpt_cold),sbr_nm//': pv tpt_cold in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,tpt_warm_id,tpt_warm),sbr_nm//': pv tpt_warm in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,tpt_std_id,tpt_std),sbr_nm//': pv tpt_std in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,abs_cff_mss_O2O2_id,abs_cff_mss_O2O2),sbr_nm//': pv abs_cff_mss_O2O2 in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,abs_xsx_O2N2_id,abs_xsx_O2N2),sbr_nm//': pv abs_xsx_O2N2 in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,abs_xsx_O2O2_dadT_id,abs_xsx_O2O2_dadT),sbr_nm//': pv abs_xsx_O2O2_dadT in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,abs_xsx_O2O2_id,abs_xsx_O2O2),sbr_nm//': pv abs_xsx_O2O2 in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,abs_xsx_O2O2_tpt_rfr_id,abs_xsx_O2O2_tpt_rfr),sbr_nm//': pv abs_xsx_O2O2_tpt_rfr in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,bnd_id,bnd),sbr_nm//': pv bnd in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,flx_bnd_dwn_TOA_id,flx_bnd_dwn_TOA),sbr_nm//': pv flx_bnd_dwn_TOA in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,flx_bnd_pht_dwn_TOA_id,flx_bnd_pht_dwn_TOA),sbr_nm//': pv flx_bnd_pht_dwn_TOA in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,flx_slr_frc_id,flx_slr_frc),sbr_nm//': pv flx_slr_frc in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,flx_spc_dwn_TOA_id,flx_spc_dwn_TOA),sbr_nm//': pv flx_spc_dwn_TOA in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,idx_rfr_air_STP_id,idx_rfr_air_STP),sbr_nm//': pv idx_rfr_air_STP in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,nrg_pht_id,nrg_pht),sbr_nm//': pv nrg_pht in '//__FILE__)
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
  ! Wrap
  rcd=nf90_wrp(nf90_put_var(nc_id,odal_O2N2_PUNC_O2_PUNP_N2_id,odal_O2N2_PUNC_O2_PUNP_N2), &
       sbr_nm//': pv odal_O2N2_PUNC_O2_PUNP_N2 in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,odal_O2O2_PUNC_O2_PUNP_O2_id,odal_O2O2_PUNC_O2_PUNP_O2), &
       sbr_nm//': pv odal_O2O2_PUNC_O2_PUNP_O2 in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,odal_O2N2_PUMC_O2_PUMP_N2_id,odal_O2N2_PUMC_O2_PUMP_N2), &
       sbr_nm//': pv odal_O2N2_PUMC_O2_PUMP_N2 in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,odal_O2O2_PUMC_O2_PUMP_O2_id,odal_O2O2_PUMC_O2_PUMP_O2), &
       sbr_nm//': pv odal_O2O2_PUMC_O2_PUMP_O2 in '//__FILE__)
  
  rcd=nf90_close(nc_id)
  write (6,'(a28,1x,a)') 'Wrote results to netCDF file',fl_out(1:ftn_strlen(fl_out))
  
  ! Get CCM wavelength grid
  call wvl_grd_CCM_SW_mk(bnd_nbr_CCM,bnd_CCM,wvl_min_CCM,wvl_max_CCM,wvl_ctr_CCM,wvl_dlt_CCM)
  
  ! Get TOA solar spectrum
  call slr_spc_get_CCM(fl_slr,wvl_min_CCM,wvl_max_CCM,bnd_nbr_CCM,flx_slr_frc_CCM,slr_spc_xtr_typ,slr_spc_xtr_typ)
  
  if (WGT_TRN) then
     ! Weight high resolution absorption cross-sections by atmospheric transmission of atmosphere with all constituents except O2O2
     call wgt_get(fl_wgt,wgt_nm,wvl_grd,bnd_nbr,wgt_spc)         
     do bnd_idx=1,bnd_nbr
        xsx_wgt_flx(bnd_idx)=(mlt_fct*abs_xsx_O2O2(bnd_idx))*wgt_spc(bnd_idx)
     enddo ! bnd_idx
  else                      ! !WGT_TRN
     ! Weight high resolution absorption cross sections by high resolution TOA solar flux
     do bnd_idx=1,bnd_nbr
        xsx_wgt_flx(bnd_idx)=(mlt_fct*abs_xsx_O2O2(bnd_idx))*flx_spc_dwn_TOA(bnd_idx)
     enddo ! bnd_idx
  endif                     ! !WGT_TRN
  
  ! Initialize default values
  xtr_typ_LHS=xtr_prt_wgt+xtr_fll_nil !+xtr_vrb
  xtr_typ_RHS=xtr_prt_wgt+xtr_fll_nil !+xtr_vrb
  
  ! Rebin solar flux
  call rbn_vec_CCM(bnd_nbr,wvl_grd,flx_spc_dwn_TOA, &
       bnd_nbr_CCM,wvl_min_CCM,wvl_max_CCM,foo_CCM, &
       xtr_typ_LHS,xtr_typ_RHS)
  ! Rebin flux-weighted absorption cross sections
  call rbn_vec_CCM(bnd_nbr,wvl_grd,xsx_wgt_flx, &
       bnd_nbr_CCM,wvl_min_CCM,wvl_max_CCM,xsx_wgt_flx_CCM, &
       xtr_typ_LHS,xtr_typ_RHS)
  if (WGT_TRN) then
     ! Rebin cross-section weights
     call rbn_vec_CCM(bnd_nbr,wvl_grd,wgt_spc, &
             bnd_nbr_CCM,wvl_min_CCM,wvl_max_CCM,wgt_spc_CCM, &
             xtr_typ_LHS,xtr_typ_RHS)
  endif                     ! endif WGT_TRN
  
  ! Multiply double (above) then divide single by large number to avoid underflow
  do bnd_idx=1,bnd_nbr_CCM
     xsx_wgt_flx_CCM(bnd_idx)=xsx_wgt_flx_CCM(bnd_idx)/mlt_fct
     ! NB: next line not necessary since xsx_wgt_flx is not used again
     ! xsx_wgt_flx(bnd_idx)=xsx_wgt_flx(bnd_idx)/mlt_fct
  enddo ! bnd_idx

  if (dbg_lvl == dbg_crr) then
     ! Examine weights
     write (6,'(5(a,1x))') '#  ','wvl        ','xsx_abs    ','flx TOA    ',wgt_nm(1:ftn_strlen(wgt_nm))
     write (6,'(5(a,1x))') '   ','m          ','m-2 mlc-1  ','W m-2 m-1  ','W m-2 m-1  '
     do bnd_idx=1,bnd_nbr
        write (6,'(i4,1x,4(es10.3,1x))') &
             bnd_idx,wvl_ctr(bnd_idx),abs_xsx_O2O2(bnd_idx),flx_spc_dwn_TOA(bnd_idx),wgt_spc(bnd_idx)
     enddo ! bnd_idx
  endif                     ! endif dbg
  
  ! Normalize flux-weighted absorption cross sections
  do bnd_idx=1,bnd_nbr_CCM
     flx_spc_dwn_TOA_CCM(bnd_idx)=flx_slr_frc_CCM(bnd_idx)*slr_cst_CCM/wvl_dlt_CCM(bnd_idx)
     if (WGT_TRN) then
        if (wgt_spc_CCM(bnd_idx) /= 0.0) then
           abs_xsx_O2O2_CCM(bnd_idx)=xsx_wgt_flx_CCM(bnd_idx)/wgt_spc_CCM(bnd_idx)
        else
           write (6,'(a,a,i2,a)') prg_nm(1:ftn_strlen(prg_nm)), &
                ': WARNING wgt_spc_CCM(',bnd_idx,') = 0.0 Setting cross section equal to zero.'
           abs_xsx_O2O2_CCM(bnd_idx)=0.0
        endif               ! endif not degenerate
     else
        abs_xsx_O2O2_CCM(bnd_idx)=xsx_wgt_flx_CCM(bnd_idx)/flx_spc_dwn_TOA_CCM(bnd_idx)
     endif                  ! endif WGT_TRN
     !     abs_cff_mss_O2O2_CCM(bnd_idx)=abs_xsx_O2O2_CCM(bnd_idx)*Avagadro/mmw_O2O2 ! 20181210 Old method
     abs_cff_mss_O2O2_CCM(bnd_idx)=abs_xsx_O2O2_CCM(bnd_idx)*(dble(Avagadro)/dble(mmw_O2))**2 ! [m5 kg-2]
     !     odal_O2O2_PUNC_O2_PUNP_O2_CCM(bnd_idx)=dble(abs_xsx_O2O2_CCM(bnd_idx))*dble(k_O2_O2) ! NB: Double precision arithmetic is required here
     odal_O2O2_PUNC_O2_PUNP_O2_CCM(bnd_idx)=abs_xsx_O2O2_CCM(bnd_idx)
     double_foo= &
             odal_O2O2_PUNC_O2_PUNP_O2_CCM(bnd_idx)*(dble(Avagadro)/dble(mmw_O2))**2 ! [m5 kg-2]
     odal_O2O2_PUMC_O2_PUMP_O2_CCM(bnd_idx)=double_foo ! [m5 kg-2] Defensive programming
  enddo ! bnd_idx
  
  ! Do the same thing for O2-N2
  if (WGT_TRN) then
     ! Weight high resolution absorption cross-sections by atmospheric transmission of atmosphere with all constituents except H2OH2O
     call wgt_get(fl_wgt,wgt_nm,wvl_grd,bnd_nbr,wgt_spc)         
     do bnd_idx=1,bnd_nbr
        xsx_wgt_flx(bnd_idx)=abs_xsx_O2N2(bnd_idx)*wgt_spc(bnd_idx)
     enddo ! bnd_idx
  else                      ! !WGT_TRN
     ! Weight high resolution absorption cross sections by high resolution TOA solar flux
     do bnd_idx=1,bnd_nbr
        xsx_wgt_flx(bnd_idx)=abs_xsx_O2N2(bnd_idx)*flx_spc_dwn_TOA(bnd_idx)
     enddo ! bnd_idx
  endif                     ! !WGT_TRN
  
  ! Rebin flux-weighted absorption cross sections
  call rbn_vec_CCM(bnd_nbr,wvl_grd,xsx_wgt_flx, &
       bnd_nbr_CCM,wvl_min_CCM,wvl_max_CCM,xsx_wgt_flx_CCM, &
       xtr_typ_LHS,xtr_typ_RHS)
  if (WGT_TRN) then
     ! Rebin cross-section weights
     call rbn_vec_CCM(bnd_nbr,wvl_grd,wgt_spc, &
          bnd_nbr_CCM,wvl_min_CCM,wvl_max_CCM,wgt_spc_CCM, &
          xtr_typ_LHS,xtr_typ_RHS)
  endif                     ! endif WGT_TRN
  
  ! Normalize flux-weighted absorption cross sections
  do bnd_idx=1,bnd_nbr_CCM
     flx_spc_dwn_TOA_CCM(bnd_idx)=flx_slr_frc_CCM(bnd_idx)*slr_cst_CCM/wvl_dlt_CCM(bnd_idx)
     if (WGT_TRN) then
        if (wgt_spc_CCM(bnd_idx) /= 0.0) then
           abs_xsx_O2N2_CCM(bnd_idx)=xsx_wgt_flx_CCM(bnd_idx)/wgt_spc_CCM(bnd_idx)
        else
           write (6,'(a,a,i2,a)') prg_nm(1:ftn_strlen(prg_nm)), &
                ': WARNING wgt_spc_CCM(',bnd_idx,') = 0.0 Setting cross section equal to zero.'
           abs_xsx_O2N2_CCM(bnd_idx)=0.0
        endif               ! endif not degenerate
     else
        abs_xsx_O2N2_CCM(bnd_idx)=xsx_wgt_flx_CCM(bnd_idx)/flx_spc_dwn_TOA_CCM(bnd_idx)
     endif                  ! endif WGT_TRN
     abs_cff_mss_O2N2_CCM(bnd_idx)=abs_xsx_O2N2_CCM(bnd_idx)*Avagadro/mmw_O2N2
     !     odal_O2N2_PUNC_O2_PUNP_N2_CCM(bnd_idx)=dble(abs_xsx_O2N2_CCM(bnd_idx))*dble(k_O2_O2) ! NB: Double precision arithmetic is required here
     odal_O2N2_PUNC_O2_PUNP_N2_CCM(bnd_idx)=dble(abs_xsx_O2N2_CCM(bnd_idx)) ! NB: Double precision arithmetic is required here
     double_foo= &
          odal_O2N2_PUNC_O2_PUNP_N2_CCM(bnd_idx)*dble(Avagadro)*dble(Avagadro)/(dble(mmw_O2)*dble(mmw_N2))
     odal_O2N2_PUMC_O2_PUMP_N2_CCM(bnd_idx)=double_foo ! Defensive programming
  enddo ! bnd_idx
  
  if (dbg_lvl == dbg_old) then
     ! Compare retrieved versus rebinned spectral fluxes
     write (6,'(3(a,1x))') 'idx','flx_spc rtr','flx_spc rbn'
     write (6,'(3(a,1x))') '   ','W m-2 m-1  ','W m-2 m-1  '
     do bnd_idx=1,bnd_nbr_CCM
        write (6,'(i4,1x,2(es10.3,1x))') &
             bnd_idx,flx_spc_dwn_TOA_CCM(bnd_idx),foo_CCM(bnd_idx)
     enddo ! bnd_idx
  endif                     ! endif dbg
  
  ! Add CCM grid to netCDF output file
  call nc_out_CCM_SW(fl_out,bnd_dmn_id)
  
  ! Add O2O2 data on CCM grid to netCDF output file
  rcd=nf90_wrp_open(fl_out,nf90_write,nc_id,sbr_nm)
  
  ! Put output file in define mode
  rcd=nf90_redef(nc_id)
  if (rcd /= nf90_noerr) call nf90_err_exit(rcd,fl_out)
  
  ! Variable definitions
  rcd=nf90_def_var(nc_id,'abs_cff_mss_O2O2_CCM',nf90_float,bnd_dmn_id,abs_cff_mss_O2O2_id)
  rcd=nf90_def_var(nc_id,'abs_xsx_O2O2_CCM',nf90_float,bnd_dmn_id,abs_xsx_O2O2_id)
  rcd=nf90_def_var(nc_id,'flx_slr_frc_CCM',nf90_float,bnd_dmn_id,flx_slr_frc_id)
  rcd=nf90_def_var(nc_id,'flx_spc_dwn_TOA_CCM',nf90_float,bnd_dmn_id,flx_spc_dwn_TOA_id)
  rcd=nf90_def_var(nc_id,'odal_O2N2_PUMC_O2_PUMP_N2_CCM',nf90_float,bnd_dmn_id,odal_O2N2_PUMC_O2_PUMP_N2_id)
  rcd=nf90_def_var(nc_id,'odal_O2O2_PUMC_O2_PUMP_O2_CCM',nf90_float,bnd_dmn_id,odal_O2O2_PUMC_O2_PUMP_O2_id)
  rcd=nf90_def_var(nc_id,'odal_O2O2_PUNC_O2_PUNP_O2_CCM',nf90_double,bnd_dmn_id,odal_O2O2_PUNC_O2_PUNP_O2_id)
  
  ! Add global attributes
  
  ! Add english text descriptions
  rcd=nf90_wrp(nf90_put_att(nc_id,abs_cff_mss_O2O2_id,'long_name','O2-O2 mass absorption coefficient'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,flx_slr_frc_id,'long_name','Fraction of solar flux in band: ' // fl_slr),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,flx_spc_dwn_TOA_id,'long_name','Spectral solar insolation at TOA'),sbr_nm)
  ! Wrap
  rcd=nf90_wrp(nf90_put_att(nc_id,odal_O2O2_PUMC_O2_PUMP_O2_id,'long_name', &
       'O2-O2 layer absorption optical depth (per kg m-3 O2 per kg m-2 O2)'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,odal_O2N2_PUMC_O2_PUMP_N2_id,'long_name', &
       'O2-N2 layer absorption optical depth (per kg m-3 O2 per kg m-2 N2)'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,odal_O2O2_PUNC_O2_PUNP_O2_id,'long_name', &
       'O2-O2 layer absorption optical depth (per mlc m-3 O2 per mlc m-2 O2)'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,abs_xsx_O2O2_id,'long_name', &
       'O2-O2 collision-induced absorption cross section'),sbr_nm)
  
  ! Add units
  rcd=nf90_wrp(nf90_put_att(nc_id,abs_cff_mss_O2O2_id,'units','meter5 kilogram-2'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,abs_xsx_O2O2_id,'units','meter5 molecule-2'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,flx_slr_frc_id,'units','fraction'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,flx_spc_dwn_TOA_id,'units','watt meter-2 meter-1'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,odal_O2N2_PUMC_O2_PUMP_N2_id,'units','meter5 kilogram-2'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,odal_O2O2_PUMC_O2_PUMP_O2_id,'units','meter5 kilogram-2'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,odal_O2O2_PUNC_O2_PUNP_O2_id,'units','meter5 molecule-2'),sbr_nm)
  
  ! All dimensions, variables, and attributes have been defined
  rcd=nf90_wrp(nf90_enddef(nc_id),sbr_nm//': enddef in '//__FILE__) 
  
  ! Write out data
  rcd=nf90_wrp(nf90_put_var(nc_id,abs_cff_mss_O2O2_id,abs_cff_mss_O2O2_CCM),sbr_nm//': pv abs_cff_mss_O2O2 in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,abs_xsx_O2O2_id,abs_xsx_O2O2_CCM),sbr_nm//': pv abs_xsx_O2O2 in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,flx_slr_frc_id,flx_slr_frc_CCM),sbr_nm//': pv flx_slr_frc in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,flx_spc_dwn_TOA_id,flx_spc_dwn_TOA_CCM),sbr_nm//': pv flx_spc_dwn_TOA in '//__FILE__)
  ! Wrap
  rcd=nf90_wrp(nf90_put_var(nc_id,odal_O2O2_PUNC_O2_PUNP_O2_id,odal_O2O2_PUNC_O2_PUNP_O2_CCM), &
       sbr_nm//': pv odal_O2O2_PUNC_O2_PUNP_O2 in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,odal_O2N2_PUMC_O2_PUMP_N2_id,odal_O2N2_PUMC_O2_PUMP_N2_CCM), &
       sbr_nm//': pv odal_O2N2_PUMC_O2_PUMP_N2 in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,odal_O2O2_PUMC_O2_PUMP_O2_id,odal_O2O2_PUMC_O2_PUMP_O2_CCM), &
       sbr_nm//': pv odal_O2O2_PUMC_O2_PUMP_O2 in '//__FILE__)
  
  ! Close output file
  rcd=nf90_wrp_close(nc_id,fl_out,'Wrote results to') ! [fnc] Close file
  
  ! Convert absorption coefficients to CGS and output block for radcsw
  do bnd_idx=1,bnd_nbr_CCM
     odal_O2O2_PUMC_O2_PUMP_O2_CCMCG(bnd_idx)=odal_O2O2_PUMC_O2_PUMP_O2_CCM(bnd_idx)*(100.0**5)/(1000.0**2) ! m5 kg-2 = m3 kg-1 O2 m2 kg-1 O2 --> cm5 g-2 = cm3 g-1 O2 cm2 g-1 O2 
     odal_O2N2_PUMC_O2_PUMP_N2_CCMCG(bnd_idx)=odal_O2N2_PUMC_O2_PUMP_N2_CCM(bnd_idx)*(100.0**5)/(1000.0**2) ! m5 kg-2 = m3 kg-1 O2 m2 kg-1 O2 --> cm5 g-2 = cm3 g-1 O2 cm2 g-1 O2 
  enddo ! bnd_idx
  nbr_dat_per_ln=4
  write (0,'(a)') 'c     O2X physics and forcing documented in Zender (1999), JGR, 104, D20, 24471--24484'
  write (0,'(a,a,a,a)') 'c     Data generated by ',prg_nm(1:ftn_strlen(prg_nm)),' on ',lcl_date_time
  write (0,'(a,a)') 'c     ',prg_ID(1:ftn_strlen(prg_ID))
  write (0,'(a,a)') 'c     Command line: ',cmd_ln(1:ftn_strlen(cmd_ln))
  write (0,'(a)') 'c     O2-O2 absorption optical depth in cm5 g-2 (per g cm-3 O2 per g cm-2 O2)'
  write (0,'(a,f7.3,a)') 'c     Valid conditions: Temperature = ',tpt_std,' K'
  write (0,'(a,a)') 'c     ',src_rfr_sng(1:ftn_strlen(src_rfr_sng))
  write (0,'(a,a)') 'c     ',src_fl_sng(1:ftn_strlen(src_fl_sng))
  write (0,'(a,a)') 'c     ',src_wgt_sng(1:ftn_strlen(src_wgt_sng))
  call dat_prn_f77(bnd_nbr_CCM,odal_O2O2_PUMC_O2_PUMP_O2_CCMCG,nbr_dat_per_ln,'abs_cff_mss_pair_O2O2'//char(0))
  call dat_prn_f90(bnd_nbr_CCM,odal_O2O2_PUMC_O2_PUMP_O2_CCMCG,nbr_dat_per_ln,'abs_cff_mss_pair_O2O2'//char(0))
  
  write (0,'(a,a,a,a)') 'c     Data generated by ',prg_nm(1:ftn_strlen(prg_nm)),' on ',lcl_date_time
  write (0,'(a,a)') 'c     ',prg_ID(1:ftn_strlen(prg_ID))
  write (0,'(a,a)') 'c     Command line: ',cmd_ln(1:ftn_strlen(cmd_ln))
  write (0,'(a)') 'c     O2-N2 absorption optical depth in cm5 g-2 (per g cm-3 O2 per g cm-2 N2)'
  write (0,'(a,f7.3,a)') 'c     Valid conditions: Temperature = ',tpt_std,' K'
  write (0,'(a,f7.3)') 'c     N2 collision partner efficiency = ',N2_cll_prt_fsh 
  write (0,'(a,a)') 'c     ',src_rfr_sng(1:ftn_strlen(src_rfr_sng))
  write (0,'(a,a)') 'c     ',src_fl_sng(1:ftn_strlen(src_fl_sng))
  write (0,'(a,a)') 'c     ',src_wgt_sng(1:ftn_strlen(src_wgt_sng))
  call dat_prn_f77(bnd_nbr_CCM,odal_O2N2_PUMC_O2_PUMP_N2_CCMCG,nbr_dat_per_ln,'abs_cff_mss_pair_O2N2'//char(0))
  call dat_prn_f90(bnd_nbr_CCM,odal_O2N2_PUMC_O2_PUMP_N2_CCMCG,nbr_dat_per_ln,'abs_cff_mss_pair_O2N2'//char(0))
  
  ! De-allocate dynamic variables
  if (allocated(abs_cff_mss_O2O2)) deallocate(abs_cff_mss_O2O2,stat=rcd)
  if (allocated(abs_xsx_O2N2)) deallocate(abs_xsx_O2N2,stat=rcd)
  if (allocated(abs_xsx_O2O2)) deallocate(abs_xsx_O2O2,stat=rcd)
  if (allocated(abs_xsx_O2O2_dadT)) deallocate(abs_xsx_O2O2_dadT,stat=rcd)
  if (allocated(abs_xsx_O2O2_tpt_rfr)) deallocate(abs_xsx_O2O2_tpt_rfr,stat=rcd)
  if (allocated(bnd)) deallocate(bnd,stat=rcd)     ! coordinate variable
  if (allocated(flx_bnd_dwn_TOA)) deallocate(flx_bnd_dwn_TOA,stat=rcd)
  if (allocated(flx_bnd_pht_dwn_TOA)) deallocate(flx_bnd_pht_dwn_TOA,stat=rcd)
  if (allocated(flx_slr_frc)) deallocate(flx_slr_frc,stat=rcd)
  if (allocated(flx_spc_dwn_TOA)) deallocate(flx_spc_dwn_TOA,stat=rcd)
  if (allocated(idx_rfr_air_STP)) deallocate(idx_rfr_air_STP,stat=rcd)
  if (allocated(nrg_pht)) deallocate(nrg_pht,stat=rcd)
  if (allocated(odal_O2N2_PUMC_O2_PUMP_N2)) deallocate(odal_O2N2_PUMC_O2_PUMP_N2,stat=rcd)
  if (allocated(odal_O2N2_PUNC_O2_PUNP_N2)) deallocate(odal_O2N2_PUNC_O2_PUNP_N2,stat=rcd)
  if (allocated(odal_O2O2_PUMC_O2_PUMP_O2)) deallocate(odal_O2O2_PUMC_O2_PUMP_O2,stat=rcd)
  if (allocated(odal_O2O2_PUNC_O2_PUNP_O2)) deallocate(odal_O2O2_PUNC_O2_PUNP_O2,stat=rcd)
  if (allocated(wgt_spc)) deallocate(wgt_spc,stat=rcd)
  if (allocated(wvl)) deallocate(wvl,stat=rcd)     ! coordinate variable
  if (allocated(wvl_ctr)) deallocate(wvl_ctr,stat=rcd)
  if (allocated(wvl_dlt)) deallocate(wvl_dlt,stat=rcd)
  if (allocated(wvl_grd)) deallocate(wvl_grd,stat=rcd)
  if (allocated(wvl_max)) deallocate(wvl_max,stat=rcd)
  if (allocated(wvl_min)) deallocate(wvl_min,stat=rcd)
  if (allocated(wvn)) deallocate(wvn,stat=rcd)     ! coordinate variable
  if (allocated(wvn_ctr)) deallocate(wvn_ctr,stat=rcd)
  if (allocated(wvn_dlt)) deallocate(wvn_dlt,stat=rcd)
  if (allocated(wvn_grd)) deallocate(wvn_grd,stat=rcd)
  if (allocated(wvn_max)) deallocate(wvn_max,stat=rcd)
  if (allocated(wvn_min)) deallocate(wvn_min,stat=rcd)
  if (allocated(xsx_wgt_flx)) deallocate(xsx_wgt_flx,stat=rcd)
  
1000 continue
  
  !call exit(exit_status)
end program O2X

