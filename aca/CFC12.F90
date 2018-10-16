! $Id$

program CFC12
  
  ! Purpose: Convert CFC12 absorption cross section data to netCDF format
  ! CCl2F2 is Dichlorodifluoromethane (R-12) is a colorless gas usually sold under the brand name Freon-12 ... Freon 12. R-12. CFC-12. P-12. Propellant 12. 
  
  ! Compilation:
  ! cd ${HOME}/sw/aca; make OPTS=D CFC12; cd -
  
  ! Use HITRAN data
  ! CFC12 --HTR16
  ! CFC12 --HTR16 --drc_out=${DATA}/aca -o abs_xsx_CFC12_HTR16.nc
  ! CFC12 -i ${DATA}/aca/CCl3F_278.1K-760.0Torr_570.0-6500.0_0.11_N2_50_43.xsc -o ${DATA}/aca/abs_xsx_CFC12.nc # Use HITRAN16 data
  ! CFC12 -S ${DATA}/aca/spc_LaN68.nc # Weight WMO85 data by LaN68

  ! Process input .xsc ASCII data files from HITRAN that look like:
  !             CCl2F2  799.9950 1270.0094 187181  294.2 350.5 1.162E-170.030         CFC-12    air 54
  ! 0.000E+00 0.000E+00 0.000E+00 0.000E+00 0.000E+00 0.000E+00 0.000E+00 0.000E+00 1.319E-21 2.788E-21
  ! 4.001E-21 4.842E-21 5.301E-21 5.451E-21 5.403E-21 5.261E-21 5.105E-21 4.995E-21 4.977E-21 5.092E-21

#if 0
  From http://hitran.org/docs/cross-sections-definitions
  "In the HITRAN FTP site, the data are presented as separate files for each individual molecule. Each portion of the file corresponding to a particular temperature-pressure pair begins with a header (see Table 1) that contains information on the wavenumber (cm−1) range, number of cross-section data in this set, temperature (K), and pressure (Torr). The maximum value of the absorption cross sections (cm2/molecule) and additional information containing the reference to that observation are also presented in each header. The cross sections have been cast into an equal wavenumber interval grid. It should be noted that the initial and final wavenumbers, νmin and νmax, respectively, of each temperature–pressure set for a given wavenumber region are not always identical. They have been taken from the analysis of the observations. The sampling intervals are also not necessarily identical for each temperature–pressure set. The wavenumber interval of the grid is obtained by taking the difference of the initial and final wavenumber and dividing this quantity by the number of points, N, minus one, i.e., Δν=(νmax−νmin)/(N−1).

  This value of N is provided so that a user’s personal program can read the table of cross- sections that follows the header. Note that the use of the features of HITRANonline makes much of this discussion transparent.

  The table below illustrates the format of each header record. Following the header, the cross-section values are arranged in records containing ten values of fields of ten for each cross-section. In other words, each record contains 100 bytes (the trailing bytes on the last line may not be meaningful if N is not a multiple of 10).

  Quantity	Field length	Type	Comment
  Molecule	20	Character	Chemical formula (right-justified)
  Minimum wavenumber, νmin	10	Real	Start of range (cm−1)
  Maximum wavenumber, νmax	10	Real	End of range (cm−1)
  Number of points, N	7	Integer	Number of cross-sections in set
  Temperature, T	7	Real	Temperature (K) of set
  Pressure, P	6	Real	Pressure of set in Torr
  Maximum cross-section value in set, σmax	10	Real	Useful for scaling plots (cm2/molecule)
  Instrument resolution	5	Real	See note
  Common name	15	Character	Familiar name of molecule
  Not currently used	4		Reserved for future use
  Broadener	3	Character	Air, N2, or self-broadened (if left blank)
  Reference	3	Integer	Index pointing to source of data
  Note: Most cross sections have been taken from Fourier transform spectrometer (FTS) measurements. In that case the resolution is given in cm−1. There are some cross-sections taken from grating spectrometer measurements in the UV. In those cases, the resolution is given in milli-Ångströms in the form xxx mÅ, where xxx are up to three digits.

  In the FTP site, for the IR cross-sections, the data on each molecule (chemical compound) are stored in separate files, which are labeled with the chemical symbol followed by an underscore and IRxx.xsc, where xx stands for the HITRAN edition that the data were originally introduced or later updated and the file extension xsc signifies that it is a list of cross-sections. For example, the file with the name C2H6_IR10.xsc contains ethane (C2H6) infrared cross-sections that were obtained in 2010.

  It is to be noted that the files may have many temperature–pressure sets for different spectral regions, as indicated by headers throughout the file. While the temperature–pressure (T,p) sets are reasonably complete for many species for an adequate simulation of atmospheric transmission in the spectral regions where those species are active, for other species an insufficiency of the (T,p) sets may become apparent. It is hoped that future measurements at extended sets of (T,p) combinations may help broaden the coverage in the database."
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
  character(len=*),parameter::sbr_nm='CFC12' ! [sng] Subroutine name
  character(*),parameter::fl_in_HTR16_cold='CCl2F2_216.2K-379.5Torr_800.0-1270.0_00.xsc'
  character(*),parameter::fl_in_HTR16_warm='CCl2F2_294.2K-350.5Torr_800.0-1270.0_00.xsc'
  character(*),parameter::fl_out_HTR16='abs_xsx_CFC12_HTR16.nc'
  character(*),parameter::fl_slr_dfl='spc_Kur95_01wvn.nc'
  character(*),parameter::nlc=char(0) ! [sng] NUL character = ASCII 0 = char(0)
  
  ! Locals with simple initialization and no command-line override
  ! integer::exit_status=0 ! [enm] Program exit status (non-standard Fortran)
  integer::rcd=nf90_noerr ! [rcd] Return success code

  integer,parameter::fl_in_unit=73
  integer,parameter::bnd_nbr_HTR16=187180 ! [nbr] 187181 is # of interfaces, 187180 is # of bands
  integer,parameter::sng_lng_dfl_fl=80 ! [nbr] Default filename string length
  integer,parameter::sng_lng_dfl_stt=200 ! [nbr] Default statement string length
  real,parameter::tpt_cold_HTR16=216.2
  real,parameter::tpt_warm_HTR16=294.2
  real,parameter::mss_val=nf90_fill_float ! Missing value = missing_value and/or _FillValue

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

  ! Local workspace
  character(sng_lng_dfl_fl)::src_fl_sng
  character(sng_lng_dfl_fl)::src_rfr_sng
  
  ! HITRAN XSC format
  character::brd_nm_htrn*3 ! [sng] Broadener (Air, N2, or self-broadened (if left blank))
  character::chr_foo_htrn*4 ! [sng] Reserved for future use
  character::mA_sng_htrn*2 ! [sng] mA
  character::mlc_frm_htrn*20 ! [sng] Molecule chemical formula (right-justified)
  character::mlc_nm_htrn*15 ! [sng] Molecule common name
  integer::bnd_nbr_htrn ! [nbr] Number of wavenumber bins
  integer::rfr_nbr_htrn ! [idx] Reference number (# of bibliographic reference in HITRAN GRH17 reference list)
  real::prs_htrn ! [Pa] Pressure
  real::rsn_nst_htrn ! [m] Instrument resolution
  real::tpt_htrn ! [K] Temperature
  real::wvn_max_htrn ! [cm-1] Wavenumber at end of range
  real::wvn_min_htrn ! [cm-1] Wavenumber at start of range
  real::xsx_max_htrn ! [cm2 mlc-1] Maximum cross-section
  real::wvn_rsn ! [cm-1] Wavenumber resolution
  
  integer::bnd_dmn_id        ! Dimension ID for bands
  integer::grd_dmn_id        ! Dimension ID for grid
  integer::bnd_idx           ! Counting index
  integer::nc_id             ! File handle
  integer::bnd_nbr ! [nbr] Dimension size
  integer::tpt_cold_id
  integer::tpt_warm_id
  integer::tpt_std_id
  integer::Rayleigh_sca_xsx_id
  integer::slr_spc_xtr_typ
  integer::abs_cff_mss_id
  integer::abs_xsx_cold_id
  integer::abs_xsx_dadT_id
  integer::abs_xsx_id
  integer::abs_xsx_tpt_rfr_id
  integer::abs_xsx_warm_id
  integer::bnd_id            ! Coordinate ID
  integer::flx_bnd_dwn_TOA_id
  integer::flx_bnd_pht_dwn_TOA_id
  integer::flx_slr_frc_id
  integer::flx_spc_dwn_TOA_id
  integer::idx_rfr_air_STP_id
  integer::nrg_pht_id
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

  ! netCDF4 
  integer::dfl_lvl=0 ! [enm] Deflate level
  integer::flg_shf=1 ! [flg] Turn on netCDF4 shuffle filter
  integer::flg_dfl=1 ! [flg] Turn on netCDF4 deflate filter
  integer::fl_out_fmt=nco_format_undefined ! [enm] Output file format
  integer::nf90_create_mode=nf90_clobber ! [enm] Mode flag for nf90_create() call
  
  ! Allocatable variables
  real,dimension(:),allocatable::Rayleigh_sca_xsx
  real,dimension(:),allocatable::abs_cff_mss
  real,dimension(:),allocatable::abs_xsx_O2
  real,dimension(:),allocatable::abs_xsx
  real,dimension(:),allocatable::abs_xsx_cold
  real,dimension(:),allocatable::abs_xsx_dadT
  real,dimension(:),allocatable::abs_xsx_tpt_rfr
  real,dimension(:),allocatable::abs_xsx_warm
  real,dimension(:),allocatable::bnd     ! coordinate variable
  real,dimension(:),allocatable::flx_bnd_dwn_TOA
  real,dimension(:),allocatable::flx_bnd_pht_dwn_TOA
  real,dimension(:),allocatable::flx_slr_frc
  real,dimension(:),allocatable::flx_spc_dwn_TOA
  real,dimension(:),allocatable::idx_rfr_air_STP
  real,dimension(:),allocatable::nrg_pht
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
  
  logical::cmd_ln_fl_in=.false.
  logical::cmd_ln_fl_out=.false.
  logical::HTR16=.true.
  
  real::tpt_cold=tpt_cold_HTR16
  real::tpt_std=250.0 ! Temperature at which generic CFC12 cross sections will be archived
  real::tpt_warm=tpt_warm_HTR16
  
  ! Main code

  ! Initialize default values
  dbg_lvl=dbg_off
  fl_slr=fl_slr_dfl
  
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
        else if (opt_sng == 'input' .or. opt_sng == 'fl_CFC12' .or. opt_sng == 'CFC12') then
           call ftn_arg_get(arg_idx,arg_val,fl_in) ! [sng] CFC12 file
        else if (opt_sng == 'HTR16') then
           HTR16=.true.
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
        fl_out_fmt=nf90_format_netcdf4 ! [enm] Output file format
     else if (dsh_key == '-D') then
        call ftn_arg_get(arg_idx,arg_val,dbg_lvl)
     else if (dsh_key == '-i') then
        call ftn_arg_get(arg_idx,arg_val,fl_in)
        cmd_ln_fl_in=.true.
     else if (dsh_key == '-H') then
        HTR16=.true.
     else if (dsh_key == '-o') then
        call ftn_arg_get(arg_idx,arg_val,fl_out)
        cmd_ln_fl_out=.true.
     else if (dsh_key == '-S') then
        call ftn_arg_get(arg_idx,arg_val,fl_slr)
     else if (dsh_key == '-T') then
        call ftn_arg_get(arg_idx,arg_val,tpt_std)
     else if (dsh_key == '-v') then
        write (6,'(a)') CVS_Id
        goto 1000
     else                   ! Option not recognized
        arg_idx=arg_idx-1   ! [idx] Counting index
        call ftn_getarg_err(arg_idx,arg_val) ! [sbr] Error handler for getarg()
     endif if_sgl_dsh       ! endif arg_val
  end do loop_while_options ! end while (arg_idx <= arg_nbr)
  
  ! Compute quantities that depend on command line input
  call ftn_strnul(fl_in)
  call ftn_strnul(fl_out)
  call ftn_strnul(fl_slr)
  call ftn_strcpy(src_fl_sng,'Original data file is ' // fl_in)
  if (HTR16) then
     bnd_nbr=bnd_nbr_HTR16
     tpt_cold=tpt_cold_HTR16
     tpt_warm=tpt_warm_HTR16
     if (.not.cmd_ln_fl_in) fl_in=fl_in_HTR16_cold//nlc
     if (.not.cmd_ln_fl_out) fl_out=fl_out_HTR16//nlc
     call ftn_strcpy(src_rfr_sng,'Data reference is HITRAN (2017) (HTR16)')
  endif ! HTR16

  ! Compute quantities that may depend on command line input
  ! Prepend user-specified path, if any, to input data file names
  if (ftn_strlen(drc_in) > 0) then
     call ftn_drcpfx(drc_in,fl_in) ! [sng] Input file
     call ftn_drcpfx(drc_in,fl_slr) ! [sng] Solar spectrum file
  endif                     ! endif drc_in
  ! Prepend user-specified path, if any, to output data file names
  if (ftn_strlen(drc_out) > 0) call ftn_drcpfx(drc_out,fl_out) ! [sng] Output file
  
  open (fl_in_unit,file=fl_in,status='old',iostat=rcd)

  if (HTR16) then            ! HTR16 data

     ! HITRAN data are in .xsc format described above
     ! First, read-in bnd_nbr in order to allocate memory
     read (fl_in_unit,'(a20,f10.3,f10.3,i7,f7.3,f6.3,e10.3,f3.0,a2,a15,a4,a3,a3)') &
          mlc_frm_htrn,wvn_min_htrn,wvn_max_htrn,bnd_nbr_htrn,tpt_htrn,prs_htrn, &
          xsx_max_htrn,rsn_nst_htrn,mA_sng_htrn,mlc_nm_htrn,chr_foo_htrn,brd_nm_htrn,rfr_nbr_htrn
     tpt_cold=tpt_htrn
     bnd_nbr=bnd_nbr_htrn-1
     write (6,'(a16,1x,a)') 'Read header from',fl_in(1:ftn_strlen(fl_in))

  endif                     ! HTR16 data

  ! Allocate space for dynamic arrays
  allocate(Rayleigh_sca_xsx(bnd_nbr),stat=rcd)
  allocate(abs_cff_mss(bnd_nbr),stat=rcd)
  allocate(abs_xsx_O2(bnd_nbr),stat=rcd)
  allocate(abs_xsx(bnd_nbr),stat=rcd)
  allocate(abs_xsx_cold(bnd_nbr),stat=rcd)
  allocate(abs_xsx_dadT(bnd_nbr),stat=rcd)
  allocate(abs_xsx_tpt_rfr(bnd_nbr),stat=rcd)
  allocate(abs_xsx_warm(bnd_nbr),stat=rcd)
  allocate(bnd(bnd_nbr),stat=rcd)     ! coordinate variable
  allocate(flx_bnd_dwn_TOA(bnd_nbr),stat=rcd)
  allocate(flx_bnd_pht_dwn_TOA(bnd_nbr),stat=rcd)
  allocate(flx_slr_frc(bnd_nbr),stat=rcd)
  allocate(flx_spc_dwn_TOA(bnd_nbr),stat=rcd)
  allocate(idx_rfr_air_STP(bnd_nbr),stat=rcd)
  allocate(nrg_pht(bnd_nbr),stat=rcd)
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

  if (HTR16) then            ! HTR16 data
  
     read (fl_in_unit,*) (abs_xsx_cold(bnd_idx),bnd_idx=1,bnd_nbr)
     close (fl_in_unit)
     write (6,'(a20,1x,a)') 'Read input data from',fl_in(1:ftn_strlen(fl_in))

     ! HITRAN data are in .xsc format described above
     fl_in=fl_in_HTR16_warm//nlc
     if (ftn_strlen(drc_in) > 0) then
        call ftn_drcpfx(drc_in,fl_in) ! [sng] Input file
     endif                     ! endif drc_in
     open (fl_in_unit,file=fl_in,status='old',iostat=rcd)
     read (fl_in_unit,'(a20,f10.3,f10.3,i7,f7.3,f6.3,e10.3,f3.0,a2,a15,a4,a3,a3)') &
          mlc_frm_htrn,wvn_min_htrn,wvn_max_htrn,bnd_nbr_htrn,tpt_htrn,prs_htrn, &
          xsx_max_htrn,rsn_nst_htrn,mA_sng_htrn,mlc_nm_htrn,chr_foo_htrn,brd_nm_htrn,rfr_nbr_htrn
     tpt_warm=tpt_htrn
     bnd_nbr=bnd_nbr_htrn-1
     read (fl_in_unit,*) (abs_xsx_warm(bnd_idx),bnd_idx=1,bnd_nbr)
     close (fl_in_unit)
     write (6,'(a20,1x,a)') 'Read input data from',fl_in(1:ftn_strlen(fl_in))

     ! Sanity check
     if (dbg_lvl >= dbg_fl) then
     
        write (6,'(a20,f10.3,f10.3,i7,f7.3,f6.3,e10.3,f3.0,a2,a15,a4,a3,a3)') &
             mlc_frm_htrn,wvn_min_htrn,wvn_max_htrn,bnd_nbr_htrn,tpt_htrn,prs_htrn, &
             xsx_max_htrn,rsn_nst_htrn,mA_sng_htrn,mlc_nm_htrn,chr_foo_htrn,brd_nm_htrn,rfr_nbr_htrn
        write (6,*) (abs_xsx_cold(bnd_idx),bnd_idx=1,bnd_nbr)
        
     endif                     ! endif dbg

     ! Convert input data to SI units where necessary
     wvn_rsn=(wvn_max_htrn-wvn_min_htrn)/(bnd_nbr-1)
     do bnd_idx=1,bnd_nbr
        wvn_grd(bnd_idx)=wvn_min_htrn+(bnd_idx-1)*wvn_rsn
     enddo
     wvn_grd(bnd_nbr+1)=wvn_min_htrn+bnd_nbr*wvn_rsn
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
        abs_xsx_cold(bnd_idx)=abs_xsx_cold(bnd_idx)*1.0e-4 ! cm2 -> m2
        abs_xsx_warm(bnd_idx)=abs_xsx_warm(bnd_idx)*1.0e-4 ! cm2 -> m2
        abs_xsx_O2(bnd_idx)=mss_val ! cm2 -> m2
        flx_bnd_pht_dwn_TOA(bnd_idx)=mss_val ! #/cm2/s -> #/m2/s
        Rayleigh_sca_xsx(bnd_idx)=mss_val ! cm2 -> m2
     enddo
     
  endif                     ! HTR16 data
  
  ! Define temperature dependance
  ! Temperature dependence is strongest around 
  do bnd_idx=1,bnd_nbr
     abs_xsx_tpt_rfr(bnd_idx)=tpt_std ! Temperature at which generic abs_xsx array will be stored
     abs_xsx_dadT(bnd_idx)= &
          (abs_xsx_warm(bnd_idx)-abs_xsx_cold(bnd_idx))/ &
          (tpt_warm-tpt_cold) 
     abs_xsx(bnd_idx)=abs_xsx_cold(bnd_idx)+ &
          (tpt_std-tpt_cold)*abs_xsx_dadT(bnd_idx)
  enddo                  ! end loop over bnd
  
  ! Get TOA solar spectrum
  slr_spc_xtr_typ=xtr_fll_ngh+xtr_vrb
  call slr_spc_get(fl_slr,wvl_min,wvl_max,bnd_nbr,flx_slr_frc,slr_spc_xtr_typ,slr_spc_xtr_typ)
  
  ! Compute diagnostic variables
  do bnd_idx=1,bnd_nbr
     abs_cff_mss(bnd_idx)=abs_xsx(bnd_idx)*Avagadro/mmw_CFC12
     bnd(bnd_idx)=0.5*(wvl_min(bnd_idx)+wvl_max(bnd_idx))
     wvl(bnd_idx)=bnd(bnd_idx)
     wvl_ctr(bnd_idx)=bnd(bnd_idx)
     wvl_dlt(bnd_idx)=wvl_max(bnd_idx)-wvl_min(bnd_idx)
     wvl_grd(bnd_idx)=wvl_max(bnd_idx) ! NB: because IR cross-sections are on increasing wavenumber grid
     nrg_pht(bnd_idx)=Planck*speed_of_light/wvl(bnd_idx)
     flx_bnd_dwn_TOA(bnd_idx)=flx_slr_frc(bnd_idx)*slr_cst_CCM
     flx_spc_dwn_TOA(bnd_idx)=flx_slr_frc(bnd_idx)*slr_cst_CCM/wvl_dlt(bnd_idx)
     ! NB: flx_bnd_pht_dwn_TOA is actually supplied in the WMO85 data set
     ! What to do, what to store?
     if(flx_bnd_pht_dwn_TOA(bnd_idx) /= mss_val) flx_bnd_pht_dwn_TOA(bnd_idx)=flx_bnd_dwn_TOA(bnd_idx)/nrg_pht(bnd_idx)
  enddo                     ! end loop over bnd
  wvl_grd(bnd_nbr+1)=wvl_max(bnd_nbr)
  
  ! Compute index of refraction through dry air at STP (Len93 p. 155)
  do bnd_idx=1,bnd_nbr
     idx_rfr_air_STP(bnd_idx)= &
          1.0+ &
          1.0e-6*(77.46+0.459/(1.0e12*wvl(bnd_idx)**2))* &
          prs_STP*0.01/tpt_STP
  enddo
  
  ! Sanity check
  if (dbg_lvl >= dbg_fl) then
     write (6,'(5(a,1x))') 'idx','wvl_ctr','abs_xsx','abs_cff_mss','flx_slr_frc'
     do bnd_idx=1,bnd_nbr
        write (6,'(i4,1x,es15.8,1x,es15.8,1x,es15.8,1x,es15.8)')  &
             bnd_idx,wvl_ctr(bnd_idx),abs_xsx(bnd_idx),abs_cff_mss(bnd_idx),flx_slr_frc(bnd_idx)
     enddo                  ! end loop over bnd
  endif                     ! endif dbg

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
  rcd=nf90_wrp(nf90_def_var(nc_id,'Rayleigh_sca_xsx',nf90_float,bnd_dmn_id,Rayleigh_sca_xsx_id),sbr_nm//': dv Rayleigh_sca_xsx')
  rcd=nf90_wrp(nf90_def_var(nc_id,'bnd',nf90_float,bnd_dmn_id,bnd_id),sbr_nm//': dv bnd')
  rcd=nf90_wrp(nf90_def_var(nc_id,'flx_bnd_dwn_TOA',nf90_float,bnd_dmn_id,flx_bnd_dwn_TOA_id),sbr_nm//': dv flx_bnd_dwn_TOA')
  rcd=nf90_wrp(nf90_def_var(nc_id,'flx_slr_frc',nf90_float,bnd_dmn_id,flx_slr_frc_id),sbr_nm//': dv flx_slr_frc')
  rcd=nf90_wrp(nf90_def_var(nc_id,'flx_spc_dwn_TOA',nf90_float,bnd_dmn_id,flx_spc_dwn_TOA_id),sbr_nm//': dv flx_spc_dwn_TOA')
  rcd=nf90_wrp(nf90_def_var(nc_id,'idx_rfr_air_STP',nf90_float,bnd_dmn_id,idx_rfr_air_STP_id),sbr_nm//': dv idx_rfr_air_STP')
  rcd=nf90_wrp(nf90_def_var(nc_id,'nrg_pht',nf90_float,bnd_dmn_id,nrg_pht_id),sbr_nm//': dv nrg_pht')
  rcd=nf90_wrp(nf90_def_var(nc_id,'wvl_ctr',nf90_float,bnd_dmn_id,wvl_ctr_id),sbr_nm//': dv wvl_ctr')
  rcd=nf90_wrp(nf90_def_var(nc_id,'wvl_grd',nf90_float,grd_dmn_id,wvl_grd_id),sbr_nm//': dv wvl_grd')
  rcd=nf90_wrp(nf90_def_var(nc_id,'wvl_dlt',nf90_float,bnd_dmn_id,wvl_dlt_id),sbr_nm//': dv wvl_dlt')
  rcd=nf90_wrp(nf90_def_var(nc_id,'wvl_max',nf90_float,bnd_dmn_id,wvl_max_id),sbr_nm//': dv wvl_max')
  rcd=nf90_wrp(nf90_def_var(nc_id,'wvl_min',nf90_float,bnd_dmn_id,wvl_min_id),sbr_nm//': dv wvl_min')
  rcd=nf90_wrp(nf90_def_var(nc_id,'wvn_ctr',nf90_float,bnd_dmn_id,wvn_ctr_id),sbr_nm//': dv wvn_ctr')
  rcd=nf90_wrp(nf90_def_var(nc_id,'wvn_grd',nf90_float,grd_dmn_id,wvn_grd_id),sbr_nm//': dv wvn_grd')
  rcd=nf90_wrp(nf90_def_var(nc_id,'wvn_dlt',nf90_float,bnd_dmn_id,wvn_dlt_id),sbr_nm//': dv wvn_dlt')
  rcd=nf90_wrp(nf90_def_var(nc_id,'wvn_max',nf90_float,bnd_dmn_id,wvn_max_id),sbr_nm//': dv wvn_max')
  rcd=nf90_wrp(nf90_def_var(nc_id,'wvn_min',nf90_float,bnd_dmn_id,wvn_min_id),sbr_nm//': dv wvn_min')
  ! Wrap
  rcd=nf90_wrp(nf90_def_var(nc_id,'abs_cff_mss',nf90_float,bnd_dmn_id,abs_cff_mss_id), &
       sbr_nm//': dv abs_cff_mss')
  rcd=nf90_wrp(nf90_def_var(nc_id,'abs_xsx',nf90_float,bnd_dmn_id,abs_xsx_id), &
       sbr_nm//': dv abs_xsx')
  rcd=nf90_wrp(nf90_def_var(nc_id,'abs_xsx_cold',nf90_float,bnd_dmn_id,abs_xsx_cold_id), &
       sbr_nm//': dv abs_xsx_cold')
  rcd=nf90_wrp(nf90_def_var(nc_id,'abs_xsx_dadT',nf90_float,bnd_dmn_id,abs_xsx_dadT_id), &
       sbr_nm//': dv abs_xsx_dadT')
  rcd=nf90_wrp(nf90_def_var(nc_id,'abs_xsx_warm',nf90_float,bnd_dmn_id,abs_xsx_warm_id), &
       sbr_nm//': dv abs_xsx_warm')
  rcd=nf90_wrp(nf90_def_var(nc_id,'flx_bnd_pht_dwn_TOA',nf90_float,bnd_dmn_id,flx_bnd_pht_dwn_TOA_id), &
       sbr_nm//': dv flx_bnd_pht_dwn_TOA')
  rcd=nf90_wrp(nf90_def_var(nc_id,'abs_xsx_tpt_rfr',nf90_float,bnd_dmn_id,abs_xsx_tpt_rfr_id), &
       sbr_nm//': dv abs_xsx_tpt_rfr')
  
  ! Add global attributes
  rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'CVS_Id',CVS_Id),sbr_nm//': pa CVS_Id')
  rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'creation_date',lcl_date_time),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'prg_ID',prg_ID(1:ftn_strlen(prg_ID))),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'cmd_ln',cmd_ln(1:ftn_strlen(cmd_ln))),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'src_rfr_sng',src_rfr_sng(1:ftn_strlen(src_rfr_sng))),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'src_fl_sng',src_fl_sng(1:ftn_strlen(src_fl_sng))),sbr_nm)
  
  ! Add english text descriptions
  rcd=nf90_wrp(nf90_put_att(nc_id,tpt_cold_id,'long_name','Temperature of coldest CFC12 measurements'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,tpt_warm_id,'long_name','Temperature of warmest CFC12 measurements'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,tpt_std_id,'long_name','Temperature at which interpolated cross sections are archived'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,Rayleigh_sca_xsx_id,'long_name','Rayleigh scattering cross section'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,abs_cff_mss_id,'long_name','CFC12 mass absorption coefficient'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,abs_xsx_cold_id,'long_name','CFC12 absorption cross section at tpt_cold'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,abs_xsx_dadT_id,'long_name','Absorption cross section temperature-dependence'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,abs_xsx_id,'long_name','CFC12 absorption cross section at tpt_rfr'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,abs_xsx_tpt_rfr_id,'long_name','Valid temperature for absorption cross section'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,abs_xsx_warm_id,'long_name','CFC12 absorption cross section at tpt_warm'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,bnd_id,'long_name','Band center wavelength'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,flx_bnd_dwn_TOA_id,'long_name','Solar Energy flux in band'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,flx_bnd_pht_dwn_TOA_id,'long_name','Photon flux in band'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,flx_slr_frc_id,'long_name','Fraction of solar flux in band: ' // fl_slr),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,flx_spc_dwn_TOA_id,'long_name','Spectral solar insolation at TOA'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,idx_rfr_air_STP_id,'long_name','Index of refraction at band center at STP'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,nrg_pht_id,'long_name','Energy of photon at band center'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,wvl_ctr_id,'long_name','Band center wavelength'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,wvl_dlt_id,'long_name','Bandwidth'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,wvl_grd_id,'long_name','Wavelength grid'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,wvl_max_id,'long_name','Band maximum wavelength'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,wvl_min_id,'long_name','Band minimum wavelength'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,wvn_ctr_id,'long_name','Band center wavenumber'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,wvn_dlt_id,'long_name','Bandwidth'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,wvn_grd_id,'long_name','Wavenumber grid'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,wvn_max_id,'long_name','Band maximum wavenumber'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,wvn_min_id,'long_name','Band minimum wavenumber'),sbr_nm)
  
  ! Add units
  rcd=nf90_wrp(nf90_put_att(nc_id,tpt_cold_id,'units','kelvin'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,tpt_warm_id,'units','kelvin'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,tpt_std_id,'units','kelvin'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,Rayleigh_sca_xsx_id,'units','meter2'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,abs_cff_mss_id,'units','meter2 kilogram-1'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,abs_xsx_cold_id,'units','meter2'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,abs_xsx_dadT_id,'units','meter2 kelvin-1'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,abs_xsx_id,'units','meter2'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,abs_xsx_tpt_rfr_id,'units','kelvin'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,abs_xsx_warm_id,'units','meter2'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,bnd_id,'units','meter'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,flx_bnd_dwn_TOA_id,'units','watt meter-2'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,flx_bnd_pht_dwn_TOA_id,'units','photon meter-2 second-1'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,flx_slr_frc_id,'units','fraction'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,flx_spc_dwn_TOA_id,'units','watt meter-2 meter-1'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,idx_rfr_air_STP_id,'units','fraction'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,nrg_pht_id,'units','joule photon-1'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,wvl_ctr_id,'units','meter'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,wvl_dlt_id,'units','meter'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,wvl_grd_id,'units','meter'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,wvl_max_id,'units','meter'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,wvl_min_id,'units','meter'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,wvn_ctr_id,'units','centimeter-1'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,wvn_dlt_id,'units','centimeter-1'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,wvn_grd_id,'units','centimeter-1'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,wvn_max_id,'units','centimeter-1'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,wvn_min_id,'units','centimeter-1'),sbr_nm)
  
  ! Add _FillValue where necessary
!  rcd=nf90_wrp(nf90_put_att(nc_id,Rayleigh_sca_xsx_id,'_FillValue',mss_val),sbr_nm)
  !rcd=nf90_wrp(nf90_put_att(nc_id,flx_bnd_pht_dwn_TOA_id,'_FillValue',mss_val),sbr_nm)

  ! All dimensions, variables, and attributes have been defined
  rcd=nf90_wrp(nf90_enddef(nc_id),sbr_nm//': enddef in '//__FILE__) 
  
  ! Write data
  rcd=nf90_wrp(nf90_put_var(nc_id,tpt_cold_id,tpt_cold),sbr_nm//': pv tpt_cold in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,tpt_warm_id,tpt_warm),sbr_nm//': pv tpt_warm in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,tpt_std_id,tpt_std),sbr_nm//': pv tpt_std in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,Rayleigh_sca_xsx_id,Rayleigh_sca_xsx),sbr_nm//': pv Rayleigh_sca_xsx in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,abs_cff_mss_id,abs_cff_mss),sbr_nm//': pv abs_cff_mss in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,abs_xsx_cold_id,abs_xsx_cold),sbr_nm//': pv abs_xsx_cold in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,abs_xsx_dadT_id,abs_xsx_dadT),sbr_nm//': pv abs_xsx_dadT in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,abs_xsx_id,abs_xsx),sbr_nm//': pv abs_xsx in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,abs_xsx_warm_id,abs_xsx_warm),sbr_nm//': pv abs_xsx_warm in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,bnd_id,bnd),sbr_nm//': pv bnd in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,flx_bnd_dwn_TOA_id,flx_bnd_dwn_TOA),sbr_nm//': pv flx_bnd_dwn_TOA in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,flx_bnd_pht_dwn_TOA_id,flx_bnd_pht_dwn_TOA),sbr_nm//': pv flx_bnd_pht_dwn_TOA in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,flx_slr_frc_id,flx_slr_frc),sbr_nm//': pv flx_slr_frc in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,flx_spc_dwn_TOA_id,flx_spc_dwn_TOA),sbr_nm//': pv flx_spc_dwn_TOA in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,idx_rfr_air_STP_id,idx_rfr_air_STP),sbr_nm//': pv idx_rfr_air_STP in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,nrg_pht_id,nrg_pht),sbr_nm//': pv nrg_pht in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,wvl_ctr_id,wvl_ctr),sbr_nm//': pv wvl_ctr in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,wvl_dlt_id,wvl_dlt),sbr_nm//': pv wvl_dlt in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,wvl_grd_id,wvl_grd),sbr_nm//': pv wvl_grd in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,wvl_max_id,wvl_max),sbr_nm//': pv wvl_max in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,wvl_min_id,wvl_min),sbr_nm//': pv wvl_min in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,wvn_ctr_id,wvn_ctr),sbr_nm//': pv wvn_ctr in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,wvn_dlt_id,wvn_dlt),sbr_nm//': pv wvn_dlt in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,wvn_grd_id,wvn_grd),sbr_nm//': pv wvn_grd in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,wvn_max_id,wvn_max),sbr_nm//': pv wvn_max in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,wvn_min_id,wvn_min),sbr_nm//': pv wvn_min in '//__FILE__)
  ! Wrap
  rcd=nf90_wrp(nf90_put_var(nc_id,abs_xsx_tpt_rfr_id,abs_xsx_tpt_rfr), &
       sbr_nm//': pv abs_xsx_tpt_rfr in '//__FILE__)
  
  rcd=nf90_wrp_close(nc_id,fl_out,'Wrote results to') ! [fnc] Close file
  
  ! De-allocate dynamic variables
  if (allocated(Rayleigh_sca_xsx)) deallocate(Rayleigh_sca_xsx,stat=rcd)
  if (allocated(abs_cff_mss)) deallocate(abs_cff_mss,stat=rcd)
  if (allocated(abs_xsx_O2)) deallocate(abs_xsx_O2,stat=rcd)
  if (allocated(abs_xsx)) deallocate(abs_xsx,stat=rcd)
  if (allocated(abs_xsx_cold)) deallocate(abs_xsx_cold,stat=rcd)
  if (allocated(abs_xsx_dadT)) deallocate(abs_xsx_dadT,stat=rcd)
  if (allocated(abs_xsx_tpt_rfr)) deallocate(abs_xsx_tpt_rfr,stat=rcd)
  if (allocated(abs_xsx_warm)) deallocate(abs_xsx_warm,stat=rcd)
  if (allocated(bnd)) deallocate(bnd,stat=rcd)     ! coordinate variable
  if (allocated(flx_bnd_dwn_TOA)) deallocate(flx_bnd_dwn_TOA,stat=rcd)
  if (allocated(flx_bnd_pht_dwn_TOA)) deallocate(flx_bnd_pht_dwn_TOA,stat=rcd)
  if (allocated(flx_slr_frc)) deallocate(flx_slr_frc,stat=rcd)
  if (allocated(flx_spc_dwn_TOA)) deallocate(flx_spc_dwn_TOA,stat=rcd)
  if (allocated(idx_rfr_air_STP)) deallocate(idx_rfr_air_STP,stat=rcd)
  if (allocated(nrg_pht)) deallocate(nrg_pht,stat=rcd)
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
  
  ! call exit(exit_status)    ! [enm] Exit with current exit status (non-standard Fortran)

end program CFC12

