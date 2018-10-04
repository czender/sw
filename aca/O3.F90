! $Id$

program O3
  
  ! Purpose: Convert O3 absorption cross section data to netCDF format
  
  ! Compilation:
  ! cd ${HOME}/sw/aca; make -W O3.F90 OPTS=D O3; cd -
  ! cd ${HOME}/sw/aca; make -W O3.F90 O3; cd -
  ! cd ${HOME}/sw/aca; make OPTS=D O3; cd -
  
  ! Use WMO85 or JPL15 data
  ! O3 --JPL15
  ! O3 --WMO85
  ! O3 --JPL15 -o ${DATA}/aca/abs_xsx_O3_JPL15.nc
  ! O3 --WMO85 -o ${DATA}/aca/abs_xsx_O3_WMO85.nc
  ! Use HITRAN16 data
  ! O3 -i ${DATA}/aca/absO3_200.0_0.0_29164.0-40798.0_04.xsc -o ${DATA}/aca/abs_xsx_O3.nc
  ! Weight WMO85 data by LaN68
  ! O3 -S ${DATA}/aca/spc_LaN68.nc
  
  ! Usage:
  ! ncks -H -C -F -v wvl_ctr_CCM,flx_slr_frc_CCM,abs_xsx_O3_CCM ${DATA}/aca/abs_xsx_O3.nc
  ! ncks -H -C -F -d bnd,0.3e-6 -d bnd_CCM,0.3e-6 -v wvl_ctr_CCM,flx_slr_frc_CCM,abs_xsx_O3_CCM,abs_xsx_O3 ${DATA}/aca/abs_xsx_O3.nc
  ! ncwa -a bnd -d bnd,0.295e-6,0.305e-6 ${DATA}/aca/abs_xsx_O3.nc ${DATA}/aca/O3.nc
  ! ncks -H -C -F -d bnd,0.3e-6 -v abs_xsx_O3 ${DATA}/aca/abs_xsx_O3.nc
  ! ncks -H -C -F -v abs_xsx_O3 ${DATA}/aca/O3.nc
  ! ncks -H -C -F -d bnd_CCM,6 -v abs_xsx_O3_CCM ${DATA}/aca/abs_xsx_O3.nc
  
  ! Process input ASCII data files from WMO85 that look like:
  
  ! ==   Reference Solar Irradiance, Rayleigh Scattering, O2 and O3 Xsection==
  ! =============================WMO, 1985, P355-362=======================
  ! spec.  wavelength     irrad.    ray.scat.  O2 Herz.    xs(cm-2) ozone
  ! inte.  range (nm)   (ph/cm2/s)   xs(cm2)    xs(cm2)   T=203K     T=273K
  ! -----------------------------------------------------------------------
  ! 1   175.4 177.0   1.74e+11   6.79e-25   4.61e-24   8.11e-19  8.11e-19 
  ! 2   177.0 178.6   2.10e+11   6.49e-25   5.03e-24   7.99e-19  7.99e-19 
  ! 3   178.6 180.2   2.38e+11   6.20e-25   5.46e-24   7.86e-19  7.89e-19 
  
  ! NB: WMO85 tabular data has a typo reported by Madronich in his TUV model:
  ! "A typo was found in the wmo85 grid (file DATA0/wmo85).  The wavelength
  ! intervals # 7 and # 8 should read
  ! 7   185.185   186.196   3.62E+11   5.16E-25   7.04E-24   6.40E-19   6.40E-19
  ! 8   186.196   188.679   4.73E+11   4.93E-25   7.36E-24   5.88E-19   5.88E-19
  ! instead of 
  ! 7   185.185   186.916   3.62E+11   5.16E-25   7.04E-24   6.40E-19   6.40E-19
  ! 8   186.916   188.679   4.73E+11   4.93E-25   7.36E-24   5.88E-19   5.88E-19
  ! I corrected my copy of the WMO85 data (${DATA}/aca/abs_xsx_WMO85.txt) on 96/12/27
  
  ! or process input .dat ASCII data files from JPL15 that look like:
  ! 20181003 Received from M. Prather, his file XO3_JPL11X.dat
  ! https://jpldataeval.jpl.nasa.gov/pdf/JPL_Publication_15-10.pdf
  ! Table 4A-4. Absorption Cross Sections of O3 at 293-298K & 218K (*1e-20)
  ! 174 (f7.3,2x,f7.3,2f10.3)
  ! range (nm) xs(cm2)  T=293-298K T=218K
  ! -----------------------------------------------------------------------
  ! 175.4    177.0       81.1      81.1
  ! 177.0    178.6       79.9      79.9

  ! or process input .xsc ASCII data files from HITRAN that look like:
  !                  O3    29164.    40798.   5818   200.    0. 1.181E-17250mA          ozone        19
  ! 1.464E-22 2.062E-22 2.295E-22 2.011E-22 1.566E-22 1.346E-22 1.380E-22 1.576E-22 1.756E-22 1.735E-22
  ! 1.627E-22 1.712E-22 1.825E-22 1.423E-22 5.071E-23 0.000E+00 0.000E+00 1.677E-22 4.810E-22 6.631E-22

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
  character(len=*),parameter::CVS_Id='$Id$' ! [sng] CVS Identification
  character(len=*),parameter::sbr_nm='O3' ! [sng] Subroutine name
  character(*),parameter::fl_in_WMO85='abs_xsx_WMO85.txt'
  character(*),parameter::fl_in_JPL15='abs_xsx_O3_JPL15.txt'
  character(*),parameter::fl_out_WMO85='abs_xsx_O3_WMO85.nc'
  character(*),parameter::fl_out_JPL15='abs_xsx_O3_JPL15.nc'
  character(*),parameter::fl_slr_dfl='spc_Kur95_01wvn.nc'
  character(*),parameter::nlc=char(0) ! [sng] NUL character = ASCII 0 = char(0)
  
  integer,parameter::fl_in_unit=73
  integer,parameter::bnd_nbr_JPL15=174
  integer,parameter::bnd_nbr_WMO85=158
  integer,parameter::sng_lng_dfl_fl=80 ! [nbr] Default filename string length
  integer,parameter::sng_lng_dfl_stt=200 ! [nbr] Default statement string length
  real,parameter::tpt_cold_JPL15=218.0
  real,parameter::tpt_warm_JPL15=295.5
  real,parameter::tpt_cold_WMO85=203.0
  real,parameter::tpt_warm_WMO85=273.0
  real,parameter::mss_val=nf90_fill_float ! Missing value = missing_value and/or _FillValue

  ! Input Arguments
  ! Input/Output Arguments
  ! Output Arguments
  ! Local workspace
  character(sng_lng_dfl_fl)::arg_val=nlc      ! [sng] Command line argument value
  character cmd_ln*500      ! [sng] Command line
  character dsh_key*2       ! [sng] Command line dash and switch
  character(sng_lng_dfl_fl)::drc_in=nlc       ! [sng] Input directory
  character(sng_lng_dfl_fl)::drc_out=nlc      ! [sng] Output directory
  character(sng_lng_dfl_fl)::opt_sng=nlc      ! [sng] Option string
  character fl_in*80
  character fl_out*80
  character fl_slr*80
  character lbl*80
  character*26::lcl_date_time ! Time formatted as Day Mth DD HH:MM:SS TZ YYYY
  character prg_ID*200
  character CVS_Date*28
  character CVS_Revision*16
  character src_fl_sng*200
  character src_rfr_sng*200
  
  integer arg_idx           ! [idx] Counting index
  integer arg_nbr           ! [nbr] Number of command line arguments
  integer exit_status       ! [enm] Program exit status
  integer int_foo           ! [nbr] Integer
  integer opt_lng           ! [nbr] Length of option
  integer rcd               ! [rcd] Return success code
  integer idx
  
  logical JPL15
  logical WMO85
  
  integer bnd_dmn_id        ! dimension ID for bands
  integer grd_dmn_id        ! dimension ID for grid
  integer bnd_idx           ! counting index
  integer nc_id             ! file handle
  integer bnd_nbr ! dimension size
  integer tpt_cold_id
  integer tpt_std_id
  integer Rayleigh_sca_xsx_id
  integer abs_cff_mss_O3_id
  integer abs_xsx_O3_cold_id
  integer abs_xsx_O3_dadT_id
  integer abs_xsx_O3_id
  integer abs_xsx_O3_tpt_rfr_id
  integer abs_xsx_O3_warm_id
  integer bnd_id            ! coordinate ID
  integer flx_bnd_dwn_TOA_id
  integer flx_bnd_pht_dwn_TOA_id
  integer flx_slr_frc_id
  integer flx_spc_dwn_TOA_id
  integer idx_rfr_air_STP_id
  integer nrg_pht_id
  integer wvl_ctr_id
  integer wvl_grd_id
  integer wvl_max_id
  integer wvl_min_id
  integer wvl_dlt_id
  ! netCDF4 
  integer::dfl_lvl=0 ! [enm] Deflate level
  integer::flg_shf=1 ! [flg] Turn on netCDF4 shuffle filter
  integer::flg_dfl=1 ! [flg] Turn on netCDF4 deflate filter
  integer::fl_out_fmt=nco_format_undefined ! [enm] Output file format
  integer::nf90_create_mode=nf90_clobber ! [enm] Mode flag for nf90_create() call
  
  real tpt_cold
  real tpt_std
  real tpt_warm
  
  ! Allocatable variables
  real,dimension(:),allocatable::Rayleigh_sca_xsx
  real,dimension(:),allocatable::abs_cff_mss_O3
  real,dimension(:),allocatable::abs_xsx_O2
  real,dimension(:),allocatable::abs_xsx_O3
  real,dimension(:),allocatable::abs_xsx_O3_cold
  real,dimension(:),allocatable::abs_xsx_O3_dadT
  real,dimension(:),allocatable::abs_xsx_O3_tpt_rfr
  real,dimension(:),allocatable::abs_xsx_O3_warm
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
  real,dimension(:),allocatable::xsx_wgt_flx
  
  integer bnd_nbr_CCM
  integer nbr_dat_per_ln
  integer slr_spc_xtr_typ
  integer xtr_typ_LHS
  integer xtr_typ_RHS
  
  real abs_cff_mss_O3_CCM(bnd_nbr_CCM_SW_max)
  real abs_cff_mss_O3_CCM_CGS(bnd_nbr_CCM_SW_max)
  real abs_xsx_O3_CCM(bnd_nbr_CCM_SW_max)
  real bnd_CCM(bnd_nbr_CCM_SW_max)
  real flx_slr_frc_CCM(bnd_nbr_CCM_SW_max)
  real flx_spc_dwn_TOA_CCM(bnd_nbr_CCM_SW_max)
  real foo_CCM(bnd_nbr_CCM_SW_max)
  real wvl_ctr_CCM(bnd_nbr_CCM_SW_max)
  real wvl_max_CCM(bnd_nbr_CCM_SW_max)
  real wvl_min_CCM(bnd_nbr_CCM_SW_max)
  real wvl_dlt_CCM(bnd_nbr_CCM_SW_max)
  real xsx_wgt_flx_CCM(bnd_nbr_CCM_SW_max)

  logical cmd_ln_fl_in
  logical cmd_ln_fl_out

  ! Main code
  
  ! Initialize default values
  drc_in='/data/zender/aca'//nlc ! [sng] Input directory
  drc_out=''                ! [sng] Output directory

  JPL15=.false.
  WMO85=.true.
  cmd_ln_fl_in=.false.
  cmd_ln_fl_out=.false.
  dbg_lvl=dbg_off
  exit_status=0
  fl_slr=fl_slr_dfl
  rcd=nf90_noerr              ! nf90_noerr == 0
  CVS_Date='$Date$'
  CVS_Revision='$Revision$'
  tpt_std=250.0 ! Temperature at which generic O3 cross sections will be archived
  
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
        else if (opt_sng == 'input' .or. opt_sng == 'fl_O3' .or. opt_sng == 'O3') then
           call ftn_arg_get(arg_idx,arg_val,fl_in) ! [sng] Ozone file
        else if (opt_sng == 'JPL15') then
           JPL15=.true.
           WMO85=.false.
        else if (opt_sng == 'WMO85') then
           WMO85=.true.
           JPL15=.false.
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
     else if (dsh_key == '-J') then
        JPL15=.true.
        WMO85=.false.
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
     else if (dsh_key == '-W') then
        WMO85=.true.
        JPL15=.false.
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
  if (WMO85) then
     bnd_nbr=bnd_nbr_WMO85
     tpt_cold=tpt_cold_WMO85
     tpt_warm=tpt_warm_WMO85
     if (.not.cmd_ln_fl_in) fl_in=fl_in_WMO85//nlc
     if (.not.cmd_ln_fl_out) fl_out=fl_out_WMO85//nlc
     call ftn_strcpy(src_rfr_sng,'Data reference is WMO (1985) (WMO85)')
  else if (JPL15) then
     bnd_nbr=bnd_nbr_JPL15
     tpt_cold=tpt_cold_JPL15
     tpt_warm=tpt_warm_JPL15
     if (.not.cmd_ln_fl_in) fl_in=fl_in_JPL15//nlc
     if (.not.cmd_ln_fl_out) fl_out=fl_out_JPL15//nlc
     call ftn_strcpy(src_rfr_sng,'Data reference is JPL (2015) (JPL15)')
  endif

  ! Compute quantities that may depend on command line input
  ! Prepend user-specified path, if any, to input data file names
  if (ftn_strlen(drc_in) > 0) then
     call ftn_drcpfx(drc_in,fl_in) ! [sng] Input file
     call ftn_drcpfx(drc_in,fl_slr) ! [sng] Solar spectrum file
  endif                     ! endif drc_in
  ! Prepend user-specified path, if any, to output data file names
  if (ftn_strlen(drc_out) > 0) call ftn_drcpfx(drc_out,fl_out) ! [sng] Output file

  ! Allocate space for dynamic arrays
  allocate(Rayleigh_sca_xsx(bnd_nbr),stat=rcd)
  allocate(abs_cff_mss_O3(bnd_nbr),stat=rcd)
  allocate(abs_xsx_O2(bnd_nbr),stat=rcd)
  allocate(abs_xsx_O3(bnd_nbr),stat=rcd)
  allocate(abs_xsx_O3_cold(bnd_nbr),stat=rcd)
  allocate(abs_xsx_O3_dadT(bnd_nbr),stat=rcd)
  allocate(abs_xsx_O3_tpt_rfr(bnd_nbr),stat=rcd)
  allocate(abs_xsx_O3_warm(bnd_nbr),stat=rcd)
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
  allocate(xsx_wgt_flx(bnd_nbr),stat=rcd)
  
  open (fl_in_unit,file=fl_in,status='old',iostat=rcd)

  if (WMO85) then            ! WMO85 data
     do idx=1,5
        read (fl_in_unit,'(a80)') lbl
     enddo
     lbl(1:1)=lbl(1:1) ! CEWI
     do bnd_idx=1,bnd_nbr
        read (fl_in_unit,*)  &
             int_foo, &
             wvl_min(bnd_idx), &
             wvl_max(bnd_idx), &
             flx_bnd_pht_dwn_TOA(bnd_idx), &
             Rayleigh_sca_xsx(bnd_idx), &
             abs_xsx_O2(bnd_idx), &
             abs_xsx_O3_cold(bnd_idx), &
             abs_xsx_O3_warm(bnd_idx)
     enddo
     int_foo=int_foo ! CEWI
     
     ! Convert input data to SI units where necessary
     do bnd_idx=1,bnd_nbr
        wvl_min(bnd_idx)=wvl_min(bnd_idx)*1.0e-9 ! nm -> m
        wvl_max(bnd_idx)=wvl_max(bnd_idx)*1.0e-9 ! nm -> m
        flx_bnd_pht_dwn_TOA(bnd_idx)=flx_bnd_pht_dwn_TOA(bnd_idx)*1.0e4 ! #/cm2/s -> #/m2/s
        Rayleigh_sca_xsx(bnd_idx)=Rayleigh_sca_xsx(bnd_idx)*1.0e-4 ! cm2 -> m2
        abs_xsx_O2(bnd_idx)=abs_xsx_O2(bnd_idx)*1.0e-4 ! cm2 -> m2
        abs_xsx_O3_cold(bnd_idx)=abs_xsx_O3_cold(bnd_idx)*1.0e-4 ! cm2 -> m2
        abs_xsx_O3_warm(bnd_idx)=abs_xsx_O3_warm(bnd_idx)*1.0e-4 ! cm2 -> m2
     enddo
     
  else if (JPL15) then            ! JPL15 data

     do idx=1,11
        read (fl_in_unit,'(a80)') lbl
     enddo
     lbl(1:1)=lbl(1:1) ! CEWI
     do bnd_idx=1,bnd_nbr
        read (fl_in_unit,*)  &
             wvl_min(bnd_idx), &
             wvl_max(bnd_idx), &
             abs_xsx_O3_warm(bnd_idx), &
             abs_xsx_O3_cold(bnd_idx)
     enddo
     int_foo=int_foo ! CEWI
     
     ! Convert input data to SI units where necessary
     do bnd_idx=1,bnd_nbr
        wvl_min(bnd_idx)=wvl_min(bnd_idx)*1.0e-9 ! nm -> m
        wvl_max(bnd_idx)=wvl_max(bnd_idx)*1.0e-9 ! nm -> m
        abs_xsx_O3_cold(bnd_idx)=abs_xsx_O3_cold(bnd_idx)*1.0e-20 ! raw -> cm2
        abs_xsx_O3_warm(bnd_idx)=abs_xsx_O3_warm(bnd_idx)*1.0e-20 ! raw -> cm2
        abs_xsx_O3_cold(bnd_idx)=abs_xsx_O3_cold(bnd_idx)*1.0e-4 ! cm2 -> m2
        abs_xsx_O3_warm(bnd_idx)=abs_xsx_O3_warm(bnd_idx)*1.0e-4 ! cm2 -> m2
        abs_xsx_O2(bnd_idx)=mss_val ! cm2 -> m2
        flx_bnd_pht_dwn_TOA(bnd_idx)=mss_val ! #/cm2/s -> #/m2/s
        Rayleigh_sca_xsx(bnd_idx)=mss_val ! cm2 -> m2
     enddo
     
  endif                     ! JPL15 data
  
  close (fl_in_unit)
  write (6,'(a20,1x,a)') 'Read input data from',fl_in(1:ftn_strlen(fl_in))
  
  ! Define temperature dependance
  ! Temperature dependence is strongest around 
  do bnd_idx=1,bnd_nbr
     abs_xsx_O3_tpt_rfr(bnd_idx)=tpt_std ! Temperature at which generic abs_xsx_O3 array will be stored
     abs_xsx_O3_dadT(bnd_idx)= &
          (abs_xsx_O3_warm(bnd_idx)-abs_xsx_O3_cold(bnd_idx))/ &
          (tpt_warm-tpt_cold) 
     abs_xsx_O3(bnd_idx)=abs_xsx_O3_cold(bnd_idx)+ &
          (tpt_std-tpt_cold)*abs_xsx_O3_dadT(bnd_idx)
  enddo                  ! end loop over bnd
  
  ! Get TOA solar spectrum
  slr_spc_xtr_typ=xtr_fll_ngh+xtr_vrb
  call slr_spc_get(fl_slr,wvl_min,wvl_max,bnd_nbr,flx_slr_frc,slr_spc_xtr_typ,slr_spc_xtr_typ)
  
  ! Compute diagnostic variables
  do bnd_idx=1,bnd_nbr
     abs_cff_mss_O3(bnd_idx)=abs_xsx_O3(bnd_idx)*Avagadro/mmw_O3
     bnd(bnd_idx)=0.5*(wvl_min(bnd_idx)+wvl_max(bnd_idx))
     wvl(bnd_idx)=bnd(bnd_idx)
     wvl_ctr(bnd_idx)=bnd(bnd_idx)
     wvl_dlt(bnd_idx)=wvl_max(bnd_idx)-wvl_min(bnd_idx)
     wvl_grd(bnd_idx)=wvl_min(bnd_idx)
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
     write (6,'(5(a,1x))') 'idx','wvl_ctr','abs_xsx_O3','abs_cff_mss_O3','flx_slr_frc'
     do bnd_idx=1,bnd_nbr
        write (6,'(i4,1x,es15.8,1x,es15.8,1x,es15.8,1x,es15.8)')  &
             bnd_idx,wvl_ctr(bnd_idx),abs_xsx_O3(bnd_idx),abs_cff_mss_O3(bnd_idx),flx_slr_frc(bnd_idx)
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
  rcd=nf90_wrp(nf90_def_var(nc_id,'tpt_std',nf90_float,tpt_std_id),sbr_nm//': dv tpt_std')
  rcd=nf90_wrp(nf90_def_var(nc_id,'Rayleigh_sca_xsx',nf90_float,bnd_dmn_id,Rayleigh_sca_xsx_id),sbr_nm//': dv Rayleigh_sca_xsx')
  rcd=nf90_wrp(nf90_def_var(nc_id,'abs_cff_mss_O3',nf90_float,bnd_dmn_id,abs_cff_mss_O3_id),sbr_nm//': dv abs_cff_mss_O3')
  rcd=nf90_wrp(nf90_def_var(nc_id,'abs_xsx_O3',nf90_float,bnd_dmn_id,abs_xsx_O3_id),sbr_nm//': dv abs_xsx_O3')
  rcd=nf90_wrp(nf90_def_var(nc_id,'abs_xsx_O3_cold',nf90_float,bnd_dmn_id,abs_xsx_O3_cold_id),sbr_nm//': dv abs_xsx_O3_cold')
  rcd=nf90_wrp(nf90_def_var(nc_id,'abs_xsx_O3_dadT',nf90_float,bnd_dmn_id,abs_xsx_O3_dadT_id),sbr_nm//': dv abs_xsx_O3_dadT')
  rcd=nf90_wrp(nf90_def_var(nc_id,'abs_xsx_O3_warm',nf90_float,bnd_dmn_id,abs_xsx_O3_warm_id),sbr_nm//': dv abs_xsx_O3_warm')
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
  ! Wrap
  rcd=nf90_wrp(nf90_def_var(nc_id,'flx_bnd_pht_dwn_TOA',nf90_float,bnd_dmn_id,flx_bnd_pht_dwn_TOA_id), &
       sbr_nm//': dv flx_bnd_pht_dwn_TOA')
  rcd=nf90_wrp(nf90_def_var(nc_id,'abs_xsx_O3_tpt_rfr',nf90_float,bnd_dmn_id,abs_xsx_O3_tpt_rfr_id), &
       sbr_nm//': dv abs_xsx_O3_tpt_rfr')
  
  ! Add global attributes
  rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'CVS_Id',CVS_Id),sbr_nm//': pa CVS_Id in '//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'creation_date',lcl_date_time),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'prg_ID',prg_ID(1:ftn_strlen(prg_ID))),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'cmd_ln',cmd_ln(1:ftn_strlen(cmd_ln))),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'src_rfr_sng',src_rfr_sng(1:ftn_strlen(src_rfr_sng))),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,nf90_global,'src_fl_sng',src_fl_sng(1:ftn_strlen(src_fl_sng))),sbr_nm)
  
  ! Add english text descriptions
  rcd=nf90_wrp(nf90_put_att(nc_id,tpt_cold_id,'long_name','Temperature of coldest O3 measurements'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,tpt_std_id,'long_name','Temperature at which interpolated O3 cross sections are archived'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,Rayleigh_sca_xsx_id,'long_name','Rayleigh scattering cross section'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,abs_cff_mss_O3_id,'long_name','Ozone mass absorption coefficient'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,abs_xsx_O3_cold_id,'long_name','Ozone absorption cross section at 203 K'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,abs_xsx_O3_dadT_id,'long_name','Slope of absorption cross section temperature dependence'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,abs_xsx_O3_id,'long_name','Ozone absorption cross section at tpt_rfr'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,abs_xsx_O3_tpt_rfr_id,'long_name','Valid temperature for absorption cross section'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,abs_xsx_O3_warm_id,'long_name','Ozone absorption cross section at 273 K'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,bnd_id,'long_name','Band center wavelength'),sbr_nm)
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
  
  ! Add units
  rcd=nf90_wrp(nf90_put_att(nc_id,tpt_cold_id,'units','kelvin'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,tpt_std_id,'units','kelvin'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,Rayleigh_sca_xsx_id,'units','meter2'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,abs_cff_mss_O3_id,'units','meter2 kilogram-1'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,abs_xsx_O3_cold_id,'units','meter2'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,abs_xsx_O3_dadT_id,'units','meter2 kelvin-1'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,abs_xsx_O3_id,'units','meter2'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,abs_xsx_O3_tpt_rfr_id,'units','kelvin'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,abs_xsx_O3_warm_id,'units','meter2'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,bnd_id,'units','meter'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,flx_bnd_dwn_TOA_id,'units','watt meter-2'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,flx_bnd_pht_dwn_TOA_id,'units','photon meter-2 second-1'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,flx_slr_frc_id,'units','fraction'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,flx_spc_dwn_TOA_id,'units','watt meter-2 meter-1'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,idx_rfr_air_STP_id,'units','fraction'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,nrg_pht_id,'units','joule photon-1'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,wvl_ctr_id,'units','meter'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,wvl_grd_id,'units','meter'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,wvl_max_id,'units','meter'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,wvl_min_id,'units','meter'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,wvl_dlt_id,'units','meter'),sbr_nm)
  
  ! Add _FillValue where necessary
!  rcd=nf90_wrp(nf90_put_att(nc_id,Rayleigh_sca_xsx_id,'_FillValue',mss_val),sbr_nm)
  !rcd=nf90_wrp(nf90_put_att(nc_id,flx_bnd_pht_dwn_TOA_id,'_FillValue',mss_val),sbr_nm)

  ! All dimensions, variables, and attributes have been defined
  rcd=nf90_wrp(nf90_enddef(nc_id),sbr_nm//': enddef in '//__FILE__) 
  
  ! Write data
  rcd=nf90_wrp(nf90_put_var(nc_id,tpt_cold_id,tpt_cold),sbr_nm//': pv tpt_cold in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,tpt_std_id,tpt_std),sbr_nm//': pv tpt_std in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,Rayleigh_sca_xsx_id,Rayleigh_sca_xsx),sbr_nm//': pv Rayleigh_sca_xsx in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,abs_cff_mss_O3_id,abs_cff_mss_O3),sbr_nm//': pv abs_cff_mss_O3 in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,abs_xsx_O3_cold_id,abs_xsx_O3_cold),sbr_nm//': pv abs_xsx_O3_cold in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,abs_xsx_O3_dadT_id,abs_xsx_O3_dadT),sbr_nm//': pv abs_xsx_O3_dadT in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,abs_xsx_O3_id,abs_xsx_O3),sbr_nm//': pv abs_xsx_O3 in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,abs_xsx_O3_tpt_rfr_id,abs_xsx_O3_tpt_rfr),sbr_nm//': pv abs_xsx_O3_tpt_rfr in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,abs_xsx_O3_warm_id,abs_xsx_O3_warm),sbr_nm//': pv abs_xsx_O3_warm in '//__FILE__)
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
  
  rcd=nf90_wrp_close(nc_id,fl_out,'Wrote results to') ! [fnc] Close file
  
  ! Get CCM wavelength grid
  call wvl_grd_CCM_SW_mk(bnd_nbr_CCM,bnd_CCM,wvl_min_CCM,wvl_max_CCM,wvl_ctr_CCM,wvl_dlt_CCM)
  
  ! Get TOA solar spectrum
  call slr_spc_get_CCM(fl_slr,wvl_min_CCM,wvl_max_CCM,bnd_nbr_CCM,flx_slr_frc_CCM,slr_spc_xtr_typ,slr_spc_xtr_typ)
  
  ! Weight high resolution absorption cross sections by high resolution TOA solar flux
  do idx=1,bnd_nbr
     xsx_wgt_flx(idx)=abs_xsx_O3(idx)*flx_spc_dwn_TOA(idx)
  enddo                     ! end loop over O3 bnd
  
  ! Initialize default values
  xtr_typ_LHS=xtr_prt_wgt+xtr_fll_nil+xtr_vrb
  xtr_typ_RHS=xtr_prt_wgt+xtr_fll_nil+xtr_vrb
  
  ! Rebin solar flux
  call rbn_vec_CCM(bnd_nbr,wvl_grd,flx_spc_dwn_TOA, &
       bnd_nbr_CCM,wvl_min_CCM,wvl_max_CCM,foo_CCM, &
       xtr_typ_LHS,xtr_typ_RHS)
  ! Rebin flux-weighted absorption cross sections
  call rbn_vec_CCM(bnd_nbr,wvl_grd,xsx_wgt_flx,  &
       bnd_nbr_CCM,wvl_min_CCM,wvl_max_CCM,xsx_wgt_flx_CCM, &
       xtr_typ_LHS,xtr_typ_RHS)
  
  ! Normalize flux-weighted absorption cross section by solar flux
  do idx=1,bnd_nbr_CCM
     flx_spc_dwn_TOA_CCM(idx)=flx_slr_frc_CCM(idx)*slr_cst_CCM/wvl_dlt_CCM(idx)
     abs_xsx_O3_CCM(idx)=xsx_wgt_flx_CCM(idx)/flx_spc_dwn_TOA_CCM(idx)
     abs_cff_mss_O3_CCM(idx)=abs_xsx_O3_CCM(idx)*Avagadro/mmw_O3
  enddo                     ! end loop over CCM bnd
  
  if (dbg_lvl == dbg_crr) then
     ! Compare retrieved versus rebinned spectral fluxes
     write (6,'(3(a,1x))') 'idx','flx_spc rtr','flx_spc rbn'
     write (6,'(3(a,1x))') '   ','W m-2 m-1  ','W m-2 m-1  '
     do idx=1,bnd_nbr_CCM
        write (6,'(i4,1x,2(es10.3,1x))') &
             idx,flx_spc_dwn_TOA_CCM(idx),foo_CCM(idx)
     enddo                  ! end loop over CCM bnd
  endif                     ! endif dbg
  
  ! Add CCM grid to netCDF output file
  call nc_out_CCM_SW(fl_out,bnd_dmn_id)
  
  ! Add O3 data on CCM grid to netCDF output file
  rcd=nf90_wrp_open(fl_out,nf90_write,nc_id)
  if (rcd /= nf90_noerr) call nf90_err_exit(rcd,fl_out)
  
  ! Put output file in define mode
  rcd=nf90_redef(nc_id)
  if (rcd /= nf90_noerr) call nf90_err_exit(rcd,fl_out)
  
  ! Variable definitions.
  rcd=nf90_wrp(nf90_def_var(nc_id,'abs_xsx_O3_CCM',nf90_float,bnd_dmn_id,abs_xsx_O3_id),sbr_nm//': dv abs_xsx_O3')
  rcd=nf90_wrp(nf90_def_var(nc_id,'abs_cff_mss_O3_CCM',nf90_float,bnd_dmn_id,abs_cff_mss_O3_id),sbr_nm//': dv abs_cff_mss_O3')
  rcd=nf90_wrp(nf90_def_var(nc_id,'flx_slr_frc_CCM',nf90_float,bnd_dmn_id,flx_slr_frc_id),sbr_nm//': dv flx_slr_frc')
  rcd=nf90_wrp(nf90_def_var(nc_id,'flx_spc_dwn_TOA_CCM',nf90_float,bnd_dmn_id,flx_spc_dwn_TOA_id),sbr_nm//': dv flx_spc_dwn_TOA')
  
  ! Add global attributes
  
  ! Add english text descriptions
  rcd=nf90_wrp(nf90_put_att(nc_id,abs_cff_mss_O3_id,'long_name','Ozone mass absorption coefficient'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,abs_xsx_O3_id,'long_name','Ozone continuum absorption cross section'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,flx_slr_frc_id,'long_name','Fraction of solar flux in band: ' // fl_slr),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,flx_spc_dwn_TOA_id,'long_name','Spectral solar insolation at TOA'),sbr_nm)
  
  ! Add units
  rcd=nf90_wrp(nf90_put_att(nc_id,abs_cff_mss_O3_id,'units','meter2 kilogram-1'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,abs_xsx_O3_id,'units','meter2'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,flx_slr_frc_id,'units','fraction'),sbr_nm)
  rcd=nf90_wrp(nf90_put_att(nc_id,flx_spc_dwn_TOA_id,'units','watt meter-2 meter-1'),sbr_nm)
  
  ! Now that all dimensions, variables, and attributes have been defined, make call to end define mode
  rcd=nf90_wrp(nf90_enddef(nc_id),sbr_nm//': enddef in '//__FILE__)
  
  ! Write out data
  rcd=nf90_wrp(nf90_put_var(nc_id,abs_cff_mss_O3_id,abs_cff_mss_O3_CCM),sbr_nm//': pv abs_cff_mss_O3 in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,abs_xsx_O3_id,abs_xsx_O3_CCM),sbr_nm//': pv abs_xsx_O3 in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,flx_slr_frc_id,flx_slr_frc_CCM),sbr_nm//': pv flx_slr_frc in '//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,flx_spc_dwn_TOA_id,flx_spc_dwn_TOA_CCM),sbr_nm//': pv flx_spc_dwn_TOA in '//__FILE__)
  
  ! Close output file
  rcd=nf90_wrp_close(nc_id,fl_out,'Wrote O3 data on CCM grid to')
  
  ! Convert absorption coefficients to cm2/gm and output block for radcsw
  do idx=1,bnd_nbr_CCM
     abs_cff_mss_O3_CCM_CGS(idx)=abs_cff_mss_O3_CCM(idx)*10. !  m2 kg-1 --> cm2 g-1 
  enddo                     ! end loop over CCM bnd
  nbr_dat_per_ln=4
  write (0,'(a,a,a,a)') 'c     Data generated by ',prg_nm(1:ftn_strlen(prg_nm)),' on ',lcl_date_time
  write (0,'(a,a)') 'c     ',prg_ID(1:ftn_strlen(prg_ID))
  write (0,'(a,a)') 'c     Command line: ',cmd_ln(1:ftn_strlen(cmd_ln))
  write (0,'(a,f7.3,a)') 'c     O3 mass absorption coefficients in cm2 g-1 for T = ',tpt_std,' K'
  write (0,'(a,a)') 'c     ',src_rfr_sng(1:ftn_strlen(src_rfr_sng))
  write (0,'(a,a)') 'c     ',src_fl_sng(1:ftn_strlen(src_fl_sng))
  call dat_prn_f77(bnd_nbr_CCM,abs_cff_mss_O3_CCM_CGS,nbr_dat_per_ln,'abO3'//char(0))
  call dat_prn_f90(bnd_nbr_CCM,abs_cff_mss_O3_CCM_CGS,nbr_dat_per_ln,'abO3'//char(0))
  
  ! De-allocate dynamic variables
  if (allocated(Rayleigh_sca_xsx)) deallocate(Rayleigh_sca_xsx,stat=rcd)
  if (allocated(abs_cff_mss_O3)) deallocate(abs_cff_mss_O3,stat=rcd)
  if (allocated(abs_xsx_O2)) deallocate(abs_xsx_O2,stat=rcd)
  if (allocated(abs_xsx_O3)) deallocate(abs_xsx_O3,stat=rcd)
  if (allocated(abs_xsx_O3_cold)) deallocate(abs_xsx_O3_cold,stat=rcd)
  if (allocated(abs_xsx_O3_dadT)) deallocate(abs_xsx_O3_dadT,stat=rcd)
  if (allocated(abs_xsx_O3_tpt_rfr)) deallocate(abs_xsx_O3_tpt_rfr,stat=rcd)
  if (allocated(abs_xsx_O3_warm)) deallocate(abs_xsx_O3_warm,stat=rcd)
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
  if (allocated(xsx_wgt_flx)) deallocate(xsx_wgt_flx,stat=rcd)
  
1000 continue
  
  call exit(exit_status)
end program O3

