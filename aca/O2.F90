! $Id$

program O2
  
  ! Purpose: Convert O2 absorption cross section data to netCDF format
  
  ! Compilation:
  ! cd ${HOME}/sw/aca; make -W O2.F OPTS=D O2; cd -
  ! cd ${HOME}/sw/aca; make -W O2.F O2; cd -
  ! cd ${HOME}/sw/aca; make OPTS=D O2; cd -
  
  ! Usage:
  ! ncks -H -C -F -d bnd,0.3e-6 -v abs_xsx_O2 ${DATA}/aca/abs_xsx_O2.nc
  
  ! Use WMO85 or JPL15 data
  ! O2 -i ${DATA}/aca/abs_xsx_WMO85.txt -o ${DATA}/aca/abs_xsx_O2.nc
  ! Use HITRAN16 data
  ! O2 -i ${DATA}/aca/absO2_200.0_0.0_29164.0-40798.0_04.xsc -o ${DATA}/aca/abs_xsx_O2.nc
  ! Weight WMO85 data by LaN68
  ! O2 -S ${DATA}/aca/spc_LaN68.nc
  
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

  ! or process input .xsc ASCII data files from HITRAN that look like:

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
  character(len=*),parameter::CVS_Id="$Id$" ! [sng] CVS Identification
  character(len=*),parameter::sbr_nm='O2' ! [sng] Subroutine name
  
  ! Arrays are allocatable, but die if size exceeds corresponding *_nbr_max
  integer fl_in_unit
  integer bnd_nbr_O2_HC_JPL15
  integer bnd_nbr_O2_HC_WMO85
  parameter(fl_in_unit=73, &
       bnd_nbr_O2_HC_WMO85=158, &
       bnd_nbr_O2_HC_JPL15=58) &
  ! Input Arguments
  ! Input/Output Arguments
  ! Output Arguments
  ! Local workspace
  character argv*80
  character cmd_ln*200
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
  
  integer arg
  integer arg_nbr
  integer int_foo
  integer exit_status       ! program exit status
  integer idx
  integer rcd               ! return success code
  
  logical JPL15
  logical WMO85
  
  integer bnd_dim_id        ! dimension ID for bands
  integer grd_dim_id        ! dimension ID for grid
  integer bnd_idx           ! counting index
  integer nc_id             ! file handle
  integer bnd_nbr_O2_HC
  integer bnd_nbr_O3_HHCWC
  integer bnd_nbr ! dimension size
  integer tpt_cold_id
  integer tpt_std_id
  integer Rayleigh_sca_xsx_id
  integer abs_cff_mss_O3_id
  integer abs_xsx_O2_id
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
  
  ! Allocatable variables
  real,dimension(:),allocatable::Rayleigh_sca_xsx
  real,dimension(:),allocatable::abs_cff_mss_O3
  real,dimension(:),allocatable::abs_xsx_O2
  real,dimension(:),allocatable::abs_xsx_O3
  real,dimension(:),allocatable::abs_xsx_O3_cold
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
  
  ! Main code
  
  ! Initialize default values
  JPL15=.false.
  WMO85=.true.
  dbg_lvl=dbg_off
  exit_status=0
  fl_in='/data/zender/aca/abs_xsx_WMO85.txt'
  fl_out='/data/zender/aca/abs_xsx_O2.nc'
  fl_slr='/data/zender/aca/spc_Kur95_01wvn.nc'
  rcd=nf90_noerr              ! nf90_noerr == 0
  CVS_Date='$Date$'
  CVS_Revision='$Revision$'
  
  ! Retrieve command line arguments
  call date_time_get(lcl_date_time)
  call ftn_cmd_ln_sng(cmd_ln)
  call ftn_prg_ID_mk(CVS_Id,CVS_Revision,CVS_Date,prg_ID)
  write (6,'(a)') prg_ID(1:ftn_strlen(prg_ID))
  arg_nbr=command_argument_count()
  do arg=1,arg_nbr
     call getarg(arg,argv)
     if (argv(1:2) == '-D') then
        call getarg(arg+1,argv)
        read (argv,'(i4)') dbg_lvl
     endif
     if (argv(1:2) == '-i') then
        call getarg(arg+1,argv)
        read (argv,'(a)') fl_in
     endif
     if (argv(1:2) == '-o') then
        call getarg(arg+1,argv)
        read (argv,'(a)') fl_out
     endif
     if (argv(1:2) == '-S') then
        call getarg(arg+1,argv)
        read (argv,'(a)') fl_slr
     endif
     if (argv(1:2) == '-T') then
        call getarg(arg+1,argv)
        read (argv,'(f8.3)') tpt_std
     endif
     if (argv(1:2) == '-v') then
        write (6,'(a)') CVS_Id
        goto 1000
     endif
     if (argv(1:2) == '-W') then
        WMO85=.true.
     endif
  end do
  
  ! Compute quantities that depend on command line input
  call ftn_strnul(fl_in)
  call ftn_strnul(fl_out)
  call ftn_strnul(fl_slr)
  call ftn_strcpy(src_fl_sng,'Original data file is ' // fl_in)
  if (WMO85) then
     bnd_nbr=bnd_nbr_O2_HC_WMO85
     call ftn_strcpy(src_rfr_sng,'Data reference is WMO (1985) (WMO85)')
  endif
  if (JPL15) then
     bnd_nbr=bnd_nbr_O2_HC_JPL15
     call ftn_strcpy(src_rfr_sng,'Data reference is JPL (2015) (JPL15)')
  endif

  ! Allocate space for dynamic arrays
  allocate(Rayleigh_sca_xsx(bnd_nbr),stat=rcd)
  allocate(abs_cff_mss_O3(bnd_nbr),stat=rcd)
  allocate(abs_xsx_O2(bnd_nbr),stat=rcd)
  allocate(abs_xsx_O3(bnd_nbr),stat=rcd)
  allocate(abs_xsx_O3_cold(bnd_nbr),stat=rcd)
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
     
  endif                     ! WMO85 data

  if (JPL15) then            ! JPL15 data
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
     
  endif                     ! JPL15 data
  
  close (fl_in_unit)
  write (6,'(a20,1x,a)') 'Read input data from',fl_in(1:ftn_strlen(fl_in))
  
  ! Get TOA solar spectrum
  slr_spc_xtr_typ=xtr_fll_ngh+xtr_vrb
  call slr_spc_get(fl_slr,wvl_min,wvl_max,bnd_nbr,flx_slr_frc,slr_spc_xtr_typ,slr_spc_xtr_typ)
  
  ! Compute diagnostic variables
  do bnd_idx=1,bnd_nbr
     abs_cff_mss_O2(bnd_idx)=abs_xsx_O2(bnd_idx)*Avagadro/mmw_O2
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
     flx_bnd_pht_dwn_TOA(bnd_idx)=flx_bnd_dwn_TOA(bnd_idx)/nrg_pht(bnd_idx)
  enddo                     ! end loop over bnd
  wvl_grd(bnd_nbr+1)=wvl_max(bnd_nbr)
  
  ! Compute index of refraction through dry air at STP (Len93 p. 155)
  do bnd_idx=1,bnd_nbr
     idx_rfr_air_STP(bnd_idx)= &
          1.0+ &
          1.0e-6*(77.46+0.459/(1.0e12*wvl(bnd_idx)**2))* &
          prs_STP*0.01/tpt_STP
  enddo
  
  ! Begin netCDF output routines
  rcd=rcd+nf90_create(fl_out,nf90_clobber,nc_id)
  if (rcd /= nf90_noerr) call nf90_err_exit(rcd,fl_out)
  
  ! Define dimension IDs
  rcd=rcd+nf90_def_dim(nc_id,'bnd',bnd_nbr,bnd_dim_id)
  if (rcd /= nf90_noerr) call nf90_err_exit(rcd,fl_out)
  rcd=rcd+nf90_def_dim(nc_id,'grd',bnd_nbr+1,grd_dim_id)
  if (rcd /= nf90_noerr) call nf90_err_exit(rcd,fl_out)
  
  ! Variable definitions
  rcd=nf90_wrp(nf90_def_var(nc_id,'Rayleigh_sca_xsx',nf90_float,bnd_dim_id,Rayleigh_sca_xsx_id),sbr_nm//": dv Rayleigh_sca_xsx")
  rcd=nf90_wrp(nf90_def_var(nc_id,'abs_cff_mss_O2',nf90_float,bnd_dim_id,abs_cff_mss_O2_id),sbr_nm//": dv abs_cff_mss_O2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'abs_xsx_O2',nf90_float,bnd_dim_id,abs_xsx_O2_id),sbr_nm//": dv abs_xsx_O2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'bnd',nf90_float,bnd_dim_id,bnd_id),sbr_nm//": dv bnd")
  rcd=nf90_wrp(nf90_def_var(nc_id,'flx_bnd_dwn_TOA',nf90_float,bnd_dim_id,flx_bnd_dwn_TOA_id),sbr_nm//": dv flx_bnd_dwn_TOA")
  rcd=nf90_wrp(nf90_def_var(nc_id,'flx_slr_frc',nf90_float,bnd_dim_id,flx_slr_frc_id),sbr_nm//": dv flx_slr_frc")
  rcd=nf90_wrp(nf90_def_var(nc_id,'flx_spc_dwn_TOA',nf90_float,bnd_dim_id,flx_spc_dwn_TOA_id),sbr_nm//": dv flx_spc_dwn_TOA")
  rcd=nf90_wrp(nf90_def_var(nc_id,'idx_rfr_air_STP',nf90_float,bnd_dim_id,idx_rfr_air_STP_id),sbr_nm//": dv idx_rfr_air_STP")
  rcd=nf90_wrp(nf90_def_var(nc_id,'nrg_pht',nf90_float,bnd_dim_id,nrg_pht_id),sbr_nm//": dv nrg_pht")
  rcd=nf90_wrp(nf90_def_var(nc_id,'wvl_ctr',nf90_float,bnd_dim_id,wvl_ctr_id),sbr_nm//": dv wvl_ctr")
  rcd=nf90_wrp(nf90_def_var(nc_id,'wvl_grd',nf90_float,grd_dim_id,wvl_grd_id),sbr_nm//": dv wvl_grd")
  rcd=nf90_wrp(nf90_def_var(nc_id,'wvl_max',nf90_float,bnd_dim_id,wvl_max_id),sbr_nm//": dv wvl_max")
  rcd=nf90_wrp(nf90_def_var(nc_id,'wvl_min',nf90_float,bnd_dim_id,wvl_min_id),sbr_nm//": dv wvl_min")
  rcd=nf90_wrp(nf90_def_var(nc_id,'wvl_dlt',nf90_float,bnd_dim_id,wvl_dlt_id),sbr_nm//": dv wvl_dlt")
  ! Wrap
  rcd=nf90_wrp(nf90_def_var(nc_id,'flx_bnd_pht_dwn_TOA',nf90_float,bnd_dim_id,flx_bnd_pht_dwn_TOA_id), &
       sbr_nm//": dv flx_bnd_pht_dwn_TOA")
  
  ! Add global attributes
  rcd=rcd+nf90_put_att(nc_id,nf90_global,'CVS_Id',CVS_Id)
  rcd=rcd+nf90_put_att(nc_id,nf90_global,'creation_date',lcl_date_time)
  rcd=rcd+nf90_put_att(nc_id,nf90_global,'prg_ID',prg_ID(1:ftn_strlen(prg_ID)))
  rcd=rcd+nf90_put_att(nc_id,nf90_global,'cmd_ln',cmd_ln(1:ftn_strlen(cmd_ln)))
  rcd=rcd+nf90_put_att(nc_id,nf90_global,'src_rfr_sng',src_rfr_sng(1:ftn_strlen(src_rfr_sng)))
  rcd=rcd+nf90_put_att(nc_id,nf90_global,'src_fl_sng',src_fl_sng(1:ftn_strlen(src_fl_sng)))
  
  ! Add english text descriptions
  rcd=rcd+nf90_put_att(nc_id,Rayleigh_sca_xsx_id,'long_name','Rayleigh scattering cross section')
  rcd=rcd+nf90_put_att(nc_id,abs_cff_mss_O2_id,'long_name','Oxygen mass absorption coefficient')
  rcd=rcd+nf90_put_att(nc_id,abs_xsx_O2_id,'long_name','Molecular Oxygen Herzberg continuum absorption cross section')
  rcd=rcd+nf90_put_att(nc_id,bnd_id,'long_name','Band center wavelength')
  rcd=rcd+nf90_put_att(nc_id,flx_bnd_dwn_TOA_id,'long_name','Solar Energy flux in band')
  rcd=rcd+nf90_put_att(nc_id,flx_bnd_pht_dwn_TOA_id,'long_name','Photon flux in band')
  rcd=rcd+nf90_put_att(nc_id,flx_slr_frc_id,'long_name','Fraction of solar flux in band: ' // fl_slr)
  rcd=rcd+nf90_put_att(nc_id,flx_spc_dwn_TOA_id,'long_name','Spectral solar insolation at TOA')
  rcd=rcd+nf90_put_att(nc_id,idx_rfr_air_STP_id,'long_name','Index of refraction at band center at STP')
  rcd=rcd+nf90_put_att(nc_id,nrg_pht_id,'long_name','Energy of photon at band center')
  rcd=rcd+nf90_put_att(nc_id,wvl_ctr_id,'long_name','Band center wavelength')
  rcd=rcd+nf90_put_att(nc_id,wvl_grd_id,'long_name','Wavelength grid')
  rcd=rcd+nf90_put_att(nc_id,wvl_max_id,'long_name','Band maximum wavelength')
  rcd=rcd+nf90_put_att(nc_id,wvl_min_id,'long_name','Band minimum wavelength')
  rcd=rcd+nf90_put_att(nc_id,wvl_dlt_id,'long_name','Bandwidth')
  
  ! Add units
  rcd=rcd+nf90_put_att(nc_id,tpt_cold_id,'units','kelvin')
  rcd=rcd+nf90_put_att(nc_id,tpt_std_id,'units','kelvin')
  rcd=rcd+nf90_put_att(nc_id,Rayleigh_sca_xsx_id,'units','meter2')
  rcd=rcd+nf90_put_att(nc_id,abs_cff_mss_O2_id,'units','meter2 kilogram-1')
  rcd=rcd+nf90_put_att(nc_id,abs_xsx_O2_id,'units','meter2')
  rcd=rcd+nf90_put_att(nc_id,bnd_id,'units','meter')
  rcd=rcd+nf90_put_att(nc_id,flx_bnd_dwn_TOA_id,'units','watt meter-2')
  rcd=rcd+nf90_put_att(nc_id,flx_bnd_pht_dwn_TOA_id,'units','photon meter-2 second-1')
  rcd=rcd+nf90_put_att(nc_id,flx_slr_frc_id,'units','fraction')
  rcd=rcd+nf90_put_att(nc_id,flx_spc_dwn_TOA_id,'units','watt meter-2 meter-1')
  rcd=rcd+nf90_put_att(nc_id,idx_rfr_air_STP_id,'units','fraction')
  rcd=rcd+nf90_put_att(nc_id,nrg_pht_id,'units','joule photon-1')
  rcd=rcd+nf90_put_att(nc_id,wvl_ctr_id,'units','meter')
  rcd=rcd+nf90_put_att(nc_id,wvl_grd_id,'units','meter')
  rcd=rcd+nf90_put_att(nc_id,wvl_max_id,'units','meter')
  rcd=rcd+nf90_put_att(nc_id,wvl_min_id,'units','meter')
  rcd=rcd+nf90_put_att(nc_id,wvl_dlt_id,'units','meter')
  
  ! All dimensions, variables, and attributes have been defined
  rcd=rcd+nf90_enddef(nc_id)
  
  ! Write data
  rcd=rcd+nf90_put_var(nc_id,tpt_cold_id,tpt_cold)
  rcd=rcd+nf90_put_var(nc_id,tpt_std_id,tpt_std)
  rcd=rcd+nf90_put_var(nc_id,Rayleigh_sca_xsx_id,Rayleigh_sca_xsx)
  rcd=rcd+nf90_put_var(nc_id,abs_cff_mss_O2_id,abs_cff_mss_O2)
  rcd=rcd+nf90_put_var(nc_id,abs_xsx_O2_id,abs_xsx_O2)
  rcd=rcd+nf90_put_var(nc_id,bnd_id,bnd)
  rcd=rcd+nf90_put_var(nc_id,flx_bnd_dwn_TOA_id,flx_bnd_dwn_TOA)
  rcd=rcd+nf90_put_var(nc_id,flx_bnd_pht_dwn_TOA_id,flx_bnd_pht_dwn_TOA)
  rcd=rcd+nf90_put_var(nc_id,flx_slr_frc_id,flx_slr_frc)
  rcd=rcd+nf90_put_var(nc_id,flx_spc_dwn_TOA_id,flx_spc_dwn_TOA)
  rcd=rcd+nf90_put_var(nc_id,idx_rfr_air_STP_id,idx_rfr_air_STP)
  rcd=rcd+nf90_put_var(nc_id,nrg_pht_id,nrg_pht)
  rcd=rcd+nf90_put_var(nc_id,wvl_ctr_id,wvl_ctr)
  rcd=rcd+nf90_put_var(nc_id,wvl_grd_id,wvl_grd)
  rcd=rcd+nf90_put_var(nc_id,wvl_max_id,wvl_max)
  rcd=rcd+nf90_put_var(nc_id,wvl_min_id,wvl_min)
  rcd=rcd+nf90_put_var(nc_id,wvl_dlt_id,wvl_dlt)
  
  rcd=rcd+nf90_close(nc_id)
  if (rcd /= nf90_noerr) call nf90_err_exit(rcd,fl_out)
  write (6,'(a28,1x,a)') 'Wrote results to netCDF file',fl_out(1:ftn_strlen(fl_out))
  
  ! De-allocate dynamic variables
  if (allocated(Rayleigh_sca_xsx)) deallocate(Rayleigh_sca_xsx,stat=rcd)
  if (allocated(abs_cff_mss_O2)) deallocate(abs_cff_mss_O2,stat=rcd)
  if (allocated(abs_xsx_O2)) deallocate(abs_xsx_O2,stat=rcd)
  if (allocated(abs_xsx_O3)) deallocate(abs_xsx_O3,stat=rcd)
  if (allocated(abs_xsx_O3_cold)) deallocate(abs_xsx_O3_cold,stat=rcd)
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
  
1000 continue
  
  call exit(exit_status)
end program O2

