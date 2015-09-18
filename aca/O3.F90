! $Id$

program O3
  
  ! Purpose: Converts O3 and O2 absorption cross section data to netCDF format.
  
  ! Compilation:
  ! cd ${HOME}/aca; make -W O3.F OPTS=D O3; cd -
  ! cd ${HOME}/aca; make -W O3.F O3; cd -
  ! cd ${HOME}/aca; make OPTS=D O3; cd -
  
  ! Usage:
  ! ncks -H -C -F -v wvl_ctr_CCM,flx_slr_frc_CCM,abs_xsx_O3_CCM ${DATA}/aca/abs_xsx_O3.nc
  ! ncks -H -C -F -d bnd,0.3e-6 -d bnd_CCM,0.3e-6 -v wvl_ctr_CCM,flx_slr_frc_CCM,abs_xsx_O3_CCM,abs_xsx_O3 ${DATA}/aca/abs_xsx_O3.nc
  ! ncwa -a bnd -d bnd,0.295e-6,0.305e-6 ${DATA}/aca/abs_xsx_O3.nc ${DATA}/aca/O3.nc
  ! ncks -H -C -F -d bnd,0.3e-6 -v abs_xsx_O3 ${DATA}/aca/abs_xsx_O3.nc
  ! ncks -H -C -F -v abs_xsx_O3 ${DATA}/aca/O3.nc
  ! ncks -H -C -F -d bnd_CCM,6 -v abs_xsx_O3_CCM ${DATA}/aca/abs_xsx_O3.nc
  
  ! Use WMO85 data
  ! O3 -i ${DATA}/aca/abs_xsx_WMO85.txt -o ${DATA}/aca/abs_xsx_O3.nc
  ! Weight WMO85 data by LaN68
  ! O3 -S ${DATA}/aca/spc_LaN68.nc
  
  ! Currently the code processes input ASCII data files that look like:
  
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
  character(len=*),parameter::sbr_nm='O3' ! [sng] Subroutine name
  
  ! Arrays are allocatable, but die if size exceeds corresponding *_nbr_max
  integer fl_in_unit
  integer bnd_nbr_WMO85
  integer bnd_nbr_max
  parameter(fl_in_unit=73, &
       bnd_nbr_WMO85=158, &
       bnd_nbr_max=bnd_nbr_WMO85) ! max(bnd_nbr_WMO85)
  
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
  
  logical WMO85
  logical std_tpt
  
  integer bnd_dim_id        ! dimension ID for bands
  integer grd_dim_id        ! dimension ID for grid
  integer bnd_idx           ! counting index
  integer nc_id             ! file handle
  integer bnd_nbr ! dimension size
  
  integer Rayleigh_sca_xsx_id
  integer abs_cff_mss_O3_id
  integer abs_xsx_O2_id
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

  ! Main code
  
  ! Initialize default values
  WMO85=.true.
  dbg_lvl=dbg_off
  exit_status=0
  fl_in='/data/zender/aca/abs_xsx_WMO85.txt'
  fl_out='/data/zender/aca/abs_xsx_O3.nc'
  fl_slr='/data/zender/aca/spc_Kur95_01wvn.nc'
  rcd=nf90_noerr              ! nf90_noerr == 0
  CVS_Date='$Date$'
  CVS_Revision='$Revision$'
  std_tpt=.true.
  tpt_cold=203.0
  tpt_std=250.0 ! Temperature at which generic O3 cross sections are archived
  tpt_warm=273.0
  
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
        std_tpt=.not.std_tpt
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
  
  if (WMO85) then           ! WMO85 data
     bnd_nbr=bnd_nbr_WMO85
  else
     stop 'WMO85 is .false. but is only O3 data source'   
  endif                     ! endif WMO85
  
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
  
  ! Compute any quantities that might depend on command line input
  call ftn_strnul(fl_in)
  call ftn_strnul(fl_out)
  call ftn_strnul(fl_slr)
  call ftn_strcpy(src_fl_sng,'Original data file is ' // fl_in)
  if (WMO85) call ftn_strcpy(src_rfr_sng,'Data reference is WMO (1985) (WMO85)')
  
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
     
     ! Define temperature dependance
     ! Temperature dependence is strongest around 
     do bnd_idx=1,bnd_nbr
        abs_xsx_O3_tpt_rfr(bnd_idx)=tpt_std ! Temperature at which generic abs_xsx_O3 array will be valid
        ! For WMO85 data, cold = 203 K, warm = 273 K
        abs_xsx_O3_dadT(bnd_idx)= &
             (abs_xsx_O3_warm(bnd_idx)-abs_xsx_O3_cold(bnd_idx))/ &
             (tpt_warm-tpt_cold) 
        abs_xsx_O3(bnd_idx)=abs_xsx_O3_cold(bnd_idx)+ &
             (tpt_std-tpt_cold)*abs_xsx_O3_dadT(bnd_idx)
     enddo                  ! end loop over bnd
     
  endif                     ! WMO85 data
  
  close (fl_in_unit)
  write (6,'(a20,1x,a)') 'Read input data from',fl_in(1:ftn_strlen(fl_in))
  
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
  rcd=nf90_wrp(nf90_def_var(nc_id,'abs_cff_mss_O3',nf90_float,bnd_dim_id,abs_cff_mss_O3_id),sbr_nm//": dv abs_cff_mss_O3")
  rcd=nf90_wrp(nf90_def_var(nc_id,'abs_xsx_O2',nf90_float,bnd_dim_id,abs_xsx_O2_id),sbr_nm//": dv abs_xsx_O2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'abs_xsx_O3',nf90_float,bnd_dim_id,abs_xsx_O3_id),sbr_nm//": dv abs_xsx_O3")
  rcd=nf90_wrp(nf90_def_var(nc_id,'abs_xsx_O3_cold',nf90_float,bnd_dim_id,abs_xsx_O3_cold_id),sbr_nm//": dv abs_xsx_O3_cold")
  rcd=nf90_wrp(nf90_def_var(nc_id,'abs_xsx_O3_dadT',nf90_float,bnd_dim_id,abs_xsx_O3_dadT_id),sbr_nm//": dv abs_xsx_O3_dadT")
  rcd=nf90_wrp(nf90_def_var(nc_id,'abs_xsx_O3_warm',nf90_float,bnd_dim_id,abs_xsx_O3_warm_id),sbr_nm//": dv abs_xsx_O3_warm")
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
  rcd=nf90_wrp(nf90_def_var(nc_id,'abs_xsx_O3_tpt_rfr',nf90_float,bnd_dim_id,abs_xsx_O3_tpt_rfr_id), &
       sbr_nm//": dv abs_xsx_O3_tpt_rfr")
  
  ! Add global attributes
  rcd=rcd+nf90_put_att(nc_id,nf90_global,'CVS_Id',CVS_Id)
  rcd=rcd+nf90_put_att(nc_id,nf90_global,'creation_date',lcl_date_time)
  rcd=rcd+nf90_put_att(nc_id,nf90_global,'prg_ID',prg_ID(1:ftn_strlen(prg_ID)))
  rcd=rcd+nf90_put_att(nc_id,nf90_global,'cmd_ln',cmd_ln(1:ftn_strlen(cmd_ln)))
  rcd=rcd+nf90_put_att(nc_id,nf90_global,'src_rfr_sng',src_rfr_sng(1:ftn_strlen(src_rfr_sng)))
  rcd=rcd+nf90_put_att(nc_id,nf90_global,'src_fl_sng',src_fl_sng(1:ftn_strlen(src_fl_sng)))
  
  ! Add english text descriptions
  rcd=rcd+nf90_put_att(nc_id,Rayleigh_sca_xsx_id,'long_name','Rayleigh scattering cross section')
  rcd=rcd+nf90_put_att(nc_id,abs_cff_mss_O3_id,'long_name','Ozone mass absorption coefficient')
  rcd=rcd+nf90_put_att(nc_id,abs_xsx_O2_id,'long_name','Molecular Oxygen Herzberg continuum absorption cross section')
  rcd=rcd+nf90_put_att(nc_id,abs_xsx_O3_cold_id,'long_name','Ozone absorption cross section at 203 K')
  rcd=rcd+nf90_put_att(nc_id,abs_xsx_O3_dadT_id,'long_name','Slope of absorption cross section temperature dependence')
  rcd=rcd+nf90_put_att(nc_id,abs_xsx_O3_id,'long_name','Ozone absorption cross section at tpt_rfr')
  rcd=rcd+nf90_put_att(nc_id,abs_xsx_O3_tpt_rfr_id,'long_name','Valid temperature for absorption cross section')
  rcd=rcd+nf90_put_att(nc_id,abs_xsx_O3_warm_id,'long_name','Ozone absorption cross section at 273 K')
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
  rcd=rcd+nf90_put_att(nc_id,Rayleigh_sca_xsx_id,'units','meter2')
  rcd=rcd+nf90_put_att(nc_id,abs_cff_mss_O3_id,'units','meter2 kilogram-1')
  rcd=rcd+nf90_put_att(nc_id,abs_xsx_O2_id,'units','meter2')
  rcd=rcd+nf90_put_att(nc_id,abs_xsx_O3_cold_id,'units','meter2')
  rcd=rcd+nf90_put_att(nc_id,abs_xsx_O3_dadT_id,'units','meter2 kelvin-1')
  rcd=rcd+nf90_put_att(nc_id,abs_xsx_O3_id,'units','meter2')
  rcd=rcd+nf90_put_att(nc_id,abs_xsx_O3_tpt_rfr_id,'units','kelvin')
  rcd=rcd+nf90_put_att(nc_id,abs_xsx_O3_warm_id,'units','meter2')
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
  rcd=rcd+nf90_put_var(nc_id,Rayleigh_sca_xsx_id,Rayleigh_sca_xsx)
  rcd=rcd+nf90_put_var(nc_id,abs_cff_mss_O3_id,abs_cff_mss_O3)
  rcd=rcd+nf90_put_var(nc_id,abs_xsx_O2_id,abs_xsx_O2)
  rcd=rcd+nf90_put_var(nc_id,abs_xsx_O3_cold_id,abs_xsx_O3_cold)
  rcd=rcd+nf90_put_var(nc_id,abs_xsx_O3_dadT_id,abs_xsx_O3_dadT)
  rcd=rcd+nf90_put_var(nc_id,abs_xsx_O3_id,abs_xsx_O3)
  rcd=rcd+nf90_put_var(nc_id,abs_xsx_O3_tpt_rfr_id,abs_xsx_O3_tpt_rfr)
  rcd=rcd+nf90_put_var(nc_id,abs_xsx_O3_warm_id,abs_xsx_O3_warm)
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
  write (6,'(a28,1x,a)') 'Wrote results to netCDF file',fl_out(1:ftn_strlen(fl_out))
  
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
  call nc_out_CCM_SW(fl_out,bnd_dim_id)
  
  ! Add O3 data on CCM grid to netCDF output file
  rcd=rcd+nf90_wrp_open(fl_out,nf90_write,nc_id)
  if (rcd /= nf90_noerr) call nf90_err_exit(rcd,fl_out)
  
  ! Put output file in define mode
  rcd=rcd+nf90_redef(nc_id)
  if (rcd /= nf90_noerr) call nf90_err_exit(rcd,fl_out)
  
  ! Variable definitions.
  rcd=nf90_wrp(nf90_def_var(nc_id,'abs_xsx_O3_CCM',nf90_float,bnd_dim_id,abs_xsx_O3_id),sbr_nm//": dv abs_xsx_O3")
  rcd=nf90_wrp(nf90_def_var(nc_id,'abs_cff_mss_O3_CCM',nf90_float,bnd_dim_id,abs_cff_mss_O3_id),sbr_nm//": dv abs_cff_mss_O3")
  rcd=nf90_wrp(nf90_def_var(nc_id,'flx_slr_frc_CCM',nf90_float,bnd_dim_id,flx_slr_frc_id),sbr_nm//": dv flx_slr_frc")
  rcd=nf90_wrp(nf90_def_var(nc_id,'flx_spc_dwn_TOA_CCM',nf90_float,bnd_dim_id,flx_spc_dwn_TOA_id),sbr_nm//": dv flx_spc_dwn_TOA")
  
  ! Add global attributes
  
  ! Add english text descriptions
  rcd=rcd+nf90_put_att(nc_id,abs_cff_mss_O3_id,'long_name','Ozone mass absorption coefficient')
  rcd=rcd+nf90_put_att(nc_id,abs_xsx_O3_id,'long_name','Ozone continuum absorption cross section')
  rcd=rcd+nf90_put_att(nc_id,flx_slr_frc_id,'long_name','Fraction of solar flux in band: ' // fl_slr)
  rcd=rcd+nf90_put_att(nc_id,flx_spc_dwn_TOA_id,'long_name','Spectral solar insolation at TOA')
  
  ! Add units
  rcd=rcd+nf90_put_att(nc_id,abs_cff_mss_O3_id,'units','meter2 kilogram-1')
  rcd=rcd+nf90_put_att(nc_id,abs_xsx_O3_id,'units','meter2')
  rcd=rcd+nf90_put_att(nc_id,flx_slr_frc_id,'units','fraction')
  rcd=rcd+nf90_put_att(nc_id,flx_spc_dwn_TOA_id,'units','watt meter-2 meter-1')
  
  ! All dimensions, variables, and attributes have been defined
  rcd=rcd+nf90_enddef(nc_id)
  if (rcd /= nf90_noerr) call nf90_err_exit(rcd,fl_out)
  
  ! Write out data
  rcd=rcd+nf90_put_var(nc_id,abs_cff_mss_O3_id,abs_cff_mss_O3_CCM)
  rcd=rcd+nf90_put_var(nc_id,abs_xsx_O3_id,abs_xsx_O3_CCM)
  rcd=rcd+nf90_put_var(nc_id,flx_slr_frc_id,flx_slr_frc_CCM)
  rcd=rcd+nf90_put_var(nc_id,flx_spc_dwn_TOA_id,flx_spc_dwn_TOA_CCM)
  
  ! Close output file
  rcd=rcd+nf90_close(nc_id)
  if (rcd /= nf90_noerr) call nf90_err_exit(rcd,fl_out)
  write (6,'(a28,1x,a)') 'Wrote O3 data on CCM grid to',fl_out(1:ftn_strlen(fl_out))
  
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

