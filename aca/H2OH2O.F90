! $Id$

program H2OH2O
  
  ! Purpose: Convert water collision complex (H2O)2 absorption cross section data to netCDF format
  
  ! Compilation:
  ! cd ${HOME}/aca; make -W H2OH2O.F OPTS=D H2OH2O; cd -
  ! cd ${HOME}/aca; make -W H2OH2O.F H2OH2O; cd -
  ! cd ${HOME}/aca; make OPTS=D H2OH2O; cd -
  
  ! Usage:
  ! ncks -H -C -F -v wvl_ctr_CCM,flx_slr_frc_CCM,abs_xsx_H2OH2O_CCM ${DATA}/aca/abs_xsx_H2OH2O.nc
  ! ncks -H -C -F -d bnd,1.2485e-6 -d bnd_CCM,1.2485e-6 -v wvl_ctr_CCM,flx_slr_frc_CCM,abs_xsx_H2OH2O_CCM,abs_xsx_H2OH2O ${DATA}/aca/abs_xsx_H2OH2O.nc
  ! ncwa -a bnd -d bnd,1.0890e-6,1.4080e-6 ${DATA}/aca/abs_xsx_H2OH2O.nc ${DATA}/aca/foo_H2OH2O.nc
  ! ncks -H -C -F -d bnd,1.2485e-6 -v abs_xsx_H2OH2O ${DATA}/aca/abs_xsx_H2OH2O.nc
  ! ncks -H -C -F -v abs_xsx_H2OH2O ${DATA}/aca/foo_H2OH2O.nc
  ! ncks -H -C -F -d bnd_CCM,10 -v abs_xsx_H2OH2O_CCM ${DATA}/aca/abs_xsx_H2OH2O.nc
  
  ! Currently the code processes input ASCII data files that look like:
  
  ! wave-number   absorption
  !  (cm^{-1})    cross-secton
  !           (10^{-20} cm^2)
  !
  ! 3.007518768 0.0000000000E+00
  ! 6.015037537 0.0000000000E+00
  ! 9.022556305 0.5747804791E-01
  ! 12.03007507 0.7174251229E-01
  
  ! Default to TGC98
  ! H2OH2O 
  ! H2OH2O -T -i ${DATA}/aca/abs_xsx_H2OH2O_TGC98.txt -o ${DATA}/aca/abs_xsx_H2OH2O.nc
  ! Weight CCM cross sections by variable flx_spc_dwn_sfc:
  ! H2OH2O -W 
  ! H2OH2O -W -w ${DATA}/tmp/swnb_trp.nc
  ! H2OH2O -W -w ${DATA}/tmp/swnb_trp.nc -f flx_spc_dwn_sfc
  ! H2OH2O -W -w ${DATA}/tmp/swnb_trp.nc -f trn_spc_atm_ttl
  ! Scale NIR bands by specified factor and weight by flx_spc_dwn_sfc:
  ! H2OH2O -W -N -n 0.6 
  
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
  character*(*) CVS_Id
  
  ! Arrays are allocatable, but die if size exceeds corresponding *_nbr_max
  integer bnd_nbr_GVD98
  integer bnd_nbr_TGC98
  integer bnd_nbr_Chy98
  integer fl_in_unit
  integer bnd_nbr_max
  parameter(fl_in_unit=73, &
       bnd_nbr_GVD98=6, &
       bnd_nbr_TGC98=6650,  & ! Newer, validated bands from Tso
       bnd_nbr_Chy98=8192,  & ! Original bands from Chylek
       bnd_nbr_max=8192,    & ! max(bnd_nbr_GVD98,bnd_nbr_TGC98,bnd_nbr_Chy98)
       CVS_Id='$Id$')
  
  ! Input Arguments
  ! Input/Output Arguments
  ! Output Arguments
  ! Local workspace
  character argv*80
  character cmd_ln*200
  character fl_in*80
  character fl_out*80
  character fl_slr*80
  character fl_wgt*80
  character lbl*80
  character*26::lcl_date_time ! Time formatted as Day Mth DD HH:MM:SS TZ YYYY
  character prg_ID*200
  character CVS_Date*28
  character CVS_Revision*16
  character src_fl_sng*200
  character src_rfr_sng*200
  character src_scl_sng*300
  character src_wgt_sng*300
  character wgt_nm*80
  
  integer arg
  integer arg_nbr
  integer exit_status       ! program exit status
  integer idx
  integer rcd               ! return success code
  
  logical GVD98
  logical TGC98
  logical Chy98
  logical std_tpt
  logical SCL_NIR
  logical WGT_TRN
  
  integer bnd_dim_id        ! dimension ID for bands
  integer grd_dim_id        ! dimension ID for grid
  integer bnd_idx           ! counting index
  integer bnd_nbr           ! dimension size
  integer nc_id             ! file handle
  
  integer abs_cff_mss_H2OH2O_id
  integer abs_xsx_H2OH2O_dadT_id
  integer abs_xsx_H2OH2O_id
  integer abs_xsx_H2OH2O_tpt_rfr_id
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
  integer wvn_ctr_id
  integer wvn_grd_id
  integer wvn_max_id
  integer wvn_min_id
  integer wvn_dlt_id
  
  real scl_NIR_fct
  real tpt_std
  
  ! Allocatable variables
  real,dimension(:),allocatable::abs_cff_mss_H2OH2O
  real,dimension(:),allocatable::abs_xsx_H2OH2O
  real,dimension(:),allocatable::abs_xsx_H2OH2O_dadT
  real,dimension(:),allocatable::abs_xsx_H2OH2O_tpt_rfr
  real,dimension(:),allocatable::bnd     ! coordinate variable
  real,dimension(:),allocatable::flx_bnd_dwn_TOA
  real,dimension(:),allocatable::flx_bnd_pht_dwn_TOA
  real,dimension(:),allocatable::flx_slr_frc
  real,dimension(:),allocatable::flx_spc_dwn_TOA
  real,dimension(:),allocatable::idx_rfr_air_STP
  real,dimension(:),allocatable::nrg_pht
  real,dimension(:),allocatable::wgt_spc
  real,dimension(:),allocatable::wvl     ! coordinate variable
  real,dimension(:),allocatable::wvl_ctr
  real,dimension(:),allocatable::wvl_dlt
  real,dimension(:),allocatable::wvl_grd
  real,dimension(:),allocatable::wvl_max
  real,dimension(:),allocatable::wvl_min
  real,dimension(:),allocatable::wvn_ctr
  real,dimension(:),allocatable::wvn_dlt
  real,dimension(:),allocatable::wvn_grd
  real,dimension(:),allocatable::wvn_max
  real,dimension(:),allocatable::wvn_min
  real,dimension(:),allocatable::xsx_wgt_flx
  
  integer bnd_nbr_CCM
  integer nbr_dat_per_ln
  integer slr_spc_xtr_typ
  integer xtr_typ_LHS
  integer xtr_typ_RHS
  
  real abs_cff_mss_H2OH2O_CCM(bnd_nbr_CCM_SW_max)
  real abs_cff_mss_H2OH2O_CCM_CGS(bnd_nbr_CCM_SW_max)
  real abs_xsx_H2OH2O_CCM(bnd_nbr_CCM_SW_max)
  real bnd_CCM(bnd_nbr_CCM_SW_max)
  real flx_slr_frc_CCM(bnd_nbr_CCM_SW_max)
  real flx_spc_dwn_TOA_CCM(bnd_nbr_CCM_SW_max)
  real foo_CCM(bnd_nbr_CCM_SW_max)
  real wgt_spc_CCM(bnd_nbr_CCM_SW_max)
  real wvl_ctr_CCM(bnd_nbr_CCM_SW_max)
  real wvl_max_CCM(bnd_nbr_CCM_SW_max)
  real wvl_min_CCM(bnd_nbr_CCM_SW_max)
  real wvl_dlt_CCM(bnd_nbr_CCM_SW_max)
  real xsx_wgt_flx_CCM(bnd_nbr_CCM_SW_max)
  ! Main code
  
  ! Initialize default values
  Chy98=.false.
  GVD98=.false.
  TGC98=.true.
  SCL_NIR=.false.
  WGT_TRN=.false.
  dbg_lvl=dbg_off
  exit_status=0
  fl_in='/data/zender/aca/abs_xsx_H2OH2O_TGC98.txt'
  fl_out='/data/zender/aca/abs_xsx_H2OH2O.nc'
  fl_slr='/data/zender/aca/spc_Kur95_01wvn.nc'
  fl_wgt='/data/zender/tmp/swnb_trp.nc'
  rcd=nf90_noerr              ! nf90_noerr == 0
  CVS_Date='$Date$'
  CVS_Revision='$Revision$'
  scl_NIR_fct=1.0
  std_tpt=.true.
  tpt_std=tpt_frz_pnt
  wgt_nm='flx_spc_dwn_sfc'
  
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
     if (argv(1:2) == '-C') then
        Chy98=.true.
        GVD98=.not.Chy98
        TGC98=.not.Chy98
     endif
     if (argv(1:2) == '-f') then
        call getarg(arg+1,argv)
        read (argv,'(a)') wgt_nm
     endif
     if (argv(1:2) == '-G') then
        GVD98=.true.
        TGC98=.not.GVD98
        Chy98=.not.GVD98
     endif
     if (argv(1:2) == '-i') then
        call getarg(arg+1,argv)
        read (argv,'(a)') fl_in
     endif
     if (argv(1:2) == '-n') then
        call getarg(arg+1,argv)
        read (argv,'(f8.3)') scl_NIR_fct
     endif
     if (argv(1:2) == '-N') then
        SCL_NIR=.not.SCL_NIR
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
        TGC98=.true.
        GVD98=.not.TGC98
        Chy98=.not.TGC98
     endif
     !     if (argv(1:2) == '-T') then
     !        std_tpt=.not.std_tpt
     !     endif
     if (argv(1:2) == '-t') then
        call getarg(arg+1,argv)
        read (argv,'(f8.3)') tpt_std
     endif
     if (argv(1:2) == '-v') then
        write (6,'(a)') CVS_Id
        goto 1000
     endif
     if (argv(1:2) == '-w') then
        call getarg(arg+1,argv)
        read (argv,'(a)') fl_wgt
     endif
     if (argv(1:2) == '-W') then
        WGT_TRN=.true.
     endif
  end do
  
  if (GVD98) then
     bnd_nbr=bnd_nbr_GVD98    ! [nbr] Number of bands
  else if (TGC98) then
     bnd_nbr=bnd_nbr_TGC98   ! [nbr] Number of bands
  else if (Chy98) then
     bnd_nbr=bnd_nbr_Chy98   ! [nbr] Number of bands
  else
     error stop 'One of GVD98, TGC98, Chy98 must be .true.'   
  endif                     ! endif
  
  ! Allocate space for dynamic arrays
  allocate(abs_cff_mss_H2OH2O(bnd_nbr),stat=rcd)
  allocate(abs_xsx_H2OH2O(bnd_nbr),stat=rcd)
  allocate(abs_xsx_H2OH2O_dadT(bnd_nbr),stat=rcd)
  allocate(abs_xsx_H2OH2O_tpt_rfr(bnd_nbr),stat=rcd)
  allocate(bnd(bnd_nbr),stat=rcd)     ! coordinate variable
  allocate(flx_bnd_dwn_TOA(bnd_nbr),stat=rcd)
  allocate(flx_bnd_pht_dwn_TOA(bnd_nbr),stat=rcd)
  allocate(flx_slr_frc(bnd_nbr),stat=rcd)
  allocate(flx_spc_dwn_TOA(bnd_nbr),stat=rcd)
  allocate(idx_rfr_air_STP(bnd_nbr),stat=rcd)
  allocate(nrg_pht(bnd_nbr),stat=rcd)
  allocate(wgt_spc(bnd_nbr),stat=rcd)
  allocate(wvl(bnd_nbr),stat=rcd)     ! coordinate variable
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
  allocate(xsx_wgt_flx(bnd_nbr),stat=rcd)
  
  ! Compute any quantities that might depend on command line input
  call ftn_strnul(fl_in)
  call ftn_strnul(fl_out)
  call ftn_strnul(fl_slr)
  call ftn_strnul(fl_wgt)
  call ftn_strnul(wgt_nm)
  call ftn_strcpy(src_fl_sng,'Original data file is '//fl_in)
  if (WGT_TRN) then
     call ftn_strcpy(src_wgt_sng,'CCM cross sections are averages of high resolution cross-sections ' &
          //'weighted by variable '//wgt_nm(1:ftn_strlen(wgt_nm))// &
          ' from model atmosphere in '//fl_wgt(1:ftn_strlen(fl_wgt)))
  else
     call ftn_strcpy(src_wgt_sng,'CCM cross sections are averages of high resolution cross-sections ' &
          //'weighted by TOA solar spectral flux')
  endif                     ! endif WGT_TRN
  if (SCL_NIR) then
     write (src_scl_sng,'(a,f9.6)') 'CCM cross sections are tuned to SWNB by scaling NIR bands by ',scl_NIR_fct
     call ftn_strnul(src_scl_sng)
  else
     call ftn_strcpy(src_scl_sng,'CCM cross sections contain no arbitrary scaling')
  endif                     ! endif SCL_NIR
  if (GVD98) then
     call ftn_strcpy(src_rfr_sng,'Data source is Goss et al. (1998) (GVD98)')
     error stop 'GVD98 not implemented yet'
  endif
  if (TGC98) then
     call ftn_strcpy(src_rfr_sng,'Data source is Tso, Geldart, and Chylek (1998) (TGC98)'// &
          ' [6650 bins at 3 cm-1 resolution from 3 cm-1 to 20000 cm-1]')
  endif
  if (Chy98) then
     call ftn_strcpy(src_rfr_sng,'Data source is Chylek (1998) (Chy98)'// &
          ' [8192 bins at 3 cm-1 resolution from 3 cm-1 to 24637 cm-1]')
  endif
  if (.not.std_tpt) then
     write (6,'(a)') 'WARNING: No data have temperature dependence yet, std_tpt ignored'
  endif
  
  open (fl_in_unit,file=fl_in,status='old',iostat=rcd)
  
  if (TGC98) then       ! TGC98 data
     
     bnd_nbr=bnd_nbr_TGC98
     do idx=1,18
        read (fl_in_unit,'(a80)') lbl
     enddo
     lbl(1:1)=lbl(1:1) ! CEWI
     do bnd_idx=1,bnd_nbr
        read (fl_in_unit,*)  &
             wvn_ctr(bnd_idx), &
             abs_xsx_H2OH2O(bnd_idx)
     enddo
     
     ! Convert input data to SI units where necessary
     do bnd_idx=1,bnd_nbr
        abs_xsx_H2OH2O(bnd_idx)=abs_xsx_H2OH2O(bnd_idx)*1.0e-24 ! 10^{-20} cm2 mlc-1 --> m2 mlc-1
     enddo
     
     ! Compute diagnostic variables
     wvn_min(1)=wvn_ctr(1)-0.5*(wvn_ctr(2)-wvn_ctr(1))
     wvn_max(bnd_nbr)=wvn_ctr(bnd_nbr)+0.5*(wvn_ctr(bnd_nbr)-wvn_ctr(bnd_nbr-1))
     do bnd_idx=2,bnd_nbr
        wvn_min(bnd_idx)=0.5*(wvn_ctr(bnd_idx-1)+wvn_ctr(bnd_idx))
     enddo
     do bnd_idx=1,bnd_nbr-1
        wvn_max(bnd_idx)=0.5*(wvn_ctr(bnd_idx)+wvn_ctr(bnd_idx+1))
     enddo
     
     ! Switch from monotonically increasing wavenumber grid to monotonically increasing wavelength grid
     call rvr_vec(wvn_ctr,bnd_nbr)
     call rvr_vec(wvn_min,bnd_nbr)
     call rvr_vec(wvn_max,bnd_nbr)
     call rvr_vec(abs_xsx_H2OH2O,bnd_nbr)
     
     ! Define wavelength grid
     do bnd_idx=1,bnd_nbr
        wvl_min(bnd_idx)=1.0/(100.0*wvn_max(bnd_idx))
        wvl_max(bnd_idx)=1.0/(100.0*wvn_min(bnd_idx))
        wvl(bnd_idx)=1.0/(100.0*wvn_ctr(bnd_idx))
        wvl_ctr(bnd_idx)=0.5*(wvl_min(bnd_idx)+wvl_max(bnd_idx))
        bnd(bnd_idx)=wvl(bnd_idx)
     enddo
     
     do bnd_idx=1,bnd_nbr
        abs_xsx_H2OH2O_tpt_rfr(bnd_idx)=tpt_frz_pnt ! all TGC98 data taken at 0 C
        abs_xsx_H2OH2O_dadT(bnd_idx)=0. ! no temperature dependence given for TGC98 data
     enddo
     
  endif                     ! Chy98 data
  
  close (fl_in_unit)
  write (6,'(a20,1x,a)') 'Read input data from',fl_in(1:ftn_strlen(fl_in))
  
  ! Get TOA solar spectrum
  slr_spc_xtr_typ=xtr_fll_ngh !+xtr_vrb
  call slr_spc_get(fl_slr,wvl_min,wvl_max,bnd_nbr,flx_slr_frc,slr_spc_xtr_typ,slr_spc_xtr_typ)
  
  ! Compute diagnostic variables
  do bnd_idx=1,bnd_nbr
     abs_cff_mss_H2OH2O(bnd_idx)=abs_xsx_H2OH2O(bnd_idx)*Avagadro/mmw_H2OH2O
     wvl_dlt(bnd_idx)=wvl_max(bnd_idx)-wvl_min(bnd_idx)
     wvl_grd(bnd_idx)=wvl_min(bnd_idx)
     nrg_pht(bnd_idx)=Planck*speed_of_light/wvl(bnd_idx)
     flx_bnd_dwn_TOA(bnd_idx)=flx_slr_frc(bnd_idx)*slr_cst_CCM
     flx_spc_dwn_TOA(bnd_idx)=flx_slr_frc(bnd_idx)*slr_cst_CCM/wvl_dlt(bnd_idx)
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
  rcd=rcd+nf90_def_var(nc_id,'abs_cff_mss_H2OH2O',nf90_float,bnd_dim_id,abs_cff_mss_H2OH2O_id)
  rcd=rcd+nf90_def_var(nc_id,'abs_xsx_H2OH2O',nf90_float,bnd_dim_id,abs_xsx_H2OH2O_id)
  rcd=rcd+nf90_def_var(nc_id,'abs_xsx_H2OH2O_dadT',nf90_float,bnd_dim_id,abs_xsx_H2OH2O_dadT_id)
  rcd=rcd+nf90_def_var(nc_id,'abs_xsx_H2OH2O_tpt_rfr',nf90_float,bnd_dim_id,abs_xsx_H2OH2O_tpt_rfr_id)
  rcd=rcd+nf90_def_var(nc_id,'bnd',nf90_float,bnd_dim_id,bnd_id)
  rcd=rcd+nf90_def_var(nc_id,'flx_bnd_dwn_TOA',nf90_float,bnd_dim_id,flx_bnd_dwn_TOA_id)
  rcd=rcd+nf90_def_var(nc_id,'flx_bnd_pht_dwn_TOA',nf90_float,bnd_dim_id,flx_bnd_pht_dwn_TOA_id)
  rcd=rcd+nf90_def_var(nc_id,'flx_slr_frc',nf90_float,bnd_dim_id,flx_slr_frc_id)
  rcd=rcd+nf90_def_var(nc_id,'flx_spc_dwn_TOA',nf90_float,bnd_dim_id,flx_spc_dwn_TOA_id)
  rcd=rcd+nf90_def_var(nc_id,'idx_rfr_air_STP',nf90_float,bnd_dim_id,idx_rfr_air_STP_id)
  rcd=rcd+nf90_def_var(nc_id,'nrg_pht',nf90_float,bnd_dim_id,nrg_pht_id)
  rcd=rcd+nf90_def_var(nc_id,'wvl_ctr',nf90_float,bnd_dim_id,wvl_ctr_id)
  rcd=rcd+nf90_def_var(nc_id,'wvl_grd',nf90_float,grd_dim_id,wvl_grd_id)
  rcd=rcd+nf90_def_var(nc_id,'wvl_max',nf90_float,bnd_dim_id,wvl_max_id)
  rcd=rcd+nf90_def_var(nc_id,'wvl_min',nf90_float,bnd_dim_id,wvl_min_id)
  rcd=rcd+nf90_def_var(nc_id,'wvl_dlt',nf90_float,bnd_dim_id,wvl_dlt_id)
  rcd=rcd+nf90_def_var(nc_id,'wvn_ctr',nf90_float,bnd_dim_id,wvn_ctr_id)
  rcd=rcd+nf90_def_var(nc_id,'wvn_grd',nf90_float,grd_dim_id,wvn_grd_id)
  rcd=rcd+nf90_def_var(nc_id,'wvn_max',nf90_float,bnd_dim_id,wvn_max_id)
  rcd=rcd+nf90_def_var(nc_id,'wvn_min',nf90_float,bnd_dim_id,wvn_min_id)
  rcd=rcd+nf90_def_var(nc_id,'wvn_dlt',nf90_float,bnd_dim_id,wvn_dlt_id)
  
  ! Add global attributes
  rcd=rcd+nf90_put_att(nc_id,nf90_global,'CVS_Id',CVS_Id)
  rcd=rcd+nf90_put_att(nc_id,nf90_global,'creation_date',lcl_date_time)
  rcd=rcd+nf90_put_att(nc_id,nf90_global,'prg_ID',prg_ID(1:ftn_strlen(prg_ID)))
  rcd=rcd+nf90_put_att(nc_id,nf90_global,'cmd_ln',cmd_ln(1:ftn_strlen(cmd_ln)))
  rcd=rcd+nf90_put_att(nc_id,nf90_global,'src_rfr_sng',src_rfr_sng(1:ftn_strlen(src_rfr_sng)))
  rcd=rcd+nf90_put_att(nc_id,nf90_global,'src_fl_sng',src_fl_sng(1:ftn_strlen(src_fl_sng)))
  rcd=rcd+nf90_put_att(nc_id,nf90_global,'src_wgt_sng',src_wgt_sng(1:ftn_strlen(src_wgt_sng)))
  rcd=rcd+nf90_put_att(nc_id,nf90_global,'src_scl_sng',src_scl_sng(1:ftn_strlen(src_scl_sng)))
  
  ! Add english text descriptions
  rcd=rcd+nf90_put_att(nc_id,abs_cff_mss_H2OH2O_id,'long_name','Water dimer mass absorption coefficient')
  rcd=rcd+nf90_put_att(nc_id,abs_xsx_H2OH2O_dadT_id,'long_name','Temperature dependence of absorption cross section')
  rcd=rcd+nf90_put_att(nc_id,abs_xsx_H2OH2O_id,'long_name','Water dimer continuum absorption cross section')
  rcd=rcd+nf90_put_att(nc_id,abs_xsx_H2OH2O_tpt_rfr_id,'long_name','Valid temperature for absorption cross section')
  rcd=rcd+nf90_put_att(nc_id,bnd_id,'long_name','Band nominal wavelength')
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
  rcd=rcd+nf90_put_att(nc_id,wvn_ctr_id,'long_name','Band center wavenumber')
  rcd=rcd+nf90_put_att(nc_id,wvn_grd_id,'long_name','Wavenumber grid')
  rcd=rcd+nf90_put_att(nc_id,wvn_max_id,'long_name','Band maximum wavenumber')
  rcd=rcd+nf90_put_att(nc_id,wvn_min_id,'long_name','Band minimum wavenumber')
  rcd=rcd+nf90_put_att(nc_id,wvn_dlt_id,'long_name','Bandwidth')
  
  ! Add units
  rcd=rcd+nf90_put_att(nc_id,abs_cff_mss_H2OH2O_id,'units','meter2 kilogram-1')
  rcd=rcd+nf90_put_att(nc_id,abs_xsx_H2OH2O_dadT_id,'units','meter2 kelvin-1')
  rcd=rcd+nf90_put_att(nc_id,abs_xsx_H2OH2O_id,'units','meter2')
  rcd=rcd+nf90_put_att(nc_id,abs_xsx_H2OH2O_tpt_rfr_id,'units','kelvin')
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
  rcd=rcd+nf90_put_att(nc_id,wvn_ctr_id,'units','centimeter-1')
  rcd=rcd+nf90_put_att(nc_id,wvn_grd_id,'units','centimeter-1')
  rcd=rcd+nf90_put_att(nc_id,wvn_max_id,'units','centimeter-1')
  rcd=rcd+nf90_put_att(nc_id,wvn_min_id,'units','centimeter-1')
  rcd=rcd+nf90_put_att(nc_id,wvn_dlt_id,'units','centimeter-1')
  
  ! All dimensions, variables, and attributes have been defined
  rcd=rcd+nf90_enddef(nc_id)
  
  ! Write data
  rcd=rcd+nf90_put_var(nc_id,abs_cff_mss_H2OH2O_id,abs_cff_mss_H2OH2O)
  rcd=rcd+nf90_put_var(nc_id,abs_xsx_H2OH2O_dadT_id,abs_xsx_H2OH2O_dadT)
  rcd=rcd+nf90_put_var(nc_id,abs_xsx_H2OH2O_id,abs_xsx_H2OH2O)
  rcd=rcd+nf90_put_var(nc_id,abs_xsx_H2OH2O_tpt_rfr_id,abs_xsx_H2OH2O_tpt_rfr)
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
  rcd=rcd+nf90_put_var(nc_id,wvn_ctr_id,wvn_ctr)
  rcd=rcd+nf90_put_var(nc_id,wvn_grd_id,wvn_grd)
  rcd=rcd+nf90_put_var(nc_id,wvn_max_id,wvn_max)
  rcd=rcd+nf90_put_var(nc_id,wvn_min_id,wvn_min)
  rcd=rcd+nf90_put_var(nc_id,wvn_dlt_id,wvn_dlt)
  
  rcd=rcd+nf90_close(nc_id)
  write (6,'(a28,1x,a)') 'Wrote results to netCDF file',fl_out(1:ftn_strlen(fl_out))
  
  ! Get CCM wavelength grid
  call wvl_grd_CCM_SW_mk(bnd_nbr_CCM,bnd_CCM,wvl_min_CCM,wvl_max_CCM,wvl_ctr_CCM,wvl_dlt_CCM)
  
  ! Get TOA solar spectrum
  call slr_spc_get_CCM(fl_slr,wvl_min_CCM,wvl_max_CCM,bnd_nbr_CCM,flx_slr_frc_CCM,slr_spc_xtr_typ,slr_spc_xtr_typ)
  
  if (WGT_TRN) then
     ! Weight high resolution absorption cross-sections by atmospheric transmission of atmosphere with all constituents except H2OH2O
     call wgt_get(fl_wgt,wgt_nm,wvl_grd,bnd_nbr,wgt_spc)         
     do idx=1,bnd_nbr
        xsx_wgt_flx(idx)=abs_xsx_H2OH2O(idx)*wgt_spc(idx)
     enddo                  ! end loop over bnd
  else                      ! !WGT_TRN
     ! Weight high resolution absorption cross sections by high resolution TOA solar flux
     do idx=1,bnd_nbr
        xsx_wgt_flx(idx)=abs_xsx_H2OH2O(idx)*flx_spc_dwn_TOA(idx)
     enddo                  ! end loop over bnd
  endif                     ! !WGT_TRN
  
  ! Initialize default values
  xtr_typ_LHS=xtr_prt_wgt+xtr_fll_nil !+xtr_vrb
  xtr_typ_RHS=xtr_prt_wgt+xtr_fll_nil !+xtr_vrb
  
  ! Rebin solar flux
  call rbn_vec_CCM(bnd_nbr,wvl_grd,flx_spc_dwn_TOA, &
       bnd_nbr_CCM,wvl_min_CCM,wvl_max_CCM,foo_CCM, &
       xtr_typ_LHS,xtr_typ_RHS)
  ! Rebin flux-weighted absorption cross sections
  call rbn_vec_CCM(bnd_nbr,wvl_grd,xsx_wgt_flx,  &
       bnd_nbr_CCM,wvl_min_CCM,wvl_max_CCM,xsx_wgt_flx_CCM, &
       xtr_typ_LHS,xtr_typ_RHS)
  if (WGT_TRN) then
     ! Rebin cross-section weights
     call rbn_vec_CCM(bnd_nbr,wvl_grd,wgt_spc, &
          bnd_nbr_CCM,wvl_min_CCM,wvl_max_CCM,wgt_spc_CCM, &
          xtr_typ_LHS,xtr_typ_RHS)
  endif                     ! endif WGT_TRN
  
  if (dbg_lvl == dbg_crr) then
     ! Examine weights
     write (6,'(5(a,1x))') '#  ','wvl        ','xsx_abs    ','flx TOA    ',wgt_nm(1:ftn_strlen(wgt_nm))
     write (6,'(5(a,1x))') '   ','m          ','m-2 mlc-1  ','W m-2 m-1  ','W m-2 m-1  '
     do idx=1,bnd_nbr
        write (6,'(i4,1x,4(es10.3,1x))') &
             idx,wvl_ctr(idx),abs_xsx_H2OH2O(idx),flx_spc_dwn_TOA(idx),wgt_spc(idx)
     enddo                  ! end loop over bnd
  endif                     ! endif dbg
  
  ! Normalize flux-weighted absorption cross sections
  do idx=1,bnd_nbr_CCM
     flx_spc_dwn_TOA_CCM(idx)=flx_slr_frc_CCM(idx)*slr_cst_CCM/wvl_dlt_CCM(idx)
     if (WGT_TRN) then
        if (wgt_spc_CCM(idx) /= 0.) then
           abs_xsx_H2OH2O_CCM(idx)=xsx_wgt_flx_CCM(idx)/wgt_spc_CCM(idx)
        else
           write (6,'(a,a,i2,a)') prg_nm,': WARNING wgt_spc_CCM(',idx,') = 0. Setting cross section equal to zero.'
           abs_xsx_H2OH2O_CCM(idx)=0.
        endif               ! endif not degenerate
     else
        abs_xsx_H2OH2O_CCM(idx)=xsx_wgt_flx_CCM(idx)/flx_spc_dwn_TOA_CCM(idx)
     endif                  ! endif WGT_TRN
     if (SCL_NIR) then
        ! Empirical fudge factor to remove 3 W m-2 MLS and TRP bias
        if (idx >= 9) then
           abs_xsx_H2OH2O_CCM(idx)=abs_xsx_H2OH2O_CCM(idx)*scl_NIR_fct
        endif               ! endif NIR bnd
     endif                  ! endif WGT_TRN
     abs_cff_mss_H2OH2O_CCM(idx)=abs_xsx_H2OH2O_CCM(idx)*Avagadro/mmw_H2OH2O
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
  
  ! Add H2OH2O data on CCM grid to netCDF output file
  rcd=rcd+nf90_wrp_open(fl_out,nf90_write,nc_id)
  if (rcd /= nf90_noerr) call nf90_err_exit(rcd,fl_out)
  
  ! Put output file in define mode
  rcd=rcd+nf90_redef(nc_id)
  if (rcd /= nf90_noerr) call nf90_err_exit(rcd,fl_out)
  
  ! Variable definitions.
  rcd=rcd+nf90_def_var(nc_id,'abs_xsx_H2OH2O_CCM',nf90_float,bnd_dim_id,abs_xsx_H2OH2O_id)
  rcd=rcd+nf90_def_var(nc_id,'abs_cff_mss_H2OH2O_CCM',nf90_float,bnd_dim_id,abs_cff_mss_H2OH2O_id)
  rcd=rcd+nf90_def_var(nc_id,'flx_slr_frc_CCM',nf90_float,bnd_dim_id,flx_slr_frc_id)
  rcd=rcd+nf90_def_var(nc_id,'flx_spc_dwn_TOA_CCM',nf90_float,bnd_dim_id,flx_spc_dwn_TOA_id)
  
  ! Add global attributes
  
  ! Add english text descriptions
  rcd=rcd+nf90_put_att(nc_id,abs_cff_mss_H2OH2O_id,'long_name','Water dimer mass absorption coefficient')
  rcd=rcd+nf90_put_att(nc_id,abs_xsx_H2OH2O_id,'long_name','Water dimer continuum absorption cross section')
  rcd=rcd+nf90_put_att(nc_id,flx_slr_frc_id,'long_name','Fraction of solar flux in band: ' // fl_slr)
  rcd=rcd+nf90_put_att(nc_id,flx_spc_dwn_TOA_id,'long_name','Spectral solar insolation at TOA')
  
  ! Add units
  rcd=rcd+nf90_put_att(nc_id,abs_cff_mss_H2OH2O_id,'units','meter2 kilogram-1')
  rcd=rcd+nf90_put_att(nc_id,abs_xsx_H2OH2O_id,'units','meter2')
  rcd=rcd+nf90_put_att(nc_id,flx_slr_frc_id,'units','fraction')
  rcd=rcd+nf90_put_att(nc_id,flx_spc_dwn_TOA_id,'units','watt meter-2 meter-1')
  
  ! All dimensions, variables, and attributes have been defined
  rcd=rcd+nf90_enddef(nc_id)
  if (rcd /= nf90_noerr) call nf90_err_exit(rcd,fl_out)
  
  ! Write out data
  rcd=rcd+nf90_put_var(nc_id,abs_cff_mss_H2OH2O_id,abs_cff_mss_H2OH2O_CCM)
  rcd=rcd+nf90_put_var(nc_id,abs_xsx_H2OH2O_id,abs_xsx_H2OH2O_CCM)
  rcd=rcd+nf90_put_var(nc_id,flx_slr_frc_id,flx_slr_frc_CCM)
  rcd=rcd+nf90_put_var(nc_id,flx_spc_dwn_TOA_id,flx_spc_dwn_TOA_CCM)
  
  ! Close output file
  rcd=rcd+nf90_close(nc_id)
  write (6,'(a32,1x,a)') 'Wrote H2OH2O data on CCM grid to',fl_out(1:ftn_strlen(fl_out))
  if (rcd /= nf90_noerr) write (6,'(a,a,i4,a)') prg_nm,': ERROR rcd = ',rcd,' on exit'
  
  ! Convert absorption coefficients to [cm2 g-1] and output block for radcsw
  do idx=1,bnd_nbr_CCM
     abs_cff_mss_H2OH2O_CCM_CGS(idx)=abs_cff_mss_H2OH2O_CCM(idx)*10. !  m2 kg-1 --> [cm2 g-1] 
  enddo                     ! end loop over CCM bnd
  nbr_dat_per_ln=4
  write (0,'(a,a,a,a)') 'c     Data generated by ',prg_nm(1:ftn_strlen(prg_nm)),' on ',lcl_date_time
  write (0,'(a,a)') 'c     ',prg_ID(1:ftn_strlen(prg_ID))
  write (0,'(a,a)') 'c     Command line: ',cmd_ln(1:ftn_strlen(cmd_ln))
  write (0,'(a,f7.3,a)') 'c     H2OH2O mass absorption coefficients in [cm2 g-1] for T = ',tpt_std,' K'
  write (0,'(a,a)') 'c     ',src_rfr_sng(1:ftn_strlen(src_rfr_sng))
  write (0,'(a,a)') 'c     ',src_fl_sng(1:ftn_strlen(src_fl_sng))
  write (0,'(a,a)') 'c     ',src_wgt_sng(1:ftn_strlen(src_wgt_sng))
  write (0,'(a,a)') 'c     ',src_scl_sng(1:ftn_strlen(src_scl_sng))
  call dat_prn_f77(bnd_nbr_CCM,abs_cff_mss_H2OH2O_CCM_CGS,nbr_dat_per_ln,'abs_cff_mss_H2OH2O'//char(0))
  call dat_prn_f90(bnd_nbr_CCM,abs_cff_mss_H2OH2O_CCM_CGS,nbr_dat_per_ln,'abs_cff_mss_H2OH2O'//char(0))
  
  ! De-allocate dynamic variables
  if (allocated(abs_cff_mss_H2OH2O)) deallocate(abs_cff_mss_H2OH2O,stat=rcd)
  if (allocated(abs_xsx_H2OH2O)) deallocate(abs_xsx_H2OH2O,stat=rcd)
  if (allocated(abs_xsx_H2OH2O_dadT)) deallocate(abs_xsx_H2OH2O_dadT,stat=rcd)
  if (allocated(abs_xsx_H2OH2O_tpt_rfr)) deallocate(abs_xsx_H2OH2O_tpt_rfr,stat=rcd)
  if (allocated(bnd)) deallocate(bnd,stat=rcd)     ! coordinate variable
  if (allocated(flx_bnd_dwn_TOA)) deallocate(flx_bnd_dwn_TOA,stat=rcd)
  if (allocated(flx_bnd_pht_dwn_TOA)) deallocate(flx_bnd_pht_dwn_TOA,stat=rcd)
  if (allocated(flx_slr_frc)) deallocate(flx_slr_frc,stat=rcd)
  if (allocated(flx_spc_dwn_TOA)) deallocate(flx_spc_dwn_TOA,stat=rcd)
  if (allocated(idx_rfr_air_STP)) deallocate(idx_rfr_air_STP,stat=rcd)
  if (allocated(nrg_pht)) deallocate(nrg_pht,stat=rcd)
  if (allocated(wgt_spc)) deallocate(wgt_spc,stat=rcd)
  if (allocated(wvl)) deallocate(wvl,stat=rcd)     ! coordinate variable
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
  if (allocated(xsx_wgt_flx)) deallocate(xsx_wgt_flx,stat=rcd)
  
1000 continue
  
  call exit(exit_status)
end program H2OH2O
