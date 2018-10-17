! $Id$

program NO2
  
  ! Purpose: Convert NO2 absorption cross section data to netCDF format
  
  ! Compilation:
  ! cd ${HOME}/aca; make -W NO2.F OPTS=D NO2; cd -
  ! cd ${HOME}/aca; make -W NO2.F NO2; cd -
  ! cd ${HOME}/aca; make OPTS=D NO2; cd -
  
  ! Usage:
  ! ncks -H -C -F -v wvl_ctr_CCM,flx_slr_frc_CCM,abs_xsx_CCM ${DATA}/aca/abs_xsx_NO2.nc
  ! ncks -H -C -F -d bnd,0.3e-6 -d bnd_CCM,0.3e-6 -v wvl_ctr_CCM,flx_slr_frc_CCM,abs_xsx_CCM,abs_xsx ${DATA}/aca/abs_xsx_NO2.nc
  ! ncwa -a bnd -d bnd,0.295e-6,0.305e-6 ${DATA}/aca/abs_xsx_NO2.nc ${DATA}/aca/NO2.nc
  ! ncks -H -C -F -d bnd,0.3e-6 -v abs_xsx ${DATA}/aca/abs_xsx_NO2.nc
  ! ncks -H -C -F -v abs_xsx ${DATA}/aca/NO2.nc
  ! ncks -H -C -F -d bnd_CCM,6 -v abs_xsx_CCM ${DATA}/aca/abs_xsx_NO2.nc
  
  ! Default to JPL data, renormalize cross sections to 0 C
  ! NO2 
  ! NO2 -J -t 273.15 -i ${DATA}/tuv/DATAE1/NO2/NO2_jpl94.abs -o ${DATA}/aca/abs_xsx_NO2.nc
  
  ! Use JPL data, renormalize cross sections to 300 K
  ! NO2 -t 300
  
  ! Use NCAR data, turn off temperature dependence
  ! NO2 -N -T -i ${DATA}/tuv/DATAE1/NO2/NO2_ncar_00.abs -o ${DATA}/aca/abs_xsx_NO2.nc
  
  ! Use NOAA data
  ! NO2 -H -t 294 -i ${DATA}/aca/NO2/abs_xsx_NO2_294K.txt -o ${DATA}/aca/NO2/abs_xsx_NO2_294K.nc
  
  ! Currently the code processes input ASCII data files that look like:
  
  ! Table 13.  Absorption Cross Sections of NO2
  ! (a* in units of 1E22 cm2 mol-1 deg-1)
  !-----------------------------------------------------------------------------
  ! lambda1         lambda2         sigma @ T0              T0      a*
  ! (nm)            (nm)            (1E-20 cm2 mol-1)       (C)     (see above)
  !-----------------------------------------------------------------------------
  ! 202.02          204.08          41.45                   25       0.
  ! 204.08          206.19          44.78                   25       0.
  ! 206.19          208.33          44.54                   25       0.
  
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
  integer bnd_nbr_Cal87
  integer bnd_nbr_JPL
  integer bnd_nbr_NCAR
  integer bnd_nbr_NOAA
  integer bnd_nbr_max
  integer fl_in_unit
  integer fl_yld_unit
  parameter(fl_in_unit=73, &
       fl_yld_unit=74, &
       bnd_nbr_Cal87=65, &
       bnd_nbr_JPL=57, &
       bnd_nbr_NCAR=750, &
       bnd_nbr_NOAA=205001, &
       bnd_nbr_max=205001,     & ! max(bnd_nbr_JPL,bnd_nbr_NCAR,bnd_nbr_NOAA)
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
  character fl_yld*80
  character lbl*80
  character*26::lcl_date_time ! Time formatted as Day Mth DD HH:MM:SS TZ YYYY
  character prg_ID*200
  character CVS_Date*28
  character CVS_Revision*16
  character src_fl_sng*200
  character src_rfr_sng*200
  
  integer arg
  integer arg_nbr
  integer exit_status       ! program exit status
  integer idx
  integer rcd               ! return success code
  
  logical Cal87
  logical JPL
  logical NCAR
  logical NOAA
  logical std_tpt
  
  real float_foo
  real tpt_std
  
  integer bnd_dim_id        ! dimension ID for bands
  integer grd_dim_id        ! dimension ID for grid
  integer bnd_idx           ! counting index
  integer bnd_nbr           ! dimension size
  integer nc_id             ! file handle
  
  integer abs_cff_mss_id
  integer abs_xsx_dadT_id
  integer abs_xsx_id
  integer abs_xsx_tpt_rfr_id
  integer bnd_id            ! coordinate ID
  integer flx_bnd_dwn_TOA_id
  integer flx_bnd_pht_dwn_TOA_id
  integer flx_slr_frc_id
  integer flx_spc_dwn_TOA_id
  integer idx_rfr_air_STP_id
  integer nrg_pht_id
  integer qnt_yld_id
  integer wvl_ctr_id
  integer wvl_grd_id
  integer wvl_max_id
  integer wvl_min_id
  integer wvl_dlt_id
  
  ! Allocatable variables
  real,dimension(:),allocatable::abs_cff_mss
  real,dimension(:),allocatable::abs_xsx
  real,dimension(:),allocatable::abs_xsx_dadT
  real,dimension(:),allocatable::abs_xsx_tpt_rfr
  real,dimension(:),allocatable::bnd     ! coordinate variable
  real,dimension(:),allocatable::flx_bnd_dwn_TOA
  real,dimension(:),allocatable::flx_bnd_pht_dwn_TOA
  real,dimension(:),allocatable::flx_slr_frc
  real,dimension(:),allocatable::flx_spc_dwn_TOA
  real,dimension(:),allocatable::idx_rfr_air_STP
  real,dimension(:),allocatable::nrg_pht
  real,dimension(:),allocatable::qnt_yld
  real,dimension(:),allocatable::wvl     ! coordinate variable
  real,dimension(:),allocatable::wvl_ctr
  real,dimension(:),allocatable::wvl_dlt
  real,dimension(:),allocatable::wvl_grd
  real,dimension(:),allocatable::wvl_max
  real,dimension(:),allocatable::wvl_min
  real,dimension(:),allocatable::xsx_wgt_flx
  
  real wvl_ctr_yyy(bnd_nbr_Cal87)
  real wvl_grd_yyy(bnd_nbr_Cal87+1)
  real wvl_max_yyy(bnd_nbr_Cal87)
  real wvl_min_yyy(bnd_nbr_Cal87)
  real qnt_yld_yyy(bnd_nbr_Cal87)
  
  integer bnd_nbr_CCM
  integer nbr_dat_per_ln
  integer slr_spc_xtr_typ
  integer xtr_typ_LHS
  integer xtr_typ_RHS
  
  real abs_cff_mss_CCM(bnd_nbr_CCM_SW_max)
  real abs_cff_mss_CCM_CGS(bnd_nbr_CCM_SW_max)
  real foo_CCM(bnd_nbr_CCM_SW_max)
  real abs_xsx_CCM(bnd_nbr_CCM_SW_max)
  real flx_spc_dwn_TOA_CCM(bnd_nbr_CCM_SW_max)
  real flx_slr_frc_CCM(bnd_nbr_CCM_SW_max)
  real bnd_CCM(bnd_nbr_CCM_SW_max)
  real wvl_ctr_CCM(bnd_nbr_CCM_SW_max)
  real wvl_max_CCM(bnd_nbr_CCM_SW_max)
  real wvl_min_CCM(bnd_nbr_CCM_SW_max)
  real wvl_dlt_CCM(bnd_nbr_CCM_SW_max)
  real xsx_wgt_flx_CCM(bnd_nbr_CCM_SW_max)

  ! Main code
  
  ! Initialize default values
  Cal87=.true.
  JPL=.true.
  NCAR=.false.
  NOAA=.false.
  dbg_lvl=dbg_off
  exit_status=0
  fl_in='/data/zender/tuv/DATAE1/NO2/NO2_jpl94.abs'
  fl_out='/data/zender/aca/abs_xsx_NO2.nc'
  fl_slr='/data/zender/aca/spc_Kur95_01wvn.nc'
  fl_yld='/data/zender/tuv/DATAJ1/YLD/NO2_calvert.yld'
  float_foo=0
  rcd=nf90_noerr              ! nf90_noerr == 0
  CVS_Date='$Date$'
  CVS_Revision='$Revision$'
  std_tpt=.true.
  tpt_std=tpt_frz_pnt
  
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
     if (argv(1:2) == '-f') then
        call getarg(arg+1,argv)
        read (argv,'(f8.3)') float_foo
     endif
     if (argv(1:2) == '-H') then
        NOAA=.true.
        JPL=.not.NOAA
        NCAR=.not.NOAA
     endif
     if (argv(1:2) == '-i') then
        call getarg(arg+1,argv)
        read (argv,'(a)') fl_in
     endif
     if (argv(1:2) == '-J') then
        JPL=.true.
        NCAR=.not.JPL
        NOAA=.not.JPL
     endif
     if (argv(1:2) == '-N') then
        NCAR=.true.
        JPL=.not.NCAR
        NOAA=.not.NCAR
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
     endif
     if (argv(1:2) == '-t') then
        call getarg(arg+1,argv)
        read (argv,'(f8.3)') tpt_std
     endif
     if (argv(1:2) == '-v') then
        write (6,'(a)') CVS_Id
        goto 1000
     endif
     if (argv(1:2) == '-Y') then
        call getarg(arg+1,argv)
        read (argv,'(a)') fl_yld
     endif
  end do
  
  if (JPL) then
     bnd_nbr=bnd_nbr_JPL    ! [nbr] Number of bands
  else if (NCAR) then
     bnd_nbr=bnd_nbr_NCAR   ! [nbr] Number of bands
  else if (NOAA) then
     bnd_nbr=bnd_nbr_NOAA   ! [nbr] Number of bands
  else
     stop 'One of JPL, NCAR, NOAA must be .true.'   
  endif                     ! endif
  
  ! Allocate space for dynamic arrays
  allocate(abs_cff_mss(bnd_nbr),stat=rcd)
  allocate(abs_xsx(bnd_nbr),stat=rcd)
  allocate(abs_xsx_dadT(bnd_nbr),stat=rcd)
  allocate(abs_xsx_tpt_rfr(bnd_nbr),stat=rcd)
  allocate(bnd(bnd_nbr),stat=rcd)     ! coordinate variable
  allocate(flx_bnd_dwn_TOA(bnd_nbr),stat=rcd)
  allocate(flx_bnd_pht_dwn_TOA(bnd_nbr),stat=rcd)
  allocate(flx_slr_frc(bnd_nbr),stat=rcd)
  allocate(flx_spc_dwn_TOA(bnd_nbr),stat=rcd)
  allocate(idx_rfr_air_STP(bnd_nbr),stat=rcd)
  allocate(nrg_pht(bnd_nbr),stat=rcd)
  allocate(qnt_yld(bnd_nbr),stat=rcd)
  allocate(wvl(bnd_nbr),stat=rcd)     ! coordinate variable
  allocate(wvl_ctr(bnd_nbr),stat=rcd)
  allocate(wvl_dlt(bnd_nbr),stat=rcd)
  allocate(wvl_grd(bnd_nbr+1),stat=rcd)
  allocate(wvl_max(bnd_nbr),stat=rcd)
  allocate(wvl_min(bnd_nbr),stat=rcd)
  allocate(xsx_wgt_flx(bnd_nbr),stat=rcd)
  
  ! Compute any quantities that might depend on command line input.
  call ftn_strnul(fl_in)
  call ftn_strnul(fl_out)
  call ftn_strnul(fl_slr)
  call ftn_strnul(fl_yld)
  call ftn_strcpy(src_fl_sng,'Original data file is ' // fl_in)
  if (JPL) call ftn_strcpy(src_rfr_sng,'Data source is Schneider et al. (1987) (JPL)')
  if (NCAR) call ftn_strcpy(src_rfr_sng,'Data source is Davidson et al. (198x) (NCAR) [263.8--648.8 nm (0.5nm), 273K]')
  if (NOAA) call ftn_strcpy(src_rfr_sng,'Data source is Harder, Brault, Johnston and Mount (1997) (NOAA) ' // &
       '[341.5--550.0 nm (0.001nm)]')
  if (NCAR.and.std_tpt) then
     write (6,'(a)') 'WARNING: NCAR data does not yet have temperature dependence, std_tpt ignored'
  endif
  
  open (fl_yld_unit,file=fl_yld,status='old',iostat=rcd)
  if (Cal87) then           ! Cal87 data
     bnd_nbr=bnd_nbr_Cal87
     do idx=1,9
        read (fl_yld_unit,'(a80)') lbl
     enddo
     lbl(1:1)=lbl(1:1) ! CEWI
     do bnd_idx=1,bnd_nbr
        read (fl_yld_unit,*)  &
             wvl_ctr_yyy(bnd_idx), &
             qnt_yld_yyy(bnd_idx)
     enddo
     
     ! Convert input data to SI units where necessary
     do bnd_idx=1,bnd_nbr
        wvl_ctr_yyy(bnd_idx)=wvl_ctr_yyy(bnd_idx)*1.0e-9 ! nm -> m
     enddo                  ! end loop over bnd
     
     ! Sanity check
     do bnd_idx=1,bnd_nbr
        if (qnt_yld_yyy(bnd_idx) < 0.0.or.qnt_yld_yyy(bnd_idx) > 1.0) then
           write (6,'(a,2(a,i4,a,es15.8))')  &
                prg_nm(1:ftn_strlen(prg_nm)), &
                ': ERROR wvl_ctr_yyy(',bnd_idx,') = ',wvl_ctr_yyy(bnd_idx), &
                ', qnt_yld_yyy(',bnd_idx,') = ',qnt_yld_yyy(bnd_idx)
           stop 'EXIT_FAILURE from main()'
        endif               ! endif error
     enddo
     
     ! Compute diagnostic variables
     wvl_min_yyy(1)=wvl_ctr_yyy(1)-0.5*(wvl_ctr_yyy(2)-wvl_ctr_yyy(1))
     do bnd_idx=2,bnd_nbr
        wvl_min_yyy(bnd_idx)=0.5*(wvl_ctr_yyy(bnd_idx-1)+wvl_ctr_yyy(bnd_idx))
     enddo
     do bnd_idx=1,bnd_nbr-1
        wvl_max_yyy(bnd_idx)=0.5*(wvl_ctr_yyy(bnd_idx)+wvl_ctr_yyy(bnd_idx+1))
     enddo
     wvl_max_yyy(bnd_nbr)=wvl_ctr_yyy(bnd_nbr)+0.5*(wvl_ctr_yyy(bnd_nbr)-wvl_ctr_yyy(bnd_nbr-1))
     do bnd_idx=1,bnd_nbr
        wvl_grd_yyy(bnd_idx)=wvl_min_yyy(bnd_idx)
     enddo
     wvl_grd_yyy(bnd_nbr+1)=wvl_max_yyy(bnd_nbr)
     
  endif                     ! Cal87
  
  open (fl_in_unit,file=fl_in,status='old',iostat=rcd)
  if (JPL) then             ! JPL data
     bnd_nbr=bnd_nbr_JPL
     do idx=1,44
        read (fl_in_unit,'(a80)') lbl
     enddo
     do bnd_idx=1,bnd_nbr
        read (fl_in_unit,*)  &
             wvl_min(bnd_idx), &
             wvl_max(bnd_idx), &
             abs_xsx(bnd_idx), &
             abs_xsx_tpt_rfr(bnd_idx), &
             abs_xsx_dadT(bnd_idx)
     enddo
     
     ! Convert input data to SI units where necessary
     do bnd_idx=1,bnd_nbr
        wvl_min(bnd_idx)=wvl_min(bnd_idx)*1.0e-9 ! nm -> m
        wvl_max(bnd_idx)=wvl_max(bnd_idx)*1.0e-9 ! nm -> m
        abs_xsx(bnd_idx)=abs_xsx(bnd_idx)*1.0e-24 ! 1.0e20 cm2 mlc-1 -> m2 mlc-1
        abs_xsx_tpt_rfr(bnd_idx)=abs_xsx_tpt_rfr(bnd_idx)+tpt_frz_pnt ! C -> K
        abs_xsx_dadT(bnd_idx)=abs_xsx_dadT(bnd_idx)*1.0e-26 ! 1.0e22 cm2 mlc-1 deg-1 -> m2 mlc-1 deg-1 
     enddo                  ! end loop over bnd
     
     ! Compute diagnostic variables.  
     do bnd_idx=1,bnd_nbr
        bnd(bnd_idx)=0.5*(wvl_min(bnd_idx)+wvl_max(bnd_idx))
        wvl_ctr(bnd_idx)=bnd(bnd_idx)
     enddo                  ! end loop over bnd
     
     ! Recompute and save cross sections at a standard temperature 
     ! Currently, only JPL supplies temperature dependence
     if (std_tpt) then
        ! a(T)=a(T0)+(T-T0)*dadT
        do bnd_idx=1,bnd_nbr
           abs_xsx(bnd_idx)=abs_xsx(bnd_idx)+abs_xsx_dadT(bnd_idx)*(tpt_std-abs_xsx_tpt_rfr(bnd_idx))
           abs_xsx_tpt_rfr(bnd_idx)=tpt_std ! Stamp data with new temperature
        enddo
        write (6,'(a30,1x,f7.3,a2)') 'Renormalized cross sections to',tpt_std,' K'
     endif
     
  else if (NCAR) then       ! NCAR data
     
     bnd_nbr=bnd_nbr_NCAR
     do bnd_idx=1,bnd_nbr
        read (fl_in_unit,*)  &
             wvl_ctr(bnd_idx), &
             abs_xsx(bnd_idx), &
             float_foo,     & ! NO2 cross-section at unknown temperature
             float_foo      ! NO2 cross-section at unknown temperature
     enddo
     float_foo=float_foo ! CEWI
     
     ! Convert input data to SI units where necessary
     do bnd_idx=1,bnd_nbr
        wvl_ctr(bnd_idx)=wvl_ctr(bnd_idx)*1.0e-9 ! nm -> m
        abs_xsx(bnd_idx)=abs_xsx(bnd_idx)*1.0e-4 ! cm2 mlc-1 -> m2 mlc-1
     enddo
     
     ! Compute diagnostic variables
     wvl_min(1)=wvl_ctr(1)-0.5*(wvl_ctr(2)-wvl_ctr(1))
     do bnd_idx=2,bnd_nbr
        wvl_min(bnd_idx)=0.5*(wvl_ctr(bnd_idx-1)+wvl_ctr(bnd_idx))
     enddo
     do bnd_idx=1,bnd_nbr-1
        wvl_max(bnd_idx)=0.5*(wvl_ctr(bnd_idx)+wvl_ctr(bnd_idx+1))
     enddo
     wvl_max(bnd_nbr)=wvl_ctr(bnd_nbr)+0.5*(wvl_ctr(bnd_nbr)-wvl_ctr(bnd_nbr-1))
     
     do bnd_idx=1,bnd_nbr
        bnd(bnd_idx)=wvl_ctr(bnd_idx)
        wvl(bnd_idx)=wvl_ctr(bnd_idx)
        abs_xsx_tpt_rfr(bnd_idx)=tpt_frz_pnt ! all NCAR data taken at 0 C
        abs_xsx_dadT(bnd_idx)=0.0 ! no temperature dependence given for NCAR data
     enddo
     
  else if (NOAA) then       ! NOAA data
     
     bnd_nbr=bnd_nbr_NOAA
     do idx=1,10
        read (fl_in_unit,'(a80)') lbl
     enddo
     do bnd_idx=1,bnd_nbr
        read (fl_in_unit,*)  &
             abs_xsx(bnd_idx)
        wvl_ctr(bnd_idx)=345.0+(bnd_idx-1)*0.001
     enddo
     
     ! Convert input data to SI units where necessary.
     do bnd_idx=1,bnd_nbr
        wvl_ctr(bnd_idx)=wvl_ctr(bnd_idx)*1.0e-9 ! nm -> m
        abs_xsx(bnd_idx)=abs_xsx(bnd_idx)*1.0e-19*1.0e-4 ! 1.0e-19 cm2 mlc-1 -> m2 mlc-1
     enddo
     
     ! Compute diagnostic variables.  
     wvl_min(1)=wvl_ctr(1)-0.5*(wvl_ctr(2)-wvl_ctr(1))
     do bnd_idx=2,bnd_nbr
        wvl_min(bnd_idx)=0.5*(wvl_ctr(bnd_idx-1)+wvl_ctr(bnd_idx))
     enddo
     do bnd_idx=1,bnd_nbr-1
        wvl_max(bnd_idx)=0.5*(wvl_ctr(bnd_idx)+wvl_ctr(bnd_idx+1))
     enddo
     wvl_max(bnd_nbr)=wvl_ctr(bnd_nbr)+0.5*(wvl_ctr(bnd_nbr)-wvl_ctr(bnd_nbr-1))
     
     ! Nominal temperature of data is currently taken from input filename.
     do bnd_idx=1,bnd_nbr
        bnd(bnd_idx)=wvl_ctr(bnd_idx)
        ! NOAA data taken at uniform temperature
        abs_xsx_tpt_rfr(bnd_idx)=tpt_std 
        ! Temperature dependence for NOAA data not currently implemented
        abs_xsx_dadT(bnd_idx)=0.0
     enddo
     
  endif                     ! NOAA data
  
  close (fl_in_unit)
  write (6,'(a20,1x,a)') 'Read input data from',fl_in(1:ftn_strlen(fl_in))
  
  ! Get TOA solar spectrum
  slr_spc_xtr_typ=xtr_fll_ngh+xtr_vrb
  call slr_spc_get(fl_slr,wvl_min,wvl_max,bnd_nbr,flx_slr_frc,slr_spc_xtr_typ,slr_spc_xtr_typ)
  
  ! Compute diagnostic variables
  do bnd_idx=1,bnd_nbr
     abs_cff_mss(bnd_idx)=abs_xsx(bnd_idx)*Avagadro/mmw_NO2
     wvl_dlt(bnd_idx)=wvl_max(bnd_idx)-wvl_min(bnd_idx)
     wvl_grd(bnd_idx)=wvl_min(bnd_idx)
     nrg_pht(bnd_idx)=Planck*speed_of_light/wvl_ctr(bnd_idx)
     flx_bnd_dwn_TOA(bnd_idx)=flx_slr_frc(bnd_idx)*slr_cst_CCM
     flx_spc_dwn_TOA(bnd_idx)=flx_slr_frc(bnd_idx)*slr_cst_CCM/wvl_dlt(bnd_idx)
     flx_bnd_pht_dwn_TOA(bnd_idx)=flx_bnd_dwn_TOA(bnd_idx)/nrg_pht(bnd_idx)
  enddo                     ! end loop over bnd
  wvl_grd(bnd_nbr+1)=wvl_max(bnd_nbr)
  
  ! Compute index of refraction through dry air at STP (Len93 p. 155)
  do bnd_idx=1,bnd_nbr
     idx_rfr_air_STP(bnd_idx)= &
          1.0+ &
          1.0e-6*(77.46+0.459/(1.0e12*wvl_ctr(bnd_idx)**2))* &
          prs_STP*0.01/tpt_STP
  enddo
  
  ! Rebin quantum yield using nearest neighbor extrapolation
  xtr_typ_LHS=xtr_prt_ngh+xtr_fll_ngh+xtr_vrb
  xtr_typ_RHS=xtr_prt_ngh+xtr_fll_ngh+xtr_vrb
  call rbn_vec(bnd_nbr_Cal87,wvl_grd_yyy,qnt_yld_yyy, &
       bnd_nbr,wvl_grd,qnt_yld, &
       xtr_typ_LHS,xtr_typ_RHS)
  ! Sanity check
  do bnd_idx=1,bnd_nbr
     if (qnt_yld(bnd_idx) < 0.0.or.qnt_yld(bnd_idx) > 1.0) then
        write (6,'(a,2(a,i4,a,es15.8))')  &
             prg_nm(1:ftn_strlen(prg_nm)), &
             ': WARNING bnd(',bnd_idx,') = ',bnd(bnd_idx), &
             ', qnt_yld(',bnd_idx,') = ',qnt_yld(bnd_idx)
        !        stop 'EXIT_FAILURE from main()'
        if (qnt_yld(bnd_idx) < 0.0) qnt_yld(bnd_idx)=0.0
        if (qnt_yld(bnd_idx) > 1.0) qnt_yld(bnd_idx)=1.0
        write (6,'(a,1(a,i4,a,es15.8))')  &
             prg_nm(1:ftn_strlen(prg_nm)), &
             ': WARNING Reset value to qnt_yld(',bnd_idx,') = ',qnt_yld(bnd_idx)
     endif                  ! endif error
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
  rcd=rcd+nf90_def_var(nc_id,'abs_cff_mss',nf90_float,bnd_dim_id,abs_cff_mss_id)
  rcd=rcd+nf90_def_var(nc_id,'abs_xsx',nf90_float,bnd_dim_id,abs_xsx_id)
  rcd=rcd+nf90_def_var(nc_id,'abs_xsx_dadT',nf90_float,bnd_dim_id,abs_xsx_dadT_id)
  rcd=rcd+nf90_def_var(nc_id,'abs_xsx_tpt_rfr',nf90_float,bnd_dim_id,abs_xsx_tpt_rfr_id)
  rcd=rcd+nf90_def_var(nc_id,'bnd',nf90_float,bnd_dim_id,bnd_id)
  rcd=rcd+nf90_def_var(nc_id,'flx_bnd_dwn_TOA',nf90_float,bnd_dim_id,flx_bnd_dwn_TOA_id)
  rcd=rcd+nf90_def_var(nc_id,'flx_bnd_pht_dwn_TOA',nf90_float,bnd_dim_id,flx_bnd_pht_dwn_TOA_id)
  rcd=rcd+nf90_def_var(nc_id,'flx_slr_frc',nf90_float,bnd_dim_id,flx_slr_frc_id)
  rcd=rcd+nf90_def_var(nc_id,'flx_spc_dwn_TOA',nf90_float,bnd_dim_id,flx_spc_dwn_TOA_id)
  rcd=rcd+nf90_def_var(nc_id,'idx_rfr_air_STP',nf90_float,bnd_dim_id,idx_rfr_air_STP_id)
  rcd=rcd+nf90_def_var(nc_id,'nrg_pht',nf90_float,bnd_dim_id,nrg_pht_id)
  rcd=rcd+nf90_def_var(nc_id,'qnt_yld',nf90_float,bnd_dim_id,qnt_yld_id)
  rcd=rcd+nf90_def_var(nc_id,'wvl_ctr',nf90_float,bnd_dim_id,wvl_ctr_id)
  rcd=rcd+nf90_def_var(nc_id,'wvl_grd',nf90_float,grd_dim_id,wvl_grd_id)
  rcd=rcd+nf90_def_var(nc_id,'wvl_max',nf90_float,bnd_dim_id,wvl_max_id)
  rcd=rcd+nf90_def_var(nc_id,'wvl_min',nf90_float,bnd_dim_id,wvl_min_id)
  rcd=rcd+nf90_def_var(nc_id,'wvl_dlt',nf90_float,bnd_dim_id,wvl_dlt_id)
  
  ! Add global attributes
  rcd=rcd+nf90_put_att(nc_id,nf90_global,'CVS_Id',CVS_Id)
  rcd=rcd+nf90_put_att(nc_id,nf90_global,'creation_date',lcl_date_time)
  rcd=rcd+nf90_put_att(nc_id,nf90_global,'prg_ID',prg_ID(1:ftn_strlen(prg_ID)))
  rcd=rcd+nf90_put_att(nc_id,nf90_global,'cmd_ln',cmd_ln(1:ftn_strlen(cmd_ln)))
  rcd=rcd+nf90_put_att(nc_id,nf90_global,'src_rfr_sng',src_rfr_sng(1:ftn_strlen(src_rfr_sng)))
  rcd=rcd+nf90_put_att(nc_id,nf90_global,'src_fl_sng',src_fl_sng(1:ftn_strlen(src_fl_sng)))
  ! Add english text descriptions
  rcd=rcd+nf90_put_att(nc_id,abs_cff_mss_id,'long_name','Nitrogen Dioxide mass absorption coefficient')
  rcd=rcd+nf90_put_att(nc_id,abs_xsx_dadT_id,'long_name','Slope of absorption cross section temperature dependence')
  rcd=rcd+nf90_put_att(nc_id,abs_xsx_id,'long_name','Nitrogen Dioxide continuum absorption cross section')
  rcd=rcd+nf90_put_att(nc_id,abs_xsx_tpt_rfr_id,'long_name','Valid temperature for absorption cross section')
  rcd=rcd+nf90_put_att(nc_id,bnd_id,'long_name','Band nominal wavelength')
  rcd=rcd+nf90_put_att(nc_id,flx_bnd_dwn_TOA_id,'long_name','Solar Energy flux in band')
  rcd=rcd+nf90_put_att(nc_id,flx_bnd_pht_dwn_TOA_id,'long_name','Photon flux in band')
  rcd=rcd+nf90_put_att(nc_id,flx_slr_frc_id,'long_name','Fraction of solar flux in band: ' // fl_slr)
  rcd=rcd+nf90_put_att(nc_id,flx_spc_dwn_TOA_id,'long_name','Spectral solar insolation at TOA')
  rcd=rcd+nf90_put_att(nc_id,idx_rfr_air_STP_id,'long_name','Index of refraction at band center at STP')
  rcd=rcd+nf90_put_att(nc_id,nrg_pht_id,'long_name','Energy of photon at band center')
  rcd=rcd+nf90_put_att(nc_id,qnt_yld_id,'long_name','Quantum yield for NO2 + hv --> O + NO')
  rcd=rcd+nf90_put_att(nc_id,wvl_ctr_id,'long_name','Band center wavelength')
  rcd=rcd+nf90_put_att(nc_id,wvl_grd_id,'long_name','Wavelength grid')
  rcd=rcd+nf90_put_att(nc_id,wvl_max_id,'long_name','Band maximum wavelength')
  rcd=rcd+nf90_put_att(nc_id,wvl_min_id,'long_name','Band minimum wavelength')
  rcd=rcd+nf90_put_att(nc_id,wvl_dlt_id,'long_name','Bandwidth')
  
  ! Add units
  rcd=rcd+nf90_put_att(nc_id,abs_cff_mss_id,'units','meter2 kilogram-1')
  rcd=rcd+nf90_put_att(nc_id,abs_xsx_dadT_id,'units','meter2 kelvin-1')
  rcd=rcd+nf90_put_att(nc_id,abs_xsx_id,'units','meter2')
  rcd=rcd+nf90_put_att(nc_id,abs_xsx_tpt_rfr_id,'units','kelvin')
  rcd=rcd+nf90_put_att(nc_id,bnd_id,'units','meter')
  rcd=rcd+nf90_put_att(nc_id,flx_bnd_dwn_TOA_id,'units','watt meter-2')
  rcd=rcd+nf90_put_att(nc_id,flx_bnd_pht_dwn_TOA_id,'units','photon meter-2 second-1')
  rcd=rcd+nf90_put_att(nc_id,flx_slr_frc_id,'units','fraction')
  rcd=rcd+nf90_put_att(nc_id,flx_spc_dwn_TOA_id,'units','watt meter-2 meter-1')
  rcd=rcd+nf90_put_att(nc_id,idx_rfr_air_STP_id,'units','fraction')
  rcd=rcd+nf90_put_att(nc_id,nrg_pht_id,'units','joule photon-1')
  rcd=rcd+nf90_put_att(nc_id,qnt_yld_id,'units','fraction')
  rcd=rcd+nf90_put_att(nc_id,wvl_ctr_id,'units','meter')
  rcd=rcd+nf90_put_att(nc_id,wvl_grd_id,'units','meter')
  rcd=rcd+nf90_put_att(nc_id,wvl_max_id,'units','meter')
  rcd=rcd+nf90_put_att(nc_id,wvl_min_id,'units','meter')
  rcd=rcd+nf90_put_att(nc_id,wvl_dlt_id,'units','meter')
  
  ! All dimensions, variables, and attributes have been defined
  rcd=rcd+nf90_enddef(nc_id)
  
  ! Write data
  rcd=rcd+nf90_put_var(nc_id,abs_cff_mss_id,abs_cff_mss)
  rcd=rcd+nf90_put_var(nc_id,abs_xsx_dadT_id,abs_xsx_dadT)
  rcd=rcd+nf90_put_var(nc_id,abs_xsx_id,abs_xsx)
  rcd=rcd+nf90_put_var(nc_id,abs_xsx_tpt_rfr_id,abs_xsx_tpt_rfr)
  rcd=rcd+nf90_put_var(nc_id,bnd_id,bnd)
  rcd=rcd+nf90_put_var(nc_id,flx_bnd_dwn_TOA_id,flx_bnd_dwn_TOA)
  rcd=rcd+nf90_put_var(nc_id,flx_bnd_pht_dwn_TOA_id,flx_bnd_pht_dwn_TOA)
  rcd=rcd+nf90_put_var(nc_id,flx_slr_frc_id,flx_slr_frc)
  rcd=rcd+nf90_put_var(nc_id,flx_spc_dwn_TOA_id,flx_spc_dwn_TOA)
  rcd=rcd+nf90_put_var(nc_id,idx_rfr_air_STP_id,idx_rfr_air_STP)
  rcd=rcd+nf90_put_var(nc_id,nrg_pht_id,nrg_pht)
  rcd=rcd+nf90_put_var(nc_id,qnt_yld_id,qnt_yld)
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
     xsx_wgt_flx(idx)=abs_xsx(idx)*flx_spc_dwn_TOA(idx)
  enddo                     ! end loop over NO2 bnd
  
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
     abs_xsx_CCM(idx)=xsx_wgt_flx_CCM(idx)/flx_spc_dwn_TOA_CCM(idx)
     abs_cff_mss_CCM(idx)=abs_xsx_CCM(idx)*Avagadro/mmw_NO2
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
  
  ! Add NO2 data on CCM grid to netCDF output file
  rcd=rcd+nf90_wrp_open(fl_out,nf90_write,nc_id)
  if (rcd /= nf90_noerr) call nf90_err_exit(rcd,fl_out)
  
  ! Put output file in define mode
  rcd=rcd+nf90_redef(nc_id)
  if (rcd /= nf90_noerr) call nf90_err_exit(rcd,fl_out)
  
  ! Variable definitions
  rcd=rcd+nf90_def_var(nc_id,'abs_cff_mss_CCM',nf90_float,bnd_dim_id,abs_cff_mss_id)
  rcd=rcd+nf90_def_var(nc_id,'abs_xsx_CCM',nf90_float,bnd_dim_id,abs_xsx_id)
  rcd=rcd+nf90_def_var(nc_id,'flx_slr_frc_CCM',nf90_float,bnd_dim_id,flx_slr_frc_id)
  rcd=rcd+nf90_def_var(nc_id,'flx_spc_dwn_TOA_CCM',nf90_float,bnd_dim_id,flx_spc_dwn_TOA_id)
  
  ! Add global attributes
  
  ! Add english text descriptions
  rcd=rcd+nf90_put_att(nc_id,abs_cff_mss_id,'long_name','Nitrogen Dioxide mass absorption coefficient')
  rcd=rcd+nf90_put_att(nc_id,abs_xsx_id,'long_name','Nitrogen Dioxide continuum absorption cross section')
  rcd=rcd+nf90_put_att(nc_id,flx_slr_frc_id,'long_name','Fraction of solar flux in band: ' // fl_slr)
  rcd=rcd+nf90_put_att(nc_id,flx_spc_dwn_TOA_id,'long_name','Spectral solar insolation at TOA')
  
  ! Add units
  rcd=rcd+nf90_put_att(nc_id,abs_cff_mss_id,'units','meter2 kilogram-1')
  rcd=rcd+nf90_put_att(nc_id,abs_xsx_id,'units','meter2')
  rcd=rcd+nf90_put_att(nc_id,flx_slr_frc_id,'units','fraction')
  rcd=rcd+nf90_put_att(nc_id,flx_spc_dwn_TOA_id,'units','watt meter-2 meter-1')
  
  ! All dimensions, variables, and attributes have been defined
  rcd=rcd+nf90_enddef(nc_id)
  if (rcd /= nf90_noerr) call nf90_err_exit(rcd,fl_out)
  
  ! Write out data
  rcd=rcd+nf90_put_var(nc_id,abs_cff_mss_id,abs_cff_mss_CCM)
  rcd=rcd+nf90_put_var(nc_id,abs_xsx_id,abs_xsx_CCM)
  rcd=rcd+nf90_put_var(nc_id,flx_slr_frc_id,flx_slr_frc_CCM)
  rcd=rcd+nf90_put_var(nc_id,flx_spc_dwn_TOA_id,flx_spc_dwn_TOA_CCM)
  
  ! Close output file
  rcd=rcd+nf90_close(nc_id)
  if (rcd /= nf90_noerr) call nf90_err_exit(rcd,fl_out)
  write (6,'(a29,1x,a)') 'Wrote NO2 data on CCM grid to',fl_out(1:ftn_strlen(fl_out))
  
  ! Convert absorption coefficients to cm2/gm and output block for radcsw
  do idx=1,bnd_nbr_CCM
     abs_cff_mss_CCM_CGS(idx)=abs_cff_mss_CCM(idx)*10.0 !  m2 kg-1 --> [cm2 g-1] 
  enddo                     ! end loop over CCM bnd
  nbr_dat_per_ln=4
  write (0,'(a,a,a,a)') 'c     Data generated by ',prg_nm(1:ftn_strlen(prg_nm)),' on ',lcl_date_time
  write (0,'(a,a)') 'c     ',prg_ID(1:ftn_strlen(prg_ID))
  write (0,'(a,a)') 'c     Command line: ',cmd_ln(1:ftn_strlen(cmd_ln))
  write (0,'(a,f7.3,a)') 'c     NO2 mass absorption coefficients in [cm2 g-1] for T = ',tpt_std,' K'
  write (0,'(a,a)') 'c     ',src_rfr_sng(1:ftn_strlen(src_rfr_sng))
  write (0,'(a,a)') 'c     ',src_fl_sng(1:ftn_strlen(src_fl_sng))
  call dat_prn_f77(bnd_nbr_CCM,abs_cff_mss_CCM_CGS,nbr_dat_per_ln,'abNO2'//char(0))
  call dat_prn_f90(bnd_nbr_CCM,abs_cff_mss_CCM_CGS,nbr_dat_per_ln,'abNO2'//char(0))
  
  ! De-allocate dynamic variables
  if (allocated(abs_cff_mss)) deallocate(abs_cff_mss,stat=rcd)
  if (allocated(abs_xsx)) deallocate(abs_xsx,stat=rcd)
  if (allocated(abs_xsx_dadT)) deallocate(abs_xsx_dadT,stat=rcd)
  if (allocated(abs_xsx_tpt_rfr)) deallocate(abs_xsx_tpt_rfr,stat=rcd)
  if (allocated(bnd)) deallocate(bnd,stat=rcd)     ! coordinate variable
  if (allocated(flx_bnd_dwn_TOA)) deallocate(flx_bnd_dwn_TOA,stat=rcd)
  if (allocated(flx_bnd_pht_dwn_TOA)) deallocate(flx_bnd_pht_dwn_TOA,stat=rcd)
  if (allocated(flx_slr_frc)) deallocate(flx_slr_frc,stat=rcd)
  if (allocated(flx_spc_dwn_TOA)) deallocate(flx_spc_dwn_TOA,stat=rcd)
  if (allocated(idx_rfr_air_STP)) deallocate(idx_rfr_air_STP,stat=rcd)
  if (allocated(nrg_pht)) deallocate(nrg_pht,stat=rcd)
  if (allocated(qnt_yld)) deallocate(qnt_yld,stat=rcd)
  if (allocated(wvl)) deallocate(wvl,stat=rcd)     ! coordinate variable
  if (allocated(wvl_ctr)) deallocate(wvl_ctr,stat=rcd)
  if (allocated(wvl_dlt)) deallocate(wvl_dlt,stat=rcd)
  if (allocated(wvl_grd)) deallocate(wvl_grd,stat=rcd)
  if (allocated(wvl_max)) deallocate(wvl_max,stat=rcd)
  if (allocated(wvl_min)) deallocate(wvl_min,stat=rcd)
  if (allocated(xsx_wgt_flx)) deallocate(xsx_wgt_flx,stat=rcd)
  
1000 continue
  
  call exit(exit_status)
end program NO2

