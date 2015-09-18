! $Id$

program O2X
  
  ! Purpose: Convert O2-O2 absorption cross section data to netCDF format
  ! NB: O2X should be compiled/run in double precision
  
  ! CCMCG stands for CCM_CGS, i.e., quantity is in CGS units not SI
  ! Abbreviation is required to keep symbol names <= 31 characters in some Fortran compilers
  
  ! Compilation:
  ! cd ${HOME}/aca; make -W O2X.F OPTS=D O2X; cd -
  ! cd ${HOME}/aca; make -W O2X.F O2X; cd -
  ! cd ${HOME}/aca; make OPTS=D O2X; cd -
  
  ! Usage:
  ! ncks -H -C -F -v wvl_ctr_CCM,flx_slr_frc_CCM,abs_xsx_O2O2_CCM ${DATA}/aca/abs_xsx_O2O2.nc
  ! ncks -H -C -F -d bnd,1.26e-6 -d bnd_CCM,1.26 -v wvl_ctr_CCM,flx_slr_frc_CCM,abs_xsx_O2O2_CCM,abs_xsx_O2O2 ${DATA}/aca/abs_xsx_O2O2.nc
  ! ncwa -a bnd -d bnd,0.295e-6,0.305e-6 ${DATA}/aca/abs_xsx_O2O2.nc ${DATA}/aca/O2O2.nc
  ! ncks -H -C -F -d bnd,0.3e-6 -v abs_xsx_O2O2 ${DATA}/aca/abs_xsx_O2O2.nc
  ! ncks -H -C -F -d bnd,1.26e-6 -v odal_O2O2_PUMC_O2_PUMP_O2,odal_O2N2_PUMC_O2_PUMP_N2 ${DATA}/aca/abs_xsx_O2O2.nc
  ! ncks -H -C -F -d bnd,1.26e-6 -v odal_O2O2_PUMC_O2_PUMP_O2,odal_O2O2_PUNC_O2_PUNP_O2,odal_O2O2_PUMC_O2_PUMP_O2_CCM ${DATA}/aca/abs_xsx_O2O2.nc
  ! ncks -H -C -F -v abs_xsx_O2O2 ${DATA}/aca/O2O2.nc
  ! ncks -H -C -F -d bnd_CCM,6 -v abs_xsx_O2O2_CCM ${DATA}/aca/abs_xsx_O2O2.nc
  
  ! Default
  ! O2X 
  ! O2X -i ${DATA}/aca/abs_xsx_O2O2.txt -o ${DATA}/aca/abs_xsx_O2O2.nc
  
  ! Exclude 1.26 micron band
  ! O2X -x -i ${DATA}/aca/abs_xsx_O2O2.txt -o ${DATA}/aca/abs_xsx_O2O2_xcl_1260nm.nc
  
  ! NB: data were measured and are recorded as absorption cross section excess
  ! (relative to the absorption cross section of O2) per unit concentration of O2.
  ! Thus typical numbers are O(1.0e-46), and must be handled in double precision.
  
  ! Currently the code processes input ASCII data files that look like:
  
  !#Greenblatt, et al. O2-O2 absorption (ref: JGR, 95, 18577, 1990)
  !#3182 points (note: 335.2-663.4 nm & 1005.0-1136.8 nm)
  !#Column 1: wavelength in nm
  !#column 2: "O4" absorption in cm^5 mol^-2
  !335.2   2.079e-47
  !335.3   1.975e-47
  !...
  
  use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
  use netcdf ! [mdl] netCDF interface
  use nf90_utl ! [mdl] netCDF utilities
  use phys_cst_mdl ! [mdl] Fundamental and derived physical constants
  use prn_mdl,only:dat_prn_f77,dat_prn_f90 ! [mdl] Print formatting & utilities
  use rt_mdl ! [mdl] Radiative transfer utilities
  use sng_mdl ! [mdl] String manipulation
  use utl_mdl ! [mdl] Utility functions (date_time_get,mnt_chk...)
  use wvl_mdl ! [mdl] Wavelength grid parameters
  use xtr_mdl ! [mdl] Extrapolation/interpolation handling
  
  implicit none
  ! Parameters
  character*(*) CVS_Id
  
  integer fl_in_unit
  integer bnd_nbr_GOB90
  integer bnd_nbr_max
  parameter(fl_in_unit=73, &
       bnd_nbr_GOB90=4086, &    
       bnd_nbr_max=bnd_nbr_GOB90, & ! max(bnd_nbr_GOB90)
       CVS_Id='$Id$')
  
  real,parameter::mlt_fct=1.0e20 ! [frc]

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
  character src_wgt_sng*300
  character wgt_nm*80
  
  integer arg
  integer arg_nbr
  integer exit_status       ! program exit status
  integer idx
  integer rcd               ! return success code
  
  logical GOB90
  logical flg_ncl_1260nm_bnd
  logical flg_ncl_1530nm_bnd
  logical std_tpt
  logical WGT_TRN
  
  integer bnd_dim_id        ! dimension ID for bands
  integer grd_dim_id        ! dimension ID for grid
  integer bnd_idx           ! counting index
  integer nc_id             ! file handle
  integer bnd_nbr           ! dimension size
  
  integer abs_cff_mss_O2O2_id
  integer abs_xsx_O2N2_id
  integer abs_xsx_O2O2_dadT_id
  integer abs_xsx_O2O2_id
  integer abs_xsx_O2O2_tpt_rfr_id
  integer bnd_id            ! coordinate ID
  integer flx_bnd_dwn_TOA_id
  integer flx_bnd_pht_dwn_TOA_id
  integer flx_slr_frc_id
  integer flx_spc_dwn_TOA_id
  integer idx_rfr_air_STP_id
  integer nrg_pht_id
  integer odal_O2N2_PUMC_O2_PUMP_N2_id
  integer odal_O2N2_PUNC_O2_PUNP_N2_id
  integer odal_O2O2_PUMC_O2_PUMP_O2_id
  integer odal_O2O2_PUNC_O2_PUNP_O2_id
  integer wvl_ctr_id
  integer wvl_grd_id
  integer wvl_max_id
  integer wvl_min_id
  integer wvl_dlt_id
  
  real N2_cll_prt_fsh       ! Efficiency of N2 (relative to O2) as a collision partner for O2 in the 1.26 micron band
  real tpt_std
  real(selected_real_kind(p=12))::double_foo
  
  ! Allocatable variables
  real(selected_real_kind(p=12)),dimension(:),allocatable::odal_O2O2_PUNC_O2_PUNP_O2
  real(selected_real_kind(p=12)),dimension(:),allocatable::odal_O2N2_PUNC_O2_PUNP_N2
  real,dimension(:),allocatable::abs_cff_mss_O2O2
  real,dimension(:),allocatable::abs_xsx_O2N2
  real(selected_real_kind(p=12)),dimension(:),allocatable::abs_xsx_O2O2
  real,dimension(:),allocatable::abs_xsx_O2O2_dadT
  real,dimension(:),allocatable::abs_xsx_O2O2_tpt_rfr
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
  real,dimension(:),allocatable::xsx_wgt_flx
  
  integer bnd_nbr_CCM
  integer nbr_dat_per_ln
  integer slr_spc_xtr_typ
  integer xtr_typ_LHS
  integer xtr_typ_RHS
  
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
  GOB90=.true.
  N2_cll_prt_fsh=0.2
  WGT_TRN=.false.
  dbg_lvl=dbg_off
  exit_status=0
  fl_in='/data/zender/aca/abs_xsx_O2O2.txt'//char(0)
  fl_out='/data/zender/aca/abs_xsx_O2O2.nc'//char(0)
  fl_slr='/data/zender/aca/spc_Kur95_01wvn.nc'//char(0)
  fl_wgt='/data/zender/tmp/swnb_trp.nc'//char(0)
  flg_ncl_1260nm_bnd=.true.
  flg_ncl_1530nm_bnd=.false.
  rcd=nf90_noerr              ! nf90_noerr == 0
  CVS_Date='$Date$'
  CVS_Revision='$Revision$'
  std_tpt=.true.
  tpt_std=296.0
  wgt_nm='flx_spc_dwn_sfc'//char(0)
  
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
        read (argv,'(a)') wgt_nm
     endif
     if (argv(1:2) == '-G') then
        GOB90=.true.
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
     endif
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
     if (argv(1:2) == '-x') then
        flg_ncl_1260nm_bnd=.not.flg_ncl_1260nm_bnd
     endif
     if (argv(1:2) == '-X') then
        flg_ncl_1530nm_bnd=.not.flg_ncl_1530nm_bnd
     endif
  end do
  
  if (GOB90) then
     bnd_nbr=bnd_nbr_GOB90
  else
     stop 'GOB90 is .false. but is only O2-O2 data source'
  endif                     ! endif GOB90
  
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
  allocate(xsx_wgt_flx(bnd_nbr),stat=rcd)
  
  ! Compute any quantities that might depend on command line input.
  call ftn_strnul(fl_in)
  call ftn_strnul(fl_out)
  call ftn_strnul(fl_slr)
  call ftn_strnul(fl_wgt)
  call ftn_strnul(wgt_nm)
  call ftn_strcpy(src_fl_sng,'Original data file is ' // fl_in)
  if (WGT_TRN) then
     call ftn_strcpy(src_wgt_sng,'CCM cross sections are averages of high resolution cross-sections ' &
             //'weighted by variable '//wgt_nm(1:ftn_strlen(wgt_nm))// &
             ' from model atmosphere in '//fl_wgt(1:ftn_strlen(fl_wgt))//char(0))
  else
     call ftn_strcpy(src_wgt_sng,'CCM cross sections are averages of high resolution cross-sections ' &
             //'weighted by TOA solar spectral flux'//char(0))
  endif                     ! endif WGT_TRN
  if (GOB90) call ftn_strcpy(src_rfr_sng,'Data source is Greenblatt et al. (1990) (GOB90)'//char(0))
  if (flg_ncl_1260nm_bnd) then
     write (6,'(a)') 'Including 1.26 micron absorption band'
  else
     write (6,'(a)') 'Excluding 1.26 micron absorption band'
  endif                     ! endif flg_ncl_1260nm_bnd
  if (flg_ncl_1530nm_bnd) then
     write (6,'(a)') 'Including 1.53 micron absorption band seen by MCB98'
     stop 'ERROR: 1.53 micron absorption band not supported yet'
  else
     write (6,'(a)') 'Excluding 1.53 micron absorption band seen by MCB98'
  endif                     ! endif flg_ncl_1530nm_bnd
  
  open (fl_in_unit,file=fl_in,status='old',iostat=rcd)
  
  if (GOB90) then
     do idx=1,7
        read (fl_in_unit,'(a80)') lbl
     enddo                  ! end loop over idx
     lbl(1:1)=lbl(1:1) ! CEWI
     do bnd_idx=1,bnd_nbr
        read (fl_in_unit,*) &
                   wvl_ctr(bnd_idx), &
                   odal_O2O2_PUNC_O2_PUNP_O2(bnd_idx)
     enddo                  ! end loop over bnd
     
     ! Convert input data to SI units where necessary
     do bnd_idx=1,bnd_nbr
        wvl_ctr(bnd_idx)=wvl_ctr(bnd_idx)*1.0e-9 ! nm -> m
        odal_O2O2_PUNC_O2_PUNP_O2(bnd_idx)=odal_O2O2_PUNC_O2_PUNP_O2(bnd_idx)*1.0e-10 ! cm5 mlc-2 -> m5 mlc-2
     enddo                  ! end loop over bnd
     
     ! Zero the 1260 nm band only if requested
     if (.not.flg_ncl_1260nm_bnd) then
        do bnd_idx=1,bnd_nbr
           ! The 1.26 micron band is currently the only O2-O2 absorption beyond 1.137 microns
           if (wvl_ctr(bnd_idx) >= 1.137e-6.and.wvl_ctr(bnd_idx) <= 1.350e-6) then
              odal_O2O2_PUNC_O2_PUNP_O2(bnd_idx)=0.0
           endif            ! endif wvl
        enddo               ! end loop over bnd
     endif                  ! endif flg_ncl_1260nm_bnd
     
     ! NB: The measurement of O2-O2 cross sections leads to non-standard nomenclature.
     
     ! Shardanand (1969) is convinced O2-O2 absorption is attributable to covalently bonded O4. 
     ! Shardanand (1977) is convinced O2-O2, O2-N2, O2-Ar absorption is attributable to van der Waals complexes.
     ! GOB90 and Blake and McCoy (1987) believe O2-O2 absorption is attributable to collision-induced absorption.
     ! Short lived collision pairs, not bound dimers, are conclusively shown to be the cause of O2-O2 and O2-N2 absorption in SPS98.
     ! Thus concentration of O2-O2 the concentration of O2-O2 collision pairs.
     ! Unfortunately, no one knows how to directly calculate the concentration of O2-O2.
     ! The physical meaningfulness of such a concentration is not clear in any case.
     ! Labs measure total monochromatic absorption cross section at given O2 pressure.
     ! Collision-pair induced absorption is inferred from pressure variation of this absorption.
     ! O2 concentration is increased by increasing O2 partial pressure.
     ! Change in measured absorption cross section is then correlated with square of O2 concentration.
     ! GOB90 show once absorption that varies linearly with O2 concentration is removed, residual absorption correlates nearly exactly with square of O2 concentration.
     ! This residual absorption is due to O2-O2 collision complexes.
     ! Unfortunately, only convolution of O2-O2 concentration and O2-O2 absorption is known with certainty.
     ! O2-O2 concentration and O2-O2 absorption cross sections are not known individually.
     ! How to implement this absorption in RT codes?
     
     ! abs_xsx_O2O2: the standard measurement of molecular cross-section. 
     ! Multiply abs_xsx_O2O2 by column path of O2-O2 complexes (# m-2) to obtain absorption optical depth of O2-O2.
     ! abs_xsx_O2O2 has the advantage of physical simplicity; it measures O2-O2 absorption in terms of O2-O2 concentration. 
     ! This allows, e.g., O2-O2 absorption cross-sections to be directly compared to, say, O3 cross-sections.
     ! Disadvantage is the uncertainty of O2-O2 concentration.
     ! Shardanand (Sha77 p. 436, 527) quotes an interaction coefficient for O2 + O2 <-> O2-O2 of K = 2.6e-23 cm3 mlc-1.
     ! GOB90 are skeptical of applying any equilibrium reaction coefficients to back out O2-O2 concentration from O2 concentration.
     ! Since SPS98 and MCB98 agree with GOB90, I disavow using O2-O2 concentrations on physical principle.
     ! I continue to use faux O2-O2 concentrations in swnb because changing clm and swnb to use either of the following more physical quantities will be time-consuming.
     
     ! The first physically based method is based on O2 number concentration.
     ! odal_O2O2_PUNC_O2_PUNP_O2: measured absorption optical depth of O2-O2 per unit number concentration of O2 per unit number path of O2.
     ! This is the quantity measured by GOB90.
     ! PUNC stands for "per unit number concentration".
     ! PUNP stands for "per unit number path".
     ! Multiply odal_O2O2_PUNC_O2_PUNP_O2 by number concentration of O2 molecules (# m-3) to obtain absorption cross-section of O2-O2 per molecule O2. 
     ! Multiply odal_O2O2_PUNC_O2_PUNP_O2 by O2 number concentration (# m-3) then by O2 number path (# m-2) to obtain absorption optical depth of O2-O2.
     ! Advantage of this method is it does not rely on uncertain equilibrium rate coefficients for O2 + O2 <-> O2-O2.
     ! There are two disadvantages to this method:
     ! 1. Computation of O2-O2 absorption optical depth in RT program is slightly non-standard because it requires knowing O2 concentration and path in the O2-O2 loop.
     ! 2. odal_O2O2_PUNC_O2_PUNP_O2 is O(1.0e-46) and thus requires double precision (eight byte) representation,
     ! Remainder of RT code has (barely) managed to avoid double precision until now
     
     ! The second, and ultimately better, physically based method is based on O2 mass concentration.
     ! odal_O2O2_PUMC_O2_PUMP_O2: measured absorption optical depth of O2-O2 per unit mass concentration of O2 per unit mass path of O2.
     ! odal_O2O2_PUMC_O2_PUMP_O2 translates to "layer absorption optical depth of O2-O2 per unit mass concentration of O2 per unit mass path of O2".
     ! PUMC stands for "per unit mass concentration".
     ! PUMP stands for "per unit mass path".
     ! Multiply odal_O2O2_PUMC_O2_PUMP_O2 by mass concentration (kg m-3) of O2 (not O2-O2) to obtain absorption optical depth of O2-O2 per unit mass path of O2.
     ! Multiply odal_O2O2_PUMC_O2_PUMP_O2 by mass concentration (kg m-3) of O2 (not O2-O2) then by mass path (kg m-2) of O2 (not O2-O2) to obtain absorption optical depth of O2-O2.
     ! A main advantage of odal_O2O2_PUMC_O2_PUMP_O2 is that using it requires no assumptions about O2-O2 concentration.
     ! All one needs to know to use odal_O2O2_PUMC_O2_PUMP_O2 is the abundance of O2, which is trivial to derive.
     ! Another advantage is that odal_O2O2_PUMC_O2_PUMP_O2 is always representable in single precision (4 byte) storage.
     
     ! Compute diagnostic variables
     wvl_min(1)=wvl_ctr(1)-0.5*(wvl_ctr(2)-wvl_ctr(1))
     do bnd_idx=2,bnd_nbr
        wvl_min(bnd_idx)=0.5*(wvl_ctr(bnd_idx-1)+wvl_ctr(bnd_idx))
     enddo                  ! end loop over bnd
     do bnd_idx=1,bnd_nbr-1
        wvl_max(bnd_idx)=0.5*(wvl_ctr(bnd_idx)+wvl_ctr(bnd_idx+1))
     enddo                  ! end loop over bnd
     wvl_max(bnd_nbr)=wvl_ctr(bnd_nbr)+0.5*(wvl_ctr(bnd_nbr)-wvl_ctr(bnd_nbr-1))
     
     do bnd_idx=1,bnd_nbr
        wvl(bnd_idx)=wvl_ctr(bnd_idx)
        bnd(bnd_idx)=wvl_ctr(bnd_idx)
        double_foo=odal_O2O2_PUNC_O2_PUNP_O2(bnd_idx)/dble(k_O2_O2)
        abs_xsx_O2O2(bnd_idx)=double_foo ! defensive programming
        abs_xsx_O2O2_tpt_rfr(bnd_idx)=tpt_std ! All GOB90 data taken at 296 K
        ! Temperature dependence for GOB90 data not currently implemented
        abs_xsx_O2O2_dadT(bnd_idx)=0.0 ! GOB90 showed very little temperature dependence 
     enddo                  ! end loop over bnd
     
     ! Set O2-N2 collision-induced absorption cross sections
     do bnd_idx=1,bnd_nbr
        if (wvl_ctr(bnd_idx) >= 1.137e-6.and.wvl_ctr(bnd_idx) <= 1.350e-6) then
           odal_O2N2_PUNC_O2_PUNP_N2(bnd_idx)=N2_cll_prt_fsh*odal_O2O2_PUNC_O2_PUNP_O2(bnd_idx)
           abs_xsx_O2N2(bnd_idx)=N2_cll_prt_fsh*abs_xsx_O2O2(bnd_idx)
        else                ! endif wvl
           odal_O2N2_PUNC_O2_PUNP_N2(bnd_idx)=0.0
           abs_xsx_O2N2(bnd_idx)=0.0
        endif               ! endelse wvl
     enddo                  ! end loop over bnd
     
  endif                     ! GOB90 data
  
  close (fl_in_unit)
  write (6,'(a20,1x,a)') 'Read input data from',fl_in(1:ftn_strlen(fl_in))
  
  ! Get TOA solar spectrum
  slr_spc_xtr_typ=xtr_fll_ngh !+xtr_vrb
  call slr_spc_get(fl_slr,wvl_min,wvl_max,bnd_nbr,flx_slr_frc,slr_spc_xtr_typ,slr_spc_xtr_typ)
  
  ! Compute diagnostic variables
  do bnd_idx=1,bnd_nbr
     abs_cff_mss_O2O2(bnd_idx)=abs_xsx_O2O2(bnd_idx)*Avagadro/mmw_O2O2
     double_foo= &
          odal_O2O2_PUNC_O2_PUNP_O2(bnd_idx)*(dble(Avagadro)/dble(mmw_O2))**2
     odal_O2O2_PUMC_O2_PUMP_O2(bnd_idx)=double_foo ! Defensive programming
     double_foo= &
          odal_O2N2_PUNC_O2_PUNP_N2(bnd_idx)*dble(Avagadro)*dble(Avagadro)/(dble(mmw_O2)*dble(mmw_N2))
     odal_O2N2_PUMC_O2_PUMP_N2(bnd_idx)=double_foo ! Defensive programming
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
  
  ! Begin netCDF output routines
  rcd=nf90_wrp_create(fl_out,nf90_clobber,nc_id)
  ! Define dimension IDs
  rcd=rcd+nf90_def_dim(nc_id,'bnd',bnd_nbr,bnd_dim_id)
  if (rcd /= nf90_noerr) call nf90_err_exit(rcd,fl_out)
  rcd=rcd+nf90_def_dim(nc_id,'grd',bnd_nbr+1,grd_dim_id)
  if (rcd /= nf90_noerr) call nf90_err_exit(rcd,fl_out)
  
  ! Variable definitions
  rcd=rcd+nf90_def_var(nc_id,'abs_cff_mss_O2O2',nf90_float,bnd_dim_id,abs_cff_mss_O2O2_id)
  rcd=rcd+nf90_def_var(nc_id,'abs_xsx_O2N2',nf90_float,bnd_dim_id,abs_xsx_O2N2_id)
  rcd=rcd+nf90_def_var(nc_id,'abs_xsx_O2O2',nf90_double,bnd_dim_id,abs_xsx_O2O2_id)
  rcd=rcd+nf90_def_var(nc_id,'abs_xsx_O2O2_dadT',nf90_float,bnd_dim_id,abs_xsx_O2O2_dadT_id)
  rcd=rcd+nf90_def_var(nc_id,'abs_xsx_O2O2_tpt_rfr',nf90_float,bnd_dim_id,abs_xsx_O2O2_tpt_rfr_id)
  rcd=rcd+nf90_def_var(nc_id,'bnd',nf90_float,bnd_dim_id,bnd_id)
  rcd=rcd+nf90_def_var(nc_id,'flx_bnd_dwn_TOA',nf90_float,bnd_dim_id,flx_bnd_dwn_TOA_id)
  rcd=rcd+nf90_def_var(nc_id,'flx_bnd_pht_dwn_TOA',nf90_float,bnd_dim_id,flx_bnd_pht_dwn_TOA_id)
  rcd=rcd+nf90_def_var(nc_id,'flx_slr_frc',nf90_float,bnd_dim_id,flx_slr_frc_id)
  rcd=rcd+nf90_def_var(nc_id,'flx_spc_dwn_TOA',nf90_float,bnd_dim_id,flx_spc_dwn_TOA_id)
  rcd=rcd+nf90_def_var(nc_id,'idx_rfr_air_STP',nf90_float,bnd_dim_id,idx_rfr_air_STP_id)
  rcd=rcd+nf90_def_var(nc_id,'nrg_pht',nf90_float,bnd_dim_id,nrg_pht_id)
  rcd=rcd+nf90_def_var(nc_id,'odal_O2N2_PUMC_O2_PUMP_N2',nf90_float,bnd_dim_id,odal_O2N2_PUMC_O2_PUMP_N2_id)
  rcd=rcd+nf90_def_var(nc_id,'odal_O2O2_PUMC_O2_PUMP_O2',nf90_float,bnd_dim_id,odal_O2O2_PUMC_O2_PUMP_O2_id)
  rcd=rcd+nf90_def_var(nc_id,'odal_O2N2_PUNC_O2_PUNP_N2',nf90_double,bnd_dim_id,odal_O2N2_PUNC_O2_PUNP_N2_id)
  rcd=rcd+nf90_def_var(nc_id,'odal_O2O2_PUNC_O2_PUNP_O2',nf90_double,bnd_dim_id,odal_O2O2_PUNC_O2_PUNP_O2_id)
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
  rcd=rcd+nf90_put_att(nc_id,nf90_global,'src_wgt_sng',src_wgt_sng(1:ftn_strlen(src_wgt_sng)))
  
  ! Add english text descriptions
  rcd=rcd+nf90_put_att(nc_id,abs_cff_mss_O2O2_id,'long_name','O2-O2 mass absorption coefficient (per kg O2-O2)')
  rcd=rcd+nf90_put_att(nc_id,abs_xsx_O2O2_dadT_id,'long_name','Slope of absorption cross section temperature dependence')
  rcd=rcd+nf90_put_att(nc_id,abs_xsx_O2O2_tpt_rfr_id,'long_name','Valid temperature for absorption cross section')
  rcd=rcd+nf90_put_att(nc_id,abs_xsx_O2O2_id,'long_name', &
       'O2-O2 collision-induced absorption cross section (per O2-O2 complex)')
  rcd=rcd+nf90_put_att(nc_id,abs_xsx_O2N2_id,'long_name', &
   'O2-N2 collision-induced absorption cross section (per O2-N2 complex)')
  rcd=rcd+nf90_put_att(nc_id,odal_O2O2_PUMC_O2_PUMP_O2_id,'long_name', &
       'O2-O2 layer absorption optical depth (per kg m-3 O2 per kg m-2 O2)')
  rcd=rcd+nf90_put_att(nc_id,odal_O2N2_PUMC_O2_PUMP_N2_id,'long_name', &
       'O2-N2 layer absorption optical depth (per kg m-3 O2 per kg m-2 N2)')
  rcd=rcd+nf90_put_att(nc_id,odal_O2O2_PUNC_O2_PUNP_O2_id,'long_name', &
       'O2-O2 layer absorption optical depth (per mlc m-3 O2 per mlc m-2 O2)')
  rcd=rcd+nf90_put_att(nc_id,odal_O2N2_PUNC_O2_PUNP_N2_id,'long_name', &
       'O2-N2 layer absorption optical depth (per mlc m-3 O2 per mlc m-2 N2)')
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
  
  ! Add units
  rcd=rcd+nf90_put_att(nc_id,abs_cff_mss_O2O2_id,'units','meter2 kilogram-1')
  rcd=rcd+nf90_put_att(nc_id,abs_xsx_O2N2_id,'units','meter2 molecule-1')
  rcd=rcd+nf90_put_att(nc_id,abs_xsx_O2O2_dadT_id,'units','meter2 kelvin-1')
  rcd=rcd+nf90_put_att(nc_id,abs_xsx_O2O2_id,'units','meter2 molecule-1')
  rcd=rcd+nf90_put_att(nc_id,abs_xsx_O2O2_tpt_rfr_id,'units','kelvin')
  rcd=rcd+nf90_put_att(nc_id,bnd_id,'units','meter')
  rcd=rcd+nf90_put_att(nc_id,flx_bnd_dwn_TOA_id,'units','watt meter-2')
  rcd=rcd+nf90_put_att(nc_id,flx_bnd_pht_dwn_TOA_id,'units','photon meter-2 second-1')
  rcd=rcd+nf90_put_att(nc_id,flx_slr_frc_id,'units','fraction')
  rcd=rcd+nf90_put_att(nc_id,flx_spc_dwn_TOA_id,'units','watt meter-2 meter-1')
  rcd=rcd+nf90_put_att(nc_id,idx_rfr_air_STP_id,'units','fraction')
  rcd=rcd+nf90_put_att(nc_id,nrg_pht_id,'units','joule photon-1')
  rcd=rcd+nf90_put_att(nc_id,odal_O2N2_PUMC_O2_PUMP_N2_id,'units','meter5 kilogram-2')
  rcd=rcd+nf90_put_att(nc_id,odal_O2N2_PUNC_O2_PUNP_N2_id,'units','meter5 molecule-2')
  rcd=rcd+nf90_put_att(nc_id,odal_O2O2_PUMC_O2_PUMP_O2_id,'units','meter5 kilogram-2')
  rcd=rcd+nf90_put_att(nc_id,odal_O2O2_PUNC_O2_PUNP_O2_id,'units','meter5 molecule-2')
  rcd=rcd+nf90_put_att(nc_id,wvl_ctr_id,'units','meter')
  rcd=rcd+nf90_put_att(nc_id,wvl_grd_id,'units','meter')
  rcd=rcd+nf90_put_att(nc_id,wvl_max_id,'units','meter')
  rcd=rcd+nf90_put_att(nc_id,wvl_min_id,'units','meter')
  rcd=rcd+nf90_put_att(nc_id,wvl_dlt_id,'units','meter')
  
  ! All dimensions, variables, and attributes have been defined
  rcd=rcd+nf90_enddef(nc_id)
  
  ! Write data
  rcd=rcd+nf90_put_var(nc_id,odal_O2N2_PUNC_O2_PUNP_N2_id,odal_O2N2_PUNC_O2_PUNP_N2)
  rcd=rcd+nf90_put_var(nc_id,odal_O2O2_PUNC_O2_PUNP_O2_id,odal_O2O2_PUNC_O2_PUNP_O2)
  rcd=rcd+nf90_put_var(nc_id,abs_cff_mss_O2O2_id,abs_cff_mss_O2O2)
  rcd=rcd+nf90_put_var(nc_id,abs_xsx_O2N2_id,abs_xsx_O2N2)
  rcd=rcd+nf90_put_var(nc_id,abs_xsx_O2O2_dadT_id,abs_xsx_O2O2_dadT)
  rcd=rcd+nf90_put_var(nc_id,abs_xsx_O2O2_id,abs_xsx_O2O2)
  rcd=rcd+nf90_put_var(nc_id,abs_xsx_O2O2_tpt_rfr_id,abs_xsx_O2O2_tpt_rfr)
  rcd=rcd+nf90_put_var(nc_id,bnd_id,bnd)
  rcd=rcd+nf90_put_var(nc_id,flx_bnd_dwn_TOA_id,flx_bnd_dwn_TOA)
  rcd=rcd+nf90_put_var(nc_id,flx_bnd_pht_dwn_TOA_id,flx_bnd_pht_dwn_TOA)
  rcd=rcd+nf90_put_var(nc_id,flx_slr_frc_id,flx_slr_frc)
  rcd=rcd+nf90_put_var(nc_id,flx_spc_dwn_TOA_id,flx_spc_dwn_TOA)
  rcd=rcd+nf90_put_var(nc_id,idx_rfr_air_STP_id,idx_rfr_air_STP)
  rcd=rcd+nf90_put_var(nc_id,nrg_pht_id,nrg_pht)
  rcd=rcd+nf90_put_var(nc_id,odal_O2N2_PUMC_O2_PUMP_N2_id,odal_O2N2_PUMC_O2_PUMP_N2)
  rcd=rcd+nf90_put_var(nc_id,odal_O2O2_PUMC_O2_PUMP_O2_id,odal_O2O2_PUMC_O2_PUMP_O2)
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
  
  if (WGT_TRN) then
     ! Weight high resolution absorption cross-sections by atmospheric transmission of atmosphere with all constituents except O2O2
     call wgt_get(fl_wgt,wgt_nm,wvl_grd,bnd_nbr,wgt_spc)         
     do idx=1,bnd_nbr
        xsx_wgt_flx(idx)=(mlt_fct*abs_xsx_O2O2(idx))*wgt_spc(idx)

     enddo                  ! end loop over bnd
  else                      ! !WGT_TRN
     ! Weight high resolution absorption cross sections by high resolution TOA solar flux
     do idx=1,bnd_nbr
        xsx_wgt_flx(idx)=(mlt_fct*abs_xsx_O2O2(idx))*flx_spc_dwn_TOA(idx)
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
  do idx=1,bnd_nbr_CCM
     xsx_wgt_flx_CCM(idx)=xsx_wgt_flx_CCM(idx)/mlt_fct
     ! NB: next line not necessary since xsx_wgt_flx is not used again
     ! xsx_wgt_flx(idx)=xsx_wgt_flx(idx)/mlt_fct
  enddo                  ! end loop over bnd

  if (dbg_lvl == dbg_crr) then
     ! Examine weights
     write (6,'(5(a,1x))') '#  ','wvl        ','xsx_abs    ','flx TOA    ',wgt_nm(1:ftn_strlen(wgt_nm))
     write (6,'(5(a,1x))') '   ','m          ','m-2 mlc-1  ','W m-2 m-1  ','W m-2 m-1  '
     do idx=1,bnd_nbr
        write (6,'(i4,1x,4(es10.3,1x))') &
             idx,wvl_ctr(idx),abs_xsx_O2O2(idx),flx_spc_dwn_TOA(idx),wgt_spc(idx)
     enddo                  ! end loop over bnd
  endif                     ! endif dbg
  
  ! Normalize flux-weighted absorption cross sections
  do idx=1,bnd_nbr_CCM
     flx_spc_dwn_TOA_CCM(idx)=flx_slr_frc_CCM(idx)*slr_cst_CCM/wvl_dlt_CCM(idx)
     if (WGT_TRN) then
        if (wgt_spc_CCM(idx) /= 0.0) then
           abs_xsx_O2O2_CCM(idx)=xsx_wgt_flx_CCM(idx)/wgt_spc_CCM(idx)
        else
           write (6,'(a,a,i2,a)') prg_nm(1:ftn_strlen(prg_nm)), &
                ': WARNING wgt_spc_CCM(',idx,') = 0.0 Setting cross section equal to zero.'
           abs_xsx_O2O2_CCM(idx)=0.0
        endif               ! endif not degenerate
     else
        abs_xsx_O2O2_CCM(idx)=xsx_wgt_flx_CCM(idx)/flx_spc_dwn_TOA_CCM(idx)
     endif                  ! endif WGT_TRN
     abs_cff_mss_O2O2_CCM(idx)=abs_xsx_O2O2_CCM(idx)*Avagadro/mmw_O2O2
     odal_O2O2_PUNC_O2_PUNP_O2_CCM(idx)=dble(abs_xsx_O2O2_CCM(idx))*dble(k_O2_O2) ! NB: Double precision arithmetic is required here
     double_foo= &
             odal_O2O2_PUNC_O2_PUNP_O2_CCM(idx)*(dble(Avagadro)/dble(mmw_O2))**2
     odal_O2O2_PUMC_O2_PUMP_O2_CCM(idx)=double_foo ! Defensive programming
  enddo                     ! end loop over CCM bnd
  
  ! Do the same thing for O2-N2
  if (WGT_TRN) then
     ! Weight high resolution absorption cross-sections by atmospheric transmission of atmosphere with all constituents except H2OH2O
     call wgt_get(fl_wgt,wgt_nm,wvl_grd,bnd_nbr,wgt_spc)         
     do idx=1,bnd_nbr
        xsx_wgt_flx(idx)=abs_xsx_O2N2(idx)*wgt_spc(idx)
     enddo                  ! end loop over bnd
  else                      ! !WGT_TRN
     ! Weight high resolution absorption cross sections by high resolution TOA solar flux
     do idx=1,bnd_nbr
        xsx_wgt_flx(idx)=abs_xsx_O2N2(idx)*flx_spc_dwn_TOA(idx)
     enddo                  ! end loop over bnd
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
  do idx=1,bnd_nbr_CCM
     flx_spc_dwn_TOA_CCM(idx)=flx_slr_frc_CCM(idx)*slr_cst_CCM/wvl_dlt_CCM(idx)
     if (WGT_TRN) then
        if (wgt_spc_CCM(idx) /= 0.0) then
           abs_xsx_O2N2_CCM(idx)=xsx_wgt_flx_CCM(idx)/wgt_spc_CCM(idx)
        else
           write (6,'(a,a,i2,a)') prg_nm(1:ftn_strlen(prg_nm)), &
                ': WARNING wgt_spc_CCM(',idx,') = 0.0 Setting cross section equal to zero.'
           abs_xsx_O2N2_CCM(idx)=0.0
        endif               ! endif not degenerate
     else
        abs_xsx_O2N2_CCM(idx)=xsx_wgt_flx_CCM(idx)/flx_spc_dwn_TOA_CCM(idx)
     endif                  ! endif WGT_TRN
     abs_cff_mss_O2N2_CCM(idx)=abs_xsx_O2N2_CCM(idx)*Avagadro/mmw_O2N2
     odal_O2N2_PUNC_O2_PUNP_N2_CCM(idx)=dble(abs_xsx_O2N2_CCM(idx))*dble(k_O2_O2) ! NB: Double precision arithmetic is required here
     double_foo= &
          odal_O2N2_PUNC_O2_PUNP_N2_CCM(idx)*dble(Avagadro)*dble(Avagadro)/(dble(mmw_O2)*dble(mmw_N2))
     odal_O2N2_PUMC_O2_PUMP_N2_CCM(idx)=double_foo ! Defensive programming
  enddo                     ! end loop over CCM bnd
  
  if (dbg_lvl == dbg_old) then
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
  
  ! Add O2O2 data on CCM grid to netCDF output file
  rcd=rcd+nf90_wrp_open(fl_out,nf90_write,nc_id)
  if (rcd /= nf90_noerr) call nf90_err_exit(rcd,fl_out)
  
  ! Put output file in define mode
  rcd=rcd+nf90_redef(nc_id)
  if (rcd /= nf90_noerr) call nf90_err_exit(rcd,fl_out)
  
  ! Variable definitions
  rcd=rcd+nf90_def_var(nc_id,'abs_cff_mss_O2O2_CCM',nf90_float,bnd_dim_id,abs_cff_mss_O2O2_id)
  rcd=rcd+nf90_def_var(nc_id,'abs_xsx_O2O2_CCM',nf90_float,bnd_dim_id,abs_xsx_O2O2_id)
  rcd=rcd+nf90_def_var(nc_id,'flx_slr_frc_CCM',nf90_float,bnd_dim_id,flx_slr_frc_id)
  rcd=rcd+nf90_def_var(nc_id,'flx_spc_dwn_TOA_CCM',nf90_float,bnd_dim_id,flx_spc_dwn_TOA_id)
  rcd=rcd+nf90_def_var(nc_id,'odal_O2N2_PUMC_O2_PUMP_N2_CCM',nf90_float,bnd_dim_id,odal_O2N2_PUMC_O2_PUMP_N2_id)
  rcd=rcd+nf90_def_var(nc_id,'odal_O2O2_PUMC_O2_PUMP_O2_CCM',nf90_float,bnd_dim_id,odal_O2O2_PUMC_O2_PUMP_O2_id)
  rcd=rcd+nf90_def_var(nc_id,'odal_O2O2_PUNC_O2_PUNP_O2_CCM',nf90_double,bnd_dim_id,odal_O2O2_PUNC_O2_PUNP_O2_id)
  
  ! Add global attributes
  
  ! Add english text descriptions
  rcd=rcd+nf90_put_att(nc_id,odal_O2O2_PUMC_O2_PUMP_O2_id,'long_name', &
       'O2-O2 layer absorption optical depth (per kg m-3 O2 per kg m-2 O2)')
  rcd=rcd+nf90_put_att(nc_id,odal_O2N2_PUMC_O2_PUMP_N2_id,'long_name', &
       'O2-N2 layer absorption optical depth (per kg m-3 O2 per kg m-2 N2)')
  rcd=rcd+nf90_put_att(nc_id,odal_O2O2_PUNC_O2_PUNP_O2_id,'long_name', &
       'O2-O2 layer absorption optical depth (per mlc m-3 O2 per mlc m-2 O2)')
  rcd=rcd+nf90_put_att(nc_id,abs_xsx_O2O2_id,'long_name', &
       'O2-O2 collision-induced absorption cross section (per O2-O2 complex)')
  rcd=rcd+nf90_put_att(nc_id,abs_cff_mss_O2O2_id,'long_name','O2-O2 mass absorption coefficient')
  rcd=rcd+nf90_put_att(nc_id,flx_slr_frc_id,'long_name','Fraction of solar flux in band: ' // fl_slr)
  rcd=rcd+nf90_put_att(nc_id,flx_spc_dwn_TOA_id,'long_name','Spectral solar insolation at TOA')
  
  ! Add units
  rcd=rcd+nf90_put_att(nc_id,abs_cff_mss_O2O2_id,'units','meter2 kilogram-1')
  rcd=rcd+nf90_put_att(nc_id,abs_xsx_O2O2_id,'units','meter2')
  rcd=rcd+nf90_put_att(nc_id,flx_slr_frc_id,'units','fraction')
  rcd=rcd+nf90_put_att(nc_id,flx_spc_dwn_TOA_id,'units','watt meter-2 meter-1')
  rcd=rcd+nf90_put_att(nc_id,odal_O2N2_PUMC_O2_PUMP_N2_id,'units','meter5 kilogram-2')
  rcd=rcd+nf90_put_att(nc_id,odal_O2O2_PUMC_O2_PUMP_O2_id,'units','meter5 kilogram-2')
  rcd=rcd+nf90_put_att(nc_id,odal_O2O2_PUNC_O2_PUNP_O2_id,'units','meter5 molecule-2')
  
  ! All dimensions, variables, and attributes have been defined
  rcd=rcd+nf90_enddef(nc_id)
  if (rcd /= nf90_noerr) call nf90_err_exit(rcd,fl_out)
  
  ! Write out data
  rcd=rcd+nf90_put_var(nc_id,odal_O2O2_PUNC_O2_PUNP_O2_id,odal_O2O2_PUNC_O2_PUNP_O2_CCM)
  rcd=rcd+nf90_put_var(nc_id,abs_cff_mss_O2O2_id,abs_cff_mss_O2O2_CCM)
  rcd=rcd+nf90_put_var(nc_id,abs_xsx_O2O2_id,abs_xsx_O2O2_CCM)
  rcd=rcd+nf90_put_var(nc_id,flx_slr_frc_id,flx_slr_frc_CCM)
  rcd=rcd+nf90_put_var(nc_id,flx_spc_dwn_TOA_id,flx_spc_dwn_TOA_CCM)
  rcd=rcd+nf90_put_var(nc_id,odal_O2N2_PUMC_O2_PUMP_N2_id,odal_O2N2_PUMC_O2_PUMP_N2_CCM)
  rcd=rcd+nf90_put_var(nc_id,odal_O2O2_PUMC_O2_PUMP_O2_id,odal_O2O2_PUMC_O2_PUMP_O2_CCM)
  
  ! Close output file
  rcd=rcd+nf90_close(nc_id)
  write (6,'(a28,1x,a)') 'Wrote O2O2 data on CCM grid to',fl_out(1:ftn_strlen(fl_out))
  if (rcd /= nf90_noerr) write (6,'(a,a,i4,a)') prg_nm(1:ftn_strlen(prg_nm)),': ERROR rcd = ',rcd,' on exit'
  
  ! Convert absorption coefficients to CGS and output block for radcsw
  do idx=1,bnd_nbr_CCM
     odal_O2O2_PUMC_O2_PUMP_O2_CCMCG(idx)=odal_O2O2_PUMC_O2_PUMP_O2_CCM(idx)*(100.0**5)/(1000.0**2) ! m5 kg-2 = m3 kg-1 O2 m2 kg-1 O2 --> cm5 g-2 = cm3 g-1 O2 cm2 g-1 O2 
     odal_O2N2_PUMC_O2_PUMP_N2_CCMCG(idx)=odal_O2N2_PUMC_O2_PUMP_N2_CCM(idx)*(100.0**5)/(1000.0**2) ! m5 kg-2 = m3 kg-1 O2 m2 kg-1 O2 --> cm5 g-2 = cm3 g-1 O2 cm2 g-1 O2 
  enddo                     ! end loop over CCM bnd
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
  if (allocated(xsx_wgt_flx)) deallocate(xsx_wgt_flx,stat=rcd)
  
1000 continue
  
  call exit(exit_status)
end program O2X

