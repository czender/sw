! $Id$

program clm
  
  ! Purpose: Turn atmospheric profile input files (CRM text or netCDF CLM),   
  ! or NWS radiosonde input files, or AFGL standard atmosphere files, 
  ! into netCDF CLM output files. 

  ! Copyright (C) 1994--2016 Charlie Zender
  ! License: GNU General Public License (GPL) Version 3
  ! See http://www.gnu.org/copyleft/gpl.html for full license text
  ! The original author of this software, Charlie Zender, seeks to improve
  ! it with your suggestions, contributions, bug-reports, and patches.
  ! Charlie Zender <zender at uci dot edu>
  ! Department of Earth System Science
  ! University of California, Irvine
  ! Irvine, CA 92697-3100
  
  ! Compilation:
  ! cd ${HOME}/aca; make -W clm.F90 OPTS=D NETCDF4=Y clm; cd -
  ! cd ${HOME}/aca; make -W clm.F90 OPTS=D clm; cd -
  ! cd ${HOME}/aca; make -W clm.F90 clm; cd -
  ! cd ${HOME}/aca; make clm; cd -
  ! cd ${HOME}/aca; make OPTS=D clm; cd -
  ! scp ~/aca/clm.F90 givre.ess.uci.edu:aca
  
  ! Output format is called CoLuMn (CLM) file because it describes properties of atmospheric column 
  ! CLM files evolved from CCM Column Radiation Model (CRM) text input files
  ! The correct format of a CRM input text file is whatever works with the CCM CRM
  ! However, CRM files are not flexible enough to include information necessary for more
  ! complex, state of the art radiative transfer and thermodynamic modeling.
  ! Thus CLM files are designed to be backward compatible with CRM files, but also provide 
  ! the flexibility and extensibility required for more detailed work.
  ! netCDF input files must have (at least) lev, p, T, and q_H2O stored in SI units
  
  ! NB: Reading from and writing to the same netCDF file is allowed because fl_in is
  ! closed before fl_out is opened. This may sound dangerous but serves a useful purpose:
  ! If netCDF in.nc was created by time interpolating between in_1.nc and in_2.nc,
  ! then none of the non-linear fields in in.nc, e.g., svp_H2O_lqd, is correct.
  ! clm is designed to fix the non-linear derived fields of in.nc 
  ! This is accomplished by only reading in primary fields from in.nc and recomputing
  ! all derived fields.
  
  ! NB: Setting trace gas concentrations to exactly zero is dangerous when the gas
  ! mixing ratio is used in the narrow band model. A divide by zero error may occur
  ! when you normalize by the pressure weighted mass path (which is zero).
  
  ! NB: Write non-fatal errors to stderr so that output CRM text files, if any, are not corrupted 
  
  ! Usage:
  ! clm -i ${DATA}/aca/mls_cld.txt -o ${DATA}/aca/mls_cld.nc
  ! clm -i ${DATA}/aca/mls_clr.txt -o ${DATA}/aca/mls_clr.nc
  ! clm -i ${DATA}/aca/mls_clr_2xCO2.txt -o ${DATA}/aca/mls_clr_2xCO2.nc
  ! clm --drc_in ${DATA}/aca --drc_out ${DATA}/aca
  ! clm --drc_in ${DATA}/aca -i mls_cld.txt --drc_out ${DATA}/aca -o mls_cld.nc
  ! clm --drc_in ${DATA}/aca -i mls_cld.txt --drc_out ${DATA}/aca -o mls_cld.nc -t clm.txt
  ! clm -l 92 -i ${DATA}/aca/mls_icrccm_92lvl.txt -o ${DATA}/aca/mls_icrccm_92lvl.nc
  ! clm -l 18 -i ${DATA}/aca/trp_icrccm_18lvl.txt -o ${DATA}/aca/trp_icrccm_18lvl.nc
  ! clm -l 35 -i ${DATA}/aca/trp_icrccm_35lvl.txt -o ${DATA}/aca/trp_icrccm_35lvl.nc
  ! clm -l 92 -i ${DATA}/aca/trp_icrccm_92lvl.txt -o ${DATA}/aca/trp_icrccm_92lvl.nc
  ! Ingest wind speeds:
  ! clm -U -l 20 -D 5 --drc_in ${HOME}/aca -i clm_wnd.txt --drc_out ${HOME}/aca -o mls_cld.nc
  ! Ingest snowpack structure:
  ! clm --lev_snw=5 -D 1 --drc_in ${DATA}/aca -i mls_snw.txt --drc_out ${DATA}/aca -o mls_snw.nc
  ! ncks -u -C -H -v '._snw' ${DATA}/aca/mls_snw.nc

  ! Print out CLM txt version of existing netCDF CLM profile:
  ! clm -n -i ${DATA}/aca/mls_clr.nc -o foo.nc
  
  ! Create saturated column for dimer studies:
  ! clm -s -D 5 -i ${HOME}/aca/dmr_clm.txt -o ${DATA}/tmp/dmr.nc
  
  ! Converting AFGL profiles (from Stefan Kinne) to CLM files is handled by
  ! ${HOME}/afgl/afgl.sh
  
  ! ARESE datasets: 
  ! "-J" (mnemonic: ) distributes command line aerosol column extinction optical depth among levels with aerosol
  ! "-A" (mnemonic: Aerosol) turns on aerosol column for text I/O, not neccessary for netCDF I/O
  ! "-Z" (mnemonic: Zenith) imposes observed SIROS CART zenith angle parameterization on surface albedo
  ! clm -Z -A -J 0.12 -n -i ${DATA}/arese/clm/951011_1200_arese_clm.nc -t ${DATA}/aca/arese_tst.txt -o ${DATA}/aca/arese_tst.nc
  ! clm -Z -A -J 0.12 -l 110 -i ${DATA}/aca/arese_tst.txt -o ${DATA}/aca/arese_tst.nc
  ! clm -n -S -Z -J 0.12 -i ${DATA}/aca/arese_tst.nc -o ${DATA}/aca/arese_tst.nc
  ! ncks -a -u -C -H -v wvl_obs_aer,wvl_obs_bga,odxc_obs_bga,odxc_obs_aer,mpc_aer,mpc_bga,dns_aer,dns_bga,ext_cff_mss_aer,ext_cff_mss_bga ${DATA}/aca/arese_tst.nc
  
  ! clm -n -W 0.5e-6 -J 0.1 -w 0.5e-6 -j 0.1 -a ${DATA}/aca/aer_mineral_dust.nc -b ${DATA}/aca/aer_mineral_dust.nc -i ${DATA}/aca/arese_tst.nc -o ${DATA}/aca/arese_tst.nc
  
  ! j_NO2:
  ! clm -n -J 1.006 -j 0.006 -a ${DATA}/aca/aer_mineral_dust.nc -b ${DATA}/aca/aer_h2so4_215K.nc -i ${DATA}/arese/clm/951011_1200_arese_clm.nc -o ${DATA}/swnb2/arese_951011_1200_clm_clr_aer_dst.nc
  ! clm -n -J 1.006 -j 0.006 -a ${DATA}/aca/aer_h2so4_300K.nc -b ${DATA}/aca/aer_h2so4_215K.nc -i ${DATA}/arese/clm/951011_1200_arese_clm.nc -o ${DATA}/swnb2/arese_951011_1200_clm_clr_aer_slf.nc
  
  ! O2-O2:
  ! ncks -a -u -F -C -H -v dns_O2_mpl_O2_clm,dns_O2_mpl_N2_clm ${DATA}/aca/mls_clr.nc
  ! ncks -a -u -F -C -H -v mpc_H2OH2O ${DATA}/aca/trp_icrccm_18lvl.nc
  
  ! Print all information about a gas:
  ! ncks -a -u -F -C -H -d lev,100000.0 -v cnc_N2O,dns_N2O,mpc_N2O,mpl_N2O,npc_N2O,npl_N2O,ppr_N2O,q_N2O,r_N2O,vmr_N2O ${DATA}/aca/mls_clr.nc
  
  use clm_mdl,only:aer_odxc_get,aer_info_get,q_o3_ntp,q_no2_ntp,q_oh_ntp,slr_crd_Bri92 ! [mdl] Column (CLM) processing
  use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
  use dmr_mdl ! [mdl] Dimers, collision complexes
  use netcdf ! [mdl] netCDF interface
  use nf90_utl ! [mdl] netCDF utilities
  use phys_cst_mdl ! [mdl] Fundamental and derived physical constants
  use slr_crd_mdl,only:slr_crd_Mic88 ! [mdl] Solar coordinate geometry
  use sng_mdl ! [mdl] String manipulation
  use tdy_mdl,only:svp_H2O_lqd_PrK78,svp_H2O_ice_PrK78 ! [mdl] Thermodynamics
  use utl_mdl,only:date_time_get ! [mdl] Utility functions (date_time_get,mnt_chk...)
  use cape_mdl ! [mdl] Parcel Profile and Column Convective Potential
  use wnd_pdf_mdl ! [mdl] Wind Speed PDF module
  use abl_trb_mdl ! [mdl] Atmospheric turbulence module
  
  implicit none
  ! Parameters
  character(len=*),parameter::CVS_Date='$Date$' ! [sng] Date string
  character(len=*),parameter::CVS_Header='$Id$' ! [sng] Full CVS Header
  character(len=*),parameter::CVS_Id='$Id$' ! [sng] CVS Identification
  character(len=*),parameter::CVS_Name='$HeadURL$' ! [sng] File name string
  character(len=*),parameter::CVS_Revision='$Revision$' ! [sng] File revision string
  character(len=*),parameter::nlc=char(0) ! [sng] NUL character = ASCII 0 = char(0)
  character(len=*),parameter::sbr_nm='clm' ! [sng] Subroutine name
  integer,parameter::fl_in_unit=73
  integer,parameter::fl_txt_unit=74
  ! Arrays are allocatable, but die if size exceeds corresponding *_nbr_max
  integer,parameter::lev_nbr_max=150
  ! Derived parameters
  integer,parameter::levp_nbr_max=lev_nbr_max+1
  ! Input
  ! Input/Output
  ! Output
  ! Local workspace
  character(80) arg_val      ! [sng] Command line argument value
  character(400) cmd_ln      ! [sng] Command line
  character(2) dsh_key       ! [sng] Command line dash and switch
  character(80) drc_in       ! [sng] Input directory
  character(80) drc_out      ! [sng] Output directory
  character(80) fl_aer
  character(80) fl_bga
  character(80) fl_in
  character(80) fl_out
  character(80) fl_txt       ! [fl] CLM text file
  character(80) opt_sng      ! [sng] Option string
  character(80) lbl
  character(26)::lcl_date_time ! Time formatted as Day Mth DD HH:MM:SS TZ YYYY
  character(80) prg_ID
  
  integer arg_idx           ! [idx] Counting index
  integer arg_nbr           ! [nbr] Number of command line arguments
  integer exit_status       ! [enm] Program exit status
  integer int_foo           ! [nbr] Integer
  integer int_foo2
  integer opt_lng           ! [nbr] Length of option
  integer rcd               ! [rcd] Return success code

  ! Wind Speed PDF bin count
  integer::PDF_bin_nbr      ! [nbr] dimension size
  integer::PDF_bin_dmn_id   ! dimension ID for PDF
  
  ! Model grid sizes
  integer::lev_nbr=18       ! [nbr] dimension size
  integer::lev_snw_nbr=0    ! [nbr] dimension size
  integer levp_nbr          ! [nbr] dimension size
  integer levp_snw_nbr      ! [nbr] dimension size
  
  integer idx               ! [idx] Counting index
  integer lev_snw_idx       ! [idx] Counting index for lev_snw
  integer levp_snw_idx      ! [idx] Counting index for levp_snw
  integer levp_idx          ! counting index
  integer lev_dmn_id        ! dimension ID for lev
  integer lev_snw_dmn_id    ! dimension ID for lev_snw
  integer levp_dmn_id       ! dimension ID for levp
  integer levp_snw_dmn_id   ! dimension ID for levp_snw
  integer nc_id             ! file handle
  ! netCDF4 
  integer::dfl_lvl=0 ! [enm] Deflate level
  integer::flg_shf=1 ! [flg] Turn on netCDF4 shuffle filter
  integer::flg_dfl=1 ! [flg] Turn on netCDF4 deflate filter
  integer::fl_out_fmt=nco_format_undefined ! [enm] Output file format
  integer::nf90_create_mode=nf90_clobber ! [enm] Mode flag for nf90_create() call
  
  ! Warning: changing prf_sng length may cause swnb problems down the line...
  character(80) prf_sng
  character(80) prf_snw_sng
  
  ! Allocatable variables
  real(selected_real_kind(p=12)),dimension(:),allocatable::cnc_O2_cnc_N2 ! [mlc2 m-6] O2 number concentration times N2 number concentration
  real(selected_real_kind(p=12)),dimension(:),allocatable::cnc_O2_cnc_O2 ! [mlc2 m-6] O2 number concentration squared
  real(selected_real_kind(p=12)),dimension(:),allocatable::cnc_O2_npl_N2 ! [mlc2 m-5] O2 number concentration times N2 number path
  real(selected_real_kind(p=12)),dimension(:),allocatable::cnc_O2_npl_O2 ! [mlc2 m-5] O2 number concentration times O2 number path
  real(selected_real_kind(p=12)),dimension(:),allocatable::cnc_O2_npl_O2_clm_frc ! [frc] Fraction of column O2-O2 at or above each level
  logical,dimension(:),allocatable::bga_lvl_flg
  real,dimension(:),allocatable::PDF_bin_spd
  real,dimension(:),allocatable::PDF_bin_wgt
  real,dimension(:),allocatable::RH
  real,dimension(:),allocatable::RH_ice
  real,dimension(:),allocatable::RH_lqd
  real,dimension(:),allocatable::alt
  real,dimension(:),allocatable::alt_dlt
  real,dimension(:),allocatable::alt_ntf
  real,dimension(:),allocatable::cld_frc
  real,dimension(:),allocatable::cnc_CFC11
  real,dimension(:),allocatable::cnc_CFC12
  real,dimension(:),allocatable::cnc_CH4
  real,dimension(:),allocatable::cnc_CO2
  real,dimension(:),allocatable::cnc_H2O
  real,dimension(:),allocatable::cnc_H2OH2O
  real,dimension(:),allocatable::cnc_N2
  real,dimension(:),allocatable::cnc_N2O
  real,dimension(:),allocatable::cnc_NO2
  real,dimension(:),allocatable::cnc_O2
  real,dimension(:),allocatable::cnc_O2O2
  real,dimension(:),allocatable::cnc_O3
  real,dimension(:),allocatable::cnc_OH
  real,dimension(:),allocatable::cnc_dry_air
  real,dimension(:),allocatable::cnc_mst_air
  real,dimension(:),allocatable::dns_CFC11
  real,dimension(:),allocatable::dns_CFC12
  real,dimension(:),allocatable::dns_CH4
  real,dimension(:),allocatable::dns_CO2
  real,dimension(:),allocatable::dns_H2O
  real,dimension(:),allocatable::dns_H2OH2O
  real,dimension(:),allocatable::dns_N2
  real,dimension(:),allocatable::dns_N2O
  real,dimension(:),allocatable::dns_NO2
  real,dimension(:),allocatable::dns_O2
  real,dimension(:),allocatable::dns_O2O2
  real,dimension(:),allocatable::dns_O2_dns_N2 ! [kg2 m-6] O2 mass concentration times N2 mass concentration
  real,dimension(:),allocatable::dns_O2_dns_O2 ! [kg2 m-6] O2 mass concentration squared
  real,dimension(:),allocatable::dns_O2_mpl_N2 ! [kg2 m-5] O2 mass concentration times N2 mass path
  real,dimension(:),allocatable::dns_O2_mpl_O2  ! [kg2 m-5] O2 mass concentration times O2 mass path
  real,dimension(:),allocatable::dns_O3
  real,dimension(:),allocatable::dns_OH
  real,dimension(:),allocatable::dns_dry_air
  real,dimension(:),allocatable::dns_mst_air
  real,dimension(:),allocatable::dns_snw ! [kg m-3] Snow density
  real,dimension(:),allocatable::dpt_dlt_snw ! [m] Snow layer thickness
  real,dimension(:),allocatable::dpt_ntf_snw ! [m] Snow depth interfaces
  real,dimension(:),allocatable::dpt_snw ! [m] Snow depth
  real,dimension(:),allocatable::foo_snw ! [K] Snow foo
  real,dimension(:),allocatable::frc_ice
  real,dimension(:),allocatable::gas_cst_mst_air
  real,dimension(:),allocatable::grv
  real,dimension(:),allocatable::lev     ! coordinate variable
  real,dimension(:),allocatable::lev_snw ! [m] Snow depth
  real,dimension(:),allocatable::levp    ! coordinate variable
  real,dimension(:),allocatable::levp_snw ! [m] Snow depth interfaces
  real,dimension(:),allocatable::mmr_mpr_snw ! [kg mpr/kg snow] Mass mixing ratio of impurities in snow
  real,dimension(:),allocatable::mmw_mst_air
  real,dimension(:),allocatable::mpl_CFC11
  real,dimension(:),allocatable::mpl_CFC12
  real,dimension(:),allocatable::mpl_CH4
  real,dimension(:),allocatable::mpl_CO2
  real,dimension(:),allocatable::mpl_CWP
  real,dimension(:),allocatable::mpl_H2O
  real,dimension(:),allocatable::mpl_H2OH2O
  real,dimension(:),allocatable::mpl_IWP
  real,dimension(:),allocatable::mpl_LWP
  real,dimension(:),allocatable::mpl_N2
  real,dimension(:),allocatable::mpl_N2O
  real,dimension(:),allocatable::mpl_NO2
  real,dimension(:),allocatable::mpl_O2
  real,dimension(:),allocatable::mpl_O2O2
  real,dimension(:),allocatable::mpl_O3
  real,dimension(:),allocatable::mpl_OH
  real,dimension(:),allocatable::mpl_aer
  real,dimension(:),allocatable::mpl_bga
  real,dimension(:),allocatable::mpl_dry_air
  real,dimension(:),allocatable::mpl_mst_air
  real,dimension(:),allocatable::npl_CFC11
  real,dimension(:),allocatable::npl_CFC12
  real,dimension(:),allocatable::npl_CH4
  real,dimension(:),allocatable::npl_CO2
  real,dimension(:),allocatable::npl_H2O
  real,dimension(:),allocatable::npl_H2OH2O
  real,dimension(:),allocatable::npl_N2
  real,dimension(:),allocatable::npl_N2O
  real,dimension(:),allocatable::npl_NO2
  real,dimension(:),allocatable::npl_O2
  real,dimension(:),allocatable::npl_O2O2
  real,dimension(:),allocatable::npl_O3
  real,dimension(:),allocatable::npl_OH
  real,dimension(:),allocatable::npl_dry_air
  real,dimension(:),allocatable::npl_mst_air
  real,dimension(:),allocatable::nrg_dry ! [J kg-1] Dry static energy
  real,dimension(:),allocatable::nrg_mst ! [J kg-1] Moist static energy
  real,dimension(:),allocatable::odxl_obs_aer
  real,dimension(:),allocatable::odxl_obs_bga
  real,dimension(:),allocatable::oneD_foo
  real,dimension(:),allocatable::pcl_DLR ![K m-1] Dry Adiabatic Lapse Rate
  real,dimension(:),allocatable::pcl_MLR ![K m-1] Moist Adiabatic Lapse Rate
  real,dimension(:),allocatable::pcl_PLR ![K m-1] Pseudo-Adiabatic Lapse Rate
  real,dimension(:),allocatable::ppr_CFC11
  real,dimension(:),allocatable::ppr_CFC12
  real,dimension(:),allocatable::ppr_CH4
  real,dimension(:),allocatable::ppr_CO2
  real,dimension(:),allocatable::ppr_H2O
  real,dimension(:),allocatable::ppr_H2OH2O
  real,dimension(:),allocatable::ppr_N2
  real,dimension(:),allocatable::ppr_N2O
  real,dimension(:),allocatable::ppr_NO2
  real,dimension(:),allocatable::ppr_O2
  real,dimension(:),allocatable::ppr_O2O2
  real,dimension(:),allocatable::ppr_O3
  real,dimension(:),allocatable::ppr_OH
  real,dimension(:),allocatable::ppr_dry_air
  real,dimension(:),allocatable::prs
  real,dimension(:),allocatable::prs_dlt
  real,dimension(:),allocatable::prs_ntf
  real,dimension(:),allocatable::prs_ttl
  real,dimension(:),allocatable::q_CFC11
  real,dimension(:),allocatable::q_CFC12
  real,dimension(:),allocatable::q_CH4
  real,dimension(:),allocatable::q_CO2
  real,dimension(:),allocatable::q_H2O
  real,dimension(:),allocatable::q_H2OH2O
  real,dimension(:),allocatable::q_H2OH2O_rcp_q_H2O
  real,dimension(:),allocatable::q_N2
  real,dimension(:),allocatable::q_N2O
  real,dimension(:),allocatable::q_NO2
  real,dimension(:),allocatable::q_O2
  real,dimension(:),allocatable::q_O2O2
  real,dimension(:),allocatable::q_O3
  real,dimension(:),allocatable::q_OH
  real,dimension(:),allocatable::qst_H2O_ice
  real,dimension(:),allocatable::qst_H2O_lqd
  real,dimension(:),allocatable::r_CFC11
  real,dimension(:),allocatable::r_CFC12
  real,dimension(:),allocatable::r_CH4
  real,dimension(:),allocatable::r_CO2
  real,dimension(:),allocatable::r_H2O
  real,dimension(:),allocatable::r_H2OH2O
  real,dimension(:),allocatable::r_N2
  real,dimension(:),allocatable::r_N2O
  real,dimension(:),allocatable::r_NO2
  real,dimension(:),allocatable::r_O2
  real,dimension(:),allocatable::r_O2O2
  real,dimension(:),allocatable::r_O3
  real,dimension(:),allocatable::r_OH
  real,dimension(:),allocatable::rds_fct_ice
  real,dimension(:),allocatable::rds_fct_lqd
  real,dimension(:),allocatable::rds_ffc_snw ! [m] Snow effective radius
  real,dimension(:),allocatable::scl_hgt
  real,dimension(:),allocatable::spc_heat_mst_air
  real,dimension(:),allocatable::tpt
  real,dimension(:),allocatable::tpt_cls
  real,dimension(:),allocatable::tpt_cls_ntf
  real,dimension(:),allocatable::tpt_dwp
  real,dimension(:),allocatable::tpt_ntf
  real,dimension(:),allocatable::tpt_ntf_snw ! [K] Snow temperature interfaces
  real,dimension(:),allocatable::tpt_pcl
  real,dimension(:),allocatable::tpt_ptn ! [K] Potential temperature
  real,dimension(:),allocatable::tpt_ptn_vrt ! [K] Virtual potential temperature
  real,dimension(:),allocatable::tpt_snw ! [K] Snow temperature
  real,dimension(:),allocatable::tpt_vrt
  real,dimension(:),allocatable::vmr_CFC11
  real,dimension(:),allocatable::vmr_CFC12
  real,dimension(:),allocatable::vmr_CH4
  real,dimension(:),allocatable::vmr_CO2
  real,dimension(:),allocatable::vmr_H2O
  real,dimension(:),allocatable::vmr_H2OH2O
  real,dimension(:),allocatable::vmr_N2
  real,dimension(:),allocatable::vmr_N2O
  real,dimension(:),allocatable::vmr_NO2
  real,dimension(:),allocatable::vmr_O2
  real,dimension(:),allocatable::vmr_O2O2
  real,dimension(:),allocatable::vmr_O3
  real,dimension(:),allocatable::vmr_OH
  ! sbc Optional wind variables
  real,dimension(:),allocatable::wnd_dir ! [deg] Wind direction
  real,dimension(:),allocatable::wnd_spd ! [m s-1] Wind speed
  real,dimension(:),allocatable::wnd_shr ! [s-1] Wind shear
  real,dimension(:),allocatable::ric_nbr ! [s-1] Wind shear
  real,dimension(:),allocatable::tke     ! [m2 s-2] TKE
  real :: q_H2O_sfc       ! [kg/kg] Surface mixing ratio
  real :: khfs            ! [mK/s]  Surface kinematic heat flux
  real :: kbfs            ! [m2/s3] Surface kinematic buoyancy flux
  real :: kqfs            ! [m s-1] Kinematic moisture flux
  real :: shflx           ! [W/m2] Sensible heat flux
  real :: lhflx           ! [W/m2] Latent heat flux
  real :: ustar           ! [m/s] Friction velocity
  real :: obk_len         ! [m] Obukhov length
  real :: wstar           ! [m/s] Convective velocity scale
  real :: abl_hgt         ! [m] Boundary layer height
  real :: wspd10m         ! [m/s] 10-m surface wind speed
  real :: tau             ! [N/m2] Surface stress
  real :: kappa
  ! !sbc
  ! Local variables for AFGL_TXT_INPUT only
  real,dimension(:),allocatable::CH4_vmr
  real,dimension(:),allocatable::CH4_vmr_ntf
  real,dimension(:),allocatable::CO_vmr
  real,dimension(:),allocatable::CO_vmr_ntf
  real,dimension(:),allocatable::H2O_vmr
  real,dimension(:),allocatable::H2O_vmr_ntf
  real,dimension(:),allocatable::N2O_vmr
  real,dimension(:),allocatable::N2O_vmr_ntf
  real,dimension(:),allocatable::O3_vmr
  real,dimension(:),allocatable::O3_vmr_ntf
  real,dimension(:),allocatable::cnc_air
  real,dimension(:),allocatable::cnc_air_ntf
  
  ! Local
  logical AFGL_TXT_INPUT
  logical CLM_NC_INPUT
  logical CLM_TXT_INPUT
  logical CRM2              ! Does CRM text file have CRM 2.x fields?
  logical cmd_ln_mpc_CWP
  logical cmd_ln_odxc_obs_aer
  logical cmd_ln_odxc_obs_bga
  logical flg_aer
  logical flg_MFRSR_obs
  logical flg_bga
  logical flg_snw ! [flg] Incorporate (read-in/write-out) multi-layer snow model
  logical flg_wnd
  logical force_cld_lvl_2b_sat
  logical force_sat_lvl_2b_cld
  logical force_cld_lvl_2b_overcast
  logical zen_ngl_dpn_alb
  integer CAPE_id
  integer CFC11_vmr_clm_id
  integer CFC12_vmr_clm_id
  integer CH4_vmr_clm_id
  integer CINE_id
  integer CO2_vmr_clm_id
  integer N2O_vmr_clm_id
  integer PDF_bin_spd_id
  integer PDF_bin_wgt_id
  integer RH_ice_id
  integer RH_id
  integer RH_lqd_id
  integer alb_sfc_NIR_dff_id
  integer alb_sfc_NIR_drc_id
  integer alb_sfc_id
  integer alb_sfc_vsb_dff_id
  integer alb_sfc_vsb_drc_id
  integer alt_cld_btm_id
  integer alt_cld_mid_id
  integer alt_cld_thick_id
  integer alt_cld_top_id
  integer alt_dlt_id
  integer alt_id
  integer alt_ntf_id
  integer cld_frc_id
  integer cnc_CFC11_id
  integer cnc_CFC12_id
  integer cnc_CH4_id
  integer cnc_CO2_id
  integer cnc_H2OH2O_id
  integer cnc_H2O_id
  integer cnc_N2O_id
  integer cnc_N2_id
  integer cnc_NO2_id
  integer cnc_O2O2_id
  integer cnc_O2_cnc_N2_id
  integer cnc_O2_cnc_O2_id
  integer cnc_O2_id
  integer cnc_O2_npl_N2_clm_id
  integer cnc_O2_npl_N2_id
  integer cnc_O2_npl_O2_clm_frc_id
  integer cnc_O2_npl_O2_clm_id
  integer cnc_O2_npl_O2_id
  integer cnc_O3_id
  integer cnc_OH_id
  integer cnc_dry_air_id
  integer cnc_mst_air_id
  integer dns_CFC11_id
  integer dns_CFC12_id
  integer dns_CH4_id
  integer dns_CO2_id
  integer dns_H2OH2O_id
  integer dns_H2O_id
  integer dns_N2O_id
  integer dns_N2_id
  integer dns_NO2_id
  integer dns_O2O2_id
  integer dns_O2_dns_N2_id
  integer dns_O2_dns_O2_id
  integer dns_O2_id
  integer dns_O2_mpl_N2_clm_id
  integer dns_O2_mpl_N2_id
  integer dns_O2_mpl_O2_clm_id
  integer dns_O2_mpl_O2_id
  integer dns_O3_id
  integer dns_OH_id
  integer dns_aer_id
  integer dns_bga_id
  integer dns_dry_air_id
  integer dns_mst_air_id
  integer dns_snw_id
  integer dpt_dlt_snw_id
  integer dpt_ntf_snw_id
  integer dpt_snw_id
  integer ext_cff_mss_aer_id
  integer ext_cff_mss_bga_id
  integer foo_snw_id
  integer frc_ice_id
  integer frc_ice_ttl_id
  integer frc_str_zen_ngl_sfc_id
  integer gas_cst_mst_air_id
  integer grv_id
  integer lat_cos_id
  integer lat_dgr_id
  integer lat_id
  integer lcl_time_hr_id
  integer lcl_yr_day_id
  integer lev_LCL_id
  integer lev_LFC_id
  integer lev_LNB_id
  integer lev_id            ! coordinate ID
  integer lev_snw_id        ! coordinate ID
  integer levp_id           ! coordinate ID
  integer levp_snw_id       ! coordinate ID
  integer mmr_mpr_snw_id
  integer mmw_mst_air_id
  integer mpc_CFC11_id
  integer mpc_CFC12_id
  integer mpc_CH4_id
  integer mpc_CO2_id
  integer mpc_CWP_id
  integer mpc_H2OH2O_id
  integer mpc_H2O_id
  integer mpc_IWP_id
  integer mpc_LWP_id
  integer mpc_N2O_id
  integer mpc_N2_id
  integer mpc_NO2_id
  integer mpc_O2O2_id
  integer mpc_O2_id
  integer mpc_O3_DU_id
  integer mpc_O3_id
  integer mpc_OH_id
  integer mpc_aer_id
  integer mpc_bga_id
  integer mpc_dry_air_id
  integer mpc_mst_air_id
  integer mpl_CFC11_id
  integer mpl_CFC12_id
  integer mpl_CH4_id
  integer mpl_CO2_id
  integer mpl_CWP_id
  integer mpl_H2OH2O_id
  integer mpl_H2O_id
  integer mpl_IWP_id
  integer mpl_LWP_id
  integer mpl_N2O_id
  integer mpl_N2_id
  integer mpl_NO2_id
  integer mpl_O2O2_id
  integer mpl_O2_id
  integer mpl_O3_id
  integer mpl_OH_id
  integer mpl_aer_id
  integer mpl_bga_id
  integer mpl_dry_air_id
  integer mpl_mst_air_id
  integer npc_CFC11_id
  integer npc_CFC12_id
  integer npc_CH4_id
  integer npc_CO2_id
  integer npc_H2OH2O_id
  integer npc_H2O_id
  integer npc_N2O_id
  integer npc_N2_id
  integer npc_NO2_id
  integer npc_O2O2_id
  integer npc_O2_id
  integer npc_O3_id
  integer npc_OH_id
  integer npc_dry_air_id
  integer npc_mst_air_id
  integer npl_CFC11_id
  integer npl_CFC12_id
  integer npl_CH4_id
  integer npl_CO2_id
  integer npl_H2OH2O_id
  integer npl_H2O_id
  integer npl_N2O_id
  integer npl_N2_id
  integer npl_NO2_id
  integer npl_O2O2_id
  integer npl_O2_id
  integer npl_O3_id
  integer npl_OH_id
  integer npl_dry_air_id
  integer npl_mst_air_id
  integer nrg_dry_id
  integer nrg_mst_id
  integer odxc_obs_aer_id
  integer odxc_obs_bga_id
  integer odxl_obs_aer_id
  integer odxl_obs_bga_id
  integer oneD_foo_id
  integer oro_id
  integer pcl_DLR_id
  integer pcl_MLR_id
  integer pcl_PLR_id
  integer ppr_CFC11_id
  integer ppr_CFC12_id
  integer ppr_CH4_id
  integer ppr_CO2_id
  integer ppr_H2OH2O_id
  integer ppr_H2O_id
  integer ppr_N2O_id
  integer ppr_N2_id
  integer ppr_NO2_id
  integer ppr_O2O2_id
  integer ppr_O2_id
  integer ppr_O3_id
  integer ppr_OH_id
  integer ppr_dry_air_id
  integer prs_cld_btm_id
  integer prs_cld_mid_id
  integer prs_cld_thick_id
  integer prs_cld_top_id
  integer prs_dlt_id
  integer prs_id
  integer prs_ntf_id
  integer prs_sfc_id
  integer q_CFC11_id
  integer q_CFC12_id
  integer q_CH4_id
  integer q_CO2_id
  integer q_H2OH2O_id
  integer q_H2OH2O_rcp_q_H2O_id
  integer q_H2O_id
  integer q_N2O_id
  integer q_N2_id
  integer q_NO2_id
  integer q_O2O2_id
  integer q_O2_id
  integer q_O3_id
  integer q_OH_id
  integer qst_H2O_ice_id
  integer qst_H2O_lqd_id
  integer r_CFC11_id
  integer r_CFC12_id
  integer r_CH4_id
  integer r_CO2_id
  integer r_H2OH2O_id
  integer r_H2O_id
  integer r_N2O_id
  integer r_N2_id
  integer r_NO2_id
  integer r_O2O2_id
  integer r_O2_id
  integer r_O3_id
  integer r_OH_id
  integer rds_fct_ice_id
  integer rds_fct_lqd_id
  integer rds_ffc_snw_id
  integer rgh_len_id
  integer scl_hgt_id
  integer sfc_ems_id
  integer slr_zen_ngl_cos_id
  integer slr_zen_ngl_dgr_id
  integer slr_zen_ngl_id
  integer snow_depth_id
  integer spc_heat_mst_air_id
  integer tpt_cls_id
  integer tpt_cls_ntf_id
  integer tpt_dwp_id
  integer tpt_id
  integer tpt_ntf_id
  integer tpt_ntf_snw_id
  integer tpt_pcl_id
  integer tpt_pcl_sfc_id
  integer tpt_ptn_id
  integer tpt_ptn_vrt_id
  integer tpt_sfc_id
  integer tpt_skn_id
  integer tpt_snw_id
  integer tpt_vrt_id
  integer vmr_CFC11_id
  integer vmr_CFC12_id
  integer vmr_CH4_id
  integer vmr_CO2_id
  integer vmr_H2OH2O_id
  integer vmr_H2O_id
  integer vmr_N2O_id
  integer vmr_N2_id
  integer vmr_NO2_id
  integer vmr_O2O2_id
  integer vmr_O2_id
  integer vmr_O3_id
  integer vmr_OH_id
  integer wvl_obs_aer_id
  integer wvl_obs_bga_id
  integer xnt_fac_id
  ! sbc
  integer abl_hgt_id 
  integer kbfs_id     
  integer khfs_id
  integer kqfs_id     
  integer lhflx_id    
  integer obk_len_id  
  integer ric_nbr_id
  integer shflx_id    
  integer tau_id
  integer tke_id
  integer ustar_id    
  integer wnd_dir_id
  integer wnd_shr_id
  integer wnd_spd_id
  integer wspd10m_id
  integer wstar_id 
  ! !sbc
  ! Keep coordinate and time variables in double precision
  real(selected_real_kind(p=12))::cnc_O2_npl_N2_clm ! [mlc2 m-5] Column total number concentration times N2 number path
  real(selected_real_kind(p=12))::cnc_O2_npl_O2_clm ! [mlc2 m-5] Column total O2 number concentration times O2 number path
  real(selected_real_kind(p=12))::lat      ! [rdn] Latitude
  real(selected_real_kind(p=12))::lat_cos  ! [frc] Cosine of latitude
  real(selected_real_kind(p=12))::lat_dgr  ! [dgr] Latitude
  real(selected_real_kind(p=12))::lcl_time_hr ! [hr] Local time hour
  real(selected_real_kind(p=12))::lcl_yr_day ! [day] Local year day
  real(selected_real_kind(p=12))::slr_zen_ngl ! [rdn] Solar zenith angle
  real(selected_real_kind(p=12))::slr_zen_ngl_cos ! [frc] Solar zenith angle cosine
  real(selected_real_kind(p=12))::slr_zen_ngl_dgr ! [dgr] Solar zenith angle 

  real CAPE ![J kg-1] Convective Available Potential Energy
  real CINE ![J kg-1] Convective Inhibition Energy
  real CO2_vmr_clm
  real N2O_vmr_clm
  real CH4_vmr_clm
  real CFC11_vmr_clm
  real CFC12_vmr_clm
  real alb_sfc_NIR_drc
  real alb_sfc_NIR_dff
  real alb_sfc
  real alb_sfc_vsb_drc
  real alb_sfc_vsb_dff
  real alt_cld_btm
  real alt_cld_mid
  real alt_cld_thick
  real alt_cld_top
  real dns_O2_mpl_N2_clm    ! [kg2 m-5] Column total O2 mass concentration times N2 mass path
  real dns_O2_mpl_O2_clm    ! [kg2 m-5] Column total O2 mass concentration times O2 mass path
  real dns_aer
  real dns_bga
  real ext_cff_mss_aer
  real ext_cff_mss_bga
  real frc_ice_ttl
  real frc_str_zen_ngl_sfc
  real lev_LCL ![m] Lifting Condensation Level
  real lev_LFC ![m] Level of Free Convection
  real lev_LNB ![m] Level of Neutral Buoyancy
  real mpc_CO2
  real mpc_CH4
  real mpc_N2O
  real mpc_CFC11
  real mpc_CFC12
  real mpc_CWP
  real mpc_H2O
  real mpc_H2OH2O
  real mpc_IWP
  real mpc_LWP
  real mpc_N2
  real mpc_NO2
  real mpc_O2
  real mpc_O2O2
  real mpc_O3
  real mpc_O3_DU
  real mpc_OH
  real mpc_aer
  real mpc_bga
  real mpc_dry_air
  real mpc_mst_air
  real npc_CO2
  real npc_CH4
  real npc_N2O
  real npc_CFC11
  real npc_CFC12
  real npc_H2O
  real npc_H2OH2O
  real npc_N2
  real npc_NO2
  real npc_O2
  real npc_O2O2
  real npc_O3
  real npc_OH
  real npc_dry_air
  real npc_mst_air
  real odxc_obs_aer
  real odxc_obs_bga
  real oro                  ! land/ocean/sea ice flag
  real prs_cld_btm
  real prs_cld_mid
  real prs_cld_thick
  real prs_cld_top
  real prs_sfc
  real rgh_len
  real sfc_ems
  real snow_depth
  real tpt_pcl_sfc ! [K] Air Parcel Surface Temp (lowest 100mb)
  real tpt_sfc
  real tpt_skn
  real wvl_obs_aer
  real wvl_obs_bga
  real xnt_fac
  
  ! Local variables
  real(selected_real_kind(p=12))::pi
  !  real(selected_real_kind(p=12))::double_foo
  
  integer aer_lvl_nbr
  integer bga_lvl_nbr
  integer cld_lvl_nbr
  integer cld_btm_lvl
  integer cld_top_lvl
  
  real prs_top ! [Pa] Pressure at top of "atmosphere"
  real H2O_tune_factor
  real NO2_vmr
  real O3_tune_factor
  real OH_vmr
  real alb_NIR_VIS_rat
  real float_foo
  real mpc_CWP_cmd_ln
  real lon_dgr_m180_p180    ! [dgr]
  real odxc_obs_aer_cmd_ln
  real odxc_obs_bga_cmd_ln
  real PDF_avg_wnd_spd      ! PDF Average wind speed
  real PDF_spd_var          ! PDF Wind spd variability low=1.05;avg=0.94;high=0.83
  real prs_btm_bga
  real prs_top_bga
  real sat_vpr_ice
  real sat_vpr_lqd
  
  !***********************************************************************
  ! netCDF time variables
  integer gmt_day_id
  integer gmt_doy_id
  integer gmt_hr_id
  integer gmt_mnt_id
  integer gmt_mth_id
  integer gmt_sec_id
  integer gmt_ydy_id
  integer gmt_yr_id
  integer lon_sec_id
  integer lmt_day_id
!  integer lmt_doy_id
  integer lmt_hr_id
  integer lmt_mnt_id
  integer lmt_mth_id
  integer lmt_sec_id
  integer lmt_ydy_id
  integer lmt_yr_id
  integer lon_dgr_id
  integer lon_id
  integer ltst_day_id
  integer ltst_doy_id
  integer ltst_hr_id
  integer ltst_mnt_id
  integer ltst_mth_id
  integer ltst_sec_id
  integer ltst_ydy_id
  integer ltst_yr_id
  integer slr_crd_gmm_dgr_id
  integer time_lmt_id
  integer time_ltst_id
  integer time_unix_id
  integer eqn_time_sec_id
  
  character(32) gmt_sng      ! [sng] GMT formatted as Day Mth DD HH:MM:SS TZ YYYY
  character(32) lmt_sng      ! [sng] LMT formatted as Day Mth DD HH:MM:SS TZ YYYY
  character(32) ltst_sng     ! [sng] LTST formatted as Day Mth DD HH:MM:SS TZ YYYY
  
  real(selected_real_kind(p=12))::gmt_doy
  real(selected_real_kind(p=12))::lon_sec
!  real(selected_real_kind(p=12))::lmt_doy
  real(selected_real_kind(p=12))::lon
  real(selected_real_kind(p=12))::lon_dgr
  real(selected_real_kind(p=12))::ltst_doy
  real(selected_real_kind(p=12))::time_lmt ! [s] Seconds between 1969 and LMT of simulation
  real(selected_real_kind(p=12))::time_ltst ! [s] Seconds between 1969 and LTST of simulation
  real(selected_real_kind(p=12))::time_unix ! [s] Seconds between 1969 and GMT of simulation
  
  integer gmt_day
  integer gmt_hr
  integer gmt_mnt
  integer gmt_mth
  integer gmt_sec
  integer gmt_ydy
  integer gmt_yr
  integer lmt_day
  integer lmt_hr
  integer lmt_mnt
  integer lmt_mth
  integer lmt_sec
  integer lmt_ydy
  integer lmt_yr
  integer ltst_day
  integer ltst_hr
  integer ltst_mnt
  integer ltst_mth
  integer ltst_sec
  integer ltst_ydy
  integer ltst_yr
  
  real eqn_time_sec
  real slr_crd_gmm_dgr
  
  ! JJM slr_crd variables
  integer slr_azi_dgr_id
  integer slr_cst_id
  integer slr_dcl_dgr_id
  integer slr_dmt_dgr_id
  integer slr_dst_au_id
  integer slr_elv_dgr_id
  integer slr_flx_TOA_id
  integer slr_flx_nrm_TOA_id
  integer slr_hr_ngl_dgr_id
  integer slr_rfr_ngl_dgr_id
  integer slr_rgt_asc_dgr_id
  
  real slr_azi_dgr
  real slr_cst
  real slr_dcl_dgr
  real slr_dmt_dgr
  real slr_dst_au
  real slr_elv_dgr
  real slr_flx_TOA
  real slr_flx_nrm_TOA
  real slr_hr_ngl_dgr
  real slr_rfr_ngl_dgr
  real slr_rgt_asc_dgr
  
  ! Local time variables
  real gmt_hr_dcm
  real slr_crd_gmm
  character(4) tz_sng        ! [sng] Current time representation, e.g., GMT
  real(selected_real_kind(p=12))::eqn_time_day
  real(selected_real_kind(p=12))::eqn_time_mnt
  real(selected_real_kind(p=12))::gtst_doy
  integer time_unix_yr_srt  ! [s] Seconds between 1969 and beginning of year
  interface
     subroutine gmt2unix( & ! [fnc] Compute UNIX time (seconds since 1969) of given GMT date-time
          yr, & ! I [yr] Christian Year
          mth, & ! I [mth] 1-based month of year [1..12]
          day, & ! I [day] 1-based day of month [1..31]
          hr, & ! I [hr] 0-based, hour of day [0..23]
          mnt, & ! I [mnt] 0-based minute of second [0..59]
          sec, & ! I [s] 0-based second of minute [0..59]
          time_unix_long) ! O [s] Pointer to UNIX time
       integer,intent(in)::yr ! I [yr] Christian Year
       integer,intent(in)::mth ! I [mth] 1-based month of year [1..12]
       integer,intent(in)::day ! I [day] 1-based day of month [1..31]
       integer,intent(in)::hr ! I [hr] 0-based, hour of day [0..23]
       integer,intent(in)::mnt ! I [mnt] 0-based minute of second [0..59]
       integer,intent(in)::sec ! I [s] 0-based second of minute [0..59]
       integer,intent(out)::time_unix_long ! O [s] Pointer to UNIX time
     end subroutine gmt2unix
     subroutine unix2gmt( & ! [fnc] Convert UNIX time (seconds since 1969) to GMT date-time
          time_unix, & ! I [s] Pointer to UNIX time stored in double precision
          yr, & ! O [yr] Christian Year
          ydy, & ! O [day] 1-based day of year [1..366]
          mth, & ! O [mth] 1-based month of year [1..12]
          day, & ! O [day] 1-based day of month [1..31]
          hr, & ! O [hr] 0-based, hour of day [0..23]
          mnt, & ! O [mnt] 0-based minute of hour [0..59]
          sec) ! O [s] 0-based second of minute [0..59]
       real(selected_real_kind(p=12)),intent(in)::time_unix ! I [s] Pointer to UNIX time stored in double precision
       integer,intent(out)::yr ! O [yr] Christian Year
       integer,intent(out)::ydy ! O [day] 1-based day of year [1..366]
       integer,intent(out)::mth ! O [mth] 1-based month of year [1..12]
       integer,intent(out)::day ! O [day] 1-based day of month [1..31]
       integer,intent(out)::hr ! O [hr] 0-based, hour of day [0..23]
       integer,intent(out)::mnt ! O [mnt] 0-based minute of hour [0..59]
       integer,intent(out)::sec ! O [s] 0-based second of minute [0..59]
     end subroutine unix2gmt
     subroutine unix2gmt_sng( & ! [fnc] Create standard date-time string for given decimal time and time-zone
          time_unix, & ! I [s] Pointer to UNIX time_t
          bfr, & ! O [sng] Character buffer to hold output date-time string
          tz_nm) ! I/O [sng] String to use for time-zone in output date-time string
          real(selected_real_kind(p=12)),intent(in)::time_unix ! I [s] Pointer to UNIX time_t
          character(len=*),intent(out)::bfr ! O [sng] Character buffer to hold output date-time string
          character(len=*),intent(inout)::tz_nm ! I/O [sng] String to use for time-zone in output date-time string
     end subroutine unix2gmt_sng
  end interface

  ! Snow property defaults
  real::dpt_dlt_snw_dfl=0.02 ! [m] Snow layer thickness
  real::tpt_snw_dfl=273.0 ! [K] Snow temperature
  real::dns_snw_dfl=0.1 ! [kg m-3] Snow density
  real::mmr_mpr_snw_dfl=0.0 ! [kg mpr/kg snow] Mass mixing ratio of impurities in snow, default
  real::rds_ffc_snw_dfl=100.0e-6 ! [m] Snow effective radius, default
  
  ! Main code
  
  ! Initialize default values
  AFGL_TXT_INPUT=.false.
  CLM_NC_INPUT=.false.
  CLM_TXT_INPUT=.true.
  CRM2=.true.               ! Does CRM text file have CRM 2.x fields?
  H2O_tune_factor=1.0
  O3_tune_factor=1.0
  cmd_ln_mpc_CWP=.false.
  cmd_ln_odxc_obs_aer=.false.
  cmd_ln_odxc_obs_bga=.false.
  dbg_lvl=0
  drc_in='/data/zender/aca'//nlc ! [sng] Input directory
  drc_out='/data/zender/aca'//nlc ! [sng] Output directory
  exit_status=0
  fl_aer='aer_mineral_dust.nc'//nlc
  fl_bga='aer_h2so4_215K.nc'//nlc
  fl_in='mls_cld.txt'//nlc
  fl_out='clm.nc'//nlc
  fl_txt='clm.txt'//nlc
  flg_aer=.false.
  flg_MFRSR_obs=.false.
  flg_bga=.false.
  flg_snw=.false. ! [flg] Incorporate multi-layer snow model
  flg_wnd=.false.
  force_cld_lvl_2b_sat=.true.  ! [flg] Force cloudy layers to be saturated
  force_sat_lvl_2b_cld=.false.
  force_cld_lvl_2b_overcast=.false. ! [flg] Force cloudy layers to be overcast (100% cloud fraction)
  odxc_obs_bga=0.006        ! Optical depth of background aerosol at specified wavelength
  PDF_avg_wnd_spd=10.0      ! Average wind speed for PDF
  PDF_bin_nbr=25            ! Number of PDF bins
  PDF_spd_var=0.94          ! PDF Wind spd variability low=1.05;avg=0.94;high=0.83
  pi=4.0*atan(1.0)
  prs_top=0.0               ! [Pa] Pressure at top of "atmosphere"
  rcd=nf90_noerr              ! nf90_noerr == 0
  wvl_obs_aer=0.5e-6        ! [m] wavelength at which column optical depth of aer is specified
  wvl_obs_bga=0.5e-6        ! [m] wavelength at which column optical depth of bga is specified
  zen_ngl_dpn_alb=.false.

  ! Retrieve command line arguments
  call date_time_get(lcl_date_time)
  call ftn_cmd_ln_sng(cmd_ln)
  call ftn_prg_ID_mk(CVS_Id,CVS_Revision,CVS_Date,prg_ID)
  write (0,'(a)') prg_ID(1:ftn_strlen(prg_ID))
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
        else if (opt_sng == 'lev_snw') then ! [nbr] Number of snow layers
           call ftn_arg_get(arg_idx,arg_val,lev_snw_nbr)
           if (lev_snw_nbr+lev_nbr > lev_nbr_max) stop 'lev_snw_nbr+lev_nbr > lev_nbr_max'
           if (lev_snw_nbr > 0) flg_snw=.true. ! [flg] Multi-layer snow model
        else if (opt_sng == 'lev_atm') then ! [nbr] Number of atmosphere layers
           call ftn_arg_get(arg_idx,arg_val,lev_nbr)
           if (lev_nbr > lev_nbr_max) stop 'lev_nbr > lev_nbr_max'
        else if (opt_sng == 'prs_top') then ! [Pa] Pressure at top of "atmosphere"
           call ftn_arg_get(arg_idx,arg_val,prs_top) 
        else if (opt_sng == 'snw') then ! [flg] Multi-layer snow model
           flg_snw=.true.
        else if (opt_sng == 'wbl_c') then
           call ftn_arg_get(arg_idx,arg_val,PDF_spd_var) ! PDF Wind spd variability low=1.05;avg=0.94;high=0.83
        else if (opt_sng == 'wnd_nbr') then
           call ftn_arg_get(arg_idx,arg_val,PDF_bin_nbr) !
        else if (opt_sng == 'wnd_znl_mdp') then
           call ftn_arg_get(arg_idx,arg_val,PDF_avg_wnd_spd) !
        else                ! Option not recognized
           arg_idx=arg_idx-1 ! [idx] Counting index
           call ftn_getarg_err(arg_idx,arg_val) ! [sbr] Error handler for getarg()
        endif               ! endif option is recognized
        ! Jump to top of while loop
        cycle loop_while_options ! C, F77, and F90 use "continue", "goto", and "cycle"
     endif if_dbl_dsh                  ! endif long option
     ! Handle short options
     if_sgl_dsh: if (dsh_key == '-3') then
        fl_out_fmt=nf90_format_classic ! [enm] Output file format
     else if (dsh_key == '-4') then
        fl_out_fmt=nf90_format_netcdf4 ! [enm] Output file format
     else if (dsh_key == '-A') then
        flg_aer=.not.flg_aer
     else if(dsh_key == '-a') then
        call ftn_arg_get(arg_idx,arg_val,fl_aer)
     else if(dsh_key == '-B') then
        flg_bga=.not.flg_bga
     else if(dsh_key == '-b') then
        call ftn_arg_get(arg_idx,arg_val,fl_bga)
     else if(dsh_key == '-C') then
        CRM2=.not.CRM2
     else if(dsh_key == '-D') then
        call ftn_arg_get(arg_idx,arg_val,dbg_lvl) ! [enm] Debugging level
     else if(dsh_key == '-f') then
        call ftn_arg_get(arg_idx,arg_val,float_foo)
     else if(dsh_key == '-h') then
        call ftn_arg_get(arg_idx,arg_val,H2O_tune_factor)
     else if(dsh_key == '-H') then
        call ftn_arg_get(arg_idx,arg_val,OH_vmr)
     else if(dsh_key == '-I') then
        AFGL_TXT_INPUT=.true.
        CLM_NC_INPUT=.false.
        CLM_TXT_INPUT=.false.
     else if(dsh_key == '-i') then
        call ftn_arg_get(arg_idx,arg_val,fl_in)
     else if(dsh_key == '-J') then
        cmd_ln_odxc_obs_aer=.not.cmd_ln_odxc_obs_aer
        call ftn_arg_get(arg_idx,arg_val,odxc_obs_aer_cmd_ln)
     else if(dsh_key == '-j') then
        cmd_ln_odxc_obs_bga=.not.cmd_ln_odxc_obs_bga
        call ftn_arg_get(arg_idx,arg_val,odxc_obs_bga_cmd_ln)
     else if(dsh_key == '-l') then
        call ftn_arg_get(arg_idx,arg_val,lev_nbr)
        if (lev_nbr > lev_nbr_max) stop 'lev_nbr > lev_nbr_max'
     else if(dsh_key == '-m') then
        cmd_ln_mpc_CWP=.not.cmd_ln_mpc_CWP
        call ftn_arg_get(arg_idx,arg_val,mpc_CWP_cmd_ln)
     else if(dsh_key == '-M') then
        flg_MFRSR_obs=.not.flg_MFRSR_obs
     else if(dsh_key == '-N') then
        call ftn_arg_get(arg_idx,arg_val,NO2_vmr)
     else if(dsh_key == '-n') then
        CLM_NC_INPUT=.true.
        AFGL_TXT_INPUT=.false.
        CLM_TXT_INPUT=.false.
     else if(dsh_key == '-o') then
        call ftn_arg_get(arg_idx,arg_val,fl_out)
     else if(dsh_key == '-O') then
        force_cld_lvl_2b_overcast=.true.
     else if(dsh_key == '-s') then
        ! force_cld_lvl_2b_sat and force_sat_lvl_2b_cld must not both be true
        force_cld_lvl_2b_sat=.true.
     else if(dsh_key == '-S') then
        force_sat_lvl_2b_cld=.true.
     else if(dsh_key == '-t') then
        call ftn_arg_get(arg_idx,arg_val,fl_txt)
     else if(dsh_key == '-U') then
        flg_wnd=.not.flg_wnd   
     else if(dsh_key == '-v') then
        write (6,'(a)') CVS_Id
        goto 1000           ! Goto exit with error status
     else if(dsh_key == '-W') then
        call ftn_arg_get(arg_idx,arg_val,wvl_obs_aer)
     else if(dsh_key == '-w') then
        call ftn_arg_get(arg_idx,arg_val,wvl_obs_bga)
     else if(dsh_key == '-z') then
        call ftn_arg_get(arg_idx,arg_val,O3_tune_factor)
     else if(dsh_key == '-Z') then
        zen_ngl_dpn_alb=.not.zen_ngl_dpn_alb
     else                   ! Option not recognized
        arg_idx=arg_idx-1   ! [idx] Counting index
        call ftn_getarg_err(arg_idx,arg_val) ! [sbr] Error handler for getarg()
     endif if_sgl_dsh       ! endif arg_val
  end do loop_while_options ! end while (arg_idx <= arg_nbr)
  
  if (CLM_NC_INPUT) then ! Input file is in CLM_NC_INPUT format
     ! Assume lev,prs,T,q_H2O are available in input file in SI units
     ! Ingest fl_in
     rcd=nf90_wrp_open(fl_in,nf90_nowrite,nc_id)
     ! Get dimension IDs
     rcd=nf90_wrp_inq_dimid(nc_id,'lev',lev_dmn_id)
     ! Get dimension sizes
     rcd=nf90_wrp(nf90_inquire_dimension(nc_id,lev_dmn_id,len=lev_nbr),sbr_nm//": inquire_dim lev")
  end if ! !CLM_NC_INPUT

  ! Compute any quantities that might depend on command line input
  levp_nbr=lev_nbr+1 ! [nbr] dimension size
  levp_snw_nbr=lev_snw_nbr+1 ! [nbr] dimension size
  
  ! Allocate space for dynamic arrays
  allocate(PDF_bin_spd(PDF_bin_nbr),stat=rcd)
  allocate(PDF_bin_wgt(PDF_bin_nbr),stat=rcd)
  allocate(RH(lev_nbr),stat=rcd)
  allocate(RH_ice(lev_nbr),stat=rcd)
  allocate(RH_lqd(lev_nbr),stat=rcd)
  allocate(alt(lev_nbr),stat=rcd)
  allocate(alt_dlt(lev_nbr),stat=rcd)
  allocate(alt_ntf(levp_nbr),stat=rcd)
  allocate(bga_lvl_flg(lev_nbr),stat=rcd)
  allocate(cld_frc(lev_nbr),stat=rcd)
  allocate(cnc_CFC11(lev_nbr),stat=rcd)
  allocate(cnc_CFC12(lev_nbr),stat=rcd)
  allocate(cnc_CH4(lev_nbr),stat=rcd)
  allocate(cnc_CO2(lev_nbr),stat=rcd)
  allocate(cnc_H2O(lev_nbr),stat=rcd)
  allocate(cnc_H2OH2O(lev_nbr),stat=rcd)
  allocate(cnc_N2(lev_nbr),stat=rcd)
  allocate(cnc_N2O(lev_nbr),stat=rcd)
  allocate(cnc_NO2(lev_nbr),stat=rcd)
  allocate(cnc_O2(lev_nbr),stat=rcd)
  allocate(cnc_O2O2(lev_nbr),stat=rcd)
  allocate(cnc_O2_cnc_N2(lev_nbr),stat=rcd) ! [mlc2 m-6] O2 number concentration times N2 number concentration
  allocate(cnc_O2_cnc_O2(lev_nbr),stat=rcd) ! [mlc2 m-6] O2 number concentration squared
  allocate(cnc_O2_npl_N2(lev_nbr),stat=rcd) ! [mlc2 m-5] O2 number concentration times N2 number path
  allocate(cnc_O2_npl_O2(lev_nbr),stat=rcd) ! [mlc2 m-5] O2 number concentration times O2 number path
  allocate(cnc_O2_npl_O2_clm_frc(lev_nbr),stat=rcd) ! [frc] Fraction of column O2-O2 at or above each level
  allocate(cnc_O3(lev_nbr),stat=rcd)
  allocate(cnc_OH(lev_nbr),stat=rcd)
  allocate(cnc_dry_air(lev_nbr),stat=rcd)
  allocate(cnc_mst_air(lev_nbr),stat=rcd)
  allocate(dns_CFC11(lev_nbr),stat=rcd)
  allocate(dns_CFC12(lev_nbr),stat=rcd)
  allocate(dns_CH4(lev_nbr),stat=rcd)
  allocate(dns_CO2(lev_nbr),stat=rcd)
  allocate(dns_H2O(lev_nbr),stat=rcd)
  allocate(dns_H2OH2O(lev_nbr),stat=rcd)
  allocate(dns_N2(lev_nbr),stat=rcd)
  allocate(dns_N2O(lev_nbr),stat=rcd)
  allocate(dns_NO2(lev_nbr),stat=rcd)
  allocate(dns_O2(lev_nbr),stat=rcd)
  allocate(dns_O2O2(lev_nbr),stat=rcd)
  allocate(dns_O2_dns_N2(lev_nbr),stat=rcd) ! [kg2 m-6] O2 mass concentration times N2 mass concentration
  allocate(dns_O2_dns_O2(lev_nbr),stat=rcd) ! [kg2 m-6] O2 mass concentration squared
  allocate(dns_O2_mpl_N2(lev_nbr),stat=rcd) ! [kg2 m-5] O2 mass concentration times N2 mass path
  allocate(dns_O2_mpl_O2 (lev_nbr),stat=rcd) ! [kg2 m-5] O2 mass concentration times O2 mass path
  allocate(dns_O3(lev_nbr),stat=rcd)
  allocate(dns_OH(lev_nbr),stat=rcd)
  allocate(dns_dry_air(lev_nbr),stat=rcd)
  allocate(dns_mst_air(lev_nbr),stat=rcd)
  allocate(frc_ice(lev_nbr),stat=rcd)
  allocate(gas_cst_mst_air(lev_nbr),stat=rcd)
  allocate(grv(lev_nbr),stat=rcd)
  allocate(lev(lev_nbr),stat=rcd)     ! coordinate variable
  allocate(levp(levp_nbr),stat=rcd)   ! coordinate variable
  allocate(mmw_mst_air(lev_nbr),stat=rcd)
  allocate(mpl_CFC11(lev_nbr),stat=rcd)
  allocate(mpl_CFC12(lev_nbr),stat=rcd)
  allocate(mpl_CH4(lev_nbr),stat=rcd)
  allocate(mpl_CO2(lev_nbr),stat=rcd)
  allocate(mpl_CWP(lev_nbr),stat=rcd)
  allocate(mpl_H2O(lev_nbr),stat=rcd)
  allocate(mpl_H2OH2O(lev_nbr),stat=rcd)
  allocate(mpl_IWP(lev_nbr),stat=rcd)
  allocate(mpl_LWP(lev_nbr),stat=rcd)
  allocate(mpl_N2(lev_nbr),stat=rcd)
  allocate(mpl_N2O(lev_nbr),stat=rcd)
  allocate(mpl_NO2(lev_nbr),stat=rcd)
  allocate(mpl_O2(lev_nbr),stat=rcd)
  allocate(mpl_O2O2(lev_nbr),stat=rcd)
  allocate(mpl_O3(lev_nbr),stat=rcd)
  allocate(mpl_OH(lev_nbr),stat=rcd)
  allocate(mpl_aer(lev_nbr),stat=rcd)
  allocate(mpl_bga(lev_nbr),stat=rcd)
  allocate(mpl_dry_air(lev_nbr),stat=rcd)
  allocate(mpl_mst_air(lev_nbr),stat=rcd)
  allocate(npl_CFC11(lev_nbr),stat=rcd)
  allocate(npl_CFC12(lev_nbr),stat=rcd)
  allocate(npl_CH4(lev_nbr),stat=rcd)
  allocate(npl_CO2(lev_nbr),stat=rcd)
  allocate(npl_H2O(lev_nbr),stat=rcd)
  allocate(npl_H2OH2O(lev_nbr),stat=rcd)
  allocate(npl_N2(lev_nbr),stat=rcd)
  allocate(npl_N2O(lev_nbr),stat=rcd)
  allocate(npl_NO2(lev_nbr),stat=rcd)
  allocate(npl_O2(lev_nbr),stat=rcd)
  allocate(npl_O2O2(lev_nbr),stat=rcd)
  allocate(npl_O3(lev_nbr),stat=rcd)
  allocate(npl_OH(lev_nbr),stat=rcd)
  allocate(npl_dry_air(lev_nbr),stat=rcd)
  allocate(npl_mst_air(lev_nbr),stat=rcd)
  allocate(nrg_dry(lev_nbr),stat=rcd)
  allocate(nrg_mst(lev_nbr),stat=rcd)
  allocate(odxl_obs_aer(lev_nbr),stat=rcd)
  allocate(odxl_obs_bga(lev_nbr),stat=rcd)
  allocate(oneD_foo(lev_nbr),stat=rcd)
  allocate(pcl_DLR(lev_nbr),stat=rcd)
  allocate(pcl_MLR(lev_nbr),stat=rcd)
  allocate(pcl_PLR(lev_nbr),stat=rcd)
  allocate(ppr_CFC11(lev_nbr),stat=rcd)
  allocate(ppr_CFC12(lev_nbr),stat=rcd)
  allocate(ppr_CH4(lev_nbr),stat=rcd)
  allocate(ppr_CO2(lev_nbr),stat=rcd)
  allocate(ppr_H2O(lev_nbr),stat=rcd)
  allocate(ppr_H2OH2O(lev_nbr),stat=rcd)
  allocate(ppr_N2(lev_nbr),stat=rcd)
  allocate(ppr_N2O(lev_nbr),stat=rcd)
  allocate(ppr_NO2(lev_nbr),stat=rcd)
  allocate(ppr_O2(lev_nbr),stat=rcd)
  allocate(ppr_O2O2(lev_nbr),stat=rcd)
  allocate(ppr_O3(lev_nbr),stat=rcd)
  allocate(ppr_OH(lev_nbr),stat=rcd)
  allocate(ppr_dry_air(lev_nbr),stat=rcd)
  allocate(prs(lev_nbr),stat=rcd)
  allocate(prs_dlt(lev_nbr),stat=rcd)
  allocate(prs_ntf(levp_nbr),stat=rcd)
  allocate(prs_ttl(lev_nbr),stat=rcd)
  allocate(q_CFC11(lev_nbr),stat=rcd)
  allocate(q_CFC12(lev_nbr),stat=rcd)
  allocate(q_CH4(lev_nbr),stat=rcd)
  allocate(q_CO2(lev_nbr),stat=rcd)
  allocate(q_H2O(lev_nbr),stat=rcd)
  allocate(q_H2OH2O(lev_nbr),stat=rcd)
  allocate(q_H2OH2O_rcp_q_H2O(lev_nbr),stat=rcd)
  allocate(q_N2(lev_nbr),stat=rcd)
  allocate(q_N2O(lev_nbr),stat=rcd)
  allocate(q_NO2(lev_nbr),stat=rcd)
  allocate(q_O2(lev_nbr),stat=rcd)
  allocate(q_O2O2(lev_nbr),stat=rcd)
  allocate(q_O3(lev_nbr),stat=rcd)
  allocate(q_OH(lev_nbr),stat=rcd)
  allocate(qst_H2O_ice(lev_nbr),stat=rcd)
  allocate(qst_H2O_lqd(lev_nbr),stat=rcd)
  allocate(r_CFC11(lev_nbr),stat=rcd)
  allocate(r_CFC12(lev_nbr),stat=rcd)
  allocate(r_CH4(lev_nbr),stat=rcd)
  allocate(r_CO2(lev_nbr),stat=rcd)
  allocate(r_H2O(lev_nbr),stat=rcd)
  allocate(r_H2OH2O(lev_nbr),stat=rcd)
  allocate(r_N2(lev_nbr),stat=rcd)
  allocate(r_N2O(lev_nbr),stat=rcd)
  allocate(r_NO2(lev_nbr),stat=rcd)
  allocate(r_O2(lev_nbr),stat=rcd)
  allocate(r_O2O2(lev_nbr),stat=rcd)
  allocate(r_O3(lev_nbr),stat=rcd)
  allocate(r_OH(lev_nbr),stat=rcd)
  allocate(rds_fct_ice(lev_nbr),stat=rcd)
  allocate(rds_fct_lqd(lev_nbr),stat=rcd)
  allocate(scl_hgt(lev_nbr),stat=rcd)
  allocate(spc_heat_mst_air(lev_nbr),stat=rcd)
  allocate(tpt(lev_nbr),stat=rcd)
  allocate(tpt_cls(lev_nbr),stat=rcd)
  allocate(tpt_cls_ntf(levp_nbr),stat=rcd)
  allocate(tpt_dwp(lev_nbr),stat=rcd)
  allocate(tpt_ntf(levp_nbr),stat=rcd)
  allocate(tpt_pcl(lev_nbr),stat=rcd)
  allocate(tpt_ptn(lev_nbr),stat=rcd)
  allocate(tpt_ptn_vrt(lev_nbr),stat=rcd)
  allocate(tpt_vrt(lev_nbr),stat=rcd)
  allocate(vmr_CFC11(lev_nbr),stat=rcd)
  allocate(vmr_CFC12(lev_nbr),stat=rcd)
  allocate(vmr_CH4(lev_nbr),stat=rcd)
  allocate(vmr_CO2(lev_nbr),stat=rcd)
  allocate(vmr_H2O(lev_nbr),stat=rcd)
  allocate(vmr_H2OH2O(lev_nbr),stat=rcd)
  allocate(vmr_N2(lev_nbr),stat=rcd)
  allocate(vmr_N2O(lev_nbr),stat=rcd)
  allocate(vmr_NO2(lev_nbr),stat=rcd)
  allocate(vmr_O2(lev_nbr),stat=rcd)
  allocate(vmr_O2O2(lev_nbr),stat=rcd)
  allocate(vmr_O3(lev_nbr),stat=rcd)
  allocate(vmr_OH(lev_nbr),stat=rcd)
  if (flg_snw) then
     allocate(dns_snw(lev_snw_nbr),stat=rcd) ! [kg m-3] Snow density
     allocate(dpt_dlt_snw(lev_snw_nbr),stat=rcd) ! [m] Snow layer thickness
     allocate(dpt_ntf_snw(levp_snw_nbr),stat=rcd) ! [m] Snow depth interfaces
     allocate(dpt_snw(lev_snw_nbr),stat=rcd) ! [m] Snow depth
     allocate(foo_snw(lev_snw_nbr),stat=rcd) ! [K] Snow foo
     allocate(lev_snw(lev_snw_nbr),stat=rcd) ! [m] Snow depth
     allocate(levp_snw(levp_snw_nbr),stat=rcd) ! [m] Snow interface depth
     allocate(mmr_mpr_snw(lev_snw_nbr),stat=rcd) ! [kg mpr/kg snow] Mass mixing ratio of impurities in snow
     allocate(rds_ffc_snw(lev_snw_nbr),stat=rcd) ! [m] Snow effective radius
     allocate(tpt_ntf_snw(levp_snw_nbr),stat=rcd) ! [K] Snow temperature interfaces
     allocate(tpt_snw(lev_snw_nbr),stat=rcd) ! [K] Snow temperature
  endif ! !flg_snw
  ! sbc optional wind variables
  allocate(ric_nbr(lev_nbr),stat=rcd)
  allocate(tke(levp_nbr),stat=rcd)  
  allocate(wnd_dir(lev_nbr),stat=rcd)
  allocate(wnd_shr(lev_nbr),stat=rcd)
  allocate(wnd_spd(lev_nbr),stat=rcd)
  ! !sbc
  ! Compute quantities that may depend on command line input
  ! Prepend user-specified path, if any, to input data file names
  if (ftn_strlen(drc_in) > 0) then
     call ftn_drcpfx(drc_in,fl_aer) ! [sng] Aerosol file
     call ftn_drcpfx(drc_in,fl_bga) ! [sng] Background aerosol file
     call ftn_drcpfx(drc_in,fl_in) ! [sng] Input file
  endif                     ! endif drc_in
  ! Prepend user-specified path, if any, to output data file names
  if (ftn_strlen(drc_out) > 0) call ftn_drcpfx(drc_out,fl_out) ! [sng] Output file
  if (ftn_strlen(drc_out) > 0) call ftn_drcpfx(drc_out,fl_txt) ! [sng] CLM text output file
  
  ! force_cld_lvl_2b_sat and force_sat_lvl_2b_cld must not both be true
  if (force_sat_lvl_2b_cld) force_cld_lvl_2b_sat=.false.
  
  if (CLM_TXT_INPUT) then 
     write(0,'(a)') 'Input presumed to be CLM file in text format' 
     if (CRM2) then 
        write(0,'(a)') 'CRM 2.x fields read from input file' 
     else                   ! not CRM2
        write(0,'(a)') 'CRM 2.x fields set to defaults'
     endif                 ! not CRM2
     if (flg_wnd) then
        write(0,'(a)') 'Extra text columns containing layer winds are expected' 
     else                 ! not flg_wnd
        write(0,'(a)') 'Text input presumed to omit layer winds'
     endif             ! !flg_wnd
     if (flg_aer) then 
        write(0,'(a)') 'Extra text column containing layer aerosol burden is expected' 
     else                   ! not flg_aer
        write(0,'(a)') 'Text input presumed to omit layer aerosol burden'
     endif                  ! not flg_aer
  else if (CLM_NC_INPUT) then
     write(0,'(a)') 'Input presumed to be CLM file in netCDF format'
  else if (AFGL_TXT_INPUT) then
     write(0,'(a)') 'Input presumed to be AFGL file in ASCII format'
  endif                     ! not CLM_NC_INPUT
  
  if (cmd_ln_odxc_obs_bga) then 
     write(0,'(a)') 'Imposing background aerosol' 
  else                      ! not cmd_ln_odxc_obs_bga
     write(0,'(a)') 'Not imposing background aerosol'
  endif                     ! not cmd_ln_odxc_obs_bga
  
  if (force_sat_lvl_2b_cld) then
     write(0,'(a)') 'Forcing saturated levels to be cloudy' 
  else                      ! not force_sat_lvl_2b_cld
     write(0,'(a)') 'Not forcing saturated levels to be cloudy'
  endif                     ! not force_sat_lvl_2b_cld
  
  if (force_cld_lvl_2b_sat) then
     write(0,'(a)') 'Forcing cloudy and snowy levels to be saturated' 
  else                      ! not force_cld_lvl_2b_sat
     write(0,'(a)') 'Not forcing cloudy levels to be saturated' 
  endif                     ! not force_cld_lvl_2b_sat
  
  ! Set column profile defaults:
  ! CRM 1.x input files overwrite most of these
  ! CRM 2.x input files overwrite all of these
  ! AFGL 2.x input files overwrite a few of these
  ! netCDF CLM input files overwrite some of these, depending on the origin of the netCDF CLM file 
  prf_sng='Default profile string'//char(0)
  prf_snw_sng='Default snow profile string'//char(0)
  lcl_yr_day=172.5          ! [day] Summer solstice
  lat=pi*40.0/180.0         ! [rdn] Baseline in Boulder 
  oro=1.0                   ! [flg] Surface type flag (0.0=ocn, 1.0=lnd, 2.0=sea ice)
  sfc_ems=1.0
  rgh_len=0.01              ! [m] Surface aerodynamic roughness (obsolete)
  snow_depth=0.0            ! [m] Snow depth (liq. equiv.)
  alb_sfc_vsb_drc=0.1       ! [frc] Albedo (Vis, direct)
  alb_sfc_vsb_dff=0.1       ! [frc] Albedo (Vis, diffuse)
  alb_sfc_NIR_drc=0.1       ! [frc] Albedo (NIR, direct)
  alb_sfc_NIR_dff=0.1       ! [frc] Albedo (NIR, diffuse)
  frc_str_zen_ngl_sfc=0.0   ! [frc] Fraction strong zenith angle dep. sfc. (obsolete)
  CO2_vmr_clm=3.55e-4       ! [frc] CCM: physics/comvmr.h: co2vmr set in control/preset()
  ! CRM2 inputs
  N2O_vmr_clm=0.311e-6      ! [frc] CCM: physics/comvmr.h: n2ovmr set in control/preset()
  CH4_vmr_clm=1.714e-6      ! [frc] CCM: physics/comvmr.h: ch4vmr set in control/preset()
  CFC11_vmr_clm=0.280e-9    ! [frc] CCM: physics/comvmr.h: f11vmr set in control/preset()
  CFC12_vmr_clm=0.503e-9    ! [frc] CCM: physics/comvmr.h: f12vmr set in control/preset()
  odxc_obs_aer=0.14         ! [frc] Aerosol visible opt. depth
  slr_cst=slr_cst_CCM       ! [W m-2] Solar constant 
  gmt_yr=1995               ! [yr] Year AD
  lon_dgr=262.5159          ! [dgr] Longitude (degrees East, from 0.0 to 360.0)
  
  ! Snowpack defaults
  if (flg_snw) then
     dpt_dlt_snw(:)=dpt_dlt_snw_dfl ! [m] Snow layer thickness
     tpt_snw(:)=tpt_snw_dfl ! [K] Snow temperature
     dns_snw(:)=dns_snw_dfl ! [kg m-3] Snow density
     mmr_mpr_snw(:)=mmr_mpr_snw_dfl ! [kg mpr/kg snow] Mass mixing ratio of impurities in snow
     rds_ffc_snw(:)=rds_ffc_snw_dfl ! [m] Snow effective radius
     foo_snw(:)=0.0 ! [K] Snow foo
  endif ! !flg_snw

  if (CLM_TXT_INPUT) then
     
     ! Read input
     open (fl_in_unit,file=fl_in,status='old',iostat=rcd)
     if (rcd /= 0 ) write (6,'(a23,1x,a)') 'ERROR opening text file',fl_in(1:ftn_strlen(fl_in))
     
     read (fl_in_unit,'(a80)') lbl
     call ftn_strini(prf_sng) ! [sng] sng(1:len)=NUL
     read (fl_in_unit,'(a)') prf_sng
     call ftn_strnul(prf_sng) ! [sbr] NUL-initialize all characters after LSC
     read (fl_in_unit,'(a80)') lbl
     read (fl_in_unit,*) lcl_yr_day ! [day] Julian day of year (1.5 = Noon, Jan 1; from 1 to 365)
     read (fl_in_unit,*) lat_dgr ! [dgr] Latitude (degrees North, from -90.0 to +90.0)
     read (fl_in_unit,'(a80)') lbl
     
     if (dbg_lvl == dbg_crr) write (6,*) lbl
     
     do idx=1,lev_nbr
        if (flg_aer.AND..NOT.flg_wnd) then 
           read (fl_in_unit,*) &
                lev(idx), &
                prs(idx), &
                tpt(idx), &
                q_H2O(idx), &
                q_O3(idx), &
                cld_frc(idx), &
                mpl_CWP(idx), &
                odxl_obs_aer(idx)
        else if (flg_aer.AND.flg_wnd) then
           read (fl_in_unit,*) &
                lev(idx), &
                prs(idx), &
                tpt(idx), &
                q_H2O(idx), &
                q_O3(idx), &
                cld_frc(idx), &
                mpl_CWP(idx), &
                odxl_obs_aer(idx), &
                wnd_dir(idx), &
                wnd_spd(idx)
        else if (.NOT.flg_aer.AND.flg_wnd) then
           read (fl_in_unit,*) &
                lev(idx), &
                prs(idx), &
                tpt(idx), &
                q_H2O(idx), &
                q_O3(idx), &
                cld_frc(idx), &
                mpl_CWP(idx), &
                wnd_dir(idx), &
                wnd_spd(idx)
           odxl_obs_aer(idx)=0.0 
        else
           read (fl_in_unit,*) &
                lev(idx), &
                prs(idx), &
                tpt(idx), &
                q_H2O(idx), &
                q_O3(idx), &
                cld_frc(idx), &
                mpl_CWP(idx)
           odxl_obs_aer(idx)=0.0
        endif
        odxl_obs_bga(idx)=0.0
     enddo ! end loop over layer
     
     read (fl_in_unit,*) prs_sfc ! Surface pressure [mb]
     read (fl_in_unit,*) tpt_sfc ! Surface air temperature [K]
     read (fl_in_unit,*) tpt_skn ! Ground (skin) temperature [K]
     if (flg_wnd) then 
        read (fl_in_unit,*) q_H2O_sfc ! Ground (skin) mixing ratio [kg kg-1]
     end if ! !flg_wnd
     read (fl_in_unit,*) oro ! Surface type flag (0=ocn, 1=lnd, 2=sea ice)
     read (fl_in_unit,*) rgh_len ! Surface aerodynamic roughness [m] (obsolete)
     read (fl_in_unit,*) snow_depth ! Snow cover liquid water equivalent [m]
     read (fl_in_unit,*) alb_sfc_vsb_drc ! Albedo (Vis, direct)
     read (fl_in_unit,*) alb_sfc_vsb_dff ! Albedo (Vis, diffuse)
     read (fl_in_unit,*) alb_sfc_NIR_drc ! Albedo (NIR, direct)
     read (fl_in_unit,*) alb_sfc_NIR_dff ! Albedo (NIR, diffuse)
     read (fl_in_unit,*) frc_str_zen_ngl_sfc ! Fraction strong zenith angle dep. sfc. (obsolete)
     read (fl_in_unit,*) CO2_vmr_clm ! CO2 volume mixing ratio
     if (CRM2) then
        read (fl_in_unit,*) N2O_vmr_clm ! N2O volume mixing ratio
        read (fl_in_unit,*) CH4_vmr_clm ! CH4 volume mixing ratio
        read (fl_in_unit,*) CFC11_vmr_clm ! CFC11 volume mixing ratio
        read (fl_in_unit,*) CFC12_vmr_clm ! CFC12 volume mixing ratio
        read (fl_in_unit,*) odxc_obs_aer ! Aerosol visible opt. depth
        read (fl_in_unit,*) slr_cst ! Solar constant [W m-2]
        read (fl_in_unit,*) gmt_yr ! Year AD
        read (fl_in_unit,*) lon_dgr ! Longitude (degrees East, from 0.0 to 360.0)
     endif                  ! not CRM2

     if (flg_snw) then
        read (fl_in_unit,'(a80)') lbl ! Label (SCRM ID)
        call ftn_strini(prf_snw_sng) ! [sng] sng(1:len)=NUL
        read (fl_in_unit,'(a)') prf_snw_sng
        call ftn_strnul(prf_snw_sng) ! [sbr] NUL-initialize all characters after LSC
        read (fl_in_unit,'(a80)') lbl ! Label (snow field names)
        read (fl_in_unit,'(a80)') lbl ! Label (snow units)
        do lev_snw_idx=1,lev_snw_nbr
           read (fl_in_unit,*) &
                lev_snw(lev_snw_idx), & ! [idx] Snow layer number (for now)
                dpt_dlt_snw(lev_snw_idx), & ! [cm] Snow layer thickness
                tpt_snw(lev_snw_idx), & ! [K] Snow temperature
                dns_snw(lev_snw_idx), & ! [g cm-3] Snow density
                mmr_mpr_snw(lev_snw_idx), & ! [kg mpr/kg snow] Mass mixing ratio of impurities in snow
                rds_ffc_snw(lev_snw_idx) ! [um] Snow effective radius
        enddo ! end loop over layer
     endif ! !flg_snw
     
     close (fl_in_unit)
     write (0,'(a18,1x,a)') 'Ingested text file',fl_in(1:ftn_strlen(fl_in))
     
     ! Convert input data to SI units where necessary
     do idx=1,lev_nbr
        prs(idx)=prs(idx)*100.0 ! [mb] -> Pa
        mpl_CWP(idx)=mpl_CWP(idx)/1000.0 ! [g m-2] -> [kg m-2]
        ! if (cld_frc(idx) > 0.99999) cld_frc(idx)=.99999
     enddo ! end loop over lev
     if (flg_snw) then
        do lev_snw_idx=1,lev_snw_nbr
           dpt_dlt_snw(lev_snw_idx)=dpt_dlt_snw(lev_snw_idx)*0.01 ! [cm]->[m] Snow layer thickness
           dns_snw(lev_snw_idx)=dns_snw(lev_snw_idx)*1000.0 ! [g cm-3]->[kg m-3] Snow density
           rds_ffc_snw(lev_snw_idx)=rds_ffc_snw(lev_snw_idx)*1.0e-6 ! [um]->[m] Snow effective radius
        enddo ! end loop over lev
     endif ! !flg_snw
     prs_sfc=prs_sfc*100.0  ! [mb] -> Pa
     lat=pi*lat_dgr/180.0   ! [dgr] -> [rdn]
     
  else if (AFGL_TXT_INPUT) then ! input file is in AFGL_TXT_INPUT format
     
     ! Allocate space for dynamic arrays
     allocate(CH4_vmr(lev_nbr),stat=rcd)
     allocate(CH4_vmr_ntf(levp_nbr),stat=rcd)
     allocate(CO_vmr(lev_nbr),stat=rcd)
     allocate(CO_vmr_ntf(levp_nbr),stat=rcd)
     allocate(H2O_vmr(lev_nbr),stat=rcd)
     allocate(H2O_vmr_ntf(levp_nbr),stat=rcd)
     allocate(N2O_vmr(lev_nbr),stat=rcd)
     allocate(N2O_vmr_ntf(levp_nbr),stat=rcd)
     allocate(O3_vmr(lev_nbr),stat=rcd)
     allocate(O3_vmr_ntf(levp_nbr),stat=rcd)
     allocate(cnc_air(lev_nbr),stat=rcd)
     allocate(cnc_air_ntf(levp_nbr),stat=rcd)
     
     ! Read input quantities
     open (fl_in_unit,file=fl_in,status='old',iostat=rcd)
     
     call ftn_strini(prf_sng) ! [sng] sng(1:len)=NUL
     read (fl_in_unit,'(a18,a20)') lbl,prf_sng
     call ftn_strnul(prf_sng) ! [sbr] NUL-initialize all characters after LSC
     read (fl_in_unit,'(a80)') lbl
     read (fl_in_unit,'(a80)') lbl
     read (fl_in_unit,'(a80)') lbl
     read (fl_in_unit,'(a80)') lbl
     
     do levp_idx=1,levp_nbr
        int_foo=levp_nbr-levp_idx+1
        read (fl_in_unit,*) &
             alt_ntf(int_foo), &
             prs_ntf(int_foo), &
             tpt_ntf(int_foo), &
             cnc_air_ntf(int_foo), &
             H2O_vmr_ntf(int_foo), &
             O3_vmr_ntf(int_foo), &
             N2O_vmr_ntf(int_foo), &
             CO_vmr_ntf(int_foo), &
             CH4_vmr_ntf(int_foo)
     enddo
     
     close (fl_in_unit)
     write (0,'(a20,1x,a)') 'Read input data from',fl_in(1:ftn_strlen(fl_in))
     
     ! Convert input data to SI units where necessary
     do idx=1,levp_nbr
        alt_ntf(idx)=alt_ntf(idx)*1000.0 ! [km] -> [m]
        prs_ntf(idx)=prs_ntf(idx)*100.0 ! [mb] -> [Pa]
        cnc_air_ntf(idx)=cnc_air_ntf(idx)/1.0e6 ! [cm-3] -> [m-3]
        H2O_vmr_ntf(idx)=H2O_vmr_ntf(idx)*1.0e-6 ! [ppmv] -> [vmr]
        O3_vmr_ntf(idx)=O3_vmr_ntf(idx)*1.0e-6 ! [ppmv] -> [vmr]
        N2O_vmr_ntf(idx)=N2O_vmr_ntf(idx)*1.0e-6 ! [ppmv] -> [vmr]
        CO_vmr_ntf(idx)=CO_vmr_ntf(idx)*1.0e-6 ! [ppmv] -> [vmr]
        CH4_vmr_ntf(idx)=CH4_vmr_ntf(idx)*1.0e-6 ! [ppmv] -> [vmr]
     enddo                  ! end loop over levp
     
     ! Place input data on layer midpoint surfaces by linear interpolation
     do idx=1,lev_nbr
        alt(idx)=0.5*(alt_ntf(idx)+alt_ntf(idx+1)) ! [m]
        tpt(idx)=0.5*(tpt_ntf(idx)+tpt_ntf(idx+1)) ! [K]
        prs(idx)=0.5*(prs_ntf(idx)+prs_ntf(idx+1)) ! [Pa]
        cnc_air(idx)=0.5*(cnc_air(idx)+cnc_air(idx+1)) ! [# m-3]
        H2O_vmr(idx)=0.5*(H2O_vmr_ntf(idx)+H2O_vmr_ntf(idx+1)) ! [mlc mlc-1]
        O3_vmr(idx)=0.5*(O3_vmr_ntf(idx)+O3_vmr_ntf(idx+1)) ! [mlc mlc-1]
        N2O_vmr(idx)=0.5*(N2O_vmr_ntf(idx)+N2O_vmr_ntf(idx+1)) ! [mlc mlc-1]
        CO_vmr(idx)=0.5*(CO_vmr_ntf(idx)+CO_vmr_ntf(idx+1)) ! [mlc mlc-1]
        CH4_vmr(idx)=0.5*(CH4_vmr_ntf(idx)+CH4_vmr_ntf(idx+1)) ! [mlc mlc-1]
     enddo                  ! end loop over lev
     
     ! Derive required CLM fields from AFGL input
     do idx=1,lev_nbr
        q_H2O(idx)=H2O_vmr(idx)*(mmw_H2O/mmw_dry_air) ! [kg kg-1]
        q_O3(idx)=O3_vmr(idx)*(mmw_O3/mmw_dry_air) ! [kg kg-1]
     enddo                  ! end loop over lev
     
     prf_sng='AFGL ' // prf_sng(1:len(prf_sng)-5)
     call ftn_strnul(prf_sng) ! [sbr] NUL-initialize all characters after LSC
     prs_sfc=prs_ntf(levp_nbr)
     tpt_sfc=tpt_ntf(levp_nbr)
     ! tpt_sfc is temperature at aerodynamic roughness height
     ! tpt_skn is skin (emitting) temperature
     ! It is possible for there to be a discontinuity between tpt_sfc and tpt_skn
     tpt_skn=tpt_sfc
     
     ! Set non-AFGL fields to defaults
     do idx=1,lev_nbr
        cld_frc(idx)=0.0
        mpl_CWP(idx)=0.0
        odxl_obs_aer(idx)=0.0
        odxl_obs_bga(idx)=0.0
     enddo
     
     ! De-allocate dynamic variables
     if (allocated(CH4_vmr)) deallocate(CH4_vmr,stat=rcd)
     if (allocated(CH4_vmr_ntf)) deallocate(CH4_vmr_ntf,stat=rcd)
     if (allocated(CO_vmr)) deallocate(CO_vmr,stat=rcd)
     if (allocated(CO_vmr_ntf)) deallocate(CO_vmr_ntf,stat=rcd)
     if (allocated(H2O_vmr)) deallocate(H2O_vmr,stat=rcd)
     if (allocated(H2O_vmr_ntf)) deallocate(H2O_vmr_ntf,stat=rcd)
     if (allocated(N2O_vmr)) deallocate(N2O_vmr,stat=rcd)
     if (allocated(N2O_vmr_ntf)) deallocate(N2O_vmr_ntf,stat=rcd)
     if (allocated(O3_vmr)) deallocate(O3_vmr,stat=rcd)
     if (allocated(O3_vmr_ntf)) deallocate(O3_vmr_ntf,stat=rcd)
     if (allocated(cnc_air)) deallocate(cnc_air,stat=rcd)
     if (allocated(cnc_air_ntf)) deallocate(cnc_air_ntf,stat=rcd)
     
  else if (CLM_NC_INPUT) then ! Input file is in CLM_NC_INPUT format
     
     ! Assume lev,prs,T,q_H2O are available in input file in SI units
     ! File is already open because we opened it to obtain lev_nbr

     ! Required fields: Get variable IDs and data
     rcd=nf90_wrp_inq_varid(nc_id,'lev',lev_id)
     if (rcd == nf90_noerr) rcd=nf90_wrp(nf90_get_var(nc_id,lev_id,lev),sbr_nm//": gv lev")
     rcd=nf90_wrp_inq_varid(nc_id,'prs',prs_id)
     if (rcd == nf90_noerr) rcd=nf90_wrp(nf90_get_var(nc_id,prs_id,prs),sbr_nm//": gv prs")
     rcd=nf90_wrp_inq_varid(nc_id,'q_H2O',q_H2O_id)
     if (rcd == nf90_noerr) rcd=nf90_wrp(nf90_get_var(nc_id,q_H2O_id,q_H2O),sbr_nm//": gv q_H2O")
     rcd=nf90_wrp_inq_varid(nc_id,'tpt',tpt_id)
     if (rcd == nf90_noerr) rcd=nf90_wrp(nf90_get_var(nc_id,tpt_id,tpt),sbr_nm//": gv tpt")
     
     ! Fill all other standard CLM fields (q_O3, cld_frc, mpl_CWP, odxl_obs_aer, and scalar fields)
     ! with defaults which MAY be overridden by the input file.
     ! Proceed in the order of a standard CLM input file
     do idx=1,lev_nbr
        if (dbg_lvl > 3) then
           write (6,'(a,i4,a,f9.3,a)') 'Interpolating O3 to prs(',idx,') = ',prs(idx)/100.0,' mb'
        endif                     ! end if dbg
        q_O3(idx)=q_O3_ntp(prs(idx))
        cld_frc(idx)=0.0
        mpl_CWP(idx)=0.0
        odxl_obs_aer(idx)=0.0
        odxl_obs_bga(idx)=0.0
     enddo
     
     ! Set derived fields which depend on input
     prs_sfc=prs(lev_nbr)+100.0 ! [Pa]
     tpt_sfc=tpt(lev_nbr)+1.0 ! [K]
     tpt_skn=tpt(lev_nbr)+1.0 ! [K]
     
     ! Override defaults with file data where possible
     
     ! Get global attributes
     rcd=nf90_wrp(nf90_inquire_attribute(nc_id,nf90_global,'prf_sng',attnum=int_foo),sbr_nm//": inquire_att prf_sng")
     if (rcd == nf90_noerr) rcd=nf90_wrp(nf90_get_att(nc_id,nf90_global,'prf_sng',prf_sng),sbr_nm//": get_att prf_sng")
     rcd=nf90_inquire_attribute(nc_id,nf90_global,'prf_snw_sng',attnum=int_foo)
     if (rcd == nf90_noerr) rcd=nf90_wrp(nf90_get_att(nc_id,nf90_global,'prf_snw_sng',prf_snw_sng),sbr_nm//": get_att prf_snw_sng")
     ! Get arrays
     rcd=nf90_wrp_inq_varid(nc_id,'q_O3',q_O3_id,rcd_opt=NF90_ENOTVAR)
     if (rcd == nf90_noerr) rcd=nf90_wrp(nf90_get_var(nc_id,q_O3_id,q_O3),sbr_nm//": gv q_O3")
     rcd=nf90_wrp_inq_varid(nc_id,'q_NO2',q_NO2_id,rcd_opt=NF90_ENOTVAR)
     if (rcd == nf90_noerr) rcd=nf90_wrp(nf90_get_var(nc_id,q_NO2_id,q_NO2),sbr_nm//": gv q_NO2")
     rcd=nf90_wrp_inq_varid(nc_id,'q_OH',q_OH_id,rcd_opt=NF90_ENOTVAR)
     if (rcd == nf90_noerr) rcd=nf90_wrp(nf90_get_var(nc_id,q_OH_id,q_OH),sbr_nm//": gv q_OH")
     rcd=nf90_wrp_inq_varid(nc_id,'cld_frc',cld_frc_id,rcd_opt=NF90_ENOTVAR)
     if (rcd == nf90_noerr) rcd=nf90_wrp(nf90_get_var(nc_id,cld_frc_id,cld_frc),sbr_nm//": gv cld_frc")
     rcd=nf90_wrp_inq_varid(nc_id,'mpl_CWP',mpl_CWP_id,rcd_opt=NF90_ENOTVAR)
     if (rcd == nf90_noerr) rcd=nf90_wrp(nf90_get_var(nc_id,mpl_CWP_id,mpl_CWP),sbr_nm//": gv mpl_CWP")
     rcd=nf90_wrp_inq_varid(nc_id,'odxl_obs_aer',odxl_obs_aer_id,rcd_opt=NF90_ENOTVAR)
     if (rcd == nf90_noerr) rcd=nf90_wrp(nf90_get_var(nc_id,odxl_obs_aer_id,odxl_obs_aer),sbr_nm//": gv odxl_obs_aer")
     rcd=nf90_wrp_inq_varid(nc_id,'odxl_obs_bga',odxl_obs_bga_id,rcd_opt=NF90_ENOTVAR)
     if (rcd == nf90_noerr) rcd=nf90_wrp(nf90_get_var(nc_id,odxl_obs_bga_id,odxl_obs_bga),sbr_nm//": gv odxl_obs_bga")
     if (flg_snw) then
        rcd=nf90_wrp_inq_varid(nc_id,'dpt_dlt_snw',dpt_dlt_snw_id,rcd_opt=NF90_ENOTVAR)
        if (rcd == nf90_noerr) rcd=nf90_wrp(nf90_get_var(nc_id,dpt_dlt_snw_id,dpt_dlt_snw),sbr_nm//": gv dpt_dlt_snw")
        rcd=nf90_wrp_inq_varid(nc_id,'tpt_snw',tpt_snw_id,rcd_opt=NF90_ENOTVAR)
        if (rcd == nf90_noerr) rcd=nf90_wrp(nf90_get_var(nc_id,tpt_snw_id,tpt_snw),sbr_nm//": gv tpt_snw")
        rcd=nf90_wrp_inq_varid(nc_id,'dns_snw',dns_snw_id,rcd_opt=NF90_ENOTVAR)
        if (rcd == nf90_noerr) rcd=nf90_wrp(nf90_get_var(nc_id,dns_snw_id,dns_snw),sbr_nm//": gv dns_snw")
        rcd=nf90_wrp_inq_varid(nc_id,'mmr_mpr_snw',mmr_mpr_snw_id,rcd_opt=NF90_ENOTVAR)
        if (rcd == nf90_noerr) rcd=nf90_wrp(nf90_get_var(nc_id,mmr_mpr_snw_id,mmr_mpr_snw),sbr_nm//": gv mmr_mpr_snw")
        rcd=nf90_wrp_inq_varid(nc_id,'rds_ffc_snw',rds_ffc_snw_id,rcd_opt=NF90_ENOTVAR)
        if (rcd == nf90_noerr) rcd=nf90_wrp(nf90_get_var(nc_id,rds_ffc_snw_id,rds_ffc_snw),sbr_nm//": gv rds_ffc_snw")
        
     endif ! !flg_snw
     
     ! Get scalars
     rcd=nf90_wrp_inq_varid(nc_id,'lcl_yr_day',lcl_yr_day_id,rcd_opt=NF90_ENOTVAR)
     if (rcd == nf90_noerr) rcd=nf90_wrp(nf90_get_var(nc_id,lcl_yr_day_id,lcl_yr_day),sbr_nm//": gv lcl_yr_day")
     rcd=nf90_wrp_inq_varid(nc_id,'lat',lat_id,rcd_opt=NF90_ENOTVAR)
     if (rcd == nf90_noerr) rcd=nf90_wrp(nf90_get_var(nc_id,lat_id,lat),sbr_nm//": gv lat")
     rcd=nf90_wrp_inq_varid(nc_id,'oro',oro_id,rcd_opt=NF90_ENOTVAR)
     if (rcd == nf90_noerr) rcd=nf90_wrp(nf90_get_var(nc_id,oro_id,oro),sbr_nm//": gv oro")
     rcd=nf90_wrp_inq_varid(nc_id,'prs_sfc',prs_sfc_id,rcd_opt=NF90_ENOTVAR)
     if (rcd == nf90_noerr) rcd=nf90_wrp(nf90_get_var(nc_id,prs_sfc_id,prs_sfc),sbr_nm//": gv prs_sfc")
     rcd=nf90_wrp_inq_varid(nc_id,'tpt_sfc',tpt_sfc_id,rcd_opt=NF90_ENOTVAR)
     if (rcd == nf90_noerr) rcd=nf90_wrp(nf90_get_var(nc_id,tpt_sfc_id,tpt_sfc),sbr_nm//": gv tpt_sfc")
     rcd=nf90_wrp_inq_varid(nc_id,'tpt_skn',tpt_skn_id,rcd_opt=NF90_ENOTVAR)
     if (rcd == nf90_noerr) rcd=nf90_wrp(nf90_get_var(nc_id,tpt_skn_id,tpt_skn),sbr_nm//": gv tpt_skn")
     rcd=nf90_wrp_inq_varid(nc_id,'sfc_ems',sfc_ems_id,rcd_opt=NF90_ENOTVAR)
     if (rcd == nf90_noerr) rcd=nf90_wrp(nf90_get_var(nc_id,sfc_ems_id,sfc_ems),sbr_nm//": gv sfc_ems")
     rcd=nf90_wrp_inq_varid(nc_id,'rgh_len',rgh_len_id,rcd_opt=NF90_ENOTVAR)
     if (rcd == nf90_noerr) rcd=nf90_wrp(nf90_get_var(nc_id,rgh_len_id,rgh_len),sbr_nm//": gv rgh_len")
     rcd=nf90_wrp_inq_varid(nc_id,'snow_depth',snow_depth_id,rcd_opt=NF90_ENOTVAR)
     if (rcd == nf90_noerr) rcd=nf90_wrp(nf90_get_var(nc_id,snow_depth_id,snow_depth),sbr_nm//": gv snow_depth")
     rcd=nf90_wrp_inq_varid(nc_id,'alb_sfc_vsb_drc',alb_sfc_vsb_drc_id,rcd_opt=NF90_ENOTVAR)
     if (rcd == nf90_noerr) rcd=nf90_wrp(nf90_get_var(nc_id,alb_sfc_vsb_drc_id,alb_sfc_vsb_drc),sbr_nm//": gv alb_sfc_vsb_drc")
     rcd=nf90_wrp_inq_varid(nc_id,'alb_sfc_vsb_dff',alb_sfc_vsb_dff_id,rcd_opt=NF90_ENOTVAR)
     if (rcd == nf90_noerr) rcd=nf90_wrp(nf90_get_var(nc_id,alb_sfc_vsb_dff_id,alb_sfc_vsb_dff),sbr_nm//": gv alb_sfc_vsb_dff")
     rcd=nf90_wrp_inq_varid(nc_id,'alb_sfc_NIR_drc',alb_sfc_NIR_drc_id,rcd_opt=NF90_ENOTVAR)
     if (rcd == nf90_noerr) rcd=nf90_wrp(nf90_get_var(nc_id,alb_sfc_NIR_drc_id,alb_sfc_NIR_drc),sbr_nm//": gv alb_sfc_NIR_drc")
     rcd=nf90_wrp_inq_varid(nc_id,'alb_sfc_NIR_dff',alb_sfc_NIR_dff_id,rcd_opt=NF90_ENOTVAR)
     if (rcd == nf90_noerr) rcd=nf90_wrp(nf90_get_var(nc_id,alb_sfc_NIR_dff_id,alb_sfc_NIR_dff),sbr_nm//": gv alb_sfc_NIR_dff")
     rcd=nf90_wrp_inq_varid(nc_id,'CO2_vmr_clm',CO2_vmr_clm_id,rcd_opt=NF90_ENOTVAR)
     if (rcd == nf90_noerr) rcd=nf90_wrp(nf90_get_var(nc_id,CO2_vmr_clm_id,CO2_vmr_clm),sbr_nm//": gv CO2_vmr_clm")
     rcd=nf90_wrp_inq_varid(nc_id,'N2O_vmr_clm',N2O_vmr_clm_id,rcd_opt=NF90_ENOTVAR)
     if (rcd == nf90_noerr) rcd=nf90_wrp(nf90_get_var(nc_id,N2O_vmr_clm_id,N2O_vmr_clm),sbr_nm//": gv N2O_vmr_clm")
     rcd=nf90_wrp_inq_varid(nc_id,'CH4_vmr_clm',CH4_vmr_clm_id,rcd_opt=NF90_ENOTVAR)
     if (rcd == nf90_noerr) rcd=nf90_wrp(nf90_get_var(nc_id,CH4_vmr_clm_id,CH4_vmr_clm),sbr_nm//": gv CH4_vmr_clm")
     rcd=nf90_wrp_inq_varid(nc_id,'CFC11_vmr_clm',CFC11_vmr_clm_id,rcd_opt=NF90_ENOTVAR)
     if (rcd == nf90_noerr) rcd=nf90_wrp(nf90_get_var(nc_id,CFC11_vmr_clm_id,CFC11_vmr_clm),sbr_nm//": gv CFC11_vmr_clm")
     rcd=nf90_wrp_inq_varid(nc_id,'CFC12_vmr_clm',CFC12_vmr_clm_id,rcd_opt=NF90_ENOTVAR)
     if (rcd == nf90_noerr) rcd=nf90_wrp(nf90_get_var(nc_id,CFC12_vmr_clm_id,CFC12_vmr_clm),sbr_nm//": gv CFC12_vmr_clm")
     rcd=nf90_wrp_inq_varid(nc_id,'odxc_obs_aer',odxc_obs_aer_id,rcd_opt=NF90_ENOTVAR)
     if (rcd == nf90_noerr) rcd=nf90_wrp(nf90_get_var(nc_id,odxc_obs_aer_id,odxc_obs_aer),sbr_nm//": gv odxc_obs_aer")
     rcd=nf90_wrp_inq_varid(nc_id,'slr_cst',slr_cst_id,rcd_opt=NF90_ENOTVAR)
     if (rcd == nf90_noerr) rcd=nf90_wrp(nf90_get_var(nc_id,slr_cst_id,slr_cst),sbr_nm//": gv slr_cst")
     rcd=nf90_wrp_inq_varid(nc_id,'gmt_yr',gmt_yr_id,rcd_opt=NF90_ENOTVAR)
     if (rcd == nf90_noerr) rcd=nf90_wrp(nf90_get_var(nc_id,gmt_yr_id,gmt_yr),sbr_nm//": gv gmt_yr")
     rcd=nf90_wrp_inq_varid(nc_id,'lon_dgr',lon_dgr_id,rcd_opt=NF90_ENOTVAR)
     if (rcd == nf90_noerr) rcd=nf90_wrp(nf90_get_var(nc_id,lon_dgr_id,lon_dgr),sbr_nm//": gv lon_dgr")
     ! Wrap
     rcd=nf90_wrp_inq_varid(nc_id,'frc_str_zen_ngl_sfc',frc_str_zen_ngl_sfc_id,rcd_opt=NF90_ENOTVAR)
     if (rcd == nf90_noerr) rcd=nf90_wrp(nf90_get_var(nc_id,frc_str_zen_ngl_sfc_id,frc_str_zen_ngl_sfc), &
          sbr_nm//": gv frc_str_zen_ngl_sfc")
     ! Close file
     rcd=nf90_wrp_close(nc_id,fl_in,'Ingested') ! [fnc] Close file
  endif                     ! endif input file is in netCDF format

  ! All necessary input data have been read
  ! Data not stored or input as SI have been converted to SI
  ! Initialize any input data that may overridden by command line switches
  
  ! Determine solar geometry
  call slr_crd_Bri92( &
       lat,                 & ! I [rdn] Latitude
       lcl_yr_day,          & ! I [day] Local year day
       slr_zen_ngl_cos,     & ! O [frc] Solar zenith angle cosine
       xnt_fac)             ! O [frc] Eccentricity factor
  
  lcl_time_hr=(lcl_yr_day-int(lcl_yr_day))*24.0 ! [hr] Local time hour
  slr_dst_au=1.0/sqrt(xnt_fac)
  slr_zen_ngl=acos(slr_zen_ngl_cos)
  slr_zen_ngl_dgr=180.0*slr_zen_ngl/pi
  slr_cst=slr_cst_CCM
  slr_flx_nrm_TOA=slr_cst*xnt_fac
  slr_flx_TOA=slr_cst*xnt_fac*slr_zen_ngl_cos
  
  ! Determine time coordinates necessary for albedo and aerosol calculations
  ! Time/date variables correspond to output of UNIX gmtime() function
  
  ! There are three typical scenarios we encounter in computing time coordinates:
  ! 1. Given LTST, lon, and year, compute GMT
  ! 2. Given GMT or UNIX time and lon, compute LTST
  ! 3. Given slr_zen_ngl_cos 
  ! Scenario 1 corresponds to typical CRM input files where ltst_doy is specified
  ! Scenario 1 requires longitude and year be input before proceeding
  ! Scenario 2 occurs for high accuracy forward RT computations like ARESE
  ! In scenario 2, the exact geographic location and GMT to be simulated are known
  ! Scenario 3 is used to run RT codes with known solar zenith angles
  ! In scenario 3, we do not care about the exact time or location, and attempting to guess them is ill-constrained 
  ! What should be stored in the time variables in scenario 3?
  
  ! The following implements scenario 1, TST->GMT:
  ltst_doy=lcl_yr_day
  do while(lon_dgr < 0.0)
     lon_dgr=lon_dgr+360.0  ! [dgr]
  end do                    ! end while
  do while(lon_dgr > 360.0)
     lon_dgr=lon_dgr-360.0  ! [dgr]
  end do                    ! end while
  lon=pi*lon_dgr/180.0      ! [dgr] -> [rdn]
  
  ! Given ltst_doy, lon, and gmt_yr, we can derive all other times 
  ! Sign convention for longitude is extremely important here
  ! Do not change this without looking carefully at implications
  if (lon_dgr >= 0.0.and.lon_dgr <= 180.0) then
     lon_dgr_m180_p180=lon_dgr ! [dgr]
  else if (lon_dgr >= 180.0.and.lon_dgr <= 360.0) then
     lon_dgr_m180_p180=lon_dgr-360.0 ! [dgr]
  else
     write (6,'(a,a,f9.3)') prg_nm(1:ftn_strlen(prg_nm)),': ERROR lon_dgr = ',lon_dgr
     stop
  endif                     ! endif
  gtst_doy=ltst_doy-lon_dgr_m180_p180/360.0
  
  ! Estimate equation of time correction by assuming gtst = gmt
  ! Estimated error of this approximation follows:
  ! eqn_time_mnt is <~ 20 minutes all year
  ! This propogates relative error of < +/- 4.e-5 days to doy used to compute slr_crd_gmm
  ! This results in absolute error of < +/- 6.e-7 radians in slr_crd_gmm
  ! This results in relative error of < +/- 1.0e-4 in eqn_tim_mnt
  ! This results in an absolute error of < +/- 0.1 seconds in gmt_doy
  slr_crd_gmm=2.0*pi*(gtst_doy-1.0)/365.0
  slr_crd_gmm_dgr=180.0*slr_crd_gmm/pi
  eqn_time_mnt=229.18*( &
       0.000075+0.001868*cos(slr_crd_gmm)- &
       0.032077*sin(slr_crd_gmm)- &
       0.014615*cos(2.0*slr_crd_gmm)- &
       0.040849*sin(2.0*slr_crd_gmm) &
       )
  eqn_time_day=eqn_time_mnt/1440.0
  ! Apply equation of time correction to get gmt from gtst
  gmt_doy=gtst_doy-eqn_time_day
  ! Recomputing slr_crd_gmm and eqn_time_mnt correction based on new gmt would introduce 
  ! slight inconsistency between ltst_doy and gmt, so do not iterate the above.
  eqn_time_sec=eqn_time_mnt*60.0 ! Equation of time in seconds of time
  lon_sec=86400.0*lon_dgr_m180_p180/360.0 ! lmt minus gmt in seconds of time
  
  ! Get UNIX time offset of beginning of year
  int_foo=0
  int_foo2=1 ! Recall gmt2unix uses 0-based months and days of month
  call gmt2unix(gmt_yr,int_foo2,int_foo2,int_foo,int_foo,int_foo,time_unix_yr_srt)
  
  ! Add seconds since beginning of year to get complete UNIX time representation
  time_unix=time_unix_yr_srt+(gmt_doy-1.0)*86400.0 ! [s] Seconds between 1969 and GMT of simulation
  time_lmt=time_unix+lon_sec ! [s] Seconds between 1969 and LMT of simulation
  time_ltst=time_lmt+eqn_time_sec ! [s] Seconds between 1969 and LTST of simulation
  
  ! Get parsed date information and string for GMT, LMT and LTST
  ! Call to unix2gmt() gets numerical values which are archived in netCDF
  ! Call to unix2gmt_sng() gets string, e.g., Fri Sep 16 13:11:08 1994 GMT
  ! Second call 
  ! GMT
  tz_sng='GMT'
  call unix2gmt(time_unix,gmt_yr,gmt_ydy,gmt_mth,gmt_day,gmt_hr,gmt_mnt,gmt_sec)
  call unix2gmt_sng(time_unix,gmt_sng,tz_sng)
  
  ! LMT
  tz_sng='LMT'
  call unix2gmt(time_lmt,lmt_yr,lmt_ydy,lmt_mth,lmt_day,lmt_hr,lmt_mnt,lmt_sec)
  call unix2gmt_sng(time_lmt,lmt_sng,tz_sng)
  
  ! LTST
  tz_sng='TST'
  call unix2gmt(time_ltst,ltst_yr,ltst_ydy,ltst_mth,ltst_day,ltst_hr,ltst_mnt,ltst_sec)
  call unix2gmt_sng(time_ltst,ltst_sng,tz_sng)
  
  if (dbg_lvl == 3) then
     write (0,*) 'Input ltst_doy  = ',ltst_doy
     write (0,*) 'Input gmt_doy  = ',gmt_doy
     write (0,*) 'Output gmt_doy = ',gmt_ydy+dble(gmt_hr/24.0+gmt_mnt/1440.0+gmt_sec/86400.0)
     write (0,*) 'time_unix_yr_srt = ',time_unix_yr_srt
     write (0,'(1x,a,f16.6)') 'time_unix = ',time_unix
     write (0,'(1x,a,f16.6)') 'time_lmt = ',time_lmt
     write (0,'(1x,a,f16.6)') 'time_ltst = ',time_ltst
     write (0,*) 'gmt_yr = ',gmt_yr
     write (0,*) 'gmt_ydy = ',gmt_ydy
     write (0,*) 'gmt_mth = ',gmt_mth
     write (0,*) 'gmt_day = ',gmt_day
     write (0,*) 'gmt_hr = ',gmt_hr
     write (0,*) 'gmt_mnt = ',gmt_mnt
     write (0,*) 'gmt_sec = ',gmt_sec
     write (0,*) 'eqn_time_sec  = ',eqn_time_sec
     write (0,*) 'lon_sec  = ',lon_sec
     write (0,*) 'gmt_sng  = ',gmt_sng
     write (0,*) 'lmt_sng  = ',lmt_sng
     write (0,*) 'ltst_sng = ',ltst_sng
     write (0,*) 'Bri92 slr_zen_ngl_dgr  = ',slr_zen_ngl_dgr
     write (0,*) 'Bri92 slr_dst_au       = ',slr_dst_au
     write (0,*) 'Bri92 slr_flx_nrm_TOA  = ',slr_flx_nrm_TOA
     write (0,*) 'Bri92 slr_flx_TOA      = ',slr_flx_TOA
  endif                     ! endif dbg
  
  ! The following implements scenario 2, GMT->TST:
  ! Currently, this overrides Bri92 slr_crd fields with more accurate Mic88 fields
  ! Doing this, however, requires we know GMT
  ! Thus Bri92 is used for TST->GMT, and Mic88 is used for GMT->slr_crd
  
  ! Offline calculations show Bri92 (scenario 1) and Mic88 (scenario 2) predict different GMT for LTST noon
  ! For LTST noon, Oct 11, 1995 at the SGP CART site:
  ! 1. LTST -> GMT method with Bri92 yields GMT 18:16:17
  ! 2. GMT -> LTST method with Mic88 yields GMT 18:16:43
  ! The difference is 26 s, suspiciously close to the number of leap seconds from 1970--1995.
  ! The 26 s causes a TOA flux difference of 0.002 W m-2
  ! This small flux bias is much smaller than the bias in TOA fluxes between Bri92 and Mic88 for the same GMT.
  ! For Oct. 11 1995 noon TST, the two routines differ in sun-earth distance by 0.0005 AU ( = 0.05%), 
  ! and in solar zenith angle by 0.25 degrees, and in TOA flux by 3.25 W m-2.
  ! Thus ~99% of the 3.25 W m-2 difference is due to Bri92 errors in solar geometry, not in specification of time.
  ! The main Bri92 error seems to be due to neglecting the year (i.e., 1995) in the geometry calculation.
  
  if (.true.) then
     
     lat_dgr=180.0*lat/pi
     gmt_hr_dcm=24.0*(gmt_doy-gmt_ydy)
     call slr_crd_Mic88( & ! Input
          real(gmt_hr_dcm), &
          real(gmt_ydy), &
          real(gmt_yr), &
          real(lat_dgr), &
          real(lon_dgr_m180_p180), &
          ! Output
          eqn_time_sec,     & ! Equation of time (second)
          slr_azi_dgr,      & ! Solar azimuth (degree)
          slr_dcl_dgr,      & ! Solar declination (degree)
          slr_dmt_dgr,      & ! Diameter of solar disc (degree)
          slr_dst_au,       & ! Solar distance (AU)
          slr_elv_dgr,      & ! Solar elevation (degree)
          slr_hr_ngl_dgr,   & ! Solar hour angle (degree)
          slr_rgt_asc_dgr,  & ! Solar right ascension (degree)
          slr_rfr_ngl_dgr   & ! Solar refraction angle (degree)
          )
     
     slr_zen_ngl_dgr=90.0-slr_elv_dgr
     slr_zen_ngl=pi*slr_zen_ngl_dgr/180.0
     slr_zen_ngl_cos=cos(slr_zen_ngl)
     xnt_fac=1.0/(slr_dst_au*slr_dst_au)
     slr_cst=slr_cst_CCM
     slr_flx_nrm_TOA=slr_cst*xnt_fac
     slr_flx_TOA=slr_cst*xnt_fac*slr_zen_ngl_cos
     
     if (dbg_lvl == 3) then
        write (0,*) 'lat_dgr                = ',lat_dgr
        write (0,*) 'lon_dgr                = ',lon_dgr
        write (0,*) 'Mic88 slr_zen_ngl_dgr  = ',slr_zen_ngl_dgr
        write (0,*) 'Mic88 slr_dst_au       = ',slr_dst_au
        write (0,*) 'Mic88 slr_flx_nrm_TOA  = ',slr_flx_nrm_TOA
        write (0,*) 'Mic88 slr_flx_TOA      = ',slr_flx_TOA
        write (0,*) 'Mic88 eqn_time_sec          = ',eqn_time_sec
     endif                  ! endif dbg
     
  endif                     ! endif Mic88
  
  ! Do background aerosol first, because it may be subtracted from observed aerosol
  if (cmd_ln_odxc_obs_bga.or.flg_bga) then
     ! Initialize background aerosol parameters
     prs_top_bga=10.0e2     ! [Pa] top level of background aerosol
     prs_btm_bga=100.0e2    ! [Pa] bottom level of background aerosol
     
     ! Background aerosol differs from regular aerosol in that, if cmd_ln_odxc_obs_bga is set, 
     ! background aerosol will be placed in all appropriate levels up to the specified
     ! optical depth, whereas regular tropospheric aerosol is only placed in levels
     ! specified in the input file.
     ! Count levels containing background aerosol 
     bga_lvl_nbr=0
     do idx=1,lev_nbr
        if (prs(idx) > prs_top_bga.and.prs(idx) < prs_btm_bga) then
           bga_lvl_nbr=bga_lvl_nbr+1
           bga_lvl_flg(idx)=.true.
        else
           bga_lvl_flg(idx)=.false.
        endif
     enddo                  ! end loop over lev
     ! Apportion background aerosol equally amongst these levels
     if (bga_lvl_nbr <= 0) stop 'bga_lvl_nbr == 0'
     
     if (cmd_ln_odxc_obs_bga) odxc_obs_bga=odxc_obs_bga_cmd_ln
     do idx=1,lev_nbr
        if (bga_lvl_flg(idx)) then
           odxl_obs_bga(idx)=odxc_obs_bga/bga_lvl_nbr
        else
           odxl_obs_bga(idx)=0.0
        endif
     enddo                  ! end loop over lev
  endif                     ! endif cmd_ln_odxc_obs_bga
  
  ! Use GMT and wvl_obs to interpolate odxc_obs_aer from MFRSR data
  if (flg_MFRSR_obs) then
     call aer_odxc_get(wvl_obs_aer,gmt_doy,odxc_obs_aer)
     write(0,'(a,f6.4)') 'MFRSR Observed total aerosol = ',odxc_obs_aer
  endif                     ! endif MFRSR
  
  ! Aerosol optical depth is specified on command line or in external file
  if (cmd_ln_odxc_obs_aer.or.flg_MFRSR_obs) then
     aer_lvl_nbr=0
     do idx=1,lev_nbr
        if (odxl_obs_aer(idx) > 0.0) aer_lvl_nbr=aer_lvl_nbr+1
     enddo                  ! end loop over lev
     if (aer_lvl_nbr == 0) stop 'ERROR: aer_lvl_nbr == 0 but aerosol is expected'
     if (cmd_ln_odxc_obs_aer) odxc_obs_aer=odxc_obs_aer_cmd_ln
     ! Subtract observed background aerosol from observed total
     if (cmd_ln_odxc_obs_bga.or.flg_bga) then
        ! Sanity check
        if (wvl_obs_aer /= wvl_obs_bga) stop 'Attempting to subtract optical depths when wvl_obs_aer /= wvl_obs_bga'
        write(0,'(a,f6.4,a3,f6.4)') 'Subtracting aerosols: Tropospheric aerosol = ',odxc_obs_aer,' - ',odxc_obs_bga
        odxc_obs_aer=odxc_obs_aer-odxc_obs_bga
     endif                  ! endif subtracting background aerosol from total
     do idx=1,lev_nbr
        if (odxl_obs_aer(idx) > 0.0) odxl_obs_aer(idx)=odxc_obs_aer/aer_lvl_nbr
     enddo                  ! end loop over lev
  endif                     ! end if overriding CLM profile odxl_obs_aer
  
  ! Get aerosol-specific information
  call aer_info_get(fl_aer,wvl_obs_aer,dns_aer,ext_cff_mss_aer)
  call aer_info_get(fl_bga,wvl_obs_bga,dns_bga,ext_cff_mss_bga)
  ! Convert optical depths to mass paths by using the extinction efficiency per unit mass
  do idx=1,lev_nbr
     mpl_aer(idx)=odxl_obs_aer(idx)/ext_cff_mss_aer ! [kg m-2]
     mpl_bga(idx)=odxl_obs_bga(idx)/ext_cff_mss_bga ! [kg m-2]
  enddo                     ! end loop over lev
  
  ! Tune water vapor and ozone column amounts 
  ! Ozone can be simply tuned to match observations using this method
  ! Water vapor, however, subsequently gets adjusted to saturation in 
  ! cloudy layers (to prevent supersaturation).
  ! Thus, matching H2O exactly can be an iterative process
  if (H2O_tune_factor /= 1.0) then
     do idx=1,lev_nbr
        q_H2O(idx)=q_H2O(idx)*H2O_tune_factor ! [kg kg-1]
     enddo
  endif
  if (O3_tune_factor /= 1.0) then
     do idx=1,lev_nbr
        q_O3(idx)=q_O3(idx)*O3_tune_factor ! [kg kg-1]
     enddo
  endif
  
  ! Compute ice fraction
  do idx=1,lev_nbr
     ! Define fractional amount of cloud that is ice per strategy of CCM3:
     ! Warmer than -10 degrees C -> liquid phase
     ! Colder than -10 degrees C but warmer than -40 C -> mixed phase
     ! Colder than -40 degrees C -> ice phase
     if (tpt(idx) > 263.16) frc_ice(idx)=0.0
     if ((tpt(idx) <= 263.16).and.(tpt(idx) >= 233.16)) then
        frc_ice(idx)=(263.16-tpt(idx))/30.0
     endif
     if (tpt(idx) < 233.16) frc_ice(idx)=1.0
  enddo                     ! end loop over lev
  
  ! Compute saturation properties of profile: 
  ! If force_cld_lvl_2b_sat (force cloudy levels to be saturated) switch is on, 
  ! then alter moisture profile so all levels that contain suspended condensate are 
  ! saturated relative to prevailing phases of water.
  ! If force_sat_lvl_2b_cld (force saturated levels to be cloudy) switch is on, 
  ! then alter cloud profile so all saturated levels contain condensate such that
  ! total condensate in profile matches input CWP.
  do idx=1,lev_nbr
     tpt_cls(idx)=tpt(idx)-tpt_frz_pnt ! Bol80 p. 1047
     
     ! Compute saturation with respect to liquid 
     ! sat_vpr_lqd=svp_H2O_lqd_Bol80(tpt(idx))
     sat_vpr_lqd=svp_H2O_lqd_PrK78(tpt(idx))
     ! sat_vpr_lqd=100.0*exp(21.6-5420.0/tpt(idx)) ! [Pa] Salby
     ! NB: Do not use svp_H2O_lqd_Wex() until it is fixed
     ! sat_vpr_lqd=svp_H2O_lqd_Wex(tpt(idx)) ! [Pa] Bol80 (Wexler)
     if (prs(idx) < one_mns_eps_H2O*sat_vpr_lqd) then
        write (0,'(2a,2(a,f9.3,a),2(a,i4,a,f9.3,a),a)')  &
             prg_nm(1:ftn_strlen(prg_nm)),': WARNING ', &
             'sat_vpr_lqd = ',sat_vpr_lqd/100.0,' mb and ', &
             'one_mns_eps_H2O*sat_vpr_lqd = ',one_mns_eps_H2O*sat_vpr_lqd/100.0,' mb at ', &
             'prs(',idx,') = ',prs(idx)/100.0,' mb, and ', &
             'tpt_cls(',idx,') = ',tpt_cls(idx),' C. ', &
             'Setting qst_H2O_lqd = RH_lqd = 0.0'
        qst_H2O_lqd(idx)=0.0
     else
        qst_H2O_lqd(idx)=eps_H2O*sat_vpr_lqd/(prs(idx)-one_mns_eps_H2O*sat_vpr_lqd)
     endif                  ! endif dbg
     
     ! Compute saturation with respect to ice
     ! double_foo=100.0*exp(24.3-6148./tpt(idx)) ! [mb] --> [Pa] Salby
     ! sat_vpr_ice=real(double_foo)
     ! NB: svp_H2O_ice_GoG() has not been thoroughly tested
     ! sat_vpr_ice=svp_H2O_ice_GoG(tpt(idx))
     sat_vpr_ice=svp_H2O_ice_PrK78(tpt(idx))
     if (prs(idx) < one_mns_eps_H2O*sat_vpr_ice) then
        write (0,'(2a,2(a,f9.3,a),2(a,i4,a,f9.3,a),a)')  &
             prg_nm(1:ftn_strlen(prg_nm)),': WARNING ', &
             'sat_vpr_ice = ',sat_vpr_ice/100.0,' mb and ', &
             'one_mns_eps_H2O*sat_vpr_ice = ',one_mns_eps_H2O*sat_vpr_ice/100.0,' mb at ', &
             'prs(',idx,') = ',prs(idx)/100.0,' mb, and ', &
             'tpt_cls(',idx,') = ',tpt_cls(idx),' C. ', &
             'Setting qst_H2O_ice = RH_ice = 0.0'
        qst_H2O_ice(idx)=0.0
     else
        qst_H2O_ice(idx)=eps_H2O*sat_vpr_ice/(prs(idx)-one_mns_eps_H2O*sat_vpr_ice)
     endif                  ! endif dbg
     
     ! Sanity check
     if (sat_vpr_lqd <= 0.0) stop 'sat_vpr_lqd <= 0.0 in main()'
     if (sat_vpr_ice <= 0.0) stop 'sat_vpr_ice <= 0.0 in main()'
     if (qst_H2O_lqd(idx) < 0.0) stop 'qst_H2O_lqd(idx) < 0.0 in main()'
     if (qst_H2O_ice(idx) < 0.0) stop 'qst_H2O_ice(idx) < 0.0 in main()'
     if (force_cld_lvl_2b_sat) then
        if (mpl_CWP(idx) > 0.0) then
           ! NB: Whether layers containing ice clouds should be only ice-saturated
           ! is debatable since air near ice clouds is often ice-supersaturated
           q_H2O(idx)= &
                ! (1-frc_ice(idx))*qst_H2O_lqd(idx)+ &
                !    frc_ice(idx)*qst_H2O_ice(idx) 
                qst_H2O_lqd(idx)
        endif
     endif
  enddo
  
  if (cmd_ln_mpc_CWP.and..not.force_sat_lvl_2b_cld) then
     cld_lvl_nbr=0
     do idx=1,lev_nbr
        if (mpl_CWP(idx) > 0.0) cld_lvl_nbr=cld_lvl_nbr+1
     enddo                  ! end loop over lev
     if (cld_lvl_nbr == 0) stop 'cld_lvl_nbr == 0'
     do idx=1,lev_nbr
        if (mpl_CWP(idx) > 0.0) mpl_CWP(idx)=mpc_CWP_cmd_ln/cld_lvl_nbr
     enddo                  ! end loop over lev
  endif                     ! end if overriding CLM profile mpl_CWP
  
  ! When processing sonde data, it may be more useful to constrain the saturated 
  ! layers to comprise a given CWP (e.g., provided by radar).
  if (force_sat_lvl_2b_cld) then
     cld_lvl_nbr=0
     do idx=1,lev_nbr
        if (q_H2O(idx) > qst_H2O_lqd(idx)) cld_lvl_nbr=cld_lvl_nbr+1
     enddo                  ! end loop over lev
     do idx=1,lev_nbr
        if (q_H2O(idx) > qst_H2O_lqd(idx)) then
           if (cmd_ln_mpc_CWP) mpl_CWP(idx)=mpc_CWP_cmd_ln/cld_lvl_nbr 
           cld_frc(idx)=1.0
        else
           cld_frc(idx)=0.0
        endif
     enddo                  ! end loop over lev
  endif
  
  ! Once layer constituents (CWP, aer) have been processed, compute column statistics
  mpc_CWP=0.0
  mpc_IWP=0.0
  mpc_LWP=0.0
  odxc_obs_aer=0.0
  mpc_aer=0.0
  mpc_bga=0.0
  do idx=1,lev_nbr
     mpl_IWP(idx)=frc_ice(idx)*mpl_CWP(idx)
     mpl_LWP(idx)=max(0.0,mpl_CWP(idx)-mpl_IWP(idx))
     mpc_CWP=mpc_CWP+mpl_CWP(idx)
     mpc_IWP=mpc_IWP+mpl_IWP(idx)
     mpc_LWP=mpc_LWP+mpl_LWP(idx)
     odxc_obs_aer=odxc_obs_aer+odxl_obs_aer(idx)
     mpc_aer=mpc_aer+mpl_aer(idx)
     mpc_bga=mpc_bga+mpl_bga(idx)
  enddo                     ! end loop over lev
  if (mpc_CWP /= 0.0) then
     frc_ice_ttl=mpc_IWP/mpc_CWP
  else
     frc_ice_ttl=0.0
  endif
  
  ! Compute diagnostic variables  
  lat_dgr=180.0*lat/pi      ! [dgr] Latitude
  lat_cos=cos(lat)          ! [frc] Latitude cosine
  if (prs_top <= 0.0) prs_top=prs(1)/2.0 ! [Pa] Pressure at top of "atmosphere"
  prs_ntf(1)=prs_top ! [Pa] Pressure at top of "atmosphere"
  prs_ntf(levp_nbr)=prs_sfc
  do idx=2,lev_nbr
     prs_ntf(idx)=0.5*(prs(idx-1)+prs(idx))
  enddo
  do idx=1,lev_nbr
     lev(idx)=prs(idx)
     prs_dlt(idx)=prs_ntf(idx+1)-prs_ntf(idx)
     tpt_vrt(idx)=(1.0+eps_H2O_rcp_m1*q_H2O(idx))*tpt(idx)
     q_CO2(idx)=CO2_vmr_clm*(mmw_CO2/mmw_dry_air)
     q_CH4(idx)=CH4_vmr_clm*(mmw_CH4/mmw_dry_air)
     q_N2O(idx)=N2O_vmr_clm*(mmw_N2O/mmw_dry_air)
     q_CFC11(idx)=CFC11_vmr_clm*(mmw_CFC11/mmw_dry_air)
     q_CFC12(idx)=CFC12_vmr_clm*(mmw_CFC12/mmw_dry_air)
     q_NO2(idx)=q_NO2_ntp(prs(idx))
     q_OH(idx)=q_OH_ntp(prs(idx))
  enddo
  do idx=1,levp_nbr
     levp(idx)=prs_ntf(idx)
  enddo
  if (flg_snw) then
     dpt_ntf_snw(1)=0.0
     tpt_ntf_snw(1)=tpt_skn
     do levp_snw_idx=2,levp_snw_nbr
        lev_snw_idx=levp_snw_idx-1
        dpt_ntf_snw(levp_snw_idx)=dpt_ntf_snw(levp_snw_idx-1)+dpt_dlt_snw(lev_snw_idx) ! [m] Snow depth interfaces
        levp_snw(levp_snw_idx)=dpt_ntf_snw(levp_snw_idx) ! [m] Snow depth interfaces
     enddo ! end loop over levp_snw
     do lev_snw_idx=1,lev_snw_nbr
        dpt_snw(lev_snw_idx)=0.5*(dpt_ntf_snw(lev_snw_idx)+dpt_ntf_snw(lev_snw_idx+1)) ! [m] Snow depth
        lev_snw(lev_snw_idx)=dpt_snw(lev_snw_idx) ! [m] Snow depth
     enddo ! end loop over snow layer
     do levp_snw_idx=2,levp_snw_nbr-1
        lev_snw_idx=levp_snw_idx-1
        tpt_ntf_snw(levp_snw_idx)=0.5*(tpt_snw(lev_snw_idx)+tpt_snw(lev_snw_idx+1)) ! [K] Snow temperature interfaces
     enddo ! end loop over levp_snw
     tpt_ntf_snw(levp_snw_nbr)=tpt_snw(lev_snw_nbr)+ &
          (dpt_ntf_snw(levp_snw_nbr)-dpt_snw(lev_snw_nbr))* &
          (tpt_snw(lev_snw_nbr)-tpt_ntf_snw(levp_snw_nbr-1))/ &
          (dpt_snw(lev_snw_nbr)-dpt_ntf_snw(levp_snw_nbr-1))
     if (dbg_lvl == dbg_old) then
        do lev_snw_idx=1,lev_snw_nbr
           write (6,'(a8,i2,a4,es8.1)') 'lev_snw(',lev_snw_idx,') = ',lev_snw(lev_snw_idx)
        enddo ! !flg_snw
     endif ! end if dbg
  endif ! !flg_snw
  
  ! CCM: physics/radtpl() sets tpt_ntf(1)=t(1). This approximation can be improved
  ! with a hypsometric integration, or an assumption about the lapse rate.
  
  tpt_ntf(1)=tpt(1)
  tpt_ntf(levp_nbr)=tpt_skn ! not tpt_sfc
  do idx=2,lev_nbr
     tpt_ntf(idx)=tpt(idx)- &
          (tpt(idx)-tpt(idx-1))* &
          (log(prs_ntf(idx))-log(prs(idx)))/ &
          (log(prs(idx-1))-log(prs(idx)))
  enddo
  do idx=1,levp_nbr
     tpt_cls_ntf(idx)=tpt_ntf(idx)-tpt_frz_pnt
  enddo
  
  ! Use surface gravity and midpoint temperature, gas constant
  ! to bootstrap scale height and find first lowest z, alt_ntf.
  alt_ntf(levp_nbr)=0.0
  grv(lev_nbr)=grv_mean_sfc
  scl_hgt(lev_nbr)= &
       gas_cst_dry_air*tpt_vrt(lev_nbr)/grv_mean_sfc
  alt_ntf(lev_nbr)=alt_ntf(levp_nbr)+ &
       scl_hgt(lev_nbr)*log(prs_ntf(levp_nbr)/prs_ntf(lev_nbr))
  alt(lev_nbr)=scl_hgt(lev_nbr)*log(prs_ntf(levp_nbr)/prs(lev_nbr))
  
  ! Compute rest of interface heights and layer heights
  ! that gravity at base of layer gives a valid scale height,
  ! and then refining that gravity to be a mid-layer quantity once
  ! we have mid-layer height.
  do idx=lev_nbr-1,1,-1
     ! Compute gravity at lower interface of layer
     grv(idx)=grv_mean_sfc/ &
          (1.0+alt_ntf(idx+1)/rds_earth)**2
     
     ! Compute scale height with mid-layer tpt_vrt and lower interface gravity
     scl_hgt(idx)= &
          gas_cst_dry_air*tpt_vrt(idx)/grv(idx)
     alt_ntf(idx)=alt_ntf(idx+1)+ &
          scl_hgt(idx)*log(prs_ntf(idx+1)/prs_ntf(idx))
     alt(idx)=alt_ntf(idx+1)+ &
          scl_hgt(idx)*log(prs_ntf(idx+1)/prs(idx))
  enddo
  do idx=1,lev_nbr
     ! Recompute and store gravity and scale height at mid-layers
     grv(idx)= &
          grv_mean_sfc/(1.0+alt(idx)/rds_earth)**2
     scl_hgt(idx)= &
          gas_cst_dry_air*tpt_vrt(idx)/grv(idx)
  enddo
  do idx=1,lev_nbr
     alt_dlt(idx)=alt_ntf(idx)-alt_ntf(idx+1)
  enddo                     ! end loop over lev
  
  ! Compute thermodynamic quantities
  do idx=1,lev_nbr
     mmw_mst_air(idx)= &
          1.0/((1.0-q_H2O(idx))/mmw_dry_air+q_H2O(idx)/mmw_H2O)
     gas_cst_mst_air(idx)= &
          gas_cst_unv/mmw_mst_air(idx)
     dns_mst_air(idx)= &
          prs(idx)/ &
          (gas_cst_dry_air*tpt_vrt(idx))
     if (qst_H2O_ice(idx) > 0.0) then 
        RH_ice(idx)=q_H2O(idx)/qst_H2O_ice(idx)
     else
        RH_ice(idx)=0.0
     endif                  ! endif
     if (qst_H2O_lqd(idx) > 0.0) then 
        RH_lqd(idx)=q_H2O(idx)/qst_H2O_lqd(idx)
     else
        RH_lqd(idx)=0.0
     endif                  ! endif
     if (tpt(idx) > tpt_frz_pnt) then 
        RH(idx)=RH_lqd(idx) 
     else 
        RH(idx)=RH_ice(idx)
     endif
     
     ! Only trust r_H2O until verifying this relationship holds for non-vapor constituents
     r_H2O(idx)=q_H2O(idx)/(1.0-q_H2O(idx)) ! Mixing ratio
     r_CO2(idx)=q_CO2(idx)/(1.0-q_CO2(idx))
     r_CH4(idx)=q_CH4(idx)/(1.0-q_CH4(idx))
     r_N2O(idx)=q_N2O(idx)/(1.0-q_N2O(idx))
     r_CFC11(idx)=q_CFC11(idx)/(1.0-q_CFC11(idx))
     r_CFC12(idx)=q_CFC12(idx)/(1.0-q_CFC12(idx))
     r_O3(idx)=q_O3(idx)/(1.0-q_O3(idx))
     r_O2(idx)=vmr_std_O2*mmw_O2/mmw_dry_air
     r_N2(idx)=vmr_std_N2*mmw_N2/mmw_dry_air
     r_NO2(idx)=q_NO2(idx)/(1.0-q_NO2(idx))
     r_OH(idx)=q_OH(idx)/(1.0-q_OH(idx))
     spc_heat_mst_air(idx)= &
          spc_heat_dry_air* &
          (1.0+cp_vpr_rcp_cp_dry_m1*q_H2O(idx)) ! (IrG81 pp. 77)
     mpl_mst_air(idx)=prs_dlt(idx)/grv(idx)
     mpl_H2O(idx)=q_H2O(idx)*mpl_mst_air(idx)
     mpl_CO2(idx)=q_CO2(idx)*mpl_mst_air(idx)
     mpl_CH4(idx)=q_CH4(idx)*mpl_mst_air(idx)
     mpl_N2O(idx)=q_N2O(idx)*mpl_mst_air(idx)
     mpl_CFC11(idx)=q_CFC11(idx)*mpl_mst_air(idx)
     mpl_CFC12(idx)=q_CFC12(idx)*mpl_mst_air(idx)
     mpl_O3(idx)=q_O3(idx)*mpl_mst_air(idx)
     mpl_NO2(idx)=q_NO2(idx)*mpl_mst_air(idx)
     mpl_OH(idx)=q_OH(idx)*mpl_mst_air(idx)
     mpl_dry_air(idx)=mpl_mst_air(idx)-mpl_H2O(idx)
     mpl_O2(idx)=r_O2(idx)*mpl_dry_air(idx)
     mpl_N2(idx)=r_N2(idx)*mpl_dry_air(idx)
     dns_CO2(idx)=mpl_CO2(idx)/alt_dlt(idx)
     dns_CH4(idx)=mpl_CH4(idx)/alt_dlt(idx)
     dns_N2O(idx)=mpl_N2O(idx)/alt_dlt(idx)
     dns_CFC11(idx)=mpl_CFC11(idx)/alt_dlt(idx)
     dns_CFC12(idx)=mpl_CFC12(idx)/alt_dlt(idx)
     dns_H2O(idx)=mpl_H2O(idx)/alt_dlt(idx)
     dns_N2(idx)=mpl_N2(idx)/alt_dlt(idx)
     dns_NO2(idx)=mpl_NO2(idx)/alt_dlt(idx)
     dns_O2(idx)=mpl_O2(idx)/alt_dlt(idx)
     dns_O3(idx)=mpl_O3(idx)/alt_dlt(idx)
     dns_OH(idx)=mpl_OH(idx)/alt_dlt(idx)
     dns_dry_air(idx)=mpl_dry_air(idx)/alt_dlt(idx)
     q_O2(idx)=mpl_O2(idx)/mpl_mst_air(idx)
     q_N2(idx)=mpl_N2(idx)/mpl_mst_air(idx)
     ppr_H2O(idx)=prs(idx)*r_H2O(idx)/(eps_H2O+r_H2O(idx))
     ppr_dry_air(idx)=prs(idx)-ppr_H2O(idx)
     ppr_CO2(idx)=dns_CO2(idx)*gas_cst_CO2*tpt(idx)
     ppr_CH4(idx)=dns_CH4(idx)*gas_cst_CH4*tpt(idx)
     ppr_N2O(idx)=dns_N2O(idx)*gas_cst_N2O*tpt(idx)
     ppr_CFC11(idx)=dns_CFC11(idx)*gas_cst_CFC11*tpt(idx)
     ppr_CFC12(idx)=dns_CFC12(idx)*gas_cst_CFC12*tpt(idx)
     ppr_O3(idx)=dns_O3(idx)*gas_cst_O3*tpt(idx)
     ppr_NO2(idx)=dns_NO2(idx)*gas_cst_NO2*tpt(idx)
     ppr_OH(idx)=dns_OH(idx)*gas_cst_OH*tpt(idx)
     ppr_O2(idx)=dns_O2(idx)*gas_cst_O2*tpt(idx)
     ppr_N2(idx)=dns_N2(idx)*gas_cst_N2*tpt(idx)
     npl_mst_air(idx)=mpl_mst_air(idx)*Avagadro/mmw_mst_air(idx)
     npl_dry_air(idx)=mpl_dry_air(idx)*Avagadro/mmw_dry_air
     npl_H2O(idx)=mpl_H2O(idx)*Avagadro/mmw_H2O
     npl_CO2(idx)=mpl_CO2(idx)*Avagadro/mmw_CO2
     npl_CH4(idx)=mpl_CH4(idx)*Avagadro/mmw_CH4
     npl_N2O(idx)=mpl_N2O(idx)*Avagadro/mmw_N2O
     npl_CFC11(idx)=mpl_CFC11(idx)*Avagadro/mmw_CFC11
     npl_CFC12(idx)=mpl_CFC12(idx)*Avagadro/mmw_CFC12
     npl_O3(idx)=mpl_O3(idx)*Avagadro/mmw_O3
     npl_NO2(idx)=mpl_NO2(idx)*Avagadro/mmw_NO2
     npl_OH(idx)=mpl_OH(idx)*Avagadro/mmw_OH
     npl_O2(idx)=mpl_O2(idx)*Avagadro/mmw_O2
     npl_N2(idx)=mpl_N2(idx)*Avagadro/mmw_N2
     vmr_H2O(idx)=q_H2O(idx)*mmw_dry_air/mmw_H2O
     vmr_CO2(idx)=q_CO2(idx)*mmw_dry_air/mmw_CO2
     vmr_CH4(idx)=q_CH4(idx)*mmw_dry_air/mmw_CH4
     vmr_N2O(idx)=q_N2O(idx)*mmw_dry_air/mmw_N2O
     vmr_CFC11(idx)=q_CFC11(idx)*mmw_dry_air/mmw_CFC11
     vmr_CFC12(idx)=q_CFC12(idx)*mmw_dry_air/mmw_CFC12
     vmr_O2(idx)=q_O2(idx)*mmw_dry_air/mmw_O2
     vmr_N2(idx)=q_N2(idx)*mmw_dry_air/mmw_N2
     vmr_O3(idx)=q_O3(idx)*mmw_dry_air/mmw_O3
     vmr_NO2(idx)=q_NO2(idx)*mmw_dry_air/mmw_NO2
     vmr_OH(idx)=q_OH(idx)*mmw_dry_air/mmw_OH
     cnc_mst_air(idx)=npl_mst_air(idx)/alt_dlt(idx)
     cnc_dry_air(idx)=npl_dry_air(idx)/alt_dlt(idx)
     cnc_H2O(idx)=npl_H2O(idx)/alt_dlt(idx)
     cnc_CO2(idx)=npl_CO2(idx)/alt_dlt(idx)
     cnc_CH4(idx)=npl_CH4(idx)/alt_dlt(idx)
     cnc_N2O(idx)=npl_N2O(idx)/alt_dlt(idx)
     cnc_CFC11(idx)=npl_CFC11(idx)/alt_dlt(idx)
     cnc_CFC12(idx)=npl_CFC12(idx)/alt_dlt(idx)
     cnc_O2(idx)=npl_O2(idx)/alt_dlt(idx)
     cnc_N2(idx)=npl_N2(idx)/alt_dlt(idx)
     cnc_O3(idx)=npl_O3(idx)/alt_dlt(idx)
     cnc_NO2(idx)=npl_NO2(idx)/alt_dlt(idx)
     cnc_OH(idx)=npl_OH(idx)/alt_dlt(idx)
     ! sbc : Interested in adiabatic heating due to mixing down to surface rather
     ! than 1000mb (which could be above or below) surface --> use prs_0=prs_sfc
     ! This formulation takes into consideration the small influence of moisture 
     ! on the adiabatic cooling/warming of the parcel
     ! From: Emanuel, 1994 pg111 and Smith, 1997 pg34
     kappa=(gas_cst_dry_air/spc_heat_dry_air)*((1.0+(q_H2O(idx)/eps_H2O))/ &
          (1.0+(q_H2O(idx)*spc_heat_H2O_vpr/spc_heat_dry_air)))
     tpt_ptn(idx)=tpt(idx)*(prs_sfc/prs(idx))**kappa
     tpt_ptn_vrt(idx)=tpt_vrt(idx)*(prs_sfc/prs(idx))**(gas_cst_dry_air/spc_heat_dry_air)
     ! !sbc
  enddo
  
  ! Compute column totals
  mpc_CO2=0.0               ! [kg m-2]
  mpc_CH4=0.0               ! [kg m-2]
  mpc_N2O=0.0               ! [kg m-2]
  mpc_CFC11=0.0             ! [kg m-2]
  mpc_CFC12=0.0             ! [kg m-2]
  mpc_H2O=0.0               ! [kg m-2]
  mpc_OH=0.0                ! [kg m-2]
  mpc_O2=0.0                ! [kg m-2]
  mpc_N2=0.0                ! [kg m-2]
  mpc_O3=0.0                ! [kg m-2]
  mpc_NO2=0.0               ! [kg m-2]
  mpc_dry_air=0.0           ! [kg m-2]
  mpc_mst_air=0.0           ! [kg m-2]
  do idx=1,lev_nbr
     mpc_CO2=mpc_CO2+mpl_CO2(idx) ! [kg m-2]
     mpc_CH4=mpc_CH4+mpl_CH4(idx) ! [kg m-2]
     mpc_N2O=mpc_N2O+mpl_N2O(idx) ! [kg m-2]
     mpc_CFC11=mpc_CFC11+mpl_CFC11(idx) ! [kg m-2]
     mpc_CFC12=mpc_CFC12+mpl_CFC12(idx) ! [kg m-2]
     mpc_H2O=mpc_H2O+mpl_H2O(idx) ! [kg m-2]
     mpc_OH=mpc_OH+mpl_OH(idx) ! [kg m-2]
     mpc_O2=mpc_O2+mpl_O2(idx) ! [kg m-2]
     mpc_N2=mpc_N2+mpl_N2(idx) ! [kg m-2]
     mpc_O3=mpc_O3+mpl_O3(idx) ! [kg m-2]
     mpc_NO2=mpc_NO2+mpl_NO2(idx) ! [kg m-2]
     mpc_dry_air=mpc_dry_air+mpl_dry_air(idx) ! [kg m-2]
     mpc_mst_air=mpc_mst_air+mpl_mst_air(idx) ! [kg m-2]
  enddo
  
  npc_CO2=0.0               ! [mlc m-2]
  npc_CH4=0.0               ! [mlc m-2]
  npc_N2O=0.0               ! [mlc m-2]
  npc_CFC11=0.0             ! [mlc m-2]
  npc_CFC12=0.0             ! [mlc m-2]
  npc_H2O=0.0               ! [mlc m-2]
  npc_OH=0.0                ! [mlc m-2]
  npc_O2=0.0                ! [mlc m-2]
  npc_N2=0.0                ! [mlc m-2]
  npc_O3=0.0                ! [mlc m-2]
  npc_NO2=0.0               ! [mlc m-2]
  npc_dry_air=0.0           ! [mlc m-2]
  npc_mst_air=0.0           ! [mlc m-2]
  do idx=1,lev_nbr
     npc_CO2=npc_CO2+npl_CO2(idx) ! [mlc m-2]
     npc_CH4=npc_CH4+npl_CH4(idx) ! [mlc m-2]
     npc_N2O=npc_N2O+npl_N2O(idx) ! [mlc m-2]
     npc_CFC11=npc_CFC11+npl_CFC11(idx) ! [mlc m-2]
     npc_CFC12=npc_CFC12+npl_CFC12(idx) ! [mlc m-2]
     npc_H2O=npc_H2O+npl_H2O(idx) ! [mlc m-2]
     npc_OH=npc_OH+npl_OH(idx) ! [mlc m-2]
     npc_O2=npc_O2+npl_O2(idx) ! [mlc m-2]
     npc_N2=npc_N2+npl_N2(idx) ! [mlc m-2]
     npc_O3=npc_O3+npl_O3(idx) ! [mlc m-2]
     npc_NO2=npc_NO2+npl_NO2(idx) ! [mlc m-2]
     npc_dry_air=npc_dry_air+npl_dry_air(idx) ! [mlc m-2]
     npc_mst_air=npc_mst_air+npl_mst_air(idx) ! [mlc m-2]
  enddo                     ! end loop over lev
  
  ! Compute chemical equilibria separately, since they involve rate coefficients
  ! and may someday be computed in a separate, simple box model chemistry routine
  mpc_O2O2=0.0              ! [kg m-2]
  npc_O2O2=0.0              ! [mlc m-2]
  do idx=1,lev_nbr
     cnc_O2O2(idx)=(k_O2_O2*cnc_O2(idx))*cnc_O2(idx) ! [mlc m-3] force multiplication to avoid overflow in single precision
     npl_O2O2(idx)=cnc_O2O2(idx)*alt_dlt(idx) ! [mlc m-2]
     mpl_O2O2(idx)=npl_O2O2(idx)*mmw_O2O2/Avagadro ! [kg m-2]
     dns_O2O2(idx)=mpl_O2O2(idx)/alt_dlt(idx) ! [kg m-3]
     q_O2O2(idx)=mpl_O2O2(idx)/mpl_mst_air(idx) ! [kg kg-1]
     r_O2O2(idx)=q_O2O2(idx)/(1.0-q_O2O2(idx)) ! [kg kg-1]
     ppr_O2O2(idx)=q_O2O2(idx)*gas_cst_O2O2*tpt(idx) ! [Pa]
     vmr_O2O2(idx)=q_O2O2(idx)*mmw_dry_air/mmw_O2O2 ! [mlc mlc-1]
     mpc_O2O2=mpc_O2O2+mpl_O2O2(idx) ! [kg m-2]
     npc_O2O2=npc_O2O2+npl_O2O2(idx) ! [mlc m-2]
  enddo                     ! end loop over lev
  
  ! Quadratic number concentrations
  cnc_O2_npl_O2_clm=0.0     ! [mlc2 m-5]
  cnc_O2_npl_N2_clm=0.0     ! [mlc2 m-5]
  do idx=1,lev_nbr
     cnc_O2_cnc_O2(idx)=dble(cnc_O2(idx))*dble(cnc_O2(idx)) ! [mlc2 m-6]
     cnc_O2_cnc_N2(idx)=dble(cnc_O2(idx))*dble(cnc_N2(idx)) ! [mlc2 m-6]
     cnc_O2_npl_O2(idx)=cnc_O2_cnc_O2(idx)*alt_dlt(idx) ! [mlc2 m-5]
     cnc_O2_npl_N2(idx)=cnc_O2_cnc_N2(idx)*alt_dlt(idx) ! [mlc2 m-5]
     cnc_O2_npl_O2_clm=cnc_O2_npl_O2_clm+cnc_O2_npl_O2(idx) ! [mlc2 m-5]
     cnc_O2_npl_N2_clm=cnc_O2_npl_N2_clm+cnc_O2_npl_N2(idx) ! [mlc2 m-5]
  enddo                     ! end loop over lev
  ! Partial quadratic column paths
  cnc_O2_npl_O2_clm_frc(1)=cnc_O2_npl_O2(1) ! [mlc2 m-5]
  do idx=2,lev_nbr          ! NB: loop starts at 2
     cnc_O2_npl_O2_clm_frc(idx)=cnc_O2_npl_O2_clm_frc(idx-1)+cnc_O2_npl_O2(idx) ! [mlc2 m-5]
  enddo                     ! end loop over lev
  do idx=1,lev_nbr
     cnc_O2_npl_O2_clm_frc(idx)=cnc_O2_npl_O2_clm_frc(idx)/cnc_O2_npl_O2_clm ! [mlc2 m-5] -> [frc]
  enddo                     ! end loop over lev
  
  ! Quadratic mass concentrations
  dns_O2_mpl_O2_clm=0.0     ! [kg2 m-5]
  dns_O2_mpl_N2_clm=0.0     ! [kg2 m-5]
  do idx=1,lev_nbr
     dns_O2_dns_O2(idx)=dns_O2(idx)*dns_O2(idx) ! [kg2 m-6]
     dns_O2_dns_N2(idx)=dns_O2(idx)*dns_N2(idx) ! [kg2 m-6]
     dns_O2_mpl_O2(idx)=dns_O2_dns_O2(idx)*alt_dlt(idx) ! [kg2 m-5]
     dns_O2_mpl_N2(idx)=dns_O2_dns_N2(idx)*alt_dlt(idx) ! [kg2 m-5]
     dns_O2_mpl_O2_clm=dns_O2_mpl_O2_clm+dns_O2_mpl_O2(idx) ! [kg2 m-5]
     dns_O2_mpl_N2_clm=dns_O2_mpl_N2_clm+dns_O2_mpl_N2(idx) ! [kg2 m-5]
  enddo                     ! end loop over lev
  
  ! Compute chemical equilibria separately, since they involve rate coefficients
  ! and may someday be computed in a separate, simple box model chemistry routine
  ! call cnc_H2OH2O_CFT99(lev_nbr,alt_dlt,tpt,RH_lqd,mpl_mst_air,mpl_H2O,cnc_H2OH2O) 
  call cnc_H2OH2O_Chy97(lev_nbr,alt_dlt,tpt,RH_lqd,mpl_H2O,cnc_H2OH2O)
  ! call cnc_H2OH2O_Mun97(lev_nbr,alt_dlt,q_H2O,prs,tpt,RH_lqd,mpl_mst_air,cnc_H2OH2O)
  ! call cnc_H2OH2O_Zen97(RH_lqd,cnc_H2OH2O,cnc_mst_air,alt_dlt,lev_nbr,mpl_CWP,mpl_mst_air,prs,q_H2O,tpt)
  
  mpc_H2OH2O=0.0
  npc_H2OH2O=0.0
  do idx=1,lev_nbr
     npl_H2OH2O(idx)=cnc_H2OH2O(idx)*alt_dlt(idx)
     mpl_H2OH2O(idx)=npl_H2OH2O(idx)*mmw_H2OH2O/Avagadro
     dns_H2OH2O(idx)=mpl_H2OH2O(idx)/alt_dlt(idx)
     q_H2OH2O(idx)=mpl_H2OH2O(idx)/mpl_mst_air(idx)
     q_H2OH2O_rcp_q_H2O(idx)=q_H2OH2O(idx)/q_H2O(idx)
     r_H2OH2O(idx)=q_H2OH2O(idx)/(1.0-q_H2OH2O(idx))
     ppr_H2OH2O(idx)=q_H2OH2O(idx)*gas_cst_H2OH2O*tpt(idx)
     vmr_H2OH2O(idx)=q_H2OH2O(idx)*mmw_dry_air/mmw_H2OH2O
     mpc_H2OH2O=mpc_H2OH2O+mpl_H2OH2O(idx)
     npc_H2OH2O=npc_H2OH2O+npl_H2OH2O(idx)
  enddo                     ! end loop over lvl
  
  ! Integral sanity checks
  do idx=1,lev_nbr
     prs_ttl(idx)=ppr_O2(idx)+ppr_N2(idx)+ & !ppr_dry_air(idx)+
          ppr_H2O(idx)+ppr_CO2(idx)+ppr_CH4(idx)+ppr_N2O(idx)+ppr_CFC11(idx)+ppr_CFC12(idx)+ &
          ppr_O3(idx)+ppr_NO2(idx)+ &
          ppr_OH(idx)
     if (prs_ttl(idx) > prs(idx)) then
        write (0,'(2a,2(a,i4,a,f9.3,a))')  &
             prg_nm(1:ftn_strlen(prg_nm)),': WARNING ', &
             'prs_ttl(',idx,') = ',prs_ttl(idx)/100.0,' mb > ', &
             'prs(',idx,') = ',prs(idx)/100.0,' mb'
     endif                  ! endif dbg
  enddo                     ! end loop over lvl
  
  if (dbg_lvl == dbg_crr) then
     write (6,'(a,i3)') 'Integral partial pressure check:'
     write (6,'(a3,1x,14(a8,1x))') 'idx', &
          'prs','prs_ttl','ppr_dry','ppr_N2','ppr_O2', &
          'ppr_H2O','ppr_CO2','ppr_O3','ppr_NO2','ppr_OH', &
          'ppr_CH4','ppr_N2O','ppr_CFC11','ppr_CFC12'
     write (6,'(a3,1x,14(a8,1x))') '', &
          'mb','mb','mb','mb','mb', &
          'mb','mb','mb','mb','mb', &
          'mb','mb','mb','mb'
     do idx=1,lev_nbr
        write (6,'(i3,1x,14(f8.3,1x))') idx, &
             prs(idx)/100.0,prs_ttl(idx)/100.0,ppr_dry_air(idx)/100.0,ppr_N2(idx)/100.0,ppr_O2(idx)/100.0, &
             ppr_H2O(idx)/100.0,ppr_CO2(idx)/100.0,ppr_O3(idx)/100.0,ppr_NO2(idx)/100.0,ppr_OH(idx)/100.0, &
             ppr_CH4(idx)/100.0,ppr_N2O(idx)/100.0,ppr_CFC11(idx)/100.0,ppr_CFC12(idx)/100.0
     enddo                  ! end loop over lvl
  endif                     ! end if dbg
  
  ! Compute total mass paths in Dobson units (DU), which are defined as the
  ! number of milli-centimeters (i.e., 1.0e-5 m) which gas column would have
  ! if it were at standard temperature and pressure. For O3, column total
  ! is usually O(100) DU.
  mpc_O3_DU=mpc_O3*gas_cst_O3*tpt_STP/prs_STP ! [m]
  mpc_O3_DU=mpc_O3_DU*1.0e5 ! [m] --> [millicm] = [m-5]
  
  ! cld_top and cld_btm are defined on interface levels but cld_mid is defined at layer midpoint
  cld_top_lvl=0
  cld_btm_lvl=0
  do idx=1,lev_nbr
     if (mpl_CWP(idx) > 0.0) then
        if (force_cld_lvl_2b_overcast) then
           cld_frc(idx)=1.0
        endif ! endif 
        cld_btm_lvl=idx
        if (cld_top_lvl == 0) then 
           cld_top_lvl=idx
        endif               ! end if cloud top doesn't exist yet
     endif                  ! end if layer is cloudy
  enddo                     ! end loop over lev
  if (cld_top_lvl /= 0) then
     prs_cld_btm=prs_ntf(cld_btm_lvl+1)
     alt_cld_btm=alt_ntf(cld_btm_lvl+1)
     prs_cld_top=prs_ntf(cld_top_lvl)
     alt_cld_top=alt_ntf(cld_top_lvl)
  else                      ! endif clouds were found
     prs_cld_btm=0.0
     alt_cld_btm=0.0
     prs_cld_top=0.0
     alt_cld_top=0.0
  endif                     ! endif no clouds were found
  alt_cld_thick=alt_cld_top-alt_cld_btm
  prs_cld_thick=prs_cld_btm-prs_cld_top
  alt_cld_mid=alt_cld_btm+alt_cld_thick/2.0
  
  ! Should prs_cld_mid be pressure at geometric cloud midpoint (alt_cld_mid), or pressure midpoint of cloud?  
  prs_cld_mid=prs_cld_top+prs_cld_thick/2.0
  
  if (zen_ngl_dpn_alb) then
     ! Fit to ARM SGP CART SIROS data of 95/10/11 and 1995/10/15 by BPB in Bri96:
     ! alb_sfc=0.215*1.6/(1.0+1.2*slr_zen_ngl_cos)
     ! Modified to produce alb_sfc = 0.17 at 1200 TST 1995/10/11 by CSZ in ZBP97:
     alb_sfc=0.2*1.6/(1.0+1.2*slr_zen_ngl_cos)
     ! Added 19961030: Partition BB (broadband) albedo between VIS and NIR:
     ! Bon96 prescribes alb(NIR)=M*alb(VIS) where M > 1 depends on land surface
     ! For grassland in October, LAI ~ SAI (Leaf Area Index, Stem Area Index) and M=2.47
     ! Assuming flx_VIS_dwn_sfc ~ flx_NIR_dwn_sfc, we can partition a given BB albedo
     ! into VIS and NIR (knowing M) and be fairly sure resulting model will have same BB albedo.
     alb_NIR_vis_rat=2.47
     alb_sfc_vsb_drc=(2.0/(alb_NIR_vis_rat+1.0))*alb_sfc
     alb_sfc_vsb_dff=(2.0/(alb_NIR_vis_rat+1.0))*alb_sfc
     alb_sfc_NIR_drc=(2.0*alb_NIR_vis_rat/(alb_NIR_vis_rat+1.0))*alb_sfc
     alb_sfc_NIR_dff=(2.0*alb_NIR_vis_rat/(alb_NIR_vis_rat+1.0))*alb_sfc
  endif                     ! end if overriding CLM albedo
  alb_sfc=0.5*(alb_sfc_vsb_drc+alb_sfc_NIR_drc) ! [frc] Surface albedo
  
  ! sbc
  ! Calculate wind speed bin weights and speeds
  call wnd_spd_PDF(PDF_avg_wnd_spd,PDF_spd_var,PDF_bin_nbr,PDF_bin_wgt,PDF_bin_spd) ! Weibull PDF
  
  ! Calculate lifted parcel profile and other column thermodynamic quantities 
  call pcl_tdy_pfl(alt, &
       grv, &
       spc_heat_mst_air, &
       tpt, &
       tpt_sfc, &
       CAPE, &
       CINE, &
       lev_LCL, &
       lev_LFC, &
       lev_LNB, &
       lev_nbr, &
       nrg_dry, &
       nrg_mst, &
       pcl_DLR, &
       pcl_MLR, &
       pcl_PLR, &
       prs, &
       prs_sfc, &
       r_H2O, &
       RH, &
       tpt_dwp, &
       tpt_pcl, &
       tpt_pcl_sfc)
  if (flg_wnd) then
     ! Diagnose ABL turbulence
     call abl_trb_pfl(dbg_lvl,lev_nbr,oro,alt,tpt,tpt_ptn,tpt_ptn_vrt,grv, & ! I
          tpt_skn,q_H2O,prs,prs_sfc,dns_mst_air,wnd_dir,wnd_spd,q_H2O_sfc, & ! I
          wnd_shr,ric_nbr,tke,khfs,kbfs,kqfs,shflx,lhflx,ustar, & ! O
          obk_len,wstar,abl_hgt,wspd10m,tau) ! O
  end if ! !flg_wnd 
  ! !sbc
  
  ! Print data in CRM text style
  open (fl_txt_unit,file=fl_txt,status='unknown',iostat=rcd)
  if (rcd /= 0) write (6,'(4a,i4)') prg_nm(1:ftn_strlen(prg_nm)),': ERROR unable to open ',fl_txt(1:ftn_strlen(fl_txt))
  write (fl_txt_unit,*) 'CCM Column Radiation Model (CRM) Input File -*-text-*-'
  write (fl_txt_unit,*) prf_sng(1:ftn_strlen(prf_sng))
  write (fl_txt_unit,*) 'CRM homepage URL is http://www.cgd.ucar.edu/cms/crm'
  ! 19990703: PGI compiler breaks lcl_yr_day line into two lines because it exceeds 80 characters
  ! Solution use explicitly formatted output rather than asterisks (list-directed output)
  write (fl_txt_unit,'(a,f15.11,1x,a)') '   ',lcl_yr_day,'Julian day of year (1.5 = Noon, Jan 1; from 1 to 365)'
  write (fl_txt_unit,*) lat*180.0/pi,' Latitude (degrees North, from -90.0 to +90.0)'
  if (flg_aer) then
     write (fl_txt_unit,*) 'Level   p [mb]       T [K]    H2O mmr    O3 mmr     Cld frc.  Cld CWP    Aer (odx)  '
  else
     write (fl_txt_unit,*) 'Level   p [mb]       T [K]    H2O mmr    O3 mmr     Cld frc.  Cld CWP               '
  endif
  do idx=1,lev_nbr
     if (flg_aer) then
        write (fl_txt_unit,98) &
             idx, &
             prs(idx)/100.0, & ! [Pa] -> [mb]
             tpt(idx), &
             q_H2O(idx), &
             q_O3(idx), &
             cld_frc(idx), &
             mpl_CWP(idx)*1000.0, & ! [kg m-2] -> [g m-2]
             odxl_obs_aer(idx)
98      format(1x,i3,1x,es12.5,1x,6(es10.3,1x))
     else
        write (fl_txt_unit,99) &
             idx, &
             prs(idx)/100.0, & ! [Pa] -> [mb]
             tpt(idx), &
             q_H2O(idx), &
             q_O3(idx), &
             cld_frc(idx), &
             mpl_CWP(idx)*1000.0 ! [kg m-2] -> [g m-2]
99      format(1x,i3,1x,es12.5,1x,5(es10.3,1x))
     endif
  enddo ! end loop over lev
  write (fl_txt_unit,*) prs_sfc/100.0,' Surface pressure [mb]' ! [Pa] -> [mb]
  write (fl_txt_unit,*) tpt_sfc,' Surface air temperature [K]' ! [K]
  write (fl_txt_unit,*) tpt_skn,' Ground (skin) temperature [K]' ! [K]
  write (fl_txt_unit,*) oro,' Surface type flag (0=ocn, 1=lnd, 2=sea ice)'
  write (fl_txt_unit,*) rgh_len,' Surface aerodynamic roughness [m] (obsolete)'
  write (fl_txt_unit,*) snow_depth,' Snow cover liquid water equivalent [m]'
  write (fl_txt_unit,*) alb_sfc_vsb_drc,' Albedo (Vis, direct)'
  write (fl_txt_unit,*) alb_sfc_vsb_dff,' Albedo (Vis, diffuse)'
  write (fl_txt_unit,*) alb_sfc_NIR_drc,' Albedo (NIR, direct)'
  write (fl_txt_unit,*) alb_sfc_NIR_dff,' Albedo (NIR, diffuse)'
  write (fl_txt_unit,*) frc_str_zen_ngl_sfc,' Fraction strong zenith angle dep. sfc. (obsolete)'
  write (fl_txt_unit,*) CO2_vmr_clm,' CO2 volume mixing ratio'
  write (fl_txt_unit,*) N2O_vmr_clm,' N2O volume mixing ratio'
  write (fl_txt_unit,*) CH4_vmr_clm,' CH4 volume mixing ratio'
  write (fl_txt_unit,*) CFC11_vmr_clm,' CFC11 volume mixing ratio'
  write (fl_txt_unit,*) CFC12_vmr_clm,' CFC12 volume mixing ratio'
  write (fl_txt_unit,*) odxc_obs_aer,' Aerosol visible opt. depth'
  write (fl_txt_unit,*) slr_cst,' Solar constant [W m-2]'
  write (fl_txt_unit,*) gmt_yr,' Year AD'
  write (fl_txt_unit,*) lon_dgr,' Longitude (degrees East, from 0.0 to 360.0)'
  if (flg_snw) then
     write (fl_txt_unit,*) 'CLM 2.0 SWNB 2.0 Snow Column Radiation Model (SCRM) Input File -*-text-*-'
     write (fl_txt_unit,*) prf_snw_sng(1:ftn_strlen(prf_snw_sng))
     write (fl_txt_unit,*) '  Level Thick	   Temper.  Density    BC mmr	  Eff Rds'
     write (fl_txt_unit,*) '    #     cm	   K	    g cm-3     mmr	  um'
     do lev_snw_idx=1,lev_snw_nbr
        write (fl_txt_unit,97) &
             lev_snw_idx, &
             dpt_dlt_snw(lev_snw_idx)*100.0, & ! [m] -> [cm]
             tpt(lev_snw_idx), &
             dns_snw(lev_snw_idx), &
             mmr_mpr_snw(lev_snw_idx), &
             rds_ffc_snw(lev_snw_idx)*1.0e6 ! [m] -> [um]
97      format(1x,i3,1x,6(es10.3,1x))
     enddo ! end loop over lev
  endif ! !flg_snw
  close (fl_txt_unit)
  ! Finish writing CRM style text output file
  
  ! Begin netCDF output routines
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
  
  ! First, define dimension IDs
  rcd=nf90_wrp(nf90_def_dim(nc_id,'lev',lev_nbr,lev_dmn_id),sbr_nm//": def_dim lev in "//__FILE__)
  rcd=nf90_wrp(nf90_def_dim(nc_id,'levp',levp_nbr,levp_dmn_id),sbr_nm//": def_dim levp in "//__FILE__)
  if (flg_snw) then
     rcd=nf90_wrp(nf90_def_dim(nc_id,'lev_snw',lev_snw_nbr,lev_snw_dmn_id),sbr_nm//": def_dim lev_snw in "//__FILE__)
     rcd=nf90_wrp(nf90_def_dim(nc_id,'levp_snw',levp_snw_nbr,levp_snw_dmn_id),sbr_nm//": def_dim levp_snw in "//__FILE__)
  endif ! !flg_snw
  rcd=nf90_wrp(nf90_def_dim(nc_id,'PDF',PDF_bin_nbr,PDF_bin_dmn_id),sbr_nm//": def_dim PDF in "//__FILE__)
  
  ! Variable definitions
  rcd=nf90_wrp(nf90_def_var(nc_id,'CFC11_vmr_clm',nf90_float,CFC11_vmr_clm_id),sbr_nm//": dv CFC11_vmr_clm")
#ifdef ENABLE_NETCDF4
  ! Set HDF Lempel-Ziv compression level, if requested
  if (dfl_lvl > 0) rcd=nf90_def_var_deflate(nc_id,CFC11_vmr_clm_id,flg_shf,flg_dfl,dfl_lvl)
#endif /* !ENABLE_NETCDF4 */
  rcd=nf90_wrp(nf90_def_var(nc_id,'CAPE',nf90_float,CAPE_id),sbr_nm//": dv CAPE")
  rcd=nf90_wrp(nf90_def_var(nc_id,'CFC12_vmr_clm',nf90_float,CFC12_vmr_clm_id),sbr_nm//": dv CFC12_vmr_clm")
  rcd=nf90_wrp(nf90_def_var(nc_id,'CH4_vmr_clm',nf90_float,CH4_vmr_clm_id),sbr_nm//": dv CH4_vmr_clm")
  rcd=nf90_wrp(nf90_def_var(nc_id,'CINE',nf90_float,CINE_id),sbr_nm//": dv CINE")
  rcd=nf90_wrp(nf90_def_var(nc_id,'CO2_vmr_clm',nf90_float,CO2_vmr_clm_id),sbr_nm//": dv CO2_vmr_clm")
  rcd=nf90_wrp(nf90_def_var(nc_id,'N2O_vmr_clm',nf90_float,N2O_vmr_clm_id),sbr_nm//": dv N2O_vmr_clm")
  rcd=nf90_wrp(nf90_def_var(nc_id,'PDF_bin_spd',nf90_float,PDF_bin_dmn_id,PDF_bin_spd_id),sbr_nm//": dv PDF_bin_spd")
  rcd=nf90_wrp(nf90_def_var(nc_id,'PDF_bin_wgt',nf90_float,PDF_bin_dmn_id,PDF_bin_wgt_id),sbr_nm//": dv PDF_bin_wgt")
  rcd=nf90_wrp(nf90_def_var(nc_id,'RH',nf90_float,lev_dmn_id,RH_id),sbr_nm//": dv RH")
  rcd=nf90_wrp(nf90_def_var(nc_id,'RH_ice',nf90_float,lev_dmn_id,RH_ice_id),sbr_nm//": dv RH_ice")
  rcd=nf90_wrp(nf90_def_var(nc_id,'RH_lqd',nf90_float,lev_dmn_id,RH_lqd_id),sbr_nm//": dv RH_lqd")
  rcd=nf90_wrp(nf90_def_var(nc_id,'alb_sfc',nf90_float,alb_sfc_id),sbr_nm//": dv alb_sfc")
  rcd=nf90_wrp(nf90_def_var(nc_id,'alb_sfc_NIR_dff',nf90_float,alb_sfc_NIR_dff_id),sbr_nm//": dv alb_sfc_NIR_dff")
  rcd=nf90_wrp(nf90_def_var(nc_id,'alb_sfc_NIR_drc',nf90_float,alb_sfc_NIR_drc_id),sbr_nm//": dv alb_sfc_NIR_drc")
  rcd=nf90_wrp(nf90_def_var(nc_id,'alb_sfc_vsb_dff',nf90_float,alb_sfc_vsb_dff_id),sbr_nm//": dv alb_sfc_vsb_dff")
  rcd=nf90_wrp(nf90_def_var(nc_id,'alb_sfc_vsb_drc',nf90_float,alb_sfc_vsb_drc_id),sbr_nm//": dv alb_sfc_vsb_drc")
  rcd=nf90_wrp(nf90_def_var(nc_id,'alt',nf90_float,lev_dmn_id,alt_id),sbr_nm//": dv alt")
  rcd=nf90_wrp(nf90_def_var(nc_id,'alt_cld_btm',nf90_float,alt_cld_btm_id),sbr_nm//": dv alt_cld_btm")
  rcd=nf90_wrp(nf90_def_var(nc_id,'alt_cld_mid',nf90_float,alt_cld_mid_id),sbr_nm//": dv alt_cld_mid")
  rcd=nf90_wrp(nf90_def_var(nc_id,'alt_cld_thick',nf90_float,alt_cld_thick_id),sbr_nm//": dv alt_cld_thick")
  rcd=nf90_wrp(nf90_def_var(nc_id,'alt_cld_top',nf90_float,alt_cld_top_id),sbr_nm//": dv alt_cld_top")
  rcd=nf90_wrp(nf90_def_var(nc_id,'alt_dlt',nf90_float,lev_dmn_id,alt_dlt_id),sbr_nm//": dv alt_dlt")
  rcd=nf90_wrp(nf90_def_var(nc_id,'alt_ntf',nf90_float,levp_dmn_id,alt_ntf_id),sbr_nm//": dv alt_ntf")
  rcd=nf90_wrp(nf90_def_var(nc_id,'cld_frc',nf90_float,lev_dmn_id,cld_frc_id),sbr_nm//": dv cld_frc")
  rcd=nf90_wrp(nf90_def_var(nc_id,'cnc_CFC11',nf90_float,lev_dmn_id,cnc_CFC11_id),sbr_nm//": dv cnc_CFC11")
  rcd=nf90_wrp(nf90_def_var(nc_id,'cnc_CFC12',nf90_float,lev_dmn_id,cnc_CFC12_id),sbr_nm//": dv cnc_CFC12")
  rcd=nf90_wrp(nf90_def_var(nc_id,'cnc_CH4',nf90_float,lev_dmn_id,cnc_CH4_id),sbr_nm//": dv cnc_CH4")
  rcd=nf90_wrp(nf90_def_var(nc_id,'cnc_CO2',nf90_float,lev_dmn_id,cnc_CO2_id),sbr_nm//": dv cnc_CO2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'cnc_H2O',nf90_float,lev_dmn_id,cnc_H2O_id),sbr_nm//": dv cnc_H2O")
  rcd=nf90_wrp(nf90_def_var(nc_id,'cnc_H2OH2O',nf90_float,lev_dmn_id,cnc_H2OH2O_id),sbr_nm//": dv cnc_H2OH2O")
  rcd=nf90_wrp(nf90_def_var(nc_id,'cnc_N2',nf90_float,lev_dmn_id,cnc_N2_id),sbr_nm//": dv cnc_N2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'cnc_N2O',nf90_float,lev_dmn_id,cnc_N2O_id),sbr_nm//": dv cnc_N2O")
  rcd=nf90_wrp(nf90_def_var(nc_id,'cnc_NO2',nf90_float,lev_dmn_id,cnc_NO2_id),sbr_nm//": dv cnc_NO2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'cnc_O2',nf90_float,lev_dmn_id,cnc_O2_id),sbr_nm//": dv cnc_O2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'cnc_O2O2',nf90_float,lev_dmn_id,cnc_O2O2_id),sbr_nm//": dv cnc_O2O2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'cnc_O2_cnc_N2',nf90_double,lev_dmn_id,cnc_O2_cnc_N2_id),sbr_nm//": dv cnc_O2_cnc_N2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'cnc_O2_cnc_O2',nf90_double,lev_dmn_id,cnc_O2_cnc_O2_id),sbr_nm//": dv cnc_O2_cnc_O2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'cnc_O2_npl_N2',nf90_double,lev_dmn_id,cnc_O2_npl_N2_id),sbr_nm//": dv cnc_O2_npl_N2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'cnc_O2_npl_N2_clm',nf90_double,cnc_O2_npl_N2_clm_id),sbr_nm//": dv cnc_O2_npl_N2_clm")
  rcd=nf90_wrp(nf90_def_var(nc_id,'cnc_O2_npl_O2',nf90_double,lev_dmn_id,cnc_O2_npl_O2_id),sbr_nm//": dv cnc_O2_npl_O2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'cnc_O2_npl_O2_clm',nf90_double,cnc_O2_npl_O2_clm_id),sbr_nm//": dv cnc_O2_npl_O2_clm")
  rcd=nf90_wrp(nf90_def_var(nc_id,'cnc_O3',nf90_float,lev_dmn_id,cnc_O3_id),sbr_nm//": dv cnc_O3")
  rcd=nf90_wrp(nf90_def_var(nc_id,'cnc_OH',nf90_float,lev_dmn_id,cnc_OH_id),sbr_nm//": dv cnc_OH")
  rcd=nf90_wrp(nf90_def_var(nc_id,'cnc_dry_air',nf90_float,lev_dmn_id,cnc_dry_air_id),sbr_nm//": dv cnc_dry_air")
  rcd=nf90_wrp(nf90_def_var(nc_id,'cnc_mst_air',nf90_float,lev_dmn_id,cnc_mst_air_id),sbr_nm//": dv cnc_mst_air")
  rcd=nf90_wrp(nf90_def_var(nc_id,'dns_CFC11',nf90_float,lev_dmn_id,dns_CFC11_id),sbr_nm//": dv dns_CFC11")
  rcd=nf90_wrp(nf90_def_var(nc_id,'dns_CFC12',nf90_float,lev_dmn_id,dns_CFC12_id),sbr_nm//": dv dns_CFC12")
  rcd=nf90_wrp(nf90_def_var(nc_id,'dns_CH4',nf90_float,lev_dmn_id,dns_CH4_id),sbr_nm//": dv dns_CH4")
  rcd=nf90_wrp(nf90_def_var(nc_id,'dns_CO2',nf90_float,lev_dmn_id,dns_CO2_id),sbr_nm//": dv dns_CO2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'dns_H2O',nf90_float,lev_dmn_id,dns_H2O_id),sbr_nm//": dv dns_H2O")
  rcd=nf90_wrp(nf90_def_var(nc_id,'dns_H2OH2O',nf90_float,lev_dmn_id,dns_H2OH2O_id),sbr_nm//": dv dns_H2OH2O")
  rcd=nf90_wrp(nf90_def_var(nc_id,'dns_N2',nf90_float,lev_dmn_id,dns_N2_id),sbr_nm//": dv dns_N2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'dns_N2O',nf90_float,lev_dmn_id,dns_N2O_id),sbr_nm//": dv dns_N2O")
  rcd=nf90_wrp(nf90_def_var(nc_id,'dns_NO2',nf90_float,lev_dmn_id,dns_NO2_id),sbr_nm//": dv dns_NO2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'dns_O2',nf90_float,lev_dmn_id,dns_O2_id),sbr_nm//": dv dns_O2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'dns_O2O2',nf90_float,lev_dmn_id,dns_O2O2_id),sbr_nm//": dv dns_O2O2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'dns_O2_dns_N2',nf90_float,lev_dmn_id,dns_O2_dns_N2_id),sbr_nm//": dv dns_O2_dns_N2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'dns_O2_dns_O2',nf90_float,lev_dmn_id,dns_O2_dns_O2_id),sbr_nm//": dv dns_O2_dns_O2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'dns_O2_mpl_N2',nf90_float,lev_dmn_id,dns_O2_mpl_N2_id),sbr_nm//": dv dns_O2_mpl_N2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'dns_O2_mpl_N2_clm',nf90_float,dns_O2_mpl_N2_clm_id),sbr_nm//": dv dns_O2_mpl_N2_clm")
  rcd=nf90_wrp(nf90_def_var(nc_id,'dns_O2_mpl_O2',nf90_float,lev_dmn_id,dns_O2_mpl_O2_id),sbr_nm//": dv dns_O2_mpl_O2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'dns_O2_mpl_O2_clm',nf90_float,dns_O2_mpl_O2_clm_id),sbr_nm//": dv dns_O2_mpl_O2_clm")
  rcd=nf90_wrp(nf90_def_var(nc_id,'dns_O3',nf90_float,lev_dmn_id,dns_O3_id),sbr_nm//": dv dns_O3")
  rcd=nf90_wrp(nf90_def_var(nc_id,'dns_OH',nf90_float,lev_dmn_id,dns_OH_id),sbr_nm//": dv dns_OH")
  rcd=nf90_wrp(nf90_def_var(nc_id,'dns_aer',nf90_float,dns_aer_id),sbr_nm//": dv dns_aer")
  rcd=nf90_wrp(nf90_def_var(nc_id,'dns_bga',nf90_float,dns_bga_id),sbr_nm//": dv dns_bga")
  rcd=nf90_wrp(nf90_def_var(nc_id,'dns_dry_air',nf90_float,lev_dmn_id,dns_dry_air_id),sbr_nm//": dv dns_dry_air")
  rcd=nf90_wrp(nf90_def_var(nc_id,'dns_mst_air',nf90_float,lev_dmn_id,dns_mst_air_id),sbr_nm//": dv dns_mst_air")
  rcd=nf90_wrp(nf90_def_var(nc_id,'eqn_time_sec',nf90_double,eqn_time_sec_id),sbr_nm//": dv eqn_time_sec")
  rcd=nf90_wrp(nf90_def_var(nc_id,'ext_cff_mss_aer',nf90_float,ext_cff_mss_aer_id),sbr_nm//": dv ext_cff_mss_aer")
  rcd=nf90_wrp(nf90_def_var(nc_id,'ext_cff_mss_bga',nf90_float,ext_cff_mss_bga_id),sbr_nm//": dv ext_cff_mss_bga")
  rcd=nf90_wrp(nf90_def_var(nc_id,'frc_ice',nf90_float,lev_dmn_id,frc_ice_id),sbr_nm//": dv frc_ice")
  rcd=nf90_wrp(nf90_def_var(nc_id,'frc_ice_ttl',nf90_float,frc_ice_ttl_id),sbr_nm//": dv frc_ice_ttl")
  rcd=nf90_wrp(nf90_def_var(nc_id,'frc_str_zen_ngl_sfc',nf90_float,frc_str_zen_ngl_sfc_id),sbr_nm//": dv frc_str_zen_ngl_sfc")
  rcd=nf90_wrp(nf90_def_var(nc_id,'gas_cst_mst_air',nf90_float,lev_dmn_id,gas_cst_mst_air_id),sbr_nm//": dv gas_cst_mst_air")
  rcd=nf90_wrp(nf90_def_var(nc_id,'gmt_day',nf90_int,gmt_day_id),sbr_nm//": dv gmt_day")
  rcd=nf90_wrp(nf90_def_var(nc_id,'gmt_doy',nf90_double,gmt_doy_id),sbr_nm//": dv gmt_doy")
  rcd=nf90_wrp(nf90_def_var(nc_id,'gmt_hr',nf90_int,gmt_hr_id),sbr_nm//": dv gmt_hr")
  rcd=nf90_wrp(nf90_def_var(nc_id,'gmt_mnt',nf90_int,gmt_mnt_id),sbr_nm//": dv gmt_mnt")
  rcd=nf90_wrp(nf90_def_var(nc_id,'gmt_mth',nf90_int,gmt_mth_id),sbr_nm//": dv gmt_mth")
  rcd=nf90_wrp(nf90_def_var(nc_id,'gmt_sec',nf90_int,gmt_sec_id),sbr_nm//": dv gmt_sec")
  rcd=nf90_wrp(nf90_def_var(nc_id,'gmt_ydy',nf90_int,gmt_ydy_id),sbr_nm//": dv gmt_ydy")
  rcd=nf90_wrp(nf90_def_var(nc_id,'gmt_yr',nf90_int,gmt_yr_id),sbr_nm//": dv gmt_yr")
  rcd=nf90_wrp(nf90_def_var(nc_id,'grv',nf90_float,lev_dmn_id,grv_id),sbr_nm//": dv grv")
  rcd=nf90_wrp(nf90_def_var(nc_id,'lat',nf90_double,lat_id),sbr_nm//": dv lat")
  rcd=nf90_wrp(nf90_def_var(nc_id,'lat_cos',nf90_double,lat_cos_id),sbr_nm//": dv lat_cos")
  rcd=nf90_wrp(nf90_def_var(nc_id,'lat_dgr',nf90_double,lat_dgr_id),sbr_nm//": dv lat_dgr")
  rcd=nf90_wrp(nf90_def_var(nc_id,'lcl_time_hr',nf90_double,lcl_time_hr_id),sbr_nm//": dv lcl_time_hr")
  rcd=nf90_wrp(nf90_def_var(nc_id,'lcl_yr_day',nf90_double,lcl_yr_day_id),sbr_nm//": dv lcl_yr_day")
  rcd=nf90_wrp(nf90_def_var(nc_id,'lev',nf90_float,lev_dmn_id,lev_id),sbr_nm//": dv lev")
  rcd=nf90_wrp(nf90_def_var(nc_id,'lev_LCL',nf90_float,lev_LCL_id),sbr_nm//": dv lev_LCL")
  rcd=nf90_wrp(nf90_def_var(nc_id,'lev_LFC',nf90_float,lev_LFC_id),sbr_nm//": dv lev_LFC")
  rcd=nf90_wrp(nf90_def_var(nc_id,'lev_LNB',nf90_float,lev_LNB_id),sbr_nm//": dv lev_LNB")
  rcd=nf90_wrp(nf90_def_var(nc_id,'levp',nf90_float,levp_dmn_id,levp_id),sbr_nm//": dv levp")
  rcd=nf90_wrp(nf90_def_var(nc_id,'lmt_day',nf90_int,lmt_day_id),sbr_nm//": dv lmt_day")
  rcd=nf90_wrp(nf90_def_var(nc_id,'lmt_hr',nf90_int,lmt_hr_id),sbr_nm//": dv lmt_hr")
  rcd=nf90_wrp(nf90_def_var(nc_id,'lmt_mnt',nf90_int,lmt_mnt_id),sbr_nm//": dv lmt_mnt")
  rcd=nf90_wrp(nf90_def_var(nc_id,'lmt_mth',nf90_int,lmt_mth_id),sbr_nm//": dv lmt_mth")
  rcd=nf90_wrp(nf90_def_var(nc_id,'lmt_sec',nf90_int,lmt_sec_id),sbr_nm//": dv lmt_sec")
  rcd=nf90_wrp(nf90_def_var(nc_id,'lmt_ydy',nf90_int,lmt_ydy_id),sbr_nm//": dv lmt_ydy")
  rcd=nf90_wrp(nf90_def_var(nc_id,'lmt_yr',nf90_int,lmt_yr_id),sbr_nm//": dv lmt_yr")
  rcd=nf90_wrp(nf90_def_var(nc_id,'lon',nf90_double,lon_id),sbr_nm//": dv lon")
  rcd=nf90_wrp(nf90_def_var(nc_id,'lon_dgr',nf90_double,lon_dgr_id),sbr_nm//": dv lon_dgr")
  rcd=nf90_wrp(nf90_def_var(nc_id,'lon_sec',nf90_double,lon_sec_id),sbr_nm//": dv lon_sec")
  rcd=nf90_wrp(nf90_def_var(nc_id,'ltst_day',nf90_int,ltst_day_id),sbr_nm//": dv ltst_day")
  rcd=nf90_wrp(nf90_def_var(nc_id,'ltst_doy',nf90_double,ltst_doy_id),sbr_nm//": dv ltst_doy")
  rcd=nf90_wrp(nf90_def_var(nc_id,'ltst_hr',nf90_int,ltst_hr_id),sbr_nm//": dv ltst_hr")
  rcd=nf90_wrp(nf90_def_var(nc_id,'ltst_mnt',nf90_int,ltst_mnt_id),sbr_nm//": dv ltst_mnt")
  rcd=nf90_wrp(nf90_def_var(nc_id,'ltst_mth',nf90_int,ltst_mth_id),sbr_nm//": dv ltst_mth")
  rcd=nf90_wrp(nf90_def_var(nc_id,'ltst_sec',nf90_int,ltst_sec_id),sbr_nm//": dv ltst_sec")
  rcd=nf90_wrp(nf90_def_var(nc_id,'ltst_ydy',nf90_int,ltst_ydy_id),sbr_nm//": dv ltst_ydy")
  rcd=nf90_wrp(nf90_def_var(nc_id,'ltst_yr',nf90_int,ltst_yr_id),sbr_nm//": dv ltst_yr")
  rcd=nf90_wrp(nf90_def_var(nc_id,'mmw_mst_air',nf90_float,lev_dmn_id,mmw_mst_air_id),sbr_nm//": dv mmw_mst_air")
  rcd=nf90_wrp(nf90_def_var(nc_id,'mpc_CFC11',nf90_float,mpc_CFC11_id),sbr_nm//": dv mpc_CFC11")
  rcd=nf90_wrp(nf90_def_var(nc_id,'mpc_CFC12',nf90_float,mpc_CFC12_id),sbr_nm//": dv mpc_CFC12")
  rcd=nf90_wrp(nf90_def_var(nc_id,'mpc_CH4',nf90_float,mpc_CH4_id),sbr_nm//": dv mpc_CH4")
  rcd=nf90_wrp(nf90_def_var(nc_id,'mpc_CO2',nf90_float,mpc_CO2_id),sbr_nm//": dv mpc_CO2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'mpc_CWP',nf90_float,mpc_CWP_id),sbr_nm//": dv mpc_CWP")
  rcd=nf90_wrp(nf90_def_var(nc_id,'mpc_H2O',nf90_float,mpc_H2O_id),sbr_nm//": dv mpc_H2O")
  rcd=nf90_wrp(nf90_def_var(nc_id,'mpc_H2OH2O',nf90_float,mpc_H2OH2O_id),sbr_nm//": dv mpc_H2OH2O")
  rcd=nf90_wrp(nf90_def_var(nc_id,'mpc_IWP',nf90_float,mpc_IWP_id),sbr_nm//": dv mpc_IWP")
  rcd=nf90_wrp(nf90_def_var(nc_id,'mpc_LWP',nf90_float,mpc_LWP_id),sbr_nm//": dv mpc_LWP")
  rcd=nf90_wrp(nf90_def_var(nc_id,'mpc_N2',nf90_float,mpc_N2_id),sbr_nm//": dv mpc_N2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'mpc_N2O',nf90_float,mpc_N2O_id),sbr_nm//": dv mpc_N2O")
  rcd=nf90_wrp(nf90_def_var(nc_id,'mpc_NO2',nf90_float,mpc_NO2_id),sbr_nm//": dv mpc_NO2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'mpc_O2',nf90_float,mpc_O2_id),sbr_nm//": dv mpc_O2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'mpc_O2O2',nf90_float,mpc_O2O2_id),sbr_nm//": dv mpc_O2O2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'mpc_O3',nf90_float,mpc_O3_id),sbr_nm//": dv mpc_O3")
  rcd=nf90_wrp(nf90_def_var(nc_id,'mpc_O3_DU',nf90_float,mpc_O3_DU_id),sbr_nm//": dv mpc_O3_DU")
  rcd=nf90_wrp(nf90_def_var(nc_id,'mpc_OH',nf90_float,mpc_OH_id),sbr_nm//": dv mpc_OH")
  rcd=nf90_wrp(nf90_def_var(nc_id,'mpc_aer',nf90_float,mpc_aer_id),sbr_nm//": dv mpc_aer")
  rcd=nf90_wrp(nf90_def_var(nc_id,'mpc_bga',nf90_float,mpc_bga_id),sbr_nm//": dv mpc_bga")
  rcd=nf90_wrp(nf90_def_var(nc_id,'mpc_dry_air',nf90_float,mpc_dry_air_id),sbr_nm//": dv mpc_dry_air")
  rcd=nf90_wrp(nf90_def_var(nc_id,'mpc_mst_air',nf90_float,mpc_mst_air_id),sbr_nm//": dv mpc_mst_air")
  rcd=nf90_wrp(nf90_def_var(nc_id,'mpl_CFC11',nf90_float,lev_dmn_id,mpl_CFC11_id),sbr_nm//": dv mpl_CFC11")
  rcd=nf90_wrp(nf90_def_var(nc_id,'mpl_CFC12',nf90_float,lev_dmn_id,mpl_CFC12_id),sbr_nm//": dv mpl_CFC12")
  rcd=nf90_wrp(nf90_def_var(nc_id,'mpl_CH4',nf90_float,lev_dmn_id,mpl_CH4_id),sbr_nm//": dv mpl_CH4")
  rcd=nf90_wrp(nf90_def_var(nc_id,'mpl_CO2',nf90_float,lev_dmn_id,mpl_CO2_id),sbr_nm//": dv mpl_CO2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'mpl_CWP',nf90_float,lev_dmn_id,mpl_CWP_id),sbr_nm//": dv mpl_CWP")
  rcd=nf90_wrp(nf90_def_var(nc_id,'mpl_H2O',nf90_float,lev_dmn_id,mpl_H2O_id),sbr_nm//": dv mpl_H2O")
  rcd=nf90_wrp(nf90_def_var(nc_id,'mpl_H2OH2O',nf90_float,lev_dmn_id,mpl_H2OH2O_id),sbr_nm//": dv mpl_H2OH2O")
  rcd=nf90_wrp(nf90_def_var(nc_id,'mpl_IWP',nf90_float,lev_dmn_id,mpl_IWP_id),sbr_nm//": dv mpl_IWP")
  rcd=nf90_wrp(nf90_def_var(nc_id,'mpl_LWP',nf90_float,lev_dmn_id,mpl_LWP_id),sbr_nm//": dv mpl_LWP")
  rcd=nf90_wrp(nf90_def_var(nc_id,'mpl_N2',nf90_float,lev_dmn_id,mpl_N2_id),sbr_nm//": dv mpl_N2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'mpl_N2O',nf90_float,lev_dmn_id,mpl_N2O_id),sbr_nm//": dv mpl_N2O")
  rcd=nf90_wrp(nf90_def_var(nc_id,'mpl_NO2',nf90_float,lev_dmn_id,mpl_NO2_id),sbr_nm//": dv mpl_NO2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'mpl_O2',nf90_float,lev_dmn_id,mpl_O2_id),sbr_nm//": dv mpl_O2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'mpl_O2O2',nf90_float,lev_dmn_id,mpl_O2O2_id),sbr_nm//": dv mpl_O2O2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'mpl_O3',nf90_float,lev_dmn_id,mpl_O3_id),sbr_nm//": dv mpl_O3")
  rcd=nf90_wrp(nf90_def_var(nc_id,'mpl_OH',nf90_float,lev_dmn_id,mpl_OH_id),sbr_nm//": dv mpl_OH")
  rcd=nf90_wrp(nf90_def_var(nc_id,'mpl_aer',nf90_float,lev_dmn_id,mpl_aer_id),sbr_nm//": dv mpl_aer")
  rcd=nf90_wrp(nf90_def_var(nc_id,'mpl_bga',nf90_float,lev_dmn_id,mpl_bga_id),sbr_nm//": dv mpl_bga")
  rcd=nf90_wrp(nf90_def_var(nc_id,'mpl_dry_air',nf90_float,lev_dmn_id,mpl_dry_air_id),sbr_nm//": dv mpl_dry_air")
  rcd=nf90_wrp(nf90_def_var(nc_id,'mpl_mst_air',nf90_float,lev_dmn_id,mpl_mst_air_id),sbr_nm//": dv mpl_mst_air")
  rcd=nf90_wrp(nf90_def_var(nc_id,'npc_CFC11',nf90_float,npc_CFC11_id),sbr_nm//": dv npc_CFC11")
  rcd=nf90_wrp(nf90_def_var(nc_id,'npc_CFC12',nf90_float,npc_CFC12_id),sbr_nm//": dv npc_CFC12")
  rcd=nf90_wrp(nf90_def_var(nc_id,'npc_CH4',nf90_float,npc_CH4_id),sbr_nm//": dv npc_CH4")
  rcd=nf90_wrp(nf90_def_var(nc_id,'npc_CO2',nf90_float,npc_CO2_id),sbr_nm//": dv npc_CO2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'npc_H2O',nf90_float,npc_H2O_id),sbr_nm//": dv npc_H2O")
  rcd=nf90_wrp(nf90_def_var(nc_id,'npc_H2OH2O',nf90_float,npc_H2OH2O_id),sbr_nm//": dv npc_H2OH2O")
  rcd=nf90_wrp(nf90_def_var(nc_id,'npc_N2',nf90_float,npc_N2_id),sbr_nm//": dv npc_N2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'npc_N2O',nf90_float,npc_N2O_id),sbr_nm//": dv npc_N2O")
  rcd=nf90_wrp(nf90_def_var(nc_id,'npc_NO2',nf90_float,npc_NO2_id),sbr_nm//": dv npc_NO2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'npc_O2',nf90_float,npc_O2_id),sbr_nm//": dv npc_O2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'npc_O2O2',nf90_float,npc_O2O2_id),sbr_nm//": dv npc_O2O2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'npc_O3',nf90_float,npc_O3_id),sbr_nm//": dv npc_O3")
  rcd=nf90_wrp(nf90_def_var(nc_id,'npc_OH',nf90_float,npc_OH_id),sbr_nm//": dv npc_OH")
  rcd=nf90_wrp(nf90_def_var(nc_id,'npc_dry_air',nf90_float,npc_dry_air_id),sbr_nm//": dv npc_dry_air")
  rcd=nf90_wrp(nf90_def_var(nc_id,'npc_mst_air',nf90_float,npc_mst_air_id),sbr_nm//": dv npc_mst_air")
  rcd=nf90_wrp(nf90_def_var(nc_id,'npl_CFC11',nf90_float,lev_dmn_id,npl_CFC11_id),sbr_nm//": dv npl_CFC11")
  rcd=nf90_wrp(nf90_def_var(nc_id,'npl_CFC12',nf90_float,lev_dmn_id,npl_CFC12_id),sbr_nm//": dv npl_CFC12")
  rcd=nf90_wrp(nf90_def_var(nc_id,'npl_CH4',nf90_float,lev_dmn_id,npl_CH4_id),sbr_nm//": dv npl_CH4")
  rcd=nf90_wrp(nf90_def_var(nc_id,'npl_CO2',nf90_float,lev_dmn_id,npl_CO2_id),sbr_nm//": dv npl_CO2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'npl_H2O',nf90_float,lev_dmn_id,npl_H2O_id),sbr_nm//": dv npl_H2O")
  rcd=nf90_wrp(nf90_def_var(nc_id,'npl_H2OH2O',nf90_float,lev_dmn_id,npl_H2OH2O_id),sbr_nm//": dv npl_H2OH2O")
  rcd=nf90_wrp(nf90_def_var(nc_id,'npl_N2',nf90_float,lev_dmn_id,npl_N2_id),sbr_nm//": dv npl_N2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'npl_N2O',nf90_float,lev_dmn_id,npl_N2O_id),sbr_nm//": dv npl_N2O")
  rcd=nf90_wrp(nf90_def_var(nc_id,'npl_NO2',nf90_float,lev_dmn_id,npl_NO2_id),sbr_nm//": dv npl_NO2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'npl_O2',nf90_float,lev_dmn_id,npl_O2_id),sbr_nm//": dv npl_O2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'npl_O2O2',nf90_float,lev_dmn_id,npl_O2O2_id),sbr_nm//": dv npl_O2O2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'npl_O3',nf90_float,lev_dmn_id,npl_O3_id),sbr_nm//": dv npl_O3")
  rcd=nf90_wrp(nf90_def_var(nc_id,'npl_OH',nf90_float,lev_dmn_id,npl_OH_id),sbr_nm//": dv npl_OH")
  rcd=nf90_wrp(nf90_def_var(nc_id,'npl_dry_air',nf90_float,lev_dmn_id,npl_dry_air_id),sbr_nm//": dv npl_dry_air")
  rcd=nf90_wrp(nf90_def_var(nc_id,'npl_mst_air',nf90_float,lev_dmn_id,npl_mst_air_id),sbr_nm//": dv npl_mst_air")
  rcd=nf90_wrp(nf90_def_var(nc_id,'nrg_dry',nf90_float,lev_dmn_id,nrg_dry_id),sbr_nm//": dv nrg_dry")
  rcd=nf90_wrp(nf90_def_var(nc_id,'nrg_mst',nf90_float,lev_dmn_id,nrg_mst_id),sbr_nm//": dv nrg_mst")
  rcd=nf90_wrp(nf90_def_var(nc_id,'odxc_obs_aer',nf90_float,odxc_obs_aer_id),sbr_nm//": dv odxc_obs_aer")
  rcd=nf90_wrp(nf90_def_var(nc_id,'odxc_obs_bga',nf90_float,odxc_obs_bga_id),sbr_nm//": dv odxc_obs_bga")
  rcd=nf90_wrp(nf90_def_var(nc_id,'odxl_obs_aer',nf90_float,lev_dmn_id,odxl_obs_aer_id),sbr_nm//": dv odxl_obs_aer")
  rcd=nf90_wrp(nf90_def_var(nc_id,'odxl_obs_bga',nf90_float,lev_dmn_id,odxl_obs_bga_id),sbr_nm//": dv odxl_obs_bga")
  rcd=nf90_wrp(nf90_def_var(nc_id,'oneD_foo',nf90_float,lev_dmn_id,oneD_foo_id),sbr_nm//": dv oneD_foo")
  rcd=nf90_wrp(nf90_def_var(nc_id,'oro',nf90_float,oro_id),sbr_nm//": dv oro")
  rcd=nf90_wrp(nf90_def_var(nc_id,'pcl_DLR',nf90_float,lev_dmn_id,pcl_DLR_id),sbr_nm//": dv pcl_DLR")
  rcd=nf90_wrp(nf90_def_var(nc_id,'pcl_MLR',nf90_float,lev_dmn_id,pcl_MLR_id),sbr_nm//": dv pcl_MLR")
  rcd=nf90_wrp(nf90_def_var(nc_id,'pcl_PLR',nf90_float,lev_dmn_id,pcl_PLR_id),sbr_nm//": dv pcl_PLR")
  rcd=nf90_wrp(nf90_def_var(nc_id,'ppr_CFC11',nf90_float,lev_dmn_id,ppr_CFC11_id),sbr_nm//": dv ppr_CFC11")
  rcd=nf90_wrp(nf90_def_var(nc_id,'ppr_CFC12',nf90_float,lev_dmn_id,ppr_CFC12_id),sbr_nm//": dv ppr_CFC12")
  rcd=nf90_wrp(nf90_def_var(nc_id,'ppr_CH4',nf90_float,lev_dmn_id,ppr_CH4_id),sbr_nm//": dv ppr_CH4")
  rcd=nf90_wrp(nf90_def_var(nc_id,'ppr_CO2',nf90_float,lev_dmn_id,ppr_CO2_id),sbr_nm//": dv ppr_CO2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'ppr_H2O',nf90_float,lev_dmn_id,ppr_H2O_id),sbr_nm//": dv ppr_H2O")
  rcd=nf90_wrp(nf90_def_var(nc_id,'ppr_H2OH2O',nf90_float,lev_dmn_id,ppr_H2OH2O_id),sbr_nm//": dv ppr_H2OH2O")
  rcd=nf90_wrp(nf90_def_var(nc_id,'ppr_N2',nf90_float,lev_dmn_id,ppr_N2_id),sbr_nm//": dv ppr_N2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'ppr_N2O',nf90_float,lev_dmn_id,ppr_N2O_id),sbr_nm//": dv ppr_N2O")
  rcd=nf90_wrp(nf90_def_var(nc_id,'ppr_NO2',nf90_float,lev_dmn_id,ppr_NO2_id),sbr_nm//": dv ppr_NO2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'ppr_O2',nf90_float,lev_dmn_id,ppr_O2_id),sbr_nm//": dv ppr_O2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'ppr_O2O2',nf90_float,lev_dmn_id,ppr_O2O2_id),sbr_nm//": dv ppr_O2O2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'ppr_O3',nf90_float,lev_dmn_id,ppr_O3_id),sbr_nm//": dv ppr_O3")
  rcd=nf90_wrp(nf90_def_var(nc_id,'ppr_OH',nf90_float,lev_dmn_id,ppr_OH_id),sbr_nm//": dv ppr_OH")
  rcd=nf90_wrp(nf90_def_var(nc_id,'ppr_dry_air',nf90_float,lev_dmn_id,ppr_dry_air_id),sbr_nm//": dv ppr_dry_air")
  rcd=nf90_wrp(nf90_def_var(nc_id,'prs',nf90_float,lev_dmn_id,prs_id),sbr_nm//": dv prs")
  rcd=nf90_wrp(nf90_def_var(nc_id,'prs_cld_btm',nf90_float,prs_cld_btm_id),sbr_nm//": dv prs_cld_btm")
  rcd=nf90_wrp(nf90_def_var(nc_id,'prs_cld_mid',nf90_float,prs_cld_mid_id),sbr_nm//": dv prs_cld_mid")
  rcd=nf90_wrp(nf90_def_var(nc_id,'prs_cld_thick',nf90_float,prs_cld_thick_id),sbr_nm//": dv prs_cld_thick")
  rcd=nf90_wrp(nf90_def_var(nc_id,'prs_cld_top',nf90_float,prs_cld_top_id),sbr_nm//": dv prs_cld_top")
  rcd=nf90_wrp(nf90_def_var(nc_id,'prs_dlt',nf90_float,lev_dmn_id,prs_dlt_id),sbr_nm//": dv prs_dlt")
  rcd=nf90_wrp(nf90_def_var(nc_id,'prs_ntf',nf90_float,levp_dmn_id,prs_ntf_id),sbr_nm//": dv prs_ntf")
  rcd=nf90_wrp(nf90_def_var(nc_id,'prs_sfc',nf90_float,prs_sfc_id),sbr_nm//": dv prs_sfc")
  rcd=nf90_wrp(nf90_def_var(nc_id,'q_CFC11',nf90_float,lev_dmn_id,q_CFC11_id),sbr_nm//": dv q_CFC11")
  rcd=nf90_wrp(nf90_def_var(nc_id,'q_CFC12',nf90_float,lev_dmn_id,q_CFC12_id),sbr_nm//": dv q_CFC12")
  rcd=nf90_wrp(nf90_def_var(nc_id,'q_CH4',nf90_float,lev_dmn_id,q_CH4_id),sbr_nm//": dv q_CH4")
  rcd=nf90_wrp(nf90_def_var(nc_id,'q_CO2',nf90_float,lev_dmn_id,q_CO2_id),sbr_nm//": dv q_CO2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'q_H2O',nf90_float,lev_dmn_id,q_H2O_id),sbr_nm//": dv q_H2O")
  rcd=nf90_wrp(nf90_def_var(nc_id,'q_H2OH2O',nf90_float,lev_dmn_id,q_H2OH2O_id),sbr_nm//": dv q_H2OH2O")
  rcd=nf90_wrp(nf90_def_var(nc_id,'q_N2',nf90_float,lev_dmn_id,q_N2_id),sbr_nm//": dv q_N2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'q_N2O',nf90_float,lev_dmn_id,q_N2O_id),sbr_nm//": dv q_N2O")
  rcd=nf90_wrp(nf90_def_var(nc_id,'q_NO2',nf90_float,lev_dmn_id,q_NO2_id),sbr_nm//": dv q_NO2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'q_O2',nf90_float,lev_dmn_id,q_O2_id),sbr_nm//": dv q_O2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'q_O2O2',nf90_float,lev_dmn_id,q_O2O2_id),sbr_nm//": dv q_O2O2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'q_O3',nf90_float,lev_dmn_id,q_O3_id),sbr_nm//": dv q_O3")
  rcd=nf90_wrp(nf90_def_var(nc_id,'q_OH',nf90_float,lev_dmn_id,q_OH_id),sbr_nm//": dv q_OH")
  rcd=nf90_wrp(nf90_def_var(nc_id,'qst_H2O_ice',nf90_float,lev_dmn_id,qst_H2O_ice_id),sbr_nm//": dv qst_H2O_ice")
  rcd=nf90_wrp(nf90_def_var(nc_id,'qst_H2O_lqd',nf90_float,lev_dmn_id,qst_H2O_lqd_id),sbr_nm//": dv qst_H2O_lqd")
  rcd=nf90_wrp(nf90_def_var(nc_id,'r_CFC11',nf90_float,lev_dmn_id,r_CFC11_id),sbr_nm//": dv r_CFC11")
  rcd=nf90_wrp(nf90_def_var(nc_id,'r_CFC12',nf90_float,lev_dmn_id,r_CFC12_id),sbr_nm//": dv r_CFC12")
  rcd=nf90_wrp(nf90_def_var(nc_id,'r_CH4',nf90_float,lev_dmn_id,r_CH4_id),sbr_nm//": dv r_CH4")
  rcd=nf90_wrp(nf90_def_var(nc_id,'r_CO2',nf90_float,lev_dmn_id,r_CO2_id),sbr_nm//": dv r_CO2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'r_H2O',nf90_float,lev_dmn_id,r_H2O_id),sbr_nm//": dv r_H2O")
  rcd=nf90_wrp(nf90_def_var(nc_id,'r_H2OH2O',nf90_float,lev_dmn_id,r_H2OH2O_id),sbr_nm//": dv r_H2OH2O")
  rcd=nf90_wrp(nf90_def_var(nc_id,'r_N2',nf90_float,lev_dmn_id,r_N2_id),sbr_nm//": dv r_N2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'r_N2O',nf90_float,lev_dmn_id,r_N2O_id),sbr_nm//": dv r_N2O")
  rcd=nf90_wrp(nf90_def_var(nc_id,'r_NO2',nf90_float,lev_dmn_id,r_NO2_id),sbr_nm//": dv r_NO2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'r_O2',nf90_float,lev_dmn_id,r_O2_id),sbr_nm//": dv r_O2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'r_O2O2',nf90_float,lev_dmn_id,r_O2O2_id),sbr_nm//": dv r_O2O2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'r_O3',nf90_float,lev_dmn_id,r_O3_id),sbr_nm//": dv r_O3")
  rcd=nf90_wrp(nf90_def_var(nc_id,'r_OH',nf90_float,lev_dmn_id,r_OH_id),sbr_nm//": dv r_OH")
  rcd=nf90_wrp(nf90_def_var(nc_id,'rds_fct_ice',nf90_float,lev_dmn_id,rds_fct_ice_id),sbr_nm//": dv rds_fct_ice")
  rcd=nf90_wrp(nf90_def_var(nc_id,'rds_fct_lqd',nf90_float,lev_dmn_id,rds_fct_lqd_id),sbr_nm//": dv rds_fct_lqd")
  rcd=nf90_wrp(nf90_def_var(nc_id,'rgh_len',nf90_float,rgh_len_id),sbr_nm//": dv rgh_len")
  rcd=nf90_wrp(nf90_def_var(nc_id,'scl_hgt',nf90_float,lev_dmn_id,scl_hgt_id),sbr_nm//": dv scl_hgt")
  rcd=nf90_wrp(nf90_def_var(nc_id,'sfc_ems',nf90_float,sfc_ems_id),sbr_nm//": dv sfc_ems")
  rcd=nf90_wrp(nf90_def_var(nc_id,'slr_azi_dgr',nf90_float,slr_azi_dgr_id),sbr_nm//": dv slr_azi_dgr")
  rcd=nf90_wrp(nf90_def_var(nc_id,'slr_crd_gmm_dgr',nf90_float,slr_crd_gmm_dgr_id),sbr_nm//": dv slr_crd_gmm_dgr")
  rcd=nf90_wrp(nf90_def_var(nc_id,'slr_cst',nf90_float,slr_cst_id),sbr_nm//": dv slr_cst")
  rcd=nf90_wrp(nf90_def_var(nc_id,'slr_dcl_dgr',nf90_float,slr_dcl_dgr_id),sbr_nm//": dv slr_dcl_dgr")
  rcd=nf90_wrp(nf90_def_var(nc_id,'slr_dmt_dgr',nf90_float,slr_dmt_dgr_id),sbr_nm//": dv slr_dmt_dgr")
  rcd=nf90_wrp(nf90_def_var(nc_id,'slr_dst_au',nf90_float,slr_dst_au_id),sbr_nm//": dv slr_dst_au")
  rcd=nf90_wrp(nf90_def_var(nc_id,'slr_elv_dgr',nf90_float,slr_elv_dgr_id),sbr_nm//": dv slr_elv_dgr")
  rcd=nf90_wrp(nf90_def_var(nc_id,'slr_flx_TOA',nf90_float,slr_flx_TOA_id),sbr_nm//": dv slr_flx_TOA")
  rcd=nf90_wrp(nf90_def_var(nc_id,'slr_flx_nrm_TOA',nf90_float,slr_flx_nrm_TOA_id),sbr_nm//": dv slr_flx_nrm_TOA")
  rcd=nf90_wrp(nf90_def_var(nc_id,'slr_hr_ngl_dgr',nf90_float,slr_hr_ngl_dgr_id),sbr_nm//": dv slr_hr_ngl_dgr")
  rcd=nf90_wrp(nf90_def_var(nc_id,'slr_rfr_ngl_dgr',nf90_float,slr_rfr_ngl_dgr_id),sbr_nm//": dv slr_rfr_ngl_dgr")
  rcd=nf90_wrp(nf90_def_var(nc_id,'slr_rgt_asc_dgr',nf90_float,slr_rgt_asc_dgr_id),sbr_nm//": dv slr_rgt_asc_dgr")
  rcd=nf90_wrp(nf90_def_var(nc_id,'slr_zen_ngl',nf90_double,slr_zen_ngl_id),sbr_nm//": dv slr_zen_ngl")
  rcd=nf90_wrp(nf90_def_var(nc_id,'slr_zen_ngl_cos',nf90_double,slr_zen_ngl_cos_id),sbr_nm//": dv slr_zen_ngl_cos")
  rcd=nf90_wrp(nf90_def_var(nc_id,'slr_zen_ngl_dgr',nf90_double,slr_zen_ngl_dgr_id),sbr_nm//": dv slr_zen_ngl_dgr")
  rcd=nf90_wrp(nf90_def_var(nc_id,'snow_depth',nf90_float,snow_depth_id),sbr_nm//": dv snow_depth")
  rcd=nf90_wrp(nf90_def_var(nc_id,'spc_heat_mst_air',nf90_float,lev_dmn_id,spc_heat_mst_air_id),sbr_nm//": dv spc_heat_mst_air")
  rcd=nf90_wrp(nf90_def_var(nc_id,'time_lmt',nf90_double,time_lmt_id),sbr_nm//": dv time_lmt")
  rcd=nf90_wrp(nf90_def_var(nc_id,'time_ltst',nf90_double,time_ltst_id),sbr_nm//": dv time_ltst")
  rcd=nf90_wrp(nf90_def_var(nc_id,'time_unix',nf90_double,time_unix_id),sbr_nm//": dv time_unix")
  rcd=nf90_wrp(nf90_def_var(nc_id,'tpt',nf90_float,lev_dmn_id,tpt_id),sbr_nm//": dv tpt")
  rcd=nf90_wrp(nf90_def_var(nc_id,'tpt_cls',nf90_float,lev_dmn_id,tpt_cls_id),sbr_nm//": dv tpt_cls")
  rcd=nf90_wrp(nf90_def_var(nc_id,'tpt_cls_ntf',nf90_float,levp_dmn_id,tpt_cls_ntf_id),sbr_nm//": dv tpt_cls_ntf")
  rcd=nf90_wrp(nf90_def_var(nc_id,'tpt_dwp',nf90_float,lev_dmn_id,tpt_dwp_id),sbr_nm//": dv tpt_dwp")
  rcd=nf90_wrp(nf90_def_var(nc_id,'tpt_ntf',nf90_float,levp_dmn_id,tpt_ntf_id),sbr_nm//": dv tpt_ntf")
  rcd=nf90_wrp(nf90_def_var(nc_id,'tpt_pcl',nf90_float,lev_dmn_id,tpt_pcl_id),sbr_nm//": dv tpt_pcl")
  rcd=nf90_wrp(nf90_def_var(nc_id,'tpt_pcl_sfc',nf90_float,tpt_pcl_sfc_id),sbr_nm//": dv tpt_pcl_sfc")
  rcd=nf90_wrp(nf90_def_var(nc_id,'tpt_ptn',nf90_float,tpt_ptn_id),sbr_nm//": dv tpt_ptn")
  rcd=nf90_wrp(nf90_def_var(nc_id,'tpt_ptn_vrt',nf90_float,tpt_ptn_vrt_id),sbr_nm//": dv tpt_ptn_vrt")
  rcd=nf90_wrp(nf90_def_var(nc_id,'tpt_sfc',nf90_float,tpt_sfc_id),sbr_nm//": dv tpt_sfc")
  rcd=nf90_wrp(nf90_def_var(nc_id,'tpt_skn',nf90_float,tpt_skn_id),sbr_nm//": dv tpt_skn")
  rcd=nf90_wrp(nf90_def_var(nc_id,'tpt_vrt',nf90_float,lev_dmn_id,tpt_vrt_id),sbr_nm//": dv tpt_vrt")
  rcd=nf90_wrp(nf90_def_var(nc_id,'vmr_CFC11',nf90_float,lev_dmn_id,vmr_CFC11_id),sbr_nm//": dv vmr_CFC11")
  rcd=nf90_wrp(nf90_def_var(nc_id,'vmr_CFC12',nf90_float,lev_dmn_id,vmr_CFC12_id),sbr_nm//": dv vmr_CFC12")
  rcd=nf90_wrp(nf90_def_var(nc_id,'vmr_CH4',nf90_float,lev_dmn_id,vmr_CH4_id),sbr_nm//": dv vmr_CH4")
  rcd=nf90_wrp(nf90_def_var(nc_id,'vmr_CO2',nf90_float,lev_dmn_id,vmr_CO2_id),sbr_nm//": dv vmr_CO2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'vmr_H2O',nf90_float,lev_dmn_id,vmr_H2O_id),sbr_nm//": dv vmr_H2O")
  rcd=nf90_wrp(nf90_def_var(nc_id,'vmr_H2OH2O',nf90_float,lev_dmn_id,vmr_H2OH2O_id),sbr_nm//": dv vmr_H2OH2O")
  rcd=nf90_wrp(nf90_def_var(nc_id,'vmr_N2',nf90_float,lev_dmn_id,vmr_N2_id),sbr_nm//": dv vmr_N2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'vmr_N2O',nf90_float,lev_dmn_id,vmr_N2O_id),sbr_nm//": dv vmr_N2O")
  rcd=nf90_wrp(nf90_def_var(nc_id,'vmr_NO2',nf90_float,lev_dmn_id,vmr_NO2_id),sbr_nm//": dv vmr_NO2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'vmr_O2',nf90_float,lev_dmn_id,vmr_O2_id),sbr_nm//": dv vmr_O2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'vmr_O2O2',nf90_float,lev_dmn_id,vmr_O2O2_id),sbr_nm//": dv vmr_O2O2")
  rcd=nf90_wrp(nf90_def_var(nc_id,'vmr_O3',nf90_float,lev_dmn_id,vmr_O3_id),sbr_nm//": dv vmr_O3")
  rcd=nf90_wrp(nf90_def_var(nc_id,'vmr_OH',nf90_float,lev_dmn_id,vmr_OH_id),sbr_nm//": dv vmr_OH")
  rcd=nf90_wrp(nf90_def_var(nc_id,'wvl_obs_aer',nf90_float,wvl_obs_aer_id),sbr_nm//": dv wvl_obs_aer")
  rcd=nf90_wrp(nf90_def_var(nc_id,'wvl_obs_bga',nf90_float,wvl_obs_bga_id),sbr_nm//": dv wvl_obs_bga")
  rcd=nf90_wrp(nf90_def_var(nc_id,'xnt_fac',nf90_float,xnt_fac_id),sbr_nm//": dv xnt_fac")
  ! Wrap
  rcd=nf90_wrp(nf90_def_var(nc_id,'q_H2OH2O_rcp_q_H2O',nf90_float,lev_dmn_id,q_H2OH2O_rcp_q_H2O_id), &
       sbr_nm//": dv q_H2OH2O_rcp_q_H2O")
  rcd=nf90_wrp(nf90_def_var(nc_id,'cnc_O2_npl_O2_clm_frc',nf90_double,lev_dmn_id,cnc_O2_npl_O2_clm_frc_id), &
       sbr_nm//": dv cnc_O2_npl_O2_clm_frc")
  if (flg_snw) then
     rcd=nf90_wrp(nf90_def_var(nc_id,'dns_snw',nf90_float,lev_snw_dmn_id,dns_snw_id),sbr_nm//": dv dns_snw")
     rcd=nf90_wrp(nf90_def_var(nc_id,'dpt_dlt_snw',nf90_float,lev_snw_dmn_id,dpt_dlt_snw_id),sbr_nm//": dv dpt_dlt_snw")
     rcd=nf90_wrp(nf90_def_var(nc_id,'dpt_ntf_snw',nf90_float,levp_snw_dmn_id,dpt_ntf_snw_id),sbr_nm//": dv dpt_ntf_snw")
     rcd=nf90_wrp(nf90_def_var(nc_id,'dpt_snw',nf90_float,lev_snw_dmn_id,dpt_snw_id),sbr_nm//": dv dpt_snw")
     rcd=nf90_wrp(nf90_def_var(nc_id,'foo_snw',nf90_float,lev_snw_dmn_id,foo_snw_id),sbr_nm//": dv foo_snw")
     rcd=nf90_wrp(nf90_def_var(nc_id,'lev_snw',nf90_float,lev_snw_dmn_id,lev_snw_id),sbr_nm//": dv lev_snw")
     rcd=nf90_wrp(nf90_def_var(nc_id,'levp_snw',nf90_float,levp_snw_dmn_id,levp_snw_id),sbr_nm//": dv levp_snw")
     rcd=nf90_wrp(nf90_def_var(nc_id,'mmr_mpr_snw',nf90_float,lev_snw_dmn_id,mmr_mpr_snw_id),sbr_nm//": dv mmr_mpr_snw")
     rcd=nf90_wrp(nf90_def_var(nc_id,'rds_ffc_snw',nf90_float,lev_snw_dmn_id,rds_ffc_snw_id),sbr_nm//": dv rds_ffc_snw")
     rcd=nf90_wrp(nf90_def_var(nc_id,'tpt_ntf_snw',nf90_float,levp_snw_dmn_id,tpt_ntf_snw_id),sbr_nm//": dv tpt_ntf_snw")
     rcd=nf90_wrp(nf90_def_var(nc_id,'tpt_snw',nf90_float,lev_snw_dmn_id,tpt_snw_id),sbr_nm//": dv tpt_snw")
  endif ! !flg_snw
  if (flg_wnd) then
     rcd=nf90_wrp(nf90_def_var(nc_id,'abl_hgt',nf90_float,abl_hgt_id),sbr_nm//": dv abl_hgt")
     rcd=nf90_wrp(nf90_def_var(nc_id,'kbfs',nf90_float,kbfs_id),sbr_nm//": dv kbfs")
     rcd=nf90_wrp(nf90_def_var(nc_id,'khfs',nf90_float,khfs_id),sbr_nm//": dv khfs")
     rcd=nf90_wrp(nf90_def_var(nc_id,'kqfs',nf90_float,kqfs_id),sbr_nm//": dv kqfs")
     rcd=nf90_wrp(nf90_def_var(nc_id,'lhflx',nf90_float,lhflx_id),sbr_nm//": dv lhflx")
     rcd=nf90_wrp(nf90_def_var(nc_id,'obk_len',nf90_float,obk_len_id),sbr_nm//": dv obk_len")
     rcd=nf90_wrp(nf90_def_var(nc_id,'ric_nbr',nf90_float,lev_dmn_id,ric_nbr_id),sbr_nm//": dv ric_nbr")
     rcd=nf90_wrp(nf90_def_var(nc_id,'shflx',nf90_float,shflx_id),sbr_nm//": dv shflx")
     rcd=nf90_wrp(nf90_def_var(nc_id,'tau',nf90_float,tau_id),sbr_nm//": dv tau")
     rcd=nf90_wrp(nf90_def_var(nc_id,'tke',nf90_float,levp_dmn_id,tke_id),sbr_nm//": dv tke")
     rcd=nf90_wrp(nf90_def_var(nc_id,'ustar',nf90_float,ustar_id),sbr_nm//": dv ustar")
     rcd=nf90_wrp(nf90_def_var(nc_id,'wnd_dir',nf90_float,lev_dmn_id,wnd_dir_id),sbr_nm//": dv wnd_dir")
     rcd=nf90_wrp(nf90_def_var(nc_id,'wnd_shr',nf90_float,lev_dmn_id,wnd_shr_id),sbr_nm//": dv wnd_shr")
     rcd=nf90_wrp(nf90_def_var(nc_id,'wnd_spd',nf90_float,lev_dmn_id,wnd_spd_id),sbr_nm//": dv wnd_spd")
     rcd=nf90_wrp(nf90_def_var(nc_id,'wspd10m',nf90_float,wspd10m_id),sbr_nm//": dv wspd10m")
     rcd=nf90_wrp(nf90_def_var(nc_id,'wstar',nf90_float,wstar_id),sbr_nm//": dv wstar")
  end if
  
  ! Add global attributes
  rcd=rcd+nf90_put_att(nc_id,nf90_global,'CVS_Id',CVS_Id)
  rcd=rcd+nf90_put_att(nc_id,nf90_global,'creation_date',lcl_date_time)
  rcd=rcd+nf90_put_att(nc_id,nf90_global,'prg_ID',prg_ID(1:ftn_strlen(prg_ID)))
  rcd=rcd+nf90_put_att(nc_id,nf90_global,'cmd_ln',cmd_ln(1:ftn_strlen(cmd_ln)))
  rcd=rcd+nf90_put_att(nc_id,nf90_global,'prf_sng',prf_sng(1:ftn_strlen(prf_sng)))
  if (flg_snw) then
     rcd=rcd+nf90_put_att(nc_id,nf90_global,'prf_snw_sng',prf_snw_sng(1:ftn_strlen(prf_snw_sng)))
  endif ! !flg_snw
  rcd=rcd+nf90_put_att(nc_id,nf90_global,'gmt_sng',gmt_sng(1:ftn_strlen(gmt_sng)))
  rcd=rcd+nf90_put_att(nc_id,nf90_global,'lmt_sng',lmt_sng(1:ftn_strlen(lmt_sng)))
  rcd=rcd+nf90_put_att(nc_id,nf90_global,'ltst_sng',ltst_sng(1:ftn_strlen(ltst_sng)))
  
  ! Add English text descriptions
  rcd=nf90_wrp(nf90_put_att(nc_id,CO2_vmr_clm_id,'long_name','Carbon Dioxide volume mixing ratio'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,N2O_vmr_clm_id,'long_name','Nitrous Oxide volume mixing ratio'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,CH4_vmr_clm_id,'long_name','Methane volume mixing ratio'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,CFC11_vmr_clm_id,'long_name','CFC11 volume mixing ratio'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,CFC12_vmr_clm_id,'long_name','CFC12 volume mixing ratio'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,RH_ice_id,'long_name','Relative humidity w/r/t ice'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,RH_id,'long_name','Relative humidity'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,PDF_bin_spd_id,'long_name','Wind speed PDF bin wind speed'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,PDF_bin_wgt_id,'long_name','Wind Speed PDF bin weight'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,RH_lqd_id,'long_name','Relative humidity w/r/t liquid'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,alb_sfc_NIR_drc_id,'long_name','NIR reflectance small zenith angles'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,alb_sfc_NIR_dff_id,'long_name','NIR reflectance large zenith angles'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,alb_sfc_id,'long_name','Broadband surface albedo'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,alb_sfc_vsb_drc_id,'long_name','Visible reflectance at small zenith angles'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,alb_sfc_vsb_dff_id,'long_name','Visible reflectance at large zenith angles'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,alt_cld_btm_id,'long_name','Highest interface beneath all clouds in column'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,alt_cld_mid_id,'long_name','Altitude at midpoint of all clouds in column'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,alt_cld_thick_id,'long_name','Thickness of region containing all clouds'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,alt_cld_top_id,'long_name','Lowest interface above all clouds in column'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,alt_dlt_id,'long_name','Layer altitude thickness'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,alt_id,'long_name','Altitude'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,alt_ntf_id,'long_name','Interface altitude'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,CAPE_id,'long_name','Convective available potential energy'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,CINE_id,'long_name','Convective inhibition energy'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,cld_frc_id,'long_name','Cloud fraction'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,cnc_CO2_id,'long_name','CO2 concentration'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,cnc_CH4_id,'long_name','CH4 concentration'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,cnc_N2O_id,'long_name','N2O concentration'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,cnc_CFC11_id,'long_name','CFC11 concentration'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,cnc_CFC12_id,'long_name','CFC12 concentration'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,cnc_H2OH2O_id,'long_name','H2O dimer concentration'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,cnc_H2O_id,'long_name','H2O concentration'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,cnc_N2_id,'long_name','N2 concentration'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,cnc_NO2_id,'long_name','NO2 concentration'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,cnc_O2O2_id,'long_name','O2O2 concentration'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,cnc_O2_cnc_N2_id,'long_name','O2 number concentration times N2 number concentration'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,cnc_O2_cnc_O2_id,'long_name','O2 number concentration squared'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,cnc_O2_id,'long_name','O2 concentration'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,cnc_O2_npl_N2_clm_id,'long_name','Column total O2 number concentration times N2 number path'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,cnc_O2_npl_N2_id,'long_name','O2 number concentration times N2 number path'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,cnc_O2_npl_O2_clm_id,'long_name','Column total O2 number concentration times O2 number path'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,cnc_O2_npl_O2_clm_frc_id,'long_name','Fraction of column total O2-O2 at or above each layer'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,cnc_O2_npl_O2_id,'long_name','O2 number concentration times O2 number path'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,cnc_O3_id,'long_name','O3 concentration'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,cnc_OH_id,'long_name','OH concentration'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,cnc_dry_air_id,'long_name','Dry concentration'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,cnc_mst_air_id,'long_name','Moist air concentration'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,nrg_dry_id,'long_name','Dry static energy'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,nrg_mst_id,'long_name','Moist static energy'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,pcl_DLR_id,'long_name','Dry adiabatic lapse rate'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,dns_CO2_id,'long_name','Density of CO2'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,dns_CH4_id,'long_name','Density of CH4'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,dns_N2O_id,'long_name','Density of N2O'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,dns_CFC11_id,'long_name','Density of CFC11'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,dns_CFC12_id,'long_name','Density of CFC12'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,dns_H2OH2O_id,'long_name','Density of H20H2O'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,dns_H2O_id,'long_name','Density of H2O'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,dns_N2_id,'long_name','Density of N2'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,dns_NO2_id,'long_name','Density of NO2'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,dns_O2O2_id,'long_name','Density of O2-O2'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,dns_O2_dns_N2_id,'long_name','O2 mass concentration times N2 mass concentration'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,dns_O2_dns_O2_id,'long_name','O2 mass concentration squared'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,dns_O2_id,'long_name','Density of O2'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,dns_O2_mpl_N2_clm_id,'long_name','Column total O2 mass concentration times N2 mass path'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,dns_O2_mpl_N2_id,'long_name','O2 mass concentration times N2 mass path'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,dns_O2_mpl_O2_clm_id,'long_name','Column total O2 mass concentration times O2 mass path'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,dns_O2_mpl_O2_id,'long_name','O2 mass concentration times O2 mass path'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,dns_O3_id,'long_name','Density of O3'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,dns_OH_id,'long_name','Density of OH'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,dns_aer_id,'long_name','Aerosol density'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,dns_bga_id,'long_name','Background aerosol density'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,dns_dry_air_id,'long_name','Density of dry air'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,dns_mst_air_id,'long_name','Density of moist air'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,eqn_time_sec_id,'long_name','Equation of time'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,ext_cff_mss_aer_id,'long_name','Aerosol mass extinction coefficient'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,ext_cff_mss_bga_id,'long_name','Background aerosol mass extinction coefficient'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,frc_ice_id,'long_name','Fraction of condensate that is ice'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,frc_ice_ttl_id,'long_name','Fraction of column condensate that is ice'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,frc_str_zen_ngl_sfc_id,'long_name','Surface fraction of strong zenith angle dependence'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,gas_cst_mst_air_id,'long_name','Specific gas constant for moist air'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,gmt_day_id,'long_name','GMT day'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,gmt_doy_id,'long_name','GMT day-of-year'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,gmt_hr_id,'long_name','GMT hour'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,gmt_mnt_id,'long_name','GMT minute'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,gmt_mth_id,'long_name','GMT month'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,gmt_sec_id,'long_name','GMT second'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,gmt_ydy_id,'long_name','GMT year-day'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,gmt_yr_id,'long_name','GMT year'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,grv_id,'long_name','Gravity'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,oro_id,'long_name','Orography flag'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,lat_cos_id,'long_name','Cosine of latitude'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,lat_dgr_id,'long_name','Latitude (degrees)'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,lat_id,'long_name','Latitude (radians)'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,lev_LCL_id,'long_name','Lifting condensation level'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,lev_LFC_id,'long_name','Level of free convection'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,lcl_time_hr_id,'long_name','Local day hour'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,lcl_yr_day_id,'long_name','Day of year in local time'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,lev_id,'long_name','Layer pressure'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,levp_id,'long_name','Interface pressure'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,lmt_day_id,'long_name','LMT day'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,lmt_hr_id,'long_name','LMT hour'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,lmt_mnt_id,'long_name','LMT minute'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,lmt_mth_id,'long_name','LMT month'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,lmt_sec_id,'long_name','LMT second'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,lmt_ydy_id,'long_name','LMT year-day'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,lmt_yr_id,'long_name','LMT year'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,lev_LNB_id,'long_name','Level of neutral buoyancy'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,lon_dgr_id,'long_name','Longitude'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,lon_id,'long_name','Longitude'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,lon_sec_id,'long_name','Longitude'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,ltst_day_id,'long_name','LTST day'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,ltst_doy_id,'long_name','LTST day-of-year'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,ltst_hr_id,'long_name','LTST hour'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,ltst_mnt_id,'long_name','LTST minute'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,ltst_mth_id,'long_name','LTST month'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,ltst_sec_id,'long_name','LTST second'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,ltst_ydy_id,'long_name','LTST year-day'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,ltst_yr_id,'long_name','LTST year'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,pcl_MLR_id,'long_name','Moist adiabatic lapse rate'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mmw_mst_air_id,'long_name','Mean molecular weight of moist air'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpc_CO2_id,'long_name','Mass path of CO2 in column'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpc_CH4_id,'long_name','Mass path of CH4 in column'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpc_N2O_id,'long_name','Mass path of N2O in column'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpc_CFC11_id,'long_name','Mass path of CFC11 in column'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpc_CFC12_id,'long_name','Mass path of CFC12 in column'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpc_CWP_id,'long_name','Total column Condensed Water Path'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpc_H2OH2O_id,'long_name','Mass path of H2O dimer in column'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpc_H2O_id,'long_name','Mass path of H2O in column'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpc_IWP_id,'long_name','Total column Ice Water Path'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpc_LWP_id,'long_name','Total column Liquid Water Path'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpc_N2_id,'long_name','Mass path of N2 in column'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpc_NO2_id,'long_name','Mass path of NO2 in column'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpc_O2O2_id,'long_name','Mass path of O2-O2 in column'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpc_O2_id,'long_name','Mass path of O2 in column'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpc_O3_DU_id,'long_name','Mass path of O3 in column'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpc_O3_id,'long_name','Mass path of O3 in column'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpc_OH_id,'long_name','Mass path of OH in column'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpc_aer_id,'long_name','Total column mass path of aerosol'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpc_bga_id,'long_name','Total column mass path of background aerosol'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpc_dry_air_id,'long_name','Mass path of dry air in column'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpc_mst_air_id,'long_name','Mass path of moist air in column'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpl_CO2_id,'long_name','Mass path of CO2 in layer'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpl_CH4_id,'long_name','Mass path of CH4 in layer'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpl_N2O_id,'long_name','Mass path of N2O in layer'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpl_CFC11_id,'long_name','Mass path of CFC11 in layer'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpl_CFC12_id,'long_name','Mass path of CFC12 in layer'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpl_CWP_id,'long_name','Layer Condensed Water Path'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpl_H2OH2O_id,'long_name','Mass path of H2O dimer in layer'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpl_H2O_id,'long_name','Mass path of H2O in layer'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpl_IWP_id,'long_name','Layer Ice Water Path'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpl_LWP_id,'long_name','Layer Liquid Water Path'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpl_N2_id,'long_name','Mass path of N2 in layer'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpl_NO2_id,'long_name','Mass path of NO2 in layer'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpl_O2O2_id,'long_name','Mass path of O2-O2 in layer'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpl_O2_id,'long_name','Mass path of O2 in layer'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpl_O3_id,'long_name','Mass path of O3 in layer'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpl_OH_id,'long_name','Mass path of OH in layer'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpl_aer_id,'long_name','Layer mass path of aerosol'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpl_bga_id,'long_name','Layer mass path of aerosol'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpl_dry_air_id,'long_name','Mass path of dry air in layer'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpl_mst_air_id,'long_name','Mass path of moist air in layer'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npc_CO2_id,'long_name','Column number path of CO2'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npc_CH4_id,'long_name','Column number path of CH4'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npc_N2O_id,'long_name','Column number path of N2O'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npc_CFC11_id,'long_name','Column number path of CFC11'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npc_CFC12_id,'long_name','Column number path of CFC12'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npc_H2OH2O_id,'long_name','Column number path of H2O dimer'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npc_H2O_id,'long_name','Column number path of H2O'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npc_N2_id,'long_name','Column number path of O2'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npc_NO2_id,'long_name','Column number path of NO2'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npc_O2O2_id,'long_name','Column number path of O2O2'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npc_O2_id,'long_name','Column number path of O2'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npc_O3_id,'long_name','Column number path of O3'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npc_OH_id,'long_name','Column number path of OH'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npc_dry_air_id,'long_name','Column number path of dry air'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npc_mst_air_id,'long_name','Column number path of moist air'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npl_CO2_id,'long_name','Number path of CO2 in layer'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npl_CH4_id,'long_name','Number path of CH4 in layer'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npl_N2O_id,'long_name','Number path of N2O in layer'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npl_CFC11_id,'long_name','Number path of CFC11 in layer'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npl_CFC12_id,'long_name','Number path of CFC12 in layer'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npl_H2OH2O_id,'long_name','Number path of H2O dimer in layer'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npl_H2O_id,'long_name','Number path of H2O in layer'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npl_N2_id,'long_name','Number path of N2 in layer'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npl_NO2_id,'long_name','Number path of NO2 in layer'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npl_O2O2_id,'long_name','Number path of O2-O2 in layer'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npl_O2_id,'long_name','Number path of O2 in layer'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npl_O3_id,'long_name','Number path of O3 in layer'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npl_OH_id,'long_name','Number path of OH in layer'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npl_dry_air_id,'long_name','Number path of dry air in layer'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npl_mst_air_id,'long_name','Number path of moist air in layer'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,odxc_obs_aer_id,'long_name','Column aerosol extinction optical depth'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,odxc_obs_bga_id,'long_name','Column background aerosol extinction optical depth'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,odxl_obs_aer_id,'long_name','Layer aerosol extinction optical depth'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,odxl_obs_bga_id,'long_name','Layer background aerosol extinction optical depth'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,oneD_foo_id,'long_name',''), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,pcl_PLR_id,'long_name','Pseudo adiabatic lapse rate'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,ppr_CO2_id,'long_name','Partial pressure of CO2'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,ppr_CH4_id,'long_name','Partial pressure of CH4'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,ppr_N2O_id,'long_name','Partial pressure of N2O'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,ppr_CFC11_id,'long_name','Partial pressure of CFC11'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,ppr_CFC12_id,'long_name','Partial pressure of CFC12'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,ppr_H2OH2O_id,'long_name','Partial pressure of H2O dimer'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,ppr_H2O_id,'long_name','Partial pressure of H2O'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,ppr_N2_id,'long_name','Partial pressure of N2'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,ppr_NO2_id,'long_name','Partial pressure of NO2'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,ppr_O2O2_id,'long_name','Partial pressure of O2O2'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,ppr_O2_id,'long_name','Partial pressure of O2'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,ppr_O3_id,'long_name','Partial pressure of O3'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,ppr_OH_id,'long_name','Partial pressure of OH'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,ppr_dry_air_id,'long_name','Partial pressure of dry air'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,prs_cld_btm_id,'long_name','Highest interface beneath all clouds in column'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,prs_cld_mid_id,'long_name','Pressure at midpoint of all clouds in column'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,prs_cld_thick_id,'long_name','Thickness of region containing all clouds'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,prs_cld_top_id,'long_name','Lowest interface above all clouds in column'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,prs_dlt_id,'long_name','Layer pressure thickness'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,prs_id,'long_name','Pressure'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,prs_ntf_id,'long_name','Interface pressure'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,prs_sfc_id,'long_name','Surface pressure'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,q_CO2_id,'long_name','Mass mixing ratio of CO2'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,q_CH4_id,'long_name','Mass mixing ratio of CH4'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,q_N2O_id,'long_name','Mass mixing ratio of N2O'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,q_CFC11_id,'long_name','Mass mixing ratio of CFC11'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,q_CFC12_id,'long_name','Mass mixing ratio of CFC12'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,q_H2OH2O_id,'long_name','Water vapor dimer mass mixing ratio'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,q_H2OH2O_rcp_q_H2O_id,'long_name','Ratio of dimer mmr to monomer mmr'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,q_H2O_id,'long_name','Water vapor mass mixing ratio'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,q_N2_id,'long_name','Mass mixing ratio of N2'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,q_NO2_id,'long_name','Mass mixing ratio of NO2'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,q_O2O2_id,'long_name','Ozone mass mixing ratio'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,q_O2_id,'long_name','Mass mixing ratio of O2'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,q_O3_id,'long_name','Ozone mass mixing ratio'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,q_OH_id,'long_name','Mass mixing ratio of OH'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,qst_H2O_ice_id,'long_name','Saturation mixing ratio w/r/t ice'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,qst_H2O_lqd_id,'long_name','Saturation mixing ratio w/r/t liquid'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,r_CO2_id,'long_name','Dry-mass mixing ratio (r) of CO2'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,r_CH4_id,'long_name','Dry-mass mixing ratio (r) of CH4'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,r_N2O_id,'long_name','Dry-mass mixing ratio (r) of N2O'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,r_CFC11_id,'long_name','Dry-mass mixing ratio (r) of CFC11'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,r_CFC12_id,'long_name','Dry-mass mixing ratio (r) of CFC12'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,r_H2OH2O_id,'long_name','Dry-mass mixing ratio (r) of H2O dimer'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,r_H2O_id,'long_name','Dry-mass mixing ratio (r) of H2O'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,r_N2_id,'long_name','Dry-mass mixing ratio (r) of N2'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,r_NO2_id,'long_name','Dry-mass mixing ratio (r) of NO2'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,r_O2O2_id,'long_name','Dry-mass mixing ratio (r) of O2O2'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,r_O2_id,'long_name','Dry-mass mixing ratio (r) of O2'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,r_O3_id,'long_name','Dry-mass mixing ratio (r) of O3'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,r_OH_id,'long_name','Dry-mass mixing ratio (r) of OH'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,rds_fct_ice_id,'long_name','Effective radius of ice crystals'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,rds_fct_lqd_id,'long_name','Effective radius of liquid droplets'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,rgh_len_id,'long_name','Aerodynamic roughness length'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,scl_hgt_id,'long_name','Local scale height'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,sfc_ems_id,'long_name','Surface emissivity'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,slr_azi_dgr_id,'long_name','Solar azimuth angle'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,slr_crd_gmm_dgr_id,'long_name','Solar coordinate gamma'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,slr_cst_id,'long_name','Solar constant'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,slr_dcl_dgr_id,'long_name','Solar declination'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,slr_dmt_dgr_id,'long_name','Diameter of solar disc'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,slr_dst_au_id,'long_name','Earth-sun distance'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,slr_elv_dgr_id,'long_name','Solar elevation'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,slr_flx_TOA_id,'long_name','Solar flux at TOA'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,slr_flx_nrm_TOA_id,'long_name','Solar constant corrected for orbital position'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,slr_hr_ngl_dgr_id,'long_name','Solar hour angle'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,slr_rfr_ngl_dgr_id,'long_name','Solar refraction angle'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,slr_rgt_asc_dgr_id,'long_name','Solar right ascension'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,slr_zen_ngl_cos_id,'long_name','Cosine solar zenith angle'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,slr_zen_ngl_dgr_id,'long_name','Solar zenith angle in degrees'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,slr_zen_ngl_id,'long_name','Solar zenith angle'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,snow_depth_id,'long_name','Snow depth'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,spc_heat_mst_air_id,'long_name','Specific heat at constant pressure of moist air'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,time_lmt_id,'long_name','Seconds between 1969 and LMT of simulation'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,time_ltst_id,'long_name','Seconds between 1969 and LTST of simulation'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,time_unix_id,'long_name','Seconds between 1969 and GMT of simulation'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,tpt_cls_id,'long_name','Layer temperature (Celsius)'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,tpt_cls_ntf_id,'long_name','Interface temperature (Celsius)'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,tpt_id,'long_name','Layer Temperature'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,tpt_ntf_id,'long_name','Interface temperature'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,tpt_pcl_id,'long_name','Lifted parcel temperature'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,tpt_dwp_id,'long_name','Dewpoint temperature'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,tpt_ptn_id,'long_name','Potential temperature'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,tpt_ptn_vrt_id,'long_name','Virtual potential temperature'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,tpt_sfc_id,'long_name','Temperature of air in contact with surface'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,tpt_pcl_sfc_id,'long_name','Air parcel surface temperature'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,tpt_skn_id,'long_name','Temperature of surface'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,tpt_vrt_id,'long_name','Virtual temperature'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,vmr_CO2_id,'long_name','Volume mixing ratio of CO2'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,vmr_CH4_id,'long_name','Volume mixing ratio of CH4'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,vmr_N2O_id,'long_name','Volume mixing ratio of N2O'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,vmr_CFC11_id,'long_name','Volume mixing ratio of CFC11'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,vmr_CFC12_id,'long_name','Volume mixing ratio of CFC12'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,vmr_H2OH2O_id,'long_name','Volume mixing ratio of H2O dimer'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,vmr_H2O_id,'long_name','Volume mixing ratio of H2O'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,vmr_N2_id,'long_name','Volume mixing ratio of N2'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,vmr_NO2_id,'long_name','Volume mixing ratio of NO2'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,vmr_O2O2_id,'long_name','Volume mixing ratio of O2O2'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,vmr_O2_id,'long_name','Volume mixing ratio of O2'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,vmr_O3_id,'long_name','Volume mixing ratio of O3'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,vmr_OH_id,'long_name','Volume mixing ratio of OH'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,wvl_obs_aer_id,'long_name','Wavelength of aerosol optical depth specification'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,wvl_obs_bga_id,'long_name','Wavelength of background aerosol optical depth specification'), &
       sbr_nm//": pa long_name in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,xnt_fac_id,'long_name','Eccentricity factor'), &
       sbr_nm//": pa long_name in "//__FILE__)
  if (flg_snw) then
     rcd=nf90_wrp(nf90_put_att(nc_id,lev_snw_id,'long_name','Snow depth'), &
          sbr_nm//": pa long_name in "//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,levp_snw_id,'long_name','Snow depth interfaces'), &
          sbr_nm//": pa long_name in "//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,dns_snw_id,'long_name','Snow density'), &
          sbr_nm//": pa long_name in "//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,dpt_dlt_snw_id,'long_name','Snow layer thickness'), &
          sbr_nm//": pa long_name in "//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,dpt_ntf_snw_id,'long_name','Snow depth interfaces'), &
          sbr_nm//": pa long_name in "//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,dpt_snw_id,'long_name','Snow depth'), &
          sbr_nm//": pa long_name in "//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,mmr_mpr_snw_id,'long_name','Mass mixing ratio of impurities in snow'), &
          sbr_nm//": pa long_name in "//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,rds_ffc_snw_id,'long_name','Snow effective radius'), &
          sbr_nm//": pa long_name in "//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,tpt_snw_id,'long_name','Snow temperature'), &
          sbr_nm//": pa long_name in "//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,tpt_ntf_snw_id,'long_name','Snow temperature interfaces'), &
          sbr_nm//": pa long_name in "//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,foo_snw_id,'long_name','foo snow'), &
          sbr_nm//": pa long_name in "//__FILE__)
  endif ! !flg_snw
  if (flg_wnd) then
     rcd=nf90_wrp(nf90_put_att(nc_id,abl_hgt_id,'long_name','Boundary layer height'), &
          sbr_nm//": pa long_name in "//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,kbfs_id,'long_name','Surface kinematic buoyancy flux'), &
          sbr_nm//": pa long_name in "//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,khfs_id,'long_name','Surface kinematic heat flux'), &
          sbr_nm//": pa long_name in "//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,kqfs_id,'long_name','Kinematic moisture flux'), &
          sbr_nm//": pa long_name in "//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,lhflx_id,'long_name','Latent heat flux'), &
          sbr_nm//": pa long_name in "//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,obk_len_id,'long_name','Obukhov length'), &
          sbr_nm//": pa long_name in "//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,ric_nbr_id,'long_name','Richardson number'), &
          sbr_nm//": pa long_name in "//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,shflx_id,'long_name','Sensible heat flux'), &
          sbr_nm//": pa long_name in "//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,tau_id,'long_name','Surface stress'), &
          sbr_nm//": pa long_name in "//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,tke_id,'long_name','Turbulence kinetic energy'), &
          sbr_nm//": pa long_name in "//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,ustar_id,'long_name','Friction velocity'), &
          sbr_nm//": pa long_name in "//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,wnd_dir_id,'long_name','Wind direction'), &
          sbr_nm//": pa long_name in "//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,wnd_shr_id,'long_name','Wind shear'), &
          sbr_nm//": pa long_name in "//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,wnd_spd_id,'long_name','Wind speed'), &
          sbr_nm//": pa long_name in "//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,wspd10m_id,'long_name','10-m surface wind speed'), &
          sbr_nm//": pa long_name in "//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,wstar_id,'long_name','Convective velocity scale'), &
          sbr_nm//": pa long_name in "//__FILE__)
  end if
  
  ! Add units
  rcd=nf90_wrp(nf90_put_att(nc_id,CO2_vmr_clm_id,'units','fraction'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,N2O_vmr_clm_id,'units','fraction'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,CH4_vmr_clm_id,'units','fraction'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,CFC11_vmr_clm_id,'units','fraction'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,CFC12_vmr_clm_id,'units','fraction'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,RH_ice_id,'units','fraction'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,RH_id,'units','fraction'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,PDF_bin_spd_id,'units','m s-1'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,PDF_bin_wgt_id,'units','fraction'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,RH_lqd_id,'units','fraction'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,alb_sfc_NIR_drc_id,'units','fraction'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,alb_sfc_NIR_dff_id,'units','fraction'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,alb_sfc_id,'units','fraction'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,alb_sfc_vsb_drc_id,'units','fraction'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,alb_sfc_vsb_dff_id,'units','fraction'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,alt_cld_btm_id,'units','meter'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,alt_cld_mid_id,'units','meter'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,alt_cld_thick_id,'units','meter'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,alt_cld_top_id,'units','meter'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,alt_dlt_id,'units','meter'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,alt_id,'units','meter'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,alt_ntf_id,'units','meter'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,CAPE_id,'units','joule kilogram-1'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,CINE_id,'units','joule kilogram-1'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,cld_frc_id,'units','fraction'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,cnc_CO2_id,'units','molecule meter-3'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,cnc_CH4_id,'units','molecule meter-3'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,cnc_N2O_id,'units','molecule meter-3'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,cnc_CFC11_id,'units','molecule meter-3'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,cnc_CFC12_id,'units','molecule meter-3'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,cnc_H2OH2O_id,'units','molecule meter-3'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,cnc_H2O_id,'units','molecule meter-3'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,cnc_N2_id,'units','molecule meter-3'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,cnc_NO2_id,'units','molecule meter-3'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,cnc_O2O2_id,'units','molecule meter-3'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,cnc_O2_cnc_N2_id,'units','molecule2 meter-6'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,cnc_O2_cnc_O2_id,'units','molecule2 meter-6'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,cnc_O2_id,'units','molecule meter-3'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,cnc_O2_npl_N2_clm_id,'units','molecule2 meter-5'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,cnc_O2_npl_N2_id,'units','molecule2 meter-5'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,cnc_O2_npl_O2_clm_id,'units','molecule2 meter-5'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,cnc_O2_npl_O2_clm_frc_id,'units','fraction'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,cnc_O2_npl_O2_id,'units','molecule2 meter-5'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,cnc_O3_id,'units','molecule meter-3'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,cnc_OH_id,'units','molecule meter-3'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,cnc_dry_air_id,'units','molecule meter-3'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,cnc_mst_air_id,'units','molecule meter-3'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,nrg_dry_id,'units','joule kilogram-1'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,nrg_mst_id,'units','joule kilogram-1'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,pcl_DLR_id,'units','kelvin meter-1'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,dns_CO2_id,'units','kilogram meter-3'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,dns_CH4_id,'units','kilogram meter-3'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,dns_N2O_id,'units','kilogram meter-3'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,dns_CFC11_id,'units','kilogram meter-3'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,dns_CFC12_id,'units','kilogram meter-3'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,dns_H2OH2O_id,'units','kilogram meter-3'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,dns_H2O_id,'units','kilogram meter-3'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,dns_N2_id,'units','kilogram meter-3'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,dns_NO2_id,'units','kilogram meter-3'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,dns_O2O2_id,'units','kilogram meter-3'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,dns_O2_dns_N2_id,'units','kilogram2 meter-6'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,dns_O2_dns_O2_id,'units','kilogram2 meter-6'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,dns_O2_id,'units','kilogram meter-3'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,dns_O2_mpl_N2_clm_id,'units','kilogram2 meter-5'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,dns_O2_mpl_N2_id,'units','kilogram2 meter-5'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,dns_O2_mpl_O2_clm_id,'units','kilogram2 meter-5'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,dns_O2_mpl_O2_id,'units','kilogram2 meter-5'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,dns_O3_id,'units','kilogram meter-3'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,dns_OH_id,'units','kilogram meter-3'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,dns_aer_id,'units','kilogram meter-3'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,dns_bga_id,'units','kilogram meter-3'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,dns_dry_air_id,'units','kilogram meter-3'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,dns_mst_air_id,'units','kilogram meter-3'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,eqn_time_sec_id,'units','second'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,ext_cff_mss_aer_id,'units','meter2 kilogram-1'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,ext_cff_mss_bga_id,'units','meter2 kilogram-1'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,frc_ice_id,'units','fraction'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,frc_ice_ttl_id,'units','fraction'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,frc_str_zen_ngl_sfc_id,'units','fraction'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,gas_cst_mst_air_id,'units','joule kilogram-1 kelvin-1'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,gmt_day_id,'units','day'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,gmt_doy_id,'units','day'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,gmt_hr_id,'units','hour'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,gmt_mnt_id,'units','minute'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,gmt_mth_id,'units','month'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,gmt_sec_id,'units','second'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,gmt_ydy_id,'units','day'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,gmt_yr_id,'units','year'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,grv_id,'units','meter second-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,oro_id,'units','flag'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,lat_cos_id,'units','fraction'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,lat_dgr_id,'units','degree'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,lat_id,'units','radian'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,lev_LCL_id,'units','meter'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,lcl_time_hr_id,'units','hour'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,lcl_yr_day_id,'units','day'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,lev_id,'units','pascal'),sbr_nm//": pa units in "//__FILE__)
  if (flg_snw) then
     rcd=nf90_wrp(nf90_put_att(nc_id,lev_snw_id,'units','meter'),sbr_nm//": pa units in "//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,levp_snw_id,'units','meter'),sbr_nm//": pa units in "//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,dns_snw_id,'units','kilogram meter-3'),sbr_nm//": pa units in "//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,dpt_dlt_snw_id,'units','meter'),sbr_nm//": pa units in "//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,dpt_ntf_snw_id,'units','meter'),sbr_nm//": pa units in "//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,dpt_snw_id,'units','meter'),sbr_nm//": pa units in "//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,mmr_mpr_snw_id,'units','kilogram kilogram-1'),sbr_nm//": pa units in "//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,rds_ffc_snw_id,'units','meter'),sbr_nm//": pa units in "//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,tpt_snw_id,'units','kelvin'),sbr_nm//": pa units in "//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,tpt_ntf_snw_id,'units','kelvin'),sbr_nm//": pa units in "//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,foo_snw_id,'units','meter'),sbr_nm//": pa units in "//__FILE__)
  endif ! !flg_snw
  rcd=nf90_wrp(nf90_put_att(nc_id,levp_id,'units','pascal'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,lev_LFC_id, 'units','meter'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,lmt_day_id,'units','day'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,lmt_hr_id,'units','hour'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,lmt_mnt_id,'units','minute'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,lmt_mth_id,'units','month'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,lmt_sec_id,'units','second'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,lmt_ydy_id,'units','day'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,lmt_yr_id,'units','year'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,lev_LNB_id, 'units','meter'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,lon_dgr_id,'units','degree'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,lon_id,'units','radian'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,lon_sec_id,'units','second'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,ltst_day_id,'units','day'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,ltst_doy_id,'units','day'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,ltst_hr_id,'units','hour'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,ltst_mnt_id,'units','minute'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,ltst_mth_id,'units','month'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,ltst_sec_id,'units','second'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,ltst_ydy_id,'units','day'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,ltst_yr_id,'units','year'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mmw_mst_air_id,'units','kilogram mole-1'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,pcl_MLR_id,'units','kelvin meter-1'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpc_CO2_id,'units','kilogram meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpc_CH4_id,'units','kilogram meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpc_N2O_id,'units','kilogram meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpc_CFC11_id,'units','kilogram meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpc_CFC12_id,'units','kilogram meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpc_CWP_id,'units','kilogram meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpc_H2OH2O_id,'units','kilogram meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpc_H2O_id,'units','kilogram meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpc_IWP_id,'units','kilogram meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpc_LWP_id,'units','kilogram meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpc_N2_id,'units','kilogram meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpc_NO2_id,'units','kilogram meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpc_O2O2_id,'units','kilogram meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpc_O2_id,'units','kilogram meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpc_O3_DU_id,'units','Dobson'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpc_O3_id,'units','kilogram meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpc_OH_id,'units','kilogram meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpc_aer_id,'units','kilogram meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpc_bga_id,'units','kilogram meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpc_dry_air_id,'units','kilogram meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpc_mst_air_id,'units','kilogram meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpl_CO2_id,'units','kilogram meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpl_CH4_id,'units','kilogram meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpl_N2O_id,'units','kilogram meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpl_CFC11_id,'units','kilogram meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpl_CFC12_id,'units','kilogram meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpl_CWP_id,'units','kilogram meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpl_H2OH2O_id,'units','kilogram meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpl_H2O_id,'units','kilogram meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpl_IWP_id,'units','kilogram meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpl_LWP_id,'units','kilogram meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpl_N2_id,'units','kilogram meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpl_NO2_id,'units','kilogram meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpl_O2O2_id,'units','kilogram meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpl_O2_id,'units','kilogram meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpl_O3_id,'units','kilogram meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpl_OH_id,'units','kilogram meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpl_aer_id,'units','kilogram meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpl_bga_id,'units','kilogram meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpl_dry_air_id,'units','kilogram meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,mpl_mst_air_id,'units','kilogram meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npc_CO2_id,'units','molecule meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npc_CH4_id,'units','molecule meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npc_N2O_id,'units','molecule meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npc_CFC11_id,'units','molecule meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npc_CFC12_id,'units','molecule meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npc_H2OH2O_id,'units','molecule meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npc_H2O_id,'units','molecule meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npc_N2_id,'units','molecule meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npc_NO2_id,'units','molecule meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npc_O2O2_id,'units','molecule meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npc_O2_id,'units','molecule meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npc_O3_id,'units','molecule meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npc_OH_id,'units','molecule meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npc_dry_air_id,'units','molecule meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npc_mst_air_id,'units','molecule meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npl_CO2_id,'units','molecule meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npl_CH4_id,'units','molecule meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npl_N2O_id,'units','molecule meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npl_CFC11_id,'units','molecule meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npl_CFC12_id,'units','molecule meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npl_H2OH2O_id,'units','molecule meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npl_H2O_id,'units','molecule meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npl_N2_id,'units','molecule meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npl_NO2_id,'units','molecule meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npl_O2O2_id,'units','molecule meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npl_O2_id,'units','molecule meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npl_O3_id,'units','molecule meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npl_OH_id,'units','molecule meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npl_dry_air_id,'units','molecule meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,npl_mst_air_id,'units','molecule meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,odxc_obs_aer_id,'units','fraction'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,odxc_obs_bga_id,'units','fraction'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,odxl_obs_aer_id,'units','fraction'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,odxl_obs_bga_id,'units','fraction'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,oneD_foo_id,'units',''),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,pcl_PLR_id,'units','kelvin meter-1'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,ppr_CO2_id,'units','pascal'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,ppr_CH4_id,'units','pascal'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,ppr_N2O_id,'units','pascal'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,ppr_CFC11_id,'units','pascal'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,ppr_CFC12_id,'units','pascal'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,ppr_H2OH2O_id,'units','pascal'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,ppr_H2O_id,'units','pascal'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,ppr_N2_id,'units','pascal'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,ppr_NO2_id,'units','pascal'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,ppr_O2O2_id,'units','pascal'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,ppr_O2_id,'units','pascal'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,ppr_O3_id,'units','pascal'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,ppr_OH_id,'units','pascal'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,ppr_dry_air_id,'units','pascal'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,prs_cld_btm_id,'units','pascal'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,prs_cld_mid_id,'units','pascal'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,prs_cld_thick_id,'units','meter'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,prs_cld_top_id,'units','pascal'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,prs_dlt_id,'units','pascal'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,prs_id,'units','pascal'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,prs_ntf_id,'units','pascal'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,prs_sfc_id,'units','pascal'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,q_CO2_id,'units','kilogram kilogram-1'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,q_CH4_id,'units','kilogram kilogram-1'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,q_N2O_id,'units','kilogram kilogram-1'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,q_CFC11_id,'units','kilogram kilogram-1'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,q_CFC12_id,'units','kilogram kilogram-1'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,q_H2OH2O_id,'units','kilogram kilogram-1'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,q_H2OH2O_rcp_q_H2O_id,'units','fraction'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,q_H2O_id,'units','kilogram kilogram-1'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,q_N2_id,'units','kilogram kilogram-1'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,q_NO2_id,'units','kilogram kilogram-1'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,q_O2O2_id,'units','kilogram kilogram-1'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,q_O2_id,'units','kilogram kilogram-1'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,q_O3_id,'units','kilogram kilogram-1'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,q_OH_id,'units','kilogram kilogram-1'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,qst_H2O_ice_id,'units','kilogram kilogram-1'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,qst_H2O_lqd_id,'units','kilogram kilogram-1'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,r_CO2_id,'units','kilogram kilogram-1'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,r_CH4_id,'units','kilogram kilogram-1'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,r_N2O_id,'units','kilogram kilogram-1'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,r_CFC11_id,'units','kilogram kilogram-1'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,r_CFC12_id,'units','kilogram kilogram-1'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,r_H2OH2O_id,'units','kilogram kilogram-1'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,r_H2O_id,'units','kilogram kilogram-1'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,r_N2_id,'units','kilogram kilogram-1'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,r_NO2_id,'units','kilogram kilogram-1'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,r_O2O2_id,'units','kilogram kilogram-1'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,r_O2_id,'units','kilogram kilogram-1'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,r_O3_id,'units','kilogram kilogram-1'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,r_OH_id,'units','kilogram kilogram-1'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,rds_fct_ice_id,'units','micron'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,rds_fct_lqd_id,'units','micron'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,rgh_len_id,'units','meter'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,scl_hgt_id,'units','meter'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,sfc_ems_id,'units','fraction'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,slr_azi_dgr_id,'units','degree'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,slr_crd_gmm_dgr_id,'units','degree'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,slr_cst_id,'units','watt meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,slr_dcl_dgr_id,'units','degree'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,slr_dmt_dgr_id,'units','degree'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,slr_dst_au_id,'units','astronomical units'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,slr_elv_dgr_id,'units','degree'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,slr_flx_TOA_id,'units','watt meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,slr_flx_nrm_TOA_id,'units','watt meter-2'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,slr_hr_ngl_dgr_id,'units','degree'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,slr_rfr_ngl_dgr_id,'units','degree'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,slr_rgt_asc_dgr_id,'units','degree'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,slr_zen_ngl_cos_id,'units','fraction'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,slr_zen_ngl_dgr_id,'units','degree'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,slr_zen_ngl_id,'units','radian'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,snow_depth_id,'units','meter'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,spc_heat_mst_air_id,'units','joule kilogram-1 kelvin-1'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,time_lmt_id,'units','second'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,time_ltst_id,'units','second'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,time_unix_id,'units','second'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,tpt_cls_id,'units','celsius'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,tpt_cls_ntf_id,'units','celsius'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,tpt_id,'units','kelvin'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,tpt_ntf_id,'units','kelvin'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,tpt_pcl_id,'units','kelvin'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,tpt_dwp_id,'units','kelvin'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,tpt_sfc_id,'units','kelvin'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,tpt_pcl_sfc_id,'units','kelvin'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,tpt_skn_id,'units','kelvin'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,tpt_vrt_id,'units','kelvin'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,vmr_CO2_id,'units','number number-1'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,vmr_CH4_id,'units','number number-1'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,vmr_N2O_id,'units','number number-1'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,vmr_CFC11_id,'units','number number-1'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,vmr_CFC12_id,'units','number number-1'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,vmr_H2OH2O_id,'units','number number-1'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,vmr_H2O_id,'units','number number-1'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,vmr_N2_id,'units','number number-1'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,vmr_NO2_id,'units','number number-1'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,vmr_O2O2_id,'units','number number-1'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,vmr_O2_id,'units','number number-1'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,vmr_O3_id,'units','number number-1'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,vmr_OH_id,'units','number number-1'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,wvl_obs_aer_id,'units','meter'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,wvl_obs_bga_id,'units','meter'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,xnt_fac_id,'units','fraction'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,tpt_ptn_id,'units','kelvin'),sbr_nm//": pa units in "//__FILE__)
  rcd=nf90_wrp(nf90_put_att(nc_id,tpt_ptn_vrt_id,'units','kelvin'),sbr_nm//": pa units in "//__FILE__)
  if (flg_wnd) then
     rcd=nf90_wrp(nf90_put_att(nc_id,abl_hgt_id,'units','meter'),sbr_nm//": pa units in "//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,kbfs_id,'units','meter2 second-3'),sbr_nm//": pa units in "//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,khfs_id,'units','meter kelvin second-1'),sbr_nm//": pa units in "//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,kqfs_id,'units','meter second-1'),sbr_nm//": pa units in "//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,lhflx_id,'units','watt meter-2'),sbr_nm//": pa units in "//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,obk_len_id,'units','meter'),sbr_nm//": pa units in "//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,ric_nbr_id,'units',''),sbr_nm//": pa units in "//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,shflx_id,'units','watt meter-2'),sbr_nm//": pa units in "//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,tau_id,'units','Newton meter-2'),sbr_nm//": pa units in "//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,tke_id,'units','meter2 second-2'),sbr_nm//": pa units in "//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,ustar_id,'units','meter second-1'),sbr_nm//": pa units in "//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,wnd_dir_id,'units','degree'),sbr_nm//": pa units in "//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,wnd_shr_id,'units','second-1'),sbr_nm//": pa units in "//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,wnd_spd_id,'units','meter second-1'),sbr_nm//": pa units in "//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,wspd10m_id,'units','meter second-1'),sbr_nm//": pa units in "//__FILE__)
     rcd=nf90_wrp(nf90_put_att(nc_id,wstar_id,'units','meter second-1'),sbr_nm//": pa units in "//__FILE__)
  end if
  
  ! Now that all dimensions, variables, and attributes have been defined, make call to end define mode.
  rcd=nf90_wrp(nf90_enddef(nc_id),sbr_nm//": enddef in "//__FILE__) 

  ! Write data
  rcd=nf90_wrp(nf90_put_var(nc_id,CFC11_vmr_clm_id,CFC11_vmr_clm),sbr_nm//": pv CFC11_vmr_clm"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,CFC12_vmr_clm_id,CFC12_vmr_clm),sbr_nm//": pv CFC12_vmr_clm"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,CH4_vmr_clm_id,CH4_vmr_clm),sbr_nm//": pv CH4_vmr_clm"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,CO2_vmr_clm_id,CO2_vmr_clm),sbr_nm//": pv CO2_vmr_clm"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,N2O_vmr_clm_id,N2O_vmr_clm),sbr_nm//": pv N2O_vmr_clm"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,RH_ice_id,RH_ice),sbr_nm//": pv RH_ice"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,RH_id,RH),sbr_nm//": pv RH"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,PDF_bin_spd_id,PDF_bin_spd),sbr_nm//": pv PDF_bin_spd"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,PDF_bin_wgt_id,PDF_bin_wgt),sbr_nm//": pv PDF_bin_wgt"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,RH_lqd_id,RH_lqd),sbr_nm//": pv RH_lqd"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,alb_sfc_NIR_dff_id,alb_sfc_NIR_dff),sbr_nm//": pv alb_sfc_NIR_dff"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,alb_sfc_NIR_drc_id,alb_sfc_NIR_drc),sbr_nm//": pv alb_sfc_NIR_drc"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,alb_sfc_id,alb_sfc),sbr_nm//": pv alb_sfc"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,alb_sfc_vsb_dff_id,alb_sfc_vsb_dff),sbr_nm//": pv alb_sfc_vsb_dff"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,alb_sfc_vsb_drc_id,alb_sfc_vsb_drc),sbr_nm//": pv alb_sfc_vsb_drc"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,alt_cld_btm_id,alt_cld_btm),sbr_nm//": pv alt_cld_btm"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,alt_cld_mid_id,alt_cld_mid),sbr_nm//": pv alt_cld_mid"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,alt_cld_thick_id,alt_cld_thick),sbr_nm//": pv alt_cld_thick"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,alt_cld_top_id,alt_cld_top),sbr_nm//": pv alt_cld_top"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,alt_dlt_id,alt_dlt),sbr_nm//": pv alt_dlt"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,alt_id,alt),sbr_nm//": pv alt"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,alt_ntf_id,alt_ntf),sbr_nm//": pv alt_ntf"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,CAPE_id,CAPE),sbr_nm//": pv CAPE"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,CINE_id,CINE),sbr_nm//": pv CINE"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,cld_frc_id,cld_frc),sbr_nm//": pv cld_frc"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,cnc_CFC11_id,cnc_CFC11),sbr_nm//": pv cnc_CFC11"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,cnc_CFC12_id,cnc_CFC12),sbr_nm//": pv cnc_CFC12"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,cnc_CH4_id,cnc_CH4),sbr_nm//": pv cnc_CH4"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,cnc_CO2_id,cnc_CO2),sbr_nm//": pv cnc_CO2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,cnc_H2OH2O_id,cnc_H2OH2O),sbr_nm//": pv cnc_H2OH2O"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,cnc_H2O_id,cnc_H2O),sbr_nm//": pv cnc_H2O"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,cnc_N2O_id,cnc_N2O),sbr_nm//": pv cnc_N2O"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,cnc_N2_id,cnc_N2),sbr_nm//": pv cnc_N2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,cnc_NO2_id,cnc_NO2),sbr_nm//": pv cnc_NO2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,cnc_O2O2_id,cnc_O2O2),sbr_nm//": pv cnc_O2O2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,cnc_O2_cnc_N2_id,cnc_O2_cnc_N2),sbr_nm//": pv cnc_O2_cnc_N2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,cnc_O2_cnc_O2_id,cnc_O2_cnc_O2),sbr_nm//": pv cnc_O2_cnc_O2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,cnc_O2_id,cnc_O2),sbr_nm//": pv cnc_O2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,cnc_O2_npl_N2_clm_id,cnc_O2_npl_N2_clm),sbr_nm//": pv cnc_O2_npl_N2_clm"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,cnc_O2_npl_N2_id,cnc_O2_npl_N2),sbr_nm//": pv cnc_O2_npl_N2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,cnc_O2_npl_O2_clm_frc_id,cnc_O2_npl_O2_clm_frc),sbr_nm//": pv cnc_O2_npl_O2_clm_frc"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,cnc_O2_npl_O2_clm_id,cnc_O2_npl_O2_clm),sbr_nm//": pv cnc_O2_npl_O2_clm"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,cnc_O2_npl_O2_id,cnc_O2_npl_O2),sbr_nm//": pv cnc_O2_npl_O2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,cnc_O3_id,cnc_O3),sbr_nm//": pv cnc_O3"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,cnc_OH_id,cnc_OH),sbr_nm//": pv cnc_OH"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,cnc_dry_air_id,cnc_dry_air),sbr_nm//": pv cnc_dry_air"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,cnc_mst_air_id,cnc_mst_air),sbr_nm//": pv cnc_mst_air"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,nrg_dry_id,nrg_dry),sbr_nm//": pv nrg_dry"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,nrg_mst_id,nrg_mst),sbr_nm//": pv nrg_mst"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,pcl_DLR_id,pcl_DLR),sbr_nm//": pv pcl_DLR"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,dns_CFC11_id,dns_CFC11),sbr_nm//": pv dns_CFC11"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,dns_CFC12_id,dns_CFC12),sbr_nm//": pv dns_CFC12"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,dns_CH4_id,dns_CH4),sbr_nm//": pv dns_CH4"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,dns_CO2_id,dns_CO2),sbr_nm//": pv dns_CO2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,dns_H2OH2O_id,dns_H2OH2O),sbr_nm//": pv dns_H2OH2O"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,dns_H2O_id,dns_H2O),sbr_nm//": pv dns_H2O"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,dns_N2O_id,dns_N2O),sbr_nm//": pv dns_N2O"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,dns_N2_id,dns_N2),sbr_nm//": pv dns_N2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,dns_NO2_id,dns_NO2),sbr_nm//": pv dns_NO2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,dns_O2O2_id,dns_O2O2),sbr_nm//": pv dns_O2O2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,dns_O2_dns_N2_id,dns_O2_dns_N2),sbr_nm//": pv dns_O2_dns_N2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,dns_O2_dns_O2_id,dns_O2_dns_O2),sbr_nm//": pv dns_O2_dns_O2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,dns_O2_id,dns_O2),sbr_nm//": pv dns_O2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,dns_O2_mpl_N2_clm_id,dns_O2_mpl_N2_clm),sbr_nm//": pv dns_O2_mpl_N2_clm"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,dns_O2_mpl_N2_id,dns_O2_mpl_N2),sbr_nm//": pv dns_O2_mpl_N2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,dns_O2_mpl_O2_clm_id,dns_O2_mpl_O2_clm),sbr_nm//": pv dns_O2_mpl_O2_clm"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,dns_O2_mpl_O2_id,dns_O2_mpl_O2),sbr_nm//": pv dns_O2_mpl_O2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,dns_O3_id,dns_O3),sbr_nm//": pv dns_O3"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,dns_OH_id,dns_OH),sbr_nm//": pv dns_OH"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,dns_aer_id,dns_aer),sbr_nm//": pv dns_aer"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,dns_bga_id,dns_bga),sbr_nm//": pv dns_bga"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,dns_dry_air_id,dns_dry_air),sbr_nm//": pv dns_dry_air"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,dns_mst_air_id,dns_mst_air),sbr_nm//": pv dns_mst_air"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,eqn_time_sec_id,eqn_time_sec),sbr_nm//": pv eqn_time_sec"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,ext_cff_mss_aer_id,ext_cff_mss_aer),sbr_nm//": pv ext_cff_mss_aer"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,ext_cff_mss_bga_id,ext_cff_mss_bga),sbr_nm//": pv ext_cff_mss_bga"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,frc_ice_id,frc_ice),sbr_nm//": pv frc_ice"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,frc_ice_ttl_id,frc_ice_ttl),sbr_nm//": pv frc_ice_ttl"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,frc_str_zen_ngl_sfc_id,frc_str_zen_ngl_sfc),sbr_nm//": pv frc_str_zen_ngl_sfc"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,gas_cst_mst_air_id,gas_cst_mst_air),sbr_nm//": pv gas_cst_mst_air"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,gmt_day_id,gmt_day),sbr_nm//": pv gmt_day"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,gmt_doy_id,gmt_doy),sbr_nm//": pv gmt_doy"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,gmt_hr_id,gmt_hr),sbr_nm//": pv gmt_hr"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,gmt_mnt_id,gmt_mnt),sbr_nm//": pv gmt_mnt"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,gmt_mth_id,gmt_mth),sbr_nm//": pv gmt_mth"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,gmt_sec_id,gmt_sec),sbr_nm//": pv gmt_sec"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,gmt_ydy_id,gmt_ydy),sbr_nm//": pv gmt_ydy"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,gmt_yr_id,gmt_yr),sbr_nm//": pv gmt_yr"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,grv_id,grv),sbr_nm//": pv grv"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,lat_cos_id,lat_cos),sbr_nm//": pv lat_cos"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,lat_dgr_id,lat_dgr),sbr_nm//": pv lat_dgr"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,lat_id,lat),sbr_nm//": pv lat"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,lev_LCL_id,lev_LCL),sbr_nm//": pv lev_LCL"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,lcl_time_hr_id,lcl_time_hr),sbr_nm//": pv lcl_time_hr"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,lcl_yr_day_id,lcl_yr_day),sbr_nm//": pv lcl_yr_day"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,lev_id,lev),sbr_nm//": pv lev"//__FILE__)
  if (flg_snw) then
     rcd=nf90_wrp(nf90_put_var(nc_id,dns_snw_id,dns_snw),sbr_nm//": pv dns_snw"//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,dpt_dlt_snw_id,dpt_dlt_snw),sbr_nm//": pv dpt_dlt_snw"//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,dpt_ntf_snw_id,dpt_ntf_snw),sbr_nm//": pv dpt_ntf_snw"//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,dpt_snw_id,dpt_snw),sbr_nm//": pv dpt_snw"//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,foo_snw_id,foo_snw),sbr_nm//": pv foo_snw"//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,lev_snw_id,lev_snw),sbr_nm//": pv lev_snw"//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,levp_snw_id,levp_snw),sbr_nm//": pv levp_snw"//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,mmr_mpr_snw_id,mmr_mpr_snw),sbr_nm//": pv mmr_mpr_snw"//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,rds_ffc_snw_id,rds_ffc_snw),sbr_nm//": pv rds_ffc_snw"//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,tpt_ntf_snw_id,tpt_ntf_snw),sbr_nm//": pv tpt_ntf_snw"//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,tpt_snw_id,tpt_snw),sbr_nm//": pv tpt_snw"//__FILE__)
  endif ! !flg_snw
  rcd=nf90_wrp(nf90_put_var(nc_id,levp_id,levp),sbr_nm//": pv levp"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,lev_LFC_id,lev_LFC),sbr_nm//": pv lev_LFC"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,lmt_day_id,lmt_day),sbr_nm//": pv lmt_day"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,lmt_hr_id,lmt_hr),sbr_nm//": pv lmt_hr"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,lmt_mnt_id,lmt_mnt),sbr_nm//": pv lmt_mnt"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,lmt_mth_id,lmt_mth),sbr_nm//": pv lmt_mth"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,lmt_sec_id,lmt_sec),sbr_nm//": pv lmt_sec"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,lmt_ydy_id,lmt_ydy),sbr_nm//": pv lmt_ydy"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,lmt_yr_id,lmt_yr),sbr_nm//": pv lmt_yr"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,lev_LNB_id,lev_LNB),sbr_nm//": pv lev_LNB"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,lon_dgr_id,lon_dgr),sbr_nm//": pv lon_dgr"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,lon_id,lon),sbr_nm//": pv lon"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,lon_sec_id,lon_sec),sbr_nm//": pv lon_sec"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,ltst_day_id,ltst_day),sbr_nm//": pv ltst_day"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,ltst_doy_id,ltst_doy),sbr_nm//": pv ltst_doy"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,ltst_hr_id,ltst_hr),sbr_nm//": pv ltst_hr"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,ltst_mnt_id,ltst_mnt),sbr_nm//": pv ltst_mnt"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,ltst_mth_id,ltst_mth),sbr_nm//": pv ltst_mth"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,ltst_sec_id,ltst_sec),sbr_nm//": pv ltst_sec"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,ltst_ydy_id,ltst_ydy),sbr_nm//": pv ltst_ydy"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,ltst_yr_id,ltst_yr),sbr_nm//": pv ltst_yr"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,mmw_mst_air_id,mmw_mst_air),sbr_nm//": pv mmw_mst_air"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,pcl_MLR_id,pcl_MLR),sbr_nm//": pv pcl_MLR"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,mpc_CFC11_id,mpc_CFC11),sbr_nm//": pv mpc_CFC11"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,mpc_CFC12_id,mpc_CFC12),sbr_nm//": pv mpc_CFC12"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,mpc_CH4_id,mpc_CH4),sbr_nm//": pv mpc_CH4"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,mpc_CO2_id,mpc_CO2),sbr_nm//": pv mpc_CO2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,mpc_CWP_id,mpc_CWP),sbr_nm//": pv mpc_CWP"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,mpc_H2OH2O_id,mpc_H2OH2O),sbr_nm//": pv mpc_H2OH2O"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,mpc_H2O_id,mpc_H2O),sbr_nm//": pv mpc_H2O"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,mpc_IWP_id,mpc_IWP),sbr_nm//": pv mpc_IWP"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,mpc_LWP_id,mpc_LWP),sbr_nm//": pv mpc_LWP"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,mpc_N2O_id,mpc_N2O),sbr_nm//": pv mpc_N2O"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,mpc_N2_id,mpc_N2),sbr_nm//": pv mpc_N2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,mpc_NO2_id,mpc_NO2),sbr_nm//": pv mpc_NO2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,mpc_O2O2_id,mpc_O2O2),sbr_nm//": pv mpc_O2O2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,mpc_O2_id,mpc_O2),sbr_nm//": pv mpc_O2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,mpc_O3_DU_id,mpc_O3_DU),sbr_nm//": pv mpc_O3_DU"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,mpc_O3_id,mpc_O3),sbr_nm//": pv mpc_O3"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,mpc_OH_id,mpc_OH),sbr_nm//": pv mpc_OH"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,mpc_aer_id,mpc_aer),sbr_nm//": pv mpc_aer"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,mpc_bga_id,mpc_bga),sbr_nm//": pv mpc_bga"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,mpc_dry_air_id,mpc_dry_air),sbr_nm//": pv mpc_dry_air"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,mpc_mst_air_id,mpc_mst_air),sbr_nm//": pv mpc_mst_air"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,mpl_CFC11_id,mpl_CFC11),sbr_nm//": pv mpl_CFC11"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,mpl_CFC12_id,mpl_CFC12),sbr_nm//": pv mpl_CFC12"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,mpl_CH4_id,mpl_CH4),sbr_nm//": pv mpl_CH4"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,mpl_CO2_id,mpl_CO2),sbr_nm//": pv mpl_CO2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,mpl_CWP_id,mpl_CWP),sbr_nm//": pv mpl_CWP"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,mpl_H2OH2O_id,mpl_H2OH2O),sbr_nm//": pv mpl_H2OH2O"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,mpl_H2O_id,mpl_H2O),sbr_nm//": pv mpl_H2O"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,mpl_IWP_id,mpl_IWP),sbr_nm//": pv mpl_IWP"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,mpl_LWP_id,mpl_LWP),sbr_nm//": pv mpl_LWP"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,mpl_N2O_id,mpl_N2O),sbr_nm//": pv mpl_N2O"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,mpl_N2_id,mpl_N2),sbr_nm//": pv mpl_N2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,mpl_NO2_id,mpl_NO2),sbr_nm//": pv mpl_NO2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,mpl_O2O2_id,mpl_O2O2),sbr_nm//": pv mpl_O2O2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,mpl_O2_id,mpl_O2),sbr_nm//": pv mpl_O2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,mpl_O3_id,mpl_O3),sbr_nm//": pv mpl_O3"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,mpl_OH_id,mpl_OH),sbr_nm//": pv mpl_OH"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,mpl_aer_id,mpl_aer),sbr_nm//": pv mpl_aer"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,mpl_bga_id,mpl_bga),sbr_nm//": pv mpl_bga"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,mpl_dry_air_id,mpl_dry_air),sbr_nm//": pv mpl_dry_air"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,mpl_mst_air_id,mpl_mst_air),sbr_nm//": pv mpl_mst_air"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,npc_CFC11_id,npc_CFC11),sbr_nm//": pv npc_CFC11"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,npc_CFC12_id,npc_CFC12),sbr_nm//": pv npc_CFC12"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,npc_CH4_id,npc_CH4),sbr_nm//": pv npc_CH4"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,npc_CO2_id,npc_CO2),sbr_nm//": pv npc_CO2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,npc_H2OH2O_id,npc_H2OH2O),sbr_nm//": pv npc_H2OH2O"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,npc_H2O_id,npc_H2O),sbr_nm//": pv npc_H2O"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,npc_N2O_id,npc_N2O),sbr_nm//": pv npc_N2O"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,npc_N2_id,npc_N2),sbr_nm//": pv npc_N2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,npc_NO2_id,npc_NO2),sbr_nm//": pv npc_NO2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,npc_O2O2_id,npc_O2O2),sbr_nm//": pv npc_O2O2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,npc_O2_id,npc_O2),sbr_nm//": pv npc_O2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,npc_O3_id,npc_O3),sbr_nm//": pv npc_O3"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,npc_OH_id,npc_OH),sbr_nm//": pv npc_OH"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,npc_dry_air_id,npc_dry_air),sbr_nm//": pv npc_dry_air"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,npc_mst_air_id,npc_mst_air),sbr_nm//": pv npc_mst_air"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,npl_CFC11_id,npl_CFC11),sbr_nm//": pv npl_CFC11"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,npl_CFC12_id,npl_CFC12),sbr_nm//": pv npl_CFC12"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,npl_CH4_id,npl_CH4),sbr_nm//": pv npl_CH4"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,npl_CO2_id,npl_CO2),sbr_nm//": pv npl_CO2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,npl_H2OH2O_id,npl_H2OH2O),sbr_nm//": pv npl_H2OH2O"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,npl_H2O_id,npl_H2O),sbr_nm//": pv npl_H2O"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,npl_N2O_id,npl_N2O),sbr_nm//": pv npl_N2O"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,npl_N2_id,npl_N2),sbr_nm//": pv npl_N2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,npl_NO2_id,npl_NO2),sbr_nm//": pv npl_NO2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,npl_O2O2_id,npl_O2O2),sbr_nm//": pv npl_O2O2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,npl_O2_id,npl_O2),sbr_nm//": pv npl_O2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,npl_O3_id,npl_O3),sbr_nm//": pv npl_O3"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,npl_OH_id,npl_OH),sbr_nm//": pv npl_OH"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,npl_dry_air_id,npl_dry_air),sbr_nm//": pv npl_dry_air"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,npl_mst_air_id,npl_mst_air),sbr_nm//": pv npl_mst_air"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,odxc_obs_aer_id,odxc_obs_aer),sbr_nm//": pv odxc_obs_aer"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,odxc_obs_bga_id,odxc_obs_bga),sbr_nm//": pv odxc_obs_bga"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,odxl_obs_aer_id,odxl_obs_aer),sbr_nm//": pv odxl_obs_aer"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,odxl_obs_bga_id,odxl_obs_bga),sbr_nm//": pv odxl_obs_bga"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,oneD_foo_id,oneD_foo),sbr_nm//": pv oneD_foo"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,oro_id,oro),sbr_nm//": pv oro"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,pcl_PLR_id,pcl_PLR),sbr_nm//": pv pcl_PLR"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,ppr_CFC11_id,ppr_CFC11),sbr_nm//": pv ppr_CFC11"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,ppr_CFC12_id,ppr_CFC12),sbr_nm//": pv ppr_CFC12"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,ppr_CH4_id,ppr_CH4),sbr_nm//": pv ppr_CH4"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,ppr_CO2_id,ppr_CO2),sbr_nm//": pv ppr_CO2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,ppr_H2OH2O_id,ppr_H2OH2O),sbr_nm//": pv ppr_H2OH2O"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,ppr_H2O_id,ppr_H2O),sbr_nm//": pv ppr_H2O"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,ppr_N2O_id,ppr_N2O),sbr_nm//": pv ppr_N2O"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,ppr_N2_id,ppr_N2),sbr_nm//": pv ppr_N2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,ppr_NO2_id,ppr_NO2),sbr_nm//": pv ppr_NO2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,ppr_O2O2_id,ppr_O2O2),sbr_nm//": pv ppr_O2O2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,ppr_O2_id,ppr_O2),sbr_nm//": pv ppr_O2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,ppr_O3_id,ppr_O3),sbr_nm//": pv ppr_O3"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,ppr_OH_id,ppr_OH),sbr_nm//": pv ppr_OH"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,ppr_dry_air_id,ppr_dry_air),sbr_nm//": pv ppr_dry_air"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,prs_cld_btm_id,prs_cld_btm),sbr_nm//": pv prs_cld_btm"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,prs_cld_mid_id,prs_cld_mid),sbr_nm//": pv prs_cld_mid"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,prs_cld_thick_id,prs_cld_thick),sbr_nm//": pv prs_cld_thick"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,prs_cld_top_id,prs_cld_top),sbr_nm//": pv prs_cld_top"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,prs_dlt_id,prs_dlt),sbr_nm//": pv prs_dlt"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,prs_id,prs),sbr_nm//": pv prs"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,prs_ntf_id,prs_ntf),sbr_nm//": pv prs_ntf"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,prs_sfc_id,prs_sfc),sbr_nm//": pv prs_sfc"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,q_CFC11_id,q_CFC11),sbr_nm//": pv q_CFC11"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,q_CFC12_id,q_CFC12),sbr_nm//": pv q_CFC12"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,q_CH4_id,q_CH4),sbr_nm//": pv q_CH4"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,q_CO2_id,q_CO2),sbr_nm//": pv q_CO2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,q_H2OH2O_id,q_H2OH2O),sbr_nm//": pv q_H2OH2O"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,q_H2OH2O_rcp_q_H2O_id,q_H2OH2O_rcp_q_H2O),sbr_nm//": pv q_H2OH2O_rcp_q_H2O"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,q_H2O_id,q_H2O),sbr_nm//": pv q_H2O"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,q_N2O_id,q_N2O),sbr_nm//": pv q_N2O"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,q_N2_id,q_N2),sbr_nm//": pv q_N2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,q_NO2_id,q_NO2),sbr_nm//": pv q_NO2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,q_O2O2_id,q_O2O2),sbr_nm//": pv q_O2O2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,q_O2_id,q_O2),sbr_nm//": pv q_O2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,q_O3_id,q_O3),sbr_nm//": pv q_O3"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,q_OH_id,q_OH),sbr_nm//": pv q_OH"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,qst_H2O_ice_id,qst_H2O_ice),sbr_nm//": pv qst_H2O_ice"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,qst_H2O_lqd_id,qst_H2O_lqd),sbr_nm//": pv qst_H2O_lqd"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,r_CFC11_id,r_CFC11),sbr_nm//": pv r_CFC11"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,r_CFC12_id,r_CFC12),sbr_nm//": pv r_CFC12"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,r_CH4_id,r_CH4),sbr_nm//": pv r_CH4"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,r_CO2_id,r_CO2),sbr_nm//": pv r_CO2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,r_H2OH2O_id,r_H2OH2O),sbr_nm//": pv r_H2OH2O"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,r_H2O_id,r_H2O),sbr_nm//": pv r_H2O"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,r_N2O_id,r_N2O),sbr_nm//": pv r_N2O"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,r_N2_id,r_N2),sbr_nm//": pv r_N2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,r_NO2_id,r_NO2),sbr_nm//": pv r_NO2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,r_O2O2_id,r_O2O2),sbr_nm//": pv r_O2O2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,r_O2_id,r_O2),sbr_nm//": pv r_O2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,r_O3_id,r_O3),sbr_nm//": pv r_O3"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,r_OH_id,r_OH),sbr_nm//": pv r_OH"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,rds_fct_ice_id,rds_fct_ice),sbr_nm//": pv rds_fct_ice"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,rds_fct_lqd_id,rds_fct_lqd),sbr_nm//": pv rds_fct_lqd"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,rgh_len_id,rgh_len),sbr_nm//": pv rgh_len"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,scl_hgt_id,scl_hgt),sbr_nm//": pv scl_hgt"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,sfc_ems_id,sfc_ems),sbr_nm//": pv sfc_ems"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,slr_azi_dgr_id,slr_azi_dgr),sbr_nm//": pv slr_azi_dgr"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,slr_crd_gmm_dgr_id,slr_crd_gmm_dgr),sbr_nm//": pv slr_crd_gmm_dgr"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,slr_cst_id,slr_cst),sbr_nm//": pv slr_cst"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,slr_dcl_dgr_id,slr_dcl_dgr),sbr_nm//": pv slr_dcl_dgr"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,slr_dmt_dgr_id,slr_dmt_dgr),sbr_nm//": pv slr_dmt_dgr"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,slr_dst_au_id,slr_dst_au),sbr_nm//": pv slr_dst_au"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,slr_elv_dgr_id,slr_elv_dgr),sbr_nm//": pv slr_elv_dgr"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,slr_flx_TOA_id,slr_flx_TOA),sbr_nm//": pv slr_flx_TOA"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,slr_flx_nrm_TOA_id,slr_flx_nrm_TOA),sbr_nm//": pv slr_flx_nrm_TOA"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,slr_hr_ngl_dgr_id,slr_hr_ngl_dgr),sbr_nm//": pv slr_hr_ngl_dgr"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,slr_rfr_ngl_dgr_id,slr_rfr_ngl_dgr),sbr_nm//": pv slr_rfr_ngl_dgr"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,slr_rgt_asc_dgr_id,slr_rgt_asc_dgr),sbr_nm//": pv slr_rgt_asc_dgr"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,slr_zen_ngl_cos_id,slr_zen_ngl_cos),sbr_nm//": pv slr_zen_ngl_cos"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,slr_zen_ngl_dgr_id,slr_zen_ngl_dgr),sbr_nm//": pv slr_zen_ngl_dgr"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,slr_zen_ngl_id,slr_zen_ngl),sbr_nm//": pv slr_zen_ngl"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,snow_depth_id,snow_depth),sbr_nm//": pv snow_depth"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,spc_heat_mst_air_id,spc_heat_mst_air),sbr_nm//": pv spc_heat_mst_air"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,time_lmt_id,time_lmt),sbr_nm//": pv time_lmt"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,time_ltst_id,time_ltst),sbr_nm//": pv time_ltst"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,time_unix_id,time_unix),sbr_nm//": pv time_unix"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,tpt_cls_id,tpt_cls),sbr_nm//": pv tpt_cls"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,tpt_cls_ntf_id,tpt_cls_ntf),sbr_nm//": pv tpt_cls_ntf"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,tpt_dwp_id,tpt_dwp),sbr_nm//": pv tpt_dwp"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,tpt_id,tpt),sbr_nm//": pv tpt"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,tpt_ntf_id,tpt_ntf),sbr_nm//": pv tpt_ntf"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,tpt_pcl_id,tpt_pcl),sbr_nm//": pv tpt_pcl"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,tpt_pcl_sfc_id,tpt_pcl_sfc),sbr_nm//": pv tpt_pcl_sfc"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,tpt_ptn_id,tpt_ptn),sbr_nm//": pv tpt_ptn"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,tpt_ptn_vrt_id,tpt_ptn_vrt),sbr_nm//": pv tpt_ptn_vrt"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,tpt_sfc_id,tpt_sfc),sbr_nm//": pv tpt_sfc"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,tpt_skn_id,tpt_skn),sbr_nm//": pv tpt_skn"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,tpt_vrt_id,tpt_vrt),sbr_nm//": pv tpt_vrt"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,vmr_CFC11_id,vmr_CFC11),sbr_nm//": pv vmr_CFC11"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,vmr_CFC12_id,vmr_CFC12),sbr_nm//": pv vmr_CFC12"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,vmr_CH4_id,vmr_CH4),sbr_nm//": pv vmr_CH4"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,vmr_CO2_id,vmr_CO2),sbr_nm//": pv vmr_CO2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,vmr_H2OH2O_id,vmr_H2OH2O),sbr_nm//": pv vmr_H2OH2O"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,vmr_H2O_id,vmr_H2O),sbr_nm//": pv vmr_H2O"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,vmr_N2O_id,vmr_N2O),sbr_nm//": pv vmr_N2O"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,vmr_N2_id,vmr_N2),sbr_nm//": pv vmr_N2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,vmr_NO2_id,vmr_NO2),sbr_nm//": pv vmr_NO2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,vmr_O2O2_id,vmr_O2O2),sbr_nm//": pv vmr_O2O2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,vmr_O2_id,vmr_O2),sbr_nm//": pv vmr_O2"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,vmr_O3_id,vmr_O3),sbr_nm//": pv vmr_O3"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,vmr_OH_id,vmr_OH),sbr_nm//": pv vmr_OH"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,wvl_obs_aer_id,wvl_obs_aer),sbr_nm//": pv wvl_obs_aer"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,wvl_obs_bga_id,wvl_obs_bga),sbr_nm//": pv wvl_obs_bga"//__FILE__)
  rcd=nf90_wrp(nf90_put_var(nc_id,xnt_fac_id,xnt_fac),sbr_nm//": pv xnt_fac"//__FILE__)
  if (flg_wnd) then
     rcd=nf90_wrp(nf90_put_var(nc_id,abl_hgt_id,abl_hgt),sbr_nm//": pv abl_hgt"//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,kbfs_id,kbfs),sbr_nm//": pv kbfs"//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,khfs_id,khfs),sbr_nm//": pv khfs"//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,kqfs_id,kqfs),sbr_nm//": pv kqfs"//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,lhflx_id,lhflx),sbr_nm//": pv lhflx"//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,obk_len_id,obk_len),sbr_nm//": pv obk_len"//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,ric_nbr_id,ric_nbr),sbr_nm//": pv ric_nbr"//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,shflx_id,shflx),sbr_nm//": pv shflx"//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,tau_id,tau),sbr_nm//": pv tau"//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,tke_id,tke),sbr_nm//": pv tke"//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,ustar_id,ustar),sbr_nm//": pv ustar"//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,wnd_dir_id,wnd_dir),sbr_nm//": pv wnd_dir"//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,wnd_shr_id,wnd_shr),sbr_nm//": pv wnd_shr"//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,wnd_spd_id,wnd_spd),sbr_nm//": pv wnd_spd"//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,wspd10m_id,wspd10m),sbr_nm//": pv wspd10m"//__FILE__)
     rcd=nf90_wrp(nf90_put_var(nc_id,wstar_id,wstar),sbr_nm//": pv wstar"//__FILE__)
  end if
  ! Close file
  rcd=nf90_wrp_close(nc_id,fl_out,'Wrote results to') ! [fnc] Close file
  
  ! De-allocate dynamic variables
  if (allocated(PDF_bin_spd)) deallocate(PDF_bin_spd,stat=rcd)
  if (allocated(PDF_bin_wgt)) deallocate(PDF_bin_wgt,stat=rcd)
  if (allocated(RH)) deallocate(RH,stat=rcd)
  if (allocated(RH_ice)) deallocate(RH_ice,stat=rcd)
  if (allocated(RH_lqd)) deallocate(RH_lqd,stat=rcd)
  if (allocated(alt)) deallocate(alt,stat=rcd)
  if (allocated(alt_dlt)) deallocate(alt_dlt,stat=rcd)
  if (allocated(alt_ntf)) deallocate(alt_ntf,stat=rcd)
  if (allocated(bga_lvl_flg)) deallocate(bga_lvl_flg,stat=rcd)
  if (allocated(cld_frc)) deallocate(cld_frc,stat=rcd)
  if (allocated(cnc_CFC11)) deallocate(cnc_CFC11,stat=rcd)
  if (allocated(cnc_CFC12)) deallocate(cnc_CFC12,stat=rcd)
  if (allocated(cnc_CH4)) deallocate(cnc_CH4,stat=rcd)
  if (allocated(cnc_CO2)) deallocate(cnc_CO2,stat=rcd)
  if (allocated(cnc_H2O)) deallocate(cnc_H2O,stat=rcd)
  if (allocated(cnc_H2OH2O)) deallocate(cnc_H2OH2O,stat=rcd)
  if (allocated(cnc_N2)) deallocate(cnc_N2,stat=rcd)
  if (allocated(cnc_N2O)) deallocate(cnc_N2O,stat=rcd)
  if (allocated(cnc_NO2)) deallocate(cnc_NO2,stat=rcd)
  if (allocated(cnc_O2)) deallocate(cnc_O2,stat=rcd)
  if (allocated(cnc_O2O2)) deallocate(cnc_O2O2,stat=rcd)
  if (allocated(cnc_O2_cnc_N2)) deallocate(cnc_O2_cnc_N2,stat=rcd) ! [mlc2 m-6] O2 number concentration times N2 number concentration
  if (allocated(cnc_O2_cnc_O2)) deallocate(cnc_O2_cnc_O2,stat=rcd) ! [mlc2 m-6] O2 number concentration squared
  if (allocated(cnc_O2_npl_N2)) deallocate(cnc_O2_npl_N2,stat=rcd) ! [mlc2 m-5] O2 number concentration times N2 number path
  if (allocated(cnc_O2_npl_O2)) deallocate(cnc_O2_npl_O2,stat=rcd) ! [mlc2 m-5] O2 number concentration times O2 number path
  if (allocated(cnc_O2_npl_O2_clm_frc)) deallocate(cnc_O2_npl_O2_clm_frc,stat=rcd) ! [frc] Fraction of column O2-O2 at or above each level
  if (allocated(cnc_O3)) deallocate(cnc_O3,stat=rcd)
  if (allocated(cnc_OH)) deallocate(cnc_OH,stat=rcd)
  if (allocated(cnc_dry_air)) deallocate(cnc_dry_air,stat=rcd)
  if (allocated(cnc_mst_air)) deallocate(cnc_mst_air,stat=rcd)
  if (allocated(dns_CFC11)) deallocate(dns_CFC11,stat=rcd)
  if (allocated(dns_CFC12)) deallocate(dns_CFC12,stat=rcd)
  if (allocated(dns_CH4)) deallocate(dns_CH4,stat=rcd)
  if (allocated(dns_CO2)) deallocate(dns_CO2,stat=rcd)
  if (allocated(dns_H2O)) deallocate(dns_H2O,stat=rcd)
  if (allocated(dns_H2OH2O)) deallocate(dns_H2OH2O,stat=rcd)
  if (allocated(dns_N2)) deallocate(dns_N2,stat=rcd)
  if (allocated(dns_N2O)) deallocate(dns_N2O,stat=rcd)
  if (allocated(dns_NO2)) deallocate(dns_NO2,stat=rcd)
  if (allocated(dns_O2)) deallocate(dns_O2,stat=rcd)
  if (allocated(dns_O2O2)) deallocate(dns_O2O2,stat=rcd)
  if (allocated(dns_O2_dns_N2)) deallocate(dns_O2_dns_N2,stat=rcd) ! [kg2 m-6] O2 mass concentration times N2 mass concentration
  if (allocated(dns_O2_dns_O2)) deallocate(dns_O2_dns_O2,stat=rcd) ! [kg2 m-6] O2 mass concentration squared
  if (allocated(dns_O2_mpl_N2)) deallocate(dns_O2_mpl_N2,stat=rcd) ! [kg2 m-5] O2 mass concentration times N2 mass path
  if (allocated(dns_O2_mpl_O2 )) deallocate(dns_O2_mpl_O2 ,stat=rcd) ! [kg2 m-5] O2 mass concentration times O2 mass path
  if (allocated(dns_O3)) deallocate(dns_O3,stat=rcd)
  if (allocated(dns_OH)) deallocate(dns_OH,stat=rcd)
  if (allocated(dns_dry_air)) deallocate(dns_dry_air,stat=rcd)
  if (allocated(dns_mst_air)) deallocate(dns_mst_air,stat=rcd)
  if (allocated(frc_ice)) deallocate(frc_ice,stat=rcd)
  if (allocated(gas_cst_mst_air)) deallocate(gas_cst_mst_air,stat=rcd)
  if (allocated(grv)) deallocate(grv,stat=rcd)
  if (allocated(lev)) deallocate(lev,stat=rcd)     ! coordinate variable
  if (allocated(mmw_mst_air)) deallocate(mmw_mst_air,stat=rcd)
  if (allocated(mpl_CFC11)) deallocate(mpl_CFC11,stat=rcd)
  if (allocated(mpl_CFC12)) deallocate(mpl_CFC12,stat=rcd)
  if (allocated(mpl_CH4)) deallocate(mpl_CH4,stat=rcd)
  if (allocated(mpl_CO2)) deallocate(mpl_CO2,stat=rcd)
  if (allocated(mpl_CWP)) deallocate(mpl_CWP,stat=rcd)
  if (allocated(mpl_H2O)) deallocate(mpl_H2O,stat=rcd)
  if (allocated(mpl_H2OH2O)) deallocate(mpl_H2OH2O,stat=rcd)
  if (allocated(mpl_IWP)) deallocate(mpl_IWP,stat=rcd)
  if (allocated(mpl_LWP)) deallocate(mpl_LWP,stat=rcd)
  if (allocated(mpl_N2)) deallocate(mpl_N2,stat=rcd)
  if (allocated(mpl_N2O)) deallocate(mpl_N2O,stat=rcd)
  if (allocated(mpl_NO2)) deallocate(mpl_NO2,stat=rcd)
  if (allocated(mpl_O2)) deallocate(mpl_O2,stat=rcd)
  if (allocated(mpl_O2O2)) deallocate(mpl_O2O2,stat=rcd)
  if (allocated(mpl_O3)) deallocate(mpl_O3,stat=rcd)
  if (allocated(mpl_OH)) deallocate(mpl_OH,stat=rcd)
  if (allocated(mpl_aer)) deallocate(mpl_aer,stat=rcd)
  if (allocated(mpl_bga)) deallocate(mpl_bga,stat=rcd)
  if (allocated(mpl_dry_air)) deallocate(mpl_dry_air,stat=rcd)
  if (allocated(mpl_mst_air)) deallocate(mpl_mst_air,stat=rcd)
  if (allocated(npl_CFC11)) deallocate(npl_CFC11,stat=rcd)
  if (allocated(npl_CFC12)) deallocate(npl_CFC12,stat=rcd)
  if (allocated(npl_CH4)) deallocate(npl_CH4,stat=rcd)
  if (allocated(npl_CO2)) deallocate(npl_CO2,stat=rcd)
  if (allocated(npl_H2O)) deallocate(npl_H2O,stat=rcd)
  if (allocated(npl_H2OH2O)) deallocate(npl_H2OH2O,stat=rcd)
  if (allocated(npl_N2)) deallocate(npl_N2,stat=rcd)
  if (allocated(npl_N2O)) deallocate(npl_N2O,stat=rcd)
  if (allocated(npl_NO2)) deallocate(npl_NO2,stat=rcd)
  if (allocated(npl_O2)) deallocate(npl_O2,stat=rcd)
  if (allocated(npl_O2O2)) deallocate(npl_O2O2,stat=rcd)
  if (allocated(npl_O3)) deallocate(npl_O3,stat=rcd)
  if (allocated(npl_OH)) deallocate(npl_OH,stat=rcd)
  if (allocated(npl_dry_air)) deallocate(npl_dry_air,stat=rcd)
  if (allocated(npl_mst_air)) deallocate(npl_mst_air,stat=rcd)
  if (allocated(nrg_dry)) deallocate(nrg_dry,stat=rcd)
  if (allocated(nrg_mst)) deallocate(nrg_mst,stat=rcd)
  if (allocated(odxl_obs_aer)) deallocate(odxl_obs_aer,stat=rcd)
  if (allocated(odxl_obs_bga)) deallocate(odxl_obs_bga,stat=rcd)
  if (allocated(oneD_foo)) deallocate(oneD_foo,stat=rcd)
  if (allocated(pcl_DLR)) deallocate(pcl_DLR,stat=rcd)
  if (allocated(pcl_MLR)) deallocate(pcl_MLR,stat=rcd)
  if (allocated(pcl_PLR)) deallocate(pcl_PLR,stat=rcd)
  if (allocated(ppr_CFC11)) deallocate(ppr_CFC11,stat=rcd)
  if (allocated(ppr_CFC12)) deallocate(ppr_CFC12,stat=rcd)
  if (allocated(ppr_CH4)) deallocate(ppr_CH4,stat=rcd)
  if (allocated(ppr_CO2)) deallocate(ppr_CO2,stat=rcd)
  if (allocated(ppr_H2O)) deallocate(ppr_H2O,stat=rcd)
  if (allocated(ppr_H2OH2O)) deallocate(ppr_H2OH2O,stat=rcd)
  if (allocated(ppr_N2)) deallocate(ppr_N2,stat=rcd)
  if (allocated(ppr_N2O)) deallocate(ppr_N2O,stat=rcd)
  if (allocated(ppr_NO2)) deallocate(ppr_NO2,stat=rcd)
  if (allocated(ppr_O2)) deallocate(ppr_O2,stat=rcd)
  if (allocated(ppr_O2O2)) deallocate(ppr_O2O2,stat=rcd)
  if (allocated(ppr_O3)) deallocate(ppr_O3,stat=rcd)
  if (allocated(ppr_OH)) deallocate(ppr_OH,stat=rcd)
  if (allocated(ppr_dry_air)) deallocate(ppr_dry_air,stat=rcd)
  if (allocated(prs)) deallocate(prs,stat=rcd)
  if (allocated(prs_dlt)) deallocate(prs_dlt,stat=rcd)
  if (allocated(prs_ntf)) deallocate(prs_ntf,stat=rcd)
  if (allocated(prs_ttl)) deallocate(prs_ttl,stat=rcd)
  if (allocated(q_CFC11)) deallocate(q_CFC11,stat=rcd)
  if (allocated(q_CFC12)) deallocate(q_CFC12,stat=rcd)
  if (allocated(q_CH4)) deallocate(q_CH4,stat=rcd)
  if (allocated(q_CO2)) deallocate(q_CO2,stat=rcd)
  if (allocated(q_H2O)) deallocate(q_H2O,stat=rcd)
  if (allocated(q_H2OH2O)) deallocate(q_H2OH2O,stat=rcd)
  if (allocated(q_H2OH2O_rcp_q_H2O)) deallocate(q_H2OH2O_rcp_q_H2O,stat=rcd)
  if (allocated(q_N2)) deallocate(q_N2,stat=rcd)
  if (allocated(q_N2O)) deallocate(q_N2O,stat=rcd)
  if (allocated(q_NO2)) deallocate(q_NO2,stat=rcd)
  if (allocated(q_O2)) deallocate(q_O2,stat=rcd)
  if (allocated(q_O2O2)) deallocate(q_O2O2,stat=rcd)
  if (allocated(q_O3)) deallocate(q_O3,stat=rcd)
  if (allocated(q_OH)) deallocate(q_OH,stat=rcd)
  if (allocated(qst_H2O_ice)) deallocate(qst_H2O_ice,stat=rcd)
  if (allocated(qst_H2O_lqd)) deallocate(qst_H2O_lqd,stat=rcd)
  if (allocated(r_CFC11)) deallocate(r_CFC11,stat=rcd)
  if (allocated(r_CFC12)) deallocate(r_CFC12,stat=rcd)
  if (allocated(r_CH4)) deallocate(r_CH4,stat=rcd)
  if (allocated(r_CO2)) deallocate(r_CO2,stat=rcd)
  if (allocated(r_H2O)) deallocate(r_H2O,stat=rcd)
  if (allocated(r_H2OH2O)) deallocate(r_H2OH2O,stat=rcd)
  if (allocated(r_N2)) deallocate(r_N2,stat=rcd)
  if (allocated(r_N2O)) deallocate(r_N2O,stat=rcd)
  if (allocated(r_NO2)) deallocate(r_NO2,stat=rcd)
  if (allocated(r_O2)) deallocate(r_O2,stat=rcd)
  if (allocated(r_O2O2)) deallocate(r_O2O2,stat=rcd)
  if (allocated(r_O3)) deallocate(r_O3,stat=rcd)
  if (allocated(r_OH)) deallocate(r_OH,stat=rcd)
  if (allocated(rds_fct_ice)) deallocate(rds_fct_ice,stat=rcd)
  if (allocated(rds_fct_lqd)) deallocate(rds_fct_lqd,stat=rcd)
  if (allocated(scl_hgt)) deallocate(scl_hgt,stat=rcd)
  if (allocated(spc_heat_mst_air)) deallocate(spc_heat_mst_air,stat=rcd)
  if (allocated(tpt)) deallocate(tpt,stat=rcd)
  if (allocated(tpt_cls)) deallocate(tpt_cls,stat=rcd)
  if (allocated(tpt_cls_ntf)) deallocate(tpt_cls_ntf,stat=rcd)
  if (allocated(tpt_dwp)) deallocate(tpt_dwp,stat=rcd)
  if (allocated(tpt_ntf)) deallocate(tpt_ntf,stat=rcd)
  if (allocated(tpt_pcl)) deallocate(tpt_pcl,stat=rcd)
  if (allocated(tpt_ptn)) deallocate(tpt_ptn,stat=rcd)
  if (allocated(tpt_ptn_vrt)) deallocate(tpt_ptn_vrt,stat=rcd)
  if (allocated(tpt_vrt)) deallocate(tpt_vrt,stat=rcd)
  if (allocated(vmr_CFC11)) deallocate(vmr_CFC11,stat=rcd)
  if (allocated(vmr_CFC12)) deallocate(vmr_CFC12,stat=rcd)
  if (allocated(vmr_CH4)) deallocate(vmr_CH4,stat=rcd)
  if (allocated(vmr_CO2)) deallocate(vmr_CO2,stat=rcd)
  if (allocated(vmr_H2O)) deallocate(vmr_H2O,stat=rcd)
  if (allocated(vmr_H2OH2O)) deallocate(vmr_H2OH2O,stat=rcd)
  if (allocated(vmr_N2)) deallocate(vmr_N2,stat=rcd)
  if (allocated(vmr_N2O)) deallocate(vmr_N2O,stat=rcd)
  if (allocated(vmr_NO2)) deallocate(vmr_NO2,stat=rcd)
  if (allocated(vmr_O2)) deallocate(vmr_O2,stat=rcd)
  if (allocated(vmr_O2O2)) deallocate(vmr_O2O2,stat=rcd)
  if (allocated(vmr_O3)) deallocate(vmr_O3,stat=rcd)
  if (allocated(vmr_OH)) deallocate(vmr_OH,stat=rcd)
  if (flg_snw) then
     if (allocated(dns_snw)) deallocate(dns_snw,stat=rcd)
     if (allocated(dpt_dlt_snw)) deallocate(dpt_dlt_snw,stat=rcd)
     if (allocated(dpt_ntf_snw)) deallocate(dpt_ntf_snw,stat=rcd)
     if (allocated(dpt_snw)) deallocate(dpt_snw,stat=rcd)
     if (allocated(foo_snw)) deallocate(foo_snw,stat=rcd)
     if (allocated(lev_snw)) deallocate(lev_snw,stat=rcd)     ! coordinate variable
     if (allocated(levp)) deallocate(levp,stat=rcd)   ! coordinate variable
     if (allocated(levp_snw)) deallocate(levp_snw,stat=rcd)     ! coordinate variable
     if (allocated(mmr_mpr_snw)) deallocate(mmr_mpr_snw,stat=rcd)
     if (allocated(rds_ffc_snw)) deallocate(rds_ffc_snw,stat=rcd)
     if (allocated(tpt_ntf_snw)) deallocate(tpt_ntf_snw,stat=rcd)
     if (allocated(tpt_snw)) deallocate(tpt_snw,stat=rcd)
  endif ! !flg_snw
  if (flg_wnd) then
     if (allocated(ric_nbr)) deallocate(ric_nbr,stat=rcd)
     if (allocated(tke)) deallocate(tke,stat=rcd)
     if (allocated(wnd_dir)) deallocate(wnd_dir,stat=rcd)
     if (allocated(wnd_shr)) deallocate(wnd_shr,stat=rcd)
     if (allocated(wnd_spd)) deallocate(wnd_spd,stat=rcd)
  end if ! !flg_wnd
1000 continue
  
  call exit(exit_status)
end program clm                       ! end clm()

