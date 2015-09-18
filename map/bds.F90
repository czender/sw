! $Id$ -*-f90-*-

! Purpose: Generate Boundary DataSets datasets for use by global models
!          Identify and parameterize dust sources
  
! Copyright (C) 1998--2012 Charlie Zender
! This software is distributed under the terms of the GNU General Public License
! See http://www.gnu.org/copyleft/gpl.html for full license text
! The original author of this software, Charlie Zender, seeks to improve
! it with your suggestions, contributions, bug-reports, and patches.
! Charlie Zender <zender at uci dot edu>
! Department of Earth System Science
! University of California, Irvine
! Irvine, CA 92697-3100

! This program creates boundary dataset (bds) files used by the 
! Dust Entrainment And Deposition (DEAD) model (http://dust.ess.uci.edu/dead)
! This boundary dataset file contains surface type information, TOMS
! climatology, and model-specific tuning parameters. 
! Typically, this program produces files named, e.g., dst_T42.nc which 
! is placed in a directory accessible to the Dust model at run-time.
! For high resolution datasets (e.g., > T106) this program may fail with
! a segmentation violation because it is a memory hog.
! Easiest solution is to run on machines with lots of RAM
! 19990128: T170 currently requires 123 MB RAM

! Known problems:
! Fails in hgt_sfc_get() when available memory is insufficient
! Use "ulimit -s unlimited" to fix (not allowed in many interactive shells)
! Solaris machines inexplicably exit on error reading mineralogy files

#if 0 /* Compilation */
! bds is intimately linked to the fortran library libcsz_f90
! These two pieces must be compiled with the same intrinsic precisions
! Debugging compilation:
export NETCDF_ROOT='/usr';export NETCDF4_ROOT='/usr'
cd ${HOME}/f  ;make cln;make OPTS=X OMP=N PRC=S
cd ${HOME}/map;make cln;make OPTS=X OMP=N PRC=S
! Production compilation:
cd ${HOME}/f  ;make cln;make OPTS=O OMP=Y PRC=D
cd ${HOME}/map;make cln;make OPTS=O OMP=Y PRC=D
#endif /* !0 */ 

#if 0 /* Debugging usage */
bds -D 4 -t 19801980 -o ${DATA}/map/map.nc # Read one year of TOMS data (1980)
bds -D 4 -t 19801989 -o ${DATA}/map/map.nc # Read ten years of TOMS data (1980--1989)
bds -D 4 -x 3 -y 2 -m 13 -o ${DATA}/map/map.nc
bds -D 4 -x 90 -y 90 -m 13 -o ${DATA}/map/map.nc
! TOMS 1 x 1.5 degree grid:
bds -D 4 -x 288 -y 180 -m 13 -o ${DATA}/map/map.nc
! CCM T42 grid:
bds -D 4 -x 128 -y 64 -m 23 -o ${DATA}/map/map.nc
! ESMF
cd ~/map
bds -D 1 --drc_in=${DATA}/map -x 128 -y 64 -m 23 -o ${DATA}/map/map.nc
${HOME}/map/bds.sh
#endif /* !0 */

#if 0 /* Test that input data equals output data when grids are identical */
! Output equals IBIS input 1 x 1 degree grid:
bds -x 360 -y 180 -m 11 -o /tmp/zender/map/dst_1x1.nc
ncks -H -C -u -F -v tex -d lon,90.0 -d lat,40.5 /fs/cgd/data0/zender/data/soil_sfc_IBIS.nc
ncks -H -C -u -F -v mss_frc_snd,mss_frc_slt,mss_frc_cly -d lon,90.0 -d lat,40.5 /tmp/zender/map/dst_1x1.nc
#endif /* !0 */

#if 0 /* Production to CGD */
! ${HOME}/map/bds.sh produces and prunes most frequently used datasets
! bds
! bds -x 288 -y 180 -t 19801992 -o ${DATA}/map/map_19851989_0112.nc;${HOME}/map/bds.sh
! bds -x 128 -y 64 -m 23 -t 19801992 -o ${DATA}/map/map_19851989_0112.nc;${HOME}/map/bds.sh
#endif /* !0 */

#if 0 /* Production to SCD */
! WARNING: Do not overwrite these files while DEAD is using them!
! Global model triangular truncation Gaussian grids:
! Convention for grid name ordering is (x,y,z)=(lon,lat,lev)
bds -x 1440 -y 720 -m 13 -o ${DATA}/data/dst_0.25x0.25.nc
bds -x 720 -y 360 -m 13 -o ${DATA}/data/dst_0.5x0.5.nc
bds -x 360 -y 180 -m 13 -o ${DATA}/data/dst_1x1.nc
bds -x 180 -y 090 -m 13 -o ${DATA}/data/dst_2x2.nc
bds -x 120 -y 060 -m 13 -o ${DATA}/data/dst_3x3.nc
bds -x 016 -y 008 -m 23 -o ${DATA}/data/dst_T5.nc
bds -x 064 -y 032 -m 23 -o ${DATA}/data/dst_T21.nc
bds -x 096 -y 048 -m 23 -o ${DATA}/data/dst_T31.nc
bds -x 128 -y 064 -m 23 -o ${DATA}/data/dst_T42.nc
bds -x 192 -y 094 -m 23 -o ${DATA}/data/dst_T62.nc
bds -x 192 -y 096 -m 23 -o ${DATA}/data/dst_T63.nc
bds -x 256 -y 128 -m 23 -o ${DATA}/data/dst_T85.nc
bds -x 320 -y 160 -m 23 -o ${DATA}/data/dst_T106.nc
bds -x 384 -y 190 -m 23 -o ${DATA}/tmp/zender/dst_T126.nc
bds -x 512 -y 256 -m 23 -o ${DATA}/tmp/zender/dst_T170.nc

! IPCC datasets:
for mdl in 'cccma_cgcm3_1 cccma_cgcm3_1_t63 cnrm_cm3 csiro_mk3_0 \
gfdl_cm2_0 gfdl_cm2_1 giss_aom giss_model_e_h giss_model_e_r \
iap_fgoals1_0_g inmcm3_0 ipsl_cm4 miroc3_2_hires miroc3_2_medres \
miub_echo_g mpi_echam5 mri_cgcm2_3_2a ncar_ccsm3_0 ncar_pcm1 \
ukmo_hadcm3 ukmo_hadgem1'; do
  cd /data/brownmc/sresa1b/atm/mo/tas/${mdl}/run1
  ncks -H -M -m -v lat,lon,lat_bnds,lon_bnds tas_A1.nc | m
done
! NB: "c" suffix to latitude in output file name indicates Arakawa C-grid
bds -x 096 -y 048 -m 23 -o ${DATA}/data/dst_T31.nc # cccma_cgcm3_1
bds -x 128 -y 064 -m 23 -o ${DATA}/data/dst_T42.nc # cccma_cgcm3_1_t63
bds -x 128 -y 064 -m 23 -o ${DATA}/data/dst_T42.nc # cnrm_cm3
bds -x 192 -y 096 -m 23 -o ${DATA}/data/dst_T63.nc # csiro_mk3_0
bds -x 144 -y 090 -m 12 -o ${DATA}/data/dst_2.5x2.nc # gfdl_cm2_0
bds -x 144 -y 090 -m ?2 -o ${DATA}/data/dst_fxm.nc # gfdl_cm2_1 (not Gaussian!)
bds -x 090 -y 060 -m 12 -o ${DATA}/data/dst_4x3.nc # giss_aom
bds -x 072 -y 046 -m 32 -o ${DATA}/data/dst_5x4c.nc # giss_model_e_h
bds -x 072 -y 046 -m 32 -o ${DATA}/data/dst_5x4c.nc # giss_model_e_r
bds -x 128 -y 060 -m ?3 -o ${DATA}/data/dst_fxm.nc # iap_fgoals1_0_g (not Gaussian!)
bds -x 072 -y 045 -m 13 -o ${DATA}/data/dst_5x4.nc # inmcm3_0
bds -x 096 -y 072 -m ?3 -o ${DATA}/data/dst_.nc # ipsl_cm4 (not Gaussian!)
bds -x 000 -y 000 -m fxm -o ${DATA}/data/dst_fxm.nc # miroc3_2_hires
bds -x 128 -y 064 -m 23 -o ${DATA}/data/dst_T42.nc # miroc3_2_medres
bds -x 128 -y 064 -m 23 -o ${DATA}/data/dst_T42.nc # miub_echo_g
bds -x 192 -y 096 -m 23 -o ${DATA}/data/dst_T63.nc # mpi_echam5
bds -x 128 -y 064 -m 23 -o ${DATA}/data/dst_T42.nc # mri_cgcm2_3_2a
bds -x 256 -y 128 -m 23 -o ${DATA}/data/dst_T85.nc # ncar_ccsm3_0
bds -x 128 -y 064 -m 23 -o ${DATA}/data/dst_T42.nc # ncar_pcm1
bds -x 096 -y 073 -m 33 -o ${DATA}/data/dst_3.75x2.5c.nc # ukmo_hadcm3
bds -x 192 -y 145 -m 33 -o ${DATA}/data/dst_1.875x1.25c.nc # ukmo_hadgem1

! For dust version 1:
ncks -v gw,src_str_mdl,csn_cyc_mdl,mss_frc_snd,mss_frc_slt,mss_frc_cly,lnd_frc_dry,sfc_typ,lai_lsm,src_odxc ${DATA}/tmp/zender/dst_T126.nc ${DATA}/data/dst_T126.nc
! For dust version 2:
ncks -O -v gw,mss_frc_snd,mss_frc_cly,lnd_frc_dry,sfc_typ,time,mbl_bsn_fct,lai_lsm,mss_frc_CaCO3 ${DATA}/data/dst_${rsn}.nc ${DATA}/data/dst_${rsn}.nc
! For dust version 3:
ncks -O -v gw,mss_frc_snd,mss_frc_cly,lnd_frc_dry,sfc_typ,time,mbl_bsn_fct,vai_lsm,vai_ttl_clm,mss_frc_CaCO3,rdb_fct_unity,rdb_fct_tpg,rdb_fct_gmr,rdb_fct_hyd,rdb_fct_rfl_mds_lnr,rdb_fct_rfl_mds_sqr,sfc_frc_bln,sfc_acm_fct,flw_acm_fct ${DATA}/data/dst_${rsn}.nc ${DATA}/tmp/dst_${rsn}.nc
! Last command above excises all diagnostic variables not required by Dust model
! For T126 datasets, this command shrinks files from 50 to 15 MB
! NB: vai_lsm is currently computed interactively and is only accessed as 
! a placeholder for future time varying surface fields
scp ${DATA}/tmp/dst_${rsn}.nc dust.ess.uci.edu:/var/www/html/dead/data
#endif /* !0 */

#if 0 /* Erodibility factors */
! Scale erodibility to give reasonable agreement in, say, Barbados
! For topographic (GCT01) erodibility, ZBN03 and ZNT03 set flx_mss_fdg_fct in dstmbl.F90 to 7.0e-4
! Set other basin factors relative to this:
!                Absolute Relative (to Topographic):
! Topographic:     7.0e-4    1.0
! Unity:          1.13e-4    6.195
! Geomorphologic: 79.9e-4   11.414 (number-based)
! Geomorphologic:       ?    5.707 (area-based)
! Hydrologic:    157.0e-4   22.429
export rsn='T62' # Following factors are tuned for T62
ncap2 -O -s 'rdb_fct_unity=0.0*mbl_bsn_fct+0.1614' -s 'rdb_fct_tpg=mbl_bsn_fct' -s 'rdb_fct_gmr=5.707*sfc_acm_fct' -s 'rdb_fct_hyd=22.429*flw_acm_fct' -s 'rdb_fct_rfl_mds_lnr=rfl_sfc_lnd_nrm_mds_lnr' -s 'rdb_fct_rfl_mds_sqr=rfl_sfc_lnd_nrm_mds_sqr' ${DATA}/data/dst_${rsn}.nc ${DATA}/data/dst_${rsn}.nc
ncap2 -O -s 'mbl_bsn_fct=rdb_fct_gmr' ${DATA}/data/dst_${rsn}.nc ${DATA}/data/dst_${rsn}_gmr.nc
scp ${DATA}/data/dst_${rsn}.nc esmf.ess.uci.edu:/data/zender/data
scp ${DATA}/data/dst_${rsn}_gmr.nc esmf.ess.uci.edu:/data/zender/data
scp ${DATA}/data/dst_${rsn}.nc dataproc.ucar.edu:/fs/cgd/data0/zender/data
scp ${DATA}/data/dst_${rsn}_gmr.nc dataproc.ucar.edu:/fs/cgd/data0/zender/data
#endif /* !0 */

#if 0 /* Figures for prp_gtcp */
ln -s -f ${DATA}/data/dst_1x1_gtcp.nc ${DATA}/map/map_clm.nc
bds --flg_cst_tpr=F --flg_rgn_xcl=F --bkf_itr_nbr=5 -x 360 -y 180 -m 13 -o ${DATA}/data/dst_1x1_gtcp.nc
#endif /* !0 */

#if 0 /* Regressions to ZBN03 and ZNT03 */
! Aerosol source classification routines keep improving
! To reproduce geomorphic erodibility factor used in ZBN03 and ZNT03:
! 1. Use bds_rgn_xcl.txt to exclude regions north of 60N
! 2. Use bds_cst_tpr.txt to taper coastal regions globally
! 3. Generate dataset: coastal tapering and region exclusion on, no iterative backfilling
! bds --bsn_fct_hrz=F --flg_cst_tpr=T --flg_rgn_xcl=T --bkf_itr_nbr=1 -x 192 -y 94 -m 23 -o ${DATA}/data/dst_T62_ZBN03.nc
! 4. Use erodibility factor normalization procedure documented above
! 5. Test for regressions
ncdiff -O -v rdb_fct_gmr ${DATA}/data/dst_T62_ZBN03.nc ${DATA}/data/dst_T62_gmr.nc foo.nc
ncwa -O -a lat,lon -v rdb_fct_gmr foo.nc foo.nc
ncks -C -F -H -u -v rdb_fct_gmr foo.nc
#endif /* !0 */

#if 0 /* GEOS-CHEM: Sent to Duncan Fairlie 20030716 and 20040411 */
bds --bkf_itr_nbr=1 --flg_cst_tpr=T --flg_rgn_xcl=T -x 072 -y 046 -m 34 -o ${DATA}/data/dst_5x4_GSC.nc
bds --bkf_itr_nbr=1 --flg_cst_tpr=T --flg_rgn_xcl=T -x 144 -y 091 -m 34 -o ${DATA}/data/dst_2.5x2_GSC.nc
bds --bkf_itr_nbr=1 --flg_cst_tpr=T --flg_rgn_xcl=T -x 360 -y 181 -m 34 -o ${DATA}/data/dst_1x1_GSC.nc
export rsn='5x4_GSC'
export rsn='2.5x2_GSC'
export rsn='1x1_GSC'
#endif /* !0 */

#if 0 /* Import alternative surface type datasets, e.g., LGM boundary conditions */
bds -x 96 -y 48 -m 23 --flg_sfc_xtr=T --fl_sfc_xtr=${DATA}/data/sfc_typ_2x2_adams.nc -o ${DATA}/data/dst_T31_adams.nc
bds -x 96 -y 48 -m 23 --flg_sfc_xtr=T --fl_sfc_xtr=${DATA}/data/sfc_typ_2x2_crowley.nc -o ${DATA}/data/dst_T31_crowley.nc
! LGM boundary conditions at T42: Default for CCM/LGM is --flg_cst_tpr=T --flg_rgn_xcl=F
bds -x 128 -y 64 -m 23 --bkf_itr_nbr=1 --flg_cst_tpr=T --flg_rgn_xcl=F -o ${DATA}/data/dst_T42.nc
bds -x 128 -y 64 -m 23 --bkf_itr_nbr=1 --flg_cst_tpr=T --flg_rgn_xcl=F --flg_sfc_xtr=T --fl_sfc_xtr=${DATA}/data/sfc_typ_2x2_adams.nc -o ${DATA}/data/dst_T42_adams.nc
bds -x 128 -y 64 -m 23 --bkf_itr_nbr=1 --flg_cst_tpr=T --flg_rgn_xcl=F --flg_sfc_xtr=T --fl_sfc_xtr=${DATA}/data/sfc_typ_2x2_crowley.nc -o ${DATA}/data/dst_T42_crowley.nc
! Process with ZNT03 normalizations, switch to geomorphic erodibility, copy to krein:
export rsn='T42'
export rsn='T42_adams'
export rsn='T42_crowley'
#endif /* !0 */

! This program computes many diagnositic variables:
! fsh_fct(lon_nbr,lat_nbr): Source strength during peak month, sources only. A value of zero means the gridcell is not a source. fsh_fct(lon_idx,lat_idx)=src_str(lon_idx,lat_idx,time_max_idx)
! src_str(lon_nbr,lat_nbr,time_out_nbr): Mean monthly optical depth during years when gridcell is identified as a source.
! src_odxc(lon_nbr,lat_nbr,time_out_nbr): Mean monthly extinction optical depth during all years but only for sources. src_odxc = odxc except for non sources, where src_odxc=0.0.
! mth_max(lon_nbr,lat_nbr): Peak month (1--12) in annual cycle, sources only. A value of zero means the gridcell is not a source.
! src_flg(lon_nbr,lat_nbr): Source flag. 0 if gridcell is never a source in any month of any year in input data. 1 if gridcell is a source in any month of any year in input data.
! src_frq(lon_nbr,lat_nbr,time_out_nbr): Frequency that the gridcell is identified as a source in a given month. src_frq is always unity for sources in one year samples... 
! ...Will only be unity in a fifteen year sample if gridcell is a source in that month in all fifteen years.
! odxc_tms(lon_nbr,lat_nbr,time_out_nbr): Monthly mean extinction optical depth proxy at gridpoint, currently mean TOMS aerosol absorption index
! csn_cyc(lon_idx,lat_idx,time_out_idx): Monthly source strength normalized by the peak month, sources only.  csn_cyc = 1.0 during the peak month, csn_cyc = 0.0 during months with no emission from the source... 
! ...csn_cyc(lon_idx,lat_idx,time_out_idx)=src_str(lon_idx,lat_idx,time_out_idx)/src_str(lon_idx,lat_idx,time_max_idx)

!++alfgr
! Usage of changetext.py:
! changetext.py "oldtext" "newtext" indented for changing between 
! /fs/cgd/ (UCI use) and /fs/cgd/ (UiO use)

program bds
  ! Purpose: Generation of off-line datasets for use by dust model
  !          Identify and parameterize dust sources
  use bds_ctl ! [mdl] Control variables drc_in,drc_out,hgt_dlt_msl
  use bds_drv,only:bds_prc ! [mdl] BDS driver, processor, output routines
  use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
  use map_cst ! [mdl] Constants used in map routines
  use netcdf ! [mdl] netCDF interface
  use nf90_utl ! [mdl] netCDF utilities
  use sng_mdl ! [mdl] String manipulation
  use rgr_mdl,only:map_rgr
  use utl_mdl,only:date_time_get ! [mdl] Utility functions (date_time_get,mnt_chk...)
  implicit none
  ! Parameters
  character(len=*),parameter::CVS_Date="$Date$" ! [sng] Date string
  character(len=*),parameter::CVS_Header="$Id$" ! [sng] Full CVS Header
  character(len=*),parameter::CVS_Id="$Id$" ! [sng] CVS Identification
  character(len=*),parameter::CVS_Name="$HeadURL$" ! [sng] File name string
  character(len=*),parameter::CVS_Revision="$Revision$" ! [sng] File revision string
  character(len=*),parameter::nlc=char(0) ! [sng] NUL character = ASCII 0 = char(0)
  character(len=*),parameter::sbr_nm="bds" ! [sng] Subroutine name
  integer,parameter::lat_nbr_max=180 ! [nbr] Maximum number of latitudes
  integer,parameter::len_sng_max=100 ! [nbr] Maximum length of string
  integer,parameter::lon_nbr_max=288 ! [nbr] Maximum number of longitudes
  integer,parameter::time_out_nbr_max=12 ! [nbr] Number of months in a seasonal cycle
  
  ! Input
  ! Input/Output
  ! Output
  ! Local

  ! set parameters for regridding,need change for different resolution and map type
  integer,parameter::lat_rgr_in_nbr = 64
  integer,parameter::lon_rgr_in_nbr = 128
  integer,parameter::lat_rgr_out_nbr = 94
  integer,parameter::lon_rgr_out_nbr = 192
  integer,parameter::map_rgr_typ_in=map_lat_Gss_lon_Grn_ctr
  integer,parameter::map_rgr_typ_out=map_lat_Gss_lon_Grn_ctr
  real  dat_in (lon_rgr_in_nbr,lat_rgr_in_nbr)
  real  dat_out (lon_rgr_out_nbr,lat_rgr_out_nbr)

  ! Set defaults for command line options 
  character(80)::fl_in="in.nc"//nlc ! [sng] Input file
  integer::bln_nbr=4 ! [nbr] Number of tri-modal soil blends
  integer::lat_nbr=64 ! [nbr] Number of latitudes
  integer::lon_nbr=128 ! [nbr] Number of longitudes
  integer::sgs_nbr=4 ! [nbr] Number of sub-gridscale patches per gridcell
  integer::map_typ_out=map_lat_Gss_lon_Grn_ctr ! [enm] Output map grid type
  integer::yr_end=1980 ! [yr] Ending year in YYYY format
  integer::yr_srt=1980 ! [yr] Starting year in YYYY format

  ! Derived fields

  ! Locals with simple initialization and no command line override
  integer::exit_status=0 ! [enm] Program exit status
  integer::rcd=nf90_noerr ! [rcd] Return success code
  integer::time_out_nbr=12 ! [nbr] Monthly averages in a climatological seasonal cycle
  logical::GRD_MTH=.false. ! [flg] Reserved for future use

  character(26)::lcl_date_time ! [sng] Time formatted as Day Mth DD HH:MM:SS TZ YYYY
  character(80)::opt_sng ! [sng] Option string
  character(80)::arg_val ! [sng] Command line argument value
  character(500)::cmd_ln ! [sng] Command line
  character(2)::dsh_key ! [sng] Command line dash and switch
  character(80)::lbl ! [sng] Label
  character(200)::prg_ID ! [sng] Program ID
  character(80)::sng_foo ! [sng] String
  character(2)::mth_sng ! [sng] Month in MM format
  character(200)::src_fl_sng ! [sng] Source file string
  character(200)::src_mth_sng ! [sng] Source method string
  character(8)::yr_abb ! [sng] Starting and ending years in YYYYYYYY format

  integer arg_idx ! [idx] Counting index
  integer arg_nbr ! [nbr] Number of command line arguments
  integer opt_lng ! [nbr] Length of option

  integer nc_id ! [id] File handle
  integer yr_crr ! [yr] Current year in YYYY format
  integer yr_idx ! [idx] Counting index
  integer yr_nbr ! [nbr] Number of years

  ! Allocatables

  ! Main code
  dbg_lvl=dbg_off ! [enm] dbg_lvl allocated in dbg.F90
  
  ! Locals requiring complex initialization expressions
  ! Initialize defaults
  
  ! Retrieve command line arguments
  call ftn_cmd_ln_sng(cmd_ln) ! [sng] Re-construct command line into single string
  call ftn_prg_ID_mk(CVS_Id,CVS_Revision,CVS_Date,prg_ID) ! [sng] Program ID
  call date_time_get(lcl_date_time) ! Time formatted as Day Mth DD HH:MM:SS TZ YYYY  
  write (6,"(a)") prg_ID(1:ftn_strlen(prg_ID))
  arg_nbr=command_argument_count() ! [nbr] Number of command line arguments
  arg_idx=1 ! [idx] Counting index
  do while (arg_idx <= arg_nbr)
     call ftn_getarg_wrp(arg_idx,arg_val) ! [sbr] Call getarg(), increment arg_idx
     dsh_key=arg_val(1:2) ! [sng] First two characters of option
     if (dsh_key == "--") then
        opt_lng=ftn_opt_lng_get(arg_val) ! [nbr] Length of option
        if (opt_lng <= 0) stop "Long option has no name"
        opt_sng=arg_val(3:2+opt_lng) ! [sng] Option string
        if (opt_sng == "atm") then
           call ftn_arg_get(arg_idx,arg_val,caseid_atm) ! [sng] caseid for atmospheric model input
        else if (opt_sng == "bkf_itr_nbr") then ! [nbr] Number of backfill iterations (converges in 5 iterations)
           call ftn_arg_get(arg_idx,arg_val,bkf_itr_nbr)
        else if (opt_sng == "bsn_fct_hrz") then ! [flg] Compute basin factor on high-resolution topography grid
           call ftn_arg_get(arg_idx,arg_val,bsn_fct_hrz)
        else if (opt_sng == "dbg" .or. opt_sng == "dbg_lvl" ) then
           call ftn_arg_get(arg_idx,arg_val,dbg_lvl) ! [enm] Debugging level
        else if (opt_sng == "drc_in") then
           call ftn_arg_get(arg_idx,arg_val,drc_in) ! [sng] Input directory
        else if (opt_sng == "drc_out") then
           call ftn_arg_get(arg_idx,arg_val,drc_out) ! [sng] Output directory
        else if (opt_sng == "drc_sfc_xtr") then
           call ftn_arg_get(arg_idx,arg_val,drc_sfc_xtr) ! [sng] External surface type dataset directory
        else if (opt_sng == "drc_tms") then
           call ftn_arg_get(arg_idx,arg_val,drc_tms) ! [sng] TOMS input directory
        else if (opt_sng == "fl_sfc_xtr") then ! [sng] Surface type dataset, external
           call ftn_arg_get(arg_idx,arg_val,fl_sfc_xtr)
        else if (opt_sng == "flg_cst_tpr") then ! [flg] Taper erodibility factors along coasts 
           call ftn_arg_get(arg_idx,arg_val,flg_cst_tpr)
        else if (opt_sng == "flg_rdb_only") then ! [flg] Only compute basin factors
           call ftn_arg_get(arg_idx,arg_val,flg_rdb_only)
        else if (opt_sng == "flg_rgn_msk") then ! [flg] Mask region
           call ftn_arg_get(arg_idx,arg_val,flg_rgn_msk)
        else if (opt_sng == "flg_rgn_xcl") then ! [flg] Zero erodibility factors based on location
           call ftn_arg_get(arg_idx,arg_val,flg_rgn_xcl)
        else if (opt_sng == "flg_sfc_xtr") then ! [flg] Read surface type directly from external file
           call ftn_arg_get(arg_idx,arg_val,flg_sfc_xtr)
        else if (opt_sng == "flg_soi_nn") then ! [flg] Use nearest neighbor algorithm to define soil texture at gridcells where needed
           call ftn_arg_get(arg_idx,arg_val,flg_soi_nn)
        else if (opt_sng == "hgt_dlt_msl") then
           call ftn_arg_get(arg_idx,arg_val,hgt_dlt_msl) ! [m] Mean sea level change relative to current climate
        else if (opt_sng == "lat_nbr") then ! [nbr] Number of latitudes
           call ftn_arg_get(arg_idx,arg_val,lat_nbr)
        else if (opt_sng == "lnd") then ! [sng] caseid for land model input
           call ftn_arg_get(arg_idx,arg_val,caseid_lnd)
        else if (opt_sng == "lon_nbr") then ! [nbr] Number of longitudes
           call ftn_arg_get(arg_idx,arg_val,lon_nbr)
        else if (opt_sng == "tms" .or. opt_sng == "toms") then ! [sng] caseid for TOMS input
           call ftn_arg_get(arg_idx,arg_val,caseid_tms)
        else ! Option not recognized
           arg_idx=arg_idx-1 ! [idx] Counting index
           call ftn_getarg_err(arg_idx,arg_val) ! [sbr] Error handler for getarg()
        endif               ! endif option is recognized
        ! Jump to top of while loop
        cycle ! C, F77, and F90 use "continue", "goto", and "cycle"
     endif ! endif long option
     ! Handle short options
     if (dsh_key == "-D") then
        call ftn_arg_get(arg_idx,arg_val,dbg_lvl) ! [enm] Debugging level
     else if (dsh_key == "-i") then
        call ftn_arg_get(arg_idx,arg_val,fl_in) ! [sng] Input file
     else if (dsh_key == "-m") then ! [enm] Output map grid type
        call ftn_arg_get(arg_idx,arg_val,map_typ_out)
     else if (dsh_key == "-o") then
        call ftn_arg_get(arg_idx,arg_val,fl_out) ! [sng] Output file
     else if (dsh_key == "-s") then ! [frc] Source strength threshold
        call ftn_arg_get(arg_idx,arg_val,src_thr)
     else if (dsh_key == "-t") then
        call ftn_arg_get(arg_idx,arg_val,yr_abb)
        read (yr_abb,"(2(i4.4))") yr_srt,yr_end
     else if (dsh_key == "-v") then
        goto 1000
     else if (dsh_key == "-x") then ! [nbr] Number of longitudes
        call ftn_arg_get(arg_idx,arg_val,lon_nbr)
     else if (dsh_key == "-y") then ! [nbr] Number of latitudes
        call ftn_arg_get(arg_idx,arg_val,lat_nbr)
     else ! Option not recognized
        arg_idx=arg_idx-1 ! [idx] Counting index
        call ftn_getarg_err(arg_idx,arg_val) ! [sbr] Error handler for getarg()
     endif ! endif arg_val
  end do ! end while (arg_idx <= arg_nbr)
  
  ! Compute any quantities that might depend on command line input
  ! Prepend user-specified path, if any, to input data file names
  if (ftn_strlen(drc_in) > 0) then
     ! Files expected to be in map data directory have paths set here
     call ftn_drcpfx(drc_in,fl_atm) ! [sng] Atmospheric model climatology file
     call ftn_drcpfx(drc_in,fl_chm) ! [sng] Mineralogy file
     call ftn_drcpfx(drc_in,fl_grd) ! [sng] Grid file
     call ftn_drcpfx(drc_in,fl_hgt) ! [sng] Surface height dataset
     call ftn_drcpfx(drc_in,fl_in) ! [sng] Input file
     call ftn_drcpfx(drc_in,fl_lak) ! [sng] Lake file
     call ftn_drcpfx(drc_in,fl_pcp) ! [sng] Precipitation file
     call ftn_drcpfx(drc_in,fl_rfl_sfc) ! [sng] Surface reflectance file
     call ftn_drcpfx(drc_in,fl_sfc) ! [sng] Surface type dataset
     call ftn_drcpfx(drc_in,fl_soi_txt) ! [sng] Soil texture dataset
     call ftn_drcpfx(drc_in,fl_wtl) ! [sng] Wetland file
  endif ! endif drc_in > 0
  ! Prepend user-specified path, if any, to surface type dataset names
  if (ftn_strlen(drc_sfc_xtr) > 0) call ftn_drcpfx(drc_sfc_xtr,fl_sfc_xtr) ! [sng] Surface type dataset, external
  ! Prepend user-specified path, if any, to output data file names
  if (ftn_strlen(drc_out) > 0) call ftn_drcpfx(drc_out,fl_out) ! [sng] Output file
  ! Compute any quantities that might depend on command line input
  call ftn_strcpy(src_fl_sng,"Grid file is "//fl_grd)
  if (GRD_MTH) then
     call ftn_strcpy(src_mth_sng,"Method is ") ! [sng] Source method string
  else
     call ftn_strcpy(src_mth_sng,"Method is ") ! [sng] Source method string
  endif                     ! endif GRD_MTH
  yr_nbr=yr_end-yr_srt+1 ! [nbr] Number of years

  ! Begin netCDF output routines
  rcd=nf90_wrp_create(fl_out,nf90_clobber,nc_id,sbr_nm=sbr_nm)
  ! Add global attributes
  rcd=rcd+nf90_put_att(nc_id,nf90_global,"CVS_Id",CVS_Id)
  rcd=rcd+nf90_put_att(nc_id,nf90_global,"creation_date",lcl_date_time)
  rcd=rcd+nf90_put_att(nc_id,nf90_global,"prg_ID",prg_ID(1:ftn_strlen(prg_ID)))
  rcd=rcd+nf90_put_att(nc_id,nf90_global,"cmd_ln",cmd_ln(1:ftn_strlen(cmd_ln)))
  ! Now that global attributes have been defined, end define mode
  rcd=rcd+nf90_enddef(nc_id)
  rcd=nf90_wrp_close(nc_id,fl_out,"Wrote global header to",sbr_nm=sbr_nm)
  
  if (dbg_lvl /= dbg_off) then
     write(6,"(a)") "Dust Boundary Data Processor Initial State"
     write(6,"(a)") "I/O:"
     write(6,"(a,a)") "fl_atm = ",fl_atm(1:ftn_strlen(fl_atm))
     write(6,"(a,a)") "fl_grd = ",fl_grd(1:ftn_strlen(fl_grd))
     write(6,"(a,a)") "fl_out = ",fl_out(1:ftn_strlen(fl_out))
     write(6,"(a,a)") "fl_pcp = ",fl_pcp(1:ftn_strlen(fl_pcp))
     write(6,"(a,a)") "drc_in = ",drc_in(1:ftn_strlen(drc_in))
     write(6,"(a)") "Output grid:"
     write(6,"(a,i3)") "bln_nbr = ",bln_nbr
     write(6,"(a,i3)") "lat_nbr = ",lat_nbr
     write(6,"(a,i3)") "lon_nbr = ",lon_nbr
     write(6,"(a,i3)") "sgs_nbr = ",sgs_nbr
     write(6,"(a,i2)") "time_out_nbr = ",time_out_nbr
     write(6,"(a,i2)") "map_typ_out = ",map_typ_out
     write(6,"(a)") "Data:"
     write(6,"(a,i4)") "yr_srt = ",yr_srt
     write(6,"(a,i4)") "yr_end = ",yr_end
     write(6,"(a,f4.3)") "src_thr = ",src_thr
  endif ! endif dbg
  
  call bds_prc( &
       bln_nbr, & ! [nbr] Number of tri-modal soil blends
       lat_nbr, & ! [nbr] Number of latitudes
       lon_nbr, & ! [nbr] Number of longitudes
       map_typ_out, & ! [enm] Output map grid type
       sgs_nbr, & ! [nbr] Number of sub-gridscale patches per gridcell
       time_out_nbr, & ! [nbr] Number of times
       yr_nbr, & ! [nbr] Number of years
       yr_srt) ! Starting year in YYYY format
  
1000 continue

  if(.not.flg_rgr) then
     ! Interpolate data to different grid
     call map_rgr( &
          dat_in,         & ! I Input file 
          map_rgr_typ_in, & ! I Input map grid type
          lon_rgr_in_nbr, & ! I Number of longitudes
          lat_rgr_in_nbr, & ! I Number of latitudes
          map_rgr_typ_out, & ! I Output map grid typ
          lon_rgr_out_nbr, & ! I Number of longitudes
          lat_rgr_out_nbr, & ! I Number of latitudes
          dat_out) ! O output data
     endif    ! endif flg_rgr
  
  call exit(exit_status)
end program bds                       ! end program bds()
