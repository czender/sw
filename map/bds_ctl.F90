! $Id$ -*-f90-*- 

! Purpose: Control variables for Boundary DataSet processor bds
! Separate this module from bds_drv to avoid circular dependency
! build problems

! Usage:
!use bds_ctl ! [mdl] Control variables drc_in,drc_out,hgt_dlt_msl

module bds_ctl ! [mdl] Control variables drc_in,drc_out,hgt_dlt_msl
  implicit none
  public ! [stt] Symbols are public unless individually qualified as private

  character(80)::caseid_atm="dstccm25"//char(0) ! [sng] caseid for atmospheric model input
  character(80)::caseid_lnd="lsm"//char(0) ! [sng] caseid for land model input
  character(80)::caseid_tms="toms"//char(0) ! [sng] caseid for TOMS input
  character(80)::drc_in="/data/zender/map"//char(0) ! [sng] Input directory
  character(80)::drc_sfc_xtr="" ! [sng] External surface type dataset directory
  character(80)::drc_out="" ! [sng] Output directory
  character(80)::drc_tms="/data/zender/map"//char(0) ! [sng] TOMS input directory
  ! character(80)::drc_tms="/data/zender/toms"//char(0) ! [sng] TOMS input directory
  ! character(80)::fl_hgt="hgt_sfc_DS754.0"//char(0) ! [sng] Surface height dataset
  character(80)::fl_atm="dstccm25_clm_0112.nc"//char(0) ! [sng] Atmospheric model climatology file
  character(80)::fl_chm="igbp_"//char(0) ! [sng] Mineralogy file
  character(80)::fl_grd="toms_198001.nc"//char(0) ! [sng] Grid file
  character(80)::fl_hgt="hgt_sfc_ngdc_tbase_gfdl.nc"//char(0) ! [sng] Surface height dataset
  character(80)::fl_hyd="clmigbp_clm_0112.nc"//char(0) ! [sng] Surface hydrology dataset
  ! character(80)::fl_lai="/datashare/inputdata/csm/lnd/clm2/rawdata/mksrf_lai.nc"//char(0) ! [sng] Leaf Area Index dataset
  character(80)::fl_lai="/data/zender/csm/inputdata/lnd/clm2/rawdata/mksrf_lai.nc"//char(0) ! [sng] Leaf Area Index dataset
  character(80)::fl_lak="lnd_frc_flak.1x1"//char(0) ! [sng] Lake file
  character(80)::fl_out="/tmp/zender/map/map.nc"//char(0) ! [sng] Output file
  character(80)::fl_rfl_sfc="rfl_sfc_MOD09_B7_2001_annual.nc"//char(0) ! [sng] Surface reflectance file
  character(80)::fl_pcp="pcp_cmap_7901_9612.nc"//char(0) ! [sng] Precipitation file
  ! character(80)::fl_pft="/datashare/inputdata/csm/lnd/clm2/rawdata/mksrf_pft.nc"//char(0) ! [sng] Plant Functional Type dataset
  character(80)::fl_pft="/data/zender/csm/inputdata/lnd/clm2/rawdata/mksrf_pft.nc"//char(0) ! [sng] Plant Functional Type dataset
  character(80)::fl_rgn_msk="bds_rgn_msk.txt"//char(0) ! [sng] Mask region
  character(80)::fl_rgn_xcl="bds_rgn_xcl.txt"//char(0) ! [sng] Surface regions to exclude as dust sources
  character(80)::fl_sfc="sfc_typ_olson.data"//char(0) ! [sng] Surface type dataset
  character(80)::fl_sfc_xtr="/data/zender/data/dst_T31_adams.nc"//char(0) ! [sng] Surface type dataset, external
  !  character(80)::fl_soi_txt="soi_txt_webb.1x1"//char(0) ! [sng] Soil Texture
  character(80)::fl_soi_txt="soi_txt_IGBP.1x1"//char(0) ! [sng] Soil Texture
  !  character(80)::fl_soi_txt="/data/mflanner/adams_sea/IGBP_adams_sea.1x1"//char(0) ! [sng] Soil Texture
  character(80)::fl_tpr="bds_cst_tpr.txt"//char(0) ! [sng] Coastal regions to taper as dust sources
  character(80)::fl_wtl="lnd_frc_swmp.1x1"//char(0) ! [sng] Wetland file

  ! character(80)::fl_sfc="/home/zender/dst/sfc_typ_ds769.txt"//char(0) ! [sng] Surface type dataset
  integer::bkf_itr_nbr=1 ! [nbr] Number of backfill iterations (converges in 5 iterations)
  logical::bsn_fct_hrz=.false. ! [flg] Compute basin factor on high-resolution topographic grid
  logical::flg_rdb_only=.false. ! [flg] Only compute basin factors
  logical::flg_cst_tpr=.true. ! [flg] Taper erodibility factors along coasts
  logical::flg_rgn_msk=.true. ! [flg] Mask region
  logical::flg_rgn_xcl=.true. ! [flg] Zero erodibility factors based on location
  logical::flg_sfc_xtr=.false. ! [flg] Read surface type directly from external file
  logical::flg_soi_nn=.true. ! [flg] Use nearest neighbor algorithm to define soil texture where needed
  logical::flg_rgr=.true. ![flg] Use regridding
  real::hgt_dlt_msl=0.0 ! [m] Mean sea level change relative to current climate
  real::src_thr=0.9 ! [frc] Threshold for source maxima

end module bds_ctl ! [mdl] Control variables drc_in,drc_out,hgt_dlt_msl
