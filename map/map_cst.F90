! $Id$ -*-f90-*-

! Purpose: Stores constants used in map routines

! Usage:
!use map_cst ! [mdl] Constants used in map routines

module map_cst ! [mdl] Constants used in map routines
  implicit none
  
  real(selected_real_kind(p=12)),parameter::earth_rds=6.370e+06 ! [m] Earth radius
  real(selected_real_kind(p=12)),parameter::earth_rds_sqr=earth_rds*earth_rds ! [m2] Earth radius squared
  integer,parameter::map_lon_srt_Grn=1 ! [enm] Longitude arrays begin at Greenwich
  integer,parameter::map_lon_srt_180=2 ! [enm] Longitude arrays begin at date line
  integer,parameter::map_lat_grd_rgl=3 ! [enm] Latitude grid is regularly spaced
  integer,parameter::map_lat_grd_Gss=4 ! [enm] Latitude grid is Gaussian
  integer,parameter::map_lon_ctr_wst=5 ! [enm] West edge of first gridcell is on Greenwich/Date line
  integer,parameter::map_lon_ctr_ctr=6 ! [enm] Center first gridcell is on Greenwich/Date line
  integer,parameter::map_lat_grd_GSC=7 ! [enm] Latitude grid is regularly spaced except at Poles (GEOS-CHEM grid)
  integer,parameter::map_lat_grd_cos=8 ! [enm] Latitude grid is cosine-weighted (equal-area)
  integer,parameter::map_lat_grd_ease_ml=9 ! [enm] Latitude grid is EASE ML
  integer,parameter::map_lat_ncr=10  ! [enm] Latitude grid increases monotonically (-90 -> +90) (unused?)
  integer,parameter::map_lat_dcr=11  ! [enm] Latitude grid decreases monotonically (+90 -> -90) (unused?)

  ! Summary types:
  ! map_lat_Gss_lon_Grn_ctr: CCM 1--3, CAM1--3, LSM, MATCH, lnd_frc_flak.1x1, lnd_frc_swmp.1x1
  ! map_lat_rgl_lon_Grn_wst: hgt_sfc_topo.nc
  ! map_lat_rgl_lon_180_wst: sfc_typ_olson.data, sfc_typ_ds769.txt, soi_txt_webb.1x1, soil_sfc_IBIS.nc, IGBP-DIS, CIESIN/SEDAC
  ! map_lat_GSC_lon_Grn_ctr: UKMO, CAM FV
  ! map_lat_GSC_lon_180_ctr: UCICTM (4x5), GEOS-CHEM
  ! map_lat_rgl_lon_Grn_ctr: TOMS AAI
  ! map_lat_cos_lon_Grn_ctr: 
  ! map_lat_Gss_lon_180_ctr: UCICTM (T42)
  ! map_lat_EML_lon_180_ctr: AMSR-E (http://nsidc.org/data/ease/ease_grid.html)

  ! The longitude type is determined by the answer to the question:
  ! Where is Greenwich or the date line in terms of the first longitude cell?
  ! Only valid answers are: [West edge]/[Center] of first longitude cell
  integer,parameter::map_lat_rgl_lon_180_wst=11 ! [enm] Latitudes are regular, Date line at West edge of first longitude bin
  integer,parameter::map_lat_rgl_lon_Grn_wst=12 ! [enm] Latitudes are regular, Greenwich at West edge of first longitude bin
  integer,parameter::map_lat_rgl_lon_Grn_ctr=13 ! [enm] Latitudes are regular, Greenwich at center of first longitude bin
  integer,parameter::map_lat_rgl_lon_180_ctr=14 ! [enm] Latitudes are regular, Date line at center of first longitude bin
  integer,parameter::map_lat_Gss_lon_180_wst=21 ! [enm] Latitudes are Gaussian, Date line at West edge of first longitude bin
  integer,parameter::map_lat_Gss_lon_Grn_wst=22 ! [enm] Latitudes are Gaussian, Greenwich at West edge of first longitude bin
  integer,parameter::map_lat_Gss_lon_Grn_ctr=23 ! [enm] Latitudes are Gaussian, Greenwich at center of first longitude bin
  integer,parameter::map_lat_Gss_lon_180_ctr=24 ! [enm] Latitudes are Gaussian, Date line at center of first longitude bin
  integer,parameter::map_lat_GSC_lon_180_wst=31 ! [enm] Latitudes are GEOS-CHEM, Date line at West edge of first longitude bin
  integer,parameter::map_lat_GSC_lon_Grn_wst=32 ! [enm] Latitudes are GEOS-CHEM, Greenwich at West edge of first longitude bin
  integer,parameter::map_lat_GSC_lon_Grn_ctr=33 ! [enm] Latitudes are GEOS-CHEM, Greenwich at center of first longitude bin
  integer,parameter::map_lat_GSC_lon_180_ctr=34 ! [enm] Latitudes are GEOS-CHEM, Date line at center of first longitude bin
  integer,parameter::map_lat_EML_lon_180_wst=41 ! [enm] Latitudes are EASE-ML, Date line at West edge of first longitude bin
  integer,parameter::map_lat_EML_lon_Grn_wst=42 ! [enm] Latitudes are EASE-ML, Greenwich at West edge of first longitude bin
  integer,parameter::map_lat_EML_lon_Grn_ctr=43 ! [enm] Latitudes are EASE-ML, Greenwich at center of first longitude bin
  integer,parameter::map_lat_EML_lon_180_ctr=44 ! [enm] Latitudes are EASE-ML, Date line at center of first longitude bin
  
end module map_cst

