// $Id$ 

// Purpose: Boundary layer meteorology utilities for C++ programs

/* Copyright (C) 1997--present Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

#include <blm.hh> // Boundary layer meteorology

int // O [rcd] Return success code
lnd_frc_mbl_get
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp doy, // I [day] Day of year [1.0..367.0)
 const prc_cmp lat_rdn, // I [rdn] Latitude
 const prc_cmp *lnd_frc_dry, // I [frc] Dry land fraction
 prc_cmp *lnd_frc_mbl, // O [frc] Bare ground fraction
 const prc_cmp *oro, // I [frc] Orography
 const long *sfc_typ, // I [idx] LSM surface type (0..28)
 const prc_cmp *snw_frc) // I [frc] Fraction of surface covered by snow
{
  /* Purpose: Return fraction of each gridcell suitable for dust mobilization
     lnd_frc_mbl_get() is called by dst_mbl()
     The date is used to obtained the time-varying vegetation cover
     Routine currently computes latitude slice of LAI based on surface type, date, and latitude
     LAI algorithm is from CCM:lsm/phenol() Bon96
     The LSM data are mid-month values, i.e., valid on the 15th of the month */
  // Requires <lsm.hh>,<cln.hh>,<cmath>,<cassert>
  using lsm::pnt_typ; // [idx] LSM plant type (1..14 = nbr_LSM_pnt_typ)
  using lsm::pnt_frc; // [frc] Weight of corresponding plant type (sums to 1.0)
  using lsm::tai; // [m2 m-2] Monthly leaf area index + stem area index, one-sided
  using cln::dpy_365; // (365) [nbr] Model days per year
  using cln::ldom_doy_365; // [day] Day of year of last day of month
  using cln::mpy; // (12) [nbr] Model months per year
  using cln::fdom_doy_365; // [day] Day of year of first day of month 

  using mth::cst_M_PIl; // (3.1415926535897932384626433832795029L) [frc] 3
  const prc_cmp pi(static_cast<prc_cmp>(mth::cst_M_PIl)); // [frc] 3
  const prc_cmp vai_mbl_thr(0.2); // [m2 m-2] VAI threshold quenching dust mobilization
  // Local
  prc_cmp day_grw; // [day] Days since Jan 1 in NH or Jul 1 in SH (1.0..365.0)
  prc_cmp dom; // [day] Current day of month [1.0..32.0)
  prc_cmp lat_dgr; // [dgr] Latitude
  prc_cmp mth_grw; // [mth] Months since Jan (NH) or Jul (SH) (1.0..12.0)
  prc_cmp pnt_frc_mbl; // [frc] "Bare ground" fraction of sub-gridscale cell
  prc_cmp pnt_frc_sgs; // [frc] Plant fraction of current sub-gridscale cell
  prc_cmp vai_sgs; // [m2 m-2] Leaf + stem area index, one-sided
  prc_cmp wgt_glb; // [frc] Interpolation weight
  prc_cmp wgt_lub; // [frc] Interpolation weight
  int rcd(0); // O [rcd] Return success code
  long iday_NH; // [day] Day number Northern Hemisphere [1..365]
  long iday_SH; // [day] iday_NH shifted 6 mth for SH [1..365]
  long idoy; // [day] Current day of year [1..365]
  long idx_idx; // [idx] Counting index
  long idx_mth_glb; // [idx] Interpolation month, future
  long idx_mth_lub; // [idx] Interpolation month, past
  long imoy; // [mth] Current month of year [1..12]
  std::valarray<long int> lnd_idx(lon_nbr); // [idx] Longitude index array (land)
  long lnd_nbr(0L); // [nbr] Number of land points
  long lon_idx; // [idx] Counting index for longitude
  long pnt_typ_idx; // [idx] Plant type index
  long sfc_typ_idx; // [idx] Surface type index
  long sgs_idx; // [idx] Surface sub-gridscale index
  
  // Main Code
  // Sanity check
  assert(doy >= 1.0 && doy < dpy_365+1);
  assert(vai_mbl_thr > 0.0);
  // Initialize defaults
  lat_dgr=180.0*lat_rdn/pi; // [dgr] Latitude
  vec_set(lnd_frc_mbl,lon_nbr,0.0); // [frc]
  
  // Land ahoy!
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    if(oro_is_lnd(oro[lon_idx])){
      lnd_idx[lnd_nbr]=lon_idx;
      lnd_nbr++;
    }  // endif
  }  // end loop over lon
  
  // Much ado about nothing
  if(lnd_nbr == 0) return rcd;

  // Convert doy to moy and dom
  imoy=1;
  idoy=static_cast<long>(doy);
  while(ldom_doy_365[imoy] < idoy && imoy < mpy) imoy++;
  dom=doy-fdom_doy_365[imoy]+1.0;
  // Sanity check
  assert(imoy >= 1 && imoy <= mpy);
  assert(dom >= 1.0 && dom <= 32.0);
  /* Compute iday_NH and iday_SH:
     iday_NH = day of current year since Jan 0 [1..365]
     iday_SH = iday_NH shifted 6 mth for SH: 1->183, 183->365, 184->1, 365->182 */
  iday_NH=static_cast<long>(doy); // Day of year [1..365]
  iday_SH=1L+((iday_NH-1+dpy_365/2)%dpy_365);
  // Compute cyclic indices into monthly data idx_mth_glb,idx_mth_lub
  if(lat_dgr >= 0.0) day_grw=iday_NH; else day_grw=iday_SH;
  mth_grw=12.0*(day_grw-0.5)/dpy_365;
  idx_mth_glb=static_cast<long>(mth_grw+0.5);
  idx_mth_lub=idx_mth_glb+1;
  wgt_glb=(idx_mth_glb+0.5)-mth_grw;
  wgt_lub=1.0-wgt_glb;
  if(idx_mth_glb < 1) idx_mth_glb=mpy;
  if(idx_mth_lub > mpy) idx_mth_lub=1;

  // Land points
  for(idx_idx=0;idx_idx<lnd_nbr;idx_idx++){
    lon_idx=lnd_idx[idx_idx];
    // Store surface blend of current gridpoint
    sfc_typ_idx=sfc_typ[lon_idx];
    if(sfc_typ_idx <= 1 || sfc_typ_idx >= 27){ // Inland lakes and land ice or wetlands
      lnd_frc_mbl[lon_idx]=0.0; // [frc] Bare ground fraction
    }else{ // Normal land
      for(sgs_idx=0;sgs_idx<3;sgs_idx++){
	// Non-vegetated surfaces have lsm::pnt_typ=14
	pnt_typ_idx=lsm::pnt_typ[sgs_idx][sfc_typ_idx]; // NB: C pnt_typ is transposed from Fortran 
	pnt_frc_sgs=pnt_frc[sgs_idx][sfc_typ_idx]; // [frc]
	// Current LAI is weighted past LAI...
	vai_sgs= // [m2 m-2] Leaf + stem area index, one-sided
	  +pnt_frc_sgs*tai[pnt_typ_idx][idx_mth_glb]*wgt_glb
	  +pnt_frc_sgs*tai[pnt_typ_idx][idx_mth_lub]*wgt_lub;
	// The "bare ground" fraction of the current sub-gridscale cell decreases
	// linearly from 1.0 to 0.0 as VAI increases from 0.0 to vai_mbl_thr
	pnt_frc_mbl= // [frc] "Bare ground" fraction of sub-gridscale cell
	  1.0-min_cpv(vai_sgs,vai_mbl_thr)/vai_mbl_thr;
	if(dbg_lvl_get() == dbg_old){
	  std::cerr << "lnd_frc_mbl_get(): vai_sgs = " << vai_sgs << " m-2 m-2, pnt_frc = " << pnt_frc_sgs << " , pnt_frc_mbl = " << pnt_frc_mbl << std::endl;
	} // endif dbg
	lnd_frc_mbl[lon_idx]+=pnt_frc_sgs*pnt_frc_mbl; // [frc]
      } // end loop over number of possible plant types
    } // endif normal land

    // Adjust for factors which constrin entire gridcell
    lnd_frc_mbl[lon_idx]=  
      lnd_frc_mbl[lon_idx] // [frc] Bare ground fraction
      *lnd_frc_dry[lon_idx] // Account for wetlands, lake, ocean, ice
      *(1.0-snw_frc[lon_idx]); // Account for snow coverage
    assert(lnd_frc_mbl[lon_idx] <= 1.0 && lnd_frc_mbl[lon_idx] >= 0.0);
  } // end loop over land
  return rcd;
} // end lnd_frc_mbl_get()

std::string soi_typ_sng_get(const long soi_typ_idx)
{
  // Purpose: Convert LSM soil type index into descriptive string
  using lsm::nbr_LSM_soi_typ; // dry land, glacier, deep lake, shallow lake, wetland
  // Local
  const std::string soi_typ_sng[nbr_LSM_soi_typ+1]= // 
  {"Undefined (ERROR---should not be used)", // soi_typ 0
   "Dry land", // soi_typ 1
   "Land ice", // soi_typ 2
   "Deep lake", // soi_typ 3
   "Shallow lake", // soi_typ 4
   "Wetland" // soi_typ 5
  }; // end soi_typ_sng
  return soi_typ_sng[soi_typ_idx];
} // end soi_typ_sng_get()

std::string oro_sng_get(const prc_cmp oro)
{
  /* Purpose: Convert orography value into descriptive string
     Requires access to round() template in mth.hh */
  // Local
  const std::string oro_sng[3]= // 
  {"Ocean", // oro 0
   "Land", // oro 1
   "Sea ice" // oro 2
  }; // end oro_sng
  return oro_sng[mth_fnc::round(oro)];
} // end oro_sng_get()

std::string pnt_typ_sng_get(const long pnt_typ_idx)
{
  // Purpose: Convert LSM plant type index into descriptive string
  using lsm::nbr_LSM_pnt_typ; // [nbr] Number of LSM plant types
  // Local
  const std::string pnt_typ_sng[nbr_LSM_pnt_typ+1]= // 
  {"Ocean (ERROR---should not be used)", // pnt_typ  0
   "Needleleaf evergreen tree", // pnt_typ  1
   "Needleleaf deciduous tree", // pnt_typ  2
   "Broadleaf evergreen tree", // pnt_typ  3
   "Broadleaf deciduous tree", // pnt_typ  4
   "Tropical seasonal tree", // pnt_typ  5
   "Cool grass (C3)", // pnt_typ  6
   "Evergreen shrub", // pnt_typ  7
   "Deciduous shrub", // pnt_typ  8
   "Arctic deciduous shrub", // pnt_typ  9
   "Arctic grass", // pnt_typ 10
   "Crop", // pnt_typ 11
   "Irrigated crop", // pnt_typ 12
   "Warm grass (C4)", // pnt_typ 13
   "Not vegetated" // pnt_typ 14
  }; // end pnt_typ_sng
  return pnt_typ_sng[pnt_typ_idx];
} // end pnt_typ_sng_get()

std::string pft_sng_get(const long pft_idx)
{
  // Purpose: Convert CLM plant functional type index into descriptive string
  using lsm::pft_nbr_CLM; // [nbr] Number of CLM plant functional types
  // Local
  const std::string pft_sng[pft_nbr_CLM+1]= // 
  {"Not vegetated (ERROR---should not be used)", // pft  0
   "Needleleaf evergreen temperate tree", // pft  1
   "Needleleaf evergreen boreal tree", // pft  2
   "Needleleaf deciduous boreal tree", // pft  3
   "Broadleaf evergreen tropical tree", // pft  4
   "Broadleaf evergreen temperate tree", // pft  5
   "Broadleaf deciduous tropical tree", // pft  6
   "Broadleaf deciduous temperate tree", // pft  7
   "Broadleaf deciduous boreal tree", // pft  8
   "Broadleaf evergreen shrub", // pft  9
   "Broadleaf deciduous temperate shrub", // pft 10
   "Broadleaf deciduous boreal shrub", // pft 11
   "C3 Arctic grass", // pft 12
   "C3 non-Arctic grass", // pft 13
   "C4 grass", // pft 14
   "Corn", // pft 15
   "Wheat" // pft 16
  }; // end pft_sng
  return pft_sng[pft_idx];
} // end pft_sng_get()

std::string sfc_typ_dsc_get(const long sfc_typ_idx)
{
  // Purpose: Construct a string describing the given surface type
  // Requires: <blm.hh>,<lsm.hh>,<cstdio>
  using lsm::nbr_LSM_sgs_sfc; // [nbr] Number of LSM sub-gridscale patches
  using lsm::pnt_frc; // [frc] Weight of corresponding plant type (sums to 1.0)
  using lsm::pnt_typ; // [idx] LSM plant type (1..14 = nbr_LSM_pnt_typ)

  long sgs_idx;
  long sgs_nbr(0L); // Number of subgrid plant types
  char c_sng[400]; // Output string, C format
  std::string sng; // Output string
  
  while(pnt_frc[sgs_nbr][sfc_typ_idx] != 0.0 && sgs_nbr < nbr_LSM_sgs_sfc) sgs_nbr++;
  assert(sgs_nbr <= nbr_LSM_sgs_sfc);
  (void)std::sprintf(c_sng,"Surface type %ld is %s comprised of %ld sub-grid plant types: ",sfc_typ_idx,sfc_typ_sng_get(sfc_typ_idx).c_str(),sgs_nbr); 
  sng=c_sng;
  
  for(sgs_idx=0;sgs_idx<sgs_nbr;sgs_idx++){
    (void)std::sprintf(c_sng,"%g%% %s (plant type %ld)",pnt_frc[sgs_idx][sfc_typ_idx]*100.0,pnt_typ_sng_get(pnt_typ[sgs_idx][sfc_typ_idx]).c_str(),pnt_typ[sgs_idx][sfc_typ_idx]);
    sng+=c_sng;
    if(sgs_idx < sgs_nbr-1) sng+=", "; 
  } // end loop over sgs_idx
  
  return sng;
} // end sfc_typ_dsc_get()

std::string sfc_typ_sng_get(const long sfc_typ_idx)
{
  // Purpose: Convert LSM surface type index into descriptive string
  using lsm::nbr_LSM_sfc_typ; // [nbr] Number of LSM surface types
  // Local
  const std::string sfc_typ_sng[nbr_LSM_sfc_typ]= // 
  {"Ocean", // sfc_typ 0
   "Land ice (glacier)", // sfc_typ 1
   "Desert", // sfc_typ 2
   "Cool needleleaf evergreen forest", // sfc_typ 3
   "Cool needleleaf deciduous forest", // sfc_typ 4
   "Cool broadleaf deciduous forest", // sfc_typ 5
   "Cool mixed needleleaf evergreen and broadleaf deciduous tree", // sfc_typ 6
   "Warm needleleaf evergreen forest", // sfc_typ 7
   "Warm broadleaf deciduous forest", // sfc_typ 8
   "Warm mixed needleleaf evergreen and broadleaf deciduous tree", // sfc_typ 9
   "Tropical broadleaf evergreen forest", // sfc_typ 10
   "Tropical seasonal deciduous forest", // sfc_typ 11
   "Savanna", // sfc_typ 12
   "Evergreen forest tundra", // sfc_typ 13
   "Deciduous forest tundra", // sfc_typ 14
   "Cool forest crop", // sfc_typ 15
   "Warm forest crop", // sfc_typ 16
   "Cool grassland", // sfc_typ 17
   "Warm grassland", // sfc_typ 18
   "Tundra", // sfc_typ 19
   "Evergreen shrubland", // sfc_typ 20
   "Deciduous shrubland", // sfc_typ 21
   "Semi-desert", // sfc_typ 22
   "Cool irrigated crop", // sfc_typ 23
   "Cool crop", // sfc_typ 24
   "Warm irrigated crop", // sfc_typ 25
   "Warm crop", // sfc_typ 26
   "Forest wetland (mangrove)", // sfc_typ 27
   "Non-forest wetland" // sfc_typ 28
  }; // end sfc_typ_sng
  return sfc_typ_sng[sfc_typ_idx];
} // end sfc_typ_sng_get()

int // Return success code
lnd_is_vgt
(const long lon_nbr, // I [nbr] Size of arrays
 const long *pnt_typ_idx, // I [idx] Plant type index 
 const prc_cmp *snw_hgt, // I [m] Geometric bulk thickness of snow
 bool *vgt) // O [flg] "Vegetated" flag
{
  // Purpose: Determine whether surface cover is "vegetated" for purposes of surface temperature calculation
  /* Notes on vegetation:
     Bon96 sets a flag called veg(k) to determine whether a surface is "vegetated"
     veg(k) is false for bare ground, and for vegetated surfaces in months when LAI+SAI = 0, otherwise veg(k) is true
     veg(k) is also false when vegetation is covered by snow
     Thus, plant type index alone is insufficient for determining vegetation flag
  */
  // Currently only uses pnt_typ_idx criterion 
  // Taken from Bon96
  // Dependencies: <lsm.hh>
  // Local
  int rcd(0); // O [rcd] Return success code
  long lon_idx; // Counting index for lon
  
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    if(pnt_typ_idx[lon_idx] != 14) vgt[lon_idx]=true; else vgt[lon_idx]=false;
    // fxm: Add in time interpolation of lai, sai, and account
    // for burying by snow
    //    vgt_area_frc=elai[pnt_typ_idx[lon_idx]]+esai[pnt_typ_idx[lon_idx]]
    if(vgt[lon_idx]){
      if(snw_hgt[lon_idx] > 0.0){
      } // endif
    } // endif vegetated
  }  // end loop over lon
  return rcd;
} // end lnd_is_vgt()

int mmn_dlt_evl
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *hgt_lwr, // I [m] Height of lower level
 const prc_cmp *hgt_upr, // I [m] Height of upper level
 const prc_cmp *hgt_zpd, // I [m] Zero plane displacement
 prc_cmp *mmn_dlt, // O [m s-1] Upper level - lower level windspeed
 const prc_cmp *mno_lng, // I [m] Monin-Obukhov length
 const prc_cmp *wnd_frc) // I [m s-1] Surface friction velocity
{
  // Purpose: Compute difference between wind speed at upper and lower heights
  // Taken from Bon96
  // Dependencies: phys_cst.hh
  using phc::cst_von_krm_rcp; // (2.5) [frc] Reciprocal of Von Karman's constant
  // Local
  a2d_cls<prc_cmp> mno_stb_crc_mmn(lon_nbr,2); // [frc] Monin-Obukhov stability correction term
  a2d_cls<prc_cmp> mno_stb_prm(lon_nbr,2); // [frc] Monin-Obukhov stability parameter 
  prc_cmp sml_fnc_mmn_uns_rcp; // Reciprocal of similarity function for momentum, unstable atmosphere
  prc_cmp tmp2; // Term in stability correction computation
  prc_cmp tmp3; // Term in stability correction computation
  prc_cmp tmp4; // Term in stability correction computation
  std::valarray<prc_cmp> wnd_crc_fct(lon_nbr); // [frc] Wind correction factor
  int rcd(0); // O [rcd] Return success code
  long idx_idx; // Counting index for lon_idx
  long lon_idx; // Counting index for lon
  long lvl_idx; // Stability computation loop index
  // Code uses notation of Bon96 p. 50, where lvl_idx=0 is lower height, lvl_idx=1 is upper height
  const long hgt_lwr_idx(0L); // Named index for lower height
  const long hgt_upr_idx(1); // Named index for upper height
  
  // Initialize some variables
  vec_set(mmn_dlt,lon_nbr,0.0); // [m s-1]
  
  std::valarray<long> vld_idx(lon_nbr); // Valid indices
  long vld_nbr; // Number of valid indices
  whenfvlt(lon_nbr,hgt_zpd,hgt_lwr,&vld_idx[0],vld_nbr);
  
  // Compute horizontal wind speed at reference height
  for(idx_idx=0;idx_idx<vld_nbr;idx_idx++){
    lon_idx=vld_idx[idx_idx];
    // Stability functions only valid for stability parameter zeta <= 1.0 Bon96 p. 52, Bru82
    mno_stb_prm(lon_idx,hgt_lwr_idx)=min_cpv((hgt_lwr[lon_idx]-hgt_zpd[lon_idx])/mno_lng[lon_idx],1.0); // [frc]
    mno_stb_prm(lon_idx,hgt_upr_idx)=min_cpv((hgt_upr[lon_idx]-hgt_zpd[lon_idx])/mno_lng[lon_idx],1.0); // [frc]
    for(lvl_idx=0;lvl_idx<2;lvl_idx++){
      if(mno_stb_prm(lon_idx,lvl_idx) < 0.0){
	sml_fnc_mmn_uns_rcp=std::pow(1.0-16.0*mno_stb_prm(lon_idx,lvl_idx),0.25);
	tmp2=std::log((1.0+sml_fnc_mmn_uns_rcp*sml_fnc_mmn_uns_rcp)/2.0);
	tmp3=std::log((1.0+sml_fnc_mmn_uns_rcp)/2.0);
	mno_stb_crc_mmn(lon_idx,lvl_idx)=2.0*tmp3+tmp2-2.0*std::atan(sml_fnc_mmn_uns_rcp)+1.5707963; // [frc]
      }else{ // not stable
	mno_stb_crc_mmn(lon_idx,lvl_idx)=-5.0*mno_stb_prm(lon_idx,lvl_idx); // [frc]
      } // not stable
    } // end loop over lvl_idx
    tmp4=std::log((hgt_upr[lon_idx]-hgt_zpd[lon_idx])/(hgt_lwr[lon_idx]-hgt_zpd[lon_idx]));
    wnd_crc_fct[lon_idx]=tmp4-mno_stb_crc_mmn(lon_idx,hgt_upr_idx)+mno_stb_crc_mmn(lon_idx,hgt_lwr_idx); // [frc]
    // wnd_crc_fct[lon_idx]=tmp4; // Neutral stability assumption
    mmn_dlt[lon_idx]=wnd_frc[lon_idx]*cst_von_krm_rcp*wnd_crc_fct[lon_idx]; // [m s-1]
  } // end loop over lon
  
  return rcd;
} // end mmn_dlt_get()

int rss_aer_get
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *hgt_mdp, // I [m] Midlayer height above surface
 const prc_cmp *hgt_zpd, // I [m] Zero plane displacement height
 const prc_cmp *mno_lng, // I [m] Monin-Obukhov length
 const prc_cmp *rgh_mmn, // I [m] Roughness length momentum
 prc_cmp *rss_aer, // O [s m-1] Aerodynamic resistance
 const prc_cmp *wnd_frc) // I [m s-1] Surface friction velocity
{
  // Purpose: Given the surface layer structure and dynamics properties, compute and return the aerodynamic resistance 
  // Dependencies: <phys_cst.hh>
  using phc::cst_von_krm; // (0.4) [frc] Von Karman's constant

  /* Notes on algorithm:
     Routine is written using nomenclature of SeP97 p. 963
     SeP97 obtains stability correction to resistance by subtracting correction
     at roughness height from correction at midlayer height.
     Bon96 assumes correction at roughness height is near zero,
     which is a good assumption for z0m << 1.0, z0m ~ 0.0.
  */
  // Local
  int rcd(0); // O [rcd] Return success code
  long lon_idx; // Counting index for lon
  prc_cmp eta_sqr_rlm; // Eta squared at roughness height
  prc_cmp eta_sqr_gcm; // Eta squared at GCM layer height
  prc_cmp eta_rlm; // Eta at roughness height
  prc_cmp eta_gcm; // Eta at GCM layer height
  prc_cmp nmr_rlm; // Numerator
  prc_cmp dnm_gcm; // Denominator
  prc_cmp tmp4; // Correction to neutral atmosphere factor
  prc_cmp tmp5; // Neutral atmosphere factor
  prc_cmp mno_prm_rlm; // [frc] Monin-Obukhov parameter at roughness height
  prc_cmp mno_prm_gcm; // [frc] Monin-Obukhov parameter at GCM layer height

  // Main Code
  //Compute stability parameter
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    // Maximum value is 1 because stability correction function is valid only for zeta < 1, e.g., Bon96 p. 52, Bru82 p. 71, SeP97 p. 963
    mno_prm_rlm=min_cpv(rgh_mmn[lon_idx]/mno_lng[lon_idx],1.0); // [frc]
    mno_prm_gcm=min_cpv((hgt_mdp[lon_idx]-hgt_zpd[lon_idx])/mno_lng[lon_idx],1.0); // [frc]
    if(mno_lng[lon_idx] < 0.0){
      // Difference between unstable corrections (pi/2 cancels out)
      eta_sqr_rlm=std::sqrt(1.0-16.0*mno_prm_rlm);
      eta_sqr_gcm=std::sqrt(1.0-16.0*mno_prm_gcm);
      eta_rlm=std::sqrt(eta_sqr_rlm); // [frc] Reciprocal of similarity function for momentum, unstable atmosphere
      eta_gcm=std::sqrt(eta_sqr_gcm); // [frc] Reciprocal of similarity function for momentum, unstable atmosphere
      nmr_rlm=(eta_sqr_rlm+1.0)*(eta_rlm+1.0)*(eta_rlm+1.0);
      dnm_gcm=(eta_sqr_gcm+1.0)*(eta_gcm+1.0)*(eta_gcm+1.0);
      tmp4=std::log(nmr_rlm/dnm_gcm)+2.0*(std::atan(eta_gcm)-std::atan(eta_rlm)); // [frc]
    }else{ // not stable
      // Difference between stable corrections
      tmp4=5.0*(mno_prm_gcm-mno_prm_rlm); // [frc]
    } // not stable
    tmp5=std::log((hgt_mdp[lon_idx]-hgt_zpd[lon_idx])/rgh_mmn[lon_idx]);
    rss_aer[lon_idx]=(tmp4+tmp5)/(cst_von_krm*wnd_frc[lon_idx]); // [s m-1] Bon96 p. 54
  } // end loop over lon

  return rcd;
} // end rss_aer_get()
      
int wnd_rfr_get
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *hgt_mdp, // I [m] Midpoint height above surface
 const prc_cmp hgt_rfr, // I [m] Reference height
 const prc_cmp *hgt_zpd, // I [m] Zero plane displacement
 const prc_cmp *mno_lng, // I [m] Monin-Obukhov length
 const prc_cmp *wnd_frc, // I [m s-1] Surface friction velocity
 const prc_cmp *wnd_mdp, // I [m s-1] Surface layer mean wind speed
 prc_cmp *wnd_rfr) // O [m s-1] Wind speed at reference height
{
  // Purpose: Convert wind speed at given height to wind speed at reference height
  // Taken from Bon96
  // Dependencies: phys_cst.hh
  using phc::cst_von_krm_rcp; // (2.5) [frc] Reciprocal of Von Karman's constant
  using blm::wnd_min_dps; // (1.0) [m s-1] Minimum windspeed used for deposition Bon96 p. 51
  // Local
  a2d_cls<prc_cmp> mno_stb_crc_mmn(lon_nbr,2); // [frc] Monin-Obukhov stability correction term
  a2d_cls<prc_cmp> mno_stb_prm(lon_nbr,2); // [frc] Monin-Obukhov stability parameter 
  prc_cmp sml_fnc_mmn_uns_rcp; // Reciprocal of similarity function for momentum, unstable atmosphere
  prc_cmp tmp2; // Term in stability correction computation
  prc_cmp tmp3; // Term in stability correction computation
  prc_cmp tmp4; // Term in stability correction computation
  std::valarray<prc_cmp> wnd_crc_fct(lon_nbr); // [frc] Wind correction factor
  int rcd(0); // O [rcd] Return success code
  long idx_idx; // Counting index for lon_idx
  long lon_idx; // Counting index for lon
  long lvl_idx; // Stability computation loop index
  
  // Initialize some variables
  const long rfr_hgt_idx(0L); // Named index for lower (target, reference) height
  const long gcm_hgt_idx(1); // Named index for upper (known, GCM) height
  vec_set(wnd_rfr,lon_nbr,0.0); // [m s-1]
  
  std::valarray<long> vld_idx(lon_nbr); // Valid indices
  long vld_nbr; // Number of valid indices
  whenflt(lon_nbr,hgt_zpd,hgt_rfr,&vld_idx[0],vld_nbr);
  
  // Compute horizontal wind speed at reference height
  for(idx_idx=0;idx_idx<vld_nbr;idx_idx++){
    lon_idx=vld_idx[idx_idx];
    // Code uses notation of Bon96 p. 50, where lvl_idx=1 is 10 m ref. hgt, lvl_idx=2 is atm. hgt.
    // Stability functions only valid for stability parameter zeta <= 1.0 Bon96 p. 52, Bru82
    mno_stb_prm(lon_idx,rfr_hgt_idx)=min_cpv((hgt_rfr-hgt_zpd[lon_idx])/mno_lng[lon_idx],1.0); // [frc]
    mno_stb_prm(lon_idx,gcm_hgt_idx)=min_cpv((hgt_mdp[lon_idx]-hgt_zpd[lon_idx])/mno_lng[lon_idx],1.0); // [frc]
    for(lvl_idx=0;lvl_idx<2;lvl_idx++){
      if(mno_stb_prm(lon_idx,lvl_idx) < 0.0){
	sml_fnc_mmn_uns_rcp=std::pow(1.0-16.0*mno_stb_prm(lon_idx,lvl_idx),0.25);
	tmp2=std::log((1.0+sml_fnc_mmn_uns_rcp*sml_fnc_mmn_uns_rcp)/2.0);
	tmp3=std::log((1.0+sml_fnc_mmn_uns_rcp)/2.0);
	mno_stb_crc_mmn(lon_idx,lvl_idx)=2.0*tmp3+tmp2-2.0*std::atan(sml_fnc_mmn_uns_rcp)+1.5707963; // [frc]
      }else{ // not stable
	mno_stb_crc_mmn(lon_idx,lvl_idx)=-5.0*mno_stb_prm(lon_idx,lvl_idx); // [frc]
      } // not stable
    } // end loop over lvl_idx
    tmp4=std::log((hgt_mdp[lon_idx]-hgt_zpd[lon_idx])/(hgt_rfr-hgt_zpd[lon_idx]));
    wnd_crc_fct[lon_idx]=tmp4-mno_stb_crc_mmn(lon_idx,gcm_hgt_idx)+mno_stb_crc_mmn(lon_idx,rfr_hgt_idx); // [frc]
    // Correct neutral stability assumption
    wnd_rfr[lon_idx]=wnd_mdp[lon_idx]-wnd_frc[lon_idx]*cst_von_krm_rcp*wnd_crc_fct[lon_idx]; // [m s-1]
    // Using neutral stability assumption in calm winds does not seem to help
    if(wnd_rfr[lon_idx] < 0.0){
      prc_cmp wnd_rfr_tmp(wnd_rfr[lon_idx]);
      wnd_rfr[lon_idx]=wnd_mdp[lon_idx]-wnd_frc[lon_idx]*cst_von_krm_rcp*tmp4;
      std::cerr << prg_nm_get() << ": WARNING wnd_rfr_get(): Using neutral stability assumption changes U(10 m) from " << wnd_rfr_tmp << " m s-1 to " << wnd_rfr[lon_idx] << " m s-1" << std::endl;
    } // endif 
    // Set reference wind speed to minimum-allowed value
    wnd_rfr[lon_idx]=max_cpv(wnd_rfr[lon_idx],wnd_min_dps); // [m s-1]
  } // end loop over lon
  
  return rcd;
} // end wnd_rfr_get()

int snw_frc_get
(const long lon_nbr, // I [nbr] Size of arrays
 prc_cmp *snw_frc, // O [frc] Fraction of surface covered by snow
 prc_cmp *snw_hgt, // O [m] Geometric bulk thickness of snow
 const prc_cmp *snw_hgt_lqd) // I [m] Equivalent liquid water snow depth
{
  // Purpose: Convert equivalent liquid water snow depth to fractional snow cover
  // Use snow thickness -> fraction algorithm of Bon96
  //  const prc_cmp snw_dns(250.0); // [kg m-3] Bulk density of snow CCM:lsm/snoconi.F (NB: CCM:physics/tsinti() uses 100.0)
  using phc::dns_H2O_snw_gnd_std; // (100.0) [kg m-3] Standard bulk density of snow on ground WiW80 p. 2724, 2725, CCM:physics/tsinti()
  using phc::dns_H2O_lqd_std; // (1000.0) [kg m-3] Density of liquid water
  const prc_cmp snw_hgt_thr(0.05); // [m] Geometric snow thickness for 100% coverage CCM:lsm/snoconi.F 
  // Derived
  const prc_cmp hgt_lqd_snw_cnv(dns_H2O_lqd_std/dns_H2O_snw_gnd_std); // [frc] Conversion factor from liquid water depth to geometric snow thickness

  // Output
  int rcd(0); // O [rcd] Return success code
  // Local
  long lon_idx; // Counting index for lon
  // Main code
  // Fractional snow cover
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    snw_hgt[lon_idx]=snw_hgt_lqd[lon_idx]*hgt_lqd_snw_cnv; // [m] NB: CCM3 and LSM disagree on this
    snw_frc[lon_idx]=min_cpv(snw_hgt[lon_idx]/snw_hgt_thr,1.0); // [frc]
  } // end loop over lon
  return rcd;
} // end snw_frc_get()

int rgh_zpd_get
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *oro,  // I [frc] Orography
 prc_cmp *rgh_mmn, // O [m] Roughness length momentum
 const long *sfc_typ, // I [idx] LSM surface type (0..28)
 const prc_cmp *snw_frc, // I [frc] Fraction of surface covered by snow
 const prc_cmp *wnd_10m, // I [m s-1] 10 m wind speed
 prc_cmp hgt_zpd[]) // O [m] Zero plane displacement
{     
  // Purpose: Set roughness length and zero plane displacement 
  // NB: SeP97 apparently uses rlm but rlh may be more correct for aerosol
  // NB: Currently all lakes are treated as unfrozen
  // Dependencies: <phys_cst.hh>, <lsm.hh>
  using phc::cst_von_krm; // (0.4) [frc] Von Karman's constant
  using lsm::z0mvt; // [m] Momentum roughness length
  using lsm::zpdvt; // [m] Zero plane displacement height
  using lsm::pnt_frc; // [frc] Weight of corresponding plant type (sums to 1.0)
  using lsm::pnt_typ; // [idx] LSM plant type (1..14 = nbr_LSM_pnt_typ)
  using blm::rgh_mmn_ice_lnd; // (0.05) [m] Roughness length over ice, bare ground, wetlands Bon96 p. 59
  using blm::rgh_mmn_ice_ocn; // (0.0005) [m] Roughness length over sea ice BKL97 p. F-3 (updated)
  using blm::rgh_mmn_lak_wrm; // (0.001) [m] Roughness length over unfrozen lakes Bon96 p. 59
  using blm::rgh_mmn_snw; // (0.04) [m] Roughness length over snow Bon96 p. 59 CCM:lsm/snoconi.F
  using blm::wnd_min_dps; // (1.0) [m s-1] Minimum windspeed used for deposition Bon96 p. 51
  // Output
  int rcd(0); // O [rcd] Return success code
  // Local
  prc_cmp rlm_crr(0.0); // [m] Roughness length of current sub-gridscale cell
  prc_cmp wnd_10m_dps; // [m s-1] Bounded wind speed at 10 m
  prc_cmp xch_cff_mmn_ocn_ntr; // [frc] Neutral 10 m drag coefficient over ocean
  std::valarray<long> ice_idx(lon_nbr); // Longitude index array (sea ice)
  long ice_nbr(0L); // Number of sea ice points
  long idx_idx; // Counting index
  std::valarray<long> lnd_idx(lon_nbr); // Longitude index array (land)
  long lnd_nbr(0L); // Number of land points
  long lon_idx; // Counting index
  std::valarray<long> ocn_idx(lon_nbr); // Longitude index array (ocean)
  long ocn_nbr(0L); // Number of ocean points
  long pnt_typ_idx; // [idx] Plant type index NB: pnt_typ
  long sfc_typ_idx; // [idx] Surface type index
  long sgs_idx; // Surface sub-gridscale index
  
  // Main Code
  std::string sbr_nm("rgh_zpd_get");
  
  // Initialize arrays
  vec_set(rgh_mmn,lon_nbr,0.0); // [m]
  vec_set(hgt_zpd,lon_nbr,0.0); // [m]
  
  // Construct ocean, sea-ice, and land vectors
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    if(oro_is_ocn(oro[lon_idx])){
      ocn_idx[ocn_nbr]=lon_idx;
      ocn_nbr++;
    }  // endif
  }  // end loop over lon
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    if(oro_is_ice(oro[lon_idx])){
      ice_idx[ice_nbr]=lon_idx;
      ice_nbr++;
    }  // endif
  }  // end loop over lon
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    if(oro_is_lnd(oro[lon_idx])){
      lnd_idx[lnd_nbr]=lon_idx;
      lnd_nbr++;
    }  // endif
  }  // end loop over lon
  
  // Ocean points
  for(idx_idx=0;idx_idx<ocn_nbr;idx_idx++){
    lon_idx=ocn_idx[idx_idx];
    // Convert wind speed to roughness length over ocean
    wnd_10m_dps=max_cpv(wnd_min_dps,wnd_10m[lon_idx]); // [m s-1]
    xch_cff_mmn_ocn_ntr=xch_cff_mmn_ocn_ntr_get(wnd_10m_dps); // [frc]
    rgh_mmn[lon_idx]=10.0*std::exp(-cst_von_krm/std::sqrt(xch_cff_mmn_ocn_ntr)); // [m] BKL97 p. F-4, LaP81 p. 327 (14) 
    hgt_zpd[lon_idx]=0.0; // [m]
  }  // end loop over ocn lon
  
  // Sea ice points
  for(idx_idx=0;idx_idx<ice_nbr;idx_idx++){
    lon_idx=ice_idx[idx_idx];
    rgh_mmn[lon_idx]=snw_frc[lon_idx]*rgh_mmn_snw+(1.0-snw_frc[lon_idx])*rgh_mmn_ice_ocn; // [m] Bon96 p. 59
    hgt_zpd[lon_idx]=0.0; // [m]
  }  // end loop over ice lon
  
  // Land points
  for(idx_idx=0;idx_idx<lnd_nbr;idx_idx++){
    lon_idx=lnd_idx[idx_idx];
    // Store surface blend for current gridpoint
    sfc_typ_idx=sfc_typ[lon_idx];
    if(sfc_typ_idx == 0){ // Inland lake
      // NB: DBG XXX Need to add temperature input and so ability to discriminate warm from frozen lakes here
      rgh_mmn[lon_idx]=rgh_mmn_lak_wrm; // [m] Bon96 p. 59
      hgt_zpd[lon_idx]=0.0; // [m]
    }else if(sfc_typ_idx == 1){ // Land ice
      rgh_mmn[lon_idx]=snw_frc[lon_idx]*rgh_mmn_snw+(1.0-snw_frc[lon_idx])*rgh_mmn_ice_lnd; // [m] Bon96 p. 59
      hgt_zpd[lon_idx]=0.0; // [m]
    }else{ // Normal land
      for(sgs_idx=0;sgs_idx<3;sgs_idx++){
	// Bare ground is pnt_typ=14, ocean is pnt_typ=0
	pnt_typ_idx=pnt_typ[sgs_idx][sfc_typ_idx]; // NB: C pnt_typ is transposed from Fortran 
	if(pnt_typ_idx == 14){ // Bare ground
	  rlm_crr=snw_frc[lon_idx]*rgh_mmn_snw+(1.0-snw_frc[lon_idx])*rgh_mmn_ice_lnd; // [m] Bon96 p. 59 (glacial ice is same as bare ground)
	}else if(pnt_typ_idx > 0){ // Regular plant type
	  rlm_crr=snw_frc[lon_idx]*rgh_mmn_snw+(1.0-snw_frc[lon_idx])*z0mvt[pnt_typ_idx]; // [m] Bon96 p. 59
	}else{ // Presumably ocean snuck through
	  err_prn(sbr_nm,"pnt_typ_idx == 0");
	} // endif
	rgh_mmn[lon_idx]+=pnt_frc[sgs_idx][sfc_typ_idx]*rlm_crr; // [m] NB: C pnt_frc is transposed from Fortran
	hgt_zpd[lon_idx]+=pnt_frc[sgs_idx][sfc_typ_idx]*zpdvt[pnt_typ_idx]; // [m] NB: C pnt_frc is transposed from Fortran
      }// end loop over number of possible plant types
    } // endif normal land
  } // end loop over lon
  return rcd;
} // end rgh_zpd_get()

int rgh_mmn_get
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp oro[],  // I [frc] Orography
 prc_cmp rgh_mmn[], // O [m] Roughness length momentum
 const long sfc_typ[], // I [idx] LSM surface type (0..28)
 const prc_cmp snw_frc[], // I [frc] Fraction of surface covered by snow
 const prc_cmp wnd_10m[]) // I [m s-1] 10 m wind speed
{     
  // Purpose: Set roughness length
  // NB: SeP97 apparently uses rlm but rlh may be more correct for aerosol
  // NB: Currently all lakes are treated as unfrozen
  // Dependencies: <phys_cst.hh>, <lsm.hh>
  using phc::cst_von_krm; // (0.4) [frc] Von Karman's constant
  using lsm::pnt_frc; // [frc] Weight of corresponding plant type (sums to 1.0)
  using lsm::pnt_typ; // [idx] LSM plant type (1..14 = nbr_LSM_pnt_typ)
  using lsm::z0mvt; // [m] Momentum roughness length
  using blm::rgh_mmn_ice_lnd; // (0.05) [m] Roughness length over ice, bare ground, wetlands Bon96 p. 59
  using blm::rgh_mmn_ice_ocn; // (0.0005) [m] Roughness length over sea ice BKL97 p. F-3 CCM:dom/parpbl.h
  using blm::rgh_mmn_lak_wrm; // (0.001) [m] Roughness length over unfrozen lakes Bon96 p. 59
  using blm::rgh_mmn_snw; // (0.04) [m] Roughness length over snow Bon96 p. 59 CCM:lsm/snoconi.F
  using blm::wnd_min_dps; // (1.0) [m s-1] Minimum windspeed used for deposition Bon96 p. 51
  
  // Output
  int rcd(0); // O [rcd] Return success code
  // Local
  std::valarray<long> ice_idx(lon_nbr); // Longitude index array (sea ice)
  long ice_nbr(0L); // Number of sea ice points
  long idx_idx; // Counting index
  std::valarray<long> lnd_idx(lon_nbr); // Longitude index array (land)
  long lnd_nbr(0L); // Number of land points
  long lon_idx; // Counting index
  std::valarray<long> ocn_idx(lon_nbr); // Longitude index array (ocean)
  long ocn_nbr(0L); // Number of ocean points
  long pnt_typ_idx; // [idx] Plant type index NB: pnt_typ
  long sfc_typ_idx; // [idx] Surface type index
  long sgs_idx; // Surface sub-gridscale index
  prc_cmp rlm_crr(0.0); // [m] Roughness length of current sub-gridscale cell
  prc_cmp wnd_10m_dps; // [m s-1] Bounded wind speed at 10 m
  prc_cmp xch_cff_mmn_ocn_ntr; // [frc] Neutral 10 m drag coefficient over ocean
  
  // Main Code
  std::string sbr_nm("rgh_mmn_get");
  
  // Initialize arrays
  vec_set(rgh_mmn,lon_nbr,0.0); // [m]
  
  // Construct ocean, sea-ice, and land vectors
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    if(oro_is_ocn(oro[lon_idx])){
      ocn_idx[ocn_nbr]=lon_idx;
      ocn_nbr++;
    }  // endif
  }  // end loop over lon
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    if(oro_is_ice(oro[lon_idx])){
      ice_idx[ice_nbr]=lon_idx;
      ice_nbr++;
    }  // endif
  }  // end loop over lon
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    if(oro_is_lnd(oro[lon_idx])){
      lnd_idx[lnd_nbr]=lon_idx;
      lnd_nbr++;
    }  // endif
  }  // end loop over lon
  
  // Ocean points
  for(idx_idx=0;idx_idx<ocn_nbr;idx_idx++){
    lon_idx=ocn_idx[idx_idx];
    // Convert wind speed to roughness length over ocean
    wnd_10m_dps=max_cpv(wnd_min_dps,wnd_10m[lon_idx]); // [m s-1]
    xch_cff_mmn_ocn_ntr=xch_cff_mmn_ocn_ntr_get(wnd_10m_dps); // [frc]
    rgh_mmn[lon_idx]=10.0*std::exp(-cst_von_krm/std::sqrt(xch_cff_mmn_ocn_ntr)); // [m] BKL97 p. F-4, LaP81 p. 327 (14) 
  }  // end loop over ocn lon
  
  // Sea ice points
  for(idx_idx=0;idx_idx<ice_nbr;idx_idx++){
    lon_idx=ice_idx[idx_idx];
    rgh_mmn[lon_idx]=snw_frc[lon_idx]*rgh_mmn_snw+(1.0-snw_frc[lon_idx])*rgh_mmn_ice_ocn; // [m] Bon96 p. 59
  }  // end loop over ice lon
  
  // Land points
  for(idx_idx=0;idx_idx<lnd_nbr;idx_idx++){
    lon_idx=lnd_idx[idx_idx];
    // Store surface blend for current gridpoint
    sfc_typ_idx=sfc_typ[lon_idx];
    if(sfc_typ_idx == 0){ // Inland lake
      // NB: DBG XXX Need to add temperature input and so ability to discriminate warm from frozen lakes here
      rgh_mmn[lon_idx]=rgh_mmn_lak_wrm; // [m] Bon96 p. 59
    }else if(sfc_typ_idx == 1){ // Land ice
      rgh_mmn[lon_idx]=snw_frc[lon_idx]*rgh_mmn_snw+(1.0-snw_frc[lon_idx])*rgh_mmn_ice_lnd; // [m] Bon96 p. 59
    }else{ // Normal land
      for(sgs_idx=0;sgs_idx<3;sgs_idx++){
	// Bare ground is pnt_typ=14, ocean is pnt_typ=0
	pnt_typ_idx=pnt_typ[sgs_idx][sfc_typ_idx]; // NB: C pnt_typ is transposed from Fortran 
	if(pnt_typ_idx == 14){ // Bare ground
	  rlm_crr=snw_frc[lon_idx]*rgh_mmn_snw+(1.0-snw_frc[lon_idx])*rgh_mmn_ice_lnd; // [m] Bon96 p. 59 (glacial ice is same as bare ground)
	}else if(pnt_typ_idx > 0){ // Regular plant type
	  rlm_crr=snw_frc[lon_idx]*rgh_mmn_snw+(1.0-snw_frc[lon_idx])*z0mvt[pnt_typ_idx]; // [m] Bon96 p. 59
	}else{ // Presumably ocean snuck through
	  err_prn(sbr_nm,"dst: ERROR: pnt_typ_idx == 0 in rgh_mmn_get())");
	} // endif
	rgh_mmn[lon_idx]+=pnt_frc[sgs_idx][sfc_typ_idx]*rlm_crr; // [m] NB: C pnt_frc is transposed from Fortran
      }// end loop over number of possible plant types
    } // endif normal land
  } // end loop over lon
  return rcd;
} // end rgh_mmn_get()

int rgh_get
(const long lon_nbr, // I [nbr] Size of arrays
 const bool *vgt, // I [flg] "Vegetated" flag
 const bool *lak, // I [flg] Lake flag
 const prc_cmp *tpt_gnd, // I [K] Ground temperature
 const prc_cmp *oro,  // I [frc] Orography
 prc_cmp *rgh_mmn_gnd, // O [m] Roughness length momentum ground
 prc_cmp *rgh_mmn, // O [m] Roughness length momentum
 prc_cmp *rgh_heat, // O [m] Roughness length heat
 prc_cmp *rgh_heat_gnd, // O [m] Roughness length heat ground
 const long *sfc_typ, // I [idx] LSM surface type (0..28)
 const prc_cmp *snw_frc, // I [frc] Fraction of surface covered by snow
 const prc_cmp *wnd_10m) // I [m s-1] 10 m wind speed
{     
  // Purpose: Set roughness length
  // NB: SeP97 apparently uses rlm but rlh may be more correct for aerosol
  // Dependencies: <phys_cst.hh>, <lsm.hh>
  using phc::tpt_frz_pnt; // (273.15) [K] Kelvin--Celsius scale offset, see, e.g., Bol80
  using phc::cst_von_krm; // (0.4) [frc] Von Karman's constant
  using lsm::pnt_frc; // [frc] Weight of corresponding plant type (sums to 1.0)
  using lsm::pnt_typ; // [idx] LSM plant type (1..14 = nbr_LSM_pnt_typ)
  using lsm::z0mvt; // [m] Momentum roughness length
  using blm::rgh_mmn_ice_lak; // (0.04) [m] Roughness length over frozen lakes Bon96 p. 59
  using blm::rgh_mmn_ice_lnd; // (0.05) [m] Roughness length over ice, bare ground, wetlands Bon96 p. 59
  using blm::rgh_mmn_ice_ocn; // (0.0005) [m] Roughness length over sea ice BKL97 p. F-3 CCM:dom/parpbl.h
  using blm::rgh_mmn_lak_wrm; // (0.001) [m] Roughness length over unfrozen lakes Bon96 p. 59
  using blm::rgh_mmn_snw; // (0.04) [m] Roughness length over snow Bon96 p. 59 CCM:lsm/snoconi.F
  using blm::wnd_min_dps; // (1.0) [m s-1] Minimum windspeed used for deposition Bon96 p. 51
  
  // Output
  int rcd(0); // O [rcd] Return success code
  // Local
  std::valarray<long> ice_idx(lon_nbr); // Longitude index array (sea ice)
  long ice_nbr(0L); // Number of sea ice points
  long idx_idx; // Counting index
  std::valarray<long> lnd_idx(lon_nbr); // Longitude index array (land)
  long lnd_nbr(0L); // Number of land points
  long lon_idx; // Counting index
  std::valarray<long> ocn_idx(lon_nbr); // Longitude index array (ocean)
  long ocn_nbr(0L); // Number of ocean points
  long sfc_typ_idx; // [idx] Surface type index
  long sgs_idx; // Surface sub-gridscale index
  prc_cmp rlm_crr(0.0); // [m] Roughness length of current sub-gridscale cell
  prc_cmp wnd_10m_dps; // [m s-1] Bounded wind speed at 10 m
  prc_cmp xch_cff_mmn_ocn_ntr; // [frc] Neutral 10 m drag coefficient over ocean
  long pnt_typ_idx; // [idx] Plant type index 
  
  // Main Code
  std::string sbr_nm("rgh_get");
  
  // Initialize arrays
  vec_set(rgh_mmn,lon_nbr,0.0); // [m]
  
  // Construct ocean, sea-ice, and land vectors
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    if(oro_is_ocn(oro[lon_idx])){
      ocn_idx[ocn_nbr]=lon_idx;
      ocn_nbr++;
    }  // endif
  }  // end loop over lon
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    if(oro_is_ice(oro[lon_idx])){
      ice_idx[ice_nbr]=lon_idx;
      ice_nbr++;
    }  // endif
  }  // end loop over lon
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    if(oro_is_lnd(oro[lon_idx])){
      lnd_idx[lnd_nbr]=lon_idx;
      lnd_nbr++;
    }  // endif
  }  // end loop over lon
  
  // Ocean points
  for(idx_idx=0;idx_idx<ocn_nbr;idx_idx++){
    lon_idx=ocn_idx[idx_idx];
    // Convert wind speed to roughness length over ocean
    wnd_10m_dps=max_cpv(wnd_min_dps,wnd_10m[lon_idx]); // [m s-1]
    xch_cff_mmn_ocn_ntr=xch_cff_mmn_ocn_ntr_get(wnd_10m_dps); // [frc]
    rgh_mmn[lon_idx]=10.0*std::exp(-cst_von_krm/std::sqrt(xch_cff_mmn_ocn_ntr)); // [m] BKL97 p. F-4, LaP81 p. 327 (14) 
  }  // end loop over ocn lon
  
  // Sea ice points
  for(idx_idx=0;idx_idx<ice_nbr;idx_idx++){
    lon_idx=ice_idx[idx_idx];
    rgh_mmn[lon_idx]=snw_frc[lon_idx]*rgh_mmn_snw+(1.0-snw_frc[lon_idx])*rgh_mmn_ice_ocn; // [m] Bon96 p. 59
  }  // end loop over ice lon
  
  // Land points
  for(idx_idx=0;idx_idx<lnd_nbr;idx_idx++){
    lon_idx=lnd_idx[idx_idx];
    // Store surface blend for current gridpoint
    sfc_typ_idx=sfc_typ[lon_idx];
    if(lak[lon_idx] && tpt_gnd[lon_idx] < tpt_frz_pnt){
      rgh_mmn_gnd[lon_idx]=snw_frc[lon_idx]*rgh_mmn_snw+(1.0-snw_frc[lon_idx])*rgh_mmn_ice_lak; // [m] Bon96 p. 59
    }else{
      //      rgh_mmn_gnd[lon_idx]=snw_frc[lon_idx]*rgh_mmn_snw+(1.0-snw_frc[lon_idx])*z0msoi[soi_typ_idx]; // [m] Bon96 p. 59; // [m] Bon96 p. 59
    } // endif
    if(sfc_typ_idx == 0){ // Inland lake
      rgh_mmn[lon_idx]=rgh_mmn_lak_wrm; // [m] Bon96 p. 59
    }else if(sfc_typ_idx == 1){ // Land ice
      rgh_mmn[lon_idx]=snw_frc[lon_idx]*rgh_mmn_snw+(1.0-snw_frc[lon_idx])*rgh_mmn_ice_lnd; // [m] Bon96 p. 59
    }else{ // Normal land
      for(sgs_idx=0;sgs_idx<3;sgs_idx++){
	// Bare ground is pnt_typ=14, ocean is pnt_typ=0
	pnt_typ_idx=pnt_typ[sgs_idx][sfc_typ_idx]; // NB: C pnt_typ is transposed from Fortran 
	if(pnt_typ_idx == 14){ // Bare ground
	  rlm_crr=snw_frc[lon_idx]*rgh_mmn_snw+(1.0-snw_frc[lon_idx])*rgh_mmn_ice_lnd; // [m] Bon96 p. 59 (glacial ice is same as bare ground)
	}else if(pnt_typ_idx > 0){ // Regular plant type
	  rlm_crr=snw_frc[lon_idx]*rgh_mmn_snw+(1.0-snw_frc[lon_idx])*z0mvt[pnt_typ_idx]; // [m] Bon96 p. 59
	}else{ // Presumably ocean snuck through
	  err_prn(sbr_nm,"dst: ERROR: pnt_typ_idx == 0 in rgh_get())");
	} // endif
	rgh_mmn[lon_idx]+=pnt_frc[sgs_idx][sfc_typ_idx]*rlm_crr; // [m] NB: C pnt_frc is transposed from Fortran
      }// end loop over number of possible plant types
    } // endif normal land
  } // end loop over lon

  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    rgh_heat[lon_idx]=rgh_mmn[lon_idx]; // O [m] Roughness length heat
    rgh_heat_gnd[lon_idx]=rgh_mmn[lon_idx]; // O [m] Roughness length heat ground
  } // end loop over lon

  if(vgt[0]){;} // CEWU Compiler Error Warning Usage

  return rcd;
} // end rgh_get()

int zpd_get
(const long lon_nbr, // I [nbr] Size of arrays
 prc_cmp *hgt_zpd, // O [m] Zero plane displacement
 const prc_cmp *oro,  // I [frc] Orography
 const long *sfc_typ) // I [idx] LSM surface type (0..28)
{     
  // Purpose: Set zero plane displacement 
  // Dependencies: <lsm.hh>
  using lsm::pnt_frc; // [frc] Weight of corresponding plant type (sums to 1.0)
  using lsm::pnt_typ; // [idx] LSM plant type (1..14 = nbr_LSM_pnt_typ)
  using lsm::zpdvt; // [m] Zero plane displacement height
  // Output
  int rcd(0); // O [rcd] Return success code
  // Local
  long idx_idx; // Counting index
  std::valarray<long> lnd_idx(lon_nbr); // Longitude index array (land)
  long lnd_nbr(0L); // Number of land points
  long lon_idx; // Counting index
  long pnt_typ_idx; // [idx] Plant type index NB: pnt_typ
  long sfc_typ_idx; // [idx] Surface type index
  long sgs_idx; // Surface sub-gridscale index
  
  // Main Code
  std::string sbr_nm("zpd_get");
  
  // Initialize arrays
  // Zero plane displacement is identically 0.0 everywhere except land
  vec_set(hgt_zpd,lon_nbr,0.0); // [m]
  
  // Land ahoy!
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    if(oro_is_lnd(oro[lon_idx])){
      lnd_idx[lnd_nbr]=lon_idx;
      lnd_nbr++;
    }  // endif
  }  // end loop over lon
  
  // Land points
  for(idx_idx=0;idx_idx<lnd_nbr;idx_idx++){
    lon_idx=lnd_idx[idx_idx];
    // Store surface blend for current gridpoint
    sfc_typ_idx=sfc_typ[lon_idx];
    if(sfc_typ_idx == 0){ // Inland lake
      hgt_zpd[lon_idx]=0.0; // [m]
    }else if(sfc_typ_idx == 1){ // Land ice
      hgt_zpd[lon_idx]=0.0; // [m]
    }else{ // Normal land
      for(sgs_idx=0;sgs_idx<3;sgs_idx++){
	// Bare ground is pnt_typ=14, ocean is pnt_typ=0
	pnt_typ_idx=pnt_typ[sgs_idx][sfc_typ_idx]; // NB: C pnt_typ is transposed from Fortran 
	hgt_zpd[lon_idx]+=pnt_frc[sgs_idx][sfc_typ_idx]*zpdvt[pnt_typ_idx]; // [m] NB: C pnt_frc is transposed from Fortran
      }// end loop over number of possible plant types
    } // endif normal land
  } // end loop over lon
  return rcd;
} // end zpd_get()

int sfc_ems_get
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *oro,  // I [frc] Orography
 const long *sfc_typ, // I [idx] LSM surface type (0..28)
 prc_cmp *msv_sfc) // O [m] Zero plane displacement
{     
  // Purpose: Set zero plane displacement 
  // Dependencies: <lsm.hh>
  using lsm::zpdvt; // [m] Zero plane displacement height
  using lsm::pnt_typ; // [idx] LSM plant type (1..14 = nbr_LSM_pnt_typ)
  using lsm::pnt_frc; // [frc] Weight of corresponding plant type (sums to 1.0)
  // Output
  int rcd(0); // O [rcd] Return success code
  // Local
  long idx_idx; // Counting index
  std::valarray<long> lnd_idx(lon_nbr); // Longitude index array (land)
  long lnd_nbr(0L); // Number of land points
  long lon_idx; // Counting index
  long pnt_typ_idx; // [idx] Plant type index NB: pnt_typ
  long sfc_typ_idx; // [idx] Surface type index
  long sgs_idx; // Surface sub-gridscale index
  
  // Main Code
  std::string sbr_nm("sfc_ems_get");
  
  // Initialize arrays
  // Zero plane displacement is identically 0.0 everywhere except land
  vec_set(msv_sfc,lon_nbr,0.0); // [m]
  
  // Land ahoy!
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    if(oro_is_lnd(oro[lon_idx])){
      lnd_idx[lnd_nbr]=lon_idx;
      lnd_nbr++;
    }  // endif
  }  // end loop over lon
  
  // Land points
  for(idx_idx=0;idx_idx<lnd_nbr;idx_idx++){
    lon_idx=lnd_idx[idx_idx];
    // Store surface blend for current gridpoint
    sfc_typ_idx=sfc_typ[lon_idx];
    if(sfc_typ_idx == 0){ // Inland lake
      msv_sfc[lon_idx]=0.0; // [m]
    }else if(sfc_typ_idx == 1){ // Land ice
      msv_sfc[lon_idx]=0.0; // [m]
    }else{ // Normal land
      for(sgs_idx=0;sgs_idx<3;sgs_idx++){
	// Bare ground is pnt_typ=14, ocean is pnt_typ=0
	pnt_typ_idx=pnt_typ[sgs_idx][sfc_typ_idx]; // NB: C pnt_typ is transposed from Fortran 
	msv_sfc[lon_idx]+=pnt_frc[sgs_idx][sfc_typ_idx]*zpdvt[pnt_typ_idx]; // [m] NB: C pnt_frc is transposed from Fortran
      }// end loop over number of possible plant types
    } // endif normal land
  } // end loop over lon
  return rcd;
} // end sfc_ems_get()

int 
trn_fsh_vpr_soi_atm_get
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *tpt_soi, // I [K] Soil temperature
 const prc_cmp tpt_soi_frz, // I [K] Temperature of frozen soil
 prc_cmp *trn_fsh_vpr_soi_atm, // O [frc] Transfer efficiency of vapor from soil to atmosphere
 const prc_cmp *vwc_dry, // I [m3 m-3] Dry volumetric water content (no E-T)
 const prc_cmp *vwc_opt, // I [m3 m-3] E-T optimal volumetric water content
 const prc_cmp *vwc_sfc) // I [m3 m-3] Volumetric water content
{
  // Purpose: Compute factor describing effects of soil texture and moisture on vapor transfer between soil and atmosphere
  /* The trn_fsh_vpr_soi_atm efficiency factor attempts to tie soil texture and 
     moisture properties to the vapor conductance of the soil-atmosphere system.
     When the soil temperature is sub-freezing, the conductance describes the 
     resistance to vapor sublimation (or deposition) and transport through the 
     open soil pores to the atmosphere.
     For warm soils, vapor transfer is most efficient at the optimal VWC for E-T
     Thus when vwc_sfc = vwc_opt, soil vapor transfer is perfectly efficient 
     (trn_fsh_vpr_soi_atm = 1.0) so the soil does not contribute any resistance  
     to the surface vapor transfer.
     When vwc_sfc > vwc_opt, the soil has an excess of moisture and, again,
     vapor transfer is not limited by soil characteristics.
     In fact, according to Bon96 p. 98, vwc_dry is only slightly smaller than
     vwc_opt, so trn_fsh_vpr_soi_atm is usually either 0 or 1 and intermediate
     efficiencies occur over only a relatively small range of VWC.
     When vwc_sfc < vwc_dry, the soil matrix is subsaturated and acts as a 
     one-way sink for vapor through osmotic and capillary potentials.
     In this case trn_fsh_vpr_soi_atm = 0, which would cause the surface resistance
     rss_vpr_sfc to blow up, but this is guarded against and rss_sfc_vpr 
     is set to ~1.0e6*rss_aer_vpr instead. 
     Note that this formulation does not seem to allow vapor transfer from
     the atmosphere to the soil when vwc_sfc < vwc_dry, even when 
     e_atm > esat(Tg).
     Air at the apparent sink for moisture is has vapor pressure es = e_sfc
     e_atm = Vapor pressure of ambient air at z = hgt_mdp
     e_sfc = Vapor pressure at apparent sink for moisture at z = zpd + rgh_vpr 
     e_gnd = Vapor pressure at air/ground interface temperature 
     Air at the soil interface is assumed saturated, i.e., e_gnd = esat(Tg)
  */
  // Requires: <phys_cst.hh>, <utl.hh>
  // Taken from Bon96 p. 59, CCM:lsm/surphys
  // Output
  int rcd(0); // O [rcd] Return success code
  // Local
  long lon_idx; // Counting index
  // Main code
  const prc_cmp trn_fsh_vpr_soi_atm_frz(0.01); // [frc] Soil vapor transfer efficiency of frozen soil
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    if(tpt_soi[lon_idx] > tpt_soi_frz){
      trn_fsh_vpr_soi_atm[lon_idx]= // [frc] Soil vapor transfer efficiency
	min_cpv(max_cpv(vwc_sfc[lon_idx]-vwc_dry[lon_idx],0.0)/(vwc_opt[lon_idx]-vwc_dry[lon_idx]),1.0); // CCM: lsm/surphys Bon96 p. 59
    }else{
      trn_fsh_vpr_soi_atm[lon_idx]=trn_fsh_vpr_soi_atm_frz; // [frc] Bon96 p. 59
    } // endelse      
  }  // end loop over lon
  return rcd;
} // end trn_fsh_vpr_soi_atm_get()


