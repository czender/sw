// $Id$ 

// Purpose: Mineral dust physics utilities for C++ programs

/* Copyright (C) 1997--present Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

#include <mnr_dst.hh> // Mineral dust physics

int flx_mss_hrz_slt_ttl_Whi79_get
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *dns_mdp, // I [kg m-3] Midlayer density
 prc_cmp *flx_mss_hrz_slt_ttl, // O [kg m-1 s-1] Vertically integrated streamwise mass flux
 const prc_cmp *wnd_frc, // I [m s-1] Surface friction velocity
 const prc_cmp *wnd_frc_thr_slt) // I [m s-1] Threshold friction speed for saltation
{     
  // Purpose: Compute vertically integrated streamwise mass flux of particles
  // Theory: Uses method proposed by White (1979)
  /* Validation: 
     Whi79 p. 4650 Fig. 9 shows the streamwise mass flux for Earth and Mars
     Simulating a strong wind on Earth with, e.g., mie --wnd_znl_mdp=20.0,
     yields good agreement with (x,y) = (0.32,0.844) g cm-1 s-1 */
  // Dependencies: <phys_cst.hh>
  using phc::grv_sfc_mean; // (9.80665) [m s-2] Mean gravitational acceleration at Earth's surface
  // Output
  int rcd(0); // O [rcd] Return success code
  // Local
  const prc_cmp cst_slt(2.61); // [frc] Whi79 p. 4648, MaB97 p. 16422 
  prc_cmp wnd_frc_rat; // [frc] Ratio of wind friction threshold to wind friction
  long idx_idx; // Counting index for lon_idx
  long lon_idx; // Counting index for lon
  std::valarray<long> vld_idx(lon_nbr); // Valid indices
  long vld_nbr(0L); // Number of valid indices
  
  // Initialize some variables
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    if(wnd_frc[lon_idx] > wnd_frc_thr_slt[lon_idx]){
      vld_idx[vld_nbr++]=lon_idx;
    }else{
      flx_mss_hrz_slt_ttl[lon_idx]=0.0; // [kg m-1 s-1]
    } // endif
  }  // end loop over lon
  
  for(idx_idx=0;idx_idx<vld_nbr;idx_idx++){
    lon_idx=vld_idx[idx_idx];
    wnd_frc_rat=wnd_frc_thr_slt[lon_idx]/wnd_frc[lon_idx]; // [frc]
    flx_mss_hrz_slt_ttl[lon_idx]= // [kg m-1 s-1] 
      cst_slt*dns_mdp[lon_idx]*std::pow(wnd_frc[lon_idx],PRC_CMP(3.0))*(1.0-wnd_frc_rat)*(1.0+wnd_frc_rat)*(1.0+wnd_frc_rat)/grv_sfc_mean; // Whi79 p. 4648 (19), MaB97 p. 16422 (28)
  }  // end loop over lon
  
  return rcd;
} // end flx_mss_hrz_slt_ttl_Whi79_get()

int flx_mss_vrt_dst_ttl_MaB95_get
(const long lon_nbr, // I [nbr] Size of arrays
 prc_cmp *dst_slt_flx_rat_ttl, // O [m-1] Ratio of vertical dust flux to streamwise mass flux
 const prc_cmp *flx_mss_hrz_slt_ttl, // I [kg m-1 s-1] Vertically integrated streamwise mass flux
 prc_cmp *flx_mss_vrt_dst_ttl, // O [kg m-2 s-1] Total vertical mass flux of dust
 const prc_cmp *mss_frc_cly) // I [frc] Mass fraction clay 
{     
  // Purpose: Diagnose total vertical mass flux of dust from vertically integrated streamwise mass flux
  /* Theory: Uses clay-based method proposed by Marticorena & Bergametti (1997)
     Their parameterization is based only on data for mss_frc_cly < 0.20
     For clayier soils, dst_slt_flx_rat_ttl may behave dramatically different
     How this behavior changes when mss_frc_cly > 0.20 is unknown, however
     Thus we use min[mss_frc_cly,0.20] in their parameterization */
  // Dependencies: <phys_cst.hh>
  // Output
  int rcd(0); // O [rcd] Return success code
  // Local
  long lon_idx; // Counting index for lon
  prc_cmp mss_frc_cly_vld; // [frc] Mass fraction clay limited to 0.20
  
  // Initialize some variables
  const prc_cmp ln10(std::log(10.0)); // Natural log of 10
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    mss_frc_cly_vld=min_cpv(mss_frc_cly[lon_idx],0.20); // [frc]
    dst_slt_flx_rat_ttl[lon_idx]= // [m-1]
      100.0*std::exp(ln10*(13.4*mss_frc_cly_vld-6.0)); // MaB95 p. 16423 (47)
    flx_mss_vrt_dst_ttl[lon_idx]=flx_mss_hrz_slt_ttl[lon_idx]*dst_slt_flx_rat_ttl[lon_idx]; // [kg m-1 s-1] 
  }  // end loop over lon
  
  return rcd;
} // end flx_mss_vrt_dst_ttl_MaB95_get()

int // O [rcd] Return success code
wnd_frc_thr_get
(const long sz_nbr, // I [nbr] Size of arrays
 const prc_cmp *dmt_ctr, // I [m] Diameter at bin center
 const prc_cmp dns_prt, // I [kg m-3] Density of particle
 const prc_cmp dns_mdp, // I [kg m-3] Midlayer density
 const prc_cmp vsc_knm_atm, // I [m2 s-1] Kinematic viscosity of atmosphere 
 prc_cmp *ryn_nbr_frc_thr, // O [frc] Threshold friction Reynolds number
 prc_cmp *ryn_nbr_frc_thr_prx, // O [frc] Threshold friction Reynolds number approximation
 prc_cmp *wnd_frc_thr, // O [m s-1] Threshold friction speed
 prc_cmp *wnd_frc_thr_prx) // O [m s-1] Threshold friction speed approximation
{
  /* Purpose: */
  using phc::grv_sfc_mean; // (9.80665) [m s-2] Mean gravitational acceleration at Earth's surface
  // Output
  int rcd(0); // O [rcd] Return success code
  // Local
  const prc_cmp eps_max(1.0e-4); // [frc] Relative accuracy for convergence
  prc_cmp K_MaB95_4; // Factor in threshold speed calculation
  prc_cmp eps_crr; // [frc] Current relative accuracy
  prc_cmp wnd_frc_thr_old; // [m s-1] Previous threshold friction speed
  long itr_idx; // [idx] Counting index
  long idx; // Counting index

  // Main code
  const std::string sbr_nm("wnd_frc_thr_get"); // Subroutine name
  const std::string prg_nm(prg_nm_get()); // Program name
  const unsigned short dbg_lvl(dbg_lvl_get()); // Debugging level

  // Iterative solution for threshold friction speed
  for(idx=0;idx<sz_nbr;idx++){
    // Initialize accuracy and counter
    eps_crr=eps_max+1.0; // [frc] Current relative accuracy
    itr_idx=1; // Counting index
    // Initial guess for threshold friction speed
    if(dmt_ctr[idx] < 10.0e-6) wnd_frc_thr[idx]=1.0; else wnd_frc_thr[idx]=0.2; // [m s-1] 
    // Properties of K do not require iteration
    K_MaB95_4=std::sqrt(1.0+6.0e-07/(dns_prt*grv_sfc_mean*std::pow(dmt_ctr[idx],PRC_CMP(2.5)))); // IvW82 p. 115 (6) MaB95 p. 16417 (4)
    K_MaB95_4*=std::sqrt(dns_prt*grv_sfc_mean*dmt_ctr[idx]/dns_mdp); // IvW82 p. 115 (6) MaB95 p. 16417 (4)
    if(dbg_lvl == dbg_old){
      (void)std::fprintf(stderr,"%4s %4s %7s %7s %8s %8s %8s %8s %8s\n",
		    "idx ","itr ","dmt_ctr","   K   "," Rprx  "," Re*t  ","  uprx ","  u*t  ","  eps  ");
      (void)std::fprintf(stderr,"%4s %4s %7s %7s %8s %8s %8s %8s %8s\n",
		    "    ","    ","  um   ","       ","  frc  ","  frc  "," m s-1 "," m s-1 ","  frc  ");
    } // end if dbg
    ryn_nbr_frc_thr_prx[idx]=0.38+1331.0*std::pow(PRC_CMP(100.0)*dmt_ctr[idx],PRC_CMP(1.56)); // [frc] "B" MaB95 p. 16417 (5)
    if(ryn_nbr_frc_thr_prx[idx] < 10.0){
      if(ryn_nbr_frc_thr_prx[idx] < 0.03) wrn_prn(prg_nm,sbr_nm,"ryn_nbr_frc_thr_prx[idx] < 0.03");
      wnd_frc_thr_prx[idx]=0.1291*K_MaB95_4/std::sqrt(-1.0+1.928*std::pow(ryn_nbr_frc_thr_prx[idx],PRC_CMP(0.0922))); // [m s-1] IvW82 p. 114 (3), MaB95 p. 16417 (6)
    }else{ // ryn_nbr_frc_thr_prx[idx] > 10.0
      wnd_frc_thr_prx[idx]=0.120*K_MaB95_4*(1.0-0.0858*std::exp(-0.0617*(ryn_nbr_frc_thr_prx[idx]-10.0))); // [m s-1] IvW82 p. 115 (6), MaB95 p. 16417 (7)
    } // end if
    while(eps_crr > eps_max){
      // Save threshold speed for convergence test
      wnd_frc_thr_old=wnd_frc_thr[idx]; // [m s-1] 
      // Set threshold friction Reynolds number
      ryn_nbr_frc_thr[idx]=wnd_frc_thr[idx]*dmt_ctr[idx]/vsc_knm_atm; // [frc] "B" IvW82 p. 111
      // Update threshold speed based on new friction Reynolds number
      if(ryn_nbr_frc_thr[idx] < 10.0){
	wnd_frc_thr[idx]=0.1291*K_MaB95_4/std::sqrt(1.928*std::pow(ryn_nbr_frc_thr[idx],PRC_CMP(0.0922))-1.0); // [m s-1] IvW82 p. 114 (3), MaB95 p. 16417 (3)
      }else{ // ryn_nbr_frc_thr[idx] > 10.0
	wnd_frc_thr[idx]=0.120*K_MaB95_4*(1.0-0.0858*std::exp(-0.0617*(ryn_nbr_frc_thr[idx]-10.0))); // [m s-1] IvW82 p. 115 (6), MaB95 p. 16417 (4)
      } // end if
      eps_crr=std::fabs((wnd_frc_thr[idx]-wnd_frc_thr_old)/wnd_frc_thr[idx]); // Relative convergence
      if(dbg_lvl == dbg_old){
	(void)std::fprintf(stderr,"%4ld %4ld %7.3f %7.5f %8.6f %8.6f %8.6f %8.6f %8.6f\n",
		      idx,itr_idx,dmt_ctr[idx]*1.0e6,K_MaB95_4,ryn_nbr_frc_thr_prx[idx],ryn_nbr_frc_thr[idx],wnd_frc_thr_prx[idx],wnd_frc_thr[idx],eps_crr);
      } // end if dbg
      if(itr_idx > 10){
	if(eps_crr > 1.0e-3) err_prn(prg_nm,sbr_nm,"Threshold friction speed not converging"); 
	// Do not die if ping pong occurs but convergence is fairly good (eps < 1.0e-3), e.g., for dmt_ctr=455.275
	wrn_prn(prg_nm,sbr_nm,"Threshold friction relative convergence only "+nbr2sng(eps_crr)+" for dmt_ctr = "+nbr2sng(dmt_ctr[idx]));
	eps_crr=0.0;
      } // endif
      itr_idx++;
    } // end loop over itr
    // Threshold friction Reynolds number may lie outside range of validity of IvW82 when D <~ 0.1 um
    if(ryn_nbr_frc_thr[idx] < 0.03 && dmt_ctr[idx] > 0.1e-6) std::cerr << prg_nm << ": WARNING " << sbr_nm << ": idx = " << idx << ", dmt_ctr = " << dmt_ctr[idx]*1.0e6 << " um, ryn_nbr_frc_thr = " << ryn_nbr_frc_thr[idx] << " < 0.03" << std::endl;
  } // end loop over sz
  return rcd;
} // end wnd_frc_thr_get()

int // O [rcd] Return success code
wnd_frc_thr_get_ShL00
(const long sz_nbr, // I [nbr] Size of arrays
 const prc_cmp *dmt_ctr, // I [m] Diameter at bin center
 const prc_cmp dns_prt, // I [kg m-3] Density of particle
 const prc_cmp dns_mdp, // I [kg m-3] Midlayer density
 prc_cmp *ryn_nbr_frc_thr, // O [frc] Threshold friction Reynolds number
 prc_cmp *ryn_nbr_frc_thr_prx, // O [frc] Threshold friction Reynolds number approximation
 prc_cmp *wnd_frc_thr, // O [m s-1] Threshold friction speed
 prc_cmp *wnd_frc_thr_prx) // O [m s-1] Threshold friction speed approximation
{
  /* Purpose: Compute threshold friction speed using method of Shao and Lu (2000)
     Routine also returns implied threshold friction Reynolds number for diagnostic
     fxm: finish this routine */
  // Output
  int rcd(0); // O [rcd] Return success code
  // Local
  long idx; // Counting index

  // Main code
  const std::string sbr_nm("wnd_frc_thr_get_ShL00"); // Subroutine name
  const std::string prg_nm(prg_nm_get()); // Program name
  const unsigned short dbg_lvl(dbg_lvl_get()); // Debugging level

  // Solution for threshold friction speed
  for(idx=0;idx<sz_nbr;idx++){
    
  } // end loop over idx
  if(dbg_lvl > 0){rcd+=rcd;}
  return rcd;
} // end wnd_frc_thr_get_ShL00()

int wnd_frc_thr_slt_get
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *dns_mdp, // I [kg m-3] Midlayer density
 prc_cmp *wnd_frc_thr_slt) // O [m s-1] Threshold friction speed for saltation
{
  // Purpose: Return dry threshold friction velocity for saltation
  // Dependencies: <phys_cst.hh>
  using phc::grv_sfc_mean; // (9.80665) [m s-2] Mean gravitational acceleration at Earth's surface
  // Local
  const prc_cmp dmt_slt_opt(75.0e-6); // [m] Optimal diameter for saltation IvW82 p. 117 Fgr. 8, Pye87 p. 31, MBA97 p. 4388, SRL96 (2)
  const prc_cmp dns_slt(2650.0); // [kg m-3] Density of optimal saltation particles MBA97 p. 4388
  int rcd(0); // O [rcd] Return success code
  long lon_idx; // Counting index for lon
  prc_cmp ryn_nbr_frc_thr_prx_opt; // [frc] Threshold friction Reynolds number approximation for optimal size
  prc_cmp ryn_nbr_frc_thr_opt_fnc; // [frc] Threshold friction Reynolds factor for saltation calculation
  prc_cmp dns_fct; // Density ratio factor for saltation calculation
  prc_cmp icf_fct; // Interpartical cohesive forces factor for saltation calculation
  prc_cmp tmp1; // Factor in saltation computation
  // Main Code
  // Initialize some variables
  
  // MaB95 pzn. for Re*t(D_opt) circumvents iterative solution
  ryn_nbr_frc_thr_prx_opt=0.38+1331.0*std::pow(PRC_CMP(100.0)*dmt_slt_opt,PRC_CMP(1.56)); // [frc] "B" MaB95 p. 16417 (5)
  // Given Re*t(D_opt), compute time independent factors contributing to u*t
  icf_fct=1.0+6.0e-07/(dns_slt*grv_sfc_mean*std::pow(dmt_slt_opt,PRC_CMP(2.5))); // IvW82 p. 115 (6) MaB95 p. 16417 (4) Interparticle cohesive forces
  dns_fct=dns_slt*grv_sfc_mean*dmt_slt_opt; // IvW82 p. 115 (6) MaB95 p. 16417 (4)
  assert(ryn_nbr_frc_thr_prx_opt >= 0.03);
  if(ryn_nbr_frc_thr_prx_opt < 10.0){
    ryn_nbr_frc_thr_opt_fnc=-1.0+1.928*std::pow(ryn_nbr_frc_thr_prx_opt,PRC_CMP(0.0922)); // IvW82 p. 114 (3), MaB95 p. 16417 (6)
    ryn_nbr_frc_thr_opt_fnc=0.1291*0.1291/ryn_nbr_frc_thr_opt_fnc;
  }else{
    ryn_nbr_frc_thr_opt_fnc=1.0-0.0858*std::exp(-0.0617*(ryn_nbr_frc_thr_prx_opt-10.0)); // IvW82 p. 114 (3), MaB95 p. 16417 (7)
    ryn_nbr_frc_thr_opt_fnc=0.120*0.120*ryn_nbr_frc_thr_opt_fnc*ryn_nbr_frc_thr_opt_fnc; // [frc]
  } // endif
  // This method minimizes number of square root computations performed
  tmp1=std::sqrt(icf_fct*dns_fct*ryn_nbr_frc_thr_opt_fnc);
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    wnd_frc_thr_slt[lon_idx]=tmp1/std::sqrt(dns_mdp[lon_idx]); // [m s-1] Threshold friction velocity for saltation dry ground 
  } // end loop over lon
  
  return rcd;
} // end wnd_frc_thr_slt_get()

int wnd_rfr_thr_slt_get
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *wnd_frc, // I [m s-1] Surface friction velocity
 const prc_cmp *wnd_frc_thr_slt, // I [m s-1] Threshold friction speed for saltation
 const prc_cmp *wnd_mdp, // I [m s-1] Surface layer mean wind speed
 const prc_cmp *wnd_rfr, // I [m s-1] Wind speed at reference height
 prc_cmp *wnd_rfr_thr_slt) // O [m s-1] Threshold 10 m wind speed for saltation
{
  // Purpose: Infer threshold 10 m wind speed for saltation from threshold friction speed for saltation and 10 m wind speed
  // Dependencies: 
  // Local
  int rcd(0); // O [rcd] Return success code
  long lon_idx; // Counting index for lon
  // Main Code
  // Initialize some variables
  
  // Compute threshold horizontal wind speed at reference height
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    // A more complicated procedure would recompute mno_lng for wnd_frc_thr,
    // and then integrate vertically from rgh_mmn+hgt_zpd to hgt_rfr
    // wnd_crc_fct is (1/k)*[ln(z-D)/z0 - psi(zeta2) + psi(zeta1)]
    wnd_rfr_thr_slt[lon_idx]=wnd_frc_thr_slt[lon_idx]*wnd_rfr[lon_idx]/wnd_frc[lon_idx]; // [m s-1]
  } // end loop over lon
  
  if(wnd_mdp[0] == wnd_mdp[0]){;} // CEWU Compiler Error Warning Usage
  return rcd;
} // end wnd_rfr_thr_slt_get()

int wnd_frc_slt_get
(const long lon_nbr, // I [nbr] Size of arrays
 const prc_cmp *wnd_frc_mbl, // I [m s-1] Surface friction velocity
 prc_cmp *wnd_frc_slt, // O [m s-1] Saltating friction velocity
 prc_cmp *wnd_frc_slt_dlt, // O [m s-1] Friction velocity increase from saltation
 const prc_cmp *wnd_rfr_mbl, // I [m s-1] Wind speed at reference height
 const prc_cmp *wnd_rfr_thr_slt) // I [m s-1] Threshold 10 m wind speed for saltation
{
  /* Purpose: Compute the saltating friction velocity
     Saltation increases friction speed by roughening surface, AKA "Owen's effect"
     This acts as a positive feedback to the friction speed
     GMB98 parameterized this feedback in terms of 10 m windspeeds */
  // Dependencies: 
  // Local
  int rcd(0); // O [rcd] Return success code
  long lon_idx; // Counting index for lon
  // Main Code
  // Initialize some variables
  
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    // Saltating friction velocity equals friction velocity if no saltation
    wnd_frc_slt[lon_idx]=wnd_frc_mbl[lon_idx]; // [m s-1] Saltating friction velocity
    wnd_frc_slt_dlt[lon_idx]=0.0; // [m s-1] Friction velocity increase from saltation
    if(wnd_rfr_mbl[lon_idx] > wnd_rfr_thr_slt[lon_idx]){
      /* Saltation roughens the boundary layer, AKA "Owen's effect"
	 GMB98 p. 6206 Fig. 1 shows observed/computed u* dependence on observed U(1 m)
	 GMB98 eqn(12) p. 6209 has u* in cm s-1 and U, Ut in m s-1, personal communication, D. Gillette, 19990529
	 With everything in MKS, the 0.3 coefficient in GMB98 (12) becomes 0.003 
	 Increase in friction velocity due to saltation varies as square of 
	 difference between reference wind speed and reference threshold speed 
      */
      wnd_frc_slt_dlt[lon_idx]=0.003*(wnd_rfr_mbl[lon_idx]-wnd_rfr_thr_slt[lon_idx])*(wnd_rfr_mbl[lon_idx]-wnd_rfr_thr_slt[lon_idx]); // [m s-1] Friction velocity increase from saltation GMB98 p. 6209
      wnd_frc_slt[lon_idx]=wnd_frc_mbl[lon_idx]+wnd_frc_slt_dlt[lon_idx]; // [m s-1] Saltating friction velocity
    } // endif wnd_frc_mbl > wnd_frc_thr_slt
  } // end loop over lon
  
  return rcd;
} // end wnd_frc_slt_get()

int frc_thr_ncr_drg_get
(const long lon_nbr, // I [nbr] Size of arrays
 prc_cmp *frc_thr_ncr_drg, // O [frc] Factor by which surface roughness increases threshold friction velocity
 prc_cmp *wnd_frc_fsh_frc, // O [frc] Efficient fraction of wind friction
 const prc_cmp rgh_mmn_mbl, // I [m] Roughness length momentum for erodible surfaces
 const prc_cmp rgh_mmn_smt) // I [m] Smooth roughness length
{
  /* Purpose: Compute factor by which surface roughness increases threshold friction velocity
     This parameterization is based on MaB95 and GMB98
     Input rgh_mmn_mbl will be changed only if it exceeds maximum permissible value */

  /* Drag partition formula was only evaluated for rgh_mmn_mbl is < 0.1 cm
     Drag partition formula fails for rgh_mmn_mbl > rgh_mmn_mbl_max = 0.5 cm
     Hence use max(rgh_mmn_mbl,rgh_mmn_mbl_max) in formula and manually set
     output wnd_frc_fsh_frc to an arbitrarily small value so that threshold
     velocity is unattainable.
     rgh_mmn_mbl > rgh_mmn_mbl_max is reasonable over many sheltered surfaces
     such as forests, e.g., where rgh_mmn ~ 1.0 m.
     Hence changing rgh_mmn_mbl would alter other turbulent fluxes (SH, LH). */

  /* The ratio of smooth to bulk friction velocity, AKA "efficient fraction", accounts for partition of drag over rough surfaces
     Non-erodible objects, e.g., pebbles, sticks, bushes, grasses, generally determine the measured roughness length rgh_mmn_dps
     The roughness length determined by the erodible particles themselves, called the "smooth roughness length", z0ms, is about D/30 (MaB95 p. 16425, SeP97 p. 858)
     For coarse sand particles D ~ 200 um, so z0ms ~5 um = 5.0e-4 cm = 5.0e-6 m
     Values of z0ms for natural conditions:
     MaB95 p. 16426 recommend z0ms ~ 10 um = 1.0e-3 cm = 10.0e-6 m
     GMB98 p. 6207 found z0ms ~ 5 um = 5.0e-4 cm = 5.0e-6 m */

  /* Agreement with GMB98 data is an important test of boundary layer processes
     For their Site 5010 at Owen's Lake California:
     Measured z0NS = 0.000098 m, u*t = 0.42 m s-1
     Modeled u*t = 0.3827 m s-1
     19990521: Simulations of u* and u*t agree well with these values using
     mie --rgh_mmn_mbl=98.0e-6 --rgh_mmn_smt=5.0e-6 --prs_mdp=95000.0 --prs_ntf=97500.0 */

  // Dependencies: 
  // Local
  const std::string prg_nm(prg_nm_get()); // Program name
  const std::string sbr_nm("frc_thr_ncr_drg_get"); // Subroutine name
  int rcd(0); // O [rcd] Return success code
  long lon_idx; // Counting index for lon
  // Main Code
  const prc_cmp rgh_mmn_mbl_max(0.005); // [m] Roughness length momentum for erodible surfaces, maximum value allowed MaB95 p. 16420, GMB98 p. 6205
  const prc_cmp wnd_frc_fsh_frc_min(1.0e-3); // O [frc] Efficient fraction of wind friction
  if(rgh_mmn_mbl >= rgh_mmn_mbl_max){
    wrn_prn(prg_nm,sbr_nm,"Reducing rgh_mmn_mbl from user-specified value of "+nbr2sng(rgh_mmn_mbl)+" m to maximum permissible value of "+nbr2sng(rgh_mmn_mbl_max)+" m for drag-partition computation only. Change should only affect dust production threshold velocities and turbulent fluxes on surfaces completely dominated by non-erodible surfaces.");
  } // endif
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    wnd_frc_fsh_frc[lon_idx]= // [frc] Efficient fraction of wind friction MaB95 p. 16420, GMB98 p. 6207
      +1.0-std::log(rgh_mmn_mbl/rgh_mmn_smt)/std::log(0.35*std::pow(PRC_CMP(0.1)/rgh_mmn_smt,PRC_CMP(0.8)));
    if(wnd_frc_fsh_frc[lon_idx] <= 0.0){
      wrn_prn(prg_nm,sbr_nm,"Setting wnd_frc_fsh_frc to "+nbr2sng(wnd_frc_fsh_frc_min)+" instead of computed value of "+nbr2sng(wnd_frc_fsh_frc)+" to avoid invalid range of MaB95 parameterization.");
      wnd_frc_fsh_frc[lon_idx]=wnd_frc_fsh_frc_min; // [frc] Efficient fraction of wind friction MaB95 p. 16420, GMB98 p. 6207
    } // endif
    assert(wnd_frc_fsh_frc[lon_idx] <= 1.0);
    frc_thr_ncr_drg[lon_idx]=1.0/wnd_frc_fsh_frc[lon_idx]; // [frc] Factor by which surface roughness increases threshold friction velocity
  } // end loop over lon
  
  return rcd;
} // end frc_thr_ncr_drg_get()

int frc_thr_ncr_wtr_get
(const long lon_nbr, // I [nbr] Size of arrays
 prc_cmp *frc_thr_ncr_wtr, // O [frc] Factor by which moisture increases threshold friction velocity
 prc_cmp *vwc_thr, // O [m3 m-3] Threshold volumetric water content to affect mobilization
 const prc_cmp *mss_frc_cly, // I [frc] Mass fraction of clay
 const prc_cmp *vwc_sfc) // I [frc] Volumetric water content
{
  /* Purpose: Compute factor by which soil moisture increases threshold friction velocity 
     This parameterization is based on FMB99 */
  // Dependencies: 
  // Local
  int rcd(0); // O [rcd] Return success code
  long lon_idx; // Counting index for lon
  // Main Code
  for(lon_idx=0;lon_idx<lon_nbr;lon_idx++){
    vwc_thr[lon_idx]= // [m3 m-3] Threshold volumetric water content to affect mobilization 
      mss_frc_cly[lon_idx]*(0.17+0.14*mss_frc_cly[lon_idx]); // FMB99 p. 154 (14)
    if(vwc_sfc[lon_idx] < vwc_thr[lon_idx]) frc_thr_ncr_wtr[lon_idx]=1.0; else frc_thr_ncr_wtr[lon_idx]=std::sqrt(1.0+1.21*std::pow(PRC_CMP(100.0)*(vwc_sfc[lon_idx]-vwc_thr[lon_idx]),PRC_CMP(0.68))); // FMB99 p. 155 (15)
  } // end loop over lon
  
  return rcd;
} // end frc_thr_ncr_wtr_get()
