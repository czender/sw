// $Id$ 

// Implementation (declaration) of phase function utilities

/* Copyright (C) 1997--present Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

#include <phz_fnc.hh> // Phase functions

// Friendly functions begin

// Friendly functions end
// Static members begin

// Static members end
// Static member functions begin

// Static member functions end
// Public member functions begin

// Public member functions end
// Private member functions begin

// Private member functions end
// Global functions with C++ linkages begin

prc_cmp // O [frc] Legendre polynomial of order n+1
lgn_Pnp1 // [fnc] Legendre polynomial of order n+1
(int lgn_rdr_n, // [nbr] Legendre polynomial order (n >= -1)
 prc_cmp abc) // [frc] Abscissa [-1,1]
{
  /* Purpose: Return Legendre polynomial of order n+1 (not n)
     Method follows Mic63 p. 238 (9)
     See also free sample code at http://www.storage-b.com/c/18
     That code uses namespaces and inlines and is more efficient */
  prc_cmp nsw; // [frc] P_{n+1}
  prc_cmp trm_Pn; // [frc] Term containing P_{n}
  prc_cmp trm_Pnm1; // [frc] Term containing P_{n-1}
  assert(lgn_rdr_n >= -1);
  if(lgn_rdr_n == -1) return 0.0; // P_{0} is constant 1.0
  if(lgn_rdr_n == 0) return abc; // P_{1} is x
  trm_Pnm1=lgn_rdr_n*lgn_Pnp1(lgn_rdr_n-2,abc);
  trm_Pn=(2*lgn_rdr_n+1.0)*lgn_Pnp1(lgn_rdr_n-1,abc)*abc;
  nsw=(trm_Pn-trm_Pnm1)/(lgn_rdr_n+1.0);
  return nsw;
} // end lgn_Pnp1()

prc_cmp // O [frc] First derivative of Legendre polynomial of order n+1
lgn_frs_drv_Pnp1 // [fnc] First derivative of Legendre polynomial of order n+1
(int lgn_rdr_n, // [nbr] Legendre polynomial order (n >= -1)
 prc_cmp abc) // [frc] Abscissa [-1,1]
{
  // Purpose: Return first derivative of Legendre polynomial of order n+1 (not n)
  // Method follows Mic63 p. 238 (10)
  prc_cmp nsw; // [frc] First derivative of P_{n+1}
  prc_cmp trm_Pn; // [frc] Term containing first derivative of P_{n}
  prc_cmp trm_Pnm1; // [frc] Term containing first derivative of P_{n-1}
  assert(lgn_rdr_n >= -1);
  if(lgn_rdr_n == -1) return 0.0; // First derivative of P_{0} is constant 0.0
  if(lgn_rdr_n == 0) return 1.0; // First derivative of P_{1} is constant 1.0
  trm_Pnm1=(lgn_rdr_n+1.0)*lgn_frs_drv_Pnp1(lgn_rdr_n-2,abc);
  trm_Pn=(2*lgn_rdr_n+1.0)*lgn_frs_drv_Pnp1(lgn_rdr_n-1,abc)*abc;
  nsw=(trm_Pn-trm_Pnm1)/lgn_rdr_n;
  return nsw;
} // end lgn_frs_drv_Pnp1()

prc_cmp // O [frc] Second derivative of Legendre polynomial of order n+1
lgn_scn_drv_Pnp1 // [fnc] Second derivative of Legendre polynomial of order n+1
(int lgn_rdr_n, // [nbr] Legendre polynomial order (n >= 0)
 prc_cmp abc) // [frc] Abscissa [-1,1]
{
  // Purpose: Return second derivative of Legendre polynomial of order n+1 (not n)
  // Method follows Mic63 p. 238 (11)
  prc_cmp nsw; // [frc] Second derivative of P_{n+1}
  prc_cmp trm_Pn; // [frc] Term containing second derivative of P_{n}
  prc_cmp trm_Pnm1; // [frc] Term containing second derivative of P_{n-1}
  assert(lgn_rdr_n >= 0);
  if(lgn_rdr_n == 0) return 0.0; // Second derivative of P_{1} is constant 0.0
  if(lgn_rdr_n == 1) return 3.0; // Second derivative of P_{2} is constant 3.0
  trm_Pnm1=(lgn_rdr_n+2.0)*lgn_scn_drv_Pnp1(lgn_rdr_n-2,abc);
  trm_Pn=(2*lgn_rdr_n+1.0)*lgn_scn_drv_Pnp1(lgn_rdr_n-1,abc)*abc;
  nsw=(trm_Pn-trm_Pnm1)/(lgn_rdr_n-1.0);
  return nsw;
} // end lgn_scn_drv_Pnp1()

int // O [enm] Return success code
ngl_grd_get // [fnc] Create angular grid
(const long ngl_nbr, // I [nbr] Number of angles
 const std::string ngl_sng, // I [sng] Angle grid type
 prc_cmp * const ngl, // O [rdn] Angle
 prc_cmp * const ngl_dgr, // O [dgr] Angle degrees
 prc_cmp * const ngl_dlt, // O [rdn] Width of angle bin
 prc_cmp * const ngl_wgt) // O [rdn] Weight of angle bin
{
  /* Purpose: Create angular grid 
     ngl_wgt is _not_ normalized to unity 
     ngl_wgt units are [rdn] not [frc] (so integrals have correct dimensions)
     ngl_wgt is d(theta) for linear grids and \sum ngl_dlt = \sum ngl_wgt = pi
     ngl_wgt is quadrature weight Lobatto grids and \sum wgt = pi (or 2)
     ngl_wgt and ngl_dlt are interchange-able for linear grids
     For quadrature grids, ngl_dlt = spacing between angle interfaces
     Quadrature spacing is for diagnostic purposes only---not for integration 

     NB: Indices printed with 'F' or 'ftn' notation are Fortran-convention (1-based)
     This helps maintain consistency with Lobatto quadrature literature 
     Output arrays are all C-style so angle one has subscript zero */

  long ngl_idx; // [idx] Counting index for angle
  const std::string sbr_nm("ngl_grd_get"); // [sng] Subroutine name
  const unsigned short dbg_lvl(dbg_lvl_get()); // Debugging level
  int rcd(0); // [enm] Return code

  if(ngl_sng == "regular"){
    /* First and last angles are 0 and pi, respectively 
       They are on left and right edges of their respective bins
       (Equivalently, they are centered at 0 and pi and weigh half as much)
       Remaining angles are centered in bins and weigh twice as much 
       C-grid? */
    const double ngl_dlt_scl=mth::cst_M_PIl/(2.0*(ngl_nbr-1)); // [rdn] Width of angle bin
    ngl[0]=0.0; // [rdn] Scattering angle
    ngl[2*ngl_nbr-2]=mth::cst_M_PIl; // [rdn] Scattering angle
    ngl_dlt[0]=ngl_dlt[2*ngl_nbr-2]=0.5*ngl_dlt_scl; // [rdn] Width of angle bin
    for(ngl_idx=1;ngl_idx<2*ngl_nbr-2;ngl_idx++){ // NB: Starts at one
      ngl[ngl_idx]=ngl_dlt_scl*ngl_idx; // [rdn] Scattering angle
      ngl_dlt[ngl_idx]=ngl_dlt_scl; // [rdn] Width of angle bin
    } // end loop over ngl
    for(ngl_idx=0;ngl_idx<2*ngl_nbr-1;ngl_idx++){
      // Weight is bin-width for regular (linear) grids
      ngl_wgt[ngl_idx]=ngl_dlt[ngl_idx]; // [rdn] Weight of angle bin
    } // end loop over ngl
  }else if(ngl_sng == "lobatto"){ // ngl_grd!=linear
    /* Compute and use Gauss-Lobatto quadrature angles and weights
       Michels, H. H. (1963): Abscissas and weight coefficients for Lobatto quadrature, Math. Comput., 17(83), 237--244
       Mic63 method is fairly straightforward
       Wis771 recommend performing quadrature on [0,pi]
       This requires some thought and weight re-normalization, done at the end
       Lobatto abscissae and weights for L<=7: http://mathworld.wolfram.com/LobattoQuadrature.html
       NB: Lobatto abscissae and weights are symmetric about the origin
       Gauss-Lobatto and Gauss-Legendre quadrature have weighting functions W(x)=1
       Testing:
       mie --dbg=3 --ngl_sng=lobatto > ~/foo 2>&1 */
    using mth::cst_M_PIl; // (3.1415926535897932384626433832795029L) [frc] 3
    // NB: Lobatto quadrature abscissa converge remarkably smoothly and rapidly
    const prc_cmp eps_max(1.0e-14); // [frc] Relative accuracy for convergence
    const long lbb_nbr(2*ngl_nbr-1); // [nbr] Lobatto quadrature order

    double lgn_Lm1; // [frc] Legendre polynomial of degree lbb_nbr-1 evaluated at Lobatto abscissae
    gsl_sf_result nsw_dbl; // [frc] GSL result structure
    long itr_idx; // [idx] Counting index
    long lbb_idx; // [idx] Lobatto index [0..lbb_nbr-1]
    long lbbp1_idx; // [idx] Lobatto index [1..lbb_nbr]
    prc_cmp eps_crr; // [frc] Current relative accuracy
    prc_cmp ngl_cos_old; // [frc] Old angle cosine
    std::valarray<prc_cmp> ngl_cos(2*ngl_nbr-1); // [frc] Angle cosine

    // Intialize constant factors
    const double pi_sqr(cst_M_PIl*cst_M_PIl); // [frc] Pi-squared
    const double fct_tmp1((pi_sqr-4.0)/(4.0*pi_sqr)); // [frc] Factor Mic63 p. 238 (7)
    const double fct_tmp2((lbb_nbr-0.5)*(lbb_nbr-0.5)); // [frc] Factor Mic63 p. 238 (7)
    const double fct_tmp3(std::sqrt(fct_tmp1+fct_tmp2)); // [frc] Factor Mic63 p. 238 (7)

    for(lbb_idx=1;lbb_idx<=lbb_nbr-2;lbb_idx++){ // NB: loop over interior angles
      // lbbp1_idx is index 'k' in Mic63
      lbbp1_idx=lbb_idx+1; // [idx] Legendre index [1..lbb_nbr]
      /* Find roots of Bessel function of order 1 
	 20080121: GSL answers agree to all decimals with those at
	 http://webcomputing.bio.bas.bg/webMathematica/webComputing/BesselZeros.jsp */
      rcd=gsl_sf_bessel_zero_J1_e(static_cast<unsigned int>(lbb_idx),&nsw_dbl);
      if(dbg_lvl == dbg_old) std::cout << sbr_nm+" reports Lobatto j(ftn" << lbbp1_idx << ") = " << nsw_dbl.val << std::endl;
      if(rcd != 0) std::cout << "WARNING: rcd = " << rcd << ", nsw_dbl.err = " << nsw_dbl.err << std::endl;
      // Initialize accuracy and counter
      eps_crr=eps_max+1.0; // [frc] Current relative accuracy
      itr_idx=1; // Counting index
      // Initial guess, usually good to 2-3+ decimal places
      ngl_cos[lbb_idx]=std::cos(nsw_dbl.val/fct_tmp3); // [frc] Angle cosine
      // Iterate until convergence
      if(dbg_lvl == dbg_vrb){
	(void)std::fprintf(stderr,"%4s %4s %16s %16s %16s %16s\n",
			   "idx ","itr ","    ngl_cos    ","    ngl_rdn    ","    ngl_dgr    ","      eps      ");
	(void)std::fprintf(stderr,"%4s %4s %16s %16s %16s %16s\n",
			   "ftn ","    ","      frc      ","      rdn      ","      dgr      ","      frc      ");
      } // end if dbg
      while(eps_crr > eps_max){
	// Save angle for convergence test
	ngl_cos_old=ngl_cos[lbb_idx]; // [rdn] 
	// Update angle using Newton-Raphson
	ngl_cos[lbb_idx]=ngl_cos_old-Legendre::frs_drv_Pn(lbb_nbr-1,ngl_cos_old)/Legendre::scn_drv_Pn(lbb_nbr-1,ngl_cos_old); // [frc] Mic63 p. 238 (8)
	eps_crr=PRC_CMP_ABS((ngl_cos[lbb_idx]-ngl_cos_old)/max_cpv(std::fabs(ngl_cos[lbb_idx]),1.0e-12)); // Relative convergence
	if(dbg_lvl >= dbg_vrb){
	  (void)std::fprintf(stderr,"%4ld %4ld %16.12f %16.12f %16.12f %16.12f\n",
			     lbbp1_idx,itr_idx,ngl_cos[lbb_idx],std::acos(ngl_cos[lbb_idx]),180.0*std::acos(ngl_cos[lbb_idx])/cst_M_PIl,eps_crr);
	} // end if dbg
	itr_idx++;
      } // end loop over itr

      // Lobatto quadrature abscissa has converged for this angle, continue

    } // end loop over quadrature angles
    
    // Choose Lobatto quadrature type
    const bool lobatto_quadrature_Wis771_A1(false);
    const bool lobatto_quadrature_Wis771_A2(true);
    double lobatto_interval_span;
    if (lobatto_quadrature_Wis771_A1){
      // Wis771 (A1): Lobatto quadrature in cos(theta) on [-1,1]: ngl=acos(ngl_cos), \sum wgt = 2
      // Not recommended for strongly peaked phase functions
      // Recommended for symmetric phase functions, e.g., Rayleigh scattering
      lobatto_interval_span=2.0;
    }else if(lobatto_quadrature_Wis771_A2){
      // Wis771 (A2): Lobatto quadrature in theta on [0,pi]: ngl=-pi*(ngl_cos-1)/2, \sum wgt = pi
      // Highly recommended for strongly peaked phase functions
      lobatto_interval_span=cst_M_PIl;
    }else{
      err_prn(sbr_nm,"Unknown quadrature type");
    } // endif

    // Endpoint angles always 0 and 180 degrees for Lobatto quadrature
    ngl[0]=0.0; // [rdn] Angle
    ngl[lbb_nbr-1]=cst_M_PIl; // [rdn] Angle

    // Endpoint weights depend on quadrature span 
    ngl_wgt[0]=ngl_wgt[lbb_nbr-1]=lobatto_interval_span/(lbb_nbr*(lbb_nbr-1)); // [rdn] Weight of angle bin Mic63 (4)

    // Interior angles and weights depend on Lobatto quadrature type
    for(lbb_idx=1;lbb_idx<=lbb_nbr-2;lbb_idx++){ // NB: Loop over interior angles
      if(lobatto_quadrature_Wis771_A1){
	// Forward clustering in cos(theta) so angles evenly distributed in theta
	ngl[lbb_idx]=acos(ngl_cos[lbb_idx]); // [rdn] Angle
      }else if(lobatto_quadrature_Wis771_A2){
	// Forward clustering in theta
	ngl[lbb_idx]=-cst_M_PIl*(ngl_cos[lbb_idx]-1.0)/2.0; // [rdn] Angle
      }else{
	err_prn(sbr_nm,"Unknown quadrature type");
      } // endif

      // Mic63 p. 237 (3)
      // Evaluate Legendre polynomial of degree L-1 at abscissae
      lgn_Lm1=Legendre::Pn(lbb_nbr-1,ngl_cos[lbb_idx]);
      ngl_wgt[lbb_idx]=lobatto_interval_span/(lbb_nbr*(lbb_nbr-1.0)*lgn_Lm1*lgn_Lm1); // [rdn] Weight of angle bin Mic63 (4)

    } // end loop over quadrature angles

    ngl_dlt[0]=0.5*ngl[1]; // [rdn] Width of angle bin
    ngl_dlt[lbb_nbr-1]=cst_M_PIl-0.5*(cst_M_PIl+ngl[lbb_nbr-2]); // [rdn] Width of angle bin
    for(lbb_idx=1;lbb_idx<lbb_nbr-1;lbb_idx++){ // NB: Short grid
      // Define diagnostic ngl_dlt for quadrature angles
      ngl_dlt[lbb_idx]=0.5*(ngl[lbb_idx+1]+ngl[lbb_idx])-0.5*(ngl[lbb_idx]+ngl[lbb_idx-1]); // [rdn] Width of angle bin
    } // end loop over quadrature angles

  }else{ // ngl_grd!=linear
    err_prn(sbr_nm,"Angular grid ngl_sng == "+ngl_sng+" not implemented yet");
  } // ngl_grd!=linear

  // Compute diagnostics for angle grid
  prc_cmp ngl_dlt_ttl(0.0); // [rdn] Angle grid span
  prc_cmp ngl_wgt_ttl(0.0); // [frc] Angle grid weight sum
  using mth::cst_M_PIl; // (3.1415926535897932384626433832795029L) [frc] 3
  for(ngl_idx=0;ngl_idx<2*ngl_nbr-1;ngl_idx++){
    ngl_dlt_ttl+=ngl_dlt[ngl_idx]; // [rdn] Angle grid span
    ngl_wgt_ttl+=ngl_wgt[ngl_idx]; // [rdn] Angle grid weight sum
    ngl_dgr[ngl_idx]=ngl[ngl_idx]*180.0/mth::cst_M_PIl; // [dgr] Scattering angle
  } // end loop over ngl

  if(dbg_lvl == dbg_vrb){
    /* 20080123 Verified Lobatto abscissae and weights (for L=5) exactly match
       http://mathworld.wolfram.com/LobattoQuadrature.html */
    std::cout << "ngl_dlt_ttl = " << ngl_dlt_ttl << std::endl; // [rdn] Angle grid span
    std::cout << "ngl_wgt_ttl = " << ngl_wgt_ttl << std::endl; // [rdn] Angle grid weight sum
    (void)std::fprintf(stderr,"%4s %16s %16s %16s %16s %16s\n",
		       "idx ","    ngl_cos    ","    ngl_rdn    ","    ngl_dgr    ","    ngl_dlt    ","      wgt      ");
    (void)std::fprintf(stderr,"%4s %16s %16s %16s %16s %16s\n",
		       "ftn ","      frc      ","      rdn      ","      dgr      ","      dgr      ","      rdn      ");
    for(ngl_idx=0;ngl_idx<2*ngl_nbr-1;ngl_idx++){
      // lbbp1_idx is index 'k' in Mic63
      (void)std::fprintf(stderr,"%4ld %16.12f %16.12f %16.12f %16.12f %16.12f\n",
			 ngl_idx+1,std::cos(ngl[ngl_idx]),ngl[ngl_idx],ngl_dgr[ngl_idx],180.0*ngl_dlt[ngl_idx]/cst_M_PIl,ngl_wgt[ngl_idx]);
    } // end loop over quadrature angles

  } // end if dbg
  
  return rcd; // [enm] Return code
} // end ngl_grd_get()

int // O [enm] Return success code
phz_fnc_mdl // [fnc] Phase function diagnostics module
(const int dmn_nbr_max, // I [nbr] Maximum number of dimensions allowed in single variable in output file 
 const int nc_out, // I [fl] netCDF file for output 
 const long lgn_nbr, // I [nbr] Order of phase function Legendre expansion
 const long ngl_nbr, // I [nbr] Number of polar angles in one hemisphere
 const long wvl_idx_dbg, // I [idx] Debugging wavelength bin
 const long wvl_nbr, // I [nbr] Number of output wavelength bands
 const prc_cmp ngl_dbg_dgr, // I [dgr] Debugging angle
 const prc_cmp * const ngl, // I [rdn] Scattering angle
 const prc_cmp * const ngl_dgr, // I [dgr] Angle degrees
 const prc_cmp * const ngl_dlt, // I [rdn] Width of angle bin
 const prc_cmp * const ngl_wgt, // I [rdn] Weight of angle bin
 const a2d_cls<prc_cmp> &phz_fnc_ffc, // I [sr-1] Effective phase function (weighted over size distribution and sub-bands)
 const prc_cmp * const phz_fnc_dgn, // I [sr-1] Phase function at diagnostic size, wavelength
 const prc_cmp * const plz_dgn) // I [frc] Degree of linear polarization at diagnostic size, wavelength
{
  /* Purpose: Verify correct normalization of phase function
     Compute phase function Legendre expansion

     NB: lgn_nbr=8 means expand out to eighth order Legendre polynomials
     Number of terms in algebraic expansion is phase function order plus one
     Algebraic expansion always includes a zeroth term which equals 1.0.
     This routine does not need to store the zeroth term, since it should be 1.0
     The difference between the zeroth term and 1.0 indicates the expansion accuracy
     As does the difference between the first term and the asymmetry parameter
     DISORT numbers phase function moments from 0:mmn_nbr, i.e., 0:lgn_nbr
     So most intuitive procedure is to store moments 0:mmn_nbr, i.e., 0:lgn_nbr
     Subscripts appear correctly in C-notation, not in (default) Fortran notation

     This diagnostic procedure assumes azimuthally invariant phase functions
     fxm: this procedure delete()'s memory allocated in main (e.g., ngl) even though pointers are passed as const
     Mie program has somewhat haphazard definition of angular grid:
     Currently, angular grid is defined to have 2*ngl_nbr-1 angles, because...
     Mie solver uses angular symmetry to reduce computational loops to [1..ngl_nbr]
     However, this makes interpretation of ngl_nbr confusing
     Angular grid should be made similar to longitude grid, well-defined interfaces, midpoints, etc.
     Currently accepts linear and Lobatto angular grids
     Improvements that could be made:
     1. Differentiate between ngl_nbr used in Mie solvers, and ngl_nbr 

     Generate phase function diagnostics on 1 degree grid
     mie -dbg --ngl_nbr=91 --lgn_nbr=16 --cmp_prt=sulfate --psd_typ=lognormal --wvl_mnm=0.495 --wvl_mxm=0.505 --wvl_nbr=2 --bnd_nbr=1 --sz_grd=log --sz_mnm=1.0e-3 --sz_mxm=5.0 --sz_nbr=400 --rds_nma=0.4 --gsd_anl=1.6 ${DATA}/mie/mie.nc
     ncks -C -d wvl,0 -v asm_prm,asm_prm_lgn,lgn_xpn_cff,ngl_dgr,phz_fnc_dgn,phz_fnc_ffc,phz_fnc_bsf,phz_fnc_fsf,phz_fnc_lgn,phz_fnc_ntg,phz_fnc_lgn_ntg ${DATA}/mie/mie.nc | m

     Phase function moments for single particle/wavelength:
     mie -dbg --ngl_nbr=91 --lgn_nbr=16 --cmp_prt=sulfate --psd_typ=lognormal --wvl_mnm=0.495 --wvl_mxm=0.505 --wvl_nbr=1 --bnd_nbr=1 --sz_grd=lnr --sz_mnm=0.999 --sz_mxm=1.001 --sz_nbr=1 --rds_nma=1.0 --gsd_anl=1.6 ${DATA}/mie/mie.nc

     View phase function: ${HOME}/idl/mie.pro:phz_fnc(),aer_gph()
     aer_gph,fl_nm='/data/zender/mie/mie.nc' */

  long rcd(0); // Return success code
  // Local
  long ngl_idx; // [idx] Counting index for angle
  long lgn_idx; // [idx] Counting index for Legendre polynomial
  long wvl_idx; // [idx] Counting index for wavelength
  const std::string sbr_nm("phz_fnc_mdl"); // [sng] Name of subroutine
  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Entering...");

  // Compute diagnostics for angle grid
  std::valarray<prc_cmp> ngl_min_dgr(2*ngl_nbr-1); // [dgr] Minimum scattering angle in bin
  std::valarray<prc_cmp> ngl_max_dgr(2*ngl_nbr-1); // [dgr] Maximum scattering angle in bin
  ngl_min_dgr[0]=0.0;
  ngl_max_dgr[0]=0.5*(ngl_dgr[0]+ngl_dgr[1]);
  if(ngl_nbr > 1) ngl_max_dgr[2*ngl_nbr-2]=180.0;
  if(ngl_nbr > 1) ngl_min_dgr[2*ngl_nbr-2]=0.5*(ngl_dgr[2*ngl_nbr-3]+ngl_dgr[2*ngl_nbr-2]);
  for(ngl_idx=1;ngl_idx<2*ngl_nbr-2;ngl_idx++){ // NB: loop ends at 2*ngl_nbr-3
    ngl_min_dgr[ngl_idx]=0.5*(ngl_dgr[ngl_idx-1]+ngl_dgr[ngl_idx]); // [dgr] Minimum scattering angle in bin
    ngl_max_dgr[ngl_idx]=0.5*(ngl_dgr[ngl_idx]+ngl_dgr[ngl_idx+1]); // [dgr] Minimum scattering angle in bin
  } // end loop over ngl

  const long ngl_idx_dbg=vec_val2idx(ngl_dgr,2*ngl_nbr-1,ngl_dbg_dgr); // [idx] Angle bin for debugging

  /* Decompose phase function into Legendre polynomial basis
     Assume even number of terms in Legendre polynomial expansion
     Legendre polynomials are defined by their order, n
     Associated Legendre polynomials also have a degree, m, where m < n
     P_0(u) = 1
     P_1(u) = u
     P_2(u) = 0.5*(3u^2 - 1)
     P_3(u) = 0.5*(5u^3 - 3u)
     P_4(u) = 0.125*(35u^4 - 30 u^2 + 3) */

  a2d_cls<prc_cmp> lgn_xpn_cff(wvl_nbr,lgn_nbr+1); // [frc] Legendre polynomial expansion coefficients
  a2d_cls<prc_cmp> phz_fnc_lgn(wvl_nbr,2*ngl_nbr-1); // [sr-1] phase function Legendre expansion
  std::valarray<prc_cmp> asm_prm_lgn(CEWI_cpv,wvl_nbr); // [frc] Asymmetry parameter estimated from Legendre expansion

  /* Initialize expansion coefficients
     lgn_val, ngl_cos must be type double for GSL compatibility */
  std::valarray<int> lgn(lgn_nbr+1); // [idx] Order of Legendre expansion term
  std::valarray<double> lgn_val(lgn_nbr+1); // [frc] Legendre polynomial values
  std::valarray<double> ngl_cos(2*ngl_nbr-1); // [frc] Cosine of scattering angle
  std::valarray<prc_cmp> ngl_sin(2*ngl_nbr-1); // [frc] Sine of scattering angle
  lgn_xpn_cff=0.0; // [frc] Legendre polynomial expansion coefficients
  for(ngl_idx=0;ngl_idx<2*ngl_nbr-1;ngl_idx++){
    ngl_sin[ngl_idx]=static_cast<prc_cmp>(std::sin(ngl[ngl_idx])); // [frc] Sine of scattering angle
    ngl_cos[ngl_idx]=std::cos(ngl[ngl_idx]); // [frc] Cosine of scattering angle
  } // end loop over angle

  prc_cmp ngl_fct; // [frc] Factor common to given angle
  // For each angle where phase function is known...
  for(ngl_idx=0;ngl_idx<2*ngl_nbr-1;ngl_idx++){
    // Evaluate all orders of Legendre polynomial for this angle
    rcd+=gsl_sf_legendre_Pl_array(lgn_nbr,ngl_cos[ngl_idx],&lgn_val[0]); // [fnc] Legendre function
    ngl_fct=0.5*ngl_sin[ngl_idx]*ngl_wgt[ngl_idx]; // [frc] Factor common to given angle
    // For each wavelength...
    for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
      // For each order of Legendre polynomial for this wavelength...
      for(lgn_idx=0;lgn_idx<=lgn_nbr;lgn_idx++){
	// Wis771 shows this requires fancy quadrature (e.g., Lobatto)
	lgn_xpn_cff(wvl_idx,lgn_idx)+=phz_fnc_ffc(wvl_idx,ngl_idx)*static_cast<prc_cmp>(lgn_val[lgn_idx])*ngl_fct; // [frc] Legendre polynomial expansion coefficients ThS99 p. 177 (6.29)
      } // end loop over lgn
    } // end loop over wvl
  } // end loop over ngl
  for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
    if(lgn_nbr >= 1) asm_prm_lgn[wvl_idx]=lgn_xpn_cff(wvl_idx,1); // [frc] Asymmetry parameter estimated from Legendre expansion ThS99 p. 178
  } // end loop over wvl

  if(true){
    /* 20080113: (re-)Normalize Legendre expansion coefficients to zeroth moment
       Re-normalization attempts to represent the non-resolved phase function
       This is primarily (only?) the forward-scattered peak
       The re-normalization assumes that the un-resolved phase function, 
       determined by the zeroth moment (which should be 1.0 when fully resolved), 
       behaves the same as the resolved phase function.
       Hence all moments are scaled by the same amount
       This properly normalizes the resolved phase function to unity
       The error in the phase function approximation is still detectable by 
       the difference between the first moment and the asymmetry parameter */

    /* Wiscombe's notes on integrating over size distributions
       From MIEV.doc that accompanies Wis79/Wis96 distribution
       Further discussed in Van De Hulst (1982)

       INTEGRATING OVER SIZES
       ----------------------
       The normalized phase function for a single size parameter is
       
       P(one size) = 4 / ( XX**2 * QSCA ) * ( i1 + i2 ) / 2
       
       where  i1 + i2 = CABS(S1)**2 + CABS(S2)**2.  But it is
       ( i1 + i2 ), not  P(one size), that must be integrated
       over sizes when a size distribution is involved.
       (Physically, this means that intensities are added,
       not probabilities. )  An a posteriori normalization
       then gives the correct size-averaged phase function.
       
       Similarly, it is the CROSS-SECTIONS, proportional to
       XX**2  times QEXT,QSCA,QPR, which should be integrated
       over sizes, not QEXT,QSCA,QPR themselves.
       
       Similar remarks apply to PMOM.   The normalized moments are
       4 / ( XX**2 * QSCA ) * PMOM,  but it is PMOM itself, not
       these normalized moments, which should be integrated over
       a size distribution. */
    
    prc_cmp lgn_xpn_cff_zro; // [frc] Zeroth moment
    for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
      // For each order of Legendre polynomial for this wavelength...
      lgn_xpn_cff_zro=lgn_xpn_cff(wvl_idx,0);
      lgn_xpn_cff(wvl_idx,0)=1.0;
      for(lgn_idx=1;lgn_idx<=lgn_nbr;lgn_idx++){ // NB: loop starts from one
	lgn_xpn_cff(wvl_idx,lgn_idx)*=1.0/lgn_xpn_cff_zro;
      } // end loop over lgn
    } // end loop over wvl
  } // !true
  
  // Reconstruct phase function using Legendre expansion coefficients
  // fxm: Use fancier quadrature, e.g., trapezoidal or Simpson
  phz_fnc_lgn=0.0; // [sr-1] phase function Legendre expansion
  for(ngl_idx=0;ngl_idx<2*ngl_nbr-1;ngl_idx++){
    // Evaluate all orders of Legendre polynomial for this angle
    rcd+=gsl_sf_legendre_Pl_array(lgn_nbr,ngl_cos[ngl_idx],&lgn_val[0]); // [fnc] Legendre function
    for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
      for(lgn_idx=0;lgn_idx<=lgn_nbr;lgn_idx++){
	phz_fnc_lgn(wvl_idx,ngl_idx)+=(2*lgn_idx+1)*lgn_xpn_cff(wvl_idx,lgn_idx)*static_cast<prc_cmp>(lgn_val[lgn_idx]); // [sr-1] phase function Legendre expansion ThS99 p. 177 (6.28)
      } // end loop over lgn
    } // end loop over wvl
  } // end loop over ngl
  
  // Given effective phase function, compute integrated diagnostics that describe distribution of energy in forward-scattered direction
  
  // Compute partial sums
  std::valarray<prc_cmp> phz_fnc_bsf(2*ngl_nbr-1); // [frc] Phase function backward scattered fraction (fraction of energy scattered behind given angle, phz_fnc_bsf = 1.0 at 0 degrees)
  std::valarray<prc_cmp> phz_fnc_fsf(2*ngl_nbr-1); // [frc] Phase function forward scattered fraction (fraction of energy scattered forward of given angle, phz_fnc_fsf = 1.0 at 180 degrees)
  prc_cmp phz_fnc_nrg_frc; // [frc] Phase function fractional energy in current angular bin
  prc_cmp phz_fnc_lgn_nrg_frc; // [frc] Phase function Legendre expansion fractional energy in current angular bin
  // Initialize phase function integration
  std::valarray<prc_cmp> phz_fnc_ntg(0.0,wvl_nbr); // [frc] Spherically integrated phase function
  std::valarray<prc_cmp> phz_fnc_lgn_ntg(0.0,wvl_nbr); // [frc] Spherically integrated phase function Legendre expansion
  for(ngl_idx=0;ngl_idx<2*ngl_nbr-1;ngl_idx++){
    /* Factor of one-half multiplies phase function integral
       1/2 = 2*pi/(4*pi) = azimuthally symmetric phi integral / phase function normalization */
    ngl_fct=0.5*ngl_sin[ngl_idx]*ngl_wgt[ngl_idx]; // [frc] Factor common to given angle
    for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){

      phz_fnc_nrg_frc=phz_fnc_ffc(wvl_idx,ngl_idx)*ngl_fct;
      phz_fnc_lgn_nrg_frc=phz_fnc_lgn(wvl_idx,ngl_idx)*ngl_fct;
      
      if(wvl_idx == wvl_idx_dbg){
	if(ngl_idx == 0) phz_fnc_fsf[ngl_idx]=phz_fnc_nrg_frc; else phz_fnc_fsf[ngl_idx]=phz_fnc_fsf[ngl_idx-1]+phz_fnc_nrg_frc;
	phz_fnc_bsf[ngl_idx]=1.0-phz_fnc_fsf[ngl_idx];
      } // endif dbg
      
      // Zen01b eqn:phz_fnc_nrm_plr
      phz_fnc_ntg[wvl_idx]+=phz_fnc_nrg_frc; // [frc] Spherically integrated phase function
      phz_fnc_lgn_ntg[wvl_idx]+=phz_fnc_lgn_nrg_frc; // [frc] Spherically integrated phase function Legendre expansion
    } // end loop over wvl
  } // end loop over ngl
  
  // Normalize backward and forward scattered fractions so results are exact at 0 and 180 degrees, respectively
  for(ngl_idx=0;ngl_idx<2*ngl_nbr-1;ngl_idx++){
    phz_fnc_bsf[ngl_idx]/=phz_fnc_ntg[wvl_idx_dbg];
    phz_fnc_fsf[ngl_idx]/=phz_fnc_ntg[wvl_idx_dbg];
  } // end loop over ngl

  if(dbg_lvl_get() > dbg_fl){
    // Sanity check: Correctly normalized phase function should integrate to unity
    bool apx_eql; // [flg] Arguments are indistinguishable
    for(wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
      apx_eql=apx_eql_chk // [fnc] Determine whether arguments are indistinguishable
	(sbr_nm, // I [sng] Subroutine name of calling routine
	 true, // I [flg] Verbose output
	 PRC_CMP(1.0), // I [frc] Target argument
	 phz_fnc_ntg[wvl_idx], // I [frc] Approximation to target argument
	 PRC_CMP(1.0e-2), // I [frc] Relative precision (currently crude)
	 "Vetting phase function normalization"); // I [sng] Descriptive message of context
    } // end loop over wvl
  } // endif dbg
  
  if(dbg_lvl_get() == dbg_crr){
    (void)std::fprintf(stdout,"idx\tngl\tngl_dgr\tphz_fnc_ffc\tphz_fnc_lgn\tphz_fnc_bsf\tphz_fnc_fsf\n");
    for(ngl_idx=0;ngl_idx<2*ngl_nbr-1;ngl_idx++)
      (void)std::fprintf
	(stdout,"%4ld\t%.4f\t%.4f\t%.1e\t%.4f\t%.4f\t%.4f\n",
	 ngl_idx,ngl[ngl_idx],ngl_dgr[ngl_idx],phz_fnc_ffc(wvl_idx_dbg,ngl_idx),phz_fnc_lgn(wvl_idx_dbg,ngl_idx),phz_fnc_bsf[ngl_idx],phz_fnc_fsf[ngl_idx]);
  } // endif dbg
  
  std::cout << "Phase Function Diagnostics:" << std::endl;
  std::cout << "  Phase function: Archived at " << 2*ngl_nbr-1 << " angles" << std::endl;
  std::cout << "  Diagnostic angle is bin " << ngl_idx_dbg << " from " << ngl_min_dgr[ngl_idx_dbg] << "--" << ngl_max_dgr[ngl_idx_dbg] << " degrees" << std::endl;
  std::cout << "  Scattering fraction resolved by phase function = " << phz_fnc_ntg[wvl_idx_dbg] << " (must be very near 1.0 for accurate phase function expansion)" << std::endl;
  std::cout << "  Scattering fraction resolved by phase function Legendre expansion = " << phz_fnc_lgn_ntg[wvl_idx_dbg] << std::endl;
  std::cout << "  Asymmetry parameter estimated from Legendre expansion: " << asm_prm_lgn[wvl_idx_dbg] << std::endl;
  std::cout << "  " << lgn_nbr+1 << " Legendre expansion coefficients of order [0.." << lgn_nbr << "] = ";
  for(lgn_idx=0;lgn_idx<=lgn_nbr;lgn_idx++) std::cout << lgn_xpn_cff(wvl_idx_dbg,lgn_idx) << (lgn_idx != lgn_nbr ? ", " : "");
  std::cout << std::endl;

  // Populate diagnostic arrays at diagnostic angle
  const double ngl_cos_dbg=std::cos(ngl[ngl_idx_dbg]); // [frc] Cosine of scattering angle
  for(lgn_idx=0;lgn_idx<=lgn_nbr;lgn_idx++){
    rcd+=gsl_sf_legendre_Pl_array(lgn_nbr,ngl_cos_dbg,&lgn_val[0]); // [fnc] Legendre function
    lgn[lgn_idx]=lgn_idx; // [idx] Order of Legendre expansion term
  } // end loop over lgn

  // NB: number of returned angles is 2*ngl_nbr-1
  if(rcd != 0) err_prn(sbr_nm," rcd != 0 after Legendre computations in phz_fnc_mdl() (perhaps gsl_sf_legendre_Pl_array() failed?)");
  rcd=nco_redef(nc_out,NC_EINDEFINE); // [fnc] Put open netCDF dataset into define mode
  const int wvl_dmn(nco_inq_dimid(nc_out,static_cast<std::string>("wvl"))); // [dmn] Wavelength dimension
  const int ngl_dmn=nco_def_dim(nc_out,static_cast<std::string>("ngl"),2*ngl_nbr-1); // [dmn] Angle dimension
  const int lgn_dmn=nco_def_dim(nc_out,static_cast<std::string>("lgn"),lgn_nbr+1); // [dmn] Legendre dimension
  const int *dmn_scl((int *)NULL); // [dmn] Dummy dimension for scalars, should not be accessed
  // Derived dimensions
  const int *dmn_lgn(&lgn_dmn); // [dmn] Pointer to Legendre dimension
  const int *dmn_ngl(&ngl_dmn); // [dmn] Pointer to angle dimension
  const int *dmn_wvl(&wvl_dmn); // [dmn] Pointer to wavelength dimension
  const int dmn_wvl_lgn[2]={wvl_dmn,lgn_dmn}; // [dmn] Dimension ID array: wavelength, Legendre
  const int dmn_wvl_ngl[2]={wvl_dmn,ngl_dmn}; // [dmn] Dimension ID array: wavelength, angle
  const nc_type nco_xtyp(nco_get_xtype(static_cast<prc_cmp>(1.0))); // [enm] External netCDF type

  var_mtd_sct var_mtd[]={
    {0,"asm_prm_lgn",NC_FLOAT,1,dmn_wvl,"long_name","Asymmetry parameter estimated from Legendre expansion","units","fraction"},
    {0,"lgn_nbr",NC_FLOAT,0,dmn_scl,"long_name","Order of phase function Legendre expansion","units","number"},
    {0,"lgn",NC_INT,1,dmn_lgn,"long_name","Order of Legendre expansion term","units","index"},
    {0,"lgn_val",NC_FLOAT,1,dmn_lgn,"long_name","Legendre polynomial values","units","fraction"},
    {0,"lgn_xpn_cff",nco_xtyp,2,dmn_wvl_lgn,"long_name","Legendre polynomial expansion coefficients","units","fraction"},
    {0,"ngl",NC_FLOAT,1,dmn_ngl,"long_name","Scattering angle","units","radian"},
    {0,"ngl_min_dgr",NC_FLOAT,1,dmn_ngl,"long_name","Smallest scattering angle in bin","units","degree"},
    {0,"ngl_max_dgr",NC_FLOAT,1,dmn_ngl,"long_name","Largest scattering angle in bin","units","degree"},
    {0,"ngl_dgr",NC_FLOAT,1,dmn_ngl,"long_name","Scattering angle","units","degree"},
    {0,"ngl_dlt",NC_FLOAT,1,dmn_ngl,"long_name","Width of angle bin","units","radian"},
    {0,"ngl_wgt",NC_FLOAT,1,dmn_ngl,"long_name","Weight of angle bin","units","radian"},
    {0,"ngl_nbr",NC_FLOAT,0,dmn_scl,"long_name","Number of angle bins in phase function","units","number"},
    {0,"phz_fnc_bsf",NC_FLOAT,1,dmn_ngl,"long_name","Phase function backward scattered fraction (fraction of energy scattered behind given angle, phz_fnc_fsf = 1.0 at 0 degrees)","units","fraction"},
    {0,"phz_fnc_ffc",NC_FLOAT,2,dmn_wvl_ngl,"long_name","Effective phase function (weighted over size distribution)","units","steradian-1"},
    {0,"phz_fnc_fsf",NC_FLOAT,1,dmn_ngl,"long_name","Phase function forward scattered fraction (fraction of energy scattered forward of given angle, phz_fnc_fsf = 1.0 at 180 degrees)","units","fraction"},
    {0,"phz_fnc_lgn",NC_FLOAT,2,dmn_wvl_ngl,"long_name","phase function Legendre expansion","units","steradian-1"},
    {0,"phz_fnc_lgn_ntg",NC_FLOAT,1,dmn_wvl,"long_name","Spherically integrated phase function Legendre expansion","units","fraction"},
    {0,"phz_fnc_ntg",NC_FLOAT,1,dmn_wvl,"long_name","Spherically integrated phase function","units","fraction"},
    {0,"phz_fnc_dgn",NC_FLOAT,1,dmn_ngl,"long_name","Phase function at diagnostic size, wavelength","units","steradian-1"},
    {0,"plz_dgn",NC_FLOAT,1,dmn_ngl,"long_name","Degree of linear polarization at diagnostic size, wavelength","units","fraction"}
    //            {0,"foo",NC_FLOAT,1,dmn_ngl,"long_name","foo","units","fraction"},
  }; // end var_mtd_sct var_mtd[]
  const int var_mtd_nbr(sizeof(var_mtd)/sizeof(var_mtd_sct)); // [nbr] Number of variables in array
  rcd=nco_var_dfn(nc_out,var_mtd,var_mtd_nbr,dmn_nbr_max); // [fnc] Define variables in output netCDF file
  
  // Handle multi-dimensional arrays separately to conserve space if requested
  if(dmn_nbr_max >= 2){
    rcd=nco_put_var(nc_out,static_cast<std::string>("lgn_xpn_cff"),&lgn_xpn_cff(0,0));
    rcd=nco_put_var(nc_out,static_cast<std::string>("phz_fnc_ffc"),&phz_fnc_ffc(0,0));
    rcd=nco_put_var(nc_out,static_cast<std::string>("phz_fnc_lgn"),&phz_fnc_lgn(0,0));
  } // endif dmn_nbr_max >= 2
  
    // After writing, delete arrays 
  rcd=nco_put_var(nc_out,static_cast<std::string>("asm_prm_lgn"),&asm_prm_lgn[0]);
  rcd=nco_put_var(nc_out,static_cast<std::string>("lgn"),&lgn[0]);
  rcd=nco_put_var(nc_out,static_cast<std::string>("lgn_nbr"),lgn_nbr);
  rcd=nco_put_var(nc_out,static_cast<std::string>("lgn_val"),&lgn_val[0]);
  rcd=nco_put_var(nc_out,static_cast<std::string>("ngl"),ngl); delete []ngl;
  rcd=nco_put_var(nc_out,static_cast<std::string>("ngl_dgr"),&ngl_dgr[0]); delete []ngl_dgr;
  rcd=nco_put_var(nc_out,static_cast<std::string>("ngl_dlt"),ngl_dlt); delete []ngl_dlt;
  rcd=nco_put_var(nc_out,static_cast<std::string>("ngl_wgt"),ngl_wgt); delete []ngl_wgt;
  rcd=nco_put_var(nc_out,static_cast<std::string>("ngl_nbr"),ngl_nbr);
  rcd=nco_put_var(nc_out,static_cast<std::string>("phz_fnc_bsf"),&phz_fnc_bsf[0]);
  rcd=nco_put_var(nc_out,static_cast<std::string>("phz_fnc_fsf"),&phz_fnc_fsf[0]);
  rcd=nco_put_var(nc_out,static_cast<std::string>("phz_fnc_lgn_ntg"),&phz_fnc_lgn_ntg[0]);
  rcd=nco_put_var(nc_out,static_cast<std::string>("phz_fnc_ntg"),&phz_fnc_ntg[0]);
  rcd=nco_put_var(nc_out,static_cast<std::string>("phz_fnc_dgn"),phz_fnc_dgn); delete []phz_fnc_dgn;
  rcd=nco_put_var(nc_out,static_cast<std::string>("plz_dgn"),plz_dgn); delete []plz_dgn;
  
  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Exiting...");
  
  return rcd; // O [enm] Return success code
} // end phz_fnc_mdl()

// Global functions with C++ linkages end
