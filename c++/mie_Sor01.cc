// $Id$

// Purpose: Implementation of expressions for absorption and scattering cross sections of electromagnetic radiation by fractal aggregates.

/* Copyright (C) 2005--present Charlie Zender and Jorge Talamantes
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

/* History:
   200506: Jorge Talamantes implemented equations 30 and 33 in Sorensen (2001).
   
   Reference articles: Sor01
   
   Sorensen, Christopher, M., Light Scattering by Fractal Aggregates: A Review, Aerosol Science and Technology, 35 (2), pp. 648-687, 2001. */

/* Comments and equation numbers below refer to the Sorensen Review */

#include <mie_Sor01.hh> // Mie scattering solutions from MaS99

prc_cmp // O [m2] Aggregate cross-section for absorption (sigma_abs^agg)
xsx_abs_grg_Sor01 // [fnc] Compute aggregate cross-section for absorption (sigma_abs^agg)
(int mnm_nbr, // [nbr] Number of monomers in aggregate (N)
 prc_cmp rds_mnm, // [m] Monomer radius (a)
 prc_cmp wvn, // [m-1] Radiation wavenumber (k)
 std::complex<prc_cmp> idx_rfr) // [frc] Refractive index (m)
{
  const std::complex<prc_cmp> m2(idx_rfr*idx_rfr);
  prc_cmp E=imag((m2-1.0)/(m2+2.0)); // Sor01 p. 654 (31)
  
  prc_cmp xsx_abs_mnm=4.0*mth::cst_M_PIl*wvn*pow(rds_mnm,3)*E; // Sor01 p. 654 (30)
  return mnm_nbr*xsx_abs_mnm; // Sor01 p. 654 (29)
} // end xsx_abs_grg_Sor01()

prc_cmp // [m2] Aggregate cross-section for scattering (sigma_sca^agg)
xsx_sct_grg_Sor01 // [fnc] Compute aggregate cross-section for scattering (sigma_sca^agg)
(int mnm_nbr, // [nbr] Number of monomers in aggregate (N)
 prc_cmp dmn_frc, // [frc] Fractal dimension of aggregate (D)
 prc_cmp rds_gyr, // [m] Radius of gyration of aggregate (R_g)
 prc_cmp rds_mnm, // [m] Monomer radius (a)
 prc_cmp wvn, // [m-1] Radiation wavenumber (k)
 std::complex<prc_cmp> idx_rfr) // [frc] Refractive index (m)
{
  const std::complex<prc_cmp> m2(idx_rfr*idx_rfr);
  prc_cmp F=std::abs((m2-1.0)/(m2+2.0)); // Sor01 p. 654 (27)
  F*=F; // Sor01 p. (27)
  
  prc_cmp G=1.0+4.0*wvn*wvn*rds_gyr*rds_gyr/(3.0*dmn_frc);
  G=std::pow(G,-(dmn_frc/2.0)); // Sor01 p. 655 (35)
  
  prc_cmp xsx_sct_mnm=(8.0*mth::cst_M_PIl/3.0)*std::pow(wvn,4)*std::pow(rds_mnm,6)*F; // monomer scattering cross-section Sor01 p. 655 (34)

  return mnm_nbr*mnm_nbr*xsx_sct_mnm*G; // Sor01 p. 654 (33)
} // end xsx_sct_grg_Sor01()

int // O [enm] Return success code
mie_Sor01_tst // [fnc] Test function for Sor01 physics
(void)
{
  // Purpose: Test function for Sor01 physics

  int rcd(0); // [enm] Return success code

  prc_cmp lambda=1.0; // wavelength of light in microns
  
  prc_cmp wvn=2.0*mth::cst_M_PIl/lambda; // wavenumber
  prc_cmp rds_mnm=0.016; // monomer radius in microns
  prc_cmp dmn_frc=1.8; // [frc] Fractal dimension of aggregate
  prc_cmp rds_gyr=1.0/wvn; // radius of gyration of aggregate
  std::complex<prc_cmp> idx_rfr(1.7,0.7); // index of refraction of soot
  
  for(int mnm_nbr=1;mnm_nbr<=200;mnm_nbr++)
    std::cout << "N = " << mnm_nbr << std::endl
	      << "sigma abs agg = " << xsx_abs_grg_Sor01(mnm_nbr,rds_mnm,wvn,idx_rfr) << std::endl
	      << "sigma sca agg = " << xsx_sct_grg_Sor01(mnm_nbr,dmn_frc,rds_gyr,rds_mnm,wvn,idx_rfr) << std::endl
	      << std::endl;
  
  return rcd;
} // end mie_Sor01_tst()
