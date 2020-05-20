// $Id$ 

// Purpose: Implementation (declaration) of line-by-line utilities

/* Copyright (C) 1997--present Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

#include <lbl.hh> // Line-by-line radiative transfer utilites

// line-by-line classes

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

int // O [enm] Return success code
lnshp_lrn_vct // [fnc] Lorentzian line shape profile
(const long frq_nbr, // I [nbr] Size of arrays
 const double &HWHM_lrn, // I [Hz,cm-1] Lorentz half-width at half-maximum
 const double *frq_dlt, // I [Hz,cm-1] Distance from line center
 double *lnshp_lrn_val) // O [Hz-1,cm] Lorentzian line shape profile
{
  /* Purpose: Return Lorentzian line shape profile for array of frequencies
     Single scalar computations should use lnshp_lrn_fst_scl() */
  int rcd(0); // Return success code
  // Local
  using mth::cst_M_PIl; // (3.1415926535897932384626433832795029L) [frc] 3
  long idx; // [idx] Counting index
  const std::string sbr_nm("lnshp_lrn_vct"); // [sng] Name of subroutine
  const double lnshp_fct(HWHM_lrn/mth::cst_M_PIl); // [Hz,cm-1] Line shape factor
  const double HWHM_lrn_sqr(HWHM_lrn*HWHM_lrn); // [Hz2,cm-2] Line shape factor
  // Main code
  for(idx=0;idx<frq_nbr;idx++){
    lnshp_lrn_val[idx]=lnshp_fct/(frq_dlt[idx]*frq_dlt[idx]+HWHM_lrn_sqr); // O [Hz-1,cm] Lorentzian line shape profile
  } // end loop over frq
  return rcd; // [enm] Return success code
} // end lnshp_lrn_vct()

int // O [enm] Return success code
lnshp_dpp_vct // [fnc] Doppler line shape profile
(const long frq_nbr, // I [nbr] Size of arrays
 const double &HWHM_dpp, // I [Hz,cm-1] Doppler width (half-width at (1/e)-maximum)
 const double *frq_dlt, // I [Hz,cm-1] Distance from line center
 double *lnshp_dpp_val) // O [Hz-1,cm] Doppler line shape profile
{
  /* Purpose: Return Doppler line shape profile for array of frequencies
     Single scalar computations should use lnshp_lrn_fst_scl() */
  int rcd(0); // Return success code
  // Local
  long idx; // [idx] Counting index
  const std::string sbr_nm("lnshp_dpp_vct"); // [sng] Name of subroutine
  double xpn_fct; // [frc] Factor in exponent
  using mth::cst_M_PIl; // (3.1415926535897932384626433832795029L) [frc] 3
  const double pre_xpn_fct(1.0/(HWHM_dpp*std::sqrt(mth::cst_M_PIl))); // [Hz-1,cm] Factor preceding exponential
  // Main code
  for(idx=0;idx<frq_nbr;idx++){
    xpn_fct=frq_dlt[idx]/HWHM_dpp; // [frc] Factor preceding exponential
    lnshp_dpp_val[idx]=pre_xpn_fct*std::exp(-xpn_fct*xpn_fct); // O [Hz-1,cm] Doppler line shape profile
  } // end loop over frq
  return rcd; // [enm] Return success code
} // end lnshp_dpp_vct()

int // O [enm] Return success code
lnshp_vgt_vct // [fnc] Doppler line shape profile
(const long frq_nbr, // I [nbr] Size of arrays
 const double &HWEM_dpp, // I [Hz,cm-1] Doppler width (half-width at (1/e)-maximum)
 const double &HWHM_lrn, // I [Hz,cm-1] Lorentz half-width at half-maximum
 const double *frq_dlt, // I [Hz,cm-1] Distance from line center
 double *lnshp_vgt_val) // O [Hz-1,cm] Doppler line shape profile
{
  // Purpose: Return Voigt line shape profile
  int rcd(0); // Return success code
  // Local
  long idx; // [idx] Counting index
  const std::string sbr_nm("lnshp_vgt_vct"); // [sng] Name of subroutine
  /* Voigt half-width parameter aka damping ratio 
   ThS99 p. 71 (3.17), GoY89 p. 112 (3.83), Lio92 p. 31 (2.2.10)
   GoY89 note that most tabulations are in terms of 2*HWHM_lrn/HWHM_dpp not 1*HWHM_lrn/HWHM_dpp */

  // Main code
  // Order is important!
  // Convert Doppler half-width at 1/e maximum to HWHM
  using mth::cst_M_PIl; // (3.1415926535897932384626433832795029L) [frc] 3
  using mth::cst_M_LN2l; // (0.6931471805599453094172321214581766L) [frc] log_e(2.0)
  const double sqrt_ln_2(std::sqrt(cst_M_LN2l)); // [Hz,cm-1] Doppler half-width at half-maximum
  const double HWHM_dpp(sqrt_ln_2*HWEM_dpp); // [Hz,cm-1] Doppler half-width at half-maximum

  // Compute frequency-invariant factors
  const double vgt_pre_fct(1.0/HWHM_dpp); // [Hz-1,cm] Voigt function pre-factor
  const double arg_rl_fct(sqrt_ln_2/HWHM_dpp); // [Hz-1,cm] Factor in real part of complex error function argument
  const double err_fnc_cpx_arg_img(sqrt_ln_2*HWHM_lrn/HWHM_dpp); // [frc] Complex error function imaginary argument

  // Frequency varying factors 
  double err_fnc_cpx_arg_rl; // [frc] Complex error function real argument
  std::complex<double> err_fnc_cpx_arg; // [frc] Complex error function argument
  std::complex<double> err_fnc_cpx_val; // [frc] Complex error function value
  for(idx=0;idx<frq_nbr;idx++){
    err_fnc_cpx_arg_rl=arg_rl_fct*frq_dlt[idx]; // [frc] Complex error function real argument
    err_fnc_cpx_arg=std::complex<double>(err_fnc_cpx_arg_rl,err_fnc_cpx_arg_img); // [frc] Complex error function argument
    err_fnc_cpx_val=err_fnc_cpx_Hum82(err_fnc_cpx_arg); // [frc] Complex error function
    lnshp_vgt_val[idx]=vgt_pre_fct*err_fnc_cpx_val.real(); // O [Hz-1,cm] Voigt line shape profile
  } // end loop over frq

  // Compute frequency-invariant factors
  const double HWHM_vgt= // [Hz,cm-1] "Half-width" of Voigt profile Lio92 p. 31 (2.2.17)
    0.5*(HWHM_lrn+std::sqrt(HWHM_lrn*HWHM_lrn+4.0*cst_M_LN2l*HWEM_dpp*HWEM_dpp))+
    0.05*HWHM_lrn*(1.0-2.0*HWHM_lrn/std::sqrt(HWHM_lrn*HWHM_lrn+4.0*cst_M_LN2l*HWEM_dpp*HWEM_dpp));

  // Frequency varying factors 
  double eta; // [frc] 
  double xi; // [frc] 
  for(idx=0;idx<frq_nbr;idx++){
    eta=frq_dlt[idx]/HWHM_vgt; // [frc]
    xi=HWHM_lrn/HWHM_vgt; // [frc]
    lnshp_vgt_val[idx]= // O [Hz-1,cm] Voigt line shape profile Lio92 p. 31 (2.2.16)
      std::sqrt(cst_M_LN2l/mth::cst_M_PIl)*(1.0/HWHM_vgt)*(1.0-xi)*std::exp(-cst_M_LN2l*eta*eta)+
      (1.0/(mth::cst_M_PIl*HWHM_vgt))*xi*(1.0/(1.0+eta*eta))-
      (1.0/(mth::cst_M_PIl*HWHM_vgt)*xi*(1.0-xi)*(xi+1.0+1.5/cst_M_LN2l))*
      (0.066*std::exp(-0.4*eta*eta)-1.0/(40.0-5.5*eta*eta+eta*eta*eta*eta)); 
  } // end loop over frq

  return rcd; // [enm] Return success code
} // end lnshp_vgt_vct()

int // O [enm] Return success code
rt_lbl // [fnc] Single line radiative transfer
(const int &nc_out, // I [fl] netCDF file for output 
 const int &dmn_nbr_max, // I [nbr] Maximum number of dimensions allowed in single variable in output file
 const prc_cmp &prs_dlt, // I [Pa] Pressure thickness
 const prc_cmp &prs_mdp, // I [Pa] Midlayer pressure
 const prc_cmp &slr_zen_ngl_cos, // I [frc] Cosine solar zenith angle
 const prc_cmp &tpt_mdp, // I [K] Midlayer temperature
 const prc_cmp &vmr_gas, // I [mlc mlc-1] Volume mixing ratio of gas
 const std::string &lbl_tst, // I [sng] Name of line-by-line test
 const wvl_grd_cls &wvlgrd) // I [m] Wavelength grid
{
  /* Purpose: Radiative transfer of single transition */

  /* HITRAN test:
     mie --dbg=1 --tst=lbl
     mie --dbg=1 --no_wrn_ntp --wvn_mdp=936.804 --wvn_dlt=0.2 --wvn_nbr=11 --tst=lbl --lbl_tst=CO2 --prs_mdp=10000.0
     ncks -F -C -v mmw_mlc,mmw_iso ${DATA}/mie/mie.nc
     ncks -F -C -v trn_clm,lnshp_dpp_wvn,lnshp_prs_wvn,lnshp_vgt_wvn,abs_cff_mss_gas ${DATA}/mie/mie.nc
   */

  int rcd(0); // Return success code
  // Local
  long idx; // [idx] Counting index
  long mlc_idx; // [idx] Counting index for molecule
  long iso_idx; // [idx] Counting index for isotopomer
  const std::string sbr_nm("rt_lbl"); // [sng] Name of subroutine
  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Entering...");
  
  // Open HITRAN line database file
  std::string fl_htrn("/data/zender/hitran/CO2.nc"); // [sng] HITRAN File
  const int nc_in(nco_open(fl_htrn,NC_NOWRITE)); // [fnc] Open netCDF file
  // Input line statistics
  const long ln_nbr(nco_inq_dimlen(nc_in,static_cast<std::string>("ln_ctr"))); // [nbr] Number of lines in input file
  // Close HITRAN line database file
  rcd=nco_close(nc_in);
  std::cerr << prg_nm_get() << ": INFO Ingested " << fl_htrn << std::endl;
  if(dbg_lvl_get() > dbg_off) std::cerr << "Number of lines in "+fl_htrn+" is "<< ln_nbr << std::endl;

  /* HITRAN line statistics
     Suffix _htrn indicates variable is HITRAN-tabulated value (CGS units, 296 K), not SI value
     _htrn variables are */
  double HWHM_air_htrn(CEWI_dbl); // [cm-1 atm-1] Air-broadened HWHM, HITRAN
  double HWHM_slf_htrn(CEWI_dbl); // [cm-1 atm-1] Self-broadened HWHM, HITRAN
  double HWHM_tpt_dpn_xpn(CEWI_dbl); // [frc] Exponent of temperature dependence of air-broadened HWHM
  double ln_ctr_htrn(CEWI_dbl); // [cm-1] Wavenumber at line center, HITRAN
  double ln_nrg_lwr_htrn(CEWI_dbl); // [cm-1] Lower state energy of transition, HITRAN
  double ln_sft_air_htrn(CEWI_dbl); // [cm-1 atm-1] Air-broadened pressure shift of line transition, HITRAN
  double ln_str_htrn(CEWI_dbl); // [cm-1 cm2 mlc-1] Line strength, HITRAN
  short ln_iso(CEWI_sht); // [enm] HITRAN isotope number (1..9)
  short ln_mlc(CEWI_sht); // [enm] HITRAN molecule number (1..37)
  short iso_id(CEWI_sht); // [enm] HITRAN Isotopomer ID
  short mlc_id(CEWI_sht); // [enm] HITRAN Molecule ID
  short qnt_glb_lwr(CEWI_sht); // [enm] Lower state global quanta index (vibrational quanta)
  short qnt_glb_upr(CEWI_sht); // [enm] Upper state global quanta index (vibrational quanta)
  std::string qnt_lcl_lwr; // [sng] Lower state local quanta index (rotational quanta)
  std::string qnt_lcl_upr; // [sng] Upper state local quanta index (rotational quanta)

  /* Notes about specific line-by-line tests:
     Documented examples of specific lines in literature are rare
     Lio92 p. 32 Figure 2.3 gives nice H2O, O3, and CO2 examples */
  if(lbl_tst=="CO2"){
    // 936.8 cm-1 CO2 transition pictured in ThS99 p. 92 Fgr 4.5
    // ncks -F -d ln_ctr,936.8 ${DATA}/hitran/CO2.nc
    HWHM_air_htrn=0.069; // [cm-1 atm-1] Air-broadened HWHM, HITRAN
    HWHM_slf_htrn=0.0897; // [cm-1 atm-1] Self-broadened HWHM, HITRAN
    HWHM_tpt_dpn_xpn=0.78; // [frc] Exponent of temperature dependence of air-broadened HWHM
    ln_ctr_htrn=936.804; // [cm-1] Wavenumber at line center, HITRAN
    ln_nrg_lwr_htrn=1704.94; // [cm-1] Lower state energy of transition, HITRAN
    ln_sft_air_htrn=0.0; // [cm-1 atm-1] Air-broadened pressure shift of line transition, HITRAN
    ln_str_htrn=1.382e-23; // [cm-1 cm2 mlc-1] Line strength, HITRAN
    ln_iso=1; // [enm] HITRAN isotope number (1..9)
    ln_mlc=2; // [enm] HITRAN molecule number (1..37)
    iso_id=1; // [enm] HITRAN Isotopomer ID
    mlc_id=2; // [enm] HITRAN Molecule ID
    qnt_glb_lwr=5; // [enm] Lower state global quanta index (vibrational quanta)
    qnt_glb_upr=9; // [enm] Upper state global quanta index (vibrational quanta)
    qnt_lcl_lwr="    P 28 "; // [sng] Lower state local quanta index (rotational quanta)
    qnt_lcl_upr="         "; // [sng] Upper state local quanta index (rotational quanta)
  }else if(lbl_tst=="H2O"){
    // H2O transition nearest to 1.0 um
    // ncks -F -d ln_ctr,10000.0 ${DATA}/hitran/H2O.nc
    HWHM_air_htrn=0.044; // [cm-1 atm-1] Air-broadened HWHM, HITRAN
    HWHM_slf_htrn=0.285; // [cm-1 atm-1] Self-broadened HWHM, HITRAN
    HWHM_tpt_dpn_xpn=0.3; // [frc] Exponent of temperature dependence of air-broadened HWHM
    ln_ctr_htrn=9999.49; // [cm-1] Wavenumber at line center, HITRAN
    ln_nrg_lwr_htrn=744.163; // [cm-1] Lower state energy of transition, HITRAN
    ln_sft_air_htrn=-0.0253; // [cm-1 atm-1] Air-broadened pressure shift of line transition, HITRAN
    ln_str_htrn=1.471e-25; // [cm-1 cm2 mlc-1] Line strength, HITRAN
    ln_iso=1; // [enm] HITRAN isotope number (1..9)
    ln_mlc=1; // [enm] HITRAN molecule number (1..37)
    iso_id=1; // [enm] HITRAN Isotopomer ID (1..90)
    mlc_id=1; // [enm] HITRAN Molecule ID (1..37)
    qnt_glb_lwr=1; // [enm] Lower state global quanta index (vibrational quanta)
    qnt_glb_upr=20; // [enm] Upper state global quanta index (vibrational quanta)
    qnt_lcl_lwr=" 8 1 8   "; // [sng] Lower state local quanta index (rotational quanta)
    qnt_lcl_upr=" 9 1 9   "; // [sng] Upper state local quanta index (rotational quanta)
  }else{
    err_prn(sbr_nm,"Unknown line to test");
  } // endif

  if(false){
    std::cout << "qnt_glb_lwr = " << qnt_glb_lwr << std::endl; // CEWI
    std::cout << "qnt_glb_upr = " << qnt_glb_upr << std::endl; // CEWI
  } // endif false

    // Load database of names and mean molecular weights of all molecules and isotopomers
  using htrn::mlc_nbr_max_htrn; // (38) [nbr] Number of molecules in HITRAN database
  using htrn::iso_nbr_max_htrn; // (90) [nbr] Number of isotopomers in HITRAN database
  using htrn::iso_per_mlc_nbr_max_htrn; // (8) [nbr] Maximum number of isotopes of a molecule in HITRAN database
  std::string *mlc_sng=new std::string[mlc_nbr_max_htrn+1]; // [sng] HITRAN molecule names
  std::string *iso_sng=new std::string[iso_nbr_max_htrn+1]; // [sng] HITRAN isotopomer names
  double *mmw_mlc=new double[mlc_nbr_max_htrn+1]; // [kg mol-1] Mean molecular weight of molecule
  double *mmw_iso=new double[iso_nbr_max_htrn+1]; // [kg mol-1] Mean molecular weight of isotopomer
  double *xpn_mlc=new double[mlc_nbr_max_htrn+1]; // [frc] Exponent in temperature dependence of partition function
  short iso_idx_map[mlc_nbr_max_htrn+1][iso_per_mlc_nbr_max_htrn+1]; // [frc] Exponent in temperature dependence of partition function
  rcd+=mmw_mlc_get(mmw_mlc); // [kg mol-1] Mean molecular weight of molecule
  rcd+=mmw_iso_get(mmw_iso); // [kg mol-1] Mean molecular weight of isotopomer
  rcd+=mlc_sng_get(mlc_sng); // [sng] HITRAN molecule names
  rcd+=iso_sng_get(iso_sng); // [sng] HITRAN isotopomer names
  rcd+=rtl_fnc_tpt_xpn_get(xpn_mlc); // [sng] HITRAN isotopomer names
  rcd+=iso_idx_map_get(iso_idx_map); // [map] Map [mlc_id,iso_id]->istpmr_id
  const short istpmr_idx(iso_idx_map[mlc_id][iso_id]); // [kg mol-1] Isotopomer index
  const double mmw_gas(mmw_iso[istpmr_idx]); // [kg mol-1] Mean molecular weight of gas
  const double xpn_gas(xpn_mlc[mlc_id]); // [frc] Exponent in temperature dependence of partition function

  if(true){
    std::cout << "Single line radiative transfer:" << std::endl;
    std::cout << "  ln_mlc = " << ln_mlc << " = "+mlc_sng[ln_mlc]+", ln_iso = " << ln_iso << " = "+iso_sng[istpmr_idx] << std::endl;
    std::cout << "  HITRAN input:" << std::endl;
    std::cout << "  HWHM_air_htrn = " << HWHM_air_htrn << " cm-1 atm-1" << std::endl;
    std::cout << "  HWHM_slf_htrn = " << HWHM_slf_htrn << " cm-1 atm-1" << std::endl;
    std::cout << "  ln_ctr_htrn = " << ln_ctr_htrn << " cm-1" << std::endl;
    std::cout << "  ln_nrg_lwr_htrn = " << ln_nrg_lwr_htrn << " cm-1" << std::endl;
    std::cout << "  ln_sft_air_htrn = " << ln_sft_air_htrn << " cm-1" << std::endl;
    std::cout << "  ln_str_htrn = " << ln_str_htrn << " cm-1 cm2 mlc-1" << std::endl;
  } // endif true

  /* Store HITRAN variables in wavenumber, frequency, and wavelength units
     Frequency-based SI equivalents of HITRAN variables have _frq appended
     Wavelength-based SI equivalents of HITRAN variables have _wvl appended
     Append _wvn to wavenumber-based variables scaled to ambient T and p conditions */

  using phc::cst_Avagadro; // (6.022045e+23) [mlc mol-1] Avagadro's number
  using phc::cst_Planck; // (6.62620e-34) [J s] Planck's constant
  using phc::grv_sfc_mean; // (9.80665) [m s-2] Mean gravitational acceleration at Earth's surface
  using phc::mmw_dry_air; // (28.9644e-3) [kg mol-1] (Source: radcsw.F in CCM2/3)
  using phc::prs_HITRAN; // (101325.0) [Pa] Reference pressure for HITRAN database
  using phc::speed_of_light; // (2.99793e+08) [m s-1] Speed of light in vacuo
 
  // Thermodynamic environment
  const prc_cmp mmr_gas(vmr_gas*mmw_gas/mmw_dry_air); // [kg kg-1] Mass mixing ratio of gas
  const prc_cmp mpc_gas(mmr_gas*prs_mdp/grv_sfc_mean); // [kg m-2] Mass path column of gas
  const prc_cmp mpl_gas(mmr_gas*prs_dlt/grv_sfc_mean); // [kg m-2] Mass path layer of gas

  // Air-broadened HWHM, HITRAN
  const double HWHM_air_htrn_wvn(HWHM_air_htrn); // [cm-1] Air-broadened HWHM, HITRAN
  const double HWHM_air_htrn_frq(HWHM_air_htrn*100.0*speed_of_light); // [cm-1]->[Hz] Air-broadened HWHM, HITRAN
  const double HWHM_air_htrn_wvl(0.01/(ln_ctr_htrn-HWHM_air_htrn)-0.01/(ln_ctr_htrn+HWHM_air_htrn)); // [cm-1]->[m] Air-broadened HWHM, HITRAN

  // Self-broadened HWHM, HITRAN
  const double HWHM_slf_htrn_wvn(HWHM_slf_htrn); // [cm-1] Self-broadened HWHM, HITRAN
  const double HWHM_slf_htrn_frq(HWHM_slf_htrn*100.0*speed_of_light); // [cm-1]->[Hz] Self-broadened HWHM, HITRAN
  const double HWHM_slf_htrn_wvl(0.01/(ln_ctr_htrn-HWHM_slf_htrn)-0.01/(ln_ctr_htrn+HWHM_slf_htrn)); // [cm-1]->[m] Self-broadened HWHM, HITRAN

  // Unshifted line center
  const double ln_ctr_htrn_wvn(ln_ctr_htrn); // [cm-1] Unshifted line center, HITRAN
  const double ln_ctr_htrn_frq(ln_ctr_htrn*100.0*speed_of_light); // [cm-1]->[Hz] Unshifted line center, HITRAN
  const double ln_ctr_htrn_wvl(0.01/ln_ctr_htrn); // [cm-1]->[m] Unshifted line center, HITRAN

  // Lower state energy of transition
  const double ln_nrg_lwr_wvn(ln_nrg_lwr_htrn); // [cm-1] Lower state energy of transition
  const double ln_nrg_lwr_frq(100.0*speed_of_light*ln_nrg_lwr_htrn); // [cm-1]->[Hz] Lower state energy of transition
  const double ln_nrg_lwr_wvl(0.01/ln_nrg_lwr_htrn); // [cm-1]->[m] Lower state energy of transition
  const double ln_nrg_lwr_J(100.0*speed_of_light*cst_Planck*ln_nrg_lwr_htrn); // [cm-1]->[J] Lower state energy of transition

  // Pressure shift
  const double prs_fct(prs_mdp/prs_HITRAN); // [frc] Fraction of atmosphere
  const double ln_sft_air_wvn(ln_sft_air_htrn*prs_fct); // [cm-1 atm-1]->[cm-1] Air-broadened pressure shift of line transition
  const double ln_sft_air_frq(ln_sft_air_htrn*prs_fct*100.0*speed_of_light); // [cm-1 atm-1]->[Hz] Air-broadened pressure shift of line transition
  const double ln_sft_air_wvl(ln_sft_air_htrn*prs_fct*0.01/ln_ctr_htrn); // [cm-1 atm-1]->[m] Air-broadened pressure shift of line transition

  // Line strength
  const double ln_str_htrn_wvn(ln_str_htrn*1.0e-4); // [cm-1 cm2 mlc-1]->[cm-1 m2 mlc-1] Line strength, HITRAN
  const double ln_str_htrn_frq(ln_str_htrn*1.0e-4*100.0*speed_of_light); // [cm-1 cm2 mlc-1]->[Hz m2 mlc-1] Line strength, HITRAN
  // fxm: ln_str_wvl may be bogus
  const double ln_str_htrn_wvl(ln_str_htrn*1.0e-4*0.01/(ln_ctr_htrn)); // [cm-1 cm2 mlc-1]->[m m2 mlc-1] Line strength, HITRAN

  // Pressure-shifted line center
  const double ln_ctr_wvn(ln_ctr_htrn_wvn+ln_sft_air_wvn); // [cm-1] Pressure-shifted line center
  const double ln_ctr_frq(ln_ctr_htrn_frq+ln_sft_air_frq); // [Hz] Pressure-shifted line center
  const double ln_ctr_wvl(ln_ctr_htrn_wvl+ln_sft_air_wvl); // [m] Pressure-shifted line center

  // Determine partial pressure of gas from volume mixing ratio, temperature, and pressure
  using phc::gas_cst_unv; // (8.31441) [J mol-1 K-1] Universal gas constant
  const prc_cmp gas_cst_gas(gas_cst_unv/mmw_gas); // [J kg-1 K-1] Gas constant of gas
  const prc_cmp cnc_air(prs_mdp/(gas_cst_unv*tpt_mdp)); // [# m-3] Number concentration of air
  const prc_cmp cnc_gas(vmr_gas*cnc_air); // [# m-3] Number concentration of gas
  const prc_cmp dns_gas(cnc_gas*mmw_gas); // [kg m-3] Mass density of gas
  const prc_cmp prs_gas(dns_gas*gas_cst_gas*tpt_mdp); // [Pa] Partial pressure of gas

  /* Scale tabulated HWHM to ambient p and T
     Pressure scaling: HWHM increases linearly with pressure because collision frequency
     depends linearly on number concentration of collision partners, i.e., air molecules
     Temperature scaling: HWHM decreases as T increases
     This is because although thermal speed increases as sqrt(T), 
     the mean free path between collisions increases as T.
     Net result of increasing T is decreased collision frequency, thus decreased HWHM */
  using phc::tpt_HITRAN; // (296.0) [K] Reference temperature for HITRAN database
  using phc::cst_Boltzmann; // (1.38063e-23) [J K-1] Boltzmann's constant
  using phc::joules_per_eV; // (1.60217733e-19) [J ev-1] Joules per electron volt CRC95 inside back cover
  const double tpt_fct(std::pow(tpt_HITRAN/tpt_mdp,HWHM_tpt_dpn_xpn)); // [frc] Temperature dependent HWHM factor

  // Foreign broadening
  const double HWHM_air_wvn(tpt_fct*HWHM_air_htrn_wvn*(prs_mdp-prs_gas)/prs_HITRAN); // [cm-1] Air-broadened HWHM
  const double HWHM_air_frq(tpt_fct*HWHM_air_htrn_frq*(prs_mdp-prs_gas)/prs_HITRAN); // [Hz] Air-broadened HWHM
  const double HWHM_air_wvl(tpt_fct*HWHM_air_htrn_wvl*(prs_mdp-prs_gas)/prs_HITRAN); // [m] Air-broadened HWHM

  // Self broadening
  const double HWHM_slf_wvn(tpt_fct*HWHM_slf_htrn_wvn*prs_gas/prs_HITRAN); // [cm-1] Self-broadened HWHM
  const double HWHM_slf_frq(tpt_fct*HWHM_slf_htrn_frq*prs_gas/prs_HITRAN); // [Hz] Self-broadened HWHM
  const double HWHM_slf_wvl(tpt_fct*HWHM_slf_htrn_wvl*prs_gas/prs_HITRAN); // [m] Self-broadened HWHM

  // Total HWHM is sum of foreign- and self-broadening
  const double HWHM_prs_wvn(HWHM_air_wvn+HWHM_slf_wvn); // [cm-1] Pressure-broadened HWHM
  const double HWHM_prs_frq(HWHM_air_frq+HWHM_slf_frq); // [Hz] Pressure-broadened HWHM
  const double HWHM_prs_wvl(HWHM_air_wvl+HWHM_slf_wvl); // [m] Pressure-broadened HWHM

  // Doppler width depends only on pressure
  // Factor of mmw_gas/cst_Avagadro is simply the molecular mass
  const double HWEM_dpp_wvn(ln_ctr_wvn*std::sqrt(2.0*cst_Boltzmann*cst_Avagadro*tpt_mdp/mmw_gas)/speed_of_light); // [cm-1] Doppler width (half-width at (1/e)-maximum)
  const double HWEM_dpp_frq(ln_ctr_frq*std::sqrt(2.0*cst_Boltzmann*cst_Avagadro*tpt_mdp/mmw_gas)/speed_of_light); // [Hz] Doppler width (half-width at (1/e)-maximum)
  const double HWEM_dpp_wvl(ln_ctr_wvl*std::sqrt(2.0*cst_Boltzmann*cst_Avagadro*tpt_mdp/mmw_gas)/speed_of_light); // [m] Doppler width (half-width at (1/e)-maximum)
  if(HWEM_dpp_wvn > HWHM_prs_wvn) err_prn(sbr_nm,"HWEM_dpp > HWHM_prs so approximation in inhomogeneous path to space breaks down");

  // Doppler half-width at half-maximum is useful diagnostic
  // Note: HWEM > HWHM since sqrt(log(2.0)) < 1.0
  using mth::cst_M_LN2l; // (0.6931471805599453094172321214581766L) [frc] log_e(2.0)
  const double HWHM_dpp_wvn(HWEM_dpp_wvn*std::sqrt(cst_M_LN2l)); // [cm-1] Doppler half-width at half-maximum
  const double HWHM_dpp_frq(HWEM_dpp_frq*std::sqrt(cst_M_LN2l)); // [Hz] Doppler half-width at half-maximum
  const double HWHM_dpp_wvl(HWEM_dpp_wvl*std::sqrt(cst_M_LN2l)); // [m] Doppler half-width at half-maximum

  // Voigt line profile parameter
  const double HWHM_vgt_wvn(HWHM_prs_wvn/HWHM_dpp_wvn); // [cm-1] Voigt line profile parameter
  const double HWHM_vgt_frq(HWHM_prs_frq/HWHM_dpp_frq); // [Hz] Voigt line profile parameter
  const double HWHM_vgt_wvl(HWHM_prs_wvl/HWHM_dpp_wvl); // [m] Voigt line profile parameter

  /* Scale tabulated line strength to T
     Line strength must be scaled by temperature-sensitivity of 
     1. Partition function
     2. Boltzmann factor for lower state energy
     3. Stimulated emission factor 
     Scaling is described in detail in RRG98 p. 710 and htrn2nb.F */
  double ln_str_scl_wvn; // [cm-1 m2 mlc-1] Line strength at ambient temperature
  double ln_str_scl_frq; // [Hz m2 mlc-1] Line strength at ambient temperature
  double ln_str_scl_wvl; // [m m2 mlc-1] Line strength at ambient temperature
  // Initialize to HITRAN value at 296 K
  ln_str_scl_wvn=ln_str_htrn_wvn; // [cm-1 m2 mlc-1] Line strength at ambient temperature
  ln_str_scl_frq=ln_str_htrn_frq; // [Hz m2 mlc-1] Line strength at ambient temperature
  ln_str_scl_wvl=ln_str_htrn_wvl; // [m m2 mlc-1] Line strength at ambient temperature

  // Scale line strength by temperature dependence of partition function
  const double prt_fnc_tpt_scl(std::pow(tpt_HITRAN/tpt_mdp,xpn_gas)); // [frc] ! Lio92 p. 33 (2.2.21),
  ln_str_scl_wvn*=prt_fnc_tpt_scl; // [cm-1 m2 mlc-1] Line strength at ambient temperature
  ln_str_scl_frq*=prt_fnc_tpt_scl; // [Hz m2 mlc-1] Line strength at ambient temperature
  ln_str_scl_wvl*=prt_fnc_tpt_scl; // [m m2 mlc-1] Line strength at ambient temperature

  // Scale line strength by temperature dependence of Boltzmann factor for lower state energy
  const double blt_fct_tpt_scl(std::exp(ln_nrg_lwr_J*((1.0/tpt_HITRAN)-(1.0/tpt_mdp))/cst_Boltzmann)); // [frc] Temperatur dependence of Boltzmann factor for lower state energy
  ln_str_scl_wvn*=blt_fct_tpt_scl; // [cm-1 m2 mlc-1] Line strength at ambient temperature
  ln_str_scl_frq*=blt_fct_tpt_scl; // [Hz m2 mlc-1] Line strength at ambient temperature
  ln_str_scl_wvl*=blt_fct_tpt_scl; // [m m2 mlc-1] Line strength at ambient temperature

  // Scale line strength by temperature dependence of stimulated emission correction factor
  const double nrg_rcp_blt(-100.0*cst_Planck*speed_of_light*ln_ctr_wvn/cst_Boltzmann); // [K] Stimulated emission correction factor = h*frq/k = h*c*wvn/k = E/k
  const double stm_msn_tpt_scl((1.0-std::exp(nrg_rcp_blt/tpt_mdp))/(1.0-std::exp(nrg_rcp_blt/tpt_HITRAN))); // [frc] Stimulated emission correction factor
  ln_str_scl_wvn*=stm_msn_tpt_scl; // [cm-1 m2 mlc-1] Line strength at ambient temperature
  ln_str_scl_frq*=stm_msn_tpt_scl; // [Hz m2 mlc-1] Line strength at ambient temperature
  ln_str_scl_wvl*=stm_msn_tpt_scl; // [m m2 mlc-1] Line strength at ambient temperature

  if(true){
    std::cout << "  Molecule: " << mlc_sng[mlc_id] << ", Isotopomer = " << iso_sng[istpmr_idx] << std::endl;
    std::cout << "  MMW: " << mmw_gas << " kg mol-1, xpn_gas = " << xpn_gas << " " << std::endl;
    std::cout << "  prs_mdp = " << prs_mdp/100.0 << " mb, prs_dlt = " << prs_dlt/100.0 << " mb, T = " << tpt_mdp << " K" << std::endl;
    std::cout << "  vmr_gas = " << vmr_gas*1.0e6 << " ppmv, prs_gas = " << prs_gas << " Pa = " << prs_gas/100.0 << " mb" << std::endl;
    std::cout << "  mmr_gas = " << mmr_gas*1.0e6 << " mg kg-1, mpl_gas = " << mpl_gas*1.0e3 << " g m-2, mpc_gas = " << mpc_gas*1.0e3 << " g m-2" << std::endl;
    std::cout << "  Scaled to ambient p, T:" << std::endl;
    std::cout << "  ln_ctr: " << ln_ctr_wvn << " cm-1, " << ln_ctr_frq/1.0e12 << " THz, " << ln_ctr_wvl*1.0e6 << " um" << std::endl;
    std::cout << "  ln_sft_air: " << ln_sft_air_wvn << " cm-1, " << ln_sft_air_frq/1.0e9 << " MHz, " << ln_sft_air_wvl*1.0e9 << " nm" << std::endl;
    std::cout << "  ln_nrg_lwr: " << ln_nrg_lwr_wvn << " cm-1, " << ln_nrg_lwr_frq/1.0e12 << " THz, " << ln_nrg_lwr_wvl*1.0e6 << " um, " << ln_nrg_lwr_J << " J, " << ln_nrg_lwr_J*1.0e7 << " erg, " << ln_nrg_lwr_J/joules_per_eV << " eV"<< std::endl;
    std::cout << "  HWHM_air: " << HWHM_air_wvn << " cm-1, " << HWHM_air_frq/1.0e9 << " GHz, " << HWHM_air_wvl*1.0e9 << " nm" << std::endl;
    std::cout << "  HWHM_slf: " << HWHM_slf_wvn << " cm-1, " << HWHM_slf_frq/1.0e9 << " GHz, " << HWHM_slf_wvl*1.0e9 << " nm" << std::endl;
    std::cout << "  HWHM_prs: " << HWHM_prs_wvn << " cm-1, " << HWHM_prs_frq/1.0e9 << " GHz, " << HWHM_prs_wvl*1.0e9 << " nm" << std::endl;
    std::cout << "  HWEM_dpp: " << HWEM_dpp_wvn << " cm-1, " << HWEM_dpp_frq/1.0e9 << " GHz, " << HWEM_dpp_wvl*1.0e9 << " nm" << std::endl;
    std::cout << "  HWHM_dpp: " << HWHM_dpp_wvn << " cm-1, " << HWHM_dpp_frq/1.0e9 << " GHz, " << HWHM_dpp_wvl*1.0e9 << " nm" << std::endl;
    std::cout << "  HWHM_vgt: " << HWHM_vgt_wvn << " cm-1, " << HWHM_vgt_frq/1.0e9 << " GHz, " << HWHM_vgt_wvl*1.0e9 << " nm" << std::endl;
    std::cout << "  tpt_fct = " << tpt_fct << " frc" << std::endl;
    std::cout << "  ln_str = " << ln_str_htrn << " cm-1 cm2 mlc-1, " << ln_str_htrn_wvn << " cm-1 m2 mlc-1, " << ln_str_htrn_frq << " Hz m2 mlc-1, " << ln_str_htrn_wvl << " m m2 mlc-1" << std::endl;
    std::cout << "  Scaling of line strength from T = " << tpt_HITRAN << " K to T = " << tpt_mdp << " K:" << std::endl;
    std::cout << "  Partition function = " << prt_fnc_tpt_scl << ", Boltzmann factor = " << blt_fct_tpt_scl << ", stimulated emission = " << stm_msn_tpt_scl << std::endl;
    std::cout << "  Scaled line strength: " << ln_str_scl_wvn << " cm-1 m2 mlc-1, " << ln_str_scl_frq << " Hz m2 mlc-1, " << ln_str_scl_wvl << " m m2 mlc-1" << std::endl;
  } // endif true

  /* fxm: Define these once standardized on particular units
  // Ambient values of HITRAN variables (scaled to ambient T and p)
  double HWHM_air; // [m] Air-broadened HWHM
  double HWHM_slf; // [m] Self-broadened HWHM
  double ln_ctr; // [m] Wavelength at line center
  double ln_nrg_lwr; // [J] Lower state energy of transition
  double ln_sft_air; // [m] Air-broadened pressure shift of line transition
  double ln_str; // [Hz m2 mlc-1] Line strength
  */

  double chp_fnc; // [m] Chapman function of solar zenith angle
  // fxm: Adopt more sophisticated Chapman scheme a la David Huestis' chp_Hue00.f
  if(slr_zen_ngl_cos > 0.0) chp_fnc=1.0/slr_zen_ngl_cos; else chp_fnc=0.0; // [m] Chapman function of solar zenith angle
  if(true){
    std::cout << "  Cosine zenith angle = " << slr_zen_ngl_cos << ", Chapman function = " << chp_fnc << std::endl;
  } // endif true

  // Line transmittances
  const long wvl_nbr=wvlgrd.wvl_nbr_get(); // [nbr] Number of wavelength bins
  double *ln_ctr_dlt_frq=new double[wvl_nbr]; // [Hz] Distance from line center
  double *ln_ctr_dlt_wvn=new double[wvl_nbr]; // [cm-1] Distance from line center
  double *ln_ctr_dlt_wvl=new double[wvl_nbr]; // [m] Distance from line center
  double *lnshp_dpp_frq=new double[wvl_nbr]; // [Hz-1] Doppler line shape profile
  double *lnshp_dpp_wvn=new double[wvl_nbr]; // [cm] Doppler line shape profile
  double *lnshp_prs_frq=new double[wvl_nbr]; // [Hz-1] Pressure-broadened line shape profile
  double *lnshp_prs_wvn=new double[wvl_nbr]; // [cm] Pressure-broadened line shape profile
  double *lnshp_vgt_frq=new double[wvl_nbr]; // [Hz-1] Voigt line shape profile
  double *lnshp_vgt_wvn=new double[wvl_nbr]; // [cm] Voigt line shape profile
  double frq_ctr; // [Hz] Frequency at band center
  double wvl_ctr; // [m] Wavelength at band center
  double wvn_ctr; // [cm-1] Wavenumber at band center
  double *tau_abs_clm=new double[wvl_nbr]; // [frc] Column absorption optical depth 
  double *trn_clm=new double[wvl_nbr]; // [frc] Column transmittance of line
  using mth::cst_M_PIl; // (3.1415926535897932384626433832795029L) [frc] 3
  for(idx=0;idx<wvl_nbr;idx++){
    // Assemble current wavenumber, frequency, wavelength information
    wvn_ctr=wvlgrd.wvn_ctr_get()[idx]; // [cm-1] Wavenumber at band center
    frq_ctr=wvlgrd.frq_ctr_get()[idx]; // [Hz] Frequency at band center
    wvl_ctr=wvlgrd.wvl_ctr_get()[idx]; // [m] Wavelength at band center
    ln_ctr_dlt_wvn[idx]=wvn_ctr-ln_ctr_wvn; // [cm-1] Distance from line center
    ln_ctr_dlt_frq[idx]=frq_ctr-ln_ctr_frq; // [Hz] Distance from line center
    ln_ctr_dlt_wvl[idx]=wvl_ctr-ln_ctr_wvl; // [m] Distance from line center

    // Doppler profile
    lnshp_dpp_wvn[idx]=lnshp_dpp_fst_scl(HWEM_dpp_wvn,ln_ctr_dlt_wvn[idx]); // [cm] Doppler line shape profile
    lnshp_dpp_frq[idx]=lnshp_dpp_fst_scl(HWEM_dpp_frq,ln_ctr_dlt_frq[idx]); // [Hz-1] Doppler line shape profile

    // Pressure-broadening
    lnshp_prs_wvn[idx]=lnshp_lrn_fst_scl(HWHM_prs_wvn,ln_ctr_dlt_wvn[idx]); // [cm] Pressure-broadened line shape profile
    lnshp_prs_frq[idx]=lnshp_lrn_fst_scl(HWHM_prs_frq,ln_ctr_dlt_frq[idx]); // [Hz-1] Pressure-broadened line shape profile

    // Voigt profile
    lnshp_vgt_wvn[idx]=lnshp_vgt_fst_scl(HWEM_dpp_wvn,HWHM_prs_wvn,ln_ctr_dlt_wvn[idx]); // [cm] Voigt line shape profile
    lnshp_vgt_frq[idx]=lnshp_vgt_fst_scl(HWEM_dpp_frq,HWHM_prs_frq,ln_ctr_dlt_frq[idx]); // [Hz-1] Voigt line shape profile

    // Column properties
    tau_abs_clm[idx]=ln_str_scl_wvn*mpc_gas*cst_Avagadro/(2.0*mth::cst_M_PIl*HWHM_prs_wvn*mmw_gas); // [frc] Absorption optical depth
    tau_abs_clm[idx]*=std::log((wvn_ctr*wvn_ctr+HWHM_prs_wvn*HWHM_prs_wvn)/(wvn_ctr*wvn_ctr+HWEM_dpp_wvn*HWEM_dpp_wvn)); // [frc] Absorption optical depth ThS99 p. 400 (10.50)
    trn_clm[idx]=std::exp(-tau_abs_clm[idx]); // [frc] Transmittance of line
  } // end loop over wvn

  // Vector Lorentz profile
  rcd+=lnshp_lrn_vct // [fnc] Lorentz line shape profile
    (wvl_nbr, // I [nbr] Size of arrays
     HWHM_prs_wvn, // I [Hz-1,cm] Lorentz half-width at half-maximum
     ln_ctr_dlt_wvn, // I [Hz-1,cm] Distance from line center
     lnshp_prs_wvn); // O [Hz-1,cm] Lorentz line shape profile

  // Vector Doppler profile
  rcd+=lnshp_dpp_vct // [fnc] Doppler line shape profile
    (wvl_nbr, // I [nbr] Size of arrays
     HWEM_dpp_wvn, // I [Hz-1,cm] Doppler width (half-width at (1/e)-maximum)
     ln_ctr_dlt_wvn, // I [Hz-1,cm] Distance from line center
     lnshp_dpp_wvn); // O [Hz-1,cm] Doppler line shape profile

  // Vector Voigt profile
  rcd+=lnshp_vgt_vct // [fnc] Voigt line shape profile
    (wvl_nbr, // I [nbr] Size of arrays
     HWEM_dpp_wvn, // I [Hz,cm-1] Doppler width (half-width at (1/e)-maximum)
     HWHM_prs_wvn, // I [Hz,cm-1] Lorentz half-width at half-maximum
     ln_ctr_dlt_wvn, // I [Hz,cm-1] Distance from line center
     lnshp_vgt_wvn); // O [Hz-1,cm] Doppler line shape profile

  // Layer properties
  double *tau_abs_lyr=new double[wvl_nbr]; // [frc] Layer absorption optical depth
  double *abs_cff_mss_gas=new double[wvl_nbr]; // [m2 kg-1] Mass absorption coefficient
  double *trn_lyr=new double[wvl_nbr]; // [frc] Layer transmittance of line
  for(idx=0;idx<wvl_nbr;idx++){
    abs_cff_mss_gas[idx]=ln_str_scl_wvn*lnshp_vgt_wvn[idx]; // [m2 kg-1] Absorption cross-section
    tau_abs_lyr[idx]=abs_cff_mss_gas[idx]*mpl_gas; // [frc] Absorption optical depth
    trn_lyr[idx]=std::exp(-tau_abs_lyr[idx]); // [frc] Transmittance of line
  } // end loop over wvn

  if(true){
    // fxm: Add absorption cross section
    (void)std::fprintf(stdout,"idx\twvn_ctr\twvl_ctr\ttau_abs\ttrn_clm\tphi_prs\tphi_dpp\tphi_vgt\n");
    (void)std::fprintf(stdout,"   \t cm-1  \t  um   \t  frc  \t  frc \t  cm   \t  cm   \t  cm   \n");
    for(idx=0;idx<wvl_nbr;idx++)
      (void)std::fprintf(stdout,
		    "%4ld\t%.3f\t%.4f\t%.1e\t%.4f\t%.3f\t%.3f\t%.3f\n",
		    idx,wvlgrd.wvn_ctr_get()[idx],wvlgrd.wvl_ctr_get()[idx]*1.0e6,tau_abs_clm[idx],trn_clm[idx],lnshp_prs_wvn[idx],lnshp_dpp_wvn[idx],lnshp_vgt_wvn[idx]);
  } // endif true

  if(true){
    (void)std::fprintf(stdout,"idx\twvn_ctr\twvl_ctr\ttau_abs\ttrn_lyr\tphi_prs\tphi_dpp\tphi_vgt\n");
    (void)std::fprintf(stdout,"   \t cm-1  \t  um   \t  frc  \t  frc \t  cm   \t  cm   \t  cm   \n");
    for(idx=0;idx<wvl_nbr;idx++)
      (void)std::fprintf(stdout,
		    "%4ld\t%.3f\t%.4f\t%.1e\t%.4f\t%.4f\t%.4f\t%.4f\n",
		    idx,wvlgrd.wvn_ctr_get()[idx],wvlgrd.wvl_ctr_get()[idx]*1.0e6,tau_abs_lyr[idx],trn_lyr[idx],lnshp_prs_wvn[idx],lnshp_dpp_wvn[idx],lnshp_vgt_wvn[idx]);
  } // endif true

  const long mlc_nbr(mlc_nbr_max_htrn); // [nbr] Number of molecules
  const long iso_nbr(iso_nbr_max_htrn); // [nbr] Number of isotopomers
  size_t mlc_nm_lng_max(0); // [nbr] Maximum length of molecule name
  size_t iso_nm_lng_max(0); // [nbr] Maximum length of isotopomer name
  for(mlc_idx=0;mlc_idx<mlc_nbr;mlc_idx++){
    if(mlc_sng[mlc_idx].length() > mlc_nm_lng_max) mlc_nm_lng_max=mlc_sng[mlc_idx].length(); // [nbr] Maximum length of molecule name
  } /* end loop over molecules */
  for(iso_idx=0;iso_idx<iso_nbr;iso_idx++){
    if(iso_sng[iso_idx].length() > iso_nm_lng_max) iso_nm_lng_max=iso_sng[iso_idx].length(); // [nbr] Maximum length of molecule name
  } /* end loop over isotopomers */

  // Output results of module
  // Allow file to already be in define mode
  rcd=nco_redef(nc_out,NC_EINDEFINE); // [fnc] Put open netCDF dataset into define mode
  // Get existing dimensions
  const int wvl_dmn(nco_inq_dimid(nc_out,static_cast<std::string>("wvl"))); // [dmn] Wavelength dimension
  const int *dmn_wvl(&wvl_dmn); // Pointer to wavelength dimension
  const int *dmn_scl(&wvl_dmn); // Dummy argument, not used
  // Define new dimensions
  const int mlc_dmn(nco_def_dim(nc_out,static_cast<std::string>("mlc"),mlc_nbr)); // [dmn] Molecule dimension
  const int iso_dmn(nco_def_dim(nc_out,static_cast<std::string>("iso"),iso_nbr)); // [dmn] Isotopomer dimension
  const int *mlc_dmn_ptr(&mlc_dmn); // Pointer to molecule dimension
  const int *iso_dmn_ptr(&iso_dmn); // Pointer to isotopomer dimension
  // fxm: must define sngarr2ptrarr before writing 2-dimensional string arrays with C++ netCDF interface
  //  const int mlc_sng_dmn(nco_def_dim(nc_out,static_cast<std::string>("mlc_sng"),mlc_nm_lng_max)); // [dmn] Molecule string length dimension
  // const int iso_sng_dmn(nco_def_dim(nc_out,static_cast<std::string>("iso_sng"),iso_nm_lng_max)); // [dmn] Isotopomer string length dimension
  //  const int dmn_mlc_mlc_sng[2]={mlc_dmn,mlc_sng_dmn};
  //  const int dmn_iso_iso_sng[2]={iso_dmn,iso_sng_dmn};

  //  char **mlc_sng_chr_ptr=sngarr2ptrarr(mlc_sng,sizeof(mlc_sng)/sizeof(mlc_sng[0])); // [ptr] List of pointers to C-strings
  //  char **iso_sng_chr_ptr=sngarr2ptrarr(iso_sng,sizeof(iso_sng)/sizeof(iso_sng[0])); // [ptr] List of pointers to C-strings

  var_mtd_sct var_mtd[]={
    //    {0,"iso_sng",NC_CHAR,2,dmn_iso_iso_sng,"long_name","HITRAN isotopomer names","units","string"},
    //    {0,"mlc_sng",NC_CHAR,2,dmn_mlc_mlc_sng,"long_name","HITRAN molecule names","units","string"},
    {0,"abs_cff_mss_gas",NC_DOUBLE,1,dmn_wvl,"long_name","Mass absorption coefficient","units","meter2 kilogram-1"},
    {0,"chp_fnc",NC_DOUBLE,0,dmn_scl,"long_name","Chapman function","units","fraction"},
    {0,"ln_str_htrn_frq",NC_DOUBLE,0,dmn_scl,"long_name","Line Strength","units","hertz meter2 molecule-1"},
    {0,"ln_str_htrn_wvn",NC_DOUBLE,0,dmn_scl,"long_name","Line Strength","units","centimeter-1 meter2 molecule-1"},
    {0,"ln_str_scl_frq",NC_DOUBLE,0,dmn_scl,"long_name","Line Strength","units","hertz meter2 molecule-1"},
    {0,"ln_str_scl_wvn",NC_DOUBLE,0,dmn_scl,"long_name","Line Strength","units","centimeter-1 meter2 molecul-1"},
    {0,"lnshp_dpp_wvn",NC_DOUBLE,1,dmn_wvl,"long_name","Doppler line shape profile","units","centimeter"},
    {0,"lnshp_prs_wvn",NC_DOUBLE,1,dmn_wvl,"long_name","Pressure-broadened line shape profile","units","centimeter"},
    {0,"lnshp_vgt_wvn",NC_DOUBLE,1,dmn_wvl,"long_name","Voigt line shape profile","units","centimeter"},
    {0,"mmw_iso",NC_DOUBLE,1,iso_dmn_ptr,"long_name","Mean molecular weight of isotopomer","units","kilogram mole-1"},
    {0,"mmw_mlc",NC_DOUBLE,1,mlc_dmn_ptr,"long_name","Mean molecular weight of molecule","units","kilogram mole-1"},
    {0,"tau_abs_clm",NC_DOUBLE,1,dmn_wvl,"long_name","Column absorption optical depth","units","fraction"},
    {0,"tau_abs_lyr",NC_DOUBLE,1,dmn_wvl,"long_name","Layer absorption optical depth","units","fraction"},
    {0,"trn_clm",NC_DOUBLE,1,dmn_wvl,"long_name","Column transmittance of line","units","fraction"},
    {0,"trn_lyr",NC_DOUBLE,1,dmn_wvl,"long_name","Layer transmittance of line","units","fraction"},
    {0,"xpn_mlc",NC_DOUBLE,1,mlc_dmn_ptr,"long_name","Exponent defining temperature dependence of rotational partition function","units","fraction"},
  }; // end var_mtd_sct var_mtd[]
  const int var_mtd_nbr(sizeof(var_mtd)/sizeof(var_mtd_sct)); // [nbr] Number of variables in array
  rcd+=nco_var_dfn(nc_out,var_mtd,var_mtd_nbr,dmn_nbr_max); // [fnc] Define variables in output netCDF file

  // After writing, delete arrays 

  /* cvt: Offset following four arrays by one to skip zeroth value when writing
     This keeps HITRAN nomenclature intact with Fortran indexing */
  //  rcd=nco_put_var(nc_out,const_cast<std::string>("mlc_sng"),mlc_sng_chr_ptr+1); delete []mlc_sng;
  //  rcd=nco_put_var(nc_out,const_cast<std::string>("iso_sng"),iso_sng_chr_ptr+1); delete []iso_sng;
  rcd=nco_put_var(nc_out,static_cast<std::string>("xpn_mlc"),xpn_mlc+1); delete []xpn_mlc;
  rcd=nco_put_var(nc_out,static_cast<std::string>("mmw_mlc"),mmw_mlc+1); delete []mmw_mlc;
  rcd=nco_put_var(nc_out,static_cast<std::string>("mmw_iso"),mmw_iso+1); delete []mmw_iso;
  // end offset arrays
  rcd=nco_put_var(nc_out,static_cast<std::string>("abs_cff_mss_gas"),abs_cff_mss_gas); delete []abs_cff_mss_gas;
  rcd=nco_put_var(nc_out,static_cast<std::string>("chp_fnc"),chp_fnc);
  rcd=nco_put_var(nc_out,static_cast<std::string>("ln_str_htrn_frq"),ln_str_htrn_frq);
  rcd=nco_put_var(nc_out,static_cast<std::string>("ln_str_htrn_wvn"),ln_str_htrn_wvn);
  rcd=nco_put_var(nc_out,static_cast<std::string>("ln_str_scl_frq"),ln_str_scl_frq);
  rcd=nco_put_var(nc_out,static_cast<std::string>("ln_str_scl_wvn"),ln_str_scl_wvn);
  rcd=nco_put_var(nc_out,static_cast<std::string>("lnshp_dpp_wvn"),lnshp_dpp_wvn); delete []lnshp_dpp_wvn;
  rcd=nco_put_var(nc_out,static_cast<std::string>("lnshp_prs_wvn"),lnshp_prs_wvn); delete []lnshp_prs_wvn;
  rcd=nco_put_var(nc_out,static_cast<std::string>("lnshp_vgt_wvn"),lnshp_vgt_wvn); delete []lnshp_vgt_wvn;
  rcd=nco_put_var(nc_out,static_cast<std::string>("tau_abs_clm"),tau_abs_clm); delete []tau_abs_clm;
  rcd=nco_put_var(nc_out,static_cast<std::string>("tau_abs_lyr"),tau_abs_lyr); delete []tau_abs_lyr;
  rcd=nco_put_var(nc_out,static_cast<std::string>("trn_clm"),trn_clm); delete []trn_clm;
  rcd=nco_put_var(nc_out,static_cast<std::string>("trn_lyr"),trn_lyr); delete []trn_lyr;

  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Exiting...");
  return rcd; // [enm] Return success code
} // end rt_lbl()

// Global functions with C++ linkages end
