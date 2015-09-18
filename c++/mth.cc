// $Id$ 

// Purpose: Mathematical utilities, constants for C++ programs

/* Copyright (C) 1997--2014 Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

#include <mth.hh> // Mathematical utilities, constants

int // O [rcd] Return success code
eqn_qdr_slv // [fnc] Solve quadratic equation
(const prc_cmp prm_a, // I [frc] Parameter a of quadratic equation
 const prc_cmp prm_b, // I [frc] Parameter b of quadratic equation
 const prc_cmp prm_c, // I [frc] Parameter c of quadratic equation
 prc_cmp *sln_1, // O [frc] First root of quadratic equation
 prc_cmp *sln_2) // O [frc] Second root of quadratic equation
{
  // Purpose: Solve quadratic equation a*x^2+b*x+c=0
  // Output
  int rcd(0); // O [rcd] Return success code
  // Local
  const std::string sbr_nm("eqn_qdr_slv"); // [sng] Name of subroutine
  prc_cmp qdr_fct; // [frc] Factor in solution of quadratic equation
  prc_cmp qdr_dsc; // [frc] Discriminant of quadratic equation
  // Main code
  if(prm_a == 0.0){
    if(prm_b == 0.0){
      err_prn(sbr_nm,"a = b = 0.0 in quadratic equation");
    } // endif prm_a
    *sln_1=*sln_2=-prm_b/prm_c;
    return rcd;
  } // endif prm_a
  qdr_dsc=std::sqrt(prm_b*prm_b-4.0*prm_a*prm_c); // [frc] Discriminant of quadratic equation
  qdr_fct=-0.5*(prm_b+sgn(prm_b)*qdr_dsc); // [frc] Factor in solution of quadratic equation PTV96 p. 178 (5.6.4)
  *sln_1=prm_c/qdr_fct; // O [frc] First root of quadratic equation
  *sln_2=qdr_fct/prm_a; // O [frc] Second root of quadratic equation
  return rcd;
} // end eqn_qdr_slv()

std::complex<prc_cmp> // [frc] Complex value
sng2cpx // [fnc] Convert string to complex value
(const std::string &sng) // [sng] Input string
{
  /* Purpose: Convert input string to complex number
     Following formats are correctly parsed as complex numbers:
     Let u = units digits, d = decimal digits, i = imaginary unit
     [+-]u.d [+-]u.di [+-]u.d[+-]u.di (r,i)
     Following formats are not yet correctly parsed:
     [+-]u.d[+-]u.m[eEdD]xpni (r,i)
  */

  const std::string sbr_nm("sng2cpx"); // [sng] Name of subroutine
  std::string rl_sng("0.0");
  std::string img_sng("0.0");

  unsigned long i_psn=sng.find("i"); // [nbr] Position of "i"
  unsigned long xpn_frs_psn=sng.find_first_of("eEdD"); // [nbr] Position of first exponent indicator
  unsigned long xpn_scn_psn(std::string::npos); // [nbr] Position of second exponent indicator
  unsigned long pm_frs_psn=sng.find_first_of("+-"); // [nbr] Position of first "+" or "-"
  unsigned long pm_scn_psn(std::string::npos); // [nbr] Position of second "+" or "-"
  unsigned long pm_thr_psn(std::string::npos); // [nbr] Position of third "+" or "-"
  unsigned long pm_frt_psn(std::string::npos); // [nbr] Position of fourth "+" or "-"
  unsigned long nbr_frs_psn=sng.find_first_of(".1234567890"); // [nbr] Position of first digit or decimal point
  unsigned long lft_prn_psn=sng.find_first_of("("); // [nbr] Position of "("
  unsigned long rgt_prn_psn=sng.find_first_of(")"); // [nbr] Position of ")"
  unsigned long cmm_psn=sng.find_first_of(","); // [nbr] Position of ","

  if(pm_frs_psn != std::string::npos) pm_scn_psn=sng.find_first_of("+-",pm_frs_psn+1UL);
  if(pm_scn_psn != std::string::npos) pm_thr_psn=sng.find_first_of("+-",pm_scn_psn+1UL);
  if(pm_thr_psn != std::string::npos) pm_frt_psn=sng.find_first_of("+-",pm_thr_psn+1UL);
  if(xpn_frs_psn != std::string::npos) xpn_scn_psn=sng.find_first_of("eEdD",xpn_frs_psn+1UL);

  if(dbg_lvl_get() == dbg_old){
    std::cerr << "i_psn = " << i_psn << std::endl;
    std::cerr << "xpn_frs_psn = " << xpn_frs_psn << std::endl;
    std::cerr << "pm_frs_psn = " << pm_frs_psn << std::endl;
    std::cerr << "pm_scn_psn = " << pm_scn_psn << std::endl;
    std::cerr << "pm_thr_psn = " << pm_thr_psn << std::endl;
    std::cerr << "nbr_frs_psn = " << nbr_frs_psn << std::endl;
  } // endif dbg

  // Sanity check
  if(nbr_frs_psn == std::string::npos){ 
    err_prn(sbr_nm,"String \""+sng+"\" is not a number");
  } // endif
  if(lft_prn_psn != std::string::npos && cmm_psn != std::string::npos && rgt_prn_psn != std::string::npos){ // (real,complex) format, e.g., (0.5,0.5)
    rl_sng=sng.substr(lft_prn_psn+1,cmm_psn-lft_prn_psn-1);
    img_sng=sng.substr(cmm_psn+1,rgt_prn_psn-cmm_psn-1);
  }else if(i_psn == std::string::npos){ // No "i"
    // e.g., 1.33
    rl_sng=sng;
  }else{ // There is an "i"
    // In following calls to substr(), count is simply end position + 1 minus start
    if(pm_frt_psn != std::string::npos){ // There are two unary +/-'s and two +/-'s for the exponents
      // e.g., +3.0e-37-4.0e-37i
      rl_sng=sng.substr(0,pm_thr_psn);
      img_sng=sng.substr(pm_thr_psn,i_psn-pm_thr_psn);
    }else if(pm_thr_psn != std::string::npos){ // There are three +/-'s
      if(xpn_scn_psn != std::string::npos){ // There is one unary +/- and two +/-'s for exponents
	// e.g., 3.0e-37-4.0e-37i
	rl_sng=sng.substr(0,pm_scn_psn);
	img_sng=sng.substr(pm_scn_psn,i_psn-pm_scn_psn);
      }else{
	if(xpn_frs_psn == pm_thr_psn-1UL){ // There are two unary +/-'s and a +/- for the imaginary exponent
	  // e.g., +3-4.0e-37i
	  rl_sng=sng.substr(pm_frs_psn,pm_scn_psn-pm_frs_psn);
	  img_sng=sng.substr(pm_scn_psn,i_psn-pm_scn_psn);
	}else if(xpn_frs_psn == pm_thr_psn-1UL){ // There are two unary +/-'s and a +/- for the real exponent
	  // e.g., +3.0e-37-4.0i
	  rl_sng=sng.substr(pm_frs_psn,pm_thr_psn-pm_frs_psn);
	  img_sng=sng.substr(pm_thr_psn,i_psn-pm_thr_psn);
	} // endelse
      } // endif three +/-'s
    }else if(pm_scn_psn != std::string::npos){ // There are two +/-'s
      if(xpn_frs_psn == pm_scn_psn-1UL || xpn_scn_psn == pm_scn_psn-1UL){ // There is one unary +/- and a +/- for the imaginary exponent
	// e.g., 3-4.0e-37i
	rl_sng=sng.substr(0,pm_frs_psn);
	img_sng=sng.substr(pm_frs_psn,i_psn-pm_frs_psn);
      }else if(xpn_frs_psn == pm_frs_psn-1UL){ // There is one unary +/- and a +/- for the real exponent
	// e.g., 3.0e-37-4.0i
	rl_sng=sng.substr(0,pm_frs_psn);
	img_sng=sng.substr(pm_frs_psn,i_psn-pm_frs_psn);
      }else if(pm_frs_psn != pm_thr_psn){ // There are two unary +/-'s
	// e.g., +3-4i
	rl_sng=sng.substr(pm_frs_psn,pm_thr_psn-pm_frs_psn);
	img_sng=sng.substr(pm_thr_psn,i_psn-pm_thr_psn);
      } // endif two +/-'s
    }else if(pm_frs_psn != std::string::npos){ // One unary +/-
      if (nbr_frs_psn < pm_frs_psn){
	// e.g., 3+4i
	rl_sng=sng.substr(0,pm_frs_psn);
	img_sng=sng.substr(pm_frs_psn,i_psn-pm_frs_psn);
      }else{
	// e.g., +4i
	img_sng=sng.substr(pm_frs_psn,i_psn-pm_frs_psn);
      } // endif
    }else{ // Zero unary +/-
      // e.g., .1i
      img_sng=sng.substr(nbr_frs_psn,i_psn-nbr_frs_psn);
    } // endelse
  } // endif there is an "i"
  // C++ static_cast operator 
  std::complex<double> dcx(std::strtod(rl_sng.c_str(),(char **)NULL),std::strtod(img_sng.c_str(),(char **)NULL));     
  std::complex<prc_cmp> rcx=static_cast<std::complex<prc_cmp> >(dcx);
  if(dbg_lvl_get() == dbg_old) std::cerr << sbr_nm << "() translates " << sng << " into " << rcx << std::endl;
  
  return rcx;
} // end sng2cpx()

prc_cmp // O [frc] Return function value
fnc_nwt_rph // [fnc] Function value and the overloaded solver
(const std::string eqn_sng, // [sng] Equation to solve
 const prc_cmp abc)
{
  /* Purpose: Function computer for Newton-Raphson root finder
     Usage: Define function so that LHS = 0, e.g.,
     x^2 - 1 = 0 */
  // Local
  prc_cmp wrk_tmp_1; // [frc] Workspace for function evaluation
  prc_cmp wrk_tmp_2; // [frc] Workspace for function evaluation
  // Output
  prc_cmp fnc_val(CEWI_cpv); // [frc] Return function value
  
  const std::string sbr_nm("fnc_nwt_rph"); // [sng] Name of subroutine
  if(eqn_sng == "tst"){ // Test equation x^2 - 1 = 0
    fnc_val=abc*abc-1; // [frc] Function value
  }else if(eqn_sng == "lnstr_rat"){
    /* Transcendental equation for line strength ratio in Malkmus distribution
       abc represents lnstrrat */
    if(abc < 0.0) err_prn(prg_nm_get(),sbr_nm,"Square root of negative argument prevention");
    wrk_tmp_1=std::sqrt(abc); // [frc]
    if(wrk_tmp_1 == -1.0) err_prn(prg_nm_get(),sbr_nm,"Divide by zero prevention 1");
    wrk_tmp_2=std::exp(4.0*(wrk_tmp_1-1)/(wrk_tmp_1+1)); // [frc]
    fnc_val=abc-wrk_tmp_2; // [frc] Function value
  }else if(eqn_sng == "sci_prg"){ // x*exp(x) - 10 = 0
    if(abc == -1.0) err_prn(prg_nm_get(),sbr_nm,"Divide by zero prevention 1");
    fnc_val=abc*std::exp(abc)-10;
  }else if(eqn_sng == "wien"){
    /* Transcendental equation for h*nu/(k*tpt) in deriving Wien's displacement law */
    wrk_tmp_1=std::exp(-abc); // [frc] Workspace for function evaluation
    fnc_val=abc-3.0*(1-wrk_tmp_1); // [frc] Function value
  }else{
    err_prn(prg_nm_get(),sbr_nm,"Unknown equation "+eqn_sng);
  } // endelse

  return fnc_val; // [enm] Return success code
} // end fnc_nwt_rph()

prc_cmp // Return derivative value
drv_fnc_nwt_rph // [fnc] Derivative value
(const std::string eqn_sng, // [sng] Equation to solve
 const prc_cmp abc) // [frc] I Current value of root guess
{
  /* Purpose: Derivative computer for Newton-Raphson root finder
     Usage: Define function so that LHS = 0, e.g., x^2 - 1 = 0 */
  // Output
  prc_cmp drv_val(CEWI_cpv); // [frc] Derivative value

  // Local
  prc_cmp wrk_tmp_1; // [frc] Workspace for function evaluation
  prc_cmp wrk_tmp_2; // [frc] Workspace for function evaluation
  prc_cmp wrk_tmp_3; // [frc] Workspace for function evaluation

  const std::string sbr_nm("drv_fnc_nwt_rph"); // [sng] Name of subroutine
  if(eqn_sng == "tst"){
    /* Test equation x^2 -1 = 0 */
    drv_val=2.0*abc; // [frc] Derivative value

  }else if(eqn_sng == "lnstr_rat"){
    /* Transcendental equation for line strength ratio in Malkmus distribution
       abc represents lnstrrat */
    if(abc < 0.0) err_prn(prg_nm_get(),sbr_nm,"square root of negative argument prevention");
    wrk_tmp_1=std::sqrt(abc); // [frc] Workspace for function evaluation
    // if(wrk_tmp_1 == -1.0) err_prn(prg_nm_get(),sbr_nm,"Divide by zero prevention 1");
    // wrk_tmp_1 was checked in fnc_wt_rph
    wrk_tmp_2=std::exp(4.0*(wrk_tmp_1-1)/(wrk_tmp_1+1)); // [frc] Workspace for function evaluation
    if(wrk_tmp_1 == 0.0) err_prn(prg_nm_get(),sbr_nm,"Divide by zero prevention 2");
    wrk_tmp_3=4.0/(wrk_tmp_1*(wrk_tmp_1+1.0)*(wrk_tmp_1+1.0)); // [frc] Workspace for function evaluation
    drv_val=1.0-4.0*wrk_tmp_3*wrk_tmp_2; // [frc] Derivative value
  }else if(eqn_sng == "wien"){
   /* Transcendental equation for h*nu/(k*tpt) in deriving Wien's displacement law */
    wrk_tmp_1=std::exp(-abc); // [frc] Workspace for function evaluation
    drv_val=1.0-3.0*abc*wrk_tmp_1; // [frc] Derivative value
  }else if(eqn_sng == "sci_prg"){ // solve x*exp(x)=10
    if(abc == -1.0) err_prn(prg_nm_get(),sbr_nm,"Divide by zero prevention 1");
    drv_val=std::exp(abc)*(1+abc); // [frc] Derivative value
  }else{
    err_prn(prg_nm_get(),sbr_nm,"Unknown equation "+eqn_sng);
  } // endelse

  return drv_val; //  Return derivative value
} // end drv_fnc_nwt_rph()
