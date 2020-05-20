// $Id$ 

// Purpose: Mathematical utilities, constants

/* Copyright (C) 1997--present Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

/* Header includes mathematical constants in mth_cst.hh
   If routine requires utility functions AND constants, then #include mth.hh
   If routine requires only constants, then #include mth_cst.hh
   Usage:
   #include <mth.hh> // Mathematical utilities, constants */

#ifndef MTH_HH // Contents have not yet been inserted in current source file  
#define MTH_HH

// C++ headers
#include <iostream> // Standard C++ I/O streams cout, cin
#include <string> // Standard C++ string class
#include <typeinfo> // typeid()

// Standard C headers
#include <cmath> // sin cos cos sin 3.14159 
#include <cstdio> // stderr, FILE, NULL, etc.
#include <cstdlib> // abort, exit, getopt, malloc, strtod, strtol

// Personal headers
#include <dbg.hh> // Debugging constants
#include <utl.hh> // err_prn(), wrn_prn(), Exit_gracefully()
#include <mth_cst.hh> // Mathematical constants, cst_M_PIl, cst_M_El...

/* Template definitions must be included in every file (translation unit)
   unless the compiler supports the "export" keyword (gcc 2.96 does not). 
   Thus actual function definitions must appear in .hh file being #include'd */

// Forward declarations

// Typedefs

// Declare functions with C++ linkages

// Neither of these appears to work
//template <> double sgn<double> ( double );
//template<> double sgn(double);
// fxm: C++ standard allows inlined templates, but does g++?
template<class prc_T>int sgn(const prc_T val); // O [fnc] +1 if positive, -1 if negative
template<class prc_T>int // O [sgn] +1 if positive, -1 if negative
sgn // [fnc] sgn() function: +1 if positive, -1 if negative
(const prc_T val) // Value whose sign will be tested
{
  /* Purpose: Signum function
     sgn() function returns +1 if argument is positive, -1 if argument is negative
     Note that http://planetmath.org/encyclopedia/SignumFunction.html defines 
     sgn(x) = 0 when x = 0.
     This should make no difference in any analysis since functional analysis
     should not depend on the value of a function at a single point.
     Hence, definition of sgn(0) is arbitrary and may be chosen for other reasons
     However, it does _appear_ to make a difference in quadratic equation solutions
     PTV96 have sgn(b) in what seems to be an important location
     This has never caused a problem in reality, though, so I assume that my
     understanding of sgn() is poor and the books are correct.
     For generality, use implicit conversion rules to promote input variable to complex type
     Using std::complex<prc_cmp> instead of std::complex<double> is fine since we are not concerned with precision */
  std::complex<prc_cmp> val_cpx(val); // Complex representation of value whose sign will be tested
  return (val_cpx.real() >= 0 ? 1 : -1); // [sgn] +1 if positive, -1 if negative
} // end sgn()

namespace mth_fnc{ // [nms] Mathematical function namespace
  template<class prc_T>int round(const prc_T val); // O [fnc] Round to nearest integer
  template<class prc_T>int // O [idx] Nearest integer
  round // [fnc] Round to nearest integer
  (const prc_T val) // [frc] Value to be rounded
  {
    /* Purpose: round() function returns the nearest integer to argument 
       Technique is taken from the USENET C FAQ
       Function will not work for complex input types */
    return static_cast<int>(val<0 ? val-0.5 : val+0.5); // [idx] Nearest integer
  } // end round()
} // Mathematical function namespace mth_fnc

std::complex<prc_cmp> // [frc] Complex value
sng2cpx // [fnc] Convert string to complex value
(const std::string &sng); // [sng] Input string

// Define inline'd functions in header so source is visible to calling files
inline int sgn_cpv(const prc_cmp val){return (val >= 0 ? 1 : -1);}

int // O [rcd] Return success code
eqn_qdr_slv // [fnc] Solve quadratic equation
(const prc_cmp prm_a, // I [frc] Parameter a of quadratic equation
 const prc_cmp prm_b, // I [frc] Parameter b of quadratic equation
 const prc_cmp prm_c, // I [frc] Parameter c of quadratic equation
 prc_cmp *sln_1, // O [frc] First root of quadratic equation
 prc_cmp *sln_2); // O [frc] Second root of quadratic equation

template<class prc_T>int // O [enm] Return success code
nwt_rph_slvr // [fnc] Newton-Raphson root finder
(const std::string eqn_sng, // [sng] Equation to solve
 const prc_T bnd_LHS, // I [frc] Left hand side boundary
 const prc_T bnd_RHS, // I [frc] Right hand side boundary
 const prc_T cnv_thr, // I [frc] Convergence threshold
 prc_T *sln); // O [frc] Solution of equation

template<class prc_T>int // O [enm] Return success code
nwt_rph_fnc // [fnc] Function value and derivative value
(const std::string eqn_sng, // [sng] Equation to solve
 const prc_T abc, // [frc] I Current value of root guess
 prc_T *fnc_val, // [frc] O Function value
 prc_T *drv_val); // [frc] O Derivative value

template<class prc_T>int // O [enm] Return success code
nwt_rph_fnc // [fnc] Function value and derivative value
(const std::string eqn_sng, // [sng] Equation to solve
 const prc_T abc, // [frc] I Current value of root guess
 prc_T *fnc_val, // [frc] O Function value
 prc_T *drv_val)  // [frc] O Derivative value
{
  /* Purpose: Function and derivative computer for Newton-Raphson root finder
     Usage: Define function so that LHS = 0, e.g., x^2 - 1 = 0 */
  // Output
  int rcd(0); // O [rcd] Return success code
  // Local
  prc_T wrk_tmp_1; // [frc] Workspace for function evaluation
  prc_T wrk_tmp_2; // [frc] Workspace for function evaluation
  prc_T wrk_tmp_3; // [frc] Workspace for function evaluation
  const std::string sbr_nm("nwt_rph_fnc"); // [sng] Name of subroutine
  if(eqn_sng == "tst"){
    /* Test equation x^2 - 1 = 0 */
    *fnc_val=abc*abc-1.0; // [frc] Function value
    *drv_val=2.0*abc; // [frc] Derivative value
  }else if(eqn_sng == "lnstr_rat"){
    /* Transcendental equation for line strength ratio in Malkmus distribution
       abc represents lnstrrat */
    if(abc < 0.0) err_prn(prg_nm_get(),sbr_nm,"Square root of negative argument prevention");
    wrk_tmp_1=std::sqrt(abc); // [frc] Workspace for function evaluation
    if(wrk_tmp_1 == -1.0) err_prn(prg_nm_get(),sbr_nm,"Divide by zero prevention 1");
    wrk_tmp_2=std::exp(4.0*(wrk_tmp_1-1)/(wrk_tmp_1+1)); // [frc] Workspace for function evaluation
    *fnc_val=abc-wrk_tmp_2; // [frc] Function value
    if(wrk_tmp_1 == 0.0) err_prn(prg_nm_get(),sbr_nm,"Divide by zero prevention 2");
    wrk_tmp_3=4.0/(wrk_tmp_1*(wrk_tmp_1+1.0)*(wrk_tmp_1+1.0)); // [frc] Workspace for function evaluation
    *drv_val=1.0-4.0*wrk_tmp_3*wrk_tmp_2; // [frc] Derivative value
  }else if(eqn_sng == "wien"){
    /* Transcendental equation for h*nu/(k*tpt) in deriving Wien's displacement law */
    wrk_tmp_1=std::exp(-abc); // [frc] Workspace for function evaluation
    *fnc_val=abc-3.0*(1-wrk_tmp_1); // [frc] Function value
    *drv_val=1.0-3.0*abc*wrk_tmp_1; // [frc] Derivative value
  }else if(eqn_sng == "sci_prg"){
    /* Test equation x*exp(x)-10=0 */ 
    wrk_tmp_1=std::exp(abc); // [frc] Workspace for function evaluation
    *fnc_val=abc*wrk_tmp_1-10.0; // [frc] Function value
    *drv_val=wrk_tmp_1+abc*wrk_tmp_1; // [frc] Derivative value
  }else{
    err_prn(prg_nm_get(),sbr_nm,"Unknown equation "+eqn_sng);
  } // endelse

  return rcd; // [enm] Return success code
} // end nwt_rph_fnc()

template<class prc_T>int // O [enm] Return success code
nwt_rph_fnc // [fnc] Function value and derivative value  - overload the original
(prc_cmp (*fnc_ptr)(const prc_cmp), // I [fnc] Pointer to function
 prc_cmp (*drv_ptr)(const prc_cmp), // I [fnc] Pointer to function derivative
 const prc_T abc, // I [frc] Current value of root guess
 prc_T *fnc_val, // O [frc] Function value
 prc_T *drv_val) // O [frc] Derivative value
{
  /* Purpose: Function and derivative computer for Newton-Raphson root finder
     Usage: Define function so that LHS = 0 */
  // Output
  int rcd(0); // O [rcd] Return success code
  // Local
  const std::string sbr_nm("nwt_rph_fnc"); // [sng] Name of subroutine
  
  *fnc_val=(*fnc_ptr)(abc); // [frc] Function value
  *drv_val=(*drv_ptr)(abc); // [frc] Derivative value
  
  return rcd; // [enm] Return success code
} // end overloaded nwt_rph_fnc()

template<class prc_T>int // O [enm] Return success code
nwt_rph_slvr // [fnc] Driver for Newton-Raphson root finder
(const std::string eqn_sng, // [sng] Equation to solve
 const prc_T bnd_LHS, // I [frc] Left hand side boundary
 const prc_T bnd_RHS, // I [frc] Right hand side boundary
 const prc_T cnv_thr, // I [frc] Convergence threshold
 prc_T *sln) // O [frc] Solution of equation
{
  /* Purpose: Driver for Newton-Raphson root finding
     Usage: Define function so that LHS = 0 */
  // Output
  int rcd(0); // O [rcd] Return success code
  const size_t itr_nbr(200); // [nbr] Maximum number of iterations to attempt
  // Local
  const std::string sbr_nm("nwt_rph_slvr"); // [sng] Name of subroutine
  size_t itr_idx; // [idx] Iteration index
  prc_T sln_dlt; // [frc] Adjustment for next guess
  prc_T fnc_val; // [frc] Function value
  prc_T drv_val; // [frc] Derivative value
  prc_T sln_old; // [frc] Solution of equation, old
  prc_T sln_dlt_old(CEWI_flt); // [frc] Adjustment for next guess, old
  prc_T fnc_val_old(CEWI_flt); // [frc] Function value, old

  if(dbg_lvl_get() > dbg_scl) std::cout << sbr_nm+"(): LHS Boundary = " << bnd_LHS << ", RHS Boundary = " << bnd_RHS << ", cnv_thr = " << cnv_thr << std::endl;
  *sln=0.5*(bnd_LHS+bnd_RHS); // [frc] Current value of root guess
  for(itr_idx=0;itr_idx<itr_nbr;itr_idx++){
    rcd+=nwt_rph_fnc // [fnc] Function value and derivative value
      (eqn_sng, // I [sng] Equation to solve
#if 0
       // Pass function pointers instead of string
       fnc_ptr, // I [fnc] Pointer to function
       drv_ptr, // I [fnc] Pointer to function derivative
#endif // endif 0
       *sln, // I [frc] Current value of root guess
       &fnc_val, // O [frc] Function value
       &drv_val); // O [frc] Derivative value
    if(drv_val == 0.0) err_prn(prg_nm_get(),sbr_nm,"Derivative equals zero will cause divide by zero exception");
    sln_dlt=fnc_val/drv_val; // [frc] Adjustment for next guess

    if(dbg_lvl_get() >= dbg_io) std::cout << sbr_nm+"(): itr_idx = " << itr_idx << ", sln = " << *sln << ", fnc_val = " << fnc_val << ", drv_val = " << drv_val << ", sln_dlt = " << sln_dlt << std::endl;

    if(itr_idx != 0 && std::fabs(fnc_val) > std::fabs(fnc_val_old)){
      // Not converging, re-set root guess back half step
      (*sln)+=sln_dlt_old/2.0; // [frc] Current value of root guess
      sln_dlt_old=sln_dlt_old/2.0;  ///[frc] Split adjustment step by two
      std::cout << sbr_nm+"(): New root = " << *sln << ", adjustment = " << sln_dlt_old << std::endl;
    }else{
      // Converging, adjust root guess
      sln_dlt_old=sln_dlt; // [frc] Adjustment of last iteration
      fnc_val_old=fnc_val; // [frc] Function value of last iteration
      sln_old=*sln; // [frc] Current value of root guess
      (*sln)-=sln_dlt; // [frc] Current value of root guess

      // If new guess goes out of range
      if((bnd_LHS-*sln)*(*sln-bnd_RHS) < 0.0){
	std::cout << sbr_nm+"(): New guess sln = " << *sln << " is out of range" << std::endl;
	err_prn(prg_nm_get(),sbr_nm,"Range exceeded");
      } // endif

      // Absolute convergence test
      if(sln_dlt*sgn(sln_dlt) < cnv_thr){
	if(dbg_lvl_get() > dbg_scl) std::cout << sbr_nm+"(): Passed absolute convergence test abs(sln_dlt) < " << cnv_thr << std::endl;
	return rcd;
      } // endif

      // Relative convergence test
      if(*sln == 0.0) err_prn(prg_nm_get(),sbr_nm,"Solution equals zero will cause divide by zero exception in relative convergence test");
      if(std::fabs((*sln) - sln_old)/(*sln) < cnv_thr){
	if(dbg_lvl_get() > dbg_scl) std::cout << sbr_nm+"(): Passed relative convergence test abs((sln-sln_old)/sln) < " << cnv_thr << std::endl;
	return rcd;
      } // endif
    } //endif

  } // end loop over itr
  if(itr_idx == itr_nbr) err_prn(prg_nm_get(),sbr_nm,"Did not converge within itr_nbr = "+nbr2sng(itr_nbr)+" iterations");
  return rcd; // [enm] Return success code
} // end nwt_rph_slvr()

template<class prm_T,class sln_T>int // O [enm] Return success code
eqn_qdr_slvr // [fnc] Solve quadratic equation
(const prm_T prm_a, // I [frc] Parameter a of quadratic equation
 const prm_T prm_b, // I [frc] Parameter b of quadratic equation
 const prm_T prm_c, // I [frc] Parameter c of quadratic equation
 sln_T *sln_1, // O [frc] First root of quadratic equation
 sln_T *sln_2); // O [frc] Second root of quadratic equation

template<class prm_T,class sln_T>int // O [enm] Return success code
eqn_qdr_slvr // [fnc] Solve quadratic equation
(const prm_T prm_a, // I [frc] Parameter a of quadratic equation
 const prm_T prm_b, // I [frc] Parameter b of quadratic equation
 const prm_T prm_c, // I [frc] Parameter c of quadratic equation
 sln_T *sln_1, // O [frc] First root of quadratic equation
 sln_T *sln_2) // O [frc] Second root of quadratic equation
{
  // Purpose: Solve quadratic equation a*x^2+b*x+c=0 with real or complex coefficients
  // Output
  int rcd(0); // O [rcd] Return success code
  // Local
  const std::string sbr_nm("eqn_qdr_slvr"); // [sng] Name of subroutine
  // Main code
  if(prm_a == static_cast<prm_T>(0.0)){
    if(prm_b == static_cast<prm_T>(0.0)){
      err_prn(sbr_nm,"a = b = 0.0 in quadratic equation ax^2 + bx + c = 0");
      rcd=-1;
    } // endif prm_a
    *sln_1=*sln_2=-prm_b/prm_c;
    return rcd;
  } // endif prm_a
  const prm_T qdr_dsc(prm_b*prm_b-static_cast<prm_T>(4.0)*prm_a*prm_c); // [frc] Discriminant of quadratic equation
  const sln_T qdr_dsc_sqrt(std::sqrt(static_cast<sln_T>(qdr_dsc))); // [frc] Square root of discriminant of quadratic equation
  // Pick sign 
  sln_T sgn_dsc_sqrt; // [sgn] Signum of square root of discriminant
  const std::string typ_nm_prm(typeid(prm_T).name());
  const std::string typ_nm_sln(typeid(sln_T).name());
  const std::complex<prc_cmp> prm_b_cpx(prm_b); // CEWI [frc] Complex representation of prm_b
  if(typ_nm_prm.find("complex") != std::string::npos){
    // Complex parameters
    if((conj(prm_b_cpx)*qdr_dsc).real() >= 0.0) sgn_dsc_sqrt=static_cast<sln_T>(1); else sgn_dsc_sqrt=static_cast<sln_T>(-1); // PTV96 p. 178
  }else{
    // Real parameters
    // Make sure solution type is complex if solutions are complex
    if(qdr_dsc_sqrt.imag() != 0.0 && typ_nm_sln.find("complex") == std::string::npos) err_prn(sbr_nm,"Solutions to quadratic equation are complex but calling routine passed reals to hold answers");
    sgn_dsc_sqrt=static_cast<sln_T>(sgn(prm_b));
  } // endif
  const sln_T qdr_fct(static_cast<prm_T>(-0.5)*(prm_b+sgn_dsc_sqrt*qdr_dsc_sqrt)); // [frc] Factor in solution of quadratic equation PTV96 p. 178 (5.6.4)
  assert(qdr_fct != static_cast<sln_T>(0.0));
  *sln_1=prm_c/qdr_fct; // O [frc] First root of quadratic equation
  *sln_2=qdr_fct/prm_a; // O [frc] Second root of quadratic equation

  if(dbg_lvl_get() == dbg_old) std::cout << "Solution(s) to " << prm_a << "*x^2 + " << prm_b << "*x + " << prm_c << " = 0 are x1 = " << *sln_1 << ", x2 = " << *sln_2 << std::endl;

  return rcd;
} // end eqn_qdr_slvr()

template<class prm_T,class sln_T>int // O [enm] Return success code
eqn_cbc_slvr // [fnc] Solve cubic equation
(const prm_T prm_a, // I [frc] Parameter a of cubic equation
 const prm_T prm_b, // I [frc] Parameter b of cubic equation
 const prm_T prm_c, // I [frc] Parameter c of cubic equation
 sln_T *sln_1, // O [frc] First root of cubic equation
 sln_T *sln_2, // O [frc] Second root of cubic equation
 sln_T *sln_3); // O [frc] Third root of cubic equation

template<class prm_T,class sln_T>int // O [enm] Return success code
eqn_cbc_slvr // [fnc] Solve cubic equation
(const prm_T prm_a, // I [frc] Parameter a of cubic equation
 const prm_T prm_b, // I [frc] Parameter b of cubic equation
 const prm_T prm_c, // I [frc] Parameter c of cubic equation
 sln_T *sln_1, // O [frc] First root of cubic equation
 sln_T *sln_2, // O [frc] Second root of cubic equation
 sln_T *sln_3) // O [frc] Third root of cubic equation
{
  /* Purpose: Solve cubic equation x^3+a*x^2+b*x+c=0 with real or complex coefficients
     Parameter type prm_T is std::complex<prc_cmp>, i.e., std::complex<float>, std::complex<double>
     Solution type sln_T is std::complex<prc_cmp>, i.e., std::complex<float>, std::complex<double>
     Usage:
     Sample polynomials
     (z-1)(z-1)(z-1)=0=z^3-3z^2+3z-1
     ccc --dbg=3 --tst=mth --dbg=3 --prm_a=-3 --prm_b=3 --prm_c=-1
     (z-i)(z-i)(z-i)=0=z^3-3iz^2+3z-i
     ccc --dbg=3 --tst=mth --prm_a=(0,-3) --prm_b=(3,0) --prm_c=(0,-1)
     ccc --dbg=3 --tst=mth --prm_a="0-3i" --prm_b="3+0i" --prm_c="0-1i"
     (z-1)(z-2)(z-3)=0=z^3-6z^2+11z-6
     ccc --dbg=3 --tst=mth --prm_a=-6 --prm_b=11 --prm_c=-6
     ccc --dbg=3 --tst=cpx
     ccc --dbg=3 --tst=mth
     fxm: This routine has never worked */
  // Output
  int rcd(0); // O [rcd] Return success code
  // Local
  const std::string typ_nm_sln(typeid(sln_T).name()); // [sng] Type name of solutions
  const std::string sbr_nm("eqn_cbc_slvr"); // [sng] Name of subroutine
  const std::complex<prc_cmp> prm_a_over_3(prm_a/static_cast<prm_T>(3.0)); // [frc] a/3
  // Main code
  if(typ_nm_sln.find("complex") == std::string::npos) err_prn(sbr_nm," solutions to cubic equation must be of a complex type");
  const sln_T fct_Q((prm_a*prm_a-static_cast<prm_T>(3.0)*prm_b)/static_cast<prm_T>(9.0)); // [frc] Factor in solution of cubic equation
  const sln_T fct_R((static_cast<prm_T>(2.0)*prm_a*prm_a*prm_a-PRC_CMP(9.0)*prm_a*prm_b+PRC_CMP(27.0)*prm_c)/PRC_CMP(54.0)); // [frc] Factor in solution of cubic equation
  const sln_T fct_R_sqr(fct_R*fct_R); // [frc] Square of factor in solution of cubic equation
  if(dbg_lvl_get() == dbg_old) std::cout << sbr_nm << "() calling pow(fct_Q,3) with " << fct_Q << " and 3" << std::endl;
  const sln_T fct_Q_cbc(std::pow(fct_Q,static_cast<sln_T>(3))); // [frc] Cube of factor in solution of cubic equation
  if(dbg_lvl_get() == dbg_old) std::cout << sbr_nm << "() called pow(fct_Q,3) with " << fct_Q << " and 3" << std::endl;
  using mth::cst_M_PIl; // (3.1415926535897932384626433832795029L) [frc] 3
  if(fct_Q.imag() == 0.0 && fct_R.imag() == 0.0 && fct_R_sqr.real() < fct_Q_cbc.real()){
    assert(fct_Q_cbc.real() != 0.0);
    const prm_T sln_ngl(std::acos(fct_R.real()/std::sqrt(fct_Q_cbc.real()))); // [rdn] Solution angle
    const sln_T fct_Q_sqrt(-2.0*std::sqrt(fct_Q_cbc.real())); // [frc] Square root of factor
    *sln_1=fct_Q_sqrt*static_cast<sln_T>(std::cos(sln_ngl/PRC_CMP(3.0)))-prm_a_over_3.real(); // [frc] First root of cubic equation
    *sln_2=fct_Q_sqrt*static_cast<sln_T>(std::cos((sln_ngl+static_cast<prm_T>(2.0*mth::cst_M_PIl))/PRC_CMP(3.0)))-static_cast<sln_T>(prm_a_over_3.real()); // [frc] Second root of cubic equation
    *sln_3=fct_Q_sqrt*static_cast<sln_T>(std::cos((sln_ngl-static_cast<prm_T>(2.0*mth::cst_M_PIl))/PRC_CMP(3.0)))-static_cast<sln_T>(prm_a_over_3.real()); // [frc] Third root of cubic equation
  }else{ // endif real roots
    const sln_T cbc_dsc(std::sqrt(fct_R_sqr-fct_Q_cbc)); // [frc] Discriminant of cubic equation PTV96 p. 179 (5.6.13)
    // Pick sign 
    sln_T sgn_dsc_sqrt; // [sgn] Sign of square root of discriminant
    std::complex<prc_cmp> fct_A; // [frc] Factor in solution of cubic equation
    std::complex<prc_cmp> fct_B; // [frc] Factor in solution of cubic equation
    const prm_T one_third(1.0/3.0); // [frc] 1/3
    if((conj(fct_R)*cbc_dsc).real() >= 0.0) sgn_dsc_sqrt=static_cast<sln_T>(1); else sgn_dsc_sqrt=static_cast<sln_T>(-1); // PTV96 p. 178
    if(dbg_lvl_get() == dbg_old) std::cout << sbr_nm << "() calling pow() with " << fct_R+sgn_dsc_sqrt*cbc_dsc << " and " << one_third << std::endl;
    fct_A=-std::pow(fct_R+sgn_dsc_sqrt*cbc_dsc,one_third); // [frc] Factor in solution of cubic equation PTV96 p. 179 (5.6.15)
    if(dbg_lvl_get() == dbg_old) std::cout << sbr_nm << "() called pow() with " << fct_R+sgn_dsc_sqrt*cbc_dsc << " and " << one_third << std::endl;
    if(fct_A == static_cast<std::complex<prc_cmp> >(0.0)) fct_B=static_cast<std::complex<prc_cmp> >(0.0); else fct_B=fct_Q/fct_A; // PTV96 p. 179 (5.6.16)
    const sln_T ApB(fct_A+fct_B); // [frc] A+B
    const sln_T AmB(fct_A-fct_B); // [frc] A-B
    const sln_T i_sqrt_3_over_2(0.0,0.5*std::sqrt(3.0)); // [frc] 0.5*sqrt(-3.0)
    *sln_1=static_cast<sln_T>(-0.5)*ApB-prm_a_over_3+i_sqrt_3_over_2*AmB; // [frc] Second root of cubic equation PTV96 p. 178 (5.6.18)
    // PTV96 p. 178 state that this is the real root when there just a single real root
    *sln_2=ApB-prm_a_over_3; // [frc] First root of cubic equation PTV96 p. 178 (5.6.17)
    *sln_3=static_cast<sln_T>(-0.5)*ApB-prm_a_over_3-i_sqrt_3_over_2*AmB; // [frc] Third root of cubic equation PTV96 p. 178 (5.6.18)
  } // endif some imaginary roots
    
  if(dbg_lvl_get() == dbg_old) std::cout << "Solution(s) to z^3 + " << prm_a << "*z^2 + " << prm_b << "*z + " << prm_c << " = 0 are z1 = " << *sln_1 << ", z2 = " << *sln_2 << ", z3 = " << *sln_3 << std::endl;
    
  return rcd;
} // end eqn_cbc_slvr()

template<class prc_T>prc_cmp // O [frc] Function value
fnc_ptr // [fnc] Function to be solved
(const prc_T abc); // I [frc] Function argument
template<class prc_T>prc_cmp // O [frc] Function value
fnc_ptr // [fnc] Function to be solved
(const prc_T abc) // I [frc] Function argument
{
  // Purpose: Evaluate function to be solved
  return std::sin(abc);
} // end fnc_ptr()

template<class prc_T>prc_cmp // O [frc] Derivative value
drv_ptr // [fnc] Derivative of function to be solved
(const prc_T abc); // I [frc] Derivative argument
template<class prc_T>prc_cmp // O [frc] Derivative value
drv_ptr // [fnc] Derivative of function to be solved
(const prc_T abc) // I [frc] Derivative argument
{
  // Purpose: Evaluate derivative of function to be solved
  return std::sin(abc);
} // end drv_ptr()

#endif // MTH_HH  






