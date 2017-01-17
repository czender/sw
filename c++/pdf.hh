// $Id$ 

// Purpose: Description (definition) of probability density function classes

/* Copyright (C) 1997--2017 Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

// Usage:
// #include <pdf.hh> // Probability density functions

#ifndef PDF_HH // Contents have not yet been inserted in current source file
#define PDF_HH

// C++ headers
#include <iostream> // Standard C++ I/O streams cout, cin
#include <string> // Standard C++ string class
#include <map> // STL map and multimap class templates

// Standard C headers
#include <cmath> // sin cos cos sin 3.14159 
#include <cstdio> // stderr, FILE, NULL, etc. (required by pdf_prn())

// 3rd party vendors
#include <gsl/gsl_sf_gamma.h> // GNU Scientific Library special functions gamma functions

// Personal headers
#include <dbg.hh> // Debugging constants
#include <mth.hh> // Mathematical utilities, constants
#include <utl.hh> // err_prn(), wrn_prn(), Exit_gracefully()
#include <psd.hh> // Particle size distributions
#include <spc_fnc.hh> // Special functions

// Typedefs

// Forward declarations
template<class prc_T> class LogNormalFunction; 
template<class prc_T> class GammaFunction;

/* Terminology of FTV89 p. 17 for lognormal class:
   FTV89 uses two parameters, D_n and sigma, and one variable, D
   FTV89 denotes by D_n the scaling diameter
   FTV89 denotes by sigma what most other sources call ln(gsd), 
   where the other sources say gsd is the geometric standard deviation.
   It is not hard to show that gsd is NOT the geometric standard deviation.
   FTV89 uses the variable D, the diameter in the probability density function.

   My implementation of lognormal distributions:
   uses r, radius, rather than D, diameter, in the PDF
   uses gsd in the same sense SeP97 and PaG77
   uses gsd in the same sense as FTV89 uses exp(sigma)
   uses lngsd in the same sense as FTV89 uses sigma
   uses radius instead of diameter
   labels the median of the lognormal distribution the 
   "number median radius analytic" instead of the "scaling radius"
   Thus, a function from my LogNormal class is created with two parameters:
   rds_nma = r_n, the scaling radius, and gsd = exp(lngsd)
   Some sources, e.g., WKB96, quote ln(gsd), while others, e.g., PaG77, quote gsd */

// Define LogNormalFunction class template
template<class prc_T>
class LogNormalFunction{
public:
  static int nst_nbr_get(); // [nbr] Number of instantiated objects
  // Public member functions
  LogNormalFunction // Default constructor
  (const prc_T& sz_nma_arg=1.0e-6, // [m] Number median size analytic
   const prc_T& gsd_anl_arg=2.0, // [frc] Geometric standard deviation
   const prc_T& cnc_nbr_anl_arg=1.0, // [# m-3] Number concentration analytic
   const std::string &abc_sng_arg="rds"); // [sng] Abscissa type (radius or diameter)

  ~LogNormalFunction(){ // [fnc] Destructor
    nst_nbr--; // [nbr] Number of instantiated class members
  }; // end LogNormalFunction<prc_T>::~LogNormalFunction()

  int // [enm] Return success code
  prs_abc_sng(const std::string &abc_sng); // [fnc] Set abc_typ to rds_typ or dmt_typ

  int // [enm] Return success code
  prm_rst(const psd_cls &psd); // [fnc] Reset parameters

  int // [enm] Return success code
  rds_nma_set(const prc_T rds_nma); // [m] Number median radius analytic
  prc_T& rds_nma_get(){return rds_nma;}; // [m] Number median radius analytic

  int // [enm] Return success code
  dmt_nma_set(const prc_T dmt_nma); // [m] Number median diameter analytic
  prc_T& dmt_nma_get(){return dmt_nma;}; // [m] Number median diameter analytic

  int // [enm] Return success code
  gsd_anl_set(const prc_T gsd_anl); // [frc] Geometric standard deviation
  prc_T& gsd_anl_get(){return gsd_anl;}; // [frc] Geometric standard deviation

  int // [enm] Return success code
  abc_typ_set(const std::string &abc_sng_arg); // [sng] Abscissa type (radius or diameter)
  int& abc_typ_get(){return abc_typ;}; // [enm] Abscissa type (radius or diameter)
  std::string& abc_sng_get(){return abc_sng;}; // [sng] Abscissa type (radius or diameter)

  int // [enm] Return success code
  cnc_nbr_anl_set(const prc_T cnc_nbr_anl); // [# m-3] Number concentration analytic
  prc_T& cnc_nbr_anl_get(){return cnc_nbr_anl;}; // [# m-3] Number concentration analytic

  prc_T& xsa_get(){return xsa;}; // [m2 m-3] Cross-sectional area concentration analytic
  prc_T& sfc_get(){return sfc;}; // [m2 m-3] Surface area concentration analytic
  prc_T& vlm_get(){return vlm;}; // [m3 m-3] Volume concentration analytic
  //  prc_T& mss_get(){return mss;}; // [kg m-3] Mass concentration analytic
  prc_T& rds_nwa_get(){return rds_nwa;}; // [m] Number weighted mean radius analytic
  prc_T& rds_sma_get(){return rds_sma;}; // [m] Surface median radius analytic
  prc_T& rds_vma_get(){return rds_vma;}; // [m] Volume median radius analytic
  prc_T& rds_swa_get(){return rds_swa;}; // [m] Surface area weighted mean radius analytic
  prc_T& rds_vwa_get(){return rds_vwa;}; // [m] Volume weighted mean radius analytic

  prc_T operator()(const prc_T& abc)const; // [fnc] Evaluate distribution function
protected: // fxm: Not sure why I made this protected rather than private
  // Publicly set independent parameters
  int rcd; // [enm] Return success code
  prc_T rds_nma; // [m] Number median radius analytic
  prc_T gsd_anl; // [frc] Geometric standard deviation 
  prc_T cnc_nbr_anl; // # [m3] Total number concentration of aerosol
  prc_T xsa; // [m2 m-3] Cross-sectional area of aerosol per unit volume of air
  prc_T sfc; // [m2 m-3] Surface area of aerosol per unit volume of air
  prc_T vlm; // [m3 m-3] Volume of aerosol per unit volume of air
  //  prc_T mss; // [kg m-3] Mass of aerosol per unit volume of air
  prc_T rds_nwa; // [m] Number weighted mean radius analytic
  prc_T rds_sma; // [m] Surface median radius analytic
  prc_T rds_vma; // [m] Volume median radius analytic
  prc_T rds_swa; // [m] Surface area weighted mean radius analytic
  prc_T rds_vwa; // [m] Volume weighted mean radius analytic

  prc_T dmt_nma; // [m] Number median diameter analytic
  prc_T dmt_nwa; // [m] Number weighted mean diameter analytic
  prc_T dmt_sma; // [m] Surface median diameter analytic
  prc_T dmt_vma; // [m] Volume median diameter analytic
  prc_T dmt_swa; // [m] Surface area weighted mean diameter analytic
  prc_T dmt_vwa; // [m] Volume weighted mean diameter analytic
  // Privately set dependent parameters
  int abc_typ; // [enm] Abscissa type (radius or diameter)
  std::string abc_sng; // [sng] Abscissa type (radius or diameter)
  double lngsd; // [frc] Log of geometric standard deviation ln(gsd)
  double lngsd_sqrt_two_pi_rcp; // [frc] Normalization constant [2*ln(gsd)*sqrt(2*pi)]^{-1}
  double two_lngsd_sqr; // [frc] 2*[ln(gsd)*ln(gsd)]
  double ln_abc_nma; // [frc] Log of number median size analytic ln(abc_nma)
  double abc_nma; // [m] Number median abscissa analytic
  // Private member functions
  int recompute();
private:
  bool rcm_flg; // [flg] Invoke recompute() on set() calls
  // Static private members
  static const int rds_typ; // [enm] Independent variable is radius
  static const int dmt_typ; // [enm] Independent variable is diameter
  static int nst_nbr; // [nbr] Number of instantiated class members
}; // end class LogNormalFunction

// Initialize static member data of template classes _after_ class definition DeD01 p. 717
template<class prc_T>const int LogNormalFunction<prc_T>::rds_typ(0); // [enm] Independent variable is radius
template<class prc_T>const int LogNormalFunction<prc_T>::dmt_typ(1); // [enm] Independent variable is diameter
template<class prc_T>int LogNormalFunction<prc_T>::nst_nbr(0); // [nbr] Number of instantiated class members
// Initialize static member functions 
template<class prc_T>int LogNormalFunction<prc_T>::nst_nbr_get(){return nst_nbr;} // [nbr] Number of instantiated class members

////////////////////////////////////////////////////////////////////////

// Define GammaFunction class template
template<class prc_T>
class GammaFunction{
public:
  // Public member functions
  GammaFunction // Default constructor
  (const prc_T& rds_ffc_arg=50.0, // [m] Effective radius
   const prc_T& var_ffc_arg=0.99, // [frc] Effective variance
   const prc_T& cnc_nbr_anl_arg=1.0); // [# m-3] Number concentration analytic
  
  prc_T& rds_ffc_get(){return rds_ffc;};
  void rds_ffc_set(const prc_T rds_ffc_arg){
    rds_ffc=rds_ffc_arg;
    rcd+=recompute();};// [fnc] Recompute properties of object
  
  prc_T& var_ffc_get(){return var_ffc;};
  void var_ffc_set(const prc_T var_ffc_arg){
    var_ffc=var_ffc_arg;
    rcd+=recompute();}; // [fnc] Recompute properties of object
  
  prc_T& cnc_nbr_anl_get(){return cnc_nbr_anl;};
  void cnc_nbr_anl_set(const prc_T cnc_nbr_anl_arg){
    cnc_nbr_anl=cnc_nbr_anl_arg;
    rcd+=recompute();}; // [fnc] Recompute properties of object
  
  prc_T& xsa_get(){return xsa;};
  prc_T& sfc_get(){return sfc;};
  prc_T& vlm_get(){return vlm;};
  //  prc_T& mss_get();
  
  prc_T operator()(const prc_T& rds)const{
    /* Purpose: Compute and return n(r)=dN(r)/dr
       Reference FTV89
       Return gamma_const*std::pow(rds,(1.0-3.0*b)/b)*std::exp(-rds/(a*b));
       Return cnc_nbr_anl*normalization*std::pow(rds,one_m3bob)*std::exp(-rds/ab);
       Use logarithms to avoid overflow */
    double foo=log_gmm_nrm_cnc_nbr_anl+one_m3bob*std::log(rds)-rds/ab;
    return std::exp(foo);
  }; // end GammaFunction<prc_T>::operator()
  
protected:
  // Publicly set independent parameters
  int rcd; // [enm] Return success code
  prc_T rds_ffc;
  prc_T var_ffc;
  prc_T cnc_nbr_anl; // [# m-3] Total number concentration of aerosol
  prc_T xsa; // [m2 m-3] Cross-sectional area of aerosol per unit volume of air
  prc_T sfc; // [m2 m-3] Surface area of aerosol per unit volume of air
  prc_T vlm; // [m3 m-3] Volume of aerosol per unit volume of air
  //  prc_T mss; // [kg m-3] Mass of aerosol per unit volume of air
  // Privately set dependent parameters
  double ab;
  double two_bm1ob;
  double one_m2bob;
  double one_m3bob;
  double gmm_nrm;
  double log_gmm_nrm_cnc_nbr_anl;
  // Private member functions
  int recompute();
}; // end class GammaFunction

// Class LogNormalFunction Member Function Implementation
template<class prc_T> LogNormalFunction<prc_T>::LogNormalFunction // [fnc] Default constructor
(const prc_T& sz_nma_arg, // [m] Number median size analytic
 const prc_T& gsd_anl_arg, // [frc] Geometric standard deviation
 const prc_T& cnc_nbr_anl_arg, // [# m-3] Number concentration analytic
 const std::string &abc_sng_arg) // [sng] Abscissa type (radius or diameter)
{
  // Purpose: Default constructor for LogNormalFunction<prc_T> objects
  
  std::string sbr_nm("LogNormalFunction<prc_T>::LogNormalFunction"); // [sng] Name of subroutine
  
  // Error status of this object should remain zero
  rcd=0; // [enm] Return success code
  
  // Set rcm_flg to false until essential members have been set
  // set() functions will still perform range-checking and other diagnostics
  // but overhead of numeric computations will be avoided
  rcm_flg=false; // [flg] Invoke recompute() on set() calls
  
  rcd+=prs_abc_sng(abc_sng_arg); // [fnc] Set abc_typ to rds_typ or dmt_typ
  if(abc_typ == rds_typ){
    rcd+=rds_nma_set(sz_nma_arg); // [m] Number median radius analytic
  }else if(abc_typ == dmt_typ){
    rcd+=dmt_nma_set(sz_nma_arg); // [m] Number median diameter analytic
  } // end if
  
  // Use set() functions so range-checking and other diagnostics are performed
  rcd+=gsd_anl_set(gsd_anl_arg); // [frc] Geometric standard deviation
  rcd+=cnc_nbr_anl_set(cnc_nbr_anl_arg); // [# m-3] Number concentration analytic
  
  // Set rcm_flg to true for public access
  rcm_flg=true; // [flg] Invoke recompute() on set() calls
  
  rcd+=recompute(); // [fnc] Recompute properties of object
  
  if(rcd != 0) err_prn(sbr_nm,"Error constructing object");
  
  nst_nbr++; // [nbr] Number of instantiated class members
} // end LogNormalFunction<prc_T>::LogNormalFunction()

template<class prc_T>int // [enm] Return success code
LogNormalFunction<prc_T>::recompute() // [fnc] Recompute properties of object
{
  // Purpose: Recompute properties of lognormal size distribution
  
  const double sqrt_two_pi_rcp=1.0/std::sqrt(2.0*mth::cst_M_PIl);
  
  // Compute parameters which operator() uses in evaluating distribution
  // The parameters depend upon median size, GSD, and number concentration
  lngsd=std::log(gsd_anl); // [frc] Log of geometric standard deviation ln(gsd)
  lngsd_sqrt_two_pi_rcp=sqrt_two_pi_rcp/lngsd; // [frc] Normalization constant [2*ln(gsd)*sqrt(2*pi)]^{-1}
  two_lngsd_sqr=2.0*lngsd*lngsd; // [frc] 2*[ln(gsd)*ln(gsd)]
  
  // abc_typ determines whether operator() returns values per unit radius or per unit diameter:
  if(abc_typ == rds_typ) abc_nma=rds_nma; else abc_nma=dmt_nma; // [frc] Number median size analytic
  ln_abc_nma=std::log(abc_nma); // [frc] Log of number median size analytic ln(abc_nma)
  
  // Compute analytic properties of entire distribution in terms of number median radius analytic
  // Assume spherical particles
  xsa=mth::cst_M_PIl*std::pow(rds_nma,static_cast<prc_T>(2.0))*cnc_nbr_anl*std::exp(2.0*lngsd*lngsd); // [m2 m-3] Cross-sectional area concentration analytic
  sfc=4.0*xsa; // [m2 m-3] Surface area concentration analytic
  vlm=(4.0/3.0)*mth::cst_M_PIl*std::pow(rds_nma,static_cast<prc_T>(3.0))*cnc_nbr_anl*std::exp(4.5*lngsd*lngsd); // [m3 m-3] Volume concentration analytic
  
  rds_nwa=rds_nma*std::exp(1.0*lngsd*lngsd/2.0); // [m] Number weighted mean radius analytic
  rds_sma=rds_nma*std::exp(2.0*lngsd*lngsd); // [m] Surface median radius analytic
  rds_vma=rds_nma*std::exp(3.0*lngsd*lngsd); // [m] Volume median radius analytic
  rds_swa=rds_nma*std::exp(5.0*lngsd*lngsd/2.0); // [m] Surface area weighted mean radius analytic
  rds_vwa=rds_nma*std::exp(7.0*lngsd*lngsd/2.0); // [m] Volume weighted mean radius analytic
  
  dmt_nwa=2.0*rds_nwa; // [m] Number weighted mean diameter analytic
  dmt_sma=2.0*rds_sma; // [m] Surface median diameter analytic
  dmt_vma=2.0*rds_vma; // [m] Volume median diameter analytic
  dmt_swa=2.0*rds_swa; // [m] Surface area weighted mean diameter analytic
  dmt_vwa=2.0*rds_vwa; // [m] Volume weighted mean diameter analytic
  
  //  mss=vlm*dns_prt; // [kg m-3] Mass concentration analytic
  return rcd; // [enm] Return success code
} // end LogNormalFunction<prc_T>::recompute()

template<class prc_T>int LogNormalFunction<prc_T>::rds_nma_set(const prc_T rds_nma_arg){ // [m] Number median radius analytic
  // Purpose: Set size distribution number median radius
  rds_nma=rds_nma_arg; // [m] Number median radius analytic
  if(rds_nma < 0.0 || rds_nma > 1.0) std::cerr << prg_nm_get() << ": ERROR LogNormalFunction<prc_T>::rds_nma_set() reports out of range error with rds_nma = " << rds_nma << std::endl;
  dmt_nma=rds_nma*2.0; // [m] Number median diameter analytic
  if(rcm_flg) rcd+=recompute(); // [fnc] Recompute properties of object
  return rcd; // [enm] Return success code
} // end LogNormalFunction<prc_T>::rds_nma_set()

template<class prc_T>int LogNormalFunction<prc_T>::dmt_nma_set(const prc_T dmt_nma_arg){ // [m] Number median diameter analytic
  // Purpose: Set size distribution number median diameter
  dmt_nma=dmt_nma_arg; // [m] Number median diameter analytic
  if(dmt_nma < 0.0 || dmt_nma > 1.0) std::cerr << prg_nm_get() << ": ERROR LogNormalFunction<prc_T>::dmt_nma_set() reports out of range error with dmt_nma = " << dmt_nma << std::endl;
  rds_nma=dmt_nma/2.0; // [m] Number median radius analytic
  if(rcm_flg) rcd+=recompute(); // [fnc] Recompute properties of object
  return rcd; // [enm] Return success code
} // end LogNormalFunction<prc_T>::dmt_nma_set()

template<class prc_T>int LogNormalFunction<prc_T>::gsd_anl_set(const prc_T gsd_anl_arg){ // [frc] Geometric standard deviation
  gsd_anl=gsd_anl_arg; // [frc] Geometric standard deviation
  if(rcm_flg) rcd+=recompute(); // [fnc] Recompute properties of object
  return rcd; // [enm] Return success code
} // end LogNormalFunction<prc_T>::gsd_anl_set()

template<class prc_T>int LogNormalFunction<prc_T>::abc_typ_set(const std::string &abc_sng_arg){ // [fnc] Set abc_typ to rds_typ or dmt_typ
  rcd+=prs_abc_sng(abc_sng_arg); // [fnc] Set abc_typ to rds_typ or dmt_typ
  return rcd; // [enm] Return success code
} // end LogNormalFunction<prc_T>::abc_typ_set()

template<class prc_T>int LogNormalFunction<prc_T>::cnc_nbr_anl_set(const prc_T cnc_nbr_anl_arg){ // [# m-3] Number concentration analytic
  cnc_nbr_anl=cnc_nbr_anl_arg; // [# m-3] Number concentration analytic
  if(rcm_flg) rcd+=recompute(); // [fnc] Recompute properties of object
  return rcd; // [enm] Return success code
} // end LogNormalFunction<prc_T>::cnc_nbr_anl_set()

template<class prc_T>int LogNormalFunction<prc_T>::prm_rst(const psd_cls &psd){ // Reset parameters
  // Purpose: Reset all parameters of a lognormal distribution
  
  // Set rcm_flg to false until essential members have been set
  rcm_flg=false; // [flg] Invoke recompute() on set() calls
  // NB: Could call either rds_nma_set() or dmt_nma_set() here
  rcd+=rds_nma_set(psd.rds_nma_get()); // [m] Number median radius analytic
  rcd+=gsd_anl_set(psd.gsd_anl_get()); // [frc] Geometric standard deviation
  rcd+=cnc_nbr_anl_set(psd.cnc_nbr_anl_get()); // [# m-3] Number concentration analytic
  // Set rcm_flg to true for public access
  rcm_flg=true; // [flg] Invoke recompute() on set() calls
  
  rcd+=recompute(); // [fnc] Recompute properties of object
  return rcd; // [enm] Return success code
} // end LogNormalFunction<prc_T>::prm_rst()

template<class prc_T>int LogNormalFunction<prc_T>::prs_abc_sng(const std::string &abc_sng_arg){ // [fnc] Set abc_typ to rds_typ or dmt_typ
  /* Purpose: Parse string describing abscissa type and set abc_typ accordingly
     Abscissa type is independent variable to be used by distribution function
     If abc_typ==rds then the distribution computed is per unit radius
     If abc_typ==dmt then the distribution computed is per unit diameter */
  std::string sbr_nm("LogNormalFunction<prc_T>::prs_abc_sng"); // [sng] Name of subroutine
  if(abc_sng_arg.find("dmt") != std::string::npos){
    abc_typ=dmt_typ; // [enm] Abscissa type (radius or diameter)
    abc_sng="dmt"; // [sng] Abscissa type (radius or diameter)
  }else if(abc_sng_arg.find("rds") != std::string::npos){
    abc_typ=rds_typ; // [enm] Abscissa type (radius or diameter)
    abc_sng="rds"; // [sng] Abscissa type (radius or diameter)
  }else{
    err_prn(sbr_nm,"Unknown abc_sng_arg");
    rcd+=1;
  } // endif
  return rcd; // [enm] Return success code
} // end LogNormalFunction<prc_T>::prs_abc_sng()

template<class prc_T>prc_T LogNormalFunction<prc_T>::operator()(const prc_T& abc)const{ // [fnc] Evaluate distribution function
  /* Purpose: Compute and return n(r)=dN(r)/dr or n(D)=dN(D)/dD
     Value of abc_typ determines whether to return n(r) or n(D)
     If abc_typ == rds_typ assume input size argument is radius
     If abc_typ == dmt_typ assume input size argument is diameter
     This is all handled transparently when parameters of size distribution are recomputed
     Thus input argument is named generically abc */
  
  double ln_abc=std::log(abc);
  double ln_abcmln_abc_nma_sqr=(ln_abc-ln_abc_nma)*(ln_abc-ln_abc_nma);
  //  double ln_r_over_rg=std::log(abc/abc_nma);
  //  double ln_r_over_rg_sqr=lnln_r_over_rg*ln_r_over_rg;
  
  return cnc_nbr_anl*lngsd_sqrt_two_pi_rcp*std::exp(-ln_abcmln_abc_nma_sqr/two_lngsd_sqr)/abc;
} // end LogNormalFunction<prc_T>::operator()

// Class GammaFunction Member Function Implementation
template<class prc_T> GammaFunction<prc_T>::GammaFunction
(const prc_T& rds_ffc_arg, // [m] Effective radius
 const prc_T& var_ffc_arg, // [frc] Effective variance
 const prc_T& cnc_nbr_anl_arg) // [# m-3] Number concentration analytic
{
  rcd=0; // [enm] Return success code
  rds_ffc=rds_ffc_arg; // [m] Effective radius
  var_ffc=var_ffc_arg; // [frc] Effective variance
  cnc_nbr_anl=cnc_nbr_anl_arg; // [# m-3] Number concentration analytic
  rcd+=recompute(); // [fnc] Recompute properties of object
} // end GammaFunction<prc_T> constructor

template<class prc_T>int // [enm] Return success code
GammaFunction<prc_T>::recompute() // [fnc] Recompute properties of object
{
  const std::string sbr_nm("GammaFunction<prc_T>::recompute"); // [sng] Name of subroutine
  
  ab=rds_ffc*var_ffc;
  two_bm1ob=(2.0*var_ffc-1)/var_ffc;
  one_m2bob=(1.0-2.0*var_ffc)/var_ffc;
  if(one_m2bob == mth_fnc::round(one_m2bob) && one_m2bob <= 0){
    wrn_prn(sbr_nm,"one_m2bob = "+nbr2sng(one_m2bob)+" is a negative integer or zero.\nrds_ffc = "+nbr2sng(rds_ffc)+", var_ffc = "+nbr2sng(var_ffc)+", (1.0-2.0*var_ffc)/var_ffc = "+nbr2sng(one_m2bob)+"\nThe ln(gamma) function is undefined for negative integer arguments, so the program will crash.\nHINT: Choose var_ffc such that (1.0-2.0*var_ffc)/var_ffc is not a negative integer or zero");
  } // endif
  one_m3bob=(1.0-3.0*var_ffc)/var_ffc;
  /* Compute logarithm of normalization to avoid overflows
     Argument (one_m2bob) must not be a negative integer
     Deprecated Numerical Recipes solution:
     double log_gmm_nrm=two_bm1ob*std::log(ab)-nr_gammln(one_m2bob); */
  double log_gmm_nrm=two_bm1ob*std::log(ab)-gsl_sf_lngamma(one_m2bob);
  gmm_nrm=std::exp(log_gmm_nrm);
  //  gmm_nrm=std::pow(ab,two_bm1ob)/gamma(one_m2bob);
  log_gmm_nrm_cnc_nbr_anl=std::log(cnc_nbr_anl)+log_gmm_nrm;
  
  xsa=mth::cst_M_PIl*std::pow(rds_ffc,static_cast<prc_T>(2.0))*(1.0-var_ffc)*(1.0-2.0*var_ffc)*cnc_nbr_anl;
  sfc=4.0*xsa; // Assumes spherical particles
  vlm=(4.0/3.0)*mth::cst_M_PIl*std::pow(rds_ffc,static_cast<prc_T>(3.0))*(1.0-var_ffc)*(1.0-2.0*var_ffc)*cnc_nbr_anl;
  //  mss=vlm*dns_prt;
  return rcd; // [enm] Return success code
} // end GammaFunction<prc_T>::recompute()

////// Prototype functions with C++ linkages

// Define inline'd functions in header so source is visible to calling files

template<class prc_T>prc_T swa2nma(const prc_T sz_swa,const prc_T gsd); // O [fnc] Convert surface area weighted mean size to number median size
template<class prc_T>prc_T // O [m] Number median particle size analytic
swa2nma // O [fnc] Convert surface area weighted mean size to number median size
(const prc_T sz_swa, // I [m] Surface area weighted mean size analytic
 const prc_T gsd) // I [frc] Geometric standard deviation
{
  /* Purpose: Given the surface area weighted mean particle size and geometric standard deviation
     for a lognormal size distribution, compute and return the number median size of each distribution
     Definition of gsd conforms to SeP97 */
  prc_T ln_gsd(std::log(gsd)); // [frc] Log of geometric standard deviation
  return sz_swa*std::exp(-2.5*ln_gsd*ln_gsd);
} // end swa2nma()

template<class prc_T>prc_T vma2nma(const prc_T sz_vma,const prc_T gsd); // O [fnc] Convert volume median size to number median size
template<class prc_T>prc_T // O [m] Number median particle size analytic
vma2nma // O [fnc] Convert volume median size to number median size
(const prc_T sz_vma, // I [m] Volume median size analytic
 const prc_T gsd) // I [frc] Geometric standard deviation
{
  /* Purpose: Given the volume median particle size and geometric standard deviation
     for a lognormal size distribution, compute and return the number median size of each distribution
     Definition of gsd conforms to SeP97 */
  prc_T ln_gsd(std::log(gsd)); // [frc] Log of geometric standard deviation
  return sz_vma*std::exp(-3.0*ln_gsd*ln_gsd);
} // end vma2nma()

#endif // PDF_HH  
