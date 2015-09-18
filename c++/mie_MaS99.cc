// $Id$ 

// Purpose: Implementation (declaration) of Markel and Shalaev (1999) Mie scattering solutions

/* Copyright (C) 2005--2014 Charlie Zender, Jorge Talamantes, Vadim Markel, and Vladimir Shalaev
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

/* History:
   200412--200702: Jorge Talamantes translated Markel's routines from Fortran77 to C++
   200502--200703: Charlie Zender merged Jorge's codes into libcsm_c++ 
   200509--200511: Jorge Talamantes updated code to fix Markel's programming errors and to handle arbitrary D. Code uses quadrature for D!=3, and exact solution for D=3 (D = fractal dimesionality of soot aggregates)
   200511--      : Charlie merges update code into libcsm_c++
   
   Reference articles: MaS99, MaS00, Mar02
   Markel, V. A., V. M. Shalaev, Absorption of light by soot particles in micro-droplets of water, J. Quant. Spectrosc. Radiat. Transfer, 63(2-6), pp. 321-339, 1999.
   
   Markel, V. A., V. M. Shalaev, Erratum on "Absorption of light by soot particles in micro-droplets of water", J. Quant. Spectrosc. Radiat. Transfer, 66(6), p. 591, 2000.
   
   Markel, V. A., The effects of averaging on the enhancement factor for absorption of light by carbon particles in microdroplets of water, J. Quant. Spectrosc. Radiat. Transfer, 72(6), pp. 765-774, 2002. 
   
   Markel, V. A. Erratum to "The effects of averaging on the enhancement factor for absorption of light by carbon particles in microdroplets of water", to appear in J. Quant. Spectrosc. Radiat. Transfer */

#include <mie_MaS99.hh> // Mie scattering solutions from MaS99

// structure definition to pass parameters to gsl_integration_XXX
struct M{
  double dmn_frc;
  int alpha;
  int n;
}; // end struct M

/* function defined to compute the values of the integrand in Markel's
   paper I Erratum (JQSRT 66 (2000) p. 591):
   I_markel == x^(dmn_frc-alpha)*[j_n(x)]^2
   (implementation for dmn_frc != 3) */
double I_markel
(double x, void*p)
{
  struct M * params=(struct M *)p;
  double dmn_frc=(params->dmn_frc);
  int alpha=(params->alpha);
  int n=(params->n);
  
  double N=0.5+static_cast<double>(n);
  double j;
  
  if ((N-x) > 70.0){
    j=0.0; // Prevent underflow error in Bessel function
  }else{
    j=gsl_sf_bessel_Jnu(N,x);
  } // endif
  
  double f=std::pow(x,dmn_frc-static_cast<double>(alpha)-1)*j*j;
  return f;
} // end I_markel()

double QAGP /* [fnc] Adaptive stepsize integration with known singular points */
(double x1,
 double dmn_frc,
 int alpha,
 int n)
{
  /* Purpose: Implementation of GSL's adaptive integration with known singular points 
     (dmn_frc != 3) */
  
  gsl_integration_workspace *w=gsl_integration_workspace_alloc(100000);
  
  double result,error;
  using mth::cst_M_PIl; // (3.1415926535897932384626433832795029L) [frc] 3
  double pi=mth::cst_M_PIl;
  double max_rel_error=1.0e-12;
  
  struct M VM;
  VM.dmn_frc=dmn_frc;
  VM.alpha=alpha;
  VM.n=n;
  
  gsl_function F;
  F.function=&I_markel;
  F.params=&VM;
  
  double pts[2];
  pts[0]=0.0;
  pts[1]=x1;
  
  size_t npts=2;
  
  gsl_integration_qagp(&F,pts,npts,0,max_rel_error,1000,w,&result,&error); 
  
  return (pi/2.0)*result;
} // end QAGP()

double I_i_1
(const double sz_prm_idx_rfr,
 const int i)
{
  // Purpose: I_markel == x^(dmn_frc-alpha) * [j_n(x)]^2 with dmn_frc = 3, alpha = 1
  using mth::cst_M_PIl; // (3.1415926535897932384626433832795029L) [frc] 3
  double pi=mth::cst_M_PIl;
  double N=0.5+static_cast<double>(i);
  double j=gsl_sf_bessel_Jnu(N,sz_prm_idx_rfr);
  double dj=0.5*( gsl_sf_bessel_Jnu(N-1.0,sz_prm_idx_rfr)-gsl_sf_bessel_Jnu(N+1.0,sz_prm_idx_rfr));
  double x1sqr=sz_prm_idx_rfr *sz_prm_idx_rfr;
  double result=0.5*x1sqr*dj*dj+0.5*(x1sqr-N*N)*j*j;
  return (pi/2.0)*result;
} // end I_i_1()

double // O [fnc] Absorption enhancement for inclusions in weakly-absorbing spheres
mie_sph_abs_fct_MaS99 // [fnc] Absorption enhancement for inclusions in weakly-absorbing spheres
(const double sz_prm, // I [m m-1] Size parameter
 const double idx_rfr_mdm_rl, // I [frc] Refractive index of weakly-absorbing sphere, real component
 const double dmn_frc) // I [frc] Fractal dimensionality of inclusions
{
  /* Purpose: Vadim Markel's method to compute an enhancement factor (G) 
     of energy absorption by absorbing spheres in weakly-absorbing media.
     This function was initially largely a direct translation of Markel's f77 code
     into C++, with modifications to improve convergence. That f77 code only
     did dmn_frc = 3 (a very special case from the point of view of computation).
     This new version has been corrected to account for the error in VM's
     second Erratum; also, this version has been changed to do dmn_frc=1.8
     
     Program tested against matlab results for x=1,1000, delta x=1
     
     Input:
     size parameter sz_prm = 2*pi*a/lambda
     idx_rfr_mdm_rl = real part of index of refraction of water
     (imaginary part is taken to be zero)
     
     Notes:
     e - small constant. Summation over n stops when both
     coefficients c_n^2 and d_n^2 become less than e */
  
  const std::string sbr_nm("mie_sph_abs_fct_MaS99"); // [sng] Subroutine name
  
  double bj[5001],by[5001];
  double dbj[5001],dby[5001];
  double bj1[5001],by1[5001];
  double dbj1[5001],dby1[5001];
  bool   flag_cc,flag_dd;
  
  const double pi(4.0*std::atan(1.0));
  const double eps_max(1.0e-12); // [frc] Relative accuracy for convergence
  const double coeff=sqrt(pi/2.0);
  
  const double idx_rfr_mdm_rl_sqr(idx_rfr_mdm_rl*idx_rfr_mdm_rl);
  const double idx_rfr_mdm_rl_sqr_m1(idx_rfr_mdm_rl_sqr-1.0);
  
  const double sz_prm_idx_rfr(idx_rfr_mdm_rl*sz_prm);
  
  const int trm_nbr(3); // Markel had trm_nbr = xnfact*sz_prm
  
  const bool R3=std::abs(dmn_frc-3.0) < 0.001; // [flg] True if dmn_frc=3
  
  /* Compute first trm_nbr Bessel functions plus the 0th and 
     trm_nbr+first Bessel functions (needed for the derivatives) */
  for (int i=0;i<=trm_nbr+1;i++){
    bj[i] =gsl_sf_bessel_Jnu(i+0.5,sz_prm);
    by[i] =gsl_sf_bessel_Ynu(i+0.5,sz_prm);
    bj1[i]=gsl_sf_bessel_Jnu(i+0.5,sz_prm_idx_rfr);
    by1[i]=gsl_sf_bessel_Ynu(i+0.5,sz_prm_idx_rfr);
  } // end loop over i
  
  for(int i=0;i<=trm_nbr+1;i++){
    bj[i] =bj[i] *coeff/sqrt(sz_prm);
    bj1[i]=bj1[i]*coeff/sqrt(sz_prm_idx_rfr);
    by[i] =by[i] *coeff/sqrt(sz_prm);
    by1[i]=by1[i]*coeff/sqrt(sz_prm_idx_rfr);
  } // end loop over i
  
  // Compute first trm_nbr derivatives of the Bessel functions 
  for(int i=1;i<=trm_nbr;i++){
    dbj[i] =(i*bj [i-1]-(i+1)*bj [i+1])/(2*i+1);
    dbj1[i]=(i*bj1[i-1]-(i+1)*bj1[i+1])/(2*i+1);
    dby[i] =(i*by [i-1]-(i+1)*by [i+1])/(2*i+1);
    dby1[i]=(i*by1[i-1]-(i+1)*by1[i+1])/(2*i+1);
  } // end loop over i
  
  // Initialize enhancement absorption factor (G in Markel's papers)
  double abs_fct_MaS99=0.0;
  
  // Initialize booleans to help decide convergence
  flag_cc=false;
  flag_dd=false;
  
  // Compute first trm_nbr terms in series for G
  for(int i=1;i<=trm_nbr;i++){
    
    double xi=static_cast<double>(i);
    
    // a1, ..., a5 are intermediate variables
    double a1=bj[i]*dby[i]-by[i]*dbj[i];
    double a2=bj1[i]*dbj[i]-idx_rfr_mdm_rl*bj[i]*dbj1[i];// VM's first error fixed
    double a3=bj1[i]*dby[i]-idx_rfr_mdm_rl*dbj1[i]*by[i];
    double a4=idx_rfr_mdm_rl_sqr_m1*bj[i]*bj1[i]
     +sz_prm_idx_rfr*(idx_rfr_mdm_rl*bj1[i]*dbj[i]-bj[i]*dbj1[i]);
    double a5=idx_rfr_mdm_rl_sqr_m1*bj1[i]*by[i]
     +sz_prm_idx_rfr*(idx_rfr_mdm_rl*bj1[i]*dby[i]-by[i]*dbj1[i]);
    
    // cc and dd are |c_i|^2 and |d_i|^2 respectively 
    double cc=a1*a1/(a2*a2+a3*a3);
    double dd=idx_rfr_mdm_rl_sqr*sz_prm*sz_prm*a1*a1/(a4*a4+a5*a5);
    
    /* Test for convergence---done if both |c_i|^2 < eps_max and |d_i|^2 < eps_max */
    if(cc < eps_max) flag_cc=true;
    if(dd < eps_max) flag_dd=true;
    if(flag_cc && flag_dd) goto label3;
    
    /* If convergence has not been achieved, but either
       |c_i|^2< eps_max or |d_i|^2 < eps_max, 
       then re-set that one to zero */
    double ccc=cc;
    double ddd=dd;
    if(flag_cc) ccc=0.0;
    if(flag_dd) ddd=0.0;
    
    // Compute integrals \int_0^x1 x^(dmn_frc-alpha) [j_n(x)]^2, alpha=1, 3
    double v,w;
    if(R3){
      v=I_i_1(sz_prm_idx_rfr,i);
      w=0;
    }else{ // !R3
      v=QAGP(sz_prm_idx_rfr,dmn_frc,1,i);
      w=QAGP(sz_prm_idx_rfr,dmn_frc,3,i);
    } // !R3
    
    // Compute G value after i terms 
    abs_fct_MaS99=abs_fct_MaS99+(2.*xi+1.0)
      *(ccc*v+ddd*((5.0-dmn_frc)*std::pow(sz_prm_idx_rfr,dmn_frc-2.0)*bj1[i]*bj1[i]/2.0
		   +std::pow(sz_prm_idx_rfr,dmn_frc-1.0)*dbj1[i]*bj1[i]+v
		   +w*(4.0-dmn_frc)*(3.0-dmn_frc)/2.0));
  }
  
  /* Pursue convergence one term at a time. Patch begins
     The following section is essentially a repeat of the previous one.
     Same comments apply as above */
  
  for(int k=1;k<=(5000-(trm_nbr+1));k++){
    
    int i=(trm_nbr+1)+k;
    bj [i]=gsl_sf_bessel_Jnu(i+0.5,sz_prm);
    by [i]=gsl_sf_bessel_Ynu(i+0.5,sz_prm);
    bj1[i]=gsl_sf_bessel_Jnu(i+0.5,sz_prm_idx_rfr);
    by1[i]=gsl_sf_bessel_Ynu(i+0.5,sz_prm_idx_rfr);
    
    bj[i] =bj[i] *coeff/sqrt(sz_prm);
    bj1[i]=bj1[i]*coeff/sqrt(sz_prm_idx_rfr);
    by[i] =by[i] *coeff/sqrt(sz_prm);
    by1[i]=by1[i]*coeff/sqrt(sz_prm_idx_rfr);
    
    i=trm_nbr+k;
    dbj[i] =(i*bj [i-1]-(i+1)*bj [i+1])/(2*i+1);
    dbj1[i]=(i*bj1[i-1]-(i+1)*bj1[i+1])/(2*i+1);
    dby[i] =(i*by [i-1]-(i+1)*by [i+1])/(2*i+1);
    dby1[i]=(i*by1[i-1]-(i+1)*by1[i+1])/(2*i+1);
    
    flag_cc=false;
    flag_dd=false;
    
    double xi=static_cast<double>(i);
    
    double a1=bj[i]*dby[i]-by[i]*dbj[i];
    double a2=bj1[i]*dbj[i]-idx_rfr_mdm_rl*bj[i]*dbj1[i];//VM's first error fixd
    double a3=bj1[i]*dby[i]-idx_rfr_mdm_rl*dbj1[i]*by[i];
    double a4=idx_rfr_mdm_rl_sqr_m1*bj[i]*bj1[i]
     +sz_prm_idx_rfr*(idx_rfr_mdm_rl*bj1[i]*dbj[i]-bj[i]*dbj1[i]);
    double a5=idx_rfr_mdm_rl_sqr_m1*bj1[i]*by[i]
     +sz_prm_idx_rfr*(idx_rfr_mdm_rl*bj1[i]*dby[i]-by[i]*dbj1[i]);
    
    double cc=a1*a1/(a2*a2+a3*a3);
    double dd=idx_rfr_mdm_rl_sqr*sz_prm*sz_prm*a1*a1/(a4*a4+a5*a5);
    if(cc < eps_max) flag_cc=true;
    if(dd < eps_max) flag_dd=true;
    if(flag_cc && flag_dd) goto label3;
    double ccc=cc;
    double ddd=dd;
    if(flag_cc) ccc=0.0;
    if(flag_dd) ddd=0.0;
    
    // Compute integrals \int_0^x1 x^(dmn_frc-alpha) [j_n(x)]^2, alpha=1, 3
    double v,w;
    if(R3){
      v=I_i_1(sz_prm_idx_rfr,i);
      w=0;
    }else{ // !R3
      v=QAGP(sz_prm_idx_rfr,dmn_frc,1,i);
      w=QAGP(sz_prm_idx_rfr,dmn_frc,3,i);
    } // !R3
    
    // Compute value of G after i terms 
    abs_fct_MaS99=abs_fct_MaS99+(2.0*xi+1.0)*
      (ccc*v+ddd*((5.0-dmn_frc)*std::pow(sz_prm_idx_rfr,dmn_frc-2.0)*bj1[i]*bj1[i]/2.0
		  +std::pow(sz_prm_idx_rfr,dmn_frc-1.0)*dbj1[i]*bj1[i]+v
		  +w*(4.0-dmn_frc)*(3.0-dmn_frc)/2.0));
    
  } // end of patch
  
  wrn_prn(sbr_nm,"convergence failed after 5000 iterations for size parameter sz_prm = "+nbr2sng(sz_prm)+" and real refractive index of medium = "+nbr2sng(idx_rfr_mdm_rl));
  
 label3:
  abs_fct_MaS99=abs_fct_MaS99*dmn_frc/(2.0*std::pow(sz_prm_idx_rfr,dmn_frc)); // VM's second error fixed
  if (abs_fct_MaS99 < 0.0) abs_fct_MaS99=0.0; 
  
  return abs_fct_MaS99;
  
} // end mie_sph_abs_fct_MaS99()

int // O [enm] Return success code
mie_sph_abs_fct_MaS99_tst // [fnc] Test mie_sph_abs_fct_MaS99()
(void)
{
  /* Purpose: Test mie_sph_abs_fct_MaS99() */
  
  double idx_rfr_mdm_rl=1.33;
  double dmn_frc=3.0;

  int rcd(0); // [enm] Return success code
  
  for (int l=1;l<=10;l++){
    double sz_prm=static_cast<double>(l);
    double abs_fct_MaS99=mie_sph_abs_fct_MaS99(sz_prm,idx_rfr_mdm_rl,dmn_frc);
    std::cout << sz_prm << "  " << abs_fct_MaS99 << std::endl;
  } // end for
  
  return rcd;
} // end mie_sph_abs_fct_MaS99()
