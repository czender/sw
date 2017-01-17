// $Id$ 

// Purpose: Implementation (declaration) of water refractive indices

/* Copyright (C) 1997--2017 Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

/* Source:
   Converted refraction index routines (from f77 to C++)
   idx_rfr_H2O_lqd_bpb.f and idx_rfr_H2O_ice_sgw.f to idx_rfr_H2O.cc
   D. Newman Nov 2001
   Extensive Cleanup:
   C. Zender Sep 2002
   C. Zender Feb 2005
   C. Zender Jan 2008 add Kou et al. AO (1993) KLC93 ice refractive indices
   C. Zender Mar 2008 add Gosse et al. AO (1995) GLC95 ice refractive indices in comments
   C. Zender Aug 2008 notes that WaB08 adopt GLC95 for most of NIR: 1.4--2.9 um and 3.4--7.8 um. WaB08 uses War84 which used Schaaf and Williams for 2.9--3.4 um
   
   Notes:
   1. All one-dimensional arrays (except tempRef) have Fortran (1-based) indexing
   Never access element zero, x[0]
   2. Two-dimensional arrays tabReT and tabImT are dimensioned
   tabReT[ntemp][nwlT1].  
   Access first index starting at 0, second index starting at 1  
   Also remember reversal of indices from fortran to C++
   3. Can test with idx_rfr_H2O_tst() */  

#include <idx_rfr_H2O.hh>

// Namespaces
using namespace idx_rfr_H2O; // [nms] Water refractive index namespace

std::complex<prc_cmp> // O [frc] Refractive index of liquid water
RayWater // [fnc] Refractive index of liquid water in microwave regime
(const prc_cmp wvl_ctr, // I [m] Wavelength at band center
 const prc_cmp tpt) // I [K] Temperature
{
  /* Purpose: Calculates complex refractive index of pure liquid water for
     wavelengths between 2.0 microns and 10.0 m.
     
     Input:  wvl = wavelength (2 microns to 1.0e7 microns)
     tpt = temperature [K]
     
     Output: raywat = complex refractive index
     (with positive imaginary part)
     
     Method : Ray's analytic fits based on some theories of Debye
     
     Reference : Ray, p., 1972:  Broadband complex refractive indices
     of ice and water, Appl. Opt. 11, 1836--1844
     
     Specifications of local variables
     
     epsils :  static dielectric constant
     (epsilon-sub-s, Ray eq. 4)
     epsinf :  high-frequency dielectric const
     (epsilon-sub-s, Ray eq. 7a)
     epsre, :  real and imaginary parts of dielectric constant
     epsim      (Ray eqs. 5,6)
     pi     :  3.14159...
     ab     :  summation terms, Ray eq. 8
     arre   :  correction to real index for wvl < 6 micron
     brre   :  correction to imag. index
     rre    :  real part of refractive index
     rim    :  imag part of refractive index
     tpt_cls     :  centigrade temperature
     term   :  summation terms, Ray eq. 9
     tbarp1 :  tbar+ 1, coefficient in Ray eq. 9
     wvlcm   :  wavelength in cm
     wvldeby :  wavelength above which to apply debye theory
     wvl_min_mcr  :  minimum tabulated wavelength
     wvl_max_mcr  :  maximum tabulated wavelength
     wvls    :  relaxation wavelength in cm
     (lambda-sub-s, Ray eq. 7c)
     alpha  :  temporary storage variables */
  
  const std::string sbr_nm("RayWater"); // [sng] Subroutine name
  prc_cmp alpha,epsils,epsinf,rim,rre,tpt_cls,wvlcm,wvls;
  std::complex<prc_cmp> ci(0.0,1.0),ceps,raywat;
  const prc_cmp wvl_min_mcr=1.0,wvl_max_mcr=1.0e+8;
  const prc_cmp wvl_ctr_mcr(wvl_ctr*1.0e6); // [m]->[um] Wavelength at band center
  
  if(wvl_ctr_mcr < wvl_min_mcr || wvl_ctr_mcr > wvl_max_mcr) std::cout << "RayWater--wavelength out of table range" << std::endl;
  
  // Use Ray's fits to debye theory expressions
  tpt_cls=tpt-273;
  epsinf=5.27137+tpt_cls*(0.0216474-0.00131198*tpt_cls);
  epsils=78.54*(1.0+(tpt_cls-25.0)*(-4.579e-3+(tpt_cls-25.0)*(1.19e-5-2.8e-8*(tpt_cls-25.0))));
  
  // Ray eq. 7b
  alpha=-16.8129/tpt+0.0609265;
  wvlcm=1.0e-4*wvl_ctr_mcr;
  wvls =3.3836e-4*std::exp(2514.0/tpt);
  using mth::cst_M_PIl; // (3.1415926535897932384626433832795029L) [frc] 3
  ceps=epsinf+(epsils-epsinf)/
    (PRC_CMP(1.0)+
     (static_cast<prc_cmp>(std::cos(0.5*mth::cst_M_PIl*(1.0-alpha)))+
      ci*static_cast<prc_cmp>(std::sin(0.5*mth::cst_M_PIl*(1.0-alpha))))*
     static_cast<prc_cmp>(std::pow(wvls/wvlcm,PRC_CMP(1.0)-alpha)))
    -PRC_CMP(0.00666667)*ci*wvlcm;
  
  // Complex refractive index from cole-cole extension to debye theory
  raywat=std::sqrt(ceps);
  rim=-raywat.imag();
  rre=raywat.real();
  
  // Corrections to imaginary index to account for absorption bands
  // (Ray eq. 8 + table 2)
  if(wvl_ctr_mcr < 3000 && wvl_ctr_mcr > 300)
    rim=rim+ab(wvl_ctr_mcr,0.25,300.0,0.47,3.0)
      +ab(wvl_ctr_mcr,0.41,62.0,0.35,1.7)
      +ab(wvl_ctr_mcr,0.39,17.0,0.45,1.3);
  else if(wvl_ctr_mcr <= 300 && wvl_ctr_mcr > 62)
    rim=rim+ab(wvl_ctr_mcr,0.25,300.0,0.40,2.0)
      +ab(wvl_ctr_mcr,0.41,62.0,0.35,1.7)
      +ab(wvl_ctr_mcr,0.39,17.0,0.45,1.3);
  else if(wvl_ctr_mcr <= 62 && wvl_ctr_mcr > 17)
    rim=rim+ab(wvl_ctr_mcr,0.25,300.0,0.40,2.0)
      +ab(wvl_ctr_mcr,0.41,62.0,0.22,1.8)
      +ab(wvl_ctr_mcr,0.39,17.0,0.45,1.3);
  else if(wvl_ctr_mcr <= 17 && wvl_ctr_mcr > 6.1)
    rim=rim+ab(wvl_ctr_mcr,0.12,6.1,0.042,0.6)
      +ab(wvl_ctr_mcr,0.39,17.0,0.165,2.4)
      +ab(wvl_ctr_mcr,0.41,62.0,0.22,1.8);
  else if(wvl_ctr_mcr <= 6.1 && wvl_ctr_mcr > 4.95)
    rim=rim+ab(wvl_ctr_mcr,0.12,6.1,0.009,2.0)
      +ab(wvl_ctr_mcr,0.01,4.95,0.05,1.0);
  else if(wvl_ctr_mcr <= 4.95 && wvl_ctr_mcr > 2.95)
    rim=rim+ab(wvl_ctr_mcr,0.27,2.97,0.04,2.0)
      +ab(wvl_ctr_mcr,0.01,4.95,0.06,1.0);
  else if(wvl_ctr_mcr <= 2.95)
    rim=rim+ab(wvl_ctr_mcr,0.27,2.97,0.025,2.0)
      +ab(wvl_ctr_mcr,0.01,4.95,0.06,1.0);
  else
    wrn_prn(prg_nm_get(),sbr_nm," No correction available for imaginary index of refraction of liquid H2O at this wavelength = "+nbr2sng(wvl_ctr_mcr)+" um.");
  
  // Corrections to real index
  // (Ray eq. 9 + table 3)
  if(wvl_ctr_mcr < 1000 && wvl_ctr_mcr > 340)
    rre=rre*((wvl_ctr_mcr-340.0)/660.0)+brre(tpt_cls,wvl_ctr_mcr)*((1000.0-wvl_ctr_mcr)/660.0);
  else if(wvl_ctr_mcr <= 340 && wvl_ctr_mcr > 7)
    rre=brre(tpt_cls,wvl_ctr_mcr);
  else if(wvl_ctr_mcr <= 7 && wvl_ctr_mcr > 6)
    rre=arre(tpt_cls,wvl_ctr_mcr)*(7.0-wvl_ctr_mcr)+brre(tpt_cls,wvl_ctr_mcr)*(wvl_ctr_mcr-6.0);
  else if(wvl_ctr_mcr <= 6) 
    rre=arre(tpt_cls,wvl_ctr_mcr);
  else
    wrn_prn(prg_nm_get(),sbr_nm," No correction available for real index of refraction of liquid H2O at this wavelenth = "+nbr2sng(wvl_ctr_mcr)+" um.");
  
  raywat=std::complex<prc_cmp>(rre,rim);
  return raywat;
} // end Raywater()

int // O [enm] Return success code
idx_rfr_H2O_lqd_get // [fnc] Refractive index of liquid water
(prc_cmp wvl_ctr, // I [m] Wavelength at band center
 prc_cmp tpt, // I [K] Temperature
 std::complex<prc_cmp> *idx_rfr) // O [frc] Refractive index of liquid water
{
  /* Purpose: Refractive index of liquid water
     Author: Bruce P. Briegleb, NCAR
     Modifications: Charlie Zender, NCAR
     
     Calculates complex refractive index of pure liquid water for
     wavelengths between 0.01 microns and 10.0 m.  Temperature
     dependence is only considered for wavelengths greater
     than 10.0 microns.
     
     I N P U T :  wvl_ctr = wavelength (microns)
     tpt = temperature [K] -- for wvl_ctr > 10 only
     
     O U T P U T :  REFWAT = complex refractive index
     (with positive-definite imaginary part)
     
     METHOD :
     
     for wavelengths less then 10.0 microns :
     tabular interpolation assuming real index and 
     log(imaginary index) linear in log(wavelength).
     
     for wavelengths 10.0-20.0 microns :
     weighted data correction using Ray's model to account
     for temperature dependence
     
     for wavelengths 20.0-1.0e7 microns :
     data correction using Ray's model to account for 
     temperature dependence
     
     for wavelengths greater then 1.0e7 microns :
     Ray's analytic fits based on some theories of Debye
     
     REFERENCES :
     
     (1) for 0.01-1.0e7 microns:  Segelstein, D., 1981:
     "The Complex Refractive Index of Water", M.S. Thesis,
     University of Missouri--Kansas City
     
     (2) for 10.0-1.0e7 microns:  Ray, P., 1972:  Broadband Complex
     Refractive Indices of Ice and Water, Appl. Opt. 11,
     1836-1844
     
     (There is a new reference, WWQ89, Wieliczka, D. et al., Appl. Opt.
     28, 1714-1719, 1989, with some updated data for the IR)
     ---------------------------------------------------------------------
     
     *** specifications of local variables
     
     frac   :  quantity between 0 and 1 used for interpolation
     rre    :  real part of refractive index
     rim    :  imag part of refractive index
     t3im   :  temporay imag part of refractive index at 300k
     t3re   :  temporary real part of refractive index at 300k
     corre  :  correction to real index
     corim  :  correction to imag. index
     wd     :  weight factor for data points
     wr     :  weight factor for ray model
     wvl_min_mcr  :  minimum tabulated wavelength
     wvl_max_mcr  :  maximum tabulated wavelength */
  
  int rcd(0); // [enm] Return success code
  int i;
  prc_cmp corim,corre,frac,rim=CEWI_cpv,rimd,rre=CEWI_cpv,rred,t3im,t3re;
  std::complex<prc_cmp> ray3,raytem,refwat;
  const prc_cmp wvl_min_mcr=0.01,wvl_max_mcr=1.0e+8;  
  const prc_cmp wvl_ctr_mcr(wvl_ctr*1.0e6); // [m]->[um] Wavelength at band center
  
  if(wvl_ctr_mcr < wvl_min_mcr || wvl_ctr_mcr > wvl_max_mcr) std::cout << "refwat--wavelength out of table range" << std::endl;
  
  if(wvl_ctr_mcr < 10){
    for(i=2;i<=idx_rfr_h2o_lqd_tbl_nbr;i++)
      if(wvl_ctr_mcr <= idx_rfr_h2o_lqd_tbl_wvl[i]) break;
    frac=std::log(wvl_ctr_mcr/idx_rfr_h2o_lqd_tbl_wvl[i-1])/std::log(idx_rfr_h2o_lqd_tbl_wvl[i]/idx_rfr_h2o_lqd_tbl_wvl[i-1]);
    rred=idx_rfr_h2o_lqd_tbl_rl[i-1]+frac*(idx_rfr_h2o_lqd_tbl_rl[i]-idx_rfr_h2o_lqd_tbl_rl[i-1]);
    rimd=idx_rfr_h2o_lqd_tbl_img[i-1]*std::pow((idx_rfr_h2o_lqd_tbl_img[i]/idx_rfr_h2o_lqd_tbl_img[i-1]),frac);
    refwat=std::complex<prc_cmp>(rred,rimd);
  }else if(wvl_ctr_mcr >= 10 && wvl_ctr_mcr <= 1e7){
    // Temperature correction for data using Ray model
    ray3=RayWater(wvl_ctr,298.0);
    t3re=ray3.real();
    t3im=ray3.imag();
    raytem=RayWater(wvl_ctr,tpt);
    corre=raytem.real()-t3re;
    corim=raytem.imag()-t3im;
    for(i=2;i<=idx_rfr_H2O::idx_rfr_h2o_lqd_tbl_nbr;i++)
      if(wvl_ctr_mcr <= idx_rfr_h2o_lqd_tbl_wvl[i]) break;
    frac=std::log(wvl_ctr_mcr/idx_rfr_h2o_lqd_tbl_wvl[i-1])/std::log(idx_rfr_h2o_lqd_tbl_wvl[i]/idx_rfr_h2o_lqd_tbl_wvl[i-1]);
    rred=idx_rfr_h2o_lqd_tbl_rl[i-1]+frac*(idx_rfr_h2o_lqd_tbl_rl[i]-idx_rfr_h2o_lqd_tbl_rl[i-1]);
    rimd=idx_rfr_h2o_lqd_tbl_img[i-1]*std::pow((idx_rfr_h2o_lqd_tbl_img[i]/idx_rfr_h2o_lqd_tbl_img[i-1]),frac);
    if(wvl_ctr_mcr >= 10 && wvl_ctr_mcr < 20){
      // weighted average
      rre=((wvl_ctr_mcr-10)/10.0)*corre+rred;
      rim=((wvl_ctr_mcr-10)/10.0)*corim+rimd;
    } else if(wvl_ctr_mcr >= 20 && wvl_ctr_mcr <= 1e7){
      rre=rred+corre;
      rim=rimd+corim;
    }
    refwat=std::complex<prc_cmp>(rre,rim);   
  } else if(wvl_ctr_mcr > 1e7){
    // use Ray model
    refwat=RayWater(wvl_ctr,tpt);
  }
  
  *idx_rfr=refwat;
  return rcd;
} // end idx_rfr_H2O_lqd_get()

int // O [enm] Return success code
idx_rfr_H2O_ice_get // [fnc] Refractive index of liquid water
(const prc_cmp wvl_ctr, // I [m] Wavelength at band center
 const prc_cmp tpt, // I [K] Temperature
 std::complex<prc_cmp> *idx_rfr) // O [frc] Refractive index of ice water
{
  /* Purpose: Refractive index of ice
     Calculates complex refractive index of ice between 45 nm and 8.6 m
     Temperature dependence is included for 213 < T < 272 K for wvl_ctr > 167 um
     
     Method :  tabular interpolation
     (1) real index linearly in log(wavelength)
     (2) real index linearly in temperature
     (3) log(imag. index) linearly in log(wavelength)
     (4) log(imag. index) linearly in temperature
     
     Author:  Prof. Stephen Warren, Univ. of Washington (Sept., 1983)
     
     Input parameters
     wvl_ctr = Wavelength [m]
     tpt     = Temperature [K]
     
     Output
     idx_rfr = Complex index of refraction (positive imaginary part)
     
     References: 
     War84: 
     Warren, S.G., Optical Constants of Ice from the Ultraviolet to the Microwave,
     Applied Optics, vol. 23, pp. 1206--1225, 1984.

     PeG91:
     Perovich, Donald K. and John W. Govoni,
     Absorption Coefficients Of Ice From 250 To 400 nm,
     Geophysical Research Letters, 18:7, pp. 1233--1235, July 1991.
     
     Copied from prime directory "patty" on 1987/12/15 J Tillman
     
     19920905: Included Perovich and Govoni, Charlie Zender zender@ncar.ucar.edu 
     
     This should bring the values up to date except for the imaginary 
     indices from 180-250 nm, which remain too low, according to Warren.

     20080125: Add two new references, GLC95 and RGC01
     GLC95:
     Gosse, S., D. Labrie, and P. Chylek (1995), 
     Refractive index of ice in the 1.4-7.8 um spectral range, 
     Appl. Opt., 34(28), 6582-6586.

     RGC01:
     Rajaram, B., D. L. Glandorf, D. B. Curtis, M. A. Tolbert, O. B. Toon, and N. Ockman (2001), 
     Temperature-dependent optical constants of water ice in the near infrared: new results and a critical review of the available measurements, 
     Appl. Opt., 40(25), 4449-4462.

     20080125: Updated to Warren's 1995 REFICE 
     Local source: /data/zender/libRadtran-1.3/src_f/REFICE.f
     Header notes from that file

        Calculates complex refractive index of Ice 1H for wavelengths
        between 45 nm and 8.6 m.  For wavelengths above 167 microns,
        temperature dependence is included for temperatures between
        213 and 272K.  Mainly intended for applications in Earth ice
        clouds and snow, not other planets or interstellar space;
        the temperature dependence or crystalline form of ice may be
        incorrect for these latter applications.


      I N P U T :  WAVLEN = wavelength (microns)
                            (range:  0.0443 to 8.600E+06)

                   TEMP   = temperature (K) ( for WAVLEN.GT.167 only )
                            (range:  213.16 to 272.16)

      O U T P U T :  REFICE = complex refractive index
                              ( with positive imaginary part )

      (WARNING:  input out of range will print a warning message and
                 return REFICE=(0.,0.) in order not to unnecessarily 
                 halt the calling program;  the calling program should
                 test the real part of REFICE to catch these errors)


      METHOD :  Tabular interpolation, assuming

                (1) real index is linear in log(wavelength)
                    and linear in temperature

                (2) log(imag. index) is linear in log(wavelength)
                    and linear in temperature


     AUTHORS : Stephen Warren, Univ. of Washington (1983)
               (sgw@cloudy.atmos.washington.edu)

               Bo-Cai Gao, JCESS, Univ. of Maryland (1995)
               (gao@imagecube.gsfc.nasa.gov)

               Warren Wiscombe, NASA Goddard (1995)
               (wiscombe@climate.gsfc.nasa.gov)

     MODIFICATIONS IN 1995 :

       Gao, Warren, and (to a small extent) Wiscombe modified the
       original Warren REFICE program from 1984 to change values of
       imaginary refractive index in the 0.161-0.410 and 1.445-2.50
       micron regions.  The values in 0.161-0.410 were incorrect and
       the values in 1.445-2.50 were among the most uncertain in 1984.
       New measurements have made it possible to improve both regions.

       No changes were made to real refractive indices (by re-doing a
       Kramers-Kronig analysis), because the values of imaginary
       index MIM involved are so small (below 0.001) that the
       resulting changes in real index MRE would be in the third
       decimal place at most.  (MIM has negligible effect on MRE
       when MIM << MRE.)

       The 0.161-0.410 micron region was changed using data provided
       by Warren, which correct his misinterpretation of Minton's
       measurements for 0.181-0.185 micron, and incorporate new
       measurements of Perovich and Govoni (1991) for 0.250-0.400
       micron.  Warren (1984) correctly represented UV measurements
       of Seki et al. and visible measurements of Grenfell/Perovich,
       but he plotted Minton's measurements a factor of 2.3 too low
       because he misinterpreted base-10 as base-e.  (The UV
       measurements of Dressler/Schnepp and Shibaguchi et al are also
       probably expressed as absorption coefficients on base-10;
       therefore those values also were probably plotted a factor of
       2.3 too low in Warren's (1984) Figure 2.)

       The details of how the present imaginary index data for
       0.161-0.410 micron is obtained are as follows.  Point A in
       Warren's Figure 2 at 161 nm is joined with a straight line to
       Minton's corrected point B at 181 nm.  Minton's reported
       values for 181-185 nm have been smoothed within his stated
       uncertainty.  Now a smooth curve is drawn to join Minton at
       185 nm to Perovich/Govoni (PG) at 250 nm.  PG's values from
       their Table 1 show some unrealistic wiggles that are smaller
       than their error bars, so a smooth curve was fitted through
       them and values were taken from the smoothed curve at 10-nm
       intervals.  PG ends at 400 nm, where Grenfell/Perovich (GP)
       starts.  At 400 nm we take imaginary index=2.82E-9, the
       average of PG (2.93E-9) and GP (2.71E-9).

       The Warren (1984) values of imaginary index in the 1.445-2.50
       micron region were replaced by those of Kou et al.(1993).  In
       order to remove the resulting discontinuities near 1.445 and
       2.5 micron, the Warren values at 1.43 and 1.44 micron were
       changed to 0.9E-04 and 1.3E-04 respectively, and his values at
       2.52, 2.55, and 2.565 micron were changed to 8.255E-04,
       8.578E-04, and 8.739E-04, respectively. The latter change
       eliminated a small local maximum at 2.5 micron which was not
       realistic and has never been seen in spectra of snow bracketing
       that wavelength.

     REFERENCES :

       Warren, S., 1984: Optical Constants of Ice from the Ultraviolet
          to the Microwave, Appl. Opt. 23, 1206-1225

       Kou, L., D. Labrie, and P. Chylek, 1993: Refractive indices
          of water and ice in the 0.65- to 2.5-micron spectral range,
          Appl. Opt. 32, 3531-3540

       Perovich, D., and J. Govoni, 1991: Absorption Coefficients
          of Ice from 250 to 400 nm, Geophys. Res. Lett. 18, 1233-1235 */
  int rcd(0); // [enm] Return success code
  int idx,idx_temp;
  prc_cmp mRe,mIm,yLo,yHi,frac;
  const prc_cmp wvl_ctr_mcr(wvl_ctr*1.0e6); // [m]->[um] Wavelength at band center
  
  if(wvl_ctr_mcr < 0.045 || wvl_ctr_mcr > 8.6e6){
    std::cerr << "ERROR: No refractive index data for wavelength " << wvl_ctr_mcr << std::endl;
    std::abort();
  } // endif
  
  if(wvl_ctr_mcr <= 167){
    // Temperature-independent for wavelengths 0.045--167.0 um
    for(idx=2;idx <= idx_rfr_H2O::nwvl;idx++)
      if(wvl_ctr_mcr <= wvl_ice_tbl_mcr[idx]) break;
    frac=std::log(wvl_ctr_mcr/wvl_ice_tbl_mcr[idx-1])/std::log(wvl_ice_tbl_mcr[idx]/wvl_ice_tbl_mcr[idx-1]);
    mRe=tabRe[idx-1]+frac*(tabRe[idx]-tabRe[idx-1]);
    mIm=tabIm[idx-1]*std::pow((tabIm[idx]/tabIm[idx-1]),frac);
  }else{
    // Temperature-dependent for wavelengths longer than 167.0 um
    if(tpt < tempRef[ntemp-1] || tpt > tempRef[0]){
      std::cerr << "ERROR: No refractive index data for ice at temperature"
		<< tpt << std::endl;
      std::abort();
    }
    // find position in temperature array
    for(idx_temp=1;idx_temp < ntemp;idx_temp++)
      if(tpt >= tempRef[idx_temp]) break;
    // find position in wavelength array
    for(idx=2;idx <= idx_rfr_H2O::nwvlT;idx++)
      if(wvl_ctr_mcr <= wvlT[idx]) break;
    frac= std::log(wvl_ctr_mcr/wvlT[idx-1])/ std::log(wvlT[idx]/wvlT[idx-1]);
    yLo=tabReT[idx_temp][idx-1]+frac*(tabReT[idx_temp][idx]-tabReT[idx_temp][idx-1]);
    yHi=tabReT[idx_temp-1][idx-1]+frac*(tabReT[idx_temp-1][idx]-tabReT[idx_temp-1][idx-1]);
    mRe=yLo+(yHi-yLo)*(tpt-tempRef[idx_temp])/(tempRef[idx_temp-1]-tempRef[idx_temp]);
    yLo=std::log(tabImT[idx_temp][idx-1])+ frac*std::log(tabImT[idx_temp][idx]/tabImT[idx_temp][idx-1]);
    yHi=std::log(tabImT[idx_temp-1][idx-1])+frac*std::log(tabImT[idx_temp-1][idx]/tabImT[idx_temp-1][idx-1]);
    mIm=std::exp(yLo+(yHi-yLo)*(tpt-tempRef[idx_temp])/(tempRef[idx_temp-1]-tempRef[idx_temp]));
  } // end if wavelength longer than 167.0 um
  
  *idx_rfr=std::complex<prc_cmp>(mRe,mIm);
  return rcd;
} // end idx_rfr_H2O_ice_get()

prc_cmp ab(prc_cmp wvl,prc_cmp bet,prc_cmp wvlcen,prc_cmp del,prc_cmp gam)
{
  prc_cmp result;
  result=bet*std::exp(-1*std::pow(std::fabs(std::log10(wvl/wvlcen)/del),gam));
  return result;
} // end ab()

prc_cmp tbarp1(prc_cmp tpt_cls,prc_cmp wvl)
{
  prc_cmp result;
  result= 1.0+1.0e-4*(tpt_cls-25.0)*std::exp(std::pow((wvl/4.0),0.25));
  return result;
} // end tbarp1()

prc_cmp term(prc_cmp wvl,prc_cmp frec,prc_cmp bet,prc_cmp gam)
{
  prc_cmp result;
  prc_cmp wvn_sqr((PRC_CMP(1.0e4)/wvl)*(PRC_CMP(1.0e4)/wvl));
  prc_cmp fct(frec*frec-wvn_sqr);
  result=bet*fct/(fct*fct+gam*wvn_sqr);
  return result;
} // end term()

prc_cmp brre(prc_cmp tpt_cls,prc_cmp wvl)
{
  prc_cmp result;
  result=
    tbarp1(tpt_cls,wvl)*
    std::sqrt(1.83899+term(wvl,1639.0,52340.4,10399.2)
	      +term(wvl,588.24,345005.0,259913.0)
	      +term(wvl,161.29,43319.7,27661.2));
  return result;
} // end brre()

prc_cmp arre(prc_cmp tpt_cls,prc_cmp wvl)
{
  prc_cmp result;
  result=
    tbarp1(tpt_cls,wvl)*
    std::sqrt(1.79907+term(wvl,3352.27,999140.0,151963.0) 
	      +term(wvl,1639.00,50483.5,9246.27) 
	      +term(wvl,588.24,844697.0,1076150.0));
  return result;
} // end arre()

void idx_rfr_H2O_tst()
{
  /* Purpose: Test correct conversion of refractive index routines
     Usage: Simply call idx_rfr_H2O_tst() to test RayWater, idx_rfr_H2O_lqd_get, idx_rfr_H2O_ice_get */
  int rcd(0); // [enm] Return success code
  
  const int wvl_nbr(10);
  const prc_cmp wvl_ctr[wvl_nbr]={0.2e-6,0.5e-6,1.0e-6,2.0e-6,4.0e-6,8.0e-6,16.0e-6,32.0e-6,64.0e-6,128.0e-6};
  const prc_cmp tpt[wvl_nbr]={200,210,220,230,240,250,260,270,280,290};
  std::complex<prc_cmp> idx_rfr;
  
  std::cout << "Program to test idx_rfr_H2O.cc" << std::endl;
  
  std::cout << ">>> Check RayWater()" << std::endl;
  for(int wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
    std::cout << "Wavelength  = " << wvl_ctr[wvl_idx] << std::endl;
    std::cout << "Temperature = " << tpt[wvl_idx] << std::endl;
    std::cout << "RayWater    = " << RayWater(wvl_ctr[wvl_idx],tpt[wvl_idx]) << std::endl;
    std::cout << std::endl;
  } // end loop over wvl
  
  std::cout << ">>> Check idx_rfr_H2O_lqd_get()" << std::endl;
  for(int wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
    rcd=idx_rfr_H2O_lqd_get // [fnc] Refractive index of liquid water
      (wvl_ctr[wvl_idx], // I [m] Wavelength at band center
       tpt[wvl_idx], // I [K] Temperature
       &idx_rfr); // O [frc] Refractive index of liquid water
    std::cout << "Wavelength  = " << wvl_ctr[wvl_idx] << std::endl;
    std::cout << "Temperature = " << tpt[wvl_idx] << std::endl;
    std::cout << "idx_rfr_H2O_lqd_get (re) = " << idx_rfr.real() << std::endl;
    std::cout << "idx_rfr_H2O_lqd_get (im) = " << idx_rfr.imag() << std::endl;
    std::cout << std::endl;
  } // end loop over wvl
  
  std::cout << ">>> Check idx_rfr_H2O_ice_get()" << std::endl;
  for(int wvl_idx=0;wvl_idx<wvl_nbr;wvl_idx++){
    rcd=idx_rfr_H2O_ice_get // [fnc] Refractive index of ice water
      (wvl_ctr[wvl_idx], // I [m] Wavelength at band center
       tpt[wvl_idx], // I [K] Temperature
       &idx_rfr); // O [frc] Refractive index of ice water
    std::cout << "Wavelength  = " << wvl_ctr[wvl_idx] << std::endl;
    std::cout << "Temperature = " << tpt[wvl_idx] << std::endl;
    std::cout << "idx_rfr_H2O_ice_get (re) = " << idx_rfr.real() << std::endl;
    std::cout << "idx_rfr_H2O_ice_get (im) = " << idx_rfr.imag() << std::endl;
    std::cout << std::endl;
  } // end loop over wvl
} // end idx_rfr_H2O_ice_get_tst()

