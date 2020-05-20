// $Id$ 

// Purpose: Implementation (declaration) of Mie scattering solutions and utilities

/* Copyright (C) 1997--present Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

/* Library of functions for Bohren & Huffman Mie scattering solvers
   Source: 
   mie_sph_BoH83(), mie_ngl_BoH83(), mie_bck_hms_Chy73() 
   Translated from Fortran77 into C++ by Jeff G. Wong 19960319
   Modifications by Charlie Zender:
   1996 Add additional diagnostics
   1998 Adapted to ANSI C++
   1999 Cleaned up: Return -1 on error, 0 on success
   2000 Add coated sphere solution
   2001 Improved casting for ISO compliance, std::pow(int,int)
   2004 Add Wis79 solution
   2005 Add MaS99 solution */
   
#include <mie_sln.hh> // Mie scattering solutions and utilities

int // O [enm] Return success code
mie_prc // [fnc] Mie processor
(const bool abs_ncl_wk_mdm_flg, // I [flg] Absorbing inclusion in weakly-absorbing sphere (MaS99)
 const bool coat_flg, // I [flg] Assume coated spheres
 const bool slf_tst_flg, // I [frg] Self-test flag
 const long bnd_idx, // I [idx] Counting index for band
 const long lgn_nbr, // I [nbr] Order of phase function Legendre expansion
 const long ngl_nbr, // I [nbr] Angle number (number of angles is 2*ngl_nbr-1)
 const long sz_idx, // I [idx] Counting index for size
 const long wvl_idx, // I [idx] Counting index for wavelength
 const long wvl_idx_dbg, // I [idx] Debugging wavelength bin
 const prc_cmp bnd_ctr, // I [m] Wavelength at band center
 const prc_cmp dmn_frc, // I [frc] Fractal dimensionality of inclusions
 const prc_cmp rds_cor, // I [m] Radius of core
 const prc_cmp rds_mnt, // I [m] Radius of mantle
 const prc_cmp sz_ctr_sph, // I [m] Equivalent sphere size at bin center
 const prc_cmp sz_prm_rsn_usr_spc, // I [m m-1] Size parameter resolution, user specified
 const std::complex<prc_cmp> idx_rfr_cor, // I [frc] Refractive index of core
 const std::complex<prc_cmp> idx_rfr_ffc, // I [frc] Effective refractive index of particle
 const std::complex<prc_cmp> idx_rfr_mdm, // I [frc] Refractive index of medium
 const std::complex<prc_cmp> idx_rfr_mnt, // I [frc] Refractive index of mantle
 const std::complex<prc_cmp> idx_rfr_mtx, // I [frc] Refractive index of matrix
 const std::complex<prc_cmp> idx_rfr_ncl, // I [frc] Refractive index of inclusion
 const std::complex<prc_cmp> idx_rfr_prt, // I [frc] Refractive index of particle
 const std::string slv_sng, // I [sng] Mie solver to use (BoH83 or Wis79)
 double &abs_fct_MaS99, // O [frc] Absorption enhancement for inclusions in weakly-absorbing spheres
 double &asm_prm, // O [frc] Asymmetry parameter
 double &bck_hms, // O [frc] Hemispheric backscatter
 double &q_abs, // O [frc] Absorption efficiency
 double &q_bck, // O [frc] Backscattering efficiency
 double &q_ext, // O [frc] Extinction efficiency
 double &q_sct, // O [frc] Scattering efficiency
 double &spk_val, // O [frc] Mie coefficient of spike
 const prc_cmp * const ngl, // I [rdn] Scattering angle
 const prc_cmp * const ngl_dlt, // I [rdn] Width of angle bin
 prc_cmp * const phz_fnc, // O [sr-1] Phase function
 prc_cmp * const plz) // O [frc] Degree of linear polarization
{
  /* Purpose: Driver for optical efficiencies solvers
     Call Mie scattering solution for spheres or coated spheres, or
     call anomalous diffraction theory as appropriate
  
     Conventions: 
     Seven different refractive indices are passed as input
     At most three are actually used, in most cases only two are used
     idx_rfr_mdm is always used for refractive index of medium
     idx_rfr_ffc is always used for refractive index of particle unless
     1. coat_flg == true -> idx_rfr_cor and idx_rfr_mnt are used
     2. abs_ncl_wk_mdm_flg == true -> idx_rfr_mtx is used in mie_MaS99
     3. 
     Hence routine assumes all pre-processing of effective medium approximations (EMAs)
     has already occurred and correct EMA, if any, is stored in idx_rfr_ffc
     These assumption can and will be improved someday, probably by adding an
     enumerated type with more fine-grained control of components
     This would make it simple, for example, to specify that idx_rfr_ffc
     (rather than the default idx_rfr_cor) be used for for cores of coated particles.
     It is intended that idx_rfr_mtr and idx_rfr_ncl be used for homogeneously distributed
     absorbing inclusions in water spheres once Markel's codes are implemented. */

  int rcd(0); // [enm] Return success code
  unsigned short dbg_lvl(dbg_lvl_get()); // [enm] Debugging level
  const std::string sbr_nm("mie_prc"); // [sng] Name of subroutine
  // Work variables for Mie routines
  long ngl_idx; // [idx] Counting index for angle
  prc_cmp sz_prm_wrn(30.0); // [frc] Largest size parameter before warning
  std::complex<double> *s1=new std::complex<double>[2*ngl_nbr]; // [frc] Scalar amplitude scattering matrix element S1
  std::complex<double> *s2=new std::complex<double>[2*ngl_nbr]; // [frc] Scalar amplitude scattering matrix element S2
  // fxm: Why make these valarrays? Are valarrays guaranteed to be consecutive?
  //  std::valarray<std::complex<double> > s1(2*ngl_nbr); // [frc] Scalar amplitude scattering matrix element S1
  //  std::valarray<std::complex<double> > s2(2*ngl_nbr); // [frc] Scalar amplitude scattering matrix element S2
  double *XMU=new double[2*ngl_nbr]; // [frc] Cosine polar scattering angle  Wis79 p. 16 (9B)

  /* Size parameter is 2*pi*r*idx_rfr_mdm/lambda_vacuum, where
     lambda_vacuum is wavelength in vacuum. See van de Hulst (1981, p. 130) */
  // fxm: Verify factor is real(idx_rfr_mdm) and not, e.g., norm(idx_rfr_mdm)
  const double sz_prm(2.0*mth::cst_M_PIl*sz_ctr_sph*real(idx_rfr_mdm)/bnd_ctr); // [m m-1] Size parameter
  
  if(dbg_lvl >= dbg_io) std::cerr << "Calling mie_sph_BoH83() for wvl_idx = " << wvl_idx << ", bnd_idx = " << bnd_idx << ", sz_idx = " << sz_idx << ", sz_prm = " << sz_prm << std::endl;

  // Standard spherical solution unless coated sphere is indicated
  if(!coat_flg){
    const std::complex<double> idx_rfr_rlt=static_cast<std::complex<double> >(idx_rfr_ffc/idx_rfr_mdm); // [frc] Effective refractive index of particle (e.g., aerosol) relative to medium (e.g., air)

    if(slv_sng == "BoH83"){
      // Method of Bohren & Huffman (1983)
      rcd+=mie_sph_BoH83 // [fnc] Mie solution for homogeneous spheres, BoH83
      (sz_prm, // I [m m-1] Size parameter
       idx_rfr_rlt, // I [frc] Refractive index of particle (aerosol) relative to medium (air)
       ngl_nbr, // I [nbr] Angle number (number of angles is 2*ngl_nbr-1)
       ngl, // I [rdn] Scattering angle
       asm_prm, // O [frc] Asymmetry parameter
       bck_hms, // O [frc] Hemispheric backscatter
       q_bck, // O [frc] Backscattering efficiency
       q_ext, // O [frc] Extinction efficiency
       q_sct, // O [frc] Scattering efficiency
       &s1[0], // O [frc] Scalar amplitude scattering matrix element S1
       &s2[0]); // O [frc] Scalar amplitude scattering matrix element S2

    }else if(slv_sng == "Wis79"){
      /* Wis79 Algorithm interface
	 Input variables (actually, all _dimension variables are input):
	 double XX, std::complex<double> CREFIN, bool PERFCT, double MIMCUT,
	 bool ANYANG, int NUMANG, double XMU[], int XMU_dimension,
	 int NMOM, int IPOLZN, int MOMDIM, bool PRNT[], int PRNT_dimension,
	 Output variables:
	 double &QEXT, double &QSCA, double &GQSC,
	 double PMOM[][5], int PMOM_dimension, std::complex<double> &SFORW,
	 std::complex<double> &SBACK, std::complex<double> S1[], int S1_dimension, std::complex<double> S2[], int S2_dimension,
	 std::complex<double> TFORW[], int TFORW_dimension, std::complex<double> TBACK[], int TBACK_dimension,
	 double &SPIKE ) */
      
      // Primary dimensions
      // "+1"'s indicate dimension is expanded by one relative to Fortran to equivalence subscripting
      const bool NOPMOM(true); // [flg] Do not compute Legendre moments
      const int MAXANG(2*static_cast<int>(ngl_nbr)-1); // [nbr] Maximum number of angles in angular integration
      /* fxm: Dynamically dimensioning MOMDIM will be difficult since there is a
	 dependency between it and MAXANG. Leave MOMDIM at 200 for now. */
      const int MOMDIM(200); // [nbr] Maximum number of Legendre moments
      //      const int MOMDIM(lgn_nbr); // [nbr] Maximum number of Legendre moments
      const int PRNT_dimension(2+1); // [nbr] Rank of print code
      const int TFORW_dimension(2+1); // [nbr] Rank of scattering amplitude ratio arrays
      const int TBACK_dimension(2+1); // [nbr] Rank of scattering amplitude ratio arrays

      // Derived dimensions
      const int S1_dimension(MAXANG+1); // [frc] Rank of scattering amplitude matrix S1
      const int S2_dimension(MAXANG+1); // [frc] Rank of scattering amplitude matrix S2
      const int PMOM_dimension(MOMDIM+1); // [nbr] Rank of Legendre moment arrays
      const int XMU_dimension(MAXANG+1); // [nbr] Rank of scattering angle array

      // Input variables:
      double XX; // [m m-1] Size parameter Wis79 p. 15 (1)
      std::complex<double> CREFIN; // [frc] Complex refractive index of particle relative to medium Wis79 p. 15 (2)
      bool PERFCT; // [flg] Perfectly conducting/reflecting sphere Wis79 p. 58 (A.1)
      double MIMCUT; // [frc] Threshold imaginary index of refraction for non-absorption 
      bool ANYANG; // [flg] Pick arbitrary angles in angular integration
      int NUMANG; // [nbr] Number of angles in angular integration
      int NMOM; // [nbr] Number of Legendre moments [~2*sz_prm]
      int IPOLZN; // [enm] Polarization code
      bool PRNT[PRNT_dimension]; // [enm] Code for printing diagnostics
      // Output variables:
      double QEXT; // [frc] Extinction efficiency Wis79 p. 15 (6), p. 58 (A.1), Wis80 p. 1505 (1a)
      double QSCA; // [frc] Scattering efficiency Wis79 p. 15 (7), p. 58 (A.1), Wis80 p. 1505 (1b)
      double GQSC; // [frc] Asymmetry factor times scattering efficiency Wis79 p. 58 (A.1), Wis80 p. 1505 (1c)
      double PMOM[PMOM_dimension][PMOM_second_dimension]; // [frc] Legendre polynomial expansion coefficients (moments) 
      std::complex<double> SFORW; // [frc] Forward scattering amplitude Wis79 p. 16 (9), p. 60 (B.2)
      std::complex<double> SBACK; // [frc] Backscattering amplitude Wis79 p. 16 (9), p. 60 (B.1)
      std::complex<double> *S1; // [frc] Scalar amplitude scattering matrix element S1 Wis79 p. 16 (9a, 11a), Wis80 p. 1505 (1d)
      std::complex<double> *S2; // [frc] Scalar amplitude scattering matrix element S2 Wis79 p. 16 (9b, 11b), Wis80 p. 1505 (1e)
      //      std::complex<double> S1[S1_dimension]; // [frc] Scalar amplitude scattering matrix element S1 Wis79 p. 16 (9a, 11a), Wis80 p. 1505 (1d)
      //      std::complex<double> S2[S2_dimension]; // [frc] Scalar amplitude scattering matrix element S2 Wis79 p. 16 (9b, 11b), Wis80 p. 1505 (1e)
      // TFORW=[C filler, T1(1), T2(1)] and TBACK=[C filler, T1(-1), T2(-1)] on Wis79 p. 60 (B.5--B.8)
      std::complex<double> TFORW[TFORW_dimension]; // [frc] Forward ratio of scattering amplitudes useful for polarized RT Wis79 p. 60 (B.5, B.6)
      std::complex<double> TBACK[TBACK_dimension]; // [frc] Backward ratio of scattering amplitudes useful for polarized RT Wis79 p. 60 (B.7, B.8)
      double SPIKE; // [frc] Spike diagnostic Wis79 p. 62 (B.10)

      // Initialize Wis79 input variables with appropriate values from driver
      XX=sz_prm; // [m m-1] Size parameter Wis79 p. 15 (1)
      /* Wis79 and BoH83 use opposite conventions for imaginary part of index of refraction so amplitude scattering matrices returned are complex conjugates of eachother 
	 Wis79 uses m=n_r - n_i*i while BoH83 uses m=n_r + n_i*i */
      CREFIN=conj(idx_rfr_rlt); // [frc] Complex refractive index of particle relative to medium Wis79 p. 15 (2)
      PERFCT=false; // [flg] Perfectly conducting/reflecting sphere Wis79 p. 58 (A.1)
      MIMCUT=1.0e-12; // [frc] Threshold imaginary index of refraction for non-absorption 
      ANYANG=true; // [flg] Pick arbitrary angles in angular integration
      NUMANG=MAXANG; // [nbr] Number of angles in angular integration

      // Once method works with default (7) angles, use angles specified by host mie program
      for(ngl_idx=1;ngl_idx<=2*ngl_nbr-1;ngl_idx++){ // NB: Loop index starts at 1
	XMU[ngl_idx]=std::cos(ngl[ngl_idx-1]); // [frc] Cosine polar scattering angle Wis79 p. 16 (9B);
      } // end loop over angle */

      NMOM=static_cast<int>(2.0*sz_prm); // [nbr] Number of Legendre moments [~2*sz_prm]
      if((XX < 1.0) || (XX > 100.0)) NMOM=1; // [nbr] Number of Legendre moments [~2*sz_prm]
      if(XX > 1000.0 ) NMOM=0; // [nbr] Number of Legendre moments [~2*sz_prm]
      if(NOPMOM) NMOM=0; // [nbr] Number of Legendre moments [~2*sz_prm]

      IPOLZN=1234; // [enm] Polarization code (fxm: Meaning?)
      PRNT[1]=true; // [enm] Code for printing diagnostics (fxm: Meaning?)
      PRNT[2]=true; // [enm] Code for printing diagnostics (fxm: Meaning?)
      
      // Initialize output buffers shared between various Mie Solvers
      S1=s1; // [frc] Scalar amplitude scattering matrix element S1 Wis79 p. 16 (9a, 11a), Wis80 p. 1505 (1d)
      S2=s2; // [frc] Scalar amplitude scattering matrix element S2 Wis79 p. 16 (9b, 11b), Wis80 p. 1505 (1e)

      // Method of Wiscombe (1979)
      rcd+=mie_sph_Wis79 // [fnc] Mie solution for homogeneous spheres (MIEV0), Wis79
	( XX, CREFIN, PERFCT, MIMCUT,
	  ANYANG, NUMANG, XMU, XMU_dimension,
	  NMOM, IPOLZN, MOMDIM, PRNT, PRNT_dimension,
	  QEXT, QSCA, GQSC,
	  PMOM, PMOM_dimension, SFORW,
	  SBACK, S1, S1_dimension, S2, S2_dimension,
	  TFORW, TFORW_dimension, TBACK, TBACK_dimension,
	  SPIKE );

      if(dbg_lvl == dbg_old){
	int tbl_clm_wdt(8); // [nbr] Width of table columns

	std::cerr << prg_nm_get() << ": INFO Wis79 diagnostics from " << sbr_nm << "():" << std::endl;
	std::cerr << "MAXANG = " << MAXANG << " = Maximum number of angles in angular integration" << std::endl;
	std::cerr << "MOMDIM = " << MOMDIM << " = Maximum number of Legendre moments" << std::endl;
	std::cerr << "NMOM = " << NMOM << " = Number of Legendre moments [~2*sz_prm]" << std::endl;
	std::cerr << "S1_dimension = MAXANG+1 = " << MAXANG+1 << std::endl;
	std::cerr << "PMOM_dimension = MOMDIM+1 = " << MOMDIM+1 << std::endl;
	for(ngl_idx=1;ngl_idx<=NUMANG;ngl_idx++){
	  std::cerr << "S1[" << std::setw(2) << ngl_idx << "] = " << std::setw(tbl_clm_wdt) << S1[ngl_idx] << ", S2[" << ngl_idx << "] = " << std::setw(tbl_clm_wdt) << S2[ngl_idx] << std::endl;
	} // end loop over ngl

	// fxm: Understand and implement PMOM diagnostics
	for(int lgn_idx=1;lgn_idx<=PMOM_dimension;lgn_idx++){
	  for(int lgn2_idx=1;lgn2_idx<=PMOM_second_dimension;lgn2_idx++){
	    std::cerr << "PMOM[" << std::setw(2) << lgn_idx << "][ " << lgn2_idx << "] = " << std::setw(tbl_clm_wdt) << PMOM[lgn_idx][lgn2_idx] << std::endl;
	    ;
	  } // end loop over lgn2
	} // end loop over lgn
      } // end if dbg

      /* Read Wis79 output variables into appropriate values from driver
	 Variables from driver that are not set by Wis79 should be zeroed
	 Otherwise program will fail in netCDF write of unrepresentable data */
      q_ext=QEXT; // [frc] Extinction efficiency Wis79 p. 15 (6), p. 58 (A.1), Wis80 p. 1505 (1a)
      q_sct=QSCA; // [frc] Scattering efficiency Wis79 p. 15 (7), p. 58 (A.1), Wis80 p. 1505 (1b)
      // Verified 20040531 that S1[NUMANG]=s1[2*ngl_nbr-1]
      q_bck=4.0/(sz_prm*sz_prm)*norm(S1[NUMANG]); // O [frc] Backscattering efficiency evaluated by optical theorem BoH83 p. 121 (4.82)

      // Derived output
      if(q_sct == 0.0) asm_prm=0.0; else asm_prm=GQSC/q_sct; // [frc] Asymmetry factor times scattering efficiency Wis79 p. 58 (A.1), Wis80 p. 1505 (1c)

      for(ngl_idx=1;ngl_idx<=2*ngl_nbr-1;ngl_idx++){ // NB: Loop index starts at 1
	s1[ngl_idx]=S1[ngl_idx]; // [frc] Scalar amplitude scattering matrix element S1 Wis79 p. 16 (9a, 11a), Wis80 p. 1505 (1d)
	s2[ngl_idx]=S2[ngl_idx]; // [frc] Scalar amplitude scattering matrix element S2 Wis79 p. 16 (9b, 11b), Wis80 p. 1505 (1e)
      } // end loop over ngl

      PMOM[0][0]=PMOM[0][0]; // [frc] Legendre polynomial expansion coefficients (moments) 

      SFORW=SFORW; // [frc] Forward scattering amplitude Wis79 p. 16 (9), p. 60 (B.2)
      SBACK=SBACK; // [frc] Backscattering amplitude Wis79 p. 16 (9), p. 60 (B.1)

      // TFORW=[C filler, T1(1), T2(1)] and TBACK=[C filler, T1(-1), T2(-1)] on Wis79 p. 60 (B.5--B.8)
      TFORW[0]=TFORW[0]; // [frc] Forward ratio of scattering amplitudes useful for polarized RT Wis79 p. 60 (B.5, B.6)
      TBACK[0]=TBACK[0]; // [frc] Backward ratio of scattering amplitudes useful for polarized RT Wis79 p. 60 (B.7, B.8)
      SPIKE+=0.; // [frc] CEWI Spike diagnostic Wis79 p. 62 (B.10)
    } // end if slv_sng

    if(dbg_lvl >= dbg_io){
      int tbl_clm_wdt(8); // [nbr] Width of table columns
      
      std::cerr << prg_nm_get() << ": INFO Wis79 diagnostics from " << sbr_nm << "():" << std::endl;
      for(ngl_idx=1;ngl_idx<=2*ngl_nbr-1;ngl_idx++){ // NB: Loop index starts at 1
	std::cerr << "s1[" << std::setw(2) << ngl_idx << "] = " << std::setw(tbl_clm_wdt) << s1[ngl_idx] << ", s2[" << ngl_idx << "] = " << std::setw(tbl_clm_wdt) << s2[ngl_idx] << std::endl;
      } // end loop over ngl
    } // end if dbg
    
    if(slf_tst_flg){
      std::cerr << "\nCompare these results to Bohren & Huffman (1983), \n\"Absorption and Scattering of Light by Small Particles\" from Wiley Interscience Press\n" << std::endl;
      std::cerr << "p. 482 quotes values QSCA = 3.10543  QEXT =  3.10543 QBACK = 2.92534\n" << std::endl;
      
      std::cerr << "SPHERE SCATTERING PROGRAM\n\n" << std::endl;
      std::cerr << "\tREFMED = " << idx_rfr_mdm << "\tREFRE = " << idx_rfr_ffc.real() << "\tREFIM = " << idx_rfr_ffc.imag() << std::endl;
      std::cerr << "\tSPHERE RADIUS = " << sz_ctr_sph*1.0e6 << " \tWAVELENGTH = " << bnd_ctr*1.0e6 << std::endl;
      std::cerr << "\tSIZE PARAMETER = " << sz_prm << "\n" << std::endl;
      std::cerr << "\tQSCA = " << q_sct << "\tQEXT =  " << q_ext << "\tQBACK = " << q_bck << "\n" << std::endl;
    } // endif slf_tst_flg
    
    if(abs_ncl_wk_mdm_flg){
      /* Assume moderately-to-strongly absorbing inclusions are homogeneously mixed
	 in non-to-weakly absorbing medium.
	 Absorption enhancement relative to externally mixed inclusions as per MaS99
	 Usage of MaS99 for soot in cloud droplets:
	 Assign mdm=air, mtx=h2o_lqd, ncl=soot, prt=h2o_lqd
	 Mie program sets idx_rfr_ffc with prt unless ffc_mdm_typ != nil 
	 Best to set prt=h2o_lqd so MOPs and sz_prm sent to abs_fct_MaS99() represent cloud droplets
	 Absorption enhancement factors from these simulations should then be applied to soot absorption cross-sections from separate simulations */
      abs_fct_MaS99= // O [frc] Absorption enhancement for inclusions in weakly-absorbing spheres
	mie_sph_abs_fct_MaS99 // [fnc] Absorption enhancement for inclusions in weakly-absorbing spheres
	(sz_prm, // I [m m-1] Size parameter
	 idx_rfr_mtx.real(), // I [frc] Refractive index of weakly-absorbing sphere, real component
	 dmn_frc); // I [frc] Fractal dimensionality of inclusions
    } // endif abs_ncl_wk_mdm_flg

  } // endif !coat_flg
  
  // Coated sphere when indicated
  if(coat_flg){
    const double sz_prm_cor=2.0*mth::cst_M_PIl*rds_cor*real(idx_rfr_mdm)/bnd_ctr; // [m m-1] Size parameter of core
    const double sz_prm_mnt=2.0*mth::cst_M_PIl*rds_mnt*real(idx_rfr_mdm)/bnd_ctr; // [m m-1] Size parameter of mantle
    // Sanity check
    if(sz_prm_cor*idx_rfr_cor.imag() > 30.0 ||
       sz_prm_cor*idx_rfr_mnt.imag() > 30.0 ||
       sz_prm_mnt*idx_rfr_mnt.imag() > 30.0){
      err_prn(sbr_nm,"mie_sph_coat_BoH83() fails for large, highly absorbing spheres, see BoH83 p. 485");
      err_prn(sbr_nm,"Failing as safeguard");
    } // endif outside range of validity of mie_sph_coat_BoH83()

    const std::complex<double> idx_rfr_rlt_cor=std::complex<double>(idx_rfr_cor/idx_rfr_mdm); // [frc] Core refractive index relative to medium
    const std::complex<double> idx_rfr_rlt_mnt=std::complex<double>(idx_rfr_mnt/idx_rfr_mdm); // [frc] Mantle refractive index relative to medium
    rcd+=mie_sph_coat_BoH83 // [fnc] Mie solution for coated spheres
      (sz_prm_cor, // [m m-1] Size parameter of core
       sz_prm_mnt, // [m m-1] Size parameter of mantle
       idx_rfr_rlt_cor, // [frc] Core refractive index relative to medium
       idx_rfr_rlt_mnt, // [frc] Mantle refractive index relative to medium
       q_ext, // O [frc] Extinction efficiency
       q_sct, // O [frc] Scattering efficiency
       q_bck, // O [frc] Backscattering efficiency
       asm_prm, // O [frc] Asymmetry parameter
       bck_hms); // O [frc] Hemispheric backscatter

    if(slf_tst_flg){
      std::cerr << "p. 489 quotes values QSCA = 1.14341  QEXT =  2.32803 QBACK = 0.0285099\n" << std::endl;
      std::cerr << "COATED SPHERE SCATTERING PROGRAM\n\n" << std::endl;
      std::cerr << "\tREFMED = " << idx_rfr_mdm << std::endl;
      std::cerr << "\tREFRE1 = " << idx_rfr_cor.real() << "\tREFIM1 = " << idx_rfr_cor.imag() << std::endl;
      std::cerr << "\tREFRE2 = " << idx_rfr_mnt.real() << "\tREFIM2 = " << idx_rfr_mnt.imag() << std::endl;
      std::cerr << "\tCORE RADIUS = " << rds_cor*1.0e6 << "\tCOAT RADIUS = " << rds_mnt*1.0e6 << std::endl;
      std::cerr << "\tWAVELENGTH = " << bnd_ctr*1.0e6 << std::endl;
      std::cerr << "\tCORE SIZE PARAMETER = " << sz_prm_cor << "\tCOAT SIZE PARAMETER = " << sz_prm_mnt << std::endl;
      std::cerr << "\tQSCA = " << q_sct << "\tQEXT =  " << q_ext << "\tQBACK = " << q_bck << "\n" << std::endl;
    } // endif self-test
    
  } // !coat_flg endif coated spheres 
  
  // fxm: Switch to anomalous diffraction when sz_prm > sz_prm_wrn? How to handle g?
  if(sz_prm > sz_prm_wrn && dbg_lvl >= dbg_scl) std::cerr << "WARNING: size parameter = " << sz_prm << " > sz_prm_wrn = " << sz_prm_wrn << ", consider finishing Anomalous Diffraction Theory (ADT) implementation" << std::endl;
  
  // Anomalous diffraction when warranted
  if(sz_prm > sz_prm_wrn && false){
    rcd+=adt_apx // [fnc] Anomalous Diffraction Theory approximation
      (sz_prm, // I [m m-1] Size parameter
       idx_rfr_ffc, // I [frc] Effective refractive index of particle
       idx_rfr_mdm, // I [frc] Refractive index of medium
       q_abs, // I/O [frc] Absorption efficiency
       q_ext, // I/O [frc] Extinction efficiency
       q_sct); // I/O [frc] Scattering efficiency
  } // endif anomalous diffraction approximation
  
  // Diagnostic information

  // Angle routine required to diagnose size-integrated phase function
  // Accumulate size-distribution-integrated phase function for each band and size
  rcd+=mie_ngl_BoH83_csz
    (q_sct, // I [frc] Scattering efficiency
     sz_prm, // I [m m-1] Size parameter
     ngl_nbr, // I [nbr] Angle number (number of angles is 2*ngl_nbr-1)
     &s1[0], // I [frc] Scalar amplitude scattering matrix element S1
     &s2[0], // I [frc] Scalar amplitude scattering matrix element S2
     ngl, // I [rdn] Scattering angle
     ngl_dlt, // I [rdn] Width of angle bin
     phz_fnc, // O [sr-1] Phase function
     plz); // O [frc] Degree of linear polarization
  
  if(slf_tst_flg){
    (void)std::fprintf(stdout,"Note: S11 normalization differs from BoH83 by 4*pi but polarization should be identical\n");
    (void)std::fprintf(stdout,"  ANGLE       S11             POL             S33             S34\n");
    for(ngl_idx=0;ngl_idx<2*ngl_nbr-1;ngl_idx++)
      (void)std::fprintf(stdout," %6.2f  %13.6e  %13.6e\n",ngl[ngl_idx]*180.0/mth::cst_M_PIl,phz_fnc[ngl_idx],plz[ngl_idx]);
  } // end self-test
  
  // Sanity check
  q_abs=q_ext-q_sct; // [frc] Absorption efficiency
  if(q_abs < 0.0){
    // Warning prints once if mie_flg is false
    std::cerr << "WARNING: q_abs = " << q_abs << ", re-setting to 0.0" << std::endl;
    q_sct=q_ext; // [frc] Scattering efficiency
    q_abs=0.0; // [frc] Absorption efficiency
  } // end if q_abs < 0
  
  // Free amplitude scattering matrices
  if(s1) delete[]s1; // [frc] Scalar amplitude scattering matrix element S1
  if(s2) delete[]s2; // [frc] Scalar amplitude scattering matrix element S1
  if(XMU) delete[]XMU; // [frc] Cosine polar scattering angle Wis79 p. 16 (9B)

  if(dbg_lvl >= dbg_sbr) dbg_prn(sbr_nm,"Exiting...");
  return rcd; // [rcd] Return code
} // end mie_prc()

int // O [enm] Return success code
mie_sph_BoH83 // [fnc] Mie solution for homogeneous spheres
(const double sz_prm, // I [m m-1] Size parameter
 const std::complex<double> &idx_rfr_rlt, // I [frc] Refractive index of particle (aerosol) relative to medium (air)
 const long ngl_nbr, // I [nbr] Angle number (number of angles is 2*ngl_nbr-1)
 const double * const ngl, // I [rdn] Scattering angle
 double &asm_prm, // O [frc] Asymmetry parameter
 double &bck_hms, // O [frc] Hemispheric backscatter
 double &q_bck, // O [frc] Backscattering efficiency
 double &q_ext, // O [frc] Extinction efficiency
 double &q_sct, // O [frc] Scattering efficiency
 std::complex<double> * const s1, // O [frc] Scalar amplitude scattering matrix element S1
 std::complex<double> * const s2) // O [frc] Scalar amplitude scattering matrix element S2
{
  /* Purpose: Calculate amplitude scattering matrix elements and 
     efficiencies for extinction, total scattering and backscattering 
     for a given size parameter and relative refractive index 
     Inputs: sz_prm is actually size parameter times index of refraction of medium */

  // Local
  int rcd(0); // [rcd] Return code
  std::complex<double> xi,xi0,xi1;
  std::complex<double> y2;
  /* 20010727 csz++: JGW used long double type for variables on following line
     Original BoH83 bhmie() routine does not use long doubles
     Long doubles cause g++ compiler with strict flags to emit warnings like 
     "no match for `complex<double> * long double &'"
     Since long doubles are non-standard, and may simply be typedef'd to double,
     we revert these variables to plain doubles.
     ISO C99 does require long double type, so we may
     re-promote these variables to long double in the future. */
  /*  long double psi,psi0,psi1,aPsi,aPsi0,aPsi1;*/
  double psi,psi0,psi1,aPsi,aPsi0,aPsi1;
  double chi,chi0,chi1;
  //  double ngl_dlt_scl;
  double fn,P,T,trm_nbr_dbl,y2_mds,trm_idx_dbl;
  double trm1,trm2;
  long trm_nbr_max,trm_nbr;
  long jj;
  double rn_dbl;
  long trm_idx; // [idx] Term index
  long ngl_idx; // [idx] Angle index

  // Check input parameters
  if(sz_prm < 0.0) return -1;
  if(ngl_nbr <= 1) return -1;

  // Preliminary calculations
  bck_hms=0.0;
  y2=sz_prm*idx_rfr_rlt;
  
  // trm_nbr_max is based on Wiscombe's formulation
  trm_nbr_dbl=sz_prm+4.0*std::pow(sz_prm,1.0/3.0)+2.0; // [nbr] Series terminates after trm_nbr terms
  // Series terminates after trm_nbr terms
  trm_nbr=static_cast<long>(trm_nbr_dbl); // [nbr] Number of terms in Ricatti-Bessel functions
  y2_mds=std::sqrt(norm(y2));
  using mth::cst_M_PIl; // (3.1415926535897932384626433832795029L) [frc] 3
  //ngl_dlt_scl=mth::cst_M_PIl/(2.0*(ngl_nbr-1)); // [rdn] Width of angle bin
  // fxm: Not sure why trm_nbr_max differs from trm_nbr
  trm_nbr_max=15+((trm_nbr_dbl > y2_mds) ? static_cast<long>(trm_nbr_dbl) : static_cast<long>(y2_mds));

  // Memory allocation
  double *mu=new double[ngl_nbr+1];
  double *pi0=new double[ngl_nbr+1];
  double *pi1=new double[ngl_nbr+1];
  double *pi=new double[ngl_nbr+1];
  double *tau=new double[ngl_nbr+1];
  double *theta=new double[ngl_nbr+1];
  std::complex<double> *D=new std::complex<double>[trm_nbr_max+1];
  std::complex<double> *an=new std::complex<double>[trm_nbr_max+1];
  std::complex<double> *bn=new std::complex<double>[trm_nbr_max+1];

  for(trm_idx=0;trm_idx<=trm_nbr_max;trm_idx++) D[trm_idx]=std::complex<double>(0.0,0.0);
  // fxm: Allow non-uniform angular grids for, e.g., Lobatto quadrature
  for(ngl_idx=1;ngl_idx<=ngl_nbr;ngl_idx++){ // NB: Loop index starts at 1
    //    theta[ngl_idx]=(ngl_idx-1)*ngl_dlt_scl; // [rdn] Polar scattering angle
    theta[ngl_idx]=ngl[ngl_idx-1]; // [rdn] Polar scattering angle
    mu[ngl_idx]=std::cos(theta[ngl_idx]); // [frc] Cosine polar scattering angle
  } // end loop over ngl
  // Logarithmic derivative D[ngl_idx] calculated by downward
  // recurrence beginning with initial value i at ngl_idx=trm_nbr_max
  const long idx_end_trm(trm_nbr_max-1);
  for(ngl_idx=1;ngl_idx<=idx_end_trm;ngl_idx++){ // NB: Loop index starts at 1
    // 20010728: rn takes integer values but store as double to avoid compiler warnings
    rn_dbl=static_cast<double>(trm_nbr_max-ngl_idx+1);
    D[trm_nbr_max-ngl_idx]=rn_dbl/y2-1.0/(D[trm_nbr_max-ngl_idx+1]+rn_dbl/y2);
  } // end loop over ngl
  for(ngl_idx=1;ngl_idx<=ngl_nbr;ngl_idx++){ // NB: Loop index starts at 1
    pi0[ngl_idx]=0.0;
    pi1[ngl_idx]=1.0;
  } // end loop over ngl
  const long idx_end_ngl(2*ngl_nbr-1);
  // fxm 20040127: s1, s1 are dimension 2*ngl_nbr but only have 2*ngl_nbr-1 initialized!!!
  for(ngl_idx=1;ngl_idx<=idx_end_ngl;ngl_idx++){ // NB: Loop index starts at 1
    s1[ngl_idx]=std::complex<double>(0.0,0.0); // [frc] Scalar amplitude scattering matrix element S1 BoH83 p. 63 (3.12)
    s2[ngl_idx]=std::complex<double>(0.0,0.0); // [frc] Scalar amplitude scattering matrix element S2 BoH83 p. 63 (3.12)
  } // end loop over ngl
  // Riccati-Bessel Functions with real argument x
  // calculated by upward recurrence
  psi0=std::cos(sz_prm);
  psi1=std::sin(sz_prm);
  chi0=-std::sin(sz_prm);
  chi1=std::cos(sz_prm);
  aPsi0=psi0;
  aPsi1=psi1;
  xi0=std::complex<double>(aPsi0,-chi0);
  xi1=std::complex<double>(aPsi1,-chi1);
  q_sct=0.0;
  asm_prm=0.0;
  trm_idx=1L;
  an[trm_idx]=std::complex<double>(0.0,0.0);
  bn[trm_idx]=std::complex<double>(0.0,0.0);
  do{ // begin while (trm_idx-1-trm_nbr < 0);
    fn=(2.0*trm_idx+1.0)/(trm_idx*(trm_idx+1.0));
    psi=(2.0*trm_idx-1.0)*psi1/sz_prm-psi0;
    aPsi=psi;
    chi=(2.0*trm_idx-1.0)*chi1/sz_prm-chi0;
    xi=std::complex<double>(aPsi,-chi);
    /* csz++
       Mixing std::complex<double> and long double is not quite kosher in GCC
       Next four lines */
    an[trm_idx]=((D[trm_idx]/idx_rfr_rlt)+(trm_idx/sz_prm))*aPsi-aPsi1;
    an[trm_idx]/=(((D[trm_idx]/idx_rfr_rlt)+(trm_idx/sz_prm))*xi-xi1);
    bn[trm_idx]=((D[trm_idx]*idx_rfr_rlt)+(trm_idx/sz_prm))*aPsi-aPsi1;
    bn[trm_idx]/=(((D[trm_idx]*idx_rfr_rlt)+(trm_idx/sz_prm))*xi-xi1);
    // csz--
    q_sct+=((2.0*trm_idx+1.0)*(norm(an[trm_idx])+norm(bn[trm_idx]))); // ThS99 p. 78 (3.27)
    // Asymmetry parameter
    if(trm_idx > 1L){
      trm_idx_dbl=static_cast<double>(trm_idx);
      trm1=(trm_idx_dbl-1.0)*(trm_idx_dbl+1.0)/trm_idx_dbl*real(an[trm_idx-1]*conj(an[trm_idx])+bn[trm_idx-1]*conj(bn[trm_idx]));
      trm2=(2.0*trm_idx_dbl-1.0)/(trm_idx_dbl*(trm_idx_dbl-1.0))*real(an[trm_idx-1]*conj(bn[trm_idx-1]));
      asm_prm+=(trm1+trm2);
    } // endif trm_idx > 1
    
    for(ngl_idx=1;ngl_idx<=ngl_nbr;ngl_idx++){ // NB: Loop index starts at 1
      jj=2*ngl_nbr-ngl_idx;
      pi[ngl_idx]=pi1[ngl_idx];
      tau[ngl_idx]=trm_idx*mu[ngl_idx]*pi[ngl_idx]-(trm_idx+1.0)*pi0[ngl_idx];
      // [frc] Scalar amplitude scattering matrix elements BoH83 p. 63 (3.12)
      s1[ngl_idx]+=(fn*(an[trm_idx]*pi[ngl_idx]+bn[trm_idx]*tau[ngl_idx]));
      s2[ngl_idx]+=(fn*(an[trm_idx]*tau[ngl_idx]+bn[trm_idx]*pi[ngl_idx]));
      /* 20000804: Cast to int to ensure fast std::pow() function
	 Unfortunately, AIX compiler has no fast int^int version of pow 
	 Thus casting to double appears to be best solution */
      P=std::pow(-1.0,static_cast<double>(trm_idx-1));
      T=std::pow(-1.0,static_cast<double>(trm_idx));
      if(ngl_idx != jj){
	// [frc] Scalar amplitude scattering matrix elements BoH83 p. 63 (3.12)
	s1[jj]+=(fn*(an[trm_idx]*pi[ngl_idx]*P+bn[trm_idx]*tau[ngl_idx]*T));
	s2[jj]+=(fn*(an[trm_idx]*tau[ngl_idx]*T+bn[trm_idx]*pi[ngl_idx]*P));
      } // endif
    } // end loop over ngl
    psi0=psi1;
    psi1=psi;
    aPsi1=psi1;
    chi0=chi1;
    chi1=chi;
    xi1=std::complex<double>(aPsi1,-chi1);
    trm_idx++;
    for(ngl_idx=1;ngl_idx<=ngl_nbr;ngl_idx++){ // NB: Loop index starts at 1
      pi1[ngl_idx]=(2.0*trm_idx-1.0)/(trm_idx-1.0)*mu[ngl_idx]*pi[ngl_idx]-trm_idx*pi0[ngl_idx]/(trm_idx-1.0);
      pi0[ngl_idx]=pi[ngl_idx];
    } // end loop over ngl
  }while(trm_idx-1-trm_nbr < 0);
  /* Evaluate optical efficiencies using different methods */
  // Scattering efficiency evaluated by series expansion
  q_sct*=2.0/(sz_prm*sz_prm); // ThS99 p. 78 (3.27)
  // Extinction efficiency evaluated by optical theorem
  q_ext=4.0/(sz_prm*sz_prm)*real(s1[1]);
  // Backscattering efficiency evaluated by optical theorem BoH83 p. 121 (4.82)
  q_bck=4.0/(sz_prm*sz_prm)*norm(s1[2*ngl_nbr-1]);
  // Asymmetry parameter evaluated by series expansion BoH83 p. 120
  asm_prm*=4.0/(sz_prm*sz_prm*q_sct);
  rcd+=mie_bck_hms_Chy73(sz_prm,trm_nbr,an,bn,bck_hms);
  
  // Free dynamically allocated memory
  if(D) delete []D;
  if(an) delete []an;
  if(bn) delete []bn;
  if(mu) delete []mu;
  if(pi) delete []pi;
  if(pi0) delete []pi0;
  if(pi1) delete []pi1;
  if(tau) delete []tau;
  if(theta) delete []theta;
  
  return rcd; // [rcd] Return code
} // end mie_sph_BoH83()

int // O [rcd] Return success code
mie_ngl_BoH83_csz // [fnc] Mie scattering matrix for homogeneous spheres
(const double q_sct, // I [frc] Scattering efficiency
 const double sz_prm, // I [m m-1] Size parameter
 const long ngl_nbr, // I [nbr] Angle number (number of angles is 2*ngl_nbr-1)
 const std::complex<double> * const s1, // I [frc] Scalar amplitude scattering matrix element S1
 const std::complex<double> * const s2, // I [frc] Scalar amplitude scattering matrix element S2
 const prc_cmp * const ngl, // I [rdn] Scattering angle
 const prc_cmp * const ngl_dlt, // I [rdn] Width of angle bin
 prc_cmp *phz_fnc, // O [sr-1] Phase function
 prc_cmp *plz) // O [frc] Degree of linear polarization
{
  /* Purpose: Same as mie_ngl_BoH83() only customized to return Mie properties 
     of greatest interest for archiving---phase function and polarization */
  int rcd(0); // [rcd] Return code
  double phz_fnc_BoH83_nrm; // [sr-1] Phase function with BoH83 normalization
  double plr_lnr; // [frc] Degree of linear polarization
  double s11; // [] Mueller matrix: Ratio of scattered irradiance to incident irradiance
  double s12; // [] Mueller matrix BoH83 p. 67
  double s33; // [] Mueller matrix 
  double s34; // [] Mueller matrix 
  long ngl_idx;
  
  // Vet input
  if(!s1 || !s2) err_prn("!s1 || !s2");
  if(ngl_nbr <= 1) err_prn("ngl_nbr <= 1");
  
  //using mth::cst_M_PIl; // (3.1415926535897932384626433832795029L) [frc] 3
  //const double ngl_dlt_scl=mth::cst_M_PIl/(2.0*(ngl_nbr-1)); // [rdn] Width of angle bin
  const double s11_fwd_rcp=2.0/(norm(s1[1])+norm(s2[1])); // [] Inverse of s11 in forward direction BoH83 p. 112 (4.77)
  // Define normalization for p, the phase function, as in BoH83 p. 384 (13.3)
  // phz_fnc_nrm=1.0/(cst_M_PIl*sz_prm*sz_prm*q_sct); // [sr-1] Phase function normalization factor
  /* Define p normalization as in HaT74, ThS99 definition 
     This is 4*pi times BoH83 definition p. 384 (13.3) 
     4*pi factor required so integral of p over sphere = 1 */
  const double phz_fnc_nrm=4.0/(sz_prm*sz_prm*q_sct); // [] Phase function normalization factor
  for(ngl_idx=1;ngl_idx<=2*ngl_nbr-1;ngl_idx++){ // NB: Loop index starts at 1
    /* BoH83 p. 65 has s11=0.5*(norm(s1)+norm(s2)+norm(s3)+norm(s4))
       BoH83 p. 112 (4.77) and actual BHMIE code has following s11: */
    s11=0.5*(norm(s1[ngl_idx])+norm(s2[ngl_idx])); // BoH83 p. 65 (3.16)
    /* BoH83 p. 65 has s12=0.5*(norm(s2)-norm(s1)+norm(s3)-norm(s3))
       BoH83 p. 112 (4.77) and actual BHMIE code has following s12: */
    s12=0.5*(norm(s2[ngl_idx])-norm(s1[ngl_idx])); // BoH83 p. 65 (3.16)
    /* BoH83 p. 65 has plz==sqrt(((s21*s21)+(s31*s31)+(s41*s41))/(s11*s11))
       U=0-->s31=0. V=0-->s41=0. Thus if both U, V are 0, then plz=s21/s11 
       When is s21 == -s12 ? When incident light is unpolarized...
       BoH83 p. 112 (4.78) and p. 382 show plz=-s12/s11 */
    plr_lnr=-s12/s11; // BoH83 p. 113 (4.78), p. 112 (4.77)
    // BoH83 p. 65 has s33=real(s1*conj(s2)+s3*conj(s4)))
    s33=real(s2[ngl_idx]*conj(s1[ngl_idx]))/s11; // BoH83 p. 65 (3.16), p. 112 (4.77)
    // BoH83 p. 65 has s34=imag(s2*conj(s1)+s4*conj(s3)))
    s34=imag(s2[ngl_idx]*conj(s1[ngl_idx]))/s11; // BoH83 p. 65 (3.16), p. 112 (4.77)
    /* BoH83 codes (arbitrarily) normalizes s11 to be unity in forward direction
       Although convenient, this is not correct normalization for phase function
       See BoH83 p. 383 (13.2) */
    phz_fnc_BoH83_nrm=s11*s11_fwd_rcp; // [sr-1] Phase function with BoH83 normalization
    
    /* Store for diagnostics and netCDF output 
       Note we assign to ngl_idx-1 in this loop */
    // ngl[ngl_idx-1]=ngl_dlt_scl*(ngl_idx-1); // [rdn] Scattering angle
    plz[ngl_idx-1]=plr_lnr; // [frc] Degree of linear polarization
    phz_fnc[ngl_idx-1]=s11*phz_fnc_nrm; // [sr-1] Phase function BoH83 p. 383 (13.2)
    // ngl_dlt[ngl_idx-1]=ngl_dlt_scl; // [rdn] Width of angle bin
  } // end loop over ngl

  return rcd; // [rcd] Return code
} // end mie_ngl_BoH83_csz()

int // O [rcd] Return success code
mie_bck_hms_Chy73 // [fnc] Mie solution for hemispheric backscatter
(const double sz_prm, // I [m m-1] Size parameter
 const long trm_nbr, // I [nbr] Number of terms
 const std::complex<double> * const an, // I [frc] an coefficients
 const std::complex<double> * const bn, // I [frc] bn coefficients
 double &bck_hms) // O [frc] Hemispheric backscatter
{
  /* Purpose: Calculate hemispheric backscattering
     Method is as per Chylek et al. (1973)
     Coded by Peter Damiano
     Translated to C++ by Jeff Wong */
  unsigned long idx_k;
  unsigned long idx_l;
  double A,B,C,D,E,F,EF;
  double rk,rl,sb1,sb2,sb3;
  double temp;

  // Vet input parameters
  if(trm_nbr < 1) return -1; 
  if(!an || !bn || !bck_hms) return -1;
  const unsigned int trm_nbr_unsgn(static_cast<unsigned int>(trm_nbr));

  // sb1 is first term in equation for sb:(2k+1)(|ak|^2+|bk|^2)
  sb1=0.0;
  sb2=0.0;
  sb3=0.0;
  for(idx_k=1u;idx_k<=trm_nbr_unsgn;idx_k++){
    rk=idx_k;
    sb1+=(2.0*rk+1.0)*(norm(an[idx_k])+norm(bn[idx_k]));
    /* Following loop computes the second term
       Only increment over even k and odd l for this part
       Indices start at 1,k-1 -> k and k -> k+1 here */
    if((idx_k % 2u) == 0u){
      for(idx_l=1u;idx_l<=trm_nbr_unsgn;idx_l+=2u){
	rl=idx_l;
	/* nr_lnfactrl2 computes ln(x!!) = log of x double factorial
	   GSL gsl_sf_lndoublefact GDT01 p. 51 
	   Prototype is double gsl_sf_lndoublefact(const unsigned int n) */
	A=gsl_sf_lndoublefact(static_cast<unsigned int>(idx_k)-1u);
	B=gsl_sf_lndoublefact(static_cast<unsigned int>(idx_l));
	C=gsl_sf_lndoublefact(static_cast<unsigned int>(idx_k));
	D=gsl_sf_lndoublefact(static_cast<unsigned int>(idx_l)-1u);
	/* nr_lnfactrl2 is deprecated Numerical Recipes ln(x!!) solver
	   A=nr_lnfactrl2(static_cast<int>(idx_k-1));
	   B=nr_lnfactrl2(static_cast<int>(idx_l));
	   C=nr_lnfactrl2(static_cast<int>(idx_k));
	   D=nr_lnfactrl2(static_cast<int>(idx_l-1)); */
	E=(2.0*rk+1.0)*(2.0*rl+1.0)
	  /((rk-rl)*(rk+rl+1.0));
	F=A+B-C-D;
	EF=std::exp(F);
	// 20000804: fxm: change this to int^int
	temp=std::pow(-1.0,(rk+rl-1.0)/2.0)*E*EF;
	sb2 += temp*real(an[idx_k]*conj(an[idx_l])
			  +bn[idx_k]*conj(bn[idx_l]));
      } // end loop over idx_l
    } // endif
    /* Compute third term
       Increment over odd l and odd k only */
    if((idx_k % 2u) != 0u) {
      for(idx_l=1u;idx_l<=trm_nbr_unsgn;idx_l+=2u){
	rl=idx_l;
	/* nr_lnfactrl2 computes ln(x!!) = log of x double factorial
	   GSL gsl_sf_lndoublefact GDT01 p. 51 */
	A=gsl_sf_lndoublefact(idx_k);
	B=gsl_sf_lndoublefact(idx_l);
	C=gsl_sf_lndoublefact(idx_k-1u);
	D=gsl_sf_lndoublefact(idx_l-1u);
	/* nr_lnfactrl2 is deprecated Numerical Recipes ln(x!!) solver
	   A=nr_lnfactrl2(static_cast<int>(idx_k));
	   B=nr_lnfactrl2(static_cast<int>(idx_l));
	   C=nr_lnfactrl2(static_cast<int>(idx_k-1));
	   D=nr_lnfactrl2(static_cast<int>(idx_l-1)); */
	E=(2.0*rk+1.0)*(2.0*rl+1.0)
	 /(rk*(rk+1.0)*rl*(rl+1.0));
	F=A+B-C-D;
	EF=std::exp(F);
	// 20000804: fxm: change this to int^int
	sb3+=std::pow(-1.0,((rk+rl)/2.0))*E*EF
	 *real(an[idx_k]*conj(bn[idx_l]));
      } // end loop over idx_l
    } // endif
  } // end loop over idx_k
  bck_hms=sb1+2.0*sb2+2.0*sb3;
  bck_hms*=1.0/(sz_prm*sz_prm);

  return 0;
} // end mie_bck_hms_Chy73()

int // O [enm] Return success code
mie_sph_coat_BoH83 // [fnc] Mie solution for coated spheres
(const double sz_prm_cor, // I [m m-1] Size parameter of core
 const double sz_prm_mnt, // I [m m-1] Size parameter of mantle
 const std::complex<double> idx_rfr_rlt_cor, // I [frc] Core index of refraction relative to medium
 const std::complex<double> idx_rfr_rlt_mnt, // I [frc] Mantle index of refraction relative to medium
 double &q_ext, // O [frc] Extinction efficiency
 double &q_sct, // O [frc] Scattering efficiency
 double &q_bck, // O [frc] Backscattering efficiency
 double &asm_prm, // O [frc] Asymmetry parameter
 double &bck_hms) // O [frc] Hemispheric backscatter
{
  /* Subroutine mie_sph_coat_BoH83() calculates Q_ext, Q_sct, Q_bck for coated spheres
     All bessel functions computed by upward recurrence
     Input:
     sz_prm_cor=2*pi*rds_cor*idx_rfr_mdm/wvl
     sz_prm_mnt=2*pi*rds_mnt*idx_rfr_mdm/wvl
     idx_rfr_rlt_cor=idx_rfr_cor/idx_rfr_mdm
     idx_rfr_rlt_mnt=idx_rfr_mnt/idx_rfr_mdm 
     where 
     idx_rfr_cor = Complex refractive index of core
     idx_rfr_mnt = Complex refractive index of mantle
     idx_rfr_mdm = Real refractive index of medium
     rds_cor = Radius of core
     rds_mnt = Radius of mantle
     wvl = Wavelength of light in ambient medium
     
     Routine mie_sph_coat_BoH83() is bhcoat() from Bohren & Huffman (1983) 
     Fortran version obtained from website of Piotr Flatau 
     Original Fortran version from C. L. Joseph
     
     History:
     1992/11/24 (BTD) Explicit declaration of all variables
     1999/11/10 CSZ Translated to C++
     2000/02/20 CSZ Adding explicit calculation of asymmetry parameter and hemispheric backscattering
  */
  // Local variables:
  int rcd(0); // Return code
  long iflag;
  long trm_idx;
  long trm_nbr,trm_nbr_max;
  // qext and q_sct are for q_sct and q_ext but do not need de-referencing
  double chi0y,chi1y,chiy,psi0y,psi1y,psiy,qext,qsca,trm_nbr_dbl;
  double trm_idx_dbl; // [nbr] Term index double precision
  std::complex<double> amess1,amess2,amess3,amess4,an,ancap;
  std::complex<double> bn,bncap,brack;
  std::complex<double> chi0x2,chi0y2,chi1x2,chi1y2,chix2,chipx2,chipy2,chiy2,crack;
  std::complex<double> d0x1,d0x2,d0y2,d1x1,d1x2,d1y2,dnbar,gnbar;
  std::complex<double> refrel;
  std::complex<double> xback,xi0y,xi1y,xiy;
  std::complex<double> x1,x2,y2;
  double del(1.0e-8); // [frc] Inner sphere convergence criterion
  std::complex<double> ii(0.0e0,1.0e0);
  double y2_mds; // [] Modulus of y2

  // Terms in asymmetry parameter calculation
  double trm1; // [frc] Term in asymmetry parameter
  double trm2; // [frc] Term in asymmetry parameter

  // Initialize variables which accumulate
  bck_hms=0.0;
  asm_prm=0.0;

  // Main Code
  x1=idx_rfr_rlt_cor*sz_prm_cor;
  x2=idx_rfr_rlt_mnt*sz_prm_cor;
  y2=idx_rfr_rlt_mnt*sz_prm_mnt;
  trm_nbr_dbl=sz_prm_mnt+4.0*std::pow(sz_prm_mnt,1.0/3.0)+2.0; // [nbr] Series terminates after trm_nbr terms
  trm_nbr=static_cast<long>(trm_nbr_dbl); // [nbr] Series terminated after trm_nbr terms
  y2_mds=std::sqrt(norm(y2)); // [] Modulus of y2
  trm_nbr_max=15+((trm_nbr_dbl > y2_mds) ? static_cast<long>(trm_nbr_dbl) : static_cast<long>(y2_mds));

  // fxm: Not sure why trm_nbr_max differs from trm_nbr 
  std::complex<double> *an_arr=new std::complex<double>[trm_nbr_max+1];
  std::complex<double> *bn_arr=new std::complex<double>[trm_nbr_max+1];

  refrel=idx_rfr_rlt_mnt/idx_rfr_rlt_cor;
  d0x1=std::cos(x1)/std::sin(x1);
  d0x2=std::cos(x2)/std::sin(x2);
  d0y2=std::cos(y2)/std::sin(y2);
  psi0y=std::cos(sz_prm_mnt);
  psi1y=std::sin(sz_prm_mnt);
  chi0y=-std::sin(sz_prm_mnt);
  chi1y=std::cos(sz_prm_mnt);
  xi0y=psi0y-ii*chi0y;
  xi1y=psi1y-ii*chi1y;
  chi0y2=-std::sin(y2);
  chi1y2=std::cos(y2);
  chi0x2=-std::sin(x2);
  chi1x2=std::cos(x2);
  qsca=0.0;
  qext=0.0;
  xback=std::complex<double>(0.0,0.0);
  trm_idx=1L;
  iflag=0;
  do{ // begin do while(trm_idx-1L-trm_nbr < 0)
    trm_idx_dbl=static_cast<double>(trm_idx);
    psiy=(2.0*trm_idx_dbl-1.0)*psi1y/sz_prm_mnt-psi0y;
    chiy=(2.0*trm_idx_dbl-1.0)*chi1y/sz_prm_mnt-chi0y;
    xiy=psiy-ii*chiy;
    d1y2=1.0/(trm_idx_dbl/y2-d0y2)-trm_idx_dbl/y2;
    if(iflag == 1) goto goto_lbl;
    d1x1=1.0/(trm_idx_dbl/x1-d0x1)-trm_idx_dbl/x1;
    d1x2=1.0/(trm_idx_dbl/x2-d0x2)-trm_idx_dbl/x2;
    chix2=(2.0*trm_idx_dbl-1.0)*chi1x2/x2-chi0x2;
    chiy2=(2.0*trm_idx_dbl-1.0)*chi1y2/y2-chi0y2;
    chipx2=chi1x2-trm_idx_dbl*chix2/x2;
    chipy2=chi1y2-trm_idx_dbl*chiy2/y2;
    ancap=refrel*d1x1-d1x2;
    ancap/=(refrel*d1x1*chix2-chipx2);
    ancap/=(chix2*d1x2-chipx2);
    brack=ancap*(chiy2*d1y2-chipy2);
    bncap=refrel*d1x2-d1x1;
    bncap/=(refrel*chipx2-d1x1*chix2);
    bncap/=(chix2*d1x2-chipx2);
    crack=bncap*(chiy2*d1y2-chipy2);
    amess1=brack*chipy2;
    amess2=brack*chiy2;
    amess3=crack*chipy2;
    amess4=crack*chiy2;
    /* z = x + iy 
       sqrt(x*x+y*y) = |z| is called the complex modulus of z
       The complex modulus is always representable by a real number
       x^2+y^2 = |z|^2 is called the norm of z
       The norm is always representable by a real number
       Fortran abs(complex) intrinsic is equivalent to C++ sqrt(norm(complex)) 
       PTV96 p. 171 explains why sqrt(norm(complex)) is a bad thing!
       fxm: sqrt(norm()) should be replaced by intrinsic function */
    if(std::sqrt(norm(amess1)) > del*std::sqrt(norm(d1y2))) goto goto_lbl;
    if(std::sqrt(norm(amess2)) > del) goto goto_lbl;
    if(std::sqrt(norm(amess3)) > del*std::sqrt(norm(d1y2))) goto goto_lbl;
    if(std::sqrt(norm(amess4)) > del) goto goto_lbl;
    brack=std::complex<double>(0.0,0.0);
    crack=std::complex<double>(0.0,0.0);
    iflag=1;
  goto_lbl: dnbar=d1y2-brack*chipy2;
    dnbar/=(1.0-brack*chiy2);
    gnbar=d1y2-crack*chipy2;
    gnbar/=(1.0-crack*chiy2);
    an=(dnbar/idx_rfr_rlt_mnt+trm_idx_dbl/sz_prm_mnt)*psiy-psi1y;
    an/=((dnbar/idx_rfr_rlt_mnt+trm_idx_dbl/sz_prm_mnt)*xiy-xi1y);
    bn=(idx_rfr_rlt_mnt*gnbar+trm_idx_dbl/sz_prm_mnt)*psiy-psi1y;
    bn/=((idx_rfr_rlt_mnt*gnbar+trm_idx_dbl/sz_prm_mnt)*xiy-xi1y);
    qsca+=(2.0*trm_idx_dbl+1.0)*(norm(an)+norm(bn));

    /* Store complete an[] and bn[] arrays so asymmetry parameter and hemispheric backscatter may be computed analytically
       Analytical form of asymmetry parameter mixes terms of index n with n-1
       This requires keeping at least two terms of an[] and bn[] 
       We store entire an[] and bn[] arrays for simplicity
       Arrays of an[] and bn[] also required for analytical solution for hemispheric backscatter
       fxm: Should compute s1 and s2 scattering matrix elements so phase function can be computed */

    an_arr[trm_idx]=an;
    bn_arr[trm_idx]=bn;
    // Asymmetry parameter BoH83 p. 120
    if(trm_idx > 1){
      trm_idx_dbl=static_cast<double>(trm_idx);
      trm1=(trm_idx_dbl-1.0)*(trm_idx_dbl+1.0)/trm_idx_dbl*real(an_arr[trm_idx-1]*conj(an_arr[trm_idx])+bn_arr[trm_idx-1]*conj(bn_arr[trm_idx]));
      trm2=(2.0*trm_idx_dbl-1.0)/(trm_idx_dbl*(trm_idx_dbl-1.0))*real(an_arr[trm_idx-1]*conj(bn_arr[trm_idx-1]));
      asm_prm+=(trm1+trm2);
    } // endif trm_idx > 1
    
    xback+=(2.0*trm_idx_dbl+1.0)*std::pow(-1.0,static_cast<double>(trm_idx))*(an-bn);
    qext+=(2.0*trm_idx_dbl+1.0)*(real(an)+real(bn));
    psi0y=psi1y;
    psi1y=psiy;
    chi0y=chi1y;
    chi1y=chiy;
    xi1y=psi1y-ii*chi1y;
    chi0x2=chi1x2;
    chi1x2=chix2;
    chi0y2=chi1y2;
    chi1y2=chiy2;
    d0x1=d1x1;
    d0x2=d1x2;
    d0y2=d1y2;
    trm_idx++;
  }while(trm_idx-1-trm_nbr < 0);   // end do while(trm_idx-1-trm_nbr < 0)
  q_sct=(2.0/(sz_prm_mnt*sz_prm_mnt))*qsca;
  q_ext=(2.0/(sz_prm_mnt*sz_prm_mnt))*qext;
  q_bck=norm(xback);
  q_bck*=(1.0/(sz_prm_mnt*sz_prm_mnt));
  asm_prm*=4.0/(sz_prm_mnt*sz_prm_mnt*q_sct);
  rcd+=mie_bck_hms_Chy73(sz_prm_mnt,trm_nbr,an_arr,bn_arr,bck_hms);

  return rcd;
} // end mie_sph_coat_BoH83()

int // O [enm] Return success code
adt_apx // [fnc] Anomalous Diffraction Theory approximation
(const double &sz_prm, // I [m m-1] Size parameter
 const std::complex<prc_cmp> &idx_rfr_ffc, // I [frc] Effective refractive index of particle
 const std::complex<prc_cmp> &idx_rfr_mdm, // I [frc] Refractive index of medium
 double &q_abs, // I/O [frc] Absorption efficiency
 double &q_ext, // I/O [frc] Extinction efficiency
 double &q_sct) // I/O [frc] Scattering efficiency
{
  /* Purpose: Determine optical efficiencies using Anomalous Diffraction Theory (ADT) approximation
     This routine never really worked and needs re-implementation
     ADT descriptions in Van57 p. 179, Ste94, BoH83

     mie --dbg=3 --sz_nbr=1 --wvl_mnm=0.99 --wvl_mxm=1.01 --sz_mnm=99.99 --sz_mxm=100.01 --idx_rfr_prt="1.33+1.0e-8" ${DATA}/mie/mie.nc > ~/mie/foo 2>&1
  */
  const std::string sbr_nm("adt_apx"); // [sng] Name of subroutine
  const std::complex<prc_cmp> i_cpx(0.0,1.0); // [frc] Imaginary number
  int rcd(0); // [enm] Return code
  prc_cmp beta;
  prc_cmp cpv_foo2;
  prc_cmp cpv_foo;
  prc_cmp phz_dly_ctr; // [frc] Phase delay center
  prc_cmp q_ext_adt_Van57_non_abs_apx; // [frc] Non-absorbing aerosol approximation
  prc_cmp q_abs_adt_BoH83; // [frc] Absorption efficiency, anomalous diffraction theory, BoH83
  prc_cmp q_ext_adt_BoH83; // [frc] Extinction efficiency, anomalous diffraction theory, BoH83
  prc_cmp q_sct_adt_BoH83; // [frc] Scattering efficiency, anomalous diffraction theory, BoH83
  prc_cmp q_ext_adt_Van57; // [frc] Extinction efficiency, anomalous diffraction theory, Van57
  prc_cmp rho; // [frc] Phase delay center
  std::complex<prc_cmp> exp_mns_phz_dly_ctr_i; // [frc] Exponential of rho_i
  std::complex<prc_cmp> KKK; // [frc] Area-integrated anomalous phase delay center
  std::complex<prc_cmp> phz_dly_ctr_i; // [frc] Phase delay center times i
  std::complex<prc_cmp> rho_i; // [frc] Phase delay center times i
  
  /* fxmadt1: Is rho intendeded to be a complex parameter?
     Or, as a geometric parameter, must it be real?
     fxmadt2: eqn:phz_dly_ctr shows make assumption n_mdm=1 before defining rho 
     Must we make that exact assumption? */

  /* Implement BoH83(?) algorithm for extinction efficiency */
  phz_dly_ctr=2.0*sz_prm*(idx_rfr_ffc.real()-idx_rfr_mdm.real()); // [frc] Phase delay center
  phz_dly_ctr_i=phz_dly_ctr*i_cpx; // [frc] Phase delay center times i_cpx
  exp_mns_phz_dly_ctr_i=std::exp(-phz_dly_ctr_i); // [frc] Exponential of negative phz_dly_ctr_i
  KKK=PRC_CMP(0.5)+exp_mns_phz_dly_ctr_i/phz_dly_ctr_i+(exp_mns_phz_dly_ctr_i-PRC_CMP(1.0))/(phz_dly_ctr_i*phz_dly_ctr_i); // [frc] Area-integrated anomalous phase delay center

  q_ext_adt_BoH83=4.0*KKK.real(); // [frc] Extinction efficiency, ADT

  /* Implement BoH83(?) algorithm for absorption efficiency */
  if(idx_rfr_ffc.imag() != 0.0){
    // Absorbing approximation
    prc_cmp fct_b; // [frc] Four times size parameter times imaginary index of refraction (Lio02)
    prc_cmp exp_mns_fct_b; // [frc] Four times size parameter times imaginary index of refraction (Lio02)
    fct_b=4.0*sz_prm*idx_rfr_ffc.imag(); // [frc] Lio02 p. 101 (3.3.33)
    exp_mns_fct_b=std::exp(-fct_b); // [frc] Lio02 p. 101 (3.3.33)
    q_abs_adt_BoH83=1.0+(2.0/fct_b)*exp_mns_fct_b+2.0*(exp_mns_fct_b-1.0)/(fct_b*fct_b);
  }else{
    // Non-absorbing approximation
    q_abs_adt_BoH83=0.0;
  } // endif

  // Sanity check
  q_sct_adt_BoH83=q_ext_adt_BoH83-q_abs_adt_BoH83; // [frc] Scattering efficiency, ADT
  if(q_sct_adt_BoH83 < 0.0){
    // Warning prints once if mie_flg is false
    std::cerr << "WARNING: q_abs_adt_BoH83 = " << q_abs_adt_BoH83 << ", re-setting to 0.0" << std::endl;
    q_sct_adt_BoH83=0.0; // [frc] Scattering efficiency, ADT
    q_abs_adt_BoH83=q_ext_adt_BoH83; // [frc] Absorption efficiency, ADT
    err_prn(sbr_nm,"Exiting...");
  } // end if q_abs_adt_BoH83 < 0
  
  /* Implement Van57(?) */
  rho=2.0*sz_prm*(idx_rfr_ffc.real()-idx_rfr_mdm.real()); // [frc] Phase delay center
  beta=std::atan(idx_rfr_ffc.imag()/(idx_rfr_ffc.real()-real(idx_rfr_mdm)));
  cpv_foo=std::cos(beta)/rho;
  cpv_foo2=4.0*std::exp(-rho*std::tan(beta));
  /* Use anomalous diffraction theory to compute extinction efficiency
     Absorbing aerosol solution (also behaves correctly in non-absorbing limit): */
  q_ext_adt_Van57=2.0; // [frc] Extinction efficiency, ADT
  q_ext_adt_Van57+=-cpv_foo2*cpv_foo*std::sin(rho-beta); // [frc] Extinction efficiency, ADT
  q_ext_adt_Van57+=-cpv_foo2*cpv_foo*cpv_foo*std::cos(rho-2.0*beta); // [frc] Extinction efficiency, ADT
  q_ext_adt_Van57+=4.0*cpv_foo*cpv_foo*std::cos(2.0*beta); // [frc] Extinction efficiency, ADT
  
  // Non-absorbing aerosol approximation:
  q_ext_adt_Van57_non_abs_apx=2.0;
  q_ext_adt_Van57_non_abs_apx+=-4.0*std::sin(rho)/rho;
  q_ext_adt_Van57_non_abs_apx+=4.0*(1.0-std::cos(rho))/(rho*rho);

  // Print diagnostics
  if(dbg_lvl_get() == dbg_crr) std::cerr << "q_ext = " << q_ext << ", q_ext_adt_BoH83 = " << q_ext_adt_BoH83 << ", q_ext_adt_Van57 = " << q_ext_adt_Van57 << ", q_ext_adt_Van57_non_abs_apx = " << q_ext_adt_Van57_non_abs_apx << std::endl;
  if(dbg_lvl_get() == dbg_crr) std::cerr << "q_sct = " << q_sct << ", q_sct_adt_BoH83 = " << q_sct_adt_BoH83 << std::endl;
  if(dbg_lvl_get() == dbg_crr) std::cerr << "q_abs = " << q_abs << ", q_abs_adt_BoH83 = " << q_abs_adt_BoH83 << std::endl;
  
  return rcd;
} // end adt_apx()

