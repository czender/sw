// $Id$ 

// Implementation (declaration) of netCDF utilities

/* Copyright (C) 1997--2017 Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

/* File contains a grab bag of routines which are self-contained modules and perform their own netCDF I/O
   All these routines require NCO netCDF interfaces 
   Eventually all modules here except truly generic netCDF utilities like nco_var_dfn()
   should be moved into separate files or turned into classes so the
   I/O work can be handled by class methods. 
   nco.cc contains netCDF dependencies, so it is not appropriate to bundle nco.o into libcsz_c++
   Instead, nco.o is bundled into libcsm_c++
   Ideally, nco.cc should be completely generic to netCDF, and have no model-specific functions
   In the meantime, nco.cc is a convenient location to prototype "stand-alone" modules so it belongs in libcsm_c++
*/

#include <nco.hh> // netCDF utilities

// netCDF classes

// Friendly functions begin

std::ostream & // [srm] Reference to output stream for cascading
operator<< // [fnc] Stream insertion operator
(std::ostream &srm_out, // [srm] Output stream
 const gsl_complex &gsl_cpx) // [obj] Object to insert in stream
{
  /* Purpose: Overloaded stream insertion operator for gsl_complex objects
     This operator attempts to imitate the behavior of the default stream
     insertion operator implemented by g++ for std::complex<float> types.
     Overloaded stream operators discussed on DeD01 p. 529
     Usage: 
     std::cout << gsl_cpx;
  */
  srm_out << "(" << GSL_REAL(gsl_cpx) << "," << GSL_IMAG(gsl_cpx) << ")" << std::endl;
  // Following definition prints in real +/- imag i format
  //  srm_out << GSL_REAL(gsl_cpx) << (GSL_IMAG(gsl_cpx) > 0.0 ? "+" : "") << GSL_IMAG(gsl_cpx) << "i" << std::endl;

  return srm_out; // [srm] Reference to output stream for cascading
} // end operator<<()

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
aer_htg // [fnc] Determine aerosol heating characteristics
(const int &nc_out, // I [fl] netCDF file for output 
 const int &dmn_nbr_max, // I [nbr] Maximum number of dimensions allowed in single variable in output file
 const prc_cmp &dns_mdp, // I [kg m-3] Midlayer density
 const prc_cmp &prs_mdp, // I [Pa] Midlayer pressure
 const prc_cmp &tpt_mdp, // I [K] Midlayer temperature
 const prc_cmp &wnd_znl_mdp, // I [m s-1] Surface layer zonal wind speed
 const prc_cmp &dmt_dtc, // I [m] Diameter of detector
 const spc_slr_cls &flx_slr_src, // I [obj] Solar spectrum
 const aer_cls &aer, // I [obj] Aerosol
 const psd_cls *psd_lst, // I [obj] Particle size distribution
 const prc_cmp *abs_fsh, // I [frc] Absorption efficiency
 const prc_cmp *cnc, // I [# m-3] Number concentration 
 const prc_cmp *mss, // I [kg] Mass 
 const prc_cmp *rds_ctr, // I [m] Radius at bin center
 const prc_cmp *xsa, // I [m2] Cross-sectional area
 const prc_cmp *ss_co_alb_fsh, // I [frc] Single scattering co-albedo
 const long sz_nbr, // I [nbr] Number of size bins
 const prc_cmp &abs_cff_mss, // I [m2 kg-1] Mass absorption coefficient
 const prc_cmp &mss_rsl, // I [kg m-3] Mass concentration resolved
 const prc_cmp &ss_co_alb) // I [frc] Single scattering co-albedo
{
  /* Purpose: Determine aerosol heating characteristics
     nc_out is modified as variables are defined
     fxm: mss, cnc, sz_nbr should be containerized in psd objects/list
     fxm: psd_lst is currently unused
     fxm: Verify that all heating quantities scale correctly with boundary conditions
     Routine is template for module that works with class objects and performs own I/O
     Declare routine in nco.cc (rather than, say, aer.cc) to reduce number of files which depend on netCDF
  */

  /* Saltzman Laser: This routine is to validate Eric Saltzman's omega-ometer
   mie --no_wrn_ntp --ss_alb=0.9 --slr_spc=lsr --slr_cst=10000.0 --wnd_znl=1.0 --dmt_dtc=0.01 --dbg=1 --cmp_prt=sulfate --dist=lognormal --wvl_grd=reg --wvl_nbr=1 --wvl_mnm=0.495 --wvl_mxm=0.505 --bnd_nbr=1 --sz_grd=log --sz_mnm=0.01 --sz_mxm=1.0 --sz_nbr=200 --sz_dbg=0.05 --rds_nma=0.275 --gsd_anl=2.0 --cnc_nbr_anl=1.0e9 ${DATA}/mie/aer_sulfate_lsr.nc
   ncks -C -F -q -d sz,0.5e-6 -v htg_aer_prt,htg_aer_dgn,mss,cnc,slr_cst,abs_cff_mss,spc_heat_aer,mmw_aer,nrg_dps_vlm,tpt_dlt_aer_max,tpt_dlt_aer_ttl,cnc_nbr_rsl,mss_rsl,tpt_dlt_air,ss_co_alb,ss_co_alb_fsh,nrg_cnd_ttl,tpt_dff,tau_htg_prt ${DATA}/mie/aer_sulfate_lsr.nc
  */

  // Namespaces
  using phc::spc_heat_dry_air; // (1004.697) [J kg-1 K-1] IrG81 p. 25

  int rcd(0); // Return success code
  // Local
  bool flg_dry(true); // [flg] Particle is dry
  long idx; // [idx] Counting index
  const std::string sbr_nm("aer_htg"); // [sng] Name of subroutine
  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Entering...");
  
  const prc_cmp slr_cst(flx_slr_src.slr_cst_get()); // [W m-2] Solar constant
  const prc_cmp spc_heat_aer(aer.spc_heat_get()); // [J kg-1 K-1] Specific heat capacity
  const prc_cmp dff_H2O_air(dff_H2O_air_fst_scl(tpt_mdp,prs_mdp)); // [m2 s-1] Binary diffusivity of H2O in air
  // fxm: should use cnd_trm of moist air
  const prc_cmp cnd_trm_air(cnd_trm_dry_air_fst_scl(tpt_mdp)); // [W m-1 K-1] Thermal conductivity of air
  // fxm: dff_trm_air should use spc_heat of moist air
  const prc_cmp dff_trm_air(cnd_trm_air/(dns_mdp*spc_heat_dry_air)); // [m2 s-1] Thermal diffusivity of air PrK98 p. 507
  const prc_cmp vsc_dyn_atm(vsc_dyn_atm_fst_scl(tpt_mdp)); // [kg m-1 s-1] Dynamic viscosity of atmosphere
  const prc_cmp vsc_knm_atm(vsc_dyn_atm/dns_mdp); // [m2 s-1] Kinematic viscosity of atmosphere 
  const prc_cmp shm_nbr_trm(vsc_knm_atm/dff_trm_air); // [frc] Schmidt number for thermal diffusion in air PrK98 p. 541
  const prc_cmp shm_nbr_vpr(vsc_knm_atm/dff_H2O_air); // [frc] Schmidt number for vapor diffusion in air PrK98 p. 538

  /* Convolve aerosol mass concentration, absorbed radiant energy, and specific heat capacity to determine aerosol heating
     Aerosol heating is complicated
     PrK98 p. 542 (13-65) describes the balance of condensation and heating
     PrK78 p. 447 (13-67) does same
     RoY94 p. 103 covers the same physics
     SeP97 p. 685 discusses diffusional growth as one of several limiting growth processes, others being surface and volume reactions
     Aerosol heating is inherently a two dimensional problem since absorption depends on aerosol size and spectral flux
     The "laser problem" of Saltzman is to determine absorption properties (i.e., single scattering albedo) indirectly by measuring temperature change due to aerosol absorption of a laser beam
     Instruments are sensitive to temperature changes of O(0.01 K)
     For laser problem, our assumptions include:
     1. Spectral absorption is constant over domain of laser power (very reasonable)
     2. All aerosol absorption is immediately converted to thermal energy, none goes into chemical dissociation (also fine)
     4. Aerosol absorption is determined by the convolution of aerosol mass concentration, absorption efficiency, laser wavelength, and laser power (by definition)
     Absorbed thermal energy immediately changes temperature of aerosol particle
     3. If relative humidity is <~ 0.99 (depends on size), aerosol tends to evaporate
     3. Evaporation from aerosol surface takes place with timescale tau_1
     Solid particles do not evaporate and so must expel excess heat through conduction
     Conduction removes heat from aerosol with timescale tau_2 (tau_2 >> tau_1)
     
     The experimental setup is as follows:
     Powerful multi-channel laser operates at visible wavelengths of interest
     These wavelengths determine where single scattering albedo is measured
     Reasonable values would be, e.g., 0.5 um, 0.4 um, 0.63 um
     Collimated laser beam, order 1 square millimeter in diameter,
     is directed into cell orthogonal to flow of ambient air/aerosol.
     Laser is mounted across cell from an imaging thermal detector
     Thermal detector measures temperature in cell indirectly by measuring infrared radiation in some optimal bandpass to be determined.
     As laser is turned on and off, air temperature in cell changes due to aerosol heating. 
     Signal is difference between air temperature when laser is on and off

     Questions to address include:
     1. What is actually being measured? Air temperature changes or aerosol temperature changes?
     The aerosol heating is quickly diffused to the surrounding air
     The timescale for this conversion is determined in this subroutine
     2. How are the two related?
     3. Aerosol emission (IR) changes as temperature changes, how will this effect net heating? (very minor)
     4. What is timescale of heat conduction from aerosol to air?
     5. What is optimal bandpass for thermal detector? 
     Optimal bandpass should have linear radiation reponse to heating
  */
  prc_cmp *htg_aer_prt=new prc_cmp[sz_nbr]; // [K s-1] Aerosol heating rate per particle
  prc_cmp htg_aer_dgn; // [K s-1] Aerosol heating rate, total
  prc_cmp *nrg_dps_prt=new prc_cmp[sz_nbr]; // [J s-1] Energy deposition rate per particle
  prc_cmp nrg_dps_vlm(0.0); // [J s-1 m-3] Energy deposition rate
  prc_cmp nrg_dps_dgn; // [J s-1 m-3] Energy deposition rate, diagnostic
  prc_cmp *tpt_dlt_aer_max=new prc_cmp[sz_nbr]; // [K] Temperature change, maximum
  prc_cmp tpt_dlt_aer_ttl(0.0); // [K] Temperature change, total
  const prc_cmp tm_rrd(dmt_dtc/wnd_znl_mdp); // [s] Duration of irradiation
  // fxm: following computes absolute rates in each size bin, probably more kosher to compute spectral rates (per unit particle size), then multiply by bin width later
  for(idx=0;idx<sz_nbr;idx++){
    nrg_dps_prt[idx]=slr_cst*xsa[idx]*abs_fsh[idx]; // [J s-1] Energy deposition rate per particle
    nrg_dps_vlm+=cnc[idx]*nrg_dps_prt[idx]; // [J s-1 m-3] Energy deposition rate, volumetric
    htg_aer_prt[idx]=nrg_dps_prt[idx]/(spc_heat_aer*mss[idx]); // [K s-1] Aerosol heating rate per particle
    // Maximum temperature change assumes no evaporation or diffusion of heat
    tpt_dlt_aer_max[idx]=htg_aer_prt[idx]*tm_rrd; // [K] Temperature change, maximum
  } // end loop over sz
  nrg_dps_dgn=slr_cst*mss_rsl*abs_cff_mss; // [J s-1 m-3] Energy deposition rate, diagnostic
  htg_aer_dgn=slr_cst*mss_rsl*abs_cff_mss/(spc_heat_aer*mss_rsl); // [K s-1] Aerosol heating rate, total
  tpt_dlt_aer_ttl=htg_aer_dgn*tm_rrd; // [K] Temperature change, total

    /* Estimating timescale for heat conduction from aerosol to air:
       Relation between evaporating drop falling in subsaturated air is PrK98 p. 544 (13-65)
       Radiation plays same role (but opposite sign) as evaporation in particle heat balance
       Radiative absorption will maintain a constant temperature difference between particle and environment, tpt_dlt
       Quasi-equilibrium of tpt_dlt is reached when heat conduction from particle to air equals radiative heating
       Heat conduction increases linearly with tpt_dlt
       Thus, for dry particles, equilibrium tpt_dlt is determined by radiative absorption
       Inversely, a measured tpt_dlt can be inverted to obtain radiative absorption assuming heat conduction is known
       Of course heat conduction depends on size distribution so this is not trivial
       For wet particles, equilibrium tpt_dlt must account for both radiative absorption and evaporation
     */
  prc_cmp *nrg_cnd_prt=new prc_cmp[sz_nbr]; // [J s-1] Energy conduction rate per particle
  prc_cmp nrg_cnd_ttl(0.0); // [J s-1 m-3] Energy conduction rate, total
  prc_cmp *tau_htg_prt=new prc_cmp[sz_nbr]; // [s] Timescale for temperature difference to reach steady state
  prc_cmp *tpt_dff=new prc_cmp[sz_nbr]; // [K] Temperature difference, aerosol - environment
  prc_cmp *vnt_trm_aer=new prc_cmp[sz_nbr]; // [frc] Ventilation coefficient for heat transfer
  prc_cmp *vnt_mss_aer=new prc_cmp[sz_nbr]; // [frc] Ventilation coefficient for mass transfer
  prc_cmp cnd_trm_fct; // [J s-1 K-1] Factor in heat conductance
  using mth::cst_M_PIl; // (3.1415926535897932384626433832795029L) [frc] 3
  for(idx=0;idx<sz_nbr;idx++){
    // fxm: ventilation coefficients are WRONG!!!
    vnt_trm_aer[idx]=1.0; // [frc] Ventilation coefficient for heat transfer PrK98 p. 541 (13-60)
    vnt_mss_aer[idx]=1.0; // [frc] Ventilation coefficient for mass transfer PrK98 p. 541 (13-60)
    // Place heat conductance on per unit volume basis (i.e., multiply by cnc) for consistancy with nrg_dps_prt
    cnd_trm_fct=4.0*mth::cst_M_PIl*vnt_trm_aer[idx]*rds_ctr[idx]*cnd_trm_air; // [J s-1 K-1] Factor in heat conductance PrK98 p. 544 (13-65)
    // Equilibrium conditions for dry particles
    // NB: Studies are for radiant heat so use opposite convention of PrK98 and let tpt_dff represent aerosol-environment temperature which is positive definite for dry particles
    tpt_dff[idx]=nrg_dps_prt[idx]/cnd_trm_fct; // [K] Temperature difference, aerosol - environment PrK98 p. 544 (13-65)
    nrg_cnd_prt[idx]=cnd_trm_fct*tpt_dff[idx]; // [J s-1] Energy conduction rate per particle PrK98 p. 544 (13-65)
    nrg_cnd_ttl+=cnc[idx]*nrg_cnd_prt[idx]; // [J s-1 m-3] Energy conduction rate, total
    if(flg_dry){ // [flg] Particle is dry
      // Neglecting latent heating
      tau_htg_prt[idx]=mss[idx]*spc_heat_aer/cnd_trm_fct; // [s] Timescale for temperature difference to reach steady state PrK98 p. 544 (13-71)
    }else{
      prc_cmp dsvddt_H2O; // [kg m-3 K-1] Derivative of saturation vapor density over planar water
      dsvddt_H2O=dsvddt_H2O_PrK78_fst_scl(tpt_mdp); // [kg m-3 K-1] Derivative of saturation vapor density over planar liquid water
      err_prn(sbr_nm,"Wet heating physics not implemented yet...");
    } // endelse flg_dry
  } // end loop over sz

  /* Determine temperature change of air in limit of complete conduction from particles to air
     Aerosol mass mixing ratios are typically O(10^-9) (micrograms per kilogram)
     Thus we expect air temperature change to be O(10^-9) less than aerosol temperature change */ 
  prc_cmp htg_air; // [K s-1] Air heating rate
  prc_cmp tpt_dlt_air; // [K] Temperature change of air
  htg_air=nrg_dps_vlm/(dns_mdp*spc_heat_dry_air); // [K s-1] Air heating rate
  tpt_dlt_air=htg_air*tm_rrd; // [K] Temperature change of air

  // Sanity check: Sum of size resolved aerosol absorption should equal total energy deposition
  bool apx_eql; // [flg] Arguments are indistinguishable
  // fxm: use generic apx_eql_chk instead
  apx_eql=apx_eql_chk // [fnc] Determine whether arguments are indistinguishable
    (sbr_nm, // I [sng] Subroutine name of calling routine
     true, // I [flg] Verbose output
     nrg_dps_dgn, // I [frc] Target argument
     nrg_dps_vlm, // I [frc] Approximation to target argument
     PRC_CMP(1.0e-4), // I [frc] Relative precision
     "Vetting energy deposition"); // I [sng] Descriptive message of context

  // fxm: find way of inserting this to logfile output stream 
  if(true){
    std::cout << "Aerosol Heating:" << std::endl;
    // fxm: Use resolved concentration of multimodal distribution
    std::cout << "  Aerosol number concentration: " << psd_lst[0].cnc_nbr_anl_get() << " # m-3, Mass concentration = " << mss_rsl << " kg m-3" << std::endl;
    std::cout << "  Single scattering albedo = " << 1.0-ss_co_alb << ", Bulk mass absorption coefficient = " << abs_cff_mss << " m2 kg-1, Absorption optical depth to 1 m = " << mss_rsl*abs_cff_mss*1.0 << std::endl;
    std::cout << "  Absorption coefficient = " << mss_rsl*abs_cff_mss*1000.0 << " km-1, Extinction coefficient = " << (ss_co_alb == 0.0 ? 0.0 : (mss_rsl*abs_cff_mss*1000.0)/ss_co_alb) << " km-1" << std::endl;
    std::cout << "  Irradiance = " << slr_cst << " J m-2 s-1, Duration = " << tm_rrd << " s" << std::endl;
    std::cout << "  Energy deposition rate, volumetric = " << nrg_dps_vlm << " J m-3 s-1" << std::endl;
    std::cout << "  Aerosol heating rate (dry) = " << htg_aer_dgn << " K s-1, Air heating rate = " << htg_air << " K s-1" << std::endl;
    std::cout << "  Temperature change, aerosol = " << tpt_dlt_aer_ttl << " K, air = " << tpt_dlt_air << " K" << std::endl;
    std::cout << "  Timescale for steady state = " << tau_htg_prt[sz_nbr-1] << " s" << std::endl;
    std::cout << "  Steady state temperature difference = " << tpt_dff[sz_nbr-1] << " K" << std::endl;
  } // endif true

  // Output results of module
  const int sz_dmn(nco_inq_dimid(nc_out,static_cast<std::string>("sz"))); // [dmn] Size dimension
  const int *dmn_sz(&sz_dmn); // Pointer to size dimension
  const int *dmn_scl((int *)NULL); // [dmn] Dummy dimension for scalars CLIP

  var_mtd_sct var_mtd[]={
    {0,"htg_aer_prt",NC_FLOAT,1,dmn_sz,"long_name","Aerosol heating rate per particle","units","kelvin second-1"},
    {0,"htg_aer_dgn",NC_FLOAT,0,dmn_scl,"long_name","Aerosol heating rate, total","units","kelvin second-1"},
    {0,"htg_air",NC_FLOAT,0,dmn_scl,"long_name","Air heating rate","units","kelvin second-1"},
    {0,"nrg_cnd_prt",NC_FLOAT,1,dmn_sz,"long_name","Energy conduction rate per particle","units","joule second-1"},
    {0,"nrg_cnd_ttl",NC_FLOAT,0,dmn_scl,"long_name","Energy conduction rate, total","units","joule second-1 meter-3"},
    {0,"nrg_dps_prt",NC_FLOAT,1,dmn_sz,"long_name","Energy deposition rate per particle","units","joule second-1"},
    {0,"nrg_dps_vlm",NC_FLOAT,0,dmn_scl,"long_name","Energy deposition rate, total","units","joule second-1 meter-3"},
    {0,"tm_rrd",NC_FLOAT,0,dmn_scl,"long_name","Duration of irradiation","units","second"},
    {0,"tau_htg_prt",NC_FLOAT,1,dmn_sz,"long_name","Timescale for temperature difference to reach steady state","units","second"},
    {0,"tpt_dff",NC_FLOAT,1,dmn_sz,"long_name","Temperature difference, aerosol - environment","units","kelvin"},
    {0,"tpt_dlt_aer_max",NC_FLOAT,1,dmn_sz,"long_name","Temperature change, maximum","units","kelvin"},
    {0,"tpt_dlt_aer_ttl",NC_FLOAT,0,dmn_scl,"long_name","Temperature change, total","units","kelvin"},
    {0,"dff_H2O_air",NC_FLOAT,0,dmn_scl,"long_name","Binary diffusivity of H2O in air","units","meter2 second-1"},
    {0,"dff_trm_air",NC_FLOAT,0,dmn_scl,"long_name","Thermal diffusivity of air","units","meter2 second-1"},
    {0,"shm_nbr_trm",NC_FLOAT,0,dmn_scl,"long_name","Schmidt number for thermal diffusion in air","units","fraction"},
    {0,"shm_nbr_vpr",NC_FLOAT,0,dmn_scl,"long_name","Schmidt number for vapor diffusion in air","units","fraction"},
    {0,"vnt_mss_aer",NC_FLOAT,1,dmn_sz,"long_name","Ventilation coefficient for mass transfer","units","fraction"},
    {0,"vnt_trm_aer",NC_FLOAT,1,dmn_sz,"long_name","Ventilation coefficient for heat transfer","units","fraction"},
    {0,"cnd_trm_air",NC_FLOAT,0,dmn_scl,"long_name","Thermal conductivity of air","units","watt meter-1 kelvin-1"},
    {0,"tpt_dlt_air",NC_FLOAT,0,dmn_scl,"long_name","Temperature change, air","units","kelvin"}
  }; // end var_mtd_sct var_mtd[]
  const int var_mtd_nbr(sizeof(var_mtd)/sizeof(var_mtd_sct)); // [nbr] Number of variables in array
  rcd+=nco_var_dfn(nc_out,var_mtd,var_mtd_nbr,dmn_nbr_max); // [fnc] Define variables in output netCDF file

  // After writing, delete arrays 
  rcd=nco_put_var(nc_out,static_cast<std::string>("htg_aer_prt"),htg_aer_prt); delete []htg_aer_prt;
  rcd=nco_put_var(nc_out,static_cast<std::string>("htg_aer_dgn"),htg_aer_dgn);
  rcd=nco_put_var(nc_out,static_cast<std::string>("htg_air"),htg_air);
  rcd=nco_put_var(nc_out,static_cast<std::string>("nrg_cnd_prt"),nrg_cnd_prt); delete []nrg_cnd_prt;
  rcd=nco_put_var(nc_out,static_cast<std::string>("nrg_cnd_ttl"),nrg_cnd_ttl);
  rcd=nco_put_var(nc_out,static_cast<std::string>("nrg_dps_prt"),nrg_dps_prt); delete []nrg_dps_prt;
  rcd=nco_put_var(nc_out,static_cast<std::string>("nrg_dps_vlm"),nrg_dps_vlm);
  rcd=nco_put_var(nc_out,static_cast<std::string>("tm_rrd"),tm_rrd);
  rcd=nco_put_var(nc_out,static_cast<std::string>("tau_htg_prt"),tau_htg_prt); delete []tau_htg_prt;
  rcd=nco_put_var(nc_out,static_cast<std::string>("tpt_dff"),tpt_dff); delete []tpt_dff;
  rcd=nco_put_var(nc_out,static_cast<std::string>("tpt_dlt_aer_max"),tpt_dlt_aer_max); delete []tpt_dlt_aer_max;
  rcd=nco_put_var(nc_out,static_cast<std::string>("tpt_dlt_aer_ttl"),tpt_dlt_aer_ttl);
  rcd=nco_put_var(nc_out,static_cast<std::string>("dff_H2O_air"),dff_H2O_air);
  rcd=nco_put_var(nc_out,static_cast<std::string>("dff_trm_air"),dff_trm_air);
  rcd=nco_put_var(nc_out,static_cast<std::string>("shm_nbr_trm"),shm_nbr_trm);
  rcd=nco_put_var(nc_out,static_cast<std::string>("shm_nbr_vpr"),shm_nbr_vpr);
  rcd=nco_put_var(nc_out,static_cast<std::string>("vnt_mss_aer"),vnt_mss_aer); delete []vnt_mss_aer;
  rcd=nco_put_var(nc_out,static_cast<std::string>("vnt_trm_aer"),vnt_trm_aer); delete []vnt_trm_aer;
  rcd=nco_put_var(nc_out,static_cast<std::string>("cnd_trm_air"),cnd_trm_air);
  rcd=nco_put_var(nc_out,static_cast<std::string>("tpt_dlt_air"),tpt_dlt_air);

  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Exiting...");

  if(ss_co_alb_fsh[0] == ss_co_alb_fsh[0]){;} // CEWU Compiler Error Warning Usage

  return rcd; // [enm] Return success code
} // end aer_htg()

int // O [enm] Return success code
rnd_chm // [fnc] Raindrop chemistry
(const int &nc_out, // I [fl] netCDF file for output 
 const int &dmn_nbr_max, // I [nbr] Maximum number of dimensions allowed in single variable in output file
 const prc_cmp &tpt_mdp, // I [K] Midlayer temperature
 const long sz_nbr, // I [nbr] Number of size bins
 const prc_cmp *cnc, // I [# m-3] Number concentration 
 const prc_cmp *mss, // I [kg] Mass 
 const prc_cmp *rds_ctr, // I [m] Radius at bin center
 const prc_cmp &vmr_CO2) // [mlc mlc-1] Volume mixing ratio of CO2
{
  /* Purpose: Determine raindrop chemistry
     nc_out is modified as variables are defined */

  // Namespaces
  using phc::rxn_ntp_CO2_H2O_298K; // [J mol-1] Reaction enthalpy (heat of dissolution) at 298K inferred from E/R quoted by Sep97 p. 391 Tbl. 6.A.1
  using phc::tpt_298; // (298.00) [K] Standard temperature for Henry's Law
  using phc::gas_cst_unv; // (8.31441) [J mol-1 K-1] Universal gas constant
  using phc::cst_Avagadro; // (6.02214199e+23) [mlc mol-1] Avagadro's number (CODATA)
  using phc::dns_H2O_lqd_std; // (1000.0) [kg m-3] Density of liquid water
  using phc::mmw_H2O; // (1.8015259e-02) [kg mol-1] Mean molecular weight of H2O HITRAN96

  int rcd(0); // Return success code
  // Local
  long idx; // [idx] Counting index
  const std::string sbr_nm("rnd_chm"); // [sng] Name of subroutine
  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Entering...");
  
  const prc_cmp cff_hnr_CO2_H2O_298K(3.4e-2); // [mol ltr-1 atm-1] Henry's Law coefficient for CO2(g)<-->H2O(l) equilibrium at 298K SeP97 p. 341 Table 6.2
  prc_cmp cff_hnr_CO2_H2O; // [mol ltr-1 atm-1] Henry's Law coefficient for CO2(g)<-->H2O(l) equilibrium
  cff_hnr_CO2_H2O=cff_hnr_CO2_H2O_298K*std::exp(rxn_ntp_CO2_H2O_298K*(1.0/tpt_298-1.0/tpt_mdp)/gas_cst_unv); // [mol ltr-1 atm-1] SeP97 p. 342 (6.5)

  const prc_cmp eqm_cst_H2O_298K(1.0e-14); // [mol2 ltr-2] Equilibrium constant for H2O(l)<-->H+(l) + OH-(l) at 298K SeP97 p. 345 (6.12)
  // fxm: Account for temperature dependence of K_H2O
  const prc_cmp eqm_cst_H2O(eqm_cst_H2O_298K); // [mol2 ltr-2] Equilibrium constant for H2O(l)<-->H+(l) + OH-(l) SeP97 p. 345 (6.12)

  // fxm: Account for temperature dependence
  const prc_cmp eqm_cst_CO2H2O_HppHCO3m_298K(4.3e-7); // [mol ltr-1] Equilibrium constant for CO2*H2O<-->H+ + HCO3- at 298K SeP97 p. 345 Table (6.4), p. 346 (6.18)
  const prc_cmp eqm_cst_HCO3m_HppCO3mm_298K(4.7e-11); // [mol ltr-1] Equilibrium constant for HCO3-<-->H+ + CO3-- at 298K SeP97 p. 345 Table (6.4), p. 346 (6.19)

  // Compute coefficients for cubic equation in H+: [H+]^3 + a[H+]^2 + b[H+] + c = 0
  const prc_cmp prm_a(0.0); // [mol ltr-1] SeP97 p. 348 (6.29)
  const prc_cmp prm_b(-(eqm_cst_H2O+vmr_CO2*cff_hnr_CO2_H2O*eqm_cst_CO2H2O_HppHCO3m_298K)); // [mol2 ltr-2] SeP97 p. 348 (6.29)
  const prc_cmp prm_c(-2.0*vmr_CO2*cff_hnr_CO2_H2O*eqm_cst_CO2H2O_HppHCO3m_298K*eqm_cst_HCO3m_HppCO3mm_298K); // [mol3 ltr-3] SeP97 p. 348 (6.29)

  prc_cmp *mss_C_aqs=new prc_cmp[sz_nbr]; // [kg] Aqueous mass of C per particle
  prc_cmp cnc_mss_C_aqs_ttl(0.0); // [kg m-3] Mass concentration of aqueous C

  std::complex<prc_cmp> sln_1; // [frc] First root of complex equation
  std::complex<prc_cmp> sln_2; // [frc] Second root of complex equation
  std::complex<prc_cmp> sln_3; // [frc] Third root of complex equation

  // fxm: cubic solver in mth.hh is broken
  rcd+=eqn_cbc_slvr(prm_a,prm_b,prm_c,&sln_1,&sln_2,&sln_3);

  prc_cmp rnd_pH(-std::log10(sln_1.real())); // [frc] Raindrop pH
  prc_cmp rnd_cnc_Hp_mol_ltr(sln_1.real()); // [M] = [mol l-1] Raindrop H+ concentration
  prc_cmp rnd_cnc_Hp_vlm(rnd_cnc_Hp_mol_ltr*cst_Avagadro*1.0e6); // [# m-3] Raindrop H+ volume concentration
  prc_cmp rnd_cnc_Hp_mss(rnd_cnc_Hp_mol_ltr*cst_Avagadro*1.0e6*dns_H2O_lqd_std); // [# kg-1] Raindrop H+ mass concentration
  prc_cmp rnd_mmr_Hp(rnd_cnc_Hp_mol_ltr*cst_Avagadro*1.0e6*dns_H2O_lqd_std*mmw_H2O); // [kg kg-1] Raindrop H+ mass mixing ratio
  if(dbg_lvl_get() == dbg_crr){
    // fxm: find way of inserting this to logfile output stream 
    std::cout << "Raindrop chemistry:" << std::endl;
    std::cout << prg_nm_get() << ": "+sbr_nm+"() calling eqn_cbc_slvr() with prm_a = " << prm_a << ", prm_b = " << prm_b << ", prm_c = " << prm_c << std::endl;
    std::cout << "Solution(s) to z^3 + " << prm_a << "*z^2 + " << prm_b << "*z + " << prm_c << " = 0 are z1 = " << sln_1 << ", z2 = " << sln_2 << ", z3 = " << sln_3 << std::endl;
    std::cout << "  CO2 vmr = " << vmr_CO2*1.0e6 << " ppm" << std::endl;
    std::cout << "  Raindrop pH = " << rnd_pH << std::endl;
    std::cout << "  H+ concentration = " << rnd_cnc_Hp_mol_ltr << " mol l-1 = " << rnd_cnc_Hp_vlm << " # m-3 = " << rnd_cnc_Hp_mss << " # kg-1 = " << rnd_mmr_Hp << " kg kg-1" << std::endl;
  } // endif dbg

  for(idx=0;idx<sz_nbr;idx++){
    mss_C_aqs[idx]=0.0; // [kg] Aqueous mass of C per particle
    cnc_mss_C_aqs_ttl+=mss_C_aqs[idx]; // [kg m-3] Mass concentration of aqueous C
  } // end loop over sz

  // Output module results
  const int sz_dmn(nco_inq_dimid(nc_out,static_cast<std::string>("sz"))); // [dmn] Size dimension
  const int *dmn_sz(&sz_dmn); // Pointer to size dimension
  const int *dmn_scl((int *)NULL); // [dmn] Dummy dimension for scalars CLIP

  var_mtd_sct var_mtd[]={
    {0,"cnc_mss_C_aqs_ttl",NC_FLOAT,0,dmn_scl,"long_name","Mass concentration of aqueous C","units","kilogram meter-3"},
    {0,"mss_C_aqs",NC_FLOAT,1,dmn_sz,"long_name","Aqueous mass of C per particle","units","kilogram"},
    {0,"rnd_cnc_Hp_mol_ltr",NC_FLOAT,0,dmn_scl,"long_name","Raindrop H+ concentration","units","mole liter-1"},
    {0,"rnd_cnc_Hp_mss",NC_FLOAT,0,dmn_scl,"long_name","Raindrop H+ mass concentration","units","kilogram-1"},
    {0,"rnd_cnc_Hp_vlm",NC_FLOAT,0,dmn_scl,"long_name","Raindrop H+ volume concentration","units","meter-3"},
    {0,"rnd_mmr_Hp",NC_FLOAT,0,dmn_scl,"long_name","Raindrop H+ mass mixing ratio","units","kilogram kilogram-1"},
    {0,"rnd_pH",NC_FLOAT,0,dmn_scl,"long_name","Raindrop pH","units","fraction"},
    {0,"vmr_CO2",NC_FLOAT,0,dmn_scl,"long_name","Volume mixing ratio of CO2","units","molecule molecule-1"}
  }; // end var_mtd_sct var_mtd[]
  const int var_mtd_nbr(sizeof(var_mtd)/sizeof(var_mtd_sct)); // [nbr] Number of variables in array
  rcd+=nco_var_dfn(nc_out,var_mtd,var_mtd_nbr,dmn_nbr_max); // [fnc] Define variables in output netCDF file

  // After writing, delete arrays 
  rcd=nco_put_var(nc_out,static_cast<std::string>("rnd_cnc_Hp_mol_ltr"),rnd_cnc_Hp_mol_ltr);
  rcd=nco_put_var(nc_out,static_cast<std::string>("rnd_cnc_Hp_vlm"),rnd_cnc_Hp_vlm);
  rcd=nco_put_var(nc_out,static_cast<std::string>("rnd_cnc_Hp_mss"),rnd_cnc_Hp_mss);
  rcd=nco_put_var(nc_out,static_cast<std::string>("rnd_mmr_Hp"),rnd_mmr_Hp);
  rcd=nco_put_var(nc_out,static_cast<std::string>("rnd_pH"),rnd_pH);
  rcd=nco_put_var(nc_out,static_cast<std::string>("cnc_mss_C_aqs_ttl"),cnc_mss_C_aqs_ttl);
  rcd=nco_put_var(nc_out,static_cast<std::string>("mss_C_aqs"),mss_C_aqs); delete []mss_C_aqs;
  rcd=nco_put_var(nc_out,static_cast<std::string>("vmr_CO2"),vmr_CO2);

  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Exiting...");

  if(cnc[0] == (mss[0] == rds_ctr[0])){;} // CEWU Compiler Errror Warning Usage

  return rcd; // [enm] Return success code
} // end rnd_chm()

int // O [enm] Return success code
rfl_frs // [fnc] Fresnel reflectance
(const int &nc_out, // I [fl] netCDF file for output 
 const int &dmn_nbr_max, // I [nbr] Maximum number of dimensions allowed in single variable in output file
 const long wvl_nbr, // I [nbr] Number of wavelength bins
 const std::complex<prc_cmp> *idx_rfr_1, // I [frc] Refractive index of transmitted medium
 const std::complex<prc_cmp> *idx_rfr_2, // I [frc] Refractive index of incident medium
 const prc_cmp &slr_zen_ngl_cos) // I [frc] Cosine solar zenith angle
{
  /* Purpose: Determine Fresnel reflectance and transmittance
     nc_out is modified as variables are defined
     Fresnel reflectance and transmittance describe light propogation past an 
     interface separating two homogeneous substances with differing indices of 
     refraction. 
     Nomenclature of medium 1 and medium 2 follows convention of BoH83:
     Medium 1 = "transmitted medium" = "lower slab" = destination of radiation     
     Medium 1 defaults to air, indices of refraction from idx_rfr_mdm in mie()
     Medium 2 = "incident medium" = "upper slab" = source of radiation
     Medium 2 defaults to air, indices of refraction from idx_rfr_prt in mie()

     Definitions:
     Plane of incidence is contains ray of incident light and normal to surface
     Polarization components are parallel and perpindicular to this plane
     Parallel component travels along ray in plane of incdence 
     Perpindicular component travels along ray but orthogonal to plane of incdence 
     Angle of incidence is smallest angle from surface normal to indident ray
     Angle of reflection is smallest angle from surface normal to reflected ray
     Angle of reflection equals angle of incidence
     Angle of refraction is smallest angle from surface normal on other side of 
     surface to refracted ray
     Thus both angle of incidence and angle or refraction are less than pi/2
     Brewster angle is measured from surface normal like angle of incidence

     Routine gracefully handles most general case where absorption in 
     transmitting medium causes refracted angle to be complex.

     Simulate reflection at air/water interface:
     mie --idx_rfr_prt=1.01 --cmp_mdm=h2o_lqd
     ncks -C -F -q -d wvl,0.5e-6 -v ngl_inc_dgr,ngl_rfl_dgr,ngl_rfr_dgr,rfl_prl,rfl_prp,trn_prl,trn_prp,rfl_flx_frs,trn_flx_frs ${DATA}/mie/out.nc
  */

  int rcd(0); // Return success code
  // Local
  long idx; // [idx] Counting index
  const std::string sbr_nm("rfl_frs"); // [sng] Name of subroutine
  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Entering...");
  
  // Append suffix "_frs" to variables which might otherwise be confused with aerosol "slab" or "layer" properties in output file
  prc_cmp *trn_flx_frs=new prc_cmp[wvl_nbr]; // [frc] Flux transmittance of unpolarized light, Fresnel
  prc_cmp *rfl_flx_frs=new prc_cmp[wvl_nbr]; // [frc] Flux reflectance of unpolarized light, Fresnel

  // Compute angles
  prc_cmp ngl_inc(std::acos(slr_zen_ngl_cos)); // [rdn] Angle of incidence
  prc_cmp ngl_rfl(ngl_inc); // [rdn] Angle of reflection
  using mth::cst_M_PIl; // (3.1415926535897932384626433832795029L) [frc] 3
  prc_cmp ngl_inc_dgr(180.0*ngl_inc/mth::cst_M_PIl); // [dgr] Angle of incidence
  prc_cmp ngl_rfl_dgr(180.0*ngl_rfl/mth::cst_M_PIl); // [dgr] Angle of reflection
  prc_cmp cos_ngl_inc(std::cos(ngl_inc)); // [frc] Cosine of angle of incidence
  std::complex<prc_cmp> sin_ngl_inc(static_cast<std::complex<prc_cmp> >(std::sin(ngl_inc))); // [rdn] Sine of angle of incidence

  prc_cmp *ngl_brw=new prc_cmp[wvl_nbr]; // [rdn] Brewster angle
  prc_cmp *ngl_brw_dgr=new prc_cmp[wvl_nbr]; // [dgr] Brewster angle
  std::complex<prc_cmp> *ngl_rfr=new std::complex<prc_cmp>[wvl_nbr]; // [rdn] Complex angle of refraction
  prc_cmp *ngl_rfr_dgr=new prc_cmp[wvl_nbr]; // [dgr] Angle of refraction, real component
  std::complex<prc_cmp> *idx_rfr_rlt=new std::complex<prc_cmp>[wvl_nbr]; // [frc] Relative refractive index
  std::complex<prc_cmp> arc_sin_arg; // [frc] Argument to arcsin function
  gsl_complex ngl_rfr_gsl; // [rdn] Complex angle of refraction
  gsl_complex arc_sin_arg_gsl; // [frc] Argument to arcsin function

  for(idx=0;idx<wvl_nbr;idx++){
    // Relative index of refraction is transmitted index over incident index
    assert(idx_rfr_2[idx] != static_cast<std::complex<prc_cmp> >(0.0));
    idx_rfr_rlt[idx]=idx_rfr_1[idx]/idx_rfr_2[idx]; // [frc] Relative refractive index
    /* If transmitting medium is absorbing, refracting angle is complex and does not
       equal, in any simple sense, geometric angle of refraction.
       If transmitting medium is non-absorbing, refracting angle is real and equals 
       geometric angle of refraction. 
       Computing complex angle of refraction requires taking complex arcsine of complex argument, which is not as easy in C as in Fortran */
    arc_sin_arg=sin_ngl_inc/idx_rfr_rlt[idx]; // [frc] Argument to arcsin function
    arc_sin_arg_gsl=gsl_complex_rect(arc_sin_arg.real(),arc_sin_arg.imag()); // [frc] Argument to arcsin function
    ngl_rfr_gsl=gsl_complex_arcsin(arc_sin_arg_gsl); // [frc] Complex angle of refraction
    ngl_rfr[idx]=std::complex<prc_cmp>(GSL_REAL(ngl_rfr_gsl),GSL_IMAG(ngl_rfr_gsl)); // [rdn] Complex angle of refraction
    // Archive only real part of complex angle of refraction
    ngl_rfr_dgr[idx]=180.0*ngl_rfr[idx].real()/mth::cst_M_PIl; // [dgr] Angle of refraction, real component
    ngl_brw[idx]=std::atan(idx_rfr_rlt[idx].real()); // [rdn] Brewster angle BoH83 p. 36 (2.71.5)
    ngl_brw_dgr[idx]=180.0*ngl_brw[idx]/mth::cst_M_PIl; // [dgr] Brewster angle
  } // end loop over wvl

  /* Compute reflectance properties for all wavelengths for each polarization state
     Use relative index of refraction to save some operations
     Reflectance and transmission are complex quantities */
  std::complex<prc_cmp> *rfl_prl=new std::complex<prc_cmp>[wvl_nbr]; // [frc] Reflectance of electric field amplitude parallel to plane of incidence
  std::complex<prc_cmp> *rfl_prp=new std::complex<prc_cmp>[wvl_nbr]; // [frc] Reflectance of electric field amplitude orthogonal to plane of incidence
  std::complex<prc_cmp> *trn_prl=new std::complex<prc_cmp>[wvl_nbr]; // [frc] Transmittance of electric field amplitude parallel to plane of incidence
  std::complex<prc_cmp> *trn_prp=new std::complex<prc_cmp>[wvl_nbr]; // [frc] Transmittance of electric field amplitude orthogonal to plane of incidence
  gsl_complex cos_ngl_rfr_gsl; // [frc] Cosine of angle of refraction
  std::complex<prc_cmp> cos_ngl_rfr; // [frc] Cosine of angle of refraction
  gsl_complex ngl_rfr_crr_gsl; // [frc] Current angle of refraction
  for(idx=0;idx<wvl_nbr;idx++){
    // Fresnel formulae BoH83 p. 34, ThS99 p. 508
    // fxm: 20010726 Following line "should" work but broke in gsl-0.9
    //    cos_ngl_rfr_gsl=gsl_complex_cos(gsl_complex_rect(ngl_rfr[idx].real(),ngl_rfr[idx].imag())); // [frc] Cosine of angle of refraction
    ngl_rfr_crr_gsl=gsl_complex_rect(ngl_rfr[idx].real(),ngl_rfr[idx].imag()); // [frc] Current angle of refraction
    cos_ngl_rfr_gsl=gsl_complex_cos(ngl_rfr_crr_gsl); // [frc] Cosine of angle of refraction
    assert(gsl_complex_abs(cos_ngl_rfr_gsl) <= 1.0);
    cos_ngl_rfr=std::complex<prc_cmp>(GSL_REAL(cos_ngl_rfr_gsl),GSL_IMAG(cos_ngl_rfr_gsl)); // [frc] Cosine of angle of refraction
    assert(norm(cos_ngl_rfr) <= 1.0);
    rfl_prl[idx]=(cos_ngl_rfr-idx_rfr_rlt[idx]*cos_ngl_inc)/(cos_ngl_rfr+idx_rfr_rlt[idx]*cos_ngl_inc); // [frc] Reflectance of electric field amplitude parallel to plane of incidence BoH83 p. 34 (2.67), ThS99 p. 508 (E.17.5)
    rfl_prp[idx]=(cos_ngl_inc-idx_rfr_rlt[idx]*cos_ngl_rfr)/(cos_ngl_inc+idx_rfr_rlt[idx]*cos_ngl_rfr); // [frc] Reflectance of electric field amplitude orthogonal to plane of incidence BoH83 p. 35 (2.69), ThS99 p. 508 (E.17.5)
    trn_prl[idx]=PRC_CMP(2.0)*cos_ngl_inc/(cos_ngl_rfr+idx_rfr_rlt[idx]*cos_ngl_inc); // [frc] Transmittance of electric field amplitude parallel to plane of incidence BoH83 p. 35 (2.68), ThS99 p. 508 (E.17.5)
    trn_prp[idx]=PRC_CMP(2.0)*cos_ngl_inc/(cos_ngl_inc+idx_rfr_rlt[idx]*cos_ngl_rfr); // [frc] Transmittance of electric field amplitude orthogonal to plane of incidence BoH83 p. 35 (2.70), ThS99 p. 508 (E.17.5)
  } // end loop over wvl

  // Compute flux reflectance and flux transmittance
  prc_cmp *rfl_flx_prl=new prc_cmp[wvl_nbr]; // [frc] Flux reflectance of parallel polarization state
  prc_cmp *rfl_flx_prp=new prc_cmp[wvl_nbr]; // [frc] Flux reflectance of perpindicular polarization state
  prc_cmp *trn_flx_prl=new prc_cmp[wvl_nbr]; // [frc] Flux transmittance of parallel polarization state
  prc_cmp *trn_flx_prp=new prc_cmp[wvl_nbr]; // [frc] Flux transmittance of perpindicular polarization state
  std::complex<prc_cmp> trn_flx_prl_cpx; // [frc] Flux transmittance of parallel polarization state
  std::complex<prc_cmp> trn_flx_prp_cpx; // [frc] Flux transmittance of perpindicular polarization state
  for(idx=0;idx<wvl_nbr;idx++){
    /* Flux reflectance or plane albedo is ratio of reflected to incident irradiance
       Incident irradiance is m_i norm(A) cos(theta_i) / (2 magnetic_permeability cst_spd_lgt_vcm)
       Reflected irradiance is m_i norm(R) cos(theta_r) / (2 magnetic_permeability cst_spd_lgt_vcm)
       Since theta_i = theta_r, ratio of reflected to incident irradiance is 
       norm(R)/norm(A)
       Flux reflectance is norm of reflectance of electric field amplitudes */
    rfl_flx_prl[idx]=norm(rfl_prl[idx]); // [frc] Flux reflectance of parallel polarization state
    rfl_flx_prp[idx]=norm(rfl_prp[idx]); // [frc] Flux reflectance of perpindicular polarization state
    /* Flux transmittance is more, well, complex
       Flux transmittance is ratio of transmitted irradiance to incident irradiance
       Transmitted irradiance is m_t norm(T) cos(theta_t) / (2 magnetic_permeability cst_spd_lgt_vcm)
       Incident irradiance is m_i norm(A) cos(theta_i) / (2 magnetic_permeability cst_spd_lgt_vcm)
       Ratio of transmitted to incident irradiance is
       m_t norm(T) cos(theta_t) / (m_i norm(A) cos(theta_i))
       Thus flux transmittance is complex if m_t is complex
       We archive only real part of flux transmittance
     */
    // fxm: same factors are used twice, could be precomputed or moved into one loop
    // fxm: 20010726 Following line "should" work but broke in gsl-0.9
    //   cos_ngl_rfr_gsl=gsl_complex_cos(gsl_complex_rect(ngl_rfr[idx].real(),ngl_rfr[idx].imag())); // [frc] Cosine of angle of refraction
    ngl_rfr_crr_gsl=gsl_complex_rect(ngl_rfr[idx].real(),ngl_rfr[idx].imag()); // [frc] Current angle of refraction
    cos_ngl_rfr_gsl=gsl_complex_cos(ngl_rfr_crr_gsl); // [frc] Cosine of angle of refraction
    cos_ngl_rfr=std::complex<prc_cmp>(GSL_REAL(cos_ngl_rfr_gsl),GSL_IMAG(cos_ngl_rfr_gsl)); // [frc] Cosine of angle of refraction
    trn_flx_prl_cpx=idx_rfr_rlt[idx]*norm(trn_prl[idx])*cos_ngl_rfr/cos_ngl_inc; // [frc] Flux transmittance of parallel polarization state
    trn_flx_prp_cpx=idx_rfr_rlt[idx]*norm(trn_prp[idx])*cos_ngl_rfr/cos_ngl_inc; // [frc] Flux transmittance of perpindicular polarization state
    trn_flx_prl[idx]=trn_flx_prl_cpx.real(); // [frc] Flux transmittance of parallel polarization state
    trn_flx_prp[idx]=trn_flx_prp_cpx.real(); // [frc] Flux transmittance of perpindicular polarization state
  } // end loop over wvl

  // Sanity check
  // In order to conserve energy, flux transmittance should be one minus corresponding flux reflectance

  // Flux reflectance or transmittance of unpolarized light is mean of polarized quantities
  for(idx=0;idx<wvl_nbr;idx++){
    rfl_flx_frs[idx]=0.5*(rfl_flx_prl[idx]+rfl_flx_prp[idx]); // [frc] Flux reflectance of unpolarized light, Fresnel ThS99 p. 509
    trn_flx_frs[idx]=0.5*(trn_flx_prl[idx]+trn_flx_prp[idx]); // [frc] Flux transmittance of unpolarized light, Fresnel ThS99 p. 509
  } // end loop over wvl

  if(true){
    std::cout << "Fresnel reflectance:" << std::endl;
    std::cout << "  Refractive indices: incident medium = " << idx_rfr_2[0] << ", transmitted medium = " << idx_rfr_1[0] << std::endl;
    std::cout << "  ngl_inc_dgr = " << ngl_inc_dgr <<  ", ngl_rfl_dgr = " << ngl_rfl_dgr << ", ngl_rfr_dgr = " << ngl_rfr[0]*PRC_CMP(180.0)/static_cast<prc_cmp>(mth::cst_M_PIl) << ", ngl_brw_dgr = " << ngl_brw_dgr[0]<< std::endl;
    std::cout << "  trn_prl = " << trn_prl[0] << ", trn_prp = " << trn_prp[0] << std::endl;
    std::cout << "  rfl_prl = " << rfl_prl[0] << ", rfl_prp = " << rfl_prp[0] << std::endl;
    std::cout << "  trn_flx_prl = " << trn_flx_prl[0] << ", trn_flx_prp = " << trn_flx_prp[0] << std::endl;
    std::cout << "  rfl_flx_prl = " << rfl_flx_prl[0] << ", rfl_flx_prp = " << rfl_flx_prp[0] << std::endl;
    std::cout << "  rfl_flx_frs = " << rfl_flx_frs[0] << ", trn_flx_frs = " << trn_flx_frs[0] << std::endl;
    std::cout << "  Energy conservation: rfl_flx_frs + trn_flx_frs = " << rfl_flx_frs[0]+trn_flx_frs[0] << std::endl;
  } // endif true

  // Sanity check
  if(ngl_inc == 0.0){
    /* Reflection of normal incident light is independent of polarization and may
       be written as simple function of relative index of refraction */
    prc_cmp rfl_flx_nrm; // [frc] Flux reflectance of normally incident light
    bool apx_eql; // [flg] Arguments are indistinguishable
    for(idx=0;idx<wvl_nbr;idx++){
      rfl_flx_nrm=norm((PRC_CMP(1.0)-idx_rfr_rlt[idx])/(PRC_CMP(1.0)+idx_rfr_rlt[idx])); // [frc] Flux reflectance of normally incident light BoH83 p. 32 (2.58)
      // fxm: use generic apx_eql_chk instead
      apx_eql=apx_eql_chk // [fnc] Determine whether arguments are indistinguishable
	(sbr_nm, // I [sng] Subroutine name of calling routine
	 true, // I [flg] Verbose output
	 rfl_flx_nrm, // I [frc] Target argument
	 rfl_flx_frs[idx], // I [frc] Approximation to target argument
	 PRC_CMP(1.0e-4), // I [frc] Relative precision
	 "Vetting flux reflectance"); // I [sng] Descriptive message of context
    } // end loop over wvl
  } // endif

  // Delete obsolete arrays
  delete []idx_rfr_rlt; // [frc] Relative refractive index
  delete []ngl_rfr; // [rdn] Angle of refraction, real component
  delete []ngl_brw; // [rdn] Brewster angle
  delete []rfl_prl; // [frc] Reflectance of electric field amplitude parallel to plane of incidence
  delete []rfl_prp; // [frc] Reflectance of electric field amplitude orthogonal to plane of incidence
  delete []trn_prl; // [frc] Transmittance of electric field amplitude parallel to plane of incidence
  delete []trn_prp; // [frc] Transmittance of electric field amplitude orthogonal to plane of incidence

  // Output results of module
  const int wvl_dmn(nco_inq_dimid(nc_out,static_cast<std::string>("wvl"))); // [dmn] Wavelength dimension
  const int *dmn_wvl(&wvl_dmn); // Pointer to wavelength dimension
  const int *dmn_scl((int *)NULL); // [dmn] Dummy dimension for scalars CLIP

  var_mtd_sct var_mtd[]={
    {0,"ngl_inc_dgr",NC_FLOAT,0,dmn_scl,"long_name","Angle of incidence","units","degree"},
    {0,"ngl_rfl_dgr",NC_FLOAT,0,dmn_scl,"long_name","Angle of reflection","units","degree"},
    {0,"ngl_rfr_dgr",NC_FLOAT,1,dmn_wvl,"long_name","Angle of refraction, real component","units","degree"},
    {0,"ngl_brw_dgr",NC_FLOAT,1,dmn_wvl,"long_name","Angle of refraction, real component","units","degree"},
    {0,"rfl_flx_frs",NC_FLOAT,1,dmn_wvl,"long_name","Flux reflectance of unpolarized light, Fresnel","units","fraction"},
    {0,"rfl_flx_prl",NC_FLOAT,1,dmn_wvl,"long_name","Flux reflectance of parallel polarization state","units","fraction"},
    {0,"rfl_flx_prp",NC_FLOAT,1,dmn_wvl,"long_name","Flux reflectance of perpindicular polarization state","units","fraction"},
    {0,"trn_flx_frs",NC_FLOAT,1,dmn_wvl,"long_name","Flux transmittance of unpolarized light, Fresnel","units","fraction"},
    {0,"trn_flx_prl",NC_FLOAT,1,dmn_wvl,"long_name","Flux transmittance of parallel polarization state","units","fraction"},
    {0,"trn_flx_prp",NC_FLOAT,1,dmn_wvl,"long_name","Flux transmittance of perpindicular polarization state","units","fraction"}
  }; // end var_mtd_sct var_mtd[]
  const int var_mtd_nbr(sizeof(var_mtd)/sizeof(var_mtd_sct)); // [nbr] Number of variables in array
  rcd+=nco_var_dfn(nc_out,var_mtd,var_mtd_nbr,dmn_nbr_max); // [fnc] Define variables in output netCDF file

  // After writing, delete arrays 
  rcd=nco_put_var(nc_out,static_cast<std::string>("ngl_rfr_dgr"),ngl_rfr_dgr); delete []ngl_rfr_dgr;
  rcd=nco_put_var(nc_out,static_cast<std::string>("ngl_brw_dgr"),ngl_brw_dgr); delete []ngl_brw_dgr;
  rcd=nco_put_var(nc_out,static_cast<std::string>("ngl_inc_dgr"),ngl_inc_dgr);
  rcd=nco_put_var(nc_out,static_cast<std::string>("ngl_rfl_dgr"),ngl_rfl_dgr);
  rcd=nco_put_var(nc_out,static_cast<std::string>("rfl_flx_prl"),rfl_flx_prl); delete []rfl_flx_prl;
  rcd=nco_put_var(nc_out,static_cast<std::string>("rfl_flx_prp"),rfl_flx_prp); delete []rfl_flx_prp;
  rcd=nco_put_var(nc_out,static_cast<std::string>("trn_flx_prl"),trn_flx_prl); delete []trn_flx_prl;
  rcd=nco_put_var(nc_out,static_cast<std::string>("trn_flx_prp"),trn_flx_prp); delete []trn_flx_prp;
  rcd=nco_put_var(nc_out,static_cast<std::string>("trn_flx_frs"),trn_flx_frs); delete []trn_flx_frs;
  rcd=nco_put_var(nc_out,static_cast<std::string>("rfl_flx_frs"),rfl_flx_frs); delete []rfl_flx_frs;

  if(dbg_lvl_get() >= dbg_sbr) dbg_prn(sbr_nm,"Exiting...");
  return rcd; // [enm] Return success code
} // end rfl_frs()

// Global functions with C++ linkages end
