// $Id$ 

// Purpose: Physical constants for C++ in SI units

/* Copyright (C) 1997--present Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

/* Source: Based on 19980726 ${HOME}/f/parameter.com, v. 1.4
   Try to keep phys_cst.hh and parameter.com synchronized
   Values should be consistent with 1998 CODATA adjustment reported in MoT00
   Where necessary for compatibility with models like CCM, less precise values
   may be used in place of CODATA, but this should be noted in the source data 
   for the constant 

   CODATA physical constants are available from 
   NIST Constants, Units and Uncertainty: http://www.physics.nist.gov/cuu/Constants
   Constants taken from this site are denote (CODATA) in source attribution */

/* Usage: Defines physical constants as "const double" rather than preprocessor tokens
   These constants often have global (program) scope
#include <phys_cst.hh> // Physical constants NB: Live code, not a header
using namespace phc; // [nms] Physical constant namespace */

/* NB: Defining mean molecular weight in SI units [kg mol-1] instead of [g mol-1]
   keeps Avagadro's number the same and still gets rid of all those pesky factors of 1000.0 */

// High precision mmw's come from HITRAN isotopic distribution in hitran.com

#ifndef PHYS_CST_HH
#define PHYS_CST_HH

namespace phc{ // [nms] Physical constant namespace

// Fundamental constants

  // 2018 SI NIST definition refers to
  // doi:10.1088/1681-7575/aa950a (see also https://www.nist.gov/si-redefinition/meet-constants)
  //  const double cst_Avagadro(6.02214199e+23); // (6.02214199e+23) [mlc mol-1] Avagadro's number
  //const double cst_Avagadro(6.022045e+23); // (6.022045e+23) [mlc mol-1] Avagadro's number (CCM?)
  const double cst_Avagadro(6.02214076e+23); // (6.02214076e+23) [mlc mol-1] Avagadro's number (2022 CODATA) 2018 SI NIST
  // const double cst_Boltzmann(1.3806503e-23); // (1.3806503e-23) [J K-1] Boltzmann's constant MoT00 p. BG11
  //const double cst_Boltzmann(1.38063e-23); // (1.38063e-23) [J K-1] Boltzmann's constant (CCM?)
  const double cst_Boltzmann(1.380649e-23); // (1.380649e-23) [J K-1] Boltzmann's constant (2018 SI NIST, 2022 CODATA 2019, exact)
  // const double cst_Planck(6.62620e-34); // (6.62620e-34) [J s] Planck's constant (CCM?)
  //  const double cst_Planck(6.62606876e-34); // (6.62606876e-34) [J s] Planck's constant (CODATA, 2018 SI NIST) (exact)
  const double cst_Planck(6.62607015e-34); // (6.62607015e-34) [J s] Planck's constant (2022 CODATA) (exact)
  //  const double speed_of_light(2.99793e+08); // (2.99792458e+08) [m s-1] Speed of light in vacuo (???)
  const double speed_of_light(2.99792458e+08); // (2.99792458e+08) [m s-1] Speed of light in vacuo (2022 CODATA, 2018 SI NIST)

  //  const double cst_Gravitation(6.673e-11); // (6.673e-11) [N m2 kg-2] Universal gravitational constant (CODATA)
  const double cst_Gravitation(6.67430e-11); // (6.67430e-11) [N m2 kg-2] Universal gravitational constant (2022 CODATA)
  const double cff_hnr_HNO3_H2O_298K(2.1e5); // (2.1e5) [mol ltr-1 atm-1] Henry's Law coefficient of HNO3 in liquid water at 298K SeP97 p. 341 Table 6.2
  const double cst_Loschmidt(2.6871e+25); // (2.6871e+25) [mlc m-3] Loschmidt's number (molecules of air per cubic meter at STP) (NB: this is derivable in principle)
  // const double cst_Stefan_Boltzmann_GoY89(5.67e-8); // (5.67e-8) [W m-2 K-4] Stefan-Boltzmann constant E3SM used this until at least 20250508
  // const double cst_Stefan_Boltzmann_GoY89(5.67032e-8); // (5.67032e-8) [W m-2 K-4] Stefan-Boltzmann constant GoY89 p. 462 NB: I used this value until 20250505
  // const double cst_Stefan_Boltzmann(5.670374419e-8); // (5.670374419e-8) [W m-2 K-4] Stefan-Boltzmann constant CODATA exact https://physics.nist.gov/cgi-bin/cuu/Value?sigma
  const double cst_Stefan_Boltzmann(5.670374419184431e-8); // (5.670374419184431e-8) [W m-2 K-4] Stefan-Boltzmann constant CODATA exact double precision from ncap2 -O -s 'M_PI=3.14159265358979323846264338327950288;cst_Boltzmann=1.380649e-23;cst_Planck=6.62607015e-34;speed_of_light=2.99792458e+08;cst_Stefan_Boltzmann=2*M_PI^5*cst_Boltzmann^4/(15*speed_of_light^2*cst_Planck^3);print(cst_Stefan_Boltzmann);' ~/foo.nc

  const double cst_von_krm(0.4); // (0.4) [frc] Von Karman's constant
  const double dmt_cll_CO2(3.34e-10); // (3.34e-10) [m] Collision diameter of CO2 CRC95 p. 6-244
  const double dmt_cll_HNO3(3.5e-10); // (3.5e-10) [m] fxm: Pure guess
  const double dmt_cll_N2(3.15e-10); // (3.15e-10) [m] Collision diameter of N2 CRC95 p. 6-244 
  const double dmt_cll_O2(2.98e-10); // (2.98e-10) [m] Collision diameter of O2 CRC95 p. 6-244
  const double dmt_cll_air(3.46e-10); // (3.46e-10) [m] Mean collision diameter of air SeP97 p. 1292 Table A.7
  const double gas_cst_unv(8.314472); // (8.314472) [J mol-1 K-1] Universal gas constant (CODATA)
  //const double gas_cst_unv(8.31441); // (8.31441) [J mol-1 K-1] Universal gas constant (Used until 20250601)
  const double grv_sfc_mean(9.80665); // (9.80665) [m s-2] Mean gravitational acceleration at Earth's surface (RRTMGP)
  const double hc(1.986488377e-25); // (1.986488377e-25) [J m] Planck constant times speed of light = hc
  const double hc2(5.9553531e-17); // (5.9553531e-17) [J m2 s-1] Planck constant times speed of light squared = hc2
  const double joules_per_calorie(4.1855); // (4.1855) [J cal-1] Cal = energy to heat 1g H20 1C @ 15C
  const double joules_per_eV(1.60217733e-19); // (1.60217733e-19) [J ev-1] Joules per electron volt CRC95 inside back cover
  const double k_O2_O2(2.6e-29); // (2.6e-29) [m3 mlc-1] Equilibrium rate constant for O2 + O2 <-> O2-O2 (Sha77 p. 436, 527)
  const double ltn_heat_fsn_H2O_0C(0.334e6); // (0.334e6) [J kg-1] Latent heat of fusion of H2O at 0 C Tre922 p. xxix
  const double ltn_heat_fsn_H2O_old(0.3336e06); // (0.3336e06) [J kg-1] Latent heat of fusion of H2O at 0 C, standard CCM:lsm/phyconi.F  
  const double ltn_heat_fsn_H2O_std(0.3337e06); // (0.3336e06) [J kg-1] Latent heat of fusion of H2O at 0 C, standard CCSM3:shr_const_latice
  const double ltn_heat_sbl_H2O_std(2.8440e06); // (2.8440e06) [J kg-1] Latent heat of sublimation of H2O, standard CCM:lsm/phyconi.F
  const double ltn_heat_vpr_H2O_100C(2.25e6); // (2.25e6) [J kg-1] Latent heat of vaporization of H2O at 100 C Tre922 p. xxix
  const double ltn_heat_vpr_H2O_old(2.5104e06); // (2.5104e06) [J kg-1] Latent heat of vaporization of H2O, standard CCM:lsm/phyconi.F
  const double ltn_heat_vpr_H2O_std(2.501e06); // (2.501e06) [J kg-1] Latent heat of vaporization of H2O, standard CCSM3:shr_const_latvap
  const double mmw_Al(26.981539e-03); // (26.981539e-03) [kg mol-1] Mean molecular weight of Al IUPAC
  const double mmw_C(12.011e-03); // (12.011e-03) [kg mol-1] Mean molecular weight of C IUPAC
  const double mmw_CO2(4.4009743e-02); // (4.4009743e-02) [kg mol-1] Mean molecular weight of CO2 HITRAN96
  const double mmw_Ca(40.078e-03); // (40.078e-03) [kg mol-1] Mean molecular weight of Ca IUPAC
  const double mmw_Cl(35.4527e-03); // (35.4527e-03) [kg mol-1] Mean molecular weight of Cl IUPAC
  const double mmw_Fe(55.847e-03); // (55.847e-03) [kg mol-1] Mean molecular weight of Fe IUPAC
  const double mmw_H(1.00794e-03); // (1.00794e-03) [kg mol-1] Mean molecular weight of H IUPAC
  //  const double mmw_H2O(1.8016-02); // (1.8016e-02) [kg mol-1] Mean molecular weight of H2O 2018 SI NIST
  const double mmw_H2O(1.8015259e-02); // (1.8015259e-02) [kg mol-1] Mean molecular weight of H2O HITRAN96
  const double mmw_H2O2(3.4005480e-02); // (3.4005480e-02) [kg mol-1] Mean molecular weight of H2O2 HITRAN96
  const double mmw_H2OH2O(36.03e-03); // (36.03e-3) [kg mol-1] Mean molecular weight of H2OH2O
  const double mmw_HNO3(6.2995644e-02); // (6.2995644e-02) [kg mol-1] Mean molecular weight of HNO3 HITRAN96
  const double mmw_K(39.0983e-03); // (39.0983e-03) [kg mol-1] Mean molecular weight of K IUPAC
  const double mmw_Mg(24.3050e-03); // (24.3050e-03) [kg mol-1] Mean molecular weight of Mg IUPAC
  const double mmw_Mn(54.93805e-03); // (54.93805e-03) [kg mol-1] Mean molecular weight of Mn IUPAC
  const double mmw_N(14.00674e-03); // (14.00674e-03) [kg mol-1] Mean molecular weight of N IUPAC
  const double mmw_N2(2.8006147e-02); // (2.8006147e-02) [kg mol-1] Mean molecular weight of N2 HITRAN96
  const double mmw_N2O(4.4012674e-02); // (4.4012674e-02) [kg mol-1] Mean molecular weight of N2O HITRAN96
  const double mmw_NH3(1.7030201e-02); // (1.7030201e-02) [kg mol-1] Mean molecular weight of NH3 HITRAN96
  const double mmw_NO(3.0005630e-02); // (3.0005630e-02) [kg mol-1] Mean molecular weight of NO HITRAN96
  const double mmw_NO2(4.5992904e-02); // (4.5992904e-02) [kg mol-1] Mean molecular weight of NO2 HITRAN96
  const double mmw_Na(22.989768e-03); // (22.989768e-03) [kg mol-1] Mean molecular weight of Na IUPAC
  const double mmw_O(15.9994e-03); // (15.9994e-03) [kg mol-1] Mean molecular weight of O IUPAC
  const double mmw_O2(3.1998575e-02); // (3.1998575e-02) [kg mol-1] Mean molecular weight of O2 HITRAN96
  const double mmw_O2N2(60.0e-03); // (60.0e-3) [kg mol-1] (this is a guess)
  const double mmw_O2O2(64.0e-03); // (64.0e-3) [kg mol-1] (this is a guess)
  const double mmw_O3(4.7997832e-02); // (4.7997832e-02) [kg mol-1 ] Mean molecular weight of O3 HITRAN96
  const double mmw_OH(1.7006907e-02); // (1.7006907e-02) [kg mol-1] Mean molecular weight of OH HITRAN96
  const double mmw_P(30.973762e-03); // (30.973762e-03) [kg mol-1] Mean molecular weight of P IUPAC
  const double mmw_S(32.066e-03); // (32.066e-03) [kg mol-1] Mean molecular weight of S IUPAC
  const double mmw_SO2(6.4046674e-02); // (6.4046674e-02) [kg mol-1] Mean molecular weight of SO2 HITRAN96
  const double mmw_Si(28.0855e-03); // (28.0855) [kg mol-1] Mean molecular weight of Si IUPAC
  // const double mmw_dry_air(28.964e-03); // (28.964e-3) [kg mol-1] (2018 SI NIST RRTMGP)
  const double mmw_dry_air(28.9644e-03); // (28.9644e-3) [kg mol-1] (Source: radcsw.F in CCM2/3)
  const double mmw_illite(389.34e-03); // (389.34e-03) [kg mol-1] Mean molecular weight of Illite http://webmineral.com/data/Illite.shtml
  const double mmw_kaolinite(258.16e-03); // (258.16e-03) [kg mol-1] Mean molecular weight of Kaolinite http://webmineral.com/data/Kaolinite.shtml
  const double mmw_montmorillonite(549.07e-03); // (549.07e-03) [kg mol-1] Mean molecular weight of Montmorillonite http://webmineral.com/data/Montmorillonite.shtml
  const double mss_earth(5.98e+24); // (5.98e+24) [kg] Mass of the Earth
  const prc_cmp msv_snw_std(0.97); // [frc] Emissivity of snow CCM:lsm/snoconi
  const double prs_1000(100000.0); // (100000.0) [Pa] Reference pressure for potential temperature
  const double prs_HITRAN(101325.0); // (101325.0) [Pa] Reference pressure for HITRAN database
  const double prs_STP(101325.0); // (101325.0) [Pa] Standard mean sea-level pressure
  const double r_O2(0.23143); // (0.23143) [kg kg-1] (Source: radcsw.F in CCM2/3)7
  const double rds_earth(6.370e+06); // (6.370e+06) [m] Radius of sphere of same volume as Earth
  const double sfc_tns_wtr_lqd_STP(7.610e-3); // [J m-2] Surface tension of liquid water at STP
  const double slr_cst_CCM(1367.0); // (1367.0) [W m-2] Solar constant used in CCM
  const double slr_cst_FDE24_hst(1361.353); // (1361.353) [W m-2] Solar constant 1850-2023 daily mean
  const double slr_cst_FDE24_PI_old(1361.37); // (1361.37) [W m-2] Solar constant "pre-industrial" 18500101-18730128 daily mean HEPPA 4-3
  const double slr_cst_FDE24_PI(1361.617); // (1361.617) [W m-2] Solar constant "pre-industrial" 18500101-18730128 daily mean HEPPA 4-5, 4-6
  const double slr_cst_FDE24_mnm(1360.8); // (1360.8) [W m-2] Solar constant "quiet" during solar minimum used in CMIP7 (FDE24 p. 1218)
  const double slr_cst_MFA17_hst(1360.86); // (1360.86) [W m-2] Solar constant 1850-2014 monthly mean
  const double slr_cst_MFA17_PI(1360.747); // (1360.747) [W m-2] Solar constant "pre-industrial" 18500101-18730128 daily mean used in CMIP6
  const double spc_heat_H2O_ice_vlm(2.094e06); // (2.094e06) [J m-3 K-1] Volumetric specific heat capacity of ice water CCM:lsm/phyconi.F
  const double spc_heat_H2O_ice(2108.0); // (2108.0) [J kg-1 K-1] Specific heat capacity of ice water (quora.com)
  const double spc_heat_H2O_lqd(4187.0); // (4187.0) [J kg-1 K-1] Specific heat capacity of liquid water RoY94 p. 15 (fxm: better to derive this from density and volumetric spec. heat?)
  const double spc_heat_H2O_lqd_vlm(4.188e06); // (4.188e06) [J m-3 K-1] Volumetric specific heat capacity of liquid water CCM:lsm/phyconi.F
  const double spc_heat_H2SO4_lqd(1416.2); // (1416.2) [J kg-1 K-1] Specific heat capacity of sulfuric acid CRC95 p. 6-145
  const double spc_heat_HNO3_lqd(1744.1); // (1744.1) [J kg-1 K-1] Specific heat capacity of nitric acid CRC95 p. 6-145
  const double spc_heat_NaCl_sld(854.0); // (854.0) [J kg-1 K-1] Specific heat capacity of sodium chloride http://www.crystran.co.uk/nacldata.htm
  const double ssh_H2O_sln_crc(0.98); // (0.98) [frc] Salinity correction to saturation specific humidity of H2O LaP81 p. 328 CCM:dom/flxoce()
  const double tpt_298(298.00); // (298.00) [K] Standard temperature for Henry's Law
  const double tpt_HITRAN(296.0); // (296.0) [K] Reference temperature for HITRAN database
  const double tpt_STP(273.16); // (273.16) [K] Standard temperature
  const double tpt_frz_pnt(273.15); // (273.15) [K] Kelvin--Celsius scale offset, see, e.g., Bol80
  const double tpt_ocn_sea_ice_cls(-1.7999); // (-1.7999) [C] Standard temperature at which sea water freezes to sea ice CCM:dom/parsst.h/tsice
  const double tpt_trp_pnt(273.16); // (273.16) [K] Standard temperature, see, e.g., Bol80
  const double vmr_std_A(0.00934); // (0.00934) [mlc mlc-1] Volume mixing ratio of A in standard dry air GoY89 p. 9 Tbl 1.1
  const double vmr_std_CH4(1.6e-6); // (1.6e-6) [mlc mlc-1] Volume mixing ratio of CH4 in standard dry airGoY89 p. 9 Tbl 1.1 
  const double vmr_std_CO2(355.0e-6); // (355.0e-6) [mlc mlc-1] Volume mixing ratio of CO2 in standard dry air GoY89 p. 9 Tbl 1.1 (updated to present)
  const double vmr_std_N2(0.78084); // (0.78084) [mlc mlc-1] Volume mixing ratio of N2 in standard dry air GoY89 p. 9 Tbl 1.1
  const double vmr_std_O2(0.20946); // (0.20946) [mlc mlc-1] Volume mixing ratio of O2 in standard dry air GoY89 p. 9 Tbl 1.1
  const double vsc_dyn_H2O(0.001); // (0.001) [kg m-1 s-1] Dynamic viscosity of liquid H2O PrK98 p. 387
  const prc_cmp cnd_trm_H2O_ice(2.2); // (2.2) [W m-1 K-1] Thermal conductivity of ice water CCM:lsm/phyconi
  const prc_cmp cnd_trm_H2O_lqd(0.6); // (0.6) [W m-1 K-1] Thermal conductivity of liquid water CCM:lsm/phyconi
  const prc_cmp cnd_trm_snw(0.34); // (0.34) [W m-1 K-1] Thermal conductivity of snow CCM:lsm/snoconi fxm: discrepancy between CCM and CRC95
  //  const prc_cmp cnd_trm_snw(0.16); // (0.16) [W m-1 K-1] Thermal conductivity of snow at 0 C CRC95 p. 12-180 fxm: discrepancy between CCM and CRC95
  const prc_cmp cnd_trm_H2O_ocn(0.5984); // (0.5984) [W m-1 K-1] Thermal conductivity of liquid water at 20 C CRC95 p. 6-11 NB: varies by 10% from T=0--30 C fxm
  const prc_cmp cnd_trm_NaCl_273K(1.15); // (1.15) [W m-1 K-1] Thermal conductivity of sodium chloride at 273 K http://www.crystran.co.uk/nacldata.htm
  const prc_cmp cnd_trm_dry_air_300K(0.0262); // (0.0262) [W m-1 K-1] Thermal conductivity of air at 300 K CRC95 p. 6-251, RoY94 p. 103
  const prc_cmp cnd_trm_snd(0.33); // (0.33) [W m-1 K-1] Thermal conductivity of dry sand at 20 C CRC95 p. 12-180
  const prc_cmp dns_CaCO3(2710.0); // (2710.0) [kg m-3] Density of calcite (calcium carbonate) (http://webmineral.com/data/Calcite.shtml)
  const prc_cmp dns_CaSO4_2H2O(2300.0); // (2300.0) [kg m-3] Density of gypsum (hydrous calcium sulfate) (http://webmineral.com/data/Gypsum.shtml)
  const prc_cmp dns_Fe2O3(5260.0); // (5260.0) [kg m-3] Density of hematite crystals (http://webmineral.com/data/Hematite.shtml)
  const prc_cmp dns_FeOOH(3800.0); // (3800.0) [kg m-3] Density of goethite crystals (http://webmineral.com/data/Goethite.shtml)
  const prc_cmp dns_SiO2(2620.0); // (2620.0) [kg m-3] Density of SiO2 (http://webmineral.com/data/Quartz.shtml)
  const prc_cmp dns_H2O_ice_MYM02(910.0); // (910.0) [kg m-3] Density of ice crystals (MYM02 p. 7)
  const prc_cmp dns_H2O_ice_NGW03(917.0); // (917.0) [kg m-3] Density of ice crystals (NGW03 p. 6)
  const prc_cmp dns_H2O_ice_std(917.0); // (917.0) [kg m-3] Density of ice crystals (NGW03 p. 6)
  const prc_cmp dns_H2O_lqd_std(1000.0); // (1000.0) [kg m-3] Density of liquid water
  const prc_cmp dns_H2O_snw_gnd_std(100.0); // (100.0) [kg m-3] Standard bulk density of snow on ground WiW80 p. 2724, 2725, CCM:physics/tsinti()
  const prc_cmp dns_H2O_snw_gnd_LSM(250.0); // (250.0) [kg m-3] Bulk density of packed snow on ground CCM:lsm/snoconi.F
  const prc_cmp dns_H2SO4_61pct_273K_std(1526.16); // (1526.16) [kg m-3] Density of 61% H2SO4 solution (by weight) at 273 K (Timmermans in $HOME/sw/idx_rfr/hitran/timmerma.dat)
  const prc_cmp dns_H2SO4_72pct_273K_std(1652.95); // (1652.95) [kg m-3] Density of 72% H2SO4 solution (by weight) at 273 K (Timmermans in $HOME/sw/idx_rfr/hitran/timmerma.dat)
  const prc_cmp dns_H2SO4_75pct_273K_std(1688.83); // (1688.83) [kg m-3] Density of 75% H2SO4 solution (by weight) at 273 K (Timmermans in $HOME/sw/idx_rfr/hitran/timmerma.dat)
  const prc_cmp dns_H2SO4_75pct_303K_std(1659.65); // (1659.65) [kg m-3] Density of 75% H2SO4 solution (by weight) at 303 K (Timmermans in $HOME/sw/idx_rfr/hitran/timmerma.dat)
  const prc_cmp dns_H2SO4_std(1851.69); // (1851.69) [kg m-3] Density of 100% H2SO4 at 273 K (Timmermans in $HOME/sw/idx_rfr/hitran/timmerma.dat)
  const prc_cmp dns_MgSO4_std(2660.0); // (2660.0) [kg m-3] Density of MgSO4 Yu's Thesis Ch. 4
  const prc_cmp dns_NH4HSO4_std(1780.0); // (1780.0) [kg m-3] Density of NH4HSO4 Yu PhD Tbl. 4.2
  const prc_cmp dns_NH4NH4SO4_std(1769.0); // (1769.0) [kg m-3] Density of NH4NH4SO4 Yu PhD Tbl. 4.2
  const prc_cmp dns_NH4NO3_std(1725.0); // (1725.0) [kg m-3] Density of NH4NO3 Yu PhD Tbl. 4.2
  const prc_cmp dns_NaCl_std(2170.0); // (2170.0) [kg m-3] Density of NaCl http://www.crystran.co.uk/nacldata.htm
  const prc_cmp dns_dst_MCK09(2060.0); // (2060.0) [kg m-3] Density of dust MCK09 p. 8478 (dry bulk Chinese Loess from Liu, Loess in China, 1985)
  const prc_cmp dns_dst_RJM03(2000.0); // (2000.0) [kg m-3] Density of dust RJM03 ("dust effective density")
  const prc_cmp dns_dst_Cra97(2560.0); // (2560.0) [kg m-3] Density of crustal material Cra97 (Soil Mechanics book)
  const prc_cmp dns_dst_DKS91(1600.0); // (1600.0) [kg m-3] Density of dust DKS91 p. 118
  const prc_cmp dns_dst_MSI03(2000.0); // (2000.0) [kg m-3] Density of dust MSI03 p. 5
  const prc_cmp dns_dst_PaG77(2500.0); // (2500.0) [kg m-3] Density of dust PaG77 p. 2076
  const prc_cmp dns_dst_Vol73(2650.0); // (2650.0) [kg m-3] Density of dust Vol73
  const prc_cmp dns_dst_mtr_DKS91(2500.0); // (2500.0) [kg m-3] Density of meteoric dust DKS91 p. 118
  const prc_cmp dns_dst_mtr_HTT80(2000.0); // (2000.0) [kg m-3] Density of meteoric dust HTT80 p. 1348
  const prc_cmp dns_dst_std(2500.0); // (2500.0) [kg m-3] Standard density of dust
  const prc_cmp dns_mtr_dst_std(2000.0); // (2000.0) [kg m-3] Standard density of meteoric dust HTT80 p. 1348
  const prc_cmp dns_illite(2750.0); // (2750.0) [kg m-3] Density of illite http://webmineral.com/data/Illite.shtml
  const prc_cmp dns_kaolinite(2600.0); // (2600.0) [kg m-3] Density of kaolinite http://webmineral.com/data/Kaolinite.shtml
  const prc_cmp dns_montmorillonite(2350.0); // (2350.0) [kg m-3] Density of montmorillonite http://webmineral.com/data/Montmorillonite.shtml
  const prc_cmp dns_mdp_max(2000.0); // (2000.0) [kg m-3] Minimum valid density of particle
  const prc_cmp dns_prt_min(0.1); // (0.1) [kg m-3] Minimum valid density of particle
  const prc_cmp dns_lac_DKS91(2300.0); // (2300.0) [kg m-3] Density of soot DKS91 p. 118
  const prc_cmp dns_lac_CLV96(2000.0); // (2000.0) [kg m-3] Density of soot CLV96 p. 23366
  const prc_cmp dns_lac_JaC04(1500.0); // (1500.0) [kg m-3] Density of soot Jac04, based on Maricq et al. (2000), Fgr. 5
  const prc_cmp dns_lac_HKS98(1000.0); // (1000.0) [kg m-3] Density of soot HKS98 p. 836
  const prc_cmp dns_wsoc_HKS98(1800.0); // (1800.0) [kg m-3] Density of water soluble organic/inorganic mixture HKS98 p. 836
  const prc_cmp dns_lac_WaW80(2050.0); // (2050.0) [kg m-3] Density of soot WaW80 p. 2739
  const prc_cmp dns_lac_YZS00(1860.0); // (1860.0) [kg m-3] Density of soot YZS01 p. 3971
  const prc_cmp dns_lac_BoB05(1800.0); // (1800.0) [kg m-3] Density of soot BoB05 p. 52
  const prc_cmp dns_lac_BoB052(2000.0); // (2000.0) [kg m-3] Density of soot BoB05 p. 50 Tbl. 6 
  const prc_cmp dns_lac_ZGD09(1322.0); // (1322.0) [kg m-3] Density of soot ZGD09 (tuned so ChC90 MAC = 7500 m2 kg-1 at 550 nm)
  const prc_cmp dns_lac_std(1800.0); // (1800.0) [kg m-3] Standard density of soot

  const prc_cmp rhd_NH4HSO4_std(0.40); // (0.40) [frc] Relative humidity at deliquescence of NH4HSO4 Yu PhD Tbl. 4.2
  const prc_cmp rhd_NH4NH4SO4_std(0.80); // (0.80) [frc] Relative humidity at deliquescence of NH4NH4SO4 Yu PhD Tbl. 4.2
  const prc_cmp rhd_NH4NO3_std(0.62); // (0.62) [frc] Relative humidity at deliquescence of NH4NO3 Yu PhD Tbl. 4.2
  const prc_cmp rhd_NaCl_std(0.753); // (0.753) [frc] Relative humidity at deliquescence of NaCl Yu PhD Tbl. 4.2
  const std::complex<prc_cmp> idx_rfr_MgSO4_std(1.560,0.0); // (1.560,0.0) [frc] Standard refractive index of MgSO4 Yu PhD Tbl. 4.2
  const std::complex<prc_cmp> idx_rfr_NH4HSO4_std(1.473,0.0); // (1.473,0.0) [frc] Standard refractive index of NH4HSO4 Yu PhD Tbl. 4.2
  const std::complex<prc_cmp> idx_rfr_NH4NH4SO4_std(1.521,0.0); // (1.521,0.0) [frc] Standard refractive index of NH4NH4SO4 Yu PhD Tbl. 4.2
  const std::complex<prc_cmp> idx_rfr_NH4NO3_std(0.0,0.0); // (0.0,0.0) [frc] Standard refractive index of NH4NO3 Yu PhD Tbl. 4.2 fxm
  const std::complex<prc_cmp> idx_rfr_NaCl_std(1.544,0.0); // (1.544,0.0) [frc] Standard refractive index of NaCl Yu PhD Tbl. 4.2

  // Derived constants
  const double calories_per_joule(1.0/joules_per_calorie); // (0.23892) [cal J-1] Calorie = energy to heat 1g H20 1C @ 15C
  const double dmt_earth(2.0*rds_earth); // (1.27e+07) [m] Diameter of sphere of same volume as Earth
  const double N2_per_O2(vmr_std_N2/vmr_std_O2); // (3.72787) [mlc mlc-1] Tropospheric N2 to O2 vmr ratio 
  const double N_STP(prs_STP/(cst_Boltzmann*tpt_STP)); // (2.68666e25) [mlc m-3] Air molecule concentration at STP (used for index of refraction pzn)
  const double cst_von_krm_rcp(1.0/cst_von_krm); // (2.5) [frc] Reciprocal of Von Karman's constant
  const double eps_H2O(mmw_H2O/mmw_dry_air); // (0.622) [frc] molec wgt vapor/molec wgt dry air
  const double gas_cst_CO2(gas_cst_unv/mmw_CO2); // (188.92) [J kg-1 K-1] Gas constant of CO2
  const double gas_cst_H2O(gas_cst_unv/mmw_H2O); // (461.65) [J kg-1 K-1] Gas constant of H2O
  const double gas_cst_H2OH2O(gas_cst_unv/mmw_H2OH2O); // (230.76) [J kg-1 K-1] Gas constant of H2OH2O
  const double gas_cst_HNO3(gas_cst_unv/mmw_HNO3); // (131.984) [J kg-1 K-1] Gas constant of HNO3
  const double gas_cst_N2(gas_cst_unv/mmw_N2); // (296.88) [J kg-1 K-1] Gas constant of N2 
  const double gas_cst_NO2(gas_cst_unv/mmw_NO2); // (180.7) [J kg-1 K-1] Gas constant of NO2 (a guess)
  const double gas_cst_O2(gas_cst_unv/mmw_O2); // (259.83) [J kg-1 K-1] Gas constant of O2 
  const double gas_cst_O2O2(gas_cst_unv/mmw_O2O2); // (129.9) [J kg-1 K-1] Gas constant of O2O2 (a guess)
  const double gas_cst_O3(gas_cst_unv/mmw_O3); // (173.2) [J kg-1 K-1] Gas constant of O3 (a guess, because mmw(O3) is always changing)
  const double gas_cst_OH(gas_cst_unv/mmw_OH); // (488.65) [J kg-1 K-1] Gas constant of OH 
  const double gas_cst_SO2(gas_cst_unv/mmw_SO2); // (129.818) [J kg-1 K-1] Gas constant of SO2
  const double gas_cst_dry_air(gas_cst_unv/mmw_dry_air); // (287.05) [J kg-1 K-1] Gas constant of dry_air IrG81 p. 25, p. 245
  const double mmw_CaCO3(mmw_Ca+mmw_C+3.0*mmw_O); // (100.087e-03) [kg mol-1] Mean molecular weight of CaCO3
  const double mmw_CaSO4(mmw_Ca+mmw_S+4.0*mmw_O); // [kg mol-1] Mean molecular weight of CaSO4
  const double mmw_CaSO4_2H2O(mmw_Ca+mmw_S+4.0*mmw_O+2*mmw_H+2.0*mmw_O); // (172.17e-03) [kg mol-1] Mean molecular weight of CaSO4_2H2O
  const double mmw_Fe2O3(2.0*mmw_Fe+3.0*mmw_O); // (159.692e-03) [kg mol-1] Mean molecular weight of Fe2O3
  const double mmw_FeOOH(mmw_Fe+2.0*mmw_O+mmw_H); // (88.85e-03) [kg mol-1] Mean molecular weight of FeOOH
  const double mmw_H2SO4(2.0*mmw_H+mmw_S+4.0*mmw_O); // [kg mol-1] Mean molecular weight of H2SO4
  const double mmw_HO2(mmw_H+2.0*mmw_O); // (33.0067e-03) [kg mol-1] Mean molecular weight of HO2
  const double mmw_N2O5(2.0*mmw_N+5.0*mmw_O); // (108.01e-03) [kg mol-1] Mean molecular weight of N2O5
  const double mmw_NH4(mmw_N+4.0*mmw_H); // (18.0385e-03) [kg mol-1] Mean molecular weight of NH4
  const double mmw_NH4NH4SO4(2*(mmw_N+4.0*mmw_H)+mmw_S+4.0*mmw_O); // (132.141e-03) [kg mol-1] Mean molecular weight of NH4NH4SO4
  const double mmw_MgSO4(mmw_Mg+mmw_S+4.0*mmw_O); // (120.369e-03) [kg mol-1] Mean molecular weight of MgSO4
  const double mmw_NO3(mmw_N+3.0*mmw_O); // (62.0049e-03) [kg mol-1] Mean molecular weight of NO3
  const double mmw_NaCl(mmw_Na+mmw_Cl); // (58.4425e-03) [kg mol-1] Mean molecular weight of NaCl
  const double mmw_SO4(mmw_S+4.0*mmw_O); // (96.0636e-03) [kg mol-1] Mean molecular weight of SO4
  const double mmw_SiO2(mmw_Si+2.0*mmw_O); // (60.0843e-03) [kg mol-1] Mean molecular weight of SiO2

  const double rxn_ntp_HNO3_H2O_298K(8700.0*gas_cst_unv); // [J mol-1] Reaction enthalpy (heat of dissolution) at 298K inferred from E/R quoted by Sep97 p. 391 Tbl. 6.A.1
  const double rxn_ntp_CO2_H2O_298K(-1000.0*gas_cst_unv); // [J mol-1] Reaction enthalpy (heat of dissolution) at 298K inferred from E/R quoted by Sep97 p. 391 Tbl. 6.A.1
  const double rxn_ntp_H2O_H2O_298K(-6710.0*gas_cst_unv); // [J mol-1] Reaction enthalpy (heat of dissolution) at 298K inferred from E/R quoted by Sep97 p. 391 Tbl. 6.A.1
  
  // Doubly derived constants
  const double eps_H2O_rcp_m1(-1.0+1.0/eps_H2O); // (0.60777) [frc] Constant for virtual temperature
  const double g_rcp_R_dry(grv_sfc_mean/gas_cst_dry_air); // (0.03416356) [m K-1] for computing lapse rates
  const double one_mns_eps_H2O(1.0-eps_H2O); // (0.378) [frc] 1 - eps_H2O Constant for saturation specific humidity
  const double spc_heat_H2O_vpr(4.0*gas_cst_H2O); // (1850.0) [J kg-1 K-1] IrG81 pp. 77, 245
  //const double spc_heat_dry_air(1004.64); // (1004.64) [J kg-1 K-1] RRTMGP
  const double spc_heat_dry_air(7.0*gas_cst_dry_air/2.0); // (1004.697) [J kg-1 K-1] IrG81 p. 25
  const double spc_heat_SiO2_sld(1000.0*44.4/mmw_SiO2); // (738962.0) [J kg-1 K-1] Specific heat capacity of quartz CRC95 p. 5-21
  const double spc_heat_Fe2O3_sld(1000.0*103.9/mmw_Fe2O3); // (650627.0) [J kg-1 K-1] Specific heat capacity of hematite (iron oxide) CRC95 p. 5-15
  const double spc_heat_CaSO4_sld(1000.0*99.7/mmw_CaSO4); // () [J kg-1 K-1] Specific heat capacity of gypsum (hydrous calcium sulfate) CRC95 p. 5-9
  const double spc_heat_CaCO3_sld(1000.0*83.5/mmw_CaCO3); // () [J kg-1 K-1] Specific heat capacity of calcite (calcium carbonate) CRC95 p. 5-24
  const double spc_heat_C_sld(1000.0*8.5/mmw_C); // () [J kg-1 K-1] Specific heat capacity of graphite (carbon) CRC95 p. 5-23
  
  // Trebly derived constants
  const double cp_vpr_rcp_cp_dry_m1(spc_heat_H2O_vpr/spc_heat_dry_air-1.0); // (0.83745) Constant for moist specific heat IrG81 p. 77
  const double kappa_dry_air(gas_cst_dry_air/spc_heat_dry_air); // (0.286 = 2/7) [frc] Constant in potential temperature IrG81 p. 25, Tre922 p. 72 
  
} // end physical constant namespace phc

#endif // PHYS_CST_HH
