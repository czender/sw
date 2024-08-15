! $Id$

! Purpose: Store fundamental physical constants in SI units

! Copyright (C) 1994--present Charlie Zender
! This software is distributed under the terms of the GNU General Public License
! See http://www.gnu.org/copyleft/gpl.html for full license text

! These constants should be the same as found in fnd_cst.com and drv_cst.com
! This file serves as the master copy

! Usage:
! Non-HITRAN programs:
!use phys_cst_mdl ! [mdl] Fundamental and derived physical constants

! HITRAN programs:
!use fnd_cst ! Fundamental physical constants
!use htrn_mdl ! [mdl] HITRAN constants, subroutines
! fxm: module drv_cst actually exports fnd_cst and htrn_mdl?
!use drv_cst ! Derived physical constants

module phys_cst_mdl ! [mdl] Fundamental and derived physical constants

  ! Fundamental constants
  real,parameter::Avagadro=6.022045e+23 ! [mlc mol-1] Avagadro's number
  real,parameter::Boltzmann=1.38063e-23 ! [J K-1] Boltzmann's constant
  real,parameter::Joules_per_calorie=4.1855 ! [J cal-1] Cal = energy to heat 1g H20 1C @ 15C
  real,parameter::Loschmidt=2.6871e+25 ! [mlc m-3] Loschmidt's number (molecules of air per cubic meter at STP) (NB: this is derivable in principle)
  real,parameter::Planck=6.62620e-34  ! [J s] Planck's constant
  real,parameter::dns_H2O_ice_std=917.0 ! (917.0) [kg m-3] Density of ice crystals (NGW03 p. 6)
  real,parameter::gas_cst_unv=8.31441 ! [J mol-1 K-1] Universal gas constant
  real,parameter::grv_cst=6.672e-11 ! [N m2 kg-2] Universal gravitational constant
  real,parameter::k_O2_O2=2.6e-29     ! [m3 mlc-1] Equilibrium rate constant for O2 + O2 <-> O2-O2 (Sha77 p. 436, 527)
  real,parameter::grv_mean_sfc=9.80665 ! [m s-2] Mean gravitational acceleration at Earth's surface
  real,parameter::lumens_per_Watt_555nm=683.002 ! [lm W-1] Lumens per Watt at 555 nm (683 exactly at 540e12 Hz) for photopic response (FCE11 p. 2719 (1))
  real,parameter::lumens_per_Watt_scotopic=1699 ! [lm W-1] Lumens per Watt for scotopic response (FCE11 p. 2719 (2))
  ! NB: Defining mean molecular weight in SI units instead of [g mol-1]
  ! keeps Avagadro's number the same and still gets rid of all those pesky factors of 1000.0
  ! High precision mmw's come from HITRAN isotopic distribution in hitran.com
  real,parameter::mmw_CO2=4.4009743e-02 ! [kg mol-1]
  real,parameter::mmw_CO=2.8010445e-02 ! [kg mol-1]
  real,parameter::mmw_CH4=1.6043060e-02 ! [kg mol-1]
  real,parameter::mmw_CFC11=136.0e-03 ! [kg mol-1] (CCM: physics/comvmr.h)
  real,parameter::mmw_CFC12=120.0e-03 ! [kg mol-1] (CCM: physics/comvmr.h)
  real,parameter::mmw_dry_air=28.9644e-03 ! [kg mol-1] (CCM: physics/radcsw.F)
  real,parameter::mmw_O2=3.1998575e-02 ! [kg mol-1]
  real,parameter::mmw_H2O=1.8015259e-02 ! [kg mol-1]
  real,parameter::mmw_H2OH2O=36.03e-03 ! [kg mol-1]
  real,parameter::mmw_OH=1.7006907e-02 ! [kg mol-1]
  real,parameter::mmw_O3=4.7997832e-02 ! [kg mol-1] 
  real,parameter::mmw_O2O2=64.0e-03    ! [kg mol-1] (this is a guess)
  real,parameter::mmw_N2=2.8006147e-02 ! [kg mol-1]
  real,parameter::mmw_O2N2=60.0e-03    ! [kg mol-1] (this is a guess)
  real,parameter::mmw_N2O=4.4012674e-02 ! [kg mol-1]
  real,parameter::mmw_NO=3.0005630e-02 ! [kg mol-1]
  real,parameter::mmw_SO2=6.4046674e-02 ! [kg mol-1]
  real,parameter::mmw_NO2=4.5992904e-02 ! [kg mol-1]
  real,parameter::mmw_H2O2=3.4005480e-02 ! [kg mol-1]
  real,parameter::mss_earth=5.98e+24  ! [kg] Mass of the Earth
  real,parameter::out_of_bounds=1.0e+36     
  real,parameter::prs_1000=100000.0   ! [Pa] Ref. pres. for potential temperature
  real,parameter::prs_HITRAN=101325.0 ! [Pa] Ref. pres. for Malkmus band model
  real,parameter::prs_STP=101325.0    ! [Pa] Ref. pres. for index of refraction
  !  real,parameter::r_O2=.23143         ! [kg kg-1] r = mass fraction of dry air (Source: radcsw.F in CCM2/3)
  real,parameter::rds_earth_ease_ml=6371.228e+03 ! [m] Radius of Earth, authalic sphere, used in EASE ML grids
  real,parameter::rds_earth=6.370e+06 ! [m] Radius of sphere of same volume as Earth
  real,parameter::slr_cst_CCM=1367.0  ! [W m-2] Solar constant used in CCM
  real,parameter::slr_cst_FDE24=1361.353 ! [W m-2] Solar constant 1850-2023 mean used in CMIP7
  real,parameter::speed_of_light=2.99793e+08 ! [m s-] Speed of light in vacuo
  real,parameter::tpt_frz_pnt=273.15  ! [K] Kelvin--Celsius scale offset Bol80
  real,parameter::tpt_trp_pnt=273.16  ! [K] Standard temperature See, e.g., Bol80
  real,parameter::tpt_HITRAN=296.0    ! [K] Ref. temp. for HITRAN database
  real,parameter::tpt_STP=273.16      ! [K] Ref. temp. for index of refraction
  real,parameter::vmr_std_A=0.00934   ! [mlc mlc-1]  Volume mixing ratio of A in standard dry air GoY89 p. 9 Tbl 1.1
  real,parameter::vmr_std_CH4=1.6e-6  ! [mlc mlc-1]  Volume mixing ratio of CH4 in standard dry air GoY89 p. 9 Tbl 1.1 
  real,parameter::vmr_std_CO2=355.0e-6 ! [mlc mlc-1]  Volume mixing ratio of CO2 in standard dry air GoY89 p. 9 Tbl 1.1 (updated to present)
  real,parameter::vmr_std_N2=0.78084  ! [mlc mlc-1]  Volume mixing ratio of N2 in standard dry air GoY89 p. 9 Tbl 1.1
  real,parameter::vmr_std_O2=0.20946  ! [mlc mlc-1  Volume mixing ratio of O2 in standard dry air GoY89 p. 9 Tbl 1.1
  real,parameter::cst_von_krm=0.4     ! [frc] Von Karman's constant
  
  ! Derived constants
  real,parameter::N2_per_O2=vmr_std_N2/vmr_std_O2 ! 3.72787
  real,parameter::eps_H2O=mmw_H2O/mmw_dry_air ! (0.622) molec wgt vapor/molec wgt dry air
  real,parameter::gas_cst_CO2=gas_cst_unv/mmw_CO2 ! (188.92) [J kg-1 K-1]
  real,parameter::gas_cst_CO=gas_cst_unv/mmw_CO ! () [J kg-1 K-1]
  real,parameter::gas_cst_CH4=gas_cst_unv/mmw_CH4 ! () [J kg-1 K-1]
  real,parameter::gas_cst_N2O=gas_cst_unv/mmw_N2O ! () [J kg-1 K-1]
  real,parameter::gas_cst_CFC11=gas_cst_unv/mmw_CFC11 ! () [J kg-1 K-1]
  real,parameter::gas_cst_CFC12=gas_cst_unv/mmw_CFC12 ! () [J kg-1 K-1]
  real,parameter::gas_cst_H2O=gas_cst_unv/mmw_H2O ! (461.65) [J kg-1 K-1]
  real,parameter::gas_cst_H2OH2O=gas_cst_unv/mmw_H2OH2O ! (230.76) [J kg-1 K-1]
  real,parameter::gas_cst_OH=gas_cst_unv/mmw_OH ! (488.65) [J kg-1 K-1] 
  real,parameter::gas_cst_O2=gas_cst_unv/mmw_O2 ! (259.83) [J kg-1 K-1] 
  real,parameter::gas_cst_N2=gas_cst_unv/mmw_N2 ! (296.88) [J kg-1 K-1] 
  real,parameter::gas_cst_O3=gas_cst_unv/mmw_O3 ! (173.2) [J kg-1 K-1] (a guess, because mmw(O3) is always changing)
  real,parameter::gas_cst_O2O2=gas_cst_unv/mmw_O2O2 ! (129.9) [J kg-1 K-1] (a guess)
  real,parameter::gas_cst_NO2=gas_cst_unv/mmw_NO2 ! (180.7) [J kg-1 K-1] (a guess)
  real,parameter::gas_cst_dry_air=gas_cst_unv/mmw_dry_air ! (287.05) [J kg-1 K-1] IrG81 p. 245
  real,parameter::g_rcp_R_dry=grv_mean_sfc/gas_cst_dry_air ! (0.03416356) [m K-1]
  real,parameter::mlr_vlm_STP=gas_cst_unv*tpt_STP/prs_STP ! (22.414) [m3 mol-1]
  real,parameter::N_STP=prs_STP/(Boltzmann*tpt_STP) ! (2.68666e25) [mlc m-3] for computing STP index of refr.
  real,parameter::one_mns_eps_H2O=1.0-eps_H2O ! (0.378) 1 - eps_H2O for computing sat. spec. humidity 
  real,parameter::eps_H2O_rcp_m1=-1.0+1.0/eps_H2O ! (0.6077) [frc] for computing virtual Temperature
  real,parameter::spc_heat_H2O_vpr=4.0*gas_cst_H2O ! (1850.0) [J kg-1 K-1] IrG81 pp. 77, 245
  real,parameter::spc_heat_dry_air=7.0*gas_cst_dry_air/2.0 ! (1004.675) [J kg-1 K-1] IrG81 p. 25
  real,parameter::cst_von_krm_rcp=1.0/cst_von_krm ! (2.5) [frc] Reciprocal of Von Karman's constant
  
  ! Trebly derived constants
  real,parameter::kappa_dry_air=gas_cst_dry_air/spc_heat_dry_air ! (0.286 = 2/7) [frc] Constant in potential temperature IrG81 p. 25, Tre922 p. 72 
  real,parameter::cp_vpr_rcp_cp_dry_m1=spc_heat_H2O_vpr/spc_heat_dry_air-1.0 ! 0.87 [frc] for computing moist specific heatIrG81 pp. 77
  
end module phys_cst_mdl
