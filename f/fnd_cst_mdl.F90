! $Id$

! Purpose: Store fundamental physical constants in SI units

! Copyright (C) 1994--2018 Charlie Zender
! This software is distributed under the terms of the GNU General Public License
! See http://www.gnu.org/copyleft/gpl.html for full license text

! Usage:
!use fnd_cst_mdl ! [mdl] Fundamental physical constants

module fnd_cst_mdl ! [mdl] Fundamental physical constants
  
  ! Fundamental constants
  real,parameter::Avagadro=6.022045e+23 ! [mlc mol-1] Avagadro's number
  real,parameter::Boltzmann=1.38063e-23 ! [J K-1] Boltzmann's constant
  real,parameter::Joules_per_calorie=4.1855 ! [J cal-1] Cal = energy to heat 1g H20 1C @ 15C
  real,parameter::Loschmidt=2.6871e+25 ! [mlc m-3] Loschmidt's number (molecules of air per cubic meter at STP) (NB: this is derivable in principle)
  real,parameter::Planck=6.62620e-34  ! [J s] Planck's constant
  real,parameter::gas_cst_unv=8.31441 ! [J mol-1 K-1] Universal gas constant
  real,parameter::grv_cst=6.672e-11 ! [N m2 kg-2] Universal gravitational constant
  real,parameter::k_O2_O2=2.6e-29     ! [m3 mlc-1] Equilibrium rate constant for O2 + O2 <-> O2-O2 (Sha77 p. 436, 527)
  real,parameter::grv_mean_sfc=9.80665 ! [m s-2] Mean gravitational acceleration at Earth's surface
  ! NB: Defining mean molecular weight in SI units instead of [g mol-1]
  ! keeps Avagadro's number the same and still gets rid of all those pesky factors of 1000.0
  ! High precision mmw's come from HITRAN isotopic distribution in hitran.com
  real,parameter::mmw_CFC11=136.0e-03 ! [kg mol-1] (CCM: physics/comvmr.h)
  real,parameter::mmw_CFC12=120.0e-03 ! [kg mol-1] (CCM: physics/comvmr.h)
  real,parameter::mmw_dry_air=28.9644e-03 ! [kg mol-1] (CCM: physics/radcsw.F)
  real,parameter::mmw_H2OH2O=36.03e-03 ! [kg mol-1]
  real,parameter::mmw_O2O2=64.0e-03    ! [kg mol-1] (this is a guess)
  real,parameter::mmw_O2N2=60.0e-03    ! [kg mol-1] (this is a guess)
  real,parameter::mss_earth=5.98e+24  ! [kg] Mass of the Earth
  real,parameter::out_of_bounds=1.0e+36     
  real,parameter::prs_1000=100000.0   ! [Pa] Ref. pres. for potential temperature
  real,parameter::prs_HITRAN=101325.0 ! [Pa] Ref. pres. for Malkmus band model
  real,parameter::prs_STP=101325.0    ! [Pa] Ref. pres. for index of refraction
  !  real,parameter::r_O2=.23143         ! [kg kg-1] r = mass fraction of dry air (Source: radcsw.F in CCM2/3)
  real,parameter::rds_earth=6.370e+06 ! [m] Radius of sphere of same volume as Earth
  real,parameter::slr_cst_CCM=1367.0  ! [W m-2] Solar constant used in CCM
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
  
end module fnd_cst_mdl
