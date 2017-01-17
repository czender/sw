! $Id$

! Purpose: Store derived physical constants in SI units

! Copyright (C) 1994--2017 Charlie Zender
! This software is distributed under the terms of the GNU General Public License
! See http://www.gnu.org/copyleft/gpl.html for full license text

! Usage:
!use drv_cst_mdl ! [mdl] Derived physical constants

module drv_cst_mdl ! [mdl] Derived physical constants
  
  use fnd_cst_mdl ! [mdl] Fundamental physical constants
  use htrn_mdl ! [mdl] HITRAN constants, subroutines

  ! Derived constants
  real,parameter::N2_per_O2=vmr_std_N2/vmr_std_O2 ! 3.72787
  real,parameter::eps_H2O=mmw_H2O/mmw_dry_air ! (0.622) molec wgt vapor/molec wgt dry air
  real,parameter::gas_cst_CO2=gas_cst_unv/mmw_CO2 ! (188.92) [J kg-1 K-1]
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
  
end module drv_cst_mdl
