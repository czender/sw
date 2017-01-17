! $Id$

! Purpose: Description, parameters, and nomenclature for various wavelength grids

! Copyright (C) 1994--2017 Charlie Zender
! This software is distributed under the terms of the GNU General Public License
! See http://www.gnu.org/copyleft/gpl.html for full license text

! Usage:
! use wvl_mdl ! [mdl] Wavelength grid parameters

module wvl_mdl

  ! Wavelength grid types are as follows:
  integer,parameter::wvl_grd_dfl=0 ! Default grid
  integer,parameter::wvl_grd_CCM_SW=1 ! CCM SW grid
  integer,parameter::wvl_grd_CCM_LW=2 ! CCM LW grid
  integer,parameter::wvl_grd_NBM_SW=3 ! NBM SW grid
  integer,parameter::wvl_grd_NBM_LW=4 ! NBM LW grid

  ! CCM grid parameters for non-CCM programs are used in 
  ! wvl_grd_CCM_SW_mk(), wvl_grd_CCM_LW_mk(), slr_spc_get_CCM(), and rbn_vec_CCM()
  integer,parameter::bnd_nbr_CCM_max=19 ! [nbr] Number of CCM bands
  integer,parameter::bnd_nbr_CCM_SW_max=19 ! [nbr] Number of CCM3 SW bands
  integer,parameter::bnd_nbr_CCM_LW_max=6 ! [nbr] Number of CCM3 LW bands

end module wvl_mdl


