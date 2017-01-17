! $Id$

! Purpose: Physical constants and anthropogenic parameters/assumptions related to skyglow

! Copyright (C) 2016--2017 Charlie Zender
! This software is distributed under the terms of the GNU General Public License
! See http://www.gnu.org/copyleft/gpl.html for full license text

! Usage:
!use skg_mdl ! [mdl] Skyglow parameters

module skg_mdl ! [mdl] Physical constants and anthropogenic parameters/assumptions related to skyglow

  ! Pupil factor parameters
  real,parameter::ppl_age_yr_std=23.0 ! [yr] Mean age of participants in threshold brightness tests

  ! Garstang model parameters Gar00 p. 84 (3), CFE01 p. 37 (19)
  real,parameter::cst_one=3.451e-9 ! 
  real,parameter::cst_two=4.276e-8 ! 
  real,parameter::kst_one=0.109 ! 
  real,parameter::kst_two=1.51e-3 ! 
  real,parameter::yyy_one=2.0e-5 ! 
  real,parameter::yyy_two=1.29e-3 ! 
  real,parameter::zzz_one=0.174 ! 
  real,parameter::zzz_two=0.0587 ! 
  real,parameter::alpha_one=2.35e-4 ! 
  real,parameter::alpha_two=5.81e-3 ! 

end module skg_mdl
