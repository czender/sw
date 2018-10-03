! $Id$ 

! Purpose: Thermodynamic utilities for Fortran programs

! Copyright (C) 1994--2018 Charlie Zender
! This software is distributed under the terms of the GNU General Public License
! See http://www.gnu.org/copyleft/gpl.html for full license text

! Usage:
! use tdy_mdl ! [mdl] Thermodynamics

module tdy_mdl ! [mdl] Thermodynamics
  implicit none
  private ! [stt] Symbols are private unless individually qualified as public
  public::svp_H2O_ice_GoG ! [fnc] Saturation vapor pressure over planar ice water, GoG
  public::svp_H2O_ice_PrK78 ! [fnc] Saturation vapor pressure over planar ice water, PrK78
  public::svp_H2O_ice_PrK98 ! [fnc] Saturation vapor pressure over planar ice water, PrK98
  public::svp_H2O_lqd_Bol80 ! [fnc] Saturation vapor pressure over planar liquid water, Bol80
  public::svp_H2O_lqd_PrK78 ! [fnc] Saturation vapor pressure over planar liquid water, PrK78
  public::svp_H2O_lqd_Wex ! [fnc] Saturation vapor pressure over planar liquid water, Wex
  public::tpt_vrt_get ! [fnc] Virtual temperature
  
contains
  
  real function svp_H2O_lqd_Wex( & ! [fnc] Saturation vapor pressure over planar liquid water
       tpt) ! [K] Temperature
    ! Purpose: Compute saturation vapor pressure over planar liquid water
    ! Input temperature in degrees kelvin 
    ! Saturation vapor pressure returned in [Pa]
    ! See Bol80 for details
    ! fxm: This function currently does not work, for unknown reasons
    implicit none
    real,intent(in)::tpt ! [K]
    integer idx
    real(selected_real_kind(p=12))::svp_H2O_lqd_dbl
    real(selected_real_kind(p=12))::tpt_dbl
    real(selected_real_kind(p=12))::cff(0:7)= &
         (/-2.9912729e3,-6.0170128e3,1.887643854e1,-2.8354721e-2,1.7838301e-5,-8.4150417e-10,4.4412543e-13,2.858487e0/)
    ! Main code
    tpt_dbl=dble(tpt)
    svp_H2O_lqd_dbl=0.0
    do idx=0,6
       svp_H2O_lqd_dbl=svp_H2O_lqd_dbl+cff(idx)*tpt_dbl**(idx-2)
    enddo
    svp_H2O_lqd_dbl=svp_H2O_lqd_dbl+cff(7)*log(tpt_dbl)
    svp_H2O_lqd_dbl=100.0*exp(svp_H2O_lqd_dbl) ! [mb] --> [Pa]
    svp_H2O_lqd_Wex=real(svp_H2O_lqd_dbl)
    ! Sanity check
    if (svp_H2O_lqd_Wex < 0.0) stop 'svp_H2O_lqd_Wex < 0.0 in svp_H2O_lqd_Wex()'
    return
  end function svp_H2O_lqd_Wex                       ! end svp_H2O_lqd_Wex()
  
  real function svp_H2O_ice_GoG( & ! [fnc] Saturation vapor pressure over planar ice water, GoG
       tpt) ! [K] Temperature
    ! Purpose:
    ! Compute saturation vapor pressure over planar ice water
    ! Input temperature in degrees kelvin 
    ! Saturation vapor pressure returned in [Pa]
    ! Stolen from CCM2 gffgch.F by J.J. Hack and G.T. Taylor
    ! Source reference unknown
    implicit none
    real,intent(in)::tpt                  ! [K]
    real tpt0
    real term1
    real term2
    real term3
    ! Main code
    tpt0=273.16
    term1=2.01889049/(tpt0/tpt)
    term2=3.56654*log(tpt0/tpt)
    term3=20.947031*(tpt0/tpt)
    svp_H2O_ice_GoG=575.185606e10*exp(-(term1+term2+term3))
    return
  end function svp_H2O_ice_GoG                       ! end svp_H2O_ice_GoG()
  
  real function svp_H2O_ice_PrK78( & ! [fnc] Saturation vapor pressure over planar ice water, PrK78
       tpt) ! [K] Temperature
    ! Purpose:
    ! Compute saturation vapor pressure over planar ice water
    ! Input temperature in degrees kelvin 
    ! Saturation vapor pressure returned in [Pa]
    ! Taken from Lowe and Ficke (1974) as reported in PrK78 p. 625
    ! Range of validity is -50 C < T < 0 C
    use phys_cst_mdl,only:tpt_frz_pnt ! [mdl] Fundamental and derived physical constants
    implicit none
    real,intent(in)::tpt ! [K]
    real(selected_real_kind(p=12))::tpt_cls ! [C]
    real(selected_real_kind(p=12))::svp_H2O_ice_dbl
    real(selected_real_kind(p=12))::cff(0:6)= &
         (/6.109177956,5.034698970e-1,1.886013408e-2,4.176223716e-4,5.824720280e-6,4.838803174e-8,1.838826904e-10/)
    ! Main code
    tpt_cls=tpt-tpt_frz_pnt ! [C]
    svp_H2O_ice_dbl= & ! [mb]
         cff(0)+tpt_cls*(cff(1)+tpt_cls*(cff(2)+tpt_cls*(cff(3)+tpt_cls*(cff(4)+tpt_cls*(cff(5)+cff(6)*tpt_cls)))))
    svp_H2O_ice_PrK78=real(svp_H2O_ice_dbl)*100.0 ! [mb] --> [Pa]
    return
  end function svp_H2O_ice_PrK78                       ! end svp_H2O_ice_PrK78()
  
  real function svp_H2O_ice_PrK98( & ! [fnc] Saturation vapor pressure over planar ice water, PrK98
       tpt) ! [K] Temperature
    ! Purpose:
    ! Compute saturation vapor pressure over planar ice water
    ! Input temperature in degrees kelvin 
    ! Saturation vapor pressure returned in [Pa]
    ! Taken from Lowe and Ficke (1974, 1981) as reported in PrK98 p. 854
    ! Range of validity is -50 C < T < 0 C
    use phys_cst_mdl,only:tpt_frz_pnt ! [mdl] Fundamental and derived physical constants
    implicit none
    real,intent(in)::tpt                  ! [K]
    real(selected_real_kind(p=12))::tpt_cls  ! [C]
    real(selected_real_kind(p=12))::svp_H2O_ice_dbl
    real(selected_real_kind(p=12))::cff(0:6)= &
         (/6.10690449,5.02660639e-1,1.87743264e-2,4.13476180e-4,5.72333773e-6,4.71651246e-8,1.78086695e-10/)
    ! Main code
    tpt_cls=tpt-tpt_frz_pnt ! [C]
    svp_H2O_ice_dbl= & ! [mb]
         cff(0)+tpt_cls*(cff(1)+tpt_cls*(cff(2)+tpt_cls*(cff(3)+tpt_cls*(cff(4)+tpt_cls*(cff(5)+cff(6)*tpt_cls)))))
    svp_H2O_ice_PrK98=real(svp_H2O_ice_dbl)*100.0 ! [mb] --> [Pa]
    return
  end function svp_H2O_ice_PrK98                       ! end svp_H2O_ice_PrK98()
  
  real function svp_H2O_lqd_PrK78( & ! [fnc] Saturation vapor pressure over planar liquid water, PrK78
       tpt) ! [K] Temperature
    ! Purpose:
    ! Compute saturation vapor pressure over planar liquid water
    ! Input temperature in degrees kelvin 
    ! Saturation vapor pressure returned in [Pa]
    ! Taken from PrK78 p. 625
    ! Range of validity is -50 C < T < 50 C
    use phys_cst_mdl,only:tpt_frz_pnt ! [mdl] Fundamental and derived physical constants
    implicit none
    real,intent(in)::tpt                  ! [K]
    real(selected_real_kind(p=12))::tpt_cls  ! [C]
    real(selected_real_kind(p=12))::svp_H2O_lqd_dbl
    real(selected_real_kind(p=12))::cff(0:6)= &
         (/6.107799961,4.436518521e-1,1.428945805e-2,2.650648471e-4,3.031240396e-6,2.034080948e-8,6.136820929e-11/)
    ! Main code
    tpt_cls=tpt-tpt_frz_pnt   ! [C]
    if (tpt_cls > -50.0) then
       svp_H2O_lqd_dbl= & ! [mb]
            cff(0)+tpt_cls*(cff(1)+tpt_cls*(cff(2)+tpt_cls*(cff(3)+tpt_cls*(cff(4)+tpt_cls*(cff(5)+cff(6)*tpt_cls)))))
       svp_H2O_lqd_PrK78=real(svp_H2O_lqd_dbl)*100.0 ! [mb] --> [Pa]
    else                      ! 
       svp_H2O_lqd_PrK78=svp_H2O_lqd_Bol80(tpt)
    endif                     ! endif
    return
  end function svp_H2O_lqd_PrK78                       ! end svp_H2O_lqd_PrK78()
  
  real function svp_H2O_lqd_Bol80( & ! [fnc] Saturation vapor pressure over planar liquid water, Bol80
       tpt) ! [K] Temperature
    ! Purpose:
    ! Compute saturation vapor pressure over planar liquid water
    ! Input temperature in degrees kelvin 
    ! Saturation vapor pressure returned in [Pa]
    ! Taken from Bol80
    ! This Bol80 parameterization has an accuracy of 0.1% for -30 < tpt_cls < 35 C
    use phys_cst_mdl,only:tpt_frz_pnt ! [mdl] Fundamental and derived physical constants
    implicit none
    real,intent(in)::tpt ! [K]
    real(selected_real_kind(p=12))::tpt_cls ! [C]
    real(selected_real_kind(p=12))::svp_H2O_lqd_dbl
    ! Main code
    tpt_cls=tpt-tpt_frz_pnt ! [C]
    svp_H2O_lqd_dbl=6.112*exp(17.67*tpt_cls/(tpt_cls+243.5)) ! [mb] Bol80
    svp_H2O_lqd_Bol80=real(svp_H2O_lqd_dbl)*100.0 ! [mb] --> [Pa]
    return
  end function svp_H2O_lqd_Bol80                       ! end svp_H2O_lqd_Bol80()
  
  real function tpt_vrt_get( & ! [fnc] Virtual temperature
       tpt, & ! [K] Temperature
       q_H2O_vpr) ! [kg kg-1] Water vapor mixing ratio
    ! Purpose: Return virtual temperature
    ! Prototype:
    ! float tpt_vrt_get       ! [K] Virtual temperature (libcsz_f90)
    use phys_cst_mdl,only:eps_H2O_rcp_m1 ! [mdl] Fundamental and derived physical constants
    implicit none
    real,intent(in)::tpt                  ! [K] Temperature
    real,intent(in)::q_H2O_vpr            ! [kg kg-1] Water vapor mixing ratio
    ! Main Code
    tpt_vrt_get=tpt*(1.0+eps_H2O_rcp_m1*q_H2O_vpr) ! [K]
    return
  end function tpt_vrt_get                       ! end tpt_vrt_get()
  
end module tdy_mdl ! [mdl] Thermodynamics
