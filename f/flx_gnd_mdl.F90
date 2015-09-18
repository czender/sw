! The routines in this module were gathered from
! the CCSM 3.0
! ../models/lnd/clm/src/biogeophys/baregroundfluxes.F90

module flx_gnd_mdl
  use abl_typ_mdl
  implicit none
    
contains
  subroutine flx_gnd_get(pfl_typ,sfc_typ,lev_nbr)
    implicit none
    ! Calculate surface sensible and latent heat fluxes [W m-2]
    ! Input
    integer, intent(in)::lev_nbr
    type(pfl_mdp), intent(inout), dimension(lev_nbr) :: pfl_typ
    type(sfc_var), intent(inout)                     :: sfc_typ
    
    ! locals
    real dqdz      ! [kg kg-1 m-1] Spec hmd gradient
    real qstar     ! []      Spec hmd scale
    real theta_star! []      Potential temp scale 
    real dtheta    ! [K m-1] Potential temp difference
    real theta_avg ! [K]     Avg theta in the lowest layer
    real thetav    ! [K]     Virtual Potential Temp
    real dthetav   ! [K m-1] Thetav diff between sfc and lowest layer
    real thvstar   ! []      Virt pot temp scale
    real obu_lng   ! [m]     Monin-Obukhov length
    real um        ! [m s-1] wind speed including the stability effect
    real rgh_mmn   ! [m]     Momentum Roughness length
    real zldis     ! [m]     Ref hght "minus" zero displacement height
    real displa    ! [m]     Displacement height
    integer iter   ! [idx]   Iteration index
    real rgh_heat
    real rgh_vpr
    real temp1
    real temp2
    real temp12m
    real temp22m
    real ustar     ! [m s-1] Friction velocity
    real wstar     ! [m s-1] Convective velocity scale
    real zeta 
    real ram       ! [s m-1] Aerodynamic res for mmn xsfer
    real rah       ! [s m-1] Aerodynamic res for sensible heat xsfer
    real raw       ! [s m-1] Aerodynamic res for water vapor xsfer
    real raih
    real raiw
    !real t_ref2m   ! [K] 2 m height specific humidity
    !real q_ref2m   ! [kg kg-1] 2 m height specific humidity
    ! Parameters
    real   ,parameter::zii                = 1000. ! Convective abl hght
    real   ,parameter::rgh_mmn_snw        = 0.0024! Snow roughness len
    real   ,parameter::rgh_mmn_glc_sl_wtl = 0.01  ! soil, glaciers and wetlands
    integer,parameter::niters             = 3     ! max nbr of iterations
    !---------------------------------------------------
    ! Determine surface fluxes based on roughness length
    !   sfc wind speed, sfc temp
        
    ! Potential temperature
    dtheta      = pfl_typ(lev_nbr)%tpt_ptn - sfc_typ%tpt_skn
    theta_avg   = (pfl_typ(lev_nbr)%tpt_ptn + sfc_typ%tpt_skn)/2.
        
    ! Assumption: sfc spc hmd is 1.05*mid-layer spc hmd
    dqdz  = pfl_typ(lev_nbr)%spc_hmd - sfc_typ%q_H2O_sfc
    ! Determine ref hgt and sfc theta_v gradient
    ! diff of theta_v between ref. hgt and sfc
    dthetav   = dtheta*(1.+0.61*pfl_typ(lev_nbr)%spc_hmd)+ &
         (0.61*theta_avg*dqdz)
    thetav = pfl_typ(lev_nbr)%tpt_ptn*(1.+0.61*pfl_typ(lev_nbr)%spc_hmd)
    
    ! Roughness lengths 
    ! for soil, glaciers and wetlands=0.01 snow=0.0024
    rgh_mmn = rgh_mmn_glc_sl_wtl
    ! Displacement height = 0 for bare ground
    displa  = 0.
    zldis   = pfl_typ(lev_nbr)%alt - displa
    
    call MoninObukovInit(dthetav,pfl_typ(lev_nbr)%wnd_spd,zldis,thetav, &
         pfl_typ(lev_nbr)%grv,rgh_mmn,um,obu_lng)
    
    ! Perform stability iteration to converge on
    ! ustar, qstar and thetastar values
    do iter = 1, niters
           call get_scl_prm(pfl_typ(lev_nbr)%alt,displa, rgh_mmn, rgh_heat, rgh_vpr, &
               obu_lng, um, ustar, &
               temp1, temp2, temp12m, temp22m)
           
           theta_star = temp1*dtheta
           qstar      = temp2*dqdz
           rgh_heat   = rgh_mmn/exp(0.13 * (ustar*rgh_mmn/1.5e-5)**0.45)
           rgh_vpr    = rgh_heat   ! Temp and Humidity Rough lgths are the same
           
           thvstar    = theta_star*(1.+0.61*pfl_typ(lev_nbr)%spc_hmd) + 0.61*pfl_typ(lev_nbr)%tpt_ptn*qstar
           zeta       = zldis*SHR_CONST_KARMAN*pfl_typ(lev_nbr)%grv*thvstar/(ustar**2*thetav)
           if (zeta >= 0.) then !stable
              zeta = min(2.,max(zeta,0.01))
              um   = max(pfl_typ(lev_nbr)%wnd_spd,0.1)
           else                 !unstable
              zeta = max(-100.,min(zeta,-0.01))
              wstar   = (-pfl_typ(lev_nbr)%grv*ustar*thvstar*zii/thetav)**0.333
              um   = sqrt(pfl_typ(lev_nbr)%wnd_spd*pfl_typ(lev_nbr)%wnd_spd + &
                   wstar*wstar)
           end if
           obu_lng = zldis/zeta
    end do ! END stability iteration
    
    ram  = 1./(ustar*ustar/um) ! Aerodyn resistance for mmn transfer
    rah  = 1./(temp1*ustar)    ! Thermal resistance
    raw  = 1./(temp2*ustar)    ! Water vapor resistance
    raih = pfl_typ(lev_nbr)%dns_mst_air*SHR_CONST_CPDAIR/rah
    raiw = pfl_typ(lev_nbr)%dns_mst_air/raw
    
    ! Output to pft-level data structures
    ! Derivative of fluxes with respect to ground temperature
    !cgrnds(p) = raih
    !cgrndl(p) = raiw*dqgdT(c)
    
    ! Surface fluxes of momentum, sensible and latent heat
    ! using ground temperatures from previous time step
    sfc_typ%taux_bare = -pfl_typ(lev_nbr)%dns_mst_air*pfl_typ(lev_nbr)%wnd_u/ram
    sfc_typ%tauy_bare = -pfl_typ(lev_nbr)%dns_mst_air*pfl_typ(lev_nbr)%wnd_v/ram
    sfc_typ%tau_bare  = sqrt(sfc_typ%taux_bare**2+sfc_typ%tauy_bare**2)
    sfc_typ%shf_bare  = -raih*dtheta
    sfc_typ%lhf_bare  = -raiw*dqdz
    
    ! 2 m height air temperature
    !t_ref2m = theta + temp1*dtheta*(1./temp12m - 1./temp1)
    
    ! 2 m height specific humidity
    !q_ref2m = pfl_typ(lev_nbr)%spc_hmd + temp2*dqdz*(1./temp22m - 1./temp2)
    return
  end subroutine flx_gnd_get

  subroutine get_scl_prm(alt,displa,rgh_mmn,rgh_heat,rgh_vpr, &
       obu_lng,um, &                  ! I
       ustar, temp1,temp2,temp12m,temp22m) ! O
    implicit none

    real, intent(in)  :: alt     ! [m] Height of measurements
    real, intent(in)  :: displa  ! [m] displacement height
    real, intent(in)  :: rgh_mmn ! [m] roughness length over vegetation, momentum
    real, intent(in)  :: rgh_heat! [m] roughness length over vegetation, sensible heat
    real, intent(in)  :: rgh_vpr ! [m] roughness length over vegetation, latent heat
    real, intent(in)  :: obu_lng ! [m] monin-obukhov length
    real, intent(in)  :: um      ! [m/s] wind speed including the stablity effect
    real, intent(out) :: ustar   ! [m/s] friction velocity
    real, intent(out) :: temp1   ! relation for potential temperature profile
    real, intent(out) :: temp12m ! relation for potential temperature profile applied at 2-m
    real, intent(out) :: temp2   ! relation for specific humidity profile
    real, intent(out) :: temp22m ! relation for specific humidity profile applied at 2-m
    
    real              :: zldis   ! reference height "minus" zero displacement heght [m]
    real              :: zeta    ! dimensionless height used in Monin-Obukhov theory
    real, parameter   :: zetam = 1.574 ! transition point of flux-gradient relation (wind profile)
    real, parameter   :: zetat = 0.465 ! transition point of flux-gradient relation (temp. profile)

    ! Wind profile
    zldis = alt-displa
    zeta  = zldis/obu_lng
    if (zeta < -zetam) then
       ustar = SHR_CONST_KARMAN*um/(log(-zetam*obu_lng/rgh_mmn)&
            - StabilityFunc1(-zetam) &
            + StabilityFunc1(rgh_mmn/obu_lng) &
            + 1.14*((-zeta)**0.333-(zetam)**0.333))
    else if (zeta < 0.) then
       ustar = SHR_CONST_KARMAN*um/(log(zldis/rgh_mmn)&
            - StabilityFunc1(zeta)&
            + StabilityFunc1(rgh_mmn/obu_lng))
    else if (zeta <=  1.) then
       ustar = SHR_CONST_KARMAN*um/(log(zldis/rgh_mmn) + 5.*zeta -5.*rgh_mmn/obu_lng)
    else
       ustar = SHR_CONST_KARMAN*um/(log(obu_lng/rgh_mmn)+5.-5.*rgh_mmn/obu_lng &
            +(5.*log(zeta)+zeta-1.))
    end if
    
    ! Temperature profile
    zldis = alt-displa
    zeta = zldis/obu_lng
    if (zeta < -zetat) then
       temp1 = SHR_CONST_KARMAN/(log(-zetat*obu_lng/rgh_heat)&
            - StabilityFunc2(-zetat) &
            + StabilityFunc2(rgh_heat/obu_lng) &
            + 0.8*((zetat)**(-0.333)-(-zeta)**(-0.333)))
    else if (zeta < 0.) then
       temp1 = SHR_CONST_KARMAN/(log(zldis/rgh_heat) &
            - StabilityFunc2(zeta) &
            + StabilityFunc2(rgh_heat/obu_lng))
    else if (zeta <=  1.) then
       temp1 = SHR_CONST_KARMAN/(log(zldis/rgh_heat) + 5.*zeta - 5.*rgh_heat/obu_lng)
    else
       temp1 = SHR_CONST_KARMAN/(log(obu_lng/rgh_heat) + 5. - 5.*rgh_heat/obu_lng &
            + (5.*log(zeta)+zeta-1.))
    end if
    
    ! Humidity profile
    if (rgh_vpr == rgh_heat) then
       temp2 = temp1
    else
       zldis = alt-displa
       zeta = zldis/obu_lng
       if (zeta < -zetat) then
          temp2 = SHR_CONST_KARMAN/(log(-zetat*obu_lng/rgh_vpr) &
               - StabilityFunc2(-zetat) &
               + StabilityFunc2(rgh_vpr/obu_lng) &
               + 0.8*((zetat)**(-0.333)-(-zeta)**(-0.333)))
       else if (zeta < 0.) then
          temp2 = SHR_CONST_KARMAN/(log(zldis/rgh_vpr) &
               - StabilityFunc2(zeta) &
               + StabilityFunc2(rgh_vpr/obu_lng))
       else if (zeta <=  1.) then
          temp2 = SHR_CONST_KARMAN/(log(zldis/rgh_vpr) + 5.*zeta-5.*rgh_vpr/obu_lng)
       else
          temp2 = SHR_CONST_KARMAN/(log(obu_lng/rgh_vpr) + 5. - 5.*rgh_vpr/obu_lng &
               + (5.*log(zeta)+zeta-1.))
       end if
    endif
    
    ! Temperature profile applied at 2-m
    zldis = 2.0 + rgh_heat
    zeta = zldis/obu_lng
    if (zeta < -zetat) then
       temp12m = SHR_CONST_KARMAN/(log(-zetat*obu_lng/rgh_heat)&
            - StabilityFunc2(-zetat) &
            + StabilityFunc2(rgh_heat/obu_lng) &
            + 0.8*((zetat)**(-0.333)-(-zeta)**(-0.333)))
    else if (zeta < 0.) then
       temp12m = SHR_CONST_KARMAN/(log(zldis/rgh_heat) &
            - StabilityFunc2(zeta)  &
            + StabilityFunc2(rgh_heat/obu_lng))
    else if (zeta <=  1.) then
       temp12m = SHR_CONST_KARMAN/(log(zldis/rgh_heat) + 5.*zeta - 5.*rgh_heat/obu_lng)
    else
       temp12m = SHR_CONST_KARMAN/(log(obu_lng/rgh_heat) + 5. - 5.*rgh_heat/obu_lng &
            + (5.*log(zeta)+zeta-1.))
    end if
    
    ! Humidity profile applied at 2-m
    if (rgh_vpr == rgh_heat) then
       temp22m = temp12m
    else
       zldis = 2.0 + rgh_vpr
       zeta = zldis/obu_lng
       if (zeta < -zetat) then
          temp22m = SHR_CONST_KARMAN/(log(-zetat*obu_lng/rgh_vpr) - &
               StabilityFunc2(-zetat) + StabilityFunc2(rgh_vpr/obu_lng) &
               + 0.8*((zetat)**(-0.333)-(-zeta)**(-0.333)))
       else if (zeta < 0.) then
          temp22m = SHR_CONST_KARMAN/(log(zldis/rgh_vpr) - &
               StabilityFunc2(zeta)+StabilityFunc2(rgh_vpr/obu_lng))
       else if (zeta <=  1.) then
          temp22m = SHR_CONST_KARMAN/(log(zldis/rgh_vpr) + 5.*zeta-5.*rgh_vpr/obu_lng)
       else
          temp22m = SHR_CONST_KARMAN/(log(obu_lng/rgh_vpr) + 5. - 5.*rgh_vpr/obu_lng &
               + (5.*log(zeta)+zeta-1.))
       end if
    end if
    return
  end subroutine get_scl_prm
  
  real function StabilityFunc1(zeta)
    implicit none

    real, intent(in) :: zeta  ! dimensionless height used in Monin-Obukhov theory
    real :: chik, chik2
        
    chik2 = sqrt(1.-16.*zeta)
    chik = sqrt(chik2)
    StabilityFunc1 = 2.*log((1.+chik)*0.5) &
         + log((1.+chik2)*0.5)-2.*atan(chik)+SHR_CONST_PI*0.5
  end function StabilityFunc1
  
  real function StabilityFunc2(zeta)
    implicit none

    real, intent(in) :: zeta  ! dimensionless height used in Monin-Obukhov theory
    real :: chik2
    
    chik2 = sqrt(1.-16.*zeta)
    StabilityFunc2 = 2.*log((1.+chik2)*0.5)
  end function StabilityFunc2
  
  subroutine MoninObukovInit(dthetav,wnd_spd,zldis,thetav,grv,rgh_mmn,um,obu_lng)
    implicit none
    
    real, intent(in)  :: wnd_spd ! wind speed at reference height [m/s]
    real, intent(in)  :: thetav  ! virtual potential temperature (kelvin)
    real, intent(in)  :: dthetav ! diff of vir. poten. temp. between ref. height and surface
    real, intent(in)  :: zldis   ! reference height "minus" zero displacement heght [m]
    real, intent(in)  :: grv     ! lowest level g
    real, intent(in)  :: rgh_mmn ! roughness length, momentum [m]
    real, intent(out) :: um      ! wind speed including the stability effect [m/s]
    real, intent(out) :: obu_lng ! Monin-Obukhov length (m)

    real,parameter :: wstar = 0.5  ! convective velocity [m/s]
    real :: rib   ! bulk Richardson number
    real :: zeta  ! dimensionless height used in Monin-Obukhov theory
    real,parameter :: ustar = 0.06 ! friction velocity [m/s]

    ! Initial guess for the Monin-Obukhov length L
    ! using ustar=0.06 and convective velocity scale = 0.5
    if (dthetav >= 0.) then ! Stable
       um=max(wnd_spd,0.1)
    else                 ! unstable
       um=sqrt(wnd_spd*wnd_spd+wstar*wstar)
    endif

    rib=grv*zldis*dthetav/(thetav*um*um)
    ! See Arya 2001 pg. 220
    if (rib >= 0.) then      ! neutral or stable
       zeta = rib*log(zldis/rgh_mmn)/(1.-5.*min(rib,0.19))
       zeta = min(2.,max(zeta,0.01))
    else                     ! unstable
       zeta = rib*log(zldis/rgh_mmn)
       zeta = max(-100.,min(zeta,-0.01 ))
    endif

    obu_lng=zldis/zeta
    return
  end subroutine MoninObukovInit

end module flx_gnd_mdl
