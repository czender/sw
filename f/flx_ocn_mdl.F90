!-----------------------------------------------------------------------
!
! !MODULE: flx_ocn_mdl
! !Author: Routines from the CCSM 3.0
!          Modified by: Scott Capps
! !DESCRIPTION:	Calculate ocean/atm fluxes
!          	
!-----------------------------------------------------------------------
module flx_ocn_mdl
  use abl_typ_mdl ! Shared constants

  implicit none
  
contains
  subroutine flx_ocn_mdl_get(dbg_lvl,pfl_typ,sfc_typ,lev_nbr)
    implicit none
    ! Input
    integer, intent(in)::dbg_lvl
    integer, intent(in)::lev_nbr
    type(pfl_mdp), intent(inout), dimension(lev_nbr) :: pfl_typ
    type(sfc_var), intent(inout)                     :: sfc_typ
        
    ! Locals
    real     theta_bot ! [K] Bottom level pot temp
    real     dtheta    ! [K] Pot temp difference
    real     rho_bot   ! [kg m-3] Bottom level air density
    real     wnd_spd   ! [m s-1] Wind speed
    real     delq      ! [kg kg-1] Specific hum difference
    real     thetav_bot! [K] Bottom level virtual pot temp
    real     ssq       ! [kg kg-1] Surface saturation specific humidity
    real     cpvir     ! 
    real     cp_mst    ! Specific heat of moist air
    real     rdn       ! sqrt of neutral exchange coeff (momentum)
    real     rhn       ! sqrt of neutral exchange coeff (heat)
    real     ren       ! sqrt of neutral exchange coeff (tracers)
    real     stable    ! Stability factor
    real     ustar     ! [m s-1] Friction velocity
    real     tstar     ! tstar
    real     qstar     ! qstar
    real     zeta      ! Ref hgt (10m) / monin-obukhov length
    real     psimh     ! Stability funct at ref lev (momentum)
    real     psixh     ! Stability funct at ref lev (heat & tracers) 
    real     rd        ! sqrt of exchange coeff (momentum)
    real     rh        ! sqrt of exchange coeff (heat)
    real     re        ! sqrt of exchange coeff (tracers)
    real     alz       ! ln(zbot/z10)
    real     xsq       ! Temporary variable
    real     xqq       ! Temporary variable
    real     rgh_mmn   ! Momentum roughness length 
    
    ! Parameters
    real    ,parameter :: umin             = 1.0 ! Minimum wind speed at bottom level
    real    ,parameter :: alt_ref          = 10.0! 10m reference height
    
    ! Functions ---------------------------
    real     xd             ! Dummy argument
    real     Tk             ! Temperature (K)
    real     Umps           ! Wind velocity (m/sec)
    real     psimhu         ! Unstable part of psimh
    real     psixhu         ! Unstable part of psixh
    real     qsat           ! Saturation specific humidty of air
    real     cdn            ! Neutral drag coeff at bottom model level
    
    ! Functions ---------------------------
    qsat(Tk)  = 640380. / exp(5107.4/Tk)
    cdn(Umps) = 0.0027 / Umps + .000142 + .0000764 * Umps
    psimhu(xd)= log((1.+xd*(2.+xd))*(1.+xd*xd)/8.) - 2.*atan(xd) + 1.571
    psixhu(xd)= 2. * log((1. + xd*xd)/2.)
    ! -------------------------------------
        
    cpvir = SHR_CONST_CPWV/SHR_CONST_CPDAIR - 1.
    ! Bottom level pot temp
    theta_bot = pfl_typ(lev_nbr)%tpt_ptn
    ! Bottom level density
    rho_bot   = pfl_typ(lev_nbr)%prs/(SHR_CONST_RDAIR*pfl_typ(lev_nbr)%tpt)

    ! Bottom level wind magnitude
    wnd_spd   = max(umin, pfl_typ(lev_nbr)%wnd_spd)
    
    ! Virtual pot temp
    thetav_bot= theta_bot * (1.0 + SHR_CONST_ZVIR*pfl_typ(lev_nbr)%spc_hmd)
    ! Surface saturation specific humidity
    ssq       = 0.98 * qsat(sfc_typ%tpt_skn)/rho_bot
    ! Potential T difference
    dtheta    = theta_bot - sfc_typ%tpt_skn
    ! Specific humidity difference
    delq      = pfl_typ(lev_nbr)%spc_hmd - ssq 
    ! ln(zbot/z10)
    alz       = log(pfl_typ(lev_nbr)%alt/alt_ref)
    ! Specific heat of moist air
    cp_mst     = SHR_CONST_CPDAIR*(1. + cpvir*ssq) 

    ! Diagnostic output
    if (dbg_lvl > 3) then
       write(6,'(3(a,F13.6))')' Bot-layer height=',pfl_typ(lev_nbr)%alt, &
            ' Ref height=',alt_ref,' Bot-layer density =',rho_bot
       write(6,'(5(a,F13.6))')' bot-layer prs=',pfl_typ(lev_nbr)%prs,' bot-layer tpt=',pfl_typ(lev_nbr)%tpt, &
            ' ubot=',pfl_typ(lev_nbr)%wnd_u,' vbot=',pfl_typ(lev_nbr)%wnd_v,' wspd=',wnd_spd
       write(6,'(5(a,F13.6))')' bot-layer theta=',theta_bot,' bot-layer theta_v=', &
            thetav_bot,' bot-layer spc_hmd=',pfl_typ(lev_nbr)%spc_hmd, &
            ' sfc sat spc hmd=',ssq,' Skin temp=',sfc_typ%tpt_skn
       write(6,'(a,F13.6)')'ln(zbot/z10)=',alz
    end if
    !---------------------------------------------------------------
    ! First iteration to converge on Z/L and hence the fluxes
    !---------------------------------------------------------------
    ! Initial guess for roots of neutral exchange coefficients, 
    ! assume z/L=0. and u10n is approximated by vmag.
    ! Stable if (thbot > ts ).
    !
    stable = 0.5 + sign(0.5, dtheta)
    rdn    = sqrt(cdn(wnd_spd))
    rhn    = (1.-stable) * 0.0327 + stable * 0.018 
    ren    = 0.0346 
    !
    ! Initial guess of ustar,tstar and qstar
    !
    ustar = rdn*wnd_spd
    tstar = rhn*dtheta
    qstar = ren*delq
    ! Diagnostic output
    if (dbg_lvl > 3) then
       write(6,'(a)')'Initial guess'
       write(6,'(5(a,F13.6))')' dtheta=',dtheta,' Stability factor=',stable, &
            ' rdn=',rdn,' rhn=',rhn,' ren=',ren
       write(6,'(3(a,F13.6))')' ustar=',ustar,' tstar=',tstar,' qstar=',qstar
    end if

    ! Compute stability and evaluate all stability functions
    ! Stable if (thbot > ts or zeta > 0 )
    !
    zeta = SHR_CONST_KARMAN*pfl_typ(lev_nbr)%grv*pfl_typ(lev_nbr)%alt*(tstar/thetav_bot + &
         qstar/(1./SHR_CONST_ZVIR+pfl_typ(lev_nbr)%spc_hmd)) / ustar**2
    zeta = sign( min(abs(zeta),10.), zeta )
    stable = 0.5 + sign(0.5, zeta)
    xsq   = max(sqrt(abs(1. - 16.*zeta)) , 1.)
    xqq   = sqrt(xsq)
    psimh = -5. * zeta * stable + (1.-stable)*psimhu(xqq)
    psixh = -5. * zeta * stable + (1.-stable)*psixhu(xqq)
    
    ! Shift 10m neutral wind speed using old rdn coefficient
    !
    rd   = rdn / (1.+rdn/SHR_CONST_KARMAN*(alz-psimh))
    sfc_typ%u10n = wnd_spd * rd / rdn
    !
    ! Update the neutral transfer coefficients at 10m and neutral stability
    !
    rdn = sqrt(cdn(sfc_typ%u10n))
    ren = 0.0346
    rhn = (1.-stable) * 0.0327 + stable * 0.018 
    !
    ! Shift all coeffs to measurement height and stability
    !
    rd = rdn / (1.+rdn/SHR_CONST_KARMAN*(alz-psimh)) 
    rh = rhn / (1.+rhn/SHR_CONST_KARMAN*(alz-psixh)) 
    re = ren / (1.+ren/SHR_CONST_KARMAN*(alz-psixh))
    
    ! Update ustar, tstar, qstar using updated, shifted coeffs 
    !
    ustar = rd * wnd_spd 
    tstar = rh * dtheta 
    qstar = re * delq 
    !
    !---------------------------------------------------------------
    ! Second iteration to converge on Z/L and hence the fluxes
    !---------------------------------------------------------------
    !
    ! Recompute stability & evaluate all stability functions  
    ! Stable if (thbot > ts or hol > 0 )
    ! 
    zeta = SHR_CONST_KARMAN*pfl_typ(lev_nbr)%grv*pfl_typ(lev_nbr)%alt*(tstar/thetav_bot + &
         qstar/(1./SHR_CONST_ZVIR+pfl_typ(lev_nbr)%spc_hmd))/ustar**2
    zeta = sign( min(abs(zeta),10.), zeta )
    stable = 0.5 + sign(0.5, zeta)
    xsq   = max(sqrt(abs(1. - 16.*zeta)) , 1.)
    xqq   = sqrt(xsq)
    psimh = -5. * zeta * stable + (1.-stable)*psimhu(xqq)
    psixh = -5. * zeta * stable + (1.-stable)*psixhu(xqq)
    !
    ! Shift 10m neutral wind speed using old rdn coefficient
    !
    rd   = rdn / (1.+rdn/SHR_CONST_KARMAN*(alz-psimh))
    sfc_typ%u10n = wnd_spd * rd / rdn
    rgh_mmn = 10.*exp(-SHR_CONST_KARMAN/cdn(sfc_typ%u10n))
    
    ! Update the neutral transfer coefficients at 10m and neutral stability
    !
    rdn = sqrt(cdn(sfc_typ%u10n))
    ren = 0.0346
    rhn = (1.-stable) * 0.0327 + stable * 0.018 
    !
    ! Shift all coeffs to measurement height and stability
    
    rd = rdn / (1.+rdn/SHR_CONST_KARMAN*(alz-psimh)) 
    rh = rhn / (1.+rhn/SHR_CONST_KARMAN*(alz-psixh)) 
    re = ren / (1.+ren/SHR_CONST_KARMAN*(alz-psixh))
    
    !---------------------------------------------------------------
    ! Compute the fluxes
    !---------------------------------------------------------------
    !
    ! Update ustar, tstar, qstar using updated, shifted coeffs 
    !
    ustar = rd * wnd_spd 
    tstar = rh * dtheta 
    qstar = re * delq 
    !
    ! Compute surface stress components
    !
    sfc_typ%tau_ocn  =  rho_bot * ustar * ustar 
    sfc_typ%taux_ocn = -sfc_typ%tau_ocn * pfl_typ(lev_nbr)%wnd_u / wnd_spd 
    sfc_typ%tauy_ocn = -sfc_typ%tau_ocn * pfl_typ(lev_nbr)%wnd_v / wnd_spd 
    
    ! Compute heat flux components at current surface temperature
    ! (Define positive latent and sensible heat as upwards into the atm)
    !
    sfc_typ%shf_ocn  = -cp_mst * sfc_typ%tau_ocn * tstar / ustar 
    sfc_typ%lhf_ocn  = -SHR_CONST_LATVAP * sfc_typ%tau_ocn * qstar / ustar
    !
    if (dbg_lvl > 3) then
       write(6,'(4(a,F13.6))')' shf=',sfc_typ%shf_ocn,' lhf=',sfc_typ%lhf_ocn, &
            ' ustar=',ustar,' tstar=',tstar
       write(6,'(2(a,F13.6))')' qstar=',qstar,' tau=',sfc_typ%tau_ocn
       write(6,'(2(a,F13.6))')" 10-m wind=",sfc_typ%u10n," Rough len(mmn)=",rgh_mmn
    end if
    return
  end subroutine flx_ocn_mdl_get
end module flx_ocn_mdl
