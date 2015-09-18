!-----------------------------------------------------------------------
!
! !MODULE: abl_trb_mdl
! !Author: Scott Capps
!          The University of California, Irvine
! !Date:   11/2006
!
! !DESCRIPTION:	Calculate atmospheric turbulence variables
!          	
!-----------------------------------------------------------------------
module abl_trb_mdl
  use flx_gnd_mdl ! [mdl] Bareground fluxes
  use flx_ocn_mdl ! [mdl] Ocean surface fluxes
  use abl_typ_mdl ! [mdl] ABL derived types and constants
  use dbg_mdl     ! [mdl] Debugging constants, prg_nm, dbg_lvl
  
  implicit none
  
contains
  
  subroutine abl_trb_pfl(dbg_lvl,lev_nbr,oro,alt,tpt,tpt_ptn,tpt_ptn_vrt, & ! I
       grv,tpt_skn,q_H2O,prs,prs_sfc,dns_mst_air,wnd_dir,wnd_spd,q_H2O_sfc, & ! I
       wnd_shr,ri_nbr,tke,khfs,kbfs,kqfs,shflx,lhflx,ustar, &       ! O
       obk_len,wstar,abl_hgt,u10n,tau) ! O
    implicit none
    
    ! Input
    integer,intent(in)::dbg_lvl          ! [nbr] Debug level
    integer,intent(in)::lev_nbr          ! [nbr] Number of levels in the sounding file
    real   ,intent(in)::oro              ! Surface type flag (0=ocn, 1=lnd, 2=sea ice)
    real   ,intent(in)::alt(lev_nbr)     ! [m] Midlayer Height of Measurement
    real   ,intent(in)::tpt(lev_nbr)     ! [K] Midlayer Environment Air Temperature
    real   ,intent(in)::tpt_ptn(lev_nbr) ! [K] Midlayer potential temperature
    real   ,intent(in)::tpt_ptn_vrt(lev_nbr) ! [K] Virtual potential temperature
    real   ,intent(in)::grv(lev_nbr)     ! [m s-2] Midlayer Acceleration Due to Gravity
    real   ,intent(in)::tpt_skn          ! [K] Surface Temperature
    real   ,intent(in)::q_H2O(lev_nbr)   ! [kg kg-1] Specific humidity
    real   ,intent(in)::prs(lev_nbr)     ! [Pa] Mid-layer pressure
    real   ,intent(in)::prs_sfc          ! [Pa] Surface Pressure
    real   ,intent(in)::dns_mst_air(lev_nbr)! [kg m-3] Density of moist air
    real   ,intent(in)::wnd_dir(lev_nbr) ! [deg] Wind direction
    real   ,intent(in)::wnd_spd(lev_nbr) ! [m s-1] Wind speed
    real   ,intent(in)::q_H2O_sfc        ! [kg kg-1] Surface H2O mixing ratio
    ! Output
    real , intent(out)::wnd_shr(lev_nbr)! [s-1] Wind shear
    real , intent(out)::ri_nbr(lev_nbr) ! [] Richardson number
    real , intent(out)::tke(lev_nbr+1)  ! [m2 s-2] TKE
    real , intent(out)::khfs            ! [mK/s]  Surface kinematic heat flux
    real , intent(out)::kbfs            ! [m2/s3] Surface kinematic buoyancy flux
    real , intent(out)::kqfs            ! [m s-1] Kinematic moisture flux
    real , intent(out)::shflx           ! [W m-2] Sensible heat flux
    real , intent(out)::lhflx           ! [W m-2] Latent heat flux
    real , intent(out)::ustar           ! [m s-1] Friction velocity
    real , intent(out)::obk_len         ! [m] Obukhov length
    real , intent(out)::wstar           ! [m s-1] Convective velocity scale
    real , intent(out)::abl_hgt         ! [m] Boundary layer height
    real , intent(out)::u10n            ! [m/s] 10-m surface wind speed
    real , intent(out)::tau             ! [N/m2] Surface stress
    ! Locals     
    real   , dimension(lev_nbr)::wnd_dir_rad
    ! Derived types
    type(sfc_var)              :: sfc_typ
    type(pfl_mdp), allocatable :: pfl_mdp_typ(:)
    type(abl_var)              :: abl_typ
    type(pfl_ntf), allocatable :: pfl_ntf_typ(:)
    ! -------------------------------------------------
    
    allocate(pfl_mdp_typ(lev_nbr))    ! Mid-layer profile variables
    allocate(pfl_ntf_typ(lev_nbr+1))  ! Interface variables
    call init_var(pfl_mdp_typ,sfc_typ,abl_typ,pfl_ntf_typ)! Initialize all to 0
    
    call init_clm_data(lev_nbr,oro,pfl_mdp_typ,sfc_typ,alt,tpt,tpt_ptn,tpt_ptn_vrt, &
         grv,q_H2O,prs,dns_mst_air,wnd_dir,wnd_spd,prs_sfc,tpt_skn,q_H2O_sfc)
    
    ! Decompose wind magnitude
    ! Convert to radians
    wnd_dir_rad       = pfl_mdp_typ%wnd_dir*SHR_CONST_PI/180.
    pfl_mdp_typ%wnd_u = -pfl_mdp_typ%wnd_spd*sin(wnd_dir_rad)
    pfl_mdp_typ%wnd_v = -pfl_mdp_typ%wnd_spd*cos(wnd_dir_rad)
    
    if (oro==0.0) then
       ! Ocean to atmosphere fluxes
       call flx_ocn_mdl_get(dbg_lvl,pfl_mdp_typ,sfc_typ,lev_nbr)
       if (dbg_lvl > 3) then
          write(6,*)"ocn shf=",sfc_typ%shf_ocn," W m-2; ocn lhf=",sfc_typ%lhf_ocn," W m-2"
          write(6,*)"ocn taux=",sfc_typ%taux_ocn," N m-2; ocn tauy=",sfc_typ%tauy_ocn," N m-2"
       end if
       shflx   = sfc_typ%shf_ocn   ! [W/m2] Sensible heat flux
       lhflx   = sfc_typ%lhf_ocn   ! [W/m2] Latent heat flux
       u10n    = sfc_typ%u10n      ! [m/s] 10-m surface wind speed
       tau     = sfc_typ%tau_ocn   ! 
    elseif (oro==1.0) then
       ! Bare Ground fluxes
       call flx_gnd_get(pfl_mdp_typ,sfc_typ,lev_nbr)
       if (dbg_lvl > 3) then
          write(6,*)"bare shf=",sfc_typ%shf_bare," W m-2; bare lhf=",sfc_typ%lhf_bare," W m-2"
          write(6,*)"bare taux=",sfc_typ%taux_bare," N m-2; bare tauy=",sfc_typ%tauy_bare," N m-2"
       end if
       shflx   = sfc_typ%shf_bare    ! [W/m2] Sensible heat flux
       lhflx   = sfc_typ%lhf_bare    ! [W/m2] Latent heat flux
       tau     = sfc_typ%tau_bare    ! 
    else
       ! TODO: other sfc types
    end if
            
    ! Once the surface fluxes are computed, we can model turbulent transfer within
    ! the boundary layer
    call trb_init(pfl_mdp_typ,sfc_typ,abl_typ,lev_nbr)

    call pblintd(abl_typ,pfl_mdp_typ,sfc_typ,lev_nbr)
        
    call austausch_atm(abl_typ,pfl_mdp_typ,pfl_ntf_typ,lev_nbr)

    call austausch_pbl(abl_typ,sfc_typ,pfl_mdp_typ,pfl_ntf_typ,lev_nbr)

    call phys_shape_get(abl_typ,pfl_mdp_typ,pfl_ntf_typ,sfc_typ,lev_nbr)

    ! Output variables
    wnd_shr = pfl_mdp_typ%wnd_shr
    ri_nbr  = pfl_mdp_typ%ri_nbr
    tke     = pfl_ntf_typ%tke
    khfs    = sfc_typ%khfs        ! [mK/s]  Surface kinematic heat flux
    kbfs    = sfc_typ%kbfs        ! [m2/s3] Surface kinematic buoyancy flux
    kqfs    = sfc_typ%kqfs        ! [m s-1] Kinematic moisture flux
    ustar   = sfc_typ%ustar       ! [m/s] Friction velocity
    obk_len = sfc_typ%obk_len     ! [m] Obukhov length
    wstar   = sfc_typ%wstar       ! [m/s] Convective velocity scale
    abl_hgt = abl_typ%abl_hgt     ! [m] ABL height
    
    return
  end subroutine abl_trb_pfl
  
  subroutine trb_init(pfl_mdp_typ,sfc_typ,abl_typ,lev_nbr)
    implicit none
    ! Input
    integer,       intent(in) :: lev_nbr
    type(pfl_mdp), target,intent(inout), dimension(lev_nbr) :: pfl_mdp_typ
    type(sfc_var), target,intent(inout)                     :: sfc_typ
    type(abl_var), target,intent(inout)                     :: abl_typ
    
    ! Locals ----------------------
    integer  idx        ! [idx] Counting index
    real     dudz,dvdz
    real     rrho_bot   ! [kg m-3] 1./Bottom level air density
    real     dvdz2
    real     dz
    ! Parmeters --------------------
    real    , parameter :: ustar_min = 0.01        ! min permitted value of ustar
    ! Pointers
    real,    pointer :: ustar
    real,    pointer :: oro
    real,    pointer :: taux_ocn
    real,    pointer :: tauy_ocn
    real,    pointer :: tpt(:)
    real,    pointer :: prs(:)
    integer, pointer :: nbot_turb
    real,    pointer :: khfs ! [K m s-1]  Surface kinematic heat flux
    real,    pointer :: kqfs ! [m s-1] Kinematic moisture flux
    real,    pointer :: shf_ocn
    real,    pointer :: lhf_ocn
    real,    pointer :: shf_bare
    real,    pointer :: lhf_bare
    real,    pointer :: taux_bare
    real,    pointer :: tauy_bare
    real,    pointer :: kbfs
    real,    pointer :: tpt_ptn(:)
    real,    pointer :: obk_len ! [m] Obukhov length
    real,    pointer :: tpt_ptn_vrt(:)
    real,    pointer :: grv(:)
    real,    pointer :: wnd_u(:)
    real,    pointer :: wnd_v(:)
    real,    pointer :: alt(:)
    real,    pointer :: wnd_shr2(:)
    real,    pointer :: wnd_shr(:)
    real,    pointer :: ri_nbr(:)
    real,    pointer :: brt_vsl_frq2(:)
    ! Pointer assignments
    ustar       => sfc_typ%ustar
    oro         => sfc_typ%oro
    taux_ocn    => sfc_typ%taux_ocn
    tauy_ocn    => sfc_typ%tauy_ocn
    tpt         => pfl_mdp_typ%tpt
    prs         => pfl_mdp_typ%prs
    nbot_turb   => abl_typ%nbot_turb
    khfs        => sfc_typ%khfs
    kqfs        => sfc_typ%kqfs
    shf_ocn     => sfc_typ%shf_ocn
    shf_bare    => sfc_typ%shf_bare
    lhf_bare    => sfc_typ%lhf_bare
    lhf_ocn     => sfc_typ%lhf_ocn
    taux_bare   => sfc_typ%taux_bare
    tauy_bare   => sfc_typ%tauy_bare 
    kbfs        => sfc_typ%kbfs 
    tpt_ptn     => pfl_mdp_typ%tpt_ptn
    obk_len     => sfc_typ%obk_len
    tpt_ptn_vrt => pfl_mdp_typ%tpt_ptn_vrt
    grv         => pfl_mdp_typ%grv
    wnd_u       => pfl_mdp_typ%wnd_u
    wnd_v       => pfl_mdp_typ%wnd_v
    alt         => pfl_mdp_typ%alt
    wnd_shr2    => pfl_mdp_typ%wnd_shr2
    wnd_shr     => pfl_mdp_typ%wnd_shr
    ri_nbr      => pfl_mdp_typ%ri_nbr
    brt_vsl_frq2=> pfl_mdp_typ%brt_vsl_frq2
    ! END variable declaration -----------------------------
    
    nbot_turb = lev_nbr
        
    rrho_bot       = (SHR_CONST_RDAIR*tpt(lev_nbr))/prs(lev_nbr)
    if (oro==0.0) then
       ustar = max(sqrt(sqrt(taux_ocn**2 + tauy_ocn**2)*rrho_bot),ustar_min)
       khfs  = shf_ocn*rrho_bot/SHR_CONST_CPDAIR ! Kinematic heat flux = F/(rho*cpair)
       kqfs  = lhf_ocn/SHR_CONST_LATVAP ! Moisture flux =LH[W/m2]/l_vap[J/kg]=kg/(m2 s)
    elseif (oro==1.0) then
       ustar = max(sqrt(sqrt(taux_bare**2 + tauy_bare**2)*rrho_bot),ustar_min)
       khfs  = shf_bare*rrho_bot/SHR_CONST_CPDAIR ! Kinematic heat flux = F/(rho*cpair)
       kqfs  = lhf_bare/SHR_CONST_LATVAP ! Moisture flux =LH[W/m2]/l_vap[J/kg]=kg/(m2 s)
    end if
    
    kqfs     = kqfs*rrho_bot ! Kinematic moisture flux [m/s]
    kbfs     = khfs + (0.61*tpt_ptn(lev_nbr)*kqfs)
            
    obk_len  = -tpt_ptn_vrt(lev_nbr)*ustar**3/ &
         (grv(lev_nbr)*SHR_CONST_KARMAN*(kbfs + sign(1.e-10,kbfs)))
    
    do idx=ntop_turb,nbot_turb-1 ! Coming down from TOM
       ! Determine wind shear at each level
       dudz                 = wnd_u(idx) - wnd_u(idx+1)
       dvdz                 = wnd_v(idx) - wnd_v(idx+1)
       dvdz2                = dudz**2 + dvdz**2
       dvdz2                = max(dvdz2,1.e-36)
       dz                   = alt(idx)-alt(idx+1)
       wnd_shr(idx) = sqrt(dvdz2)/dz
       wnd_shr2(idx) = dvdz2/(dz**2)
       ! Calculate Brunt-Vaisala frequency^2:undefined for unstable environs
       pfl_mdp_typ(idx)%brt_vsl_frq2 = grv(lev_nbr)*2.0*(tpt_ptn_vrt(idx)- &
            tpt_ptn_vrt(idx+1))/((tpt_ptn_vrt(idx)+ &
            tpt_ptn_vrt(idx+1))*dz)
       ! Calculate bulk Richardson Number
       ri_nbr(idx) = brt_vsl_frq2(idx)/wnd_shr2(idx)
    end do
    wnd_shr(lev_nbr) = 0. ! No shear btwn bottom level and ground
    ! Write out diagnostics
    if (dbg_lvl > 3) then
       write(6,'(2(a,I3),2(a,F13.6))')' nbot turb=',nbot_turb,' ntop turb=',ntop_turb, &
            ' rrho bot=',rrho_bot,' obk len=',obk_len
       write(6,'(4(a,F13.6))')' ustar=',ustar,' khfs=',khfs,' kqfs=',kqfs, &
            ' kbfs=',kbfs
    end if
    return
  end subroutine trb_init

  subroutine pblintd(abl_typ,pfl_mdp_typ,sfc_typ,lev_nbr)
    implicit none
    ! Input
    integer,       intent(in)::lev_nbr
    type(pfl_mdp), target, intent(inout), dimension(lev_nbr) :: pfl_mdp_typ
    type(sfc_var), target, intent(inout)                     :: sfc_typ
    type(abl_var), target, intent(inout)                     :: abl_typ
    ! Pointers
    real,    pointer :: abl_hgt
    integer, pointer :: abl_nbr
    real,    pointer :: alt(:)
    real,    pointer :: prs(:)
    integer, pointer :: nbot_turb
    real,    pointer :: ustar
    real,    pointer :: wnd_u(:)
    real,    pointer :: wnd_v(:)
    real,    pointer :: tpt_ptn_vrt(:)
    real,    pointer :: grv(:)
    real,    pointer :: kbfs
    real,    pointer :: obk_len
    real,    pointer :: wstar
    
    ! Locals
    real     :: rino(lev_nbr) ! bulk Rich nbr from level to ref lev
    integer  :: idx
    real     :: vvk           ! velocity magnitude squared
    logical  :: check         ! True=>chk if Richardson no.>critcal
    real     :: phiminv       ! inverse phi function for momentum
    real     :: tlv           ! ref. level pot tmp + tmp excess
    
    ! Pointer assignments
    abl_hgt     => abl_typ%abl_hgt
    alt         => pfl_mdp_typ%alt
    prs         => pfl_mdp_typ%prs
    nbot_turb   => abl_typ%nbot_turb
    abl_nbr     => abl_typ%abl_nbr
    ustar       => sfc_typ%ustar
    wnd_u       => pfl_mdp_typ%wnd_u
    wnd_v       => pfl_mdp_typ%wnd_v
    tpt_ptn_vrt => pfl_mdp_typ%tpt_ptn_vrt
    grv         => pfl_mdp_typ%grv
    kbfs        => sfc_typ%kbfs
    obk_len     => sfc_typ%obk_len
    wstar       => sfc_typ%wstar
    
    rino(lev_nbr) = 0.0
    abl_hgt       = alt(lev_nbr)
    
    check = .true.
    tlv   = 0.
    
    ! Limit pbl height to regions below 400 mb
    ! abl_nbr = max number of levels (from bottom) in abl
    do idx=nbot_turb,ntop_turb,-1
       if (prs(idx) >= pblmaxp) then
          abl_nbr = abl_nbr + 1
       end if
    end do
    abl_nbr = max(abl_nbr,1)
    ! Write out diagnostics
    if (dbg_lvl > 3) then
       write(6,'(a,I4,a,F13.6)')'ABL height will be limited to bottom ',abl_nbr, &
         ' levels. Top is ',prs(lev_nbr+1-abl_nbr),' pascals'
    end if
    ! ABL height determination:  Scan upward until the Richardson number between
    ! the first level and the current level exceeds the "critical" value.
    do idx=lev_nbr-1,lev_nbr-abl_nbr+1,-1
       if (check) then
          vvk       = (wnd_u(idx)-wnd_u(lev_nbr))**2+(wnd_v(idx)- &
               wnd_v(lev_nbr))**2+fac*ustar**2
          vvk       = max(vvk,tiny)
          rino(idx) = grv(idx)*(tpt_ptn_vrt(idx) - tpt_ptn_vrt(lev_nbr))* &
               (alt(idx)-alt(lev_nbr))/(tpt_ptn_vrt(lev_nbr)*vvk)
          if (rino(idx) >= ricr) then
             
             abl_hgt= alt(idx+1) + (ricr - rino(idx+1))/(rino(idx) - rino(idx+1)) * &
                  (alt(idx) - alt(idx+1))
             ! Write out diagnostics
             if (dbg_lvl > 3) then
                write(6,'(3(a,F13.6),a,I4)') 'Initial ABL hgt=',abl_hgt, &
                     'Critical RiNum Exceeded=',ricr,' Bulk Num=',rino(idx),' idx=',idx
             end if
             check  = .false.
          end if
       end if
    end do
        
    ! Estimate an effective surface temperature to account for surface fluctuations
    if (check) abl_hgt = alt((lev_nbr+1)-abl_nbr)
    ! Write out diagnostics
    if (dbg_lvl > 3) then
       write(6,'((a,F13.6))')'abl hgt after check=',abl_hgt
    end if
    check  = (kbfs > 0.)
    if (check) then
       phiminv       = (1. - binm*abl_hgt/obk_len)**onet
       rino(lev_nbr) = 0.0
       tlv           = tpt_ptn_vrt(lev_nbr) + kbfs*fak/(ustar*phiminv)
    end if
    
    ! Improve pblh estimate for unstable conditions using the convective temperature excess:
    do idx=lev_nbr-1,lev_nbr-abl_nbr+1,-1
       if (check) then
          vvk = (wnd_u(idx)-wnd_u(lev_nbr))**2+(wnd_v(idx)- &
               wnd_v(lev_nbr))**2+fac*ustar**2
          vvk = max(vvk,tiny)
          rino(idx) = grv(idx)*(tpt_ptn_vrt(idx) - tlv)*(alt(idx)- &
               alt(lev_nbr))/(tpt_ptn_vrt(lev_nbr)*vvk)
          if (rino(idx) >= ricr) then
             abl_hgt = alt(idx+1) + (ricr - rino(idx+1))/(rino(idx) - rino(idx+1))* &
                  (alt(idx) - alt(idx+1))
             if (dbg_lvl > 3) then
                write(6,'(3(a,F13.6),a,I4)') 'ABL hgt improvement=',abl_hgt, &
                  'Critical RiNum Exceeded=',ricr,' Bulk Num=',rino(idx),' idx=',idx
             end if
             check = .false.
          end if
       end if
    end do

    ! PBL height must be greater than some minimum mechanical mixing depth
    ! Several investigators have proposed minimum mechanical mixing depth
    ! relationships as a function of the local friction velocity, u*.  We
    ! make use of a linear relationship of the form h = c u* where c=700.
    ! The scaling arguments that give rise to this relationship most often
    ! represent the coefficient c as some constant over the local coriolis
    ! parameter.  Here we make use of the experimental results of Koracin
    ! and Berkowicz (1988) [BLM, Vol 43] for wich they recommend 0.07/f
    ! where f was evaluated at 39.5 N and 52 N.  Thus we use a typical mid
    ! latitude value for f so that c = 0.07/f = 700.  Also, do not allow 
    ! PBL to exceed some maximum (npbl) number of allowable points
    if (check) abl_hgt = alt((lev_nbr+1)-abl_nbr)
    abl_hgt = max(abl_hgt,700.0*ustar)
    wstar   = (max(0.,kbfs)*grv(idx)*abl_hgt/tpt_ptn_vrt(lev_nbr))**onet
    
    ! Final requirement on PBL height is that it must be greater than the depth
    ! of the lowest model level over ocean if there is any cloud diagnosed in 
    ! the lowest model level.  This is to deal with the inadequacies of the 
    ! current "dry" formulation of the boundary layer, where this test is 
    ! used to identify circumstances where there is marine stratus in the 
    ! lowest level, and to provide a weak ventilation of the layer to avoid
    ! a pathology in the cloud scheme (locking in low-level stratiform cloud)
    ! If over an ocean surface, and any cloud is diagnosed in the 
    ! lowest level, set pblh to 50 meters higher than top interface of lowest level
    !
    !  jrm This is being applied everywhere (not just ocean)!
    !ocncldcheck = .false.
    !if (cldn(pver).ge.0.0) ocncldcheck = .true.
    !if (ocncldcheck) abl_hgt = max(abl_hgt,zi(pver) + 50.)
    
    ! Write out diagnostics
    if (dbg_lvl > 3) then
       write(6,'(a,F13.6)')' wstar=',wstar
       write(6,'(a,F13.6)')' ABL Height=',abl_hgt
    end if
    return
  end subroutine pblintd

  subroutine austausch_atm(abl_typ,pfl_mdp_typ,pfl_ntf_typ,lev_nbr)
    implicit none
    ! Input
    integer,       intent(in)                                :: lev_nbr
    type(pfl_mdp), target, intent(inout), dimension(lev_nbr) :: pfl_mdp_typ
    type(pfl_ntf), target, intent(inout), dimension(lev_nbr+1) :: pfl_ntf_typ
    type(abl_var), target, intent(inout)                     :: abl_typ
    ! Pointers
    real,    pointer :: ri_nbr(:)
    real,    pointer :: wnd_shr2(:)
    integer, pointer :: nbot_turb
    real,    pointer :: kvf(:)
    ! Locals
    integer  :: idx
    real     :: fofri   ! f(ri)
    real     :: kvn     ! neutral Kv
    real     :: ml2(lev_nbr+1) ! Mixing lengths squared

    
    ri_nbr      => pfl_mdp_typ%ri_nbr
    wnd_shr2    => pfl_mdp_typ%wnd_shr2
    nbot_turb   => abl_typ%nbot_turb
    kvf         => pfl_ntf_typ%kvf
    
    ! Set the square of the mixing lengths.
    !
    ml2(ntop_turb) = 0.
    do idx = ntop_turb+1, nbot_turb
       ml2(idx) = 30.0**2
    end do
    ml2(nbot_turb+1) = 0.

    ! Set the vertical diffusion coefficient above the top diffusion level
    ! Note that nbot_turb != pver is not supported
    kvf(1:ntop_turb) = 0.0
    !
    ! Compute the free atmosphere vertical diffusion coefficients: kvh = kvq = kvm. 
    do idx = ntop_turb,nbot_turb-1
       if (ri_nbr(idx) < 0.0) then
          fofri = sqrt(max(1. - 18.*ri_nbr(idx),0.))
       else 
          fofri = 1.0/(1.0 + 10.0*ri_nbr(idx)*(1.0 + 8.0*ri_nbr(idx)))    
       end if
       kvn = ml2(idx)*sqrt(wnd_shr2(idx))
       kvf(idx+1) = max(zkmin,kvn*fofri)
    end do

    return
  end subroutine austausch_atm

  subroutine austausch_pbl(abl_typ,sfc_typ,pfl_mdp_typ,pfl_ntf_typ,lev_nbr)
    implicit none
    ! Nonlocal scheme that determines eddy diffusivities based on a
    ! specified boundary layer height and a turbulent velocity scale;
    ! also, countergradient effects for heat and moisture, and constituents
    ! are included, along with temperature and humidity perturbations which
    ! measure the strength of convective thermals in the lower part of the
    ! atmospheric boundary layer.
    !
    ! For more information, see Holtslag, A.A.M., and B.A. Boville, 1993:
    ! Local versus Nonlocal Boundary-Layer Diffusion in a Global Climate
    ! Model. J. Clim., vol. 6., p. 1825--1842.
    
    ! Input
    integer,       intent(in)::lev_nbr
    type(pfl_mdp), target, intent(inout), dimension(lev_nbr) :: pfl_mdp_typ
    type(pfl_ntf), target, intent(inout), dimension(lev_nbr+1):: pfl_ntf_typ
    type(sfc_var), target, intent(inout)                     :: sfc_typ
    type(abl_var), target, intent(inout)                     :: abl_typ
    logical  :: unstbl ! pts w/unstbl pbl (positive virtual ht flx)
    ! pointers
    real,    pointer :: kbfs
    real,    pointer :: abl_hgt
    integer, pointer :: abl_nbr
    real,    pointer :: ustar
    real,    pointer :: obk_len
    real,    pointer :: wstar
    real,    pointer :: kvf(:)
    real,    pointer :: kvm(:)
    real,    pointer :: kvh(:)
    real,    pointer :: tke(:)
    real,    pointer :: alt(:)
    ! Locals
    integer  :: idx
    real     :: pblk     ! level eddy diffusivity for momentum
    real     :: fak1     ! k*ustar*pblh     
    real     :: fak2     ! k*wm*pblh
    real     :: fak3     ! fakn*wstar/wm
    real     :: phiminv  ! inverse phi function for momentum
    real     :: phihinv  ! inverse phi function for heat
    real     :: wm       ! turbulent velocity scale for momentum
    integer  :: ktopbl   ! index of first midpoint inside pbl
    integer  :: ktopblmn ! min value of ktopbl
    logical  :: pblpt    ! pts within pbl
    real     :: zp       ! current level height + one level up
    real     :: zl       ! zmzp / Obukhov length
    real     :: zh       ! zmzp / pblh
    real     :: zzh      ! (1-(zmzp/pblh))**2
    real     :: zmzp     ! level height halfway between zm and zp
    real     :: term     ! intermediate calculation
    real     :: pr       ! Prandtl number for eddy diffusivities
    real     :: ccon     ! fak * sffrac * vk
    ! Pointer assignments
    abl_hgt     => abl_typ%abl_hgt
    abl_nbr     => abl_typ%abl_nbr
    kbfs        => sfc_typ%kbfs
    ustar       => sfc_typ%ustar
    obk_len     => sfc_typ%obk_len
    wstar       => sfc_typ%wstar
    kvf         => pfl_ntf_typ%kvf
    tke         => pfl_ntf_typ%tke
    kvm         => pfl_ntf_typ%kvm
    kvh         => pfl_ntf_typ%kvh
    alt         => pfl_mdp_typ%alt
    
    ! CEWI (compiler-error-warning-initializer) 
    unstbl = .false.
    fak3   = 0.
    fak2   = 0.
    phihinv= 0.
    phiminv= 0.
    ktopbl = 0
    ! END CEWI
    
    ccon   = fak*sffrac*SHR_CONST_KARMAN
    unstbl = (kbfs > 0.)
    pblk   = 0.0
    fak1   = ustar*abl_hgt*SHR_CONST_KARMAN
    if (unstbl) then
       phiminv = (1. - binm*abl_hgt/obk_len)**onet
       phihinv = sqrt(1. - binh*abl_hgt/obk_len)
       wm      = ustar*phiminv
       fak2    = wm*abl_hgt*SHR_CONST_KARMAN
       fak3    = fakn*wstar/wm
       !tpert   = max(khfs*fak/wm,0.)
       !qpert   = max(kqfs*fak/wm,0.)
    else
       !tpert   = max(khfs*fak/ustar,0.)
       !qpert   = max(kqfs*fak/ustar,0.)
    end if
        
    ! Initialize output arrays with free atmosphere values
    do idx=1,lev_nbr+1
       kvm(idx) = kvf(idx)
       kvh(idx) = kvf(idx)
    end do
    ! Main level loop to compute the diffusivities and counter-gradient terms. These terms are 
    ! only calculated at points determined to be in the interior of the pbl (pblpt(i)==.true.),
    ! and then calculations are directed toward regime: stable vs unstable, surface vs outer 
    ! layer.
    do idx=lev_nbr,lev_nbr-abl_nbr+2,-1
       pblpt = (alt(idx) < abl_hgt)
       if (pblpt) then
          ktopbl = idx
          zp  = alt(idx-1)
          if (zkmin == 0.0 .and. zp > abl_hgt) zp = abl_hgt
          zmzp    = 0.5*(alt(idx) + zp)
          zh   = zmzp/abl_hgt
          zl   = zmzp/obk_len
          zzh  = zh*max(0.,(1. - zh))**2
          if (unstbl) then
             if (zh < sffrac) then
                term     = (1. - betam*zl)**onet
                pblk  = fak1*zzh*term
                pr    = term/sqrt(1. - betah*zl)
             else
                pblk  = fak2*zzh
                pr    = phiminv/phihinv + ccon*fak3/fak
             end if
          else
             if (zl <= 1.) then
                pblk = fak1*zzh/(1. + betas*zl)
             else
                pblk = fak1*zzh/(betas + zl)
             end if
             pr    = 1.
          end if
          kvm(idx) = max(pblk,kvf(idx))
          kvh(idx) = max(pblk/pr,kvf(idx))
       end if
    end do
    
    ! Check whether last allowed midpoint is within pbl, determine ktopblmn
    !
    ktopblmn = lev_nbr
    idx = lev_nbr-abl_nbr+1
    if (alt(idx) < abl_hgt) ktopbl = idx
    ktopblmn = min(ktopblmn, ktopbl)
    
    ! Crude estimate of tke (tke=0 above boundary layer)
    ! kvm = eddy diffusivity for momentum [m2/s] 
    ! abl hgt in meters
    ! TKE [m2/s2]
    do idx = ktopblmn,lev_nbr+1
       !tke(idx) = (kvm(idx)/abl_hgt)**2
       tke(idx) = 340.*((kvm(idx)/abl_hgt)**2)
       if (dbg_lvl > 3) then
          write(6,'(a,F13.6)')'TKE=',tke(idx)
       end if
    end do
    
    return
  end subroutine austausch_pbl

  subroutine phys_shape_get (abl_typ,pfl_mdp_typ,pfl_ntf_typ,sfc_typ,lev_nbr)
    implicit none
    ! Input
    integer,       intent(in) :: lev_nbr
    type(pfl_mdp), target, intent(inout), dimension(lev_nbr)  :: pfl_mdp_typ
    type(pfl_ntf), target, intent(inout), dimension(lev_nbr+1):: pfl_ntf_typ
    type(sfc_var), target, intent(inout)                      :: sfc_typ
    type(abl_var), target, intent(inout)                      :: abl_typ
    ! locals
    integer idx
    real    dz
    real    tke_tot
    real    buoy_tot
    real    wnd_gst(lev_nbr)
    real    wnd_gst_max
    real    ri_nbr_tot
    real    ri_nbr_avg
    real    dz_wgt
    
    ! Pointers
    real,    pointer :: abl_hgt
    real,    pointer :: wstar
    real,    pointer :: alt(:)
    real,    pointer :: grv(:)
    integer, pointer :: nbot_turb
    real,    pointer :: tke(:)
    real,    pointer :: tpt_ptn_vrt(:)
    real,    pointer :: wnd_u(:)
    real,    pointer :: wnd_v(:)
    real,    pointer :: ri_nbr(:)
    
    ! Pointer assignments
    abl_hgt     => abl_typ%abl_hgt
    alt         => pfl_mdp_typ%alt
    grv         => pfl_mdp_typ%grv
    nbot_turb   => abl_typ%nbot_turb
    tke         => pfl_ntf_typ%tke
    tpt_ptn_vrt => pfl_mdp_typ%tpt_ptn_vrt
    wnd_u       => pfl_mdp_typ%wnd_u
    wnd_v       => pfl_mdp_typ%wnd_v
    wstar       => sfc_typ%wstar
    ri_nbr      => pfl_mdp_typ%ri_nbr
    
    tke_tot     = 0.
    buoy_tot    = 0.
    wnd_gst(:)  = 0.
    wnd_gst_max = 0.
    ri_nbr_tot  = 0.

!    ---------- ntf
!
!    ++++++++++ ntop_turb
!
!    ---------- ntf
!    ********** ABL hgt
!    ++++++++++ mid lyr
!
!    ---------- ntf
!
!    ++++++++++ nbot_turb
!
!    ---------- surface
        
    ! From bottom mid-layer level to ABL hgt
    do idx = nbot_turb-1,ntop_turb,-1
              
       if (alt(idx)>abl_hgt) then
          write(6,*)'ABL hgt exceeded, exiting...'
          exit
       end if
       
       dz  = alt(idx-1)-alt(idx)
       ! weight
       dz_wgt = dz/abl_hgt
       if (dbg_lvl > 3) then
          write(6,'(4(a,F13.6))')'height=',alt(idx),' abl_hgt=',abl_hgt, &
               ' lyr thick=',dz,' wgt=',dz_wgt
       end if
       
       ! Integrate TKE below parcel height
       tke_tot = tke_tot + (tke(idx)*dz)
       
       ! Integrate buoyant energy below parcel height
       buoy_tot = buoy_tot + (dz*grv(idx)*(tpt_ptn_vrt(idx-1)- &
            tpt_ptn_vrt(idx))/tpt_ptn_vrt(idx))
       ! Calculate RI_NUM at each level within ABL
       ri_nbr_tot = ri_nbr_tot + (ri_nbr(idx)*dz_wgt)
       
       if (dbg_lvl > 3) then
          write(6,'(2(a,F13.6))')' tke below this lyr=',tke_tot,' buoyancy consumption=',buoy_tot
          write(6,'(2(a,F13.6))')' wind speed this lyr=',wnd_gst(idx),' RI_NUM=',ri_nbr(idx)
       end if
       
       if (tke_tot/alt(idx)>=buoy_tot) then
          wnd_gst(idx) = sqrt(wnd_u(idx)**2 + wnd_v(idx)**2)
       else
          wnd_gst(idx) = 0. ! No wind mixed down to sfc
       end if
       
    end do
    
    write(6,*)'Ending wind variability computation...'
    
    ! Calculate avg ri_nbr for ABL
    ri_nbr_avg = ri_nbr_tot
    
    ! Chose the fastest wind speed throughout the ABL
    wnd_gst_max = maxval(wnd_gst)
    if (dbg_lvl > 3) then
       write(6,'(a,F13.6)')'Max Gust=',wnd_gst_max,' ri_nbr_avg=',ri_nbr_avg, &
            ' wstar=',wstar
    end if
    
    return
  end subroutine phys_shape_get
end module abl_trb_mdl
