!-----------------------------------------------------------------------
!
! !MODULE: abl_typ_mdl
! !Author: Scott Capps
!          The University of California, Irvine
! !Date:   11/2006
!
! !DESCRIPTION:	Module to declare derived types and their respective
!        maintenance routines.  ABL turbulence and surface level  
!        variables are declared and stored here.  
!	
!-----------------------------------------------------------------------
module abl_typ_mdl

  implicit none
  
  !----------------------------------------------------------------------------
  ! physical constants (all data public)
  !----------------------------------------------------------------------------
  real    ,parameter :: SHR_CONST_PI     = 3.14159265358979323846  ! pi
  real    ,parameter :: SHR_CONST_BOLTZ  = 1.38065e-23  ! Boltzmann's constant ~ J/K/molecule
  real    ,parameter :: SHR_CONST_AVOGAD = 6.02214e26   ! Avogadro's number ~ molecules/kmole
  real    ,parameter :: SHR_CONST_RGAS   = SHR_CONST_AVOGAD*SHR_CONST_BOLTZ ! Universal gas constant ~ J/K/kmole
  real    ,parameter :: SHR_CONST_MWDAIR = 28.966       ! molecular weight dry air ~ kg/kmole
  real    ,parameter :: SHR_CONST_MWWV   = 18.016       ! molecular weight water vapor
  real    ,parameter :: SHR_CONST_RDAIR  = SHR_CONST_RGAS/SHR_CONST_MWDAIR  ! Dry air gas constant ~ J/K/kg
  real    ,parameter :: SHR_CONST_RWV    = SHR_CONST_RGAS/SHR_CONST_MWWV    ! Water vapor gas constant ~ J/K/kg
  real    ,parameter :: SHR_CONST_ZVIR   = (SHR_CONST_RWV/SHR_CONST_RDAIR)-1.0   ! RWV/RDAIR - 1.0
  real    ,parameter :: SHR_CONST_KARMAN = 0.4          ! Von Karman constant
  real    ,parameter :: SHR_CONST_CPDAIR = 1.00464e3    ! specific heat of dry air ~ J/kg/K
  real    ,parameter :: SHR_CONST_CPWV   = 1.810e3      ! specific heat of water vap ~ J/kg/K
  real    ,parameter :: SHR_CONST_LATVAP = 2.501e6      ! latent heat of evaporation ~ J/kg
  ! ABL
  integer ,parameter :: ntop_turb = 1      ! Top lvl to which turblnt vertical diffusion is applied.
  real    ,parameter :: pblmaxp   = 4.e4   ! pbl max depth in pressure units
  real    ,parameter :: tiny      = 1.e-36 ! lower bound for wind magnitude
  real    ,parameter :: fac       = 100.   ! ustar parameter in height diagnosis
  real    ,parameter :: ricr      =  0.3   ! Critical richardson number

  real    ,parameter :: onet  = 1./3. ! 1/3 power in wind gradient expression
  real    ,parameter :: betam = 15.0  ! Constant in wind gradient expression
  real    ,parameter :: betas =  5.0  ! Constant in surface layer gradient expression
  real    ,parameter :: betah = 15.0  ! Constant in temperature gradient expression 
  real    ,parameter :: fakn  =  7.2  ! Constant in turbulent prandtl number
  real    ,parameter :: fak   =  8.5  ! Constant in surface temperature excess         
  real    ,parameter :: sffrac=  0.1  ! Surface layer fraction of boundary layer
  real    ,parameter :: binm  = betam*sffrac ! betam * sffrac
  real    ,parameter :: binh  = betah*sffrac ! betah * sffrac
  real    ,parameter :: zkmin = 0.01  ! Minimum kneutral*f(ri)

  ! Public routines
  public init_var
  public init_clm_data
  ! Public dervied types
  public sfc_var
  public pfl_mdp
  public pfl_ntf
  public abl_var

  ! Surface flux and state variables
  type sfc_var
     real     tpt_skn ! [K] Ground Temperature
     real     prs_sfc ! [Pa] Surface Pressure
     real     khfs    ! [mK/s]  Surface kinematic heat flux
     real     kbfs    ! [m2/s3] Surface kinematic buoyancy flux
     real     kqfs    ! [m s-1] Kinematic moisture flux
     real     shf_ocn
     real     lhf_ocn
     real     tau_ocn
     real     taux_ocn
     real     tauy_ocn
     real     shf_bare
     real     lhf_bare
     real     tau_bare
     real     taux_bare
     real     tauy_bare 
     real     ustar   ! [m/s] Friction velocity
     real     obk_len ! [m] Obukhov length
     real     wstar   ! [m/s] Convective velocity scale
     real     oro     ! Surface type flag (0=ocn, 1=lnd, 2=sea ice)
     real     q_H2O_sfc ! [kg kg-1] Surface H2O mixing ratio
     real     u10n    ! 10-m surface wind speed magnitude
  end type sfc_var
  
  ! Profile state variables at layer mid-point
  type pfl_mdp
     real     grv         ! [m s-2] Midlayer Acceleration Due to Gravity
     real     prs         ! [Pa] Mid-layer pressure
     real     dns_mst_air ! [kg m-3] Density of moist air
     real     alt
     real     tpt
     real     tpt_ptn     ! [C] Potential temp
     real     wnd_dir
     real     wnd_spd
     real     wnd_shr
     real     wnd_shr2    ! [] Wind shear squared
     real     tpt_ptn_vrt
     real     wnd_u
     real     wnd_v
     real     ri_nbr      ! [] Bulk Richardson number
     real     brt_vsl_frq2!    Brunt-Vaisaila frequency squared
     real     spc_hmd     ! [kg kg-1] Specific humidity
  end type pfl_mdp
  
  ! Profile state variables at interface
  type pfl_ntf
     real     kvf  ! coefficient for heat and tracers
     real     kvm  ! eddy diffusivity for momentum [m2/s]
     real     kvh  ! eddy diffusivity for heat [m2/s]
     real     cgh  ! counter-gradient term for heat [J/kg/m]
     real     cgs  ! counter-gradient star (cg/flux)
     real     tke  ! [m2/s2] Turbulence kinetic energy
  end type pfl_ntf
  
  type abl_var
     real     abl_hgt  ! ABL height
     integer  abl_nbr  ! Maximum number of levels in ABL from surface
     integer  nbot_turb! Bot lvl to which turbulent vrtcl diff is applied
  end type abl_var

contains
  !======================================================================
  ! PUBLIC ROUTINES: Following routines are publically accessable
  !======================================================================
  subroutine init_var(typ_in,typ_in2,typ_in3,typ_in4)
    implicit none
    
    type(pfl_mdp), intent(inout) :: typ_in(:)
    type(sfc_var), intent(inout) :: typ_in2
    type(abl_var), intent(inout) :: typ_in3
    type(pfl_ntf), intent(inout) :: typ_in4(:)
    
    typ_in%grv        = 0.  
    typ_in%prs        = 0.
    typ_in%dns_mst_air= 0.
    typ_in%alt        = 0.
    typ_in%tpt        = 0.
    typ_in%tpt_ptn    = 0.
    typ_in%wnd_dir    = 0.
    typ_in%wnd_spd    = 0.
    typ_in%wnd_shr    = 0.
    typ_in%wnd_shr2   = 0.
    typ_in%tpt_ptn_vrt= 0.
    typ_in%wnd_u      = 0.
    typ_in%wnd_v      = 0.
    typ_in%ri_nbr     = 0.
    typ_in%brt_vsl_frq2 = 0.
    typ_in%spc_hmd    = 0.
    
    typ_in2%tpt_skn   = 0. ! [K] Ground Temperature
    typ_in2%prs_sfc   = 0. ! [Pa] Surface Pressure
    typ_in2%khfs      = 0. ! [mK/s]  Surface kinematic heat flux
    typ_in2%kbfs      = 0. ! [m2/s3] Surface kinematic buoyancy flux
    typ_in2%kqfs      = 0. ! [kg m-2 s-1] Moisture flux
    typ_in2%shf_ocn   = 0.
    typ_in2%lhf_ocn   = 0.
    typ_in2%tau_ocn   = 0.
    typ_in2%taux_ocn  = 0.
    typ_in2%tauy_ocn  = 0.
    typ_in2%shf_bare  = 0.
    typ_in2%lhf_bare  = 0.
    typ_in2%tau_bare  = 0.
    typ_in2%taux_bare = 0.
    typ_in2%tauy_bare = 0.
    typ_in2%ustar     = 0.  ! [m/s] Friction velocity
    typ_in2%obk_len   = 0.
    typ_in2%wstar     = 0.
    typ_in2%oro       = 0.
    typ_in2%q_H2O_sfc = 0.
    typ_in2%u10n      = 0.
    
    typ_in3%abl_hgt   = 0.
    typ_in3%abl_nbr   = 0
    typ_in3%nbot_turb = 0

    typ_in4%kvf       = 0.
    typ_in4%kvm       = 0.
    typ_in4%kvh       = 0.
    typ_in4%cgh       = 0.
    typ_in4%cgs       = 0.
    typ_in4%tke       = 0.
    return
  end subroutine init_var

  subroutine init_clm_data(lev_nbr,oro,typ_in_pfl,typ_in_sfc,alt,tpt,tpt_ptn, &
       tpt_ptn_vrt,grv,q_H2O,prs,dns_mst_air,wnd_dir,wnd_spd,prs_sfc,tpt_skn, &
       q_H2O_sfc)
    implicit none
    ! ONLY initialize input data from CLM
    type(pfl_mdp), intent(inout) :: typ_in_pfl(:)
    type(sfc_var), intent(inout) :: typ_in_sfc
    
    integer,intent(in)::lev_nbr
    real   ,intent(in)::oro              ! [flg] (0=ocn, 1=lnd, 2=sea ice)
    real   ,intent(in)::alt(lev_nbr)     ! [m] Midlayer Height of Measurement
    real   ,intent(in)::tpt(lev_nbr)     ! [K] Midlayer Environment Air Temperature
    real   ,intent(in)::tpt_ptn(lev_nbr) ! [K] Midlayer potential temperature
    real   ,intent(in)::tpt_ptn_vrt(lev_nbr) ! [K] Virtual potential temperature
    real   ,intent(in)::grv(lev_nbr)     ! [m s-2] Midlayer Acceleration Due to Gravity
    real   ,intent(in)::q_H2O(lev_nbr)   ! [kg kg-1] Specific humidity
    real   ,intent(in)::prs(lev_nbr)     ! [Pa] Mid-layer pressure
    real   ,intent(in)::dns_mst_air(lev_nbr)! [kg m-3] Density of moist air
    real   ,intent(in)::wnd_dir(lev_nbr) ! [deg] Wind direction
    real   ,intent(in)::wnd_spd(lev_nbr) ! [m s-1] Wind speed 
    real   ,intent(in)::prs_sfc          ! [Pa] Surface Pressure
    real   ,intent(in)::tpt_skn          ! [K] Ground Temperature
    real   ,intent(in)::q_H2O_sfc        ! [kg kg-1] Surface H2O mixing ratio
    ! Profile variables
    typ_in_pfl%grv        = grv  
    typ_in_pfl%prs        = prs
    typ_in_pfl%dns_mst_air= dns_mst_air
    typ_in_pfl%alt        = alt
    typ_in_pfl%tpt        = tpt
    typ_in_pfl%tpt_ptn    = tpt_ptn
    typ_in_pfl%tpt_ptn_vrt= tpt_ptn_vrt
    typ_in_pfl%wnd_dir    = wnd_dir
    typ_in_pfl%wnd_spd    = wnd_spd
    typ_in_pfl%spc_hmd    = q_H2O
    ! Surface variables
    typ_in_sfc%tpt_skn    = tpt_skn
    typ_in_sfc%prs_sfc    = prs_sfc
    typ_in_sfc%oro        = oro
    typ_in_sfc%q_H2O_sfc  = q_H2O_sfc
    return
  end subroutine init_clm_data
end module abl_typ_mdl
