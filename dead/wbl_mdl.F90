! $Id$

! Purpose: Utilities for Weibull PDF of wind speeds
! Author: Alf Grini, alf.grini@geofysikk.uio.no
! References:
! C. G. Justus et al., J. Appl. Met., 17(3), 350-354 (JHM78)
! C. S. Zender, Natural Aerosols in the Climate System, http://dust.ess.uci.edu/facts/aer/aer.pdf (Zen01c)
! Notes: 
! JHM78 & Zen01c use "c" and "k", respectively, for Weibull scale and shape parameters
! History:
! 20030120 A. Grini  Original Version
! 20030218 C. Zender Clean up
! 20030427 C. Zender Patch to work for wind speed = 0.0, more clean up
! Usage: 
! use wbl_mdl ! [mdl] Weibull wind speed distribution

module wbl_mdl ! [mdl] Weibull wind speed distribution
  use shr_kind_mod,only:r8=>shr_kind_r8 ! [mdl] Precision r8, i8, ...
  implicit none
  private ! [stt] Symbols are private unless individually qualified as public
  private::wbl_prm_rfr_get ! [sbr] Compute Weibull c & k parameters at reference height
  private::wbl_rfr2mdp ! [sbr] Move Weibull parameters c & k to midpoint height
  private::wnd_mdp_wgt_wbl_get ! [sbr] Determine weights of discrete wind bins
  private::wnd_min_max_wbl_get ! [sbr] Determin minimum and maximum discretized wind speeds
  public::wbl_wnd ! [sbr] Generate Weibull wind speed PDF abscissae and weights
  public::wnd_mdp_wbn_get ! [sbr] Determine nominal wind speed in bin
  
  real(r8),parameter::wnd_min_wbl=0.001 ! [m s-1] Minimum windspeed used for Weibull definitions to prevent divide by zero

contains
  
  subroutine wbl_wnd( & ! [sbr] Generate Weibull wind speed PDF abscissae and weights
       hgt_mdp, & ! I [m] Height of layer midpoint
       hgt_rfr, & ! I [m] Reference height
       wnd_frc_rsl, & ! I [frc] Fraction of wind PDF to resolve
       wnd_mdp, & ! I [m s-1] Wind speed 
       wnd_mdp_max, & ! O [m s-1] Maximum discretized wind speed  
       wnd_mdp_min, & ! O [m s-1] Minimum discretized wind speed
       wnd_mdp_nbr, & ! I [nbr] Number of wind speed bins
       wnd_mdp_wgt, & ! O [frc] Wind speed bin weight
       wnd_rfr) ! I [m s-1] Wind speed at reference height
    ! Purpose: Generate Weibull wind speed PDF abscissae and weights
    ! wbl_wnd() is called by dst_mbl()
    ! Discretize Weibull wind speed PDF into wbl_wnd_nbr bins with wbl_wnd_nbr+1 limits
    use pmgrid,only:plon,plond ! [mdl] Spatial resolution parameters
    implicit none
    ! Input
    real(r8),intent(in)::hgt_mdp(plond) ! [m] Midpoint height
    real(r8),intent(in)::hgt_rfr ! [m] Reference height
    real(r8),intent(in)::wnd_frc_rsl ! [frc] Fraction of wind PDF to resolve
    integer,intent(in)::wnd_mdp_nbr ! [nbr] Number of wind speed bins
    real(r8),intent(in)::wnd_rfr(plond) ! [m s-1] Wind speed at reference height
    real(r8),intent(in)::wnd_mdp(plond) ! [m s-1] Wind speed 
    
    ! Output
    real(r8),intent(out)::wnd_mdp_max(plond) ! [m s-1] Maximum discretized wind speed
    real(r8),intent(out)::wnd_mdp_min(plond) ! [m s-1] Minimum discretized wind speed
    real(r8),intent(out)::wnd_mdp_wgt(plond,wnd_mdp_nbr) ! [frc] Wind speed bin weight
    
    ! Local
    real(r8)::wbl_prm_scl_mdp(plond) ! [m s-1] Weibull scale parameter (midpoint)
    real(r8)::wbl_prm_scl_rfr(plond) ! [m s-1] Weibull scale parameter (reference)
    real(r8)::wbl_prm_shp_mdp(plond) ! [frc] Weibull shape parameter (midpoint)
    real(r8)::wbl_prm_shp_rfr(plond) ! [frc] Weibull shape parameter (reference)
    
    ! Compute Weibull c & k parameters at reference height
    call wbl_prm_rfr_get( &
         wbl_prm_scl_rfr, & ! O [m s-1] Weibull scale parameter (reference)
         wbl_prm_shp_rfr, & ! O [frc] Weibull shape parameter (reference)
         wnd_rfr) ! I [m s-1] Wind speed at reference
    
    ! Move Weibull c & k parameters from reference to midpoint height
    call wbl_rfr2mdp( &
         hgt_mdp, & ! I [m] Midpoint height
         hgt_rfr, & ! I [m] Reference height
         wbl_prm_scl_mdp, & ! O [m s-1] Weibull scale parameter (midpoint)
         wbl_prm_scl_rfr, & ! I [m s-1] Weibull scale parameter (reference)
         wbl_prm_shp_mdp, & ! O [frc] Weibull shape parameter (midpoint)
         wbl_prm_shp_rfr, & ! I [frc] Weibull shape parameter (reference)
         wnd_mdp) ! I [m s-1] Wind speed
    
    ! Determine minimum and maximum discretized wind speeds at layer midpoint
    call wnd_min_max_wbl_get( &
         wbl_prm_scl_mdp, & ! I [m s-1] Weibull scale parameter (midpoint)
         wbl_prm_shp_mdp, & ! I [frc] Weibull shape parameter (midpoint)
         wnd_frc_rsl, & ! I [frc] Fraction of wind PDF to resolve
         wnd_mdp_max, & ! O [m s-1] Maximum discretized wind speed
         wnd_mdp_min) ! O [m s-1] Minimum discretized wind speed
    
    ! Determine weights of discrete wind bins
    call wnd_mdp_wgt_wbl_get( &
         wbl_prm_scl_mdp, & ! I [m s-1] Weibull scale parameter (midpoint)
         wbl_prm_shp_mdp, & ! I [frc] Weibull shape parameter (midpoint)
         wnd_mdp_max, & ! I [m s-1] Maximum discretized wind speed
         wnd_mdp_min, & ! I [m s-1] Minimum discretized wind speed
         wnd_mdp_nbr, & ! I [nbr] Number of wind speed bins
         wnd_mdp_wgt) ! O [frc] weight given to each wind speed "bin"
    
  end subroutine wbl_wnd
  
  subroutine wbl_prm_rfr_get( & ! [sbr] Compute Weibull c & k parameters at reference height
       wbl_prm_scl_rfr, & ! O [m s-1] Weibull scale parameter
       wbl_prm_shp_rfr, & ! O [frc] Weibull shape parameter
       wnd_rfr) ! I [m s-1] Mean wind speed at reference height
    ! Purpose: Compute c & k of Weibull distribution winds at anemometer height (10 m)
    ! Author: Alf Grini (alf.grini@geofysikk.uio.no)
    use pmgrid,only:plon,plond ! [mdl] Spatial resolution parameters
    use gmm_mdl,only:gamma ! [mdl] Gamma function gamma()
    implicit none
    ! Input:
    real(r8),intent(in)::wnd_rfr(plond) ! I [m s-1] Reference height wind speed
    ! Output:
    real(r8),intent(out)::wbl_prm_scl_rfr(plond) ! O [m s-1] Weibull scale parameter
    real(r8),intent(out)::wbl_prm_shp_rfr(plond) ! O [frc] Weibull shape parameter 
    ! Local:
    ! fxm: Adopt relationship from high variability sites? Can it be assumed that time-periods of dustiness are high-variability time periods? Or would this be "double-counting" variabiltity since wind PDF is already accounted for?
    ! real(r8),parameter::cst_wbl_shp_var=0.83_r8 ! Weibull shape-mean wind constant, high variablility sites JHM78 p. 352 (20) 
    real(r8),parameter::cst_wbl_shp_var=0.94_r8 ! Weibull shape-mean wind constant, average variablility sites JHM78 p. 352 (20) 
    ! real(r8),parameter::cst_wbl_shp_var=1.05_r8 ! Weibull shape-mean wind constant, low variablility sites JHM78 p. 352 (20) 
    real(r8)::gmm_prm ! Parameter related to gamma distribution
    integer::lon_idx ! [idx] Counting index for longitude
    
    do lon_idx=1,plon
       ! Weibull shape parameter at reference height
       ! fxm csz: 20030427 ensure non-zero shape parameter
       wbl_prm_shp_rfr(lon_idx)=cst_wbl_shp_var*sqrt(max(wnd_rfr(lon_idx),wnd_min_wbl)) ! [frc] JHM78, p. 352 (20)
       ! write(6,*)'Reference height wind = ',wnd_rfr(lon_idx)
       ! write(6,*)'wbl_prm_shp_rfr = ',wbl_prm_shp_rfr(lon_idx)
       
       ! Get Gamma function of (1+1/k) 
       gmm_prm=gamma((1.0_r8+1.0_r8/wbl_prm_shp_rfr(lon_idx))) ! JHM78, p. 351 (16)
       ! write(6,*)'gmm_prm = ',gmm_prm
       
       ! Get Weibull scale parameter from mean wind speed
       wbl_prm_scl_rfr(lon_idx)=wnd_rfr(lon_idx)/gmm_prm ! [m s-1] JHM78, p. 351 (16)
       ! write(6,*)'wbl_prm_scl_rfr = ',wbl_prm_scl_rfr
    enddo ! end loop over lon
    
    ! write(6,*)'gamma 0.5',gamma(0.5_r8)
    ! write(6,*)'sqrt (pi)',sqrt(3.141592654_r8)
    ! do i=1,24
    !   xx=real(i,r8)/6.0_r8
    !   write(6,*)'gamma',xx,gamma(xx)
    ! enddo
    
  end subroutine wbl_prm_rfr_get
  
  subroutine wbl_rfr2mdp( &
       hgt_mdp, & ! I [m] Midpoint height
       hgt_rfr, & ! I [m] Reference height
       wbl_prm_scl_mdp, & ! O [m s-1] Weibull scale parameter (midpoint)
       wbl_prm_scl_rfr, & ! I [m s-1] Weibull scale parameter (reference)
       wbl_prm_shp_mdp, & ! O [frc] Weibull shape parameter (midpoint)
       wbl_prm_shp_rfr, & ! I [frc] Weibull shape parameter (reference)
       wnd_mdp) ! I [m s-1] wind at layer midpoint
    ! Purpose: 
    ! Interpolate Weibull shape and scale parameters from reference height 
    ! (where they are calculated) to midpoint height (where they are used)
    ! Author: Alf Grini, alf.grini@geofysikk.uio.no
    use pmgrid,only:plon,plond ! [mdl] Spatial resolution parameters
    use gmm_mdl
    implicit none
    ! Input
    real(r8),intent(in)::wbl_prm_scl_rfr(plond) ! I [m s-1] Weibull scale parameter (reference)
    real(r8),intent(in)::hgt_mdp(plond) ! I [m] Height of layer midpoint
    real(r8),intent(in)::hgt_rfr ! I [m] Reference height
    real(r8),intent(in)::wbl_prm_shp_rfr(plond) ! I [frc] Weibull shape parameter (reference)
    real(r8),intent(in)::wnd_mdp(plond) ! I [m s-1] wind speed at layer midpoint
    ! Output
    real(r8),intent(out)::wbl_prm_scl_mdp(plond) ! O [m s-1] Weibull scale parameter (midpoint)
    real(r8),intent(out)::wbl_prm_shp_mdp(plond) ! O [frc] Weibull distribution shale factor (midpoint)
    ! Local
    integer::i ! [idx] Counting index for longitude
    !    real(r8)::n_wbl ! [-] JHM78 p. 351 (9)
    real(r8) CEWI_flt ! CEWI
    
    ! Main Code
    CEWI_flt=sum(wbl_prm_scl_rfr) ! CEWI
    
    do i=1,plon
       ! JHM78 Method A proposes
       ! n_wbl= & ! JHM78 p. 351 (9)
       ! (0.37_r8-0.088_r8*log(wbl_prm_scl_rfr(i))) &
       !  /(1.0_r8-0.088_r8*log(hgt_rfr/10.0_r8)) 
       ! ...and...
       ! wbl_prm_scl_mdp(i)=wbl_prm_scl_rfr(i)*(hgt_mdp(i)/hgt_rfr)**n_wbl ! JHM78 p. 351 (8)
       ! ...but Method A over-determines problem since:
       ! 1. We know wind at layer midpoint
       ! 2. We assume wind obeys Weibull distribution
       ! Instead shape of wind distribution to another height as per Method A...
       wbl_prm_shp_mdp(i)= & ! [frc] JHM78 p. 351 (8)
            wbl_prm_shp_rfr(i)*(1.0_r8-0.088_r8*log(hgt_rfr/10.0_r8)) &
            /(1.0_r8-0.088_r8*log(hgt_mdp(i)/10.0_r8))
       
       ! ...then use Method C which relates scale parameter to mean wind and shape parameter
       ! fxm csz: 20030427 ensure non-zero scale parameter to prevent divide by zero later
       wbl_prm_scl_mdp(i)=max(wnd_mdp(i),wnd_min_wbl)/gamma(1.0_r8+1.0_r8/wbl_prm_shp_mdp(i)) ! [m s-1] JHM78 p. 351 (16), Zen01c p. 170 (17.20)
       ! write(6,*)'c mdp /rfr ',wbl_prm_scl_mdp,wbl_prm_scl_rfr
       ! write(6,*)'k mdp /rfr ',wbl_prm_shp_mdp,wbl_prm_shp_rfr
    enddo
    
  end subroutine wbl_rfr2mdp
  
  subroutine wnd_min_max_wbl_get( &
       wbl_prm_scl_mdp, & ! I [m s-1] Weibull scale parameter
       wbl_prm_shp_mdp, & ! I [frc] Weibull shape parameter
       wnd_frc_rsl, & ! I [frc] Fraction of wind PDF to resolve
       wnd_mdp_max, & ! O [m s-1] Maximum discretized wind speed
       wnd_mdp_min) ! O [m s-1] Minimum discretized wind speed
    ! Purpose: 
    ! Determine minimum and maximum discretized winds
    ! If, e.g., wnd_frc_rsl=0.95, the maxiumu wind is the wind speed
    ! at which 95% of winds speeds are smaller.
    ! Similarly, the minimum wind speed would be greater than 5% of the wind PDF
    ! Thus wnd_frc_rsl is a number between zero and one
    ! Author: Alf Grini, alf.grini@geofysikk.uio.no
    use pmgrid,only:plon,plond ! [mdl] Spatial resolution parameters
    implicit none
    ! Input
    real(r8),intent(in)::wbl_prm_scl_mdp(plond) ! [m s-1] Weibull scale parameter
    real(r8),intent(in)::wbl_prm_shp_mdp(plond) ! [frc] Weibull shape parameter
    real(r8),intent(in)::wnd_frc_rsl ! [frc] Fraction of wind PDF to resolve
    ! Output
    real(r8),intent(out)::wnd_mdp_max(plond) ! [m s-1] Maximum discretized wind speed
    real(r8),intent(out)::wnd_mdp_min(plond) ! [m s-1] Minimum discretized wind speed
    ! Local
    integer::lon_idx ! [idx] Counting index for longitude
    real(r8)::wnd_frc_min_lss ! [frc] Fraction of winds slower than wnd_mdp_min_wbl
    real(r8)::wnd_frc_max_lss ! [frc] Fraction of winds slower than wnd_mdp_max_wbl
    real(r8)::xpn_fct_max ! [frc] Factor (U_max/c)**k in cumulative distribution
    real(r8)::xpn_fct_min ! [frc] Factor (U_min/c)**k in cumulative distribution
    ! Begin code
    wnd_frc_min_lss=1.0_r8-wnd_frc_rsl ! [frc] Fraction of winds slower than wnd_mdp_min_wbl 
    wnd_frc_max_lss=wnd_frc_rsl ! [frc] Fraction of winds slower than wnd_mdp_max_wbl
    
    if(wnd_frc_min_lss >= wnd_frc_max_lss) then
       write(6,*)'ERROR in wnd_min_max_wbl_get(): '
       write(6,*)'Specify fraction larger than 0.5 to resolve'
       write(6,*)'Doing otherwise is no better than using mean wind speed'
       error stop
    endif ! endif err
    
    do lon_idx=1,plon
       ! Cumulative distribution with parameters c & k
       ! p(U < Ut) = 1-exp(-(Ut/c)**k) (where Ut is cutoff wind speed)
       ! Solve for factor (Ut/c)**k first, then invert that to get Ut
       xpn_fct_min=-1.0_r8*log(1.0_r8-wnd_frc_min_lss) ! [frc] Zen01c p. 170 (17.18)
       ! Get the factor (wnd_max/c)**k
       xpn_fct_max=-1.0_r8*log(1.0_r8-wnd_frc_max_lss) ! [frc] Zen01c p. 170 (17.18)
       wnd_mdp_max(lon_idx)= & ! [m s-1] Maximum discretized wind speed
            wbl_prm_scl_mdp(lon_idx)*xpn_fct_max**(1.0_r8/wbl_prm_shp_mdp(lon_idx))
       wnd_mdp_min(lon_idx)= & ! [m s-1] Minimum discretized wind speed
            wbl_prm_scl_mdp(lon_idx)*xpn_fct_min**(1.0_r8/wbl_prm_shp_mdp(lon_idx))
       ! write(6,*)'wnd_mdp_min ',wnd_mdp_min
       ! write(6,*)'wnd_mdp_max ',wnd_mdp_max
    enddo
    
  end subroutine wnd_min_max_wbl_get
  
  subroutine wnd_mdp_wgt_wbl_get( &
       wbl_prm_scl_mdp, & ! I [m s-1] Weibull scale parameter
       wbl_prm_shp_mdp, & ! I [frc] Weibull shape parameter
       wnd_mdp_max, & ! I [m s-1] Maximum discretized wind speed
       wnd_mdp_min, & ! I [m s-1] Minimum discretized wind speed
       wnd_mdp_nbr, & ! I [nbr] Number of wind speed bins
       wnd_mdp_wgt) ! O [frc] Wind speed bin weight
    ! Purpose: Given Weibull shape and scale parameters and minimum and maximum 
    ! wind speeds to discretize, compute weights for wind speeds distribution
    ! Author: Alf Grini (alf.grini@geofysikk.uio.no)
    use pmgrid,only:plon,plond ! [mdl] Spatial resolution parameters
    implicit none
    ! Input
    real(r8),intent(in)::wbl_prm_scl_mdp(plond) ! I [m s-1] Weibull scale parameter (midpoint)
    real(r8),intent(in)::wbl_prm_shp_mdp(plond) ! I [frc] Weibull shape parameter (midpoint)
    integer,intent(in)::wnd_mdp_nbr ! I [nbr] Number of wind speed bins
    real(r8),intent(in)::wnd_mdp_max(plond) ! I [m s-1] Maximum discretized wind speed
    real(r8),intent(in)::wnd_mdp_min(plond) ! I [m s-1] Minimum discretized wind speed
    ! Output
    real(r8),intent(out)::wnd_mdp_wgt(plond,wnd_mdp_nbr) ! O [frc] Wind speed bin weight
    ! Local
    integer::i ! [idx] Counting index for longitude
    integer::wnd_mdp_idx ! [idx] Counting index for wind midpoint
    real(r8)::wgt_ttl ! [frc] Total weight of all bins
    real(r8)::wnd_frc_max_lss ! [frc] Fraction of winds slower than maximum wind speed in bin
    real(r8)::wnd_frc_min_lss ! [frc] Fraction of winds slower than minimum wind speed in bin
    real(r8)::wnd_max_wbn ! [m s-1] Maximum wind speed for which we calculate
    real(r8)::wnd_min_wbn ! [m s-1] Minimum wind in a bin
    real(r8)::wnd_ncr(plond) ! [m s-1] Wind increments
    real(r8)::wnd_srt(plond) ! [m s-1] Start wind speed for a bin
    ! Begin code
    
    ! Initialize
    wnd_srt(:)=wnd_mdp_min(:) ! [m s-1] Minimum discretized wind speed
    ! Wind increments at all longitudes
    wnd_ncr(:)=(wnd_mdp_max(:)-wnd_mdp_min(:))/real(wnd_mdp_nbr) ! [m s-1]
    do wnd_mdp_idx=1,wnd_mdp_nbr
       do i=1,plon
          wnd_min_wbn=wnd_srt(i) ! [m s-1] Minimum wind speed in bin
          wnd_max_wbn=wnd_srt(i)+wnd_ncr(i) ! [m s-1] Maximum wind speed in bin
          ! Fraction of winds slower than minimum wind speed in bin
          wnd_frc_min_lss=1.0_r8-exp(-1.0_r8*(wnd_min_wbn/wbl_prm_scl_mdp(i))**wbl_prm_shp_mdp(i)) ! [frc] Zen01c p. 170 (17.18)
          ! Fraction of winds slower than maximum wind speed in bin
          wnd_frc_max_lss=1.0_r8-exp(-1.0_r8*(wnd_max_wbn/wbl_prm_scl_mdp(i))**wbl_prm_shp_mdp(i)) ! [frc] Zen01c p. 170 (17.18)
          ! Weight of current wind speed bin
          wnd_mdp_wgt(i,wnd_mdp_idx)=wnd_frc_max_lss-wnd_frc_min_lss ! [frc] Wind speed bin weight
          ! New start point for next wind speed bin
          wnd_srt(i)=wnd_max_wbn ! [m s-1] Minimum wind speed in bin
       enddo ! end loop over longitude
       
    enddo ! end loop over winds
    
    ! Normalize discretized weights to sum to one
    do i=1,plon
       wgt_ttl=sum(wnd_mdp_wgt(i,:)) ! [frc] Total weight of all bins
       wnd_mdp_wgt(i,:)=wnd_mdp_wgt(i,:)/wgt_ttl ! [frc] Wind speed bin weight
    enddo
    
  end subroutine wnd_mdp_wgt_wbl_get
  
  subroutine wnd_mdp_wbn_get( & ! [fnc] Determine nominal wind speed in bin
       wnd_mdp_idx, & ! I [idx] Counting index for wind bin
       wnd_mdp_max, & ! I [m s-1] Maximum discretized wind speed 
       wnd_mdp_min, & ! I [m s-1] Minimum discretized wind speed 
       wnd_mdp_nbr, & ! I [nbr] Number of wind speed bins
       wnd_mdp_wbn) ! O [m s-1] Nominal wind speed in bin
    ! Purpose: Return nominal wind speed in bin
    ! Given minimum and maximum discretized wind speeds, number of bins,
    ! and bin index, compute nominal wind speed in bin.
    ! Nominal windspeed in bin is currently defined as mean wind speed in bin
    use pmgrid,only:plon,plond ! [mdl] Spatial resolution parameters
    implicit none
    ! Input
    integer,intent(in)::wnd_mdp_idx ! [idx] Counting index for wind bin
    integer,intent(in)::wnd_mdp_nbr ! [nbr] Number of winds to calculate
    real(r8),intent(in)::wnd_mdp_max(plond) ! [m s-1] Maximum discretized wind speed 
    real(r8),intent(in)::wnd_mdp_min(plond) ! [m s-1] Minimum discretized wind speed 
    ! Output
    real(r8),intent(out)::wnd_mdp_wbn(plond) ! [m s-1] Nominal wind speed in bin
    ! Local 
    integer::i ! [idx] Counting index for longitude
    real(r8)::wnd_rng_wbn ! [m s-1] Wind speed range in current bin
    real(r8)::wnd_rng_ttl ! [m s-1] Total range of discretized wind speeds 
    real(r8)::wnd_max_wbn ! [m s-1] Maximum wind speed in bin
    real(r8)::wnd_min_wbn ! [m s-1] Minimum wind speed in bin
    ! Begin code
    do i=1,plon
       wnd_rng_ttl=wnd_mdp_max(i)-wnd_mdp_min(i) ! [m s-1] Total range of discretized wind speeds
       wnd_rng_wbn=wnd_rng_ttl/real(wnd_mdp_nbr,r8) ! [m s-1] Wind speed range in current bin
       wnd_min_wbn=real(wnd_mdp_idx-1,r8)*wnd_rng_wbn+wnd_mdp_min(i) ! [m s-1] Minimum wind speed in bin
       wnd_max_wbn=wnd_min_wbn+wnd_rng_wbn ! [m s-1] Maximum wind speed in bin
       wnd_mdp_wbn(i)=0.5_r8*(wnd_min_wbn+wnd_max_wbn) ! [m s-1] Nominal wind speed in bin
    enddo
    
  end subroutine wnd_mdp_wbn_get
  
end module wbl_mdl ! [mdl] Weibull wind speed distribution
