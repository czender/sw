module wnd_PDF_mdl
  use shr_kind_mod,only:r8=>shr_kind_r8 ! [mdl] Precision r8, i8, ...
  use gmm_mdl,only:gamma,gratio ! [mdl] Gamma function gamma()
  
  implicit none
contains
  
  subroutine wnd_spd_PDF(wnd_spd_avg,wnd_spd_var,bin_nbr,bin_wgt,bin_spd)
    !Computes a PDF for surface wind speed passed in by CLM module via
    !the command-line (default=10m s-1)
    !command-line inputs into CLM=mean wind speed, nbr of bins, wind spd variability
    
    implicit none

    !Input
    real(r8),intent(in)::wnd_spd_avg ! [m s-1] Average wind speed
    real(r8),intent(in)::wnd_spd_var ! [frc] Wind spd variability:low=1.05;avg=0.94;high=0.83 Justus et al.(1978)
    integer, intent(in)::bin_nbr     ! [nbr] Number of wind speed bins
    
    !Output
    real(r8), intent(out)::bin_wgt(bin_nbr) ! [frc] Weight assigned to each bin
    real(r8), intent(out)::bin_spd(bin_nbr) ! [m s-1] Bin average wind speed
    
    !Static
    real(r8)::bin_srt        ! [m s-1] Starting wind speed
    integer   idx            ! [idx]
    real(r8)::wbl_alpha      ! (wbl_shp_prm + 1.0_r8)/wbl_shp_prm
    real(r8)::wbl_shp_prm    ! [frc] Weibull shape parameter determines variance of PDF
    real(r8)::wbl_scl_prm    ! [m s-1] Weibull scale parameter determines max value of PDF
        
    !-----------------------------------------------
    ! Compute equal probability wind speed bin size
    bin_wgt(:) = 1.0_r8/real(bin_nbr,r8)
    
    !Calculate scale and shape parameters
    !The mean wind speed is defined by the first moment of the Weibull Distribution
    !Ubar=W(Ut=0,n=1)=c(gamma(1+1/k)) -> c=Ubar/(gamma(1+1/k))
    if (wnd_spd_avg == 0.0_r8) then
       write(6,*) "Wind speed is zero and will cause a divide by zero condition...exiting"
       error stop
    end if
    ! [frc] Weibull shape parameter
    wbl_shp_prm = wnd_spd_var*sqrt(wnd_spd_avg)
    ! [m s-1] Weibull scale parameter
    wbl_scl_prm = wnd_spd_avg/(gamma(1.0_r8+1.0_r8/wbl_shp_prm))
    wbl_alpha   = (wbl_shp_prm + 1.0_r8)/wbl_shp_prm
    
    bin_srt = 0.0_r8
    do idx=1,bin_nbr
       call get_wbn_spd(idx,bin_srt,bin_nbr,wbl_alpha, & ! I
            wbl_scl_prm,wbl_shp_prm,bin_wgt(1), & ! I
            bin_spd(idx)) ! O
    end do
    return
  end subroutine wnd_spd_PDF
  
  subroutine get_wbn_spd(idx,bin_srt,bin_nbr,wbl_alpha, & ! I
       wbl_scl_prm,wbl_shp_prm,bin_wgt, & ! I
       bin_spd) ! O
        
    implicit none
    
    ! Out
    real(r8),intent(out)  :: bin_spd
    
    ! Input
    real(r8), intent(inout):: bin_srt
    integer,  intent(in)   :: bin_nbr
    integer,  intent(in)   :: idx
    real(r8), intent(in)   :: wbl_alpha
    real(r8), intent(in)   :: wbl_scl_prm
    real(r8), intent(in)   :: wbl_shp_prm
    real(r8), intent(in)   :: bin_wgt
    
    ! static variables
    real(r8)              :: min_bin_spd
    real(r8)              :: max_bin_spd
    real(r8)              :: max_inc_gmm_p
    real(r8)              :: min_inc_gmm_q
    real(r8)              :: mmt_0           ! Weibull moments
    real(r8)              :: mmt_1
    real(r8)              :: prb_u_lss
    real(r8)              :: xt_max
    real(r8)              :: xt_min

    if (idx == 1) then
       prb_u_lss = 0.0_r8
       min_inc_gmm_q = 1.0_r8
    else
       prb_u_lss = bin_wgt*(real(idx-1,r8))
       min_inc_gmm_q = bin_srt
    end if

    if (idx == bin_nbr) then
       min_bin_spd = wbl_scl_prm*(-1.0_r8*log(1.0_r8 - (prb_u_lss)))**(1.0_r8/wbl_shp_prm)
       xt_min = (min_bin_spd/wbl_scl_prm)**wbl_shp_prm
       mmt_0 = exp(-1.0_r8*xt_min)
       mmt_1 = wbl_scl_prm*gamma(wbl_alpha)*min_inc_gmm_q
       bin_spd = mmt_1/mmt_0
    else
       min_bin_spd = wbl_scl_prm*(-1.0_r8*log(1.0_r8 - (prb_u_lss)))**(1.0_r8/wbl_shp_prm)
       max_bin_spd = wbl_scl_prm*(-1.0_r8*log(1.0_r8 - (prb_u_lss+bin_wgt)))**(1.0_r8/wbl_shp_prm)
       xt_min = (min_bin_spd/wbl_scl_prm)**wbl_shp_prm
       xt_max = (max_bin_spd/wbl_scl_prm)**wbl_shp_prm
       call gratio(wbl_alpha,xt_max,max_inc_gmm_p,bin_srt,0)
       mmt_0 = exp(-1.0_r8*xt_min)-exp(-1.0_r8*xt_max)
       mmt_1 = wbl_scl_prm*gamma(wbl_alpha)*(min_inc_gmm_q - bin_srt)
       bin_spd = mmt_1/mmt_0
    end if  
    
    return
  end subroutine get_wbn_spd
end module wnd_PDF_mdl
