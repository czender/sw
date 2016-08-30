      REAL FUNCTION BDREF( WVNMLO, WVNMHI, MU, MUP, DPHI )
c      Supplies surface bi-directional reflectivity.
c      NOTE 1: Bidirectional reflectivity in DISORT is defined
c              by Eq. 39 in STWL.
c      NOTE 2: Both MU and MU0 (cosines of reflection and incidence
c              angles) are positive.
c  INPUT:
c    WVNMLO : Lower wavenumber (inv cm) of spectral interval
c    WVNMHI : Upper wavenumber (inv cm) of spectral interval
c    MU     : Cosine of angle of reflection (positive)
c    MUP    : Cosine of angle of incidence (positive)
c    DPHI   : Difference of azimuth angles of incidence and reflection
c                (radians)
c   Called by- DREF, SURFAC
c +-------------------------------------------------------------------+
c     .. Scalar Arguments ..
      REAL      DPHI, MU, MUP, WVNMHI, WVNMLO
      common /brdf_com/ brdf_typ
      integer brdf_typ
      common /brdf_EOS/ f_iso, f_vol, f_geo
c++   csz 20041109 fix g95 warning
      real brdf_CoxMunk
      real brdf_Hapke
      real brdf_Minnaert
      real brdf_RahmanPinty
      real brdf_LommelSeeliger
      real brdf_LiSparse
      real brdf_RossThick
c--   csz
      real f_iso
      real f_vol
      real f_geo
c     local variables
      real pi
      parameter(pi = 3.141592741)
      if (brdf_typ .eq. 0) then
         bdref = f_iso * pi
     &        + f_vol * brdf_RossThick(mu, mup, dphi)
     &        + f_geo * brdf_LiSparse(mu, mup, dphi)
      elseif (brdf_typ .eq. 1) then
         bdref = brdf_CoxMunk(mu, mup, dphi)
      elseif (brdf_typ .eq. 2) then
         bdref = brdf_Hapke(mu, mup, dphi)
      elseif (brdf_typ .eq. 3) then
         bdref = brdf_Minnaert(mu, mup)
      elseif (brdf_typ .eq. 4) then
         bdref = brdf_RahmanPinty(mu, mup, dphi)
      elseif (brdf_typ .eq. 5) then
         bdref = brdf_LommelSeeliger(mu, mup)
      elseif (brdf_typ .eq. 6) then
         bdref = brdf_Gar86(mu)
      elseif (brdf_typ .eq. 7) then
         bdref = brdf_CiF12(mu)
      else
         stop 'Invalid brdf_typ in BDREF()'
      endif
c++   csz 20070217 fix gfortran warnings
      wvnmlo=wvnmlo
      wvnmhi=wvnmhi
c--   csz
c     print*, 'BDREF: ', bdref, mu, mup, dphi
      RETURN
      END
      
      real function brdf_CoxMunk(mu, mup, dphi)
      implicit none
c     Cox and Munk BDRF
c     input variables
      real mu, mup, dphi
c     Cox and Munk parameters
      common /brdf_CoxMunk_prm/ nrm_cff_CM, wnd_spd, idx_rfr_srf
      real nrm_cff_CM           ! empirical normalization constant
      real wnd_spd              ! suface windspeed m / s
      real idx_rfr_srf          ! surface index of refraction
c     local variables
      real pi
      parameter(pi = 3.141592741)
      real theta                ! scattering angle between incidence and reflection
      real costheta             ! cos scattering angle
      real omega                ! angle between incidence/reflection and facet normal
      real cosomega             ! cos omega
      real beta                 ! angle between sloped facet normal and vertical
      real cosbeta              ! cos beta
      real r_Fresnel            ! Fresnel reflection coefficient
      real srf_slp_CoxMunk      ! facet slope distribution
      real sigma_sq             ! slope variance
c++   csz 200070217
c     calculate brdf
      costheta = - mu * mup +
     &     (1.0 - mu**2)**0.5 * (1.0 - mup**2)**0.5 * cos(dphi)
      theta = acos(costheta)
      omega = 0.5 * (pi - theta)
      cosomega = cos(omega)
      cosbeta = 0.5 * (mu + mup) / cosomega
      beta = acos(cosbeta)
      sigma_sq = 0.003 + 0.00512 * wnd_spd
      r_Fresnel = 0.5 *
     &     (cosomega - sqrt(idx_rfr_srf**2 - 1 + cosomega**2))**2
     &     / (cosomega + sqrt(idx_rfr_srf**2 - 1 + cosomega**2))**2
     &     + 0.5 *
     &     (- idx_rfr_srf**2 * cosomega
     &     + sqrt(idx_rfr_srf**2 - 1 + cosomega**2))**2
     &     / (idx_rfr_srf**2 * cosomega
     &     + sqrt(idx_rfr_srf**2 - 1 + cosomega**2))**2
      srf_slp_CoxMunk = 1.0 / (pi * sigma_sq)
     &     * exp(- (tan(beta))**2 / sigma_sq)
      brdf_CoxMunk = srf_slp_CoxMunk * r_Fresnel
     &     * nrm_cff_CM / (mup * mu * cosbeta**4)
      return
      end
      
      real function brdf_Hapke(mu, mup, dphi)
      implicit none
c     Hapke BDRF
c     ** Hapke's BRDF model (times Pi/Mu0)
c     ** (Hapke, B., Theory of reflectance
c     ** and emittance spectroscopy, Cambridge
c     ** University Press, 1993, Eq. 8.89 on
c     ** page 233. Parameters are from
c     ** Fig. 8.15 on page 231, expect for w.)
c     input variables
      real mu, mup, dphi
c     Hapke parameters
      common /brdf_Hapke_prm/ b0, hh, w
      real b0                   ! empirical factor for the finite size of particles
      real hh                   ! angular width parameter of opposition effect
      real w                    ! single scattering albedo
c     local variables
      real theta                ! phase angle (pi - scattering angle)!
      real costheta             ! cos phase angle
      real b                    ! term that accounts for the opposition effect
      real gamma                ! albedo factor
      real P                    ! Hapke scattering phase function
      real H0                   ! H(mu0)
      real H                    ! H(mu)
c     calculate brdf
      costheta = mu * mup +
     &     (1.0 - mu**2)**0.5 * (1.0 - mup**2)**0.5 * cos(dphi)
      theta = acos(costheta)
      P = 1.0 + 0.5 * costheta
      b = b0 * hh / (hh + tan(theta / 2.0))
      gamma = sqrt(1.0 - w)
      H0 = (1.0 + 2.0 * mup) / (1.0 + 2.0 * mup * gamma)
      H = (1.0 + 2.0 * mu) / (1.0 + 2.0 * mu * gamma)
      brdf_Hapke = w / 4.0 / (mu + mup)
     &     * ((1.0 + b) * P + H0 * H - 1.0)
      return
      end
      
      real function brdf_Minnaert(mu, mup)
      implicit none
c     Minnaert BDRF
c     input variables
      real mu, mup
c     Minnaert parameters
      common /brdf_Minnaert_prm/ nrm_rfl_M, k_cff_M
      real nrm_rfl_M            ! normal reflectance
      real k_cff_M              ! coefficient for power dependence on mu, mup
      brdf_Minnaert = nrm_rfl_M
     &     * mu**(k_cff_M - 1.0) * mup**(k_cff_M - 1.0)
      return
      end
      
      real function brdf_RahmanPinty(mu, mup, dphi)
      implicit none
c     Rahman, Pinty and Verstraete BDRF
c     input variables
      real mu, mup, dphi
c     RahmanPinty parameters
      common /brdf_RahmanPinty_prm/ nrm_cff_RP, k_cff_RP, g_phs
      real nrm_cff_RP           ! normalization coefficient
      real k_cff_RP             ! coefficient for power dependence on mu, mup
      real g_phs                ! asymmetry parameter in Henyey Greenstein function
c     local variables
      real pi
      parameter(pi = 3.141592741)
      real costheta             ! cos scattering angle
      real tngnt                ! tangent of outgoing angle
      real tngntp               ! tangent of incident angle
      real P                    ! Henyey Greenstein phase function
      real G
      costheta = - mu * mup +
     &     (1.0 - mu**2)**0.5 * (1.0 - mup**2)**0.5 * cos(dphi)
      P = (1.0 - g_phs**2)
     &     / (1.0 + g_phs**2 - 2.0 * g_phs * costheta)**1.5
      tngnt = (1.0 - mu**2)**0.5 / mu
      tngntp = (1.0 - mup**2)**0.5 / mup
      G = (tngnt**2 + tngntp**2 - 2.0 * tngnt * tngntp * cos(dphi))**0.5
      brdf_RahmanPinty
     &     = pi * nrm_cff_RP * mu**(k_cff_RP - 1.0) * mup**(k_cff_RP - 1.0)
     &     / (mu + mup)**(k_cff_RP - 1.0)
     &     * P * (1.0 + (1.0 - nrm_cff_RP) / (1.0 + G))
      return
      end
      
      real function brdf_LommelSeeliger(mu, mup)
      implicit none
c     Lommel and Seeliger BDRF
c     input variables
      real mu, mup
c     LommelSeeliger parameters
      common /brdf_LommelSeeliger_prm/ nrm_rfl_LS
      real nrm_rfl_LS           ! normal reflectance
      brdf_LommelSeeliger = 2.0 * nrm_rfl_LS / (mu + mup)
      return
      end
      
      real function brdf_RossThick(mu, mup, dphi)
      implicit none
c     Ross-Thick BDRF
c     ** EOS MODIS BRDF/Albedo Product
c     ** Algorithm Theoretical Basis Document Version 5.0
c     ** Strahler A.H., Muller J.P., et al., April 1999
c     ** eq. 38 note: extra factor of pi due to different
c     BDRF normalization conventions
c     input variables
      real mu, mup, dphi
c     local variables
      real pi
      parameter(pi = 3.141592741)
      real xi                   ! phase angle
      real cosxi                ! cos phase angle
      real sinxi                ! sin phase angle
      cosxi = mu * mup +
     &     (1.0 - mu**2)**0.5 * (1.0 - mup**2)**0.5 * cos(dphi)
      xi = acos(cosxi)
      sinxi = sin(xi)
      brdf_RossThick = (pi * (pi / 2.0 - xi) * cosxi + pi * sinxi)
     &     / (mu + mup) - pi * pi / 4.0
      return
      end
      
      real function brdf_LiSparse(mu, mup, dphi)
      implicit none
c     LiSparse BDRF
c     ** EOS MODIS BRDF/Albedo Product
c     ** Algorithm Theoretical Basis Document Version 5.0
c     ** Strahler A.H., Muller J.P., et al., April 1999
c     ** eq. 39 - 44 note: extra factor of pi due to different
c     BDRF normalization conventions
c     input variables
      real mu, mup, dphi
      common /brdf_LiSparse_prm/ hb, br
      real hb                   ! height to vertical radius ratio
      real br                   ! vertical radius to horizontal radius ratio
c     local variables
      real pi
      parameter(pi = 3.141592741)
c     various intermediate results, in order of appearance
c     see equations cited above
      real tn, tnp, cs, csp, cosxi, dsq, secplsecp, cost, t, opi
      tn = br * tan(acos(mu))
      tnp = br * tan(acos(mup))
      cs = 1.0 / sqrt(1.0 + tn**2)
      csp = 1.0 / sqrt(1.0 + tnp**2)
      cosxi = cs * csp +
     &     (1.0 - cs**2)**0.5 * (1.0 - csp**2)**0.5 * cos(dphi)
      dsq = tn**2 + tnp**2 - 2.0 * tn * tnp * cos(dphi)
      secplsecp = 1.0 / cs + 1.0 / csp
      cost = hb * sqrt(dsq + (tn * tnp * sin(dphi))**2) / secplsecp
      if (cost .gt. 1.0) then
         cost = 1.0
      endif
      if (cost .lt. -1.0) then
         cost = -1.0
      endif
      t = acos(cost)
      opi = (t - sin(t) * cost) * secplsecp
      brdf_LiSparse = opi - pi * secplsecp
     &     + 0.5 * pi * (1.0 + cosxi) / (cs * csp)
      return
      end
      
      real function brdf_Gar86(mu)
      implicit none
c     Gar86 light pollution emission
c     input variables
      real mu
      common /brdf_Gar86_prm/ fff, ggg
      real fff
      real ggg
c     local variables
      real pi
      real theta
      real theta_sqr
      parameter(pi = 3.141592741)
      theta=acos(mu)
      theta_sqr=theta*theta
      brdf_Gar86 = 2.0*ggg*(1.0-fff)*mu+0.554*fff*theta_sqr*theta_sqr
      return
      end                       ! Gar86
      
      real function brdf_CiF12(mu)
      implicit none
c     CiF12 p. 3349 (61) light pollution emission
c     input variables
      real mu
c     CiF12 parameters
      common /brdf_CiF12_prm/ u1, u2, u3
      real u1
      real u2
      real u3
c     local variables
      real q
      real pi
      real theta
      real theta_sqr
      parameter(pi = 3.141592741)
      theta=acos(mu)
      theta_sqr=theta*theta
      if(theta < pi/6.0) then
         q=0.0
      else 
         q=1.0
      endif
!     Received from Fabio Falchi 20160827
!     thank you for your appreciation of our work.
!     Yes, q=1 for theta>30°
!     Here is the upward function obtained from the calibration:
!     http://www.wolframalpha.com/input/?i=polarplot+%5B+%7B(2.*190+Cos%5Bx%5D+%2B+52*0.5543155+x%5E4+%2B+1.777777*7.6*((sqrt((x+-+%5C%5BPi%5D%2F6.)*(x+-+%5C%5BPi%5D%2F6.))%2B(x+-+%5C%5BPi%5D%2F6.))%2F(2.*(x+-+%5C%5BPi%5D%2F6.)))*+Cos%5B3.+(x+-+%5C%5BPi%5D%2F3.)%5D)%2F(2+%5C%5BPi%5D)%7D,++%7Bx,+0,+%5BPi%5D%2F2.%7D%5D
!     I multiplied the parameters (Wa = 1.9 × 10−3, Wb = 5.2 × 10−4, Wc = 7.6 × 10−5) by 10^5 in the equation, but the function's shape is the correct one.
!     Best regards,
!     Fabio
      brdf_CiF12=2.0*mu+0.5543*u2*theta_sqr*theta_sqr+1.778*q*u3*cos(3.0*theta-pi)/(2.0*pi*(1.0+u2+u3))
      return
      end ! CiF12
      
