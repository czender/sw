      COMPLEX FUNCTION RAYWAT( WL, TEMP)

c       Calculates complex refractive index of pure liquid water for
c       wavelengths between 2.0 microns and 10.0 m.

c    I N P U T :  WL   = wavelength ( 2 microns to 1.e7 microns )
c                 TEMP = temperature (K)

c    O U T P U T :  RAYWAT = complex refractive index
c                            ( with positive imag. part )

c    METHOD : Ray's analytic fits based on some theories of Debye

c    REFERENCE : Ray, P., 1972:  Broadband Complex Refractive Indices
c                  of Ice and Water, Appl. Opt. 11, 1836-1844

c ---------------------------------------------------------------------

c            *** Specifications of local variables

c      EPSILS :  static dielectric constant
c                   ( epsilon-sub-s, Ray Eq. 4 )
c      EPSINF :  high-frequency dielectric const
c                   ( epsilon-sub-s, Ray Eq. 7a )
c      EPSRE, :  real and imaginary parts of dielectric constant
c        EPSIM      ( Ray Eqs. 5,6 )
c      PI     :  3.14159...
c      AB     :  summation terms, Ray Eq. 8
c      ARRE   :  correction to real index for WL < 6 micron
c      BRRE   :  correction to imag. index
c      RRE    :  real part of refractive index
c      RIM    :  imag part of refractive index
c      TC     :  Centigrade temperature
c      TERM   :  summation terms, Ray Eq. 9
c      TBARP1 :  TBAR + 1, coefficient in Ray Eq. 9
c      WLCM   :  wavelength in cm
c      WLDEBY :  wavelength above which to apply Debye theory
c      WLMIN  :  minimum tabulated wavelength
c      WLMAX  :  maximum tabulated wavelength
c      WLS    :  relaxation wavelength in cm
c                   ( lambda-sub-s, Ray Eq. 7C )
c      ALPHA  :  temporary storage variables
c ---------------------------------------------------------------------

c     .. Scalar Arguments ..

      REAL      TEMP, WL
c     ..
c     .. Local Scalars ..

      REAL      ALPHA, BET, DEL, EPSILS, EPSINF, FREC, GAM, PI, RIM,
     &          RRE, TC, WLCEN, WLCM, WLMAX, WLMIN, WLS, WVL
      COMPLEX   CEPS, CI
c     ..
c     .. External Subroutines ..

      EXTERNAL  ERRMSG
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ABS, AIMAG, LOG10, ASIN, CMPLX, COS, EXP,
     &          REAL, SIN, SQRT
c     ..
c     .. Statement Functions ..

      REAL      AB, ARRE, BRRE, TBARP1, TERM
c     ..
c     .. Save statement ..

      SAVE      PI
c     ..
c     .. Data statements ..

      DATA      PI / 0.0 /, CI / ( 0.0, 1.0) /
      DATA      WLMIN, WLMAX / 1.0, 1.E+8 /
c     ..
c     .. Statement Function definitions ..

      AB( WVL, BET, WLCEN, DEL, GAM) = BET*
     &   EXP( - ABS(LOG10(WVL/WLCEN)/DEL)**GAM)

      TBARP1( TC, WVL) = 1.+ 1.E-4*( TC - 25.)*EXP( (WVL/4.)**0.25)

      TERM( WVL, FREC, BET, GAM) = 
     &    BET*( FREC**2 - (1.E4/WVL)**2) /
     &       ( (FREC**2 - (1.E4/WVL)**2)**2 + GAM*(1.E4/WVL)**2)

      BRRE( TC, WVL) = TBARP1( TC, WVL) *
     &       SQRT( 1.83899 + TERM(WVL, 1639.0,  52340.4,  10399.2)
     &                     + TERM(WVL, 588.24, 345005.0, 259913.0)
     &                     + TERM(WVL, 161.29,  43319.7,  27661.2) )

      ARRE( TC, WVL) = TBARP1( TC, WVL)*
     &       SQRT( 1.79907 + TERM(WVL, 3352.27, 999140.0,  151963.0) 
     &                     + TERM(WVL, 1639.00,  50483.5,    9246.27) 
     &                     + TERM(WVL,  588.24, 844697.0, 1076150.0) )
c     ..


      IF( WL.LT.WLMIN .OR. WL.GT.WLMAX )
     &    CALL ERRMSG( 'RAYWAT--wavelength out of table range',.TRUE.)

      IF( PI.EQ.0.0 ) PI = 2.*ASIN( 1.0)

c                       *** use Ray's fits to Debye theory expressions

      TC     = TEMP - 273.
      EPSINF = 5.27137 + TC*( 0.0216474 - 0.00131198*TC)
      EPSILS = 78.54*( 1. + (TC - 25.)*
     &         (- 4.579E-3 + (TC - 25.)*(1.19E-5 - 2.8E-8*(TC -
     &         25.))))

c                                               *** Ray Eq. 7B
      ALPHA  = - 16.8129 / TEMP + 0.0609265
      WLCM   = 1.E-4*WL
      WLS    = 3.3836E-4*EXP( 2514. / TEMP)
      CEPS   = EPSINF + ( EPSILS - EPSINF) /
     &         ( 1.0 + (COS(0.5*PI*(1. - ALPHA)) + CI*SIN(0.5*PI*(1. -
     &         ALPHA)))*(WLS/WLCM)**(1. - ALPHA)) - 0.00666667*CI*WLCM

c                             ** complex refractive index from Cole-
c                             ** Cole extension to Debye theory

      RAYWAT = SQRT( CEPS)
      RIM    = - AIMAG( RAYWAT)
      RRE    = REAL( RAYWAT)
c                                  ** corrections to imag. index to
c                                  ** account for absorption bands
c                                  ** ( Ray Eq. 8 + Table 2 )

      IF( WL.LT.3000. .AND. WL.GT.300.) THEN

         RIM    = RIM + AB( WL, 0.25, 300., 0.47, 3.0) +
     &            AB( WL, 0.41, 62., 0.35, 1.7) +
     &            AB( WL, 0.39, 17., 0.45, 1.3)

      ELSE IF( WL.LE.300. .AND. WL.GT.62.) THEN

         RIM    = RIM + AB( WL, 0.25, 300., 0.40, 2.0) +
     &            AB( WL, 0.41, 62., 0.35, 1.7) +
     &            AB( WL, 0.39, 17., 0.45, 1.3)

      ELSE IF( WL.LE.62. .AND. WL.GT.17.) THEN

         RIM    = RIM + AB( WL, 0.25, 300., 0.40, 2.0) +
     &            AB( WL, 0.41, 62., 0.22, 1.8) +
     &            AB( WL, 0.39, 17., 0.45, 1.3)

      ELSE IF( WL.LE.17. .AND. WL.GT.6.1) THEN

         RIM    = RIM + AB( WL, 0.12, 6.1, 0.042, 0.6) +
     &            AB( WL, 0.39, 17., 0.165, 2.4) +
     &            AB( WL, 0.41, 62., 0.22, 1.8)

      ELSE IF( WL.LE.6.1 .AND. WL.GT.4.95) THEN

         RIM    = RIM + AB( WL, 0.12, 6.1, 0.009, 2.0) +
     &            AB( WL, 0.01, 4.95, 0.05, 1.0)

      ELSE IF( WL.LE.4.95 .AND. WL.GT.2.95) THEN

         RIM    = RIM + AB( WL, 0.27, 2.97, 0.04, 2.0) +
     &            AB( WL, 0.01, 4.95, 0.06, 1.0)

      ELSE IF( WL.LE.2.95) THEN

         RIM    = RIM + AB( WL, 0.27, 2.97, 0.025, 2.0) +
     &            AB( WL, 0.01, 4.95, 0.06, 1.0)

      END IF
c                                  ** corrections to real index
c                                  ** ( Ray Eq. 9 + Table 3 )

      IF( WL.LT.1000. .AND. WL.GT.340) THEN

         RRE    = RRE*( (WL - 340.) / 660.) +
     &            BRRE( TC, WL)*( (1000. - WL) / 660.)

      ELSE IF( WL.LE.340. .AND. WL.GT.7.0) THEN

         RRE    = BRRE( TC, WL)

      ELSE IF( WL.LE.7. .AND. WL.GT.6.) THEN

         RRE    = ARRE( TC, WL)*( 7. - WL) + BRRE( TC, WL)*( WL - 6.)

      ELSE IF( WL.LE.6.0) THEN

         RRE    = ARRE( TC, WL)

      END IF

      RAYWAT = CMPLX( RRE, RIM)

      END
