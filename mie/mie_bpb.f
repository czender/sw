      PROGRAM MIEBPB
c
c logical for test
c
      logical modp
      data modp / .false.    /
C
      parameter ( maxmu = 119 ,maxmup2 = maxmu+2)
      real gauswt(maxmup2),gauspt(maxmup2)
      real xzero(maxmup2),pzero(maxmup2)
      COMMON /RES/ MXANG,MXANG2
C
      PARAMETER (NMBRAD=1000, NMBWAV=100, numbre=1)
      DIMENSION P(MAXMUP2),SCTAN(MAXMUP2),PEFF(MAXMUP2,nmbwav)
      real qeex(nmbwav),qesc(nmbwav),weff(nmbwav)
      real geff(nmbwav),feff(nmbwav)
      real sigsc(nmbwav),sigex(nmbwav),renorm(nmbwav)
      real qe(nmbwav,nmbrad),qs(nmbwav,nmbrad)
      real we(nmbwav,nmbrad),ge(nmbwav,nmbrad)
c
      real ger(nmbrad),fer(nmbrad)
      real gr(nmbwav),fr(nmbwav)
c
      real fe(nmbwav,nmbrad)
      real sge(nmbwav,nmbrad),sgs(nmbwav,nmbrad)
      DIMENSION RAD(NMBRAD),WAVEL(NMBWAV),re(numbre)
      real s(nmbwav),refre(nmbwav),refim(nmbwav)
      DIMENSION QSCAT(NMBRAD),QEXT(NMBRAD),SIZEP(NMBRAD)
      REAL      NMBPRT(NMBRAD)
C
      logical sulfate,water,ice,gamma
c
      data sulfate / .true.  /
      data water   / .false. /
      data ice     / .false. /
      data gamma   / .false. /
c
      LOGICAL PRNT,PLOT
      REAL N0,NMIN,NMAX,MUMIN,MUMAX,MUMN
      DATA    PRNT / .FALSE./
      DATA    PLOT / .true. /
      INTEGER FRQPLT,FRQPRT
      DATA    FRQPLT /  1000 /
      DATA    FRQPRT /  1000 /
c
c      data re /1., 4., 6.,8.,10., 12., 14., 16., 18., 20., 25.,
c     $           30., 35., 40. , 45., 50. /
c
       data re /.1/
c
c ..... gamma distribution
c
c      RNMDEN(R,req,B) = (R**((1.-3*B)/B))*EXP(-R/(B*req))
c
c
c ..... log normal distribution
c
      rnmden(r,req,sig) = (1./r)*
     $                    exp(-.5*((alog(r/req)/alog(sig))**2))
c
c--------------------------------------------------------------------
c..... setup gks:
c
      call setgks
c
      write(6,1234) nmbrad
 1234 format(/' start mie computation '/
     $        ' number of radii in size distribtuion = ',i4)
C
C...... SET DIMENSIONS:
C
      call gauss(gauspt,gauswt,xzero,pzero,maxmu)
      MXAO2  = ((maxmu+2)/2) + 1
      MXANG2 =  maxmu + 2
C
      SIGRAD    =  50.
      N0        =  1.
      BVAL      =  1./6.
      sig       =  2.0
C
      if( sulfate ) then
       nwav     =    1
       wavel(1) = .250
       refre(1) = 1.450
       refim(1) = 2.e-8
       write(6,9876) wavel(1),refre(1),refim(1)
 9876  format(' sulfate aerosol for wavelength = ',f8.3,' microns'/
     $        ' real index of refraction      = ',1pe11.4/
     $        ' imaginary index of refraction = ',1pe11.4)
      endif
c
      if( water ) then
       call optwat(nwav,wavel,refre,refim)
      endif
c
      if( ice   ) then
       call optice(nwav,wavel,refre,refim)
      endif       
c
c.......... loop over equivalent radius:
c
      do 2 m=1,numbre
c      
       if( gamma ) then
         reval  = re(m)
         remax  = (1.-3*bval)*reval
         rad0   = .001*remax
         sigrad = 50.
       else
         reval  = re(m)
         rad0   = re(m) / 100.
         remax  = 30. * re(m)
         sigrad = (real(nmbrad)-n0)/alog10(remax/rad0)
       endif
C
       SUMDEN    =  0.0
       PIE       = 4. * ATAN(1.)
C
       if( nmbrad .gt. 1 ) then
        DO 10 N=1,NMBRAD
         RAD(N) = RAD0 * (10.**((REAL(N)-N0)/SIGRAD))
         RADIUS = RAD(N)
         if( gamma ) then
           NMBPRT(N) = RNMDEN(RADIUS,REVAL,BVAL)
         else
           nmbprt(n) = rnmden(radius,reval,sig)
         endif
   10   CONTINUE
       else
        rad(1)    = re(m)
        nmbprt(1) = 1.
       endif
C
       DO 15 N=2,NMBRAD
C
        RMIN = RAD(N-1)
        SMIN = NMBPRT(N-1)
C
        RMAX = RAD(N)
        SMAX = NMBPRT(N)
C
        SUMDEN    = SUMDEN + 0.5*(SMIN + SMAX) *( RMAX-RMIN )
C
   15  CONTINUE
       if (nmbrad .eq. 1) sumden = 1.0
C
       PRINT 17,SUMDEN
   17  FORMAT(/' ... TRPZD QUADRATURE OF NUMBER DEN = ',1PE11.4)
C
       write(6,19) 
   19  format(/'    size distribution  '/
     $         '  n     radius      dn/dr ')
       DO 20 N=1,NMBRAD
         NMBPRT(N) = NMBPRT(N) / SUMDEN
         write(6,21) n,rad(n),nmbprt(n)
   21    format(2x,i5,2x,2(1pe11.4,2x))
   20  CONTINUE
C
       SUM = 0.0
       DO 30 N=2,NMBRAD
C
        RMIN = RAD(N-1)
        SMIN = NMBPRT(N-1)
C
        RMAX = RAD(N)
        SMAX = NMBPRT(N)
C
        SUM    = SUM + 0.5 * (RMAX*SMAX + RMIN*SMIN) * (RMAX-RMIN)
C
   30  CONTINUE
       RADMN = SUM
C
       SUM1 = 0.0
       SUM2 = 0.0
       SUM3 = 0.0
C
       DO 40 N=2,NMBRAD
C
        RMIN = RAD(N-1) - RADMN
        SMIN = NMBPRT(N-1)
C
        RMAX = RAD(N)   - RADMN
        SMAX = NMBPRT(N)
C
        SUM1 = SUM1 + 0.5*(RMAX*RMAX*SMAX + RMIN*RMIN*SMIN)*(RMAX-RMIN)
C
C
        RMIN = RAD(N-1)
        RMAX = RAD(N)
C
        SUM2 = SUM2 + 0.5*(RMAX*RMAX*SMAX + RMIN*RMIN*SMIN)*(RMAX-RMIN)
        SUM3 = SUM3 + 0.5*(RMAX*RMAX*RMAX*SMAX +
     +                     RMIN*RMIN*RMIN*SMIN)*(RMAX-RMIN)
C
   40  CONTINUE
C
       SIGMN = SQRT(SUM1)
       AMEAN = PIE * SUM2
       RAREA = SQRT(SUM2)
       VMEAN = (4.*PIE/3.) * SUM3
c
       if( nmbrad .eq. 1 ) then
         sigmn = 0.0
         amean = pie*rad(1)**2
         rarea = rad(1)
         radmn = rad(1)
         vmean = (4.*pie/3.)*rad(1)**3
       endif
C
       if( gamma ) then 
c
       REMAX = REVAL * (1.-3.*BVAL)
       PRINT 45,REVAL,BVAL,REMAX,RADMN,SIGMN,RAREA,AMEAN,VMEAN
   45  FORMAT(/'  .... INFO ON SIZE DISTRIBUTION ....'/
     + ' ....   RE        B          RMAX        RMEAN        SIGMA '/
     +        2X,5(1PE11.4,2X)/
     + ' ....                        RAREA        AMEAN        VMEAN '/
     +       28X,4(1PE11.4,2X))
c
       else
c
       PRINT 46,rad0,REVAL,remax,sig,RADMN,RAREA,AMEAN,VMEAN
   46  FORMAT(/'  .... info on size distribution ....'/
     + ' ....   log normal distribution '/
     + ' ....   rmin        re       rmax        sig'/
     +        2X,4(1PE11.4,2X)/
     + ' ....   rmean       rarea    amean       vmean '/
     +        2X,4(1PE11.4,2X))
       endif
C
C..... PERFORM MIE COMPUTATION ....
C
      if( plot ) then
C
        CALL PLTR  (NMBPRT,NMBRAD,RAD,REVAL,BVAL,
     +              '   RADIUS (MICRONS)    $',
     +              ' NUMBER DENSITY / UNIT RADIUS$')
c
      endif
c
      DO 100 NW=1,nwav
C
        PRINT 101,WAVEL(nw),REFRE(nw),REFIM(nw)
  101   FORMAT(/' .... WAVELENGTH = ',1PE11.4,' MICRO-METERS'/
     +   ' .... REAL INDEX OF REFRACTION = ',1PE11.4/
     +   ' .... IMG  INDEX OF REFRACTION = ',1PE11.4/)
C
        do 190 na=1,mxang2
          peff(na,nw) = 0.
  190   continue
c
        DO 200 NR=1,NMBRAD
c
          write(6,201) nr
  201     format('  200 loop nr = ',i5)
C
          RADIUS = RAD(NR)
C
          PLOT = .FALSE.
          IF( MOD(NR,FRQPLT) .EQ. 0 ) PLOT = .TRUE.
          IF(  NR .EQ. 1 .OR. NR .EQ. NMBRAD ) PLOT = .TRUE.
C
          PRNT = .FALSE.
          IF( MOD(NR,FRQPRT) .EQ. 0 ) PRNT = .TRUE.
          IF(  NMBRAD .EQ. 1 ) PRNT = .TRUE.
          IF(  NR .EQ. 1 .OR. NR .EQ. NMBRAD ) prnt = .true.
C
          CALL CALLBH(MXAO2,gauspt,WAVEL(NW),RAD(NR),
     +                REFRE(nw),REFIM(nw),P,SCTAN,
     +                X,QSCA,QEXTN,prnt,plot)
C
          QSCAT(NR) = QSCA
          QEXT (NR) = QEXTN
          SIZEP(NR) = X
c
          if( modp ) then 
            do 1777 na=1,mxang2
              p(na) = 1.
 1777       continue
          endif
C
C....... COMPUTE CONTRIBUTION TO EFFECTIVE
C....... SCATTERING PHASE FUNCTION .......
C
      DO 150 NA=1,MXANG2
        PEFF(NA,nw) = PEFF(NA,nw) +
     +             P(NA)*QSCAT(NR)*PIE*RADIUS*RADIUS*NMBPRT(NR)
  150 CONTINUE
c
c perform angular integration using gaussian quadrature
c
        sumg  = 0.0
        sumg2 = 0.0
        sum   = 0.0
c
        do 170 na=1,mxang2
c
          wtg   =         p(na)*gauspt(na)*gauswt(na)
          sumg  = sumg +  p(na)*gauspt(na)*gauswt(na)
c
          p2    = 0.5 * (3.*gauspt(na)*gauspt(na) - 1.)
c
          wtg2  =         p(na)*p2*gauswt(na)
          sumg2 = sumg2 + p(na)*p2*gauswt(na)
          sum   = sum  +  p(na)*gauswt(na)
c
c          write(6,1770) na,gauspt(na),p(na),p2
c 1770     format(2x,i3,1x,3(1pe16.9,1x))
c
c          write(6,1772) wtg,wtg2
c 1772     format(' wtg,wtg2 = ',2(1pe16.9,1x))
c
c          write(6,1771) sum,sumg,sumg2
c 1771     format(' sum sumg sumg2 = ',3(1pe16.9,1x))
c
  170   continue
C
        ger(nr) =  ABS( SUMG  / SUM )
        fer(nr) =  ABS( SUMG2 / SUM )
C
  200   CONTINUE
C
      if( plot .and. nmbrad .gt. 1 ) then
        CALL PLTQ  (QSCAT,NMBRAD,SIZEP,WAVEL(nw),
     +              REFRE(nw),REFIM(nw),
     +              '    SIZE PARAMETER   $',
     +              ' SCATTERING EFFICIENCY    $')
      endif
C
C..... COMPUTE EFFECTIVE SCATTERING EFFICIENCY, ........
C..... SINGLE PARTICAL SCATTERING ALBEDO, AND ..........
C..... SCATTERING PHASE FUNCTION .....
C
      sumsc= 0.0
      sumex= 0.0
      SUM  = 0.0
      sumg = 0.0
      sumf = 0.0
C
      DO 300 N=1,NMBRAD-1
C
        RMIN  = RAD(N)
        RMAX  = RAD(N+1)
C
        NMIN  = NMBPRT(N)
        NMAX  = NMBPRT(N+1)
C
C
        SCMIN = QSCAT(N)  * PIE * (RMIN*RMIN) * NMIN
        SCMAX = QSCAT(N+1)* PIE * (RMAX*RMAX) * NMAX
C
        EXMIN = QEXT (N)  * PIE * (RMIN*RMIN) * NMIN
        EXMAX = QEXT (N+1)* PIE * (RMAX*RMAX) * NMAX
C
C
        sumsc = sumsc + ((SCMIN+SCMAX)/2.) * (RMAX-RMIN)
C
        sumex = sumex + ((EXMIN+EXMAX)/2.) * (RMAX-RMIN)
c
        sumg  = sumg  + ((SCMIN*ger(n)+SCMAX*ger(n+1))/2.) 
     $                * (RMAX-RMIN)
c
        sumf  = sumf  + ((SCMIN*fer(n)+SCMAX*fer(n+1))/2.) 
     $                * (RMAX-RMIN)
C
        SUM   = SUM   + ((NMIN  + NMAX )/2.) * (RMAX-RMIN)
C
  300 CONTINUE
c
C
      if( nmbrad .gt. 1 ) then
c
        SIGSC(nw) = sumsc / SUM
        qesc(nw)  = (sumsc / AMEAN)
C
        SIGEX(nw) = sumex / SUM
        qeex(nw)  = (sumex/ AMEAN)
c
        weff(nw) = qesc(nw) / qeex(nw)
c
        gr(nw)   =  ABS( SUMG  / sumsc )
        fr(nw)   =  ABS( SUMF  / sumsc )
c
      else
        sigsc(nw) = qscat(1)*pie*rad(1)**2
        sigex(nw) = qext (1)*pie*rad(1)**2
        qesc(nw)  = sigsc(nw)/amean
        qeex(nw)  = sigex(nw)/amean
        weff(nw)  = sigsc(nw)/sigex(nw)
        gr(nw)    = ger(1)
        fr(nw)    = fer(1)
      endif
C
C.... CONVERT CROSS SECTIONS FROM MICRO-METERS TO CM2:
C
      sigsc(nw) = sigsc(nw) * 1.E-8
      sigex(nw) = sigex(nw) * 1.E-8
C
      RENORM(nw) = PEFF(1,nw) / SIGSC(nw)
      PEFF(1,nw) = 1.0
      DO 350 NA=2,MXANG2
        PEFF(NA,nw) = PEFF(NA,nw) / sigsc(nw)
        PEFF(NA,nw) = PEFF(NA,nw) / RENORM(nw)
  350 CONTINUE
c
c perform angular integration using gaussian quadrature
c
      sumg  = 0.0
      sumg2 = 0.0
      sum   = 0.0
c
      do 370 na=1,mxang2
c
        sumg  = sumg  + peff(na,nw)*gauspt(na)*gauswt(na)
c
        p2    = 0.5 * (3.*gauspt(na)*gauspt(na) - 1.)
c
        sumg2 = sumg2 + peff(na,nw)*p2*gauswt(na)
c
        sum   = sum   + peff(na,nw)*gauswt(na)
c
  370 continue
C
      geff(nw) =  ABS( SUMG  / SUM )
      feff(nw) =  ABS( SUMG2 / SUM )
C
      if( gamma ) then
      PRINT 310,BVAL,REVAL,QESC(nw),SIGSC(nw),
     $          QEEX(nw),SIGEX(nw),WEFF(nw),geff(nw),feff(nw),
     $          gr(nw),fr(nw)
  310 FORMAT(/' ... SIZE DISTRIBUTION  B AND RE = ',2(1PE11.4,2X)/
     +        ' ... EFFECTIVE SCATTERING EFFICIENCY    = ',1PE11.4/
     + ' ... SCATTERING CROSS SECTION           = ',1PE11.4,' CM2'/
     +        ' ... EFFECTIVE EXTINCTION EFFICIENCY    = ',1PE11.4/
     + ' ... EXTINCTION CROSS SECTION           = ',1PE11.4,' CM2'/
     +        ' ... EFFECTIVE SINGLE SCATTERING ALBEDO = ',1PE11.4/
     +        ' ...    radius first, then angle '/
     +        ' ... EFFECTIVE assymetry parameter      = ',1PE15.6/
     +        ' ... EFFECTIVE forward scattering frctn = ',1PE15.6/
     +        ' ...    angle first, then radius '/
     +        ' ... EFFECTIVE assymetry parameter      = ',1PE15.6/
     +        ' ... EFFECTIVE forward scattering frctn = ',1PE15.6)
      else
      PRINT 311,sig,REVAL,QESC(nw),SIGSC(nw),
     $          QEEX(nw),SIGEX(nw),WEFF(nw),geff(nw),feff(nw),
     $          gr(nw),fr(nw)
  311 FORMAT
     $(/' ... log normal SIZE DSTRBTN  sig and re= ',2(1PE11.4,2X)/
     +        ' ... EFFECTIVE SCATTERING EFFICIENCY    = ',1PE11.4/
     + ' ... SCATTERING CROSS SECTION           = ',1PE11.4,' CM2'/
     +        ' ... EFFECTIVE EXTINCTION EFFICIENCY    = ',1PE11.4/
     + ' ... EXTINCTION CROSS SECTION           = ',1PE11.4,' CM2'/
     +        ' ... EFFECTIVE SINGLE SCATTERING ALBEDO = ',1PE11.4/
     +        ' ...    radius first, then angle '/
     +        ' ... EFFECTIVE assymetry parameter      = ',1PE15.6/
     +        ' ... EFFECTIVE forward scattering frctn = ',1PE15.6/
     +        ' ...    angle first, then radius '/
     +        ' ... EFFECTIVE assymetry parameter      = ',1PE15.6/
     +        ' ... EFFECTIVE forward scattering frctn = ',1PE15.6)
      endif
c
c print out mean scattering phase function
c
      write(6,349) 1,sctan(1),peff(1,nw)
  349 format('  n    angle      peff  '/
     $       2x,i5,2x,f7.2,1x,1pe11.4)
      DO 352 NA=2,MXANG2
        write(6,351) na,sctan(na),peff(na,nw)
  351   format(2x,i5,2x,f7.2,1x,1pe11.4)
  352 CONTINUE
c
c print out g,f for each radius r
c
      write(6,1349) 1,rad(1),nmbprt(1),ger(1),fer(1)
 1349 format('  n    rad      n(r)     g(r)      f(r)'/
     $       2x,i5,2x,4(1pe11.4,2x))
      DO 1352 nr=2,nmbrad
        write(6,1351) nr,rad(nr),nmbprt(nr),ger(nr),fer(nr)
 1351   format(2x,i5,2x,4(1pe11.4,2x))
 1352 CONTINUE
C
      qe(nw,m) = qeex(nw)
      qs(nw,m) = qesc(nw)
      we(nw,m) = weff(nw)
      ge(nw,m) = geff(nw)
c
      sgs(nw,m) = sigsc(nw)
      sge(nw,m) = sigex(nw)
c
  100 continue
C
      if( plot ) then
       if( gamma ) then
        CALL PLTQEF(PEFF,MXANG2,SCTAN,WAVEL,REVAL,BVAL,
     +              '  SCATTERING ANGLE     $',
     +              ' MEAN SCAT PHASE FNCTN  gamma $')
       else
        CALL PLTQEF(PEFF,MXANG2,SCTAN,WAVEL,REVAL,sig,
     +              '  SCATTERING ANGLE     $',
     +              ' MEAN SCAT PHASE FNCTN  log normal $')
       endif
      endif
C
    2 continue
C
      write(6,1235) 
 1235 format(/' .... final list of mie properties ...'/
     $ ' .... re      qe       qs      we       ge',
c++csz
c     $ '    sg ext     sg sct')
     $ '    sg sct     sg ext')
c--csz
      do 800 nw=1,nwav
       write(6,801) wavel(nw)
  801  format(14x,f8.3,' micro-meters')
       do 900 nr=1,numbre
        write(6,901) re(nr),qe(nw,nr),qs(nw,nr),
     $               we(nw,nr),ge(nw,nr),sgs(nw,nr),
     $               sge(nw,nr)
  901   format(3x,f5.1,1x,4(f8.6,1x),2(1pe11.4,1x))
  900  continue
  800 continue
c
      PRINT 1000
 1000 FORMAT(/' .... END OF MIE PROGRAM ....'/)
c
      call gclwk(9)
      call clsgks
c
      END
C-------------------------------------------------------------------
      SUBROUTINE CALLBH(NANG,amu,WAVEL,RAD,
     +                  REFRE,REFIM,P,SCTAN,
     +                  X,QSCA,QEXT,PRNT,PLOT)
C       ---------------------------------------------------------------
C       CALLBH CALCULATES THE SIZE PARAMETER (X) AND RELATIVE
C       REFRACTIVE INDEX (REFREL) FOR A GIVEN SPHERE REFRACTIVE
C       INDEX, MEDIUM REFRACTIVE INDEX, RADIUS, AND FREE SPACE
C       WAVELENGTH.  IT THEN CALLS BHMIE, THE SUBROUTINE THAT COMPUTES
C       AMPLITUDE SCATTERING MATRIX ELEMENTS AND EFFICIENCIES
C       ---------------------------------------------------------------
      COMPLEX REFREL,S1(10000),S2(10000)
      REAL amu(1),P(1),SCTAN(1)
      LOGICAL PRNT,PLOT
      IF( PRNT ) WRITE (6,11)
C       ------------------------------------------------------
C       REFMED = (REAL) REFRACTIVE INDEX OF SURROUNDING MEDIUM
C       ------------------------------------------------------
      REFMED = 1.0
C       --------------------------------------------
C       REFRACTIVE INDEX OF SPHERE = REFRE + I*REFIM
C       --------------------------------------------
      REFREL = CMPLX(REFRE,REFIM)/REFMED
      IF( PRNT ) WRITE(6,12) REFMED,REFRE,REFIM
C       ----------------------------------------------
C       RADIUS (RAD) AND WAVELENGTH (WAVEL) SAME UNITS
C       ----------------------------------------------
      X     = 2. * 3.14159265*RAD*REFMED/WAVEL
      IF( PRNT ) WRITE(6,13) RAD,WAVEL
      IF( PRNT ) WRITE(6,14) X
C       ------------------------------------------------
      CALL BHMIE(X,REFREL,NANG,amu,S1,S2,QEXT,QSCA,QBACK,PRNT)
      IF( PRNT ) WRITE (6,65) QSCA,QEXT,QBACK
      IF( PRNT ) WRITE (6,17)
C       --------------------------------------------------
C       S33 AND S34 MATRIX ELEMNTS NORMALIZED BY S11.
C       S11 IS NORMALIZED TO 1.0 IN THE FORMWARD DIRECTION
C       POL=DEGREE OF POLARIZATION (INCIDENT UNPOLARIZED LIGHT)
C       --------------------------------------------------
      S11NOR = 0.5*(CABS(S2(1))**2+CABS(S1(1))**2)
      NAN    = 2*NANG - 1
      pie    = 4. * atan(1.)
      DO 355 J=1,NAN
        AJ = J
        S11=0.5*CABS(S2(J))*CABS(S2(J))
        S11=S11+0.5*CABS(S1(J))*CABS(S1(J))
        S12=0.5*CABS(S2(J))*CABS(S2(J))
        S12=S12-0.5*CABS(S1(J))*CABS(S1(J))
        POL=-S12/S11
        S33=REAL(S2(J)*CONJG(S1(J)))
        S33=S33/S11
        S34=AIMAG(S2(J)*CONJG(S1(J)))
        S34=S34/S11
        S11=S11/S11NOR
        ANG = acos(amu(j))*(180./pie)
        SCTAN(J) = ANG
C
        IF( PRNT ) WRITE(6,75) ANG,S11,POL,S33,S34
        P(J)      = S11
  355 CONTINUE
      IF( PLOT ) CALL PLTPSC(P,NAN,SCTAN,X,WAVEL,QSCA,QEXT,
     +           '    SCATTERING ANGLE $',
     +           '    RELATIVE SCATTERING   $')
C
   65 FORMAT(//,1X,' QSCA =',E13.6,3X,' QEXT = ',E13.6,3X,
     + ' QBACK = ',E13.6)
   75 FORMAT(1X,F6.2,2X,E13.6,2X,E13.6,2X,E13.6,2X,E13.6)
   11 FORMAT(/ ' SPHERE SCATTERING PROGRAM'//)
   12 FORMAT(5X,' REFMED = ',F8.4,3X,' REFRE = ',E14.6,3X,
     + ' REFIM = ',E14.6)
   13 FORMAT(5X,' SPHERE RADIUS = ',F7.3,3X,' WAVELENGTH = ',F7.4)
   14 FORMAT(5X,' SIZE PARAMETER = ',F8.3/)
   17 FORMAT(//,2X,'ANGLE',7X,'S11',13X,'POL',13X,'S33',13X,'S34'//)
      END
C       -------------------------------------------------------
C       SUBROUTINE BHMIE CALCULATES AMPLITUDE SCATTERING MATRIX
C       ELEMENTS AND EFFICIENCIES FOR EXTINCTION, TOTAL SCATTERING
C       AND BACKSCATTERING FOR A GIVEN SIZE PARAMETER AND
C       RELATIVE REFRACTIVE INDEX
C       -------------------------------------------------------
      SUBROUTINE BHMIE (X,REFREL,NANG,amu,S1,S2,
     +                  QEXT,QSCA,QBACK,PRNT)
      PARAMETER (MXA=10000,MXA2=20000)
      DIMENSION AMU(1),THETA(MXA),PI(MXA),TAU(MXA),PI0(MXA),PI1(MXA)
      COMPLEX D(20000),Y,REFREL,XI,XI0,XI1,AN,BN,S1(MXA2),S2(MXA2)
      DOUBLE PRECISION PSI0,PSI1,PSI,DN,DX
      LOGICAL PRNT
      DX = X
      Y  = X*REFREL
C       -----------------------------------
C       SERIES TERMINATED AFTER NSTOP TERMS
C       -----------------------------------
      XSTOP = X + 4.*X**.3333+2.0
      NSTOP = XSTOP
      YMOD  = CABS(Y)
      NMX   = AMAX1(XSTOP,YMOD)+15
C
C     PRINT 554,NMX
C 554 FORMAT(' ..... NMX = ',I10)
C
C       ---------------------------------------------------
C       LOGARITHMIC DERIVATIVE D(J) CALCULATED BY DOWNWARD
C       RECURRENCE BEGINNING WITH INITIAL VALUE 0.0 + I*0.0
C       AT J = NMX
C       ---------------------------------------------------
      D(NMX) = CMPLX(0.0,0.0)
      NN     = NMX - 1
      DO 120 N=1,NN
        RN = NMX - N + 1
        D(NMX-N)=(RN/Y)-(1./(D(NMX-N+1)+RN/Y))
  120 CONTINUE
      DO 666 J=1,NANG
        PI0(J) = 0.0
        PI1(J) = 1.0
  666 CONTINUE
      NN = 2*NANG -1
      DO 777 J=1,NN
        S1(J) = CMPLX(0.0,0.0)
        S2(J) = CMPLX(0.0,0.0)
  777 CONTINUE
C       ---------------------------------------------
C       RICCATI-BESSEL FUNCTIONS WITH REAL ARGUMENT X
C       CALCULATED BY UPWAR RECURRENCE
C       ---------------------------------------------
      PSI0 = DCOS(DX)
      PSI1 = DSIN(DX)
      CHI0 = -SIN(X)
      CHI1 =  COS(X)
      APSI0 = PSI0
      APSI1 = PSI1
      XI0 = CMPLX(APSI0,-CHI0)
      XI1 = CMPLX(APSI1,-CHI1)
      QSCA = 0.0
      N = 1
 200  DN = N
      RN = N
      FN = (2.*RN+1.)/(RN*(RN+1.))
      PSI = (2.*DN-1.)*PSI1/DX-PSI0
      APSI = PSI
      CHI = (2.*RN-1.)*CHI1/X - CHI0
      XI = CMPLX(APSI,-CHI)
C
      AN = (D(N)/REFREL+RN/X)*APSI - APSI1
      AN = AN/((D(N)/REFREL+RN/X)*XI - XI1)
C
      BN = (REFREL*D(N)+RN/X)*APSI - APSI1
      BN = BN/((REFREL*D(N)+RN/X)*XI - XI1)
C
      QSCA = QSCA+(2.*RN+1.)*(CABS(AN)*CABS(AN)+CABS(BN)*CABS(BN))
      DO 789 J=1,NANG
        JJ = 2*NANG - J
        PI(J) = PI1(J)
        TAU(J) = RN*AMU(J)*PI(J) - (RN+1.)*PI0(J)
        P = (-1.)**(N-1)
        S1(J) = S1(J) + FN*(AN*PI(J)+BN*TAU(J))
        T = (-1.)**N
        S2(J) = S2(J) + FN*(AN*TAU(J)+BN*PI(J))
        IF( J .EQ. JJ ) GOTO 789
        S1(JJ) = S1(JJ)+FN*(AN*PI(J)*P+BN*TAU(J)*T)
        S2(JJ) = S2(JJ)+FN*(AN*TAU(J)*T+BN*PI(J)*P)
  789 CONTINUE
      PSI0 = PSI1
      PSI1 = PSI
      APSI1 = PSI1
      CHI0 = CHI1
      CHI1 = CHI
      XI1 = CMPLX(APSI1,-CHI1)
      N = N + 1
      RN = N
      DO 999 J=1,NANG
        PI1(J) = ((2.*RN-1.)/(RN-1.))*AMU(J)*PI(J)
        PI1(J) = PI1(J) -RN*PI0(J)/(RN-1.)
        PI0(J) = PI(J)
  999 CONTINUE
      IF( N-1-NSTOP) 200,300,300
  300 QSCA = (2./(X*X))*QSCA
      QEXT = (4./(X*X))*REAL(S1(1))
      QBACK = (4./(X*X))*CABS(S1(2*NANG-1))*CABS(S1(2*NANG-1))
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE PLTPSC(P,NMAX,X,SIZP,WAVEL,QSCA,QEXT,
     +                  XLAB,YLAB)
C
      DIMENSION P(NMAX),X(NMAX)
      CHARACTER*80 LABEL
      LOGICAL PRNT
      DATA PRNT / .TRUE. /
C
C     SET GRID
C
      CALL AGSETF('GRID/LEFT.'  ,.10)
      CALL AGSETF('GRID/RIGHT.' ,.70)
      CALL AGSETF('GRID/BOTTOM.',.10)
      CALL AGSETF('GRID/TOP.'   ,.80)
C
C     CHANGE MAXIMUM LINE LENGTH.
C
      CALL AGSETI('LINE/MAX.',80)
C
C     SET X MINIMUM AND MAXIMUM
C
      CALL AGSETF('X/MIN.',  0.0 )
      CALL AGSETF('X/MAX.',180.0 )
      CALL AGSETI('X/ORD.',    0 )
      CALL AGSETI('X/NICE.',   0 )
      CALL AGSETI('X/LOG.',0)
C
C     SET Y MINIMUM AND MAXIMUM
C
C
C     SET Y MINIMUM AND MAXIMUM
C
      CALL AGSETI('Y/LOG.',1)
      YMIN = P(1)
C
      DO 10 N=2,NMAX
C
        IF( P(N) .LT. YMIN ) YMIN = P(N)
C
  10  CONTINUE
      YMIN = YMIN * 0.5
      write(6,3333) nmax,ymin
 3333 format(' nmax and ymin = ',i5,2x,1pe11.4)
c
      CALL AGSETF('Y/MIN.', YMIN )
C
      YMAX = P(1)
C
      DO 20 N=2,NMAX
C
        IF( P(N) .GT. YMAX ) YMAX = P(N)
C
  20  CONTINUE
      YMAX = YMAX * 1.5
      write(6,3334) nmax,ymax
 3334 format(' nmax and ymax = ',i5,2x,1pe11.4)
      CALL AGSETF('Y/MAX.', YMAX )
      CALL AGSETI('Y/ORD.',    0 )
      CALL AGSETI('Y/NICE.',   0 )
C
C
C       SET LABELS:
C
      CALL AGSETC('LABEL/NAME.','L')
      CALL AGSETI('LINE/NUMBER.',100)
      CALL AGSETC('LINE/TEXT.',YLAB)
C
      CALL AGSETC('LABEL/NAME.','B')
      CALL AGSETI('LINE/NUMBER.',-100)
      CALL AGSETC('LINE/TEXT.',XLAB)
C
      CALL AGSETC('LABEL/NAME.','T')
C
      CALL AGSETI('LINE/NUMBER.', 100)
      WRITE(UNIT=LABEL,FMT=100)
  100 FORMAT(' ... SCATTERING PHASE FUNCTION  ......$')
      CALL AGSETC('LINE/TEXT.',LABEL)
C
      CALL AGSETI('LINE/NUMBER.',  70)
      WRITE(UNIT=LABEL,FMT=101) WAVEL,SIZP
  101 FORMAT(' WAVELENGTH = ',F7.4,' MICRONS  X = ',1PE11.3,'$')
      CALL AGSETC('LINE/TEXT.',LABEL)
C
      CALL AGSETI('LINE/NUMBER.',  40)
      WRITE(UNIT=LABEL,FMT=102) QSCA,QEXT
  102 FORMAT(' QSCA =',1PE11.3,' QEXT = ',1PE11.3,'$')
      CALL AGSETC('LINE/TEXT.',LABEL)
C
C         SET BACKGROUND FRAME AND PARAMETERS.
C
      CALL AGSTUP
      CALL AGBACK
      CALL AGBACK
C
C
      CALL setusv('IN',2000)
      CALL setusv('LW', 200)
C
      CALL CURVE(X,P,NMAX)
C
      IF( PRNT ) THEN
        WRITE(6,FMT=110)
  110   FORMAT(' ... SCATTERING PHASE FUNCTION  ......$')
C
        WRITE(6,FMT=111) WAVEL,SIZP
  111   FORMAT(' WAVELENGTH = ',F7.4,' MICRONS  X = ',1PE11.3,'$')
C
        WRITE(6,FMT=112) QSCA,QEXT
  112   FORMAT(' QSCA =',1PE11.3,' QEXT = ',1PE11.3,'$')
C
        WRITE(6,FMT=113)
  113   FORMAT(' INDEX     ANGLE         P ')
C
        DO 115 N=1,NMAX
         WRITE(6,116) N,X(N),P(N)
  116    FORMAT(1X,I3,1X,F6.2,2X,1PE11.4)
  115   CONTINUE
C
      ENDIF
C
      CALL FRAME
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE PLTQ  (P,NMAX,X,WAVEL,REFRE,REFIM,
     +                  XLAB,YLAB)
C
      DIMENSION P(NMAX),X(NMAX)
      CHARACTER*100 LABEL
C
C     SET GRID
C
      CALL AGSETF('GRID/LEFT.'  ,.10)
      CALL AGSETF('GRID/RIGHT.' ,.70)
      CALL AGSETF('GRID/BOTTOM.',.10)
      CALL AGSETF('GRID/TOP.'   ,.80)
C
C     CHANGE MAXIMUM LINE LENGTH.
C
      CALL AGSETI('LINE/MAX.',80)
C
C     SET X MINIMUM AND MAXIMUM
C
      CALL AGSETF('X/MIN.',X(1)   )
      CALL AGSETF('X/MAX.',X(NMAX))
      CALL AGSETI('X/LOG.',1)
      CALL AGSETI('X/ORD.',    0  )
      CALL AGSETI('X/NICE.',   0  )
C
      CALL AGSETI('Y/LOG.',0)
      CALL AGSETF('Y/MIN.',  0.0 )
      CALL AGSETF('Y/MAX.',  4.0 )
C
C       SET LABELS:
C
      CALL AGSETC('LABEL/NAME.','L')
      CALL AGSETI('LINE/NUMBER.',100)
      CALL AGSETC('LINE/TEXT.',YLAB)
C
      CALL AGSETC('LABEL/NAME.','B')
      CALL AGSETI('LINE/NUMBER.',-100)
      CALL AGSETC('LINE/TEXT.',XLAB)
C
      CALL AGSETC('LABEL/NAME.','T')
      CALL AGSETI('LINE/NUMBER.', 100)
      WRITE(UNIT=LABEL,FMT=100)
  100 FORMAT(' ... SCATTERING EFFICIENCY ....... $')
      CALL AGSETC('LINE/TEXT.',LABEL)
C
      CALL AGSETI('LINE/NUMBER.',  70)
      WRITE(UNIT=LABEL,FMT=101) WAVEL
  101 FORMAT(' WAVELENGTH = ',F7.4,'$')
      CALL AGSETC('LINE/TEXT.',LABEL)
C
      CALL AGSETI('LINE/NUMBER.',  40)
      WRITE(UNIT=LABEL,FMT=102) REFRE,REFIM
  102 FORMAT('  RFINX REAL = ',1PE11.4,
     +       '  RFINX  IMG = ',1PE11.4,'$')
      CALL AGSETC('LINE/TEXT.',LABEL)
C
C
      CALL AGSETI('Y/ORD.',    0 )
      CALL AGSETI('Y/NICE.',   0 )
C
C         SET BACKGROUND FRAME AND PARAMETERS.
C
      CALL AGSTUP
      CALL AGBACK
      CALL AGBACK
C
      call setusv('LW',200 )
      call setusv('IN',2000)
C
      CALL CURVE(X,P,NMAX)
C
      CALL FRAME
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE PLTR  (P,NMAX,X,RE,B,
     +                  XLAB,YLAB)
C
      DIMENSION P(NMAX),X(NMAX)
      CHARACTER*80 LABEL
C
C     SET GRID
C
      CALL AGSETF('GRID/LEFT.'  ,.10)
      CALL AGSETF('GRID/RIGHT.' ,.70)
      CALL AGSETF('GRID/BOTTOM.',.10)
      CALL AGSETF('GRID/TOP.'   ,.80)
C
C     CHANGE MAXIMUM LINE LENGTH.
C
      CALL AGSETI('LINE/MAX.',80)
C
C     SET X MINIMUM AND MAXIMUM
C
      DO 11 N=1,NMAX
         IF( P(N) .LE. 1.E-10) P(N) = 1.E-10
   11 CONTINUE
C
      CALL AGSETF('X/MIN.',X(1)   )
      CALL AGSETF('X/MAX.',X(NMAX))
C
      CALL AGSETI('X/LOG.',1)
C
C     SET Y MINIMUM AND MAXIMUM
C
      CALL AGSETI('Y/LOG.',1)
      YMIN = P(1)
C
      DO 20 N=2,NMAX
C
        IF( P(N) .LT. YMIN ) YMIN = P(N)
C
  20  CONTINUE
      YMIN = YMIN * 0.5
      CALL AGSETF('Y/MIN.', YMIN )
C
      YMAX = P(1)
C
      DO 30 N=2,NMAX
C
        IF( P(N) .GT. YMAX ) YMAX = P(N)
C
  30  CONTINUE
      YMAX = YMAX * 1.5
      CALL AGSETF('Y/MAX.', YMAX )
      CALL AGSETI('Y/ORD.',    0 )
      CALL AGSETI('Y/NICE.',   0 )
C
C       SET LABELS:
C
      CALL AGSETC('LABEL/NAME.','L')
      CALL AGSETI('LINE/NUMBER.',100)
      CALL AGSETC('LINE/TEXT.',YLAB)
C
      CALL AGSETC('LABEL/NAME.','B')
      CALL AGSETI('LINE/NUMBER.',-100)
      CALL AGSETC('LINE/TEXT.',XLAB)
C
      CALL AGSETC('LABEL/NAME.','T')
      CALL AGSETI('LINE/NUMBER.', 100)
      WRITE(UNIT=LABEL,FMT=100)
  100 FORMAT(' ... NORMALIZED SIZE DISTRIBUTION ... $')
      CALL AGSETC('LINE/TEXT.',LABEL)
C
      CALL AGSETI('LINE/NUMBER.',  70)
      WRITE(UNIT=LABEL,FMT=101) RE,B
  101 FORMAT(' RE = ',F10.6,'  B = ',1PE11.4,'$')
      CALL AGSETC('LINE/TEXT.',LABEL)
C
      CALL AGSETI('LINE/NUMBER.',  40)
      WRITE(UNIT=LABEL,FMT=102)
  102 FORMAT('$')
      CALL AGSETC('LINE/TEXT.',LABEL)
C
C         SET BACKGROUND FRAME AND PARAMETERS.
C
      CALL AGSTUP
      CALL AGBACK
      CALL AGBACK
C
      call setusv('LW',200 )
      call setusv('IN',2000)
C
      CALL CURVE(X,P,NMAX)
C
      CALL FRAME
C
      PRINT 500
  500 FORMAT(' ... DONE WITH PLTR '/)
C
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE PLTPEF(P,PHG,NMAX,X,WAVEL,GF,GB,F,RE,B,
     +                  XLAB,YLAB)
C
      DIMENSION P(NMAX),PHG(NMAX),X(NMAX)
      CHARACTER*80 LABEL
C
C     SET GRID
C
      CALL AGSETF('GRID/LEFT.'  ,.10)
      CALL AGSETF('GRID/RIGHT.' ,.70)
      CALL AGSETF('GRID/BOTTOM.',.10)
      CALL AGSETF('GRID/TOP.'   ,.80)
C
C     CHANGE MAXIMUM LINE LENGTH.
C
      CALL AGSETI('LINE/MAX.',80)
C
C     SET X MINIMUM AND MAXIMUM
C
      CALL AGSETF('X/MIN.',X(1)   )
      CALL AGSETF('X/MAX.',X(NMAX))
      CALL AGSETI('X/LOG.',0)
      CALL AGSETI('X/ORD.',    0  )
      CALL AGSETI('X/NICE.',   0  )
C
C     SET Y MINIMUM AND MAXIMUM
C
      CALL AGSETI('Y/LOG.',1)
      YMIN = P(1)
      IF( PHG(1) .LT. YMIN ) YMIN = PHG(1)
C
      DO 10 N=2,NMAX
C
        IF( P(N) .LT. YMIN ) YMIN = P(N)
        IF( PHG(N) .LT. YMIN ) YMIN = PHG(N)
C
  10  CONTINUE
      YMIN = YMIN * 0.5
      CALL AGSETF('Y/MIN.', YMIN )
C
      YMAX = P(1)
      IF( PHG(1) .GT. YMAX ) YMAX = PHG(1)
C
      DO 20 N=2,NMAX
C
        IF( P(N) .GT. YMAX ) YMAX = P(N)
        IF( PHG(N) .GT. YMAX ) YMAX = PHG(N)
C
  20  CONTINUE
      YMAX = YMAX * 1.5
      CALL AGSETF('Y/MAX.', YMAX )
      CALL AGSETI('Y/ORD.',    0 )
      CALL AGSETI('Y/NICE.',   0 )
C
C       SET LABELS:
C
      CALL AGSETC('LABEL/NAME.','L')
      CALL AGSETI('LINE/NUMBER.',100)
      CALL AGSETC('LINE/TEXT.',YLAB)
C
      CALL AGSETC('LABEL/NAME.','B')
      CALL AGSETI('LINE/NUMBER.',-100)
      CALL AGSETC('LINE/TEXT.',XLAB)
C
      CALL AGSETC('LABEL/NAME.','T')
      CALL AGSETI('LINE/NUMBER.', 100)
      WRITE(UNIT=LABEL,FMT=100)
  100 FORMAT(' ... MEAN SCAT AND HENYEY-GREENSTEIN PH FNC... $')
      CALL AGSETC('LINE/TEXT.',LABEL)
C
      CALL AGSETI('LINE/NUMBER.',  70)
      WRITE(UNIT=LABEL,FMT=101) GF,GB,F
  101 FORMAT('  HG GF =',F8.6,' HG GB = ',F8.6,'  F = ',F8.6,'$')
      CALL AGSETC('LINE/TEXT.',LABEL)
C
      CALL AGSETI('LINE/NUMBER.',  40)
      WRITE(UNIT=LABEL,FMT=102) WAVEL,RE,B
  102 FORMAT(' WAVELENGTH = ',F7.4,' RE = ',F10.6,'  B = ',1PE11.4,'$')
      CALL AGSETC('LINE/TEXT.',LABEL)
C
C
      CALL AGSETI('Y/ORD.',    0 )
      CALL AGSETI('Y/NICE.',   0 )
C
C         SET BACKGROUND FRAME AND PARAMETERS.
C
      CALL AGSTUP
      CALL AGBACK
      CALL AGBACK
C
      call setusv('LW',200 )
      call setusv('IN',2000)
C
      CALL CURVE(X,P,NMAX)
C
      CALL POINTS(X,PHG,NMAX,'+',0)
C
      CALL FRAME
C
      PRINT 500
  500 FORMAT(' ... DONE WITH PLTPEF '/)
C
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE PLTQEF(P,NMAX,X,WAVEL,RE,B,
     +                  XLAB,YLAB)
C
      DIMENSION P(NMAX),X(NMAX)
      CHARACTER*80 LABEL
C
C     SET GRID
C
      CALL AGSETF('GRID/LEFT.'  ,.10)
      CALL AGSETF('GRID/RIGHT.' ,.70)
      CALL AGSETF('GRID/BOTTOM.',.10)
      CALL AGSETF('GRID/TOP.'   ,.80)
C
C     CHANGE MAXIMUM LINE LENGTH.
C
      CALL AGSETI('LINE/MAX.',80)
C
C     SET X MINIMUM AND MAXIMUM
C
      CALL AGSETF('X/MIN.',X(1)   )
      CALL AGSETF('X/MAX.',X(NMAX))
      CALL AGSETI('X/LOG.',0)
      CALL AGSETI('X/ORD.',    0  )
      CALL AGSETI('X/NICE.',   0  )
C
C     SET Y MINIMUM AND MAXIMUM
C
      CALL AGSETI('Y/LOG.',1)
      YMIN = P(1)
C
      DO 10 N=2,NMAX
C
        IF( P(N) .LT. YMIN ) YMIN = P(N)
C
  10  CONTINUE
      YMIN = YMIN * 0.5
      CALL AGSETF('Y/MIN.', YMIN )
C
      YMAX = P(1)
C
      DO 20 N=2,NMAX
C
        IF( P(N) .GT. YMAX ) YMAX = P(N)
C
  20  CONTINUE
      YMAX = YMAX * 1.5
      CALL AGSETF('Y/MAX.', YMAX )
      CALL AGSETI('Y/ORD.',    0 )
      CALL AGSETI('Y/NICE.',   0 )
C
C       SET LABELS:
C
      CALL AGSETC('LABEL/NAME.','L')
      CALL AGSETI('LINE/NUMBER.',100)
      CALL AGSETC('LINE/TEXT.',YLAB)
C
      CALL AGSETC('LABEL/NAME.','B')
      CALL AGSETI('LINE/NUMBER.',-100)
      CALL AGSETC('LINE/TEXT.',XLAB)
C
      CALL AGSETC('LABEL/NAME.','T')
      CALL AGSETI('LINE/NUMBER.', 100)
      WRITE(UNIT=LABEL,FMT=100)
  100 FORMAT(' ... MEAN SCATTERING PHASE FUNCTION ... $')
      CALL AGSETC('LINE/TEXT.',LABEL)
C
      CALL AGSETI('LINE/NUMBER.',  70)
      WRITE(UNIT=LABEL,FMT=101) WAVEL
  101 FORMAT(' WAVELENGTH = ',F7.4,'$')
      CALL AGSETC('LINE/TEXT.',LABEL)
C
      CALL AGSETI('LINE/NUMBER.',  40)
      WRITE(UNIT=LABEL,FMT=102) RE,B
  102 FORMAT(' RE = ',F10.6,'  B = ',1PE11.4,'$')
      CALL AGSETC('LINE/TEXT.',LABEL)
C
C
      CALL AGSETI('Y/ORD.',    0 )
      CALL AGSETI('Y/NICE.',   0 )
C
C         SET BACKGROUND FRAME AND PARAMETERS.
C
      CALL AGSTUP
      CALL AGBACK
      CALL AGBACK
C
      call setusv('LW',200 )
      call setusv('IN',2000)
C
      CALL CURVE(X,P,NMAX)
C
      CALL FRAME
C
      PRINT 500
  500 FORMAT(' ... DONE WITH PLTQEF '/)
C
      RETURN
      END
      subroutine setgks
      dimension iasf(13)
      data iasf /13*1/
C
C OPEN GKS
C
      call opngks
C
C TURN OFF THE CLIPPING INDICATOR.
C
      call gsclip(0)
C
C SET ALL ASPECT SOURCE FLAGS TO "INDIVIDUAL".
C
      call gsasf(iasf)
C
C FORCE SOLID FILL.
C
      call gsfais(1)
C
C SET UP FLASHBUFFER(S)
C
      call gopwk(9,3,3)
C
C RETURN TO CALLER
C
      return 
      end
      subroutine gauss(gauspt,gauswt,xzero,pzero,maxmu)
c
      parameter (nmax =10001)
c
c  ......... compute gaussian points and weights
c
      real gauspt(maxmu),gauswt(maxmu)
      real xzero(maxmu) ,pzero(maxmu)
c
      real mu(nmax),p(nmax),pscrt(nmax)
      real gptrev(nmax),gwtrev(nmax)
c
c                  define evenly spaced mesh points over (-1,1) for
c                  evaluation of legendre polynomials:
c
      mu(1)    = -1.0
      mu(nmax) = +1.0
      deltmu   = (mu(nmax) - mu(1)) / real(nmax-1)
      do 100 n=2,nmax-1
        mu(n) = mu(n-1) + deltmu
  100 continue
c
c         compute legendre polynomials of the evenly spaced
c         grid using recursion relations:
c
       do 200 n=1,nmax
c
         x = mu(n)
c
         p(n) = plgndr(maxmu,x)
c
  200  continue
c
c
       print 300,maxmu
  300  format(/'  compute gaussian points and weights for maxord = ',i4)
c
       call fndzro(maxmu,p,nmax,mu,gptrev,pzero,pscrt)
       call fndgwt(maxmu,p,pscrt,nmax,mu,gptrev,gwtrev)
c
c reverse order of points and weights:
c
       gauspt(1) = 1.0
       gauswt(1) = 0.0
c
       do 400 n=1,maxmu
         gauspt(n+1) = gptrev(maxmu-n+1)
         gauswt(n+1) = gwtrev(maxmu-n+1)
  400  continue
c
       gauspt(maxmu+2) = -1.0
       gauswt(maxmu+2) =  0.0
c
       write(6,499) maxmu+2
  499  format(/' quadrature point and weight number = ',i5/
     $         ' n       point       weight    ')
       do 500 n=1,maxmu+2
         write(6,501) n,gauspt(n),gauswt(n)
  501    format(2x,i3,2x,2(f12.8,3x))
  500  continue
c
      return
c
      end
      function plgndr(maxord,x)
c
c
        if( maxord  .eq. 0 ) then
           plgndr = 1.0
        else if( maxord .eq. 1 ) then
           plgndr =  x
        else
           p0 = 1.0
           p1 =  x
c
           do 100 l=2,maxord
c
c .......   use pure recursion relation to compute legendre polynomials
c
              rl       = real(l-1)
              plgndr   = ((2.*rl+1.)*x*p1 - rl*p0) / (rl+1.)
              p0       = p1
              p1       = plgndr
c
  100      continue
c
        endif
c
      return
      end
      subroutine fndzro(maxmu,p,nmax,x,xzero,pzero,ppzero)
c
      dimension p(nmax),x(nmax),xzero(maxmu),pzero(maxmu),ppzero(1)
      parameter (numdiv = 1000, maxitr = 1000, difmax = 1.e-14)
      logical incrs,decrs,prnt
      data prnt / .false./
c
c     find zeros of legendre polynomial by marching through domain;
c     when x axis crossing detected, stop and use newton-raphson
c     technique to find zeros:
c
      izero = 0
      do 100 n=1,nmax-1
        xmin = x(n)
        xmax = x(n+1)
        pmin = p(n)
        pmax = p(n+1)
        incrs = (pmin .lt. 0.) .and. (pmax .gt. 0.)
        decrs = (pmin .gt. 0.) .and. (pmax .lt. 0.)
        if( incrs .or. decrs ) then
          iter  = 1
          xbar  = (xmin + xmax) / 2.
          pmid  = plgndr(maxmu,xbar)
          deltx = (xmax - xmin) / real(numdiv)
  150     if( abs(pmid) .gt. difmax ) then
            xm = xbar - deltx
            xp = xbar + deltx
            ppmid = (plgndr(maxmu,xp) - plgndr(maxmu,xm)) / (2.*deltx)
            xbar  = (xm + xp) / 2.
            pmid  = plgndr(maxmu,xbar)
            iter  = iter + 1
            if( iter .gt. maxitr ) go to 200
            xbar = xbar - (pmid/ppmid)
            goto 150
          endif
  200     izero = izero + 1
          xzero(izero)  = xbar
          pzero(izero)  = pmid
          ppzero(izero) = ppmid
        endif
c
  100 continue
c
      if( prnt ) then
         print 300,maxmu,izero,iter,maxitr
  300    format(/' legendre polynomial of order = ',i3/
     +           ' number of zeros found = ',i3,
     +           ' iterations needed = ',i3,'  out of possible = ',i3)
         do 350 i=1,izero
           print 355,i,xzero(i),pzero(i),ppzero(i)
  355      format('  zero = ',i3,'  mu value = ',f10.8,'  p value = ',
     +            f10.8,'  pp value = ',f10.5)
  350    continue
      endif
c
      return
      end
      subroutine fndgwt(maxmu,p,pp,nmax,x,xzero,gaussw)
      logical prnt,bad
      data prnt / .true.  /
      data diftol / 1.e-14  /
c
      dimension p(nmax),x(nmax),pp(1),xzero(1),gaussw(1)
c
      do 100 iw=1,maxmu
        sum = 0.
        do 200 n=2,nmax-1,2
          difmin = x(n-1) - xzero(iw)
          difmid = x(n)   - xzero(iw)
          difmax = x(n+1) - xzero(iw)
          bad = .false.
          if( abs(difmin) .lt. diftol ) bad = .true.
          if( abs(difmid) .lt. diftol ) bad = .true.
          if( abs(difmax) .lt. diftol ) bad = .true.
          if( .not. bad ) then
            sum = sum + ((x(n+1) - x(n-1))/6.)
     +                * (    (p(n-1) / difmin)
     +                +   4.*(p(n)   / difmid)
     +                +      (p(n+1) / difmax)  )
          endif
  200   continue
        gaussw(iw) = sum / pp(iw)
  100 continue
c
c
      if( prnt ) then
        print 208
  208     format('         num              zero value   ',16x,
     +            '    gaussian weight ')
        do 209 iw=1,maxmu
          print 210,iw,xzero(iw),gaussw(iw)
  210     format(9x,i2,10x,f17.14,19x,f17.14)
  209   continue
      endif
      if( prnt ) then
        sum = 0.0
        do 215 iw=1,maxmu
          sum = sum + gaussw(iw)
  215   continue
        print 220,sum
  220   format(/' sum of gaussian weights  = ',f17.14)
      endif
c
      return
      end
