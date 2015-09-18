      PROGRAM MIEBPB
C
      PARAMETER (MAXMU=200)
      COMMON /RES/ MXANG,MXANG2
      DIMENSION P(MAXMU),SCTAN(MAXMU),PEFF(MAXMU)
      DIMENSION PHGFOR(MAXMU),PHGBCK(MAXMU),PHGTOT(MAXMU)
C
      PARAMETER (NMBRAD=1000, NMBWAV=1)
      DIMENSION RAD(NMBRAD),WAVEL(NMBWAV)
      DIMENSION QSCAT(NMBRAD),QEXT(NMBRAD),SIZEP(NMBRAD)
      REAL      NMBPRT(NMBRAD)
C
      LOGICAL PRNT,PLOT
      REAL N0,NMIN,NMAX,MUMIN,MUMAX,MUMN
      DATA    PRNT / .FALSE./
      DATA    PLOT / .FALSE./
      INTEGER FRQPLT,FRQPRT
      DATA    FRQPLT /  50 /
      DATA    FRQPRT /1000 /
C
      RNMDEN(R,RE,B) = (R**((1.-3*B)/B))*EXP(-R/(B*RE))
      PHG(G,U) = (1.-(G*G))/((1.+(G*G)-(2.*G*U))**1.5)
c
c..... setup gks:
c
      call setgks
C
C...... SET DIMENSIONS:
C
      MXANG  = MAXMU
      MXAO2  = 91
      MXANG2 = (MXAO2*2) - 1
      DO 5 NA=1,MXANG2
        PEFF(NA) = 0.0
    5 CONTINUE
C
C...... SET RADIUS, REAL AND IMAG INDEX REFRACTION:
C...... (FOR .500MICRO-METERS, RAD0=.00256;
C......  FOR .625MICRO-METERS, RAD0=.0032;
C......  FOR  3.7MICRO-METERS, RAD0=.018944)
C......  FOR 10.0MICRO-METERS, RAD0=.051200)
C .....
C .....  ??? RAD0:
C
      RAD0      = .002560
      SIGRAD    =  200.
      N0        =  100.
C
C..... ?????? RE AND B:
C
      BVAL      =  1./6.
      REVAL     =  1.0
C
      SUMDEN    =  0.0
      PIE       = 4. * ATAN(1.)
C
C...... IF MORE THAN ONE RADIUS ....
C
      IF( NMBRAD .GT. 1 )  THEN
C
      DO 10 N=1,NMBRAD
        RAD(N) = RAD0 * (10.**((REAL(N)-N0)/SIGRAD))
        RADIUS = RAD(N)
        NMBPRT(N) = RNMDEN(RADIUS,REVAL,BVAL)
   10 CONTINUE
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
   15 CONTINUE
C
      PRINT 17,SUMDEN
   17 FORMAT(/' ... TRPZD QUADRATURE OF NUMBER DEN = ',1PE11.4)
C
      DO 20 N=1,NMBRAD
        NMBPRT(N) = NMBPRT(N) / SUMDEN
   20 CONTINUE
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
   30 CONTINUE
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
   40 CONTINUE
C
      SIGMN = SQRT(SUM1)
      AMEAN = PIE * SUM2
      RAREA = SQRT(SUM2)
      VMEAN = (4.*PIE/3.) * SUM3
C
      REMAX = REVAL * (1.-3.*BVAL)
C
      PRINT 45,REVAL,BVAL,REMAX,RADMN,SIGMN,RAREA,AMEAN,VMEAN
   45 FORMAT(/'  .... INFO ON SIZE DISTRIBUTION ....'/
     + ' ....   RE        B          RMAX        RMEAN        SIGMA '/
     +        2X,5(1PE11.4,2X)/
     + ' ....                        RAREA        AMEAN        VMEAN '/
     +       28X,3(1PE11.4,2X))
C
C...... IF ONLY ONE RADIUS ....
C
      ELSE IF( NMBRAD .EQ. 1 )  THEN
C
        RAD(1)    = RAD0
        NMBPRT(1) =  1.0
        PRINT 50,RAD0
   50   FORMAT(/' .... ONLY ONE RADIUS DESIRED = ',1PE11.4,' MICRONS'/)
C
      ENDIF
C
C...... ??? SET WAVELENGTH IN MICRONS
C
      WAVEL(1)  =   .550
C
C.....  ??? SET OPTICAL INDEX OF REFRACTION (REAL AND IMAGINARY):
C
C...          0.55 MICRONS
C
      REFRE  = 1.333
      REFIM  = 1.96E-9
C
C...         10.00 MICRONS
C
C     REFRE  = 1.218
C     REFIM  = 0.0508
C
C...          3.70 MICRONS
C
C     REFRE  = 1.374
C     REFIM  = 3.60E-3
C
C...           .625 MICRONS
C
C     REFRE  = 1.332
C     REFIM  = 1.39E-8
C
C..... PERFORM MIE COMPUTATION ....
C
      DO 100 NW=1,NMBWAV
C
        PRINT 101,WAVEL,REFRE,REFIM
  101   FORMAT(/' .... WAVELENGTH = ',1PE11.4,' MICRO-METERS'/
     +   ' .... REAL INDEX OF REFRACTION = ',1PE11.4/
     +   ' .... IMG  INDEX OF REFRACTION = ',1PE11.4/)
C
        DO 200 NR=1,NMBRAD
C
          RADIUS = RAD(NR)
C
          PLOT = .FALSE.
          IF( MOD(NR,FRQPLT) .EQ. 0 ) PLOT = .TRUE.
          IF(  NR .EQ. 1 .OR. NR .EQ. NMBRAD ) PLOT = .TRUE.
C
          PRNT = .FALSE.
C
          IF( MOD(NR,FRQPRT) .EQ. 0 ) PRNT = .TRUE.
          IF(  NMBRAD .EQ. 1 ) PRNT = .TRUE.
C
          CALL CALLBH(MXAO2,WAVEL(NW),RAD(NR),
     +                REFRE,REFIM,P,SCTAN,
     +                X,QSCA,QEXTN,PRNT,PLOT)
C
          QSCAT(NR) = QSCA
          QEXT (NR) = QEXTN
          SIZEP(NR) = X
C
          IF( NMBRAD .EQ. 1 ) THEN
            PRINT 199,SIZEP(NR),QSCAT(NR),QEXT(NR)
  199       FORMAT('  SIZEP   QSCAT  QEXT =',2X,3(1PE12.5,2X))
          ENDIF
C
C....... COMPUTE CONTRIBUTION TO EFFECTIVE
C....... SCATTERING PHASE FUNCTION .......
C
      DO 150 NA=1,MXANG2
        PEFF(NA) = PEFF(NA) +
     +             P(NA)*QSCAT(NR)*PIE*RADIUS*RADIUS*NMBPRT(NR)
  150 CONTINUE
C
  200   CONTINUE
  100 CONTINUE
C
C..... STOP IF ONLY ONE RADIUS
C
      IF( NMBRAD .EQ. 1 ) THEN
        PRINT 201
  201   FORMAT(/' .... END OF ONE RADIUS MIE PROGRAM ....'/)
        STOP
      ENDIF
C
C..... CONTINUE WITH MORE THAN ONE RADIUS ...
C
      CALL PLTQ  (QSCAT,NMBRAD,SIZEP,WAVEL,REFRE,REFIM,
     +            '    SIZE PARAMETER   $',
     +            ' SCATTERING EFFICIENCY    $')
C
C..... COMPUTE EFFECTIVE SCATTERING EFFICIENCY, ........
C..... SINGLE PARTICAL SCATTERING ALBEDO, AND ..........
C..... SCATTERING PHASE FUNCTION .....
C
      SIGSC= 0.0
      SIGEX= 0.0
      SUM  = 0.0
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
        SIGSC = SIGSC + ((SCMIN+SCMAX)/2.) * (RMAX-RMIN)
C
        SIGEX = SIGEX + ((EXMIN+EXMAX)/2.) * (RMAX-RMIN)
C
        SUM  = SUM  + (( NMIN+ NMAX)/2.) * (RMAX-RMIN)
C
  300 CONTINUE
C
      SIGSC = SIGSC / SUM
      QESC  = SIGSC / AMEAN
C
      SIGEX = SIGEX / SUM
      QEEX  = SIGEX / AMEAN
C
      WEFF  = SIGSC / SIGEX
C
C.... CONVERT CROSS SECTIONS FROM MICRO-METERS TO CM2:
C
      SIGSC = SIGSC * 1.E-8
      SIGEX = SIGEX * 1.E-8
C
      PRINT 310,BVAL,REVAL,QESC,SIGSC,QEEX,SIGEX,WEFF
  310 FORMAT(/' ... SIZE DISTRIBUTION  B AND RE = ',2(1PE11.4,2X)/
     +        ' ... EFFECTIVE SCATTERING EFFICIENCY    = ',1PE11.4/
     + ' ... SCATTERING CROSS SECTION           = ',1PE11.4,' CM2'/
     +        ' ... EFFECTIVE EXTINCTION EFFICIENCY    = ',1PE11.4/
     + ' ... EXTINCTION CROSS SECTION           = ',1PE11.4,' CM2'/
     +        ' ... EFFECTIVE SINGLE SCATTERING ALBEDO = ',1PE11.4)
C
      RENORM  = PEFF(1) / SIGSC
      PEFF(1) = 1.0
      DO 350 NA=2,MXANG2
        PEFF(NA) = PEFF(NA) / SIGSC
        PEFF(NA) = PEFF(NA) / RENORM
  350 CONTINUE
C
      CALL PLTR  (NMBPRT,NMBRAD,SIZEP,REVAL,BVAL,
     +            '   SIZE PARAMETER      $',
     +            ' NUMBER DENSITY / UNIT RADIUS$')
C
      CALL PLTR  (NMBPRT,NMBRAD,RAD,REVAL,BVAL,
     +            '   RADIUS (MICRONS)    $',
     +            ' NUMBER DENSITY / UNIT RADIUS$')
C
      CALL PLTQEF(PEFF,MXANG2,SCTAN,WAVEL,REVAL,BVAL,
     +            '  SCATTERING ANGLE     $',
     +            ' MEAN SCAT PHASE FUNCTION $')
C
C.... DETERMINE HENYEY-GREENSTEIN PROPERTIES:
C
C.... FIND MINIMUM IN EFFECTIVE SCATTERING PHASE FUNCTION:
C
      PMIN  = PEFF(1)
      NAMIN = 1
      DO 360 NA=2,MXANG2
        IF( PEFF(NA) .LT. PMIN ) THEN
           PMIN  = PEFF(NA)
           NAMIN =  NA
        ENDIF
  360 CONTINUE
C
C.... DETERMINE GF AND GB TO BE THE G VALUES COMPUTED
C.... FOR THE SMOOTHED MIE PHASE FUNCTION FOR THE FORWARD HEMISPHERE
C.... (FROM 0 TO SPECIAL ANGLE) AND FOR THE BACKWARD HEMISPHERE
C.... (FROM SPECIAL ANGLE TO 180) RESPECTIVELY:
C
      SUMG = 0.0
      SUM  = 0.0
      MAXG = 1
      SPCANG = 80.
      DO 365 NA=2,MXANG2
        IF( SCTAN(NA-1) .LE. SPCANG .AND. SPCANG .LT. SCTAN(NA) ) THEN
           MAXG = NA
        ENDIF
  365 CONTINUE
C
      PRINT 366,SPCANG,MAXG
  366 FORMAT(/' .... SPCANG AND MAXG = ',F7.2,2X,I6)
C
      DO 370 NA=2,MAXG,2
C
        PMIN  = PEFF(NA-1)
        MUMIN = COS(PIE*SCTAN(NA-1)/180.)
C
        PMN   = PEFF(NA)
        MUMN  = COS(PIE*SCTAN(NA)/180.)
C
        PMAX  = PEFF(NA+1)
        MUMAX = COS(PIE*SCTAN(NA+1)/180.)
C
        SUMG  = SUMG + ((MUMAX-MUMIN)/6.) *
     +                 (PMAX*MUMAX + PMIN*MUMIN + 4.*PMN*MUMN)
        SUM   = SUM  + ((MUMAX-MUMIN)/6.) *
     +                 (PMAX       + PMIN       + 4.*PMN     )
C
  370 CONTINUE
C
      GF     =  ABS( SUMG / SUM )
C
C
      SUMG = 0.0
      SUM  = 0.0
      DO 380 NA=MAXG+1,MXANG2-1,2
C
        PMIN  = PEFF(NA-1)
        MUMIN = COS(PIE*SCTAN(NA-1)/180.)
C
        PMN   = PEFF(NA)
        MUMN  = COS(PIE*SCTAN(NA)/180.)
C
        PMAX  = PEFF(NA+1)
        MUMAX = COS(PIE*SCTAN(NA+1)/180.)
C
        SUMG  = SUMG + ((MUMAX-MUMIN)/6.) *
     +                 (PMAX*MUMAX + PMIN*MUMIN + 4.*PMN*MUMN)
        SUM   = SUM  + ((MUMAX-MUMIN)/6.) *
     +                 (PMAX       + PMIN       + 4.*PMN     )
C
  380 CONTINUE
C
      GB     = ABS( SUMG / SUM )
C
C
      PRINT 362,SCTAN(NAMIN)
  362 FORMAT(/' ... MINIMUM SCT PHS FUNCTION ANGLE = ',F8.2,' DEGREES'/)
      UMIN = COS(PIE*SCTAN(NAMIN)/180.)
C
      ALPHA = (3.*GF*(1.-GF*GF)) /
     +        ((1. + (GF*GF) - (2.*GF*UMIN))**2.5)
C
      BETA  = (3.*GB*(1.-GB*GB)) /
     +        ((1. + (GB*GB) + (2.*GB*UMIN))**2.5)
C
      F  = BETA / (ALPHA + BETA)
      GM = F*GF + (1.-F)*GB
C
      PRINT 363,UMIN,ALPHA,BETA,F,GM
  363 FORMAT(/' .... FRACTION FORWARD COMPUTED BASED ON MINIMUM'/
     +        ' .... OF DOUBLE HG PHS FUNCTION '/
     +        ' .... UMIN  = ',F10.6/
     +        ' .... ALPHA = ',F10.6/
     +        ' .... BETA  = ',F10.6/
     +        ' .... F     = ',F10.8/
     +        ' .... GMEAN = ',F10.8)
C
      PRINT 399,GF,GB,F
  399 FORMAT(/' .... HENYEY GREENSTEIN PHASE FUNCTIONS ....'/
     +        ' .... G FORWARD = ',F10.7,' .... G BACKWARDS = ',F10.7/
     +        ' .... AND FOR F = ',F10.7/
     +      ' .... NA U  PHGFOR(NA)   PHGBCK(NA)   PHGTOT(NA) '/)
      SUMG = 0.0
      SUM  = 0.0
      DO 400 NA=1,MXANG2
        U = COS(PIE*SCTAN(NA)/180.)
        PHGFOR(NA) = PHG(GF,U)
        PHGBCK(NA) = PHG(GB,-U)
        PHGTOT(NA) = F*PHGFOR(NA) + (1.-F)*PHGBCK(NA)
C
        PRINT 401,NA,U,PHGFOR(NA),PHGBCK(NA),PHGTOT(NA)
  401   FORMAT(2X,I3,2X,F8.5,2X,3(1PE11.4,2X))
C
        IF( NA .EQ. MXANG2 ) GOTO 400
        PMIN  = PEFF(NA)
        MUMIN = COS(PIE*SCTAN(NA)/180.)
C
        PMAX  = PEFF(NA+1)
        MUMAX = COS(PIE*SCTAN(NA+1)/180.)
C
        SUMG  = SUMG + 0.5 * (PMAX*MUMAX + PMIN*MUMIN) * (MUMAX-MUMIN)
        SUM   = SUM  + 0.5 * (PMAX       + PMIN      ) * (MUMAX-MUMIN)
C
  400 CONTINUE
      GEFTRP = SUMG / SUM
C
      SUMG = 0.0
      SUM  = 0.0
      DO 500 NA=2,MXANG2-1,2
C
        PMIN  = PEFF(NA-1)
        MUMIN = COS(PIE*SCTAN(NA-1)/180.)
C
        PMN   = PEFF(NA)
        MUMN  = COS(PIE*SCTAN(NA)/180.)
C
        PMAX  = PEFF(NA+1)
        MUMAX = COS(PIE*SCTAN(NA+1)/180.)
C
        SUMG  = SUMG + ((MUMAX-MUMIN)/6.) *
     +                 (PMAX*MUMAX + PMIN*MUMIN + 4.*PMN*MUMN)
        SUM   = SUM  + ((MUMAX-MUMIN)/6.) *
     +                 (PMAX       + PMIN       + 4.*PMN     )
C
  500 CONTINUE
C
      GEFSMP = SUMG / SUM
C
C.... RENORMALIZE MIE SCATTERING PHASE FUNCTION:
C
      DO 450 NA=1,MXANG2
C
        PEFF(NA) = PEFF(NA) * (-2./SUM)
C
  450 CONTINUE
C
      PRINT 402,GEFTRP,GEFSMP,SUM
  402 FORMAT(/' .... ACTUAL G FOR SMOOTHED MIE SCT PHS FUNC'/
     +        ' .... USING TRAPEZOIDAL RULE  = ',F8.6/
     +        ' .... USING SIMPSONS    RULE  = ',F8.6/
     +        ' .... RENORMILIZATION VALUE     = ',1PE11.4)
C
      PRINT 458
  458 FORMAT(
     +  ' .... NA SCTANG  PHGTOT(NA)   PEFF(NA)   % DIFF ARE '/)
      DO 460 NA=1,MXANG2
C
        PCTDIF = 100. * (PHGTOT(NA) - PEFF(NA)) / PEFF(NA)
C
        PRINT 459,NA,SCTAN(NA),PHGTOT(NA),PEFF(NA),PCTDIF
  459   FORMAT(2X,I3,2X,F8.2,2X,3(1PE11.4,2X))
C
  460 CONTINUE
C
      CALL PLTPEF(PEFF,PHGTOT,MXANG2,SCTAN,
     +            WAVEL,GF,GB,F,REVAL,BVAL,
     +            '  SCATTERING ANGLE     $',
     +            '   SCAT PHASE FUNCTION $')
C
      PRINT 1000
 1000 FORMAT(/' .... END OF MIE PROGRAM ....'/)
c
      call gclwk(9)
      call clsgks
c
      END
C-------------------------------------------------------------------
      SUBROUTINE CALLBH(NANG,WAVEL,RAD,
     +                  REFRE,REFIM,P,SCTAN,
     +                  X,QSCA,QEXT,PRNT,PLOT)
C       ---------------------------------------------------------------
C       CALLBH CALCULATES THE SIZE PARAMETER (X) AND RELATIVE
C       REFRACTIVE INDEX (REFREL) FOR A GIVEN SPHERE REFRACTIVE
C       INDEX, MEDIUM REFRACTIVE INDEX, RADIUS, AND FREE SPACE
C       WAVELENGTH.  IT THEN CALLS BHMIE, THE SUBROUTINE THAT COMPUTES
C       AMPLITUDE SCATTERING MATRIX ELEMENTS AND EFFICIENCIES
C       ---------------------------------------------------------------
      COMPLEX REFREL,S1(400),S2(400)
      REAL P(1),SCTAN(1)
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
      CALL BHMIE(X,REFREL,NANG,S1,S2,QEXT,QSCA,QBACK,PRNT)
      IF( PRNT ) WRITE (6,65) QSCA,QEXT,QBACK
      IF( PRNT ) WRITE (6,17)
C       --------------------------------------------------
C       S33 AND S34 MATRIX ELEMNTS NORMALIZED BY S11.
C       S11 IS NORMALIZED TO 1.0 IN THE FORMWARD DIRECTION
C       POL=DEGREE OF POLARIZATION (INCIDENT UNPOLARIZED LIGHT)
C       --------------------------------------------------
      S11NOR = 0.5*(CABS(S2(1))**2+CABS(S1(1))**2)
      NAN    = 2*NANG - 1
      DANG   = 1.570796327/FLOAT(NANG-1)
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
        ANG = DANG*(AJ-1.)*57.2958
        SCTAN(J) = ANG
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
      SUBROUTINE BHMIE (X,REFREL,NANG,S1,S2,
     +                  QEXT,QSCA,QBACK,PRNT)
      PARAMETER (MXA=200,MXA2=400)
      DIMENSION AMU(MXA),THETA(MXA),PI(MXA),TAU(MXA),PI0(MXA),PI1(MXA)
      COMPLEX D(10000),Y,REFREL,XI,XI0,XI1,AN,BN,S1(MXA2),S2(MXA2)
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
      DANG  = 1.570796327/FLOAT(NANG-1)
      DO 555 J=1,NANG
        THETA(J) = (FLOAT(J)-1.)*DANG
        AMU(J)   = COS(THETA(J))
  555 CONTINUE
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
      INTEGER INTEN,SSIZE
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
      INTEN = 6500
      SSIZE =  200
C
      CALL setusv('IN',INTEN)
      CALL setusv('LW',SSIZE)
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
      INTEGER INTEN,SSIZE
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
C
      INTEN = 6500
      SSIZE =  200
C
      CALL setusv('IN',INTEN)
      CALL setusv('LW',SSIZE)
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
      INTEGER INTEN,SSIZE
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
C
      INTEN = 6500
      SSIZE =  200
C
      CALL setusv('IN',INTEN)
      CALL setusv('LW',SSIZE)
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
      INTEGER INTEN,SSIZE
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
C
      INTEN = 6500
      SSIZE =  200
C
      CALL setusv('IN',INTEN)
      CALL setusv('LW',SSIZE)
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
      INTEGER INTEN,SSIZE
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
C
      INTEN = 6500
      SSIZE =  200
C
      CALL setusv('IN',INTEN)
      CALL setusv('LW',SSIZE)
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
