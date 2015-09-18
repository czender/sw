      program mie_csz
      LOGICAL prnt
      data prnt /.true./ 
      NANG=11
      REFRE=1.3726
      REFIM=3.6641E-03
      RAD=2.
      WAVEL=3.7
c
      call CALLBH(NANG,WAVEL,RAD,
     +REFRE,REFIM,PRNT,
     +X,QSCA,QEXT)
c
      end
C-------------------------------------------------------------------
      SUBROUTINE CALLBH(NANG,WAVEL,RAD,
     +                  REFRE,REFIM,PRNT,
     +                  X,QSCA,QEXT)
C       ---------------------------------------------------------------
C       CALLBH CALCULATES THE SIZE PARAMETER (X) AND RELATIVE
C       REFRACTIVE INDEX (REFREL) FOR A GIVEN SPHERE REFRACTIVE
C       INDEX, MEDIUM REFRACTIVE INDEX, RADIUS, AND FREE SPACE
C       WAVELENGTH.  IT THEN CALLS BHMIE, THE SUBROUTINE THAT COMPUTES
C       AMPLITUDE SCATTERING MATRIX ELEMENTS AND EFFICIENCIES
C       ---------------------------------------------------------------
      PARAMETER (MXA=200,MXA2=400)
      COMPLEX REFREL,S1(MXA2),S2(MXA2)
      LOGICAL PRNT
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
        IF( PRNT ) WRITE(6,75) ANG,S11,POL,S33,S34
  355 CONTINUE
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
     +                  QEXT,QSCA,QBACK)
      PARAMETER (MXA=200,MXA2=400)
      DIMENSION AMU(MXA),THETA(MXA),PI(MXA),TAU(MXA),PI0(MXA),PI1(MXA)
      COMPLEX D(10000),Y,REFREL,XI,XI0,XI1,AN,BN,S1(MXA2),S2(MXA2)
      DOUBLE PRECISION PSI0,PSI1,PSI,DN,DX
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
