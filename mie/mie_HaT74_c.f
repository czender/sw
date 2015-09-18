C********1*********2*********3*********4*********5*********6*********7**
      subroutine mie_c(NSIZEPZ,NGZ,IPARTZ,NUMANGZ,STEPSZ,ANGLUSZ,
     &NUMSDZ,NSDZ,AZ,BZ,CZ,DZ,EZ,NRZ,NIZ,ALAMZ,R1Z,R2Z)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 NR,NI,NUMPAR(10),KSCA(10),KEXT(10)
      REAL*8 NRZ,NIZ
      COMMON R1,R2,NR,NI,ALAM,PI,A(10),B(10),C(10),D(10),G(10),REFF(10),
     &       VEFF(10),SEFF(10),NUMPAR(10),KSCA(10),KEXT(10),PIZERO(10),
     &       COSBAR(10),REXT(10),RSCA(10),RBAR(10),RR(96),WW(96),E(10),
     &       R(2000),W(2000),SIZDIS(10,2000),TD(160),CT(80),S2(80),
     &       P1(10,160),P2(10,160),P3(10,160),P4(10,160),P1F(80),
     &       P2F(80),P3F(80),P4F(80),P1B(80),P2B(80),P3B(80),P4B(80),
     &       NSD(10),JX,NGAUSS,IPART,NGT,NUMSD
      DIMENSION STEPS(20),ANGLUS(20)
      OPEN(UNIT=6,STATUS='UNKNOWN',FILE='mie.out')

      NSIZEP=NSIZEPZ
      NG=NGZ
      IPART=IPARTZ
      NUMANG=NUMANGZ
      STEPS(1)=STEPSZ
      ANGLUS(1)=ANGLUSZ
      NUMSD=NUMSDZ
      NSD(1)=NSDZ
      A(1)=AZ
      B(1)=BZ
      C(1)=CZ
      D(1)=DZ
      E(1)=EZ
      NR=NRZ
      NI=NIZ
      WRITE (6,*) 'nrz=',nrz,'nr=',nr,'alam=',alam
      ALAM=ALAMZ
      R1=R1Z
      R2=R2Z

      PI = 3.1415926535897932D0
      TWOPI = PI * 2.D0
      WRITE (6,10) NSIZEP,NG,IPART
   10 FORMAT(' NSIZEP=',I2,' NG=',I2,' IPART=',I3)
      NGT = NG * IPART
C********1*********2*********3*********4*********5*********6*********7**
C COMPUTE SCATTERING ANGLES
C********1*********2*********3*********4*********5*********6*********7**
      ANGLEL = 0.D0
      WRITE (6,20) NUMANG
   20 FORMAT(' THE NUMBER OF INTERVALS USED IN SCATTERING ANGLE DETER
     &MINATION IS NUMANG =',I3)
      K = 1
   30 FORMAT(' INTERVAL NO.',I3,2X,'IS FROM',F7.2,2X,'TO',F7.2,2X,
     &'IN STEPS OF',F7.3)
      DO 50 J = 1,NUMANG
         STEP = STEPS(J) 
         ANGLEU = ANGLUS(J)
         WRITE (6,30) J,ANGLEL,ANGLEU,STEP
         INDEX = (((ANGLEU - ANGLEL) + 10.D-8) / STEP) + 1
         K = K - 1
         DO 40 I = 1,INDEX
            K = K + 1
            TD(K) = ANGLEL + (I - 1) * STEP
   40    CONTINUE
         ANGLEL = ANGLEU
   50 CONTINUE
      JX = K
      JXF = 2 * JX - 1
      DO 60 J = 1,JX
         CS = DCOS(PI * TD(J) / 180.D0)
         CT(J) = CS
         S2(J) = 1.D0 - CS * CS
         L = JXF - (J - 1)
         TD(L) = 180.D0 - TD(J)
   60 CONTINUE
C********1*********2*********3*********4*********5*********6*********7**
C READ SIZE DISTRIBUTION DATA
C********1*********2*********3*********4*********5*********6*********7**
      WRITE (6,*) NUMSD
   70 FORMAT(' NUMSD=',I4)
      DO 80 N = 1,NUMSD
         WRITE (6,90) NSD(N),A(N),B(N),C(N),D(N),E(N)
   80 CONTINUE
   90 FORMAT(' NSD=',I3,'  A=',E15.8,'  B=',E15.8,'  C=',E15.8,
     &' D=',E15.8,'  E=',E15.8)
C********1*********2*********3*********4*********5*********6*********7**
C READ OPTICAL CONSTANTS, ETC.
C********1*********2*********3*********4*********5*********6*********7**
C  100 READ (5,*,END = 290) NR,NI,ALAM,R1,R2
C  100 CONTINUE
      IF(NSIZEP .EQ. 7) ALAM = TWOPI
      WRITE (6,110) NR,NI,ALAM,R1,R2
  110 FORMAT(' NR=',F12.5,' NI=',F12.5,' ALAM=',F10.5,
     &' R1=',F12.5,' R2=',F12.5)
      COEFIN = ALAM * ALAM / PI
C********1*********2*********3*********4*********5*********6*********7**
C OBTAIN INTEGRATION POINTS AND WEIGHTS
C********1*********2*********3*********4*********5*********6*********7**
      RNGE = R2 - R1
      DIV = RNGE / IPART
      DO 120 I = 1,IPART
         XBOT = (I - 1) * DIV + R1
         XTOP = XBOT + DIV
         LOWER = (I - 1) * NG
         CALL GAUSST(NG,XBOT,XTOP,RR,WW)
         DO 125 KK = 1,NG
            KKK = LOWER + KK
            R(KKK) = RR(KK)
            W(KKK) = WW(KK)
  125    CONTINUE
  120 CONTINUE
C********1*********2*********3*********4*********5*********6*********7**
C COMPUTE SIZE DISTRIBUTION AT INTEGRATION POINTS
C********1*********2*********3*********4*********5*********6*********7**
      CALL SIZEDS
C********1*********2*********3*********4*********5*********6*********7**
C INITIALIZATIONS
C********1*********2*********3*********4*********5*********6*********7**
      DO 130 N = 1,NUMSD
         NUMPAR(N) = 0.D0
         KSCA(N) = 0.D0
         KEXT(N) = 0.D0
         COSBAR(N) = 0.D0
         RSCA(N) = 0.D0
         REXT(N) = 0.D0
         RBAR(N) = 0.D0
         REFF(N) = 0.D0
         VEFF(N) = 0.D0
         SEFF(N) = 0.D0
         G(N) = 0.D0
         DO 135 I = 1,JXF
            P1(N,I) = 0.D0
            P2(N,I) = 0.D0
            P3(N,I) = 0.D0
            P4(N,I) = 0.D0
  135    CONTINUE
  130 CONTINUE
      JCO = 0
C********1*********2*********3*********4*********5*********6*********7**
C JCO LOOP (NEW PARTICLE) BEGINS HERE
C********1*********2*********3*********4*********5*********6*********7**
  140 JCO = JCO + 1
      X = R(JCO) * TWOPI / ALAM
      CALL SINPAR(X,QEXT,QSCA,COSB)
      DO 170 N = 1,NUMSD
         SW = SIZDIS(N,JCO) * W(JCO)
         SWPI = SW * PI
         DO 150 I = 1,JX
            P1(N,I) = P1F(I) * SW + P1(N,I)
            P2(N,I) = P2F(I) * SW + P2(N,I)
            P3(N,I) = P3F(I) * SW + P3(N,I)
            P4(N,I) = P4F(I) * SW + P4(N,I)
  150    CONTINUE
         JXX = JX + 1
         JXF = 2 * JX - 1
         KK = 0
         DO 160 I = JXX,JXF
            KK = KK + 1
            K = I - 2 * KK
            P1(N,I) = P1B(K) * SW + P1(N,I)
            P2(N,I) = P2B(K) * SW + P2(N,I)
            P3(N,I) = P3B(K) * SW + P3(N,I)
            P4(N,I) = P4B(K) * SW + P4(N,I)
  160    CONTINUE
         Z1 = R(JCO)
         Z2 = Z1 * Z1
         Z3 = Z2 * Z1
         PIR2QS = Z2 * QSCA * SWPI
         KSCA(N) = KSCA(N) + PIR2QS
         KEXT(N) = KEXT(N) + Z2 * QEXT * SWPI
         RSCA(N) = RSCA(N) + Z1 * PIR2QS
         REXT(N) = REXT(N) + Z3 * QEXT * SWPI
         NUMPAR(N) = NUMPAR(N) + SW
         RBAR(N) = RBAR(N) + Z1 * SW
         REFF(N) = REFF(N) + Z3 * SWPI
         G(N) = G(N) + Z2 * SWPI
         COSBAR(N) = COSBAR(N) + COSB * PIR2QS
  170 CONTINUE
      IF(JCO .NE. NGT) GO TO 140
C********1*********2*********3*********4*********5*********6*********7**
C JCO LOOP (NEW PARTICLE) ENDS HERE
C********1*********2*********3*********4*********5*********6*********7**
      DO 260 N = 1,NUMSD
         COSBAR(N) = COSBAR(N) / KSCA(N)
         PIZERO(N) = KSCA(N) / KEXT(N)
         RSCA(N) = RSCA(N) / KSCA(N)
         REXT(N) = REXT(N) / KEXT(N)
         RBAR(N) = RBAR(N) / NUMPAR(N)
         REFF(N) = REFF(N) / G(N)
         RF = REFF(N)
         Z2 = 0.D0
         Z3 = 0.D0
         DO 180 JCO = 1,NGT
            RMRF = R(JCO) - RF
            RM2 = RMRF * RMRF * SIZDIS(N,JCO) * W(JCO) * R(JCO) * R(JCO)
            RM3 = RM2 * RMRF
            Z2 = Z2 + RM2
            Z3 = Z3 + RM3
  180    CONTINUE
         VEFF(N) = Z2 * PI / G(N) / RF / RF
         SEFF(N) = Z3 * PI / G(N) / RF / RF / RF / VEFF(N) ** 1.5
         DO 190 I = 1,JXF
            CS = COEFIN / KSCA(N)
            P1CS = P1(N,I) * CS
            P2CS = P2(N,I) * CS
            P1(N,I) = 0.5D0 * (P1CS + P2CS)
            P2(N,I) = 0.5D0 * (P1CS - P2CS)
            P3(N,I) = P3(N,I) * CS
            P4(N,I) = P4(N,I) * CS
  190    CONTINUE
C********1*********2*********3*********4*********5*********6*********7**
C PRINTED OUTPUT
C********1*********2*********3*********4*********5*********6*********7**
         WRITE (6,200) NSD(N),A(N),B(N),C(N),D(N),E(N)
  200    FORMAT(' FOLLOWING FOR NSD=',I3,' A=',E15.8,' B=',E15.8,
     &   ' C=',E15.8,' D=',E15.8,' E=',E15.8)
         WRITE (6,210)
  210    FORMAT(' SCAT. ANGLE',3X,'P11',12X,'P21',12X,'P33',
     &   12X,'P43',6X,'POLARIZATION',/)
         DO 220 I = 1,JXF
            XP = -100.D0 * P2(N,I) / P1(N,I)
            WRITE (6,230) TD(I),P1(N,I),P2(N,I),P3(N,I),P4(N,I),XP
  220    CONTINUE
  230    FORMAT(1X,F8.3,4D15.8,D12.5)
         WRITE (6,240) REXT(N),RSCA(N),RBAR(N),REFF(N)
  240    FORMAT(' REXT='E12.5,' RSCA=',E12.5,' RBAR=',E12.5,' REFF=',
     &   E12.5)
         WRITE (6,250) VEFF(N),SEFF(N),G(N)
  250    FORMAT(' VEFF=',E12.5,' SEFF=',E12.5,'   G=',E12.5)
         QSCA = KSCA(N) / G(N)
         QEXT = KEXT(N) / G(N)
         WRITE (6,270) COSBAR(N),PIZERO(N),NUMPAR(N)
         WRITE (6,280) KSCA(N),KEXT(N),QSCA,QEXT
  260 CONTINUE
  270 FORMAT(' COSBAR=',F12.6,' PIZERO=',F12.6,' NUMPAR=',F12.6)
  280 FORMAT(' KSCA=',F12.6,' KEXT=',F12.6,' QSCA=',F12.6,
     &' QEXT=',F12.6)
C      GO TO 100
  290 STOP
      END
C********1*********2*********3*********4*********5*********6*********7**
C********1*********2*********3*********4*********5*********6*********7**
      SUBROUTINE SINPAR(X,QEXT,QSCA,COSB)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 NR,NI,NUMPAR(10),KSCA(10),KEXT(10)
      COMPLEX*16 WN,WN1,WN2,SMAN,SMAN1,SMBN,SMBN1,NC,Z
      COMPLEX*16 RRF,RRFX,TC1,TC2,CDUM1,CDUM2
      DIMENSION TAUN(80),TAUN1(80),TAUN2(80)
      DIMENSION PIX1(80),PIX2(80),PIX3(80),Z(7000)
      COMMON R1,R2,NR,NI,ALAM,PI,A(10),B(10),C(10),D(10),G(10),REFF(10),
     &       VEFF(10),SEFF(10),NUMPAR(10),KSCA(10),KEXT(10),PIZERO(10),
     &       COSBAR(10),REXT(10),RSCA(10),RBAR(10),RR(96),WW(96),E(10),
     &       R(2000),W(2000),SIZDIS(10,2000),TD(160),CT(80),S2(80),
     &       P1(10,160),P2(10,160),P3(10,160),P4(10,160),P1F(80),
     &       P2F(80),P3F(80),P4F(80),P1B(80),P2B(80),P3B(80),P4B(80),
     &       NSD(10),JX,NGAUSS,IPART,NGT,NUMSD
      EQUIVALENCE (WN,WNR),(WN1,WN1R)
   10 NC = DCMPLX(NR,-NI)
      RRF = 1.D0 / NC
      RX = 1.D0 / X
      RRFX = RRF * RX
      T1 = DSQRT(X * X * (NR * NR + NI * NI))
      NMX1 = 1.1D0 * T1
      IF(NMX1 .LE. 6999) GO TO 30
      WRITE (6,20) NMX1
   20 FORMAT(' DIMENSION OF Z AND LIMIT FOR NMX1 SHOULD BE INCREASED,',
     & ' NMX1 =',I10)
      STOP
   30 NMX2 = T1
      IF(NMX1 .GT. 150) GO TO 40
      NMX1 = 150
      NMX2 = 135
   40 Z(NMX1 + 1) = (0.D0,0.D0)
      DO 50 N = 1,NMX1
         NN = NMX1 - N + 1
         Z(NN) = (NN + 1) * RRFX - 1.D0 / ((NN + 1) * RRFX + Z(NN+1))
   50 CONTINUE
      DO 60 J = 1,JX
         PIX1(J) = 0.D0
         PIX2(J) = 1.D0
         TAUN2(J) = 0.D0
         TAUN1(J) = CT(J)
   60 CONTINUE
      T1 = DCOS(X)
      T2 = DSIN(X)
      WN2 = DCMPLX(T1,-T2)
      WN1 = DCMPLX(T2,T1)
      WN = RX * WN1 - WN2
      TC1 = Z(1) * RRF + RX
      TC2 = Z(1) * NC + RX
      SMAN = (TC1 * WNR - WN1R) / (TC1 * WN - WN1)
      SMBN = (TC2 * WNR - WN1R) / (TC2 * WN - WN1)
      SMANI = DIMAG(SMAN)
      SMBNI = DIMAG(SMBN)
      SMANR = DREAL(SMAN)
      SMBNR = DREAL(SMBN)
      IF(DABS(SMANR) .LT. 1.D-30) SMANR = 0.D0
      IF(DABS(SMANI) .LT. 1.D-30) SMANI = 0.D0
      IF(DABS(SMBNR) .LT. 1.D-30) SMBNR = 0.D0
      IF(DABS(SMBNI) .LT. 1.D-30) SMBNI = 0.D0
      SMBN1 = SMBN
      SMAN1 = SMAN
      SMAN1I = DIMAG(SMAN1)
      SMBN1I = DIMAG(SMBN1)
      SMAN1R = DREAL(SMAN1)
      SMBN1R = DREAL(SMBN1)
      SMANR = 1.5 * SMANR
      SMANI = 1.5 * SMANI
      SMBNR = 1.5 * SMBNR
      SMBNI = 1.5 * SMBNI
      DO 70 J = 1,JX
         P = PIX2(J)
         T = TAUN1(J)
         SMANRP = SMANR * P
         SMANIP = SMANI * P
         SMBNRP = SMBNR * P
         SMBNIP = SMBNI * P
         SMANRT = SMANR * T
         SMANIT = SMANI * T
         SMBNRT = SMBNR * T
         SMBNIT = SMBNI * T
         P1F(J) = SMANRP + SMBNRT
         P2F(J) = SMANIP + SMBNIT
         P3F(J) = SMBNRP + SMANRT
         P4F(J) = SMBNIP + SMANIT
         P1B(J) = SMANRP - SMBNRT
         P2B(J) = SMANIP - SMBNIT
         P3B(J) = SMBNRP - SMANRT
         P4B(J) = SMBNIP - SMANIT
   70 CONTINUE
      QEXT = 2.D0 * (SMANR + SMBNR)
      QSCA = (SMANR * SMANR + SMANI * SMANI + SMBNR * SMBNR +
     &        SMBNI * SMBNI) / 0.75D0
      COSBQS = 0.D0
      N = 2
C********1*********2*********3*********4*********5*********6*********7**
C MAJOR LOOP BEGINS HERE
C********1*********2*********3*********4*********5*********6*********7**
   80 T1 = 2 * N - 1
      T3 = T1 + 2
      WN2 = WN1
      WN1 = WN
      WN = T1 * RX * WN1 - WN2
      CDUM1 = Z(N)
      CDUM2 = N * RX
      TC1 = CDUM1 * RRF + CDUM2
      TC2 = CDUM1 * NC + CDUM2
      SMAN = (TC1 * WNR - WN1R) / (TC1 * WN - WN1)
      SMBN = (TC2 * WNR - WN1R) / (TC2 * WN - WN1)
      SMANI = DIMAG(SMAN)
      SMBNI = DIMAG(SMBN)
      SMANR = DREAL(SMAN)
      SMBNR = DREAL(SMBN)
      IF(DABS(SMANR) .LT. 1.D-30) SMANR = 0.D0
      IF(DABS(SMANI) .LT. 1.D-30) SMANI = 0.D0
      IF(DABS(SMBNR) .LT. 1.D-30) SMBNR = 0.D0
      IF(DABS(SMBNI) .LT. 1.D-30) SMBNI = 0.D0
      QEXT = QEXT + T3 * (SMANR + SMBNR)
      TX = SMANR * SMANR + SMANI * SMANI + SMBNR * SMBNR +
     &     SMBNI * SMBNI
      QSCA = QSCA + T3 * TX
      T2 = N - 1
      DO 90 J = 1,JX
         T1PIX2 = T1 * PIX2(J)
         CTJ = CT(J)
         PIX1J = PIX1(J)
         PIX3(J) = (T1PIX2 * CTJ - N * PIX1J) / T2
         TAUN(J) = CTJ * (PIX3(J) - PIX1J) - T1PIX2 * S2(J) + TAUN2(J)
   90 CONTINUE
      T5 = N
      T4 = T1 / (T5 * T2)
      T2 = (T2 * (T5 + 1.D0)) / T5
      COSBQS = COSBQS + T2 * (SMAN1R * SMANR + SMAN1I * SMANI +
     &         SMBN1R * SMBNR + SMBN1I * SMBNI) + T4 * (SMAN1R *
     &         SMBN1R + SMAN1I * SMBN1I)
      T2 = N * (N + 1)
      T1 = T3 / T2
      K = (N / 2) * 2
      DO 110 J = 1,JX
         P = T1 * PIX3(J)
         T = T1 * TAUN(J)
         SMANRP = SMANR * P
         SMANIP = SMANI * P
         SMBNRP = SMBNR * P
         SMBNIP = SMBNI * P
         SMANRT = SMANR * T
         SMANIT = SMANI * T
         SMBNRT = SMBNR * T
         SMBNIT = SMBNI * T
         P1F(J) = P1F(J) + SMANRP + SMBNRT
         P2F(J) = P2F(J) + SMANIP + SMBNIT
         P3F(J) = P3F(J) + SMBNRP + SMANRT
         P4F(J) = P4F(J) + SMBNIP + SMANIT
         IF(K .EQ. N) GO TO 100
         P1B(J) = P1B(J) + SMANRP - SMBNRT
         P2B(J) = P2B(J) + SMANIP - SMBNIT
         P3B(J) = P3B(J) + SMBNRP - SMANRT
         P4B(J) = P4B(J) + SMBNIP - SMANIT
         GO TO 110
  100    P1B(J) = P1B(J) - SMANRP + SMBNRT
         P2B(J) = P2B(J) - SMANIP + SMBNIT
         P3B(J) = P3B(J) - SMBNRP + SMANRT
         P4B(J) = P4B(J) - SMBNIP + SMANIT
  110 CONTINUE
      IF(TX .LT. 1.D-14) GO TO 130
      DO 120 J = 1,JX
         PIX1(J) = PIX2(J)
         PIX2(J) = PIX3(J)
         TAUN2(J) = TAUN1(J)
         TAUN1(J) = TAUN(J)
  120 CONTINUE
      SMBN1 = SMBN
      SMAN1 = SMAN
      SMAN1I = DIMAG(SMAN1)
      SMBN1I = DIMAG(SMBN1)
      SMAN1R = DREAL(SMAN1)
      SMBN1R = DREAL(SMBN1)
      N = N + 1
      IF(N .LE. NMX2) GO TO 80
C********1*********2*********3*********4*********5*********6*********7**
C MAJOR LOOP ENDS HERE
C********1*********2*********3*********4*********5*********6*********7**
      WRITE (6,20)
      STOP
  130 DO 140 J = 1,JX
         T1 = P1F(J)
         T2 = P2F(J)
         T3 = P3F(J)
         T4 = P4F(J)
         P1F(J) = T3 * T3 + T4 * T4
         P2F(J) = T1 * T1 + T2 * T2
         P3F(J) = T1 * T3 + T2 * T4
         P4F(J) = T2 * T3 - T4 * T1
         T1 = P1B(J)
         T2 = P2B(J)
         T3 = P3B(J)
         T4 = P4B(J)
         P1B(J) = T3 * T3 + T4 * T4
         P2B(J) = T1 * T1 + T2 * T2
         P3B(J) = T1 * T3 + T2 * T4
         P4B(J) = T2 * T3 - T4 * T1
  140 CONTINUE
      T1 = 2.D0 * RX * RX
      QEXT = QEXT * T1
      QSCA = QSCA * T1
      COSB = 2.D0 * COSBQS * T1 / QSCA
      RETURN
      END
C********1*********2*********3*********4*********5*********6*********7**
C********1*********2*********3*********4*********5*********6*********7**
      SUBROUTINE SIZEDS
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 NR,NI,NUMPAR(10),KSCA(10),KEXT(10)
      COMMON R1,R2,NR,NI,ALAM,PI,A(10),B(10),C(10),D(10),G(10),REFF(10),
     &       VEFF(10),SEFF(10),NUMPAR(10),KSCA(10),KEXT(10),PIZERO(10),
     &       COSBAR(10),REXT(10),RSCA(10),RBAR(10),RR(96),WW(96),E(10),
     &       R(2000),W(2000),SIZDIS(10,2000),TD(160),CT(80),S2(80),
     &       P1(10,160),P2(10,160),P3(10,160),P4(10,160),P1F(80),
     &       P2F(80),P3F(80),P4F(80),P1B(80),P2B(80),P3B(80),P4B(80),
     &       NSD(10),JX,NGAUSS,IPART,NGT,NUMSD
      DO 190 N = 1,NUMSD
         NSDN = NSD(N)
         GO TO (10,50,70,100,120,115,111,112),NSDN
C********1*********2*********3*********4*********5*********6*********7**
C THE FOLLOWING (NSD=1) IS THE TWO PARAMETER GAMMA DISTRIBUTION (2.56)
C (STANDARD)
C********1*********2*********3*********4*********5*********6*********7**
   10    RB = 1.D0 / B(N)
         IF(RB .GT. 2.D0) GO TO 30
         WRITE (6,20) B(N)
   20    FORMAT(' B FOR THE GAMMA DISTRIBUTION MUST BE LESS THAN',
     &   ' 0.5, BUT HAS BEEN READ IN AS',1PD11.4)
         STOP
   30    AB = RB / A(N)
         GM = RB - 2.D0
         DLAB = GM * DLOG(AB)
         DLGM = DLGAMA(GM)
         DO 40 I =1,NGT
            SIZDIS(N,I) = DEXP(DLAB + (GM - 1.D0) *
     &      DLOG(R(I)) - AB * R(I) - DLGM)
   40    CONTINUE
         GO TO 190
C********1*********2*********3*********4*********5*********6*********7**
C THE FOLLOWING (NSD=2) IS THE THREE PARAMETER GAMMA DISTRIBUTION
C (MODIFIED) DEIRMENDJIAN
C********1*********2*********3*********4*********5*********6*********7**
   50    RAC = (A(N) + 1.D0) / C(N)
         DLGN = DLGAMA(RAC)
         DLB = RAC * DLOG(B(N)) - DLGN
         DO 60 I = 1,NGT
            SIZDIS(N,I) = DEXP(DLB + A(N) * DLOG(R(I)) -
     &                    B(N) * R(I) ** C(N)) * C(N)
   60    CONTINUE
         GO TO 190
C********1*********2*********3*********4*********5*********6*********7**
C THE FOLLOWING (NSD=3) IS THE BIMODAL GAMMA DISTRIBUTION (2.59,C=A2)
C********1*********2*********3*********4*********5*********6*********7**
   70    RB = 1.D0 / B(N)
         IF(RB .GT. 2.D0) GO TO 80
         WRITE (6,20) B(N)
         STOP
   80    A1B = RB / A(N)
         A2B = RB / C(N)
         GM = RB - 2.D0
         DLA1B = GM * DLOG(A1B)
         DLA2B = GM * DLOG(A2B)
         DLGM = DLGAMA(GM)
         DO 90 I = 1,NGT
            D1 = DEXP(DLA1B + (GM - 1.D0) * DLOG(R(I)) -
     &           A1B * R(I) - DLGM)
            D2 = DEXP(DLA2B + (GM - 1.D0) * DLOG(R(I)) -
     &           A2B * R(I) - DLGM)
            SIZDIS(N,I) = 0.5D0 * D1 + 0.5D0 * D2
   90    CONTINUE
         GO TO 190
C********1*********2*********3*********4*********5*********6*********7**
CTHE FOLLOWING (NSD=4) IS THE LOG-NORMAL DISTRIBUTION(2.60,A=RG,B=SIGMA)
C********1*********2*********3*********4*********5*********6*********7**
  100    TEMP = DSQRT(2.D0 * PI) * B(N)
         A2 = DLOG(A(N))
         A3 = 2.D0 * B(N) * B(N)
         DO 110 I = 1,NGT
            A1 = DLOG(R(I))
            SIZDIS(N,I) = DEXP(-(A1 - A2) ** 2 / A3) / TEMP / R(I)
  110    CONTINUE
         GO TO 190
C********1*********2*********3*********4*********5*********6*********7**
CTHE FOLLOWING (NSD=7) IS THE LOG-NORMAL DISTRIBUTION(2.60,A=RG,B=SIGMA)
Cbut not 2.60 exactly, note log sigma in denominator, etc. zender
C********1*********2*********3*********4*********5*********6*********7**
 111     TEMP = DSQRT(2.D0 * PI) * DLOG(B(N))
         WRITE (6,*) 'you are using a hack option 7: log-normal'
         A2 = DLOG(A(N))
         A3 = 2.D0 * DLOG(B(N)) * DLOG(B(N))
         DO 114 I = 1,NGT
            A1 = DLOG(R(I))
            SIZDIS(N,I) = DEXP(-(A1 - A2) ** 2 / A3) / TEMP / R(I)
 114        CONTINUE
         GO TO 190
C********1*********2*********3*********4*********5*********6*********7**
CTHE FOLLOWING (NSD=8) IS THE delta fn DISTRIBUTION
Cthat is all particles are the same size. make sure rmin and rmax 
Cenclose the radius bin desired. zender
C********1*********2*********3*********4*********5*********6*********7**
 112     WRITE (6,*) 'you are using a hack option 8: delta function'
         DO 113 I = 1,NGT
            SIZDIS(N,I) = 1./(R2-R1)
 113     CONTINUE
         GO TO 190
C********1*********2*********3*********4*********5*********6*********7**
CTHE FOLLOWING (NSD=6) IS THE BIMODAL LOG-NORMAL DISTRIBUTION
C(2.60,A=RG0,B=SIGMA0,C=RG1,D=SIGMA1,E=frac in mode 0. zender
C********1*********2*********3*********4*********5*********6*********7**
 115     TEMP = DSQRT(2.D0 * PI) * DLOG(B(N))
         TEMP2 = DSQRT(2.D0 * PI) * DLOG(D(N))
         WRITE (6,117) 
 117     FORMAT(' You are using an untested software hack',
     &   ' the bi-modal log-normal distribution')
         A2 = DLOG(A(N))
         A3 = 2.D0 * DLOG(B(N)) * DLOG(B(N))
         B2 = DLOG(C(N))
         B3 = 2.D0 * DLOG(D(N)) * DLOG(D(N))
         DO 116 I = 1,NGT
            A1 = DLOG(R(I))
            B1 = DLOG(R(I))
            SIZDIS(N,I) = E(N)*DEXP(-(A1 - A2)**2 / A3) /TEMP/R(I)
     & +(1.D0-E(N))*DEXP(-(B1 - B2) ** 2 / B3)/TEMP2/R(I)
 116     CONTINUE
         GO TO 190
C********1*********2*********3*********4*********5*********6*********7**
C THE FOLLOWING (NSD=5) IS THE POWER LAW DISTRIBUTION (2.61,
C N(R)=CONST*R**(-A) FROM RMIN=B TO RMAX=C AND ZERO OTHERWISE
C********1*********2*********3*********4*********5*********6*********7**
  120    A1 = 1.D0 - A(N)
         IF(B(N) .GT. 0.D0) GO TO 140
         WRITE (6,130)
  130    FORMAT(' RMIN=B MUST BE GREATER THAN ZERO FOR THE POWER',
     &   ' LAW - EXECUTION TERMINATED')
         STOP
  140    IF(DABS(A1) .GT. 1.D-15) GO TO 150
         CONST = DLOG(C(N) / B(N))
         GO TO 160
  150    CONST = A1 / (C(N) ** A1 - B(N) ** A1)
  160    CONTINUE
         DO 170 I = 1,NGT
            SIZDIS(N,I) = 0.D0
            IF(R(I) .LT. B(N) .OR. R(I) .GT. C(N)) GO TO 170
            SIZDIS(N,I) = CONST * R(I) ** (-A(N))
  170    CONTINUE
  180    CONTINUE
         DO 185 I = 1,NGT
            IF(SIZDIS(N,I) .LT. 1.D-30) SIZDIS(N,I) = 0.D0
  185    CONTINUE
  190 CONTINUE
      RETURN
      END
C********1*********2*********3*********4*********5*********6*********7**
C********1*********2*********3*********4*********5*********6*********7**
      SUBROUTINE GAUSST(NG,X1,X2,XP,WT)
C********1*********2*********3*********4*********5*********6*********7**
C NG=TOTAL NUMBER OF GAUSS POINTS INCLUDING BOTH SIDES
C X1=START X2=END XP=POINTS WT=WEIGHT
C********1*********2*********3*********4*********5*********6*********7**
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION XP(1),WT(1)
      DATA PI,PS,DXL/3.141592653589793D0,1.013211836423378D-01,1.D-16/
      XMID = (X2 + X1) / 2.D0
      XHAF = (X2 - X1) / 2.D0
      DNG = NG
      NN = NG / 2
      N2 = NN * 2
      IF(N2 .EQ. NG) GO TO 20
      XP(NN + 1) = XMID
      WT(NN + 1) = 1.D0
      IF(NG .LT. 2) RETURN
      PN = 1.D0
      N = 0
   10 N = N + 2
      DN = N
      DM = DN - 1.D0
      PN = PN * (DM / DN)
      IF(N. LT .N2) GO TO 10
      WT(NN + 1) = 2.D0 * XHAF / (DNG * PN) ** 2
   20 I = 0
      C = PI / DSQRT(DNG * (DNG + 1.D0) + 0.5D0 - PS) / 105.D0
   30 I = I + 1
      DI = I
      Z = PS/(4.D0 * DI - 1.D0) ** 2
      ZZ = (105.D0 + Z * (210.D0 - Z * (2170.D0 - Z * (105812.D0
     &     - 12554474.D0 * Z))))
      X = DCOS(ZZ * C * (DI - 0.25D0))
   40 N = 1
      DM = 1.D0
      PNI = 1.D0
      PNJ = X
   50 N = N + 1
      DN = N
      PNK = ((DM + DN) * X * PNJ - DM * PNI) / DN
      PNI = PNJ
      PNJ = PNK
      DM = DN
      IF(N .LT. NG) GO TO 50
      DX = PNJ * (1.D0 - X * X) / DNG / (PNI - X * PNJ)
      X = X - DX
      IF(DABS(DX) .GT. DXL) GO TO 40
      J = NG + 1 - I
      XP(I) = XMID - XHAF * X
      XP(J) = XMID + XHAF * X
      WT(I) = 2.D0 * XHAF * (1.D0 - X * X) / (DNG * PNI) ** 2
      WT(J) = WT(I)
      IF(I .LT. NN) GO TO 30
      RETURN
      END
C********1*********2*********3*********4*********5*********6*********7**
      FUNCTION DLGAMA(X)
C********1*********2*********3*********4*********5*********6*********7**
C THIS FUNCTION CALCULATES THE LOG OF THE GAMMA FUNCTION FOR USE IN
C THE MIE PARTICLE SIZE DISTRIBUTION SUBROUTINE SIZEDS
C********1*********2*********3*********4*********5*********6*********7**
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NTERMS = 26)
      DOUBLE PRECISION C(NTERMS)
      DATA C /   1.D0,       5.772156649015329D-1,-6.558780715202538D-1,
     & -4.20026350340952D-2, 1.665386113822915D-1,-4.21977345555443D-2,
     & -9.621971527877D-3,   7.218943246663D-3,   -1.1651675918591D-3,
     & -2.152416741149D-4,   1.280502823882D-4,   -2.01348547807D-5,
     & -1.2504934821D-6,     1.133027232D-6,      -2.056338417D-7,
     &  6.116095D-9,         5.0020075D-9,        -1.1812746D-9,
     &  1.043427D-10,        7.7823D-12,          -3.6968D-12,
     &  5.1D-13,            -2.06D-14,            -5.4D-15,
     &  1.4D-15,             1.D-16       /
      PARAMETER (XTOL = -30.D0)
C********1*********2*********3*********4*********5*********6*********7**
      IF(X .LE. 0.D0) THEN
        WRITE(*,'(''0'',''> ERROR : INVALID ARGUMENT IN <DLGAMA>'')') X
        STOP
      END IF
      IF(DABS(DNINT(X)-X) .LT. 1.D-3) THEN
C********1*********2*********3*********4*********5*********6*********7**
C INTEGER VALUES OF X
C********1*********2*********3*********4*********5*********6*********7**
        SUM = 0.D0
        DO 10 Y = 2.D0,X - 1.D0,1.D0
          SUM = SUM + DLOG(Y)
   10   CONTINUE
        DLGAMA = SUM
      ELSE
C********1*********2*********3*********4*********5*********6*********7**
C ARBITRARY X  -  REDUCE TO INTERVAL [0,1] IF NECESSARY
C********1*********2*********3*********4*********5*********6*********7**
        IF(X .LE. 1.D0) THEN
          Z = X
        ELSE
          N = X
          Z = X - N
        END IF
C********1*********2*********3*********4*********5*********6*********7**
C SERIES EXPANSION OVER [0,1]
C********1*********2*********3*********4*********5*********6*********7**
        SUM = 0.D0
        DO 20 K = NTERMS,1,-1
          IF(K * DLOG(Z) .GT. XTOL) THEN
            SUM = SUM + C(K) * Z ** K
          END IF
   20   CONTINUE
        DLGAMA = -DLOG(SUM)
        IF(X .GT. 1.D0) THEN
C********1*********2*********3*********4*********5*********6*********7**
C RECURSION FORMULA FOR X > 1
C********1*********2*********3*********4*********5*********6*********7**
          SUM = 0.D0
          DO 30 I = 1,N - 1
            SUM = SUM + DLOG(N - I + Z)
   30     CONTINUE
          DLGAMA = DLGAMA + SUM + DLOG(Z)
        END IF
      END IF
      RETURN
      END
C********1*********2*********3*********4*********5*********6*********7**
