/*	fall_speed[size][layer]=*/
/*	  vt_(*/
/*	      crystal_length+size,*/
/*	      prism_radius+size,*/
/*	      crystal_mass+size,*/
/*	      equiv_rad_squared+size,*/
/*	      env_pressure+layer,*/
/*	      env_temperature+layer,*/
/*	      env_density+layer,*/
/*	      );*/
C
C     Fall speed routines from Heymsfield's microphysics program. 
C     Edited and sent by Joanne Parish
C
Cwhat you need:
CDIAM diameter (cm -- I think)
CTHICK thickness 
CM     habit
C      1 -- DENDRITES
C      2 -- PLANAR AGGREGATES (NOT CURRENTLY IN USE)
C      3 -- SPATIAL AGGREGATES
C      4 -- GRAUPEL
C      5 -- GRAUPEL-2 (high
C      6 -- WATER DROPS
C      7 -- FROZEN DROPS
C      8 -- COLUMNS
C      9 -- BULLET ROSETTES
C      for 8 and 9, DIAM IS LENGTH (LARGEST DIAM), THICK IS ICE DIAMETER
CAMASS mass
CAREA  cross-sectional area
CANU   KINEMATIC VISCOSITY
C      ANU=(.043/PS*ABT**2.5)/(ABT+120.)
CPS    atmospheric pressure (bars)
CABT   absolute temperature
CRHOF  AIR DENSITY                                                               
C      RHOF=348.38*PS/ABT                                                        
      SUBROUTINE VT(crys_length,prism_rad,crys_mass,
     &     equiv_r2,env_pres,env_temp,env_dens)
C
C     THIS SUBROUTINE CALLS XRE TO GET THE CURRENT REYNOLDS NUMBER,
C        THEN CALCULATES TERMINAL VELOCITY.
C
cz      COMMON /PARAM/ PS,RFT,ABT,ANU,RHOF
cz      COMMON /TOTDIM/ DIAM,VT0,THICK,AMASS,VOL,AREA,SHAPE,TH,M,ISPHER
cz      COMMON /XRE/ RE,X,F,KTYP
C
      real crys_length,prism_rad,crys_mass,
     &     equiv_r2,env_pres,env_temp,env_dens

      real diam,thick,amass,area,ps,abt,rhof,
cconvert the inputs to the correct units
      diam=crys_length*100.     !microns --> cm
      thick=2.*prism_rad.*100.  !microns --> cm
      m=8                       !for columns 
      amass                     !kg --> g
      area                      !m^2 --> cm^2
      ps,                       !Pa --> bars
      abt,                
      rhof,                     !kg/m^3 --> ??

      ANU=(.043/PS*ABT**2.5)/(ABT+120.)
      CALL XRE(0,M,DIAM,THICK,AMASS,AREA,X,RE,F,KTYP)
      IF(DIAM.LE.0.0)THEN
          VT0=0.0 
          GOTO 100
      END IF
      IF (M.GE.8) THEN
         VT0=ANU*RE/THICK
      ELSE
         VT0=ANU*RE/DIAM
      END IF
  100 RETURN
      END
      SUBROUTINE XRE(IFLAG,M,DIAM,THICK,AMASS,AREA,X,RE,F,KTYP)
C
C     THIS SUBROUTINE CALCULATES X AND RE FOR DETERMINING
C        TERMINAL VELOCITY.
C
      COMMON /PARAM/ PS,RFT,ABT,ANU,RHOF
      DIMENSION IKR(9)
      DATA IKR /1,2,2,3,4,5,2,6,6/
      DATA RMS /1.414214/
      DATA PI /3.1415926/
C
      KTYP=1
      X=2495.6*AMASS/(RHOF*1.E-6*ANU**2)
      IF(X.LT.0.00001)X=0.00001
C
         IR=IKR(M)
         IF (M.EQ.4 .AND. X.LT.1138.) IR=5
         IF (M.EQ.5 .AND. X.LT.1580.) IR=5
C
      GO TO (10,20,30,40,50,60) IR
C
10    CONTINUE
C
C        DENDRITES
C        THIS IS CHANGE IN RE ACCORDING TO JAYAWEERA
         AX=ALOG10(X)
         IF (X.GT.2.2E4) THEN
            RE=0.5*X**.558
         ELSE IF (X.GT.5.) THEN
            RE=10.**(-1.17758+0.84874*AX+0.02248*AX**2-0.00998*AX**3)
         ELSE
C
C           ACCOUNT FOR SMALL PLATES   DY1 W/L=.1   DY2 W/L=.2
            DY1=-1.140514+1.021225*AX-.09244896*AX**2-.07297163*AX**3
            DY2=-1.297247+.8369597*AX-.1217281*AX**2+.3278072*AX**3
            RE1=10.**DY1
            RE2=10.**DY2
            E=DIAM/THICK
            IF (E.LT.2.) THEN
               RE=RE2
            ELSE IF (E.LE.10.) THEN
               RE=RE1+(RE2-RE1)*(1.80886-.5721/SQRT(1./E))
            ELSE
               RE=RE1
            END IF
         END IF
         F=ALOG(RE)
         RETURN
C
20    CONTINUE
C
C           SPHERES AND AGGREGATES(SMOOTH SPHERE CD-RE)
         IF(AMASS.LT.0.004)THEN
            AX=ALOG(X)
            F=-3.18657+.992696*AX-.153193E-2*AX**2-.987059E-3*AX**3
     *        -.578878E-3*AX**4+.855176E-4*AX**5-.327815E-5*AX**6
            IF(DIAM.LT.0.00001)DIAM=0.00001
            CSC=1.+2.51*RL/DIAM
            RE=CSC*EXP(F)
         ELSE IF(AMASS.LT. 0.40)THEN
           AX=ALOG10(X)
           F=-1.81391+1.34671*AX-0.12427*AX**2+0.006344*AX**3
           RE=10.0**F
         ELSE
           AX=ALOG10(X)
           F=5.33283-1.21728*AX+0.19007*AX**2-0.007005*AX**3
           RE=10.0**F
         END IF
C
C        CHECK FOR MINIMUM DRAG COEFFICIENT OF 0.6
         CD=X/(RE*RE)
         IF (CD.LT.0.6) RE=SQRT(X/0.6)
         GO TO 200
C
30    CONTINUE
C
C           GRAUPEL (KNIGHT AND HEYMSFIELD CD-RE)
            RE=0.4487*X**0.5536
            IF (X.LE.1.75E3) THEN
               AX=ALOG10(X)
               F=-1.81391+1.34671*AX-0.12427*AX**2+0.006344*AX**3
               RE=10.0**F
            END IF
         GO TO 200
C
40    CONTINUE
C
C          MATSON AND HUGGINS CD-RE
C           GRAUPEL-2
            RE=0.4742*X**0.5461
         GO TO 200
C
50    CONTINUE
C
C***********************************************************************
C
C         BEARD CD-RE FOR WATER DROPS
C
        I300=0
        IF(DIAM.LT. 0.1)THEN
            AX=ALOG(X)
            Y=-3.18657+0.992696*AX-0.00153193*AX**2-0.000987059*AX**3
     1      -0.578878E-3*AX**4+0.855176E-4*AX**5-0.327815E-5*AX**6
            RE=EXP(Y)
            I300=1
         END IF
         IF (I300.EQ.0 .OR. RE.GT.300.) THEN
               ANP=484.23/(ANU**4*RHOF**2*1.E-12)
               ANP6=ANP**(1./6.)
               X=ALOG(16.752*DIAM**2*ANP6)
               Y=-5.00015+5.23778*X-2.04914*X**2+.475294*X**3
     *             -.0542819*X**4+.00238449*X**5
               RE=ANP6*EXP(Y)
         END IF
        IF (DIAM.GE.1. .AND. IFLAG.EQ.0) PRINT 700,DIAM,RE
  700 FORMAT(' ','WATER DROP TOO LARGE, DIAM=',F10.4,'RE=',E12.4)
         GO TO 200
C
60    CONTINUE
C
C***************************************************************************
C
C           COLUMNS AND BULLET ROSETTES
            AREAX=AREA
            IF (M.EQ.9) AREAX=AREA*RMS/2.
            X=2.*AMASS*THICK*THICK*980./(RHOF*1.E-6*ANU**2*AREAX)
            AX=ALOG10(X)
            E=THICK/DIAM
            IF (X.GT.25.) THEN
               RE=10.**(-1.10114+1.05687*AX-.09244*AX*AX+.00535*AX**3)
               DEL=10.**(-.90186+1.0034*AX-.10142*AX*AX+.0083*AX**3)-RE
               RE=RE+DEL*(1.81-2.56*SQRT(E))
            ELSE
C
C              ACCOUNT FOR SMALL COLUMNS -- CY1  W/L=.5     CY2  W/L=.1
              CY1=-1.141252+.9574312*AX-.03779655*AX*AX+.006800489*AX**3
               CY2=-.8463266+.9076694*AX-.0631194*AX*AX-.00392794*AX**3
               RE=10.**CY1
               IF (E.LE.0.5) RE=RE+(10.**CY2-RE)*(1.809-2.56*SQRT(E))
            END IF
         GO TO 200
C
C
200   CONTINUE
      IF (RE.LE.0.) THEN
         PRINT 101, X,RE,M,IR
         RE=1.E-4
      END IF
      IF (RE.GT.1.E6) PRINT 102, RE,X,M,IR
      F=ALOG(RE)
      IF (RE.GT.400.) KTYP=2
      RETURN
 100   FORMAT (10X,*RE NEGATIVE -- X,RE,A,B=*,4E11.3)
101   FORMAT (10X,'RE NEGATIVE -- X,RE,M,IR=',2E11.3,2I5)
102   FORMAT (10X,'RE TOO HIGH -- X,RE,M,IR=',2E11.3,2I5)
      END


