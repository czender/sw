C     $Id$
      
C     Purpose: This file contains four routines to return the fraction
C     of the solar flux contained within a give wavelength interval.
C     Source code is all from Bruce P. Briegleb. There are 4 functions,
C     the first two give data from Thekeakara and Drummond (ThD71), 
C     the last two give data from Labs and Neckel (LaN68). 
      
C     SLFFTD(LMIN,LMAX): lmin,lmax are wavelengths in microns.
C     SLFWTD(WAVHGH,WAVLOW): wavhgh,wavlow are higest,lowest wavenumbers
C     SLFFLN(LMIN,LMAX): lmin,lmax are wavelengths in microns.
C     SLFWLN(WAVHGH,WAVLOW): wavhgh,wavlow are higest,lowest wavenumbers
      
C     ThD71 binned their data every 0.01 microns from 0.2 up to 1.0 micron 
C     then by 0.2 microns up to 5.0 microns for a total of 120 intervals. 
C     Valid wavelengths are between 0.15 and 4.95 microns inclusive.
      
C     LaN68 binned their data every 0.01 microns from 0.2 up to 0.6 microns 
C     then by 0.1 microns up to 5 microns for a total of 120 intervals. 
C     Valid wavelengths are between 0.2 and 5.0 microns inclusive.
      
C     Users of ThD71: TaL89, Ste78, Sli89
C     Users of LaN68: CCM2, 
C     Users of NeL84: EbC92
      
      FUNCTION SLFFTD(LMIN,LMAX)
      PARAMETER (MAXWAV=67)
      DIMENSION DLAMBD(MAXWAV)
      REAL LMIN,LMAX,LAMBDA(MAXWAV)
      
C     SOLAR FLUX FRACTIONS ACCORDING TO:
C     THEKEAKARA AND DRUMMOND (1971)
      
C     INTEGRATING VALUES FROM 0.200 MICRONS   TO  4.950 MICRONS
C     GIVES A SOLAR CONSTANT OF: 1351.3 WATTS/M2
      
C     VALUES BELOW ARE ACTUALLY THE FRACTION OF SOLAR FLUX SHORTWAVE
C     OF THE PARTICULAR WAVELENGTH.
      
      DATA LAMBDA/ .15 ,
     +     .20 , .22 , .23 , .24 , .25 , .26  ,  .27 ,
     +     .28 ,
     +     .29 , .30 , .31 , .32 , .33 , .34  ,  .35 ,
     +     .36 , .37 , .38 , .39 , .40 , .41  ,  .42 ,
     +     .43 , .44 , .45 , .46 , .47 , .48  ,  .49 ,
     +     .50 , .51 , .52 , .53 , .54 , .55  ,  .56 ,
     +     .57 ,
     +     .58 , .59 , .60 , .62 , .64 , .66  ,  .68 ,
     +     .70 , .72 , .75 , .80 , .90 ,1.00  , 1.20 ,
     +     1.40 ,1.60 ,1.80 ,2.00 ,2.20 ,2.40  , 2.60 ,
     +     2.80 ,3.00 ,3.20 ,3.40 ,3.60 ,3.80  , 4.00 ,
     +     5.00 /
      DATA DLAMBD/.0000,
     +     .0001,.0005,.0010,.0014,.0019,.0027 ,.0041 ,
     +     .0056,
     +     .0081,.0121,.0165,.0222,.0293,.0372 ,.0452 ,
     +     .0532,.0615,.0700,.0782,.0873,.0992 ,.1122 ,
     +     .1247,.1373,.1514,.1665,.1817,.1968 ,.2115 ,
     +     .2260,.2401,.2538,.2674,.2808,.2938 ,.3065 ,
     +     .3191,
     +     .3318,.3444,.3568,.3810,.4042,.4266 ,.4481 ,
     +     .4688,.4886,.5169,.5602,.6336,.6946 ,.7839 ,
     +     .8434,.8861,.9159,.9349,.9483,.9589 ,.9667 ,
     +     .9731,.9783,.9822,.9850,.9872,.9891 ,.9906 ,
     +     .9951/
      
      IF( LMIN .GE. 0.150 .AND. LMIN .LT. 4.950 ) THEN
         
C     FIND LOWER LIMIT
         N = 1
 100     IF( LMIN .GE. LAMBDA(N) ) THEN
            N = N + 1
            GO TO 100
         ENDIF
         DMIN = DLAMBD(N-1) + ((DLAMBD(N) - DLAMBD(N-1)) / (LAMBDA(N)
     1        - LAMBDA(N-1))) * ( LMIN - LAMBDA(N-1))
      ELSE IF( LMIN .LT. 0.150 ) THEN
         DMIN = 0.0
      ELSE IF( LMIN .GT. 4.950 ) THEN
         SLFFTD =  0.0
         RETURN
      ENDIF
      
C     FIND UPPER LIMIT
      IF( LMAX .GT. 0.200 .AND. LMAX .LE. 4.950 ) THEN
         N = 1
 200     IF( LMAX .GE. LAMBDA(N) ) THEN
            N = N + 1
            GO TO 200
         ENDIF
         DMAX = DLAMBD(N-1) + ((DLAMBD(N) - DLAMBD(N-1)) / (LAMBDA(N)
     1        - LAMBDA(N-1))) * ( LMAX - LAMBDA(N-1))
      ELSE IF( LMAX .GT. 4.950 ) THEN
         DMAX = 1.0000
c++csz
c     ELSE IF( LMAX .LT. 0.150 ) THEN
      ELSE IF( LMAX .LE. 0.2 ) THEN
c--csz
         SLFFTD =  0.0
         RETURN
      ENDIF
      
C     SOLAR FLUX FRACTION
      SLFFTD =  DMAX - DMIN
      RETURN
      END
      
      FUNCTION SLFWTD(WAVHGH,WAVLOW)
      
C     COMPUTES SOLAR FLUX FRACTION BETWEEN HIGHEST (WAVGHG) AND LOWEST
C     (WAVLOW) WAVENUMBER:
      
      PARAMETER (MAXWAV= 67)
      real lamdat,wavnum,solwav
      common /stddat / LAMDAT(MAXWAV),WAVNUM(MAXWAV),SOLWAV(MAXWAV)
      
C     SOLAR FLUX FRACTIONS ACCORDING TO:
C     THEKEAKARA AND DRUMMOND (1971)
      
C     INTEGRATING VALUES FROM 0.200 MICRONS   TO  4.950 MICRONS
C     GIVES A SOLAR CONSTANT OF: 1351.3 WATTS/M2
      
C     VALUES BELOW ARE ACTUALLY THE FRACTION OF SOLAR FLUX SHORTWAVE
C     OF THE PARTICULAR WAVELENGTH.
      
      DATA LAMDAT/ .15 ,
     +     .20 , .22 , .23 , .24 , .25 , .26  ,  .27 ,
     +     .28 ,
     +     .29 , .30 , .31 , .32 , .33 , .34  ,  .35 ,
     +     .36 , .37 , .38 , .39 , .40 , .41  ,  .42 ,
     +     .43 , .44 , .45 , .46 , .47 , .48  ,  .49 ,
     +     .50 , .51 , .52 , .53 , .54 , .55  ,  .56 ,
     +     .57 ,
     +     .58 , .59 , .60 , .62 , .64 , .66  ,  .68 ,
     +     .70 , .72 , .75 , .80 , .90 ,1.00  , 1.20 ,
     +     1.40 ,1.60 ,1.80 ,2.00 ,2.20 ,2.40  , 2.60 ,
     +     2.80 ,3.00 ,3.20 ,3.40 ,3.60 ,3.80  , 4.00 ,
     +     5.00 /
      DATA SOLWAV/.0000,
     +     .0001,.0005,.0010,.0014,.0019,.0027 ,.0041 ,
     +     .0056,
     +     .0081,.0121,.0165,.0222,.0293,.0372 ,.0452 ,
     +     .0532,.0615,.0700,.0782,.0873,.0992 ,.1122 ,
     +     .1247,.1373,.1514,.1665,.1817,.1968 ,.2115 ,
     +     .2260,.2401,.2538,.2674,.2808,.2938 ,.3065 ,
     +     .3191,
     +     .3318,.3444,.3568,.3810,.4042,.4266 ,.4481 ,
     +     .4688,.4886,.5169,.5602,.6336,.6946 ,.7839 ,
     +     .8434,.8861,.9159,.9349,.9483,.9589 ,.9667 ,
     +     .9731,.9783,.9822,.9850,.9872,.9891 ,.9906 ,
     +     .9951/
      
      DATA ITER / 0 /
      
      IF( ITER .EQ. 0 ) THEN
         DO 10 N=1,MAXWAV
            WAVNUM(N) = 10000. / LAMDAT(N)
 10      CONTINUE
      ENDIF
      ITER = ITER + 1
      
      IF( WAVLOW .GT. 50000. ) THEN
         SLFWTD = 0.0
         RETURN
      ENDIF
      
      IF( WAVLOW .GE. 2000. ) THEN
         M = 1
         DO 100 N=2,MAXWAV
            
C     WAVENUMBERS ORDERED IN ARRAY 'WAVNUM' FROM LARGE TO SMALL:
            IF( WAVNUM(N-1) .GE. WAVLOW .AND. WAVLOW .GE. WAVNUM(N)) THEN
               M = N
               GOTO 101
            ENDIF
 100     CONTINUE
         
 101     CONTINUE
         DMAX = SOLWAV(M-1) + ((SOLWAV(M) - SOLWAV(M-1)) / (WAVNUM(M)
     1        - WAVNUM(M-1))) * ( WAVLOW - WAVNUM(M-1))
      ELSE
         DMAX = 1.0
      ENDIF
      
      IF( WAVHGH .LT. 2000. ) THEN
         SLFWTD = 0.0
         RETURN
      ENDIF
      
      IF( WAVHGH .LE. 50000. ) THEN
         M = 1
         DO 200 N=2,MAXWAV
            
C     WAVENUMBERS ORDERED IN ARRAY 'WAVNUM' FROM LARGE TO SMALL:
            IF( WAVNUM(N-1) .GE. WAVHGH .AND. WAVHGH .GE. WAVNUM(N)) THEN
               M = N
               GOTO 201
            ENDIF
 200     CONTINUE
         
 201     CONTINUE
         DMIN = SOLWAV(M-1) + ((SOLWAV(M) - SOLWAV(M-1)) / (WAVNUM(M)
     1        - WAVNUM(M-1))) * ( WAVHGH - WAVNUM(M-1))
      ELSE
         DMIN = 0.0
      ENDIF
      
C     SOLAR FLUX FRACTION
      SLFWTD =  DMAX - DMIN
      RETURN
      END
      
      FUNCTION SLFFLN(LMIN,LMAX)
      PARAMETER (MAXWAV=122)
      REAL LAMBDA(MAXWAV),DLAMBD(MAXWAV)
      REAL LMIN,LMAX
      
C     SOLAR FLUX FRACTIONS ACCORDING TO:
C     LABS AND NECKEL (1968)
      
C     INTEGRATING VALUES FROM 0.200 MICRONS   TO  4.950 MICRONS
C     GIVES A SOLAR CONSTANT OF: 1360.39 WATTS/M2
      
C     VALUES BELOW ARE ACTUALLY THE FRACTION OF SOLAR FLUX SHORTWAVE
C     OF THE PARTICULAR WAVELENGTH.
      
      DATA LAMBDA/ .200,.205,.215,.225,.235,.245,.255,.265,.275,
     +     .285,.295,
     +     .305,.315,.325,.335,.345,.355,.365,.375,.385,.395,
     +     .405,.415,.425,.435,.445,.455,.465,.475,.485,.495,
     +     .505,.515,.525,.535,.545,.555,.565,.575,.585,.595,
     +     .605,.615,.625,.635,.645,.655,.665,.675,.685,.695,
     +     .705,.715,.725,.735,.745,.755,.765,.775,.785,.795,
     +     .805,.815,.825,.835,.845,.855,.865,.875,.885,.895,
     +     .905,.915,.925,.935,.945,.955,.965,.975,.985,.995,
     +     1.05,1.15,1.25,1.35,1.45,1.55,1.65,1.75,1.85,1.95,
     +     2.05,2.15,2.25,2.35,2.45,2.55,2.65,2.75,2.85,2.95,
     +     3.05,3.15,3.25,3.35,3.45,3.55,3.65,3.75,3.85,3.95,
     +     4.05,4.15,4.25,4.35,4.45,4.55,4.65,4.75,4.85,4.95,
     +     5.00/
      
      DATA DLAMBD  / .0000,
     +     .0001,.0003,.0006,.0010,.0015,.0020,.0029,.0042,.0059,.0088,
     +     .0127,.0172,.0227,.0291,.0358,.0427,.0502,.0580,.0654,.0732,
     +     .0835,.0959,.1085,.1209,.1343,.1489,.1638,.1786,.1930,.2073,
     +     .2217,.2356,.2493,.2633,.2774,.2911,.3047,.3183,.3318,.3451,
     +     .3581,.3709,.3834,.3956,.4076,.4191,.4304,.4417,.4528,.4636,
     +     .4741,.4844,.4945,.5043,.5139,.5232,.5324,.5414,.5502,.5588,
     +     .5672,.5755,.5835,.5913,.5989,.6062,.6134,.6204,.6273,.6340,
     +     .6407,.6472,.6536,.6598,.6660,.6719,.6778,.6835,.6892,.6947,
     +     .7230,.7671,.8034,.8336,.8590,.8808,.8992,.9144,.9269,.9372,
     +     .9457,.9529,.9590,.9642,.9686,.9724,.9757,.9786,.9811,.9834,
     +     .9853,.9870,.9885,.9899,.9911,.9922,.9931,.9940,.9948,.9955,
     +     .9962,.9968,.9973,.9978,.9982,.9987,.9990,.9994,.9997,.9999,
     +     1.0/
      
      IF( LMIN .GE. 0.200 .AND. LMIN .LT. 5.000 ) THEN
         
C     FIND LOWER LIMIT
         N = 1
 100     IF( LMIN .GE. LAMBDA(N) ) THEN
            N = N + 1
            GO TO 100
         ENDIF
         DMIN = DLAMBD(N-1) + ((DLAMBD(N) - DLAMBD(N-1)) / (LAMBDA(N)
     1        - LAMBDA(N-1))) * ( LMIN - LAMBDA(N-1))
      ELSE IF( LMIN .LT. 0.200 ) THEN
         DMIN = 0.0
c++csz
c     ELSE IF( LMIN .GT. 5.000 ) THEN
      ELSE IF( LMIN .GE. 5.000 ) THEN
c--csz
         SLFFLN =  0.0
         RETURN
      ENDIF
      
C     FIND UPPER LIMIT
      
c++csz
c     IF( LMAX .LE. 5.000 ) THEN
      IF( LMAX .GT. 0.200 .AND. LMAX .LT. 5. ) THEN
c--csz
         N = 1
 200     IF( LMAX .GE. LAMBDA(N) ) THEN
            N = N + 1
            GO TO 200
         ENDIF
         DMAX = DLAMBD(N-1) + ((DLAMBD(N) - DLAMBD(N-1)) / (LAMBDA(N)
     1        - LAMBDA(N-1))) * ( LMAX - LAMBDA(N-1))
      ELSE IF( LMAX .GE. 5.000 ) THEN
         DMAX = 1.0000
c++csz
c     ELSE IF( LMAX .LT. 0.200 ) THEN
      ELSE IF( LMAX .LE. 0.200 ) THEN
c--csz
         SLFFLN =  0.0
         RETURN
      ENDIF
      
C     SOLAR FLUX FRACTION
      SLFFLN =  DMAX - DMIN
      RETURN
      END
      
      FUNCTION SLFWLN(WAVHGH,WAVLOW)
      
C     COMPUTES SOLAR FLUX FRACTION BETWEEN HIGHEST (WAVGHG) AND LOWEST
C     (WAVLOW) WAVENUMBER:
      PARAMETER (MAXWAV=122)
      real lamdat,wavnum,solwav
      common /slndat / LAMDAT(MAXWAV),WAVNUM(MAXWAV),SOLWAV(MAXWAV)
      
C     SOLAR FLUX FRACTIONS ACCORDING TO:
C     LABS AND NECKEL (1968)
      
C     INTEGRATING VALUES FROM 0.200 MICRONS   TO  4.950 MICRONS
C     GIVES A SOLAR CONSTANT OF: 1360.39 WATTS/M2
      
C     VALUES BELOW ARE ACTUALLY THE FRACTION OF SOLAR FLUX SHORTWAVE
C     OF THE PARTICULAR WAVELENGTH.
      
      DATA LAMDAT/ .200,.205,.215,.225,.235,.245,.255,.265,.275,
     +     .285,.295,
     +     .305,.315,.325,.335,.345,.355,.365,.375,.385,.395,
     +     .405,.415,.425,.435,.445,.455,.465,.475,.485,.495,
     +     .505,.515,.525,.535,.545,.555,.565,.575,.585,.595,
     +     .605,.615,.625,.635,.645,.655,.665,.675,.685,.695,
     +     .705,.715,.725,.735,.745,.755,.765,.775,.785,.795,
     +     .805,.815,.825,.835,.845,.855,.865,.875,.885,.895,
     +     .905,.915,.925,.935,.945,.955,.965,.975,.985,.995,
     +     1.05,1.15,1.25,1.35,1.45,1.55,1.65,1.75,1.85,1.95,
     +     2.05,2.15,2.25,2.35,2.45,2.55,2.65,2.75,2.85,2.95,
     +     3.05,3.15,3.25,3.35,3.45,3.55,3.65,3.75,3.85,3.95,
     +     4.05,4.15,4.25,4.35,4.45,4.55,4.65,4.75,4.85,4.95,
     +     5.00/
      
      DATA SOLWAV  / .0000,
     +     .0001,.0003,.0006,.0010,.0015,.0020,.0029,.0042,.0059,.0088,
     +     .0127,.0172,.0227,.0291,.0358,.0427,.0502,.0580,.0654,.0732,
     +     .0835,.0959,.1085,.1209,.1343,.1489,.1638,.1786,.1930,.2073,
     +     .2217,.2356,.2493,.2633,.2774,.2911,.3047,.3183,.3318,.3451,
     +     .3581,.3709,.3834,.3956,.4076,.4191,.4304,.4417,.4528,.4636,
     +     .4741,.4844,.4945,.5043,.5139,.5232,.5324,.5414,.5502,.5588,
     +     .5672,.5755,.5835,.5913,.5989,.6062,.6134,.6204,.6273,.6340,
     +     .6407,.6472,.6536,.6598,.6660,.6719,.6778,.6835,.6892,.6947,
     +     .7230,.7671,.8034,.8336,.8590,.8808,.8992,.9144,.9269,.9372,
     +     .9457,.9529,.9590,.9642,.9686,.9724,.9757,.9786,.9811,.9834,
     +     .9853,.9870,.9885,.9899,.9911,.9922,.9931,.9940,.9948,.9955,
     +     .9962,.9968,.9973,.9978,.9982,.9987,.9990,.9994,.9997,.9999,
     +     1.0000/
      
      DATA ITER / 0 /
      
      IF( ITER .EQ. 0 ) THEN
         DO 10 N=1,MAXWAV
            WAVNUM(N) = 10000. / LAMDAT(N)
 10      CONTINUE
      ENDIF
      ITER = ITER + 1
      
      IF( WAVLOW .GT. 50000. ) THEN
         SLFWLN = 0.0
         RETURN
      ENDIF
      
      IF( WAVLOW .GE. 2000. ) THEN
         M = 1
         DO 100 N=2,MAXWAV
            
C     WAVENUMBERS ORDERED IN ARRAY 'WAVNUM' FROM LARGE TO SMALL:
            IF( WAVNUM(N-1) .GE. WAVLOW .AND. WAVLOW .GE. WAVNUM(N)) THEN
               M = N
            ENDIF
 100     CONTINUE
         
 101     CONTINUE
         DMAX = SOLWAV(M-1) + ((SOLWAV(M) - SOLWAV(M-1)) / (WAVNUM(M)
     1        - WAVNUM(M-1))) * ( WAVLOW - WAVNUM(M-1))
      ELSE
         DMAX = 1.0
      ENDIF
      
      IF( WAVHGH .LT. 2000. ) THEN
         SLFWLN = 0.0
         RETURN
      ENDIF
      
      IF( WAVHGH .LE. 50000. ) THEN
         M = 1
         DO 200 N=2,MAXWAV
            
C     WAVENUMBERS ORDERED IN ARRAY 'WAVNUM' FROM LARGE TO SMALL:
            IF( WAVNUM(N-1) .GE. WAVHGH .AND. WAVHGH .GE. WAVNUM(N)) THEN
               M = N
            ENDIF
 200     CONTINUE
         
 201     CONTINUE
         DMIN = SOLWAV(M-1) + ((SOLWAV(M) - SOLWAV(M-1)) / (WAVNUM(M)
     1        - WAVNUM(M-1))) * ( WAVHGH - WAVNUM(M-1))
      ELSE
         DMIN = 0.0
      ENDIF
      
C     SOLAR FLUX FRACTION
      SLFWLN =  DMAX - DMIN
      RETURN
      END
      
      
