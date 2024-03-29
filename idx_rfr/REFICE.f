
      COMPLEX FUNCTION REFICE( WAVLEN, TEMP)

c        Calculates complex refractive index of Ice 1H for wavelengths
c        between 45 nm and 8.6 m.  For wavelengths above 167 microns,
c        temperature dependence is included for temperatures between
c        213 and 272K.

c      I N P U T :  WAVLEN = wavelength (microns)
c                   TEMP   = temperature (K) ( for WAVLEN.GT.167 only )

c      O U T P U T :  REFICE = complex refractive index
c                              ( with positive imaginary part )

c      METHOD :  Tabular interpolation, assuming

c                (1) real index is linear in log(wavelength)
c                    and linear in temperature

c                (2) log(imag. index) is linear in log(wavelength)
c                    and linear in temperature

c     AUTHOR :  Prof. Stephen Warren, Univ. of Washington (Sept., 1983)

c     REFERENCE :  Warren, S., 1984:  Optical Constants of Ice from the
c                     Ultraviolet to the Microwave, Appl. Opt. 23,
c                     1206-1225


c     .. Parameters ..

      INTEGER   NWL, NWLT
      PARAMETER ( NWL = 468, NWLT = 62)
c     ..
c     .. Scalar Arguments ..

      REAL      TEMP, WAVLEN
c     ..
c     .. Local Scalars ..

      INTEGER   I, L
      REAL      FRAC, MIM, MRE, YHI, YLO
c     ..
c     .. External Subroutines ..

      EXTERNAL  ERRMSG
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC LOG, CMPLX, EXP
c     ..
c     .. Common blocks ..

c                                        ** Refractive index table

      COMMON / ICEREF / WL, WLT, TABRE, TABRET, TABIM, TABIMT, TEMREF

      REAL      TABIM( NWL), TABIMT( NWLT, 4), TABRE( NWL),
     &          TABRET( NWLT, 4), TEMREF( 4), WL( NWL), WLT( NWLT)
c     ..


      IF( WAVLEN.LT.WL(1) .OR. WAVLEN.GT.WLT(NWLT)) THEN

         CALL ERRMSG( 'REFICE--WAVLEN outside table boundaries',.TRUE.)

         STOP

      END IF


      IF( WAVLEN.LE.167.) THEN
c                                  ** Wavelength between 0.045 and 167
c                                  ** microns. No temperature dependence
         DO 10 I = 2, NWL

            IF( WAVLEN.LE.WL(I)) GO TO 20

   10    CONTINUE

   20    CONTINUE
         FRAC   = LOG( WAVLEN / WL(I-1)) / LOG( WL(I) / WL(I-1))
         MRE    = TABRE(I-1) + FRAC*( TABRE(I) - TABRE(I-1))
         MIM    = TABIM(I-1) * ( TABIM(I) / TABIM(I-1) )**FRAC

      ELSE
c               ** Wavelength greater than 167 microns
c               ** (temperature-dependent case)

         IF( TEMP.LT.TEMREF(4) .OR. TEMP.GT.TEMREF(1)) THEN

            CALL ERRMSG( 'REFICE--TEMP OUTSIDE TABLE BOUNDARIES',.TRUE.)

            STOP

         END IF
c                         ** Find position in temperature array

         DO 30 L = 2, 4

            IF( TEMP.GE.TEMREF(L)) GO TO 40

   30    CONTINUE
c                         ** Find position in wavelength array

   40    CONTINUE
         DO 50 I = 2, NWLT

            IF( WAVLEN.LE.WLT(I)) GO TO 60

   50    CONTINUE


   60    CONTINUE
         FRAC   = LOG( WAVLEN / WLT(I-1)) /
     &            LOG( WLT(I) / WLT(I-1))

         YLO    = TABRET( I-1, L) + 
     &            FRAC*( TABRET(I, L) - TABRET(I-1, L))
         YHI    = TABRET( I-1, L-1) +
     &            FRAC*( TABRET(I, L-1) - TABRET(I-1, L-1))
         MRE    = YLO + ( YHI - YLO)*( TEMP - TEMREF(L)) /
     &                               ( TEMREF(L-1) - TEMREF(L))

         YLO    = LOG( TABIMT(I-1, L)) +
     &            FRAC*LOG( TABIMT(I, L) / TABIMT(I-1, L))
         YHI    = LOG( TABIMT(I-1, L-1)) +
     &            FRAC*LOG( TABIMT(I, L-1) / TABIMT(I-1, L-1))
         MIM    = EXP( YLO + (YHI - YLO)*(TEMP - TEMREF(L)) /
     &                                   (TEMREF(L-1) - TEMREF(L)))

      END IF


      REFICE = CMPLX( MRE, MIM)

      END

      BLOCK DATA ICECON

c        Ice-refractive-index vs. wavelength table for function REFICE

c     .. Parameters ..

      INTEGER   NWL, NWLT
      PARAMETER ( NWL = 468, NWLT = 62)
c     ..
c     .. Local Scalars ..

      INTEGER   I
c     ..
c     .. Common blocks ..

      COMMON / ICEREF / WL, WLT, TABRE, TABRET, TABIM, TABIMT, TEMREF

      REAL      TABIM( NWL), TABIMT( NWLT, 4), TABRE( NWL),
     &          TABRET( NWLT, 4), TEMREF( 4), WL( NWL), WLT( NWLT)
c     ..
c     .. Data statements ..

c        WL, WLT          wavelengths (microns) for temperature-
c                         independent and temperature-dependent
c                         regimes, respectively

c       TABRE, TABRET     real refractive indices for temperature-
c                         independent and temperature-dependent
c                         regimes, respectively

c       TABIM, TABIMT     imaginary refractive indices for temperature-
c                         independent and temperature-dependent
c                         regimes, respectively

c       TEMREF            reference temperatures (-1,-5,-20,-60 deg C)
c                         for TABRET,TABIMT

      DATA  TEMREF / 272.16, 268.16, 253.16, 213.16 /

      DATA ( WL(I), I = 1,114 ) /
     &0.4430E-01,0.4510E-01,0.4590E-01,0.4680E-01,0.4770E-01,0.4860E-01,
     &0.4960E-01,0.5060E-01,0.5170E-01,0.5280E-01,0.5390E-01,0.5510E-01,
     &0.5640E-01,0.5770E-01,0.5900E-01,0.6050E-01,0.6200E-01,0.6360E-01,
     &0.6530E-01,0.6700E-01,0.6890E-01,0.7080E-01,0.7290E-01,0.7380E-01,
     &0.7510E-01,0.7750E-01,0.8000E-01,0.8270E-01,0.8550E-01,0.8860E-01,
     &0.9180E-01,0.9300E-01,0.9540E-01,0.9920E-01,0.1033E+00,0.1078E+00,
     &0.1100E+00,0.1127E+00,0.1140E+00,0.1181E+00,0.1210E+00,0.1240E+00,
     &0.1272E+00,0.1295E+00,0.1305E+00,0.1319E+00,0.1333E+00,0.1348E+00,
     &0.1362E+00,0.1370E+00,0.1378E+00,0.1387E+00,0.1393E+00,0.1409E+00,
     &0.1425E+00,0.1435E+00,0.1442E+00,0.1450E+00,0.1459E+00,0.1468E+00,
     &0.1476E+00,0.1480E+00,0.1485E+00,0.1494E+00,0.1512E+00,0.1531E+00,
     &0.1540E+00,0.1550E+00,0.1569E+00,0.1580E+00,0.1589E+00,0.1610E+00,
     &0.1625E+00,0.1648E+00,0.1669E+00,0.1692E+00,0.1713E+00,0.1737E+00,
     &0.1757E+00,0.1779E+00,0.1802E+00,0.1809E+00,0.1821E+00,0.1833E+00,
     &0.1843E+00,0.1850E+00,0.1860E+00,0.1870E+00,0.1880E+00,0.1890E+00,
     &0.1900E+00,0.1910E+00,0.1930E+00,0.1950E+00,0.2100E+00,0.2500E+00,
     &0.3000E+00,0.3500E+00,0.4000E+00,0.4100E+00,0.4200E+00,0.4300E+00,
     &0.4400E+00,0.4500E+00,0.4600E+00,0.4700E+00,0.4800E+00,0.4900E+00,
     &0.5000E+00,0.5100E+00,0.5200E+00,0.5300E+00,0.5400E+00,0.5500E+00/
      DATA ( WL(I), I = 115,228 ) /
     &0.5600E+00,0.5700E+00,0.5800E+00,0.5900E+00,0.6000E+00,0.6100E+00,
     &0.6200E+00,0.6300E+00,0.6400E+00,0.6500E+00,0.6600E+00,0.6700E+00,
     &0.6800E+00,0.6900E+00,0.7000E+00,0.7100E+00,0.7200E+00,0.7300E+00,
     &0.7400E+00,0.7500E+00,0.7600E+00,0.7700E+00,0.7800E+00,0.7900E+00,
     &0.8000E+00,0.8100E+00,0.8200E+00,0.8300E+00,0.8400E+00,0.8500E+00,
     &0.8600E+00,0.8700E+00,0.8800E+00,0.8900E+00,0.9000E+00,0.9100E+00,
     &0.9200E+00,0.9300E+00,0.9400E+00,0.9500E+00,0.9600E+00,0.9700E+00,
     &0.9800E+00,0.9900E+00,0.1000E+01,0.1010E+01,0.1020E+01,0.1030E+01,
     &0.1040E+01,0.1050E+01,0.1060E+01,0.1070E+01,0.1080E+01,0.1090E+01,
     &0.1100E+01,0.1110E+01,0.1120E+01,0.1130E+01,0.1140E+01,0.1150E+01,
     &0.1160E+01,0.1170E+01,0.1180E+01,0.1190E+01,0.1200E+01,0.1210E+01,
     &0.1220E+01,0.1230E+01,0.1240E+01,0.1250E+01,0.1260E+01,0.1270E+01,
     &0.1280E+01,0.1290E+01,0.1300E+01,0.1310E+01,0.1320E+01,0.1330E+01,
     &0.1340E+01,0.1350E+01,0.1360E+01,0.1370E+01,0.1380E+01,0.1390E+01,
     &0.1400E+01,0.1410E+01,0.1420E+01,0.1430E+01,0.1440E+01,0.1449E+01,
     &0.1460E+01,0.1471E+01,0.1481E+01,0.1493E+01,0.1504E+01,0.1515E+01,
     &0.1527E+01,0.1538E+01,0.1563E+01,0.1587E+01,0.1613E+01,0.1650E+01,
     &0.1680E+01,0.1700E+01,0.1730E+01,0.1760E+01,0.1800E+01,0.1830E+01,
     &0.1840E+01,0.1850E+01,0.1855E+01,0.1860E+01,0.1870E+01,0.1890E+01/
      DATA ( WL(I), I = 229,342 ) /
     &0.1905E+01,0.1923E+01,0.1942E+01,0.1961E+01,0.1980E+01,0.2000E+01,
     &0.2020E+01,0.2041E+01,0.2062E+01,0.2083E+01,0.2105E+01,0.2130E+01,
     &0.2150E+01,0.2170E+01,0.2190E+01,0.2220E+01,0.2240E+01,0.2245E+01,
     &0.2250E+01,0.2260E+01,0.2270E+01,0.2290E+01,0.2310E+01,0.2330E+01,
     &0.2350E+01,0.2370E+01,0.2390E+01,0.2410E+01,0.2430E+01,0.2460E+01,
     &0.2500E+01,0.2520E+01,0.2550E+01,0.2565E+01,0.2580E+01,0.2590E+01,
     &0.2600E+01,0.2620E+01,0.2675E+01,0.2725E+01,0.2778E+01,0.2817E+01,
     &0.2833E+01,0.2849E+01,0.2865E+01,0.2882E+01,0.2899E+01,0.2915E+01,
     &0.2933E+01,0.2950E+01,0.2967E+01,0.2985E+01,0.3003E+01,0.3021E+01,
     &0.3040E+01,0.3058E+01,0.3077E+01,0.3096E+01,0.3115E+01,0.3135E+01,
     &0.3155E+01,0.3175E+01,0.3195E+01,0.3215E+01,0.3236E+01,0.3257E+01,
     &0.3279E+01,0.3300E+01,0.3322E+01,0.3345E+01,0.3367E+01,0.3390E+01,
     &0.3413E+01,0.3436E+01,0.3460E+01,0.3484E+01,0.3509E+01,0.3534E+01,
     &0.3559E+01,0.3624E+01,0.3732E+01,0.3775E+01,0.3847E+01,0.3969E+01,
     &0.4099E+01,0.4239E+01,0.4348E+01,0.4387E+01,0.4444E+01,0.4505E+01,
     &0.4547E+01,0.4560E+01,0.4580E+01,0.4719E+01,0.4904E+01,0.5000E+01,
     &0.5100E+01,0.5200E+01,0.5263E+01,0.5400E+01,0.5556E+01,0.5714E+01,
     &0.5747E+01,0.5780E+01,0.5814E+01,0.5848E+01,0.5882E+01,0.6061E+01,
     &0.6135E+01,0.6250E+01,0.6289E+01,0.6329E+01,0.6369E+01,0.6410E+01/
      DATA ( WL(I), I = 343,456 ) /
     &0.6452E+01,0.6494E+01,0.6579E+01,0.6667E+01,0.6757E+01,0.6897E+01,
     &0.7042E+01,0.7143E+01,0.7246E+01,0.7353E+01,0.7463E+01,0.7576E+01,
     &0.7692E+01,0.7812E+01,0.7937E+01,0.8065E+01,0.8197E+01,0.8333E+01,
     &0.8475E+01,0.8696E+01,0.8929E+01,0.9091E+01,0.9259E+01,0.9524E+01,
     &0.9804E+01,0.1000E+02,0.1020E+02,0.1031E+02,0.1042E+02,0.1053E+02,
     &0.1064E+02,0.1075E+02,0.1087E+02,0.1100E+02,0.1111E+02,0.1136E+02,
     &0.1163E+02,0.1190E+02,0.1220E+02,0.1250E+02,0.1282E+02,0.1299E+02,
     &0.1316E+02,0.1333E+02,0.1351E+02,0.1370E+02,0.1389E+02,0.1408E+02,
     &0.1429E+02,0.1471E+02,0.1515E+02,0.1538E+02,0.1563E+02,0.1613E+02,
     &0.1639E+02,0.1667E+02,0.1695E+02,0.1724E+02,0.1818E+02,0.1887E+02,
     &0.1923E+02,0.1961E+02,0.2000E+02,0.2041E+02,0.2083E+02,0.2222E+02,
     &0.2260E+02,0.2305E+02,0.2360E+02,0.2460E+02,0.2500E+02,0.2600E+02,
     &0.2857E+02,0.3100E+02,0.3333E+02,0.3448E+02,0.3564E+02,0.3700E+02,
     &0.3824E+02,0.3960E+02,0.4114E+02,0.4276E+02,0.4358E+02,0.4458E+02,
     &0.4550E+02,0.4615E+02,0.4671E+02,0.4736E+02,0.4800E+02,0.4878E+02,
     &0.5003E+02,0.5128E+02,0.5275E+02,0.5350E+02,0.5424E+02,0.5500E+02,
     &0.5574E+02,0.5640E+02,0.5700E+02,0.5746E+02,0.5840E+02,0.5929E+02,
     &0.6000E+02,0.6100E+02,0.6125E+02,0.6250E+02,0.6378E+02,0.6467E+02,
     &0.6558E+02,0.6655E+02,0.6760E+02,0.6900E+02,0.7053E+02,0.7300E+02/
      DATA ( WL(I), I = 457,468 ) /
     &0.7500E+02,0.7629E+02,0.8000E+02,0.8297E+02,0.8500E+02,0.8680E+02,
     &0.9080E+02,0.9517E+02,0.1000E+03,0.1200E+03,0.1500E+03,0.1670E+03/

      DATA  WLT /                      0.1670E+03,0.1778E+03,0.1884E+03,
     &0.1995E+03,0.2113E+03,0.2239E+03,0.2371E+03,0.2512E+03,0.2661E+03,
     &0.2818E+03,0.2985E+03,0.3162E+03,0.3548E+03,0.3981E+03,0.4467E+03,
     &0.5012E+03,0.5623E+03,0.6310E+03,0.7943E+03,0.1000E+04,0.1259E+04,
     &0.2500E+04,0.5000E+04,0.1000E+05,0.2000E+05,0.3200E+05,0.3500E+05,
     &0.4000E+05,0.4500E+05,0.5000E+05,0.6000E+05,0.7000E+05,0.9000E+05,
     &0.1110E+06,0.1200E+06,0.1300E+06,0.1400E+06,0.1500E+06,0.1600E+06,
     &0.1700E+06,0.1800E+06,0.2000E+06,0.2500E+06,0.2900E+06,0.3200E+06,
     &0.3500E+06,0.3800E+06,0.4000E+06,0.4500E+06,0.5000E+06,0.6000E+06,
     &0.6400E+06,0.6800E+06,0.7200E+06,0.7600E+06,0.8000E+06,0.8400E+06,
     &0.9000E+06,0.1000E+07,0.2000E+07,0.5000E+07,0.8600E+07/

      DATA ( TABRE(I), I = 1,114 ) /
     &   0.83441,   0.83676,   0.83729,   0.83771,   0.83827,   0.84038,
     &   0.84719,   0.85522,   0.86047,   0.86248,   0.86157,   0.86093,
     &   0.86419,   0.86916,   0.87764,   0.89296,   0.91041,   0.93089,
     &   0.95373,   0.98188,   1.02334,   1.06735,   1.11197,   1.13134,
     &   1.15747,   1.20045,   1.23840,   1.27325,   1.32157,   1.38958,
     &   1.41644,   1.40906,   1.40063,   1.40169,   1.40934,   1.40221,
     &   1.39240,   1.38424,   1.38075,   1.38186,   1.39634,   1.40918,
     &   1.40256,   1.38013,   1.36303,   1.34144,   1.32377,   1.30605,
     &   1.29054,   1.28890,   1.28931,   1.30190,   1.32025,   1.36302,
     &   1.41872,   1.45834,   1.49028,   1.52128,   1.55376,   1.57782,
     &   1.59636,   1.60652,   1.61172,   1.61919,   1.62522,   1.63404,
     &   1.63689,   1.63833,   1.63720,   1.63233,   1.62222,   1.58269,
     &   1.55635,   1.52453,   1.50320,   1.48498,   1.47226,   1.45991,
     &   1.45115,   1.44272,   1.43498,   1.43280,   1.42924,   1.42602,
     &   1.42323,   1.42143,   1.41897,   1.41660,   1.41434,   1.41216,
     &   1.41006,   1.40805,   1.40423,   1.40067,   1.38004,   1.35085,
     &   1.33394,   1.32492,   1.31940,   1.31854,   1.31775,   1.31702,
     &   1.31633,   1.31569,   1.31509,   1.31452,   1.31399,   1.31349,
     &   1.31302,   1.31257,   1.31215,   1.31175,   1.31136,   1.31099/
      DATA ( TABRE(I), I = 115,228 ) /
     &   1.31064,   1.31031,   1.30999,   1.30968,   1.30938,   1.30909,
     &   1.30882,   1.30855,   1.30829,   1.30804,   1.30780,   1.30756,
     &   1.30733,   1.30710,   1.30688,   1.30667,   1.30646,   1.30625,
     &   1.30605,   1.30585,   1.30566,   1.30547,   1.30528,   1.30509,
     &   1.30491,   1.30473,   1.30455,   1.30437,   1.30419,   1.30402,
     &   1.30385,   1.30367,   1.30350,   1.30333,   1.30316,   1.30299,
     &   1.30283,   1.30266,   1.30249,   1.30232,   1.30216,   1.30199,
     &   1.30182,   1.30166,   1.30149,   1.30132,   1.30116,   1.30099,
     &   1.30082,   1.30065,   1.30048,   1.30031,   1.30014,   1.29997,
     &   1.29979,   1.29962,   1.29945,   1.29927,   1.29909,   1.29891,
     &   1.29873,   1.29855,   1.29837,   1.29818,   1.29800,   1.29781,
     &   1.29762,   1.29743,   1.29724,   1.29705,   1.29686,   1.29666,
     &   1.29646,   1.29626,   1.29605,   1.29584,   1.29563,   1.29542,
     &   1.29521,   1.29499,   1.29476,   1.29453,   1.29430,   1.29406,
     &   1.29381,   1.29355,   1.29327,   1.29299,   1.29272,   1.29252,
     &   1.29228,   1.29205,   1.29186,   1.29167,   1.29150,   1.29130,
     &   1.29106,   1.29083,   1.29025,   1.28962,   1.28891,   1.28784,
     &   1.28689,   1.28623,   1.28521,   1.28413,   1.28261,   1.28137,
     &   1.28093,   1.28047,   1.28022,   1.27998,   1.27948,   1.27849/
      DATA ( TABRE(I), I = 229,342 ) /
     &   1.27774,   1.27691,   1.27610,   1.27535,   1.27471,   1.27404,
     &   1.27329,   1.27240,   1.27139,   1.27029,   1.26901,   1.26736,
     &   1.26591,   1.26441,   1.26284,   1.26036,   1.25860,   1.25815,
     &   1.25768,   1.25675,   1.25579,   1.25383,   1.25179,   1.24967,
     &   1.24745,   1.24512,   1.24266,   1.24004,   1.23725,   1.23270,
     &   1.22583,   1.22198,   1.21548,   1.21184,   1.20790,   1.20507,
     &   1.20209,   1.19566,   1.17411,   1.14734,   1.10766,   1.06739,
     &   1.04762,   1.02650,   1.00357,   0.98197,   0.96503,   0.95962,
     &   0.97269,   0.99172,   1.00668,   1.02186,   1.04270,   1.07597,
     &   1.12954,   1.21267,   1.32509,   1.42599,   1.49656,   1.55095,
     &   1.59988,   1.63631,   1.65024,   1.64278,   1.62691,   1.61284,
     &   1.59245,   1.57329,   1.55770,   1.54129,   1.52654,   1.51139,
     &   1.49725,   1.48453,   1.47209,   1.46125,   1.45132,   1.44215,
     &   1.43366,   1.41553,   1.39417,   1.38732,   1.37735,   1.36448,
     &   1.35414,   1.34456,   1.33882,   1.33807,   1.33847,   1.34053,
     &   1.34287,   1.34418,   1.34634,   1.34422,   1.33453,   1.32897,
     &   1.32333,   1.31800,   1.31432,   1.30623,   1.29722,   1.28898,
     &   1.28730,   1.28603,   1.28509,   1.28535,   1.28813,   1.30156,
     &   1.30901,   1.31720,   1.31893,   1.32039,   1.32201,   1.32239/
      DATA ( TABRE(I), I = 343,456 ) /
     &   1.32149,   1.32036,   1.31814,   1.31705,   1.31807,   1.31953,
     &   1.31933,   1.31896,   1.31909,   1.31796,   1.31631,   1.31542,
     &   1.31540,   1.31552,   1.31455,   1.31193,   1.30677,   1.29934,
     &   1.29253,   1.28389,   1.27401,   1.26724,   1.25990,   1.24510,
     &   1.22241,   1.19913,   1.17150,   1.15528,   1.13700,   1.11808,
     &   1.10134,   1.09083,   1.08734,   1.09254,   1.10654,   1.14779,
     &   1.20202,   1.25825,   1.32305,   1.38574,   1.44478,   1.47170,
     &   1.49619,   1.51652,   1.53328,   1.54900,   1.56276,   1.57317,
     &   1.58028,   1.57918,   1.56672,   1.55869,   1.55081,   1.53807,
     &   1.53296,   1.53220,   1.53340,   1.53289,   1.51705,   1.50097,
     &   1.49681,   1.49928,   1.50153,   1.49856,   1.49053,   1.46070,
     &   1.45182,   1.44223,   1.43158,   1.41385,   1.40676,   1.38955,
     &   1.34894,   1.31039,   1.26420,   1.23656,   1.21663,   1.20233,
     &   1.19640,   1.19969,   1.20860,   1.22173,   1.24166,   1.28175,
     &   1.32784,   1.38657,   1.46486,   1.55323,   1.60379,   1.61877,
     &   1.62963,   1.65712,   1.69810,   1.72065,   1.74865,   1.76736,
     &   1.76476,   1.75011,   1.72327,   1.68490,   1.62398,   1.59596,
     &   1.58514,   1.59917,   1.61405,   1.66625,   1.70663,   1.73713,
     &   1.76860,   1.80343,   1.83296,   1.85682,   1.87411,   1.89110/
      DATA ( TABRE(I), I = 457,468 ) /
     &   1.89918,   1.90432,   1.90329,   1.88744,   1.87499,   1.86702,
     &   1.85361,   1.84250,   1.83225,   1.81914,   1.82268,   1.82961/

      DATA ( TABRET(I,1), I = 1,NWLT ) /  1.82961,   1.83258,   1.83149,
     &   1.82748,   1.82224,   1.81718,   1.81204,   1.80704,   1.80250,
     &   1.79834,   1.79482,   1.79214,   1.78843,   1.78601,   1.78434,
     &   1.78322,   1.78248,   1.78201,   1.78170,   1.78160,   1.78190,
     &   1.78300,   1.78430,   1.78520,   1.78620,   1.78660,   1.78680,
     &   1.78690,   1.78700,   1.78700,   1.78710,   1.78710,   1.78720,
     &   1.78720,   1.78720,   1.78720,   1.78720,   1.78720,   1.78720,
     &   1.78720,   1.78720,   1.78720,   1.78720,   1.78720,   1.78720,
     &   1.78720,   1.78720,   1.78720,   1.78720,   1.78720,   1.78720,
     &   1.78720,   1.78720,   1.78720,   1.78720,   1.78720,   1.78720,
     &   1.78720,   1.78720,   1.78720,   1.78720,   1.78800/
      DATA ( TABRET(I,2), I = 1,NWLT ) /
     &                         1.82961,   1.83258,   1.83149,   1.82748,
     &   1.82224,   1.81718,   1.81204,   1.80704,   1.80250,   1.79834,
     &   1.79482,   1.79214,   1.78843,   1.78601,   1.78434,   1.78322,
     &   1.78248,   1.78201,   1.78170,   1.78160,   1.78190,   1.78300,
     &   1.78430,   1.78520,   1.78610,   1.78630,   1.78640,   1.78650,
     &   1.78650,   1.78650,   1.78650,   1.78650,   1.78650,   1.78650,
     &   1.78650,   1.78650,   1.78650,   1.78650,   1.78650,   1.78650,
     &   1.78650,   1.78650,   1.78650,   1.78650,   1.78650,   1.78650,
     &   1.78650,   1.78650,   1.78650,   1.78650,   1.78650,   1.78650,
     &   1.78650,   1.78650,   1.78650,   1.78650,   1.78650,   1.78650,
     &   1.78650,   1.78650,   1.78650,   1.78720/
      DATA ( TABRET(I,3), I = 1,NWLT ) /
     &              1.82961,   1.83258,   1.83149,   1.82748,   1.82224,
     &   1.81718,   1.81204,   1.80704,   1.80250,   1.79834,   1.79482,
     &   1.79214,   1.78843,   1.78601,   1.78434,   1.78322,   1.78248,
     &   1.78201,   1.78160,   1.78140,   1.78160,   1.78220,   1.78310,
     &   1.78380,   1.78390,   1.78400,   1.78400,   1.78400,   1.78400,
     &   1.78400,   1.78390,   1.78380,   1.78370,   1.78370,   1.78370,
     &   1.78370,   1.78370,   1.78370,   1.78370,   1.78370,   1.78370,
     &   1.78370,   1.78370,   1.78370,   1.78370,   1.78370,   1.78370,
     &   1.78370,   1.78370,   1.78370,   1.78370,   1.78370,   1.78370,
     &   1.78370,   1.78370,   1.78370,   1.78370,   1.78370,   1.78370,
     &   1.78370,   1.78400,   1.78450/
      DATA ( TABRET(I,4), I = 1,NWLT ) /
     &   1.82961,   1.83258,   1.83149,   1.82748,   1.82224,   1.81718,
     &   1.81204,   1.80704,   1.80250,   1.79834,   1.79482,   1.79214,
     &   1.78843,   1.78601,   1.78434,   1.78322,   1.78248,   1.78201,
     &   1.78150,   1.78070,   1.78010,   1.77890,   1.77790,   1.77730,
     &   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,
     &   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,
     &   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,
     &   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,
     &   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,
     &   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,   1.77720,
     &   1.77720,   1.77800/

      DATA ( TABIM(I), I = 1,114 ) /
     &0.1640E+00,0.1730E+00,0.1830E+00,0.1950E+00,0.2080E+00,0.2230E+00,
     &0.2400E+00,0.2500E+00,0.2590E+00,0.2680E+00,0.2790E+00,0.2970E+00,
     &0.3190E+00,0.3400E+00,0.3660E+00,0.3920E+00,0.4160E+00,0.4400E+00,
     &0.4640E+00,0.4920E+00,0.5170E+00,0.5280E+00,0.5330E+00,0.5340E+00,
     &0.5310E+00,0.5240E+00,0.5100E+00,0.5000E+00,0.4990E+00,0.4680E+00,
     &0.3800E+00,0.3600E+00,0.3390E+00,0.3180E+00,0.2910E+00,0.2510E+00,
     &0.2440E+00,0.2390E+00,0.2390E+00,0.2440E+00,0.2470E+00,0.2240E+00,
     &0.1950E+00,0.1740E+00,0.1720E+00,0.1800E+00,0.1940E+00,0.2130E+00,
     &0.2430E+00,0.2710E+00,0.2890E+00,0.3340E+00,0.3440E+00,0.3820E+00,
     &0.4010E+00,0.4065E+00,0.4050E+00,0.3890E+00,0.3770E+00,0.3450E+00,
     &0.3320E+00,0.3150E+00,0.2980E+00,0.2740E+00,0.2280E+00,0.1980E+00,
     &0.1720E+00,0.1560E+00,0.1100E+00,0.8300E-01,0.5800E-01,0.2200E-01,
     &0.1000E-01,0.3000E-02,0.1000E-02,0.3000E-03,0.1000E-03,0.3000E-04,
     &0.1000E-04,0.3000E-05,0.1000E-05,0.7000E-06,0.4000E-06,0.2000E-06,
     &0.1000E-06,0.6377E-07,0.3750E-07,0.2800E-07,0.2400E-07,0.2200E-07,
     &0.1900E-07,0.1750E-07,0.1640E-07,0.1590E-07,0.1325E-07,0.8623E-08,
     &0.5504E-08,0.3765E-08,0.2710E-08,0.2510E-08,0.2260E-08,0.2080E-08,
     &0.1910E-08,0.1540E-08,0.1530E-08,0.1550E-08,0.1640E-08,0.1780E-08,
     &0.1910E-08,0.2140E-08,0.2260E-08,0.2540E-08,0.2930E-08,0.3110E-08/
      DATA ( TABIM(I), I = 115,228 ) /
     &0.3290E-08,0.3520E-08,0.4040E-08,0.4880E-08,0.5730E-08,0.6890E-08,
     &0.8580E-08,0.1040E-07,0.1220E-07,0.1430E-07,0.1660E-07,0.1890E-07,
     &0.2090E-07,0.2400E-07,0.2900E-07,0.3440E-07,0.4030E-07,0.4300E-07,
     &0.4920E-07,0.5870E-07,0.7080E-07,0.8580E-07,0.1020E-06,0.1180E-06,
     &0.1340E-06,0.1400E-06,0.1430E-06,0.1450E-06,0.1510E-06,0.1830E-06,
     &0.2150E-06,0.2650E-06,0.3350E-06,0.3920E-06,0.4200E-06,0.4440E-06,
     &0.4740E-06,0.5110E-06,0.5530E-06,0.6020E-06,0.7550E-06,0.9260E-06,
     &0.1120E-05,0.1330E-05,0.1620E-05,0.2000E-05,0.2250E-05,0.2330E-05,
     &0.2330E-05,0.2170E-05,0.1960E-05,0.1810E-05,0.1740E-05,0.1730E-05,
     &0.1700E-05,0.1760E-05,0.1820E-05,0.2040E-05,0.2250E-05,0.2290E-05,
     &0.3040E-05,0.3840E-05,0.4770E-05,0.5760E-05,0.6710E-05,0.8660E-05,
     &0.1020E-04,0.1130E-04,0.1220E-04,0.1290E-04,0.1320E-04,0.1350E-04,
     &0.1330E-04,0.1320E-04,0.1320E-04,0.1310E-04,0.1320E-04,0.1320E-04,
     &0.1340E-04,0.1390E-04,0.1420E-04,0.1480E-04,0.1580E-04,0.1740E-04,
     &0.1980E-04,0.2500E-04,0.5400E-04,0.1040E-03,0.2030E-03,0.2708E-03,
     &0.3511E-03,0.4299E-03,0.5181E-03,0.5855E-03,0.5899E-03,0.5635E-03,
     &0.5480E-03,0.5266E-03,0.4394E-03,0.3701E-03,0.3372E-03,0.2410E-03,
     &0.1890E-03,0.1660E-03,0.1450E-03,0.1280E-03,0.1030E-03,0.8600E-04,
     &0.8220E-04,0.8030E-04,0.8500E-04,0.9900E-04,0.1500E-03,0.2950E-03/
      DATA ( TABIM(I), I = 229,342 ) /
     &0.4687E-03,0.7615E-03,0.1010E-02,0.1313E-02,0.1539E-02,0.1588E-02,
     &0.1540E-02,0.1412E-02,0.1244E-02,0.1068E-02,0.8414E-03,0.5650E-03,
     &0.4320E-03,0.3500E-03,0.2870E-03,0.2210E-03,0.2030E-03,0.2010E-03,
     &0.2030E-03,0.2140E-03,0.2320E-03,0.2890E-03,0.3810E-03,0.4620E-03,
     &0.5480E-03,0.6180E-03,0.6800E-03,0.7300E-03,0.7820E-03,0.8480E-03,
     &0.9250E-03,0.9200E-03,0.8920E-03,0.8700E-03,0.8900E-03,0.9300E-03,
     &0.1010E-02,0.1350E-02,0.3420E-02,0.7920E-02,0.2000E-01,0.3800E-01,
     &0.5200E-01,0.6800E-01,0.9230E-01,0.1270E+00,0.1690E+00,0.2210E+00,
     &0.2760E+00,0.3120E+00,0.3470E+00,0.3880E+00,0.4380E+00,0.4930E+00,
     &0.5540E+00,0.6120E+00,0.6250E+00,0.5930E+00,0.5390E+00,0.4910E+00,
     &0.4380E+00,0.3720E+00,0.3000E+00,0.2380E+00,0.1930E+00,0.1580E+00,
     &0.1210E+00,0.1030E+00,0.8360E-01,0.6680E-01,0.5400E-01,0.4220E-01,
     &0.3420E-01,0.2740E-01,0.2200E-01,0.1860E-01,0.1520E-01,0.1260E-01,
     &0.1060E-01,0.8020E-02,0.6850E-02,0.6600E-02,0.6960E-02,0.9160E-02,
     &0.1110E-01,0.1450E-01,0.2000E-01,0.2300E-01,0.2600E-01,0.2900E-01,
     &0.2930E-01,0.3000E-01,0.2850E-01,0.1730E-01,0.1290E-01,0.1200E-01,
     &0.1250E-01,0.1340E-01,0.1400E-01,0.1750E-01,0.2400E-01,0.3500E-01,
     &0.3800E-01,0.4200E-01,0.4600E-01,0.5200E-01,0.5700E-01,0.6900E-01,
     &0.7000E-01,0.6700E-01,0.6500E-01,0.6400E-01,0.6200E-01,0.5900E-01/
      DATA ( TABIM(I), I = 343,456 ) /
     &0.5700E-01,0.5600E-01,0.5500E-01,0.5700E-01,0.5800E-01,0.5700E-01,
     &0.5500E-01,0.5500E-01,0.5400E-01,0.5200E-01,0.5200E-01,0.5200E-01,
     &0.5200E-01,0.5000E-01,0.4700E-01,0.4300E-01,0.3900E-01,0.3700E-01,
     &0.3900E-01,0.4000E-01,0.4200E-01,0.4400E-01,0.4500E-01,0.4600E-01,
     &0.4700E-01,0.5100E-01,0.6500E-01,0.7500E-01,0.8800E-01,0.1080E+00,
     &0.1340E+00,0.1680E+00,0.2040E+00,0.2480E+00,0.2800E+00,0.3410E+00,
     &0.3790E+00,0.4090E+00,0.4220E+00,0.4220E+00,0.4030E+00,0.3890E+00,
     &0.3740E+00,0.3540E+00,0.3350E+00,0.3150E+00,0.2940E+00,0.2710E+00,
     &0.2460E+00,0.1980E+00,0.1640E+00,0.1520E+00,0.1420E+00,0.1280E+00,
     &0.1250E+00,0.1230E+00,0.1160E+00,0.1070E+00,0.7900E-01,0.7200E-01,
     &0.7600E-01,0.7500E-01,0.6700E-01,0.5500E-01,0.4500E-01,0.2900E-01,
     &0.2750E-01,0.2700E-01,0.2730E-01,0.2890E-01,0.3000E-01,0.3400E-01,
     &0.5300E-01,0.7550E-01,0.1060E+00,0.1350E+00,0.1761E+00,0.2229E+00,
     &0.2746E+00,0.3280E+00,0.3906E+00,0.4642E+00,0.5247E+00,0.5731E+00,
     &0.6362E+00,0.6839E+00,0.7091E+00,0.6790E+00,0.6250E+00,0.5654E+00,
     &0.5433E+00,0.5292E+00,0.5070E+00,0.4883E+00,0.4707E+00,0.4203E+00,
     &0.3771E+00,0.3376E+00,0.3056E+00,0.2835E+00,0.3170E+00,0.3517E+00,
     &0.3902E+00,0.4509E+00,0.4671E+00,0.4779E+00,0.4890E+00,0.4899E+00,
     &0.4873E+00,0.4766E+00,0.4508E+00,0.4193E+00,0.3880E+00,0.3433E+00/
      DATA ( TABIM(I), I = 457,468 ) /
     &0.3118E+00,0.2935E+00,0.2350E+00,0.1981E+00,0.1865E+00,0.1771E+00,
     &0.1620E+00,0.1490E+00,0.1390E+00,0.1200E+00,0.9620E-01,0.8300E-01/

      DATA ( TABIMT(I,1), I = 1,NWLT ) /
     &                                 0.8300E-01,0.6900E-01,0.5700E-01,
     &0.4560E-01,0.3790E-01,0.3140E-01,0.2620E-01,0.2240E-01,0.1960E-01,
     &0.1760E-01,0.1665E-01,0.1620E-01,0.1550E-01,0.1470E-01,0.1390E-01,
     &0.1320E-01,0.1250E-01,0.1180E-01,0.1060E-01,0.9540E-02,0.8560E-02,
     &0.6210E-02,0.4490E-02,0.3240E-02,0.2340E-02,0.1880E-02,0.1740E-02,
     &0.1500E-02,0.1320E-02,0.1160E-02,0.8800E-03,0.6950E-03,0.4640E-03,
     &0.3400E-03,0.3110E-03,0.2940E-03,0.2790E-03,0.2700E-03,0.2640E-03,
     &0.2580E-03,0.2520E-03,0.2490E-03,0.2540E-03,0.2640E-03,0.2740E-03,
     &0.2890E-03,0.3050E-03,0.3150E-03,0.3460E-03,0.3820E-03,0.4620E-03,
     &0.5000E-03,0.5500E-03,0.5950E-03,0.6470E-03,0.6920E-03,0.7420E-03,
     &0.8200E-03,0.9700E-03,0.1950E-02,0.5780E-02,0.9700E-02/
      DATA ( TABIMT(I,2), I = 1,NWLT ) /
     &                      0.8300E-01,0.6900E-01,0.5700E-01,0.4560E-01,
     &0.3790E-01,0.3140E-01,0.2620E-01,0.2240E-01,0.1960E-01,0.1760E-01,
     &0.1665E-01,0.1600E-01,0.1500E-01,0.1400E-01,0.1310E-01,0.1230E-01,
     &0.1150E-01,0.1080E-01,0.9460E-02,0.8290E-02,0.7270E-02,0.4910E-02,
     &0.3300E-02,0.2220E-02,0.1490E-02,0.1140E-02,0.1060E-02,0.9480E-03,
     &0.8500E-03,0.7660E-03,0.6300E-03,0.5200E-03,0.3840E-03,0.2960E-03,
     &0.2700E-03,0.2520E-03,0.2440E-03,0.2360E-03,0.2300E-03,0.2280E-03,
     &0.2250E-03,0.2200E-03,0.2160E-03,0.2170E-03,0.2200E-03,0.2250E-03,
     &0.2320E-03,0.2390E-03,0.2600E-03,0.2860E-03,0.3560E-03,0.3830E-03,
     &0.4150E-03,0.4450E-03,0.4760E-03,0.5080E-03,0.5400E-03,0.5860E-03,
     &0.6780E-03,0.1280E-02,0.3550E-02,0.5600E-02/
      DATA ( TABIMT(I,3), I = 1,NWLT ) /
     &           0.8300E-01,0.6900E-01,0.5700E-01,0.4560E-01,0.3790E-01,
     &0.3140E-01,0.2620E-01,0.2190E-01,0.1880E-01,0.1660E-01,0.1540E-01,
     &0.1470E-01,0.1350E-01,0.1250E-01,0.1150E-01,0.1060E-01,0.9770E-02,
     &0.9010E-02,0.7660E-02,0.6520E-02,0.5540E-02,0.3420E-02,0.2100E-02,
     &0.1290E-02,0.7930E-03,0.5700E-03,0.5350E-03,0.4820E-03,0.4380E-03,
     &0.4080E-03,0.3500E-03,0.3200E-03,0.2550E-03,0.2120E-03,0.2000E-03,
     &0.1860E-03,0.1750E-03,0.1660E-03,0.1560E-03,0.1490E-03,0.1440E-03,
     &0.1350E-03,0.1210E-03,0.1160E-03,0.1160E-03,0.1170E-03,0.1200E-03,
     &0.1230E-03,0.1320E-03,0.1440E-03,0.1680E-03,0.1800E-03,0.1900E-03,
     &0.2090E-03,0.2160E-03,0.2290E-03,0.2400E-03,0.2600E-03,0.2920E-03,
     &0.6100E-03,0.1020E-02,0.1810E-02/
      DATA ( TABIMT(I,4), I = 1,NWLT ) /
     &0.8300E-01,0.6900E-01,0.5700E-01,0.4450E-01,0.3550E-01,0.2910E-01,
     &0.2440E-01,0.1970E-01,0.1670E-01,0.1400E-01,0.1235E-01,0.1080E-01,
     &0.8900E-02,0.7340E-02,0.6400E-02,0.5600E-02,0.5000E-02,0.4520E-02,
     &0.3680E-02,0.2990E-02,0.2490E-02,0.1550E-02,0.9610E-03,0.5950E-03,
     &0.3690E-03,0.2670E-03,0.2510E-03,0.2290E-03,0.2110E-03,0.1960E-03,
     &0.1730E-03,0.1550E-03,0.1310E-03,0.1130E-03,0.1060E-03,0.9900E-04,
     &0.9300E-04,0.8730E-04,0.8300E-04,0.7870E-04,0.7500E-04,0.6830E-04,
     &0.5600E-04,0.4960E-04,0.4550E-04,0.4210E-04,0.3910E-04,0.3760E-04,
     &0.3400E-04,0.3100E-04,0.2640E-04,0.2510E-04,0.2430E-04,0.2390E-04,
     &0.2370E-04,0.2380E-04,0.2400E-04,0.2460E-04,0.2660E-04,0.4450E-04,
     &0.8700E-04,0.1320E-03/
c     ..

      END
