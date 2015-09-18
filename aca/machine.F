      REAL FUNCTION R1MACH(I)

c        Single-precision machine constants

c  Assume floating-point numbers are represented in the t-digit,
c  base-b form

c         sign (b**e)*( (x(1)/b) + ... + (x(t)/b**t) )

c  where 0.le.x(i).lt.b  for  i = 1,...,t,
c  0.lt.x(1), and  emin.LE.e.LE.emax.  then

c  R1MACH(1) = b**(emin-1), the smallest positive magnitude
c              (use TINY(R) in Fortran 90, where R is a single
c              precision variable)

c  R1MACH(2) = b**emax*(1 - b**(-t)), the largest magnitude
c              (use HUGE(R) in Fortran 90, where R is a single
c              precision variable))

c  R1MACH(3) = b**(-t), the smallest relative spacing.

c  R1MACH(4) = b**(1-t), the largest relative spacing.  i.e.,
c              smallest positive eps such that  1+eps .ne. 1
c              (use EPSILON(R) in Fortran 90, where R is a single
c              precision variable))

c  R1MACH(5) = LOG10(b)

c  Reference: Fox P.A., Hall A.D., Schryer N.L.,'Framework For A
c               Portable Library', ACM Transactions On Mathematical
c               Software, Vol. 4, No. 2, June 1978, pp. 177-188.

c  To alter this function for a particular environment,
c  the desired set of data statements should be activated by
c  deleting the C from column 1.

c  For IEEE-arithmetic machines (binary standard), the first
c  set of constants below should be appropriate.

c  Where possible, decimal, octal or hexadecimal constants are used
c  to specify the constants exactly.  Sometimes this requires using
c  equivalent integer arrays.  If your compiler uses half-word
c  integers by default (sometimes called INTEGER*2), you may need to
c  change INTEGER to INTEGER*4 or otherwise instruct your compiler
c  to use full-word integers in the next 5 declarations.

      INTEGER SMALL(2), LARGE(2), RIGHT(2), DIVER(2), LOG10(2), SC
      REAL RMACH(5)

      EQUIVALENCE (RMACH(1),SMALL(1)), (RMACH(2),LARGE(1)),
     $            (RMACH(3),RIGHT(1)), (RMACH(4),DIVER(1)),
     $            (RMACH(5),LOG10(1))

      LOGICAL  PASS1
      SAVE     PASS1
      DATA     PASS1/.TRUE./

c IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T
c 3B SERIES, MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T
c PC 7300), AND 8087 BASED MICROS (E.G. IBM PC AND AT&T 6300).

c++csz
      DATA SMALL(1)/8388608/, LARGE(1)/2139095039/,
     $     RIGHT(1)/864026624/, DIVER(1)/872415232/,
     $     LOG10(1)/ 1050288283/, SC/987/
c--csz

c AMDAHL MACHINES.

c      DATA SMALL(1)/1048576/, LARGE(1)/2147483647/,
c     $     RIGHT(1)/990904320/, DIVER(1)/1007681536/,
c     $     LOG10(1)/1091781651/, SC/987/

c BURROUGHS 1700 SYSTEM.

c      DATA RMACH/Z400800000,Z5FFFFFFFF,Z4E9800000,Z4EA800000,
c     $             Z500E730E8/, SC/987/

c BURROUGHS 5700/6700/7700 SYSTEMS.

c      DATA RMACH/O1771000000000000,O0777777777777777,O1311000000000000,
c     $             O1301000000000000,O1157163034761675/, SC/987/

c FTN4 ON CDC 6000/7000 SERIES.

c      DATA RMACH/00564000000000000000B,37767777777777777776B,
c     $ 16414000000000000000B,16424000000000000000B,
c     $ 17164642023241175720B/, SC/987/

c FTN5 ON CDC 6000/7000 SERIES.

c      DATA RMACH/O"00564000000000000000",O"37767777777777777776",
c     $ O"16414000000000000000",O"16424000000000000000",
c     $ O"17164642023241175720"/, SC/987/

c CONVEX C-1.

c      DATA RMACH/'00800000'X,'7FFFFFFF'X,'34800000'X,
c     $ '35000000'X,'3F9A209B'X/, SC/987/

c CRAY 1, XMP, 2, AND 3.

c      DATA RMACH/200034000000000000000B,577767777777777777776B,
c     $ 377224000000000000000B,377234000000000000000B,
c     $ 377774642023241175720B/, SC/987/

c DATA GENERAL ECLIPSE S/200.
c NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING LINE -
c STATIC RMACH(5)

c      DATA SMALL/20K,0/, LARGE/77777K,177777K/, RIGHT/35420K,0/,
c     $  DIVER/36020K,0/, LOG10/40423K,42023K/, SC/987/

c HARRIS SLASH 6 AND SLASH 7.

c      DATA SMALL/'20000000,'00000201/, LARGE/'37777777,'00000177/,
c     $  RIGHT/'20000000,'00000352/, DIVER/'20000000,'00000353/,
c     $  LOG10/'23210115,'00000377/, SC/987/

c HONEYWELL DPS 8/70 SERIES.

c      DATA RMACH/O402400000000,O376777777777,O714400000000,
c     $ O716400000000,O776464202324/, SC/987/

c IBM 360/370 SERIES,
c XEROX SIGMA 5/7/9 AND SEL SYSTEMS 85/86.

c      DATA RMACH/Z00100000,Z7FFFFFFF,Z3B100000,Z3C100000,
c     $ Z41134413/, SC/987/

c INTERDATA 8/32 WITH UNIX SYSTEM FORTRAN 77 COMPILER.
c FOR INTERDATA FORTRAN VII COMPILER REPLACE
c Z'S SPECIFYING HEX CONSTANTS WITH Y'S.

c      DATA RMACH/Z'00100000',Z'7EFFFFFF',Z'3B100000',Z'3C100000',
c     $ Z'41134413'/, SC/987/

c PDP-10 (KA OR KI PROCESSOR).

c      DATA RMACH/"000400000000,"377777777777,"146400000000,
c     $ "147400000000,"177464202324/, SC/987/

c PDP-11 FORTRANS SUPPORTING 32-BIT INTEGERS
c (EXPRESSED IN INTEGER AND OCTAL).

c      DATA SMALL(1)/8388608/, LARGE(1)/2147483647/,
c     $  RIGHT(1)/880803840/, DIVER(1)/889192448/,
c     $  LOG10(1)/1067065499/, SC/987/

c      DATA RMACH/O00040000000,O17777777777,O06440000000,
c     $ O06500000000,O07746420233/, SC/987/

c PDP-11 FORTRANS SUPPORTING 16-BIT INTEGERS
c (EXPRESSED IN INTEGER AND OCTAL).

c      DATA SMALL/128,0/, LARGE/32767,-1/, RIGHT/13440,0/,
c     $  DIVER/13568,0/, LOG10/16282,8347/, SC/987/

c      DATA SMALL/O000200,O000000/, LARGE/O077777,O177777/,
c     $  RIGHT/O032200,O000000/, DIVER/O032400,O000000/,
c     $  LOG10/O037632,O020233/, SC/987/

c SEQUENT BALANCE 8000.

c      DATA SMALL(1)/$00800000/, LARGE(1)/$7F7FFFFF/,
c     $  RIGHT(1)/$33800000/, DIVER(1)/$34000000/,
c     $  LOG10(1)/$3E9A209B/, SC/987/

c UNIVAC 1100 SERIES.

c      DATA RMACH/O000400000000,O377777777777,O146400000000,
c     $ O147400000000,O177464202324/, SC/987/

c VAX UNIX F77 COMPILER.

c      DATA SMALL(1)/128/, LARGE(1)/-32769/, RIGHT(1)/13440/,
c     $  DIVER(1)/13568/, LOG10(1)/547045274/, SC/987/

c VAX-11 WITH FORTRAN IV-PLUS COMPILER.

c      DATA RMACH/Z00000080,ZFFFF7FFF,Z00003480,Z00003500,
c     $ Z209B3F9A/, SC/987/

c VAX/VMS VERSION 2.2.

c      DATA RMACH/'80'X,'FFFF7FFF'X,'3480'X,'3500'X,
c     $ '209B3F9A'X/, SC/987/


      IF( PASS1 )  THEN

         PASS1 = .FALSE.
         IF (SC.NE.987)
     $       CALL ErrMsg('R1MACH--no DATA statements active',.TRUE.)

c                      ** Calculate machine precision
         EPSNEW = 0.01
   10    EPS = EPSNEW
            EPSNEW = EPSNEW / 1.1
c                                 ** This may force 'S' to be stored
c                                    but there is no guarantee;  if it
c                                    is kept in a register, it may be
c                                    kept in higher precision
            S = 1.0 + EPSNEW
            IF( S.GT.1.0 ) GO TO 10

         RATIO = EPS / RMACH(4)

c++csz
         IF( RATIO.LT.0.5 .OR. RATIO.GT.2.0 ) 
     $        then
         write (6,'(a,e15.8)') 'r1mach() reports ratio = ',ratio
         write (6,'(a)') 'HINT: LINUX g77 fails here only when optimization is invoked'
         CALL ErrMsg('R1MACH--tabulated precision wrong',.TRUE.)
      endif                     ! endif
c--csz 

      END IF

      IF (I.LT.1.OR.I.GT.5)
     $    CALL ErrMsg('R1MACH--argument out of bounds',.TRUE.)
      R1MACH = RMACH(I)
      RETURN
      END
      DOUBLE PRECISION FUNCTION D1MACH(I)

c  Double-precision machine constants (see R1MACH for documentation)

c  For IEEE-arithmetic machines (binary standard), one of the first
c  two sets of constants below should be appropriate.

      INTEGER SMALL(4), LARGE(4), RIGHT(4), DIVER(4), LOG10(4), SC
      DOUBLE PRECISION DMACH(5), EPS, EPSNEW, S

      EQUIVALENCE (DMACH(1),SMALL(1)), (DMACH(2),LARGE(1)),
     $            (DMACH(3),RIGHT(1)), (DMACH(4),DIVER(1)),
     $            (DMACH(5),LOG10(1))

      LOGICAL  PASS1
      SAVE     PASS1
      DATA     PASS1/.TRUE./

c IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T 3B SERIES AND
c MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T PC 7300),
c IN WHICH THE MOST SIGNIFICANT BYTE IS STORED FIRST.

c++csz
#ifndef LINUX
      DATA (SMALL(N),N=1,2)/1048576,0/, (LARGE(N),N=1,2)/2146435071,-1/,
     $  (RIGHT(N),N=1,2)/1017118720,0/, (DIVER(N),N=1,2)/1018167296,0/,
     $  (LOG10(N),N=1,2)/1070810131,1352628735/, SC/987/
#endif /* LINUX */ 
c--csz

c IEEE ARITHMETIC MACHINES AND 8087-BASED MICROS, SUCH AS THE IBM PC
c AND AT&T 6300, IN WHICH THE LEAST SIGNIFICANT BYTE IS STORED FIRST.

c++csz
#ifdef LINUX
      DATA (SMALL(N),N=1,2)/0,1048576/, (LARGE(N),N=1,2)/-1,2146435071/,
     $  (RIGHT(N),N=1,2)/0,1017118720/, (DIVER(N),N=1,2)/0,1018167296/,
     $  (LOG10(N),N=1,2)/1352628735,1070810131/, SC/987/
#endif /* not LINUX */ 
c--csz

c AMDAHL MACHINES.

c      DATA (SMALL(N),N=1,2)/1048576,0/, (LARGE(N),N=1,2)/2147483647,-1/,
c     $ (RIGHT(N),N=1,2)/856686592,0/, (DIVER(N),N=1,2)/ 873463808,0/,
c     $ (LOG10(N),N=1,2)/1091781651,1352628735/, SC/987/

c BURROUGHS 1700 SYSTEM.

c      DATA (SMALL(N),N=1,2)/ZC00800000,Z000000000/,
c     $ (LARGE(N),N=1,2)/ZDFFFFFFFF,ZFFFFFFFFF/,
c     $ (RIGHT(N),N=1,2)/ZCC5800000,Z000000000/,
c     $ (DIVER(N),N=1,2)/ZCC6800000,Z000000000/,
c     $ (LOG10(N),N=1,2)/ZD00E730E7,ZC77800DC0/, SC/987/

c BURROUGHS 5700 SYSTEM.

c      DATA (SMALL(N),N=1,2)/O1771000000000000,O0000000000000000/,
c     $  (LARGE(N),N=1,2)/O0777777777777777,O0007777777777777/,
c     $  (RIGHT(N),N=1,2)/O1461000000000000,O0000000000000000/,
c     $  (DIVER(N),N=1,2)/O1451000000000000,O0000000000000000/,
c     $  (LOG10(N),N=1,2)/O1157163034761674,O0006677466732724/, SC/987/

c BURROUGHS 6700/7700 SYSTEMS.

c      DATA (SMALL(N),N=1,2)/O1771000000000000,O7770000000000000/,
c     $  (LARGE(N),N=1,2)/O0777777777777777,O7777777777777777/,
c     $  (RIGHT(N),N=1,2)/O1461000000000000,O0000000000000000/,
c     $  (DIVER(N),N=1,2)/O1451000000000000,O0000000000000000/,
c     $  (LOG10(N),N=1,2)/O1157163034761674,O0006677466732724/, SC/987/

c FTN4 ON THE CDC 6000/7000 SERIES.

c      DATA
c     $  (SMALL(N),N=1,2)/00564000000000000000B,00000000000000000000B/,
c     $  (LARGE(N),N=1,2)/37757777777777777777B,37157777777777777774B/,
c     $  (RIGHT(N),N=1,2)/15624000000000000000B,00000000000000000000B/,
c     $  (DIVER(N),N=1,2)/15634000000000000000B,00000000000000000000B/,
c     $  (LOG10(N),N=1,2)/17164642023241175717B,16367571421742254654B/,
c     $  SC/987/

c FTN5 ON THE CDC 6000/7000 SERIES.

c      DATA
c     $(SMALL(N),N=1,2)/O"00564000000000000000",O"00000000000000000000"/,
c     $(LARGE(N),N=1,2)/O"37757777777777777777",O"37157777777777777774"/,
c     $(RIGHT(N),N=1,2)/O"15624000000000000000",O"00000000000000000000"/,
c     $(DIVER(N),N=1,2)/O"15634000000000000000",O"00000000000000000000"/,
c     $(LOG10(N),N=1,2)/O"17164642023241175717",O"16367571421742254654"/,
c     $ SC/987/

c CONVEX C-1

c      DATA (SMALL(N),N=1,2)/'00100000'X,'00000000'X/,
c     $  (LARGE(N),N=1,2)/'7FFFFFFF'X,'FFFFFFFF'X/,
c     $  (RIGHT(N),N=1,2)/'3CC00000'X,'00000000'X/,
c     $  (DIVER(N),N=1,2)/'3CD00000'X,'00000000'X/,
c     $  (LOG10(N),N=1,2)/'3FF34413'X,'509F79FF'X/, SC/987/

c CRAY 1, XMP, 2, AND 3.

c      DATA
c     $ (SMALL(N),N=1,2)/201354000000000000000B,000000000000000000000B/,
c     $ (LARGE(N),N=1,2)/577767777777777777777B,000007777777777777776B/,
c     $ (RIGHT(N),N=1,2)/376434000000000000000B,000000000000000000000B/,
c     $ (DIVER(N),N=1,2)/376444000000000000000B,000000000000000000000B/,
c     $ (LOG10(N),N=1,2)/377774642023241175717B,000007571421742254654B/,
c     $ SC/987/

c DATA GENERAL ECLIPSE S/200
c NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING LINE -
c STATIC DMACH(5)

c      DATA SMALL/20K,3*0/, LARGE/77777K,3*177777K/,
c     $  RIGHT/31420K,3*0/, DIVER/32020K,3*0/,
c     $  LOG10/40423K,42023K,50237K,74776K/, SC/987/

c HARRIS SLASH 6 AND SLASH 7

c      DATA (SMALL(N),N=1,2)/'20000000,'00000201/,
c     $  (LARGE(N),N=1,2)/'37777777,'37777577/,
c     $  (RIGHT(N),N=1,2)/'20000000,'00000333/,
c     $  (DIVER(N),N=1,2)/'20000000,'00000334/,
c     $  (LOG10(N),N=1,2)/'23210115,'10237777/, SC/987/

c HONEYWELL DPS 8/70 SERIES.

c      DATA (SMALL(N),N=1,2)/O402400000000,O000000000000/,
c     $  (LARGE(N),N=1,2)/O376777777777,O777777777777/,
c     $  (RIGHT(N),N=1,2)/O604400000000,O000000000000/,
c     $  (DIVER(N),N=1,2)/O606400000000,O000000000000/,
c     $  (LOG10(N),N=1,2)/O776464202324,O117571775714/, SC/987/

c IBM 360/370 SERIES, XEROX SIGMA 5/7/9 AND THE SEL SYSTEMS 85/86.

c      DATA (SMALL(N),N=1,2)/Z00100000,Z00000000/,
c     $  (LARGE(N),N=1,2)/Z7FFFFFFF,ZFFFFFFFF/,
c     $  (RIGHT(N),N=1,2)/Z33100000,Z00000000/,
c     $  (DIVER(N),N=1,2)/Z34100000,Z00000000/,
c     $  (LOG10(N),N=1,2)/Z41134413,Z509F79FF/, SC/987/

c INTERDATA 8/32 WITH THE UNIX SYSTEM FORTRAN 77 COMPILER.
c FOR THE INTERDATA FORTRAN VII COMPILER REPLACE
c THE Z'S SPECIFYING HEX CONSTANTS WITH Y'S.

c      DATA (SMALL(N),N=1,2)/Z'00100000',Z'00000000'/,
c     $  (LARGE(N),N=1,2)/Z'7EFFFFFF',Z'FFFFFFFF'/,
c     $  (RIGHT(N),N=1,2)/Z'33100000',Z'00000000'/,
c     $  (DIVER(N),N=1,2)/Z'34100000',Z'00000000'/,
c     $  (LOG10(N),N=1,2)/Z'41134413',Z'509F79FF'/, SC/987/

c PDP-10 (KA PROCESSOR).

c      DATA (SMALL(N),N=1,2)/"033400000000,"000000000000/,
c     $  (LARGE(N),N=1,2)/"377777777777,"344777777777/,
c     $  (RIGHT(N),N=1,2)/"113400000000,"000000000000/,
c     $  (DIVER(N),N=1,2)/"114400000000,"000000000000/,
c     $  (LOG10(N),N=1,2)/"177464202324,"144117571776/, SC/987/

c PDP-10 (KI PROCESSOR).

c      DATA (SMALL(N),N=1,2)/"000400000000,"000000000000/,
c     $  (LARGE(N),N=1,2)/"377777777777,"377777777777/,
c     $  (RIGHT(N),N=1,2)/"103400000000,"000000000000/,
c     $  (DIVER(N),N=1,2)/"104400000000,"000000000000/,
c     $  (LOG10(N),N=1,2)/"177464202324,"047674776746/, SC/987/

c PDP-11 FORTRANS SUPPORTING 32-BIT INTEGERS
c (EXPRESSED IN INTEGER AND OCTAL).

c      DATA (SMALL(N),N=1,2)/8388608,0/, (LARGE(N),N=1,2)/2147483647,-1/,
c     $  (RIGHT(N),N=1,2)/612368384,0/, (DIVER(N),N=1,2)/620756992,0/,
c     $  (LOG10(N),N=1,2)/1067065498,-2063872008/, SC/987/

c      DATA (SMALL(N),N=1,2)/O00040000000,O00000000000/,
c     $  (LARGE(N),N=1,2)/O17777777777,O37777777777/,
c     $  (RIGHT(N),N=1,2)/O04440000000,O00000000000/,
c     $  (DIVER(N),N=1,2)/O04500000000,O00000000000/,
c     $  (LOG10(N),N=1,2)/O07746420232,O20476747770/, SC/987/

c PDP-11 FORTRANS SUPPORTING 16-BIT INTEGERS
c (EXPRESSED IN INTEGER AND OCTAL).

c      DATA SMALL/128,3*0/, LARGE/32767,3*-1/, RIGHT/9344,3*0/,
c     $  DIVER/9472,3*0/, LOG10/16282,8346,-31493,-12296/, SC/987/

c      DATA SMALL/O000200,3*O000000/, LARGE/O077777,3*O177777/,
c     $  RIGHT/O022200,3*O000000/, DIVER/O022400,3*O000000/,
c     $  LOG10/O037632,O020232,O102373,O147770/, SC/987/

c PRIME 50 SERIES SYSTEMS WITH 32-BIT INTEGERS AND 64V MODE
c INSTRUCTIONS, SUPPLIED BY IGOR BRAY.

c      DATA (SMALL(N),N=1,2)/:10000000000,:00000100001/,
c     $  (LARGE(N),N=1,2)/:17777777777,:37777677775/,
c     $  (RIGHT(N),N=1,2)/:10000000000,:00000000122/,
c     $  (DIVER(N),N=1,2)/:10000000000,:00000000123/,
c     $  (LOG10(N),N=1,2)/:11504046501,:07674600177/, SC/987/

c SEQUENT BALANCE 8000

c      DATA (SMALL(N),N=1,2)/$00000000, $00100000/,
c     $  (LARGE(N),N=1,2)/$FFFFFFFF, $7FEFFFFF/,
c     $  (RIGHT(N),N=1,2)/$00000000, $3CA00000/,
c     $  (DIVER(N),N=1,2)/$00000000, $3CB00000/,
c     $  (LOG10(N),N=1,2)/$509F79FF, $3FD34413/, SC/987/

c UNIVAC 1100 SERIES.

c      DATA (SMALL(N),N=1,2)/O000040000000,O000000000000/,
c     $  (LARGE(N),N=1,2)/O377777777777,O777777777777/,
c     $  (RIGHT(N),N=1,2)/O170540000000,O000000000000/,
c     $  (DIVER(N),N=1,2)/O170640000000,O000000000000/,
c     $  (LOG10(N),N=1,2)/O177746420232,O411757177572/, SC/987/

c VAX UNIX F77 COMPILER

c      DATA (SMALL(N),N=1,2)/128,0/, (LARGE(N),N=1,2)/-32769,-1/,
c     $  (RIGHT(N),N=1,2)/9344,0/, (DIVER(N),N=1,2)/9472,0/,
c     $  (LOG10(N),N=1,2)/546979738,-805796613/, SC/987/

c VAX-11 WITH FORTRAN IV-PLUS COMPILER

c      DATA (SMALL(N),N=1,2)/Z00000080,Z00000000/,
c     $  (LARGE(N),N=1,2)/ZFFFF7FFF,ZFFFFFFFF/,
c     $  (RIGHT(N),N=1,2)/Z00002480,Z00000000/,
c     $  (DIVER(N),N=1,2)/Z00002500,Z00000000/,
c     $  (LOG10(N),N=1,2)/Z209A3F9A,ZCFF884FB/, SC/987/

c VAX/VMS VERSION 2.2

c      DATA (SMALL(N),N=1,2)/'80'X,'0'X/,
c     $  (LARGE(N),N=1,2)/'FFFF7FFF'X,'FFFFFFFF'X/,
c     $  (RIGHT(N),N=1,2)/'2480'X,'0'X/, (DIVER(N),N=1,2)/'2500'X,'0'X/,
c     $  (LOG10(N),N=1,2)/'209A3F9A'X,'CFF884FB'X/, SC/987/


      IF( PASS1 )  THEN

         PASS1 = .FALSE.
         IF (SC.NE.987)
     $       CALL ErrMsg( 'D1MACH--no DATA statements active',.TRUE.)

c                        ** Calculate machine precision
         EPSNEW = 0.01D0
   10    EPS = EPSNEW
            EPSNEW = EPSNEW / 1.1D0
c                                 ** This may force 'S' to be stored
c                                    but there is no guarantee;  if it
c                                    is kept in a register, it may be
c                                    kept in higher precision
            S = 1.D0 + EPSNEW
            IF( S.GT.1.D0 ) GO TO 10

         IF( EPS/DMACH(4).LT.0.5D0 .OR. EPS/DMACH(4).GT.2.D0 )
     $       CALL ErrMsg( 'D1MACH--tabulated precision wrong',.TRUE.)

      END IF

      IF (I.LT.1.OR.I.GT.5)
     $    CALL ERRMSG( 'D1MACH--argument out of bounds',.TRUE.)
      D1MACH = DMACH(I)
      RETURN
      END
      SUBROUTINE  ErrMsg( MESSAG, FATAL )

c        Print out a warning or error message;  abort if error
c        after making symbolic dump (machine-specific)

      LOGICAL       FATAL, MsgLim, Cray
      CHARACTER*(*) MESSAG
      INTEGER       MaxMsg, NumMsg
      SAVE          MaxMsg, NumMsg, MsgLim
      DATA NumMsg / 0 /,  MaxMsg / 100 /,  MsgLim / .FALSE. /


      IF ( FATAL )  THEN
         WRITE ( *, '(//,2A,//)' )  ' ******* ERROR >>>>>>  ', MESSAG
         STOP
      END IF

      NumMsg = NumMsg + 1
      IF( MsgLim )  RETURN

      IF ( NumMsg.LE.MaxMsg )  THEN
         WRITE ( *, '(/,2A,/)' )  ' ******* WARNING >>>>>>  ', MESSAG
      ELSE
         WRITE ( *,99 )
         MsgLim = .True.
      ENDIF

      RETURN

   99 FORMAT( //,' >>>>>>  TOO MANY WARNING MESSAGES --  ',
     $   'They will no longer be printed  <<<<<<<', // )
      END

      LOGICAL FUNCTION  WrtBad ( VarNam )

c          Write names of erroneous variables and return 'TRUE'

c      INPUT :   VarNam = Name of erroneous variable to be written
c                         ( CHARACTER, any length )

      CHARACTER*(*)  VarNam
      INTEGER        MaxMsg, NumMsg
      SAVE  NumMsg, MaxMsg
      DATA  NumMsg / 0 /,  MaxMsg / 50 /


      WrtBad = .TRUE.
      NumMsg = NumMsg + 1
      WRITE ( *, '(3A)' )  ' ****  Input variable  ', VarNam,
     $                     '  in error  ****'
      IF ( NumMsg.EQ.MaxMsg )
     $   CALL  ErrMsg ( 'Too many input errors.  Aborting...', .TRUE. )
      RETURN
      END

      LOGICAL FUNCTION  WrtDim ( DimNam, MinVal )

c          Write name of too-small symbolic dimension and
c          the value it should be increased to;  return 'TRUE'

c      INPUT :  DimNam = Name of symbolic dimension which is too small
c                        ( CHARACTER, any length )
c               Minval = Value to which that dimension should be
c                        increased (at least)

      CHARACTER*(*)  DimNam
      INTEGER        MinVal


      WRITE ( *, '(3A,I7)' )  ' ****  Symbolic dimension  ', DimNam,
     $                     '  should be increased to at least ', MinVal
      WrtDim = .TRUE.
      RETURN
      END

      LOGICAL FUNCTION  TstBad( VarNam, RelErr )

c       Write name (VarNam) of variable failing self-test and its
c       percent error from the correct value;  return  'FALSE'.

      CHARACTER*(*)  VarNam
      REAL           RelErr


      TstBad = .FALSE.
      WRITE( *, '(/,3A,E11.2,A)' )
     $       ' Output variable ', VarNam,' differed by ', 100.*RelErr,
     $       ' per cent from correct value.  Self-test failed.'
      RETURN
      END
