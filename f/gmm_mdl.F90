! $Id$
! Purpose: Portable implementation of gamma functions
! Routines are external subroutines so that host program may
! link to machine's native erf functions if available
! fxm: Add incomplete gamma function gammq() to gmm_mdl

! Copyright (C) 1994--2018 Charlie Zender
! You may copy, distribute, and/or modify this software under the terms of the GNU General Public License (GPL) Version 2
! See http://www.gnu.org/copyleft/gpl.html for full license text

! Usage:
! use gmm_mdl ! [mdl] Gamma function gamma()

! Source: libspecfun.a specfun/src.dp/gamma.f
! http://scicomp.ewha.ac.kr/netlib/specfun/
! Modifications by Alf Grini and Charlie Zender
! 20030121: AG Turn into Fortran 90 free-format module
! 20030121: AG Use generic precision
! 20030219: CZ Clean up
! 20030219: CZ Inline use of real()
! 20061125: CZ Fix all all g95 'implicit precision' and 'uninitialized' warnings with CEWIs
! 20070314: CZ Add precision to some naked constants
! 20070330: CZ Change from precision module to shr_kind_mod module

module gmm_mdl ! [mdl] Gamma function gamma()
  use shr_kind_mod,only:r8=>shr_kind_r8 ! [mdl] Precision r8, i8, ...
  
  implicit none
  private ! [stt] Symbols are private unless individually qualified as public
  public::gamma ! [fnc] Compute gamma function gamma(x)
  public::gratio! [sub] Compute incomplete gamma function
  
  INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(15, 60)

contains
  
  real(r8) function gamma(x) ! [fnc] Compute gamma function gamma(x)
    ! This routine calculates the GAMMA function for a real argument X.
    ! Computation is based on an algorithm outlined in reference 1.
    ! The program uses rational functions that approximate the GAMMA
    ! function to at least 20 significant decimal digits.  Coefficients
    ! for the approximation over the interval (1,2) are unpublished.
    ! Those for the approximation for X .GE. 12 are from reference 2.
    ! The accuracy achieved depends on the arithmetic system, the
    ! compiler, the intrinsic functions, and proper selection of the
    ! machine-dependent constants.
    
    ! Explanation of machine-dependent constants
    ! beta   - radix for the floating-point representation
    ! maxexp - the smallest positive power of beta that overflows
    ! XBIG   - the largest argument for which GAMMA(X) is representable
    !          in the machine, i.e., the solution to the equation
    !                  GAMMA(XBIG) = beta**maxexp
    ! XINF   - the largest machine representable floating-point number;
    !          approximately beta**maxexp
    ! EPS    - the smallest positive floating-point number such that
    !          1.0+EPS .GT. 1.0
    ! XMININ - the smallest positive floating-point number such that
    !          1/XMININ is machine representable
    
    !     Approximate values for some important machines are:
    !                            beta       maxexp        XBIG
    ! CRAY-1         (S.P.)        2         8191        966.961
    ! Cyber 180/855
    !   under NOS    (S.P.)        2         1070        177.803
    ! IEEE (IBM/XT,
    !   SUN, etc.)   (S.P.)        2          128        35.040
    ! IEEE (IBM/XT,
    !   SUN, etc.)   (D.P.)        2         1024        171.624
    ! IBM 3033       (D.P.)       16           63        57.574
    ! VAX D-Format   (D.P.)        2          127        34.844
    ! VAX G-Format   (D.P.)        2         1023        171.489
    !                            XINF         EPS        XMININ
    ! CRAY-1         (S.P.)   5.45E+2465   7.11E-15    1.84E-2466
    ! Cyber 180/855
    !   under NOS    (S.P.)   1.26E+322    3.55E-15    3.14E-294
    ! IEEE (IBM/XT,
    !   SUN, etc.)   (S.P.)   3.40E+38     1.19E-7     1.18E-38
    ! IEEE (IBM/XT,
    !   SUN, etc.)   (D.P.)   1.79D+308    2.22D-16    2.23D-308
    ! IBM 3033       (D.P.)   7.23D+75     2.22D-16    1.39D-76
    ! VAX D-Format   (D.P.)   1.70D+38     1.39D-17    5.88D-39
    ! VAX G-Format   (D.P.)   8.98D+307    1.11D-16    1.12D-308
    
    ! Error returns
    ! The program returns the value XINF for singularities or
    ! when overflow would occur.  The computation is believed
    ! to be free of underflow and overflow.
    
    ! Intrinsic functions required are:
    ! INT, DBLE, EXP, LOG, REAL, SIN
    
    ! References: "An Overview of Software Development for Special
    ! Functions", W. J. Cody, Lecture Notes in Mathematics,
    ! 506, Numerical Analysis Dundee, 1975, G. A. Watson
    ! (ed.), Springer Verlag, Berlin, 1976.
    ! Computer Approximations, Hart, Et. Al., Wiley and
    ! sons, New York, 1968.
    ! Latest modification: October 12, 1989
    ! Authors: W. J. Cody and L. Stoltz
    ! Applied Mathematics Division
    ! Argonne National Laboratory
    ! Argonne, IL 60439
    integer I,N
    logical PARITY
    real(r8) &
         ! C,CONV,EPS,FACT,HALF,ONE,P,PI,Q,RES,SQRTPI,SUM,TWELVE, &
         C,EPS,FACT,HALF,ONE,P,PI,Q,RES,SQRTPI,SUM,TWELVE, &
         TWO,X,XBIG,XDEN,XINF,XMININ,XNUM,Y,Y1,YSQ,Z,ZERO
    
    dimension C(7),P(8),Q(8)
    ! Mathematical constants
    data ONE,HALF,TWELVE,TWO,ZERO/1.0E0_r8,0.5E0_r8,12.0E0_r8,2.0E0_r8,0.0E0_r8/, &
         SQRTPI/0.9189385332046727417803297E0_r8/, &
         PI/3.1415926535897932384626434E0_r8/
    ! Machine dependent parameters
    data XBIG,XMININ,EPS/35.040E0_r8,1.18E-38_r8,1.19E-7_r8/,  &
         XINF/3.4E38_r8/
    ! Numerator and denominator coefficients for rational minimax
    ! approximation over (1,2).
    data P/-1.71618513886549492533811E+0_r8,2.47656508055759199108314E+1_r8,  &
         -3.79804256470945635097577E+2_r8,6.29331155312818442661052E+2_r8,    &
         8.66966202790413211295064E+2_r8,-3.14512729688483675254357E+4_r8,    &
         -3.61444134186911729807069E+4_r8,6.64561438202405440627855E+4_r8/  
    data Q/-3.08402300119738975254353E+1_r8,3.15350626979604161529144E+2_r8, & 
         -1.01515636749021914166146E+3_r8,-3.10777167157231109440444E+3_r8, &
         2.25381184209801510330112E+4_r8,4.75584627752788110767815E+3_r8, &
         -1.34659959864969306392456E+5_r8,-1.15132259675553483497211E+5_r8/
    ! Coefficients for minimax approximation over (12, INF).
    data C/-1.910444077728E-03_r8,8.4171387781295E-04_r8,                  &
         -5.952379913043012E-04_r8,7.93650793500350248E-04_r8,              &
         -2.777777777777681622553E-03_r8,8.333333333333333331554247E-02_r8, &
         5.7083835261E-03_r8/
    ! csz++ Get rid of statement function by placing directly inline
    ! Statement functions for conversion between integer and float
    ! CONV(I) = REAL(I,kind=r8)
    ! CD    CONV(I) = DBLE(I)
    ! csz-- 
    PARITY = .FALSE.
    FACT = ONE
    N = 0
    Y = X
    IF (Y .LE. ZERO) THEN
       !  Argument is negative
       Y = -X
       Y1 = AINT(Y,kind=r8)
       RES = Y - Y1
       IF (RES .NE. ZERO) THEN
          IF (Y1 .NE. AINT(Y1*HALF,kind=r8)*TWO) PARITY = .TRUE.
          FACT = -PI / SIN(PI*RES)
          Y = Y + ONE
       ELSE
          RES = XINF
          GO TO 900
       END IF
    END IF
    !  Argument is positive
    IF (Y .LT. EPS) THEN
       !  Argument .LT. EPS
       IF (Y .GE. XMININ) THEN
          RES = ONE / Y
       ELSE
          RES = XINF
          GO TO 900
       END IF
    ELSE IF (Y .LT. TWELVE) THEN
       Y1 = Y
       IF (Y .LT. ONE) THEN
          !  0.0 .LT. argument .LT. 1.0
          Z = Y
          Y = Y + ONE
       ELSE
          !  1.0 .LT. argument .LT. 12.0, reduce argument if necessary
          N = INT(Y) - 1
          Y = Y - REAL(N,kind=r8)
          Z = Y - ONE
       END IF
       !  Evaluate approximation for 1.0 .LT. argument .LT. 2.0
       XNUM = ZERO
       XDEN = ONE
       DO 260 I = 1, 8
          XNUM = (XNUM + P(I)) * Z
          XDEN = XDEN * Z + Q(I)
260       CONTINUE
          RES = XNUM / XDEN + ONE
          IF (Y1 .LT. Y) THEN
             !  Adjust result for case  0.0 .LT. argument .LT. 1.0
             RES = RES / Y1
          ELSE IF (Y1 .GT. Y) THEN
             !  Adjust result for case  2.0 .LT. argument .LT. 12.0
             DO 290 I = 1, N
                RES = RES * Y
                Y = Y + ONE
290             CONTINUE
             END IF
          ELSE
             !  Evaluate for argument .GE. 12.0,
             IF (Y .LE. XBIG) THEN
                YSQ = Y * Y
                SUM = C(7)
                DO 350 I = 1, 6
                   SUM = SUM / YSQ + C(I)
350                CONTINUE
                   SUM = SUM/Y - Y + SQRTPI
                   SUM = SUM + (Y-HALF)*LOG(Y)
                   RES = EXP(SUM)
                ELSE
                   RES = XINF
                   GO TO 900
                END IF
             END IF
             !  Final adjustments and return
             IF (PARITY) RES = -RES
             IF (FACT .NE. ONE) RES = FACT / RES
900          GAMMA = RES
             !CD900 DGAMMA = RES
             RETURN
           end function gamma
           
           !sbc++ 12152005 Added incomplete gamma function subroutine
SUBROUTINE gratio(a, x, ans, qans, ind)
!-----------------------------------------------------------------------

!     EVALUATION OF THE INCOMPLETE GAMMA RATIO FUNCTIONS
!                   P(A,X) AND Q(A,X)

!                     ----------

!  IT IS ASSUMED THAT A AND X ARE NONNEGATIVE, WHERE A AND X ARE NOT BOTH 0.

!  ANS AND QANS ARE VARIABLES.  GRATIO ASSIGNS ANS THE VALUE P(A,X)
!  AND QANS THE VALUE Q(A,X). IND MAY BE ANY INTEGER.
!  IF IND = 0 THEN THE USER IS REQUESTING AS MUCH ACCURACY AS POSSIBLE
!  (UP TO 14 SIGNIFICANT DIGITS).  OTHERWISE, IF IND = 1 THEN ACCURACY
!  IND = 1 THEN ACCURACY IS REQUESTED TO WITHIN 1 UNIT OF THE 6-TH SIGNIFICANT
!  DIGIT, AND IF IND .NE. 0,1 THEN ACCURACY IS REQUESTED TO WITHIN 1 UNIT
!  OF THE 3RD SIGNIFICANT DIGIT.

!  ERROR RETURN ...

!     ANS IS ASSIGNED THE VALUE 2 WHEN A OR X IS NEGATIVE,
!  WHEN A*X = 0, OR WHEN P(A,X) AND Q(A,X) ARE INDETERMINANT.
!  P(A,X) AND Q(A,X) ARE COMPUTATIONALLY INDETERMINANT WHEN
!  X IS EXCEEDINGLY CLOSE TO A AND A IS EXTREMELY LARGE.

!--------------------------------------------------------------------
!  WRITTEN BY ALFRED H. MORRIS, JR.
!     NAVAL SURFACE WARFARE CENTER
!     DAHLGREN, VIRGINIA
!  REVISED ... DEC 1991
!-------------------------
REAL(r8), INTENT(IN)   :: a, x
INTEGER, INTENT(IN)    :: ind
REAL(r8), INTENT(OUT)  :: ans, qans

! Local variables
REAL(r8), PARAMETER :: acc0(3) = (/ 5.E-15_r8, 5.E-7_r8, 5.E-4_r8 /),  &
        big(3) = (/ 25.0_r8, 14.0_r8 , 10.0_r8 /),  &
        e0(3) = (/ .25E-3_r8, .25E-1_r8, 0.14_r8 /),  &
        x0(3) = (/ 31.0_r8, 17.0_r8, 9.7_r8 /), ALOG10 = 2.30258509299405_r8,  &
        rt2pin = .398942280401433_r8, rtpi = 1.77245385090552_r8,  &
        d00 = -.333333333333333_r8, d10 = -.185185185185185E-02_r8,   &
        d20 =  .413359788359788E-02_r8, d30 = .649434156378601E-03_r8,   &
        d40 = -.861888290916712E-03_r8, d50 = -.336798553366358E-03_r8,   &
        d60 =  .531307936463992E-03_r8, d70 = .344367606892378E-03_r8,   &
        d80 = -.652623918595309E-03_r8
REAL(r8), PARAMETER :: a0(4) = (/ -.231272501940775E-02_r8, -.335378520024220E-01_r8,  &
                           -.159840143443990_r8, -.333333333333333_r8 /),  &
        a1(4) = (/ -.398783924370770E-05_r8, -.587926036018402E-03_r8,  &
                   -.491687131726920E-02_r8, -.185185185184291E-02_r8 /),  &
        a2(2) = (/ .669564126155663E-03_r8, .413359788442192E-02_r8 /),  &
        a3(2) = (/ .810586158563431E-03_r8, .649434157619770E-03_r8 /),  &
        a4(2) = (/ -.105014537920131E-03_r8, -.861888301199388E-03_r8 /),  &
        a5(2) = (/ -.435211415445014E-03_r8, -.336806989710598E-03_r8 /),  &
        a6(2) = (/ -.182503596367782E-03_r8, .531279816209452E-03_r8 /),  &
        a7(2) = (/ .443219646726422E-03_r8, .344430064306926E-03_r8 /),  &
        a8(2) = (/ .878371203603888E-03_r8, -.686013280418038E-03_r8 /)
REAL(r8), PARAMETER :: b0(6) = (/ .633763414209504E-06_r8, -.939001940478355E-05_r8,  &
        .239521354917408E-02_r8, .376245718289389E-01_r8, .238549219145773_r8,  &
        .729520430331981_r8 /),  &
        b1(4) = (/ .386325038602125E-02_r8,  &
        .506042559238939E-01_r8, .283344278023803_r8, .780110511677243_r8 /),  &
        b2(5) = (/ -.421924263980656E-03_r8, .650837693041777E-02_r8,  &
        .682034997401259E-01_r8, .339173452092224_r8, .810647620703045_r8 /),  &
        b3(5) = (/ -.632276587352120E-03_r8, .905375887385478E-02_r8,  &
        .906610359762969E-01_r8, .406288930253881_r8, .894800593794972_r8 /),  &
        b4(4) = (/ .322609381345173E-01_r8, .178295773562970_r8,  &
        .591353097931237_r8, .103151890792185E+01_r8 /),  &
        b5(3) = (/ .178716720452422_r8, .600380376956324_r8,  &
        .108515217314415E+01_r8 /),  &
        b6(2) = (/ .345608222411837_r8, .770341682526774_r8 /),  &
        b7(2) = (/ .821824741357866_r8, .115029088777769E+01_r8 /)
REAL(r8) :: d0(6) = (/ .833333333333333E-01_r8, -.148148148148148E-01_r8,  &
        .115740740740741E-02_r8, .352733686067019E-03_r8, -.178755144032922E-03_r8,  &
        .391926317852244E-04_r8 /), d1(4) = (/ -.347222222222222E-02_r8,  &
        .264550264550265E-02_r8, -.990226337448560E-03_r8, .205761316872428E-03_r8 /),  &
        d2(2) = (/ -.268132716049383E-02_r8, .771604938271605E-03_r8 /),  &
        d3(2) = (/ .229472093621399E-03_r8, -.469189494395256E-03_r8 /),  &
        d4(1) = (/ .784039221720067E-03_r8 /), d5(1) = (/ -.697281375836586E-04_r8 /),  &
        d6(1) = (/ -.592166437353694E-03_r8 /)
REAL(r8) :: acc, amn, apn, a2n, a2nm1, b2n, b2nm1, c, c0, c1, c2, c3, c4, c5, &
            c6, c7, c8, e, g, h, j, l, r, rta, rtx, s, sum, t, tol, twoa, u,  &
            w, wk(20), y, z
INTEGER  :: i, iop, m, n, nl1
!-------------------------

!     ****** E IS A MACHINE DEPENDENT CONSTANT. E IS THE SMALLEST
!            FLOATING POINT NUMBER FOR WHICH 1.0 + E > 1.0 .

e = EPSILON(1.0_r8)

!-------------------------
IF (a < 0.0 .OR. x < 0.0) GO TO 320
IF (a == 0.0 .AND. x == 0.0) GO TO 320
IF (a*x == 0.0) GO TO 310

iop = ind + 1
IF (iop /= 1 .AND. iop /= 2) iop = 3
acc = MAX(acc0(iop),e)

!            SELECT THE APPROPRIATE ALGORITHM

IF (a < 1.0) THEN
  IF (a == 0.5) GO TO 290
  IF (x < 1.1) GO TO 100
  r = drcomp(a,x)
  IF (r == 0.0) GO TO 280
  GO TO 160
END IF

IF (a < big(iop)) THEN
  IF (a > x .OR. x >= x0(iop)) GO TO 10
  twoa = a + a
  m = INT(twoa)
  IF (twoa /= REAL(m)) GO TO 10
  i = m / 2
  IF (a == REAL(i)) GO TO 130
  GO TO 140
END IF

l = x / a
IF (l == 0.0) GO TO 270
s = 0.5 + (0.5-l)
z = drlog(l)
IF (z >= 700.0/a) GO TO 300
y = a * z
rta = SQRT(a)
IF (ABS(s) <= e0(iop)/rta) GO TO 230
IF (ABS(s) <= 0.4) GO TO 180

10 r = drcomp(a,x)
IF (r == 0.0) GO TO 310
IF (x > MAX(a, ALOG10)) THEN
  IF (x < x0(iop)) GO TO 160
ELSE

!                 TAYLOR SERIES FOR P/R

  apn = a + 1.0
  t = x / apn
  wk(1) = t
  DO n = 2, 20
    apn = apn + 1.0
    t = t * (x/apn)
    IF (t <= 1.E-3_r8) GO TO 30
    wk(n) = t
  END DO
  n = 20

  30 sum = t
  tol = 0.5 * acc
  40 apn = apn + 1.0
  t = t * (x/apn)
  sum = sum + t
  IF (t > tol) GO TO 40

  nl1 = n - 1
  DO m = 1, nl1
    n = n - 1
    sum = sum + wk(n)
  END DO
  ans = (r/a) * (1.0+sum)
  qans = 0.5 + (0.5-ans)
  RETURN
END IF

!                 ASYMPTOTIC EXPANSION

amn = a - 1.0
t = amn / x
wk(1) = t
DO n = 2, 20
  amn = amn - 1.0
  t = t * (amn/x)
  IF (ABS(t) <= 1.E-3_r8) GO TO 70
  wk(n) = t
END DO
n = 20

70 sum = t
80 IF (ABS(t) >= acc) THEN
  amn = amn - 1.0
  t = t * (amn/x)
  sum = sum + t
  GO TO 80
END IF

nl1 = n - 1
DO m = 1, nl1
  n = n - 1
  sum = sum + wk(n)
END DO
qans = (r/x) * (1.0+sum)
ans = 0.5 + (0.5-qans)
RETURN

!             TAYLOR SERIES FOR P(A,X)/X**A

100 l = 3.0
c = x
sum = x / (a+3.0)
tol = 3.0 * acc / (a+1.0)
110 l = l + 1.0
c = -c * (x/l)
t = c / (a+l)
sum = sum + t
IF (ABS(t) > tol) GO TO 110
j = a * x * ((sum/6.0 - 0.5/(a+2.0))*x + 1.0/(a+1.0))

z = a * LOG(x)
h = dgam1(a)
g = 1.0 + h
IF (x >= 0.25) THEN
  IF (a < x/2.59) GO TO 120
ELSE
  IF (z > -.13394) GO TO 120
END IF

w = EXP(z)
ans = w * g * (0.5+(0.5-j))
qans = 0.5 + (0.5-ans)
RETURN

120 l = drexp(z)
w = 0.5 + (0.5+l)
qans = (w*j-l) * g - h
IF (qans < 0.0) GO TO 280
ans = 0.5 + (0.5-qans)
RETURN

!             FINITE SUMS FOR Q WHEN A >= 1
!                 AND 2*A IS AN INTEGER

130 sum = EXP(-x)
t = sum
n = 1
c = 0.0
GO TO 150

140 rtx = SQRT(x)
sum = derfc1(0,rtx)
t = EXP(-x) / (rtpi*rtx)
n = 0
c = -0.5

150 IF (n /= i) THEN
  n = n + 1
  c = c + 1.0
  t = (x*t) / c
  sum = sum + t
  GO TO 150
END IF
qans = sum
ans = 0.5 + (0.5-qans)
RETURN

!              CONTINUED FRACTION EXPANSION

160 tol = MAX(8.0*e,4.0*acc)
a2nm1 = 1.0
a2n = 1.0
b2nm1 = x
b2n = x + (1.0-a)
c = 1.0
170 a2nm1 = x * a2n + c * a2nm1
b2nm1 = x * b2n + c * b2nm1
c = c + 1.0
t = c - a
a2n = a2nm1 + t * a2n
b2n = b2nm1 + t * b2n

a2nm1 = a2nm1 / b2n
b2nm1 = b2nm1 / b2n
a2n = a2n / b2n
b2n = 1.0
IF (ABS(a2n-a2nm1/b2nm1) >= tol*a2n) GO TO 170

qans = r * a2n
ans = 0.5 + (0.5-qans)
RETURN

180 IF (ABS(s) <= 2.0*e .AND. a*e*e > 3.28E-3_r8) GO TO 320
c = EXP(-y)
w = 0.5 * derfc1(1,SQRT(y))
u = 1.0 / a
z = SQRT(z+z)
IF (l < 1.0) z = -z
IF (iop == 2) THEN
  GO TO 200
ELSE IF (iop > 2) THEN
  GO TO 210
END IF

IF (ABS(s) <= 1.E-3_r8) GO TO 240

!            USING THE MINIMAX APPROXIMATIONS

c0 = (((a0(1)*z + a0(2))*z + a0(3))*z + a0(4)) / ((((((b0(1)*z + b0(2))*z +  &
     b0(3))*z + b0(4))*z + b0(5))*z + b0(6))*z + 1.0)
c1 = (((a1(1)*z + a1(2))*z + a1(3))*z + a1(4)) / ((((b1(1)*z + b1(2))*z +   &
     b1(3))*z + b1(4))*z + 1.0)
c2 = (a2(1)*z + a2(2)) / (((((b2(1)*z + b2(2))*z + b2(3))*z + b2(4))*z +  &
     b2(5))*z + 1.0)
c3 = (a3(1)*z + a3(2)) / (((((b3(1)*z + b3(2))*z + b3(3))*z + b3(4))*z +  &
     b3(5))*z + 1.0)
c4 = (a4(1)*z + a4(2)) / ((((b4(1)*z + b4(2))*z + b4(3))*z + b4(4))*z + 1.0)
c5 = (a5(1)*z + a5(2)) / (((b5(1)*z + b5(2))*z + b5(3))*z + 1.0)
c6 = (a6(1)*z + a6(2)) / ((b6(1)*z + b6(2))*z + 1.0)
c7 = (a7(1)*z + a7(2)) / ((b7(1)*z + b7(2))*z + 1.0)
c8 = a8(1) * z + a8(2)
t = (((((((c8*u + c7)*u + c6)*u + c5)*u + c4)*u + c3)*u + c2)*u + c1)*u + c0
GO TO 220

!                    TEMME EXPANSION

200 c0 = (((((d0(6)*z+d0(5))*z+d0(4))*z+d0(3))*z+d0(2))*z+d0(1)) * z +d00
c1 = (((d1(4)*z+d1(3))*z+d1(2))*z+d1(1)) * z + d10
c2 = d2(1) * z + d20
t = (c2*u+c1) * u + c0
GO TO 220

210 t = ((d0(3)*z+d0(2))*z+d0(1)) * z + d00

220 IF (l >= 1.0) THEN
  qans = c * (w+rt2pin*t/rta)
  ans = 0.5 + (0.5-qans)
  RETURN
END IF
ans = c * (w-rt2pin*t/rta)
qans = 0.5 + (0.5-ans)
RETURN

!               TEMME EXPANSION FOR L = 1

230 IF (a*e*e > 3.28E-3_r8) GO TO 320
c = 0.5 + (0.5-y)
w = (0.5 - SQRT(y)*(0.5 + (0.5 - y/3.0))/rtpi) / c
u = 1.0 / a
z = SQRT(z+z)
IF (l < 1.0) z = -z
IF (iop == 2) THEN
  GO TO 250
ELSE IF (iop > 2) THEN
  GO TO 260
END IF

240 c0 = ((d0(3)*z+d0(2))*z+d0(1)) * z + d00
c1 = ((d1(3)*z+d1(2))*z+d1(1)) * z + d10
c2 = (d2(2)*z+d2(1)) * z + d20
c3 = (d3(2)*z+d3(1)) * z + d30
c4 = d4(1) * z + d40
c5 = d5(1) * z + d50
c6 = d6(1) * z + d60
t = (((((((d80*u+d70)*u+c6)*u+c5)*u+c4)*u+c3)*u+c2)*u+c1) * u + c0
GO TO 220

250 c0 = (d0(2)*z+d0(1)) * z + d00
c1 = d1(1) * z + d10
t = (d20*u+c1) * u + c0
GO TO 220

260 t = d0(1) * z + d00
GO TO 220

!                     SPECIAL CASES

270 ans = 0.0
qans = 1.0
RETURN

280 ans = 1.0
qans = 0.0
RETURN

290 IF (x < 0.25) THEN
  ans = derf(SQRT(x))
  qans = 0.5 + (0.5-ans)
  RETURN
END IF
qans = derfc1(0,SQRT(x))
ans = 0.5 + (0.5-qans)
RETURN

300 IF (ABS(s) <= 2.0*e) GO TO 320
310 IF (x <= a) GO TO 270
GO TO 280

!                     ERROR RETURN

320 ans = 2.0
RETURN
END SUBROUTINE gratio



SUBROUTINE gaminv(a, x, x0, p, q, ierr)
!--------------------------------------------------------------------

!          INVERSE INCOMPLETE GAMMA RATIO FUNCTION

!  GIVEN POSITIVE A, AND NONEGATIVE P AND Q WHERE P + Q = 1.
!  THEN X IS COMPUTED WHERE P(A,X) = P AND Q(A,X) = Q. SCHRODER
!  ITERATION IS EMPLOYED. THE ROUTINE ATTEMPTS TO COMPUTE X
!  TO 10 SIGNIFICANT DIGITS IF THIS IS POSSIBLE FOR THE
!  PARTICULAR COMPUTER ARITHMETIC BEING USED.

!                     ------------

!  X IS A VARIABLE. IF P = 0 THEN X IS ASSIGNED THE VALUE 0,
!  AND IF Q = 0 THEN X IS SET TO THE LARGEST FLOATING POINT
!  NUMBER AVAILABLE. OTHERWISE, GAMINV ATTEMPTS TO OBTAIN
!  A SOLUTION FOR P(A,X) = P AND Q(A,X) = Q. IF THE ROUTINE
!  IS SUCCESSFUL THEN THE SOLUTION IS STORED IN X.

!  X0 IS AN OPTIONAL INITIAL APPROXIMATION FOR X. IF THE USER DOES NOT
!  WISH TO SUPPLY AN INITIAL APPROXIMATION, THEN SET X0 <= 0.

!  IERR IS A VARIABLE THAT REPORTS THE STATUS OF THE RESULTS.
!  WHEN THE ROUTINE TERMINATES, IERR HAS ONE OF THE FOLLOWING
!  VALUES ...

!    IERR =  0    THE SOLUTION WAS OBTAINED.  ITERATION WAS NOT USED.
!    IERR >  0    THE SOLUTION WAS OBTAINED.  IERR ITERATIONS WERE PERFORMED.
!    IERR = -2    (INPUT ERROR) A <= 0
!    IERR = -3    NO SOLUTION WAS OBTAINED. THE RATIO Q/A IS TOO LARGE.
!    IERR = -4    (INPUT ERROR) P OR Q IS NEGATIVE, OR P + Q .NE. 1.
!    IERR = -6    20 ITERATIONS WERE PERFORMED. THE MOST
!                 RECENT VALUE OBTAINED FOR X IS GIVEN.
!                 THIS CANNOT OCCUR IF X0 <= 0.
!    IERR = -7    ITERATION FAILED. NO VALUE IS GIVEN FOR X.
!                 THIS MAY OCCUR WHEN X IS APPROXIMATELY 0.
!    IERR = -8    A VALUE FOR X HAS BEEN OBTAINED, BUT THE ROUTINE IS NOT
!                 CERTAIN OF ITS ACCURACY.
!                 ITERATION CANNOT BE PERFORMED IN THIS CASE.
!                 IF X0 <= 0, THIS CAN OCCUR ONLY WHEN P OR Q IS
!                 APPROXIMATELY 0.  IF X0 IS POSITIVE THEN THIS CAN OCCUR
!                 WHEN A IS EXCEEDINGLY CLOSE TO X AND A IS EXTREMELY
!                 LARGE (SAY A >= 1.E20).

!--------------------------------------------------------------------
!  WRITTEN BY ALFRED H. MORRIS, JR.
!     NAVAL SURFACE WARFARE CENTER
!     DAHLGREN, VIRGINIA
!  REVISED ... JANUARY 1992
!------------------------
REAL(r8), INTENT(IN)   :: a, x0, p, q
REAL(r8), INTENT(OUT)  :: x
INTEGER, INTENT(OUT)   :: ierr

! Local variables
REAL(r8), PARAMETER :: ln10 = 2.302585_r8, bmin(2) = (/ 1.E-28_r8, 1.E-13_r8 /),  &
        emin(2) = (/ 2.E-03_r8, 6.E-03_r8 /), c = .577215664901533_r8, tol = 1.E-5_r8
REAL(r8)  :: amax, amin, am1, ap1, ap2, ap3, apn, b, c1, c2, c3, c4, c5, d,  &
             e, eps, e2, g, h, pn, qg, qn, r, rta, s, sum, s2, t, u, w, xmin, &
             xn, y, z
INTEGER   :: ier, iop
!------------------------
!     LN10 = LN(10)
!     C = EULER CONSTANT
!------------------------

!     ****** E AND XMIN ARE MACHINE DEPENDENT CONSTANTS. E IS THE
!            SMALLEST NUMBER FOR WHICH 1.0 + E > 1.0, AND XMIN
!            IS THE SMALLEST POSITIVE NUMBER.

b=0.0_r8 ! csz 20061125 CEWI
e = EPSILON(1.0_r8)
xmin = TINY(1.0_r8)

!------------------------
x = 0.0
IF (a > 0.0) THEN
  IF (p < 0.0 .OR. q < 0.0) GO TO 120
  t = ((p+q)-0.5) - 0.5
  !SCAPPS IF (ABS(t) > 5.0*MAX(e,1.D-15)) GO TO 120
  IF (ABS(t) > 5.0*MAX(e,1.E-15_r8)) GO TO 120

  ierr = 0
  xmin = xmin / e
  IF ((p/e) > xmin) THEN
    IF ((q/e) <= xmin) GO TO 160
    IF (a == 1.0) GO TO 100

    e2 = e + e
    amax = 0.4E-10_r8 / (e*e)
    !SCAPPS eps = MAX(100.0*e,1.D-10)
    eps = MAX(100.0*e,1.E-10_r8)
    iop = 1
    IF (e > 1.E-10_r8) iop = 2
    xn = x0
    IF (x0 <= 0.0) THEN

!        SELECTION OF THE INITIAL APPROXIMATION XN OF X
!                       WHEN A < 1

      IF (a > 1.0) GO TO 30
      g = dgamma(a+1.0)
      qg = q * g
      IF (qg == 0.0) GO TO 160
      b = qg / a
      IF (qg > 0.6*a) GO TO 20
      IF (a < 0.30 .AND. b >= 0.35) THEN
        t = EXP(-(b+c))
        u = t * EXP(t)
        xn = t * EXP(u)
        GO TO 50
      END IF

      IF (b >= 0.45) GO TO 20
      IF (b == 0.0) GO TO 160
      y = -LOG(b)
      s = 0.5 + (0.5-a)
      z = LOG(y)
      t = y - s * z
      IF (b >= 0.15) THEN
        xn = y - s * LOG(t) - LOG(1.0 + s/(t+1.0))
        GO TO 80
      END IF
      IF (b > 1.E-2_r8) THEN
        u = ((t + 2.0*(3.0-a))*t + (2.0-a)*(3.0-a)) / ((t + (5.0-a))*t + 2.0)
        xn = y - s * LOG(t) - LOG(u)
        GO TO 80
      END IF
      10 c1 = -s * z
      c2 = -s * (1.0+c1)
      c3 = s * ((0.5*c1+(2.0-a))*c1+(2.5-1.5*a))
      c4 = -s * (((c1/3.0+(2.5-1.5*a))*c1+((a-6.0)*a+7.0))*c1 +  &
                ((11.0*a-46.0)*a+47.0)/6.0)
      c5 = -s * ((((-c1/4.0+(11.0*a-17.0)/6.0)*c1+((-3.0*a+13.0)*a  &
           -13.0))*c1+0.5*(((2.0*a-25.0)*a+72.0)*a-61.0))*c1 +  &
           (((25.0*a-195.0)*a+477.0)*a-379.0)/12.0)
      xn = ((((c5/y+c4)/y+c3)/y+c2)/y+c1) + y
      IF (a > 1.0) GO TO 80
      IF (b > bmin(iop)) GO TO 80
      x = xn
      RETURN

      20 IF (b*q <= 1.E-8_r8) THEN
        xn = EXP(-(q/a+c))
      ELSE
        IF (p > 0.9) THEN
          xn = EXP((dlnrel(-q) + dgmln1(a))/a)
        ELSE
          xn = EXP(LOG(p*g)/a)
        END IF
      END IF

      IF (xn == 0.0) GO TO 110
      t = 0.5 + (0.5 - xn/(a+1.0))
      xn = xn / t
      GO TO 50

!        SELECTION OF THE INITIAL APPROXIMATION XN OF X
!                       WHEN A > 1

      30 t = p - 0.5
      IF (q < 0.5) t = 0.5 - q
      CALL dpni(p, q, t, s, ier)
      IF (ier /= 0) WRITE(*, *) '** Error in call to PNI from GAMINV **'

      rta = SQRT(a)
      s2 = s * s
      xn = (((12.0*s2 - 243.0)*s2 - 923.0)*s2 + 1472.0) / 204120.0
      xn = (xn/a + s*((9.0*s2 + 256.0)*s2 - 433.0)/(38880.0*rta)) -   &
           ((3.0*s2 + 7.0)*s2 - 16.0) / 810.0
      xn = a + s * rta + (s2-1.0) / 3.0 + s * (s2-7.0) / (36.0*rta) + xn / a
! csz++ add precision to naked constant
      xn = MAX(xn, 0.0_r8)

      amin = 20.0
      IF (e < 1.E-8_r8) amin = 250.0
      IF (a >= amin) THEN
        x = xn
        d = 0.5 + (0.5-x/a)
        IF (ABS(d) <= 1.E-1_r8) RETURN
      END IF

      IF (p > 0.5) THEN
        IF (xn < 3.0*a) GO TO 80
        w = LOG(q)
        y = -(w + dgamln(a))
! csz++ add precision to naked constant
        d = MAX(2.0_r8,a*(a-1.0))
        IF (y >= ln10*d) THEN
          s = 1.0 - a
          z = LOG(y)
          GO TO 10
        END IF
        t = a - 1.0
        xn = y + t * LOG(xn) - dlnrel(-t/(xn+1.0))
        xn = y + t * LOG(xn) - dlnrel(-t/(xn+1.0))
        GO TO 80
      END IF

      ap1 = a + 1.0
      IF (xn > 0.70*ap1) GO TO 60
      w = LOG(p) + dgamln(ap1)
      IF (xn <= 0.15*ap1) THEN
        ap2 = a + 2.0
        ap3 = a + 3.0
        x = EXP((w+x)/a)
        x = EXP((w+x-LOG(1.0+(x/ap1)*(1.0+x/ap2)))/a)
        x = EXP((w+x-LOG(1.0+(x/ap1)*(1.0+x/ap2)))/a)
        x = EXP((w+x-LOG(1.0+(x/ap1)*(1.0+(x/ap2)*(1.0+x/ap3))))/a)
        xn = x
        IF (xn <= 1.E-2_r8*ap1) THEN
          IF (xn <= emin(iop)*ap1) RETURN
          GO TO 60
        END IF
      END IF

      apn = ap1
      t = xn / apn
      sum = 1.0 + t
      40 apn = apn + 1.0
      t = t * (xn/apn)
      sum = sum + t
      IF (t > 1.E-4_r8) GO TO 40
      t = w - LOG(sum)
      xn = EXP((xn+t)/a)
      xn = xn * (1.0-(a*LOG(xn)-xn-t)/(a-xn))
      GO TO 60
    END IF

!                 SCHRODER ITERATION USING P

    50 IF (p > 0.5) GO TO 80
    60 IF (p <= xmin) GO TO 150
    am1 = (a-0.5) - 0.5
    70 IF (a > amax) THEN
      d = 0.5 + (0.5-xn/a)
      IF (ABS(d) <= e2) GO TO 150
    END IF

    IF (ierr >= 20) GO TO 130
    ierr = ierr + 1
    CALL gratio(a,xn,pn,qn,0)
    IF (pn == 0.0.OR.qn == 0.0) GO TO 150
    r = drcomp(a,xn)
    IF (r < xmin) GO TO 150
    t = (pn-p) / r
    w = 0.5 * (am1-xn)
    IF (ABS(t) > 0.1.OR.ABS(w*t) > 0.1) THEN
      x = xn * (1.0-t)
      IF (x <= 0.0) GO TO 140
      d = ABS(t)
    ELSE

      h = t * (1.0+w*t)
      x = xn * (1.0-h)
      IF (x <= 0.0) GO TO 140
      IF (ABS(w) >= 1.0 .AND. ABS(w)*t*t <= eps) RETURN
      d = ABS(h)
    END IF
    xn = x
    IF (d > tol) GO TO 70
    IF (d <= eps) RETURN
    IF (ABS(p-pn) <= tol*p) RETURN
    GO TO 70

!                 SCHRODER ITERATION USING Q

    80 IF (q <= xmin) GO TO 150
    am1 = (a-0.5) - 0.5
    90 IF (a > amax) THEN
      d = 0.5 + (0.5-xn/a)
      IF (ABS(d) <= e2) GO TO 150
    END IF

    IF (ierr >= 20) GO TO 130
    ierr = ierr + 1
    CALL gratio(a,xn,pn,qn,0)
    IF (pn == 0.0 .OR. qn == 0.0) GO TO 150
    r = drcomp(a,xn)
    IF (r < xmin) GO TO 150
    t = (q-qn) / r
    w = 0.5 * (am1-xn)
    IF (ABS(t) > 0.1 .OR. ABS(w*t) > 0.1) THEN
      x = xn * (1.0-t)
      IF (x <= 0.0) GO TO 140
      d = ABS(t)
    ELSE

      h = t * (1.0+w*t)
      x = xn * (1.0-h)
      IF (x <= 0.0) GO TO 140
      IF (ABS(w) >= 1.0 .AND. ABS(w)*t*t <= eps) RETURN
      d = ABS(h)
    END IF
    xn = x
    IF (d > tol) GO TO 90
    IF (d <= eps) RETURN
    IF (ABS(q-qn) <= tol*q) RETURN
    GO TO 90
  END IF

!                       SPECIAL CASES

  ierr = -8
  RETURN

  100 IF (q >= 0.9) THEN
    x = -dlnrel(-p)
    RETURN
  END IF
  x = -LOG(q)
  RETURN
END IF

!                       ERROR RETURN

ierr = -2
RETURN

110 ierr = -3
RETURN

120 ierr = -4
RETURN

130 ierr = -6
RETURN

140 ierr = -7
RETURN

150 x = xn
ierr = -8
RETURN

160 x = HUGE(1.0)
ierr = -8
RETURN
END SUBROUTINE gaminv



FUNCTION derf(x) RESULT(fn_val)
!-----------------------------------------------------------------------
!        REAL (dp) EVALUATION OF THE ERROR FUNCTION
!-----------------------------------------------------------------------
REAL (r8), INTENT(IN)  :: x
REAL (r8)              :: fn_val

! Local variables
REAL (r8)  :: ax, t, w
INTEGER    :: i, k
REAL (r8), PARAMETER :: a(21) = (/ .1283791670955125738961589031215_r8,  &
        -.3761263890318375246320529677070_r8,  &
        .1128379167095512573896158902931_r8,  &
        -.2686617064513125175943235372542E-01_r8,  &
        .5223977625442187842111812447877E-02_r8,  &
        -.8548327023450852832540164081187E-03_r8,  &
        .1205533298178966425020717182498E-03_r8,  &
        -.1492565035840625090430728526820E-04_r8,  &
        .1646211436588924261080723578109E-05_r8,  &
        -.1636584469123468757408968429674E-06_r8,  &
        .1480719281587021715400818627811E-07_r8,  &
        -.1229055530145120140800510155331E-08_r8,  &
        .9422759058437197017313055084212E-10_r8,  &
        -.6711366740969385085896257227159E-11_r8,  &
        .4463222608295664017461758843550E-12_r8,  &
        -.2783497395542995487275065856998E-13_r8,  &
        .1634095572365337143933023780777E-14_r8,  &
        -.9052845786901123985710019387938E-16_r8,  &
        .4708274559689744439341671426731E-17_r8,  &
        -.2187159356685015949749948252160E-18_r8,  &
        .7043407712019701609635599701333E-20_r8 /)
!-------------------------------

!                     ABS(X) <= 1

ax = ABS(x)
IF (ax <= 1._r8) THEN
  t = x * x
  w = a(21)
  DO i = 1, 20
    k = 21 - i
    w = t * w + a(k)
  END DO
  fn_val = x * (1._r8+w)
  RETURN
END IF

!                     ABS(X) > 1

IF (ax < 8.5_r8) THEN
  fn_val = 0.5_r8 + (0.5_r8 - EXP(-x*x)*derfc0(ax))
  IF (x < 0._r8) fn_val = -fn_val
  RETURN
END IF

!                 LIMIT VALUE FOR LARGE X

fn_val = SIGN(1._r8,x)
RETURN
END FUNCTION derf



FUNCTION derfc1(ind, x) RESULT(fn_val)
!--------------------------------------------------------------------

!      EVALUATION OF THE COMPLEMENTARY ERROR FUNCTION

!       DERFC1(IND,X) = ERFC(X)           IF IND = 0
!       DERFC1(IND,X) = EXP(X*X)*ERFC(X)  OTHERWISE

!--------------------------------------------------------------------
INTEGER, INTENT(IN)    :: ind
REAL (r8), INTENT(IN)  :: x
REAL (r8)              :: fn_val

! Local variables
REAL (r8)  :: ax, t, w
INTEGER    :: i, k
REAL (r8), PARAMETER :: a(21) = (/ .1283791670955125738961589031215_r8,  &
        -.3761263890318375246320529677070_r8,  &
        .1128379167095512573896158902931_r8,  &
        -.2686617064513125175943235372542E-01_r8,  &
        .5223977625442187842111812447877E-02_r8,  &
        -.8548327023450852832540164081187E-03_r8,  &
        .1205533298178966425020717182498E-03_r8,  &
        -.1492565035840625090430728526820E-04_r8,  &
        .1646211436588924261080723578109E-05_r8,  &
        -.1636584469123468757408968429674E-06_r8,  &
        .1480719281587021715400818627811E-07_r8,  &
        -.1229055530145120140800510155331E-08_r8,  &
        .9422759058437197017313055084212E-10_r8,  &
        -.6711366740969385085896257227159E-11_r8,  &
        .4463222608295664017461758843550E-12_r8,  &
        -.2783497395542995487275065856998E-13_r8,  &
        .1634095572365337143933023780777E-14_r8,  &
        -.9052845786901123985710019387938E-16_r8,  &
        .4708274559689744439341671426731E-17_r8,  &
        -.2187159356685015949749948252160E-18_r8,  &
        .7043407712019701609635599701333E-20_r8 /)
!-------------------------------

!                     ABS(X) <= 1

ax = ABS(x)
IF (ax <= 1._r8) THEN
  t = x * x
  w = a(21)
  DO i = 1, 20
    k = 21 - i
    w = t * w + a(k)
  END DO
  fn_val = 0.5_r8 + (0.5_r8-x*(1._r8+w))
  IF (ind /= 0) fn_val = EXP(t) * fn_val
  RETURN
END IF

!                       X < -1

IF (x <= 0._r8) THEN
  IF (x < -8.3_r8) GO TO 20
  IF (ind /= 0) THEN
    fn_val = 2._r8 * EXP(x*x) - derfc0(ax)
    RETURN
  END IF
  fn_val = 2._r8 - EXP(-x*x) * derfc0(ax)
  RETURN
END IF

!                       X > 1

IF (ind /= 0) THEN
  fn_val = derfc0(x)
  RETURN
END IF
fn_val = 0._r8
IF (x > 100._r8) RETURN
t = x * x
IF (t > -dxparg(1)) RETURN
fn_val = EXP(-t) * derfc0(x)
RETURN

!             LIMIT VALUE FOR LARGE NEGATIVE X

20 fn_val = 2._r8
IF (ind /= 0) fn_val = 2._r8 * EXP(x*x)
RETURN
END FUNCTION derfc1


FUNCTION derfc0(x) RESULT(fn_val)
!-----------------------------------------------------------------------
REAL (r8), INTENT(IN)  :: x
REAL (r8)              :: fn_val

!        EVALUATION OF EXP(X**2)*ERFC(X) FOR X >= 1

!--------------------------------------------------------------------
!  WRITTEN BY ALFRED H. MORRIS, JR.
!     NAVAL SURFACE WARFARE CENTER
!     DAHLGREN, VIRGINIA
!     APRIL 1992
!-------------------------------
REAL (r8)            :: t, u, v, z
REAL (r8), PARAMETER :: rpinv = .56418958354775628694807945156077259
REAL (r8), PARAMETER :: p0 = .16506148041280876191828601E-03,  &
                        p1 =  .15471455377139313353998665E-03_r8,  &
                        p2 =  .44852548090298868465196794E-04_r8,  &
                        p3 = -.49177280017226285450486205E-05_r8,  &
                        p4 = -.69353602078656412367801676E-05_r8,  &
                        p5 = -.20508667787746282746857743E-05_r8,  &
                        p6 = -.28982842617824971177267380E-06_r8,  &
                        p7 = -.17272433544836633301127174E-07_r8,  &
                        q1 =  .16272656776533322859856317E+01_r8,  &
                        q2 =  .12040996037066026106794322E+01_r8,  &
                        q3 =  .52400246352158386907601472_r8,  &
                        q4 =  .14497345252798672362384241_r8,  &
                        q5 =  .25592517111042546492590736E-01_r8,  &
                        q6 =  .26869088293991371028123158E-02_r8,  &
                        q7 =  .13133767840925681614496481E-03_r8
REAL (r8), PARAMETER :: r0 =  .145589721275038539045668824025_r8,  &
                        r1 = -.273421931495426482902320421863_r8,  &
                        r2 =  .226008066916621506788789064272_r8,  &
                        r3 = -.163571895523923805648814425592_r8,  &
                        r4 =  .102604312032193978662297299832_r8,  &
                        r5 = -.548023266949835519254211506880E-01_r8,  &
                        r6 =  .241432239725390106956523668160E-01_r8,  &
                        r7 = -.822062115403915116036874169600E-02_r8,  &
! SCAPPS: Compiler error when:r8 =  .180296241564687154310619200000E-02_r8
!REPLACED WITH: r8_gmm =  .180296241564687154310619200000E-02_r8
                        r8_gmm =  .180296241564687154310619200000E-02_r8
REAL (r8), PARAMETER :: a0 = -.45894433406309678202825375E-03_r8,   &
                        a1 = -.12281298722544724287816236E-01_r8,  &
                        a2 = -.91144359512342900801764781E-01_r8,  &
                        a3 = -.28412489223839285652511367E-01_r8,  &
                        a4 =  .14083827189977123530129812E+01_r8,  &
                        a5 =  .11532175281537044570477189E+01_r8,  &
                        a6 = -.72170903389442152112483632E+01_r8,  &
                        a7 = -.19685597805218214001309225E+01_r8,  &
                        a8 =  .93846891504541841150916038E+01_r8,  &
                        b1 =  .25136329960926527692263725E+02_r8,  &
                        b2 =  .15349442087145759184067981E+03_r8,  &
                        b3 = -.29971215958498680905476402E+03_r8,  &
                        b4 = -.33876477506888115226730368E+04_r8,  &
                        b5 =  .28301829314924804988873701E+04_r8,  &
                        b6 =  .22979620942196507068034887E+05_r8,  &
                        b7 = -.24280681522998071562462041E+05_r8,  &
                        b8 = -.36680620673264731899504580E+05_r8,  &
                        b9 =  .42278731622295627627042436E+05_r8,  &
                        b10=  .28834257644413614344549790E+03_r8,  &
                        b11=  .70226293775648358646587341E+03_r8
REAL (r8), PARAMETER :: c0 = -.7040906288250128001000086E-04_r8,   &
                        c1 = -.3858822461760510359506941E-02_r8,  &
                        c2 = -.7708202127512212359395078E-01_r8,  &
                        c3 = -.6713655014557429480440263_r8,  &
                        c4 = -.2081992124162995545731882E+01_r8,  &
                        c5 =  .2898831421475282558867888E+01_r8,  &
                        c6 =  .2199509380600429331650192E+02_r8,  &
                        c7 =  .2907064664404115316722996E+01_r8,  &
                        c8 = -.4766208741588182425380950E+02_r8,  &
                        d1 =  .5238852785508439144747174E+02_r8,  &
                        d2 =  .9646843357714742409535148E+03_r8,  &
                        d3 =  .7007152775135939601804416E+04_r8,  &
                        d4 =  .8515386792259821780601162E+04_r8,  &
                        d5 = -.1002360095177164564992134E+06_r8,  &
                        d6 = -.2065250031331232815791912E+06_r8,  &
                        d7 =  .5695324805290370358175984E+06_r8,  &
                        d8 =  .6589752493461331195697873E+06_r8,  &
                        d9 = -.1192930193156561957631462E+07_r8
REAL (r8), PARAMETER :: e0 = .540464821348814822409610122136_r8,  &
                        e1 = -.261515522487415653487049835220E-01_r8, &
                        e2 = -.288573438386338758794591212600E-02_r8, &
                        e3 = -.529353396945788057720258856000E-03_r8
REAL (r8), PARAMETER :: s1 = .75000000000000000000_r8,   &
        s2  = -.18750000000000000000E+01_r8, s3  = .65625000000000000000E+01_r8,  &
        s4  = -.29531250000000000000E+02_r8, s5  = .16242187500000000000E+03_r8,  &
        s6  = -.10557421875000000000E+04_r8, s7  = .79180664062500000000E+04_r8,  &
        s8  = -.67303564453125000000E+05_r8, s9  = .63938386230468750000E+06_r8,  &
        s10 = -.67135305541992187500E+07_r8, s11 = .77205601373291015625E+08_r8
!-------------------------------
!     RPINV = 1/SQRT(PI)
!-------------------------------

!                     1 <= X <= 2

IF (x <= 2._r8) THEN
  u = ((((((p7*x + p6)*x + p5)*x + p4)*x + p3)*x + p2)*x + p1) * x + p0
  v = ((((((q7*x + q6)*x + q5)*x + q4)*x + q3)*x + q2)*x + q1) * x + 1._r8
  t = (x-3.75_r8) / (x+3.75_r8)
  fn_val = (((((((((u/v)*t + r8_gmm)*t + r7)*t + r6)*t + r5)*t + r4)*t + r3)*t + &
           r2)*t + r1) * t + r0
  RETURN
END IF

!                     2 < X <= 4

IF (x <= 4._r8) THEN
  z = 1._r8 / (2.5_r8 + x*x)
  u = (((((((a8*z + a7)*z + a6)*z + a5)*z + a4)*z + a3)*z + a2)*z + a1) * z + a0
  v = ((((((((((b11*z + b10)*z + b9)*z + b8)*z + b7)*z + b6)*z + b5)*z +  &
      b4)*z + b3)*z + b2)*z + b1) * z + 1._r8
  t = 13._r8 * z - 1._r8
  fn_val = ((((u/v)*t + e2)*t + e1)*t + e0) / x
  RETURN
END IF

!                     4 < X < 50

IF (x < 50._r8) THEN
  z = 1._r8 / (2.5_r8 + x*x)
  u = (((((((c8*z + c7)*z + c6)*z + c5)*z + c4)*z + c3)*z + c2)*z + c1) * z + &
      c0
  v = ((((((((d9*z + d8)*z + d7)*z + d6)*z + d5)*z + d4)*z + d3)*z + d2)*z +  &
      d1)*z + 1._r8
  t = 13._r8 * z - 1._r8
  fn_val = (((((u/v)*t + e3)*t + e2)*t + e1)*t + e0) / x
  RETURN
END IF

!                        X >= 50

t = (1._r8/x) ** 2
z = (((((((((((s11*t + s10)*t + s9)*t + s8)*t + s7)*t + s6)*t + s5)*t +  &
    s4)*t + s3)*t + s2)*t + s1)*t - 0.5_r8) * t + 1._r8
fn_val = rpinv * (z/x)
RETURN
END FUNCTION derfc0



FUNCTION drexp(x) RESULT(fn_val)
!-----------------------------------------------------------------------
!            EVALUATION OF THE FUNCTION EXP(X) - 1
!-----------------------------------------------------------------------
REAL (r8), INTENT(IN)  :: x
REAL (r8)              :: fn_val

! Local variables
REAL (r8)  :: e, w, z
REAL (r8)  :: a0 = .248015873015873015873016E-04_r8,   &
    a1 = -.344452080605731005808147E-05_r8, a2 = .206664230430046597475413E-06_r8,  &
    a3 = -.447300111094328162971036E-08_r8, a4 = .114734027080634968083920E-11_r8,  &
    b1 = -.249994190011341852652396_r8, b2 = .249987228833107957725728E-01_r8,  &
    b3 = -.119037506846942249362528E-02_r8, b4 = .228908693387350391768682E-04_r8
REAL (r8) :: c1 = .1666666666666666666666666666666667_r8,   &
             c2 = .4166666666666666666666666666666667E-01_r8,   &
             c3 = .8333333333333333333333333333333333E-02_r8,   &
             c4 = .1388888888888888888888888888888889E-02_r8,   &
             c5 = .1984126984126984126984126984126984E-03_r8
!---------------------------
IF (ABS(x) <= 0.15_r8) THEN

!     Z IS A MINIMAX APPROXIMATION OF THE SERIES

!             C6 + C7*X + C8*X**2 + ....

!     THIS APPROXIMATION IS ACCURATE TO WITHIN
!     1 UNIT OF THE 23-RD SIGNIFICANT DIGIT.
!     THE RESULTING VALUE FOR W IS ACCURATE TO
!     WITHIN 1 UNIT OF THE 33-RD SIGNIFICANT DIGIT.

  z = ((((a4*x + a3)*x + a2)*x + a1)*x + a0) /  &
      ((((b4*x + b3)*x + b2)*x + b1)*x + 1._r8)
  w = ((((((z*x + c5)*x + c4)*x + c3)*x + c2)*x + c1)*x + 0.5_r8)*x + 1._r8
  fn_val = x * w
  RETURN
END IF

IF (x >= 0._r8) THEN
  e = EXP(x)
  fn_val = e * (0.5_r8 + (0.5_r8 - 1._r8/e))
  RETURN
END IF
IF (x >= -77._r8) THEN
  fn_val = (EXP(x) - 0.5_r8) - 0.5_r8
  RETURN
END IF
fn_val = -1._r8
RETURN
END FUNCTION drexp



FUNCTION dlnrel(a) RESULT(fn_val)
!-----------------------------------------------------------------------
!            EVALUATION OF THE FUNCTION LN(1 + A)
!-----------------------------------------------------------------------
REAL (r8), INTENT(IN)  :: a
REAL (r8)              :: fn_val

! Local variables
REAL (r8) :: t, t2, w, z
REAL (r8) :: p0 = .7692307692307692307680E-01_r8,   &
       p1 = -.1505958055914600184836_r8, p2 =  .9302355725278521726,   &
       p3 = -.1787900022182327735804E-01_r8, q1 = -.2824412139355646910683E+01_r8,  &
       q2 =  .2892424216041495392509E+01_r8, q3 = -.1263560605948009364422E+01_r8,  &
       q4 =  .1966769435894561313526_r8
REAL (r8) :: c1 = .3333333333333333333333333333333_r8,   &
             c2 = .2000000000000000000000000000000_r8,   &
             c3 = .1428571428571428571428571428571_r8,   &
             c4 = .1111111111111111111111111111111_r8,   &
             c5 = .9090909090909090909090909090909E-01_r8
!-------------------------
IF (ABS(a) >= 0.375_r8) THEN
  t = 1._r8 + a
  IF (a < 0._r8) t = 0.5_r8 + (0.5_r8 + a)
  fn_val = LOG(t)
  RETURN
END IF

!     W IS A MINIMAX APPROXIMATION OF THE SERIES

!            C6 + C7*T**2 + C8*T**4 + ...

!     THIS APPROXIMATION IS ACCURATE TO WITHIN 1.6 UNITS OF THE 21-ST
!     SIGNIFICANT DIGIT.
!     THE RESULTING VALUE FOR 1._dp + T2*Z IS ACCURATE TO WITHIN 1 UNIT OF
!     THE 30-TH SIGNIFICANT DIGIT.

t = a / (a + 2._r8)
t2 = t * t
w = (((p3*t2 + p2)*t2 + p1)*t2 + p0) /  &
    ((((q4*t2 + q3)*t2 + q2)*t2 + q1)*t2 + 1._r8)
z = ((((w*t2 + c5)*t2 + c4)*t2 + c3)*t2 + c2)*t2 + c1
fn_val = 2._r8 * t * (1._r8 + t2*z)
RETURN
END FUNCTION dlnrel



FUNCTION drlog(x) RESULT(fn_val)
!--------------------------------------------------------------------
!          EVALUATION OF THE FUNCTION X - 1 - LN(X)
!--------------------------------------------------------------------

REAL (r8), INTENT(IN)  :: x
REAL (r8)              :: fn_val

! Local variables
REAL (r8), PARAMETER :: a = .566749439387323789126387112411845E-01_r8,   &
                        b = .456512608815524058941143273395059E-01_r8
REAL (r8), PARAMETER :: p0 = .7692307692307692307680E-01_r8,   &
       p1 = -.1505958055914600184836_r8, p2 = .9302355725278521726,   &
       p3 = -.1787900022182327735804E-01_r8, q1 = -.2824412139355646910683E+01_r8, &
       q2 = .2892424216041495392509E+01_r8, q3 = -.1263560605948009364422E+01_r8,  &
       q4 = .1966769435894561313526_r8
REAL (r8), PARAMETER :: c1 = .333333333333333333333333333333333_r8,   &
                        c2 = .200000000000000000000000000000000_r8,   &
                        c3 = .142857142857142857142857142857143_r8,   &
                        c4 = .111111111111111111111111111111111_r8,   &
                        c5 = .909090909090909090909090909090909E-01_r8
REAL (r8)  :: r, t, u, up2, w, w1, z
!-------------------------
!     A = DRLOG (0.7)
!     B = DRLOG (4/3)
!-------------------------
IF (x >= 0.61_r8 .AND. x <= 1.57_r8) THEN
  IF (x >= 0.82_r8) THEN
    IF (x > 1.18_r8) GO TO 10

!                 ARGUMENT REDUCTION

    u = (x-0.5_r8) - 0.5_r8
    up2 = u + 2._r8
    w1 = 0._r8
    GO TO 20
  END IF

  u = (x-0.7_r8) / 0.7_r8
  up2 = u + 2._r8
  w1 = a - u * 0.3_r8
  GO TO 20

  10 t = 0.75_r8 * (x-1._r8)
  u = t - 0.25_r8
  up2 = t + 1.75_r8
  w1 = b + u / 3._r8

!                  SERIES EXPANSION

  20 r = u / up2
  t = r * r

!     Z IS A MINIMAX APPROXIMATION OF THE SERIES

!            C6 + C7*R**2 + C8*R**4 + ...

!     FOR THE INTERVAL (0.0, 0.375). THE APPROXIMATION IS ACCURATE
!     TO WITHIN 1.6 UNITS OF THE 21-ST SIGNIFICANT DIGIT.

  z = (((p3*t + p2)*t + p1)*t + p0) / ((((q4*t + q3)*t + q2)*t + q1)*t + 1._r8)

  w = ((((z*t + c5)*t + c4)*t + c3)*t + c2) * t + c1
  fn_val = r * (u-2._r8*t*w) + w1
  RETURN
END IF

r = (x-0.5_r8) - 0.5_r8
fn_val = r - LOG(x)
RETURN
END FUNCTION drlog



FUNCTION dpdel(x) RESULT(fn_val)
!--------------------------------------------------------------------

!  COMPUTATION OF THE FUNCTION DEL(X) FOR  X >= 10  WHERE
!  LN(GAMMA(X)) = (X - 0.5)*LN(X) - X + 0.5*LN(2*PI) + DEL(X)

!                      --------

!  THE SERIES FOR DPDEL ON THE INTERVAL 0.0 TO 1.0 DERIVED BY
!  A.H. MORRIS FROM THE CHEBYSHEV SERIES IN THE SLATEC LIBRARY
!  OBTAINED BY WAYNE FULLERTON (LOS ALAMOS).

!--------------------------------------------------------------------
REAL (r8), INTENT(IN)  :: x
REAL (r8)              :: fn_val

! Local variables
REAL (r8), PARAMETER :: a(15) = (/ .833333333333333333333333333333E-01_r8,  &
        -.277777777777777777777777752282E-04_r8,  &
         .793650793650793650791732130419E-07_r8,  &
        -.595238095238095232389839236182E-09_r8,  &
         .841750841750832853294451671990E-11_r8,  &
        -.191752691751854612334149171243E-12_r8,  &
         .641025640510325475730918472625E-14_r8,  &
        -.295506514125338232839867823991E-15_r8,  &
         .179643716359402238723287696452E-16_r8,  &
        -.139228964661627791231203060395E-17_r8,  &
         .133802855014020915603275339093E-18_r8,  &
        -.154246009867966094273710216533E-19_r8,  &
         .197701992980957427278370133333E-20_r8,  &
        -.234065664793997056856992426667E-21_r8,  &
         .171348014966398575409015466667E-22_r8 /)
REAL (r8) :: t, w
INTEGER   :: i, k
!-----------------------------------------------------------------------
t = (10._r8/x) ** 2
w = a(15)
DO i = 1, 14
  k = 15 - i
  w = t * w + a(k)
END DO
fn_val = w / x
RETURN
END FUNCTION dpdel



FUNCTION dsin1(x) RESULT(fn_val)
!--------------------------------------------------------------------

!             REAL (dp) EVALUATION OF SIN(PI*X)

!                          --------------

!  THE EXPANSION FOR SIN(PI*A) (ABS(A) <= PI/4) USING A1,...,A13
!  IS ACCURATE TO WITHIN 2 UNITS OF THE 40-TH SIGNIFICANT DIGIT, AND
!  THE EXPANSION FOR COS(PI*A) (ABS(A) <= PI/4) USING B1,...,B13
!  IS ACCURATE TO WITHIN 4 UNITS OF THE 40-TH SIGNIFICANT DIGIT.

!--------------------------------------------------------------------
REAL (r8), INTENT(IN)  :: x
REAL (r8)              :: fn_val

! Local variables
REAL (r8)            :: a, t, w
REAL (r8), PARAMETER :: pi = 3.141592653589793238462643383279502884197_r8
REAL (r8), PARAMETER :: a1 = -.1028083791780141522795259479153765743002_r8,   &
      a2  = .3170868848763100170457042079710451905600E-02_r8,   &
      a3  = -.4657026956105571623449026167864697920000E-04_r8,  &
      a4  = .3989844942879455643410226655783424000000E-06_r8,   &
      a5  = -.2237397227721999776371894030796800000000E-08_r8,  &
      a6  = .8847045483056962709715066675200000000000E-11_r8,   &
      a7  = -.2598715447506450292885585920000000000000E-13_r8,  &
      a8  = .5893449774331011070033920000000000000000E-16_r8 ,  &
      a9  = -.1062975472045522550784000000000000000000E-18_r8,   &
      a10 = .1561182648301780992000000000000000000000E-21_r8,    &
      a11 = -.1903193516670976000000000000000000000000E-24_r8,   &
      a12 = .1956617650176000000000000000000000000000E-27_r8,    &
      a13 = -.1711276032000000000000000000000000000000E-30_r8
REAL (r8), PARAMETER :: b1 = -.3084251375340424568385778437461297229882_r8, &
      b2  = .1585434424381550085228521039855226435920E-01_r8,   &
      b3  = -.3259918869273900136414318317506279360000E-03_r8,  &
      b4  = .3590860448591510079069203991239232000000E-05_r8,   &
      b5  = -.2461136950494199754009084061808640000000E-07_r8,  &
      b6  = .1150115912797405152263195572224000000000E-09_r8,   &
      b7  = -.3898073171259675439899172864000000000000E-12_r8,  &
      b8  = .1001886461636271969091584000000000000000E-14_r8,   &
      b9  = -.2019653396886572027084800000000000000000E-17_r8,  &
      b10 = .3278483561466560512000000000000000000000E-20_r8,   &
      b11 = -.4377345082051788800000000000000000000000E-23_r8,  &
      b12 = .4891532381388800000000000000000000000000E-26_r8,   &
      b13 = -.4617089843200000000000000000000000000000E-29_r8
INTEGER  :: max, n
!------------------------

!     ****** MAX IS A MACHINE DEPENDENT CONSTANT. MAX IS THE
!            LARGEST POSITIVE INTEGER THAT MAY BE USED.

!                       MAX = IPMPAR(3)
max = HUGE(3)

!------------------------
a = ABS(x)
t = real(MAX) ! csz 20090804 type-cast int to real
IF (a >= t) THEN
  fn_val = 0._r8
  RETURN
END IF

n = int(a)
t = real(n) ! csz 20090804 type-cast int to real
a = a - t
IF (a <= 0.75_r8) THEN
  IF (a < 0.25_r8) GO TO 10

!                    0.25 <= A <= 0.75

  a = 0.25_r8 + (0.25_r8-a)
  t = 16._r8 * a * a
  fn_val = (((((((((((((b13*t + b12)*t + b11)*t + b10)*t + b9)*t + b8)*t  &
           + b7)*t + b6)*t + b5)*t + b4)*t + b3)*t + b2)*t + b1)*t +  &
           0.5_r8) + 0.5_r8
  GO TO 20
END IF

!                 A < 0.25  OR  A > 0.75

a = 0.25_r8 + (0.75_r8-a)
10 t = 16._r8 * a * a
w = (((((((((((((a13*t + a12)*t + a11)*t + a10)*t + a9)*t + a8)*t + a7)*t  &
    + a6)*t + a5)*t + a4)*t + a3)*t + a2)*t + a1)*t + 0.5_r8) + 0.5_r8
fn_val = pi * a * w

!                        TERMINATION

20 IF (x < 0.0) fn_val = -fn_val
IF (MOD(n,2) /= 0) fn_val = -fn_val
RETURN
END FUNCTION dsin1



FUNCTION dgamma(a) RESULT(fn_val)
!--------------------------------------------------------------------

!             EVALUATION OF THE GAMMA FUNCTION FOR
!                  REAL (dp) ARGUMENTS

!                        -----------

!  DGAMMA(A) IS ASSIGNED THE VALUE 0 WHEN THE GAMMA FUNCTION CANNOT
!  BE COMPUTED.

!--------------------------------------------------------------------
!  WRITTEN BY ALFRED H. MORRIS, JR.
!       NAVAL SURFACE WEAPONS CENTER
!       DAHLGREN, VIRGINIA
!--------------------------------------------------------------------
REAL (r8), INTENT(IN)  :: a
REAL (r8)              :: fn_val

! Local variables
REAL (r8), PARAMETER :: d = 0.41893853320467274178032973640562_r8,  &
                        pi = 3.14159265358979323846264338327950_r8
REAL (r8) :: s, t, x, w
INTEGER   :: j, n
!-----------------------------------------------------------------------
!     D = 0.5*(LN(2*PI) - 1)
!-----------------------------------------------------------------------
s=0.0_r8 ! csz 20061125 CEWI
fn_val = 0._r8
x = a
IF (ABS(a) <= 20._r8) THEN
!-----------------------------------------------------------------------
!             EVALUATION OF DGAMMA(A) FOR ABS(A) <= 20
!-----------------------------------------------------------------------
  t = 1._r8
  n = int(x)
  n = n - 1

!     LET T BE THE PRODUCT OF A-J WHEN A >= 2

  IF (n < 0) THEN
    GO TO 40
  ELSE IF (n == 0) THEN
    GO TO 30
  END IF

  DO j = 1, n
    x = x - 1._r8
    t = x * t
  END DO
  30 x = x - 1._r8
  GO TO 60

!     LET T BE THE PRODUCT OF A+J WHEN A < 1

  40 t = a
  IF (a <= 0._r8) THEN
    n = -n - 1
    IF (n /= 0) THEN
      DO j = 1, n
        x = x + 1._r8
        t = x * t
      END DO
    END IF
    x = (x+0.5_r8) + 0.5_r8
    t = x * t
    IF (t == 0._r8) RETURN
  END IF


!     THE FOLLOWING CODE CHECKS IF 1/T CAN OVERFLOW. THIS
!     CODE MAY BE OMITTED IF DESIRED.

  IF (ABS(t) < 1.E-33_r8) THEN
    IF (ABS(t)*HUGE(1.0_r8) <= 1.000000001_r8) RETURN
    fn_val = 1._r8 / t
    RETURN
  END IF

!     COMPUTE DGAMMA(1 + X) FOR 0 <= X < 1

  60 fn_val = 1._r8 / (1._r8 + dgam1(x))

!     TERMINATION

  IF (a >= 1._r8) THEN
    fn_val = fn_val * t
    RETURN
  END IF
  fn_val = fn_val / t
  RETURN
END IF
!-----------------------------------------------------------------------
!           EVALUATION OF DGAMMA(A) FOR ABS(A) > 20
!-----------------------------------------------------------------------
IF (ABS(a) >= 1.E3_r8) RETURN
IF (a <= 0._r8) THEN
  s = dsin1(a) / pi
  IF (s == 0._r8) RETURN
  x = -a
END IF

!     COMPUTE THE MODIFIED ASYMPTOTIC SUM

w = dpdel(x)

!     FINAL ASSEMBLY

w = (d+w) + (x-0.5_r8) * (LOG(x)-1._r8)
IF (w > dxparg(0)) RETURN
fn_val = EXP(w)
IF (a < 0._r8) fn_val = (1._r8/(fn_val*s)) / x

RETURN
END FUNCTION dgamma



FUNCTION dgam1(x) RESULT(fn_val)
!--------------------------------------------------------------------
!  EVALUATION OF 1/GAMMA(1 + X) - 1  FOR -0.5 <= X <= 1.5
!--------------------------------------------------------------------

!  THE FOLLOWING ARE THE FIRST 49 COEFFICIENTS OF THE MACLAURIN
!  EXPANSION FOR 1/GAMMA(1 + X) - 1. THE COEFFICIENTS ARE
!  CORRECT TO 40 DIGITS.  THE COEFFICIENTS WERE OBTAINED BY
!  ALFRED H. MORRIS JR. (NAVAL SURFACE WARFARE CENTER) AND ARE
!  GIVEN HERE FOR REFERENCE.  ONLY THE FIRST 14 COEFFICIENTS ARE
!  USED IN THIS CODE.

!                        -----------

!  DATA A(1)  / .5772156649015328606065120900824024310422_dp/,
! *     A(2)  /-.6558780715202538810770195151453904812798_dp/,
! *     A(3)  /-.4200263503409523552900393487542981871139E-01_r8/,
! *     A(4)  / .1665386113822914895017007951021052357178_dp/,
! *     A(5)  /-.4219773455554433674820830128918739130165E-01_r8/,
! *     A(6)  /-.9621971527876973562114921672348198975363E-02_r8/,
! *     A(7)  / .7218943246663099542395010340446572709905E-02_r8/,
! *     A(8)  /-.1165167591859065112113971084018388666809E-02_r8/,
! *     A(9)  /-.2152416741149509728157299630536478064782E-03_r8/,
! *     A(10) / .1280502823881161861531986263281643233949E-03_r8/
!  DATA A(11) /-.2013485478078823865568939142102181838229E-04_r8/,
! *     A(12) /-.1250493482142670657345359473833092242323E-05_r8/,
! *     A(13) / .1133027231981695882374129620330744943324E-05_r8/,
! *     A(14) /-.2056338416977607103450154130020572836513E-06_r8/,
! *     A(15) / .6116095104481415817862498682855342867276E-08_r8/,
! *     A(16) / .5002007644469222930055665048059991303045E-08_r8/,
! *     A(17) /-.1181274570487020144588126565436505577739E-08_r8/,
! *     A(18) / .1043426711691100510491540332312250191401E-09_r8/,
! *     A(19) / .7782263439905071254049937311360777226068E-11_r8/,
! *     A(20) /-.3696805618642205708187815878085766236571E-11_r8/
!  DATA A(21) / .5100370287454475979015481322863231802727E-12_r8/,
! *     A(22) /-.2058326053566506783222429544855237419746E-13_r8/,
! *     A(23) /-.5348122539423017982370017318727939948990E-14_r8/,
! *     A(24) / .1226778628238260790158893846622422428165E-14_r8/,
! *     A(25) /-.1181259301697458769513764586842297831212E-15_r8/,
! *     A(26) / .1186692254751600332579777242928674071088E-17_r8/,
! *     A(27) / .1412380655318031781555803947566709037086E-17_r8/,
! *     A(28) /-.2298745684435370206592478580633699260285E-18_r8/,
! *     A(29) / .1714406321927337433383963370267257066813E-19_r8/,
! *     A(30) / .1337351730493693114864781395122268022875E-21_r8/
!  DATA A(31) /-.2054233551766672789325025351355733796682E-21_r8/,
! *     A(32) / .2736030048607999844831509904330982014865E-22_r8/,
! *     A(33) /-.1732356445910516639057428451564779799070E-23_r8/,
! *     A(34) /-.2360619024499287287343450735427531007926E-25_r8/,
! *     A(35) / .1864982941717294430718413161878666898946E-25_r8/,
! *     A(36) /-.2218095624207197204399716913626860379732E-26_r8/,
! *     A(37) / .1297781974947993668824414486330594165619E-27_r8/,
! *     A(38) / .1180697474966528406222745415509971518560E-29_r8/,
! *     A(39) /-.1124584349277088090293654674261439512119E-29_r8/,
! *     A(40) / .1277085175140866203990206677751124647749E-30_r8/
!  DATA A(41) /-.7391451169615140823461289330108552823711E-32_r8/,
! *     A(42) / .1134750257554215760954165259469306393009E-34_r8/,
! *     A(43) / .4639134641058722029944804907952228463058E-34_r8/,
! *     A(44) /-.5347336818439198875077418196709893320905E-35_r8/,
! *     A(45) / .3207995923613352622861237279082794391090E-36_r8/,
! *     A(46) /-.4445829736550756882101590352124643637401E-38_r8/,
! *     A(47) /-.1311174518881988712901058494389922190237E-38_r8/,
! *     A(48) / .1647033352543813886818259327906394145400E-39_r8/,
! *     A(49) /-.1056233178503581218600561071538285049997E-40_r8/

!                        -----------

!  C = A(1) - 1 IS ALSO FREQUENTLY NEEDED. C HAS THE VALUE ...

!  DATA C /-.4227843350984671393934879099175975689578_dp/

!--------------------------------------------------------------------
REAL (r8), INTENT(IN) :: x
REAL (r8)             :: fn_val

! Local variables
REAL (r8) :: d, t, w, z
REAL (r8), PARAMETER :: a0 = .611609510448141581788E-08_r8, a1  &
        = .624730830116465516210E-08_r8, b1 = .203610414066806987300_r8, b2  &
        = .266205348428949217746E-01_r8, b3 = .493944979382446875238E-03_r8, b4  &
        = -.851419432440314906588E-05_r8, b5 = -.643045481779353022248E-05_r8, b6  &
        = .992641840672773722196E-06_r8, b7 = -.607761895722825260739E-07_r8, b8  &
        = .195755836614639731882E-09_r8
REAL (r8), PARAMETER :: p0 = .6116095104481415817861E-08_r8, p1  &
        = .6871674113067198736152E-08_r8, p2 = .6820161668496170657, p3  &
        = .4686843322948848031080E-10_r8, p4 = .1572833027710446286995E-11_r8, p5  &
        = -.1249441572276366213222E-12_r8, p6 = .4343529937408594255178E-14_r8, q1  &
        = .3056961078365221025009_r8, q2 = .5464213086042296536016E-01_r8, q3  &
        = .4956830093825887312, q4 = .2692369466186361192876E-03_r8
REAL (r8), PARAMETER :: c = -.422784335098467139393487909917598_r8, c0  &
        = .577215664901532860606512090082402_r8, c1  &
        = -.655878071520253881077019515145390_r8, c2  &
        = -.420026350340952355290039348754298E-01_r8, c3  &
        = .166538611382291489501700795102105_r8, c4  &
        = -.421977345555443367482083012891874E-01_r8, c5  &
        = -.962197152787697356211492167234820E-02_r8, c6  &
        = .721894324666309954239501034044657E-02_r8, c7  &
        = -.116516759185906511211397108401839E-02_r8, c8  &
        = -.215241674114950972815729963053648E-03_r8, c9  &
        = .128050282388116186153198626328164E-03_r8, c10  &
        = -.201348547807882386556893914210218E-04_r8, c11  &
        = -.125049348214267065734535947383309E-05_r8, c12  &
        = .113302723198169588237412962033074E-05_r8, c13  &
        = -.205633841697760710345015413002057E-06_r8
!----------------------------
t = x
d = x - 0.5_r8
IF (d > 0._r8) t = d - 0.5_r8
IF (t < 0.0_r8) THEN
  GO TO 30
ELSE IF (t > 0.0_r8) THEN
  GO TO 20
END IF

fn_val = 0._r8
RETURN
!------------

!             CASE WHEN 0 < T <= 0.5

!           W IS A MINIMAX APPROXIMATION FOR
!           THE SERIES A(15) + A(16)*T + ...

!------------
20 w = ((((((p6*t + p5)*t + p4)*t + p3)*t + p2)*t + p1)*t + p0) /   &
       ((((q4*t+q3)*t + q2)*t + q1)*t + 1._r8)
z = (((((((((((((w*t + c13)*t + c12)*t + c11)*t + c10)*t + c9)*t + c8)*t + c7)*t  &
    + c6)*t + c5)*t + c4)*t + c3)*t + c2)*t + c1) * t + c0

IF (d <= 0._r8) THEN
  fn_val = x * z
  RETURN
END IF
fn_val = (t/x) * ((z-0.5_r8)-0.5_r8)
RETURN
!------------

!             CASE WHEN -0.5 <= T < 0

!           W IS A MINIMAX APPROXIMATION FOR
!           THE SERIES A(15) + A(16)*T + ...

!------------
30 w = (a1*t + a0) / ((((((((b8*t + b7)*t + b6)*t + b5)*t + b4)*t + b3)*t + b2)*t + b1)*t+1._r8)
z = (((((((((((((w*t + c13)*t + c12)*t + c11)*t + c10)*t + c9)*t + c8)*t + c7)*t  &
    + c6)*t + c5)*t + c4)*t + c3)*t + c2)*t + c1) * t + c

IF (d <= 0._r8) THEN
  fn_val = x * ((z+0.5_r8)+0.5_r8)
  RETURN
END IF
fn_val = t * z / x
RETURN
END FUNCTION dgam1



FUNCTION dgamln(a) RESULT(fn_val)
!--------------------------------------------------------------------

!        EVALUATION OF LN(GAMMA(A)) FOR POSITIVE A

!--------------------------------------------------------------------
!  WRITTEN BY ALFRED H. MORRIS
!       NAVAL SURFACE WEAPONS CENTER
!       DAHLGREN, VIRGINIA
!--------------------------------------------------------------------
!  D = 0.5*(LN(2*PI) - 1)
!--------------------------------------------------------------------
REAL (r8), INTENT(IN)  :: a
REAL (r8)              :: fn_val

! Local variables
REAL (r8), PARAMETER  :: d = 0.41893853320467274178032973640562_r8
REAL (r8) :: w, x
INTEGER   :: i, n
!--------------------------
IF (a < 0.5_r8) THEN
  fn_val = dgmln1(a) - LOG(a)
  RETURN
END IF
IF (a <= 2.5_r8) THEN
  x = a - 1._r8
  IF (a < 1._r8) x = (a-0.5_r8) - 0.5_r8
  fn_val = dgmln1(x)
  RETURN
END IF

IF (a < 10._r8) THEN
  n = int(a - 1.5_r8)
  x = a
  w = 1._r8
  DO i = 1, n
    x = x - 1._r8
    w = x * w
  END DO
  fn_val = dgmln1(x-1._r8) + LOG(w)
  RETURN
END IF

w = dpdel(a)
fn_val = (d+w) + (a-0.5_r8) * (LOG(a)-1._r8)
RETURN
END FUNCTION dgamln



FUNCTION dgmln1(x) RESULT(fn_val)
!-----------------------------------------------------------------------
!     EVALUATION OF LN(GAMMA(1 + X)) FOR -0.5 <= X <= 1.5
!-----------------------------------------------------------------------
REAL (r8), INTENT(IN)  :: x
REAL (r8)              :: fn_val

! Local variables
REAL (r8)  :: w
!-----------------------
w = dgam1(x)
fn_val = -dlnrel(w)
RETURN
END FUNCTION dgmln1



FUNCTION drcomp(a, x) RESULT(fn_val)
!-----------------------------------------------------------------------
!              EVALUATION OF EXP(-X)*X**A/GAMMA(A)
!-----------------------------------------------------------------------
REAL (r8), INTENT(IN)  :: a, x
REAL (r8)              :: fn_val

! Local variables
REAL (r8), PARAMETER  :: c = .398942280401432677939946059934_r8
REAL (r8)             :: t, w
!--------------------------
!     C = 1/SQRT(2*PI)
!--------------------------
fn_val = 0._r8
IF (x == 0._r8) RETURN
IF (a <= 20._r8) THEN

  t = a * LOG(x) - x
  IF (t < dxparg(1)) RETURN
  IF (a < 1._r8) THEN
    fn_val = (a*EXP(t)) * (1._r8 + dgam1(a))
    RETURN
  END IF
  fn_val = EXP(t) / dgamma(a)
  RETURN
END IF

t = x / a
IF (t == 0._r8) RETURN
w = -(dpdel(a) + a*drlog(t))
IF (w >= dxparg(1)) fn_val = c * SQRT(a) * EXP(w)
RETURN
END FUNCTION drcomp



SUBROUTINE dpni(p, q, d, w, ierr)
!--------------------------------------------------------------------

!      EVALUATION OF THE INVERSE NORMAL DISTRIBUTION FUNCTION

!                        ------------

!  LET F(T) = 1/(SQRT(2*PI)*EXP(-T*T/2)). THEN THE FUNCTION

!     PROB(X) = INTEGRAL FROM MINUS INFINITY TO X OF F(T)

!  IS THE NORMAL DISTRIBUTION FUNCTION OF ZERO MEAN AND UNIT
!  VARIANCE. IT IS ASSUMED THAT P > 0, Q > 0, P + Q = 1,
!  AND D = P - 0.5. THE VALUE W IS COMPUTED WHERE PROB(W) = P.

!  IERR IS A VARIABLE THAT REPORTS THE STATUS OF THE RESULTS.

!    IERR = 0  NO INPUT ERRORS WERE DETECTED. W WAS COMPUTED.
!    IERR = 1  EITHER P OR Q IS INCORRECT.
!    IERR = 2  D IS INCORRECT.

!--------------------------------------------------------------------
REAL (r8), INTENT(IN)   :: p, q, d
REAL (r8), INTENT(OUT)  :: w
INTEGER, INTENT(OUT)    :: ierr

! Local variables
REAL (r8), PARAMETER  :: rt2 = 1.4142135623730950488016887242097_r8
REAL (r8)             :: eps, t, u, v
!------------------------
!     RT2 = SQRT(2)
!------------------------
t = MIN(p,q)
IF (t > 0._r8) THEN
  eps = EPSILON(1.0_r8)
  w = 0.5_r8 + (0.5_r8 - (p+q))
  IF (ABS(w) <= 2._r8*eps) THEN

    u = ABS(d+d)
    v = t + t
    w = derfi(u,v)
    IF (w < 0._r8) GO TO 10

    ierr = 0
    w = rt2 * w
    IF (d < 0._r8) w = -w
    RETURN
  END IF
END IF

!                         ERROR RETURN

ierr = 1
RETURN
10 ierr = 2
RETURN
END SUBROUTINE dpni



FUNCTION derfi(p, q) RESULT(fn_val)
!-----------------------------------------------------------------------
REAL (r8), INTENT(IN)  :: p, q
REAL (r8)              :: fn_val

!               REAL (dp) COMPUTATION OF
!                 THE INVERSE ERROR FUNCTION

!                      ----------------

!  FOR 0 <= P <= 1,  W = DERFI(P,Q) WHERE ERF(W) = P. IT
!  IS ASSUMED THAT Q = 1 - P. IF P < 0, Q <= 0, OR P + Q
!  IS NOT 1, THEN DERFI(P,Q) IS SET TO A NEGATIVE VALUE.

!--------------------------------------------------------------------
!  REFERENCE. MATHEMATICS OF COMPUTATION,OCT.1976,PP.827-830.
!               J.M.BLAIR,C.A.EDWARDS,J.H.JOHNSON
!--------------------------------------------------------------------
REAL (r8) :: c = .5625_r8, c1 = .87890625_r8, c2  &
        = -.2302585092994045684017991454684364E+03_r8, r  &
        = .8862269254527580136490837416705726E+00_r8, eps, f, lnq, s, t, x
REAL (r8) :: a(7) = (/ .841467547194693616E-01_r8,  &
        .160499904248262200E+01_r8, .809451641478547505E+01_r8,  &
        .164273396973002581E+02_r8, .154297507839223692E+02_r8,  &
        .669584134660994039E+01_r8, .108455979679682472E+01_r8 /), a1(7)  &
        = (/ .552755110179178015E+2_r8, .657347545992519152E+3_r8,  &
        .124276851197202733E+4_r8, .818859792456464820E+3_r8,  &
        .234425632359410093E+3_r8, .299942187305427917E+2_r8,  &
        .140496035731853946E+1_r8 /), a2(7) = (/ .500926197430588206E+1_r8,  &
        .111349802614499199E+3_r8, .353872732756132161E+3_r8,  &
        .356000407341490731E+3_r8, .143264457509959760E+3_r8,  &
        .240823237485307567E+2_r8, .140496035273226366E+1_r8 /), a3(11)  &
        = (/ .237121026548776092E4_r8, .732899958728969905E6_r8,  &
        .182063754893444775E7_r8, .269191299062422172E7_r8, .304817224671614253E7_r8  &
       , .130643103351072345E7_r8, .296799076241952125E6_r8,  &
        .457006532030955554E5_r8, .373449801680687213E4_r8, .118062255483596543E3_r8  &
       , .100000329157954960E1_r8 /), a4(9) = (/ .154269429680540807E12_r8,  &
        .430207405012067454E12_r8, .182623446525965017E12_r8,  &
        .248740194409838713E11_r8, .133506080294978121E10_r8,  &
        .302446226073105850E08_r8, .285909602878724425E06_r8,  &
        .101789226017835707E04_r8, .100000004821118676E01_r8 /),  &
        b(7) = (/ .352281538790042405E-02_r8, .293409069065309557E+00_r8,  &
        .326709873508963100E+01_r8, .123611641257633210E+02_r8,  &
        .207984023857547070E+02_r8, .170791197367677668E+02_r8,  &
        .669253523595376683E+01_r8 /), b1(6) = (/ .179209835890172156E+3_r8,  &
        .991315839349539886E+3_r8, .138271033653003487E+4_r8, .764020340925985926E+3_r8,  &
        .194354053300991923E+3_r8, .228139510050586581E+2_r8 /),  &
        b2(6) = (/ .209004294324106981E+2_r8, .198607335199741185E+3_r8,  &
        .439311287748524270E+3_r8, .355415991280861051E+3_r8, .123303672628828521E+3_r8,  &
        .186060775181898848E+2_r8 /), b3(10) = (/ .851911109952055378E6_r8,  &
        .194746720192729966E7_r8, .373640079258593694E7_r8, .397271370110424145E7_r8,  &
        .339457682064283712E7_r8 , .136888294898155938E7_r8, .303357770911491406E6_r8,  &
        .459721480357533823E5_r8, .373762573565814355E4_r8, .118064334590001264E3_r8 /),  &
        b4(9) = (/ .220533001293836387E12_r8, .347822938010402687E12_r8,  &
        .468373326975152250E12_r8, .185251723580351631E12_r8,  &
        .249464490520921771E11_r8, .133587491840784926E10_r8,  &
        .302480682561295591E08_r8, .285913799407861384E06_r8,  &
        .101789250893050230E04_r8 /)
!-----------------------------------------------------------------------
!     C2 = LN(1.E-100)
!     R  = SQRT(PI)/2
!-----------------------------------------------------------------------
IF (p >= 0._r8 .AND. q > 0._r8) THEN
  eps = EPSILON(1.0_r8)
  t = 0.5_r8 + (0.5_r8-(p+q))
  IF (ABS(t) > 3._r8*eps) GO TO 10

!                      0 <= P <= 0.75

  IF (p <= 0.75_r8) THEN
    x = c - p * p
    s = (((((a(1)*x+a(2))*x+a(3))*x+a(4))*x+a(5))*x+a(6)) * x +a(7)
    t = ((((((b(1)*x+b(2))*x+b(3))*x+b(4))*x+b(5))*x+b(6))*x+b(7)) * x + 1._r8
    fn_val = p * (s/t)
    IF (eps > 1.E-19_r8) RETURN

    x = fn_val
    f = derf(x) - p
    fn_val = x - r * EXP(x*x) * f
    RETURN
  END IF

!                    0.75 < P <= 0.9375

  IF (p <= 0.9375_r8) THEN
    x = c1 - p * p
    IF (x <= 0.1_r8) THEN
      s = ((((((a1(1)*x+a1(2))*x+a1(3))*x+a1(4))*x+a1(5))*x+a1(6))*x+a1(7))
      t = ((((((b1(1)*x+b1(2))*x+b1(3))*x+b1(4))*x+b1(5))*x+b1(6))*x+1._r8)
    ELSE

      s = ((((((a2(1)*x+a2(2))*x+a2(3))*x+a2(4))*x+a2(5))*x+a2(6))*x+a2(7))
      t = ((((((b2(1)*x+b2(2))*x+b2(3))*x+b2(4))*x+b2(5))*x+b2(6))*x+1._r8)
    END IF

    fn_val = p * (s/t)
    IF (eps > 1.E-19_r8) RETURN

    x = fn_val
    t = derfc1(1,x) - EXP(x*x) * q
    fn_val = x + r * t
    RETURN
  END IF

!                  1.E-100 <= Q < 0.0625

  lnq = LOG(q)
  x = 1._r8 / SQRT(-lnq)
  IF (lnq >= c2) THEN
    s = (((((((((a3(1)*x+a3(2))*x+a3(3))*x+a3(4))*x+a3(5))*x+  &
    a3(6))*x+a3(7))*x+a3(8))*x+a3(9))*x+a3(10)) * x + a3(11)
    t = (((((((((b3(1)*x+b3(2))*x+b3(3))*x+b3(4))*x+b3(5))*x+  &
    b3(6))*x+b3(7))*x+b3(8))*x+b3(9))*x+b3(10)) * x + 1._r8
  ELSE

!                 1.E-10000 <= Q < 1.E-100

    s = (((((((a4(1)*x+a4(2))*x+a4(3))*x+a4(4))*x+a4(5))*x+a4(6))*  &
    x+a4(7))*x+a4(8)) * x + a4(9)
    t = ((((((((b4(1)*x+b4(2))*x+b4(3))*x+b4(4))*x+b4(5))*x+  &
    b4(6))*x+b4(7))*x+b4(8))*x+b4(9)) * x + 1._r8
  END IF

  fn_val = s / (x*t)
  IF (eps > 5.E-20_r8) RETURN

  x = fn_val
  t = derfc1(1,x)
  f = (LOG(t)-lnq) - x * x
  fn_val = x + r * t * f
  RETURN
END IF

!                         ERROR RETURN

fn_val = -1._r8
RETURN
10 fn_val = -2._r8
RETURN
END FUNCTION derfi

! Contains the NSWC functions IPMPAR, SPMPAR, DPMPAR, EPSLN, DEPSLN,
! EXPARG & DXPARG
!-----------------------------------------------------------------------
!     WRITTEN using F90 intrinsics by
!        Alan Miller
!        CSIRO Mathematical & Information Sciences
!        CLAYTON, VICTORIA, AUSTRALIA 3169
!     Latest revision - 1 February 1997
!-----------------------------------------------------------------------

FUNCTION ipmpar (i) RESULT(fn_val)
!-----------------------------------------------------------------------

!     IPMPAR PROVIDES THE INTEGER MACHINE CONSTANTS FOR THE COMPUTER
!     THAT IS USED. IT IS ASSUMED THAT THE ARGUMENT I IS AN INTEGER
!     HAVING ONE OF THE VALUES 1-10. IPMPAR(I) HAS THE VALUE ...

!  INTEGERS.

!     ASSUME INTEGERS ARE REPRESENTED IN THE N-DIGIT, BASE-A FORM

!               SIGN ( X(N-1)*A**(N-1) + ... + X(1)*A + X(0) )

!               WHERE 0 .LE. X(I) .LT. A FOR I=0,...,N-1.

!     IPMPAR(1) = A, THE BASE (radix).

!     IPMPAR(2) = N, THE NUMBER OF BASE-A DIGITS (digits).

!     IPMPAR(3) = A**N - 1, THE LARGEST MAGNITUDE (huge).

!  FLOATING-POINT NUMBERS.

!     IT IS ASSUMED THAT THE SINGLE AND DOUBLE PRECISION FLOATING
!     POINT ARITHMETICS HAVE THE SAME BASE, SAY B, AND THAT THE
!     NONZERO NUMBERS ARE REPRESENTED IN THE FORM

!               SIGN (B**E) * (X(1)/B + ... + X(M)/B**M)

!               WHERE X(I) = 0,1,...,B-1 FOR I=1,...,M,
!               X(1) .GE. 1, AND EMIN .LE. E .LE. EMAX.

!     IPMPAR(4) = B, THE BASE.

!  SINGLE-PRECISION

!     IPMPAR(5) = M, THE NUMBER OF BASE-B DIGITS.

!     IPMPAR(6) = EMIN, THE SMALLEST EXPONENT E.

!     IPMPAR(7) = EMAX, THE LARGEST EXPONENT E.

!  DOUBLE-PRECISION

!     IPMPAR(8) = M, THE NUMBER OF BASE-B DIGITS.

!     IPMPAR(9) = EMIN, THE SMALLEST EXPONENT E.

!     IPMPAR(10) = EMAX, THE LARGEST EXPONENT E.

!-----------------------------------------------------------------------

IMPLICIT NONE
INTEGER, INTENT(IN) :: i
INTEGER             :: fn_val

fn_val=0 ! csz 20061125 CEWI
SELECT CASE(i)
  CASE( 1)
    fn_val = RADIX(i)
  CASE( 2)
    fn_val = DIGITS(i)
  CASE( 3)
    fn_val = HUGE(i)
  CASE( 4)
    fn_val = RADIX(1.0)
  CASE( 5)
    fn_val = DIGITS(1.0)
  CASE( 6)
    fn_val = MINEXPONENT(1.0)
  CASE( 7)
    fn_val = MAXEXPONENT(1.0)
  CASE( 8)
    fn_val = DIGITS(1.0E0_r8)
  CASE( 9)
    fn_val = MINEXPONENT(1.0E0_r8)
  CASE(10)
    fn_val = MAXEXPONENT(1.0E0_r8)
  CASE DEFAULT
    RETURN
END SELECT

RETURN
END FUNCTION ipmpar



FUNCTION spmpar (i) RESULT(fn_val)
!-----------------------------------------------------------------------

!     SPMPAR PROVIDES THE SINGLE PRECISION MACHINE CONSTANTS FOR
!     THE COMPUTER BEING USED. IT IS ASSUMED THAT THE ARGUMENT
!     I IS AN INTEGER HAVING ONE OF THE VALUES 1, 2, OR 3. IF THE
!     SINGLE PRECISION ARITHMETIC BEING USED HAS M BASE B DIGITS AND
!     ITS SMALLEST AND LARGEST EXPONENTS ARE EMIN AND EMAX, THEN

!        SPMPAR(1) = B**(1 - M), THE MACHINE PRECISION,

!        SPMPAR(2) = B**(EMIN - 1), THE SMALLEST MAGNITUDE,

!        SPMPAR(3) = B**EMAX*(1 - B**(-M)), THE LARGEST MAGNITUDE.
!-----------------------------------------------------------------------

IMPLICIT NONE
INTEGER, INTENT(IN) :: i
REAL                :: fn_val

! Local variable
REAL                :: one = 1.0

fn_val=0.0 ! csz 20061125 CEWI
SELECT CASE (i)
  CASE (1)
    fn_val = EPSILON(one)
  CASE (2)
    fn_val = TINY(one)
  CASE (3)
    fn_val = HUGE(one)
END SELECT

RETURN
END FUNCTION spmpar



FUNCTION dpmpar (i) RESULT(fn_val)
!-----------------------------------------------------------------------

!     DPMPAR PROVIDES THE DOUBLE PRECISION MACHINE CONSTANTS FOR
!     THE COMPUTER BEING USED. IT IS ASSUMED THAT THE ARGUMENT
!     I IS AN INTEGER HAVING ONE OF THE VALUES 1, 2, OR 3. IF THE
!     DOUBLE PRECISION ARITHMETIC BEING USED HAS M BASE B DIGITS AND
!     ITS SMALLEST AND LARGEST EXPONENTS ARE EMIN AND EMAX, THEN

!        DPMPAR(1) = B**(1 - M), THE MACHINE PRECISION,

!        DPMPAR(2) = B**(EMIN - 1), THE SMALLEST MAGNITUDE,

!        DPMPAR(3) = B**EMAX*(1 - B**(-M)), THE LARGEST MAGNITUDE.
!-----------------------------------------------------------------------

IMPLICIT NONE
INTEGER, INTENT(IN) :: i
REAL (r8)           :: fn_val

! Local variable
REAL (r8)    :: one = 1._r8

fn_val=0.0_r8 ! csz 20061125 CEWI
SELECT CASE (i)
  CASE (1)
    fn_val = EPSILON(one)
  CASE (2)
    fn_val = TINY(one)
  CASE (3)
    fn_val = HUGE(one)
END SELECT

RETURN
END FUNCTION dpmpar


FUNCTION epsln () RESULT(fn_val)
!--------------------------------------------------------------------
!     THE EVALUATION OF LN(EPS) WHERE EPS IS THE SMALLEST NUMBER
!     SUCH THAT 1.0 + EPS .GT. 1.0 .  L IS A DUMMY ARGUMENT.
!--------------------------------------------------------------------
IMPLICIT NONE
REAL                :: fn_val

! Local variable
REAL                :: one = 1.0

fn_val = LOG( EPSILON(one) )
RETURN
END FUNCTION epsln


FUNCTION exparg (l) RESULT(fn_val)
!--------------------------------------------------------------------
!     IF L = 0 THEN  EXPARG(L) = THE LARGEST POSITIVE W FOR WHICH
!     EXP(W) CAN BE COMPUTED.
!
!     IF L IS NONZERO THEN  EXPARG(L) = THE LARGEST NEGATIVE W FOR
!     WHICH THE COMPUTED VALUE OF EXP(W) IS NONZERO.
!
!     NOTE... ONLY AN APPROXIMATE VALUE FOR EXPARG(L) IS NEEDED.
!--------------------------------------------------------------------
IMPLICIT NONE
INTEGER, INTENT(IN) :: l
REAL                :: fn_val

! Local variable
REAL                :: one = 1.0

IF (l == 0) THEN
  fn_val = LOG( HUGE(one) )
ELSE
  fn_val = LOG( TINY(one) )
END IF
RETURN
END FUNCTION exparg


FUNCTION depsln () RESULT(fn_val)
!--------------------------------------------------------------------
!     THE EVALUATION OF LN(EPS) WHERE EPS IS THE SMALLEST NUMBER
!     SUCH THAT 1.D0 + EPS .GT. 1.D0 .  L IS A DUMMY ARGUMENT.
!--------------------------------------------------------------------
IMPLICIT NONE
REAL (r8)           :: fn_val

! Local variable
REAL (r8)    :: one = 1._r8

fn_val = LOG( EPSILON(one) )
RETURN
END FUNCTION depsln


FUNCTION dxparg (l) RESULT(fn_val)
!--------------------------------------------------------------------
!     IF L = 0 THEN  DXPARG(L) = THE LARGEST POSITIVE W FOR WHICH
!     DEXP(W) CAN BE COMPUTED.
!
!     IF L IS NONZERO THEN  DXPARG(L) = THE LARGEST NEGATIVE W FOR
!     WHICH THE COMPUTED VALUE OF DEXP(W) IS NONZERO.
!
!     NOTE... ONLY AN APPROXIMATE VALUE FOR DXPARG(L) IS NEEDED.
!--------------------------------------------------------------------
IMPLICIT NONE
INTEGER, INTENT(IN) :: l
REAL (r8)           :: fn_val

! Local variable
REAL (r8)    :: one = 1._r8

IF (l == 0) THEN
  fn_val = LOG( HUGE(one) )
ELSE
  fn_val = LOG( TINY(one) )
END IF
RETURN
END FUNCTION dxparg
! sbc--
end module gmm_mdl ! [mdl] Gamma function gamma()
