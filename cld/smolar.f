codes below are three different MPDATA routines. first goes 1-d, then 2-d on
cartesian domain and finally 2-d on the sphere. the spherical routine is very
close to the cartesian one, except that it explicitly incorporates boundary
conditions characteristic for the problems on the sphere. cartesian routine
can be used in arbitrary coordinates, as it admitts extra field H that may be
considered for the jacobian of coordinate transformation. in truly cartesian
coordinates H should be preset to ``1''. these routines were pulled from the
codes that deal with nontrivial dynamic applications (i.e., they work). some
comments are included in MPDATA2 and they apply to all 3 routines. NOTE: both
codes MPDATA1 and MPDATA2 require some boundary routines (search for ``call b''
call to see example. these may be very simple in the spirit of
cyclic, e.g., x(1,j)=x(n-1,j), or zero-gradient, e.g., x(1,j)=x(2,j),
conditions. MPDATAS calls xbc which is provided and requires nothing else.

confirm this e-mail,please.


      SUBROUTINE MPDATA1(U,X,M,IORD,ISOR,IFCT)
C THIS SUBROUTINE SOLVES 1-D ADVECTION IN CARTESIAN GEOMETRY ONLY
C SEE MPDATA2 FOR SOME COMMENTS
      PARAMETER(NM=161+16,N=NM+1)
      DIMENSION X(M),U(M+1),F(N),V(N),CP(NM),CN(NM)
      REAL MX(NM),MN(NM)
      DATA EP /1.E-15/  
      DATA IDIV/1/
      DONOR(Y1,Y2,A)=AMAX1(0.,A)*Y1+AMIN1(0.,A)*Y2
      VDYF(X1,X2,A)=(ABS(A)-A**2)*(ABS(X2)-ABS(X1))
     1                           /(ABS(X2)+ABS(X1)+EP)
      VCOR3(X0,X1,X2,X3,A)= -A*(1.-3*ABS(A)+2.*A**2)
     1                        *(ABS(X0)+ABS(X3)-ABS(X1)-ABS(X2))
     2                        /(ABS(X0)+ABS(X3)+ABS(X1)+ABS(X2)+EP)/3.
      VCU(A1,A2,A3)=0.25*A2*(A3-A1)
      PP(Y)=AMAX1(0.,Y)
      PN(Y)=AMIN1(0.,Y)
C
C
      DO 1 I=1,N
    1 V(I)=U(I) 
C
      IF(IFCT.EQ.1) THEN
      DO 400 I=2,NM-1
      MX(I)=AMAX1(X(I-1),X(I),X(I+1))
  400 MN(I)=AMIN1(X(I-1),X(I),X(I+1))
      MX(1)=AMAX1(X(1),X(2))
      MX(NM)=AMAX1(X(NM-1),X(NM))
      MN(1)=AMIN1(X(1),X(2))
      MN(NM)=AMIN1(X(NM-1),X(NM))
      ENDIF
C    
C 
                         DO 3 K=1,IORD
C
      DO 331 I=2,NM
  331 F(I)=DONOR(X(I-1),X(I),V(I))
      DO 333 I=2,NM-1
  333 X(I)=X(I)-(F(I+1)-F(I))
      CALL BOUNDARY(X,NM,0)   
C
      IF(K.EQ.IORD) GO TO 6
      DO 501 I=2,NM
  501 F(I)=V(I)
      DO 51 I=2,NM
   51 V(I)=VDYF(X(I-1),X(I),F(I))
C
      IF(IDIV.EQ.1) THEN
      DO 50 I=3,NM-1
      V(I)=V(I)-VCU(F(I-1),F(I),F(I+1))
   50 CONTINUE
      ENDIF
C
      IF(ISOR.EQ.3) THEN
      DO 503 I=3,NM-1
  503 V(I)=V(I)+VCOR3(X(I-2),X(I-1),X(I),X(I+1),F(I))
      ENDIF 
      IF(IFCT.EQ.0) GO TO 3
C                    FLUX LIMITER
C
      DO 401 I=2,NM-1
      MX(I)=AMAX1(X(I-1),X(I),X(I+1),MX(I))
  401 MN(I)=AMIN1(X(I-1),X(I),X(I+1),MN(I))
      MX(1)=AMAX1(X(1),X(2),MX(1))
      MX(NM)=AMAX1(X(NM-1),X(NM),MX(NM))
      MN(1)=AMIN1(X(1),X(2),MN(1))
      MN(NM)=AMIN1(X(NM-1),X(NM),MX(NM))
      DO 399 I=2,NM
  399 F(I)=DONOR(X(I-1),X(I),V(I))
C
      DO 402 I=2,NM-1
      CN(I)=AMIN1(1., (X(I)-MN(I))/(PP(F(I+1))-PN(F(I))+EP))
  402 CP(I)=AMIN1(1., (MX(I)-X(I))/(-PN(F(I+1))+PP(F(I))+EP))
      DO 404 I=3,NM-1
      V(I)=PP(V(I))*
     1 ( AMIN1(CP(I),CN(I-1))*PP(SIGN(1., X(I-1)))
     1  +AMIN1(CP(I-1),CN(I))*PP(SIGN(1.,-X(I-1))) )
     2    +PN(V(I))*
     2 ( AMIN1(CP(I-1),CN(I))*PP(SIGN(1., X(I  )))
     2  +AMIN1(CP(I),CN(I-1))*PP(SIGN(1.,-X(I  ))) )
  404 CONTINUE
      V(2)=PP(V(2))*(CP(2)*PP(SIGN(1.,X(1)))+CN(2)*PP(SIGN(1.,-X(1))))
     2    +PN(V(2))*(CN(2)*PP(SIGN(1.,X(2)))+CP(2)*PP(SIGN(1.,-X(2))))
      V(NM)=PP(V(NM))*
     1 (CN(NM-1)*PP(SIGN(1.,X(NM-1)))+CP(NM-1)*PP(SIGN(1.,-X(NM-1))))
     2     +PN(V(NM))*
     2 (CP(NM-1)*PP(SIGN(1.,X(NM  )))+CN(NM-1)*PP(SIGN(1.,-X(NM  ))))
    3                      CONTINUE
    6 CONTINUE
      RETURN
      END



      SUBROUTINE MPDATA2(U1,U2,X,H,N,M,IORD,ISOR,NONOS,IDIV)
C THIS SUBROUTINE SOLVES 2-D ADVECTIVE TRANSPORT IN CARTESIAN GEOMETRY
C HOWEVER ANY COORDINATE TRANSFORMATION MAY BE INCORPORATED THROUGH MATRIX H
C THAT MAY BE EMPLOYED FOR THE JACOBIAN. FOR TRULY CARTESIAN GEOMETRY H(I,J)=1.
C SHOULD BE PRESET IN THE CALLING PROGRAM
C
C
C U1, U2 ARE ADVECTIVE COURANT NUMBERS (VELOCITY_i*DT/DX_i) DEFINED
C AT STAGGERED POSITIONS SUCH THAT THEIR INDEX PRECEDS THE INDEX OF THE
C TRANSPORTED FIELD X (AS WELL AS THAT OF H WHICH COINCIDES WITH X)
C
C
C               U2(I,J+1)
C
C      U1(I,J)   X(I,J)   U1(I+1,J)
C
C               U2(I,J)
C
C N,M ARE FIRST AND SECOND DIMENSION OF X AND H; NOTE DIMENSION OF VELOCITIES
C U1(N+1,M), U2(N,M+1) THAT INCLUDE EXTRA POINT IN THE STAGGERED DIRECTION
C
C IORD IS A NUMBER OF ANTIDIFFUSIVE ITERATIONS: IORD=1 GIVES CLASSIC UPWIND
C SCHEME, IORD=2 GIVES BASIC SECOND-ORDER MPDATA. IN PRACTICE, THERE IS NO
C GAIN IN USING IORD LARGER THAN 3. IORD=3 AND ISOR=3 GIVE SPECIAL OPTION
C OF THE THIRD ORDER ACCURATE MPDATA. NONOS=0 IS A SIGN PRESERVING MPDATA
C (SAY, REGULAR) NONOS=1 GIVES STRICTLY MONOTONE SCHEME (MORE EXPENSIVE).
C IDIV=0 IS FOR NONDIVERGENT FLOWS, I.E., DU1/DX1+DU2/DX2=0. IF THIS ROUTINE
C IS SPLITTED WITH MPDATA1 TO GENERATE SOLUTION TO 3-D PROBLEM IDIV=1 SHOULD
C BE USED EVEN IF THE FULL 3-D FLOW IS NONDIVERGENT (THIS IS RELATIVELY CHEAP
C OPTION).  IN MPDATAS THERE IS ONE MORE PARAMETER IBC, WHICH INFORMS ROUTINE
C WHETHER THE FIELD X CHANGES SIGN (IBC=-1) OR NOT (IBC=1) WHILE PASSING OVER
C THE POLES
C
C N1 AND N2 ARE ADJUSTABLE PARAMETERS AND THEY ARE EQUAL TO N+1 AND M+1,
C RESPECTIVELY.

      PARAMETER(N1=102,N2=102)
      PARAMETER(N1M=N1-1,N2M=N2-1)
      DIMENSION U1(N+1,M),U2(N,M+1),X(N,M),H(N,M)
C FOLLOWING VARIABLES ARE LOCAL TO THE ROUTINE
      DIMENSION V1(N1,N2M),V2(N1M,N2),F1(N1,N2M),F2(N1M,N2)
     *         ,CP(N1M,N2M),CN(N1M,N2M)
      REAL MX(N1M,N2M),MN(N1M,N2M)
      DATA EP/1.E-10/
C
      DONOR(Y1,Y2,A)=CVMGM(Y2,Y1,A)*A
      VDYF(X1,X2,A,R)=(ABS(A)-A**2/R)*(ABS(X2)-ABS(X1))
     1                               /(ABS(X2)+ABS(X1)+EP)
      VCORR(A,B,Y1,Y2,R)=-0.125*A*B*Y1/(Y2*R)
      VCOR31(A,X0,X1,X2,X3,R)= -(A -3.*ABS(A)*A/R+2.*A**3/R**2)/3.
     1                         *(ABS(X0)+ABS(X3)-ABS(X1)-ABS(X2))
     2                         /(ABS(X0)+ABS(X3)+ABS(X1)+ABS(X2)+EP)
      VCOR32(A,B,Y1,Y2,R)=0.25*B/R*(ABS(A)-2.*A**2/R)*Y1/Y2
      VDIV1(A1,A2,A3,R)=0.25*A2*(A3-A1)/R
      VDIV2(A,B1,B2,B3,B4,R)=0.25*A*(B1+B2-B3-B4)/R
      PP(Y)= AMAX1(0.,Y)
      PN(Y)=-AMIN1(0.,Y)
C
      IF(ISOR.EQ.3) IORD=MAX0(IORD,3)
C
      DO 1 J=1,N2-1
      DO 1 I=1,N1
    1 V1(I,J)=U1(I,J)
      DO 2 J=1,N2
      DO 2 I=1,N1-1
    2 V2(I,J)=U2(I,J)
C
C     
      IF(NONOS.EQ.1) THEN
      DO 400 J=2,N2-2
      DO 400 I=2,N1-2
      MX(I,J)=AMAX1(X(I-1,J),X(I,J),X(I+1,J),X(I,J-1),X(I,J+1))
  400 MN(I,J)=AMIN1(X(I-1,J),X(I,J),X(I+1,J),X(I,J-1),X(I,J+1))
      ENDIF
C
                         DO 3 K=1,IORD
C
      DO 331 J=2,N2-2
      DO 331 I=2,N1-1
  331 F1(I,J)=DONOR(X(I-1,J),X(I,J),V1(I,J))
      DO 332 J=2,N2-1
      DO 332 I=2,N1-2
  332 F2(I,J)=DONOR(X(I,J-1),X(I,J),V2(I,J))
C 
      DO 333 J=2,N2-2
      DO 333 I=2,N1-2
  333 X(I,J)=X(I,J)-(F1(I+1,J)-F1(I,J)+F2(I,J+1)-F2(I,J))/H(I,J)
C
      CALL BOUNDARY(X,N1-1,N2-1,3,0,0)   
C
      IF(K.EQ.IORD) GO TO 6
      DO 49 J=1,N2-1
      DO 49 I=1,N1
   49 F1(I,J)=V1(I,J)
      DO 50 J=1,N2
      DO 50 I=1,N1-1
   50 F2(I,J)=V2(I,J)
      DO 51 J=2,N2-2
      DO 51 I=2,N1-1
   51 V1(I,J)=VDYF(X(I-1,J),X(I,J),V1(I,J),.5*(H(I-1,J)+H(I,J)))
     *       +VCORR(V1(I,J), F2(I-1,J)+F2(I-1,J+1)+F2(I,J+1)+F2(I,J),
     *   ABS(X(I-1,J+1))+ABS(X(I,J+1))-ABS(X(I-1,J-1))-ABS(X(I,J-1)),
     *   ABS(X(I-1,J+1))+ABS(X(I,J+1))+ABS(X(I-1,J-1))+ABS(X(I,J-1))+EP,
     *                 .5*(H(I-1,J)+H(I,J)))
      IF(IDIV.EQ.1) THEN
      DO 511 J=2,N2-2
      DO 511 I=2,N1-1
  511 V1(I,J)=V1(I,J)
     *    -VDIV1(F1(I-1,J),F1(I,J),F1(I+1,J),.5*(H(I-1,J)+H(I,J)))
     *    -VDIV2(F1(I,J),F2(I-1,J+1),F2(I,J+1),F2(I-1,J),F2(I,J),
     *                 .5*(H(I-1,J)+H(I,J)))
      ENDIF
      DO 52 J=2,N2-1
      DO 52 I=2,N1-2
   52 V2(I,J)=VDYF(X(I,J-1),X(I,J),V2(I,J),.5*(H(I,J-1)+H(I,J)))
     *       +VCORR(V2(I,J), F1(I,J-1)+F1(I,J)+F1(I+1,J)+F1(I+1,J-1),
     *   ABS(X(I+1,J-1))+ABS(X(I+1,J))-ABS(X(I-1,J-1))-ABS(X(I-1,J)),
     *   ABS(X(I+1,J-1))+ABS(X(I+1,J))+ABS(X(I-1,J-1))+ABS(X(I-1,J))+EP,
     *                 .5*(H(I,J-1)+H(I,J)))
      IF(IDIV.EQ.1) THEN
      DO 521 J=2,N2-1
      DO 521 I=2,N1-2
  521 V2(I,J)=V2(I,J)
     *    -VDIV1(F2(I,J-1),F2(I,J),F2(I,J+1),.5*(H(I,J-1)+H(I,J)))
     *    -VDIV2(F2(I,J),F1(I+1,J),F1(I+1,J-1),F1(I,J-1),F1(I,J),
     *                 .5*(H(I,J-1)+H(I,J)))
      ENDIF
      IF(ISOR.EQ.3) THEN
      DO 61 J=2,N2-2
      DO 61 I=3,N1-2
   61 V1(I,J)=V1(I,J)     +VCOR31(F1(I,J),
     1        X(I-2,J),X(I-1,J),X(I,J),X(I+1,J),.5*(H(I-1,J)+H(I,J)))
      DO 62 J=2,N2-2
      DO 62 I=3,N1-2
   62 V1(I,J)=V1(I,J)
     1 +VCOR32(F1(I,J),F2(I-1,J)+F2(I-1,J+1)+F2(I,J+1)+F2(I,J),
     *   ABS(X(I,J+1))-ABS(X(I,J-1))-ABS(X(I-1,J+1))+ABS(X(I-1,J-1)),
     *   ABS(X(I,J+1))+ABS(X(I,J-1))+ABS(X(I-1,J+1))+ABS(X(I-1,J-1))+EP,
     *                   .5*(H(I-1,J)+H(I,J)))
      DO 63 J=3,N2-2
      DO 63 I=2,N1-2
   63 V2(I,J)=V2(I,J)     +VCOR31(F2(I,J),
     1        X(I,J-2),X(I,J-1),X(I,J),X(I,J+1),.5*(H(I,J-1)+H(I,J)))
      DO 64 J=3,N2-2
      DO 64 I=2,N1-2
   64 V2(I,J)=V2(I,J)
     1 +VCOR32(F2(I,J),F1(I,J-1)+F1(I+1,J-1)+F1(I+1,J)+F1(I,J),
     *   ABS(X(I+1,J))-ABS(X(I-1,J))-ABS(X(I+1,J-1))+ABS(X(I-1,J-1)),
     *   ABS(X(I+1,J))+ABS(X(I-1,J))+ABS(X(I+1,J-1))+ABS(X(I-1,J-1))+EP,
     *                   .5*(H(I,J-1)+H(I,J)))
      ENDIF
C
C
      IF(NONOS.EQ.0) GO TO 3
C                 NON-OSSCILATORY OPTION
      DO 401 J=2,N2-2
      DO 401 I=2,N1-2
      MX(I,J)=AMAX1(X(I-1,J),X(I,J),X(I+1,J),X(I,J-1),X(I,J+1),MX(I,J))
  401 MN(I,J)=AMIN1(X(I-1,J),X(I,J),X(I+1,J),X(I,J-1),X(I,J+1),MN(I,J))
C
      DO 402 J=2,N2-2 
      DO 402 I=2,N1-1
  402 F1(I,J)=DONOR(X(I-1,J),X(I,J),V1(I,J))
      DO 403 J=2,N2-1
      DO 403 I=2,N1-2
  403 F2(I,J)=DONOR(X(I,J-1),X(I,J),V2(I,J))
      DO 404 J=2,N2-2   
      DO 404 I=2,N1-2
      CP(I,J)=(MX(I,J)-X(I,J))*H(I,J)/
     1(PN(F1(I+1,J))+PP(F1(I,J))+PN(F2(I,J+1))+PP(F2(I,J))+EP)
      CN(I,J)=(X(I,J)-MN(I,J))*H(I,J)/
     1(PP(F1(I+1,J))+PN(F1(I,J))+PP(F2(I,J+1))+PN(F2(I,J))+EP)
  404 CONTINUE
      DO 405 J=3,N2-2 
      DO 405 I=3,N1-2 
      V1(I,J)=PP(V1(I,J))*
     1  ( AMIN1(1.,CP(I,J),CN(I-1,J))*PP(SIGN(1., X(I-1,J)))
     1   +AMIN1(1.,CP(I-1,J),CN(I,J))*PP(SIGN(1.,-X(I-1,J))) )
     2       -PN(V1(I,J))*
     2  ( AMIN1(1.,CP(I-1,J),CN(I,J))*PP(SIGN(1., X(I ,J )))
     2   +AMIN1(1.,CP(I,J),CN(I-1,J))*PP(SIGN(1.,-X(I ,J ))) )
  405 V2(I,J)=PP(V2(I,J))*
     1  ( AMIN1(1.,CP(I,J),CN(I,J-1))*PP(SIGN(1., X(I,J-1)))
     1   +AMIN1(1.,CP(I,J-1),CN(I,J))*PP(SIGN(1.,-X(I,J-1))) )
     1       -PN(V2(I,J))*
     2  ( AMIN1(1.,CP(I,J-1),CN(I,J))*PP(SIGN(1., X(I ,J )))
     2   +AMIN1(1.,CP(I,J),CN(I,J-1))*PP(SIGN(1.,-X(I ,J ))) )
C
    3                      CONTINUE
    6 CONTINUE
      RETURN
      END   

      SUBROUTINE MPDATAS(U1,U2,X,H,N,M,IORD,ISOR,NONOS,IDIV,IBC)
      PARAMETER(N1=130,N2=64)
      DIMENSION U1(N,M),U2(N,M+1),X(N,M),H(N,M)
      DIMENSION V1(N1,N2),V2(N1,N2+1),F1(N1,N2),F2(N1,N2+1)
     *         ,CP(N1,N2),CN(N1,N2),IP(N1)
      REAL MX(N1,N2),MN(N1,N2)
      DATA EP/1.E-15/
C           
      DONOR(Y1,Y2,A)=AMAX1(0.,A)*Y1+AMIN1(0.,A)*Y2
      VDYF(X1,X2,A,R)=(ABS(A)-A**2/R)*(ABS(X2)-ABS(X1))
     1                               /(ABS(X2)+ABS(X1)+EP)
      VCORR(A,B,Y1,Y2,R)=-0.125*A*B*Y1/(Y2*R)
      VCOR31(A,X0,X1,X2,X3,R)= -(A -3.*ABS(A)*A/R+2.*A**3/R**2)/3.
     1                         *(ABS(X0)+ABS(X3)-ABS(X1)-ABS(X2))
     2                         /(ABS(X0)+ABS(X3)+ABS(X1)+ABS(X2)+EP)
      VCOR32(A,B,Y1,Y2,R)=0.25*B/R*(ABS(A)-2.*A**2/R)*Y1/Y2
      VDIV1(A1,A2,A3,R)=0.25*A2*(A3-A1)/R
      VDIV2(A,B1,B2,B3,B4,R)=0.25*A*(B1+B2-B3-B4)/R
      PP(Y)= AMAX1(0.,Y)
      PN(Y)=-AMIN1(0.,Y)              
C
      IF(ISOR.EQ.3) IORD=MAX0(IORD,3)
      DO 287 I=1,N1
  287 IP(I)=MOD(I+(N1-2)/2-1,N1-2)+1
C
      DO 1 J=1,N2
      DO 1 I=1,N1
    1 V1(I,J)=U1(I,J)
      DO 2 J=1,N2+1
      DO 2 I=1,N1
    2 V2(I,J)=U2(I,J)

C                 
      IF(NONOS.EQ.1) THEN           
      DO 400 J=2,N2-1
      DO 400 I=2,N1-1
      MX(I,J)=AMAX1(X(I-1,J),X(I,J),X(I+1,J),X(I,J-1),X(I,J+1))
  400 MN(I,J)=AMIN1(X(I-1,J),X(I,J),X(I+1,J),X(I,J-1),X(I,J+1))
      DO 398 I=2,N1-1
      MX(I,1)=AMAX1(X(I-1,1),X(I,1),X(I+1,1),IBC*X(IP(I),1),X(I,2))
      MN(I,1)=AMIN1(X(I-1,1),X(I,1),X(I+1,1),IBC*X(IP(I),1),X(I,2))
      MX(I,N2)=AMAX1(X(I-1,N2),X(I,N2),X(I+1,N2),X(I,N2-1),
     1                                       IBC*X(IP(I),N2))
  398 MN(I,N2)=AMIN1(X(I-1,N2),X(I,N2),X(I+1,N2),X(I,N2-1),
     1                                  IBC*X(IP(I),N2))
      CALL XBC(MX,N1,N2)
      CALL XBC(MN,N1,N2)
      ENDIF
C    
                        DO 3 K=1,IORD
COMPUTE DONOR-CELL FLUXES
      DO 330 J=1,N2
      DO 330 I=2,N1-1
  330 F1(I,J)=DONOR(X(I-1,J),X(I,J),V1(I,J))
      DO 331 J=2,N2
      DO 331 I=2,N1-1
  331 F2(I,J)=DONOR(X(I,J-1),X(I,J),V2(I,J))
      DO 332 I=2,N1-1
      F2(I,N2+1)=DONOR(X(I,N2),IBC*X(IP(I),N2),V2(I,N2+1))
  332 F2(I,1)=DONOR(IBC*X(IP(I),1),X(I,1),V2(I,1))
      CALL XBC(F1,N1,N2)
      CALL XBC(F2,N1,N2+1)
COMPUTE NEW UPWIND-SOLUTION
      DO 333 J=1,N2
      DO 333 I=2,N1-1
  333 X(I,J)=X(I,J)-(F1(I+1,J)-F1(I,J)+F2(I,J+1)-F2(I,J))/H(I,J)
      CALL XBC(X,N1,N2)
C
      IF(K.EQ.IORD) GO TO 6
CONVERT VELOCITIES TO LOCAL STORAGE
      DO 49 J=1,N2
      DO 49 I=1,N1
   49 F1(I,J)=V1(I,J)
      DO 50 J=1,N2+1
      DO 50 I=1,N1
   50 F2(I,J)=V2(I,J) 

CALCULATE PSEUDO VELOCITIES      
COMPUTE FIRST DIRECTION    
      DO 51 J=2,N2-1
      DO 51 I=2,N1-1
   51 V1(I,J)=VDYF(X(I-1,J),X(I,J),F1(I,J),.5*(H(I-1,J)+H(I,J)))
     * +VCORR(F1(I,J), F2(I-1,J)+F2(I-1,J+1)+F2(I,J+1)+F2(I,J),
     *   ABS(X(I-1,J+1))+ABS(X(I,J+1))-ABS(X(I-1,J-1))-ABS(X(I,J-1)),
     *   ABS(X(I-1,J+1))+ABS(X(I,J+1))+ABS(X(I-1,J-1))+ABS(X(I,J-1))+EP,
     *                 .5*(H(I-1,J)+H(I,J)))
COMPUTE B.C IN Y DIRECTION
      DO 510 I=2,N1-1
      V1(I,1)=VDYF(X(I-1,1),X(I,1),F1(I,1),.5*(H(I-1,1)+H(I,1)))
     * +VCORR(F1(I,1), F2(I-1,1)+F2(I-1,2)+F2(I,2)+F2(I,1),
     *   ABS(X(I-1,2))+ABS(X(I,2))-ABS(X(IP(I-1),1))-ABS(X(IP(I),1)),
     *   ABS(X(I-1,2))+ABS(X(I,2))+ABS(X(IP(I-1),1))+ABS(X(IP(I),1))+EP,
     *                 .5*(H(I-1,1)+H(I,1)))
  510 V1(I,N2)=VDYF(X(I-1,N2),X(I,N2),F1(I,N2),.5*(H(I-1,N2)+H(I,N2)))
     * +VCORR(F1(I,N2), F2(I-1,N2)+F2(I-1,N2+1)+F2(I,N2+1)+F2(I,N2),
     *           ABS(X(IP(I-1),N2))+ABS(X(IP(I),N2))
     *          -ABS(X(I-1,N2-1))-ABS(X(I,N2-1)),
     *           ABS(X(IP(I-1),N2))+ABS(X(IP(I),N2))
     *          +ABS(X(I-1,N2-1))+ABS(X(I,N2-1))+EP,
     *                 .5*(H(I-1,N2)+H(I,N2)))        
      IF(IDIV.EQ.1) THEN
COMPUTE FLOW-DIVERGENCE CORRECTION
      DO 5101 J=1,N2
      DO 5101 I=2,N1-1
 5101 V1(I,J)=V1(I,J)
     *    -VDIV1(F1(I-1,J),F1(I,J),F1(I+1,J),.5*(H(I-1,J)+H(I,J)))
     *    -VDIV2(F1(I,J),F2(I-1,J+1),F2(I,J+1),F2(I-1,J),F2(I,J),
     *                 .5*(H(I-1,J)+H(I,J)))
      ENDIF   

COMPUTE SECOND DIRECTION
      DO 52 J=2,N2
      DO 52 I=2,N1-1
   52 V2(I,J)=VDYF(X(I,J-1),X(I,J),F2(I,J),.5*(H(I,J-1)+H(I,J)))
     * +VCORR(F2(I,J), F1(I,J-1)+F1(I,J)+F1(I+1,J)+F1(I+1,J-1),
     *   ABS(X(I+1,J-1))+ABS(X(I+1,J))-ABS(X(I-1,J-1))-ABS(X(I-1,J)),
     *   ABS(X(I+1,J-1))+ABS(X(I+1,J))+ABS(X(I-1,J-1))+ABS(X(I-1,J))+EP,
     *                 .5*(H(I,J-1)+H(I,J)))
COMPUTE B.C IN Y-DIRECTION
      DO 520 I=2,N1-1
      V2(I,1)=VDYF(X(IP(I),1),X(I,1),F2(I,1),.5*(H(IP(I),1)+H(I,1)))
     * +VCORR(F2(I,1), F1(IP(I),1)+F1(I,1)+F1(I+1,1)+F1(IP(I+1),1),
     *    ABS(X(IP(I+1),1))+ABS(X(I+1,1))
     *   -ABS(X(IP(I-1),1))-ABS(X(I-1,1)),
     *    ABS(X(IP(I+1),1))+ABS(X(I+1,1))
     *   +ABS(X(IP(I-1),1))+ABS(X(I-1,1))+EP,
     *                 .5*(H(IP(I),1)+H(I,1)))
  520 V2(I,N2+1)=VDYF(X(I,N2),X(IP(I),N2),F2(I,N2+1),
     1                        .5*(H(I,N2)+H(IP(I),N2)))
     *+VCORR(F2(I,N2+1),F1(I,N2)+F1(IP(I),N2)+F1(IP(I+1),N2)+F1(I+1,N2),
     *   ABS(X(I+1,N2))+ABS(X(IP(I+1),N2))
     *  -ABS(X(I-1,N2))-ABS(X(IP(I-1),N2)),
     *   ABS(X(I+1,N2))+ABS(X(IP(I+1),N2))
     *  +ABS(X(I-1,N2))+ABS(X(IP(I-1),N2))+EP,
     *                 .5*(H(I,J-1)+H(I,J)))
      IF(IDIV.EQ.1) THEN
      DO 5201 J=2,N2
      DO 5201 I=2,N1-1
 5201 V2(I,J)=V2(I,J)
     *    -VDIV1(F2(I,J-1),F2(I,J),F2(I,J+1),.5*(H(I,J-1)+H(I,J)))
     *    -VDIV2(F2(I,J),F1(I+1,J-1),F1(I+1,J),F1(I,J-1),F1(I,J),
     *                 .5*(H(I,J-1)+H(I,J)))
      DO 5202 I=2,N1-1
      V2(I,1)=V2(I,1)
     *    -VDIV1(-F2(IP(I),1),F2(I,1),F2(I,2),.5*(H(IP(I),1)+H(I,1)))
     *    -VDIV2( F2(I,1),F1(IP(I+1),1),F1(I+1,1),F1(IP(I),1),F1(I,1),
     *                 .5*(H(IP(I),1)+H(I,1)))
 5202 V2(I,N2+1)=V2(I,N2+1)
     *    -VDIV1(F2(I,N2),F2(I,N2+1),-F2(IP(I),N2+1)
     *                              ,.5*(H(I,N2)+H(IP(I),N2)))
     *    -VDIV2(F2(I,N2+1),
     *           F1(I+1,N2),F1(IP(I+1),N2),F1(I,N2),F1(IP(I),N2),
     *                 .5*(H(I,N2)+H(IP(I),N2)))
      ENDIF
C
C THIRD ORDER CORRECTION
      IF(ISOR.EQ.3) THEN
C  FIRST DIRECTION
      DO 61 J=1,N2
      DO 610 I=3,N1-1
  610 V1(I,J)=V1(I,J)     +VCOR31(F1(I,J),
     1        X(I-2,J),X(I-1,J),X(I,J),X(I+1,J),.5*(H(I-1,J)+H(I,J)))
   61 V1(2,J)=V1(2,J)     +VCOR31(F1(2,J),
     1        X(N1-2,J),X(1,J),X(2,J),X(3,J),.5*(H(1,J)+H(2,J)))
      DO 62 J=2,N2-1
      DO 62 I=2,N1-1               
   62 V1(I,J)=V1(I,J)
     1 +VCOR32(F1(I,J),F2(I-1,J)+F2(I-1,J+1)+F2(I,J+1)+F2(I,J),
     *   ABS(X(I,J+1))-ABS(X(I,J-1))-ABS(X(I-1,J+1))+ABS(X(I-1,J-1)),
     *   ABS(X(I,J+1))+ABS(X(I,J-1))+ABS(X(I-1,J+1))+ABS(X(I-1,J-1))+EP,
     *                   .5*(H(I-1,J)+H(I,J)))
C B.C. FOLLOW
      DO 620 I=2,N1-1
      V1(I,1)=V1(I,1)
     1 +VCOR32(F1(I,1),F2(I-1,1)+F2(I-1,2)+F2(I,2)+F2(I,1),
     *   ABS(X(I,2))-ABS(X(IP(I),1))-ABS(X(I-1,2))+ABS(X(IP(I-1),1)),
     *   ABS(X(I,2))+ABS(X(IP(I),1))+ABS(X(I-1,2))+ABS(X(IP(I-1),1))+EP,
     *                   .5*(H(I-1,1)+H(I,1)))
  620 V1(I,N2)=V1(I,N2)
     1 +VCOR32(F1(I,N2),F2(I-1,N2)+F2(I-1,N2+1)+F2(I,N2+1)+F2(I,N2),
     *     ABS(X(IP(I),N2))-ABS(X(I,N2-1))
     *                     -ABS(X(IP(I-1),N2))+ABS(X(I-1,N2-1)),
     *     ABS(X(IP(I),N2))+ABS(X(I,N2-1))
     *                     +ABS(X(IP(I-1),N2))+ABS(X(I-1,N2-1))+EP,
     *                   .5*(H(I-1,N2)+H(I,N2)))
      DO 621 J=1,N2
      V1(1,J)=V1(N1-1,J)
  621 V1(N1,J)=V1(2,J)
C
C  SECOND DIRECTION
      DO 63 J=3,N2-1
      DO 63 I=2,N1-1
   63 V2(I,J)=V2(I,J)     +VCOR31(F2(I,J),
     1        X(I,J-2),X(I,J-1),X(I,J),X(I,J+1),.5*(H(I,J-1)+H(I,J)))  
      DO 630 I=2,N1-1
      V2(I,1)=V2(I,1)     +VCOR31(F2(I,1),
     1    X(IP(I),2),X(IP(I),1),X(I,1),X(I,2),.5*(H(IP(I),1)+H(I,1)))  
      V2(I,2)=V2(I,2)     +VCOR31(F2(I,2),
     1        X(IP(I),1),X(I,1),X(I,2),X(I,3),.5*(H(I,1)+H(I,2)))  
      V2(I,N2)=V2(I,N2)     +VCOR31(F2(I,N2),
     1 X(I,N2-2),X(I,N2-1),X(I,N2),X(IP(I),N2),.5*(H(I,N2-1)+H(I,N2)))  
      V2(I,N2+1)=V2(I,N2+1) +VCOR31(F2(I,N2+1), X(I,N2-1),X(I,N2),
     1    X(IP(I),N2),X(IP(I),N2-1),.5*(H(I,N2)+H(IP(I),N2)))  
  630 CONTINUE
      DO 64 J=2,N2
      DO 64 I=2,N1-1
   64 V2(I,J)=V2(I,J)
     1 +VCOR32(F2(I,J),F1(I,J-1)+F1(I+1,J-1)+F1(I+1,J)+F1(I,J),
     *  ABS(X(I+1,J))-ABS(X(I-1,J))-ABS(X(I+1,J-1))+ABS(X(I-1,J-1)),
     *  ABS(X(I+1,J))+ABS(X(I-1,J))+ABS(X(I+1,J-1))+ABS(X(I-1,J-1))+EP,
     *                  .5*(H(I,J-1)+H(I,J)))
      ENDIF

CALL B.C IN X DIRECTION
      CALL XBC(V1,N1,N2)
      CALL XBC(V2,N1,N2+1)
C
      IF(NONOS.EQ.0) GO TO 3
C                 NON-OSSCILATORY OPTION
      DO 401 J=2,N2-1
      DO 401 I=2,N1-1
      MX(I,J)=AMAX1(X(I-1,J),X(I,J),X(I+1,J),X(I,J-1),X(I,J+1),MX(I,J))
  401 MN(I,J)=AMIN1(X(I-1,J),X(I,J),X(I+1,J),X(I,J-1),X(I,J+1),MN(I,J))
      DO 3981 I=2,N1-1
      MX(I,1)=AMAX1(X(I-1,1),X(I,1),X(I+1,1),IBC*X(IP(I),1),
     1                                        X(I,2),MX(I,1))
      MN(I,1)=AMIN1(X(I-1,1),X(I,1),X(I+1,1),IBC*X(IP(I),1),
     1                                        X(I,2),MN(I,1))
      MX(I,N2)=AMAX1(X(I-1,N2),X(I,N2),X(I+1,N2),X(I,N2-1),
     1                              IBC*X(IP(I),N2),MX(I,N2))
 3981 MN(I,N2)=AMIN1(X(I-1,N2),X(I,N2),X(I+1,N2),X(I,N2-1),
     1                              IBC*X(IP(I),N2),MN(I,N2))
      CALL XBC(MX,N1,N2)
      CALL XBC(MN,N1,N2)
C
      DO 4021 J=1,N2
      DO 4021 I=2,N1-1
 4021 F1(I,J)=DONOR(X(I-1,J),X(I,J),V1(I,J))
      DO 403 J=2,N2
      DO 403 I=2,N1-1
  403 F2(I,J)=DONOR(X(I,J-1),X(I,J),V2(I,J))
      DO 410 I=2,N1-1    
      F2(I,N2+1)=DONOR(X(I,N2),IBC*X(IP(I),N2),V2(I,N2+1))
  410 F2(I,1)=DONOR(IBC*X(IP(I),1),X(I,1),V2(I,1))
      CALL XBC(F1,N1,N2)
      CALL XBC(F2,N1,N2+1)

      DO 404 J=1,N2
      DO 404 I=2,N1-1
      CP(I,J)=(MX(I,J)-X(I,J))*H(I,J)/
     1(PN(F1(I+1,J))+PP(F1(I,J))+PN(F2(I,J+1))+PP(F2(I,J))+EP)
      CN(I,J)=(X(I,J)-MN(I,J))*H(I,J)/
     1(PP(F1(I+1,J))+PN(F1(I,J))+PP(F2(I,J+1))+PN(F2(I,J))+EP)
  404 CONTINUE
      CALL XBC(CP,N1,N2)
      CALL XBC(CN,N1,N2)

      DO 405 J=2,N2
      DO 405 I=2,N1-1
      V1(I,J)=PP(V1(I,J))*
     1  ( AMIN1(1.,CP(I,J),CN(I-1,J))*PP(SIGN(1., X(I-1,J)))
     1   +AMIN1(1.,CP(I-1,J),CN(I,J))*PP(SIGN(1.,-X(I-1,J))) )
     2       -PN(V1(I,J))*
     2  ( AMIN1(1.,CP(I-1,J),CN(I,J))*PP(SIGN(1., X(I ,J )))
     2   +AMIN1(1.,CP(I,J),CN(I-1,J))*PP(SIGN(1.,-X(I ,J ))) )
  405 V2(I,J)=PP(V2(I,J))*
     1  ( AMIN1(1.,CP(I,J),CN(I,J-1))*PP(SIGN(1., X(I,J-1)))
     1   +AMIN1(1.,CP(I,J-1),CN(I,J))*PP(SIGN(1.,-X(I,J-1))) )
     1       -PN(V2(I,J))*
     2  ( AMIN1(1.,CP(I,J-1),CN(I,J))*PP(SIGN(1., X(I ,J )))
     2   +AMIN1(1.,CP(I,J),CN(I,J-1))*PP(SIGN(1.,-X(I ,J ))) )
C B.C. FOLLOW
      DO 4051 I=2,N1-1
      V2(I,1)=PP(V2(I,1))*
     1  ( AMIN1(1.,CP(I,1),CN(IP(I),1))*PP(SIGN(1., IBC*X(IP(I),1)))
     1   +AMIN1(1.,CP(IP(I),1),CN(I,1))*PP(SIGN(1.,-IBC*X(IP(I),1))) )
     1       -PN(V2(I,1))*
     2  ( AMIN1(1.,CP(IP(I),1),CN(I,1))*PP(SIGN(1., X(I ,1 )))
     2   +AMIN1(1.,CP(I,1),CN(IP(I),1))*PP(SIGN(1.,-X(I ,1 ))) )
      V2(I,N2+1)=PP(V2(I,N2+1))*
     1  ( AMIN1(1.,CP(IP(I),N2),CN(I,N2))*PP(SIGN(1., X(I,N2)))
     1   +AMIN1(1.,CP(I,N2),CN(IP(I),N2))*PP(SIGN(1.,-X(I,N2))) )
     1       -PN(V2(I,N2+1))*
     2( AMIN1(1.,CP(I,N2),CN(IP(I),N2))*PP(SIGN(1., IBC*X(IP(I),N2)))
     2 +AMIN1(1.,CP(IP(I),N2),CN(I,N2))*PP(SIGN(1.,-IBC*X(IP(I),N2))) )
 4051 V1(I,1)=PP(V1(I,1))*
     1  ( AMIN1(1.,CP(I,1),CN(I-1,1))*PP(SIGN(1., X(I-1,1)))
     1   +AMIN1(1.,CP(I-1,1),CN(I,1))*PP(SIGN(1.,-X(I-1,1))) )
     2       -PN(V1(I,1))*
     2  ( AMIN1(1.,CP(I-1,1),CN(I,1))*PP(SIGN(1., X(I ,1 )))
     2   +AMIN1(1.,CP(I,1),CN(I-1,1))*PP(SIGN(1.,-X(I ,1 ))) )
      CALL XBC(V1,N1,N2)
      CALL XBC(V2,N1,N2+1)
C
C         END OF NONOSCILLATORY OPTION
C
    3                      CONTINUE
    6 CONTINUE
      RETURN
      END   

      SUBROUTINE XBC(X,N,M)
      DIMENSION X(N,M)
      DO 1 J=1,M
      X(1,J)=X(N-1,J)
    1 X(N,J)=X(2,J)
      RETURN
      END



