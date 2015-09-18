// $Id$ 

// Purpose: Implementation (declaration) of Wiscombe (1979, 1980, 1996) Mie scattering solutions and utilities 

/* Copyright (C) 2004--2014 Charlie Zender, Jorge Talamantes, Warren Wiscombe
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */


/* History:
   2003--2004: Jorge Talamantes translated Wiscombe's routines from Fortran77 to C++
   200402--: Charlie Zender merged Jorge's codes to libcsz_c++ 

   Reference articles: Wis79 (=Wis96), Wis80 

   Wiscombe, W. J., Mie Scattering Calculations: Advances in Technique and Fast, Vector-Speed Computer Codes, NCAR/TN-140+STR, National Center for Atmospheric Research, Boulder, Colo., 1979, edited/revised 1996.

   Wiscombe, W. J., Improved Mie scattering algorithms, Appl. Opt., 19(9), pp. 1505-1509, 1980. */

#include <mie_Wis79.hh> // Mie scattering solutions from Wis79

void BIGA ( dcomp CIOR, double XX, int NTRM, bool NOABS, bool YESANG,
            double RBIGA[], int RBIGA_dimension,
            dcomp CBIGA[], int CBIGA_dimension,
            bool *MSGLIM_thr,int *NUMMSG_errmsg_thr) // csz: Additional arguments for threading
{
/*   Calculate logarithmic derivatives of J-Bessel-function

     Input :  CIOR, XX, NTRM, NOABS, YESANG  (defined in MIEV0)

     Output :  RBIGA or CBIGA  (defined in MIEV0)

     Routines called :  CONFRA

     INTERNAL VARIABLES :

        CONFRA     Value of Lentz continued fraction for CBIGA(NTRM),
                      used to initialize downward recurrence

        DOWN       = True, use down-recurrence.  False, do not.

        F1,F2,F3   Arithmetic statement functions used in determining
                      whether to use up-  or down-recurrence
                      ( Ref. 2, Eqs. 6-8 )

        MRE        Real refractive index
        MIM        Imaginary refractive index

        REZINV     1 / ( MRE * XX ); temporary variable for recurrence
        ZINV       1 / ( CIOR * XX ); temporary variable for recurrence
*/

//    .. Local Scalars ..

  bool     DOWN;
  int      N;
  double   MIM, MRE, REZINV, RTMP;
  dcomp    CTMP, ZINV;

                               //  ** Decide whether BigA can be
                               //  ** calculated by up-recurrence
  MRE = real( CIOR );
  MIM = abs( imag( CIOR ) );

  if ( (MRE < 1.0) || (MRE > 10.0) || (MIM > 10.0) )

    DOWN = true;

  else if ( YESANG )
    {
      DOWN = true;
                                                 // ** Eq. R48
      if ( (MIM*XX) < F2( MRE ) ) DOWN = false;
    }
  else
    {
      DOWN = true;
                                                 // ** Eq. R48
      if ( (MIM*XX) < F1( MRE ) ) DOWN = false;
    }


  ZINV   = 1.0 / ( CIOR*XX );
  REZINV = 1.0 / ( MRE*XX );


  if ( DOWN )     //       ** Compute initial high-order BigA using
    {             //       ** Lentz method ( Ref. 1, pp. 17-20 )

      CTMP = CONFRA( NTRM, ZINV,
		     MSGLIM_thr,NUMMSG_errmsg_thr); // csz: Additional arguments for threading
		     

                                 // *** Downward recurrence for BigA
      if ( NOABS )
        {                        //      ** No-absorption case; Eq (R23)
          RBIGA[ NTRM ] = real( CTMP );

          for (N = NTRM; N >= 2; N--)
            RBIGA[ N - 1 ] = ( N*REZINV ) -
                             1.0 / ( ( N*REZINV ) + RBIGA[ N ] );
        }
      else
        {                        //      ** Absorptive case; Eq (R23)
          CBIGA[ NTRM ] = CTMP;

          for (N = NTRM; N >= 2; N--)
             CBIGA[ N-1 ] = (dcomp(N,0.)*ZINV) - 1.0 /
                             ( (dcomp(N,0.)*ZINV) + CBIGA[ N ] );
        }

    }
  else
    {             //                    *** Upward recurrence for BigA
      if ( NOABS )
        {         //               ** No-absorption case; Eq (R20,21)
          RTMP = sin ( MRE*XX );
          RBIGA[ 1 ] = - REZINV + RTMP /
                       ( RTMP*REZINV - cos ( MRE*XX ) );

          for (N = 2; N <= NTRM; N++)
            RBIGA[ N ] = -( N*REZINV ) +
                         1.0 / ( ( N*REZINV ) - RBIGA[ N - 1 ] );
         }
      else
         {        //                  ** Absorptive case; Eq (R20,22)
            CTMP = exp ( - dcomp(0.,2.)*CIOR*XX );
            CBIGA[ 1 ] = - ZINV + (1.-CTMP) /
                           ( ZINV * (1.-CTMP) - dcomp (0.,1.)*(1.+CTMP) );

            for (N = 2; N <= NTRM; N++)
               CBIGA[ N ] = - (dcomp(N,0.)*ZINV) +
                            1.0 / ((dcomp(N,0.)*ZINV) - CBIGA[ N-1 ]);
         }

    }

  return;
} /* end BIGA() */

void CKINMI ( // input variables:
              int NUMANG, int MAXANG, double XX, bool PERFCT,
              dcomp CREFIN, int MOMDIM, int NMOM, int IPOLZN,
              int ANYANG, double XMU[], int XMU_dimension,
              // output variables:
              bool CALCMO[], int CALCMO_dimension, int &NPQUAN,
              bool *MSGLIM_thr,int *NUMMSG_errmsg_thr,int *NUMMSG_wrtbad_thr) // csz: Additional arguments for threading
{
//       Check for bad input to MIEV0 and calculate CALCMO, NPQUAN

//       Routines called :  ERRMSG, WRTBAD, WRTDIM

//    .. Local Scalars ..

  int   N [ 4 ];
  bool  INPERR;
  int   I, IP, J, L;
  char  STRING [ 4 ];

  INPERR = false;

  if ( NUMANG > MAXANG ) INPERR = WRTDIM((char *)"MaxAng", NUMANG );
  if ( NUMANG < 0 ) INPERR = WRTBAD( (char *)"NUMANG",MSGLIM_thr,NUMMSG_errmsg_thr,NUMMSG_wrtbad_thr); // csz: Additional arguments for threading

  if ( XX < 0. ) INPERR = WRTBAD( (char *)"XX",MSGLIM_thr,NUMMSG_errmsg_thr,NUMMSG_wrtbad_thr); // csz: Additional arguments for threading

  if ( (! PERFCT) && ( real ( CREFIN ) <= 0.) )
                                  INPERR = WRTBAD( (char *)"CREFIN",MSGLIM_thr,NUMMSG_errmsg_thr,NUMMSG_wrtbad_thr); // csz: Additional arguments for threading

  if ( MOMDIM < 0 ) INPERR = WRTBAD( (char *)"MOMDIM",MSGLIM_thr,NUMMSG_errmsg_thr,NUMMSG_wrtbad_thr); // csz: Additional arguments for threading


  if ( NMOM != 0 )
    {

     if ( (NMOM < 0) || (NMOM > MOMDIM) ) INPERR = WRTBAD( (char *)"NMOM",MSGLIM_thr,NUMMSG_errmsg_thr,NUMMSG_wrtbad_thr); // csz: Additional arguments for threading

     if ( abs( IPOLZN ) > 4444 ) INPERR = WRTBAD( (char *)"IPOLZN",MSGLIM_thr,NUMMSG_errmsg_thr,NUMMSG_wrtbad_thr); // csz: Additional arguments for threading

     NPQUAN = 0;

     for (int L = 1; L <= 4; L++) CALCMO[ L ] = false;

     if ( IPOLZN != 0 )
       {                   //     ** Parse out IPOLZN into its digits
                           //     ** to find which phase quantities are
                           //     ** to have their moments calculated

         PARSEI ( IPOLZN, N, 4);

         for (int J = 0; J <= 3; J++)
           {
             IP = N [ J ];

             if ( (IP >= 1) &&  (IP <= 4) ) CALCMO[ IP ] = true;

             if ( (IP == 0) || ( (IP >= 5) && (IP <= 9) ) ) 
                                    INPERR = WRTBAD( (char *)"IPOLZN",MSGLIM_thr,NUMMSG_errmsg_thr,NUMMSG_wrtbad_thr); // csz: Additional arguments for threading

             if (NPQUAN < IP) NPQUAN = IP;
           }

       }

    }


  if ( ANYANG )
    {                        //  ** Allow for slight imperfections in
                             //  ** computation of cosine
      for (int I = 1; I <= NUMANG; I++)
        if ( (XMU[ I ] < -1.00001) || (XMU[ I ] > 1.00001) )
                INPERR = WRTBAD( (char *)"XMU",MSGLIM_thr,NUMMSG_errmsg_thr,NUMMSG_wrtbad_thr); // csz: Additional arguments for threading

    }
  else
    {
      for (int I = 1; I <= ( NUMANG + 1 )/2; I++)
        if ( (XMU[ I ] < -0.00001) || (XMU[ I ] > 1.00001) )
                INPERR = WRTBAD( (char *)"XMU",MSGLIM_thr,NUMMSG_errmsg_thr,NUMMSG_wrtbad_thr); // csz: Additional arguments for threading

    }


  if ( INPERR ) ERRMSG( (char *)"MIEV0--Input error(S).  Aborting...",
			true ,
			MSGLIM_thr,NUMMSG_errmsg_thr); // csz: Additional arguments for threading

  /* Split warnings into different categories to ease diagnosis */
  if(XX > 20000.0) ERRMSG( (char *)"MIEV0--size parameter > 20000, outside tested range",false,MSGLIM_thr,NUMMSG_errmsg_thr);

  if(real(CREFIN) > 10.0) ERRMSG( (char *)"MIEV0--idx_rfr_rl > 10, outside tested range",false,MSGLIM_thr,NUMMSG_errmsg_thr);

  if(abs(imag(CREFIN)) > 10.0) ERRMSG( (char *)"MIEV0--idx_rfr_img > 10, outside tested range",false,MSGLIM_thr,NUMMSG_errmsg_thr);

  /* if ( (XX > 20000.0) || (real( CREFIN ) > 10.0) ||
          (abs( imag( CREFIN ) ) > 10.0) )
          ERRMSG( (char *)"MIEV0--XX or CREFIN outside tested range",
		  false,
		  MSGLIM_thr,NUMMSG_errmsg_thr); */

  return;
} /* end CKINMI() */

dcomp CONFRA ( int N, dcomp ZINV,
	       bool *MSGLIM_thr,int *NUMMSG_errmsg_thr) // csz: Additional arguments for threading
{
/*        Compute Bessel function ratio A-sub-N from its
          continued fraction using Lentz method

          ZINV = Reciprocal of argument of A


     I N T E R N A L    V A R I A B L E S
     ------------------------------------

     CAK      Term in continued fraction expansion of A (Eq. R25)

     CAPT     Factor used in Lentz iteration for A (Eq. R27)

     CNUMER   Numerator   in capT  ( Eq. R28A )
     CDENOM   Denominator in capT  ( Eq. R28B )

     CDTD     Product of two successive denominators of capT factors
                  ( Eq. R34C )
     CNTN     Product of two successive numerators of capT factors
                  ( Eq. R34B )

     EPS1     Ill-conditioning criterion
     EPS2     Convergence criterion

     KK       Subscript k of cAk  ( Eq. R25B )

     KOUNT    Iteration counter ( used to prevent infinite looping )

     MAXIT    Max. allowed no. of iterations

     MM       + 1  and - 1, alternately
  --------------------------------------------------------------------
*/

//    .. Local Scalars ..

  int      KK, KOUNT, MM;
  dcomp    CAK, CAPT, CDENOM, CDTD, CNTN, CNUMER;

  double   EPS1 = 1.e-2, EPS2 = 1.E-8;
//double   EPS1 = 1.e-2, EPS2 = 1.E-10;
  int      MAXIT = 10000;

  dcomp    confra;

                               // ** Eq. R25a
  confra = ( dcomp ( N + 1 ) ) * ZINV;
  MM     = - 1;
  KK     = 2*N + 3;
                               // ** Eq. R25b, k=2
  CAK    = ( dcomp ( MM*KK ) ) * ZINV;
  CDENOM = CAK;
  CNUMER = CDENOM + 1.0 / confra;
  KOUNT  = 1;

  label10: KOUNT = KOUNT + 1;

  if ( KOUNT > MAXIT )
    ERRMSG((char *)"ConFra--Iteration failed to converge", true,
	   MSGLIM_thr,NUMMSG_errmsg_thr); // csz: Additional arguments for threading

  MM  = - MM;
  KK  = KK + 2;
                                               // ** Eq. R25b
  CAK = ( dcomp ( MM*KK ) ) * ZINV;
                                               // ** Eq. R32
  if ( ( abs( CNUMER / CAK ) <= EPS1 ) ||
       ( abs( CDENOM / CAK ) <= EPS1 ) )
    {
                                // ** Ill-conditioned case -- stride
                                // ** two terms instead of one

                                //      ** Eq. R34
       CNTN   = CAK * CNUMER + 1.0;
       CDTD   = CAK * CDENOM + 1.0;
                                //      ** Eq. R33
       confra = ( CNTN / CDTD ) * confra;

       MM  = - MM;
       KK  = KK + 2;
                                //      ** Eq. R25b
       CAK = ( dcomp ( MM*KK ) ) * ZINV;
                                //      ** Eq. R35
       CNUMER = CAK + CNUMER / CNTN;
       CDENOM = CAK + CDENOM / CDTD;
       KOUNT  = KOUNT + 1;
       goto label10;
    }
  else
    {           //                      *** Well-conditioned case

                                //       ** Eq. R27
       CAPT   = CNUMER / CDENOM;
                                //       ** Eq. R26
       confra = CAPT * confra;
                                //      ** Check for convergence; Eq. R31

       if (    ( abs( real (CAPT) - 1.0 ) >= EPS2 )  ||
               ( abs( imag (CAPT)       ) >= EPS2 ) )
         {                      //       ** Eq. R30
            CNUMER = CAK + 1.0 / CNUMER;
            CDENOM = CAK + 1.0 / CDENOM;

            goto label10;
          }

    }


  return confra;
} /* end CONFRA() */

void ERRMSG( char *MESSAG, bool FATAL,
	     bool *MSGLIM_thr,int *NUMMSG_errmsg_thr) // csz: Additional arguments for threading
{
//       Print out a warning or error message;  abort if error
//       after making symbolic dump (machine-specific)

//    .. Local Scalars ..

  // static bool  MSGLIM = false; // csz: Replace static value with local argument
  // static int NUMMSG = 0; // csz: Replace static value with local argument
  const int MAXMSG=100; // csz: Replace static with const

  if ( FATAL )
    {
      cout << endl << endl << " ****** ERROR *****  " << MESSAG << endl;
      exit(0);
    }


  *NUMMSG_errmsg_thr = *NUMMSG_errmsg_thr + 1;

  if ( *MSGLIM_thr ) return;

  if ( *NUMMSG_errmsg_thr <= MAXMSG )

    cout << " ****** WARNING *****  #" << *NUMMSG_errmsg_thr << ": " << MESSAG << endl;

  else
    {
      cout << endl << endl << " ****** TOO MANY WARNING MESSAGES -- ";
      cout << "They will no longer be printed *******" << endl;
      *MSGLIM_thr = true;
    }


  return;
} /* end ERRMSG() */

double F1 ( double MRE)
  {                                                  // ** Eq. R47c
    return ( -8.0 + pow(MRE, 2.)*( 26.22 +
               MRE*( -0.4474 + pow(MRE, 3.)*( 0.00204 - 0.000175*MRE ) ) ) );
  } /* end F1() */

double F2 ( double MRE)
  {                                                  // ** Eq. R47b
    return ( 3.9 + MRE*( -10.8 + 13.78*MRE ) );
  } /* end F2() */

double F3 ( double MRE)
  {                                                  // ** Eq. R47a
    return ( -15.04 + MRE*( 8.42 + 16.35*MRE ) );
  } /* end F3() */

void LPCO1T ( int NMOM, int IPOLZN, int MOMDIM,
              bool CALCMO[], int CALCMO_dimension,
              dcomp A[], int A_dimension, dcomp B[], int B_dimension,
              double PMOM[][PMOM_second_dimension], int PMOM_dimension )
{
/*        Calculate Legendre polynomial expansion coefficients (also
          called moments) for phase quantities in special case where
          no. terms in Mie series = 1

         INPUT:  NMOM, IPOLZN, MOMDIM     MIEV0 arguments
                 CALCMO                   Flags calculated from IPOLZN
                 A(1), B(1)               Mie series coefficients

         OUTPUT: PMOM                     Legendre moments
*/

//    .. Local Scalars ..

  int     L, NUMMOM;
  double  A1SQ, B1SQ;
  dcomp   A1B1C, CTMP;


  A1SQ   = SQ( A[ 1 ] );
  B1SQ   = SQ( B[ 1 ] );
  A1B1C  = A[ 1 ] * conj( B[ 1 ] );

  if ( IPOLZN < 0 )
    { 
    if ( CALCMO[ 1 ] ) PMOM[ 0 ][ 1 ] = 2.25*B1SQ;

    if ( CALCMO[ 2 ] ) PMOM[ 0 ][ 2 ] = 2.25*A1SQ;

    if ( CALCMO[ 3 ] ) PMOM[ 0 ][ 3 ] = 2.25*real( A1B1C );

    if ( CALCMO[ 4 ] ) PMOM[ 0 ][ 4 ] = 2.25*imag( A1B1C );
    }
  else
    {
    if (NMOM >= 2) NUMMOM = 2; else NUMMOM = NMOM;

    label10: for (int L = 0; L <= NUMMOM; L++)  // ** Loop over moments
      {
            if ( IPOLZN == 0 )
              {
               if ( L == 0 ) PMOM [ L ][ 1 ] = 1.5*( A1SQ + B1SQ );

               if ( L == 1 ) PMOM [ L ][ 1 ] = 1.5*real( A1B1C );

               if ( L == 2 ) PMOM [ L ][ 1 ] = 0.15*( A1SQ + B1SQ );

              // continue; // use this instead of backwards goto
              goto label10;
              }


            if ( CALCMO[ 1 ] )
              {
               if ( L == 0 ) PMOM [ L ][ 1 ] = 2.25*( A1SQ + B1SQ / 3.);

               if ( L == 1 ) PMOM [ L ][ 1 ] = 1.5*real( A1B1C );

               if ( L == 2 ) PMOM [ L ][ 1 ] = 0.3*B1SQ;
              }


            if ( CALCMO[ 2 ] )
              {
               if ( L == 0 ) PMOM [ L ][ 2 ] = 2.25*( B1SQ + A1SQ / 3. );

               if ( L == 1 ) PMOM [ L ][ 2 ] = 1.5*real( A1B1C );

               if ( L == 2 ) PMOM [ L ][ 2 ] = 0.3*A1SQ;
              }


            if ( CALCMO[ 3 ] )
              {
               if ( L == 0 ) PMOM [ L ][ 3 ] = 3.0*real( A1B1C );

               if ( L == 1 ) PMOM [ L ][ 3 ] = 0.75*( A1SQ + B1SQ );

               if ( L == 2 ) PMOM [ L ][ 3 ] = 0.3*real( A1B1C );
              }


            if ( CALCMO[ 4 ] )
              {
               if ( L == 0 ) PMOM [ L ][ 4 ] = -1.5*imag( A1B1C );

               if ( L == 1 ) PMOM [ L ][ 4 ] = 0.0;

               if ( L == 2 ) PMOM [ L ][ 4 ] = 0.3*imag( A1B1C );
              }

    }

  }

  return;
} /* end LPCO1T() */

void LPCO2T ( int NMOM, int IPOLZN, int MOMDIM,
              bool CALCMO[], int CALCMO_dimension,
              dcomp A[], int A_dimension, dcomp B[], int B_dimension,
              double PMOM[][PMOM_second_dimension], int PMOM_dimension )
{
/*        Calculate Legendre polynomial expansion coefficients (also
          called moments) for phase quantities in special case where
          no. terms in Mie series = 2

         INPUT:  NMOM, IPOLZN, MOMDIM     MIEV0 arguments
                 CALCMO                   Flags calculated from IPOLZN
                 A(1-2), B(1-2)           Mie series coefficients

         OUTPUT: PMOM                     Legendre moments
*/

//    .. Local Scalars ..

  int     L, NUMMOM;
  double  A2SQ, B2SQ, PM1, PM2;
  dcomp   A2C, B2C, CA, CAC, CAT, CB, CBC, CBT, CG, CH, CTMP;


  CA   = 3.*A[ 1 ] - 5.*B[ 2 ];
  CAT  = 3.*B[ 1 ] - 5.*A[ 2 ];
  CAC  = conj( CA );
  A2SQ = SQ( A[ 2 ] );
  B2SQ = SQ( B[ 2 ] );
  A2C  = conj( A[ 2 ] );
  B2C  = conj( B[ 2 ] );

  if ( IPOLZN < 0 )
    {                                   //** Loop over Sekera moments
      if (NMOM < 2) NUMMOM = NMOM;
      else NUMMOM = 2;

      for (int L = 0; L <= NUMMOM; L++)
        {
          if ( CALCMO[ 1 ] )
            {
               if ( L == 0 ) PMOM [ L ][ 1 ] = 0.25 * ( SQ( CAT )
                                              + (100./3.)* B2SQ );

               if ( L == 1 ) PMOM [ L ][ 1 ] = (5./3.)*real( CAT*B2C );

               if ( L == 2 ) PMOM [ L ][ 1 ] = (10./3.)*B2SQ;
            }

          if ( CALCMO[ 2 ] )
            {
               if ( L == 0 ) PMOM [ L ][ 2 ] = 0.25 * ( SQ( CA )
                                              + (100./3.) * A2SQ );

               if ( L == 1 ) PMOM [ L ][ 2 ] = (5./3.)*real( CA*A2C );

               if ( L == 2 ) PMOM [ L ][ 2 ] = (10./3.)*A2SQ;
            }

          if ( CALCMO[ 3 ] )
            {
               if ( L == 0 ) PMOM [ L ][ 3 ] = 0.25 * real( CAT * CAC
                                              + (100./3.) * B[2] * A2C );

               if ( L == 1 ) PMOM [ L ][ 3 ] = 5./6.*
                                              real( B[2]*CAC + CAT*A2C );

               if ( L == 2 ) PMOM [ L ][ 3 ] = 10./3.* real( B[2]*A2C );
            }

          if ( CALCMO[ 4 ] )
            {
               if ( L == 0 ) PMOM [ L ][ 4 ] = -0.25 * imag( CAT * CAC
                                              + (100./3.)* B[2] * A2C );

               if ( L == 1 ) PMOM [ L ][ 4 ] = -5./ 6.*
                                               imag( B[2]*CAC + CAT*A2C );

	       if ( L == 2 ) PMOM [ L ][ 4 ] = -10./ 3.* imag( B[2]*A2C );
            }

        }
    }
  
  else

    {

      CB  = 3.*B[ 1 ] + 5.*A[ 2 ];
      CBT = 3.*A[ 1 ] + 5.*B[ 2 ];
      CBC = conj( CB );
      CG  = ( CBC*CBT + 10.*( CAC*A[ 2 ] + B2C*CAT ) ) / 3.;
      CH  = 2.*( CBC*A[ 2 ] + B2C*CBT );

      if (NMOM < 4) NUMMOM = NMOM;
      else NUMMOM = 4;

                                        // ** Loop over Mueller moments

      for (int L = 0; L <=NUMMOM; L++)
        {

          if ( IPOLZN == 0 || CALCMO[ 1 ] )
            {
               if ( L == 0 ) PM1 = 0.25*SQ( CA ) + SQ( CB ) / 12.
                                  + (5./3.)*real( CA*B2C ) + 5.*B2SQ;

               if ( L == 1 ) PM1 = real( CB * ( CAC / 6.+ B2C ) );

               if ( L == 2 ) PM1 = SQ( CB ) / 30.+ (20./7.)*B2SQ
                                  + (2./3.)*real( CA*B2C );

               if ( L == 3 ) PM1 = (2./7.) * real( CB*B2C );

               if ( L == 4 ) PM1 = (40./63.) * B2SQ;

               if ( CALCMO[ 1 ] ) PMOM [ L ][ 1 ] = PM1;
            }

          if ( IPOLZN == 0 || CALCMO[ 2 ] )
            {
               if ( L == 0 ) PM2 = 0.25*SQ( CAT ) + SQ( CBT ) / 12.
                                  + ( 5./ 3.) * real( CAT*A2C )
                                  + 5.*A2SQ;

               if ( L == 1 ) PM2 = real( CBT *
                                       ( conj( CAT ) / 6.+ A2C ) );

               if ( L == 2 ) PM2 = SQ( CBT ) / 30.
                                  + ( 20./7.) * A2SQ
                                  + ( 2./3.) * real( CAT*A2C );

               if ( L == 3 ) PM2 = (2./7.) * real( CBT*A2C );

               if ( L == 4 ) PM2 = (40./63.) * A2SQ;

               if ( CALCMO[ 2 ] ) PMOM [ L ][ 2 ] = PM2;
            }

          if ( IPOLZN == 0 )
            {
              PMOM [ L ][ 1 ] = 0.5*( PM1 + PM2 );
              continue;
            }

          if ( CALCMO[ 3 ] )
            {
               if ( L == 0 ) PMOM [ L ][ 3 ] = 0.25 * real( CAC*CAT + CG
                                               + 20.* B2C * A[2] );

               if ( L == 1 ) PMOM [ L ][ 3 ] = real( CAC*CBT + CBC*CAT
                                                + 3.*CH ) / 12.;

               if ( L == 2 ) PMOM [ L ][ 3 ] = 0.1 * real( CG
                                                + (200./7.) * B2C * A[2] );

               if ( L == 3 ) PMOM [ L ][ 3 ] = real( CH ) / 14.;

               if ( L == 4 ) PMOM [ L ][ 3 ] = 40./63.* real( B2C*A[2] );
            }

          if ( CALCMO[ 4 ] )
            {
               if ( L == 0 ) PMOM [ L ][ 4 ] = 0.25 * imag( CAC*CAT + CG
                                              + 20.* B2C * A[2] );

               if ( L == 1 ) PMOM [ L ][ 4 ] = imag( CAC*CBT + CBC*CAT
                                                 + 3.*CH ) / 12.;

               if ( L == 2 ) PMOM [ L ][ 4 ] = 0.1 * imag( CG
                                              + (200./7.) * B2C * A[2] );

               if ( L == 3 ) PMOM [ L ][ 4 ] = imag( CH ) / 14.;

               if ( L == 4 ) PMOM [ L ][ 4 ] = 40./63.* imag( B2C*A[2] );
            }

      }

    }

  return;
} /* end LPCO2T() */

void LPCOEF( int NTRM, int NMOM, int IPOLZN, int MOMDIM, bool CALCMO[],
             int CALCMO_dimension, int NPQUAN,
             dcomp A[], int A_dimension, dcomp B[], int B_dimension,
             double PMOM[][PMOM_second_dimension], int PMOM_dimension,
	     bool *PASS1_lpcoef_thr,double *RECIP_thr, // csz: Additional arguments for threading
             bool *MSGLIM_thr,int *NUMMSG_errmsg_thr) // csz: Additional arguments for threading
{

/*        Calculate Legendre polynomial expansion coefficients (also
          called moments) for phase quantities ( Ref. 5 formulation )

      INPUT:  NTRM                    Number terms in Mie series
              NMOM, IPOLZN, MOMDIM    MIEV0 arguments
              CALCMO                  Flags calculated from IPOLZN
              NPQUAN                  Defined in MIEV0
              A, B                    Mie series coefficients

      OUTPUT: PMOM                   Legendre moments (MIEV0 argument)

      Routines called :  ERRMSG, LPCO1T, LPCO2T

      *** NOTES ***

          (1)  Eqs. 2-5 are in error in Dave, Appl. Opt. 9,
          1888 (1970).  Eq. 2 refers to M1, not M2;  eq. 3 refers to
          M2, not M1.  In eqs. 4 and 5, the subscripts on the second
          term in square brackets should be interchanged.

          (2)  The general-case logic in this subroutine works correctly
          in the two-term Mie series case, but subroutine LPCO2T
          is called instead, for speed.

          (3)  Subroutine  LPCO1T, to do the one-term case, is never
          called within the context of MIEV0, but is included for
          complete generality.

          (4)  Some improvement in speed is obtainable by combining the
          310- and 410-loops, if moments for both the third and fourth
          phase quantities are desired, because the third phase quantity
          is the real part of a complex series, while the fourth phase
          quantity is the imaginary part of that very same series.  But
          most users are not interested in the fourth phase quantity,
          which is related to circular polarization, so the present
          scheme is usually more efficient.


            ** Definitions of local variables ***

       AM(M)       Numerical coefficients  a-sub-m-super-l
                      in Dave, Eqs. 1-15, as simplified in Ref. 5.

       BI(I)       Numerical coefficients  b-sub-i-super-l
                      in Dave, Eqs. 1-15, as simplified in Ref. 5.

       BIDEL(I)    1/2 Bi(I) times factor capital-del in Dave

       CM,DM()     Arrays C and D in Dave, Eqs. 16-17 (Mueller form),
                      calculated using recurrence derived in Ref. 5

       CS,DS()     Arrays C and D in Ref. 4, Eqs. A5-A6 (Sekera form),
                      calculated using recurrence derived in Ref. 5

       C,D()       Either CM,DM or CS,DS, depending on IPOLZN

       EVENL       True for even-numbered moments;  false otherwise

       IDEL        1 + little-del  in Dave

       MAXTRM      Max. no. of terms in Mie series

       MAXMOM      Max. no. of non-zero moments

       NUMMOM      Number of non-zero moments

       RECIP(K)    1 / K
*/

//    .. Parameters ..

  const int    MAXTRM = 50000; // csz 50000 tried for large snow particles
  // const int    MAXTRM = 10100; // csz 10100 used for Wiscombe's test cases
  // const int MAXTRM = 1102; // csz MAXTRM should be consistent, or, more precisely, less than or equal to MAXTRM in MIEV0()=mie_sph_Wis79()
  const int MAXMOM = 2*MAXTRM;
  const int MXMOM2 = MAXMOM / 2;
  const int MAXRCP = 4*MAXTRM + 2;

//    .. Local Scalars ..

  bool         EVENL;
  // static bool  PASS1 = true; // csz: Replace static value with local argument
  int          I, IDEL, IMAX, J, K, L, LD2, M, MMAX, NUMMOM;
  double       SUM;
  int          comp;

//    .. Local Arrays ..

  double          AM [ MAXTRM + 1 ], BI [ MXMOM2 + 1 ], BIDEL [ MXMOM2  + 1 ];
  // fxm: got to here in threading
  //  static double   RECIP [ MAXRCP + 1 ]; // csz: Replace static value with local argument
  dcomp           C [ MAXTRM + 1 ], *CM = C, // C & CM are EQUIVALENT
                  D [ MAXTRM + 1 ], *DM = D,
                  CS [ MAXTRM + 1 ], DS[ MAXTRM + 1 ];


  if ( *PASS1_lpcoef_thr )
    {
      for (int K = 1; K <= MAXRCP; K++) RECIP_thr [ K ] = 1.0 / K;

      *PASS1_lpcoef_thr  = false;
    }

  if (1 > NPQUAN) comp = 1;
  else comp = NPQUAN;

  for (int J = 1; J <=  comp; J++)
    {
      for (int L = 0; L <= NMOM; L++) PMOM [L][J] = 0.0;
    }

  if ( NTRM == 1 )
    {
      LPCO1T ( NMOM, IPOLZN, MOMDIM, CALCMO, CALCMO_dimension,
              A, A_dimension, B, B_dimension, PMOM, PMOM_dimension );
      return;
    }
  else if ( NTRM == 2 )
    {
      LPCO2T ( NMOM, IPOLZN, MOMDIM, CALCMO, CALCMO_dimension,
              A, A_dimension, B, B_dimension, PMOM, PMOM_dimension );
      return;
    }

    if ( (NTRM + 2) > MAXTRM )
      ERRMSG((char *)"LPCoef--PARAMETER MaxTrm too small (csz: HINT: set MaxTrm >~ sz_prm)", true,
	     MSGLIM_thr,NUMMSG_errmsg_thr); // csz: Additional arguments for threading

                                   // ** Calculate Mueller C, D arrays
      CM[ NTRM + 2 ] = dcomp ( 0., 0. );
      DM[ NTRM + 2 ] = dcomp ( 0., 0. );
      CM[ NTRM + 1 ] = ( 1. - RECIP_thr[ NTRM+1 ] ) * B[ NTRM ];
      DM[ NTRM + 1 ] = ( 1. - RECIP_thr[ NTRM+1 ] ) * A[ NTRM ];
      CM[ NTRM ] = ( RECIP_thr[ NTRM ] + RECIP_thr[ NTRM+1 ] ) * A[ NTRM ] +
                   ( 1. - RECIP_thr[ NTRM ] )*B[ NTRM-1 ];
      DM[ NTRM ] = ( RECIP_thr[ NTRM ] + RECIP_thr[ NTRM+1 ] ) * B[ NTRM ] +
                   ( 1. - RECIP_thr[ NTRM ] )*A[ NTRM-1 ];

      for (int K = NTRM-1; K >=2; K--)
        {
         CM[ K ] = CM[ K+2 ] - ( 1. + RECIP_thr[K+1] ) * B[ K+1 ]
                             + ( RECIP_thr[K] + RECIP_thr[K+1] ) * A[ K ]
                             + ( 1. - RECIP_thr[K] ) * B[ K-1 ];
         DM[ K ] = DM[ K+2 ] - ( 1. + RECIP_thr[K+1] ) * A[ K+1 ]
                             + ( RECIP_thr[K] + RECIP_thr[K+1] ) * B[ K ]
                             + ( 1. - RECIP_thr[K] ) * A[ K-1 ];
        }

      CM[ 1 ] = CM[ 3 ] + 1.5 * ( A[ 1 ] - B[ 2 ] );
      DM[ 1 ] = DM[ 3 ] + 1.5 * ( B[ 1 ] - A[ 2 ] );


      if ( IPOLZN >= 0 )
        {
         for (int K = 1; K <= (NTRM + 2); K++)
           {
            C[ K ] = static_cast <double> ( 2*K - 1 ) * CM[ K ];
            D[ K ] = static_cast <double> ( 2*K - 1 ) * DM[ K ];
           }
        }
      else
        {                       //   ** Compute Sekera C and D arrays
         CS[ NTRM + 2 ] = dcomp ( 0., 0. );
         DS[ NTRM + 2 ] = dcomp ( 0., 0. );
         CS[ NTRM + 1 ] = dcomp ( 0., 0. );
         DS[ NTRM + 1 ] = dcomp ( 0., 0. );

         for (int K = NTRM; K >= 1; K--)
           {
            CS[ K ] = CS[ K+2 ] +
                      static_cast <double> ( 2*K + 1 ) * ( CM[ K+1 ] - B[ K ] );
            DS[ K ] = DS[ K+2 ] +
                      static_cast <double> ( 2*K + 1 ) * ( DM[ K+1 ] - A[ K ] );
           }

         for (int K = 1; K <= (NTRM + 2); K++)
           {
            C[ K ] = static_cast <double> ( 2*K - 1 ) * CS[ K ];
            D[ K ] = static_cast <double> ( 2*K - 1 ) * DS[ K ];
           }
        }


      if ( IPOLZN <  0 ) if (NMOM < (2*NTRM - 2)) NUMMOM = NMOM;
      else NUMMOM = 2*NTRM -2;

      if ( IPOLZN >= 0 ) if (NMOM < 2*NTRM ) NUMMOM = NMOM;
      else NUMMOM = 2*NTRM;

      if ( NUMMOM > MAXMOM )
	ERRMSG((char *)"LPCoef--PARAMETER MaxTrm too small", true,
	       MSGLIM_thr,NUMMSG_errmsg_thr); // csz: Additional arguments for threading


                       //  ** Loop over moments

      for (int L = 0; L <= NUMMOM; L++)
        {
         LD2 = L / 2;
         EVENL  = ( L % 2 ) == 0;
                                  // ** Calculate numerical coefficients
                                  // ** a-sub-m and b-sub-i in Dave
                                  // ** double-sums for moments
         if ( L == 0 )
           {
            IDEL = 1;

            for (int M = 0; M <= NTRM; M++)
               AM[ M ] = 2.0 * RECIP_thr[ 2*M + 1 ];

            BI[ 0 ] = 1.0;
           }
         else if ( EVENL )
           {
            IDEL = 1;

            for ( int M = LD2; M <= NTRM; M++)
               AM[ M ] = ( 1. + RECIP_thr[ 2*M - L + 1 ] ) * AM[ M ];

            for (int I = 0; I <= (LD2 - 1); I++)
               BI[ I ] = ( 1. - RECIP_thr[ L - 2*I ] ) * BI[ I ];

            BI[ LD2 ] = ( 2. - RECIP_thr[ L ] ) * BI[ LD2 - 1 ];
           }
         else
           {
            IDEL = 2;

            for ( int M = LD2; M <= NTRM; M++)
               AM[ M ] = ( 1. - RECIP_thr[ 2*M + L + 2 ] ) * AM[ M ];

            for (int I = 0; I <= LD2; I++)
               BI[ I ] = ( 1. - RECIP_thr[ L + 2*I + 1 ] ) * BI[ I ];

           }
                                  //  ** Establish upper limits for sums
                                  //  ** and incorporate factor capital-
                                  //  ** del into b-sub-i
         MMAX = NTRM - IDEL;
         if ( IPOLZN >= 0 ) MMAX = MMAX + 1;

         if (LD2 < MMAX - LD2) IMAX = LD2;
         else IMAX = MMAX - LD2;

         if ( IMAX < 0 ) goto label250;

         for (int I = 0; I <= IMAX; I++) BIDEL[ I ] = BI[ I ];

         if ( EVENL ) BIDEL[ 0 ] = 0.5*BIDEL[ 0 ];

                                 //  ** Perform double sums just for
                                 //  ** phase quantities desired by user
         if ( IPOLZN == 0 )
           {
            for (int I = 0; I <=IMAX; I++)
              {                        //   ** vectorizable loop

               SUM = 0.0;

               for (int M = LD2; M <= (MMAX - I); M++)
                 {
                  SUM = SUM + AM[ M ] *
                            ( real( C[M-I+1] * conj( C[M+I+IDEL] ) )
                            + real( D[M-I+1] * conj( D[M+I+IDEL] ) ) );
                 }

               PMOM [ L ][ 1 ] = PMOM [ L ][ 1 ] + BIDEL[ I ] * SUM;

              }

            PMOM [ L ][ 1 ] = 0.5 * PMOM [ L ][ 1 ];
            continue;

           }


         if ( CALCMO[ 1 ] )
           {
            for (int I = 0; I <=IMAX; I++)
              {
               SUM = 0.0;
                                      //    ** vectorizable loop
               for (int M = LD2; M <= (MMAX - I); M++)
                  SUM = SUM + AM[ M ] *
                              real( C[M-I+1] * conj( C[M+I+IDEL] ) );

               PMOM [ L ][ 1 ] = PMOM [ L ][ 1 ] + BIDEL[ I ] * SUM;
              }

           }

         if ( CALCMO[ 2 ] )
           {
            for (int I = 0; I <= IMAX; I++)
              {
               SUM = 0.0;
                                       //   ** vectorizable loop
               for (int M = LD2; M <= (MMAX - I); M++)
                  SUM = SUM + AM[ M ] *
                              real( D[M-I+1] * conj( D[M+I+IDEL] ) );

               PMOM [ L ][ 2 ] = PMOM [ L ][ 2 ] + BIDEL[ I ] * SUM;
              }
           }

         if ( CALCMO[ 3 ] )
           {
            for (int I = 0; I <= IMAX; I++)
              {
               SUM = 0.0;
                                      //    ** vectorizable loop
               for (int M = LD2; M <= (MMAX - I); M++)
                  SUM = SUM + AM[ M ] *
                            ( real( C[M-I+1] * conj( D[M+I+IDEL] ) )
                            + real( C[M+I+IDEL] * conj( D[M-I+1] ) ) );

               PMOM [ L ][ 3 ] = PMOM [ L ][ 3 ] + BIDEL[ I ] * SUM;
              }

            PMOM [ L ][ 3 ] = 0.5 * PMOM [ L ][ 3 ];
           }

         if ( CALCMO[ 4 ] )
           {
            for (int I = 0; I <= IMAX; I++)
              {
               SUM= 0.0;
                                      //    ** vectorizable loop
               for (int M = LD2; M <= (MMAX - I); M++)
                  SUM = SUM + AM[ M ] *
                            ( imag( C[M-I+1] * conj( D[M+I+IDEL] ) )
                            + imag( C[M+I+IDEL] * conj( D[M-I+1] ) ) );

               PMOM [ L ][ 4 ] = PMOM [ L ][ 4 ] + BIDEL[ I ] * SUM;
              }

            PMOM [ L ][ 4 ] = - 0.5 * PMOM [ L ][ 4 ];

           }

        }

  label250: return;
} /* end LPCOEF() */

double MAX (double first, double second)
{
  if (first >= second) return first; else return second;
} /* end MAX() */

int mie_sph_Wis79(  // [fnc] Mie solution for homogeneous spheres (MIEV0), Wis79
// input variables:
       double XX, dcomp CREFIN, bool PERFCT, double MIMCUT,
       bool ANYANG, int NUMANG, double XMU[], int XMU_dimension,
       int NMOM, int IPOLZN, int MOMDIM, bool PRNT[], int PRNT_dimension,
// output variables:
       double &QEXT, double &QSCA, double &GQSC,
       double PMOM[][PMOM_second_dimension], int PMOM_dimension, dcomp &SFORW,
       dcomp &SBACK, dcomp S1[], int S1_dimension, dcomp S2[], int S2_dimension,
       dcomp TFORW[], int TFORW_dimension, dcomp TBACK[], int TBACK_dimension,
       double &SPIKE )
{

/*   Computes Mie scattering and extinction efficiencies; asymmetry
     factor;  forward- and backscatter amplitude;  scattering
     amplitudes vs. scattering angle for incident polarization parallel
     and perpendicular to the plane of scattering;
     coefficients in the Legendre polynomial expansions of either the
     unpolarized phase function or the polarized phase matrix;
     some quantities needed in polarized radiative transfer;  and
     information about whether or not a resonance has been hit.

     Input and output variables are described in file MIEV.doc. 
     Many statements are accompanied by comments referring to 
     references in MIEV.doc, notably the NCAR Mie report which is now
     available electronically and which is referred to using the
     shorthand (Rn), meaning Eq. (n) of the report.

     CALLING TREE:

         MIEV0
             TESTMI
                 TSTBAD
                 MIPRNT
                 ERRMSG
             CKINMI
                 WRTBAD
                 WRTDIM
                 ERRMSG
             SMALL1
             SMALL2
             ERRMSG
             BIGA
                 CONFRA
                     ERRMSG
             LPCOEF
                 LPCO1T
                 LPCO2T
                 ERRMSG
             MIPRNT


       I N T E R N A L   V A R I A B L E S
       -----------------------------------

   AN,BN           Mie coefficients a-sub-n, b-sub-n ( Ref. 1, Eq. 16 )

   ANM1,BNM1       Mie coefficients  a-sub-(n-1),
                      b-sub-(n-1);  used in GQSC sum

   ANP             Coeffs. in S+ expansion ( Ref. 2, p. 1507 )
   BNP             Coeffs. in S- expansion ( Ref. 2, p. 1507 )
   ANPM            Coeffs. in S+ expansion ( Ref. 2, p. 1507 )
                      when  MU  is replaced by  - MU
   BNPM            Coeffs. in S- expansion ( Ref. 2, p. 1507 )
                      when  MU  is replaced by  - MU

   CALCMO(K)       TRUE, calculate moments for K-th phase quantity
                      (derived from IPOLZN)

   CBIGA(N)        Bessel function ratio A-sub-N (Ref. 2, Eq. 2)
                      ( COMPLEX version )

   CDENAN,         (COMPLEX) denominators of An,Bn
    CDENBN

   CIOR            Complex index of refraction with negative
                      imaginary part (Van de Hulst convention)
   CIORIV          1 / cIoR

   COEFF           ( 2N + 1 ) / ( N ( N + 1 ) )

   CSUM1,2         temporary sum variables for TFORW, TBACK

   FN              Floating point version of loop index for
                      Mie series summation

   LITA,LITB(N)    Mie coefficients An, Bn, saved in arrays for
                      use in calculating Legendre moments PMOM

   MAXTRM          Max. possible no. of terms in Mie series

   MM              (-1)^(n+1), where n is Mie series sum index 

   MIM             Magnitude of imaginary refractive index
   MRE             Real part of refractive index

   MAXANG          Max. possible value of input variable NUMANG
   NANGD2          (NUMANG+1)/2 ( no. of angles in 0-90 deg; ANYANG=F )

   NOABS           TRUE, sphere non-absorbing (determined by MIMCUT)

   NP1DN           ( N + 1 ) / N

   NPQUAN          Highest-numbered phase quantity for which moments are
                      to be calculated (the largest digit in IPOLZN
                      if  IPOLZN .NE. 0)

   NTRM            No. of terms in Mie series

   PASS1           TRUE on first entry, FALSE thereafter; for self-test

   PIN(J)          Angular function pi-sub-n ( Ref. 2, Eq. 3 )
                      at J-th angle
   PINM1(J)        pi-sub-(n-1) ( see PIn ) at J-th angle

   PSINM1          Ricatti-Bessel function psi-sub-(n-1), argument XX
   PSIN            Ricatti-Bessel function psi-sub-n of argument XX
                      ( Ref. 1, p. 11 ff. )

   RBIGA(N)        Bessel function ratio A-sub-N (Ref. 2, Eq. 2)
                      ( REAL version, for when imag refrac index = 0 )

   RIORIV          1 / Mre

   RN              1 / N

   RTMP            (REAL) temporary variable

   SP(J)           S+  for J-th angle  ( Ref. 2, p. 1507 )
   SM(J)           S-  for J-TH angle  ( Ref. 2, p. 1507 )
   SPS(J)          S+  for (NUMANG+1-J)-th angle ( ANYANG=FALSE )
   SMS(J)          S-  for (NUMANG+1-J)-th angle ( ANYANG=FALSE )

   TAUN            Angular function tau-sub-n ( Ref. 2, Eq. 4 )
                      at J-th angle

   TCOEF           N ( N+1 ) ( 2N+1 ) (for summing TFORW,TBACK series)

   TWONP1          2N + 1

   YESANG          TRUE if scattering amplitudes are to be calculated

   ZETNM1          Ricatti-Bessel function  zeta-sub-(n-1) of argument
                      XX  ( Ref. 2, Eq. 17 )
   ZETN            Ricatti-Bessel function  zeta-sub-n of argument XX
  ----------------------------------------------------------------------
*/
/*                                 ** NOTE --  MAXTRM = 10100  is neces-
                                   ** sary to do some of the test probs,
                                   ** but 1100 is sufficient for most
                                   ** conceivable applications
*/

  // csz++
  int rcd(0); // [rcd] Return code
  int NUMMSG_wrtbad_thr(0); // csz: Replace static value with local argument
  // fxm: 20070924 static NUMMSG_errmsg_thr may break threading
  static int NUMMSG_errmsg_thr(0); // csz: Replace static value with local argument
  static bool MSGLIM_thr(false); // csz: Replace static value with local argument
  bool PASS1_lpcoef_thr(true); // csz: Replace static value with local argument
  // LPCOEF will likely fail if sz_prm > MAXTRM
  /* 20040711: 1102 terms is too small for large snow particles */
  const int    MAXTRM(50000); // csz 50000 tried for large snow particles
  // const int    MAXTRM(10100); // csz 10100 used for Wiscombe's test cases
  // const int MAXTRM(1102); // csz LPCOEF uses 1102 by default---be consistent!
  const int MAXRCP(4*MAXTRM+2);
  double   RECIP_thr[MAXRCP+1]; // csz: Replace static value with local argument
  // Values for TESTMI
  bool PERSAV_thr, ANYSAV_thr; // csz: Replace static value with local argument
  double XXSAV_thr, MIMSAV_thr, XMUSAV_thr; // csz: Replace static value with local argument
  dcomp CRESAV_thr; // csz: Replace static value with local argument
  int NMOSAV_thr, IPOSAV_thr, NUMSAV_thr; // csz: Replace static value with local argument
  // csz--

//    .. Parameters ..

//  const int    MAXANG = 501;
// csz Increase MAXANG to 721 for quarter-degree resolution
  const int    MAXANG = 721;
  const int    MXANG2 = MAXANG / 2 + 1;
  const double ONETHR = 1. / 3.;

//    .. Local Scalars ..

  bool        NOABS, YESANG;
  bool PASS1_miev0_thr = true; // csz: Replace static value with local argument
  int         I, J, N, NANGD2, NPQUAN, NTRM;
  double      CHIN, CHINM1, COEFF, DENAN, DENBN, FN, MIM, MM, MRE,
              NP1DN, PSIN, PSINM1, RATIO, RIORIV, RN, RTMP, TAUN,
              TCOEF, TWONP1, XINV;
  dcomp       AN, ANM1, ANP, ANPM, BN, BNM1, BNP, BNPM, CDENAN,
              CDENBN, CIOR, CIORIV, CSUM1, CSUM2, CTMP, ZET, 
              ZETN, ZETNM1;

//    .. Local Arrays ..

  const int     CALCMO_dimension = 4 + 1;

  const int     PIN_dimension    = MAXANG + 1;
  const int     PINM1_dimension  = MAXANG + 1;
  const int     SM_dimension     = MAXANG + 1;
  const int     SP_dimension     = MAXANG + 1;

  const int     RBIGA_dimension = MAXTRM + 1;
  const int     CBIGA_dimension = MAXTRM + 1;
  const int     LITA_dimension  = MAXTRM + 1;
  const int     LITB_dimension  = MAXTRM + 1;

  const int     SMS_dimension = MXANG2 + 1;
  const int     SPS_dimension = MXANG2 + 1;

  bool    CALCMO [ CALCMO_dimension ];

  double  PIN [ PIN_dimension ], PINM1 [ PINM1_dimension ],
          RBIGA [ RBIGA_dimension ];

  dcomp   CBIGA [ CBIGA_dimension ], LITA [ LITA_dimension ],
          LITB [ LITB_dimension ], SM [ SM_dimension ], SMS [SMS_dimension ],
          SP [ SP_dimension ], SPS [ SPS_dimension ];

//      ..........................................................

                      // Save some input variables and replace them
                      // with values needed to do the self-test

      if ( PASS1_miev0_thr ) TESTMI( false, XX, CREFIN, MIMCUT,
             PERFCT, ANYANG, NMOM, IPOLZN, NUMANG,
             XMU, XMU_dimension, QEXT, QSCA,
             GQSC, SFORW, SBACK,
             S1, S1_dimension, S2, S2_dimension, 
             TFORW, TFORW_dimension, TBACK,
             TBACK_dimension, PMOM, PMOM_dimension,
	     MOMDIM,
	     &ANYSAV_thr,&CRESAV_thr,&IPOSAV_thr,&MIMSAV_thr,&NMOSAV_thr,&NUMSAV_thr,&PERSAV_thr,&XMUSAV_thr,&XXSAV_thr, // csz: Additional arguments for threading
             &MSGLIM_thr,&NUMMSG_errmsg_thr,&NUMMSG_wrtbad_thr); // csz: Additional arguments for threading
	
   label10:
                                         // Check input and calculate
                                         // certain variables from input

      CKINMI ( NUMANG, MAXANG, XX, PERFCT, CREFIN, MOMDIM, NMOM, IPOLZN,
	       ANYANG, XMU, XMU_dimension, CALCMO, CALCMO_dimension, NPQUAN, 
               &MSGLIM_thr,&NUMMSG_errmsg_thr,&NUMMSG_wrtbad_thr); // csz: Additional arguments for threading

      if ( (PERFCT) && (XX <= 0.1) )         // Use totally-reflecting
                                             // small-particle limit
         {
           SMALL1( XX, NUMANG, XMU, XMU_dimension, QEXT, QSCA, GQSC, SFORW,
             SBACK, S1, S1_dimension, S2, S2_dimension, TFORW, TFORW_dimension,
             TBACK, TBACK_dimension, LITA, LITA_dimension,
             LITB, LITB_dimension );

           NTRM = 2;
           goto label100;
         }


      NOABS = true;

      if ( ! PERFCT )
        {
         CIOR = CREFIN;

         if ( imag(CIOR) > 0.0 ) CIOR = conj( CIOR );

         MRE    = real( CIOR );
         MIM    = - imag( CIOR );
         NOABS  = (MIM <= MIMCUT);
         CIORIV = 1.0 / CIOR;
         RIORIV = 1.0 / MRE;

         if ( XX*MAX( 1.0, abs(CIOR) ) <= 0.1 )
           {                              // ** Use general-refractive-index
                                          // ** small-particle limit
             SMALL2( XX, CIOR, ( MIM > MIMCUT ), NUMANG,
                 XMU, XMU_dimension, QEXT, QSCA, GQSC, SFORW, SBACK,
                 S1, S1_dimension, S2, S2_dimension, TFORW, TFORW_dimension,
                 TBACK, TBACK_dimension, LITA, LITA_dimension,
                 LITB, LITB_dimension );
            NTRM = 2;
            goto label100;

           }
       }


      NANGD2 = ( NUMANG + 1 ) / 2;
      YESANG = ( NUMANG > 0 );

                                 //  ** Number of terms in Mie series; Eq R50
      if ( XX <= 8.0 )

         NTRM = static_cast <int> (XX + 4.*pow(XX, ONETHR) + 1.);

      else if ( XX < 4200. )
       {
         NTRM = static_cast <int> (XX + 4.05*pow(XX, ONETHR) + 2.);
       }
      else

         NTRM = static_cast <int> (XX + 4.*pow(XX, ONETHR) + 2.);


      if ( (NTRM+1) > MAXTRM )
                 ERRMSG((char *)"MIEV0--PARAMETER MaxTrm TOO SMALL", true,
		 &MSGLIM_thr,&NUMMSG_errmsg_thr); // csz: Additional arguments for threading

                        //   ** Calculate logarithmic derivatives of
                        //   ** J-Bessel-fcn., A-sub-(1 to NTrm)
      if ( ! PERFCT ) 
        {
          BIGA ( CIOR, XX, NTRM, NOABS, YESANG, RBIGA, RBIGA_dimension,
                 CBIGA, CBIGA_dimension,
	         &MSGLIM_thr,&NUMMSG_errmsg_thr); // csz: Additional arguments for threading
        }


                        //   ** Initialize Ricatti-Bessel functions
                        //   ** (psi,chi,zeta)-sub-(0,1) for upward
                        //   ** recurrence ( Eq. R19 )
      XINV   = 1.0 / XX;
      PSINM1 = sin( XX );
      CHINM1 = cos( XX );
      PSIN   = PSINM1*XINV - CHINM1;
      CHIN   = CHINM1*XINV + PSINM1;
      ZETNM1 = dcomp ( PSINM1, CHINM1 );
      ZETN   = dcomp ( PSIN, CHIN );
                                 //   ** Initialize previous coeffi-
                                 //   ** cients for GQSC series
      ANM1 = dcomp ( 0.0, 0.0 );
      BNM1 = dcomp ( 0.0, 0.0 );
                         //   ** Initialize angular function  pi
                         //   ** and sums for S+, S- ( Ref. 2, p. 1507 )
      if ( ANYANG )
        {
          for (int J = 1; J <= NUMANG; J++)
            {                           //  ** Eq. R39
              PINM1 [ J ] = 0.0;
              PIN [ J ] = 1.0;

              SP [ J ] = dcomp ( 0.0, 0.0 );
              SM [ J ] = dcomp ( 0.0, 0.0 );
            }
        }
      else
        {
          for (int J = 1; J <= NANGD2; J++)
            {                          ///  ** Eq. R39
              PINM1 [ J ] = 0.0;
              PIN [ J ] = 1.0;

              SP [ J ] = dcomp ( 0.0, 0.0 );
              SM [ J ] = dcomp ( 0.0, 0.0 );
              SPS [ J ] = dcomp ( 0.0, 0.0 );
              SMS [ J ] = dcomp ( 0.0, 0.0 );
            }
        }

                   //   ** Initialize Mie sums for efficiencies, etc.
      QSCA  = 0.0;
      GQSC  = 0.0;
      SFORW = dcomp ( 0., 0. );
      SBACK = dcomp ( 0., 0. );
      CSUM1 = dcomp ( 0., 0. );
      CSUM2 = dcomp ( 0., 0. );

// ---------  LOOP TO SUM MIE SERIES  -----------------------------------

      MM     = +1.0;
      SPIKE  = 1.0;

      for (int N = 1; N <= NTRM; N++)
        {           //           ** Compute various numerical coefficients
          FN     = static_cast <double> ( N );
          RN     = 1.0 / FN;
          NP1DN  = 1.0 + RN;
          TWONP1 = static_cast <double > (2*N + 1);
          COEFF  = TWONP1 / ( FN * static_cast <double> ( N + 1 ) );
          TCOEF  = TWONP1 * ( FN * static_cast <double> ( N + 1 ) );
  
                    //      ** Calculate Mie series coefficients
         if ( PERFCT )
           {        //            ** Totally-reflecting case; Eq R/A.1,2

              AN = ( ( FN*XINV )*PSIN - PSINM1 ) /
                   ( ( FN*XINV )*ZETN - ZETNM1 );
              BN = PSIN / ZETN;
           }
         else if ( NOABS )
           {        //                 ** No-absorption case; Eq (R16)

              CDENAN = ( RIORIV*RBIGA[N] + ( FN*XINV ) ) * ZETN - ZETNM1;
              AN   = ( ( RIORIV*RBIGA[N] + ( FN*XINV ) ) * PSIN - PSINM1 )
                     / CDENAN;
              CDENBN = ( MRE*RBIGA[N] + ( FN*XINV ) ) * ZETN - ZETNM1;
              BN   = ( ( MRE*RBIGA[N] + ( FN*XINV ) ) * PSIN - PSINM1 )
                     / CDENBN;
           }
         else
           {         //                 ** Absorptive case; Eq (R16)

              CDENAN = ( CIORIV*CBIGA[ N ] + ( FN*XINV ) )*ZETN - ZETNM1;
              CDENBN =   ( CIOR*CBIGA[ N ] + ( FN*XINV ) )*ZETN - ZETNM1;
              AN   = ( ( CIORIV*CBIGA[ N ] + ( FN*XINV ) )*PSIN - PSINM1 )
                       / CDENAN;
              BN     = ( ( CIOR*CBIGA[ N ] + ( FN*XINV ) )*PSIN - PSINM1 )
                       / CDENBN;
                      //                  ** Eq (R7)

              QSCA   = QSCA + TWONP1*( SQ( AN ) + SQ( BN ) );
           }

                      //** Save Mie coefficients for PMOM calculation

         LITA[ N ] = AN;
         LITB[ N ] = BN;


         if ( ( ! PERFCT ) && ( N > XX ) )
           {                                //  ** Flag resonance spikes
            DENAN  = abs( CDENAN );
            DENBN  = abs( CDENBN );
                                            //      ** Eq. R/B.9
            RATIO  = DENAN / DENBN;
            if ( ( RATIO <= 0.2 ) || ( RATIO >= 5.0 ) )
                { 
                  if (SPIKE > DENAN) SPIKE = DENAN;
                  if (SPIKE > DENBN) SPIKE = DENBN;
                }

           }
                               //  ** Increment Mie sums for non-angle-
                               //  ** dependent quantities

                               //                   ** Eq. R/B.2
         SFORW = SFORW + TWONP1*( AN + BN );
                               //                   ** Eq. R/B.5,6
         CSUM1 = CSUM1 + TCOEF *( AN - BN );
                               //                   ** Eq. R/B.1
         SBACK = SBACK + ( MM*TWONP1 )*( AN - BN );
                               //                   ** Eq. R/B.7,8
         CSUM2 = CSUM2 + ( MM*TCOEF ) *( AN + BN );

                               //         ** Eq (R8)

         GQSC  = GQSC  + ( FN - RN ) * real( ANM1 * conj ( AN ) +
                                             BNM1 * conj ( BN ) )
                 + COEFF * real( AN * conj ( BN ) );


         if ( YESANG )
           {                       //  ** Put Mie coefficients in form
                                   //  ** needed for computing S+, S-
                                   //  ** ( Eq R10 )
            ANP = COEFF*( AN + BN );
            BNP = COEFF*( AN - BN );

                                    // ** Increment Mie sums for S+, S-
                                    // ** while upward recursing
                                    // ** angular functions pi and tau
            if ( ANYANG )
              {                     //    ** Arbitrary angles

                                    //         ** vectorizable loop
               for (int J = 1; J <= NUMANG; J++)
                 {                  //            ** Eq. (R37b)

                  RTMP = ( XMU[J] * PIN[J] ) - PINM1[ J ];

                                    //            ** Eq. (R38b)
                  TAUN   = FN * RTMP - PINM1[ J ];

                                    //              ** Eq (R10)

                  SP[ J ] = SP[ J ] + ANP * ( PIN[ J ] + TAUN );
                  SM[ J ] = SM[ J ] + BNP * ( PIN[ J ] - TAUN );

                  PINM1[ J ] = PIN[ J ];
                                    //            ** Eq. R37c

                  PIN[ J ] = ( XMU[ J ] * PIN[ J ] ) + NP1DN * RTMP;
                 }
              }
            else
              {
                                // ** Angles symmetric about 90 degrees
               ANPM = MM*ANP;
               BNPM = MM*BNP;
                                //         ** vectorizable loop
               for (int J = 1; J <= NANGD2; J++)
                 {              //                ** Eq. (R37b)

                  RTMP = ( XMU[J] * PIN[J] ) - PINM1[ J ];

                                //                ** Eq. (R38b)
                  TAUN = FN * RTMP - PINM1[ J ];

                                //                ** Eq (R10,12)

                  SP [ J ] = SP [ J ] + ANP * ( PIN[ J ] + TAUN );
                  SMS[ J ] = SMS[ J ] + BNPM *( PIN[ J ] + TAUN );
                  SM [ J ] = SM [ J ] + BNP * ( PIN[ J ] - TAUN );
                  SPS[ J ] = SPS[ J ] + ANPM *( PIN[ J ] - TAUN );

                  PINM1[ J ] = PIN[ J ];
                               //                 ** Eq. R37c

                  PIN[ J ] = ( XMU[J] * PIN[J] ) + NP1DN * RTMP;
                 }

              }

           }
                        // ** Update relevant quantities for next
                        // ** pass through loop
         MM   = - MM;
         ANM1 = AN;
         BNM1 = BN;
                        //  ** Upward recurrence for Ricatti-Bessel
                        //  ** functions ( Eq. R17 )

         ZET    = dcomp(( TWONP1*XINV ), 0.) * ZETN - ZETNM1;
         ZETNM1 = ZETN;
         ZETN   = ZET;
         PSINM1 = PSIN;
         PSIN   = real( ZETN );
        }

//---------- END LOOP TO SUM MIE SERIES --------------------------------


                                     //   ** Eq (R6)
      QEXT = 2. / pow ( XX, 2.) * real( SFORW );

      if ( PERFCT || NOABS )

         QSCA = QEXT;

      else

         QSCA = 2./ pow ( XX, 2.) * QSCA;

      GQSC   = 4./ pow ( XX, 2.) * GQSC;
      SFORW  = 0.5*SFORW;
      SBACK  = 0.5*SBACK;
      TFORW[ 1 ] =  0.5*SFORW - 0.125*CSUM1;
      TFORW[ 2 ] =  0.5*SFORW + 0.125*CSUM1;
      TBACK[ 1 ] = -0.5*SBACK + 0.125*CSUM2;
      TBACK[ 2 ] =  0.5*SBACK + 0.125*CSUM2;


      if ( YESANG ) 
         {                   //  ** Recover scattering amplitudes
                             //  ** from S+, S- ( Eq (R11) )

         if ( ANYANG )
           {                 //           ** vectorizable loop
            for (int J = 1; J <= NUMANG; J++)
              {              //                     ** Eq (R11)
               S1[ J ] = 0.5*( SP[ J ] + SM[ J ] );
               S2[ J ] = 0.5*( SP[ J ] - SM[ J ] );
              }
           }
         else
           {                 //           ** vectorizable loop
            for (int J = 1; J <= NANGD2 ; J++)
              {              //                     ** Eq (R11)
               S1[ J ] = 0.5*( SP[ J ] + SM[ J ] );
               S2[ J ] = 0.5*( SP[ J ] - SM[ J ] );
              }
                             //           ** vectorizable loop
            for (int J = 1; J <= NANGD2; J++)
              {
               S1[ NUMANG + 1 - J ] = 0.5*( SPS[ J ] + SMS[ J ] );
               S2[ NUMANG + 1 - J ] = 0.5*( SPS[ J ] - SMS[ J ] );
              }
           }

         }
                             //   ** Calculate Legendre moments

  label100:


      if ( NMOM > 0 ) LPCOEF( NTRM, NMOM, IPOLZN, MOMDIM,
                              CALCMO, CALCMO_dimension, NPQUAN,
                              LITA, LITA_dimension, LITB, LITB_dimension,
                              PMOM, PMOM_dimension,
			      &PASS1_lpcoef_thr,RECIP_thr, // csz: Additional arguments for threading
			      &MSGLIM_thr,&NUMMSG_errmsg_thr); // csz: Additional arguments for threading

      if ( imag ( CREFIN ) > 0.0 )
        {                             //  ** Take complex conjugates
                                      //  ** of scattering amplitudes

         SFORW = conj( SFORW );
         SBACK = conj( SBACK );

         for (int I = 1; I <= 2; I++)
           {
            TFORW[ I ] = conj( TFORW[ I ] );
            TBACK[ I ] = conj( TBACK[ I ] );
           }

         for (int J = 1; J <= NUMANG; J++)
           {
            S1[ J ] = conj( S1[ J ] );
            S2[ J ] = conj( S2[ J ] );
           }

        }


        if ( PASS1_miev0_thr )
          {            //   ** Compare test case results with
                       //   ** correct answers and abort if bad;
                       //   ** otherwise restore user input and proceed

            TESTMI( true, XX, CREFIN, MIMCUT,
                   PERFCT, ANYANG, NMOM, IPOLZN, NUMANG,
                   XMU, XMU_dimension, QEXT, QSCA,
                   GQSC, SFORW, SBACK,
                   S1, S1_dimension, S2, S2_dimension, 
                   TFORW, TFORW_dimension, TBACK,
                   TBACK_dimension, PMOM, PMOM_dimension,
		   MOMDIM,
		   &ANYSAV_thr,&CRESAV_thr,&IPOSAV_thr,&MIMSAV_thr,&NMOSAV_thr,&NUMSAV_thr,&PERSAV_thr,&XMUSAV_thr,&XXSAV_thr, // csz: Additional arguments for threading
	           &MSGLIM_thr,&NUMMSG_errmsg_thr,&NUMMSG_wrtbad_thr); // csz: Additional arguments for threading
  
         PASS1_miev0_thr  = false;
         goto label10;

          }


//    IF( PRNT( 1 ) .OR. PRNT( 2 ) ) 
//   &  CALL MIPRNT( PRNT, XX, PERFCT, CREFIN, NUMANG, XMU, QEXT,
//   &               QSCA, GQSC, NMOM, IPOLZN, MOMDIM, CALCMO, PMOM,
//   &               SFORW, SBACK, TFORW, TBACK, S1, S2 )

  return rcd;
} /* end MIEV0() */

int MIEV0_drv() {
/* -----------------------------------------------------------------
  Note about the original FORTRAN code:

  NOTE:  R1MACH is invoked through a USE RDI1MACH_f90 statement,
         which is Fortran-90.  The reasons are given in the file
         README.RDI1MACH.  One can do things the old way by 
         removing the USE statement and getting R1MACH from netlib.
         Remember that MODULEs must precede the program unit USEing
         them in Fortran-90, or be separately compiled and linked.
  ----------------------------------------------------------------- */

/*        Run 19 test cases for MIEV0, testing all the pathways
          through the code.  In parentheses beneath each result is
          printed its ratio to the answer obtained on a Cray
          computer using 14-digit precision.  If these ratios drift
          away from unity as size parameter increases, your computer
          probably is of lower precision than the Cray, OR it may
          handle arithmetic differently (truncating instead of
          rounding, for example).  Before becoming overly concerned
          about non-unit ratios, re-run in double precision
          (most compilers have an auto-doubling option) and see if
          your results improve.

      NOTES:

         ** Set NoPMOM = True at the beginning of the executable
            statements below if using NoPMOM version of MIEV0

         ** Temporarily set PARAMETER MAXTRM in MIEV0 to 10,100 to
            run these test problems.  (Be sure to lower it again to
            some reasonable value when you are finished.)

         ** Timing is done by calls to the Fortran-90 intrinsic
            function SYSTEM_CLOCK; if you absolutely can't find an
            f90 compiler, these calls can just be deleted, or replaced
            using some local system-dependent routine.

         ** To keep storage requirements reasonable (i.e. to keep from
            having to set  MAXTRM = 10100 in LPCOEF also), and to avoid
            enormous DATA statements to hold the correct answers, 
            Legendre moments are not checked for the size parameter = 
            10,000 cases.  Also, only the first two moments are checked 
            for size parameter = 1000 cases.

         ** Pay especial attention to the cases where two successive
            values of size parameter are small and close.  The lower
            value uses small-particle approximations, the upper value
            uses the full Mie series.  The results should be close.

         ** TBACK is the quantity most sensitive to precision.
            This is because the coefficients in the TBACK series
            are large, increasing, and alternating in sign,
            making the order in which the terms are summed
            important in finite-precision arithmetic.  SBACK
            and S1, S2 near 180 deg are also sensitive.

         ** High-order Legendre moments, being near zero, are also
            subject to a lot of cancellation in their defining series,
            and hence may also be sensitive to computer precision.

      CALLING TREE:

           MITEST
              RATIO
                 R1MACH
              MIEV0
           CHEKMI
*/

  /*--------------------------------------------------------------------
  -----------  I/O SPECIFICATIONS FOR SUBROUTINE  MIEV0  ---------------
  ----------------------------------------------------------------------*/

  const int MAXANG = 7;
  const int MOMDIM = 200;

  const int XMU_dimension  = MAXANG + 1;
  const int PRNT_dimension = 2 + 1;
  const int PMOM_dimension = MOMDIM + 1;
  const int S1_dimension   = MAXANG + 1;
  const int S2_dimension   = MAXANG + 1;

  const int TFORW_dimension = 2 + 1;
  const int TBACK_dimension = 2 + 1;

  bool ANYANG, PERFCT;

  bool PRNT[3];  // I intend not to use prnt[0], and leave instructions alone
               // for prnt[1], and prnt[2].

  int IPOLZN, NUMANG, NMOM;

  double GQSC, MIMCUT, PMOM [ PMOM_dimension ][ PMOM_second_dimension ], QEXT, QSCA, SPIKE,
         XMU [ XMU_dimension ], XX;

  dcomp CREFIN, SFORW, SBACK, S1 [ S1_dimension ], S2 [ S2_dimension ],
        TFORW [ TFORW_dimension ], TBACK [ TBACK_dimension ];

/*----------------------------------------------------------------------

  --------------- LOCAL VARIABLES --------------------------------------

      .. Local Scalars ..     */

  bool NOPMOM;
  int I, J, K, NCAS, NPQUAN, time0, time1, cntrat, maxcnt;
  double DEGPOL, FNORM, I1, I2, INTEN, PI, QABS, TESTIN;

/*    ..
      .. Local Arrays .. */

  double ANGLE [ MAXANG + 1];

/*----------------------------------------------------------------------
         Input specifications and 'correct' answers to test problems
  ----------------------------------------------------------------------*/

  const int NCASES = 19 + 1;
  bool     TESTAN [ NCASES ];
  int      TESTIP [ NCASES ];
  double   TESTXX [ NCASES ];
  dcomp    TESTCR [ NCASES ];

      PI     = 2.* asin( 1.0 );
      NOPMOM = false;
                             // Set MIEV0 input variables that are
                             // the same for every test case
      MIMCUT = 1.e-6;
      NUMANG = MAXANG;

      for (int I = 1; I <= NUMANG; I++){
         ANGLE [ I ] = ( I - 1 )*180. / ( NUMANG - 1 );
         XMU [ I ] = cos( PI / 180.*ANGLE [ I ] );
      }

                /* Call once for very small perfectly conducting sphere;
                ** this does MIEV0 self-test, so that self-test does
                ** not affect timing results, and tests print option */

      XX     = 0.001;
      CREFIN = dcomp( 0., 0. );
      PERFCT = true;
      ANYANG = true;
      IPOLZN = +1234;
      NMOM   = 1;

      if ( NOPMOM ) NMOM = 0;

      PRNT [ 1 ] = true;
      PRNT [ 2 ] = true;

      for (int I = 1; I <= PMOM_dimension; I++)
        for (int J = 1; J <= 5; J++)
          PMOM [I][J] = 0.;

      /* csz++
	 Rename main solution routine for compatibility with BoH83
	 MIEV0( XX, CREFIN, PERFCT, MIMCUT, */
      mie_sph_Wis79( XX, CREFIN, PERFCT, MIMCUT,
        ANYANG, NUMANG, XMU, XMU_dimension,
        NMOM, IPOLZN, MOMDIM, PRNT, PRNT_dimension,
        QEXT, QSCA, GQSC,
        PMOM, PMOM_dimension, SFORW,
        SBACK, S1, S1_dimension, S2, S2_dimension,
        TFORW, TFORW_dimension, TBACK, TBACK_dimension,
        SPIKE );

/*    cout << endl << " **** C output" << endl;
      cout << "XX    = " << XX <<endl;
      cout << "QEXT  = " << QEXT <<endl;
      cout << "QSCA  = " << QSCA  <<endl;
      cout << "GQSC  = " << GQSC  <<endl;
      cout << "SPIKE = " << SPIKE  <<endl;
**/
      cout.setf(ios::scientific, ios::floatfield);
      cout << setprecision(20) << XX <<endl;
      cout << QEXT <<endl;
      cout << QSCA  <<endl;
      cout << GQSC  <<endl;
      cout << SPIKE  <<endl;

      PRNT [ 1 ] = false;
      PRNT [ 2 ] = false;

  // Initialize for Test cases
  
  // TESTXX [] :
  
  TESTXX [  1 ] =     0.09900;
  TESTXX [  2 ] =     0.10100;
  TESTXX [  3 ] =   100.00000;
  TESTXX [  4 ] = 10000.00000;
  TESTXX [  5 ] =     0.09900;
  TESTXX [  6 ] =     0.10100;
  TESTXX [  7 ] =    10.00000;
  TESTXX [  8 ] =  1000.00000;
  TESTXX [  9 ] =     1.00000;
  TESTXX [ 10 ] =   100.00000;
  TESTXX [ 11 ] = 10000.00000;
  TESTXX [ 12 ] =     0.05500;
  TESTXX [ 13 ] =     0.05600;
  TESTXX [ 14 ] =     1.00000;
  TESTXX [ 15 ] =   100.00000;
  TESTXX [ 16 ] = 10000.00000;
  TESTXX [ 17 ] =     1.00000;
  TESTXX [ 18 ] =   100.00000;
  TESTXX [ 19 ] = 10000.00000;
  
  // TESTCR [] :
  
  TESTCR [  1 ] = dcomp ( 0.0000E+00 ,  0.0000E+00);
  TESTCR [  2 ] = dcomp ( 0.0000E+00 ,  0.0000E+00);
  TESTCR [  3 ] = dcomp ( 0.0000E+00 ,  0.0000E+00);
  TESTCR [  4 ] = dcomp ( 0.0000E+00 ,  0.0000E+00);
  TESTCR [  5 ] = dcomp ( 0.7500E+00 ,  0.0000E+00);
  TESTCR [  6 ] = dcomp ( 0.7500E+00 ,  0.0000E+00);
  TESTCR [  7 ] = dcomp ( 0.7500E+00 ,  0.0000E+00);
  TESTCR [  8 ] = dcomp ( 0.7500E+00 ,  0.0000E+00);
  TESTCR [  9 ] = dcomp ( 0.1330E+01 , -0.1000E-04);
  TESTCR [ 10 ] = dcomp ( 0.1330E+01 , -0.1000E-04);
  TESTCR [ 11 ] = dcomp ( 0.1330E+01 , -0.1000E-04);
  TESTCR [ 12 ] = dcomp ( 0.1500E+01 , -0.1000E+01);
  TESTCR [ 13 ] = dcomp ( 0.1500E+01 , -0.1000E+01);
  TESTCR [ 14 ] = dcomp ( 0.1500E+01 , -0.1000E+01);
  TESTCR [ 15 ] = dcomp ( 0.1500E+01 , -0.1000E+01);
  TESTCR [ 16 ] = dcomp ( 0.1500E+01 , -0.1000E+01);
  TESTCR [ 17 ] = dcomp ( 0.1000E+02 , -0.1000E+02);
  TESTCR [ 18 ] = dcomp ( 0.1000E+02 , -0.1000E+02);
  TESTCR [ 19 ] = dcomp ( 0.1000E+02 , -0.1000E+02);
  
  // TESTAN [] :
  
  TESTAN [  1 ] = 0;
  TESTAN [  2 ] = 0;
  TESTAN [  3 ] = 0;
  TESTAN [  4 ] = 0;
  TESTAN [  5 ] = 1;
  TESTAN [  6 ] = 1;
  TESTAN [  7 ] = 1;
  TESTAN [  8 ] = 1;
  TESTAN [  9 ] = 1;
  TESTAN [ 10 ] = 1;
  TESTAN [ 11 ] = 1;
  TESTAN [ 12 ] = 0;
  TESTAN [ 13 ] = 0;
  TESTAN [ 14 ] = 0;
  TESTAN [ 15 ] = 0;
  TESTAN [ 16 ] = 0;
  TESTAN [ 17 ] = 1;
  TESTAN [ 18 ] = 1;
  TESTAN [ 19 ] = 1;
  
  // TESTIP [] :
  
  TESTIP [  1 ] = -1234;
  TESTIP [  2 ] = -1234;
  TESTIP [  3 ] = -1234;
  TESTIP [  4 ] = -1234;
  TESTIP [  5 ] =     0;
  TESTIP [  6 ] =     0;
  TESTIP [  7 ] =     0;
  TESTIP [  8 ] =     0;
  TESTIP [  9 ] =  1234;
  TESTIP [ 10 ] =  1234;
  TESTIP [ 11 ] =  1234;
  TESTIP [ 12 ] = -1234;
  TESTIP [ 13 ] = -1234;
  TESTIP [ 14 ] = -1234;
  TESTIP [ 15 ] = -1234;
  TESTIP [ 16 ] = -1234;
  TESTIP [ 17 ] =     0;
  TESTIP [ 18 ] =     0;
  TESTIP [ 19 ] =     0;

//                                ** Loop over test cases

      for (int NCAS = 1; NCAS <= NCASES - 1; NCAS++)
        {                                 // TTT
         XX     = TESTXX [ NCAS ];
         CREFIN = TESTCR [ NCAS ];

         if ( NCAS <= 4 ) PERFCT = true;
         if ( NCAS >  4 ) PERFCT = false;

         ANYANG = TESTAN [ NCAS ];
         IPOLZN = TESTIP [ NCAS ];
         NMOM   = static_cast <int> (2.* XX);

         if ( (XX < 1.) || (XX > 100.) ) NMOM   = 1;
         if ( XX > 1000. )               NMOM   = 0;
         if ( NOPMOM )                   NMOM   = 0;

         mie_sph_Wis79( XX, CREFIN, PERFCT, MIMCUT,
           ANYANG, NUMANG, XMU, XMU_dimension,
           NMOM, IPOLZN, MOMDIM, PRNT, PRNT_dimension,
           QEXT, QSCA, GQSC,
           PMOM, PMOM_dimension, SFORW,
           SBACK, S1, S1_dimension, S2, S2_dimension,
           TFORW, TFORW_dimension, TBACK, TBACK_dimension,
           SPIKE );

/*    cout << " " << endl;
      cout << " Test case : " << NCAS  << endl;
      cout << "XX    = " << XX <<endl;
      cout << "QEXT  = " << QEXT <<endl;
      cout << "QSCA  = " << QSCA  <<endl;
      cout << "GQSC  = " << GQSC  <<endl;
      cout << "SPIKE = " << SPIKE  <<endl;
*/
      cout << XX <<endl;
      cout << QEXT <<endl;
      cout << QSCA  <<endl;
      cout << GQSC  <<endl;
      cout << SPIKE  <<endl;

        }   // TTT

  return 0;
} /* end MIEV0_drv() */

void PARSEI ( int IPOLZN, int N[], int N_dimension )
{
/* The original FORTRAN code has a variable POLZN whose digits have 
   individual meaning -- but the value of POLZN itself does not. Thus,
   POLZN has to be 'taken apart'. The original code accomplished this
   differently from the way done here. The task does not have to be 
   accomplished as it is done here, but we wanted to output exactly
   the same as Wiscombe's program.
*/
  int M [ 4 ];

  M [ 2 ] = abs(IPOLZN) % 1000;
  M [ 1 ] =     M [ 2 ] % 100 ;
  M [ 0 ] =     M [ 1 ] % 10  ;
  M [ 3 ] = abs(IPOLZN) / 1000; if (abs (IPOLZN) < 1000) M [ 3 ] = -16;
  M [ 2 ] =     M [ 2 ] / 100 ; if (abs (IPOLZN) < 100 ) M [ 2 ] = -16;
  M [ 1 ] =     M [ 1 ] / 10  ; if (abs (IPOLZN) < 10  ) M [ 1 ] = -16;

  int index = 0;
  for (int k = 3; k >= 0; k--)
    { N [ index ] = M [ k ];
      index++;
    }
  

  return;
} /* end PARSEI() */

void SMALL1( double XX, int NUMANG, double XMU[], int XMU_dimension,
             double &QEXT, double &QSCA, double &GQSC, dcomp &SFORW,
             dcomp &SBACK, dcomp S1[], int S1_dimension, dcomp S2[],
             int S2_dimension, dcomp TFORW[], int TFORW_dimension,
             dcomp TBACK[], int TBACK_dimension,
             dcomp A[], int LITA_dimension,
             dcomp B[], int LITB_dimension )
{
/*      Small-particle limit of Mie quantities in totally reflecting
        limit ( Mie series truncated after 2 terms )

         A,B       First two Mie coefficients, with numerator and
                   denominator expanded in powers of XX ( a factor
                   of XX**3 is missing but is restored before return
                   to calling program )  ( Ref. 2, p. 1508 )
*/

//    .. Parameters ..

      const double  TWOTHR = 2./3.;
      const double  FIVTHR = 5./3.;
      const double  FIVNIN = 5./9.;


//    .. Local Scalars ..

      int     J;
      double  RTMP;
      dcomp   CTMP;

//                                                      ** Eq. R/A.5
      A[ 1 ] = dcomp( 0., TWOTHR * ( 1. - 0.2 * pow (  XX, 2.) ) ) /
               dcomp( 1. - 0.5 * pow( XX, 2.), TWOTHR * pow ( XX ,3.) );
//                                                      ** Eq. R/A.6
      B[ 1 ] = dcomp( 0., - ( 1. - 0.1 * pow ( XX, 2.) ) / 3.) /
               dcomp( 1. + 0.5 * pow ( XX, 2.), - pow ( XX, 3.) / 3.);
//                                                      ** Eq. R/A.7,8
      A[ 2 ] = dcomp( 0.,   pow ( XX, 2.) / 30.);
      B[ 2 ] = dcomp( 0., - pow ( XX, 2.) / 45.);
//                                                      ** Eq. R/A.9
      QSCA = 6.* pow ( XX, 4.) * ( SQ( A[1] ) + SQ( B[1] ) +
                 FIVTHR*( SQ( A[2] ) + SQ( B[2] ) ) );
      QEXT = QSCA;
//                                                      ** Eq. R/A.10
      GQSC = 6.* pow (XX, 4.) * real( A[1] * conj( A[2] + B[1] ) +
               ( B[1] + FIVNIN*A[2] ) * conj( B[2] ) );

      RTMP   = 1.5 * pow ( XX, 3.);
      SFORW  = RTMP*( A[1] + B[1] + FIVTHR*( A[2] + B[2] ) );
      SBACK  = RTMP*( A[1] - B[1] - FIVTHR*( A[2] - B[2] ) );
      TFORW[ 1 ] = RTMP*( B[1] + FIVTHR*( 2.*B[2] - A[2] ) );
      TFORW[ 2 ] = RTMP*( A[1] + FIVTHR*( 2.*A[2] - B[2] ) );
      TBACK[ 1 ] = RTMP*( B[1] - FIVTHR*( 2.*B[2] + A[2] ) );
      TBACK[ 2 ] = RTMP*( A[1] - FIVTHR*( 2.*A[2] + B[2] ) );


      for ( int J = 1; J<=NUMANG; J++)
        {                                         // ** Eq. R/A.11,12
         S1[ J ] = RTMP*( A[1] + B[1]*XMU[ J ] +
                          FIVTHR*( A[2]*XMU[ J ] + 
                                   B[2]*( 2.*pow (XMU[ J ], 2.) - 1.) ) );
         S2[ J ] = RTMP*( B[1] + A[1]*XMU[ J ] +
                          FIVTHR*( B[2]*XMU[ J ] + 
                                   A[2]*( 2.*pow (XMU[ J ], 2.) - 1.) ) );
        }

                                   // ** Recover actual Mie coefficients
      A[ 1 ] = A[1] * pow (XX, 3.);
      A[ 2 ] = A[2] * pow (XX, 3.);
      B[ 1 ] = B[1] * pow (XX, 3.);
      B[ 2 ] = B[2] * pow (XX, 3.);

  return;
} /* end SMALL1() */

void SMALL2( double XX, dcomp CIOR, bool CALCQE, int NUMANG,
             double XMU[], int XMU_dimension,
             double &QEXT, double &QSCA, double &GQSC, dcomp &SFORW,
             dcomp &SBACK, dcomp S1[], int S1_dimension, dcomp S2[],
             int S2_dimension, dcomp TFORW[], int TFORW_dimension,
             dcomp TBACK[], int TBACK_dimension,
             dcomp A[], int LITA_dimension,
             dcomp B[], int LITB_dimension )
{

/*      Small-particle limit of Mie quantities for general refractive
        index ( Mie series truncated after 2 terms )

         A,B       First two Mie coefficients, with numerator and
                   denominator expanded in powers of XX ( a factor
                   of XX**3 is missing but is restored before return
                   to calling program )

         CIORSQ    Square of refractive index
*/

//    .. Parameters ..

  const double  TWOTHR = 2./3., FIVTHR = 5./3.;

//    .. Local Scalars ..

  int     J;
  double  RTMP;
  dcomp   CIORSQ, CTMP;


  CIORSQ = pow (CIOR, 2.);
  CTMP   = dcomp ( 0., TWOTHR )*( CIORSQ - 1.0 );

//                                          ** Eq. R42a
  A[ 1 ] = CTMP*( 1.- 0.1*pow(XX, 2.) +
              ( CIORSQ / 350. + 1./280.)*pow(XX, 4.) ) /
              ( CIORSQ + 2.+ ( 1.- 0.7*CIORSQ )*pow(XX, 2.) -
              ( pow(CIORSQ, 2.) / 175.- 0.275*CIORSQ + 0.25 )*pow(XX, 4.) +
              pow(XX, 3.) * CTMP * ( 1.- 0.1*pow(XX, 2.) ) );

//                                          ** Eq. R42b
  B[ 1 ] = ( pow(XX, 2.) / 30. )*CTMP*( 1.+
              ( CIORSQ / 35. - 1./ 14.)*pow(XX, 2.) ) /
              ( 1.- ( CIORSQ / 15. - 1./6.)*pow(XX, 2.) );

//                                          ** Eq. R42c

  A[ 2 ] = ( 0.1*pow(XX, 2.) )*CTMP*( 1.- pow(XX, 2.) / 14. ) /
              ( 2.*CIORSQ + 3.- ( CIORSQ / 7.- 0.5 ) * pow(XX, 2.) );

//                                          ** Eq. R40a

  QSCA = (6.*pow(XX, 4.)) * ( SQ( A[1] ) + SQ( B[1] ) +
              FIVTHR * SQ( A[2] ) );

//                                          ** Eq. R40b
  QEXT = QSCA;
  if ( CALCQE ) QEXT = 6.*XX * real( A[1] + B[1] + FIVTHR*A[2] );

//                                          ** Eq. R40c

  GQSC = (6.*pow(XX, 4.)) * real( A[1]*conj( A[2] + B[1] ) );

  RTMP   = 1.5 * pow(XX, 3.);
  SFORW  = RTMP*( A[1] + B[1] + FIVTHR*A[2] );
  SBACK  = RTMP*( A[1] - B[1] - FIVTHR*A[2] );
  TFORW[ 1 ] = RTMP*( B[1] - FIVTHR*A[2] );
  TFORW[ 2 ] = RTMP*( A[1] + 2.*FIVTHR*A[2] );
  TBACK[ 1 ] = TFORW[1];
  TBACK[ 2 ] = RTMP*( A[1] - 2.*FIVTHR*A[2] );


  for (J = 1; J <= NUMANG; J++)
    {                                       //  ** Eq. R40d,e
         S1[ J ] = RTMP*( A[1] + ( B[1] + FIVTHR*A[2] )*XMU[ J ] );
         S2[ J ] = RTMP*( B[1] + A[1]*XMU[ J ] +
                          FIVTHR*A[2]*( 2.*pow(XMU[ J ], 2.) - 1.) );
    }

//                                    ** Recover actual Mie coefficients
  A[ 1 ] = pow(XX, 3.) * A[1];
  A[ 2 ] = pow(XX, 3.) * A[2];
  B[ 1 ] = pow(XX, 3.) * B[1];
  B[ 2 ] = dcomp( 0., 0.);

  return;
} /* end SMALL2() */

double SQ ( dcomp CTMP )
{
  return real ( (CTMP)*conj(CTMP) );
} /* end SQ() */

void TESTMI( bool COMPAR, double &XX, dcomp &CREFIN, double &MIMCUT,
             bool &PERFCT, bool &ANYANG, int &NMOM, int &IPOLZN, int &NUMANG,
             double XMU[], int XMU_dimension, double &QEXT, double &QSCA,
             double &GQSC, dcomp &SFORW, dcomp &SBACK,
             dcomp S1[], int S1_dimension, dcomp S2[], int S2_dimension, 
             dcomp TFORW[], int TFORW_dimension, dcomp TBACK[],
             int TBACK_dimension, double PMOM[][PMOM_second_dimension], int PMOM_dimension,
             int MOMDIM,
	     bool *ANYSAV_thr,dcomp *CRESAV_thr,int *IPOSAV_thr,double *MIMSAV_thr,int *NMOSAV_thr,int *NUMSAV_thr,bool *PERSAV_thr,double *XMUSAV_thr,double *XXSAV_thr, // csz: Additional arguments for threading
             bool *MSGLIM_thr,int *NUMMSG_errmsg_thr,int *NUMMSG_wrtbad_thr) // csz: Additional arguments for threading
{
/*        Set up to run test case when  COMPAR = False;  when  = True,
          compare Mie code test case results with correct answers
          and abort if even one result is inaccurate.

          The test case is :  Mie size parameter = 10
                              refractive index   = 1.5 - 0.1 i
                              scattering angle = 140 degrees
                              1 Sekera moment

          Results for this case may be found among the test cases
          at the end of reference (1).

          *** NOTE *** When running on some computers, esp. in single
          precision, the Accur criterion below may have to be relaxed.
          However, if Accur must be set larger than 10**-3 for some
          size parameters, your computer is probably not accurate
          enough to do Mie computations for those size parameters.
*/


//    .. Local Scalars ..

      bool   OK;
      int    M, N;
      double ACCUR, CALC, EXACT, TESTGQ, TESTQE, TESTQS;
      //static bool PERSAV, ANYSAV; // csz: Replace static value with local argument
      //static double XXSAV, MIMSAV, XMUSAV; // csz: Replace static value with local argument
      //static dcomp CRESAV; // csz: Replace static value with local argument
      //static int NMOSAV, IPOSAV, NUMSAV; // csz: Replace static value with local argument
      dcomp  TESTS1, TESTS2, TESTSB, TESTSF;

//    .. Local Arrays ..

      bool   CALCMO [ 4 + 1 ], PRNT [ 2 + 1];
      double TESTPM [ 2 ];
      dcomp  TESTTB [ 2 + 1 ], TESTTF [ 2 + 1 ];

// Parameters

      TESTQE       = 2.459791;
      TESTQS       = 1.235144;
      TESTGQ       = 1.139235;
      TESTSF       = dcomp ( 61.49476, -3.177994 );
      TESTSB       = dcomp ( 1.493434, 0.2963657 );
      TESTS1       = dcomp ( -0.1548380, -1.128972 );
      TESTS2       = dcomp ( 0.05669755, 0.5425681 );
      TESTTF [ 1 ] = dcomp ( 12.95238, -136.6436 );
      TESTTF [ 2 ] = dcomp ( 48.54238, 133.4656 );
      TESTTB [ 1 ] = dcomp ( 41.88414, -15.57833 );
      TESTTB [ 2 ] = dcomp ( 43.37758, -15.28196 );
      TESTPM [ 0 ] = 227.1975;
      TESTPM [ 1 ] = 183.6898;

      ACCUR  = 1.e-4;

      if ( !COMPAR )
        {                        // ** Save certain user input values
         *XXSAV_thr  = XX;
         *CRESAV_thr = CREFIN;
         *MIMSAV_thr = MIMCUT;
         *PERSAV_thr = PERFCT;
         *ANYSAV_thr = ANYANG;
         *NMOSAV_thr = NMOM;
         *IPOSAV_thr = IPOLZN;
         *NUMSAV_thr = NUMANG;
         *XMUSAV_thr = XMU[ 1 ];
                                 // ** Reset input values for test case
         XX     = 10.0;
         CREFIN = dcomp ( 1.5, -0.1 );
         MIMCUT = 0.0;
         PERFCT = false;
         ANYANG = true;
         NMOM   = 1;
         IPOLZN = -1;
         NUMANG = 1;
         XMU[ 1 ] = -0.7660444;
        }
      else
        {                        //  ** Compare test case results with
                                 //  ** correct answers and abort if bad
         OK = true;

         if ( WRONG( QEXT,TESTQE ) )
             OK = TSTBAD( (char *)"QEXT", abs( ( QEXT - TESTQE ) / TESTQE ) );

         if ( WRONG( QSCA,TESTQS ) )
             OK = TSTBAD( (char *)"QSCA", abs( ( QSCA - TESTQS ) / TESTQS ) );

         if ( WRONG( GQSC,TESTGQ ) )
             OK = TSTBAD( (char *)"GQSC", abs( ( GQSC - TESTGQ ) / TESTGQ ) );

         if ( WRONG( real( SFORW ),real( TESTSF ) ) ||
              WRONG( imag( SFORW ),imag( TESTSF ) ) )
             OK = TSTBAD( (char *)"SFORW", abs( ( SFORW - TESTSF ) / TESTSF ) );

         if ( WRONG( real( SBACK ),real( TESTSB ) ) ||
              WRONG( imag( SBACK ),imag( TESTSB ) ) )
             OK = TSTBAD( (char *)"SBACK", abs( ( SBACK - TESTSB ) / TESTSB ) );

         if ( WRONG( real( S1[1] ),real( TESTS1 ) ) ||
              WRONG( imag( S1[1] ),imag( TESTS1 ) ) )
             OK = TSTBAD( (char *)"S1", abs( ( S1[1] - TESTS1 ) / TESTS1 ) );

         if ( WRONG( real( S2[1] ),real( TESTS2 ) ) ||
              WRONG( imag( S2[1] ),imag( TESTS2 ) ) )
             OK = TSTBAD( (char *)"S2", abs( ( S2[1] - TESTS2 ) / TESTS2 ) );


         for ( N = 1; N <= 2; N++)
           {
            if ( WRONG( real( TFORW[N] ),real( TESTTF[N] ) ) ||
                 WRONG( imag( TFORW[N] ),
                 imag( TESTTF[N] ) ) ) OK = TSTBAD( (char *)"TFORW",
                 abs( ( TFORW[N] - TESTTF[N] ) / TESTTF[N] ) );

            if ( WRONG( real( TBACK[N] ),real( TESTTB[N] ) ) ||
                 WRONG( imag( TBACK[N] ),
                 imag( TESTTB[N] ) ) ) OK = TSTBAD( (char *)"TBACK",
                 abs( ( TBACK[N] - TESTTB[N] ) / TESTTB[N] ) );
           }


         for ( M = 0; M <= 1; M++)
            if ( WRONG( PMOM[M][1], TESTPM[M] ) )
                 OK =  TSTBAD( (char *)"PMOM", abs( (PMOM[M][1]-TESTPM[M]) /
                                            TESTPM[M] ) );


         if ( !OK )
           {
            PRNT[ 1 ]   = true;
            PRNT[ 2 ]   = true;
            CALCMO[ 1 ] = true;
            CALCMO[ 2 ] = false;
            CALCMO[ 3 ] = false;
            CALCMO[ 4 ] = false;

/*          CALL MIPRNT( PRNT, XX, PERFCT, CREFIN, NUMANG, XMU, QEXT,
     &                   QSCA, GQSC, NMOM, IPOLZN, MOMDIM, CALCMO, PMOM,
     &                   SFORW, SBACK, TFORW, TBACK, S1, S2 )
*/
            ERRMSG( (char *)"mie_sph_Wis79 -- Self-test failed", true,
		    MSGLIM_thr,NUMMSG_errmsg_thr); // csz: Additional arguments for threading
           }
//                                      ** Restore user input values
         XX     = *XXSAV_thr;
         CREFIN = *CRESAV_thr;
         MIMCUT = *MIMSAV_thr;
         PERFCT = *PERSAV_thr;
         ANYANG = *ANYSAV_thr;
         NMOM   = *NMOSAV_thr;
         IPOLZN = *IPOSAV_thr;
         NUMANG = *NUMSAV_thr;
         XMU[ 1 ] = *XMUSAV_thr;

        }

  return;
} /* end TESTMI() */

bool TSTBAD( char *VarNam, double RelErr )
{
/*      Write name (VarNam) of variable failing self-test and its
        percent error from the correct value;  return  'FALSE'.
*/

  cout << " *** Output variable " << VarNam;
  cout << " differed by " << 100.*RelErr;
  cout << " per cent from correct value.  Self-test failed." << endl;

  return false;
} /* end TSTBAD() */

bool WRONG( double CALC, double EXACT )
{
  double ACCUR  = 1.e-4;
  return abs ( ( CALC - EXACT ) / EXACT ) > ACCUR;
} /* end WRONG() */

bool WRTBAD (char *VarNam,
             bool *MSGLIM_thr,int *NUMMSG_errmsg_thr,int *NUMMSG_wrtbad_thr) // csz: Additional arguments for threading
{
/*         Write names of erroneous variables and return 'TRUE'

       INPUT :   VarNam = Name of erroneous variable to be written
                          ( CHARACTER, any length )
*/

//    .. Local Scalars ..

  // static int NUMMSG = 0; // csz: Replace static value with local argument
    const int MAXMSG = 50; // csz: Replace static with const

  *NUMMSG_wrtbad_thr = *NUMMSG_wrtbad_thr + 1; // csz
  cout << " ****  Input variable  " << VarNam;
  cout << "  in error  ****" << endl;

  if ( *NUMMSG_wrtbad_thr == MAXMSG )
    ERRMSG( (char *)"Too many input errors.  Aborting...", true,
	    MSGLIM_thr,NUMMSG_errmsg_thr); // csz: Additional arguments for threading

  return true;
} /* end WRTBAD() */

bool WRTDIM (char *DimNam, int Minval)
{
/*         Write name of too-small symbolic dimension and
           the value it should be increased to;  return 'TRUE'

       INPUT :  DimNam = Name of symbolic dimension which is too small
                         ( CHARACTER, any length )
                Minval = Value to which that dimension should be
                         increased (at least)
*/

  cout << " ****  Symbolic dimension  " << DimNam;
  cout << "  should be increased to at least " << Minval << endl;

  return true;
} /* end WRTDIM() */

