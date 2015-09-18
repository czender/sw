// $Id$ 

// Purpose: Description (definition) of Wiscombe (1979, 1980, 1996) Mie scattering solutions and utilities 

/* Copyright (C) 2004--2014 Charlie Zender, Jorge Talamantes, Warren Wiscombe
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

// Usage:
// #include <mie_Wis79.hh> // Mie scattering solutions from Wis79

#ifndef MIE_WIS79_HH // Contents have not yet been inserted in current source file
#define MIE_WIS79_HH

// Standard C++ headers 
#include <complex> // Standard C++ complex class
#include <iomanip> // Standard C++ I/O manipulation: setw()
#include <iostream> // Standard C++ I/O streams: cout, cin, cerr

// Standard C headers 
#include <cmath> // sin cos cos sin 3.14159 
#include <cstdlib> // strtod, strtol, malloc, getopt  

// 3rd party vendors

// Personal headers
#include <dbg.hh> // Debugging constants
#include <mth.hh> // Mathematical utilities, constants */
#include <utl.hh> // err_prn(), wrn_prn(), Exit_gracefully()

// Namespaces
using namespace std;

// Typedefs
typedef complex<double> dcomp;

// Forward declarations
const int PMOM_second_dimension(1); // [nbr] Second dimension of Legendre moment arrays (defined in MIEV.doc, set to 1 for unpolarized light, must agree with sizes PMOM_second_dimension size defined in mie_Wis79.cc)

// Prototype functions with C++ linkages
void BIGA ( dcomp CIOR, double XX, int NTRM, bool NOABS, bool YESANG,
            double RBIGA[], int RBIGA_dimension,
            dcomp CBIGA[], int CBIGA_dimension,
            bool *MSGLIM_thr,int *NUMMSG_errmsg_thr); // csz: Additional arguments for threading
void CKINMI ( // input variables:
              int NUMANG, int MAXANG, double XX, bool PERFCT,
              dcomp CREFIN, int MOMDIM, int NMOM, int IPOLZN,
              int ANYANG, double XMU[], int XMU_dimension,
              // output variables:
              bool CALCMO[], int CALCMO_dimension, int &NPQUAN,
             bool *MSGLIM_thr,int *NUMMSG_errmsg_thr,int *NUMMSG_wrtbad_thr); // csz: Additional arguments for threading


dcomp CONFRA ( int N, dcomp ZINV,
            bool *MSGLIM_thr,int *NUMMSG_errmsg_thr); // csz: Additional arguments for threading

void ERRMSG( char *MESSAG, bool FATAL,
	     bool *MSGLIM_thr,int *NUMMSG_errmsg_thr); // csz: Additional arguments for threading

double F1 ( double MRE);

double F2 ( double MRE);

double F3 ( double MRE);

void LPCO1T ( int NMOM, int IPOLZN, int MOMDIM,
              bool CALCMO[], int CALCMO_dimension,
              dcomp A[], int A_dimension, dcomp B[], int B_dimension,
              double PMOM[][PMOM_second_dimension], int PMOM_dimension );

void LPCO2T ( int NMOM, int IPOLZN, int MOMDIM,
              bool CALCMO[], int CALCMO_dimension,
              dcomp A[], int A_dimension, dcomp B[], int B_dimension,
              double PMOM[][PMOM_second_dimension], int PMOM_dimension );

void LPCOEF( int NTRM, int NMOM, int IPOLZN, int MOMDIM, bool CALCMO[],
             int CALCMO_dimension, int NPQUAN,
             dcomp A[], int A_dimension, dcomp B[], int B_dimension,
             double PMOM[][PMOM_second_dimension], int PMOM_dimension,
	     bool *PASS1_lpcoef_thr,double *RECIP_thr, // csz: Additional arguments for threading
             bool *MSGLIM_thr,int *NUMMSG_errmsg_thr,int *NUMMSG_wrtbad_thr); // csz: Additional arguments for threading

int // O [enm] Return success code
mie_sph_Wis79(  // [fnc] Mie solution for homogeneous spheres (MIEV0), Wis79
// input variables:
       double XX, dcomp CREFIN, bool PERFCT, double MIMCUT,
       bool ANYANG, int NUMANG, double XMU[], int XMU_dimension,
       int NMOM, int IPOLZN, int MOMDIM, bool PRNT[], int PRNT_dimension,
// output variables:
       double &QEXT, double &QSCA, double &GQSC,
       double PMOM[][PMOM_second_dimension], int PMOM_dimension, dcomp &SFORW,
       dcomp &SBACK, dcomp S1[], int S1_dimension, dcomp S2[], int S2_dimension,
       dcomp TFORW[], int TFORW_dimension, dcomp TBACK[], int TBACK_dimension,
       double &SPIKE );

int MIEV0_drv();

void PARSEI ( int IPOLZN, int N[], int N_dimension);

void SMALL1( double XX, int NUMANG, double XMU[], int XMU_dimension,
             double &QEXT, double &QSCA, double &GQSC, dcomp &SFORW,
             dcomp &SBACK, dcomp S1[], int S1_dimension, dcomp S2[],
             int S2_dimension, dcomp TFORW[], int TFORW_dimension,
             dcomp TBACK[], int TBACK_dimension,
             dcomp A[], int LITA_dimension,
             dcomp B[], int LITB_dimension );

void SMALL2( double XX, dcomp CIOR, bool CALCQE, int NUMANG,
             double XMU[], int XMU_dimension,
             double &QEXT, double &QSCA, double &GQSC, dcomp &SFORW,
             dcomp &SBACK, dcomp S1[], int S1_dimension, dcomp S2[],
             int S2_dimension, dcomp TFORW[], int TFORW_dimension,
             dcomp TBACK[], int TBACK_dimension,
             dcomp A[], int LITA_dimension,
             dcomp B[], int LITB_dimension );

double SQ ( dcomp );

void TESTMI( bool COMPAR, double &XX, dcomp &CREFIN, double &MIMCUT,
             bool &PERFCT, bool &ANYANG, int &NMOM, int &IPOLZN, int &NUMANG,
             double XMU[], int XMU_dimension, double &QEXT, double &QSCA,
             double &GQSC, dcomp &SFORW, dcomp &SBACK,
             dcomp S1[], int S1_dimension, dcomp S2[], int S2_dimension, 
             dcomp TFORW[], int TFORW_dimension, dcomp TBACK[],
             int TBACK_dimension, double PMOM[][PMOM_second_dimension], int PMOM_dimension,
             int MOMDIM, 
	     bool *ANYSAV_thr,dcomp *CRESAV_thr,int *IPOSAV_thr,double *MIMSAV_thr,int *NMOSAV_thr,int *NUMSAV_thr,bool *PERSAV_thr,double *XMUSAV_thr,double *XXSAV_thr, // csz: Additional arguments for threading
             bool *MSGLIM_thr,int *NUMMSG_errmsg_thr,int *NUMMSG_wrtbad_thr); // csz: Additional arguments for threading

bool TSTBAD( char *VarNam, double RelErr );

bool WRONG( double CALC, double EXACT );

bool WRTBAD (char *VarNam,
	     bool *MSGLIM_thr,int *NUMMSG_errmsg_thr,int *NUMMSG_wrtbad_thr); // csz: Additional arguments for threading

bool WRTDIM (char *DimNam, int Minval);


#endif // MIE_WIS79_HH  
