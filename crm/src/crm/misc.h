c     $Header$ -*-fortran-*-

#ifndef MISC_SET
# define MISC_SET
c#define NCPREC NF_FLOAT
#define NCPREC NF_DOUBLE
#define PVP 
#ifdef CRAY
#define REALTYPE MPI_REAL
#undef  SHELL_MSS
#undef  FORTFFT
#else /* not CRAY */
#define REALTYPE MPI_DOUBLE_PRECISION
#define SHELL_MSS
#define FORTFFT
#endif /* not CRAY */
#undef  COUP_SOM
#undef  COUP_CSM
#undef  SPMD
#if ( ! defined CRAY ) && ( ! defined LINUX ) && ( ! defined RS6K ) && ( ! defined SUN ) && ( ! defined SGI )
You must define one of CRAY, LINUX, RS6K, SUN, or SGI 
#endif /* not CCM_ARCH */ 
#endif /* not MISC_SET */


