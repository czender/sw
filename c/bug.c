#include <stdio.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h> /* OpenMP pragmas */
#endif /* not _OPENMP */

int main()
{
  float fabsf(float);
  float cosf(float);
  
  float minus_one_f=-1.0;
  double minus_one_d=-1.0;
  float float_cosf_minus_one_f,float_fabsf_minus_one_f;
  double double_cos_minus_one_d,double_fabs_minus_one_d;
  
  /* Produces correct results */
  double_fabs_minus_one_d=fabs(minus_one_d);
  double_cos_minus_one_d=cos(minus_one_d);
  fprintf(stdout,"abs(-1.0d) = double_fabs_minus_one_d = %g\n",double_fabs_minus_one_d);
  fprintf(stdout,"cos(-1.0d) = double_cos_minus_one_d = %g\n",double_cos_minus_one_d);
  
  /* Produces incorrect results on AIX */
  float_fabsf_minus_one_f=fabsf(minus_one_f);
  float_cosf_minus_one_f=cosf(minus_one_f);
  fprintf(stdout,"fabsf(-1.0f) = float_fabsf_minus_one_f = %g\n",float_fabsf_minus_one_f);
  fprintf(stdout,"cosf(-1.0f) = float_cosf_minus_one_f = %g\n",float_cosf_minus_one_f);
  
  FILE * const fp_stderr=stderr; /* [fl] stderr filehandle CEWI */
#ifdef _OPENMP /* OpenMP-compliant compilers define _OPENMP=YYYYMM = year and month of OpenMP specification */
  int thr_nbr_max; /* [nbr] Maximum number of threads system/user allow program to use */  
  thr_nbr_max=omp_get_max_threads(); /* [nbr] Maximum number of threads system/user allow program to use */
#pragma omp parallel default(none) shared(fp_stderr)
  { /* begin OpenMP parallel */
#pragma omp single nowait
    { /* begin OpenMP single */
      (void)fprintf(fp_stderr,"c: INFO OpenMP multi-threading using %d threads\n",omp_get_num_threads());
    } /* end OpenMP single */
  } /* end OpenMP parallel */
#else /* not _OPENMP */
  (void)fprintf(fp_stderr,"c: INFO Not attempting OpenMP multi-threading\n");
#endif /* not _OPENMP */

} /* end main() */

/* Provide substitutes */
/*float fabsf(float x){return (float)(fabs((double)x));}*/
/*float cosf(float x){return (float)(cos((double)x));}*/

/* Purpose: Demonstrate handling of float-valued math intrinsics 
   scp ~/c/bug.c esmf.ess.uci.edu:c
   
   AIX/xlc: 
   System (libC.a) versions work iff fabsf/cosf are explicitly prototyped in main()
   Private versions work iff fabsf/cosf are explicitly prototyped in main()
   xlc_r -o bug bug.c -lm
   xlc_r -o bug bug.c -lm -lC
   xlc_r -qlanglvl=extended -o bug bug.c -lm -lC
   xlc_r -qlanglvl=extended -qsmp=omp -o bug bug.c -lm -lC
   
   Linux/GCC: 
   Compiling with
   gcc -Wconversion -o bug bug.c -lm
   give this warning: 
   bug.c:19: warning: passing arg 1 of `fabsf' as `float' rather than `double' due to prototype
   Adding (redundant, since it's in math.h) prototype 
   float fabsf(float);
   does not help
*/
