/*
 * elapse.c
 *
 * compile on Sun with : cc -O -v -o elapse elapse.c -lm
 *            IBM with : cc -O -v -DSVR3 -o elapse elapse.c -lm
 *           Cray with : cc -O -## -DSVR3 -o elapse elapse.c -lm
 *
 */

#include        <sys/types.h>
#include	<sys/param.h>

#ifdef SVR3
#include	<time.h>
#include	<sys/time.h>
#else
#include	<sys/timeb.h>
#endif

#include	<sys/times.h>
#include	<math.h>

struct cpu_elapse {
#ifdef SVR3
	time_t		elapse;
#else
	struct timeb	elapse;	/* defined in <sys/timeb.h> */
#endif
	struct tms	cpu;	/* defined in <sys/times.h> */
}	time_buffer[2][64];	/* 64 pairs; the [0][n] is set by
				 * cpu_start, the [1][n] is set by
				 * cpu_end which then difference
				 * both; 64 is an arbitrary value
				 */

#define		IMAX		1000000

main()
{
	double	elt;	/* elapse (wall clock) time */
	double	put;	/* user time, parent */
	double	pst;	/* system time, parent */
	double	cut;	/* user time, childrens */
	double	cst;	/* system time, childrens */
	double	cpu;	/* total CPU time */
	double	rate;	/* total CPU time / elapse time */
	float	var[IMAX + 1];
	int	i;
	int	counter;

	printf( "\n" );
#ifdef SVR3
	printf( "the time granularity is CLK_TCK = %d\n", CLK_TCK );
	printf( "the smallest time increment is %es\n", 1.0 / (float)CLK_TCK );
#else
	printf( "the time granularity is HZ = %d\n", HZ );
	printf( "the smallest time increment is %es\n", 1.0 / (float)HZ );
#endif
	printf( "\n" );

	counter = 13;	/* for example */
	if( cpu_start( &counter ) < 0 ) {
		printf( "cpu_start failed.\n" );
		exit( -1 );
	}

	for( i = 0; i <= IMAX; i++ ) {
		var[i] = (float)sqrt( (double)i );
	}
	printf( "var[%d] is: %e\n", IMAX, var[IMAX] );

	if( cpu_end( &counter, &elt, &put, &pst, &cut, &cst ) < 0 ) {
		printf( "cpu_end failed.\n" );
		exit( -1 );
	}

	printf( "\n" );
	printf( "elapse (wall clock) time : %f\n", elt );
	printf( "user time, parent        : %f\n", put );
	printf( "system time, parent      : %f\n", pst );
	printf( "user time, childrens     : %f\n", cut );
	printf( "system time, childrens   : %f\n", cst );
	cpu = put + pst + cut + cst;
	printf( "total CPU time           : %f\n", cpu );
	rate = 100.0 * ( cpu / elt );
	printf( "rate (cpu/elapse)        : %f%%\n", rate );
	printf( "\n" );

	exit( 0 );
}

/*
 * cpu_start initialize wall clock (elapse) time and CPU time counter
 *
 * the argument is made a pointer just to be able to call the
 * routine from Fortran without headache ... 8=)
 */
int
cpu_start( counter )
int	*counter;
{
#ifdef SVR3
	if ( ( time_buffer[0][*counter].elapse =
	    times( &(time_buffer[0][*counter].cpu) ) ) < 0 ) {
		return( -1 );
	}
#else
	if ( ftime( &(time_buffer[0][*counter].elapse) ) < 0 ) {
		return( -1 );
	}
	if ( times( &(time_buffer[0][*counter].cpu) ) < 0 ) {
		return( -1 );
	}
#endif
	return( 0 );
}

/*
 * cpu_end used the initialized wall clock (elapse) time
 * and CPU time counter to difference then with the current
 * state and return the differences as double
 *
 * the argument are made a pointers just to be able to call the
 * routine from Fortran without headache ... 8=)
 */
int
cpu_end( counter, elapse_time, parent_user, parent_syst, childs_user,
    childs_syst )
int	*counter;	/* one might want to have more then one cpu
			 * and elapse time counter ...
			 */
double	*elapse_time;	/* wall clock (or elapse time) */
double	*parent_user;	/* user time */
double	*parent_syst;	/* system time */
double	*childs_user;	/* user time, children */
double	*childs_syst;	/* system time, children */
{
#ifdef SVR3
	if ( ( time_buffer[1][*counter].elapse =
	    times( &(time_buffer[1][*counter].cpu) ) ) < 0 ) {
		return( -1 );
	}

        *elapse_time = (double)(time_buffer[1][*counter].elapse -
            time_buffer[0][*counter].elapse) / (double)CLK_TCK;

	*parent_user = (double)(time_buffer[1][*counter].cpu.tms_utime -
	    time_buffer[0][*counter].cpu.tms_utime) / (double)CLK_TCK;
	*parent_syst = (double)(time_buffer[1][*counter].cpu.tms_stime -
	    time_buffer[0][*counter].cpu.tms_stime) / (double)CLK_TCK;
	*childs_user = (double)(time_buffer[1][*counter].cpu.tms_cutime -
	    time_buffer[0][*counter].cpu.tms_cutime) / (double)CLK_TCK;
	*childs_syst = (double)(time_buffer[1][*counter].cpu.tms_cstime -
	    time_buffer[0][*counter].cpu.tms_cstime) / (double)CLK_TCK;
#else
	if ( ftime( &(time_buffer[1][*counter].elapse) ) < 0 ) {
		return( -1 );
	}
	if ( times( &(time_buffer[1][*counter].cpu) ) < 0 ) {
		return( -1 );
	}

	*elapse_time = 0.001 * (double)(
	      (long)time_buffer[1][*counter].elapse.time * 1000L
	    + (long)time_buffer[1][*counter].elapse.millitm
	    - (long)time_buffer[0][*counter].elapse.time * 1000L
	    - (long)time_buffer[0][*counter].elapse.millitm );

	*parent_user = (double)(time_buffer[1][*counter].cpu.tms_utime -
	    time_buffer[0][*counter].cpu.tms_utime) / HZ;
	*parent_syst = (double)(time_buffer[1][*counter].cpu.tms_stime -
	    time_buffer[0][*counter].cpu.tms_stime) / HZ;
	*childs_user = (double)(time_buffer[1][*counter].cpu.tms_cutime -
	    time_buffer[0][*counter].cpu.tms_cutime) / HZ;
	*childs_syst = (double)(time_buffer[1][*counter].cpu.tms_cstime -
	    time_buffer[0][*counter].cpu.tms_cstime) / HZ;
#endif
	return( 0 );
}
