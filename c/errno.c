/* Purpose: Example of using errno from François P. Thibaud <thibaud@ucar.edu> */
#include <errno.h>
#include <stdio.h>
main(){
    extern char *sys_errlist[];		/* translation of error codes */
    extern int errno;			/* error code # */
    extern int sys_nerr;		/* # of system error codes */
    
    int i=1;
    int j;

    j=i/0;
    fprintf( stderr, "XXX failed: %d, %s\n", errno, sys_errlist[errno] );
} /* end errno */


