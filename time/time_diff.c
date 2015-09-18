/*
The output from this program is...

The difference between the two dates...
Thu Jan  1 00:00:00 1970

...and...
Sun Jan  1 00:00:00 1995

...is 788918400 seconds, or 9131 days.

This result obviously does not include leap seconds.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

main()
{
 struct tm ts1, ts2;
 time_t t1, t2;

 ts1.tm_year = 70;   /* means the year 1970 */
 ts1.tm_mon = 0;     /* note that the months are 0, 1,..., 11 */
 ts1.tm_mday = 1;
 ts1.tm_hour = 0;
 ts1.tm_min = 0;
 ts1.tm_sec = 0;
 
 ts2.tm_year = 95;  /* means the year 1995 */
 ts2.tm_mon = 0;    /* note that the months are 0, 1,..., 11 */
 ts2.tm_mday = 1;
 ts2.tm_hour = 0;
 ts2.tm_min = 0;
 ts2.tm_sec = 0;
 
 t1 = mktime(&ts1);
 t2 = mktime(&ts2);
 
 printf("The difference between the two dates...\n");
 puts(ctime(&t1));
 printf("...and...\n");
 puts(ctime(&t2));
 printf("...is %.0f seconds, or %.0f days.\n", difftime(t2, t1), 
                                             difftime(t2, t1)/(float)86400);
}
