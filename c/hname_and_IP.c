/*
 *----------------------------------------------------------------------
 * hname_and_IP.c
 *----------------------------------------------------------------------
 * Last written:
 * Time-stamp: "1995/02/16 17:08:21 thibaud@kether"
 *----------------------------------------------------------------------
 */

static char hname_and_IP_Id[] = "$Id$";

#include <sys/types.h>
#include <errno.h>
#include <malloc.h>
#include <netdb.h>
#include <netinet/in.h>
#include <stdio.h>
#include <string.h>
#include <sys/param.h>
#include <sys/socket.h>

#include <arpa/nameser.h>
#include <arpa/inet.h>
#include <resolv.h>

main( argc, argv, envp )
    int argc;				/* # of arguments on command line */
    char *argv[];			/* list arguments on command line */
    char *envp[];			/* list of environment variables */
{
    char *ch;				/* loop index */
    char *host_name;			/* name of the host */
    extern char *sys_errlist[];		/* translation of error codes */
    extern int errno;			/* error code # */
    extern int sys_nerr;		/* # of system error codes */
    int i;				/* loop index */
    struct hostent *host_str;		/* host structure */
    struct in_addr *aptr;		/* Internet address structure */

    if ( ( ( host_name = (char *)malloc( MAXHOSTNAMELEN ) ) == NULL ) ) {
        fprintf( stderr, "malloc failed\n" );
        return( -1  );
    }
    if ( ( host_str = (struct hostent *)malloc( sizeof( host_str ) ) )
        == NULL ) {
        fprintf( stderr, "malloc failed: %d, %s\n",
                errno, sys_errlist[errno] );
        exit( -1  );
    }

    /* get host name */
    if ( gethostname( host_name, MAXHOSTNAMELEN ) < 0 ) {
        fprintf( stderr, "gethostname failed: %d, %s\n",
                errno, sys_errlist[errno] );
        return( -1 );
    }

    /* seach for a `.' in host_name */
    for ( ch = host_name, i = 0;
	 *ch != '.' && i < (int)strlen( host_name );
	 ch++, i++ );
    if ( i < (int)strlen( host_name ) ) {
	/* there is a `.', assume host name include domain name */
	printf( "host name: %s\n", host_name );
    } else {
	/* there is NO `.', assume host name is the short one */
	/* initialize name server routine */
	res_init();
	strcat( host_name, "." );
	/* default IP domain */
	strcat( host_name, _res.defdname );
	printf( "full host name: %s\n", host_name );
    }

    /* get host structure */
    if ( ( host_str = gethostbyname( host_name ) ) == NULL ) {
        fprintf( stderr, "gethostbyname failed for %s: %d, %s\n",
                host_name, errno, sys_errlist[errno] );
        return( -1 );
    }
    /* official name */
    printf( "official host name: %s\n", host_str->h_name );
    /* aliases */
    while ( ( ch = *(host_str->h_aliases) ) != NULL ) {
        printf( "aliases: %s\n", ch );
        host_str->h_aliases++;
    }
    /* address type */
    printf( "address type: %d\n", host_str->h_addrtype );
    /* address length */
    printf( "address length: %d\n", host_str->h_length );
    /* IP addresses */
    while ( ( aptr = (struct in_addr *) *(host_str->h_addr_list) ) != NULL ) {
        printf( "IP address: %s\n", inet_ntoa( *aptr ) );
        host_str->h_addr_list++;
    }

    return( 0 );
}


/*
 *----------------------------------------------------------------------
 * RCS identification
 *----------------------------------------------------------------------
 * $Author: zender $
 * $Date$
 * $Locker:  $
 * $Revision$
 * $Source: /home/zender/cvs/c/hname_and_IP.c,v $
 * $State: Exp $
 * $Log: not supported by cvs2svn $
 * Revision 1.1.1.1  1998-08-31 01:25:20  zender
 * Imported sources
 *
 *----------------------------------------------------------------------
 * For GNU Emacs:
 *----------------------------------------------------------------------
 * Local Variables:
 * mode: C
 * abbrev-mode: t
 * comment-column: 40
 * version-control: t
 * End:
 *----------------------------------------------------------------------
 */
