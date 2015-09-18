/* $Id$ */

/* Purpose:  Print network information about the host */

/* Parts of utsname.c are based on hname_and_IP.c by Francois Thibaud */

/* Example Usage (place mousable command lines here):
   utsname
 */

#define bool int

/* Standard header files */
#include <math.h> /* sin cos cos sin 3.14159 */
#include <stdio.h> /* stderr, FILE, NULL, etc. */
#include <stdlib.h> /* atof, atoi, malloc, getopt */
#include <string.h> /* strcmp. . . */
#include <time.h> /* machine time */
#include <unistd.h> /* all sorts of POSIX stuff */
#include <sys/utsname.h> /* defines struct utsname */

#include <pwd.h> /* password structures for getpwuid() */
#include <sys/types.h> /* needed for _res */
#include <netinet/in.h> /* needed for _res */
#if (defined ALPHA) || (defined OSF1)
#include <sys/param.h> /* needed for MAXHOSTNAMELEN */
#endif /* not ALPHA || OSF1 */
#ifndef WIN32
#include <arpa/nameser.h> /* needed for _res */
#include <netdb.h> /* needed for MAXHOSTNAMELEN */
#include <resolv.h> /* Internet structures for _res */
#endif /* WIN32 */

#if (defined LINUX) || (defined LINUXALPHA)
#include <getopt.h> /* GNU getopt() is standard on Linux */
#else /* not LINUX */
#include "getopt.h" /* GNU getopt() */
#endif /* not LINUX */

/* 3rd party vendors */

/* Personal headers */
#include <dbg.h> /* Debugging constants */

/* Global variables declared here */
unsigned short int dbg_lvl=0; /* Option D */
unsigned short int dbg_val=0; /* Option d */

int main(int argc,char **argv)
{
  char *cmd_ln_sng(int,char **);
  void Exit_gracefully(void);
  void usg_prn(char *);

  char *cmd_ln; /* [sng] Parsed command line */
  char *dmn_nm; /* [sng] Domain name */
  char *fl_in; /* [sng] Input file */
  char *fl_out; /* [sng] Output file */
  char *hst_nm; /* [sng] Host name */
  char *prg_nm; /* [sng] Name of program */
  char *time_bfr_srt; /* [sng] Current date and time */

  char CVS_Id[]="$Id$"; /* [sng] CVS identification string */
  char CVS_Revision[]="$Revision$"; /* [sng] CVS revision string */
  char sbr_nm[]="main()"; /* [sng] Subroutine name */

  double dbl_foo=0.0; /* [frc] Intrinsic double temporary variable */

  float flt_foo=0.0; /* [frc] Intrinsic float temporary variable */

  int idx; /* Counting index */

  short sht_foo=73; /* [frc] Intrinsic short temporary variable */

  struct hostent *hst_ntr; /* Host entry */
  struct utsname unm;
  struct passwd *usr_pwd;

  time_t time_crr_time_t; /* [tm] Current date and time */

  uid_t usr_uid;
  uid_t usr_euid;

  /* Variables for option processing */
  char *opt_crr; /* String representation of current long-option name */
  char *opt_sng=NULL; /* String representation of current optarg, if any */ /* CEWI */
  extern char *optarg; /* char * representation of current optarg, if any (this memory is owned by system) */
  extern int optind; /* extern enumerating cardinal of current option */
  int opt; /* Value is zero if current argument is long type, else value contains single letter version of command line argument */
  int opt_idx=0; /* Index of current long option into opt_lng array */
  /* Short options: no colon = no arg, one colon = required arg, two colons = optional arg */
  char opt_sht_lst[]="D:f:o:v"; /* List of single-letter (C-style) option abbreviations */
  static struct option opt_lng[]=
  {
    /* Option structure is {char *name,int has_arg,int *flag,int val} 
       has_arg is compared to enum _argtype{no_argument,required_argument,optional_argument}, 
       flag points to a variable that gets set to val whenever the name option is set.
       For long options that have zero flag field, getopt() returns contents of val
    */
    /* Long options with no argument */
    /* Long options with argument */
    {"dbl_foo",1,0,0}, /* [frc] Intrinsic double temporary variable */
    {"sht_foo",1,0,0}, /* [frc] Intrinsic short temporary variable */
    /* Long options with optional argument */
    /* Long options with short counterparts */
    {"dbg_lvl",1,0,'D'}, /* [enm] Debugging level */
    {"debug",1,0,'D'}, /* [enm] Debugging level */
    {"flt_foo",1,0,'f'}, /* [frc] Intrinsic float temporary variable */
    {"output",1,0,'o'},
    {"version",0,0,'v'},
    {"verbose",0,0,'D'},
    /* Last option named "0" to signal getopt_long() to stop processing */
    {0,0,0,0}
  }; /* end opt_lng */

  /* Set defaults */

  /* Start clock and save command line */
  cmd_ln=cmd_ln_sng(argc,argv); /* [sng] Parsed command line */
  time_crr_time_t=time((time_t *)NULL); /* [tm] Current date and time */
  time_bfr_srt=ctime(&time_crr_time_t); /* [tm] Current date and time */
  (void)fprintf(stderr,"\tStart = %s",time_bfr_srt);
  prg_nm=((prg_nm=strrchr(argv[0],'/')) == NULL) ? argv[0] : prg_nm++; /* [sng] Name of program */

  /* Parse command line arguments */
  while(1){
    /* getopt_long_only() allows a single dash '-' to prefix long options as well */
    opt=getopt_long_only(argc,argv,opt_sht_lst,opt_lng,&opt_idx);
    /* NB: access to opt_crr is only valid when long_opt was detected */
    opt_crr=(char *)strdup(opt_lng[opt_idx].name);
    if(optarg) opt_sng=optarg; /* Copy system memory into local string for safer operations */
    if(opt == EOF) break; /* Parse positional arguments once getopt_long_only() returns EOF */
    /* Process long options without short option counterparts */
    if(opt == 0){
      if(dbg_lvl >= dbg_io){
	(void)fprintf(stdout,"Long option name: %s, ",opt_lng[opt_idx].name);
	if(optarg) (void)fprintf(stdout," Argument: %s\n",opt_sng); else (void)fprintf(stdout," No argument\n");
      } /* end if dbg */
      if(!strcmp(opt_crr,"dbl_foo")) dbl_foo=strtod(opt_sng,(char **)NULL);
      if(!strcmp(opt_crr,"flt_foo")) flt_foo=(float)strtod(opt_sng,(char **)NULL);
      if(!strcmp(opt_crr,"sht_foo")) sht_foo=(short)strtol(opt_sng,(char **)NULL,10);
    } /* opt != 0 */
    switch(opt){
    case 0: /* Long options have already been processed, return */
      break;
    case 'D': /* Debugging level.  Default is 0. */
      if(optarg) dbg_lvl=(unsigned short int)strtol(optarg,(char **)NULL,10); else dbg_lvl=dbg_fl;
      break;
    case 'f': /* Generic tuning parameter.  Default is 0. */
      flt_foo=(float)strtod(optarg,(char **)NULL);
      break;
    case 'o': /* Get the output file name. Default is stdout */
      fl_out=opt_sng;
      break;
    case 'v': /* CVS program info */
      (void)fprintf(stderr,"%s %s\n",CVS_Revision,CVS_Id);
      exit(EXIT_SUCCESS);
      break;
    default: /* Print proper usage */
      (void)usg_prn(opt_sht_lst);
      exit(EXIT_FAILURE);
    } /* end switch */
  } /* end while loop */
  
  /* Process positional arguments */
  if(optind < argc){
    int psn_arg_nbr=argc-optind;
    if(psn_arg_nbr > 2){
      (void)fprintf(stdout,"%s: ERROR %s reports too many positional arguments\n",prg_nm,sbr_nm);
      exit(EXIT_FAILURE);
    }else if(psn_arg_nbr == 1){
      fl_out=argv[optind++];
    }else if(psn_arg_nbr == 2){
      fl_in=argv[optind++];
      fl_out=argv[optind++];
    } /* end else */
  } /* end if */

  /* Main body of code */
  (void)fprintf(stdout,"\n1. Contents of uname(2) system call with utsname structure:\n");
#ifdef LINUX
  (void)fprintf(stdout,"Size limit of system name: SYS_NMLN = _UTSNAME_LENGTH = %d\n",_UTSNAME_LENGTH);
  (void)fprintf(stdout,"Size limit of system names: _UTSNAME_NODENAME_LENGTH = %d\n",_UTSNAME_NODENAME_LENGTH);
  (void)fprintf(stdout,"Size limit of system names: _UTSNAME_DOMAIN_LENGTH = %d\n",_UTSNAME_DOMAIN_LENGTH);
#endif /* LINUX */
#if ( defined SUN4 ) || ( defined SUN4SOL2 ) || ( defined SUNMP ) || ( defined SGI64 )
  (void)fprintf(stdout,"Size limit of system name: SYS_NMLN = %d\n",SYS_NMLN);
#endif /* Solaris, IRIX64 */
  (void)uname(&unm);
  (void)fprintf(stdout,"utsname.sysname = %s\n",unm.sysname);
  (void)fprintf(stdout,"utsname.nodename = %s\n",unm.nodename);
  (void)fprintf(stdout,"utsname.release = %s\n",unm.release);
  (void)fprintf(stdout,"utsname.version = %s\n",unm.version);
  (void)fprintf(stdout,"utsname.machine = %s\n",unm.machine);
  /* NB: domainname element is not in Solaris utsname struct as of 1998/01/13 */
#ifdef LINUX
#ifdef __USE_GNU
  (void)fprintf(stdout,"utsname.domainname = %s\n",unm.domainname);
#else /* __USE_GNU */
  (void)fprintf(stdout,"utsname.__domainname = %s\n",unm.__domainname);
#endif /* __USE_GNU */
#endif /* LINUX */

  /* Get information about current user */
  usr_uid=getuid();
  usr_euid=geteuid();
  usr_pwd=getpwuid(usr_uid);
  (void)fprintf(stdout,"\n1. getuid(2), geteuid(2), and getpwuid(2) system calls:\n");
  (void)fprintf(stdout,"getuid() returns %d\n",usr_uid);
  (void)fprintf(stdout,"geteuid() returns %d\n",usr_euid);
  (void)fprintf(stdout,"user name = getpwuid()->passwd.pw_name: %s\n",usr_pwd->pw_name);
  (void)fprintf(stdout,"user password = getpwuid()->passwd.pw_passwd: %s\n",usr_pwd->pw_passwd);
  (void)fprintf(stdout,"user ID = getpwuid()->passwd.pw_uid: %d\n",usr_pwd->pw_uid);
  (void)fprintf(stdout,"group ID = getpwuid()->passwd.pw_gid: %d\n",usr_pwd->pw_gid);
  (void)fprintf(stdout,"real name = getpwuid()->passwd.pw_gecos: %s\n",usr_pwd->pw_gecos);
  (void)fprintf(stdout,"home directory = getpwuid()->passwd.pw_dir: %s\n",usr_pwd->pw_dir);
  (void)fprintf(stdout,"login shell = getpwuid()->passwd.pw_shell: %s\n",usr_pwd->pw_shell);

  /* Get information about current machine */
  dmn_nm=(char *)malloc(MAXHOSTNAMELEN*sizeof(char));
  hst_nm=(char *)malloc(MAXHOSTNAMELEN*sizeof(char));
  (void)getdomainname(dmn_nm,MAXHOSTNAMELEN);
  (void)gethostname(hst_nm,MAXHOSTNAMELEN);
  (void)fprintf(stdout,"\n2. getdomainname(2) and gethostname(2) system calls:\n");
  (void)fprintf(stdout,"Size limit of host name: MAXHOSTNAMELEN = %d\n",MAXHOSTNAMELEN);
  (void)fprintf(stdout,"gethostname() returns %s\n",hst_nm);
  (void)fprintf(stdout,"getdomainname() returns %s\n",dmn_nm);

  /* Get information about current machine from hostent */
  (void)fprintf(stdout,"\n3. hostent structure:\n");
  hst_ntr=gethostbyname(hst_nm); /* [sct] Host entry */
  (void)fprintf(stdout,"hostname = hst_ntr->h_name = %s\n",hst_ntr->h_name);
  (void)fprintf(stdout,"address length = hst_ntr->h_length = %i\n",hst_ntr->h_length);
  (void)fprintf(stdout,"address type = hst_ntr->h_addrtype = %i\n",hst_ntr->h_addrtype);
  idx=0;
  /* Aliases list is zero-terminated meaning end of list is marked by two successive NULs */
  while(hst_ntr->h_aliases[idx][0] != '\0'){ /* struct hostent in <netdb.h> specifies char **h_addr_list as zero-terminated a= NUL-terminated array of network addresses */
    /*    int chr_idx=0;
    while(hst_ntr->h_aliases[idx][chr_idx] != '\0'){
      (void)fprintf(stdout,"h_aliases[%i][%i] = %c = %i\n",idx,chr_idx,hst_ntr->h_aliases[idx][chr_idx],(signed int)hst_ntr->h_aliases[idx][chr_idx]);
      chr_idx++;
      } */ /* end loop over idx */
    /* "zero-terminated" means a hex 0 = '\0' = NUL not character '0' 
       Terminating with '0' easily causes string formatting routines to core dump */
    /* fxm: LINUX 20010730 core dumps on following line */
    /* (void)fprintf(stdout,"alias #%d = hst_ntr->h_aliases.[%d] = %s\n",(unsigned int)idx+1,(unsigned int)idx,hst_ntr->h_aliases[idx]); */
    idx++;
  } /* end loop over idx */

  /* Get information about current machine from _res */
  (void)fprintf(stdout,"\n4. _res system calls:\n");
  (void)res_init();
  (void)fprintf(stdout,"default domain = _res.defdname = %s\n",_res.defdname); /* Default domain name */
  (void)fprintf(stdout,"# of name servers = _res.nscount = %d\n",_res.nscount); /* Number of name servers */
  idx=0;
  for(idx=0;idx<_res.nscount;idx++){
    /* _res member nsaddr_list is a list of struct sockaddr_in (defined in <netinet/in.h>) */
    /* fxm: LINUX 20010730 unable to compile following line */
    /* (void)fprintf(stdout,"name server #%d = _res.nsaddr_list.[%d] = %i\n",(unsigned int)idx+1,(unsigned int)idx,(int)(_res.nsaddr_list[idx].sin_addr)); */
    idx++;
  } /* end loop over idx */

  Exit_gracefully();
  return EXIT_SUCCESS;
} /* end main() */

void usg_prn(char *opt_sng)
{
  (void)fprintf(stderr,"\nusage: c_template [-options] where options are one or more of:\n\n");
  (void)fprintf(stderr,"%s\n\n",opt_sng);
  (void)fprintf(stderr,"-D dbg_lvl The debugging level.  Default is 0\n");
  (void)fprintf(stderr,"-E STDERR Toggle stderr stream.  Default is True\n");
  (void)fprintf(stderr,"-I STDIN Toggle stdin stream.  Default is True\n");
  (void)fprintf(stderr,"-O STDOUT Toggle stdout stream.  Default is True\n");
  (void)fprintf(stderr,"-d dbg_val The second debug level.  Default is 0\n");
  (void)fprintf(stderr,"-f float_foo Set the generic float. Default is 0.\n");
  (void)fprintf(stderr,"-v print the CVS program version\n");
  (void)fprintf(stderr,"\n");
} /* end usg_prn() */

void Exit_gracefully(void)
{
  char *time_bfr_end;
  time_t clock;

  /* end the clock */
  
  clock=time((time_t *)NULL);
  time_bfr_end=ctime(&clock);
  (void)fprintf(stderr,"\tfinish = %s\n",time_bfr_end);

  (void)fclose(stderr);
  (void)fclose(stdin);
  (void)fclose(stdout);

  exit(0);
} /* end Exit_gracefully() */

char *
cmd_ln_sng(int argc,char **argv)
/* 
   int argc: input argument count
   char **argv: input argument list
   char *cmd_ln_sng(): output command line
*/
{
  char *cmd_ln;
  
  int cmd_ln_sz=0;
  int idx;

  for(idx=0;idx<argc;idx++){
    cmd_ln_sz+=(int)strlen(argv[idx])+1;
  } /* end loop over args */
  cmd_ln=(char *)malloc(cmd_ln_sz*sizeof(char));
  if(argc <= 0){
    cmd_ln=(char *)malloc(sizeof(char));
    cmd_ln[0]='\0';
  }else{
    (void)strcpy(cmd_ln,argv[0]);
    for(idx=1;idx<argc;idx++){
      (void)strcat(cmd_ln," ");
      (void)strcat(cmd_ln,argv[idx]);
    } /* end loop over args */
  } /* end else */

  return cmd_ln;
} /* end cmd_ln_sng() */
      




