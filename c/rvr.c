/* $Id$ */

/* Purpose: Reverse a file byte-by-byte
   Applying operator twice restores original file
   This operator is needed to overcome downloading problems on some executables 
   The need for this program became evident when I was trying to download Netscape binaries
   The binaries are very long and the FTP was always stalling in the same place
   Even after splitting the file with split (1), the FTP would stop on the exact same byte
   Presumably there is a code sequence which is causing my modem to stall
   I did not know this was possible with UNIX, because I thought the binary data would only be
   seen by the FTP protocol and be successfully hidden from the modem hardware.
   It must be that the executable has an unfortunate string of bytes which matches a modem control sequence
   The fact that I can always FTP files successfullly from NCAR, but not from home, 
   convinces me this must be a modem-related problem, and not some form of copy protection.
 */

/* Example Usage (place mousable command lines here):
cat >&! foo.tar.gz
1234567890
rvr -D 1 foo.tar.gz oof.tar.gz
split -b 1m oof.tar.gz 
cat x?? >! oof.tar.gz
rvr -D 1 oof.tar.gz foo.tar.gz

cat xaa xab >! oof.tar.gz
rvr -D 1 xaf fax
rvr -D 1 fax xaf
 */

#define bool int

/* Standard header files */
#include <stdio.h> /* stderr, FILE, NULL, etc. */
#include <time.h> /* machine time */
#include <string.h> /* strcmp. . . */
#include <sys/stat.h> /* stat() */
#include <stdlib.h> /* atof, atoi, malloc, getopt */
#include <unistd.h> /* all sorts of POSIX stuff */

#if (defined LINUX) || (defined LINUXALPHA)
#include <getopt.h> /* GNU getopt() is standard on Linux */
#else /* not LINUX */
#include "getopt.h" /* GNU getopt() */
#endif /* not LINUX */

/* Global variables declared here */
char *prg_nm;

unsigned short int dbg_lvl=0; /* Option D */
unsigned short int dbg_val=0; /* Option d */

int main(int argc,char **argv)
{
  char *cmd_ln_sng(int,char **);
  void Exit_gracefully(void);
  void usg_prn(char *);

  FILE *fp_in;
  FILE *fp_out;

  char *bff_in;
  char *bff_out;
  char *cmd_ln;
  char *fl_in=NULL;
  char *fl_out=NULL;
  char *opt_sng;
  char *time_buf_srt;

  char rcs_Id[] = "$Id$\n";
  char rcs_Revision[] = "$Revision$\n";
  
  extern char *optarg;
  extern int optind;
  
  long bff_sz;
  long idx;

  int opt;
  int psn_arg_nbr;
  int rcd=0;
  
  struct stat stat_sct;

  time_t clock;

  /* set defaults */
  bff_sz=-1L; /* Option b */

  /* Start the clock and save the command line */
  cmd_ln=cmd_ln_sng(argc,argv);
  clock=time((time_t *)NULL);
  time_buf_srt=ctime(&clock);
  (void)fprintf(stderr,"\tStart = %s",time_buf_srt);
  (void)fprintf(stderr,"\tCommand line = %s\n",cmd_ln);

  /* Get the program name */
  prg_nm=argv[0];

  /* Parse command line arguments */
  while(1){
    int opt_idx = 0;
    static struct option opt_long[] =
      {
	/* The option structure is {char *name,int has_arg,int *flag,int val} 
	   has_arg is compared to enum _argtype{no_argument,required_argument,optional_argument}, 
	   flag points to a variable that gets set to val whenever the name option is set.
	   For long options that have a zero flag field, getopt() returns the contents of val.
	 */
	{"debug",1,0,'D'},
	{"output",1,0,'o'},
	{"version",0,0,'v'},
	{"verbose",0,0,'D'},
	/* The last option must have a name of "0" to signal to getopt_long() to stop processing */
	{0,0,0,0}
      };
    
    opt_sng="b:D:o:v";
    opt=getopt_long_only(argc,argv,opt_sng,opt_long,&opt_idx);
    
    /* Is it time to parse the positional arguments yet? */
    if(opt == EOF) break;
    
    switch(opt){
    case 0:
      (void)fprintf(stdout,"option %s",opt_long[opt_idx].name);
      if(optarg) (void)fprintf(stdout," with arg %s",optarg);
      (void)fprintf(stdout,"\n");

      if(!strcmp(opt_long[opt_idx].name,"debug")){
	if(optarg) dbg_lvl=(unsigned short int)atoi(optarg);
	 (void)fprintf(stdout,"dbg_lvl = %d\n",dbg_lvl);
       } /* end if */
      break;
    case 'D': /* The debugging level.  Default is 0. */
      if(optarg) dbg_lvl=(unsigned short int)atoi(optarg); else dbg_lvl=1;
      (void)fprintf(stdout,"dbg_lvl = %d\n",dbg_lvl);
      break;
    case 'b': /* Set the generic tuning parameter.  Default is 0. */
      bff_sz=atol(optarg);
      break;
    case 'i': /* Get the input file name. Default is stdout */
      fl_in=optarg;
      break;
    case 'o': /* Get the output file name. Default is stdout */
      fl_out=optarg;
      break;
    case 'v': /* Print the RCS program info */
      (void)fprintf(stderr,rcs_Id);
      (void)fprintf(stderr,rcs_Revision);
      (void)fprintf(stderr,"$Author: zender $\n");
      (void)fprintf(stderr,"$Date$\n");
      (void)fprintf(stderr,"$Locker:  $\n");
      (void)fprintf(stderr,"$RCSfile: rvr.c,v $\n");
      (void)fprintf(stderr,"$Source: /home/zender/cvs/c/rvr.c,v $\n");
      (void)fprintf(stderr,"$Id$\n");
      (void)fprintf(stderr,"$State: Exp $\n");
      exit(EXIT_SUCCESS);
      break;
    default: /* Print proper usage */
      (void)usg_prn(opt_sng);
      exit(EXIT_FAILURE);
    } /* end switch */
  } /* end while loop */
  
  /* Process the positional arguments */
  psn_arg_nbr=argc-optind;
  if(dbg_lvl > 1){
    int tmp_idx=optind;

    (void)fprintf(stdout,"There are %d non-option ARGV-elements:\n",psn_arg_nbr);
    while(optind < argc) (void)fprintf(stdout,"%s ",argv[tmp_idx++]);
    (void)fprintf(stdout,"\n");
  } /* endif dbg */

  if(psn_arg_nbr > 2){
    (void)fprintf(stdout,"%s: ERROR Too many (%d) non-option ARGV-elements, exiting...\n",prg_nm,psn_arg_nbr);
    exit(EXIT_FAILURE);
  } /* endif */
  
  if(psn_arg_nbr > 0)
    if(fl_out == NULL)
      fl_out=argv[argc-1];
  
  if(psn_arg_nbr == 2)
    if(fl_in == NULL)
      fl_in=argv[argc-2];
  
  fp_in=fopen(fl_in,"r"); 
  if(fp_in == NULL){
    (void)fprintf(stderr,"%s: ERROR Could not open input file %s\n",prg_nm,fl_in);
      exit(EXIT_FAILURE);
  } /* endif */

  fp_out=fopen(fl_out,"w"); 
  if(fp_out == NULL){
    (void)fprintf(stderr,"%s: ERROR Could not open output file %s\n",prg_nm,fl_out);
      exit(EXIT_FAILURE);
  } /* endif */

  /* If buffer size was not specified on command line, set it to the size of the input file */
  if(bff_sz <= 0){
    rcd=stat(fl_in,&stat_sct);
    bff_sz=stat_sct.st_size;
  } /* endif */

  if(dbg_lvl > 0){
    (void)fprintf(stderr,"%s: fl_in = %s, fl_out = %s, bff_sz = %ld\n",prg_nm,fl_in,fl_out,bff_sz);
  } /* end if dbg */

  bff_in=(char *)malloc(bff_sz);
  bff_out=(char *)malloc(bff_sz);
  rcd=fread((void *)bff_in,1,bff_sz,fp_in);
  (void)fclose(fp_in); 
  for(idx=0;idx<bff_sz;idx++){
    bff_out[bff_sz-idx-1]=bff_in[idx];
  } /* end loop over bff */
  rcd=fwrite((void *)bff_out,1,bff_sz,fp_out);
  (void)fclose(fp_out); 

  if(rcd != 1){
    /* This condition may not actually be a fatal error. Linux fwrite returns nitems, but
       Solaris appears to return nitems*bff_sz */
    (void)fprintf(stderr,"%s: ERROR: Attempt to write 1 item(s) of size %ld failed, only %d items written.\n",prg_nm,bff_sz,rcd);
    exit(EXIT_FAILURE);
  } /* endif */

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
  (void)fprintf(stderr,"-b bff_sz Manually set the file size (not recommended). Default is 0.\n");
  (void)fprintf(stderr,"-v print the RCS program version\n");
  (void)fprintf(stderr,"\n");
} /* end usg_prn() */

void Exit_gracefully(void)
{
  char *time_buf_finish;
  time_t clock;

  /* end the clock */
  
  clock=time((time_t *)NULL);
  time_buf_finish=ctime(&clock);
  (void)fprintf(stderr,"\tfinish = %s\n",time_buf_finish);

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





