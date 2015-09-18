static char rcs_Id[] = "$Id$\n";
static char rcs_Revision[] = "$Revision$\n";

/* $Author: zender $
 * $Date$
 * $Locker:  $
 * $RCSfile: pvm_test.c,v $
 * $Source: /home/zender/cvs/c/pvm_test.c,v $
 * $Id$
 * $State: Exp $
 * */

/* Purpose:  Template for PVM. */

/* Notes: 

   PVM log file from remote processes may be browsed with
   less /tmp/pvml.3555

   PVM error codes may be browsed with
   less /opt/pvm3/include/pvm3.h
   less /usr/local/pvm3/include/pvm3.h
   less $PVM_ROOT/include/pvm3.h
   
   I use some PVM codes for debugging purposes so here is 
   a list of the ones i might use from PVM 3.3.2:
   Option value        MEANING
   PvmTaskDefault 0    PVM can choose any machine to start task
   PvmTaskHost    1    where specifies a particular host
   PvmTaskArch    2    where specifies a type of architecture
   PvmTaskDebug   4    Start up processes under debugger
   PvmTaskTrace   8    Processes will generate PVM trace data. *
   PvmMppFront    16   Start process on MPP front-end.
   PvmHostCompl   32   Use complement host set
   */

/* I also define my own debug values: */
#define foo_debug 1

/* Example Usage (place mousable command lines here):
   run 2 tasks, one each on the host and eagle7, in debug mode
   pvm_test -t 2 -D 5 -s eagle7
   run 2 tasks, both on the host, in debug mode
   pvm_test -t 2 -D 4 -s hostname
   run 2 tasks, on any machine, in debug mode
   pvm_test -t 2 -D 4
   run 5 tasks, on any machine, in production mode
   pvm_test -t 5
   */

/* $Log: not supported by cvs2svn $
/* Revision 1.2  2000-09-13 18:57:55  zender
/* Fixed Makefile to use .PHONY on binaries
/*
/* Revision 1.1.1.1  1998/08/31 01:25:20  zender
/* Imported sources
/*
 * Revision 1.2  1994/09/14  20:46:17  zender
 * this version works well under SunOS and AIX. about to add
 * pvm_config() and improve task distribution.
 *
 * Revision 1.1  1994/09/13  08:30:25  zender
 * Initial revision
 *
 * */

#define Boolean int
#define True 1
#define False 0
#define YES 1
#define NO 0
#define FILE_NAME_SIZE 80
#define CMD_LINE_SIZE 200
#define MACH_NAME_SIZE 30

/* Standard header files */
#include <stdio.h>              /* stderr, FILE, NULL, etc. */
#include <time.h>               /* machine time */
#include <string.h>             /* strcmp. . . */
#include <math.h>               /* sin cos cos sin 3.14159 */
#include <stdlib.h>             /* atof, atoi, malloc, getopt */
#include <unistd.h>             /* all sorts of POSIX stuff */

#include <pvm3.h>
#include <pvmsdpro.h>
#include <pvmtev.h>

#if defined (AIX)
#include "/usr/local/pvm3/pvmgs/pvmgs.h"
#include "/usr/local/pvm3/pvmgs/pvmgdef.h"
#endif

#if defined (Solaris)
#include "/opt/SUNWspro/SC3.0/include/cc/sunmath.h"            /* IEEE signal handling */
#include <siginfo.h>            /* IEEE stuff, je ne sais pas */
#include <ucontext.h>           /* IEEE stuff, je ne sais pas */
#endif

#if defined (SunOS)
#include <signal.h>            /* IEEE stuff, je ne sais pas */
#endif

/* Global variables declared here */
FILE *fp_err,*fp_in,*fp_out;

unsigned short int debug=0; /* Option D */
unsigned short int debug_value=0; /* Option d */

int main(int argc,char **argv)
{
  extern char *ctime(const time_t *);  
  extern time_t time(time_t *);
  extern void exit(int);

  int gethostname
    (char *,
     int);
  int pvm_config
    (int *,
     int *,
     struct pvmhostinfo **);

  void Exit_gracefully(void);
#if defined (SunOS)
  void ieee_exception_handler
    (int,
     int,
     struct sigcontext,
     char *);
#endif
#if defined (Solaris)
  void ieee_exception_handler
    (int,
     siginfo_t *, 
     ucontext_t *);
#endif
  void print_usage(char *);
  void string_cmd_line
    (int,
     char *[],
     char *,
     int);

  Boolean STDERR;
  Boolean STDIN;
  Boolean STDOUT;
  Boolean VERBOSE;

  char *time_buf_start;
  char *option_string;

  char cmd_line[CMD_LINE_SIZE];
  char out_file[FILE_NAME_SIZE];
  char in_file[FILE_NAME_SIZE];
  char err_file[FILE_NAME_SIZE];
  char slave_name[FILE_NAME_SIZE];
  
  extern char *optarg;
  extern int optind;
  
  float float_foo;
  float float_foo_input;

  int int_foo;
  int opt;
  
  time_t clock;

  /* PVM names and IDentifiers  */

  struct pvmhostinfo{
    int  task_id;
    char *name;
    char *arch;
    int  speed;
  };
  struct pvmhostinfo *host;

  struct pvmtaskinfo{
    int task_id;
    int parent_task_id;
    int pvmd_host_task_id;
    int task_status;
    char *task_name;
  }; 
  struct pvmtaskinfo *task;

  int *task_id;
  int left_task_id;
  int right_task_id;
  int msg_len;
  int msg_tag;
  int sender_task_id;
  int pvm_spawn_code;

  /* Communication  data buffer layouts */
  double *data;
	
  /* "The Great Unwashed" */
  char **mach_name;

  int buffer_size;
  int idx;
  int loop_idx;
  int mach_idx;
  int my_task_id;
  int num_loop;
  int num_pvmd;
  int num_arch;
  int num_mach;
  int num_task;
  int parent_task_id;
  int rcode;
  int task_idx;

  /* set defaults */
  STDERR = True; /* Option E */
  STDIN = True; /* Option I */
  STDOUT = True; /* Option O */
  
  buffer_size = 2000; /* Option b */
  float_foo_input = 0.; /* Option f */
  num_loop = 4; /* Option l */
  num_mach = 2; /* Option m */
  num_task = 2; /* Option t */

  (void)strcpy(in_file,"stdin"); /* Option i */
  (void)strcpy(out_file,"stdout"); /* Option o */
  (void)strcpy(err_file,"stderr"); /* Option e */
#if defined (AIX)
  (void)strcpy(slave_name,"eagle7"); /* Option s */
#endif
#if defined (SunOS)
  (void)strcpy(slave_name,"francke"); /* Option s */
#endif
  
  /* parse command line arguments */
  option_string="b:D:d:Ee:f:Ii:l:m:Oo:s:t:Vv";
  while((opt=getopt(argc,argv,option_string)) != EOF){
    switch(opt){
    case 'b':
      /* Set the data buffer size.  Default is 2000 */
      buffer_size=atoi(optarg);
      break;
    case 'D':
      /* The debugging level.  Default is 0. */
      debug=(unsigned short int)atoi(optarg);
      break;
    case 'd':
      /* The second debugging level.  Default is 0. */
      debug_value=(unsigned short int)atoi(optarg);
      break;
    case 'E':
      /* Toggle the error file stream. Default is True */
      STDERR=!STDERR;
      break;
    case 'e':
      /* get the error file name. Default is stderr */
      (void)strcpy(err_file,optarg);
      break;
    case 'f':
      /* Set the generic tuning parameter.  Default is 0. */
      float_foo_input=atof(optarg);
      break;
    case 'I':
      /* Toggle the input file stream. Default is True */
      STDIN=!STDIN;
      break;
    case 'i':
      /* get the input file name. Default is stdin */
      (void)strcpy(in_file,optarg);
      break;
    case 'l':
      /* Set the number of loops.  Default is 4 */
      num_loop=atoi(optarg);
      break;
    case 'm':
      /* Set the number of machines.  Default is 1 */
      num_mach=atoi(optarg);
      break;
    case 'O':
      /* Toggle the output file stream. Default is True */
      STDOUT=!STDOUT;
      break;
    case 'o':
      /* get the output file name. Default is stdout */
      (void)strcpy(out_file,optarg);
      break;
    case 's':
      /* get the name of the slave machine. Default is starfury/wildhorse */
      (void)strcpy(slave_name,optarg);
      break;
    case 't':
      /* Set the number of tasks.  Default is 2 */
      num_task=atoi(optarg);
      break;
    case 'v':
      /* print the RCS program info */
      (void)fprintf(stderr,rcs_Id);
      (void)fprintf(stderr,rcs_Revision);
      (void)fprintf(stderr,"$Author: zender $\n");
      (void)fprintf(stderr,"$Date$\n");
      (void)fprintf(stderr,"$Locker:  $\n");
      (void)fprintf(stderr,"$RCSfile: pvm_test.c,v $\n");
      (void)fprintf(stderr,"$Source: /home/zender/cvs/c/pvm_test.c,v $\n");
      (void)fprintf(stderr,"$Id$\n");
      (void)fprintf(stderr,"$State: Exp $\n");
      exit(0);
      break;
    case 'V':
      /* toggle verbose printing out of WARNINGS. Default is True */
      VERBOSE=!VERBOSE;
      break;
    case '?':
      /* print proper usage */
      (void)print_usage(option_string);
      exit(1);
    } /* end switch */
  } /* end while loop */
  
  if(STDERR){
    fp_err = stderr;
  }else{
    if( (fp_err = fopen( err_file, "w")) == NULL) {
      (void)fprintf(stderr,"\nError in opening error file %s\n",err_file);
      exit(1);
    } /* end if */
  } /* end else */
  
  if(STDIN){
    fp_in = stdin;
  }else{
    if( (fp_in = fopen( in_file, "r")) == NULL) {
      (void)fprintf(stderr,"\nError in opening input file %s\n",in_file);
      exit(1);
    } /* end if */
  } /* end else */
  
  if(STDOUT){
    fp_out = stdout;
  }else{
    if( (fp_out = fopen( out_file, "w")) == NULL) {
      (void)fprintf(stderr,"\nError in opening output file %s\n",out_file);
      exit(1);
    } /* end if */
  } /* end else */
  
  /* start the clock and save the command line */
  string_cmd_line(argc,argv,cmd_line,CMD_LINE_SIZE);
  (void)fprintf(fp_err,"Command Line: %s\n",cmd_line);
  clock=time((time_t *)NULL);
  time_buf_start=ctime(&clock);
  (void)fprintf(fp_err,"\tstart = %s",time_buf_start);

  /* Attempt to set up IEEE signal handling */
#if defined (sun)
  if(ieee_handler("set","common",(sigfpe_handler_type)ieee_exception_handler)!=0){
    (void)fprintf(fp_err,"IEEE trapping not supported here.\n");
  } /* end if */
#endif

  /* Main body of code */
  if(
     ((task_id=(int *)malloc((num_task)*sizeof(int))) == NULL ) ||
     ((data=(double *)malloc((buffer_size)*sizeof(double))) == NULL ) ||     
     ((mach_name=(char **)malloc((num_task)*sizeof(char *))) == NULL ) ||
     False ){
    (void)fprintf(fp_err,"Unable to allocate array in main\n");
    exit(1);
  } /* end if */
  
  for(task_idx=0;task_idx<=num_task-1;task_idx++){
    if(
       ((mach_name[task_idx]=(char *)malloc((MACH_NAME_SIZE)*sizeof(char))) == NULL ) ||
       False ){
      (void)fprintf(fp_err,"Unable to allocate mach_name array in main\n");
      exit(1);
    } /* end if */
  } /* end for */

  my_task_id=pvm_mytid();
  parent_task_id=pvm_parent();
    
  rcode=gethostname(mach_name[0],MACH_NAME_SIZE);

  if(debug == 0){
    (void)fprintf(fp_err,"Machine: %s. My task ID: %d. ",mach_name[0],my_task_id);
    if (parent_task_id == PvmNoParent){
      (void)fprintf(fp_err,"No Parent Task.\n");
    }else{ /* end if this task was spawned */
      (void)fprintf(fp_err,"Parent Task ID: %d\n",parent_task_id);
    } /* end if this task was not spawned */
  } /* end debug */
  
  if(parent_task_id == PvmNoParent){
    task_id[0]=my_task_id;
    pvm_spawn_code=PvmTaskDefault;
    if((int)debug%(2*PvmTaskHost) >= (int)PvmTaskHost){
      pvm_spawn_code=pvm_spawn_code+PvmTaskHost;
    } /* end if */
    if((int)debug%(2*PvmTaskDebug) >= PvmTaskDebug){
      pvm_spawn_code=pvm_spawn_code+PvmTaskDebug;
    } /* end if */
    rcode=pvm_spawn
      ("pvm_test", /* name of executable */
       argv+1, /* argv for executable */
       pvm_spawn_code, /* pvm_spawn_code */
       slave_name, /* machine or architecture */
       num_task-1, /* # of copies of executable to start */
       task_id+1); /* space for returned task id's */
    if(rcode == num_task-1){
      (void)fprintf(fp_err,"%d tasks started by pvm_spawn():\n",rcode);
      /* get the environment */
      rcode=pvm_config
	(&num_pvmd,
	 &num_arch,
	 (struct pvmhostinfo **)&host);
      if(rcode != 0){
	(void)fprintf(fp_err,"Error returning from pvm_config()\n");
      } /* end if */
      (void)fprintf(fp_err,"Returning from pvm_config(). # pvmds = %d, # data formats is %d.\n",num_pvmd,num_arch);
      (void)fprintf(fp_err,"%10s%10s%30s%10s%10s\n",
		    "Ord. idx","Host ID","Host","Arch","Speed");
      for(task_idx=0;task_idx<=num_pvmd-1;task_idx++){
	(void)fprintf
	  (fp_err,"%10d%10d%30s%10s%10d\n",
	   task_idx,
	   host[task_idx].task_id,
	   host[task_idx].name,
	   host[task_idx].arch,
	   host[task_idx].speed
	   );
      } /* end loop over task */
      rcode=pvm_tasks
	(0, /* 0 means entire virtual machine */
	 &int_foo, /* # of tasks being reported on */
	 (struct pvmtaskinfo **)&task);
      if(int_foo != num_task){
	(void)fprintf(fp_err,"Discrepency in # of tasks: pvm_tasks() says there are %d, i wanted %d. Trying to kill old tasks with.\n",int_foo,num_task);
      } /* end if */
      if(rcode != 0){
	(void)fprintf(fp_err,"Error returning from pvm_tasks()\n");
      } /* end if */
      (void)fprintf(fp_err,"Returning from pvm_tasks(). # tasks = %d\n",int_foo);
      (void)fprintf(fp_err,"%10s%10s%10s%10s%10s%20s\n",
		    "Ord. idx","Task ID","Parent ID","Host ID","Status","Task Name");
      for(task_idx=0;task_idx<=int_foo-1;task_idx++){
/*	(void)fprintf
	  (fp_err,"%10d%10d%10d%10d%10d%20s\n",
	   task_idx,
	   task[task_idx].task_id,
	   task[task_idx].parent_task_id,
	   task[task_idx].pvmd_host_task_id,
	   task[task_idx].task_status,
	   task[task_idx].task_name
	   ); */
      } /* end loop over task */
      for(task_idx=0;task_idx<=num_task-1;task_idx++){
	(void)fprintf(fp_err,"\tTask # %d ID is %d.\n",task_idx,task_id[task_idx]);
      } /* end loop over task */
    }else{ /* end if pvm_spawn() was successful */
      if(rcode < 0){
	(void)fprintf(fp_err,"System error. Attempted to spawn %d tasks. Error code is %d. Shutting down.\n",num_task-1,rcode);
	pvm_exit();
	exit(1);
      }else{ /* end if pvm_spawn() caused system error */
	(void)fprintf(fp_err,"Attempted to spawn %d tasks, only spawned %d.\n",num_task-1,rcode);
	for(task_idx=rcode+1;task_idx<=num_task-1;task_idx++){
	  (void)fprintf(fp_err,"Error code for ordinal task %d is %d.\n",task_idx,task_id[task_idx]);
	} /* end loop over task */
	(void)fprintf(fp_err,"Shutting down.\n");
      } /* endif */
      for(task_idx=1;task_idx<=rcode;task_idx++){
	if(task_id[task_idx]!=my_task_id){
	  pvm_kill(task_id[task_idx]);
	} /* endif */
      } /* end loop over task */
      pvm_exit();
      exit(1);
    } /* end if pvm_spawn() was unsuccessful */
  }else{ /* end if this is the parent task */
    task_id[1]=my_task_id;
    task_id[0]=pvm_parent();
  } /* end if this is a child task */

  for(idx=0;idx<buffer_size;idx++){
    data[idx]=sin((float)idx) ;
  } /* endfor */

  for(loop_idx=1;loop_idx<=num_loop;loop_idx++){

    /* Sending to the "left" and "right" comes from the token ring paradigm */
    if(my_task_id == task_id[0]){
      left_task_id=right_task_id=task_id[1];
    } /* endif */
    if(my_task_id == task_id[1]){
      left_task_id=right_task_id=task_id[0];
    } /* endif */

    if(debug >= 0){
      (void)fprintf(fp_err,"%s (%d): Beginning loop # %d. Sending left to %d and right to %d.\n",mach_name[0],my_task_id,loop_idx,left_task_id,right_task_id);
    } /* end debug */
    
    if(
       pvm_psend
       (left_task_id, /* destination process */
	loop_idx, /* message tag */
	data, /* buffer */
	buffer_size, /* length of buffer in multiple of data type size */
	PVM_DOUBLE) /* type of data in buffer */
       < 0){
      /* An error occurred */
      (void)fprintf(fp_err,"Send error in pvm_psend()\n");
      for(task_idx=0;task_idx<num_task;task_idx++){
	if(task_id[task_idx]!=my_task_id){
	  pvm_kill(task_id[task_idx]);
	} /* endif */
      } /* end loop over task */
      pvm_exit();
      exit(1);
    } /* endif pvm_psend() didn't work*/
    
    if(
       pvm_psend
       (right_task_id, /* destination process */
	loop_idx, /* message tag */
	data, /* buffer */
	buffer_size, /* length of buffer in multiple of data type size */
	PVM_DOUBLE) /* type of data in buffer */
       < 0){
      /* An error occurred */
      (void)fprintf(fp_err,"Send error in pvm_psend()\n");
      for(task_idx=0;task_idx<num_task;task_idx++){
	if(task_id[task_idx]!=my_task_id){
	  pvm_kill(task_id[task_idx]);
	} /* endif */
      } /* end loop over task */
      pvm_exit();
      exit(1);
    } /* endif pvm_psend() didn't work*/
    
    if(debug >= 0){
      (void)fprintf(fp_err,"%s (%d): Sent messages to left and right. Waiting for recvs.\n",mach_name[0],my_task_id);
    } /* end debug */
    
    if(
       (rcode=
	pvm_precv
	(left_task_id, /* task id of sending process (to match) */
	 loop_idx, /* message tag (to match) */
	 data, /* buffer to receive data */
	 buffer_size, /* length of buffer */
	 PVM_DOUBLE, /* data type of buffer */
	 &sender_task_id, /* actual task id of sender */
	 &msg_tag, /* actual message tag */
	 &msg_len) /* actual message length */
	)<0){
      /* An error occurred */
      (void)fprintf(fp_err,"Receive error on pvm_precv()\n");
      for(task_idx=0;task_idx<num_task;task_idx++){
	if(task_id[task_idx]!=my_task_id){
	  pvm_kill(task_id[task_idx]);
	} /* endif */
      } /* end loop over task */
      pvm_exit();
      exit(1);
    }else{ /* endif receive failed */
      (void)fprintf(fp_err,"%s (%d): Received data of length %d and tag %d from %d.\n",mach_name[0],my_task_id,msg_len,msg_tag,sender_task_id);
    } /* end if receive succeeded */
    
    if(
       (rcode=
	pvm_precv
	(right_task_id, /* task id of sending process (to match) */
	 loop_idx, /* message tag (to match) */
	 data, /* buffer to receive data */
	 buffer_size, /* length of buffer */
	 PVM_DOUBLE, /* data type of buffer */
	 &sender_task_id, /* actual task id of sender */
	 &msg_tag, /* actual message tag */
	 &msg_len) /* actual message length */
	)<0){
      /* An error occurred */
      (void)fprintf(fp_err,"Receive error on pvm_precv()\n");
      for(task_idx=0;task_idx<num_task;task_idx++){
	if(task_id[task_idx]!=my_task_id){
	  pvm_kill(task_id[task_idx]);
	} /* endif */
      } /* end loop over task */
      pvm_exit(); 
      exit(1);
    }else{ /* endif receive failed */
      (void)fprintf(fp_err,"%s (%d): Received data of length %d and tag %d from %d.\n",mach_name[0],my_task_id,msg_len,msg_tag,sender_task_id);
    } /* end if receive succeeded */
    
    if(debug >= 0){
      (void)fprintf(fp_err,"%s (%d): Ending loop # %d.\n\n",mach_name[0],my_task_id,loop_idx);
    } /* end debug */
  } /* end loop over iterations */
  
  /* Leave PVM */
  pvm_exit();
  
  /* Find out if there are any outstanding exceptions */
#if defined (sun)
  ieee_retrospective(fp_err);
#endif
  
  /* Clear all signal handling traps before exiting */
#if defined (sun)
  ieee_handler("clear","common",(sigfpe_handler_type)ieee_exception_handler);
#endif
  
  Exit_gracefully();
} /* end main() */

#if defined (SunOS)
/* See the Sun Numerical Computation Guide p. 145 */
void ieee_exception_handler
(int signal,
 int code,
 struct sigcontext *scp, 
 char *address)
#endif
#if defined (Solaris)
void ieee_exception_handler
(int signal,
 siginfo_t *sip, 
 ucontext_t *uap)
#endif
#if defined (AIX)
void ieee_exception_handler()
#endif
{
  ;
  /* See /usr/include/signal.h for SIGFPE codes.
     This example is taken from the Sun Numerical Computation Guide, p. 146 */
/*  (void)fprintf(fp_err,"floating point exception %x at address %x\n",*/
/*                code,*/
/*                address);*/
} /* end ieee_exception_handler() */

void print_usage(char *option_string)
{
  (void)fprintf(stderr,"\nusage: clouds [-options] where options are one or more of:\n\n");
  (void)fprintf(stderr,"%s\n\n",option_string);
  (void)fprintf(stderr,"-b buffer_size Set the data buffer size. Default is 2000\n");
  (void)fprintf(stderr,"-D debug The debugging level.  Default is 0\n");
  (void)fprintf(stderr,"-E STDERR Toggle stderr stream.  Default is True\n");
  (void)fprintf(stderr,"-I STDIN Toggle stdin stream.  Default is True\n");
  (void)fprintf(stderr,"-O STDOUT Toggle stdout stream.  Default is True\n");
  (void)fprintf(stderr,"-V toggle verbose printing of WARNINGS. Default is True\n");
  (void)fprintf(stderr,"-d debug_value The second debug level.  Default is 0\n");
  (void)fprintf(stderr,"-e err_file Set the error file name. Default is stderr\n");
  (void)fprintf(stderr,"-f float_foo Set the generic float. Default is 0.\n");
  (void)fprintf(stderr,"-i in_file Set the input file name. Default is stdin\n");
  (void)fprintf(stderr,"-l num_loop Set the number of loops. Default is 4\n");
  (void)fprintf(stderr,"-m num_mach Set the number of machines. Default is 1\n");
  (void)fprintf(stderr,"-o out_file Set the output file name. Default is stdout\n");
  (void)fprintf(stderr,"-s slave_name Set the slave machine name. Default is starfury/wildhorse\n");
  (void)fprintf(stderr,"-t num_task Set the number of tasks. Default is 2\n");
  (void)fprintf(stderr,"-v print the RCS program version\n");
  (void)fprintf(stderr,"\n");
} /* end print_usage() */

void Exit_gracefully(void)
{
  extern char *ctime(const time_t *);  
  extern time_t time(time_t *);

  char *time_buf_finish;
  time_t clock;

  /* end the clock */
  
  clock=time((time_t *)NULL);
  time_buf_finish=ctime(&clock);
  (void)fprintf(fp_err,"\tfinish = %s\n",time_buf_finish);

  (void)fclose(fp_err);
  (void)fclose(fp_in);
  (void)fclose(fp_out);

  exit(0);
} /* end Exit_gracefully() */

/* Module:	cmdparse.c (Command Parse)
 * Purpose:	Set options from command line
 * Subroutine:	parse_cmd_line()			returns: int
 * Subroutine:	string_cmd_line()		returns: void
 * Subroutine:	usage()				returns: int
 * Xlib calls:	none
 * Copyright:	1989 Smithsonian Astrophysical Observatory
 *		You may do anything you like with this file except remove
 *		this copyright.  The Smithsonian Astrophysical Observatory
 *		makes no representations about the suitability of this
 *		software for any purpose.  It is provided "as is" without
 *		express or implied warranty.
 * Modified:	{0} Michael VanHilst	initial version	       5 January 1989
 *              {1} MVH BSDonly strings.h compatability           19 Feb 1990
 *              {2} CSZ ANSI C compatability                      29 May 1993
 *		{n} <who> -- <does what> -- <when>
 */

void string_cmd_line(int argc,char *argv[],char *cmd_line,int linemax)
{
  int current_arg;

  if(argc <= 0){
    cmd_line[0]='\0';
  }else{
    (void)strcpy(cmd_line, argv[0]);
    for(current_arg=1;current_arg<argc;current_arg++){
      (void)strncat(cmd_line," ",linemax);
      (void)strncat(cmd_line,argv[current_arg],linemax);
    } /* end loop over args */
  } /* end else */
} /* end string_cmd_line() */
