static char rcsid[] = "$Id$\n";
static char rcsrev[] = "$Revision$";

/* $Log: not supported by cvs2svn $
/* Revision 1.1.1.1  1998-09-15 02:06:44  zender
/* Imported sources
/*
 * Revision 0.3  1993/03/19  21:28:03  zender
 * going to try to convert to motif
 *
 * Revision 0.2  1993/02/16  06:02:04  zender
 * this version runs cloud, but gets compiler errors with the
 * istat and qstat casting in the execute() function.
 *
 * Revision 0.1  1993/02/16  04:28:18  zender
 * using Athena widgets to drive cloud program. initial version
 * just opens some windows, exit_gracefully isn't working as a
 * callback function yet.
 * */ 

/*#define Boolean int*/ /* this crashes with Intrinsic.h line 146 */ 
#define True 1
#define False 0
#define YES 1
#define NO 0
#define FILESIZE 80
#define TITLE_SIZE 80
#define CMDLINE_SIZE 80

/* UNIX standard include files */
#include <stdio.h>              /* stderr, FILE, NULL, etc. */
#include <time.h>               /* machine time */
#include <string.h>             /* strcmp. . . */
#include <math.h>               /* sin cos cos sin 3.14159 */
#include <vfork.h>              /* vfork */
#include <signal.h>             /* signal stuff */

/* Xlib include files */
/*#include <X11/Xlib.h>*/           /* X window stuff */
/*#include <X11/Xutil.h>*/          /* X window manager stuff */
/*#include <X11/Xos.h>*/            /* for portability */
/*#include <X11/Xatom.h>*/          /* ??? */

/* Xt include files */
#include <X11/Intrinsic.h>     /* Intrinsics Definitions*/
#include <X11/StringDefs.h>    /* Standard Name-String definitions*/
#include <X11/Shell.h>         /* Shell widgets */

/* Athena widget include files */
#include <X11/Xaw/Box.h>   /* Athena Label Widget */
#include <X11/Xaw/Command.h>   /* Athena Label Widget */
#include <X11/Xaw/Dialog.h>    /* Athena Label Widget */
#include <X11/Xaw/Form.h>      /* Athena Label Widget */
#include <X11/Xaw/Label.h>     /* Athena Label Widget */

/*Display *display;*/

int debug=0; /* Option D */
int debug_value=0; /* Option d */

Widget command;
Widget pshell;

void PressMeCB(w, client_data, call_data)
     Widget w;
     XtPointer client_data, call_data;
{ 
  fprintf(stderr, "Thankyou!\n"); 
}

void PopupDialogCB(w, client_data, call_data)
     Widget w;
     XtPointer client_data, call_data;
{ 
  Dimension height;
  Dimension width;

  Position x;
  Position y;

  String string;

  Widget topLevel = (Widget) client_data;

  /* get the coordinates of the middle of topLevel widget. */ 
  XtVaGetValues(topLevel, /* widget */ 
		XtNwidth, &width,
		XtNheight, &height,
		NULL); /* terminate varargs list */

  /* translate coordinates in application top-level window into 
   coordinates from root window origin. */ 
  XtTranslateCoords(topLevel, /* widget */ 
		    (Position) width/2, /* x */ 
		    (Position) height/2, /* y */ 
		    &x, &y); /* coords on root window */ 

  /* move popup shell to this position (it's not visible yet) */ 
  XtVaSetValues(pshell, /* widget */ 
		XtNx, x,
		XtNy, y,
		NULL); /* terminate varargs list */

  /* indicate to user that no other application functions are 
     valid while dialog is popped up  */ 
  XtSetSensitive(command, FALSE);

  XtPopup(pshell, XtGrabNonexclusive);
}

void DialogDoneCB(w, client_data, call_data)
     Widget w;
     XtPointer client_data, call_data;
{ 
  Widget dialog = (Widget) client_data;
  String string;

  XtPopdown(pshell);

  XtSetSensitive(command, TRUE);

  string = XawDialogGetValueString(dialog);

  (void)fprintf(stdout,"New command: %s\n",string);
}

void RunCB(w, client_data, call_data)
     Widget w;
     XtPointer client_data, call_data;
{ 
  String command[CMDLINE_SIZE];
  int i;
  int execute();
  int status;

  /* code from O'Reilly winman, Xlib Programming Manual (vol 1) p. 438 */
  /* Force child processes to disinherit the TCP file . descriptor; this helps the
     shell command (creating new xterm) forked and executed from the menu to work
     properly */
/*  if((fcntl(ConnectionNumber(display), F_SETFD, 1)) == -1){*/
/*    fprintf(stdout,"xcloud: child cannot disinherit TCP fd");*/
/*    exit(-1);*/
/*  }*/
  sprintf(command,"/cgd/home/zender/cloud/cloud -n 1 -x 40 -l 80 -R -o xcloud.dat &");
  status=execute(command);

  if(debug == 0){
    (void)fprintf(stdout,"executed %s, status is %d\n",command,status);
  } /* end debug */
}

void QuitCB(w, client_data, call_data)
     Widget w;
     XtPointer client_data, call_data;
{ 
  void (*exit_function)();
  exit_function = (void (*)())call_data;
/*  (*exit_function)();*/
  exit(0);
}

main(argc, argv)
     int argc;
     char **argv;
{
  void Exit_gracefully();
  void print_usage();
  void string_cmdline();

  String string_foo;

  Widget topLevel;
  Widget titlebar;
  Widget form;
  Widget quit;
  Widget pressme;
  Widget run;
  Widget dialog;
  Widget dialogdone;

  XtAppContext app_context;

  Boolean STDERR;
  Boolean STDIN;
  Boolean STDOUT;
  Boolean VERBOSE;

  char *char_ptr_foo;
  char *char_ptr_foo2;
  char *time_buf_start;
  
  char title[TITLE_SIZE];
  char cmdline[CMDLINE_SIZE];
  char run_command[CMDLINE_SIZE];
  char out_file[FILESIZE];
  char in_file[FILESIZE];
  char err_file[FILESIZE];
  
  extern char *optarg;
  extern int optind;
  
  float float_foo;

  int int_foo;
  int opt;
  
  time_t clock;

  /* set defaults */
  STDERR = True; /* Option E */
  STDIN = True; /* Option I */
  STDOUT = True; /* Option O */
  
  (void)strcpy(run_command,"cloud -n 10 -l 80 -x 40 -R"); /* Option i */
  string_foo=run_command;
  (void)strcpy(in_file,"stdin"); /* Option i */
  (void)strcpy(out_file,"stdout"); /* Option o */
  (void)strcpy(err_file,"stderr"); /* Option e */

  /* parse command line arguments */
  while((opt = getopt(argc, argv, "D:d:Ee:Ii:Oo:Vv")) != EOF){
    switch(opt){
    case 'D':
      /* The debugging level.  Default is 0. */
      debug = (unsigned short int)atoi(optarg);
      break;
    case 'd':
      /* The second debugging level.  Default is 0. */
      debug_value = (unsigned short int)atoi(optarg);
      break;
    case 'E':
      /* Toggle the error file stream. Default is True */
      STDERR = !STDERR;
      break;
    case 'e':
      /* get the error file name. Default is stderr */
      (void)strcpy(err_file,optarg);
      break;
    case 'I':
      /* Toggle the input file stream. Default is True */
      STDIN = !STDIN;
      break;
    case 'i':
      /* get the input file name. Default is stdin */
      (void)strcpy(in_file,optarg);
      break;
    case 'O':
      /* Toggle the output file stream. Default is True */
      STDOUT = !STDOUT;
      break;
    case 'o':
      /* get the output file name. Default is stdout */
      (void)strcpy(out_file,optarg);
      break;
    case 'v':
      /* print the RCS program version */
      (void)fprintf(stderr,rcsid);
      exit(0);
      break;
    case 'V':
      /* toggle verbose printing out of WARNINGS. Default is True */ 
      VERBOSE=!VERBOSE;
      break;
    case '?':
      /* print proper usage */
      (void)print_usage();
      exit(1);
    } /* end switch */
  } /* end while loop */
  
  /* start the clock and save the command line */  
  string_cmdline( argc, argv, cmdline, CMDLINE_SIZE );
  (void)fprintf(stdout,"Command Line: %s\n",cmdline);
  clock=time((time_t *)NULL);
  time_buf_start=ctime(&clock);
  (void)fprintf(stderr,"\tstart = %s",time_buf_start);
  (void)sprintf(title,"%s v. ",argv[0]);
  char_ptr_foo2=strchr(rcsrev,':');
  char_ptr_foo=strncat(title,char_ptr_foo2+2,4);

  topLevel = XtVaAppInitialize(
			       &app_context, /* Application context */
			       "XCloud", /* Application class */
			       NULL, 0, /* command line option list */
			       &argc, argv, /* command line args */
			       NULL, /* for missing app-defaults file */
			       XtNx,100,
			       XtNy,100,
/*			       XtNwidth,400,*/
/*			       XtNheight,100,*/
			       NULL); /* terminate varargs list */

  /* Check to see that all arguments were processed, and if not then
     report an error and exit. */
/*  if (argc != 1)*/
/*    Syntax(app_context, argv[0]);*/

  form = XtVaCreateManagedWidget(
				 "form", /* widget name */
				 formWidgetClass, /* widget class */
				 topLevel, /* parent widget*/
				 NULL); /* argument list*/
  
  titlebar = XtVaCreateManagedWidget(
				     "titlebar", /* widget name */
				     labelWidgetClass, /* widget class */
				     form, /* parent widget */
				     XtNlabel,char_ptr_foo,
				     NULL); /* terminate varargs list */
  
  quit = XtVaCreateManagedWidget(
				 "quit", /* widget name */
				 commandWidgetClass, /* widget class */
				 form, /* parent widget*/
				 XtNfromVert,titlebar,
				 XtNlabel,"Quit",
				 NULL); /* argument list*/
  
  pressme = XtVaCreateManagedWidget(
				    "pressme", /* widget name   */
				    commandWidgetClass, /* widget class */
				    form, /* parent widget*/
				    XtNfromVert,titlebar,
				    NULL); /* argument list*/

  run = XtVaCreateManagedWidget(
				"run", /* widget name   */
				commandWidgetClass, /* widget class */
				form, /* parent widget*/
				XtNfromHoriz,pressme, /* resource setting */ 
				XtNfromVert,titlebar,
				NULL); /* argument list*/

  command = XtVaCreateManagedWidget(
				    "command", /* widget name   */
				    commandWidgetClass, /* widget class */
				    form, /* parent widget*/
				    XtNfromVert,quit,
				    NULL); /* argument list*/
  
  pshell = XtVaCreatePopupShell(
				"pshell", /* widget name   */
				transientShellWidgetClass, /* widget class */
				topLevel, /* parent widget*/
				NULL); /* argument list*/
  
  dialog = XtVaCreateManagedWidget(
				   "dialog", /* widget name   */
				   dialogWidgetClass, /* widget class */
				   pshell, /* parent widget*/
				   XtNfromVert,quit,
				   XtNlabel,"Current command line:",
				   XtNvalue,string_foo, 
				   XtNwidth,400,
				   NULL); /* argument list*/

  dialogdone = XtVaCreateManagedWidget(
				       "dialog_done", /* widget name   */
				       commandWidgetClass, /* widget class */
				       dialog, /* parent widget*/
				       XtNlabel,"Done",
				       NULL); /* argument list*/
  
  XtVaSetValues(pressme,
		XtNfromHoriz,quit, /* resource setting */ 
		NULL); /* terminate varargs list */ 

  XtAddCallback(quit, XtNcallback, QuitCB, Exit_gracefully);
  XtAddCallback(pressme, XtNcallback, PressMeCB, 0);
  XtAddCallback(run, XtNcallback, RunCB, 0);
  XtAddCallback(command, XtNcallback, PopupDialogCB, topLevel);
  XtAddCallback(dialogdone, XtNcallback, DialogDoneCB, dialog);

  /* Create windows for widgets and map them. */
  XtRealizeWidget(topLevel);
  
  /* Loop for events. */
  XtAppMainLoop(app_context);
}

void print_usage()
{
  (void)fprintf(stderr,"\nusage: cloud [-options] where options are one or more of:\n\n");
  (void)fprintf(stderr,"D:d:Ee:Ii:Oo:vV\n\n");
  (void)fprintf(stderr,"-D debug The debugging level.  Default is 0\n");
  (void)fprintf(stderr,"-d debug_value The second debug level.  Default is 0\n");
  (void)fprintf(stderr,"-E STDERR Toggle stderr stream.  Default is True\n");
  (void)fprintf(stderr,"-I STDIN Toggle stdin stream.  Default is True\n");
  (void)fprintf(stderr,"-O STDOUT Toggle stdout stream.  Default is True\n");
  (void)fprintf(stderr,"-e err_file Get the error file name. Default is stderr\n");
  (void)fprintf(stderr,"-i in_file get the input file name. Default is stdin\n");
  (void)fprintf(stderr,"-o out_file get the output file name. Default is stdout\n");
  (void)fprintf(stderr,"-v print the RCS program version\n");
  (void)fprintf(stderr,"-V toggle verbose printing of WARNINGS. Default is True\n");
  (void)fprintf(stderr,"\n");
}

void Exit_gracefully()
{
  char *time_buf_finish;
  time_t clock;

  /* end the clock */  
  
  clock=time((time_t *)NULL);
  time_buf_finish=ctime(&clock);
  (void)fprintf(stdout,"\tfinish = %s\n",time_buf_finish);

  exit(0);
}

/* Module:	cmdparse.c (Command Parse)
 * Purpose:	Set options from command line
 * Subroutine:	parse_cmdline()			returns: int
 * Subroutine:	string_cmdline()		returns: void
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
 *		{n} <who> -- <does what> -- <when>
 */

void string_cmdline( argc, argv, cmdline, linemax )
     int argc;
     char *argv[];
     char *cmdline;
     int linemax;
{
  int i;

  if( argc <= 0 ) {
    cmdline[0] = '\0';
  } else {
    (void)strcpy(cmdline, argv[0]);
    for( i=1; i<argc; i++ ) {
      (void)strncat(cmdline, " ", linemax);
      (void)strncat(cmdline, argv[i], linemax);
    }
  }
}

int execute(s)
     char *s;
{
  int status, pid, w;
  register int (*istat)(),(*qstat)();
  
  /* code from O'Reilly winman, book 1 p. 438 */
  if((pid = vfork()) == 0){
    signal(SIGINT, SIG_DFL);
    signal(SIGQUIT, SIG_DFL);
    signal(SIGHUP, SIG_DFL);
    execl("/bin/sh", "sh", "-c", s, (char *)0);
    _exit(127);
  }
  istat = signal(SIGINT, SIG_IGN);
  qstat = signal(SIGQUIT, SIG_IGN);
  while((w = wait(&status)) != pid && w != -1)
    ;
  if(w == -1)
    status = -1;
  signal(SIGINT, istat);
  signal(SIGQUIT, qstat);
  
  return(status);
}
