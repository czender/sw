static char rcs_Id[] = "$Id$\n";
static char rcs_Revision[] = "$Revision$\n";

/* $Author: zender $
 * $Date$
 * $Locker:  $
 * $RCSfile: xcld.c,v $
 * $Source: /home/zender/cvs/cld/xcld.c,v $
 * $Id$
 * $State: Exp $
 * */

/* Purpose: Provide an easy windowing interface to the cloud program. */ 

/* $Log: not supported by cvs2svn $
/* Revision 1.1.1.1  1998-09-15 02:06:43  zender
/* Imported sources
/*
 * Revision 0.5  1995/05/24  16:51:38  zender
 * synchronization checkin before updating to heterogeneous cpu
 * environment especially solaris.
 *
 * Revision 0.4  1993/04/05  03:03:06  zender
 * a pretty good template Motif version.  all previous Athena
 * widget are stripped out.  Accelerators and icon functionality
 * is good.  Needs more bells and whistles for running "cloud".
 *
 * Revision 0.3  1993/03/19  21:28:03  zender
 * going to try to convert to motif
 *
 * Revision 0.2  1993/02/16  06:02:04  zender
 * this version runs cloud, but gets compiler errors with the
 * istat and qstat casting in the execute() function.
 *
 * Revision 0.1  1993/02/16  04:28:18  zender
 * using Athena widgets to drive cloud program. initial version
 * just opens some windows, Exit_gracefully isn't working as a
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

/* standard header files */
#include <stdio.h>              /* stderr, FILE, NULL, etc. */
#include <time.h>               /* machine time */
#include <string.h>             /* strcmp. . . */
#include <math.h>               /* sin cos cos sin 3.14159 */
#include <stdlib.h>             /* atof, atoi, malloc, getopt */ 
#include <unistd.h>             /* all sorts of POSIX stuff */ 
#include <netcdf.h>             /* netCDF routines */
#if ( defined Solaris )
#include "/opt/SUNWspro/SC3.0/include/cc/sunmath.h"            /* IEEE signal handling */
#include <siginfo.h>            /* IEEE stuff, je ne sais pas */ 
#include <ucontext.h>           /* IEEE stuff, je ne sais pas */ 
/* #include <ieeefp.h>*/
#endif
#if ( defined SunOS )
#include <signal.h>            /* IEEE stuff, je ne sais pas */ 
#endif

/* Motif widget include files */
#include <Xm/Xm.h>                 /* Standard Motif definitions */
#include <Xm/CascadeB.h>           /* CascadeButton (for menubar labels) */
#include <Xm/Command.h>            /*  */
#include <Xm/FileSB.h>             /*  */ 
#include <Xm/Form.h>               /* Constraint widget */
#include <Xm/Frame.h>              /* Frame (simulated custom widget) */
#include <Xm/Label.h>              /* Label widget */
#include <Xm/MainW.h>              /* MainWindow */
#include <Xm/MessageB.h>           /* MessageBox dialog (for help) */
#include <Xm/PushB.h>              /* Motif PushButton Widget */
#include <Xm/RowColumn.h>          /* for MenuBar (actually a RowColumn) */

/*Display *display;*/

int debug=0; /* Option D */
int debug_value=0; /* Option d */

Widget command;
Widget pshell;
Widget toplevel;

char current_netCDF_file[1024]="cloud.nc";

#define MAX_NUM_XT_ARGS 10
#define HELP_MSG \
"Use the FileSelection dialog to find netCDF files to\n\
display in the scrolling area in the main window.  Use\n\
the edit menu to create new command lines."
  
int main(int argc,char **argv)
{
  void Exit_gracefully(void);
  void SetIconPixmap(void);
  void print_usage(char *);
  void string_cmdline(int,char *[],char *,int);

  void LoadNetCDFFileCB(Widget,XtPointer,XtPointer);
  void HelpCB(Widget,XtPointer,XtPointer);
  void EditMenuCB(Widget,XtPointer,XtPointer);
  void FileMenuCB(Widget,XtPointer,XtPointer);
  void CommandWindowCB(Widget,XtPointer,XtPointer);
  void QuitCB(Widget,XtPointer,XtPointer);

  Arg args[MAX_NUM_XT_ARGS];

  String string_foo;

  Widget commandwindow;
  Widget editmenu;
  Widget filemenu;
  Widget helpmenu;
  Widget mainwindow;
  Widget menubar;
  Widget textwindow;
  Widget widget;

  XmString compound_string_foo;
  XmString compound_string_foo2;
  XmString compound_string_foo3;

  XtAppContext app_context;

  Boolean STDERR;
  Boolean STDIN;
  Boolean STDOUT;
  Boolean VERBOSE;

  char *char_ptr_foo;
  char *char_ptr_foo2;
  char *option_string;
  char *time_buf_start;

  char *help_accelerator_text="Ctrl<Key>H";
  
  char cmdline[CMDLINE_SIZE];
  char err_file[FILESIZE];
  char icon_title[TITLE_SIZE];
  char in_file[FILESIZE];
  char out_file[FILESIZE];
  char run_command[CMDLINE_SIZE];
  char titlebar[TITLE_SIZE];
  
  extern char *optarg;
  extern int optind;
  
  float float_foo;

  int int_foo;
  int num_args;
  int opt;
  
  time_t clock;

  /* set defaults */
  STDERR = True; /* Option E */
  STDIN = True; /* Option I */
  STDOUT = True; /* Option O */
  
  (void)strcpy(icon_title,"Xcloud"); /* Option i */
  (void)strcpy(run_command,"cloud -n 10 -l 80 -x 40 -R"); /* Option i */
  (void)strcpy(in_file,"stdin"); /* Option i */
  (void)strcpy(out_file,"stdout"); /* Option o */
  (void)strcpy(err_file,"stderr"); /* Option e */

  /* parse command line arguments */
  option_string="D:d:Ee:Ii:Oo:Vv";
  while((opt = getopt(argc,argv,option_string)) != EOF){
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
      /* print the RCS program info */
      (void)fprintf(stderr,rcs_Id);
      (void)fprintf(stderr,rcs_Revision);
      (void)fprintf(stderr,"$Author: zender $\n");
      (void)fprintf(stderr,"$Date$\n");
      (void)fprintf(stderr,"$Locker:  $\n");
      (void)fprintf(stderr,"$RCSfile: xcld.c,v $\n");
      (void)fprintf(stderr,"$Source: /home/zender/cvs/cld/xcld.c,v $\n");
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
  
  /* start the clock and save the command line */  
  string_cmdline( argc,argv,cmdline,CMDLINE_SIZE );
  (void)fprintf(stdout,"Command Line: %s\n",cmdline);
  clock=time((time_t *)NULL);
  time_buf_start=ctime(&clock);
  (void)fprintf(stderr,"\tstart = %s",time_buf_start);
  (void)sprintf(titlebar,"%s v. ",argv[0]);
  char_ptr_foo2=strchr(rcs_Revision,':');
  char_ptr_foo=strncat(titlebar,char_ptr_foo2+2,4);

  toplevel=XtVaAppInitialize
    (&app_context, /* Application context */
     "XCloud", /* Application class */
     NULL, 0, /* command line option list */
     &argc, argv, /* command line args */
     NULL, /* for missing app-defaults file */
     XmNx,100, /* this WM hint isn't usually obeyed */ 
     XmNy,25, /* this WM hint isn't usually obeyed */ 
     XmNiconX,635, /* this WM hint is usually obeyed */ 
     XmNiconY,0, /* this WM hint is usually obeyed */ 
     XmNtitle,titlebar, /* text in titlebar created by the window manager */ 
     XmNiconName,icon_title, /* text beneath icon */ 
/*     XmNwidth,400,*/
/*     XmNheight,100,*/
     NULL); /* terminate varargs list */

  /* set the icon */ 
  SetIconPixmap();

  /* Check to see that all arguments were processed, and if not then
     report an error and exit. */
/*  if (argc != 1)*/
/*    Syntax(app_context,argv[0]);*/

  mainwindow=XtVaCreateManagedWidget
    ("mainwindow",
     xmMainWindowWidgetClass,
     toplevel,
/* when you use a text widget as the workwindow then you can use the text
widget's scrollbars, not the ones that come automatically with a MainWindow */ 
/*     XmNscrollBarDisplayPolicy,XmAS_NEEDED,*/
/*     XmNscrollingPolicy,XmAUTOMATIC,*/
     NULL); /* terminate varargs list */
  
  /* create the menu bar */ 
  compound_string_foo=XmStringCreateSimple("File");
  compound_string_foo2=XmStringCreateSimple("Edit");
  compound_string_foo3=XmStringCreateSimple("Help");
  menubar=XmVaCreateSimpleMenuBar
    (mainwindow,
     "menubar",
     XmVaCASCADEBUTTON,compound_string_foo,'F',
     XmVaCASCADEBUTTON,compound_string_foo2,'E',
     XmVaCASCADEBUTTON,compound_string_foo3,'H',
     NULL); /* terminate varargs list */
  XmStringFree(compound_string_foo);
  XmStringFree(compound_string_foo2);
  /* don't free "help" compound string yet -- reuse it later */ 

  /* Tell the menubar which button is the help menu  */
  if (widget=XtNameToWidget(menubar,"button_2"))
    XtVaSetValues(menubar,XmNmenuHelpWidget,widget,NULL);

  /* First menu is the File menu -- callback is FileMenuCB() */
  compound_string_foo=XmStringCreateSimple("New ...");
  compound_string_foo2=XmStringCreateSimple("Quit");
  filemenu=XmVaCreateSimplePulldownMenu
    (menubar, /* widget ID of the MenuShell's parent */ 
     "filemenu", /* name of the newly created widget */ 
     0, /* parent's cascade button that posts the menu */ 
     FileMenuCB, /* procedure to call when a button is activated */ 
     XmVaPUSHBUTTON,compound_string_foo,'N',NULL,NULL,
     XmVaSEPARATOR,
     XmVaPUSHBUTTON,compound_string_foo2,'Q',NULL,NULL,
     NULL); /* terminate varargs list */
  XmStringFree(compound_string_foo);
  XmStringFree(compound_string_foo2);

  /* Second menu is the Edit menu -- callback is EditMenuCB() */
  compound_string_foo=XmStringCreateSimple("Cray YMP");
  compound_string_foo2=XmStringCreateSimple("Sun");
  editmenu=XmVaCreateSimplePulldownMenu
    (menubar, /* widget ID of the MenuShell's parent */
     "editmenu", /* name of the newly created widget */
     1, /* parent's cascade button that posts the menu */ 
     EditMenuCB, /* procedure to call when a button is activated */
     XmVaRADIOBUTTON,compound_string_foo,'C',NULL,NULL,
     XmVaRADIOBUTTON,compound_string_foo2,'S',NULL,NULL,
     XmNradioBehavior,True, /* RowColumn resources to enforce */
     XmNradioAlwaysOne,True, /* radio behavior in Menu */
     NULL); /* terminate varargs list */
  XmStringFree(compound_string_foo);
  XmStringFree(compound_string_foo2);
  
  /* Initialize menu so that "sun" is selected. */
  if(widget=XtNameToWidget(editmenu,"button_0"))
    XtVaSetValues(widget,XmNset,True,NULL);
  
  /* Third menu is the help menu -- callback is HelpCB() */
  helpmenu=XmVaCreateSimplePulldownMenu
    (menubar, /* widget ID of the MenuShell's parent */ 
     "helpmenu", /* name of the newly created widget */ 
     2, /* parent's cascade button that posts the menu */ 
     HelpCB, /* procedure to call when a button is activated */ 
     XmVaPUSHBUTTON,compound_string_foo3,'H',NULL,NULL,
     NULL); /* terminate varargs list */
  XmStringFree(compound_string_foo3); /* we're done with it; now we can free it */

  /* Install a keyboard accelerator on the help menu  */
  compound_string_foo=XmStringCreateSimple("Ctrl-h");
  if(widget=XtNameToWidget(helpmenu,"button_0"))
    XtVaSetValues(widget,
		  XmNaccelerator,help_accelerator_text,
		  XmNacceleratorText,compound_string_foo,
		  NULL);
  
  XtManageChild(menubar);
  
  /* Create ScrolledText -- this is work area for the MainWindow */
  num_args=0;
  XtSetArg(args[0], XmNrows,      24); num_args++;
  XtSetArg(args[1], XmNcolumns,   80); num_args++;
  XtSetArg(args[2], XmNeditable,  False); num_args++;
  XtSetArg(args[3], XmNeditMode,  XmMULTI_LINE_EDIT); num_args++;
  textwindow=XmCreateScrolledText(mainwindow,"textwindow",args,num_args);
  XtManageChild(textwindow);
  
  /* store textwindow as user data in "File" menu for file_cb() callback */
/*  XtVaSetValues(menu, XmNuserData, textwindow, NULL);*/
  
  /* Create the command area -- this must be a Command class widget */
  compound_string_foo=XmStringCreateSimple("Command:");
  compound_string_foo2=XmStringCreateSimple(run_command);
  commandwindow=XtVaCreateWidget
    ("commandwindow", 
     xmCommandWidgetClass, 
     mainwindow,
     XmNpromptString, compound_string_foo,
     XmNcommand, compound_string_foo2,
     NULL);
  XmStringFree(compound_string_foo);
  XmStringFree(compound_string_foo2);
  XtAddCallback(commandwindow, XmNcommandEnteredCallback, 
		CommandWindowCB, textwindow);
  XtManageChild(commandwindow);
  
  /* set the textwindow as the "work area" of the main window */
  /*  XtVaSetValues(mainwindow,XmNmenuBar,menubar,NULL);*/
  XmMainWindowSetAreas
    (mainwindow, /* widget ID of MainWindow widet */ 
     menubar, /* widget ID for MenuBar */ 
     commandwindow, /* widget ID for Command window */ 
     NULL, /* widget ID for horizontal ScrollBar */ 
     NULL, /* widget ID for vertical ScrollBar */ 
     XtParent(textwindow)); /* widget ID for the work window */ 
  
  /* Create windows for widgets and map them. */
  XtRealizeWidget(toplevel);
  
  /* Loop for events. */
  XtAppMainLoop(app_context);
} /* end main() */ 

void print_usage(char *option_string)
{
  (void)fprintf(stderr,"\nusage: cloud [-options] where options are one or more of:\n\n");
  (void)fprintf(stderr,"%s\n\n",option_string);
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
} /* end print_usage() */ 

void Exit_gracefully(void)
{
  char *time_buf_finish;
  time_t clock;
  
  /* end the clock */  
  
  clock=time((time_t *)NULL);
  time_buf_finish=ctime(&clock);
  (void)fprintf(stdout,"\tfinish = %s\n",time_buf_finish);
  
  exit(0);
} /* end Exit_gracefully() */ 

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

void string_cmdline(int argc,char *argv[],char *cmdline,int linemax )
{
  int i;

  if( argc <= 0 ) {
    cmdline[0] = '\0';
  } else {
    (void)strcpy(cmdline,argv[0]);
    for( i=1; i<argc; i++ ) {
      (void)strncat(cmdline," ",linemax);
      (void)strncat(cmdline,argv[i],linemax);
    }
  }
} /* end string_cmdline() */ 

  /* execute the command and redirect output to the ScrolledText window */
  void CommandWindowCB(Widget commandwindow,XtPointer client_data,
		       XtPointer call_data)
{
  char *cmd,buf[BUFSIZ];
  XmTextPosition pos;
  FILE *pp,*popen();
  
  Widget textwindow=(Widget)client_data; /* passed the textwindow as client_data */
  XmCommandCallbackStruct *cbs=(XmCommandCallbackStruct *)call_data;
  
  XmStringGetLtoR(cbs->value,XmSTRING_DEFAULT_CHARSET,&cmd);
  
  if(!cmd || !*cmd){ /* nothing typed? */
    if (cmd)
      XtFree(cmd);
    return;
  } /* end if */ 
  
  /* make sure the file is a regular text file and open it */
  if (!(pp=popen(cmd,"r")))
    perror(cmd);
  XtFree(cmd);
  if (!pp)
    return;
  
  /* put the output of the command in the Text widget by reading
     until EOF (meaning that the command has terminated). */
  for (pos=0; fgets(buf,sizeof buf,pp); pos += strlen(buf))
    XmTextReplace(textwindow,pos,pos,buf);
  
  pclose(pp);
} /* end CommandWindowCB() */ 

/* The help button in the help menu from the menubar was selected.
   Display help information defined above for how to use the program.
   This is done by creating a Motif information dialog box.  Again,
   make the dialog static so we can reuse it. */

void HelpCB(Widget w,XtPointer client_data,XtPointer call_data)
{
  static Widget helpdialog;
  
  if(!helpdialog){
    Arg args[1];
    XmString help_msg=XmStringCreateLtoR(HELP_MSG,XmSTRING_DEFAULT_CHARSET);
    XtSetArg(args[0],XmNmessageString,help_msg);
    helpdialog=XmCreateInformationDialog(toplevel,"helpdialog",args,1);
  } /* end if */ 
  XtManageChild(helpdialog);
  XtPopup(XtParent(helpdialog),XtGrabNone);
} /* end HelpCB() */ 

void EditMenuCB(Widget w,XtPointer client_data,XtPointer call_data)
{
  ;
} /* end EditMenuCB() */ 

/* Any item the user selects from the File menu calls this function.
   It will either be "New" (item_no == 0) or "Quit" (item_no == 1). */
void FileMenuCB(Widget w,XtPointer client_data,XtPointer call_data)
{
    static Widget fileselect; /* make it static for reuse */
    extern void LoadNetCDFFileCB();
    extern void Exit_gracefully();

    int item_no;  /* the index into the menu */

    item_no=(int)client_data;

    if(item_no==1){
      Exit_gracefully(); /* the "quit" item */
    } /* end if */

    /* "New" was selected.  Create a Motif FileSelectionDialog w/callback,
       but make sure it doesn't already exist first! */
    if(!fileselect){
      fileselect=XmCreateFileSelectionDialog(toplevel,"fileselect",NULL,0);
      XtAddCallback(fileselect,XmNokCallback,LoadNetCDFFileCB,NULL);
      XtAddCallback(fileselect,XmNcancelCallback,
		    (XtCallbackProc)XtUnmanageChild,NULL);
    } /* end if */ 

    XtManageChild(fileselect);
    XtPopup(XtParent(fileselect),XtGrabNone);
} /* end EditMenuCB() */ 

/* The Ok button was selected from the FileSelectionDialog (or, the user
   double-clicked on a file selection).  Try to read the file as a netCDF
   data file. */
void LoadNetCDFFileCB(Widget fileselect,XtPointer client_data,XtPointer call_data)
{
  XmFileSelectionBoxCallbackStruct *cbs; /* NULL if called from change_color() */

  char *file=NULL;

  /* cbs is NULL if called from (dynapix) change_color() */
  cbs=(XmFileSelectionBoxCallbackStruct *)call_data; 
  
  if(cbs){
    if (!XmStringGetLtoR(cbs->value,XmSTRING_DEFAULT_CHARSET,&file))
      return; /* internal error */
    (void)strcpy(current_netCDF_file,file);
    XtFree(file); /* free allocated data from XmStringGetLtoR() */
  } /* end if */
} /* end LoadNetCDFFileCB() */ 

void QuitCB(Widget w,XtPointer client_data,XtPointer call_data)
{ 
  void (*exit_function)(void);

  exit_function=client_data;
  (*exit_function)();
} /* end QuitCB() */ 

void SetIconPixmap()
{
  Screen *screen;
  Pixmap pixmap;

  screen = XtScreen(toplevel);
  pixmap = XmGetPixmap
    (screen, /* screen on which pixmap will be drawn */ 
     "xcloudicon.h", /* string name of the image */ 
     BlackPixelOfScreen(screen), /* foreground pixel to combine with image */ 
     WhitePixelOfScreen(screen)); /* background pixel to combine with image */ 
  
  XtVaSetValues(toplevel,
		XmNiconPixmap, pixmap,
		NULL);
} /* end SetIconPixmap() */ 

