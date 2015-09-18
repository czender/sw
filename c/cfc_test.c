#ifndef lint
  static char rcsid[] = "$Id$";
#endif

/* $Log: not supported by cvs2svn $
/* Revision 1.2  2000-09-13 18:57:54  zender
/* Fixed Makefile to use .PHONY on binaries
/*
/* Revision 1.1.1.1  1998/08/31 01:25:20  zender
/* Imported sources
/* */

#define YES 1
#define NO 0
#define FILESIZE 80

#define Boolean int
#define True 1
#define False 0

#include <stdio.h>              /* stderr, FILE, NULL, etc. */
#include <time.h>               /* machine time */
#include <string.h>             /* strcmp. . . */
#include <math.h>               /* sin cos cos sin 3.14159 */

FILE *fp_out,*fp_err,*fp_in;

int debug=0; /* Option D */

main(argc,argv)
     int argc;
     char * argv[];
{

#ifdef CRAY2
  void TEST_CF();
  void TEST_C();
  void TEST_TRANS();
#else
  void test_cf_();
  void test_c_();
  void test_trans_();
#endif

  void Exit_gracefully();
  void print_usage();
  void transpose_conc();

  Boolean STDERR;
  Boolean STDIN;
  Boolean STDOUT;
  Boolean VERBOSE;

  char *time_buf_start;
  char out_file[FILESIZE];
  char in_file[FILESIZE];
  char err_file[FILESIZE];
  
  extern char *optarg;
  extern int optind;
  
  float **matrix_c;
  float **matrix_c_trans;
  float **matrix_cf;
  float **matrix_f;

  float * matrix_c_ptr;
  float * matrix_c_trans_ptr;
  float * matrix_cf_ptr;
  float * matrix_f_ptr;

  int layer;
  int num_layers=7;
  int num_sizes=5;
  int opt;
  int size;

  time_t clock;

  /* set defaults */
  STDERR = True; /* Option E */
  STDIN = True; /* Option I */
  STDOUT = False; /* Option O */
  
  (void)strcpy(in_file,"stdin"); /* Option i */
  (void)strcpy(out_file,"cfc.out"); /* Option o */
  (void)strcpy(err_file,"stderr"); /* Option e */
  
  /* parse command line arguments */
  while((opt = getopt(argc, argv, "D:Ee:Ii:Oo:Vv")) != EOF){
    switch(opt){
    case 'D':
      /* The debugging level.  Default is 0. */
      debug = (unsigned short int)atoi(optarg);
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
  
  /* start the clock */
  clock=time((time_t *)NULL);
  time_buf_start=ctime(&clock);
  (void)fprintf(stderr,"\tstart = %s",time_buf_start);
  
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
  
  if(
     ((matrix_c_ptr=
       (float *)malloc((num_layers+2)*(num_sizes+2)*sizeof(float))) == NULL ) ||
     ((matrix_c_trans_ptr=
       (float *)malloc((num_layers+2)*(num_sizes+2)*sizeof(float))) == NULL ) ||
     ((matrix_cf_ptr=
       (float *)malloc((num_layers+2)*(num_sizes+2)*sizeof(float))) == NULL ) ||
     ((matrix_f_ptr=
       (float *)malloc((num_layers+2)*(num_sizes+2)*sizeof(float))) == NULL ) ||
     False ){
    (void)fprintf(fp_err,"Unable to allocate array in main\n");
    exit(1);
  } /* end if */

  if(
     ((matrix_c=
       (float **)malloc((num_layers+2)*sizeof(float *))) == NULL ) ||
     ((matrix_c_trans=
       (float **)malloc((num_sizes+2)*sizeof(float *))) == NULL ) ||
     ((matrix_cf=
       (float **)malloc((num_layers+2)*sizeof(float *))) == NULL ) ||
     ((matrix_f=
       (float **)malloc((num_layers+2)*sizeof(float *))) == NULL ) ||
     False ){
    (void)fprintf(fp_err,"Unable to allocate array in main\n");
    exit(1);
  } /* end if */

  /* REMEMBER: these matrices are indexed as [layer][size] */
  for(layer=0;layer<num_layers+2;layer++){
    matrix_c[layer]=matrix_c_ptr+(num_sizes+2)*layer;    
    matrix_cf[layer]=matrix_cf_ptr+(num_sizes+2)*layer;    
    matrix_f[layer]=matrix_f_ptr+(num_sizes+2)*layer;    
  } /* end loop over layers */

  /* REMEMBER: these matrices are indexed as [size][layer] */
  for(size=0;size<num_sizes+2;size++){
    matrix_c_trans[size]=matrix_c_trans_ptr+(num_layers+2)*size;    
  } /* end loop over sizes */

  for(layer=1;layer<=num_layers;layer++){
    for(size=1;size<=num_sizes;size++){
      matrix_c[layer][size]=layer+.1*size;
      matrix_cf[layer][size]=0.;
      matrix_f[layer][size]=0.;
    } /* end loop over sizes */
  } /* end loop over layers */

  for(layer=1;layer<=num_layers;layer++){
    for(size=1;size<=num_sizes;size++){
      matrix_cf[size][layer]=matrix_c[layer][size];
    } /* end loop over sizes */
  } /* end loop over layers */

  transpose_conc(matrix_c,matrix_c_trans,num_layers,num_sizes);

#ifdef CRAY2
/*  TEST_TRANS(matrix_c_trans_ptr,matrix_f_ptr,&num_layers,&num_sizes);*/
/*  TEST_CF(matrix_cf_ptr,matrix_f_ptr,&num_layers,&num_sizes);*/
  (void)fprintf(fp_err,"using CRAY2 commands\n");
  TEST_C(matrix_c_ptr,matrix_f_ptr,&num_layers,&num_sizes);
#elif defined(PEACE)
/*  test_trans_(matrix_c_trans_ptr,matrix_f_ptr,&num_layers,&num_sizes);*/
/*  test_cf_(matrix_cf_ptr,matrix_f_ptr,&num_layers,&num_sizes);*/
  test_c_(matrix_c_ptr,matrix_f_ptr,&num_layers,&num_sizes);
  (void)fprintf(fp_err,"using peace ( = PEACE) commands\n");
#else
  (void)fprintf(fp_err,"unknown machine.\n");
#endif

  if(True){
    if(fp_out == stdout){
      (void)fprintf
	(fp_err,"Not writing binary data to stdout (screen). Continuing ...\n");
      return 0 ;
    } /* end if */
    
    if(
       (fwrite(matrix_c_ptr,(num_layers+2)*(num_sizes+2)*sizeof(float),1,fp_out) != 1) ||
       False){
      (void)fprintf(fp_err,"Unable to write correct data array in main\n");
      exit(1);
    }
    (void)fflush(fp_out);
    (void)fprintf(fp_err,"Wrote out cloud data to %s\n",out_file);
  } /* end if True */
  (void)fclose(fp_out);

  if(True){
    if(fp_in == stdin){
      (void)fprintf
	(fp_err,"Not reading binary data to stdin (screen). Continuing ...\n");
      return 0 ;
    } /* end if */
    
    if(
       (fread(matrix_c_ptr,(num_layers+2)*(num_sizes+2)*sizeof(float),1,fp_in) != 1) ||
       False){
      (void)fprintf(fp_err,"Unable to read correct data array in main\n");
      exit(1);
    }
    (void)fprintf(fp_err,"Read in cloud data from %s\n",in_file);
  } /* end if True */
  (void)fclose(fp_in);

  if(debug == 1){
    for(layer=1;layer<=num_layers;layer++){
      for(size=1;size<=num_sizes;size++){
	(void)fprintf(fp_err,"matrix_c[%i][%i] = %g, matrix_cf[%i][%i] = %g, matrix_f[%i][%i] = %g, matrix_f[%i][%i] = %g\n",
		      layer,size,matrix_c[layer][size],
		      size,layer,matrix_cf[size][layer],
		      layer,size,matrix_f[layer][size],
		      size,layer,matrix_f[size][layer]);
      } /* end loop over sizes */
    } /* end loop over layers */
  } /* end debug */
  
  if(debug == 2){
    for(layer=1;layer<=num_layers;layer++){
      for(size=1;size<=num_sizes;size++){
	(void)fprintf(fp_err,"matrix_c[%i][%i] = %g, matrix_c_trans[%i][%i] = %g, matrix_f[%i][%i] = %g, matrix_f[%i][%i] = %g\n",
		      layer,size,matrix_c[layer][size],
		      size,layer,matrix_c_trans[size][layer],
		      layer,size,matrix_f[layer][size],
		      size,layer,matrix_f[size][layer]);
      } /* end loop over sizes */
    } /* end loop over layers */
  } /* end debug */
  
  Exit_gracefully();
} /* end main() */

void print_usage()
{
  (void)fprintf(stderr,"\nusage: clouds [-options] where options are one or more of:\n\n");
  (void)fprintf(stderr,"D:Ee:Ii:Oo:vV\n\n");
  (void)fprintf(stderr,"-D debug The debugging level.  Default is 0.\n");
  (void)fprintf(stderr,"-E STDERR Toggle stderr stream.  Default is True\n");
  (void)fprintf(stderr,"-I STDIN Toggle stdin stream.  Default is True\n");
  (void)fprintf(stderr,"-O STDOUT Toggle stdout stream.  Default is False\n");
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
  (void)fprintf(fp_err,"\tfinish = %s\n",time_buf_finish);

  (void)fclose(fp_err);
  (void)fclose(fp_in);
  (void)fclose(fp_out);

  exit(1);
}

void transpose_conc(concentration,conc_transpose,num_layers,num_sizes)
     float **concentration,**conc_transpose;
     int num_layers,num_sizes;
{
  /* Fill the conc_transpose matrix with the transpose of the 
     concentration matrix */

  int layer,size;

  for(layer=0;layer<=num_layers+1;layer++){
    for(size=0;size<=num_sizes+1;size++){
      conc_transpose[size][layer]=concentration[layer][size];
    } /* end loop over sizes */
  } /* end loop over layers */
}




