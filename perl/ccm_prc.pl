#!/usr/local/bin/perl 
########################################################################
# RCS Identification
########################################################################
# $Author: zender $
# $Date$
# $Id$
# $Revision$
# $Locker:  $
# $RCSfile: ccm_prc.pl,v $
# $Source: /home/zender/cvs/perl/ccm_prc.pl,v $
# $Id$
# $State: Exp $
#
# NB: get RCS formatting in Perl files by using rcs -U -c"# " proc.pl
#
# Purpose: Generate CCM processor scripts from Perl, using some of 
# Perl's batch job capabilities.
#
# Perl Cautions: 
# The usage of " is distinct from that of '. Use " for strings,
# like filenames, but use ' for passing those strings to functions.
# Integer format string to use is %d NOT %i.
# Use "elsif" instead of "else if".
#
# Splitting lines into words is done with @split_line=split(/ +/,$line);
# the / +/ eliminates groups of one or more spaces.
# It may be neccessary in some circumstances to separate the division 
# operator, /, from the operands by spaces or else syntax checkers may
# mistake it for the beginning of a regexp.
# 				
# Hazards: comments within the UNIX script must begin with a '#' and
# comments within the cat -- END command must begin with a ';' and
# there must not be any comments within the actual namelist.
#
# Processor peccadilloes: ICP's are not allowed in column 1. Do not
# place tab characters at the end of lines. Do not allow ICP 
# continuation lines to begin before the preceding lines, i.e., 
# keep the hanging indent.
#
# Example usage:
# proc.pl -D 1 | m
# proc.pl -d 2
# qsub proc.pl
# setenv DAY_OF_RUN 01; qsub -x -o omega.$DAY_OF_RUN.out proc.pl
# setenv DAY_OF_RUN 02; qsub -x -o omega.$DAY_OF_RUN.out proc.pl
# setenv DAY_OF_RUN 03; qsub -x -o omega.$DAY_OF_RUN.out proc.pl
# setenv DAY_OF_RUN 04; qsub -x -o omega.$DAY_OF_RUN.out proc.pl
# setenv DAY_OF_RUN 05; qsub -x -o omega.$DAY_OF_RUN.out proc.pl
# setenv DAY_OF_RUN 06; qsub -x -o omega.$DAY_OF_RUN.out proc.pl
# setenv DAY_OF_RUN 07; qsub -x -o omega.$DAY_OF_RUN.out proc.pl
# setenv DAY_OF_RUN 08; qsub -x -o omega.$DAY_OF_RUN.out proc.pl
# setenv DAY_OF_RUN 09; qsub -x -o omega.$DAY_OF_RUN.out proc.pl
# setenv DAY_OF_RUN 10; qsub -x -o omega.$DAY_OF_RUN.out proc.pl
# setenv DAY_OF_RUN 11; qsub -x -o omega.$DAY_OF_RUN.out proc.pl
# setenv DAY_OF_RUN 12; qsub -x -o omega.$DAY_OF_RUN.out proc.pl
# setenv DAY_OF_RUN 13; qsub -x -o omega.$DAY_OF_RUN.out proc.pl
# setenv DAY_OF_RUN 14; qsub -x -o omega.$DAY_OF_RUN.out proc.pl
#
# $Log: not supported by cvs2svn $
# Revision 1.1.1.1  1998-09-06 06:11:13  zender
# Imported sources
#
# Revision 1.1  1994/10/04  01:21:10  zender
# Initial revision
#
########################################################################

# Hard code the QSUB directives---these can be overridden from the 
# command line.
## QSUB -q reg -A 35071220 -lM 10Mw -lm 10Mw -lT 10000 -lt 10000 -eo
# QSUB -q reg -A 35071220 -lM 10Mw -lm 10Mw -lT 2000 -lt 2000 -eo

require <getopts.pl>;
require <importenv.pl>;
require <pwd.pl>;

# Required by pwd.pl
&initpwd;

# Set output flushing to help debugging on hard crashes. 
# These options update the filehandle after every output statement.
# See the Perl manual, p. 110.
select((select(STDOUT),$|=1)[0]);
select((select(STDERR),$|=1)[0]);

# Timing information
($user_time,$system_time,$child_user_time,$child_system_time)=times;
$start_user_time=$user_time;
$start_system_time=$system_time;
$start_child_user_time=$child_user_time;
$start_child_system_time=$child_system_time;
$local_date_time=&date_time();
printf STDERR ("Start time: %s\n",$local_date_time);

# Set defaults 
$False=0;
$True=1;
#
$CCM_dir="/usr/tmp/ccm";
$FOUND_MSS_FILE=$False;
$MSS_dir="/ZENDER/proc/omega";
$MSS_retention_period=365;
$NEED_TO_MSWRITE=$False;
$RCP_path="ra.cgd.ucar.edu:/d2/zender";
$cmd_line=$False;
if (defined($DAY_OF_RUN)){
    $day_of_run=$DAY_OF_RUN;
}else{ # end if $DAY_OF_RUN was an environment variable
    $day_of_run=1;
} # end if $DAY_OF_RUN was not an environment variable
$debug=0;
$in_file="STDIN";
$out_file="STDOUT";
$qsub_A="35071220";
$qsub_G="short 1";
$qsub_N="omega_processor";
$qsub_eo=$True;
$qsub_lM=10;
$qsub_lT=2000;
$qsub_o="proc.pl.out";
$qsub_q="reg";
$qsub_s="/usr/local/bin/perl";
				 
# Parse command line arguments 
&Getopts('A:D:d:e:f:G:i:M:m:N:O:o:q:Qs:T:t:');
if (defined($opt_A)) {$qsub_A=$opt_A};
if (defined($opt_D)) {$debug=$opt_D};
if (defined($opt_G)) {$qsub_G=$opt_G};
if (defined($opt_M)) {$qsub_lM=$opt_M};
if (defined($opt_N)) {$qsub_N=$opt_N};
if (defined($opt_Q)) {$cmd_line=!$cmd_line};
if (defined($opt_T)) {$qsub_lT=$opt_T};
if (defined($opt_d)) {$day_of_run=$opt_d};
if (defined($opt_e)) {$qsub_eo=!$qsub_eo};
if (defined($opt_f)) {$float_foo=$opt_f};
if (defined($opt_i)) {$in_file=$opt_i};
if (defined($opt_m)) {$qsub_lm=$opt_m};
if (defined($opt_o)) {$qsub_o=$opt_o};
if (defined($opt_q)) {$qsub_q=$opt_q};
if (defined($opt_s)) {$qsub_s=$opt_s};
if (defined($opt_t)) {$qsub_lt=$opt_t};
#if (defined($opt_f)) {$qsub_=$opt_f};

# Some options best depend on others...
$qsub_lm=10;
$qsub_lt=2000;

# Assemble the NQS directives.
if($OS == "UNICOS"){		 
    $qsub_cmd="";
    $qsub_cmd=$qsub_cmd." -q ".$qsub_q;
    $qsub_cmd=$qsub_cmd." -A ".$qsub_A;
    $qsub_cmd=$qsub_cmd." -lM ".$qsub_lM."Mw";
    $qsub_cmd=$qsub_cmd." -lm ".$qsub_lm."Mw";
    $qsub_cmd=$qsub_cmd." -lT ".$qsub_lT;
    $qsub_cmd=$qsub_cmd." -lt ".$qsub_lt;
# don't give /usr/local/bin/perl to qsub--it barfs.
#    $qsub_cmd=$qsub_cmd." -s ".$qsub_s;
    if($qsub_eo){
	$qsub_cmd=$qsub_cmd." -eo"
    } # endif -eo command
    $qsub_cmd=$qsub_cmd." -o ".$qsub_o;
}elsif($OS == "AIX"){  # endif OS is UNICOS
    $qsub_cmd="";
    $qsub_cmd=$qsub_cmd." -A ".$qsub_A;
    $qsub_cmd=$qsub_cmd." -G ".$qsub_G;
    $qsub_cmd=$qsub_cmd." -N ".$qsub_N;
    $qsub_cmd=$qsub_cmd." -s ".$qsub_s;
    if($qsub_eo){
	$qsub_cmd=$qsub_cmd." -eo"
    } # endif -eo command
    $qsub_cmd=$qsub_cmd." -o ".$qsub_o;
} # endif OS is AIX

if($cmd_line){
    printf STDERR ("qsub%s %s\n",$qsub_cmd,$0);
    exit;
} # endif -eo command

# Execute any preliminary shell commands
&chdir($TMPDIR) || die "Can't chdir to $TMPDIR, stopped";
printf STDERR ("Working directory is now %s\n",$ENV{"PWD"});

# Figure out the task-dependent ICP lines
$base_day=303.0139;
$start_day=$day_of_run-1+$base_day;
$end_day=$start_day+71./72.;
$day_incr=1./72.;
$start_hist_tape_num=$day_of_run*6-5;
$out_file_stub=sprintf("omega.%02d",$day_of_run);
			    
# Form a string for the Processor input
$icp=
    "C JOBSTEP 1. Get the data from the daily mean average tapes.\nC Condense the history save tapes to include only the fields\nC of interest, and get rid of stratospheric levels and levels not in\nC all post-1989 ECMWF tapes.\n".
    " MSPFXI  = \'/ZENDER/ccm2/OMEGA_1/hist/\'\n".
    sprintf(" TAPESA  = \'h%04d\'\n",$start_hist_tape_num).
    sprintf(" DAYSA   = %.7f,%.7f,%.7f\n",$start_day,$end_day,$day_incr).
    " NINTAPA = 6,4,1\n".
    " NDYHSTA = 72\n".
    " TYPEA   = \'CCM1\'\n".
    " FIELDA1 = \'CMFMC\',\'PRECC\',\'PRECT\',\'ORO\'\n".
    "C\n".
    " PRESSLE = 1000.,850.,700.,500.,400.,300.,250.\n".
    "          ,200.,150.,100.\n".
    " LBTDP   =  4\n".
    " PINTXL  = \'CMFMC\',2,4,0\n".
    "C\n".
    " TITLEA  = \'OMEGA Grid Point Fields Pressure Levels\'\n".
    sprintf(" MSRTOA  = \'%d\'\n",$MSS_retention_period).
    " BPHSTA  = \'YES\'\n".
    " OFTHSTA = \'CCM2\'\n".
    sprintf(" SAVHSTA = \'%s/%s.hst\'\n",$MSS_dir,$out_file_stub).
    " PRNTHD  = \'FULL\'\n".
    "ENDOFDATA --------------------------------------------------------------\n".
    "\n";

# Create the processor ICP file
$icp_file = "parms";
open(ICP,">".$icp_file) || die "cannot open $icp_file, stopped";
printf STDERR ("proc ICPs in %s are:\n%s",$icp_file,$icp);
printf ICP ("%s",$icp);
close ICP;

# Execute the processor
$proc_cmd = "/ccm/proc/Processor";
printf STDERR ("Executing %s\n",$proc_cmd);
open(PROC,$proc_cmd." |") || die "cannot execute $proc_cmd, stopped";
while(<PROC>){
    printf STDOUT ("%s",$_);
} # endwhile
close PROC;

# Form a string for the Cond input
$icp=
    ";Comments after semicolons are OK\n".
    "E\$HIST\n".
    " FORMAT  = \'netcdf\'\n".
    " DELINP  = .false.\n".
    sprintf(" MSDIRI  = \'%s\'\n",$MSS_dir).
    sprintf(" FILESI  = \'%s.hst\'\n",$out_file_stub).
#    sprintf(" RCPDIR  = \'%s\'\n",$RCP_path).
    sprintf(" MSDIRO  = \'%s\'\n",$MSS_dir).
    sprintf(" FILESO  = \'%s.nc\'\n",$out_file_stub).
    " DTYPE   = \'DAY\'\n".
    " NTHIST  = 50\n".
    sprintf(" RETPD   = %d\n",$MSS_retention_period).
    " WRPSWD  = \' \' \$\n".
    "\n";

# Create the cond ICP file
$icp_file = "cond.$$";
open(ICP,">".$icp_file) || die "cannot open $icp_file, stopped";
printf STDERR ("cond ICPs in %s are:\n%s",$icp_file,$icp);
printf ICP ("%s",$icp);
close ICP;

# Execute the NetCDF conversion program
$cond_cmd = "/crestone/u1/jet/craybin/condnew < ".$icp_file;
printf STDERR ("Executing %s\n",$cond_cmd);
open(COND,$cond_cmd." |") || die "cannot execute $cond_cmd, stopped";
while(<COND>){
    printf STDOUT ("%s",$_);
} # endwhile
close COND;

if($NEED_TO_MSWRITE){
# Copy the NetCDF file to the MSS
    @in_file_dirs=($TMPDIR,$CCM_dir.$MSS_dir);
    foreach $in_dir (@in_file_dirs) {
	$foo=$in_dir."/".$out_file_stub.".nc";
	if(-e $foo){
	    $in_file=$foo;
	    $FOUND_MSS_FILE=$True;
	    last;
	} # end if the required file was found on disk
    } # end loop over input directories
    if($FOUND_MSS_FILE){
	$mss_cmd=sprintf("mswrite -nowait -t %d %s %s/%s.nc",$MSS_retention_period,$in_file,$MSS_dir,$out_file_stub);
	printf STDERR ("Executing %s\n",$mss_cmd);
	open(MSS,$mss_cmd." |") || die  "cannot perform $mss_cmd";
	while(<MSS>){
	    printf STDOUT ("%s",$_);
	} # endwhile
	close MSS;
	$FOUND_MSS_FILE=$False;
    }else{ # end if the required file was found on disk
	printf STDERR ("Couldn't find %s.nc, looked in...\n%s\n",$out_file_stub,join("\n",@in_file_dirs));
    } # end if the required file couln't be found
    
# Copy the History tape file to the MSS
    @in_file_dirs=($TMPDIR,$CCM_dir.$MSS_dir);
    foreach $in_dir (@in_file_dirs) {
	$foo=$in_dir."/".$out_file_stub.".hst";
	if(-e $foo){
	    $in_file=$foo;
	    $FOUND_MSS_FILE=$True;
	    last;
	} # end if the required file was found on disk
    } # end loop over input directories
    if($FOUND_MSS_FILE){
	$mss_cmd=sprintf("mswrite -nowait -t %d %s %s/%s.hst",$MSS_retention_period,$in_file,$MSS_dir,$out_file_stub);
	printf STDERR ("Executing %s\n",$mss_cmd);
	open(MSS,$mss_cmd." |") || die  "cannot perform $mss_cmd";
	while(<MSS>){
	    printf STDOUT ("%s",$_);
	} # endwhile
	close MSS;
	$FOUND_MSS_FILE=$False;
    }else{ # end if the required file was found on disk
	printf STDERR ("Couldn't find %s.hst, looked in...\n%s\n",$out_file_stub,join("\n",@in_file_dirs));
    } # end if the required file couln't be found
} # end if NEED_TO_MSWRITE

# Timing information
($user_time,$system_time,$child_user_time,$child_system_time)=times;
$end_user_time=$user_time;
$end_system_time=$system_time;
$end_child_user_time=$child_user_time;
$end_child_system_time=$child_system_time;
printf STDERR ("Elapsed time: User CPU %.2f s, System CPU %.2f s\n",$end_user_time-$start_user_time,$end_system_time-$start_system_time);
printf STDERR ("Elapsed time: Child user CPU %.2f s, Child system CPU %.2f s\n",$end_child_user_time-$start_child_user_time,$end_child_system_time-$start_child_system_time);
$local_date_time=&date_time();
printf STDERR ("End time: %s\n",$local_date_time);

sub date_time{
# Returns local time string in "Sun Sep 18 20:14:49 1994" format
    local($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst);
    ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
    $date_time_string=
	sprintf
	    ("%s %s %02d %02d:%02d:%02d 19%02d",
	     (Sun,Mon,Tue,Wed,Thu,Fri,Sat)[$wday],
	     (Jan,Feb,Mar,Apr,May,Jun,Jul,Aug,Sep,Oct,Nov,Dec)[$mon],
	     $mday,
	     $hour,
	     $min,
	     $sec,
	     $year);

    return $date_time_string;
} # end sub date_time()

