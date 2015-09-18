#!/usr/local/bin/perl 
'di';
'ig00';
# --  # -*-Perl-*- Automatically loads Perl mode in Emacs
# $Id$

# Purpose: Takes a CCM run script output file as input. 
# Output is a netCDF file containing the global statistics.

# Usage:
# /home/zender/perl/stat2nc.pl /fs/cgd/home0/zender/dst/cray/dst02.txt /data/zender/dst/dgn/dst02_dgn.nc
# /home/zender/perl/stat2nc.pl /data/zender/dst/dgn/dst04.txt /data/zender/dst/dgn/dst04_dgn.nc

require <getopts.pl>;
require <importenv.pl>;
require <pwd.pl>;

# Required by pwd.pl
&initpwd;

# Set output flushing to help debugging on hard crashes
# These options update the filehandle after every output statement
# See the Perl manual, p. 110.
select((select(STDOUT),$|=1)[0]);
select((select(STDERR),$|=1)[0]);

# Timing information
($start_user_time,$start_system_time,$start_child_user_time,$start_child_system_time)=times;
$local_date_time=&date_time();
#printf STDERR ("Start time: %s\n",$local_date_time);

# Set defaults 
$False=0;
$True=1;

$Boolean=$True;
$HELP=$False;
$RCS_Header="$Id$";
$RCS_Revision="$Revision$";
$REMOVE_CDL_FILE=$True;
$ZERO_CFL_vrt_max=$False;
$debug=0;
$float_foo=102.543;
$pi=atan2(1,1)*4.;
				 
# Parse command line arguments 
$cmd_line=join(" ",$0,@ARGV);
&Getopts('BD:f:hrz');
if (defined($opt_B)) {$Boolean=!$Boolean};
if (defined($opt_D)) {$debug=$opt_D};
if (defined($opt_h)) {$HELP=!$HELP};
if (defined($opt_f)) {$float_foo=$opt_f};
if (defined($opt_r)) {$REMOVE_CDL_FILE=!$REMOVE_CDL_FILE};
if (defined($opt_z)) {$ZERO_CFL_vrt_max=!$ZERO_CFL_vrt_max};

if($HELP){		 
    &prn_usg();
    exit 0;
} # end HELP

if($#ARGV < 1){
    &prn_usg();
    exit 0;
} # end if

# the input file name should be the second to last positional argument
$in_file=$ARGV[$#ARGV-1];

# the output file name should be the last positional argument
$out_file=$ARGV[$#ARGV];

if($debug == 1){		 
    printf STDERR ("Reading input CCM data from %s\n",$in_file);
    printf STDERR ("Writing output CDL data to %s.cdl\n",$out_file);
    printf STDERR ("Writing output netCDF data to %s\n",$out_file);
} # end debug 

open(IN_FILE,$in_file) || die "cannot open $in_file";

printf STDERR ("Reading %s\n",$in_file);
while(<IN_FILE>){
    if (/^ NSTEP =/){

# Example line: See the camel book p. 28 for parsing examples.
# NSTEP = 53 8.541418100014336E-05 1.827422388945251E-06 252.895 9.85095E+04 2.925252412325957E+01 0.08 0.02
# This method doesn't work for some reason.
# 	/^ NSTEP = (\S+) (\S+) (\S+) (\S+) (\S+) (\S+) (\S+) (\S+)$/;
#	push(@time_step,$1);
#	push(@vrt_rms,$2);
#	push(@dvr_rms,$3);
#	push(@tpt_rms,$4);
#	push(@prs_sfc,$5);
#	push(@mpc_H2O,$6);
#	push(@CFL_hrz_max,$7);
#	push(@CFL_vrt_max,$8);
#
	chop;
	s/NSTEP \=//g;
	@line=split(/ +/,$_);
	push(@time_step,$line[1]);	# $time_step[++$#time_step]=$line[3]
	push(@vrt_rms,$line[2]);
	push(@dvr_rms,$line[3]);
	push(@tpt_rms,$line[4]);
	push(@prs_sfc,$line[5]);
	push(@mpc_H2O,$line[6]);
	push(@CFL_hrz_max,$line[7]);
	push(@CFL_vrt_max,$line[8]);
#
    } # end if
} # end while

close IN_FILE;

# Backwards compatability to when CFL_vrt_max wasn't written
if ($ZERO_CFL_vrt_max){
#    grep($_ *= 0.,@CFL_vrt_max); # camel p. 
    @CFL_vrt_max= (0.) x @CFL_vrt_max # camel p. 83
} # end if

# Verify we found all the time steps
if($#time_step != $time_step[$#time_step]){
    printf STDERR ("\nWarning:\n Array length = %d doesn't match last time step value = %d.\n This is fine if you are analyzing a restart, regeneration or branch\n (NSREST= 1, 2, or 3) run. Otherwise suspect that the timestep\n line-matching algorithm isn't matching all the right lines.\n\n",
		  $#time_step,
		  $time_step[$#time_step]);
}

if($debug >= 1){
    printf STDERR ("Size of RMS arrays = %d, last value = %d\n",$#time_step-$[+1,$time_step[$#time_step]);
} # end debug 

if($debug == 1){		 
    printf STDERR ("%4s %10s %10s %10s %10s %10s %10s %10s\n",
                   "Step","Vorticity","Divergence","Temperatur",
                   "Moisture","Sfc. Pres.","Hor. CFL","Vert. CFL");
    printf STDERR ("%4s %10s %10s %10s %10s %10s %10s %10s\n",
                   "","s-1","s-1","K",
                   "kg m-2","Pa","fraction","fraction");
    for ($idx=0;$idx<=$#time_step;$idx++){
	printf STDERR ("%4d %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n",
		       $time_step[$idx],
		       $vrt_rms[$idx],
		       $dvr_rms[$idx],
		       $tpt_rms[$idx],
		       $mpc_H2O[$idx],
		       $prs_sfc[$idx],
		       $CFL_hrz_max[$idx],
		       $CFL_vrt_max[$idx]
		       );
    } # endfor
} # end debug 

open(OUT_FILE,">".$out_file.".cdl") || die "cannot open $out_file.cdl";
printf STDERR ("Creating %s.cdl\n",$out_file);

$cdl=
    sprintf("netcdf %s {\n",$out_file).
    "dimensions:\n".
#    sprintf("\ttime = %d;\n",$#time_step-$[+1).
    sprintf("\ttime = UNLIMITED;\n",$#time_step-$[+1).
    "variables:\n".
    sprintf("\t\t:cmd_line=\"%s\";\n",$cmd_line).
    sprintf("\t\t:creation_date=\"%s\";\n",$local_date_time).
    sprintf("\t\t:RCS_Header=\"%s\";\n",$RCS_Header).
    "\tint\ttime_step(time);\n".
    "\t\ttime_step:long_name=\"Model timestep\";\n".
    "\t\ttime_step:units=\"# steps\";\n".
    "\tfloat vrt_rms(time);\n".
    "\t\tvrt_rms:long_name=\"Global root mean square vorticity\";\n".
    "\t\tvrt_rms:units=\"second-1\";\n".
    "\tfloat dvr_rms(time);\n".
    "\t\tdvr_rms:long_name=\"Global root mean square divergence\";\n".
    "\t\tdvr_rms:units=\"second-1\";\n".
    "\tfloat tpt_rms(time);\n".
    "\t\ttpt_rms:long_name=\"Global root mean square temperature\";\n".
    "\t\ttpt_rms:units=\"kelvin\";\n".
    "\tfloat mpc_H2O(time);\n".
    "\t\tmpc_H2O:long_name=\"Global average column water vapor\";\n".
    "\t\tmpc_H2O:units=\"kilogram meter-2\";\n".
    "\tfloat prs_sfc(time);\n".
    "\t\tprs_sfc:long_name=\"Global average surface pressure\";\n".
    "\t\tprs_sfc:units=\"pascal\";\n".
    "\tfloat CFL_hrz_max(time);\n".
    "\t\tCFL_hrz_max:long_name=\"Global maximum horizontal CFL condition\";\n".
    "\t\tCFL_hrz_max:units=\"fraction\";\n".
    "\tfloat CFL_vrt_max(time);\n".
    "\t\tCFL_vrt_max:long_name=\"Global maximum vertical CFL condition\";\n".
    "\t\tCFL_vrt_max:units=\"fraction\";\n".
    "data:\n";
#    "\n".
#    sprintf("\n").

# Print out the CDL file from which the netCDF file is constructable.
print OUT_FILE $cdl;
print OUT_FILE "\ttime_step\t= ".join(",",@time_step).";\n";
print OUT_FILE "\tvrt_rms\t= ".join(",",@vrt_rms).";\n";
print OUT_FILE "\tdvr_rms\t= ".join(",",@dvr_rms).";\n";
print OUT_FILE "\ttpt_rms\t= ".join(",",@tpt_rms).";\n";
print OUT_FILE "\tmpc_H2O\t= ".join(",",@mpc_H2O).";\n";
print OUT_FILE "\tprs_sfc\t= ".join(",",@prs_sfc).";\n";
print OUT_FILE "\tCFL_hrz_max\t= ".join(",",@CFL_hrz_max).";\n";
print OUT_FILE "\tCFL_vrt_max\t= ".join(",",@CFL_vrt_max).";\n";
print OUT_FILE "\n}\n";
close OUT_FILE;

$ncgen_cmd = "ncgen -b -o ".$out_file." ".$out_file.".cdl";
printf STDERR ("Executing %s\n",$ncgen_cmd);
open(NCGEN,$ncgen_cmd." |") || die "cannot execute $ncgen_cmd, stopped";
while(<NCGEN>){
    printf STDERR ("%s",$_);
} # endwhile
close NCGEN;

# Removing the CDL file is recommended because CDL files are large
# and can always be reconstructed with the ncgen command.
if ($REMOVE_CDL_FILE){
    $rm_cmd = "rm -f ".$out_file.".cdl";
    printf STDERR ("Executing %s\n",$rm_cmd);
    open(RM,$rm_cmd." |") || die "cannot execute $rm_cmd, stopped";
    while(<RM>){
	printf STDERR ("%s",$_);
    } # endwhile
    close RM;
} # endif

# Timing information
($end_user_time,$end_system_time,$end_child_user_time,$end_child_system_time)=times;
#printf STDERR ("Elapsed time: User CPU %.2f s, System CPU %.2f s\n",$end_user_time-$start_user_time,$end_system_time-$start_system_time);
#printf STDERR ("Elapsed time: Child user CPU %.2f s, Child system CPU %.2f s\n",$end_child_user_time-$start_child_user_time,$end_child_system_time-$start_child_system_time);
$local_date_time=&date_time();
#printf STDERR ("End time: %s\n",$local_date_time);

sub prn_usg{
    print STDERR "\
Usage:  stat2nc.pl [options] in_file out_file\

   -D dbg_lvl      Sets the debugging level\
   -h              Display this help message\
   -r              Do not remove intermediate netCDF CDL file\
   -z              Fill the CFL_vrt_max field with zeroes\
   in_file         Name of the ASCII input file to process\
   out_file        Name of the binary netCDF output file\
                   The intermediate CDL file will be named out_file.cdl\n";
} # end sub prn_usg()

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

###############################################################

# These next few lines are legal in both Perl and nroff.

.00;                            # finish .ig
 
'di           \" finish diversion--previous line must be blank
.nr nl 0-1    \" fake up transition to first page again
.nr % 0       \" start at page 1
'; __END__    ##### From here on it's a standard manual page #####

.de SB          \" small and bold
.if !"\\$1"" \\s-2\\fB\&\\$1\\s0\\fR\\$2 \\$3 \\$4 \\$5
..
.\"
.de T           \" switch to typewriter font
.ft CW          \" probably want CW if you don't have TA font
..
.\"
.de TY          \" put $1 in typewriter font
.if t .T
.if n ``\c
\\$1\c
.if t .ft P
.if n \&''\c
\\$2
..
.\"
.de M           \" man page reference
\\fI\\$1\\fR\\|(\\$2\)\\$3
..
.TH STAT2NC.PL 1 "March 6, 1995"
.AT 3
.SH NAME
stat2nc.pl \- build netCDF files from CCM2 text output file
.SH SYNOPSIS
.nf
.B stat2nc.pl [ \-D dbg_lvl ] [ \-r ] [ \-z ] [[\-i] input_file] [[\-o] output_file]
.fi
.SH DESCRIPTION
stat2nc.pl is a filter that converts the input stream (input_file) from a
standard ASCII CCM2 output file into a binary netCDF file format
(output_file) containing all the global diagnostics from each
timestep.  stat2nc.pl works by searching the input file for lines of the
sort

NSTEP = 53 8.541418100014336E-05 1.827422388945251E-06 252.895
9.85095E+04 2.925252412325957E+01 0.08 0.02 

stat2nc.pl should work when you concatenate many restart files together,
as long as the time steps are in sequential order. 
stat2nc.pl takes as input the ASCII textual output file from a CCM2 run
(e.g., sample.input). stat2nc.pl assembles an ASCII netCDF file (e.g.,
sample.output.cdl) using the CDL language, then stat2nc.pl invokes ncgen
to convert the file to binary netCDF format (e.g., sample.output.nc),
then stat2nc.pl deletes the CDL output file because it is often very
large, and it can always be recreated with ncgen.  Each CCM2 global is
stored in its own one dimensional array, so the output netCDF file has
eight variables and one dimension:

.nf
long time_step(time)            Model timestep
float vrt_rms(time)            Global root mean square vorticity
float dvr_rms(time)             Global root mean square divergence
float tpt_rms(time)               Global root mean square temperature
float mpc_H2O(time)       Global average column water vapor
float prs_sfc(time)             Global average surface pressure
float CFL_hrz_max(time)         Global maximum horizontal CFL condition
float CFL_vrt_max(time)        Global maximum vertical CFL condition
.fi

Use a netCDF viewer to examine the output.  The textual viewer
available at all netCDF sites is ncdump. For a graphical viewer, write
your own, or, if your site licenses IDL, you can use the publicly
available IDL netCDF viewer called ncbrowse, available via WWW at
ftp://ftp.unidata.ucar.edu:/pub/gopherd/gopher-data/anon.ftp/netcdf/contrib/ncbrowse.pro

.SH OPTIONS
.IP "\fB\-h\fP"
Display small help message.
.IP "\fB\-d dbg_lvl\fP"
Set the verbosity of the diagnostic information reported to stderr.
.IP "\fB\-r\fP"
Turn off automatic removal of the intermediate netCDF CDL meta-language file.
.IP "\fB\-z\fP"
Automatcially fill the CFL_vrt_max field with zeroes. Use this option with output from earlier versions of CCM which did not report the CFL_vrt_max field in the diagnostic output.
.IP "\fB\-i input_file\fP"
Name of the ASCII input file to process. Default is to read from stdin.
.IP "\fB\-o output_file\fP"
Name stem of the netCDF output files. The intermediate CDL language file will be name output_file.cdl. The binary netCDF file will be named output_file.nc. When no output_file is given, the output files will be written to input_file.cdl and input_file.nc, respectively.
.SH EXAMPLE
The sample input file is just the output file from a one month restart
run of the CCM2 (version 2.1), so sample.output has 2160 timesteps.

.nf
/home/zender/web/ccm: stat2nc.pl -i sample.input -o sample.output
Start time: Tue Feb 21 20:56:08 1995

Warning:
 Array length = 2159 doesn't match last time step value = 8784.
 This is fine if you are analyzing a restart, regeneration or branch
 (NSREST= 1, 2, or 3) run. Otherwise suspect that the timestep line-matching algorithm isn't matching all the right lines.

Executing ncgen -b -o sample.output.nc sample.output.cdl
Executing rm -f sample.output.cdl
Elapsed time: User CPU 2.88 s, System CPU 0.77 s
Elapsed time: Child user CPU 8.50 s, Child system CPU 0.33 s
End time: Tue Feb 21 20:56:23 1995
/home/zender/web/ccm: ncdump -h sample.output.nc
netcdf sample.output {
dimensions:
        time = 2160 ;

variables:
        long time_step(time) ;
                time_step:long_name = "Model timestep" ;
                time_step:units = "# steps" ;
        float vrt_rms(time) ;
                vrt_rms:long_name = "Global root mean square vorticity" ;
                vrt_rms:units = "second-1" ;
        float dvr_rms(time) ;
                dvr_rms:long_name = "Global root mean square divergence" ;
                dvr_rms:units = "second-1" ;
        float tpt_rms(time) ;
                tpt_rms:long_name = "Global root mean square temperature" ;
                tpt_rms:units = "kelvin" ;
        float mpc_H2O(time) ;
                mpc_H2O:long_name = "Global average column water vapor" ;
                mpc_H2O:units = "kilogram meter-2" ;
        float prs_sfc(time) ;
                prs_sfc:long_name = "Global average surface pressure" ;
                prs_sfc:units = "pascal" ;
        float CFL_hrz_max(time) ;
                CFL_hrz_max:long_name = "Global maximum horizontal CFL condition" ;
                CFL_hrz_max:units = "fraction" ;
        float CFL_vrt_max(time) ;
                CFL_vrt_max:long_name = "Global maximum vertical CFL condition" ;
                CFL_vrt_max:units = "fraction" ;

// global attributes:
                :cmd_line = "stat2nc.pl -r -i sample.input -o sample.output" ; 
                :creation_date = "Tue Feb 21 21:00:16 1995" ;
                :RCS_Header = ": /home/zender/perl/RCS/stat2nc.pl,v 1.3 1995/02/22 03:22:41 zender Exp $" ;
}
/home/zender/web/ccm:  
.fi
.SH CAVEATS
Some early
versions of CCM2 did not print out the vertical CFL condition each
timestep. If you know perl, then stat2nc.pl is easy to modify to
handle situations with missing/extra diagnostics; see the -z option
above.  
.SH ENVIRONMENT
No environment variables are used.
.SH FILES
None.
.SH AUTHOR
Charlie Zender
.SH "SEE ALSO"
ncgen(1), ncdump(1), netcdf(3), perl(1)
.SH DIAGNOSTICS
stat2nc.pl will report instances where the CCM2 timestep number is not in sync with the number of lines meeting the search criteria ("NSTEP = .."). 
.SH BUGS
stat2nc.pl has not been tested on very long CCM2 runs ( > 6000 timesteps).
