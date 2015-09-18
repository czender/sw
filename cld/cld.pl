#!/usr/local/bin/perl 
#! /data/zender/bin/SUNMP/perl 
# --  # -*-Perl-*- Automatically loads Perl mode in Emacs
# $Id$

# Purpose: Batch run the cld program

# Example usage:
# cld.pl >&! cld.txt &

# $Log: not supported by cvs2svn $
# Revision 1.1.1.1  1998-09-15 02:06:40  zender
# Imported sources
#
# Revision 1.2  1996/01/05  19:49:47  zender
# version used to generate updraft sensitivity study for 10 km
# high cloud with updrafts from 1-20 cm/s.
#
# Revision 1.1  1995/06/18  02:12:47  zender
# Initial revision

require "getopts.pl";
require "importenv.pl";
require "pwd.pl";

# Required by pwd.pl
&initpwd;

# Set output flushing to help debugging on hard crashes. 
# These options update the filehandle after every output statement.
# See the Perl manual, p. 110.
select((select(STDOUT),$|=1)[0]);
select((select(STDERR),$|=1)[0]);

# Timing information
($user_time,$system_time,$child_user_time,$child_system_time)=times;
$srt_user_time=$user_time;
$srt_system_time=$system_time;
$srt_child_user_time=$child_user_time;
$srt_child_system_time=$child_system_time;
$local_date_time=&date_time();
printf STDERR ("Start time: %s\n",$local_date_time);

# Set defaults 
$False=0;
$True=1;
#
$Boolean=$True;
$RCS_Header="\$Id$";
$RCS_Revision="\$Revision$";
$debug=0;
$float_foo=102.543;
$hgt_srt=5000.;
$hgt_end=16000.;
$hgt_incr=1000.;
$pi=atan2(1,1)*4.;
				 
# Parse command line arguments 
$cmd_line=join(" ",$0,@ARGV);
&Getopts('BD:e:i:f:s:v?');
if (defined($opt_B)) {$Boolean=!$Boolean};
if (defined($opt_D)) {$debug=$opt_D};
if (defined($opt_f)) {$float_foo=$opt_f};
if (defined($opt_s)) {$hgt_srt=$opt_s};
if (defined($opt_e)) {$hgt_end=$opt_e};
if (defined($opt_i)) {$hgt_incr=$opt_i};
if (defined($opt_v)) {print STDERR $RCS_Header."\n"};

if($debug == 1){		 
    print STDERR "Processing data...\n";
    print STDERR "Command Line = ".$cmd_line."\n";
} # end debug 

$hgt_srt=5000.;
$hgt_end=16000.;
$hgt_incr=1000.;
$w_srt=.0;
$w_end=.0;
$w_incr=.02;
$RH_ice_srt=1.;
$RH_ice_end=1.;
$RH_ice_incr=.1;
$RH_liq=.8;
for ($hgt=$hgt_srt;$hgt<=$hgt_end;$hgt+=$hgt_incr){
    for ($w=$w_srt;$w<=$w_end;$w+=$w_incr){
	for ($RH_ice=$RH_ice_srt;$RH_ice<=$RH_ice_end;$RH_ice+=$RH_ice_incr){

# Filenames are cld_HH_WW_RR where HH = hgt in km, WW = updraft in cm/s, RR = RH_ice in percent
	    $out_fl=sprintf("cld_%02d_%02d_%03d.nc",&round($hgt/1000.),&round($w*100.),&round($RH_ice*100.));
	    $cld_cmd=sprintf("cld -C -D 55 -G 5000 -g 13000 -H %05d -h 2000 -k 3. -l 130 -M 1 -m 3.e-6 -N 1.e6 -n 1200 -o %s -p 20 -r 0 -R -S %04.2f -s %04.2f -u -w %04.3f -X -x 50",$hgt,$out_fl,$RH_ice,$RH_liq,$w);
	    printf STDOUT ("%s\n",$cld_cmd);
	    `$cld_cmd`;
#
	    if($PVM_ARCH eq "SUN4" || $PVM_ARCH eq "SUN4SOL2" || $PVM_ARCH eq "SUNMP" || $PVM_ARCH eq "SGI5"){		 
#		$mv_cmd="mv -f ".$out_fl." /data2/zender/cld";
#		printf STDOUT ("%s\n",$mv_cmd);
#		`$mv_cmd`;
	    }elsif($PVM_ARCH eq "RS6K" ||$PVM_ARCH eq "CRAY"){		 
		$rcp_cmd="rcp ".$out_fl." ".$err_fl." heinlein.cgd.ucar.edu:/data2/zender/cld";
		printf STDOUT ("%s\n",$rcp_cmd);
		`$rcp_cmd`;
		
		$rm_cmd="/bin/rm -f ".$out_fl." ".$err_fl;
		printf STDOUT ("%s\n",$rm_cmd);
		`$rm_cmd`;
	    } # end else	
	    
	} # end for
    } # end for
} # end for

# Timing information
($user_time,$system_time,$child_user_time,$child_system_time)=times;
$end_user_time=$user_time;
$end_system_time=$system_time;
$end_child_user_time=$child_user_time;
$end_child_system_time=$child_system_time;
printf STDERR ("Elapsed time: User CPU %.2f s, System CPU %.2f s\n",$end_user_time-$srt_user_time,$end_system_time-$srt_system_time);
printf STDERR ("Elapsed time: Child user CPU %.2f s, Child system CPU %.2f s\n",$end_child_user_time-$srt_child_user_time,$end_child_system_time-$srt_child_system_time);
$local_date_time=&date_time();
printf STDERR ("End time: %s\n",$local_date_time);

sub date_time{
# Returns local time string in "Sun Sep 18 20:14:49 1994" format
# Since working in GMT is a pain, adding the time zone for which
# this time is valid would be nice and remove ambiguity, i.e.,
# "Sun Sep 18 20:14:49 1994 MDT"
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

sub round{
# rounds to nearest integer. from the perl FAQ.
    $number=shift;
    return int($number+.5*($number <=> 0));
} # end sub round()
