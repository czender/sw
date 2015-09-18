#!/usr/local/bin/perl 
# --  # -*-Perl-*- Automatically loads Perl mode in Emacs
########################################################################
# RCS Identification
########################################################################
# $Author: zender $
# $Date$
# $Id$
# $Revision$
# $Locker:  $
# $RCSfile: data2nc.pl,v $
# $Source: /home/zender/cvs/perl/data2nc.pl,v $
# $Id$
# $State: Exp $
#
# NB: get RCS formatting in Perl files by using rcs -U -c"# " stat2nc.pl
#
# Purpose: Convert arbitrary input text data files to netCDF.
#
# Cautions: 
# The usage of " is distinct from that of '. Use " for strings,
# like filenames, but use ' for passing those strings to functions.
# Integer format string to use is %d NOT %i.
# Use "elsif" instead of "else if".
# Splitting lines into words is done with @split_line=split(/ +/,$line);
# the / +/ eliminates groups of one or more spaces.
# It may be neccessary in some circumstances to separate the division 
# operator, /, from the operands by spaces or else syntax checkers may
# mistake it for the beginning of a regexp.
# 				
# Example usage:
# /home/zender/perl/data2nc /data/zender/aca/nrel/nrel_obs_01.txt /data/zender/aca/nrel/nrel_obs_01.nc
# /home/zender/perl/data2nc /data/zender/aca/nrel/nrel_obs_02.txt /data/zender/aca/nrel/nrel_obs_02.nc
#
# $Log: not supported by cvs2svn $
# Revision 1.1.1.1  1998-09-06 06:11:13  zender
# Imported sources
#
########################################################################
#
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
printf STDOUT ("Start time: %s\n",$local_date_time);

# Set defaults 
$False=0;
$True=1;
#
$Boolean=$True;
$HELP=$False;
$RCS_Header="$Id$";
$RCS_Revision="$Revision$";
$REMOVE_CDL_FILE=$True;
$debug=0;
$float_foo=102.543;
$pi=atan2(1,1)*4.;
				 
# Parse command line arguments 
$cmd_line=join(" ",$0,@ARGV);
&Getopts('BD:f:hr');
if (defined($opt_B)) {$Boolean=!$Boolean};
if (defined($opt_D)) {$debug=$opt_D};
if (defined($opt_h)) {$HELP=!$HELP};
if (defined($opt_f)) {$float_foo=$opt_f};
if (defined($opt_r)) {$REMOVE_CDL_FILE=!$REMOVE_CDL_FILE};

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
    printf STDOUT ("Reading input data from %s\n",$in_file);
    printf STDOUT ("Writing output CDL data to %s.cdl\n",$out_file);
    printf STDOUT ("Writing output netCDF data to %s\n",$out_file);
} # end debug 

open(IN_FILE,$in_file) || die "cannot open $in_file";

printf STDOUT ("Reading %s\n",$in_file);

while(<IN_FILE>){

# Throw away the first line, which is a text header
    if($. <= 1){next;}

# Process the rest of the lines.
# Convert the data to SI as we go.
    chop;
    @line=split(/ +/,$_);
    push(@band,$line[1]*1.e-9);	
    push(@flux_spec_down_sfc_cld_1,$line[2]*1.e9);
    push(@flux_spec_down_sfc_cld_2,$line[3]*1.e9);
    push(@flux_spec_down_sfc_cld_3,$line[4]*1.e9);
    push(@flux_spec_down_sfc_clr,$line[5]*1.e9);

} # end while

close IN_FILE;

if($debug == 1){		 
    printf STDOUT ("%10s %10s %10s %10s %10s %10s %10s %10s\n",
                   "Band","Flux Cld","Flux Cld","Flux Clr");
    printf STDOUT ("%10s %10s %10s %10s %10s %10s %10s %10s\n",
                   "meter","w/m2/m","w/m2/m","w/m2/m");
    for ($idx=0;$idx<=$#band;$idx++){
	printf STDOUT ("%10.4e %10.4e %10.4e %10.4e\n",
		       $band[$idx],
		       $flux_spec_down_sfc_cld_1[$idx],
		       $flux_spec_down_sfc_cld_2[$idx],
		       $flux_spec_down_sfc_clr[$idx]
		       );
    } # endfor
} # end debug 

@wavelen=@band;
@wavelen_width=@band;
for ($idx=0;$idx<=$#band-1;$idx++){
    $wavelen_width[$idx]=$band[$idx+1]-$band[$idx];
} # endfor
$wavelen_width[$#band]=$wavelen_width[$#band-1];

$flux_bb_down_sfc_cld_1=0.;
$flux_bb_down_sfc_cld_2=0.;
$flux_bb_down_sfc_cld_3=0.;
$flux_bb_down_sfc_clr=0.;
for ($idx=0;$idx<=$#band;$idx++){
    $flux_bb_down_sfc_cld_1+=$flux_spec_down_sfc_cld_1[$idx]*$wavelen_width[$idx];
    $flux_bb_down_sfc_cld_2+=$flux_spec_down_sfc_cld_2[$idx]*$wavelen_width[$idx];
    $flux_bb_down_sfc_cld_3+=$flux_spec_down_sfc_cld_3[$idx]*$wavelen_width[$idx];
    $flux_bb_down_sfc_clr+=$flux_spec_down_sfc_clr[$idx]*$wavelen_width[$idx];
} # endfor

if($debug == 2){
    printf STDOUT ("Size of arrays = %d, last value = %d\n",$#band-$[+1,$band[$#band]);
} # end debug 

open(OUT_FILE,">".$out_file.".cdl") || die "cannot open $out_file.cdl";
printf STDOUT ("Creating %s.cdl\n",$out_file);

$cdl=
    sprintf("netcdf %s {\n",$out_file).
    "dimensions:\n".
    sprintf("\tband = %d;\n",$#band-$[+1).
    "variables:\n".
    sprintf("\t\t:cmd_line=\"%s\";\n",$cmd_line).
    sprintf("\t\t:creation_date=\"%s\";\n",$local_date_time).
    sprintf("\t\t:RCS_Header=\"%s\";\n",$RCS_Header).
    "\tfloat\tflux_bb_down_sfc_cld_1;\n".
    "\tfloat\tflux_bb_down_sfc_cld_2;\n".
    "\tfloat\tflux_bb_down_sfc_cld_3;\n".
    "\tfloat\tflux_bb_down_sfc_clr;\n".
    "\tfloat\tband(band);\n".
    "\t\tband:long_name=\"Wavelength at band center\";\n".
    "\t\tband:units=\"meter\";\n".
    "\tfloat\twavelen(band);\n".
    "\t\twavelen:long_name=\"Wavelength at band center\";\n".
    "\t\twavelen:units=\"meter\";\n".
    "\tfloat\twavelen_width(band);\n".
    "\t\twavelen_width:long_name=\"bandpass of each spectral interval\";\n".
    "\t\twavelen_width:units=\"meter\";\n".
    "\tfloat flux_spec_down_sfc_cld_1(band);\n".
    "\t\tflux_spec_down_sfc_cld_1:long_name=\"cloudy sky broadband downwelling spectral flux at the surface\";\n".
    "\t\tflux_spec_down_sfc_cld_1:units=\"watt meter-2 meter-1\";\n".
    "\tfloat flux_spec_down_sfc_cld_2(band);\n".
    "\t\tflux_spec_down_sfc_cld_2:long_name=\"cloudy sky broadband downwelling spectral flux at the surface\";\n".
    "\t\tflux_spec_down_sfc_cld_2:units=\"watt meter-2 meter-1\";\n".
    "\tfloat flux_spec_down_sfc_cld_3(band);\n".
    "\t\tflux_spec_down_sfc_cld_3:long_name=\"cloudy sky broadband downwelling spectral flux at the surface\";\n".
    "\t\tflux_spec_down_sfc_cld_3:units=\"watt meter-2 meter-1\";\n".
    "\tfloat flux_spec_down_sfc_clr(band);\n".
    "\t\tflux_spec_down_sfc_clr:long_name=\"clear sky broadband downwelling spectral flux at the surface\";\n".
    "\t\tflux_spec_down_sfc_clr:units=\"watt meter-2 meter-1\";\n".
    "data:\n";
#    "\n".
#    sprintf("\n").

# Print out the CDL file from which the netCDF file is constructable.
print OUT_FILE $cdl;
print OUT_FILE "\tflux_bb_down_sfc_cld_1\t= ".join(",",$flux_bb_down_sfc_cld_1).";\n";
print OUT_FILE "\tflux_bb_down_sfc_cld_2\t= ".join(",",$flux_bb_down_sfc_cld_2).";\n";
print OUT_FILE "\tflux_bb_down_sfc_cld_3\t= ".join(",",$flux_bb_down_sfc_cld_3).";\n";
print OUT_FILE "\tflux_bb_down_sfc_clr\t= ".join(",",$flux_bb_down_sfc_clr).";\n";
print OUT_FILE "\tband\t= ".join(",",@band).";\n";
print OUT_FILE "\twavelen\t= ".join(",",@wavelen).";\n";
print OUT_FILE "\twavelen_width\t= ".join(",",@wavelen_width).";\n";
print OUT_FILE "\tflux_spec_down_sfc_cld_1\t= ".join(",",@flux_spec_down_sfc_cld_1).";\n";
print OUT_FILE "\tflux_spec_down_sfc_cld_2\t= ".join(",",@flux_spec_down_sfc_cld_2).";\n";
print OUT_FILE "\tflux_spec_down_sfc_cld_3\t= ".join(",",@flux_spec_down_sfc_cld_3).";\n";
print OUT_FILE "\tflux_spec_down_sfc_clr\t= ".join(",",@flux_spec_down_sfc_clr).";\n";
print OUT_FILE "\n}\n";
close OUT_FILE;

$ncgen_cmd = "ncgen -b -o ".$out_file." ".$out_file.".cdl";
printf STDOUT ("Executing %s\n",$ncgen_cmd);
open(NCGEN,$ncgen_cmd." |") || die "cannot execute $ncgen_cmd, stopped";
while(<NCGEN>){
    printf STDOUT ("%s",$_);
} # endwhile
close NCGEN;

# Removing the CDL file is recommended because CDL files are large
# and can always be reconstructed with the ncgen command.
if ($REMOVE_CDL_FILE){
    $rm_cmd = "rm -f ".$out_file.".cdl";
    printf STDOUT ("Executing %s\n",$rm_cmd);
    open(RM,$rm_cmd." |") || die "cannot execute $rm_cmd, stopped";
    while(<RM>){
	printf STDOUT ("%s",$_);
    } # endwhile
    close RM;
} # endif

# Timing information
($user_time,$system_time,$child_user_time,$child_system_time)=times;
$end_user_time=$user_time;
$end_system_time=$system_time;
$end_child_user_time=$child_user_time;
$end_child_system_time=$child_system_time;
printf STDOUT ("Elapsed time: User CPU %.2f s, System CPU %.2f s\n",$end_user_time-$start_user_time,$end_system_time-$start_system_time);
printf STDOUT ("Elapsed time: Child user CPU %.2f s, Child system CPU %.2f s\n",$end_child_user_time-$start_child_user_time,$end_child_system_time-$start_child_system_time);
$local_date_time=&date_time();
printf STDOUT ("End time: %s\n",$local_date_time);

sub prn_usg{
    print STDOUT "\
Usage:  data2nc [options] in_file out_file\

   -D dbg_lvl      Sets the debugging level\
   -h              Display this help message\
   -r              Do not remove intermediate netCDF CDL file\
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
