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
# $RCSfile: dss3532nc.pl,v $
# $Source: /home/zender/cvs/perl/dss3532nc.pl,v $
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
# /home/zender/perl/dss3532nc /data/zender/aca/dss/raobs.93oct19-23 foo.nc
# /home/zender/perl/dss3532nc /data/zender/aca/dss/raobs.95apr16-19 foo.nc
# /home/zender/perl/dss3532nc /data/zender/aca/dss/raobs.95apr15 foo.nc
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
$CURRENTLY_IN_SONDE_SECTION_OF_PROFILE=$False;
$HELP=$False;
$RCS_Header="$Id$";
$RCS_Revision="$Revision$";
$REMOVE_CDL_FILE=$True;
$debug=0;
$epsilon=.622; # molec wgt vapor/molec wgt dry air
$float_foo=102.543;
$gas_const_dry_air=287.05; # J/kg/K
$lat_heat_vap=2.501e6; # J/kg
$mean_sfc_gravity=9.80665; # m/s^2
$pi=atan2(1,1)*4.;

# Initialize derived fields
$l_over_R=$lat_heat_vap/$gas_const_dry_air;

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

    if (/^ BLOCK /){
	chop;
	@line=split(/ +/,$_);
	$year=$line[7];
	$month=$line[8];
	$day=$line[9];
	$hour=$line[11];
	$lat=$line[13];
	$lon=$line[14];
	$elevation=$line[15];

	$out_file=sprintf("%02d%02d%02d_%02d00_dia_rsd.nc",$year,$month,$day,$hour);

	$CURRENTLY_IN_SONDE_SECTION_OF_PROFILE=$True;
	next;
    } # end if

    next if(/BLOCK/);
    next if(/LEVEL/);
    next if(/STACKED DATA/);

    if (/WIND BY PRESSURE/){
# We've reached the end of the sonde profile, so output the file and skip to the next profile.
	$CURRENTLY_IN_SONDE_SECTION_OF_PROFILE=$False;

	@lev=@p;
	$num_lev=$#p+1;
	
# Compute some quantities involving the whole column
	@dp=@p;
	@mass_path_lvl_H2O=@p;
	for ($idx=1;$idx<=$#p;$idx++){
	    $dp[$idx]=$p[$idx]-$p[$idx-1];
	} # endfor
	$dp[0]=$dp[1];
        $mass_path_tot_H2O=0.;
	for ($idx=0;$idx<=$#p;$idx++){
	    $mass_path_lvl_H2O[$idx]=$q_H2O[$idx]*$dp[$idx]/$mean_sfc_gravity;
            $mass_path_tot_H2O+=$mass_path_lvl_H2O[$idx];
	} # endfor
	
	open(OUT_FILE,">".$out_file.".cdl") || die "cannot open $out_file.cdl";
	printf STDOUT ("Creating %s.cdl\n",$out_file);
	
	$cdl=
    sprintf("netcdf %s {\n",$out_file).
    "dimensions:\n".
    sprintf("\tlev = %d;\n",$#lev-$[+1).
    "variables:\n".
    sprintf("\t\t:cmd_line=\"%s\";\n",$cmd_line).
    sprintf("\t\t:creation_date=\"%s\";\n",$local_date_time).
    sprintf("\t\t:RCS_Header=\"%s\";\n",$RCS_Header).
    "\tint\tnum_lev;\n".
    "\t\tnum_lev:long_name=\"number of vertical layers\";\n".
    "\t\tnum_lev:units=\"number\";\n".
    "\tfloat mass_path_tot_H2O;\n".
    "\t\tmass_path_tot_H2O:long_name=\"Mass path of H2O in column\";\n".
    "\t\tmass_path_tot_H2O:units=\"kilogram meter-2\";\n".
    "\tfloat\tlev(lev);\n".
    "\t\tlev:long_name=\"pressure at layer midpoint\";\n".
    "\t\tlev:units=\"pascal\";\n".
    "\tfloat\tp(lev);\n".
    "\t\tp:long_name=\"pressure at layer midpoint\";\n".
    "\t\tp:units=\"pascal\";\n".
    "\tfloat t(lev);\n".
    "\t\tt:long_name=\"temperature at layer midpoint\";\n".
    "\t\tt:units=\"kelvin\";\n".
    "\tfloat dew_dep(lev);\n".
    "\t\tdew_dep:long_name=\"dewpoint depression\";\n".
    "\t\tdew_dep:units=\"kelvin\";\n".
    "\tfloat RH_liq(lev);\n".
    "\t\tRH_liq:long_name=\"relative humidity w/r/t liquid water\";\n".
    "\t\tRH_liq:units=\"fraction\";\n".
    "\tfloat r_H2O(lev);\n".
    "\t\tr_H2O:long_name=\"mass mixing ratio\";\n".
    "\t\tr_H2O:units=\"kilogram kilogram-1\";\n".
    "\tfloat q_H2O(lev);\n".
    "\t\tq_H2O:long_name=\"specific humidity\";\n".
    "\t\tq_H2O:units=\"kilogram kilogram-1\";\n".
    "\tfloat mass_path_lvl_H2O(lev);\n".
    "\t\tmass_path_lvl_H2O:long_name=\"Mass path of H2O in layer\";\n".
    "\t\tmass_path_lvl_H2O:units=\"kilogram meter-2\";\n".
#    "\tfloat foo(lev);\n".
#    "\t\tfoo:long_name=\"dewpoint depression\";\n".
#    "\t\tfoo:units=\"kelvin\";\n".
    "data:\n";
#    "\n".
#    sprintf("\n").

# Print out the CDL file from which the netCDF file is constructable.
	print OUT_FILE $cdl;
	print OUT_FILE "\tnum_lev\t= ".join(",",$num_lev).";\n";
	print OUT_FILE "\tmass_path_tot_H2O\t= ".join(",",$mass_path_tot_H2O).";\n";
	print OUT_FILE "\tlev\t= ".join(",",@p).";\n";
	print OUT_FILE "\tp\t= ".join(",",@p).";\n";
	print OUT_FILE "\tt\t= ".join(",",@t).";\n";
	print OUT_FILE "\tdew_dep\t= ".join(",",@dew_dep).";\n";
	print OUT_FILE "\tRH_liq\t= ".join(",",@RH_liq).";\n";
	print OUT_FILE "\tr_H2O\t= ".join(",",@r_H2O).";\n";
	print OUT_FILE "\tq_H2O\t= ".join(",",@q_H2O).";\n";
	print OUT_FILE "\tmass_path_lvl_H2O\t= ".join(",",@mass_path_lvl_H2O).";\n";
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

# Now that the data has been written out, zero-the arrays beforea the next profile.	
	undef(@p);
	undef(@t);
	undef(@dew_dep);
	undef(@r_H2O);
	undef(@q_H2O);
	undef(@mass_path_lvl_H2O);
	undef(@RH_liq);
#	undef(@foo);

    } # end if
    
    if ($CURRENTLY_IN_SONDE_SECTION_OF_PROFILE){
# Process the line.
	chop;
	@line=split(/ +/,$_);

	$p_tmp=$line[2];
	$t_tmp=$line[4];
	$dew_dep_tmp=$line[5];

# Skip over the bad data, which is usually just the mandatory pressure levels
# being blocked because the ground elevation was above 1000 mb.
	next if($t_tmp == 999.9);

# Convert the data to SI as we go.
	$p_tmp*=100.;
	if($t_tmp != 999.9){$t_tmp+=273.15;}else{$t_tmp=1.e36;}
	if($dew_dep_tmp != 99.9){$dew_dep_tmp=$dew_dep_tmp;}else{$dew_dep_tmp=1.e36;}
	if(($dew_dep_tmp < 1.e36) && ($t_tmp < 1.e36)){
	    $t_dew=$t_tmp-$dew_dep_tmp;
            $RH_liq_tmp=exp(-$l_over_R*$dew_dep_tmp/($t_tmp*$t_dew));
            $sat_vap_liq=100.*exp(21.6-5420./$t_tmp);
            $pp_H2O=$RH_liq_tmp*$sat_vap_liq;
            $r_H2O_tmp=$epsilon*$pp_H2O/($p_tmp-$pp_H2O);
            $q_H2O_tmp=$r_H2O_tmp/(1.+$r_H2O_tmp);
	}else{
            $RH_liq_tmp=1.e36;
            $r_H2O_tmp=1.e36;
            $q_H2O_tmp=1.e36;
        } # end else

# So far we are working only with radiosonde data, which has the first level
# at the surface, but we want to generate CCM RT style coordinate arrays 
# where the first level is at TOA. Therefore we use unshift() rather than
# push() to prepend the new elements to the array, rather than append them.
	unshift(@p,$p_tmp);
	unshift(@t,$t_tmp);
	unshift(@dew_dep,$dew_dep_tmp);
	unshift(@r_H2O,$r_H2O_tmp);
	unshift(@q_H2O,$q_H2O_tmp);
	unshift(@RH_liq,$RH_liq_tmp);
#	unshift(@foo,$foo_tmp);
    } # end if

} # end while

close IN_FILE;

if($debug == 1){		 
    printf STDOUT ("%10s %10s %10s %10s\n",
		   "idx","p","T","Dew Dep");
    printf STDOUT ("%10s %10s %10s %10s\n",
		   "#","pascal","kelvin","kelvin");
    for ($idx=0;$idx<=$#lev;$idx++){
	printf STDOUT ("%10d %10.4e %10.4e %10.4e\n",
		       $idx,
		       $p[$idx],
		       $t[$idx],
		       $dew_dep[$idx]
		       );
    } # endfor
} # end debug 

if($debug == 2){
    printf STDOUT ("Size of arrays = %d, last value = %d\n",$#lev-$[+1,$lev[$#lev]);
} # end debug 

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
Usage:  dss3532nc [options] in_file out_file\

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
