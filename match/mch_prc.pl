#!/usr/bin/perl

# $Id$

# Purpose: Process data related to the MATCH experiment

# Usage:
# ${HOME}/match/mch_prc.pl
# ${HOME}/match/mch_prc.pl --yyyymmdd_srt=19990208 --yyyymmdd_end=19990322 
# ${HOME}/match/mch_prc.pl --yyyymmdd_srt=19990208 --yyyymmdd_end=19990322 ${DATA}/dstmch74/dstmch74_19990208_19990322.nc

BEGIN{
    unshift @INC,$ENV{'HOME'}.'/perl'; # Location of DBG.pm HaS98 p. 170
} # end BEGIN

my $CVS_Header='$Id$';

# Specify modules
use strict; # Protect all namespaces
use Getopt::Long; # GNU-style getopt
use File::Basename; # For parsing filenames

# 3rd party modules

# Personal modules
use DBG; # Debugging constants
require 'csz.pl'; # Personal library, &date_time(), &cmd_prc(), ...

# Declare local variables
my ($idx,$rcd,$foo);
my ($prg_nm,$prg_dsc,$prg_vrs,$prg_date);
my ($pth_in,$pth_out);
my ($fl_nm,$fl_pth,$fl_sfx);
my ($fl_in,$fl_out);
my ($drc_in,$drc_out);
my ($yr,$mth,$day);
my ($yr_srt,$mth_srt,$day_srt);
my ($yr_end,$mth_end,$day_end);
my ($day_srt_crr,$day_end_crr);
my ($mth_srt_crr,$mth_end_crr);
my ($yyyymmdd_sng);

# Set defaults 
my $False=0;
my $True=1;

my $CVS_Date='$Date$';
my $CVS_Id='$Id$';
my $CVS_Revision='$Revision$';
my $HELP=$False;
my $DATA=$ENV{'DATA'};
my $HOME=$ENV{'HOME'};
my $PVM_ARCH=$ENV{'PVM_ARCH'};
my $usr_nm=$ENV{'USER'};

# Set defaults for command line arguments
my $Boolean=$True;
my $caseid='dstmch74';
my $dbg_lvl=0; # Debugging level
my $fl_in_usr=$False; # [flg] True if user specifies fl_in
my $float_foo=102.543;
my $nco_flg=$False; # [flg] Postprocess data with NCO
my $yyyymmdd_srt=19990208; # Start date in YYYYMMDD format
my $yyyymmdd_end=19990322; # End date in YYYYMMDD format
my $mfilt=10; # Days per tape in simulation
my $yyyymmdd_sml_srt=19980101; # Start date of simulation in YYYYMMDD format
my $drc_in='/data/'+$usr_nm+'/'+'ZENDER'+'/match/'+$caseid+'/hist' # Location of input files

# Derived fields

$prg_dsc='MATCH to netCDF converter'; # Program description
($prg_nm,$prg_vrs)=$CVS_Id =~ /: (.+).pl,v ([\d.]+)/; # Program name and version
$prg_vrs.='*' if length('$Locker:  $ ') > 12; # Tack '*' if it is not checked in into CVS.
($prg_nm,$pth_in,$fl_sfx)=fileparse($0,''); # $0 is program name Camel p. 136.
if(length($CVS_Date) > 6){($prg_date)=unpack '@7 a19',$CVS_Date;}else{$prg_date='Unknown';}

# Parse command line arguments: '!' means Boolean, '|' is OR, '=' specifies required argument: 'i' is integer, 'f' is float, 's' is string
$rcd=Getopt::Long::Configure('no_ignore_case'); # Turn on case-sensitivity
$rcd=GetOptions( # man Getopt::GetoptLong
		 'Boolean!' => \$Boolean,
		 'caseid=s' => \$caseid,
		 'dbg_lvl=i' => \$dbg_lvl,
		 'yyyymmdd_end=i', => \$yyyymmdd_end, # [day] Maximum day number
		 'yyyymmdd_srt=i' => \$yyyymmdd_srt, # [day] Minimum day number
		 'float_foo=f' => \$float_foo, 
		 'help|usage', => \$HELP,
		 'mfilt', => $mfilt, # Days per tape in simulation
		 'yyyymmdd_sml_srt', => $yyyymmdd_sml_srt, # Start date of simulation in YYYYMMDD format
		 'nco_flg!' => \$nco_flg, # [flg] Postprocess data with NCO
		 ); # end GetOptions arguments

if($HELP){&usg_prn();exit 0;} # end HELP

# Parse positional arguments, if any
if($#ARGV+1 > 1){die "$prg_nm: ERROR Called with $#ARGV+1 positional arguments, need no more than 1\n";}
elsif($#ARGV+1 == 2){
# Input file name is first positional argument, if any
# Output file name is last positional argument, if any
    $fl_in=$ARGV[0];
    $fl_out=$ARGV[1];
    $fl_in_usr=$True;
}elsif($#ARGV+1 == 1){
    $fl_in=$ARGV[0];
    $fl_in_usr=$True;
}elsif($#ARGV+1 == 0){
    print STDOUT "$prg_nm: INFO Using default filenames\n";
} # end else

# Definitions that depend on command line input
$drc_in="$DATA/$caseid"; # Directory for input files
$drc_out="$DATA/$caseid"; # Directory for output files

if(!$fl_in_usr){
# Create input filename
    $fl_in=$drc_in.'/'.$caseid.'/'.$caseid.'_'.$yyyymmdd_srt.'_'.$yyyymmdd_end.'.nc';
} # endif

($yr_srt,$mth_srt,$day_srt)=&yyyymmdd_prs($yyyymmdd_srt);
($yr_end,$mth_end,$day_end)=&yyyymmdd_prs($yyyymmdd_end);

# Print initialization state
if($dbg_lvl > 1){print ("$prg_nm: $prg_dsc, version $prg_vrs of $prg_date\n");} # endif dbg
if($dbg_lvl >= 2){print ("$prg_nm: \$dbg_lvl = $dbg_lvl\n");} # endif dbg
if($dbg_lvl >= 2){print ("$prg_nm: \$caseid = $caseid\n");} # endif dbg
if($dbg_lvl >= 2){print ("$prg_nm: \$yyyymmdd_end = $yyyymmdd_end\n");} # endif dbg
if($dbg_lvl >= 2){print ("$prg_nm: \$yyyymmdd_srt = $yyyymmdd_srt\n");} # endif dbg
if($dbg_lvl >= 2){print ("$prg_nm: \$fl_in = $fl_in\n");} # endif dbg
if($dbg_lvl >= 2){print ("$prg_nm: \$yr_srt = $yr_srt\n");} # endif dbg
if($dbg_lvl >= 2){print ("$prg_nm: \$mth_srt = $mth_srt\n");} # endif dbg
if($dbg_lvl >= 2){print ("$prg_nm: \$day_srt = $day_srt\n");} # endif dbg
if($dbg_lvl >= 2){print ("$prg_nm: \$yr_end = $yr_end\n");} # endif dbg
if($dbg_lvl >= 2){print ("$prg_nm: \$mth_end = $mth_end\n");} # endif dbg
if($dbg_lvl >= 2){print ("$prg_nm: \$day_end = $day_end\n");} # endif dbg

&cmd_prc("ncks -H -F -u -C -v time $fl_in");

foreach $yr ($yr_srt..$yr_end){
    if($yr == $yr_srt){$mth_srt_crr=$mth_srt;}else{$mth_srt_crr=1}
    if($yr == $yr_end){$mth_end_crr=$mth_end;}else{$mth_end_crr=12}
    foreach $mth ($mth_srt_crr..$mth_end_crr){
	if($yr == $yr_srt && $mth == $mth_srt){$day_srt_crr=$day_srt;}else{$day_srt_crr=1}
	if($yr == $yr_end && $mth == $mth_end){$day_end_crr=$day_end;}else{$day_end_crr=&dpmyr($mth,$yr);}
	foreach $day ($day_srt_crr..$day_end_crr){

	    $yyyymmdd_sng=sprintf('%04d%02d%02d',$yr,$mth,$day);

	    &cmd_prc("ncks -O -F -d time,$day $fl_in ${DATA}/${caseid}/${caseid}_${yyyymmdd_sng}.nc");

# C-indices of month boundaries, inclusive
# J   F   M    A    M    J    J    A    S    O    N    D
# 0, 31, 59,  90, 120, 151, 181, 212, 243, 273, 304, 334
#30, 58, 89, 119, 150, 180, 211, 242, 272, 303, 333, 364
my $doy_mth_srt_C=()
my $doy_mth_end_C=

# Fortran-indices of month boundaries, inclusive
# J   F   M    A    M    J    J    A    S    O    N    D
# 1, 32, 60,  91, 121, 152, 182, 213, 244, 274, 305, 335
#31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365

# Quicker to average on krein and then copy to dust
# Following time coordinates work when simulation starts on January 1
export caseid='dstmch77'
export yyyy='1994'
ncra -O -D 3 -d time,001.0,031.0 /data/zender/ZENDER/match/${caseid}/hist/h00[0-3]?.nc ${caseid}_${yyyy}01.nc &
ncra -O -D 3 -d time,032.0,059.0 /data/zender/ZENDER/match/${caseid}/hist/h00[0-3]?.nc ${caseid}_${yyyy}02.nc &
ncra -O -D 3 -d time,060.0,090.0 /data/zender/ZENDER/match/${caseid}/hist/h00[0-3]?.nc ${caseid}_${yyyy}03.nc &
ncra -O -D 3 -d time,091.0,120.0 /data/zender/ZENDER/match/${caseid}/hist/h00[0-3]?.nc ${caseid}_${yyyy}04.nc &
ncra -O -D 3 -d time,121.0,151.0 /data/zender/ZENDER/match/${caseid}/hist/h00[0-3]?.nc ${caseid}_${yyyy}05.nc &
ncra -O -D 3 -d time,152.0,181.0 /data/zender/ZENDER/match/${caseid}/hist/h00[0-3]?.nc ${caseid}_${yyyy}06.nc &
ncra -O -D 3 -d time,182.0,212.0 /data/zender/ZENDER/match/${caseid}/hist/h00[0-3]?.nc ${caseid}_${yyyy}07.nc &
ncra -O -D 3 -d time,213.0,243.0 /data/zender/ZENDER/match/${caseid}/hist/h00[0-3]?.nc ${caseid}_${yyyy}08.nc &
ncra -O -D 3 -d time,244.0,273.0 /data/zender/ZENDER/match/${caseid}/hist/h00[0-3]?.nc ${caseid}_${yyyy}09.nc &
ncra -O -D 3 -d time,274.0,304.0 /data/zender/ZENDER/match/${caseid}/hist/h00[0-3]?.nc ${caseid}_${yyyy}10.nc &
ncra -O -D 3 -d time,305.0,334.0 /data/zender/ZENDER/match/${caseid}/hist/h00[0-3]?.nc ${caseid}_${yyyy}11.nc &
ncra -O -D 3 -d time,335.0,365.0 /data/zender/ZENDER/match/${caseid}/hist/h00[0-3]?.nc ${caseid}_${yyyy}12.nc &

export yyyy='1995'
ncra -O -D 3 -d time,366.0,396.0 /data/zender/ZENDER/match/${caseid}/hist/h00[3-7]?.nc ${caseid}_${yyyy}01.nc &
ncra -O -D 3 -d time,397.0,424.0 /data/zender/ZENDER/match/${caseid}/hist/h00[3-7]?.nc ${caseid}_${yyyy}02.nc &
ncra -O -D 3 -d time,425.0,455.0 /data/zender/ZENDER/match/${caseid}/hist/h00[3-7]?.nc ${caseid}_${yyyy}03.nc &
ncra -O -D 3 -d time,456.0,485.0 /data/zender/ZENDER/match/${caseid}/hist/h00[3-7]?.nc ${caseid}_${yyyy}04.nc &
ncra -O -D 3 -d time,486.0,516.0 /data/zender/ZENDER/match/${caseid}/hist/h00[3-7]?.nc ${caseid}_${yyyy}05.nc &
ncra -O -D 3 -d time,517.0,546.0 /data/zender/ZENDER/match/${caseid}/hist/h00[3-7]?.nc ${caseid}_${yyyy}06.nc &
ncra -O -D 3 -d time,547.0,577.0 /data/zender/ZENDER/match/${caseid}/hist/h00[3-7]?.nc ${caseid}_${yyyy}07.nc &
ncra -O -D 3 -d time,578.0,608.0 /data/zender/ZENDER/match/${caseid}/hist/h00[3-7]?.nc ${caseid}_${yyyy}08.nc &
ncra -O -D 3 -d time,609.0,638.0 /data/zender/ZENDER/match/${caseid}/hist/h00[3-7]?.nc ${caseid}_${yyyy}09.nc &
ncra -O -D 3 -d time,639.0,669.0 /data/zender/ZENDER/match/${caseid}/hist/h00[3-7]?.nc ${caseid}_${yyyy}10.nc &
ncra -O -D 3 -d time,670.0,699.0 /data/zender/ZENDER/match/${caseid}/hist/h00[3-7]?.nc ${caseid}_${yyyy}11.nc &
ncra -O -D 3 -d time,700.0,730.0 /data/zender/ZENDER/match/${caseid}/hist/h00[3-7]?.nc ${caseid}_${yyyy}12.nc &

export yyyy='1996'
ncra -O -D 3 -d time,731.0,761.0 /data/zender/ZENDER/match/${caseid}/hist/h0???.nc ${caseid}_${yyyy}01.nc &
ncra -O -D 3 -d time,762.0,790.0 /data/zender/ZENDER/match/${caseid}/hist/h0???.nc ${caseid}_${yyyy}02.nc &
ncra -O -D 3 -d time,791.0,821.0 /data/zender/ZENDER/match/${caseid}/hist/h0???.nc ${caseid}_${yyyy}03.nc &
ncra -O -D 3 -d time,822.0,851.0 /data/zender/ZENDER/match/${caseid}/hist/h0???.nc ${caseid}_${yyyy}04.nc &
ncra -O -D 3 -d time,852.0,882.0 /data/zender/ZENDER/match/${caseid}/hist/h0???.nc ${caseid}_${yyyy}05.nc &
ncra -O -D 3 -d time,883.0,912.0 /data/zender/ZENDER/match/${caseid}/hist/h0???.nc ${caseid}_${yyyy}06.nc &
ncra -O -D 3 -d time,913.0,943.0 /data/zender/ZENDER/match/${caseid}/hist/h0???.nc ${caseid}_${yyyy}07.nc &
ncra -O -D 3 -d time,944.0,974.0 /data/zender/ZENDER/match/${caseid}/hist/h0???.nc ${caseid}_${yyyy}08.nc &
ncra -O -D 3 -d time,975.0,1004.0 /data/zender/ZENDER/match/${caseid}/hist/h0???.nc ${caseid}_${yyyy}09.nc &
ncra -O -D 3 -d time,1005.0,1035.0 /data/zender/ZENDER/match/${caseid}/hist/h0???.nc ${caseid}_${yyyy}10.nc &
ncra -O -D 3 -d time,1036.0,1065.0 /data/zender/ZENDER/match/${caseid}/hist/h0???.nc ${caseid}_${yyyy}11.nc &
ncra -O -D 3 -d time,1066.0,1096.0 /data/zender/ZENDER/match/${caseid}/hist/h0???.nc ${caseid}_${yyyy}12.nc &

export yyyy='1997'
ncra -O -D 3 -d time,1097.0,1127.0 /data/zender/ZENDER/match/${caseid}/hist/h0???.nc ${caseid}_${yyyy}01.nc &
ncra -O -D 3 -d time,1128.0,1155.0 /data/zender/ZENDER/match/${caseid}/hist/h0???.nc ${caseid}_${yyyy}02.nc &
ncra -O -D 3 -d time,1156.0,1186.0 /data/zender/ZENDER/match/${caseid}/hist/h0???.nc ${caseid}_${yyyy}03.nc &
ncra -O -D 3 -d time,1187.0,1216.0 /data/zender/ZENDER/match/${caseid}/hist/h0???.nc ${caseid}_${yyyy}04.nc &
ncra -O -D 3 -d time,1217.0,1247.0 /data/zender/ZENDER/match/${caseid}/hist/h0???.nc ${caseid}_${yyyy}05.nc &
ncra -O -D 3 -d time,1248.0,1277.0 /data/zender/ZENDER/match/${caseid}/hist/h0???.nc ${caseid}_${yyyy}06.nc &
ncra -O -D 3 -d time,1278.0,1308.0 /data/zender/ZENDER/match/${caseid}/hist/h0???.nc ${caseid}_${yyyy}07.nc &
ncra -O -D 3 -d time,1309.0,1339.0 /data/zender/ZENDER/match/${caseid}/hist/h0???.nc ${caseid}_${yyyy}08.nc &
ncra -O -D 3 -d time,1340.0,1369.0 /data/zender/ZENDER/match/${caseid}/hist/h0???.nc ${caseid}_${yyyy}09.nc &
ncra -O -D 3 -d time,1370.0,1400.0 /data/zender/ZENDER/match/${caseid}/hist/h0???.nc ${caseid}_${yyyy}10.nc &
ncra -O -D 3 -d time,1401.0,1430.0 /data/zender/ZENDER/match/${caseid}/hist/h0???.nc ${caseid}_${yyyy}11.nc &
ncra -O -D 3 -d time,1431.0,1461.0 /data/zender/ZENDER/match/${caseid}/hist/h0???.nc ${caseid}_${yyyy}12.nc &

export yyyy='1998'
ncra -O -D 3 -d time,1462.0,1492.0 /data/zender/ZENDER/match/${caseid}/hist/h0???.nc ${caseid}_${yyyy}01.nc &
ncra -O -D 3 -d time,1493.0,1520.0 /data/zender/ZENDER/match/${caseid}/hist/h0???.nc ${caseid}_${yyyy}02.nc &
ncra -O -D 3 -d time,1521.0,1551.0 /data/zender/ZENDER/match/${caseid}/hist/h0???.nc ${caseid}_${yyyy}03.nc &
ncra -O -D 3 -d time,1552.0,1581.0 /data/zender/ZENDER/match/${caseid}/hist/h0???.nc ${caseid}_${yyyy}04.nc &
ncra -O -D 3 -d time,1582.0,1612.0 /data/zender/ZENDER/match/${caseid}/hist/h0???.nc ${caseid}_${yyyy}05.nc &
ncra -O -D 3 -d time,1613.0,1642.0 /data/zender/ZENDER/match/${caseid}/hist/h0???.nc ${caseid}_${yyyy}06.nc &
ncra -O -D 3 -d time,1643.0,1673.0 /data/zender/ZENDER/match/${caseid}/hist/h0???.nc ${caseid}_${yyyy}07.nc &
ncra -O -D 3 -d time,1674.0,1704.0 /data/zender/ZENDER/match/${caseid}/hist/h0???.nc ${caseid}_${yyyy}08.nc &
ncra -O -D 3 -d time,1705.0,1734.0 /data/zender/ZENDER/match/${caseid}/hist/h0???.nc ${caseid}_${yyyy}09.nc &
ncra -O -D 3 -d time,1735.0,1765.0 /data/zender/ZENDER/match/${caseid}/hist/h0???.nc ${caseid}_${yyyy}10.nc &
ncra -O -D 3 -d time,1766.0,1795.0 /data/zender/ZENDER/match/${caseid}/hist/h0???.nc ${caseid}_${yyyy}11.nc &
ncra -O -D 3 -d time,1796.0,1826.0 /data/zender/ZENDER/match/${caseid}/hist/h0???.nc ${caseid}_${yyyy}12.nc &

# Following time coordinates work when simulation starts on May 1

export caseid='dstmch25'
export date_sng='19980101_19981231'
export date_sng='19950101_19970729'
# On krein:
ncrcat -O -D 5 -d lev,1000.0, -v date,time,gw,DSTQ,DSTODXC,ORO /data/zender/ZENDER/match/${caseid}/hist/h00??.nc /data/zender/ZENDER/match/${caseid}/hist/${caseid}_${date_sng}.nc
# Move files to dust:
mv /data/zender/ZENDER/match/${caseid}/hist/${caseid}_1998??.nc /data/zender/${caseid}
mv /data/zender/ZENDER/match/${caseid}/hist/${caseid}_${date_sng}.nc /data/zender/${caseid}/${caseid}_${date_sng}.nc
# Append time coordinate from Prospero files to concentration timeseries
ncrename -v time,doy /data/zender/${caseid}/${caseid}_${date_sng}.nc
ncks -C -A -v time /data/zender/prspr/prspr_${date_sng}_xy_rgn_Brb.nc /data/zender/${caseid}/${caseid}_${date_sng}.nc
# Now run region hyperslabber
${HOME}/dst/rgn.sh dly

	} # endfor $day
    } # endfor $mth
} # endfor $yr
