#!/data/zender/bin/SUNMP/perl
##!/data/zender/bin/SUNMP/perl -w -d
##!/data/zender/bin/SUNMP/perl -w
##!/usr/local/bin/perl

# $Id$

# Purpose: Create seasonal cycle files of regionally averaged anomalies

# Example usage:

# Debugging
# /home/zender/aca/reg_avg.pl -d 1

# Production
# /home/zender/aca/reg_avg.pl -d 1

# $Log: not supported by cvs2svn $
# Revision 1.1.1.1  2000-07-13 04:47:51  zender
# Imported sources
#

# Specify modules
use Getopt::Long;
use NetCDF;

# Set output flushing to help debugging on hard crashes. 
# These options update the filehandle after every output statement.
# See the camel book, p. 110.
select((select(STDOUT),$|=1)[0]);
select((select(STDERR),$|=1)[0]);

# Timing information
($start_user_time,$start_system_time,$start_child_user_time,$start_child_system_time)=times;
$local_date_time=date_time();
printf STDOUT ("Start time: %s\n",$local_date_time);

# Set defaults 
$False=0;
$True=1;

$Boolean=$True;
$combine_isotopes=$False;
$HELP=$False;
$RCS_Header='$Id$';
$RCS_Revision='$Revision$';
$dbg_lvl=0;
$ln_lo=0.;
$ln_hi=50000.;
$mlc_nbr=1;
$iso_nbr=1;
$float_foo=102.543;
$pi=atan2(1,1)*4.;
				 
# Parse command line arguments 
$result=GetOptions(
		   'Boolean!',
                   'combine_isotopes!',
		   'debug=i',
		   'float_foo=f',
		   'help|usage',
                   'high_wavenumber=f',
		   'isotope=i',
                   'low_wavenumber=f',
		   'molecule=i',
		   );
if (defined($opt_Boolean)) {$Boolean=!$Boolean};
if (defined($opt_combine_isotopes)) {$combine_isotopes=!$combine_isotopes};
if (defined($opt_debug)) {$dbg_lvl=$opt_debug};
if (defined($opt_float_foo)) {$float_foo=$opt_float_foo};
if (defined($opt_help)) {$HELP=$opt_help};
if (defined($opt_high_wavenumber)) {$ln_hi=$opt_high_wavenumber};
if (defined($opt_isotope)) {$iso_nbr=$opt_isotope};
if (defined($opt_low_wavenumber)) {$ln_lo=$opt_low_wavenumber};
if (defined($opt_molecule)) {$mlc_nbr=$opt_molecule};

if($HELP){		 
    &prn_usg();
    exit 0;
} # end HELP

%case=('erbe_b','amip5','spcp_81');
@var=('LWCF','SWCF');

if($dbg_lvl >= 1){		 
    printf STDOUT ("result = %d, Boolean = %d, combine_isotopes = %d, dbg_lvl = %d, float_foo = %f, mlc_nbr = %d, iso_nbr=%d, ln_lo = %f, ln_hi = %f, HELP = %d, mlc_sng = %s, iso = %s\n",$result,$Boolean,$combine_isotopes,$dbg_lvl,$float_foo,$mlc_nbr,$iso_nbr,$ln_lo,$ln_hi,$HELP,$mlc_sng{$mlc_nbr},${$mlc_sng{$mlc_nbr}}{$iso_nbr});
} # end debug 

%reg=(
      'Africa_South',[['lat',-35.,-10.],['lon',10.,40.],['ORO','eq',1.]],
      'Alaska_NW_Canada',[['lat',50.,70.],['lon',190.,260.],['ORO','eq',1.]],
      'Amazon_Basin',[['lat',-10.,0.],['lon',290.,310.],['ORO','eq',1.]],
      'America_Central',[['lat',10.,25.],['lon',250.,280.],['ORO','eq',1.]],
      'America_South_South',[['lat',-60.,-25.],['lon',280.,310.],['ORO','eq',1.]],
      'Antarctica',[['lat',-90.,-65.],['lon',0.,360.],['ORO','eq',1.]],
      'Atlantic_North',[['lat',30.,60.],['lon',300.,330.],['ORO','eq',0.]],
      'Australia',[['lat',-40.,-10.],['lon',110.,160.],['ORO','eq',1.]],
      'Congo',[['lat',-10.,5.],['lon',10.,30.],['ORO','eq',1.]],
#      'Europe_Central',[['lat',40.,55.],['lon',-10.,40.],['ORO','eq',1.]],
      'Europe_Northern',[['lat',55.,70.],['lon',5.,60.],['ORO','eq',1.]],
      'Greenland',[['lat',60.,90.],['lon',300.,340.],['ORO','eq',1.]],
      'India',[['lat',10.,30.],['lon',90.,110.],['ORO','eq',1.]],
      'Indian_Equatorial',[['lat',-10.,10.],['lon',50.,100.],['ORO','eq',0.]],
      'Indian_Tropical',[['lat',-20.,10.],['lon',50.,100.],['ORO','eq',0.]],
      'Indochina',[['lat',10.,30.],['lon',90.,120.],['ORO','eq',1.]],
      'Indonesia_Land',[['lat',-10.,10.],['lon',90.,150.],['ORO','eq',1.]],
      'Indonesia_Ocean',[['lat',-10.,10.],['lon',90.,150.],['ORO','eq',0.]],
      'Pacific_South',[['lat',-60.,-30.],['lon',150.,270.],['ORO','eq',0.]],
      'Pacific_North',[['lat',30.,50.],['lon',150.,210.],['ORO','eq',0.]],
      'Pacific_Equatorial',[['lat',-10.,10.],['lon',140.,270.],['ORO','eq',0.]],
      'Pacific_Tropical',[['lat',-20.,20.],['lon',140.,270.],['ORO','eq',0.]],
      'Pacific_Equatorial_Central',[['lat',-10.,10.],['lon',180.,230.],['ORO','eq',0.]],
      'Pacific_Equatorial_Eastern',[['lat',-10.,10.],['lon',230.,270.],['ORO','eq',0.]],
      'Pacific_Equatorial_Western',[['lat',-10.,10.],['lon',140.,180.],['ORO','eq',0.]],
      'Pacific_Tropical_Central',[['lat',-20.,20.],['lon',180.,230.],['ORO','eq',0.]],
      'Pacific_Tropical_Eastern',[['lat',-20.,20.],['lon',230.,270.],['ORO','eq',0.]],
      'Pacific_Tropical_Western',[['lat',-20.,20.],['lon',140.,180.],['ORO','eq',0.]],
      'Pacific_Tropical_Northwest',[['lat',0.,30.],['lon',120.,180.],['ORO','eq',0.]],
      'Pacific_Western_Warm_Pool',[['lat',-10.,10.],['lon',140.,170.],['ORO','eq',0.]],
#      'Sahara_Arabian_Peninsula',[['lat',10.,30.],['lon',-20.,50.],['ORO','eq',1.]],
      'Siberia_Eastern',[['lat',50.,70.],['lon',90.,140.],['ORO','eq',1.]],
      'Siberia_Western',[['lat',50.,70.],['lon',60.,90.],['ORO','eq',1.]],
      'Tibetan_Plateau',[['lat',30.,40.],['lon',80.,100.],['ORO','eq',1.]],
      'US_Central',[['lat',30.,50.],['lon',250.,270.],['ORO','eq',1.]],
      'US_Eastern',[['lat',30.,50.],['lon',270.,290.],['ORO','eq',1.]],
      'US_Western',[['lat',30.,50.],['lon',230.,250.],['ORO','eq',1.]],
      );

foreach $case (@case){

    foreach $region (@region){

# Produce a command like the following:
# ncwa -O -a lat,lon -m ORO -o eq -M 1. -w gw -d lat,10.,30. -d lon,90.,110. -v LWCF,SWCF /data2/zender/amip5/amip5_anom_8589x87_0112.nc /data2/zender/amip5/amip5_anom_reg_India_8589x87_0112.nc
	$cmd_sng="ncwa -O -a lat,lon -m ORO -o ".$op_type." -M ".$op_val." -w gw -d ";

	if($dbg_lvl > 0){		 
	    printf STDOUT ("case=%s reg_nm=%s lon_wst=%g lon_est=%g lat_sth=%g lat_nth=%g\n",$case,$reg_nm,$lon_wst,$lon_est,$lat_sth,$lat_nth);
	} # end debug 

    } # end loop over regions

} # end loop over cases

# Timing information
($end_user_time,$end_system_time,$end_child_user_time,$end_child_system_time)=times;
printf STDOUT ("Elapsed time: User CPU %.2f s, System CPU %.2f s\n",$end_user_time-$start_user_time,$end_system_time-$start_system_time);
printf STDOUT ("Elapsed time: Child user CPU %.2f s, Child system CPU %.2f s\n",$end_child_user_time-$start_child_user_time,$end_child_system_time-$start_child_system_time);
$local_date_time=date_time();
printf STDOUT ("End time: %s\n",$local_date_time);

sub date_time{
# Returns local time string in "Sun Sep 18 20:14:49 1994" format
    local($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst);
    ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
    $date_time_string=
	sprintf
	    ('%s %s %02d %02d:%02d:%02d 19%02d',
	     (Sun,Mon,Tue,Wed,Thu,Fri,Sat)[$wday],
	     (Jan,Feb,Mar,Apr,May,Jun,Jul,Aug,Sep,Oct,Nov,Dec)[$mon],
	     $mday,
	     $hour,
	     $min,
	     $sec,
	     $year);

    return $date_time_string;
} # end sub date_time()
