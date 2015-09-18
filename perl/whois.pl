#!/contrib/bin/perl 

my $CVS_Header='$Id$';

# Purpose: Search internic for valid domain names

# Usage:
# /home/zender/perl/whois.pl
# /home/zender/perl/whois.pl >! whois.txt &

use strict; # Protect all namespaces; this produces many errors with non-expert code
use File::Basename; # Parses filenames

# Specify modules
use Getopt::Long; # GNU-style getopt
use strict; # Protect all namespaces
use lib qw(/home/zender/perl);
use DBG; # Debugging constants
require 'csz.pl'; # Personal library: date_time(), YYYYMMDD(), ...

# Set output flushing to help debugging on hard crashes. 
# These options update the filehandle after every output statement.
# See Camel book, p. 110.
select((select(STDOUT),$|=1)[0]);
select((select(STDERR),$|=1)[0]);

# Timing information
my $lcl_date_time=&time_srt();

# Declare local variables
my ($prg_nm,$prg_dsc,$prg_vrs,$prg_date);
my ($pth_in,$fl_sfx,$fl_out);
my ($rcd);

# Set defaults 
my $False=0;
my $True=1;
my $CVS_Date='$Date$';
my $CVS_Id='$Id$';
my $CVS_Revision='$Revision$';

# Set defaults for command line arguments
my $dbg_lvl=0;

$prg_dsc='whois searcher'; # Program description
($prg_nm,$prg_vrs)=$CVS_Id =~ /: (.+).pl,v ([\d.]+)/; # Program name and version
$prg_vrs.='*' if length('$Locker:  $ ') > 12; # Tack '*' if it is not checked in into CVS.
($prg_nm,$pth_in,$fl_sfx)=fileparse($0,''); # $0 is program name Camel p. 136.
($prg_date)=unpack '@7 a19',$CVS_Date if length($CVS_Date) > 7;

# Parse command line arguments: '!' means Boolean, '|' is OR, '=' specifies required argument: 'i' is integer, 'f' is float
$rcd=Getopt::Long::Configure('no_ignore_case'); # Turn on case-sensitivity
$rcd=GetOptions(
		'dbg_lvl=i' => \$dbg_lvl,
		);
if($dbg_lvl != $dbg_off){print ("$prg_nm: $prg_dsc, version $prg_vrs of $prg_date\n");} # endif dbg
# Output file name should be only positional argument
die "$prg_nm: Too many positional arguments" if($#ARGV > 0);
$fl_out=$ARGV[$#ARGV] if($#ARGV == 0);

# Main Code
my ($dmn_nm,$whois_cmd,$rx_sng);
my ($sng1,$sng2,$sng3);
my ($rgs_flg);
$rx_sng="No match for";
if(defined($fl_out)){open(FL_OUT,">".$fl_out) || die "$prg_nm: Error opening $fl_out";}
foreach $sng1 (1..10000){
#    foreach $sng2 ('a'..'z'){
#	foreach $sng3 ('a'..'z'){
	    $rgs_flg=$True;
	    $dmn_nm=$sng1.'.com';
	    $whois_cmd='whois '.$dmn_nm;
	    open(WHOIS,$whois_cmd." |") || die "Cannot execute $whois_cmd, stopped";
	    printf STDERR ("Trying %s...",$dmn_nm);
	    while(<WHOIS>){
		if(/$rx_sng/){
		    $rgs_flg=$False;
		    printf STDERR ("we have a winner!\n");
		    if(defined($fl_out)){
			printf FL_OUT ("Unregistered domain: %s\n",$dmn_nm);
		    }else{
			printf STDOUT ("Unregistered domain: %s\n",$dmn_nm);
		    } # endif
		    last;
		} # endif
	    } # endwhile
	    close WHOIS;
	    if($rgs_flg){printf STDERR ("registered\n");}
#	} # endfor $sng3
#    } # endfor $sng2
} # endfor $sng1





