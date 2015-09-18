#!/usr/bin/perl

# Purpose: Get NCEP/DSS files from NCAR
#
# Example: get_grb.pl yyyymmdd_srt yyyymmdd_end
# Example: get_grb.pl 19931201 19941231
#   yyyymmdd_srt   = start date
#   yyyymmdd_end   = end date
#

use strict;
use Getopt::Long;
sub usage { die "usage: @ARGV$0 yyyymmdd_srt yyyymmdd_end\n"; }

# Declare local variables
my ($yyyymmdd_srt, $yyyymmdd_end);
my ($yyyymmdd);
my ($cmd);
my ($file, $size, $yyyymmdd_string);
my ($file_last, $size_last, $yyyymmdd_last, $got_first);

GetOptions(
	   'yyyymmdd_srt=i'   => \$yyyymmdd_srt,   # [day] Start date
	   'yyyymmdd_end=i'   => \$yyyymmdd_end    # [day] End date
	   );
if ($#ARGV != 1) {usage();}
$yyyymmdd_srt=$ARGV[0];
$yyyymmdd_end=$ARGV[1];
$got_first=0;

print "Getting grib files:\n";
print "  yyyymmdd_srt = ", $yyyymmdd_srt, "\n";
print "  yyyymmdd_end = ", $yyyymmdd_end, "\n";

# fxm csz
get_grb_file("$ENV{\"HOME\"}/match_dst/readers/ncep/data/grb2d.list");
get_grb_file("$ENV{\"HOME\"}/match_dst/readers/ncep/data/grbsanl.list");
#get_grb_file("$ENV{\"HOME\"}/match_dst/readers/ncep/data/grb2d.list.dat_old_fmt_new");
#get_grb_file("$ENV{\"HOME\"}/match_dst/readers/ncep/data/grbsanl.list.dat_old_fmt_new");

exit;

#--------------------------------------------------------------
# subroutine to get files from list
#--------------------------------------------------------------
sub get_grb_file ($) {
    my $file_list = shift;
    my $link;
    print "file_list = ", $file_list, "\n";
    open(FP,$file_list);
    if (!<FP>) {die "ERROR: couldn't open $file_list";}
    if ($file_list =~ /grb2d/)   { $link = ".grb2d"; }
    if ($file_list =~ /grbsanl/) { $link = ".grbsanl"; }
    $got_first=0;
    while (<FP>) {
	$file_last     = $file;
	$size_last     = $size;
	$yyyymmdd_last = $yyyymmdd;
	($file, $size, $yyyymmdd_string) = split / /,$_;
	($yyyymmdd) = ($yyyymmdd_string =~ /(.{8})/);
	if ($yyyymmdd >= $yyyymmdd_srt) {
	    if ($got_first == 0) {
		print "...getting ", $file_last, " ", $yyyymmdd_last, "\n";
		$cmd="msread -f BI -R /data/zender/DSS/".$file_last." /DSS/".$file_last;
		print $cmd, "\n";
		system($cmd);
		$cmd="ln -s /data/zender/DSS/".$file_last." /data/zender/DSS/".$yyyymmdd_last.$link;
		print $cmd, "\n";
		system($cmd);
		$got_first = 1;
	    }
	    while ($yyyymmdd < $yyyymmdd_end) {
		print "...getting ", $file, " ", $yyyymmdd, "\n";	    
		$cmd="msread -f BI -R /data/zender/DSS/".$file." /DSS/".$file;
		print $cmd, "\n";
		system($cmd);
		$cmd="ln -s /data/zender/DSS/".$file." /data/zender/DSS/".$yyyymmdd.$link;
		print $cmd, "\n";
		system($cmd);
		$_ = <FP>;
		($file, $size, $yyyymmdd_string) = split / /,$_;
		($yyyymmdd) = ($yyyymmdd_string =~ /(.{8})/);
	    }
	}
    }
    close(FP);
}
