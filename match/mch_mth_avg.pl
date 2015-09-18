#!/usr/bin/perl

# Purpose: Compute monthly average files from h0??? files

# Example: ~/match/mch_mth_avg.pl --yyyymm_begin=199312 --yyyymm_srt=199401 --yyyymm_end=199411 dstmch90

# scp ${HOME}/match/mch_mth_avg.pl dust.ps.uci.edu:${HOME}/match

use strict;
use Getopt::Long;
use Date::Calc qw(Days_in_Month);
use Date::Calc qw(Delta_Days);
sub usage { die "usage: @ARGV$0 --yyyymm_begin=yyyymm_begin --yyyymm_srt=yyyymm_srt --yyyymm_end=yyyymm_end caseid\n"; }

# Declare local variables
my ($caseid);
my ($yyyymm_begin, $yyyymm_srt, $yyyymm_end);
my ($yyyy0, $mm0, $dd0, $yyyy, $mm, $dd);
my ($yyyymm, $day_srt, $day_end);
my ($cmd);

GetOptions(
	   'yyyymm_begin=i' => \$yyyymm_begin, # [day] Start of run
	   'yyyymm_srt=i'   => \$yyyymm_srt,   # [day] Start of averages
	   'yyyymm_end=i'   => \$yyyymm_end    # [day] End of averages
	   );
if ($#ARGV != 0) {usage();}
$caseid=$ARGV[0];

print "caseid       = ", $caseid, "\n";
print "yyyymm_begin = ", $yyyymm_begin, "\n";
print "yyyymm_srt   = ", $yyyymm_srt, "\n";
print "yyyymm_end   = ", $yyyymm_end, "\n";

$yyyy0  = int($yyyymm_begin/100);
$mm0    =     $yyyymm_begin%100;
$dd0    = 01;
$dd     = 01;

# loop
$yyyymm=$yyyymm_srt;
while ($yyyymm <= $yyyymm_end) {

    # get day boundaries
    $yyyy     = int($yyyymm/100);
    $mm       =     $yyyymm%100;
    $day_srt  = Delta_Days( $yyyy0,$mm0,$dd0,$yyyy,$mm,$dd);
    $day_end  = $day_srt + Days_in_Month( $yyyy,$mm);
    $day_srt += 1;
    #print "day_srt = ", $day_srt, "\n";
    #print "day_end = ", $day_end, "\n";

    # compute monthly average
#    $cmd="ncra -O -D 3 -d time,".$day_srt.".0,".$day_end.".0 /data/zender/ZENDER/match/".$caseid."/hist/h0???.nc ".$caseid."_".$yyyymm.".nc";
    $cmd="ncra -O -D 3 -d time,".$day_srt.".0,".$day_end.".0 /data/zender/".$caseid."/h0???.nc /data/zender/".$caseid."/".$caseid."_".$yyyymm.".nc";
    print $cmd, "\n";
    system($cmd);
    print "processed ", $yyyymm,  "\n";
    
    # increment month
    $yyyymm += 1;
    if ($mm == 12) {
	$yyyy  += 1;
	$yyyymm = 100*$yyyy+1;
    }
}

