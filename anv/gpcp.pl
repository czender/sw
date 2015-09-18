#!/contrib/bin/perl
##!/usr/local/bin/perl 
				
# Purpose: Process the GPCP data

# Example usage: /home/zender/anv/gpcp.pl

foreach $case ("gpcp") {
#    foreach $yr (87..93) {
#	foreach $mth (07) {
    foreach $yr (88..94) {
	foreach $mth (01) {
		$case_sng=sprintf("%s",$case);
		$yr_sng=sprintf("%02d",$yr);
		$mth_sng=sprintf("%02d",$mth);

# Monthly: 
		&cmd_prc("ncks -O -D 3 -d time,".$yr_sng.$mth_sng."."." /data2/zender/".$case_sng."/".$case_sng."_8707_9406.nc /data2/zender/".$case_sng."/".$case_sng."_".$yr_sng."_".$mth_sng.".nc"); 
	} # endfor $mth
    } # endfor $yr
} # endfor $case

#ncra -O -D 3 gpcp_??_07.nc gpcp_8793_07.nc
#ncwa -O -D 3 -a lon gpcp_8793_07.nc gpcp_xavg_8793_07.nc
#ncra -O -D 3 gpcp_??_01.nc gpcp_8894_01.nc
#ncwa -O -D 3 -a lon gpcp_8894_01.nc gpcp_xavg_8894_01.nc

sub cmd_prc{
    my($cmd)=@_;
    print STDOUT $cmd."\n";
    $status=open(CMD,$cmd." |");
    if(!$status){print STDOUT "Unable to open() on ".$cmd;}
    while(<CMD>){print STDOUT $_;}
    close CMD
    } # end cmd_prc()

