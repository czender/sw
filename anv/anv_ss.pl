#!/contrib/bin/perl
##!/data/zender/bin/SUNMP/perl
##!/usr/local/bin/perl 

# Purpose: reduces results from cloud sensitivity studies into a form usable
# by graphics routines in spcp.pro

# Specify modules
use Env;

# Set output flushing to help debugging on hard crashes. 
# These options update the filehandle after every output statement.
# See the Camel book, p. 110.
select((select(STDOUT),$|=1)[0]);
select((select(STDERR),$|=1)[0]);

#%ENV{'PATH'}=$PATH;

# There are currently 3 sensitivity studies: updraft, sedimentation, and temperature
# Activate the data reduction you want by turning on the flag for that sensitivity study.
$updraft_ss=1;
$temperature_ss=1;
$decay_ss=0;

if($decay_ss){
#-j resets the courant safety factor for ice crystal transport
#-U turns off VAPOR_ADVECT
#cld -j 1. -D 55 -G 5000 -g 13000 -H 15000 -h 2000 -k 10. -l 130 -M 1 -m 3.e-6 -N 1.e6 -n 2160 -o cld_15_00_100_decay.nc -p 60 -r 0 -R -S 1.00 -s .40 -u -U -w 0.0 -X -x 50
    foreach $hgt_base (5000,10000,15000) {
	$fl_nm_arr="";
	$hgt_top=$hgt_base+2000.;
	for ($snapshot=0;$snapshot<=36;$snapshot+=1) {
	    $fl_nm_tmp=sprintf("/data2/zender/cld/cld_%02d_00_100_decay.nc_%02d_foo",$hgt_base/1000.,$snapshot);
	    &cmd_prc(sprintf("ncwa -O -N -w dz -v IWC_snapshot -a altitude -d altitude,%05d.,%05d. -d num_framep1,%02d /data2/zender/cld/cld_%02d_00_100_decay.nc $fl_nm_tmp",$hgt_base,$hgt_base+2000.,$snapshot,$hgt_base/1000.,$fl_nm_tmp));
	    $fl_nm_arr=$fl_nm_arr." ".$fl_nm_tmp;
	} # endfor
	&cmd_prc(sprintf("ncecat -O -v IWC_snapshot %s /data2/zender/cld/cld_%02d_00_100_decay_ss.nc",$fl_nm_arr,$hgt_base/1000.));
	&cmd_prc(sprintf("ncrename -O -v IWC_snapshot,IWP_cld /data2/zender/cld/cld_%02d_00_100_decay_ss.nc",$hgt_base/1000.));
    } # endfor

    `/bin/rm /data2/zender/cld/*foo*`;

} # $decay_ss

if($updraft_ss){
    $hgt_base_min=5000;
    $hgt_base_max=15000;
    $hgt_base_incr=5000;
    $w_min=.0;
    $w_max=.20;
    $w_incr=.02;
    &cmd_prc(sprintf("/bin/rm /data2/zender/cld/cld_??_w_???_ss.nc"));

    for ($hgt_base=$hgt_base_min;$hgt_base<=$hgt_base_max;$hgt_base+=$hgt_base_incr) {
	$hgt_top=$hgt_base+2000.;
	$hgt_mid=$hgt_base+1000.;

	# Assemble the final IWP
	$fl_nm_arr="";
	$hgt_top=$hgt_base+2000.;
	for ($w=$w_min;$w<=$w_max;$w+=$w_incr) {
	    $fl_nm_tmp=sprintf("/data2/zender/cld/cld_%02d_%02d_100.nc_foo",$hgt_base/1000.,$w*100.);
	    &cmd_prc(sprintf("ncwa -O -N -w dz -v IWC_snapshot -a altitude -d altitude,%05d.,%05d. -d num_framep1,60 /data2/zender/cld/cld_%02d_%02d_100.nc %s",$hgt_base,$hgt_base+2000.,$hgt_base/1000.,$w*100.,$fl_nm_tmp));
	    $fl_nm_arr=$fl_nm_arr." ".$fl_nm_tmp;
	} # endfor
	&cmd_prc(sprintf("ncecat -O -v IWC_snapshot %s /data2/zender/cld/cld_%02d_w_100_ss.nc",$fl_nm_arr,$hgt_base/1000.));
	&cmd_prc(sprintf("ncrename -O -v IWC_snapshot,IWP_end /data2/zender/cld/cld_%02d_w_100_ss.nc",$hgt_base/1000.));

	# Append the initial IWP
        &cmd_prc(sprintf("ncwa -O -N -w dz -v IWC_snapshot -a altitude -d altitude,%05d.,%05d. -d num_framep1,0 /data2/zender/cld/cld_%02d_00_100.nc /data2/zender/cld/foo.nc",$hgt_base,$hgt_top,$hgt_base/1000.));
        &cmd_prc(sprintf("ncrename -O -v IWC_snapshot,IWP_srt /data2/zender/cld/foo.nc"));
        &cmd_prc(sprintf("ncwa -A -a num_framep1 -v IWP_srt /data2/zender/cld/foo.nc /data2/zender/cld/cld_%02d_w_100_ss.nc",$hgt_base/1000.));

	# Append the initial IWC
        &cmd_prc(sprintf("ncwa -O -v IWC -a altitude -d altitude,%05d. /data2/zender/cld/cld_%02d_00_100.nc /data2/zender/cld/foo.nc",$hgt_mid,$hgt_base/1000.));
        &cmd_prc(sprintf("ncrename -O -v IWC,IWC_srt /data2/zender/cld/foo.nc"));
        &cmd_prc(sprintf("ncks -A -v IWC_srt /data2/zender/cld/foo.nc /data2/zender/cld/cld_%02d_w_100_ss.nc",$hgt_base/1000.));

	# Assemble the updraft,temperature,density,pressure
	$fl_nm_arr="";
	for ($w=$w_min;$w<=$w_max;$w+=$w_incr) {
	    $fl_nm_tmp=sprintf("/data2/zender/cld/cld_%02d_%02d_100.nc_foo",$hgt_base/1000.,$w*100.);
	    &cmd_prc(sprintf("ncwa -O -C -a altitude -v wind_speed,env_temperature,orig_env_temp,env_density,env_pressure,mmr_vapor,orig_mmr_vapor,pp_vapor,saturation_ice,vapor_density,IWC -d altitude,%05d. /data2/zender/cld/cld_%02d_%02d_100.nc %s",$hgt_mid,$hgt_base/1000.,$w*100.,$fl_nm_tmp));
	    $fl_nm_arr=$fl_nm_arr." ".$fl_nm_tmp;
	} # endfor
	&cmd_prc(sprintf("ncecat -O -C %s /data2/zender/cld/foo.nc",$fl_nm_arr));
	&cmd_prc(sprintf("ncrename -O -v env_temperature,temperature_end -v orig_env_temp,temperature_srt -v env_density,density -v env_pressure,pressure -v mmr_vapor,mmr_vapor_end -v orig_mmr_vapor,mmr_vapor_srt -v IWC,IWC_end /data2/zender/cld/foo.nc"));
	&cmd_prc(sprintf("ncks -A -C /data2/zender/cld/foo.nc /data2/zender/cld/cld_%02d_w_100_ss.nc",$hgt_base/1000.));
	&cmd_prc(sprintf("/bin/rm /data2/zender/cld/*foo*"));
	
    } # endfor

} # $updraft_ss

if($temperature_ss){
    $hgt_base_min=5000;
    $hgt_base_max=16000;
    $hgt_base_incr=1000;
    $w_min=.0;
    $w_max=.10;
    $w_incr=.02;
    &cmd_prc(sprintf("/bin/rm /data2/zender/cld/cld_z_??_???_ss.nc"));

    for ($w=$w_min;$w<=$w_max;$w+=$w_incr) {

	# Assemble the final IWP
	$fl_nm_arr="";
	for ($hgt_base=$hgt_base_min;$hgt_base<=$hgt_base_max;$hgt_base+=$hgt_base_incr) {
	    $hgt_top=$hgt_base+2000.;
	    $hgt_mid=$hgt_base+1000.;
	    $fl_nm_tmp=sprintf("/data2/zender/cld/cld_%02d_%02d_100.nc_foo",$hgt_base/1000.,$w*100.);
	    &cmd_prc(sprintf("ncwa -O -C -N -w dz -v IWC_snapshot -a altitude,num_framep1 -d altitude,%05d.,%05d. -d num_framep1,60 /data2/zender/cld/cld_%02d_%02d_100.nc %s",$hgt_base,$hgt_base+2000.,$hgt_base/1000.,$w*100.,$fl_nm_tmp));
	    $fl_nm_arr=$fl_nm_arr." ".$fl_nm_tmp;
	} # endfor
	&cmd_prc(sprintf("ncecat -O -C -v IWC_snapshot %s /data2/zender/cld/cld_z_%02d_100_ss.nc",$fl_nm_arr,$w*100.));
	&cmd_prc(sprintf("ncrename -O -v IWC_snapshot,IWP_end /data2/zender/cld/cld_z_%02d_100_ss.nc",$w*100.));

	# Assemble the final updraft,temperature,density,pressure
	$fl_nm_arr="";
	for ($hgt_base=$hgt_base_min;$hgt_base<=$hgt_base_max;$hgt_base+=$hgt_base_incr) {
	    $hgt_top=$hgt_base+2000.;
	    $hgt_mid=$hgt_base+1000.;
	    $fl_nm_tmp=sprintf("/data2/zender/cld/cld_%02d_%02d_100.nc_foo",$hgt_base/1000.,$w_min*100.);
	    &cmd_prc(sprintf("ncwa -O -C -a altitude -v wind_speed,env_temperature,orig_env_temp,env_density,env_pressure,mmr_vapor,orig_mmr_vapor,pp_vapor,saturation_ice,vapor_density,IWC -d altitude,%05d. /data2/zender/cld/cld_%02d_%02d_100.nc %s",$hgt_mid,$hgt_base/1000.,$w*100.,$fl_nm_tmp));
	    $fl_nm_arr=$fl_nm_arr." ".$fl_nm_tmp;
	} # endfor
	&cmd_prc(sprintf("ncecat -O -C %s /data2/zender/cld/foo.nc",$fl_nm_arr));
	&cmd_prc(sprintf("ncrename -O -v env_temperature,temperature_end -v orig_env_temp,temperature_srt -v env_density,density -v env_pressure,pressure -v mmr_vapor,mmr_vapor_end -v orig_mmr_vapor,mmr_vapor_srt -v IWC,IWC_end /data2/zender/cld/foo.nc"));
	&cmd_prc(sprintf("ncks -A /data2/zender/cld/foo.nc /data2/zender/cld/cld_z_%02d_100_ss.nc",$w*100.));
    
    } # endfor
    &cmd_prc(sprintf("/bin/rm /data2/zender/cld/*foo*"));

    # Assemble the initial IWP
    $fl_nm_arr="";
    for ($hgt_base=$hgt_base_min;$hgt_base<=$hgt_base_max;$hgt_base+=$hgt_base_incr) {
	$hgt_top=$hgt_base+2000.;
	$hgt_mid=$hgt_base+1000.;
	$fl_nm_tmp=sprintf("/data2/zender/cld/cld_%02d_%02d_100.nc_foo",$hgt_base/1000.,$w_min*100.);
	&cmd_prc(sprintf("ncwa -O -C -N -w dz -v IWC_snapshot -a altitude -d altitude,%05d.,%05d. -d num_framep1,0 /data2/zender/cld/cld_%02d_%02d_100.nc %s",$hgt_base,$hgt_top,$hgt_base/1000.,$w_min*100.,$fl_nm_tmp));
	$fl_nm_arr=$fl_nm_arr." ".$fl_nm_tmp;
    } # endfor
    &cmd_prc(sprintf("ncecat -O -C -v IWC_snapshot %s /data2/zender/cld/foo.nc",$fl_nm_arr));
    &cmd_prc(sprintf("ncrename -O -v IWC_snapshot,IWP_srt /data2/zender/cld/foo.nc"));
    
    # Append the initial IWP to all the updraft files
    for ($w=$w_min;$w<=$w_max;$w+=$w_incr) {
	&cmd_prc(sprintf("ncwa -A -C -a num_framep1 -v IWP_srt /data2/zender/cld/foo.nc /data2/zender/cld/cld_z_%02d_100_ss.nc",$w*100.));
    } # endfor
    &cmd_prc(sprintf("/bin/rm /data2/zender/cld/*foo*"));

    # Assemble the initial IWC
    $fl_nm_arr="";
    for ($hgt_base=$hgt_base_min;$hgt_base<=$hgt_base_max;$hgt_base+=$hgt_base_incr) {
	$hgt_top=$hgt_base+2000.;
	$hgt_mid=$hgt_base+1000.;
	$fl_nm_tmp=sprintf("/data2/zender/cld/cld_%02d_%02d_100.nc_foo",$hgt_base/1000.,$w_min*100.);
	&cmd_prc(sprintf("ncwa -O -C -v IWC -a altitude -d altitude,%05d. /data2/zender/cld/cld_%02d_%02d_100.nc %s",$hgt_mid,$hgt_base/1000.,$w_min*100.,$fl_nm_tmp));
	$fl_nm_arr=$fl_nm_arr." ".$fl_nm_tmp;
    } # endfor
    &cmd_prc(sprintf("ncecat -O -C -v IWC %s /data2/zender/cld/foo.nc",$fl_nm_arr));
    &cmd_prc(sprintf("ncrename -O -v IWC,IWC_srt /data2/zender/cld/foo.nc"));
    
    # Append the initial IWC to all the updraft files
    for ($w=$w_min;$w<=$w_max;$w+=$w_incr) {
	&cmd_prc(sprintf("ncks -A -C -v IWC_srt /data2/zender/cld/foo.nc /data2/zender/cld/cld_z_%02d_100_ss.nc",$w*100.));
    } # endfor
    &cmd_prc(sprintf("/bin/rm /data2/zender/cld/*foo*"));


} # $temperature_ss

sub cmd_prc{
    my($cmd)=@_;
    print STDOUT $cmd."\n";
    $status=open(CMD,$cmd." |");
    if(!$status){print STDOUT "Unable to open() on ".$cmd;}
    while(<CMD>){print STDOUT $_;}
    close CMD
    } # end cmd_prc()
