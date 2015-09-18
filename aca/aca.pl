#!/usr/bin/perl 

# $Id$

# Purpose: Process lots of ACA data

# Usage:
# Debugging:
# ${HOME}/aca/aca.pl -RS6K -dbg 2 -srt 5 -end 1 & # recompile under RS6K
# ${HOME}/aca/aca.pl -dbg 5 -RS6K -str 16 -srt 420 -end 1020 -ncr 5 > foo
# ${HOME}/aca/aca.pl -sng -bnd 1 -dbg 2 &
# ${HOME}/aca/aca.pl -sng -bnd 1 -dbg 2 -srt 420 -end 1020 -ncr 60 
# ${HOME}/aca/aca.pl -sng -bnd 1 -dbg 5 -srt 600 -end 670 -ncr 10 &

# Production:

# ER2
# ${HOME}/aca/aca.pl -clr -dbg 2 -day 06 -str 4 -srt 670 -end 680 -ncr 5 > foo.06.clr.aer 2>&1 &
# ${HOME}/aca/aca.pl -clr -cln -dbg 2 -day 06 -str 4 -srt 670 -end 680 -ncr 5 > foo.06.clr.cln 2>&1 &

# cloudy turbid collocated with small droplets
# ${HOME}/aca/aca.pl --rds_lqd=07 -dbg 2 -day 30 -str 4 -srt 660 -end 660 -ncr 5 > foo1 2>&1 &
# ${HOME}/aca/aca.pl --rds_lqd=07 -dbg 2 -day 30 -str 4 -srt 725 -end 725 -ncr 5 > foo2 2>&1 &
# ${HOME}/aca/aca.pl --rds_lqd=07 -dbg 2 -day 30 -str 4 -srt 755 -end 755 -ncr 5 > foo3 2>&1 &
# ${HOME}/aca/aca.pl --rds_lqd=07 -dbg 2 -day 30 -str 4 -srt 775 -end 775 -ncr 5 > foo4 2>&1 &
# ${HOME}/aca/aca.pl --rds_lqd=07 -dbg 2 -day 30 -str 4 -srt 780 -end 780 -ncr 5 > foo5 2>&1 &

# cloudy turbid
# ${HOME}/aca/aca.pl -dbg 2 -day 30 -str 4 -srt 420 -end 1020 -ncr 5 > foo.30.cld.aer 2>&1 &

# cloudy clean
# ${HOME}/aca/aca.pl -cln -dbg 2 -day 30 -str 4 -srt 420 -end 1020 -ncr 5 > foo.30.cld.cln 2>&1 &

# clear turbid
# ${HOME}/aca/aca.pl -clr -dbg 2 -day 11 -str 4 -srt 420 -end 1020 -ncr 15 > foo.11.clr.aer 2>&1 &
# ${HOME}/aca/aca.pl -clr -dbg 2 -day 13 -str 4 -srt 420 -end 1020 -ncr 30 > foo.13.clr.aer 2>&1 &
# ${HOME}/aca/aca.pl -clr -dbg 2 -day 15 -str 4 -srt 420 -end 1020 -ncr 30 > foo.15.clr.aer 2>&1 &
# ${HOME}/aca/aca.pl -clr -dbg 2 -day 17 -str 4 -srt 420 -end 1020 -ncr 30 > foo.17.clr.aer 2>&1 &
# ${HOME}/aca/aca.pl -clr -dbg 2 -day 19 -str 4 -srt 420 -end 1020 -ncr 30 > foo.19.clr.aer 2>&1 &
# ${HOME}/aca/aca.pl -clr -dbg 2 -day 30 -str 4 -srt 420 -end 1020 -ncr 30 > foo.30.clr.aer 2>&1 &

# clear clean
# ${HOME}/aca/aca.pl -clr -cln -dbg 2 -day 11 -str 4 -srt 420 -end 1020 -ncr 15 > foo.11.clr.cln 2>&1 &
# ${HOME}/aca/aca.pl -clr -cln -dbg 2 -day 13 -str 4 -srt 420 -end 1020 -ncr 30 > foo.13.clr.cln 2>&1 &
# ${HOME}/aca/aca.pl -clr -cln -dbg 2 -day 15 -str 4 -srt 420 -end 1020 -ncr 30 > foo.15.clr.cln 2>&1 &
# ${HOME}/aca/aca.pl -clr -cln -dbg 2 -day 17 -str 4 -srt 420 -end 1020 -ncr 30 > foo.17.clr.cln 2>&1 &
# ${HOME}/aca/aca.pl -clr -cln -dbg 2 -day 19 -str 4 -srt 420 -end 1020 -ncr 30 > foo.19.clr.cln 2>&1 &
# ${HOME}/aca/aca.pl -clr -cln -dbg 2 -day 30 -str 4 -srt 420 -end 1020 -ncr 30 > foo.30.clr.cln 2>&1 &
                           
# Specify modules
use Getopt::Long;
$Getopt::Long::ignorecase=0; # Turn on case-sensitivity
use POSIX qw(ceil floor);

require '/home/zender/perl/csz.pl'; # Personal library, &date_time(), &cmd_prc(), ...

# Set output flushing to help debugging on hard crashes. 
# These options update the filehandle after every output statement.
# See the Perl manual, p. 110.
select((select(STDOUT),$|=1)[0]);
select((select(STDERR),$|=1)[0]);

# Set defaults 
$False=0;
$True=1;

$recompile=$False;
$CRR_FL_IN_WAS_RCP=$False;
$PVM_ARCH=$ENV{'PVM_ARCH'};
$HOST=$ENV{'HOST'};
$aer=$True;
$cln=!$aer;
$cld=$True;
$clr=!$cld;
$ThD71=$False;
$day=11;
$dbg_bnd=1;
$dbg_lvl=0;
$end_idx=480; # min, or just a cardinal #
@fl_in_stbs=('arese_clm_tst.nc','mls_clm_cld.nc','H2O.nc','CO2.nc','O2.nc','O3.nc','O4.nc','OH.nc','NO2.nc','lqd_10.nc','ice_20.nc','aer_dust_like.nc','aer_mineral.nc','aer_sulphate','nst_FSBR.nc','nst_TSBR.nc');
$ncr_idx=5; # min, or just a cardinal #
$str_nbr=4;
$rds_lqd=10; # microns
$sng_bnd_cmp=$False;
$srt_idx=480; # min, or just a cardinal #

# Parse command line arguments 
$result=GetOptions(
		   'recompile!',
		   'aer!',
		   'cln!',
		   'cld!',
		   'clr!',
		   'sng!',
		   'bnd=i',
		   'dbg=i',
		   'day=i',
		   'end=i',
		   'ncr=i',
		   'rds_lqd=i',
		   'str=i',
		   'srt=i',
		   'ThD71!',
		   'plr=i',
		   );
				# 
if (defined($opt_recompile)) {$recompile=!$recompile};
if (defined($opt_aer)) {$aer=!$aer};
if (defined($opt_cln)) {$cln=!$cln};
if (defined($opt_clr)) {$clr=!$clr};
if (defined($opt_cld)) {$cld=!$cld};
if (defined($opt_sng)) {$sng_bnd_cmp=!$sng_bnd_cmp};
if (defined($opt_bnd)) {$dbg_bnd=$opt_bnd};
if (defined($opt_dbg)) {$dbg_lvl=$opt_dbg};
if (defined($opt_day)) {$day=$opt_day};
if (defined($opt_end)) {$end_idx=$opt_end};
if (defined($opt_ncr)) {$ncr_idx=$opt_ncr};
if (defined($opt_rds_lqd)) {$rds_lqd=$opt_rds_lqd};
if (defined($opt_str)) {$str_nbr=$opt_str};
if (defined($opt_srt)) {$srt_idx=$opt_srt};
if (defined($opt_ThD71)) {$ThD71=!$ThD71};
				 
# Ensure switches are in exclusive state after parsing user options
$aer=!$cln;
$cld=!$clr;

# Do some OS-dependent housekeeping to speed things up 
# and avoid possible ambiguities.
if($PVM_ARCH eq 'SUNMP' || $PVM_ARCH eq 'SUN4' || $PVM_ARCH eq 'SUN4SOL2' || $PVM_ARCH eq 'LINUX' ){
    $NEED_TO_RCP_FL_OUT=$False;
    @fl_in_dirs=('/data/zender/aca');
    $fl_prf_in_dir='/data/zender/arese/clm';
    $fl_out_dir='/data/zender/arese/mdl';
    $rcp_out_dir='/data/zender/arese/mdl';
}elsif($PVM_ARCH eq 'RS6K' || $PVM_ARCH eq 'SGI64'){		 

    $NEED_TO_RCP_FL_OUT=$True;
    @fl_in_dirs=('/usr/tmp/zender');
    $fl_out_dir='/usr/tmp/zender';
    $fl_prf_in_dir='/usr/tmp/zender';
    $home_mch='flagstaff.cgd.ucar.edu';
    $home_mch_data_dir='/data/zender/arese/mdl';
} # endif PVM_ARCH is RS6K

if($dbg_lvl >= 1){		 
    printf STDERR ("HOST is %s\n",$HOST);
    printf STDERR ("PVM_ARCH is %s\n",$PVM_ARCH);
    printf STDERR ("Debug level is %d\n",$dbg_lvl);
    printf STDERR ("Search paths for any pre-existing files:\n");
    for($idx=0;$idx<=$#fl_in_dirs;$idx++){
	printf STDERR ("%s\n",$fl_in_dirs[$idx]);
    } # end loop over idx
    printf STDERR ("Output path for storage of processed data: %s\n",$fl_out_dir);
    print STDERR "Processing ACA data...\n";
} # endif dbg
				 
if($PVM_ARCH eq 'RS6K' || $PVM_ARCH eq 'SGI64'){		 

# Remove data files that might have changed

    `/bin/rm -f mls.cld.nc`;
    print STDOUT ("Removed some old data files.\n");

    if($recompile){		 
        
# Copy the source from flagstaff and compile it

        `/bin/rm -f date_time.c`;
        `/bin/rm -f parameter.com`;
        `/bin/rm -f swnb.F`;
        `/bin/rm -f disort.f`;
        `/bin/rm -f linpak.f`;
        `/bin/rm -f machine.f`;
        `/bin/rm -f flx_slr_frc.f`;
        print STDOUT ("Removed any old source files.\n");
        
        `rcp flagstaff.cgd.ucar.edu:/home/zender/c/date_time.c .`;
        `rcp flagstaff.cgd.ucar.edu:/home/zender/f/parameter.com .`;
        `rcp flagstaff.cgd.ucar.edu:/home/zender/aca/swnb.F .`;
        `rcp flagstaff.cgd.ucar.edu:/home/zender/aca/disort.f .`;
        `rcp flagstaff.cgd.ucar.edu:/home/zender/aca/linpak.f .`;
        `rcp flagstaff.cgd.ucar.edu:/home/zender/aca/machine.f .`;
        `rcp flagstaff.cgd.ucar.edu:/home/zender/aca/flx_slr_frc.f .`;
        print STDOUT ("Copied over a new set of source files.\n");
        
        `gcc -c -O -DRS6K date_time.c -o date_time.o`;
        `/lib/cpp -P -DRS6K -I./ -I/usr/local/include swnb.F swnb.f`;
        `xlf -c -O -NS2000 -qfixed=132 -qmaxmem=16000 swnb.f -o swnb.o`;
        `xlf -c -O -NS2000 -qfixed=132 disort.f -o disort.o`;
        `xlf -c -O -NS2000 -qfixed=132 linpak.f -o linpak.o`;
        `xlf -c -O -NS2000 -qfixed=132 machine.f -o machine.o`;
        `xlf -c -O -NS2000 -qfixed=132 flx_slr_frc.f -o flx_slr_frc.o`;
        
        `xlf -o /crestone/u1/zender/bin/RS6K/swnb -O -NS2000 -qfixed=132 swnb.o date_time.o disort.o linpak.o machine.o flx_slr_frc.o -L/usr/local/lib -lm -lnetcdf -lncaru`;
        print STDOUT ("Compiled a new version of swnb.\n");
        
    } # endif recompile
} # endif PVM_ARCH is RS6K or SGI64

# Assemble all the input files
foreach $fl_in_stb (@fl_in_stbs) {
    
# Start with the presumption that the desired input file is not
# available locally, then search for it, then go to the appropriate
# source retrieval mechanism.
    $CRR_FL_IN_ALREADY_ON_DISK=$False;

    foreach $in_dir (@fl_in_dirs) {
	$foo=$in_dir.'/'.$fl_in_stb;
	if(-e $foo){
	    $fl_in=$foo;
	    $CRR_FL_IN_ALREADY_ON_DISK=$True;
	    $NEED_TO_RCP_FL_IN=$False;
	    last;
	} # end if the required file is already on disk
    } # end loop over input directories

    if(!$CRR_FL_IN_ALREADY_ON_DISK){
        if($PVM_ARCH eq 'RS6K' || $PVM_ARCH eq 'SGI64'){		 
            $NEED_TO_RCP_FL_IN=$True;
        } # endelse
    } # end if the required file is not already on disk
    
    if ($NEED_TO_RCP_FL_IN){
        $in_dir=$fl_out_dir;
        &cmd_prc('rcp flagstaff.cgd.ucar.edu:/data/zender/aca/'.$fl_in_stb.' '.$in_dir.'/'.$fl_in_stb);
        $fl_in=$in_dir.'/'.$fl_in_stb;
        $CRR_FL_IN_WAS_RCP=$True;
    } # endif
} # end loop over input files

$cnt=0;
for ($day_min=$srt_idx;$day_min<=$end_idx;$day_min+=$ncr_idx){
    $min=$day_min % 60;
    $hr=$day_min/60.;
    $HHMM=sprintf("%02d%02d",$hr,$min);

# NB: this works, but using reals in the loop may cause it to skip a time
#    $hr=$day_min/60.;
#    $min=&round(($hr-floor($hr))*60.);
#    if($min == 60){$min=0;$hr=int($hr+1);}
#    $HHMM=sprintf("%02d%02d",$hr,$min);

# NB: this does not work
#    $hr=$day_min/60.;
#    $min=($hr-floor($hr))*60.;
#    if(&round($min) == 60){$HHMM=sprintf("%02d00",int($hr+1));}else{$HHMM=sprintf("%02d%02d",$hr,$min)};

#    if($dbg_lvl > 3){printf STDOUT ("day_min = %f, round(day_min) = %f, hr = %f, min = %f, HHMM = %s\n",$day_min,&round($day_min),$hr,$min,$min,$HHMM);}

    $cnt=$cnt+1;
    $fl_out_stb=sprintf("9510%02d_%s_arese_mdl_",$day,$HHMM);
    if($cld){$fl_out_stb=$fl_out_stb.'cld_';}
    if($clr){$fl_out_stb=$fl_out_stb.'clr_';}
    if($aer){$fl_out_stb=$fl_out_stb.'aer';}
    if($cln){$fl_out_stb=$fl_out_stb.'cln';}
    if(defined($opt_rds_lqd)){
	$rds_lqd_stb=sprintf("_%02d",$rds_lqd);
	$fl_out_stb=$fl_out_stb.$rds_lqd_stb;
	$fl_lqd='/data/zender/aca/lqd'.$rds_lqd_stb.'.nc';
    };
    $fl_out_stb=$fl_out_stb.'.nc';
    $fl_out=$fl_out_dir.'/'.$fl_out_stb;
    $fl_prf_stb=sprintf("9510%02d_%s_arese_clm.nc",$day,$HHMM);
    $fl_prf=$fl_prf_in_dir.'/'.$fl_prf_stb;
    
    if($dbg_lvl == 1){
	printf STDOUT ("cnt = %d, idx = %f, fl_prf = %s fl_out = %s\n",
		       $cnt,
		       $idx,
		       $fl_prf,
		       $fl_out_stb);
    } # end if dbg 

    $swnb_cmd='swnb2';
# -J -K -Q -X -Y removes Planck source function, O2O2, H2OH2O, NO2, OH 
#    $swnb_cmd=$swnb_cmd.' -J -K -Q -X -Y';
# -Q -K -U removes H2OH2O but adds O2O2 and O2N2
    $swnb_cmd=$swnb_cmd.' -Q -K -U';
    if($cln){$swnb_cmd=$swnb_cmd.' -A -B';}
    if($clr){$swnb_cmd=$swnb_cmd.' -I -L';}
    $swnb_cmd=sprintf("%s -s %d",$swnb_cmd,$str_nbr);
    if($sng_bnd_cmp){$swnb_cmd=sprintf("%s -E -e %d",$swnb_cmd,$dbg_bnd);}
    if($ThD71){$swnb_cmd=$swnb_cmd.' -T /data/zender/aca/spc_ThD71.nc';}
    if($NeL84){$swnb_cmd=$swnb_cmd.' -T /data/zender/aca/spc_NeL84.nc';}
    if($Kur95_01wvn){$swnb_cmd=$swnb_cmd.' -T /data/zender/aca/spc_Kur95_01wvn.nc';}
    if($Kur95_20wvn){$swnb_cmd=$swnb_cmd.' -T /data/zender/aca/spc_Kur95_20wvn.nc';}
    if(defined($opt_rds_lqd)){$swnb_cmd=sprintf("%s -l %s",$swnb_cmd,$fl_lqd);}
    $swnb_cmd=sprintf("%s -p %s -d %s",$swnb_cmd,$fl_prf,$fl_out);

    if($dbg_lvl < 3){&cmd_prc($swnb_cmd);}else{printf STDOUT ("\n%s\n",$swnb_cmd)};

    if ($NEED_TO_RCP_FL_OUT){
	&cmd_prc('rcp '.$fl_out.' '.$home_mch.':'.$home_mch_data_dir.'/'.$fl_out_stb);
	&cmd_prc('/bin/rm -f '.$fl_out);
    } # endif

} # end outer loop over $day_min
