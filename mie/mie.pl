#!/usr/bin/perl -w
# $Id$ 

# Purpose: Generate, process aerosol optical properties for enhanced absorption studies
# Call mie program for specified size distributions and wavelength grids
# View output size distributions with ${HOME}/idl/mie.pro:aer_gph(),psd_bch()

# Usage:
# Set at least one of the six Booleans CCM_SW, CCM_LW, GSFC, RGL, SWNB, or TOMS 
# Calls to CCM_SW, CCM_LW, GSFC, SWNB, and TOMS can be identical in all other parameters
# Calls to RGL (regular grid) should be identical except for additional specification of wvl_rgl_mnm and wvl_rgl_mxm
# Calls to SWNB should set bnd_nbr to a small number since resolution is so high
# fxm: handling of size bin-independent parameters is hinkey
# Goal is to allow them to be passed as scalars or arrays
# fxm: Allow choice of alternate solar spectra, e.g., Kur95 for SWNB
# NB: Most size and wavelength arguments passed in microns rather than meters

# Usage:
# scp ~/mie/mie.pl ~/mie/mie.sh esmf.ess.uci.edu:mie
# scp ~/mie/mie.pl ~/mie/mie.sh goldhill.cgd.ucar.edu:/fs/cgd/home0/zender/mie

# Debugging
# ${HOME}/mie/mie.pl --dbg=5
# ${HOME}/mie/mie.pl --dbg=8 --foo={1,2,3,4}
# ${HOME}/mie/mie.pl --dbg=1 
# ${HOME}/mie/mie.pl --dbg=1 --xpt_dsc=Resonance_experiment --fl_lbl=mie01 --fl_nbr=10 --aer_typ=h2o_lqd --wvl_rgl_mnm=0.3 --wvl_rgl_mxm=2.0 --wvl_nbr_dfl=10 --bnd_nbr_dfl=1 --sz_nbr_dfl=1 --sz_mnm_dfl=0.1 --sz_mxm_dfl=30.0 --rds_swa_dfl=10.0 --gsd_anl_dfl=1.6 --drc_out=/data/zender/tmp > ${DATA}/mie/mie_RGL.txt 2>&1

# Production liquid cloud droplet optics: 
# Fine wavelength grid:
# ${HOME}/mie/mie.pl --dbg=1 --fl_lbl=${CASEID} --fl_nbr=480 --aer_typ=h2o_lqd --wvl_rgl_mnm=0.2 --wvl_rgl_mxm=5.0 --wvl_nbr_dfl=100000 --bnd_nbr_dfl=1 --sz_nbr_dfl=30 --sz_mnm_dfl=0.1 --sz_mxm_dfl=30.0 --rds_swa_dfl=10.0 --gsd_anl_dfl=1.6 > ${DATA}/mie/mie_RGL.txt 2>&1
# "Coarse" wavelength grid:
# ${HOME}/mie/mie.pl --dbg=1 --fl_lbl=${CASEID} --fl_nbr=480 --aer_typ=h2o_lqd --wvl_rgl_mnm=0.2 --wvl_rgl_mxm=5.0 --wvl_nbr_dfl=1 --bnd_nbr_dfl=100000 --sz_nbr_dfl=30 --sz_mnm_dfl=0.1 --sz_mxm_dfl=30.0 --rds_swa_dfl=10.0 --gsd_anl_dfl=1.6 > ${DATA}/mie/mie_RGL.txt 2>&1

BEGIN{
    unshift @INC,$ENV{'HOME'}.'/perl'; # Location of csz.pl and DBG.pm HaS98 p. 170
} # end BEGIN

my $CVS_Header='$Id$';

# Specify modules
use strict; # Protect all namespaces
use Getopt::Long; # GNU-style getopt
use File::Basename; # For parsing filenames

# 3rd party modules

# Personal modules
use DBG; # Debugging constants
require 'csz.pl'; # Contains date_time()

# Set output flushing to help debugging on hard crashes 
# These options update the filehandle after every output statement
select((select(STDOUT),$|=1)[0]); # Camel book, p. 110
select((select(STDERR),$|=1)[0]); # Camel book, p. 110

# Timing information
my ($lcl_date_time,$srt_usr_tm,$srt_sys_tm,$srt_child_usr_tm,$srt_child_sys_tm);
($lcl_date_time,$srt_usr_tm,$srt_sys_tm,$srt_child_usr_tm,$srt_child_sys_tm)=time_srt();
printf STDOUT ("Start user time %f\n",$srt_usr_tm);

# Declare local variables
my ($idx,$rcd);
my ($prg_nm,$prg_dsc,$prg_vrs,$prg_date);
my ($fl_out,$fl_out_nc,$fl_out_txt,$fl_sfx);
my ($fl_idx,$psd_idx,$wvl_fl_sng);
my ($cmd_sng,$drc_rgs_sng,$pcp_rgs_sng);
my ($rst_cmd); # [sng] Command to write restart file
my ($wvl_dlt_mcr,$wvl_rgs_sng,$sz_rgs_sng);
my (@sz_mcr); # [um] Size measurements to pass

# Set defaults 
my $False=0;
my $True=1;

my $CVS_Id='$Id$';
my $CVS_Revision='$Revision$';
my $CVS_Date='$Date$';
my $pi=atan2(1,1)*4.0;
my $nvr_data=$ENV{'DATA'} ? $ENV{'DATA'} : '';
my $nvr_data_out=$ENV{'DATA_OUT'} ? $ENV{'DATA_OUT'} : '';
my $nvr_data_rt=$ENV{'DATA_RT'} ? $ENV{'DATA_RT'} : '';
my $wvl_rgl_mnm_mcr=0.495; # [um] Wavelength minimum for regular grids
my $wvl_rgl_mxm_mcr=0.505; # [um] Wavelength maximum for regular grids
my $dmt_pcp_nma_cnv_mcr=1000.0; # [um] Diameter number median analytic, raindrop, convective
my $dmt_pcp_nma_str_mcr=400.0; # [um] Diameter number median analytic, raindrop, stratiform
my $flx_vlm_pcp_cnv=2.77e-6; # [m s-1]=[m3 m-2 s-1] Precipitation volume flux, convective [10 mm hr-1]
my $flx_vlm_pcp_str=0.277e-6; # [m s-1]=[m3 m-2 s-1] Precipitation volume flux, stratiform [1 mm hr-1]

# Set defaults for command line arguments
my $Boolean=$True;
my $CCM_LW=$False;
my $CCM_SW=$False;
my $CNV=$False; # [flg] Convective precipitation characteristics
my $GSFC=$False;
my $HELP=$False;
my $RGL=$False; # [flg] Regular wavelength grid
my $RST=$False; # [flg] Current shell is child of restart command
my $STR=$False; # [flg] Stratiform precipitation characteristics
my $SWNB=$False;
my $TOMS=$False;
my $aer_typ='h2o_lqd'; # [sng] Aerosol type
my $drc_dat="$nvr_data_rt"; # [sng] Data directory
my $drc_out="$nvr_data_out"; # [sng] Output directory
my $fl_rst='restart'; # [sng] Restart file
my $sz_typ='rds_swa'; # [sng] Size measurement to pass: "rds_nma", "rds_swa", "dmt_vma"
my $bch_flg=$True; # [flg] Batch behavior
my $bnd_SW_LW=5.0e-6; # [m] Boundary between SW and LW weighting
my $bnd_nbr_dfl=1; # [frc] Number of sub-bands per output band
my $wvl_nbr_dfl=1; # [frc] Number of output wavelength bands
my $dbg_lvl=0;
my $dsd_mnm_mcr=1.0; # [um] Minimum diameter in raindrop distribution
my $dsd_mxm_mcr=4000.0; # [um] Maximum diameter in raindrop distribution
my $dsd_nbr=200; # [nbr] Number of raindrop size bins
my $fdg_val=1.0; # [frc] Tuning factor for all bands
my $flt_foo=73.37;
my $foo=73;
my $ftn_fxd_flg=$False; # [flg] Fortran fixed format
my $gsd_anl_dfl=2.0; # [frc] Geometric standard deviation analytic
my $gsd_pcp_anl=1.86; # [frc] Geometric standard deviation, raindrop NGD94 p. 2337 Table 1 
my $hrz_flg=$False; # [flg] Print size-resolved optical properties at debug wavelength
my $idx_rfr_prt_dfl=''; # [sng] Index of refraction of particle
my $mss_val=1.0e36; # [frc] Missing value
my $psd_nbr=1; # [nbr] Number of particle size distributions
my $fl_lbl=''; # [sng] Optional file label (e.g., ${CASEID})
my $xpt_dsc=''; # [sng] Experiment description
my $fl_nbr=1; # [nbr] Number of files
my $psd_typ='lognormal'; # [sng] Particle size distribution type
my $rds_nma_dfl=0.2986; # [um] Number median radius analytic
my $dmt_vma_dfl=3.5; # [um] Volume median diameter analytic
my $rds_swa_dfl=10.0; # [um] Surface area weighted radius analytic
my $spc_idx_dfl=''; # [sng] Species index for Fortran data
my $sz_grd='log';
my $sz_mnm_dfl=0.01; # [um] Minimum radius in aerosol distribution
my $sz_mxm_dfl=1.0; # [um] Maximum radius in aerosol distribution
my $sz_nbr_dfl=2; # [nbr] Number of aerosol size bins
my $thr_nbr=0; # [nbr] Thread number
my $wvl_grd_sng='RGL'; # [sng] Type of wavelength grid
my $wrn_ntp_flg=$False; # [flg] Print WARNINGs from ntp_vec()

my @bnd_nbr_arr=(); # [frc] Number of sub-bands per output band
my @spc_idx_arr=(); # [sng] Species index for Fortran data
my @foo=();
my @gsd_anl_arr=(); # [frc] Geometric standard deviation analytic
my @idx_rfr_prt_arr=(); # [sng] Index of refraction of particle
my @rds_nma_arr=(); # [um] Number median radius analytic
my @dmt_vma_arr=(); # [um] Volume median diameter analytic
my @rds_swa_arr=(); # [um] Surface area weighted radius analytic
my @sz_mnm_arr=(); # [um] Minimum radius in aerosol distribution
my @sz_mxm_arr=(); # [um] Maximum radius in aerosol distribution
my @sz_nbr_arr=(); # [nbr] Number of aerosol size bins
my @wvl_nbr_arr=(); # [frc] Number of output wavelength bands
my @wvl_mnm_mcr=(); # [um] Minimum wavelength
my @wvl_mxm_mcr=(); # [um] Maximum wavelength

# Derived fields
$prg_dsc='Mie Resonance Processor'; # Program description
if(length($CVS_Id) > 4){($prg_nm,$prg_vrs)=$CVS_Id =~ /: (.+).pl,v ([\d.]+)/;}else{$prg_vrs='Unknown';}
$prg_vrs.='*' if length('$Locker:  $ ') > 12; # Tack '*' if it is not checked in into CVS
$prg_nm=$0; # $0 is program name Camel p. 136
if(length($CVS_Date) > 6){($prg_date)=unpack '@7 a19',$CVS_Date;}else{$prg_date='Unknown';}

my ($arg_idx,$cmd_ln);
$cmd_ln=$0.' ';
for($arg_idx=0;$arg_idx<=$#ARGV;$arg_idx++){
    $cmd_ln.=$ARGV[$arg_idx].' ';
} # end loop over arg
if($dbg_lvl != $dbg_off){printf ("$prg_nm: \$cmd_ln = $cmd_ln\n");} # endif dbg

# Parse command line arguments: '!' means Boolean, '|' is OR, '=' specifies required argument: 'i' is integer, 'f' is float, 's' is string
$rcd=Getopt::Long::Configure('no_ignore_case'); # Turn on case-sensitivity
$rcd=GetOptions( # man Getopt::Long
		'Boolean!' => \$Boolean,
		'CNV!' => \$CNV,
		'RST!' => \$RST, # [flg] Current shell is child of restart command
		'STR!' => \$STR,
		'aer_typ=s' => \$aer_typ, # [sng] Aerosol type
		'bch_flg!' => \$bch_flg, # [flg] Batch behavior
		'bnd_nbr_arr=i' => \@bnd_nbr_arr, # [frc] Number of sub-bands per output band
		'bnd_nbr_dfl=i' => \$bnd_nbr_dfl, # [frc] Number of sub-bands per output band
		'dbg_lvl=i' => \$dbg_lvl,
		'dmt_vma_arr=f' => \@dmt_vma_arr, # [um] Volume median diameter analytic
		'dmt_vma_dfl=f' => \$dmt_vma_dfl, # [um] Volume median diameter analytic
		'drc_dat=s' => \$drc_dat, # [sng] Data directory
		'drc_out=s' => \$drc_out, # [sng] Output directory
		'fdg_val=f' => \$fdg_val, 
		'fl_lbl=s' => \$fl_lbl, # [sng] Optional file label (e.g., ${CASEID})
		'xpt_dsc=s' => \$xpt_dsc, # [sng] Experiment description
		'fl_nbr=i' => \$fl_nbr, # [nbr] Number of files
		'flt_foo=f' => \$flt_foo, 
		'foo=i' => \@foo,
		'ftn_fxd!' => \$ftn_fxd_flg, # [flg] Fortran fixed format
		'gsd_anl_arr=f' => \@gsd_anl_arr, # [frc] Geometric standard deviation analytic
		'gsd_anl_dfl=f' => \$gsd_anl_dfl, # [frc] Geometric standard deviation analytic
		'idx_rfr_prt_arr=s' => \@idx_rfr_prt_arr, # [sng] Index of refraction of particle
		'idx_rfr_prt_dfl=s' => \$idx_rfr_prt_dfl, # [sng] Index of refraction of particle
		'gsd_pcp_anl=f' => \$gsd_pcp_anl, # [frc] Geometric standard deviation, raindrop
		'help|usage|?', => \$HELP,
		'hrz_flg!' => \$hrz_flg, # [flg] Print size-resolved optical properties at debug wavelength
		'psd_nbr=i' => \$psd_nbr, # [nbr] Number of particle size distributions
		'rds_nma_arr=f' => \@rds_nma_arr, # [um] Number median radius analytic
		'rds_nma_dfl=f' => \$rds_nma_dfl, # [um] Number median radius analytic
		'rds_swa_arr=f' => \@rds_swa_arr, # [um] Surface area weighted radius analytic
		'rds_swa_dfl=f' => \$rds_swa_dfl, # [um] Surface area weighted radius analytic
		'spc_idx_arr=s' => \@spc_idx_arr, # [sng] Species index for Fortran data
		'spc_idx_dfl=s' => \$spc_idx_dfl, # [sng] Species index for Fortran data
		'sz_mnm_arr=f' => \@sz_mnm_arr, # [um] Minimum radius in aerosol distribution
		'sz_mnm_dfl=f' => \$sz_mnm_dfl, # [um] Minimum radius in aerosol distribution
		'sz_mxm_arr=f' => \@sz_mxm_arr, # [um] Maximum radius in aerosol distribution
		'sz_mxm_dfl=f' => \$sz_mxm_dfl, # [um] Maximum radius in aerosol distribution
		'sz_nbr_arr=i' => \@sz_nbr_arr, # [nbr] Number of aerosol size bins
		'sz_nbr_dfl=i' => \$sz_nbr_dfl, # [nbr] Number of aerosol size bins
		'sz_typ=s' => \$sz_typ, # [sng] Size measurement to pass: "rds_nma", "rds_swa", "dmt_vma"
		'thr_nbr=f' => \$thr_nbr, # [nbr] Thread number
		'wrn_ntp_flg!' => \$wrn_ntp_flg, # [flg] Print WARNINGs from ntp_vec()
		'wvl_grd_sng=s' => \$wvl_grd_sng, # [sng] Type of wavelength grid
		'wvl_nbr_arr=i' => \@wvl_nbr_arr, # [nbr] Number of output wavelength bands
		'wvl_nbr_dfl=i' => \$wvl_nbr_dfl, # [nbr] Number of output wavelength bands
		'wvl_rgl_mnm_mcr=f' => \$wvl_rgl_mnm_mcr, # [um] Wavelength minimum for regular grids
		'wvl_rgl_mxm_mcr=f' => \$wvl_rgl_mxm_mcr, # [um] Wavelength maximum for regular grids
		 ); # end GetOptions arguments

# Automatically determine if restart run should be performed
my $RST_READ=$False; # [flg] Read restart files
my $RST_WRT=$True; # [flg] Write restart files
if($RST){printf "Current batch script is child spawned by restart command\n";}else{printf "Current batch script is parent, will search for valid restart file...\n";}
# Avoid infinite recursion---read restart files only in parent shell mode
if(!$RST){
    $RST_READ=$True if -x $fl_rst; # if file is executable by effective uid/gid
    if($RST_READ){
	printf "Restart file $fl_rst is present and executable...will spawn restart child\n";
# Perform manual restart rather than initial run logic 
	open(FL_RST,$fl_rst) || die "cannot open $fl_rst";
	my ($cmd_old,$rst_idx,$rst_nbr);
	$rst_idx=0; # 
	while(<FL_RST>){
	    $rst_cmd=$_;
	    print STDOUT "Old restart command = $rst_cmd\n";
	} # end while <FL_RST>
	close FL_RST;
	$rst_cmd.=' --RST'; # [sng] Command to invoke restart child
# Spawn restart child
	print "Spawning restart child with command = $rst_cmd\n";
	if($dbg_lvl < $dbg_io){&cmd_prc($rst_cmd);}else{print $rst_cmd."\n";}
    } # endif RST_READ
} # endif RST

# Here follows the real program...
# Implement this only when present shell is not a restart child
if(!$RST_READ){
    if(!$RST){print "No valid restart file found, continuing as initial invocation...\n";}
    
    if(!$STR && !$CNV){$CNV=$True;}elsif($STR && $CNV){die "$prg_nm: ERROR Both CNV and STR are set\n";}
    if($HELP){&usg_prn();exit 0;} # end HELP
    
# Lists not specified on command line get default values here
# Initialize arrays with repetition operator, Camel p. 82
    if(!@bnd_nbr_arr){@bnd_nbr_arr=($bnd_nbr_dfl)x$fl_nbr;} # [frc] Number of sub-bands per output band
    if(!@gsd_anl_arr){@gsd_anl_arr=($gsd_anl_dfl)x$fl_nbr;} # [frc] Geometric standard deviation analytic
    if(!@idx_rfr_prt_arr){@idx_rfr_prt_arr=($idx_rfr_prt_dfl)x$fl_nbr;} # [sng] Index of refraction of particle
    if(!@spc_idx_arr){@spc_idx_arr=($spc_idx_dfl)x$fl_nbr;} # [sng] Species index for Fortran data
    if(!@sz_mnm_arr){@sz_mnm_arr=($sz_mnm_dfl)x$fl_nbr;} # [um] Minimum radius in aerosol distribution
    if(!@sz_mxm_arr){@sz_mxm_arr=($sz_mxm_dfl)x$fl_nbr;} # [um] Maximum radius in aerosol distribution
    if(!@sz_nbr_arr){@sz_nbr_arr=($sz_nbr_dfl)x$fl_nbr;} # [nbr] Number of aerosol size bins
    if(!@wvl_nbr_arr){@wvl_nbr_arr=($wvl_nbr_dfl)x$fl_nbr;} # [frc] Number of output wavelength bands
    if(!@foo){@foo=($foo)x$fl_nbr;}
    
# Expand default into full array if necessary
    if(!@dmt_vma_arr){
	for ($fl_idx=0;$fl_idx<$fl_nbr;$fl_idx++){
	    $dmt_vma_arr[$fl_idx]=$dmt_vma_dfl; # [um] Volume median diameter analytic
	} # end loop over $fl
    } # endif !defined(@dmt_vma_arr)
    if(!@rds_nma_arr){
	for ($fl_idx=0;$fl_idx<$fl_nbr;$fl_idx++){
	    $rds_nma_arr[$fl_idx]=$rds_nma_dfl; # [um] Number median radius analytic
	} # end loop over $fl
    } # endif !defined(@rds_nma_arr)
    if(!@rds_swa_arr){
	for ($fl_idx=0;$fl_idx<$fl_nbr;$fl_idx++){
	    $rds_swa_arr[$fl_idx]=$rds_swa_dfl; # [um] Surface area weighted radius analytic
	} # end loop over $fl
    } # endif !defined(@rds_swa_arr)
    if(!@wvl_mnm_mcr || !@wvl_mxm_mcr){
	$wvl_dlt_mcr=($wvl_rgl_mxm_mcr-$wvl_rgl_mnm_mcr)/$fl_nbr; # [um] Bandwidth
	for ($fl_idx=0;$fl_idx<$fl_nbr;$fl_idx++){
	    $wvl_mnm_mcr[$fl_idx]=$wvl_rgl_mnm_mcr+$fl_idx*$wvl_dlt_mcr; # [um] Minimum wavelength
	    $wvl_mxm_mcr[$fl_idx]=$wvl_rgl_mnm_mcr+($fl_idx+1)*$wvl_dlt_mcr; # [um] Maximum wavelength
	} # end loop over $fl
	$wvl_mxm_mcr[$fl_nbr-1]=$wvl_rgl_mxm_mcr; # [um] Maximum wavelength
    } # endif !defined(@wvl_mnm_mcr)
    
    if($wvl_grd_sng eq 'CCM_SW'){$CCM_SW=$True;$bnd_SW_LW=5.0e-6;}
    elsif($wvl_grd_sng eq 'CCM_LW'){$CCM_LW=$True;$bnd_SW_LW=3.0e-6;$fdg_val=1.22;}
    elsif($wvl_grd_sng eq 'GSFC'){$GSFC=$True;$bnd_SW_LW=3.0e-6;}
    elsif($wvl_grd_sng eq 'RGL'){$RGL=$True;$bnd_SW_LW=5.0e-6;}
    elsif($wvl_grd_sng eq 'SWNB'){$SWNB=$True;$bnd_SW_LW=5.0e-6;}
    elsif($wvl_grd_sng eq 'TOMS'){$TOMS=$True;$bnd_SW_LW=5.0e-6;}
    else{die "$prg_nm: ERROR Unknown wavelength grid string \$wvl_grd_sng = $wvl_grd_sng\n";}
    
# [sng] Size measurements that are passed
    if($sz_typ =~ 'rds_nma'){@sz_mcr=@rds_nma_arr;} 
    elsif($sz_typ =~ 'rds_swa'){@sz_mcr=@rds_swa_arr;}
    elsif($sz_typ =~ 'dmt_vma'){@sz_mcr=@dmt_vma_arr;}
    else{die "$prg_nm: ERROR Incorrect sz_typ string\n";}
    
# Print initialization state
    if($dbg_lvl > 1){print ("$prg_nm: $prg_dsc, version $prg_vrs of $prg_date\n");} # endif dbg
    if($dbg_lvl >= $dbg_fl){
	printf ("Initialization State:\n");
	printf ("cmd_ln = $cmd_ln\n");
	printf ("aer_typ = $aer_typ\n");
	printf ("wvl_grd=$wvl_grd_sng");
	if($wvl_grd_sng =~ 'rgl'){printf (", wvl_rgl_mnm = $wvl_rgl_mnm_mcr um, wvl_rgl_mxm = $wvl_rgl_mxm_mcr um");}
	printf ("\n");
	if($thr_nbr != 0){printf ("thr_nbr = $thr_nbr, ");}
	printf ("fl_nbr = $fl_nbr, wvl_nbr = $wvl_nbr_dfl, bnd_nbr = $bnd_nbr_dfl, sz_nbr = $sz_nbr_dfl\n");
	printf ("Mie calls per file = wvl_nbr*bnd_nbr*sz_nbr = %d\n",$wvl_nbr_dfl*$bnd_nbr_dfl*$sz_nbr_dfl);
	printf ("sz_typ = $sz_typ\n");
	printf ("spc_idx = $spc_idx_dfl\n");
	printf ("idx_rfr_prt = $idx_rfr_prt_dfl\n");
	printf ("bnd_SW_LW = $bnd_SW_LW um\n");
	printf ("drc_dat = $drc_dat, drc_out = $drc_out\n");
	printf ("fl_idx\twvl_mnm\twvl_mxm\tsz_mnm\tsz_mxm\t${sz_typ}\tgsd_anl\twvl_nbr\tbnd_nbr\tsz_nbr\n");
	printf ("#\t\tum\tum\tum\tum\tum\t\t\t\t\n");
	for ($fl_idx=0;$fl_idx<$fl_nbr;$fl_idx++){
	    print ("$fl_idx\t$wvl_mnm_mcr[$fl_idx]\t$wvl_mxm_mcr[$fl_idx]\t$sz_mnm_arr[$fl_idx]\t$sz_mxm_arr[$fl_idx]\t$sz_mcr[$fl_idx]\t$gsd_anl_arr[$fl_idx]\t$wvl_nbr_arr[$fl_idx]\t$bnd_nbr_arr[$fl_idx]\t$sz_nbr_arr[$fl_idx]\n");
	} # end loop over $aer
	print ("\n");
    } # endif dbg
    
    if($fdg_val != 1.0){print ("$prg_nm: WARNING Using non-unity \$fdg_val = $fdg_val. This value will multiply all computed extinction coefficients. This is intended to tune optical properties on coarse grids (e.g., CCM LW) to obtain reasonable agreement with high resolution grids. If current optical properties are not for CCM LW in CCM/CRM or this \$fdg_val was not tuned by you personally then be careful using these results!\n");}
    
    if($dbg_lvl == $dbg_crr){
	print ('@foo = ');
	foreach $foo (@foo){
	    print ("$foo, ");
	} # end loop over $foo
	print ("\n");
    } # endif dbg
    
    for ($fl_idx=0;$fl_idx<$fl_nbr;$fl_idx++){
	
# Construct loop-independent portion of command line arguments
	$cmd_sng="mie"; # [sng] Initialize to program name
	$cmd_sng.=" --dbg=$dbg_lvl --tst=nsz --aer_typ=$aer_typ"; # [sng] 
	if($thr_nbr != 0){$cmd_sng.=" --thr_nbr=$thr_nbr";}
	if($xpt_dsc){$cmd_sng.=" --xpt_dsc=$xpt_dsc";}
	if($bch_flg){$cmd_sng.=" --bch";} # [sng] 
	if($ftn_fxd_flg){$cmd_sng.=" --ftn_fxd";} # [sng] 
	if($hrz_flg){$cmd_sng.=" --hrz";} # [sng] 
	if(!$wrn_ntp_flg){$cmd_sng.=" --no_wrn_ntp";} # [sng] 
	if($fdg_val != 1.0){$cmd_sng.=" --fdg_val=$fdg_val";}
	$drc_rgs_sng=''; # [sng] Directory arguments
	if($drc_dat ne ''){$drc_rgs_sng.=" --drc_dat=$drc_dat";} # [sng] Directory arguments
	if($drc_out ne ''){$drc_rgs_sng.=" --drc_out=$drc_out";} # [sng] Directory arguments
	$wvl_rgs_sng=" --bnd_SW_LW=$bnd_SW_LW"; # [sng] Wavelength distribution arguments
	$sz_rgs_sng=" --psd_typ=$psd_typ --sz_grd=$sz_grd"; # [sng] Size distribution arguments
	$pcp_rgs_sng=' --dsd_mnm='.$dsd_mnm_mcr.' --dsd_mxm='.$dsd_mxm_mcr.' --dsd_nbr='.$dsd_nbr.' --gsd_pcp_anl='.$gsd_pcp_anl; # [sng] Precipitation distribution arguments
	if($CNV){$pcp_rgs_sng.=' --dmt_pcp_nma='.$dmt_pcp_nma_cnv_mcr.' --flx_vlm_pcp='.$flx_vlm_pcp_cnv;}else{$pcp_rgs_sng.=' --dmt_pcp_nma='.$dmt_pcp_nma_str_mcr.' --flx_vlm_pcp='.$flx_vlm_pcp_str;} # [sng] Precipitation distribution arguments
	
# Construct loop-dependent portion of command line arguments
	$fl_out='aer_'.$aer_typ.'_';
	$fl_out.=$sz_typ.'_'.sprintf('%02d',$sz_mcr[$fl_idx]).'_';
	if($fl_lbl){$fl_out.=$fl_lbl.'_';}
	if($spc_idx_arr[$fl_idx]){
	    $drc_rgs_sng.=" --spc_idx=$spc_idx_arr[$fl_idx]"; # [sng] Directory arguments
	    $fl_out.=$spc_idx_arr[$fl_idx].'_';
	} # endif spc_idx_arr
	if($RGL){
# Compute optical properties over specified wavelength grid
	    $wvl_fl_sng=sprintf('%010.7f_%010.7f',$wvl_mnm_mcr[$fl_idx],$wvl_mxm_mcr[$fl_idx]);
	    $fl_out.=$wvl_fl_sng;
	} # endif RGL
# Compute broadband optical properties for solar or terrestrial spectral regimes
	if($CCM_SW || $CCM_LW || $GSFC || $SWNB || $TOMS){$fl_out.=$wvl_grd_sng;}
	$fl_out_nc=$fl_out.'.nc';
	$fl_out_txt=$fl_out.'.txt';
	
	$drc_rgs_sng.=" $fl_out_nc"; # [sng] Directory arguments
	$sz_rgs_sng.=" --sz_nbr=$sz_nbr_arr[$fl_idx] --sz_mnm=$sz_mnm_arr[$fl_idx] --sz_mxm=$sz_mxm_arr[$fl_idx] --${sz_typ}=$sz_mcr[$fl_idx] --gsd_anl=$gsd_anl_arr[$fl_idx]"; # [sng] Size distribution arguments
	if($RGL){$wvl_rgs_sng.=" --wvl_grd=$wvl_grd_sng --wvl_nbr=$wvl_nbr_arr[$fl_idx] --bnd_nbr=$bnd_nbr_arr[$fl_idx] --wvl_mnm=$wvl_mnm_mcr[$fl_idx] --wvl_mxm=$wvl_mxm_mcr[$fl_idx]";} # [sng] Wavelength arguments
	if($CCM_SW || $CCM_LW || $GSFC || $SWNB || $TOMS){$wvl_rgs_sng.=" --wvl_grd=$wvl_grd_sng --bnd_nbr=$bnd_nbr_arr[$fl_idx]";} # [sng] Wavelength arguments
	if($idx_rfr_prt_arr[$fl_idx]){
	    $wvl_rgs_sng.=" --idx_rfr_prt=$idx_rfr_prt_arr[$fl_idx]"; # [sng] Wavelength arguments	    
	} # endif idx_rfr_prt_arr

	# Assemble final execution command from component substrings
	$cmd_sng.=$wvl_rgs_sng.$sz_rgs_sng.$pcp_rgs_sng.$drc_rgs_sng;
	
# Construct/update restart command for mie.pl for current 
	if($RST_WRT){
	    my $RST_rx=' --RST'; # [sng] Regular expression for fl_nbr of batch command
	    my $fl_nbr_rx='--fl_nbr=\d{1,6}'; # [sng] Regular expression for fl_nbr of batch command
	    my $fl_nbr_rpl=''; # [sng] Replacement for fl_nbr portion of batch command
	    my $fl_rmn=$fl_nbr-$fl_idx; # [nbr] Number of files remaining
	    my $wvl_mnm_rx='--wvl_rgl_mnm=-?([0-9]+(\.[0-9]*)?|\.[0-9]+)'; # [sng] Regular expression for wvl_mnm of batch command
	    my $wvl_mnm_rpl=''; # [sng] Replacement for wvl_mnm portion of batch command
	    my $xpt_dsc_rx='--xpt_dsc=(.*) --fl_lbl'; # [sng] Regular expression for experiment description
	    my $xpt_dsc_rpl=''; # [sng] Replacement for experiment description

	    $rst_cmd=$cmd_ln; # [sng] Batch command to restart processing
	    # Strip child invocation from batch command
	    if($rst_cmd =~ m/($RST_rx)/){
		if($dbg_lvl > 1){print STDOUT "Pattern is $RST_rx, first match is $1\n";} 
		$rst_cmd =~ s/$RST_rx//g;
	    }
	    # Update starting wavelength
	    $wvl_mnm_rpl="--wvl_rgl_mnm=$wvl_mnm_mcr[$fl_idx]"; # [sng] wvl_mnm portion of batch command
	    if($rst_cmd =~ m/($wvl_mnm_rx)/){
		if($dbg_lvl > 1){print STDOUT "Pattern is $wvl_mnm_rx, first match is $1\n";}
		$rst_cmd =~ s/$wvl_mnm_rx/$wvl_mnm_rpl/g;
	    }  # endif
	    # Update number of files
	    $fl_nbr_rpl="--fl_nbr=$fl_rmn"; # [sng] fl_nbr portion of batch command
	    if($rst_cmd =~ m/($fl_nbr_rx)/){
		if($dbg_lvl > 1){print STDOUT "Pattern is $fl_nbr_rx, first match is $1\n";}
		$rst_cmd =~ s/$fl_nbr_rx/$fl_nbr_rpl/g;
	    }  # endif
	    # Update experiment description
	    if($rst_cmd =~ m/($xpt_dsc_rx)/){
		if($dbg_lvl > 1){print STDOUT "Pattern is $xpt_dsc_rx, first match is $1\n";}
		$xpt_dsc_rpl='\''.$1.'\'';
#		$rst_cmd =~ s/$xpt_dsc_rx/$xpt_dsc_rpl/g;
	    } # endif
	    print STDOUT "Updated batch restart command: $rst_cmd\n";
	    open FL_RST,'>'.$fl_rst;
	    print FL_RST "$rst_cmd";
# Change mode to executable
	    chmod 0755, $fl_rst;
	    printf "Wrote restart information to file $fl_rst...\n";
	} # endif if $RST_WRT
	
	# Execute computation
	if($dbg_lvl < $dbg_io){&cmd_prc($cmd_sng);}else{print $cmd_sng."\n";}
	
    } # end loop over fl
    
} # endif $RST_READ

# Get elapsed times
time_end($srt_usr_tm,$srt_sys_tm,$srt_child_usr_tm,$srt_child_sys_tm);

