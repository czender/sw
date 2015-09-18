#!/usr/bin/perl -w
# $Id$ 

# Purpose: Call mie program for specified size distributions and wavelength grids
# This script generates all aerosol optical properties used in dust model studies
# View output size distributions with ${HOME}/idl/mie.pro:psd_bch()
# View output optical properties with ${HOME}/anl/spc_ln.ncl

# Usage:
# Set one the eight Booleans CAM_SW, CAM_LW, GSFC, RGL, RRTM_SW, RRTM_LW, SNICAR, SWNB, or TOMS
# Calls to SNICAR are shorthand for specific RGL configurations
# Calls to CAM_SW, CAM_LW, GSFC, RRTM_SW, RRTM_LW, SWNB, and TOMS can be identical in all other parameters
# Calls to RGL should be identical except for additional specification of wvl_mnm, wvl_mxm, and wvl_nbr
# Calls to SWNB should set bnd_nbr to a small number since resolution is so high
# fxm: handling of size bin-independent parameters is hinkey
# Goal is to allow them to be passed as scalars or arrays
# fxm: Allow choice of alternate solar spectra, e.g., Kur95 for SWNB
# NB: Most size and wavelength arguments passed in microns rather than meters

# Debugging
# chmod a+x ~/mie/psd.pl
# ${HOME}/mie/psd.pl --dbg=5
# ${HOME}/mie/psd.pl --dbg=8 --foo={1,2,3,4}
# ${HOME}/mie/psd.pl --dbg=1 
# Production:
# rsync tephra.ess.uci.edu:/data/zender/dst/mie /data/zender/dst

# Dust version 3 optics: Default as of 20030717
# ${HOME}/mie/psd.pl --dbg=1 --CAM_SW --ftn_fxd --psd_nbr=4 --spc_idx={01,02,03,04} --sz_mnm={0.05,0.5,1.25,2.5} --sz_mxm={0.5,1.25,2.5,5.0} --sz_nbr={200,25,25,25} --dmt_vma_dfl=3.5 > ${DATA}/dst/mie/psd_CAM_SW.txt 2>&1 &
# ${HOME}/mie/psd.pl --dbg=1 --CAM_LW --ftn_fxd --psd_nbr=4 --spc_idx={01,02,03,04} --sz_mnm={0.05,0.5,1.25,2.5} --sz_mxm={0.5,1.25,2.5,5.0} --sz_nbr={200,25,25,25} --dmt_vma_dfl=3.5 > ${DATA}/dst/mie/psd_CAM_LW.txt 2>&1 &
# ${HOME}/mie/psd.pl --dbg=1 --RGL --ftn_fxd --psd_nbr=4 --spc_idx={01,02,03,04} --sz_mnm={0.05,0.5,1.25,2.5} --sz_mxm={0.5,1.25,2.5,5.0} --sz_nbr={200,25,25,25} --dmt_vma_dfl=3.5 --wvl_mnm=0.625 --wvl_mxm=0.635 > ${DATA}/dst/mie/psd_ODX.txt 2>&1 &

# Dust version 4 optics: Default as of 20060502
# ${HOME}/mie/psd.pl --dbg=1 --CAM_SW --EMA=20060904 --vts_flg --dns_prt=2500.0 --psd_nbr=4 --spc_idx={01,02,03,04} --sz_mnm={0.05,0.5,1.25,2.5} --sz_mxm={0.5,1.25,2.5,5.0} --sz_nbr={200,25,25,25} --dmt_vma_dfl=3.5 > ${DATA}/dst/mie/psd_CAM_SW.txt 2>&1 &
# ${HOME}/mie/psd.pl --dbg=1 --CAM_LW --EMA=20060904 --vts_flg --dns_prt=2500.0 --psd_nbr=4 --spc_idx={01,02,03,04} --sz_mnm={0.05,0.5,1.25,2.5} --sz_mxm={0.5,1.25,2.5,5.0} --sz_nbr={200,25,25,25} --dmt_vma_dfl=3.5 > ${DATA}/dst/mie/psd_CAM_LW.txt 2>&1 &
# ${HOME}/mie/psd.pl --dbg=1 --SNICAR --EMA=20060904 --vts_flg --dns_prt=2500.0 --psd_nbr=4 --spc_idx={01,02,03,04} --sz_mnm={0.05,0.5,1.25,2.5} --sz_mxm={0.5,1.25,2.5,5.0} --sz_nbr={200,25,25,25} --dmt_vma_dfl=3.5 > ${DATA}/dst/mie/psd_SNICAR.txt 2>&1 &
# ${HOME}/mie/psd.pl --dbg=1 --SNICAR --EMA=1 --cmp_mtx='SiO2_avg_hitran96' --cmp_ncl='Fe2O3_doccd' --mss_frc_ncl=0.01 --psd_nbr=4 --spc_idx={01,02,03,04} --sz_mnm={0.05,0.5,1.25,2.5} --sz_mxm={0.5,1.25,2.5,5.0} --sz_nbr={200,25,25,25} --dmt_vma_dfl=3.5 > ${DATA}/dst/mie/psd_SNICAR.txt 2>&1 &
# ${HOME}/mie/psd.pl --dbg=1 --CAM_SW --cmp_prt='Fe2O3_doccd' --psd_nbr=4 --spc_idx={01,02,03,04} --sz_mnm={0.05,0.5,1.25,2.5} --sz_mxm={0.5,1.25,2.5,5.0} --sz_nbr={200,25,25,25} --dmt_vma_dfl=3.5 > ${DATA}/dst/mie/psd_CAM_SW.txt 2>&1 &
# ${HOME}/mie/psd.pl --dbg=1 --CAM_LW --psd_nbr=4 --spc_idx={01,02,03,04} --sz_mnm={0.05,0.5,1.25,2.5} --sz_mxm={0.5,1.25,2.5,5.0} --sz_nbr={200,25,25,25} --dmt_vma_dfl=3.5 > ${DATA}/dst/mie/psd_CAM_LW.txt 2>&1 &
# ${HOME}/mie/psd.pl --dbg=1 --RGL --psd_nbr=4 --spc_idx={01,02,03,04} --sz_mnm={0.05,0.5,1.25,2.5} --sz_mxm={0.5,1.25,2.5,5.0} --sz_nbr={200,25,25,25} --dmt_vma_dfl=3.5 --wvl_mnm=0.625 --wvl_mxm=0.635 > ${DATA}/dst/mie/psd_ODX.txt 2>&1 &

# Dust version 5 optics: Never really used for anything but narrow band spectral studies
# ${HOME}/mie/psd.pl --dbg=1 --RGL --EMA=20060904 --vts_flg --dns_prt=2500.0 --psd_nbr=4 --spc_idx={01,02,03,04} --sz_mnm={0.05,0.5,1.25,5.0} --sz_mxm={0.5,1.25,5.0,10.0} --sz_nbr={200,25,25,25} --dmt_vma_dfl=3.5 --wvl_mnm=0.625 --wvl_mxm=0.635 > ${DATA}/dst/mie/psd_ODX.txt 2>&1 &

# Dust version 6 optics: Default as of ... (20130119: fxm still experimental)
# ${HOME}/mie/psd.pl --dbg=1 --RRTM_SW --EMA=20130119 --vts_flg --dns_prt=2500.0 --psd_nbr=4 --spc_idx={01,02,03,04} --sz_mnm={0.05,0.5,1.25,2.5} --sz_mxm={0.5,1.25,2.5,5.0} --sz_nbr={200,25,25,25} --dmt_vma_dfl=3.5 > ${DATA}/dst/mie/psd_RRTM_SW.txt 2>&1 &
# ${HOME}/mie/psd.pl --dbg=1 --RRTM_LW --EMA=20130119 --vts_flg --dns_prt=2500.0 --psd_nbr=4 --spc_idx={01,02,03,04} --sz_mnm={0.05,0.5,1.25,2.5} --sz_mxm={0.5,1.25,2.5,5.0} --sz_nbr={200,25,25,25} --dmt_vma_dfl=3.5 > ${DATA}/dst/mie/psd_RRTM_LW.txt 2>&1 &

# Stratiform scavenging with SeaWiFs channel
# ${HOME}/mie/psd.pl --dbg=1 --RGL --STR --ftn_fxd --psd_nbr=4 --spc_idx={01,02,03,04} --sz_mnm={0.05,0.5,1.25,2.5} --sz_mxm={0.5,1.25,2.5,5.0} --sz_nbr={200,25,25,25} --dmt_vma_dfl=3.5 --wvl_mnm=0.860 --wvl_mxm=0.870 > ${DATA}/dst/mie/psd_ODX.txt 2>&1 &

BEGIN{
    unshift @INC,$ENV{'HOME'}.'/perl'; # Location of csz.pl and DBG.pm HaS98 p. 170
} # end BEGIN

my $CVS_Header='$Id$';

# Specify modules
use strict; # Protect all namespaces
use Getopt::Long; # GNU-style getopt
use File::Basename; # For parsing filenames

# 3rd party modules
#use NetCDF; # man netCDFPerl

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
my ($fl_out_nc,$fl_out_txt,$fl_sfx);
my ($psd_idx,$wvl_grd,$wvl_sng);
my ($cmd_sng,$pcp_sng);
my ($cmp_sng_fl,$cmp_sng_arg,$dns_prt_arg,$mca_bln_arg);

# Set defaults 
my $False=0;
my $True=1;

my $CVS_Id='$Id$';
my $CVS_Revision='$Revision$';
my $CVS_Date='$Date$';
my $pi=atan2(1,1)*4.0;
my $data_nm=$ENV{'DATA'};
my $wvl_nbr=1; # [nbr] Number of output wavelengths (regular grids only)
my $wvl_mnm_mcr=0.495; # [um] Minimum wavelength (regular grids only)
my $wvl_mxm_mcr=0.505; # [um] Maximum wavelength (regular grids only)
my $dmt_pcp_nma_cnv_mcr=1000.0; # [um] Diameter number median analytic, raindrop, convective
my $dmt_pcp_nma_str_mcr=400.0; # [um] Diameter number median analytic, raindrop, stratiform
my $flx_vlm_pcp_cnv=2.77e-6; # [m s-1]=[m3 m-2 s-1] Precipitation volume flux, convective [10 mm hr-1]
my $flx_vlm_pcp_str=0.277e-6; # [m s-1]=[m3 m-2 s-1] Precipitation volume flux, stratiform [1 mm hr-1]

# Set defaults for command line arguments
#my $gsd_anl=2.2; # [frc] Geometric standard deviation analytic PaG77 p. 2080 Table 1 (Old dust version 1 optics)
my $Boolean=$True;
my $CAM_LW=$False;
my $CAM_SW=$False;
my $CNV=$False; # [flg] Convective precipitation characteristics
my $EMA=0; # [enm] Effective Medium Approximation (0=None,1=Binary, >1 = MCA)
my $GSFC=$False;
my $HELP=$False;
my $ODX=$False;
my $RGL=$False; # [flg] Regular wavelength grid
my $RRTM_LW=$False;
my $RRTM_SW=$False;
my $SNICAR=$False; # [flg] SNICAR wavelength grid
my $STR=$False; # Stratiform precipitation characteristics
my $SWNB=$False;
my $TOMS=$False;
my $asp_rat_hxg=1.5; # [frc] Hexagonal prism aspect ratio
my $bch_flg=$True; # [flg] Batch behavior
my $bnd_SW_LW=5.0e-6; # [m] Boundary between SW and LW weighting
my $bnd_nbr=20; # [frc] Number of sub-bands per output band
my $bnd_nbr_rgl=2; # [frc] Number of sub-bands per output band
my $cmp_mtx='SiO2_avg_hitran96'; # [sng] Composition of matrix
my $cmp_ncl='Fe2O3_doccd'; # [sng] Composition of inclusion
my $cmp_prt='saharan_dust'; # [sng] Composition of particle
my $dbg_lvl=0;
my $dmt_vma_dfl=3.5; # [um] Volume median diameter analytic, default RJM03 Table 1
my $dns_prt=0.0; # [kg m-3] Density of particle
#my $dmt_vma_dfl=2.524; # [um] Volume median diameter analytic, default She84 p. 75 Table 1 
#my $rds_nma_dfl=0.2986; # [um] Number median radius analytic, default She84 p. 75 Table 1 
my $dsd_mnm_mcr=1.0; # [um] Minimum diameter in raindrop distribution
my $dsd_mxm_mcr=4000.0; # [um] Maximum diameter in raindrop distribution
my $dsd_nbr=200; # [nbr] Number of raindrop size bins
my $fdg_val=1.0; # [frc] Tuning factor for all bands
my $ffc_mdm_typ='ffc_mdm_mxg'; # [enm] Effective medium type
my $flt_foo=73.37;
my $foo=73;
my $ftn_fxd_flg=$False; # [flg] Fortran fixed format
my $gsd_anl=2.0; # [frc] Geometric standard deviation analytic SBG98 p. 10581 Table 1, She84 p. 75 Table 1
my $gsd_pcp_anl=1.86; # [frc] Geometric standard deviation, raindrop NGD94 p. 2337 Table 1 
my $hrz_flg=$False; # [flg] Print size-resolved optical properties at debug wavelength
my $mss_frc_ncl=0.01; # [frc] Mass fraction in inclusion
my $mss_val=1.0e36; # [frc] Missing value
my $psd_nbr=1;
my $psd_typ='lognormal'; # [sng] Particle size distribution type
my $spc_idx='01'; # [sng] Species index for Fortran data
my $sz_grd='log';
my $sz_mnm=0.01; # [um] Minimum radius in aerosol distribution
my $sz_mxm=1.0; # [um] Maximum radius in aerosol distribution
my $sz_nbr=200; # [nbr] Number of aerosol size bins
my $vts_flg=$False; # [flg] Apply equal-V/S approximation for aspherical optical properties
my $wrn_ntp_flg=$False; # Print WARNINGs from ntp_vec()

my @bnd_nbr=(); # [frc] Number of sub-bands per output band
my @spc_idx=(); # [sng] Species index for Fortran data
my @foo=();
my @gsd_anl=(); # [frc] Geometric standard deviation analytic
my @rds_nma=(); # [um] Number median radius analytic
my @dmt_vma=(); # [um] Volume median diameter analytic
my @sz_mnm=(); # [um]
my @sz_mxm=(); # [um]
my @sz_nbr=(); # [nbr]

# Derived fields
my $pth_out=$data_nm.'/dst/mie'; # Path to output netCDF files

$prg_dsc='Dust Mie Processor'; # Program description
if(length($CVS_Id) > 4){($prg_nm,$prg_vrs)=$CVS_Id =~ /: (.+).pl,v ([\d.]+)/;}else{$prg_vrs='Unknown';}
$prg_vrs.='*' if length('$Locker:  $ ') > 12; # Tack '*' if it is not checked in into CVS
$prg_nm=$0; # $0 is program name Camel p. 136
if(length($CVS_Date) > 6){($prg_date)=unpack '@7 a19',$CVS_Date;}else{$prg_date='Unknown';}

# Parse command line arguments: '!' means Boolean, '|' is OR, '=' specifies required argument: 'i' is integer, 'f' is float, 's' is string
$rcd=Getopt::Long::Configure('no_ignore_case'); # Turn on case-sensitivity
$rcd=GetOptions( # man Getopt::GetoptLong
		'Boolean!' => \$Boolean,
		'CAM_LW!' => \$CAM_LW,
		'CAM_SW!' => \$CAM_SW,
		'RGL!' => \$RGL, # [flg] Regular wavelength grid
		'RRTM_LW!' => \$RRTM_LW,
		'RRTM_SW!' => \$RRTM_SW,
		'SNICAR!' => \$SNICAR, # [flg] SNICAR wavelength grid
		'CNV!' => \$CNV,
		'EMA=i' => \$EMA, # [flg] Effective Medium Approximation
		'GSFC!' => \$GSFC,
		'STR!' => \$STR,
		'SWNB!' => \$SWNB,
		'TOMS!' => \$TOMS,
		'asp_rat_hxg=f' => \$asp_rat_hxg, # [frc] Hexagonal prism aspect ratio
		'bch_flg!' => \$bch_flg, # [flg] Batch behavior
		'bnd_nbr=i' => \$bnd_nbr, # [frc] Number of sub-bands per output band
		'cmp_mtx=s' => \$cmp_mtx, # [sng] Composition of matrix
		'cmp_ncl=s' => \$cmp_ncl, # [sng] Composition of inclusion
		'cmp_prt=s' => \$cmp_prt, # [sng] Composition of particle
		'dbg_lvl=i' => \$dbg_lvl,
		'dmt_vma=f' => \@dmt_vma, # [um] Volume median diameter analytic
		'dmt_vma_dfl=f' => \$dmt_vma_dfl, # [um] Volume median diameter analytic
		'dns_prt=f' => \$dns_prt, # [kg m-3] Density of particle
		'fdg_val=f' => \$fdg_val, 
		'flt_foo=f' => \$flt_foo, 
		'ffc_mdm_typ=s' => \$ffc_mdm_typ, # [enm] Effective medium type
		'foo=i' => \@foo,
		'ftn_fxd!' => \$ftn_fxd_flg, # [flg] Fortran fixed format
		'gsd_anl=f' => \@gsd_anl, # [frc] Geometric standard deviation analytic
		'gsd_pcp_anl=f' => \$gsd_pcp_anl, # [frc] Geometric standard deviation, raindrop
		'help|usage|?', => \$HELP,
		'hrz_flg!' => \$hrz_flg, # [flg] Print size-resolved optical properties at debug wavelength
		'mss_frc_ncl=f' => \$mss_frc_ncl, # [frc] Mass fraction in inclusion
		'psd_nbr=i' => \$psd_nbr,
		'spc_idx=s' => \@spc_idx, # [sng] Species index for Fortran data
		'sz_mnm=f' => \@sz_mnm,
		'sz_mxm=f' => \@sz_mxm,
		'sz_nbr=i' => \@sz_nbr,
		'vts_flg!' => \$vts_flg, # [flg] Apply equal-V/S approximation for aspherical optical properties
		'wrn_ntp_flg!' => \$wrn_ntp_flg,
		'wvl_nbr=i' => \$wvl_nbr, # [nbr] Number of output wavelengths (regular grids only)
		'wvl_mnm_mcr=f' => \$wvl_mnm_mcr, # [um] Minimum wavelength (regular grids only)
		'wvl_mxm_mcr=f' => \$wvl_mxm_mcr, # [um] Maximum wavelength (regular grids only)
		 ); # end GetOptions arguments

if(!$STR && !$CNV){$CNV=$True;}elsif($STR && $CNV){die "$prg_nm: ERROR Both CNV and STR are set\n";}
if($HELP){&usg_prn();exit 0;} # end HELP

# Lists not specified on command line get default values here
# Initialize arrays with repetition operator, Camel p. 82
if(!@bnd_nbr){@bnd_nbr=($bnd_nbr)x$psd_nbr;} # [frc] Number of sub-bands per output band
if(!@gsd_anl){@gsd_anl=($gsd_anl)x$psd_nbr;} # [frc] Geometric standard deviation analytic
if(!@spc_idx){@spc_idx=($spc_idx)x$psd_nbr;} # [sng] Species index for Fortran data
if(!@sz_mnm){@sz_mnm=($sz_mnm)x$psd_nbr;}
if(!@sz_mxm){@sz_mxm=($sz_mxm)x$psd_nbr;}
if(!@sz_nbr){@sz_nbr=($sz_nbr)x$psd_nbr;}
if(!@foo){@foo=($foo)x$psd_nbr;}

# Derive rds_nma from mnm and mxm
if(!@dmt_vma){
    for ($psd_idx=0;$psd_idx<$psd_nbr;$psd_idx++){
#	$rds_nma[$psd_idx]=0.5*($sz_mnm[$psd_idx]+$sz_mxm[$psd_idx]); # [um] Number median radius analytic (Old dust version 2 optics)
	$dmt_vma[$psd_idx]=$dmt_vma_dfl; # [um] Volume median diameter analytic (New dust version 3 optics)
    } # end loop over $aer
} # endif !defined(@dmt_vma)

if($SNICAR){
    # SNICAR flag may be used as short-hand to specify following grid
    $RGL=$True;
    $wvl_nbr=470;
    $wvl_mnm_mcr=0.3;
    $wvl_mxm_mcr=5.0;
} # !SNICAR

if($CAM_SW){$wvl_sng='CAM_SW';$bnd_SW_LW=5.0e-6;}
# NB: I think I used fdg_val=1.22 for volcanic aerosol in LW code
# No known reason to use it for other aerosols
# It should probably go away as it is intrinsically bad
# elsif($CAM_LW){$wvl_sng='CAM_LW';$bnd_SW_LW=3.0e-6;$fdg_val=1.22;}
elsif($CAM_LW){$wvl_sng='CAM_LW';$bnd_SW_LW=3.0e-6;}
elsif($GSFC){$wvl_sng='GSFC_LW';$bnd_SW_LW=3.0e-6;}
elsif($RGL){$wvl_sng='rgl';$bnd_SW_LW=5.0e-6;}
elsif($RRTM_SW){$wvl_sng='RRTM_SW';$bnd_SW_LW=5.0e-6;}
elsif($RRTM_LW){$wvl_sng='RRTM_LW';$bnd_SW_LW=3.0e-6;}
elsif($SWNB){$wvl_sng='SWNB';$bnd_SW_LW=5.0e-6;}
elsif($TOMS){$wvl_sng='TOMS';$bnd_SW_LW=5.0e-6;}

$pcp_sng='--dsd_mnm='.$dsd_mnm_mcr.' --dsd_mxm='.$dsd_mxm_mcr.' --dsd_nbr='.$dsd_nbr.' --gsd_pcp_anl='.$gsd_pcp_anl;
if($CNV){$pcp_sng.=' --dmt_pcp_nma='.$dmt_pcp_nma_cnv_mcr.' --flx_vlm_pcp='.$flx_vlm_pcp_cnv;}else
{$pcp_sng.=' --dmt_pcp_nma='.$dmt_pcp_nma_str_mcr.' --flx_vlm_pcp='.$flx_vlm_pcp_str;}

# Print initialization state
if($dbg_lvl > 1){print ("$prg_nm: $prg_dsc, version $prg_vrs of $prg_date\n");} # endif dbg
if($dbg_lvl >= $dbg_fl){
    printf ("Initialization State:\n");
    if($EMA == 1){printf ("cmp_mtx = $cmp_mtx and cmp_ncl = $cmp_ncl\n");}else{printf ("cmp_prt = $cmp_prt\n");}
    printf ("wvl_mnm = $wvl_mnm_mcr um, wvl_mxm = $wvl_mxm_mcr um\n");
    printf ("bnd_SW_LW = $bnd_SW_LW um\n");
    printf ("fdg_val = $fdg_val\n");
    printf ("psd_idx\tspc_idx\tsz_mnm\tsz_mxm\tdmt_vma\tgsd_anl\tsz_nbr\tbnd_nbr\n");
    printf ("#\t\tum\tum\tum\t\t\t\n");
    for ($psd_idx=0;$psd_idx<$psd_nbr;$psd_idx++){
	print ("$psd_idx\t$spc_idx[$psd_idx]\t$sz_mnm[$psd_idx]\t$sz_mxm[$psd_idx]\t$dmt_vma[$psd_idx]\t$gsd_anl[$psd_idx]\t$sz_nbr[$psd_idx]\t$bnd_nbr[$psd_idx]\n");
    } # end loop over $aer
    print ("\n");
} # endif dbg

if($fdg_val != 1.0){print ("$prg_nm: WARNING Using non-unity \$fdg_val = $fdg_val. This value will multiply all computed extinction coefficients. This is intended to tune optical properties on coarse grids (e.g., CAM LW) to obtain reasonable agreement with high resolution grids. If current optical properties are not for CAM LW in CAM/CRM or this \$fdg_val was not tuned by you personally then be careful using these results!\n");}

if($dbg_lvl == $dbg_crr){
    print ('@foo = ');
    foreach $foo (@foo){
	print ("$foo, ");
    } # end loop over $foo
    print ("\n");
} # endif dbg

if($EMA == 0){
    # Single component aerosols
    $cmp_sng_fl=$cmp_prt;
    $cmp_sng_arg='--cmp_prt='.$cmp_prt;
}else{ # EMA != 0
    # Multi-component aerosols handled in two stages
    # First, particular "named" blends, which are usually fragile MCAs
    # Then easier binary mixtures specifiable by, e.g., matrix/inclusion terminology
    if($EMA == 1){
	$cmp_sng_fl=$cmp_mtx.'_'.$cmp_ncl;
	$cmp_sng_arg='--ffc_mdm_typ='.$ffc_mdm_typ.' --cmp_mtx='.$cmp_mtx.' --cmp_ncl='.$cmp_ncl.' --mss_frc_ncl='.$mss_frc_ncl;
    }elsif($EMA == 20060621){
	$mca_bln_arg='--ncl=0.01 --ncl=0.25 --ncl=0.25 --ncl=0.005 --ffc_mdm_typ=ffc_mdm_mxg --cmp_mtx=SiO2_avg_hitran96 --cmp_cor=limestone_krk --cmp_mnt=montmorillonite --cmp_ncl=illite --cmp_prt=Fe2O3_doccd';
	$cmp_sng_fl='dst_bln_20060621';
	$cmp_sng_arg='--ffc_mdm_typ='.$ffc_mdm_typ.' '.$mca_bln_arg;
    }elsif($EMA == 20060904){
	$mca_bln_arg='--ncl=0.02 --ncl=0.25 --ncl=0.25 --ncl=0.004 --ffc_mdm_typ=ffc_mdm_mxg --cmp_mtx=SiO2_avg_hitran96 --cmp_cor=limestone_krk --cmp_mnt=montmorillonite --cmp_ncl=illite --cmp_prt=Fe2O3_avg_roush';
	$cmp_sng_fl='dst_bln_20060904';
	$cmp_sng_arg='--ffc_mdm_typ='.$ffc_mdm_typ.' '.$mca_bln_arg;
    }elsif($EMA == 20060907){
	$mca_bln_arg='--ncl=0.0014 --ncl=0.25 --ncl=0.25 --ncl=0.004 --ffc_mdm_typ=ffc_mdm_mxg --cmp_mtx=SiO2_avg_hitran96 --cmp_cor=lac_ChC90 --cmp_mnt=montmorillonite --cmp_ncl=illite --cmp_prt=Fe2O3_avg_roush';
	$cmp_sng_fl='dst_bln_20060907';
	$cmp_sng_arg='--ffc_mdm_typ='.$ffc_mdm_typ.' '.$mca_bln_arg;
    }elsif($EMA == 20130119){
	$mca_bln_arg='--ncl=0.02 --ncl=0.25 --ncl=0.25 --ncl=0.000 --ffc_mdm_typ=ffc_mdm_mxg --cmp_mtx=SiO2_avg_hitran96 --cmp_cor=limestone_krk --cmp_mnt=montmorillonite --cmp_ncl=illite --cmp_prt=Fe2O3_avg_roush';
	$cmp_sng_fl='dst_bln_20130119';
	$cmp_sng_arg='--ffc_mdm_typ='.$ffc_mdm_typ.' '.$mca_bln_arg;
    }else{ # !EMA_mlt
	die "$prg_nm: ERROR EMA type unknown\n";
    } # !EMA_mlt
} # EMA != 0

# Allow manual density specification to override computed EMA density
if($dns_prt == 0.0){$dns_prt_arg='';}else{$dns_prt_arg="--dns_prt=$dns_prt ";}
$cmp_sng_arg=$dns_prt_arg.$cmp_sng_arg;

for ($psd_idx=0;$psd_idx<$psd_nbr;$psd_idx++){

# Number of size bins used for each distribution decreases with increasing size 
# This is due to computational expense of large size parameters (2*pi*rds/lambda)
# Accuracy is not compromised because optical importance decreases even faster than discretization error

    if($CAM_SW || $CAM_LW || $GSFC || $RRTM_SW || $RRTM_LW || $SWNB || $TOMS){
# Compute broadband optical properties for solar or terrestrial spectral regimes
	$fl_out_nc=$pth_out.'/aer_'.$cmp_sng_fl.'_'.$spc_idx[$psd_idx].'_'.$wvl_sng.'.nc';
	$fl_out_txt=$pth_out.'/aer_'.$cmp_sng_fl.'_'.$spc_idx[$psd_idx].'_'.$wvl_sng.'.txt';
	$cmd_sng="--dbg=$dbg_lvl --spc_idx=$spc_idx[$psd_idx] $cmp_sng_arg --psd_typ=$psd_typ --bnd_nbr=$bnd_nbr[$psd_idx] --wvl_grd=$wvl_sng --bnd_SW_LW=$bnd_SW_LW --sz_grd=$sz_grd --sz_mnm=$sz_mnm[$psd_idx] --sz_mxm=$sz_mxm[$psd_idx] --sz_nbr=$sz_nbr[$psd_idx] --dmt_vma=$dmt_vma[$psd_idx] --gsd_anl=$gsd_anl[$psd_idx] $pcp_sng $fl_out_nc";
        if($bch_flg){$cmd_sng="--bch ".$cmd_sng;}
        if($ftn_fxd_flg){$cmd_sng="--ftn_fxd ".$cmd_sng;}
        if($hrz_flg){$cmd_sng="--hrz ".$cmd_sng;}
        if($vts_flg){$cmd_sng="--vts_flg --asp_rat_hxg=$asp_rat_hxg ".$cmd_sng;}
        if(!$wrn_ntp_flg){$cmd_sng="--no_wrn_ntp ".$cmd_sng;}
	if($fdg_val != 1.0){$cmd_sng="--fdg_val=$fdg_val ".$cmd_sng;}
	if($dbg_lvl < $dbg_io){&cmd_prc('mie '.$cmd_sng);}else{print 'mie '.$cmd_sng."\n";}
    } # endif $CAM_SW || $CAM_LW || $GSFC || $RRTM_SW || $RRTM_LW || || $SWNB || $TOMS

    if($RGL){ # Regular (i.e., evenly spaced) Wavelength Grids
# Compute narrowband optical properties for diagnostic optical depths
# The wavelength grid for this command differs from the above commands
	if($wvl_nbr==470 && $wvl_mnm_mcr==0.3 && $wvl_mxm_mcr==5.0){$SNICAR=$True;}
	if($wvl_nbr==1){$ODX=$True;}
	if($ODX){
	    $bnd_nbr_rgl=2;
	    $wvl_sng=sprintf('%03d',50*($wvl_mnm_mcr+$wvl_mxm_mcr));
	} # !ODX
	if($SNICAR){
	    $bnd_nbr_rgl=1;
	    $wvl_sng='snicar';
	} # !SNICAR
	$fl_out_nc=$pth_out.'/aer_'.$cmp_sng_fl.'_'.$spc_idx[$psd_idx].'_'.$wvl_sng.'.nc';
	$fl_out_txt=$pth_out.'/aer_'.$cmp_sng_fl.'_'.$spc_idx[$psd_idx].'_'.$wvl_sng.'.txt';
	$cmd_sng="--dbg=$dbg_lvl --spc_idx=$spc_idx[$psd_idx] $cmp_sng_arg --psd_typ=$psd_typ --wvl_grd=reg --wvl_nbr=$wvl_nbr --wvl_mnm=$wvl_mnm_mcr --wvl_mxm=$wvl_mxm_mcr --bnd_nbr=$bnd_nbr_rgl --sz_grd=$sz_grd --sz_mnm=$sz_mnm[$psd_idx] --sz_mxm=$sz_mxm[$psd_idx] --sz_nbr=$sz_nbr[$psd_idx] --dmt_vma=$dmt_vma[$psd_idx] --gsd_anl=$gsd_anl[$psd_idx] $pcp_sng $fl_out_nc";
	if($bch_flg){$cmd_sng="--bch ".$cmd_sng;}
	if($ftn_fxd_flg){$cmd_sng="--ftn_fxd ".$cmd_sng;}
	if($hrz_flg){$cmd_sng="--hrz ".$cmd_sng;}
        if($vts_flg){$cmd_sng="--vts_flg --asp_rat_hxg=$asp_rat_hxg ".$cmd_sng;}
	if(!$wrn_ntp_flg){$cmd_sng="--no_wrn_ntp ".$cmd_sng;}
	if($fdg_val != 1.0){$cmd_sng="--fdg_val=$fdg_val ".$cmd_sng;}
	if($dbg_lvl < $dbg_io){&cmd_prc('mie '.$cmd_sng);}else{print 'mie '.$cmd_sng."\n";}
    } # endif $RGL

} # end loop over aer

# Get elapsed times
time_end($srt_usr_tm,$srt_sys_tm,$srt_child_usr_tm,$srt_child_sys_tm);

