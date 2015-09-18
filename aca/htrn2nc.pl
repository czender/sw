#!/usr/bin/perl

# Purpose: Extract line info from HITRAN database, save in netCDF format

# Usage:

# Production:
# Shell script $HOME/aca/htrn2nc.sh processes most useful gases

# Debugging:
# head -50 ${DATA}/hitran/hitran12.txt > ${DATA}/hitran/foo.txt
# $HOME/aca/htrn2nc.pl --mlc={1,7} --dbg=1 ${DATA}/hitran/foo.txt ${DATA}/hitran/foo.nc
# $HOME/aca/htrn2nc.pl --mlc=7 --iso=1 --dbg=1 ${DATA}/hitran/foo.txt ${DATA}/hitran/foo.nc
# $HOME/aca/htrn2nc.pl --dbg=1 --wvn_min=17900.0 --wvn_max=30000.0 ${DATA}/hitran/hitran12.txt ${DATA}/hitran/foo.nc
# $HOME/aca/htrn2nc.pl --dbg=1 --wvn_min=1931.0 --wvn_max=1939.0 ${DATA}/hitran/hitran12.txt ${DATA}/hitran/foo.nc

# LBL validation:
# GoY89 p. 120 fig 3.25, the Q-branch of the 3v2 band of CO2:
# $HOME/aca/htrn2nc.pl --mlc=2 --wvn_min=1931.0 --wvn_max=1939.0 ${DATA}/hitran/hitran12.txt ${DATA}/hitran/lbl_CO2_1931_1939.nc

# Water Vapor
# $HOME/aca/htrn2nc.pl --mlc=1 --wvn_min=2000.0 --wvn_max=17900.0 ${DATA}/hitran/hitran12.txt ${DATA}/hitran/H2O.nc
# $HOME/aca/htrn2nc.pl --mlc=1 --wvn_min=1.0e-4 --wvn_max=50000.0 /cdrom/hitran96/hitran96.par ${DATA}/hitran/H2O_all.nc
# $HOME/aca/htrn2nc.pl --mlc=1 --iso=1 --wvn_min=2000.0 --wvn_max=17900.0 ${DATA}/hitran/hitran12.txt ${DATA}/hitran/1H2_16O.nc
# $HOME/aca/htrn2nc.pl --mlc=1 --iso=2 --wvn_min=2000.0 --wvn_max=17900.0 ${DATA}/hitran/hitran12.txt ${DATA}/hitran/1H2_18O.nc
# $HOME/aca/htrn2nc.pl --mlc=1 --iso=3 --wvn_min=2000.0 --wvn_max=17900.0 ${DATA}/hitran/hitran12.txt ${DATA}/hitran/1H2_17O.nc
# $HOME/aca/htrn2nc.pl --mlc=1 --iso=4 --wvn_min=2000.0 --wvn_max=17900.0 ${DATA}/hitran/hitran12.txt ${DATA}/hitran/1H_2H_16O.nc

BEGIN{
    unshift @INC,$ENV{'HOME'}.'/perl'; # Location of DBG.pm HaS98 p. 170
    unshift @INC,$ENV{'HOME'}.'/ck'; # Location of HITRAN.pm HaS98 p. 170
} # end BEGIN

my $CVS_Header='$Id$';

# Specify modules
use strict; # Protect all namespaces
use Getopt::Long; # GNU-style getopt
use File::Basename; # For parsing filenames

# 3rd party modules
use NetCDF; # man netCDFPerl

# Personal modules
use HITRAN; # HITRAN database
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
my ($fl_in,$fl_out,$fl_sfx);
my ($pth_in,$pth_out);
my ($mlc_nm,@foo,$mlc_nbr,$iso_nbr,$mlc_sng);

# Declare netCDF output fields
my ($ln_mlc,$ln_iso,$ln_ctr,$ln_str,$HWHM_air,$HWHM_slf,$ln_nrg_lwr,$HWHM_tpt_dpn_xpn,$ln_sft_air,$qnt_glb_lwr,$qnt_glb_upr,$qnt_lcl_lwr,$qnt_lcl_upr);

# Declare netCDF metadata variables
my ($nc_id,@scalar);
my ($dim_nm,%dim_cnt,%dim_srt,%dim_id,%dim_vec);
my (%var_dim,%var_srt,%var_cnt,%var_id,$var_nm,%var_att,%var_type);
my (%att_id,$att_type,$att_nm,$att_idx,$att_val,$att_id);

# Set defaults 
my $False=0;
my $True=1;

my $CVS_Date='$Date$';
my $CVS_Id='$Id$';
my $CVS_Revision='$Revision$';
my $HELP=$False;
my $data_nm=$ENV{'DATA'};
my $mss_val_out=1.0e36; # Missing value in output data
my $pi=atan2(1.0,1.0)*4.0;
my $usr_nm=$ENV{'USER'};
my ($mlc,$mlc_idx,$mlc_id,$iso,$iso_idx);
my (@crr_iso,$crr_iso,$iso_sng,$iso_nm);
my ($ln_cnt,$rx_sng,@swap);
my (@iso_crr_htrn,$iso_nbr_crr_htrn,$iso_nbr_crr_usr,@iso_crr_usr);

# Set defaults for command line arguments
my $Boolean=$True;
my $dbg_lvl=0;
my $wvn_min=1.0e-36; # [cm-1] Minimum wavenumber
my $wvn_max=50000.0; # [cm-1] Maximum wavenumber
my @mlc=();
my @iso=();
my $float_foo=102.543;
				 
# Derived fields

# Package name
my $my_package='HITRAN to netCDF converter';
# Program name and version
my ($my_name,$my_version)= $CVS_Id =~ /: (.+).pl,v ([\d.]+)/;
# Tack on '*' if it is not checked in into CVS
$my_version.='*' if length('$Locker:  $ ') > 12;

$prg_dsc='HITRAN to netCDF converter'; # Program description
($prg_nm,$prg_vrs)=$CVS_Id =~ /: (.+).pl,v ([\d.]+)/; # Program name and version
$prg_vrs.='*' if length('$Locker:  $ ') > 12; # Tack '*' if it is not checked in into CVS.
($prg_nm,$pth_in,$fl_sfx)=fileparse($0,''); # $0 is program name Camel p. 136.
if(length($CVS_Date) > 6){($prg_date)=unpack '@7 a19',$CVS_Date;}else{$prg_date='Unknown';}

# Parse command line arguments: '!' means Boolean, '|' is OR, '=' specifies required argument: 'i' is integer, 'f' is float, 's' is string
$rcd=Getopt::Long::Configure('no_ignore_case'); # Turn on case-sensitivity
$rcd=GetOptions( # man Getopt::GetoptLong
		   'Boolean!' => \$Boolean,
		   'dbg_lvl=i' => \$dbg_lvl,
		   'float_foo=f' => \$float_foo, 
		   'help|usage', => \$HELP,
                   'wvn_max=f', => \$wvn_max,
                   'wvn_min=f' => \$wvn_min,
		   'iso=i' => \@iso,
		   'mlc=i' => \@mlc,
                   'foo=i' => \@foo,
		 ); # end GetOptions arguments

if($HELP){&usg_prn();exit 0;} # end HELP

# Parse positional arguments, if any
if($#ARGV > 1){die "$prg_nm: ERROR Called with $#ARGV+1 positional arguments, need no more than 2\n";}
elsif($#ARGV == 1){
# Input file name is first positional argument, if any
# Output file name is last positional argument, if any
    $fl_in=$ARGV[0];
    $fl_out=$ARGV[1];}
elsif($#ARGV == 0){$fl_out=$ARGV[0];}

# Definitions that depend on command line input
if($wvn_min == 0.0){
    print ("$prg_nm: HINT Use, e.g., --wvn_min=1.0e-36 to avoid divide by zero\n");
    die "$prg_nm: ERROR wvn_min = $wvn_min cm-1\n";
} # endif dbg
my $wvl_min=1.0/(100.0*$wvn_max); # [m] Minimum wavelength
my $wvl_max=1.0/(100.0*$wvn_min); # [m] Maximum wavelength

# Print initialization state
if($dbg_lvl > 1){print ("$prg_nm: $prg_dsc, version $prg_vrs of $prg_date\n");} # endif dbg

# fxm: Should add a function to unique'ize iso and mlc lists

# Check lists for errors
if(defined(@mlc)){
    foreach $mlc (@mlc){
	die "ERROR Invalid molecule ID $mlc\n" if($mlc > $mlc_nbr_max_htrn || $mlc < 1);
    } # end loop over $mlc
} # endif
if(defined(@iso)){
    foreach $iso (@iso){
	die "ERROR Invalid isotopomer ID $iso\n" if($iso > $iso_nbr_max_htrn || $iso < 1);
    } # end loop over $iso
} # endif

if(defined(@mlc) && defined(@iso)){
# Specifying both molecules and isotopes is currently only OK when mlc_nbr == 1
    if($#mlc+1 == 1){
# The isotope list is assumed to contain isotopic abundance indices (1..8) when only one molecule is specified
# Change to the isotope list into an isotopomer list
	@swap=@iso;
	@iso=();
	foreach $iso (@swap){
	    die "ERROR Invalid isotope ID $iso for molecule $mlc_sng{$mlc[0]}\n" if($iso >  $#{$mlc_iso{$mlc_sng{$mlc[0]}}}+1 || $iso < 1);
	    push @iso,$mlc_iso{$mlc_sng{$mlc[0]}}[$iso-1];
	} # end loop over swap
    }else{
	die "Unable to parse arbitrary molecule and isotope combinations at this time\n";
    } # endif
} # endif

# Ensure molecule list is defined
if(!defined(@mlc) && !defined(@iso)){@mlc=(1..$mlc_nbr_max_htrn);}
if(!defined(@mlc) && defined(@iso)){
# Isotope list is assumed to contain isotopomer indices (1..90) when there is no molecule list
    foreach $iso (@iso){
	for($mlc_idx=1;$mlc_idx<=$mlc_nbr_max_htrn;$mlc_idx++){
            $mlc_nm=$mlc_sng{$mlc_idx};
            for $idx (0..$#{$mlc_iso{$mlc_nm}}){
		if($iso == $mlc_iso{$mlc_nm}[$idx]){push @mlc,$mlc_idx;}
            } # end loop over $idx
	} # end loop over $mlc
    } # end loop over $iso
} # endif
$mlc_nbr=$#mlc+1;

# Ensure isotope list is defined
if(!defined(@iso)){
# Not specifying any isotopes is same as requesting all isotopomers of each molecule in molecule list
    foreach $mlc_nm (@mlc_sng{@mlc}){
	push @iso,@{$mlc_iso{$mlc_nm}}; # Put braces around list to extract it from an HoL
    } # end loop over $mlc_nm
} # endif
$iso_nbr=$#iso+1;

# Now that molecule and isotope lists are complete, sort them
@iso=sort numerically @iso; # Camel book, p. 218
@mlc=sort numerically @mlc; # Camel book, p. 218

# Print list of all HITRAN isotopes for each user-specified molecule
if($dbg_lvl >= 1){		 
    print ("All HITRAN isotopes for each user-specified molecule:\n");
    foreach $mlc_nm (@mlc_sng{@mlc}){
	@crr_iso=@{$mlc_iso{$mlc_nm}}; # Put braces around list to extract it from an HoL
	print ("Molecule $mlc_nm: ");
	foreach $crr_iso (@crr_iso){print ("$iso_sng{$crr_iso} ($crr_iso), ");}
	print ("\n");
    } # end loop over $mlc_nm
    print ("\n");
} # endif dbg

# Print user-specified molecule and isotope lists
if($dbg_lvl >= 1){		 
    print ("User-specified molecule and isotope lists:\n");
    print ("\$mlc_nbr = $mlc_nbr, \$\#mlc = $#mlc\n");
    for $idx (0..$#mlc){print ("\$mlc[$idx] = $mlc[$idx] = $mlc_sng{$mlc[$idx]}\n");}
    print ("\$iso_nbr = $iso_nbr, \$\#iso = $#iso\n");
    for $idx (0..$#iso){print ("\$iso[$idx] = $iso[$idx] = $iso_sng{$iso[$idx]}\n");}
    print ("\n");
} # endif dbg

# Create strings containing all molecules and isotopes requested
if($mlc_nbr == $mlc_nbr_max_htrn){
    $mlc_sng='All';
    $iso_sng='All';
}else{
    $iso_sng='';
    $mlc_sng='';
    foreach $mlc_nm (@mlc_sng{@mlc}){
	$mlc_sng.=$mlc_nm;
	if($mlc_nm ne $mlc_sng{$mlc[$#mlc]}){$mlc_sng.=' ';}
    } # end loop over $mlc
    foreach $iso_nm (@iso_sng{@iso}){
	$iso_sng.=$iso_nm;
	if($iso_nm ne $iso_sng{$iso[$#iso]}){$iso_sng.=' ';}
    } # end loop over $iso
} # endelse
print ("\$mlc_sng = $mlc_sng\n");
print ("\$iso_sng = $iso_sng\n\n");

# Print initialization state
if($dbg_lvl >= 1){		 
    printf STDOUT ("Reading input HITRAN data from %s\n",$fl_in);
    printf STDOUT ("Writing output netCDF data to %s\n",$fl_out);
    printf STDOUT ("rcd = %d, Boolean = %d, dbg_lvl = %d, float_foo = %f, mlc_nbr = %d, iso_nbr=%d, wvn_min = %f, wvn_max = %f, HELP = %d, mlc_sng = %s, iso_sng = %s\n",$rcd,$Boolean,$dbg_lvl,$float_foo,$mlc_nbr,$iso_nbr,$wvn_min,$wvn_max,$HELP,$mlc_sng,$iso_sng);
} # end if dbg 

# Dimensions
%dim_cnt=( # Hash of scalars
	   'ln_ctr',NetCDF::UNLIMITED,
	   'mlc',$#mlc+1,
	   'iso',$#iso+1,
	   'qnt_lcl',9,
	   );
%dim_srt=( # Hash of scalars
	   'ln_ctr',0,
	   'mlc',0,
	   'iso',0,
	   'qnt_lcl',0,
	   );

# Variables
%var_dim=( # Hash of lists
	   'HWHM_air',['ln_ctr'],
	   'HWHM_slf',['ln_ctr'],
	   'HWHM_tpt_dpn_xpn',['ln_ctr'],
	   'ln_ctr',['ln_ctr'],
	   'ln_iso',['ln_ctr'],
	   'ln_mlc',['ln_ctr'],
	   'ln_nrg_lwr',['ln_ctr'],
	   'ln_sft_air',['ln_ctr'],
	   'ln_str',['ln_ctr'],
	   'qnt_glb_lwr',['ln_ctr'],
	   'qnt_glb_upr',['ln_ctr'],
	   'qnt_lcl_lwr',['ln_ctr','qnt_lcl'],
	   'qnt_lcl_upr',['ln_ctr','qnt_lcl'],
	   'mlc_id',['mlc'],
	   'iso_id',['iso'],
	   );
%var_type=( # Hash of scalars
	    'HWHM_air',NetCDF::DOUBLE,
	    'HWHM_slf',NetCDF::DOUBLE,
	    'HWHM_tpt_dpn_xpn',NetCDF::DOUBLE,
	    'ln_ctr',NetCDF::DOUBLE,
	    'ln_iso',NetCDF::SHORT,
	    'ln_mlc',NetCDF::SHORT,
	    'ln_nrg_lwr',NetCDF::DOUBLE,
	    'ln_sft_air',NetCDF::DOUBLE,
	    'ln_str',NetCDF::DOUBLE,
	    'qnt_glb_lwr',NetCDF::SHORT,
	    'qnt_glb_upr',NetCDF::SHORT,
	    'qnt_lcl_lwr',NetCDF::CHAR,
	    'qnt_lcl_upr',NetCDF::CHAR,
	    'mlc_id',NetCDF::SHORT,
	    'iso_id',NetCDF::SHORT,
	    );
%var_att=( # Hash of lists of lists
	   'HWHM_air',[
		       ['long_name',NetCDF::CHAR,'Air-broadened HWHM'],
		       ['units',NetCDF::CHAR,'cm-1 atm-1 @ 296 K'],
		       ],
	   'HWHM_slf',[
		       ['long_name',NetCDF::CHAR,'Self-broadened HWHM'],
		       ['units',NetCDF::CHAR,'cm-1 atm-1 @ 296 K'],
		       ],
	   'HWHM_tpt_dpn_xpn',[
			       ['long_name',NetCDF::CHAR,'Exponent of temperature dependence of air-broadened HWHM'],
			       ['units',NetCDF::CHAR,'fraction'],
			       ],
	   'ln_ctr',[
		     ['long_name',NetCDF::CHAR,'Wavenumber at line center'],
		     ['units',NetCDF::CHAR,'cm-1'],
		     ],
	   'ln_iso',[
		     ['long_name',NetCDF::CHAR,'HITRAN isotope number (1..9)'],
		     ['units',NetCDF::CHAR,'number'],
		     ],
	   'ln_mlc',[
		     ['long_name',NetCDF::CHAR,'HITRAN molecule number (1..37)'],
		     ['units',NetCDF::CHAR,'number'],
		     ],
	   'ln_nrg_lwr',[
			['long_name',NetCDF::CHAR,'Lower state energy of transition'],
			['units',NetCDF::CHAR,'cm-1'],
			],
	   'ln_sft_air',[
			 ['long_name',NetCDF::CHAR,'Air-broadened pressure shift of line transition'],
			 ['units',NetCDF::CHAR,'cm-1 atm-1 @ 296 K'],
			 ],
	   'ln_str',[
		     ['long_name',NetCDF::CHAR,'Line strength'],
		     ['units',NetCDF::CHAR,'cm-1 molecule-1 cm2 @ 296 K'],
		     ],
	   'qnt_glb_lwr',[
			  ['long_name',NetCDF::CHAR,'Lower state global quanta index (vibrational quanta)'],
			  ['units',NetCDF::CHAR,'key'],
			  ],
	   'qnt_glb_upr',[
			  ['long_name',NetCDF::CHAR,'Upper state global quanta index (vibrational quanta)'],
			  ['units',NetCDF::CHAR,'key'],
			  ],
	   'qnt_lcl_lwr',[
		     ['long_name',NetCDF::CHAR,'Lower state local quanta index (rotational quanta)'],
		     ['units',NetCDF::CHAR,'key'],
		     ],
	   'qnt_lcl_upr',[
		     ['long_name',NetCDF::CHAR,'Upper state local quanta index (rotational quanta)'],
		     ['units',NetCDF::CHAR,'key'],
		     ],
	   'mlc_id',[
		     ['long_name',NetCDF::CHAR,'HITRAN Molecule IDs in this file'],
		     ['units',NetCDF::CHAR,'number'],
		     ['valid_range',NetCDF::SHORT,[1,$mlc_nbr_max_htrn]],
		     ],
	   'iso_id',[
		     ['long_name',NetCDF::CHAR,'HITRAN Isotopomer IDs in this file'],
		     ['units',NetCDF::CHAR,'number'],
		     ['valid_range',NetCDF::SHORT,[1,$iso_nbr_max_htrn]],
		     ],
	   ); # end var_att

# Create netCDF output file
print STDOUT 'Creating netCDF file...';
$nc_id=NetCDF::create($fl_out,NetCDF::CLOBBER);
die "Could not open netCDF file\n" if $nc_id < 0;
print STDOUT "ok\n";

# Put global attributes
print STDERR "Writing global attributes...";
$att_id=NetCDF::attput($nc_id,NetCDF::GLOBAL,"creation_date",NetCDF::CHAR,$lcl_date_time);
$att_id=NetCDF::attput($nc_id,NetCDF::GLOBAL,"CVS_Header",NetCDF::CHAR,$CVS_Header);
$att_id=NetCDF::attput($nc_id,NetCDF::GLOBAL,"molecule",NetCDF::CHAR,$mlc_sng);
$att_id=NetCDF::attput($nc_id,NetCDF::GLOBAL,"isotope",NetCDF::CHAR,$iso_sng);
$att_id=NetCDF::attput($nc_id,NetCDF::GLOBAL,"wvn_min_rgn",NetCDF::DOUBLE,$wvn_min);
$att_id=NetCDF::attput($nc_id,NetCDF::GLOBAL,"wvn_max_rgn",NetCDF::DOUBLE,$wvn_max);
$att_id=NetCDF::attput($nc_id,NetCDF::GLOBAL,"wvl_min_rgn",NetCDF::DOUBLE,$wvl_min);
$att_id=NetCDF::attput($nc_id,NetCDF::GLOBAL,"wvl_max_rgn",NetCDF::DOUBLE,$wvl_max);
print STDERR "ok\n";

# Set fill mode
#print STDOUT 'Setting fill mode...';
#NetCDF::setfill($nc_id,NetCDF::NOFILL) == 0 || die "Could not set fill mode\n";
#print STDOUT "ok\n";

if($dbg_lvl >= 3){		 
    print STDOUT "var_dim:\n";
    foreach $var_nm (keys %var_dim){
	print STDOUT "\t$var_nm: @{$var_dim{$var_nm}}\n";
    } # end loop over variables
    print "\n";
} # end if dbg 

# Define dimensions
print STDOUT 'Defining dimensions...';
foreach $dim_nm (keys(%dim_cnt)){
    $dim_id{$dim_nm}=NetCDF::dimdef($nc_id,$dim_nm,$dim_cnt{$dim_nm});
    die "Could not define dimension\n" if $dim_id{$dim_nm} < 0;
} # end loop over dimensions
print STDOUT "ok\n";

if($dbg_lvl >= 3){		 
    print STDOUT "dim_id:\n";
    foreach $dim_nm (keys %dim_id){
	print STDOUT "\t$dim_nm: $dim_id{$dim_nm}\n";
    } # end loop over dimensions
} # end if dbg 

# Now that we have dimension IDs, form dimension ID vectors
@scalar=();
foreach $var_nm (sort keys(%var_dim)){
    if(defined($var_dim{$var_nm}[0])){
	foreach $dim_nm (@{$var_dim{$var_nm}}){
	    push @{$dim_vec{$var_nm}},$dim_id{$dim_nm};
	    push @{$var_srt{$var_nm}},$dim_srt{$dim_nm};
	    push @{$var_cnt{$var_nm}},$dim_cnt{$dim_nm};
	} # end loop over dimensions
    }else{
	$dim_vec{$var_nm}=(\@scalar);
    } # endif
} # end loop over variables

if($dbg_lvl >= 3){		 
    print STDOUT "dim_vec:\n";
    foreach $var_nm (keys %dim_vec){
	print STDOUT "\t$var_nm: @{$dim_vec{$var_nm}}\n";
	print STDOUT "\t$var_nm: @{$var_srt{$var_nm}}\n";
	print STDOUT "\t$var_nm: @{$var_cnt{$var_nm}}\n";
    } # end loop over variables
    print "\n";
} # end if dbg 

# Define variables in alphabetical order so we know order to use in recput() call without using recinq()
print STDOUT 'Defining variables...';
foreach $var_nm (sort keys(%var_dim)){
    
    $var_id{$var_nm}=NetCDF::vardef($nc_id,$var_nm,$var_type{$var_nm},\$dim_vec{$var_nm});
    die "Could not define variable\n" if $var_id{$var_nm} < 0;
    
    for $att_idx (0..$#{$var_att{$var_nm}}){
	$att_nm=$var_att{$var_nm}[$att_idx][0];
	$att_type=$var_att{$var_nm}[$att_idx][1];
	$att_val=$var_att{$var_nm}[$att_idx][2];
	push @{$att_id{$var_nm}},NetCDF::attput($nc_id,$var_id{$var_nm},$att_nm,$att_type,$att_val).',';
	die "Could not put attribute\n" if $att_id{$var_nm} < 0;
    } # end loop over attributes of variable
    
} # end loop over variables
print STDOUT "ok\n";

if($dbg_lvl >= 3){		 
    print STDOUT "att_id:\n";
    foreach $var_nm (keys %att_id){
	print STDOUT "\t$var_nm: @{$att_id{$var_nm}}\n";
    } # end loop over variables
    print "\n";
} # end if dbg 

# End definition
print STDOUT 'Ending definition mode...';
$rcd=NetCDF::endef($nc_id);
die "Could not end definition mode\n" if $rcd < 0;
print STDOUT "ok\n";

# Output variables whose values are already known
$rcd=NetCDF::varput($nc_id,$var_id{'mlc_id'},\@{$var_srt{'mlc_id'}},\@{$var_cnt{'mlc_id'}},\@mlc);
$rcd=NetCDF::varput($nc_id,$var_id{'iso_id'},\@{$var_srt{'iso_id'}},\@{$var_cnt{'iso_id'}},\@iso);

# Assemble regular expression string
if($mlc_nbr == $mlc_nbr_max_htrn && $iso_nbr == $iso_nbr_max_htrn){
    $rx_sng='.'; # Each line is due to a valid isotopomer
}else{ # Get only specified molecules and isotopes of specified molecules
    $rx_sng='^';
    if($mlc_nbr > 1){$rx_sng.='(';}
    foreach $mlc_id (@mlc){
	$mlc_nm=$mlc_sng{$mlc_id};
	@iso_crr_htrn=@{$mlc_iso{$mlc_nm}}; # Put braces around list to extract it from an HoL
	$iso_nbr_crr_htrn=$#{$mlc_iso{$mlc_nm}}+1; # Put braces around list to extract it from an HoL
	$iso_nbr_crr_usr=0;
	@iso_crr_usr=(); # This array will hold the isotopic indices relative to the current molecule
	for($idx=0;$idx<$iso_nbr_crr_htrn;$idx++){
	    for($iso_idx=0;$iso_idx<$iso_nbr;$iso_idx++){
		if($iso_crr_htrn[$idx] == $iso[$iso_idx]){push @iso_crr_usr,$idx+1;}
	    } # end loop over iso_idx
	} # end loop over idx
	$iso_nbr_crr_usr=$#iso_crr_usr+1;
# Only specify isotopes in search string if user did not request all isotopes for a given molecule
	if($mlc_nbr > 1){$rx_sng.='(';}
	if($iso_nbr_crr_htrn == $iso_nbr_crr_usr){
	    $rx_sng.=sprintf('%2d',$mlc_id);
	}else{
	    if($iso_nbr_crr_usr == 1){
		$rx_sng.=sprintf('%2d%1d',$mlc_id,$iso_crr_usr[0]);
	    }else{
		$rx_sng.=sprintf('%2d(',$mlc_id);
		for($idx=0;$idx<$iso_nbr_crr_usr;$idx++){
		    $rx_sng.=sprintf('%1d',$iso_crr_usr[$idx]);
		    $rx_sng.='|' if $idx+1 != $iso_nbr_crr_usr;
		} # endfor        
		$rx_sng.=')';
	    } # endelse
	} # end else
	if($mlc_nbr > 1){$rx_sng.=')';}
	$rx_sng.='|' if $mlc_id != $mlc[$#mlc];
    } # end loop over $mlc
    if($mlc_nbr > 1){$rx_sng.=')';}
} # end else
printf STDOUT ("rx_sng = %s\n",$rx_sng);

# Open ASCII input data file
open FL_IN,$fl_in or die "$prg_nm: ERROR unable to open $fl_in: $!\n"; # $! is the system error sng

# Following is a one-based ruler, a zero-based ruler, and a "typical" HITRAN line
# Perl indexing of characters and strings is zero-based
#1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
#0--------10--------20--------30--------40--------50--------60--------70--------80--------90---------
#0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
# 71     .000001 1.850E-49 0.000E+00.0320.0000 4383.2507 .50 .000000  1  1         Q55Q55   000 0 0 0

$ln_cnt=0; # [nbr] Number of records written
while(<FL_IN>){
    if (/$rx_sng/){
# Parse input line
	($ln_mlc,$ln_iso,$ln_ctr,$ln_str,$HWHM_air,$HWHM_slf,$ln_nrg_lwr,$HWHM_tpt_dpn_xpn,$ln_sft_air,$qnt_glb_upr,$qnt_glb_lwr,$qnt_lcl_upr,$qnt_lcl_lwr)=unpack '@0 a2 @2 a1 @3 a12 @15 a10 @35 a5 @40 a5 @46 a10 @55 a4 @59 a8 @67 a3 @70 a3 @73 a9 @82 a9',$_;
	if($ln_ctr >= $wvn_min && $ln_ctr <= $wvn_max){
# Write record, making sure record variables are listed in order of creation
	    $rcd=NetCDF::recput($nc_id,$ln_cnt,[\$HWHM_air,\$HWHM_slf,\$HWHM_tpt_dpn_xpn,\$ln_ctr,\$ln_iso,\$ln_mlc,\$ln_nrg_lwr,\$ln_sft_air,\$ln_str,\$qnt_glb_lwr,\$qnt_glb_upr,\$qnt_lcl_lwr,\$qnt_lcl_upr]);
	    die "Could not write record variables\n" if $rcd < 0;
	    $ln_cnt++; # [nbr] Number of records written
	} # end if line is in correct region
    } # end if isotopomer matches desired list
} # end while <fl_in>
close fl_in;
print STDOUT "Ingested $fl_in\n";
print STDOUT "Wrote $ln_cnt records\n";

# Close netCDF file
print STDOUT 'Closing netCDF file...';
$rcd=NetCDF::close($nc_id);
die "Could not close netCDF file\n" if $rcd < 0;
print STDOUT "ok\n";

if($dbg_lvl == 1){		 
    if($ln_cnt != 0){		 
	printf STDOUT ("ln_mlc=%s ln_iso=%s ln_ctr=%12.6f ln_str=%10.3e HWHM_air=%5.4f HWHM_slf=%5.4f ln_nrg_lwr=%10.4f HWHM_tpt_dpn_xpn=%4.2f ln_sft_air=%8.6f qnt_glb_lwr=%3d qnt_glb_upr=%3d qnt_lcl_lwr=%s qnt_lcl_upr=%s\n",$ln_mlc,$ln_iso,$ln_ctr,$ln_str,$HWHM_air,$HWHM_slf,$ln_nrg_lwr,$HWHM_tpt_dpn_xpn,$ln_sft_air,$qnt_glb_lwr,$qnt_glb_upr,$qnt_lcl_lwr,$qnt_lcl_upr);
    } # end if $ln_cnt
} # end if dbg 

# Get elapsed times
&time_end($srt_usr_tm,$srt_sys_tm,$srt_child_usr_tm,$srt_child_sys_tm);

# Exit gracefully
exit 0;

sub usg_prn($){
    my $exit=@_;

    print STDOUT <<EndOfUsage;
Usage:  htrn2nc.pl [options] fl_in fl_out

   --dbg integer            Sets the debugging level (1..9)
   --help                   Display this help message
   --iso integer            HITRAN isotope number (1..9)
   --mlc integer            HITRAN molecule number (1..37)
   --mlt                    Toggles multiple molecules
   --wvn_max float          Maximum wavenumber to extract
   --wvn_min float          Minimum wavenumber to extract
   fl_in                    Name of the HITRAN ASCII input file to process
   fl_out                   Name of the binary netCDF output file

and unambiguous option abbreviations are permitted.
EndOfUsage
    exit $exit if $exit != 0;
} # end sub usg_prn()

