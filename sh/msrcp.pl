#!/usr/bin/perl 
# -*-Perl-*-

# Purpose: Emulate msrcp, msread, mswrite, and rdsp call on non-NCAR machines
#          Copy files between NCAR MSS and UCI via intermediary (e.g., dataproc.ucar.edu)

# Usage: Program determines which behavior to emulate (msrcp, msread, mswrite, or rdsp)
# based on command it is invoked with. Thus Perl script msrcp.pl should be
# symbolically linked to msrcp, msread, mswrite, and rdsp.

# Features:
# 1. msread (and msrcp in read mode) use cosconvert on intermediate machine to remove COS-blocking from files when "-f BI" option is used
# 2. Configured to redirect msread and mswrite commands to local disk for particular hosts to facilitate running NCAR models at remote locations
# 3. Automatically attempts retrieval from NCAR MSS if search of local disk fails
# 4. Able to override searches of local disk in favor of always using NCAR MSS

# Known Problems:
# 1. Code handles all cases where source and destination are files
# 2. Does not always work when source and/or destination are directories
# 3. Does not pass through all valid msrcp, msread, and mswrite switches
#    A. Some msread and mswrite switches are obsolete (do not translate to msrcp),
#    B. Some valid msrcp switches are simply not implemented (e.g., -maxtrans).
# 4. Unable to verify asynchronous writes were successful
# 5. Options have not been thoroughly checked to avoid ambiguities between msread, mswrite, and msrcp (e.g., -R)
# 6. Unable to use `.' for destination (but this should be easy to fix)

# msrcp mss:/DSS/A03797 /ptmp/zender/NCEP/A03797 
# Backgrounding jobs requires dummy tty input from /dev/null
# msrcp mss:/ZENDER/match/dstmch29/hist/h0001 /data/zender/dstmch29/h0001 < /dev/null &
# This /dev/null redirection is required for all three front ends: msrcp, msread, and mswrite

# Testing:
# msrcp --dbg=2 mss:/ZENDER/tmp/in.nc in.nc 
# msrcp --dbg=2 mss:/ZENDER/tmp/in.nc ./in.nc 
# msrcp --dbg=2 ~/nco/data/in.nc mss:/ZENDER/tmp/in.nc
# msrcp --dbg=2 --async ~/nco/data/in.nc mss:/ZENDER/tmp/in.nc

# Redirect calls to msread and mswrite to local disk 
# mswrite --dbg=2 -t 365 -v ~/nco/data/in.nc /ZENDER/tmp/in.nc
# msread --dbg=2 in.nc /ZENDER/tmp/in.nc

# Verify success on all three forms of mswrite command issued by MATCH
# mswrite -t 1000 -v -c "DAYS: 0.027-0.111 DATES: 0.666Z 930501 - 2.666Z 930501" h0001 /ZENDER/match/dstmch04/hist/h0001
# mswrite -t  365 -v rstrt /ZENDER/match/dstmch04/rest/rstrt
# mswrite -t  365 -v -c "DATE: 2.666Z 930501" r0001 /ZENDER/match/dstmch04/rest/r0001

# Verify success on all three forms of rdsp command issued by MATCH
# rdsp -p 1000 -d -w passwd -c "DAYS: 0.027-0.027 DATES: 0.666Z 19980101 - 0.666Z 19980101" h0001 /ZENDER/match/dstmch73/hist/h0001 &
# rdsp -p  365 rstrt /ZENDER/match/dstmch73/rest/rstrt &
# rdsp -p  365 -d -c "DATE: 0.666Z 19980101" r0001 /ZENDER/match/dstmch73/rest/r0001 &

# Verify success on all forms of msread command issued by MATCH
# Reading data tapes:
# msread -f BI -R /ptmp/zender/NCEP/A03797 /DSS/A03797
# Reading restart tapes:
# msread r0010 /ZENDER/match/dstmch32/rest/r0010
# msread in.nc /ZENDER/tmp/in.nc

# mswrite commands issued by CCM3/LSM:
# 19991228: Unlike MATCH, CCM and LSM background shell mswrite commands by default
# mswrite -c "time constant history file" -t  365 ./lsmh_timcon.nc /ZENDER/csm/dstccm80/lnd/hist/lsmh_timcon.nc
# mswrite -c "initial lsm data" -t  365 ./lsmi_00000601_00000.nc /ZENDER/csm/dstccm80/lnd/hist/lsmi_00000601_00000.nc; /bin/rm ./lsmi_00000601_00000.nc
# mswrite -c "MONTH MEAN DAYS: 000001.000-000031.000 DATES: 0.000Z 00000502 - 0.000Z 00000601" -t  365 ./lsmh_0000-05.nc /ZENDER/csm/dstccm80/lnd/hist/lsmh_0000-05.nc
# mswrite  -t  365 ./lsmr_00000601_00000 /ZENDER/csm/dstccm80/lnd/rest/lsmr_00000601_00000
# mswrite -c "DAYS: 00000000.000-00000000.042 DATES: 0.000Z 00000101 - 1.000Z 00000101" -t  365 -w passwd h0001 /ZENDER/csm/dstccm80/ccm3/hist/h0001
# mswrite -c "MONTHLY-AVERAGE FILE - DATES:  0.000Z 00000502 - 0.000Z 00000601" -t  365 -w passwd 0000-05 ZENDER/csm/dstccm80/ccm3/hist/0000-05
# mswrite -c "DAYS: 00000031.000-00000031.000 DATES: 0.000Z 00000601 - 0.000Z 00000601" -t  365 -w passwd ccmi_00000601_00000 /ZENDER/csm/dstccm80/ccm3/hist/ccmi_00000601_00000
# mswrite  -t  365 -w passwd ccmr_00000601_00000 /ZENDER/csm/dstccm80/ccm3/rest/ccmr_00000601_00000
# mswrite  -t  365 -w passwd ccmr_00000601_00000.A /ZENDER/csm/dstccm80/ccm3/rest/ccmr_00000601_00000.A

BEGIN{unshift @INC,$ENV{'HOME'}.'/perl'} # Location of DBG.pm HaS98 p. 170

my $CVS_Header='$Id$';

# Standard modules
use strict; # Protect all namespaces
use File::Basename; # Parses filenames
use Getopt::Long; # GNU-style getopt

# 3rd party modules

# Personal modules
use DBG; # Debugging constants
require 'csz.pl'; # Personal library: date_time(), YYYYMMDD(), ...

# Set output flushing to help debugging on hard crashes 
# These options update filehandle after every output statement
# See Camel book, p. 110
select((select(STDOUT),$|=1)[0]);
select((select(STDERR),$|=1)[0]);

# Timing information
my ($lcl_date_time,$srt_usr_tm,$srt_sys_tm,$srt_child_usr_tm,$srt_child_sys_tm);
($lcl_date_time,$srt_usr_tm,$srt_sys_tm,$srt_child_usr_tm,$srt_child_sys_tm)=time_srt();
#printf STDOUT ("Start time: %s\n",$lcl_date_tm);

# Declare local variables
my ($prg_dsc,$prg_vrs,$prg_date);
my ($cmd_ln,$prg_nm,$prg_pth,$prg_sfx);
my ($idx,$rcd);
my ($ssh_cmd,$scp_cmd,$scp_flg,$cp_cmd,$mv_cmd,$rm_cmd,$ccm2nc_cmd,$usr_nm);
my ($drc_lcl,$drc_rmt,$drc_mss);
my ($READ,$WRITE,$opt_cmd);
my ($fl_src,$fl_dst);
my ($fl_lcl,$fl_lcl_nm,$fl_lcl_pth,$fl_lcl_sfx);
my ($fl_rmt,$fl_rmt_nm,$fl_rmt_pth,$fl_rmt_sfx);
my ($drc_lcl_pfx); # Local disk replacement prefix for mswrite
my ($fl_rmt_lcl); # Look for local copy here before trying MSS

# Set defaults 
my $False=0;
my $True=1;

my $CVS_Date='$Date$';
my $CVS_Id='$Id$';
my $CVS_Revision='$Revision$';
my $READ=$False; # [flg] Operation is read
my $WRITE=$False; # [flg] Operation is write
my $ccm2nc_cmd='ccm2nc'; # [cmd] Command that behaves like ccm2nc (add -xz option here if desired)
my $cp_cmd='cp -p -f'; # [cmd] Command that behaves like cp
my $mv_cmd='mv -f'; # [cmd] Command that behaves like mv
my $data=$ENV{'DATA'};
my $fl_is_CCM_history_tape=$False; # [flg] CCM history tapes may be treated specially
my $home=$ENV{'HOME'}; # [sng] Home directory
my $host=$ENV{'HOST'}; # [sng] Host name
my $mkdir_cmd='mkdir -p'; # [cmd] Command that behaves like mkdir
my $mss_pfx='mss:'; # [sng] Prefix which identifies file on NCAR MSS
my $rm_cmd='rm -f'; # [cmd] Command that behaves like rm
my $scp_cmd='scp -p'; # [cmd] Command that behaves like scp
my $ssh_cmd='ssh'; # [cmd] Command that behaves like ssh
my $uci_dmn_sng='uci.edu'; # [sng] UCI domain string
my $uci_flg=$False; # [flg] Machine is at UCI
my $usr_nm=$ENV{'USER'}; # [sng] User name

# Set defaults for command line arguments
my $async=$False; # Perform writes asynchronously
my $dbg_lvl=0; # Debugging level
my $cos_cnv=$False; # [flg] cosconvert file (for msrcp.pl, same as msread -f BI)
my $msrcp=$False; # [flg] Emulate msrcp
my $msrcp_cmd='msrcp'; # Command that behaves like msrcp
my $msread=$False; # [flg] Emulate msread
my $msread_R=$False; # [flg] Tells msread to overwrite local file Obsolete
my $msread_fmt_arg=''; # Option f Format argument for msread: "BI" means COS-blocked binary
my $mswrite=$False; # [flg] Emulate mswrite or rdsp
my $mswrite_cmt_arg=''; # Comment argument for mswrite Obsolete
my $mswrite_verbose=$False; # [flg] Verbose Obsolete
my $mss_lcl=$True; # [flg] msread and mswrite will look to local disk first
my $pwd_read=''; # Read password
my $pwd_write=''; # Write password
my $rm_rmt_mch_fl=$True; # [flg] Remove file on remote machine after use (saves space)
my $rm_lcl_mch_fl=$False; # [flg] Remove local file after rdsp
#my $rmt_mch_nm='dataproc.ucar.edu'; # Remote machine name
my $rmt_mch_nm='goldhill.cgd.ucar.edu'; # Remote machine name
my $rtn_prd=365; # Retention period in days
my $scp_flg=$False; # [flg] Use scp when required
my $use_ccm2nc=$True; # [flg] Create and store netCDF version of CCM history tapes

# Derived fields
#my $drc_rmt_pfx='/ptmp/'.$usr_nm; # Remote directory (i.e., intermediate directory on remote server)
my $drc_rmt_pfx='/tmp/'.$usr_nm; # Remote directory (i.e., intermediate directory on remote server)
if($rm_cmd =~ m/( -r)|( -R)|( --recursive)/){die "$prg_nm: ERROR Dangerous setting \$rm_cmd eq $rm_cmd";} # This would be disastrous

$prg_dsc='Clone of NCAR msrcp, msread, mswrite, and rdsp protocols'; # Program description
($prg_nm,$prg_vrs)=$CVS_Id =~ /: (.+).pl,v ([\d.]+)/; # Program name and version
$prg_vrs.='*' if length('$Locker:  $ ') > 12; # Tack on '*' if it is not checked in into CVS
($prg_nm,$prg_pth,$prg_sfx)=fileparse($0,''); # $0 is program name Camel p. 136
if(length($CVS_Date) > 6){($prg_date)=unpack '@7 a19',$CVS_Date;}else{$prg_date='Unknown';}

if($prg_nm =~ m/^msrcp/){$msrcp=$True;}
elsif($prg_nm =~ m/^msread/){$msread=$True;}
elsif($prg_nm =~ m/^mswrite/){$mswrite=$True;}
elsif($prg_nm =~ m/^rdsp/){$mswrite=$True;} # Treat rdsp same as mswrite for now
else{die "$prg_nm: ERROR executable name does not contain msrcp, msread, mswrite, or rdsp\n";}

# Save command line before digesting options
$cmd_ln=$prg_nm.' '; # [sng] Command line
for $idx (0..$#ARGV){
    $cmd_ln.=$ARGV[$idx].' ';
} # end loop over arg

# Parse command line arguments: '!' means Boolean (default asserts truth, 'no' prefix asserts false), '|' is OR, '=' specifies required argument: 'i' is integer, 'f' is float, 's' is string
$Getopt::Long::ignorecase=0; # Turn on case-sensitivity (old method)
#$rcd=Getopt::Long::Configure('no_ignore_case'); # Turn on case-sensitivity  (new method)
$rcd=GetOptions( # man Getopt::GetoptLong
		 'R!' => \$msread_R, # Tells msread to overwrite local file Obsolete
		 'async!' => \$async, # Perform writes asynchronously
		 'comment|c=s' => \$mswrite_cmt_arg, # Comment argument for mswrite Obsolete
		 'cos_cnv!' => \$cos_cnv, # cosconvert file (for msrcp.pl, same as msread -f BI)
		 'dbg_lvl|D=i' => \$dbg_lvl, # Debugging level
		 'f=s' => \$msread_fmt_arg, # Option f Format argument for msread: "BI" means COS-blocked binary
		 'mss_lcl!' => \$mss_lcl, # msread and mswrite will look to local disk first
		 'period|t|p=i' => \$rtn_prd, # Retention period in days
		 'rm_rmt_mch_fl!' => \$rm_rmt_mch_fl, # Remove file on remote machine after use (saves space)
		 'rm_lcl_mch_fl|d!' => \$rm_lcl_mch_fl, # Remove local file after rdsp
		 'rmt_mch_nm=s' => \$rmt_mch_nm, # Remote machine name
		 'rpwd|r=s' => \$pwd_read, # Read password
		 'scp_flg!' => \$scp_flg, # scp flag
		 'use_ccm2nc!' => \$use_ccm2nc, # Create and store netCDF version of CCM history tapes
		 'verbose!' => \$mswrite_verbose, # Verbose Obsolete
		 'wpwd|w=s' => \$pwd_write, # Write password
		); # end GetOptions arguments

if($dbg_lvl >= 2){
    print ("prg_nm = $prg_nm\n");
    print ("cmd_ln = $cmd_ln\n");
    print ("Positional arguments remaining after GetOptions() parsing:\n");
    for $idx (0..$#ARGV){
	print ("ARGV[$idx] = $ARGV[$idx]\n");
    } # end loop over arg
} # endif debug

# Definitions that depend on command line input
if($msrcp && $msread_R){die "$prg_nm: ERROR Script does not support -R = Recursive switch\n";}

# Parse positional arguments
if($#ARGV+1 < 2){die "$prg_nm: ERROR Called with $#ARGV+1 positional arguments, need at least 2\n";}
else{
# Source file name is penultimate positional argument
# Destination file name is final positional argument
    if($msrcp || $mswrite){
	$fl_src=$ARGV[$#ARGV-1];
	$fl_dst=$ARGV[$#ARGV];
    }elsif($msread){
	$fl_dst=$ARGV[$#ARGV-1];
	$fl_src=$ARGV[$#ARGV];
    } # endif    
# Decide which file is local and which is remote
# Add MSS prefix to remote file if calling from msread or mswrite
    if($msread || $fl_src =~ m/$mss_pfx/){
	if($msread){$fl_src=$mss_pfx.$fl_src;}
	$fl_rmt=$fl_src;
	$fl_lcl=$fl_dst;
	$READ=$True;
    } # endif
    if($mswrite || $fl_dst =~ m/$mss_pfx/){
	if($mswrite){$fl_dst=$mss_pfx.$fl_dst;}
	$fl_rmt=$fl_dst;
	$fl_lcl=$fl_src;
	$WRITE=$True;
    } # endif
    if($READ && $WRITE){die "$prg_nm: ERROR both READ and WRITE are True\n";}
    if(!$READ && !$WRITE){die "$prg_nm: ERROR both READ and WRITE are False\n$prg_nm: HINT Make sure $mss_pfx prefixes either source or destination file\n";}
} # endelse

if($#ARGV+1 > 2){
    for $idx (1..$ARGV[$#ARGV-2]){
	$opt_cmd.=" $ARGV[$idx];"
    } # end loop over idx
} # endif

# Decide if MSS is local or remote
if($host =~ m/$uci_dmn_sng/){$uci_flg=$True;}
# Machines not at UCI should always write to real MSS 
if(!$uci_flg && $WRITE){$mss_lcl=$False;}
# Machines at UCI may look locally first if invoked with mswrite, msread, or rdsp
if($uci_flg && $mss_lcl && ($mswrite || $msread)){
    if($data eq ''){$drc_lcl_pfx=$home;}else{$drc_lcl_pfx=$data;}
} # endif $uci_flg

# Decide if writing CCM history tape
if($fl_lcl =~ m/^h[0-9]{4}/){$fl_is_CCM_history_tape=$True;} # Right name format
if($fl_lcl =~ m/.nc$/){$fl_is_CCM_history_tape=$False;} # No, already converted

# Set optional commands for WRITE mode 
if($WRITE){
# Retention period
$opt_cmd.=" -period $rtn_prd";
# -async turns on asynchronous mode
if($async){$opt_cmd.=' -async';}
if($pwd_write){$opt_cmd.=" -wpwd $pwd_write";}
#if($mswrite_cmt_arg){$opt_cmd.=" -c $mswrite_cmt_arg";} # Obsolete
#if($mswrite_verbose){$opt_cmd.=" -v";} # Obsolete
} # endif WRITE

# Set optional commands for READ mode 
if($READ){
if($pwd_read){$opt_cmd.=" -rpwd $pwd_read";}
#if($msread_R)=''; #  # Obsolete
} # endif READ

# Print initialization state
if($dbg_lvl != $dbg_off){print ("$prg_nm: $prg_dsc, version $prg_vrs of $prg_date\n");}
if($dbg_lvl >= 2){print ("$prg_nm: Called with $#ARGV+1 positional arguments\n");}
if($dbg_lvl >= 2){print ("$prg_nm: \$#ARGV = $#ARGV\n");}
if($dbg_lvl >= 2){print ("$prg_nm: \$async = $async\n");}
if($dbg_lvl >= 2){print ("$prg_nm: \$dbg_lvl = $dbg_lvl\n");}
if($dbg_lvl >= 2){print ("$prg_nm: \$msread_R = $msread_R\n");}
if($dbg_lvl >= 2){print ("$prg_nm: \$msread_fmt_arg = $msread_fmt_arg\n");}
if($dbg_lvl >= 2){print ("$prg_nm: \$mswrite_cmt_arg = $mswrite_cmt_arg\n");}
if($dbg_lvl >= 2){print ("$prg_nm: \$pwd_read = $pwd_read\n");}
if($dbg_lvl >= 2){print ("$prg_nm: \$pwd_write = $pwd_write\n");}
if($dbg_lvl >= 2){print ("$prg_nm: \$rm_lcl_mch_fl = $rm_lcl_mch_fl\n");}
if($dbg_lvl >= 2){print ("$prg_nm: \$rm_rmt_mch_fl = $rm_rmt_mch_fl\n");}
if($dbg_lvl >= 2){print ("$prg_nm: \$rtn_prd = $rtn_prd\n");}
if($dbg_lvl >= 2){print ("$prg_nm: \$use_ccm2nc = $use_ccm2nc\n");}
if($dbg_lvl >= 2){print ("$prg_nm: \$msrcp = $msrcp\n");}
if($dbg_lvl >= 2){print ("$prg_nm: \$msread = $msread\n");}
if($dbg_lvl >= 2){print ("$prg_nm: \$mswrite = $mswrite\n");}
if($dbg_lvl >= 2){print ("$prg_nm: \$READ = $READ\n");}
if($dbg_lvl >= 2){print ("$prg_nm: \$WRITE = $WRITE\n");}
if($dbg_lvl >= 2){print ("$prg_nm: \$fl_is_CCM_history_tape = $fl_is_CCM_history_tape\n");}
if($dbg_lvl >= 2){print ("$prg_nm: \$usr_nm = $usr_nm\n");}
if($dbg_lvl >= 2){print ("$prg_nm: \$host = $host\n");}
if($dbg_lvl >= 2){print ("$prg_nm: \$uci_flg = $uci_flg\n");}
if($dbg_lvl >= 2){print ("$prg_nm: \$ssh_cmd = $ssh_cmd\n");}
if($dbg_lvl >= 2){print ("$prg_nm: \$scp_cmd = $scp_cmd\n");}
if($dbg_lvl >= 2){print ("$prg_nm: \$scp_flg = $scp_flg\n");}
if($dbg_lvl >= 2){print ("$prg_nm: \$mss_pfx = $mss_pfx\n");}
if($dbg_lvl >= 2){print ("$prg_nm: \$rmt_mch_nm = $rmt_mch_nm\n");}
if($dbg_lvl >= 2){print ("$prg_nm: \$drc_rmt_pfx = $drc_rmt_pfx\n");}
if($dbg_lvl >= 2){print ("$prg_nm: \$fl_src = $fl_src\n");}
if($dbg_lvl >= 2){print ("$prg_nm: \$fl_dst = $fl_dst\n");}
if($dbg_lvl >= 2){print ("$prg_nm: \$fl_lcl = $fl_lcl\n");}
if($dbg_lvl >= 2){print ("$prg_nm: \$fl_rmt = $fl_rmt\n");}
if($dbg_lvl >= 2){print ("$prg_nm: \$opt_cmd = $opt_cmd\n");}
if($dbg_lvl >= 2){print ("$prg_nm: \$mss_lcl = $mss_lcl\n");}
if($dbg_lvl >= 2){print ("$prg_nm: \$drc_lcl_pfx = $drc_lcl_pfx\n");}

# Main Code
# Temporary remote storage directory is $drc_rmt_pfx plus the MSS path
($fl_lcl_nm,$fl_lcl_pth,$fl_lcl_sfx)=fileparse($fl_lcl,'');
($fl_rmt_nm,$fl_rmt_pth,$fl_rmt_sfx)=fileparse($fl_rmt,'');
$drc_mss=$fl_rmt_pth; # $fl_rmt_pth has a trailing /
$drc_mss =~ s#/$##g; # Remove trailing /
$drc_mss =~ s#^$mss_pfx##g; # Remove leading mss:
$drc_rmt=$drc_rmt_pfx.$drc_mss; # Append MSS directory to temporary root directory
$drc_lcl=$drc_lcl_pfx.$drc_mss; # Append MSS directory to local directory
$fl_rmt_lcl=$drc_lcl_pfx.$drc_mss.'/'.$fl_rmt_nm; # Look for local copy here before trying MSS

if($dbg_lvl >= 2){print ("$prg_nm: \$fl_lcl = $fl_lcl\n");}
if($dbg_lvl >= 2){print ("$prg_nm: \$fl_lcl_nm = $fl_lcl_nm\n");}
if($dbg_lvl >= 2){print ("$prg_nm: \$fl_lcl_pth = $fl_lcl_pth\n");}
if($dbg_lvl >= 2){print ("$prg_nm: \$fl_lcl_sfx = $fl_lcl_sfx\n");}
if($dbg_lvl >= 2){print ("$prg_nm: \$fl_rmt = $fl_rmt\n");}
if($dbg_lvl >= 2){print ("$prg_nm: \$fl_rmt_nm = $fl_rmt_nm\n");}
if($dbg_lvl >= 2){print ("$prg_nm: \$fl_rmt_pth = $fl_rmt_pth\n");}
if($dbg_lvl >= 2){print ("$prg_nm: \$fl_rmt_sfx = $fl_rmt_sfx\n");}
if($dbg_lvl >= 2){print ("$prg_nm: \$fl_rmt_lcl = $fl_rmt_lcl\n");}
if($dbg_lvl >= 2){print ("$prg_nm: \$drc_mss = $drc_mss\n");}
if($dbg_lvl >= 2){print ("$prg_nm: \$drc_rmt_pfx = $drc_rmt_pfx\n");}
if($dbg_lvl >= 2){print ("$prg_nm: \$drc_rmt = $drc_rmt\n");}
if($dbg_lvl >= 2){print ("$prg_nm: \$drc_lcl = $drc_lcl\n");}

if($READ){

    if($mss_lcl){

	if($dbg_lvl >= 0){print ("$prg_nm: Attempting to find $fl_rmt on local disk in $drc_lcl/$fl_rmt_nm\n");} # endif dbg

# Copy file from local disk to specified directory
	$rcd=1; 
	if(-f $fl_rmt_lcl){ # WCS96 p. 86
	    if(-r $fl_rmt_lcl){ # WCS96 p. 86 Perhaps should use "_" ?
		$rcd=cmd_prc("$cp_cmd $drc_lcl_pfx$drc_mss/$fl_rmt_nm $fl_lcl");
		if(-f $fl_lcl && -r $fl_lcl){$rcd=0;}else{print ("$prg_nm: Unable to copy $fl_rmt_lcl\n");}
	    }else{print ("$prg_nm: Unable to read $fl_rmt_lcl\n");} # endelse
	}else{print ("$prg_nm: Unable to find $fl_rmt_lcl\n");} # endelse
	if($rcd != 0){print ("$prg_nm: Will attempt to read from NCAR MSS instead\n")} # endif

    } # endif mss_lcl

# Fall through to MSS if necessary
    if(!$mss_lcl || $rcd != 0){
# Read from MSS to local file
	if($dbg_lvl >= 0){print ("$prg_nm: Reading $fl_rmt from NCAR MSS to $fl_lcl via $rmt_mch_nm\n");} # endif dbg
	
# Make remote directory
	if($scp_flg){
	    cmd_prc("$ssh_cmd -l $usr_nm $rmt_mch_nm $mkdir_cmd $drc_rmt");
	    
# Copy file from MSS to remote disk
	    cmd_prc("$ssh_cmd -l $usr_nm $rmt_mch_nm $msrcp_cmd $fl_rmt $drc_rmt/$fl_lcl_nm");
	    
# Convert COS-blocked files to binary files
	    if($msread_fmt_arg eq 'BI' || $cos_cnv){
		cmd_prc("$ssh_cmd -l $usr_nm $rmt_mch_nm cosconvert -b $drc_rmt/$fl_lcl_nm");
	    } # endif COS-blocked
	} # endif scp_flg
	
# Copy file from remote disk to local disk
	if($dbg_lvl >= 0){print ("$prg_nm: WARNING Attempt to scp file from remote machine aborted. Finding the file on a local disk apparently failed. Since the NCAR security perimeter was installed in ~2003, attempts to read/write automatically to NCAR usually fail. Only workaround is to find file on the local disk.\n");} # endif dbg
	
	# Following two commands require password-less authentication
	if($scp_flg){
	    cmd_prc("$scp_cmd $usr_nm\@$rmt_mch_nm:$drc_rmt/$fl_lcl_nm $fl_lcl");
# Remove file from remote disk when transfer is complete
	    if($rm_rmt_mch_fl){cmd_prc("$ssh_cmd -l $usr_nm $rmt_mch_nm $rm_cmd $drc_rmt/$fl_lcl_nm");}
	} # endif scp enabled
	
    } # endif read from MSS
    
# endif $READ
}elsif($WRITE){

    if($mss_lcl){
# Make local directory
	if($dbg_lvl >= 0){print ("$prg_nm: Writing $fl_lcl to local disk $drc_lcl/$fl_lcl_nm\n");} # endif dbg
	-d $drc_lcl or cmd_prc("$mkdir_cmd $drc_lcl"); # -d test on Camel p. 20

# Redirect mswrite to store file locally rather than on MSS
# Do not bother moving/copying CCM history tapes which have been converted to netCDF
# Experience shows that the primary cause of errors with data archival is full disks
# Be conservative and do not archive redundant copies
	if($fl_is_CCM_history_tape && $use_ccm2nc){
# Convert local file to netCDF if it looks like a history tape
	    cmd_prc("$ccm2nc_cmd $fl_lcl $drc_lcl/$fl_lcl_nm.nc");
# Simply remove original since we have netCDF archive
	    if($rm_lcl_mch_fl){cmd_prc("$rm_cmd $fl_lcl")};
	}else{
# If instructed to delete, then use mv instead of sequential cp,rm
# This eliminates risk of using rm on $fl_lcl when $fl_lcl == $drc_lcl/$fl_lcl_nm
	    if($rm_lcl_mch_fl){cmd_prc("$mv_cmd $fl_lcl $drc_lcl/$fl_lcl_nm")}else{cmd_prc("$cp_cmd $fl_lcl $drc_lcl/$fl_lcl_nm")};
	} # end if 

    }else{ # endif mss_lcl
# Write local file to MSS
	if($dbg_lvl >= 0){print ("$prg_nm: Writing $fl_lcl to NCAR MSS file $fl_rmt via $rmt_mch_nm\n");} # endif dbg
	
# Make remote directory
	if($scp_flg){
	    cmd_prc("$ssh_cmd -l $usr_nm $rmt_mch_nm $mkdir_cmd $drc_rmt");
	    
# Copy local file to remote disk
	    cmd_prc("$scp_cmd $fl_lcl $usr_nm\@$rmt_mch_nm:$drc_rmt/$fl_lcl_nm");
	    
# Convert remote file to netCDF if it looks like a history tape
	    if($fl_is_CCM_history_tape && $use_ccm2nc){
		cmd_prc("$ssh_cmd -l $usr_nm $rmt_mch_nm $ccm2nc_cmd $drc_rmt/$fl_lcl_nm $drc_rmt/$fl_lcl_nm.nc");
		cmd_prc("$ssh_cmd -l $usr_nm $rmt_mch_nm $msrcp_cmd $opt_cmd $drc_rmt/$fl_lcl_nm.nc $fl_rmt.nc");
	    } # end if $fl_is_CCM_history_tape
	    
# Copy remote file to MSS
	    cmd_prc("$ssh_cmd -l $usr_nm $rmt_mch_nm $msrcp_cmd $opt_cmd $drc_rmt/$fl_lcl_nm $fl_rmt");
	    
	    if(!$async){
# Remove remote file -- never do this in asynchronous mode
		if($rm_rmt_mch_fl){cmd_prc("$ssh_cmd -l $usr_nm $rmt_mch_nm $rm_cmd $drc_rmt/$fl_lcl_nm");}
		
# Remove remote netCDF file -- never do this in asynchronous mode
		if($rm_rmt_mch_fl && $fl_is_CCM_history_tape && $use_ccm2nc){
		    cmd_prc("$ssh_cmd -l $usr_nm $rmt_mch_nm $rm_cmd $drc_rmt/$fl_lcl_nm.nc");
		} # end if $fl_is_CCM_history_tape
	    } # endif not async
	} # endif scp_flg
    } # endif write to MSS
} # endelse $WRITE

# Get elapsed times
#&time_end($srt_usr_tm,$srt_sys_tm,$srt_child_usr_tm,$srt_child_sys_tm);
