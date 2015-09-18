#!/usr/bin/perl 

# Purpose: Backup files to NCAR MSS or to Zip or LS120 drive or to tape

# Usage: Run as root so that read permission is available on all files
# /home/zender/bin/sh/bck.pl --zip

# Backup of dust.ps.uci.edu to zip drive:
# sudo mount /zip
# cd /

# sudo tar --preserve-permissions -cvzf /zip/dust-home-20001028.tar.gz ./home/zender/*
# sudo tar --preserve-permissions -cvzf /zip/dust-etc-20001028.tar.gz ./etc/*
# sudo tar --preserve-permissions -cvzf /zip/dust-data-20001028.tar.gz ./data/zender/aca/* ./data/zender/dst/* ./data/zender/igbp/* ./data/zender/mie/* ./data/zender/nr/* ./data/zender/pix/* ./data/zender/specfun/* ./data/zender/tex/* 

# Backup of lanina.ps.uci.edu to ls120 drive:
# sudo mount /ls120
# cd /
# sudo tar --preserve-permissions -cvzf /ls120/lanina-home-20001028.tar.gz ./home/zender/*
# sudo tar --preserve-permissions -cvzf /ls120/lanina-etc-20001028.tar.gz ./etc/*
# sudo tar --preserve-permissions -cvzf /ls120/lanina-data-20001028.tar.gz ./data/zender/aca/* ./data/zender/dst/* ./data/zender/igbp/* ./data/zender/mie/* ./data/zender/nr/* ./data/zender/pix/* ./data/zender/specfun/* ./data/zender/tex/* 

# Debugging:

# Production:

# Normal backups of home and NCAR and UCI
# bck.pl

# Backup z.ppp.ucar.edu to tape (must be root)
# bck.pl --tape

# Backup odin.cgd.ucar.edu to MSS (must be root)
# bck.pl --tape

BEGIN{
    unshift @INC,$ENV{'HOME'}.'/perl'; # Location of csz.pl and DBG.pm HaS98 p. 170
} # end BEGIN

my $CVS_Header='$Id$';

# Specify modules
use strict; # Protect all namespaces
use Getopt::Long; # GNU-style getopt

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
&time_srt($lcl_date_time,$srt_usr_tm,$srt_sys_tm,$srt_child_usr_tm,$srt_child_sys_tm);
printf STDOUT ("Start user time %f\n",$srt_usr_tm);

# Declare local variables
$HELP=$False;
$RCS_Header='$Id$';
$RCS_Revision='$Revision$';
$tape=$False;
$Boolean=$False;
$zip=$False;
$MSS=$False;
$CRR_FL_IN_WAS_RCP=$False;
$NEED_TO_RCP_FL_OUT=$False;
$HOST=$ENV{'HOST'};
$HOSTNAME=$ENV{'HOSTNAME'};
$PVM_ARCH=$ENV{'PVM_ARCH'};
$dbg_lvl=0;
$float_foo=102.543;
$int_foo=73;

# Parse command line arguments: '!' means Boolean, '|' is OR, '=' specifies required argument: 'i' is integer, 'f' is float, 's' is string
$rcd=Getopt::Long::Configure('no_ignore_case'); # Turn on case-sensitivity
$rcd=GetOptions( # man Getopt::GetoptLong
		   'tape!' => \$tape,
		   'Boolean!' => \$Boolean,
		   'dbg_lvl=i' => \$dbg_lvl,
		   'float_foo=f' => \$float_foo, 
		   'help|usage', => \$HELP,
		   'foo=i' => \@foo,
		   'zip!' => \$zip,
		   'MSS!' => \$MSS,
		 ); # end GetOptions arguments

if($HELP){		 
#    &usg_prn();
    exit 0;
} # end HELP

# Do some OS-dependent housekeeping to speed things up and avoid possible ambiguities
if($PVM_ARCH =~ m/SUN/){ # See Camel p. 81 for =~ and m//
    $tar_cmd='gtar';
}elsif($PVM_ARCH eq 'LINUX'){		 
    $tar_cmd='tar';
} # endif PVM_ARCH

$fl_lst=''; $fl_bck='bck_'.$HOST.'_'.&YYYYMMDD_HHMM().'.tar.gz';
$home_fl=''; $fl_bck_home='bck_'.&YYYYMMDD_HHMM().'_'.$HOST.'_home.tar.gz';
$data_fl=''; $fl_bck_data='bck_'.&YYYYMMDD_HHMM().'_'.$HOST.'_data.tar.gz';
$nt_fl=''; $fl_bck_nt='bck_'.&YYYYMMDD_HHMM().'_'.$HOST.'_nt.tar.gz';
$linux_fl=''; $fl_bck_linux='bck_'.&YYYYMMDD_HHMM().'_'.$HOST.'_linux.tar.gz';
$usrlcl_fl=''; $fl_bck_usrlcl='bck_'.&YYYYMMDD_HHMM().'_'.$HOST.'_usrlcl.tar.gz';
$dell_fl=''; $fl_bck_dell='bck_'.&YYYYMMDD_HHMM().'_'.$HOST.'_dell.tar.gz';
$etc_fl=''; $fl_bck_etc='bck_'.&YYYYMMDD_HHMM().'_'.$HOST.'_etc.tar.gz';

if($HOST eq 'dust.ps.uci.edu'){		 
    $rootdot_fl="/root/.[a-zA-Z]*"; $fl_lst=$fl_lst.$rootdot_fl.' ';
} # endif LINUX
$dot_fl=$HOME."/.[a-zA-Z]*"; $fl_lst=$fl_lst.$dot_fl.' ';

$home_fl='/home';
$data_fl="/data/zender/aca/* /data/zender/dst/* /data/zender/igbp/* /data/zender/mie/* /data/zender/nr/* /data/zender/pix/* /data/zender/specfun/* /data/zender/tex/* /data/zender/aeroce/* /data/zender/data/* /data/zender/ps/*";
$usrlcl_fl='/usr/local';
$dell_fl='/dell';
$etc_fl='/etc';

if($MSS){ # MSS backup
    &cmd_prc(sprintf("/bin/rm -f /data/zender/tmp/%s",$fl_bck));
    &cmd_prc(sprintf("%s --absolute-names -cvzf /data/zender/tmp/%s %s",$tar_cmd,$fl_bck,$fl_lst));
    &cmd_prc(sprintf("ls -l /data/zender/tmp/%s",$fl_bck));
    &cmd_prc(sprintf("ssh dataproc.ucar.edu mkdir -p /tmp/zender"));
    &cmd_prc(sprintf("scp -p /data/zender/tmp/%s dataproc.ucar.edu:/tmp/zender",$fl_bck));
    &cmd_prc(sprintf("ssh dataproc.ucar.edu msrcp -period 365 -wpwd zender /tmp/zender/%s mss:/ZENDER/bck/%s",$fl_bck,$fl_bck));
    if($dbg_lvl > 0){&cmd_prc(sprintf("ssh dataproc.ucar.edu ls -lg /tmp/zender/%s",$fl_bck))};
    if($dbg_lvl > 0){&cmd_prc(sprintf("ssh dataproc.ucar.edu %s -tvzf /tmp/zender/%s",$tar_cmd,$fl_bck))};
    &cmd_prc(sprintf("ssh dataproc.ucar.edu /bin/rm -f /tmp/zender/%s",$fl_bck));
#&cmd_prc(sprintf("/bin/rm -f /data/zender/tmp/%s",$fl_bck));

}elsif($zip){ # Zip drive backup

    if($HOST eq 'dust.ps.uci.edu'){
	$linux_fl='/etc /usr/lib/texmf/texmf/tex/latex/csz /usr/local/info/dir /root /usr/local/src /usr/local/netscape'; $fl_lst=$fl_lst.$linux_fl.' ';
	$nt_fl='/nt/home/zender /nt/boot.ini /nt/autoexec.bat /nt/home/r /nt/bin /nt/gnuwin32 /nt/netcdf-3.3.1 /nt/emacs-19.34'; $fl_lst=$fl_lst.$nt_fl.' ';
    } # endif dust

    &cmd_prc(sprintf("/bin/rm -f /zip/bck*tar*gz"));
    &cmd_prc(sprintf("%s --absolute-names --preserve-permissions -cvzf /zip/%s %s",$tar_cmd,$fl_bck_home,$home_fl));
    &cmd_prc(sprintf("%s --absolute-names --preserve-permissions -cvzf /zip/%s %s",$tar_cmd,$fl_bck_data,$data_fl));
    &cmd_prc(sprintf("%s --absolute-names --preserve-permissions -cvzf /zip/%s %s",$tar_cmd,$fl_bck_etc,$etc_fl));
    &cmd_prc(sprintf("%s --absolute-names --preserve-permissions -cvzf /zip/%s %s",$tar_cmd,$fl_bck_usrlcl,$usrlcl_fl));
    &cmd_prc(sprintf("%s --absolute-names --preserve-permissions -cvzf /zip/%s %s",$tar_cmd,$fl_bck_dell,$dell_fl));

}elsif($tape){ # Tape backup

    if($HOST eq 'dust.ps.uci.edu'){		 
	$linux_fl='/etc /usr/lib/texmf/texmf/tex/latex/csz /usr/local/info/dir /root /usr/local/src /usr/local/netscape'; $fl_lst=$fl_lst.$linux_fl.' ';
	$nt_fl='/nt/home/zender /nt/boot.ini /nt/autoexec.bat /nt/home/r /nt/bin /nt/gnuwin32 /nt/netcdf-3.3.1 /nt/emacs-19.34'; $fl_lst=$fl_lst.$nt_fl.' ';
    } # endif dust
    &cmd_prc(sprintf("/bin/rm -f /usr/tmp/zender/bck*z*tar*gz"));
    &cmd_prc(sprintf("%s --absolute-names --preserve-permissions -cvzf /usr/tmp/zender/%s %s",$tar_cmd,$fl_bck,$fl_lst));

    &cmd_prc(sprintf("mt -f /dev/ftape rewind",));
    &cmd_prc(sprintf("%s --absolute-names --preserve-permissions -cvf /dev/ftape /usr/tmp/zender/%s",$tar_cmd,$fl_bck));

}else{
    die 'Backup type unknown';
} # endelse
