#!/contrib/bin/perl -w

# $Id$

# Purpose: Transfer files among NCAR, home, tape backup

# Usage:
# rcp -p odin.cgd.ucar.edu:/home/zender/bin/sh/rcp.pl z.ppp.ucar.edu:/home/zender/bin/sh/rcp.pl
# rcp -p /home/zender/bin/sh/rcp.pl odin.cgd.ucar.edu:/home/zender/bin/sh/rcp.pl 

$foo='
cd /home/zender/c++;cvs update
cd /home/zender/c;cvs update
cd /home/zender/crr;cvs update
cd /home/zender/dot;cvs update
cd /home/zender/dst;cvs update
cd /home/zender/f;cvs update
cd /home/zender/mie;cvs update
cd /home/zender/tex;cvs update
cd /fs/cgd/data0/zender/ccm_dst;cvs update
cd /fs/cgd/data0/zender/match_dst;cvs update
';

# Debugging:
# rcp.pl -CRM -dbg 2 -day 1

# Production:

# Specify modules
use Getopt::Long;
# Keep long options from being lowercased (bug before v. 2.4 of Getopt::Long)
$Getopt::Long::ignorecase=0;
use POSIX qw(ceil floor);
require '/home/zender/perl/csz.pl'; # Contains date_time(), cmd_prc(), round()

# Set output flushing to help debugging on hard crashes. 
# These options update filehandle after every output statement.
# See Camel book, p. 110.
select((select(STDOUT),$|=1)[0]);
select((select(STDERR),$|=1)[0]);

# Set defaults 
$False=0;
$True=1;

$HELP=$False;
$RCS_Header='$Id$';
$RCS_Revision='$Revision$';
$Boolean=$False;
$HOME=$ENV{'HOME'};
$HOST=$ENV{'HOST'};
$HOSTNAME=$ENV{'HOSTNAME'};
$PVM_ARCH=$ENV{'PVM_ARCH'};
$day=11;
$dbg_lvl=0;
$float_foo=102.543;
$int_foo=73;

# Parse command line arguments 
$result=GetOptions(
		   'Boolean!',
		   'dbg=i',
		   'float_foo=f',
		   'help|usage',
		   'int_foo=i',
		   'day=i',
		   );
if(defined($opt_Boolean)) {$Boolean=!$Boolean};
if(defined($opt_dbg)) {$dbg_lvl=$opt_dbg};
if(defined($opt_float_foo)) {$float_foo=$opt_float_foo};
if(defined($opt_help)) {$HELP=$opt_help};
if(defined($opt_int_foo)) {$int_foo=$opt_int_foo};
if(defined($opt_day)) {$day=$opt_day};

if($HELP){exit 0;} # endif HELP

# Do OS-dependent housekeeping to speed things up and avoid possible ambiguities
if($PVM_ARCH =~ m/SUN/){ # See Camel p. 81 for =~ and m//
    $host_lcl='odin.cgd.ucar.edu';
    $host_rmt='z.ppp.ucar.edu';
    $tar_cmd_lcl='gtar';
    $tar_cmd_rmt='tar';
}elsif($PVM_ARCH eq 'LINUX'){		 
    $host_lcl='z.ppp.ucar.edu';
    $host_rmt='odin.cgd.ucar.edu';
    $tar_cmd_lcl='tar';
    $tar_cmd_rmt='gtar';
} # endif PVM_ARCH

$fl_lst='';
$fl_rcp='rcp.tar.gz';

#$aca_fl='/home/zender/aca/*.F /home/zender/aca/*.pl /home/zender/aca/*.sh /home/zender/aca/Makefile'; $fl_lst=$fl_lst.$aca_fl.' ';
#$aer_fl='/home/zender/aer/*.dat /home/zender/aer/*.doc /home/zender/aer/*.cdl /home/zender/aer/*.nc'; $fl_lst=$fl_lst.$aer_fl.' ';
#$afgl_fl='/home/zender/afgl/*.txt'; $fl_lst=$fl_lst.$afgl_fl.' ';
#$anv_fl='/home/zender/anv/anv.tex /home/zender/anv/jgr_anv_02.tex /home/zender/anv/jgr_anv_03.tex '; $fl_lst=$fl_lst.$anv_fl.' ';
#$arese_fl='/home/zender/arese/arese.sh /home/zender/arese/arese.pl'; $fl_lst=$fl_lst.$arese_fl.' ';
#$avhrr_fl='/home/zender/avhrr/*.pl /home/zender/avhrr/*.pm'; $fl_lst=$fl_lst.$avhrr_fl.' ';
#$binsh_fl='/home/zender/bin/sh/rcp.pl /home/zender/bin/sh/bck.pl'; $fl_lst=$fl_lst.$binsh_fl.' ';
#$c_fl='/home/zender/c/*.c /home/zender/c/*.cc /home/zender/c/*.hh /home/zender/c/Makefile'; $fl_lst=$fl_lst.$c_fl.' ';
#$ck_fl='/home/zender/ck/*.F /home/zender/ck/*.c /home/zender/ck/Makefile'; $fl_lst=$fl_lst.$ck_fl.' ';
#$csz_fl='/home/zender/csz/*.F /home/zender/csz/*.c /home/zender/csz/*.h /home/zender/csz/Makefile'; $fl_lst=$fl_lst.$csz_fl.' ';
$dot_fl='/home/zender/.login /home/zender/.Xdefaults /home/zender/.xinitrc /home/zender/.rhosts /home/zender/.fvwm2rc /home/zender/.cshrc /home/zender/.mailrc'; $fl_lst=$fl_lst.$dot_fl.' ';
#$dmr_fl='/home/zender/dmr/*'; $fl_lst=$fl_lst.$dmr_fl.' ';
#$dmr_fl='/fs/cgd/home0/zender/dmr/dmr/*.F /fs/cgd/home0/zender/dmr/dmr/*.h'; $fl_lst=$fl_lst.$dmr_fl.' ';
#$dst_fl='/home/zender/dst/Makefile /home/zender/dst/*.cdl /home/zender/dst/*.F /home/zender/dst/*.h /home/zender/dst/*.f'; $fl_lst=$fl_lst.$dst_fl.' ';
#$dst_fl='/fs/cgd/home0/zender/dst/dst/*.F /fs/cgd/home0/zender/dst/dst/*.h /fs/cgd/home0/zender/dst/dst/*.com'; $fl_lst=$fl_lst.$dst_fl.' ';
#$elisp_fl='/home/zender/elisp/my_macros.el'; $fl_lst=$fl_lst.$elisp_fl.' ';
#$f_fl='/home/zender/f/*.F /home/zender/f/*.com /home/zender/f/*.h /home/zender/f/Makefile'; $fl_lst=$fl_lst.$f_fl.' ';
#$frc_fl='/home/zender/frc/*'; $fl_lst=$fl_lst.$frc_fl.' ';
#$job_fl='/home/zender/job/*.tex'; $fl_lst=$fl_lst.$job_fl.' ';
#$linux_fl='/home/zender/linux/*'; $fl_lst=$fl_lst.$linux_fl.' ';
$mk_fl='/home/zender/mk/Makefile*'; $fl_lst=$fl_lst.$mk_fl.' ';
#$mie_fl='/home/zender/mie/*.cc /home/zender/mie/*.hh /home/zender/mie/*.sh /home/zender/mie/Makefile'; $fl_lst=$fl_lst.$mie_fl.' ';
#$nc_fl='/home/zender/nc/nco.texi /home/zender/nc/Makefile /home/zender/nc/csz.c /home/zender/nc/ncap* /home/zender/nc/nc*.c /home/zender/nc/nc*.h /home/zender/nc/TODO /home/zender/nc/ChangeLog /home/zender/nc/in.cdl /home/zender/nc/*.F'; $fl_lst=$fl_lst.$nc_fl.' ';
#$nc_fl='/home/zender/nc/*'; $fl_lst=$fl_lst.$nc_fl.' ';
#$root_fl='/root/*'; $fl_lst=$fl_lst.$root_fl.' ';
#$rvw_fl='/home/zender/rvw/DSP97.tex /home/zender/rvw/jas_DSP97.tex'; $fl_lst=$fl_lst.$rvw_fl.' ';
#$ncl_fl='/home/zender/ncl/tst.ncl /home/zender/ncl/tddr_abs.ncl'; $fl_lst=$fl_lst.$ncl_fl.' ';
#$perl_fl='/home/zender/perl/csz.pl'; $fl_lst=$fl_lst.$perl_fl.' ';
#$ps_fl='/data2/zender/ps/arese_trn_alb.eps /data2/zender/ps/arese_9510[13][150]*.eps'; $fl_lst=$fl_lst.$ps_fl.' ';
#$slr_spc_fl='/home/zender/slr_spc/slr_spc* /home/zender/slr_spc/spc*.txt'; $fl_lst=$fl_lst.$slr_spc_fl.' ';
#$sage_fl='/home/zender/sage/*.i'; $fl_lst=$fl_lst.$sage_fl.' ';
#$tex_fl='/home/zender/tex/bib.bib /home/zender/tex/csz.sty'; $fl_lst=$fl_lst.$tex_fl.' ';
#$time_fl='/home/zender/time/Makefile /home/zender/time/slr_crd_jjm.F /home/zender/time/slr_crd_utl.F'; $fl_lst=$fl_lst.$time_fl.' ';
#$toms_fl='/home/zender/toms/*.pl /home/zender/toms/*.pm'; $fl_lst=$fl_lst.$toms_fl.' ';
#$www_fl='/home/zender/www/formtest*'; $fl_lst=$fl_lst.$www_fl.' ';

&cmd_prc(sprintf("/bin/rm -f /data/zender/tmp/%s",$fl_rcp));
&cmd_prc(sprintf("%s --absolute-names -cvzf /data/zender/tmp/%s %s",$tar_cmd_lcl,$fl_rcp,$fl_lst));
&cmd_prc(sprintf("ls -l /data/zender/tmp/%s",$fl_rcp));
&cmd_prc(sprintf("rcp -p /data/zender/tmp/%s %s:/data/zender/tmp/%s",$fl_rcp,$host_rmt,$fl_rcp));
&cmd_prc(sprintf("rsh %s %s --absolute-names -xvzf /data/zender/tmp/%s",$host_rmt,$tar_cmd_rmt,$fl_rcp));
&cmd_prc(sprintf("rsh %s /bin/rm -f /data/zender/tmp/%s",$host_rmt,$fl_rcp));
&cmd_prc(sprintf("/bin/rm -f /data/zender/tmp/%s",$fl_rcp));

