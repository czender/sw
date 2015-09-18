#!/contrib/bin/perl 

my $RCS_Header='$Id$';

# Purpose: Download and install egcs

# Usage:
# /home/zender/perl/egcs_nst.pl

use strict; # Protect all namespaces; this produces many errors with non-expert code
use File::Basename; # Parses filenames

# Specify modules
use Getopt::Long; # GNU-style getopt
require '/home/zender/perl/csz.pl'; # Personal library: date_time(), YYYYMMDD(), ...

# Set output flushing to help debugging on hard crashes. 
# These options update the filehandle after every output statement.
# See Camel book, p. 110.
select((select(STDOUT),$|=1)[0]);
select((select(STDERR),$|=1)[0]);

# Timing information
my $lcl_date_time=&time_srt();

# Declare local variables
my ($prg_nm,$prg_dsc,$prg_vrs,$prg_date);
my ($pth_in,$fl_sfx);

my ($mk_cmd,$tar_cmd,$rmt_mch);

# Set defaults 
my $False=0;
my $True=1;

my $RCS_Date='$Date$';
my $RCS_Id='$Id$';
my $RCS_Revision='$Revision$';
my $HOME=$ENV{'HOME'};
my $HOST=$ENV{'HOST'};
my $HOSTNAME=$ENV{'HOSTNAME'};
my $PVM_ARCH=$ENV{'PVM_ARCH'};
my $src_pth='/l9/zender/data/egcs-1.1b'; # Where the distribution source will be
my $obj_pth='/l9/zender/data/foo'; # Where the distribution will be built
my $nst_pth='/usr/local'; # Where the distribution will be installed
my $dst_fl='egcs-1.1.1.tar.gz'; # Name of gzipped tarfile distribution

if($PVM_ARCH =~ m/SUN/){ # See Camel p. 81 for =~ and m//
    $tar_cmd='gtar';
    $mk_cmd='gmake';
}elsif($PVM_ARCH =~ m/CRAY/){
    $tar_cmd='tar';
    $mk_cmd='gnumake';
}else{
    $tar_cmd='tar';
    $mk_cmd='make';
} # endelse
if($src_pth eq '/home/zender'){die "$prg_nm: ERROR \$src_pth eq $src_pth";} # This would be disasterous

$prg_dsc='egcs installer'; # Program description
($prg_nm,$prg_vrs)=$RCS_Id =~ /: (.+).pl,v ([\d.]+)/; # Program name and version
$prg_vrs.='*' if length('$Locker:  $ ') > 12; # Tack '*' if it is not checked in into RCS.
($prg_nm,$pth_in,$fl_sfx)=fileparse($0,''); # $0 is program name Camel p. 136.
($prg_date)=unpack '@7 a19',$RCS_Date;

&cmd_prc("mkdir -p $src_pth"); # Make distribution directory if necessary
&cmd_prc("mv $dst_fl $src_pth"); # Move tarfile
chdir $src_pth or die "$prg_nm: ERROR unable to chdir to $src_pth: $!\n"; # $! is the system error sng
&cmd_prc("$tar_cmd -xvzf $dst_fl"); # Unpack tarfile

&cmd_prc("mkdir -p $obj_pth"); # Make build directory if necessary (should not be beneath src or nst)
&cmd_prc("mkdir -p $nst_pth"); # Make install directory if necessary
chdir $obj_pth or die "$prg_nm: ERROR unable to chdir to $obj_pth: $!\n"; # $! is the system error sng
&cmd_prc("$src_pth/configure --prefix=$nst_pth");
&cmd_prc("make bootstrap");
&cmd_prc("make install");
# NB: egcs seems to build with -O2 -g both on by default. Should I get rid of -g?

# Currently egcs seems to fail building something called ginfo:
# /usr/foo/texinfo/info: make
# terminal.o: In function `terminal_begin_using_terminal':
# /usr/egcs-1.1.1/texinfo/info/terminal.c:139: undefined reference to `tputs', `tgoto', `tgetstr' ...
