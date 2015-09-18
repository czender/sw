# $CVS_Header='$Id$';

# Purpose: Library of Perl subroutines

# Usage:
# do '/home/zender/perl/csz.pl'; # Personal library, cmd_prc(), date_time(), YYYYMMDD()...
# use lib qw(/home/zender/perl);
# require 'csz.pl'; # Personal library, cmd_prc(), date_time(), YYYYMMDD()...

sub numerically{$a <=> $b;} # Camel book, p. 218

sub round($){
# Purpose: Rounds to nearest integer. From the Perl FAQ.
# Usage: 
    my $flt=shift;
    return int($flt+0.5*($flt <=> 0));
} # end sub round()

sub yday(){
# Purpose: Returns Julian day of year
# Usage: 
    my($sec,$min,$hr,$mday,$mon,$yr,$wday,$yday,$isdst);
    ($sec,$min,$hr,$mday,$mon,$yr,$wday,$yday,$isdst)=localtime(time);
    return $yday;
} # end doy()

sub cmd_prc($){
# Purpose: Execute shell commands and pipe their standard output to stdout
# Usage: Newlines are automatically attached to incoming command before printing
# &cmd_prc("echo foo");
    my ($cmd)=@_;
    print STDOUT $cmd."\n";
    my $pid=open(CMD,$cmd.' |');
# It would be very nice to check command return status
# This does not seem possible with open(), which returns process PID
# Perhaps system() call is required, but can that return STDOUT?
#    if(defined($pid)){print "Successful open\n";}else{print "Unsuccessful open\n";} # Camel p. 193
    if(!$pid){print STDOUT 'Unable to open() on '.$cmd;}
    while(<CMD>){print STDOUT $_;}
    close CMD;
    return $pid;
} # end cmd_prc()

sub typeof{
# Purpose: Print type of a Perl variable
# Usage: printf "Type of \$var is ".&typeof($var)."\n";
    my $arg=shift(@_);
    my $typ=ref($arg);
    if(not $typ){$typ_sng='not a reference';}
    elsif($typ eq 'HASH'){$typ_sng='reference to HASH';}
    elsif($typ eq 'REF'){$typ_sng='reference to REF';}
    elsif($typ eq 'SCALAR'){$typ_sng='reference to SCALAR';}
    elsif($typ eq 'ARRAY'){$typ_sng='reference to ARRAY';}
    elsif($typ eq 'CODE'){$typ_sng='reference to CODE';}
    elsif($typ eq 'GLOB'){$typ_sng='reference to GLOB';}
    else{$typ_sng='reference to unknown type';}
    return $typ_sng;
} # end typeof()

sub time_srt(){
# Purpose: Return date and timing information for use by calling program
# Usage: 
# ($lcl_date_time,$srt_usr_tm,$srt_sys_tm,$srt_child_usr_tm,$srt_child_sys_tm)=time_srt();
# Values of arguments are all overwritten on output
# Copy passed arguments into named variables
    my ($lcl_date_tm,$srt_usr_tm,$srt_sys_tm,$srt_child_usr_tm,$srt_child_sys_tm);
    ($srt_usr_tm,$srt_sys_tm,$srt_child_usr_tm,$srt_child_sys_tm)=times;
    $lcl_date_tm=date_time();
    return ($lcl_date_tm,$srt_usr_tm,$srt_sys_tm,$srt_child_usr_tm,$srt_child_sys_tm);
} # end time_srt()

sub time_end{
# Purpose: Return date and timing information for use by calling program
# Usage: Call with arguments previously set using time_srt()
# time_end($srt_usr_tm,$srt_sys_tm,$srt_child_usr_tm,$srt_child_sys_tm);
    my ($srt_usr_tm,$srt_sys_tm,$srt_child_usr_tm,$srt_child_sys_tm);
    ($srt_usr_tm,$srt_sys_tm,$srt_child_usr_tm,$srt_child_sys_tm)=@_;
    my ($end_usr_tm,$end_sys_tm,$end_child_usr_tm,$end_child_sys_tm);
    ($end_usr_tm,$end_sys_tm,$end_child_usr_tm,$end_child_sys_tm)=times;
#    printf STDOUT ("Elapsed time: User CPU %.2f s, System CPU %.2f s\n",$end_usr_tm-$srt_usr_tm,$end_sys_tm-$srt_sys_tm);
    printf STDOUT ("Start user time %f\n",$srt_usr_tm);
    printf STDOUT ("End user time %f\n",$end_usr_tm);
    printf STDOUT ("Elapsed time: Child user CPU %.2f s, Child system CPU %.2f s\n",$end_child_usr_tm-$srt_child_usr_tm,$end_child_sys_tm-$srt_child_sys_tm);
    my $lcl_date_tm=date_time();
    printf STDOUT ("End time: %s\n",$lcl_date_tm);
} # end time_end()

sub mth2mm($){
# Purpose: Converts month name (or three letter abbreviation) to two digit month 
# "January" and "jan" return 01, "December" and "Dec" return 12
# Routine is inverse of mm2mmm()
# Usage: $mm=mth2mm("jan");
    my ($mth_abb)=@_;
    my $mm;
    if(!$mth_abb){die "mth2mm() reports invalid \$mth_abb = $mth_abb\n";}
    if($mth_abb =~ m/jan/i){$mm=1;} # i = Case-insensitive matching Camel p. 69
    elsif($mth_abb =~ m/feb/i){$mm=2;} # i = Case-insensitive matching Camel p. 69
    elsif($mth_abb =~ m/mar/i){$mm=3;} # i = Case-insensitive matching Camel p. 69
    elsif($mth_abb =~ m/apr/i){$mm=4;} # i = Case-insensitive matching Camel p. 69
    elsif($mth_abb =~ m/may/i){$mm=5;} # i = Case-insensitive matching Camel p. 69
    elsif($mth_abb =~ m/jun/i){$mm=6;} # i = Case-insensitive matching Camel p. 69
    elsif($mth_abb =~ m/jul/i){$mm=7;} # i = Case-insensitive matching Camel p. 69
    elsif($mth_abb =~ m/aug/i){$mm=8;} # i = Case-insensitive matching Camel p. 69
    elsif($mth_abb =~ m/sep/i){$mm=9;} # i = Case-insensitive matching Camel p. 69
    elsif($mth_abb =~ m/oct/i){$mm=10;} # i = Case-insensitive matching Camel p. 69
    elsif($mth_abb =~ m/nov/i){$mm=11;} # i = Case-insensitive matching Camel p. 69
    elsif($mth_abb =~ m/dec/i){$mm=12;} # i = Case-insensitive matching Camel p. 69
    else{die "ERROR Invalid month name $mth_abb in csz.pl\n";}
    return $mm;
} # end mth2mm()

sub mm2mmm($){
# Purpose: Converts month index [1..12] to three letter month abbreviation
# 01 returns "Jan", 12 returns "Dec"
# Routine is inverse of mth2mm()
# Usage: $mmm=mm2mmm($mm);
    my ($mm)=@_;
    if(!$mm){die "mm2mmm() reports invalid \$mm = $mm\n";}
    @mmm=('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec');
    if($mm < 1 || $mm > 12){die "ERROR mm2mmm() reports invalid month index \$mm = $mm\n";}
    return $mmm[$mm-1]; # -1 is for 0-based indices
} # end mm2mmm()

sub mm2mmm_lc($){
# Purpose: Converts month index [1..12] to three letter month abbreviation
# Lowercase months as in TOMS AOD files
# 01 returns "jan", 12 returns "dec"
# Routine is inverse of mth2mm()
# Usage: $mmm=mm2mmm($mm);
    my ($mm)=@_;
    if(!$mm){die "mm2mmm() reports invalid \$mm = $mm\n";}
    @mmm=('jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec');
    if($mm < 1 || $mm > 12){die "ERROR mm2mmm() reports invalid month index \$mm = $mm\n";}
    return $mmm[$mm-1]; # -1 is for 0-based indices
} # end mm2mmm()

sub yyyymmdd_prs($){
# Purpose: Parses date (in yyyymmdd format) into constituent year, month, and day
# Usage: ($yyyy,$mm,$dd)=yyyymmdd_prs($yyyymmdd);
    my ($yyyymmdd)=@_;
    if(length($yyyymmdd) != 8){
	printf STDOUT ("ERROR yyyymmdd_prs() reports length(\$yyyymmdd) = length($yyyymmdd) = %d\n",length($yyyymmdd));
	die;
    } # endif
    my ($yyyy,$mm,$dd);
    $dd=$yyyymmdd%100; 
    $mm=(($yyyymmdd%10000)-$dd)/100;
    $yyyy=int($yyyymmdd/10000);
    return ($yyyy,$mm,$dd);
} # end yyyymmdd_prs()

sub lst_add(\@$){
# Purpose: Sum values of list which may contain missing values
# Values in list equal to $mss_val are ignored
# List is passed by reference to avoid unnecessary copying
# Usage: lst_ttl=lst_add(@lst,$mss_val)
    my $mss_val; # [frc] Missing value
    my $lst; # [frc] Reference to array
    my $lst_nbr; # [nbr] Number of elements in list
    my $dbg_lvl=0; # [enm] Debugging level
    my $lst_idx=0; # [idx] Index for list
    my $lst_ttl=0.0; # [frc] Total of elements in list
# Last element in input is mss_val
    $mss_val=pop @_; # [frc] Missing value
# Syntax: Evaluate $lst in list context
    ($lst)=@_; # [frc] Reference to array
# Syntax: $lst is reference to @lst, and $# operates on it
    $lst_nbr=$#$lst; # [nbr] Number of elements in list
    for my $lst_lmn (@$lst){ # WCO00 p. 224
	if($dbg_lvl){printf STDOUT ("lst_add: lst_idx = %d, val= %g\n",$lst_idx,$lst_lmn);}
	if($lst_lmn != $mss_val){$lst_ttl+=$lst_lmn;}
	$lst_idx++
    } # end loop over lst
    if($dbg_lvl){printf STDOUT ("lst_add: lst_nbr = %d, mss_val = %g, lst_ttl = %g\n",$lst_nbr,$mss_val,$lst_ttl);}
    return $lst_ttl;
} # end lst_add()

sub yyyymmdd_ncr($$){
# Purpose: Increment a date by a given number of days
# Usage: $yyyymmdd_new=yyyymmdd_ncr($yyyymmdd_srt,$day_nbr);
# $day_nbr may be positive (days in the future) or negative (days in the past)
# $day_nbr should be an integer
    my ($yyyymmdd_srt,$day_nbr);
    ($yyyymmdd_srt,$day_nbr)=@_;
    my ($yyyy_srt,$mm_srt,$dd_srt);
    my ($yyyymmdd_new);
    my ($yyyy_new,$mm_new,$dd_new);
# Parse current date
    ($yyyy_srt,$mm_srt,$dd_srt)=&yyyymmdd_prs($yyyymmdd_srt);
# Set defaults for new date
    $yyyy_new=$yyyy_srt;
    $mm_new=$mm_srt;
    $dd_new=$dd_srt;
# Add days to current date
    $dpm_srt=dpmyr($mm_srt,$yyyy_srt);
    $dd_prv=$dd_srt+$day_nbr; # [day] Provisional day
    if($dd_prv <= $dpm_srt && $dd_prv >= 1){
	$dd_new=$dd_prv;
	$yyyymmdd_new=sprintf("%04d%02d%02d",$yyyy_new,$mm_new,$dd_new);
    }else{
	$doy_srt=yyyymmdd2doy($yyyymmdd_srt);
	$doy_new=$doy_srt+$day_nbr;
	$dpy_new=dpy_get($yyyy_new);
	if($day_nbr >= 0){
	    while($doy_new > $dpy_new){
		$doy_new-=$dpy_new;
		$yyyy_new++;
		$dpy_new=dpy_get($yyyy_new);
	    } # endwhile
	}else{ # endif $day_nbr >= 0
	    while($doy_new < 1){
		$yyyy_new--;
		$dpy_new=dpy_get($yyyy_new);
		$doy_new+=$dpy_new;
	    } # endwhile
	} # endif $day_nbr < 0
	$yyyymmdd_new=&yrdoy2yyyymmdd($yyyy_new,$doy_new);
    } # endelse
    return $yyyymmdd_new;
} # end yyyymmdd_ncr()

sub yrdoy2yyyymmdd($$){
# Purpose: Convert year and day of year [1.0..366.0) into yyyymmdd format
# yrdoy2yyyymmdd() is inverse of yyyymmdd2doy()
# Usage: $yyyymmdd=&yrdoy2yyyymmdd($yyyy,$doy);
    my ($yyyy,$doy); # I [day] Day of year
    my ($yyyymmdd); # O Date determined by yyyy and doy inputs
    my ($mm,$dd); # 
    my (@day_ttl_prc_mth,@day_ttl_prc_mth_non_leap_yr,@day_ttl_prc_mth_leap_yr); # Number of days in all preceding months
    ($yyyy,$doy)=@_;
# Days in preceding months
    @day_ttl_prc_mth_non_leap_yr=(0,31,59,90,120,151,181,212,243,273,304,334); # Number of days in all preceding months of a non-leap year
    @day_ttl_prc_mth_leap_yr=(0,31,60,91,121,152,182,213,244,274,305,335); # Number of days in all preceding months of a leap year
# Vet input
    if($doy > dpy_get($yyyy)){
	printf STDOUT "ERROR yrdoy2yyyymmdd() reports \$doy = $doy in \$yyyy = $yyyy\n";
	die;
    } # endif
# Compute calendar for this year
    my $leap_yr_flg=is_leap_yr($yyyy);
    if($leap_yr_flg){@day_ttl_prc_mth=@day_ttl_prc_mth_leap_yr;}else{@day_ttl_prc_mth=@day_ttl_prc_mth_non_leap_yr;}
    $mm=1; # Start at January
    while($mm < 13 && $doy > $day_ttl_prc_mth[$mm-1]){ # -1 is for 0-based indices
	$mm++;
    } # end while mm
    $mm--; # Compensate for overshoot
    $dd=$doy-$day_ttl_prc_mth[$mm-1]; 
    if($dd < 1 || $dd > dpmyr($mm,$yyyy)){printf STDOUT "ERROR yrdoy2yyyymmdd() reports \$yyyy = $yyyy, \$doy = $doy, \$mm = $mm, \$dd = $dd, \$day_ttl_prc_mth[\$mm-1] = $day_ttl_prc_mth[$mm-1]\n";}
    $yyyymmdd=sprintf("%04d%02d%02d",$yyyy,$mm,$dd);
    return $yyyymmdd;
} # end yrdoy2yyyymmdd()

sub yyyymmdd2doy($){
# Purpose: Convert date in yyyymmdd format to day of year [1.0..366.0)
# yyyymmdd2doy() is inverse of yrdoy2yyyymmdd()
# Usage: $doy=&yyyymmdd2doy($yyyymmdd);
    my ($yyyymmdd);
    my ($yyyy,$mm,$dd);
    my ($doy);
    ($yyyymmdd)=@_;
# Parse current date
    ($yyyy,$mm,$dd)=&yyyymmdd_prs($yyyymmdd);
# Add up days in preceding months
    @day_ttl_prc_mth_non_leap_yr=(0,31,59,90,120,151,181,212,243,273,304,334); # Number of days in all preceding months of a non-leap year
    $doy=$dd+$day_ttl_prc_mth_non_leap_yr[$mm-1];
# Correct for leap years
    if($mm > 2){
	if(is_leap_yr($yyyy)){$doy+=1;}
    } # endif leap year
#    print STDOUT "yyyymmdd2doy() reports \$yyyymmdd = $yyyymmdd, \$doy = $doy\n";
    return $doy;
} # end yyyymmdd2doy()

sub date_time(){
# Purpose: Returns local time string in "Sun Sep 18 20:14:49 1994" format
# Usage: 
    my ($sec,$min,$hr,$mday,$mon,$yr,$wday,$yday,$isdst);
    ($sec,$min,$hr,$mday,$mon,$yr,$wday,$yday,$isdst)=localtime(time);
    my $date_tm_sng=
	sprintf
	    ('%s %s %02d %02d:%02d:%02d 19%02d',
	     (Sun,Mon,Tue,Wed,Thu,Fri,Sat)[$wday],
	     (Jan,Feb,Mar,Apr,May,Jun,Jul,Aug,Sep,Oct,Nov,Dec)[$mon],
	     $mday,
	     $hr,
	     $min,
	     $sec,
	     $yr);
    return $date_tm_sng;
} # end date_time()

sub dpm($){
# Purpose: Returns # of days in given month (no leap days)
# Assume $mm is a two digit integer ranging from 1..12
# Usage: $dpm=dpm($mm)
    my ($mm)=@_; # Assume $mm is a two digit integer ranging from 1..12
    @dpm=(31,28,31,30,31,30,31,31,30,31,30,31);
    return $dpm[$mm-1]; # Subtract one for base-0 indexing
} # end dpm()

sub dpy_get($){
# Purpose: Returns # of days in given year (365 or 366)
# Usage: $dpy=&dpy_get($yyyy)
    my ($yyyy)=@_;
    my $dpy;
    return is_leap_yr($yyyy) ? 366 : 365;
} # end dpy_get()

sub is_leap_yr($){
# Purpose: Returns true (1) if given year is a leap year else returns false (0)
# Usage: $leap_yr_flg=is_leap_yr($yyyy)
    my ($yyyy)=@_;
    $leap_yr_flg=0; 
    if($yyyy%4 == 0 && $yyyy%100 != 0){$leap_yr_flg=1;} # 1900 is not a leap year
    if($yyyy%400 == 0){$leap_yr_flg=1;} # 2000 is a leap year
    return $leap_yr_flg;
} # end is_leap_yr()

sub dpmyr($$){
# Purpose: Returns # of days in given month (accounts for leap days)
# Usage: dpm=dpmyr($mm,$yyyy)
# Assume $mm is a two digit integer ranging from 1..12
# Assume $yyyy is in YYYY format
    my ($mm,$yyyy);
    ($mm,$yyyy)=@_;
    @dpm=(31,28,31,30,31,30,31,31,30,31,30,31);
    my $dpm=$dpm[$mm-1]; # Subtract one for base-0 indexing
    if($mm == 2){
	if(is_leap_yr($yyyy)){$dpm+=1;}
    } # endif February
    return $dpm;
} # end dpmyr()

sub YYYYMMDD_HHMM(){
# Purpose: Returns local time sng in "19940312_1209" format
# Usage: 
    local($sec,$min,$hr,$mday,$mth,$yr,$wday,$yday,$isdst);
    ($sec,$min,$hr,$mday,$mth,$yr,$wday,$yday,$isdst)=localtime(time);
    $YYYYMMDD_HHMM_sng=
        sprintf
            ("19%02d%02d%02d_%02d%02d",
	     $yr,
             $mth+1,
             $mday,
             $hr,
             $min);
    return $YYYYMMDD_HHMM_sng;
} # end sub YYYYMMDD_HHMM()

sub YYYYMMDD(){
# Purpose: Returns local time sng in "19640312" format
# Usage: 
    local($sec,$min,$hr,$mday,$mth,$yr,$wday,$yday,$isdst);
    ($sec,$min,$hr,$mday,$mth,$yr,$wday,$yday,$isdst)=localtime(time);
    $YYYYMMDD_sng=
        sprintf
            ("19%02d%02d%02d",
	     $yr,
             $mth+1,
             $mday);
    return $YYYYMMDD_sng;
} # end sub YYYYMMDD()

# Files included with 'require' must return TRUE as the last statement.
# Not doing this causes errors like "foo.pl did not return a true value"
1; 
