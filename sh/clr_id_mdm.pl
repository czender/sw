#!/usr/local/bin/perl5 -w

# Purpose: Get computer to perform caller ID on incoming calls

# History: Began life as cid-modem by jwz
# Caller-ID-logger, by jwz (19-Jan-97)
# Modified: 22-Sep-99
# Downloaded and renamed clr_id_mdm.pl by CSZ on 20001022

# Opens the modem and waits for it to print caller-ID data.  When it does,
# it logs it to a file, parses it, and pops up a window using "xmessage".
# If the number is present in your .bbdb file, it shows the name (or company)
# associated with it.  

# Todo:
#  Currently BBDB only gets the phone number, not the name.  It should look
#  up the name if the number doesn't match.

#  Modems other than ZyXELs have different caller-ID formats, and this doesn't
#  deal with those.

# Set this to the device that your modem is attached to
$modem_device = "/dev/ttyd1";

# kludge...
if ( ! $ENV{HOME} ) { $ENV{HOME} = "/home/zender"; }

# This is your .bbdb file.  (Set it to null if you don't want to do BBDB
# lookups at all, but why would you want to go and do a thing like that?)
$bbdb_file    = "$ENV{HOME}/.bbdb";

# A shell command to use to cause emacs to pop up the BBDB buffer
# (bbdb-srv.pl is a good choice, so it defaults to the value of the
# shell environment variable $NS_MSG_DISPLAY_HOOK.)
$bbdb_cmd     = $ENV{NS_MSG_DISPLAY_HOOK};

# If you want the $bbdb_cmd to be run on a different host, set it here
$bbdb_host    = "grendel";
$rsh_command  = "ssh -n";

# If you want each call to be logged to a file as well, name it here
$logfile      = "/usr/local/etc/clr_id_log.txt";

# For verbosity...
$debug        = 1;

# How to pop up a dialog box
$xmessage_cmd	= "xmessage";
@xmessage_args	= ("-display",	":0",
		   "-name",	"Caller ID",
		   # SGI lossage
		   "-xrm",	"*useSchemes: none",
		   # roughly centered on my screen; YMMV.
#		   "-geometry",	"+400+400",
		   "-geometry",	"+260+320",
#		   "-xrm",	"*Font: -*-new cent*-bold-r-normal-*-240-*",
		   "-xrm",
		       "*Font: -*-new cent*-bold-r-normal-*-*-512-*-*-*-*-*-*",
		   "-xrm",	"*Foreground: black",
		   "-xrm",	"*Background: lightgreen",
		   # no buttons on the window: dismiss it by clicking in it.
		   "-button",	"",
		   "-xrm", "*form.Translations: #override <BtnDown>: exit(0)",
#		   "-xrm", "*Command.Font: -*-new cent*-bold-r-normal-*-120-*",
		   "-xrm",
	       "*Command.Font: -*-new cent*-bold-r-normal-*-*-240-*-*-*-*-*-*",
		   "-xrm", "*Command.horizDistance: 200",
		   "-xrm", "*internalWidth: 25",
		   "-xrm", "*internalHeight:25"
		   );

# Do you want to turn off the screensaver before popping up the window?
# I've decided that I don't
#$pre_dialog_cmd = "xscreensaver-command -deactivate";
$pre_dialog_cmd = 0;

# To make the computer do intelligent ringing instead of the phone
# (This could work better.)
$ring_cmd = "/u/jwz/src/cid/ring.sh";


# commands (and their expected responses) used to initialize the modem
@modem_init   = ( "AT",		"OK",		# ping
		  "ATZ",	"OK",		# reset
		  "ATE0",	"OK",		# don't echo commands
		  "ATM0",	"OK",		# turn off speaker
		  "ATN0",	"OK",		# turn off ringer
		  "ATS40.2=1",	"OK",		# turn on caller ID
	        );


# for diagnostics: if the modem ever asynchronously prints something that
# doesn't match this, we issue a warning.
$expected_responses = "^CALLER NUMBER"		. "|" .
    		      "^CALLER NAME"		. "|" .
    		      "^REASON FOR NO CALLER "	. "|" .
    		      "^RING\$"			. "|" .
    		      "^TIME: [-0-9: ]+\$";


##############################################################################
# Talking to the serial port...
if ( $ debug ) {
    use diagnostics;
}


sub DEBUG {
    if ( $debug ) {
        my ($str) = @_;
        my $d = localtime;
        print STDERR "$d:   $str\n";
    }
}

sub open_modem {
    use IPC::Open2;

    # Close the terminal streams before forking `cu', because otherwise
    # it fucks around with the stty settings.
    #
# 8-Dec-98 -- for some reason this doesn't work any more; it makes the 
# forked "cu" no longer function.  What changed?  All I can think of is that
# I installed a newer version of perl.

#    open(SAVEIN,  "<&STDIN")  || die("can't dup stdin");
#    open(SAVEOUT, ">&STDOUT") || die("can't dup stdout");
#    open(SAVEERR, ">&STDERR") || die("can't dup stderr");
#    close(STDIN);
#    close(STDOUT);
#    close(STDERR);

    my $cu_pid = open2( \*MODEM_IN, \*MODEM_OUT,
		       "cu -l$modem_device -s2400 2>&1");

    # Now that cu has been launched, we can restore them.
    #
#    open(STDIN,  "<&SAVEIN")  || die("can't restore stdin");
#    open(STDOUT, ">&SAVEOUT") || die("can't restore stdout");
#    open(STDERR, ">&SAVEERR") || die("can't restore stderr");
#    close(SAVEIN);
#    close(SAVEOUT);
#    close(SAVEERR);

    # The following doesn't seem to work and I don't know why...
    #
    # Set up a signal handler to try and kill off the cu process
    # when we exit, instead of waiting ~30 seconds for it to notice
    # that the pipe is gone...
    #
#    $SIG{INT} = sub { my $signame = shift;
#		      DEBUG "sending $signame to $cu_pid";
#		      }
#		      print MODEM_OUT "\r\n~.\r\n";
#		      close MODEM_OUT;
#		      close MODEM_IN;
#		      kill ($signame, $cu_pid);
#		      exit (1);
#		    };

    $_ = <MODEM_IN>;
    chop;
    if ( !m/^Connected/ ) {
	DEBUG "$0: cu printed `$_' instead of `Connected'";
    }
}

sub read_line {
    $_ = <MODEM_IN>;
    $_ || die("got eof on modem");
    s/[\r\n]+$//;
    if ( $_ eq "" ) {
	$_ = <MODEM_IN>;
	$_ || die("got eof on modem");
	s/[\r\n]+$//;
    }
    return $_;
}

sub command {
    my ( $command, $expected_response) = @_;

    if ( $debug ) {
	DEBUG "sending: $command";
    }

    print MODEM_OUT "$command\r\n";
    my $line = read_line();

    if ( $line eq $command ) {
	if ( $debug ) {
	    DEBUG "    got echo: reading next line too...";
	}
	$line = read_line();
    }

    if ( $line ne $expected_response ) {
	DEBUG "    got: $line ; expected: $expected_response";
    } elsif ( $debug ) {
	DEBUG "    got: $line";
    }
}

sub init_modem {
    open_modem;

    my $len = $#modem_init + 1;
    my $i;
    for ($i = 0; $i < $len; $i += 2) {
	command($modem_init[$i], $modem_init[$i+1]);
    }
}

sub handle_async_line {
    local ( $_ ) = @_;

    if (!m/$expected_responses/) {
	DEBUG "modem turd:   $_";

    } elsif (m/CALLER NUMBER/) {
	if ( $debug ) {
	    DEBUG "number: $_";
	}
	handle_cid_line($_);

    } elsif (m/CALLER NAME/) {
	if ( $debug ) {
	    DEBUG "name:   $_";
	}
	handle_cname_line($_);

    } elsif (m/^RING$/) {
	if ( $debug ) {
	    DEBUG "ring: $_";
	}
	handle_ring_line($_);

    } elsif ( $debug ) {
	if ( $_ eq '' ) {
	    DEBUG "ignored: blank line";
	} else {
	    DEBUG "ignored: $_";
	}
    }
}


##############################################################################
#
# Parsing BBDB and CID data...
#

sub find_bbdb_record {
    my ( $area, $exchange, $suffix ) = @_;

    if ( ! $bbdb_file ) {
	return undef;
    }

    # strip off leading 0's, to match the way it's stored in .bbdb.
    $area     =~ s/^0+(.)/$1/;
    $exchange =~ s/^0+(.)/$1/;
    $suffix   =~ s/^0+(.)/$1/;

    my $which = undef;
    my $bbdb_rec = undef;
    my $pat = "\\[\"([^\"]+)\" $area $exchange $suffix (nil|[0-9]+)\\]";

    open(BBDB, "<$bbdb_file") || die("error opening $bbdb_file: $!\n");

    while (<BBDB>) {
	if ( m/$pat/ ) {
	    $bbdb_rec = $_;
            $which = $1;
	    last;
	}
    }
    close(BBDB);
    return ($bbdb_rec, $which);
}


# note: global (kludge!)
$pretty_number = 0;

sub make_message_string {
    my ( $number, $date, $fn, $ln, $co, $which, $caller_name, $error ) = @_;

    my $msg_date;
    my $msg_name;
    my $msg_name2;
    my $msg_number;

    my $line_prefix = " ";
    my $line_suffix = " ";


    if ( $caller_name && $caller_name eq "" ) {
        $caller_name = 0;
    }

    # First reformat the date.
    #
    $_ = $date;
    my ( $dotw, $mon, $day, $hr, $min, $sec, $year ) =
	m/^([^ ]+) +([^ ]+) +([^ ]+) +([^:]+):([^:]+):([^:]+) +([^ ]+) *$/;
    $year =~ s/^..(..)/$1/;
    $day  =~ s/^0//;
    $hr   =~ s/^0//;
    if ($hr < 12) {
	$ampm = "am";
        if ($hr == 0) { $hr = 12; };
    } else {
	$ampm = "pm";
	if ($hr > 12) { $hr -= 12; };
    }
    $msg_date = "$hr:$min$ampm, $day $mon ($dotw)";


    # Next format the caller name, company, or error message.
    # (This is the BBDB version of the name.)
    #
    if ( $error ) {
	$msg_name = "$error";
    } elsif ( $co && !$fn && !$ln ) {
	$msg_name = "$co";
    } elsif ( $fn || $ln ) {
	$msg_name = "$fn $ln";
    }

    # Next format the telephone company's version of the name.
    # (This is the version we read off the line, if any.)
    #
    $msg_name2 = "";
    if ( $caller_name ) {

        $_ = $caller_name;

        s/(\w+)/\u\L$1/g;		# capitalize line.
        s/^([^ ]+) (.+)$/$2 $1/;	# move first word to end.

        if ( ! $fn ) { $fn = ""; }
        if ( ! $ln ) { $ln = ""; }
        if ( ! m/^$fn $ln$/i ) {	# skip if exact match
            $msg_name2 = "$_";
        }
    }

    # Next format the phone number.
    #
    $pretty_number = 0;
    if ( $number ) {
	my $area = 0;
	my $exchange = 0;
	my $suffix = 0;
	$_ = $number;
	( $area, $exchange, $suffix ) =
	    m/^([0-9][0-9][0-9])([0-9][0-9][0-9])([0-9][0-9][0-9][0-9]+)/;

	# note: global (kludge!)
	$pretty_number = "($area) $exchange-$suffix";
    }

    if ($which) {
        if (length($which) >= 8) {
            $pretty_number .= "$line_suffix\n$line_prefix($which)";
        } else {
            $pretty_number .= " ($which)";
        }
    }

    my $msg;
    my $short_msg;

    $msg = $line_prefix . $msg_date . $line_suffix;
    $msg .= "\n" . $line_prefix . $msg_name . $line_suffix;

    $short_msg = "$msg_date: $msg_name";

    if ( $msg_name2 && $msg_name2 ne "" ) {
	$msg .= "\n" . $line_prefix . "($msg_name2)" . $line_suffix;
	$short_msg .= " ($msg_name2)";
    }

    if ( $pretty_number && $pretty_number ne "" ) {
	$msg .= "\n" . $line_prefix . $pretty_number . $line_suffix;
	$short_msg .= "; $pretty_number";
    }

    # Log the short message too...
    #
    if ( $logfile ) {
	if (open(LOG, ">>$logfile")) {
	    print LOG "** $short_msg\n";
	    close LOG;
	} else {
	    DEBUG "error opening $logfile: $!";
	}
    }


    return $msg;
}

use POSIX;
sub reaper {
    $SIG{CHLD} = \&reaper;  # loathe sysV
    my $signame = shift;
    if ( $debug >= 2 ) {
	DEBUG "  (got SIG$signame...)";
    }
    my $child;
    while ( ( $child = waitpid(-1,WNOHANG) ),
	    $child > 0 ) {
	if ( $debug >= 2 ) {
	    DEBUG "    (pid $child exited with $?)";
	}
    }
}

sub fork_and_exec {
    my @cmd_list = @_;

    $SIG{CHLD} = \&reaper;

    if ( $debug >= 2 ) {
	$_ = $cmd_list[0];
	s/ .*//;
	DEBUG "forking for " . $_ . " at " . (localtime) . ".";
    }

    my $pid;
    if ($pid = fork()) {
	# parent
    } elsif (!defined $pid) {
	DEBUG "$0: fork failed: $!";
    } else {
	# child

	if ( $debug ) {
	    $_ = $cmd_list[0];
	    s/ .*//;
	    DEBUG "exec'ing $_ in pid $$.";
	}
	close(STDIN);
#	close(STDOUT);
#	close(STDERR);
	exec @cmd_list;
    }
}


sub fork_and_exec_for_bbdb {
    my @cmd_list = @_;
    my $number = shift @cmd_list;

    $SIG{CHLD} = \&reaper;

    if ( $debug >= 2 ) {
	$_ = $cmd_list[0];
	s/ .*//;
	DEBUG "forking for " . $_ . " at " . (localtime) . ".";
    }

    my $pid;
    if ($pid = fork()) {
	# parent
    } elsif (!defined $pid) {
	DEBUG "$0: fork failed: $!";
	exit (1);
    } else {
	# child

	if ( $debug ) {
	    $_ = $cmd_list[0];
	    s/ .*//;
	    DEBUG "exec'ing $_ in pid $$.";
	}
	if ( system @cmd_list ) {
	    my $cmd = "gnudoit -q '(bbdb-srv-add-phone \"$pretty_number\")'";
	    if ( $bbdb_host ) {
		$cmd =~ s/([()\"])/\\$1/g;
		$cmd = "$rsh_command $bbdb_host $cmd";
	    }
	    exec $cmd;
	}
	exit (0);
    }
}


sub pop_up_dialog {
    my ( $msg, $buttonp, $number ) = @_;

    if ( $pre_dialog_cmd ) {
	fork_and_exec $pre_dialog_cmd;
    }

#    $msg = "\n$msg\n\n";
    if ( ! $buttonp ) {
	fork_and_exec $xmessage_cmd, @xmessage_args, $msg;
    } else {
	my @args = ( @xmessage_args, "-button", "Add To BBDB" );
	fork_and_exec_for_bbdb $number, $xmessage_cmd, @args, $msg;
    }
}

sub pop_up_bbdb_buffer {
    my ( $caller ) = @_;
    if ( $bbdb_cmd ) {
	my $cmd = $bbdb_cmd;
	if ( $bbdb_host ) {
	    $cmd = "$rsh_command $bbdb_host $cmd";
	}
	$caller =~ s/\\/\\\\/g;
	$caller =~ s/\"/\\\\\"/g;
	`echo "Path:\nFrom: \\\"$caller\\\" <>" | $cmd >&- 2>&- &`;
    }
}


sub handle_cid_line {
    my($line) = @_;

    my $date = localtime;

    # Log the call...
    #
    if ( $logfile ) {
	if (open(LOG, ">>$logfile")) {
	    print LOG "$date\t$line\n";
	    close LOG;
	} else {
	    DEBUG "error opening $logfile: $!";
	}
    }

    # Pull the phone number out of the message...
    #
    my $number = "";
    my $error = "";

    $_ = $line;
    if ( m/^CALLER NUMBER/ ) {
	( $number ) = m/^[^:]+: *(.*) *$/;
    } else {
	$error = $line;
    }

    my $caller = undef;

    my $fn = undef;
    my $ln = undef;
    my $co = undef;
    my $which = undef;
    my $buttonp = 0;

    if ( !$number || $number eq "" ) {
       	$error =~ tr#A-Z#a-z#;
	$error =~ s/^REASON FOR NO CALLER (NUMBER|NAME)/Caller unknown/i;

    } else {
	$_ = $number;

	my $area = 0;
	my $exchange = 0;
	my $suffix = 0;
	( $area, $exchange, $suffix ) =
	    m/^([0-9][0-9][0-9])([0-9][0-9][0-9])([0-9][0-9][0-9][0-9]+)/;

	my $bbdb_rec;
        ($bbdb_rec, $which) = find_bbdb_record($area, $exchange, $suffix);

	if ( $bbdb_rec ) {
	    my $junk = 0;
	    $_ = $bbdb_rec;
	    # This will lose if names or aliases have double-quotes in them.
	    # No doubt there's some hairier regexp magic that handles that...
	    ( $fn, $ln ) = m/^[\[]\"([^\"]*)\" *\"([^\"]*)\"/;

            # handle folks with only one name
            if ( !$fn && !$ln ) {
                ( $fn ) = m/^[\[]\"([^\"]*)\" nil /;
                if ($fn) { $ln = ""; }
            }

	    ( $junk, $junk, $junk, $co ) =
      m/^[[](nil|\"[^\"]*\") *(nil|\"[^\"]*\") (nil|[(][^)]*[)]) \"([^\"]*)\"/;

	    if ( $fn || $ln ) {
		$caller = "$fn $ln";
	    }

            if ( $which ) {
                $_ = $which;
                s/[ \t][ \t]+/ /g;       # compress whitespace
                s/\(cid\)//g;            # lose "cid" in parens
                s/\bcid\b//g;            # lose "cid" as a word
                s/^ *\((.*)\) *$/$1/;    # "(foo)" => "foo"
                s/\(.*?\)//g;            # lose anything in parens
                s/[ \t][ \t]+/ /g;       # compress whitespace
                s/^ //; s/ $//;          # lose leading/trailing whitespace
                $which = $_;
                if ($which =~ /^$/) {
                    $which = undef;
                }
            }

	} else {
	    $buttonp = 1;
	}
    }

    my $msg = make_message_string($number, $date, $fn, $ln, $co, $which,
                                  $caller_name, $error);
    pop_up_dialog($msg, $buttonp, $pretty_number);

    if ( $caller ) {
	pop_up_bbdb_buffer($caller);
    }

    # clear it for next time
    $caller_name = "";
}

sub handle_cname_line {
    my($line) = @_;

    my $date = localtime;

    # Log the call...
    #
    if ( $logfile ) {
	if (open(LOG, ">>$logfile")) {
	    print LOG "$date\t$line\n";
	    close LOG;
	} else {
	    DEBUG "error opening $logfile: $!";
	}
    }

    # Pull the caller name out of the message...
    #
    $caller_name = "";

    $_ = $line;
    if ( m/^CALLER NAME/ ) {
	( $caller_name ) = m/^[^:]+: *(.*) *$/;
        $caller_name =~ s/[ \r\n]*$//;
    }
}


# Make a phone-ringing sound.
#
# Rather than just letting ring.sh decide how loud to ring it each time,
# ask ring.sh how loud it would do it the first time, then do all the
# rings at that volume.  That's because the answer will change when we
# turn off the screensaver...
#
# Wait N minutes before asking again, where N is hopefully longer than
# the screensaver timeout value -- and therefore, if the console is idle,
# the screensaver will have reactivated before we ask again.
#
# The other way in which this sucks is that the rings that the computer
# makes are delayed by about two seconds from the rings that the phone
# would make.  This is is because the modem doesn't actually print a
# "RING" line until the ring has *completed*.
#
# TODO: Given that we don't end up processing the first RING line until just 
#       before we process the caller ID info (they come within 1/8 second of
#       each other) we might as well ring differently for different cid info,
#       right?  Like, "don't ring if privacy or out-of-area" or something.
#
#       The hard part here is, how do you detect "end of call"?  Maybe by
#       noticing that more than N seconds has passed between RING lines?
#       That would only do the wrong thing if one call came right on the
#       tail of another (but the second cid line would be the clue there.)

my $be_quiet = 0;
my $be_quiet_time = 0;

sub handle_ring_line {
    my($line) = @_;

    if ( $pre_dialog_cmd ) {

	my $minutes = 6;	# duration of cache
	if ( $be_quiet_time < (time - ($minutes * 60)) ) {
	    $be_quiet_time = time;
	    my $q = `$ring_cmd -query`;
	    chop $q;
	    $be_quiet = ( "$q" eq "quiet" );
	}

	if ( $ring_cmd ) {
	    if ( $be_quiet ) {
		fork_and_exec "$ring_cmd -quiet";
	    } else {
		fork_and_exec "$ring_cmd -loud";
	    }
	}

    } else {		# the simple way...
	fork_and_exec "$ring_cmd";
    }
}



##############################################################################
#
# hey ho.  let's go.
#
sub main {
    init_modem();
    while (1) {
	handle_async_line(read_line());
    }
    exit (1);
}

main();
