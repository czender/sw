# $CVS_Header='$Id$';

# Purpose: Systematic nomenclature for debugging/verbosity levels

# Example usage:
# use DBG; # Debugging constants

package DBG; # NB: Module name must be uppercase version of filename stub
require Exporter;
@ISA=qw(Exporter);
@EXPORT=qw($dbg_off $dbg_fl $dbg_scl $dbg_crr $dbg_sbr $dbg_io $dbg_vec $dbg_vrb $dbg_old); # Export by default
@EXPORT_OK=qw($dbg_nbr_max); # Export on request

$dbg_nbr_max=9; # Number of different debugging levels

$dbg_off=0; # Production mode. Debugging is turned off.
$dbg_fl=1; # Filenames
$dbg_scl=2; # Scalars
$dbg_crr=3; # Current task
$dbg_sbr=4; # Subroutine names on entry and exit
$dbg_io=5; # Subroutine I/O
$dbg_vec=6; # Entire vectors
$dbg_vrb=7; # Everything
$dbg_old=8; # Old debugging blocks not used anymore

# NB: This database was made by hand
%dbg_db=( # Hash of scalars
	 dbg_off => 'Production mode. Debugging is turned off.',
	 dbg_fl => 'Filenames',
	 dbg_scl => 'Scalars',
	 dbg_crr => 'Current task',
	 dbg_sbr => 'Subroutine names on entry and exit',
	 dbg_io => 'Subroutine I/O',
	 dbg_vec => 'Entire vectors',
	 dbg_vrb => 'Everything',
	 dbg_old => 'Old debugging blocks not used anymore',
); # end %dbg_db
