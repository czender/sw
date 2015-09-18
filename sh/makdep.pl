#!/usr/local/bin/perl

# Generate dependencies in a form suitable for inclusion into a Makefile.
# The source filenames are provided in a file, one per line.  Directories
# to be searched for the source files and for their dependencies are provided
# in another file, one per line.  Output is written to STDOUT.
#
# For CPP type dependencies (lines beginning with #include) the dependency
# search is recursive.  Only dependencies that are found in the specified
# directories are included.  So, for example, the standard include file
# stdio.h would not be included as a dependency unless /usr/include were
# one of the specified directories to be searched.
#
# For Fortran module USE dependencies (lines beginning with a case
# insensitive "USE", possibly preceded by whitespace) the Fortran compiler
# must be able to access the .mod file associated with the .o file that
# contains the module.  In order to correctly generate these dependencies
# two restrictions must be observed.
# 1) All modules must be contained in files that have the same base name as
#    the module, in a case insensitive sense.
# 2) All modules that are to be contained in the dependency list must be
#    contained in one of the source files in the list provided on the command
#    line.
# The reason for the second restriction is that since the makefile doesn't
# contain rules to build .mod files the dependency takes the form of the .o
# file that contains the module.  If a module is being used for which the
# source code is not available (e.g., a module from a library), then adding
# a .o dependency for that module is a mistake because make will attempt to
# build that .o file, and will fail if the source code is not available.
#
# Author: B. Eaton
#         Climate Modelling Section, NCAR
#         Feb 2001

use Getopt::Std;
use File::Basename;

# Check for usage request.
@ARGV >= 2                          or usage();

# Process command line.
%opt = ();
getopts( "t:", \%opt )            or usage();
$filepaths = shift()              or usage();
$srcfiles = shift()               or usage();
@ARGV == 0                        or usage();  # Check that all args were processed.

if ( defined $opt{'t'} ) { $fmt = $opt{'t'}; }
open(FILEPATH, $filepaths) or die "Can't open $filepaths: $!\n";
open(SRCFILES, $srcfiles) or die "Can't open $srcfiles: $!\n";

# Make list of paths to use when looking for files.  Remove newline characters.
# Prepend "." so search starts in current directory.  This default is for
# consistency with the way GNU Make searches for dependencies.
@paths = <FILEPATH>;
chomp @paths;
unshift(@paths,'.');

# Make list of files containing source code.
@src = <SRCFILES>;
chomp @src;

# Generate a hash of the source filenames in order to efficiently search for the 
# file that contains a module in a case insensitive way.
@suffixes = ('\.[fF]90', '\.[fF]' );
foreach $f (@src) {
    ($name, $path, $suffix) = fileparse($f, @suffixes);
    ($mod = $name) =~ tr/a-z/A-Z/;  # convert file's basename to uppercase
    $srchash{$mod} = $name;
}
#while ( ($k,$v) = each %srchash ) {
#    print "$k => $v\n";
#}

# Write dependencies formatted for inclusion into a Makefile to STDOUT.
FILE:
foreach $file ( @src ) {

    local( $list, $x, $y, $target, @dep );

    @dep = (); # make @dep the null list

    @list = dependents( $file );

    # Check the list returned from dependents for error flags and remove
    # redundant dependencies.

  LIST_MEMBER:
    for ( $i = 0; $i <= $#list; ++$i ) {
	$x = $list[$i];
	if ( $x == -1 ) {                            # error return code
	    if ( $file eq $list[$i+1] ) {
		print STDERR "$file not found\n";
		next FILE;
	    } else {
		++$i;
		print STDERR "$file dependency $list[$i] not found\n";
		# Remove dependency not found from the dependency list
		# (or make will try to create it!).
		# It is last one in @dep since it was added to the list before the
		# function "dependents" tried to find it.
		pop( @dep );
		next LIST_MEMBER;
	    }
	} else {        # add filename to dependency list - check for duplication.
	    # If the dependency is a ".o" file, add the appropriate formatting
	    if ( $x =~ /\w+\.o/ ) {
		if ( $fmt eq "OBJ" ) { $x = "OBJ/$x"; }
	    }
	    foreach $y ( @dep ) {
		if ( $y eq $x ) { 
		    next LIST_MEMBER;
		}
	    }
	    push( @dep, $x );
	}
    }

    # No errors were encountered... format the list for inclusion in makefile.
    $file =~ /\s*(\w[^.]*)/;
    $target = "$1.o";
    if ( $fmt eq "OBJ" ) { $target = "OBJ/$1.o"; }

    print "$target : $file @dep\n";

}

#--------------------------------------------------------------------------------------

sub findsrc {

# Search for the specified file in the list of directories in the global
# array @path.  Return the first occurance found, or the null string if
# the file is not found.

    local( $file ) = @_;
    local( $dir, $fname );

    foreach $dir (@paths) {

	if( $dir =~ m#/$# ) {             # allow directory name to end with /
	    $fname = $dir . $file;
	} else {
	    $fname = $dir . '/' . $file;
	}

	if ( -e  $fname ) {
	    return $fname;
	}
      
    }
    return '';  # file not found
}

#--------------------------------------------------------------------------------------

sub dependents {

    # Search recursively for all files that are "#include"d by the cpp
    # preprocessor.  

    local( $file ) = @_;
    local( @out, $fh );

    # Find the file.

    if( ! ($absname = findsrc( $file )) ) {
	return -1, $file; # file not found
    }

    # Make a unique filehandle.

    $fh = $file . 'FH';
    $fh =~ tr/a-z/A-Z/;

    open( $fh, $absname );

    while ( <$fh> ) {
	if ( /^#include\s+[<"](.*)[>"]/ ) {
            # Search for "#include" and strip filename when found.  Append the file to
	    # the output array and descend into that file to continue the search.
            push @out, $1;
	    push @out, &dependents( $1 );
	} 
	elsif ( /^\s*[Uu][Ss][Ee]\s+(\w+)/ ) {
	    # Search for module dependencies.
	    # Since Fortran is case insensitive must check the source filenames to
	    # find the file that contains the module in a case insensitive manner.
	    # Don't need to do recursive search on modules because the dependency is
	    # on a .o file and we assume that that .o file is created by one on the
	    # files in the Srcfiles list, and hence has its own dependencies generated.
	    ($module = $1) =~ tr/a-z/A-Z/;
	    if ( defined $srchash{$module} ) { push @out, "$srchash{$module}.o"; }
	}
    }
    close( $fh );
    return @out;
}

#--------------------------------------------------------------------------------------

sub usage {
    ($ProgName = $0) =~ s!.*/!!;            # name of program
    die <<EOF
SYNOPSIS
     $ProgName [-t fmt] Filepath Srcfiles
OPTIONS
     -t fmt
          Target format.  The only valid value for fmt is OBJ.  If this
          option is set the .o files that are targets in the dependency
          rules have the form OBJ/file.o.
ARGUMENTS
     Filepath is the name of a file containing the directories (one per 
     line) to be searched for dependencies.  Srcfiles is the name of a
     file containing the names of files (one per line) for which
     dependencies will be generated.
EOF
}
