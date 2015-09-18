#!/bin/sh

# $Id$

# Purpose: Test bash scripts

# Bash login shells execute one of (in order) .bash_profile, .bash_login, .profile
# .bashrc is executed at the start of any subshell
# .bash_logout is executed every time a Bash login shell exits
# /etc/bash_completion contains shortcute commands on Debian systems

# Bash implements many features not supported by the standard shell
# bash and ksh functions may be defined with either of two syntaxes:
# 1. fnc_nm () { foo ; }
# 2. function fnc_nm { foo ; }
# whereas sh on, e.g., Solaris, only supports the former method 
# Bash understands 'local' keyword, but ksh & sh do not

# Set debugging options to bash with NeR98 p. 221
# cmd line method | script method | Purpose
# bash -n           set -o noexec   Do not run commands, just look for syntax errors
# bash -v           set -o verbose  Echo commands before running them
# bash -x           set -o xtrace   Echo commands after command-line processing

# Usage: 
# ${DATA}/bashdb/bashdb ${HOME}/sh/tst_bsh.sh
# chmod a+x ~/sh/tst_bsh.sh;~/sh/tst_bsh.sh

# stdout and stderr redirection:
# ${HOME}/sh/tst_bsh.sh > ~/foo 2>&1 <-- Overwrite file foo if it exists
# ${HOME}/sh/tst_bsh.sh >> ~/foo 2>&1 <-- Append to file foo if it exists

fl_is_exe () {
# function fl_is_exe { # function fnc_nm { foo ; } works in bash & ksh but not sh
# Purpose: Check if file exists and is executable
# Usage: fl_is_exe $fl
    local fl
    fl=$1 # String substition operator NeR98 p. 95, 98
    fl=${fl:?"ERROR: fl_is_exe() requires filename as argument"}
    test -f "$fl" -a -x "$fl"
    # Functions return exit status of last command executed
} # end fl_is_exe()

if fl_is_exe ${HOME}/sh/pvmgetarch ; then
	export PVM_ARCH=`${HOME}/sh/pvmgetarch`
elif fl_is_exe /fs/cgd/home0/zender/sh/pvmgetarch ; then
	export PVM_ARCH=`/fs/cgd/home0/zender/sh/pvmgetarch`
else
	export PVM_ARCH='SUNMP'
fi # endif pvmgetarch

printf "Testing determination of \$HOST:\n"
if [ -z "${HOST}" ]; then
    if fl_is_exe /bin/hostname ; then
	export HOST=`/bin/hostname`
    elif fl_is_exe /usr/bin/hostname ; then
	export HOST=`/usr/bin/hostname`
    fi # endif hostname exists
fi # endif HOST
printf "HOST = ${HOST}\n"

printf "Testing if [ 1 ] then continue; fi:\n"
if [ 1 ]; then 
    continue
fi # endif

printf "Testing if [ -f file ] then...else...fi:\n"
if [ -f ${HOME}/.bashrc ]; then 
    printf "~/.bashrc is a file (correct)\n"
else 
    printf "~/.bashrc is not a file (incorrect)"
fi # endif

printf "Testing if [ ! -f file ] then...else...fi:\n"
if [ ! -f ${HOME}/.bashrc ]; then 
    printf "~/.bashrc is not a file (incorrect)"
else 
    printf "~/.bashrc is a file (correct)\n"
fi # endif

if [ "${PVM_ARCH}" = "SUN*" ]; then
    printf "PVM_ARCH = ${PVM_ARCH} matches string comparison to SUN*\n"
else 
    printf "PVM_ARCH = ${PVM_ARCH} does not match string comparison to SUN*\n"
fi # endif

# NB: Warning: Bash does not match string patterns in logical string comparisons
# Following block proves this and should not be copied as a template
if [ "${PVM_ARCH}" = "LIN*" ]; then
    printf "PVM_ARCH = ${PVM_ARCH} matches string equality comparison to LIN*\n"
else 
    printf "PVM_ARCH = ${PVM_ARCH} does not match string equality comparison to LIN*\n"
fi # endif

# NeR98 p. 99
if [ ${PVM_ARCH%LIN} ]; then
    printf "PVM_ARCH = ${PVM_ARCH} pattern-matches to LIN, i.e., ${PVM_ARCH} contains LIN\n"
else 
    printf "PVM_ARCH = ${PVM_ARCH} does not pattern match to LIN\n"
fi # endif

# NeR98 p. 99
CVS_DRC=':ext:zender@esmf.ess.uci.edu:/u/zender/cvs'
CVS_DRC=':ext:zender@esmf.ess.uci.edu:/u/scapps/cvs'
pattern='zender/cvs'
#pattern='foo'
if [ ${CVS_DRC%${pattern}} ]; then
    printf "if-statement pattern-matches ${CVS_DRC} to ${pattern}, i.e., ${CVS_DRC} contains ${pattern}\n"
else 
    printf "if-statement does not pattern match ${CVS_DRC} to ${pattern}\n"
fi # endif

case "${CVS_DRC}" in 
    *${pattern}* ) 
    printf "case-statement pattern-matches ${CVS_DRC} to ${pattern}, i.e., ${CVS_DRC} contains ${pattern}\n"
    ;; # endif match*
    * )
	printf "case-statement does not pattern match ${CVS_DRC} to ${pattern}\n"
    ;; # endif default
esac # endcase ${CVS_DRC}

printf "Testing a simple loop: for foo in foo bar ; do ; ... ; done:\n"
for foo in foo bar ; do 
    printf "Traversing loop element ${foo}\n"
done # end loop over foo

printf "Testing a while loop: while ; do ; ... ; done:\n"
while [ $# -gt 0 ]; do
    sleep 1
done # end loop over foo

case "${PVM_ARCH}" in 
    LIN* ) 
	printf "PVM_ARCH = ${PVM_ARCH} matches case LIN*\n"
    ;; # endif LIN*
    SUN* ) 
	printf "PVM_ARCH = ${PVM_ARCH} matches case SUN*\n"
    ;; # endif SUN*
    * )
	printf "PVM_ARCH = ${PVM_ARCH} had no case matches\n"
    ;; # endif default
esac # endcase ${PVM_ARCH}

if [ -n "${PVM_ARCH}" ]; then
    printf "'-n' test shows variable PVM_ARCH = ${PVM_ARCH} exists in environment\n"
else 
    printf "'-n' test shows variable PVM_ARCH does not exist in environment\n"
fi # endif

if [ '0' = '1' ]; then
printf "CVS Root files beneath ${PWD} are:\n"
for fl in `find . -name Root -print`; do # Space before semi-colon is unnecessary
    printf "${fl}\n"
done # end loop over fl
fi # endif

# Use path name manipulation NeR98 p. 100 to isolate file stub from path
if [ '0' = '1' ]; then
    for fl in `ls *.F90 *.c`; do
       makdep ${fl} > ${fl%.*}.d
    done
fi # endif

if [ '0' = '1' ]; then
    printf "Timing speeds of file transfer methods for 315 MB file\n"
    /bin/rm /data/zender/dstmch41/h00??.nc
    printf "Using /bin/cp across NFS cross-mounted partitions...\n"
    echo `date`
    /bin/cp /data/zender/ZENDER/match/dstmch41/hist/h0034.nc /data/zender/dstmch41/h0034.nc
    echo `date`
    printf "Using scp...\n"
    echo `date`
    scp krein.math.uci.edu:/data/zender/ZENDER/match/dstmch41/hist/h0033.nc /data/zender/dstmch41/h0033.nc
    echo `date`
fi # endif

