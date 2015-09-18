#!/bin/csh -f
# -f = fast start, ignore .cshrc and .login
# -x = echo commands after substitutions and just before execution
# -fx = fast with echoing (use only one dash)

# $Id$

# Purpose: Test C-shell scripts

# Set debugging options to bash with NeR98 p. 221
# cmd line method | script method | Purpose
# csh -f           set              Fast start, ignore .cshrc and .login
# csh -n           set              Parse commands but do not execute them
# csh -v           set verbose      Echo commands before running them
# csh -x           set echo         Echo commands after command-line processing

# Usage: ~/sh/tst_csh.sh
# scp ~/sh/tst_csh.sh dataproc.ucar.edu:/fs/cgd/home0/zender/sh/tst_csh.sh

# Save stdout and stderr in C-shell scripts with, e.g., tst_csh.sh >&! foo

# NB: Never never never place comments (e.g., # foo) on same line as setenv
# It works in bash, but not in csh

echo
cat >! ~/foo.txt << EOF
Testing file catenation from ${0}
EOF

if ('0' == '1') then
    printf "ERROR: '0' == '1'\n"
else if ('1' == '1') then
    printf "'1' == '1'\n"
else
    printf "End of string comparison block\n"
endif

if (${HOST} =~ *.uci.edu) then
    printf "HOST = ${HOST} is in UCI domain\n"
else if (${HOST} =~ *.ucar.edu) then
    printf "HOST = ${HOST} is in UCAR domain\n"
else
    printf "ERROR: \${HOST} = ${HOST} has unknown domain\n"
endif # endif $HOST

if (${HOST} =~ esmf* || ${HOST} =~ b[bfs]*en) then
    printf "HOST = ${HOST} matches esmf* or b[bfs]*en\n"
else if (${HOST} =~ ipcc* || ${HOST} =~ pbs*) then
    printf "HOST = ${HOST} matches ipcc* or pbs*\n"
else
    printf "ERROR: \${HOST} = ${HOST} matches nothing\n"
endif # endif $HOST

if (${PVM_ARCH} =~ *AIX*) then
    printf "PVM_ARCH = ${PVM_ARCH} matches *AIX*\n"
else if (${PVM_ARCH} =~ *AMD64) then
    printf "PVM_ARCH = ${PVM_ARCH} matches *AMD64\n"
else if (${PVM_ARCH} =~ *LINUX) then
    printf "PVM_ARCH = ${PVM_ARCH} matches *LINUX\n"
else
    printf "ERROR: \${PVM_ARCH} = ${PVM_ARCH} does not match options\n"
    exit
endif # endif $PVM_ARCH

printf "0 = $0\n"
printf "{0} = ${0}\n"
printf "argv = $argv\n"
printf "argv[0] = ${argv[0]}\n"
printf "argv[0] = ${argv}[0]\n"
printf "argv[1] = ${argv}[1]\n"

set foo = /home/zender/cvs/cvs.txt
printf "${foo:t}\n" # ${0} is argv[0], :t extracts filename component only
