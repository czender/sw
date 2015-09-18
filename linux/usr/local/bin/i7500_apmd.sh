#!/bin/bash
# Purpose: apmcontinue script is executed when i7500 resumes
# This script is called by /etc/sysconfig/apmscript when laptop resumes
# Source: Suggested by dmg@csg.uwaterloo.ca
# Fix clock which will be wrong after suspends
# Fix slow keyboard rate after suspends
(
/sbin/hwclock --hctosys
/sbin/kbdrate -r 30 -d 250
\echo -n >/dev/tty0
)
