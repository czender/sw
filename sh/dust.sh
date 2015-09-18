# $Id$

# Purpose: seasalt-specific startup commands

# Set the clock over the network
sudo rdate -s time.nist.gov
sudo /sbin/hwclock --systohc
# Override the timezone
export TZ='PST8PDT' # Colorado = MST7MDT, California = PST8PDT
