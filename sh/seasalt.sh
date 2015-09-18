# $Id$

# Purpose: seasalt-specific startup commands

# 20001002: seasalt does not pay attention to XF86Config swapkeys directory
xmodmap ~/.xmodmaprc

# Set the clock over the network
sudo rdate -s time.nist.gov
sudo /sbin/hwclock --systohc
# Override the timezone
export TZ='PST8PDT' # Colorado = MST7MDT, California = PST8PDT
