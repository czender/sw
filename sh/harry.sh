# $Id$

# Purpose: Attach lanina to network over Harry's @home server

# Usage: Script must be run as root
# Must have @home gateway specified in resolv.conf

# 20011203: Following lines to restart network necessary after RH Linux 7.2 upgrade
/etc/rc.d/init.d/network restart
/sbin/ifconfig eth0 down 
/etc/rc.d/init.d/pcmcia restart
/sbin/ifconfig eth0 up
# It should be possible to remove the preceding two lines by modifying the
# appropriate files in /etc (which were working fine with RHL 7.0)

# Delete existing Ethernet interface(s), if any
/sbin/ifconfig eth0 down 
sleep 1
/sbin/ifconfig eth0 add address 192.168.1.73 add netmask 255.255.255.0 

# Remove static route that may prevent accessing subnet 14 from home
# This works interactively (with sudo) but not in script for some reason
# sudo /sbin/route add default gw 192.168.1.1;sudo /sbin/route del -net 128.200.14.0 netmask 255.255.255.0 eth0;sudo /sbin/route del default gw 128.200.14.1
/sbin/route add default gw 192.168.1.1 
/sbin/route del -net 128.200.14.0 netmask 255.255.255.0 eth0
/sbin/route del default gw 128.200.14.1

# Deprecate this section once ntpd works
# Only synchronize with machines guaranteed to be up, or boots may fail when timekeeper is down
# ntp.ucsd.edu = 132.239.254.49 is a nearby Level 2 server
tm_srv='ntp.ucsd.edu'
echo "Synchronizing system time and hardware clock time with ${tm_srv}..."
ntpdate -s ${tm_srv}
/sbin/hwclock --systohc

# Override the timezone
export TZ='PST8PDT' # Colorado = MST7MDT, California = PST8PDT
