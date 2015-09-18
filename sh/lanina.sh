# $Id$

# Purpose: Attach lanina to network over PPP

# Usage: Script must be run as root

# lanina is PPP dialin host and gateway for home LAN
# Execute following commands on lanina
# Delete existing Ethernet interface(s), if any
hostname lanina.zender.org
hostname -i -v
/sbin/ifconfig eth0 down 
/sbin/ifconfig eth0:1 down 
# Ensure default route interface is not preset to eth0 before dialing
/sbin/route del default gw 192.168.1.1 metric 1
ppp-go
#ppp-go -NCAR
sleep 60
/sbin/ifconfig eth0 lanina.zender.org # Connect IP address with Ethernet interface
#/etc/ppp/chain start # Turn on IP masquerading

# Set the clock over the network
rdate -s time.nist.gov
/sbin/hwclock --systohc
# Override the timezone
#export TZ='PST8PDT' # Colorado = MST7MDT, California = PST8PDT
export TZ='MST7MDT' # Colorado = MST7MDT, California = PST8PDT
