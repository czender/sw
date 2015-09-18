# $Id$

# Purpose: Attach lanina to UIO network as fixed IP fullstorm.uio.no

# Usage: Script must be run as root

hostname fullstorm.uio.no
hostname -i -v
# Delete existing Ethernet interface(s), if any
/sbin/ifconfig eth0 down
/sbin/ifconfig eth0:1 down 
/sbin/ifconfig eth0 inet 129.240.20.81 # Connect IP address with Ethernet interface
/sbin/ifconfig eth0 netmask 255.255.254.0 broadcast 129.240.21.255
# Route to gateway host all datagrams bound for hosts outside intranet
# fxm: "metric 1" argument appears necessary
/sbin/route add default gw 129.240.20.1 metric 1
