#! /bin/sh
# Restart pcmcia
/etc/init.d/pcmcia restart
# Connect to Harry's home network
ifconfig eth0 down
ifconfig eth1 down
ifconfig eth1 up
route add default gw 192.168.1.1
ping -c 2 mangogw
ping -c 2 www.uci.edu
