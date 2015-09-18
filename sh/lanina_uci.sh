# $Id$

# Purpose: Attach lanina to fixed IP network at UCI

# Usage: Script must be run as root

# lanina is PPP dialin host and gateway for home LAN
# Execute following commands on lanina
# Delete existing Ethernet interface(s), if any
sudo route add default gw 128.200.14.1
