#!/usr/bin/perl -w

# Purpose: Filter dangerous input
# bxm_mch.pl is called by bxm_run.sh to get hgt an mydate input
# Script accepts only numeric input to protect against malicious attacks using web forms

# Usage:
# sudo cp ~/bxm/bxm_mch.pl /tmp/bxm/bxm_mch.pl

# Permissions:
# Keep file out of HTML-accessible directory so crackers may not examine it for weaknesses
# sudo chmod 744 /tmp/bxm/bxm_mch.pl
# sudo chown apache /tmp/bxm/bxm_mch.pl
# sudo chgrp apache /tmp/bxm/bxm_mch.pl

use strict;
while(<>){
    $_ =~ s/:/::/g;
    $_ =~ s/\//-/g;
    if(/(.*)=(.*).\d\d\d/){print $2;}
    else{
	if(/(.*)=(.*)/){print $2;}
    } # end else
} # end while
