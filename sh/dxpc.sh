#!/bin/csh -fv

# From dxpc man page:
#       Assume  that you're running a real X server on the console
#       of a local workstation called homepc, and that you want to
#       run  some  X  applications on a remote system called work-
#       server and have them display on the console of  the  local
#       system.
#       On homepc, run
#           $ export DISPLAY=unix:0
#           $ dxpc -s1

#       On workserver, run
#           $ dxpc -s1 homepc
#           $ export DISPLAY=unix:8
#           $ xterm&
#           $ xemacs&

# Following is what really works, from JMR:
# home: xhost dust.acd.ucar.edu
# ncar: setenv DISPLAY dust:8
# ncar: dxpc -s1 &
# home: dxpc dust.acd.ucar.edu &
# ncar: emacs &

# home: xhost sanitas.cgd.ucar.edu
# ncar: setenv DISPLAY sanitas:8
# ncar: dxpc -s1 &
# home: dxpc sanitas.cgd.ucar.edu &
# ncar: emacs &

# home: xhost dataproc.ucar.edu
# ncar: setenv DISPLAY dataproc:8
# ncar: dxpc -s1 &
# home: dxpc dataproc &
# ncar: emacs &

# Be careful not to let this shell script kill itself
#killproc dxpc
rsh odin "killproc dxpc"

xhost odin.cgd.ucar.edu
rsh odin "setenv DISPLAY odin:8"

rsh odin "dxpc -s1 &"
dxpc odin &

rsh odin "emacs &"

exit


