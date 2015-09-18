#!/bin/bash

# Purpose: contour plot of one variable, one to four netCDF files
#
# Example:
#   plotvar.sh < plotvar.in

# Todo:
# 1. Option for auto contour levels
# 2. Option to omit x and y axis labels depending on number of plots
# 3. Option of not printing avg
# 4. Documentation

#---------------------------------------------------------------
# parse args
#---------------------------------------------------------------
usage="usage: $0 [-xph] < plotvar.in"
function help ()
{
    echo " "
    echo $usage
    echo " "
    echo "  -x  send output to screen (x11)"
    echo "  -p  create postscript file"
    echo "  -h  display this help"
    echo " "
    echo "note: must specify contour level min, max and delta "
    echo " "
    exit 1
}
dev="eps"
while getopts "xph" opt; do
    case $opt in
        x  ) dev="x11" ;;
        p  ) dev="ps"  ;;
        h  ) help      ;;
        \? ) echo $usage
             exit 1    ;;
    esac
done
shift $(($OPTIND - 1))
if [ $# -ne 0 ]; then
    echo $usage
    exit 1
fi

#function crop
#{
#  line=$*
#  line=${line%%#*}; line=${line// /}; 
#  return $line; # can't return strings
#}

#---------------------------------------------------------------
# read input
#---------------------------------------------------------------
while read line; do
    if [ ${line:0:1} != '#' ]; then break; fi
done
declare -i nplots

           line=${line%%#*}; line=${line// /}; nplots=$line; echo $nplots
read line; line=${line%%#*};                   title=$line;  echo $title
read line; line=${line%%#*}; line=${line// /}; var=$line;    echo $var
read line; line=${line%%#*}; line=${line// /}; min=$line;    echo $min
read line; line=${line%%#*}; line=${line// /}; max=$line;    echo $max
read line; line=${line%%#*}; line=${line// /}; delta=$line;  echo $delta
read line; line=${line%%#*}; line=${line// /}; file1=$line;  echo $file1
read line; line=${line%%#*}; line=${line// /}; title1=$line; echo $title1
read line; line=${line%%#*}; line=${line// /}; file2=$line;  echo $file2
read line; line=${line%%#*}; line=${line// /}; title2=$line; echo $title2
read line; line=${line%%#*}; line=${line// /}; file3=$line;  echo $file3
read line; line=${line%%#*}; line=${line// /}; title3=$line; echo $title3
read line; line=${line%%#*}; line=${line// /}; file4=$line;  echo $file4
read line; line=${line%%#*}; line=${line// /}; title4=$line; echo $title4

if [ $nplots -lt 1 ]; then
    echo "$0 error: nplots < 1"
    exit
fi
if [ $nplots -gt 4 ]; then
    echo "$0 error: nplots > 4"
    exit
fi

#---------------------------------------------------------------
# variables
#---------------------------------------------------------------
psfile=${file1%.nc}.${var}.${nplots}plot

#---------------------------------------------------------------
# create ncl script
#---------------------------------------------------------------
rm -f tmp.ncl

cat > tmp.ncl <<EOF
load "/opt/local/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "/opt/local/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
begin

  wks                      = gsn_open_wks("${dev}","${psfile}")            
                             gsn_define_colormap(wks,"gui_default")
  res                      = True               ; plot mods desired
  res@cnFillOn             = True               ; turn on color fill
  res@gsnDraw              = False              ; do not draw yet
  res@gsnFrame             = False              ; do notadvance frame yet
  res@gsnSpreadColors      = True               ; use full colormap
  res@mpFillOn             = False              ; turn off continent gray
  res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
  res@cnMinLevelValF       = $min               ; set min contour level
  res@cnMaxLevelValF       = $max               ; set max contour level
  res@cnLevelSpacingF      = $delta             ; set contour spacing
  res@lbLabelBarOn         = False              ; turn off individual label bars
  nplots                   = $nplots
  plot                     = new($nplots,graphic)     ; create graphic array for panel plot

  ;#---------------------------------------------------------------
  ;# plot(0)
  ;#---------------------------------------------------------------
  a = addfile("${file1}","r")
  u = a->${var}
  n = dimsizes(dimsizes(u))
  if (n .eq. 2) then
    v = u
  end if
  if (n .eq. 3) then
    v = u(0,:,:)
  end if
  if (n .eq. 4) then
    v = u(0,0,:,:)
  end if
  res@tiMainString    = "${title1}"
  res@gsnCenterString = "avg = " + sprintf("%.3g", avg(v))
  plot(0)             = gsn_csm_contour_map(wks,v,res)    

  ;#---------------------------------------------------------------
  ;# plot(1)
  ;#---------------------------------------------------------------
  if ( nplots .gt. 1) then
    a = addfile("${file2}","r")
    u = a->${var}
    n = dimsizes(dimsizes(u))
    if (n .eq. 2) then
      v = u
    end if
    if (n .eq. 3) then
      v = u(0,:,:)
    end if
    if (n .eq. 4) then
      v = u(0,0,:,:)
    end if
    res@tiMainString    = "${title2}"
    res@gsnCenterString = "avg = " + sprintf("%.3g", avg(v))
    plot(1)             = gsn_csm_contour_map(wks,v,res)    
  end if

  ;#---------------------------------------------------------------
  ;# plot(2)
  ;#---------------------------------------------------------------
  if ( nplots .gt. 2) then
    a = addfile("${file3}","r")
    u = a->${var}
    n = dimsizes(dimsizes(u))
    if (n .eq. 2) then
      v = u
    end if
    if (n .eq. 3) then
      v = u(0,:,:)
    end if
    if (n .eq. 4) then
      v = u(0,0,:,:)
    end if
    res@tiMainString    = "${title3}"
    res@gsnCenterString = "avg = " + sprintf("%.3g", avg(v))
    plot(2)             = gsn_csm_contour_map(wks,v,res)    
  end if

  ;#---------------------------------------------------------------
  ;# plot(3)
  ;#---------------------------------------------------------------
  if ( nplots .gt. 3) then
    a = addfile("${file4}","r")
    u = a->${var}
    n = dimsizes(dimsizes(u))
    if (n .eq. 2) then
      v = u
    end if
    if (n .eq. 3) then
      v = u(0,:,:)
    end if
    if (n .eq. 4) then
      v = u(0,0,:,:)
    end if
    res@tiMainString    = "${title4}"
    res@gsnCenterString = "avg = " + sprintf("%.3g", avg(v))
    plot(3)             = gsn_csm_contour_map(wks,v,res)    
  end if

  ;#---------------------------------------------------------------
  ;# create panel plot
  ;#---------------------------------------------------------------
  resPanel                    = True               ; panel mods desired
  resPanel@txString           = "${title}"
  resPanel@gsnPanelLabelBar   = True               ; label bar on panel
  resPanel@lbLabelFontHeightF = .01                ; cb label font size
  resPanel@gsnMaximize        = True               ; blow up plot

  ;resPanel@gsnPanelBottom    = 0.05
  ;resPanel@pmLabelBarWidthF  = 10.6               ; control size of colorbar
  ;resPanel@pmLabelBarHeightF = 0.1
  ;resPanel@pmLabelBarOrthogonalPosF = -0.02       ; position cb wrt plot
  ;resPanel@lbAutoManage      = False              ; dont let plot manager contr


  if ( nplots .lt. 4) then
    gsn_panel(wks,plot,(/$nplots,1/),resPanel)
  else
    gsn_panel(wks,plot,(/2,2/),resPanel)
  end if

end
EOF

# run ncl with tmp.ncl
/opt/local/bin/ncl < tmp.ncl

# rm -f tmp.ncl
echo ">>> created file: ${psfile}.${dev}"
