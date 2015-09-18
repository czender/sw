; $Id$

; Purpose: Template NCL for drawing multi-XY plots

; Example usage:

; ncl < /home/zender/ncl/odxc.ncl
; load "/home/zender/ncl/odxc.ncl"
; tst("x11")

begin ; begin tst()

; Initialize some fields
if(.not.isvar("dev")) then
	dev="x11"; "x11", "ps", or "ncgm"
;	dev="ps"; "x11", "ps", or "ncgm"
end if
mss_val=1.e36

; Get data
fl_in=addfile("/data/zender/tmp/951011_1200_arese_mdl_clr_aer.nc","r")
ttl_sng="Extinction Optical Depth 1200 TST October 11, 1995"
;fl_in=addfile("/data/zender/tmp/951015_1200_arese_mdl_clr_aer.nc","r")
;ttl_sng="Extinction Optical Depth 1200 TST October 15, 1995"
;fl_in=addfile("/data/zender/tmp/951030_1200_arese_mdl_cld_aer.nc","r")
;ttl_sng="Extinction Optical Depth 1200 TST October 30, 1995"

bnd_nbr=dimsizes(fl_in->bnd)
var_nm=(/\
	"odxc_spc_O2",\
	"odxc_spc_CO2",\
;	"odxc_spc_H2O",\
	"odxc_spc_O3",\
	"odxc_spc_O2O2",\
	"odxc_spc_O2N2",\
	"odxc_spc_NO2",\
	"odxc_spc_ice",\
	"odxc_spc_liq",\
	"odxc_spc_aer",\
	"odxc_spc_bga",\
	"odxc_spc_Ray",\
	"odxc_spc_H2OH2O",\
	"odxc_spc_ttl"\
	/)
var_sng=(/\
	"O:B:2:N:",\
	"CO:B:2:N:",\
;	"H:B:2:N:O",\
	"O:B:3:N:",\
	"O:B:2:N:-O:B:2:N:",\
	"O:B:2:N:-N:B:2:N:",\
	"NO:B:2:N:",\
	"Ice",\
	"Liq",\
	"Aer",\
	"BGA",\
	"Ray",\
	"(H:B:2:N:O):B:2:N:",\
	"Total"\
	/)
cst_nbr=dimsizes(var_nm)
odxc=new((/cst_nbr,bnd_nbr/),float,mss_val)

do cst_idx=0,cst_nbr-1
	odxc(cst_idx,:)=fl_in->$var_nm(cst_idx)$
end do; end loop over cst_idx

; Convert from SI to nice units
bnd_mcr=1.e6*fl_in->bnd

; Create application object. 
; Application name determines name of resource file ("odxc.res")
app_id = create "tst" appClass defaultapp
   "appDefaultParent": "True"
   "appUsrDir": "/home/zender/ncl"
end create

; NCAR Graphics commands
; Create color map.  
nbr_sys_clr=2
clr=new((/cst_nbr/),integer)
clr_map=new((/cst_nbr+nbr_sys_clr,3/),float)
; Color indices '0' and '1' are background and foreground colors respectively.
clr_map(0,:)=(/1.00,1.00,1.00/)
clr_map(1,:)=(/0.00,0.00,0.00/)
if(dev.eq."x11") then
	clr_map(0,:)=(/0.00,0.00,0.00/)
	clr_map(1,:)=(/1.00,1.00,1.00/)
end if ; endif "ps"
do cst_idx=0,cst_nbr-1
	clr_idx=cst_idx+nbr_sys_clr
	clr(cst_idx)=clr_idx
	red=cst_idx/(cst_nbr-1.)
	green=1.-red
	blue=.5
	clr_map(clr_idx,:)=(/red,green,blue/)
end do; end loop over cst_idx
;print(clr_map)
;print(clr)

if(dev.eq."x11") then
	wks_id=create "wks" xWorkstationClass noparent
	end create
end if ; endif "x11"
if(dev.eq."ps") then
	wks_id=create "wks" psWorkstationClass defaultapp
		"wkPSFileName": "/data2/zender/ps/odxc.ps"
        end create
end if ; endif "ps"
setvalues wks_id
	"wkColorMap": "default"
;	"wkColorMap": clr_map
end setvalues

; ID of data object becomes resource value for "xyCoordData".
odxc_id=create "xyData" coordArraysClass noparent
; CoordArrays resources
    "caXArray": bnd_mcr
    "caXMaxV": 2.
    "caXMinV": .2
    "caYArray": odxc
    "caYMaxV": .12
    "caYMinV": 0
    "caYMissingV": mss_val    
end create

vp=(/.15,.95,.3,.9/)
plot_id=create "xyPlot_01" xyPlotClass wks_id
; XYPlot resources
	"xyCoordData": odxc_id
; Legend resources
	"lgAutoManage": "False"
	"lgPerimOn": "True"
	"lgOrientation": "Horizontal"
; PlotManager resources
	"pmLegendDisplayMode": "always"
; Title resources
	"tiMainFontHeightF": 0.02
	"tiMainString": ttl_sng
	"tiXAxisFontHeightF": 0.02
	"tiXAxisString": "Wavelength (microns)"
	"tiYAxisFontHeightF": 0.02
	"tiYAxisString": "Optical Depth"
; TickMark resources
	"tmXBMinorPerMajor": 2
; ViewPort resources
	"vpHeightF" : vp(3)-vp(2)
	"vpWidthF" : vp(1)-vp(0)
	"vpXF" : vp(0)
	"vpYF" : vp(3)
end create

getvalues plot_id
	"xyCoordDataSpec": data_spec_id
end getvalues
setvalues data_spec_id
; XyDataSpec resources
	"xyExplicitLabels": var_sng
	"xyExplicitLegendLabels": var_sng
	"xyLabelMode": "NOLABELS"
	"xyLineColors":	clr
	"xyLineLabelFontHeightF": 0.02
	"xyMarkLineMode": "LINES"
	"xyMarkerColors": clr
;	"xyLineLabelFontColors": (/transparent,3,4,5/)
;	"xyDashPatterns": (/0,5,10,15/)
;	"xyLineThicknesses": (/1.0,2.0,3.0,4.0/)
end setvalues

draw(plot_id)
frame(wks_id)

end ; end tst()

