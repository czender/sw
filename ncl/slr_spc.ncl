; $Id$

; Purpose: Template NCL for drawing multi-XY plots

; Example usage:

; ncl < ${HOME}/sw/ncl/slr_spc.ncl
; load "${HOME}/sw/ncl/slr_spc.ncl"
; spc("x11")

begin ; begin spc()

; Initialize some fields
if(.not.isvar("dev")) then
;	dev="x11"; "x11", "ps", or "ncgm"
	dev="ps"; "x11", "ps", or "ncgm"
end if

; Create application object. 
; Application name determines name of resource file ("slr_spc.res")
app_id = create "spc" appClass defaultapp
   "appDefaultParent": "True"
   "appUsrDir": "${HOME}/sw/ncl"
end create

mss_val=1.e36
spc_nm=(/\
	"Kur95_01wvn",\
	"Kur95_20wvn",\
;	"NeL84",\
	"ThD71",\
	"LaN68"\
	/)
cst_nbr=dimsizes(spc_nm)
fl_nm=new((/cst_nbr/),string)
data_id=new(cst_nbr,graphic)
;fl_in=new((/cst_nbr/),file)
do idx=0,cst_nbr-1
	fl_nm(idx)="$DATA/aca/spc_"+spc_nm(idx)+".nc"
;	fl_in(idx)=addfile(fl_nm(idx),"r")
	fl_in=addfile(fl_nm(idx),"r")
; ID of data object becomes resource value for "xyCoordData".
	data_id(idx)=create "xyData"+idx coordArraysClass defaultapp
	; CoordArrays resources
		"caXArray": 1.e6*fl_in->wvl
		"caXMaxV": 2.0
		"caXMinV": 0.
		"caYArray": 1.e-9*fl_in->flx_spc_wvl
;		"caYMaxV": 1.
		"caYMinV": 0.
		"caYMissingV": mss_val    
	end create ; end create data_id
end do; end loop over idx

; NCAR Graphics commands
; Create color map.  
nbr_sys_clr=2
clr=new((/cst_nbr/),integer)
clr_map=new((/cst_nbr+nbr_sys_clr,3/),float)
; Color indices '0' and '1' are background and foreground colors respectively.
clr_map(0,:)=(/1.00,1.00,1.00/)
clr_map(1,:)=(/0.00,0.00,0.00/)
do idx=0,cst_nbr-1
	clr_idx=idx+nbr_sys_clr
	clr(idx)=clr_idx
	red=idx/(cst_nbr-1.)
	green=1.-red
	blue=.5
	clr_map(clr_idx,:)=(/red,green,blue/)
end do; end loop over idx
;print(clr_map)
;print(clr)

if(dev.eq."x11") then
	wks_id=create "wks" xWorkstationClass noparent
		"wkColorMap": "default"
		"wkPause": "True"
		"wkXColorMode": "mixed"
	end create; wks_id
end if ; endif "x11"
if(dev.eq."ps") then
	wks_id=create "wks" psWorkstationClass defaultapp
		"wkPSFileName": "/data2/zender/ps/slr_spc.ps"
        end create; wks_id
end if ; endif "ps"
setvalues wks_id
	"wkColorMap": "default"
;	"wkColorMap": clr_map
end setvalues; wks_id

vp=(/.15,.95,.3,.9/)
plot_id=create "xyPlot_01" xyPlotClass wks_id
; XYPlot resources
	"xyCoordData": data_id
; Legend resources
	"lgAutoManage": "False"
	"lgPerimOn": "True"
; PlotManager resources
	"pmLegendDisplayMode": "always"
; Title resources
	"tiMainFont": 21
	"tiMainFontHeightF": 0.02
	"tiMainString": "TOA Solar Spectrum (1367 W m:S:-2:N:)"
	"tiXAxisFontHeightF": 0.02
	"tiXAxisString": "Wavelength (microns)"
	"tiYAxisFontHeightF": 0.02
	"tiYAxisString": "Spectral Flux (W m:S:-2:N: nm:S:-1:N:)"
; TickMark resources
	"tmXBMinorPerMajor": 3
	"tmYLMinorPerMajor": 2
; ViewPort resources
	"vpHeightF" : vp(3)-vp(2)
	"vpWidthF" : vp(1)-vp(0)
	"vpXF" : vp(0)
	"vpYF" : vp(3)
end create

getvalues plot_id
	"xyCoordDataSpec": data_spec_id
end getvalues
do idx=0,cst_nbr-1
setvalues data_spec_id(idx)
; XyDataSpec resources
	"xyExplicitLabels": spc_nm(idx)
	"xyExplicitLegendLabels": spc_nm(idx)
	"xyLabelMode": "NOLABELS"
	"xyLineColors":	clr(idx)
	"xyLineLabelFontHeightF": 0.02
	"xyMarkLineMode": "LINES"
	"xyMarkerColors": clr(idx)
;	"xyLineLabelFontColor": (/transparent,3,4,5/)
;	"xyDashPattern": (/0,5,10,15/)
;	"xyLineThickness": (/1.0,2.0,3.0,4.0/)
end setvalues; end setvalues data_spec_id(idx)
end do; end loop over idx

draw(plot_id)
frame(wks_id)

end ; end spc()

