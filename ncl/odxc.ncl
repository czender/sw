; $Id$

; Purpose: Plot extinction optical depth from swnb model

; Example usage:

; ncl < ${HOME}/sw/ncl/odxc.ncl
; load "${HOME}/sw/ncl/odxc.ncl"

begin ; begin odxc()

; These files are loaded by default in NCL V6.2.0 and newer
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

; Initialize some fields
if(.not.isvar("dvc")) then
;	dvc="x11"; "x11", "ps", or "ncgm"
;	dvc="ps"; "x11", "ps", or "ncgm"
	dvc="eps"; "x11", "ps", or "ncgm"
end if
mss_val=1.e36

; Get data
;ttl_sng="Extinction Optical Depth for Mid-Latitude Summer"
;fl_in=addfile("$DATA/tmp/swnb_mls_O2O2_O2N2.nc","r")
;ttl_sng="Extinction Optical Depth 1200 TST October 11, 1995"
;fl_in=addfile("$DATA/arese/mdl/951011_1200_arese_mdl_clr_aer.nc","r")
;ttl_sng="Sub-Arctic Summer"
;fl_in=addfile("$DATA/tmp/swnb_sas_O2O2_O2N2.nc","r")
;ttl_sng="Sub-Arctic Winter"
;fl_in=addfile("$DATA/tmp/swnb_saw_O2O2_O2N2.nc","r")
;ttl_sng="Extinction Optical Depth 1200 TST October 30, 1995"
;fl_in=addfile("$DATA/arese/mdl/951030_1200_arese_mdl_cld_aer.nc","r")
ttl_sng="Absorption Optical Depth" ; Use in OdX paper
fl_in=addfile("$DATA/arese/mdl/951015_1200_arese_mdl_clr_aer.nc","r")
;fl_in=addfile("$DATA/arese/mdl/951015_1200_arese_mdl_clr_aer_nbm.nc","r")
bnd_nbr=dimsizes(fl_in->bnd)
var_nm=(/\
	"odxc_spc_H2O",\
	"odxc_spc_O3",\
	"odxc_spc_O2",\
	"odxc_spc_CO2",\
	"odxc_spc_NO2",\
	"odxc_spc_O2O2",\
	"odxc_spc_O2N2"\
;	"odac_spc_aer"\
;	"odac_spc_bga"\
;	"odxc_spc_ice",\
;	"odxc_spc_liq",\
;	"odxc_spc_Ray"\
;	"odxc_spc_H2OH2O"\
;	"odxc_spc_ttl"\
	/)
var_sng=(/\
	"H~B~2~N~O",\
	"O~B~3~N~",\
	"O~B~2~N~",\
	"CO~B~2~N~",\
	"NO~B~2~N~",\
	"O~B~2~N~O~B~2~N~",\
	"O~B~2~N~N~B~2~N~"\
;	"Dust"\
;	"SO4"\
;	"Ice",\
;	"Liq",\
;	"Ray"\
;	"(H~B~2~N~O)~B~2~N~"\
;	"Total"\
	/)
cst_nbr=dimsizes(var_nm)
odxc=new((/cst_nbr,bnd_nbr/),float,mss_val)

do cst_idx=0,cst_nbr-1
	odxc(cst_idx,:)=fl_in->$var_nm(cst_idx)$
end do; end loop over cst_idx

; Convert from SI to nice units
bnd_mcr=1.0e6*fl_in->bnd

; Create application object 
; Application name determines name of resource file ("odxc.res")
app_id = create "tst" appClass defaultapp
   "appDefaultParent": "True"
   "appUsrDir": "${HOME}/sw/ncl"
end create

; NCAR Graphics commands
; Create color map  
nbr_sys_clr=2
clr=new((/cst_nbr/),integer)
clr_map=new((/cst_nbr+nbr_sys_clr,3/),float)
; Color indices '0' and '1' are background and foreground colors respectively.
clr_map(0,:)=(/1.00,1.00,1.00/)
clr_map(1,:)=(/0.00,0.00,0.00/)
if(dvc.eq."x11") then
	clr_map(0,:)=(/0.00,0.00,0.00/)
	clr_map(1,:)=(/1.00,1.00,1.00/)
end if ; endif "ps"
idx_yll=0
idx_blu=1
idx_blu_lgt=2
idx_gry=3
idx_grn_lgt=4
idx_red=5
idx_grn=6
do cst_idx=0,cst_nbr-1
	clr_idx=cst_idx+nbr_sys_clr
	clr(cst_idx)=clr_idx
	red=cst_idx/(cst_nbr-1.0)
	green=1.0-red
	blue=0.5
	clr_map(clr_idx,:)=(/red,green,blue/)
	if(cst_idx.eq.idx_yll) then clr_map(clr_idx,:)=(/1.0,1.0,0.0/) end if
	if(cst_idx.eq.idx_red) then clr_map(clr_idx,:)=(/1.0,0.0,0.0/) end if
	if(cst_idx.eq.idx_grn) then clr_map(clr_idx,:)=(/0.0,1.0,0.0/) end if
	if(cst_idx.eq.idx_blu) then clr_map(clr_idx,:)=(/0.0,0.0,1.0/) end if
	if(cst_idx.eq.idx_blu_lgt) then clr_map(clr_idx,:)=(/0.0,1.0,1.0/) end if
	if(cst_idx.eq.idx_grn_lgt) then clr_map(clr_idx,:)=(/0.5,1.0,0.0/) end if
	if(cst_idx.eq.idx_gry) then clr_map(clr_idx,:)=(/0.3,0.3,0.3/) end if
end do; end loop over cst_idx
;print(clr_map)
;print(clr)

fl_out="$DATA/ps/odxc"
wks=gsn_open_wks(dvc,fl_out)
setvalues wks
;	"wkColorMap": "default"
	"wkColorMap": clr_map
end setvalues

; ID of data object becomes resource value for "xyCoordData".
odxc_id=create "xyData" coordArraysClass noparent
; CoordArrays resources
    "caXArray": bnd_mcr
    "caXMaxV": 1.4
    "caXMinV": 0.3
    "caYArray": odxc
    "caYMaxV": 0.04
    "caYMinV": 0.0
    "caYMissingV": mss_val    
end create

vp=(/0.1,0.975,0.3,0.9/)
plot_id=create "xyPlot_01" xyPlotClass wks
; XYPlot resources
	"xyCoordData": odxc_id
; Legend resources
	"lgAutoManage": "False"
	"lgPerimOn": "False"
	"lgOrientation": "Horizontal"
; PlotManager resources
	"pmLegendDisplayMode": "always"
; Title resources
	"tiMainFontHeightF": 0.03
	"tiMainString": ttl_sng
	"tiXAxisFontHeightF": 0.03
;	"tiXAxisString": "Wavelength (~F33~m~F21~m)"
	"tiXAxisString": "Wavelength (~0714~m)"
	"tiYAxisFontHeightF": 0.03
	"tiYAxisString": "Optical Depth"
; TickMark resources
;	"tmXBMinorPerMajor": 5
; Tranformation  resources
;	"trYLog": "True"
;	"trYMinF": 0.01
;	"trYMaxF": 0.2
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
	"xyLineThicknesses": (/3.0,3.0,3.0,3.0,3.0,3.0,3.0/)
end setvalues

draw(plot_id)
frame(wks)

end ; end odxc()

