; $Id$

; Purpose: Plot AMSR-E soil moisture data

; Example usage:

; ncl < ~/ncl/amsre_gph.ncl
; gv ${DATA}/ps/amsre_200406_mdl_obs_dlt_vwc_sfc.ps &
; kpdf ${DATA}/ps/amsre_200406_mdl_obs_dlt_vwc_sfc.pdf &
; load "/home/zender/ncl/amsre_gph.ncl"
; tst("x11")

begin

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

fl_mdl = addfile("/data/zender/amsre/amsre_200406_1x1.nc","r")
fl_obs = addfile("/data/zender/amsre/amsre_200406_1x1.nc","r")
fl_dlt = addfile("/data/zender/amsre/amsre_200406_1x1.nc","r")
vwc_sfc_mdl = fl_mdl->vwc_sfc(0,:,:)
vwc_sfc_obs = fl_obs->vwc_sfc(0,:,:)
vwc_sfc_dlt = fl_dlt->vwc_sfc(0,:,:)

;wks_id = gsn_open_wks("x11","/data/zender/ps/amsre_200406_mdl_obs_dlt_vwc_sfc")
; Initialize graphics
if(.not.isvar("dev")) then
;	dev="x11"; "x11", "ps", or "ncgm"
	dev="ps"; "x11", "ps", or "ncgm"
end if
if(dev.eq."x11") then
	wks_id=create "wks" xWorkstationClass noparent
	end create
end if ; endif "x11"
if(dev.eq."ps") then
	wks_id=create "wks" psWorkstationClass defaultapp
		"wkPSFileName": "/data/zender/ps/amsre_200406_mdl_obs_dlt_vwc_sfc.ps"
        end create
end if ; endif "ps"
;setvalues wks_id
;	"wkColorMap": "gui_default"
;;	"wkColorMap": clr_map
;end setvalues

gsn_define_colormap(wks_id,"gui_default")
plot = new(2,graphic) 

res                     = True
res@gsnDraw             = False           ; don't draw
res@gsnFrame            = False           ; don't advance frame
res@cnInfoLabelOn       = False           ; turn off cn info label
res@cnFillOn            = True            ; turn on color
res@gsnSpreadColors     = True            ; spread out color table
res@gsnSpreadColorStart = 2 
res@gsnSpreadColorEnd   = 21 
;res@gsnAddCyclic        = True            ; fxm: longitude
res@lbLabelBarOn        = False           ; turn off individual cb's
res@cnLinesOn = False
res@mpProjection = "Mollweide"
res@mpPerimOn = True
res@mpEllipticalBoundary = True
res@cnInfoLabelFont = "simplex_roman"
res@gsnCenterStringParallelPosF = 0.4
res@cnLevelSelectionMode =  "ManualLevels"   
res@cnMinLevelValF       = 0.0
res@cnMaxLevelValF       = 0.2
res@cnLevelSpacingF     = 0.02

res@gsnStringFont = "simplex_roman"
res@gsnStringFont = "simplex_roman"
res@gsnLeftString = "vwc_sfc"
res@tiMainString = "AMSR-E June 2004"
res@gsnRightString = "Mean=fxm :V1:m:S:3:H10::N:m:S:-3:N:"
plot(0) = gsn_csm_contour_map(wks_id,vwc_sfc_obs,res) 
res@gsnStringFont = "simplex_roman"
res@gsnLeftString = "vwc_sfc"
res@tiMainString = "CLM/SNICAR June Avg."
res@gsnRightString = "Mean=fxm :V1:m:S:3:H10::N:m:S:-3:N:"
plot(1) = gsn_csm_contour_map(wks_id,vwc_sfc_mdl,res) 

pres1                  = True                   ; modify the panel plot
pres1@gsnPanelLabelBar = True                   ; add common colorbar
pres1@lbLabelBarOn = True
pres1@gsnFrame = False
pres1@lbAutoManage = True
pres1@vpHeightF = .55
pres1@lbOrientation = "vertical"
pres1@lbLabelFont = "simplex_roman"
pres1@gsnPanelBottom = 0.4
gsn_panel(wks_id,plot,(/2,1/),pres1)               ; now draw as one plot

res@gsnStringFont = "simplex_roman"
res@cnMinLevelValF = -0.2
res@cnMaxLevelValF = 0.2
res@cnLevelSpacingF = 0.05
res@gsnLeftString = "vwc_sfc"
res@tiMainString = "Difference"
res@gsnRightString = "Mean=fxm :V1:m:S:3:H10::N:m:S:-3:N:"
res@gsnStringFontHeight = .013
gsn_define_colormap(wks_id,"nrl_sirkes_3white")
plot2 = gsn_csm_contour_map(wks_id,vwc_sfc_mdl-vwc_sfc_obs,res)

pres2 = True
pres2@gsnPanelLabelBar = True
pres2@gsnPanelTop = 0.4
pres2@gsnPanelBottom = 0.1
pres2@lbOrientation = "vertical"
pres2@vpHeightF = .325
pres2@lbLabelFont = "simplex_roman"
pres2@lbLabelFontHeightF = .0115
pres2@gsnFrame = False
gsn_panel(wks_id,plot2,(/1,1/),pres2)
frame(wks_id)

end
