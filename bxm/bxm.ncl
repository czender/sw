; $Id$

; Purpose: Plot results of box model run
; bxm.ncl is called by bxm_run.sh
; bxm.ncl creates bxm.ps from box model run (dead.nc)

; Usage:
; scp ~/bxm/bxm.ncl dust.ess.uci.edu:/var/www/html/dead/bxm.ncl

load "/usr/local/ncarg/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "/usr/local/ncarg/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin
  wks  = gsn_open_wks("ps","bxm")
  plot = new(5,graphic)                          ; create graphics array

  res                 = True                     ; Use plot options
  res@gsnDraw         = False                    ; do not draw plot (just yet)
  res@gsnFrame        = False                    ; do not advance frame
  res@vpHeightF       = 0.5                      ; h, w, loc of plot
  res@vpWidthF        = 2.5
  ;res@vpXF           = 0.13
  ;res@vpYF           = 0.85
  res@tiXAxisString   = ""                       ; null string till last plot
  ;res@tmXBLabelFontHeightF = 0.014	         ; resize tick labels
  res@tmYLLabelFontHeightF  = 0.04
  res@tiMainFontHeightF     = 0.05               ; resize title
  ;res@txFontHeightF        = 0.01               ; resize the text
  res@tmXBLabelsOn  = False
  res@tmXTOn        = False
  res@tmYROn        = False
  res@tmXTBorderOn  = False
  res@tmYRBorderOn  = False
  res@pmLegendDisplayMode    = "Always"
  res@pmLegendSide           = "Top"                        ; Change location of 
  res@pmLegendParallelPosF   = 1.3                          ; move units right
  res@pmLegendOrthogonalPosF = -1.2                         ; move units down
  res@pmLegendWidthF         = 0.6                          ; Change width and
  res@pmLegendHeightF        = 0.6                          ; height of legend.  
  res@lgPerimOn              = False                        ; turn off box around
  res@tiXAxisString = ""
  res@tiYAxisString = ""
  res@xyDashPatterns  = (/0,0,0,0,0/)

;----------------------------------------------------------
; plot 1: Wind Speed
;----------------------------------------------------------

  a = addfile("dead.nc","r")
  t = a->time
  time_nbr = dimsizes(t)

  a1 = addfile("tvbds.nc","r")
  y1 = new((/3, time_nbr/), float)
  tmp1 = a1->wnd_znl_mdp
  tmp2 = a1->wnd_mrd_mdp
  do n=0,time_nbr-1
	y1(2,n) = sqrt(tmp1(n)*tmp1(n) + tmp2(n)*tmp2(n))
  end do
  y1(1,:) = a->wnd_rfr_mbl(:,0)
  y1(0,:) = a->wnd_rfr_thr_slt(:,0)
 
  ;res@tiYAxisString = "wnd-mdp"
  res@tiMainString = "Wind Speed (m s:S:-1:N:)"
  res@xyExplicitLegendLabels = (/"10 m Slt. Thr.", "10 m Wind", "sed_txt_1 m Wind","",""/)
  res@xyLineColors  = (/3,2,1,5,4/)
  plot(0) = gsn_xy(wks,t,y1,res)

;----------------------------------------------------------
; plot 2: Wind Friction Speed
;----------------------------------------------------------

  y2 = new((/3, time_nbr/), float)
  y2(2,:) = a->wnd_frc_mbl(:,0)*100
  y2(1,:) = a->wnd_frc_slt(:,0)*100
  y2(0,:) = a->wnd_frc_thr_slt(:,0)*100
  res@tiMainString = "Wind Friction Speed (cm s:S:-1:N:)"
  res@xyExplicitLegendLabels = (/"u:B3:*,t:S1:s", "u:B1:*:S1:s", "u:B1:*","",""/)
  plot(1) = gsn_xy(wks,t,y2,res)

;----------------------------------------------------------
; plot 3: Dust Mobilization Flux
;----------------------------------------------------------

  y3 = new((/5, time_nbr/), float)
  y3(4,:) = a->flx_mss_vrt_dst_ttl(:,0)*1e9
  y3(3,:) = a->flx_mss_vrt_dst(:,0,0)*1e9
  y3(2,:) = a->flx_mss_vrt_dst(:,1,0)*1e9
  y3(1,:) = a->flx_mss_vrt_dst(:,2,0)*1e9
  y3(0,:) = a->flx_mss_vrt_dst(:,3,0)*1e9
  res@tiMainString = "Dust Mobilization Flux (:F33:m:F:g m:S:-2:N: s:S:-1:N:)"
  res@xyExplicitLegendLabels = (/"5-10 :F33:m:F:m","2.5-5 :F33:m:F:m","1-2.5 :F33:m:F:m","0-1 :F33:m:F:m","Total"/)
  res@xyLineColors  = (/5,4,3,2,1/)
  plot(2) = gsn_xy(wks,t,y3,res)

;----------------------------------------------------------
; plot 4: Dry Deposition Flux
;----------------------------------------------------------

  y4 = new((/5, time_nbr/), float)
  y4(4,:) = a->flx_mss_dry_sfc_ttl(:,0)*1e9
  y4(3,:) = a->flx_mss_dry_sfc(:,0,0)*1e9
  y4(2,:) = a->flx_mss_dry_sfc(:,1,0)*1e9
  y4(1,:) = a->flx_mss_dry_sfc(:,2,0)*1e9
  y4(0,:) = a->flx_mss_dry_sfc(:,3,0)*1e9
  res@tiMainString = "Dry Deposition Flux (:F33:m:F:g m:S:-2:N: s:S:-1:N:)"
  plot(3) = gsn_xy(wks,t,y4,res)

;----------------------------------------------------------
; plot 5: Column Mass Path
;----------------------------------------------------------

  y5 = new((/5, time_nbr/), float)
  y5(4,:) = a->mpc_dst_ttl(:,0)*1e3
  y5(3,:) = a->mpc_dst(:,0,0)*1e3
  y5(2,:) = a->mpc_dst(:,1,0)*1e3
  y5(1,:) = a->mpc_dst(:,2,0)*1e3
  y5(0,:) = a->mpc_dst(:,3,0)*1e3

  res@tmXBLabelsOn  = True
  res@tiMainString  = "Column Mass Path (g m:S:-2:N:)"
  res@tiXAxisString = "Time (days)"
  plot(4) = gsn_xy(wks,t,y5,res)


;************************************************************
; now create panel plot with its own resources
;************************************************************

  resPanel                  = True 
  ;resPanel@gsnMaximize     = True                 ; fill up the page
  ;resPanel@gsnPanelLeft    = 0.001
  ;resPanel@gsnPanelRight   = 0.999
  ;resPanel@gsnPanelBottom  = 0.001
  ;resPanel@gsnPanelTop     = 0.999
  ;resPanel@txString        = "DEAD Box Model"
  resPanel@gsnPanelLabelBar = True
  resPanel@gsnFrame         = False                  ; do not advance panel plot
  resPanel@gsnPanelBottom   = 0.05
  gsn_panel(wks,plot,(/5,1/),resPanel)

  txres               = True
  txres@txFontHeightF = 0.01
  gsn_text_ndc(wks,"sed_txt_2 DEAD Box Model (http:://dust.ess.uci.edu/dead) sed_txt_3",0.4,0.02,txres)
  frame(wks)

end
