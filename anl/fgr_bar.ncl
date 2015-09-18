; $Id$ -*-f90-*-

; Purpose: Plot bar charts
; Based on: 
; http://www.ncl.ucar.edu/Applications/Images/unique_5_lg.png
; http://www.ncl.ucar.edu/Applications/Scripts/unique_5.ncl
; NB: Intended for use with defaults set in ~/.hluresfile

; Example usage:
; cd ~/anl;ncl < fgr_bar.ncl
; cd ~/anl;ncl 'dvc="x11"' fgr_bar.ncl
; cd ~/anl;ncl flg_clr=False fgr_bar.ncl
; cd ~/anl;ncl abc_nbr_plt=3 ord_max=0.18 'dvc="x11"' fgr_bar.ncl
; cd ~/anl;ncl abc_nbr_plt=3 ord_nbr_plt=3 ord_max=0.18 'dvc="x11"' fgr_bar.ncl

begin

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

; Initialize some fields
prg_nm="fgr_bar" ; [sng] Program name
fgr_stb="fgr_bar" ; [sng] Figure name stub
dbg_lvl=0
mss_val_flt=1.0e36

if(.not.isvar("dvc")) then dvc="eps" end if
if(.not.isvar("plt_ffc_frc_rsp")) then plt_ffc_frc_rsp=False end if
if(.not.isvar("plt_atm_snw")) then plt_atm_snw=False end if
if(.not.isvar("plt_tpt_rsp")) then plt_tpt_rsp=False end if
if(.not.isvar("fl_out_nm")) then fl_out_flg=False else fl_out_flg=True end if
if(.not.isvar("flg_clr")) then flg_clr=True end if
if(.not.isvar("flg_rng")) then flg_rng=False end if
if(.not.isvar("fld_nm")) then fld_nm="tpt_rsp" end if
if(.not.isvar("ord_max")) then ord_max=2.8 end if
if(.not.isvar("arc_scl")) then arc_scl=1.0 end if ; To plot Arctic response set Arctic amplification as (mean of 1997+2001 responses)/(global response) = (0.5*(1.61+0.5))/0.15 ~ 7.0
if(.not.isvar("nrm_flg")) then nrm_flg=False end if

if(plt_ffc_frc_rsp.eq.False.and.plt_atm_snw.eq.False.and.plt_tpt_rsp.eq.False) then print("ERROR: Set one of plt_ffc_frc_rsp, plt_atm_snw, and plt_tpt_rsp to True") end if

if (flg_clr) then
   fl_sfx_clr="_clr"
else ; !flg_clr
   fl_sfx_clr="_bw"
end if; !flg_clr

; Markers: http://www.ncl.ucar.edu/Document/Graphics/Images/markers.gif
idx_mrk_ast=03 ; Asterisk
idx_mrk_crc_crs=15 ; Circle with cross
idx_mrk_crc_dot=14 ; Circle with dot
idx_mrk_crc_sld=16 ; Circle solid
idx_mrk_crc=04 ; Circle 
idx_mrk_crs=05 ; Cross
idx_mrk_dmn=09 ; Diamond
idx_mrk_dot=01 ; Dot
idx_mrk_pls=02 ; Plus
idx_mrk_sqr=06 ; Square
idx_mrk_str=12 ; Star
idx_mrk_tri_dwn=08 ; Triangle down
idx_mrk_tri_up=07 ; Triangle up

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Figure fgr_bar
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if(.not.isvar("ttl_sng")) then ttl_sng="Warming by Dirty Snow" end if
if(.not.isvar("x_ttl_sng")) then x_ttl_sng="Climate" end if
if(.not.isvar("y_ttl_sng")) then y_ttl_sng="Warming ~F33~D~F10~T~F21~ [~S~o~N~C]" end if

abc_nbr=4 ; [nbr] Number of timeperiods
; Actual number of abscissae to plot: abc_nbr_plt <= abc_nbr
if(.not.isvar("abc_nbr_plt")) then abc_nbr_plt=abc_nbr end if
abc=ispan(1,abc_nbr_plt,1) ; srt,end,srd
abc_dlt=0.15 ; Bar width (offset between successive x-locations)
abc_max=max(abc) ; First abscissa numeric value
abc_min=min(abc) ; Last abscissa numeric value
ord_nbr=3 ; [nbr] Number of aerosols
data=new((/abc_nbr,ord_nbr/),float,mss_val_flt)
; Actual number of ordinates to plot: ord_nbr_plt <= ord_nbr
if(.not.isvar("ord_nbr_plt")) then ord_nbr_plt=ord_nbr end if

; Order of climates along horizontal axis
idx_pd=1 ; [idx] PD index
idx_pi=0 ; [idx] PI index
idx_2050a2=2 ; [idx] 2050 A2 index
idx_lgm=3 ; [idx] LGM index (abc_nbr_plt discards abscissae with highest indices)

; Forcing data
frc_dst=new((/abc_nbr/),float,mss_val_flt) ; [W m-2] Snowpack heating by dust impurity
frc_bc=new((/abc_nbr/),float,mss_val_flt) ; [W m-2] Snowpack heating by soot impurity
frc_aer=new((/abc_nbr/),float,mss_val_flt) ; [W m-2] Snowpack heating by aerosol impurity
frc_dst(idx_lgm)=0.2015 ; [W m-2] Snowpack heating by dust impurity
frc_dst(idx_pi)=0.0363 ; [W m-2] Snowpack heating by dust impurity
frc_dst(idx_pd)=0.0165 ; [W m-2] Snowpack heating by dust impurity
frc_dst(idx_2050a2)=0.0105 ; [W m-2] Snowpack heating by dust impurity
frc_bc(idx_lgm)=0.0152 ; [W m-2] Snowpack heating by soot impurity
frc_bc(idx_pi)=0.0099 ; [W m-2] Snowpack heating by soot impurity
frc_bc(idx_pd)=0.0347 ; [W m-2] Snowpack heating by soot impurity
frc_bc(idx_2050a2)=0.0220 ; [W m-2] Snowpack heating by soot impurity
frc_aer(idx_lgm)=0.2202 ; [W m-2] Snowpack heating by aerosol impurity
frc_aer(idx_pi)=0.0482 ; [W m-2] Snowpack heating by aerosol impurity
frc_aer(idx_pd)=0.0578 ; [W m-2] Snowpack heating by aerosol impurity
frc_aer(idx_2050a2)=0.0357 ; [W m-2] Snowpack heating by aerosol impurity

; Temperature Response data
tpt_rsp_dst=new((/abc_nbr/),float,mss_val_flt) ; [K] Temperature response to dust impurity
tpt_rsp_bc=new((/abc_nbr/),float,mss_val_flt) ; [K] Temperature response to soot impurity
tpt_rsp_aer=new((/abc_nbr/),float,mss_val_flt) ; [K] Temperature response to aerosol impurity

ffc_aer=new((/abc_nbr/),float,mss_val_flt) ; [frc] Efficacy of aerosol impurity
ffc_aer(idx_lgm)=17.39 ; [K] Efficacy of aerosol impurity
ffc_aer(idx_pi)=1.112 ; [K] Efficacy of aerosol impurity
ffc_aer(idx_pd)=3.871 ; [K] Efficacy of aerosol impurity
ffc_aer(idx_2050a2)=5.058 ; [K] Efficacy of aerosol impurity

tpt_rsp_aer=new((/abc_nbr/),float,mss_val_flt) ; [K] Temperature response to aerosol impurity
tpt_rsp_aer(idx_lgm)=2.750 ; [K] Temperature response to aerosol impurity
tpt_rsp_aer(idx_pi)=0.0370 ; [K] Temperature response to aerosol impurity
tpt_rsp_aer(idx_pd)=0.1544 ; [K] Temperature response to aerosol impurity
tpt_rsp_aer(idx_2050a2)=0.1237 ; [K] Temperature response to aerosol impurity
; Correct for Arctic response
if (arc_scl.ne.1.0) then tpt_rsp_aer=arc_scl*tpt_rsp_aer end if
if (arc_scl.ne.1.0) then ttl_sng="Arctic "+ttl_sng end if

; Partition total into component aerosol responses
tpt_rsp_aer_rng=0.0*tpt_rsp_aer
tpt_rsp_dst=(frc_dst/frc_aer)*tpt_rsp_aer
tpt_rsp_dst_rng=0.0*tpt_rsp_dst
tpt_rsp_bc=(frc_bc/frc_aer)*tpt_rsp_aer
tpt_rsp_bc_rng=0.0*tpt_rsp_bc

; Legend places lines in reverse order
idx_aer=0
idx_bc=1
idx_dst=2

; Legend places lines in reverse order
idx_ffc=0
idx_frc=1
idx_rsp=2

; Legend places lines in reverse order
idx_atm=0
idx_snw=1

if(plt_ffc_frc_rsp) then
   data(:,idx_ffc)=ffc_aer
   data(:,idx_frc)=frc_aer
   data(:,idx_rsp)=tpt_rsp_aer
end if ; !plt_ffc_frc_rsp
if(plt_atm_snw) then
   if(fld_nm.eq."frc_toa_dst") then
      data(:,idx_atm)=frc_toa_dst_atm
      data(:,idx_snw)=frc_toa_dst_snw
   end if 
   if(fld_nm.eq."tpt_rsp_dst") then
      data(:,idx_atm)=tpt_rsp_dst_atm
      data(:,idx_snw)=tpt_rsp_dst_snw
   end if 
   if(fld_nm.eq."ffc_dst") then
      data(:,idx_atm)=ffc_dst_atm
      data(:,idx_snw)=ffc_dst_snw
   end if 
   if(fld_nm.eq."frc_ffc_dst") then
      data(:,idx_atm)=frc_ffc_dst_atm
      data(:,idx_snw)=frc_ffc_dst_snw
   end if
end if ; !plt_atm_snw
if(plt_tpt_rsp) then ;
   data(:,idx_dst)=tpt_rsp_dst
   data(:,idx_bc)=tpt_rsp_bc
   data(:,idx_aer)=tpt_rsp_aer
end if ; ! Default plot tpt_rsp

if(nrm_flg) then 
; Normalize forcing and response by same amount to be visible
   nrm_fct_frc=data(idx_pi,idx_ffc)/data(idx_pi,idx_frc)
   nrm_fct=data(idx_pi,idx_ffc)/data(idx_pi,idx_frc)
   data(:,idx_frc)=nrm_fct*data(:,idx_frc)
   data(:,idx_rsp)=nrm_fct*data(:,idx_rsp)
end if ; !nrm_flg

abc_sng=new((/abc_nbr/),string)
abc_sng(idx_lgm)="LGM"
abc_sng(idx_pi)="Pre-industrial"
abc_sng(idx_pd)="Present"
abc_sng(idx_2050a2)="2050A2"

dat_sng=new((/ord_nbr/),string)
if (plt_ffc_frc_rsp) then
   dat_sng(idx_ffc)="Efficacy"
   dat_sng(idx_frc)="Forcing"
   dat_sng(idx_rsp)="Response"
else ; !plt_ffc_frc_rsp
   dat_sng(idx_dst)="Dust"
   dat_sng(idx_bc)="BC"
   dat_sng(idx_aer)="Total"
end if ; !plt_ffc_frc_rsp

clr_sng=new((/ord_nbr/),string)
clr_sng(:)="black"
if (flg_clr) then
   clr_sng(idx_aer)="black"
   clr_sng(idx_dst)="blue"
   clr_sng(idx_bc)="green"
end if; !flg_clr

; Begin plot information
if (fl_out_flg) then
   fl_out=fl_out_nm
else 
   drc_out=getenv("HOME")+"/ppr_ZFM07/fgr"
   fl_out=drc_out+"/"+fgr_stb+"_"+fld_nm+fl_sfx_clr
end if
print("Output will go to "+fl_out+"."+dvc)
  
rsr=True
;rsr@vpWidthF=0.7
;rsr@vpHeightF=0.5
;rsr@vpXF=.15
rsr@trXMinF=0.4
rsr@trXMaxF=abc_max+0.6
rsr@trYMinF=0.0
rsr@trYMaxF=ord_max ; Set this or else axes overprint in loop
rsr@gsnDraw=True
rsr@gsnFrame=False
rsr@gsnXYBarChart=True
rsr@gsnXYBarChartBarWidth=abc_dlt ; Change bar widths
rsr@tmXBMode="Explicit" ; Explicit labels
rsr@tmXBValues=ispan(1,abc_nbr_plt,1) ; srt,end,srd
;rsr@tmXBValues=fspan(1,abc_nbr_plt,abc_nbr_plt) ; srt,end,nbr

rsr@tmXBLabels=abc_sng(0:abc_nbr_plt-1)
rsr@tmXBLabelFontHeightF=0.0205
rsr@tmXTLabelFontHeightF=0.0205
rsr@tmYLLabelFontHeightF=0.0225
rsr@tiMainFontHeightF=0.025
rsr@tiMainFont="helvetica"
rsr@tiMainString=ttl_sng
rsr@gsnRightString=""
rsr@tiYAxisString=y_ttl_sng

; Open workstation and plot data
wks=gsn_open_wks(dvc,fl_out)

; Draw timeseries separately by offsetting initial point
ord_gph=new((/ord_nbr_plt/),graphic)
abc_fst=0.775 ; Left-most bar x-location (Assumes x-axis is [1..abc_nbr])
do ord_idx=0,ord_nbr_plt-1
   rsr@gsnXYBarChartColors=(/clr_sng(ord_idx)/)	
   abc_srt=abc_fst+ord_idx*abc_dlt ; Initial bar x-location
   abc_end=abc_max-1+abc_srt ; Final bar x-location
   ord_gph=gsn_csm_xy(wks,fspan(abc_srt,abc_end,abc_nbr_plt),data(0:abc_nbr_plt-1,ord_idx),rsr)
end do; end loop over ord

; Labelbar inside bar chart
rsr_lbl=True ; Labelbar resources
rsr_lbl@vpWidthF=0.3 ; Labelbar width
rsr_lbl@vpHeightF=0.10 ; Labelbar height
rsr_lbl@lbBoxMajorExtentF =0.36 ; Space between color boxes
rsr_lbl@lbMonoFillPattern =True ; Solid fill pattern
rsr_lbl@lbLabelFontHeightF=0.035 ; Font height
rsr_lbl@lbLabelJust="CenterLeft" ; Left justify labels
rsr_lbl@lbPerimOn=False
rsr_lbl@lgPerimColor="white"
rsr_lbl@lbFillColors=clr_sng(0:ord_nbr_plt-1:-1) ; Make left-most into top-most label
if(.not.isvar("lbl_lft")) then lbl_lft=0.17 end if
if(.not.isvar("lbl_top")) then lbl_top=0.75 end if
gsn_labelbar_ndc(wks,ord_nbr_plt,dat_sng(0:ord_nbr_plt-1:-1),lbl_lft,lbl_top,rsr_lbl)

frame(wks)      

delete(data)
delete(rsr)
delete(rsr_lbl)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Figure fgr_bar
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
end
