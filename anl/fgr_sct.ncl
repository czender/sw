; $Id$ -*-f90-*-

; Purpose: Scatterplots
; Based on: 
; http://www.ncl.ucar.edu/Applications/Images/unique_5_lg.png
; http://www.ncl.ucar.edu/Applications/Scripts/unique_5.ncl
; clay.ess.uci.edu:~qhan/Aldata/ncl_scripts/scatter_compare.ncl
; clay.ess.uci.edu:~qhan/ecosys/scripts/panel_model_obs_surAl.ncl

; NB: Intended for use with defaults set in ~/.hluresfile

; Example usage:
; cd ~/sw/anl;ncl < fgr_sct.ncl
; cd ~/sw/anl;ncl 'dvc="x11"' fgr_sct.ncl
; cd ~/sw/anl;ncl flg_ZGD09=True fgr_sct.ncl
; cd ~/sw/anl;ncl flg_crc=True flg_ZGD09=True fgr_sct.ncl
; cd ~/sw/anl;ncl flg_crc=True flg_ZGD09=True 'dvc="x11"' fgr_sct.ncl
; cd ~/sw/anl;ncl flg_crc=True flg_ZGD09=True flg_ssa=True 'dvc="x11"' fgr_sct.ncl
; cd ~/sw/anl;ncl flg_crc=True flg_ZGD09=True flg_shp_wvl=False abc_min=0.75 ord_min=0.75 abc_max=0.90 ord_max=0.90 'dvc="x11"' fgr_sct.ncl
; cd ~/sw/anl;ncl flg_crc=True flg_ZGD09=True 'fl_in_nm="/data/zender/icr/swnb_lgge_bee_nsm_mpr.nc"' 'dvc="x11"' fgr_sct.ncl

begin

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

; Personal headers
if(.not.isvar("csz_ncl")) then load "${HOME}/sw/ncl/csz.ncl" end if

; Initialize some fields
prg_nm="fgr_sct" ; [sng] Program name
fgr_stb="fgr_sct" ; [sng] Figure name stub
dbg_lvl=0
campaign_defined=False

if(.not.isvar("abc_max")) then abc_max=0.9 end if
if(.not.isvar("abc_min")) then abc_min=0.1 end if
if(.not.isvar("dvc")) then dvc="x11" end if
if(.not.isvar("fl_in_nm")) then fl_in_flg=False else fl_in_flg=True end if 
if(.not.isvar("fl_out_nm")) then fl_out_flg=False else fl_out_flg=True end if
if(.not.isvar("fld_nm")) then fld_nm="alb_spc_snw" end if
if(.not.isvar("flg_ZGD09")) then flg_ZGD09=False end if ; Plot ZGD09 mdl/msr
if(.not.isvar("flg_bee")) then flg_bee=True end if ; Best estimate of everything
if(.not.isvar("flg_clr")) then flg_clr=True end if
if(.not.isvar("flg_crc")) then flg_crc=True end if ; Apply geometry corrections
if(.not.isvar("flg_log_x")) then flg_log_x=False end if ; Logarithmic x-axis
if(.not.isvar("flg_log_y")) then flg_log_y=False end if ; Logarithmic y-axis
if(.not.isvar("flg_mop")) then flg_mop=False end if ; Apply MOP range
if(.not.isvar("flg_nir")) then flg_nir=True end if
if(.not.isvar("flg_nmn")) then flg_nmn=False end if ; Plot nominal snow
if(.not.isvar("flg_phb_he")) then flg_phb_he=False end if ; Plot phb_he results
if(.not.isvar("flg_phb_low")) then flg_phb_low=False end if ; Plot phb_low results
if(.not.isvar("flg_rms")) then flg_rms=True end if ; Print RMS statistic on graph
if(.not.isvar("flg_rng")) then flg_rng=False end if
if(.not.isvar("flg_shp_wvl")) then flg_shp_wvl=True end if ; Marker shape denotes wavelength
if(.not.isvar("flg_ssa")) then flg_ssa_dns=False end if ; Plot experiments predicted with measured SSA and density
if(.not.isvar("flg_vsb")) then flg_vsb=True end if
if(.not.isvar("ord_max")) then ord_max=0.9 end if
if(.not.isvar("ord_min")) then ord_min=0.1 end if
if(.not.isvar("plt_snw_lgge")) then plt_snw_lgge=True end if
if(.not.isvar("wvl_max")) then wvl_max=0.630e-6 end if
if(.not.isvar("wvl_min")) then wvl_min=0.640e-6 end if

if(plt_snw_lgge.eq.False) then print("ERROR: Set one of plt_snw_lgge to True") end if

if (flg_clr) then
   fl_sfx_clr="_clr"
else ; !flg_clr
   fl_sfx_clr="_bw"
end if; !flg_clr

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Figure fgr_sct
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if(.not.isvar("ttl_sng")) then ttl_sng="" end if
if(.not.isvar("x_ttl_sng")) then x_ttl_sng="Measured Albedo" end if
if(.not.isvar("y_ttl_sng")) then y_ttl_sng="Modeled Albedo" end if

if(.not.isvar("ord_nbr")) then ord_nbr=7 end if ; [nbr] Groups read

; Handle pre-defined dataset groups ("campaigns")
if (flg_ZGD09) then
   if(.not.isvar("snw_rfl_ncl")) then load "$HOME/anl/snw_rfl.ncl" end if
;   print(dns_sal_0635nm_20071206e)
   fl_nbr=1
end if ; end if flg_ZGD09

; Assemble file names for non-campaigns
if (.not.campaign_defined) then
   ord_nbr=7
   fl_nbr=0   
   data_xpt=fspan(abc_min,abc_max,ord_nbr)
   data_mdl=fspan(abc_min,abc_max,ord_nbr)
end if ; endif not campaign

; Input file for model data
if(isvar("fl_in")) then delete(fl_in) end if 
fl_in=new(1,string)
fl_in(0)=getenv("DATA")+"/icr/"+"swnb_lgge_bee_nsm_mpr.nc"
if (fl_in_flg) then fl_in(0)=fl_in_nm end if; !fl_in_flg
if (fl_nbr.gt.0) then print("Input data file "+idx_zero+" = "+fl_in(idx_zero)) end if

; Legend places data bottom-to-top (small-to-large index)
; ord_nbr_plt discards ordinates with smallest indices
idx_top  =ord_nbr-1 ; [idx] top index
idx_topm1=ord_nbr-2 ; [idx] 
idx_topm2=ord_nbr-3 ; [idx] 
idx_topm3=ord_nbr-4 ; [idx] 
idx_topm4=ord_nbr-5 ; [idx] 
idx_topm5=ord_nbr-6 ; [idx] 
idx_topm6=ord_nbr-7 ; [idx] 

clr_sng=new((/ord_nbr/),string)
clr_sng(:)="black"
if (flg_clr) then
   clr_sng(idx_top)="black"
   clr_sng(idx_topm1)="blue"
   clr_sng(idx_topm2)="green"
   clr_sng(idx_topm3)="orange"
   clr_sng(idx_topm4)="black"
   clr_sng(idx_topm5)="blue"
   clr_sng(idx_topm6)="green"
end if; !flg_clr

; Line properties
ln_sty=0*new(ord_nbr,integer,"No_FillValue")
ln_thk=2+0*new(ord_nbr,integer,"No_FillValue")

; Marker properties
mrk_clr=clr_sng ; Marker color
mrk_ln_md=new((/ord_nbr/),string)
mrk_ln_md(:)="Markers" ; Line/marker combinations
mrk_sz=new(ord_nbr,float)
mrk_sz(:)=0.02
mrk_typ=new(ord_nbr,integer)
mrk_typ(idx_top)=idx_mrk_crc_sld
mrk_typ(idx_topm1)=idx_mrk_crs
mrk_typ(idx_topm2)=idx_mrk_crc
mrk_typ(idx_topm3)=idx_mrk_dmn
mrk_typ(idx_topm4)=idx_mrk_sqr
mrk_typ(idx_topm5)=idx_mrk_str
mrk_typ(idx_topm6)=idx_mrk_tri_dwn

; Begin plot information
if (fl_out_flg) then
   fl_out=fl_out_nm
else 
   drc_out=getenv("HOME")+"/ppr_ZFM07/fgr"
   fl_out=drc_out+"/"+fgr_stb+"_"+fld_nm+fl_sfx_clr
end if
if(dvc.ne."x11") then print("Print output to "+fl_out+"."+dvc) end if

if (flg_ZGD09) then
   if (flg_vsb.and.flg_nir) then 
      data_obs=data_xpt
   else ; !both
      if (flg_vsb) then 
	 data_obs=data_xpt_vsb
      end if ; !flg_vsb
      if (flg_nir) then 
         data_obs=data_xpt_nir
      end if ; !flg_nir
   end if ; !both
   data_mdl=data_obs
   nsm_nbr=dimsizes(data_obs(:,0))

   fl_inq_ptr=addfile(fl_in(0),"r")
   var_in=fl_inq_ptr->$fld_nm$
   dat_in=var_in
   dmn_sz=dimsizes(var_in)
   mpr_nbr_in=dmn_sz(0)
   wvl_nbr_in=dmn_sz(1)
   if (mpr_nbr_in.ne.mpr_nbr) then print("ERROR: mpr_nbr_in.ne.mpr_nbr") end if
   if (wvl_nbr_in.ne.wvl_nbr) then print("ERROR: wvl_nbr_in.ne.wvl_nbr") end if
   foo=new((/1,1/),float,mss_val_flt)
   do mpr_idx=0,mpr_nbr-1
      do wvl_idx=0,wvl_nbr-1
;       print("mmr_mpr("+mpr_idx+")="+mmr_mpr(mpr_idx))
	foo=var_in({mmr_mpr(mpr_idx)*0.99:mmr_mpr(mpr_idx)*1.01},{wvl_obs(wvl_idx)*0.99:wvl_obs(wvl_idx)*1.01})
        dat_in(mpr_idx,wvl_idx)=foo(0,0)
      end do; end loop over wvl
   end do; end loop over mpr

   do xpt_idx=0,xpt_nbr-1
       data_mdl(xpt_idx,:)=dat_in(xpt2mpr(xpt_idx),xpt2wvl(xpt_idx))
       ; Apply geometrical correction if requested
       if (flg_crc) then data_mdl(xpt_idx,:)=data_mdl(xpt_idx,:)*alb_spc_snw_gmt_crc(xpt2wvl(xpt_idx)) end if ; Correct modeled values to represent 3D geometry
       ; fxm: should derive geometrical correction factor
   end do; end loop over fl
   
;   if(flg_rng) then
;      ; fxm
;   end if ; !flg_rng

end if ; !flg_ZGD09

; Scale data
if(isvar("scl_fct")) then data=scl_fct*data end if

rsr=True 
;rsr@tiMainFontHeightF=.020
;rsr@tmXBMode="Explicit"
;rsr@tmXBValues=ind_xaxis
rsr@gsnDraw=False
rsr@gsnFrame=False
rsr@tiMainString=ttl_sng
rsr@tiXAxisString=x_ttl_sng
rsr@tiYAxisFontHeightF=.020
rsr@tiYAxisOn=True
rsr@tiYAxisString=y_ttl_sng
rsr@tmXTOn=True
rsr@tmYROn=True
rsr@tmXBMajorOutwardLengthF=0.0 ; Setting to zero turns ticks inwards
rsr@tmXBMinorOutwardLengthF=0.0 ; Setting to zero turns ticks inwards
rsr@tmYLMajorOutwardLengthF=0.0 ; Setting to zero turns ticks inwards
rsr@tmYLMinorOutwardLengthF=0.0 ; Setting to zero turns ticks inwards
rsr@tmLabelAutoStride=True      ; Nice tick mark labels
if (flg_log_x) then rsr@trXLog=True end if; set X-Axis coordinate system logarithmic
if (flg_log_y) then rsr@trYLog=True end if; set Y-Axis coordinate system logarithmic
;rsr@trXMinF=1.0
rsr@trXMinF=abc_min
rsr@trXMaxF=abc_max
rsr@trYMinF=ord_min
rsr@trYMaxF=ord_max
rsr@xyDashPatterns=ln_sty(idx_top)
rsr@xyLineColors=clr_sng(idx_top)
rsr@xyLineThicknesses=ln_thk(idx_top)
rsr@xyMarkLineModes=mrk_ln_md(idx_top)
rsr@xyMarkerColors=clr_sng(idx_top)
rsr@xyMarkerSizeF=mrk_sz(idx_top)
rsr@xyMarkers=mrk_typ(idx_top)

; Open workstation and plot data
wks=gsn_open_wks(dvc,fl_out)
;plot=gsn_csm_xy(wks,data_obs,data_mdl,rsr)
plot=gsn_csm_xy(wks,(/0.0,1.0/),(/0.0,1.0/),rsr) 

rsr_xpt=True ; Polyline modifications desired
rsr_xpt@xyMarkLineMode="Markers" ; Line/marker combinations
do xpt_idx=0,xpt_nbr-1
;   print("xpt2mpr("+xpt_idx+")="+xpt2mpr(xpt_idx))
   rsr_xpt@gsMarkerColor=clr_sng(xpt2mpr(xpt_idx)) ; Color denotes impurity concentration
   if (flg_shp_wvl) then 
      rsr_xpt@gsMarkerIndex=mrk_typ(xpt2wvl(xpt_idx)) ; Marker shape denotes wavelength
   else ; !flg_shp_wvl		      
      rsr_xpt@gsMarkerIndex=mrk_typ(xpt2pss(xpt_idx)) ; Marker shape denotes date
   end if ; !flg_shp_wvl
   gph_xpt(xpt_idx)=gsn_add_polymarker(wks,plot,data_obs(xpt_idx,:),data_mdl(xpt_idx,:),rsr_xpt)
end do ; end loop over xpt_idx
rsr_ln=True ; Polyline modifications desired
rsr_ln@gsLineThicknessF=1.0 ; Line thicker
gsn_polyline(wks,plot,(/0.0,1.0/),(/0.0,1.0/),rsr_ln)
rsr_ln@gsLineDashPattern=50
gsn_polyline(wks,plot,(/0.05,1.0/),(/0.0,0.95/),rsr_ln)
gsn_polyline(wks,plot,(/0.0,0.95/),(/0.05,1.0/),rsr_ln)

if (True) then
   lgn_ln_x1=0.57
   lgn_ln_dx=0.01
   lgn_ln_x2=lgn_ln_x1+lgn_ln_dx
   lgn_ln_x1_x2=(/lgn_ln_x1,lgn_ln_x2/)
   lgn_txt_x1=lgn_ln_x1+lgn_ln_dx
   lgn_y_top=0.5
   lgn_dy=0.03
   lgn_y=lgn_y_top-ispan(0,mpr_nbr-1,1)*lgn_dy
   lgn_y=lgn_y(::-1)

   rsr_lgn=True ; Resources for legend
   rsr_lgn@txFont="helvetica"
   rsr_lgn@txFontHeightF=0.017 ; Change text size
   rsr_lgn@txJust="CenterLeft" ; Text justification
   do mpr_idx=0,mpr_nbr-1
      rsr_lgn@gsMarkerIndex=idx_mrk_crc_sld ; Circle solid
      rsr_lgn@gsMarkerColor=clr_sng(mpr_idx) ; Color denotes impurity concentration
      gsn_polymarker_ndc(wks,lgn_ln_x1,lgn_y(mpr_idx),rsr_lgn)
      gsn_text_ndc(wks,dat_sng_sht(mpr_idx),lgn_txt_x1,lgn_y(mpr_idx),rsr_lgn)
   end do ; end if
   rsr_lgn@gsMarkerColor=clr_sng(idx_top)
   if (flg_shp_wvl) then 
      ; Marker shape denotes wavelength      
      do wvl_idx=0,wvl_nbr-1
      	 rsr_lgn@gsMarkerIndex=mrk_typ(wvl_idx)
	 gsn_polymarker_ndc(wks,lgn_ln_x1,lgn_y(0)-(wvl_idx+1)*lgn_dy,rsr_lgn)
	 gsn_text_ndc(wks,wvl_sng(wvl_idx),lgn_txt_x1,lgn_y(0)-(wvl_idx+1)*lgn_dy,rsr_lgn)
      end do ; end loop over wvl
   else ; !flg_shp_wvl		      
      ; Marker shape denotes pss date
      do pss_idx=0,pss_nbr-1
      	 rsr_lgn@gsMarkerIndex=mrk_typ(pss_idx)
	 gsn_polymarker_ndc(wks,lgn_ln_x1,lgn_y(0)-(pss_idx+1)*lgn_dy,rsr_lgn)
	 gsn_text_ndc(wks,pss_sng(pss_idx),lgn_txt_x1,lgn_y(0)-(pss_idx+1)*lgn_dy,rsr_lgn)
      end do ; end loop over pss
   end if ; !flg_shp_wvl
end if ; endif add_legend

dmn_sz=dimsizes(data_obs)
lmn_nbr=dmn_sz(0)*dmn_sz(1)
obsmmdl=data_obs-data_mdl
obsmmdl_rel=(data_obs-data_mdl)/data_obs
obsmmdl_sqr=obsmmdl*obsmmdl
obsmmdl_rel_sqr=obsmmdl_rel*obsmmdl_rel
obsmmdl_sqr_ttl=dim_sum(dim_sum(obsmmdl_sqr))
obsmmdl_rel_sqr_ttl=dim_sum(dim_sum(obsmmdl_rel_sqr))
obsmmdl_rms_abs=sqrt(obsmmdl_sqr_ttl/lmn_nbr)
obsmmdl_rms_rel=sqrt(obsmmdl_rel_sqr_ttl/lmn_nbr)
print("RMS absolute bias = "+obsmmdl_rms_abs*100+"%")
print("RMS relative bias = "+obsmmdl_rms_rel*100+"%")
if(flg_rms) then
   rsr_lgn=True ; Resources for legend
   rsr_lgn@txFont="helvetica"
   rsr_lgn@txFontHeightF=0.017 ; Change text size
   rsr_lgn@txJust="CenterLeft" ; Text justification
   lgn_rms_x1=0.3
   lgn_rms_y1=0.7
   gsn_text_ndc(wks,"RMS abs. bias = "+sprintf("%5.2f",obsmmdl_rms_abs*100)+"%",lgn_rms_x1,lgn_rms_y1,rsr_lgn)
   gsn_text_ndc(wks,"RMS rel. bias = "+sprintf("%5.2f",obsmmdl_rms_rel*100)+"%",lgn_rms_x1,lgn_rms_y1-lgn_dy,rsr_lgn)
end if ; !flg_rms

draw(plot)
frame(wks) 

; Clean up for next plot
delete(clr_sng)
delete(ln_sty)
delete(ln_thk)
delete(mrk_clr)
delete(mrk_ln_md)
delete(mrk_sz)
delete(mrk_typ)
delete(rsr)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Figure fgr_sct
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
end
