; $Id$ -*-f90-*-

; Purpose: Plot line figures, generally spectral figures from mie/swnb2 output
; Based on zonal average plots in ~/anl/fgr_vy.ncl
; NB: Intended for use with defaults set in ~/.hluresfile

; Example usage:
; cd ~/anl;ncl < mie_xv.ncl
; cd ~/anl;ncl 'dvc="x11"' mie_xv.ncl
; cd ~/anl;ncl 'dvc="x11"' 'fld_nm="alb_spc_snw"' mie_xv.ncl
; cd ~/anl;ncl 'fl_in_nm="/data/zender/aca/aer_h2o_ice_snw.nc"' 'fl_out_nm="/data/zender/ps/mie_xv"' 'fld_nm="alb_spc_snw"' 'dvc="x11"' mie_xv.ncl
; cd ~/anl;ncl flg_GWM94="True" 'fl_out_nm="/data/zender/ps/mie_xv"' 'fld_nm="alb_spc_snw"' 'dvc="x11"' mie_xv.ncl
; cd ~/anl;ncl flg_crc=False flg_ZGD09=True 'fl_out_nm="/data/zender/ps/mie_xv"' 'fld_nm="alb_spc_snw"' 'ttl_sng="Measured + Modeled Snow Albedo"' 'y_ttl_sng="Lab Snow Diffuse Albedo"' 'dvc="x11"' mie_xv.ncl
; cd ~/anl;ncl 'fl_in_nm="/data/zender/aca/aer_h2o_ice_snw.nc"' 'fl_out_nm="/data/zender/ps/mie_xv"' 'fld_nm="alb_spc_snw"' 'dvc="x11"' mie_xv.ncl
; cd ~/anl;ncl dff=True 'dvc="x11"' mie_xv.ncl
; cd ~/anl;ncl flg_clr="False" mie_xv.ncl

; scp ${DATA}/ps/mie_xv.pdf sand.ess.uci.edu:/var/www/html/tmp

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

; Personal headers
if(.not.isvar("csz_ncl")) then load "$HOME/ncl/csz.ncl" end if

; Initialize some fields
prg_nm="mie_xv" ; [sng] Program name
fgr_stb="mie_xv" ; [sng] Figure name stub
dbg_lvl=0
mss_val_flt=1.0e36
mss_val_dbl=1.0d36
campaign_defined=False

if(.not.isvar("abc_max")) then abc_max=3.0e-6 end if
if(.not.isvar("abc_min")) then abc_min=0.3e-6 end if
if(.not.isvar("dvc")) then dvc="eps" end if
if(.not.isvar("fl_in_nm")) then fl_in_flg=False else fl_in_flg=True end if 
if(.not.isvar("fl_out_nm")) then fl_out_flg=False else fl_out_flg=True end if
if(.not.isvar("fld_nm")) then fld_nm="alb_spc_snw" end if
if(.not.isvar("flg_GWM94")) then flg_GWM94=False end if ; Plot GWM94 data
if(.not.isvar("flg_ZGD09")) then flg_ZGD09=False end if ; Plot ZGD09 mdl/msr
if(.not.isvar("flg_bee")) then flg_bee=False end if ; Best estimate of everything
if(.not.isvar("flg_tst")) then flg_tst=True end if ; New best estimate of everything (current manipulations)
if(.not.isvar("flg_clr")) then flg_clr=True end if
if(.not.isvar("flg_crc")) then flg_crc=False end if ; Apply geometry corrections
if(.not.isvar("flg_log_x")) then flg_log_x=False end if ; Logarithmic x-axis
if(.not.isvar("flg_nmn")) then flg_nmn=False end if ; Plot nominal snow
if(.not.isvar("flg_phb_he")) then flg_phb_he=False end if ; Plot phb_he results
if(.not.isvar("flg_phb_low")) then flg_phb_low=False end if ; Plot phb_low results
if(.not.isvar("flg_rng")) then flg_rng=False end if
if(.not.isvar("ord_max")) then ord_max=1.0 end if
if(.not.isvar("ord_min")) then ord_min=0.0 end if
if(.not.isvar("ssa_dns")) then ssa_dns=False end if ; Plot experiments predicted with measured SSA and density

if (flg_clr) then
   fl_sfx_clr="_clr"
else ; !flg_clr
   fl_sfx_clr="_bw"
end if; !flg_clr
if (flg_ZGD09) then
   flg_mdl_ZGD09=True
   flg_msr_ZGD09=True
else
   flg_mdl_ZGD09=False
   flg_msr_ZGD09=False
end if; !flg_ZGD09

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Figure mie_xv
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Plot defaults that may be over-ridden for specific variables
sym_sng=""
units_sng=""
; Set fld_sng, units_sng, and, maybe, cnt_lvl, scl_fct
if(.not.isvar("fld_prp_ncl")) then load "$HOME/ncl/fld_prp.ncl" end if
if(.not.isvar("fld_sng")) then fld_sng=fld_nm end if
if(.not.isvar("ttl_sng")) then ttl_sng=fld_sng end if
if(.not.isvar("y_ttl_sng")) then y_ttl_sng=fld_sng+" "+sym_sng end if
if(.not.isvar("x_ttl_sng")) then x_ttl_sng="Wavelength ~F33~l~F21~ [~F33~m~F21~m]" end if
y_ttl_sng=y_ttl_sng+" ["+units_sng+"]"

if(.not.isvar("ord_nbr")) then ord_nbr=7 end if ; [nbr] Ordinates read

; Legend places data bottom-to-top (small-to-large index)
; ord_nbr_plt discards ordinates with smallest indices
idx_top  =ord_nbr-1 ; [idx] top index
idx_topm1=ord_nbr-2 ; [idx] 
idx_topm2=ord_nbr-3 ; [idx] 
idx_topm3=ord_nbr-4 ; [idx] 
idx_topm4=ord_nbr-5 ; [idx] 
idx_topm5=ord_nbr-6 ; [idx] 
idx_topm6=ord_nbr-7 ; [idx] 

; Handle pre-defined dataset groups ("campaigns")
if (flg_ZGD09) then
;   if(.not.isvar(snw_rfl_ncl)) then load "$HOME/anl/snw_rfl.ncl" end if
   load "$HOME/anl/snw_rfl.ncl"
   print("quark 1 Input data file 0 = "+fl_in(0))
end if ; end if flg_ZGD09
if (flg_GWM94) then
   if(isvar("fl_in")) then delete(fl_in) end if 
   mdl_nbr=7
   if(.not.isvar("ord_nbr_plt")) then ord_nbr_plt=4 end if ; [nbr] Ordinates plotted: ord_nbr_plt <= ord_nbr
   fl_nbr=mdl_nbr
   fl_in=new(mdl_nbr,string)
   fl_in(idx_top)="/data/zender/icr/swnb_prp_dff_030_100.nc"
   ;      fl_in(idx_topm1)="/home/zender/swnb_prp_dff_100.nc"
   fl_in(idx_topm1)="/data/zender/icr/swnb_prp_dff_045_100.nc"
   fl_in(idx_topm2)="/data/zender/icr/swnb_sal_dff_045_100.nc"
   ;      fl_in(idx_topm2)="/data/zender/icr/swnb_prp_mie_045_100.nc"
   ;   fl_in(idx_topm2)="/data/zender/icr/swnb_prp_dff_noatm_045_100.nc"
   fl_in(idx_topm3)="/data/zender/icr/swnb_prp_drc_045_100.nc"
   fl_in(idx_topm4)=fl_in(idx_top)
   fl_in(idx_topm5)=fl_in(idx_top)
   fl_in(idx_topm6)=fl_in(idx_top)

   if(isvar("dat_sng")) then delete(dat_sng) end if 
   dat_sng=new(mdl_nbr,string)
   dat_sng(idx_top)="  Clean ~F10~r~B~e~N~ = ~F21~ 30, 100 ~F33~m~F21~m"
   ; dat_sng(idx_topm1)="  Clean ~F10~r~B~e~N~ = ~F21~ 100 ~F33~m~F21~m"
   dat_sng(idx_topm1)="  Clean ~F10~r~B~e~N~ = ~F21~ 45, 100 ~F33~m~F21~m"
         dat_sng(idx_topm2)="  Dirty ~F10~r~B~e~N~ = ~F21~ 45, 100 ~F33~m~F21~m"
   ;      dat_sng(idx_topm2)="  Clean Mie ~F10~r~B~e~N~ = ~F21~ 45, 100 ~F33~m~F21~m"
   ;dat_sng(idx_topm2)="  Clean NoAtm ~F10~r~B~e~N~ = ~F21~ 45, 100 ~F33~m~F21~m"
   dat_sng(idx_topm3)="  Clean, Direct ~F10~r~B~e~N~ = ~F21~ 45, 100 ~F33~m~F21~m"
   dat_sng(idx_topm4)=dat_sng(idx_top)
   dat_sng(idx_topm5)=dat_sng(idx_top)
   dat_sng(idx_topm6)=dat_sng(idx_top)

   ; Indicate to host program that input datasets have already been defined
   campaign_defined=True
end if

; Assemble file names for non-campaigns
if (.not.campaign_defined) then
   mdl_nbr=ord_nbr
   fl_nbr=mdl_nbr
   ord_nbr_plt=1
   xpt_nm="aer_h2o_ice_snw"
   drc_in=getenv("DATA")+"/aca"
   fl_in=new(fl_nbr,string)
   if (fl_in_flg) then
       fl_in(:)=fl_in_nm
   else ; !fl_in_flg
       fl_in(:)=drc_in+"/"+xpt_nm+".nc"
   end if; !fl_in_flg
   dat_sng=new(mdl_nbr,string)
   dat_sng(:)="Empty"
end if ; endif not campaign

do fl_idx=0,fl_nbr-1
   print("Input data file "+fl_idx+" = "+fl_in(fl_idx))
end do; end loop over fl

; Begin print output information
if (fl_out_flg) then
   fl_out=fl_out_nm
else 
   drc_out=getenv("DATA")+"/mie"
   fl_out=drc_out+"/"+fgr_stb+"_mie_"
   fl_out=fl_out+fld_nm+fl_sfx_clr
end if
if(dvc.ne."x11") then print("Print output will go to "+fl_out+"."+dvc) end if

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
; NB: 0*new(nbr,typ) behaves strangely due to _FillValue
ln_sty=0*new(ord_nbr,integer,"No_FillValue")
ln_thk=2+0*new(ord_nbr,integer,"No_FillValue")

; Marker properties
mrk_clr=clr_sng ; Marker color
mrk_ln_md=new((/ord_nbr/),string)
mrk_ln_md(:)="Lines" ; Line/marker combinations
mrk_sz=new(ord_nbr,float)
mrk_sz(:)=0.01
mrk_typ=new(ord_nbr,integer)
mrk_typ(idx_top)=idx_mrk_crc_sld
mrk_typ(idx_topm1)=idx_mrk_crs
mrk_typ(idx_topm2)=idx_mrk_crc
mrk_typ(idx_topm3)=idx_mrk_dmn
mrk_typ(idx_topm4)=idx_mrk_sqr
mrk_typ(idx_topm5)=idx_mrk_str
mrk_typ(idx_topm6)=idx_mrk_tri_dwn

; Determine grid size for data
fl_inq_ptr=addfile(fl_in(0),"r")
wvl=fl_inq_ptr->wvl({abc_min:abc_max}) ; [m] Wavelength
abc=wvl*1.0e6; [m]->[um]
abc_nbr=dimsizes(abc)
fld_typ=typeof(fl_inq_ptr->$fld_nm$)
if (fld_typ.eq."float") then
   data=new((/ord_nbr,abc_nbr/),float,mss_val_flt)
else
   data=new((/ord_nbr,abc_nbr/),double,mss_val_dbl)
end if

print("quark 2")
fl_lst_ptr=addfiles(fl_in,"r") ; Note plural "addfiles"
print("quark 3")
ListSetType(fl_lst_ptr,"join") ; "join" for ensembles, "cat" for timeseries
print("quark 4")
data=fl_lst_ptr[:]->$fld_nm$({abc_min:abc_max}) ; Read fld from all files
print("quark 5")
if (fld_typ.eq."double") then
   data_flt=new((/ord_nbr,abc_nbr/),float,mss_val_flt)
   data_flt=doubletofloat(data)
   delete(data)
   data=data_flt
   delete(data_flt)
end if ; endif double
delete(fl_lst_ptr)

print("quark 6")
; Scale data
if(isvar("scl_fct")) then data=scl_fct*data end if

; Average data
data_avg=new((/ord_nbr/),float,mss_val_flt)
do fl_idx=0,fl_nbr-1
    data_avg(fl_idx)=dim_avg(data(fl_idx,:)) ; Average over rightmost dimension (wvl) (and retain metadata)
end do; end loop over fl
;printVarSummary(data_avg)
;print(data)
;print(data_avg)
;print("min(data_avg)="+min(data_avg))
;print("max(data_avg)="+max(data_avg))

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
if (flg_log_x) then 
; Hard to get good labels with logarithmic axis
;   rsr@tmXBStyle="Log"
;   rsr@tmXBMode="Manual"
;   rsr@tmXBTickSpacingF=1.0
;   rsr@tmXBDataLeftF=min(abc)
;   rsr@tmXBDataRightF=max(abc)
end if ; end if flg_log_x
if (flg_log_x) then rsr@trXLog=True end if; set X-Axis coordinate system logarithmic
;rsr@trXMinF=1.0
rsr@trXMinF=abc_min*1.0e6 ; [m]->[um]
rsr@trXMaxF=abc_max*1.0e6 ; [m]->[um]
rsr@trYMinF=ord_min
rsr@trYMaxF=ord_max
rsr@xyDashPatterns=ln_sty(ord_nbr-ord_nbr_plt:ord_nbr-1)
rsr@xyLineColors=clr_sng(ord_nbr-ord_nbr_plt:ord_nbr-1)
rsr@xyLineThicknesses=ln_thk(ord_nbr-ord_nbr_plt:ord_nbr-1)
rsr@xyMarkLineModes=mrk_ln_md(ord_nbr-ord_nbr_plt:ord_nbr-1)
rsr@xyMarkerColors=clr_sng(ord_nbr-ord_nbr_plt:ord_nbr-1)
rsr@xyMarkerSizeF=mrk_sz(ord_nbr-ord_nbr_plt:ord_nbr-1)
rsr@xyMarkers=mrk_typ(ord_nbr-ord_nbr_plt:ord_nbr-1)

; Legend
rsr@lgLabelFontHeightF=0.020 ; Legend label font height
rsr@lgPerimOn=False ; Turn off box surrounding legend
rsr@pmLegendDisplayMode="Always" ; Turn on legend
rsr@pmLegendHeightF=0.20 ; Change legend height
rsr@pmLegendOrthogonalPosF=-0.5 ; Move legend down NDC units
rsr@pmLegendParallelPosF=0.65 ; Move legend right NDC units
rsr@pmLegendSide="Top" ; Change legend location
rsr@pmLegendWidthF=0.15 ; Change legend width 
rsr@xyExplicitLegendLabels=dat_sng(ord_nbr-ord_nbr_plt:ord_nbr-1) ; Legend labels

; Open workstation and plot data
wks=gsn_open_wks(dvc,fl_out)
plot=gsn_csm_xy(wks,abc,data(ord_nbr-ord_nbr_plt:ord_nbr-1,:),rsr)

if (flg_GWM94) then
   ; Data plotted in Grenfell et al. GWM94 Fgr 4 and in Flanner et al. FZR07 Fgr 2
   ; Received from Steve Warren 20050215
   wvl_GWM94=(/.31 , .32 , .33 , .34 , .35 , .36 , .37 , .38 , .4 , .452 , .464 , .492 , .51 , .5 , .523 , .543 , .555 , .569 , .584 , .591 , .6 , .606 , .62 , .627 , .634 , .64 , .656 , .672 , .685 , .714 , .732 , .745 , .774 , .806 , .862 , .903 , .95 , 1.004 , 1.042 , 1.101 , 1.134 , .76 , .784 , .841 , .878 , .904 , .945 , .971 , 1 , 1.031 , 1.045 , 1.075 , 1.104 , 1.119 , 1.133 , 1.145 , 1.179 , 1.212 , 1.239 , 1.297 , 1.334 , 1.361 , 1.4 , 1.5 , 1.6 , 1.7 , 1.8 , 1.9 , 2 , 2.15 , 2.25 , 2.38 , 2.45 , 2.5/)
   wvl_GWM94=wvl_GWM94*1.0e-6; [um]->[m]
   alb_GWM94=(/.972 , .973 , .973 , .975 , .978 , .977 , .974 , .973 , .977 , .997 , .996 , .991 , .99 , .975 , .992 , .992 , .986 , .987 , .984 , .985 , .97 , .982 , .981 , .978 , .98 , .976 , .974 , .97 , .967 , .963 , .954 , .936 , .93 , .917 , .899 , .869 , .855 , .766 , .747 , .778 , .765 , .939 , .903 , .903 , .903 , .872 , .871 , .843 , .78 , .754 , .769 , .788 , .785 , .775 , .777 , .737 , .672 , .585 , .569 , .561 , .572 , .543 , .552 , .194 , .193 , .279 , .362 , .342 , .079 , .161 , .25 , .155 , .141 , .126/)
   alb_sd_GWM94=(/0 , 0 , 0 , 0 , 0 , 0 , 0 , 0, 0 , .01 , .01 , 0 , 0 , 0 , .013 , .011 , 0 , .013 , 0 , .012 , 0 , .016 , 0 , .015 , 0 , .018 , 0 , .018 , 0 , .018 , 0 , .01 , 0 , .018 , 0 , .011 , 0 , .029 , 0 , .022 , 0 , .025 , 0 , .015 , 0 , .011 , 0 , .036 , 0 , .038 , 0 , .025 , 0 , .025 , 0 , .041 , .041 , .041 , 0 , .04 , 0 , .038 , .019 , .025 , .046 , .046 , .036 , .029 , .025 , .031 , .047 , .042 , .044 , .038/)
   abc_obs=wvl_GWM94
   abc_nbr=dimsizes(abc_obs)
   ord_nbr=1 ; Meaning there is one set of observations
   abc_2d=new((/ord_nbr,abc_nbr/),float,mss_val_flt)
   data_obs=new((/ord_nbr,abc_nbr/),float,mss_val_flt)
   data_obs_rng=new((/ord_nbr,abc_nbr/),float,mss_val_flt)
   idx_GWM94=0 ; GWM94 data are first ordinate in 2D arrays
   idx_rng=(/idx_GWM94/) ; Array of ordinates with array information
   rng_nbr=dimsizes(idx_rng) ; rng_nbr may be smaller than ord_nbr
   abc_2d(idx_GWM94,:)=wvl_GWM94*1.0e6; [m]->[um]
   data_obs(idx_GWM94,:)=alb_GWM94
   data_obs_rng(idx_GWM94,:)=alb_sd_GWM94
   
   ; Code for adding whiskers to points taken from ~/ppr_ZeM06/nco_fgr.ncl
   ; Draw range bars
   rsr_rng=True ; Polyline modifications desired
   rsr_rng@gsLineThicknessF=1.5 ; Line thicker
   ; fxm: convert to NDCs?
   rng_bar_hrz_lng=0.02
   rng_bar_ctr=new((/ord_nbr,abc_nbr/),graphic)
   rng_bar_vrt=new((/ord_nbr,abc_nbr/),graphic)
   rng_bar_hrz_top=new((/ord_nbr,abc_nbr/),graphic)
   rng_bar_hrz_btm=new((/ord_nbr,abc_nbr/),graphic)
   do rng_idx=0,rng_nbr-1 ; Loop over ordinates with error information
      ord_idx=idx_rng(rng_idx)
      do abc_idx=0,abc_nbr-1
         if (.not.ismissing(data_obs(ord_idx,abc_idx))) then
            ; Draw all points
            rng_bar_ctr(ord_idx,abc_idx)=gsn_add_polymarker(wks,plot,(/abc_2d(ord_idx,abc_idx)/),(/data_obs(ord_idx,abc_idx)/),rsr_rng)
            ; Draw ranges when greater than zero
            if (data_obs_rng(ord_idx,abc_idx).gt.0.0) then
               ; Draw vertical bar
               rng_bar_vrt(ord_idx,abc_idx)=gsn_add_polyline(wks,plot,(/abc_2d(ord_idx,abc_idx),abc_2d(ord_idx,abc_idx)/),(/data_obs(ord_idx,abc_idx)+data_obs_rng(ord_idx,abc_idx),data_obs(ord_idx,abc_idx)-data_obs_rng(ord_idx,abc_idx)/),rsr_rng)
               ; Draw horizontal line on top of vertical bar
               rng_bar_hrz_top(ord_idx,abc_idx)=gsn_add_polyline(wks,plot,(/abc_2d(ord_idx,abc_idx)-rng_bar_hrz_lng,abc_2d(ord_idx,abc_idx)+rng_bar_hrz_lng/),(/data_obs(ord_idx,abc_idx)+data_obs_rng(ord_idx,abc_idx),data_obs(ord_idx,abc_idx)+data_obs_rng(ord_idx,abc_idx)/),rsr_rng)
               ; Draw horizontal line beneath vertical bar
               rng_bar_hrz_btm(ord_idx,abc_idx)=gsn_add_polyline(wks,plot,(/abc_2d(ord_idx,abc_idx)-rng_bar_hrz_lng,abc_2d(ord_idx,abc_idx)+rng_bar_hrz_lng/),(/data_obs(ord_idx,abc_idx)-data_obs_rng(ord_idx,abc_idx),data_obs(ord_idx,abc_idx)-data_obs_rng(ord_idx,abc_idx)/),rsr_rng)
            end if; end if non-zero range
         end if; end if not missing
      end do; end loop over abc
   end do; end loop over rng

   delete(abc_2d);
   delete(data_obs);
   delete(data_obs_rng);
end if; end if flg_GWM94

if (flg_msr_ZGD09) then
; Resources for experimental points
   abc_2d=wvl_2d*1.0e6; [m]->[um]
   rsr_xpt=True ; Polyline modifications desired
   rsr_xpt@gsLineThicknessF=1.5 ; Line thicker
   rsr_xpt@xyMarkLineMode="Markers" ; Line/marker combinations

   ; Plot experiment data
   do xpt_idx=0,xpt_nbr-1
      rsr_xpt@gsMarkerColor=clr_sng(xpt2mpr(xpt_idx)) ; Color denotes impurity concentration
      rsr_xpt@gsMarkerIndex=mrk_typ(xpt2pss(xpt_idx)) ; Marker style denotes pss date
      if (flg_crc) then data_xpt(xpt_idx,:)=data_xpt(xpt_idx,:)*(1.0/alb_spc_snw_gmt_crc(xpt2wvl(xpt_idx))) end if ; Correct measured values to represent plane parallel geometry
      gph_xpt(xpt_idx)=gsn_add_polymarker(wks,plot,abc_2d(xpt_idx,:),data_xpt(xpt_idx,:),rsr_xpt)
   end do ; end loop over xpt_idx
end if ; end if flg_msr_ZGD09

print("Drawing frame...")
draw(plot)
frame(wks) 

; Clean up for next plot
delete(fl_in)
delete(clr_sng)
delete(dat_sng)
delete(data)
delete(data_avg)
delete(ln_sty)
delete(ln_thk)
delete(mrk_clr)
delete(mrk_ln_md)
delete(mrk_sz)
delete(mrk_typ)
delete(rsr)
if (flg_msr_ZGD09) then
   delete(rsr_xpt)
   delete(gph_xpt)
   delete(abc_2d)
   delete(data_xpt)
;   delete(foo)
end if ; end if flg_msr_ZGD09
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Figure mie_xv
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
end
