; $Id$

; Purpose: Spectral line plots for optical effects
; NB: Intended for use with defaults set in ~/.hluresfile

; /bin/rm ~/pnp/ppr_ZFM07/fgr/idx_rfr_clr.pdf ~/pnp/ppr_ZFM07/fgr/ssa_clr.pdf
; /bin/rm ${DATA}/ps/dst_idx_rfr.pdf
; /bin/rm ${DATA}/ps/dst_ssa.pdf

; Example usage:
; cd ~/sw/anl;ncl < spc_ln.ncl &
; cd ~/sw/anl;ncl 'dvc="x11"' spc_ln.ncl &
; cd ~/sw/anl;ncl 'dvc="eps"' spc_ln.ncl &
; cd ~/sw/anl;ncl 'dvc="x11"' 'fl_in_stb="idx_rfr_mlt.nc"' spc_ln.ncl &
; cd ~/sw/anl;ncl 'dvc="x11"' 'bln="20060621"' spc_ln.ncl &
; cd ~/sw/anl;ncl 'dvc="x11"' 'fl_in_stb="idx_rfr_bln_wlk1.nc"' 'bln="20060904"' spc_ln.ncl &
; cd ~/sw/anl;ncl flg_clr="True" spc_ln.ncl &

begin

; These files are loaded by default in NCL V6.2.0 and newer
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

; Initialize some fields
prg_nm="spc_ln"
dbg_lvl=0
mss_val_flt=1.0e36
mss_val_ntg=-999

if(.not.isvar("dvc")) then dvc="eps" end if ; [sng] Device
if(.not.isvar("flg_clr")) then flg_clr=True end if ; [flg] Plot in color
if(.not.isvar("flg_vlm")) then flg_vlm=True end if ; [flg] Plot volume-weighted refractive indices of components
if(.not.isvar("fl_in_stb")) then fl_in_stb="idx_rfr_mlt.nc" end if

if (flg_clr) then
	fl_sfx_clr="_clr"
else ; !flg_clr
	fl_sfx_clr="_bw"
end if ; !flg_clr

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
; Begin Figure 1
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
ttl_sng="Dust Imaginary Refractive Index"
x_ttl_sng="Wavelength ~F33~l~F21~ [~F33~m~F21~m]"
y_ttl_sng="Imaginary Refractive Index ~F10~n~B~i~F21~"

; Begin Arithmetic
drc_in=getenv("DATA")+"/mie"
fl_in_nm=drc_in+"/"+fl_in_stb
fl_in=addfile(fl_in_nm,"r")
wvl_mcr_mdl=doubletofloat(fl_in->wvl)*1.0e6 ; [m]->[um]
abc_mdl=wvl_mcr_mdl
abc_nbr_mdl=dimsizes(abc_mdl)
ord_nbr_mdl=6
data_mdl=new((/ord_nbr_mdl,abc_nbr_mdl/),float,mss_val_flt)

idx_mdl_ffc=0
idx_mdl_cor=1
idx_mdl_mnt=2
idx_mdl_mtx=3
idx_mdl_ncl=4
idx_mdl_prt=5

var_nm=new(ord_nbr_mdl,string)
var_nm(idx_mdl_ffc)="idx_rfr_ffc_img"
var_nm(idx_mdl_cor)="idx_rfr_cor_img"
var_nm(idx_mdl_mnt)="idx_rfr_mnt_img"
var_nm(idx_mdl_mtx)="idx_rfr_mtx_img"
var_nm(idx_mdl_ncl)="idx_rfr_ncl_img"
var_nm(idx_mdl_prt)="idx_rfr_prt_img"

wvl_mcr_aeronet=(/0.44, 0.67, 0.87, 1.02/)
wvl_mcr_TBH02=(/0.331,0.340,0.360,0.380,0.550/)
idx_rfr_TBH02_rl=(/1.58,1.58,1.57,1.58,1.56/)
idx_rfr_TBH02_img=(/0.00654,0.00616,0.00528,0.00440,0.00140/)
wvl_mcr_CTT02=(/0.360,0.360/) 
idx_rfr_CTT02_rl=(/1.0e36,1.0e36/)
idx_rfr_CTT02_img=(/0.004,0.005/) ; 

abc=wvl_mcr_aeronet
abc_nbr_aeronet=dimsizes(wvl_mcr_aeronet)
ord_nbr_aeronet=7
ord_nbr_TBH02=1
ord_nbr_CTT02=1
idx_fst_mdl=ord_nbr_aeronet+ord_nbr_TBH02+ord_nbr_CTT02-1 ; Store model data before obs data so it appears after obs data in legends (which do things backwards)
abc_nbr_TBH02=dimsizes(wvl_mcr_TBH02)
abc_nbr_CTT02=dimsizes(wvl_mcr_CTT02)
abc_nbr=max((/abc_nbr_aeronet,abc_nbr_TBH02,abc_nbr_CTT02,abc_nbr_mdl/))
ord_nbr=ord_nbr_aeronet+ord_nbr_TBH02+ord_nbr_CTT02+ord_nbr_mdl

idx_rfr_aeronet_Bhr_rl=(/1.55,1.55,1.55,1.55/)
idx_rfr_aeronet_Bhr_img=(/0.0025, 0.0014, 0.001,  0.001/)
idx_rfr_aeronet_Bnz_rl=(/1.505,1.505,1.505,1.505/)
idx_rfr_aeronet_Bnz_img=(/0.0015, 0.001,  0.001,  0.001/)
idx_rfr_aeronet_Brb_rl=(/1.43,1.43,1.43,1.43/)
idx_rfr_aeronet_Brb_img=(/0.0016, 0.0024, 0.0033, 0.0038/)
idx_rfr_aeronet_CpV_rl=(/1.48,1.48,1.48,1.48/)
idx_rfr_aeronet_CpV_img=(/0.0025, 0.0007, 0.0006, 0.0006/)
idx_rfr_aeronet_Mng_rl=(/1.555,1.555,1.555,1.555/)
idx_rfr_aeronet_Mng_img=(/0.0037, 0.0027, 0.0027, 0.0027/)
idx_rfr_aeronet_Ogd_rl=(/1.48,1.48,1.48,1.48/)
idx_rfr_aeronet_Ogd_img=(/0.0018, 0.0016, 0.0012, 0.0012/)
idx_rfr_aeronet_SdA_rl=(/1.56,1.56,1.56,1.56/)
idx_rfr_aeronet_SdA_img=(/0.0029, 0.0013, 0.001,  0.001/)
idx_rfr_aeronet_foo_rl=(/1.56,1.56,1.56,1.56/)
idx_rfr_aeronet_foo_img=(/0.0029, 0.0013, 0.001,  0.001/)

idx_rfr_DKS91_rl=(/1.53,1.53,1.53,1.53/)
idx_rfr_DKS91_img=(/0.0085,0.0045,0.0012,0.001/)
idx_rfr_PaG81_rl=(/1.56,1.56,1.56,1.56/)
idx_rfr_PaG81_img=(/0.010859,0.003774,0.0036,0.0039/)
idx_rfr_SAJ93_rl=(/1.5596,1.5599,1.5599,1.5598/)
idx_rfr_SAJ93_img=(/0.0031,0.0032,0.0036,0.0039/)

; Legend places lines in reverse order
idx_aeronet_Bhr=ord_nbr_mdl+8
idx_aeronet_Bnz=ord_nbr_mdl+7
idx_aeronet_Brb=ord_nbr_mdl+6
idx_aeronet_CpV=ord_nbr_mdl+5
idx_aeronet_Mng=ord_nbr_mdl+4
idx_aeronet_Ogd=ord_nbr_mdl+3
idx_aeronet_SdA=ord_nbr_mdl+2
idx_TBH02=ord_nbr_mdl+1
idx_CTT02=ord_nbr_mdl+0

abc_2d=new((/ord_nbr,abc_nbr/),float,mss_val_flt)
abc_2d(:,:)=mss_val_flt
abc_2d(idx_aeronet_Bhr,0:abc_nbr_aeronet-1)=wvl_mcr_aeronet
abc_2d(idx_aeronet_Bnz,0:abc_nbr_aeronet-1)=wvl_mcr_aeronet
abc_2d(idx_aeronet_Brb,0:abc_nbr_aeronet-1)=wvl_mcr_aeronet
abc_2d(idx_aeronet_CpV,0:abc_nbr_aeronet-1)=wvl_mcr_aeronet
abc_2d(idx_aeronet_Mng,0:abc_nbr_aeronet-1)=wvl_mcr_aeronet
abc_2d(idx_aeronet_Ogd,0:abc_nbr_aeronet-1)=wvl_mcr_aeronet
abc_2d(idx_aeronet_SdA,0:abc_nbr_aeronet-1)=wvl_mcr_aeronet
abc_2d(idx_CTT02,0:abc_nbr_CTT02-1)=wvl_mcr_CTT02
abc_2d(idx_TBH02,0:abc_nbr_TBH02-1)=wvl_mcr_TBH02
abc_2d(idx_mdl_ffc,0:abc_nbr_mdl-1)=wvl_mcr_mdl
abc_2d(idx_mdl_cor,0:abc_nbr_mdl-1)=wvl_mcr_mdl
abc_2d(idx_mdl_mnt,0:abc_nbr_mdl-1)=wvl_mcr_mdl
abc_2d(idx_mdl_mtx,0:abc_nbr_mdl-1)=wvl_mcr_mdl
abc_2d(idx_mdl_ncl,0:abc_nbr_mdl-1)=wvl_mcr_mdl
abc_2d(idx_mdl_prt,0:abc_nbr_mdl-1)=wvl_mcr_mdl

data_obs=new((/ord_nbr,abc_nbr/),float,mss_val_flt)
data_obs(:,:)=mss_val_flt
data_obs(idx_aeronet_Bhr,0:abc_nbr_aeronet-1)=idx_rfr_aeronet_Bhr_img
data_obs(idx_aeronet_Bnz,0:abc_nbr_aeronet-1)=idx_rfr_aeronet_Bnz_img
data_obs(idx_aeronet_Brb,0:abc_nbr_aeronet-1)=idx_rfr_aeronet_Brb_img
data_obs(idx_aeronet_CpV,0:abc_nbr_aeronet-1)=idx_rfr_aeronet_CpV_img
data_obs(idx_aeronet_Mng,0:abc_nbr_aeronet-1)=idx_rfr_aeronet_Mng_img
data_obs(idx_aeronet_Ogd,0:abc_nbr_aeronet-1)=idx_rfr_aeronet_Ogd_img
data_obs(idx_aeronet_SdA,0:abc_nbr_aeronet-1)=idx_rfr_aeronet_SdA_img
data_obs(idx_TBH02,0:abc_nbr_TBH02-1)=idx_rfr_TBH02_img
data_obs(idx_CTT02,0:abc_nbr_CTT02-1)=idx_rfr_CTT02_img
do ord_idx=0,ord_nbr_mdl-1
	data_obs(ord_idx,0:abc_nbr_mdl-1)=fl_in->$var_nm(ord_idx)$
end do; end loop over ord_idx
; following all defined by mie.cc:vlm_frc_ncl_mlt when ncl_nbr==4
if(flg_vlm) then 
	; Plot volume-weighted refractive indices of components
	vlm_frc_cor=fl_in->vlm_frc_cor
	vlm_frc_mnt=fl_in->vlm_frc_mnt
	vlm_frc_mtx=fl_in->vlm_frc_mtx
	vlm_frc_ncl=fl_in->vlm_frc_ncl
	vlm_frc_prt=fl_in->vlm_frc_prt
	data_obs(idx_mdl_cor,0:abc_nbr_mdl-1)=data_obs(idx_mdl_cor,0:abc_nbr_mdl-1)*vlm_frc_cor
	data_obs(idx_mdl_mnt,0:abc_nbr_mdl-1)=data_obs(idx_mdl_mnt,0:abc_nbr_mdl-1)*vlm_frc_mnt
	data_obs(idx_mdl_mtx,0:abc_nbr_mdl-1)=data_obs(idx_mdl_mtx,0:abc_nbr_mdl-1)*vlm_frc_mtx
	data_obs(idx_mdl_ncl,0:abc_nbr_mdl-1)=data_obs(idx_mdl_ncl,0:abc_nbr_mdl-1)*vlm_frc_ncl
	data_obs(idx_mdl_prt,0:abc_nbr_mdl-1)=data_obs(idx_mdl_prt,0:abc_nbr_mdl-1)*vlm_frc_prt
end if ; !flg_vlm

; dbg mdl
; ncks -v idx_rfr_ffc_img -d wvl,0.5e-6 ${DATA}/mie/idx_rfr_mlt.nc
dbg_sng="data_obs("+sprinti("%d",idx_mdl_ffc)+",0.5e-6)="+sprintf("%f",data_obs(idx_mdl_ffc,{0.5e-6}))
print(dbg_sng)

dat_sng=new((/ord_nbr/),string)
dat_sng(idx_aeronet_Bhr)="  AERONET: Bahrain"
dat_sng(idx_aeronet_Bnz)="  AERONET: Banizoumbou"
dat_sng(idx_aeronet_Brb)="  AERONET: Barbados"
dat_sng(idx_aeronet_CpV)="  AERONET: Cape Verde"
dat_sng(idx_aeronet_Mng)="  AERONET: Mongolia"
dat_sng(idx_aeronet_Ogd)="  AERONET: Ougadougou"
dat_sng(idx_aeronet_SdA)="  AERONET: Solar Village"
dat_sng(idx_TBH02)="  TOMS/AERONET (Sinyuk et al.)" ; STD03
dat_sng(idx_CTT02)="  CARMA/TOMS (Colarco et al.)" ; 
dat_sng(idx_mdl_ffc)="  Effective Medium"
dat_sng(idx_mdl_cor)="  "+fl_in@composition_core
dat_sng(idx_mdl_mnt)="  "+fl_in@composition_mantle
dat_sng(idx_mdl_mtx)="  "+fl_in@composition_matrix
dat_sng(idx_mdl_ncl)="  "+fl_in@composition_inclusion
dat_sng(idx_mdl_prt)="  "+fl_in@composition_particle

clr_sng=new((/ord_nbr/),string)
clr_sng(:)="black"
if (flg_clr) then
	clr_sng(idx_aeronet_Bnz)="black"
	clr_sng(idx_aeronet_Brb)="black"
	clr_sng(idx_aeronet_Bhr)="black"
	clr_sng(idx_aeronet_CpV)="black"
	clr_sng(idx_aeronet_Mng)="black"
	clr_sng(idx_aeronet_Ogd)="black"
	clr_sng(idx_aeronet_SdA)="black"
	clr_sng(idx_TBH02)="black"
	clr_sng(idx_CTT02)="black"
	clr_sng(idx_mdl_ffc)="blue"
	clr_sng(idx_mdl_cor)="green"
	clr_sng(idx_mdl_mnt)="green"
	clr_sng(idx_mdl_mtx)="green"
	clr_sng(idx_mdl_ncl)="green"
	clr_sng(idx_mdl_prt)="green"
end if ; !flg_clr

ln_sty=0*new(ord_nbr,integer,"No_FillValue")
ln_sty(idx_aeronet_Bhr)=0
ln_sty(idx_aeronet_Bnz)=0
ln_sty(idx_aeronet_Brb)=0
ln_sty(idx_aeronet_CpV)=0
ln_sty(idx_aeronet_Mng)=0
ln_sty(idx_aeronet_Ogd)=0
ln_sty(idx_aeronet_SdA)=0
ln_sty(idx_TBH02)=50
ln_sty(idx_CTT02)=50
ln_sty(idx_mdl_ffc)=0
ln_sty(idx_mdl_cor)=0
ln_sty(idx_mdl_mnt)=0
ln_sty(idx_mdl_mtx)=0
ln_sty(idx_mdl_ncl)=0
ln_sty(idx_mdl_prt)=0

ln_thk=new(ord_nbr,integer,mss_val_ntg)
ln_thk(:)=2

; Marker properties
mrk_clr=clr_sng ; Marker color
mrk_ln_md=new((/ord_nbr/),string)
mrk_ln_md(:)="MarkLines" ; Line/marker combinations
mrk_sz=new(ord_nbr,float,mss_val_flt)
mrk_sz(:)=0.01
mrk_typ=new(ord_nbr,integer,mss_val_ntg)
mrk_typ(idx_aeronet_Bnz)=idx_mrk_crc_sld
mrk_typ(idx_aeronet_Bhr)=idx_mrk_crs
mrk_typ(idx_aeronet_Brb)=idx_mrk_crc_crs
mrk_typ(idx_aeronet_CpV)=idx_mrk_crc_dot
mrk_typ(idx_aeronet_Mng)=idx_mrk_dmn
mrk_typ(idx_aeronet_Ogd)=idx_mrk_pls
mrk_typ(idx_aeronet_SdA)=idx_mrk_sqr
mrk_typ(idx_TBH02)=idx_mrk_tri_up
mrk_typ(idx_CTT02)=idx_mrk_tri_dwn
mrk_typ(idx_mdl_ffc)=idx_mrk_crc_sld
mrk_typ(idx_mdl_cor)=idx_mrk_crs
mrk_typ(idx_mdl_mnt)=idx_mrk_dmn
mrk_typ(idx_mdl_mtx)=idx_mrk_crc_crs
mrk_typ(idx_mdl_ncl)=idx_mrk_crc
mrk_typ(idx_mdl_prt)=idx_mrk_pls

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
rsr@tmXTOn=False
rsr@tmYROn=False
rsr@trXMinF=0.3
rsr@trXMaxF=1.1
rsr@trYMinF=0.0
rsr@trYMaxF=0.01
rsr@xyDashPatterns=ln_sty
rsr@xyLineColors=clr_sng
rsr@xyLineThicknesses=ln_thk
rsr@xyMarkLineModes=mrk_ln_md
rsr@xyMarkerColors=clr_sng
rsr@xyMarkerSizeF=mrk_sz
rsr@xyMarkers=mrk_typ

; Legend
rsr@lgLabelFontHeightF=0.015 ; Legend label font height
rsr@lgPerimOn=False ; Turn off box surrounding legend
rsr@pmLegendDisplayMode="Always" ; Turn on legend
rsr@pmLegendHeightF=0.33 ; Change legend height
rsr@pmLegendOrthogonalPosF=-0.65 ; Move legend down NDC units
rsr@pmLegendParallelPosF=0.5 ; Move legend right NDC units
rsr@pmLegendSide="Top" ; Change legend location
rsr@pmLegendWidthF=0.15 ; Change legend width 
rsr@xyExplicitLegendLabels=dat_sng ; Legend labels
;res@gsnDraw           = False
;res@gsnFrame          = False

; Open workstation and plot data
drc_out="${DATA}/ps"
;fl_out=drc_out+"/dst_idx_rfr"+fl_sfx_clr
fl_out=drc_out+"/dst_idx_rfr"
if (dvc.ne."x11") then print("Output will go to "+fl_out+"."+dvc) end if
wks=gsn_open_wks(dvc,fl_out)
plot=gsn_csm_xy(wks,abc_2d,data_obs,rsr)
;plot=gsn_csm_xy(wks,abc_mdl,data_mdl,rsr)

draw(plot)

frame(wks)

; Clean up for next plot
delete(abc)
delete(abc_2d)
delete(abc_mdl)
delete(clr_sng)
delete(dat_sng)
delete(data_mdl)
delete(data_obs)
delete(fl_in)
delete(fl_in_nm)
delete(ln_sty)
delete(ln_thk)
delete(mrk_clr)
delete(mrk_ln_md)
delete(mrk_sz)
delete(mrk_typ)
delete(rsr)
delete(var_nm)
delete(wvl_mcr_mdl)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Figure 1
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Figure 2
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
ttl_sng="Dust Single Scattering Albedo"
x_ttl_sng="Wavelength ~F33~l~F21~ [~F33~m~F21~m]"
y_ttl_sng="Single Scattering Albedo ~F33~w~F21~"

; Command line arguments for this plot
if(.not.isvar("bln")) then bln="20060904" end if

; Begin Arithmetic
fl_in_nm=(/"$DATA/dst/mie/aer_dst_bln_"+bln+"_04_snicar.nc",\
	"$DATA/dst/mie/aer_dst_bln_"+bln+"_03_snicar.nc",\
	"$DATA/dst/mie/aer_dst_bln_"+bln+"_02_snicar.nc",\
	"$DATA/dst/mie/aer_dst_bln_"+bln+"_01_snicar.nc"/)
fl_in=addfiles(fl_in_nm,"r")
wvl_mcr_mdl=doubletofloat(fl_in[0]->wvl)*1.0e6 ; [m]->[um] NB: square brackets
abc_mdl=wvl_mcr_mdl
abc_nbr_mdl=dimsizes(abc_mdl)
ord_nbr_mdl=dimsizes(fl_in_nm)
data_mdl=new((/ord_nbr_mdl,abc_nbr_mdl/),float,mss_val_flt)

; Quantities used for derived fields
ext_cff_mss_mdl=new((/ord_nbr_mdl,abc_nbr_mdl/),float,mss_val_flt)
ssa_mdl=new((/ord_nbr_mdl,abc_nbr_mdl/),float,mss_val_flt)
mpc_dst=new((/ord_nbr_mdl/),float,mss_val_flt)
;mpc_dst=(/7.54052895365e-06,2.03528735631e-05,1.6365892713e-05,2.34241234801e-05/) ; [kg m-2] Global mean column paths of bins 1--4 from cspdfdb07
;mpc_dst=(/4.22817e-06,1.13764e-05,8.86696e-06,1.28381e-05/) ; [kg m-2] Global mean column paths of bins 1--4 from esmfipcc01 year 0
;mpc_dst=(/4.0527e-05,0.000112495,0.000110711,0.000216529/) ; [kg m-2] Bahrain column paths of bins 1--4 from esmfipcc01 year 0
;mpc_dst=(/4.07515e-05,0.000109976,6.91932e-05,6.4165e-05/) ; [kg m-2] Cape Verde column paths of bins 1--4 from esmfipcc01 years 0-4
;mpc_dst=(/3.54606e-05,9.76973e-05,8.42204e-05,0.000122992/) ; [kg m-2] Solar Village column paths of bins 1--4 from esmfipcc01 years 0-4 NB: These Solar Village data were used to generate figures disseminated from 2007-2013 as dst_idx_rfr_ssa.pdf
mpc_dst=(/0.025952,0.192418,0.421066,0.360564/) ; [kg m-2] Samuel Albani supplied these Solar Village simulations 20130212
;mpc_dst=(/0.1,0.2,0.3,0.4/) ; [kg m-2] Bahrain ~observed coarse/fine fractions DHE02
;mpc_dst=(/0.05,0.2,0.3,0.5/) ; [kg m-2] Test paths
;mpc_dst=(/1.0,0.0,0.0,0.0/) ; [kg m-2] Test paths
; Arrays are stored backwards from size order (uggh, fxm)
mpc_dst=mpc_dst(::-1) ; http://www.ncl.ucar.edu/Document/Language/array_reverse.shtml
mpc_dst_frc=mpc_dst/sum(mpc_dst) ; Normalized fractional weights of each size bin

idx_mdl_bin4=0
idx_mdl_bin3=1
idx_mdl_bin2=2
idx_mdl_bin1=3

var_nm=new(ord_nbr_mdl,string)
var_nm(idx_mdl_bin4)="ss_alb"
var_nm(idx_mdl_bin3)="ss_alb"
var_nm(idx_mdl_bin2)="ss_alb"
var_nm(idx_mdl_bin1)="ss_alb"

wvl_mcr_landsat=(/0.47,0.64,0.83/) ; Read carefully KTD01 p. 1481 Landsat channels used are blue, red, and NIR: http://public.lanl.gov/kelly/CANDID/tm/index.shtml and 
wvl_mcr_aeronet=(/0.44, 0.67, 0.87, 1.02/)
; drv is a placeholder for future (non-AERONET) ssa model/observations
wvl_mcr_drv=abc_mdl ;

abc_nbr_aeronet=dimsizes(wvl_mcr_aeronet)
abc_nbr_landsat=dimsizes(wvl_mcr_landsat)
ord_nbr_aeronet=5
ord_nbr_landsat=1
ord_nbr_drv=1
idx_fst_mdl=ord_nbr_aeronet+ord_nbr_landsat+ord_nbr_drv-1 ; Store model data before obs data so it appears after obs data in legends (which do things backwards)
abc_nbr_drv=dimsizes(wvl_mcr_drv)
abc_nbr=max((/abc_nbr_aeronet,abc_nbr_landsat,abc_nbr_drv,abc_nbr_mdl/))
ord_nbr=ord_nbr_aeronet+ord_nbr_landsat+ord_nbr_drv+ord_nbr_mdl

; #s following +- are measurement standard deviations, not uncertainties
ssa_aeronet_Bhr=(/0.92,0.95,0.96,0.97/) ; +- 0.03
;ssa_aeronet_Bnz=(//)
ssa_aeronet_CpV=(/0.93,0.98,0.98,0.99/) ; +- 0.01
ssa_aeronet_CpV_KTD01_05=(/0.906, 0.930, 0.968, 0.975/) ; +- 0.02 eyeballed from KTD01 fgr. 4 for tau = 0.5
ssa_aeronet_CpV_KTD01_10=(/0.900, 0.947, 0.974, 0.982/) ; +- 0.02 eyeballed from KTD01 fgr. 4 for tau = 1.0
ssa_aeronet_SdA=(/0.92,0.96,0.97,0.97/) ; +- 0.02
ssa_landsat_Sng_KTD01=(/0.87,0.99,0.998/)

asm_aeronet_Bhr=(/0.68,0.66,0.66,0.66/) ; +- 0.04
;asm_aeronet_Bnz=(//)
asm_aeronet_CpV=(/0.73,0.71,0.71,0.71/) ; +- 0.04
;asm_aeronet_CpV_KTD01_05=(//)
;asm_aeronet_CpV_KTD01_10=(//)
asm_aeronet_SdA=(/0.69,0.66,0.65,0.65/) ; +- 0.04
;asm_landsat_Sng_KTD01=(//) ; +- 0.03,0.02

;ssa_DKS91=(//)
;ssa_PaG81=(//)
;ssa_SAJ93=(//)

; Legend places lines in reverse order
idx_aeronet_Bhr=ord_nbr_mdl+6
;idx_aeronet_Bnz=ord_nbr_mdl+6
idx_aeronet_SdA=ord_nbr_mdl+5
idx_aeronet_CpV=ord_nbr_mdl+4
idx_aeronet_CpV_KTD01_05=ord_nbr_mdl+3
idx_aeronet_CpV_KTD01_10=ord_nbr_mdl+2
idx_landsat_Sng_KTD01=ord_nbr_mdl+1
idx_drv=ord_nbr_mdl+0

abc_2d=new((/ord_nbr,abc_nbr/),float,mss_val_flt)
abc_2d(:,:)=mss_val_flt
abc_2d(idx_aeronet_Bhr,0:abc_nbr_aeronet-1)=wvl_mcr_aeronet
;abc_2d(idx_aeronet_Bnz,0:abc_nbr_aeronet-1)=wvl_mcr_aeronet
abc_2d(idx_aeronet_CpV,0:abc_nbr_aeronet-1)=wvl_mcr_aeronet
abc_2d(idx_aeronet_CpV_KTD01_05,0:abc_nbr_aeronet-1)=wvl_mcr_aeronet
abc_2d(idx_aeronet_CpV_KTD01_10,0:abc_nbr_aeronet-1)=wvl_mcr_aeronet
abc_2d(idx_aeronet_SdA,0:abc_nbr_aeronet-1)=wvl_mcr_aeronet
abc_2d(idx_landsat_Sng_KTD01,0:abc_nbr_landsat-1)=wvl_mcr_landsat
abc_2d(idx_drv,0:abc_nbr_drv-1)=wvl_mcr_drv
abc_2d(idx_mdl_bin4,0:abc_nbr_mdl-1)=wvl_mcr_mdl
abc_2d(idx_mdl_bin3,0:abc_nbr_mdl-1)=wvl_mcr_mdl
abc_2d(idx_mdl_bin2,0:abc_nbr_mdl-1)=wvl_mcr_mdl
abc_2d(idx_mdl_bin1,0:abc_nbr_mdl-1)=wvl_mcr_mdl

data_obs=new((/ord_nbr,abc_nbr/),float,mss_val_flt)
data_obs(:,:)=mss_val_flt
data_obs(idx_aeronet_Bhr,0:abc_nbr_aeronet-1)=ssa_aeronet_Bhr
;data_obs(idx_aeronet_Bnz,0:abc_nbr_aeronet-1)=ssa_aeronet_Bnz
data_obs(idx_aeronet_CpV,0:abc_nbr_aeronet-1)=ssa_aeronet_CpV
data_obs(idx_aeronet_CpV_KTD01_05,0:abc_nbr_aeronet-1)=ssa_aeronet_CpV_KTD01_05
data_obs(idx_aeronet_CpV_KTD01_10,0:abc_nbr_aeronet-1)=ssa_aeronet_CpV_KTD01_10
data_obs(idx_aeronet_SdA,0:abc_nbr_aeronet-1)=ssa_aeronet_SdA
data_obs(idx_landsat_Sng_KTD01,0:abc_nbr_landsat-1)=ssa_landsat_Sng_KTD01
do ord_idx=0,ord_nbr_mdl-1
	data_obs(ord_idx,0:abc_nbr_mdl-1)=doubletofloat(fl_in[ord_idx]->$var_nm(ord_idx)$) ; NB: square brackets
	ext_cff_mss_mdl(ord_idx,0:abc_nbr_mdl-1)=doubletofloat(fl_in[ord_idx]->ext_cff_mss) ; [m2 kg-1] Mass extinction coefficient
	ssa_mdl(ord_idx,0:abc_nbr_mdl-1)=doubletofloat(fl_in[ord_idx]->ss_alb) ; [fraction] Single scattering albedo
end do; end loop over ord_idx
data_obs!1="wvl" ; Name dimension
data_obs&wvl=doubletofloat(fl_in[0]->wvl) ; Assign coordinate

; dbg mdl
; ncks -C -v ss_alb -d wvl,0.5e-6 ${DATA}/dst/mie/aer_dst_bln_20060621_04_snicar.nc
dbg_sng="data_obs("+sprinti("%d",idx_mdl_bin4)+",0.5e-6)="+sprintf("%f",data_obs(idx_mdl_bin4,{0.5e-6}))
print(dbg_sng)
; Derive variables
; Optical depth matrix is dst_nbr x wvl_nbr
mpc_dst_frc_2d=conform(ext_cff_mss_mdl,mpc_dst_frc,0) ; Broadcast mpc_dst_frc along dimension 0 to fit shape of ext_cff_mss_mdl
odxc_mdl=ext_cff_mss_mdl*mpc_dst_frc_2d ; [frc] Spectral optical depth in each size
odxc_mdl_ttl=dim_sum(transpose(odxc_mdl)) ; [frc] Total spectral optical depth
; Weighted single scatter albedo is array of size wvl_nbr
ssa_drv=dim_sum(transpose(ssa_mdl*odxc_mdl))/odxc_mdl_ttl
;printVarSummary(ssa_drv)
copy_VarCoords(wvl_mcr_drv,ssa_drv) ; Copy named dimensions and coordinate variables from a to b
data_obs(idx_drv,0:abc_nbr_drv-1)=(/ssa_drv/) ; Use value assign to avoid writing over data_obs metadata

dat_sng=new((/ord_nbr/),string)
dat_sng(idx_aeronet_Bhr)="  AERONET: Bahrain"
;dat_sng(idx_aeronet_Bnz)="  AERONET: Banizoumbou"
dat_sng(idx_aeronet_CpV)="  AERONET: Cape Verde"
dat_sng(idx_aeronet_CpV_KTD01_05)="  AERONET: Cape Verde (~F33~t~F21~=0.5)"
dat_sng(idx_aeronet_CpV_KTD01_10)="  AERONET: Cape Verde (~F33~t~F21~=1.0)"
dat_sng(idx_aeronet_SdA)="  AERONET: Solar Village"
dat_sng(idx_landsat_Sng_KTD01)="  Landsat: Senegal (~F33~t~F21~=2.4)"
dat_sng(idx_drv)="  DEAD 0.1 < ~F10~D~F21~ < 10.0 ~F33~m~F21~m"
dat_sng(idx_mdl_bin4)="  DEAD 5.0 < ~F10~D~F21~ < 10.0 ~F33~m~F21~m"
dat_sng(idx_mdl_bin3)="  DEAD 2.5 < ~F10~D~F21~ < 5.0 ~F33~m~F21~m"
dat_sng(idx_mdl_bin2)="  DEAD 1.0 < ~F10~D~F21~ < 2.5 ~F33~m~F21~m"
dat_sng(idx_mdl_bin1)="  DEAD 0.1 < ~F10~D~F21~ < 1.0 ~F33~m~F21~m"

clr_sng=new((/ord_nbr/),string)
clr_sng(:)="black"
if (flg_clr) then
	clr_sng(idx_aeronet_Bhr)="black"
;	clr_sng(idx_aeronet_Bnz)="black"
	clr_sng(idx_aeronet_CpV)="black"
	clr_sng(idx_aeronet_CpV_KTD01_05)="black"
	clr_sng(idx_aeronet_CpV_KTD01_10)="black"
	clr_sng(idx_aeronet_SdA)="black"
	clr_sng(idx_landsat_Sng_KTD01)="black"
	clr_sng(idx_drv)="blue"
	clr_sng(idx_mdl_bin4)="green"
	clr_sng(idx_mdl_bin3)="green"
	clr_sng(idx_mdl_bin2)="green"
	clr_sng(idx_mdl_bin1)="green"
end if ; !flg_clr

ln_sty=0*new(ord_nbr,integer,"No_FillValue")
ln_sty(idx_aeronet_Bhr)=0
;ln_sty(idx_aeronet_Bnz)=0
ln_sty(idx_aeronet_CpV)=0
ln_sty(idx_aeronet_CpV_KTD01_05)=0
ln_sty(idx_aeronet_CpV_KTD01_10)=0
ln_sty(idx_aeronet_SdA)=0
ln_sty(idx_landsat_Sng_KTD01)=0
ln_sty(idx_drv)=0
ln_sty(idx_mdl_bin4)=0
ln_sty(idx_mdl_bin3)=0
ln_sty(idx_mdl_bin2)=0
ln_sty(idx_mdl_bin1)=0

ln_thk=new(ord_nbr,integer,mss_val_ntg)
ln_thk(:)=2

; Marker properties
mrk_clr=clr_sng ; Marker color
mrk_ln_md=new((/ord_nbr/),string)
mrk_ln_md(:)="MarkLines" ; Line/marker combinations
mrk_sz=new(ord_nbr,float,mss_val_flt)
mrk_sz(:)=0.01
mrk_typ=new(ord_nbr,integer,mss_val_ntg)
mrk_typ(idx_aeronet_Bhr)=idx_mrk_crs
;mrk_typ(idx_aeronet_Bnz)=idx_mrk_crc_sld
mrk_typ(idx_aeronet_CpV)=idx_mrk_crc_dot
mrk_typ(idx_aeronet_CpV_KTD01_05)=idx_mrk_dmn
mrk_typ(idx_aeronet_CpV_KTD01_10)=idx_mrk_pls
mrk_typ(idx_aeronet_SdA)=idx_mrk_sqr
mrk_typ(idx_landsat_Sng_KTD01)=idx_mrk_crc_crs
mrk_typ(idx_drv)=idx_mrk_tri_dwn
mrk_typ(idx_mdl_bin4)=idx_mrk_crc_crs
mrk_typ(idx_mdl_bin3)=idx_mrk_crs
mrk_typ(idx_mdl_bin2)=idx_mrk_dmn
mrk_typ(idx_mdl_bin1)=idx_mrk_crc_sld

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
rsr@tmXTOn=False
rsr@tmYROn=False
rsr@trXMinF=0.3
rsr@trXMaxF=1.1
rsr@trYMinF=0.8
rsr@trYMaxF=1.0
rsr@xyDashPatterns=ln_sty
rsr@xyLineColors=clr_sng
rsr@xyLineThicknesses=ln_thk
rsr@xyMarkLineModes=mrk_ln_md
rsr@xyMarkerColors=clr_sng
rsr@xyMarkerSizeF=mrk_sz
rsr@xyMarkers=mrk_typ

; Legend
rsr@lgLabelFontHeightF=0.015 ; Legend label font height
rsr@lgPerimOn=False ; Turn off box surrounding legend
rsr@pmLegendDisplayMode="Always" ; Turn on legend
rsr@pmLegendHeightF=0.30 ; Change legend height
rsr@pmLegendOrthogonalPosF=-1.1 ; Move legend down NDC units
rsr@pmLegendParallelPosF=0.66 ; Move legend right NDC units
rsr@pmLegendSide="Top" ; Change legend location
rsr@pmLegendWidthF=0.15 ; Change legend width 
rsr@xyExplicitLegendLabels=dat_sng ; Legend labels
;res@gsnDraw           = False
;res@gsnFrame          = False

; Open workstation and plot data
drc_out="${DATA}/ps"
;fl_out=drc_out+"/dst_ssa"+fl_sfx_clr
fl_out=drc_out+"/dst_ssa"
if (dvc.ne."x11") then print("Output will go to "+fl_out+"."+dvc) end if
wks=gsn_open_wks(dvc,fl_out)
plot=gsn_csm_xy(wks,abc_2d,data_obs,rsr)
;plot=gsn_csm_xy(wks,abc_mdl,data_mdl,rsr)

draw(plot)

frame(wks)

; Clean up for next plot
delete(abc_2d)
delete(abc_mdl)
delete(clr_sng)
delete(dat_sng)
delete(data_mdl)
delete(data_obs)
delete(fl_in)
delete(fl_in_nm)
delete(ln_sty)
delete(ln_thk)
delete(mrk_clr)
delete(mrk_ln_md)
delete(mrk_sz)
delete(mrk_typ)
delete(rsr)
delete(var_nm)
delete(wvl_mcr_mdl)

; Delete extra diagnostic inputs
delete(ext_cff_mss_mdl)
delete(ssa_mdl)
delete(mpc_dst)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Figure 2
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

end
