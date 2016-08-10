; $Id$

; Pupose: Compute TDDR band average absorption

; Usage:
; ncl < ~/sw/ncl/tddr_abs.ncl

begin ; begin tddr_abs()

; Initialize some fields
mss_val=1.0e36
bnd_ctr=0.5e-6
bnd_dlt=10.0e-9
bnd_lft=bnd_ctr-0.5*bnd_dlt
bnd_rgt=bnd_ctr+0.5*bnd_dlt
day=83 ; Day > 31 is a code to use SS study files instead
;day=73 ; Day > 31 is a code to use SS study files instead
;day=75 ; Day > 31 is a code to use SS study files instead
;day=80 ; Day > 31 is a code to use SS study files instead
;day=85 ; Day > 31 is a code to use SS study files instead
;day=90 ; Day > 31 is a code to use SS study files instead
;day=95 ; Day > 31 is a code to use SS study files instead
day=99 ; Day > 31 is a code to use SS study files instead
;day=11
;day=30
;aer_sng="cln"
aer_sng="aer"
prs_top=16200.0 ; Pa
if(day.eq.30) then prs_btm=93200.0 end if ; Pa
if(day.eq.11) then prs_btm=83200.0 end if ; Pa z=1456.37 m
if(day.gt.31) then prs_btm=83200.0 end if ; Pa z=1456.37 m
;if(day.eq.11) then prs_btm=81000.0 end if ; Pa z=1717.0 m
;if(day.eq.11) then prs_btm=93500.0 end if ; Pa z=811.9 m
; Diagnose column values
;prs_top=50.0 ; Pa
;prs_btm=100000.0 ; Pa
prs_top_sng=floattochar(prs_top)
prs_btm_sng=floattochar(prs_btm)

; Get data
;fl_nm="${DATA}/arese/mdl/951011_1200_arese_mdl_clr_aer.nc"
;drc_sng="${DATA}/arese/mdl/"
drc_sng="${DATA}/tmp/"
if(day.gt.31) then 
	day_sng=inttochar(day)
	if(day.eq.83) then ss_alb_sng="083" end if
	if(day.eq.73) then ss_alb_sng="000" end if
	if(day.eq.75) then ss_alb_sng="075" end if
	if(day.eq.80) then ss_alb_sng="080" end if
	if(day.eq.85) then ss_alb_sng="085" end if
	if(day.eq.90) then ss_alb_sng="090" end if
	if(day.eq.95) then ss_alb_sng="095" end if
	if(day.eq.99) then ss_alb_sng="099" end if
	if(day.eq.100) then ss_alb_sng="100" end if
	fl_nm=(/\
        drc_sng+"arese_swnb_w_"+ss_alb_sng+"_AOD_000.nc",\
        drc_sng+"arese_swnb_w_"+ss_alb_sng+"_AOD_006.nc",\
        drc_sng+"arese_swnb_w_"+ss_alb_sng+"_AOD_012.nc"\
	/) 
end if ; endif day.gt.31
if(day.eq.11) then 
	day_sng="951011"
	fl_nm=(/\
        drc_sng+day_sng+"_0900_arese_mdl_clr_"+aer_sng+".nc",\
        drc_sng+day_sng+"_1000_arese_mdl_clr_"+aer_sng+".nc",\
        drc_sng+day_sng+"_1100_arese_mdl_clr_"+aer_sng+".nc",\
        drc_sng+day_sng+"_1200_arese_mdl_clr_"+aer_sng+".nc",\
        drc_sng+day_sng+"_1300_arese_mdl_clr_"+aer_sng+".nc"\
	/) 
end if ; endif day.eq.11
if(day.eq.30) then 
	day_sng="951030"
	fl_nm=(/\
        drc_sng+day_sng+"_0900_arese_mdl_cld_"+aer_sng+".nc",\
        drc_sng+day_sng+"_1000_arese_mdl_cld_"+aer_sng+".nc",\
        drc_sng+day_sng+"_1100_arese_mdl_cld_"+aer_sng+".nc",\
        drc_sng+day_sng+"_1200_arese_mdl_cld_"+aer_sng+".nc",\
        drc_sng+day_sng+"_1300_arese_mdl_cld_"+aer_sng+".nc"\
	/) 
end if ; endif day.eq.30
fl_nbr=dimsizes(fl_nm)

print("Quantities computed for YYMMDD = "+day_sng+", conditions = "+aer_sng)
fl_in=addfile(fl_nm(0),"r")
alt_top=fl_in->alt_int({prs_top})
alt_btm=fl_in->alt_int({prs_btm})
print("Egrett level = "+fl_in->levp({prs_top})/100.+" mb = "+alt_top/1000.+" km")
print("Otter level = "+fl_in->levp({prs_btm})/100.+" mb = "+alt_btm/1000.+" km")
print("TST	flx_dwn_egr_ttl	flx_dwn_ott_ttl flx_abs_ttl	abs_clm		trn_clm")
print("HH:MM	W/m2		W/m2		W/m2		fraction	fraction")
abs_clm=new((/fl_nbr/),float)
trn_clm=new((/fl_nbr/),float)
odxc_aer=new((/fl_nbr/),float)
do fl_idx=0,fl_nbr-1
	fl_in=addfile(fl_nm(fl_idx),"r")
	odxc_aer(fl_idx)=fl_in->odxc_obs_aer
	lcl_time_hr=fl_in->lcl_time_hr
	bnd=fl_in->bnd({bnd_lft:bnd_rgt})
	bnd_nbr=dimsizes(bnd)
	lev=fl_in->lev({prs_top:prs_btm})
	lev_nbr=dimsizes(lev)
	flx_spc_abs=fl_in->flx_spc_abs({prs_top:prs_btm},{bnd_lft:bnd_rgt})
	wvl_sz=fl_in->wvl_sz({bnd_lft:bnd_rgt})
	wvl_sz_lev=onedtond(wvl_sz,(/lev_nbr,bnd_nbr/))
	flx_abs=wvl_sz_lev*flx_spc_abs
; Normalize model bands to TDDR bandwidth (assume SRF = 1)
	flx_abs_ttl=(bnd_dlt/sum(wvl_sz))*sum(flx_abs) 

; Absolute down flux in given channel at Egrett
	flx_spc_dwn_egr=fl_in->flx_spc_dwn({prs_top},{bnd_lft:bnd_rgt})
	flx_dwn_egr=wvl_sz*flx_spc_dwn_egr
	flx_dwn_egr_ttl=(bnd_dlt/sum(wvl_sz))*sum(flx_dwn_egr)

; Absolute down flux in given channel at Otter
	flx_spc_dwn_ott=fl_in->flx_spc_dwn({prs_btm},{bnd_lft:bnd_rgt})
	flx_dwn_ott=wvl_sz*flx_spc_dwn_ott
	flx_dwn_ott_ttl=(bnd_dlt/sum(wvl_sz))*sum(flx_dwn_ott)

; Derived fields absorptance and transmittance
	abs_clm(fl_idx)=flx_abs_ttl/flx_dwn_egr_ttl
	trn_clm(fl_idx)=flx_dwn_ott_ttl/flx_dwn_egr_ttl

	print(lcl_time_hr+":00	"+flx_dwn_egr_ttl+"		"+flx_dwn_ott_ttl+"		"+flx_abs_ttl+"	"+abs_clm(fl_idx)+"	"+trn_clm(fl_idx))
end do; end loop over fl_idx

print("")
print("TST	flx_dwn_dff_ott_ttl	flx_dwn_drc_ott_ttl	dff_drc_rat")
print("HH:MM	W/m2			W/m2			ratio")
do fl_idx=0,fl_nbr-1
	fl_in=addfile(fl_nm(fl_idx),"r")
	lcl_time_hr=fl_in->lcl_time_hr
	wvl_sz=fl_in->wvl_sz({bnd_lft:bnd_rgt})

; Absolute down direct and diffuse fluxes in given channel at Otter
	flx_spc_dwn_dff_ott=fl_in->flx_spc_dwn_dff({prs_btm},{bnd_lft:bnd_rgt})
	flx_dwn_dff_ott=wvl_sz*flx_spc_dwn_dff_ott
	flx_dwn_dff_ott_ttl=(bnd_dlt/sum(wvl_sz))*sum(flx_dwn_dff_ott)

	flx_spc_dwn_drc_ott=fl_in->flx_spc_dwn_drc({prs_btm},{bnd_lft:bnd_rgt})
	flx_dwn_drc_ott=wvl_sz*flx_spc_dwn_drc_ott
	flx_dwn_drc_ott_ttl=(bnd_dlt/sum(wvl_sz))*sum(flx_dwn_drc_ott)

; Derived fields direct:diffuse ratio
	dff_drc_rat=flx_dwn_dff_ott_ttl/flx_dwn_drc_ott_ttl

print(lcl_time_hr+":00	"+flx_dwn_dff_ott_ttl+"			"+flx_dwn_drc_ott_ttl+"			"+dff_drc_rat)
end do; end loop over fl_idx

if(day.gt.31) then 
fl_out_nm="/data2/zender/arese/arese_951011_w_"+ss_alb_sng+"_trn_abs.nc"
system("/bin/rm -f "+fl_out_nm)
fl_out=addfile(fl_out_nm,"c")
fl_out@history="Created by ~/sw/ncl/tddr_abs.ncl"
fl_out@usage="File is intended as input to ~/sw/idl/mie.pro:trn_abs_gph()"
ss_alb=stringtofloat(ss_alb_sng)/100.0
ss_alb@long_name="Aerosol single scattering albedo at 0.525 um"
odxc_aer!0="odxc_aer"
odxc_aer@long_name="Aerosol extinction optical depth"
tddr_trn_mdl=trn_clm
tddr_trn_mdl!0="odxc_aer"
tddr_trn_mdl@long_name="Transmittance between aircraft in TDDR 500 nm channel"
tddr_abs_mdl=abs_clm
tddr_abs_mdl!0="odxc_aer"
tddr_trn_mdl@long_name="Absorptance between aircraft in TDDR 500 nm channel"
fl_out->ss_alb=ss_alb
fl_out->odxc_aer=odxc_aer
fl_out->tddr_trn_mdl=tddr_trn_mdl
fl_out->tddr_abs_mdl=tddr_abs_mdl
print("Wrote transmittance and absorptance to "+fl_out_nm)
end if ; endif day.gt.31
 
end ; end tddr_abs()


