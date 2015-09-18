; $Id$

; Purpose: offline analysis routines for looking at ACA data 
; These routines rely heavily on netCDF

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Running averager
; Inputs:
; data: ordinate array (e.g., spectral flux) containing data to be averaged
; crd_ctr: abscissa (coordinate) array (e.g., wavelength)
; crd_sz: coordinate-width array which contains the interval of each coordinate, 
; i.e., crd_sz(i)=crd_rgt(i)-crd_lft(i)
; avg_ntv: scalar specifying duration of running average in units of crd_ctr. 
; The original array is not altered.
; The returned array is a running average of the original. 
;
; Boundary technique: boundary points are only averaged over the available
; bandwidth to one side.
;
; NB: It would be nice if the user could supply an entire response
; function along with the bandwidth. For now, the response function
; is considered to be a square window whose FWFM is bandwidth.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function rnn_avg,data,crd_ctr,crd_sz=crd_sz,avg_ntv

ntvo2=avg_ntv/2.

; Dimension output running average array.
data_nbr=n_elements(data)
avg=data
wgt=data

; Create the crd_sz array if it was not passed in
; Assume edge sizes are identical to adjacent interior sizes
if n_elements(crd_sz) ne 0 then begin
	crd_sz=0.5*(shift(data,-1)-shift(data,1))
	crd_sz(0)=crd_sz(1)
	crd_sz(data_nbr-1)=crd_sz(data_nbr-2)
endif; endif auto-generated crd_sz

; Find points where special boundary averaging rules must be applied.
for idx=0,data_nbr-1 do begin
nbr_pt=1
avg(idx)=0.
wgt(idx)=0.
rgt_idx=idx
lft_idx=idx

while (abs(crd_ctr(rgt_idx)-crd_ctr(idx)) le ntvo2) do begin
	avg(idx)=avg(idx)+crd_sz(rgt_idx)*data(rgt_idx)
	wgt(idx)=wgt(idx)+crd_sz(rgt_idx)
	if rgt_idx lt data_nbr-1 then rgt_idx=rgt_idx+1 else goto,end_rgt
endwhile ; end loop over adding points to the right
end_rgt: foo=1

while (abs(crd_ctr(lft_idx)-crd_ctr(idx)) le ntvo2) do begin
	avg(idx)=avg(idx)+crd_sz(lft_idx)*data(lft_idx)
	wgt(idx)=wgt(idx)+crd_sz(lft_idx)
	if lft_idx gt 0 then lft_idx=lft_idx-1 else goto,end_lft
endwhile ; end loop over adding points to the left
end_lft: foo=1

; Correct for (intentionally) adding data(idx) twice to accumulating arrays. 
avg(idx)=avg(idx)-data(idx)*crd_sz(idx)
wgt(idx)=wgt(idx)-crd_sz(idx)

; Normalize the accumulated sum over the averaging interval.
avg(idx)=avg(idx)/wgt(idx)

endfor ; end loop over all data points
end_of_function: foo=1
return,avg
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Running averager
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure Cess
; compute the % and absolute amount of solar absorption in a 
; cloud between two given wavelengths
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro cess, $
	cld_lev=cld_lev, $
	srt_bnd_idx=srt_bnd_idx, $
	end_bnd_idx=end_bnd_idx, $
	fl_nm=fl_nm
;
; Where are the wavelengths of interest?
; Output from netget:
; bnd[1228] = 7.00035e-07, bnd[1228] = 7.00035e-07 meter
; bnd[1642] = 3.0075e-07, bnd[1642] = 3.0075e-07 meter
;
; Recall that these are the C and IDL indices, the
; Fortran indices are obtained by adding one to these.
;
if n_elements(cld_lev) eq 0 then cld_lev=13
if n_elements(srt_bnd_idx) eq 0 then srt_bnd_idx=1228
if n_elements(end_bnd_idx) eq 0 then end_bnd_idx=1642
if n_elements(fl_nm) eq 0 then fl_nm='/data/zender/aca/aca.100.nc'
;
print,'Gathering data from '+fl_nm
nc_id=ncdf_open(fl_nm)
;
; Get the dimension sizes
;
bnd_id=ncdf_dimid(nc_id,'bnd')
ncdf_diminq,nc_id,bnd_id,foo,nbr_bnd
levp_id=ncdf_dimid(nc_id,'levp')
ncdf_diminq,nc_id,levp_id,foo,levp_nbr
lev_id=ncdf_dimid(nc_id,'lev')
ncdf_diminq,nc_id,lev_id,foo,lev_nbr
pol_id=ncdf_dimid(nc_id,'pol')
ncdf_diminq,nc_id,pol_id,foo,nbr_pol
;
; Get the scalars
;
ncdf_varget,nc_id,'alb_sfc',alb_sfc
ncdf_varget,nc_id,'slr_zen_ngl_cos',slr_zen_ngl_cos
ncdf_varget,nc_id,'CWP_ttl',CWP_ttl
ncdf_varget,nc_id,'flx_bb_dwn_TOA',flx_bb_dwn_TOA
ncdf_varget,nc_id,'frc_ice_ttl',frc_ice_ttl
;
; Get the one-dimensional arrays
;
ncdf_varget,nc_id,'wvl_dlt',wvl_dlt
ncdf_varget,nc_id,'bnd',bnd
ncdf_varget,nc_id,'lev',lev
;
; Get the two-dimensional arrays
;
ncdf_varget,nc_id,'flx_spc_dwn_TOA',flx_spc_dwn_TOA
ncdf_varget,nc_id,'flx_spc_abs',flx_spc_abs
;
; End of NetCDF commands
;
ncdf_close,nc_id
;
; Truncate the data
;
;wvl_dlt=wvl_dlt(srt_bnd_idx:end_bnd_idx)
;flx_spc_dwn_TOA=flx_spc_dwn_TOA(srt_bnd_idx:end_bnd_idx)
;flx_spc_abs=flx_spc_abs(srt_bnd_idx:end_bnd_idx,cld_lev)
;
; Find the absolute absorbed and TOA flux
;
flx_bb_abs_cld=flx_spc_abs(*,cld_lev)*wvl_dlt
flx_bb_dwn_TOA=flx_spc_dwn_TOA*wvl_dlt
flx_int_abs_cld=flx_spc_abs(srt_bnd_idx:end_bnd_idx,cld_lev)*wvl_dlt(srt_bnd_idx:end_bnd_idx)
flx_int_dwn_TOA=flx_spc_dwn_TOA(srt_bnd_idx:end_bnd_idx)*wvl_dlt(srt_bnd_idx:end_bnd_idx)
;
; Total up the diagnostics
;
flx_bb_dwn_TOA=total(flx_bb_dwn_TOA)
flx_bb_abs_cld=total(flx_bb_abs_cld)
flx_int_dwn_TOA=total(flx_int_dwn_TOA)
flx_int_abs_cld=total(flx_int_abs_cld)
;
abs_int_cld=flx_int_abs_cld/flx_int_dwn_TOA
abs_bb_cld=flx_bb_abs_cld/flx_bb_dwn_TOA
;
print,'# of streams nbr_pol: '+auto_sng(nbr_pol,1)
print,'# of levels lev_nbr: '+auto_sng(lev_nbr,1)
print,'CCM cloud layer cld_lev: '+auto_sng(cld_lev,1)
print,'Midpt. pressure lev(cld_lev): '+auto_sng(.01*lev(cld_lev),1)+' mb'
print,'Surface albedo alb_sfc: '+auto_sng(alb_sfc,2)
print,'Cosine solar zenith angle slr_zen_ngl_cos: '+auto_sng(slr_zen_ngl_cos,2)
print,'Condensed Water Path CWP_ttl: '+auto_sng(1000.*CWP_ttl,1)+' g/m2'
print,'Fraction of ice in cloud frc_ice_ttl: '+auto_sng(100.*frc_ice_ttl,1)+' %'
print,''
print,'Total broadband TOA flux flx_bb_dwn_TOA: '+auto_sng(flx_bb_dwn_TOA,1)+' W/m2'
print,'Total broadband flux absorbed in layer flx_bb_abs_cld: '+auto_sng(flx_bb_abs_cld,1)+' W/m2'
print,'Total broadband absorptance in layer abs_bb_cld: '+auto_sng(100.*abs_bb_cld,2)+' %'
print,''
print,'Interval start wavelength bnd(srt_bnd_idx): '+auto_sng(1.0e6*bnd(srt_bnd_idx),3)+' microns'
print,'Interval end wavelength bnd(end_bnd_idx): '+auto_sng(1.0e6*bnd(end_bnd_idx),3)+' microns'
print,''
print,'Interval TOA flux flx_int_dwn_TOA: '+auto_sng(flx_int_dwn_TOA,1)+' W/m2'
print,'Interval flux absorbed in layer flx_int_abs_cld: '+auto_sng(flx_int_abs_cld,1)+' W/m2'
print,'Interval absorptance in layer abs_int_cld: '+auto_sng(100.*abs_int_cld,2)+' %'
;
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Cess
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure psiphi
; Compare the exact to the least-squares fit temperature dependence
; of a band
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro psiphi, $
	fl_nm=fl_nm
;
; psiphi,fl_nm='/data/zender/aca/H2_16O.nc'
; psiphi,fl_nm='/data/zender/aca/OH.nc'
;
if n_elements(fl_nm) eq 0 then fl_nm='/data/zender/aca/H2_16O.nc'
;
print,'Gathering data from '+fl_nm
nc_id=ncdf_open(fl_nm)
;
; Get the dimension sizes
;
t_id=ncdf_dimid(nc_id,'t')
ncdf_diminq,nc_id,t_id,foo,nbr_t
;
; Get the scalars
;
;ncdf_varget,nc_id,'alb_sfc',alb_sfc
;
; Get the one-dimensional arrays
;
ncdf_varget,nc_id,'t',t
ncdf_varget,nc_id,'phi_exact',phi_exact
ncdf_varget,nc_id,'psi_exact',psi_exact
ncdf_varget,nc_id,'psi_fit',psi_fit
ncdf_varget,nc_id,'phi_fit',phi_fit
;
; End of NetCDF commands
;
ncdf_close,nc_id
;
plot, $
	t, $
	phi_exact, $
	tit='!5 Exact !7U!8!Ik!N!5 and Least squares fit !7U!8(T)!5', $
	xtit='Temperature !8T!5 (K)', $
	ytit='!5Weak line limit sensitivity !7U!8(T)!5', $
	psym=1, $
	xstyle=1, $
	ystyle=1, $
	/ynozero, $
	thick=3.0, $
	charsize=2.0, $
	linestyle=0
;
oplot, $
	t, $
	phi_fit, $
	thick=3.0, $
	linestyle=0
;
ln_lgn_x1=.3
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=0.8
lgn_dy=0.1
lgn_y=lgn_y_top-indgen(6)*lgn_dy
;
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(0)+0.013,linestyle=0,psym=1,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013,linestyle=0,thick=3.0,/NORMAL
;
xyouts,txt_lgn_x,lgn_y(0),'!5Exact',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(1),'!5Least squares fit',size=txt_lgn_sz,/NORMAL
;
print,'Hit any key to continue, or q to quit ...'
junk=get_kbrd(1)
if junk eq 'q' then goto,end_of_procedure
;
plot, $
	t, $
	psi_exact, $
	tit='!5 Exact !7W!8!Ik!N!5 and Least squares fit !7W!8(T)!5', $
	xtit='Temperature !8T!5 (K)', $
	ytit='!5Strong line limit sensitivity !7W!8(T)!5', $
	psym=1, $
	xstyle=1, $
	ystyle=1, $
	/ynozero, $
	thick=3.0, $
	charsize=2.0, $
	linestyle=0
;
oplot, $
	t, $
	psi_fit, $
	thick=3.0, $
	linestyle=0
;
ln_lgn_x1=.3
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=0.8
lgn_dy=0.1
lgn_y=lgn_y_top-indgen(6)*lgn_dy
;
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(0)+0.013,linestyle=0,psym=1,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013,linestyle=0,thick=3.0,/NORMAL
;
xyouts,txt_lgn_x,lgn_y(0),'!5Exact',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(1),'!5Least squares fit',size=txt_lgn_sz,/NORMAL
;
end_of_procedure: foo=1
;
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End psiphi
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure ACA Original
; Scatterplot the TDDR channel intensity vs. the .5 microns intensity.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro aca_orig, $
	srt_idx=srt_idx, $
	end_idx=end_idx, $
	ncr_idx=ncr_idx, $
	chn_bnd=chn_bnd, $
	fld_nbr=fld_nbr, $
	top_int=top_int, $
	btm_int=btm_int, $
	fl_in=fl_in, $
	fl_cmp=fl_cmp, $
	fzz=fzz, $
	NTN=NTN, $
	FLX=FLX, $
	TRN=TRN, $
	CLR=CLR, $
	RTA=RTA, $
	TWO=TWO, $
	SWCF=SWCF, $
	ann=ann, $
	level=level
;
; aca_orig,srt_idx=0,end_idx=0,ncr_idx=25
; aca_orig,srt_idx=100,end_idx=100,ncr_idx=0
; aca_orig,RTA=1,FLX=0,NTN=0,fzz=0,TRN=0,fl_in='/data/zender/arese/mdl/951030_1200_arese_mdl_cld_aer.nc',CLR='/data/zender/arese/mdl/951030_1200_arese_mdl_clr_aer.nc',btm_int=93,top_int=26
; aca_orig,RTA=0,FLX=1,NTN=1,fzz=1,TRN=0,fl_in='/data/zender/aca/aca.100.nc',CLR='/data/zender/aca/aca.000.nc'
; aca_orig,RTA=0,FLX=1,NTN=1,fzz=0,TRN=0,TWO=1,srt_idx=12,end_idx=18,ncr_idx=1
; aca_orig,RTA=0,FLX=0,NTN=1,fzz=0,TRN=0,fl_in='/data/zender/aca/aca.ice.nc',fl_cmp='/data/zender/aca/aca.100.nc'
; aca_orig,fl_in='/data/zender/aca/950417_1940_nrel_dsrt_cld.nc',CLR='/data/zender/aca/950417_1940_nrel_dsrt_clr.nc',btm_int=18,top_int=18
; aca_orig,fl_in='/data/zender/aca/931021_2200_nrel_dsrt_cld.nc',CLR='/data/zender/aca/931021_2200_nrel_dsrt_clr.nc',btm_int=18,top_int=18
;
True=1
False=0
pi=3.141592654
;
if n_elements(srt_idx) eq 0 then srt_idx=0
if n_elements(end_idx) eq 0 then end_idx=0
if n_elements(ncr_idx) eq 0 then ncr_idx=25
if n_elements(chn_bnd) eq 0 then chn_bnd=800
if n_elements(fld_nbr) eq 0 then fld_nbr=2
if n_elements(top_int) eq 0 then top_int=7
if n_elements(btm_int) eq 0 then btm_int=15
if n_elements(fl_in) ne 0 then cmd_ln_fl=True else cmd_ln_fl=False
if n_elements(fl_cmp) ne 0 then CMP_FLG=True else CMP_FLG=False
if n_elements(fzz) eq 0 then fzz=False
if n_elements(NTN) eq 0 then NTN=True
if n_elements(FLX) eq 0 then FLX=True
if n_elements(TRN) eq 0 then TRN=False
if n_elements(CLR) eq 0 then CLR='/data/zender/aca/aca.000.nc'
if n_elements(RTA) eq 0 then RTA=False
if n_elements(TWO) eq 0 then TWO=False
if n_elements(SWCF) eq 0 then SWCF=True
if n_elements(ann) eq 0 then ann=False
if n_elements(vis_bnd) eq 0 then vis_bnd=1603
;if n_elements(foo) eq 0 then foo=1
;
nbr_sim=(end_idx-srt_idx)/ncr_idx+1
nbr_chn=11
print,'nbr_sim = ',nbr_sim
;
chn_rad_top=fltarr(nbr_chn,nbr_sim)
chn_rad_btm=fltarr(nbr_chn,nbr_sim)
;
foo={chn_sct, $
	idx:0, $
	TDDR_wvl_ctr:0.0, $
	aca_wvl_ctr:0.0, $
	aca_ctr_idx:0, $
	aca_wvl_min:0.0, $
	aca_min_idx:0, $
	aca_wvl_max:0.0, $
	aca_max_idx:0}
chans=replicate({chn_sct},nbr_chn)	
;
; Recall that these are the C and IDL indices, the
; Fortran indices are obtained by adding one to these.
;
;chans(0)={chn_sct,0,0.0,0.0,0,0.0,0,0.0,0}
;chans(0)={chn_sct,0,.500,.500,1602,0.0,0,0.0,0} ; visible channel
;chans(1)={chn_sct,0,.380,.380,1626,0.0,0,0.0,0}
;chans(2)={chn_sct,0,.412,.410,1620,0.0,0,0.0,0}
;chans(3)={chn_sct,0,0.675,0.674992,1281,0.0,0,0.0,0}
;chans(4)={chn_sct,0,0.862,0.861698,960,0.0,0,0.0,0}
;chans(5)={chn_sct,0,1.064,1.0644,739,0.0,0,0.0,0}
;chans(6)={chn_sct,0,1.640,1.64069,409,0.0,0,0.0,0}
;chans(7)={chn_sct,0,0.91,0.910332,898,0.0,0,0.0,0}
;chans(8)={chn_sct,0,1.18,1.17994,647,0.0,0,0.0,0}
;chans(9)={chn_sct,0,1.235,1.23533,609,0.0,0,0.0,0}
;chans(10)={chn_sct,0,2.2,2.20022,254,0.0,0,0.0,0}
chans(0)={chn_sct,0,.500,.500,1602,0.0,0,0.0,0} ; visible channel
chans(1)={chn_sct,0,0.862,0.861698,960,0.0,0,0.0,0}
chans(2)={chn_sct,0,1.064,1.0644,739,0.0,0,0.0,0}
chans(3)={chn_sct,0,1.249,1.24922,600,0.0,0,0.0,0}
chans(4)={chn_sct,0,1.501,1.50038,466,0.0,0,0.0,0}
chans(5)={chn_sct,0,1.651,1.65153,405,0.0,0,0.0,0}
chans(6)={chn_sct,0,1.750,1.74978,371,0.0,0,0.0,0}
;
chans.idx=indgen(nbr_chn)
;
; Before the time loop, gather some useful information that
; does not change (too much!) with each passing file. 
;
if CMP_FLG then begin
print,'Gathering CMP_FLG data from '+fl_cmp
nc_id=ncdf_open(fl_cmp)
levp_id=ncdf_dimid(nc_id,'levp')
ncdf_diminq,nc_id,levp_id,foo,levp_nbr
ncdf_varget,nc_id,'ntn_spc_aa_ndr',ntn_spc_aa_ndr
ncdf_varget,nc_id,'ntn_spc_aa_zen',ntn_spc_aa_zen
;ncdf_varget,nc_id,'foo',foo_clr
ncdf_close,nc_id
ntn_spc_aa_ndr2=ntn_spc_aa_ndr(*,btm_int)
ntn_spc_aa_zen2=ntn_spc_aa_zen(*,top_int)
endif; endif CMP_FLG
;
; Loop over the requested simulations
;
end_idx=long(end_idx)
srt_idx=long(srt_idx)
ncr_idx=long(ncr_idx)
;
sim=0
for idx=srt_idx,end_idx,ncr_idx do begin
;
; Construct the file name
;
if TWO then begin
idx_sng=string(format='(I2.2)',idx)
endif else begin
idx_sng=string(format='(I3.3)',idx)
endelse
fl_nm='aca.'+idx_sng+'.nc'
;
; Read in binary data from the NetCDF file 
;
fl_path='/data/zender/aca'
fl_nm=fl_path+'/'+fl_nm
if TWO then clr=fl_path+'/'+'aca.'+idx_sng+'.clr.nc'
if cmd_ln_fl then fl_nm=fl_in
print,'Gathering data from '+fl_nm
nc_id=ncdf_open(fl_nm)
;
; Get the dimension sizes
;
bnd_id=ncdf_dimid(nc_id,'bnd')
ncdf_diminq,nc_id,bnd_id,foo,nbr_bnd
levp_id=ncdf_dimid(nc_id,'levp')
ncdf_diminq,nc_id,levp_id,foo,levp_nbr
;
; Get the two-dimensional arrays
;
ncdf_varget,nc_id,'ntn_spc_aa_ndr',ntn_spc_aa_ndr
ncdf_varget,nc_id,'ntn_spc_aa_zen',ntn_spc_aa_zen
ncdf_varget,nc_id,'flx_spc_dwn',flx_spc_dwn
ncdf_varget,nc_id,'flx_spc_dwn_TOA',flx_spc_dwn_TOA
ncdf_varget,nc_id,'flx_spc_dwn_sfc',flx_spc_dwn_sfc
ncdf_varget,nc_id,'flx_bb_dwn_TOA',flx_bb_dwn_TOA
ncdf_varget,nc_id,'flx_bb_dwn_sfc',flx_bb_dwn_sfc
ncdf_varget,nc_id,'flx_bb_dwn_drc',flx_bb_dwn_drc
ncdf_varget,nc_id,'flx_bb_dwn_dff',flx_bb_dwn_dff
ncdf_varget,nc_id,'flx_bb_dwn',flx_bb_dwn
ncdf_varget,nc_id,'flx_bb_up',flx_bb_up
ncdf_varget,nc_id,'wvl_ctr',wvl_ctr
ncdf_varget,nc_id,'wvl_dlt',wvl_dlt
ncdf_varget,nc_id,'trn_spc_atm_ttl',trn_spc_atm_ttl
ncdf_varget,nc_id,'rfl_spc_SAS',rfl_spc_SAS
ncdf_varget,nc_id,'abs_spc_atm',abs_spc_atm
ncdf_varget,nc_id,'CWP_ttl',CWP_ttl
ncdf_varget,nc_id,'frc_ice_ttl',frc_ice_ttl
ncdf_varget,nc_id,'alb_sfc',alb_sfc
ncdf_attget,nc_id,'prf_sng',prf_sng,/GLOBAL
ncdf_varget,nc_id,'lcl_time_hr',lcl_time_hr
ncdf_varget,nc_id,'lcl_yr_day',lcl_yr_day
ncdf_varget,nc_id,'z_cld_btm',z_cld_btm
ncdf_varget,nc_id,'lat_dgr',lat_dgr
ncdf_varget,nc_id,'z_cld_thick',z_cld_thick
ncdf_varget,nc_id,'slr_zen_ngl_cos',slr_zen_ngl_cos
ncdf_varget,nc_id,'flx_bb_abs_atm',flx_bb_abs_atm
ncdf_varget,nc_id,'flx_bb_abs_sfc',flx_bb_abs_sfc
ncdf_varget,nc_id,'rfl_bb_SAS',rfl_bb_SAS
ncdf_varget,nc_id,'flx_bb_net',flx_bb_net
ncdf_varget,nc_id,'opt_dep_ext_clm_ttl',opt_dep_ext_clm_ttl
if TRN then begin
ncdf_varget,nc_id,'trn_spc_atm_CO2',trn_spc_atm_CO2
ncdf_varget,nc_id,'trn_spc_atm_H2O',trn_spc_atm_H2O
ncdf_varget,nc_id,'trn_spc_atm_O2',trn_spc_atm_O2
ncdf_varget,nc_id,'trn_spc_atm_O3',trn_spc_atm_O3
ncdf_varget,nc_id,'trn_spc_atm_Ray',trn_spc_atm_Ray
ncdf_varget,nc_id,'trn_spc_atm_ice',trn_spc_atm_ice
ncdf_varget,nc_id,'trn_spc_atm_liq',trn_spc_atm_liq
ncdf_varget,nc_id,'trn_spc_atm_ttl',trn_spc_atm_ttl
endif; endif plotting individual transmissions
;
; End of NetCDF commands
;
ncdf_close,nc_id
;
; The clear sky file might change with each file
; like when we're varying the solar zenith angle, so this needs
; to go in the main time loop.
;
if SWCF then begin
print,'Gathering clear sky data from '+CLR
nc_id=ncdf_open(CLR)
levp_id=ncdf_dimid(nc_id,'levp')
ncdf_diminq,nc_id,levp_id,foo,levp_nbr
ncdf_varget,nc_id,'flx_bb_net',flx_bb_net_clr
;ncdf_varget,nc_id,'foo',foo_clr
ncdf_close,nc_id
flx_bb_net_clr_TOA=flx_bb_net_clr(0)
flx_bb_net_clr_sfc=flx_bb_net_clr(levp_nbr-1)
flx_bb_net_clr_top=flx_bb_net_clr(top_int)
flx_bb_net_clr_btm=flx_bb_net_clr(btm_int)
endif; endif SWCF
;
;print,'levp_nbr = ',levp_nbr
;print,'nbr_bnd = ',nbr_bnd
;
flx_spc_dwn_top=flx_spc_dwn(*,top_int)
flx_spc_dwn_btm=flx_spc_dwn(*,btm_int)
flx_bb_net_TOA=flx_bb_net(0)
flx_bb_net_sfc=flx_bb_net(levp_nbr-1)
flx_bb_net_top=flx_bb_net(top_int)
flx_bb_net_btm=flx_bb_net(btm_int)
flx_bb_dwn_top=flx_bb_dwn(top_int)
flx_bb_dwn_btm=flx_bb_dwn(btm_int)
flx_bb_dwn_drc_top=flx_bb_dwn_drc(top_int)
flx_bb_dwn_drc_btm=flx_bb_dwn_drc(btm_int)
flx_bb_dwn_dff_top=flx_bb_dwn_dff(top_int)
flx_bb_dwn_dff_btm=flx_bb_dwn_dff(btm_int)
flx_bb_up_top=flx_bb_up(top_int)
flx_bb_up_btm=flx_bb_up(btm_int)
flx_bb_abs_cld=(flx_bb_net_top-flx_bb_net_btm)/flx_bb_dwn_TOA
opt_dep_ext_cld_vis=0.01*(round(100.*(opt_dep_ext_clm_ttl(1602)-.1472))) ; subtract clear sky value
;
ntn_spc_aa_ndr=ntn_spc_aa_ndr(*,btm_int)
ntn_spc_aa_zen=ntn_spc_aa_zen(*,top_int)
;
if FZZ then begin
print,'Averaging Intensities...'
ntn_spc_aa_ndr=rnn_avg(ntn_spc_aa_ndr,wvl_ctr,crd_sz=wvl_dlt,1.0e-8)
ntn_spc_aa_zen=rnn_avg(ntn_spc_aa_zen,wvl_ctr,crd_sz=wvl_dlt,1.0e-8)
;foo=rnn_avg(foo,wvl_ctr,crd_sz=wvl_dlt,1.0e-8)
if FLX then begin
	print,'Averaging FLXes...'
	flx_spc_dwn_TOA=rnn_avg(flx_spc_dwn_TOA,wvl_ctr,crd_sz=wvl_dlt,1.0e-8)
	flx_spc_dwn_sfc=rnn_avg(flx_spc_dwn_sfc,wvl_ctr,crd_sz=wvl_dlt,1.0e-8)
	flx_spc_dwn_top=rnn_avg(flx_spc_dwn_top,wvl_ctr,crd_sz=wvl_dlt,1.0e-8)
	flx_spc_dwn_btm=rnn_avg(flx_spc_dwn_btm,wvl_ctr,crd_sz=wvl_dlt,1.0e-8)
endif; endif FLX averaging
if RTA then begin
	print,'Averaging RTA...'
	trn_spc_atm_ttl=rnn_avg(trn_spc_atm_ttl,wvl_ctr,crd_sz=wvl_dlt,1.0e-8)
	rfl_spc_SAS=rnn_avg(rfl_spc_SAS,wvl_ctr,crd_sz=wvl_dlt,1.0e-8)
	abs_spc_atm=rnn_avg(abs_spc_atm,wvl_ctr,crd_sz=wvl_dlt,1.0e-8)
endif; endif RTA averaging
if CMP_FLG then begin
	print,'Averaging CMP_FLG...'
	ntn_spc_aa_zen2=rnn_avg(ntn_spc_aa_zen2,wvl_ctr,crd_sz=wvl_dlt,1.0e-8)
	ntn_spc_aa_ndr2=rnn_avg(ntn_spc_aa_ndr2,wvl_ctr,crd_sz=wvl_dlt,1.0e-8)
endif; endif CMP_FLG averaging
if TRN then begin
	print,'Averaging TRNs...'
	trn_spc_atm_CO2=rnn_avg(trn_spc_atm_CO2,wvl_ctr,crd_sz=wvl_dlt,1.0e-8)
	trn_spc_atm_H2O=rnn_avg(trn_spc_atm_H2O,wvl_ctr,crd_sz=wvl_dlt,1.0e-8)
	trn_spc_atm_O2=rnn_avg(trn_spc_atm_O2,wvl_ctr,crd_sz=wvl_dlt,1.0e-8)
	trn_spc_atm_O3=rnn_avg(trn_spc_atm_O3,wvl_ctr,crd_sz=wvl_dlt,1.0e-8)
	trn_spc_atm_Ray=rnn_avg(trn_spc_atm_Ray,wvl_ctr,crd_sz=wvl_dlt,1.0e-8)
	trn_spc_atm_ice=rnn_avg(trn_spc_atm_ice,wvl_ctr,crd_sz=wvl_dlt,1.0e-8)
	trn_spc_atm_liq=rnn_avg(trn_spc_atm_liq,wvl_ctr,crd_sz=wvl_dlt,1.0e-8)
	trn_spc_atm_ttl=rnn_avg(trn_spc_atm_ttl,wvl_ctr,crd_sz=wvl_dlt,1.0e-8)
endif; end TRN averaging
endif ; end if applying 10 nm window
;
ntn_spc_aa_ndr=ntn_spc_aa_ndr*1.0e-6 ; W/m2/m/sr --> W/m2/micron/sr
ntn_spc_aa_zen=ntn_spc_aa_zen*1.0e-6 ; W/m2/m/sr --> W/m2/micron/sr
flx_spc_dwn_TOA=flx_spc_dwn_TOA*1.0e-6 ; W/m2/m --> W/m2/micron
flx_spc_dwn_sfc=flx_spc_dwn_sfc*1.0e-6 ; W/m2/m --> W/m2/micron
flx_spc_dwn_top=flx_spc_dwn_top*1.0e-6 ; W/m2/m --> W/m2/micron
flx_spc_dwn_btm=flx_spc_dwn_btm*1.0e-6 ; W/m2/m --> W/m2/micron
;foo=foo*1.0e-6 ; W/m2/m --> W/m2/micron
if CMP_FLG then begin
ntn_spc_aa_ndr2=ntn_spc_aa_ndr2*1.0e-6 ; W/m2/m/sr --> W/m2/micron/sr
ntn_spc_aa_zen2=ntn_spc_aa_zen2*1.0e-6 ; W/m2/m/sr --> W/m2/micron/sr
endif; endif CMP_FLG
;
wvl_ctr=wvl_ctr*1.0e6 ; m --> micron
wvl_dlt=wvl_dlt*1.0e6 ; m --> micron
;
for chan=0,nbr_chn-1 do begin
	chn_rad_top(chan,sim)= $	
	ntn_spc_aa_zen(chans(chan).aca_ctr_idx)
	chn_rad_btm(chan,sim)= $	
	ntn_spc_aa_ndr(chans(chan).aca_ctr_idx)
endfor ; end loop over current channels
;
if CWP_ttl eq 0. then begin
opt_dep_ext_cld_vis=0.
endif; endif there is no cloud
if lat_dgr ge 0. then hemi_sng='N' else hemi_sng='S'
prf_sng=strtrim(string(prf_sng),2)
info_1_sng=prf_sng+', Lat = '+auto_sng(lat_dgr,2)+'!E!12_!5!N'+hemi_sng+', Day = '+auto_sng(fix(lcl_yr_day),0)+', Hour = '+auto_sng(lcl_time_hr,2)+', !7H!5 = '+auto_sng(180.*acos(slr_zen_ngl_cos)/pi,2)+'!E!12_!5!N, cos(!7H!5) = '+auto_sng(slr_zen_ngl_cos,2)+', !8A!5!Isfc!N = '+auto_sng(alb_sfc,2)
info_2_sng='!8CWP!5 = '+auto_sng(round(CWP_ttl*1000.),0)+' g m!E-2!N, Base = '+auto_sng(z_cld_btm/1000.0,1)+' km, Thick = '+auto_sng(z_cld_thick/1000.0,1)+' km, '+' !8f!5!Iice!N = '+auto_sng(round(100.*frc_ice_ttl),0)+' %'
info_3_sng='!8A!5!ITOA!N =  '+auto_sng(rfl_bb_SAS,2)+', !7s!5!Icld!N(.5 !7l!5m) = '+auto_sng(opt_dep_ext_cld_vis,2)+', !7j!5(cld) = '+auto_sng(flx_bb_abs_cld,2)
info_3_sng=info_3_sng+', !8F!5!95!5!N(TOA) = '+auto_sng(round(flx_bb_dwn_TOA),0)+', !8F!5!95!5!N(sfc) = '+auto_sng(round(flx_bb_dwn_sfc),0)+', !8F!5!Iabs!N(atm) = '+auto_sng(round(flx_bb_abs_atm),0)+', !8F!5!Iabs!N(sfc) = '+auto_sng(round(flx_bb_abs_sfc),0)+', !8F!5!Iabs!N(cld) = '+auto_sng(round(flx_bb_net_top-flx_bb_net_btm),0)
;
info_4_sng='13 km: !8F!5!95!5!N = '+auto_sng(round(flx_bb_dwn_top),0)+', !8F!5!S!E!95!5!N!R!Idir!N = '+auto_sng(round(flx_bb_dwn_drc_top),0)+', !8F!5!97!5!N = '+auto_sng(round(flx_bb_up_top),0)+', !8F!5!INet!N = '+auto_sng(round(flx_bb_net_top),0)
info_4_sng=info_4_sng+'. 1 km: !8F!5!95!5!N = '+auto_sng(round(flx_bb_dwn_btm),0)+', !8F!5!S!E!95!5!N!R!Idir!N = '+auto_sng(round(.01*(round(100.*flx_bb_dwn_drc_btm))),0)+', !8F!5!97!5!N = '+auto_sng(round(flx_bb_up_btm),0)+', !8F!5!INet!N = '+auto_sng(round(flx_bb_net_btm),0)
if SWCF then begin
SWCF_TOA=flx_bb_net_clr_TOA-flx_bb_net_TOA
SWCF_sfc=flx_bb_net_clr_sfc-flx_bb_net_sfc
SWCF_top=flx_bb_net_clr_top-flx_bb_net_top
SWCF_btm=flx_bb_net_clr_btm-flx_bb_net_btm
if SWCF_TOA ne 0. then SWCF_rat_clm=SWCF_sfc/SWCF_TOA else SWCF_rat_clm=1.
if SWCF_top ne 0. then SWCF_rat_cld=SWCF_btm/SWCF_top else SWCF_rat_cld=1.
endif else begin; endif SWCF
SWCF_TOA=0.
SWCF_sfc=0.
SWCF_top=0.
SWCF_btm=0.
SWCF_rat_clm=1.
SWCF_rat_cld=1.
endelse; endif no clear sky file
;if not CMP_FLG then begin
info_2_sng=info_2_sng+'. CF(TOA) = '+auto_sng(round(SWCF_TOA),0)+', CF(sfc) = '+auto_sng(round(SWCF_sfc),0)+', CFR(col) = '+auto_sng(SWCF_rat_clm,2)+', CFR(cld) = '+auto_sng(SWCF_rat_cld,2)
;endif; endif not CMP_FLG
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Start cloud radiance graph
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if NTN then begin
;
ord_max=max([ntn_spc_aa_ndr,ntn_spc_aa_zen])
;
plot, $
	wvl_ctr, $
	ntn_spc_aa_zen, $
	tit='!5Radiance at Flight Levels', $
	xtit='!5Wavelength !7k!5 (!7l!5m)', $
	ytit='!5Radiance !8I!I!7k!5!N (!5W m!E-2!N !7l!5m!E-1!N sr!E-1!N)', $
	xrange=[.2,2.5], $
	yrange=[0.0,ord_max], $
	xmargin=[5.5,3.], $	;[10,3] is [left,right] default
	ymargin=[3.0,1.0], $	;[4,2] is [bottom,top] default
	xstyle=1, $
	ystyle=1, $
	thick=1.0, $
	charsize=2.0, $
	linestyle=0
;
oplot,	$
	wvl_ctr, $
	ntn_spc_aa_ndr, $
	thick=1.0, $
	linestyle=2
;
; Loop over the current and projected channels
;
for chan=0,6 do begin
oplot,	$
	[chans(chan).TDDR_wvl_ctr,chans(chan).TDDR_wvl_ctr], $
	[0.0,ord_max], $
	thick=1.0, $
	linestyle=1
endfor ; end loop over current channels
;
for chan=7,10 do begin
oplot,	$
	[chans(chan).TDDR_wvl_ctr,chans(chan).TDDR_wvl_ctr], $
	[0.0,ord_max], $
	thick=1.0, $
	linestyle=3
endfor ; end loop over new channels
;
ln_lgn_x1=0.6
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=0.8
lgn_dy=0.1
lgn_y=lgn_y_top-indgen(6)*lgn_dy
;
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(0)+0.013,linestyle=0,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013,linestyle=2,thick=3.0,/NORMAL
;plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(2)+0.013,linestyle=1,thick=3.0,/NORMAL
;plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(3)+0.013,linestyle=3,thick=3.0,/NORMAL
;
xyouts,txt_lgn_x,lgn_y(0),'!8I!S!I!7k!R!E!97!5!N !8z!5 = 13 km',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(1),'!8I!S!I!7k!R!E!95!5!N !8z!5 = 1 km',size=txt_lgn_sz,/NORMAL
;xyouts,txt_lgn_x,lgn_y(2),'Existing Channel',size=txt_lgn_sz,/NORMAL
;xyouts,txt_lgn_x,lgn_y(3),'Candidate Channel',size=txt_lgn_sz,/NORMAL
;
;xyouts,.1,0.93,info_sng,alignment=0.0,orientation=0.0,size=0.65,/NORMAL
;xyouts,0.99,.12,info_3_sng,alignment=0.0,orientation=90.0,size=0.65,/NORMAL
xyouts,1.0,0.95,info_1_sng,alignment=0.0,orientation=270.0,size=0.65,/NORMAL
xyouts,0.985,0.95,info_2_sng,alignment=0.0,orientation=270.0,size=0.65,/NORMAL
xyouts,0.97,0.95,info_3_sng,alignment=0.0,orientation=270.0,size=0.65,/NORMAL
xyouts,0.955,0.95,info_4_sng,alignment=0.0,orientation=270.0,size=0.65,/NORMAL
;
if ann then ntt
;
print,'Hit any key to continue, or q to quit ...'
junk=get_kbrd(1)
if junk eq 'q' then goto,end_of_procedure
;
endif; endif NTN
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End cloud radiance graph
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Start cloud compare graph
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if CMP_FLG then begin
;
; In this case the printed radiance graphs will compare the radiances
; of two different cloud scenarios at the same level, instead of 
; the normal, which is to show both up and downwelling radiances at
; different levels for the same cloud scenario.
;
; Replace the '_zen' string with '_ndr' and switch the commented with
; the uncommented labels to change which level for the radiances.
; Also remember to change the descriptive text by the line labels to
; discriminate between the cases, e.g., add 'Liquid' to one and 'Ice'
; the other.
;
ord_max=max([ntn_spc_aa_zen,ntn_spc_aa_zen2])
;
plot, $
	wvl_ctr, $
	ntn_spc_aa_zen, $
	xtit='!5Wavelength !7k!5 (!7l!5m)', $
	tit='!5Radiance at !8z!5 = 13 km', $
	ytit='!5Radiance !8I!S!E!97!5!N!R!I!7k!5!N (!5W m!E-2!N !7l!5m!E-1!N sr!E-1!N)', $
;	tit='!5Radiance at !8z!5 = 1 km', $
;	ytit='!5Radiance !8I!S!E!95!5!N!R!I!7k!5!N (!5W m!E-2!N !7l!5m!E-1!N sr!E-1!N)', $
	xrange=[.2,2.5], $
	yrange=[0.0,ord_max], $
	xmargin=[5.5,3.], $	;[10,3] is [left,right] default
	ymargin=[3.0,1.0], $	;[4,2] is [bottom,top] default
	xstyle=1, $
	ystyle=1, $
	thick=1.0, $
	charsize=2.0, $
	linestyle=0
;
oplot,	$
	wvl_ctr, $
	ntn_spc_aa_zen2, $
	thick=1.0, $
	linestyle=2
;
; Loop over the current and projected channels
;
for chan=0,6 do begin
oplot,	$
	[chans(chan).TDDR_wvl_ctr,chans(chan).TDDR_wvl_ctr], $
	[0.0,ord_max], $
	thick=1.0, $
	linestyle=1
endfor ; end loop over current channels
;
for chan=7,10 do begin
oplot,	$
	[chans(chan).TDDR_wvl_ctr,chans(chan).TDDR_wvl_ctr], $
	[0.0,ord_max], $
	thick=1.0, $
	linestyle=3
endfor ; end loop over new channels
;
ln_lgn_x1=0.6
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=0.8
lgn_dy=0.1
lgn_y=lgn_y_top-indgen(6)*lgn_dy
;
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(0)+0.013,linestyle=0,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013,linestyle=2,thick=3.0,/NORMAL
;plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(2)+0.013,linestyle=1,thick=3.0,/NORMAL
;plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(3)+0.013,linestyle=3,thick=3.0,/NORMAL
;
xyouts,txt_lgn_x,lgn_y(0),'!8I!S!I!7k!R!E!97!5!N Ice',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(1),'!8I!S!I!7k!R!E!97!5!N Liquid',size=txt_lgn_sz,/NORMAL
;xyouts,txt_lgn_x,lgn_y(0),'!8I!S!I!7k!R!E!95!5!N Ice',size=txt_lgn_sz,/NORMAL
;xyouts,txt_lgn_x,lgn_y(1),'!8I!S!I!7k!R!E!95!5!N Liquid',size=txt_lgn_sz,/NORMAL
;xyouts,txt_lgn_x,lgn_y(2),'Existing Channel',size=txt_lgn_sz,/NORMAL
;xyouts,txt_lgn_x,lgn_y(3),'Candidate Channel',size=txt_lgn_sz,/NORMAL
;
;xyouts,.1,0.93,info_sng,alignment=0.0,orientation=0.0,size=0.65,/NORMAL
;xyouts,0.99,.12,info_3_sng,alignment=0.0,orientation=90.0,size=0.65,/NORMAL
xyouts,1.0,0.95,info_1_sng,alignment=0.0,orientation=270.0,size=0.65,/NORMAL
xyouts,0.985,0.95,info_2_sng,alignment=0.0,orientation=270.0,size=0.65,/NORMAL
xyouts,0.97,0.95,info_3_sng,alignment=0.0,orientation=270.0,size=0.65,/NORMAL
xyouts,0.955,0.95,info_4_sng,alignment=0.0,orientation=270.0,size=0.65,/NORMAL
;
if ann then ntt
;
print,'Hit any key to continue, or q to quit ...'
junk=get_kbrd(1)
if junk eq 'q' then goto,end_of_procedure
;
endif; endif CMP_FLG
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End cloud compare graph
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Start cloud flux graph
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if FLX then begin
;
ord_max=max(flx_spc_dwn_TOA)
;
plot, $
	wvl_ctr, $
	flx_spc_dwn_TOA, $
	tit='!5Flux at TOA and Flight Levels', $
	xtit='!5Wavelength !7k!5 (!7l!5m)', $
	ytit='!5Flux !8F!I!7k!5!N (!5W m!E-2!N !7l!5m!E-1!N)', $
	xrange=[.2,2.5], $
	yrange=[0.0,ord_max], $
	xmargin=[5.5,3.], $	;[10,3] is [left,right] default
	ymargin=[3.0,1.0], $	;[4,2] is [bottom,top] default
	xstyle=1, $
	ystyle=1, $
	thick=1.0, $
	charsize=2.0, $
	linestyle=0
;
oplot,	$
	wvl_ctr, $
	flx_spc_dwn_top, $
	thick=1.0, $
	linestyle=0
;
oplot,	$
	wvl_ctr, $
	flx_spc_dwn_btm, $
	thick=1.0, $
	linestyle=2
;
; Loop over the current and projected channels
;
for chan=0,6 do begin
oplot,	$
	[chans(chan).TDDR_wvl_ctr,chans(chan).TDDR_wvl_ctr], $
	[0.0,ord_max], $
	thick=1.0, $
	linestyle=1
endfor ; end loop over current channels
;
for chan=7,10 do begin
oplot,	$
	[chans(chan).TDDR_wvl_ctr,chans(chan).TDDR_wvl_ctr], $
	[0.0,ord_max], $
	thick=1.0, $
	linestyle=3
endfor ; end loop over new channels
;
ln_lgn_x1=0.6
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=0.8
lgn_dy=0.1
lgn_y=lgn_y_top-indgen(6)*lgn_dy
;
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(0)+0.013,linestyle=0,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013,linestyle=0,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(2)+0.013,linestyle=2,thick=3.0,/NORMAL
;plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(3)+0.013,linestyle=1,thick=3.0,/NORMAL
;plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(4)+0.013,linestyle=3,thick=3.0,/NORMAL
;
xyouts,txt_lgn_x,lgn_y(0),'!8F!S!I!7k!R!E!95!5!N TOA',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(1),'!8F!S!I!7k!R!E!95!5!N 13 km',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(2),'!8F!S!I!7k!R!E!95!5!N 1 km',size=txt_lgn_sz,/NORMAL
;xyouts,txt_lgn_x,lgn_y(3),'Existing Channel',size=txt_lgn_sz,/NORMAL
;xyouts,txt_lgn_x,lgn_y(4),'Candidate Channel',size=txt_lgn_sz,/NORMAL
;
xyouts,1.0,0.95,info_1_sng,alignment=0.0,orientation=270.0,size=0.65,/NORMAL
xyouts,0.985,0.95,info_2_sng,alignment=0.0,orientation=270.0,size=0.65,/NORMAL
xyouts,0.97,0.95,info_3_sng,alignment=0.0,orientation=270.0,size=0.65,/NORMAL
xyouts,0.955,0.95,info_4_sng,alignment=0.0,orientation=270.0,size=0.65,/NORMAL
;
if ann then ntt
;
print,'Hit any key to continue, or q to quit ...'
junk=get_kbrd(1)
if junk eq 'q' then goto,end_of_procedure
;
endif; endif FLX
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End cloud flux graph
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Start transmittance graph
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if RTA then begin
;
ord_max=1.
;
plot, $
	wvl_ctr, $
	trn_spc_atm_ttl, $
	tit='!5Transmittance of Atmosphere', $
	xtit='!5Wavelength !7k!5 (!7l!5m)', $
	ytit='!5Transmittance !8T!I!7k!5!N', $
	xrange=[.2,2.5], $
;	yrange=[.2,2.5], $
	xmargin=[5.5,3.], $	;[10,3] is [left,right] default
	ymargin=[3.0,1.0], $	;[4,2] is [bottom,top] default
	xstyle=0, $
	ystyle=0, $
	thick=1.0, $
	charsize=2.0, $
	linestyle=0
;
; Loop over the current and projected channels
;
for chan=0,6 do begin
oplot,	$
	[chans(chan).TDDR_wvl_ctr,chans(chan).TDDR_wvl_ctr], $
	[0.0,ord_max], $
	thick=1.0, $
	linestyle=1
endfor ; end loop over current channels
;
for chan=7,10 do begin
oplot,	$
	[chans(chan).TDDR_wvl_ctr,chans(chan).TDDR_wvl_ctr], $
	[0.0,ord_max], $
	thick=1.0, $
	linestyle=3
endfor ; end loop over new channels
;
ln_lgn_x1=0.6
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=0.8
lgn_dy=0.1
lgn_y=lgn_y_top-indgen(6)*lgn_dy
;
;plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(0)+0.013,linestyle=1,thick=3.0,/NORMAL
;plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013,linestyle=3,thick=3.0,/NORMAL
;xyouts,txt_lgn_x,lgn_y(0),'Existing Channel',size=txt_lgn_sz,/NORMAL
;xyouts,txt_lgn_x,lgn_y(1),'Candidate Channel',size=txt_lgn_sz,/NORMAL
;
xyouts,1.0,0.95,info_1_sng,alignment=0.0,orientation=270.0,size=0.65,/NORMAL
xyouts,0.985,0.95,info_2_sng,alignment=0.0,orientation=270.0,size=0.65,/NORMAL
xyouts,0.97,0.95,info_3_sng,alignment=0.0,orientation=270.0,size=0.65,/NORMAL
xyouts,0.955,0.95,info_4_sng,alignment=0.0,orientation=270.0,size=0.65,/NORMAL
;
if ann then ntt
;
print,'Hit any key to continue, or q to quit ...'
junk=get_kbrd(1)
if junk eq 'q' then goto,end_of_procedure
;
endif; endif RTA
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of transmittance graph
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Start absorptance graph
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if RTA then begin
;
ord_max=1.
;
plot, $
	wvl_ctr, $
	abs_spc_atm, $
	tit='!5Absorptance of Atmosphere', $
	xtit='!5Wavelength !7k!5 (!7l!5m)', $
	ytit='!5Absorptance !8A!I!7k!5!N', $
	xrange=[.2,2.5], $
	yrange=[0.0,1.0], $
	xmargin=[5.5,3.], $	;[10,3] is [left,right] default
	ymargin=[3.0,1.0], $	;[4,2] is [bottom,top] default
	xstyle=0, $
	ystyle=0, $
	thick=1.0, $
	charsize=2.0, $
	linestyle=0
;
; Loop over the current and projected channels
;
for chan=0,6 do begin
oplot,	$
	[chans(chan).TDDR_wvl_ctr,chans(chan).TDDR_wvl_ctr], $
	[0.0,ord_max], $
	thick=1.0, $
	linestyle=1
endfor ; end loop over current channels
;
for chan=7,10 do begin
oplot,	$
	[chans(chan).TDDR_wvl_ctr,chans(chan).TDDR_wvl_ctr], $
	[0.0,ord_max], $
	thick=1.0, $
	linestyle=3
endfor ; end loop over new channels
;
ln_lgn_x1=0.6
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=0.8
lgn_dy=0.1
lgn_y=lgn_y_top-indgen(6)*lgn_dy
;
;plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(0)+0.013,linestyle=1,thick=3.0,/NORMAL
;plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013,linestyle=3,thick=3.0,/NORMAL
;xyouts,txt_lgn_x,lgn_y(0),'Existing Channel',size=txt_lgn_sz,/NORMAL
;xyouts,txt_lgn_x,lgn_y(1),'Candidate Channel',size=txt_lgn_sz,/NORMAL
;
xyouts,1.0,0.95,info_1_sng,alignment=0.0,orientation=270.0,size=0.65,/NORMAL
xyouts,0.985,0.95,info_2_sng,alignment=0.0,orientation=270.0,size=0.65,/NORMAL
xyouts,0.97,0.95,info_3_sng,alignment=0.0,orientation=270.0,size=0.65,/NORMAL
xyouts,0.955,0.95,info_4_sng,alignment=0.0,orientation=270.0,size=0.65,/NORMAL
;
if ann then ntt
;
print,'Hit any key to continue, or q to quit ...'
junk=get_kbrd(1)
if junk eq 'q' then goto,end_of_procedure
;
endif; endif RTA
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of absorptance graph
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Start reflectance graph
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if RTA then begin
;
ord_max=1.
;
plot, $
	wvl_ctr, $
	rfl_spc_SAS, $
	tit='!5Reflectance of Surface-Atmosphere System', $
	xtit='!5Wavelength !7k!5 (!7l!5m)', $
	ytit='!5Reflectance !8R!I!7k!5!N', $
	xrange=[.2,2.5], $
	yrange=[0.0,1.0], $
	xmargin=[5.5,3.], $	;[10,3] is [left,right] default
	ymargin=[3.0,1.0], $	;[4,2] is [bottom,top] default
	xstyle=0, $
	ystyle=0, $
	thick=1.0, $
	charsize=2.0, $
	linestyle=0
;
; Loop over the current and projected channels
;
for chan=0,6 do begin
oplot,	$
	[chans(chan).TDDR_wvl_ctr,chans(chan).TDDR_wvl_ctr], $
	[0.0,ord_max], $
	thick=1.0, $
	linestyle=1
endfor ; end loop over current channels
;
for chan=7,10 do begin
oplot,	$
	[chans(chan).TDDR_wvl_ctr,chans(chan).TDDR_wvl_ctr], $
	[0.0,ord_max], $
	thick=1.0, $
	linestyle=3
endfor ; end loop over new channels
;
ln_lgn_x1=0.6
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=0.8
lgn_dy=0.1
lgn_y=lgn_y_top-indgen(6)*lgn_dy
;
;plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(0)+0.013,linestyle=1,thick=3.0,/NORMAL
;plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013,linestyle=3,thick=3.0,/NORMAL
;xyouts,txt_lgn_x,lgn_y(0),'Existing Channel',size=txt_lgn_sz,/NORMAL
;xyouts,txt_lgn_x,lgn_y(1),'Candidate Channel',size=txt_lgn_sz,/NORMAL
;
xyouts,1.0,0.95,info_1_sng,alignment=0.0,orientation=270.0,size=0.65,/NORMAL
xyouts,0.985,0.95,info_2_sng,alignment=0.0,orientation=270.0,size=0.65,/NORMAL
xyouts,0.97,0.95,info_3_sng,alignment=0.0,orientation=270.0,size=0.65,/NORMAL
xyouts,0.955,0.95,info_4_sng,alignment=0.0,orientation=270.0,size=0.65,/NORMAL
;
if ann then ntt
;
print,'Hit any key to continue, or q to quit ...'
junk=get_kbrd(1)
if junk eq 'q' then goto,end_of_procedure
;
endif; endif RTA
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of reflectance graph
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Start multi-transmittance graphs
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if TRN then begin
;
nbr_trn=8
trn_spc_atm=fltarr(nbr_bnd,nbr_trn)
trn_spc_atm(*,0)=trn_spc_atm_CO2
trn_spc_atm(*,1)=trn_spc_atm_H2O
trn_spc_atm(*,2)=trn_spc_atm_O2
trn_spc_atm(*,3)=trn_spc_atm_O3
trn_spc_atm(*,4)=trn_spc_atm_Ray
trn_spc_atm(*,5)=trn_spc_atm_ice
trn_spc_atm(*,6)=trn_spc_atm_liq
trn_spc_atm(*,7)=trn_spc_atm_ttl
;
title=strarr(nbr_trn)
title(0)='!5CO!I2!N Transmittance'
title(1)='!5H!I2!NO Transmittance'
title(2)='!5O!I2!N Transmittance'
title(3)='!5O!I3!N Transmittance'
title(4)='!5Rayleigh Scattering Transmittance'
title(5)='!5Ice Water Transmittance'
title(6)='!5Liquid Water Transmittance'
title(7)='!5Total Transmittance'
;
for trn_gph=0,nbr_trn-1 do begin
;
ord_max=1.
;
plot, $
	wvl_ctr, $
	trn_spc_atm(*,trn_gph), $
	tit=title(trn_gph), $
	xtit='!5Wavelength !7k!5 (!7l!5m)', $
	ytit='!5Transmittance !8T!I!7k!5!N', $
	xrange=[.2,2.5], $
;	yrange=[.2,2.5], $
	xmargin=[5.5,3.], $	;[10,3] is [left,right] default
	ymargin=[3.0,1.0], $	;[4,2] is [bottom,top] default
	xstyle=0, $
	ystyle=0, $
	thick=1.0, $
	charsize=2.0, $
	linestyle=0
;
; Loop over the current and projected channels
;
for chan=0,6 do begin
oplot,	$
	[chans(chan).TDDR_wvl_ctr,chans(chan).TDDR_wvl_ctr], $
	[0.0,ord_max], $
	thick=1.0, $
	linestyle=1
endfor ; end loop over current channels
;
for chan=7,10 do begin
oplot,	$
	[chans(chan).TDDR_wvl_ctr,chans(chan).TDDR_wvl_ctr], $
	[0.0,ord_max], $
	thick=1.0, $
	linestyle=3
endfor ; end loop over new channels
;
ln_lgn_x1=0.6
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=0.8
lgn_dy=0.1
lgn_y=lgn_y_top-indgen(6)*lgn_dy
;
;plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(0)+0.013,linestyle=1,thick=3.0,/NORMAL
;plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013,linestyle=3,thick=3.0,/NORMAL
;xyouts,txt_lgn_x,lgn_y(0),'Existing Channel',size=txt_lgn_sz,/NORMAL
;xyouts,txt_lgn_x,lgn_y(1),'Candidate Channel',size=txt_lgn_sz,/NORMAL
;
xyouts,1.0,0.95,info_1_sng,alignment=0.0,orientation=270.0,size=0.65,/NORMAL
xyouts,0.985,0.95,info_2_sng,alignment=0.0,orientation=270.0,size=0.65,/NORMAL
xyouts,0.97,0.95,info_3_sng,alignment=0.0,orientation=270.0,size=0.65,/NORMAL
xyouts,0.955,0.95,info_4_sng,alignment=0.0,orientation=270.0,size=0.65,/NORMAL
;
if ann then ntt
;
print,'Hit any key to continue, or q to quit ...'
junk=get_kbrd(1)
if junk eq 'q' then goto,end_of_procedure
;
endfor; end loop over transmittance graphs
;
endif; endif TRN
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of multi-transmittance graphs
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
sim=sim+1
if srt_idx eq end_idx then goto,end_of_procedure
if cmd_ln_fl then goto,end_of_procedure
endfor ; end loop over sims
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Start radiance ratio graph
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Replace the '_top' string with '_btm' and switch the commented with
; the uncommented labels to change which level for the radiances.
;
ord_max=max(chn_rad_btm(1:nbr_chn-1,*))
chan=1
;
plot, $
	chn_rad_btm(0,*), $ ; visible channel
	chn_rad_btm(chan,*), $
;	tit='!5Radiance Ratios at !8z!5 = 13 km', $
;	xtit='!5Visible Radiance !8I!S!E!97!5!N!R!I!7k!5!N(.5 !7l!5m) (!5W m!E-2!N !7l!5m!E-1!N sr!E-1!N)', $
;	ytit='!5Channel Radiance !8I!S!E!97!5!N!R!I!7k!5!N (!5W m!E-2!N !7l!5m!E-1!N sr!E-1!N)', $
	tit='!5Radiance Ratios at !8z!5 = 1 km', $
	xtit='!5Visible Radiance !8I!S!E!95!5!N!R!I!7k!5!N(.5 !7l!5m) (!5W m!E-2!N !7l!5m!E-1!N sr!E-1!N)', $
	ytit='!5Channel Radiance !8I!S!E!95!5!N!R!I!7k!5!N (!5W m!E-2!N !7l!5m!E-1!N sr!E-1!N)', $
;	xrange=[.2,2.5], $
	yrange=[0.0,ord_max], $
	xmargin=[5.5,3.], $	;[10,3] is [left,right] default
	ymargin=[3.0,1.0], $	;[4,2] is [bottom,top] default
	xstyle=0, $
	ystyle=0, $
	thick=1.0, $
	charsize=2.0, $
	psym=-chan
;
for chan=2,6 do begin
oplot, $
	chn_rad_btm(0,*), $ ; visible channel
	chn_rad_btm(chan,*), $
	thick=1.0, $
	psym=-((chan mod 6) + 1)
endfor ; end loop over current channels
;
if ann then ntt
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End radiance ratio graph
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
end_of_procedure: foo=1
;
end; end aca_orig()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure ACA Original commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure ACA NREL
; Scatterplot the TDDR channel intensity vs. the .5 microns intensity.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro aca_nrel, $
	srt_idx=srt_idx, $
	end_idx=end_idx, $
	ncr_idx=ncr_idx, $
	chn_bnd=chn_bnd, $
	fld_nbr=fld_nbr, $
	top_int=top_int, $
	btm_int=btm_int, $
	fl_in=fl_in, $
	fl_cmp=fl_cmp, $
	fzz=fzz, $
	NTN=NTN, $
	FLX=FLX, $
	TRN=TRN, $
	CLR=CLR, $
	RTA=RTA, $
	TWO=TWO, $
	SWCF=SWCF, $
	ann=ann, $
	level=level
;
; aca_nrel,srt_idx=0,end_idx=0,ncr_idx=25
; aca_nrel,srt_idx=100,end_idx=100,ncr_idx=0
; aca_nrel,RTA=0,FLX=1,NTN=1,fzz=1,TRN=0,fl_in='/data/zender/aca/aca.100.nc',CLR='/data/zender/aca/aca.000.nc'
; aca_nrel,RTA=0,FLX=1,NTN=1,fzz=0,TRN=0,TWO=1,srt_idx=12,end_idx=18,ncr_idx=1
; aca_nrel,RTA=0,FLX=0,NTN=1,fzz=0,TRN=0,fl_in='/data/zender/aca/aca.ice.nc',fl_cmp='/data/zender/aca/aca.100.nc'
; aca_nrel,fl_in='/data/zender/aca/nrel_clr.nc',CLR='/data/zender/aca/nrel_clr.nc',btm_int=18,top_int=18,fzz=1
; aca_nrel,fl_in='/data/zender/aca/nrel_cld.nc',CLR='/data/zender/aca/nrel_cld.nc',btm_int=18,top_int=18,fzz=1
; aca_nrel,fl_in='/data/zender/aca/950417_1940_nrel_d2n_cld.nc',CLR='/data/zender/aca/950417_1940_nrel_d2n_clr.nc',btm_int=30,top_int=0,fzz=1
; aca_nrel,fl_in='/data/zender/aca/931021_2200_nrel_d2n_clr.nc',CLR='/data/zender/aca/931021_2200_nrel_d2n_clr.nc',btm_int=28,top_int=0,fzz=1
; aca_nrel,fl_in='/data/zender/aca/950513_2030_nrel_d2n_clr.nc',CLR='/data/zender/aca/950513_2030_nrel_d2n_clr.nc',btm_int=28,top_int=0,fzz=1
;
True=1
False=0
pi=3.141592654
;
if n_elements(srt_idx) eq 0 then srt_idx=0
if n_elements(end_idx) eq 0 then end_idx=0
if n_elements(ncr_idx) eq 0 then ncr_idx=25
if n_elements(chn_bnd) eq 0 then chn_bnd=800
if n_elements(fld_nbr) eq 0 then fld_nbr=2
if n_elements(top_int) eq 0 then top_int=7
if n_elements(btm_int) eq 0 then btm_int=15
if n_elements(fl_in) ne 0 then cmd_ln_fl=True else cmd_ln_fl=False
if n_elements(fl_cmp) ne 0 then CMP_FLG=True else CMP_FLG=False
if n_elements(fzz) eq 0 then fzz=False
if n_elements(NTN) eq 0 then NTN=False
if n_elements(FLX) eq 0 then FLX=True
if n_elements(TRN) eq 0 then TRN=False
if n_elements(CLR) eq 0 then CLR='/data/zender/aca/aca.000.nc'
if n_elements(RTA) eq 0 then RTA=False
if n_elements(TWO) eq 0 then TWO=False
if n_elements(SWCF) eq 0 then SWCF=True
if n_elements(ann) eq 0 then ann=False
if n_elements(vis_bnd) eq 0 then vis_bnd=1603
;if n_elements(foo) eq 0 then foo=1
;
if CMP_FLG then begin
print,'Gathering CMP_FLG data from '+fl_cmp
nc_id=ncdf_open(fl_cmp)
levp_id=ncdf_dimid(nc_id,'levp')
ncdf_diminq,nc_id,levp_id,foo,levp_nbr
ncdf_varget,nc_id,'ntn_spc_aa_ndr',ntn_spc_aa_ndr
ncdf_varget,nc_id,'ntn_spc_aa_zen',ntn_spc_aa_zen
;ncdf_varget,nc_id,'foo',foo_clr
ncdf_close,nc_id
ntn_spc_aa_ndr2=ntn_spc_aa_ndr(*,btm_int)
ntn_spc_aa_zen2=ntn_spc_aa_zen(*,top_int)
endif; endif CMP_FLG
;
;NREL='/data/zender/aca/nrel/nrel_obs_01.nc'
NREL='/data/zender/aca/nrel/nrel_obs_02.nc'
if NREL then begin
print,'Gathering NREL data from '+NREL
nc_id=ncdf_open(NREL)
nrel_bnd_id=ncdf_dimid(nc_id,'bnd')
ncdf_diminq,nc_id,nrel_bnd_id,foo,nrel_nbr_bnd
ncdf_varget,nc_id,'flx_spc_dwn_sfc_cld_1',nrel_flx_spc_dwn_sfc_cld_1
ncdf_varget,nc_id,'flx_spc_dwn_sfc_cld_2',nrel_flx_spc_dwn_sfc_cld_2
;ncdf_varget,nc_id,'flx_spc_dwn_sfc_cld_3',nrel_flx_spc_dwn_sfc_cld_3
ncdf_varget,nc_id,'flx_spc_dwn_sfc_clr',nrel_flx_spc_dwn_sfc_clr
ncdf_varget,nc_id,'wvl_ctr',nrel_wvl_ctr
;ncdf_varget,nc_id,'foo',nrel_foo
ncdf_close,nc_id
endif; endif NREL
;
; Construct the file name
;
if cmd_ln_fl then fl_nm=fl_in
print,'Gathering data from '+fl_nm
nc_id=ncdf_open(fl_nm)
;
; Get the dimension sizes
;
bnd_id=ncdf_dimid(nc_id,'bnd')
ncdf_diminq,nc_id,bnd_id,foo,nbr_bnd
levp_id=ncdf_dimid(nc_id,'levp')
ncdf_diminq,nc_id,levp_id,foo,levp_nbr
;
; Get the two-dimensional arrays
;
ncdf_varget,nc_id,'flx_spc_dwn',flx_spc_dwn
ncdf_varget,nc_id,'flx_spc_dwn_TOA',flx_spc_dwn_TOA
ncdf_varget,nc_id,'flx_spc_dwn_sfc',flx_spc_dwn_sfc
ncdf_varget,nc_id,'flx_bb_dwn_TOA',flx_bb_dwn_TOA
ncdf_varget,nc_id,'flx_bb_dwn_sfc',flx_bb_dwn_sfc
ncdf_varget,nc_id,'flx_bb_dwn_drc',flx_bb_dwn_drc
ncdf_varget,nc_id,'flx_bb_dwn_dff',flx_bb_dwn_dff
ncdf_varget,nc_id,'flx_bb_dwn',flx_bb_dwn
ncdf_varget,nc_id,'flx_bb_up',flx_bb_up
ncdf_varget,nc_id,'wvl_ctr',wvl_ctr
ncdf_varget,nc_id,'wvl_dlt',wvl_dlt
ncdf_varget,nc_id,'trn_spc_atm_ttl',trn_spc_atm_ttl
ncdf_varget,nc_id,'rfl_spc_SAS',rfl_spc_SAS
ncdf_varget,nc_id,'abs_spc_atm',abs_spc_atm
ncdf_varget,nc_id,'CWP_ttl',CWP_ttl
ncdf_varget,nc_id,'frc_ice_ttl',frc_ice_ttl
ncdf_varget,nc_id,'alb_sfc',alb_sfc
ncdf_attget,nc_id,'prf_sng',prf_sng,/GLOBAL
ncdf_varget,nc_id,'lcl_time_hr',lcl_time_hr
ncdf_varget,nc_id,'lcl_yr_day',lcl_yr_day
ncdf_varget,nc_id,'z_cld_btm',z_cld_btm
ncdf_varget,nc_id,'lat_dgr',lat_dgr
ncdf_varget,nc_id,'z_cld_thick',z_cld_thick
ncdf_varget,nc_id,'slr_zen_ngl_cos',slr_zen_ngl_cos
ncdf_varget,nc_id,'flx_bb_abs_atm',flx_bb_abs_atm
ncdf_varget,nc_id,'flx_bb_abs_sfc',flx_bb_abs_sfc
ncdf_varget,nc_id,'rfl_bb_SAS',rfl_bb_SAS
ncdf_varget,nc_id,'flx_bb_net',flx_bb_net
ncdf_varget,nc_id,'opt_dep_ext_clm_ttl',opt_dep_ext_clm_ttl
ncdf_varget,nc_id,'alt_ntf',alt_ntf
;
; End of NetCDF commands
;
ncdf_close,nc_id
;
; The clear sky file might change with each file
; like when we're varying the solar zenith angle, so this needs
; to go in the main time loop.
;
if SWCF then begin
print,'Gathering clear sky data from '+CLR
nc_id=ncdf_open(CLR)
levp_id=ncdf_dimid(nc_id,'levp')
ncdf_diminq,nc_id,levp_id,foo,levp_nbr
ncdf_varget,nc_id,'flx_bb_net',flx_bb_net_clr
;ncdf_varget,nc_id,'foo',foo_clr
ncdf_close,nc_id
flx_bb_net_clr_TOA=flx_bb_net_clr(0)
flx_bb_net_clr_sfc=flx_bb_net_clr(levp_nbr-1)
flx_bb_net_clr_top=flx_bb_net_clr(top_int)
flx_bb_net_clr_btm=flx_bb_net_clr(btm_int)
endif; endif SWCF
;
flx_spc_dwn_top=flx_spc_dwn(*,top_int)
flx_spc_dwn_btm=flx_spc_dwn(*,btm_int)
flx_bb_net_TOA=flx_bb_net(0)
flx_bb_net_sfc=flx_bb_net(levp_nbr-1)
flx_bb_net_top=flx_bb_net(top_int)
flx_bb_net_btm=flx_bb_net(btm_int)
flx_bb_dwn_top=flx_bb_dwn(top_int)
flx_bb_dwn_btm=flx_bb_dwn(btm_int)
flx_bb_dwn_drc_top=flx_bb_dwn_drc(top_int)
flx_bb_dwn_drc_btm=flx_bb_dwn_drc(btm_int)
flx_bb_dwn_dff_top=flx_bb_dwn_dff(top_int)
flx_bb_dwn_dff_btm=flx_bb_dwn_dff(btm_int)
flx_bb_up_top=flx_bb_up(top_int)
flx_bb_up_btm=flx_bb_up(btm_int)
flx_bb_abs_cld=(flx_bb_net_top-flx_bb_net_btm)/flx_bb_dwn_TOA
opt_dep_ext_cld_vis=0.01*(round(100.*(opt_dep_ext_clm_ttl(1602)-.1472))) ; subtract clear sky value
;
if FZZ then begin
if FLX then begin
	print,'Averaging FLXes...'
	flx_spc_dwn_TOA=rnn_avg(flx_spc_dwn_TOA,wvl_ctr,crd_sz=wvl_dlt,1.0e-8)
	flx_spc_dwn_sfc=rnn_avg(flx_spc_dwn_sfc,wvl_ctr,crd_sz=wvl_dlt,1.0e-8)
	flx_spc_dwn_top=rnn_avg(flx_spc_dwn_top,wvl_ctr,crd_sz=wvl_dlt,1.0e-8)
	flx_spc_dwn_btm=rnn_avg(flx_spc_dwn_btm,wvl_ctr,crd_sz=wvl_dlt,1.0e-8)
endif; endif FLX averaging
endif ; end if fzz
;
flx_spc_dwn_TOA=flx_spc_dwn_TOA*1.0e-9 ; W/m2/m --> W/m2/nm
flx_spc_dwn_sfc=flx_spc_dwn_sfc*1.0e-9 ; W/m2/m --> W/m2/nm
flx_spc_dwn_top=flx_spc_dwn_top*1.0e-9 ; W/m2/m --> W/m2/nm
flx_spc_dwn_btm=flx_spc_dwn_btm*1.0e-9 ; W/m2/m --> W/m2/nm
;foo=foo*1.0e-9 ; W/m2/m --> W/m2/nm
if NREL then begin
bad_idx=where(nrel_flx_spc_dwn_sfc_cld_1 le 0.0,cnt)
if cnt ne 0 then nrel_flx_spc_dwn_sfc_cld_1(bad_idx)=1.0e-1
bad_idx=where(nrel_flx_spc_dwn_sfc_cld_2 le 0.0,cnt)
if cnt ne 0 then nrel_flx_spc_dwn_sfc_cld_2(bad_idx)=1.0e-1
bad_idx=where(nrel_flx_spc_dwn_sfc_clr le 0.0,cnt)
if cnt ne 0 then nrel_flx_spc_dwn_sfc_clr(bad_idx)=1.0e-1
nrel_flx_spc_dwn_sfc_cld_1=nrel_flx_spc_dwn_sfc_cld_1*1.0e-9 ; W/m2/m --> W/m2/nm
nrel_flx_spc_dwn_sfc_cld_2=nrel_flx_spc_dwn_sfc_cld_2*1.0e-9 ; W/m2/m --> W/m2/nm
nrel_flx_spc_dwn_sfc_clr=nrel_flx_spc_dwn_sfc_clr*1.0e-9 ; W/m2/m --> W/m2/nm
nrel_wvl_ctr=nrel_wvl_ctr*1.0e6 ; m --> micron
endif; endif NREL
;
wvl_ctr=wvl_ctr*1.0e6 ; m --> micron
wvl_dlt=wvl_dlt*1.0e6 ; m --> micron
;
if CWP_ttl eq 0. then begin
opt_dep_ext_cld_vis=0.
endif; endif there is no cloud
if lat_dgr ge 0. then hemi_sng='N' else hemi_sng='S'
prf_sng=strtrim(string(prf_sng),2)
btm_sng=auto_sng(fix(alt_ntf(btm_int)/1000.),0)+' km'
top_sng=auto_sng(fix(alt_ntf(top_int)/1000.),0)+' km'
info_1_sng=prf_sng+', Lat = '+auto_sng(lat_dgr,2)+'!E!12_!5!N'+hemi_sng+', Day = '+auto_sng(fix(lcl_yr_day),0)+', Hour = '+auto_sng(lcl_time_hr,2)+', !7H!5 = '+auto_sng(180.*acos(slr_zen_ngl_cos)/pi,2)+'!E!12_!5!N, cos(!7H!5) = '+auto_sng(slr_zen_ngl_cos,2)+', !8A!5!Isfc!N = '+auto_sng(alb_sfc,2)
info_2_sng='!8CWP!5 = '+auto_sng(round(CWP_ttl*1000.),0)+' g m!E-2!N, Base = '+auto_sng(z_cld_btm/1000.0,1)+' km, Thick = '+auto_sng(z_cld_thick/1000.0,1)+' km, '+' !8f!5!Iice!N = '+auto_sng(round(100.*frc_ice_ttl),0)+' %'
info_3_sng='!8A!5!ITOA!N =  '+auto_sng(rfl_bb_SAS,2)+', !7s!5!Icld!N(.5 !7l!5m) = '+auto_sng(opt_dep_ext_cld_vis,2)+', !7j!5(cld) = '+auto_sng(flx_bb_abs_cld,2)
info_3_sng=info_3_sng+', !8F!5!95!5!N(TOA) = '+auto_sng(round(flx_bb_dwn_TOA),0)+', !8F!5!95!5!N(sfc) = '+auto_sng(round(flx_bb_dwn_sfc),0)+', !8F!5!Iabs!N(atm) = '+auto_sng(round(flx_bb_abs_atm),0)+', !8F!5!Iabs!N(sfc) = '+auto_sng(round(flx_bb_abs_sfc),0)+', !8F!5!Iabs!N(cld) = '+auto_sng(round(flx_bb_net_top-flx_bb_net_btm),0)
;
info_4_sng=top_sng+': !8F!5!95!5!N = '+auto_sng(round(flx_bb_dwn_top),0)+', !8F!5!S!E!95!5!N!R!Idir!N = '+auto_sng(round(flx_bb_dwn_drc_top),0)+', !8F!5!97!5!N = '+auto_sng(round(flx_bb_up_top),0)+', !8F!5!INet!N = '+auto_sng(round(flx_bb_net_top),0)
info_4_sng=info_4_sng+'. '+btm_sng+': !8F!5!95!5!N = '+auto_sng(round(flx_bb_dwn_btm),0)+', !8F!5!S!E!95!5!N!R!Idir!N = '+auto_sng(round(.01*(round(100.*flx_bb_dwn_drc_btm))),0)+', !8F!5!97!5!N = '+auto_sng(round(flx_bb_up_btm),0)+', !8F!5!INet!N = '+auto_sng(round(flx_bb_net_btm),0)
if SWCF then begin
SWCF_TOA=flx_bb_net_clr_TOA-flx_bb_net_TOA
SWCF_sfc=flx_bb_net_clr_sfc-flx_bb_net_sfc
SWCF_top=flx_bb_net_clr_top-flx_bb_net_top
SWCF_btm=flx_bb_net_clr_btm-flx_bb_net_btm
if SWCF_TOA ne 0. then SWCF_rat_clm=SWCF_sfc/SWCF_TOA else SWCF_rat_clm=1.
if SWCF_top ne 0. then SWCF_rat_cld=SWCF_btm/SWCF_top else SWCF_rat_cld=1.
endif else begin; endif SWCF
SWCF_TOA=0.
SWCF_sfc=0.
SWCF_top=0.
SWCF_btm=0.
SWCF_rat_clm=1.
SWCF_rat_cld=1.
endelse; endif no clear sky file
;if not CMP_FLG then begin
info_2_sng=info_2_sng+'. CF(TOA) = '+auto_sng(round(SWCF_TOA),0)+', CF(sfc) = '+auto_sng(round(SWCF_sfc),0)+', CFR(col) = '+auto_sng(SWCF_rat_clm,2)+', CFR(cld) = '+auto_sng(SWCF_rat_cld,2)
;endif; endif not CMP_FLG
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Start model flux graph
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if FLX then begin
;
ord_max=max(flx_spc_dwn_TOA)
;
plot_io, $
	wvl_ctr, $
	flx_spc_dwn_TOA, $
	tit='!5Flux at TOA and Surface', $
	xtit='!5Wavelength !7k!5 (!7l!5m)', $
	ytit='!5Flux !8F!I!7k!5!N (!5W m!E-2!N nm!E-1!N)', $
	xrange=[.2,2.5], $
	yrange=[1.0e-6,1.0e1], $
	xmargin=[5.75,3.], $	;[10,3] is [left,right] default
	ymargin=[3.0,1.0], $	;[4,2] is [bottom,top] default
	xstyle=1, $
	ystyle=1, $
	thick=1.0, $
	yticklen=1.0, $
	charsize=2.0, $
	linestyle=0
;
oplot,	$
	wvl_ctr, $
	flx_spc_dwn_btm, $
	thick=1.0, $
	linestyle=0
;
ln_lgn_x1=0.6
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=0.9
lgn_dy=0.1
lgn_y=lgn_y_top-indgen(6)*lgn_dy
;
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(0)+0.013,linestyle=0,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013,linestyle=0,thick=3.0,/NORMAL
;
xyouts,txt_lgn_x,lgn_y(0),'!8F!S!I!7k!R!E!95!5!N TOA',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(1),'!8F!S!I!7k!R!E!95!5!N '+btm_sng,size=txt_lgn_sz,/NORMAL
;
xyouts,1.0,0.95,info_1_sng,alignment=0.0,orientation=270.0,size=0.65,/NORMAL
xyouts,0.985,0.95,info_2_sng,alignment=0.0,orientation=270.0,size=0.65,/NORMAL
xyouts,0.97,0.95,info_3_sng,alignment=0.0,orientation=270.0,size=0.65,/NORMAL
xyouts,0.955,0.95,info_4_sng,alignment=0.0,orientation=270.0,size=0.65,/NORMAL
;
if ann then ntt
;
print,'Hit any key to continue, or q to quit ...'
junk=get_kbrd(1)
if junk eq 'q' then goto,end_of_procedure
;
endif; endif FLX
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End model flux graph
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Start clear sky graph
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if FLX then begin
;
max_flx_spc_dwn_TOA=max(flx_spc_dwn_TOA)
max_nrel_flx_spc_dwn_sfc_clr=max(nrel_flx_spc_dwn_sfc_clr)
scale=max_nrel_flx_spc_dwn_sfc_clr/max_flx_spc_dwn_TOA
;flx_spc_dwn_TOA=flx_spc_dwn_TOA*scale
;flx_spc_dwn_btm=flx_spc_dwn_btm*scale
ord_max=max([max_flx_spc_dwn_TOA,max_nrel_flx_spc_dwn_sfc_clr])
;
;plot_io, $
plot, $
	wvl_ctr, $
	flx_spc_dwn_TOA, $
	tit='!5Clear Sky Flux at TOA and Surface', $
	xtit='!5Wavelength !7k!5 (!7l!5m)', $
	ytit='!5Flux !8F!I!7k!5!N (!5W m!E-2!N nm!E-1!N)', $
	xrange=[.2,2.5], $
;	yrange=[1.0e-4,1.0e1], $
	yrange=[0.0,ord_max], $
	xmargin=[5.75,3.], $	;[10,3] is [left,right] default
	ymargin=[3.0,1.0], $	;[4,2] is [bottom,top] default
	xstyle=1, $
	ystyle=1, $
	thick=1.0, $
	yticklen=1.0, $
	charsize=2.0, $
	linestyle=0
;
oplot,	$
	nrel_wvl_ctr, $
	nrel_flx_spc_dwn_sfc_clr, $
	thick=1.0, $
	linestyle=0
;
oplot,	$
	wvl_ctr, $
	flx_spc_dwn_btm, $
	thick=1.0, $
	linestyle=3
;
ln_lgn_x1=0.6
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=0.9
lgn_dy=0.1
lgn_y=lgn_y_top-indgen(6)*lgn_dy
;
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(0)+0.013,linestyle=0,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013,linestyle=0,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(2)+0.013,linestyle=3,thick=3.0,/NORMAL
;
xyouts,txt_lgn_x,lgn_y(0),'TOA model',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(1),'Surface observed',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(2),'Surface model',size=txt_lgn_sz,/NORMAL
;
xyouts,1.0,0.95,info_1_sng,alignment=0.0,orientation=270.0,size=0.65,/NORMAL
xyouts,0.985,0.95,info_2_sng,alignment=0.0,orientation=270.0,size=0.65,/NORMAL
xyouts,0.97,0.95,info_3_sng,alignment=0.0,orientation=270.0,size=0.65,/NORMAL
xyouts,0.955,0.95,info_4_sng,alignment=0.0,orientation=270.0,size=0.65,/NORMAL
;
if ann then ntt
;
print,'Hit any key to continue, or q to quit ...'
junk=get_kbrd(1)
if junk eq 'q' then goto,end_of_procedure
;
endif; endif FLX
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End clear sky graph
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Start cloudy sky graph
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if FLX then begin
;
max_flx_spc_dwn_TOA=max(flx_spc_dwn_TOA)
max_nrel_flx_spc_dwn_sfc_cld_1=max(nrel_flx_spc_dwn_sfc_cld_1)
scale=max_nrel_flx_spc_dwn_sfc_cld_1/max_flx_spc_dwn_TOA
;flx_spc_dwn_TOA=flx_spc_dwn_TOA*scale
;flx_spc_dwn_btm=flx_spc_dwn_btm*scale
ord_max=max([max_flx_spc_dwn_TOA,max_nrel_flx_spc_dwn_sfc_cld_1])
;
;plot_io, $
plot, $
	wvl_ctr, $
	flx_spc_dwn_TOA, $
	tit='!5Cloudy Sky Flux at TOA and Surface', $
	xtit='!5Wavelength !7k!5 (!7l!5m)', $
	ytit='!5Flux !8F!I!7k!5!N (!5W m!E-2!N nm!E-1!N)', $
	xrange=[.2,2.5], $
;	yrange=[1.0e-6,1.0e1], $
	yrange=[0.0,ord_max], $
	xmargin=[5.75,3.], $	;[10,3] is [left,right] default
	ymargin=[3.0,1.0], $	;[4,2] is [bottom,top] default
	xstyle=1, $
	ystyle=1, $
	thick=1.0, $
	yticklen=1.0, $
	charsize=2.0, $
	linestyle=0
;
oplot,	$
	nrel_wvl_ctr, $
	nrel_flx_spc_dwn_sfc_cld_1, $
	thick=1.0, $
	linestyle=0
;
oplot,	$
	nrel_wvl_ctr, $
	nrel_flx_spc_dwn_sfc_cld_2, $
	thick=1.0, $
	linestyle=0
;
oplot,	$
	wvl_ctr, $
	flx_spc_dwn_btm, $
	thick=1.0, $
	linestyle=3
;
ln_lgn_x1=0.6
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=0.9
lgn_dy=0.1
lgn_y=lgn_y_top-indgen(6)*lgn_dy
;
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(0)+0.013,linestyle=0,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013,linestyle=0,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(2)+0.013,linestyle=0,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(3)+0.013,linestyle=3,thick=3.0,/NORMAL
;
xyouts,txt_lgn_x,lgn_y(0),'TOA model',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(1),'Surface observed',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(2),'Surface observed',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(3),'Surface model',size=txt_lgn_sz,/NORMAL
;
xyouts,1.0,0.95,info_1_sng,alignment=0.0,orientation=270.0,size=0.65,/NORMAL
xyouts,0.985,0.95,info_2_sng,alignment=0.0,orientation=270.0,size=0.65,/NORMAL
xyouts,0.97,0.95,info_3_sng,alignment=0.0,orientation=270.0,size=0.65,/NORMAL
xyouts,0.955,0.95,info_4_sng,alignment=0.0,orientation=270.0,size=0.65,/NORMAL
;
if ann then ntt
;
print,'Hit any key to continue, or q to quit ...'
junk=get_kbrd(1)
if junk eq 'q' then goto,end_of_procedure
;
endif; endif FLX
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End cloudy sky graph
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
end_of_procedure: foo=1
;
end; end aca_nrel()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure ACA NREL commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

