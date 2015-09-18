; $Id$

; Purpose: Holds functions which should always be loaded before procedures, 
; and procedures whose use is generic and should be in a utilities file.
; By their nature, these procedures should not require access to any common blocks.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin ncdf_varget_2D
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro ncdf_varget_2D,nc_id,var_nm,var_val,fl_idx
; Purpose: Place new array into input array and return changed array 
; fl_idx is the location in the final dimension of var_val where the new array will be placed
; This routine is useful when building 2D arrays by concatentating 1D arrays from a list of files
; For example, suppose var_val=fltarr(crd_nbr,fl_nbr) 
; Then in a loop over files, we fill in var_val by calling
; ncdf_varget_2D,nc_id,var_nm,var_val,fl_idx ; (fl_idx=fl_idx+1)
; Calling syntax is same as vanilla ncdf_varget(), but with additional index argument (fl_idx) at end
; Routine allows for extra space in var_val (but only when var_nm on disk is one dimensional)
; See also ncdf_varget_1D, which inserts netCDF data into one dimensional arrays
dmn_nbr_out=(size(var_val))(0)
data_nbr_out=(size(var_val))(dmn_nbr_out-1)
if dmn_nbr_out le 0 or dmn_nbr_out gt 3 then begin
	print,'ERROR: ncdf_varget_2D() reports var_nm = ',var_nm,', var_val has ',auto_sng(dmn_nbr_out,0),' dimensions, fl_idx = ',auto_sng(fl_idx,0)
endif; endif err
ncdf_varget,nc_id,var_nm,data_foo
dmn_nbr_in=(size(data_foo))(0)
data_nbr_in=(size(data_foo))(dmn_nbr_in+2)
if data_nbr_in gt data_nbr_out then begin
	print,'ERROR: ncdf_varget_2D() reports var_nm = ',var_nm,', dmn_nbr_in = ',auto_sng(dmn_nbr_in,0),', data_nbr_in = ',auto_sng(data_nbr_in,0),' dmn_nbr_out = ',auto_sng(dmn_nbr_out,0),' data_nbr_out = ',auto_sng(data_nbr_out,0)
endif; endif err
if dmn_nbr_out eq 1 then var_val(fl_idx)=data_foo else if dmn_nbr_out eq 2 then var_val(0:data_nbr_in-1,fl_idx)=data_foo else if dmn_nbr_out eq 3 then var_val(*,*,fl_idx)=data_foo
end; ncdf_varget_2D()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End ncdf_varget_2D()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin ncdf_varget_1D
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro ncdf_varget_1D,nc_id,var_nm,var_val,crd_idx,dmn0_idx=dmn0_idx,dbg_lvl=dbg_lvl,val_dbg=val_dbg
; Purpose: Insert new array into specified location in one dimensional input array and return changed array 
; crd_idx=[idx_srt,idx_end] are bounding indices in one dimension input array where var_nm array will be placed
; This routine is useful when building 1D arrays by concatentating 1D arrays from a list of files
; For example, suppose var_val=fltarr(sz_nbr_ttl) and idx_srt=idx_end=intarr(fl_nbr)
; Then in a loop over files, we fill in var_val by calling
; ncdf_varget_1D,nc_id,var_nm,var_val,[idx_srt(fl_idx),idx_end(fl_idx)] ; (fl_idx=fl_idx+1)
; Calling syntax is same as vanilla ncdf_varget(), but with additional index array argument (crd_idx) at end
; Optional argument dmn0_idx specifies, for two dimensional disk variables only, offset of desired array
; See also ncdf_var_get_2D, which inserts netCDF data into multidimensional arrays
if n_elements(dmn0_idx) eq 0 then dmn0_idx=0 ; [idx] Dimension index
if n_elements(dbg_lvl) eq 0 then dbg_lvl=0 ; [idx] Debugging level
if n_elements(val_dbg) eq 0 then val_dbg=0.0 ; [???] Debugging value fxm: not fully implemented yet

hyp_flg=0 ; [flg] Hyperslab input variable
dmn_nbr_out=(size(var_val))(0)
data_nbr_out=(size(var_val))(dmn_nbr_out)
if dmn_nbr_out ne 1 then begin
	print,'WARNING: ncdf_varget_1D() reports var_nm = ',var_nm,', var_val has ',auto_sng(dmn_nbr_out,0),' dimensions, fl_idx = ',auto_sng(fl_idx,0)
endif; endif err
ncdf_varget,nc_id,var_nm,data_foo
dmn_nbr_in=(size(data_foo))(0)
data_nbr_in=(size(data_foo))(dmn_nbr_in+2)
if data_nbr_in gt data_nbr_out and dmn_nbr_in eq 2 then begin
	data_nbr_in=(size(data_foo))(1)
	if dbg_lvl then print,'WARNING: ncdf_varget_1D() 2D -> 1D reading IDL ',var_nm,'(0:',auto_sng(data_nbr_in-1,0),',',auto_sng(dmn0_idx,0),') = C ',var_nm,'[',auto_sng(dmn0_idx,0),'][0:',auto_sng(data_nbr_in-1,0),']'
	hyp_flg=1 ; [flg] Hyperslab input variable
endif; endif dmn_nbr
if data_nbr_in gt data_nbr_out then begin
	print,'ERROR: ncdf_varget_1D() reports var_nm = ',var_nm,', dmn_nbr_in = ',auto_sng(dmn_nbr_in,0),', data_nbr_in = ',auto_sng(data_nbr_in,0),' dmn_nbr_out = ',auto_sng(dmn_nbr_out,0),' data_nbr_out = ',auto_sng(data_nbr_out,0)
endif; endif err
crd_fll_sz_out=crd_idx(1)-crd_idx(0)+1
if crd_fll_sz_out ne data_nbr_in then begin
	print,'WARNING: ncdf_varget_1D() reports ',var_nm,' has crd_fll_sz_out = ',auto_sng(crd_fll_sz_out,0),' but data_nbr_in = ',auto_sng(data_nbr_in)
endif; endif err
if hyp_flg then var_val(crd_idx(0):crd_idx(1))=data_foo(*,dmn0_idx) else var_val(crd_idx(0):crd_idx(1))=data_foo
end; ncdf_varget_1D()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End ncdf_varget_1D()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin mth2sng
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function mth2sng,mth_idx
; Purpose: Convert month index to string
; Usage: mth_sng=mth2sng(3) returns 'Mar'
mth_mdm=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
mth_sng=mth_mdm(mth_idx-1)
return,mth_sng
end; mth2sng()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End mth2sng()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin szgrd_mk
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro szgrd_mk,sz_mnm,sz_mxm,sz_nbr,grd_sng=grd_sng,sz_grd=sz_grd,sz_min=sz_min,sz_max=sz_max,sz_dlt=sz_dlt,sz_ctr=sz_ctr
; Purpose: Construct a size grid
; sz_mnm, sz_mxm, and sz_nbr are required input (all scalars)
; sz_grd, sz_min, sz_max, sz_dlt, sz_ctr are output arrays
; Usage: szgrd_mk,sz_mnm,sz_mxm,sz_nbr,grd_sng='logarithmic',szgrd=szgrd,sz_min=sz_min,sz_max=sz_max,sz_dlt=sz_dlt,sz_ctr=sz_ctr
; Based on mie_cls.cc:SzGrd::recompute()
if n_elements(grd_sng) eq 0 then grd_sng='linear' ; [sng] Type of size grid
if n_elements(sz_ctr) eq 0 then sz_ctr=fltarr(sz_nbr)
if n_elements(sz_dlt) eq 0 then sz_dlt=fltarr(sz_nbr)
if n_elements(sz_grd) eq 0 then sz_grd=fltarr(sz_nbr+1)
if n_elements(sz_max) eq 0 then sz_max=fltarr(sz_nbr)
if n_elements(sz_min) eq 0 then sz_min=fltarr(sz_nbr)
sz_min(0)=sz_mnm
if grd_sng eq 'logarithmic' then begin ; Space size bins logarithmically
    max_min_ratio=sz_mxm/sz_mnm
    series_ratio=max_min_ratio^(1.0/sz_nbr)
    if sz_mnm eq 0.0 then print,'ERROR: szgrd_mk() reports sz_mnm = 0.0 for size grid type = '+grd_sng
    for idx=1,sz_nbr-1 do begin
      sz_min(idx)=sz_min(idx-1)*series_ratio;
      sz_max(idx-1)=sz_min(idx);
    endfor ; end loop over idx
endif else if grd_sng eq 'linear' then begin ; Space size bins linearly
    sz_ncr=(sz_mxm-sz_mnm)/sz_nbr;
    for idx=0,sz_nbr-1 do begin
      sz_min(idx)=sz_mnm+idx*sz_ncr;
      sz_max(idx)=sz_mnm+(idx+1)*sz_ncr;
    endfor ; end loop over idx
endif else begin ; Space size bins linearly
	print,'ERROR: szgrd_mk() reports ',grd_sng,' is not a valid grid type'
endelse; end else grd_sng
sz_max(sz_nbr-1)=sz_mxm ; Ensure final max is not affected by roundoff
for idx=0,sz_nbr-1 do begin 
	sz_dlt(idx)=sz_max(idx)-sz_min(idx);
endfor; end loop over idx
for idx=0,sz_nbr-1 do begin
    sz_ctr(idx)=0.5*(sz_max(idx)+sz_min(idx));
    sz_grd(idx)=sz_min(idx);
endfor ; end loop over idx
sz_grd(sz_nbr)=sz_max(sz_nbr-1);
end; szgrd_mk()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End szgrd_mk()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin lgn_evl
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function lgn_evl,abc_ctr,mdn,gsd
; Purpose: Evaluate the lognormal size distribution for given parameters
; abc_ctr is required input array
; mdn and gsd are required input scalars
; lgn_dst is output array
; Usage: dst=lgn_evl(dmt_ctr,dmt_mdn,gsd,dst)
; Based on szdstlgn.F:lgn_evl()
; The distribution returned by this procedure is normalized
abc_nbr=n_elements(abc_ctr)
lgn_dst=fltarr(abc_nbr)
ln_gsd=alog(gsd)
lngsdsqrttwopi_rcp=1.0/(ln_gsd*sqrt(2.0*!pi))
tmp=alog(abc_ctr/mdn)/ln_gsd
xpn=exp(-0.5*tmp*tmp)
lgn_dst=lngsdsqrttwopi_rcp*xpn/abc_ctr
return,lgn_dst
end; lgn_evl()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End lgn_evl()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin ngl_dgr_2_m180_p180
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function ngl_dgr_2_m180_p180,ngl_arb
; Purpose: Convert arbitrary angle (in degrees) to range of -180 -- +180 degrees
while ngl_arb lt -180.0 do ngl_arb=ngl_arb+360.0
while ngl_arb gt +180.0 do ngl_arb=ngl_arb-360.0
return,ngl_arb
end; ngl_dgr_2_m180_p180()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End ngl_dgr_2_m180_p180()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin bytarr2sng
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function bytarr2sng,bytarr
; Purpose: Convert array of characters to string
; Useful for converting netCDF string attributes to strings
type_code=(size(bytarr))((size(bytarr))(0)+1) 
if type_code ne 7 then goto,err_exit
sng=''
byt_nbr=(size(bytarr))(1)
for byt_idx=0,byt_nbr-1 do begin
	sng=sng+bytarr[byt_idx]
endfor ; end loop over byt
return,sng
err_exit: foo=1
print,'ERROR: not a byte array in bytarr2sng()'
end; bytarr2sng()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End bytarr2sng()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin pnl_lbl_get
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function pnl_lbl_get,pnl_idx
; Purpose: Return appropriate alphabetical panel label
alphabet=['a','b','c','d','e','f','g','h','i']
;alphabet=['A','B','C','D','E','F','G','H','I']
return,'('+alphabet(pnl_idx)+')'
end; pnl_lbl_get()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End pnl_lbl_get()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin aeronet_unit_sng_get
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function aeronet_unit_sng_get,chn_idx
; Purpose: Return appropriate AERONET wavelength as string
; Indices are 0-based
wvl_sng=['1.02','0.87','0.67','0.5','0.44','0.38','0.34']
return,wvl_sng(chn_idx)+' !7l!5m'
;wvl_sng=['1020','870','670','500','440','380','340']
;return,wvl_sng(chn_idx)+' 5nm'
end; aeronet_unit_sng_get()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End aeronet_unit_sng_get()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin trc_nm_sbs
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function trc_nm_sbs,fld_sng,xpt_nm
; Purpose: Replace tracer placeholder string in field name with actual string for tracer based on experiment name
if xpt_nm eq 'dmr03' then trc_nm='!5(H!I2!NO)!I2!N'
if xpt_nm eq 'dmr04' then trc_nm='!5O!I2!N!9. !5O!I2!N + !5O!I2!N!9. !5N!I2!N'
if xpt_nm eq 'dmr05' then trc_nm='!5O!I2!N!9. !5N!I2!N'
if xpt_nm eq 'dmr06' then trc_nm='!5O!I2!N!9. !5O!I2!N'
if xpt_nm eq 'dmr08' then trc_nm='!5O!I2!N!9. !5O!I2!N + !5O!I2!N!9. !5N!I2!N'
dst_psn=strpos(xpt_nm,'dst')
if dst_psn ne -1 then trc_nm='Dust'
dmr_psn=strpos(fld_sng,'Tracer')
if dmr_psn ne -1 then begin
	head_sng=strmid(fld_sng,0,dmr_psn)
	tail_sng=strmid(fld_sng,dmr_psn+6,strlen(fld_sng))
	fld_sng=head_sng+trc_nm+tail_sng
endif; endif dmr_psn
return,fld_sng
end; trc_nm_sbs()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End trc_nm_sbs()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin mth2lbl
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function mth2lbl,mth_sng
; Purpose: Convert a monthe string, e.g., '07' or '1202', into a label string,
; 'Jul' or 'DJF', respectively
if mth_sng eq '1202' then lbl='DJF'
if mth_sng eq '0305' then lbl='MAM'
if mth_sng eq '0608' then lbl='JJA'
if mth_sng eq '0911' then lbl='SON'
if mth_sng eq '01' then lbl='Jan'
if mth_sng eq '02' then lbl='Feb'
if mth_sng eq '03' then lbl='Mar'
if mth_sng eq '04' then lbl='Apr'
if mth_sng eq '05' then lbl='May'
if mth_sng eq '06' then lbl='Jun'
if mth_sng eq '07' then lbl='Jul'
if mth_sng eq '08' then lbl='Aug'
if mth_sng eq '09' then lbl='Sep'
if mth_sng eq '10' then lbl='Oct'
if mth_sng eq '11' then lbl='Nov'
if mth_sng eq '12' then lbl='Dec'
return,lbl
end; end mth2lbl()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End mth2lbl()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin nint
; A copy of this should be in ~/.idle
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function nint,floater
if floater ge 0. then if floater-fix(floater) gt .5 then floater=floater+1 
if floater le 0. then if abs(floater-fix(floater)) gt .5 then floater=floater-1 
if abs(floater) > 32768 then return,long(floater) else return,fix(floater)
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End nint
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin auto_fmt
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function auto_fmt,nbr_to_fmt,dcm_plc_nbr,dbg_lvl=dbg_lvl

if n_elements(dbg_lvl) eq 0 then dbg_lvl=0

abs_nbr_to_fmt=abs(nbr_to_fmt)

; Leave room for decimal point and negative sign, if any
if nbr_to_fmt lt 0. then sgn_chr_nbr=1 else sgn_chr_nbr=0

; Definitions
; pfx_chr_nbr is number of characters preceding the decimal point
; mnt_chr_nbr is number of characters in the mantissa (i.e., following the decimal point) (same as dcm_plc_nbr)
; dcm_plc_nbr is number of characters following the decimal point
; sci_ntt_chr_nbr is number of characters used for scientific notation part (e.g., x 10^6)
; sgn_chr_nbr is number of characters needed for positive and negative signs (i.e., is 0 or 1)
; chr_ttl_nbr is total number of characters in formatted string
; dgt_nbr is number of digits in integer (used in integer formats only)

; IDL type_code definition:
; 0 Undefined
; 1 Byte
; 2 Integer
; 3 Longword integer
; 4 Floating-point
; 5 Double-precision floating
; 6 Complex floating
; 7 String
; 8 Structure

; Determine what we're dealing with
type_code=(size(nbr_to_fmt))((size(nbr_to_fmt))(0)+1) 

if ((type_code eq 4) or (type_code eq 5)) then begin
; Number is a float, more or less
; Pre-assign fmt_type just in case number falls through if statements 
fmt_type='f'
if abs_nbr_to_fmt ge 1.0e6 then begin
        pfx_chr_nbr=1
        sci_ntt_chr_nbr=4
        fmt_type='g'
endif else if ((abs_nbr_to_fmt lt 1.0e6) and (abs_nbr_to_fmt ge 1.0)) then begin
        pfx_chr_nbr=fix(alog10(round(abs_nbr_to_fmt)))+1
        sci_ntt_chr_nbr=0
        fmt_type='f'
endif else if ((abs_nbr_to_fmt lt 1.0) and (abs_nbr_to_fmt ge 1.0e-3)) then begin
; unfortunately, it seems impossible to get rid of leading zeroes before formatting
        pfx_chr_nbr=1
        sci_ntt_chr_nbr=0
        fmt_type='f'
endif else if (abs_nbr_to_fmt eq 0.0) then begin
        pfx_chr_nbr=1
        sci_ntt_chr_nbr=0
        fmt_type='f'
endif else begin 
        pfx_chr_nbr=1
        sci_ntt_chr_nbr=4
        fmt_type='g'
endelse

if fmt_type eq 'e' then $
        mnt_chr_nbr=dcm_plc_nbr $
else if fmt_type eq 'f' then $
        mnt_chr_nbr=dcm_plc_nbr $
else if fmt_type eq 'g' then $
        mnt_chr_nbr=dcm_plc_nbr

; Sum the characters in the numeric representation from left to right 
chr_ttl_nbr=sgn_chr_nbr+pfx_chr_nbr+1+dcm_plc_nbr+sci_ntt_chr_nbr
         
if chr_ttl_nbr lt 10 then fmt_sng=string(format='("(",A1,I1,".",I1,")")',fmt_type,chr_ttl_nbr,mnt_chr_nbr) else fmt_sng=string(format='("(",A1,I2,".",I1,")")',fmt_type,chr_ttl_nbr,mnt_chr_nbr)

endif else if ((type_code eq 2) or (type_code eq 3)) then begin
; Number is an integer
fmt_type='i'
if abs_nbr_to_fmt ne 0 then dgt_nbr=fix(alog10(abs_nbr_to_fmt))+1 else dgt_nbr=1 
	
chr_ttl_nbr=dgt_nbr+sgn_chr_nbr

if chr_ttl_nbr lt 10 then fmt_sng=string(format='("(",A1,I1,")")',fmt_type,chr_ttl_nbr) else fmt_sng=string(format='("(",A1,I2,")")',fmt_type,chr_ttl_nbr)

endif else begin; end if number is integer
; Number is currently unformattable
print,'unformattable number'
endelse;

the_sng=string(format=fmt_sng,nbr_to_fmt)
if strpos(the_sng,'*') ne -1 then begin
	print,'auto_fmt() reports nbr_to_fmt = ',nbr_to_fmt,' and fmt_sng = ',fmt_sng,' results in string '+string(format=fmt_sng,nbr_to_fmt)
	print,'chr_ttl_nbr = ',chr_ttl_nbr
	print,'mnt_chr_nbr = ',mnt_chr_nbr
    	print,'dcm_plc_nbr = ',dcm_plc_nbr
	print,'sci_ntt_chr_nbr = ',sci_ntt_chr_nbr
endif; endif dbg

return,fmt_sng
end; endif auto_fmt()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End auto_fmt
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin auto_sng
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function auto_sng,nbr_to_fmt,dcm_plc_nbr,dbg_lvl=dbg_lvl

if n_elements(dbg_lvl) eq 0 then dbg_lvl=0

if n_elements(nbr_to_fmt) eq 0 then return,''
if n_elements(nbr_to_fmt) eq 1 then begin
	the_sng=string(format=auto_fmt(nbr_to_fmt,dcm_plc_nbr),nbr_to_fmt)
	if strpos(the_sng,'*') ne -1 then begin
		print,'WARNING: auto_sng() reports ',nbr_to_fmt,' is formatted as ',the_sng
		print,'auto_sng() received:'
		print,'nbr_to_fmt = ',nbr_to_fmt
		print,'dcm_plc_nbr = ',dcm_plc_nbr
		print,'the_sng = ',the_sng
	endif; endif dbg
endif else begin
data_nbr=n_elements(nbr_to_fmt)
the_sng=strarr(data_nbr)
for idx=0,data_nbr-1 do begin
	the_sng(idx)=string(format=auto_fmt(nbr_to_fmt(idx),dcm_plc_nbr),nbr_to_fmt(idx))
endfor
endelse

if dbg_lvl ne 0 then begin
	print,'auto_sng() received:'
	print,'nbr_to_fmt = ',nbr_to_fmt
	print,'dcm_plc_nbr = ',dcm_plc_nbr
	print,'the_sng = ',the_sng
endif; endif dbg_lvl

return,the_sng
end; endif auto_sng()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End auto_sng
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin pth_nm_get
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro pth_nm_get,full_nm,path,name
path=''
name='False'
start_psn=0
while name eq 'False' do begin
	slash_psn=strpos(full_nm,'/',start_psn)
	next_slash_psn=strpos(full_nm,'/',slash_psn+1)
	if next_slash_psn eq -1 then begin
	; There is only one slash
		path=strmid(full_nm,0,slash_psn)
		name_lng=strlen(full_nm)-slash_psn-1	
		name=strmid(full_nm,slash_psn+1,name_lng)
	endif else begin
		dir_nm_lng=next_slash_psn-slash_psn
		path=path+strmid(full_nm,slash_psn,dir_nm_lng)
		start_psn=next_slash_psn
	endelse
endwhile
;
end; end pth_nm_get()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End pth_nm_get
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin xm
; A copy of this should be in ~/.idle but that file does not always get
; loaded, mysteriously enough.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro xm
xmanager
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End xm
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin open_ps
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro open_ps, $
	scale_fct=scale_fct, $
	fl_nm=fl_nm, $
	eps=eps, $
	ps=ps, $
	aspect=aspect, $
	x_sz=x_sz, $
	y_sz=y_sz, $
	color=color, $
	port=port, $
	land=land, $
	half=half, $
	full_col=full_col, $
	one_clm=one_clm, $
	two_col=two_col, $
	square=square

; Default plot to B&W postscript device with default dimensions
; IDL default is /portrait
; Caveats: 
; /port and /land refer to specific sizes on the page (as well as orientations)
; Use device,encapsulated=0 to turn off encapsulated postscript 
; Use device,font=0 to turn on postscript fonts, font=-1 to use default Hershey fonts, font=1 for truetype fonts

cmd='device'
;cmd=cmd+',font=0' ; Turn on Postscript fonts
if keyword_set(fl_nm) then cmd=cmd+',filename='''+fl_nm+'''
if keyword_set(scale_fct) then cmd=cmd+',scale_fct='+auto_sng(scale_fct,5)
if keyword_set(color) then cmd=cmd+',/color,bits_per_pixel=8'
if keyword_set(eps) then begin
	cmd=cmd+',/encapsulated,xoffset=0.0'
	if keyword_set(land) then cmd=cmd+',yoffset=10.0' else cmd=cmd+',yoffset=0.0'
; yoffsets should be same as xsizes of corresponding landscape plots.
endif; endif eps
if keyword_set(port) then cmd=cmd+',/portrait'
if keyword_set(square) then cmd=cmd+',/portrait'
if keyword_set(half) then cmd=cmd+',/portrait'
if keyword_set(full_col) then cmd=cmd+',/portrait'
if keyword_set(one_clm) then cmd=cmd+',/portrait'
if keyword_set(two_col) then cmd=cmd+',/portrait'
if keyword_set(aspect) then cmd=cmd+',/portrait'
if keyword_set(land) then cmd=cmd+',/landscape'
;if keyword_set(foo) then cmd=cmd+',/foo'

if keyword_set(x_sz) and keyword_set(y_sz) then begin
	if x_sz le 6.5 then cmd=cmd+',/portrait' else cmd=cmd+',/landscape'
	cmd=cmd+',xsize='+auto_sng(x_sz,5)+',ysize='+auto_sng(y_sz,5)+',/inches'
	y_off=(11.0-y_sz)/2.0
	x_off=(8.5-x_sz)/2.0
	if not keyword_set(eps) then cmd=cmd+',xoffset='+auto_sng(x_off,2)+',yoffset='+auto_sng(y_off,2)
endif else if keyword_set(aspect) then begin
; aspect ratio given in format x:y in decimal form
	x_sz=6.5
	y_sz=6.5/aspect
	cmd=cmd+',xsize='+auto_sng(x_sz,5)+',ysize='+auto_sng(y_sz,5)+',/inches'
	if not keyword_set(eps) then cmd=cmd+',xoffset=1.0,yoffset=1.0'
endif else if keyword_set(half) then begin
; 1.4:1 aspect ratio
	cmd=cmd+',xsize=7.,ysize=5.,/inches'
	if not keyword_set(eps) then cmd=cmd+',xoffset=1.0,yoffset=1.0'
endif else if keyword_set(full_col) then begin
; 2.5:1 aspect ratio
	cmd=cmd+',xsize=3.2,ysize=8,/inches'
	if not keyword_set(eps) then cmd=cmd+',xoffset=1.0,yoffset=1.0'
endif else if keyword_set(one_clm) then begin
; 1.25:1 aspect ratio
	cmd=cmd+',xsize=6.5,ysize=5.2,/inches'
	if not keyword_set(eps) then cmd=cmd+',xoffset=1.0,yoffset=1.0'
endif else if keyword_set(two_col) then begin
; 2.3:1 aspect ratio
	cmd=cmd+',xsize=6.5,ysize=2.82,/inches'
	if not keyword_set(eps) then cmd=cmd+',xoffset=1.0,yoffset=1.0'
endif else if keyword_set(square) then begin
; 1:1 aspect ratio
	cmd=cmd+',xsize=7.5,ysize=7.5,/inches'
	if not keyword_set(eps) then cmd=cmd+',xoffset=0.5,yoffset=2.'
endif else if keyword_set(port) then begin
	cmd=cmd+',xsize=7.5,ysize=10.0,/inches'
	if not keyword_set(eps) then cmd=cmd+',xoffset=0.5,yoffset=0.5'
endif else if keyword_set(land) then begin
	cmd=cmd+',xsize=10.0,ysize=7.5,/inches'
	if not keyword_set(eps) then cmd=cmd+',xoffset=0.5,yoffset=10.5'
endif else begin
; Default is landscape
	cmd=cmd+',xsize=10.0,ysize=7.5,/inches'
	if not keyword_set(eps) then cmd=cmd+',xoffset=0.5,yoffset=10.5,/landscape'
endelse

set_plot,'ps'
result=execute(cmd)
print,'open_ps(): cmd = ',cmd
if result ne 1 then print,'Error Printing...'

end; end open_ps()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End open_ps()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin close_ps
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro close_ps, $
	fl_nm=fl_nm, $
	printer=printer, $
	acdps4n=acdps4n, $
	dust=dust, $
	ghostview=ghostview, $
	hp4600=hp4600, $
	hp4600t=hp4600t, $
	lp=lp, $
	renoir=renoir, $
	renoirt=renoirt, $
	spectrum=spectrum, $
	tek360=tek360, $
	tek360T=tek360T, $
	tags=tags

; Script to save the current contents of the currently open graphical
; output file, as previously set by device,filename=fl_nm
; If the device is currently X windows, i.e., there is no graphical
; output file open, the script will look to see if the printer is 
; defined and, if so, print the given file to the printer.
; If the option is to send to tags, then the color 35mm slide options
; are set.
; If the option is to send to ghostview, then a ghostview process is
; opened on the host machine, from which the user may print anywhere.

if n_elements(fl_nm) eq 0 then fl_nm='idl.ps'
if keyword_set(acdps4n) then printer='acdps4n'
if keyword_set(dust) then printer='dust'
if keyword_set(lp) then printer='lp'
if keyword_set(renoir) then printer='renoir'
if keyword_set(renoirt) then printer='-o tray2 renoir'
if keyword_set(spectrum) then printer='spectrum'
if keyword_set(tek360) then printer='tek360'
if keyword_set(tek360T) then printer='tek360T'
if keyword_set(hp4600) then printer='hp4600'
if keyword_set(hp4600t) then printer='hp4600t'
;if keyword_set(foo) then cmd=cmd+',/foo'

; Close the Postscript file and reset the device to X Windows
if !d.name eq 'X' then begin
	print,'Device was already X so not saving current image'
	if ( $
	(n_elements(printer) eq 0) and $
	(not keyword_set(tags)) and $
	(not keyword_set(ghostview)) and $
	1 ) then return
endif else begin
	device,/close_file
	device,encapsulated=0
	set_plot,'X'
; 98/10/29: Backing store: 0 = None, 1 = Server-implemented, 2 = IDL-implemented
	device,retain=2
	print,'Saving image to '+fl_nm+' and returning to X'
endelse

; Print the file if requested
if n_elements(printer) ne 0 then begin
	cmd='lpr -P'+printer+' '+fl_nm
	spawn,cmd
	print,cmd
	cmd='lpq -P'+printer
	spawn,cmd
endif

if keyword_set(tags) then begin

; When using TAGS it is important to move the file to a temporary
; storage location and give it a unique name, because TAGS takes
; some time ( ~ 5-10 min. ?) to actually copy the file to downstairs.
; It is better that this is transparent to the user, so spawn a 
; dummy processs and capture it's id to generate a unique filename
; for the file to be sent to TAGS.
	usr_sng=getenv('LOGNAME')
	tmp_dir='/usr/tmp/'+usr_sng
	cmd='ls '+tmp_dir
	spawn,cmd,pid=p_id
	print,cmd
	tmp_fl_nm=tmp_dir+'/tags.'+strtrim(string(p_id),2)+'.ps'
	cmd='cp '+fl_nm+' '+tmp_fl_nm
	spawn,cmd
	print,cmd
	cmd='nrnet sendtg '+tmp_fl_nm+' r macr=slidesclsq req=ps mnt=yes'
	spawn,cmd
	print,cmd
	cmd='rsh shavano stattg'
	spawn,cmd
	print,cmd
endif

if keyword_set(ghostview) then begin
	cmd='ghostview -swap '+fl_nm
	spawn,cmd
	print,cmd
endif

end; end close_ps()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End close_ps()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Koontz
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; From velt@splinter.cgd.ucar.edu  Mon Aug 15 06:44:02 1994
; Received: from sage.cgd.ucar.EDU by ra.cgd.ucar.EDU (4.1/ NCAR Mail Server 04/10/90)
; 	id AA11111; Mon, 15 Aug 94 06:44:02 MDT
; Received: from splinter.rad.usf.edu by sage.cgd.ucar.EDU (8.6.4/ NCAR Mail Server 04/10/90)
; 	id GAA29803; Mon, 15 Aug 1994 06:44:00 -0600
; Received: from grater.rad (grater.rad.usf.edu [131.247.75.15]) by splinter.rad.usf.edu (8.6.5/8.6.5) with SMTP id IAA11378 for <zender@sage.cgd.ucar.edu>; Mon, 15 Aug 1994 08:43:57 -0400
; Received: by grater.rad (5.0/SMI-SVR4)
; 	id AA10233; Mon, 15 Aug 1994 08:43:39 +0500
; Date: Mon, 15 Aug 1994 08:43:39 +0500
; From: velt@splinter.rad.usf.edu (Robert Velthuizen)
; Message-Id: <9408151243.AA10233@grater.rad>
; To: zender@sage.cgd.ucar.EDU
; Subject: Re: high/low labels on contour plots
; Content-Length: 2044
; Status: RO
; 
; The routine below finds "peakness" in a 2D array; it looks for the relative
; height (depth) compared to immediate neighbors. If you're only interested
; in real local minima and maxima, do:
; 
; maxima = where( koontz( image) ge 8)
; minima = where( koontz( image) lt 1)
; 
; I have a bunch of variants of this routine, this one is for 2 dimensional
; images and distributes the  "peakness" over neighbors of equal height. You
; could also give peakness based on a random number. It may also be possible
; to improve on efficiency.
; 
; Hope this helps,
; Robert Velthuizen,
; Digital Medical Imaging Program of the
; H. Lee Moffitt Cancer Center and Research Institute at the
; University of South Florida.
; 
; I think bytscl is by far the simplest way to make "koontz" work with
; real data. The algorithm was designed to work on histogram data, which
; necessarily is positive integer. However, you could just do a type
; change in koontz to floats, and use min(fspace)-1 as instead of -1
; for the border, and select the seach space as everything > min(fspace).
; 
; If the neighborhood of the minima and
; maxima is very flat, then you may want to do something nonlinear to
; enhance the gradients around the candidate peaks. You could do that
; in two runs of koonz: first select the candidate peaks using bytscl,
; use the result to create a new image but with the neighborhood of
; each candidate extremum scaled between 0-255; the area's that are left
; out of this you set to 128 or so.
; 
; Good luck,
; Robert Velthuizen.
; ===============================================================================
 
Function Koontz, fspace
;+
; Name:		Koontz
;
; Purpose:	find the peaks in an image
;		written to find peaks in 2D histogram of feature space
;
; Reference:	A Khotanzad and A Bouarfa
;		Image segmentation by a parallel, non-parametric histogram
;		based clustering algorithm
;		Pattern Recognition, 23(9), pp 961-973, 1990
;-
d=size(fspace)
;
; The border can be [-1,min(fspace)-1] in which case the 
; search space should be [where(Im ge 0), where(Im ge min(fspace))]
;
border=-1
;border=min(fspace)-1
;
im=intarr(d(1)+1,d(2)+1)+border         ; +border, so boundary will never be a top
;im=fltarr(d(1)+1,d(2)+1)+border           ; +border, so boundary will never be a top
im(1:d(1),1:d(2))=fspace             ; Im holds the image in a border of +border's
Nm=n_elements(Im)
                                     ; N is the 3x3 neighborhood
N=reform([[shift(Im,-1,-1)],[shift(Im,-1, 0)],[shift(Im,-1, 1)],$
          [shift(Im, 0,-1)],[shift(Im, 0, 0)],[shift(Im, 0, 1)],$
          [shift(Im, 1,-1)],[shift(Im, 1, 0)],[shift(Im, 1, 1)]],Nm,9)
                                     ; Offsets for the neighborhood pixels
Offset=[1,1,1,0,0,0,-1,-1,-1] + [1,0,-1,1,0,-1,1,0,-1]*(d(1)+1)
Counts=0.0*Im
ss=where(Im ge 0)                    ; ss stands for searchspace
;ss=where(Im ge min(fspace))           ; ss stands for searchspace
Nm=n_elements(ss)
for i=0,Nm-1 do begin
   top=where(N(ss(i),*) eq max(N(ss(i),*)),nt)  ; Find all highest neighbors
   address=ss(i)+offset(top)                    ; pointer
   Counts(address) = Counts(address)+1.0/nt     ; distribute "peakness"
endfor
return,counts(1:d(1),1:d(2))
end
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Koontz
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Generic Power Spectra
; This routine is called by other routines. The calling routine supplies
; all of the input this routine needs to draw a nice power spectra.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro generic_pow_spec, $
	data=data, $
	time_nbr=time_nbr, $
	delta_time=delta_time, $
	x_ttl_sng=spec_x_ttl_sng, $
	y_ttl_sng=spec_y_ttl_sng
;
; Supply the defaults for lazy calling routines...
;
if n_elements(data) eq 0 then data=0
if n_elements(time_nbr) eq 0 then time_nbr=n_elements(data)
if n_elements(delta_time) eq 0 then delta_time=1.
if n_elements(x_ttl_sng) eq 0 then x_ttl_sng='Frequency !8f!5 (day!E-1!N)'
if n_elements(y_ttl_sng) eq 0 then y_ttl_sng='Normalized Power'
;if n_elements(foo) eq 0 then foo=0
;
; Analyze and display the spectral characteristics of a generic 
; time series.
;
; See discussion in Numerical Recipes in C p. 406,411 for description
; of fft output format.
;
; Zero frequency is stored in idx 0.
; Positive frequencies 0 < f < f_c correspond to 1 <= idx <= N/2 - 1
; Negative frequencies -f_c < f < 0 correspond to N/2 + 1 <= idx <= N-1
; Nyquist frequencies (both of them) are in idx = N/2
; Is that clear :-) ?
;
; I'd like to estimate the power spectral density in the N/2 + 1
; postive frequencies by the 'periodigram' method discussed in
; NR p. 439. The N/2 + 1 values of the periodigram estimated power
; have the property that their sum is equal to the mean squared 
; amplitude of the discrete sample of the original function (precip(t)).
;
; Correct normalization of the power depends upon convention.
; For the periodigram choice, dividing the power array by the sum
; of the squares of the input array (i.e., precip(t)) seems the 
; best choice.
;
; Note that IDL applies a 1/N normalization factor to the forward
; direction transform, so the 1/N^2 factor is already included 
; in the power spectrum after squaring the fft.
;
transform=fft(data(0:time_nbr-1),-1)
if time_nbr mod 2 ne 0 then print,'Array size = '+auto_sng(time_nbr,0)+' isn''t evenly divisible by 2 in fft().'
;
power=fltarr(1+time_nbr/2)
frequency=fltarr(1+time_nbr/2)
;
power(0)=abs(transform(0))^2
frequency(0)=0.
for idx=1,time_nbr/2-1 do begin
	frequency(idx)=idx/(time_nbr*delta_time)
	power(idx)=abs(transform(idx))^2+abs(transform(time_nbr-idx))^2
endfor
idx=time_nbr/2
frequency(idx)=idx/(time_nbr*delta_time)
power(idx)=abs(transform(idx))^2
;
; Normalize the power
;
power=time_nbr*power/total(data*data) ; NB: multiplying by time_nbr get the normalization of Parseval's theorem correct
;
if !d.name ne 'X' then chr_sz=2.0 else chr_sz=1.0

plot_io, $
	frequency*86400.0, $
	power, $
	xtit=x_ttl_sng, $
	ytit=y_ttl_sng, $
	xstyle=8, $
	ystyle=0, $
	xrange=[0.0,12.], $
	/ynozero, $
	thick=3.0, $
	xtick_get=x_tick_coords, $ ; [flt] Tick mark values (coordinates) at each tick mark
	charsize=chr_sz, $ ; [frc] Character size
	psym=-6, $ ; plot data as squares connected by lines
	ticklen=1, $
	xmargin=[8.75,1.5], $ ;[10,3] is [left,right] default
	ymargin=[3,2.5], $  ;[4,2] is [bottom,top] default
	linestyle=1
;
x_tck_nbr=n_elements(x_tick_coords)
xstring=strarr(x_tck_nbr)
for tick=0,x_tck_nbr-1 do begin
	period=x_tick_coords(tick)/86400. ; convert coord. to SI
	if period eq 0. then begin
		xstring(tick)='!9$!5'
	endif else begin
		period=1./period ; invert frequency to get period
		period=period/3600. ; convert period to hours
		xstring(tick)=auto_sng(period,1)
	endelse
endfor
;
axis, $
	xaxis=1, $
	xticks=x_tck_nbr-1, $ ; [nbr] Number of major intervals, 1 less than # of ticks
	xtickn=xstring, $
	xtitle='!5Period !7s!5 (hrs.)', $
	xmargin=[8.75,1.5], $ ;[10,3] is [left,right] default
	charsize=chr_sz ; [frc] Character size
;
end_of_procedure: foo=1
;
end; end generic_pow_spec()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Generic Power Spectra
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Generic Histogram
; This routine is called by other routines. The calling routine supplies
; all of the input this routine needs to draw a nice histogram.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro generic_hst, $
	abb_sng=abb_sng, $
	bin_sz=bin_sz, $
	chr_sz=chr_sz, $ ; [frc] Character size
	clean=clean, $
	data=data, $
	dbg_txt_sz=dbg_txt_sz, $
	hostname_sng=hostname_sng, $
	lat_max=lat_max, $
	lat_min=lat_min, $
	lon_max=lon_max, $
	lon_min=lon_min, $
	data_max=data_max, $
	data_min=data_min, $
	data_nbr=data_nbr, $
	rgn_sng=rgn_sng, $
	slice_sng=slice_sng, $
	sub_ttl_sng=sub_ttl_sng, $
	ttl_sng=ttl_sng, $
	sym_sng=sym_sng, $
	unit_sng=unit_sng, $
	usr_sng=usr_sng, $
	very_big_nbr=very_big_nbr, $
	rng_y=rng_y
;
; Supply the defaults for lazy calling routines...
;
if n_elements(bin_sz) eq 0 then bin_sz=1.
if n_elements(chr_sz) eq 0 then chr_sz=1.2
if n_elements(clean) eq 0 then clean=0
if n_elements(data) eq 0 then data=0
if n_elements(dbg_txt_sz) eq 0 then dbg_txt_sz=0.65
if n_elements(hostname_sng) eq 0 then hostname_sng=getenv('HOST')
if n_elements(lat_max) eq 0 then foo=1
if n_elements(lat_min) eq 0 then foo=1
if n_elements(lon_max) eq 0 then foo=1
if n_elements(lon_min) eq 0 then foo=1
if n_elements(data_max) eq 0 then data_max=max(data)
if n_elements(data_min) eq 0 then data_min=min(data)
if n_elements(data_nbr) eq 0 then data_nbr=n_elements(data)
if n_elements(rgn_sng) eq 0 then rgn_sng=''
if n_elements(slice_sng) eq 0 then slice_sng=''
if n_elements(sub_ttl_sng) eq 0 then sub_ttl_sng='Source: '
if n_elements(sym_sng) eq 0 then sym_sng=''
if n_elements(ttl_sng) eq 0 then ttl_sng=''
if n_elements(unit_sng) eq 0 then unit_sng=''
if n_elements(usr_sng) eq 0 then usr_sng=getenv('LOGNAME')
if n_elements(very_big_nbr) eq 0 then very_big_nbr=1.0e20

; find maxima after masked values have been removed
good_data=where(data lt very_big_nbr,count)
if count ne -1 then good_data=data(good_data) else good_data=data
max_good_data=max(good_data)
good_data_nbr=n_elements(good_data)

if n_elements(bin_sz) eq 0 then bin_sz=cntr_ntv
if n_elements(bin_max) eq 0 then bin_max=max_good_data
if n_elements(bin_min) eq 0 then bin_min=data_min

if bin_min ne data_min then begin
	outside=where(data lt bin_min,count)
	if count ne 0 then print,'points smaller than minimum bin: ',count
endif
if bin_max ne data_max then begin
	outside=where(data gt bin_max,count)
	if count ne 0 then print,'points larger than maximum bin: ',count
endif
;
hst=histogram(    $
		data, $
		binsize=bin_sz, $
		omax=hst_max, $
		omin=hst_min, $
		max=bin_max, $
		min=bin_min, $
		reverse_indices=inverse_data)
;
; Scale the histogram to frequency
frequency_axis=100.*hst/float(good_data_nbr)
;
print,'hst_min = ',hst_min
print,'hst_max = ',hst_max

bin_nbr=(size(hst))(1)
abc_axis=bin_sz*(findgen(bin_nbr)+.5)+hst_min
print,'There are ',bin_nbr,' bins covering from ',hst_min,' to ',hst_min+bin_nbr*bin_sz

cum_freq_axis=fltarr(bin_nbr)
cum_freq_axis(0)=frequency_axis(0)
for i=1,bin_nbr-1 do begin
	cum_freq_axis(i)=cum_freq_axis(i-1)+frequency_axis(i)
endfor

if n_elements(rng_y) le 1 then begin
	ord_min=min(ord)
	ord_max=max(ord)
endif else begin
	ord_min=rng_y(0)
	ord_max=rng_y(1)
endelse

plot, $
	abc_axis, $
	frequency_axis, $
	tit=ttl_sng, $
	xtit=slice_sng+' '+sym_sng+' ('+unit_sng+')', $
	ytit='Frequency (%)', $
	xrange=[0.0,40.0], $
	yrange=[ord_min,ord_max], $
	xstyle=0, $
	ystyle=9, $
	thick=3.0, $
	charsize=chr_sz, $ ; [frc] Character size
	xmargin=[7,7.5], $  ;[10,3] is [left,right] default
	ymargin=[3.2,2], $  ;[4,2] is [bottom,top] default
	psym=10

; Add the first half-bin
plots, $
	[hst_min, $
	hst_min, $
	abc_axis(0)], $
	[0.0, $
	frequency_axis(0), $
	frequency_axis(0)], $
	linestyle=0, $
	thick=3.0, $
	clip=[!x.crange(0),!y.crange(0),!x.crange(1),!y.crange(1)], $
	/data

plots, $
	[hst_min+bin_nbr*bin_sz, $
	hst_min+bin_nbr*bin_sz, $
	abc_axis(bin_nbr-1)], $
	[0.0, $
	frequency_axis(bin_nbr-1), $
	frequency_axis(bin_nbr-1)], $
	linestyle=0, $
	thick=3.0, $
	clip=[!x.crange(0),!y.crange(0),!x.crange(1),!y.crange(1)], $
	/data

plot, $
	abc_axis, $
	cum_freq_axis, $
	xstyle=4, $
	ystyle=4, $
;	/ynozero, $
	thick=3.0, $
	charsize=chr_sz, $ ; [frc] Character size
	linestyle=2, $
	xmargin=[7,7.5], $ ;[10,3] is [left,right] default
	ymargin=[3.2,2], $  ;[4,2] is [bottom,top] default
	/noerase

axis, $
	yaxis=1, $
	ystyle=1, $
	ytitle='!5Cumulative Frequency (%)', $
	charsize=chr_sz ; [frc] Character size

if clean then goto,end_of_procedure

possible_elements_sng=auto_sng(data_nbr,0)
plotted_elements_sng=auto_sng(round(total(hst)),0)
bin_sz_sng=auto_sng(bin_sz,1)
bin_sz_sng='!5Bin size = '+bin_sz_sng+' '+unit_sng+'.'
avg_sng=auto_sng(total(good_data)/good_data_nbr,1)
avg_sng='!5Average '+sym_sng+' = '+avg_sng+' '+unit_sng+'.'
mdn_sng=auto_sng(median(data),1)
mdn_sng='!5Median '+sym_sng+' = '+mdn_sng+' '+unit_sng+'.'

if $
	n_elements(lat_max) eq 0 or $
	n_elements(lat_min) eq 0 or $
	n_elements(lon_max) eq 0 or $
	n_elements(lon_min) eq 0 then begin

lon_sng=''	
lat_sng=''	

endif else begin; endif region is not defined

lat_min_sng=auto_sng(lat_min,1)
lat_max_sng=auto_sng(lat_max,1)
lat_sng='!5'+lat_min_sng+'!E!12_!5!N < !7u!5 < '+ $
		lat_max_sng+'!E!12_!5!N'
lon_min_sng=auto_sng(lon_min,1)
lon_max_sng=auto_sng(lon_max,1)
lon_sng='!5'+lon_min_sng+'!E!12_!5!N < !7k!5 < '+ $
		lon_max_sng+'!E!12_!5!N'

endelse; endif region is defined

coord_sng=' in '+rgn_sng+':'+lon_sng+', '+lat_sng

xyouts,0.98,0.025,sub_ttl_sng+' '+sym_sng+coord_sng,alignment=0.0,orientation=90.0,size=dbg_txt_sz,/NORMAL

xyouts,1.0,0.025,'Used '+plotted_elements_sng+'/'+possible_elements_sng+' points.  '+bin_sz_sng+'  '+avg_sng+'  '+mdn_sng,alignment=0.0,orientation=90.0,size=dbg_txt_sz,/NORMAL

sys_time=systime(0) 
xyouts,1.0,0.0,usr_sng+'@'+hostname_sng+' '+sys_time,size=dbg_txt_sz,orientation=0.0,alignment=1.0,/NORMAL

end_of_procedure: foo=1

end; end generic_hst()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Generic Histogram
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Generic Scatterplot 
; This routine is called by other routines. The calling routine supplies
; all of the input this routine needs to draw a nice scatterplot.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro generic_scat, $
	abb_sng=abb_sng, $
	bin_sz=bin_sz, $
	chr_sz=chr_sz, $ ; [frc] Character size
	clean=clean, $
	data=data, $
	dbg_txt_sz=dbg_txt_sz, $
	hostname_sng=hostname_sng, $
	info=info, $
	lat_max=lat_max, $
	lat_min=lat_min, $
	lon_max=lon_max, $
	lon_min=lon_min, $
	data_max=data_max, $
	data_min=data_min, $
	data_nbr=data_nbr, $
	rgn_sng=rgn_sng, $
	slice_sng=slice_sng, $
	sub_ttl_sng=sub_ttl_sng, $
	sym_sng=sym_sng, $
	unit_sng=unit_sng, $
	usr_sng=usr_sng, $
	very_big_nbr=very_big_nbr

; Supply defaults for lazy calling routines...
if n_elements(abb_sng) eq 0 then abb_sng=['','']
if n_elements(bin_sz) eq 0 then bin_sz=[1.0,1.0]
if n_elements(chr_sz) eq 0 then chr_sz=1.2
if n_elements(clean) eq 0 then clean=0
if n_elements(data) eq 0 then data=[0,0]
if n_elements(dbg_txt_sz) eq 0 then dbg_txt_sz=0.65
if n_elements(hostname_sng) eq 0 then hostname_sng=getenv('HOST')
if n_elements(info) eq 0 then info=1
if n_elements(lat_max) eq 0 then foo=1
if n_elements(lat_min) eq 0 then foo=1
if n_elements(lon_max) eq 0 then foo=1
if n_elements(lon_min) eq 0 then foo=1
if n_elements(data_max) eq 0 then data_max=[max(data(0,*)),max(data(1,*))]
if n_elements(data_min) eq 0 then data_min=[min(data(0,*)),min(data(1,*))]
if n_elements(data_nbr) eq 0 then data_nbr=[n_elements(data(0,*)),n_elements(data(1,*))]
if n_elements(rgn_sng) eq 0 then rgn_sng=''
if n_elements(slice_sng) eq 0 then slice_sng=['','']
if n_elements(sub_ttl_sng) eq 0 then sub_ttl_sng=['Source: ','Source: ']
if n_elements(sym_sng) eq 0 then sym_sng=['','']
if n_elements(unit_sng) eq 0 then unit_sng=['','']
if n_elements(usr_sng) eq 0 then usr_sng=getenv('LOGNAME')
if n_elements(very_big_nbr) eq 0 then very_big_nbr=1.0e20

; Make sure there are two equally sized fields selected
if data_nbr(0) eq 0 then begin
	print,'Must save a field first'
	goto, end_of_procedure
endif
if data_nbr(0) ne data_nbr(1) then begin
	print,'Fields are of differing sizes, check your data...'
	goto, end_of_procedure
endif

; Find size of good_data arrays
good_idx=where(data(0,*) lt very_big_nbr,good_data_nbr)
good_data=[data(0,good_idx),data(1,good_idx)]
good_idx=where(good_data(1,*) lt very_big_nbr,good_data_nbr)
good_data=[good_data(0,good_idx),good_data(1,good_idx)]

plot, $
	good_data(0,*), $
	good_data(1,*), $
	max_value=very_big_nbr, $
	tit='Least squares fit: '+abb_sng(1)+' & '+abb_sng(0), $
	xtit=slice_sng(0)+' '+sym_sng(0)+' ('+unit_sng(0)+')', $
	ytit=slice_sng(1)+' '+sym_sng(1)+' ('+unit_sng(1)+')', $
	xstyle=0, $
	ystyle=0, $
	thick=2.0, $
	charsize=chr_sz, $ ; [frc] Character size
	xmargin=[7,2], $ ;[10,3] is [left,right] default
	ymargin=[3.1,2], $  ;[4,2] is [bottom,top] default
	psym=1

; Perform statistical analysis of correlation between two datasets
fit_coeffs=poly_fit( $
	good_data(0,*), $
	good_data(1,*), $
	1, $
	fitted_ords)
;
oplot, $
	[data_min(0), $
	data_max(0)], $
	[fit_coeffs(0)+data_min(0)*fit_coeffs(1), $
	fit_coeffs(0)+data_max(0)*fit_coeffs(1)], $
	thick=3.0, $
	linestyle=0	
;
if $
	n_elements(lat_max) eq 0 or $
	n_elements(lat_min) eq 0 or $
	n_elements(lon_max) eq 0 or $
	n_elements(lon_min) eq 0 then begin
;
lon_sng=''	
lat_sng=''	
;
endif else begin; endif region is not defined
;
lat_min_sng=auto_sng(lat_min,1)
lat_max_sng=auto_sng(lat_max,1)
lat_sng='!5'+lat_min_sng+'!E!12_!5!N < !7u!5 < '+ $
		lat_max_sng+'!E!12_!5!N'
lon_min_sng=auto_sng(lon_min,1)
lon_max_sng=auto_sng(lon_max,1)
lon_sng='!5'+lon_min_sng+'!E!12_!5!N < !7k!5 < '+ $
		lon_max_sng+'!E!12_!5!N'
;
endelse; endif region is defined
;
if info then begin
coord_sng=' in '+rgn_sng+':'+lon_sng+', '+lat_sng
;
if not clean then begin
	xyouts,0.98,0.025,sub_ttl_sng(1)+' '+sym_sng(1)+' !8vs.!5 '+sub_ttl_sng(0)+' '+sym_sng(0)+coord_sng,alignment=0.0,orientation=90.0,size=dbg_txt_sz,/NORMAL
	sys_time=systime(0) 
	xyouts,1.0,0.0,usr_sng+'@'+hostname_sng+' '+sys_time,size=dbg_txt_sz,orientation=0.0,alignment=1.0,/NORMAL
endif; endif not clean
endif; endif info
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin linear (Pearson) regression
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pearson_mtx=correl_matrix(good_data)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End linear (Pearson) regression
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Non-parametric (Spearman) correlation
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
matrix=good_data
spearman,matrix,spearman_mtx
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Non-parametric (Spearman) correlation
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
print,'b_term = ',fit_coeffs(0)
print,'m_term = ',fit_coeffs(1)
b_sng=auto_sng(abs(fit_coeffs(0)),2)
m_sng=auto_sng(fit_coeffs(1),2)
if fit_coeffs(0) ge 0. then sign_sng=' + ' else sign_sng=' - '
fit_sng=sym_sng(1)+' = '+m_sng+' '+sym_sng(0)+sign_sng+b_sng
pearson_sng=string(format='(F5.2)',pearson_mtx(0,1))
pearson_sng='!8r!5 = '+pearson_sng
spearman_sng=string(format='(F5.2)',spearman_mtx(0,1))
spearman_sng='!8r!Is!N!5 = '+spearman_sng
;
if not clean then begin
xyouts,1.0,0.025,'Least squares fit: '+fit_sng+'; Linear '+pearson_sng+', Non-parametric '+spearman_sng,size=dbg_txt_sz,orientation=90.0,alignment=0.0,/NORMAL
endif; endif not clean
xyouts,0.25,0.80,'Correlation '+pearson_sng,size=2.0,/NORMAL
xyouts,0.25,0.70,'Slope = '+m_sng,size=2.0,/NORMAL
;xyouts,0.2,0.75,fit_sng,size=2.0,orientation=90.0,alignment=0.0,/NORMAL
;xyouts,0.2,0.75,pearson_sng,size=2.0,orientation=90.0,alignment=0.0,/NORMAL
;xyouts,0.2,0.70,spearman_sng,size=2.0,orientation=90.0,alignment=0.0,/NORMAL
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Full multiple regression
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;matrix=fltarr(1,good_data_nbr)
;matrix(0,*)=good_data(0,*)
;weights=replicate(1,good_data_nbr)
;m_term=regress( $
;		matrix, $
;		reform(good_data(1,*),good_data_nbr), $
;		weights, $
;		/relative_weight, $
;		fitted_ords, $
;		b_term, $
;		sigma, $
;		f_test, $
;		pearson, $
;		multi_pearson, $
;		chi_squared)
;
;print,'b_term = ',b_term
;print,'m_term = ',m_term
;print,'sigma = ',sigma
;print,'f_test = ',f_test
;print,'pearson = ',pearson
;print,'multi_pearson = ',multi_pearson
;print,'chi_squared = ',chi_squared
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Full multiple regression
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
end_of_procedure: foo=1
;
end; end generic_scat()
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Generic Scatterplot
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Generic Correlation 
; This routine is called by other routines. The calling routine supplies
; all of the input this routine needs to draw a nice correlation plot
; of the input data.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro corr_tester,max_mag_lag=max_mag_lag
;
; This tester should show that the first curve, sin(t), lags
; the second curve, sin(t+pi/4). Since the first curve is 
; lagging, the maximum corr() occurs for positive tau, namely
; at tau = pi/4 = .7854
;
if n_elements(max_mag_lag) eq 0 then max_mag_lag=50
data=fltarr(2,100)
data(0,*)=sin(findgen(100)*2.*!pi/99.)
data(1,*)=sin(.7854+findgen(100)*2.*!pi/99.)
delta_time=2.*!pi/99.
plot,data(0,*),linestyle=0
oplot,data(1,*),linestyle=1
print,' Hit any key to continue...'
junk = get_kbrd(1)
generic_corr, $
	data=data(0:1,*), $
	delta_time=delta_time, $
	max_mag_lag=max_mag_lag
end

pro generic_corr, $
	abb_sng=abb_sng, $
	data=data, $
	dbg_txt_sz=dbg_txt_sz, $
	delta_time=delta_time, $
	hostname_sng=hostname_sng, $
	lat_max=lat_max, $
	lat_min=lat_min, $
	lon_max=lon_max, $
	lon_min=lon_min, $
	max_mag_lag=max_mag_lag, $
	data_max=data_max, $
	data_min=data_min, $
	data_nbr=data_nbr, $
	rgn_sng=rgn_sng, $
	slice_sng=slice_sng, $
	sub_ttl_sng=sub_ttl_sng, $
	sym_sng=sym_sng, $
	unit_sng=unit_sng, $
	usr_sng=usr_sng, $
	very_big_nbr=very_big_nbr

; Supply defaults for lazy calling routines...
if n_elements(abb_sng) eq 0 then abb_sng=['','']
if n_elements(data) eq 0 then data=[0,0]
if n_elements(dbg_txt_sz) eq 0 then dbg_txt_sz=0.65
if n_elements(delta_time) eq 0 then delta_time=1.
if n_elements(hostname_sng) eq 0 then hostname_sng=getenv('HOST')
if n_elements(lat_max) eq 0 then foo=1
if n_elements(lat_min) eq 0 then foo=1
if n_elements(lon_max) eq 0 then foo=1
if n_elements(lon_min) eq 0 then foo=1
if n_elements(data_max) eq 0 then data_max=[max(data(0,*)),max(data(1,*))]
if n_elements(data_min) eq 0 then data_min=[min(data(0,*)),min(data(1,*))]
if n_elements(data_nbr) eq 0 then data_nbr=[n_elements(data(0,*)),n_elements(data(1,*))]
if n_elements(rgn_sng) eq 0 then rgn_sng=''
if n_elements(slice_sng) eq 0 then slice_sng=['','']
if n_elements(sub_ttl_sng) eq 0 then sub_ttl_sng=['Source: ','Source: ']
if n_elements(sym_sng) eq 0 then sym_sng=['','']
if n_elements(unit_sng) eq 0 then unit_sng=['','']
if n_elements(usr_sng) eq 0 then usr_sng=getenv('LOGNAME')
if n_elements(very_big_nbr) eq 0 then very_big_nbr=1.0e20

; Make sure there are two equally sized fields selected
if data_nbr(0) eq 0 then begin
	print,'Must save a field first'
	goto, end_of_procedure
endif
if data_nbr(0) ne data_nbr(1) then begin
	print,'Fields are of differing sizes, check your data...'
	goto, end_of_procedure
endif

; Time series data probably has no blocked values
; Spatial data, however, may have blocked values
; Find the size of the good_data arrays

good_idx=where(data(0,*) lt very_big_nbr,good_data_nbr)
good_data=[data(0,good_idx),data(1,good_idx)]
good_idx=where(good_data(1,*) lt very_big_nbr,good_data_nbr)
good_data=[good_data(0,good_idx),good_data(1,good_idx)]

; Correlations are very simply computed in spectral space, see
; Numerical Recipes p. 433 for discussion. 

; The correlation (at a lag T) of two time series is the 
; summation of the product of the first time series offset by 
; T elements with the second time series.

; The spectral method gets rid of the summation using 
; the correlation theorem: 

; The correlation (at a lag T) of two time series is the 
; Tth element of the (sequentially ordered) inverse transform 
; of the product of the forward transform of the first time series 
; with the complex conjugate of the forward transform of the 
; second time series. 

; The only numerical snag is that these time series aren't periodic,
; so they need to be zero padded in order to handle end effects
; correctly. From NR, "If you are interested in the correlation
; for lags as large as +/- K, then you must append a buffer zone of
; K zeros at the end of both input data sets." This is simply to
; provide a large enough defined data window to perform the 
; correlation out to the specified lag.

; M is the size of the (unpadded) raw data series.
; N is the size of the (padded) fft input data series.
; K = the number of lags to conmpute in each direction.
; The maximum possible K is M.
; The default K in this routine is M/2.
; N=M+K
; Map (unscramble) the size N=M+K output array into a size 2K+1 array. 
; i.e., when K=M/2 map the correlations into a size N+1 array.

; Put the negative lags into the first half of the lag array.
; The sequentially ordered (unscrambled) array contains the 
; correlations ranging from negative lag M/2 to positive lag M/2. 
; The disparity in sizes between the unsequential and sequentially
; ordered arrays is common to fft problems: one of the elements of 
; the sequentially ordered array is a "linearly dependent" value
; because positive and negative lag M/2 are, of course, identical.
; (check this).

; Unscrambling:
; Increasing output order when K is arbitrary, input is:
; Neg. lag K will be the 1 element of uns_corr, uns_corr(0)
; Neg. lag K-1 will be the 2 element of uns_corr, uns_corr(1)
; Neg. lag 1 will be the K element of uns_corr, uns_corr(K-1)
; Zero lag will be the K+1 (midpoint) element of uns_corr, uns_corr(K)
; Pos. lag 1 will be the K+2 element of uns_corr, uns_corr(K+1)
; Pos. lag K-1 will be the 2K element of uns_corr, uns_corr(2K-1)
; Pos. lag K will be the 2K+1 element of uns_corr, uns_corr(2K)

; Increasing input order when K is arbitrary, output is:
; Zero lag is the first element of corr, corr(0)
; Pos. lag 1 is the 2 element of corr, corr(1)
; Pos. lag K-1 is the K element of corr, corr(K-1)
; Pos. lag K is the K+1 element of corr, corr(K)
; Neg. lag K is the N-K+1 element of corr, corr(N-K)
; Neg. lag K-1 is the N-K+2 element of corr, corr(N-K+1)
; Neg. lag 2 is the N-1 element of corr, corr(N-2)
; Neg. lag 1 is the N element of corr, corr(N-1)

; Increasing output order when K=M and N=2M, input is:
; Neg. lag M will be the 1 element of uns_corr, uns_corr(0)
; Neg. lag M-1 will be the 2 element of uns_corr, uns_corr(1)
; Neg. lag 1 will be the M element of uns_corr, uns_corr(M-1)
; Zero lag will be the M+1 (midpoint) element of uns_corr, uns_corr(M)
; Pos. lag 1 will be the M+2 element of uns_corr, uns_corr(M+1)
; Pos. lag M-1 will be the 2M element of uns_corr, uns_corr(2M-1)
; Pos. lag M will be the 2M+1 element of uns_corr, uns_corr(2M)

; Increasing input order when K=M and N=2M, output is:
; Zero lag is the first element of corr, corr(0)
; Pos. lag 1 is the 2 element of corr, corr(1)
; Pos. lag M-1 is the M element of corr, corr(M-1)
; Pos. lag M is the M+1 element of corr, corr(M) = Neg. lag M
; Neg. lag M is the M+1 element of corr, corr(M) = Pos. lag M
; Neg. lag M-1 is the M+2 element of corr, corr(M+1)
; Neg. lag 1 is the 2M element of corr, corr(2M-1)

time_nbr=good_data_nbr		; time_nbr is even, hopefully
M=time_nbr			; M is even, hopefully
if n_elements(max_mag_lag) eq 0 then K=floor(M/2) else K=max_mag_lag
N=M+K
nbr_lag=2*K+1			; i.e., nbr_lag is odd

if M mod 2 ne 0 then print,'Array size M = '+auto_sng(M,0)+' isn''t evenly divisible by 2 in fft().'
if N mod 2 ne 0 then print,'Array size N = '+auto_sng(N,0)+' isn''t evenly divisible by 2 in fft().'

data_0=reform(data(0,*))
data_0=[data_0,fltarr(K)]
data_1=reform(data(1,*))
data_1=[data_1,fltarr(K)]

; Note that IDL applies a 1/N normalization factor to the forward
; direction transform, so the 1/N^2 factor is already included 
; in the power spectrum after squaring the fft.

corr=fft(fft(data_0,-1)*conj(fft(data_1,-1)),1)

lag=fltarr(nbr_lag)
uns_corr=fltarr(nbr_lag)
for idx=0,K-1 do begin
	uns_corr(idx)=corr(N-K+idx)
	uns_corr(K+idx+1)=corr(idx+1)
endfor
uns_corr(K)=corr(0)

lag=(indgen(nbr_lag)-K)*delta_time

if !d.name ne 'X' then chr_sz=2.0 else chr_sz=1.0
plot, $
	lag, $
	uns_corr, $
	max_value=very_big_nbr, $
	tit='Lag Corr. of '+abb_sng(0)+' and '+abb_sng(1), $
	xtit='Lag !7s!5 (min.)', $
	ytit='Corr('+slice_sng(0)+' '+sym_sng(0)+', '+slice_sng(1)+' '+sym_sng(1)+')', $
	xstyle=1, $
	ystyle=1, $
	xmargin=[7,2], $ ;[10,3] is [left,right] default
	ymargin=[3,2], $  ;[4,2] is [bottom,top] default
	thick=3.0, $
	charsize=chr_sz, $ ; [frc] Character size
	linestyle=0

if $
	n_elements(lat_max) eq 0 or $
	n_elements(lat_min) eq 0 or $
	n_elements(lon_max) eq 0 or $
	n_elements(lon_min) eq 0 then begin

lon_sng=''	
lat_sng=''	

endif else begin; endif region is not defined

lat_min_sng=auto_sng(lat_min,1)
lat_max_sng=auto_sng(lat_max,1)
lat_sng='!5'+lat_min_sng+'!E!12_!5!N < !7u!5 < '+ $
		lat_max_sng+'!E!12_!5!N'
lon_min_sng=auto_sng(lon_min,1)
lon_max_sng=auto_sng(lon_max,1)
lon_sng='!5'+lon_min_sng+'!E!12_!5!N < !7k!5 < '+ $
		lon_max_sng+'!E!12_!5!N'

endelse; endif region is defined

coord_sng=' in '+rgn_sng+':'+lon_sng+', '+lat_sng

xyouts,0.98,0.025,sub_ttl_sng(1)+' '+sym_sng(1)+' !8vs.!5 '+sub_ttl_sng(0)+' '+sym_sng(0)+coord_sng,alignment=0.0,orientation=90.0,size=dbg_txt_sz,/NORMAL

sys_time=systime(0) 
xyouts,1.0,0.0,usr_sng+'@'+hostname_sng+' '+sys_time,size=dbg_txt_sz,orientation=0.0,alignment=1.0,/NORMAL

end_of_procedure: foo=1

end; end generic_corr()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Generic Correlation
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Generic Multi-plot 
; This routine is called by other routines. The calling routine supplies
; all of the input this routine needs to draw multiple ords for
; the same abscissa.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro generic_multi, $
	abb_sng=abb_sng, $
	bin_sz=bin_sz, $
	data=data, $
	dbg_txt_sz=dbg_txt_sz, $
	hostname_sng=hostname_sng, $
	lat_max=lat_max, $
	lat_min=lat_min, $
	lon_max=lon_max, $
	lon_min=lon_min, $
	data_max=data_max, $
	data_min=data_min, $
	data_nbr=data_nbr, $
	rgn_sng=rgn_sng, $
	slice_sng=slice_sng, $
	sub_ttl_sng=sub_ttl_sng, $
	sym_sng=sym_sng, $
	unit_sng=unit_sng, $
	usr_sng=usr_sng, $
	very_big_nbr=very_big_nbr
;
; Supply the defaults for lazy calling routines...
;
if n_elements(abb_sng) eq 0 then abb_sng=['','']
if n_elements(bin_sz) eq 0 then bin_sz=[1.0,1.0]
if n_elements(data) eq 0 then data=[0,0]
if n_elements(dbg_txt_sz) eq 0 then dbg_txt_sz=0.65
if n_elements(hostname_sng) eq 0 then hostname_sng=getenv('HOST')
if n_elements(lat_max) eq 0 then foo=1
if n_elements(lat_min) eq 0 then foo=1
if n_elements(lon_max) eq 0 then foo=1
if n_elements(lon_min) eq 0 then foo=1
if n_elements(data_max) eq 0 then data_max=[max(data(0,*)),max(data(1,*))]
if n_elements(data_min) eq 0 then data_min=[min(data(0,*)),min(data(1,*))]
if n_elements(data_nbr) eq 0 then data_nbr=[n_elements(data(0,*)),n_elements(data(1,*))]
if n_elements(rgn_sng) eq 0 then rgn_sng=''
if n_elements(slice_sng) eq 0 then slice_sng=['','']
if n_elements(sub_ttl_sng) eq 0 then sub_ttl_sng=['Source: ','Source: ']
if n_elements(sym_sng) eq 0 then sym_sng=['','']
if n_elements(unit_sng) eq 0 then unit_sng=['','']
if n_elements(usr_sng) eq 0 then usr_sng=getenv('LOGNAME')
if n_elements(very_big_nbr) eq 0 then very_big_nbr=1.0e20

; Make sure there is data to plot
if data_nbr(0) eq 0 then begin
	print,'Error: no data to plot'
	goto, end_of_procedure
endif

; Find size of good_data arrays
dmn_nbr=(size(data))(0)

if !d.name ne 'X' then chr_sz=2.0 else chr_sz=1.0
plot, $
	abc(0,*), $
	good_data(0,*), $
	charsize=chr_sz, $ ; [frc] Character size
	max_value=very_big_nbr, $
	thick=2.0, $
	tit=abb_sng(0)+' !8vs.!5 '+abb_sng(1), $
	xmargin=[7,2], $ ;[10,3] is [left,right] default
	xstyle=1, $
	xtit=slice_sng(0)+' '+sym_sng(0)+' ('+unit_sng(0)+')', $
	ymargin=[3,2], $ ;[4,2] is [bottom,top] default
	ystyle=1, $
	ytit=slice_sng(1)+' '+sym_sng(1)+' ('+unit_sng(1)+')'
;
oplot, $
	good_data(1,*), $
	thick=3.0, $
	linestyle=0	
;
if $
	n_elements(lat_max) eq 0 or $
	n_elements(lat_min) eq 0 or $
	n_elements(lon_max) eq 0 or $
	n_elements(lon_min) eq 0 then begin
;
lon_sng=''	
lat_sng=''	
;
endif else begin; endif region is not defined
;
lat_min_sng=auto_sng(lat_min,1)
lat_max_sng=auto_sng(lat_max,1)
lat_sng='!5'+lat_min_sng+'!E!12_!5!N < !7u!5 < '+ $
		lat_max_sng+'!E!12_!5!N'
lon_min_sng=auto_sng(lon_min,1)
lon_max_sng=auto_sng(lon_max,1)
lon_sng='!5'+lon_min_sng+'!E!12_!5!N < !7k!5 < '+ $
		lon_max_sng+'!E!12_!5!N'
;
endelse; endif region is defined
;
coord_sng=' in '+rgn_sng+':'+lon_sng+', '+lat_sng
;
xyouts,0.98,0.025,sub_ttl_sng(1)+' '+sym_sng(1)+' !8vs.!5 '+sub_ttl_sng(0)+' '+sym_sng(0)+coord_sng,alignment=0.0,orientation=90.0,size=dbg_txt_sz,/NORMAL
;
sys_time=systime(0) 
xyouts,1.0,0.0,usr_sng+'@'+hostname_sng+' '+sys_time,size=dbg_txt_sz,orientation=0.0,alignment=1.0,/NORMAL
;
end_of_procedure: foo=1
;
end; end generic_multi()
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Generic Multi-plot
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin cntr_lbl_mk()
; This routine computes the labels and decorations for a set of contour levels
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro cntr_lbl_mk,cntr_lvl,cntr_lvl_min,cntr_ntv,cntr_lvl_nbr,unit_pfx_sng, $
	cntr_fll_idx,cntr_lvl_lbl,cntr_ln_sty,cntr_thk,cntr_which_lbl,dbg_lvl=dbg_lvl
if n_elements(dbg_lvl) eq 0 then dbg_lvl=0

dcm_plc_nbr=0
unit_pfx_sng=''
p_fmt_fct=1
max_abs_cntr_lvl=max(abs(cntr_lvl))

if max_abs_cntr_lvl le 1.0e-5 then begin
	dcm_plc_nbr=1
	p_fmt_fct=1.0e6
	unit_pfx_sng='!9X!5 10!E-6!N '
endif else if max_abs_cntr_lvl le 1.0e-4 then begin
	if cntr_ntv mod 0.0000025 eq 0.0 then if cntr_ntv mod 0.0000050 ne 0.0 then dcm_plc_nbr=2 else dcm_plc_nbr=1
	p_fmt_fct=1.0e5
	unit_pfx_sng='!9X!5 10!E-5!N '
endif else if max_abs_cntr_lvl le 1.0e-3 then begin
	if cntr_ntv mod 0.000025 eq 0.0 then if cntr_ntv mod 0.000050 ne 0.0 then dcm_plc_nbr=2 else dcm_plc_nbr=1
	p_fmt_fct=1.0e4
	unit_pfx_sng='!9X!5 10!E-4!N '
endif else if max_abs_cntr_lvl le 1.0e-2 then begin
	if cntr_ntv mod 0.00025 eq 0.0 then if cntr_ntv mod 0.00050 ne 0.0 then dcm_plc_nbr=2 else dcm_plc_nbr=1
	p_fmt_fct=1.0e3
	unit_pfx_sng='!9X!5 10!E-3!N '
endif else if max_abs_cntr_lvl le 1.0e-1 then begin
	if cntr_ntv mod 0.0025 eq 0.0 then if cntr_ntv mod 0.0050 ne 0.0 then dcm_plc_nbr=2 else dcm_plc_nbr=1
	p_fmt_fct=1.0e2
	unit_pfx_sng='!9X!5 10!E-2!N '
endif else if max_abs_cntr_lvl le 1.0e+1 then begin
	if cntr_ntv mod 0.25 eq 0.0 then if cntr_ntv mod 0.50 ne 0.0 then dcm_plc_nbr=2 else dcm_plc_nbr=1
endif else if max_abs_cntr_lvl le 1.0e+3 then begin
	dcm_plc_nbr=0
endif else if max_abs_cntr_lvl le 1.0e+4 then begin
	dcm_plc_nbr=1
	p_fmt_fct=1.0e-3
	unit_pfx_sng='!9X!5 10!E3!N '
endif else if max_abs_cntr_lvl le 1.0e+5 then begin
	dcm_plc_nbr=1
	p_fmt_fct=1.0e-4
	unit_pfx_sng='!9X!5 10!E4!N '
endif else if max_abs_cntr_lvl le 1.0e+6 then begin
	dcm_plc_nbr=1
	p_fmt_fct=1.0e-5
	unit_pfx_sng='!9X!5 10!E5!N '
endif else if max_abs_cntr_lvl le 1.0e+7 then begin
	dcm_plc_nbr=1
	p_fmt_fct=1.0e-6
	unit_pfx_sng='!9X!5 10!E6!N '
endif 

cntr_fll_idx=indgen(cntr_lvl_nbr-1)+2
cntr_lvl_lbl=strarr(cntr_lvl_nbr)
cntr_ln_sty=indgen(cntr_lvl_nbr)*0		;linestyle of each level
cntr_thk=replicate(2.0,cntr_lvl_nbr)		;thickness of each level
cntr_which_lbl=indgen(cntr_lvl_nbr)*0		;which levels to label

label_maker_srt: foo=1
cntr_lvl_lbl(0)=auto_sng(cntr_lvl(0)*p_fmt_fct,dcm_plc_nbr,dbg_lvl=dbg_lvl)
for lvl_idx=1,cntr_lvl_nbr-1 do begin
	if cntr_lvl(lvl_idx) lt 0.0 then cntr_ln_sty(lvl_idx)=2 
	cntr_lvl_lbl(lvl_idx)=auto_sng(cntr_lvl(lvl_idx)*p_fmt_fct,dcm_plc_nbr,dbg_lvl=dbg_lvl)
	if cntr_lvl_lbl(lvl_idx) eq cntr_lvl_lbl(lvl_idx-1) then begin
		dcm_plc_nbr=dcm_plc_nbr+1
		print,'cntr_lvl_mk(): level ',lvl_idx,' = ',cntr_lvl_lbl(lvl_idx),' equals level ',lvl_idx-1,' = ',cntr_lvl_lbl(lvl_idx),', increasing dcm_plc_nbr to ',dcm_plc_nbr
		goto, label_maker_srt
	endif
endfor

; Remove superfluous decimal points and zeroes from contour labels
for lvl_idx=0,cntr_lvl_nbr-1 do begin
        cntr_lvl_crr=cntr_lvl_lbl(lvl_idx)
	if dbg_lvl gt 0 then print,'cntr_lbl_mk: lvl = ',cntr_lvl(lvl_idx),', lbl = ',cntr_lvl_crr
	cntr_lvl_crr=strtrim(cntr_lvl_crr)
	; Remove 0's trailing decimal places
	result=strpos(cntr_lvl_crr,'.')
	if result ne -1 then begin
              ; truncate excess zeros after decimal point until no more zeros are left
		if strlen(cntr_lvl_crr) gt 1 then begin
                    truncate_trailing_zeroes_after_decimal_point:
                    if (dbg_lvl gt 2) then print,'cntr_lbl_mk: lvl_idx = ',lvl_idx,', cntr_lvl_crr = ',cntr_lvl_crr,', strlen(cntr_lvl_crr) = ',strlen(cntr_lvl_crr),', strmid(cntr_lvl_crr,strlen(cntr_lvl_crr)-1,1) = ',strmid(cntr_lvl_crr,strlen(cntr_lvl_crr)-1,1)
                    if strmid(cntr_lvl_crr,strlen(cntr_lvl_crr)-1,1) eq '0' then begin
                        cntr_lvl_crr=strmid(cntr_lvl_crr,0,strlen(cntr_lvl_crr)-1)
                        goto,truncate_trailing_zeroes_after_decimal_point
                    endif ; endif truncated a trail zero after decimal point
                endif ; string longer than one
	endif; endif
	result=strpos(cntr_lvl_crr,'0.')
	if result eq 0 and strlen(cntr_lvl_crr) gt 2 then cntr_lvl_crr=strmid(cntr_lvl_crr,1,strlen(cntr_lvl_crr)-1)
	result=strpos(cntr_lvl_crr,'-0')
	if result eq 0 then cntr_lvl_crr='-'+strmid(cntr_lvl_crr,2,strlen(cntr_lvl_crr)-2)
	result=strpos(cntr_lvl_crr,'.')
	if result eq strlen(cntr_lvl_crr)-1 then cntr_lvl_crr=strmid(cntr_lvl_crr,0,strlen(cntr_lvl_crr)-1)
	if cntr_lvl_crr eq '-' then cntr_lvl_crr='0'
	if cntr_lvl_crr eq '.0' then cntr_lvl_crr='0.'
	if dbg_lvl gt 0 then print,'cntr_lbl_mk: lvl = ',cntr_lvl(lvl_idx),', lbl = ',cntr_lvl_crr
        cntr_lvl_lbl(lvl_idx)=cntr_lvl_crr ; Reassign working copy to array
endfor
if dbg_lvl gt 0 then begin
	print,'dcm_plc_nbr = ',dcm_plc_nbr
	print,'p_fmt_fct = ',p_fmt_fct
endif; endif dbg_lvl

end; end cntr_lbl_mk()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End cntr_lbl_mk()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Draw Lat Axes
; Define the axes drawing procedure for contour and image maps.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro lat_axz_drw, $
	ntt=ntt, $ ; [flg] Annotate primary axis with labels
	axz_thk=axz_thk, $
	chr_sz=chr_sz, $ ; [frc] Character size
	chr_thk=chr_thk, $
	dbg_lvl=dbg_lvl, $
	lat_btm=lat_btm, $
	lat_top=lat_top, $
	tck_lng=tck_lng, $
	tck_lcn=tck_lcn, $	
	ttl=ttl, $
	axz_vrt=axz_vrt

if n_elements(ntt) eq 0 then ntt=1 ; [flg] Annotate primary axis with labels
if n_elements(axz_thk) eq 0 then axz_thk=1.0
if n_elements(chr_sz) eq 0 then chr_sz=!p.charsize
if n_elements(chr_thk) eq 0 then chr_thk=chr_sz ; [frc] Character size
if n_elements(dbg_lvl) eq 0 then dbg_lvl=0
if n_elements(tck_lng) eq 0 then tck_lng=0.02 ; default is 0.02, negative means ticks point outwards
if n_elements(tck_lcn) eq 0 then tck_lcn=0
if n_elements(ttl) eq 0 then ttl=''	
if n_elements(axz_vrt) eq 0 then axz_vrt=1

; Color initialization
@ibp_clr.com

if lat_top gt lat_btm then lat_spn=lat_top-lat_btm else lat_spn=360.0-lat_btm+lat_top

if lat_spn ge 120.0 then begin
	lat_lbl_ntv=30.0 ; [nbr] Units between adjacent major tick marks (and labels)
	lat_tck_mnr_nbr=3 ; [nbr] Number of minor tickmarks
endif else $
if lat_spn ge 90.0 then begin
	lat_lbl_ntv=20.0
	lat_tck_mnr_nbr=4 ; [nbr] Number of minor tickmarks
endif else $
if lat_spn ge 60.0 then begin
	lat_lbl_ntv=10.0
	lat_tck_mnr_nbr=5 ; [nbr] Number of minor tickmarks
endif else begin
	lat_lbl_ntv=5.0
	lat_tck_mnr_nbr=5 ; [nbr] Number of minor tickmarks
endelse

lat_lbl_nbr=fix(lat_spn/lat_lbl_ntv)
if lat_top eq -lat_btm and lat_lbl_nbr mod 2 ne 0 then lat_lbl_nbr=lat_lbl_nbr+1

if ntt then tck_fmt_fnc='lat_tck_fmt' else tck_fmt_fnc='null_fmt'

if axz_vrt then begin; draw the latitude labels on a vertical axis

axis, $
	yaxis=0, $
	ystyle=1, $
	ythick=axz_thk, $
	ytick_get=lat_tck, $ ; [flt] Tick mark values (coordinates) at each tick mark
	yticks=lat_lbl_nbr, $ ; [nbr] Number of major intervals, 1 less than # of ticks
	yminor=lat_tck_mnr_nbr, $ ; [nbr] Number of minor tickmarks
	ticklen=tck_lng, $
	ytit=ttl, $
	ytickformat=tck_fmt_fnc, $ ; [fmt,fnc] String format or function() returning string for each tick label
;	yrange=[lat_btm,lat_top], $
	charthick=chr_thk, $
	charsize=chr_sz, $ ; [frc] Character size
	color=clr_blk_idx, $
	/nodata

; NB: using a null string in xtickname (in order to have blank labels)
; does not work because it is over-ridden by some internal labeling 
; routine.

axis, $
	yaxis=1, $
	ystyle=1, $
	ythick=axz_thk, $
	ytick_get=lat_tck, $ ; [flt] Tick mark values (coordinates) at each tick mark
	yticks=lat_lbl_nbr, $ ; [nbr] Number of major intervals, 1 less than # of ticks
	yminor=lat_tck_mnr_nbr, $ ; [nbr] Number of minor tickmarks
	ticklen=tck_lng, $
	ytickformat='null_fmt', $ ; [fmt,fnc] String format or function() returning string for each tick label
;	ytickname=strarr(lat_lbl_nbr), $ ; [sng] Annotation string for each tickmark
;	yrange=[lat_btm,lat_top], $
	charthick=chr_thk, $
	charsize=chr_sz, $ ; [frc] Character size
	color=clr_blk_idx, $
	/nodata

endif else begin; draw the latitude labels on a horizontal axis

axis, $
	xaxis=0, $
	xstyle=1, $
	xthick=axz_thk, $
	xtick_get=lat_tck, $ ; [flt] Tick mark values (coordinates) at each tick mark
	xticks=lat_lbl_nbr, $ ; [nbr] Number of major intervals, 1 less than # of ticks
	xminor=lat_tck_mnr_nbr, $ ; [nbr] Number of minor tickmarks
	ticklen=tck_lng, $
	xtit=ttl, $
	xtickformat=tck_fmt_fnc, $ ; [fmt,fnc] String format or function() returning string for each tick label
	charthick=chr_thk, $
	charsize=chr_sz, $ ; [frc] Character size
	color=clr_blk_idx, $
	/nodata

; NB: using a null string in xtickname (in order to have blank labels)
; does not work because it is over-ridden by some internal labeling 
; routine.

axis, $
	xaxis=1, $
	xstyle=1, $
	xthick=axz_thk, $
	xtick_get=lat_tck, $ ; [flt] Tick mark values (coordinates) at each tick mark
	xticks=lat_lbl_nbr, $ ; [nbr] Number of major intervals, 1 less than # of ticks
	xminor=lat_tck_mnr_nbr, $ ; [nbr] Number of minor tickmarks
	ticklen=tck_lng, $
	xtickformat='null_fmt', $ ; [fmt,fnc] String format or function() returning string for each tick label
;	xtickname=strarr(lat_lbl_nbr), $ ; [sng] Annotation string for each tickmark
	charthick=chr_thk, $
	charsize=chr_sz, $ ; [frc] Character size
	color=clr_blk_idx, $
	/nodata

endelse

if dbg_lvl gt 0 then begin
	print,'!x.crange(0), !x.crange(1) = ',!x.crange(0),!x.crange(1)
	print,'lat_top = ',lat_top
	print,'lat_btm = ',lat_btm
;	print,'tck_lcn = ',tck_lcn
	print,'lat_spn = ',lat_spn
	print,'lat_lbl_ntv = ',lat_lbl_ntv
	print,'lat_min = ',min(lat_tck)
	print,'lat_max = ',max(lat_tck)
	print,'lat_tck_mnr_nbr = ',lat_tck_mnr_nbr
      	print,'lat_lbl_nbr = ',lat_lbl_nbr
      	print,'lat_tck = ',lat_tck
endif; endif dbg

end; end lat_axz_drw()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Draw Lat Axes
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Draw Lon Axes
; Define the axes drawing procedure for contour and image maps.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro lon_axz_drw, $
	ntt=ntt, $ ; [flg] Annotate primary axis with labels
	axz_thk=axz_thk, $
	chr_sz=chr_sz, $ ; [frc] Character size
	chr_thk=chr_thk, $ ; [frc] Character thickness
	dbg_lvl=dbg_lvl, $
	lon_btm=lon_btm, $
	lon_top=lon_top, $
	tck_lcn=tck_lcn, $	
	tck_lng=tck_lng, $
	ttl=ttl, $
	axz_vrt=axz_vrt

if n_elements(ntt) eq 0 then ntt=1 ; [flg] Annotate primary axis with labels
if n_elements(axz_thk) eq 0 then axz_thk=1.0
if n_elements(chr_sz) eq 0 then chr_sz=!p.charsize
if n_elements(chr_thk) eq 0 then chr_thk=chr_sz ; [frc] Character thickness
if n_elements(dbg_lvl) eq 0 then dbg_lvl=0
if n_elements(tck_lng) eq 0 then tck_lng=0.02 ; default is 0.02, negative means ticks point outwards
if n_elements(tck_lcn) eq 0 then tck_lcn=0
if n_elements(ttl) eq 0 then ttl=''	
if n_elements(axz_vrt) eq 0 then axz_vrt=1

if lon_top gt lon_btm then lon_spn=lon_top-lon_btm else lon_spn=360.0-lon_btm+lon_top
if lon_top gt lon_btm and lon_top gt 180.0 and lon_btm lt 180.0 then begin
; When region straddles date line axis is weird
!x.range=[lon_btm+360.0,lon_top+360.0] ; Longitude increases on x-axis
endif; endif resetting axis
if lon_top lt lon_btm then begin
; When region straddles date line axis is weird
!x.range=[lon_btm,lon_top+360.0] ; Longitude increases on x-axis
endif; endif resetting axis

if n_elements(tck_lcn) gt 1 then begin ; user-scaled ticks
	lon_lbl_ntv=tck_lcn(0)
	lon_tck_mnr_nbr=tck_lcn(1) ; [nbr] Number of minor tickmarks
endif else begin; autoscaled ticks
if lon_spn ge 180.0 then begin
	lon_lbl_ntv=60.0
	lon_tck_mnr_nbr=6 ; [nbr] Number of minor tickmarks
endif else $
if lon_spn eq 160.0 then begin
	lon_lbl_ntv=30
	lon_tck_mnr_nbr=6 ; [nbr] Number of minor tickmarks
endif else $
if lon_spn gt 120.0 then begin
	lon_lbl_ntv=45.0
	lon_tck_mnr_nbr=3 ; [nbr] Number of minor tickmarks
endif else $
if lon_spn eq 120.0 then begin
	lon_lbl_ntv=30.0
	lon_tck_mnr_nbr=3 ; [nbr] Number of minor tickmarks
endif else $
if lon_spn ge 90.0 then begin
	lon_lbl_ntv=30.0 ; [nbr] Units between adjacent major tick marks (and labels)
	lon_tck_mnr_nbr=4 ; [nbr] Number of minor tickmarks
endif else $
if lon_spn ge 60.0 then begin
	lon_lbl_ntv=10.0 ; [nbr] Units between adjacent major tick marks (and labels)
	lon_tck_mnr_nbr=5 ; [nbr] Number of minor tickmarks
endif else begin
	lon_lbl_ntv=5.0 ; [nbr] Units between adjacent major tick marks (and labels)
	lon_tck_mnr_nbr=5 ; [nbr] Number of minor tickmarks
endelse;
endelse; autoscaled ticks

lon_lbl_nbr=fix(lon_spn/lon_lbl_ntv)
if lon_top eq -lon_btm and lon_lbl_nbr mod 2 ne 0 then lon_lbl_nbr=lon_lbl_nbr-1

if ntt then tck_fmt_fnc='lon_tck_fmt' else tck_fmt_fnc='null_fmt'

if axz_vrt then begin ; Draw longitude labels on vertical axis

axis, $
	yaxis=0, $
	ystyle=1, $
	ythick=axz_thk, $
	ytick_get=lon_tck, $ ; [flt] Tick mark values (coordinates) at each tick mark
	yticks=lon_lbl_nbr, $ ; [nbr] Number of major intervals, 1 less than # of ticks
	yminor=lon_tck_mnr_nbr, $ ; [nbr] Number of minor tickmarks
	ticklen=tck_lng, $
	ytit=ttl, $
	ytickformat=tck_fmt_fnc, $ ; [fmt,fnc] String format or function() returning string for each tick label
;	yrange=[lon_btm,lon_top], $
	charthick=chr_thk, $ ; [frc] Character thickness
	charsize=chr_sz, $ ; [frc] Character size
	/nodata

; NB: using a null string in xtickname (in order to have blank labels)
; does not work because it is over-ridden by some internal labeling 
; routine.

axis, $
	yaxis=1, $
	ystyle=1, $
	ythick=axz_thk, $
	ytick_get=lon_tck, $ ; [flt] Tick mark values (coordinates) at each tick mark
	yticks=lon_lbl_nbr, $ ; [nbr] Number of major intervals, 1 less than # of ticks
	yminor=lon_tck_mnr_nbr, $ ; [nbr] Number of minor tickmarks
	ticklen=tck_lng, $
	ytickformat=tck_fmt_fnc, $ ; [fmt,fnc] String format or function() returning string for each tick label
;	ytickname=strarr(lon_lbl_nbr), $ ; [sng] Annotation string for each tickmark
;	yrange=[lon_btm,lon_top], $
	charthick=chr_thk, $ ; [frc] Character thickness
	charsize=chr_sz, $ ; [frc] Character size
	/nodata

endif else begin; end vertical, begin the longitude labels on a horizontal axis

axis, $
	xaxis=0, $
	xstyle=1, $
	xthick=axz_thk, $
	xtick_get=lon_tck, $ ; [flt] Tick mark values (coordinates) at each tick mark
	xticks=lon_lbl_nbr, $ ; [nbr] Number of major intervals, 1 less than # of ticks
	xminor=lon_tck_mnr_nbr, $ ; [nbr] Number of minor tickmarks
	ticklen=tck_lng, $
	xtit=ttl, $
	xtickformat=tck_fmt_fnc, $ ; [fmt,fnc] String format or function() returning string for each tick label
	charthick=chr_thk, $ ; [frc] Character thickness
	charsize=chr_sz, $ ; [frc] Character size
	/nodata

; NB: using a null string in xtickname (in order to have blank labels)
; does not work because it is over-ridden by some internal labeling 
; routine.

axis, $
	xaxis=1, $
	xstyle=1, $
	xthick=axz_thk, $
	xtick_get=lon_tck, $ ; [flt] Tick mark values (coordinates) at each tick mark
	xticks=lon_lbl_nbr, $ ; [nbr] Number of major intervals, 1 less than # of ticks
	xminor=lon_tck_mnr_nbr, $ ; [nbr] Number of minor tickmarks
	ticklen=tck_lng, $
	xtickformat='null_fmt', $ ; [fmt,fnc] String format or function() returning string for each tick label
;	xtickname=strarr(lon_lbl_nbr), $ ; [sng] Annotation string for each tickmark
	charthick=chr_thk, $ ; [frc] Character thickness
	charsize=chr_sz, $ ; [frc] Character size
	/nodata

endelse

if dbg_lvl gt 0 then begin
	print,'!x.crange(0), !x.crange(1) = ',!x.crange(0),!x.crange(1)
	print,'lon_top = ',lon_top
	print,'lon_btm = ',lon_btm
	print,'tck_lcn = ',tck_lcn
       	print,'lon_spn = ',lon_spn
	print,'lon_lbl_ntv = ',lon_lbl_ntv
	print,'lon_min = ',min(lon_tck)
	print,'lon_max = ',max(lon_tck)
	print,'lon_tck_mnr_nbr = ',lon_tck_mnr_nbr
      	print,'lon_lbl_nbr = ',lon_lbl_nbr
       	print,'lon_tck = ',lon_tck
endif; endif dbg

end; end lon_axz_drw()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Draw Lon Axes
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Polar Tick Format
; For use in polar coordinate plotting routines
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function plr_tck_fmt,axis,index,value

tickmark='10!E'+strtrim(string(value,format='(f5.1)'),2)+'!N'

xpn=value
;if xpn lt 0 then fmt_sng='(i'+strtrim(abs(floor(xpn))+2,2)+'.'+strtrim(abs(floor(xpn)),2)+')'
if xpn lt 0 then fmt_sng='(i1)'
if xpn ge 0 then fmt_sng='(i1)'

;tck_lbl='10!E'+string(abs(value),format=fmt_sng)+'!N' ; Label is 10^xpn
tck_lbl='!5'+string(abs(value),format=fmt_sng)+'!N' ; Label is xpn

print,'plr_tck_fmt() reports value = ',value,', fmt_sng = ',fmt_sng,', tck_lbl = ',tck_lbl

; It only makes sense to plot radial coordinate on +'ve x-axis
if value lt 0.0 then tck_lbl=''

return,tck_lbl
end; end plr_tck_fmt()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Polar Tick Format
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Positive Tick Format
; For use in polar coordinate plotting routines
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function pst_tck_fmt,axis,index,value

tck_lbl='!5'+auto_sng(abs(value),1)+'!N' ; Label is positive definite

;print,'plr_tck_fmt() reports value = ',value,', fmt_sng = ',fmt_sng,', tck_lbl = ',tck_lbl

; It only makes sense to plot radial coordinate on +'ve x-axis
if value lt 0.0 then tck_lbl=''

return,tck_lbl
end; end pst_tck_fmt()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Positive Tick Format
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Lon Tick Format
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function lon_tck_fmt,axis,index,value

while value gt 180.0 do value=value-360.0
while value lt -180.0 do value=value+360.0
nint_value=round(value)
abs_nint_value=abs(nint_value)

if abs_nint_value ge 100 then fmt_sng='(I3)' else $
if abs_nint_value ge 10  then fmt_sng='(I2)' else $
fmt_sng='(I1)' 

if abs_nint_value eq 0 then return,'!50!E!12_!N!5'
if abs_nint_value eq 180 then return,'!5180!E!12_!N!5'
if nint_value gt 0 then return,'!5'+ $
	string(format=fmt_sng,abs_nint_value)+ $
	'!E!12_!N!5E'
if nint_value lt 0 then return,'!5'+ $
	string(format=fmt_sng,abs_nint_value)+ $
	'!E!12_!N!5W'

return,tck_lbl
end; end lon_tck_fmt()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Lon Tick Format
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Lat Tick Format
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function lat_tck_fmt,axis,index,value

if value gt 180.0 then value=value-360.0
nint_value=round(value)
abs_nint_value=abs(nint_value)

if abs_nint_value ge 100 then fmt_sng='(I3)' else $
if abs_nint_value ge 10  then fmt_sng='(I2)' else $
fmt_sng='(I1)' 

if abs_nint_value eq 0 then return,'!50!E!12_!N!5'
if abs_nint_value eq 180 then return,'!5180!E!12_!N!5'
if nint_value gt 0 then return,'!5'+ $
	string(format=fmt_sng,abs_nint_value)+ $
	'!E!12_!N!5N'
if nint_value lt 0 then return,'!5'+ $
	string(format=fmt_sng,abs_nint_value)+ $
	'!E!12_!N!5S'

return,tck_lbl
end; end lat_tck_fmt()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Lat Tick Format
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Year Tick Format
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function yr_tck_fmt,axis,index,value

value=value mod 12
nint_value=round(value)
abs_nint_value=abs(nint_value)

mth_shrt=['J','F','M','A','M','J','J','A','S','O','N','D']
tck_lbl=mth_shrt(abs_nint_value)
;if abs_nint_value eq 0 then tck_lbl=year

return,tck_lbl
end; end yr_tck_fmt()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Year Tick Format
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Day Tick Format
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function day_tck_fmt,axis,index,value

nint_value=round(value)
abs_nint_value=abs(nint_value)
tck_lbl=auto_sng(abs_nint_value)

return,tck_lbl
end; end day_tck_fmt()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Year Tick Format
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Null Format
; Purpose: Dummy function to return blanks for labels
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function null_fmt,axis,index,value

return,''

end; end null_fmt()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Null Format
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Draw Time Axes
; Define axes drawing procedure for time axes
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro time_axz_drw, $
	ntt=ntt, $ ; [flg] Annotate primary axis with labels
	axz_thk=axz_thk, $
	chr_sz=chr_sz, $ ; [frc] Character size
	chr_thk=chr_thk, $ ; [frc] Character thickness
	dbg_lvl=dbg_lvl, $
	lbl_sty=lbl_sty, $ ; [day, mth_shrt, mth_mdm, mth_lng, none]
	tck_lng=tck_lng, $
	time_ncr=time_ncr, $ ; Separation of tick marks and labels in time units
	time_lbl_ntv=time_lbl_ntv, $ ; [nbr] Units between adjacent major tick marks (and labels)
	time_max=time_max, $ ; Usually, e.g., x.crange(1)
	time_min=time_min, $ ; Usually, e.g., x.crange(0)
	time_unit=time_unit, $ ; Unit of the time axis, 'day', 'mth', 'yr'
	ttl=ttl, $
	axz_vrt=axz_vrt, $
	day_srt=day_srt, $ ; Start day
	mth_srt=mth_srt, $ ; Start month
	yr_srt=yr_srt, $ ; Start year in YY format
	yyyymmdd_srt=yyyymmdd_srt ; Start date in YYYYMMDD format

; NB: [xy]style with which plot has been drawn may need to be passed and used by this axis routine
; NB: it is impossible to write this routine generally with IDL 3.6
; Main difficulties are inability to arbitrarily place the major ticks, plus static limit of tick labels to 30

if n_elements(ntt) eq 0 then ntt=1 ; [flg] Annotate primary axis with labels
if n_elements(axz_thk) eq 0 then axz_thk=1.0
if n_elements(chr_sz) eq 0 then chr_sz=1.5
if n_elements(chr_thk) eq 0 then chr_thk=chr_sz ; [frc] Character thickness
if n_elements(dbg_lvl) eq 0 then dbg_lvl=0
if n_elements(lbl_sty) eq 0 then lbl_sty='none'	; [mth_shrt, mth_mdm, mth_lng, none]
if n_elements(tck_lng) eq 0 then tck_lng=0.02 ; default is 0.02, negative means ticks point outwards
if n_elements(time_ncr) eq 0 then time_ncr=1 ; Separation of tick marks and labels in time units
if n_elements(time_lbl_ntv) eq 0 then time_lbl_ntv=1 ; [nbr] Units between adjacent major tick marks (and labels)
if n_elements(time_max) eq 0 then time_max=0 ; Usually, e.g., x.crange(1)
if n_elements(time_min) eq 0 then time_min=0 ; Usually, e.g., x.crange(0)
if n_elements(time_unit) eq 0 then time_unit='mth' ; unit of the time axis
if n_elements(ttl) eq 0 then ttl=''	
if n_elements(axz_vrt) eq 0 then axz_vrt=1
if n_elements(day_srt) eq 0 then day_srt=1
if n_elements(mth_srt) eq 0 then mth_srt=1
if n_elements(yr_srt) eq 0 then yr_srt=85
if n_elements(x_sty) eq 0 then x_sty=0
if n_elements(y_sty) eq 0 then y_sty=0
if n_elements(yyyymmdd_srt) eq 0 then yyyymmdd_srt=19990312 ; Start date in YYYYMMDD format

time_spn=time_max-time_min+1 ; [] Difference between minimum and maximum abscissa
time_nbr=1+(time_max-time_min)/time_ncr ; [nbr] Number of abscissa units in domain
tck_crd=time_min+indgen(time_nbr)*time_ncr
tck_crd_mth=tck_crd mod 12

; Defaults
time_lbl_nbr=time_nbr ; [nbr] Number of time_ncr's that fit in domain
time_tck_mnr_nbr=1 ; [nbr] Number of minor tickmarks

mth_nmr=indgen(12)+1
mth_shrt=['J','F','M','A','M','J','J','A','S','O','N','D']
mth_mdm=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
mth_lng=['January','February','March','April','May','June','July','August','Sepember','October','November','December']
none=replicate(' ',1)

; In month mode, 0 is January, 11 is December, 12 is January of year 2 ...
; Number of major tickmarks equals number of Januarys

; Label with day number
if time_unit eq 'day' then begin

	x_sty=1
	y_sty=0

	if time_spn ge 180.0 then begin
		time_lbl_ntv=28 ; [nbr] Units between adjacent major tick marks (and labels)
		time_tck_mnr_nbr=4 ; [nbr] Number of minor tickmarks
	endif else if time_spn ge 20.0 then begin
		time_lbl_ntv=7 ; [nbr] Units between adjacent major tick marks (and labels)
		time_tck_mnr_nbr=7 ; [nbr] Number of minor tickmarks
	endif else begin
		time_lbl_ntv=1 ; [nbr] Units between adjacent major tick marks (and labels)
		time_tck_mnr_nbr=1 ; [nbr] Number of minor tickmarks
	endelse ; end else

	time_lbl_nbr=fix(time_spn/time_lbl_ntv)+1 ; [nbr] Number of labels that fit in domain
	yyyymmdd_prs,yyyymmdd_srt,yr=yr_srt,mth=mth_srt,day=day_srt
	tck_lbl=strarr(time_nbr)
;	ttl=mth2sng(mth_srt)+' '+string(format='(I4.4)',yr_srt)
	jdy_srt=julday(mth_srt,day_srt,yr_srt)
	jdy_end=jdy_srt+time_spn
	for idx=0,time_nbr-1 do begin
		tck_lbl(idx)=''
		if idx mod time_lbl_ntv eq 0 then tck_lbl(idx)=string(day_srt+idx)
	endfor; end loop over idx
endif ; endif time_unit eq 'day'

; Label with month abbreviation
if time_unit eq 'mth' then begin
	if time_spn le 13 then begin 
		tck_lbl=strarr(time_nbr)
		for idx=0,time_nbr-1 do begin
			if tck_crd_mth(idx) lt 0 then tck_crd_mth(idx)=tck_crd_mth(idx)+time_nbr-1
			if lbl_sty eq 'mth_shrt' then tck_lbl(idx)=mth_shrt(tck_crd_mth(idx))
			if lbl_sty eq 'mth_mdm' then tck_lbl(idx)=mth_mdm(tck_crd_mth(idx))
			if lbl_sty eq 'mth_lng' then tck_lbl(idx)=mth_lng(tck_crd_mth(idx))
		endfor; end loop over idx
		time_lbl_nbr=time_nbr ; [nbr] Number of labels that fit in domain
		time_tck_mnr_nbr=1 ; [nbr] Number of minor tickmarks
	endif; endif time_spn lt 12
endif; endif time_unit eq 'mth'

; Label with year number
if time_unit eq 'mth' then begin
	if time_spn gt 13 then begin 
		if lbl_sty eq 'mth_shrt' then begin

			jan_nbr=0
			for idx=0,time_nbr-1 do if tck_crd_mth(idx) eq 0 then jan_nbr=jan_nbr+1
			if dbg_lvl gt 0 then print,'tck_crd_mth = ',tck_crd_mth
			if dbg_lvl gt 0 then print,'jan_nbr = ',jan_nbr
			tck_lbl=strarr(jan_nbr)
			for idx=0,jan_nbr-1 do tck_lbl(idx)=auto_sng(yr_srt+idx,0)
		endif; endif lbl_sty eq 'mth_shrt'
		time_lbl_nbr=jan_nbr ; [nbr] Number of labels that fit in domain
		time_tck_mnr_nbr=12 ; [nbr] Number of minor tickmarks
	endif; endif time_spn gt 13
endif; endif time_unit eq 'mth'

if lbl_sty eq 'mth_nmr' then begin
endif else if lbl_sty eq 'mth_shrt' then begin
endif else if lbl_sty eq 'mth_mdm' then begin
endif else if lbl_sty eq 'mth_lng' then begin
endif else if lbl_sty eq 'none' then begin
endif ; endif

; Draw only tickmarks unless axis has label 
if ntt then begin
	if time_unit eq 'mth' then tck_fmt_fnc='yr_tck_fmt'
	if time_unit eq 'day' then tck_fmt_fnc='label_date'
endif else begin
	tck_fmt_fnc='null_fmt'
	tck_lbl(*)=' '
endelse ; endelse

if axz_vrt then begin ; Draw time labels on vertical axis

axis, $
	yaxis=0, $
	ystyle=y_sty, $
	ythick=axz_thk, $
	ytick_get=time_tck, $ ; [flt] Tick mark values (coordinates) at each tick mark
	yticks=time_lbl_nbr-1, $ ; [nbr] Number of major intervals, 1 less than # of ticks
	ytickname=tck_lbl, $ ; [sng] Annotation string for each tickmark
	ytit=ttl, $
	ytickformat=null_fmt, $ ; [fmt,fnc] String format or function() returning string for each tick label
	yminor=time_tck_mnr_nbr, $ ; [nbr] Number of minor tickmarks
	ticklen=tck_lng, $
	charthick=chr_thk, $ ; [frc] Character thickness
	charsize=chr_sz, $ ; [frc] Character size
	/nodata

; NB: using a null string in xtickname (in order to have blank labels)
; does not work because it is over-ridden by some internal labeling routine

axis, $
	yaxis=1, $
	ystyle=y_sty, $
	ythick=axz_thk, $
	ytick_get=time_tck, $ ; [flt] Tick mark values (coordinates) at each tick mark
	yticks=time_lbl_nbr-1, $ ; [nbr] Number of major intervals, 1 less than # of ticks
	yminor=time_tck_mnr_nbr, $ ; [nbr] Number of minor tickmarks
	ticklen=tck_lng, $
	ytickformat=tck_fmt_fnc, $ ; [fmt,fnc] String format or function() returning string for each tick label
	charthick=chr_thk, $ ; [frc] Character thickness
	charsize=chr_sz, $ ; [frc] Character size
	/nodata

endif else begin ; Draw time labels on horizontal axis

axis, $
	xaxis=0, $
	xstyle=x_sty, $
	xthick=axz_thk, $
	xtick_get=time_tck, $ ; [flt] Tick mark values (coordinates) at each tick mark
	xticks=time_lbl_nbr-1, $ ; [nbr] Number of major intervals, 1 less than # of ticks
	xminor=time_tck_mnr_nbr, $ ; [nbr] Number of minor tickmarks
;	xtickname=tck_lbl, $ ; [sng] Annotation string for each tickmark
	xtit=ttl, $
	xtickformat=tck_fmt_fnc, $ ; [fmt,fnc] String format or function() returning string for each tick label
	ticklen=tck_lng, $
	charthick=chr_thk, $ ; [frc] Character thickness
	charsize=chr_sz, $ ; [frc] Character size
	/nodata

; NB: using a null string in xtickname (in order to have blank labels)
; does not work because it is over-ridden by some asinine internal labeling routine.

axis, $
	xaxis=1, $
	xstyle=x_sty, $
	xthick=axz_thk, $
	xtick_get=time_tck, $ ; [flt] Tick mark values (coordinates) at each tick mark
	xticks=time_lbl_nbr-1, $ ; [nbr] Number of major intervals, 1 less than # of ticks
	xminor=time_tck_mnr_nbr, $ ; [nbr] Number of minor tickmarks
	ticklen=tck_lng, $
	xtickformat='null_fmt', $ ; [fmt,fnc] String format or function() returning string for each tick label
	charthick=chr_thk, $ ; [frc] Character thickness
	charsize=chr_sz, $ ; [frc] Character size
	/nodata

endelse

if dbg_lvl gt 0 then begin
	print,'!x.crange(0) = ',!x.crange(0)
	print,'!x.crange(1) = ',!x.crange(1)
	print,'time_unit = ',time_unit,' = Unit of the time axis, day or mth'
	print,'ntt = ',ntt ; [flg] Annotate primary axis with labels
	print,'x_sty = ',x_sty
	print,'day_srt = ',day_srt
	print,'mth_srt = ',mth_srt
	print,'yr_srt = ',yr_srt
	print,'yyyymmdd_srt = ',yyyymmdd_srt
	print,'time_min = ',time_min,' =  Usually, e.g., x.crange(0)'
	print,'time_max = ',time_max,' =  Usually, e.g., x.crange(1)'
	print,'time_ncr = ',time_ncr
	print,'time_spn = ',time_spn
	print,'time_lbl_ntv = ',time_lbl_ntv,' = [nbr] Units between adjacent major tick marks (and labels)'
	print,'time_tck_mnr_nbr = ',time_tck_mnr_nbr,' = [nbr] Number of minor tickmarks'
	print,'time_lbl_nbr = ',time_lbl_nbr
	print,'time_tck_min = ',min(time_tck)
	print,'time_tck_max = ',max(time_tck)
	print,'tck_fmt_fnc = ',tck_fmt_fnc
	print,'time_tck = ',time_tck
	print,'tck_lbl = ',tck_lbl
endif; endif dbg

end; end time_axz_drw()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Draw Time Axes
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Draw Generic Axes
; Define axes drawing procedure for generic axes
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro gnr_axz_drw, $
	abc=abc, $ ; [frc] Abscissas
	ord=ord, $ ; [frc] Ordinates
	log_flg=log_flg, $ ; [flg] Axis is logarithmic
	ntt=ntt, $ ; [flg] Annotate primary axis with labels
	axz_thk=axz_thk, $
	chr_sz=chr_sz, $ ; [frc] Character size
	chr_thk=chr_thk, $ ; [frc] Character thickness
	dbg_lvl=dbg_lvl, $
	lbl_sty=lbl_sty, $
	tck_lng=tck_lng, $
	axz_clr=axz_clr, $ ; [idx] Axis color
	axz_max=axz_max, $
	axz_min=axz_min, $
	axz_unit=axz_unit, $
	ttl=ttl, $
	axz_vrt=axz_vrt

; NB: [xy]style with which plot has been drawn may need to be passed and used by this axis routine
; NB: it is impossible to write this routine generally with IDL 3.6
; Main difficulties are inability to arbitrarily place the major ticks, plus static limit of tick labels to 30

if n_elements(abc) eq 0 then abc=1 ; [frc] Abscissas
if n_elements(ord) eq 0 then ord=1 ; [frc] Ordinates 
if n_elements(ntt) eq 0 then ntt=1 ; [flg] Annotate primary axis with labels axes labels
if n_elements(axz_thk) eq 0 then axz_thk=1.0
if n_elements(chr_sz) eq 0 then chr_sz=!p.charsize
if n_elements(chr_thk) eq 0 then chr_thk=1.0 ; [frc] Character thickness
if n_elements(dbg_lvl) eq 0 then dbg_lvl=0
if n_elements(lbl_sty) eq 0 then lbl_sty='none'	
if n_elements(tck_lng) eq 0 then tck_lng=0.02 ; default is 0.02, negative means ticks point outwards
if n_elements(axz_ncr) eq 0 then axz_ncr=1 ; Separation of tick marks and labels in time units
if n_elements(axz_lbl_ntv) eq 0 then axz_lbl_ntv=1 ; [nbr] Units between adjacent major tick marks (and labels)
if n_elements(axz_unit) eq 0 then axz_unit='mth' ; unit of the time axis
if n_elements(ttl) eq 0 then ttl=''	
if n_elements(axz_vrt) eq 0 then axz_vrt=1
if n_elements(x_sty) eq 0 then x_sty=0
if n_elements(y_sty) eq 0 then y_sty=0
if n_elements(axz_clr) eq 0 then axz_clr=0 ; [idx] Axis color

if n_elements(log_flg) eq 0 then begin
	if axz_vrt eq 0 then log_flg=!x.type else log_flg=!y.type
endif ; endif log_flg
axz_typ=log_flg

if n_elements(axz_max) eq 0 then begin
	if axz_vrt eq 0 then axz_max=!x.crange(1) else axz_max=!y.crange(1)
	if axz_typ eq 1 then axz_max=10.0^axz_max
endif ; endif axz_max
if n_elements(axz_min) eq 0 then begin
	if axz_vrt eq 0 then axz_min=!x.crange(0) else axz_min=!y.crange(0)
	if axz_typ eq 1 then axz_min=10.0^axz_min
endif ; endif axz_min

axz_spn=axz_max-axz_min

; For some lame reason the following is required to set logarithmic axis coordinates
if axz_vrt then begin
	!y.range=[axz_min,axz_max]
endif else begin
	!x.range=[axz_min,axz_max]
endelse; endif axz_vrt

; Draw only tickmarks unless axis has label 
if ntt then begin
	tck_fmt_fnc='null_fmt'
endif else begin
	tck_fmt_fnc='null_fmt'
endelse ; endelse

if axz_vrt then begin ; Draw vertical axis

axis, $
	yaxis=0, $
	ylog=axz_typ, $ ; [flg] Axis is logarithmic
	ystyle=y_sty, $
	ythick=axz_thk, $
	ytick_get=axz_tck, $ ; [flt] Tick mark values (coordinates) at each tick mark
;	yticks=axz_lbl_nbr-1, $ ; [nbr] Number of major intervals, 1 less than # of ticks
;	ytickname=tck_lbl, $ ; [sng] Annotation string for each tickmark
	ytit=ttl, $
;	ytickformat=tck_fmt_fnc, $ ; [fmt,fnc] String format or function() returning string for each tick label
;	yminor=axz_tck_mnr_nbr, $ ; [nbr] Number of minor tickmarks
	ticklen=tck_lng, $
	charthick=chr_thk, $ ; [frc] Character thickness
	charsize=chr_sz, $ ; [frc] Character size
	color=axz_clr, $ ; [idx] Axis color
	/nodata

; NB: using null string in xtickname (in order to have blank labels)
; does not work because it is over-ridden by some internal labeling routine.

axis, $
	yaxis=1, $
	ylog=axz_typ, $ ; [flg] Axis is logarithmic
	ystyle=y_sty, $
	ythick=axz_thk, $
	ytick_get=axz_tck, $ ; [flt] Tick mark values (coordinates) at each tick mark
;	yticks=axz_lbl_nbr-1, $ ; [nbr] Number of major intervals, 1 less than # of ticks
;	yminor=axz_tck_mnr_nbr, $ ; [nbr] Number of minor tickmarks
	ticklen=tck_lng, $
	ytickformat='null_fmt', $ ; [fmt,fnc] String format or function() returning string for each tick label
	charthick=chr_thk, $ ; [frc] Character thickness
	charsize=chr_sz, $ ; [frc] Character size
	color=axz_clr, $ ; [idx] Axis color
	/nodata

endif else begin ; Draw horizontal axis

if ntt then begin

axis, $
	xaxis=0, $
	xlog=axz_typ, $ ; [flg] Axis is logarithmic
	xstyle=x_sty, $
	xthick=axz_thk, $
	xtick_get=axz_tck, $ ; [flt] Tick mark values (coordinates) at each tick mark
;	xticks=axz_lbl_nbr-1, $ ; [nbr] Number of major intervals, 1 less than # of ticks
;	xminor=axz_tck_mnr_nbr, $ ; [nbr] Number of minor tickmarks
;	xtickname=tck_lbl, $ ; [sng] Annotation string for each tickmark
	xtit=ttl, $
;	xtickformat=tck_fmt_fnc, $ ; [fmt,fnc] String format or function() returning string for each tick label
	ticklen=tck_lng, $
	charthick=chr_thk, $ ; [frc] Character thickness
	charsize=chr_sz, $ ; [frc] Character size
	color=axz_clr, $ ; [idx] Axis color
	/nodata

endif else begin

axis, $
	xaxis=0, $
	xlog=axz_typ, $ ; [flg] Axis is logarithmic
	xstyle=x_sty, $
	xthick=axz_thk, $
	xtick_get=axz_tck, $ ; [flt] Tick mark values (coordinates) at each tick mark
;	xticks=axz_lbl_nbr-1, $ ; [nbr] Number of major intervals, 1 less than # of ticks
;	xminor=axz_tck_mnr_nbr, $ ; [nbr] Number of minor tickmarks
;	xtickname=tck_lbl, $ ; [sng] Annotation string for each tickmark
	xtit=ttl, $
	xtickformat=tck_fmt_fnc, $ ; [fmt,fnc] String format or function() returning string for each tick label
	ticklen=tck_lng, $
	charthick=chr_thk, $ ; [frc] Character thickness
	charsize=chr_sz, $ ; [frc] Character size
	color=axz_clr, $ ; [idx] Axis color
	/nodata

endelse; endelse ntt

; NB: using a null string in xtickname (in order to have blank labels)
; does not work because it is over-ridden by some internal labeling routine.

axis, $
	xaxis=1, $
	xlog=axz_typ, $ ; [flg] Axis is logarithmic
	xstyle=x_sty, $
	xthick=axz_thk, $
	xtick_get=axz_tck, $ ; [flt] Tick mark values (coordinates) at each tick mark
;	xticks=axz_lbl_nbr-1, $ ; [nbr] Number of major intervals, 1 less than # of ticks
;	xminor=axz_tck_mnr_nbr, $ ; [nbr] Number of minor tickmarks
	ticklen=tck_lng, $
	xtickformat='null_fmt', $ ; [fmt,fnc] String format or function() returning string for each tick label
	charthick=chr_thk, $ ; [frc] Character thickness
	charsize=chr_sz, $ ; [frc] Character size
	color=axz_clr, $ ; [idx] Axis color
	/nodata

endelse

if dbg_lvl gt 0 then begin
	print,'!x.crange(0), !x.crange(1) = ',!x.crange(0),!x.crange(1)
	print,'axz_unit = ',axz_unit
	print,'ntt = ',ntt ; [flg] Annotate primary axis with labels
	print,'x_sty = ',x_sty
	print,'axz_typ = ',axz_typ
	print,'axz_min = ',axz_min
	print,'axz_max = ',axz_max
	print,'axz_ncr = ',axz_ncr
	print,'axz_spn = ',axz_spn
;	print,'axz_lbl_ntv = ',axz_lbl_ntv
;	print,'axz_tck_mnr_nbr = ',axz_tck_mnr_nbr
;	print,'axz_lbl_nbr = ',axz_lbl_nbr
	print,'axz_tck_min = ',min(axz_tck)
	print,'axz_tck_max = ',max(axz_tck)
	print,'tck_fmt_fnc = ',tck_fmt_fnc
	print,'axz_tck = ',axz_tck
;	print,'tck_lbl = ',tck_lbl
endif; endif dbg

end; end gnr_axz_drw()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Draw Generic Axes
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
