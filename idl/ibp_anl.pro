;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; RCS Identification
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; $Author: zender $
; $Date$
; $Id$
; $Revision$
; $Locker:  $
; $RCSfile: ibp_anl.pro,v $
; $Source: /home/zender/cvs/idl/ibp_anl.pro,v $
; $Id$
; $State: Exp $
;
; NB: get RCS formatting in IDL files by using rcs -U -c"; " foo.pro
;
; $Log: not supported by cvs2svn $
; Revision 1.6  2000-03-01 02:28:33  zender
; *** empty log message ***
;
; Revision 1.5  2000/01/10 23:36:31  zender
; *** empty log message ***
;
; Revision 1.3  1999/12/31 20:12:47  zender
; *** empty log message ***
;
; Revision 1.2  1999/12/31 02:09:38  zender
; *** empty log message ***
;
; Revision 1.1  1999/12/31 00:18:20  zender
; *** empty log message ***
;
; Revision 1.1.1.1  1999/05/11 22:20:47  zender
; Imported sources
;
; Revision 3.8  1998/10/19 19:01:36  zender
; maintenance check-in
;
; Revision 3.7  1995/09/22  17:43:49  zender
; removed lon_lat map from main screen, implemented adjustable
; window sizes. still have bug in color table 4.
;
; Revision 3.6  1995/06/15  17:31:35  zender
; added azimuthal polar plots, changed distribution slightly.
;
; Revision 3.5  1995/05/09  22:55:21  zender
; finished merging doetzl's mods. added a customize screen to
; each figure type. removed buttons from main window. added
; user-specified contour levels. merged contour level database
; into pre_fld structure.
;
; Revision 3.4  1995/05/08  17:14:46  zender
; finished merge with doetzl's customization mods. about to
; bring customization to all cases.
;
; Revision 3.3  1995/05/07  22:59:09  zender
; about to merge in doetzl's mods.
;
; Revision 3.2  1995/04/07  21:38:54  zender
; first version for general release. renamed everything ibp
; instead of ccm.
;
; Revision 3.1  1995/04/07  21:17:20  zender
; final version before ccm -> ibp name switch. diddled with
; color routines to add a black/white palette. switched off
; user input of cntr_lvls in ibp_lat_lev but now allow exact
; contour specs. in ibp_fgr for pre-defined fields.
;
; Revision 3.0  1995/03/27  20:30:54  zender
; added animation widgets (still unstable), fixed longitude shifting
; bug, beta-tested with BEE's and JET's netCDF files OK.
;
; Revision 1.2  1995/02/09  22:02:00  zender
; synchronization checkin. ibp_lat_lev is in and working, although it's not
; very generic. don't know what to do next....
;
; Revision 1.1  1995/02/08  23:14:42  zender
; Initial revision
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin USER HOOK analysis
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro usr_hook
;
end_of_procedure: foo=1
end; end usr_hook()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End USER HOOK analysis
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Differencing
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro difference
@ibp_cmn.com
;
fld_nm='d'+fld_nm
orig_data=data
bad_data_idx=where(data gt very_big_nbr,count)
bad_sv_data_idx=where(sv_data gt very_big_nbr,sv_count)
data=orig_data-sv_data
;
if count ne 0 then data(bad_data_idx)=orig_data(bad_data_idx)
if sv_count ne 0 then data(bad_sv_data_idx)=sv_data(bad_sv_data_idx)
;
diagnose_data
form_cntr_lvl
comp_cntr_lvl
;
; do a manipulation to ease the writing of the GIF files
;
src_nm='388_obs'
;
src_sng=src_sng+' - '+sv_src_sng
;fld_sng=fld_sng+' - '+sv_fld_sng
abb_sng=abb_sng+' - '+sv_abb_sng
ttl_sng=slice_sng+' '+fld_sng+', '+month_sng
sub_ttl_sng=src_sng
;
end_of_procedure: foo=1
end; end difference()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Differencing
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin High/Low Labels
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro high_low_lbls
@ibp_cmn.com
;
;lat_idx=[0]
;lon_idx=[0]
;max_count=0
;;
;for lat_idx=1,lat_nbr-2 do begin
;for lon_idx=1,lon_nbr-2 do begin
;	if data(lon_idx,lat_idx) lt very_big_nbr then $
;	if data(lon_idx,lat_idx) gt data(lon_idx-1,lat_idx-1) then $
;	if data(lon_idx,lat_idx) gt data(lon_idx-1,lat_idx) then $
;	if data(lon_idx,lat_idx) gt data(lon_idx-1,lat_idx+1) then $
;	if data(lon_idx,lat_idx) gt data(lon_idx,lat_idx-1) then $
;	if data(lon_idx,lat_idx) gt data(lon_idx,lat_idx+1) then $
;	if data(lon_idx,lat_idx) gt data(lon_idx+1,lat_idx-1) then $
;	if data(lon_idx,lat_idx) gt data(lon_idx+1,lat_idx) then $
;	if data(lon_idx,lat_idx) gt data(lon_idx+1,lat_idx+1) then begin
;		lat_idx=[lat_idx,lat_idx]		
;		lon_idx=[lon_idx,lon_idx]		
;		max_count=max_count+1
;	endif
;endfor
;endfor
;;
;print,'max_count = ',max_count
;lat_idx=lat_idx(1:max_count)
;lon_idx=lon_idx(1:max_count)
;;
;if max_count ne 0 then begin
;hi_sngs=strarr(max_count)
;for idx=0,max_count-1 do begin
;	hi_sngs(idx)=auto_sng(data(lon_idx(idx),lat_idx(idx)),0)
;	print,hi_sngs(idx)
;endfor
;wset,lon_lat_wdw
;outline_points,max_count,lon_idx,lat_idx
;endif; endif max_count ne 0
;
; The technique used here is taken from image processing. 
; In order to get the koontz() routine to work quickly and reliably,
; the data array is first converted to an integer array by bytscl.
; The returned indices of hi's and lo's are then reconstructed
; to latitude and longitude indices and the labels created and
; plotted.
;
byte_data=bytscl( $
		data, $
		max=very_big_nbr, $
;		min=data_min, $
		top=lat_nbr*lon_nbr)	
;
koontz_data=koontz(byte_data)
max_idx=where(koontz_data ge 1.1,max_count)
print,'max_count = ',max_count
if max_count ne 0 then begin
;
; max_idx is an array of longword indices into the original data.
; We can now invert it to find the latitude and longitude of the 
; local maxima.
;
lat_idx=indgen(max_count)
lon_idx=indgen(max_count)
hi_sngs=strarr(max_count)
;
for idx=0,max_count-1 do begin
	lon_idx(idx)=max_idx(idx) mod lon_nbr
	lat_idx(idx)=fix(max_idx(idx)/lon_nbr)
	hi_sngs(idx)=auto_sng(data(max_idx(idx)),0)
	print,hi_sngs(idx)
;
	print,'lon = '+auto_sng(lon(lon_idx(idx)),0)+ $
	', lat = '+auto_sng(lat(lat_idx(idx)),0)+ $
	', '+fld_nm+' = '+ $
	auto_sng(data(lon_idx(idx),lat_idx(idx)),0)
endfor
;
wset,lon_lat_wdw
outline_points,max_count,lon_idx,lat_idx
endif; endif max_count ne 0
;
end_of_procedure: foo=1
end; end high_low_lbls()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End High/Low Labels
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin PRECC CMFMC analysis
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro precc_cmfmc
@ibp_cmn.com
;
good_data=where(data gt 15.,good_data_nbr)
if good_data_nbr eq 0 then begin
	print,'No matches in current analysis, exiting precc_cmfmc()...'
	goto, end_of_procedure
endif
precc=data(good_data)
cmfmc=sv_data(good_data)
oro=msk_data(good_data)
;
; find the oro field and mask out those events not over ocean
;
;fld_idx=0
;while (oro_idx eq 0 and fld_idx lt fld_nbr) do begin
;	if 'ORO' eq fld_lst(fld_idx).name then begin
;		oro_idx=fld_idx
;	endif	
;	fld_idx=fld_idx+1
;	nc_id=ncdf_open(fl_in)
;	ncdf_varget,nc_id,'ORO',oro
;endwhile
;
ocean_data=where(oro eq 0.0,count)
if count ne 0 then begin
	precc=precc(ocean_data)
	cmfmc=cmfmc(ocean_data)
	oro=oro(ocean_data)
endif
;
;print,'time = ',time(time_slc)
;print,'precc = ',precc
;print,'cmfmc = ',cmfmc
;print,'oro = ',oro
;
pth_nm_get,fl_out,path,name
analyze_fl=path+'/analyze.nc'
matching_fls=findfile(analyze_fl)
print,'matching_fls = ',matching_fls
count=n_elements(matching_fls)
if (size(matching_fls))(0) eq 0 then begin
	unlim_dim_sz=0
	time_dim_sz=0
	lev_dim_sz=10
	nc_id=ncdf_create(analyze_fl,/noclobber)
	unlim_dim_id=ncdf_dimdef(nc_id,'unlim',/unlimited)
	time_dim_id=ncdf_dimdef(nc_id,'time')
	lev_dim_id=ncdf_dimdef(nc_id,'lev',lev_dim_sz)
	lev_id=ncdf_vardef(nc_id,'lev',[lev_dim_id],/float)
	cmfmc_id=ncdf_vardef(nc_id,'CMFMC',[lev_dim_id,unlim_dim_id],/float)
	time_id=ncdf_vardef(nc_id,'time',[time_dim_id],/float)
;	cmfmc_id=ncdf_vardef(nc_id,'CMFMC',[unlim_dim_id],/float)
	precc_id=ncdf_vardef(nc_id,'PRECC',[unlim_dim_id],/float)
	ncdf_control,nc_id,/endef
	ncdf_varput,nc_id,lev_id,lev
endif else begin
	nc_id=ncdf_open(analyze_fl,/write)
	foo=ncdf_inquire(nc_id)
	time_id=ncdf_varid(nc_id,'time')
	cmfmc_id=ncdf_varid(nc_id,'CMFMC')
	precc_id=ncdf_varid(nc_id,'PRECC')
	lev_dim_id=ncdf_dimid(nc_id,'lev')
	unlim_dim_id=ncdf_dimid(nc_id,'unlim')
	time_dim_id=ncdf_dimid(nc_id,'time')
	ncdf_diminq,nc_id,unlim_dim_id,foo,unlim_dim_sz
	ncdf_diminq,nc_id,time_dim_id,foo,time_dim_sz
	ncdf_diminq,nc_id,lev_dim_id,foo,lev_dim_sz
endelse
;
ncdf_varput,nc_id,time_id,time(time_slc),offset=[time_dim_sz]
ncdf_varput,nc_id,cmfmc_id,cmfmc,offset=[0,unlim_dim_sz]
ncdf_varput,nc_id,precc_id,precc,offset=[unlim_dim_sz]
ncdf_close,nc_id
;
end_of_procedure: foo=1
end; end precc_cmfmc()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End PRECC CMFMC analysis
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



