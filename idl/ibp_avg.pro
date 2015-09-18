;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Area Average
; Purpose: This routine performs an area-weighted average on
; a region of the sphere in grid space.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function area_avg, $
	data, $
	lat, $
	lon, $
	lat_nbr, $
	lon_nbr, $
	lat_max_idx, $	; index into lat array
	lat_min_idx, $	; index into lat array
	lat_wgt=lat_wgt, $
	msk_data=msk_data, $
	msk_mode=msk_mode, $
	msk_val=msk_val, $
	vld_avg_thr_pct=vld_avg_thr_pct, $
	very_big_nbr=very_big_nbr

if n_elements(very_big_nbr) eq 0 then very_big_nbr=1.0e20
if n_elements(vld_avg_thr_pct) eq 0 then vld_avg_thr_pct=66.6
if n_elements(msk_mode) eq 0 then msk_mode=0 ; 0 is False
if n_elements(msk_data) eq 0 then msk_data=0 ;
if n_elements(msk_val) eq 0 then msk_val=0 ; 0 is CCM Ocean
if n_elements(foo) eq 0 then foo=0
if n_elements(lat_wgt) eq 0 then lat_wgt=(2.0/lat_nbr)+0.0*lat ; Default is all weights are equal and sum to 2.0

if (lat_max_idx-lat_min_idx+1) ne lat_nbr then begin 
	print,'mismatched arrays in area_avg()'
	goto,end_of_procedure
endif

; This nifty matrix trick fills each of the lat_nbr rows of the result
; with the corresponding element of lat_wgt. There are lon_nbr
; elements per row in the result.
lat_wgt_mtx=replicate(1.0,lon_nbr)#lat_wgt(lat_min_idx:lat_max_idx)

if not msk_mode then begin
	good_data=data
	good_lat_wgt_mtx=lat_wgt_mtx
endif else begin ; endif not msk_mode
	if n_elements(msk_data) ne n_elements(data) then begin
		print,'Mask Region is not the same size as data, cannot perform average...'
		average=0.		
		goto,end_of_procedure	
	endif

	good_idx=where(msk_data eq msk_val,count)
	if count ne 0 then begin 
		good_data=data(good_idx)
		good_lat_wgt_mtx=lat_wgt_mtx(good_idx)
	endif else begin
		print,'Region is completely masked, cannot perform average...'
		average=0.		
		goto,end_of_procedure	
	endelse
endelse; endif mask mode

good_idx=where(good_data lt very_big_nbr,count)
vld_cnt_pct=100.0*count/(lon_nbr*lat_nbr)

if vld_cnt_pct lt vld_avg_thr_pct then begin
	print,'WARNING: Only '+auto_sng(vld_cnt_pct,1)+'% of good data in this regional average.'
endif

if count ne 0 then begin
        good_data=good_data(good_idx) 
        good_lat_wgt_mtx=good_lat_wgt_mtx(good_idx) 
        average=total(good_data*good_lat_wgt_mtx)/total(good_lat_wgt_mtx)
endif else begin
	print,'Region that isn''t completely masked is completely blocked data, cannot perform average...'
	average=0.		
	goto,end_of_procedure	
endelse

end_of_procedure: foo=1
return,average

end; end area_avg()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Area Average
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
