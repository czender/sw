; $Id$

; Purpose: offline analysis routines for looking at CCM data. 
; These routines rely heavily on netCDF.

; $Log: not supported by cvs2svn $
; Revision 1.10  2000-10-10 22:52:15  zender
; *** empty log message ***
;
; Revision 1.9  2000/01/15 02:07:49  zender
; *** empty log message ***
;
; Revision 1.7  2000/01/01 01:55:48  zender
; *** empty log message ***
;
; Revision 1.5  1999/12/31 02:09:36  zender
; *** empty log message ***
;
; Revision 1.4  1999/12/31 00:18:17  zender
; *** empty log message ***
;
; Revision 1.2  1999/10/12 07:37:41  zender
; *** empty log message ***
;
; Revision 1.1.1.1  1999/05/11 22:20:48  zender
; Imported sources
;
; Revision 1.5  1995/12/28  00:13:58  zender
; added _bch() routines
;
; Revision 1.4  1995/11/10  23:55:38  zender
; recoded calls to refer to /data/zender/_aux0_/...
;
; Revision 1.3  1994/09/21  02:13:34  zender
; fixed problem with large datasets. implemented generic histogram
; routines. spruced up precip routines. added area averaging.
;
; Revision 1.2  1994/08/18  19:57:14  zender
; this version works with smallish netcdf files but runs out of memory
; while allocating array space for my daily instantaneous (72 step)
; files of cmfmc and precc. will try only reading in neccessary data,
; i.e., reading in the region-truncated data not the whole kielbasa.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Analysis Batch
; This batch routine calls the selected analysis procedure with all
; the available data.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro analysis_bch, $
	precc_thr=precc_thr, $
	dom_avg=dom_avg, $
	foo=foo
;
; analysis_bch,dom_avg=1,precc_thr=15.
; analysis_bch,dom_avg=1,precc_thr=0.
;
if n_elements(dom_avg) eq 0 then dom_avg=0
if n_elements(foo) eq 0 then foo=0
if n_elements(precc_thr) eq 0 then precc_thr=15.
;
if DOM_AVG then begin
	spawn,'rm -f /data/zender/_aux0_/cem/dom_avg.nc'
endif else begin
	spawn,'rm -f /data/zender/_aux0_/cem/grid_point.nc'
endelse
;
for day=1,14 do begin
;
fl_in='omega.'+string(format='(I2.2)',day)+'.nc'
fl_in='/data/zender/_aux0_/omega0_2/'+fl_in
;
if DOM_AVG then begin
;
cmfmc_dom_avg, $
	lev_min_idx=1,lev_max_idx=7, $
	time_min_idx=0,time_max_idx=71, $
	precc_thr=precc_thr, $
	fl_in=fl_in, $
	fl_out='/data/zender/_aux0_/cem/dom_avg.nc'
;
endif else begin
;
cmfmc_point, $
	lev_min_idx=1,lev_max_idx=7, $
	time_min_idx=0,time_max_idx=71, $
	precc_thr=precc_thr, $
	fl_in=fl_in, $
	fl_out='/data/zender/_aux0_/cem/grid_point.nc'
;
endelse
;
endfor; end loop over days
;
end; end analysis_bch()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Analysis Batch
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Graphics Batch
; This batch routine calls the selected graphics procedure with all
; the available data.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro graphics_bch, $
	dom_avg=dom_avg, $
	foo=foo
;
if n_elements(dom_avg) eq 0 then dom_avg=0
if n_elements(foo) eq 0 then foo=0
;
open_ps,fl_nm='/data/zender/ps/ccm_mc_hst_500mb.eps',x_sz=5.1,y_sz=3.0,/eps
cmfmc_graphs,level=2,cmfmc_bin_sz=2.5,fl_in='/data/zender/_aux0_/cem/dom_avg.nc',graph_prect=graph_prect,ttl_sng='GCM',rng_y=[0,30]
close_ps,fl_nm='/data/zender/ps/ccm_mc_hst_500mb.eps'
;
open_ps,fl_nm='/data/zender/ps/ccm_mc_hst_200mb.eps',x_sz=5.1,y_sz=3.0,/eps
cmfmc_graphs,level=6,cmfmc_bin_sz=2.5,fl_in='/data/zender/_aux0_/cem/dom_avg.nc',graph_prect=graph_prect,ttl_sng='GCM',rng_y=[0,45]
close_ps,fl_nm='/data/zender/ps/ccm_mc_hst_200mb.eps'
;
;for level=0,6 do begin
;if DOM_AVG then begin
;	if level eq 6 then graph_prect=1 else graph_prect=0
;	cmfmc_graphs,level=level,cmfmc_bin_sz=2.5,fl_in='/data/zender/_aux0_/cem/dom_avg.nc',graph_prect=graph_prect
;endif else begin
;	cmfmc_graphs,level=level,fl_in='/data/zender/_aux0_/cem/grid_point.nc'
;endelse
;
;endfor; end loop over levels
;
end; end graphics_bch()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Graphics Batch
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Domain Average analysis
; In this analysis procedure, values meeting the specified masks and
; other selection criteria are domain averaged into one value. Time
; is the unlimited dimension, and level is the known dimension. 
; The philosophy behind this procedure is that the region has already
; been selected by some time-dependent criteria, i.e., its monthly 
; average convective rainfall.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro cmfmc_dom_avg, $
	time_min_idx=time_min_idx, $
	time_max_idx=time_max_idx, $
	lev_min_idx=lev_min_idx, $
	lev_max_idx=lev_max_idx, $
	precc_thr=precc_thr, $
	fl_in=fl_in, $
	fl_out=fl_out, $
	foo=foo
;
; cmfmc_dom_avg,lev_min_idx=7,lev_max_idx=9
; cmfmc_dom_avg,time_min_idx=0,time_max_idx=2
; cmfmc_dom_avg,lev_min_idx=3,lev_max_idx=7,fl_in='/data/zender/_aux0_/omega0_2/omega.01.nc'
; cmfmc_dom_avg,lev_min_idx=3,lev_max_idx=7,time_min_idx=0,time_max_idx=71,fl_in='/data/zender/_aux0_/omega0_2/omega.01.nc'
;
if n_elements(time_min_idx) eq 0 then time_min_idx=0
if n_elements(time_max_idx) eq 0 then time_max_idx=71
if n_elements(lev_min_idx) eq 0 then lev_min_idx=0
if n_elements(lev_max_idx) eq 0 then lev_max_idx=9
if n_elements(lat_min) eq 0 then lat_min=10.
if n_elements(lon_min) eq 0 then lon_min=120.
if n_elements(lat_max) eq 0 then lat_max=20.
if n_elements(lon_max) eq 0 then lon_max=130.
if n_elements(precc_thr) eq 0 then precc_thr=15.
if n_elements(fl_in) eq 0 then fl_in='/data/zender/_aux0_/omega0_2/omega.01.nc'
if n_elements(fl_out) eq 0 then fl_out='/data/zender/_aux0_/cem/dom_avg.nc'
if n_elements(foo) eq 0 then foo=0
;
; Read dimension data from the NetCDF file 
;
nc_id=ncdf_open(fl_in)
ncdf_varget,nc_id,'time',time
ncdf_varget,nc_id,'lat',lat
ncdf_varget,nc_id,'lev',lev
ncdf_varget,nc_id,'lon',lon
;
foo=min(abs(lon-lon_min),lon_min_idx)
foo=min(abs(lon-lon_max),lon_max_idx)
foo=min(abs(lat-lat_min),lat_min_idx)
foo=min(abs(lat-lat_max),lat_max_idx)
;
lat_nbr=(lat_max_idx-lat_min_idx)+1
lon_nbr=(lon_max_idx-lon_min_idx)+1
lev_nbr=(lev_max_idx-lev_min_idx)+1
time_nbr=(time_max_idx-time_min_idx)+1
;
print,'Beginning Region Truncation...'
print,'lon_nbr = '+auto_sng(lon_nbr,1)+' min(lon) = '+auto_sng(min(lon),1)+' max(lon) = '+auto_sng(max(lon),1)
print,'lon_min = '+auto_sng(lon_min,1)+' lon_min_idx = '+auto_sng(lon_min_idx,1)+' = '+auto_sng(lon(lon_min_idx),1)
print,'lon_max = '+auto_sng(lon_max,1)+' lon_max_idx = '+auto_sng(lon_max_idx,1)+' = '+auto_sng(lon(lon_max_idx),1)
print,'lat_nbr = '+auto_sng(lat_nbr,1)+' min(lat) = '+auto_sng(min(lat),1)+' max(lat) = '+auto_sng(max(lat),1)
print,'lat_min = '+auto_sng(lat_min,1)+' lat_min_idx = '+auto_sng(lat_min_idx,1)+' = '+auto_sng(lat(lat_min_idx),1)
print,'lat_max = '+auto_sng(lat_max,1)+' lat_max_idx = '+auto_sng(lat_max_idx,1)+' = '+auto_sng(lat(lat_max_idx),1)
print,''
;
; Extract the correct region from the dimensions
;
lon=lon(lon_min_idx:lon_max_idx)
lat=lat(lat_min_idx:lat_max_idx)
time=time(time_min_idx:time_max_idx)
lev=lev(lev_min_idx:lev_max_idx)
;
; Form the offset and count arrays so that only useful data is read
; into memory (memory is a severe limitation with this 4-dimensional
; dataset).
;
offset=[lon_min_idx,lev_min_idx,lat_min_idx,time_min_idx]
count=[lon_nbr,lev_nbr,lat_nbr,time_nbr]
ncdf_varget,nc_id,'CMFMC',cmfmc,count=count,offset=offset
;
offset=[lon_min_idx,lat_min_idx,time_min_idx]
count=[lon_nbr,lat_nbr,time_nbr]
ncdf_varget,nc_id,'PRECC',precc,count=count,offset=offset
ncdf_varget,nc_id,'PRECT',prect,count=count,offset=offset
ncdf_varget,nc_id,'ORO',oro,count=count,offset=offset
;
; End of NetCDF commands
;
ncdf_close,nc_id
;
; Scale the data
;
cmfmc=1.0e3*cmfmc	; cmfmc is now in g/m^2/s
precc=8.64e7*precc	; precc is now in mm/day
prect=8.64e7*prect	; prect is now in mm/day
;
ocean_idx=where(oro(*,*,0) lt .1,nbr_ocean_idx)
print,'There are '+auto_sng(lat_nbr*lon_nbr,0)+' grid points in the specified region of which '+auto_sng(nbr_ocean_idx,0)+' are over ocean.'
;
print,'Processing '+fl_in+' for '+auto_sng(time_nbr,0)+' time steps: '+auto_sng(time_min_idx,0)+' -- '+auto_sng(time_max_idx,0)
for time_step=0,time_nbr-1 do begin
	foo=where(precc(*,*,time_step) gt precc_thr,nbr_precc_points)
	print,'step = '+auto_sng(time_step,0)+', time = '+auto_sng(time(time_step),3)+', has '+auto_sng(nbr_precc_points,0)+' gridpoints meeting PRECC threshold.'
endfor
;
print,'Processing '+fl_in+' for '+auto_sng(lev_nbr,0)+' vertical levels: '+auto_sng(lev_min_idx,0)+' -- '+auto_sng(lev_max_idx,0)
for level=0,lev_nbr-1 do begin
	print,'level = '+auto_sng(level,0)+', lev = '+auto_sng(lev(level),3)
endfor
;
; If regional averaging is specified then we must carefully average
; spatially over the valid points. This collapses the two spatial
; dimensions into a point so cmfmc will be 2-D (level,time) and 
; precc will be 1-D (time) after this.
;
dom_avg_cmfmc=fltarr(lev_nbr,time_nbr)
dom_avg_precc=fltarr(time_nbr)
dom_avg_prect=fltarr(time_nbr)
for time_slc=0,time_nbr-1 do begin
	for level=0,lev_nbr-1 do begin
;
		slice=cmfmc(*,level,*,time_slc)
		slice=slice(ocean_idx)
		dom_avg_cmfmc(level,time_slc)=total(slice)/nbr_ocean_idx
;
	endfor; end loop over levels
;
	slice=precc(*,*,time_slc)
	slice=slice(ocean_idx)
	dom_avg_precc(time_slc)=total(slice)/nbr_ocean_idx
;
	slice=prect(*,*,time_slc)
	slice=slice(ocean_idx)
	dom_avg_prect(time_slc)=total(slice)/nbr_ocean_idx
;
endfor; end loop over times
;
; Find out if the output file exists.
;
matching_fls=findfile(fl_out)
print,'matching_fls = ',matching_fls
;
if (size(matching_fls))(0) eq 0 then begin
;
; The output file does not yet exist, create it.
;
	time_dim_sz=0
	lev_dim_sz=lev_nbr
	nc_id=ncdf_create(fl_out,/noclobber)
	time_dim_id=ncdf_dimdef(nc_id,'time',/unlimited)
	lev_dim_id=ncdf_dimdef(nc_id,'lev',lev_dim_sz)
	lev_id=ncdf_vardef(nc_id,'lev',[lev_dim_id],/float)
	precc_thr_id=ncdf_vardef(nc_id,'precc_thr',/float)
	lat_min_id=ncdf_vardef(nc_id,'lat_min',/float)
	lat_max_id=ncdf_vardef(nc_id,'lat_max',/float)
	lon_min_id=ncdf_vardef(nc_id,'lon_min',/float)
	lon_max_id=ncdf_vardef(nc_id,'lon_max',/float)
	cmfmc_id=ncdf_vardef(nc_id,'CMFMC',[lev_dim_id,time_dim_id],/float)
	time_id=ncdf_vardef(nc_id,'time',[time_dim_id],/float)
	precc_id=ncdf_vardef(nc_id,'PRECC',[time_dim_id],/float)
	prect_id=ncdf_vardef(nc_id,'PRECT',[time_dim_id],/float)
	ncdf_control,nc_id,/endef
	ncdf_varput,nc_id,lev_id,lev
	ncdf_varput,nc_id,precc_thr_id,precc_thr
	ncdf_varput,nc_id,lat_min_id,lat_min
	ncdf_varput,nc_id,lat_max_id,lat_max
	ncdf_varput,nc_id,lon_min_id,lon_min
	ncdf_varput,nc_id,lon_max_id,lon_max
endif else begin
;
; The output file does exist: gather the variable id's
;
	nc_id=ncdf_open(fl_out,/write)
	time_id=ncdf_varid(nc_id,'time')
	cmfmc_id=ncdf_varid(nc_id,'CMFMC')
	precc_id=ncdf_varid(nc_id,'PRECC')
	prect_id=ncdf_varid(nc_id,'PRECT')
	lev_dim_id=ncdf_dimid(nc_id,'lev')
	time_dim_id=ncdf_dimid(nc_id,'time')
	ncdf_diminq,nc_id,time_dim_id,foo,time_dim_sz
	ncdf_diminq,nc_id,lev_dim_id,foo,lev_dim_sz
	ncdf_varget,nc_id,'time',foo_time
	print,'Previous data processed for '+auto_sng(time_dim_sz,0)+' time steps: ',foo_time(0:time_dim_sz-1)
endelse
;
ncdf_varput,nc_id,time_id,time,offset=[time_dim_sz]
ncdf_varput,nc_id,precc_id,dom_avg_precc,offset=[time_dim_sz]
ncdf_varput,nc_id,prect_id,dom_avg_prect,offset=[time_dim_sz]
ncdf_varput,nc_id,cmfmc_id,dom_avg_cmfmc,offset=[0,time_dim_sz]
;
; End of NetCDF commands
;
ncdf_close,nc_id
;
end_of_procedure: foo=1
end; end cmfmc_dom_avg()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Domain Average analysis
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin CMFMC Gridpoint analysis
; In this analysis procedure, all time and horizontal space information
; is lost. Gridpoints which meet the masking and other criteria are all
; dumped onto the tape level by level.
; Keeping the time information in this procedure would be difficult,
; because there can't be two unlimited dimensions and it's not really
; known beforehand how many times there will be. Making time the 
; unlimited dimension is an option, but that would require getting
; rid of all the interactive selection criteria, i.e., we could keep
; the oro mask but not the precc criteria, since that would result in
; another unknown dimension.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro cmfmc_point, $
	time_min_idx=time_min_idx, $
	time_max_idx=time_max_idx, $
	lev_min_idx=lev_min_idx, $
	lev_max_idx=lev_max_idx, $
	precc_thr=precc_thr, $
	fl_in=fl_in, $
	fl_out=fl_out, $
	rgn_avg=rgn_avg, $
	foo=foo
;
; cmfmc_point,lev_min_idx=7,lev_max_idx=9
; cmfmc_point,time_min_idx=0,time_max_idx=2
; cmfmc_point,lev_min_idx=1,lev_max_idx=7,fl_in='/data/zender/_aux0_/omega0_2/omega.01.nc'
; cmfmc_point,lev_min_idx=1,lev_max_idx=7,time_min_idx=0,time_max_idx=71,fl_in='/data/zender/_aux0_/omega0_2/omega.01.nc'
;
if n_elements(time_min_idx) eq 0 then time_min_idx=0
if n_elements(time_max_idx) eq 0 then time_max_idx=2
if n_elements(lev_min_idx) eq 0 then lev_min_idx=0
if n_elements(lev_max_idx) eq 0 then lev_max_idx=9
if n_elements(lat_min) eq 0 then lat_min=10.
if n_elements(lon_min) eq 0 then lon_min=120.
if n_elements(lat_max) eq 0 then lat_max=20.
if n_elements(lon_max) eq 0 then lon_max=130.
if n_elements(precc_thr) eq 0 then precc_thr=15.
if n_elements(fl_in) eq 0 then fl_in='/data/zender/_aux0_/omega0_2/omega.01.nc'
if n_elements(fl_out) eq 0 then fl_out='/data/zender/_aux0_/cem/analyze.nc'
if n_elements(rgn_avg) eq 0 then rgn_avg=0
if n_elements(foo) eq 0 then foo=0
;
; Read dimension data from the NetCDF file 
;
nc_id=ncdf_open(fl_in)
ncdf_varget,nc_id,'time',time
ncdf_varget,nc_id,'lat',lat
ncdf_varget,nc_id,'lev',lev
ncdf_varget,nc_id,'lon',lon
;
foo=min(abs(lon-lon_min),lon_min_idx)
foo=min(abs(lon-lon_max),lon_max_idx)
foo=min(abs(lat-lat_min),lat_min_idx)
foo=min(abs(lat-lat_max),lat_max_idx)
;
lat_nbr=(lat_max_idx-lat_min_idx)+1
lon_nbr=(lon_max_idx-lon_min_idx)+1
lev_nbr=(lev_max_idx-lev_min_idx)+1
time_nbr=(time_max_idx-time_min_idx)+1
;
print,'Beginning Region Truncation...'
print,'lon_nbr = '+auto_sng(lon_nbr,1)+' min(lon) = '+auto_sng(min(lon),1)+' max(lon) = '+auto_sng(max(lon),1)
print,'lon_min = '+auto_sng(lon_min,1)+' lon_min_idx = '+auto_sng(lon_min_idx,1)+' = '+auto_sng(lon(lon_min_idx),1)
print,'lon_max = '+auto_sng(lon_max,1)+' lon_max_idx = '+auto_sng(lon_max_idx,1)+' = '+auto_sng(lon(lon_max_idx),1)
print,'lat_nbr = '+auto_sng(lat_nbr,1)+' min(lat) = '+auto_sng(min(lat),1)+' max(lat) = '+auto_sng(max(lat),1)
print,'lat_min = '+auto_sng(lat_min,1)+' lat_min_idx = '+auto_sng(lat_min_idx,1)+' = '+auto_sng(lat(lat_min_idx),1)
print,'lat_max = '+auto_sng(lat_max,1)+' lat_max_idx = '+auto_sng(lat_max_idx,1)+' = '+auto_sng(lat(lat_max_idx),1)
print,''
;
; Extract the correct region from the dimensions
;
lon=lon(lon_min_idx:lon_max_idx)
lat=lat(lat_min_idx:lat_max_idx)
time=time(time_min_idx:time_max_idx)
lev=lev(lev_min_idx:lev_max_idx)
;
; Form the offset and count arrays so that only useful data is read
; into memory (memory is a severe limitation with this 4-dimensional
; dataset).
;
offset=[lon_min_idx,lev_min_idx,lat_min_idx,time_min_idx]
count=[lon_nbr,lev_nbr,lat_nbr,time_nbr]
ncdf_varget,nc_id,'CMFMC',cmfmc,count=count,offset=offset
;
offset=[lon_min_idx,lat_min_idx,time_min_idx]
count=[lon_nbr,lat_nbr,time_nbr]
ncdf_varget,nc_id,'PRECC',precc,count=count,offset=offset
ncdf_varget,nc_id,'PRECT',prect,count=count,offset=offset
ncdf_varget,nc_id,'ORO',oro,count=count,offset=offset
;
; End of NetCDF commands
;
ncdf_close,nc_id
;
; Scale the data
;
cmfmc=1.0e3*cmfmc	; cmfmc is now in g/m^2/s
precc=8.64e7*precc	; precc is now in mm/day
prect=8.64e7*prect	; prect is now in mm/day
;
foo=where(oro(*,*,0) eq 0.0,nbr_ocean_points)
print,'There are '+auto_sng(lat_nbr*lon_nbr,0)+' grid points in the specified region of which '+auto_sng(nbr_ocean_points,0)+' are over ocean.'
;
print,'Processing '+fl_in+' for '+auto_sng(time_nbr,0)+' time steps: '+auto_sng(time_min_idx,0)+' -- '+auto_sng(time_max_idx,0)
for time_step=0,time_nbr-1 do begin
	foo=where(precc(*,*,time_step) gt precc_thr,nbr_precc_points)
	print,'step = '+auto_sng(time_step,0)+', time = '+auto_sng(time(time_step),3)+', has '+auto_sng(nbr_precc_points,0)+' gridpoints meeting PRECC threshold.'
endfor
;
print,'Processing '+fl_in+' for '+auto_sng(lev_nbr,0)+' vertical levels: '+auto_sng(lev_min_idx,0)+' -- '+auto_sng(lev_max_idx,0)
for level=0,lev_nbr-1 do begin
	print,'level = '+auto_sng(level,0)+', lev = '+auto_sng(lev(level),3)
endfor
;
; If regional averaging is specified then we must carefully average
; spatially over the valid points. This collapses the two spatial
; dimensions into a point so cmfmc will be 2-D (level,time) and 
; precc will be 1-D (time) after this.
;
good_precc_idx=where(precc gt precc_thr,nbr_good_precc_idx)
if nbr_good_precc_idx eq 0 then begin
	print,'No PRECC matches in current analysis, exiting cmfmc_point()...'
	goto, end_of_procedure
endif
precc=precc(good_precc_idx)
prect=prect(good_precc_idx)
oro=oro(good_precc_idx)
;
good_oro_idx=where(oro eq 0.0,nbr_good_oro_idx)
if nbr_good_oro_idx eq 0 then begin
	print,'No ORO matches in current analysis, exiting cmfmc_point()...'
	goto, end_of_procedure
endif
precc=precc(good_oro_idx)
prect=prect(good_oro_idx)
oro=oro(good_oro_idx)
;
; Find out if the output file exists.
;
matching_fls=findfile(fl_out)
print,'matching_fls = ',matching_fls
;
if (size(matching_fls))(0) eq 0 then begin
;
; The output file does not yet exist, create it.
;
	unlim_dim_sz=0
	time_dim_sz=1008
	lev_dim_sz=lev_nbr
	tot_num_time=0
	nc_id=ncdf_create(fl_out,/noclobber)
	unlim_dim_id=ncdf_dimdef(nc_id,'unlim',/unlimited)
	time_dim_id=ncdf_dimdef(nc_id,'time',time_dim_sz)
	lev_dim_id=ncdf_dimdef(nc_id,'lev',lev_dim_sz)
	lev_id=ncdf_vardef(nc_id,'lev',[lev_dim_id],/float)
	tot_num_time_id=ncdf_vardef(nc_id,'tot_num_time',/long)
	precc_thr_id=ncdf_vardef(nc_id,'precc_thr',/float)
	lat_min_id=ncdf_vardef(nc_id,'lat_min',/float)
	lat_max_id=ncdf_vardef(nc_id,'lat_max',/float)
	lon_min_id=ncdf_vardef(nc_id,'lon_min',/float)
	lon_max_id=ncdf_vardef(nc_id,'lon_max',/float)
	cmfmc_id=ncdf_vardef(nc_id,'CMFMC',[lev_dim_id,unlim_dim_id],/float)
	time_id=ncdf_vardef(nc_id,'time',[time_dim_id],/float)
	precc_id=ncdf_vardef(nc_id,'PRECC',[unlim_dim_id],/float)
	prect_id=ncdf_vardef(nc_id,'PRECT',[unlim_dim_id],/float)
	ncdf_control,nc_id,/endef
	ncdf_varput,nc_id,lev_id,lev
	ncdf_varput,nc_id,precc_thr_id,precc_thr
	ncdf_varput,nc_id,lat_min_id,lat_min
	ncdf_varput,nc_id,lat_max_id,lat_max
	ncdf_varput,nc_id,lon_min_id,lon_min
	ncdf_varput,nc_id,lon_max_id,lon_max
endif else begin
;
; The output file does exist: gather the variable id's
;
	nc_id=ncdf_open(fl_out,/write)
	time_id=ncdf_varid(nc_id,'time')
	cmfmc_id=ncdf_varid(nc_id,'CMFMC')
	precc_id=ncdf_varid(nc_id,'PRECC')
	prect_id=ncdf_varid(nc_id,'PRECT')
	tot_num_time_id=ncdf_varid(nc_id,'tot_num_time')
	lev_dim_id=ncdf_dimid(nc_id,'lev')
	unlim_dim_id=ncdf_dimid(nc_id,'unlim')
	time_dim_id=ncdf_dimid(nc_id,'time')
	ncdf_varget,nc_id,'tot_num_time',tot_num_time
	ncdf_diminq,nc_id,unlim_dim_id,foo,unlim_dim_sz
	ncdf_diminq,nc_id,time_dim_id,foo,time_dim_sz
	ncdf_diminq,nc_id,lev_dim_id,foo,lev_dim_sz
	ncdf_varget,nc_id,'tot_num_time',tot_num_time
	ncdf_varget,nc_id,'time',foo_time
	print,'Previous data processed for '+auto_sng(tot_num_time,0)+' time steps: ',foo_time(0:tot_num_time-1)
endelse
;
; These variables are level-independent, so write them before
; the level loop.
;
ncdf_varput,nc_id,tot_num_time_id,tot_num_time+time_nbr
ncdf_varput,nc_id,time_id,time,offset=[tot_num_time]
ncdf_varput,nc_id,precc_id,precc,offset=[unlim_dim_sz]
ncdf_varput,nc_id,prect_id,prect,offset=[unlim_dim_sz]
good_cmfmc=fltarr(lev_nbr,nbr_good_oro_idx)
;
; Loop over the levels for the entire time.
; This is necessary because the where() function loses the multi-
; dimensional level information, i.e., where() returns single 
; dimensional indices that are degenerate in level.
; Fortunately, we don't care about the time dimension in this
; application so we needn't loop over it as well.
;
for level=0,lev_nbr-1 do begin
;
; Sift out the good values in the same order that the good_????_idx
; were generated by the successive where() commands.
;
cmfmc_lev=cmfmc(*,level,*,*)
;
cmfmc_lev=cmfmc_lev(good_precc_idx)
good_cmfmc(level,*)=cmfmc_lev(good_oro_idx)
;
endfor ; end loop over levels
;
; These variables are level-dependent, so write them after
; the level loop.
;
ncdf_varput,nc_id,cmfmc_id,good_cmfmc,offset=[0,unlim_dim_sz]
;
; End of NetCDF commands
;
ncdf_close,nc_id
;
end_of_procedure: foo=1
end; end cmfmc_point()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End CMFMC Gridpoint analysis
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure CMFMC Graphs
; Histograph the convective mass flux from CCM output.
; If time is the unlimited dimension, then it is assumed that
; the data is domain-averaged and available once per time step.
; In this case the data is also 
; plotted as a time series and precipitation is also graphed.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro cmfmc_graphs, $
	fl_in=fl_in, $
	level=level, $
	graph_prect=graph_prect, $
	cmfmc_bin_sz=cmfmc_bin_sz, $
	ttl_sng=ttl_sng, $
	rng_y=rng_y, $
	foo=foo
;
; cmfmc_graphs,level=6,cmfmc_bin_sz=2.5,fl_in='/data/zender/_aux0_/cem/dom_avg.nc',graph_prect=1
;
if n_elements(graph_prect) eq 0 then graph_prect=0
if n_elements(cmfmc_bin_sz) eq 0 then cmfmc_bin_sz=10.
if n_elements(foo) eq 0 then foo=1
if n_elements(fl_in) eq 0 then fl_in='/data/zender/_aux0_/cem/dom_avg.nc'
if n_elements(ttl_sng) eq 0 then ttl_sng=''
if n_elements(rng_y) le 1 then rng_y=0
if n_elements(level) eq 0 then level=4
;
; Read dimension data from the NetCDF file 
;
nc_id=ncdf_open(fl_in)
foo=ncdf_inquire(nc_id)
unlim_dim_id=foo.recdim
if unlim_dim_id eq -1 then begin
	print,'There is no unlimited dimension in '+fl_in+', exiting.'
	goto,end_of_procedure
endif
;
ncdf_diminq,nc_id,unlim_dim_id,unlim_dim_nm,unlim_dim_sz
if unlim_dim_nm eq 'time' then DOM_AVG=1 else DOM_AVG=0
;
lev_dim_id=ncdf_dimid(nc_id,'lev')
ncdf_diminq,nc_id,lev_dim_id,foo,lev_dim_sz
ncdf_varget,nc_id,'lev',lev
ncdf_varget,nc_id,'time',time
if not DOM_AVG then ncdf_varget,nc_id,'tot_num_time',tot_num_time
ncdf_varget,nc_id,'precc_thr',precc_thr
ncdf_varget,nc_id,'lat_min',lat_min
ncdf_varget,nc_id,'lat_max',lat_max
ncdf_varget,nc_id,'lon_min',lon_min
ncdf_varget,nc_id,'lon_max',lon_max
;
offset=[level,0]
count=[1,unlim_dim_sz]
ncdf_varget,nc_id,'CMFMC',cmfmc,count=count,offset=offset
;
if DOM_AVG then begin
offset=[0]
count=[unlim_dim_sz]
ncdf_varget,nc_id,'PRECT',prect,count=count,offset=offset
endif; end if DOM_AVG
;
ncdf_close,nc_id
;
; End of NetCDF commands
;
data=cmfmc*1. ; data is already in g/m^2/s
if not DOM_AVG then begin
	sub_ttl_sng='Source: OMEGA July instantaneous oceanic gridpoint CMFMC where PRECC > '+auto_sng(precc_thr,0)+' mm/day.'
endif else begin
	sub_ttl_sng='Source: OMEGA July instantaneous oceanic domain average CMFMC.'
endelse
;
generic_hst, $
	abb_sng='Conv. Mass Flux', $
	bin_sz=cmfmc_bin_sz, $
	clean=1, $
	data=data, $
	lat_max=lat_max, $
	lat_min=lat_min, $
	lon_max=lon_max, $
	lon_min=lon_min, $
	data_max=max(data), $
	data_min=min(data), $
	data_nbr=n_elements(data), $
	rgn_sng='', $
	slice_sng=auto_sng(round(lev(level)),0)+' mb', $
	sub_ttl_sng=sub_ttl_sng, $
	ttl_sng=ttl_sng, $
	sym_sng='!8M!5!Ic!N', $
	unit_sng='!5g m!e-2!N s!E-1!N', $
	rng_y=rng_y
;
if DOM_AVG and GRAPH_PRECT then begin
;
print,' Hit any key to continue...'
junk = get_kbrd(1)
;
data=prect*1. ; data is already in mm/day
;
slice_sng=''
sym_sng='!8P!5!IT!N'
unit_sng='!5mm day!E-1!N'
;
generic_hst, $
	abb_sng='Domain Avg. Precip', $
	bin_sz=1.0, $
	clean=1, $
	data=data, $
	lat_max=lat_max, $
	lat_min=lat_min, $
	lon_max=lon_max, $
	lon_min=lon_min, $
	data_max=max(data), $
	data_min=min(data), $
	data_nbr=n_elements(data), $
	rgn_sng='', $
	slice_sng=slice_sng, $
	sub_ttl_sng=sub_ttl_sng, $
	sym_sng=sym_sng, $
	ttl_sng=ttl_sng, $
	unit_sng=unit_sng, $
	rng_y=rng_y
;
print,' Hit any key to continue...'
junk = get_kbrd(1)
;
plt_rgn_nrm=[ $
	.13, $ ; x_min
	.11, $ ; y_min
	.87, $ ; x_max
	.90] ; y_max
;
plot, $
	time, $
	data, $
	tit='!5Domain-Averaged Precipitation vs. Time', $
	xtit='Model Day !8t!5 (days)', $
	ytit=slice_sng+' '+sym_sng+' ('+unit_sng+')', $
	xstyle=1, $
	ystyle=0, $
	/ynozero, $
	thick=3.0, $
	charsize=2.0, $
	position=plt_rgn_nrm, $
	linestyle=0
;
print,' Hit any key to continue...'
junk = get_kbrd(1)
;
time_nbr=n_elements(time)
incr_time_min=20.
;
generic_pow_spec, $
	data=data, $
	time_nbr=time_nbr, $
	delta_time=incr_time_min*60.0, $
	x_ttl_sng=spec_x_ttl_sng, $
	y_ttl_sng=spec_y_ttl_sng
;
endif; end if DOM_AVG and GRAPH_PRECT
;
end_of_procedure: foo=1
;
end; end cmfmc_graphs()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure CMFMC Graphs commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
