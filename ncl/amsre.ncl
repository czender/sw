; $Id$

; Purpose: Read and plot AMSR-E soil moisture data
; Equal Area Scalable Earth Grid (EASE) Grid
; AMSR-E Level 2/3 data stored on ML grid http://nsidc.org/data/ease/ease_grid.html

; Example usage:

; ncl < ~/ncl/amsre.ncl
; ncview ${DATA}/amsre/amsre_200406_1x1.nc &
; load "/home/zender/ncl/amsre.ncl"
; tst("x11")

begin ; begin amsre_cnv()

; Initialize some fields
prg_nm="amsre_cnv"
dbg_lvl=0
mss_val=1.0e36

; Get data
ttl_sng="Volumetric Soil Moisture" ; 
; File suffix .heg is not recognized
; File suffix .hdfeos is recognized
; 200507 NSIDC changed file suffix changed to .hdf
fl_in=addfile("/data/zender/amsre/amsre_200406_1x1.hdf","r")
;print(fl_in)
mss_cnc_H2O_mgxcm3=fl_in->soil_moisture
mss_val_in=mss_cnc_H2O_mgxcm3@NoRetrieval
dmn_sz=dimsizes(mss_cnc_H2O_mgxcm3)
lat_nbr=dmn_sz(0)
lon_nbr=dmn_sz(1)
time_nbr=1
;lat=new((/lat_nbr/),float,mss_val_in)
lat=fspan(89.5,-89.5,lat_nbr)
lat!0="lat"
lat@long_name="Latitude"
lat@units="degrees_north"
;lon=new((/lon_nbr/),float,mss_val_in)
lon=fspan(-180,180,lon_nbr)
lon!0="lon"
lon@long_name="Longitude"
lon@units="degrees_east"
time=1
time!0="time"
time@long_name="Time"
time@units="none"

; Re-grid data
; Original data are ordered N->S and start at 180.0
; First re-order latitude from N->S to S->N
; Second, move first longitude from 180.0 to 0.0


; Create netCDF output
fl_nm_out="/data/zender/amsre/amsre_200406_1x1.nc"
xst_flg=isfilepresent(fl_nm_out)
if (xst_flg) then
	system ("/bin/rm -f "+fl_nm_out) ; Remove any pre-existing file
end if
; "c" creates without overwriting, "w" writes to existing
fl_out=addfile(fl_nm_out,"c")

; Define File Dimensions
delete(dmn_sz)
dmn_nm=(/"time","lat","lon"/)
dmn_sz=(/time_nbr,lat_nbr,lon_nbr/)
dmn_ulm=(/True,False,False/)
filedimdef(fl_out,dmn_nm,dmn_sz,dmn_ulm)

; Define Global Attributes
nl=integertochar(10) ; newline character
att_glb=True ; Variable to hold global attributes
att_glb@history=systemfunc("date")+": ncl < amsre.ncl"
att_glb@source="Received 20050405 from Siri Jodha Singh Khalsa <sjsk@nsidc.org> (303) 492-1445"
att_glb@notes="Okay, here is try #2.  The grid is still in the orientation of the input, i.e. UL is +90N, 180W.  All missing data are set to 0 (we will change this in the future to something like -9999.)"+ \
"I've attached a plot comparing the input and output grids, with a detail of S. America.  Note the large areas of the Amazon Basin where there are no retrievals because of dense vegetation. These get smoothed out in the process of resampling. We may modify the resampling method to yield more representative results in these areas."+ \
"Also note the anomalously high values along the eastern coast. This is most likely an antenna sidelobe contamination effect. In the future we could apply a mask that would remove some margin along the coasts, plus eliminate the stray island and Antarctic coastal pixels."+ \
"The units are g/cm**3 after you divide by 1000 (or percent volumetric if you divide by 10).  The range of possible values is .015 to .5 g/cm**3."+ \
"A new version of the product is being generated now at the facility in Huntsville. We expect there to be improvements, as some calibration problems with the low frequency channels have been addressed."+ \
"I checked the HDF file with the hdp utility and it looks fine."
fileattdef(fl_out,att_glb)

; Convert from weird units to SI
mss_val_out=1.0e36 ; fxm Change to 1.0e36
; Conveniently, input mass concentration is mg cm-3 = kg m3
mss_cnc_H2O=(/mss_cnc_H2O_mgxcm3/)*1.0e-6*1.0e6 ; kg m-3 = mg->kg * cm3->m3
mss_cnc_H2O!0="lat"
mss_cnc_H2O!1="lon"
; Convert old mss_val_in to mss_val_out
; Yes, NCL does "the right thing" when re-assigning missing_values
mss_cnc_H2O@missing_value=mss_val_out

; Handle weird data values
if(1.eq.1) then
; This loop does not do so well at converting missing values
;do lat_idx=0,lat_nbr-1
;	do lon_idx=0,lon_nbr-1
;		if(ismissing(mss_cnc_H2O(lat_idx,lon_idx))) then
;		      	mss_cnc_H2O(lat_idx,lon_idx)=mss_val_out
;	        end if
;	end do; end loop over lon_idx
;end do; end loop over lat_idx
mss_cnc_H2O=mss_cnc_H2O < 1000.0
mss_cnc_H2O=mss_cnc_H2O > 0.0
else
; Find and fix data less than 0.0 kg m-3
mss_cnc_H2O_1D=ndtooned(mss_cnc_H2O)
err_idx=ind(mss_cnc_H2O_1D.lt.0.0)
err_idx_2D=ind_resolve(err_idx,dimsizes(mss_cnc_H2O))
; Set these data to missing_value
mss_cnc_H2O_1D(err_idx)=mss_val_out
; Re-map corrected 1D array into 2D array
mss_cnc_H2O=onedtond(mss_cnc_H2O_1D,dimsizes(mss_cnc_H2O))
; fxm: other_lt0 shows where original values were not -9999, but were < 0
other_lt0=0.0*new((/lat_nbr,lon_nbr/),float,mss_val_out)
other_lt0_1D=ndtooned(other_lt0)
other_lt0_1D(err_idx)=1.0
other_lt0=onedtond(other_lt0_1D,dimsizes(mss_cnc_H2O))

; Find and fix data greater than 1000.0 kg m-3
if (num(mss_cnc_H2O_1D.gt.1000.0).gt.0) then
mss_cnc_H2O_1D=ndtooned(mss_cnc_H2O)
err_idx=ind(mss_cnc_H2O_1D.gt.1.0)
; Set these data to missing_value
err_idx_2D=ind_resolve(err_idx,dimsizes(mss_cnc_H2O))
; Set these data to missing_value
mss_cnc_H2O_1D(err_idx)=mss_val_out
; Re-map corrected 1D array into 2D array
mss_cnc_H2O=onedtond(mss_cnc_H2O_1D,dimsizes(mss_cnc_H2O))
; fxm: other_gt0 shows where original values were not were > 1.0
other_gt1=0.0*new((/lat_nbr,lon_nbr/),float,mss_val_out)
other_gt1_1D=ndtooned(other_gt1)
other_gt1_1D(err_idx)=1.0
other_gt1=onedtond(other_gt1_1D,dimsizes(mss_cnc_H2O))
end if
end if

; Create output variable
vwc_sfc=new((/time_nbr,lat_nbr,lon_nbr/),float,mss_val_out)
vwc_sfc!0="time"
vwc_sfc!1="lat"
vwc_sfc!2="lon"
vwc_sfc@long_name="Volumetric Water Content at surface"
vwc_sfc@units="meter3 meter-3"
vwc_sfc@missing_value=mss_val_out

; Define output variable
; Convert H2O mass concentration (kg H2O m-3 soil) to VWC (m3 H2O m-3 soil)
vwc_sfc(0,:,:)=(/mss_cnc_H2O/)*1.0e-3 ; m3 m-3 = kg H2O m-3 soil / (1000 kg H2O m-3 H2O) 

; Write output file
fl_out->time=time
fl_out->lon=lon
fl_out->lat=lat
fl_out->vwc_sfc=vwc_sfc
;fl_out->other_gt1=other_gt1
;fl_out->other_lt0=other_lt0

end ; end amsre()

