From: afl@cdc.noaa.gov (Andy Loughe)
Newsgroups: comp.lang.idl-pvwave
Subject: Re: HDF viewer available?
Date: 16 Aug 1995 16:44:16 GMT
Organization: Climate Diagnostic Center
Message-ID: <40t790$pk3@lace.Colorado.EDU>
References: <DDD2oC.2uH@midway.uchicago.edu> <40qna9$pl3@lace.Colorado.EDU> <DDDFCz.DJt@midway.uchicago.edu>


To view the contents of an HDF SDS one can use netCDF's ncdump utility, 
or hdp which I believe is a new HDF dump utility from NCSA.

Here is something I wrote for specific HDF SDS files.
I must add that I have been very impressed with RSI's 
new HDF interface.  Good job!


;
;  Function to obtain information on NASA Goddard's GEOS-1 HDF files.
;
;  Originator:  Andrew F. Loughe
;
;     History:  Written 4/28/95
;  

FUNCTION GET_HDF_INFO, filename, varname, $
         varnames=varnames, coordnames=coordnames, $
         lons=lons, lats=lats, levels=levels, times=times, $
         title=title, units=units, format=format, id=id, dims=dims

on_error, 2

;  Help string
USAGE1='err = get_hdf_info(filename, varname, ' + $
       'varnames=varnames, coordnames=coordnames,' 
USAGE2='          lons=lons, lats=lats, levels=levels, times=times,'
USAGE3='          title=title, units=units, format=format, id=id, ' + $
       'dims=dims)'

;  Need at least one parameter passed in.
if (N_params() lt 1) then begin
   print, ' '
   print, USAGE1
   print, USAGE2
   print, USAGE3
   return, -1
endif

;  Supply nonsense varname.
if ( n_elements(varname) eq 0 ) then varname='ZxYzXy'
varname = strlowcase(varname)

;  Assign dummy values to ERR and ID.
ERR = -1 
ID  = -9999

;	Check that file is in HDF.
if HDF_ISHDF( filename ) ne 1 then $
	message, "File " + filename + " is not an HDF file."

;  Open HDF file and initialize the SD interface.
sd_id = HDF_SD_START( filename, /read )

;Get number of datasets and global attributes in this file.
HDF_SD_FILEINFO, sd_id, nmfsds, nglobatts

;  Loop through all datasets until correct SDS_NAME (and SDS_ID) is found.
if nmfsds gt 0 then begin
   varnames = ''
   coordnames = ''
   for i = 0, nmfsds-1 do begin
      sds_id = HDF_SD_SELECT(sd_id, i)
      HDF_SD_GETINFO, sds_id, name=n, ndims=r, type=t, $
                      natts=nats, dims=dimensions

;  Collect all variable names.
      if ( HDF_SD_ISCOORDVAR(sds_id) ne 1 ) then $ 
         varnames = [varnames, n + ' ; '] 
 
;  If desired variable name is found, get more information.
      if ( strpos( strlowcase(n) , varname) ge 0 ) then begin   
         ERR  = 0 
         ID   = sds_id
         dims = dimensions
         for j = 0, nats-1 do begin
             HDF_SD_ATTRINFO, sds_id, j, name=n, data=d
             if (n eq 'long_name') then title=d
             if (n eq 'units')     then units=d
             if (n eq 'format')    then format=d
         endfor
      endif


;  Get information on coordinate variables.
      if ( HDF_SD_ISCOORDVAR(sds_id) eq 1 ) then begin
         coordnames = [coordnames, n + ' ; ']
         if ( strpos(strlowcase(n), 'longitude') ge 0 ) then $
            HDF_SD_GETDATA, sds_id, lons
         if ( strpos(strlowcase(n), 'latitude') ge 0 ) then $
            HDF_SD_GETDATA, sds_id, lats
         if ( strpos(strlowcase(n), 'level') ge 0 ) then $
            HDF_SD_GETDATA, sds_id, levels
         if ( strpos(strlowcase(n), 'time') ge 0 ) then $
            HDF_SD_GETDATA, sds_id, times
      endif

      HDF_SD_ENDACCESS, sds_id
   endfor
endif

;  Omit first varname and coordname which is set to a dummy value.
if ( N_elements(varnames) gt 0 ) then $
   varnames=varnames( indgen(N_elements(varnames)-1) + 1 )

if ( N_elements(coordnames) gt 0 ) then $
   coordnames=coordnames( indgen(N_elements(coordnames)-1) + 1 )


;  Close the HDF file.
HDF_SD_END, sd_id


return, err
end
