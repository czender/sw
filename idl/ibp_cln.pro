; $Id$

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Translation of calendar routines in /c++/cln.cc
; ibp_cln.pro
; Dave Newman
; Nov 2001
;
; Notes:
; 1. cc routines left in code (and commented out) for validation purposes (can delete)
; 2. Noteworthy comments marked with <djn>
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin dpmyr
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function dpmyr,mm,yyyy
; Purpose: Returns # of days in given month (accounts for leap days)
; Usage: dpm = dpmyr(mm, yyyy)
; Assume mm is a two digit integer ranging from 1..12
; Assume yyyy is in YYYY format
dpm_365 = [31,28,31,30,31,30,31,31,30,31,30,31]
dpm     = dpm_365(mm-1)
if ((mm eq 2) and (is_leap_yr(yyyy))) then dpm = dpm + 1
return,dpm
end; dpmyr()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End dpmyr()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin yyyymmdd_prs
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro yyyymmdd_prs,yyyymmdd,yr=yr,mth=mth,day=day
; Purpose: Deconstruct date in yyyymmdd format into yr, mth, and day
; yyyymmdd is input
; yr,mth,day are output
; Usage: yyyymmdd_prs,yyyymmdd_srt,yr=yr_srt,mth=mth_srt,day=day_srt

day = yyyymmdd mod 100                 ; Start day
mth = ((yyyymmdd mod 10000) - day)/100 ; Start month
yr  = yyyymmdd/10000                   ; Start year
end; yyyymmdd_prs()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End yyyymmdd_prs()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin yymmdd_prs
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro yymmdd_prs,yymmdd,yr=yr,mth=mth,day=day
; Purpose: Deconstruct date in yymmdd format into yr, mth, and day
; yymmdd is input
; yr,mth,day are output
; Usage: yymmdd_prs,yymmdd_srt,yr=yr_srt,mth=mth_srt,day=day_srt
day=yymmdd mod 100 ; Start day
mth=((yymmdd mod 10000) - day)/100 ; Start month
yr=yymmdd/10000 ; Start year
if yr eq 0 then yr=2000
if yr lt 100 then yr=yr+1900
end; yymmdd_prs()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End yymmdd_prs()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin yyyymmdd2lng
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function yyyymmdd2lng,yyyy,mm,dd
; Purpose: Combines constituent year, month, and day into single yyyymmdd
; Usage: yyyymmdd = yyyymmdd2lng(yyyy, mm, dd)
if ((yyyy lt 1) or (yyyy gt 9999)) then print,'ERROR: yyyymmdd2lng'
if ((mm   lt 1) or (mm   gt 12))   then print,'ERROR: yyyymmdd2lng'
if ((dd   lt 1) or (dd   gt 31))   then print,'ERROR: yyyymmdd2lng'
yyyy     = long(yyyy)
mm       = long(mm)
dd       = long(dd) 
yyyymmdd = 10000*yyyy + 100*mm + dd
return,yyyymmdd
end; yyyymmdd2lng()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End yyyymmdd2lng()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin is_leap_yr
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function is_leap_yr,yyyy
; Purpose: Returns 1 (true) if given year is a leap year else returns 0 (false)
; Usage: leap_yr_flg = is_leap_yr(yyyy);
leap_yr_flg = 0
if ((yyyy mod 4 eq 0) and (yyyy mod 100 ne 0)) then leap_yr_flg = 1
if (yyyy mod 400 eq 0)                         then leap_yr_flg = 1
return,leap_yr_flg
end; is_leap_yr()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End is_leap_yr()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin dpy_get
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function dpy_get,yyyy
; Purpose: Returns # of days in given year (365 or 366)
; Usage: dpy = dpy_get(yyyy)
dpy = 365
if (is_leap_yr(yyyy)) then dpy = 366
return,dpy
end; dpy_get()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End dpy_get()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin dpm_get
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function dpm_get,mm
; Purpose: Returns # of days in given month (no leap days)
; Assume mm is a two digit integer ranging from 1..12
; Usage: dpm = dpm_get(mm)
dpm = [31,28,31,30,31,30,31,31,30,31,30,31]
return,dpm(mm-1)
end; dpm_get()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End dpm_get()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin yyyymmdd_ncr
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function yyyymmdd_ncr,yyyymmdd_srt,day_nbr
; Purpose: Increment a date by a given number of days
; Usage: yyyymmdd_new = yyyymmdd_ncr(yyyymmdd_srt, day_nbr)
;        day_nbr may be positive (days in the future) or negative (days in the past)
;        day_nbr should be an integer

; Parse current date
yyyymmdd_prs,yyyymmdd_srt,yr=yyyy_srt,mth=mm_srt,day=dd_srt

; Set defaults for new date
; yyyy_srt            [yr]  Year in start date
; mm_srt              [mth] Month in start date [1..12]
; dd_srt              [day] Day in start date [1..31]
yyyy_new = yyyy_srt ; [yr]  Year in new date
mm_new   = mm_srt   ; [mth] Month in new date [1..12]
dd_new   = dd_srt   ; [day] Day in new date [1..31]

; Add days to current date
dpm_srt = dpmyr(mm_srt, yyyy_srt)
dd_prv  = dd_srt + day_nbr            ;[day] Provisional day

if (((dd_srt+day_nbr) le dpm_srt) and (dd_prv ge 1)) then begin
    dd_new = dd_srt + day_nbr         ;[day] Day in new date [1..31]
    yyyymmdd_new = yyyymmdd2lng(yyyy_new, mm_new, dd_new)
endif else begin
    doy_srt = yyyymmdd2doy(yyyymmdd_srt)
    doy_new = doy_srt + day_nbr
    dpy_new = dpy_get(yyyy_new)
    if (day_nbr ge 0) then begin
        while (doy_new gt dpy_new) do begin
            doy_new  = doy_new - dpy_new
            yyyy_new = yyyy_new + 1
            dpy_new  = dpy_get(yyyy_new)
        endwhile
    endif else begin                  ;endif day_nbr >= 0
        while (doy_new lt 1) do begin
	    yyyy_new = yyyy_new - 1
	    dpy_new  = dpy_get(yyyy_new)
	    doy_new  = doy_new + dpy_new
        endwhile
    endelse                           ;endif day_nbr < 0
    yyyymmdd_new = yrdoy2yyyymmdd(yyyy_new, doy_new)
endelse

return,yyyymmdd_new
end; yyyymmdd_ncr()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End yyyymmdd_ncr()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin yrdoy2yyyymmdd
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function yrdoy2yyyymmdd,yyyy,doy
; Purpose: Convert year and day of year [1.0..366.0) into yyyymmdd format
; yrdoy2yyyymmdd() is inverse of yyyymmdd2doy()
; Usage: yyyymmdd = yrdoy2yyyymmdd(yyyy, doy)


; Number of days in all preceding months of non-leap year
; <djn> make array length 13 because idl does not short circuit boolean logic
day_ttl_prc_mth_non_leap_yr = [0,31,59,90,120,151,181,212,243,273,304,334,365]
; Number of days in all preceding months of a leap year
day_ttl_prc_mth_leap_yr     = [0,31,60,91,121,152,182,213,244,274,305,335,366]
; Vet input
if (doy gt dpy_get(yyyy)) then print,'ERROR: yrdoy2yyyymmdd'

;<djn> can't figure pointer to array, so do long way in if ... block
; Compute calendar for this year
;if (is_leap_yr(yyyy)) then day_ttl_prc_mth = day_ttl_prc_mth_leap_yr else day_ttl_prc_mth = day_ttl_prc_mth_non_leap_yr

if (is_leap_yr(yyyy)) then begin
    mm = 1
    while ((mm lt 13) and (doy gt day_ttl_prc_mth_leap_yr(mm-1))) do mm = mm + 1
    mm = mm - 1
    dd = doy - day_ttl_prc_mth_leap_yr(mm-1)
    if ((dd lt 1) or (dd gt dpmyr(mm,yyyy))) then print,'ERROR: yrdoy2yyyymmdd'
endif else begin
    mm = 1
    while ((mm lt 13) and (doy gt day_ttl_prc_mth_non_leap_yr(mm-1))) do mm = mm + 1
    mm = mm - 1
    dd = doy - day_ttl_prc_mth_non_leap_yr(mm-1)
    if ((dd lt 1) or (dd gt dpmyr(mm,yyyy))) then print,'ERROR: yrdoy2yyyymmdd'
endelse

yyyymmdd = yyyymmdd2lng(yyyy, mm, dd)
return,yyyymmdd
end; yrdoy2yyyymmdd()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End yrdoy2yyyymmdd()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin yyyymmdd2doy
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function yyyymmdd2doy,yyyymmdd
; Purpose: Convert date in yyyymmdd format to day of year [1.0..366.0)
; yyyymmdd2doy() is inverse of yrdoy2yyyymmdd()
; Usage: doy = yyyymmdd2doy(yyyymmdd);
yyyymmdd_prs,yyyymmdd,yr=yyyy,mth=mm,day=dd
day_ttl_prc_mth = [0,31,59,90,120,151,181,212,243,273,303,334]
doy = dd + day_ttl_prc_mth(mm-1)
if ((mm gt 2) and (is_leap_yr(yyyy))) then doy = doy + 1
return,doy
end; yyyymmdd2doy()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End yyyymmdd2doy()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin date2vrs
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function date2vrs
print,'ERROR: date2vrs not implemented'
return,0
end; date2vrs()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End date2vrs()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;long // O [yyyymmdd] Date
;date2yyyymmdd // [fnc] Convert date string into yyyymmdd
;(const std::string &date_sng) // [sng] POSIX date
;{
;  /* Purpose: Convert date string into yyyymmdd format long
;     Expects date_sng in format "Sat Oct 13 19:13:08 PDT 2001"
;     Returns numberic representation in format "20011213"
;     date --iso-8601
;     date --rfc-822
;     date +%Y%M%d.%H%M%S */
;  long yyyymmdd=std::strtol(date_sng.c_str(),(char **)NULL,10); // [yyyymmdd] Date
;  return yyyymmdd;
;} // end date2vrs()

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin date2yyyymmdd
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function date2yyyymmdd
print,'ERROR: date2yyyymmdd not implemented'
return,0
end; date2yyyymmdd()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End date2yyyymmdd()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
