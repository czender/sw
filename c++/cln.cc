// $Id$ 

// Purpose: Calendar properties

/* Copyright (C) 1997--present Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

#include <cln.hh> // Calendar properties

yyyymmdd_sct // O [yyyymmdd]
yyyymmdd_prs(const long yyyymmdd){
  // Purpose: Parses date (in yyyymmdd format) into constituent year, month, and day
  // Usage: yyyymmdd_trp=yyyymmdd_prs(yyyymmdd);
  yyyymmdd_sct yyyymmdd_trp;

  //  assert(yyyymmdd < 1.0e9);
  
  yyyymmdd_trp.dd=yyyymmdd%100; // [day] Day [1..31]
  yyyymmdd_trp.mm=((yyyymmdd%10000)-yyyymmdd_trp.dd)/100; // [mth] Month [1..12]
  yyyymmdd_trp.yyyy=int(yyyymmdd/10000); // [yyyy] Year

  return yyyymmdd_trp;
} // end yyyymmdd_prs()

long // O [yyyymmdd] Date
yyyymmdd2lng
(const short yyyy, // I [yyyy] Year in yyyy format
 const short mm, // I [mth] Month [1..12]
 const short dd) // I [day] Day [1..31]
{
  // Purpose: Combines constituent year, month, and day into single yyyymmdd
  // Usage: yyyymmdd=yyyymmdd2lng(yyyy,mm,dd);
  using cln::mpy; // (12) [nbr] Model months per year
  const short yyyymmdd_sng_lng(8); // [nbr] Length of yyyymmdd string
  char yyyymmdd_sng[yyyymmdd_sng_lng+1]; // [sng] Date string
  long yyyymmdd; // O [yyyymmdd] Date
  assert(mm >= 1 && mm <= mpy);
  assert(dd >= 1 && dd <= 31);
  long chr_nbr=std::sprintf(yyyymmdd_sng,"%04d%02d%02d",yyyy,mm,dd);
  if(chr_nbr != yyyymmdd_sng_lng){
    err_prn("yyyymmdd2lng","Printed too many characters to yyyymdd_sng");
  } // endif
  yyyymmdd=std::strtol(yyyymmdd_sng,(char **)NULL,10);
  return yyyymmdd;
} // end yyyymmdd2lng()

long // O [yyyymmdd]
yyyymmdd2lng(const yyyymmdd_sct yyyymmdd) // I [yyyymmdd] Date
{
  // Purpose: Combines constituent year, month, and day into single yyyymmdd
  // Usage: yyyymmdd=yyyymmdd2lng(yyyy,mm,dd);
  // yyyymmdd2lng() is the inverse of yyyymmdd_prs()
  using cln::mpy; // (12) [nbr] Model months per year
  const short yyyymmdd_sng_lng(8); // [nbr] Length of yyyymmdd string
  char yyyymmdd_sng[yyyymmdd_sng_lng+1]; // [sng] Date string
  long yyyymmdd_lng; // O [yyyymmdd] Date
  assert(yyyymmdd.mm >= 1 && yyyymmdd.mm <= mpy);
  assert(yyyymmdd.dd >= 1 && yyyymmdd.dd <= 31);
  long chr_nbr=std::sprintf(yyyymmdd_sng,"%04d%02d%02d",yyyymmdd.yyyy,yyyymmdd.mm,yyyymmdd.dd);
  if(chr_nbr != yyyymmdd_sng_lng){
    err_prn("yyyymmdd2lng","Printed too many characters to yyyymdd_sng");
  } // endif
  yyyymmdd_lng=std::strtol(yyyymmdd_sng,(char **)NULL,10);
  return yyyymmdd_lng;
} // end yyyymmdd2lng()

bool // O [flg]
is_leap_yr(const short yyyy){
  // Purpose: Returns true if given year is a leap year else returns false
  // Usage: leap_yr_flg=is_leap_yr(yyyy);
    bool leap_yr_flg(false); 

    if(yyyy%4 == 0 && yyyy%100 != 0) leap_yr_flg=true; // 1900 is not a leap year
    if(yyyy%400 == 0) leap_yr_flg=true; // 2000 is a leap year

    return leap_yr_flg;
} // end is_leap_yr()

short 
dpy_get(const short yyyy)
{
  // Purpose: Returns # of days in given year (365 or 366)
  // Usage: dpy=dpy_get(yyyy);
    return is_leap_yr(yyyy) ? 366 : 365;
} // end dpy_get()

short
dpm_get(const short mm)
{
  // Purpose: Returns # of days in given month (no leap days)
  // Assume mm is a two digit integer ranging from 1..12
  // Usage: dpm=dpm_get(mm)
  using cln::mpy; // (12) [nbr] Model months per year
  const short dpm[mpy]={31,28,31,30,31,30,31,31,30,31,30,31};
  return dpm[mm-1]; // Subtract one for base-0 indexing
} // end dpm_get()

long // O [yyyymmdd] Date at end
yyyymmdd_ncr
(const long yyyymmdd_srt, // I [yyyymmdd] Date at start
 const long day_nbr) // I [day] Days to increment
{
  /* Purpose: Increment a date by a given number of days
     Usage: yyyymmdd_new=yyyymmdd_ncr(yyyymmdd_srt,day_nbr);
     day_nbr may be positive (days in the future) or negative (days in the past)
     day_nbr should be an integer */
  long yyyymmdd_new; // [yyyymmdd] New date
  short yyyy_srt; // [yr] Year in start date
  short mm_srt; // [mth] Month in start date [1..12]
  short dd_srt; // [day] Day in start date [1..31]
  short dd_prv; // [day] Provisional day in new date [1..31]
  short yyyy_new; // [yr] Year in new date
  short mm_new; // [mth] Month in new date [1..12]
  short dd_new; // [day] Day in new date [1..31]
  short dpm_srt; // [day] Days in current month
  yyyymmdd_sct yyyymmdd_trp_srt; // [yyyymmdd] 
  // Parse current date
  yyyymmdd_trp_srt=yyyymmdd_prs(yyyymmdd_srt);
  // Set defaults for new date
  yyyy_srt=yyyymmdd_trp_srt.yyyy; // [yr] Year in start date
  mm_srt=yyyymmdd_trp_srt.mm; // [mth] Month in start date [1..12]
  dd_srt=yyyymmdd_trp_srt.dd; // [day] Day in start date [1..31]
  yyyy_new=yyyy_srt; // [yr] Year in new date
  mm_new=mm_srt; // [mth] Month in new date [1..12]
  dd_new=dd_srt; // [day] Day in new date [1..31]
  // Add days to current date
  dpm_srt=dpmyr(mm_srt,yyyy_srt);
  dd_prv=dd_srt+day_nbr; // [day] Provisional day
  if(dd_srt+day_nbr <= dpm_srt && dd_prv >= 1){
    dd_new=dd_srt+day_nbr; // [day] Day in new date [1..31]
    yyyymmdd_new=yyyymmdd2lng(yyyy_new,mm_new,dd_new);
  }else{
    prc_cmp doy_srt=yyyymmdd2doy(yyyymmdd_srt);
    prc_cmp doy_new=doy_srt+day_nbr;
    short dpy_new=dpy_get(yyyy_new);
    if(day_nbr >= 0){
      while(doy_new > dpy_new){
	doy_new-=dpy_new;
	yyyy_new++;
	dpy_new=dpy_get(yyyy_new);
      } // end while
    }else{ // endif day_nbr >= 0
      while(doy_new < 1){
		yyyy_new--;
		dpy_new=dpy_get(yyyy_new);
		doy_new+=dpy_new;
      } // endwhile
    } // endif day_nbr < 0
    yyyymmdd_new=yrdoy2yyyymmdd(yyyy_new,doy_new);
  } // endelse
  return yyyymmdd_new;
} // end yyyymmdd_ncr()

long // O [yyyymmdd] Date at end
yrdoy2yyyymmdd
(const short yyyy, // I [yyyy] Year 
 const prc_cmp doy) // I [day] Day of year
{
  // Purpose: Convert year and day of year [1.0..367.0) into yyyymmdd format
  // yrdoy2yyyymmdd() is inverse of yyyymmdd2doy()
  // Usage: yyyymmdd=yrdoy2yyyymmdd(yyyy,doy);
  using cln::mpy; // (12) [nbr] Model months per year
  long yyyymmdd; // O [yyyymmdd] Date determined by yyyy and doy inputs
  short mm; // [mm] Month determined by yyyy and doy inputs
  short dd; // [mm] Day determined by yyyy and doy inputs
  // Number of days in all preceding months of non-leap year  
  short day_ttl_prc_mth_non_leap_yr[mpy]={0,31,59,90,120,151,181,212,243,273,304,334};
  // Number of days in all preceding months of a leap year
  short day_ttl_prc_mth_leap_yr[mpy]={0,31,60,91,121,152,182,213,244,274,305,335}; 
  short *day_ttl_prc_mth; // Number of days in all preceding months
  // Vet input
  std::string sbr_nm("yrdoy2yyyymmdd"); // [sng] Subroutine name
  if(doy > dpy_get(yyyy)){
    std::cerr << prg_nm_get() << ": ERROR " << sbr_nm << "() reports doy = " << doy << " in yyyy = " << yyyy << std::endl;
    err_prn(sbr_nm,"Exiting");
  } // endif
  // Compute calendar for this year
  bool leap_yr_flg=is_leap_yr(yyyy); // [flg] Leap year flag
  if(leap_yr_flg) day_ttl_prc_mth=day_ttl_prc_mth_leap_yr; else day_ttl_prc_mth=day_ttl_prc_mth_non_leap_yr;
  mm=1; // Start at January
  while(mm < 13 && doy > day_ttl_prc_mth[mm-1]){ // -1 is for 0-based indices
    // printf STDOUT "INFO yrdoy2yyyymmdd() reports \$yyyy = $yyyy, \$doy = $doy, \$leap_yr_flg = $leap_yr_flg, \$mm = $mm, \$day_ttl_prc_mth[$mm-1] = $day_ttl_prc_mth[$mm-1]\n";
    mm++;
  } // end while mm
  mm--; // Compensate for overshoot
  dd=static_cast<short>(doy-day_ttl_prc_mth[mm-1]); 
  if(dd < 1 || dd > dpmyr(mm,yyyy)){
    std::cerr << prg_nm_get() << ": ERROR " << sbr_nm << "() reports yyyy = " << yyyy << ", doy = " << doy << ", mm = " << mm << ", dd = " << dd << ", day_ttl_prc_mth[mm-1] = " << day_ttl_prc_mth[mm-1] << std::endl;
    err_prn(sbr_nm,"Exiting");
  } // endif err
  yyyymmdd=yyyymmdd2lng(yyyy,mm,dd);
  return yyyymmdd;
} // end yrdoy2yyyymmdd()

prc_cmp // O [day] Day of year
yyyymmdd2doy(const long yyyymmdd) // I [yyyymmdd] Date   
{
  // Purpose: Convert date in yyyymmdd format to day of year [1.0..367.0)
  // yyyymmdd2doy() is inverse of yrdoy2yyyymmdd()
  // Usage: doy=yyyymmdd2doy(yyyymmdd);
  short yyyy; // [yr] Year in date
  short mm; // [mth] Month in date [1..12]
  short dd; // [day] Day in date [1..31]
  prc_cmp doy; // O [day] Day of year
  yyyymmdd_sct yyyymmdd_trp; // [yyyymmdd] Date
  // Parse current date
  yyyymmdd_trp=yyyymmdd_prs(yyyymmdd);
  yyyy=yyyymmdd_trp.yyyy; // [yr] Year in date
  mm=yyyymmdd_trp.mm; // [mth] Month in date [1..12]
  dd=yyyymmdd_trp.dd; // [day] Day in date [1..31]
  // Add up days in preceding months
  using cln::mpy; // (12) [nbr] Model months per year
  const short day_ttl_prc_mth[mpy]={0,31,59,90,120,151,181,212,243,273,303,334}; // Number of days in all preceding months of non-leap year
  doy=dd+day_ttl_prc_mth[mm-1];
  // Correct for leap years
  if(mm > 2)
    if(is_leap_yr(yyyy)) 
      doy++;
  //  print STDOUT "yyyymmdd2doy() reports \$yyyymmdd = $yyyymmdd, \$doy = $doy\n";
  return doy;
} // end yyyymmdd2doy()

short // O [day] Days in month
dpmyr
(const short mm, // I [mth] Month [1..12]
 const short yyyy) // I [yyyy] Year
{
  // Purpose: Returns # of days in given month (accounts for leap days)
  // Usage: dpm=dpmyr(mm,yyyy)
  // Assume mm is a two digit integer ranging from 1..12
  // Assume yyyy is in YYYY format
  using cln::mpy; // (12) [nbr] Model months per year
  const short dpm_365[mpy]={31,28,31,30,31,30,31,31,30,31,30,31};
  short dpm;
  dpm=dpm_365[mm-1]; // Subtract one for base-0 indexing
  if(mm == 2){
    if(is_leap_yr(yyyy)) dpm++;
  } // endif February
  return dpm;
} // end dpmyr()

std::string // O [sng] Date in yyyymmdd format
date2vrs // [fnc] Convert date string into yyyymmdd
(const std::string &date_sng) // [sng] POSIX date
{
  /* Purpose: Convert date string into version string
     Expects date_sng in format "Sat Oct 13 19:13:08 PDT 2001"
     Returns version string in format "20011213"
     Similar to date --iso-8601
     date --rfc-822
     date +%Y%m%d.%H%M%S */
  // const long yyyymmdd(date2yyyymmdd(date_sng)); // [yyyymmdd] Date
  // const std::string sng(static_cast<std::string>(yyyymmdd));
  //  sng.assign(yyyymmdd);
  // How to convert a long to a string?
  std::string sng="1";
  return sng;
} // end date2vrs()

long // O [yyyymmdd] Date
date2yyyymmdd // [fnc] Convert date string into yyyymmdd
(const std::string &date_sng) // [sng] POSIX date
{
  /* Purpose: Convert date string into yyyymmdd format long
     Expects date_sng in format "Sat Oct 13 19:13:08 PDT 2001"
     Returns numberic representation in format "20011213"
     date --iso-8601
     date --rfc-822
     date +%Y%M%d.%H%M%S */
  long yyyymmdd=std::strtol(date_sng.c_str(),(char **)NULL,10); // [yyyymmdd] Date
  return yyyymmdd;
} // end date2yyyymmdd()
