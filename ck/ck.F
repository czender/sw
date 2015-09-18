c     $Id$

      program ck

#if ( defined DOUBLE )
#define COMPUTATIONAL_PRECISION double precision
c#define NF_REAL NF_DOUBLE
#else
c#define NF_DOUBLE NF_REAL
#define COMPUTATIONAL_PRECISION real
#endif

c     Purpose: Compute correlated-k parameters from HITRAN line data.

c     NB: Compiling and running program in single precision may cause inexplicable
c     crashes and floating point errors. This depends on gas species and wavelength
c     interval because some HITRAN line strengths are simply too small
c     ( < 10^-36) to represent in single precision. Even when line strength can be 
c     represented in single precision, the arithmatic below (performed in SI units) 
c     may cause IEEE underflow etc. Thus our strategy is to read and process the 
c     data in double precision, and then store results in single precision.
c     This requires compiling with double precision options, e.g., use -r8 -i4 on Suns.
c     Makefile implements this with
c     make -W ck.F -k precision=double ck

c     Example usage:
c     For H2O vapor (default) one can just use ck.
c     Other gases require explicitly setting # of bands and names and molecule #:

c     Debugging:
c     ck -m 7 -n -t 28 -l 1.e-35 -h 1. -b 2 -i /data2/zender/hitran/foo.nc -o /data2/zender/hitran/foo_out.nc

c     Water Vapor:

      implicit none
c     Parameters
      character*(*) 
     $     rcs_Id,
     $     rcs_Revision

      integer
     $     fl_in_unit,
     $     bnd_nbr_max,
     $     max_ln_nbr,
     $     max_t_nbr

      COMPUTATIONAL_PRECISION
     $     t_Malkmus_BPB        ! reference temperature for Malkmus parameters

      parameter(fl_in_unit=73,
     $     bnd_nbr_max=4800,    ! 2000-50000 cm-1 in 10cm-1 bands
     $     max_ln_nbr=44228,    ! O3 from .2-5 microns
     $     max_t_nbr=140,
     $     rcs_Id='$Id$',
     $     rcs_Revision='$Revision$',
     $     t_Malkmus_BPB=250.)

c     Commons
#include <parameter.com>
#include <netcdf.inc>
c     Input Arguments
c     Input/Output Arguments
c     Output Arguments
c     Local workspace

      character
     $     argv*80,
     $     char_foo*4,
     $     fl_in*80,
     $     lbl*80,
     $     fl_out*80
      character*26::lcl_date_time ! Time formatted as Day Mth DD HH:MM:SS TZ YYYY

      integer
     $     arg,
     $     dbg_lvl,
     $     bnd_dbg,
     $     exit_status,         ! program exit status
     $     int_foo,
     $     nbr_arg,
     $     rcd                ! return success code

      logical
     $     netCDF_input

      integer
     $     bnd_dim_id,          ! dimension ID for band
     $     bnd_idx,             ! counting index
     $     ln_dim_id,           ! dimension ID for ln
     $     ln_idx,              ! counting index
     $     bnd_nbr,             ! dimension size
     $     ln_nbr,              ! dimension size
     $     nc_id,               ! file handle
     $     srt(1)               ! starting offsets (long *) 

c     HITRAN netCDF input variables
      integer
     $     iso_nbr_id,
     $     ln_str_id,
     $     nrg_lo_id,
     $     mlc_nbr_id,
     $     HWHM_id,
     $     HWHM_t_dep_id,
     $     ln_id

      integer*2                 ! NF_SHORT
     $     iso_nbr,
     $     mlc_nbr

      COMPUTATIONAL_PRECISION
     $     ln_str(max_ln_nbr),
     $     nrg_lo(max_ln_nbr),
     $     HWHM(max_ln_nbr),
     $     HWHM_t_dep(max_ln_nbr),
     $     ln(max_ln_nbr)

      integer
     $     bnd_id,             ! coordinate ID
     $     bnd_sz_id,
     $     bnd_ln_nbr_id,
     $     wvl_ctr_id,
     $     wvl_max_id,
     $     wvl_min_id,
     $     wvl_sz_id,
     $     wvn_ctr_id,
     $     wvn_max_id,
     $     wvn_min_id,
     $     wvn_sz_id
     
      COMPUTATIONAL_PRECISION
     $     bnd(bnd_nbr_max),  ! coordinate variable
     $     bnd_sz,
     $     float_foo,
     $     wvl_ctr(bnd_nbr_max),
     $     wvl_max(bnd_nbr_max),
     $     wvl_min(bnd_nbr_max),
     $     wvl_sz(bnd_nbr_max),
     $     wvn_ctr(bnd_nbr_max),
     $     wvn_max(bnd_nbr_max),
     $     wvn_min(bnd_nbr_max),
     $     wvn_sz(bnd_nbr_max)

c     Single-precision netCDF output variables
      real*4                    ! NF_REAL
     $     bnd_r4(bnd_nbr_max), ! coordinate variable
     $     bnd_sz_r4,
     $     wvl_ctr_r4(bnd_nbr_max),
     $     wvl_max_r4(bnd_nbr_max),
     $     wvl_min_r4(bnd_nbr_max),
     $     wvl_sz_r4(bnd_nbr_max),
     $     wvn_ctr_r4(bnd_nbr_max),
     $     wvn_max_r4(bnd_nbr_max),
     $     wvn_min_r4(bnd_nbr_max),
     $     wvn_sz_r4(bnd_nbr_max)

c Local variables
      character*10
     $     iso_sng,
     $     mlc_sng

      COMPUTATIONAL_PRECISION
     $     ln_hi,
     $     ln_lo

      integer
     $     ln_idx_max(bnd_nbr_max),
     $     ln_idx_min(bnd_nbr_max),
     $     bnd_ln_nbr(bnd_nbr_max)

c     Externals
      integer 
     $     iargc
c     Main code

c     Initialize default values
      dbg_lvl=0
      exit_status=0
      fl_in='/data/zender/hitran/H2O.nc'
      fl_out='/data/zender/aca/ck.nc'
      ln_lo=2000.
      ln_hi=17900.
      netcdf_input=.true.
      mlc_nbr=1

c     Retrieve command line arguments
      call date_time_get(lcl_date_time)
      nbr_arg=iargc()
      do arg=1,nbr_arg
         call getarg(arg,argv)
         if (argv(1:2).eq.'-B') then
            call getarg(arg+1,argv)
            read (argv,'(i4)') bnd_dbg
         endif
         if (argv(1:2).eq.'-b') then
            call getarg(arg+1,argv)
            read(argv,'(i4)') bnd_nbr
         endif
         if (argv(1:2).eq.'-D') then
            call getarg(arg+1,argv)
            read(argv,'(i4)') dbg_lvl
         endif
         if (argv(1:2).eq.'-f') then
            call getarg(arg+1,argv)
            read(argv,'(f8.3)') float_foo
         endif
         if (argv(1:2).eq.'-h') then
            call getarg(arg+1,argv)
            read(argv,'(f8.3)') ln_hi
         endif
         if (argv(1:2).eq.'-i') then
            call getarg(arg+1,argv)
            read(argv,'(a)') fl_in
         endif
         if (argv(1:2).eq.'-l') then
            call getarg(arg+1,argv)
            read(argv,'(f8.3)') ln_lo
         endif
         if (argv(1:2).eq.'-m') then
            call getarg(arg+1,argv)
            read(argv,'(i)') mlc_nbr
         endif
         if (argv(1:2).eq.'-n') then
            netCDF_input=.not.netCDF_input
         endif
         if (argv(1:2).eq.'-o') then
            call getarg(arg+1,argv)
            read(argv,'(a)') fl_out
         endif
         if (argv(1:2).eq.'-v') then
            write(6,'(a)') rcs_Id
            goto 1000
         endif
      end do

      if (netCDF_input) then 
         write(6,'(a)') 'Input presumed to be netCDF format HITRAN line data'
      else 
         write(6,'(a)') 'Input presumed to be BPB-style binary Malkmus random band data' 
      endif

c     Set netCDF error handling
      call ncpopt(NCVERBOS+NCFATAL)
      
c     Read in netCDF data.
      nc_id=ncopn(fl_in,NCNOWRIT,rcd)
      
c     Get global attributes
      call ncagtc(nc_id,NF_GLOBAL,'molecule',mlc_sng,len(mlc_sng),rcd)
      call ncagtc(nc_id,NF_GLOBAL,'isotope',iso_sng,len(iso_sng),rcd)
      do bnd_idx=1,len(mlc_sng)
         if (mlc_sng(bnd_idx:bnd_idx).eq.char(0)) goto 600
      end do
 600  do bnd_idx=bnd_idx,len(mlc_sng)
         mlc_sng(bnd_idx:bnd_idx)=char(0)
      end do
      do bnd_idx=1,len(iso_sng)
         if (iso_sng(bnd_idx:bnd_idx).eq.char(0)) goto 700
      end do
 700  do bnd_idx=bnd_idx,len(iso_sng)
         iso_sng(bnd_idx:bnd_idx)=char(0)
      end do
      
c     Get dimension IDs
      ln_dim_id=ncdid(nc_id,'wvn_ctr',rcd)
      
c     Get dimension sizes
      call ncdinq(nc_id,ln_dim_id,lbl,ln_nbr,rcd)
      
c     Get variable IDs
      ln_id=ncvid(nc_id,'wvn_ctr',rcd)
      mlc_nbr_id=ncvid(nc_id,'mlc_nbr',rcd)
      iso_nbr_id=ncvid(nc_id,'iso_nbr',rcd)
      ln_str_id=ncvid(nc_id,'ln_str',rcd)
      HWHM_id=ncvid(nc_id,'HWHM',rcd)
      nrg_lo_id=ncvid(nc_id,'nrg_lo',rcd)
      HWHM_t_dep_id=ncvid(nc_id,'HWHM_t_dep',rcd)
      
c     Get data
      srt(1)=1
      
      call ncvgt1(nc_id,mlc_nbr_id,1,mlc_nbr,rcd)
      call ncvgt1(nc_id,iso_nbr_id,1,iso_nbr,rcd)
      call ncvgt(nc_id,ln_id,srt,ln_nbr,ln,rcd)
      call ncvgt(nc_id,ln_str_id,srt,ln_nbr,ln_str,rcd)
      call ncvgt(nc_id,HWHM_id,srt,ln_nbr,HWHM,rcd)
      call ncvgt(nc_id,nrg_lo_id,srt,ln_nbr,nrg_lo,rcd)
      call ncvgt(nc_id,HWHM_t_dep_id,srt,ln_nbr,HWHM_t_dep,rcd)
      
      call ncclos(nc_id,rcd)
      write(6,'(a20,1x,a60)') 'Ingested netCDF file',fl_in
      
c     Convert input data to SI units where necessary. Normally i never store
c     anything in units other than SI in a netCDF file, but this HITRAN input is
c     an exception. Line data is exactly as it appears in the HITRAN database,
c     which is basically CGS. The reason why ln_str should be scaled to SI, even
c     though there is a possible loss of precision involved for weak lines in 
c     single precision, is that ln_str is used in all sorts of SI algebra 
c     involving fundamental constants and it's easiest to do it all in SI.
c     Once ln_str is scaled to SI, the band parameters (S_d, S_p, ...) computed
c     from ln_str will be in SI.
      do ln_idx=1,ln_nbr
         ln_str(ln_idx)=ln_str(ln_idx)*1.e-4 ! cm2/mlc -> m2/mlc
         ln_str(ln_idx)=ln_str(ln_idx)*Avagadro ! m2/mlc -> m2/mol
      enddo                     ! end loop over lines
      
c     Compute wavenumber coordinates
      bnd_sz=(ln_hi-ln_lo)/bnd_nbr
      float_foo=bnd_sz/2.
      do bnd_idx=1,bnd_nbr
         wvn_min(bnd_idx)=ln_lo+(bnd_idx-1)*bnd_sz
      enddo                     ! end loop over bnd
      do bnd_idx=1,bnd_nbr-1
         wvn_max(bnd_idx)=wvn_min(bnd_idx+1)
      enddo                     ! end loop over bnd
      wvn_max(bnd_nbr)=ln_hi
      do bnd_idx=1,bnd_nbr
         wvn_ctr(bnd_idx)=wvn_min(bnd_idx)+float_foo
      enddo                     ! end loop over bnd
      
      if (dbg_lvl.gt.1) then
         write (6,'(2(a8,f9.2,a5,/),a10,i5,/,a10,i3)') 
     $        'ln_lo = ',ln_lo,' cm-1',
     $        'ln_hi = ',ln_hi,' cm-1',
     $        'ln_nbr = ',ln_nbr,
     $        'bnd_nbr = ',bnd_nbr
      endif                     ! endif dbg
      
c     Find first line within specified spectral interval.
c     This index bootstraps next loop.
      do ln_idx=1,ln_nbr
         if (ln(ln_idx).ge.ln_lo.and.ln(ln_idx).le.ln_hi) goto 500
         if (ln_idx.eq.ln_nbr) then
            write(6,'(a41,/,a27,g,a5,/,a26,g,a5,a8)') 
     $           'No lines found within specified interval.',
     $           'First line in input data = ',ln(1),' cm-1',
     $           'Last line in input data = ',ln(ln_nbr),' cm-1',
     $           'Exiting.'
            exit_status=-1
            call exit(exit_status)
         endif
      enddo                     ! end loop over lines
 500  continue
      
c     Find first and last lines in each band. We assume that on entry ln_idx points
c     to beginning of first band, and that ln is monotonically increasing.
      do bnd_idx=1,bnd_nbr
         do while (ln_idx.le.ln_nbr.and.ln(ln_idx).ge.wvn_min(bnd_idx).and.ln(ln_idx).lt.wvn_max(bnd_idx))
            if (ln_idx_min(bnd_idx).eq.0) ln_idx_min(bnd_idx)=ln_idx
            ln_idx_max(bnd_idx)=ln_idx
            ln_idx=ln_idx+1
         end do                 ! end loop over lines within each band
      enddo                     ! end loop over bnd
      int_foo=0
      do bnd_idx=1,bnd_nbr
         bnd_ln_nbr(bnd_idx)=ln_idx_max(bnd_idx)-ln_idx_min(bnd_idx)+1
         int_foo=int_foo+bnd_ln_nbr(bnd_idx)
      enddo                     ! end loop over bnd
      
c     Now we know which lines are in each band.
      write (6,'(a29,f12.6,a9,f12.6,a18,f12.6,a11)') 
     $     'c-k parameters computed from ',
     $     ln_lo,' cm-1 to ',ln_hi,' cm-1 for regular ',bnd_sz,' cm-1 bands'
      write (6,'(i6,a47)') int_foo,' lines fall within specified spectral interval.'
      write (6,'(a76,i4,a3,f12.6,a8,f12.6,a9)') 
     $     'Fitted and exact weak and strong line temperature dependence saved for band ',bnd_dbg,
     $     ' = ',wvn_ctr(bnd_dbg),' cm-1 = ',wvl_ctr(bnd_dbg)*1.e6,' microns.'
      
      do bnd_idx=1,bnd_nbr
         wvl_max(bnd_idx)=1./(100.*wvn_min(bnd_idx))
         wvl_min(bnd_idx)=1./(100.*wvn_max(bnd_idx))
         wvl_ctr(bnd_idx)=
     $        .5*(wvl_min(bnd_idx)+wvl_max(bnd_idx))
         wvl_sz(bnd_idx)=
     $        wvl_max(bnd_idx)-wvl_min(bnd_idx)
         wvn_sz(bnd_idx)=
     $        wvn_max(bnd_idx)-wvn_min(bnd_idx)
         bnd(bnd_idx)=wvn_ctr(bnd_idx)
      enddo                     ! end loop over bnd
      
c     Shoehorn computational-precision variables into 4-byte floats for storage as netCDF NF_REALs.
      bnd_sz_r4=bnd_sz
      do bnd_idx=1,bnd_nbr
         
         bnd_r4(bnd_idx)=bnd(bnd_idx)
         wvl_ctr_r4(bnd_idx)=wvl_ctr(bnd_idx)
         wvl_max_r4(bnd_idx)=wvl_max(bnd_idx)
         wvl_min_r4(bnd_idx)=wvl_min(bnd_idx)
         wvl_sz_r4(bnd_idx)=wvl_sz(bnd_idx)
         wvn_ctr_r4(bnd_idx)=wvn_ctr(bnd_idx)
         wvn_max_r4(bnd_idx)=wvn_max(bnd_idx)
         wvn_min_r4(bnd_idx)=wvn_min(bnd_idx)
         wvn_sz_r4(bnd_idx)=wvn_sz(bnd_idx)
      enddo                     ! end loop over bnd
      
c     Begin NetCDF output routines
      nc_id=nccre(fl_out,NF_CLOBBER,rcd)
      
c     Define dimension IDs
      bnd_dim_id=ncddef(nc_id,'bnd',bnd_nbr,rcd)
      
c     Variable definitions
      bnd_id=ncvdef(nc_id,'bnd',NF_REAL,1,bnd_dim_id,rcd)
      bnd_ln_nbr_id=ncvdef(nc_id,'bnd_ln_nbr',NF_INT,1,bnd_dim_id,rcd)
      bnd_sz_id=ncvdef(nc_id,'bnd_sz',NF_REAL,0,0,rcd)
      iso_nbr_id=ncvdef(nc_id,'iso_nbr',NF_SHORT,0,0,rcd)
      mlc_nbr_id=ncvdef(nc_id,'mlc_nbr',NF_SHORT,0,0,rcd)
      wvl_ctr_id=ncvdef(nc_id,'wvl_ctr',NF_REAL,1,bnd_dim_id,rcd)
      wvl_max_id=ncvdef(nc_id,'wvl_max',NF_REAL,1,bnd_dim_id,rcd)
      wvl_min_id=ncvdef(nc_id,'wvl_min',NF_REAL,1,bnd_dim_id,rcd)
      wvl_sz_id=ncvdef(nc_id,'wvl_sz',NF_REAL,1,bnd_dim_id,rcd)
      wvn_ctr_id=ncvdef(nc_id,'wvn_ctr',NF_REAL,1,bnd_dim_id,rcd)
      wvn_max_id=ncvdef(nc_id,'wvn_max',NF_REAL,1,bnd_dim_id,rcd)
      wvn_min_id=ncvdef(nc_id,'wvn_min',NF_REAL,1,bnd_dim_id,rcd)
      wvn_sz_id=ncvdef(nc_id,'wvn_sz',NF_REAL,1,bnd_dim_id,rcd)
      
c     Add global attibutes
      call ncaptc(nc_id,NF_GLOBAL,'rcs_Id',NF_CHAR,len(rcs_Id),rcs_Id,rcd)
      call ncaptc(nc_id,NF_GLOBAL,'creation_date',NF_CHAR,len(lcl_date_time),lcl_date_time,rcd)
      call ncaptc(nc_id,NF_GLOBAL,'molecule',NF_CHAR,len(mlc_sng),mlc_sng,rcd)
      call ncaptc(nc_id,NF_GLOBAL,'isotope',NF_CHAR,len(iso_sng),iso_sng,rcd)
      
c     Add English text descriptions
      call ncaptc(nc_id,bnd_id,'long_name',NF_CHAR,11,'Band center',rcd)
      call ncaptc(nc_id,bnd_ln_nbr_id,'long_name',NF_CHAR,30,'# of HITRAN lines in each band',rcd)
      call ncaptc(nc_id,bnd_sz_id,'long_name',NF_CHAR,22,'Uniform width of bands',rcd)
      call ncaptc(nc_id,iso_nbr_id,'long_name',NF_CHAR,28,'HITRAN isotope number (1..9)',rcd)
      call ncaptc(nc_id,mlc_nbr_id,'long_name',NF_CHAR,30,'HITRAN molecule number (1..31)',rcd)
      call ncaptc(nc_id,wvl_ctr_id,'long_name',NF_CHAR,27,'Midpoint wavelength in band',rcd)
      call ncaptc(nc_id,wvl_max_id,'long_name',NF_CHAR,26,'Maximum wavelength in band',rcd)
      call ncaptc(nc_id,wvl_min_id,'long_name',NF_CHAR,26,'Minimum wavelength in band',rcd)
      call ncaptc(nc_id,wvl_sz_id,'long_name',NF_CHAR,13,'Size of band in wavelength',rcd)
      call ncaptc(nc_id,wvn_ctr_id,'long_name',NF_CHAR,26,'Maximum wavenumber in band',rcd)
      call ncaptc(nc_id,wvn_max_id,'long_name',NF_CHAR,26,'Maximum wavenumber in band',rcd)
      call ncaptc(nc_id,wvn_min_id,'long_name',NF_CHAR,26,'Minimum wavenumber in band',rcd)
      call ncaptc(nc_id,wvn_sz_id,'long_name',NF_CHAR,24,'Size of band in wavenumbers',rcd)
      
c     Add units
      call ncaptc(nc_id,bnd_id,'units',NF_CHAR,12,'centimeter-1',rcd)
      call ncaptc(nc_id,bnd_ln_nbr_id,'units',NF_CHAR,8,'cardinal',rcd)
      call ncaptc(nc_id,bnd_sz_id,'units',NF_CHAR,12,'centimeter-1',rcd)
      call ncaptc(nc_id,wvl_ctr_id,'units',NF_CHAR,5,'meter',rcd)
      call ncaptc(nc_id,wvl_max_id,'units',NF_CHAR,5,'meter',rcd)
      call ncaptc(nc_id,wvl_min_id,'units',NF_CHAR,5,'meter',rcd)
      call ncaptc(nc_id,wvl_sz_id,'units',NF_CHAR,5,'meter',rcd)
      call ncaptc(nc_id,wvn_ctr_id,'units',NF_CHAR,12,'centimeter-1',rcd)
      call ncaptc(nc_id,wvn_max_id,'units',NF_CHAR,12,'centimeter-1',rcd)
      call ncaptc(nc_id,wvn_min_id,'units',NF_CHAR,12,'centimeter-1',rcd)
      call ncaptc(nc_id,wvn_sz_id,'units',NF_CHAR,12,'centimeter-1',rcd)
      
c     Now that all dimensions, variables, and attributes have been defined, make call to end define mode
      call ncendf(nc_id,rcd)
      
c     Write data
      srt(1)=1
      
      call ncvpt1(nc_id,bnd_sz_id,srt,bnd_sz_r4,rcd)
      call ncvpt1(nc_id,iso_nbr_id,srt,iso_nbr,rcd)
      call ncvpt1(nc_id,mlc_nbr_id,srt,mlc_nbr,rcd)
      
      call ncvpt(nc_id,bnd_id,srt,bnd_nbr,bnd_r4,rcd)
      call ncvpt(nc_id,bnd_ln_nbr_id,srt,bnd_nbr,bnd_ln_nbr,rcd)
      call ncvpt(nc_id,wvl_ctr_id,srt,bnd_nbr,wvl_ctr_r4,rcd)
      call ncvpt(nc_id,wvl_max_id,srt,bnd_nbr,wvl_max_r4,rcd)
      call ncvpt(nc_id,wvl_min_id,srt,bnd_nbr,wvl_min_r4,rcd)
      call ncvpt(nc_id,wvl_sz_id,srt,bnd_nbr,wvl_sz_r4,rcd)
      call ncvpt(nc_id,wvn_ctr_id,srt,bnd_nbr,wvn_ctr_r4,rcd)
      call ncvpt(nc_id,wvn_max_id,srt,bnd_nbr,wvn_max_r4,rcd)
      call ncvpt(nc_id,wvn_min_id,srt,bnd_nbr,wvn_min_r4,rcd)
      call ncvpt(nc_id,wvn_sz_id,srt,bnd_nbr,wvn_sz_r4,rcd)
      
      call ncclos(nc_id,rcd)
      
      write (6,'(a16,1x,a40)') 'Wrote results to',fl_out
      
 1000 continue
      
      call exit(exit_status)
      end
      
      
      
      
      
      
      
      



