C     FORTRAN TEMPLATE FOR FILE= /usr/tmp/zender/spcp_85_8589_amip5_8589_01_r15.nc                       
csz   INCLUDE'netcdf.inc'
csz   PARAMETER (NVARS=32) !NUMBER OF VARIABLES
      PARAMETER (NREC=     1)   !CHANGE THIS TO GENERALIZE
C     VARIABLE IDS RUN SEQUENTIALLY FROM 1 TO NVARS= 32
      INTEGER*4 RCODE
      INTEGER*4 RECDIM
C     ****VARIABLES FOR THIS NETCDF FILE****
C
      REAL             P0                   
      INTEGER          ntrm                 
      INTEGER          ntrn                 
      INTEGER          ntrk                 
      INTEGER          ndbase               
      INTEGER          nsbase               
      INTEGER          nbdate               
      INTEGER          nbsec                
      INTEGER          mdt                  
      INTEGER          mhisf                
      CHARACTER*1      current_mss          (  81)
      CHARACTER*1      first_mss            (  81)
      CHARACTER*1      init_mss             (  81)
      CHARACTER*1      tibnd_mss            (  81)
      CHARACTER*1      sst_mss              (  81)
      CHARACTER*1      ozone_mss            (  81)
      REAL             hyai                 (  19)
      REAL             hybi                 (  19)
      REAL             hyam                 (  18)
      REAL             hybm                 (  18)
      REAL             lev                  (  18)
      REAL             ilev                 (  19)
      REAL             lat                  (  40)
      REAL             gw                   (  40)
      REAL             lon                  (  48)
      INTEGER          date                 (NREC)
      INTEGER          datesec              (NREC)
      DOUBLE PRECISION time                 (NREC)
      REAL             QR                   (  48,  18,  40,NREC)
      REAL             DTCOND               (  48,  18,  40,NREC)
      REAL             PS                   (  48,  40,NREC)
      REAL             PRES                 (  48,  18,  40,NREC)
C*************************************
      INTEGER*4 START(10)
      INTEGER*4 COUNT(10)
      INTEGER VDIMS(10) !ALLOW UP TO 10 DIMENSIONS
      CHARACTER*31 DUMMY
c++csz
      integer
     $     in_fl_unit,
     $     out_fl_unit,
     $     max_nbr_lat,
     $     max_nbr_lev,
     $     max_nbr_lon,
     $     max_nbr_time
      parameter(in_fl_unit=73,
     $     out_fl_unit=74,
     $     max_nbr_lat=40,
     $     max_nbr_lev=18,
     $     max_nbr_lon=48,
     $     max_nbr_time=1)

      character
     $     in_fl*80,
     $     out_fl*80

      integer
     $     dbg_lvl,
     $     idx_lat,             ! counting index
     $     idx_lev,             ! counting index
     $     idx_lon,             ! counting index
     $     idx_time,            ! counting index
     $     nbr_lat,             ! dimension size
     $     nbr_lev,             ! dimension size
     $     nbr_lon,             ! dimension size
     $     nbr_time             ! dimension size

      integer
     $     idx_ctn,
     $     idx_var,
     $     idx_lev9
      real
     $     dst,
     $     dst_ctn,
     $     dst_frc_ctn

      real sigma9(9)
      data sigma9 /.009,.074,.189,.336,.500,.664,.811,.926,.991/

      real sigma18(48,18,40,NREC)
      real qr9(48,9,40,NREC)
      real dtcond9(48,9,40,NREC)
      real record(48,40)
c--csz
      NCID=NCOPN(
     +'/usr/tmp/zender/spcp_85_8589_amip5_8589_01_r15.nc           '
     +,NCNOWRIT,RCODE)
      CALL NCINQ(NCID,NDIMS,NVARS,NGATTS,RECDIM,RCODE)
      CALL NCDINQ(NCID,RECDIM,DUMMY,NRECS,RCODE)
C     !NRECS! NOW CONTAINS NUM RECORDS FOR THIS FILE
C
C    statements to fill           P0                   
C
      CALL NCVINQ(NCID, 1,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO  10 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
  10  CONTINUE
      CALL NCVGT(NCID, 1,START,COUNT,
     +          P0                   ,RCODE)
C
C    statements to fill           ntrm                 
C
      CALL NCVINQ(NCID, 2,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO  20 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
  20  CONTINUE
      CALL NCVGT(NCID, 2,START,COUNT,
     +          ntrm                 ,RCODE)
C
C    statements to fill           ntrn                 
C
      CALL NCVINQ(NCID, 3,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO  30 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
  30  CONTINUE
      CALL NCVGT(NCID, 3,START,COUNT,
     +          ntrn                 ,RCODE)
C
C    statements to fill           ntrk                 
C
      CALL NCVINQ(NCID, 4,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO  40 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
  40  CONTINUE
      CALL NCVGT(NCID, 4,START,COUNT,
     +          ntrk                 ,RCODE)
C
C    statements to fill           ndbase               
C
      CALL NCVINQ(NCID, 5,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO  50 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
  50  CONTINUE
      CALL NCVGT(NCID, 5,START,COUNT,
     +          ndbase               ,RCODE)
C
C    statements to fill           nsbase               
C
      CALL NCVINQ(NCID, 6,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO  60 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
  60  CONTINUE
      CALL NCVGT(NCID, 6,START,COUNT,
     +          nsbase               ,RCODE)
C
C    statements to fill           nbdate               
C
      CALL NCVINQ(NCID, 7,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO  70 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
  70  CONTINUE
      CALL NCVGT(NCID, 7,START,COUNT,
     +          nbdate               ,RCODE)
C
C    statements to fill           nbsec                
C
      CALL NCVINQ(NCID, 8,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO  80 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
  80  CONTINUE
      CALL NCVGT(NCID, 8,START,COUNT,
     +          nbsec                ,RCODE)
C
C    statements to fill           mdt                  
C
      CALL NCVINQ(NCID, 9,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO  90 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
  90  CONTINUE
      CALL NCVGT(NCID, 9,START,COUNT,
     +          mdt                  ,RCODE)
C
C    statements to fill           mhisf                
C
      CALL NCVINQ(NCID,10,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 100 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 100  CONTINUE
      CALL NCVGT(NCID,10,START,COUNT,
     +          mhisf                ,RCODE)
C
C    statements to fill           current_mss          
C
      CALL NCVINQ(NCID,11,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 110 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 110  CONTINUE
      CALL NCVGTC(NCID,11,START,COUNT,
     +          current_mss          ,LENSTR,RCODE)
C
C    statements to fill           first_mss            
C
      CALL NCVINQ(NCID,12,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 120 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 120  CONTINUE
      CALL NCVGTC(NCID,12,START,COUNT,
     +          first_mss            ,LENSTR,RCODE)
C
C    statements to fill           init_mss             
C
      CALL NCVINQ(NCID,13,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 130 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 130  CONTINUE
      CALL NCVGTC(NCID,13,START,COUNT,
     +          init_mss             ,LENSTR,RCODE)
C
C    statements to fill           tibnd_mss            
C
      CALL NCVINQ(NCID,14,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 140 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 140  CONTINUE
      CALL NCVGTC(NCID,14,START,COUNT,
     +          tibnd_mss            ,LENSTR,RCODE)
C
C    statements to fill           sst_mss              
C
      CALL NCVINQ(NCID,15,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 150 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 150  CONTINUE
      CALL NCVGTC(NCID,15,START,COUNT,
     +          sst_mss              ,LENSTR,RCODE)
C
C    statements to fill           ozone_mss            
C
      CALL NCVINQ(NCID,16,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 160 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 160  CONTINUE
      CALL NCVGTC(NCID,16,START,COUNT,
     +          ozone_mss            ,LENSTR,RCODE)
C
C    statements to fill           hyai                 
C
      CALL NCVINQ(NCID,17,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 170 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 170  CONTINUE
      CALL NCVGT(NCID,17,START,COUNT,
     +          hyai                 ,RCODE)
C
C    statements to fill           hybi                 
C
      CALL NCVINQ(NCID,18,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 180 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 180  CONTINUE
      CALL NCVGT(NCID,18,START,COUNT,
     +          hybi                 ,RCODE)
C
C    statements to fill           hyam                 
C
      CALL NCVINQ(NCID,19,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 190 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 190  CONTINUE
      CALL NCVGT(NCID,19,START,COUNT,
     +          hyam                 ,RCODE)
C
C    statements to fill           hybm                 
C
      CALL NCVINQ(NCID,20,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 200 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 200  CONTINUE
      CALL NCVGT(NCID,20,START,COUNT,
     +          hybm                 ,RCODE)
C
C    statements to fill           lev                  
C
      CALL NCVINQ(NCID,21,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 210 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 210  CONTINUE
      CALL NCVGT(NCID,21,START,COUNT,
     +          lev                  ,RCODE)
C
C    statements to fill           ilev                 
C
      CALL NCVINQ(NCID,22,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 220 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 220  CONTINUE
      CALL NCVGT(NCID,22,START,COUNT,
     +          ilev                 ,RCODE)
C
C    statements to fill           lat                  
C
      CALL NCVINQ(NCID,23,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 230 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 230  CONTINUE
      CALL NCVGT(NCID,23,START,COUNT,
     +          lat                  ,RCODE)
C
C    statements to fill           gw                   
C
      CALL NCVINQ(NCID,24,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 240 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 240  CONTINUE
      CALL NCVGT(NCID,24,START,COUNT,
     +          gw                   ,RCODE)
C
C    statements to fill           lon                  
C
      CALL NCVINQ(NCID,25,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 250 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 250  CONTINUE
      CALL NCVGT(NCID,25,START,COUNT,
     +          lon                  ,RCODE)
C
C    statements to fill           date                 
C
      CALL NCVINQ(NCID,26,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 260 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 260  CONTINUE
      CALL NCVGT(NCID,26,START,COUNT,
     +          date                 ,RCODE)
C
C    statements to fill           datesec              
C
      CALL NCVINQ(NCID,27,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 270 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 270  CONTINUE
      CALL NCVGT(NCID,27,START,COUNT,
     +          datesec              ,RCODE)
C
C    statements to fill           time                 
C
      CALL NCVINQ(NCID,28,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 280 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 280  CONTINUE
      CALL NCVGT(NCID,28,START,COUNT,
     +          time                 ,RCODE)
C
C    statements to fill           QR                   
C
      CALL NCVINQ(NCID,29,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 290 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 290  CONTINUE
      CALL NCVGT(NCID,29,START,COUNT,
     +          QR                   ,RCODE)
C
C    statements to fill           DTCOND               
C
      CALL NCVINQ(NCID,30,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 300 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 300  CONTINUE
      CALL NCVGT(NCID,30,START,COUNT,
     +          DTCOND               ,RCODE)
C
C    statements to fill           PS                   
C
      CALL NCVINQ(NCID,31,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 310 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 310  CONTINUE
      CALL NCVGT(NCID,31,START,COUNT,
     +          PS                   ,RCODE)
C
C    statements to fill           PRES                 
C
      CALL NCVINQ(NCID,32,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 320 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 320  CONTINUE
      CALL NCVGT(NCID,32,START,COUNT,
     +          PRES                 ,RCODE)
C
C     HERE IS WHERE YOU WRITE STATEMENTS TO USE THE DATA
C
c++csz
c     Initialize default values
      dbg_lvl=2
      in_fl='/usr/tmp/zender/spcp_85_8589_amip5_8589_01_r15.nc'
      out_fl='/usr/tmp/zender/r15.cos'
      nbr_lat=max_nbr_lat
      nbr_lev=max_nbr_lev
      nbr_lon=max_nbr_lon
      nbr_time=max_nbr_time
      open (out_fl_unit,file=out_fl,form='unformatted')

c     precompute sigma levels of 18 layer grid
      if (dbg_lvl.gt.0) write(6,'(a)') 'Computing sigma grid...'
      do idx_lon=1,nbr_lon
         do idx_lev=1,nbr_lev
            do idx_lat=1,nbr_lat
               do idx_time=1,nbr_time
                  sigma18(idx_lon,idx_lev,idx_lat,idx_time)=
     $                 pres(idx_lon,idx_lev,idx_lat,idx_time)/
     $                 ps(idx_lon,idx_lat,idx_time)
               end do           ! end loop over time
            end do              ! end loop over lat
         end do                 ! end loop over lev
      end do                    ! end loop over lon

c     vertical interpolation
      if (dbg_lvl.gt.0) write(6,'(a)') 'Vertical interpolation...'
      do idx_lon=1,nbr_lon
         do idx_lev9=1,9
            do idx_lat=1,nbr_lat
               do idx_time=1,nbr_time
c     bracket the 9 layer sigma level by two 180 layer hybrid levels
c     initialize the current top neighbor (ctn)
                  idx_ctn=1
                  dst_ctn=sigma9(idx_lev9)-
     $                 sigma18(idx_lon,idx_ctn,idx_lat,idx_time)
                  if (dst_ctn.lt.0.) stop 'Bad 1st top level'
c     find the 18 layer sigma level closest to the current 9 layer sigma level
                  do idx_lev=2,nbr_lev
                     dst=sigma9(idx_lev9)-
     $                    sigma18(idx_lon,idx_lev,idx_lat,idx_time) 
c     if candidate hybrid level is above 9 layer sigma level and closer than ctn
                     if ((dst.ge.0.).and.(dst.lt.dst_ctn)) then 
                        idx_ctn=idx_lev
                        dst_ctn=dst
                     endif
                  end do        ! end loop over lev

c     make sure everything is ok
                  if (idx_ctn.eq.18) stop 'ctn is lev 18'
                  if (
     $                 (sigma9(idx_lev9).lt.
     $                 sigma18(idx_lon,idx_ctn,idx_lat,idx_time))
     $                 .or.
     $                 (sigma9(idx_lev9).gt.
     $                 sigma18(idx_lon,idx_ctn+1,idx_lat,idx_time))
     $                 ) stop 'Bracketing failed'

c     the 9 layer sigma level is bracketed
c     find the fractional distance to the ctn and cbn
                  dst_frc_ctn=
     $                 (sigma9(idx_lev9)-
     $                 sigma18(idx_lon,idx_ctn,idx_lat,idx_time))/
     $                 (sigma18(idx_lon,idx_ctn+1,idx_lat,idx_time)-
     $                 sigma18(idx_lon,idx_ctn,idx_lat,idx_time))

                  if ((dst_frc_ctn.lt.0.).or.
     $                 (dst_frc_ctn.gt.1.)) 
     $                 stop 'dst_frc_ctn out of bounds'

c     assign the 9 level sigma values
                  dtcond9(idx_lon,idx_lev9,idx_lat,idx_time)=
     $                 (1-dst_frc_ctn)*
     $                 dtcond(idx_lon,idx_ctn,idx_lat,idx_time)+
     $                 dst_frc_ctn*
     $                 dtcond(idx_lon,idx_ctn+1,idx_lat,idx_time)

c     assign the 9 level sigma values
                  qr9(idx_lon,idx_lev9,idx_lat,idx_time)=
     $                 (1-dst_frc_ctn)*
     $                 qr(idx_lon,idx_ctn,idx_lat,idx_time)+
     $                 dst_frc_ctn*
     $                 qr(idx_lon,idx_ctn+1,idx_lat,idx_time)

               end do           ! end loop over time
            end do              ! end loop over lat
         end do                 ! end loop over lev9
      end do                    ! end loop over lon

c     output heating perturbations and switch longitudes 
      if (dbg_lvl.gt.0) write(6,'(a)') 'Output to binary file...'
      do idx_var=1,2
         do idx_lev9=1,9
            do idx_lon=1,nbr_lon
c     longitude indices of output array start at 180 W
               idx_lon9=mod(idx_lon+24,nbr_lon)
               if (idx_lon9.eq.0) idx_lon9=48
               do idx_lat=1,nbr_lat
                  if (dbg_lvl.gt.1) write(6,73)
     $                 'lon9(',idx_lon9,') = ',lon(idx_lon),', ',
     $                 'lev(',idx_lev9,') = ',sigma9(idx_lev9),', ',
     $                 'lat(',idx_lat,') = ',lat(idx_lat),', ',
     $         'dtcond9 = ',dtcond9(idx_lon,idx_lev9,idx_lat,1),', ',
     $                 'qr9 = ',qr9(idx_lon,idx_lev9,idx_lat,1),', ',
     $                 ' '
 73               format(
     $                 a,i2,a,f7.3,a,
     $                 a,i1,a,f4.3,a,
     $                 a,i2,a,f7.3,a,
     $                 a,1p,e9.2,0p,a,
     $                 a,1p,e9.2,0p,a,
     $                 a) 
                  
c     first slew of records is condensation perturbation, second is radiative perturbation
                  if (idx_var.eq.1) then 
                     record(idx_lon9,idx_lat)=
     $                    dtcond9(idx_lon,idx_lev9,idx_lat,1)
                  else
                     record(idx_lon9,idx_lat)=
     $                    qr9(idx_lon,idx_lev9,idx_lat,1)
                  endif
                  
               end do           ! end loop over lat
            end do              ! end loop over lon
            write(out_fl_unit) record
         end do                 ! end loop over lev9
      end do                    ! end loop over var
      close (out_fl_unit)
c--csz
C
C
      STOP
      END
