      PROGRAM LOOKUPTABLE
C-----------------------------------------------------------------------------------
      IMPLICIT NONE
C----------------------------------------------------------------------------------------
c     VARIABLES USED TO DESCRIBE SOIL SIZE DISTRIBUTION
      INTEGER, PARAMETER            :: max_soil_type=4           !Number of soil types
      DOUBLE PRECISION              :: mmd_soil(max_soil_type)   !mass median diameter 
      DOUBLE PRECISION              :: smd_soil                  !surface median diameter
      DOUBLE PRECISION              :: nmd_soil                  !number median diameter
      DOUBLE PRECISION              :: sigma(max_soil_type)      !sigma of soils
      DOUBLE PRECISION, PARAMETER   :: total_surface=1.d0        !Total surface
      DOUBLE PRECISION, PARAMETER   :: pi=3.141592654d0          !pi
      DOUBLE PRECISION, PARAMETER   :: cst_slt=2.61d0            !Constant in horizontal flux 
      DOUBLE PRECISION, PARAMETER   :: dns_mdp=1.2d0             !Air density 
      DOUBLE PRECISION, PARAMETER   :: grv_sfc=9.81d0            !gravity
      INTEGER                       :: Nread,N                   !Counting
      INTEGER, PARAMETER            :: max_sizes=10000             !Number of size bins along soil
      INTEGER, PARAMETER            :: max_ustr=100              !Number of friction wind speeds
      Character                     :: trash1,trash2,trash3      !read in an throw away
      DOUBLE PRECISION              :: sz(max_sizes)             !Sizes of soil
      DOUBLE PRECISION              :: logstep                   !value of dlnD
      DOUBLE PRECISION              :: wnd_frc_thr(max_sizes)    !friction wind thr
      DOUBLE PRECISION              :: wnd_frc(max_ustr)         !Friction wind
      DOUBLE PRECISION              :: wnd_frc_rat               !Ratio of wind/thr wind
      DOUBLE PRECISION              :: dS(max_sizes)             !value of dS for a given D
      DOUBLE PRECISION              :: 
     &     size_segr_hrz_flx(max_ustr,max_sizes)                 !Size segregated Horizontal flux
      DOUBLE PRECISION              ::
     &     fraction(max_ustr,max_sizes)                          !Fraction of horizontal flux
      DOUBLE PRECISION              :: total                     !Checking variable
      INTEGER                       :: wnd_frc_idx, sz_idx       !Counting variables
c
c     VARIABLES USED BY ZENDER DUST MODULE
      INTEGER, PARAMETER            :: dst_src_nbr=3                !Number of source modes
      INTEGER                       :: src_idx                      !Counting var. for source modes
      DOUBLE PRECISION              :: dmt_vma_src(dst_src_nbr)     ![m] MMD of source modes
      DOUBLE PRECISION              :: sgm_src(dst_src_nbr)         ![-] standard deviation of source modes
      DOUBLE PRECISION              :: dmt_vmm_src(dst_src_nbr)     ![m] MASS MEAN source modes
      DOUBLE PRECISION              :: egy_bnd_src(dst_src_nbr)     ![J] Binding energy for source modes
      DOUBLE PRECISION              :: egy_frc_src(dst_src_nbr)     ![-] Energy fraction of source modes
      DOUBLE PRECISION              :: egy_frc_src_tmp(dst_src_nbr) ![-] Energy fraction of source modes
      DOUBLE PRECISION              :: mss_frc_src(dst_src_nbr)     ![-] Mass fraction of source modes
      DOUBLE PRECISION              :: mss_frc_src_tmp(dst_src_nbr) ![-] Mass fraction of source modes for iteration
      DOUBLE PRECISION              :: nbr_frc_src(dst_src_nbr)     ![-] number fraction of source modes
      DOUBLE PRECISION              :: nbr_frc_src_tmp(dst_src_nbr) ![-] Number fraction of source modes for one iteration
      DOUBLE PRECISION              :: dst_slt_flx_rat_ttl          ![m-1] Alpha (vert/hor flux)
      DOUBLE PRECISION              :: dst_slt_flx_rat_ttl_tmp      ![m-1] Alpha for one iteration
      DOUBLE PRECISION              :: flx_mss_vrt_dst_ttl          ![kgm-2s-1] Vertical mass flux 
      DOUBLE PRECISION              :: flx_mss_vrt_dst_ttl_tmp      ![kgm-2s-1] Vertical mass flux for one iteration
      DOUBLE PRECISION              :: flx_nbr_vrt_dst_ttl          ![#m-2s-1] Vertical number flux
      DOUBLE PRECISION              :: flx_nbr_vrt_dst_ttl_tmp      ![#m-2s-1] Vertical number flux for one iteration
      DOUBLE PRECISION              :: flx_mss_hrz_slt_ttl          ![kgm-1s-1] Horizontal flux
      DOUBLE PRECISION              :: flx_mss_hrz_slt_ttl_tmp      ![kgm-1s-1] Horizontal flux for one iteration
      DOUBLE PRECISION,PARAMETER    :: dns_soil=2650.d0             !Density of soil      
      DOUBLE PRECISION              :: dmt_slt_bin                  ![m] Diameter of saltating bin 
c
      LOGICAL                       :: flg_mbl             !Mobilization flag
      INTEGER                       :: soil_type_idx
      INTEGER, PARAMETER            :: fileno=21           !file number of outfile
c
c     OUTPUT TO FILE
      DOUBLE PRECISION              :: 
     &     alphavalue(max_ustr,max_soil_type)                     ![m-1] Value of alpha for a given u* and soiltype
      DOUBLE PRECISION              ::
     &     alphanumbervalue(max_ustr,max_soil_type)               ![#/kg/m] Number based alpha for a given u* and soiltype
      DOUBLE PRECISION              :: 
     &     mss_frc_src_value(max_ustr,max_soil_type,dst_src_nbr)  ![-] Value of mss_frc_src for given u* and soiltype
      CHARACTER*50                   :: soilname(max_soil_type)   ![-] name of soil types
      CHARACTER*50                   :: infilename                ![-] name of file to read in
c
      DATA dmt_vma_src /1.5d-6, 6.7d-6, 14.2d-6 /
      DATA sgm_src     /1.7d0, 1.6d0, 1.5d0 /
      DATA egy_bnd_src /3.61d-7, 3.52d-7, 3.46d-7 /
      DATA mmd_soil    /125.d-6, 210d-6, 520d-6, 690d-6/
      DATA sigma       /1.8, 1.6, 1.6, 1.6 /                      !Corrected, now OK with MBA97
      DATA soilname    /'ASS ', 'FS', 'SS' , 'CS' /
      flg_mbl = .TRUE.
c
c--calculate mass mean diameters of source mode
      do src_idx=1,dst_src_nbr
         dmt_vmm_src(src_idx)=dmt_vma_src(src_idx)
     &        *exp(log(sgm_src(src_idx))*log(sgm_src(src_idx))*(-1.5d0))
      enddo
      write(6,*)'median',dmt_vma_src
      write(6,*)'mean ',dmt_vmm_src

c      stop
c---read in values from Charlie's mie program
      if(max_sizes.eq.100)then
         infilename='input_100.txt'
      elseif(max_sizes.eq.1000)then
         infilename='input_1000.txt'
      elseif(max_sizes.eq.10000)then
         infilename='input_10000.txt'
      else
         write(6,*)'check if you have input for'
         write(6,*)max_sizes ,'soil sizes'
         stop
      endif
      open(20,file=infilename)
      do N=1,max_sizes
         read(20,*)trash1,sz(N),trash2,wnd_frc_thr(N),trash3
c         write(6,*)N,trash1,sz(N),trash2,wnd_frc_thr(N),trash3
      enddo
c
c     Conversion to from radius to diameter
      sz(:) = sz(:)*2.d0      

c     Setting array of wind frictions
c     to values from 0-1 ms-1
      do wnd_frc_idx=1,max_ustr
         wnd_frc(wnd_frc_idx)=float(wnd_frc_idx)*1.d0/100.d0
      enddo
c
c---Calculating dlnD for use in size distribution function
      logstep=
     &     (log(sz(max_sizes))-log(sz(1)))
     &     /float(max_sizes)
c
c*****************************************************************
c***STARTING LOOP FOR A SOIL TYPE
c****************************************************************
c
      do soil_type_idx=1,max_soil_type
c
c****************************************************************
c     FIRST PART OF PROGRAM IS TO GENERATE THE ARRAY CALLED 
c     "FRACTION" WHICH IS THE FRACTION OF THE HORIZONTAL (SALTATING)
c     FLUX WITH A CERTAIN SIZE
c*****************************************************************
c
c     Using formulas7.50 / 7.52 in Seinfeld
c     Number median diameter of soil :
         nmd_soil = exp (log(mmd_soil(soil_type_idx))-
     &        3.d0*log(sigma(soil_type_idx))*log(sigma(soil_type_idx)))
c     Surface median diameter of soil :
         smd_soil = exp (log(nmd_soil)+
     &        2.d0*log(sigma(soil_type_idx))*log(sigma(soil_type_idx)))
c     
c     Calculating value for dS/dlnD for the classes we have
c     Calculating formula 7.33 in Sinfeld for surface
c     In the and, multiplying by logstep to get dS for a given size
         do N=1,max_sizes
            dS(N) = total_surface
     &           /(sqrt(2.d0*pi)*log(sigma(soil_type_idx)))
     &           *exp(-1.d0*((log(sz(N))-log(smd_soil))**2
     &           /(2.d0*log(sigma(soil_type_idx))
     &           *log(sigma(soil_type_idx)))))
     &           *logstep
         enddo
c     
c     copy from aerosol model (zender dust module) program
c     fraction(wnd_frc_idx,sz_idx) is an array containing how much of the 
c     horizontal (saltating) flux has that specific size
c     
         total=0.d0
         do wnd_frc_idx=1,max_ustr
c     initializing "total" for all wind friction speeds
            total=0.d0
            do sz_idx=1,max_sizes
c     
               wnd_frc_rat=wnd_frc_thr(sz_idx)/wnd_frc(wnd_frc_idx) ! [frc]
c---  Horizontal flux for each size:      
               size_segr_hrz_flx(wnd_frc_idx,sz_idx)= ! [kg m-1 s-1] 
     &              max(0.d0,
     $              cst_slt*dns_mdp*(wnd_frc(wnd_frc_idx)**3.0)*
     $              (1.0-wnd_frc_rat)*(1.0+wnd_frc_rat)*
     &              (1.0+wnd_frc_rat)
     &              /grv_sfc 
     &              *dS(sz_idx) )
c     
c     Summing up total flux for that size:
               total=total + size_segr_hrz_flx(wnd_frc_idx,sz_idx)
            enddo               !loop on size
            if(total.gt.0.d0)then
               fraction(wnd_frc_idx,:) = 
     &              size_segr_hrz_flx(wnd_frc_idx,:)
     &              /total
            else
               fraction(wnd_frc_idx,:)=0.d0
            endif
         enddo                  !loop on wind friction
c     
c     Check if we have a total of 1 for all winds 
         do wnd_frc_idx=1,max_ustr
            total=0.d0
            do sz_idx=1,max_sizes
               total=total+fraction(wnd_frc_idx,sz_idx)
            enddo
         enddo
c     
c*******************************************************************************
c     SECOND PART OF PROGRAM IS TO USE THE ARRAY FRACTION TO CALCULATE 
c     THE MASS FRACTION OF EACH SOURCE MODE AND THE VALUE OF "ALPHA" 
c     AND THEN PRINT THIS TO A LOOK UP TABLE TO USE IN THE DUST MODULE
c******************************************************************************
c     
         flx_mss_hrz_slt_ttl = 1.d0
c     
         do wnd_frc_idx=1,max_ustr
            flx_mss_vrt_dst_ttl= 0.d0
            flx_nbr_vrt_dst_ttl = 0.d0
            nbr_frc_src(:)=0.d0
            mss_frc_src(:)=0.d0
c
            write(6,*)'calculation for ustr',wnd_frc(wnd_frc_idx)
c     
            do sz_idx=1,max_sizes
c     
c     Diameter of the saltating size
               dmt_slt_bin=sz(sz_idx)
c     Set total horizontal flux to 1.
c     Then horizontal flux in loop is equal to "fraction"
               flx_mss_hrz_slt_ttl_tmp =
     &              flx_mss_hrz_slt_ttl
     &              *fraction(wnd_frc_idx,sz_idx)
c     
c     Calculate the mass fraction of each source mode
c     using the AlG01 formulation
c     
               call mss_frc_src_AlG01_get(
     &              wnd_frc(wnd_frc_idx), ! I [m/s] wind friction velocity
     &              egy_bnd_src,          ! I [J] binding energy of source aerosols
     &              mss_frc_src_tmp,      ! O [frc] Mass fraction of each source mode for specific diameter
     &              dmt_vmm_src,          ! I  [m] mass MEAN diameter of source distributions  
c     &              dmt_vma_src,          ! I [m] mass median diamter of sourc distributions
     &              egy_frc_src_tmp,      ! O [-] Energy fraction used on each mode (p(i) in AlG01) for specific diameter
     &              dmt_slt_bin,          ! I [m] Diameter of saltating bin
     &              dns_soil,             ! [kgm-3] density of soil
     &              dst_src_nbr           ! I [-] number of source modes
     &              )
c     
c     Calculate alpha and total vertical flux
c     Done for a all longitudes
               call flx_mss_vrt_dst_ttl_AlG01_get(
     &              dst_slt_flx_rat_ttl_tmp,  ! O [frc] ratio of vertical to horizontal dust flux
     &              flg_mbl,                  ! I [flg] mobilization flag
     &              mss_frc_src_tmp,          ! I [frc] mass fraction source
     &              flx_mss_hrz_slt_ttl_tmp,  ! I [kgms-1] total horizontal flux of soil
     &              flx_mss_vrt_dst_ttl_tmp,  ! O [kgm-2s-1] total vertical dust flux
     &              dmt_vmm_src,              ! I [m] mass MEAN diameter of source
c     &              dmt_vma_src,              ! I [m] mass median diamter of source
     &              egy_bnd_src,              ! I [J] bounding energy of source dust
     &              egy_frc_src_tmp,          ! I [-] Energy fraction use in each mode (p(i) in AlG01)
     &              dns_soil,                 ! I [kgm-3] density of soil
     &              dst_src_nbr,              ! I [-] number of source modes
     &              flx_nbr_vrt_dst_ttl_tmp   ! O [#m-2s-1] number flux of dust
     &              )                          
                                            
c     summing up the total dust flux
               flx_mss_vrt_dst_ttl = flx_mss_vrt_dst_ttl 
     &              + flx_mss_vrt_dst_ttl_tmp     

c     Summing up total number flux
               flx_nbr_vrt_dst_ttl = 
     &              flx_nbr_vrt_dst_ttl 
     &              + flx_nbr_vrt_dst_ttl_tmp
c
c     summing up mass fraction of source weighted by vertical dust flux FSD
               mss_frc_src(:) = 
     &              mss_frc_src_tmp(:)*flx_mss_vrt_dst_ttl_tmp
     &              + mss_frc_src(:)
c     
            enddo               !sz_idx
C
            mss_frc_src(:)=mss_frc_src(:)
     &           / max(1.d-40,flx_mss_vrt_dst_ttl)
c
            if(flx_mss_vrt_dst_ttl.gt.0.d0)then
               alphavalue(wnd_frc_idx,soil_type_idx) = 
     &              flx_mss_vrt_dst_ttl 
     &              /flx_mss_hrz_slt_ttl
c
               alphanumbervalue(wnd_frc_idx,soil_type_idx) = 
     &              flx_nbr_vrt_dst_ttl
     &              /flx_mss_hrz_slt_ttl
c     
               mss_frc_src_value(wnd_frc_idx,soil_type_idx,:)
     &              = mss_frc_src(:)         
c     
            else
               alphavalue(wnd_frc_idx,soil_type_idx) = 0.d0
               mss_frc_src_value(wnd_frc_idx,soil_type_idx,:)= 0.d0
            endif
c
            write(6,*)'total mass flux ',flx_mss_vrt_dst_ttl
            write(6,*)'total nmbr flux ',flx_nbr_vrt_dst_ttl
         enddo                  !ustar
c     
      enddo                     !SOILTYPE
c*****************************************************************
c     PART 3 OF PROGRAM IS TO WRITE THE RESULTS TO A NICE FILE
c     WHICH CAN BE USED BY THE ZENDER DUST MODULE
c*****************************************************************
      open(fileno,file='AlG01data.dat',form='FORMATTED')
      do soil_type_idx=1,max_soil_type
         write(fileno,1001)SOILNAME(soil_type_idx)
         write(fileno,1002)'ustr','alpha','numberalpha',
     &        'massfrac(1)','massfrac(2)',
     &        'massfrac(3)'
         do wnd_frc_idx=1,max_ustr
            write(fileno,1003)wnd_frc(wnd_frc_idx),
     &           alphavalue(wnd_frc_idx,soil_type_idx),
     &           alphanumbervalue(wnd_frc_idx,soil_type_idx),
     &        (mss_frc_src_value(wnd_frc_idx,soil_type_idx,src_idx)
     &           ,src_idx=1,dst_src_nbr)
         enddo  !loop on wnd_frc
      enddo     !loop on soil_type
      close(fileno)
c
 1001 format(a50)
 1002 format(6a15)
 1003 format(6e15.6)
c
c---output to matlab
c
      open(fileno, file='formatlab.m', form='FORMATTED')
      write(fileno,1010)'ustar',1,'=[',
     &     (wnd_frc(wnd_frc_idx),wnd_frc_idx=1,max_ustr),']'

      do soil_type_idx=2,2 !max_soil_type
         write(fileno,1010)'alpha',soil_type_idx,' = [',
     &        (alphavalue(wnd_frc_idx,soil_type_idx)
     &        ,wnd_frc_idx=1,max_ustr),']'
      enddo
c
c      do soil_type_idx=2,2 !max_soil_type
c         write(fileno,1011)'numberalpha',soil_type_idx,' = [',
c     &        (alphanumbervalue(wnd_frc_idx,soil_type_idx)
c     &        ,wnd_frc_idx=1,max_ustr),']'
c      enddo
 1010 format(a5,i1,a4,100e15.6,a1)
c 1011 format(a11,i1,a4,100e15.6,a1)

      open(52,file='dstsltsbl.F90',form='FORMATTED')
      do soil_type_idx=1,max_soil_type
         write(52,1022)'      ! SOIL TYPE ',soil_type_idx,
     &        'MMD=',mmd_soil(soil_type_idx)*1e6,'um',
     &        'sigma=',sigma(soil_type_idx)
         write(52,'(a6,a50)')
     &        ' ','! First index is wind friction speed in cm/s'
         write(52,'(a6,a50)')' ','! Second index is soil type '
         write(52,'(a6,a50)')
     &        ' ','! Third index is number of source mode'
      do wnd_frc_idx=1,max_ustr
       write(52,1021)
     &        '      ',
     &        'data ',
     &        '(mss_frc_src_lut(',wnd_frc_idx,',',soil_type_idx,
     &        ',src_idx),src_idx=1,dst_src_nbr)) /'
     &        ,mss_frc_src_value(wnd_frc_idx,soil_type_idx,1),','
     &        ,mss_frc_src_value(wnd_frc_idx,soil_type_idx,2),','
     &        ,mss_frc_src_value(wnd_frc_idx,soil_type_idx,3),
     &        '/'
      enddo
      enddo
 1021 format(a6,a5,a17,i3,a1,i1,a35,3(e15.6,a1))
 1022 format(a20,i4,a6,f8.1,a4,a8,f8.1)
c
      write(52,'(a6,a50)')' ','!Alpha values, vert/hor flux [m-1]'
c
      do soil_type_idx=1,max_soil_type
         write(52,'(a6,a50)')' ','! First index is u* in cm/s'
         write(52,'(a6,a50)')' ','! Second index soil type (see above)'
      do wnd_frc_idx=1,max_ustr
         write(52,1023)' ','data ','dst_slt_flx_rat_ttl_lut ('
     &        ,wnd_frc_idx,',',soil_type_idx,
     &        ') /'
     &        ,alphavalue(wnd_frc_idx,soil_type_idx),
     &        '/'
      enddo
      enddo
 1023    format(a6,a5,a25,i3,a1,i1,a3,e15.6,a1)
c
      end program lookuptable  

c********************************************************
c********************************************************
C***                      SUBROUTINES                 ***
c********************************************************

c********************************************************************
      subroutine flx_mss_vrt_dst_ttl_AlG01_get(
     &     dst_slt_flx_rat_ttl,    ! O [-] Ratio of vert/hor. flux
     &     flg_mbl,                ! I [flg] Mobilization candidate flag
     &     mss_frc_src,            ! I [-] mass fraction in source
     &     flx_mss_hrz_slt_ttl,    ! I [kgm-1s-1] vertically integrated streamwize flux
     &     flx_mss_vrt_dst_ttl,    ! O [kgm-2s-1] total vertical dust flux
     &     dmt_vma_src,            ! I [m] median diameters of source 
     &     egy_bnd_src,            ! I [J] bounding energy of source dust
     &     egy_frc_src,            ! I [-] Energy fraction used in each mode (p(i) in AlG01)           
     &     dns_soil,               ! I [kgm-3] Density of soil
     &     dst_src_nbr,            ! I [-] number of source modes
     &     flx_nbr_vrt_dst_ttl     ! O [# m-2 s-1] Number flux of dust
     &     )
c
c---Purpose:
c---Calculate alpha (vertflux/horflux) and total vertical flux of dust 
c---for a given friction velocity,
c---density and size distribution of soil and source aerosol.
c---Formulation is given by Alfaro and Gomes (2001) (eq.6).
c
c---Will be called from dstmbl. Will replace old routine called
c---flx_mss_vrt_dst_ttl_MaB95_get where alpha is calculated as a function
c---of clay content.
c
c---Coded by Alf Grini, spring 2002, UCI/UiO
c------------------------------------------------------------------
      IMPLICIT NONE
c-------------------------------------------------------
c
c---OUTPUT
      double precision, intent(out) :: dst_slt_flx_rat_ttl      ![1/m] vert / hor. flux
      double precision, intent(out) :: flx_mss_vrt_dst_ttl      ![kgm-2s-1] total vertical mass flux
      double precision, intent(out) :: flx_nbr_vrt_dst_ttl      ![#m-2s-1] total vertical number flux
c---INPUT
      integer, intent(in)  :: dst_src_nbr                      ! [-] number of source modes
      logical, intent(in)  :: flg_mbl                          ! [flg] Mobilization candidate (?)
      double precision, intent(in) :: mss_frc_src(dst_src_nbr) ! [-] Mass fraction source (p(i) in A/G))
      double precision, intent(in) :: flx_mss_hrz_slt_ttl      ! [kgm-1s-1] vertically integrated streamwize flux
      double precision, intent(in) :: dmt_vma_src(dst_src_nbr) ! [m] mass median diameter of source
      double precision, intent(in) :: egy_bnd_src(dst_src_nbr) ! [J] binding energy of source modes
      double precision, intent(in) :: egy_frc_src(dst_src_nbr) ! [-] Energy fraction used in each source mode (p(i) in AlG01)
      double precision, intent(in) :: dns_soil                 ! [kgm-3] density of soil 
c---LOCAL IN THIS ROUTINE
      integer               :: lon_idx, src_idx                ! counting variables
      double precision, parameter   :: pi = 3.141592654d0      ! pi
      double precision, parameter   :: pisixths = pi/6.d0      ! pi/6
      double precision, parameter   :: beta = 16300.d-2        ! ms-2 (given AlG01-article)
c
c---Begin code
c
c---initializing alpha
      dst_slt_flx_rat_ttl=0.d0
c---initializing number flux
      flx_nbr_vrt_dst_ttl=0.d0

      if(flg_mbl)then
c     
         do src_idx = 1, dst_src_nbr !Loop on source distributions
c     
c---solving eq. 6 in Alfaro/Gomes.
c     
c---finding alpha (= vertical flux/horizontal flux)
c---in the limit of very high kinetic energy of soil aggregates, 
c---the mss_frc_src will be close to 1 for the smallest mode,
c---and 0 for the two largest ones 
c---Then alpha will reach a limit value  
c
c---alpha is sum of contribution from 3 source distributions:
c---solving eq 6 in AlG01 :
c
            dst_slt_flx_rat_ttl  = 
     &           dst_slt_flx_rat_ttl  +
     &           pisixths * dns_soil * beta 
     &           * egy_frc_src(src_idx)
     &           * dmt_vma_src(src_idx)**3
     &           / egy_bnd_src(src_idx)
c
c     
         enddo                  !loop on source distributions
c
c---getting number flux of each mode in #/m2/s
c---eqn 3 in AlG01
         do src_idx=1,dst_src_nbr
            flx_nbr_vrt_dst_ttl =
     &           beta 
     &           * flx_mss_hrz_slt_ttl 
     &           * egy_frc_src(src_idx) 
     &           /egy_bnd_src(src_idx)
     &           +flx_nbr_vrt_dst_ttl
         enddo
c
c---  finding the total vertical dust flux from alpha and horizontal
c---  for each longitude:
c     
         flx_mss_vrt_dst_ttl  = 
     &        dst_slt_flx_rat_ttl 
     &        * flx_mss_hrz_slt_ttl
c     
c     
      endif                     !Check on mobilization flag
c     
      return
      end subroutine flx_mss_vrt_dst_ttl_AlG01_get
c********************************************************************
c********************************************************************
c**************************************************************************
      subroutine mss_frc_src_AlG01_get(
     &     wnd_frc,             ! I [m/s] wind friction velocity
     &     egy_bnd_src,         ! I [kgm2s-2] binding energy of sources
     &     mss_frc_src,         ! O [-] Mass fraction of each of the three src dst.)
     &     dmt_vma_src,         ! I [m] Mass median diamter of source functions
     &     egy_frc_src,         ! O [-] Fraction of energy used in each source mode (p(i) in AlG01)
     &     dmt_slt_bin,         ! I [m] Diameter of the saltating soil
     &     dns_soil,            ! I [kgm-3] density of soil
     &     dst_src_nbr          ! I [-] number of source modes
     &     )
c
c---Purpose:
c---Determine based on friction wind speed how much energy is available
c---for saltation. Based on this energy, we can determine the mass fraction
c---of each of the three distributions proposed by Alfaro and Gomes.
c
c---Called from dstmbl
c
c---Code by Alf Grini, spring 2002, UCI/UiO
c---------------------------------------------------------------
      implicit none
c----------------------------------------------------------------
c
c---INPUT
      integer         , intent(in)        :: dst_src_nbr              ![-] Number of source modes
      double precision, intent(in)        :: wnd_frc                  ![m/s] Wind friction vel 
      double precision, intent(in)        :: egy_bnd_src(dst_src_nbr) ![J] binding energy of soil
      double precision, intent(in)        :: dmt_vma_src(dst_src_nbr) ![m] mass median diameter of source functions
      double precision, intent(in)        :: dmt_slt_bin              ![m] diameter of saltating bin
      double precision, intent(in)        :: dns_soil                 ! [kgm-3] density of soil
c---OUTPUT
      double precision, intent(out) :: mss_frc_src(dst_src_nbr)  ![-] fraction of mass in each or three source distr.
      double precision, intent(out) :: egy_frc_src(dst_src_nbr)  ![-] Fraction of energy used in three source distr.
c---LOCAL IN ROUTINE
      double precision            :: egy_kin_agr             ! [kgm2s-2] Kinetic energy of soil aggregates
      double precision, parameter :: pi = 3.141592654d0      ! [-] The number pi
      double precision, parameter :: onetwelwth = 1.d0/12.d0 ! [-] One divided by twelve
      double precision            :: egy_mss_src_ttl         ! [m3J-1] Total energy weighted by volume of modes
      integer                     :: lon_idx                 ! [-] counting variable for longitude
      integer                     :: src_idx                 ! [-] counting variable for source distributions
      double precision            :: nbr_frc_src(dst_src_nbr)! [-] number fraction of source
c
c---THE NEW THING HERE (COMPARED TO She84) , IS THAT MSS_FRC_SRC 
c---IS CALCULATED IN EACH TIME STEP. 
c
c---initializing
      mss_frc_src(:)=0.d0     !No mass if egy_kin_agr is too small
      egy_frc_src(:)=0.d0
      nbr_frc_src(:)=0.d0
c
c      write(6,*)'starting loop on longitude'
c     
c---  kinetic energy of aggregates : (Eq. 1 in AlG01)
      egy_kin_agr =
     &     dns_soil * pi * onetwelwth
     &     *dmt_slt_bin*dmt_slt_bin*dmt_slt_bin
     &     *(20.d0*wnd_frc)
     &     *(20.d0*wnd_frc) !units [kgm2s-2 = J]
c
c         write(6,*)'dmt_slt_bin', dmt_slt_bin*1.d6,'um'
c         write(6,*)'egy_kin_agr  ',egy_kin_agr,'J'
c         write(6,*)'wnd_frc',wnd_frc,'m/s'
c         write(6,*)'egy_bnd_src',egy_bnd_src
c         write(6,*)'dns_soil ',dns_soil,'kg/m3'
c         write(6,*)'onetwelwth*12',onetwelwth*12.d0
c         write(6,*)'pi ',pi
c
c---This if-test is basically the sum of table 2. in AlG01 :
c
         if(egy_kin_agr.gt.egy_bnd_src(3)
     &        .and.egy_kin_agr.lt.egy_bnd_src(2))then
c---Only largest mode available
            egy_frc_src(3)=1.d0
c            write(6,*)'only largest mode'
c
         elseif(egy_kin_agr.gt.egy_bnd_src(2)
     &           .and.egy_kin_agr.lt.egy_bnd_src(1))then
c---Only two largest modes available
c            write(6,*)'only largest two modes'
            egy_frc_src(2)= 
     &           (egy_kin_agr - egy_bnd_src(2))
     &           /(egy_kin_agr - egy_bnd_src(3))
c
            egy_frc_src(3)=
     &           1.d0 - egy_frc_src(2)

         elseif(egy_kin_agr.gt.egy_bnd_src(1))THEN
c---All three modes available
c            write(6,*)'all three modes'
            egy_frc_src(1)=
     &           (egy_kin_agr - egy_bnd_src(1))
     &           /(egy_kin_agr - egy_bnd_src(3))
            egy_frc_src(2)=
     &           (1.d0 - egy_frc_src(1))*
     &           (egy_kin_agr - egy_bnd_src(2))
     &           /(egy_kin_agr - egy_bnd_src(3))
            egy_frc_src(3)=
     &           1.d0 
     &           - egy_frc_src(2) 
     &           - egy_frc_src(1)
         endif
c
c---initializing
         egy_mss_src_ttl = 0.d0          
c
c---summing total kinetic energy used for saltiation weighted by volume
c---this is the denominator in eq: 9 AlG01
         do src_idx = 1, dst_src_nbr
            egy_mss_src_ttl = egy_mss_src_ttl 
     &           + egy_frc_src(src_idx)
     &           * dmt_vma_src(src_idx)**3
     &           /egy_bnd_src(src_idx)
         enddo
c
c---finding mass fraction of each mode by using eq. (9) in AlG01
         IF(egy_mss_src_ttl.gt.0.d0)THEN
            do src_idx = 1,dst_src_nbr          
               mss_frc_src(src_idx)=egy_frc_src(src_idx)
     &              *dmt_vma_src(src_idx)**3
     &              /egy_bnd_src(src_idx) 
     &              / egy_mss_src_ttl
            enddo               !Loop over source
c
c_test
            do src_idx=1,dst_src_nbr
               nbr_frc_src(src_idx)=
     &              (egy_frc_src(src_idx)/egy_bnd_src(src_idx))
     &              /
     &              ( (egy_frc_src(1)/egy_bnd_src(1))
     &              + (egy_frc_src(2)/egy_bnd_src(2))
     &              + (egy_frc_src(3)/egy_bnd_src(3)) )*100.d0

            enddo
c_test
c
         ELSE
            mss_frc_src(:)=0.d0   !No saltation
         ENDIF
c
c---If ekin lower than ekin_src(3) the if-loop will not be done
c---and we will end here. Mass fractions will be zero. 
c---NOTE!!!! :It has to be like this, because with 
c---with AlG01 formulation we can have saltation but no sandblasting.
c---( if u*>u*,tr, but ec<e(3) ) 
c

      return
      end subroutine mss_frc_src_AlG01_get
c*****************************************************************************


