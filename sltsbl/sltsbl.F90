! $Id$

! Purpose: Offline driver to generate saltation-sandblasting lookup-tables for
! Dust Entrainment And Deposition (DEAD) model

! Copyright (C) 2003--2005 Charlie Zender and Alf Grini
! This software is distributed under the terms of the GNU General Public License Version 2
! See http://www.gnu.ai.mit.edu/copyleft/gpl.html for full license text
! The original author of this software, Charlie Zender, wants to improve it
! with the help of your suggestions, improvements, bug-reports, and patches.
! Charlie Zender <zender at uci dot edu>
! Department of Earth System Science
! University of California at Irvine
! Irvine, CA 92697-3100

! Compilation:
! Debugging             : cd ~/sltsbl;make OPTS=D;cd - 
! Double precision reals: cd ~/sltsbl;make OPTS=D PRC=D;cd -
! Single precision reals: cd ~/sltsbl;make OPTS=D PRC=S;cd -

! Usage:
! sltsbl

program sltsbl
  use precision ! [mdl] Precision r8, i8, ...
  implicit none
  ! Variables used to describe soil size distribution
  integer,parameter::max_soil_type=4           !number of soil types
  real(r8)::mmd_soil(max_soil_type)   !mass median diameter 
  real(r8)::smd_soil                  !surface median diameter
  real(r8)::nmd_soil                  !number median diameter
  real(r8)::sigma(max_soil_type)      !sigma of soils
  real(r8),parameter::total_surface=1.0_r8        !total surface
  real(r8),parameter::pi=3.141592654d0          !pi
  real(r8),parameter::cst_slt=2.61d0            !constant in horizontal flux 
  real(r8),parameter::dns_mdp=1.2d0             !air density 
  real(r8),parameter::grv_sfc=9.81d0            !gravity
  integer::nread,n                   !counting
  integer,parameter::max_sizes=10000             !number of size bins along soil
  integer,parameter::max_ustr=100              !number of friction wind speeds
  character::trash1,trash2,trash3      !read in an throw away
  real(r8)::sz(max_sizes)             !sizes of soil
  real(r8)::logstep                   !value of dlnd
  real(r8)::wnd_frc_thr(max_sizes)    !friction wind thr
  real(r8)::wnd_frc(max_ustr)         !friction wind
  real(r8)::wnd_frc_rat               !ratio of wind/thr wind
  real(r8)::ds(max_sizes)             !value of ds for a given d
  real(r8)::size_segr_hrz_flx(max_ustr,max_sizes)                 !size segregated horizontal flux
  real(r8)::fraction(max_ustr,max_sizes)                          !fraction of horizontal flux
  real(r8)::total                     !checking variable
  integer::wnd_frc_idx,sz_idx       !counting variables
  ! variables used by Zender dust module
  integer,parameter::dst_src_nbr=3                !number of source modes
  integer::src_idx                      !counting var. for source modes
  real(r8)::dmt_vma_src(dst_src_nbr)     ![m] mmd of source modes
  real(r8)::sgm_src(dst_src_nbr)         ![-] standard deviation of source modes
  real(r8)::dmt_vmm_src(dst_src_nbr)     ![m] mass mean source modes
  real(r8)::egy_bnd_src(dst_src_nbr)     ![j] binding energy for source modes
  real(r8)::egy_frc_src(dst_src_nbr)     ![-] energy fraction of source modes
  real(r8)::egy_frc_src_tmp(dst_src_nbr) ![-] energy fraction of source modes
  real(r8)::mss_frc_src(dst_src_nbr)     ![-] mass fraction of source modes
  real(r8)::mss_frc_src_tmp(dst_src_nbr) ![-] mass fraction of source modes for iteration
  real(r8)::nbr_frc_src(dst_src_nbr)     ![-] number fraction of source modes
  real(r8)::nbr_frc_src_tmp(dst_src_nbr) ![-] number fraction of source modes for one iteration
  real(r8)::dst_slt_flx_rat_ttl          ![m-1] alpha (vert/hor flux)
  real(r8)::dst_slt_flx_rat_ttl_tmp      ![m-1] alpha for one iteration
  real(r8)::flx_mss_vrt_dst_ttl          ![kg m-2 s-1] vertical mass flux 
  real(r8)::flx_mss_vrt_dst_ttl_tmp      ![kg m-2 s-1] vertical mass flux for one iteration
  real(r8)::flx_nbr_vrt_dst_ttl          ![# m-2 s-1] vertical number flux
  real(r8)::flx_nbr_vrt_dst_ttl_tmp      ![# m-2 s-1] vertical number flux for one iteration
  real(r8)::flx_mss_hrz_slt_ttl          ![kg m-1 s-1] horizontal flux
  real(r8)::flx_mss_hrz_slt_ttl_tmp      ![kg m-1 s-1] horizontal flux for one iteration
  real(r8),parameter::dns_soil=2650.0_r8             !density of soil      
  real(r8)::dmt_slt_bin                  ![m] diameter of saltating bin 
  logical::flg_mbl             !mobilization flag
  integer::soil_type_idx       !Counter for soil type
  integer,parameter::fileno=21           !file number of outfile
  ! output to file
  real(r8)::alphavalue(max_ustr,max_soil_type)                     ![m-1] value of alpha for a given u* and soiltype
  real(r8)::alphanumbervalue(max_ustr,max_soil_type)               ![#/kg/m] number based alpha for a given u* and soiltype
  real(r8)::mss_frc_src_value(max_ustr,max_soil_type,dst_src_nbr)  ![-] value of mss_frc_src for given u* and soiltype
  real(r8)::horizontal_saltation_flux(max_ustr,max_soil_type)      ![kg m-1 s-1] horizontal soil flux for a given soil and u*
  character(len=50)::soilname(max_soil_type)   ![-] name of soil types
  character(len=50)::infilename                ![-] name of file to read in
  data dmt_vma_src /1.5d-6, 6.7d-6, 14.2d-6 /
  data sgm_src     /1.7d0, 1.6d0, 1.5d0 /
  data egy_bnd_src /3.61d-7, 3.52d-7, 3.46d-7 /
  data mmd_soil    /125.d-6, 210d-6, 520d-6, 690d-6/
! data sigma       /1.8, 1.6, 1.6, 1.6 / ! MBA97 values have a typo in Table fxm
  data sigma       /1.6, 1.8, 1.6, 1.6 / ! Corrected values follow Chatenet96 which also agree with Alf01
  data soilname    /"ass ", "fs", "ss" , "cs" /
  flg_mbl=.true.
  
  ! calculate mass mean diameters of source mode
  do src_idx=1,dst_src_nbr
     dmt_vmm_src(src_idx)=dmt_vma_src(src_idx) &
          *exp(log(sgm_src(src_idx))*log(sgm_src(src_idx))*(-1.5d0))
  enddo
  write(6,*) "median",dmt_vma_src
  write(6,*) "mean ",dmt_vmm_src
  
  !  stop
  ! read in values from Zender's mie program
  if (max_sizes == 100) then
     infilename="input_100.txt"
  elseif (max_sizes == 1000) then
     infilename="input_1000.txt"
  elseif (max_sizes == 10000) then
     infilename="input_10000.txt"
  else
     write(6,*) "check if you have input for"
     write(6,*)max_sizes ,"soil sizes"
     stop
  endif
  open(20,file=infilename)
  do n=1,max_sizes
     read(20,*)trash1,sz(n),trash2,wnd_frc_thr(n),trash3
     !     write(6,*)n,trash1,sz(n),trash2,wnd_frc_thr(n),trash3
  enddo
  
  ! conversion to from radius to diameter
  sz(:)=sz(:)*2.0_r8      
  
  ! setting array of wind frictions
  ! to values from 0-1 ms-1
  do wnd_frc_idx=1,max_ustr
     wnd_frc(wnd_frc_idx)=float(wnd_frc_idx)*1.0_r8/100.0_r8
  enddo
  
  ! calculating dlnd for use in size distribution function
  logstep=(log(sz(max_sizes))-log(sz(1)))/float(max_sizes)
  
  ! starting loop for a soil type
  do soil_type_idx=1,max_soil_type
     ! first part of program is to generate the array called 
     ! "fraction" which is the fraction of the horizontal (saltating)
     ! flux with a certain size
     ! using formulas7.50 / 7.52 in seinfeld
     ! number median diameter of soil :
     nmd_soil=exp(log(mmd_soil(soil_type_idx))-3.0_r8*log(sigma(soil_type_idx))*log(sigma(soil_type_idx)))
     ! surface median diameter of soil :
     smd_soil=exp(log(nmd_soil)+2.0_r8*log(sigma(soil_type_idx))*log(sigma(soil_type_idx)))
     ! calculating value for ds/dlnd for the classes we have
     ! calculating formula 7.33 in seinfeld for surface
     ! in the and, multiplying by logstep to get ds for a given size
     do n=1,max_sizes
        ds(n)=total_surface &
             /(sqrt(2.0_r8*pi)*log(sigma(soil_type_idx))) &
             *exp(-1.0_r8*((log(sz(n))-log(smd_soil))**2 &
             /(2.0_r8*log(sigma(soil_type_idx)) &
             *log(sigma(soil_type_idx))))) &
             *logstep
     enddo
     ! copy from aerosol model (zender dust module) program
     ! fraction(wnd_frc_idx,sz_idx) is an array containing how much of the 
     ! horizontal (saltating) flux has that specific size
     total=0.0_r8
     do wnd_frc_idx=1,max_ustr
        ! initializing "total" for all wind friction speeds
        total=0.0_r8
        do sz_idx=1,max_sizes
           wnd_frc_rat=wnd_frc_thr(sz_idx)/wnd_frc(wnd_frc_idx) ! [frc]
           !   horizontal flux for each size:      
           size_segr_hrz_flx(wnd_frc_idx,sz_idx)= & ! [kg m-1 s-1] 
                max(0.0_r8, &
                cst_slt*dns_mdp*(wnd_frc(wnd_frc_idx)**3.0)* &
                (1.0-wnd_frc_rat)*(1.0+wnd_frc_rat)* &
                (1.0+wnd_frc_rat) &
                /grv_sfc  &
                *ds(sz_idx))
           !write(6,*)'wnd_frc',wnd_frc_thr(sz_idx),wnd_frc(wnd_frc_idx),ds(sz_idx),size_segr_hrz_flx(wnd_frc_idx,sz_idx)
           ! summing up total flux for that size:
           total=total+size_segr_hrz_flx(wnd_frc_idx,sz_idx)
        enddo               !loop on size
        if (total > 0.0_r8) then
           fraction(wnd_frc_idx,:)=size_segr_hrz_flx(wnd_frc_idx,:)/total
        else
           fraction(wnd_frc_idx,:)=0.0_r8
        endif
        !Archive the horizonal saltation flux (kg m-1 s-1)
        horizontal_saltation_flux(wnd_frc_idx,soil_type_idx) = total
        write(6,*)'TOTAL ',trim(soilname(soil_type_idx)),wnd_frc_idx,total
     enddo                  !loop on wind friction
     !stop
     ! check if we have a total of 1 for all winds 
     do wnd_frc_idx=1,max_ustr
        total=0.0_r8
        do sz_idx=1,max_sizes
           total=total+fraction(wnd_frc_idx,sz_idx)
        enddo
     enddo
     ! second part of program is to use the array fraction to calculate 
     ! the mass fraction of each source mode and the value of "alpha" 
     ! and then print this to a look up table to use in the dust module
     flx_mss_hrz_slt_ttl=1.0_r8
     do wnd_frc_idx=1,max_ustr
        flx_mss_vrt_dst_ttl= 0.0_r8
        flx_nbr_vrt_dst_ttl=0.0_r8
        nbr_frc_src(:)=0.0_r8
        mss_frc_src(:)=0.0_r8
        write(6,*) "calculation for ustr",wnd_frc(wnd_frc_idx)
        do sz_idx=1,max_sizes
           ! diameter of the saltating size
           dmt_slt_bin=sz(sz_idx)
           ! set total horizontal flux to 1.
           ! then horizontal flux in loop is equal to "fraction"
           flx_mss_hrz_slt_ttl_tmp=flx_mss_hrz_slt_ttl*fraction(wnd_frc_idx,sz_idx)
           ! calculate the mass fraction of each source mode
           ! using the alg01 formulation
           call mss_frc_src_alg01_get( &
                wnd_frc(wnd_frc_idx), & ! i [m/s] wind friction velocity
                egy_bnd_src,          & ! i [j] binding energy of source aerosols
                mss_frc_src_tmp,      & ! o [frc] mass fraction of each source mode for specific diameter
                dmt_vmm_src,          & ! i  [m] mass mean diameter of source distributions  
                !               dmt_vma_src,          & ! i [m] mass median diamter of sourc distributions
                egy_frc_src_tmp,      & ! o [-] energy fraction used on each mode (p(i) in alg01) for specific diameter
                dmt_slt_bin,          & ! i [m] diameter of saltating bin
                dns_soil,             & ! [kg m-3] density of soil
                dst_src_nbr           & ! i [-] number of source modes
                )
           ! calculate alpha and total vertical flux
           ! done for a all longitudes
           call flx_mss_vrt_dst_ttl_alg01_get( &
                dst_slt_flx_rat_ttl_tmp,  & ! o [frc] ratio of vertical to horizontal dust flux
                flg_mbl,                  & ! i [flg] mobilization flag
                mss_frc_src_tmp,          & ! i [frc] mass fraction source
                flx_mss_hrz_slt_ttl_tmp,  & ! i [kgms-1] total horizontal flux of soil
                flx_mss_vrt_dst_ttl_tmp,  & ! o [kg m-2 s-1] total vertical dust flux
                dmt_vmm_src,              & ! i [m] mass mean diameter of source
                !              dmt_vma_src,              & ! i [m] mass median diamter of source
                egy_bnd_src,              & ! i [j] bounding energy of source dust
                egy_frc_src_tmp,          & ! i [-] energy fraction use in each mode (p(i) in alg01)
                dns_soil,                 & ! i [kg m-3] density of soil
                dst_src_nbr,              & ! i [-] number of source modes
                flx_nbr_vrt_dst_ttl_tmp   & ! o [# m-2 s-1] number flux of dust
                )                          
           ! summing up the total dust flux
           flx_mss_vrt_dst_ttl=flx_mss_vrt_dst_ttl+flx_mss_vrt_dst_ttl_tmp     
           ! summing up total number flux
           flx_nbr_vrt_dst_ttl=flx_nbr_vrt_dst_ttl+flx_nbr_vrt_dst_ttl_tmp
           ! summing up mass fraction of source weighted by vertical dust flux fsd
           mss_frc_src(:)=mss_frc_src_tmp(:)*flx_mss_vrt_dst_ttl_tmp+mss_frc_src(:)
        enddo               !sz_idx
        !Make sure we don't divide by zero
        mss_frc_src(:)=mss_frc_src(:)/max(1.e-40,flx_mss_vrt_dst_ttl)
        if (flx_mss_vrt_dst_ttl > 0.0_r8) then
           alphavalue(wnd_frc_idx,soil_type_idx)=flx_mss_vrt_dst_ttl/flx_mss_hrz_slt_ttl
           alphanumbervalue(wnd_frc_idx,soil_type_idx)=flx_nbr_vrt_dst_ttl/flx_mss_hrz_slt_ttl
           mss_frc_src_value(wnd_frc_idx,soil_type_idx,:)=mss_frc_src(:)         
        else
           alphavalue(wnd_frc_idx,soil_type_idx)=0.0_r8
           alphanumbervalue(wnd_frc_idx,soil_type_idx)=0.0_r8
           mss_frc_src_value(wnd_frc_idx,soil_type_idx,:)= 0.0_r8
        endif
        write(6,*) "total mass flux ",flx_mss_vrt_dst_ttl
        write(6,*) "total nmbr flux ",flx_nbr_vrt_dst_ttl
     enddo                  !ustar
  enddo                     !soiltype
  ! part 3 of program is to write the results to a nice file
  ! which can be used by the zender dust module
  open(fileno,file="AlG01data.dat",form="formatted")
  do soil_type_idx=1,max_soil_type
     write(fileno,1001)soilname(soil_type_idx)
     write(fileno,1002) "ustr","alpha","numberalpha", &
          "massfrac(1)","massfrac(2)", &
          "massfrac(3)"
     do wnd_frc_idx=1,max_ustr
        write(fileno,1003)wnd_frc(wnd_frc_idx), &
             alphavalue(wnd_frc_idx,soil_type_idx), &
             alphanumbervalue(wnd_frc_idx,soil_type_idx), &
             (mss_frc_src_value(wnd_frc_idx,soil_type_idx,src_idx),src_idx=1,dst_src_nbr)
     enddo  !loop on wnd_frc
  enddo     !loop on soil_type
  close(fileno)
1001 format(a50)
1002 format(6a15)
1003 format(6e15.6)
  ! output to matlab
  open(fileno, file="formatlab.m", form="formatted")
  write(fileno,1010) "ustar",1,"=[", &
       (wnd_frc(wnd_frc_idx),wnd_frc_idx=1,max_ustr),"]"
  do soil_type_idx=2,2 !max_soil_type
     write(fileno,1010) "alpha",soil_type_idx," = [", &
          (alphavalue(wnd_frc_idx,soil_type_idx) &
          ,wnd_frc_idx=1,max_ustr),"]"
  enddo
  !  do soil_type_idx=2,2 !max_soil_type
  !     write(fileno,1011) "numberalpha",soil_type_idx," = [", &
  !        (alphanumbervalue(wnd_frc_idx,soil_type_idx) &
  !        ,wnd_frc_idx=1,max_ustr),"]"
  !  enddo
1010 format(a5,i1,a4,100e15.6,a1)
  ! 1011 format(a11,i1,a4,100e15.6,a1)
  open(73,file="dstsltsbl.FF",form="formatted")
  do soil_type_idx=1,max_soil_type
     write(73,1022) "    ! Soil type ",soil_type_idx, &
          ", MMD = ",mmd_soil(soil_type_idx)*1e6,"microns", &
          ", Geometric standard deviation (sigma) = ",sigma(soil_type_idx)
     write(73,"(a6,a50)") " ","! First index is wind friction speed in cm s-1 [1..100]"
     write(73,"(a6,a50)") " ","! Second index is soil type index [1..4]"
     write(73,"(a6,a50)") " ","! Third index is number of source mode [1..3]"
     do wnd_frc_idx=1,max_ustr
        write(73,1021) &
             "      ", &
             "data ", &
             "(mss_frc_src_lktx(",wnd_frc_idx,",",soil_type_idx, &
             ",src_idx),src_idx=1,dst_src_nbr) /" &
             ,mss_frc_src_value(wnd_frc_idx,soil_type_idx,1),"," &
             ,mss_frc_src_value(wnd_frc_idx,soil_type_idx,2),"," &
             ,mss_frc_src_value(wnd_frc_idx,soil_type_idx,3), &
             "/"
     enddo
  enddo
1021 format(a6,a5,a18,i3,a1,i1,a35,3(e15.6,a1))
1022 format(a20,i4,a6,f8.1,a4,a8,f8.1)
  write(73,"(a6,a50)") " ","!alpha values, vert/hor flux [m-1]"
  do soil_type_idx=1,max_soil_type
     write(73,"(a6,a50)") " ","! first index is u* in cm/s"
     write(73,"(a6,a50)") " ","! second index soil type (see above)"
     do wnd_frc_idx=1,max_ustr
        write(73,1023) " ","data ","dst_slt_flx_rat_ttl_lktx (" &
             ,wnd_frc_idx,",",soil_type_idx, &
             ") /" &
             ,alphavalue(wnd_frc_idx,soil_type_idx), &
             "/"
     enddo
  enddo
1023 format(a6,a5,a26,i3,a1,i1,a3,e15.6,a1)

  ! Write out total horizontal fluxes too. We need them to be consistent
  write(73,*)'       '
  write(73,*)"! Horizontal flux for different soils"
  do soil_type_idx=1,max_soil_type
     write(73,*)'   '
     do wnd_frc_idx=1,max_ustr
        write(73,1025) &
             "      ", &
             "data ", &
             "flx_mss_hrz_slt_ttl_lktx(",wnd_frc_idx,",",soil_type_idx,") /" &
             ,horizontal_saltation_flux(wnd_frc_idx,soil_type_idx),"/"

1025    format(a6,a5,a28,i3,a1,i1,a3,e15.6,a3)
     enddo
  enddo
end program sltsbl

subroutine flx_mss_vrt_dst_ttl_alg01_get( &
     dst_slt_flx_rat_ttl,    & ! o [-] ratio of vert/hor. flux
     flg_mbl,                & ! i [flg] mobilization candidate flag
     mss_frc_src,            & ! i [-] mass fraction in source
     flx_mss_hrz_slt_ttl,    & ! i [kg m-1 s-1] vertically integrated streamwize flux
     flx_mss_vrt_dst_ttl,    & ! o [kg m-2 s-1] total vertical dust flux
     dmt_vma_src,            & ! i [m] median diameters of source 
     egy_bnd_src,            & ! i [j] bounding energy of source dust
     egy_frc_src,            & ! i [-] energy fraction used in each mode (p(i) in alg01)           
     dns_soil,               & ! i [kg m-3] density of soil
     dst_src_nbr,            & ! i [-] number of source modes
     flx_nbr_vrt_dst_ttl     & ! o [# m-2 s-1] number flux of dust
     )
  ! Purpose:
  ! calculate alpha (vertflux/horflux) and total vertical flux of dust 
  ! for a given friction velocity,
  ! density and size distribution of soil and source aerosol.
  ! formulation is given by alfaro and gomes (2001) (eq.6).
  ! will be called from dstmbl. will replace old routine called
  ! flx_mss_vrt_dst_ttl_mab95_get where alpha is calculated as a function
  ! of clay content.
  ! coded by Alf Grini, spring 2002, UCI/UIO
  ! ---------------------------------------------------------------
  use precision ! [mdl] precision r8, i8, ...
  implicit none
  ! ----------------------------------------------------
  ! output
  real(r8),intent(out)::dst_slt_flx_rat_ttl      ![m-1] vert / hor. flux
  real(r8),intent(out)::flx_mss_vrt_dst_ttl      ![kg m-2 s-1] total vertical mass flux
  real(r8),intent(out)::flx_nbr_vrt_dst_ttl      ![# m-2 s-1] total vertical number flux
  ! input
  integer,intent(in)::dst_src_nbr                      ! [-] number of source modes
  logical,intent(in)::flg_mbl                          ! [flg] mobilization candidate (?)
  real(r8),intent(in)::mss_frc_src(dst_src_nbr) ! [-] mass fraction source (p(i) in a/g))
  real(r8),intent(in)::flx_mss_hrz_slt_ttl      ! [kg m-1 s-1] vertically integrated streamwize flux
  real(r8),intent(in)::dmt_vma_src(dst_src_nbr) ! [m] mass median diameter of source
  real(r8),intent(in)::egy_bnd_src(dst_src_nbr) ! [j] binding energy of source modes
  real(r8),intent(in)::egy_frc_src(dst_src_nbr) ! [-] energy fraction used in each source mode (p(i) in alg01)
  real(r8),intent(in)::dns_soil                 ! [kg m-3] density of soil 
  ! local in this routine
  integer::lon_idx, src_idx                ! counting variables
  real(r8),parameter::pi=3.141592654d0      ! pi
  real(r8),parameter::pisixths=pi/6.0_r8      ! pi/6
  real(r8),parameter::beta=16300.d-2        ! ms-2 (given alg01-article)
  ! begin code
  ! initializing alpha
  dst_slt_flx_rat_ttl=0.0_r8
  ! initializing number flux
  flx_nbr_vrt_dst_ttl=0.0_r8
  if (flg_mbl) then
     do src_idx=1, dst_src_nbr !loop on source distributions
        ! solving eq. 6 in alfaro/gomes.
        ! finding alpha (= vertical flux/horizontal flux)
        ! in the limit of very high kinetic energy of soil aggregates, 
        ! the mss_frc_src will be close to 1 for the smallest mode,
        ! and 0 for the two largest ones 
        ! then alpha will reach a limit value  
        ! alpha is sum of contribution from 3 source distributions:
        ! solving eq 6 in alg01 :
        dst_slt_flx_rat_ttl= &
             dst_slt_flx_rat_ttl+ &
             pisixths*dns_soil*beta  &
             *egy_frc_src(src_idx) &
             *dmt_vma_src(src_idx)**3 &
             /egy_bnd_src(src_idx)
     enddo                  !loop on source distributions
     ! getting number flux of each mode in #/m2/s
     ! eqn 3 in alg01
     do src_idx=1,dst_src_nbr
        flx_nbr_vrt_dst_ttl= &
             beta*flx_mss_hrz_slt_ttl*egy_frc_src(src_idx)/egy_bnd_src(src_idx) &
             +flx_nbr_vrt_dst_ttl
     enddo
     !   finding the total vertical dust flux from alpha and horizontal
     !   for each longitude:
     flx_mss_vrt_dst_ttl=dst_slt_flx_rat_ttl*flx_mss_hrz_slt_ttl
  endif                     !check on mobilization flag
  return
end subroutine flx_mss_vrt_dst_ttl_alg01_get

subroutine mss_frc_src_alg01_get( &
     wnd_frc,             & ! i [m/s] wind friction velocity
     egy_bnd_src,         & ! i [kgm2s-2] binding energy of sources
     mss_frc_src,         & ! o [-] mass fraction of each of the three src dst.)
     dmt_vma_src,         & ! i [m] mass median diamter of source functions
     egy_frc_src,         & ! o [-] fraction of energy used in each source mode (p(i) in alg01)
     dmt_slt_bin,         & ! i [m] diameter of the saltating soil
     dns_soil,            & ! i [kg m-3] density of soil
     dst_src_nbr          & ! i [-] number of source modes
     )
  ! Purpose:
  ! determine based on friction wind speed how much energy is available
  ! for saltation. based on this energy, we can determine the mass fraction
  ! of each of the three distributions proposed by alfaro and gomes.
  ! called from dstmbl
  ! code by Alf Grini, spring 2002, UCI/UIO
  ! ------------------------------------------------------------
  use precision ! [mdl] precision r8, i8, ...
  implicit none
  ! -------------------------------------------------------------
  ! input
  integer,intent(in)::dst_src_nbr              ![-] number of source modes
  real(r8),intent(in)::wnd_frc                  ![m/s] wind friction vel 
  real(r8),intent(in)::egy_bnd_src(dst_src_nbr) ![j] binding energy of soil
  real(r8),intent(in)::dmt_vma_src(dst_src_nbr) ![m] mass median diameter of source functions
  real(r8),intent(in)::dmt_slt_bin              ![m] diameter of saltating bin
  real(r8),intent(in)::dns_soil                 ! [kg m-3] density of soil
  ! output
  real(r8),intent(out)::mss_frc_src(dst_src_nbr)  ![-] fraction of mass in each or three source distr.
  real(r8),intent(out)::egy_frc_src(dst_src_nbr)  ![-] fraction of energy used in three source distr.
  ! local in routine
  real(r8)::egy_kin_agr             ! [kgm2s-2] kinetic energy of soil aggregates
  real(r8),parameter::pi=3.141592654d0      ! [-] the number pi
  real(r8),parameter::onetwelwth=1.0_r8/12.0_r8 ! [-] one divided by twelve
  real(r8)::egy_mss_src_ttl         ! [m3j-1] total energy weighted by volume of modes
  integer::lon_idx                 ! [-] counting variable for longitude
  integer::src_idx                 ! [-] counting variable for source distributions
  real(r8)::nbr_frc_src(dst_src_nbr)! [-] number fraction of source
  ! the new thing here (compared to she84) , is that mss_frc_src 
  ! is calculated in each time step. 
  ! initializing
  mss_frc_src(:)=0.0_r8     !no mass if egy_kin_agr is too small
  egy_frc_src(:)=0.0_r8
  nbr_frc_src(:)=0.0_r8
  !  write(6,*) "starting loop on longitude"
  !   kinetic energy of aggregates : (eq. 1 in alg01)
  egy_kin_agr= &
       dns_soil*pi*onetwelwth &
       *dmt_slt_bin*dmt_slt_bin*dmt_slt_bin &
       *(20.0_r8*wnd_frc) &
       *(20.0_r8*wnd_frc) !units [kgm2s-2 = j]
  !     write(6,*) "dmt_slt_bin", dmt_slt_bin*1.d6,"um"
  !     write(6,*) "egy_kin_agr  ",egy_kin_agr,"j"
  !     write(6,*) "wnd_frc",wnd_frc,"m/s"
  !     write(6,*) "egy_bnd_src",egy_bnd_src
  !     write(6,*) "dns_soil ",dns_soil,"kg/m3"
  !     write(6,*) "onetwelwth*12",onetwelwth*12.0_r8
  !     write(6,*) "pi ",pi
  ! this if-test is basically the sum of table 2. in alg01 :
  if (egy_kin_agr > egy_bnd_src(3).and.egy_kin_agr < egy_bnd_src(2)) then
     ! only largest mode available
     egy_frc_src(3)=1.0_r8
     !        write(6,*) "only largest mode"
  elseif (egy_kin_agr > egy_bnd_src(2).and.egy_kin_agr < egy_bnd_src(1)) then
     ! only two largest modes available
     !        write(6,*) "only largest two modes"
     egy_frc_src(2)=(egy_kin_agr-egy_bnd_src(2))/(egy_kin_agr-egy_bnd_src(3))
     egy_frc_src(3)=1.0_r8-egy_frc_src(2)
  elseif (egy_kin_agr > egy_bnd_src(1)) then
     ! all three modes available
     !        write(6,*) "all three modes"
     egy_frc_src(1)=(egy_kin_agr-egy_bnd_src(1))/(egy_kin_agr-egy_bnd_src(3))
     egy_frc_src(2)=(1.0_r8-egy_frc_src(1))*(egy_kin_agr-egy_bnd_src(2))/(egy_kin_agr-egy_bnd_src(3))
     egy_frc_src(3)=1.0_r8-egy_frc_src(2)-egy_frc_src(1)
  endif
  ! initializing
  egy_mss_src_ttl=0.0_r8          
  ! summing total kinetic energy used for saltiation weighted by volume
  ! this is the denominator in eq: 9 alg01
  do src_idx=1, dst_src_nbr
     egy_mss_src_ttl=egy_mss_src_ttl &
          +egy_frc_src(src_idx)*dmt_vma_src(src_idx)**3/egy_bnd_src(src_idx)
  enddo
  ! finding mass fraction of each mode by using eq. (9) in alg01
  if (egy_mss_src_ttl > 0.0_r8) then
     do src_idx=1,dst_src_nbr          
        mss_frc_src(src_idx)=egy_frc_src(src_idx)*dmt_vma_src(src_idx)**3 &
             /egy_bnd_src(src_idx) &
             /egy_mss_src_ttl
     enddo               !loop over source
     do src_idx=1,dst_src_nbr
        nbr_frc_src(src_idx)=(egy_frc_src(src_idx)/egy_bnd_src(src_idx))/ &
             ((egy_frc_src(1)/egy_bnd_src(1))+(egy_frc_src(2)/egy_bnd_src(2))+(egy_frc_src(3)/egy_bnd_src(3)))*100.0_r8 
     enddo
  else
     mss_frc_src(:)=0.0_r8 !no saltation
  endif
  
  ! if ekin lower than ekin_src(3) the if-loop will not be done
  ! and we will end here. mass fractions will be zero. 
  ! note!!!! :it has to be like this, because with 
  ! with alg01 formulation we can have saltation but no sandblasting.
  ! ( if u*>u*,tr, but ec<e(3)) 
  return
end subroutine mss_frc_src_alg01_get
