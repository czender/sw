! $Id$ -*-f90-*-

! Purpose: Library of geomorphological/hyrdological routines

! Copyright (C) 2002--2005 Charlie Zender

! Given a digital elevation model "hgt_sfc", and surface runoff "sfc_flw", 
! dst_src_flw computes area accumulation factor "sfc_acm"
! and flow accumulation factor "flw_acm".
! Factors are normalized to [0.0,1.0]

! Flow direction and flow accumulation algorithm is based on
! Jenson S. K., and Domingue J. O.,
! "Extracting Topographic Structure from Digital Elevation Data
! for Geographic Information System Analysis"
! Photogrammetric Engineering and Remote Sensing, Vol 54, No 11, Nov 1988 pp 1593-1600

! Author: David Newman
! newman@uci.edu
! Feb 2002


! Modified by Chao Luo, Dec. 2005
!
! Notes:
! 1. dst_src_flw does NOT first preprocess hgt_sfc to create a
!    depressionless elevation model (i.e. one where all land points
!    drain to the ocean).  Therefore internal basins may be produced.

module geo_mdl ! [mdl] Geomorphology and flow accumulation

  implicit none
  integer::geo_dbg_lvl=9
  public ::dst_src_flw ! [sbr] Compute area and flow accumulation factors
  private::initialize                ! set initilaization for some parameters
  private::flw_drc             ! flow direction calculation subroutine
  private::flw_acm_get              ! upstream area accumulation  subroutine
  private::accumulate_one              ! upstream area accumulation calculation
  private::get_neighbor              ! get longitude and latitude index of neighbor cell
  private::get_cells                 ! get neighbor cell's hight
  private::is_close                  ! check neighbour elevations 
  private::sort_elevation            ! sort elevation 
  private::myqsort                   ! set elevation in order 
  private::write_matlab_dir          ! write out the flow direction
  private::chk_flw                   ! check flow direction
  private::bsn_sz_get                ! get Basin size
contains

!-------------------------------------------------------------------
! subroutine dst_src_flw 
!-------------------------------------------------------------------
   subroutine dst_src_flw ( & ! [sbr] Compute area and flow accumulation factors
              lat_nbr, &   ! I  [nbr] number of longitude
              lon_nbr, &   ! I  [nbr] number of latitude
              hgt_sfc, &   ! I  [m] surface elevation
              oro_real, &  ! I  orography
              sfc_area, &  ! I  [m2] surface area
              sfc_flw,  &  ! I  surface flow
              bsn_enm,  &  ! O  drainage basin
              flw_acm,  &  ! O  flow acumlate
              sfc_acm,  &  ! O  [m2] surface acumulate
              bsn_sz)    ! O  [m2] Basin size
          
! Purpose: Compute area accumulation factor and flow accumulation factor
! dst_src_flw() called by bds_prc() in bds_drv.F90

  use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
  use bds_ctl,only:bkf_itr_nbr ! [mdl] Control variables drc_in,drc_out,hgt_dlt_msl
  implicit none

  integer,intent(in)::lat_nbr    ![nbr] number of latitudes
  integer,intent(in)::lon_nbr    ![nbr] number of longitudes
  real,intent(inout)::hgt_sfc(lon_nbr,lat_nbr)  ! [m] surface elevation
  real,intent(in)::oro_real(lon_nbr,lat_nbr)  ! Orography
  real,intent(in)::sfc_area(lon_nbr,lat_nbr)    ! [m2] surface area
  real,intent(in)::sfc_flw(lon_nbr,lat_nbr)   ! surface flow
  real,intent(out)::sfc_acm(lon_nbr,lat_nbr)   ! [m2] surface accumulate
  real,intent(out)::flw_acm(lon_nbr,lat_nbr)  ! flow accumulate
  integer,intent(out)::bsn_enm(lon_nbr,lat_nbr)  ! [nbr] drainage basin
  integer::hgt_order(lon_nbr*lat_nbr)     ! [m] surface elvation order
  integer::hgt_map(lon_nbr*lat_nbr,2)     ! [m] 1-D surface elevation array
  logical::oro(lon_nbr,lat_nbr)  !Orography flag
  integer::flw_dir(lon_nbr,lat_nbr)     !flow direction
  integer::is_flwd(lon_nbr,lat_nbr)     !flow direction flag
  logical::is_acm(lon_nbr,lat_nbr)      ! flag
  
  real::sfc_acm_orig(lon_nbr,lat_nbr)   ! [m2] surface accumulate
  real::dlt_hgt(lon_nbr,lat_nbr)        ! [m] surface elevation different
  integer::pointsTo(lon_nbr,lat_nbr)
  integer::sink(lon_nbr,lat_nbr)
  integer::internal(lon_nbr,lat_nbr)
  integer::flw(lon_nbr,lat_nbr)
  integer::numBasins           ! [nbr] number of basin
  integer::bck_nbr_idx         !index of Number of backfill iterations
  real bsn_sz(lon_nbr*lat_nbr)
  
  ! Main code
!  if (dbg_lvl >= dbg_sbr) write(6,*) "Entering dst_src_flw()"

! set initialization
  call initialize(lon_nbr, &    !I  [nbr] number of longitude
                  lat_nbr, &    !I  [nbr] number of latitude
                  hgt_sfc, &    !I  [m] surface elevation
                  oro_real,&    !I  [nbr] orography
                  oro,     &    !O  [nbr] orography
                  flw_dir, &    !O  flow direction
                  is_flwd, &    !O  flow direction flag
                  is_acm)       !O  flag for accumulation 
! set surface high initial
  call preprocess(lon_nbr, &    !I [nbr] number of longitude
                  lat_nbr, &    !I [nbr] number of latitude
                  oro,     &    !I [nbr] orography
                  hgt_sfc)      !O [m] surface elevation

  do bck_nbr_idx=1,bkf_itr_nbr
! set initial
     call init(lon_nbr, &    !I [nbr] number of longitude
               lat_nbr, &    !I [nbr] number of latitude
               flw,     &    !I/O flow direction
               pointsTo,&    !I/O
               sink,    &    !I/O drainage basin
               bsn_enm, &    !I/O drainage basin
               internal,&    !I/O internal basin
               dlt_hgt)      !O  [m] surface elevation difference
!set minimum surface high value
     call fillSingleCellDepression(lon_nbr, &    !I [nbr] number of longitude
                                   lat_nbr, &    !I [nbr] number of latitude
                                   oro,     &    !I [nbr] orography
                                   hgt_sfc)      !O [m] surface elevation
! calculate flow direction
     call flw_drc(lon_nbr, &    !I [nbr] number of longitude
                        lat_nbr, &    !I [nbr] number of latitude
                        oro,     &    !I [nbr] orography
                        hgt_sfc, &    !I [m] surface elevation
                        flw)          !O flow direction
! check that all flows are downhill
     call chk_flw_1(lon_nbr, &    !I [nbr] number of longitud
                    lat_nbr, &    !I [nbr] number of latitude
                    oro,     &    !I [nbr] orography
                    hgt_sfc, &    !I [m] surface elevation
                    flw)          !I/O flow direction
!check circular flow
     call chk_flw_2(lon_nbr, &    !I [nbr] number of longitud
                    lat_nbr, &    !I [nbr] number of latitude
                    oro,     &    !I [nbr] orography
                    hgt_sfc, &    !I [m] surface elevation
                    flw)          !I/O  flow direction
!finds all flat points and sees if there is a path to the surf
     call fixFlowDir(lon_nbr, &   !I [nbr] number of longitud
                     lat_nbr, &   !I [nbr] number of latitude
                     oro,     &   !I [nbr] orography
                     hgt_sfc, &   !I [m] surface elevation
                     flw)         !I/O flow direction
     !print *, "unresolved = ", unresolvedFlow(n, m, oro, flw)
     
! get points
     call getPointsTo(lon_nbr,  &  !I [nbr] number of longitud
                      lat_nbr,  &  !I [nbr] number of latitude
                      oro,      &  !I [nbr] orography
                      flw,      &  !I flow direction
                      pointsTo)    !O 
     numBasins = getSink (lon_nbr, &  !I [nbr] number of longitud
                          lat_nbr, &  !I [nbr] number of latitude
                          pointsTo,&  !I 
                          sink)       !O drainage basin
      print *, "number of sinks = ", numBasins
     
! Merge Sink
     call mergeSink(lon_nbr,  &   !I [nbr] number of longitud
                    lat_nbr,  &   !I [nbr] number of latitude
                    oro,      &   !I [nbr] orography
                    hgt_sfc,  &   !I [m] surface elevation
                    flw,      &   !I flow direction
                    sink)         !I/O drainage basin
! get Basin
     call getBasin(lon_nbr,   &   !I [nbr] number of longitud
                   lat_nbr,   &   !I [nbr] number of latitude
                   pointsTo,  &   !I 
                   sink,      &   !I drainage basin
                   bsn_enm)       !I/O drainage basin
! get internal
     call getInternal(lon_nbr, &  !I [nbr] number of longitud
                      lat_nbr, &  !I [nbr] number of latitude
                      oro,     &  !I [nbr] orography
                      sink,    &  !I drainage basin
                      internal)   !O basin
     !print *, "zzz num internal = ", countInternal(n, m, numBasins, internal)
     
     if (bck_nbr_idx < bkf_itr_nbr) then
        call fillBasin(lon_nbr,  &   !I  [nbr] number of longitud
                       lat_nbr,  &   !I  [nbr] number of latitude
                       numBasins,&   !I  [nbr] number of basins
                       oro,      &   !I  [nbr] orography
                       hgt_sfc,  &   !I  [m] surface elevation
                       bsn_enm,  &   !I  drainage basins
                       internal, &   !I  internal basin
                       dlt_hgt)      !O  [m] surface elevation difference
        !print *, "zzz dh           = ", sum(dh)
        call adjustHeight(lon_nbr, &  !I  [nbr] number of longitud
                          lat_nbr, &  !I  [nbr] number of latitude
                          hgt_sfc, &  !I  [m] surface elevation
                          dlt_hgt)    !O  [m] surface elevation difference
     end if

  end do ! end loop over bck_nbr_idx

!  calculate drainage area
   call bsn_sz_get (lon_nbr,  &   !I [nbr] number of Longitude
                      lat_nbr,  &   !I [nbr] number of latitude
                      oro,      &   !I  Orography flag
                      numBasins, & !I [nbr] number of Basins
                      sfc_area,  &  !I [m2] grid area
                      bsn_enm,  &   !I [nbr] Basin ID
                      bsn_sz)     !O flow direction


  flw_dir = flw

  flw_dir(:,1) = 0
  flw_dir(:,lat_nbr) = 0
  where (flw_dir < 0) flw_dir = 0
  
  call chk_flw_1 (lon_nbr, &    !I  [nbr] number of longitud
                  lat_nbr, &    !I  [nbr] number of latitude
                  oro,     &    !I  [nbr] orography
                  hgt_sfc, &    !I  [m] surface elevation
                  flw_dir)      !O  flow direction
  call chk_flw_2 (lon_nbr, &    !I  [nbr] number of longitud
                  lat_nbr, &    !I  [nbr] number of latitude
                  oro,     &    !I  [nbr] orography
                  hgt_sfc, &    !I  [m] surface elevation
                  flw_dir)      !O  flow direction
  call chk_flw   (lon_nbr, &    !I  [nbr] number of longitud
                  lat_nbr, &    !I  [nbr] number of latitude
                  oro,     &    !I  [nbr] orography
                  hgt_sfc, &    !I  [m] surface elevation
                  flw_dir)      !O  flow direction
!  Sort elevation
  call sort_elevation (lon_nbr, &  !I [nbr] number of longitud
                       lat_nbr, &  !I [nbr] number of latitude
                       oro,     &  !I [nbr] orography
                       hgt_sfc, &  !I [m] surface elevation
                       hgt_order, &  !O [m] 1-D surface elevation from low to high
                       hgt_map)      !O [m] 2-D surface elevation coordinator 
 
  !sfc_acm= 1.0      ! initialize sfc_acm with 1.0 for testing purposes
  sfc_acm = sfc_area ! initialize sfc_acm with sfc_area
! surface accumulate
  call flw_acm_get   (lon_nbr,   &    !I  [nbr] number of longitud
                       lat_nbr,   &    !I  [nbr] number of latitude
                       oro,       &    !I  [nbr] orography
                       hgt_sfc,   &    !I  [m] surface elevation
                       hgt_order, &    !I  [m] 1-D surface elevation from low to high
                       hgt_map,   &    !I  [m] 2-D surface elevation coordinator
                       flw_dir,   &    !I  flow direction
                       sfc_acm,   &    !O  [m2] surface accumulation
                       is_acm)         !O  flag for surface accumulation

  sfc_acm_orig = sfc_acm
  flw_acm  = sfc_flw  ! initialize flw_acm with sfc_flw
! flow accumulate
  call flw_acm_get   (lon_nbr,   &    !I  [nbr] number of longitud
                       lat_nbr,   &    !I  [nbr] number of latitude
                       oro,       &    !I  [nbr] orography
                       hgt_sfc,   &    !I  [m] surface elevation
                       hgt_order, &    !I  [m] 1-D surface elevation from low to high
                       hgt_map,   &    !I  [m] 2-D surface elevation coordinator
                       flw_dir,   &    !I  flow direction
                       flw_acm,   &    !O  surface flow accumulation
                       is_acm)         !O  flag for surface flow accumulation
  
  ! normalize factors to be in the range 0->1
  if (geo_dbg_lvl > 0) write(6,*) "INFO: max sfc_acm = ", maxval(sfc_acm)
  if (geo_dbg_lvl > 0) write(6,*) "INFO: max flw_acm = ", maxval(flw_acm)
  if (abs(maxval(sfc_acm)) > 1e-30) sfc_acm = sfc_acm/maxval(sfc_acm)
  if (abs(maxval(flw_acm))  > 1e-30) flw_acm  = flw_acm/maxval(flw_acm)

  ! use following matlab datafiles for debugging
   call write_matlab_dir (lon_nbr, lat_nbr, hgt_sfc, flw_dir, flw_acm)
  !call write_matlab     (lon_nbr, lat_nbr, hgt_sfc,     "data_hgt_sfc.m", "hgt_sfc=[")
  !call write_matlab     (lon_nbr, lat_nbr, oro_real, "data_oro.m",     "oro=[")
  !call write_matlab     (lon_nbr, lat_nbr, sfc_flw,  "data_sfc_flw.m", "sfc_flw=[")
  !call write_matlab     (lon_nbr, lat_nbr, flw_acm,  "data_flw_acm.m", "flw_acm=[")
  !call write_matlab     (lon_nbr, lat_nbr, sfc_acm, "data_sfc_acm.m","sfc_acm=[")
   call write_matlab     (lon_nbr, lat_nbr, sfc_acm_orig, "data_sfc_acm.m","sfc_acm=[")

  if (dbg_lvl >= dbg_sbr) write(6,*) "Exiting dst_src_flw()"
  
end subroutine dst_src_flw

!-------------------------------------------------------------------
! subroutine initialize 
!-------------------------------------------------------------------
subroutine initialize (lon_nbr, &    !I [nbr] Number of longitude
                       lat_nbr, &    !I [nbr] Number of latitude
                       hgt_sfc, &    !I [m] surface elevation
                       oro_real,&    !I orography
                       oro,     &    !O orography
                       flw_dir, &    !O flow direction
                       is_flwd, &    !O flow direction flag 
                       is_acm)       !O surface accumualtion flag
! Purpose: initilization of some parameters
! initialize() called by dst_src_flw in geo_mdl.F90 
  use dbg_mdl
  implicit none
  integer, intent(in)   ::lon_nbr     ! [nbr] Number of longitude
  integer, intent(in)   ::lat_nbr     ! [nbr] Number of latitude
  real,    intent(inout)::hgt_sfc (lon_nbr,lat_nbr) ![m] surface elevation
  real,    intent(in)   ::oro_real (lon_nbr,lat_nbr) !orography
  logical, intent(out)  ::oro      (lon_nbr,lat_nbr) !orography
  integer, intent(out)  ::flw_dir  (lon_nbr,lat_nbr) ! flow direction
  integer, intent(out)  ::is_flwd  (lon_nbr,lat_nbr) ! flow direction flag
  logical, intent(out)  ::is_acm   (lon_nbr,lat_nbr) ! surface accumualtion flag
  
  !if (geo_dbg_lvl > 0) write(6,*) "INFO: oro=1 and h<0 at", count (oro_real > 0.5 .and. h < 0.0), "cells"
  
  where (oro_real > 0.5)
     oro = .true.
  elsewhere
     oro = .false.
  end where
  
  where (hgt_sfc < 0.0)
     hgt_sfc = 0.0
  end where
  
  flw_dir = 0
  
  where (oro)
     is_flwd = -1
     is_acm  = .false.
  elsewhere
     is_flwd = 1
     is_acm  = .true.
  end where
  
end subroutine initialize


!-------------------------------------------------------------------
! subroutine get_neighbor 
!-------------------------------------------------------------------

subroutine get_neighbor (lon_nbr,  &   !I [nbr] number of Longitude
                         lat_nbr,  &   !I [nbr] number of latitude
                         lon_idx,  &   !I [nbr] longitude index
                         lat_idx,  &   !I [nbr] latitude index
                         dir,      &   !O flow direction
                         lon_neb_idx, &  !O [nbr] longitude index for neighbor cell
                         lat_neb_idx)    !O [nbr] latitude index for neighbor cell
! Purpose: get longitude and latitude index of neighbor cell 
! get_neighbor() called by get_cell in geo_mdl.F90
  implicit none
  integer, intent(in) ::lon_nbr
  integer, intent(in) ::lat_nbr
  integer, intent(in) ::lon_idx
  integer, intent(in) ::lat_idx
  integer, intent(in) :: dir
  integer, intent(out)::lon_neb_idx
  integer, intent(out) ::lat_neb_idx

  !------------
  ! 5 | 4 | 3 |
  !------------
  ! 6 | x | 2 |
  !------------
  ! 7 | 8 | 1 |
  !------------

! fxm: 20030720 change to select (case)  
  if      (dir==1) then
     lon_neb_idx = lon_idx+1
     lat_neb_idx = lat_idx-1
  else if (dir==2) then
     lon_neb_idx = lon_idx+1
     lat_neb_idx = lat_idx
  else if (dir==3) then
     lon_neb_idx = lon_idx+1
     lat_neb_idx = lat_idx+1
  else if (dir==4) then
     lon_neb_idx = lon_idx
     lat_neb_idx = lat_idx+1
  else if (dir==5) then
     lon_neb_idx = lon_idx-1
     lat_neb_idx = lat_idx+1
  else if (dir==6) then
     lon_neb_idx = lon_idx-1
     lat_neb_idx = lat_idx
  else if (dir==7) then
     lon_neb_idx = lon_idx-1
     lat_neb_idx = lat_idx-1
  else if (dir==8) then
     lon_neb_idx = lon_idx
     lat_neb_idx = lat_idx-1
  else
     stop "ERROR: [get_neighbor] bad dir"
  endif
  
  if (lon_neb_idx < 0 .or. lon_neb_idx > lon_nbr+1 .or. lat_neb_idx < 0 .or. lat_neb_idx > lat_nbr+1) stop "ERROR: [get_neighbor]"

  ! correct for spherical boundary conditions
  if (lon_neb_idx==0) then
     lon_neb_idx = lon_nbr
  endif
  if (lon_neb_idx==lon_nbr+1) then
     lon_neb_idx = 1
  endif
  if (lat_neb_idx==0) then
     lon_neb_idx = mod(lon_neb_idx-1+lon_nbr/2,lon_nbr) + 1
     lat_neb_idx = 1
  endif
  if (lat_neb_idx==lat_nbr+1) then
     lon_neb_idx = mod(lon_neb_idx-1+lon_nbr/2,lon_nbr) + 1
     lat_neb_idx = lat_nbr
  endif
  
  if (lon_neb_idx< 1 .or. lon_neb_idx> lon_nbr .or. lat_neb_idx < 1 .or. lat_neb_idx > lat_nbr) stop "ERROR: [get_neighbor]"

end subroutine get_neighbor


!-------------------------------------------------------------------
! subroutine get_cells 
!-------------------------------------------------------------------
subroutine get_cells (lon_nbr, &   !I [nbr] number of longitude
                      lat_nbr, &   !I [nbr] number of latitude
                      lon_idx ,&   !I [nbr] longitude index
                      lat_idx, &   !I [nbr] latitude index
                      hgt_sfc, &   !I [m] surface elevation
                      hgt_dir)     !O [nbr] surface elevation of 8 directions
! Purpose: get neighbor cell's hight
! get_cell() called by flow_direction in geo_mdl.F90
  implicit none
  integer, intent(in) ::lon_nbr
  integer, intent(in) ::lat_nbr
  integer, intent(in) ::lon_idx
  integer, intent(in) ::lat_idx
  real,    intent(in) ::hgt_sfc(lon_nbr,lat_nbr)
  real,    intent(out)::hgt_dir(8)
! local
  integer::dir
  integer::lon_neb_idx
  integer::lat_neb_idx
  do dir = 1, 8
     call get_neighbor (lon_nbr, &   ! I  [nbr] number of longitude
                        lat_nbr, &   ! I  [nbr] number of latitude
                        lon_idx, &   ! I  [nbr] longitude index
                        lat_idx, &   ! I  [nbr] latitude index
                        dir,     &   ! I  nbr] flow direction
                        lon_neb_idx,  & !O [nbr] number of longitude of neighbor cell
                        lat_neb_idx)    !O [nbr] number of latitude of neighbor cell
     hgt_dir(dir) = hgt_sfc(lon_neb_idx,lat_neb_idx)
  end do
end subroutine get_cells


!-------------------------------------------------------------------
! function is_close
!-------------------------------------------------------------------
function is_close(hgt1,hgt2)
!Purpose: check the neighbor elevation
!is_close() called by acumulate_one in geo_mdl.F90
    implicit none
    logical:: is_close
    real, intent(in):: hgt1, hgt2
    is_close = .false.
    if (abs(hgt1-hgt2) < 1.0e-6) is_close = .true.
end function is_close


!-------------------------------------------------------------------
! subroutine flw_drc 
!-------------------------------------------------------------------
subroutine flw_drc (lon_nbr, &    !I [nbr] Number of longitude
                          lat_nbr, &    !I [nbr] Number of latitude
                          oro,     &    !I Orography
                          hgt_sfc, &    !I [m] surface elevation
                          flw)          !O [nbr] flow direction
! Purpose: flow direction calculation
! flw_drc() called by dst_src_flw in geo_mdl.F90
  use dbg_mdl
  implicit none
  integer, intent(in)   ::lon_nbr  !I [nbr] Number of longitude
  integer, intent(in)   ::lat_nbr  !I [nbr] Number of latitude
  logical, intent(in)   ::oro(lon_nbr,lat_nbr) !Orography
  real,    intent(in)   ::hgt_sfc(lon_nbr,lat_nbr) ![m] surface elevation
  integer, intent(out)  ::flw(lon_nbr,lat_nbr) ![nbr] flow direction
! local
  integer::lon_idx        ! [nbr] Number of longitude index
  integer::lat_idx        ! [nbr] Number of latitude index
  integer::k_loc(1)       ! [nbr] index of location of minimum value of array
  integer::k_loc_min      ! [nbr] index of location of minimum value of array
  integer::dir            ! [nbr] flow direction
  integer::lon_neb_idx    ! [nbr] Number of longitude index of neighbor cell
  integer::lat_neb_idx    ! [nbr] Number of latitude index of neighbor cell
  integer::opp_dir        ! [nbr] opposite direction
  real ::hgt_dir(8)       ! [nbr] flow direction
  real::dlt_hgt(8)        ! [m]  surface elevation difference
  real::hgt_tmp           ! [m] surface elevation difference
  real::fac               ! [fct] 
  real::dlt_hgt_min       ! [m]  surface elevation difference
  logical::is_many_min    ! flag
  
  ! Main code
  !if (dbg_lvl >= dbg_sbr) write(6,*) "Entering flw_drc()"
  
  fac = 1.0/sqrt(2.0)
  
  do lat_idx = 1, lat_nbr
     do lon_idx = 1, lon_nbr
        if (.not. oro(lon_idx,lat_idx)) cycle
        call get_cells (lon_nbr,lat_nbr,lon_idx,lat_idx,hgt_sfc,hgt_dir)
        hgt_tmp     = hgt_sfc(lon_idx,lat_idx)
        dlt_hgt(1)  = fac * (hgt_tmp - hgt_dir(1))
        dlt_hgt(2)  =       (hgt_tmp - hgt_dir(2))
        dlt_hgt(3)  = fac * (hgt_tmp - hgt_dir(3))
        dlt_hgt(4)  =       (hgt_tmp - hgt_dir(4))
        dlt_hgt(5)  = fac * (hgt_tmp - hgt_dir(5))
        dlt_hgt(6)  =       (hgt_tmp - hgt_dir(6))
        dlt_hgt(7)  = fac * (hgt_tmp - hgt_dir(7))
        dlt_hgt(8)  =       (hgt_tmp - hgt_dir(8))
        k_loc      = maxloc (dlt_hgt)
        k_loc_min  = k_loc(1)
        dlt_hgt_min = dlt_hgt(k_loc(1))
        
	if (dlt_hgt_min > 1e-6) then
           flw(lon_idx,lat_idx) = k_loc(1)	  
        else if (dlt_hgt_min < -1e-6) then
           print *, "ERROR: found single cell depression at ", lon_idx, lat_idx
           stop
        else
           ! flat: set to -1 and process later
           !print *, "flat at ", lon_idx, lat_idx
           flw(lon_idx,lat_idx) = -1
        end if
     enddo
  enddo
end subroutine flw_drc

!-------------------------------------------------------------------
! RECURSIVE SUBROUTINE myqsort (xvec, order, left, right)
!-------------------------------------------------------------------
RECURSIVE SUBROUTINE myqsort (xvec, &     !I arry that be sorted elevation
                              order,&     !O arry that hold sorted elevation
                              left, &     !I elvation series
                              right)      !I elevation series
! Purpose: Sort elevation in order 
! myqsort() called by sort_elevation in geo_mdl.F90
  IMPLICIT NONE
  REAL,    intent(inout)::xvec(:)
  integer, intent(inout)::order(:)
  INTEGER, intent(in)   ::left,right
  REAL                  ::temp
  INTEGER               ::apu
  integer::last
  integer::idx
  integer               ::tempi

  IF(left >= right) RETURN
  apu         = (left+right)/2
  temp        = xvec(left)
  xvec(left)     = xvec(apu)
  xvec(apu)      = temp
  tempi       = order(left)
  order(left) = order(apu)
  order(apu)  = tempi
  last        = left
  DO idx=left+1,right
     IF(xvec(idx) > xvec(left)) THEN
        last        = last+1
        temp        = xvec(last)
        xvec(last)     = xvec(idx)
        xvec(idx)        = temp
        tempi       = order(last)
        order(last) = order(idx)
        order(idx)    = tempi
     ENDIF
  END DO
  temp        = xvec(left)
  xvec(left)     = xvec(last)
  xvec(last)     = temp
  tempi       = order(left)
  order(left) = order(last)
  order(last) = tempi
  CALL myqsort (xvec, order, left,  last-1)
  CALL myqsort (xvec, order, last+1, right)
  RETURN
END SUBROUTINE myqsort


!-------------------------------------------------------------------
! subroutine sort_elevation 
!-------------------------------------------------------------------
subroutine sort_elevation (lon_nbr, &    !I [nbr] Number of longitude
                           lat_nbr, &    !I [nbr] Number of latituden
                           oro,     &    !I Orographgy
                           hgt_sfc, &    !I [m] surface elevation
                           hgt_order,&   !O elevation order
                           hgt_map)      !O elevation map
! Purpose: Sort grid elevation 
! Sort_elevation() called by dst_src_flw in geo_mdl.F90
  implicit none
  integer, intent(in) ::lon_nbr
  integer, intent(in) ::lat_nbr
  logical, intent(in) ::oro(lon_nbr,lat_nbr)
  real,    intent(in) ::hgt_sfc(lon_nbr,lat_nbr)
  integer, intent(out)::hgt_order(lon_nbr*lat_nbr)
  integer, intent(out)::hgt_map(lon_nbr*lat_nbr,2)
  real                ::zvec(lon_nbr*lat_nbr)
  integer                 lon_idx, lat_idx, index
  
  hgt_order = -1
  hgt_map   = -1
  index = 1
  do lat_idx = 1, lat_nbr
     do lon_idx = 1, lon_nbr
        if (.not. oro(lon_idx,lat_idx)) cycle
        hgt_order(index) = index
        zvec(index)    = hgt_sfc(lon_idx,lat_idx)
        hgt_map(index,1) = lon_idx
        hgt_map(index,2) = lat_idx
        index          = index + 1
     end do
  end do
  index = index - 1
  if (index /= count(oro)) stop "ERROR: [sort_elevation] index /= count(oro)"
  call myqsort (zvec, hgt_order, 1, index)
end subroutine sort_elevation


!-------------------------------------------------------------------
! subroutine flw_acm_get
!-------------------------------------------------------------------
subroutine flw_acm_get (lon_nbr, &    !I [nbr] Number of longitude
                         lat_nbr, &    !I [nbr] Number of latituden
                         oro,     &    !I Orographgy
                         hgt_sfc, &    !I [m] surface elevation
                         hgt_order,&   !I elevation order
                         hgt_map,  &   !I elevation map
                         flw_dir,  &   !I flow direction
                         flw_acm,  &   !I/O [m2] flow accumulation
                         is_acm)       !O accumulation flag
! Purpose: Upatream area accumulation 
! flw_acm_get() called by dst_src_flw in geo_mdl.F90
  use dbg_mdl
  implicit none
  integer, intent(in)   ::lon_nbr
  integer, intent(in)   ::lat_nbr
  logical, intent(in)   ::oro(lon_nbr,lat_nbr)
  real,    intent(in)   ::hgt_sfc(lon_nbr,lat_nbr)
  integer, intent(in)   ::hgt_order (lon_nbr*lat_nbr)
  integer, intent(in)   ::hgt_map   (lon_nbr*lat_nbr,2)
  integer, intent(in)   ::flw_dir (lon_nbr,lat_nbr)
  real,    intent(inout)::flw_acm (lon_nbr,lat_nbr)
  logical, intent(out)  ::is_acm  (lon_nbr,lat_nbr)
!local
  integer  lon_idx,lat_idx,count_idx
  
  where (oro)
     is_acm  = .false.
  elsewhere
     is_acm  = .true.
     flw_acm = 0.0
  end where
  
  do count_idx = 1, count(oro)
     lon_idx = hgt_map(hgt_order(count_idx),1)
     lat_idx = hgt_map(hgt_order(count_idx),2)
     if (lon_idx < 1 .or. lon_idx > lon_nbr .or. lat_idx < 1 .or. lat_idx > lat_nbr) &
        stop "ERROR: [flw_acm_get] lon_idx or lat_idx out of bounds"
     if (.not. oro(lon_idx,lat_idx))                         stop "ERROR: [flw_acm_get] oro = 0"
     if (.not.(is_acm(lon_idx,lat_idx))) call accumulate_one (lon_nbr, lat_nbr, lon_idx, lat_idx, hgt_sfc, flw_dir, flw_acm, is_acm)
  end do
  
  where (.not. oro) flw_acm = 0.0

end subroutine flw_acm_get


!-------------------------------------------------------------------
! recursive subroutine accumulate_one
!-------------------------------------------------------------------
recursive subroutine accumulate_one (lon_nbr,   &  !I [nbr] number of longitude
                                   lat_nbr,   &  !I [nbr]number of latitude
                                   lon_idx,   &  !I longitude index
                                   lat_idx,   &  !I latitude index
                                   hgt_sfc,   &  !I [m] surface elevation
                                   flw_dir,   &  !I flow direction
                                   flw_acm,   &  !I/O [m2] flow accumulate
                                   is_acm)       !I/O accumulation flag
! Purpose: Upstream area accumulation 
! accumulate_one() called by flw_acm_get in geo_mdl.F90
    use dbg_mdl
    implicit none
    integer, intent(in)   ::lon_nbr
    integer, intent(in)   ::lat_nbr
    integer, intent(in)   ::lon_idx
    integer, intent(in)   ::lat_idx
    real,    intent(in)   ::hgt_sfc(lon_nbr,lat_nbr)
    integer, intent(in)   ::flw_dir (lon_nbr,lat_nbr)
    real,    intent(inout)::flw_acm (lon_nbr,lat_nbr)
    logical, intent(inout)::is_acm  (lon_nbr,lat_nbr)
!local
    integer  lon_neb_idx, lat_neb_idx
    integer  dir, opp_dir

    ! First check for special case ==>
    ! a) Is neighbor at same elevation?
    ! b) Does neighbor flw into i,j
    ! c) Has neighbor not yet been accumulated ?
    ! If all above true, then accumulate neighbor first (may need to do recursively)

!--------------------------------------------- all dirs

    do dir = 1, 8
       call get_neighbor (lon_nbr,  &   !I [nbr] number of longitude
                          lat_nbr,  &   !I [nbr] number of latitude
                          lon_idx,  &   !I longitude index
                          lat_idx,  &   !I latitude index
                          dir,      &   !I flow direction
                          lon_neb_idx, & !O longitude index of neighbor cell
                          lat_neb_idx)   !O latitude index of neighbor cell
       if (is_close(hgt_sfc(lon_idx,lat_idx),hgt_sfc(lon_neb_idx,lat_neb_idx))) then
          opp_dir = mod(dir+3,8)+1
          if (flw_dir(lon_neb_idx,lat_neb_idx)==opp_dir .and. flw_dir(lon_idx,lat_idx)/=dir) then
             if (.not.is_acm(lon_neb_idx,lat_neb_idx)) then
                !if (geo_dbg_lvl > 1) write(6,*) "...going upstream from (", i, j, ") to (", ii, jj, ")"
                call accumulate_one (lon_nbr,   &   !I [nbr] number of longitude
                                   lat_nbr,   &   !I [nbr] number of latitude
                                   lon_neb_idx,&  !I longitude index of neighbor cell
                                   lat_neb_idx, & !I latitude index of neighbor cell
                                   hgt_sfc,    &  !I [m] surface elevation
                                   flw_dir,    &  !I flow direction
                                   flw_acm,    &  !O flow accumulation
                                   is_acm)        !O flag for flow accumulation
             endif
          endif
       endif
    end do
    
!--------------------------------------------- confirm and dont accumulate single cell depressions

    do dir = 1, 8
       call get_neighbor (lon_nbr,  &   !I  [nbr] number of longitude
                          lat_nbr,  &   !I  [nbr] number of latitude
                          lon_idx,  &   !I  longitude index
                          lat_idx,  &   !I  latitude index
                          dir,      &   !I  flow direction
                          lon_neb_idx, & !O longitude index of neighbor cell
                          lat_neb_idx)   !O latitude index of neighbor cell
       if (is_close(hgt_sfc(lon_idx,lat_idx),hgt_sfc(lon_neb_idx,lat_neb_idx))) then
          opp_dir = mod(dir+3,8)+1
          if (flw_dir(lon_neb_idx,lat_neb_idx)==opp_dir .and. flw_dir(lon_idx,lat_idx)==dir) then
             !if (geo_dbg_lvl > 1) write(6,*) "... single cell depression at ", i, j
          endif
       endif
    end do
    
!--------------------------------------------- standard acmulation case

    ! Else standard acmulation case
    if (flw_dir(lon_idx,lat_idx)==0) then
       !if (geo_dbg_lvl > 1) write(6,*) "WARNING: [accumulate_one] unknown flw direction"
       is_acm(lon_idx,lat_idx)    = .true.
    else
       call get_neighbor (lon_nbr,  &   !I [nbr] number of logitude
                          lat_nbr,  &   !I [nbr] number of latitude
                          lon_idx,  &   !I logitude index
                          lat_idx,  &   !I latitude index
                          flw_dir(lon_idx,lat_idx),    &   !I flow direction
                          lon_neb_idx, & !O logitude index of neighbor cell
                          lat_neb_idx)   !O latitude index of neighbor cell
       flw_acm(lon_neb_idx,lat_neb_idx) = flw_acm(lon_neb_idx,lat_neb_idx) + flw_acm(lon_idx,lat_idx)
       is_acm(lon_idx,lat_idx)    = .true.
    endif
    
end subroutine accumulate_one

!-------------------------------------------------------------------
! subroutine chk_flw_1 
!-------------------------------------------------------------------
subroutine chk_flw_1 (lon_nbr,  &   !I [nbr] number of logitude
                      lat_nbr,  &   !I [nbr] number of latitude
                      oro,      &   !I orography
                      hgt_sfc,  &   !I [m] surface elevation
                      flw)          !I  flow direction
! Purpose: Check that all flows are downhill 
! chk_flw_1() called by dst_src_flw in geo_mdl.F90
  implicit none
  integer, intent(in)::lon_nbr
  integer, intent(in)::lat_nbr
  logical, intent(in)::oro (lon_nbr,lat_nbr)
  real,    intent(in)::hgt_sfc (lon_nbr,lat_nbr)
  integer, intent(in)::flw (lon_nbr,lat_nbr)
!local
  integer::lon_idx
  integer::lat_idx
  integer::dir
  integer::lon_neb_idx
  integer::lat_neb_idx
  integer::opp_dir
  
  do lat_idx = 1, lat_nbr
     do lon_idx = 1, lon_nbr
        if (.not. oro(lon_idx,lat_idx)) cycle
        if (flw(lon_idx,lat_idx) > 0) then
           dir = flw(lon_idx,lat_idx)
           call get_neighbor (lon_nbr,  &   !I  [nbr] number of logitude
                              lat_nbr,  &   !I  [nbr] number of latitude
                              lon_idx,  &   !I  logitude index
                              lat_idx,  &   !I  latitude index
                              dir,      &   !I  flow direction
                              lon_neb_idx, & !O longitude index of neighbor cell 
                              lat_neb_idx)   !O latitude index of neighbor cell
           if (hgt_sfc(lon_idx,lat_idx) < hgt_sfc(lon_neb_idx,lat_neb_idx)) then
              print *, "ERROR: found uphill flow ", lon_idx, lat_idx, hgt_sfc(lon_idx,lat_idx), &
                        lon_neb_idx, lat_neb_idx, hgt_sfc(lon_neb_idx,lat_neb_idx)
              stop
           end if
        end if
     enddo
  enddo
end subroutine chk_flw_1

!-------------------------------------------------------------------
! subroutine chk_flw_2 
!-------------------------------------------------------------------
subroutine chk_flw_2 (lon_nbr,  &   !I  [nbr] number of logitude
                      lat_nbr,  &   !I  [nbr] number of latitude
                      oro,      &   !I  orography
                      hgt_sfc,  &   !I  [m] surface elevation
                      flw)          !I  flow direction
! Purpose: check that all flows are downhill 
! chk_flw_2() called by dst_src_flw in geo_mdl.F90
  implicit none
  integer, intent(in)::lon_nbr
  integer, intent(in)::lat_nbr
  logical, intent(in)::oro (lon_nbr,lat_nbr)
  real,    intent(in)::hgt_sfc (lon_nbr,lat_nbr)
  integer, intent(in)::flw (lon_nbr,lat_nbr)
!local
  integer::lon_idx
  integer::lat_idx
  integer::dir
  integer::lon_neb_idx
  integer::lat_neb_idx
  integer::opp_dir
  
  do lat_idx = 1, lat_nbr
     do lon_idx = 1, lon_nbr
        if (.not. oro(lon_idx,lat_idx)) cycle
        if (flw(lon_idx,lat_idx) > 0) then
           dir = flw(lon_idx,lat_idx)
           opp_dir = mod(dir+3,8)+1
           call get_neighbor (lon_nbr,  &   !I [nbr] number of longitude
                              lat_nbr,  &   !I [nbr] number of latitude
                              lon_idx,  &   !I longitude index
                              lat_idx,  &   !I latitude index
                              dir,      &   !I flow direction
                              lon_neb_idx, & !O longitude index of neighbor cell
                              lat_neb_idx)   !O latitude index of neighbor cell 
           if (flw(lon_neb_idx,lat_neb_idx) == opp_dir) then
              print *, "ERROR: found circular flow ", lon_idx, lat_idx, flw(lon_idx,lat_idx), &
              hgt_sfc(lon_idx,lat_idx), lon_neb_idx, lat_neb_idx, flw(lon_neb_idx,lat_neb_idx),&
              hgt_sfc(lon_neb_idx,lat_neb_idx)
              stop
           end if
        end if
     enddo
  enddo
end subroutine chk_flw_2

!-------------------------------------------------------------------
! subroutine chk_flw 
!-------------------------------------------------------------------
subroutine chk_flw (lon_nbr,  &   !I  [nbr] number of longitude
                    lat_nbr,  &   !I  [nbr] number of latitude
                    oro,      &   !I  orography
                    hgt_sfc,  &   !I  [m] surface elevation
                    flw_dir)      !I  flow direction
! Purpose: check that all flows are downhill
! chk_flw() called by dst_src_flw in geo_mdl.F90
  use dbg_mdl
  implicit none
  integer, intent(in)::lon_nbr
  integer, intent(in)::lat_nbr
  logical, intent(in)::oro (lon_nbr,lat_nbr)
  real,    intent(in)::hgt_sfc (lon_nbr,lat_nbr)
  integer, intent(in)::flw_dir (lon_nbr,lat_nbr)
!local
  integer::lon_idx
  integer::lat_idx
  integer::dir
  integer::lon_neb_idx
  integer::lat_neb_idx
  integer::opp_dir

  do lat_idx = 2, lat_nbr-1
     do lon_idx = 1, lon_nbr
        if (.not. oro(lon_idx,lat_idx)) cycle
        if (flw_dir(lon_idx,lat_idx) > 0) then
           dir = flw_dir(lon_idx,lat_idx)
           call get_neighbor (lon_nbr,  &   !I [nbr] number of longitude
                              lat_nbr,  &   !I [nbr] number of latitude
                              lon_idx,  &   !I longitude index
                              lat_idx,  &   !I latitude index
                              dir,      &   !I flow direction
                              lon_neb_idx, & !O longitude index of neighbor cell
                              lat_neb_idx)   !O latitude index of neighbor cell
           opp_dir = mod(dir+3,8)+1
           if (flw_dir(lon_neb_idx,lat_neb_idx)==opp_dir) then
              write (6,*) "ERROR: found circular flow at (lon_idx,lat_idx), dir = ", lon_idx, lat_idx, dir
              print *, "lon_idx   = ", lon_idx
              print *, "lat_idx   = ", lat_idx
              print *, "dir = ", dir
              print *, "lon_neb_idx  = ", lon_neb_idx
              print *, "lat_neb_idx  = ", lat_neb_idx
              print *, "opp = ", opp_dir
              print *, "hgt_sfc   = ", hgt_sfc(lon_idx,lat_idx)
              print *, "hgt_sfc_nb  = ", hgt_sfc(lon_neb_idx,lat_neb_idx)
              print *, "oro = ", oro(lon_idx,lat_idx)
              print *, "oo  = ", oro(lon_neb_idx,lat_neb_idx)
              stop "ERROR: found circular flow"
           endif
           if (hgt_sfc(lon_idx,lat_idx) < hgt_sfc(lon_neb_idx,lat_neb_idx)) then
              write (6,*) "ERROR: flow uphill at (lon_idx,lat_idx), dir = ", lon_idx, lat_idx, dir
              stop "ERROR: flow uphill"
           endif
        endif
     end do
  end do
  
  write (6,*) "INFO: chk_flw ok"

end subroutine chk_flw

subroutine preprocess (lon_nbr, &    !I [nbr] number of longitude
                       lat_nbr, &    !I [nbr] number of latitude
                       oro,     &    !I orography
                       hgt_sfc)      !O [m] surface elevation
! Purpose: set some variables as reasonable values 
! preprocess() called by dst_src_flw() in geo_mdl.F90
  implicit none
  integer, intent(in)   ::lon_nbr, lat_nbr
  logical, intent(in)   ::oro(lon_nbr,lat_nbr)
  real,    intent(inout)::hgt_sfc(lon_nbr,lat_nbr)
  integer::lon_idx, lat_idx, count
  
  !where (oro(1,:)) h(1,:) = 6000.0
  !where (oro(m,:)) h(m,:) = 6000.0
  
  do lon_idx = 1, lon_nbr
     lat_idx = 1
     if (oro(lon_idx,lat_idx)) hgt_sfc(lon_idx,lat_idx) = 6000.0
     lat_idx = lat_nbr
     if (oro(lon_idx,lat_idx)) hgt_sfc(lon_idx,lat_idx) = 6000.0
  enddo

  count = 0
  do lat_idx = 1, lat_nbr
     do lon_idx = 1, lon_nbr
        if (oro(lon_idx,lat_idx) .and. hgt_sfc(lon_idx,lat_idx) <= 0.0) then
           hgt_sfc(lon_idx,lat_idx) = 1e-3
           count = count + 1
        end if
        if (.not. oro(lon_idx,lat_idx) .and. hgt_sfc(lon_idx,lat_idx) > 0.0) then
           hgt_sfc(lon_idx,lat_idx) = 0.0
           count = count + 1
        end if
     enddo
  enddo
  !print *, ">>> num fixed points = ", count
end subroutine preprocess

subroutine init (lon_nbr, &    !I     [nbr] number of longitude
                 lat_nbr, &    !I     [nbr] number of latitude
                 flw,     &    !I/O   flow direction
                 pointsTo,&    !I/O
                 sink,    &    !I/O   basin
                 bsn_enm, &    !I/O   drainage basin
                 internal,&    !I/O
                 dlt_hgt)      !O     [m] difference of surface elevation
! Purpose: initilization of some cariables 
! init called by dst_src_flw()  in dst_src_flw.F90
  implicit none
  integer, intent(in)   ::lon_nbr, lat_nbr
  integer, intent(inout)::flw(lon_nbr,lat_nbr)
  integer, intent(inout)::pointsTo(lon_nbr,lat_nbr)
  integer, intent(inout)::sink(lon_nbr,lat_nbr)
  integer, intent(inout)::bsn_enm(lon_nbr,lat_nbr)
  integer, intent(inout)::internal(lon_nbr,lat_nbr)
  real,    intent(out)  ::dlt_hgt(lon_nbr,lat_nbr)
  
  flw      = 0
  pointsTo = -1
  sink     = -1
  bsn_enm    = -1
  internal = -1
  dlt_hgt       = 0.0
  
end subroutine init

subroutine fillSingleCellDepression (lon_nbr, &    !I  [nbr] number of longitude
                                     lat_nbr, &    !I  [nbr] number of latitude
                                     oro,     &    !I  orography
                                     hgt_sfc)      !O  [m] surface elevation
! Purpose: set minimum elevation for some grid 
! fillSingleCellDepression called by dst_src_flw() in geo_mdl.F90
  implicit none
  integer, intent(in)   ::lon_nbr, lat_nbr
  logical, intent(in)   ::oro(lon_nbr,lat_nbr)
  real,    intent(inout)::hgt_sfc(lon_nbr,lat_nbr)
!local
  integer::lon_idx, lat_idx
  real   ::hgt_dir(8), hgt_min

  do lat_idx = 1, lat_nbr
     do lon_idx = 1, lon_nbr
        if (.not. oro(lon_idx,lat_idx)) cycle
        call get_cells (lon_nbr, &   !I [nbr] number of longitude
                        lat_nbr, &   !I [nbr] number of latitude
                        lon_idx ,&   !I longitude index
                        lat_idx, &   !I latitude index
                        hgt_sfc, &   !I [m] surface elevation
                        hgt_dir)     !O direction [1-8]
        hgt_min = minval(hgt_dir)
        if (hgt_sfc(lon_idx,lat_idx) < hgt_min) then
           !print *, "INFO: single cell depression at ", lon_idx, lat_idx
           hgt_sfc(lon_idx,lat_idx) = hgt_min
        endif
     enddo
  enddo
  
end subroutine fillSingleCellDepression

function notLoop(lon_nbr,    &  !I [nbr] number of longitude
                 lat_nbr,    &  !I [nbr] number of latitude 
                 lon_idx_in, &  !I longitude index
                 lat_idx_in, &  !I latitude index
                 oro,        &  !I orography
                 flw)           !I/O flow direction
! Purpose: get neighbor lat and lon indeces
! notLoop() called by fix_flw_drc() in geo_mdl.F90
  implicit none
  logical:: notLoop
  integer, intent(in)::lon_nbr, lat_nbr, lon_idx_in, lat_idx_in
  logical, intent(in)::oro(lon_nbr,lat_nbr)
  integer, intent(in)::flw(lon_nbr,lat_nbr)
  integer::lon_idx, lat_idx, lon_neb_idx, lat_neb_idx, dir, count
  
  count = 0
  lon_idx = lon_idx_in
  lat_idx = lat_idx_in
  do while (oro(lon_idx,lat_idx))
     dir = flw(lon_idx,lat_idx)
     if (dir < 1) exit
     call get_neighbor(lon_nbr,  &   !I [nbr] number of Longitude
                       lat_nbr,  &   !I [nbr] number of latitude
                       lon_idx,  &   !I longitude index
                       lat_idx,  &   !I latitude index
                       dir,      &   !O flow direction
                       lon_neb_idx, &  !O longitude index for neighbor cell
                       lat_neb_idx)    !O latitude index for neighbor cell
     lon_idx = lon_neb_idx
     lat_idx = lat_neb_idx
     count = count + 1
     ! fxm 20030720: What does this mean?
!     if (count > 100) then
!       print *, "ERROR: infinite loop in notLoop"
!       stop
!    end if
    ! This arbitrary limit is not necessary, since the loops will be 
    ! automatically setup by code according to the resolutions.  
  end do
  notLoop = .true.
end function notLoop


function unresolvedFlow(lon_nbr, &     !I [nbr] number of longitude
                        lat_nbr, &     !I [nbr] number of latitude
                        oro,     &     !I  orography
                        flw)           !O  flow direction
! Purpose: get numbers of unresolvedFlow
! unresolvedFlow() called by dst_src_flw() in geo_mdl.F90
  implicit none
  integer::unresolvedFlow
  integer, intent(in)::lon_nbr, lat_nbr
  logical, intent(in)::oro(lon_nbr,lat_nbr)
  integer, intent(in)::flw(lon_nbr, lat_nbr)
  integer::lon_idx, lat_idx

  unresolvedFlow = 0

  do lat_idx = 1, lat_nbr
     do lon_idx = 1, lon_nbr
        if (oro(lon_idx,lat_idx) .and. flw(lon_idx,lat_idx) == -1) then
           unresolvedFlow = unresolvedFlow + 1
        end if
     enddo
  enddo
end function unresolvedFlow
  
subroutine fixFlowDir (lon_nbr, &   !I   [nbr] number of longitude
                       lat_nbr, &   !I   [nbr] number of latitude
                       oro,     &   !I   orography
                       hgt_sfc, &   !I   [m] surface elevation
                       flw)         !I/O flow direction

! Purpose: check flat point 
! fixFlowDir called by dst_src_flw() in geo_mdl.F90

  implicit none
  integer, intent(in)   ::lon_nbr, lat_nbr
  logical, intent(in)   ::oro(lon_nbr,lat_nbr)
  real,    intent(in)   ::hgt_sfc(lon_nbr,lat_nbr)
  integer, intent(inout)::flw(lon_nbr,lat_nbr)
  integer::count

  call fix_flw_drc(lon_nbr, &   !I [nbr] number of longitude
                         lat_nbr, &   !I [nbr] number of latitude
                         oro,     &   !I orography
                         hgt_sfc, &   !I [m] surface elevation
                         flw,     &   !I/O flow direction
                         count)       !O [nbr] count number
  !print *, "count = ", count
  call chk_flw_1 (lon_nbr, &   !I  [nbr] number of longitude
                  lat_nbr, &   !I  [nbr] number of latitude
                  oro,     &   !I  orography
                  hgt_sfc, &   !I  [m] surface elevation
                  flw)         !I/O flow direction
  call chk_flw_2 (lon_nbr, &   !I  [nbr] number of longitude
                  lat_nbr, &   !I  [nbr] number of latitude
                  oro,     &   !I  orography
                  hgt_sfc, &   !I  [m] surface elevation
                  flw)         !I/O flow direction
  
  do while (count > 0) 
   call fix_flw_drc (lon_nbr, &   !I  [nbr] number of longitude
                           lat_nbr, &   !I  [nbr] number of latitude
                           oro,     &   !I  orography
                           hgt_sfc, &   !I  [m] surface elevation
                           flw,     &   !I/O flow direction
                           count)       !O count number
     !print *, "count = ", count
   call chk_flw_1 (lon_nbr, &   !I  [nbr] number of longitude
                   lat_nbr, &   !I  [nbr] number of latitude
                   oro,     &   !I  orography
                   hgt_sfc, &   !I  [m] surface elevation
                   flw)         !I/O flow direction
   call chk_flw_2 (lon_nbr, &   !I  [nbr] number of longitude
                   lat_nbr, &   !I  [nbr] number of latitude
                   oro,     &   !I  orography
                   hgt_sfc, &   !I  [m] surface elevation
                   flw)         !I/O flow direction
  end do

end subroutine fixFlowDir

subroutine fix_flw_drc (lon_nbr,  &   !I  [nbr] number of longitud
                              lat_nbr,  &   !I  [nbr] number of latitude
                              oro,      &   !I  orography
                              hgt_sfc,  &   !I  [m] surface elevatio
                              flw,      &   !I/O flow direction
                              count)        !O  [nbr] count number
! Purpose: check all flat points and sees if there is a path to the surf
! fix_flw_drc called by fixFlowDir()  in geo_mdl.F90
  implicit none
  integer, intent(in)::lon_nbr
  integer, intent(in)::lat_nbr
  logical, intent(in)::oro (lon_nbr,lat_nbr)
  real,    intent(in)::hgt_sfc (lon_nbr,lat_nbr)
  integer, intent(inout)::flw (lon_nbr,lat_nbr)
  integer, intent(out)  ::count
!Local
  integer::lon_idx
  integer::lat_idx
  integer::dir
  integer::lon_neb_idx
  integer::lat_neb_idx

  ! finds all flat points and sees if there is a path to the surf
  count = 0
  do lat_idx = 1, lat_nbr
     do lon_idx = 1, lon_nbr
        if (.not. oro(lon_idx,lat_idx)) cycle
        
        if (flw(lon_idx,lat_idx) == -1) then
           do dir = 1, 8
              
           call get_neighbor (lon_nbr,  &   !I  [nbr] number of longitude
                              lat_nbr,  &   !I  [nbr] number of latitude
                              lon_idx,  &   !I  longitude index
                              lat_idx,  &   !I  latitude index
                              dir,      &   !I  flow direction
                              lon_neb_idx, & !O longitude index of neighbor cell
                              lat_neb_idx)   !O latitude index of neighbor cell
              
              if (flw(lon_neb_idx,lat_neb_idx) /= -1) then
                 
                 if (flw(lon_neb_idx,lat_neb_idx) /= mod(dir+3,8)+1) then
                    if (hgt_sfc(lon_idx,lat_idx) >= hgt_sfc(lon_neb_idx,lat_neb_idx)) then
                       
                       ! flows to ocean

                       if (notLoop(lon_nbr, &    !I  [nbr] number of longitude
                                   lat_nbr, &    !I  [nbr] number of latitude
                                   lon_neb_idx, &    !I longitude index of neighbor cell
                                   lat_neb_idx, &    !I latitude index of neighbor cell
                                   oro,     &    !I  orography
                                   flw)) then    !O  flow direction
                          flw(lon_idx,lat_idx) = dir
                          count = count + 1
                          !print *, "fixed lon_idx, lat_idx, hgt_sfc, dir = ", lon_idx, lat_idx, hgt_sfc(lon_idx,lat_idx), dir
                       end if
                    end if
                 end if
              end if
           enddo
        end if
     enddo
  enddo
end subroutine fix_flw_drc

subroutine getPointsTo (lon_nbr,  &  !I  [nbr] number of longitude
                        lat_nbr,  &  !I  [nbr] number of latitude
                        oro,      &  !I  orography
                        flw,      &  !I  flow direction
                        pointsTo)    !O
! Purpose: set ID to each point 
! getPointsTo() called by dst_src_flw() in geo_mdl.F90
  implicit none
  integer, intent(in)   ::lon_nbr,lat_nbr
  logical, intent(in)   ::oro(lon_nbr,lat_nbr)
  integer, intent(in)   ::flw(lon_nbr,lat_nbr)
  integer, intent(inout)::pointsTo(lon_nbr,lat_nbr)
!local
  integer::lon_idx, lat_idx, lon_neb_idx, lat_neb_idx, dir
  
  do lat_idx = 1, lat_nbr
     do lon_idx = 1, lon_nbr
        if (oro(lon_idx,lat_idx)) then
           ! see where cell (i,j) points to
           dir = flw(lon_idx,lat_idx)
             if (dir > 0) then
              call get_neighbor (lon_nbr,  &   !I [nbr] number of Longitude
                                 lat_nbr,  &   !I [nbr] number of latitude
                                 lon_idx,  &   !I longitude index
                                 lat_idx,  &   !I latitude index
                                 dir,      &   !I flow direction
                                 lon_neb_idx, &  !O longitude index for neighbor cell
                                 lat_neb_idx)    !O latitude index for neighbor cell)
                  pointsTo(lon_idx,lat_idx) = (lat_neb_idx-1)*lon_nbr + (lon_neb_idx-1)
              !print *, lon_idx, lat_idx, dir, lon_neb_idx, lat_neb_idx, pointsTo(lon_idx,lat_idx)
           end if
        end if
     enddo
  enddo
end subroutine getPointsTo

function getSink(lon_nbr, &  !I  [nbr] number of Longitude
                 lat_nbr, &  !I  [nbr] number of latitude
                 pointsTo,&  !I
                 sink)       !I/O  basin
! Purpose: compute sink number 
! getSink called by dst_src_flw() in geo_mdl.F90
  implicit none
  integer::getSink
  integer, intent(in) ::lon_nbr, lat_nbr
  integer, intent(in) ::pointsTo(lon_nbr,lat_nbr)
  integer, intent(inout)::sink(lon_nbr,lat_nbr)
!local
  integer::lon_idx, lat_idx, lon_neb_idx, lat_neb_idx, point, sinkId
  
  sinkId = 0
  
  do lat_idx = 1, lat_nbr
     do lon_idx = 1, lon_nbr
        point = pointsTo(lon_idx,lat_idx)
        if (point /= -1) then
           do while (point /= -1)
              lon_neb_idx = mod(point,lon_nbr) + 1
              lat_neb_idx = point/lon_nbr + 1
              point = pointsTo(lon_neb_idx,lat_neb_idx)
              !print *, lon_neb_idx, lat_neb_idx, point
           enddo
           sink(lon_neb_idx,lat_neb_idx) = sinkId
           sinkId = sinkId + 1
        end if
     enddo
  enddo
  getSink = sinkId
end function getSink

subroutine mergeSink (lon_nbr,  &   !I  [nbr] number of longitude
                      lat_nbr,  &   !I  [nbr] number of latitude
                      oro,      &   !I  orography
                      hgt,      &   !I  [m] surface elevation
                      flw,      &   !I  flow direction
                      sink)         !I/O drainage basin
! Purpose: merge same elevation grid to same sink
! mergeSink() called by dst_src_flw() in geo_mdl.F90
  implicit none
  integer, intent(in)   ::lon_nbr, lat_nbr
  logical, intent(in)   ::oro(lon_nbr,lat_nbr)
  integer, intent(in)   ::flw(lon_nbr,lat_nbr)
  real,    intent(in)   ::hgt(lon_nbr,lat_nbr)
  integer, intent(inout)::sink(lon_nbr,lat_nbr)
!local
  integer::lon_idx, lat_idx, lon_neb_idx, lat_neb_idx, dir
  
  do lat_idx = 1, lat_nbr
     do lon_idx = 1, lon_nbr
        if (oro(lon_idx,lat_idx) .and. flw(lon_idx,lat_idx) == -1) then
           do dir = 1, 8
              call get_neighbor (lon_nbr,  &   !I [nbr] number of Longitude
                                 lat_nbr,  &   !I [nbr] number of latitude
                                 lon_idx,  &   !I longitude index
                                 lat_idx,  &   !I latitude index
                                 dir,      &   !I/O flow direction
                                 lon_neb_idx, &  !O longitude index for neighbor cell
                                 lat_neb_idx)    !O latitude index for neighbor cell
              if (oro(lon_neb_idx,lat_neb_idx) .and. flw(lon_neb_idx,lat_neb_idx) == -1) then 
                 if (hgt(lon_idx,lat_idx) == hgt(lon_neb_idx,lat_neb_idx)) then
                    if (sink(lon_idx,lat_idx) /= sink(lon_neb_idx,lat_neb_idx)) then
                       !print *, "merge ", lon_idx, lat_idx, sink(lon_idx,lat_idx), lon_neb_idx, lat_neb_idx, sink(lon_neb_idx,lat_neb_idx)
                       sink(lon_idx,lat_idx) = sink(lon_neb_idx,lat_neb_idx)
                    end if
                 end if
              end if
           enddo
        end if
     enddo
  enddo
end subroutine mergeSink

subroutine getBasin (lon_nbr,  &   !I  [nbr] number of longitude
                     lat_nbr,  &   !I  [nbr] number of latitude
                     pointsTo, &   !I
                     sink,     &   !I  basin
                     bsn_enm)      !I/O basin
! Purpose: Get Basin ID 
! getBasin() called by dst_src_flw() in geo_mdl.F90
  implicit none
  integer, intent(in)   ::lon_nbr, lat_nbr
  integer, intent(in)   ::pointsTo(lon_nbr,lat_nbr)
  integer, intent(in)   ::sink(lon_nbr, lat_nbr)
  integer, intent(inout)::bsn_enm(lon_nbr, lat_nbr)
!Local
  integer::lon_idx, lat_idx, lon_neb_idx, lat_neb_idx, point

  do lat_idx = 1, lat_nbr
     do lon_idx = 1, lon_nbr
        point = pointsTo(lon_idx,lat_idx)
        if (point /= -1) then
           do while (point /= -1)
              lon_neb_idx = mod(point,lon_nbr) + 1
              lat_neb_idx = point/lon_nbr + 1
              point = pointsTo(lon_neb_idx,lat_neb_idx)
           enddo
           bsn_enm(lon_idx,lat_idx) = sink(lon_neb_idx,lat_neb_idx)
        end if
     enddo
  enddo
  
  do lat_idx = 1, lat_nbr
     do lon_idx = 1, lon_nbr
        if (sink(lon_idx,lat_idx) /= -1) then
           bsn_enm(lon_idx,lat_idx) = sink(lon_idx,lat_idx)
        end if
     enddo
  enddo
end subroutine getBasin

subroutine getInternal(lon_nbr, &  !I  [nbr] number of longitude
                       lat_nbr, &  !I  [nbr] number of latitude
                       oro,     &  !I  orography
                       sink,    &  !I  basin
                       internal)   !I/O internal basin
! Purpose: find internal basin
! getInternal called by dst_src_flw() in geo_mdl.F90
  implicit none
  integer, intent(in)   ::lon_nbr, lat_nbr
  logical, intent(in)   ::oro(lon_nbr,lat_nbr)
  integer, intent(in)   ::sink(lon_nbr,lat_nbr)
  integer, intent(inout)::internal(lon_nbr,lat_nbr)
!Local
  integer::lon_idx, lat_idx

  do lat_idx = 1, lat_nbr
     do lon_idx = 1, lon_nbr
        if (sink(lon_idx,lat_idx) /= -1 .and. oro(lon_idx,lat_idx)) then
           internal(lon_idx,lat_idx) = sink(lon_idx,lat_idx)
        end if
     enddo
  enddo
end subroutine getInternal

function countInternal(lon_nbr,  &    !I [nbr] number of longitude
                       lat_nbr,  &    !I [nbr] number of latitude
                       numBasins,&    !I [nbr] number of basin
                       internal)      !I internal basin
! Purpose: count internal basin number 
! countInternal() called by dst_src_flw() in geo_mdl.F90
  implicit none
  integer::countInternal
  integer, intent(in)::lon_nbr, lat_nbr, numBasins
  integer, intent(in)::internal(lon_nbr,lat_nbr)
  integer::num_idx

  countInternal = 0
  
  do num_idx = 1, numBasins
     if (isInternal(lon_nbr,lat_nbr,num_idx,internal)) countInternal = countInternal+1
  enddo
end function countInternal

function isInternal(lon_nbr, &      !  longitude number
                    lat_nbr, &      !  latitude number
                    num_idx, &      !number of internal
                    internal)       ! internal basin
! Purpose: check if it is internal basin 
! isInternal() called by countInternal function in geo_mdl.F90
  implicit none
  logical::isInternal
  integer, intent(in)::lon_nbr, lat_nbr, num_idx
  integer, intent(in)::internal(lon_nbr,lat_nbr)
  integer::lon_idx, lat_idx

  isInternal = .false.
  do lat_idx = 1, lat_nbr
     do lon_idx = 1, lon_nbr
        if (internal(lon_idx, lat_idx) == num_idx) then
           isInternal = .true.
           exit
        end if
     enddo
  enddo
end function isInternal

function lon2i(lon_nbr, lon)
!Purpose: compute longitude index
!lon2i() called by isTrueInternal2 in geo_mdl.F90
  implicit none
  integer::lon2i
  integer, intent(in)::lon_nbr
  real, intent(in)::lon
  real::lon_tmp
  lon_tmp=lon
  if (lon_tmp < 0.0) lon_tmp = lon_tmp + 360.0
  lon2i = nint(lon_tmp*lon_nbr/360.0) + 1
end function lon2i

function lat2j(lat_nbr, lat)
!Purpose: compute latitude index
!lat2i() called by isTrueInternal2 in geo_mdl.F90
  implicit none
  integer::lat2j
  integer, intent(in)::lat_nbr
  real, intent(in)::lat
  lat2j = nint((lat+90.0)*lat_nbr/180.0) + 1
end function lat2j

function isTrueInternal2 (lon_nbr,   &   !I  [nbr] number of longitude
                          lat_nbr,   &   !I  [nbr] number of latitude
                          num_idx,   &   !I  [nbr] number of coordinator 
                          bsn_enm)       !I  drainage basin
! Purpose: set logical if it is internal basin 
! isTrueInternal2() called by fillBasin() in geo_mdl.F9
  implicit none
  integer,parameter::ntr_bsn_nbr=27 ! [nbr] Number of internal basins
  logical::isTrueInternal2
  integer, intent(in)::lon_nbr, lat_nbr, num_idx
  integer, intent(in)::bsn_enm(lon_nbr,lat_nbr)
  integer::lon_idx, lat_idx, bsn_idx
  real,parameter::inLon(ntr_bsn_nbr)=(/ &
       14.45,   35.13,   31.36,   36.43,   -5.83,   32.06,   36.04,   25.84,   59.84,   76.07, &
       77.39,   90.38,   92.11,   89.30,  140.04,  140.24,   50.11,   37.45, -118.02, -115.83, &
       -112.54, -117.83,  -67.49,  -67.69,  -68.27,  -68.03,  -67.45 /)
  real,parameter::inLat(ntr_bsn_nbr)=(/ &
       21.09,   -3.59,   22.87,    2.91,   21.09,   -7.83,    1.75,  -20.82,   44.86,   46.69, &
       42.41,   39.95,   42.70,   42.68,  -28.55,  -29.47,   43.02,   31.33,   39.80,   33.31, &
       41.11,   36.70,  -17.87,  -24.85,  -23.53,  -26.25,  -26.54 /) 
  
  isTrueInternal2 = .false.
  
  do lat_idx = 1, lat_nbr
     do lon_idx = 1, lon_nbr
        do bsn_idx = 1, ntr_bsn_nbr
           if (bsn_enm(lon_idx,lat_idx)==num_idx .and. lon_idx==lon2i(lon_nbr,inLon(bsn_idx)) .and. &
             lat_idx==lat2j(lat_nbr,inLat(bsn_idx))) then
              !print *, "xxxxxxxxxxxxxx internal bsn_enm, num_idx, lon_idx, lat_idx", num_idx, lon_idx, lat_idx
              isTrueInternal2 = .true.
              exit
           end if
        enddo
     enddo
  enddo
end function isTrueInternal2

subroutine fillBasin(lon_nbr,  &   !I  [nbr] number of longitude
                     lat_nbr,  &   !I  [nbr] number of latitude
                     numBasins,&   !I  [nbr] number of basin
                     oro,      &   !I  orography
                     hgt,      &   !I  [m] surface elevation
                     bsn_enm,  &   !I  drainage basin
                     internal, &   !I  internal basin
                     dlt_hgt)      !O  [m] difference of surface elevations
! Purpose: adjust elevation of neighbor grid 
! fillBasin() called by dst_src_flw() in geo_mdl.F90
  implicit none
  integer, intent(in)   ::lon_nbr, lat_nbr, numBasins
  logical, intent(in)   ::oro(lon_nbr,lat_nbr)
  real,    intent(in)   ::hgt(lon_nbr,lat_nbr)
  integer, intent(in)   ::bsn_enm(lon_nbr,lat_nbr)
  integer, intent(in)   ::internal(lon_nbr,lat_nbr)
  real,    intent(inout)::dlt_hgt(lon_nbr,lat_nbr)
  integer::lon_idx, lat_idx
  integer::lon_neb_idx, lat_neb_idx
  integer::num_idx, dir
  real   ::hgt_pour

  do num_idx = 1, numBasins
     if (isInternal(lon_nbr,lat_nbr,num_idx,internal) .and. .not.(isTrueInternal2(lon_nbr,lat_nbr,num_idx,bsn_enm))) then
        
        hgt_pour  = 1e6
        
        do lat_idx = 1, lat_nbr
           do lon_idx = 1, lon_nbr
              if (bsn_enm(lon_idx,lat_idx) == num_idx) then
                 do dir = 1, 4 !only pour in N, W, W, E directions
                    call get_neighbor (lon_nbr, &  !I [nbr] number of longitude
                                       lat_nbr, &  !I [nbr] number of latitude
                                       lon_idx, &  !I longitude index
                                       lat_idx, &  !I latitude index
                                       2*dir,   &  !I flow direction
                                       lon_neb_idx, & !I longitude index of neighbor cell
                                       lat_neb_idx)   !I latitude indx of neighbor cell
                    if (oro(lon_neb_idx,lat_neb_idx) .and. bsn_enm(lon_neb_idx,lat_neb_idx) /= num_idx) then
                       if (hgt(lon_neb_idx,lat_neb_idx) < hgt_pour) hgt_pour = hgt(lon_neb_idx,lat_neb_idx)
                    end if
                 enddo
              end if
           enddo
        enddo
        
        !print *, "num_idx, hgt_pour = ", num_idx, hgt_pour

        ! create dh
        do lat_idx = 1, lat_nbr
           do lon_idx = 1, lon_nbr
              if (bsn_enm(lon_idx,lat_idx) == num_idx) then
                 if (hgt(lon_idx,lat_idx) < hgt_pour) then
                    if (dlt_hgt(lon_idx,lat_idx) > 0.0) then
                       !print *, "dh: i, j, dh, h_pour ", i, j, dh(i,j), h_pour
                       if (dlt_hgt(lon_idx,lat_idx) > hgt_pour) dlt_hgt(lon_idx,lat_idx) = hgt_pour
                       !if (dh(i,j) < h_pour) dh(i,j) = h_pour less conservative
                    else
                       dlt_hgt(lon_idx,lat_idx) = hgt_pour
                    end if
                 end if
              end if
           enddo
        enddo
     end if
  end do
end subroutine fillBasin

subroutine adjustHeight(lon_nbr,   &   !I  [nbr] number of longitud
                        lat_nbr,   &   !I  [nbr] number of latitude
                        hgt,       &   !I/O [m] surface elevation
                        dlt_hgt)       !I  [m] difference of elevations
! Purpose: adjust neighbor height 
! adjustHeight() called by dst_src_flw() in geo_mdl.F9o
  implicit none
  integer, intent(in)   ::lon_nbr, lat_nbr
  real,    intent(inout)::hgt(lon_nbr,lat_nbr)
  real,    intent(in)   ::dlt_hgt(lon_nbr,lat_nbr)
  integer::lon_idx, lat_idx

  do lat_idx = 1, lat_nbr
     do lon_idx = 1, lon_nbr
        if (dlt_hgt(lon_idx,lat_idx) > 0.0) then
           
           if (dlt_hgt(lon_idx,lat_idx) < hgt(lon_idx,lat_idx)) then
              print *, "ERROR: lon_idx,lat_idx, hgt, dlt_hgt = ", lon_idx, lat_idx, hgt(lon_idx,lat_idx), dlt_hgt(lon_idx,lat_idx)
              stop
           end if
           
           !print *, "raising: i, j, h_old, h_new", i, j, hgt(i,j), dh(i,j)
           hgt(lon_idx,lat_idx) = dlt_hgt(lon_idx,lat_idx)

        end if
     end do
  end do
end subroutine adjustHeight


!-------------------------------------------------------------------
! subroutine write_matlab 
!-------------------------------------------------------------------
subroutine write_matlab (lon_nbr,    &     !I  [nbr] number of longitude
                         lat_nbr,    &     !I  [nbr] number of latitude
                         dat,        &     !I/O data in the file
                         filename,   &     !I/O file name
                         label)            !I  data label
  use dbg_mdl
  implicit none
  integer,          intent(in)::lon_nbr, lat_nbr
  real,             intent(in)::dat(lon_nbr,lat_nbr)
  character(len=*), intent(in)::filename
  character(len=*), intent(in)::label
  integer, parameter          ::fnum = 12
  integer                     ::OpenStatus, lon_idx, lat_idx

  open (unit=fnum, file=filename, status="unknown", iostat=OpenStatus)
  if (OpenStatus > 0) stop "ERROR: [write_matlab] cannot open file for writing"
  write (fnum, *) label
  do lat_idx = 1, lat_nbr
     write (fnum, '(1000(E13.6, 1X))') (dat(lon_idx,lat_idx), lon_idx = 1, lon_nbr)
  end do
  write (fnum, *) "];"
  close (fnum)
  if (geo_dbg_lvl > 0) write (6,*) "[write_matlab] wrote file ", filename
end subroutine write_matlab

!-------------------------------------------------------------------
! subroutine write_matlab_dir 
!-------------------------------------------------------------------
subroutine write_matlab_dir (lon_nbr,  &    !I  [nbr] number of longitude
                             lat_nbr,  &    !I  [nbr] number of latitude
                             hgt_sfc,  &    !I  [m] surface elevation
                             flw_dir,  &    !I/O flow direction
                             flw_acm)       !I/O flow accumulation
! Purpose: write out flow direction, and upstream areas 
! write_matlab_dir() called by dst_src_flw() in geo_mdl.F90
  implicit none
  integer, intent(in)::lon_nbr, lat_nbr
  real,    intent(in)::hgt_sfc(lon_nbr,lat_nbr)
  integer, intent(in)::flw_dir(lon_nbr,lat_nbr)
  real,    intent(in)::flw_acm(lon_nbr,lat_nbr)
  character          ::filename*32	
  character          ::label*32
  real               ::flw_dir_x(lon_nbr,lat_nbr)
  real               ::flw_dir_y(lon_nbr,lat_nbr)
  real               ::flw_acm_real(lon_nbr,lat_nbr)
  integer            ::lon_idx, lat_idx, dir

  filename="hgt_sfc_data.m"
  label="hgt_sfc = ["
  call write_matlab (lon_nbr, lat_nbr, hgt_sfc, filename, label)
  
  !------------
  ! 5 | 4 | 3 |
  !------------
  ! 6 | x | 2 |
  !------------
  ! 7 | 8 | 1 |
  !------------

  do lat_idx = 1, lat_nbr
     do lon_idx = 1, lon_nbr
        flw_acm_real(lon_idx,lat_idx) = flw_acm(lon_idx,lat_idx)
        dir               = flw_dir(lon_idx,lat_idx)
        if      (dir==1) then
           flw_dir_x(lon_idx,lat_idx) =  0.707
           flw_dir_y(lon_idx,lat_idx) = -0.707
        else if (dir==2) then
           flw_dir_x(lon_idx,lat_idx) =  1
           flw_dir_y(lon_idx,lat_idx) =  0
        else if (dir==3) then
           flw_dir_x(lon_idx,lat_idx) =  0.707
           flw_dir_y(lon_idx,lat_idx) =  0.707
        else if (dir==4) then
           flw_dir_x(lon_idx,lat_idx) =  0
           flw_dir_y(lon_idx,lat_idx) =  1
        else if (dir==5) then
           flw_dir_x(lon_idx,lat_idx) = -0.707
           flw_dir_y(lon_idx,lat_idx) =  0.707
        else if (dir==6) then
           flw_dir_x(lon_idx,lat_idx) = -1
           flw_dir_y(lon_idx,lat_idx) =  0
        else if (dir==7) then
           flw_dir_x(lon_idx,lat_idx) = -0.707
           flw_dir_y(lon_idx,lat_idx) = -0.707
        else if (dir==8) then
           flw_dir_x(lon_idx,lat_idx) =  0
           flw_dir_y(lon_idx,lat_idx) = -1
        else
           flw_dir_x(lon_idx,lat_idx) =  0
           flw_dir_y(lon_idx,lat_idx) =  0
        endif
     end do
  end do
  
  filename="flw_dir_x_data.m"
  label="flw_dir_x = ["
  call write_matlab (lon_nbr, lat_nbr, flw_dir_x, filename, label)
  
  filename="flw_dir_y_data.m"
  label="flw_dir_y = ["
  call write_matlab (lon_nbr, lat_nbr, flw_dir_y, filename, label)
  
! filename="flw_acm_data.m"
! label="flw_acm = ["
! call write_matlab (lon_nbr, lat_nbr, flw_acm_real, filename, label)
  
end subroutine write_matlab_dir

subroutine bsn_sz_get (lon_nbr,  &   !I [nbr] number of Longitude
                         lat_nbr,  &   !I [nbr] number of latitude
                         oro,      &   !I  Orography flag
                         numBasins, & !I [nbr] number of Basins
                         sfc_area,  &  !I [m2] grid area
                         bsn_enm,  &   !I [nbr] Basin ID
                         bsn_sz)     !O flow direction
! Purpose: compute drainage basin areas
! bsn_sz_get() called by dst_src_flw() in geo_mdl.F90
  implicit none
  integer, intent(in)   ::lon_nbr     ! [nbr] Number of longitude
  integer, intent(in)   ::lat_nbr     ! [nbr] Number of latitude
  integer, intent(in)   ::numBasins  ! [nbr] number of basins
  real,    intent(in)   ::sfc_area (lon_nbr,lat_nbr) ! [m2] surface area
  integer,    intent(in)   ::bsn_enm (lon_nbr,lat_nbr) ! [nbr] Basin ID
!  real,    intent(out)  ::bsn_sz (numBasins) ! [m2] Basin size
   real,    intent(out)  ::bsn_sz (lon_nbr*lat_nbr) ! [m2] Basin size
   logical, intent(in) ::oro(lon_nbr,lat_nbr) !  Orography flag

  integer  bsn_idx, lat_idx, lon_idx

  bsn_sz = 0.

  do  bsn_idx = 1, numBasins
     do lon_idx = 1, lon_nbr
        do lat_idx = 1, lat_nbr
            if (oro(lon_idx,lat_idx)) then
                if (bsn_enm(lon_idx,lat_idx).eq.bsn_idx) then
                   bsn_sz(bsn_idx)=bsn_sz(bsn_idx)+sfc_area(lon_idx,lat_idx)
                end if
           end if
        enddo
      enddo
  enddo 

! do  bsn_idx = 1, numBasins
!    if (bsn_sz(bsn_idx).ne.0.) then
!        write (999,'(I8,E12.4)') bsn_idx, bsn_sz(bsn_idx)
!    end if
! enddo 
end subroutine bsn_sz_get

end module geo_mdl ! [mdl] Geomorphology and flow accumulation
