! $Id$ -*-f90-*-

! Purpose: Process maps and identify aerosol sources

! Copyright (C) 1997--2002 Charlie Zender
! Portions are copyright by their respective contributors
! This software is distributed under the terms of the GNU General Public License
! See http://www.gnu.org/copyleft/gpl.html for full license text
! The original author of this software, Charlie Zender, seeks to improve
! it with your suggestions, contributions, bug-reports, and patches.
! Charlie Zender <zender at uci dot edu>
! Department of Earth System Science
! University of California at Irvine
! Irvine, CA 92697-3100
  
! Usage: 
! use src_prc ! [mdl] Source processing, identification

module src_prc ! [mdl] Source processing, identification
  implicit none
  private ! [stt] Symbols are private unless individually qualified as public
  public::aod_nmd_get ! [sbr] Process non-mineral dust (nmd) sources
  public::rdb_fct_tpg_GCT01 ! [sbr] Compute source efficiency factor from topography
  public::src_id ! [sbr] Find and parameterize dust sources
  public::oro_is_ocn ! [fnc] Point is > 50% ocean
  public::oro_is_lnd ! [fnc] Point is > 50% land
  public::cst_tpr ! [sbr] Taper erodibility factors along coasts
  public::rgn_msk ! [sbr] Mask region
  public::rgn_xcl ! [sbr] Zero field based on location
 
contains
  
  logical function oro_is_ocn(oro_val)
    ! Purpose: True if > 50% ocean
    real,intent(in)::oro_val ! [frc] Orography
    oro_is_ocn=nint(oro_val)==0
  end function oro_is_ocn
  
  logical function oro_is_lnd(oro_val)
    ! Purpose: True if > 50% land
    real,intent(in)::oro_val ! [frc] Orography
    oro_is_lnd=nint(oro_val)==1
  end function oro_is_lnd
  
  subroutine aod_nmd_get( & ! [sbr] Process non-mineral dust (nmd) sources
       lat_grd,lat_nbr,lon_grd,lon_nbr,time_nbr, & ! I
       aod_nmd_frc) ! O
    ! Purpose: Find and parameterize non-mineral dust (nmd) sources
    ! aod_nmd_get() is called by io_drv()
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    implicit none
    ! Parameters
    character(len=*),parameter::sbr_nm="aod_nmd_get" ! [sng] Subroutine name
    integer,parameter::rgn_nbr=3 ! [nbr] Dimension size
    ! Commons
    ! Input
    integer lat_nbr ! [nbr] Number of latitudes
    integer lon_nbr ! [nbr] Number of longitudes
    integer time_nbr ! [nbr] Number of times
    real lat_grd(lat_nbr+1) ! [dgr] Interface latitudes
    real lon_grd(lon_nbr+1) ! [dgr] Interface longitudes
    ! Input/Output
    ! Output
    real aod_nmd_frc(lon_nbr,lat_nbr,time_nbr) ! [frc] TOMS-observed aerosol index not due to dust
    ! Local
    integer idx_rgn_angola    ! Named index into rgn array
    integer idx_rgn_brazil    ! Named index into rgn array
    integer idx_rgn_gold_coast ! Named index into rgn array
    integer lat_idx           ! [idx] Counting index
    integer lon_idx           ! [idx] Counting index    
    integer rgn_idx           ! [idx] Counting index
    integer time_idx          ! [idx] Counting index
    real lon_est_out          ! [dgr] Eastern edge of output grid
    real lon_est_rgn          ! [dgr] Eastern edge of input (region) grid
    real rgn_aod_nmd(rgn_nbr,12) ! Seasonal modulation factor
    real rgn_lat_nrt(rgn_nbr) ! [dgr] Northmost latitude of region
    real rgn_lat_sth(rgn_nbr) ! [dgr] Southmost latitude of region
    real rgn_lon_est(rgn_nbr) ! [dgr] Eastmost longitude of region
    real rgn_lon_wst(rgn_nbr) ! [dgr] Westmost longitude of region
    
    real aod_nmd_angola(12) ! Seasonal modulation factor
    real aod_nmd_brazil(12) ! Seasonal modulation factor
    real aod_nmd_gold_coast(12) ! Seasonal modulation factor
    data aod_nmd_angola / &
         0.0,0.0,0.0, &
         0.0,0.0,0.0, &
         0.5,0.5,0.5, &
         0.0,0.0,0.5/
    data aod_nmd_brazil / &
         1.0,1.0,1.0, &
         1.0,1.0,1.0, &
         1.0,1.0,1.0, &
         1.0,1.0,1.0/
    data aod_nmd_gold_coast / &
         0.75,0.75,0.25, &
         0.0,0.0,0.0, &
         0.0,0.0,0.0, &
         0.0,0.0,0.5/
    
    ! Main code
    if (dbg_lvl >= dbg_sbr) write(6,"(a)") "Entering "//sbr_nm
    
    ! Enough memory?
    if (time_nbr /= 12) stop "aod_nmd_get(): time_nbr /= 12"
    
    ! Initialize arrays
    do lon_idx=1,lon_nbr
       do lat_idx=1,lat_nbr
          do time_idx=1,time_nbr
             aod_nmd_frc(lon_idx,lat_idx,time_idx)=0.0
          end do ! end loop over time
       end do ! end loop over lat
    end do ! end loop over lon
    
    ! Named indices must be unique, like a C enum
    idx_rgn_angola=1
    idx_rgn_brazil=2
    idx_rgn_gold_coast=3
    
    ! NB: Make sure last idx equals rgn_nbr
    if (idx_rgn_gold_coast /= rgn_nbr) stop "aod_nmd_get(): ERROR: idx_rgn_gold_coast /= rgn_nbr"
    
    ! Initialize non-mineral dust aerosol arrays
    do rgn_idx=1,rgn_nbr
       rgn_lat_nrt(rgn_idx)=0.0
       rgn_lat_sth(rgn_idx)=0.0
       rgn_lon_est(rgn_idx)=0.0
       rgn_lon_wst(rgn_idx)=0.0
       do time_idx=1,time_nbr
          rgn_aod_nmd(rgn_idx,time_idx)=0.0
       end do ! end loop over mth
    end do ! end loop over rgn
    
    ! Assign values for each non-mineral dust region
    
    ! African SH BMB sources
    rgn_lat_nrt(idx_rgn_angola)=0.0
    rgn_lat_sth(idx_rgn_angola)=-10.0
    rgn_lon_est(idx_rgn_angola)=30.0
    rgn_lon_wst(idx_rgn_angola)=10.0
    do time_idx=1,time_nbr
       rgn_aod_nmd(idx_rgn_angola,time_idx)=aod_nmd_angola(time_idx)
    end do ! end loop over mth
    
    ! Brazilian BMB sources
    rgn_lat_nrt(idx_rgn_brazil)=50.0
    rgn_lat_sth(idx_rgn_brazil)=-50.0
    rgn_lon_est(idx_rgn_brazil)=50.0
    rgn_lon_wst(idx_rgn_brazil)=300.0
    do time_idx=1,time_nbr
       rgn_aod_nmd(idx_rgn_brazil,time_idx)=aod_nmd_brazil(time_idx)
    end do ! end loop over mth
    
    ! African NH BMB sources
    rgn_lat_nrt(idx_rgn_gold_coast)=10.0
    rgn_lat_sth(idx_rgn_gold_coast)=0.0
    rgn_lon_est(idx_rgn_gold_coast)=30.0
    rgn_lon_wst(idx_rgn_gold_coast)=340.0
    do time_idx=1,time_nbr
       rgn_aod_nmd(idx_rgn_gold_coast,time_idx)=aod_nmd_gold_coast(time_idx)
    end do ! end loop over mth
    
    ! Loop over regions and combine all nmd factors into single array
    do rgn_idx=1,rgn_nbr
       ! Ensure longitudes which are compared increase monotonically from West to East
       if (rgn_lon_est(rgn_idx) > rgn_lon_wst(rgn_idx)) then ! region is not wrapped
          lon_est_rgn=rgn_lon_est(rgn_idx) 
       else                   ! region is wrapped
          lon_est_rgn=rgn_lon_est(rgn_idx)+360.0
       endif                  ! endif region is wrapped
       do lon_idx=1,lon_nbr
          if (lon_est_rgn <= 360.0) then ! region is not wrapped
             lon_est_out=lon_grd(lon_idx+1)
          else                ! region is wrapped
             lon_est_out=lon_grd(lon_idx+1)+360.0
          endif               ! endif region is wrapped
          do lat_idx=1,lat_nbr
             if ( &
                  (lat_grd(lat_idx+1) <= rgn_lat_nrt(rgn_idx)).and. &
                  (lat_grd(lat_idx) >= rgn_lat_sth(rgn_idx)).and. &
                  (lon_est_out <= lon_est_rgn).and. &
                  (lon_grd(lon_idx) >= rgn_lon_wst(rgn_idx)).and. &
                  .true.) then
                do time_idx=1,time_nbr
                   aod_nmd_frc(lon_idx,lat_idx,time_idx)= &
                        aod_nmd_frc(lon_idx,lat_idx,time_idx)+ &
                        rgn_aod_nmd(rgn_idx,time_idx)
                end do ! end loop over time
             endif            ! endif
          end do ! end loop over lat
       end do ! end loop over lon
    end do ! end loop over rgn
    
    ! Sanity check
    do lon_idx=1,lon_nbr
       do lat_idx=1,lat_nbr
          do time_idx=1,time_nbr
             if ((aod_nmd_frc(lon_idx,lat_idx,time_idx) < 0.0).or. &
                  (aod_nmd_frc(lon_idx,lat_idx,time_idx) > 1.0)) then
                write (6,"(a,f9.6)") "ERROR aod_nmd_get() reports invalid aod_nmd_frc() = ",aod_nmd_frc(lon_idx,lat_idx,time_idx)
                write (6,"(a)") "Cell edge locations:"
                write (6,"(a)") "(idx)   Lat Sth (idx)   Lat Nrt (idx)   Lon Wst (idx)   Lon Est"
                write (6,"(4(a1,i3,a1,1x,f9.4,1x))") &
                     "(",lat_idx,")",lat_grd(lat_idx),"(",lat_idx+1,")",lat_grd(lat_idx+1), &
                     "(",lon_idx,")",lon_grd(lon_idx),"(",lon_idx+1,")",lon_grd(lon_idx+1)
                stop
             endif
          end do ! end loop over time
       end do ! end loop over lat
    end do ! end loop over lon
    
    if (dbg_lvl >= dbg_sbr) write(6,"(a)") "Exiting "//sbr_nm
    return
  end subroutine aod_nmd_get
  
  subroutine rdb_fct_tpg_GCT01( & ! [sbr] Compute source efficiency factor from topography
       lat_grd,lat_nbr,lon_grd,lon_nbr, & ! I
       hgt_sfc,oro, & ! I
       mbl_bsn_fct) ! O
    ! Purpose: Compute source efficiency factor from topography 
    ! Source: Paul Ginoux ginoux@rondo.gsfc.nasa.gov ${HOME}/dst/ginoux.F:filter()
    ! Modifications by C. Zender 
    ! 1. Generalized for arbitrary horizontal resolution
    ! 2. Removed latitudinal constraints
    ! 3. Removed vertical constraints
    ! 4. Correct handling of poles and date-line
    ! 5. Some provisions for irregular or reduced grids
    ! 6. Added tapering of basin factor for coastal points/islands
    ! rdb_fct_tpg_GCT01() is called by io_drv()
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use map_cst ! [mdl] Constants used in map routines
    implicit none
    ! Parameters
    character(len=*),parameter::sbr_nm="rdb_fct_tpg_GCT01" ! [sng] Subroutine name
    integer,parameter::itr_nbr_max=1 ! [nbr] Maximum number of iterations
    real,parameter::hgt_dlt_bsn_max=8000.0 ! [m] Maximum elevation difference within a basin
    real,parameter::hgt_dlt_bsn_min=1.0 ! [m] Minimum elevation difference within a basin
    ! Ginoux code uses basins that are 15 degrees on a side at 1 degree resolution
    ! Ginoux paper GCT01 states basins are 10x10 degrees
    ! Following code implements the 10x10 basin
    real,parameter::bsn_sz_lon=1111.1e3 ! [m] Zonal size of basin at equator
    real,parameter::bsn_sz_lat=1111.1e3 ! [m] Meridional size of basin at equator
    ! fxm: cst_dst_dmp not yet used (but may be a good idea)
    real,parameter::cst_dst_dmp=100.0e3 ! [m] Distance from coast to damp sources
    ! Commons
    ! Input
    integer,intent(in)::lat_nbr ! [nbr] Dimension size
    integer,intent(in)::lon_nbr ! [nbr] Dimension size
    real,intent(in)::oro(lon_nbr,lat_nbr) ! [frc] Orography
    real,intent(in)::hgt_sfc(lon_nbr,lat_nbr) ! [m] Surface height
    real,intent(in)::lat_grd(lat_nbr+1) ! [dgr] Latitude grid
    real,intent(in)::lon_grd(lon_nbr+1) ! [dgr] Longitude grid
    ! Output
    real,intent(out)::mbl_bsn_fct(lon_nbr,lat_nbr) ! [frc] Mobilization enhancement due to basin characteristics
    ! Local
    real(selected_real_kind(p=12))::pi       ! [frc] 3
    integer bsn_pnt_nbr       ! [nbr] Number of gridpoints in basin
    integer itr_idx           ! [idx] Iteration index
    integer lat_bsn_idx       ! [idx] Counting index
    integer lat_bsn_idx_bnd   ! [idx] Bounded counting index
    integer lon_bsn_idx_bnd   ! [idx] Bounded counting index
    integer bsn_lat_nbr_hlf   ! [nbr] Half number of latitude points in basin
    integer bsn_lon_nbr_hlf   ! [nbr] Half number of longitude points in basin
    integer lat_idx           ! [idx] Counting index
    integer lat_ngh_idx       ! [idx] Latitude neighbor index
    integer lat_ngh_idx_bnd   ! [idx] Bounded latitude neighbor index
    integer lon_bsn_idx       ! [idx] Counting index
    integer lon_idx           ! [idx] Counting index
    integer lon_ngh_idx       ! [idx] Longitude neighbor index
    integer lon_ngh_idx_bnd   ! [idx] Bounded longitude neighbor index
    integer ocn_ngh_nbr       ! [nbr] Number of neighbors which are ocean
    real earth_crc            ! [m] Earth"s circumference
    real hgt_dlt_bsn          ! [m] Elevation range spanned by basin 
    real hgt_dlt_bsn_bnd      ! [m] Bounded elevation range spanned by basin
    real hgt_rng_lcl_nrm      ! [m] Basin peak minus local mean elevation normalized by basin elevation range
    real hgt_max              ! [m] Highest elevation in basin
    real hgt_min              ! [m] Lowest elevation in basin
    real hgt_sfc_bnd(lon_nbr,lat_nbr) ! [m] Surface height bounded
    real hgt_sfc_bsn_avg      ! [m] Mean elevation in basin
    real lat_dlt_dgr          ! [dgr] Cell latitude size
    real lat_dlt_m            ! [m] Cell latitude size
    real lon_dlt_dgr          ! [dgr] Cell longitude size
    real lon_dlt_m            ! [m] Cell longitude size
    real wgt_hgt_rng          ! [frc] Elevation-range-dependent weight
    real wgt_slp              ! [frc] Slope-dependent weight
    real wgt_vgt              ! [frc] Vegetation-dependent weight
    ! Main Code
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Entering "//sbr_nm
    ! Initialize scalars
    pi=4.0*atan(1.0d0) ! [frc] 3
    earth_crc=2.0*pi*earth_rds ! [m] Earth's circumference
    ! Initialize array
    mbl_bsn_fct(:,:)=0.0 ! [frc] Mobilization enhancement due to basin characteristics
    
    ! Assign minimal surface height of 0.0 m
    where (hgt_sfc >= 0.0) ! [m]
       hgt_sfc_bnd=hgt_sfc       ! [m] Surface height bounded
    elsewhere
       hgt_sfc_bnd=0.0        ! [m] Surface height bounded
    end where                 ! end where hgt_sfc >= 0.0
    
    ! Iterate to enlarge surrounding domain
    ! This reinforces weight of large basins compared to small depressions
    do itr_idx=1,itr_nbr_max
       
       ! Compute properties of local "basin" for every latitude slice
       do lat_idx=1,lat_nbr
          
          ! Number of latitude points in basin depends on latitude if grid is irregular
          lat_dlt_dgr=abs((lat_grd(lat_idx+1)-lat_grd(lat_idx))) ! [dgr]
          if (lat_dlt_dgr > 100.0) stop "ERROR: lat_dlt_dgr > 100.0 in rdb_fct_tpg_GCT01()"
          lat_dlt_m=lat_dlt_dgr*earth_crc/360.0 ! [m] Meters per latitude point
          ! Add itr_idx to half-number of points in basin so basins grow with each iteration
          ! Reason for this is unclear, ask Ginoux for rationale
          bsn_lat_nbr_hlf=int(0.5*(itr_idx+bsn_sz_lat/lat_dlt_m)) ! [nbr] Half number of latitude points in basin
          if (2*bsn_lat_nbr_hlf+1 > lat_nbr) stop "ERROR: 2*bsn_lat_nbr_hlf+1 > lat_nbr in rdb_fct_tpg_GCT01()"
          
          do lon_idx=1,lon_nbr
             
             ! Number of longitude points in basin depends on latitude for irregular grids
             ! and may even depend on longitude for very irregular (e.g., reduced) grids
             lon_dlt_dgr=abs((lon_grd(lon_idx+1)-lon_grd(lon_idx))) ! [dgr]
             if (lon_dlt_dgr > 100.0) stop "ERROR: lon_dlt_dgr > 100.0 in rdb_fct_tpg_GCT01()"
             lon_dlt_m=lon_dlt_dgr*earth_crc/360.0 ! [m] Meters per longitude point
             ! Add itr_idx to half-number of points in basin so basins grow with each iteration
             ! Reason for this is unclear, ask Ginoux for rationale
             bsn_lon_nbr_hlf=int(0.5*(itr_idx+bsn_sz_lon/lon_dlt_m)) ! [nbr] Half number of longitude points in basin
             if (2*bsn_lon_nbr_hlf+1 > lon_nbr) stop "ERROR: 2*bsn_lon_nbr_hlf+1 > lon_nbr in rdb_fct_tpg_GCT01()"
             
             ! Initialize description of current basin
             bsn_pnt_nbr=0    ! [nbr] Number of gridpoints in basin
             hgt_min=1.0e36   ! [m] Lowest elevation in basin
             hgt_max=-1.0e36  ! [m] Highest elevation in basin
             hgt_sfc_bsn_avg=0.0 ! [m] Mean elevation in basin
             ocn_ngh_nbr=0    ! [nbr] Number of neighbors which are ocean
             
             ! Count neighboring land points
             do lat_ngh_idx=lat_idx-1,lat_idx+1
                do lon_ngh_idx=lon_idx-1,lon_idx+1
                   lon_ngh_idx_bnd=min(max(1,lon_ngh_idx),lon_nbr)
                   lat_ngh_idx_bnd=min(max(1,lat_ngh_idx),lat_nbr)
                   if (oro_is_ocn(oro(lon_ngh_idx_bnd,lat_ngh_idx_bnd))) ocn_ngh_nbr=ocn_ngh_nbr+1
                end do ! end loop over neighbor lons
             end do ! end loop over neighbor lats
             
             ! Loop over grid points contributing to current basin
             do lat_bsn_idx=lat_idx-bsn_lat_nbr_hlf,lat_idx+bsn_lat_nbr_hlf
                
                ! Bound [lon,lat]_bsn_idx indices to set of valid array indices
                ! Latitude coordinate does not wrap
                lat_bsn_idx_bnd=lat_bsn_idx ! Normal point
                if (lat_bsn_idx < 1.or.lat_bsn_idx > lat_nbr) cycle ! Jump to next iteration of lat_bsn_idx loop
                
                do lon_bsn_idx=lon_idx-bsn_lon_nbr_hlf,lon_idx+bsn_lon_nbr_hlf
                   
                   ! Longitude coordinate wraps
                   lon_bsn_idx_bnd=lon_bsn_idx
                   if (lon_bsn_idx < 1) then 
                      lon_bsn_idx_bnd=lon_nbr-abs(lon_bsn_idx)
                   else if (lon_bsn_idx > lon_nbr) then 
                      lon_bsn_idx_bnd=lon_bsn_idx-lon_nbr
                   endif      ! endif
                   
                   bsn_pnt_nbr=bsn_pnt_nbr+1
                   hgt_sfc_bsn_avg=hgt_sfc_bsn_avg+hgt_sfc_bnd(lon_bsn_idx_bnd,lat_bsn_idx_bnd)
                   ! Elevation extrema in current basin
                   hgt_max=max(hgt_max,hgt_sfc_bnd(lon_bsn_idx_bnd,lat_bsn_idx_bnd))
                   hgt_min=min(hgt_min,hgt_sfc_bnd(lon_bsn_idx_bnd,lat_bsn_idx_bnd))
                end do ! end loop over lon_bsn_idx
             end do ! end loop over lat_bsn_idx
             
             ! Sanity checks
             if (bsn_pnt_nbr > (2*bsn_lat_nbr_hlf+1)*(2*bsn_lon_nbr_hlf+1)) stop "ERROR: bsn_pnt_nbr too large in rdb_fct_tpg_GCT01"
             if (bsn_pnt_nbr > 0) then
                hgt_sfc_bsn_avg=hgt_sfc_bsn_avg/bsn_pnt_nbr
             else
                hgt_sfc_bsn_avg=0.0
             endif            ! endif bsn_pnt_nbr
             hgt_dlt_bsn=hgt_max-hgt_min ! [m] Elevation range spanned by basin 
             if (hgt_dlt_bsn < 0.0) stop "ERROR: hgt_dlt_bsn < 0.0 in rdb_fct_tpg_GCT01()"
             ! Approximately 30 basins (at T42) have hgt_dlt_bsn = 0.0 m
             ! Enforce minimum hgt_dlt_bsn = 1.0 m in these regions to avoid divide-by-zero below in wgt_bsn
             hgt_dlt_bsn_bnd=max(hgt_dlt_bsn_min,hgt_dlt_bsn) ! [m] Bounded elevation range spanned by basin
             
             ! Potential sources are points above MSL with some land
             if (oro(lon_idx,lat_idx) > 0.5) then
                
                ! wgt_slp_fct2=hgt_dlt_bsn_bnd/1000.0 ! [frc] Factor used in slope weight
                ! wgt_slp=(8.0-0.25*wgt_slp_fct2)/8.0 ! [frc] Slope-dependent weight
                wgt_slp=1.0   ! [frc] Slope-dependent weight
                
                ! Basin characteristics
                ! Find elevation difference between peak of basin and local mean elevation and normalize this by 
                ! (bounded) height difference anywhere in basin
                ! Combined with following step, this ensures deeper portions of basin have more mobilization potential
                hgt_rng_lcl_nrm=(hgt_max-hgt_sfc_bnd(lon_idx,lat_idx))/hgt_dlt_bsn_bnd ! [m] Basin peak minus local mean elevation normalized by basin elevation range
                
                ! 0.0 < wgt_hgt_rng <= 1.0:
                ! wgt_hgt_rng is 0.0 when local point is highest point in basin
                ! wgt_hgt_rng is 1.0 when local point is lowest point in basin
                ! Fifth power dependence comes from Ginoux's tuning
                wgt_hgt_rng=hgt_rng_lcl_nrm**5 ! [frc] Elevation-range-dependent weight
                
                !!     Vegetation is not a topographic variable so it is counted in lnd_frc_mbl
                !              if (sfc_typ_idx==0.or. & ! Exclude ocean
                !                  sfc_typ_idx==1.or. & ! Exclude ice and sea ice 
                !!                      sfc_typ_idx==26.or. & ! Exclude warm crops (mostly Europe)
                !                  sfc_typ_idx==27.or. & ! Exclude wetlands
                !                  sfc_typ_idx==28.or. & ! Exclude wetlands
                !                  .false.) then 
                !              wgt_vgt=0.0   ! [frc] Vegetation-dependent weight
                !           else
                !              wgt_vgt=1.0   ! [frc] Vegetation-dependent weight
                !           endif
                wgt_vgt=1.0   ! [frc] Vegetation-dependent weight
                
                ! Increment mbl_bsn_fct (rather than straight assign) to take advantage of iteration capability
                mbl_bsn_fct(lon_idx,lat_idx)= & ! [frc] Mobilization enhancement due to basin characteristics
                     mbl_bsn_fct(lon_idx,lat_idx)+ &
                     wgt_slp*wgt_hgt_rng*wgt_vgt/itr_nbr_max
                
             endif            ! end if point is potential source
             
          end do ! end loop over lon
       end do ! end loop over lat
       
    end do ! end loop over itr
    
    ! Sanity check
    do lat_idx=1,lat_nbr
       do lon_idx=1,lon_nbr
          if (mbl_bsn_fct(lon_idx,lat_idx) < 0.0.or.mbl_bsn_fct(lon_idx,lat_idx) > 1.0) then
             write (6,"(a)") "rdb_fct_tpg_GCT01(): ERROR "
             write (6,"(9(a,es9.2,a))") &
                  "hgt_sfc =",hgt_sfc(lon_idx,lat_idx)," m, ", &
                  "hgt_sfc_bnd =",hgt_sfc_bnd(lon_idx,lat_idx)," m, ", &
                  "hgt_sfc_bsn_avg =",hgt_sfc_bsn_avg," m, ", &
                  "hgt_min =",hgt_min," m, ", &
                  "hgt_max =",hgt_max," m, ", &
                  "hgt_dlt_bsn =",hgt_dlt_bsn," m, ", &
                  "hgt_dlt_bsn_bnd =",hgt_dlt_bsn_bnd," m, ", &
                  "hgt_rng_lcl_nrm =",hgt_rng_lcl_nrm," m, ", &
                  "mbl_bsn_fct =",mbl_bsn_fct(lon_idx,lat_idx)," frc, "
             write (6,"(a)") "Output cell edge locations:"
             write (6,"(a)") "(idx)   Lat Sth (idx)   Lat Nrt (idx)   Lon Wst (idx)   Lon Est"
             write (6,"(4(a1,i3,a1,1x,f9.4,1x))") &
                  "(",lat_idx,")",lat_grd(lat_idx),"(",lat_idx+1,")",lat_grd(lat_idx+1), &
                  "(",lon_idx,")",lon_grd(lon_idx),"(",lon_idx+1,")",lon_grd(lon_idx+1)
          endif               ! endif err
       end do ! end loop over lon
    end do ! end loop over lat

    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Exiting "//sbr_nm
    return 
  end subroutine rdb_fct_tpg_GCT01
  
  subroutine src_id( & ! [sbr] Find and parameterize dust sources
       odxc_in,mss_flg,mss_val, & ! I
       lat_grd,lat_nbr,lon_grd,lon_nbr, & ! I
       lnd_msk, & ! I
       odxc_out,odxc_tll,src_frq,src_str, & ! I/Output
       src_flg) ! O
    ! Purpose: Find and parameterize dust sources
    ! src_id() is called by toms_get()
    use bds_ctl,only:src_thr ! [mdl] Control variables drc_in,drc_out,hgt_dlt_msl
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    implicit none
    ! Parameters
    character(len=*),parameter::sbr_nm="src_id" ! [sng] Subroutine name
    ! Commons
    ! Input
    integer lat_nbr ! [nbr] Number of latitudes
    integer lon_nbr ! [nbr] Number of longitudes
    integer lnd_msk(lon_nbr,lat_nbr) ! Land-water mask
    integer odxc_tll(lon_nbr,lat_nbr) ! Tally of valid optical depth months
    logical mss_flg           ! [flg] Check for mss_val on input data
    real lat_grd(lat_nbr+1) ! [dgr] Interface latitudes
    real lon_grd(lon_nbr+1) ! [dgr] Interface longitudes
    real mss_val              ! Missing value
    real odxc_in(lon_nbr,lat_nbr) ! [frc] Extinction optical depth
    ! Input/Output
    real odxc_out(lon_nbr,lat_nbr) ! [frc] Extinction optical depth
    real src_frq(lon_nbr,lat_nbr) ! Source frequency
    real src_str(lon_nbr,lat_nbr) ! [frc] Source strength
    ! Output
    integer src_flg(lon_nbr,lat_nbr) ! [flg] Source flag
    ! Local
    integer cnt               ! Counter
    integer cnt_thr           ! Counter source threshold
    integer lat_idx           ! [idx] Counting index
    integer lon_idx           ! [idx] Counting index    
    real odxc_pad(0:lon_nbr+1,0:lat_nbr+1) ! [frc] Padded input optical depth
    real thr_top(0:lon_nbr+1,0:lat_nbr+1) ! [frc] Padded input optical depth times threshold
    ! Main code
    if (dbg_lvl >= dbg_sbr) write(6,"(a)") "Entering "//sbr_nm
    
    ! Initialize scalars
    cnt_thr=8                 ! Pairwise source criteria must be met for at least cnt_thr neighbors
    
    ! Accumulate optical depth
    do lon_idx=1,lon_nbr
       do lat_idx=1,lat_nbr
          if (odxc_in(lon_idx,lat_idx) /= mss_val) then 
             odxc_out(lon_idx,lat_idx)=odxc_out(lon_idx,lat_idx)+odxc_in(lon_idx,lat_idx)
             odxc_tll(lon_idx,lat_idx)=odxc_tll(lon_idx,lat_idx)+1
          endif               ! mss_val
       end do ! end loop over lat
    end do ! end loop over lon
    
    ! Initialize arrays
    do lon_idx=1,lon_nbr
       do lat_idx=1,lat_nbr
          src_flg(lon_idx,lat_idx)=0
          odxc_pad(lon_idx,lat_idx)=odxc_in(lon_idx,lat_idx)
       end do ! end loop over lat
    end do ! end loop over lon
    ! Pad latitude
    do lon_idx=1,lon_nbr
       odxc_pad(lon_idx,0)=odxc_pad(lon_idx,1)
       odxc_pad(lon_idx,lat_nbr+1)=odxc_pad(lon_idx,lat_nbr)
    end do ! end loop over lon
    ! Pad longitude
    do lat_idx=0,lat_nbr+1
       odxc_pad(0,lat_idx)=odxc_pad(lon_nbr,lat_idx)
       odxc_pad(lon_nbr+1,lat_idx)=odxc_pad(1,lat_idx)
    end do ! end loop over lat
    ! Construct threshold optical depths
    do lon_idx=0,lon_nbr+1
       do lat_idx=0,lat_nbr+1
          ! NB: In TOMS data, missing values are poleward of about 60 degrees
          if (mss_flg.and.odxc_pad(lon_idx,lat_idx)==mss_val) then 
             thr_top(lon_idx,lat_idx)=mss_val
          else                ! mss_val
             thr_top(lon_idx,lat_idx)=src_thr*odxc_pad(lon_idx,lat_idx)
          endif               ! endif
       end do ! end loop over lat
    end do ! end loop over lon
    
    ! Select maxima (i.e., sources)
    do lon_idx=1,lon_nbr
       do lat_idx=1,lat_nbr
          ! Initialize flags and counters for each candidate gridpoint
          cnt=0
          ! Currently, a source must ...
          if ( &
               ! ... be a valid datapoint
               (.not.mss_flg.or.(mss_flg.and.odxc_pad(lon_idx,lat_idx) /= mss_val)).and. &
               ! ... be a land gridpoint which exceeds the TOMS threshold (has non-zero optical depth)
               (lnd_msk(lon_idx,lat_idx)==1.and.odxc_pad(lon_idx,lat_idx) > 0.0) &
               ) then
             ! ... and which exceeds a threshold (e.g., 90%) of the strength of 
             if (odxc_pad(lon_idx,lat_idx) >= thr_top(lon_idx,lat_idx+1)) cnt=cnt+1 ! N
             if (odxc_pad(lon_idx,lat_idx) >= thr_top(lon_idx+1,lat_idx+1)) cnt=cnt+1 ! NE
             if (odxc_pad(lon_idx,lat_idx) >= thr_top(lon_idx+1,lat_idx)) cnt=cnt+1 ! E
             if (odxc_pad(lon_idx,lat_idx) >= thr_top(lon_idx+1,lat_idx-1)) cnt=cnt+1 ! SE
             if (odxc_pad(lon_idx,lat_idx) >= thr_top(lon_idx,lat_idx-1)) cnt=cnt+1 ! S
             if (odxc_pad(lon_idx,lat_idx) >= thr_top(lon_idx-1,lat_idx-1)) cnt=cnt+1 ! SW
             if (odxc_pad(lon_idx,lat_idx) >= thr_top(lon_idx-1,lat_idx)) cnt=cnt+1 ! W
             if (odxc_pad(lon_idx,lat_idx) >= thr_top(lon_idx-1,lat_idx+1)) cnt=cnt+1 ! NW
             ! ... a certain number (e.g., seven) of its eight neighbors.
             if (cnt >= cnt_thr) then 
                src_flg(lon_idx,lat_idx)=1
                src_frq(lon_idx,lat_idx)=src_frq(lon_idx,lat_idx)+1.0
             endif               ! endif src
          endif                  ! endif land
       end do ! end loop over lat
    end do ! end loop over lon
    
    ! Correct absorbing aerosol index for non-dust aerosols, if necessary
    
    ! Parameterize strengths
    do lon_idx=1,lon_nbr
       do lat_idx=1,lat_nbr
          if (src_flg(lon_idx,lat_idx) > 0) then
             src_str(lon_idx,lat_idx)=src_str(lon_idx,lat_idx)+odxc_pad(lon_idx,lat_idx)
          endif               ! endif source
       end do ! end loop over lat
    end do ! end loop over lon
    
    ! Sanity check
    do lon_idx=1,lon_nbr
       do lat_idx=1,lat_nbr
          if (src_flg(lon_idx,lat_idx) > 0.and.odxc_pad(lon_idx,lat_idx) <= 0.0) then
             write(6,"(a,i3,a,i3,a,a,f6.4,a,i2)") &
                  "src_id(): (lon_idx,lat_idx) = (", &
                  lon_idx,",",lat_idx,"), ", &
                  " odxc_pad = ",odxc_pad(lon_idx,lat_idx), &
                  " src_flg = ",src_flg(lon_idx,lat_idx)
          endif               ! endif source
       end do ! end loop over lat
    end do ! end loop over lon
    
    if (dbg_lvl==dbg_old) then
       do lon_idx=1,lon_nbr
          do lat_idx=1,lat_nbr
             write(6,"(a,i3,a,i3,a,f6.4)") &
                  "src_id(): src_str(", &
                  lon_idx,",",lat_idx,") = ", &
                  src_str(lon_idx,lat_idx)
          end do ! end loop over lat
       end do ! end loop over lon
    endif                     ! endif dbg
    
    if (dbg_lvl >= dbg_sbr) write(6,"(a)") "Exiting "//sbr_nm
    return
  end subroutine src_id

  subroutine cst_tpr( & ! [sbr] Taper field at coast
       lat,lat_nbr,lon,lon_nbr, & ! I
       lat_min,lat_max,lon_min,lon_max, & ! I 
       oro, & ! I
       arr_2D) ! I/O
    ! Purpose: Set arr_2D to zero (or some fraction) if one or more costal points
    ! cst_tpr() is called by rdb_fct_tpg_GCT01()
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use map_cst ! [mdl] Constants used in map routines
    implicit none
    ! fxm: cst_dst_dmp not yet used (but may be a good idea)
    ! Parameters
    character(len=*),parameter::sbr_nm="cst_tpr" ! [sng] Subroutine name
    real,parameter::cst_dst_dmp=100.0e3 ! [m] Distance from coast to damp sources
    ! Input
    integer,intent(in)::lat_nbr   ! [nbr] Number of latitudes
    integer,intent(in)::lon_nbr   ! [nbr] Number of longitudes
    real,intent(in)::lat(lat_nbr) ! [dgr] Latitude
    real,intent(in)::lon(lon_nbr) ! [dgr] Longitude
    real,intent(in)::lat_min      ! [dgr] Min latitude  to apply filter
    real,intent(in)::lat_max      ! [dgr] Max latitude  to apply filter
    real,intent(in)::lon_min      ! [dgr] Min longitude to apply filter
    real,intent(in)::lon_max      ! [dgr] Max longitude to apply filter
    real,intent(in)::oro(lon_nbr,lat_nbr) ! [frc] Orography
    ! Output
    real,intent(inout)::arr_2D(lon_nbr,lat_nbr) ! [frc] Field to filter
    ! Local
    integer lat_idx           ! [idx] Counting index
    integer lat_ngh_idx       ! [idx] Latitude neighbor index
    integer lat_ngh_idx_bnd   ! [idx] Bounded latitude neighbor index
    integer lon_idx           ! [idx] Counting index
    integer lon_ngh_idx       ! [idx] Longitude neighbor index
    integer lon_ngh_idx_bnd   ! [idx] Bounded longitude neighbor index
    integer ocn_ngh_nbr       ! [nbr] Number of neighbors which are ocean
    integer zro_cnt           ! [nbr] Number of points zeroed
    real lon_crr              ! [dgr] Current longitude
    real lon_crr_p360         ! [dgr] Current longitude shifted 360 degrees East
    real lon_crr_m360         ! [dgr] Current longitude shifted 360 degrees West
    real wgt_ocn_ngh          ! [frc] Coastline-dependent weight

    ! Main Code
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Entering "//sbr_nm
    
    ! Loop through all latitudes and longitudes within limits
    zro_cnt=0
    do lat_idx=1,lat_nbr
       if ((lat(lat_idx) >= lat_min).and.(lat(lat_idx) <= lat_max)) then
          do lon_idx=1,lon_nbr
             lon_crr=lon(lon_idx) ! [dgr] Current longitude
             lon_crr_p360=lon_crr+360.0 ! [dgr] Current longitude shifted 360 degrees East
             lon_crr_m360=lon_crr-360.0 ! [dgr] Current longitude shifted 360 degrees West
             if ( &
                (lon_crr >= lon_min).and.(lon_crr <= lon_max).or. &
                (lon_crr_p360 >= lon_min).and.(lon_crr_p360 <= lon_max).or. &
                (lon_crr_m360 >= lon_min).and.(lon_crr_m360 <= lon_max).or. &
                .false.) then
                
                ! Only land points are tapered
                if (oro_is_lnd(oro(lon_idx,lat_idx))) then
                   
                   ! Initialize current basin counter
                   ocn_ngh_nbr=0 ! [nbr] Number of neighbors which are ocean
                   
                   ! Count neighboring land points
                   do lat_ngh_idx=lat_idx-1,lat_idx+1
                      do lon_ngh_idx=lon_idx-1,lon_idx+1
                         lon_ngh_idx_bnd=min(max(1,lon_ngh_idx),lon_nbr)
                         lat_ngh_idx_bnd=min(max(1,lat_ngh_idx),lat_nbr)
                         if (oro_is_ocn(oro(lon_ngh_idx_bnd,lat_ngh_idx_bnd))) ocn_ngh_nbr=ocn_ngh_nbr+1
                      end do ! end loop over neighbor lons
                   end do ! end loop over neighbor lats
                   
                   ! Coastal points are usually lower than surrounding inland points
                   ! However, they are not really basins, just victims of geography
                   ! Taper mobilization efficiency at coastlines by dividing by two for every neighboring non-land point
                   ! Gridpoint islands should not be basins either
                   ! 20020305: fxm: Test not tapering points with two or fewer ocean neighbors
                   if (ocn_ngh_nbr == 0) then ! No coasts allowed
                      !                if (ocn_ngh_nbr <= 2) then ! No tapering, but some coasts allowed
                      wgt_ocn_ngh=1.0 ! [frc] Coastline-dependent weight
                      ! wgt_ocn_ngh=2.0**(-ocn_ngh_nbr) ! [frc] Taper coastlines
                   else 
                      wgt_ocn_ngh=0.0 ! [frc] Coastline-dependent weight
                      zro_cnt = zro_cnt + 1
                   endif ! endif coast
                   
                   ! Sanity check
                   if (wgt_ocn_ngh < 1.0 .and. ocn_ngh_nbr == 0) then
                      write (6,"(a,f9.6,a,2i3)") &
                           "ERROR"//sbr_nm//"reports invalid coast weight = ", &
                           wgt_ocn_ngh," at (lon_idx,lat_idx) = ",lon_idx,lat_idx
                      stop
                   endif        ! endif insane
                   arr_2D(lon_idx,lat_idx)=arr_2D(lon_idx,lat_idx)*wgt_ocn_ngh
                endif ! oro_is_lnd 
             endif ! endif
          end do ! end loop over lon
       endif ! endif
    end do ! end loop over lat
    
    if (dbg_lvl >= dbg_sbr) write (6,"(a,i4)") "INFO: zro_cnt = ", zro_cnt
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Exiting "//sbr_nm
    return 
  end subroutine cst_tpr
  
  subroutine rgn_msk( & ! [sbr] Zero field based on location
       lat,lat_nbr,lon,lon_nbr, & ! I
       lat_min,lat_max,lon_min,lon_max, & ! I 
       arr_2D) ! I/O
    ! Purpose: Mask field based on location
    ! Currently this sets field true (unity) inside rectangular regions in bds_rgn_msk.txt
    ! NB: Calling routine should initialize variable to false everywhere---
    ! This routine does not set exterior regions to false (zero)!
    ! rgn_msk() is called by bds_prc()
    ! 20050726: Based on rgn_xcl
    ! fxm: Routine could/should be made more general by taking optional
    ! arguments values:
    ! val_ntr: I [frc] Value for points interior to region, default is 0.0
    ! flg_xtr: I [flg] Set points exterior to region to val_xtr, default is False
    ! val_xtr: I [frc] Value for points exterior to region, default is 1.0
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    implicit none
    ! Input
    integer,intent(in)::lat_nbr   ! [nbr] Number of latitudes
    integer,intent(in)::lon_nbr   ! [nbr] Number of longitudes
    real,intent(in)::lat(lat_nbr) ! [dgr] Latitude
    real,intent(in)::lon(lon_nbr) ! [dgr] Longitude
    real,intent(in)::lat_min      ! [dgr] Min latitude  to apply filter
    real,intent(in)::lat_max      ! [dgr] Max latitude  to apply filter
    real,intent(in)::lon_min      ! [dgr] Min longitude to apply filter
    real,intent(in)::lon_max      ! [dgr] Max longitude to apply filter
    ! Output
    integer,intent(inout)::arr_2D(lon_nbr,lat_nbr) ! [frc] Field to be filtered
    ! Local
    integer lat_idx           ! [idx] Counting index
    integer lon_idx           ! [idx] Counting index
    integer msk_cnt           ! [nbr] Number of points zeroed
    real lon_crr              ! [dgr] Current longitude
    real lon_crr_p360         ! [dgr] Current longitude shifted 360 degrees East
    real lon_crr_m360         ! [dgr] Current longitude shifted 360 degrees West
    real wgt_lcn              ! [frc] Location-dependent weight
    ! Main Code
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Entering rgn_msk()"
    
    ! Loop through all latitudes and longitudes
    msk_cnt=0
    do lat_idx=1,lat_nbr
       ! Find all points interior to specified rectangle
       if ((lat(lat_idx) >= lat_min).and.(lat(lat_idx) <= lat_max)) then
          do lon_idx=1,lon_nbr
             lon_crr=lon(lon_idx) ! [dgr] Current longitude
             lon_crr_p360=lon_crr+360.0 ! [dgr] Current longitude shifted 360 degrees East
             lon_crr_m360=lon_crr-360.0 ! [dgr] Current longitude shifted 360 degrees West
             if ( &
                  ((lon_crr >= lon_min).and.(lon_crr <= lon_max)).or. &
                  ((lon_crr_p360 >= lon_min).and.(lon_crr_p360 <= lon_max)).or. &
                  ((lon_crr_m360 >= lon_min).and.(lon_crr_m360 <= lon_max)).or. &
                  .false.) then
                ! Find points interior to specified rectangle are set to val_ntr
                wgt_lcn=0.0 ! [frc] Location-dependent weight
                msk_cnt=msk_cnt+1
                arr_2D(lon_idx,lat_idx)=1
             endif ! endif
          end do ! end loop over lon
       endif ! endif
    end do ! end loop over lat
    if (dbg_lvl >= dbg_sbr) write (6,"(a,i4)") "INFO: msk_cnt = ",msk_cnt
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Exiting rgn_msk()"
    return 
  end subroutine rgn_msk
  
  subroutine rgn_xcl( & ! [sbr] Zero field based on location
       lat,lat_nbr,lon,lon_nbr, & ! I
       lat_min,lat_max,lon_min,lon_max, & ! I 
       arr_2D) ! I/O
    ! Purpose: Zero field based on location
    ! Currently this zeros field based on rectangular regions in bds_rgn_xcl.txt
    ! rgn_xcl() is called by bds_prc()
    ! 20020522: Original use of routine is to remove excessive Canadian/Siberian emissions
    ! This is required for reasonable fluxes using CLM vegetation
    ! Atmospheric burdens are not much affected since washout is fast here
    ! fxm: Routine could/should be made more general by taking optional
    ! arguments values:
    ! val_ntr: I [frc] Value for points interior to region, default is 0.0
    ! flg_xtr: I [flg] Set points exterior to region to val_xtr, default is False
    ! val_xtr: I [frc] Value for points exterior to region, default is 1.0
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    implicit none
    ! Input
    integer,intent(in)::lat_nbr   ! [nbr] Number of latitudes
    integer,intent(in)::lon_nbr   ! [nbr] Number of longitudes
    real,intent(in)::lat(lat_nbr) ! [dgr] Latitude
    real,intent(in)::lon(lon_nbr) ! [dgr] Longitude
    real,intent(in)::lat_min      ! [dgr] Min latitude  to apply filter
    real,intent(in)::lat_max      ! [dgr] Max latitude  to apply filter
    real,intent(in)::lon_min      ! [dgr] Min longitude to apply filter
    real,intent(in)::lon_max      ! [dgr] Max longitude to apply filter
    ! Output
    real,intent(inout)::arr_2D(lon_nbr,lat_nbr) ! [frc] Field to be filtered
    ! Local
    integer lat_idx           ! [idx] Counting index
    integer lon_idx           ! [idx] Counting index
    integer zro_cnt           ! [nbr] Number of points zeroed
    real lon_crr              ! [dgr] Current longitude
    real lon_crr_p360         ! [dgr] Current longitude shifted 360 degrees East
    real lon_crr_m360         ! [dgr] Current longitude shifted 360 degrees West
    real wgt_lcn              ! [frc] Location-dependent weight
    ! Main Code
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Entering rgn_xcl()"

    ! Loop through all latitudes and longitudes
    zro_cnt=0
    do lat_idx=1,lat_nbr
       ! Find all points interior to specified rectangle
       if ((lat(lat_idx) >= lat_min).and.(lat(lat_idx) <= lat_max)) then
          do lon_idx=1,lon_nbr
             lon_crr=lon(lon_idx) ! [dgr] Current longitude
             lon_crr_p360=lon_crr+360.0 ! [dgr] Current longitude shifted 360 degrees East
             lon_crr_m360=lon_crr-360.0 ! [dgr] Current longitude shifted 360 degrees West
             if ( &
                (lon_crr >= lon_min).and.(lon_crr <= lon_max).or. &
                (lon_crr_p360 >= lon_min).and.(lon_crr_p360 <= lon_max).or. &
                (lon_crr_m360 >= lon_min).and.(lon_crr_m360 <= lon_max).or. &
                .false.) then
                ! Find points interior to specified rectangle are set to val_ntr
                
                wgt_lcn=0.0 ! [frc] Location-dependent weight
                zro_cnt=zro_cnt+1
                
                arr_2D(lon_idx,lat_idx)=arr_2D(lon_idx,lat_idx)*wgt_lcn
             endif ! endif
          end do ! end loop over lon
       endif ! endif
    end do ! end loop over lat
    if (dbg_lvl >= dbg_sbr) write (6,"(a,i4)") "INFO: zro_cnt = ",zro_cnt
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Exiting rgn_xcl()"
 return 
end subroutine rgn_xcl

end module src_prc ! [mdl] Source processing, identification
