! $Id$ -*-f90-*-

! Purpose: Routines to process maps and regridding

! Copyright (C) 2001 Charlie Zender
! Portions are copyright by their respective contributors

! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 3
! of the License, or (at your option) any later version.
! See http://www.gnu.org/copyleft/gpl.html for details.

! The author of this software, Charlie Zender, would like to receive
! your suggestions, improvements, bug-reports, and patches. 
! Charlie Zender, zender@uci.edu
! Department of Earth System Science
! University of California at Irvine
! Irvine, CA 92697-3100

! Usage:
!use map_grd ! [mdl] Map grids and regridding

module map_grd ! [mdl] Map grids and regridding
  implicit none
  public ! [stt] Symbols are public unless individually qualified as private
!  public::map_typ_sng_get ! [sbr] Print map type
  
contains

  subroutine map_rbn(dat_in,mss_flg,mss_val, & ! I
       lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
       lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
       ovr_nbr_max,ovr_lat_idx,ovr_lon_idx,ovr_nbr,ovr_wgt, & ! I
       dat_out) ! O
    ! Purpose: Rebin continuous valued field from input grid to output grid
    ! Missing values, if specified, are handled correctly
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    implicit none
    ! Parameters
    ! Commons
    ! Input
    integer,intent(in)::lat_in_nbr ! I [nbr] Number of latitudes
    integer,intent(in)::lat_out_nbr ! I [nbr] Number of latitudes
    integer,intent(in)::lon_in_nbr ! I [nbr] Number of longitudes
    integer,intent(in)::lon_out_nbr ! I [nbr] Number of longitudes
    integer,intent(in)::ovr_nbr_max ! I [nbr] Maximum number of input cells which overlap any output cell
    integer,intent(in)::ovr_lat_idx(lon_out_nbr,lat_out_nbr,ovr_nbr_max) ! I [idx] Map into input grid of latitude indices of overlap cells
    integer,intent(in)::ovr_lon_idx(lon_out_nbr,lat_out_nbr,ovr_nbr_max) ! I [idx] Map into input grid of longitude indices of overlap cells
    integer,intent(in)::ovr_nbr(lon_out_nbr,lat_out_nbr) !  I [nbr] Number of input gridcells which overlap each output gridcell
    logical,intent(in)::mss_flg ! I [flg] Check for mss_val on input data
    real,intent(in)::dat_in(lon_in_nbr,lat_in_nbr) ! I data
    real,intent(in)::lat_in_grd(lat_in_nbr+1) ! I [dgr] Interface latitudes
    real,intent(in)::lat_out_grd(lat_out_nbr+1) ! I [dgr] Interface latitudes
    real,intent(in)::lon_in_grd(lon_in_nbr+1) ! I [dgr] Interface longitudes
    real,intent(in)::lon_out_grd(lon_out_nbr+1) ! I [dgr] Interface longitudes
    real,intent(in)::mss_val  ! I [frc] Missing value
    real,intent(in)::ovr_wgt(lon_out_nbr,lat_out_nbr,ovr_nbr_max) ! I [frc] Weight of overlapping input gridcells onto each output gridcell
    ! Output
    real,intent(out)::dat_out(lon_out_nbr,lat_out_nbr) ! O Mean of input field 
    ! Local
    integer lat_out_idx       ! [idx] Counting index
    integer lon_out_idx       ! [idx] Counting index
    integer lat_in_idx        ! [idx] Counting index
    integer lon_in_idx        ! [idx] Counting index
    integer mss_cnt           ! [nbr] Number of missing values in current output gridcell
    integer ovr_idx           ! [idx] Counting index
    integer ovr_nbr_crr       ! [nbr] Current number of overlapping gridcells
    
    ! Main Code
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Entering map_rbn()"
    
    ! Initialize arrays
    do lon_out_idx=1,lon_out_nbr
       do lat_out_idx=1,lat_out_nbr
          dat_out(lon_out_idx,lat_out_idx)=0.0
       end do ! end loop over lat
    end do ! end loop over lon
    
    !$omp parallel default(none) &
    !$omp private(lat_out_idx,lon_out_idx) &
    !$omp private(lat_in_idx,lon_in_idx) &
    !$omp private(ovr_nbr_crr,mss_cnt,ovr_idx) &
    !$omp shared(lat_out_nbr,lon_out_nbr) &
    !$omp shared(ovr_nbr,ovr_lon_idx,ovr_lat_idx,ovr_wgt) &
    !$omp shared(dat_in,dat_out,mss_val,mss_flg)
    
    ! Rebin input data to output grid
    if (mss_flg) then
       !$omp do
       do lat_out_idx=1,lat_out_nbr
          do lon_out_idx=1,lon_out_nbr
             ovr_nbr_crr=ovr_nbr(lon_out_idx,lat_out_idx)
             mss_cnt=0
             do ovr_idx=1,ovr_nbr_crr
                lon_in_idx=ovr_lon_idx(lon_out_idx,lat_out_idx,ovr_idx) ! [idx] lon index (input grid) of overlap cell
                lat_in_idx=ovr_lat_idx(lon_out_idx,lat_out_idx,ovr_idx) ! [idx] lat index (input grid) of overlap cell
                if (dat_in(lon_in_idx,lat_in_idx) /= mss_val) then
                   dat_out(lon_out_idx,lat_out_idx)= &
                        dat_out(lon_out_idx,lat_out_idx)+ &
                        ovr_wgt(lon_out_idx,lat_out_idx,ovr_idx)*dat_in(lon_in_idx,lat_in_idx)
                else          ! mss_val
                   mss_cnt=mss_cnt+1
                endif         ! mss_val
             end do ! end loop over ovr
             if (mss_cnt==ovr_nbr_crr) then
                dat_out(lon_out_idx,lat_out_idx)=mss_val
             endif            ! endif
          end do ! end loop over lon
       end do ! end loop over lat
       !$omp end do
    endif                     ! not mss_flg
    if (.not.mss_flg) then
       !$omp do
       do lat_out_idx=1,lat_out_nbr
          do lon_out_idx=1,lon_out_nbr
             ovr_nbr_crr=ovr_nbr(lon_out_idx,lat_out_idx)
             do ovr_idx=1,ovr_nbr_crr
                lon_in_idx=ovr_lon_idx(lon_out_idx,lat_out_idx,ovr_idx) ! [idx] lon index (input grid) of overlap cell
                lat_in_idx=ovr_lat_idx(lon_out_idx,lat_out_idx,ovr_idx) ! [idx] lat index (input grid) of overlap cell
                dat_out(lon_out_idx,lat_out_idx)= &
                     dat_out(lon_out_idx,lat_out_idx)+ &
                     ovr_wgt(lon_out_idx,lat_out_idx,ovr_idx)*dat_in(lon_in_idx,lat_in_idx)
             end do           ! end loop over ovr
          end do              ! end loop over lon
       end do                 ! end loop over lat
       !$omp end do
    endif ! mss_flg
    !$omp end parallel
    
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Exiting map_rbn()"
    return 
  end subroutine map_rbn                       ! end map_rbn()
  
  subroutine map_rbn_var(dat_in,mss_flg,mss_val, & ! I
       lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
       lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
       ovr_nbr_max,ovr_lat_idx,ovr_lon_idx,ovr_nbr,ovr_wgt, & ! I
       dat_out, & ! I
       dat_var_out) ! O
    ! Purpose: Find variance of continuous valued field on coarse output grid
    ! Given an input grid, a field, and the mean of the field on an output grid,
    ! return the variance of the field on the output grid
    ! Missing values, if specified, are handled correctly
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    implicit none
    ! Parameters
    ! Commons
    ! Input
    integer,intent(in)::lat_in_nbr ! I [nbr] Number of latitudes
    integer,intent(in)::lat_out_nbr ! I [nbr] Number of latitudes
    integer,intent(in)::lon_in_nbr ! I [nbr] Number of longitudes
    integer,intent(in)::lon_out_nbr ! I [nbr] Number of longitudes
    integer,intent(in)::ovr_nbr_max ! I [nbr] Maximum number of input cells which overlap any output cell
    integer,intent(in)::ovr_lat_idx(lon_out_nbr,lat_out_nbr,ovr_nbr_max) ! I [idx] Map into input grid of latitude indices of overlap cells
    integer,intent(in)::ovr_lon_idx(lon_out_nbr,lat_out_nbr,ovr_nbr_max) ! I [idx] Map into input grid of longitude indices of overlap cells
    integer,intent(in)::ovr_nbr(lon_out_nbr,lat_out_nbr) !  I [nbr] Number of input gridcells which overlap each output gridcell
    logical,intent(in)::mss_flg ! I [flg] Check for mss_val on input data
    real,intent(in)::dat_in(lon_in_nbr,lat_in_nbr) ! I [frc] Data
    real,intent(in)::dat_out(lon_out_nbr,lat_out_nbr) ! [frc] Mean of input field 
    real,intent(in)::lat_in_grd(lat_in_nbr+1) ! I [dgr] Interface latitudes
    real,intent(in)::lat_out_grd(lat_out_nbr+1) ! I [dgr] Interface latitudes
    real,intent(in)::lon_in_grd(lon_in_nbr+1) ! I [dgr] Interface longitudes
    real,intent(in)::lon_out_grd(lon_out_nbr+1) ! I [dgr] Interface longitudes
    real,intent(in)::mss_val ! I [frc] Missing value
    real,intent(in)::ovr_wgt(lon_out_nbr,lat_out_nbr,ovr_nbr_max) ! I [frc] Weight of overlapping input gridcells onto each output gridcell
    ! Output
    real,intent(out)::dat_var_out(lon_out_nbr,lat_out_nbr) ! O [frc] Variance of input field 
    ! Local
    integer lat_out_idx ! [idx] Counting index
    integer lon_out_idx ! [idx] Counting index
    integer lat_in_idx ! [idx] Counting index
    integer lon_in_idx ! [idx] Counting index
    integer mss_cnt ! [nbr] Number of missing values in current output gridcell
    integer ovr_idx ! [idx] Counting index
    integer ovr_nbr_crr ! [nbr] Current number of overlapping gridcells
    real ovr_wgt_ttl(lon_out_nbr,lat_out_nbr) ! [frc] Total weight of input cells contributing to output point
    
    ! Main Code
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Entering map_rbn_var()"
    
    ! Initialize arrays
    ovr_wgt_ttl=0.0 ! [frc]
    dat_var_out=0.0 ! [frc]
    
    ! Rebin input data to output grid
    if (mss_flg) then
       do lat_out_idx=1,lat_out_nbr
          do lon_out_idx=1,lon_out_nbr
             ovr_nbr_crr=ovr_nbr(lon_out_idx,lat_out_idx)
             mss_cnt=0
             do ovr_idx=1,ovr_nbr_crr
                lon_in_idx=ovr_lon_idx(lon_out_idx,lat_out_idx,ovr_idx) ! [idx] lon index (input grid) of overlap cell
                lat_in_idx=ovr_lat_idx(lon_out_idx,lat_out_idx,ovr_idx) ! [idx] lat index (input grid) of overlap cell
                if (dat_in(lon_in_idx,lat_in_idx) /= mss_val) then
                   ! ovr_wgt_ttl will be 1.0 except if missing values are present
                   ovr_wgt_ttl(lon_out_idx,lat_out_idx)= &
                        ovr_wgt_ttl(lon_out_idx,lat_out_idx)+ &
                        ovr_wgt(lon_out_idx,lat_out_idx,ovr_idx)
                   dat_var_out(lon_out_idx,lat_out_idx)= &
                        dat_var_out(lon_out_idx,lat_out_idx)+ &
                        ovr_wgt(lon_out_idx,lat_out_idx,ovr_idx)* &
                        (dat_in(lon_in_idx,lat_in_idx)-dat_out(lon_out_idx,lat_out_idx))**2
                else          ! mss_val
                   mss_cnt=mss_cnt+1
                endif         ! mss_val
             end do           ! end loop over ovr
             if (mss_cnt==ovr_nbr_crr) then
                dat_var_out(lon_out_idx,lat_out_idx)=mss_val
             else
                dat_var_out(lon_out_idx,lat_out_idx)= &
                     dat_var_out(lon_out_idx,lat_out_idx)/ovr_wgt_ttl(lon_out_idx,lat_out_idx)
             endif ! endif
          end do ! end loop over lon
       end do ! end loop over lat
    endif ! not mss_flg
    if (.not.mss_flg) then
       do lat_out_idx=1,lat_out_nbr
          do lon_out_idx=1,lon_out_nbr
             ovr_nbr_crr=ovr_nbr(lon_out_idx,lat_out_idx)
             do ovr_idx=1,ovr_nbr_crr
                lon_in_idx=ovr_lon_idx(lon_out_idx,lat_out_idx,ovr_idx) ! [idx] lon index (input grid) of overlap cell
                lat_in_idx=ovr_lat_idx(lon_out_idx,lat_out_idx,ovr_idx) ! [idx] lat index (input grid) of overlap cell
                dat_var_out(lon_out_idx,lat_out_idx)= &
                     dat_var_out(lon_out_idx,lat_out_idx)+ &
                     ovr_wgt(lon_out_idx,lat_out_idx,ovr_idx)* &
                     (dat_in(lon_in_idx,lat_in_idx)-dat_out(lon_out_idx,lat_out_idx))**2
             end do ! end loop over ovr
             dat_var_out(lon_out_idx,lat_out_idx)= &
                  dat_var_out(lon_out_idx,lat_out_idx)/ovr_nbr_crr
          end do ! end loop over lon
       end do ! end loop over lat
    endif ! mss_flg
    
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Exiting map_rbn_var()"
    return 
  end subroutine map_rbn_var ! end map_rbn_var()
  
  subroutine map_ovr_nbr_max_get( &
       lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
       lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
       ovr_nbr_max) ! O
    ! Purpose: Compute maximum number of input cells which overlap any output cell
    ! NB: ovr_nbr==0 is NOT an error for any given gridcell, because the input grid may be staggered from the output grid
    ! We assume ovr_nbr_max would not change (increase) even if the input grid were rotated and map_ovr_nbr_max_get() were called iteratively (twice)
    ! Thus, we currently call map_ovr_nbr_max_get() once to obtain ovr_nbr_max
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use utl_mdl,only:mnt_ncr_chk ! [mdl] Utility functions (date_time_get,mnt_chk...)
    implicit none
    ! Parameters
    ! Commons
    ! Input
    integer,intent(in)::lat_in_nbr        ! [nbr] Number of latitudes
    integer,intent(in)::lat_out_nbr       ! [nbr] Number of latitudes
    integer,intent(in)::lon_in_nbr        ! [nbr] Number of longitudes
    integer,intent(in)::lon_out_nbr       ! [nbr] Number of longitudes
    real,intent(in)::lat_in_grd(lat_in_nbr+1) ! [dgr] Interface latitudes
    real,intent(in)::lat_out_grd(lat_out_nbr+1) ! [dgr] Interface latitudes
    real,intent(in)::lon_in_grd(lon_in_nbr+1) ! [dgr] Interface longitudes
    real,intent(in)::lon_out_grd(lon_out_nbr+1) ! [dgr] Interface longitudes
    ! Output
    integer,intent(out)::ovr_nbr_max       ! [nbr] Maximum number of input cells which overlap any output cell
    ! Local
    integer lat_in_idx        ! [idx] Counting index
    integer lat_out_idx       ! [idx] Counting index
    integer lon_in_idx        ! [idx] Counting index
    integer lon_out_idx       ! [idx] Counting index
    integer ovr_nbr           ! [nbr] Number of input cells which overlap current output cell
    logical mnt_ncr           ! Monotonic and increasing flag
    ! Main Code
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Entering map_ovr_nbr_max_get()"
    
    ! Check for monotonically increasing grids
    mnt_ncr=mnt_ncr_chk(lon_in_grd,lon_in_nbr+1)
    if (.not.mnt_ncr) stop "lon_in_grd not monotonically increasing in map_ovr_nbr_max_get()"
    mnt_ncr=mnt_ncr_chk(lat_in_grd,lat_in_nbr+1)
    if (.not.mnt_ncr) stop "lat_in_grd not monotonically increasing in map_ovr_nbr_max_get()"
    mnt_ncr=mnt_ncr_chk(lon_out_grd,lon_out_nbr+1)
    if (.not.mnt_ncr) stop "lon_out_grd not monotonically increasing in map_ovr_nbr_max_get()"
    mnt_ncr=mnt_ncr_chk(lat_out_grd,lat_out_nbr+1)
    if (.not.mnt_ncr) stop "lat_out_grd not monotonically increasing in map_ovr_nbr_max_get()"
    
    ! Initialize scalars
    ovr_nbr_max=0
    
    ! NB: Since we have that all coordinate arrays are monotonically increasing, and we assume the arrays are regular, the following properties are true:
    ! 1. lon_in_grd(i) is the western edge longitude of cell (i,j)
    ! 2. lat_in_grd(i) is the southern edge latitude of cell (i,j)
    do lat_out_idx=1,lat_out_nbr
       if (dbg_lvl >= dbg_sbr) write (6,"(a1)",advance="no") "."
       do lon_out_idx=1,lon_out_nbr
          ovr_nbr=0           ! Initialize counter for this output gridcell
          do lat_in_idx=1,lat_in_nbr
             ! Check if the latitudes overlap
             ! NB: This math assumes latitude grid obeys -90.0 < lat < 90.0
             if ((lat_in_grd(lat_in_idx) < lat_out_grd(lat_out_idx+1)).and. &
                  (lat_in_grd(lat_in_idx+1) > lat_out_grd(lat_out_idx))) then
                ! Check if the longitudes overlap
                do lon_in_idx=1,lon_in_nbr
                   ! Latitudes already known to overlap, so if longitudes overlap, we have a winner
                   if ((lon_in_grd(lon_in_idx) < lon_out_grd(lon_out_idx+1)).and. &
                        (lon_in_grd(lon_in_idx+1) > lon_out_grd(lon_out_idx))) ovr_nbr=ovr_nbr+1
                end do        ! end loop over lon_in
             endif            ! endif latitudes overlap
          end do              ! end loop over lat_in
          ! Loop over input grid is complete, set maximum overlap number to current overlap number if warranted
          if (ovr_nbr > ovr_nbr_max) ovr_nbr_max=ovr_nbr
       end do ! end loop over lon_out
    end do ! end loop over lat_out

    if (dbg_lvl >= dbg_sbr) write (6,"()") 
    if (ovr_nbr_max==0) then
       stop "ovr_nbr_max==0 in map_ovr_nbr_max_get()"
    endif                     ! endif
    if (dbg_lvl >= dbg_io) write (6,"(a,i3)") "map_ovr_nbr_max_get() reports ovr_nbr_max is ",ovr_nbr_max
    
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Exiting map_ovr_nbr_max_get()"
    return 
  end subroutine map_ovr_nbr_max_get                       ! end map_ovr_nbr_max_get()
  
  subroutine map_ovr_wgt_drv( &
       lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
       lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
       ovr_nbr_max,area_out, & ! I
       ovr_lat_idx,ovr_lon_idx,ovr_nbr,ovr_wgt) ! O
    ! Purpose: Driver routine for map_ovr_wgt_get()
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use utl_mdl,only:mnt_ncr_chk ! [mdl] Utility functions (date_time_get,mnt_chk...)
    implicit none
    ! Parameters
    ! Commons
    ! Input
    integer,intent(in)::lat_in_nbr        ! [nbr] Number of latitudes
    integer,intent(in)::lat_out_nbr       ! [nbr] Number of latitudes
    integer,intent(in)::lon_in_nbr        ! [nbr] Number of longitudes
    integer,intent(in)::lon_out_nbr       ! [nbr] Number of longitudes
    integer,intent(in)::ovr_nbr_max       ! [nbr] Maximum number of input cells which overlap any output cell
    real,intent(in)::area_out(lon_out_nbr,lat_out_nbr) ! [m2] Areas of all gridcells
    real,intent(in)::lat_in_grd(lat_in_nbr+1) ! [dgr] Interface latitudes
    real,intent(in)::lat_out_grd(lat_out_nbr+1) ! [dgr] Interface latitudes
    real,intent(in)::lon_in_grd(lon_in_nbr+1) ! [dgr] Interface longitudes
    real,intent(in)::lon_out_grd(lon_out_nbr+1) ! [dgr] Interface longitudes
    ! Input/Output
    ! Output
    integer,intent(out)::ovr_lat_idx(lon_out_nbr,lat_out_nbr,ovr_nbr_max) ! [idx] Map into input grid of latitude indices of overlap cells
    integer,intent(out)::ovr_lon_idx(lon_out_nbr,lat_out_nbr,ovr_nbr_max) ! [idx] Map into input grid of longitude indices of overlap cells
    integer,intent(out)::ovr_nbr(lon_out_nbr,lat_out_nbr) ! [nbr] Number of input gridcells which overlap each output gridcell
    real,intent(out)::ovr_wgt(lon_out_nbr,lat_out_nbr,ovr_nbr_max) ! [frc] Weight of overlapping input gridcells onto each output gridcell
    ! Local
    real lon_in_grd_rtt(lon_in_nbr+1) ! [dgr] Interface longitudes, rotated
    integer lat_in_idx        ! [idx] Counting index
    integer lat_out_idx       ! [idx] Counting index
    integer lon_in_idx        ! [idx] Counting index
    integer lon_out_idx       ! [idx] Counting index
    integer ovr_idx           ! [idx] Counting index
    integer ovr_nbr_crr       ! [nbr] Current number of overlapping gridcells
    logical mnt_ncr           ! Monotonic and increasing flag
    logical flg_err           ! [flg] Error flag
    real eps_rlt              ! [frc] Relative error allowed in ovr_wgt_ttl
    real lon_rtt              ! [dgr] Longitude rotation amount
    real ovr_wgt_ttl          ! [frc] Total weigt of overlapping points for current gridcell
    ! Main Code
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Entering map_ovr_wgt_drv()"
    
    ! Check for monotonically increasing grids
    mnt_ncr=mnt_ncr_chk(lon_in_grd,lon_in_nbr+1)
    if (.not.mnt_ncr) stop "lon_in_grd not monotonically increasing in map_ovr_wgt_drv()"
    mnt_ncr=mnt_ncr_chk(lat_in_grd,lat_in_nbr+1)
    if (.not.mnt_ncr) stop "lat_in_grd not monotonically increasing in map_ovr_wgt_drv()"
    mnt_ncr=mnt_ncr_chk(lon_out_grd,lon_out_nbr+1)
    if (.not.mnt_ncr) stop "lon_out_grd not monotonically increasing in map_ovr_wgt_drv()"
    mnt_ncr=mnt_ncr_chk(lat_out_grd,lat_out_nbr+1)
    if (.not.mnt_ncr) stop "lat_out_grd not monotonically increasing in map_ovr_wgt_drv()"
    
    ! Initialize scalars
    eps_rlt=1.0e-5            ! [frc] Relative error allowed in ovr_wgt_ttl
    
    ! Initialize arrays
    do lon_out_idx=1,lon_out_nbr
       do lat_out_idx=1,lat_out_nbr
          ovr_nbr(lon_out_idx,lat_out_idx)=0  
          do ovr_idx=1,ovr_nbr_max
             ovr_lat_idx(lon_out_idx,lat_out_idx,ovr_idx)=0
             ovr_lon_idx(lon_out_idx,lat_out_idx,ovr_idx)=0
             ovr_wgt(lon_out_idx,lat_out_idx,ovr_idx)=0.0
          end do ! end loop over ovr
       end do ! end loop over lat
    end do ! end loop over lon
    
    call map_ovr_wgt_get( &
         lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
         lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
         ovr_nbr_max,area_out, & ! I
         ovr_lat_idx,ovr_lon_idx,ovr_nbr,ovr_wgt) ! I/O
    
    ! map_ovr_wgt_get() does not report all overlaps whenever the input and output grids are staggered in longitude
    ! This is because lon_out_grd may contain, e.g., [0.0,360.0] and lon_in_grd may be [-180.0,180.0]
    ! When this occurs then map_ovr_wgt_get() will not find any of the overlap gridcells in [-180.0,0.0]
    ! A nice solution is to shift lon_in_grd by 360.0 degrees and call map_ovr_wgt_get() again, adding on to the counts and values in the original overlap arrays
    ! The main purpose of map_ovr_wgt_drv() is, in fact, to hide the details of this solution from the user
    ! When the input and output grids are identical, then the first call to map_ovr_wgt_get() finds all the overlaps
    ! In this case, rotating the input grid 360 degrees and calling map_ovr_wgt_get() does no harm because there will be no overlaps on the second call
    
    ! Rotate input grid
    if (lon_in_grd(1) < lon_out_grd(1)) then 
       lon_rtt=360.0          ! [dgr]
    else 
       lon_rtt=-360.0         ! [dgr]
    endif                     ! endif
    do lon_in_idx=1,lon_in_nbr+1
       lon_in_grd_rtt(lon_in_idx)=lon_in_grd(lon_in_idx)+lon_rtt ! [dgr] Interface longitudes, rotated
    end do                     ! end loop over lon
    
    call map_ovr_wgt_get( &
         lat_in_grd,lat_in_nbr,lon_in_grd_rtt,lon_in_nbr, & ! I
         lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
         ovr_nbr_max,area_out, & ! I
         ovr_lat_idx,ovr_lon_idx,ovr_nbr,ovr_wgt) !Input/Output

!  write (88,'(10I5)')((ovr_nbr(lon_out_idx,lat_out_idx),lon_out_idx=1,lon_out_nbr),lat_out_idx=1,lat_out_nbr)

    ! Sanity check
    do lon_out_idx=1,lon_out_nbr
       do lat_out_idx=1,lat_out_nbr
          flg_err=.false.
          ovr_wgt_ttl=0.0
          ovr_nbr_crr=ovr_nbr(lon_out_idx,lat_out_idx)
          if (ovr_nbr_crr <= 0) flg_err=.true.
          do ovr_idx=1,ovr_nbr_crr
             lon_in_idx=ovr_lon_idx(lon_out_idx,lat_out_idx,ovr_idx) ! [idx] lon index (input grid) of overlap cell
             lat_in_idx=ovr_lat_idx(lon_out_idx,lat_out_idx,ovr_idx) ! [idx] lat index (input grid) of overlap cell
             ovr_wgt_ttl=ovr_wgt_ttl+ovr_wgt(lon_out_idx,lat_out_idx,ovr_idx) ! weight of overlap cell
             if (lon_in_idx <= 0.or.lat_in_idx <= 0 ) flg_err=.true.
          end do ! end loop over ovr
          ! Sum of overlap weights should equal unity for each output gridcell
          if (abs(ovr_wgt_ttl-1.0) >= eps_rlt) flg_err=.true.
          if (flg_err) then
             write (6,"(a)") "map_ovr_wgt_drv(): Error flagged"
             write (6,"(a,i7)") "map_ovr_wgt_drv(): # overlapping input gridcells = ",ovr_nbr_crr
             write (6,"(a,f9.3)") "map_ovr_wgt_drv(): Weights of input cells sum to = ",ovr_wgt_ttl
             write (6,"(a)") "Output cell edge locations:"
             write (6,"(a)") "(idx)   Lat Sth (idx)   Lat Nrt (idx)   Lon Wst (idx)   Lon Est"
             write (6,"(4(a1,i3,a1,1x,f9.4,1x))") &
                  "(",lat_out_idx,")",lat_out_grd(lat_out_idx),"(",lat_out_idx+1,")",lat_out_grd(lat_out_idx+1), &
                  "(",lon_out_idx,")",lon_out_grd(lon_out_idx),"(",lon_out_idx+1,")",lon_out_grd(lon_out_idx+1)
             write (6,"(a)") "Overlapping points and weights:"
             write (6,"(a)") "(idx)   Lat Sth (idx)   Lat Nrt (idx)   Lon Wst (idx)   Lon Est"
             do ovr_idx=1,ovr_nbr_crr ! Overlap cell index
                lon_in_idx=ovr_lon_idx(lon_out_idx,lat_out_idx,ovr_idx) ! [idx] lon index (input grid) of overlap cell
                lat_in_idx=ovr_lat_idx(lon_out_idx,lat_out_idx,ovr_idx) ! [idx] lat index (input grid) of overlap cell
                write (6,"(4(a1,i3,a1,1x,f9.4,1x),f9.4)") &
                     "(",lat_in_idx,")",lat_in_grd(lat_in_idx),"(",lat_in_idx+1,")",lat_in_grd(lat_in_idx+1), &
                     "(",lon_in_idx,")",lon_in_grd(lon_in_idx),"(",lon_in_idx+1,")",lon_in_grd(lon_in_idx+1), &
                     ovr_wgt(lon_out_idx,lat_out_idx,ovr_idx)
             end do ! end loop over overlapping cells
             stop
          endif ! endif err
       end do ! end loop over lat
    end do ! end loop over lon
    
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Exiting map_ovr_wgt_drv()"
    return 
  end subroutine map_ovr_wgt_drv                       ! end map_ovr_wgt_drv()
  
  subroutine map_ovr_wgt_get( &
       lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
       lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
       ovr_nbr_max,area_out, & ! I
       ovr_lat_idx,ovr_lon_idx,ovr_nbr,ovr_wgt) ! I/O
    ! Purpose: Compute locations, numbers, and weights of input cells which overlap each output cell
    ! map_ovr_wgt_get() should normally be called by map_ovr_wgt_drv(), which performs rigorous error checking
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use map_cst ! [mdl] Constants used in map routines
    use utl_mdl,only:mnt_ncr_chk ! [mdl] Utility functions (date_time_get,mnt_chk...)
    implicit none
    ! Parameters
    ! Commons
    ! Input
    integer,intent(in)::lat_in_nbr        ! [nbr] Number of latitudes
    integer,intent(in)::lat_out_nbr       ! [nbr] Number of latitudes
    integer,intent(in)::lon_in_nbr        ! [nbr] Number of longitudes
    integer,intent(in)::lon_out_nbr       ! [nbr] Number of longitudes
    integer,intent(in)::ovr_nbr_max       ! [nbr] Maximum number of input cells which overlap any output cell
    real,intent(in)::area_out(lon_out_nbr,lat_out_nbr) ! [m2] Areas of all gridcells
    real,intent(in)::lat_in_grd(lat_in_nbr+1) ! [dgr] Interface latitudes
    real,intent(in)::lat_out_grd(lat_out_nbr+1) ! [dgr] Interface latitudes
    real,intent(in)::lon_in_grd(lon_in_nbr+1) ! [dgr] Interface longitudes
    real,intent(in)::lon_out_grd(lon_out_nbr+1) ! [dgr] Interface longitudes
    ! Input/Output
    integer,intent(inout)::ovr_lat_idx(lon_out_nbr,lat_out_nbr,ovr_nbr_max) ! [idx] Map into input grid of latitude indices of overlap cells
    integer,intent(inout)::ovr_lon_idx(lon_out_nbr,lat_out_nbr,ovr_nbr_max) ! [idx] Map into input grid of longitude indices of overlap cells
    integer,intent(inout)::ovr_nbr(lon_out_nbr,lat_out_nbr) ! [nbr] Number of input gridcells which overlap each output gridcell
    real,intent(inout)::ovr_wgt(lon_out_nbr,lat_out_nbr,ovr_nbr_max) ! [frc] Weight of overlapping input gridcells onto each output gridcell
    ! Output
    ! Local
    real(selected_real_kind(p=12))::dgr2rdn  ! Convert degrees to radians
    real(selected_real_kind(p=12))::pi       ! [frc] 3
    integer lat_in_idx        ! [idx] Counting index
    integer lat_out_idx       ! [idx] Counting index
    integer lon_in_idx        ! [idx] Counting index
    integer lon_out_idx       ! [idx] Counting index
    integer ovr_idx           ! [idx] Counting index
    logical mnt_ncr           ! Monotonic and increasing flag
    real lat_sin_dlt          ! [dgr] Latitude overlap
    real ovr_lat_nrt          ! [dgr] Northern boundary of overlap
    real ovr_lat_sth          ! [dgr] Southern boundary of overlap
    real ovr_lon_dgr          ! [dgr] Longitude overlap
    real ovr_lon_est          ! [dgr] Eastern boundary of overlap
    real ovr_lon_rdn          ! [rdn] Longitude overlap
    real ovr_lon_wst          ! [dgr] Western boundary of overlap
    ! Main Code
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Entering map_ovr_wgt_get()"
    
    ! Check for monotonically increasing grids
    mnt_ncr=mnt_ncr_chk(lon_in_grd,lon_in_nbr+1)
    if (.not.mnt_ncr) stop "lon_in_grd not monotonically increasing in map_ovr_wgt_get()"
    mnt_ncr=mnt_ncr_chk(lat_in_grd,lat_in_nbr+1)
    if (.not.mnt_ncr) stop "lat_in_grd not monotonically increasing in map_ovr_wgt_get()"
    mnt_ncr=mnt_ncr_chk(lon_out_grd,lon_out_nbr+1)
    if (.not.mnt_ncr) stop "lon_out_grd not monotonically increasing in map_ovr_wgt_get()"
    mnt_ncr=mnt_ncr_chk(lat_out_grd,lat_out_nbr+1)
    if (.not.mnt_ncr) stop "lat_out_grd not monotonically increasing in map_ovr_wgt_get()"
    
    ! Initialize scalars
    pi=4.0*atan(1.0d0) ! [frc] 3
    dgr2rdn=pi/180.0d0        ! rdn dgr-1
    
    ! All coordinate arrays are monotonically increasing, and we assume the arrays are regular, so the following properties are true:
    ! 1. lon_in_grd(i) is the western edge longitude of cell (i,j)
    ! 2. lat_in_grd(i) is the southern edge latitude of cell (i,j)
    
    !$omp parallel default(none)
    !$omp single
    !$write (6,"(a,a36,i3,a8)") prg_nm(1:ftn_strlen(prg_nm)), & 
    !$": INFO OpenMP multi-threading using ",omp_get_num_threads()," threads"
    !$omp end single
    !$omp end parallel
    
    !$omp parallel default(none) &
    !$omp private(lat_out_idx,lon_out_idx,lat_in_idx,lon_in_idx) &
    !$omp private(ovr_lat_nrt,ovr_lat_sth,ovr_lon_est,ovr_lon_wst,ovr_lon_dgr,ovr_lon_rdn,ovr_idx) &
    !$omp private(lat_sin_dlt) &
    !$omp shared(lat_out_nbr,lon_out_nbr,lat_in_nbr,lon_in_nbr) &
    !$omp shared(lat_out_grd,lon_out_grd,lat_in_grd,lon_in_grd) &
    !$omp shared(ovr_nbr,ovr_lon_idx,ovr_lat_idx,ovr_wgt,ovr_nbr_max,area_out) &
    !$omp shared(dgr2rdn) &
    !$omp shared(dbg_lvl)
    
    !$omp do

    do lat_out_idx=1,lat_out_nbr
       if (dbg_lvl >= dbg_sbr) write (6,"(a1)",advance="no") "."
       do lat_in_idx=1,lat_in_nbr
          ! Check if the latitudes overlap
          lat_sin_dlt=0.0
          ! NB: This math assumes latitude grid obeys -90.0 < lat < 90.0
          if ((lat_in_grd(lat_in_idx) < lat_out_grd(lat_out_idx+1)).and. &
               (lat_in_grd(lat_in_idx+1) > lat_out_grd(lat_out_idx))) then
             ovr_lat_nrt=min(lat_in_grd(lat_in_idx+1),lat_out_grd(lat_out_idx+1))
             ovr_lat_sth=max(lat_in_grd(lat_in_idx),lat_out_grd(lat_out_idx))
             lat_sin_dlt=sin(ovr_lat_nrt*dgr2rdn)-sin(ovr_lat_sth*dgr2rdn)
          endif               ! endif latitudes overlap
          if (lat_sin_dlt > 0.0) then
             ! Check if the longitudes overlap
             do lon_out_idx=1,lon_out_nbr
                ! Uncomment this block to look at a particular cell
                ! if (lon_out_idx==1.and.lat_out_idx==1) then
                ! write (6,"(a)") "Cell edge locations:"
                ! write (6,"(a)") "(idx)   Lat Sth (idx)   Lat Nrt (idx)   Lon Wst (idx)   Lon Est"
                ! write (6,"(4(a1,i3,a1,1x,f9.4,1x))") 
                ! $                    "(",lat_out_idx,")",lat_out_grd(lat_out_idx),"(",lat_out_idx+1,")",lat_out_grd(lat_out_idx+1), &
                ! $                    "(",lon_out_idx,")",lon_out_grd(lon_out_idx),"(",lon_out_idx+1,")",lon_out_grd(lon_out_idx+1)
                ! endif         ! endif dbg
                do lon_in_idx=1,lon_in_nbr
                   ovr_lon_dgr=0.0
                   ! Latitudes already known to overlap, so, if longitudes overlap, we have a winner
                   if ((lon_in_grd(lon_in_idx) < lon_out_grd(lon_out_idx+1)).and. &
                        (lon_in_grd(lon_in_idx+1) > lon_out_grd(lon_out_idx))) then
                      ovr_lon_est=min(lon_in_grd(lon_in_idx+1),lon_out_grd(lon_out_idx+1))
                      ovr_lon_wst=max(lon_in_grd(lon_in_idx),lon_out_grd(lon_out_idx))
                      ovr_lon_dgr=ovr_lon_est-ovr_lon_wst
                   endif      ! endif latitudes overlap
                   if (ovr_lon_dgr > 0.0) then
                      ovr_nbr(lon_out_idx,lat_out_idx)=ovr_nbr(lon_out_idx,lat_out_idx)+1
                      ovr_idx=ovr_nbr(lon_out_idx,lat_out_idx)
                      ! Sanity check OK but AIX does not compile stop"s in OMP loops
                      ! if (ovr_idx > ovr_nbr_max) stop "ovr_idx > ovr_nbr_max in map_ovr_wgt_get()"
                      ovr_lat_idx(lon_out_idx,lat_out_idx,ovr_idx)=lat_in_idx
                      ovr_lon_idx(lon_out_idx,lat_out_idx,ovr_idx)=lon_in_idx
                      ovr_lon_rdn=ovr_lon_dgr*dgr2rdn
                      ! Overlap weights are normalized to sum to unity for each output gridcell
                      ovr_wgt(lon_out_idx,lat_out_idx,ovr_idx)= &
                           earth_rds_sqr*lat_sin_dlt*ovr_lon_rdn/ &
                           area_out(lon_out_idx,lat_out_idx)
                      ! Uncomment this block to look at a particular cell
                      ! if (lon_out_idx==1.and.lat_out_idx==1) then
                      ! write (6,"(a)") "(idx)   Lat Sth (idx)   Lat Nrt (idx)   Lon Wst (idx)   Lon Est"
                      ! write (6,"(4(a1,i3,a1,1x,f9.4,1x),f9.4)") 
                      ! $                          "(",lat_in_idx,")",lat_in_grd(lat_in_idx),"(",lat_in_idx+1,")",lat_in_grd(lat_in_idx+1), &
                      ! $                          "(",lon_in_idx,")",lon_in_grd(lon_in_idx),"(",lon_in_idx+1,")",lon_in_grd(lon_in_idx+1), &
                      ! $                          ovr_wgt(lon_out_idx,lat_out_idx,ovr_idx)
                      ! endif ! endif dbg
                   endif ! endif cell overlaps
                end do ! end loop over lon_in
             end do ! end loop over lon_out
          endif ! endif latitudes overlap
       end do ! end loop over lat_in
    end do ! end loop over lat_out
    !$omp end do
    !$omp end parallel

    if (dbg_lvl >= dbg_sbr) write (6,"()") ! This puts a newline after the dots
    
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Exiting map_ovr_wgt_get()"
    return 
  end subroutine map_ovr_wgt_get ! end map_ovr_wgt_get()
  
  subroutine map_area_get(lat_grd,lat_nbr,lon_grd,lon_nbr, & ! I
       map_area) ! O
    ! Purpose: Compute map area of each cell
    ! Entire computation should be performed in double precision if possible
    ! because errors may accumulate for large input grids with O(10^6) cells
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use map_cst ! [mdl] Constants used in map routines
    use utl_mdl,only:mnt_ncr_chk ! [mdl] Utility functions (date_time_get,mnt_chk...)
    implicit none
    ! Parameters
    real,parameter::eps_rlt=1.0e-5 ! [frc] Relative error allowed in area_ttl
    ! Commons
    ! Input
    integer,intent(in)::lon_nbr           ! [nbr] Number of longitudes
    integer,intent(in)::lat_nbr           ! [nbr] Number of latitudes
    real,intent(in)::lon_grd(lon_nbr+1) ! [dgr] Interface longitudes
    real,intent(in)::lat_grd(lat_nbr+1) ! [dgr] Interface latitudes
    ! Output
    real,intent(out)::map_area(lon_nbr,lat_nbr) ! [m2] Areas of all gridcells
    ! Local
    real(selected_real_kind(p=12))::area_ttl ! [m2] Sum of map_area array
    real(selected_real_kind(p=12))::dgr2rdn  ! [frc] Convert degrees to radians
    real(selected_real_kind(p=12))::earth_area ! [m2] Area of Earth
    real(selected_real_kind(p=12))::lat_sin_dlt ! [rdn] Cell latitude size
    real(selected_real_kind(p=12))::lon_dlt_rdn ! [rdn] Cell longitude size
    real(selected_real_kind(p=12))::pi ! [frc] 3
    integer lat_idx ! [idx] Counting index
    integer lon_idx ! [idx] Counting index
    logical mnt_ncr ! [flg] Monotonicity flag (increasing)
    ! Main Code
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Entering map_area_get()"
    
    ! Check for monotonicity and increasing
    mnt_ncr=mnt_ncr_chk(lon_grd,lon_nbr+1)
    if (.not.mnt_ncr) stop "lon_grd not monotonically increasing in map_area_get()"
    mnt_ncr=mnt_ncr_chk(lat_grd,lat_nbr+1)
    if (.not.mnt_ncr) stop "lat_grd not monotonically increasing in map_area_get()"
    
    ! Initialize scalars
    pi=4.0*atan(1.0d0) ! [frc] 3
    dgr2rdn=pi/180.0d0 ! [rdn dgr-1]
    area_ttl=0.0 ! [m2]
    earth_area=4.0*pi*earth_rds_sqr ! [m2]
    
    ! Initialize arrays
    do lon_idx=1,lon_nbr
       do lat_idx=1,lat_nbr
          ! Assume lon and lat are in degrees and monotonically increasing
          lon_dlt_rdn=(lon_grd(lon_idx+1)-lon_grd(lon_idx))*dgr2rdn ! [rdn]
          lat_sin_dlt=sin(lat_grd(lat_idx+1)*dgr2rdn)-sin(lat_grd(lat_idx)*dgr2rdn)
          map_area(lon_idx,lat_idx)=earth_rds_sqr*lon_dlt_rdn*lat_sin_dlt ! [m2]
          area_ttl=area_ttl+map_area(lon_idx,lat_idx) ! [m2]
       end do ! end loop over lat
    end do ! end loop over lon
    
    ! Sanity check
    if (abs(area_ttl-earth_area) > eps_rlt*earth_area) then
       write (6,"(3(a,es14.7))") &
            "map_area_get(): area_ttl = ",area_ttl, &
            " m2 and earth_area = ",earth_area, &
            " m2 differ by a factor of ",area_ttl/earth_area
       ! Analyze map grid
       call map_grd_anl(lat_grd,lat_nbr,lon_grd,lon_nbr)
       stop "area_ttl not equal to earth_area in map_area_get()"
    endif                     ! not insane
    
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Exiting map_area_get()"
    return 
  end subroutine map_area_get                       ! end map_area_get()
  
  subroutine map_grd_anl(lat_grd,lat_nbr,lon_grd,lon_nbr) ! I
    ! Purpose: Analyze map grid
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use map_cst ! [mdl] Constants used in map routines
    use sng_mdl,only:ftn_strlen ! [mdl] String manipulation
    use utl_mdl,only:mnt_ncr_chk ! [mdl] Utility functions (date_time_get,mnt_chk...)
    implicit none
    ! Parameters
    character(*),parameter::sbr_nm="map_grd_anl()" ! [sng] Subroutine name
    ! Commons
    ! Input
    integer,intent(in)::lon_nbr ! [nbr] Number of longitudes
    integer,intent(in)::lat_nbr ! [nbr] Number of latitudes
    real,intent(in)::lon_grd(lon_nbr+1) ! [dgr] Interface longitudes
    real,intent(in)::lat_grd(lat_nbr+1) ! [dgr] Interface latitudes
    ! Output
    ! Local
    real map_area(lon_nbr,lat_nbr) ! [m2] Gridcell area
    real(selected_real_kind(p=12))::area_ttl ! [m2] Sum of map_area array
    real(selected_real_kind(p=12))::dgr2rdn  ! [frc] Convert degrees to radians
    real(selected_real_kind(p=12))::earth_area ! [m2] Area of Earth
    real(selected_real_kind(p=12))::lat_sin_dlt ! [rdn] Cell latitude size
    real(selected_real_kind(p=12))::lon_dlt_rdn ! [rdn] Cell longitude size
    real(selected_real_kind(p=12))::pi ! [frc] 3
    integer lat_idx ! [idx] Counting index
    integer lon_idx ! [idx] Counting index
    logical mnt_ncr ! [flg] Monotonicity flag (increasing)
    real eps_rlt ! [frc] Relative error allowed in ovr_wgt_ttl
    real lat_ctr(lat_nbr) ! [dgr] Latitude center
    real lat_dlt(lat_nbr) ! [dgr] Latitude span
    real lon_ctr(lon_nbr) ! [dgr] Longitude center
    real lon_dlt(lon_nbr) ! [dgr] Longitude span
    ! Main Code
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Entering "//sbr_nm
    
    ! Check for monotonicity and increasing
    mnt_ncr=mnt_ncr_chk(lon_grd,lon_nbr+1)
    if (.not.mnt_ncr) then
       write (6,"(2a)") prg_nm(1:ftn_strlen(prg_nm)), & 
            "lon_grd not monotonically increasing in "//sbr_nm
    endif ! endif err
    mnt_ncr=mnt_ncr_chk(lat_grd,lat_nbr+1)
    if (.not.mnt_ncr) then
       write (6,"(2a)") prg_nm(1:ftn_strlen(prg_nm)), & 
            "lat_grd not monotonically increasing in "//sbr_nm
    endif ! endif err
    
    ! Initialize scalars
    pi=4.0*atan(1.0d0) ! [frc] 3
    dgr2rdn=pi/180.0d0 ! [rdn dgr-1]
    area_ttl=0.0 ! [m2]
    earth_area=4.0*pi*earth_rds_sqr ! [m2]
    ! Entire computation should be performed in double precision if possible
    ! because errors may accumulate for large input grids with O(10^6) cells
    eps_rlt=1.0e-3 ! [frc] Relative error allowed in area_ttl
    
    ! Initialize arrays
    do lon_idx=1,lon_nbr
       do lat_idx=1,lat_nbr
          ! Assume lon and lat are in degrees and monotonically increasing
          lon_dlt(lon_idx)=lon_grd(lon_idx+1)-lon_grd(lon_idx) ! [dgr] Longitude span
          lat_dlt(lat_idx)=lat_grd(lat_idx+1)-lat_grd(lat_idx) ! [dgr] Latitude span
          lon_ctr(lon_idx)=0.5*(lon_grd(lon_idx)+lon_grd(lon_idx+1)) ! [dgr] Longitude center
          lat_ctr(lat_idx)=0.5*(lat_grd(lat_idx)+lat_grd(lat_idx+1)) ! [dgr] Latitude center
          lon_dlt_rdn=(lon_grd(lon_idx+1)-lon_grd(lon_idx))*dgr2rdn ! [rdn]
          lat_sin_dlt=sin(lat_grd(lat_idx+1)*dgr2rdn)-sin(lat_grd(lat_idx)*dgr2rdn)
          map_area(lon_idx,lat_idx)=earth_rds_sqr*lon_dlt_rdn*lat_sin_dlt ! [m2]
          area_ttl=area_ttl+map_area(lon_idx,lat_idx) ! [m2]
       end do ! end loop over lat
    end do ! end loop over lon
    
    ! Latitude
    write (6,"(a)") "Latitude grid properties:"
    write (6,"(a)") "Cell edge locations, center, and spacing:"
    write (6,"(a)") "(idx)   Lat Wst (idx)   Lat Est   lat_ctr   lat_dlt"
    do lat_idx=1,lat_nbr
       write (6,"(2(a1,i3,a1,1x,f9.4,1x),1(f9.4,1x,f8.5,1x))") &
            "(",lat_idx,")",lat_grd(lat_idx),"(",lat_idx+1,")",lat_grd(lat_idx+1), &
            lat_ctr(lat_idx),lat_dlt(lat_idx)
    end do ! end loop over lat
    
    ! Longitude
    write (6,"(a)") "Longitude grid properties:"
    write (6,"(a)") "Cell edge locations, center, and spacing:"
    write (6,"(a)") "(idx)   Lon Wst (idx)   Lon Est   lon_ctr   lon_dlt"
    do lon_idx=1,lon_nbr
       write (6,"(2(a1,i3,a1,1x,f9.4,1x),1(f9.4,1x,f8.5,1x))") &
            "(",lon_idx,")",lon_grd(lon_idx),"(",lon_idx+1,")",lon_grd(lon_idx+1), &
            lon_ctr(lon_idx),lon_dlt(lon_idx)
    end do ! end loop over lon
    
    ! Gridcells 
    ! Printing entire 2-D map array takes copious space/time for high resolution grids
    ! Presumably this information is most useful for irregular grids
    if (dbg_lvl >= dbg_io) then
       write (6,"(a)") "Map grid properties:"
       write (6,"(a)") "Input cell edge locations:"
       write (6,"(a)") "(idx)   Lat Sth (idx)   Lat Nrt (idx)   Lon Wst (idx)   Lon Est   lat_dlt   lon_dlt"
       do lat_idx=1,lat_nbr
          do lon_idx=1,lon_nbr
             write (6,"(4(a1,i3,a1,1x,f9.4,1x),2(f8.5,1x))") &
                  "(",lat_idx,")",lat_grd(lat_idx),"(",lat_idx+1,")",lat_grd(lat_idx+1), &
                  "(",lon_idx,")",lon_grd(lon_idx),"(",lon_idx+1,")",lon_grd(lon_idx+1), &
                  lat_dlt(lat_idx),lon_dlt(lon_idx)
          end do ! end loop over lon
       end do ! end loop over lat
    endif ! endif dbg
    
    ! Sanity check
    if (abs(area_ttl-earth_area) > eps_rlt*earth_area) then
       write (6,"(3(a,es14.7))") &
            "map_grd_anl(): area_ttl = ",area_ttl, &
            " m2 and earth_area = ",earth_area, &
            " m2 differ by a factor of ",area_ttl/earth_area
       ! Analyzed map grid
       stop "area_ttl not equal to earth_area in map_grd_anl()"
    endif ! not insane
    
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Exiting "//sbr_nm
    return 
  end subroutine map_grd_anl ! end map_grd_anl()
  
  subroutine map_grd_read(fl_in, & ! I
       lat_grd,lat_out_nbr,lon_grd,lon_out_nbr) ! O
    ! Purpose: Return the map grid to the calling program
    ! NB: This routine is intended to work on netCDF files that only have lat,lon coordinates
    ! lat_grd and lon_grd are generated from these coordinates
    use netcdf ! [mdl] netCDF interface
    use nf90_utl ! [mdl] netCDF utilities
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use sng_mdl,only:ftn_strlen ! [mdl] String manipulation
    use utl_mdl,only:mnt_chk ! [mdl] Utility functions (date_time_get,mnt_chk...)
    implicit none
    ! Parameters
    ! Commons
    ! Input
    character(*),intent(in)::fl_in       ! [sng] Name of netCDF input file
    ! Input/Output
    integer,intent(inout)::lon_out_nbr       ! [nbr] Number of longitudes
    integer,intent(inout)::lat_out_nbr       ! [nbr] Number of latitudes
    ! Output
    real,intent(out)::lat_grd(lat_out_nbr+1) ! [dgr] Interface latitudes
    real,intent(out)::lon_grd(lon_out_nbr+1) ! [dgr] Interface longitudes
    ! Local
    integer idx               ! [idx] Counting index
    integer lat_dmn_id        ! [id] Dimension ID for lat
    integer lat_id            ! [id] Coordinate ID for lat
    integer lat_in_nbr        ! [nbr] Number of latitudes
    integer lat_nbr           ! [nbr] Number of latitudes
    integer lon_dmn_id        ! [id] Dimension ID for lon
    integer lon_id            ! [id] Coordinate ID for lon
    integer lon_in_nbr        ! [nbr] Number of longitudes
    integer lon_nbr           ! [nbr] Number of longitudes
    integer nc_id             ! netCDF file ID of output file
    integer rcd               ! [enm] Return success code
    logical mnt               ! [flg] Monotonicity flag
    real lat(lat_out_nbr) ! Nominal latitude coordinate
    real lat_max(lat_out_nbr) ! Maximum latitude in bin
    real lat_min(lat_out_nbr) ! Minimum latitude in bin
    real lat_pad(0:lat_out_nbr+1) ! Padded longitude coordinate
    real lat_sz(lat_out_nbr) ! Latitudinal size of bin
    real lon(lon_out_nbr) ! Nominal longitude coordinate
    real lon_max(lon_out_nbr) ! Maximum longitude in bin
    real lon_min(lon_out_nbr) ! Minimum longitude in bin
    real lon_pad(0:lon_out_nbr+1) ! Padded longitude coordinate
    real lon_sz(lon_out_nbr) ! Longitudinal size of bin
    
    ! Main Code
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Entering map_grd_read()"
    
    ! Read netCDF data
    rcd=nf90_noerr              ! nf90_noerr == 0
    rcd=rcd+nf90_wrp_open(fl_in,nf90_nowrite,nc_id)
    ! Get dimension IDs
    rcd=nf90_wrp_inq_dimid(nc_id,"lat",lat_dmn_id)
    rcd=nf90_wrp_inq_dimid(nc_id,"lon",lon_dmn_id)
    ! Get dimension sizes
    rcd=rcd+nf90_inquire_dimension(nc_id,lat_dmn_id,len=lat_in_nbr)
    rcd=rcd+nf90_inquire_dimension(nc_id,lon_dmn_id,len=lon_in_nbr)
    ! Enough memory? 
    if (lat_in_nbr > lat_out_nbr) stop "lat_in_nbr > lat_out_nbr in map_grd_read()"
    if (lon_in_nbr > lon_out_nbr) stop "lon_in_nbr > lon_out_nbr in map_grd_read()"
    ! Get variable IDs
    rcd=nf90_wrp_inq_varid(nc_id,"lat",lat_id)
    rcd=nf90_wrp_inq_varid(nc_id,"lon",lon_id)
    ! Get data
    rcd=nf90_wrp(nf90_get_var(nc_id,lon_id,lon),"get_var lon")
    rcd=nf90_wrp(nf90_get_var(nc_id,lat_id,lat),"get_var lat")
    rcd=rcd+nf90_close(nc_id)
    if (rcd /= nf90_noerr) then
       write (6,"(4a,i4)") prg_nm(1:ftn_strlen(prg_nm)), & 
            ": ERROR ingesting ",fl_in(1:ftn_strlen(fl_in))," in map_grd_read(), rcd = ",rcd
       stop
    endif                     ! endif err
    if (dbg_lvl >= dbg_fl) write (6,"(a,1x,a)") "Read map grid from netCDF file",fl_in(1:ftn_strlen(fl_in))
    
    ! Set the output dimension sizes to the size of the input dimensions
    lon_out_nbr=lon_in_nbr
    lat_out_nbr=lat_in_nbr
    lon_nbr=lon_in_nbr
    lat_nbr=lat_in_nbr
    
    ! Pad longitude array
    do idx=1,lon_out_nbr 
       lon_pad(idx)=lon(idx)
    end do                     ! end loop over lon
    if (abs(lon_pad(1)-lon_pad(lon_nbr)) > 180.) then
       lon_pad(0)=lon_pad(lon_nbr)-360. 
       lon_pad(lon_nbr+1)=lon_pad(1)+360.
    else 
       lon_pad(0)=lon_pad(1)-(lon_pad(2)-lon_pad(1))
       lon_pad(lon_nbr+1)=lon_pad(lon_nbr)+(lon_pad(lon_nbr)-lon_pad(lon_nbr-1))
    endif                     ! endif
    ! Compute longitude grid
    do idx=1,lon_out_nbr 
       lon_min(idx)=0.5*(lon_pad(idx-1)+lon_pad(idx))
       lon_max(idx)=0.5*(lon_pad(idx)+lon_pad(idx+1))
       lon_sz(idx)=lon_max(idx)-lon_min(idx)
       lon_grd(idx)=lon_min(idx)
    end do                     ! end loop over lon
    lon_grd(lon_nbr+1)=lon_max(lon_nbr)
    
    ! Pad latitude array
    do idx=1,lat_out_nbr 
       lat_pad(idx)=lat(idx)
    end do                     ! end loop over lat
    lat_pad(0)=lat_pad(1)-2.0*(lat_pad(1)-(-90.0))
    lat_pad(lat_nbr+1)=lat_pad(lat_nbr)+2.0*(90.0-lat_pad(lat_nbr))
    ! Compute latitude grid
    do idx=1,lat_out_nbr 
       lat_min(idx)=0.5*(lat_pad(idx-1)+lat_pad(idx))
       lat_max(idx)=0.5*(lat_pad(idx)+lat_pad(idx+1))
       lat_sz(idx)=lat_max(idx)-lat_min(idx)
       lat_grd(idx)=lat_min(idx)
    end do                     ! end loop over lat
    lat_grd(lat_nbr+1)=lat_max(lat_nbr)
    
    ! Check for monotonicity
    mnt=mnt_chk(lon_grd,lon_nbr+1)
    if (.not.mnt) stop "lon_grd not monotonic in map_grd_read()"
    mnt=mnt_chk(lat_grd,lat_nbr+1)
    if (.not.mnt) stop "lat_grd not monotonic in map_grd_read()"
    
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Exiting map_grd_read()"
    return 
  end subroutine map_grd_read                       ! end map_grd_read()
  
  subroutine map_grd_mk(lat_nbr,lon_nbr,map_typ, & ! I
       lat_grd,lon_grd) ! O
    ! Purpose: Construct and return a map grid based on given grid characteristics
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use map_cst ! [mdl] Constants used in map routines
    use sng_mdl,only:ftn_strlen ! [mdl] String manipulation
    implicit none
    ! Parameters
    ! Commons
    ! Input
    integer,intent(in)::lon_nbr           ! [nbr] Number of longitudes
    integer,intent(in)::lat_nbr           ! [nbr] Number of latitudes
    integer,intent(in)::map_typ           ! [enm] Map type in LSM "mn" format
    ! Output
    real,intent(out)::lon_grd(lon_nbr+1) ! [dgr] Interface longitudes
    real,intent(out)::lat_grd(lat_nbr+1) ! [dgr] Interface latitudes
    ! Local
    real(selected_real_kind(p=12))::pi       ! [frc] 3
    real(selected_real_kind(p=12))::rdn2dgr  ! Convert radians to degrees
    integer lat_idx           ! [idx] Counting index
    integer lon_idx           ! [idx] Counting index
    integer map_grd_m         ! [enm] Latitude component of map type
    integer map_grd_n         ! [enm] Longitude component of map type
    integer map_lon_srt_typ   ! Location of first longitude in each latitude band
    integer map_lon_ctr_typ   ! Centering type of first longitude
    integer map_lat_grd_typ   ! Latitude grid type
    real lat(lat_nbr) ! Centered latitudes
    real lat_ncr              ! [dgr] Latitude increment
    real lat_sin(lat_nbr) ! Sine of Gaussian latitudes
    real lon_ncr              ! [dgr] Longitude increment
    real wgt_Gss(lat_nbr) ! Gaussian weights
    real lat_ml(lat_nbr)

    ! Main Code
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Entering map_grd_mk()"
    
    ! Initialize scalar variables
    lon_ncr=360.0/lon_nbr     ! [dgr]
    pi=4.0*atan(1.0d0) ! [frc] 3
    rdn2dgr=180.0d0/pi        ! [dgr rdn-1]
    
    ! Decode latitude code
    map_grd_m=map_typ/10
    if (map_grd_m==1) then
       map_lat_grd_typ=map_lat_grd_rgl
    else if (map_grd_m==2) then
       map_lat_grd_typ=map_lat_grd_Gss
    else if (map_grd_m==3) then
       map_lat_grd_typ=map_lat_grd_GSC
    else if (map_grd_m==4) then
       map_lat_grd_typ=map_lat_grd_ease_ml
    else
       write (6,"(2a,i4)") prg_nm(1:ftn_strlen(prg_nm)), & 
            ": ERROR map_grd_mk() reports map_grd_m = ",map_grd_m
       stop
    endif                     ! endif
    
    ! Decode longitude code
    map_grd_n=mod(map_typ,10)
    if (map_grd_n==1) then
       map_lon_srt_typ=map_lon_srt_180
       map_lon_ctr_typ=map_lon_ctr_wst
    else if (map_grd_n==2) then
       map_lon_srt_typ=map_lon_srt_Grn
       map_lon_ctr_typ=map_lon_ctr_wst
    else if (map_grd_n==3) then
       map_lon_srt_typ=map_lon_srt_Grn
       map_lon_ctr_typ=map_lon_ctr_ctr
    else if (map_grd_n==4) then
       map_lon_srt_typ=map_lon_srt_180
       map_lon_ctr_typ=map_lon_ctr_ctr
    else
       write (6,"(2a,i4)") prg_nm(1:ftn_strlen(prg_nm)), & 
            ": ERROR map_grd_mk() reports map_grd_n = ",map_grd_n
       stop
    endif                     ! endif
    
    ! All maps begin at longitude 0.0 or 180.0 degrees
    if (map_lon_srt_typ==map_lon_srt_Grn) then
       lon_grd(1)=0.0
    else if (map_lon_srt_typ==map_lon_srt_180) then
       lon_grd(1)=-180.0
    else
       write (6,"(2a,i4)") prg_nm(1:ftn_strlen(prg_nm)), & 
            ": ERROR map_grd_mk() reports map_lon_srt_typ = ",map_lon_srt_typ
       stop
    endif
    ! Whether 0.0 or 180.0 refers to cell center or Western edge is specified with map_lon_ctr_typ argument
    if (map_lon_ctr_typ==map_lon_ctr_ctr) then
       lon_grd(1)=lon_grd(1)-(lon_ncr/2.0)
    else if (map_lon_ctr_typ /= map_lon_ctr_wst) then
       write (6,"(2a,i4)") prg_nm(1:ftn_strlen(prg_nm)), & 
            ": ERROR map_grd_mk() reports map_lon_ctr_typ = ",map_lon_ctr_typ
       stop
    endif
    
    ! Create longitude array
    do lon_idx=2,lon_nbr
       lon_grd(lon_idx)=lon_grd(1)+(lon_idx-1)*lon_ncr
    end do                     ! end loop over lon
    ! This ensures rounding errors do not produce unphysical grid
    lon_grd(lon_nbr+1)=lon_grd(1)+360.0
    
    ! All maps begin at southernmost latitude
    lat_grd(1)=-90.0
    if (map_lat_grd_typ==map_lat_grd_rgl) then
       lat_ncr=180.0/lat_nbr ! [dgr] Latitude increment
       do lat_idx=2,lat_nbr
          lat_grd(lat_idx)=lat_grd(1)+(lat_idx-1)*lat_ncr
       end do                  ! end loop over lat
    else if (map_lat_grd_typ==map_lat_grd_Gss) then
       call lat_Gss_mk(lat_nbr, & ! I
            lat_sin,wgt_Gss) ! O
       lat=-asin(lat_sin)*rdn2dgr
       do lat_idx=2,lat_nbr
          lat_grd(lat_idx)=0.5*(lat(lat_idx-1)+lat(lat_idx))
       end do                  ! end loop over lat
    else if (map_lat_grd_typ==map_lat_grd_GSC) then
       lat_ncr=180.0/(lat_nbr-1) ! [dgr] Latitude increment
       lat_grd(2)=lat_grd(1)+0.5*lat_ncr
       do lat_idx=3,lat_nbr ! NB: Weird GEOS grid is uneven at Poles
          lat_grd(lat_idx)=lat_grd(2)+(lat_idx-2)*lat_ncr ! NB: lat_idx-2
       end do                  ! end loop over lat
    else if (map_lat_grd_typ==map_lat_grd_ease_ml) then
       call lat_ease_ml_mk(lat_nbr,lat_ml)    
       do lat_idx=2,lat_nbr
          lat_grd(lat_idx)=lat_ml(lat_idx) ! Low boundary from S->N
       enddo
    else
       write (6,"(2a,i4)") prg_nm(1:ftn_strlen(prg_nm)), & 
            ": ERROR map_grd_mk() reports map_lat_grd_typ = ",map_lat_grd_typ
       stop
    endif
    ! This ensures rounding errors do not produce unphysical grid
    lat_grd(lat_nbr+1)=90.0

    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Exiting map_grd_mk()"
    return 
  end subroutine map_grd_mk                       ! end map_grd_mk()
  
  subroutine map_lat_wgt_mk(lat_nbr,map_typ, & ! I
       lat_wgt) ! O
    ! Purpose: Construct and return latitude weights (integrating factors) of grid
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use map_cst ! [mdl] Constants used in map routines
    use sng_mdl,only:ftn_strlen ! [mdl] String manipulation
    implicit none
    ! Parameters
    ! Commons
    ! Input
    integer,intent(in)::lat_nbr ! [nbr] Number of latitudes
    integer,intent(in)::map_typ ! [enm] Map type in LSM "mn" format
    ! Output
    real,intent(out)::lat_wgt(lat_nbr) ! [frc] Latitude weights
    ! Local
    real lat_grd(lat_nbr+1) ! [dgr] Interface latitudes
    real lat_ml(lat_nbr)
    real(selected_real_kind(p=12))::pi       ! [frc] 3
    real(selected_real_kind(p=12))::rdn2dgr  ! Convert radians to degrees
    real(selected_real_kind(p=12))::dgr2rdn  ! Convert degrees to radians
    integer lat_idx           ! [idx] Counting index
    integer map_grd_m         ! [enm] Latitude component of map type
    integer map_lat_grd_typ   ! Latitude grid type
    real lat(lat_nbr) ! Centered latitudes
    real lat_ncr              ! [dgr] Latitude increment
    real lat_sin(lat_nbr) ! Sine of Gaussian latitudes
    ! Main Code
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Entering map_lat_wgt_mk()"
    
    ! Initialize scalar variables
    pi=4.0*atan(1.0d0) ! [frc] 3
    rdn2dgr=180.0d0/pi        ! [dgr rdn-1]
    dgr2rdn=pi/180.0d0        ! [rdn dgr-1]
    
    ! Decode latitude code
    map_grd_m=map_typ/10
    if (map_grd_m==1) then
       map_lat_grd_typ=map_lat_grd_rgl
    else if (map_grd_m==2) then
       map_lat_grd_typ=map_lat_grd_Gss
    else if (map_grd_m==3) then
       map_lat_grd_typ=map_lat_grd_GSC
    else if (map_grd_m==4) then
       map_lat_grd_typ=map_lat_grd_ease_ml
    else
       write (6,"(2a,i4)") prg_nm(1:ftn_strlen(prg_nm)), & 
            ": ERROR map_lat_wgt_mk() reports map_grd_m = ",map_grd_m
       stop
    endif                     ! endif
    
    ! Create latitude array
    ! All maps begin at southernmost latitude
    lat_grd(1)=-90.0          ! [dgr]
    if (map_lat_grd_typ==map_lat_grd_rgl) then
       lat_ncr=180.0/lat_nbr     ! [dgr] Latitude increment
       do lat_idx=2,lat_nbr
          lat_grd(lat_idx)=lat_grd(1)+(lat_idx-1)*lat_ncr ! [dgr]
       end do                  ! end loop over lat
       ! This ensures rounding errors do not produce unphysical grid
       lat_grd(lat_nbr+1)=90.0 ! [dgr]
       do lat_idx=1,lat_nbr
          lat(lat_idx)=0.5*(lat_grd(lat_idx)+lat_grd(lat_idx+1)) ! [dgr]
          lat_wgt(lat_idx)=cos(dgr2rdn*lat(lat_idx))
       end do                 ! end loop over lat
    else if (map_lat_grd_typ==map_lat_grd_Gss) then
       call lat_Gss_mk(lat_nbr, & ! I
            lat_sin,lat_wgt) ! O
    else if (map_lat_grd_typ==map_lat_grd_GSC) then
       lat_ncr=180.0/(lat_nbr-1) ! [dgr] Latitude increment
       lat_grd(2)=lat_grd(1)+0.5*lat_ncr
       do lat_idx=3,lat_nbr ! NB: Weird GEOS grid is uneven at Poles
          lat_grd(lat_idx)=lat_grd(2)+(lat_idx-2)*lat_ncr ! NB: lat_idx-2
       end do                  ! end loop over lat
       ! This ensures rounding errors do not produce unphysical grid
       lat_grd(lat_nbr+1)=90.0 ! [dgr]
       do lat_idx=1,lat_nbr
          lat(lat_idx)=0.5*(lat_grd(lat_idx)+lat_grd(lat_idx+1)) ! [dgr]
          lat_wgt(lat_idx)=cos(dgr2rdn*lat(lat_idx))
       end do                 ! end loop over lat
     else if (map_lat_grd_typ==map_lat_grd_ease_ml) then
      call lat_ease_ml_mk(lat_nbr,lat_ml)
      do lat_idx=2,lat_nbr
      lat_grd(lat_idx)=lat_ml(lat_idx)
      enddo
      do lat_idx=1,lat_nbr
          lat(lat_idx)=0.5*(lat_grd(lat_idx)+lat_grd(lat_idx+1)) ! [dgr]
          lat_wgt(lat_idx)=cos(dgr2rdn*lat(lat_idx))
      end do                 ! end loop over lat
    else
       write (6,"(2a,i4)") prg_nm(1:ftn_strlen(prg_nm)), & 
            ": ERROR map_lat_wgt_mk() reports map_lat_grd_typ = ",map_lat_grd_typ
       stop
    endif
    
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Exiting map_lat_wgt_mk()"
    return 
  end subroutine map_lat_wgt_mk                       ! end map_lat_wgt_mk()
  
  subroutine lat_Gss_mk(lat_nbr, & ! I
       lat_sin,wgt_Gss) ! O
    ! Purpose: Compute and return sine of Gaussian latitudes and their weights
    ! Source: CCM /fs/cgd/csm/models/atm/ccm3.5.8/src/ccmlsm_share/gauaw.F
    ! Converted all real*16 statements to double precision (real*8)
    ! Calculate sine of latitudes lat_sin(lat_nbr) and weights wgt_Gss(lat_nbr) for Gaussian quadrature. 
    ! Algorithm is described in Davis and Rabinowitz, Journal of Research of the NBS, V 56, Jan 1956
    ! Zeros of Bessel function j0, obtained from bsl_zro_get(), are first guess for abscissae
    ! Original version:  CCM1
    ! Standardized:      L. Bath, Jun 1992
    ! L. Buja, Feb 1996
    ! Reviewed:          D. Williamson, J. Hack, Aug 1992
    ! D. Williamson, J. Hack, Feb 1996
    ! Modified 1997/01/23 by Jim Rosinski to use real*16 arithmetic in order to 
    ! achieve (nearly) identical weights and latitudes on all machines.
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use sng_mdl,only:ftn_strlen ! [mdl] String manipulation
    implicit none
    ! Parameters
    real(selected_real_kind(p=12)),parameter::eps_rlt=1.0d-15 ! Convergence criterion (NB: Threshold was 1.0d-27 in real*16)
    real(selected_real_kind(p=12)),parameter::one=1.0d0 ! High precision
    integer,parameter::itr_nbr_max=20 ! [nbr] Maximum number of iterations
    ! Commons
    ! Input
    integer,intent(in)::lat_nbr           ! [nbr] Number of latitudes from pole to pole 
    ! Output
    real,intent(out)::lat_sin(lat_nbr) ! Sine of Gaussian latitudes
    real,intent(out)::wgt_Gss(lat_nbr) ! Gaussian weights
    ! Local
    real(selected_real_kind(p=12))::c        ! Constant combination
    real(selected_real_kind(p=12))::lat_idx_dbl ! Latitude index, double precision
    real(selected_real_kind(p=12))::lat_nnr_idx_dbl ! Inner latitude index, double precision
    real(selected_real_kind(p=12))::lat_nbr_dbl ! [nbr] Number of latitudes, double precision
    real(selected_real_kind(p=12))::lat_sin_dbl(lat_nbr) ! Sine of Gaussian latitudes double precision
    real(selected_real_kind(p=12))::pi       ! [frc] 3
    real(selected_real_kind(p=12))::pk       ! Polynomial
    real(selected_real_kind(p=12))::pkm1     ! Polynomial
    real(selected_real_kind(p=12))::pkm2     ! Polynomial
    real(selected_real_kind(p=12))::pkmrk    ! Polynomial
    real(selected_real_kind(p=12))::sp       ! Current iteration latitude increment
    real(selected_real_kind(p=12))::wgt_Gss_dbl(lat_nbr) ! Gaussian weights double precision
    real(selected_real_kind(p=12))::xz       ! Abscissa estimate
    integer itr_cnt           ! Iteration counter
    integer lat_hlf_idx       ! [idx] Counting index (half latitude index)
    integer lat_idx           ! [idx] Counting index (latitude)
    integer lat_sym_idx       ! [idx] Counting index (symmetric latitude)
    integer lat_nnr_idx       ! [idx] Counting index (inner latitude loop)
    integer lat_nbr_rcp2      ! lat_nbr/2 (number of latitudes in hemisphere)
    ! Main Code
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Entering lat_Gss_mk()"
    
    ! Initialize scalars
    pi=4.0*atan(one) ! [frc] 3
    
    ! The value eps_rlt, used for convergence tests in the iterations, can be changed.  
    ! Use Newton iteration to find the abscissas. 
    c=(1.d0-(2.d0/pi)**2)*0.25
    lat_nbr_dbl=lat_nbr
    lat_nbr_rcp2=lat_nbr/2
    call bsl_zro_get(lat_nbr_rcp2, & ! I
         lat_sin_dbl) ! O
    do lat_idx=1,lat_nbr_rcp2
       xz=cos(lat_sin_dbl(lat_idx)/sqrt((lat_nbr_dbl+0.5)**2+c))
       ! First approximation to xz
       itr_cnt=0
73     pkm2=1.d0
       pkm1=xz
       itr_cnt=itr_cnt+1
       if (itr_cnt > itr_nbr_max) then
          ! Error exit
          write (6,"(2a,i3,a,i3)") prg_nm(1:ftn_strlen(prg_nm)), & 
               ": ERROR lat_Gss_mk() reports no convergence in ", &
               itr_nbr_max," iterations for lat_idx = ",lat_idx
          stop
       end if
       ! Compute Legendre polynomial
       do lat_nnr_idx=2,lat_nbr
          lat_nnr_idx_dbl=lat_nnr_idx
          pk=((2.0*lat_nnr_idx_dbl-1.0)*xz*pkm1-(lat_nnr_idx_dbl-1.0)*pkm2)/lat_nnr_idx_dbl
          pkm2=pkm1
          pkm1=pk
       end do ! end inner loop over lat
       pkm1=pkm2
       pkmrk=(lat_nbr_dbl*(pkm1-xz*pk))/(1.0-xz**2)
       sp=pk/pkmrk
       xz=xz-sp
       if (abs(sp) > eps_rlt) go to 73
       lat_sin_dbl(lat_idx)=xz
       wgt_Gss_dbl(lat_idx)=(2.0*(1.0-xz**2))/(lat_nbr_dbl*pkm1)**2
    end do ! end outer loop over lat
    if (lat_nbr /= lat_nbr_rcp2*2) then
       ! When lat_nbr is odd, compute weight at the Equator
       lat_sin_dbl(lat_nbr_rcp2+1)=0.0
       pk=2.0/lat_nbr_dbl**2
       do lat_idx=2,lat_nbr,2
          lat_idx_dbl=lat_idx
          pk=pk*lat_idx_dbl**2/(lat_idx_dbl-1.0)**2
       end do ! end loop over lat
       wgt_Gss_dbl(lat_nbr_rcp2+1)=pk
    end if                    ! endif lat_nbr is odd
    
    ! Complete sets of abscissas and weights, using symmetry properties
    ! Also note truncation from double precision to real
    do lat_idx=1,lat_nbr_rcp2
       lat_sym_idx=lat_nbr-lat_idx+1
       lat_sin(lat_idx)=lat_sin_dbl(lat_idx)
       lat_sin(lat_sym_idx)=-lat_sin_dbl(lat_idx)
       wgt_Gss(lat_idx)=wgt_Gss_dbl(lat_idx)
       wgt_Gss(lat_sym_idx)=wgt_Gss_dbl(lat_idx)
    end do ! end loop over lat
    
    if (dbg_lvl==dbg_old) then
       write (6,"(a,i3)") "lat_Gss_mk(): lat_nbr = ",lat_nbr
       write (6,"(5(a11,1x))") "idx","asin","ngl_rad","ngl_dgr","gw"
       do lat_idx=1,lat_nbr
          write (6,"(i11,1x,4(f11.7,1x))") lat_idx,lat_sin(lat_idx),asin(lat_sin(lat_idx)), &
               180.0*asin(lat_sin(lat_idx))/pi,wgt_Gss(lat_idx)
       end do ! end loop over zeros outside table
    endif                     ! endif dbg
    
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Exiting lat_Gss_mk()"
    return 
  end subroutine lat_Gss_mk                       ! end lat_Gss_mk()
  
  subroutine bsl_zro_get(bsl_zro_nbr, & ! I
       bsl_zro) ! O
    ! Purpose: Return the zeros of the Bessel function
    ! Source: CCM code /fs/cgd/csm/models/atm/ccm3.5.8/src/ccmlsm_share/bsslzr.F
    ! Return bsl_zro_nbr zeros (or if bsl_zro_nbr > 50, approximate zeros), of the Bessel function j0,in the array bes. 
    ! The first 50 zeros will be given exactly, and the remaining zeros are computed by extrapolation, and therefore are not exact.
    ! Original version:  CCM1
    ! Standardized:      J. Rosinski, June 1992
    ! Reviewed:          J. Hack, D. Williamson, August 1992
    ! Reviewed:          J. Hack, D. Williamson, April 1996
    ! Modified 1/23/97 by Jim Rosinski to use double precision arithmetic
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    implicit none
    ! Parameters
    real(selected_real_kind(p=12)),parameter::one=1.0d0 ! High precision to ensure reproducibility
    integer,parameter::bsl_zro_tbl_nbr_max=50
    ! Commons
    ! Input
    integer,intent(in)::bsl_zro_nbr       ! [nbr] Number of zeros to return
    ! Output
    real(selected_real_kind(p=12)),intent(out)::bsl_zro(bsl_zro_nbr) ! Zeros of Bessel function j0
    ! Local
    real(selected_real_kind(p=12))::pi       ! [frc] 3
    real(selected_real_kind(p=12))::bsl_zro_tbl(bsl_zro_tbl_nbr_max) ! Table of first bsl_zro_tbl_nbr_max zeros
    real val_flt              ! Scalar for debugging
    integer bsl_idx           ! [idx] Counting index
    save bsl_zro_tbl          ! Ensures re-entrancy
    data bsl_zro_tbl           / 2.4048255577,   5.5200781103, &
         8.6537279129,  11.7915344391,  14.9309177086,  18.0710639679, &
         21.2116366299,  24.3524715308,  27.4934791320,  30.6346064684, &
         33.7758202136,  36.9170983537,  40.0584257646,  43.1997917132, &
         46.3411883717,  49.4826098974,  52.6240518411,  55.7655107550, &
         58.9069839261,  62.0484691902,  65.1899648002,  68.3314693299, &
         71.4729816036,  74.6145006437,  77.7560256304,  80.8975558711, &
         84.0390907769,  87.1806298436,  90.3221726372,  93.4637187819, &
         96.6052679510,  99.7468198587, 102.8883742542, 106.0299309165, &
         109.1714896498, 112.3130502805, 115.4546126537, 118.5961766309, &
         121.7377420880, 124.8793089132, 128.0208770059, 131.1624462752, &
         134.3040166383, 137.4455880203, 140.5871603528, 143.7287335737, &
         146.8703076258, 150.0118824570, 153.1534580192, 156.2950342685/
    
    ! Main Code
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Entering bsl_zro_get()"
    
    ! Initialize scalars
    pi=4.0*atan(one) ! [frc] 3
    
    do bsl_idx=1,bsl_zro_nbr
       if (bsl_idx <= bsl_zro_tbl_nbr_max) bsl_zro(bsl_idx)=bsl_zro_tbl(bsl_idx)
    end do ! end loop over zeros in table
    if (bsl_zro_nbr > bsl_zro_tbl_nbr_max) then
       do bsl_idx=bsl_zro_tbl_nbr_max+1,bsl_zro_nbr
          bsl_zro(bsl_idx)=bsl_zro(bsl_idx-1)+pi
       end do ! end loop over zeros outside table
    end if                    ! endif more zeros than table holds
    
    if (dbg_lvl==dbg_old) then
       write (6,"(a,i3)") "bsl_zro_get(): bsl_zro_nbr = ",bsl_zro_nbr
       write (6,"(4(a10,1x))") "idx","bsl_zro"
       do bsl_idx=1,bsl_zro_nbr
          write (6,"(i10,1x,3(f10.7,1x))") bsl_idx,bsl_zro(bsl_idx)
       end do ! end loop over zeros outside table
    endif                     ! endif dbg
    
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Exiting bsl_zro_get()"
    return
  end subroutine bsl_zro_get                       ! end bsl_zro_get()
  
  subroutine rnk_vec_LSM(dat,dat_nbr,mss_val, & ! I
       rnk_1st_idx,rnk_2nd_idx) ! O
    ! Purpose: Find the first and second largest elements of a zero-based array
    ! NB: Input array indexing is zero-based for use with LSM surface type arrays
    ! Input data are assumed to be physically meaningful, e.g., weights
    ! An error is flagged if all input data are zero
    ! Missing values are handled correctly
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    implicit none
    ! Parameters
    ! Commons
    ! Input
    integer,intent(in)::dat_nbr           ! [nbr] Number of array elements (minus one)
    integer,intent(in)::mss_val           ! Code for missing/invalid data
    real,intent(in)::dat(0:dat_nbr) ! Array to rank
    ! Output
    integer,intent(out)::rnk_1st_idx       ! Index of largest member of array
    integer,intent(out)::rnk_2nd_idx       ! Index of second largest member of array
    ! Local
    integer dat_idx           ! [idx] Counting index
    
    ! Main Code
    ! if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Entering rnk_vec_LSM()"
    
    rnk_1st_idx=mss_val
    do dat_idx=0,dat_nbr      ! NB: zero-based
       if (dat(dat_idx) > 0.0) then
          if (rnk_1st_idx==mss_val) then
             rnk_1st_idx=dat_idx ! First time array is > 0.0
          else
             if (dat(dat_idx) > dat(rnk_1st_idx)) rnk_1st_idx=dat_idx
          endif               ! endif
       endif                  ! endif
    end do ! end loop over dat
    if (rnk_1st_idx==mss_val) then
       write (6,"(a,i5)") "rnk_vec_LSM(): ERROR rnk_1st_idx = ",rnk_1st_idx
       do dat_idx=0,dat_nbr   ! NB: zero-based
          write (6,"(a,i2,a,f9.6)") "dat(",dat_idx,") = ",dat(dat_idx)
       end do ! end loop over dat
       ! NB: Allow function to return to calling routine for further diagnostics
    endif                     ! endif err
    
    rnk_2nd_idx=mss_val
    do dat_idx=0,dat_nbr      ! NB: zero-based
       if (dat(dat_idx) > 0.0.and.dat_idx /= rnk_1st_idx) then
          if (rnk_2nd_idx==mss_val) then
             rnk_2nd_idx=dat_idx ! First time array is > 0.0
          else
             if (dat(dat_idx) > dat(rnk_2nd_idx)) rnk_2nd_idx=dat_idx
          endif               ! endif
       endif                  ! endif
    end do ! end loop over dat
    ! NB: rnk_2nd_idx is undefined when there is only one valid member in list
    
    ! if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Exiting rnk_vec_LSM()"
    return 
  end subroutine rnk_vec_LSM                       ! end rnk_vec_LSM()
  
  subroutine rnk_vec_gnr(dat,dat_nbr,idx_not_rnk,mss_val,rnk_nbr, & ! I
       rnk,rnk_vld_nbr) ! O
    ! Purpose: Find rnk_nbr largest elements of a zero-based array
    ! NB: Input array indexing is zero-based for use with LSM/CLM surface type arrays
    ! NB: Negative values are considered invalid
    ! Input data are assumed to be physically meaningful, e.g., weights
    ! An error is flagged if all input data are zero
    ! Missing values are almost handled correctly
    ! fxm: pass mss_flg to determine whether to check for missing values
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    implicit none
    ! Parameters
    ! Commons
    ! Input
    integer,intent(in)::dat_nbr ! [nbr] Number of array elements (minus one)
    integer,intent(in)::idx_not_rnk ! [idx] Index denoting missing/invalid data
    real,intent(in)::mss_val ! [frc] Missing value
    real,intent(in)::dat(0:dat_nbr) ! [frc] Array to rank
    integer,intent(in)::rnk_nbr ! [nbr] Number of elements to rank
    ! Output
    integer,intent(out)::rnk(rnk_nbr) ! [idx] Ranked list input data array (largest to smallest)
    integer,intent(out)::rnk_vld_nbr ! [nbr] Number of valid, ranked entries returned
    ! Local
    integer dat_idx ! [idx] Counting index for dat
    integer rnk_idx ! [idx] Counting index for rnk
    integer rnk_idx_idx ! [idx] Counting index for rnk_idx
    logical flg_not_rnk ! [flg] Value not yet ranked
    
    ! Main Code
    ! if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Entering rnk_vec_gnr()"
    
    ! Initialize to unusable index
    rnk(:)=idx_not_rnk
    do dat_idx=0,dat_nbr ! NB: zero-based
       if (dat(dat_idx) /= mss_val) then
          if (dat(dat_idx) > 0.0) then
             if (rnk(1)==idx_not_rnk) then
                ! First time value is > 0.0
                rnk(1)=dat_idx !
             else
                if (dat(dat_idx) > dat(rnk(1))) rnk(1)=dat_idx
             endif ! endif first time
          endif ! endif > 0.0
       endif ! endif mss_val
    end do ! end loop over dat
    if (rnk(1)==idx_not_rnk) then
       write (6,"(a,i5)") "rnk_vec_gnr(): ERROR rnk(1) = ",rnk(1)
       do dat_idx=0,dat_nbr   ! NB: zero-based
          write (6,"(a,i2,a,f9.6)") "dat(",dat_idx,") = ",dat(dat_idx)
       end do ! end loop over dat
       ! NB: Allow function to return to calling routine for further diagnostics
    endif                     ! endif err
    
    ! Find remaining rnk_nbr-1 largest values
    do rnk_idx=2,rnk_nbr
       do dat_idx=0,dat_nbr ! NB: zero-based
          if (dat(dat_idx) /= mss_val) then
             flg_not_rnk=.true. ! [flg] Value not yet ranked
             do rnk_idx_idx=1,rnk_idx-1
                if (dat_idx == rnk(rnk_idx_idx)) flg_not_rnk=.false. ! [flg] Value not yet ranked
             end do ! end loop over rnk
             if (flg_not_rnk.and.dat(dat_idx) > 0.0) then
                if (rnk(rnk_idx)==idx_not_rnk) then
                   ! First time array is > 0.0
                   rnk(rnk_idx)=dat_idx ! 
                else
                   if (dat(dat_idx) > dat(rnk(rnk_idx))) rnk(rnk_idx)=dat_idx
                endif ! endif
             endif ! endif
          endif ! endif mss_val
       end do ! end loop over dat
    end do ! end loop over rnk
    ! NB: rnk(i+1:rnk_nbr) is idx_not_rnk when there are only i valid members in list
    
    ! How many valid, ranked entries were returned?
    rnk_vld_nbr=0 ! [nbr] Number of valid, ranked entries returned
    do rnk_idx=1,rnk_nbr
       if (rnk(rnk_idx) == idx_not_rnk) exit
       rnk_vld_nbr=rnk_vld_nbr+1 ! [nbr] Number of valid, ranked entries returned
    end do ! end loop over rnk
    
    ! Sanity check
    if (rnk_vld_nbr == 0) stop "rnk_vld_nbr == 0 in rnk_vec_gnr()"
    
    ! if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Exiting rnk_vec_gnr()"
    return 
  end subroutine rnk_vec_gnr                       ! end rnk_vec_gnr()
  
  subroutine map_edge_mk(lat_nbr,lon_nbr,lon_nbr_1d, & ! I
       lat_ctr_2d,lon_ctr_2d,edg_grd, & ! I 
       lat_grd,lon_grd,lon_grd_2d) ! O
    ! Purpose: Convert gridpoint centers and domain boundaries (edges),
    ! to interfaces (edges) of latitude-longitude grid.
    ! For south --> north grids, interfaces are south edges of gridcells
    ! For north --> south grids, interfaces are north edges of gridcells
    ! For east --> west grids, interfaces are west edges of gridcells
    ! For west --> east grids, interfaces are east edges of gridcells
    ! Output values:
    ! Irregular grids and reduced grids require lon_grd_2d for interface locations,
    ! lon_nbr_1d for number of valid longitudes, lon_nbr for greatest dimension size
    ! lon_grd is invalid for irregular grids
    ! Regular grids are fully described by lon_grd and lon_nbr
    ! lon_grd_2d and lon_nbr_1d are valid, but superfluous, for regular grids
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use map_cst ! [mdl] Constants used in map routines
    use sng_mdl,only:ftn_strlen ! [mdl] String manipulation
    implicit none
    ! Parameters
    character(*),parameter::sbr_nm="map_edge_mk()" ! [sng] Subroutine name
    real,parameter::mss_val=-999.0 ! [dgr] Missing value for irregular grids
    ! Commons
    ! Input
    integer,intent(in)::lon_nbr ! [nbr] Maximum number of longitudes
    integer,intent(in)::lat_nbr ! [nbr] Number of latitudes
    integer,intent(in)::lon_nbr_1d(lat_nbr) ! [nbr] Longitudes per latitude
    real,intent(in)::lon_ctr_2d(lon_nbr,lat_nbr) ! [nbr] Longitude at gridcell center
    real,intent(in)::lat_ctr_2d(lon_nbr,lat_nbr) ! [nbr] Latitude at gridcell center
    real,intent(in)::edg_grd(4) ! [dgr] Grid edges (north, east, south, west)
    ! Output
    real,intent(out)::lat_grd(lat_nbr+1) ! [dgr] Interface latitudes
    real,intent(out)::lon_grd(lon_nbr+1) ! [dgr] Interface longitudes
    real,intent(out)::lon_grd_2D(lon_nbr+1,lat_nbr) ! [dgr] Interface longitudes
    ! Local
    integer lat_idx ! [idx] Counting index
    integer lon_idx ! [idx] Counting index
    real edg_nrt ! [dgr] North edge of grid
    real edg_est ! [dgr] East edge of grid
    real edg_sth ! [dgr] South edge of grid
    real edg_wst ! [dgr] West edge of grid
    real lon_ctr(lon_nbr) ! [dgr] Longitude center
    real lon_dlt(lon_nbr) ! [dgr] Longitude span
    ! Main Code
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Entering "//sbr_nm
    
    ! Sanity check
    if (lat_nbr <= 0) stop "lat_nbr <= 0 in map_edge_mk()"
    if (lon_nbr <= 0) stop "lon_nbr <= 0 in map_edge_mk()"
    
    if (dbg_lvl == dbg_io) then
       write (6,"(a)") "Data input to "//sbr_nm//":"
       write (6,"(a,i4)") "lon_nbr = ",lon_nbr
       write (6,"(a,i4)") "lat_nbr = ",lat_nbr
       write (6,"(a,4(f12.6,1x))") "Edge grid N,E,S,W: ",edg_grd(1),edg_grd(2),edg_grd(3),edg_grd(4)
       write (6,"(a,i4)") "Longitudes per latitude:"
       do lat_idx=1,lat_nbr
          write (6,"(i4,a1,1x)",advance="no") lon_nbr_1d(lat_idx),","
       end do ! end loop over lat
       write (6,"()",advance="yes") ! Carriage return after last number
    endif ! endif dbg
    
    ! Initialize scalars
    edg_nrt=edg_grd(1) ! [dgr] North edge of grid
    edg_est=edg_grd(2) ! [dgr] East edge of grid
    edg_sth=edg_grd(3) ! [dgr] South edge of grid
    edg_wst=edg_grd(4) ! [dgr] West edge of grid
    
    ! If grid is degenerate (has one latitude)...
    if (lat_nbr == 1) then
       ! ...then cell boundaries are simply north and south edges
       lat_grd(1)=edg_sth 
       lat_grd(lat_nbr+1)=edg_nrt 
    else ! ...otherwise the grid has more than one latitude...
       ! ...and south to north grids...
       if (lat_ctr_2d(1,2) > lat_ctr_2d(1,1)) then
          ! ...have south edge first, north edge last...
          lat_grd(1)=edg_sth 
          lat_grd(lat_nbr+1)=edg_nrt 
       else ! ...while north to south grids...
          ! ...have north edge first, south edge last...
          lat_grd(1)=edg_nrt 
          lat_grd(lat_nbr+1)=edg_sth 
       end if ! endif north to south grid
    end if ! endif multiple latitudes
    
    ! Intermediate interfaces are arithmetic means of adjacent centerpoints
    do lat_idx=2,lat_nbr
       lat_grd(lat_idx)=0.5*(lat_ctr_2d(1,lat_idx-1)+lat_ctr_2d(1,lat_idx))
    end do ! end loop over lat
    
    do lat_idx=1,lat_nbr
       ! Bogus CLM method assumes regularly spaced longitudes:
       ! lon_dlt=(edg_est-edg_wst)/lon_nbr_1d(lat_idx) ! [dgr] Longitude increment
       ! Western edge of each latitude band must be identical
       lon_grd_2d(1,lat_idx)=edg_wst
       ! Interior points
       do lon_idx=2,lon_nbr_1d(lat_idx)
          ! Bogus CLM method assumes regularly spaced longitudes:
          ! lon_grd_2d(lon_idx,lat_idx)=lon_grd_2d(1,lat_idx)+(lon_idx-1)*lon_dlt
          ! Interior interfaces (edges) are arithmetic means of adjacent centerpoints
          lon_grd_2d(lon_idx,lat_idx)=0.5* &
               (lon_ctr_2d(lon_idx-1,lat_idx)+lon_ctr_2d(lon_idx,lat_idx))
       end do ! end loop over lon
       ! Eastern edge of each latitude band must be identical
       lon_grd_2d(lon_nbr_1d(lat_idx)+1,lat_idx)=edg_est
       ! Irregular ("reduced") grids have undefined points "beyond" eastern edge
       do lon_idx=lon_nbr_1d(lat_idx)+2,lon_nbr
          lon_grd_2d(lon_idx,lat_idx)=mss_val
       end do ! end loop over lon
    end do ! end loop over lat
    ! Arbitrarily choose first latitude as source of 1D longitude grid
    ! User is responsible for knowing when lon_grd is invalid (all irregular grids)
    lon_grd(:)=lon_grd_2d(:,1)
    
    if (dbg_lvl == dbg_io) then
       write (6,"(a)") "Longitude grid properties from "//sbr_nm//":"
       write (6,"(a)") "Cell edge locations, center, and spacing:"
       write (6,"(a)") "(idx)   Lon Wst (idx)   Lon Est   lon_ctr   lon_dlt"
       do lon_idx=1,lon_nbr
          lon_dlt(lon_idx)=lon_grd(lon_idx+1)-lon_grd(lon_idx) ! [dgr] Longitude span
          lon_ctr(lon_idx)=0.5*(lon_grd(lon_idx)+lon_grd(lon_idx+1)) ! [dgr] Longitude center
          write (6,"(2(a1,i3,a1,1x,f9.4,1x),1(f9.4,1x,f8.5,1x))") &
               "(",lon_idx,")",lon_grd(lon_idx),"(",lon_idx+1,")",lon_grd(lon_idx+1), &
               lon_ctr(lon_idx),lon_dlt(lon_idx)
          
       end do ! end loop over lon
    endif ! endif dbg
    
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Exiting "//sbr_nm
    return 
  end subroutine map_edge_mk ! end map_edge_mk()
  
  subroutine map_typ_sng_get( & ! [sbr] Return map type description
       map_typ, & ! O [enm] Map type in LSM "mn" format
       map_typ_sng) ! O [sng] Map type description
    ! Purpose: Return description corresponding to enumerated map_typ
    ! Usage:
    ! call map_typ_sng_get(map_typ,map_typ_sng) ! [sng] Map type description
    use map_cst ! [mdl] Constants used in map routines
    use sng_mdl,only:ftn_strcpy ! [mdl] String manipulation
    implicit none
    ! Input
    integer,intent(in)::map_typ ! O [enm] Map type in LSM "mn" format
    ! Output
    character(len=*),intent(out)::map_typ_sng ! O [sng] Map type description
    ! Local
    ! Main code
    select case (map_typ)
    case (map_lat_rgl_lon_180_wst)
       call ftn_strcpy(map_typ_sng, &
            "Latitudes are regular, Date line at West edge of first longitude bin")
    case (map_lat_rgl_lon_Grn_wst)
       call ftn_strcpy(map_typ_sng, &
            "Latitudes are regular, Greenwich at West edge of first longitude bin")
    case (map_lat_rgl_lon_Grn_ctr)
       call ftn_strcpy(map_typ_sng, &
            "Latitudes are regular, Greenwich at center of first longitude bin")
    case (map_lat_rgl_lon_180_ctr)
       call ftn_strcpy(map_typ_sng, &
            "Latitudes are regular, Date line at center of first longitude bin")
    case (map_lat_Gss_lon_180_wst)
       call ftn_strcpy(map_typ_sng, &
            "Latitudes are Gaussian, Date line at West edge of first longitude bin")
    case (map_lat_Gss_lon_Grn_wst)
       call ftn_strcpy(map_typ_sng, &
            "Latitudes are Gaussian, Greenwich at West edge of first longitude bin")
    case (map_lat_Gss_lon_Grn_ctr)
       call ftn_strcpy(map_typ_sng, &
            "Latitudes are Gaussian, Greenwich at center of first longitude bin")
    case (map_lat_GSC_lon_180_ctr)
       call ftn_strcpy(map_typ_sng, &
            "Latitudes are GEOS-CHEM, Date line at center of first longitude bin")
    case default
       stop "Unknown map_typ in map_typ_sng_get()"
    end select ! end select map_typ
    return
  end subroutine map_typ_sng_get
  
  subroutine lat_ease_ml_mk (nlat,lat_ml)
    ! Convert equal area cylindrical grid coordinates to geographic coordinates (spherical earth)
    ! Input : r, s - grid column and row coordinates
    ! Output: lat, lon - geo. coords. (decimal degrees)
    implicit none
    integer,intent(in)::nlat         ! I
    real,intent(out)::lat_ml(nlat)   ! O
    !       local
    real r, s, lat
    real lat_ctr(nlat)
    integer cols, rows, scale, j, jj
    real Rg, phi, lam, rho
    real gamma, beta, epsilon, x, y
    real sinphi1, cosphi1
    real s0
    
    real,parameter::rds_earth_ease_ml=6371.228e+03 ! [m] Radius of Earth, authalic sphere, used in EASE ML grids
    real,parameter::CELL_m = 25.067525e+03 ! [m] Nominal cell size in meters
    
    !       scale factor for standard paralles at +/-30.00 degrees
    real,parameter::COS_PHI1 = .866025403
    real(selected_real_kind(p=12))::pi ! [frc] 3
    pi=4.0*atan(1.0d0) ! [frc] 3

    rows = nlat
    scale = 1
    Rg = scale * rds_earth_ease_ml/CELL_m
    s0 = (rows-1)/2.0*scale
    do j=nlat,1,-1
       jj=nlat-j+1          !index to convert lat from 90N -> 90S
       s = j-1+0.5
       y = -(s - s0)
       
       ! Allow 0.5 cell tolerance in arcsin function
       ! so that grid coordinates which are less than 0.5 cells
       ! above 90.00N or below 90.00S are given a lat of 90.00
       
       epsilon=1+0.5/Rg
       beta=y*COS_PHI1/Rg
       if (abs(beta).gt.epsilon) return
       if (beta.le.-1.0) then
          phi=-pi/2.0
       else if (beta.ge.1.0) then
          phi=pi/2.0
       else
          phi=asin(beta)
       endif
       ! lat = deg(phi)
       lat=phi*180.0/pi
       lat_ml(jj)=lat ! convert lat from 90N -> 90S
    enddo
    return       
  end subroutine lat_ease_ml_mk         ! lat_ease_ml_mk

end module map_grd ! [mdl] Map grids and regridding
