! $Id$ -*-f90-*-

! Purpose: Routines to process regridding

! Copyright (C) 2001 Charlie Zender
! Portions are copyright by their respective contributors

! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 3
! of the License, or (at your option) any later version.
! See http://www.gnu.org/copyleft/gpl.html for details.

! Author: Chao Luo
! cluo@uci.edu
! Jan. 2006

module rgr_mdl          ! Regridding module 
  implicit none
  public:: map_rgr    !mapping data in different grid
  private::map_rbn_nbr_vld  ! rebin data and count valid number
  private::read_in    !read in data in orginal grid
  private::wrt_out_rel  !write out the data in new grid
  private::wrt_out_int  !write out the valid number for regridding
contains
  
  subroutine map_rgr (dat_in,            &      !I input file
       map_typ_in,        &      !I Input map grid type
       lon_in_nbr,        &      !I Number of longitudes
       lat_in_nbr,        &      !I Number of latitudes
       map_typ_out,       &      !I Output map grid typ
       lon_out_nbr,       &      !I Number of longitudes
       lat_out_nbr,       &      !I Number of latitudes
       dat_out)                  !O output data
    use map_cst               ! [mdl] Constants used in map routines
    use map_grd               ! [mdl] Map grids and regridding
    
    !     Parameters
    !     character(len=*),intent(in)::fl_in        ! [sng] Input file
    character(80)::fl_in="rgr_in.dat"
    character(80)::fl_out_1="rgr_out.dat"
    character(80)::fl_out_2="rgr_ovr_nbr_vld.dat"
    integer,intent(in)::lat_in_nbr ! [nbr] Number of latitudes
    integer,intent(in)::lon_in_nbr ! [nbr] Number of longitudes
    integer,intent(in)::lat_out_nbr ! [nbr] Number of latitudes
    integer,intent(in)::lon_out_nbr ! [nbr] Number of longitudes
    
    integer,intent(in)::map_typ_in  ! [enm] Input map grid type
    integer,intent(in)::map_typ_out ! [enm] Output map grid type
    
    real::lat_in_grd(lat_in_nbr+1) ! [dgr] Interface latitudes
    real::lon_in_grd(lon_in_nbr+1) ! [dgr] Interface longitudes
    real::lat_out_grd(lat_out_nbr+1) ! [dgr] Interface latitudes
    real::lon_out_grd(lon_out_nbr+1) ! [dgr] Interface longitudes
    real::area_out(lon_out_nbr,lat_out_nbr) ! [m2] Area of gridcells
    integer ovr_nbr(lon_out_nbr,lat_out_nbr) ! [nbr] Number of input gridcells which overlap each output gridcell
    !     integer,intent(out):: ovr_nbr_vld(lon_out_nbr,lat_out_nbr) ! number of input valid value gridcells which overlap each output gridcell
    integer ovr_nbr_vld(lon_out_nbr,lat_out_nbr)     !number of valid value gridcells which overlap each output gridcell
    !     Allocatables
    integer,dimension(:,:,:),allocatable::ovr_lat_idx ! [idx] Map into input grid of latitude indices of overlap cells
    integer,dimension(:,:,:),allocatable::ovr_lon_idx ! [idx] Map into input grid of longitude indices of overlap cells
    real ,dimension(:,:,:),allocatable::ovr_wgt ! [frc] Weight of overlapping input gridcells onto each output gridcell
    !     Locals
    logical::mss_flg=.true.
    integer ovr_nbr_max       ! [nbr] Maximum number of input cells which overlap any output cell
    integer::rcd=0            ! [enm] Return success code
    integer::fl_in_unit=50
    integer::fl_out_unit=60
    real,parameter::mss_val_in=1.0e36   ! [frc] Missing value
    real,parameter::mss_val_out=1.0e36  ! [frc] Missing value
    integer lon_idx,lat_idx
    real  dat_in(lon_in_nbr,lat_in_nbr)
    real, intent(out)::dat_out(lon_out_nbr,lat_out_nbr)
    
    !     Create input map grid
    call map_grd_mk(lat_in_nbr,lon_in_nbr,map_typ_in, & ! I
         lat_in_grd,lon_in_grd) ! O
    
    !     Create output map grid
    call map_grd_mk(lat_out_nbr,lon_out_nbr,map_typ_out, & ! I
         lat_out_grd,lon_out_grd) ! O
    
    !     Diagnostic area on output grid
    call map_area_get(lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
         area_out)        ! O
    
    !     Determine space required by overlap arrays
    call map_ovr_nbr_max_get(lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
         lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
         ovr_nbr_max)     ! O
    
    !              write (6,'(3(a,i))') 'lon_out_nbr = ',lon_out_nbr,'lat_out_nbr = ',lat_out_nbr,'ovr_nbr_max_nbr = ',ovr_nbr_max

    !     Allocate arrays that depend on ovr_nbr_max
    allocate(ovr_lat_idx(lon_out_nbr,lat_out_nbr,ovr_nbr_max),stat=rcd) ! [idx] Map into input grid of latitude indices of overlap cells
    if(rcd /= 0) stop "allocate() failed for ovr_lat_idx"
    allocate(ovr_lon_idx(lon_out_nbr,lat_out_nbr,ovr_nbr_max),stat=rcd) ! [idx] Map into input grid of longitude indices of overlap cells
    if(rcd /= 0) stop "allocate() failed for ovr_lon_idx"
    allocate(ovr_wgt(lon_out_nbr,lat_out_nbr,ovr_nbr_max),stat=rcd) ! [frc] Weight of overlapping input gridcells onto each output gridcell
    if(rcd /= 0) stop "allocate() failed for ovr_wgt"

    !     Get overlap locations and weights for mapping input grid to output grid
    call map_ovr_wgt_drv(&
         lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
         lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
         ovr_nbr_max,area_out, & ! I
         ovr_lat_idx,ovr_lon_idx,ovr_nbr,ovr_wgt) ! O
    
    ! read in  data
    
    
    call read_in (lon_in_nbr,lat_in_nbr,dat_in,fl_in)
    
    !     Rebin land fraction from input grid to output grid
    call map_rbn_nbr_vld(dat_in,mss_flg,mss_val_in, & ! I
         lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
         lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
         ovr_nbr_max,ovr_lat_idx,ovr_lon_idx,ovr_nbr,ovr_wgt, & ! I
         ovr_nbr_vld, &  ! O
         dat_out)       ! O
    
    ! write out dat_out
    call wrt_out_rel (lon_out_nbr,lat_out_nbr,dat_out,fl_out_1)
    ! write out ovr_nbr_vld
    call wrt_out_int (lon_out_nbr,lat_out_nbr,ovr_nbr_vld,fl_out_2)
    
  end subroutine map_rgr! end program interp
  
  subroutine map_rbn_nbr_vld (dat_in,mss_flg,mss_val, & ! I
       lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
       lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
       ovr_nbr_max,ovr_lat_idx,ovr_lon_idx,ovr_nbr,ovr_wgt, & ! I
       ovr_nbr_vld, &  !O
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
    integer,intent(out):: ovr_nbr_vld(lon_out_nbr,lat_out_nbr)
    
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
    !$omp shared(ovr_nbr,ovr_lon_idx,ovr_lat_idx,ovr_wgt,ovr_nbr_vld) &
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

!cl+++      add ovr_nbr_vld
            ovr_nbr_vld(lon_out_idx,lat_out_idx)=ovr_nbr_crr-mss_cnt
!cl+++
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
  end subroutine map_rbn_nbr_vld                       ! end map_rbn()
  
  subroutine read_in (lon_in_nbr,  &    ! I
       lat_in_nbr,  &    ! I
       dat_in,      &    ! I
       fl_in)            ! I
    implicit none
    integer,intent(in)::lat_in_nbr ! [nbr] Number of latitudes
    integer,intent(in)::lon_in_nbr ! [nbr] Number of longitudes
    character(len=*), intent(in)::fl_in
    real dat_in(lon_in_nbr,lat_in_nbr)
    integer, parameter          ::unit_in = 10
    integer                     ::OpenStatus
    integer  :: lat_idx
    integer  :: lon_idx
    open (unit=unit_in, file=fl_in, status="old", iostat=OpenStatus)
    if (OpenStatus > 0) stop "ERROR: cannot open file for reading"
    read (unit_in, *) ((dat_in(lon_idx,lat_idx),lon_idx=1,lon_in_nbr),lat_idx=1,lat_in_nbr)
    close (unit_in)
  end subroutine read_in
  
  
  subroutine wrt_out_rel (lon_out_nbr,  &    ! I
       lat_out_nbr,  &    ! I
       dat_out,      &    ! I
       fl_out)            ! I
    implicit none
    integer,intent(in)::lat_out_nbr ! [nbr] Number of latitudes
    integer,intent(in)::lon_out_nbr ! [nbr] Number of longitudes
    character(len=*), intent(in)::fl_out
    real dat_out(lon_out_nbr,lat_out_nbr)
    integer, parameter          ::unit_out = 20
    integer                     ::OpenStatus
    integer  :: lat_idx
    integer  :: lon_idx
    open (unit=unit_out, file=fl_out, status="old", iostat=OpenStatus)
    if (OpenStatus > 0) stop "ERROR: cannot open file for reading"
    write (unit_out, *) ((dat_out(lon_idx,lat_idx),lon_idx=1,lon_out_nbr),lat_idx=1,lat_out_nbr)
    close (unit_out)
  end subroutine wrt_out_rel
  
  subroutine wrt_out_int (lon_out_nbr,  &    ! I
       lat_out_nbr,  &    ! I
       dat_out,      &    ! I
       fl_out)            ! I
    implicit none
    integer,intent(in)::lat_out_nbr ! [nbr] Number of latitudes
    integer,intent(in)::lon_out_nbr ! [nbr] Number of longitudes
    character(len=*), intent(in)::fl_out
    integer dat_out(lon_out_nbr,lat_out_nbr)
    integer, parameter          ::unit_out = 30
    integer                     ::OpenStatus
    integer  :: lat_idx
    integer  :: lon_idx
    open (unit=unit_out, file=fl_out, status="old", iostat=OpenStatus)
    if (OpenStatus > 0) stop "ERROR: cannot open file for reading"
    write (unit_out, *) ((dat_out(lon_idx,lat_idx),lon_idx=1,lon_out_nbr),lat_idx=1,lat_out_nbr)
    close (unit_out)
  end subroutine wrt_out_int
  
end module rgr_mdl
