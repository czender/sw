! $Id$ -*-f90-*-

! Purpose: Lakes and wetlands

! Usage: 
! use lak_wtl ! [mdl] Lakes and wetlands

module lak_wtl ! [mdl] Lakes and wetlands
  implicit none
  private ! [stt] Symbols are private unless individually qualified as public
  public::lnd_frc_dry_get ! [sbr] Determine dry land fraction
  
contains
  
  subroutine lnd_frc_dry_get(fl_lak,fl_wtl, & ! I
       lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
       lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
       ovr_nbr_max,area_out,lnd_frc_out,sfc_typ_LSM_out, & ! I
       lak_frc_out,lnd_frc_dry_out,wtl_frc_out) ! O
    ! Purpose: Retrieve lake and wetland fractions, use these to determine dry land fraction
    ! lnd_frc_dry_out is fraction of land, if any, in gridcell that is not lake or wetland
    ! Consequently, ocean points (islands and pure ocean) have lnd_frc_dry=1.0
    ! This illustrates that lnd_frc_dry_out is only physically meaninful as a multiplier of lnd_frc
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use map_grd ! [mdl] Map grids and regridding
    use sng_mdl,only:ftn_strlen,ftn_strstr ! [mdl] String manipulation
    use utl_mdl,only:mnt_ncr_chk ! [mdl] Utility functions (date_time_get,mnt_chk...)
    implicit none
    ! Parameters
    integer,parameter::fl_in_unit=73 ! Unit for reading input data
    integer,parameter::unit_dgn=74 ! Unit for writing general Map diagnoses
    ! Commons
    ! Input
    character(80),intent(in)::fl_wtl       ! [sng] Input file, wetlands
    character(80),intent(in)::fl_lak       ! [sng] Input file, lakes
    integer,intent(in)::lat_in_nbr        ! [nbr] Number of latitudes
    integer,intent(in)::lat_out_nbr       ! [nbr] Number of latitudes
    integer,intent(in)::lon_in_nbr        ! [nbr] Number of longitudes
    integer,intent(in)::lon_out_nbr       ! [nbr] Number of longitudes
    integer,intent(in)::ovr_nbr_max       ! [nbr] Maximum number of input cells which overlap any output cell
    integer,intent(in)::sfc_typ_LSM_out(lon_out_nbr,lat_out_nbr) ! [enm] Surface type 
    real,intent(in)::area_out(lon_out_nbr,lat_out_nbr) ! [m2] Area of gridcells
    real,intent(in)::lat_in_grd(lat_in_nbr+1) ! [dgr] Interface latitudes
    real,intent(in)::lat_out_grd(lat_out_nbr+1) ! [dgr] Interface latitudes
    real,intent(in)::lnd_frc_out(lon_out_nbr,lat_out_nbr) ! Land fraction of gridcell 
    real,intent(in)::lon_in_grd(lon_in_nbr+1) ! [dgr] Interface longitudes
    real,intent(in)::lon_out_grd(lon_out_nbr+1) ! [dgr] Interface longitudes
    ! Output
    real,intent(out)::lnd_frc_dry_out(lon_out_nbr,lat_out_nbr) ! Dry fraction of gridcell 
    real,intent(out)::lak_frc_out(lon_out_nbr,lat_out_nbr) ! Lake fraction of gridcell 
    real,intent(out)::wtl_frc_out(lon_out_nbr,lat_out_nbr) ! Wetland fraction of gridcell 
    ! Local
    character(80)::fl_dgn_lak   ! Diagnostic file, lakes
    character(80)::fl_dgn_wtl   ! Diagnostic file, wetlands
    integer frc_int(lon_in_nbr,lat_in_nbr) ! Integer array for input
    integer idx               ! [idx] Counting index
    integer lat_in_idx        ! [idx] Counting index for lat
    integer lat_out_idx       ! [idx] Counting index for lat
    integer lon_in_idx        ! [idx] Counting index for lon
    integer lon_out_idx       ! [idx] Counting index for lon
    integer ovr_idx           ! [idx] Counting index
    integer ovr_lat_idx(lon_out_nbr,lat_out_nbr,ovr_nbr_max) ! [idx] Map into input grid of latitude indices of overlap cells
    integer ovr_lon_idx(lon_out_nbr,lat_out_nbr,ovr_nbr_max) ! [idx] Map into input grid of longitude indices of overlap cells
    integer ovr_nbr(lon_out_nbr,lat_out_nbr) ! [nbr] Number of input gridcells which overlap each output gridcell
    integer ovr_nbr_crr       ! [nbr] Current number of overlapping gridcells
    integer rcd               ! [enm] Return success code
    logical mnt_ncr           ! [flg] Monotonic and increasing flag
    logical mss_flg           ! [flg] Check for mss_val on input data
    real area_in(lon_in_nbr,lat_in_nbr) ! [m2] Area of gridcells
    real eps_rlt              ! [frc] Relative error allowed in frc_ttl
    real frc_ttl              ! [frc] Total fraction of gridcell accounted for by all soil textures
    real lak_frc_in(lon_in_nbr,lat_in_nbr) ! Lake fraction of gridcell 
    real mss_val              ! Missing value
    real ovr_wgt(lon_out_nbr,lat_out_nbr,ovr_nbr_max) ! [frc] Weight of overlapping input gridcells onto each output gridcell
    real ovr_wgt_crr          ! [frc] Overlap weight of current gridcell
    real wtl_frc_in(lon_in_nbr,lat_in_nbr) ! Wetland fraction of gridcell 
    
    ! Main Code
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Entering lnd_frc_dry_get()"
    
    ! Check for monotonically increasing grids
    mnt_ncr=mnt_ncr_chk(lon_out_grd,lon_out_nbr+1)
    if (.not.mnt_ncr) stop "lon_out_grd not monotonically increasing in lnd_frc_dry_get()"
    mnt_ncr=mnt_ncr_chk(lat_out_grd,lat_out_nbr+1)
    if (.not.mnt_ncr) stop "lat_out_grd not monotonically increasing in lnd_frc_dry_get()"
    
    ! Initialize scalars
    eps_rlt=1.0e-5            ! [frc] Relative error allowed in frc_ttl
    mss_flg=.false.
    mss_val=-99999
    
    ! Initialize output arrays
    do lat_out_idx=1,lat_out_nbr
       do lon_out_idx=1,lon_out_nbr
          wtl_frc_out(lon_out_idx,lat_out_idx)=0.0
          lak_frc_out(lon_out_idx,lat_out_idx)=0.0
          lnd_frc_dry_out(lon_out_idx,lat_out_idx)=0.0
       end do ! end loop over lon
    end do ! end loop over lat
    
    ! Compute gridcell area 
    call map_area_get(lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
         area_in) ! O
    
    ! Cogley's data are stored as a one column formatted text file, a separate file for lakes and for wetlands
    ! The datum, i.e., lake or wetland (swamp) percentage for each (latitude,longitude) point is written alone on one row
    ! The dataset resolution is 1 x 1 degree and the data are written from the South pole to the North pole, East to West, beginning at Greenwich
    ! Read strategy is read the entire file at once, as only Fortran can
    
    ! Read lake fraction
    open (fl_in_unit,file=fl_lak,status="old",iostat=rcd)
    if (rcd /= 0) write (6,"(4a,i4)") prg_nm(1:ftn_strlen(prg_nm)), & 
         ": ERROR lnd_frc_dry_get() unable to open ",fl_lak(1:ftn_strlen(fl_lak))
    read (fl_in_unit,*,iostat=rcd) frc_int
    if (rcd /= 0) write (6,"(3a)") prg_nm(1:ftn_strlen(prg_nm)), & 
         ": ERROR lnd_frc_dry_get() reading ",fl_lak(1:ftn_strlen(fl_lak))
    close (fl_in_unit)
    write (6,"(a,1x,a)") "Read lake fraction data from",fl_lak(1:ftn_strlen(fl_lak))
    do lat_in_idx=1,lat_in_nbr
       do lon_in_idx=1,lon_in_nbr
          lak_frc_in(lon_in_idx,lat_in_idx)=frc_int(lon_in_idx,lat_in_idx)*0.01
       end do ! end loop over lon
    end do ! end loop over lat
    
    ! Read wetland fraction
    open (fl_in_unit,file=fl_wtl,status="old",iostat=rcd)
    if (rcd /= 0) write (6,"(4a,i4)") prg_nm(1:ftn_strlen(prg_nm)), & 
         ": ERROR lnd_frc_dry_get() unable to open ",fl_wtl(1:ftn_strlen(fl_wtl))
    read (fl_in_unit,*,iostat=rcd) frc_int
    if (rcd /= 0) write (6,"(3a)") prg_nm(1:ftn_strlen(prg_nm)), & 
         ": ERROR lnd_frc_dry_get() reading ",fl_wtl(1:ftn_strlen(fl_wtl))
    close (fl_in_unit)
    write (6,"(a,1x,a)") "Read wetland fraction data from",fl_wtl(1:ftn_strlen(fl_wtl))
    do lat_in_idx=1,lat_in_nbr
       do lon_in_idx=1,lon_in_nbr
          wtl_frc_in(lon_in_idx,lat_in_idx)=frc_int(lon_in_idx,lat_in_idx)*0.01
       end do ! end loop over lon
    end do ! end loop over lat
    
    ! Sanity check
    do lat_in_idx=1,lat_in_nbr
       do lon_in_idx=1,lon_in_nbr
          frc_ttl=wtl_frc_in(lon_in_idx,lat_in_idx)+ &
               lak_frc_in(lon_in_idx,lat_in_idx)
          ! 20021002: Lahey compiled code dies here due to sum > 1
          ! fxm: real solution is to test if sum is less than 1+epsilon
          !          if (frc_ttl > 1.0.and.frc_ttl <= 1+tiny(0.0)) frc_ttl=1.0
          if (frc_ttl > 1.0.and.frc_ttl <= 1.0+1.0e-5) frc_ttl=1.0
          if ( &
               (wtl_frc_in(lon_in_idx,lat_in_idx) < 0.0.or.wtl_frc_in(lon_in_idx,lat_in_idx) > 1.0).or. &
               (lak_frc_in(lon_in_idx,lat_in_idx) < 0.0.or.lak_frc_in(lon_in_idx,lat_in_idx) > 1.0).or. &
               (frc_ttl > 1.0)) then
             write (6,"(a)") "ERROR lnd_frc_dry_get() reports wetland error during input"
             write (6,"(a)") "Wetland properties:"
             write (6,"(3(a,7x))") "Wetland","Lake","Total"
             write (6,"(3(f18.12,1x))") &
                  wtl_frc_in(lon_in_idx,lat_in_idx), &
                  lak_frc_in(lon_in_idx,lat_in_idx), &
                  frc_ttl
             write (6,"(a)") "Input cell edge locations:"
             write (6,"(a)") "(idx)   Lat Sth (idx)   Lat Nrt (idx)   Lon Wst (idx)   Lon Est"
             write (6,"(4(a1,i3,a1,1x,f9.4,1x))") &
                  "(",lat_in_idx,")",lat_in_grd(lat_in_idx),"(",lat_in_idx+1,")",lat_in_grd(lat_in_idx+1), &
                  "(",lon_in_idx,")",lon_in_grd(lon_in_idx),"(",lon_in_idx+1,")",lon_in_grd(lon_in_idx+1)
             stop
          end if              ! endif err
       end do ! end loop over lon
    end do ! end loop over lat
    
    if (dbg_lvl==dbg_crr) then
       fl_dgn_lak="/tmp/zender/map/map_lak.txt"//char(0) ! Diagnostics file
       fl_dgn_wtl="/tmp/zender/map/map_wtl.txt"//char(0) ! Diagnostics file
       ! Write lake data in CSM/CCM/LSM format (one value per line)
       ! These files should be very similar to LSM input files:
       ! diff -c -w /tmp/zender/map/map_lak.txt /fs/cgd/csm/input/lnd/lnd_frc_flak.1x1 | m
       ! diff -c -w /tmp/zender/map/map_wtl.txt /fs/cgd/csm/input/lnd/lnd_frc_swmp.1x1 | m
       
       open (unit_dgn,file=fl_dgn_lak,status="unknown",iostat=rcd)
       if (rcd /= 0) write (6,"(4a,i4)") prg_nm(1:ftn_strlen(prg_nm)), & 
            ": ERROR lnd_frc_dry_get() unable to open ",fl_dgn_lak(1:ftn_strlen(fl_dgn_lak))
       do lat_in_idx=1,lat_in_nbr ! NB: Outer loop over lat
          do lon_in_idx=1,lon_in_nbr
             frc_int(lon_in_idx,lat_in_idx)=nint(lak_frc_in(lon_in_idx,lat_in_idx)*100.0)
             write (unit_dgn,*) frc_int(lon_in_idx,lat_in_idx)
          end do ! end loop over lon
       end do ! end loop over lat
       close (unit_dgn)
       write (6,"(a,1x,a)") "Wrote lake data in LSM format to",fl_dgn_lak(1:ftn_strlen(fl_dgn_lak))
       
       open (unit_dgn,file=fl_dgn_wtl,status="unknown",iostat=rcd)
       if (rcd /= 0) write (6,"(4a,i4)") prg_nm(1:ftn_strlen(prg_nm)), & 
            ": ERROR lnd_frc_dry_get() unable to open ",fl_dgn_wtl(1:ftn_strlen(fl_dgn_wtl))
       do lat_in_idx=1,lat_in_nbr ! NB: Outer loop over lat
          do lon_in_idx=1,lon_in_nbr
             frc_int(lon_in_idx,lat_in_idx)=nint(wtl_frc_in(lon_in_idx,lat_in_idx)*100.0)
             write (unit_dgn,*) frc_int(lon_in_idx,lat_in_idx)
          end do ! end loop over lon
       end do ! end loop over lat
       close (unit_dgn)
       write (6,"(a,1x,a)") "Wrote wetland data in LSM format to",fl_dgn_wtl(1:ftn_strlen(fl_dgn_wtl))
    endif                     ! endif dbg
    
    ! Get overlap locations and weights for mapping input grid to output grid
    call map_ovr_wgt_drv( &
         lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
         lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
         ovr_nbr_max,area_out, & ! I
         ovr_lat_idx,ovr_lon_idx,ovr_nbr,ovr_wgt) ! O
    
    ! Rebin wetland fraction from input grid to output grid
    call map_rbn(wtl_frc_in,mss_flg,mss_val, & ! I
         lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
         lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
         ovr_nbr_max,ovr_lat_idx,ovr_lon_idx,ovr_nbr,ovr_wgt, & ! I
         wtl_frc_out) ! O
    
    ! Rebin lake fraction from input grid to output grid
    call map_rbn(lak_frc_in,mss_flg,mss_val, & ! I
         lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
         lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
         ovr_nbr_max,ovr_lat_idx,ovr_lon_idx,ovr_nbr,ovr_wgt, & ! I
         lak_frc_out) ! O
    
    ! Lake and wetland fractions > 1.0 are possible after re-binning due to rounding errors
    ! This shows up at T170
    where (wtl_frc_out > 1.0) ! [frc]
       wtl_frc_out=1.0   ! [frc]
    end where                 ! end where wtl_frc_out < 0.0
    where (lak_frc_out > 1.0) ! [frc]
       lak_frc_out=1.0   ! [frc]
    end where                 ! end where lak_frc_out < 0.0
    
    ! Process Cogley data in accordance with LSM algorithms
    ! NB: Following code from LSM landwat.F
    do lat_out_idx=1,lat_out_nbr
       do lon_out_idx=1,lon_out_nbr
          
          ! Exclude lake and wetland areas less than 5% of cell
          if (lak_frc_out(lon_out_idx,lat_out_idx) < 0.05) lak_frc_out(lon_out_idx,lat_out_idx)=0.0
          if (wtl_frc_out(lon_out_idx,lat_out_idx) < 0.05) wtl_frc_out(lon_out_idx,lat_out_idx)=0.0
          
          ! Set oceans to zero lake and wetland fraction
          if (sfc_typ_LSM_out(lon_out_idx,lat_out_idx)==0) then
             lak_frc_out(lon_out_idx,lat_out_idx)=0.0
             wtl_frc_out(lon_out_idx,lat_out_idx)=0.0
          endif               ! not Ocean     
          
          ! Dry land is land fraction (meaning non-ocean fraction) minus lakes and wetlands
          ! fxm: Someday it might be nice to replace 1.0 in the following by lnd_frc().
          lnd_frc_dry_out(lon_out_idx,lat_out_idx)=1.0-lak_frc_out(lon_out_idx,lat_out_idx)-wtl_frc_out(lon_out_idx,lat_out_idx)
          if (lnd_frc_dry_out(lon_out_idx,lat_out_idx) < 0.0 .and. &
               abs(lnd_frc_dry_out(lon_out_idx,lat_out_idx)) < eps_rlt) &
               lnd_frc_dry_out(lon_out_idx,lat_out_idx)=0.0
          
       end do ! end loop over lon
    end do ! end loop over lat
    
    ! Sanity check
    do lat_out_idx=1,lat_out_nbr
       do lon_out_idx=1,lon_out_nbr
          frc_ttl=wtl_frc_out(lon_out_idx,lat_out_idx)+ &
               lak_frc_out(lon_out_idx,lat_out_idx)+ &
               lnd_frc_dry_out(lon_out_idx,lat_out_idx)
          if ( &
               (wtl_frc_out(lon_out_idx,lat_out_idx) < 0.0).or. &
               (wtl_frc_out(lon_out_idx,lat_out_idx) > 1.0).or. &
               (lak_frc_out(lon_out_idx,lat_out_idx) < 0.0).or. &
               (lak_frc_out(lon_out_idx,lat_out_idx) > 1.0).or. &
               (lnd_frc_dry_out(lon_out_idx,lat_out_idx) < 0.0).or. &
               (lnd_frc_dry_out(lon_out_idx,lat_out_idx) > 1.0).or. &
               ((sfc_typ_LSM_out(lon_out_idx,lat_out_idx) > 0).and. &
               (abs(frc_ttl-1.0) >= eps_rlt)) &
               ) then
             write (6,"(a)") "lnd_frc_dry_get(): ERROR during output"
             write (6,"(a)") "Wetland properties:"
             write (6,"(4(a,7x))") "Wetland","Lake","Dry","Total"
             write (6,"(4(es14.7,1x))") &
                  wtl_frc_out(lon_out_idx,lat_out_idx), &
                  lak_frc_out(lon_out_idx,lat_out_idx), &
                  lnd_frc_dry_out(lon_out_idx,lat_out_idx), &
                  frc_ttl
             write (6,"(a)") "Output cell edge locations:"
             write (6,"(a)") "(idx)   Lat Sth (idx)   Lat Nrt (idx)   Lon Wst (idx)   Lon Est"
             write (6,"(4(a1,i3,a1,1x,f9.4,1x))") &
                  "(",lat_out_idx,")",lat_out_grd(lat_out_idx),"(",lat_out_idx+1,")",lat_out_grd(lat_out_idx+1), &
                  "(",lon_out_idx,")",lon_out_grd(lon_out_idx),"(",lon_out_idx+1,")",lon_out_grd(lon_out_idx+1)
             stop
          endif               ! endif err
       end do ! end loop over lon
    end do ! end loop over lat
    
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Exiting lnd_frc_dry_get()"
    return 
  end subroutine lnd_frc_dry_get                       ! end lnd_frc_dry_get()
  
end module lak_wtl ! [mdl] Lakes and wetlands
