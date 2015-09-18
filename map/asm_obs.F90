! $Id$ -*-f90-*-

! Purpose: Assimilate observations

! Usage: 
! use asm_obs ! [mdl] Assimilate observations

module asm_obs ! [mdl] Assimilate observations
  implicit none
  private ! [stt] Symbols are private unless individually qualified as public
  public::toms_get ! [sbr] Retrieve and process TOMS satellite data
  public::lnd_frc_get ! [sbr] Determine land fraction based on lat/lon
  
contains
  
  subroutine toms_get( & ! [sbr] Retrieve and process TOMS satellite data
       drc_in,yr_nbr,yr_srt, & ! I
       lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr,time_in_nbr, & ! I
       lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr,time_out_nbr, & ! I
       area_out,lnd_frc_out,lnd_msk_out,ovr_nbr_max, & ! I
       fsh_fct,odxc_tms_out,src_flg,src_frq,src_str) ! O 
    ! Purpose: Retrieve optical depth data from satellite climatology
    ! toms_get() is called by asm_drv()
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use bds_ctl,only:caseid_tms,src_thr ! [mdl] Control variables drc_in,drc_out,hgt_dlt_msl
    use map_cst ! [mdl] Constants used in map routines
    use map_grd ! [mdl] Map grids and regridding
    use netcdf ! [mdl] netCDF interface
    use nf90_utl ! [mdl] netCDF utilities
    use src_prc,only:src_id ! [mdl] Source processing, identification
    use sng_mdl,only:ftn_strlen,ftn_strcpy,ftn_strcat ! [mdl] String manipulation
    use utl_mdl,only:mnt_chk ! [mdl] Utility functions (date_time_get,mnt_chk...)
    implicit none
    ! Parameters
    character(len=*),parameter::sbr_nm="toms_get" ! [sng] Subroutine name
    ! Commons
    ! Input
    character(len=*) drc_in      ! [sng] Input directory
    integer lat_in_nbr        ! [nbr] Number of latitudes
    integer lat_out_nbr       ! [nbr] Number of latitudes
    integer lon_in_nbr        ! [nbr] Number of longitudes
    integer lon_out_nbr       ! [nbr] Number of longitudes
    integer lnd_msk_out(lon_out_nbr,lat_out_nbr) ! [msk] Land mask (integer 0 or 1)
    integer ovr_nbr_max       ! [nbr] Maximum number of input cells which overlap any output cell
    integer time_in_nbr       ! [nbr] Number of records (months) in input file
    integer time_out_nbr      ! [nbr] Dimension size
    integer yr_nbr            ! [nbr] Number of years
    integer yr_srt            ! [yr] Starting year in YYYY format
    real area_out(lon_out_nbr,lat_out_nbr) ! [m2] Area of gridcells
    real lat_in_grd(lat_in_nbr+1) ! [dgr] Interface latitudes
    real lat_out_grd(lat_out_nbr+1) ! [dgr] Interface latitudes
    real lnd_frc_out(lon_out_nbr,lat_out_nbr) ! Land fraction of gridcell 
    real lon_in_grd(lon_in_nbr+1) ! [dgr] Interface longitudes
    real lon_out_grd(lon_out_nbr+1) ! [dgr] Interface longitudes
    ! Output
    integer src_flg(lon_out_nbr,lat_out_nbr) ! Source flag
    real fsh_fct(lon_out_nbr,lat_out_nbr) ! Efficiency factor
    real odxc_tms_out(lon_out_nbr,lat_out_nbr,time_out_nbr) ! [frc] Optical depth
    real src_frq(lon_out_nbr,lat_out_nbr,time_out_nbr) ! [frc] Source frequency
    real src_str(lon_out_nbr,lat_out_nbr,time_out_nbr) ! [frc] Source strength
    ! Local
    character(80)::fl_in        ! [sng] Input file
    character(2)::mth_sng       ! [sng] Month in MM format
    character(4)::yr_sng        ! [sng] Year in YYYY format
    integer lat_out_idx       ! [idx] Counting index for lat
    integer lon_out_idx       ! [idx] Counting index for lon
    integer mss_val_id        ! [id] Attribute ID
    integer nc_id             ! [id] File handle
    integer odxc_id           ! [id] Variable ID
    integer odxc_tll(lon_out_nbr,lat_out_nbr,time_out_nbr) ! Tally of valid optical depth months
    integer ovr_idx           ! [idx] Counting index
    integer ovr_lat_idx(lon_out_nbr,lat_out_nbr,ovr_nbr_max) ! [idx] Map into input grid of latitude indices of overlap cells
    integer ovr_lon_idx(lon_out_nbr,lat_out_nbr,ovr_nbr_max) ! [idx] Map into input grid of longitude indices of overlap cells
    integer ovr_nbr(lon_out_nbr,lat_out_nbr) ! [nbr] Number of input gridcells which overlap each output gridcell
    integer ovr_nbr_crr       ! [nbr] Current number of overlapping gridcells
    integer rcd               ! [enm] Return success code
    integer src_flg_crr(lon_out_nbr,lat_out_nbr) ! Source flag
    integer time_dmn_id       ! [id] Dimension ID for time
    integer time_id           ! [id] Coordinate ID
    integer time_in(time_in_nbr) ! Dimension
    integer time_out_idx      ! [idx] Counting index for time
    integer yr_crr            ! [yr] Current year in YYYY format
    integer yr_idx            ! [idx] Counting index
    logical mnt               ! [flg] Monotonicity flag
    logical mss_flg           ! [flg] Check for mss_val on input data
    real mss_val              ! Missing value
    real odxc_prc(lon_out_nbr,lat_out_nbr) ! [frc] Optical depth 
    real odxc_in(lon_in_nbr,lat_in_nbr) ! [frc] Optical depth 
    real ovr_wgt(lon_out_nbr,lat_out_nbr,ovr_nbr_max) ! [frc] Weight of overlapping input gridcells onto each output gridcell
    
    ! Main Code
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Entering toms_get()"
    
    ! Check for monotonicity
    mnt=mnt_chk(lon_out_grd,lon_out_nbr+1)
    if (.not.mnt) stop "lon_out_grd not monotonic in toms_get()"
    mnt=mnt_chk(lat_out_grd,lat_out_nbr+1)
    if (.not.mnt) stop "lat_out_grd not monotonic in toms_get()"
    ! Check for increasing monotonicity
    if (lon_out_grd(2) < lon_out_grd(1)) stop "lon_out_grd not increasing in toms_get()"
    if (lat_out_grd(2) < lat_out_grd(1)) stop "lat_out_grd not increasing in toms_get()"
    
    ! Initialize output arrays
    do lon_out_idx=1,lon_out_nbr
       do lat_out_idx=1,lat_out_nbr
          src_flg(lon_out_idx,lat_out_idx)=0
          fsh_fct(lon_out_idx,lat_out_idx)=0.0
          do time_out_idx=1,time_out_nbr
             odxc_tms_out(lon_out_idx,lat_out_idx,time_out_idx)=0.0
             odxc_tll(lon_out_idx,lat_out_idx,time_out_idx)=0
             src_frq(lon_out_idx,lat_out_idx,time_out_idx)=0.0
             src_str(lon_out_idx,lat_out_idx,time_out_idx)=0.0
          end do ! end loop over time_out
       end do ! end loop over lat
    end do ! end loop over lon
    
    ! Get overlap locations and weights for mapping input grid to output grid
    call map_ovr_wgt_drv( &
         lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
         lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
         ovr_nbr_max,area_out, & ! I
         ovr_lat_idx,ovr_lon_idx,ovr_nbr,ovr_wgt) ! O
    
    ! Loops over months and years
    do yr_idx=1,yr_nbr
       yr_crr=yr_srt+yr_idx-1
       write (yr_sng,"(i4.4)") yr_crr
       do time_out_idx=1,time_out_nbr
          write (mth_sng,"(i2.2)") time_out_idx
          fl_in=char(0)
          call ftn_strcpy(fl_in,drc_in)
          call ftn_strcat(fl_in,"/"//caseid_tms)
          call ftn_strcat(fl_in,"_"//yr_sng//mth_sng//".nc")
          ! Read netCDF data
          rcd=nf90_noerr        ! nf90_noerr == 0
          rcd=nf90_wrp_open(fl_in,NF90_NOWRITE,nc_id,sbr_nm=sbr_nm)
          ! Get dimension IDs
          rcd=nf90_wrp_inq_dimid(nc_id,"time",time_dmn_id)
          ! Get dimension sizes
          rcd=rcd+nf90_inquire_dimension(nc_id,time_dmn_id,len=time_in_nbr)
          if (time_in_nbr > 1) then
             write (6,"(2a,f12.6)") prg_nm(1:ftn_strlen(prg_nm)), & 
                  ": ERROR toms_get() reports time_in_nbr = ",time_in_nbr
             stop "time_in_nbr > 1"
          endif               ! endif err
          ! Get variable IDs
          rcd=nf90_wrp_inq_varid(nc_id,"time",time_id)
          ! fxm: Read absorbing aerosol index not guesstimated optical depth at 630 nm
          rcd=nf90_wrp_inq_varid(nc_id,"aer_idx_331_360",odxc_id)
          ! Get data
          rcd=nf90_wrp(nf90_get_var(nc_id,time_id,time_in),"get_var time_in")
          rcd=nf90_wrp(nf90_get_var(nc_id,odxc_id,odxc_in),"get_var odxc_in")
          if (rcd /= nf90_noerr) write (6,"(4a,i4)") prg_nm(1:ftn_strlen(prg_nm)), & 
               ": ERROR toms_get() reports rcd = ",rcd
          rcd=rcd+nf90_inquire_attribute(nc_id,odxc_id,"missing_value",attnum=mss_val_id)
          if (rcd==nf90_noerr) then
             mss_flg=.true.
          else
             mss_flg=.false.
             write (6,"(3a)") prg_nm(1:ftn_strlen(prg_nm)), & 
                  ": WARNING toms_get() reports no missing_value in ",fl_in(1:ftn_strlen(fl_in))
             rcd=nf90_noerr
          endif               ! endif
          rcd=rcd+nf90_get_att(nc_id,odxc_id,"missing_value",mss_val)
          if (real(mss_val) /= real(1.0e36)) write (6,"(2a,es15.8)") prg_nm(1:ftn_strlen(prg_nm)), & 
               ": WARNING toms_get() reports mss_val = ",mss_val
          rcd=nf90_wrp_close(nc_id,fl_in,'',sbr_nm=sbr_nm)
          if (dbg_lvl >= dbg_fl) write (6,"(a8,1x,a)") "Ingested",fl_in(1:ftn_strlen(fl_in))
          
          ! Rebin optical depth from input grid to output grid
          call map_rbn(odxc_in,mss_flg,mss_val, & ! I
               lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
               lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
               ovr_nbr_max,ovr_lat_idx,ovr_lon_idx,ovr_nbr,ovr_wgt, & ! I
               odxc_prc) ! O
          
          ! Identify source regions
          call src_id(odxc_prc,mss_flg,mss_val, & ! I
               lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
               lnd_msk_out, & ! I
               odxc_tms_out(1,1,time_out_idx),odxc_tll(1,1,time_out_idx),src_frq(1,1,time_out_idx),src_str(1,1,time_out_idx), & ! I/O
               src_flg_crr) ! O
          
          ! Take the union of the source flags
          do lon_out_idx=1,lon_out_nbr
             do lat_out_idx=1,lat_out_nbr
                if (src_flg_crr(lon_out_idx,lat_out_idx) > 0) &
                     src_flg(lon_out_idx,lat_out_idx)=src_flg_crr(lon_out_idx,lat_out_idx)
             end do ! end loop over lat
          end do ! end loop over lon
          
       end do ! end loop over mth
    end do ! end loop over yr
    
    ! Normalize optical depths
    do lon_out_idx=1,lon_out_nbr
       do lat_out_idx=1,lat_out_nbr
          do time_out_idx=1,time_out_nbr
             ! Compute seasonal cycle of observed optical depth, regardless of source
             if (odxc_tll(lon_out_idx,lat_out_idx,time_out_idx) /= 0) then
                odxc_tms_out(lon_out_idx,lat_out_idx,time_out_idx)= & ! fraction
                     odxc_tms_out(lon_out_idx,lat_out_idx,time_out_idx)/odxc_tll(lon_out_idx,lat_out_idx,time_out_idx)
             else
                odxc_tms_out(lon_out_idx,lat_out_idx,time_out_idx)=mss_val ! fraction
             endif            ! endif
          end do ! end loop over time_out
       end do ! end loop over lat
    end do ! end loop over lon
    
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Exiting toms_get()"
    return 
  end subroutine toms_get
  
  ! [sbr] Determine land fraction based on lat/lon
  subroutine lnd_frc_get(fl_in, & ! I
       lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
       lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
       ovr_nbr_max,area_out, & ! I
       lnd_frc_out,lnd_msk_out) ! O
    ! Purpose: Retrieve land fraction and mask and rebin to requested grid
    ! NB: Routine is deprecated
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use sng_mdl,only:ftn_strlen,ftn_strstr ! [mdl] String manipulation
    use utl_mdl,only:mnt_ncr_chk ! [mdl] Utility functions (date_time_get,mnt_chk...)
    implicit none
    ! Parameters
    integer,parameter::fl_in_unit=73 ! Unit for reading input data
    ! Commons
    ! Input
    character(80)::fl_in        ! [sng] Input file
    integer lat_in_nbr        ! [nbr] Number of latitudes
    integer lat_out_nbr       ! [nbr] Number of latitudes
    integer lon_in_nbr        ! [nbr] Number of longitudes
    integer lon_out_nbr       ! [nbr] Number of longitudes
    integer ovr_nbr_max       ! [nbr] Maximum number of input cells which overlap any output cell
    real area_out(lon_out_nbr,lat_out_nbr) ! [m2] Area of gridcells
    real lat_in_grd(lat_in_nbr+1) ! [dgr] Interface latitudes
    real lat_out_grd(lat_out_nbr+1) ! [dgr] Interface latitudes
    real lon_in_grd(lon_in_nbr+1) ! [dgr] Interface longitudes
    real lon_out_grd(lon_out_nbr+1) ! [dgr] Interface longitudes
    ! Output
    real lnd_frc_out(lon_out_nbr,lat_out_nbr) ! Land fraction of gridcell 
    integer lnd_msk_out(lon_out_nbr,lat_out_nbr) ! [msk] Land mask (integer 0 or 1)
    
    ! Local
    integer int_foo
    integer lat_in_idx        ! [idx] Counting index for lat
    integer lat_out_idx       ! [idx] Counting index for lat
    integer lon_in_idx        ! [idx] Counting index for lon
    integer lon_out_idx       ! [idx] Counting index for lon
    integer mss_val           ! [frc] Code for missing/invalid data
    logical mnt_ncr           ! [flg] Monotonic and increasing flag
    real area_in(lon_in_nbr,lat_in_nbr) ! [m2] Area of gridcells
    
    ! Main Code
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Entering lnd_frc_get()"
    
    ! Check for monotonically increasing grids
    mnt_ncr=mnt_ncr_chk(lon_out_grd,lon_out_nbr+1)
    if (.not.mnt_ncr) stop "lon_out_grd not monotonically increasing in lnd_frc_get()"
    mnt_ncr=mnt_ncr_chk(lat_out_grd,lat_out_nbr+1)
    if (.not.mnt_ncr) stop "lat_out_grd not monotonically increasing in lnd_frc_get()"
    
    ! Initialize scalars
    
    ! Initialize arrays
    do lon_out_idx=1,lon_out_nbr
       do lat_out_idx=1,lat_out_nbr
          lnd_frc_out(lon_out_idx,lat_out_idx)=0.0
          lnd_msk_out(lon_out_idx,lat_out_idx)=0
       end do ! end loop over lat
    end do ! end loop over lon
    
    ! For now, just use Shea"s subroutine
    do lon_out_idx=1,lon_out_nbr
       do lat_out_idx=1,lat_out_nbr
          ! call lndsea(lon_out(lon_out_idx),lat_out(lat_out_idx),int_foo)
          ! lnd_msk_out(lon_out_idx,lat_out_idx)=int_foo
          ! lnd_frc_out(lon_out_idx,lat_out_idx)=real(int_foo)
       end do ! end loop over lat
    end do ! end loop over lon
    
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Exiting lnd_frc_get()"
    return 
  end subroutine lnd_frc_get
  
end module asm_obs ! [mdl] Assimilate observations
