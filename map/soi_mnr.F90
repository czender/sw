! $Id$ -*-f90-*-

! Purpose: Module soi_mnr stores parameters and routies related to soil mineralogy 

! Usage:
! use soi_mnr ! [mdl] Soil mineralogy

module soi_mnr ! [mdl] Soil mineralogy
  implicit none
  
  integer,parameter::chm_nbr=7 ! [nbr] Number of chemical species to process
  
  integer,parameter::idx_Al=1 ! [idx] Exchangeable aluminum
  integer,parameter::idx_C_org=2 ! [idx] Organic carbon
  integer,parameter::idx_CaCO3=3 ! [idx] Calcium carbonate
  integer,parameter::idx_K=4 ! [idx] Exchangeable potassium
  integer,parameter::idx_N=5 ! [idx] Nitrogen
  integer,parameter::idx_Na=6 ! [idx] Exchangeable sodium
  integer,parameter::idx_P2O5=7 ! [idx] Extractable phosphorus

contains

  subroutine soi_chm_IGBP_get(fl_in,chm_idx, & ! I
       lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
       lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
       ovr_nbr_max,area_out,lnd_msk_out, & ! I
       mss_frc_chm_out) ! O [kg kg-1] Generic chemical
    ! Purpose: Retrieve soil chemical composition from IGBP-DIS surface files and rebin to requested grid
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use map_grd ! [mdl] Map grids and regridding
    use sng_mdl,only:ftn_strlen,ftn_strcpy,ftn_strnul,ftn_strcat ! [mdl] String manipulation
    use utl_mdl,only:mnt_ncr_chk ! [mdl] Utility functions (date_time_get,mnt_chk...)
    implicit none
    ! Parameters
    character(len=*),parameter::nlc=char(0) ! [sng] NUL character = ASCII 0 = char(0)
    integer,parameter::fl_in_unit=73 ! Unit for reading input data
    integer,parameter::lat_rgn_in_skip=6 ! [nbr] Number of latitudes skipped before input data
    integer,parameter::lat_rgn_in_nbr=140 ! [nbr] Number of latitudes in input data region
    integer,parameter::lon_rgn_in_nbr=360 ! [nbr] Number of longitudes in input data region
    integer,parameter::unit_dgn=74 ! Unit for writing general Map diagnoses
    ! Commons
    ! Input
    character(len=*),intent(in)::fl_in ! [sng] Input file
    integer,intent(in)::chm_idx           ! [idx] Counting index for chm
    integer,intent(in)::lat_in_nbr        ! [nbr] Number of latitudes
    integer,intent(in)::lat_out_nbr       ! [nbr] Number of latitudes
    integer,intent(in)::lon_in_nbr        ! [nbr] Number of longitudes
    integer,intent(in)::lon_out_nbr       ! [nbr] Number of longitudes
    integer,intent(in)::ovr_nbr_max       ! [nbr] Maximum number of input cells which overlap any output cell
    integer,intent(in)::lnd_msk_out(lon_out_nbr,lat_out_nbr) ! [flg] Land mask (integer 0 or 1)
    real,intent(in)::area_out(lon_out_nbr,lat_out_nbr) ! [m2] Area of gridcells
    real,intent(in)::lat_in_grd(lat_in_nbr+1) ! [dgr] Interface latitudes
    real,intent(in)::lat_out_grd(lat_out_nbr+1) ! [dgr] Interface latitudes
    real,intent(in)::lon_in_grd(lon_in_nbr+1) ! [dgr] Interface longitudes
    real,intent(in)::lon_out_grd(lon_out_nbr+1) ! [dgr] Interface longitudes
    ! Output
    real,intent(out)::mss_frc_chm_out(lon_out_nbr,lat_out_nbr) ! [kg kg-1] Generic chemical
    ! Local
    character chm_sng(chm_nbr)*32 ! [sng] Array of species names
    character chm_sng_crr*32  ! [sng] Current species name
    character fl_chm*80       ! IGBP-DIS SoilData Data-Surface file
    character fl_dgn_chm*80   ! Chemical diagnostic file 
    integer idx               ! [idx] Counting index
    integer lat_rgn_in_idx    ! [idx] Counting index for lat
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
    logical mnt_ncr           ! Monotonic and increasing flag
    logical mss_flg           ! [flg] Check for mss_val on input data
    real chm_dat_in(lon_rgn_in_nbr,lat_rgn_in_nbr) ! Real array for input
    real area_in(lon_in_nbr,lat_in_nbr) ! [m2] Area of gridcells
    real eps_rlt              ! [frc] Relative error allowed in frc_ttl
    real frc_ttl              ! [frc] Total fraction of gridcell accounted for by all soil textures
    real mss_frc_chm_in(lon_in_nbr,lat_in_nbr) ! Chemical field 
    real mss_val_out          ! Missing value
    real mss_val_in_ocn       ! Missing value used for Ocean in input data
    real mss_val_in_lnd       ! Missing value used for Land in input data
    real ovr_wgt(lon_out_nbr,lat_out_nbr,ovr_nbr_max) ! [frc] Weight of overlapping input gridcells onto each output gridcell
    real ovr_wgt_crr          ! [frc] Overlap weight of current gridcell
    real val_min(chm_nbr) ! [kg kg-1] Geochemical species minimum values
    real val_max(chm_nbr) ! [kg kg-1] Geochemical species maximum values
    ! Main Code
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Entering soi_chm_IGBP_get()"
    
    ! Check for monotonically increasing grids
    mnt_ncr=mnt_ncr_chk(lon_out_grd,lon_out_nbr+1)
    if (.not.mnt_ncr) stop "lon_out_grd not monotonically increasing in soi_chm_IGBP_get()"
    mnt_ncr=mnt_ncr_chk(lat_out_grd,lat_out_nbr+1)
    if (.not.mnt_ncr) stop "lat_out_grd not monotonically increasing in soi_chm_IGBP_get()"
    
    ! Initialize scalars
    eps_rlt=1.0e-5            ! [frc] Relative error allowed in frc_ttl
    mss_flg=.true.
    ! fxm: set to zero so that global models may use value directly
    ! mss_val_out=1.0e36
    mss_val_out=0.0 !
    mss_val_in_ocn=-2.0 ! [] Missing value used for Ocean in input data
    mss_val_in_lnd=-1.0 ! [] Missing value used for Land in input data
    
    chm_sng(idx_Al)="Al"//nlc ! [sng] Exchangeable aluminum
    chm_sng(idx_C_org)="C_org"//nlc ! [sng] Organic carbon
    chm_sng(idx_CaCO3)="CaCO3"//nlc ! [sng] Calcium carbonate
    chm_sng(idx_K)="K"//nlc ! [sng] Exchangeable potassium
    chm_sng(idx_N)="N"//nlc ! [sng] Nitrogen
    chm_sng(idx_Na)="Na"//nlc ! [sng] Exchangeable sodium
    chm_sng(idx_P2O5)="P2O5"//nlc ! [sng] Extractable phosphorus
    
    ! 19990805: Following line is first use of f90 arithmetic constructs
    val_min=0.0               ! [kg kg-1]
    val_max=1.0e20            ! [kg kg-1]
    val_max(idx_C_org)=1.0    ! [kg kg-1] Organic carbon maximum value
    val_max(idx_CaCO3)=1.0    ! [kg kg-1] Calcium carbonate maximum value
    val_max(idx_N)=1.0        ! [kg kg-1] Nitrogen maximum value
    
    ! Construct filename
    call ftn_strcpy(fl_chm,fl_in)
    call ftn_strcat(fl_chm,chm_sng(chm_idx))
    call ftn_strnul(fl_chm)
    call ftn_strcat(fl_chm,".srf")
    chm_sng_crr=nlc
    call ftn_strcpy(chm_sng_crr,chm_sng(chm_idx))
    
    ! Initialize output arrays
    do lat_out_idx=1,lat_out_nbr
       do lon_out_idx=1,lon_out_nbr
          mss_frc_chm_out(lon_out_idx,lat_out_idx)=mss_val_out ! [frc]
       end do ! end loop over lon
    end do ! end loop over lat
    do lat_in_idx=1,lat_in_nbr
       do lon_in_idx=1,lon_in_nbr
          mss_frc_chm_in(lon_in_idx,lat_in_idx)=mss_val_out ! [frc]
       end do ! end loop over lon
    end do ! end loop over lat
    
    ! Compute gridcell area 
    call map_area_get(lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
         area_in) ! O
    
    ! IGBP data are stored as a one column formatted text file
    ! The dataset resolution is determined by the operator of the SoilData program
    ! The finest possible resolution is 5 arc minutes
    ! My default is 1 x 1 degree resolution
    ! Each datum, i.e., soil mass percentage of CaCO3 at a given (latitude,longitude) point, is written alone on one row
    ! Global data, i.e., data created using the world.img map, are stored in latitude bands from 84 N to 56 S
    ! Each latitude band is stored West to East, beginning at the date line (180 W)
    ! Ocean data are stored as -2
    ! Void cells (non-ocean but no data) are stored as -1
    
    ! Read chemical species information (usually input as mass percentage)
    open (fl_in_unit,file=fl_chm,status="old",iostat=rcd)
    if (rcd /= 0) then 
       write (6,"(4a,i4)") prg_nm(1:ftn_strlen(prg_nm)), & 
            ": ERROR soi_chm_IGBP_get() unable to open ",fl_chm(1:ftn_strlen(fl_chm)), &
            ", iostat return code = ",rcd
       stop
    endif ! endif err
    do lat_rgn_in_idx=1,lat_rgn_in_nbr
       do lon_in_idx=1,lon_in_nbr
          ! fxm: Sun f90 requires explicit format to read floats (not ints) correctly
          ! But I cannot figure out a generic format for the floats in these files
          ! Thus map generation does not work on Sun machines
          ! character(80) chm_dat_sng ! [sng] String to hold input data
          
          ! read (fl_in_unit,"(a)",iostat=rcd) chm_dat_sng ! [sng] String to hold input data
          ! read (chm_dat_sng,*) chm_dat_in(lon_in_idx,lat_rgn_in_idx)
          read (fl_in_unit,*,iostat=rcd) chm_dat_in(lon_in_idx,lat_rgn_in_idx)
          ! read (fl_in_unit,"(f12.7)",iostat=rcd) chm_dat_in(lon_in_idx,lat_rgn_in_idx)
#if 0
          write (6,"(2(a,i3),a,i10,a,f12.7)") "chm_dat_in(",lon_in_idx,",", &
               lat_rgn_in_idx,") (1D offset = ", &
               (lat_rgn_in_idx-1)*lon_in_nbr+lon_in_idx,") = ", &
               chm_dat_in(lon_in_idx,lat_rgn_in_idx)
#endif /* endif false */
       end do ! end loop over lon
    end do ! end loop over lat
    if (rcd /= 0) then 
       write (6,"(4a,i4)") prg_nm(1:ftn_strlen(prg_nm)), & 
            ": ERROR soi_chm_IGBP_get() reading ",fl_chm(1:ftn_strlen(fl_chm)), &
            ", iostat return code = ",rcd
       stop
    endif ! endif err
    close (fl_in_unit)
    write (6,"(4a)") "Read ",chm_sng_crr(1:ftn_strlen(chm_sng_crr))," data from ",fl_chm(1:ftn_strlen(fl_chm))
    do lat_rgn_in_idx=1,lat_rgn_in_nbr
       do lon_in_idx=1,lon_in_nbr
          lat_in_idx=lat_rgn_in_idx+lat_rgn_in_skip-1
          ! Remap SoilData grid so latitude monotonically increases from South Pole
          lat_in_idx=lat_in_nbr-lat_in_idx+1
          if (chm_dat_in(lon_in_idx,lat_rgn_in_idx) >= 0.0) then
             ! Handle valid input data differently depending on input units 
             if (chm_idx==idx_CaCO3.or.chm_idx==idx_C_org.or.chm_idx==idx_N) then
                ! Input unit is mass percentage, convert to mass fraction
                mss_frc_chm_in(lon_in_idx,lat_in_idx)=chm_dat_in(lon_in_idx,lat_rgn_in_idx)*0.01 ! [pct] -> [frc]
             else 
                ! Input unit is weird--temporarily neglect conversion
                mss_frc_chm_in(lon_in_idx,lat_in_idx)=chm_dat_in(lon_in_idx,lat_rgn_in_idx) 
             endif            ! endif input is a mass percentage
          else if (chm_dat_in(lon_in_idx,lat_rgn_in_idx)==mss_val_in_ocn) then
             mss_frc_chm_in(lon_in_idx,lat_in_idx)=mss_val_out
          else if (chm_dat_in(lon_in_idx,lat_rgn_in_idx)==mss_val_in_lnd) then
             mss_frc_chm_in(lon_in_idx,lat_in_idx)=mss_val_out
          else
             write (6,"(f9.4)") chm_dat_in(lon_in_idx,lat_in_idx)
             stop "Bad chm_dat_in in soi_chm_IGBP_get()"
          endif               ! endif
       end do ! end loop over lon
    end do ! end loop over lat
    
    ! Sanity check on input
    do lat_in_idx=1,lat_in_nbr
       do lon_in_idx=1,lon_in_nbr
          if (mss_frc_chm_in(lon_in_idx,lat_in_idx) /= mss_val_out) then
             frc_ttl=mss_frc_chm_in(lon_in_idx,lat_in_idx)
             if ( &
                  (mss_frc_chm_in(lon_in_idx,lat_in_idx) < val_min(chm_idx).or. &
                  mss_frc_chm_in(lon_in_idx,lat_in_idx) > val_max(chm_idx)).or. &
                  (frc_ttl > val_max(chm_idx))) then
                write (6,"(a)") "ERROR soi_chm_IGBP_get() reports error during input"
                write (6,"(a)") "Soil composition properties:"
                write (6,"(2(a,7x))") chm_sng_crr(1:ftn_strlen(chm_sng_crr)),"Total"
                write (6,"(2(f12.6,1x))") &
                     mss_frc_chm_in(lon_in_idx,lat_in_idx), &
                     frc_ttl
                write (6,"(a)") "Input cell edge locations:" 
                write (6,"(a)") "(idx)   Lat Sth (idx)   Lat Nrt (idx)   Lon Wst (idx)   Lon Est"
                write (6,"(4(a1,i3,a1,1x,f9.4,1x))") &
                     "(",lat_in_idx,")",lat_in_grd(lat_in_idx),"(",lat_in_idx+1,")",lat_in_grd(lat_in_idx+1), &
                     "(",lon_in_idx,")",lon_in_grd(lon_in_idx),"(",lon_in_idx+1,")",lon_in_grd(lon_in_idx+1)
                stop
             end if           ! endif err
          endif               ! endif not missing value
       end do ! end loop over lon
    end do ! end loop over lat
    
    if (dbg_lvl==dbg_crr) then
       fl_dgn_chm=nlc ! Species diagnostics file
       call ftn_strcpy(fl_dgn_chm,"/tmp/zender/map/map_")
       call ftn_strcat(fl_dgn_chm,chm_sng(chm_idx))
       call ftn_strnul(fl_dgn_chm)
       call ftn_strcat(fl_dgn_chm,".txt")
       ! Write species data in CSM/CCM/LSM format (one value per line)
       open (unit_dgn,file=fl_dgn_chm,status="unknown",iostat=rcd)
       if (rcd /= 0) write (6,"(4a,i4)") prg_nm(1:ftn_strlen(prg_nm)), & 
            ": ERROR soi_chm_IGBP_get() unable to open ",fl_dgn_chm(1:ftn_strlen(fl_dgn_chm))
       do lat_in_idx=1,lat_in_nbr ! NB: Outer loop over lat
          do lon_in_idx=1,lon_in_nbr
             write (unit_dgn,*) mss_frc_chm_in(lon_in_idx,lat_in_idx)
          end do ! end loop over lon
       end do ! end loop over lat
       close (unit_dgn)
       write (6,"(4a)") "Wrote ",chm_sng_crr," data in LSM format to",fl_dgn_chm(1:ftn_strlen(fl_dgn_chm))
    endif                     ! endif dbg
    
    ! Get overlap locations and weights for mapping input grid to output grid
    call map_ovr_wgt_drv( &
         lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
         lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
         ovr_nbr_max,area_out, & ! I
         ovr_lat_idx,ovr_lon_idx,ovr_nbr,ovr_wgt) ! O
    
    ! Rebin values from input grid to output grid
    call map_rbn(mss_frc_chm_in,mss_flg,mss_val_out, & ! I
         lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
         lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
         ovr_nbr_max,ovr_lat_idx,ovr_lon_idx,ovr_nbr,ovr_wgt, & ! I
         mss_frc_chm_out) ! O
    
    ! Surface-type specific constraints on species
    do lat_out_idx=1,lat_out_nbr
       do lon_out_idx=1,lon_out_nbr
          ! Set oceans to zero mass fraction
          if (lnd_msk_out(lon_out_idx,lat_out_idx)==0) then
             mss_frc_chm_out(lon_out_idx,lat_out_idx)=0.0 ! [kg kg-1]
          endif               ! not Ocean     
       end do ! end loop over lon
    end do ! end loop over lat
    
    ! Sanity check on output
    do lat_out_idx=1,lat_out_nbr
       do lon_out_idx=1,lon_out_nbr
          if (mss_frc_chm_out(lon_out_idx,lat_out_idx) /= mss_val_out) then
             frc_ttl=mss_frc_chm_out(lon_out_idx,lat_out_idx)
             if ( &
                  (mss_frc_chm_out(lon_out_idx,lat_out_idx) < val_min(chm_idx)).or. &
                  (mss_frc_chm_out(lon_out_idx,lat_out_idx) > val_max(chm_idx)).or. &
                  ((lnd_msk_out(lon_out_idx,lat_out_idx) > 0).and. &
                  .false.) &
                  ) then
                write (6,"(a)") "soi_chm_IGBP_get(): ERROR during output"
                write (6,"(a)") "Soil composition properties:"
                write (6,"(2(a,7x))") chm_sng_crr(1:ftn_strlen(chm_sng_crr)),"Total"
                write (6,"(2(es14.7,1x))") &
                     mss_frc_chm_out(lon_out_idx,lat_out_idx), &
                     frc_ttl
                write (6,"(a)") "Output cell edge locations:" 
                write (6,"(a)") "(idx)   Lat Sth (idx)   Lat Nrt (idx)   Lon Wst (idx)   Lon Est"
                write (6,"(4(a1,i3,a1,1x,f9.4,1x))") &
                     "(",lat_out_idx,")",lat_out_grd(lat_out_idx),"(",lat_out_idx+1,")",lat_out_grd(lat_out_idx+1), &
                     "(",lon_out_idx,")",lon_out_grd(lon_out_idx),"(",lon_out_idx+1,")",lon_out_grd(lon_out_idx+1)
                stop
             endif            ! endif err
          endif               ! endif not missing value
       end do ! end loop over lon
    end do ! end loop over lat
    
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Exiting soi_chm_IGBP_get()"
    return 
  end subroutine soi_chm_IGBP_get ! end soi_chm_IGBP_get()
  
end module soi_mnr ! [mdl] Soil mineralogy
