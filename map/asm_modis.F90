! $Id$ -*-f90-*-

! Purpose: Erodibility factors from MODIS data

! Usage: 
! use asm_modis ! [mdl] Erodibility factor from MODIS satellite

module asm_modis ! [mdl] Erodibility factor from MODIS satellite
  implicit none
  private ! [stt] Symbols are private unless individually qualified as public
  public::bsn_mds_get ! [sbr] Create erodibility factors from MODIS surface reflectivity
contains
  
  subroutine bsn_mds_get( & ! [sbr] Create erodibility factors from MODIS surface reflectivity
       fl_in, & ! I
       lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
       lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
       ovr_nbr_max,area_out,lnd_frc_out, & ! I
       rfl_sfc_lnd_nrm_mds_lnr_out,rfl_sfc_lnd_nrm_mds_sqr_out) ! O
    ! Purpose: Create erodibility factors from MODIS surface reflectivity
    ! rfl_sfc_lnd_nrm_mds_lnr is (albedo/albedo_max)
    ! rfl_sfc_lnd_nrm_mds_sqr is (albedo/albedo_max)**2
    ! Code based on soi_txt.F90, largest arrays are made allocatable to save memory
    ! Author: Alf Grini, alf.grini@geofysikk.uio.no (2003)
    use dbg_mdl  ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use map_grd  ! [mdl] Map grids and regridding
    use netcdf   ! [mdl] netCDF module
    use nf90_utl ! [mdl] netCDF utilities
    use sng_mdl,only:ftn_strlen,ftn_strstr ! [mdl] String manipulation
    use utl_mdl,only:mnt_ncr_chk ! [mdl] Utility functions (date_time_get,mnt_chk...)
    implicit none
    ! Parameters
    character(len=*),parameter::sbr_nm="bsn_mds_get" ! [sng] Subroutine name
    ! Input
    character(len=*),intent(in)::fl_in ! [sng] Input file, MODIS surface reflectance MOD09
    integer,intent(in)::lat_in_nbr ! [nbr] Number of latitudes
    integer,intent(in)::lat_out_nbr ! [nbr] Number of latitudes
    integer,intent(in)::lon_in_nbr ! [nbr] Number of longitudes
    integer,intent(in)::lon_out_nbr ! [nbr] Number of longitudes
    integer,intent(in)::ovr_nbr_max ! [nbr] Maximum number of input cells which overlap any output cell
    real,intent(in)::area_out(lon_out_nbr,lat_out_nbr) ! [m2] Area of gridcells
    real,intent(in)::lat_in_grd(lat_in_nbr+1) ! [dgr] Interface latitudes
    real,intent(in)::lat_out_grd(lat_out_nbr+1) ! [dgr] Interface latitudes
    real,intent(in)::lnd_frc_out(lon_out_nbr,lat_out_nbr) ! [frc] Land fraction of gridcell 
    real,intent(in)::lon_in_grd(lon_in_nbr+1) ! [dgr] Interface longitudes
    real,intent(in)::lon_out_grd(lon_out_nbr+1) ! [dgr] Interface longitudes
    ! Output
    real,intent(out)::rfl_sfc_lnd_nrm_mds_lnr_out(lon_out_nbr,lat_out_nbr) ! [frc] Erodibility factor from MODIS, linear
    real,intent(out)::rfl_sfc_lnd_nrm_mds_sqr_out(lon_out_nbr,lat_out_nbr) ! [frc] Erodibility factor from MODIS, squared
    ! Locals with simple initialization
    integer::rcd=nf90_noerr ! [enm] Return success code

    ! Local
    integer idx               ! [idx] Counting index
    integer lat_in_idx        ! [idx] Counting index for lat
    integer lat_out_idx       ! [idx] Counting index for lat
    integer lon_in_idx        ! [idx] Counting index for lon
    integer lon_out_idx       ! [idx] Counting index for lon
    integer ovr_idx           ! [idx] Counting index
    integer,dimension(:,:,:),allocatable::ovr_lat_idx ! [idx] Map into input grid of latitude indices of overlap cells
    integer,dimension(:,:,:),allocatable::ovr_lon_idx ! [idx] Map into input grid of longitude indices of overlap cells
    integer ovr_nbr(lon_out_nbr,lat_out_nbr) ! [nbr] Number of input gridcells which overlap each output gridcell
    integer ovr_nbr_crr       ! [nbr] Current number of overlapping gridcells
    logical mnt_ncr           ! [flg] Monotonic and increasing flag
    logical mss_flg           ! [flg] Check for mss_val on input data
    real frc_ttl              ! [frc] Total fraction of gridcell accounted for by all soil textures
    real mss_val              ! Missing value
    real,dimension(:,:,:),allocatable::ovr_wgt ! [frc] Weight of overlapping input gridcells onto each output gridcell
    real ovr_wgt_crr          ! [frc] Overlap weight of current gridcell
    real,dimension(:,:),allocatable::rfl_sfc_in  !Surface reflectance in input resolution
    real,dimension(:,:),allocatable::rfl_sfc_out !Surface reflectance in output resolution
    real max_rfl_sfc ! Maximum surface reflection in file

    ! Local netCDF read file stuff
    integer :: nc_id                 !netCDF file handle
    integer :: lat_id,lon_id         !Variable id for lon/lat
    integer :: lat_dmn_id,lon_dmn_id !Dimension id for lon/lat
    integer :: rfl_sfc_id            !Variable id for surface reflection
    real lon(lon_in_nbr)             !Longitude values in netCDF file
    real lat(lat_in_nbr)             !Latitude values in netCDF file
    integer srt_lon_lat(2)           !Start value for lon/lat
    integer cnt_lon_lat(2)           !Counting value for lon/lat
    integer :: lat_nbr               !Number of latitudes read from netCDF file
    integer :: lon_nbr               !Number of longitudes read from netCDF file

    ! Main Code
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Entering bsn_mds_get()"

    ! Check for monotonically increasing grids
    mnt_ncr=mnt_ncr_chk(lon_out_grd,lon_out_nbr+1)
    if (.not.mnt_ncr) stop "lon_out_grd not monotonically increasing in bsn_mds_get()"
    mnt_ncr=mnt_ncr_chk(lat_out_grd,lat_out_nbr+1)
    if (.not.mnt_ncr) stop "lat_out_grd not monotonically increasing in bsn_mds_get()"
    
    ! Initialize scalars
    mss_flg=.false.
    mss_val=-99999

    !Allocate memory for all arrays needed
    allocate(ovr_lat_idx(lon_out_nbr,lat_out_nbr,ovr_nbr_max),stat=rcd)
    if(rcd/=0) stop "allocate failed for ovr_lat_idx"
    allocate(ovr_lon_idx(lon_out_nbr,lat_out_nbr,ovr_nbr_max),stat=rcd)
    if(rcd/=0) stop "allocate failed for ovr_lon_idx"
    allocate(ovr_wgt(lon_out_nbr,lat_out_nbr,ovr_nbr_max),stat=rcd)
    if(rcd/=0) stop "allocate failed for ovr_wgt"
    allocate(rfl_sfc_in(lon_in_nbr,lat_in_nbr),stat=rcd)
    if(rcd/=0) stop "allocate failed for rfl_sfc_in"
    allocate(rfl_sfc_out(lon_out_nbr,lat_out_nbr),stat=rcd)
    if(rcd/=0) stop "allocate failed for rfl_sfc_out"
 
    ! Initialize output arrays
    rfl_sfc_lnd_nrm_mds_lnr_out(:,:)=0.0 ! [frc] Erodibility factor from MODIS, linear
    rfl_sfc_lnd_nrm_mds_sqr_out(:,:)=0.0 ! [frc] Erodibility factor from MODIS, squared
       
    ! Read MODIS surface reflection (0-1) annual mean, 2001 from netCDF file
    rcd=nf90_wrp_open(fl_in,nf90_nowrite,nc_id,sbr_nm=sbr_nm)
    ! Get dimension IDs
    rcd=nf90_wrp_inq_dimid(nc_id,"lat",lat_dmn_id)
    rcd=nf90_wrp_inq_dimid(nc_id,"lon",lon_dmn_id)

    ! Inquire dimensions
    rcd=nf90_wrp(nf90_inquire_dimension(nc_id,lat_dmn_id,len=lat_nbr),"bsn_mds_get: inquire_dim lat")
    if (lat_nbr /= lat_in_nbr)then
       write(6,*)"bsn_modis_get: lat_nbr /= lat_in_nbr it is",lat_nbr
       stop
    endif
    rcd=nf90_wrp(nf90_inquire_dimension(nc_id,lon_dmn_id,len=lon_nbr),"bsn_mds_get: inquire_dim lon")
    if (lon_nbr /= lon_in_nbr)then
       write(6,*)"bsn_modis_get: lon_nbr /= lon_in_nbr, it is",lon_nbr
       stop
    endif

    ! Get variable IDs
    rcd=nf90_wrp_inq_varid(nc_id,"lon",lon_id)
    rcd=nf90_wrp_inq_varid(nc_id,"lat",lat_id)
    rcd=nf90_wrp_inq_varid(nc_id,"SURFREF7",rfl_sfc_id)

    ! Get data
    rcd=nf90_wrp(nf90_get_var(nc_id,lat_id,lat),sbr_nm//": get_var lat")
    rcd=nf90_wrp(nf90_get_var(nc_id,lon_id,lon),sbr_nm//": get_var lon")
    rcd=nf90_wrp(nf90_get_var(nc_id,rfl_sfc_id,rfl_sfc_in),sbr_nm//": get_var rfl_sfc_in")

    ! Close file
    rcd=nf90_wrp_close(nc_id,fl_in,"Ingested") ! [fnc] Close file

    ! Done reading
    do lat_in_idx=1,lat_in_nbr
       do lon_in_idx=1,lon_in_nbr

          ! Sanity check
          if(rfl_sfc_in(lon_in_idx,lat_in_idx).gt.1.0.or. &
               rfl_sfc_in(lon_in_idx,lat_in_idx).lt.0.0)then
             write(6,*)'ERROR IN READ SURFACE REFLECTANCE' 
             write(6,*)'lon_in_idx,lon_lat_idx value',   &
                  lon_in_idx,lat_in_idx,rfl_sfc_in(lon_in_idx,lat_in_idx)
             stop
          endif  !sanity check
       end do ! end loop over lon
    end do ! end loop over lat
              
    ! Get overlap locations and weights for mapping input grid to output grid
    call map_ovr_wgt_drv( &
         lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
         lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
         ovr_nbr_max,area_out, & ! I
         ovr_lat_idx,ovr_lon_idx,ovr_nbr,ovr_wgt) ! O
    
    ! Rebin Surface reflectivity to output grid
    call map_rbn(rfl_sfc_in,mss_flg,mss_val, & ! I
         lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
         lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
         ovr_nbr_max,ovr_lat_idx,ovr_lon_idx,ovr_nbr,ovr_wgt, & ! I
         rfl_sfc_out) ! O
    
    ! Order or normalization is important
    ! Can rebin then normalize or normalize then rebin
    ! These are not equivalent:
    ! 1. Average albedo and then make erodibility factors or 
    ! 2. Make erodibility factors and then average them !!!

    ! Normalize surface reflectivity as fraction of maximum
    max_rfl_sfc=maxval(rfl_sfc_out)
    do lat_out_idx=1,lat_out_nbr
       do lon_out_idx=1,lon_out_nbr
          ! Linear MODIS erodibility factor
          rfl_sfc_lnd_nrm_mds_lnr_out(lon_out_idx,lat_out_idx)= &
               rfl_sfc_out(lon_out_idx,lat_out_idx)/max_rfl_sfc
          ! Quadratic MODIS erodibility factor
          rfl_sfc_lnd_nrm_mds_sqr_out(lon_out_idx,lat_out_idx)= &
               (rfl_sfc_out(lon_out_idx,lat_out_idx)/max_rfl_sfc)**2
       enddo ! end loop over lon
    enddo ! end loop over lat
        
    ! Zero basin factors over ocean
    where (lnd_frc_out == 0.0) ! [frc]
       rfl_sfc_lnd_nrm_mds_lnr_out=0.0 ! [frc] Erodibility factor from MODIS, linear
       rfl_sfc_lnd_nrm_mds_sqr_out=0.0 ! [frc] Erodibility factor from MODIS, squared
    end where ! end where lnd_frc_out

    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Exiting "//sbr_nm
    return 
  end subroutine bsn_mds_get ! end bsn_mds_get()
  
end module asm_modis ! [mdl] Erodibility factor from MODIS satellite
