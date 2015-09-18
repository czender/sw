! $Id$ -*-f90-*-

! Purpose: Routines to process surface height fields

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
! use hgt_sfc ! [mdl] Surface elevation, bathymetry

module hgt_sfc ! [mdl] Surface elevation, bathymetry
  implicit none
  public::hgt_sfc_get ! [sbr] Retrieve and regrid surface height
  
contains
  
  subroutine hgt_sfc_get( & ! [sbr] Retrieve and regrid surface height
       fl_in, & ! I
       lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
       area_out, & ! I
       hgt_sfc_gpm_out,hgt_sfc_out,hgt_sfc_std_dvn_out, & ! O
       lnd_frc_out,lnd_msk_out,oro_out, & ! O
       mbl_bsn_fct_out,bsn_enm_out,flw_acm_fct_out,sfc_acm_fct_out) ! O
    ! Purpose: Retrieve surface height and rebin to requested grid
    ! hgt_sfc_get() is called by bds_prc()
    ! hgt_sfc_get() evolved from hgt_sfc_get_old(), which has been deprecated
    ! hgt_sfc_get() can use elevation from either 
    ! NCAR SCD DSS754.0 is U.S. Navy Global Elevation 10-arcminute resolution dataset described at http://www.scd.ucar.edu/dss/catalogs/geo.html
    ! CSM topography file hgt_sfc_topo.nc was derived from DSS754.0
    ! This is the same as LSM elevation file mksrf_elev.nc
    ! Subroutine tools/getopo.F90 defines CSM method used
    ! CSM initial fields SGH, PHIS, and ORO are defined using mean and variance of height and fractional cover fields of this file
    ! Dataset resolution is 1/6 x 1/6 degree and data are written from South pole to North pole, East to West, beginning at Greenwich
    
    use bds_ctl,only:hgt_dlt_msl,bsn_fct_hrz ! [mdl] Control variables drc_in,drc_out,hgt_dlt_msl
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use geo_mdl,only:dst_src_flw ! [mdl] Compute flow accumulation factor
    use map_cst ! [mdl] Constants used in map routines
    use map_grd ! [mdl] Map grids and regridding
    use netcdf ! [mdl] netCDF interface
    use nf90_utl ! [mdl] netCDF utilities
    use sng_mdl,only:ftn_strlen ! [mdl] String manipulation
    use src_prc,only:rdb_fct_tpg_GCT01 ! [mdl] Source processing, identification
    use utl_mdl,only:mnt_ncr_chk ! [mdl] Utility functions (date_time_get,mnt_chk...)
    implicit none
    ! Parameters
    character(len=*),parameter::sbr_nm="hgt_sfc_get" ! [sng] Subroutine name
    integer,parameter::fl_in_unit=73 ! Unit for reading input data
    integer,parameter::lat_idx_dbg=50 ! [idx] Longitude for verbose output
    integer,parameter::lat_nbr_max=2160 ! [nbr] Maximum number of latitudes
    integer,parameter::lon_idx_dbg=90 ! [idx] Longitude for verbose output
    integer,parameter::lon_nbr_max=4320 ! [nbr] Maximum number of longitudes
    integer,parameter::unit_dgn=74 ! Unit for writing general Map diagnoses
    real,parameter::eps_rlt=1.0e-5 ! [frc] Relative error allowed in frc_ttl
    real,parameter::hgt_sfc_max_CLM=10000.0 ! [m] Maximum surface height
    real,parameter::hgt_dlt_msl_LGM_sis=-104.0 ! [m] Mean sea level change of LGM used by Sang-Ik Shin
    real,parameter::hgt_dlt_msl_LGM_std=-120.0 ! [m] Mean sea level change of LGM
    ! Derived parameters
    ! Commons
    ! Input
    character(len=*),intent(in)::fl_in ! [sng] Input file
    integer,intent(in)::lat_out_nbr ! [nbr] Number of latitudes
    integer,intent(in)::lon_out_nbr ! [nbr] Number of longitudes
    real,intent(in)::area_out(lon_out_nbr,lat_out_nbr) ! [m2] Area of gridcells
    real,intent(in)::lat_out_grd(lat_out_nbr+1) ! [dgr] Interface latitudes
    real,intent(in)::lon_out_grd(lon_out_nbr+1) ! [dgr] Interface longitudes
    ! Output
    integer,intent(out)::bsn_enm_out(lon_out_nbr,lat_out_nbr) ! [enm] Basin ID
    integer,intent(out)::lnd_msk_out(lon_out_nbr,lat_out_nbr) ! [flg] Land mask (integer 0 or 1)
    real,intent(out)::hgt_sfc_gpm_out(lon_out_nbr,lat_out_nbr) ! [gpm] Surface geopotential height 
    real,intent(out)::hgt_sfc_out(lon_out_nbr,lat_out_nbr) ! [m] Surface height 
    real,intent(out)::hgt_sfc_std_dvn_out(lon_out_nbr,lat_out_nbr) ! [m] Standard deviation of surface height 
    real,intent(out)::lnd_frc_out(lon_out_nbr,lat_out_nbr) ! [frc] Land fraction
    real,intent(out)::oro_out(lon_out_nbr,lat_out_nbr) ! [frc] Orography 
    real,intent(out)::mbl_bsn_fct_out(lon_out_nbr,lat_out_nbr) ! [frc] Mobilization enhancement due to basin characteristics
    real,intent(out)::flw_acm_fct_out(lon_out_nbr,lat_out_nbr) ! [frc] Flow accumulation factor
    real,intent(out)::sfc_acm_fct_out(lon_out_nbr,lat_out_nbr) ! [frc] Area accumulation factor
    
    ! Locals with simple initialization
    character(80)::fl_dgn="/tmp/zender/map/hgt_sfc.txt" ! [sng] Diagnostics file
    integer::map_typ_in=map_lat_rgl_lon_Grn_wst ! Latitudes are regular, Greenwich at West edge of first longitude bin
    integer::rcd=nf90_noerr ! [enm] Return success code
    logical::flg_read_2d_grd=.false. ! [flg] Read/use un-necessary CLM2 2D grid fields
    logical::flg_NGDC=.true. ! [flg] NGDC TerrainBase NCAR DSS 759.2 with CLM2 coding
    logical::flg_CCM=.false. ! [flg] Navy topography NCAR DSS 754.0 with CSM/CCM coding
    real::mss_val_in=1.0e36 ! [frc] Missing value
    real::mss_val_out=1.0e36 ! [frc] Missing value
    
    ! Local
    integer err_nbr ! [nbr] Number of errors processing surface types
    integer idx ! [idx] Counting index
    integer lat_in_idx ! [idx] Counting index for lat
    integer lat_in_idx_err ! [idx] Latitude of gridcell which caused error
    integer lat_in_nbr ! [nbr] Number of latitudes
    integer lat_out_idx ! [idx] Counting index for lat
    integer lat_out_idx_err ! [idx] Latitude of gridcell which caused error
    integer lon_in_idx ! [idx] Counting index for lon
    integer lon_in_idx_err ! [idx] Longitude of gridcell which caused error
    integer lon_in_nbr ! [nbr] Number of longitudes
    integer lon_out_idx ! [idx] Counting index for lon
    integer lon_out_idx_err ! [idx] Longitude of gridcell which caused error
    integer ovr_idx ! [idx] Counting index
    integer ovr_nbr(lon_out_nbr,lat_out_nbr) ! [nbr] Number of input gridcells which overlap each output gridcell
    integer ovr_nbr_crr ! [nbr] Current number of overlapping gridcells
    integer ovr_nbr_max ! [nbr] Maximum number of input cells which overlap any output cell
    logical lat_mnt_ncr ! [flg] Latitude monotonically decreases
    logical mnt_ncr ! [flg] Monotonic and increasing flag
    logical mss_flg ! [flg] Variable has missing_value attribute
    real edg_grd(4) ! [dgr] Grid edges (north, east, south, west)
    real frc_ttl ! [frc] Total fraction of gridcell accounted for by all plant functional types
    real hgt_sfc_var_out(lon_out_nbr,lat_out_nbr) ! [m2] Variance of surface height
    real ovr_wgt_crr ! [frc] Overlap weight of current gridcell
    
    ! Variables needed to read external netCDF dataset
    character(10)::hgt_sfc_nm ! [sng] Surface height variable name
    character(10)::lnd_frc_nm ! [sng] Fraction of land (not ocean) variable name
    integer edg_est_id ! [id] Variable ID
    integer edg_nrt_id ! [id] Variable ID
    integer edg_sth_id ! [id] Variable ID
    integer edg_wst_id ! [id] Variable ID
    integer hgt_sfc_dmn_nbr ! [nbr] Number of dimensions in disk hgt_sfc field
    integer lat_dmn_id ! [id] Dimension ID for lat
    integer lat_in_2d_id ! [id] Variable ID
    integer lat_in_id ! [id] Variable ID
    integer lnd_frc_id ! [id] Variable ID
    integer lon_dmn_id ! [id] Dimension ID for lon
    integer lon_in_2d_id ! [id] Variable ID
    integer mss_val_id ! [id] Attribute ID
    integer nc_id ! [id] File handle
    integer hgt_sfc_id ! [id] Variable ID
    
    ! Allocatables
    integer,dimension(:),allocatable::lon_in_nbr_1d ! [nbr] Number of longitudes per latitude
    integer,dimension(:,:,:),allocatable::ovr_lat_idx ! [idx] Map into input grid of latitude indices of overlap cells
    integer,dimension(:,:,:),allocatable::ovr_lon_idx ! [idx] Map into input grid of longitude indices of overlap cells
    real,dimension(:),allocatable::lat_in_grd ! [dgr] Interface latitudes
    real,dimension(:),allocatable::lon_in_grd ! [dgr] Interface longitudes
    real,dimension(:,:),allocatable::lon_in_grd_2d ! [dgr] Interface longitudes
    real,dimension(:,:),allocatable::area_in ! [m2] Area of gridcells
    real,dimension(:,:),allocatable::lat_in_2d ! [dgr] Latitude at gridcell center
    real,dimension(:,:),allocatable::lon_in_2d ! [dgr] Longitude at gridcell center
    real,dimension(:,:,:),allocatable::ovr_wgt ! [frc] Weight of overlapping input gridcells onto each output gridcell
    
    real,dimension(:,:),allocatable::hgt_sfc_in ! [m] Surface height 
    real,dimension(:,:),allocatable::lnd_frc_in ! [frc] Land fraction 
    real,dimension(:,:),allocatable::lnd_msk_in ! [msk] Land mask (0.0 or 1.0)
    real,dimension(:,:),allocatable::mbl_bsn_fct_in ! [frc] Mobilization enhancement due to basin characteristics
    real,dimension(:,:),allocatable::flw_acm_fct_in ! [frc] Flow accumulation factor
    real,dimension(:,:),allocatable::sfc_acm_fct_in ! [frc] Area accumulation factor
    real,dimension(:,:),allocatable::oro_in ! [frc] Orography 
    
    ! Main Code
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Entering "//sbr_nm
    
    ! Sanity check
    if(flg_NGDC.and.flg_CCM) stop "flg_NGDC and flg_CCM both true in hgt_sfc_get()"
    
    ! Check for monotonically increasing grids
    mnt_ncr=mnt_ncr_chk(lon_out_grd,lon_out_nbr+1)
    if (.not.mnt_ncr) stop "lon_out_grd not monotonically increasing in hgt_sfc_get()"
    mnt_ncr=mnt_ncr_chk(lat_out_grd,lat_out_nbr+1)
    if (.not.mnt_ncr) stop "lat_out_grd not monotonically increasing in hgt_sfc_get()"
    
    ! Initialize arrays
    bsn_enm_out(:,:)         = 0     ! [enm] Basin ID
    lnd_msk_out(:,:)         = 0     ! [flg] Land mask (integer 0 or 1)
    hgt_sfc_gpm_out(:,:)     = 0.0   ! [gpm] Surface geopotential height 
    hgt_sfc_out(:,:)         = 0.0   ! [m]   Surface height 
    hgt_sfc_std_dvn_out(:,:) = 0.0   ! [m]   Standard deviation of surface height 
    lnd_frc_out(:,:)         = 0.0   ! [frc] Land fraction 
    oro_out(:,:)             = 0.0   ! [frc] Orography 
    mbl_bsn_fct_out(:,:)     = 0.0   ! [frc] Mobilization enhancement due to basin characteristics
    flw_acm_fct_out(:,:)     = 0.0   ! [frc] Flow accumulation factor
    sfc_acm_fct_out(:,:)     = 0.0   ! [frc] Area accumulation factor
    
    ! Initialize scalars
    
    ! Read in netCDF data
    rcd=nf90_wrp_open(fl_in,nf90_nowrite,nc_id,sbr_nm=sbr_nm)
    ! Get dimension IDs
    rcd=nf90_wrp_inq_dimid(nc_id,"lat",lat_dmn_id)
    rcd=nf90_wrp_inq_dimid(nc_id,"lon",lon_dmn_id)
    ! Get dimension sizes
    rcd=rcd+nf90_inquire_dimension(nc_id,lat_dmn_id,len=lat_in_nbr)
    rcd=rcd+nf90_inquire_dimension(nc_id,lon_dmn_id,len=lon_in_nbr)
    if(rcd /= nf90_noerr) stop "Error retrieving dimension sizes"
    
    ! Enough memory? 
    if (lat_in_nbr > lat_nbr_max) stop "lat_in_nbr > lat_nbr_max in hgt_sfc_get()"
    if (lon_in_nbr > lon_nbr_max) stop "lon_in_nbr > lon_nbr_max in hgt_sfc_get()"
    
    ! Allocate based on input file dimensions
    allocate(area_in(lon_in_nbr,lat_in_nbr),stat=rcd) ! [m2] Area of gridcells
    if(rcd /= 0) stop "allocate() failed for area_in"
    allocate(lat_in_grd(lat_in_nbr+1),stat=rcd) ! [dgr] Interface latitudes
    if(rcd /= 0) stop "allocate() failed for lat_in_grd"
    allocate(lon_in_grd(lon_in_nbr+1),stat=rcd) ! [dgr] Interface longitudes
    if(rcd /= 0) stop "allocate() failed for lon_in_grd"
    allocate(lon_in_nbr_1d(lat_in_nbr),stat=rcd) ! [nbr] Number of longitudes per latitude
    if(rcd /= 0) stop "allocate() failed for lon_in_nbr_1d"
    allocate(lnd_frc_in(lon_in_nbr,lat_in_nbr),stat=rcd) ! [frc] Fraction of land (not ocean)
    if(rcd /= 0) stop "allocate() failed for lnd_frc_in"
    if(flg_NGDC.and.flg_read_2d_grd) then
       allocate(lat_in_2d(lon_in_nbr,lat_in_nbr),stat=rcd) ! [dgr] Latitude at gridcell center
       if(rcd /= 0) stop "allocate() failed for lat_in_2d"
       allocate(lon_in_grd_2d(lon_in_nbr+1,lat_in_nbr),stat=rcd) ! [dgr] Interface longitudes
       if(rcd /= 0) stop "allocate() failed for lon_in_grd_2d"
       allocate(lon_in_2d(lon_in_nbr,lat_in_nbr),stat=rcd) ! [dgr] Longitude at gridcell center
       if(rcd /= 0) stop "allocate() failed for lon_in_2d"
       allocate(lnd_msk_in(lon_in_nbr,lat_in_nbr),stat=rcd) ! [msk] Land mask (0.0 or 1.0)
       if(rcd /= 0) stop "allocate() failed for lnd_msk_in"
    endif ! endif flg_NGDC.and.flg_read_2d_grd
    
    allocate(hgt_sfc_in(lon_in_nbr,lat_in_nbr),stat=rcd) ! [m] Surface height 
    if(rcd /= 0) stop "allocate() failed for hgt_sfc_in"
    allocate(oro_in(lon_in_nbr,lat_in_nbr),stat=rcd) ! [frc] Orography
    if(rcd /= 0) stop "allocate() failed for oro_in"
    allocate(mbl_bsn_fct_in(lon_in_nbr,lat_in_nbr),stat=rcd) ! [frc] Mobilization enhancement due to basin characteristics
    if(rcd /= 0) stop "allocate() failed for mbl_bsn_fct_in"
    allocate(flw_acm_fct_in(lon_in_nbr,lat_in_nbr),stat=rcd) ! [frc] Flow accumulation factor
    if(rcd /= 0) stop "allocate() failed for flw_acm_fct_in"
    allocate(sfc_acm_fct_in(lon_in_nbr,lat_in_nbr),stat=rcd) ! [frc] Area accumulation factor
    if(rcd /= 0) stop "allocate() failed for sfc_acm_fct_in"
    
    ! Get variable IDs
    if(flg_NGDC) then
       hgt_sfc_nm="hgt_sfc" ! [sng] Surface height variable name
       lnd_frc_nm="LANDMASK" ! [sng] Fraction of land (not ocean) variable name
    endif ! endif flg_NGDC
    if(flg_CCM) then
       hgt_sfc_nm="htopo" ! [sng] Surface height variable name
       lnd_frc_nm="ftopo" ! [sng] Fraction of land (not ocean) variable name
    endif ! endif flg_CCM
    if(flg_NGDC.and.flg_read_2d_grd) then
       rcd=nf90_wrp_inq_varid(nc_id,"LATIXY",lat_in_2d_id)
       rcd=nf90_wrp_inq_varid(nc_id,"LONGXY",lon_in_2d_id)
       rcd=nf90_wrp_inq_varid(nc_id,"EDGEN",edg_nrt_id)
       rcd=nf90_wrp_inq_varid(nc_id,"EDGEE",edg_est_id)
       rcd=nf90_wrp_inq_varid(nc_id,"EDGES",edg_sth_id)
       rcd=nf90_wrp_inq_varid(nc_id,"EDGEW",edg_wst_id)
    endif ! endif flg_NGDC.and.flg_read_2d_grd
    if(flg_CCM) then
       rcd=nf90_wrp_inq_varid(nc_id,lnd_frc_nm,lnd_frc_id)
    endif ! endif flg_CCM
    rcd=nf90_wrp_inq_varid(nc_id,hgt_sfc_nm,hgt_sfc_id)
    ! Get number of dimensions
    rcd=rcd+nf90_inquire_variable(nc_id,hgt_sfc_id,ndims=hgt_sfc_dmn_nbr)
    if (hgt_sfc_dmn_nbr /= 2) then
       write (6,"(2a,i1,a)") prg_nm(1:ftn_strlen(prg_nm)), & 
            ": ERROR hgt_sfc has ",hgt_sfc_dmn_nbr," dimensions"
       stop
    endif                  ! endif err
    ! Get time-independent data
    if(flg_NGDC.and.flg_read_2d_grd) then
       rcd=nf90_wrp(nf90_get_var(nc_id,lat_in_2d_id,lat_in_2d),"get_var lat_in_2d")
       rcd=nf90_wrp(nf90_get_var(nc_id,lon_in_2d_id,lon_in_2d),"get_var lon_in_2d")
       rcd=nf90_wrp(nf90_get_var(nc_id,edg_nrt_id,edg_grd(1)),"get_var edg_nrt")
       rcd=nf90_wrp(nf90_get_var(nc_id,edg_est_id,edg_grd(2)),"get_var edg_est")
       rcd=nf90_wrp(nf90_get_var(nc_id,edg_sth_id,edg_grd(3)),"get_var edg_sth")
       rcd=nf90_wrp(nf90_get_var(nc_id,edg_wst_id,edg_grd(4)),"get_var edg_wst")
    endif ! endif flg_NGDC.and.flg_read_2d_grd
    if(flg_CCM) then
       rcd=nf90_wrp(nf90_get_var(nc_id,lnd_frc_id,lnd_frc_in),"get_var lnd_frc_in")
    endif ! endif flg_CCM
    rcd=nf90_wrp(nf90_get_var(nc_id,hgt_sfc_id,hgt_sfc_in),"get_var hgt_sfc_in")
    if (rcd /= nf90_noerr) stop "Error in nf90_get_var in hgt_sfc_get()"
    ! Get missing value
    rcd=rcd+nf90_inquire_attribute(nc_id,hgt_sfc_id,"missing_value",attnum=mss_val_id)
    if (rcd==nf90_noerr) then
       mss_flg=.true. ! [flg] Data may have missing values
       rcd=rcd+nf90_get_att(nc_id,hgt_sfc_id,"missing_value",mss_val_in)
       if (real(mss_val_in) /= real(1.0e36)) write (6,"(2a,es15.8)") prg_nm(1:ftn_strlen(prg_nm)), &
            ": WARNING "//sbr_nm//" reports mss_val_in = ",mss_val_in
    else
       mss_flg=.false. ! [flg] Data may have missing values
       rcd=nf90_noerr ! [enm] Return success code
    endif ! endif
    
    if(rcd /= nf90_noerr) stop "Error retrieving variable data"
    ! Close file
    rcd=nf90_wrp_close(nc_id,fl_in,"Ingested",sbr_nm=sbr_nm) ! [fnc] Close file
    
    ! Convert input data to SI
    if(hgt_dlt_msl /= 0.0) then
       ! Mean sea level decrease (due to, e.g., glaciation) increases elevations above sea level relative to present day
       ! Mean sea level increase (due to, e.g., warming) decreases elevations above sea level relative to present day
       hgt_sfc_in(:,:)=hgt_sfc_in(:,:)-hgt_dlt_msl ! [m] Surface height
       write (6,"(a,f9.4,a)") "Adding ",hgt_dlt_msl," m to mean sea level"
       if(flg_CCM) then
          write (6,"(2a)") prg_nm(1:ftn_strlen(prg_nm)), &
               ": WARNING "//sbr_nm//" is not accounting for MSL change in CCM landmask ftopo"
       endif ! endif flg_CCM
    endif ! endif SLC
    
    ! Determine map grid and overlap characteristics
    ! fxm: Method breaks for reduced grids, streamline, subroutinize, and generalize
    lon_in_nbr_1d(:)=lon_in_nbr ! [nbr] Number of longitudes per latitude
    
    ! Convert gridpoint centers, domain boundaries to grid interfaces
    if(flg_NGDC.and.flg_read_2d_grd) then
       call map_edge_mk(lat_in_nbr,lon_in_nbr,lon_in_nbr_1d, & ! I
            lat_in_2d,lon_in_2d,edg_grd, & ! I 
            lat_in_grd,lon_in_grd,lon_in_grd_2d) ! O
    else
       call map_grd_mk(lat_in_nbr,lon_in_nbr,map_typ_in, & ! I
            lat_in_grd,lon_in_grd) ! O
    endif ! endif flg_NGDC.and.flg_read_2d_grd
    
    mnt_ncr=mnt_ncr_chk(lat_in_grd,lat_in_nbr+1)
    if (.not.mnt_ncr) stop "ERROR: lat_in_grd not monotonically increasing in hgt_sfc_get()"
    mnt_ncr=mnt_ncr_chk(lon_in_grd,lon_in_nbr+1)
    if (.not.mnt_ncr) stop "ERROR: lon_in_grd not monotonically increasing in hgt_sfc_get()"
    
    ! Diagnostic area on input grid
    call map_area_get(lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
         area_in) ! O
    
    ! Determine space required by overlap arrays
    call map_ovr_nbr_max_get(lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
         lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
         ovr_nbr_max) ! O
    
    ! Allocate arrays that depend on ovr_nbr_max
    allocate(ovr_lat_idx(lon_out_nbr,lat_out_nbr,ovr_nbr_max),stat=rcd) ! [idx] Map into input grid of latitude indices of overlap cells
    if(rcd /= 0) stop "allocate() failed for ovr_lat_idx"
    allocate(ovr_lon_idx(lon_out_nbr,lat_out_nbr,ovr_nbr_max),stat=rcd) ! [idx] Map into input grid of longitude indices of overlap cells
    if(rcd /= 0) stop "allocate() failed for ovr_lon_idx"
    allocate(ovr_wgt(lon_out_nbr,lat_out_nbr,ovr_nbr_max),stat=rcd) ! [frc] Weight of overlapping input gridcells onto each output gridcell
    if(rcd /= 0) stop "allocate() failed for ovr_wgt"
    
    ! Derive land fraction for this subroutine
    if(flg_NGDC) then
       ! fxm: Inappropriate for basins like Mojave, Dead Sea, lower than 0 m MSL
       ! Mojave is significant dust source so this problem must be addressed
       where (hgt_sfc_in > 0.0) ! [m]
          lnd_frc_in=1.0 ! [frc] Land fraction
       elsewhere
          lnd_frc_in=0.0 ! [frc] Land fraction
       end where ! end where hgt_sfc < 0.0
    endif ! endif flg_NGDC
    
    ! Land mask is unity for LSM/CLM grid cells and 0 for ocean points
    ! Use land mask in array as weights compute land fraction on output grid
    ! For this purpose, temporarily set input land mask to unity everywhere
    if(flg_NGDC .and. flg_read_2d_grd) then
       lnd_msk_in(:,:)=1.0 ! [msk] Land mask (0.0 or 1.0)
    endif ! endif flg_NGDC.and.flg_read_2d_grd
    
    ! Archive input in human-readable format if desired
    if(.false.) then
       ! Write surface height data in text format of four columns: lat, lon, lnd_frc, and elevation
       open (unit=unit_dgn,file=fl_dgn,status="unknown",iostat=rcd)
       if (rcd /= 0) write (6,"(4a,i4)") prg_nm(1:ftn_strlen(prg_nm)), & 
            ": ERROR "//sbr_nm//" unable to open ",fl_dgn(1:ftn_strlen(fl_dgn))
       do lat_in_idx=1,lat_in_nbr ! NB: Outer loop over lat
          do lon_in_idx=1,lon_in_nbr
             write (unit_dgn,"(4(f10.4))") &
                  0.5*(lat_in_grd(lat_in_idx)+lat_in_grd(lat_in_idx+1)), & ! [dgr]
                  0.5*(lon_in_grd(lon_in_idx)+lon_in_grd(lon_in_idx+1)), & ! [dgr]
                  lnd_frc_in(lon_in_idx,lat_in_idx), & ! [frc]
                  hgt_sfc_in(lon_in_idx,lat_in_idx) ! [m]
          end do ! end loop over lon
       end do ! end loop over lat
       close (unit_dgn)
       write (6,"(a,1x,a)") "Wrote surface height data in text format to",fl_dgn(1:ftn_strlen(fl_dgn))
    endif ! endif dbg
    
    ! Sanity check
    err_nbr=0
    do lat_in_idx=1,lat_in_nbr
       do lon_in_idx=1,lon_in_nbr
          if (lnd_frc_in(lon_in_idx,lat_in_idx)==0.0.and.hgt_sfc_in(lon_in_idx,lat_in_idx) > 0.0) then
             ! According to CSM tools/getopo.F90, Caspian Sea is only region defined
             ! as ocean which should have non-zero surface height (because it is landlocked)
             if ( lon_in_grd(lon_in_idx) < 125.0.and. &
                  lon_in_grd(lon_in_idx) > 135.0.and. &
                  lat_in_grd(lat_in_idx) < 30.0.and. &
                  lat_in_grd(lat_in_idx) > 50.0) then
                err_nbr=err_nbr+1
                lon_in_idx_err=lon_in_idx
                lat_in_idx_err=lat_in_idx
             end if ! endif not Caspian Sea
          end if ! endif err
       end do ! end loop over lon
    end do ! end loop over lat
    if (err_nbr > 0) then
       frc_ttl=hgt_sfc_in(lon_in_idx_err,lat_in_idx_err)
       write (6,"(a,i3,a)") "ERROR hgt_sfc_get() reports ",err_nbr," surface height errors during input"
       write (6,"(a)") "Surface height properties:"
       write (6,"(2(a,7x))") "Land fraction","Surface height"
       write (6,"(2(f12.6,1x))") &
            lnd_frc_in(lon_in_idx_err,lat_in_idx_err), &
            hgt_sfc_in(lon_in_idx_err,lat_in_idx_err)
       write (6,"(a)") "Input cell edge locations:"
       write (6,"(a)") "(idx)   Lat Sth (idx)   Lat Nrt (idx)   Lon Wst (idx)   Lon Est"
       write (6,"(4(a1,i3,a1,1x,f9.4,1x))") &
            "(",lat_in_idx_err,")",lat_in_grd(lat_in_idx_err),"(",lat_in_idx_err+1,")",lat_in_grd(lat_in_idx_err+1), &
            "(",lon_in_idx_err,")",lon_in_grd(lon_in_idx_err),"(",lon_in_idx_err+1,")",lon_in_grd(lon_in_idx_err+1)
       stop
    end if                    ! endif err
    
    ! Get overlap locations and weights for mapping input grid to output grid
    call map_ovr_wgt_drv( &
         lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
         lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
         ovr_nbr_max,area_out, & ! I
         ovr_lat_idx,ovr_lon_idx,ovr_nbr,ovr_wgt) ! O
    
    ! Rebin land fraction from input grid to output grid
    call map_rbn(lnd_frc_in,mss_flg,mss_val_in, & ! I
         lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
         lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
         ovr_nbr_max,ovr_lat_idx,ovr_lon_idx,ovr_nbr,ovr_wgt, & ! I
         lnd_frc_out) ! O
    
    ! Rebin surface height from input grid to output grid
    call map_rbn(hgt_sfc_in,mss_flg,mss_val_in, & ! I
         lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
         lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
         ovr_nbr_max,ovr_lat_idx,ovr_lon_idx,ovr_nbr,ovr_wgt, & ! I
         hgt_sfc_out) ! O
    
    ! Find surface height variance on output grid
    ! fxm: Add 1-2-1 filter and 2x2 and 3x3 coarse cases
    call map_rbn_var(hgt_sfc_in,mss_flg,mss_val_in, & ! I
         lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
         lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
         ovr_nbr_max,ovr_lat_idx,ovr_lon_idx,ovr_nbr,ovr_wgt, & ! I
         hgt_sfc_out, & ! I
         hgt_sfc_var_out) ! O
    
    ! Process surface height data in accordance with CSM algorithms
    ! Following code from CSM tools/fmain.F90, tools/sghphis.F90
    where (hgt_sfc_var_out < 0.5) ! [m]
       hgt_sfc_std_dvn_out=0.0 ! [m]
    elsewhere
       ! fxm: 19991029 Unsure whether CCM SGH is supposed to be deviation of geometric height or of geopotential height
       hgt_sfc_std_dvn_out=sqrt(hgt_sfc_var_out) ! [m]
    end where                 ! end where hgt_sfc_var_out < 0.5
    hgt_sfc_gpm_out(:,:)=hgt_sfc_out(:,:)*9.80616 ! [gpm]
    
    ! Define orography on input and output grids
    where (lnd_frc_in >= 0.5) ! [frc]
       oro_in=1.0               ! [frc]
    elsewhere
       oro_in=0.0            ! [frc]
    end where                 ! end where lnd_frc_out >= 0.5
    
    ! Land mask (lnd_msk) determines whether point has (partial) fluxes computed by CLM/LSM
    ! Offline, CLM/LSM consider a gridcell as land only if more than 50% is land
    ! Coupled, CLM/LSM consider any gridcell with any land to be land (and then
    ! apply a land fraction correction to fluxes), as well as any points that have
    ! no land, but which the ocean model does not want to treat as having any ocean
    ! Coupled, CLM enforces conditions that south pole is land, north pole is not
    ! clm2/src/mksrfdata/mkgridMod.F90
    where (lnd_frc_out >= 0.5) ! [frc]
       lnd_msk_out=1 ! [frc]
       oro_out=1.0 ! [frc]
    elsewhere
       lnd_msk_out=0 ! [frc]
       oro_out=0.0 ! [frc]
       hgt_sfc_std_dvn_out=0.0 ! [m]
    end where ! end where lnd_frc_out >= 0.5
    
    if(bsn_fct_hrz) then
       ! Compute source efficiency factor from high-resolution topography then rebin basin factor to model grid (still experimental)
       write (6,"(a)") prg_nm(1:ftn_strlen(prg_nm)),"Computing basin factor at high resolution"
       call rdb_fct_tpg_GCT01(lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
            hgt_sfc_in,oro_in, & ! I
            mbl_bsn_fct_in) ! O
       
       ! Rebin basin factor from input grid to output grid
       call map_rbn(mbl_bsn_fct_in,mss_flg,mss_val_in, & ! I
            lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
            lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
            ovr_nbr_max,ovr_lat_idx,ovr_lon_idx,ovr_nbr,ovr_wgt, & ! I
            mbl_bsn_fct_out) ! O
       
       ! djn
       ! fxm: high-resolution calculation of flw_acm_fct and sfc_acm_fct not yet coded
       ! Need to call dst_src_flw here
       ! For now, set:
       !   flw_acm_fct_in  = 1.0
       !   sfc_acm_fct_in = 1.0
       print *, "WARNING: high-resolution calculation of flw_acm_fct and sfc_acm_fct not yet coded"
       print *, "WARNING: setting flw_acm_fct  = 1.0"
       print *, "WARNING: setting sfc_acm_fct = 1.0"
       flw_acm_fct_in(:,:)=1.0
       sfc_acm_fct_in(:,:)=1.0
       
       ! Rebin basin factor from input grid to output grid
       call map_rbn(flw_acm_fct_in,mss_flg,mss_val_in, & ! I
            lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
            lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
            ovr_nbr_max,ovr_lat_idx,ovr_lon_idx,ovr_nbr,ovr_wgt, & ! I
            flw_acm_fct_out) ! O
       
       call map_rbn(sfc_acm_fct_in,mss_flg,mss_val_in, & ! I
            lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
            lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
            ovr_nbr_max,ovr_lat_idx,ovr_lon_idx,ovr_nbr,ovr_wgt, & ! I
            sfc_acm_fct_out) ! O
       
    endif ! endif .false
    
    ! De-allocate
    if (allocated(area_in)) deallocate(area_in,stat=rcd)
    if(rcd /= 0) stop "deallocate() failed for area_in"
    if (allocated(lat_in_grd)) deallocate(lat_in_grd,stat=rcd)
    if(rcd /= 0) stop "deallocate() failed for lat_in_grd"
    if (allocated(lat_in_2d)) deallocate(lat_in_2d,stat=rcd)
    if(rcd /= 0) stop "deallocate() failed for lat_in_2d"
    if (allocated(lon_in_grd_2d)) deallocate(lon_in_grd_2d,stat=rcd)
    if(rcd /= 0) stop "deallocate() failed for lon_in_grd_2d"
    if (allocated(lon_in_2d)) deallocate(lon_in_2d,stat=rcd)
    if(rcd /= 0) stop "deallocate() failed for lon_in_2d"
    if (allocated(lon_in_grd)) deallocate(lon_in_grd,stat=rcd)
    if(rcd /= 0) stop "deallocate() failed for lon_in_grd"
    if (allocated(lon_in_nbr_1d)) deallocate(lon_in_nbr_1d,stat=rcd)
    if(rcd /= 0) stop "deallocate() failed for lon_in_nbr_1d"
    if (allocated(ovr_lat_idx)) deallocate(ovr_lat_idx,stat=rcd)
    if(rcd /= 0) stop "deallocate() failed for ovr_lat_idx"
    if (allocated(ovr_lon_idx)) deallocate(ovr_lon_idx,stat=rcd)
    if(rcd /= 0) stop "deallocate() failed for ovr_lon_idx"
    if (allocated(ovr_wgt)) deallocate(ovr_wgt,stat=rcd)
    if(rcd /= 0) stop "deallocate() failed for ovr_wgt"
    if (allocated(lnd_frc_in)) deallocate(lnd_frc_in,stat=rcd)
    if(rcd /= 0) stop "deallocate() failed for lnd_frc_in"
    if (allocated(lnd_msk_in)) deallocate(lnd_msk_in,stat=rcd)
    if(rcd /= 0) stop "deallocate() failed for lnd_msk_in"
    
    if (allocated(hgt_sfc_in)) deallocate(hgt_sfc_in,stat=rcd)
    if(rcd /= 0) stop "deallocate() failed for hgt_sfc_in"
    if (allocated(mbl_bsn_fct_in)) deallocate(mbl_bsn_fct_in,stat=rcd)
    if(rcd /= 0) stop "deallocate() failed for mbl_bsn_fct_in"
    if (allocated(flw_acm_fct_in)) deallocate(flw_acm_fct_in,stat=rcd)
    if(rcd /= 0) stop "deallocate() failed for flw_acm_fct_in"
    if (allocated(sfc_acm_fct_in)) deallocate(sfc_acm_fct_in,stat=rcd)
    if(rcd /= 0) stop "deallocate() failed for sfc_acm_fct_in"
    if (allocated(oro_in)) deallocate(oro_in,stat=rcd)
    if(rcd /= 0) stop "deallocate() failed for oro_in"
    
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Exiting "//sbr_nm
    return 
  end subroutine hgt_sfc_get ! end hgt_sfc_get()
  
  subroutine hgt_sfc_get_old( & ! 
       fl_in, & ! I
       lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
       lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
       ovr_nbr_max,area_out, & ! I
       hgt_sfc_gpm_out,hgt_sfc_out,hgt_sfc_std_dvn_out, & ! O
       lnd_frc_out,lnd_msk_out,oro_out, & ! O
       mbl_bsn_fct_out) ! O
    ! Purpose: Retrieve surface height and rebin to requested grid
    ! hgt_sfc_get_old() is called by bds_prc()
    ! hgt_sfc_get_old() has been deprecated in favor of hgt_sfc_get()
    ! hgt_sfc_get_old() does not use dynamic memory techniques and cannot use bathymetry
    use bds_ctl,only:bsn_fct_hrz ! [mdl] Control variables drc_in,drc_out,hgt_dlt_msl
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use map_grd ! [mdl] Map grids and regridding
    use netcdf ! [mdl] netCDF interface
    use nf90_utl ! [mdl] netCDF utilities
    use sng_mdl,only:ftn_strlen,ftn_strstr,ftn_strnul ! [mdl] String manipulation
    use src_prc,only:rdb_fct_tpg_GCT01 ! [mdl] Source processing, identification
    use utl_mdl,only:mnt_ncr_chk ! [mdl] Utility functions (date_time_get,mnt_chk...)
    implicit none
    ! Parameters
    character(len=*),parameter::sbr_nm="hgt_sfc_get_old" ! [sng] Subroutine name
    integer,parameter::fl_in_unit=73 ! Unit for reading input data
    integer,parameter::unit_dgn=74 ! Unit for writing general Map diagnoses
    ! Commons
    ! Input
    character(len=*),intent(in)::fl_in ! [sng] Input file
    integer,intent(in)::lat_in_nbr ! [nbr] Number of latitudes
    integer,intent(in)::lat_out_nbr ! [nbr] Number of latitudes
    integer,intent(in)::lon_in_nbr ! [nbr] Number of longitudes
    integer,intent(in)::lon_out_nbr ! [nbr] Number of longitudes
    integer,intent(in)::ovr_nbr_max ! [nbr] Maximum number of input cells which overlap any output cell
    real,intent(in)::area_out(lon_out_nbr,lat_out_nbr) ! [m2] Area of gridcells
    real,intent(in)::lat_in_grd(lat_in_nbr+1) ! [dgr] Interface latitudes
    real,intent(in)::lat_out_grd(lat_out_nbr+1) ! [dgr] Interface latitudes
    real,intent(in)::lon_in_grd(lon_in_nbr+1) ! [dgr] Interface longitudes
    real,intent(in)::lon_out_grd(lon_out_nbr+1) ! [dgr] Interface longitudes
    ! Output
    integer,intent(out)::lnd_msk_out(lon_out_nbr,lat_out_nbr) ! [flg] Land mask (integer 0 or 1)
    real,intent(out)::hgt_sfc_gpm_out(lon_out_nbr,lat_out_nbr) ! [gpm] Surface geopotential height 
    real,intent(out)::hgt_sfc_out(lon_out_nbr,lat_out_nbr) ! [m] Surface height 
    real,intent(out)::hgt_sfc_std_dvn_out(lon_out_nbr,lat_out_nbr) ! [m] Standard deviation of surface height 
    real,intent(out)::lnd_frc_out(lon_out_nbr,lat_out_nbr) ! [frc] Land fraction 
    real,intent(out)::oro_out(lon_out_nbr,lat_out_nbr) ! [frc] Orography 
    real,intent(out)::mbl_bsn_fct_out(lon_out_nbr,lat_out_nbr) ! [frc] Mobilization enhancement due to basin characteristics
    ! Local
    character fl_dgn_hgt_sfc*80 ! Diagnostics file
    integer err_nbr           ! [nbr] Number of errors processing surface types
    integer idx               ! [idx] Counting index
    integer lat_in_idx        ! [idx] Counting index for lat
    integer lat_in_idx_err    ! Latitude of gridcell which caused error
    integer lat_out_idx       ! [idx] Counting index for lat
    integer lon_in_idx        ! [idx] Counting index for lon
    integer lon_in_idx_err    ! Longitude of gridcell which caused error
    integer lon_out_idx       ! [idx] Counting index for lon
    integer ovr_idx           ! [idx] Counting index
    integer ovr_lat_idx(lon_out_nbr,lat_out_nbr,ovr_nbr_max) ! [idx] Map into input grid of latitude indices of overlap cells
    integer ovr_lon_idx(lon_out_nbr,lat_out_nbr,ovr_nbr_max) ! [idx] Map into input grid of longitude indices of overlap cells
    integer ovr_nbr(lon_out_nbr,lat_out_nbr) ! [nbr] Number of input gridcells which overlap each output gridcell
    integer ovr_nbr_crr       ! [nbr] Current number of overlapping gridcells
    integer rcd               ! [enm] Return success code
    logical mnt_ncr           ! Monotonic and increasing flag
    real area_in(lon_in_nbr,lat_in_nbr) ! [m2] Area of gridcells
    real eps_rlt              ! [frc] Relative error allowed in frc_ttl
    real frc_ttl              ! [frc] Total fraction of gridcell accounted for by all soil textures
    real oro_in(lon_in_nbr,lat_in_nbr) ! [frc] Orography 
    real ovr_wgt(lon_out_nbr,lat_out_nbr,ovr_nbr_max) ! [frc] Weight of overlapping input gridcells onto each output gridcell
    real ovr_wgt_crr          ! [frc] Overlap weight of current gridcell
    real hgt_sfc_var(lon_out_nbr,lat_out_nbr) ! [m2] Variance of surface height
    ! Variables needed for CSM_TOPO dataset
    integer lat_dmn_id        ! [id] Dimension ID for lat
    integer lat_in_id         ! [id] Variable ID
    integer lat_nbr           ! [nbr] Number of latitudes in input file
    integer lon_dmn_id        ! [id] Dimension ID for lon
    integer lon_nbr           ! [nbr] Number of longitudes in input file
    integer mss_val_id        ! [id] Attribute ID
    integer nc_id             ! [id] File handle
    integer hgt_sfc_dmn_nbr   ! [nbr] Number of dimensions in disk hgt_sfc field
    integer lnd_frc_id        ! [id] Variable ID
    integer hgt_sfc_id        ! [id] Variable ID
    logical flg_CSM_TOPO      ! Input file is CSM hgt_sfc_topo.nc
    logical lat_mnt_ncr       ! Input latitude monotonically increases
    logical mss_flg           ! Variable has missing_value attribute
    real lat_in(lat_in_nbr)   ! Midpoint latitudes
    real mss_val              ! Missing value
    real lnd_frc_in(lon_in_nbr,lat_in_nbr) ! [frc] Land fraction 
    real hgt_sfc_in(lon_in_nbr,lat_in_nbr) ! [m] Surface height 
    real mbl_bsn_fct_in(lon_in_nbr,lat_in_nbr) ! [frc] Mobilization enhancement due to basin characteristics
    ! Debugging
    integer lon_idx_dbg       ! Longitude for verbose output
    integer lat_idx_dbg       ! Longitude for verbose output
    
    ! Main Code
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Entering hgt_sfc_get_old()"
    
    ! Check for monotonically increasing grids
    mnt_ncr=mnt_ncr_chk(lon_out_grd,lon_out_nbr+1)
    if (.not.mnt_ncr) stop "lon_out_grd not monotonically increasing in hgt_sfc_get_old()"
    mnt_ncr=mnt_ncr_chk(lat_out_grd,lat_out_nbr+1)
    if (.not.mnt_ncr) stop "lat_out_grd not monotonically increasing in hgt_sfc_get_old()"
    
    ! Initialize scalars
    eps_rlt=1.0e-5            ! [frc] Relative error allowed in frc_ttl
    flg_CSM_TOPO=.true.
    lon_idx_dbg=90
    lat_idx_dbg=50
    mss_val=1.0e36
    
    ! Initialize output fields
    oro_out=0.0               ! [frc]
    hgt_sfc_gpm_out=0.0       ! [gpm]
    hgt_sfc_out=0.0           ! [m]
    hgt_sfc_std_dvn_out=0.0   ! [m]
    mbl_bsn_fct_out=0.0       ! [frc] Mobilization enhancement due to basin characteristics
    
    ! Compute gridcell area 
    call map_area_get(lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
         area_in) ! O
    
    fl_dgn_hgt_sfc="/tmp/zender/map/hgt_sfc.txt" ! Diagnostics file
    call ftn_strnul(fl_dgn_hgt_sfc)
    
    if (flg_CSM_TOPO) then
       ! NCAR SCD DSS754.0 is U.S. Navy Global Elevation 10-arcminute resolution dataset described at http://www.scd.ucar.edu/dss/catalogs/geo.html
       ! CSM topography file hgt_sfc_topo.nc was derived from DSS754.0
       ! This is the same as LSM elevation file mksrf_elev.nc
       ! Subroutine tools/getopo.F90 defines CSM method used
       ! CSM initial fields SGH, PHIS, and ORO are defined using mean and variance of height and fractional cover fields of this file
       ! Dataset resolution is 1/6 x 1/6 degree and data are written from South pole to North pole, East to West, beginning at Greenwich
       ! Read netCDF data
       rcd=nf90_noerr           ! nf90_noerr == 0
       rcd=nf90_wrp_open(fl_in,NF90_NOWRITE,nc_id,sbr_nm=sbr_nm)
       ! Get dimension IDs
       rcd=nf90_wrp_inq_dimid(nc_id,"lon",lon_dmn_id)
       rcd=nf90_wrp_inq_dimid(nc_id,"lat",lat_dmn_id)
       if (rcd /= nf90_noerr) then
          write (6,"(3a)") prg_nm(1:ftn_strlen(prg_nm)), & 
               ": ERROR Unable to find coordinate in ",fl_in(1:ftn_strlen(fl_in))
          rcd=nf90_noerr
       endif                  ! endif
       ! Get dimension sizes
       rcd=rcd+nf90_inquire_dimension(nc_id,lon_dmn_id,len=lon_nbr)
       rcd=rcd+nf90_inquire_dimension(nc_id,lat_dmn_id,len=lat_nbr)
       ! Enough memory? 
       if (lat_nbr /= lat_in_nbr) stop "lat_nbr /= lat_in_nbr in hgt_sfc_get_old()"
       if (lon_nbr /= lon_in_nbr) stop "lon_nbr /= lon_in_nbr in hgt_sfc_get_old()"
       ! Get variable IDs
       rcd=nf90_wrp_inq_varid(nc_id,"ftopo",lnd_frc_id)
       rcd=nf90_wrp_inq_varid(nc_id,"htopo",hgt_sfc_id)
       rcd=nf90_wrp_inq_varid(nc_id,"lat",lat_in_id)
       ! Get number of dimensions
       rcd=rcd+nf90_inquire_variable(nc_id,hgt_sfc_id,ndims=hgt_sfc_dmn_nbr)
       if (hgt_sfc_dmn_nbr /= 2) then
          write (6,"(2a,i1,a)") prg_nm(1:ftn_strlen(prg_nm)), & 
               ": ERROR hgt_sfc_get_old() reports htopo has ",hgt_sfc_dmn_nbr," dimensions"
          stop
       endif                  ! endif err
       rcd=rcd+nf90_inquire_attribute(nc_id,hgt_sfc_id,"missing_value",attnum=mss_val_id)
       if (rcd==nf90_noerr) then
          mss_flg=.true.
          write (6,"(3a)") prg_nm(1:ftn_strlen(prg_nm)), & 
               ": WARNING hgt_sfc_get_old() reports hgt_sfc has missing_value in ", &
               fl_in(1:ftn_strlen(fl_in))
          rcd=rcd+nf90_get_att(nc_id,hgt_sfc_id,"missing_value",mss_val)
       else
          mss_flg=.false.
          rcd=nf90_noerr
       endif                  ! endif
       ! Get data
       rcd=nf90_wrp(nf90_get_var(nc_id,lnd_frc_id,lnd_frc_in),"get_var lnd_frc_in")
       rcd=nf90_wrp(nf90_get_var(nc_id,hgt_sfc_id,hgt_sfc_in),"get_var hgt_sfc_in")
       rcd=nf90_wrp(nf90_get_var(nc_id,lat_in_id,lat_in),"get_var lat_in")
       rcd=nf90_wrp_close(nc_id,fl_in,'',sbr_nm=sbr_nm)
       
       lat_mnt_ncr=mnt_ncr_chk(lat_in,lat_in_nbr)
       if (.not.lat_mnt_ncr) stop "ERROR: CSM_TOPO lat_in decreases in hgt_sfc_get_old()"
    endif                     ! endif CSM_TOPO surface height
    write (6,"(a,1x,a)") "Read surface height data from",fl_in(1:ftn_strlen(fl_in))
    
    if (dbg_lvl==dbg_crr) then
       ! Write surface height data in text format of four columns: lat, lon, lnd_frc, and elevation
       open (unit_dgn,file=fl_dgn_hgt_sfc,status="unknown",iostat=rcd)
       if (rcd /= 0) write (6,"(4a,i4)") prg_nm(1:ftn_strlen(prg_nm)), & 
            ": ERROR hgt_sfc_get_old() unable to open ",fl_dgn_hgt_sfc(1:ftn_strlen(fl_dgn_hgt_sfc))
       do lat_in_idx=1,lat_in_nbr ! NB: Outer loop over lat
          do lon_in_idx=1,lon_in_nbr
             write (unit_dgn,"(4(f10.4))") &
                  0.5*(lat_in_grd(lat_in_idx)+lat_in_grd(lat_in_idx+1)), & ! [dgr]
                  0.5*(lon_in_grd(lon_in_idx)+lon_in_grd(lon_in_idx+1)), & ! [dgr]
                  lnd_frc_in(lon_in_idx,lat_in_idx), & ! [frc]
                  hgt_sfc_in(lon_in_idx,lat_in_idx) ! [m]
          end do ! end loop over lon
       end do ! end loop over lat
       close (unit_dgn)
       write (6,"(a,1x,a)") "Wrote surface height in text format to",fl_dgn_hgt_sfc(1:ftn_strlen(fl_dgn_hgt_sfc))
    endif                     ! endif dbg
    
    ! Sanity check
    err_nbr=0
    do lat_in_idx=1,lat_in_nbr
       do lon_in_idx=1,lon_in_nbr
          if (lnd_frc_in(lon_in_idx,lat_in_idx)==0.0.and.hgt_sfc_in(lon_in_idx,lat_in_idx) > 0.0) then
             ! According to CSM tools/getopo.F90, Caspian Sea is only region defined
             ! as ocean which should have non-zero surface height (because it is landlocked)
             if ( lon_in_grd(lon_in_idx) < 125.0.and. &
                  lon_in_grd(lon_in_idx) > 135.0.and. &
                  lat_in_grd(lat_in_idx) < 30.0.and. &
                  lat_in_grd(lat_in_idx) > 50.0) then
                err_nbr=err_nbr+1
                lon_in_idx_err=lon_in_idx
                lat_in_idx_err=lat_in_idx
             end if           ! endif not Caspian Sea
          end if              ! endif err
       end do ! end loop over lon
    end do ! end loop over lat
    if (err_nbr > 0) then
       frc_ttl=hgt_sfc_in(lon_in_idx_err,lat_in_idx_err)
       write (6,"(a,i3,a)") "ERROR hgt_sfc_get_old() reports ",err_nbr," surface height errors during input"
       write (6,"(a)") "Surface height properties:"
       write (6,"(2(a,7x))") "Land fraction","Surface height"
       write (6,"(2(f12.6,1x))") &
            lnd_frc_in(lon_in_idx_err,lat_in_idx_err), &
            hgt_sfc_in(lon_in_idx_err,lat_in_idx_err)
       write (6,"(a)") "Input cell edge locations:"
       write (6,"(a)") "(idx)   Lat Sth (idx)   Lat Nrt (idx)   Lon Wst (idx)   Lon Est"
       write (6,"(4(a1,i3,a1,1x,f9.4,1x))") &
            "(",lat_in_idx_err,")",lat_in_grd(lat_in_idx_err),"(",lat_in_idx_err+1,")",lat_in_grd(lat_in_idx_err+1), &
            "(",lon_in_idx_err,")",lon_in_grd(lon_in_idx_err),"(",lon_in_idx_err+1,")",lon_in_grd(lon_in_idx_err+1)
       stop
    end if                    ! endif err
    
    ! Get overlap locations and weights for mapping input grid to output grid
    call map_ovr_wgt_drv( &
         lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
         lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
         ovr_nbr_max,area_out, & ! I
         ovr_lat_idx,ovr_lon_idx,ovr_nbr,ovr_wgt) ! O
    
    ! Rebin land fraction from input grid to output grid
    call map_rbn(lnd_frc_in,mss_flg,mss_val, & ! I
         lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
         lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
         ovr_nbr_max,ovr_lat_idx,ovr_lon_idx,ovr_nbr,ovr_wgt, & ! I
         lnd_frc_out) ! O
    
    ! Rebin surface height from input grid to output grid
    call map_rbn(hgt_sfc_in,mss_flg,mss_val, & ! I
         lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
         lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
         ovr_nbr_max,ovr_lat_idx,ovr_lon_idx,ovr_nbr,ovr_wgt, & ! I
         hgt_sfc_out) ! O
    
    ! Find surface height variance on output grid
    ! fxm: Add 1-2-1 filter and 2x2 and 3x3 coarse cases
    call map_rbn_var(hgt_sfc_in,mss_flg,mss_val, & ! I
         lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
         lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
         ovr_nbr_max,ovr_lat_idx,ovr_lon_idx,ovr_nbr,ovr_wgt, & ! I
         hgt_sfc_out, & ! I
         hgt_sfc_var) ! O
    
    ! Process surface height data in accordance with CSM algorithms
    ! Following code from CSM tools/fmain.F90, tools/sghphis.F90
    where (hgt_sfc_var < 0.5) ! [m]
       hgt_sfc_std_dvn_out=0.0   ! [m]
    elsewhere
       ! fxm: 19991029 Unsure whether CCM SGH is supposed to be deviation of regular height or of geopotential height
       hgt_sfc_std_dvn_out=sqrt(hgt_sfc_var) ! [m]
    end where                 ! end where hgt_sfc_var < 0.5
    hgt_sfc_gpm_out=hgt_sfc_out*9.80616 ! [gpm]
    
    ! Define orography on input and output grids
    where (lnd_frc_in >= 0.5) ! [frc]
       oro_in=1.0               ! [frc]
    elsewhere
       oro_in=0.0            ! [frc]
    end where                 ! end where lnd_frc_out >= 0.5
    
    ! Land mask (lnd_msk) determines whether point has (partial) fluxes computed by CLM/LSM
    ! Offline, CLM/LSM consider a gridcell as land only if more than 50% is land
    ! Coupled, CLM/LSM consider any gridcell with any land to be land (and then
    ! apply a land fraction correction to fluxes), as well as any points that have
    ! no land, but which the ocean model does not want to treat as having any ocean
    ! Coupled, CLM enforces conditions that south pole is land, north pole is not
    ! clm2/src/mksrfdata/mkgridMod.F90
    where (lnd_frc_out >= 0.5) ! [frc]
       lnd_msk_out=1 ! [frc]
       oro_out=1.0 ! [frc]
    elsewhere
       lnd_msk_out=0 ! [frc]
       oro_out=0.0 ! [frc]
       hgt_sfc_std_dvn_out=0.0 ! [m]
    end where ! end where lnd_frc_out >= 0.5
    
    if(bsn_fct_hrz) then
       ! Compute source efficiency factor from high-resolution topography then rebin basin factor to model grid (still experimental)
       write (6,"(2a)") prg_nm(1:ftn_strlen(prg_nm)),": Computing basin factor at high resolution"
       call rdb_fct_tpg_GCT01(lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
            hgt_sfc_in,oro_in, & ! I
            mbl_bsn_fct_in) ! O
       
       ! Rebin basin factor from input grid to output grid
       call map_rbn(mbl_bsn_fct_in,mss_flg,mss_val, & ! I
            lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
            lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
            ovr_nbr_max,ovr_lat_idx,ovr_lon_idx,ovr_nbr,ovr_wgt, & ! I
            mbl_bsn_fct_out) ! O
    endif ! endif .false
    
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Exiting hgt_sfc_get_old()"
    return 
  end subroutine hgt_sfc_get_old
  
end module hgt_sfc ! [mdl] Surface elevation, bathymetry
