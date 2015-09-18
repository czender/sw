! $Id$ -*-f90-*-

! Purpose: Routines to process Community Land Model (CLM) fields

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
! use clm_fld ! [mdl] Community Land Model fields

module clm_fld ! [mdl] Community Land Model fields
  implicit none
  public::clm_get ! [sbr] CLM LAI & VAI (inelegant)
  public::lai_get ! [sbr] CLM LAI
  public::pft_get ! [sbr] CLM PFT
  
contains
  
  subroutine clm_get(fl_in, & ! I
       lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr,time_in_nbr, & ! I
       lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr,time_out_nbr, & ! I
       area_out,ovr_nbr_max, & ! I
       lai_ttl_clm_out,vai_ttl_clm_out) ! O
    ! Purpose: Assimilate vegetation data from CLM raw source data and rebin it to output grid
    ! clm_get() is called by bds_prc()
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use map_cst ! [mdl] Constants used in map routines
    use map_grd ! [mdl] Map grids and regridding
    use netcdf ! [mdl] netCDF interface
    use nf90_utl ! [mdl] netCDF utilities
    use sng_mdl,only:ftn_strlen ! [mdl] String manipulation
    use utl_mdl,only:mnt_chk ! [mdl] Utility functions (date_time_get,mnt_chk...)
    implicit none
    ! Parameters
    integer,parameter::sgs_nbr_CLM=4 ! [nbr] Number of sub-gridscale patches in CLM 
    ! Commons
    ! Input
    character,intent(in)::fl_in*80 ! [sng] Input file
    integer,intent(in)::lat_in_nbr ! [nbr] Number of latitudes
    integer,intent(in)::lat_out_nbr ! [nbr] Number of latitudes
    integer,intent(in)::lon_in_nbr ! [nbr] Number of longitudes
    integer,intent(in)::lon_out_nbr ! [nbr] Number of longitudes
    integer,intent(in)::ovr_nbr_max ! [nbr] Maximum number of input cells which overlap any output cell
    integer,intent(in)::time_in_nbr ! [nbr] Number of records (months) in input file
    integer,intent(in)::time_out_nbr ! [nbr] Dimension size
    real,intent(in)::area_out(lon_out_nbr,lat_out_nbr) ! [m2] Area of gridcells
    real,intent(in)::lat_in_grd(lat_in_nbr+1) ! [dgr] Interface latitudes
    real,intent(in)::lat_out_grd(lat_out_nbr+1) ! [dgr] Interface latitudes
    real,intent(in)::lon_in_grd(lon_in_nbr+1) ! [dgr] Interface longitudes
    real,intent(in)::lon_out_grd(lon_out_nbr+1) ! [dgr] Interface longitudes
    ! Output
    real,intent(out)::lai_ttl_clm_out(lon_out_nbr,lat_out_nbr,time_out_nbr) ! [m2 m-2] Leaf area index, one-sided
    real,intent(out)::vai_ttl_clm_out(lon_out_nbr,lat_out_nbr,time_out_nbr) ! [m2 m-2] Vegetation area index, one-sided
    
    ! Locals with simple initialization
    integer::rcd=nf90_noerr ! [enm] Return success code
    logical::mss_flg=.false. ! [flg] Data may have missing values
    real::mss_val_in=1.0e36 ! [frc] Missing value
    
    ! Local
    integer lat_in_idx ! [idx] Counting index for lat
    integer lon_in_idx ! [idx] Counting index for lon
    integer ovr_idx ! [idx] Counting index
    integer ovr_lat_idx(lon_out_nbr,lat_out_nbr,ovr_nbr_max) ! [idx] Map into input grid of latitude indices of overlap cells
    integer ovr_lon_idx(lon_out_nbr,lat_out_nbr,ovr_nbr_max) ! [idx] Map into input grid of longitude indices of overlap cells
    integer ovr_nbr(lon_out_nbr,lat_out_nbr) ! [nbr] Number of input gridcells which overlap each output gridcell
    integer ovr_nbr_crr ! [nbr] Current number of overlapping gridcells
    integer sgs_idx ! [idx] Counting index for sgs
    integer time_in_idx ! [idx] Counting index for time
    integer time_nbr ! [nbr] Number of records (months) in input file
    integer time_out_idx ! [idx] Counting index for time
    logical mnt ! [flg] Monotonicity flag
    real ovr_wgt(lon_out_nbr,lat_out_nbr,ovr_nbr_max) ! [frc] Weight of overlapping input gridcells onto each output gridcell
    real sai_clm_out(lon_out_nbr,lat_out_nbr,time_out_nbr) ! [m2 m-2] Stem area index, one-sided
    
    ! Variables needed to read external netCDF dataset
    integer lai_clm_id ! [id] Variable ID
    integer lai_dmn_nbr ! [nbr] Number of dimensions in disk lai_clm field
    integer lat_dmn_id ! [id] Dimension ID for lat
    integer lat_nbr ! [nbr] Number of latitudes in input file
    integer lon_dmn_id ! [id] Dimension ID for lon
    integer lon_nbr ! [nbr] Number of longitudes in input file
    integer mss_val_id ! [id] Attribute ID
    integer nc_id ! [id] File handle
    integer pft_frc_clm_id ! [id] Variable ID
    integer sai_clm_id ! [id] Variable ID
    integer sgs_dmn_id ! [id] Dimension ID for sgs
    integer sgs_nbr ! [nbr] Number of sub-gridscale patches per gridcell
    integer time_dmn_id ! [id] Dimension ID for time
    
    ! Allocatable
    real,dimension(:,:,:,:),allocatable::lai_clm_in ! [m2 m-2] Leaf area index, one-sided
    real,dimension(:,:,:,:),allocatable::sai_clm_in ! [m2 m-2] Stem area index, one-sided
    real,dimension(:,:,:),allocatable::lai_ttl_clm_in ! [m2 m-2] Leaf area index, one-sided
    real,dimension(:,:,:),allocatable::vai_ttl_clm_in ! [m2 m-2] Vegetation area index, one-sided
    real,dimension(:,:,:),allocatable::pft_frc_clm_in ! [frc] Fraction covered by plant functional type
    
    ! Main Code
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Entering clm_get()"
    
    ! Enough memory?
    if (time_out_nbr /= 12) stop "clm_get(): time_out_nbr /= 12"
    
    ! Check for monotonicity
    mnt=mnt_chk(lon_out_grd,lon_out_nbr+1)
    if (.not.mnt) stop "lon_out_grd not monotonic in clm_get()"
    mnt=mnt_chk(lat_out_grd,lat_out_nbr+1)
    if (.not.mnt) stop "lat_out_grd not monotonic in clm_get()"
    ! Check for increasing monotonicity
    if (lon_out_grd(2) < lon_out_grd(1)) stop "lon_out_grd not increasing in clm_get()"
    if (lat_out_grd(2) < lat_out_grd(1)) stop "lat_out_grd not increasing in clm_get()"
    
    ! Get overlap locations and weights for mapping input grid to output grid
    call map_ovr_wgt_drv( &
         lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
         lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
         ovr_nbr_max,area_out, & ! I
         ovr_lat_idx,ovr_lon_idx,ovr_nbr,ovr_wgt) ! O
    
    ! Initialize arrays
    lai_ttl_clm_out(:,:,:)=0.0 ! [m2 m-2] Leaf area index, one-sided
    vai_ttl_clm_out(:,:,:)=0.0 ! [m2 m-2] Vegetation area index, one-sided
    
    ! Read netCDF data
    rcd=rcd+nf90_wrp_open(fl_in,nf90_nowrite,nc_id)
    ! Get dimension IDs
    rcd=nf90_wrp_inq_dimid(nc_id,"lsmlon",lon_dmn_id)
    rcd=nf90_wrp_inq_dimid(nc_id,"lsmlat",lat_dmn_id)
    rcd=nf90_wrp_inq_dimid(nc_id,"lsmpft",sgs_dmn_id)
    rcd=nf90_wrp_inq_dimid(nc_id,"time",time_dmn_id)
    if(rcd /= nf90_noerr) stop "Error retrieving dimension IDs"
    ! Get dimension sizes
    rcd=rcd+nf90_inquire_dimension(nc_id,lon_dmn_id,len=lon_nbr)
    rcd=rcd+nf90_inquire_dimension(nc_id,lat_dmn_id,len=lat_nbr)
    rcd=rcd+nf90_inquire_dimension(nc_id,sgs_dmn_id,len=sgs_nbr)
    rcd=rcd+nf90_inquire_dimension(nc_id,time_dmn_id,len=time_nbr)
    if(rcd /= nf90_noerr) stop "Error retrieving dimension sizes"
    ! Enough memory? 
    if (lat_nbr /= lat_in_nbr) stop "lat_nbr /= lat_in_nbr in clm_get()"
    if (lon_nbr /= lon_in_nbr) stop "lon_nbr /= lon_in_nbr in clm_get()"
    if (sgs_nbr /= sgs_nbr_CLM) stop "sgs_nbr /= sgs_nbr_CLM in clm_get()"
    if (time_nbr /= time_in_nbr) stop "time_nbr /= time_in_nbr in clm_get()"
    ! Allocate based on input file dimensions
    allocate(lai_clm_in(lon_in_nbr,lat_in_nbr,sgs_nbr,time_in_nbr),stat=rcd) ! [m2 m-2] Leaf area index, one-sided
    if(rcd /= 0) stop "allocate() failed for lai_clm_in"
    allocate(sai_clm_in(lon_in_nbr,lat_in_nbr,sgs_nbr,time_in_nbr),stat=rcd) ! [m2 m-2] Stem area index, one-sided
    if(rcd /= 0) stop "allocate() failed for sai_clm_in"
    allocate(lai_ttl_clm_in(lon_in_nbr,lat_in_nbr,time_in_nbr),stat=rcd) ! [m2 m-2] Leaf area index, one-sided
    if(rcd /= 0) stop "allocate() failed for lai_ttl_clm_in"
    allocate(vai_ttl_clm_in(lon_in_nbr,lat_in_nbr,time_in_nbr),stat=rcd) ! [m2 m-2] Vegetation area index, one-sided
    if(rcd /= 0) stop "allocate() failed for vai_ttl_clm_in"
    allocate(pft_frc_clm_in(lon_in_nbr,lat_in_nbr,sgs_nbr),stat=rcd) ! [frc] Fraction covered by plant functional type
    if(rcd /= 0) stop "allocate() failed for pft_frc_clm_in"
    
    if (time_nbr > 1) then
       write (6,"(2a,f12.6)") prg_nm(1:ftn_strlen(prg_nm)), & 
            ": ERROR clm_get() reports time_nbr = ",time_nbr
       stop "time_nbr > 1"
    endif                  ! endif err
    ! Get variable IDs
    rcd=nf90_wrp_inq_varid(nc_id,"MONTHLY_LAI",lai_clm_id)
    rcd=nf90_wrp_inq_varid(nc_id,"MONTHLY_SAI",sai_clm_id)
    rcd=nf90_wrp_inq_varid(nc_id,"PCT_PFT",pft_frc_clm_id)
    ! Get number of dimensions
    rcd=rcd+nf90_inquire_variable(nc_id,lai_clm_id,ndims=lai_dmn_nbr)
    if (lai_dmn_nbr /= 4) then
       write (6,"(2a,i1,a)") prg_nm(1:ftn_strlen(prg_nm)), & 
            ": ERROR MONTHLY_LAI has ",lai_dmn_nbr," dimensions"
       write (6,"(2a,i1,a)") prg_nm(1:ftn_strlen(prg_nm)), & 
            ": HINT Make sure MONTHLY_LAI has phenotype and time dimensions"
       stop
    endif                  ! endif err
    ! Get data
    rcd=nf90_wrp(nf90_get_var(nc_id,lai_clm_id,lai_clm_in),"get_var lai_clm_in")
    rcd=nf90_wrp(nf90_get_var(nc_id,pft_frc_clm_id,pft_frc_clm_in),"get_var pft_frc_clm_in")
    rcd=nf90_wrp(nf90_get_var(nc_id,sai_clm_id,sai_clm_in),"get_var sai_clm_in")
    rcd=rcd+nf90_wrp_close(nc_id,fl_in,"Ingested")
    
    ! Compute derived fields on input grid 
    ! Note: Do not attempt to rebin pft_clm and pft_frc_clm using continuous methods
    ! pft_clm is enumerated integer not continuous field
    ! pft_clm can only be properly remapped using, e.g., a ranking method
    ! Ideally, we should re-generate pft_clm for output grid from raw data
    ! This is done in routine pft_get()
    do lat_in_idx=1,lat_in_nbr
       do lon_in_idx=1,lon_in_nbr
          do sgs_idx=1,sgs_nbr_CLM
             do time_in_idx=1,time_in_nbr
                lai_ttl_clm_in(lon_in_idx,lat_in_idx,time_in_idx)= & ! [m2 m-2]
                     lai_ttl_clm_in(lon_in_idx,lat_in_idx,time_in_idx)+ & ! [m2 m-2]
                     pft_frc_clm_in(lon_in_idx,lat_in_idx,sgs_idx)* & ! [frc]
                     lai_clm_in(lon_in_idx,lat_in_idx,sgs_idx,time_in_idx) ! [m2 m-2]
                vai_ttl_clm_in(lon_in_idx,lat_in_idx,time_in_idx)= & ! [m2 m-2]
                     vai_ttl_clm_in(lon_in_idx,lat_in_idx,time_in_idx)+ & ! [m2 m-2]
                     pft_frc_clm_in(lon_in_idx,lat_in_idx,sgs_idx)* & ! [frc]
                     (lai_clm_in(lon_in_idx,lat_in_idx,sgs_idx,time_in_idx)+ & ! [m2 m-2]
                     sai_clm_in(lon_in_idx,lat_in_idx,sgs_idx,time_in_idx)) ! [m2 m-2]
             end do ! end loop over time
          end do ! end loop over number of possible plant types
       end do ! end loop over lon
    end do ! end loop over lat
    
    ! De-allocate
    if (allocated(lai_clm_in)) deallocate(lai_clm_in,stat=rcd)
    if(rcd /= 0) stop "deallocate() failed for lai_clm_in"
    if (allocated(sai_clm_in)) deallocate(sai_clm_in,stat=rcd)
    if(rcd /= 0) stop "deallocate() failed for sai_clm_in"
    
    ! Loop over months because map_rbn does not yet work for 3D fields
    do time_out_idx=1,time_out_nbr
       
       ! Rebin total leaf area index from input grid to output grid
       call map_rbn(lai_ttl_clm_in(1,1,time_out_idx),mss_flg,mss_val_in, & ! I
            lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
            lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
            ovr_nbr_max,ovr_lat_idx,ovr_lon_idx,ovr_nbr,ovr_wgt, & ! I
            lai_ttl_clm_out(1,1,time_out_idx)) ! O
       
       ! Rebin total vegetation area index from input grid to output grid
       call map_rbn(vai_ttl_clm_in(1,1,time_out_idx),mss_flg,mss_val_in, & ! I
            lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
            lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
            ovr_nbr_max,ovr_lat_idx,ovr_lon_idx,ovr_nbr,ovr_wgt, & ! I
            vai_ttl_clm_out(1,1,time_out_idx)) ! O
       
    end do ! end loop over mth
    
    ! De-allocate
    if (allocated(lai_ttl_clm_in)) deallocate(lai_ttl_clm_in,stat=rcd)
    if(rcd /= 0) stop "deallocate() failed for lai_ttl_clm_in"
    if (allocated(vai_ttl_clm_in)) deallocate(vai_ttl_clm_in,stat=rcd)
    if(rcd /= 0) stop "deallocate() failed for vai_ttl_clm_in"
    
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Exiting clm_get()"
    return 
  end subroutine clm_get                       ! end clm_get()
  
  subroutine pft_get(fl_in, & ! I
       lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr,sgs_out_nbr, & ! I
       area_out,lnd_msk_out, & ! I
       lnd_frc_out,pft_clm_out,pft_frc_clm) ! O
    ! Purpose: Retrieve plant functional type and rebin to requested grid
    ! pft_get() is called by bds_prc()
    ! Plant functional type "raw" data are rebinned from PFT dataset assembled by ....
    ! Dataset resolution is 0.5 x 0.5 degree and data are written from South pole to North pole, West to East, beginning at date line
    
    ! Mixture of _in, _out, and no suffixes can be confusing:
    ! Most fields require only a var_in and a var_out array and intermediate
    ! processing may be performed on both the var_in and var_out arrays prior
    ! to returning var_out to the host model. pft_clm* is one such field. 
    ! Internally to the pft_get() subroutine:
    ! pft_clm_in is raw data on input grid
    ! pft_clm_out is processed data returned on output grid
    ! pft_clm is non-existent because it is not needed
    ! pft_frc_clm*, however, has 3 versions because its processing is more complex:
    ! pft_frc_clm_in is raw data on input grid (lon,lat,pft,time)
    ! pft_frc_clm_out is (nearly) raw data rebinned to output spatial grid (lon,lat,pft,time)
    ! pft_frc_clm is processed data returned on output spatial, sgs grid (lon,lat,sgs,time)
    ! Externally to the pft_get() subroutine:
    ! pft_clm_out is non-existent because it is not needed
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use map_grd ! [mdl] Map grids and regridding
    use netcdf ! [mdl] netCDF interface
    use nf90_utl ! [mdl] netCDF utilities
    use sng_mdl,only:ftn_strlen ! [mdl] String manipulation
    use utl_mdl,only:mnt_ncr_chk ! [mdl] Utility functions (date_time_get,mnt_chk...)
    implicit none
    ! Parameters
    character(len=*),parameter::sbr_nm="pft_get()" ! [sng] Subroutine name
    integer,parameter::lat_idx_dbg=50 ! [idx] Longitude for verbose output
    integer,parameter::lon_idx_dbg=90 ! [idx] Longitude for verbose output
    integer,parameter::idx_mss_val=99999 ! [idx] Index denoting missing/invalid data
    integer,parameter::pft_nbr_CLM=17 ! [nbr] Number of plant functional types in CLM 
    integer,parameter::pft_non_vgt_CLM=0 ! [enm] Non-vegetated plant functional type in CLM
    integer,parameter::sgs_nbr_CLM=4 ! [nbr] Number of sub-gridscale patches in CLM 
    integer,parameter::unit_dgn=74 ! Unit for writing general Map diagnoses
    real,parameter::eps_rlt=1.0e-5 ! [frc] Relative error allowed in frc_ttl
    ! Derived parameters
    integer,parameter::pft_max_CLM=pft_nbr_CLM-1 ! [enm] Maximum value of PFT in CLM
    ! Commons
    ! Input
    character,intent(in)::fl_in*80 ! [sng] Input file
    integer,intent(in)::lat_out_nbr ! [nbr] Number of latitudes
    integer,intent(in)::lon_out_nbr ! [nbr] Number of longitudes
    integer,intent(in)::sgs_out_nbr ! [nbr] Number of sub-gridscale patches per gridcell
    integer,intent(in)::lnd_msk_out(lon_out_nbr,lat_out_nbr) ! [msk] Land mask (integer 0 or 1)
    real,intent(in)::area_out(lon_out_nbr,lat_out_nbr) ! [m2] Area of gridcells
    real,intent(in)::lat_out_grd(lat_out_nbr+1) ! [dgr] Interface latitudes
    real,intent(in)::lon_out_grd(lon_out_nbr+1) ! [dgr] Interface longitudes
    ! Output
    integer,intent(out)::pft_clm_out(lon_out_nbr,lat_out_nbr,sgs_out_nbr) ! [enm] Plant functional type
    real,intent(out)::pft_frc_clm(lon_out_nbr,lat_out_nbr,sgs_out_nbr) ! [frc] Fraction covered by plant functional type
    real,intent(out)::lnd_frc_out(lon_out_nbr,lat_out_nbr) ! [frc] Fraction of land (not ocean)
    
    ! Locals with simple initialization
    character(80)::fl_dgn="/tmp/zender/map/pft.txt"//char(0) ! [sng] Diagnostics file
    character(35) pft_sng(0:pft_max_CLM) ! [sng] Description of plant functional types
    data pft_sng( 0) /"not vegetated"                      /
    data pft_sng( 1) /"needleleaf evergreen temperate tree"/
    data pft_sng( 2) /"needleleaf evergreen boreal tree"   /
    data pft_sng( 3) /"needleleaf deciduous boreal tree"   /
    data pft_sng( 4) /"broadleaf evergreen tropical tree"  /
    data pft_sng( 5) /"broadleaf evergreen temperate tree" /
    data pft_sng( 6) /"broadleaf deciduous tropical tree"  /
    data pft_sng( 7) /"broadleaf deciduous temperate tree" /
    data pft_sng( 8) /"broadleaf deciduous boreal tree"    /
    data pft_sng( 9) /"broadleaf evergreen shrub"          /
    data pft_sng(10) /"broadleaf deciduous temperate shrub"/
    data pft_sng(11) /"broadleaf deciduous boreal shrub"   /
    data pft_sng(12) /"c3 arctic grass"                    /
    data pft_sng(13) /"c3 non-arctic grass"                /
    data pft_sng(14) /"c4 grass"                           / 
    data pft_sng(15) /"corn"                               /
    data pft_sng(16) /"wheat"                              /
    integer::rcd=nf90_noerr ! [enm] Return success code
    real::mss_val_in=1.0e36 ! [frc] Missing value
    
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
    integer pft_idx ! [idx] Counting index for pft
    integer pft_in_idx ! [idx] Counting index for pft
    integer pft_max ! [enm] Maximum plant functional type
    integer pft_min ! [enm] Minimum plant functional type
    integer pft_nbr ! [nbr] Number of plant functional types
    integer pft_rnk(sgs_out_nbr) ! [idx] Ranked list of plant functional type fraction
    integer rnk_vld_nbr ! [nbr] Number of valid, ranked entries returned
    integer sgs_idx ! [idx] Counting index for sgs
    logical lat_mnt_ncr ! [flg] Latitude monotonically decreases
    logical mnt_ncr ! [flg] Monotonic and increasing flag
    logical mss_flg ! [flg] Variable has missing_value attribute
    real edg_grd(4) ! [dgr] Grid edges (north, east, south, west)
    real frc_dff ! [frc] Excess PFT fraction to redistribute
    real frc_ttl ! [frc] Total fraction of gridcell accounted for by all plant functional types
    real frc_ttl_pft ! [frc] Sum of plant functional type fractions
    real frc_ttl_sgs ! [frc] Sum of sub-gridscale fractions
    real ovr_wgt_crr ! [frc] Overlap weight of current gridcell
    real pft_frc_clm_crr(sgs_out_nbr) ! [frc] Initial PFT fraction of current gridcell
    
    ! Variables needed to read external netCDF dataset
    integer edg_est_id ! [id] Variable ID
    integer edg_nrt_id ! [id] Variable ID
    integer edg_sth_id ! [id] Variable ID
    integer edg_wst_id ! [id] Variable ID
    integer lat_dmn_id ! [id] Dimension ID for lat
    integer lat_in_2d_id ! [id] Variable ID
    integer lat_in_id ! [id] Variable ID
    integer lnd_frc_clm_id ! [id] Variable ID
    integer lon_dmn_id ! [id] Dimension ID for lon
    integer lon_in_2d_id ! [id] Variable ID
    integer mss_val_id ! [id] Attribute ID
    integer nc_id ! [id] File handle
    integer pft_dmn_id ! [id] Dimension ID for pft
    integer pft_frc_clm_dmn_nbr ! [nbr] Number of dimensions in disk pft_frc_clm field
    integer pft_frc_clm_id ! [id] Variable ID
    
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
    
    real,dimension(:,:),allocatable::lnd_frc_in ! [frc] Fraction of land (not ocean)
    real,dimension(:,:),allocatable::lnd_msk_in ! [msk] Land mask (0.0 or 1.0)
    real,dimension(:,:,:),allocatable::pft_frc_clm_in ! [frc] Fraction covered by plant functional type
    real,dimension(:,:,:),allocatable::pft_frc_clm_out ! [frc] Fraction covered by plant functional type
    real,dimension(:),allocatable::pft_frc_clm_out_crr ! [frc] Fraction covered by plant functional type in current cell
    
    ! Main Code
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Entering "//sbr_nm
    
    ! Check for monotonically increasing grids
    mnt_ncr=mnt_ncr_chk(lon_out_grd,lon_out_nbr+1)
    if (.not.mnt_ncr) stop "lon_out_grd not monotonically increasing in pft_get()"
    mnt_ncr=mnt_ncr_chk(lat_out_grd,lat_out_nbr+1)
    if (.not.mnt_ncr) stop "lat_out_grd not monotonically increasing in pft_get()"
    ! Sanity check
    if (sgs_out_nbr > pft_nbr_CLM) stop "sgs_out_nbr > pft_nbr_CLM"
    
    ! Initialize arrays
    pft_clm_out(:,:,:)=0 ! [enm] Plant functional type
    pft_frc_clm(:,:,:)=0.0 ! [frc] Fraction covered by plant functional type
    lnd_frc_out(:,:)=0.0 ! [frc] Fraction of land (not ocean)
    
    ! Read netCDF data
    rcd=rcd+nf90_wrp_open(fl_in,nf90_nowrite,nc_id)
    ! Get dimension IDs
    rcd=nf90_wrp_inq_dimid(nc_id,"lon",lon_dmn_id)
    rcd=nf90_wrp_inq_dimid(nc_id,"lat",lat_dmn_id)
    rcd=nf90_wrp_inq_dimid(nc_id,"pft",pft_dmn_id)
    if(rcd /= nf90_noerr) stop "Error retrieving dimension IDs"
    ! Get dimension sizes
    rcd=rcd+nf90_inquire_dimension(nc_id,lon_dmn_id,len=lon_in_nbr)
    rcd=rcd+nf90_inquire_dimension(nc_id,lat_dmn_id,len=lat_in_nbr)
    rcd=rcd+nf90_inquire_dimension(nc_id,pft_dmn_id,len=pft_nbr)
    if(rcd /= nf90_noerr) stop "Error retrieving dimension sizes"
    ! Enough memory? 
    if (lat_in_nbr /= 360) stop "lat_in_nbr /= 360 in pft_get()"
    if (lon_in_nbr /= 720) stop "lon_in_nbr /= 720 in pft_get()"
    if (pft_nbr /= pft_nbr_CLM) stop "pft_nbr /= pft_nbr_CLM in pft_get()"
    ! Allocate based on input file dimensions
    allocate(area_in(lon_in_nbr,lat_in_nbr),stat=rcd) ! [m2] Area of gridcells
    if(rcd /= 0) stop "allocate() failed for area_in"
    allocate(lat_in_grd(lat_in_nbr+1),stat=rcd) ! [dgr] Interface latitudes
    if(rcd /= 0) stop "allocate() failed for lat_in_grd"
    allocate(lat_in_2d(lon_in_nbr,lat_in_nbr),stat=rcd) ! [dgr] Latitude at gridcell center
    if(rcd /= 0) stop "allocate() failed for lat_in_2d"
    allocate(lon_in_grd(lon_in_nbr+1),stat=rcd) ! [dgr] Interface longitudes
    if(rcd /= 0) stop "allocate() failed for lon_in_grd"
    allocate(lon_in_grd_2d(lon_in_nbr+1,lat_in_nbr),stat=rcd) ! [dgr] Interface longitudes
    if(rcd /= 0) stop "allocate() failed for lon_in_grd_2d"
    allocate(lon_in_2d(lon_in_nbr,lat_in_nbr),stat=rcd) ! [dgr] Longitude at gridcell center
    if(rcd /= 0) stop "allocate() failed for lon_in_2d"
    allocate(lon_in_nbr_1d(lat_in_nbr),stat=rcd) ! [nbr] Number of longitudes per latitude
    if(rcd /= 0) stop "allocate() failed for lon_in_nbr_1d"
    allocate(lnd_frc_in(lon_in_nbr,lat_in_nbr),stat=rcd) ! [frc] Fraction of land (not ocean)
    if(rcd /= 0) stop "allocate() failed for lnd_frc_in"
    allocate(lnd_msk_in(lon_in_nbr,lat_in_nbr),stat=rcd) ! [msk] Land mask (0.0 or 1.0)
    if(rcd /= 0) stop "allocate() failed for lnd_msk_in"
    
    allocate(pft_frc_clm_in(lon_in_nbr,lat_in_nbr,0:pft_nbr-1),stat=rcd) ! [frc] Fraction covered by plant functional type
    if(rcd /= 0) stop "allocate() failed for pft_frc_clm_in"
    allocate(pft_frc_clm_out(lon_out_nbr,lat_out_nbr,0:pft_nbr-1),stat=rcd) ! [frc] Fraction covered by plant functional type
    if(rcd /= 0) stop "allocate() failed for pft_frc_clm_out"
    allocate(pft_frc_clm_out_crr(0:pft_nbr-1),stat=rcd) ! [frc] Fraction covered by plant functional type in current cell
    if(rcd /= 0) stop "allocate() failed for pft_frc_clm_out_crr"
    ! Get variable IDs
    rcd=nf90_wrp_inq_varid(nc_id,"LATIXY",lat_in_2d_id)
    rcd=nf90_wrp_inq_varid(nc_id,"LONGXY",lon_in_2d_id)
    rcd=nf90_wrp_inq_varid(nc_id,"EDGEN",edg_nrt_id)
    rcd=nf90_wrp_inq_varid(nc_id,"EDGEE",edg_est_id)
    rcd=nf90_wrp_inq_varid(nc_id,"EDGES",edg_sth_id)
    rcd=nf90_wrp_inq_varid(nc_id,"EDGEW",edg_wst_id)
    rcd=nf90_wrp_inq_varid(nc_id,"PCT_PFT",pft_frc_clm_id)
    rcd=nf90_wrp_inq_varid(nc_id,"LANDMASK",lnd_frc_clm_id)
    ! Get number of dimensions
    rcd=rcd+nf90_inquire_variable(nc_id,pft_frc_clm_id,ndims=pft_frc_clm_dmn_nbr)
    if (pft_frc_clm_dmn_nbr /= 3) then
       write (6,"(2a,i1,a)") prg_nm(1:ftn_strlen(prg_nm)), & 
            ": ERROR PCT_PFT has ",pft_frc_clm_dmn_nbr," dimensions"
       write (6,"(2a,i1,a)") prg_nm(1:ftn_strlen(prg_nm)), & 
            ": HINT Make sure PCT_PFT has phenotype dimension but no time dimension"
       stop
    endif                  ! endif err
    ! Get data
    rcd=nf90_wrp(nf90_get_var(nc_id,lat_in_2d_id,lat_in_2d),"get_var lat_in_2d")
    rcd=nf90_wrp(nf90_get_var(nc_id,lon_in_2d_id,lon_in_2d),"get_var lon_in_2d")
    rcd=nf90_wrp(nf90_get_var(nc_id,edg_nrt_id,edg_grd(1)),"get_var edg_nrt")
    rcd=nf90_wrp(nf90_get_var(nc_id,edg_est_id,edg_grd(2)),"get_var edg_est")
    rcd=nf90_wrp(nf90_get_var(nc_id,edg_sth_id,edg_grd(3)),"get_var edg_sth")
    rcd=nf90_wrp(nf90_get_var(nc_id,edg_wst_id,edg_grd(4)),"get_var edg_wst")
    rcd=nf90_wrp(nf90_get_var(nc_id,pft_frc_clm_id,pft_frc_clm_in),"get_var pft_frc_clm_in")
    rcd=nf90_wrp(nf90_get_var(nc_id,lnd_frc_clm_id,lnd_frc_in),"get_var lnd_frc_in")
    if (rcd /= nf90_noerr) stop "Error in nf90_get_var in pft_get()"
    ! Get missing value
    rcd=rcd+nf90_inquire_attribute(nc_id,pft_frc_clm_id,"missing_value",attnum=mss_val_id)
    if (rcd==nf90_noerr) then
       mss_flg=.true. ! [flg] Data may have missing values
       rcd=rcd+nf90_get_att(nc_id,pft_frc_clm_id,"missing_value",mss_val_in)
       if (mss_val_in /= 1.0e36) write (6,"(2a,f12.6)") prg_nm(1:ftn_strlen(prg_nm)), &
            ": WARNING "//sbr_nm//" reports mss_val_in = ",mss_val_in
    else
       mss_flg=.false. ! [flg] Data may have missing values
       rcd=nf90_noerr ! [enm] Return success code
    endif                  ! endif
    rcd=rcd+nf90_wrp_close(nc_id,fl_in,"Ingested plant functional type data from") ! [fnc] Close file
    
    ! Convert input data to SI
    do pft_in_idx=0,pft_max_CLM ! NB: Loop starts at zero
       do lat_in_idx=1,lat_in_nbr
          do lon_in_idx=1,lon_in_nbr
             ! pft_frc_clm_in is stored in percent
             if (mss_flg.and.pft_frc_clm_in(lon_in_idx,lat_in_idx,pft_in_idx)==mss_val_in) then 
                pft_frc_clm_in(lon_in_idx,lat_in_idx,pft_in_idx)=0.0 ! [frc] Fraction covered by plant functional type
             else
                pft_frc_clm_in(lon_in_idx,lat_in_idx,pft_in_idx)= & ! [frc] Fraction covered by plant functional type
                     0.01*pft_frc_clm_in(lon_in_idx,lat_in_idx,pft_in_idx)
             end if ! end if mss_flg
          end do ! end loop over lon
       end do ! end loop over lat
    end do ! end loop over pft
    
    ! Determine map grid and overlap characteristics
    ! fxm: Method breaks for reduced grids, streamline, subroutinize, and generalize
    lon_in_nbr_1d(:)=lon_in_nbr ! [nbr] Number of longitudes per latitude
    
    ! Convert gridpoint centers, domain boundaries to grid interfaces
    call map_edge_mk(lat_in_nbr,lon_in_nbr,lon_in_nbr_1d, & ! I
         lat_in_2d,lon_in_2d,edg_grd, & ! I 
         lat_in_grd,lon_in_grd,lon_in_grd_2d) ! O
    
    mnt_ncr=mnt_ncr_chk(lat_in_grd,lat_in_nbr+1)
    if (.not.mnt_ncr) stop "ERROR: lat_in_grd not monotonically increasing in pft_get()"
    mnt_ncr=mnt_ncr_chk(lon_in_grd,lon_in_nbr+1)
    if (.not.mnt_ncr) stop "ERROR: lon_in_grd not monotonically increasing in pft_get()"
    
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
    
    ! Archive input in human-readable format if desired
    if (.false.) then
       open (unit=unit_dgn,file=fl_dgn,status="unknown",iostat=rcd)
       if (rcd /= 0) write (6,"(4a,i4)") prg_nm(1:ftn_strlen(prg_nm)), & 
            ": ERROR "//sbr_nm//" unable to open ",fl_dgn(1:ftn_strlen(fl_dgn))
       do lat_in_idx=1,lat_in_nbr ! NB: Outer loop over lat
          do lon_in_idx=1,lon_in_nbr
             do pft_in_idx=0,pft_max_CLM ! NB: Loop starts at zero
                write (unit_dgn,"(f10.4)") &
                     pft_frc_clm_in(lon_in_idx,lat_in_idx,pft_in_idx)
             end do ! end loop over pft
          end do ! end loop over lon
       end do ! end loop over lat
       close (unit_dgn)
       write (6,"(a,1x,a)") "Wrote plant functional type data to",fl_dgn(1:ftn_strlen(fl_dgn))
    endif ! endif dbg
    
    ! Sanity check
    err_nbr=0
    do lat_in_idx=1,lat_in_nbr
       do lon_in_idx=1,lon_in_nbr
          frc_ttl=sum(pft_frc_clm_in(lon_in_idx,lat_in_idx,:))
          if ( &
               (lnd_frc_in(lon_in_idx,lat_in_idx) < 0.0.or.lnd_frc_in(lon_in_idx,lat_in_idx) > 1.0).or. &
               (lnd_frc_in(lon_in_idx,lat_in_idx)==1.0.and.abs(frc_ttl-1.0) > eps_rlt)) then
             err_nbr=err_nbr+1
             lon_in_idx_err=lon_in_idx
             lat_in_idx_err=lat_in_idx
          end if ! endif err
       end do ! end loop over lon
    end do ! end loop over lat
    if (err_nbr > 0) then
       frc_ttl=sum(pft_frc_clm_in(lon_in_idx_err,lat_in_idx_err,:))
       write (6,"(a,i6,a)") "ERROR "//sbr_nm//" reports ",err_nbr," errors during input"
       write (6,"(a,f15.12)") "Input land fraction: ",lnd_frc_in(lon_in_idx_err,lat_in_idx_err)
       write (6,"(a,f15.12)") "Plant functional type fractions sum to ",frc_ttl
       write (6,"(a3,1x,a35,1x,a)") "PFT","Description                        ","Fraction"
       do pft_idx=0,pft_max_CLM ! NB: Loop starts at zero
          write (6,"(i3,1x,a35,1x,f9.6)") pft_idx,pft_sng(pft_idx), &
               pft_frc_clm_in(lon_in_idx_err,lat_in_idx_err,pft_idx)
       end do ! end loop over pft
       write (6,"(a)") "Input cell edge locations:"
       write (6,"(a)") "(idx)   Lat Sth (idx)   Lat Nrt (idx)   Lon Wst (idx)   Lon Est"
       write (6,"(4(a1,i3,a1,1x,f9.4,1x))") &
            "(",lat_in_idx_err,")",lat_in_grd(lat_in_idx_err),"(",lat_in_idx_err+1,")",lat_in_grd(lat_in_idx_err+1), &
            "(",lon_in_idx_err,")",lon_in_grd(lon_in_idx_err),"(",lon_in_idx_err+1,")",lon_in_grd(lon_in_idx_err+1)
       stop
    endif ! endif err
    
    ! Land mask is unity for LSM/CLM grid cells and 0 for ocean points
    ! Use land mask in array as weights compute land fraction on output grid
    ! For this purpose, temporarily set input land mask to unity everywhere
    lnd_msk_in(:,:)=1.0 ! [msk] Land mask (0.0 or 1.0)
    
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
    
    ! lnd_msk_out is land mask of host model, two possible integer values (0, 1):
    ! 1 if there is ANY land in output gridcell (then apply lnd_frc_out to fluxes)
    ! 0 if there is NO  land in output gridcell (do not compute land fluxes)
    
    ! Land mask is unity for LSM/CLM grid cells and 0 for ocean points
    ! lnd_frc_in is valid mask because it is 1.0 or 0.0 on input grid
    ! This is only true when input file is "raw data" high resolution grid
    lnd_msk_in(:,:)=lnd_frc_in(:,:) ! [msk] Land mask (0.0 or 1.0)
    
    ! fxm: I"m not sure about purpose of second areaini() in mkpft.F90
    ! so I"ve skipped it temporarily
#if 0
    ! Get overlap locations and weights for mapping input grid to output grid
    ! Use lnd_frc_in as input mask
    call map_ovr_wgt_drv( &
         lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
         lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
         ovr_nbr_max,area_out, & ! I
         ovr_lat_idx,ovr_lon_idx,ovr_nbr,ovr_wgt) ! O
#endif /* endif 0 */
    
    do pft_idx=0,pft_max_CLM ! NB: Loop starts at zero
       ! Rebin plant functional type fractional cover from input grid to output grid
       call map_rbn(pft_frc_clm_in(1,1,pft_idx),mss_flg,mss_val_in, & ! I
            lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
            lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
            ovr_nbr_max,ovr_lat_idx,ovr_lon_idx,ovr_nbr,ovr_wgt, & ! I
            pft_frc_clm_out(1,1,pft_idx)) ! O
    end do ! end loop over pft
    
    ! Post-process data
    do lat_out_idx=1,lat_out_nbr
       do lon_out_idx=1,lon_out_nbr
          frc_ttl=0.0
          do pft_idx=0,pft_max_CLM ! NB: Loop starts at zero
             
             ! If raw data has any land in this output gridcell...
             if (lnd_frc_out(lon_out_idx,lat_out_idx) > 0.0) then
                ! ...and if host model has no land in this output gridcell...
                if (lnd_msk_out(lon_out_idx,lat_out_idx) == 0) &
                     ! ...then set fractional PFT to 0.0 for all PFTs
                     pft_frc_clm_out(lon_out_idx,lat_out_idx,pft_idx)=0.0
             else ! ...else if raw data has no land in this output gridcell...
                ! ...then set fractional PFT to 0.0 for all PFTs...
                pft_frc_clm_out(lon_out_idx,lat_out_idx,pft_idx)=0.0
                ! ...and, if host model has land any land in this output gridcell...
                if (lnd_msk_out(lon_out_idx,lat_out_idx) == 1) &
                     ! ...then set fraction of non-vegetation PFT to 1.0...
                     pft_frc_clm_out(lon_out_idx,lat_out_idx,pft_non_vgt_CLM)=1.0
             end if ! endif
             frc_ttl=frc_ttl+pft_frc_clm_out(lon_out_idx,lat_out_idx,pft_idx)
          end do ! end loop over pft
          ! fxm: Add sanity check on frc_ttl here?
       end do ! end loop over lon
    end do ! end loop over lat
    
    ! Determine dominant plant functional types
    do lat_out_idx=1,lat_out_nbr
       do lon_out_idx=1,lon_out_nbr
          pft_frc_clm_out_crr(:)=pft_frc_clm_out(lon_out_idx,lat_out_idx,:) ! [frc] Fraction covered by plant functional type in current cell
          frc_ttl_pft=sum(pft_frc_clm_out_crr(:)) ! [frc] Sum of plant functional type fractions
          
          if (lnd_msk_out(lon_out_idx,lat_out_idx) == 1) then
             call rnk_vec_gnr(pft_frc_clm_out_crr,pft_max_CLM,idx_mss_val,mss_val_in,sgs_out_nbr, & ! I
                  pft_rnk,rnk_vld_nbr) ! O
          end if ! endif lnd_msk_out
          
          ! Fill in pft_clm_out and pft_frc_clm with sgs_nbr_max PFTs 
          ! If CLM/LSM grid cell is ocean, set to no PFTs
          ! If LSM grid cell is land, there are three possibilities: 
          ! 1. If lnd_frc_out = 0, there is no PFT data from input grid 
          !    Since need land data, use bare ground
          ! 2. If lnd_frc_out > 0, there is PFT data from input grid but:
          !    a. use chosen PFT if it is not missing value
          !    b. missing value means no more PFTs with cover > 0
          
          ! If host model will compute land (non-ocean) processes on this gridpoint...
          if (lnd_msk_out(lon_out_idx,lat_out_idx) == 1) then
             ! ...then assign PFT and fractional coverage to each sub-grid patch...
             do sgs_idx=1,sgs_out_nbr
                ! ...if sgs_idx patch is unique, valid PFT...
                if (pft_rnk(sgs_idx) /= idx_mss_val) then
                   ! ...then assign next PFT from ranked list...
                   pft_clm_out(lon_out_idx,lat_out_idx,sgs_idx)=pft_rnk(sgs_idx) ! [enm] Plant functional type
                   ! ...and corresponding fractional coverage
                   pft_frc_clm(lon_out_idx,lat_out_idx,sgs_idx)= & ! [frc] Fraction covered by plant functional type
                        pft_frc_clm_out_crr(pft_rnk(sgs_idx))
                else ! ...else if sgs_idx"th patch does not exist, 
                   ! for example, whole cell is homogeneous so only one PFT exists, 
                   ! then asssign PFT for non-vegetated type and zero fraction
                   pft_clm_out(lon_out_idx,lat_out_idx,sgs_idx)=pft_non_vgt_CLM ! [enm] Plant functional type
                   pft_frc_clm(lon_out_idx,lat_out_idx,sgs_idx)=0.0 ! [frc] Fraction covered by plant functional type
                end if ! endif pft_rnk is valid
             end do ! end loop over sgs
          else ! ...else if host model has no land (non-ocean) processes on this gridpoint...
             ! ...then assign PFT and fractional coverage to non-vegetated and 0.0, respectively
             do sgs_idx=1,sgs_out_nbr
                pft_clm_out(lon_out_idx,lat_out_idx,sgs_idx)=0 ! [enm] Plant functional type
                pft_frc_clm(lon_out_idx,lat_out_idx,sgs_idx)=0.0 ! [frc] Fraction covered by plant functional type
             end do ! end loop over sgs
          end if ! endif not land
          
          ! Save initial plant function type distribution for re-normalization
          pft_frc_clm_crr(:)=pft_frc_clm(lon_out_idx,lat_out_idx,:) ! [frc] Initial PFT fraction of current gridcell
          frc_ttl_sgs=sum(pft_frc_clm_crr(:)) ! [frc] Sum of sub-gridscale fractions
          
          ! Re-normalize fractions on output grid to account for fraction occupied
          ! by PFTs not in sgs_out_nbr most important PFTs
          if (frc_ttl_sgs < frc_ttl_pft) then
             frc_dff=frc_ttl_pft-frc_ttl_sgs ! [frc] Excess PFT fraction to redistribute
             ! Sanity check
             if (frc_dff < 0.0) stop "frc_dff < 0.0 in pft_get()"
             do sgs_idx=1,sgs_out_nbr
                ! Weight redistributed fraction by existing fraction
                ! This can cause significant improvement in representation of 
                ! heterogeneous surfaces at low-resolution
                ! Note CLM redistributes excess evenly among sgs_out_nbr patches
                ! Thus this method produces differences with CLM PFT generator
                pft_frc_clm(lon_out_idx,lat_out_idx,sgs_idx)= &
                     pft_frc_clm(lon_out_idx,lat_out_idx,sgs_idx)+ &
                     pft_frc_clm_crr(sgs_idx)*frc_dff/frc_ttl_sgs
             end do ! end loop over sgs
             frc_ttl_sgs=sum(pft_frc_clm(lon_out_idx,lat_out_idx,:)) ! [frc] Total fraction of gridcell accounted for by all sub-gridscale patches
             
          end if ! endif re-normalizing
          
          ! fxm: Renormalize by fractional land in gridcell
          ! This differs from re-normalization for area in unused PFTs done above
          ! CLM has pft_frc always sum to unity, but weights contribution by lnd_frc
          ! Need to understand implications of downstream weighting better
          ! Better fix might be to write new ovr_wgt_get() routine
          if ((lnd_msk_out(lon_out_idx,lat_out_idx)==1).and. &
               (abs(lnd_frc_out(lon_out_idx,lat_out_idx)-1.0) > eps_rlt)) then
             pft_frc_clm(lon_out_idx,lat_out_idx,:)= &
                  pft_frc_clm(lon_out_idx,lat_out_idx,:)/frc_ttl_sgs
          end if ! endif re-normalizing
          
       end do ! end loop over lon
    end do ! end loop over lat
    
    ! Sanity check
    err_nbr=0
    do lat_out_idx=1,lat_out_nbr
       do lon_out_idx=1,lon_out_nbr
          frc_ttl=sum(pft_frc_clm(lon_out_idx,lat_out_idx,:))
          pft_max=maxval(pft_clm_out(lon_out_idx,lat_out_idx,:)) ! [enm] Maximum plant functional type
          pft_min=minval(pft_clm_out(lon_out_idx,lat_out_idx,:)) ! [enm] Minimum plant functional type
          ! fxm: Really we only care about regions with lnd_msk == 1
          ! Code depending on other gridcells may break in coupled or regional runs
          ! Remove checks on lnd_msk != 0 once have verified that dust model 
          ! uses consistent land mask for all input fields
          ! Until then, make sure values where lnd_msk == 0 or 1 are reasonable
          ! ...Flag error and record point for diagnostics if...
          if ( &
               ! ...vegetation fraction is completely unphysical or...
               (frc_ttl < 0.0).or. &
               ! ...land fraction is completely unphysical or...
               (lnd_frc_out(lon_out_idx,lat_out_idx) < 0.0).or. &
               ! ...valid land points have one of three things wrong...
               ((lnd_msk_out(lon_out_idx,lat_out_idx)==1).and. & ! NB: _AND_
               ! ...1. Invalid PFT...
               ((pft_min < 0.or.pft_max > pft_max_CLM).or. &
               ! ...2. Non-unity total fraction occupied by PFTs...
               (abs(frc_ttl-1.0) > eps_rlt).or. &
               ! ...3. Land fraction measurably exceeds unity...
               (lnd_frc_out(lon_out_idx,lat_out_idx)-1.0 > eps_rlt)))) then
             err_nbr=err_nbr+1
             lon_out_idx_err=lon_out_idx
             lat_out_idx_err=lat_out_idx
          end if ! endif err
       end do ! end loop over lon
    end do ! end loop over lat
    if (err_nbr > 0) then
       frc_ttl=sum(pft_frc_clm(lon_out_idx_err,lat_out_idx_err,:))
       write (6,"(a,i6,a)") "ERROR "//sbr_nm//" reports ",err_nbr," errors during rebinning"
       write (6,"(a,i2)") "Output land mask = ",lnd_msk_out(lon_out_idx_err,lat_out_idx_err)
       write (6,"(a,f15.12)") "Output land fraction = ",lnd_frc_out(lon_out_idx_err,lat_out_idx_err)
       write (6,"(a,f15.12)") "Plant functional type fractions sum to ",frc_ttl
       write (6,"(a1,1x,a3,1x,a35,1x,a)") "#","PFT","Description                       ","Fraction"
       do sgs_idx=1,sgs_out_nbr
          write (6,"(i1,1x,i3,1x,a35,1x,f9.6)") sgs_idx, &
               pft_clm_out(lon_out_idx_err,lat_out_idx_err,sgs_idx), &
               pft_sng(pft_clm_out(lon_out_idx_err,lat_out_idx_err,sgs_idx)), &
               pft_frc_clm(lon_out_idx_err,lat_out_idx_err,sgs_idx)
       end do ! end loop over sgs
       write (6,"(a)") "Output cell edge locations:"
       write (6,"(a)") "(idx)   Lat Sth (idx)   Lat Nrt (idx)   Lon Wst (idx)   Lon Est"
       write (6,"(4(a1,i3,a1,1x,f9.4,1x))") &
            "(",lat_out_idx_err,")",lat_out_grd(lat_out_idx_err),"(",lat_out_idx_err+1,")",lat_out_grd(lat_out_idx_err+1), &
            "(",lon_out_idx_err,")",lon_out_grd(lon_out_idx_err),"(",lon_out_idx_err+1,")",lon_out_grd(lon_out_idx_err+1)
       stop
    endif ! endif err
    
    ! De-allocate
    if (allocated(area_in)) deallocate(area_in,stat=rcd)
    if(rcd /= 0) stop "deallocate() failed for area_in"
    if (allocated(lat_in_grd)) deallocate(lat_in_grd,stat=rcd)
    if(rcd /= 0) stop "deallocate() failed for lat_in_grd"
    if (allocated(lat_in_2d)) deallocate(lat_in_2d,stat=rcd)
    if(rcd /= 0) stop "deallocate() failed for lat_in_2d"
    if (allocated(lon_in_grd)) deallocate(lon_in_grd,stat=rcd)
    if(rcd /= 0) stop "deallocate() failed for lon_in_grd"
    if (allocated(lon_in_grd_2d)) deallocate(lon_in_grd_2d,stat=rcd)
    if(rcd /= 0) stop "deallocate() failed for lon_in_grd_2d"
    if (allocated(lon_in_2d)) deallocate(lon_in_2d,stat=rcd)
    if(rcd /= 0) stop "deallocate() failed for lon_in_2d"
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
    
    if (allocated(pft_frc_clm_in)) deallocate(pft_frc_clm_in,stat=rcd)
    if(rcd /= 0) stop "deallocate() failed for pft_frc_clm_in"
    if (allocated(pft_frc_clm_out)) deallocate(pft_frc_clm_out,stat=rcd)
    if(rcd /= 0) stop "deallocate() failed for pft_frc_clm_out"
    if (allocated(pft_frc_clm_out_crr)) deallocate(pft_frc_clm_out_crr,stat=rcd)
    if(rcd /= 0) stop "deallocate() failed for pft_frc_clm_out_crr"
    
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Exiting "//sbr_nm
    return 
  end subroutine pft_get                       ! end pft_get()
  
  subroutine lai_get(fl_in, & ! I
       lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr,sgs_out_nbr,time_out_nbr, & ! I
       area_out,lnd_msk_out,pft_clm_out,pft_frc_clm_out, & ! I
       lai_clm,lai_ttl_clm,vai_clm,vai_ttl_clm) ! O
    ! Purpose: Retrieve leaf area index and vegetation area index and rebin to requested grid
    ! lai_get() is called by bds_prc()
    ! Leaf area index "raw" dataset are rebinned from LAI dataset assembled by ....
    ! Dataset resolution is 0.5 x 0.5 degree and data are written from South pole to North pole, West to East, beginning at date line
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use map_grd ! [mdl] Map grids and regridding
    use netcdf ! [mdl] netCDF interface
    use nf90_utl ! [mdl] netCDF utilities
    use sng_mdl,only:ftn_strlen ! [mdl] String manipulation
    use utl_mdl,only:mnt_ncr_chk ! [mdl] Utility functions (date_time_get,mnt_chk...)
    implicit none
    ! Parameters
    character(len=*),parameter::sbr_nm="lai_get()" ! [sng] Subroutine name
    integer,parameter::lat_idx_dbg=50 ! [idx] Longitude for verbose output
    integer,parameter::lon_idx_dbg=90 ! [idx] Longitude for verbose output
    integer,parameter::idx_mss_val=99999 ! [idx] Index denoting missing/invalid data
    integer,parameter::pft_nbr_CLM=17 ! [nbr] Number of plant functional types in CLM 
    integer,parameter::pft_non_vgt_CLM=0 ! [enm] Non-vegetated plant functional type in CLM
    integer,parameter::sgs_nbr_CLM=4 ! [nbr] Number of sub-gridscale patches in CLM 
    integer,parameter::time_nbr_CLM=12 ! [nbr] Number of times in CLM 
    integer,parameter::unit_dgn=74 ! Unit for writing general Map diagnoses
    real,parameter::eps_rlt=1.0e-5 ! [frc] Relative error allowed in frc_ttl
    real,parameter::lai_max_CLM=100.0 ! [m2 m-2] Maximum leaf area index, one-sided
    ! Derived parameters
    integer,parameter::pft_max_CLM=pft_nbr_CLM-1 ! [enm] Maximum value of PFT in CLM
    ! Commons
    ! Input
    character,intent(in)::fl_in*80 ! [sng] Input file
    integer,intent(in)::lat_out_nbr ! [nbr] Number of latitudes
    integer,intent(in)::lon_out_nbr ! [nbr] Number of longitudes
    integer,intent(in)::sgs_out_nbr ! [nbr] Number of sub-gridscale patches per gridcell
    integer,intent(in)::time_out_nbr ! [nbr] Number of times
    integer,intent(in)::lnd_msk_out(lon_out_nbr,lat_out_nbr) ! [msk] Land mask (integer 0 or 1)
    integer,intent(in)::pft_clm_out(lon_out_nbr,lat_out_nbr,sgs_out_nbr) ! [enm] Plant functional type
    real,intent(in)::area_out(lon_out_nbr,lat_out_nbr) ! [m2] Area of gridcells
    real,intent(in)::lat_out_grd(lat_out_nbr+1) ! [dgr] Interface latitudes
    real,intent(in)::lon_out_grd(lon_out_nbr+1) ! [dgr] Interface longitudes
    real,intent(in)::pft_frc_clm_out(lon_out_nbr,lat_out_nbr,sgs_out_nbr) ! [enm] Fraction covered by plant functional type
    
    ! Output
    real,intent(out)::lai_clm(lon_out_nbr,lat_out_nbr,sgs_out_nbr,time_out_nbr) ! [m2 m-2] Leaf area index, one-sided
    real,intent(out)::lai_ttl_clm(lon_out_nbr,lat_out_nbr,time_out_nbr) ! [m2 m-2] Total leaf area index, one-sided
    real,intent(out)::vai_clm(lon_out_nbr,lat_out_nbr,sgs_out_nbr,time_out_nbr) ! [m2 m-2] Vegetation area index, one-sided
    real,intent(out)::vai_ttl_clm(lon_out_nbr,lat_out_nbr,time_out_nbr) ! [m2 m-2] Total vegetation area index, one-sided
    
    ! Locals with simple initialization
    character(80)::fl_dgn="/tmp/zender/map/lai.txt"//char(0) ! [sng] Diagnostics file
    character(35) pft_sng(0:pft_max_CLM) ! [sng] Description of plant functional types
    data pft_sng( 0) /"not vegetated"                      /
    data pft_sng( 1) /"needleleaf evergreen temperate tree"/
    data pft_sng( 2) /"needleleaf evergreen boreal tree"   /
    data pft_sng( 3) /"needleleaf deciduous boreal tree"   /
    data pft_sng( 4) /"broadleaf evergreen tropical tree"  /
    data pft_sng( 5) /"broadleaf evergreen temperate tree" /
    data pft_sng( 6) /"broadleaf deciduous tropical tree"  /
    data pft_sng( 7) /"broadleaf deciduous temperate tree" /
    data pft_sng( 8) /"broadleaf deciduous boreal tree"    /
    data pft_sng( 9) /"broadleaf evergreen shrub"          /
    data pft_sng(10) /"broadleaf deciduous temperate shrub"/
    data pft_sng(11) /"broadleaf deciduous boreal shrub"   /
    data pft_sng(12) /"c3 arctic grass"                    /
    data pft_sng(13) /"c3 non-arctic grass"                /
    data pft_sng(14) /"c4 grass"                           / 
    data pft_sng(15) /"corn"                               /
    data pft_sng(16) /"wheat"                              /
    integer::rcd=nf90_noerr ! [enm] Return success code
    real::mss_val_in=1.0e36 ! [frc] Missing value
    
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
    integer pft_idx ! [idx] Counting index for pft
    integer pft_in_idx ! [idx] Counting index for pft
    integer pft_nbr ! [nbr] Number of plant functional types
    integer sgs_idx ! [idx] Counting index for sgs
    integer time_idx ! [idx] Counting index for time
    integer time_nbr ! [nbr] Number of times
    logical lat_mnt_ncr ! [flg] Latitude monotonically decreases
    logical mnt_ncr ! [flg] Monotonic and increasing flag
    logical mss_flg ! [flg] Variable has missing_value attribute
    real edg_grd(4) ! [dgr] Grid edges (north, east, south, west)
    real frc_ttl ! [frc] Total fraction of gridcell accounted for by all plant functional types
    real lnd_frc_out(lon_out_nbr,lat_out_nbr) ! [frc] Fraction of land (not ocean)
    real ovr_wgt_crr ! [frc] Overlap weight of current gridcell
    
    ! Variables needed to read external netCDF dataset
    integer cnt_lon_lat_pft_time(4) ! [nbr] Dimension sizes
    integer edg_est_id ! [id] Variable ID
    integer edg_nrt_id ! [id] Variable ID
    integer edg_sth_id ! [id] Variable ID
    integer edg_wst_id ! [id] Variable ID
    integer lai_clm_dmn_nbr ! [nbr] Number of dimensions in disk lai_clm field
    integer lai_clm_id ! [id] Variable ID
    integer lat_dmn_id ! [id] Dimension ID for lat
    integer lat_in_2d_id ! [id] Variable ID
    integer lat_in_id ! [id] Variable ID
    integer lnd_frc_clm_id ! [id] Variable ID
    integer lon_dmn_id ! [id] Dimension ID for lon
    integer lon_in_2d_id ! [id] Variable ID
    integer mss_val_id ! [id] Attribute ID
    integer nc_id ! [id] File handle
    integer pft_dmn_id ! [id] Dimension ID for pft
    integer srt_lon_lat_pft_time(4) ! [idx] Dimension offset indices
    integer time_dmn_id ! [id] Dimension ID for time
    integer vai_clm_id ! [id] Variable ID
    
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
    
    real,dimension(:,:),allocatable::lnd_frc_in ! [frc] Fraction of land (not ocean)
    real,dimension(:,:),allocatable::lnd_msk_in ! [msk] Land mask (0.0 or 1.0)
    real,dimension(:,:,:),allocatable::lai_clm_in ! [m2 m-2] Leaf area index, one-sided
    real,dimension(:,:,:),allocatable::vai_clm_in ! [m2 m-2] Vegetation area index, one-sided
    real,dimension(:,:,:),allocatable::lai_clm_out ! [m2 m-2] Leaf area index, one-sided
    real,dimension(:,:,:),allocatable::vai_clm_out ! [m2 m-2] Vegetation area index, one-sided
    
    ! Main Code
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Entering "//sbr_nm
    
    ! Check for monotonically increasing grids
    mnt_ncr=mnt_ncr_chk(lon_out_grd,lon_out_nbr+1)
    if (.not.mnt_ncr) stop "lon_out_grd not monotonically increasing in lai_get()"
    mnt_ncr=mnt_ncr_chk(lat_out_grd,lat_out_nbr+1)
    if (.not.mnt_ncr) stop "lat_out_grd not monotonically increasing in lai_get()"
    ! Sanity check
    if (sgs_out_nbr > pft_nbr_CLM) stop "sgs_out_nbr > pft_nbr_CLM"
    if (time_out_nbr /= time_nbr_CLM) stop "time_out_nbr /= time_nbr_CLM"
    
    ! Initialize arrays
    lai_clm(:,:,:,:)=0.0 ! [m2 m-2] Leaf area index, one-sided
    lai_ttl_clm(:,:,:)=0.0 ! [m2 m-2] Total leaf area index, one-sided
    vai_clm(:,:,:,:)=0.0 ! [m2 m-2] Vegetation area index, one-sided
    vai_ttl_clm(:,:,:)=0.0 ! [m2 m-2] Total vegetation area index, one-sided
    lnd_frc_out(:,:)=0.0 ! [frc] Fraction of land (not ocean)
    
    ! Read netCDF data
    rcd=rcd+nf90_wrp_open(fl_in,nf90_nowrite,nc_id)
    ! Get dimension IDs
    rcd=nf90_wrp_inq_dimid(nc_id,"lon",lon_dmn_id)
    rcd=nf90_wrp_inq_dimid(nc_id,"lat",lat_dmn_id)
    rcd=nf90_wrp_inq_dimid(nc_id,"pft",pft_dmn_id)
    rcd=nf90_wrp_inq_dimid(nc_id,"time",time_dmn_id)
    if(rcd /= nf90_noerr) stop "Error retrieving dimension IDs"
    ! Get dimension sizes
    rcd=rcd+nf90_inquire_dimension(nc_id,lon_dmn_id,len=lon_in_nbr)
    rcd=rcd+nf90_inquire_dimension(nc_id,lat_dmn_id,len=lat_in_nbr)
    rcd=rcd+nf90_inquire_dimension(nc_id,pft_dmn_id,len=pft_nbr)
    rcd=rcd+nf90_inquire_dimension(nc_id,time_dmn_id,len=time_nbr)
    if(rcd /= nf90_noerr) stop "Error retrieving dimension sizes"
    ! Enough memory? 
    if (lat_in_nbr /= 360) stop "lat_in_nbr /= 360 in lai_get()"
    if (lon_in_nbr /= 720) stop "lon_in_nbr /= 720 in lai_get()"
    if (pft_nbr /= pft_nbr_CLM) stop "pft_nbr /= pft_nbr_CLM in lai_get()"
    if (time_nbr /= time_nbr_CLM) stop "time_nbr /= time_nbr_CLM in lai_get()"
    ! Allocate based on input file dimensions
    allocate(area_in(lon_in_nbr,lat_in_nbr),stat=rcd) ! [m2] Area of gridcells
    if(rcd /= 0) stop "allocate() failed for area_in"
    allocate(lat_in_grd(lat_in_nbr+1),stat=rcd) ! [dgr] Interface latitudes
    if(rcd /= 0) stop "allocate() failed for lat_in_grd"
    allocate(lat_in_2d(lon_in_nbr,lat_in_nbr),stat=rcd) ! [dgr] Latitude at gridcell center
    if(rcd /= 0) stop "allocate() failed for lat_in_2d"
    allocate(lon_in_grd(lon_in_nbr+1),stat=rcd) ! [dgr] Interface longitudes
    if(rcd /= 0) stop "allocate() failed for lon_in_grd"
    allocate(lon_in_grd_2d(lon_in_nbr+1,lat_in_nbr),stat=rcd) ! [dgr] Interface longitudes
    if(rcd /= 0) stop "allocate() failed for lon_in_grd_2d"
    allocate(lon_in_2d(lon_in_nbr,lat_in_nbr),stat=rcd) ! [dgr] Longitude at gridcell center
    if(rcd /= 0) stop "allocate() failed for lon_in_2d"
    allocate(lon_in_nbr_1d(lat_in_nbr),stat=rcd) ! [nbr] Number of longitudes per latitude
    if(rcd /= 0) stop "allocate() failed for lon_in_nbr_1d"
    allocate(lnd_frc_in(lon_in_nbr,lat_in_nbr),stat=rcd) ! [frc] Fraction of land (not ocean)
    if(rcd /= 0) stop "allocate() failed for lnd_frc_in"
    allocate(lnd_msk_in(lon_in_nbr,lat_in_nbr),stat=rcd) ! [msk] Land mask (0.0 or 1.0)
    if(rcd /= 0) stop "allocate() failed for lnd_msk_in"
    
    allocate(lai_clm_out(lon_out_nbr,lat_out_nbr,0:pft_nbr-1),stat=rcd) ! [m2 m-2] Leaf area index, one-sided
    if(rcd /= 0) stop "allocate() failed for lai_clm_out"
    allocate(vai_clm_out(lon_out_nbr,lat_out_nbr,0:pft_nbr-1),stat=rcd) ! [m2 m-2] Vegetation area index, one-sided
    if(rcd /= 0) stop "allocate() failed for vai_clm_out"
    allocate(lai_clm_in(lon_in_nbr,lat_in_nbr,0:pft_nbr-1),stat=rcd) ! [m2 m-2] Leaf area index, one-sided
    if(rcd /= 0) stop "allocate() failed for lai_clm_in"
    allocate(vai_clm_in(lon_in_nbr,lat_in_nbr,0:pft_nbr-1),stat=rcd) ! [m2 m-2] Vegetation area index, one-sided
    if(rcd /= 0) stop "allocate() failed for vai_clm_in"
    
    ! Get variable IDs
    rcd=nf90_wrp_inq_varid(nc_id,"LATIXY",lat_in_2d_id)
    rcd=nf90_wrp_inq_varid(nc_id,"LONGXY",lon_in_2d_id)
    rcd=nf90_wrp_inq_varid(nc_id,"EDGEN",edg_nrt_id)
    rcd=nf90_wrp_inq_varid(nc_id,"EDGEE",edg_est_id)
    rcd=nf90_wrp_inq_varid(nc_id,"EDGES",edg_sth_id)
    rcd=nf90_wrp_inq_varid(nc_id,"EDGEW",edg_wst_id)
    rcd=nf90_wrp_inq_varid(nc_id,"MONTHLY_LAI",lai_clm_id)
    rcd=nf90_wrp_inq_varid(nc_id,"MONTHLY_SAI",vai_clm_id) ! NB: Read SAI, then add LAI
    rcd=nf90_wrp_inq_varid(nc_id,"LANDMASK",lnd_frc_clm_id)
    
    ! Get number of dimensions
    rcd=rcd+nf90_inquire_variable(nc_id,lai_clm_id,ndims=lai_clm_dmn_nbr)
    if (lai_clm_dmn_nbr /= 4) then
       write (6,"(2a,i1,a)") prg_nm(1:ftn_strlen(prg_nm)), & 
            ": ERROR MONTHLY_LAI has ",lai_clm_dmn_nbr," dimensions"
       write (6,"(2a,i1,a)") prg_nm(1:ftn_strlen(prg_nm)), & 
            ": HINT Make sure MONTHLY_LAI has phenotype and time dimensions"
       stop
    endif                  ! endif err
    ! Get time-independent data
    rcd=nf90_wrp(nf90_get_var(nc_id,lat_in_2d_id,lat_in_2d),"get_var lat_in_2d")
    rcd=nf90_wrp(nf90_get_var(nc_id,lon_in_2d_id,lon_in_2d),"get_var lon_in_2d")
    rcd=nf90_wrp(nf90_get_var(nc_id,edg_nrt_id,edg_grd(1)),"get_var edg_nrt")
    rcd=nf90_wrp(nf90_get_var(nc_id,edg_est_id,edg_grd(2)),"get_var edg_est")
    rcd=nf90_wrp(nf90_get_var(nc_id,edg_sth_id,edg_grd(3)),"get_var edg_sth")
    rcd=nf90_wrp(nf90_get_var(nc_id,edg_wst_id,edg_grd(4)),"get_var edg_wst")
    rcd=nf90_wrp(nf90_get_var(nc_id,lnd_frc_clm_id,lnd_frc_in),"get_var lnd_frc_in")
    if (rcd /= nf90_noerr) stop "Error in nf90_get_var in lai_get()"
    ! Get missing value
    rcd=rcd+nf90_inquire_attribute(nc_id,lai_clm_id,"missing_value",attnum=mss_val_id)
    if (rcd==nf90_noerr) then
       mss_flg=.true. ! [flg] Data may have missing values
       rcd=rcd+nf90_get_att(nc_id,lai_clm_id,"missing_value",mss_val_in)
       if (mss_val_in /= 1.0e36) write (6,"(2a,f12.6)") prg_nm(1:ftn_strlen(prg_nm)), &
            ": WARNING "//sbr_nm//" reports mss_val_in = ",mss_val_in
    else
       mss_flg=.false. ! [flg] Data may have missing values
       rcd=nf90_noerr ! [enm] Return success code
    endif                  ! endif
    
    ! Determine map grid and overlap characteristics
    ! fxm: Method breaks for reduced grids, streamline, subroutinize, and generalize
    lon_in_nbr_1d(:)=lon_in_nbr ! [nbr] Number of longitudes per latitude
    
    ! Convert gridpoint centers, domain boundaries to grid interfaces
    call map_edge_mk(lat_in_nbr,lon_in_nbr,lon_in_nbr_1d, & ! I
         lat_in_2d,lon_in_2d,edg_grd, & ! I 
         lat_in_grd,lon_in_grd,lon_in_grd_2d) ! O
    
    mnt_ncr=mnt_ncr_chk(lat_in_grd,lat_in_nbr+1)
    if (.not.mnt_ncr) stop "ERROR: lat_in_grd not monotonically increasing in lai_get()"
    mnt_ncr=mnt_ncr_chk(lon_in_grd,lon_in_nbr+1)
    if (.not.mnt_ncr) stop "ERROR: lon_in_grd not monotonically increasing in lai_get()"
    
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
    
    ! Land mask is unity for LSM/CLM grid cells and 0 for ocean points
    ! Use land mask in array as weights compute land fraction on output grid
    ! For this purpose, temporarily set input land mask to unity everywhere
    lnd_msk_in(:,:)=1.0 ! [msk] Land mask (0.0 or 1.0)
    
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
    
    ! lnd_msk_out is land mask of host model, two possible integer values (0, 1):
    ! 1 if there is ANY land in output gridcell (then apply lnd_frc_out to fluxes)
    ! 0 if there is NO  land in output gridcell (do not compute land fluxes)
    
    ! Land mask is unity for LSM/CLM grid cells and 0 for ocean points
    ! lnd_frc_in is valid mask because it is 1.0 or 0.0 on input grid
    ! This is only true when input file is "raw data" high resolution grid
    lnd_msk_in(:,:)=lnd_frc_in(:,:) ! [msk] Land mask (0.0 or 1.0)
    
    ! fxm: I"m not sure about purpose of second areaini() in mklai.F90
    ! so I"ve skipped it temporarily
#if 0
    ! Get overlap locations and weights for mapping input grid to output grid
    ! Use lnd_frc_in as input mask
    call map_ovr_wgt_drv( &
         lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
         lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
         ovr_nbr_max,area_out, & ! I
         ovr_lat_idx,ovr_lon_idx,ovr_nbr,ovr_wgt) ! O
#endif /* endif 0 */
    
    ! Outer loop over months reduces memory overhead
    do time_idx=1,time_nbr
       
       ! Initialize hyperslab indices for disk file
       srt_lon_lat_pft_time=(/1,1,1,time_idx/) ! [idx] Dimension offset indices
       cnt_lon_lat_pft_time=(/lon_in_nbr,lat_in_nbr,pft_nbr,1/) ! [nbr] Dimension sizes
       ! Get time-dependent data
       rcd=nf90_wrp(nf90_get_var(nc_id,lai_clm_id,lai_clm_in, &
            start=srt_lon_lat_pft_time,count=cnt_lon_lat_pft_time),"get_var lai_clm_in")
       rcd=nf90_wrp(nf90_get_var(nc_id,vai_clm_id,vai_clm_in, &
            start=srt_lon_lat_pft_time,count=cnt_lon_lat_pft_time),"get_var vai_clm_in")
       
       ! VAI (vai_clm_in) is actually stem area index (SAI) when input
       ! Convert to actual VAI by adding LAI (VAI = SAI + LAI)
       do pft_in_idx=0,pft_max_CLM ! NB: Loop starts at zero
          do lat_in_idx=1,lat_in_nbr
             do lon_in_idx=1,lon_in_nbr
                if (mss_flg) then
                   if (lai_clm_in(lon_in_idx,lat_in_idx,pft_in_idx)/=mss_val_in) then 
                      vai_clm_in(lon_in_idx,lat_in_idx,pft_in_idx)= & ! [m2 m-2] Vegetation area index, one-sided
                           vai_clm_in(lon_in_idx,lat_in_idx,pft_in_idx)+ & 
                           lai_clm_in(lon_in_idx,lat_in_idx,pft_in_idx)
                   endif ! endif mss_val
                else ! not mss_flg
                   vai_clm_in(lon_in_idx,lat_in_idx,pft_in_idx)= & ! [m2 m-2] Vegetation area index, one-sided
                        vai_clm_in(lon_in_idx,lat_in_idx,pft_in_idx)+ & 
                        lai_clm_in(lon_in_idx,lat_in_idx,pft_in_idx)
                end if ! end if mss_flg
             end do ! end loop over lon
          end do ! end loop over lat
       end do ! end loop over pft
       
       ! Sanity check
       err_nbr=0
       do lat_in_idx=1,lat_in_nbr
          do lon_in_idx=1,lon_in_nbr
             frc_ttl=maxval(lai_clm_in(lon_in_idx,lat_in_idx,:))
             ! ...Flag error and record point for diagnostics if...
             if ( &
                  ! ...land fraction is completely unphysical or...
                  (lnd_frc_in(lon_in_idx,lat_in_idx) < 0.0.or.lnd_frc_in(lon_in_idx,lat_in_idx) > 1.0).or. &
                  ! ...input land points have one of three things wrong...
                  ((lnd_frc_in(lon_in_idx,lat_in_idx)==1.0).and. & ! NB: _AND_
                  ! ...1. Unphysical maximum LAI
                  (frc_ttl < 0.0.or.frc_ttl > lai_max_CLM).or. &
                  ! ...2. Non-vegetated surface has LAI or VAI
                  (lai_clm_in(lon_in_idx,lat_in_idx,pft_non_vgt_CLM) /= 0.0.or. &
                  vai_clm_in(lon_in_idx,lat_in_idx,pft_non_vgt_CLM) /= 0.0))) then
                err_nbr=err_nbr+1
                lon_in_idx_err=lon_in_idx
                lat_in_idx_err=lat_in_idx
             end if ! endif err
          end do ! end loop over lon
       end do ! end loop over lat
       if (err_nbr > 0) then
          frc_ttl=maxval(lai_clm_in(lon_in_idx_err,lat_in_idx_err,:))
          write (6,"(a,i6,a)") "ERROR "//sbr_nm//" reports ",err_nbr," errors during input"
          write (6,"(a,i2)") "Time index (month) = ",time_idx
          write (6,"(a,f15.12)") "Input land fraction = ",lnd_frc_in(lon_in_idx_err,lat_in_idx_err)
          write (6,"(a,f15.12)") "Maximum LAI = ",frc_ttl
          write (6,"(a3,1x,a35,1x,2(a9,1x))") "PFT","Description                       ","LAI","VAI"
          do pft_idx=0,pft_max_CLM ! NB: Loop starts at zero
             write (6,"(i3,1x,a35,1x,2(f9.6,1x))") pft_idx,pft_sng(pft_idx), &
                  lai_clm_in(lon_in_idx_err,lat_in_idx_err,pft_idx), &
                  vai_clm_in(lon_in_idx_err,lat_in_idx_err,pft_idx)
          end do ! end loop over pft
          write (6,"(a)") "Input cell edge locations:"
          write (6,"(a)") "(idx)   Lat Sth (idx)   Lat Nrt (idx)   Lon Wst (idx)   Lon Est"
          write (6,"(4(a1,i3,a1,1x,f9.4,1x))") &
               "(",lat_in_idx_err,")",lat_in_grd(lat_in_idx_err),"(",lat_in_idx_err+1,")",lat_in_grd(lat_in_idx_err+1), &
               "(",lon_in_idx_err,")",lon_in_grd(lon_in_idx_err),"(",lon_in_idx_err+1,")",lon_in_grd(lon_in_idx_err+1)
          stop
       endif ! endif err
       
       do pft_idx=0,pft_max_CLM ! NB: Loop starts at zero
          ! Rebin leaf area index from input grid to output grid
          call map_rbn(lai_clm_in(1,1,pft_idx),mss_flg,mss_val_in, & ! I
               lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
               lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
               ovr_nbr_max,ovr_lat_idx,ovr_lon_idx,ovr_nbr,ovr_wgt, & ! I
               lai_clm_out(1,1,pft_idx)) ! O
          
          ! Rebin vegetation area index from input grid to output grid
          call map_rbn(vai_clm_in(1,1,pft_idx),mss_flg,mss_val_in, & ! I
               lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
               lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
               ovr_nbr_max,ovr_lat_idx,ovr_lon_idx,ovr_nbr,ovr_wgt, & ! I
               vai_clm_out(1,1,pft_idx)) ! O
       end do ! end loop over pft
       
       ! Post-process data
       ! Set input LAI and VAI to 0.0 for all PFTs where host model has no land
       do lat_out_idx=1,lat_out_nbr
          do lon_out_idx=1,lon_out_nbr
             if (lnd_msk_out(lon_out_idx,lat_out_idx)==0) then
                lai_clm_out(lon_out_idx,lat_out_idx,:)=0.0
                vai_clm_out(lon_out_idx,lat_out_idx,:)=0.0
             end if ! end if lnd
          end do ! end loop over lon
       end do ! end loop over lat
       
       ! Set output LAI/VAI for each SGS to input LAI/VAI of corresponding PFT
       do lat_out_idx=1,lat_out_nbr
          do lon_out_idx=1,lon_out_nbr
             do pft_idx=0,pft_max_CLM ! NB: Loop starts at zero
                do sgs_idx=1,sgs_out_nbr
                   if (pft_clm_out(lon_out_idx,lat_out_idx,sgs_idx) == pft_idx) then
                      lai_clm(lon_out_idx,lat_out_idx,sgs_idx,time_idx)= &
                           lai_clm_out(lon_out_idx,lat_out_idx,pft_idx)
                      vai_clm(lon_out_idx,lat_out_idx,sgs_idx,time_idx)= &
                           vai_clm_out(lon_out_idx,lat_out_idx,pft_idx)
                   endif ! endif
                end do ! end loop over sgs
             end do ! end loop over pft
          end do ! end loop over lon
       end do ! end loop over lat
       
       ! Diagnose total monthly LAI/VAI by weighting sub-gridscale contributions
       do lat_out_idx=1,lat_out_nbr
          do lon_out_idx=1,lon_out_nbr
             do sgs_idx=1,sgs_out_nbr
                lai_ttl_clm(lon_out_idx,lat_out_idx,time_idx)= & ! [m2 m-2] Total leaf area index, one-sided
                     lai_ttl_clm(lon_out_idx,lat_out_idx,time_idx)+ &
                     pft_frc_clm_out(lon_out_idx,lat_out_idx,sgs_idx)* &
                     lai_clm(lon_out_idx,lat_out_idx,sgs_idx,time_idx)
                vai_ttl_clm(lon_out_idx,lat_out_idx,time_idx)= & ! [m2 m-2] Total vegetation area index, one-sided
                     vai_ttl_clm(lon_out_idx,lat_out_idx,time_idx)+ &
                     pft_frc_clm_out(lon_out_idx,lat_out_idx,sgs_idx)* &
                     vai_clm(lon_out_idx,lat_out_idx,sgs_idx,time_idx)
             end do ! end loop over sgs
          end do ! end loop over lon
       end do ! end loop over lat
       
       ! Sanity check
       err_nbr=0
       do lat_out_idx=1,lat_out_nbr
          do lon_out_idx=1,lon_out_nbr
             frc_ttl=maxval(lai_clm_out(lon_out_idx,lat_out_idx,:))
             ! ...Flag error and record point for diagnostics if...
             if ( &
                  ! ...land fraction is completely unphysical or...
                  (lnd_frc_out(lon_out_idx,lat_out_idx) < 0.0).or. &
                  ! ...valid land points have one of two things wrong...
                  ((lnd_msk_out(lon_out_idx,lat_out_idx)==1).and. & ! NB: _AND_
                  ! ...1. Unphysical maximum LAI
                  ((frc_ttl < 0.0.or.frc_ttl > lai_max_CLM).or. &
                  ! ...2. Land fraction measurably exceeds unity...
                  (lnd_frc_out(lon_out_idx,lat_out_idx)-1.0 > eps_rlt)))) then
                err_nbr=err_nbr+1
                lon_out_idx_err=lon_out_idx
                lat_out_idx_err=lat_out_idx
             end if ! endif err
          end do ! end loop over lon
       end do ! end loop over lat
       if (err_nbr > 0) then
          frc_ttl=maxval(lai_clm_out(lon_out_idx_err,lat_out_idx_err,:))
          write (6,"(a,i6,a)") "ERROR "//sbr_nm//" reports ",err_nbr," errors during rebinning"
          write (6,"(a,i2)") "Time index (month) = ",time_idx
          write (6,"(a,i2)") "Output land mask = ",lnd_msk_out(lon_out_idx_err,lat_out_idx_err)
          write (6,"(a,f15.12)") "Output land fraction = ",lnd_frc_out(lon_out_idx_err,lat_out_idx_err)
          write (6,"(a,f15.12)") "Maximum LAI = ",frc_ttl
          write (6,"(a1,1x,a3,1x,a35,1x,3(a9,1x))") "#","PFT","Description                       ","Fraction","LAI","VAI"
          do sgs_idx=1,sgs_out_nbr
             write (6,"(i1,1x,i3,1x,a35,1x,3(f9.6,1x))") sgs_idx, &
                  pft_clm_out(lon_out_idx_err,lat_out_idx_err,sgs_idx), &
                  pft_sng(pft_clm_out(lon_out_idx_err,lat_out_idx_err,sgs_idx)), &
                  pft_frc_clm_out(lon_out_idx_err,lat_out_idx_err,sgs_idx), &
                  lai_clm_out(lon_out_idx_err,lat_out_idx_err,sgs_idx), &
                  vai_clm_out(lon_out_idx_err,lat_out_idx_err,sgs_idx)
          end do ! end loop over sgs
          write (6,"(a3,1x,a35,1x,2(a9,1x))") "PFT","Description                       ","LAI","VAI"
          do pft_idx=0,pft_max_CLM ! NB: Loop starts at zero
             write (6,"(i3,1x,a35,1x,2(f9.6,1x))") pft_idx,pft_sng(pft_idx), &
                  lai_clm_out(lon_out_idx_err,lat_out_idx_err,pft_idx), &
                  vai_clm_out(lon_out_idx_err,lat_out_idx_err,pft_idx)
          end do ! end loop over pft
          write (6,"(a)") "Output cell edge locations:"
          write (6,"(a)") "(idx)   Lat Sth (idx)   Lat Nrt (idx)   Lon Wst (idx)   Lon Est"
          write (6,"(4(a1,i3,a1,1x,f9.4,1x))") &
               "(",lat_out_idx_err,")",lat_out_grd(lat_out_idx_err),"(",lat_out_idx_err+1,")",lat_out_grd(lat_out_idx_err+1), &
               "(",lon_out_idx_err,")",lon_out_grd(lon_out_idx_err),"(",lon_out_idx_err+1,")",lon_out_grd(lon_out_idx_err+1)
          stop
       endif ! endif err
       
    end do ! end loop over time
    
    ! Close netCDF file
    rcd=rcd+nf90_wrp_close(nc_id,fl_in,"Ingested leaf area index data from") ! [fnc] Close file
    
    ! Archive input in human-readable format if desired
    if (.false.) then
       open (unit=unit_dgn,file=fl_dgn,status="unknown",iostat=rcd)
       if (rcd /= 0) write (6,"(4a,i4)") prg_nm(1:ftn_strlen(prg_nm)), & 
            ": ERROR "//sbr_nm//" unable to open ",fl_dgn(1:ftn_strlen(fl_dgn))
       do lat_in_idx=1,lat_in_nbr ! NB: Outer loop over lat
          do lon_in_idx=1,lon_in_nbr
             do pft_in_idx=0,pft_max_CLM ! NB: Loop starts at zero
                write (unit_dgn,"(2(f10.4))") &
                     lai_clm_in(lon_in_idx,lat_in_idx,pft_in_idx), & ! [m2]
                     vai_clm_in(lon_in_idx,lat_in_idx,pft_in_idx) ! [m2]
             end do ! end loop over pft
          end do ! end loop over lon
       end do ! end loop over lat
       close (unit_dgn)
       write (6,"(a,1x,a)") "Wrote leaf area index data to",fl_dgn(1:ftn_strlen(fl_dgn))
    endif ! endif dbg
    
    ! De-allocate
    if (allocated(area_in)) deallocate(area_in,stat=rcd)
    if(rcd /= 0) stop "deallocate() failed for area_in"
    if (allocated(lat_in_grd)) deallocate(lat_in_grd,stat=rcd)
    if(rcd /= 0) stop "deallocate() failed for lat_in_grd"
    if (allocated(lat_in_2d)) deallocate(lat_in_2d,stat=rcd)
    if(rcd /= 0) stop "deallocate() failed for lat_in_2d"
    if (allocated(lon_in_grd)) deallocate(lon_in_grd,stat=rcd)
    if(rcd /= 0) stop "deallocate() failed for lon_in_grd"
    if (allocated(lon_in_grd_2d)) deallocate(lon_in_grd_2d,stat=rcd)
    if(rcd /= 0) stop "deallocate() failed for lon_in_grd_2d"
    if (allocated(lon_in_2d)) deallocate(lon_in_2d,stat=rcd)
    if(rcd /= 0) stop "deallocate() failed for lon_in_2d"
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
    
    if (allocated(lai_clm_in)) deallocate(lai_clm_in,stat=rcd)
    if(rcd /= 0) stop "deallocate() failed for lai_clm_in"
    if (allocated(vai_clm_in)) deallocate(vai_clm_in,stat=rcd)
    if(rcd /= 0) stop "deallocate() failed for vai_clm_in"
    if (allocated(lai_clm_out)) deallocate(lai_clm_out,stat=rcd)
    if(rcd /= 0) stop "deallocate() failed for lai_clm_out"
    if (allocated(vai_clm_out)) deallocate(vai_clm_out,stat=rcd)
    if(rcd /= 0) stop "deallocate() failed for vai_clm_out"
    
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Exiting "//sbr_nm
    return 
  end subroutine lai_get ! end lai_get()
  
end module clm_fld ! [mdl] Community Land Model fields
