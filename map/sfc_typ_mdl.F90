! $Id$ -*-f90-*-

! Purpose: Surface type classification and derived properties of surface type

! Usage: 
! use sfc_typ_mdl ! [mdl] Surface type and derived properties

module sfc_typ_mdl ! [mdl] Surface type and derived properties
  implicit none
  private ! [stt] Symbols are private unless individually qualified as public
  public::sfc_typ_get ! [sbr] Create LSM surface type field
  public::sfc_typ_xtr_get ! [sbr] Retrieve and rebin surface type
  public::sfc2vgt ! [sbr] LAI, VAI prescribed for LSM surface type
  public::sfc2dps ! [sbr] Surface roughness prescribed for LSM surface type
  
contains
  
  subroutine sfc_typ_xtr_get( & ! [sbr] Retrieve and rebin surface type
       fl_in, & ! I
       lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
       area_out,lnd_frc_out,lnd_msk_out,flg_soi_nn, & ! I
       sfc_typ_LSM_out) ! O
    ! Purpose: Retrieve and rebin surface type
    ! sfc_typ_xtr_get() is called by lnd_bnd_cnv_rbn() in bds_drv.F90
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use map_cst ! [mdl] Constants used in map routines
    use map_grd ! [mdl] Map grids and regridding
    use netcdf ! [mdl] netCDF interface
    use nf90_utl ! [mdl] netCDF utilities
    use sng_mdl,only:ftn_strlen,ftn_strstr ! [mdl] String manipulation
    use utl_mdl,only:mnt_ncr_chk ! [mdl] Utility functions (date_time_get,mnt_chk...)
    implicit none
    ! Parameters
    character(len=*),parameter::sbr_nm="sfc_typ_xtr_get" ! [sng] Subroutine name
    !    character(len=*),parameter::fld_nm="SURF2D" ! [sng] Name of variable in input file
    character(len=*),parameter::fld_nm="sfc_typ" ! [sng] Name of variable in input file
    integer,parameter::sfc_typ_LSM_max=28 ! [enm] Maximum value of LSM surface type
    integer,parameter::sfc_typ_LSM_nbr=29 ! [nbr] Number of LSM surface types 
    integer,parameter::lat_nbr_max=180 ! [nbr] Maximum number of latitudes
    integer,parameter::lon_nbr_max=360 ! [nbr] Maximum number of longitudes
    integer,parameter::unit_dgn=74 ! Unit for writing general Map diagnoses
    ! Commons
    ! Input
    character(len=*),intent(in)::fl_in        ! [sng] Input file
    integer,intent(in)::lat_out_nbr       ! [nbr] Number of latitudes
    integer,intent(in)::lnd_msk_out(lon_out_nbr,lat_out_nbr) ! [msk] Land mask (integer 0 or 1)
    integer,intent(in)::lon_out_nbr       ! [nbr] Number of longitudes
    real,intent(in)::area_out(lon_out_nbr,lat_out_nbr) ! [m2] Area of gridcells
    real,intent(in)::lat_out_grd(lat_out_nbr+1) ! [dgr] Interface latitudes
    real,intent(in)::lnd_frc_out(lon_out_nbr,lat_out_nbr) ! [frc] Land fraction of gridcell 
    real,intent(in)::lon_out_grd(lon_out_nbr+1) ! [dgr] Interface longitudes
    logical,intent(in)::flg_soi_nn

    ! Output
    integer,intent(out)::sfc_typ_LSM_out(lon_out_nbr,lat_out_nbr) ! [enm] Surface type
    
    ! Locals with simple initialization
    character(80)::map_typ_sng=""//char(0) ! [sng] Map type description
    character(80)::fl_dgn="/tmp/zender/map/map_typ.txt"//char(0) ! [sng] Diagnostics file
    character(80)::fl_Ols="/data/mflanner/data/map_Ols.txt"//char(0) ! [sng] Diagnostics file
    integer::map_typ_in=map_lat_rgl_lon_180_ctr ! [enm] Latitudes are regular, Date line at center of first longitude bin
    integer::mss_val_ntg=-99999 ! [frc] Missing value
    integer::rcd=nf90_noerr ! [enm] Return success code
    logical::flg_read_2d_grd=.false. ! [flg] Read/use un-necessary CLM2 2D grid fields
    logical::flg_CLM=.false. ! [flg] Model input has CLM2 coding
    real::mss_val_in=1.0e36 ! [frc] Missing value
    real::mss_val_out=1.0e36 ! [frc] Missing value
    character(35) vgt_sng(0:sfc_typ_LSM_max) ! [sng] Description of surface blends
    logical::flg_set_crp_frs=.true. ! [flg] Set warm/cool crops, irrigated crops to forest crops
    integer::crp_frs_rpl_nbr=0 ! [nbr] Number of crop<-->forest replacements
    data vgt_sng( 0)/"ocean                              "/
    data vgt_sng( 1)/"land ice                           "/
    data vgt_sng( 2)/"desert                             "/
    data vgt_sng( 3)/"cool needleleaf evergreen forest   "/
    data vgt_sng( 4)/"cool needleleaf deciduous forest   "/
    data vgt_sng( 5)/"cool broadleaf deciduous forest    "/
    data vgt_sng( 6)/"cool mixed ne+bd forest            "/
    data vgt_sng( 7)/"warm needleleaf evergreen forest   "/
    data vgt_sng( 8)/"warm broadleaf deciduous forest    "/
    data vgt_sng( 9)/"warm mixed ne+bd forest            "/
    data vgt_sng(10)/"tropical broadleaf evergreen forest"/
    data vgt_sng(11)/"tropical seasonal deciduous forest "/
    data vgt_sng(12)/"savanna                            "/
    data vgt_sng(13)/"evergreen forest tundra            "/
    data vgt_sng(14)/"deciduous forest tundra            "/
    data vgt_sng(15)/"cool forest crop                   "/
    data vgt_sng(16)/"warm forest crop                   "/
    data vgt_sng(17)/"cool grassland                     "/
    data vgt_sng(18)/"warm grassland                     "/
    data vgt_sng(19)/"tundra                             "/
    data vgt_sng(20)/"evergreen shrubland                "/
    data vgt_sng(21)/"deciduous shrubland                "/
    data vgt_sng(22)/"semi-desert                        "/
    data vgt_sng(23)/"cool irrigated crop                "/
    data vgt_sng(24)/"cool crop                          "/
    data vgt_sng(25)/"warm irrigated crop                "/
    data vgt_sng(26)/"warm crop                          "/
    data vgt_sng(27)/"forest wetland                     "/
    data vgt_sng(28)/"non-forest wetland                 "/
    
    ! Local
    integer lat_in_idx        ! [idx] Counting index for lat
    integer lat_out_idx       ! [idx] Counting index for lat
    integer lon_in_idx        ! [idx] Counting index for lon
    integer lon_out_idx       ! [idx] Counting index for lon
    integer ovr_idx           ! [idx] Counting index
    integer ovr_nbr(lon_out_nbr,lat_out_nbr) ! [nbr] Number of input gridcells which overlap each output gridcell
    integer ovr_nbr_max       ! [nbr] Maximum number of input cells which overlap any output cell
    integer ovr_nbr_crr       ! [nbr] Current number of overlapping gridcells
    logical mnt_ncr           ! [flg] Monotonically increasing
    real edg_grd(4) ! [dgr] Grid edges (north, east, south, west)
    
    integer sfc_typ_LSM_crr ! [enm] Surface type of current gridcell
    integer sfc_typ_LSM_idx   ! [idx] Counting index for sfc_typ
    integer sfc_typ_rnk_1st   ! [idx] Index of most dominant LSM surface type
    integer sfc_typ_rnk_2nd   ! [idx] Index of second most dominant LSM surface type
    real area_sfc_typ_in(0:sfc_typ_LSM_max) ! [m2] Area of globe covered by each surface type
    real area_sfc_typ_out(0:sfc_typ_LSM_max) ! [m2] Area of globe covered by each surface type
    real sfc_typ_wgt(0:sfc_typ_LSM_max) ! [frc] Weight of surface type in gridcell
    integer sfc_typ_crr      ! Current surface type

    ! Local MF
    integer LSM2Ols(0:sfc_typ_LSM_max)      ! [int] LSM surface type for each Olson type
    integer temp_idx         ! [idx] Counting index for sfc_typ text output file
    integer lat_idx2         ! [idx] Counting index for nearest neighbor loop
    integer lon_idx2         ! [idx] Counting index for nearest neighbor loop
    integer lat_idx_min_dst  ! [idx] Lat index of nearest valid land gridcell for NN option
    integer lon_idx_min_dst  ! [idx] Lon index of nearest valid land gridcell for NN option
    integer nn_cnt           ! [idx] Count of how many gridcells are set set to nearest neighbor
    real distance_min        ! [km] Distance of closest valid gridcell for NN option
    real lat1_tmp            ! [nbr] Latitude1 used in NN distance calculation
    real lat2_tmp            ! [nbr] Latitude2 used in NN distance calculation
    real lon1_tmp            ! [nbr] Longitude1 used in NN distance calculation
    real lon2_tmp            ! [nbr] Longitude2 used in NN distance calculation
    real distance_tmp        ! [km] Distance between current 2 gridcells being evaluated in NN loop
    real,parameter::degtorad = 57.2958  ! [nbr] conversion factor for degrees to radians


    ! Variables needed to read external netCDF dataset
    integer cnt_lon_lat_time(3) ! [nbr] Dimension sizes
    integer edg_est_id ! [id] Variable ID
    integer edg_nrt_id ! [id] Variable ID
    integer edg_sth_id ! [id] Variable ID
    integer edg_wst_id ! [id] Variable ID
    integer lat_dmn_id ! [id] Dimension ID for lat
    integer lat_in_2d_id ! [id] Variable ID
    integer lat_in_nbr        ! [nbr] Number of latitudes
    integer lon_dmn_id ! [id] Dimension ID for lon
    integer lon_in_2d_id ! [id] Variable ID
    integer lon_in_nbr        ! [nbr] Number of longitudes
    integer mss_val_id        ! [id] Attribute ID
    integer nc_id             ! [id] File handle
    integer sfc_typ_id           ! [id] Variable ID
    integer sfc_typ_dmn_nbr ! [nbr] Number of dimensions in disk sfc_typ field
    integer srt_lon_lat_time(3) ! [idx] Dimension offset indices
    logical mss_flg           ! [flg] Check for mss_val on input data
    
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
    
    integer,dimension(:,:),allocatable::sfc_typ_LSM_in ! [enm] Surface type 
    integer,dimension(:,:),allocatable::sfc_typ_Ols_out ! [enm] Surface type output on Olson grid

    
    ! Main Code
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Entering "//sbr_nm
    
    ! Check for monotonicity
    mnt_ncr=mnt_ncr_chk(lon_out_grd,lon_out_nbr+1)
    mnt_ncr=mnt_ncr_chk(lat_out_grd,lat_out_nbr+1)
    
    ! Initialize scalars
    mss_flg=.false.
    ! Initialize arrays
    sfc_typ_LSM_out(:,:)=mss_val_ntg ! [enm] Surface type
    
    ! Read netCDF data
    rcd=nf90_noerr ! nf90_noerr == 0
    rcd=nf90_wrp_open(fl_in,NF90_NOWRITE,nc_id,sbr_nm=sbr_nm)
    ! Get dimension IDs
    rcd=nf90_wrp_inq_dimid(nc_id,"lat",lat_dmn_id)
    rcd=nf90_wrp_inq_dimid(nc_id,"lon",lon_dmn_id)
    ! Get dimension sizes
    rcd=rcd+nf90_inquire_dimension(nc_id,lon_dmn_id,len=lon_in_nbr)
    rcd=rcd+nf90_inquire_dimension(nc_id,lat_dmn_id,len=lat_in_nbr)
    if(rcd /= nf90_noerr) stop "Error retrieving dimension sizes"
    
    ! Enough memory? 
    if (lat_in_nbr > lat_nbr_max) stop "lat_in_nbr > lat_nbr_max in sfc_typ_xtr_get()"
    if (lon_in_nbr > lon_nbr_max) stop "lon_in_nbr > lon_nbr_max in sfc_typ_xtr_get()"
    
    ! Allocate based on input file dimensions
    allocate(area_in(lon_in_nbr,lat_in_nbr),stat=rcd) ! [m2] Area of gridcells
    if(rcd /= 0) stop "allocate() failed for area_in"
    allocate(lat_in_grd(lat_in_nbr+1),stat=rcd) ! [dgr] Interface latitudes
    if(rcd /= 0) stop "allocate() failed for lat_in_grd"
    allocate(lon_in_grd(lon_in_nbr+1),stat=rcd) ! [dgr] Interface longitudes
    if(rcd /= 0) stop "allocate() failed for lon_in_grd"
    allocate(lon_in_nbr_1d(lat_in_nbr),stat=rcd) ! [nbr] Number of longitudes per latitude
    if(rcd /= 0) stop "allocate() failed for lon_in_nbr_1d"
    if(flg_read_2d_grd) then
       allocate(lat_in_2d(lon_in_nbr,lat_in_nbr),stat=rcd) ! [dgr] Latitude at gridcell center
       if(rcd /= 0) stop "allocate() failed for lat_in_2d"
       allocate(lon_in_grd_2d(lon_in_nbr+1,lat_in_nbr),stat=rcd) ! [dgr] Interface longitudes
       if(rcd /= 0) stop "allocate() failed for lon_in_grd_2d"
       allocate(lon_in_2d(lon_in_nbr,lat_in_nbr),stat=rcd) ! [dgr] Longitude at gridcell center
       if(rcd /= 0) stop "allocate() failed for lon_in_2d"
    endif ! endif flg_read_2d_grd
    
    allocate(sfc_typ_LSM_in(lon_in_nbr,lat_in_nbr),stat=rcd) ! [enm] Surface type 
    if(rcd /= 0) stop "allocate() failed for sfc_typ_LSM_in"
    allocate(sfc_typ_Ols_out(lon_in_nbr*2,lat_in_nbr*2),stat=rcd) ! [enm] Surface type on Olson grid
    if(rcd /= 0) stop "allocate() failed for sfc_typ_Ols_out"


    ! Get variable IDs
    rcd=nf90_wrp_inq_varid(nc_id,fld_nm,sfc_typ_id)
    if(flg_CLM.and.flg_read_2d_grd) then
       rcd=nf90_wrp_inq_varid(nc_id,"LATIXY",lat_in_2d_id)
       rcd=nf90_wrp_inq_varid(nc_id,"LONGXY",lon_in_2d_id)
       rcd=nf90_wrp_inq_varid(nc_id,"EDGEN",edg_nrt_id)
       rcd=nf90_wrp_inq_varid(nc_id,"EDGEE",edg_est_id)
       rcd=nf90_wrp_inq_varid(nc_id,"EDGES",edg_sth_id)
       rcd=nf90_wrp_inq_varid(nc_id,"EDGEW",edg_wst_id)
    endif ! endif flg_CLM.and.flg_read_2d_grd
    ! Get number of dimensions
    rcd=rcd+nf90_inquire_variable(nc_id,sfc_typ_id,ndims=sfc_typ_dmn_nbr)
    if (sfc_typ_dmn_nbr /= 2) then
       write (6,"(2a,i1,a)") prg_nm(1:ftn_strlen(prg_nm)), & 
            ": ERROR "//fld_nm//" has ",sfc_typ_dmn_nbr," dimensions"
       write (6,"(2a,i1,a)") prg_nm(1:ftn_strlen(prg_nm)), & 
            ": HINT Make sure "//fld_nm//" has no time dimension"
       stop
    endif                  ! endif err
    ! Get time-independent data
    if(flg_CLM.and.flg_read_2d_grd) then
       rcd=nf90_wrp(nf90_get_var(nc_id,lat_in_2d_id,lat_in_2d),"get_var lat_in_2d")
       rcd=nf90_wrp(nf90_get_var(nc_id,lon_in_2d_id,lon_in_2d),"get_var lon_in_2d")
       rcd=nf90_wrp(nf90_get_var(nc_id,edg_nrt_id,edg_grd(1)),"get_var edg_nrt")
       rcd=nf90_wrp(nf90_get_var(nc_id,edg_est_id,edg_grd(2)),"get_var edg_est")
       rcd=nf90_wrp(nf90_get_var(nc_id,edg_sth_id,edg_grd(3)),"get_var edg_sth")
       rcd=nf90_wrp(nf90_get_var(nc_id,edg_wst_id,edg_grd(4)),"get_var edg_wst")
    endif ! endif flg_CLM.and.flg_read_2d_grd
    rcd=nf90_wrp(nf90_get_var(nc_id,sfc_typ_id,sfc_typ_LSM_in),"get_var sfc_typ_LSM_in")
    ! Get missing value
    rcd=rcd+nf90_inquire_attribute(nc_id,sfc_typ_id,"missing_value",attnum=mss_val_id)
    if (rcd==nf90_noerr) then
       mss_flg=.true. ! [flg] Data may have missing values
       rcd=rcd+nf90_get_att(nc_id,sfc_typ_id,"missing_value",mss_val_in)
       if (mss_val_in /= 1.0e36) write (6,"(2a,f12.6)") prg_nm(1:ftn_strlen(prg_nm)), &
            ": WARNING "//sbr_nm//" reports mss_val_in = ",mss_val_in
    else
       mss_flg=.false. ! [flg] Data may have missing values
       rcd=nf90_noerr ! [enm] Return success code
    endif ! endif
    if(rcd /= nf90_noerr) stop "Error retrieving variable data"
    
    ! Close file
    rcd=nf90_wrp_close(nc_id,fl_in,"Ingested",sbr_nm=sbr_nm) ! [fnc] Close file
    
    ! Determine map grid and overlap characteristics
    ! fxm: Method breaks for reduced grids, streamline, subroutinize, and generalize
    lon_in_nbr_1d(:)=lon_in_nbr ! [nbr] Number of longitudes per latitude
    
    ! Convert gridpoint centers, domain boundaries to grid interfaces
    if(flg_CLM.and.flg_read_2d_grd) then
       call map_edge_mk(lat_in_nbr,lon_in_nbr,lon_in_nbr_1d, & ! I
            lat_in_2d,lon_in_2d,edg_grd, & ! I 
            lat_in_grd,lon_in_grd,lon_in_grd_2d) ! O
    else
       call map_grd_mk(lat_in_nbr,lon_in_nbr,map_typ_in, & ! I
            lat_in_grd,lon_in_grd) ! O
       call map_typ_sng_get(map_typ_in, & ! I
            map_typ_sng) ! O
       write (6,"(2a,i2,2a)") prg_nm(1:ftn_strlen(prg_nm)), &  
            ": INFO "//sbr_nm//"() assumes input map type is ",map_typ_in, &
            " = ",map_typ_sng
    endif ! endif flg_CLM.and.flg_read_2d_grd
    
    mnt_ncr=mnt_ncr_chk(lat_in_grd,lat_in_nbr+1)
    if (.not.mnt_ncr) stop "ERROR: lat_in_grd not monotonically increasing in sfc_typ_xtr_get()"
    mnt_ncr=mnt_ncr_chk(lon_in_grd,lon_in_nbr+1)
    if (.not.mnt_ncr) stop "ERROR: lon_in_grd not monotonically increasing in sfc_typ_xtr_get()"
    
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
    
    ! Get overlap locations and weights for mapping input grid to output grid
    call map_ovr_wgt_drv( &
         lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
         lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
         ovr_nbr_max,area_out, & ! I
         ovr_lat_idx,ovr_lon_idx,ovr_nbr,ovr_wgt) ! O
    
    ! Rebin surface type from input grid to output grid
    ! Maximum fraction method adopted from sfc_typ_get()
    ! Process each cell on output grid

    nn_cnt=0
    do lat_out_idx=1,lat_out_nbr
       do lon_out_idx=1,lon_out_nbr
          ! Sum overlap weights by surface type
          do sfc_typ_LSM_idx=0,sfc_typ_LSM_max ! Loop starts with 0
             sfc_typ_wgt(sfc_typ_LSM_idx)=0.0 ! [frc] Weight of surface type in gridcell
          end do ! end loop over LSM surface type
          ovr_nbr_crr=ovr_nbr(lon_out_idx,lat_out_idx)
          if (ovr_nbr_crr==0) write (6,"(a,i3)") sbr_nm//"(): # overlapping input gridcells = ",ovr_nbr_crr
          do ovr_idx=1,ovr_nbr_crr ! Overlap cell index
             lon_in_idx=ovr_lon_idx(lon_out_idx,lat_out_idx,ovr_idx) ! [idx] lon index (input grid) of overlap cell
             lat_in_idx=ovr_lat_idx(lon_out_idx,lat_out_idx,ovr_idx) ! [idx] lat index (input grid) of overlap cell
             sfc_typ_LSM_crr=sfc_typ_LSM_in(lon_in_idx,lat_in_idx) ! [enm] Surface type of current gridcell
             sfc_typ_wgt(sfc_typ_LSM_crr)=sfc_typ_wgt(sfc_typ_LSM_crr)+ovr_wgt(lon_out_idx,lat_out_idx,ovr_idx) ! [frc] Weight of surface type in gridcell
          end do ! end loop over overlapping cells
          
          ! Rank non-zero weights by surface type 
          ! sfc_typ_wgt_1st is the most extensive surface type, sfc_typ_wgt_2nd is the second most extensive surface type
          ! NB: LSM code is rather ambiguous whether these ranks are indices or values
          call rnk_vec_LSM(sfc_typ_wgt,sfc_typ_LSM_max,mss_val_ntg, & ! I
               sfc_typ_rnk_1st,sfc_typ_rnk_2nd) ! O
          if (sfc_typ_rnk_1st==mss_val_ntg) then
             write (6,"(a,i3,a,f9.5)") sbr_nm//"(): ERROR occured at lon_out_grd(",lon_out_idx,") = ",lon_out_grd(lon_out_idx)
             write (6,"(a,i3,a,f9.5)") sbr_nm//"(): ERROR occured at lat_out_grd(",lat_out_idx,") = ",lat_out_grd(lat_out_idx)
             write (6,"(a,i3)") sbr_nm//"(): # overlapping input gridcells = ",ovr_nbr_crr
             write (6,"(a)") "lon_in_idx lat_in_idx sfc_typ_LSM_in ovr_wgt" 
             do ovr_idx=1,ovr_nbr_crr !overlap cell index
                lon_in_idx=ovr_lon_idx(lon_out_idx,lat_out_idx,ovr_idx) ! [idx] lon index (input grid) of overlap cell
                lat_in_idx=ovr_lat_idx(lon_out_idx,lat_out_idx,ovr_idx) ! [idx] lat index (input grid) of overlap cell
                sfc_typ_LSM_crr=sfc_typ_LSM_in(lon_in_idx,lat_in_idx) ! [enm] Surface type of current gridcell
                sfc_typ_wgt(sfc_typ_LSM_crr)=sfc_typ_wgt(sfc_typ_LSM_crr)+ovr_wgt(lon_out_idx,lat_out_idx,ovr_idx) ! [frc] Weight of surface type in gridcell
                write (6,"(i3,1x,i3,1x,i2,1x,f9.6)") lon_in_idx,lat_in_idx, &
                     sfc_typ_LSM_crr,ovr_wgt(lon_out_idx,lat_out_idx,ovr_idx)
             end do ! end loop over overlapping cells
             stop
          endif               ! endif err
          
          ! Set surface type as:
          ! If fractional land = 0: cell = ocean
          ! If fractional land > 0: cell = land
          ! a. Use most frequent surface type based on area of overlap unless this is ocean 
          ! In this case, input grid says ocean but output grid wants land
          ! b. So use next most extensive surface type so long as is not ocean
          ! c. If this is ocean, or if there is none,
          !    If flg_soi_nn is FALSE (original method) use swamp/marsh surface type
          !    This algorithm, originally used by LSM, fails for changes in mean sea 
          !    level because new land may appear that should not be swamp/marsh
          !    A nearest neighbor algorithm would be much more appropriate
          !    Same is true of lai_get(), pft_get(), soi_txt_get()
          !    If flg_soi_nn is TRUE (MF added) then set undefined land surface type 
          !    to land type of nearest (originally-defined) land gridcell

          if (lnd_frc_out(lon_out_idx,lat_out_idx)==0.0) then ! Ocean
             sfc_typ_LSM_out(lon_out_idx,lat_out_idx)=0
          else ! not ocean
             if (sfc_typ_rnk_1st /= 0) then
                sfc_typ_LSM_out(lon_out_idx,lat_out_idx)=sfc_typ_rnk_1st 
             else ! Dominant type is ocean, consider second-ranked type...
                if (sfc_typ_rnk_2nd==0.or.sfc_typ_rnk_2nd==mss_val_ntg) then
                   
                   ! Nearest Neighbor block
                   if (flg_soi_nn) then
                      
                      distance_min = 100000 ! [km] Initialized to high value
                      lat1_tmp = lat_out_grd(lat_out_idx) / degtorad ! [nbr] radians
                      lon1_tmp = lon_out_grd(lon_out_idx) / degtorad ! [nbr] radians
                      
                      ! Find nearest neighbor on input grid
                      do lat_idx2=1,lat_in_nbr
                         do lon_idx2=1,lon_in_nbr
                            if ( (sfc_typ_LSM_in(lon_idx2,lat_idx2).ne.0) ) then
                             !    (sfc_typ_LSM_in(lon_idx2,lat_idx2).ne.1) )then
                         
                               lat2_tmp = lat_in_grd(lat_idx2) / degtorad ! [nbr] radians
                               lon2_tmp = lon_in_grd(lon_idx2) / degtorad ! [nbr] radians
                            
                               ! Great Circle Distance Formula (kilometers)
                               distance_tmp=6378.7*acos((sin(lat1_tmp)*sin(lat2_tmp))+ &
                                    (cos(lat1_tmp)*cos(lat2_tmp)*cos(abs(lon2_tmp-lon1_tmp))))
                               
                               if (distance_tmp.lt.distance_min) then
                                  distance_min = distance_tmp
                                  lat_idx_min_dst = lat_idx2
                                  lon_idx_min_dst = lon_idx2
                               endif   ! endif
                            endif   ! endif
                         enddo   ! end loop over lon
                      enddo    ! end loop over lat
                      sfc_typ_LSM_out(lon_out_idx,lat_out_idx)=sfc_typ_LSM_in(lon_idx_min_dst,lat_idx_min_dst)
                      nn_cnt = nn_cnt + 1
                   else 
                      sfc_typ_LSM_out(lon_out_idx,lat_out_idx)=28   ! Original method: set gridcell to wetland/marsh
                   endif ! endif nearest neighbor

                else
                   sfc_typ_LSM_out(lon_out_idx,lat_out_idx)=sfc_typ_rnk_2nd
                end if ! endif
             end if ! endif dominant type is ocean
          end if ! endif lnd_frc_out

          ! Perform any arbitrary post-processing desired for dust data
          ! Operations in this block are not performed by LSM
          ! Thus each remapping causes discrepancy between LSM and dust datasets
          if (flg_set_crp_frs) then
             ! Set warm/cool crops to forest crops
             ! This considerably reduces bare ground fraction in northern Europe, midwest US in winter
             sfc_typ_LSM_crr=sfc_typ_LSM_out(lon_out_idx,lat_out_idx) ! [enm] Surface type of current gridcell
             if (sfc_typ_LSM_crr >= 23 .and. sfc_typ_LSM_crr <=26) crp_frs_rpl_nbr=crp_frs_rpl_nbr+1 ! [nbr] Number of crop<-->forest replacements
             ! Change cool irrigated crop to cool forest crop
             if (sfc_typ_LSM_crr==23) sfc_typ_LSM_out(lon_out_idx,lat_out_idx)=15
             ! Change cool crop to cool forest crop
             if (sfc_typ_LSM_crr==24) sfc_typ_LSM_out(lon_out_idx,lat_out_idx)=15
             ! Change warm irrigated crop to warm forest crop
             if (sfc_typ_LSM_crr==25) sfc_typ_LSM_out(lon_out_idx,lat_out_idx)=16
             ! Change warm crop to warm forest crop
             if (sfc_typ_LSM_crr==26) sfc_typ_LSM_out(lon_out_idx,lat_out_idx)=16
          endif ! endif flg_set_crp_frs
          
          ! Sanity checks
          sfc_typ_LSM_crr=sfc_typ_LSM_out(lon_out_idx,lat_out_idx) ! [enm] Surface type of current gridcell
          if (sfc_typ_LSM_crr < 0.or.sfc_typ_LSM_crr > sfc_typ_LSM_max) then
             write (6,"(a,i10,a,2(i3,a1))") &
                  "ERROR sfc_typ_xtr_get() reports invalid LSM surface type = ", &
                  sfc_typ_LSM_crr, &
                  " at (lon_out_idx,lat_out_idx) = (",lon_out_idx,",",lat_out_idx,")"
             stop
          end if ! endif
          if ( &
               (lnd_frc_out(lon_out_idx,lat_out_idx)==0.0.and.sfc_typ_LSM_crr /= 0).or. & ! Land fraction says ocean but sfc_typ_LSM says land
               (lnd_frc_out(lon_out_idx,lat_out_idx) > 0.0.and.sfc_typ_LSM_crr==0) & ! Land fraction says land but sfc_typ_LSM says ocean
               ) then
             write (6,"(a,i2,a,f9.6,a,2(i3,a1))") "ERROR sfc_type_get() reports sfc_typ_LSM_out = ", &
                  sfc_typ_LSM_crr, &
                  ", lnd_frc_out = ",lnd_frc_out(lon_out_idx,lat_out_idx), &
                  " at (lon_out_idx,lat_out_idx) = (",lon_out_idx,",",lat_out_idx,")"
             stop
          end if ! endif
          
       end do ! end loop over lon
    end do ! end loop over lat
    
    if (nn_cnt > 0) then
       write (6,"(a,i4,a)") "DEBUG: sfc_typ_xtr_get() reports ",nn_cnt," points set to nearest neighbor"
    endif         ! endif

    if (flg_set_crp_frs) then ! [flg] Set warm/cool crops, irrigated crops to forest crops
       write (6,"(2a,i6,a)") prg_nm(1:ftn_strlen(prg_nm)), & 
            ": INFO sfc_typ_xtr_get() replaced crops by forest crops in ", &
            crp_frs_rpl_nbr," gridcells"
    endif ! endif flg_set_crp_frs

    ! Compute global area of each LSM surface type 
    do sfc_typ_LSM_idx=0,sfc_typ_LSM_max ! Loop starts with 0
       area_sfc_typ_in(sfc_typ_LSM_idx)=0.0
    end do ! end loop over sfc_typ_LSM
    do lat_in_idx=1,lat_in_nbr
       do lon_in_idx=1,lon_in_nbr
          sfc_typ_LSM_crr=sfc_typ_LSM_in(lon_in_idx,lat_in_idx) ! [enm] Surface type of current gridcell
          area_sfc_typ_in(sfc_typ_LSM_crr)=area_sfc_typ_in(sfc_typ_LSM_crr)+area_in(lon_in_idx,lat_in_idx)
       end do ! end loop over lon
    end do ! end loop over lat
    
    ! Compute global area of each LSM surface type on output grid
    do sfc_typ_LSM_idx=0,sfc_typ_LSM_max ! Loop starts with 0
       area_sfc_typ_out(sfc_typ_LSM_idx)=0.0
    end do ! end loop over sfc_typ_LSM
    do lat_out_idx=1,lat_out_nbr
       do lon_out_idx=1,lon_out_nbr
          sfc_typ_LSM_crr=sfc_typ_LSM_out(lon_out_idx,lat_out_idx) ! [enm] Surface type of current gridcell
          area_sfc_typ_out(sfc_typ_LSM_crr)=area_sfc_typ_out(sfc_typ_LSM_crr)+area_out(lon_out_idx,lat_out_idx)
       end do ! end loop over lon
    end do ! end loop over lat
    
    ! Compare areas
    if (dbg_lvl==dbg_crr) then
       open (unit_dgn,file=fl_dgn,status="unknown",iostat=rcd)
       if (rcd /= 0) write (6,"(4a,i4)") prg_nm(1:ftn_strlen(prg_nm)), & 
            ": ERROR sfc_typ_xtr_get() unable to open ",fl_dgn(1:ftn_strlen(fl_dgn))
       write (unit_dgn,"(a)") sbr_nm//"() Surface Type Diagnostics File"
       write (unit_dgn,"(a,f16.3,a)") "Output grid cell area nearest the Equator: ",area_out(1,lat_out_nbr/2)*1.0e-6, " km2"
       write (unit_dgn,"(a,20x,a,1x,a)") "Vegetation type"," Input grid area","Output grid area"
       write (unit_dgn,"(33x,7x,a,8x,a)") "km2","km2"
       do sfc_typ_LSM_idx=0,sfc_typ_LSM_max ! Loop starts with 0
          write (unit_dgn,"(a35,f16.3,f17.3)") vgt_sng(sfc_typ_LSM_idx), &
               area_sfc_typ_in(sfc_typ_LSM_idx)*1.0e-6,area_sfc_typ_out(sfc_typ_LSM_idx)*1.0e-6
       end do ! end loop over sfc_typ_LSM
       close (unit_dgn)
       write (6,"(a,a)") "Wrote diagnostic land surface data to ",fl_dgn(1:ftn_strlen(fl_dgn))
    endif                     ! endif dbg

   
    ! Write sfc_typ data to text file using Olson format
    if (.false.) then
       
       LSM2Ols(0)=65
       LSM2Ols(1)=70
       LSM2Ols(2)=50
       LSM2Ols(3)=20
       LSM2Ols(4)=61
       LSM2Ols(5)=61
       LSM2Ols(6)=23
       LSM2Ols(7)=27
       LSM2Ols(8)=24
       LSM2Ols(9)=25
       LSM2Ols(10)=28
       LSM2Ols(11)=28
       LSM2Ols(12)=32
       LSM2Ols(13)=62
       LSM2Ols(14)=63
       LSM2Ols(15)=55
       LSM2Ols(16)=56
       LSM2Ols(17)=40
       LSM2Ols(18)=41
       LSM2Ols(19)=53
       LSM2Ols(20)=46
       LSM2Ols(21)=59
       LSM2Ols(22)=49
       LSM2Ols(23)=38
       LSM2Ols(24)=30
       LSM2Ols(25)=37
       LSM2Ols(26)=31
       LSM2Ols(27)=72
       LSM2Ols(28)=36


       ! Transform LSM (out) surface types to Olson surface types and write to output file in 0.5 x 0.5 resolution
       ! Assumes output resolution of 2x2 and map_typ_out=11
       open (unit_dgn,file=fl_Ols,status="unknown",iostat=rcd)
       if (rcd /= 0) write (6,"(4a,i4)") prg_nm(1:ftn_strlen(prg_nm)), & 
            ": ERROR sfc_typ_get() unable to open ",fl_Ols(1:ftn_strlen(fl_Ols))
       
       do lat_out_idx=lat_out_nbr,1,-1
          do temp_idx=1,4
             do lon_out_idx=1,lon_out_nbr   
                sfc_typ_crr = LSM2Ols(sfc_typ_LSM_out(lon_out_idx,lat_out_idx))
                write (unit_dgn,"(i2)") sfc_typ_crr
                write (unit_dgn,"(i2)") sfc_typ_crr
                write (unit_dgn,"(i2)") sfc_typ_crr
                write (unit_dgn,"(i2)") sfc_typ_crr
             end do ! end loop over lon
          end do
       end do ! end loop over lat
       close (unit_dgn)
       write (6,"(a,1x,a)") "Wrote LSM data in Olson format to ",fl_Ols(1:ftn_strlen(fl_Ols))
       
    endif ! endif write output file



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
    
    if (allocated(sfc_typ_LSM_in)) deallocate(sfc_typ_LSM_in,stat=rcd)
    if(rcd /= 0) stop "deallocate() failed for sfc_typ_LSM_in"
    
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Exiting "//sbr_nm
    return 
  end subroutine sfc_typ_xtr_get

  subroutine sfc_typ_get( & ! [sbr] Create LSM surface type field
       fl_in, & ! I
       lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
       lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
       ovr_nbr_max,area_out,lnd_frc_out,lnd_msk_out,flg_soi_nn, & ! I
       sfc_typ_LSM_out) ! O
    ! Purpose: Retrieve land fraction and mask and rebin to requested grid
    ! sfc_typ_xtr_get() is called by lnd_bnd_cnv_rbn() in bds_drv.F90
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use map_grd ! [mdl] Map grids and regridding
    use sng_mdl,only:ftn_strlen ! [mdl] String manipulation
    use utl_mdl,only:mnt_ncr_chk ! [mdl] Utility functions (date_time_get,mnt_chk...)
    implicit none
    ! Parameters
    character(len=*),parameter::sbr_nm="sfc_typ_get" ! [sng] Subroutine name
    integer,parameter::fl_in_unit=73 ! Unit for reading input data
    integer,parameter::sfc_typ_LSM_max=28 ! [enm] Maximum value of LSM surface type
    integer,parameter::sfc_typ_LSM_nbr=29 ! [nbr] Number of LSM surface types 
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
    real,intent(in)::lnd_frc_out(lon_out_nbr,lat_out_nbr) ! [frc] Land fraction of gridcell 
    integer,intent(in)::lnd_msk_out(lon_out_nbr,lat_out_nbr) ! [msk] Land mask (integer 0 or 1)
    real,intent(in)::lat_in_grd(lat_in_nbr+1) ! [dgr] Interface latitudes
    real,intent(in)::lat_out_grd(lat_out_nbr+1) ! [dgr] Interface latitudes
    real,intent(in)::lon_in_grd(lon_in_nbr+1) ! [dgr] Interface longitudes
    real,intent(in)::lon_out_grd(lon_out_nbr+1) ! [dgr] Interface longitudes
    logical,intent(in)::flg_soi_nn
    ! Output
    integer,intent(out)::sfc_typ_LSM_out(lon_out_nbr,lat_out_nbr) ! [enm] Surface type 
    ! Locals with simple initialization
    character(80)::fl_dgn="/tmp/zender/map/map_typ.txt"//char(0) ! [sng] Diagnostics file
    character(80)::fl_LSM="/tmp/zender/map/map_Ols.txt"//char(0) ! [sng] Diagnostics file
    character(80)::fl_Ols="/data/mflanner/data/map_Ols.txt"//char(0) ! [sng] Diagnostics file
    character(35) vgt_sng(0:sfc_typ_LSM_max) ! [sng] Description of surface blends
    data vgt_sng( 0)/"ocean                              "/
    data vgt_sng( 1)/"land ice                           "/
    data vgt_sng( 2)/"desert                             "/
    data vgt_sng( 3)/"cool needleleaf evergreen forest   "/
    data vgt_sng( 4)/"cool needleleaf deciduous forest   "/
    data vgt_sng( 5)/"cool broadleaf deciduous forest    "/
    data vgt_sng( 6)/"cool mixed ne+bd forest            "/
    data vgt_sng( 7)/"warm needleleaf evergreen forest   "/
    data vgt_sng( 8)/"warm broadleaf deciduous forest    "/
    data vgt_sng( 9)/"warm mixed ne+bd forest            "/
    data vgt_sng(10)/"tropical broadleaf evergreen forest"/
    data vgt_sng(11)/"tropical seasonal deciduous forest "/
    data vgt_sng(12)/"savanna                            "/
    data vgt_sng(13)/"evergreen forest tundra            "/
    data vgt_sng(14)/"deciduous forest tundra            "/
    data vgt_sng(15)/"cool forest crop                   "/
    data vgt_sng(16)/"warm forest crop                   "/
    data vgt_sng(17)/"cool grassland                     "/
    data vgt_sng(18)/"warm grassland                     "/
    data vgt_sng(19)/"tundra                             "/
    data vgt_sng(20)/"evergreen shrubland                "/
    data vgt_sng(21)/"deciduous shrubland                "/
    data vgt_sng(22)/"semi-desert                        "/
    data vgt_sng(23)/"cool irrigated crop                "/
    data vgt_sng(24)/"cool crop                          "/
    data vgt_sng(25)/"warm irrigated crop                "/
    data vgt_sng(26)/"warm crop                          "/
    data vgt_sng(27)/"forest wetland                     "/
    data vgt_sng(28)/"non-forest wetland                 "/
    logical::flg_set_crp_frs=.true. ! [flg] Set warm/cool crops, irrigated crops to forest crops
    integer::rcd=0 ! [enm] Return success code
    integer::crp_frs_rpl_nbr=0 ! [nbr] Number of crop<-->forest replacements
    integer::mss_val_ntg=-99999 ! [frc] Missing value
    
    ! Local
    integer Ols2LSM(100) ! LSM surface type for each Olson type
    integer bfr_in(15,2) ! Buffer for storing up to 15 (nbr,val) pairs
    integer err_nbr           ! [nbr] Number of errors processing surface types
    integer idx               ! [idx] Counting index
    integer lat_idx_nrt       ! [idx] Latitude index of cell to the North
    integer lat_idx_sth       ! [idx] Latitude index of cell to the South
    integer lat_in_idx        ! [idx] Counting index for lat
    integer lat_in_idx_err    ! [idx] Latitude of gridcell which caused error
    integer lat_out_idx       ! [idx] Counting index for lat
    integer lon_idx_est       ! [idx] Longitude index of cell to the East
    integer lon_idx_wst       ! [idx] Longitude index of cell to the West
    integer lon_in_idx        ! [idx] Counting index for lon
    integer lon_in_idx_err    ! [idx] Longitude of gridcell which caused error
    integer lon_out_idx       ! [idx] Counting index for lon
    integer ols_idx           ! [idx] Counting index
    integer ovr_idx           ! [idx] Counting index
    integer ovr_lat_idx(lon_out_nbr,lat_out_nbr,ovr_nbr_max) ! [idx] Map into input grid of latitude indices of overlap cells
    integer ovr_lon_idx(lon_out_nbr,lat_out_nbr,ovr_nbr_max) ! [idx] Map into input grid of longitude indices of overlap cells
    integer ovr_nbr(lon_out_nbr,lat_out_nbr) ! [nbr] Number of input gridcells which overlap each output gridcell
    integer ovr_nbr_crr       ! [nbr] Current number of overlapping gridcells
    integer pr_idx            ! [idx] Counting index for pr
    integer pr_nbr            ! [nbr] Number (nbr,val) pairs in this latitude band
    integer sfc_typ_LSM_crr ! [enm] Surface type of current gridcell
    integer sfc_typ_LSM_idx   ! [idx] Counting index for sfc_typ
    integer sfc_typ_LSM_in(lon_in_nbr,lat_in_nbr) ! LSM surface type on Olson (input) grid
    integer sfc_typ_Ols_crr   ! [enm] Surface type of current Olson gridcell
    integer sfc_typ_Ols_in(lon_in_nbr,lat_in_nbr) ! Olson surface type on Olson (input) grid
    integer sfc_typ_Ols_prc(lon_in_nbr,lat_in_nbr) ! Olson surface type on Olson (input) grid processed with intermediary LSM algorithms
    integer sfc_typ_rnk_1st   ! [idx] Index of most dominant LSM surface type
    integer sfc_typ_rnk_2nd   ! [idx] Index of second most dominant LSM surface type
    logical mnt_ncr           ! [flg] Monotonic and increasing flag
    real area_in(lon_in_nbr,lat_in_nbr) ! [m2] Area of gridcells
    real area_sfc_typ_in(0:sfc_typ_LSM_max) ! [m2] Area of globe covered by each surface type
    real area_sfc_typ_out(0:sfc_typ_LSM_max) ! [m2] Area of globe covered by each surface type
    real ovr_wgt(lon_out_nbr,lat_out_nbr,ovr_nbr_max) ! [frc] Weight of overlapping input gridcells onto each output gridcell
    real sfc_typ_wgt(0:sfc_typ_LSM_max) ! [frc] Weight of surface type in gridcell

    ! Local MF
    integer sfc_typ_crr      ! [int] Current surface type
    integer LSM2Ols(0:sfc_typ_LSM_max)      ! [int] LSM surface type for each Olson type
    integer temp_idx         ! [idx] Counting index for sfc_typ text output file
    integer lat_idx2         ! [idx] Counting index for nearest neighbor loop
    integer lon_idx2         ! [idx] Counting index for nearest neighbor loop
    integer lat_idx_min_dst  ! [idx] Lat index of nearest valid land gridcell for NN option
    integer lon_idx_min_dst  ! [idx] Lon index of nearest valid land gridcell for NN option
    integer nn_cnt           ! [idx] Count of how many gridcells are set set to nearest neighbor
    real distance_min        ! [km] Distance of closest valid gridcell for NN option
    real lat1_tmp            ! [nbr] Latitude1 used in NN distance calculation
    real lat2_tmp            ! [nbr] Latitude2 used in NN distance calculation
    real lon1_tmp            ! [nbr] Longitude1 used in NN distance calculation
    real lon2_tmp            ! [nbr] Longitude2 used in NN distance calculation
    real distance_tmp        ! [km] Distance between current 2 gridcells being evaluated in NN loop
    real,parameter::degtorad = 57.2958  ! [nbr] conversion factor for degrees to radians


    
    ! Main Code
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Entering sfc_typ_get()"
    
    ! Check for monotonically increasing grids
    mnt_ncr=mnt_ncr_chk(lon_out_grd,lon_out_nbr+1)
    if (.not.mnt_ncr) stop "lon_out_grd not monotonically increasing in sfc_typ_get()"
    mnt_ncr=mnt_ncr_chk(lat_out_grd,lat_out_nbr+1)
    if (.not.mnt_ncr) stop "lat_out_grd not monotonically increasing in sfc_typ_get()"
    
    ! Initialize output values
    sfc_typ_LSM_out(:,:)=0 ! [enm] Surface type
    
    ! Compute gridcell area 
    call map_area_get(lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
         area_in) ! O
    
    ! Olson's data, aka NCAR SCD DSS DS769.0, is stored in a packed, formatted text file
    ! Each latitude band is written separately and preceded by the number of "pairs" it contains
    ! Each pair is of the form (nbr,val) where nbr is the number of consecutive gridpoints of surface type val
    ! Up to fifteen pairs are written on each line
    ! The first line for a given latitude has a slightly different format from the succeeding lines, if any
    ! Read strategy is read up to a newline at a time, unpacking as we go
    ! This means the temporary storage buffer, bfr_in, need be only a 15x2 matrix
    ! NB: Olson grid is stored with lon starting at date line and lat starting at North Pole
    ! open (fl_in_unit,file=fl_in,status="old")
    open (fl_in_unit,file=fl_in,status="old",iostat=rcd)
    if (rcd /= 0) write (6,"(4a,i4)") prg_nm(1:ftn_strlen(prg_nm)), & 
         ": ERROR sfc_typ_get() unable to open ",fl_in(1:ftn_strlen(fl_in))
    read (fl_in_unit,*) sfc_typ_Ols_in ! Read it all in at once
    ! Old code to read DS769.0
    if (.false.) then
       do lat_in_idx=1,lat_in_nbr ! NB: Outer loop over lat
          lon_in_idx=0
          ! Ingest first line for this latitude
          read (fl_in_unit,"(i3,15(i3,i2))") pr_nbr,(bfr_in(pr_idx,1),bfr_in(pr_idx,2),pr_idx=1,15)
          pr_idx=1
          do while (pr_idx <= 15.and.lon_in_idx < 720)
             do idx=1,bfr_in(pr_idx,1)
                lon_in_idx=lon_in_idx+1
                sfc_typ_Ols_in(lon_in_idx,lat_in_idx)=bfr_in(pr_idx,2)
             end do ! end loop over cns
             pr_idx=pr_idx+1
          end do ! end while
          ! Ingest remaining lines, if any, for this latitude
          do while (lon_in_idx /= 720) 
             read (fl_in_unit,"(3x,15(i3,i2),2x)") (bfr_in(pr_idx,1),bfr_in(pr_idx,2),pr_idx=1,15)
             pr_idx=1
             do while (pr_idx <= 15.and.lon_in_idx < 720) 
                do idx=1,bfr_in(pr_idx,1)
                   lon_in_idx=lon_in_idx+1
                   sfc_typ_Ols_in(lon_in_idx,lat_in_idx)=bfr_in(pr_idx,2)
                end do ! end loop over cns
                pr_idx=pr_idx+1
             end do ! end while pairs remain in current line
          end do ! end while pairs remain in current lat
          if (lon_in_idx /= 720) stop "lon_in_idx /= 720 in sfc_typ_get()"
       end do ! end loop over lat
    endif ! endif false
    close (fl_in_unit)
    write (6,"(a,1x,a)") "Read surface type data from",fl_in(1:ftn_strlen(fl_in))
    
    ! Write surface type data in CSM/CCM/LSM format (one value per line)
    ! These data should be very similar to /fs/cgd/csm/input/lnd/sfc_typ_olson.data
    ! diff -c -w /tmp/zender/map/map_Ols.txt /fs/cgd/csm/input/lnd/sfc_typ_olson.data | m
    if (dbg_lvl==dbg_crr) then
       open (unit_dgn,file=fl_LSM,status="unknown",iostat=rcd)
       if (rcd /= 0) write (6,"(4a,i4)") prg_nm(1:ftn_strlen(prg_nm)), & 
            ": ERROR sfc_typ_get() unable to open ",fl_LSM(1:ftn_strlen(fl_LSM))
       do lat_in_idx=1,lat_in_nbr ! NB: Outer loop over lat
          do lon_in_idx=1,lon_in_nbr
             write (unit_dgn,"(i2)") sfc_typ_Ols_in(lon_in_idx,lat_in_idx)
          end do ! end loop over lon
       end do ! end loop over lat
       close (unit_dgn)
       write (6,"(a,1x,a)") "Wrote Olson data in LSM format to",fl_LSM(1:ftn_strlen(fl_LSM))
    endif ! endif dbg
    
    ! Process Olson data in accordance with LSM algorithms
    ! NB: Following code from LSM
    do lat_in_idx=1,lat_in_nbr ! NB: Outer loop over lat
       do lon_in_idx=1,lon_in_nbr
          sfc_typ_Ols_crr=sfc_typ_Ols_in(lon_in_idx,lat_in_idx)
          
          ! There are several values Olson=2, Olson=6, Olson=8 in the data set that are not defined 
          ! Use neighboring cells in this order: West, East, North, South
          lon_idx_wst=max(1,lon_in_idx-1)
          lon_idx_est=min(lon_in_nbr,lon_in_idx+1)
          lat_idx_nrt=min(1,lat_in_idx-1)
          lat_idx_sth=max(lat_in_nbr,lat_in_idx+1)
          
          if (sfc_typ_Ols_crr==2.or.sfc_typ_Ols_crr==6.or.sfc_typ_Ols_crr==8) &
               sfc_typ_Ols_crr=sfc_typ_Ols_in(lon_idx_wst,lat_in_idx)
          if (sfc_typ_Ols_crr==2.or.sfc_typ_Ols_crr==6.or.sfc_typ_Ols_crr==8) &
               sfc_typ_Ols_crr=sfc_typ_Ols_in(lon_idx_est,lat_in_idx)
          if (sfc_typ_Ols_crr==2.or.sfc_typ_Ols_crr==6.or.sfc_typ_Ols_crr==8) &
               sfc_typ_Ols_crr=sfc_typ_Ols_in(lon_in_idx,lat_idx_nrt)
          if (sfc_typ_Ols_crr==2.or.sfc_typ_Ols_crr==6.or.sfc_typ_Ols_crr==8) &
               sfc_typ_Ols_crr=sfc_typ_Ols_in(lon_in_idx,lat_idx_sth)
          
          ! Split Antarctica (17) into polar desert (69) and ice (70)
          if (sfc_typ_Ols_crr==17) then
             if (lat_in_idx <= 313) then
                sfc_typ_Ols_crr=69 ! Polar desert
             else
                sfc_typ_Ols_crr=70 ! Ice
             end if
          end if              ! endif Antarctica
          
          ! Olson=61 (eastern south taiga) will be classified as needleleaf deciduous tree
          ! Change Olson=61 to Olson=20 (main taiga = needleleaf evergreen tree) based on longitude
          if (sfc_typ_Ols_crr==61.and.lon_in_idx <= 576) sfc_typ_Ols_crr=20 
          
          ! Olson=61 (eastern south taiga) will be classified needleleaf deciduous tree 
          ! Create additional needleleaf deciduous tree from Olson=21 (main taiga) and Olson=60 (southern taiga) based on longitude
          if (sfc_typ_Ols_crr==21.and.lon_in_idx >= 555) sfc_typ_Ols_crr=61  
          if (sfc_typ_Ols_crr==60.and.lon_in_idx >= 582) sfc_typ_Ols_crr=61   
          
          ! Change Olson=26 (warm mixed) to broad-leaved humid forest based on latitude
          if (sfc_typ_Ols_crr==26.and.lat_in_idx >= 113) sfc_typ_Ols_crr=29
          
          ! Split forest tundra (62, 63) into needleleaf evergreen forest tundra (62) and needleleaf deciduous forest tundra (63) based on longitude
          if (sfc_typ_Ols_crr==63) sfc_typ_Ols_crr=62
          if (sfc_typ_Ols_crr==62.and.lon_in_idx >= 490) sfc_typ_Ols_crr=63   
          
          ! Remap Olson grid so latitude monotonically increases from South Pole
          lat_out_idx=lat_in_nbr-lat_in_idx+1
          sfc_typ_Ols_prc(lon_in_idx,lat_out_idx)=sfc_typ_Ols_crr
       end do ! end loop over lon
    end do ! end loop over lat
    
    ! Sanity check
    err_nbr=0
    do lat_in_idx=1,lat_in_nbr
       do lon_in_idx=1,lon_in_nbr
          if (sfc_typ_Ols_prc(lon_in_idx,lat_in_idx) > 100.or.sfc_typ_Ols_prc(lon_in_idx,lat_in_idx) < 0) then
             err_nbr=err_nbr+1
             lon_in_idx_err=lon_in_idx
             lat_in_idx_err=lat_in_idx
          end if              ! endif err
       end do ! end loop over lon
    end do ! end loop over lat
    if (err_nbr > 0) then
       write (6,*) err_nbr," land type errors. Last one is Olson surface type = ", &
            sfc_typ_Ols_prc(lon_in_idx_err,lat_in_idx_err)," is undef for lon,lat = ",lon_in_idx_err,lat_in_idx_err
       stop
    end if                    ! endif err
    
    ! Initialize all LSM surface types on input grid to missing value 
    sfc_typ_LSM_in(:,:)=mss_val_ntg ! [enm] Surface type
    
    ! Assign each of Olson surface type to an LSM surface type
    ! Mapping from Olson to LSM is based on BATS code
    ! NB: Ols2LSM(i) = Olson type i
    do Ols_idx=1,19
       Ols2LSM(Ols_idx)=mss_val_ntg
    end do ! end loop over idx
    Ols2LSM(20)=3                                                     
    Ols2LSM(21)=3                                                     
    Ols2LSM(22)=3                                                     
    Ols2LSM(23)=6                                                    
    Ols2LSM(24)=8                                                    
    Ols2LSM(25)=9                                                     
    Ols2LSM(26)=9                                                     
    Ols2LSM(27)=7             
    Ols2LSM(28)=10                                                     
    Ols2LSM(29)=10                                                     
    Ols2LSM(30)=24                                                     
    Ols2LSM(31)=26                                                     
    Ols2LSM(32)=12                                                     
    Ols2LSM(33)=10                                                     
    Ols2LSM(34)=mss_val_ntg                                                 
    Ols2LSM(35)=mss_val_ntg                                                 
    Ols2LSM(36)=28                                                    
    Ols2LSM(37)=25                                                    
    Ols2LSM(38)=23                                                    
    Ols2LSM(39)=23                                                    
    Ols2LSM(40)=17                                                     
    Ols2LSM(41)=18                                                     
    Ols2LSM(42)=17                                                     
    Ols2LSM(43)=12                                                     
    Ols2LSM(44)=28                                                    
    Ols2LSM(45)=28                                                    
    Ols2LSM(46)=20                                                    
    Ols2LSM(47)=20                                                    
    Ols2LSM(48)=20         
    Ols2LSM(49)=22                                                    
    Ols2LSM(50)=2                                                     
    Ols2LSM(51)=22                                                    
    Ols2LSM(52)=22         
    Ols2LSM(53)=19                                                     
    Ols2LSM(54)=19                                                     
    Ols2LSM(55)=15          
    Ols2LSM(56)=16         
    Ols2LSM(57)=15         
    Ols2LSM(58)=16          
    Ols2LSM(59)=21                                                    
    Ols2LSM(60)=6                                                    
    Ols2LSM(61)=4                                                     
    Ols2LSM(62)=13          
    Ols2LSM(63)=14                                                    
    Ols2LSM(64)=20                                                    
    Ols2LSM(65)=0         
    Ols2LSM(66)=0         
    Ols2LSM(67)=0         
    Ols2LSM(68)=0         
    Ols2LSM(69)=2                                                     
    Ols2LSM(70)=1                                                   
    Ols2LSM(71)=22                                                    
    Ols2LSM(72)=27                                                    
    Ols2LSM(73)=0         
    do Ols_idx=74,100
       Ols2LSM(Ols_idx)=mss_val_ntg
    end do ! end loop over idx
    
    ! Transform Olson surface types to LSM surface types
    err_nbr=0
    do lat_in_idx=1,lat_in_nbr                                                 
       do lon_in_idx=1,lon_in_nbr                                                 
          if (sfc_typ_Ols_prc(lon_in_idx,lat_in_idx)==0) then
             sfc_typ_LSM_in(lon_in_idx,lat_in_idx)=0
          else
             sfc_typ_LSM_in(lon_in_idx,lat_in_idx)=Ols2LSM(sfc_typ_Ols_prc(lon_in_idx,lat_in_idx))
          end if              ! endif
          if (sfc_typ_LSM_in(lon_in_idx,lat_in_idx) > sfc_typ_LSM_max.or.sfc_typ_LSM_in(lon_in_idx,lat_in_idx) < 0) then
             err_nbr=err_nbr+1
             lon_in_idx_err=lon_in_idx
             lat_in_idx_err=lat_in_idx
          end if
       end do ! end loop over lon
    end do ! end loop over lat
    if (err_nbr > 0) then
       write (6,"(a,i4,a,a,i3,a,i3,i3)") &
            sbr_nm//"(): Encountered ",err_nbr," errors translating Olson types to LSM types on Olson grid.", &
            "Undefined LSM surface type ",sfc_typ_LSM_in(lon_in_idx_err,lat_in_idx_err), &
            " at lon_in_idx,lat_in_idx = ",lon_in_idx_err,lat_in_idx_err
       stop
    end if ! endif err
    
    ! Get overlap locations and weights for mapping input grid to output grid
    call map_ovr_wgt_drv( &
         lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
         lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
         ovr_nbr_max,area_out, & ! I
         ovr_lat_idx,ovr_lon_idx,ovr_nbr,ovr_wgt) ! O

    nn_cnt=0
    ! Process each cell on output grid
    do lat_out_idx=1,lat_out_nbr
       do lon_out_idx=1,lon_out_nbr
          ! Sum overlap weights by surface type
          do sfc_typ_LSM_idx=0,sfc_typ_LSM_max ! Loop starts with 0
             sfc_typ_wgt(sfc_typ_LSM_idx)=0.0 ! [frc] Weight of surface type in gridcell
          end do ! end loop over LSM surface type
          ovr_nbr_crr=ovr_nbr(lon_out_idx,lat_out_idx)
          if (ovr_nbr_crr==0) write (6,"(a,i3)") sbr_nm//"(): # overlapping input gridcells = ",ovr_nbr_crr
          do ovr_idx=1,ovr_nbr_crr ! Overlap cell index
             lon_in_idx=ovr_lon_idx(lon_out_idx,lat_out_idx,ovr_idx) ! [idx] lon index (input grid) of overlap cell
             lat_in_idx=ovr_lat_idx(lon_out_idx,lat_out_idx,ovr_idx) ! [idx] lat index (input grid) of overlap cell
             sfc_typ_LSM_crr=sfc_typ_LSM_in(lon_in_idx,lat_in_idx) ! [enm] Surface type of current gridcell
             sfc_typ_wgt(sfc_typ_LSM_crr)=sfc_typ_wgt(sfc_typ_LSM_crr)+ovr_wgt(lon_out_idx,lat_out_idx,ovr_idx) ! [frc] Weight of surface type in gridcell
          end do ! end loop over overlapping cells
          
          ! Rank non-zero weights by surface type 
          ! sfc_typ_wgt_1st is the most extensive surface type, sfc_typ_wgt_2nd is the second most extensive surface type
          ! NB: LSM code is rather ambiguous whether these ranks are indices or values
          call rnk_vec_LSM(sfc_typ_wgt,sfc_typ_LSM_max,mss_val_ntg, & ! I
               sfc_typ_rnk_1st,sfc_typ_rnk_2nd) ! O
          if (sfc_typ_rnk_1st==mss_val_ntg) then
             write (6,"(a,i3,a,f9.5)") sbr_nm//"(): ERROR occured at lon_out_grd(",lon_out_idx,") = ",lon_out_grd(lon_out_idx)
             write (6,"(a,i3,a,f9.5)") sbr_nm//"(): ERROR occured at lat_out_grd(",lat_out_idx,") = ",lat_out_grd(lat_out_idx)
             write (6,"(a,i3)") sbr_nm//"(): # overlapping input gridcells = ",ovr_nbr_crr
             write (6,"(a)") "lon_in_idx lat_in_idx sfc_typ_LSM_in ovr_wgt" 
             do ovr_idx=1,ovr_nbr_crr !overlap cell index
                lon_in_idx=ovr_lon_idx(lon_out_idx,lat_out_idx,ovr_idx) ! [idx] lon index (input grid) of overlap cell
                lat_in_idx=ovr_lat_idx(lon_out_idx,lat_out_idx,ovr_idx) ! [idx] lat index (input grid) of overlap cell
                sfc_typ_LSM_crr=sfc_typ_LSM_in(lon_in_idx,lat_in_idx) ! [enm] Surface type of current gridcell
                sfc_typ_wgt(sfc_typ_LSM_crr)=sfc_typ_wgt(sfc_typ_LSM_crr)+ovr_wgt(lon_out_idx,lat_out_idx,ovr_idx) ! [frc] Weight of surface type in gridcell
                write (6,"(i3,1x,i3,1x,i2,1x,f9.6)") lon_in_idx,lat_in_idx, &
                     sfc_typ_LSM_crr,ovr_wgt(lon_out_idx,lat_out_idx,ovr_idx)
             end do ! end loop over overlapping cells
             stop
          endif               ! endif err
          
          ! Set surface type as:
          ! If fractional land = 0: cell = ocean
          ! If fractional land > 0: cell = land
          ! a. Use most frequent surface type based on area of overlap unless this is ocean 
          ! In this case, input grid says ocean but output grid wants land
          ! b. So use next most extensive surface type so long as is not ocean
          ! c. If this is ocean, or if there is none,
          !    If flg_soi_nn is FALSE (original method) use swamp/marsh surface type
          !    This algorithm, originally used by LSM, fails for changes in mean sea 
          !    level because new land may appear that should not be swamp/marsh
          !    A nearest neighbor algorithm would be much more appropriate
          !    Same is true of lai_get(), pft_get(), soi_txt_get()
          !    If flg_soi_nn is TRUE (MF added) then set undefined land surface type 
          !    to land type of nearest (originally-defined) land gridcell

          if (lnd_frc_out(lon_out_idx,lat_out_idx)==0.0) then ! Ocean
             sfc_typ_LSM_out(lon_out_idx,lat_out_idx)=0
          else ! not ocean
             if (sfc_typ_rnk_1st /= 0) then
                sfc_typ_LSM_out(lon_out_idx,lat_out_idx)=sfc_typ_rnk_1st 
             else ! Dominant type is ocean, consider second-ranked type...
                if (sfc_typ_rnk_2nd==0.or.sfc_typ_rnk_2nd==mss_val_ntg) then
                   
                   ! Nearest Neighbor Method
                   if (flg_soi_nn) then
                      
                      distance_min = 100000 ! [km] Initialized to high value
                      lat1_tmp = lat_out_grd(lat_out_idx) / degtorad ! [nbr] radians
                      lon1_tmp = lon_out_grd(lon_out_idx) / degtorad ! [nbr] radians
                      
                      ! Find nearest neighbor on input grid
                      do lat_idx2=1,lat_in_nbr
                         do lon_idx2=1,lon_in_nbr
                            if ( (sfc_typ_LSM_in(lon_idx2,lat_idx2).ne.0) ) then
                               !    (sfc_typ_LSM_in(lon_idx2,lat_idx2).ne.1) )then
                               
                               lat2_tmp = lat_in_grd(lat_idx2) / degtorad ! [nbr] radians
                               lon2_tmp = lon_in_grd(lon_idx2) / degtorad ! [nbr] radians
                            
                               ! Great Circle Distance Formula (kilometers)
                               distance_tmp = 6378.7*acos((sin(lat1_tmp)*sin(lat2_tmp))+ &
                                    (cos(lat1_tmp)*cos(lat2_tmp)*cos(abs(lon2_tmp-lon1_tmp))))
                               
                               if (distance_tmp.lt.distance_min) then
                                  distance_min = distance_tmp
                                  lat_idx_min_dst = lat_idx2
                                  lon_idx_min_dst = lon_idx2
                               endif ! endif
                            endif ! endif
                         enddo
                      enddo
                      sfc_typ_LSM_out(lon_out_idx,lat_out_idx)=sfc_typ_LSM_in(lon_idx_min_dst,lat_idx_min_dst)
                      nn_cnt = nn_cnt + 1
                   else 
                      sfc_typ_LSM_out(lon_out_idx,lat_out_idx)=28   ! Original method: set gridcell to wetland/marsh
                      
                   endif ! endif nearest neighbor
                   
                else
                   sfc_typ_LSM_out(lon_out_idx,lat_out_idx)=sfc_typ_rnk_2nd
                end if ! endif
             end if ! endif dominant type is ocean
          end if ! endif lnd_frc_out
          
          ! Perform any arbitrary post-processing desired for dust data
          ! Operations in this block are not performed by LSM
          ! Thus each remapping causes discrepancy between LSM and dust datasets
          if (flg_set_crp_frs) then
             ! Set warm/cool crops to forest crops
             ! This considerably reduces bare ground fraction in northern Europe, midwest US in winter
             sfc_typ_LSM_crr=sfc_typ_LSM_out(lon_out_idx,lat_out_idx) ! [enm] Surface type of current gridcell
             if (sfc_typ_LSM_crr >= 23 .and. sfc_typ_LSM_crr <=26) crp_frs_rpl_nbr=crp_frs_rpl_nbr+1 ! [nbr] Number of crop<-->forest replacements
             ! Change cool irrigated crop to cool forest crop
             if (sfc_typ_LSM_crr==23) sfc_typ_LSM_out(lon_out_idx,lat_out_idx)=15
             ! Change cool crop to cool forest crop
             if (sfc_typ_LSM_crr==24) sfc_typ_LSM_out(lon_out_idx,lat_out_idx)=15
             ! Change warm irrigated crop to warm forest crop
             if (sfc_typ_LSM_crr==25) sfc_typ_LSM_out(lon_out_idx,lat_out_idx)=16
             ! Change warm crop to warm forest crop
             if (sfc_typ_LSM_crr==26) sfc_typ_LSM_out(lon_out_idx,lat_out_idx)=16
          endif ! endif flg_set_crp_frs
          
          ! Sanity checks
          sfc_typ_LSM_crr=sfc_typ_LSM_out(lon_out_idx,lat_out_idx) ! [enm] Surface type of current gridcell
          if (sfc_typ_LSM_crr < 0.or.sfc_typ_LSM_crr > sfc_typ_LSM_max) then
             write (6,"(a,i10,a,2(i3,a1))") &
                  "ERROR sfc_typ_get() reports invalid LSM surface type = ", &
                  sfc_typ_LSM_crr, &
                  " at (lon_out_idx,lat_out_idx) = (",lon_out_idx,",",lat_out_idx,")"
             stop
          end if ! endif
          if ( &
               (lnd_frc_out(lon_out_idx,lat_out_idx)==0.0.and.sfc_typ_LSM_crr /= 0).or. & ! Land fraction says ocean but sfc_typ_LSM says land
               (lnd_frc_out(lon_out_idx,lat_out_idx) > 0.0.and.sfc_typ_LSM_crr==0) & ! Land fraction says land but sfc_typ_LSM says ocean
               ) then
             write (6,"(a,i2,a,f9.6,a,2(i3,a1))") "ERROR sfc_type_get() reports sfc_typ_LSM_out = ", &
                  sfc_typ_LSM_crr, &
                  ", lnd_frc_out = ",lnd_frc_out(lon_out_idx,lat_out_idx), &
                  " at (lon_out_idx,lat_out_idx) = (",lon_out_idx,",",lat_out_idx,")"
             stop
          end if ! endif
          
       end do ! end loop over lon
    end do ! end loop over lat
    
    if (flg_set_crp_frs) then ! [flg] Set warm/cool crops, irrigated crops to forest crops
       write (6,"(2a,i6,a)") prg_nm(1:ftn_strlen(prg_nm)), & 
            ": INFO sfc_typ_get() replaced crops by forest crops in ", &
            crp_frs_rpl_nbr," gridcells"
    endif ! endif flg_set_crp_frs

    if (nn_cnt > 0) then
       write (6,"(a,i4,a)") "DEBUG: sfc_typ_get() reports ",nn_cnt," points set to nearest neighbor"
    endif         ! endif

    
    ! Compute global area of each LSM surface type 
    do sfc_typ_LSM_idx=0,sfc_typ_LSM_max ! Loop starts with 0
       area_sfc_typ_in(sfc_typ_LSM_idx)=0.0
    end do ! end loop over sfc_typ_LSM
    do lat_in_idx=1,lat_in_nbr
       do lon_in_idx=1,lon_in_nbr
          sfc_typ_LSM_crr=sfc_typ_LSM_in(lon_in_idx,lat_in_idx) ! [enm] Surface type of current gridcell
          area_sfc_typ_in(sfc_typ_LSM_crr)=area_sfc_typ_in(sfc_typ_LSM_crr)+area_in(lon_in_idx,lat_in_idx)
       end do ! end loop over lon
    end do ! end loop over lat
    
    ! Compute global area of each LSM surface type on output grid
    do sfc_typ_LSM_idx=0,sfc_typ_LSM_max ! Loop starts with 0
       area_sfc_typ_out(sfc_typ_LSM_idx)=0.0
    end do ! end loop over sfc_typ_LSM
    do lat_out_idx=1,lat_out_nbr
       do lon_out_idx=1,lon_out_nbr
          sfc_typ_LSM_crr=sfc_typ_LSM_out(lon_out_idx,lat_out_idx) ! [enm] Surface type of current gridcell
          area_sfc_typ_out(sfc_typ_LSM_crr)=area_sfc_typ_out(sfc_typ_LSM_crr)+area_out(lon_out_idx,lat_out_idx)
       end do ! end loop over lon
    end do ! end loop over lat
    
    ! Compare areas
    if (dbg_lvl==dbg_crr) then
       open (unit_dgn,file=fl_dgn,status="unknown",iostat=rcd)
       if (rcd /= 0) write (6,"(4a,i4)") prg_nm(1:ftn_strlen(prg_nm)), & 
            ": ERROR sfc_typ_get() unable to open ",fl_dgn(1:ftn_strlen(fl_dgn))
       write (unit_dgn,"(a)") sbr_nm//"() Surface Type Diagnostics File"
       write (unit_dgn,"(a,f16.3,a)") "Output grid cell area nearest the Equator: ",area_out(1,lat_out_nbr/2)*1.0e-6, " km2"
       write (unit_dgn,"(a,20x,a,1x,a)") "Vegetation type"," Input grid area","Output grid area"
       write (unit_dgn,"(33x,7x,a,8x,a)") "km2","km2"
       do sfc_typ_LSM_idx=0,sfc_typ_LSM_max ! Loop starts with 0
          write (unit_dgn,"(a35,f16.3,f17.3)") vgt_sng(sfc_typ_LSM_idx), &
               area_sfc_typ_in(sfc_typ_LSM_idx)*1.0e-6,area_sfc_typ_out(sfc_typ_LSM_idx)*1.0e-6
       end do ! end loop over sfc_typ_LSM
       close (unit_dgn)
       write (6,"(a,a)") "Wrote diagnostic land surface data to ",fl_dgn(1:ftn_strlen(fl_dgn))
    endif                     ! endif dbg

    ! Write sfc_typ data to text file using Olson format
    if (.false.) then
       
       LSM2Ols(0)=65
       LSM2Ols(1)=70
       LSM2Ols(2)=50
       LSM2Ols(3)=20
       LSM2Ols(4)=61
       LSM2Ols(5)=61
       LSM2Ols(6)=23
       LSM2Ols(7)=27
       LSM2Ols(8)=24
       LSM2Ols(9)=25
       LSM2Ols(10)=28
       LSM2Ols(11)=28
       LSM2Ols(12)=32
       LSM2Ols(13)=62
       LSM2Ols(14)=63
       LSM2Ols(15)=55
       LSM2Ols(16)=56
       LSM2Ols(17)=40
       LSM2Ols(18)=41
       LSM2Ols(19)=53
       LSM2Ols(20)=46
       LSM2Ols(21)=59
       LSM2Ols(22)=49
       LSM2Ols(23)=38
       LSM2Ols(24)=30
       LSM2Ols(25)=37
       LSM2Ols(26)=31
       LSM2Ols(27)=72
       LSM2Ols(28)=36


       ! Transform LSM (out) surface types to Olson surface types and write to output file in 0.5 x 0.5 resolution
       ! Assumes output resolution of 2x2 and map_typ_out=11
       open (unit_dgn,file=fl_Ols,status="unknown",iostat=rcd)
       if (rcd /= 0) write (6,"(4a,i4)") prg_nm(1:ftn_strlen(prg_nm)), & 
            ": ERROR sfc_typ_get() unable to open ",fl_Ols(1:ftn_strlen(fl_Ols))
       
       do lat_out_idx=lat_out_nbr,1,-1
          do temp_idx=1,4
             do lon_out_idx=1,lon_out_nbr   
                sfc_typ_crr = LSM2Ols(sfc_typ_LSM_out(lon_out_idx,lat_out_idx))
                write (unit_dgn,"(i2)") sfc_typ_crr
                write (unit_dgn,"(i2)") sfc_typ_crr
                write (unit_dgn,"(i2)") sfc_typ_crr
                write (unit_dgn,"(i2)") sfc_typ_crr
             end do ! end loop over lon
          end do          
       end do ! end loop over lat
       close (unit_dgn)
       write (6,"(a,1x,a)") "Wrote LSM data in Olson format to ",fl_Ols(1:ftn_strlen(fl_Ols))
    endif ! endif write output file


    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Exiting sfc_typ_get()"
    return 
  end subroutine sfc_typ_get
  
  subroutine sfc2vgt( & ! [sbr] LAI, VAI prescribed for LSM surface type
       lat,lat_nbr,lon,lon_nbr,time_nbr, & ! I
       sfc_typ, & ! I
       lai_lsm,vai_lsm) ! O
    ! Purpose: Compute and return seasonal cycle of vegetation properties
    ! sfc2vgt() is called by bds_prc()
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use dstlsm,only:pln_typ,pln_frc,gai,tai ! [mdl] LSM data
    use map_grd ! [mdl] Map grids and regridding
    implicit none
    ! Parameters
    ! Commons
    ! Input
    integer,intent(in)::lat_nbr ! [nbr] Dimension size
    integer,intent(in)::lon_nbr ! [nbr] Dimension size
    integer,intent(in)::time_nbr ! [nbr] Dimension size
    integer,intent(in)::sfc_typ(lon_nbr,lat_nbr) ! [enm] Surface type code
    real,intent(in)::lat(lat_nbr) ! [dgr] Centered latitudes
    real,intent(in)::lon(lon_nbr) ! [dgr] Centered longitudes
    ! Output
    real,intent(out)::lai_lsm(lon_nbr,lat_nbr,time_nbr) ! O [m2 m-2] Leaf area index, one-sided
    real,intent(out)::vai_lsm(lon_nbr,lat_nbr,time_nbr) ! O [m2 m-2] Vegetation area index, one-sided
    ! Local
    integer lat_idx           ! [idx] Counting index
    integer lon_idx           ! [idx] Counting index
    integer pln_typ_idx       ! [idx] Plant type index
    integer sfc_typ_idx       ! [idx] Surface type index
    integer sgs_idx           ! [idx] Surface sub-gridscale index
    integer time_idx          ! [idx] Counting index
    integer time_idx_NH       ! [idx] Counting index into LSM block data
    ! Main Code
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Entering sfc2vgt()"
    ! Initialize array
    lai_lsm=0.0 ! [m2 m-2]
    vai_lsm=0.0 ! [m2 m-2]
    ! Sanity checks
    if (time_nbr /= 12) stop "ERROR: time_nbr /= 12 in sfc2vgt()"
    do lat_idx=1,lat_nbr
       do lon_idx=1,lon_nbr
          ! Store surface blend for current gridpoint
          sfc_typ_idx=sfc_typ(lon_idx,lat_idx)
          ! Vegetation information not available over ocean
          if (sfc_typ_idx /= 0) then
             do sgs_idx=1,3
                ! Non-vegetated land surfaces have pln_typ=14, ocean has pln_typ=0
                pln_typ_idx=pln_typ(sfc_typ_idx,sgs_idx)
                do time_idx_NH=1,time_nbr
                   ! LSM phenology data assumes index corresponds to NH month
                   if (lat(lat_idx) > 0.0) then
                      time_idx=time_idx_NH
                   else
                      time_idx=mod(time_idx_NH+time_nbr/2,time_nbr)
                      if (time_idx==0) time_idx=time_nbr
                   endif      ! endif SH
                   lai_lsm(lon_idx,lat_idx,time_idx)= & ! [m2 m-2]
                        lai_lsm(lon_idx,lat_idx,time_idx)+ & ! [m2 m-2]
                        pln_frc(sfc_typ_idx,sgs_idx)* & ! [frc]
                        gai(pln_typ_idx,time_idx_NH) ! [m2 m-2]
                   vai_lsm(lon_idx,lat_idx,time_idx)= & ! [m2 m-2]
                        vai_lsm(lon_idx,lat_idx,time_idx)+ & ! [m2 m-2]
                        pln_frc(sfc_typ_idx,sgs_idx)* & ! [frc]
                        tai(pln_typ_idx,time_idx_NH) ! [m2 m-2]
                end do ! end loop over time
             end do ! end loop over number of possible plant types
          endif               ! endif sfc_typ_idx /= 0
       end do ! end loop over lon
    end do ! end loop over lat
    
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Exiting sfc2vgt()"
    return 
  end subroutine sfc2vgt
  
  subroutine sfc2dps( & ! [sbr] Surface roughness prescribed for LSM surface type
       lat,lat_nbr,lon,lon_nbr, & ! I
       sfc_typ, & ! I
       rgh_mmn,zpd_lsm) ! O
    ! Purpose: Compute and return seasonal cycle of deposition properties
    ! sfc2dps() is called by bds_prc()
    ! Note that sfc2dps() is similar to sfc2vgt() except that sfc2dps() produces 
    ! gridcell averaged properties used in deposition calculation.
    ! Since deposition may occur over any surface, we must take special care
    ! to make sure the returned parameters are valid over, e.g., ocean and ice.
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use dstlsm,only:pln_typ,z0mvt,pln_frc,zpdvt ! [mdl] LSM data
    implicit none
    ! Parameters
    ! Commons
    ! Input
    integer lat_nbr ! [nbr] Number of latitudes
    integer lon_nbr ! [nbr] Number of longitudes
    integer sfc_typ(lon_nbr,lat_nbr) ! [enm] Surface type code
    real lat(lat_nbr) ! [dgr] Centered latitudes
    real lon(lon_nbr) ! [dgr] Centered longitudes
    ! Output
    real rgh_mmn(lon_nbr,lat_nbr) ! [m] Roughness length momentum
    real zpd_lsm(lon_nbr,lat_nbr) ! [m] Zero plane displacement height
    ! Local
    integer lat_idx           ! [idx] Counting index
    integer lon_idx           ! [idx] Counting index
    integer pln_typ_idx       ! [idx] Plant type index
    integer sfc_typ_idx       ! [idx] Surface type index
    integer sgs_idx           ! [idx] Surface sub-gridscale index
    real rlm_crr              ! [m] Roughness length of current sub-gridscale
    ! Main Code
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Entering sfc2dps()"
    ! Initialize array
    rgh_mmn(:,:)=0.0 ! [m]
    zpd_lsm(:,:)=0.0 ! [m]
    ! Sanity checks
    do lat_idx=1,lat_nbr
       do lon_idx=1,lon_nbr
          ! Store surface blend for current gridpoint
          sfc_typ_idx=sfc_typ(lon_idx,lat_idx)
          if (sfc_typ_idx==0) then ! Ocean
             rgh_mmn(lon_idx,lat_idx)=0.001 ! [m] Bon96 p. 59
             zpd_lsm(lon_idx,lat_idx)=0.0 ! [m]
          else if(sfc_typ_idx==1) then ! Land ice
             rgh_mmn(lon_idx,lat_idx)=0.05 ! [m] Bon96 p. 59
             zpd_lsm(lon_idx,lat_idx)=0.0 ! [m]
          else                ! Normal land
             do sgs_idx=1,3
                ! Non-vegetated land surfaces have pln_typ=14, ocean has pln_typ=0
                pln_typ_idx=pln_typ(sfc_typ_idx,sgs_idx)
                if (pln_typ_idx==14) then ! Bare ground
                   rlm_crr=0.05 ! [m] Bon96 p. 59
                else if (pln_typ_idx > 0) then ! Regular plant type
                   rlm_crr=z0mvt(pln_typ_idx) ! [m]
                else          ! Presumably ocean snuck through
                   stop "ERROR: pln_typ_idx==0 in sfc2dps()"
                endif         ! endif
                ! NB: SeP97 apparently uses rlm but rlh may be more correct for aerosol
                rgh_mmn(lon_idx,lat_idx)= & ! [m]
                     rgh_mmn(lon_idx,lat_idx)+ &
                     pln_frc(sfc_typ_idx,sgs_idx)* & ! [frc]
                     rlm_crr  ! [m]
                zpd_lsm(lon_idx,lat_idx)= & ! [m]
                     zpd_lsm(lon_idx,lat_idx)+ &
                     pln_frc(sfc_typ_idx,sgs_idx)* & ! [frc]
                     zpdvt(pln_typ_idx) ! [m]
             end do ! end loop over number of possible plant types
          endif               ! endif normal land
       end do ! end loop over lon
    end do ! end loop over lat
    
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Exiting sfc2dps()"
    return 
  end subroutine sfc2dps
  
end module sfc_typ_mdl ! [mdl] Surface type and derived properties
