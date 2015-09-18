! $Id$ -*-f90-*-

! Purpose: Assimilate model fields

! Usage: 
! use asm_mdl ! [mdl] Assimilate model fields

module asm_mdl ! [mdl] Assimilate model fields
  implicit none
  public::asm_drv ! [sbr] Assimilation driver
  private::mdl_atm_get ! [sbr] Retrieve and rebin atmospheric model data
  private::hyd_get ! [sbr] Retrieve and rebin hydrologic data 
  private::mdl_lsm_get ! [sbr] Retrieve and rebin volumetric water content
  private::lsm_get_sfc_flw ! [sbr] Retrieve and rebin surface flow data QOVER
  
contains
  
  subroutine asm_drv( & ! [sbr] Assimilation driver
       lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr,time_out_nbr, & ! I
       yr_nbr,yr_srt, & ! I
       area_out,lnd_frc_out,lnd_msk_out, & ! I
       fsh_fct,odxc_tms,src_flg,src_frq,src_str, & ! O 
       odxc_mdl,src_str_old,vwc_sfc,sfc_flw) ! O 
    ! Purpose: Assimilate data from various model and observational datasets
    ! asm_drv() is called by bds_prc()
    use asm_obs,only:toms_get ! [mdl] Assimilate observations
    use bds_ctl ! [mdl] Control variables drc_in,drc_out,hgt_dlt_msl
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use map_cst ! [mdl] Constants used in map routines
    use map_grd ! [mdl] Map grids and regridding
    use utl_mdl,only:mnt_chk ! [mdl] Utility functions (date_time_get,mnt_chk...)
    implicit none
    ! Parameters
    integer,parameter::lon_in_nbr_max=360 ! [nbr] Maximum number of longitudes
    integer,parameter::lat_in_nbr_max=360 ! [nbr] Maximum number of latitudes
    integer,parameter::time_in_nbr_max=1 ! [nbr] Maximum number of records (months) in input file
    ! Commons
    ! Input
    integer lon_out_nbr ! [nbr] Number of longitudes
    integer lat_out_nbr ! [nbr] Number of latitudes
    real area_out(lon_out_nbr,lat_out_nbr) ! [m2] Area of gridcells 
    real lon_out_grd(lon_out_nbr+1) ! [dgr] Interface longitudes
    real lat_out_grd(lat_out_nbr+1) ! [dgr] Interface latitudes
    integer time_out_nbr ! [nbr] Dimension size
    integer yr_nbr ! [nbr] Number of years
    integer yr_srt ! [yr] Starting year in YYYY format
    integer lnd_msk_out(lon_out_nbr,lat_out_nbr) ! [flg] Land mask (integer 0 or 1)
    real lnd_frc_out(lon_out_nbr,lat_out_nbr) ! [frc] Land fraction of gridcell 
    ! Output
    integer src_flg(lon_out_nbr,lat_out_nbr) ! [flg] Source flag
    real fsh_fct(lon_out_nbr,lat_out_nbr) ! [frc] Efficiency factor
    real src_frq(lon_out_nbr,lat_out_nbr,time_out_nbr) ! [frc] Source frequency
    real src_str(lon_out_nbr,lat_out_nbr,time_out_nbr) ! [frc] Source strength
    real odxc_mdl(lon_out_nbr,lat_out_nbr,time_out_nbr) ! [frc] Optical depth from model
    real odxc_tms(lon_out_nbr,lat_out_nbr,time_out_nbr) ! [frc] Optical depth from TOMS
    real src_str_old(lon_out_nbr,lat_out_nbr,time_out_nbr) ! [frc] Source strength
    real vwc_sfc(lon_out_nbr,lat_out_nbr,time_out_nbr) ! [m3 m-3] Volumetric water content
    real sfc_flw(lon_out_nbr,lat_out_nbr) ! [mm s-1] Surface runoff (QOVER)
    
    ! Local
    character(80)::fl_in        ! [sng] Input file
    character(2)::mth_sng       ! [sng] Month in MM format
    character(4)::yr_sng        ! [sng] Year in YYYY format
    integer lat_in_nbr        ! [nbr] Number of latitudes
    integer lon_in_nbr        ! [nbr] Number of longitudes
    integer time_in_nbr       ! [nbr] Number of records (months) in input file
    integer map_typ_in        ! [enm] Input map grid type
    integer ovr_nbr_max       ! [nbr] Maximum number of input cells which overlap any output cell
    logical mnt               ! [flg] Monotonicity flag
    real lat_in_grd(lat_in_nbr_max+1) ! [dgr] Interface latitudes
    real lon_in_grd(lon_in_nbr_max+1) ! [dgr] Interface longitudes
    
    ! Main Code
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Entering asm_drv()"
    
    ! Check for monotonicity
    mnt=mnt_chk(lon_out_grd,lon_out_nbr+1)
    if (.not.mnt) stop "lon_out_grd not monotonic in asm_drv()"
    mnt=mnt_chk(lat_out_grd,lat_out_nbr+1)
    if (.not.mnt) stop "lat_out_grd not monotonic in asm_drv()"
    ! Check for increasing monotonicity
    if (lon_out_grd(2) < lon_out_grd(1)) stop "lon_out_grd not increasing in asm_drv()"
    if (lat_out_grd(2) < lat_out_grd(1)) stop "lat_out_grd not increasing in asm_drv()"
    
    ! Initialize scalars
    lat_in_nbr=lat_in_nbr_max
    lon_in_nbr=lon_in_nbr_max
    time_in_nbr=time_in_nbr_max
    
    call mdl_atm_get( & ! [sbr] Retrieve and rebin atmospheric model data
         fl_atm, & ! I
         lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr,time_out_nbr, & ! I
         area_out, & ! I
         odxc_mdl,src_str_old) ! O 
    
    ! Get TOMS optical depth
    lat_in_nbr=180            ! TOMS 1 x 1.5 degree grid
    lon_in_nbr=288            ! TOMS 1 x 1.5 degree grid
    time_in_nbr=1             ! TOMS data is one timeslice (month) per file
    map_typ_in=map_lat_rgl_lon_Grn_ctr ! Latitudes are regular, Greenwich at center of first longitude bin
    ! call map_grd_read(fl_grd, & ! I
    !    lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr) ! O
    call map_grd_mk(lat_in_nbr,lon_in_nbr,map_typ_in, & ! I
         lat_in_grd,lon_in_grd) ! O
    call map_ovr_nbr_max_get(lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
         lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
         ovr_nbr_max) ! O
    call toms_get( & ! [sbr] Retrieve and process TOMS satellite data
         drc_tms,yr_nbr,yr_srt, & ! I
         lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr,time_in_nbr, & ! I
         lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr,time_out_nbr, & ! I
         area_out,lnd_frc_out,lnd_msk_out,ovr_nbr_max, & ! I
         fsh_fct,odxc_tms,src_flg,src_frq,src_str) ! O 
    
    ! Get land surface hydrology
#if 0
    ! fxm: New method---needs to be tested then activated
    ! Must add new variable to deal with 12-month sfc_flw since geomorphology wants climate mean
    call hyd_get( & ! [sbr] Retrieve and rebin hydrologic data 
         fl_hyd, & ! I
         lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr,time_out_nbr, & ! I
         area_out, & ! I
         vwc_sfc,sfc_flw) ! O 
#else /* !0 */

    ! fxm: Old method Should replace mdl_lsm_get() and lsm_get_sfc_flw() by single routine hyd_get()
    fl_in="/data/zender/map/lsm_clm01.nc"//char(0) ! [sng] Seed filename loop with January value
    lat_in_nbr=64             ! LSM T42
    lon_in_nbr=128            ! LSM T42
    time_in_nbr=1             ! LSM data is one timeslice (month) per file
    map_typ_in=map_lat_Gss_lon_Grn_ctr ! Latitudes are Gaussian, Greenwich at center of first longitude bin
    call map_grd_mk(lat_in_nbr,lon_in_nbr,map_typ_in, & ! I
         lat_in_grd,lon_in_grd) ! O
    call map_ovr_nbr_max_get(lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
         lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
         ovr_nbr_max) ! O
    call mdl_lsm_get(fl_in, & ! I
         lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr,time_in_nbr, & ! I
         lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr,time_out_nbr, & ! I
         area_out,ovr_nbr_max, & ! I
         vwc_sfc) ! O 

    ! djn: Get QOVER from different file call map_grd_mk()
    ! djn: Note: h2osoi has different # of dimension in lsmh_1x1.nc vs lsm_clm01.nc
    ! Process Sam Levis data on 1x1 grid
    fl_in="/data/zender/map/lsmh_1x1.nc"//char(0)
    lat_in_nbr=180 ! LSM 1 x 1 degree grid
    lon_in_nbr=360 ! LSM 1 x 1 degree grid
    time_in_nbr=1 ! LSM history data is one timeslice per file
    map_typ_in=map_lat_rgl_lon_180_wst ! Latitudes are Gaussian, Greenwich at center of first longitude bin
    call map_grd_mk(lat_in_nbr,lon_in_nbr,map_typ_in, & ! I
         lat_in_grd,lon_in_grd) ! O
    call map_ovr_nbr_max_get(lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
         lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
         ovr_nbr_max) ! O
    call lsm_get_sfc_flw( & ! [sbr] Retrieve and rebin surface flow data QOVER
         fl_in, & ! I
         lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr,time_in_nbr, & ! I
         lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr,time_out_nbr, & ! I
         area_out,ovr_nbr_max, & ! I
         sfc_flw) ! O 
    
#endif /* !0 */

    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Exiting asm_drv()"
    return 
  end subroutine asm_drv
  
  subroutine mdl_atm_get( & ! [sbr] Retrieve and rebin atmospheric model data
       fl_in, & ! I
       lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr,time_out_nbr, & ! I
       area_out, & ! I
       odxc_out,src_str_out) ! O 
    ! Purpose: Retrieve data from atmospheric model datasets
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use map_cst ! [mdl] Constants used in map routines
    use map_grd ! [mdl] Map grids and regridding
    use netcdf ! [mdl] netCDF interface
    use nf90_utl ! [mdl] netCDF utilities
    use sng_mdl,only:ftn_strlen,ftn_strstr ! [mdl] String manipulation
    use utl_mdl,only:mnt_ncr_chk ! [mdl] Utility functions (date_time_get,mnt_chk...)
    implicit none
    ! Parameters
    character(len=*),parameter::sbr_nm="mdl_atm_get" ! [sng] Subroutine name
    integer,parameter::lat_nbr_max=180 ! [nbr] Maximum number of latitudes
    integer,parameter::lon_nbr_max=360 ! [nbr] Maximum number of longitudes
    ! Commons
    ! Input
    character(len=*),intent(in)::fl_in        ! [sng] Input file
    integer,intent(in)::lat_out_nbr       ! [nbr] Number of latitudes
    integer,intent(in)::lon_out_nbr       ! [nbr] Number of longitudes
    integer,intent(in)::time_out_nbr      ! [nbr] Dimension size
    real,intent(in)::area_out(lon_out_nbr,lat_out_nbr) ! [m2] Area of gridcells
    real,intent(in)::lat_out_grd(lat_out_nbr+1) ! [dgr] Interface latitudes
    real,intent(in)::lon_out_grd(lon_out_nbr+1) ! [dgr] Interface longitudes
    ! Output
    real,intent(out)::odxc_out(lon_out_nbr,lat_out_nbr,time_out_nbr) ! [frc] Optical depth
    real,intent(out)::src_str_out(lon_out_nbr,lat_out_nbr,time_out_nbr) ! [frc] Source strength
    
    ! Locals with simple initialization
    integer::map_typ_in=map_lat_Gss_lon_Grn_ctr ! Latitudes are Gaussian, Greenwich at center of first longitude bin
    integer::rcd=nf90_noerr ! [enm] Return success code
    logical::flg_read_2d_grd=.false. ! [flg] Read/use un-necessary CLM2 2D grid fields
    logical::flg_CLM=.false. ! [flg] Model input has CLM2 coding
    real::mss_val_in=1.0e36 ! [frc] Missing value
    real::mss_val_out=1.0e36 ! [frc] Missing value
    
    ! Local
    integer lat_out_idx       ! [idx] Counting index for lat
    integer lon_out_idx       ! [idx] Counting index for lon
    integer ovr_idx           ! [idx] Counting index
    integer ovr_nbr(lon_out_nbr,lat_out_nbr) ! [nbr] Number of input gridcells which overlap each output gridcell
    integer ovr_nbr_max       ! [nbr] Maximum number of input cells which overlap any output cell
    integer ovr_nbr_crr       ! [nbr] Current number of overlapping gridcells
    integer time_out_idx      ! [idx] Counting index for time
    logical mnt_ncr           ! [flg] Monotonically increasing
    real edg_grd(4) ! [dgr] Grid edges (north, east, south, west)
    
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
    integer odxc_id           ! [id] Variable ID
    integer src_str_dmn_nbr ! [nbr] Number of dimensions in disk src_str field
    integer src_str_id        ! [id] Variable ID
    integer srt_lon_lat_time(3) ! [idx] Dimension offset indices
    integer time_dmn_id       ! [id] Dimension ID for time
    integer time_in_nbr       ! [nbr] Number of records (months) in input file
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
    
    real,dimension(:,:),allocatable::odxc_in ! [frc] Optical depth 
    real,dimension(:,:),allocatable::src_str_in ! [frc] Source strength 
    
    ! Main Code
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Entering "//sbr_nm
    
    ! Enough memory?
    if (time_out_nbr /= 12) stop "mdl_atm_get(): time_out_nbr /= 12"
    
    ! Check for monotonicity
    mnt_ncr=mnt_ncr_chk(lon_out_grd,lon_out_nbr+1)
    mnt_ncr=mnt_ncr_chk(lat_out_grd,lat_out_nbr+1)
    
    ! Initialize scalars
    mss_flg=.false.
    ! Initialize arrays
    odxc_out(:,:,:)=0.0 ! [frc] Optical depth
    src_str_out(:,:,:)=0.0 ! [frc] Source strength
    
    ! Read netCDF data
    rcd=nf90_noerr ! nf90_noerr == 0
    rcd=nf90_wrp_open(fl_in,NF90_NOWRITE,nc_id,sbr_nm=sbr_nm)
    ! Get dimension IDs
    rcd=nf90_wrp_inq_dimid(nc_id,"lat",lat_dmn_id)
    rcd=nf90_wrp_inq_dimid(nc_id,"lon",lon_dmn_id)
    rcd=nf90_wrp_inq_dimid(nc_id,"time",time_dmn_id)
    ! Get dimension sizes
    rcd=rcd+nf90_inquire_dimension(nc_id,lon_dmn_id,len=lon_in_nbr)
    rcd=rcd+nf90_inquire_dimension(nc_id,lat_dmn_id,len=lat_in_nbr)
    rcd=rcd+nf90_inquire_dimension(nc_id,time_dmn_id,len=time_in_nbr)
    if(rcd /= nf90_noerr) stop "Error retrieving dimension sizes"
    
    ! Enough memory? 
    if (lat_in_nbr > lat_nbr_max) stop "lat_in_nbr > lat_nbr_max in mdl_atm_get()"
    if (lon_in_nbr > lon_nbr_max) stop "lon_in_nbr > lon_nbr_max in mdl_atm_get()"
    if (time_in_nbr /= time_out_nbr) stop "time_in_nbr /= time_out_nbr in mdl_atm_get()"
    if (time_in_nbr /= 12) then
       write (6,"(2a,f12.6)") prg_nm(1:ftn_strlen(prg_nm)), & 
            ": ERROR "//sbr_nm//" reports time_in_nbr = ",time_in_nbr
       stop "time_in_nbr /= 1"
    endif                  ! endif err
    
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
    
    allocate(odxc_in(lon_in_nbr,lat_in_nbr),stat=rcd) ! [frc] Optical depth 
    if(rcd /= 0) stop "allocate() failed for odxc_in"
    allocate(src_str_in(lon_in_nbr,lat_in_nbr),stat=rcd) ! [frc] Source strength 
    if(rcd /= 0) stop "allocate() failed for src_str_in"
    
    ! Get variable IDs
    rcd=nf90_wrp_inq_varid(nc_id,"DSTODXC",odxc_id)
    rcd=nf90_wrp_inq_varid(nc_id,"mbl_bsn_fct",src_str_id)
    if(flg_CLM.and.flg_read_2d_grd) then
       rcd=nf90_wrp_inq_varid(nc_id,"LATIXY",lat_in_2d_id)
       rcd=nf90_wrp_inq_varid(nc_id,"LONGXY",lon_in_2d_id)
       rcd=nf90_wrp_inq_varid(nc_id,"EDGEN",edg_nrt_id)
       rcd=nf90_wrp_inq_varid(nc_id,"EDGEE",edg_est_id)
       rcd=nf90_wrp_inq_varid(nc_id,"EDGES",edg_sth_id)
       rcd=nf90_wrp_inq_varid(nc_id,"EDGEW",edg_wst_id)
    endif ! endif flg_CLM.and.flg_read_2d_grd
    ! Get number of dimensions
    rcd=rcd+nf90_inquire_variable(nc_id,src_str_id,ndims=src_str_dmn_nbr)
    if (src_str_dmn_nbr /= 2) then
       write (6,"(2a,i1,a)") prg_nm(1:ftn_strlen(prg_nm)), & 
            ": ERROR mbl_bsn_fct has ",src_str_dmn_nbr," dimensions"
       write (6,"(2a,i1,a)") prg_nm(1:ftn_strlen(prg_nm)), & 
            ": HINT Make sure mbl_bsn_fct has no time dimension"
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
    rcd=nf90_wrp(nf90_get_var(nc_id,src_str_id,src_str_in),"get_var src_str_in")
    ! Get missing value
    rcd=rcd+nf90_inquire_attribute(nc_id,src_str_id,"missing_value",attnum=mss_val_id)
    if (rcd==nf90_noerr) then
       mss_flg=.true. ! [flg] Data may have missing values
       rcd=rcd+nf90_get_att(nc_id,src_str_id,"missing_value",mss_val_in)
       if (mss_val_in /= 1.0e36) write (6,"(2a,f12.6)") prg_nm(1:ftn_strlen(prg_nm)), &
            ": WARNING "//sbr_nm//" reports mss_val_in = ",mss_val_in
    else
       mss_flg=.false. ! [flg] Data may have missing values
       rcd=nf90_noerr ! [enm] Return success code
    endif ! endif
    if(rcd /= nf90_noerr) stop "Error retrieving variable data"
    
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
    endif ! endif flg_CLM.and.flg_read_2d_grd
    
    mnt_ncr=mnt_ncr_chk(lat_in_grd,lat_in_nbr+1)
    if (.not.mnt_ncr) stop "ERROR: lat_in_grd not monotonically increasing in mdl_atm_get()"
    mnt_ncr=mnt_ncr_chk(lon_in_grd,lon_in_nbr+1)
    if (.not.mnt_ncr) stop "ERROR: lon_in_grd not monotonically increasing in mdl_atm_get()"
    
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
    
    ! Loop over months
    do time_out_idx=1,time_out_nbr
       srt_lon_lat_time=(/1,1,time_out_idx/) ! [idx] Dimension offset indices
       cnt_lon_lat_time=(/lon_in_nbr,lat_in_nbr,1/) ! [idx] Dimension sizes
       ! Get data
       rcd=nf90_wrp(nf90_get_var(nc_id,odxc_id,odxc_in, &
            start=srt_lon_lat_time,count=cnt_lon_lat_time),"get_var odxc_in")
       ! Rebin optical depth from input grid to output grid
       call map_rbn(odxc_in,mss_flg,mss_val_in, & ! I
            lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
            lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
            ovr_nbr_max,ovr_lat_idx,ovr_lon_idx,ovr_nbr,ovr_wgt, & ! I
            odxc_out(1,1,time_out_idx)) ! O
       
       ! Rebin source strength from input grid to output grid
       call map_rbn(src_str_in,mss_flg,mss_val_in, & ! I
            lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
            lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
            ovr_nbr_max,ovr_lat_idx,ovr_lon_idx,ovr_nbr,ovr_wgt, & ! I
            src_str_out(1,1,time_out_idx)) ! O
       
    end do ! end loop over mth
    
    ! Close file
    rcd=nf90_wrp_close(nc_id,fl_in,"Ingested",sbr_nm=sbr_nm) ! [fnc] Close file
    
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
    
    if (allocated(odxc_in)) deallocate(odxc_in,stat=rcd)
    if(rcd /= 0) stop "deallocate() failed for odxc_in"
    if (allocated(src_str_in)) deallocate(src_str_in,stat=rcd)
    if(rcd /= 0) stop "deallocate() failed for src_str_in"
    
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Exiting "//sbr_nm
    return 
  end subroutine mdl_atm_get
  
  subroutine hyd_get( & ! [sbr] Retrieve and rebin hydrologic data 
       fl_in, & ! I
       lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr,time_out_nbr, & ! I
       area_out, & ! I
       sfc_flw_out,vwc_sfc_out) ! O 
    ! Purpose: Retrieve hydrologic data from model datasets and rebin to output grid
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use map_cst ! [mdl] Constants used in map routines
    use map_grd ! [mdl] Map grids and regridding
    use netcdf ! [mdl] netCDF interface
    use nf90_utl ! [mdl] netCDF utilities
    use sng_mdl,only:ftn_strlen,ftn_strstr ! [mdl] String manipulation
    use utl_mdl,only:mnt_ncr_chk ! [mdl] Utility functions (date_time_get,mnt_chk...)
    implicit none
    ! Parameters
    character(len=*),parameter::sbr_nm="hyd_get" ! [sng] Subroutine name
    integer,parameter::lat_nbr_max=180 ! [nbr] Maximum number of latitudes
    integer,parameter::lon_nbr_max=360 ! [nbr] Maximum number of longitudes
    ! Commons
    ! Input
    character(len=*),intent(in)::fl_in        ! [sng] Input file
    integer,intent(in)::lat_out_nbr       ! [nbr] Number of latitudes
    integer,intent(in)::lon_out_nbr       ! [nbr] Number of longitudes
    integer,intent(in)::time_out_nbr      ! [nbr] Dimension size
    real,intent(in)::area_out(lon_out_nbr,lat_out_nbr) ! [m2] Area of gridcells
    real,intent(in)::lat_out_grd(lat_out_nbr+1) ! [dgr] Interface latitudes
    real,intent(in)::lon_out_grd(lon_out_nbr+1) ! [dgr] Interface longitudes
    ! Output
    real,intent(out)::sfc_flw_out(lon_out_nbr,lat_out_nbr,time_out_nbr) ! [kg m-2 s-1] Surface flow
    real,intent(out)::vwc_sfc_out(lon_out_nbr,lat_out_nbr,time_out_nbr) ! [m3 m-3] Volumetric water content
    
    ! Locals with simple initialization
    integer::map_typ_in=map_lat_Gss_lon_Grn_ctr ! Latitudes are Gaussian, Greenwich at center of first longitude bin
    integer::rcd=nf90_noerr ! [enm] Return success code
    logical::flg_read_2d_grd=.false. ! [flg] Read/use un-necessary CLM2 2D grid fields
    logical::flg_CLM=.false. ! [flg] Model input has CLM2 coding
    real::mss_val_in=1.0e36 ! [frc] Missing value
    real::mss_val_out=1.0e36 ! [frc] Missing value
    
    ! Local
    integer lat_out_idx       ! [idx] Counting index for lat
    integer lon_out_idx       ! [idx] Counting index for lon
    integer ovr_idx           ! [idx] Counting index
    integer ovr_nbr(lon_out_nbr,lat_out_nbr) ! [nbr] Number of input gridcells which overlap each output gridcell
    integer ovr_nbr_max       ! [nbr] Maximum number of input cells which overlap any output cell
    integer ovr_nbr_crr       ! [nbr] Current number of overlapping gridcells
    integer time_out_idx      ! [idx] Counting index for time
    logical mnt_ncr           ! [flg] Monotonically increasing
    real edg_grd(4) ! [dgr] Grid edges (north, east, south, west)
    
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
    integer sfc_flw_id           ! [id] Variable ID
    integer vwc_sfc_dmn_nbr ! [nbr] Number of dimensions in disk vwc_sfc field
    integer vwc_sfc_id        ! [id] Variable ID
    integer srt_lon_lat_time(3) ! [idx] Dimension offset indices
    integer time_dmn_id       ! [id] Dimension ID for time
    integer time_in_nbr       ! [nbr] Number of records (months) in input file
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
    
    real,dimension(:,:),allocatable::sfc_flw_in ! [kg m-2 s-1] Surface flow 
    real,dimension(:,:),allocatable::vwc_sfc_in ! [m3 m-3] Volumetric water content 
    
    ! Main Code
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Entering "//sbr_nm
    
    ! Enough memory?
    if (time_out_nbr /= 12) stop "hyd_get(): time_out_nbr /= 12"
    
    ! Check for monotonicity
    mnt_ncr=mnt_ncr_chk(lon_out_grd,lon_out_nbr+1)
    mnt_ncr=mnt_ncr_chk(lat_out_grd,lat_out_nbr+1)
    
    ! Initialize scalars
    mss_flg=.false.
    ! Initialize arrays
    sfc_flw_out(:,:,:)=0.0 ! [kg m-2 s-1] Surface flow
    vwc_sfc_out(:,:,:)=0.0 ! [m3 m-3] Volumetric water content
    
    ! Read netCDF data
    rcd=nf90_noerr ! nf90_noerr == 0
    rcd=nf90_wrp_open(fl_in,NF90_NOWRITE,nc_id,sbr_nm=sbr_nm)
    ! Get dimension IDs
    rcd=nf90_wrp_inq_dimid(nc_id,"lat",lat_dmn_id)
    rcd=nf90_wrp_inq_dimid(nc_id,"lon",lon_dmn_id)
    rcd=nf90_wrp_inq_dimid(nc_id,"time",time_dmn_id)
    ! Get dimension sizes
    rcd=rcd+nf90_inquire_dimension(nc_id,lon_dmn_id,len=lon_in_nbr)
    rcd=rcd+nf90_inquire_dimension(nc_id,lat_dmn_id,len=lat_in_nbr)
    rcd=rcd+nf90_inquire_dimension(nc_id,time_dmn_id,len=time_in_nbr)
    if(rcd /= nf90_noerr) stop "Error retrieving dimension sizes"
    
    ! Enough memory? 
    if (lat_in_nbr > lat_nbr_max) stop "lat_in_nbr > lat_nbr_max in hyd_get()"
    if (lon_in_nbr > lon_nbr_max) stop "lon_in_nbr > lon_nbr_max in hyd_get()"
    if (time_in_nbr /= time_out_nbr) stop "time_in_nbr /= time_out_nbr in hyd_get()"
    if (time_in_nbr /= 12) then
       write (6,"(2a,f12.6)") prg_nm(1:ftn_strlen(prg_nm)), & 
            ": ERROR "//sbr_nm//" reports time_in_nbr = ",time_in_nbr
       stop "time_in_nbr /= 1"
    endif                  ! endif err
    
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
    
    allocate(sfc_flw_in(lon_in_nbr,lat_in_nbr),stat=rcd) ! [kg m-2 s-1] Surface flow 
    if(rcd /= 0) stop "allocate() failed for sfc_flw_in"
    allocate(vwc_sfc_in(lon_in_nbr,lat_in_nbr),stat=rcd) ! [m3 m-3] Volumetric water content 
    if(rcd /= 0) stop "allocate() failed for vwc_sfc_in"
    
    ! Get variable IDs
    rcd=nf90_wrp_inq_varid(nc_id,"QOVER",sfc_flw_id)
    rcd=nf90_wrp_inq_varid(nc_id,"H2OSOI",vwc_sfc_id)
    if(flg_CLM.and.flg_read_2d_grd) then
       rcd=nf90_wrp_inq_varid(nc_id,"LATIXY",lat_in_2d_id)
       rcd=nf90_wrp_inq_varid(nc_id,"LONGXY",lon_in_2d_id)
       rcd=nf90_wrp_inq_varid(nc_id,"EDGEN",edg_nrt_id)
       rcd=nf90_wrp_inq_varid(nc_id,"EDGEE",edg_est_id)
       rcd=nf90_wrp_inq_varid(nc_id,"EDGES",edg_sth_id)
       rcd=nf90_wrp_inq_varid(nc_id,"EDGEW",edg_wst_id)
    endif ! endif flg_CLM.and.flg_read_2d_grd
    ! Get number of dimensions
    rcd=rcd+nf90_inquire_variable(nc_id,vwc_sfc_id,ndims=vwc_sfc_dmn_nbr)
    if (vwc_sfc_dmn_nbr /= 2) then
       write (6,"(2a,i1,a)") prg_nm(1:ftn_strlen(prg_nm)), & 
            ": ERROR mbl_bsn_fct has ",vwc_sfc_dmn_nbr," dimensions"
       write (6,"(2a,i1,a)") prg_nm(1:ftn_strlen(prg_nm)), & 
            ": HINT Make sure mbl_bsn_fct has no time dimension"
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
    rcd=nf90_wrp(nf90_get_var(nc_id,vwc_sfc_id,vwc_sfc_in),"get_var vwc_sfc_in")
    ! Get missing value
    rcd=rcd+nf90_inquire_attribute(nc_id,vwc_sfc_id,"missing_value",attnum=mss_val_id)
    if (rcd==nf90_noerr) then
       mss_flg=.true. ! [flg] Data may have missing values
       rcd=rcd+nf90_get_att(nc_id,vwc_sfc_id,"missing_value",mss_val_in)
       if (mss_val_in /= 1.0e36) write (6,"(2a,f12.6)") prg_nm(1:ftn_strlen(prg_nm)), &
            ": WARNING "//sbr_nm//" reports mss_val_in = ",mss_val_in
    else
       mss_flg=.false. ! [flg] Data may have missing values
       rcd=nf90_noerr ! [enm] Return success code
    endif ! endif
    if(rcd /= nf90_noerr) stop "Error retrieving variable data"
    
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
    endif ! endif flg_CLM.and.flg_read_2d_grd
    
    mnt_ncr=mnt_ncr_chk(lat_in_grd,lat_in_nbr+1)
    if (.not.mnt_ncr) stop "ERROR: lat_in_grd not monotonically increasing in hyd_get()"
    mnt_ncr=mnt_ncr_chk(lon_in_grd,lon_in_nbr+1)
    if (.not.mnt_ncr) stop "ERROR: lon_in_grd not monotonically increasing in hyd_get()"
    
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
    
    ! Loop over months
    do time_out_idx=1,time_out_nbr
       srt_lon_lat_time=(/1,1,time_out_idx/) ! [idx] Dimension offset indices
       cnt_lon_lat_time=(/lon_in_nbr,lat_in_nbr,1/) ! [idx] Dimension sizes
       ! Get data
       rcd=nf90_wrp(nf90_get_var(nc_id,sfc_flw_id,sfc_flw_in, &
            start=srt_lon_lat_time,count=cnt_lon_lat_time),"get_var sfc_flw_in")
       ! Rebin optical depth from input grid to output grid
       call map_rbn(sfc_flw_in,mss_flg,mss_val_in, & ! I
            lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
            lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
            ovr_nbr_max,ovr_lat_idx,ovr_lon_idx,ovr_nbr,ovr_wgt, & ! I
            sfc_flw_out(1,1,time_out_idx)) ! O
       
       ! Rebin source strength from input grid to output grid
       call map_rbn(vwc_sfc_in,mss_flg,mss_val_in, & ! I
            lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
            lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
            ovr_nbr_max,ovr_lat_idx,ovr_lon_idx,ovr_nbr,ovr_wgt, & ! I
            vwc_sfc_out(1,1,time_out_idx)) ! O
       
    end do ! end loop over mth
    
    ! Close file
    rcd=nf90_wrp_close(nc_id,fl_in,"Ingested",sbr_nm=sbr_nm) ! [fnc] Close file
    
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
    
    if (allocated(sfc_flw_in)) deallocate(sfc_flw_in,stat=rcd)
    if(rcd /= 0) stop "deallocate() failed for sfc_flw_in"
    if (allocated(vwc_sfc_in)) deallocate(vwc_sfc_in,stat=rcd)
    if(rcd /= 0) stop "deallocate() failed for vwc_sfc_in"
    
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Exiting "//sbr_nm
    return 
  end subroutine hyd_get

  subroutine mdl_lsm_get( & ! [sbr] Retrieve and rebin volumetric water content
       fl_in, & ! I
       lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr,time_in_nbr, & ! I
       lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr,time_out_nbr, & ! I
       area_out,ovr_nbr_max, & ! I
       vwc_sfc_out) ! O 
    ! Purpose: Retrieve and rebin land surface data
    ! mdl_lsm_get() is called by asm_drv()
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use map_cst ! [mdl] Constants used in map routines
    use map_grd ! [mdl] Map grids and regridding
    use netcdf ! [mdl] netCDF interface
    use nf90_utl ! [mdl] netCDF utilities
    use sng_mdl,only:ftn_strlen,ftn_strstr ! [mdl] String manipulation
    use utl_mdl,only:mnt_chk ! [mdl] Utility functions (date_time_get,mnt_chk...)
    implicit none
    ! Parameters
    character(len=*),parameter::sbr_nm="mdl_lsm_get" ! [sng] Subroutine name
    ! Commons
    ! Input
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
    real,intent(out)::vwc_sfc_out(lon_out_nbr,lat_out_nbr,time_out_nbr) ! [m3 m-3] Volumetric water content
    ! Input/Output
    character(len=*),intent(inout)::fl_in ! [sng] Input file
    ! Local
    character(2)::mth_sng       ! [sng] Month in MM format
    integer idx_mth_sng       ! [idx] Location of month index in input filename
    integer lat_out_idx       ! [idx] Counting index for lat
    integer lon_out_idx       ! [idx] Counting index for lon
    integer mss_val_id        ! [id] Attribute ID
    integer nc_id             ! [id] File handle
    integer vwc_sfc_id        ! [id] Variable ID
    integer ovr_idx           ! [idx] Counting index
    integer ovr_lat_idx(lon_out_nbr,lat_out_nbr,ovr_nbr_max) ! [idx] Map into input grid of latitude indices of overlap cells
    integer ovr_lon_idx(lon_out_nbr,lat_out_nbr,ovr_nbr_max) ! [idx] Map into input grid of longitude indices of overlap cells
    integer ovr_nbr(lon_out_nbr,lat_out_nbr) ! [nbr] Number of input gridcells which overlap each output gridcell
    integer ovr_nbr_crr       ! [nbr] Current number of overlapping gridcells
    integer rcd               ! [enm] Return success code
    integer time_dmn_id       ! [id] Dimension ID for time
    integer time_nbr          ! [nbr] Number of records (months) in input file
    integer time_out_idx      ! [idx] Counting index for time
    integer vwc_dmn_nbr       ! [nbr] Number of dimensions in disk vwc_sfc field
    logical mnt               ! [flg] Monotonicity flag
    logical mss_flg           ! [flg] Check for mss_val on input data
    real mss_val              ! Missing value
    real vwc_sfc_in(lon_in_nbr,lat_in_nbr) ! [m3 m-3] Volumetric water content 
    real ovr_wgt(lon_out_nbr,lat_out_nbr,ovr_nbr_max) ! [frc] Weight of overlapping input gridcells onto each output gridcell
    
    ! Main Code
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Entering mdl_lsm_get()"
    
    ! Enough memory?
    if (time_out_nbr /= 12) stop "mdl_lsm_get(): time_out_nbr /= 12"
    
    ! Check for monotonicity
    mnt=mnt_chk(lon_out_grd,lon_out_nbr+1)
    if (.not.mnt) stop "lon_out_grd not monotonic in mdl_lsm_get()"
    mnt=mnt_chk(lat_out_grd,lat_out_nbr+1)
    if (.not.mnt) stop "lat_out_grd not monotonic in mdl_lsm_get()"
    ! Check for increasing monotonicity
    if (lon_out_grd(2) < lon_out_grd(1)) stop "lon_out_grd not increasing in mdl_lsm_get()"
    if (lat_out_grd(2) < lat_out_grd(1)) stop "lat_out_grd not increasing in mdl_lsm_get()"
    
    ! Get overlap locations and weights for mapping input grid to output grid
    call map_ovr_wgt_drv( &
         lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
         lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
         ovr_nbr_max,area_out, & ! I
         ovr_lat_idx,ovr_lon_idx,ovr_nbr,ovr_wgt) ! O
    
    ! Initialize scalars
    mss_flg=.false.
    mss_val=1.0e36
    idx_mth_sng=ftn_strstr(fl_in,"01.nc")
    
    ! Loop over months
    do time_out_idx=1,time_out_nbr
       write (mth_sng,"(i2.2)") time_out_idx
       fl_in(idx_mth_sng:idx_mth_sng+1)=mth_sng
       ! Read netCDF data
       rcd=nf90_noerr           ! nf90_noerr == 0
       rcd=nf90_wrp_open(fl_in,NF90_NOWRITE,nc_id,sbr_nm=sbr_nm)
       ! Get dimension IDs
       rcd=nf90_wrp_inq_dimid(nc_id,"time",time_dmn_id)
       ! Get dimension sizes
       rcd=rcd+nf90_inquire_dimension(nc_id,time_dmn_id,len=time_nbr)
       if (time_nbr /= time_in_nbr) stop "time_nbr /= time_in_nbr in mdl_lsm_get()"
       if (time_nbr > 1) then
          write (6,"(2a,f12.6)") prg_nm(1:ftn_strlen(prg_nm)), & 
               ": ERROR mdl_lsm_get() reports time_nbr = ",time_nbr
          stop "time_nbr > 1"
       endif                  ! endif err
       ! Get variable IDs
       rcd=nf90_wrp_inq_varid(nc_id,"H2OSOI",vwc_sfc_id)
       ! Get number of dimensions
       rcd=rcd+nf90_inquire_variable(nc_id,vwc_sfc_id,ndims=vwc_dmn_nbr)
       if (vwc_dmn_nbr /= 3) then
          write (6,"(2a,i1,a)") prg_nm(1:ftn_strlen(prg_nm)), & 
               ": ERROR H2OSOI has ",vwc_dmn_nbr," dimensions"
          write (6,"(2a,i1,a)") prg_nm(1:ftn_strlen(prg_nm)), & 
               ": HINT Make sure H2OSOI has no vertical (soil depth) dimension"
          stop
       endif                  ! endif err
       ! Get data
       rcd=nf90_wrp(nf90_get_var(nc_id,vwc_sfc_id,vwc_sfc_in),"get_var vwc_sfc_in")
       rcd=nf90_wrp_close(nc_id,fl_in,'',sbr_nm=sbr_nm)
       if (dbg_lvl >= dbg_fl) write (6,"(a8,1x,a)") "Ingested",fl_in(1:ftn_strlen(fl_in))
       
       ! Rebin volumetric soil water from input grid to output grid
       call map_rbn(vwc_sfc_in,mss_flg,mss_val, & ! I
            lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
            lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
            ovr_nbr_max,ovr_lat_idx,ovr_lon_idx,ovr_nbr,ovr_wgt, & ! I
            vwc_sfc_out(1,1,time_out_idx)) ! O
       
    end do ! end loop over mth
    
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Exiting mdl_lsm_get()"
    return 
  end subroutine mdl_lsm_get
  
  subroutine lsm_get_sfc_flw( & ! [sbr] Retrieve and rebin surface flow data QOVER
       fl_in, & ! I
       lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr,time_in_nbr, & ! I
       lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr,time_out_nbr, & ! I
       area_out,ovr_nbr_max, & ! I
       sfc_flw_out) ! O 
    ! Purpose: Retrieve and rebin surface flow data QOVER
    ! lsm_get_sfc_flw() is called by asm_drv() just after lsm_get()
    ! Routine based on lsm_get()
    use netcdf ! [mdl] netCDF interface
    use nf90_utl ! [mdl] netCDF utilities
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use map_cst ! [mdl] Constants used in map routines
    use map_grd ! [mdl] Map grids and regridding
    use sng_mdl,only:ftn_strlen,ftn_strstr ! [mdl] String manipulation
    use utl_mdl,only:mnt_chk ! [mdl] Utility functions (date_time_get,mnt_chk...)
    implicit none
    ! Parameters
    character(len=*),parameter::sbr_nm="lsm_get_sfc_flw" ! [sng] Subroutine name
    ! Commons
    ! Input
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
    real,intent(out)::sfc_flw_out(lon_out_nbr,lat_out_nbr) ! [mm s-1] Surface runoff (QOVER)
    ! Input/Output
    character(80),intent(inout)::fl_in ! [sng] Input file
    ! Local
    character(2)::mth_sng       ! [sng] Month in MM format
    integer idx_mth_sng       ! [idx] Location of month index in input filename
    integer lat_out_idx       ! [idx] Counting index for lat
    integer lon_out_idx       ! [idx] Counting index for lon
    integer mss_val_id        ! [id] Attribute ID
    integer nc_id             ! [id] File handle
    integer sfc_flw_id        ! [id] Variable ID
    integer ovr_idx           ! [idx] Counting index
    integer ovr_lat_idx(lon_out_nbr,lat_out_nbr,ovr_nbr_max) ! [idx] Map into input grid of latitude indices of overlap cells
    integer ovr_lon_idx(lon_out_nbr,lat_out_nbr,ovr_nbr_max) ! [idx] Map into input grid of longitude indices of overlap cells
    integer ovr_nbr(lon_out_nbr,lat_out_nbr) ! [nbr] Number of input gridcells which overlap each output gridcell
    integer ovr_nbr_crr       ! [nbr] Current number of overlapping gridcells
    integer rcd               ! [enm] Return success code
    integer time_dmn_id       ! [id] Dimension ID for time
    integer time_nbr          ! [nbr] Number of records (months) in input file
    integer time_out_idx      ! [idx] Counting index for time
    integer sfc_flw_dmn_nbr   ! [nbr] Number of dimensions in disk sfc_flw field
    logical mnt               ! [flg] Monotonicity flag
    logical mss_flg           ! [flg] Check for mss_val on input data
    real mss_val              ! Missing value
    real sfc_flw_in(lon_in_nbr,lat_in_nbr) ! [mm s-1] Surface runoff
    real ovr_wgt(lon_out_nbr,lat_out_nbr,ovr_nbr_max) ! [frc] Weight of overlapping input gridcells onto each output gridcell
    
    ! Main Code
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Entering lsm_get_sfc_flw()"
    
    ! Enough memory?
    if (time_out_nbr /= 12) stop "lsm_get_sfc_flw(): time_out_nbr /= 12"
    
    ! Check for monotonicity
    mnt=mnt_chk(lon_out_grd,lon_out_nbr+1)
    if (.not.mnt) stop "lon_out_grd not monotonic in lsm_get_sfc_flw()"
    mnt=mnt_chk(lat_out_grd,lat_out_nbr+1)
    if (.not.mnt) stop "lat_out_grd not monotonic in lsm_get_sfc_flw()"
    ! Check for increasing monotonicity
    if (lon_out_grd(2) < lon_out_grd(1)) stop "lon_out_grd not increasing in lsm_get_sfc_flw()"
    if (lat_out_grd(2) < lat_out_grd(1)) stop "lat_out_grd not increasing in lsm_get_sfc_flw()"
    
    ! Get overlap locations and weights for mapping input grid to output grid
    call map_ovr_wgt_drv( &
         lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
         lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
         ovr_nbr_max,area_out, & ! I
         ovr_lat_idx,ovr_lon_idx,ovr_nbr,ovr_wgt) ! O
    
    ! Initialize scalars
    mss_flg=.true.
    mss_val=1.0e36
    
    ! Read netCDF data
    rcd=nf90_noerr           ! nf90_noerr == 0
    rcd=nf90_wrp_open(fl_in,NF90_NOWRITE,nc_id,sbr_nm=sbr_nm)
    ! Get dimension IDs
    rcd=nf90_wrp_inq_dimid(nc_id,"time",time_dmn_id)
    ! Get dimension sizes
    rcd=rcd+nf90_inquire_dimension(nc_id,time_dmn_id,len=time_nbr)
    if (time_nbr /= time_in_nbr) stop "time_nbr /= time_in_nbr in lsm_get_sfc_flw()"
    if (time_nbr > 1) then
       write (6,"(2a,f12.6)") prg_nm(1:ftn_strlen(prg_nm)), & 
            ": ERROR lsm_get_sfc_flw() reports time_nbr = ",time_nbr
       stop "time_nbr > 1"
    endif                  ! endif err
    ! Get variable IDs
    rcd=nf90_wrp_inq_varid(nc_id,"QOVER",sfc_flw_id)
    ! Get number of dimensions
    rcd=rcd+nf90_inquire_variable(nc_id,sfc_flw_id,ndims=sfc_flw_dmn_nbr)
    if (sfc_flw_dmn_nbr /= 3) then
       write (6,"(2a,i1,a)") prg_nm(1:ftn_strlen(prg_nm)), & 
            ": ERROR QOVER has ",sfc_flw_dmn_nbr," dimensions"
       write (6,"(2a,i1,a)") prg_nm(1:ftn_strlen(prg_nm)), & 
            ": HINT Make sure QOVER has no vertical (soil depth) dimension"
       stop
    endif                  ! endif err
    ! Get data
    rcd=nf90_wrp(nf90_get_var(nc_id,sfc_flw_id,sfc_flw_in),"get_var sfc_flw_in")
    rcd=nf90_wrp_close(nc_id,fl_in,'',sbr_nm=sbr_nm)
    if (dbg_lvl >= dbg_fl) write (6,"(a8,1x,a)") "Ingested",fl_in(1:ftn_strlen(fl_in))

    if (dbg_lvl >= dbg_fl) then
       write (6,"(2a,e14.7)") prg_nm(1:ftn_strlen(prg_nm)), &
            "DEBUG "//sbr_nm//" reports max(QOVER) = ",maxval(sfc_flw_in)
       write (6,"(2a,e14.7)") prg_nm(1:ftn_strlen(prg_nm)), &
            "DEBUG "//sbr_nm//" reports min(QOVER) = ",minval(sfc_flw_in)
    end if ! endif dbg
    
    ! Set negative QOVER to zero, and missing values to zero
    ! Do we need any of this???
    !where (sfc_flw_in > 1.0e29)
       !sfc_flw_in = 0.0
    !end where
    !where (sfc_flw_in < 0.0)
       !sfc_flw_in = 0.0
    !end where
    
    ! Rebin from input grid to output grid
    call map_rbn(sfc_flw_in,mss_flg,mss_val, & ! I
         lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
         lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
         ovr_nbr_max,ovr_lat_idx,ovr_lon_idx,ovr_nbr,ovr_wgt, & ! I
         sfc_flw_out) ! O
    
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Exiting lsm_get_sfc_flw()"
    return 
  end subroutine lsm_get_sfc_flw
  
end module asm_mdl ! [mdl] Assimilate model fields
