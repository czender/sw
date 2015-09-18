! $Id$ -*-f90-*- 

! Purpose: Control routines for Boundary DataSet processor bds

! Usage:
!use bds_drv ! [mdl] BDS driver, processor, output routines

module bds_drv ! [mdl] BDS driver, processor, output routines
  implicit none
  public::bds_prc ! [sbr] Boundary DataSet processor
  private::lnd_bnd_cnv_rbn ! [sbr] Convert, rebin external datasets
  private::bds_out ! [sbr] Write BDS output to netCDF file
 
contains 

  subroutine bds_prc( & ! [sbr] Boundary DataSet processor
       bln_nbr, & ! [nbr] Number of tri-modal soil blends
       lat_nbr, & ! [nbr] Number of latitudes
       lon_nbr, & ! [nbr] Number of longitudes
       map_typ_out, & ! [enm] Output map grid type
       sgs_nbr, & ! [nbr] Number of sub-gridscale patches per gridcell
       time_out_nbr, & ! [nbr] Number of times
       yr_nbr, & ! [nbr] Number of years
       yr_srt & ! [yr] Starting year in YYYY format
       )
    ! Purpose: Dimension all output variables then call the routines reposible for input, processing, and output
    ! bds_prc() is called by main()
    use asm_mdl,only:asm_drv ! [mdl] Assimilate model fields
    use bds_ctl,only:bsn_fct_hrz, &
         fl_chm,fl_hgt,fl_lai,fl_lak,fl_rgn_msk,fl_rgn_xcl,fl_pft,fl_rfl_sfc, &
         fl_sfc,fl_soi_txt,fl_tpr,fl_wtl, &
         flg_rdb_only,flg_cst_tpr,flg_rgn_msk,flg_rgn_xcl,src_thr! [mdl] Control variables drc_in,drc_out,hgt_dlt_msl
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use geo_mdl,only:dst_src_flw ! [mdl] Geomorphology and flow accumulation
    use map_cst ! [mdl] Constants used in map routines
    use map_grd ! [mdl] Map grids and regridding
    use sng_mdl,only:ftn_strlen ! [mdl] String manipulation
    use sfc_typ_mdl,only:sfc2dps,sfc2vgt ! [mdl] Surface type and derived properties
    use soi_txt,only:soi_hyd_prc ! [mdl] Soil texture and hydrology
    use src_prc,only:aod_nmd_get,cst_tpr,rdb_fct_tpg_GCT01,rgn_msk,rgn_xcl,src_id ! [mdl] Source processing, identification
    implicit none
    ! Parameters
    character(len=*),parameter::sbr_nm="bds_prc" ! [sng] Subroutine name
    integer,parameter::time_in_nbr_max=1 ! Input files consist of 1 time units (months)
    integer,parameter::fl_tpr_unit=73 ! [idx] Unit for reading taper rectangles
    integer,parameter::fl_rgn_msk_unit=73 ! [idx] Unit for reading mask rectangle
    integer,parameter::fl_rgn_xcl_unit=73 ! [idx] Unit for reading zero rectangle
    real,parameter::itr_fct_max=100.0 ! [frc] Arbitrary iteration cutoff
    real,parameter::itr_fct_min=0.001 ! [frc] Arbitrary iteration cutoff
    ! Commons
    ! Input
    integer,intent(in)::bln_nbr ! [nbr] Number of tri-modal soil blends
    integer,intent(in)::lat_nbr ! [nbr] Number of latitudes
    integer,intent(in)::lon_nbr ! [nbr] Number of longitudes
    integer,intent(in)::map_typ_out ! [enm] Output map grid type
    integer,intent(in)::sgs_nbr ! [nbr] Number of sub-gridscale patches per gridcell
    integer,intent(in)::time_out_nbr ! [nbr] Number of times
    integer,intent(in)::yr_nbr ! [nbr] Number of years
    integer,intent(in)::yr_srt ! [yr] Starting year in YYYY format
    ! Input/Output
    ! Output
    ! Local
    integer bsn_enm(lon_nbr,lat_nbr) ! [enm] Basin ID
    integer lat_idx           ! [idx] Counting index
    integer lnd_msk(lon_nbr,lat_nbr) ! [msk] Land mask
    integer lon_idx           ! [idx] Counting index
    integer mth_max(lon_nbr,lat_nbr) ! [idx] Month of seasonal maximum
    integer msk_rgn(lon_nbr,lat_nbr) ! [flg] Mask region
    integer pft_clm(lon_nbr,lat_nbr,sgs_nbr) ! [enm] Plant functional type
    integer sfc_typ(lon_nbr,lat_nbr) ! [enm] Surface type code
    integer src_flg(lon_nbr,lat_nbr) ! [flg] Source flag
    integer bsn_nbr ! [nbr] Number of basins
    integer slat_nbr ! [nbr] Number of staggered latitudes for FV grid
    integer slon_nbr ! [nbr] Number of staggered longitudes for FV grid
    integer time_out_idx      ! [idx] Counting index
    real aod_nmd_frc(lon_nbr,lat_nbr,time_out_nbr) ! [frc] TOMS-observed aerosol index not due to dust
    real area(lon_nbr,lat_nbr) ! [m2] Area of gridcells 
    real sfc_acm_fct(lon_nbr,lat_nbr) ! [frc] Surface area accumulation factor
    real rfl_sfc_lnd_nrm_mds_lnr(lon_nbr,lat_nbr) ! [frc] Erodibility factor from MODIS, linear
    real rfl_sfc_lnd_nrm_mds_sqr(lon_nbr,lat_nbr) ! [frc] Erodibility factor from MODIS, squared
    real csn_cyc(lon_nbr,lat_nbr,time_out_nbr) ! [frc] Seasonal cycle
    real csn_cyc_mdl(lon_nbr,lat_nbr,time_out_nbr) ! [frc] Seasonal cycle, CCM
    real flw_acm_fct(lon_nbr,lat_nbr) ! [frc] Flow accumulation factor
    real fsh_fct(lon_nbr,lat_nbr) ! [frc] Efficiency factor
    real grd_mrd_lng(lat_nbr) ! [m] Grid meridional length
    real grd_rds_rth(lat_nbr) ! [m] Distance from axis of rotation at local latitude
    real grd_znl_lng(lat_nbr) ! [m] Grid zonal length
    real gw(lat_nbr) ! Gaussian weight
    real hgt_sfc(lon_nbr,lat_nbr) ! [m] Surface height 
    real hgt_sfc_gpm(lon_nbr,lat_nbr) ! [gpm] Surface geopotential height 
    real hgt_sfc_std_dvn(lon_nbr,lat_nbr) ! [m] Standard deviation of surface height 
    real hgt_zpd_lsm(lon_nbr,lat_nbr) ! [m] Zero plane displacement height
    real lai_clm(lon_nbr,lat_nbr,sgs_nbr,time_out_nbr) ! [m2 m-2] Leaf area index, one-sided
    real lai_lsm(lon_nbr,lat_nbr,time_out_nbr) ! [m2 m-2] Leaf area index
    real lai_ttl_clm(lon_nbr,lat_nbr,time_out_nbr) ! [m2 m-2] Total leaf area index, one-sided
    real lak_frc(lon_nbr,lat_nbr) ! [frc] Lake fraction of gridcell
    real lat(lat_nbr) ! Coordinate variable
    real lat_grd(lat_nbr+1) ! Latitude grid
    real lat_sz(lat_nbr) ! Latitudinal size of bin
    real lnd_frc(lon_nbr,lat_nbr) ! [frc] Land fraction
    real lnd_frc_clm(lon_nbr,lat_nbr) ! [frc] Fraction of land (not ocean)
    real lnd_frc_dry(lon_nbr,lat_nbr) ! [frc] Dry fraction of gridcell
    real lon(lon_nbr) ! Coordinate variable
    real lon_grd(lon_nbr+1) ! Longitude grid
    real lon_sz(lon_nbr) ! Longitudinal size of bin
    real mbl_bsn_fct(lon_nbr,lat_nbr) ! [frc] Mobilization enhancement due to basin characteristics
    real mss_frc_Al(lon_nbr,lat_nbr) ! [frc] Exchangeable aluminum
    real mss_frc_C_org(lon_nbr,lat_nbr) ! [frc] Organic carbon
    real mss_frc_CaCO3(lon_nbr,lat_nbr) ! [frc] Calcium carbonate
    real mss_frc_K(lon_nbr,lat_nbr) ! [frc] Exchangeable potassium
    real mss_frc_N(lon_nbr,lat_nbr) ! [frc] Nitrogen
    real mss_frc_Na(lon_nbr,lat_nbr) ! [frc] Exchangeable sodium
    real mss_frc_P2O5(lon_nbr,lat_nbr) ! [frc] Extractable phosphorus
    real sfc_frc_bln(lon_nbr,lat_nbr,bln_nbr) ! [m2 m-2] Surface area fraction of soil blend
    real mss_frc_cly(lon_nbr,lat_nbr) ! [frc] Soil texture clay
    real mss_frc_slt(lon_nbr,lat_nbr) ! [frc] Soil texture silt
    real mss_frc_snd(lon_nbr,lat_nbr) ! [frc] Soil texture sand
    real odxc_tms(lon_nbr,lat_nbr,time_out_nbr) ! [frc] Extinction optical depth proxy from TOMS
    real oro(lon_nbr,lat_nbr) ! [frc] Orography 
    real pft_frc_clm(lon_nbr,lat_nbr,sgs_nbr) ! [enm] Fraction covered by plant functional type
    real rgh_mmn(lon_nbr,lat_nbr) ! [m] Roughness length humidity
    real smp_sat(lon_nbr,lat_nbr) ! [mm H2O] Saturated soil matric potential (sand-dependent)
    real smp_sfc(lon_nbr,lat_nbr,time_out_nbr) ! [mm H2O] Soil matric potential
    real src_frq(lon_nbr,lat_nbr,time_out_nbr) ! [frc] Source frequency
    real src_odxc(lon_nbr,lat_nbr,time_out_nbr) ! [frc] Extinction optical depth, source regions only
    real src_str(lon_nbr,lat_nbr,time_out_nbr) ! [frc] Source strength
    real src_str_mdl(lon_nbr,lat_nbr,time_out_nbr) ! [frc] Source strength for next CCM run
    real src_str_old(lon_nbr,lat_nbr,time_out_nbr) ! [frc] Source strength used by last CCM run
    real time_out(time_out_nbr) ! Coordinate variable
    real vai_clm(lon_nbr,lat_nbr,sgs_nbr,time_out_nbr) ! [m2 m-2] Vegetation area index, one-sided
    real vai_lsm(lon_nbr,lat_nbr,time_out_nbr) ! [m2 m-2] Vegetation area index
    real vai_ttl_clm(lon_nbr,lat_nbr,time_out_nbr) ! [m2 m-2] Total vegetation area index, one-sided
    real vwc_dry(lon_nbr,lat_nbr) ! [m3 m-3] Dry volumetric water content (no E-T)
    real vwc_opt(lon_nbr,lat_nbr) ! [m3 m-3] E-T optimal volumetric water content 
    real vwc_rel(lon_nbr,lat_nbr,time_out_nbr) ! [frc] Water content relative to saturation
    real vwc_sat(lon_nbr,lat_nbr) ! [m3 m-3] Saturated volumetric water content (sand-dependent)
    real vwc_sfc(lon_nbr,lat_nbr,time_out_nbr) ! [m3 m-3] Volumetric water content
    real wtl_frc(lon_nbr,lat_nbr) ! [frc] Wetland fraction of gridcell
    real(selected_real_kind(p=12))::pi ! [frc] 3
    ! Local
    integer rcd ! [enm] Return success code
    integer rct_idx ! [idx] Counting index for rectangles
    integer rct_nbr ! [nbr] Number of rectangles
    integer time_max_idx ! [idx] Time unit (month) of maximum emission
    integer lat_in_nbr ! [nbr] Number of latitudes
    integer lon_in_nbr ! [nbr] Number of longitudes
    integer lat_out_nbr ! [nbr] Number of latitudes
    integer lon_out_nbr ! [nbr] Number of longitudes

    logical itr_1st ! [flg] First iteration
    real itr_fct ! [frc] Ratio of observed to predicted optical depth
    real lat_max ! [dgr] Maximum latitude in rectangle
    real lat_min ! [dgr] Minimum latitude in rectangle
    real lon_max ! [dgr] Maximum longitude in rectangle
    real lon_min ! [dgr] Minimum longitude in rectangle
    real mss_val ! [frc] Missing value
    real mth_day_ctr(12) ! [day] Day of year at mid-month (i.e., mid-January=16)
    real odxc_mdl(lon_nbr,lat_nbr,time_out_nbr) ! [frc] Extinction optical depth from model
    real sfc_flw(lon_nbr,lat_nbr) 
    real bsn_sz(lon_nbr*lat_nbr) ! [m2] Basin area
    data mth_day_ctr/ 16, 45, 75,105,136,166,197,228,258,289,319,350/

    ! Main code
    if (dbg_lvl >= dbg_sbr) write(6,"(a)") "Entering bds_prc()"
    
    ! Initialize default values
    ! 19990128: Setting mss_val=1.0e36 is a kludge to obviate passing actual mss_val up from toms_get() 
    mss_val=1.0e36
    ! This mss_val must agree with odxc_tms mss_val
    if (mss_val /= 1.0e36) write (6,"(2a,f12.6)") prg_nm(1:ftn_strlen(prg_nm)), & 
         ": WARNING bds_prc() mss_val = ",mss_val
    
    ! Initialize output arrays
    msk_rgn(:,:)=0 ! [flg] Mask region
    mth_max(:,:)=0 ! [idx] Month of seasonal maximum
    csn_cyc(:,:,:)=0.0 ! [frc] Seasonal cycle
    src_odxc(:,:,:)=0.0 ! [frc] Extinction optical depth, source regions only
    time_out=mth_day_ctr
    
    ! Create output map grid
    call map_grd_mk(lat_nbr,lon_nbr,map_typ_out, & ! I
         lat_grd,lon_grd) ! O
    call map_lat_wgt_mk(lat_nbr,map_typ_out, & ! I
         gw) ! O
    do lon_idx=1,lon_nbr
       lon(lon_idx)=0.5*(lon_grd(lon_idx)+lon_grd(lon_idx+1))
       lon_sz(lon_idx)=lon_grd(lon_idx+1)-lon_grd(lon_idx)
    end do ! end loop over lon
    do lat_idx=1,lat_nbr
       lat(lat_idx)=0.5*(lat_grd(lat_idx)+lat_grd(lat_idx+1))
       lat_sz(lat_idx)=lat_grd(lat_idx+1)-lat_grd(lat_idx)
    end do ! end loop over lat
    
    ! Distance along each side of gridcell
    pi=4.0*atan(1.0d0) ! [frc] 3
    do lat_idx=1,lat_nbr
       grd_rds_rth(lat_idx)=earth_rds*cos(pi*lat(lat_idx)/180.0) ! [m] Distance from axis of rotation at local latitude
       ! Assumes regular longitude array
       ! Mean zonal distance across gridcell (along center latitude)
       grd_znl_lng(lat_idx)=grd_rds_rth(lat_idx)*lon_sz(1)*pi/180.0 ! [m] Grid zonal length
       ! Mean meridional distance across gridcell (along center longitude)
       grd_mrd_lng(lat_idx)=grd_rds_rth(lat_idx)*lat_sz(1)*pi/180.0 ! [m] Grid meridional length
    end do ! end loop over lat
    
    ! Compute gridcell area 
    call map_area_get(lat_grd,lat_nbr,lon_grd,lon_nbr, & ! I
         area) ! O
    
    ! Generate land surface boundary data from disparate external datasets
    call lnd_bnd_cnv_rbn(lat_grd,lat_nbr,lon_grd,lon_nbr,sgs_nbr,time_out_nbr,area, & ! I
         lak_frc,lnd_frc,lnd_frc_dry,lnd_msk,sfc_typ, & ! O
         mss_frc_Al, & ! O
         mss_frc_C_org, & ! O
         mss_frc_CaCO3, & ! O
         mss_frc_K, & ! O
         mss_frc_N, & ! O
         mss_frc_Na, & ! O
         mss_frc_P2O5, & ! O
         mss_frc_cly,mss_frc_slt,mss_frc_snd, & ! O
         hgt_sfc_gpm,hgt_sfc,hgt_sfc_std_dvn,oro, & ! O 
         mbl_bsn_fct, & ! O
         bsn_enm, & ! O [enm] Basin ID
         flw_acm_fct, & ! O
         sfc_acm_fct, & ! O
         rfl_sfc_lnd_nrm_mds_lnr, & ! O [frc] Erodibility factor from MODIS, linear
         rfl_sfc_lnd_nrm_mds_sqr, & ! O [frc] Erodibility factor from MODIS, squared
         pft_clm, & ! O [enm] Plant functional type
         pft_frc_clm, & ! O [enm] Fraction covered by plant functional type
         lnd_frc_clm, & ! O [frc] Fraction of land (not ocean)
         lai_clm, & ! O [m2 m-2] Leaf area index, one-sided
         lai_ttl_clm, & ! O [m2 m-2] Total leaf area index, one-sided
         vai_clm, & ! O [m2 m-2] Vegetation area index, one-sided
         vai_ttl_clm, & ! O [m2 m-2] Total vegetation area index, one-sided
         wtl_frc) ! O
    
    if(.not.bsn_fct_hrz) then
       ! Compute source efficiency factor from output topography grid
       write (6,"(2a)") prg_nm(1:ftn_strlen(prg_nm)),": Computing basin factor at model resolution"
       call rdb_fct_tpg_GCT01(lat_grd,lat_nbr,lon_grd,lon_nbr, & ! I
            hgt_sfc,oro, & ! I
            mbl_bsn_fct) ! O
    endif
    
    ! Assimilate modeled and observed seasonal cycle data onto output grid
    call asm_drv(lat_grd,lat_nbr,lon_grd,lon_nbr,time_out_nbr, & ! I
         yr_nbr,yr_srt, & ! I
         area,lnd_frc,lnd_msk, & ! I
         fsh_fct,odxc_tms,src_flg,src_frq,src_str, & ! O 
         odxc_mdl,src_str_old,vwc_sfc,sfc_flw) ! O 
    
    if(.not.bsn_fct_hrz) then
       ! Compute source efficiency factor from output topography grid
       write (6,"(2a)") prg_nm(1:ftn_strlen(prg_nm)),": Computing flow accumulation factors at model resolution"
       call dst_src_flw( & ! [sbr] Compute area and flow accumulation factors
            lat_nbr,lon_nbr, & ! I
            hgt_sfc,oro, &
            area,sfc_flw, & ! I
            bsn_enm,flw_acm_fct,sfc_acm_fct,bsn_sz) ! O
    endif ! endif bsn_fct_hrz
    ! Taper erodibility factors along coasts in rectangles
    if (flg_cst_tpr) then
       open (fl_tpr_unit,file=fl_tpr,status="old",iostat=rcd)
       if (rcd /= 0) then
          write (6,"(3a)") prg_nm(1:ftn_strlen(prg_nm)), & 
               ": ERROR "//sbr_nm//"() reports coastal tapering requested but unable to open ", & 
               fl_tpr(1:ftn_strlen(fl_tpr)) 
          write (6,"(6a)") prg_nm(1:ftn_strlen(prg_nm)), & 
               ": HINT ",fl_tpr(1:ftn_strlen(fl_tpr))," may be in ${HOME}/map. Run ", &
               prg_nm(1:ftn_strlen(prg_nm))," from there, or invoke with --flg_cst_tpr=F"
          stop
       endif ! endif error
       rct_idx=1 ! [idx] Counting index for rectangles
       do ! Begin loop over rectangles
          ! Input list is ordered [S,N,W,E]
          read (fl_tpr_unit,*) lat_min,lat_max,lon_min,lon_max
          ! Rectangle list ends with line of four zeroes
          if ((abs(lat_min)+abs(lat_max)+abs(lon_min)+abs(lon_max)) == 0.0) exit
          if (dbg_lvl > dbg_sbr) then
             write(6,"(2a,4(f6.2,a,1x))") prg_nm(1:ftn_strlen(prg_nm)), &
                  ": Tapering coastal region [S,N,W,E] ",lat_min,",",lat_max,",",lon_min,",",lon_max,"."
          endif ! endif dbg
          call cst_tpr(lat,lat_nbr,lon,lon_nbr, & ! I
               lat_min,lat_max,lon_min,lon_max, & ! I
               oro,flw_acm_fct)
          call cst_tpr(lat,lat_nbr,lon,lon_nbr, & ! I
               lat_min,lat_max,lon_min,lon_max, & ! I
               oro,sfc_acm_fct)
          call cst_tpr(lat,lat_nbr,lon,lon_nbr, & ! I
               lat_min,lat_max,lon_min,lon_max, & ! I
               oro,mbl_bsn_fct)
          call cst_tpr(lat,lat_nbr,lon,lon_nbr, & ! I
               lat_min,lat_max,lon_min,lon_max, & ! I
               oro,rfl_sfc_lnd_nrm_mds_lnr)
          call cst_tpr(lat,lat_nbr,lon,lon_nbr, & ! I
               lat_min,lat_max,lon_min,lon_max, & ! I
               oro,rfl_sfc_lnd_nrm_mds_sqr)
          rct_idx=rct_idx+1 ! [idx] Counting index for rectangles
       end do ! end loop over coastal regions to taper
       close (fl_tpr_unit)
       if (dbg_lvl > dbg_sbr) then
          write (6,"(a,i2,a)") "Applied cst_tpr() in ",rct_idx-1," regions"
       endif ! endif dbg
       write (6,"(a,1x,a)") "Read coastal regions to taper as dust sources from",fl_tpr(1:ftn_strlen(fl_tpr))
    end if ! endif flg_cst_tpr
 
    ! Zero erodibility factors in regions
    if (flg_rgn_xcl) then
       open (fl_rgn_xcl_unit,file=fl_rgn_xcl,status="old",iostat=rcd)
       if (rcd /= 0) then
          write (6,"(3a)") prg_nm(1:ftn_strlen(prg_nm)), & 
               ": ERROR "//sbr_nm//"() reports coastal tapering requested but unable to open ", & 
               fl_rgn_xcl(1:ftn_strlen(fl_rgn_xcl)) 
          write (6,"(6a)") prg_nm(1:ftn_strlen(prg_nm)), & 
               ": HINT ",fl_rgn_xcl(1:ftn_strlen(fl_rgn_xcl))," may be in ${HOME}/map. Run ", &
               prg_nm(1:ftn_strlen(prg_nm))," from there, or invoke with --flg_rgn_xcl=F"
          stop
       endif ! endif error
       rct_idx=1 ! [idx] Counting index for rectangles
       do ! Begin loop over rectangles
          ! Input list is ordered [S,N,W,E]
          read (fl_rgn_xcl_unit,*) lat_min,lat_max,lon_min,lon_max
          if ((abs(lat_min)+abs(lat_max)+abs(lon_min)+abs(lon_max)) == 0.0) exit
          if (dbg_lvl > dbg_sbr) then
             write(6,"(2a,4(f6.2,a,1x))") prg_nm(1:ftn_strlen(prg_nm)), &
                  ": Zeroing region [S,N,W,E] ",lat_min,",",lat_max,",",lon_min,",",lon_max,"."
          endif ! endif dbg
          call rgn_xcl(lat,lat_nbr,lon,lon_nbr, & ! I
               lat_min,lat_max,lon_min,lon_max, & ! I
               flw_acm_fct) ! I/O
          call rgn_xcl(lat,lat_nbr,lon,lon_nbr, & ! I
               lat_min,lat_max,lon_min,lon_max, & ! I
               sfc_acm_fct) ! I/O
          call rgn_xcl(lat,lat_nbr,lon,lon_nbr, & ! I
               lat_min,lat_max,lon_min,lon_max, & ! I
               mbl_bsn_fct) ! I/O
          call rgn_xcl(lat,lat_nbr,lon,lon_nbr, & ! I
               lat_min,lat_max,lon_min,lon_max, & ! I
               rfl_sfc_lnd_nrm_mds_lnr) ! I/O
          call rgn_xcl(lat,lat_nbr,lon,lon_nbr, & ! I
               lat_min,lat_max,lon_min,lon_max, & ! I
               rfl_sfc_lnd_nrm_mds_sqr) ! I/O
          rct_idx=rct_idx+1 ! [idx] Counting index for rectangles
       end do ! end loop over locations to zero
       close (fl_rgn_xcl_unit)
       if (dbg_lvl > dbg_sbr) then
          write (6,"(a,i2,a)") "Applied rgn_xcl() in ",rct_idx-1," regions"
       endif ! endif dbg
       write (6,"(a,1x,a)") "Read regions to zero as dust sources from",fl_rgn_xcl(1:ftn_strlen(fl_rgn_xcl))
    end if ! endif flg_rgn_xcl
  
    ! Mask regions
    if (flg_rgn_msk) then
       open (fl_rgn_msk_unit,file=fl_rgn_msk,status="old",iostat=rcd)
       if (rcd /= 0) then
          write (6,"(3a)") prg_nm(1:ftn_strlen(prg_nm)), & 
               ": ERROR "//sbr_nm//"() reports masking region requested but unable to open ", & 
               fl_rgn_msk(1:ftn_strlen(fl_rgn_msk)) 
          write (6,"(6a)") prg_nm(1:ftn_strlen(prg_nm)), & 
               ": HINT ",fl_rgn_msk(1:ftn_strlen(fl_rgn_msk))," may be in ${HOME}/map. Run ", &
               prg_nm(1:ftn_strlen(prg_nm))," from there, or invoke with --flg_rgn_msk=F"
          stop
       endif ! endif error
       rct_idx=1 ! [idx] Counting index for rectangles
       do ! Begin loop over rectangles
          ! Input list is ordered [S,N,W,E]
          read (fl_rgn_msk_unit,*) lat_min,lat_max,lon_min,lon_max
          if ((abs(lat_min)+abs(lat_max)+abs(lon_min)+abs(lon_max)) == 0.0) exit
          if (dbg_lvl > dbg_sbr) then
             write(6,"(2a,4(f6.2,a,1x))") prg_nm(1:ftn_strlen(prg_nm)), &
                  ": Masking true region [S,N,W,E] ",lat_min,",",lat_max,",",lon_min,",",lon_max,"."
          endif ! endif dbg
          call rgn_msk(lat,lat_nbr,lon,lon_nbr, & ! I
               lat_min,lat_max,lon_min,lon_max, & ! I
               msk_rgn) ! I/O
          rct_idx=rct_idx+1 ! [idx] Counting index for rectangles
       end do ! end loop over locations to zero
       close (fl_rgn_msk_unit)
       if (dbg_lvl > dbg_sbr) then
          write (6,"(a,i2,a)") "Applied rgn_msk() in ",rct_idx-1," regions"
       endif ! endif dbg
       write (6,"(a,1x,a)") "Read regions to mask from",fl_rgn_msk(1:ftn_strlen(fl_rgn_msk))
    end if ! endif flg_rgn_msk
  
    ! Process prescribed LSM vegetation data into seasonal cycle surface fields
    call sfc2vgt(lat,lat_nbr,lon,lon_nbr,time_out_nbr, & ! I
         sfc_typ, & ! I
         lai_lsm,vai_lsm) ! O
    
    ! Process land surface structure data into time-invariant surface fields
    call sfc2dps(lat,lat_nbr,lon,lon_nbr, & ! I
         sfc_typ, & ! I
         rgh_mmn,hgt_zpd_lsm) ! O
    
    ! Process soil hydrology data into seasonal cycle surface fields
    call soi_hyd_prc(lat,lat_nbr,lon,lon_nbr,time_out_nbr, & ! I
         mss_frc_cly,mss_frc_slt,mss_frc_snd, & ! I
         vwc_sfc, & ! I
         smp_sat,smp_sfc,vwc_dry,vwc_opt,vwc_rel,vwc_sat) ! O
    
    ! Retrieve non-mineral dust fractions
    call aod_nmd_get(lat_grd,lat_nbr,lon_grd,lon_nbr,time_out_nbr, & ! I
         aod_nmd_frc) ! O
    
    ! Normalize results by period of analyses
    do lon_idx=1,lon_nbr
       do lat_idx=1,lat_nbr
          time_max_idx=1
          do time_out_idx=1,time_out_nbr
             if (src_frq(lon_idx,lat_idx,time_out_idx) /= 0.0) then
                ! Compute mean optical depth in each source region in each month
                src_str(lon_idx,lat_idx,time_out_idx)=src_str(lon_idx,lat_idx,time_out_idx)/src_frq(lon_idx,lat_idx,time_out_idx)
                ! Search for peak month
                if (src_str(lon_idx,lat_idx,time_out_idx) > src_str(lon_idx,lat_idx,time_max_idx)) time_max_idx=time_out_idx
             endif            ! endif src
          end do ! end loop over time
          if (src_flg(lon_idx,lat_idx) > 0) then
             ! Peak monthly value is efficiency factor
             fsh_fct(lon_idx,lat_idx)=src_str(lon_idx,lat_idx,time_max_idx) ! [frc]
             mth_max(lon_idx,lat_idx)=time_max_idx ! [mth]
             do time_out_idx=1,time_out_nbr
                ! Seasonal cycle is monthly source strength normalized by peak month
                csn_cyc(lon_idx,lat_idx,time_out_idx)=src_str(lon_idx,lat_idx,time_out_idx)/src_str(lon_idx,lat_idx,time_max_idx) ! [frc]
                ! src_odxc is optical depth in source regions averaged over all years (including non-source years)
                ! 19990128: Ensure src_odxc is 0.0 during months when no TOMS data exists (even if point is known source during other months)
                if (odxc_tms(lon_idx,lat_idx,time_out_idx) /= mss_val) then
                   src_odxc(lon_idx,lat_idx,time_out_idx)=odxc_tms(lon_idx,lat_idx,time_out_idx) ! [frc]
                else
                   src_odxc(lon_idx,lat_idx,time_out_idx)=0.0 ! [frc]
                endif         ! endif
                ! Convert src_frq from ordinal to fraction
                src_frq(lon_idx,lat_idx,time_out_idx)=src_frq(lon_idx,lat_idx,time_out_idx)/yr_nbr ! [nbr] -> [frc]
             end do ! end loop over time
          endif               ! endif src
       end do ! end loop over lat
    end do ! end loop over lon
    
    ! Default CCM source strength and seasonal cycle (iteration zero)
    itr_1st=.false. ! [flg] First iteration
    do lon_idx=1,lon_nbr
       do lat_idx=1,lat_nbr
          do time_out_idx=1,time_out_nbr
             csn_cyc_mdl(lon_idx,lat_idx,time_out_idx)=csn_cyc(lon_idx,lat_idx,time_out_idx)
             if (itr_1st) then
                src_str_mdl(lon_idx,lat_idx,time_out_idx)=src_odxc(lon_idx,lat_idx,time_out_idx)
             else
                src_str_mdl(lon_idx,lat_idx,time_out_idx)=src_str_old(lon_idx,lat_idx,time_out_idx)
             endif ! endif not first iteration
          end do ! end loop over time
       end do ! end loop over lat
    end do ! end loop over lon
    
    ! Create default soil type blend for saltation/sandblasting
    if(bln_nbr /= 4) stop 'bln_nbr /= 4'
    do lon_idx=1,lon_nbr
       do lat_idx=1,lat_nbr
          ! CMS = Coarse Medium Sand = 90% Coarse sand (CS) + 10% Fine Sand (FS)
          sfc_frc_bln(lon_idx,lat_idx,1)=0.0 ! [frc] Aluminosilicated silt fraction
          sfc_frc_bln(lon_idx,lat_idx,2)=0.1 ! [frc] Fine sand fraction
          sfc_frc_bln(lon_idx,lat_idx,3)=0.0 ! [frc] Salt fraction
          sfc_frc_bln(lon_idx,lat_idx,4)=0.9 ! [frc] Coarse sand fraction
       end do ! end loop over lat
    end do ! end loop over lon
    
    ! Perform linear, iterative correction to source strengths
    do lon_idx=1,lon_nbr
       do lat_idx=1,lat_nbr
          if (src_flg(lon_idx,lat_idx) > 0) then
             do time_out_idx=1,time_out_nbr
                if ((odxc_mdl(lon_idx,lat_idx,time_out_idx) > 0.0).and. &
                     (odxc_tms(lon_idx,lat_idx,time_out_idx) > 0.0).and. &
                     ! 19990128: Check for mss_val to ensure division does not cause floating point overflow error when TOMS datum is missing
                     ! 19990129: fxm I have no idea why mss_val check does not work
                     ! (odxc_tms(lon_idx,lat_idx,time_out_idx) /= mss_val)) then 
                     (odxc_tms(lon_idx,lat_idx,time_out_idx) < 100.0)) then 
                   itr_fct=odxc_tms(lon_idx,lat_idx,time_out_idx)/odxc_mdl(lon_idx,lat_idx,time_out_idx)
                   ! fxm: Apply arbitrary cutoffs to iterative correction
                   itr_fct=min(itr_fct,itr_fct_max)
                   itr_fct=max(itr_fct,itr_fct_min)
                   if(itr_fct > 100.0) then
                      write (6,"(a,f12.6)") "WARNING: bds_prc() reports itr_fct = ",itr_fct
                      write (6,"(a)") "Cell edge locations:"
                      write (6,"(a)") "(idx)   Lat Sth (idx)   Lat Nrt (idx)   Lon Wst (idx)   Lon Est"
                      write (6,"(4(a1,i3,a1,1x,f9.4,1x))") &
                           "(",lat_idx,")",lat_grd(lat_idx),"(",lat_idx+1,")",lat_grd(lat_idx+1), &
                           "(",lon_idx,")",lon_grd(lon_idx),"(",lon_idx+1,")",lon_grd(lon_idx+1)
                   endif ! endif
                   ! New value is iterative factor times old value
                   src_str_mdl(lon_idx,lat_idx,time_out_idx)= &
                        itr_fct*src_str_mdl(lon_idx,lat_idx,time_out_idx)
                endif ! endif
             end do ! end loop over time
          endif ! endif src
       end do ! end loop over lat
    end do ! end loop over lon
    
    bsn_nbr=lon_nbr*lat_nbr ! [nbr] Number of basins
    slat_nbr=lat_nbr-1 ! [nbr] Number of staggered latitudes for FV grid
    slon_nbr=lon_nbr ! [nbr] Number of staggered longitudes for FV grid
    call bds_out( &
         aod_nmd_frc, & ! [frc] TOMS AOD not due to dust
         area, & ! [m2] Area of gridcells 
         bln_nbr, & ! [nbr] Number of tri-modal soil blends
         bsn_enm, & ! [enm] Basin ID
         bsn_nbr, &  ! [nbr] Number of basins
         csn_cyc, & ! [frc] Seasonal cycle
         csn_cyc_mdl, & ! [frc] Seasonal cycle, CCM
         flw_acm_fct, & ! [frc] Flow accumulation factor
         fsh_fct, & ! [frc] Efficiency factor
         grd_mrd_lng, & ! [m] Grid meridional length
         grd_rds_rth, & ! [m] Distance from axis of rotation at local latitude
         grd_znl_lng, & ! [m] Grid zonal length
         gw, & ! Gaussian weight
         hgt_sfc_gpm,hgt_sfc,hgt_sfc_std_dvn,oro, & ! O 
         hgt_zpd_lsm, & ! [m] Zero plane displacement height
         lai_clm, & ! [m2 m-2] Leaf area index, one-sided
         lai_lsm, & ! [m2 m-2] Leaf area index
         lai_ttl_clm, & ! [m2 m-2] Total leaf area index, one-sided
         lak_frc, & ! Lake fraction of gridcell
         lat, & ! Coordinate variable
         lat_grd, & ! Latitude grid
         lat_nbr, & ! # Latitudes
         lat_sz, & ! Latitudinal size of bin
         lnd_frc, & ! [frc] Land fraction
         lnd_frc_clm, & ! [frc] Fraction of land (not ocean)
         lnd_frc_dry, & ! Dry fraction of gridcell
         lnd_msk, & ! [msk] Land mask
         lon, & ! Coordinate variable
         lon_grd, & ! Longitude grid
         lon_nbr, & ! # Longitudes
         lon_sz, & ! Latitudinal size of bin
         mbl_bsn_fct, & ! [frc] Mobilization enhancement due to basin characteristics
         msk_rgn, & ! [flg] Mask region
         mss_frc_Al, & ! [frc] Exchangeable aluminum
         mss_frc_C_org, & ! [frc] Organic carbon
         mss_frc_CaCO3, & ! [frc] Calcium carbonate
         mss_frc_K, & ! [frc] Exchangeable potassium
         mss_frc_N, & ! [frc] Nitrogen
         mss_frc_Na, & ! [frc] Exchangeable sodium
         mss_frc_P2O5, & ! [frc] Extractable phosphorus
         mss_frc_cly, & ! Soil texture clay
         mss_frc_slt, & ! Soil texture silt
         mss_frc_snd, & ! Soil texture sand
         mth_max, & ! [idx] Month of seasonal maximum
         odxc_tms, & ! [frc] Extinction optical depth proxy from TOMS
         pft_clm, & ! [enm] Plant functional type
         pft_frc_clm, & ! [enm] Fraction covered by plant functional type
         rfl_sfc_lnd_nrm_mds_lnr, & ![frc] Erodibility factor from MODIS, linear
         rfl_sfc_lnd_nrm_mds_sqr, & ![frc] Erodibility factor from MODIS, squared
         rgh_mmn, & ! [m] Roughness length humidity
         sfc_acm_fct, & ! [frc] Surface area accumulation factor
         sfc_frc_bln, & ! [m2 m-2] Surface area fraction of soil blend
         sfc_typ, & ! [enm] Surface type code
         sgs_nbr, & ! [nbr] Number of sub-gridscale patches per gridcell
         slat_nbr, & ! [nbr] Number of staggered latitudes for FV grid
         slon_nbr, & ! [nbr] Number of staggered longitudes for FV grid
         smp_sat, & ! [mm H2O] Saturated soil matric potential (sand-dependent)
         smp_sfc, & ! [mm H2O] Soil matric potential
         src_flg, & ! [flg] Source flag
         src_frq, & ! [frc] Source frequency
         src_odxc, & ! [frc] Extinction optical depth, source regions only
         src_str, & ! [frc] Source strength
         src_str_mdl, & ! [frc] Source strength, CCM
         src_thr, & ! [frc] Source threshold
         time_out, & ! Coordinate variable
         time_out_nbr, & ! [nbr] Dimension size
         vai_clm, & ! [m2 m-2] Vegetation area index, one-sided
         vai_lsm, & ! [m2 m-2] Vegetation area index
         vai_ttl_clm, & ! [m2 m-2] Total vegetation area index, one-sided
         vwc_dry, & ! [m3 m-3] Dry volumetric water content (no E-T)
         vwc_opt, & ! [m3 m-3] E-T optimal volumetric water content 
         vwc_rel, & ! [frc] Water content relative to saturation
         vwc_sat, & ! [m3 m-3] Saturated volumetric water content (sand-dependent)
         vwc_sfc, & ! [m3 m-3] Volumetric water content
         wtl_frc, & ! [frc] Wetland fraction of gridcell
         bsn_sz & ! [m2] Basin area
         )
    
    if (dbg_lvl >= dbg_sbr) write(6,"(a)") "Exiting bds_prc()"
    return 
  end subroutine bds_prc                       ! end bds_prc()
  
  subroutine lnd_bnd_cnv_rbn( & ! [sbr] Convert, rebin external datasets
       lat_out_grd,         & ! I
       lat_out_nbr,         & ! I
       lon_out_grd,         & ! I
       lon_out_nbr,         & ! I
       sgs_out_nbr, & ! I [nbr] Number of sub-gridscale patches per gridcell
       time_out_nbr, & ! I [nbr] Number of times
       area_out,            & ! I
       lak_frc_out,         & ! O
       lnd_frc_out,         & ! O
       lnd_frc_dry_out,     & ! O
       lnd_msk_out,         & ! O
       sfc_typ_out,         & ! O
       mss_frc_Al_out,      & ! O
       mss_frc_C_org_out,   & ! O
       mss_frc_CaCO3_out,   & ! O
       mss_frc_K_out,       & ! O
       mss_frc_N_out,       & ! O
       mss_frc_Na_out,      & ! O
       mss_frc_P2O5_out,    & ! O
       mss_frc_cly_out,     & ! O
       mss_frc_slt_out,     & ! O
       mss_frc_snd_out,     & ! O
       hgt_sfc_gpm_out,hgt_sfc_out,hgt_sfc_std_dvn_out,oro_out, & ! O
       mbl_bsn_fct_out,     & ! O
       bsn_enm_out,     & ! O [enm] Basin ID
       flw_acm_fct_out,     & ! O
       sfc_acm_fct_out,    & ! O
       rfl_sfc_lnd_nrm_mds_lnr_out,     & ! O [frc] Erodibility factor from MODIS, linear
       rfl_sfc_lnd_nrm_mds_sqr_out,     & ! O [frc] Erodibility factor from MODIS, squared
       pft_clm_out, & ! O [enm] Plant functional type
       pft_frc_clm_out, & ! O [enm] Fraction covered by plant functional type
       lnd_frc_clm_out, & ! O [frc] Fraction of land (not ocean)
       lai_clm_out, & ! O [m2 m-2] Leaf area index, one-sided
       lai_ttl_clm_out, & ! O [m2 m-2] Total leaf area index, one-sided
       vai_clm_out, & ! O [m2 m-2] Vegetation area index, one-sided
       vai_ttl_clm_out, & ! O [m2 m-2] Total vegetation area index, one-sided
       wtl_frc_out          & ! O
       )
    ! Purpose: Retrieve surface data from all external datasets and rebin it to output grid
    ! lnd_bnd_cnv_rbn() is called by bds_prc()
    use bds_ctl,only:fl_chm,fl_hgt,fl_lai,fl_lak,fl_pft, &
         fl_rfl_sfc,fl_sfc,fl_sfc_xtr,fl_soi_txt,fl_wtl, &
         flg_rdb_only,flg_sfc_xtr,flg_soi_nn ! [mdl] Control variables drc_in,drc_out,hgt_dlt_msl
    use clm_fld,only:lai_get,pft_get ! [mdl] Community Land Model fields
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use hgt_sfc,only:hgt_sfc_get,hgt_sfc_get_old ! [mdl] Surface elevation, bathymetry
    use lak_wtl,only:lnd_frc_dry_get ! [mdl] Lakes and wetlands
    use map_cst ! [mdl] Constants used in map routines
    use map_grd ! [mdl] Map grids and regridding
    use sfc_typ_mdl,only:sfc_typ_get,sfc_typ_xtr_get ! [mdl] Surface type and derived properties
    use sng_mdl,only:ftn_strlen ! [mdl] String manipulation
    use soi_mnr ! [mdl] IGBP-DIS geochemical data
    use soi_txt,only:soi_txt_get ! [mdl] Soil texture and hydrology
    use utl_mdl,only:mnt_chk ! [mdl] Utility functions (date_time_get,mnt_chk...)
    use asm_modis,only:bsn_mds_get ![mdl] Erodibility factor from MODIS satellite
    use bds_ctl,only:flg_rgr
    use rgr_mdl 
    implicit none
    ! Parameters
    integer,parameter::lat_in_nbr_max=1080 ! [nbr] Maximum number of latitudes
    integer,parameter::lon_in_nbr_max=2160 ! [nbr] Maximum number of longitudes
    character(len=*),parameter::sbr_nm="lnd_bnd_cnv_rbn" ! [sng] Subroutine name
    ! Commons
    ! Input
    integer,intent(in)::lon_out_nbr ! [nbr] Number of longitudes
    integer,intent(in)::lat_out_nbr ! [nbr] Number of latitudes
    integer,intent(in)::sgs_out_nbr ! [nbr] Number of sub-gridscale patches per gridcell
    integer,intent(in)::time_out_nbr ! [nbr] Number of times
    real,intent(in)::area_out(lon_out_nbr,lat_out_nbr) ! [m2] Area of gridcells 
    real,intent(in)::lon_out_grd(lon_out_nbr+1) ! [dgr] Interface longitudes
    real,intent(in)::lat_out_grd(lat_out_nbr+1) ! [dgr] Interface latitudes
    ! Output
    integer,intent(out)::bsn_enm_out(lon_out_nbr,lat_out_nbr) ! [enm] Basin ID
    integer,intent(out)::lnd_msk_out(lon_out_nbr,lat_out_nbr) ! [msk] Land mask (integer 0 or 1)
    integer,intent(out)::sfc_typ_out(lon_out_nbr,lat_out_nbr) ! [enm] Surface type 
    integer,intent(out)::pft_clm_out(lon_out_nbr,lat_out_nbr,sgs_out_nbr) ! [enm] Plant functional type
    real,intent(out)::pft_frc_clm_out(lon_out_nbr,lat_out_nbr,sgs_out_nbr) ! [enm] Fraction covered by plant functional type
    real,intent(out)::lnd_frc_clm_out(lon_out_nbr,lat_out_nbr) ! [frc] Fraction of land (not ocean)
    real,intent(out)::lai_clm_out(lon_out_nbr,lat_out_nbr,sgs_out_nbr,time_out_nbr) ! [m2 m-2] Leaf area index, one-sided
    real,intent(out)::lai_ttl_clm_out(lon_out_nbr,lat_out_nbr,time_out_nbr) ! [m2 m-2] Total leaf area index, one-sided
    real,intent(out)::vai_clm_out(lon_out_nbr,lat_out_nbr,sgs_out_nbr,time_out_nbr) ! [m2 m-2] Vegetation area index, one-sided
    real,intent(out)::vai_ttl_clm_out(lon_out_nbr,lat_out_nbr,time_out_nbr) ! [m2 m-2] Total vegetation area index, one-sided
    real,intent(out)::lnd_frc_out(lon_out_nbr,lat_out_nbr) ! Land fraction of gridcell 
    real,intent(out)::wtl_frc_out(lon_out_nbr,lat_out_nbr) ! Wetland fraction of gridcell 
    real,intent(out)::lak_frc_out(lon_out_nbr,lat_out_nbr) ! Lake fraction of gridcell 
    real,intent(out)::lnd_frc_dry_out(lon_out_nbr,lat_out_nbr) ! Dry fraction of gridcell 
    real,intent(out)::mss_frc_cly_out(lon_out_nbr,lat_out_nbr) ! [frc] Soil texture clay
    real,intent(out)::mss_frc_snd_out(lon_out_nbr,lat_out_nbr) ! [frc] Soil texture sand
    real,intent(out)::mss_frc_slt_out(lon_out_nbr,lat_out_nbr) ! [frc] Soil texture silt
    real,intent(out)::mss_frc_Al_out(lon_out_nbr,lat_out_nbr) ! [frc] Exchangeable aluminum
    real,intent(out)::mss_frc_C_org_out(lon_out_nbr,lat_out_nbr) ! [frc] Organic carbon
    real,intent(out)::mss_frc_CaCO3_out(lon_out_nbr,lat_out_nbr) ! [frc] Calcium carbonate
    real,intent(out)::mss_frc_K_out(lon_out_nbr,lat_out_nbr) ! [frc] Exchangeable potassium
    real,intent(out)::mss_frc_N_out(lon_out_nbr,lat_out_nbr) ! [frc] Nitrogen
    real,intent(out)::mss_frc_Na_out(lon_out_nbr,lat_out_nbr) ! [frc] Exchangeable sodium
    real,intent(out)::mss_frc_P2O5_out(lon_out_nbr,lat_out_nbr) ! [frc] Extractable phosphorus
    real,intent(out)::oro_out(lon_out_nbr,lat_out_nbr) ! [frc] Orography 
    real,intent(out)::hgt_sfc_gpm_out(lon_out_nbr,lat_out_nbr) ! [gpm] Surface geopotential height 
    real,intent(out)::hgt_sfc_out(lon_out_nbr,lat_out_nbr) ! [m] Surface height 
    real,intent(out)::hgt_sfc_std_dvn_out(lon_out_nbr,lat_out_nbr) ! [m] Standard deviation of surface height 
    real,intent(out)::mbl_bsn_fct_out(lon_out_nbr,lat_out_nbr) ! [frc] Mobilization enhancement due to basin characteristics
    real,intent(out)::flw_acm_fct_out(lon_out_nbr,lat_out_nbr) ! [frc] Flow accumulation factor
    real,intent(out)::sfc_acm_fct_out(lon_out_nbr,lat_out_nbr) ! [frc] Surface area accumulation factor
    real,intent(out)::rfl_sfc_lnd_nrm_mds_lnr_out(lon_out_nbr,lat_out_nbr) ! [frc] Erodibility factor from MODIS, linear
    real,intent(out)::rfl_sfc_lnd_nrm_mds_sqr_out(lon_out_nbr,lat_out_nbr) ! [frc] Erodibility factor from MODIS, squared

    ! Local
    character(80)::fl_in        ! [sng] Input file
    integer chm_idx           ! [idx] Counting index for chm
    integer lat_in_nbr        ! [nbr] Number of latitudes
    integer lon_in_nbr        ! [nbr] Number of longitudes
    integer map_typ_in        ! [enm] Input map grid type
    integer ovr_nbr_max       ! [nbr] Maximum number of input cells which overlap any output cell
    logical mnt               ! [flg] Monotonicity flag
    real lat_in_grd(lat_in_nbr_max+1) ! [dgr] Interface latitudes
    real lon_in_grd(lon_in_nbr_max+1) ! [dgr] Interface longitudes
    real mss_frc_chm_out(lon_out_nbr,lat_out_nbr) ! [kg kg-1] Generic chemical
    real ovr_nbr_val(lon_out_nbr,lat_out_nbr)
    integer map_typ_out


    ! Main Code
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Entering"//sbr_nm//"()"
    
    ! Check for monotonicity
    mnt=mnt_chk(lon_out_grd,lon_out_nbr+1)
    if (.not.mnt) stop "lon_out_grd not monotonic in lnd_bnd_cnv_rbn()"
    mnt=mnt_chk(lat_out_grd,lat_out_nbr+1)
    if (.not.mnt) stop "lat_out_grd not monotonic in lnd_bnd_cnv_rbn()"
    ! Check for increasing monotonicity
    if (lon_out_grd(2) < lon_out_grd(1)) stop "lon_out_grd not increasing in lnd_bnd_cnv_rbn()"
    if (lat_out_grd(2) < lat_out_grd(1)) stop "lat_out_grd not increasing in lnd_bnd_cnv_rbn()"
    
    ! Surface elevation
    ! NB: hgt_sfc_get() MUST be called first since it provides inputs (lnd_frc) for following routines
    call hgt_sfc_get(fl_hgt, & ! I
         lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
         area_out, & ! I
         hgt_sfc_gpm_out,hgt_sfc_out,hgt_sfc_std_dvn_out, & ! O
         lnd_frc_out,lnd_msk_out,oro_out, & ! O
         mbl_bsn_fct_out,bsn_enm_out,flw_acm_fct_out,sfc_acm_fct_out) ! O
    ! Deprecated method for surface elevation
    if(.false.) then
       lat_in_nbr=1080 ! 1/6 x 1/6 degree hgt_sfc_topo.nc 
       lon_in_nbr=2160 ! 1/6 x 1/6 degree hgt_sfc_topo.nc
       map_typ_in=map_lat_rgl_lon_Grn_wst ! Latitudes are regular, Greenwich at West edge of first longitude bin
       call map_grd_mk(lat_in_nbr,lon_in_nbr,map_typ_in, & ! I
            lat_in_grd,lon_in_grd) ! O
       call map_ovr_nbr_max_get(lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
            lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
            ovr_nbr_max) ! O
       call hgt_sfc_get_old(fl_hgt, & ! I
            lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
            lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
            ovr_nbr_max,area_out, & ! I
            hgt_sfc_gpm_out,hgt_sfc_out,hgt_sfc_std_dvn_out, & ! O
            lnd_frc_out,lnd_msk_out,oro_out, & ! O
            mbl_bsn_fct_out) ! O
    endif ! endif .false.
    
    ! Skip remaining routines if flg_rdb_only true
    if (.not.flg_rdb_only) then
       
       ! Plant functional type
       call pft_get(fl_pft, & ! I
            lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr,sgs_out_nbr, & ! I
            area_out,lnd_msk_out, & ! I
            lnd_frc_clm_out,pft_clm_out,pft_frc_clm_out) ! O
       
       ! Leaf area index
       call lai_get(fl_lai, & ! I
            lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr,sgs_out_nbr,time_out_nbr, & ! I
            area_out,lnd_msk_out,pft_clm_out,pft_frc_clm_out, & ! I
            lai_clm_out,lai_ttl_clm_out,vai_clm_out,vai_ttl_clm_out) ! O
       
       ! Surface type
       lat_in_nbr=360 ! [nbr] 0.5 degree Olson DS769
       lon_in_nbr=720 ! [nbr] 0.5 degree Olson DS769
       map_typ_in=map_lat_rgl_lon_180_wst ! Latitudes are regular, Date line at West edge of first longitude bin
       call map_grd_mk(lat_in_nbr,lon_in_nbr,map_typ_in, & ! I
            lat_in_grd,lon_in_grd) ! O
       call map_ovr_nbr_max_get(lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
            lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
            ovr_nbr_max) ! O
       call sfc_typ_get(fl_sfc, & ! I
            lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
            lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
            ovr_nbr_max,area_out,lnd_frc_out,lnd_msk_out,flg_soi_nn, & ! I
            sfc_typ_out) ! O
       
       ! Read surface type directly from external file
       if (flg_sfc_xtr) then
          call sfc_typ_xtr_get( & ! [sbr] Retrieve and rebin surface type
               fl_sfc_xtr, & ! I
               lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
               area_out,lnd_frc_out,lnd_msk_out,flg_soi_nn, & ! I
               sfc_typ_out) ! O
       end if ! endif flg_sfc_xtr

       ! Soil texture
       ! NB: soi_txt_webb.1x1 and soi_txt_IBIS.nc are currently same resolution and map grid
       lat_in_nbr=180 ! [nbr] 1 x 1 degree
       lon_in_nbr=360 ! [nbr] 1 x 1 degree
       map_typ_in=map_lat_rgl_lon_180_wst ! Latitudes are regular, Date line at West edge of first longitude bin
       call map_grd_mk(lat_in_nbr,lon_in_nbr,map_typ_in, & ! I
            lat_in_grd,lon_in_grd) ! O
       call map_ovr_nbr_max_get(lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
            lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
            ovr_nbr_max) ! O
       call soi_txt_get(fl_soi_txt, & ! I
            lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
            lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
            ovr_nbr_max,area_out,lnd_msk_out,sfc_typ_out,flg_soi_nn, & ! I
            mss_frc_cly_out,mss_frc_slt_out,mss_frc_snd_out) ! O

       ! Basin factors from MODIS surface reflectance
       lat_in_nbr=720 ! [nbr] 0.25 x 0.25 degree
       lon_in_nbr=1440 ! [nbr] 0.25 x 0.25 degree
       map_typ_in=map_lat_rgl_lon_Grn_wst ! Latitudes are regular, Greenwich at west of first longitude bin
       call map_grd_mk(lat_in_nbr,lon_in_nbr,map_typ_in, & ! I
            lat_in_grd,lon_in_grd) ! O
       call map_ovr_nbr_max_get(lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
            lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
            ovr_nbr_max) ! O
       call bsn_mds_get(fl_rfl_sfc, & ! I
            lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
            lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
            ovr_nbr_max,area_out,lnd_frc_out, & ! I
            rfl_sfc_lnd_nrm_mds_lnr_out,rfl_sfc_lnd_nrm_mds_sqr_out) !O

       ! Lake and wetland fractions
       lat_in_nbr=180            ! 1.0 x 1.0 degree Cogley
       lon_in_nbr=360            ! 1.0 x 1.0 degree Cogley
       map_typ_in=map_lat_rgl_lon_Grn_ctr ! Latitudes are regular, Greenwich at center of first longitude bin
       call map_grd_mk(lat_in_nbr,lon_in_nbr,map_typ_in, & ! I
            lat_in_grd,lon_in_grd) ! O
       call map_ovr_nbr_max_get(lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
            lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
            ovr_nbr_max) ! O
       call lnd_frc_dry_get(fl_lak,fl_wtl, & ! I
            lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
            lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
            ovr_nbr_max,area_out,lnd_frc_out,sfc_typ_out, & ! I
            lak_frc_out,lnd_frc_dry_out,wtl_frc_out) ! O
       
       ! Soil chemical composition from IGBP-DIS
       lat_in_nbr=180            ! 1 x 1 degree
       lon_in_nbr=360            ! 1 x 1 degree
       map_typ_in=map_lat_rgl_lon_180_wst ! Latitudes are regular, Date line at West edge of first longitude bin
       call map_grd_mk(lat_in_nbr,lon_in_nbr,map_typ_in, & ! I
            lat_in_grd,lon_in_grd) ! O
       call map_ovr_nbr_max_get(lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
            lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
            ovr_nbr_max) ! O
       
       do chm_idx=1,chm_nbr
          call soi_chm_IGBP_get(fl_chm,chm_idx, & ! I
               lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
               lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
               ovr_nbr_max,area_out,lnd_msk_out, & ! I
               mss_frc_chm_out) ! O
          if (chm_idx==idx_CaCO3) then
             mss_frc_CaCO3_out=mss_frc_chm_out ! [kg kg-1] Calcium carbonate
          else if (chm_idx==idx_C_org) then
             mss_frc_C_org_out=mss_frc_chm_out ! [kg kg-1] Organic carbon
          else if (chm_idx==idx_P2O5) then
             mss_frc_P2O5_out=mss_frc_chm_out ! [kg kg-1] Extractable phosphorus
          else if (chm_idx==idx_Na) then
             mss_frc_Na_out=mss_frc_chm_out ! [kg kg-1] Exchangeable sodium
          else if (chm_idx==idx_K) then
             mss_frc_K_out=mss_frc_chm_out ! [kg kg-1] Exchangeable potassium
          else if (chm_idx==idx_Al) then
             mss_frc_Al_out=mss_frc_chm_out ! [kg kg-1] Exchangeable aluminum
          else if (chm_idx==idx_N) then
             mss_frc_N_out=mss_frc_chm_out ! [kg kg-1] Nitrogen
          else
             stop "Unknown chm_idx in lnd_bnd_cnv_rbn()"
          endif                  ! endif
       end do ! end loop over chm
       
    endif ! flg_rdb_only

    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Exiting"//sbr_nm//"()"
    return 
  end subroutine lnd_bnd_cnv_rbn
  
  subroutine bds_out( & ! [sbr] Write BDS output to netCDF file
       aod_nmd_frc, & ! [frc] TOMS AOD not due to dust
       area, & ! [m2] Area of gridcells 
       bln_nbr, & ! [nbr] Number of tri-modal soil blends
       bsn_enm, & ! [enm] Basin ID
       bsn_nbr, & ! [nbr] Number of basins
       csn_cyc, & ! [frc] Seasonal cycle
       csn_cyc_mdl, & ! [frc] Seasonal cycle, CCM
       flw_acm_fct, & ! [frc] Flow accumulation factor
       fsh_fct, & ! [frc] Efficiency factor
       grd_mrd_lng, & ! [m] Grid meridional length
       grd_rds_rth, & ! [m] Distance from axis of rotation at local latitude
       grd_znl_lng, & ! [m] Grid zonal length
       gw, & ! Gaussian weight
       hgt_sfc_gpm,hgt_sfc,hgt_sfc_std_dvn,oro, & ! O 
       hgt_zpd_lsm, & ! [m] Zero plane displacement height
       lai_clm, & ! [m2 m-2] Leaf area index, one-sided
       lai_lsm, & ! [m2 m-2] Leaf area index
       lai_ttl_clm, & ! [m2 m-2] Total leaf area index, one-sided
       lak_frc, & ! Lake fraction of gridcell
       lat, & ! Coordinate variable
       lat_grd, & ! Latitude grid
       lat_nbr, & ! # Latitudes
       lat_sz, & ! Latitudinal size of bin
       lnd_frc, & ! [frc] Land fraction
       lnd_frc_clm, & ! [frc] Fraction of land (not ocean)
       lnd_frc_dry, & ! Dry fraction of gridcell
       lnd_msk, & ! [msk] Land mask
       lon, & ! Coordinate variable
       lon_grd, & ! Longitude grid
       lon_nbr, & ! # Longitudes
       lon_sz, & ! Latitudinal size of bin
       mbl_bsn_fct, & ! [frc] Mobilization enhancement due to basin characteristics
       msk_rgn, & ! [flg] Mask region
       mss_frc_Al, & ! [frc] Exchangeable aluminum
       mss_frc_C_org, & ! [frc] Organic carbon
       mss_frc_CaCO3, & ! [frc] Calcium carbonate
       mss_frc_K, & ! [frc] Exchangeable potassium
       mss_frc_N, & ! [frc] Nitrogen
       mss_frc_Na, & ! [frc] Exchangeable sodium
       mss_frc_P2O5, & ! [frc] Extractable phosphorus
       mss_frc_cly, & ! Soil texture clay
       mss_frc_slt, & ! Soil texture silt
       mss_frc_snd, & ! Soil texture sand
       mth_max, & ! [idx] Month of seasonal maximum
       odxc_tms, & ! [frc] Extinction optical depth proxy from TOMS
       pft_clm, & ! [enm] Plant functional type
       pft_frc_clm, & ! [enm] Fraction covered by plant functional type
       rfl_sfc_lnd_nrm_mds_lnr, & ![frc] Erodibility factor from MODIS, linear
       rfl_sfc_lnd_nrm_mds_sqr, & ![frc] Erodibility factor from MODIS, squared
       rgh_mmn, & ! [m] Roughness length humidity
       sfc_acm_fct, & ! [frc] Surface area accumulation factor
       sfc_frc_bln, & ! [m2 m-2] Surface area fraction of soil blend
       sfc_typ, & ! [enm] Surface type code
       sgs_nbr, & ! [nbr] Number of sub-gridscale patches per gridcell
       slat_nbr, & ! [nbr] Number of staggered latitudes for FV grid
       slon_nbr, & ! [nbr] Number of staggered longitudes for FV grid
       smp_sat, & ! [mm H2O] Saturated soil matric potential (sand-dependent)
       smp_sfc, & ! [mm H2O] Soil matric potential
       src_flg, & ! [flg] Source flag
       src_frq, & ! [frc] Source frequency
       src_odxc, & ! [frc] Extinction optical depth, source regions only
       src_str, & ! [frc] Source strength
       src_str_mdl, & ! [frc] Source strength, CCM
       src_thr, & ! [frc] Source strength threshold
       time_out, & ! Coordinate variable
       time_out_nbr, & ! [nbr] Dimension size
       vai_clm, & ! [m2 m-2] Vegetation area index, one-sided
       vai_lsm, & ! [m2 m-2] Vegetation area index
       vai_ttl_clm, & ! [m2 m-2] Total vegetation area index, one-sided
       vwc_dry, & ! [m3 m-3] Dry volumetric water content (no E-T)
       vwc_opt, & ! [m3 m-3] E-T optimal volumetric water content 
       vwc_rel, & ! [frc] Water content relative to saturation
       vwc_sat, & ! [m3 m-3] Saturated volumetric water content (sand-dependent)
       vwc_sfc, & ! [m3 m-3] Volumetric water content
       wtl_frc, & ! [frc] Wetland fraction of gridcell
       bsn_sz & ! [m2] Basin area
       )
    ! Purpose: Write BDS output to netCDF file
    ! bds_out() is called by bds_prc()
    use bds_ctl,only:fl_out ! [mdl] Control variables drc_in,drc_out,hgt_dlt_msl
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use netcdf ! [mdl] netCDF interface
    use nf90_utl ! [mdl] netCDF utilities
    use shr_kind_mod,only:r8=>shr_kind_r8,r4=>shr_kind_r4,DBLKIND ! [mdl] Precision r8, i8, ...
    implicit none
    ! Parameters
    character(len=*),parameter::sbr_nm="bds_out" ! [sng] Subroutine name
    ! Commons
    ! Input
    integer,intent(in)::bln_nbr ! [nbr] Number of tri-modal soil blends
    integer,intent(in)::sgs_nbr ! [nbr] Number of sub-gridscale patches per gridcell
    ! Input/Output
    ! Output
    ! Local

    ! Metadata
    integer bsn_nbr ! [nbr] Number of basins
    integer dmn_bln(1) ! [id] Dimension ID for bln
    integer dmn_bsn(1) ! [id] Dimension ID for bsn
    integer dmn_lat(1) ! [id] Dimension ID for lat
    integer dmn_lat_grd(1) ! [id] Dimension ID for lat grid
    integer dmn_lon(1) ! [id] Dimension ID for lon
    integer dmn_lon_grd(1) ! [id] Dimension ID for lon grid
    integer dmn_lon_lat(2) ! [id] Dimension IDs
    integer dmn_lon_lat_bln(3) ! [id] Dimension IDs
    integer dmn_lon_lat_sgs(3) ! [id] Dimension IDs
    integer dmn_lon_lat_sgs_time(4) ! [id] Dimension IDs
    integer dmn_lon_lat_time(3) ! [id] Dimension IDs
    integer dmn_sgs(1) ! [id] Dimension ID for sgs
    integer dmn_slat(1) ! [id] Dimension ID for slat
    integer dmn_slon(1) ! [id] Dimension ID for slon
    integer dmn_time(1) ! [id] Dimension ID for time
    integer fll_mode_old ! Old fill mode
    integer lat_nbr ! [nbr] Number of latitudes
    integer lon_nbr ! [nbr] Number of longitudes
    integer nc_id ! [id] File handle
    integer rcd ! [enm] Return success code
    integer slat_nbr ! [nbr] Number of staggered latitudes for FV grid
    integer slon_nbr ! [nbr] Number of staggered longitudes for FV grid
    integer time_in_nbr ! [nbr] Dimension size
    integer time_out_nbr ! [nbr] Dimension size
    
    ! netCDF IDs
    integer aod_nmd_frc_id    ! [id] Variable ID
    integer csn_cyc_id        ! [id] Variable ID
    integer csn_cyc_mdl_id    ! [id] Variable ID
    integer lnd_frc_dry_id    ! [id] Variable ID
    integer fsh_fct_id        ! [id] Variable ID
    integer src_thr_id        ! [id] Variable ID
    integer grd_mrd_lng_id    ! [id] Variable ID
    integer grd_rds_rth_id    ! [id] Variable ID
    integer grd_znl_lng_id    ! [id] Variable ID
    integer gw_id             ! [id] Variable ID
    integer lai_lsm_id        ! [id] Variable ID
    integer lak_frc_id        ! [id] Variable ID
    integer lat_fv_id        ! [id] Coordinate ID
    integer lat_grd_id        ! [id] Coordinate ID
    integer lat_id            ! [id] Coordinate ID
    integer slon_id            ! [id] Coordinate ID
    integer slat_id            ! [id] Coordinate ID
    integer lat_sz_id         ! [id] Variable ID
    integer lnd_frc_id        ! [id] Variable ID
    integer lnd_msk_id        ! [id] Variable ID
    integer lon_grd_id        ! [id] Coordinate ID
    integer lon_id            ! [id] Coordinate ID
    integer lon_sz_id         ! [id] Variable ID
    integer mbl_bsn_fct_id    ! [id] Variable ID
    integer bsn_enm_id         ! [id] Variable ID
    integer bsn_sz_id        ! [id] Variable ID
    integer flw_acm_fct_id    ! [id] Variable ID
    integer sfc_acm_fct_id   ! [id] Variable ID
    integer rfl_sfc_lnd_nrm_mds_lnr_id    ! [id] Variable ID
    integer rfl_sfc_lnd_nrm_mds_sqr_id    ! [id] Variable ID
    integer area_id           ! [id] Variable ID
    integer msk_rgn_id        ! [id] Variable ID
    integer mss_frc_Al_id     ! [id] Variable ID
    integer mss_frc_C_org_id  ! [id] Variable ID
    integer mss_frc_CaCO3_id  ! [id] Variable ID
    integer mss_frc_K_id      ! [id] Variable ID
    integer mss_frc_N_id      ! [id] Variable ID
    integer mss_frc_Na_id     ! [id] Variable ID
    integer mss_frc_P2O5_id   ! [id] Variable ID
    integer mss_frc_cly_id    ! [id] Variable ID
    integer mss_frc_slt_id    ! [id] Variable ID
    integer mss_frc_snd_id    ! [id] Variable ID
    integer mth_max_id        ! [id] Variable ID
    integer odxc_tms_id       ! [id] Variable ID
    integer rgh_mmn_id        ! [id] Variable ID
    integer sfc_typ_id        ! [id] Variable ID
    integer smp_sat_id        ! [id] Variable ID
    integer smp_sfc_id        ! [id] Variable ID
    integer src_flg_id        ! [id] Variable ID
    integer src_frq_id        ! [id] Variable ID
    integer src_odxc_id       ! [id] Variable ID
    integer src_str_mdl_id    ! [id] Variable ID
    integer src_str_id        ! [id] Variable ID
    integer time_id           ! [id] Coordinate ID
    integer vai_lsm_id        ! [id] Variable ID
    integer vwc_dry_id        ! [id] Variable ID
    integer vwc_opt_id        ! [id] Variable ID
    integer vwc_rel_id        ! [id] Variable ID
    integer vwc_sat_id        ! [id] Variable ID
    integer vwc_sfc_id        ! [id] Variable ID
    integer wtl_frc_id        ! [id] Variable ID
    integer hgt_zpd_lsm_id    ! [id] Variable ID
    integer oro_id            ! [id] Variable ID 
    integer hgt_sfc_gpm_id    ! [id] Variable ID 
    integer hgt_sfc_id        ! [id] Variable ID 
    integer hgt_sfc_std_dvn_id ! [id] Variable ID 
    integer pft_clm_id ! [id] Variable ID
    integer lai_clm_id ! [id] Variable ID
    integer lai_ttl_clm_id ! [id] Variable ID
    integer vai_clm_id ! [id] Variable ID
    integer vai_ttl_clm_id ! [id] Variable ID
    integer sfc_frc_bln_id ! [id] Variable ID
    integer pft_frc_clm_id ! [id] Variable ID
    integer lnd_frc_clm_id ! [id] Variable ID
    
    ! Data
    integer bsn_enm(lon_nbr,lat_nbr) ! [enm] Basin ID
    integer lnd_msk(lon_nbr,lat_nbr) ! [msk] Land mask
    integer mth_max(lon_nbr,lat_nbr) ! [idx] Month of seasonal maximum
    integer msk_rgn(lon_nbr,lat_nbr) ! [flg] Mask region
    integer sfc_typ(lon_nbr,lat_nbr) ! [enm] Surface type code
    integer src_flg(lon_nbr,lat_nbr) ! [flg] Source flag
    real aod_nmd_frc(lon_nbr,lat_nbr,time_out_nbr) ! [frc] TOMS AOD not due to dust
    real area(lon_nbr,lat_nbr) ! [m2] Area of gridcells 
    real csn_cyc(lon_nbr,lat_nbr,time_out_nbr) ! [frc] Seasonal cycle
    real csn_cyc_mdl(lon_nbr,lat_nbr,time_out_nbr) ! [frc] Seasonal cycle
    real fsh_fct(lon_nbr,lat_nbr) ! [frc] Efficiency factor
    real grd_mrd_lng(lat_nbr) ! [m] Grid meridional length
    real grd_rds_rth(lat_nbr) ! [m] Distance from axis of rotation at local latitude
    real grd_znl_lng(lat_nbr) ! [m] Grid zonal length
    real gw(lat_nbr) ! Gaussian weight
    real hgt_zpd_lsm(lon_nbr,lat_nbr) ! [m] Zero plane displacement height
    real lai_lsm(lon_nbr,lat_nbr,time_out_nbr) ! [m2 m-2] Leaf area index
    real lak_frc(lon_nbr,lat_nbr) ! [frc] Lake fraction of gridcell
    real lat(lat_nbr) ! Coordinate variable
    real slat(slat_nbr) ! Coordinate variable
    real slon(slon_nbr) ! Coordinate variable
    real lat_fv(lat_nbr) ! [dgr] Dynamic latitudes for FV grid
    real lat_grd(lat_nbr+1) ! Latitude grid
    real lat_sz(lat_nbr) ! Latitudinal size of bin
    real lnd_frc(lon_nbr,lat_nbr) ! [frc] Land fraction
    real lnd_frc_dry(lon_nbr,lat_nbr) ! [frc] Dry fraction of gridcell
    real lon(lon_nbr) ! Coordinate variable
    real lon_grd(lon_nbr+1) ! Longitude grid
    real lon_sz(lon_nbr) ! Latitudinal size of bin
    real mbl_bsn_fct(lon_nbr,lat_nbr) ! [frc] Mobilization enhancement due to basin characteristics
    real flw_acm_fct(lon_nbr,lat_nbr) ! [frc] Flow accumulation factor
    real sfc_acm_fct(lon_nbr,lat_nbr) ! [frc] Surface area accumulation factor
    real rfl_sfc_lnd_nrm_mds_lnr(lon_nbr,lat_nbr) ![frc] Erodibility factor from MODIS, linear
    real rfl_sfc_lnd_nrm_mds_sqr(lon_nbr,lat_nbr) ![frc] Erodibility factor from MODIS, squared
    real mss_frc_Al(lon_nbr,lat_nbr) ! [frc] Exchangeable aluminum
    real mss_frc_C_org(lon_nbr,lat_nbr) ! [frc] Organic carbon
    real mss_frc_CaCO3(lon_nbr,lat_nbr) ! [frc] Calcium carbonate
    real mss_frc_K(lon_nbr,lat_nbr) ! [frc] Exchangeable potassium
    real mss_frc_N(lon_nbr,lat_nbr) ! [frc] Nitrogen
    real mss_frc_Na(lon_nbr,lat_nbr) ! [frc] Exchangeable sodium
    real mss_frc_P2O5(lon_nbr,lat_nbr) ! [frc] Extractable phosphorus
    real sfc_frc_bln(lon_nbr,lat_nbr,bln_nbr) ! [m2 m-2] Surface area fraction of soil blend
    real mss_frc_cly(lon_nbr,lat_nbr) ! [frc] Soil texture clay
    real mss_frc_slt(lon_nbr,lat_nbr) ! [frc] Soil texture silt
    real mss_frc_snd(lon_nbr,lat_nbr) ! [frc] Soil texture sand
    real odxc_tms(lon_nbr,lat_nbr,time_out_nbr) ! [frc] Extinction optical depth proxy from TOMS
    real oro(lon_nbr,lat_nbr) ! [frc] Orography 
    real rgh_mmn(lon_nbr,lat_nbr) ! [m] Roughness length humidity
    real hgt_sfc(lon_nbr,lat_nbr) ! [m] Surface height 
    real hgt_sfc_gpm(lon_nbr,lat_nbr) ! [gpm] Surface geopotential height 
    real hgt_sfc_std_dvn(lon_nbr,lat_nbr) ! [m] Standard deviation of surface height 
    real smp_sat(lon_nbr,lat_nbr) ! [mm H2O] Saturated soil matric potential (sand-dependent)
    real smp_sfc(lon_nbr,lat_nbr,time_out_nbr) ! [mm H2O] Soil matric potential
    real src_frq(lon_nbr,lat_nbr,time_out_nbr) ! [frc] Source frequency
    real src_odxc(lon_nbr,lat_nbr,time_out_nbr) ! [frc] Extinction optical depth, source regions only
    real src_str(lon_nbr,lat_nbr,time_out_nbr) ! [frc] Source strength
    real src_str_mdl(lon_nbr,lat_nbr,time_out_nbr) ! [frc] Source strength, CCM
    real src_thr              ! [frc] Threshold for source maxima
    real time_out(time_out_nbr) ! Coordinate variable
    real vai_lsm(lon_nbr,lat_nbr,time_out_nbr) ! [m2 m-2] Vegetation area index
    real vwc_dry(lon_nbr,lat_nbr) ! [m3 m-3] Dry volumetric water content (no E-T)
    real vwc_opt(lon_nbr,lat_nbr) ! [m3 m-3] E-T optimal volumetric water content 
    real vwc_rel(lon_nbr,lat_nbr,time_out_nbr) ! [frc] Water content relative to saturation
    real vwc_sat(lon_nbr,lat_nbr) ! [m3 m-3] Saturated volumetric water content (sand-dependent)
    real vwc_sfc(lon_nbr,lat_nbr,time_out_nbr) ! [m3 m-3] Volumetric water content
    integer pft_clm(lon_nbr,lat_nbr,sgs_nbr) ! [enm] Plant functional type
    real pft_frc_clm(lon_nbr,lat_nbr,sgs_nbr) ! [enm] Fraction covered by plant functional type
    real lnd_frc_clm(lon_nbr,lat_nbr) ! [frc] Fraction of land (not ocean)
    real lai_clm(lon_nbr,lat_nbr,sgs_nbr,time_out_nbr) ! [m2 m-2] Leaf area index, one-sided
    real lai_ttl_clm(lon_nbr,lat_nbr,time_out_nbr) ! [m2 m-2] Total leaf area index, one-sided
    real vai_clm(lon_nbr,lat_nbr,sgs_nbr,time_out_nbr) ! [m2 m-2] Vegetation area index, one-sided
    real vai_ttl_clm(lon_nbr,lat_nbr,time_out_nbr) ! [m2 m-2] Total vegetation area index, one-sided
    real wtl_frc(lon_nbr,lat_nbr) ! [frc] Wetland fraction of gridcell
    real bsn_sz(bsn_nbr) ! [m2] Basin area

    ! Local
    real(r4),parameter::mss_val=1.0e36 ! [frc] Missing value
    ! Main code
    if (dbg_lvl >= dbg_sbr) write(6,"(a)") "Entering bds_out()"
    
    ! Initialize some variables that should be initialized in grid routines eventually
    ! CAM FV grids slat does not include interface values at poles
    slat=lat_grd(2:lat_nbr) ! Coordinate variable 
    slon=lon_grd(1:lon_nbr) ! Coordinate variable
    lat_fv(2:lat_nbr-1)=lat(2:lat_nbr-1) ! [dgr] Dynamic latitudes for FV grid
    lat_fv(1)=lat_grd(1) ! [dgr] Dynamic latitudes for FV grid
    lat_fv(lat_nbr)=lat_grd(lat_nbr+1) ! [dgr] Dynamic latitudes for FV grid

    ! Begin netCDF output routines
    rcd=nf90_noerr              ! nf90_noerr == 0
    rcd=nf90_wrp_open(fl_out,nf90_write,nc_id,sbr_nm=sbr_nm)
    ! Put output file in define mode
    rcd=nf90_wrp(nf90_redef(nc_id),sbr_nm//': nf90_redef')
    rcd=nf90_wrp(nf90_set_fill(nc_id,nf90_nofill,fll_mode_old),sbr_nm//': nf90_set_fill')
    ! Define dimension IDs
    rcd=nf90_wrp(nf90_def_dim(nc_id,"time",nf90_unlimited,dmn_time(1)),sbr_nm//": def_dim time")
    rcd=nf90_wrp(nf90_def_dim(nc_id,"lat",lat_nbr,dmn_lat(1)),sbr_nm//": def_dim lat")
    rcd=nf90_wrp(nf90_def_dim(nc_id,"lon",lon_nbr,dmn_lon(1)),sbr_nm//": def_dim lon")
    rcd=nf90_wrp(nf90_def_dim(nc_id,"sgs",sgs_nbr,dmn_sgs(1)),sbr_nm//": def_dim sgs")
    rcd=nf90_wrp(nf90_def_dim(nc_id,"bln",bln_nbr,dmn_bln(1)),sbr_nm//": def_dim bln")
    rcd=nf90_wrp(nf90_def_dim(nc_id,"lat_grd",lat_nbr+1,dmn_lat_grd(1)),sbr_nm//": def_dim lat_grd")
    rcd=nf90_wrp(nf90_def_dim(nc_id,"lon_grd",lon_nbr+1,dmn_lon_grd(1)),sbr_nm//": def_dim lon_grd")
    rcd=nf90_wrp(nf90_def_dim(nc_id,"slat",slat_nbr,dmn_slat(1)),sbr_nm//": def_dim slat")
    rcd=nf90_wrp(nf90_def_dim(nc_id,"slon",slon_nbr,dmn_slon(1)),sbr_nm//": def_dim slon")
    rcd=nf90_wrp(nf90_def_dim(nc_id,"bsn",bsn_nbr,dmn_bsn(1)),sbr_nm//": def_dim bsn")

    ! Assemble ID and count vectors for each multidimensional combination of dimensions
    dmn_lon_lat=(/dmn_lon(1),dmn_lat(1)/)
    dmn_lon_lat_time=(/dmn_lon(1),dmn_lat(1),dmn_time(1)/)
    dmn_lon_lat_sgs=(/dmn_lon(1),dmn_lat(1),dmn_sgs(1)/)
    dmn_lon_lat_bln=(/dmn_lon(1),dmn_lat(1),dmn_bln(1)/)
    dmn_lon_lat_sgs_time=(/dmn_lon(1),dmn_lat(1),dmn_sgs(1),dmn_time(1)/)
    ! Variable definitions
    ! Wrapped
    rcd=nf90_wrp(nf90_def_var(nc_id,"aod_nmd_frc",nf90_float,dmn_lon_lat_time,aod_nmd_frc_id), &
         & sbr_nm//": def_var aod_nmd_frc")
    rcd=nf90_wrp(nf90_def_var(nc_id,"rfl_sfc_lnd_nrm_mds_lnr",nf90_float,dmn_lon_lat,rfl_sfc_lnd_nrm_mds_lnr_id), &
         & sbr_nm//": def_var rfl_sfc_lnd_nrm_mds_lnr")
    rcd=nf90_wrp(nf90_def_var(nc_id,"rfl_sfc_lnd_nrm_mds_sqr",nf90_float,dmn_lon_lat,rfl_sfc_lnd_nrm_mds_sqr_id), &
         & sbr_nm//": def_var rfl_sfc_lnd_nrm_mds_sqr")
    rcd=nf90_wrp(nf90_def_var(nc_id,"csn_cyc_mdl",nf90_float,dmn_lon_lat_time,csn_cyc_mdl_id), &
         & sbr_nm//": def_var csn_cyc_mdl")
    rcd=nf90_wrp(nf90_def_var(nc_id,"flw_acm_fct",nf90_float,dmn_lon_lat,flw_acm_fct_id), &
         & sbr_nm//": def_var flw_acm_fct")
    rcd=nf90_wrp(nf90_def_var(nc_id,"hgt_sfc_gpm",nf90_float,dmn_lon_lat,hgt_sfc_gpm_id), &
         & sbr_nm//": def_var hgt_sfc_gpm")
    rcd=nf90_wrp(nf90_def_var(nc_id,"hgt_sfc_std_dvn",nf90_float,dmn_lon_lat,hgt_sfc_std_dvn_id), &
         & sbr_nm//": def_var hgt_sfc_std_dvn")
    rcd=nf90_wrp(nf90_def_var(nc_id,"hgt_zpd_lsm",nf90_float,dmn_lon_lat,hgt_zpd_lsm_id), &
         & sbr_nm//": def_var hgt_zpd_lsm")
    rcd=nf90_wrp(nf90_def_var(nc_id,"lai_ttl_clm",nf90_float,dmn_lon_lat_time,lai_ttl_clm_id), &
         & sbr_nm//": def_var lai_ttl_clm")
    rcd=nf90_wrp(nf90_def_var(nc_id,"lnd_frc_clm",nf90_float,dmn_lon_lat,lnd_frc_clm_id), &
         & sbr_nm//": def_var lnd_frc_clm")
    rcd=nf90_wrp(nf90_def_var(nc_id,"lnd_frc_dry",nf90_float,dmn_lon_lat,lnd_frc_dry_id), &
         & sbr_nm//": def_var lnd_frc_dry")
    rcd=nf90_wrp(nf90_def_var(nc_id,"mbl_bsn_fct",nf90_float,dmn_lon_lat,mbl_bsn_fct_id), &
         & sbr_nm//": def_var mbl_bsn_fct")
    rcd=nf90_wrp(nf90_def_var(nc_id,"mss_frc_C_org",nf90_float,dmn_lon_lat,mss_frc_C_org_id), &
         & sbr_nm//": def_var mss_frc_C_org")
    rcd=nf90_wrp(nf90_def_var(nc_id,"mss_frc_CaCO3",nf90_float,dmn_lon_lat,mss_frc_CaCO3_id), &
         & sbr_nm//": def_var mss_frc_CaCO3")
    rcd=nf90_wrp(nf90_def_var(nc_id,"mss_frc_P2O5",nf90_float,dmn_lon_lat,mss_frc_P2O5_id), &
         & sbr_nm//": def_var mss_frc_P2O5")
    rcd=nf90_wrp(nf90_def_var(nc_id,"mss_frc_cly",nf90_float,dmn_lon_lat,mss_frc_cly_id), &
         & sbr_nm//": def_var mss_frc_cly")
    rcd=nf90_wrp(nf90_def_var(nc_id,"mss_frc_slt",nf90_float,dmn_lon_lat,mss_frc_slt_id), &
         & sbr_nm//": def_var mss_frc_slt")
    rcd=nf90_wrp(nf90_def_var(nc_id,"mss_frc_snd",nf90_float,dmn_lon_lat,mss_frc_snd_id), &
         & sbr_nm//": def_var mss_frc_snd")
    rcd=nf90_wrp(nf90_def_var(nc_id,"pft_frc_clm",nf90_float,dmn_lon_lat_sgs,pft_frc_clm_id), &
         & sbr_nm//": def_var pft_frc_clm")
    rcd=nf90_wrp(nf90_def_var(nc_id,"sfc_acm_fct",nf90_float,dmn_lon_lat,sfc_acm_fct_id), &
         & sbr_nm//": def_var sfc_acm_fct")
    rcd=nf90_wrp(nf90_def_var(nc_id,"sfc_frc_bln",nf90_float,dmn_lon_lat_bln,sfc_frc_bln_id), &
         & sbr_nm//": def_var sfc_frc_bln")
    rcd=nf90_wrp(nf90_def_var(nc_id,"src_str_mdl",nf90_float,dmn_lon_lat_time,src_str_mdl_id), &
         & sbr_nm//": def_var src_str_mdl")
    rcd=nf90_wrp(nf90_def_var(nc_id,"vai_ttl_clm",nf90_float,dmn_lon_lat_time,vai_ttl_clm_id), &
         & sbr_nm//": def_var vai_ttl_clm")

    ! Not wrapped
    rcd=nf90_wrp(nf90_def_var(nc_id,"area",nf90_double,dmn_lon_lat,area_id),sbr_nm//": def_var area")
    rcd=nf90_wrp(nf90_def_var(nc_id,"bsn_enm",nf90_int,dmn_lon_lat,bsn_enm_id),sbr_nm//": def_var bsn_enm")
    rcd=nf90_wrp(nf90_def_var(nc_id,"bsn_sz",nf90_float,dmn_bsn,bsn_sz_id),sbr_nm//": def_var bsn_sz")
    rcd=nf90_wrp(nf90_def_var(nc_id,"csn_cyc",nf90_float,dmn_lon_lat_time,csn_cyc_id),sbr_nm//": def_var csn_cyc")
    rcd=nf90_wrp(nf90_def_var(nc_id,"fsh_fct",nf90_float,dmn_lon_lat,fsh_fct_id),sbr_nm//": def_var fsh_fct")
    rcd=nf90_wrp(nf90_def_var(nc_id,"grd_mrd_lng",nf90_double,dmn_lat,grd_mrd_lng_id),sbr_nm//": def_var grd_mrd_lng")
    rcd=nf90_wrp(nf90_def_var(nc_id,"grd_rds_rth",nf90_double,dmn_lat,grd_rds_rth_id),sbr_nm//": def_var grd_rds_rth")
    rcd=nf90_wrp(nf90_def_var(nc_id,"grd_znl_lng",nf90_double,dmn_lat,grd_znl_lng_id),sbr_nm//": def_var grd_znl_lng")
    rcd=nf90_wrp(nf90_def_var(nc_id,"gw",nf90_double,dmn_lat,gw_id),sbr_nm//": def_var gw")
    rcd=nf90_wrp(nf90_def_var(nc_id,"hgt_sfc",nf90_float,dmn_lon_lat,hgt_sfc_id),sbr_nm//": def_var hgt_sfc")
    rcd=nf90_wrp(nf90_def_var(nc_id,"lai_clm",nf90_float,dmn_lon_lat_sgs_time,lai_clm_id),sbr_nm//": def_var lai_clm")
    rcd=nf90_wrp(nf90_def_var(nc_id,"lai_lsm",nf90_float,dmn_lon_lat_time,lai_lsm_id),sbr_nm//": def_var lai_lsm")
    rcd=nf90_wrp(nf90_def_var(nc_id,"lak_frc",nf90_float,dmn_lon_lat,lak_frc_id),sbr_nm//": def_var lak_frc")
    rcd=nf90_wrp(nf90_def_var(nc_id,"lat",nf90_double,dmn_lat,lat_id),sbr_nm//": def_var lat")
    rcd=nf90_wrp(nf90_def_var(nc_id,"lat_fv",nf90_double,dmn_lat,lat_fv_id),sbr_nm//": def_var lat_fv")
    rcd=nf90_wrp(nf90_def_var(nc_id,"lat_grd",nf90_double,dmn_lat_grd,lat_grd_id),sbr_nm//": def_var lat_grd")
    rcd=nf90_wrp(nf90_def_var(nc_id,"lat_sz",nf90_double,dmn_lat,lat_sz_id),sbr_nm//": def_var lat_sz")
    rcd=nf90_wrp(nf90_def_var(nc_id,"lnd_frc",nf90_float,dmn_lon_lat,lnd_frc_id),sbr_nm//": def_var lnd_frc")
    rcd=nf90_wrp(nf90_def_var(nc_id,"lnd_msk",nf90_int,dmn_lon_lat,lnd_msk_id),sbr_nm//": def_var lnd_msk")
    rcd=nf90_wrp(nf90_def_var(nc_id,"lon",nf90_double,dmn_lon,lon_id),sbr_nm//": def_var lon")
    rcd=nf90_wrp(nf90_def_var(nc_id,"lon_grd",nf90_double,dmn_lon_grd,lon_grd_id),sbr_nm//": def_var lon_grd")
    rcd=nf90_wrp(nf90_def_var(nc_id,"lon_sz",nf90_double,dmn_lon,lon_sz_id),sbr_nm//": def_var lon_sz")
    rcd=nf90_wrp(nf90_def_var(nc_id,"msk_rgn",nf90_int,dmn_lon_lat,msk_rgn_id),sbr_nm//": def_var msk_rgn")
    rcd=nf90_wrp(nf90_def_var(nc_id,"mss_frc_Al",nf90_float,dmn_lon_lat,mss_frc_Al_id),sbr_nm//": def_var mss_frc_Al")
    rcd=nf90_wrp(nf90_def_var(nc_id,"mss_frc_K",nf90_float,dmn_lon_lat,mss_frc_K_id),sbr_nm//": def_var mss_frc_K")
    rcd=nf90_wrp(nf90_def_var(nc_id,"mss_frc_N",nf90_float,dmn_lon_lat,mss_frc_N_id),sbr_nm//": def_var mss_frc_N")
    rcd=nf90_wrp(nf90_def_var(nc_id,"mss_frc_Na",nf90_float,dmn_lon_lat,mss_frc_Na_id),sbr_nm//": def_var mss_frc_Na")
    rcd=nf90_wrp(nf90_def_var(nc_id,"mth_max",nf90_int,dmn_lon_lat,mth_max_id),sbr_nm//": def_var mth_max")
    rcd=nf90_wrp(nf90_def_var(nc_id,"odxc_tms",nf90_float,dmn_lon_lat_time,odxc_tms_id),sbr_nm//": def_var odxc_tms")
    rcd=nf90_wrp(nf90_def_var(nc_id,"oro",nf90_float,dmn_lon_lat,oro_id),sbr_nm//": def_var oro")
    rcd=nf90_wrp(nf90_def_var(nc_id,"pft_clm",nf90_int,dmn_lon_lat_sgs,pft_clm_id),sbr_nm//": def_var pft_clm")
    rcd=nf90_wrp(nf90_def_var(nc_id,"rgh_mmn",nf90_float,dmn_lon_lat,rgh_mmn_id),sbr_nm//": def_var rgh_mmn")
    rcd=nf90_wrp(nf90_def_var(nc_id,"sfc_typ",nf90_int,dmn_lon_lat,sfc_typ_id),sbr_nm//": def_var sfc_typ")
    rcd=nf90_wrp(nf90_def_var(nc_id,"slat",nf90_double,dmn_slat,slat_id),sbr_nm//": def_var slat")
    rcd=nf90_wrp(nf90_def_var(nc_id,"slon",nf90_double,dmn_slon,slon_id),sbr_nm//": def_var slon")
    rcd=nf90_wrp(nf90_def_var(nc_id,"smp_sat",nf90_float,dmn_lon_lat,smp_sat_id),sbr_nm//": def_var smp_sat")
    rcd=nf90_wrp(nf90_def_var(nc_id,"smp_sfc",nf90_float,dmn_lon_lat_time,smp_sfc_id),sbr_nm//": def_var smp_sfc")
    rcd=nf90_wrp(nf90_def_var(nc_id,"src_flg",nf90_int,dmn_lon_lat,src_flg_id),sbr_nm//": def_var src_flg")
    rcd=nf90_wrp(nf90_def_var(nc_id,"src_frq",nf90_float,dmn_lon_lat_time,src_frq_id),sbr_nm//": def_var src_frq")
    rcd=nf90_wrp(nf90_def_var(nc_id,"src_odxc",nf90_float,dmn_lon_lat_time,src_odxc_id),sbr_nm//": def_var src_odxc")
    rcd=nf90_wrp(nf90_def_var(nc_id,"src_str",nf90_float,dmn_lon_lat_time,src_str_id),sbr_nm//": def_var src_str")
    rcd=nf90_wrp(nf90_def_var(nc_id,"src_thr",nf90_float,src_thr_id),sbr_nm//": def_var src_thr")
    rcd=nf90_wrp(nf90_def_var(nc_id,"time",nf90_double,dmn_time,time_id),sbr_nm//": def_var time")
    rcd=nf90_wrp(nf90_def_var(nc_id,"vai_clm",nf90_float,dmn_lon_lat_sgs_time,vai_clm_id),sbr_nm//": def_var vai_clm")
    rcd=nf90_wrp(nf90_def_var(nc_id,"vai_lsm",nf90_float,dmn_lon_lat_time,vai_lsm_id),sbr_nm//": def_var vai_lsm")
    rcd=nf90_wrp(nf90_def_var(nc_id,"vwc_dry",nf90_float,dmn_lon_lat,vwc_dry_id),sbr_nm//": def_var vwc_dry")
    rcd=nf90_wrp(nf90_def_var(nc_id,"vwc_opt",nf90_float,dmn_lon_lat,vwc_opt_id),sbr_nm//": def_var vwc_opt")
    rcd=nf90_wrp(nf90_def_var(nc_id,"vwc_rel",nf90_float,dmn_lon_lat_time,vwc_rel_id),sbr_nm//": def_var vwc_rel")
    rcd=nf90_wrp(nf90_def_var(nc_id,"vwc_sat",nf90_float,dmn_lon_lat,vwc_sat_id),sbr_nm//": def_var vwc_sat")
    rcd=nf90_wrp(nf90_def_var(nc_id,"vwc_sfc",nf90_float,dmn_lon_lat_time,vwc_sfc_id),sbr_nm//": def_var vwc_sfc")
    rcd=nf90_wrp(nf90_def_var(nc_id,"wtl_frc",nf90_float,dmn_lon_lat,wtl_frc_id),sbr_nm//": def_var wtl_frc")
    ! Add english text descriptions, long ones first
    rcd=nf90_wrp(nf90_put_att(nc_id,rfl_sfc_lnd_nrm_mds_lnr_id, &
         "long_name","Normalized annual mean land surface reflectance from MODIS, linear"), &
         sbr_nm//': put_att rfl_sfc_lnd_nrm_mds_lnr')
    rcd=nf90_wrp(nf90_put_att(nc_id,rfl_sfc_lnd_nrm_mds_sqr_id, &
         "long_name","Normalized annual mean land surface reflectance from MODIS, squared"), &
         sbr_nm//': put_att rfl_sfc_lnd_nrm_mds_sqr')
    rcd=nf90_wrp(nf90_put_att(nc_id,aod_nmd_frc_id,"long_name","Fraction of TOMS AOD not from dust"), &
         sbr_nm//': put_att aod_nmd_frc')
    rcd=nf90_wrp(nf90_put_att(nc_id,grd_rds_rth_id,"long_name","Distance from axis of rotation at local latitude"), &
         sbr_nm//': put_att grd_rds_rth')
    ! Not wrapped
    rcd=nf90_wrp(nf90_put_att(nc_id,area_id,"long_name","Surface area"),sbr_nm//': put_att area')
    rcd=nf90_wrp(nf90_put_att(nc_id,bsn_enm_id,"long_name","Basin ID"),sbr_nm//': put_att bsn_enm')
    rcd=nf90_wrp(nf90_put_att(nc_id,bsn_sz_id,"long_name","Basin area"),sbr_nm//': put_att bsn_sz')
    rcd=nf90_wrp(nf90_put_att(nc_id,csn_cyc_id,"long_name","Seasonal cycle"),sbr_nm//': put_att csn_cyc')
    rcd=nf90_wrp(nf90_put_att(nc_id,csn_cyc_mdl_id,"long_name","Seasonal cycle, CCM"),sbr_nm//': put_att csn_cyc_mdl')
    rcd=nf90_wrp(nf90_put_att(nc_id,flw_acm_fct_id,"long_name","Flow accumulation factor"),sbr_nm//': put_att flw_acm_fct')
    rcd=nf90_wrp(nf90_put_att(nc_id,fsh_fct_id,"long_name","Efficiency factor"),sbr_nm//': put_att fsh_fct')
    rcd=nf90_wrp(nf90_put_att(nc_id,grd_mrd_lng_id,"long_name","Grid meridional length"),sbr_nm//': put_att grd_mrd_lng')
    rcd=nf90_wrp(nf90_put_att(nc_id,grd_znl_lng_id,"long_name","Grid zonal length"),sbr_nm//': put_att grd_znl_lng')
    rcd=nf90_wrp(nf90_put_att(nc_id,gw_id,"long_name","Gaussian weight"),sbr_nm//': put_att gw')
    rcd=nf90_wrp(nf90_put_att(nc_id,hgt_sfc_gpm_id,"long_name","Surface geopotential height"),sbr_nm//': put_att hgt_sfc_gpm')
    rcd=nf90_wrp(nf90_put_att(nc_id,sfc_acm_fct_id,"long_name","Surface area accumulation factor"),sbr_nm//': put_att sfc_acm_fct')
    rcd=rcd+nf90_put_att(nc_id,hgt_sfc_id,"long_name","Surface height")
    rcd=rcd+nf90_put_att(nc_id,hgt_sfc_std_dvn_id,"long_name","Standard deviation of surface height")
    rcd=rcd+nf90_put_att(nc_id,hgt_zpd_lsm_id,"long_name","Zero plane displacement height")
    rcd=rcd+nf90_put_att(nc_id,lai_clm_id,"long_name","Leaf area index, one-sided")
    rcd=rcd+nf90_put_att(nc_id,lai_lsm_id,"long_name","Leaf area index, one-sided")
    rcd=rcd+nf90_put_att(nc_id,lai_ttl_clm_id,"long_name","Total leaf area index, one-sided")
    rcd=rcd+nf90_put_att(nc_id,lak_frc_id,"long_name","Lake Fraction")
    rcd=rcd+nf90_put_att(nc_id,lat_fv_id,"long_name","Latitude grid")
    rcd=rcd+nf90_put_att(nc_id,lat_grd_id,"long_name","Latitude grid")
    rcd=rcd+nf90_put_att(nc_id,lat_id,"long_name","Latitude")
    rcd=rcd+nf90_put_att(nc_id,lat_sz_id,"long_name","Latitudinal size")
    rcd=rcd+nf90_put_att(nc_id,lnd_frc_clm_id,"long_name","Fraction of land (not ocean)")
    rcd=rcd+nf90_put_att(nc_id,lnd_frc_dry_id,"long_name","Dry Land Fraction")
    rcd=rcd+nf90_put_att(nc_id,lnd_frc_id,"long_name","Land Fraction")
    rcd=rcd+nf90_put_att(nc_id,lnd_msk_id,"long_name","Land Mask")
    rcd=rcd+nf90_put_att(nc_id,lon_grd_id,"long_name","Longitude grid")
    rcd=rcd+nf90_put_att(nc_id,lon_id,"long_name","Longitude")
    rcd=rcd+nf90_put_att(nc_id,lon_sz_id,"long_name","Longitudinal size")
    rcd=rcd+nf90_put_att(nc_id,mbl_bsn_fct_id,"long_name","Mobilization enhancement due to basin characteristics")
    rcd=rcd+nf90_put_att(nc_id,msk_rgn_id,"long_name","Mask region")
    rcd=rcd+nf90_put_att(nc_id,mss_frc_Al_id,"long_name","Exchangeable aluminum")
    rcd=rcd+nf90_put_att(nc_id,mss_frc_C_org_id,"long_name","Organic carbon")
    rcd=rcd+nf90_put_att(nc_id,mss_frc_CaCO3_id,"long_name","Calcium carbonate")
    rcd=rcd+nf90_put_att(nc_id,mss_frc_K_id,"long_name","Exchangeable potassium")
    rcd=rcd+nf90_put_att(nc_id,mss_frc_N_id,"long_name","Nitrogen")
    rcd=rcd+nf90_put_att(nc_id,mss_frc_Na_id,"long_name","Exchangeable sodium")
    rcd=rcd+nf90_put_att(nc_id,mss_frc_P2O5_id,"long_name","Extractable phosphorus")
    rcd=rcd+nf90_put_att(nc_id,mss_frc_cly_id,"long_name","Clay texture fraction")
    rcd=rcd+nf90_put_att(nc_id,mss_frc_slt_id,"long_name","Clay texture fraction")
    rcd=rcd+nf90_put_att(nc_id,mss_frc_snd_id,"long_name","Clay texture fraction")
    rcd=rcd+nf90_put_att(nc_id,mth_max_id,"long_name","Peak month of annual cycle")
    rcd=rcd+nf90_put_att(nc_id,odxc_tms_id,"long_name","Extinction optical depth proxy from TOMS")
    rcd=rcd+nf90_put_att(nc_id,oro_id,"long_name","Orography")
    rcd=rcd+nf90_put_att(nc_id,pft_clm_id,"long_name","Plant functional type")
    rcd=rcd+nf90_put_att(nc_id,pft_frc_clm_id,"long_name","Fraction covered by plant functional type")
    rcd=rcd+nf90_put_att(nc_id,rgh_mmn_id,"long_name","Roughness length for humidity")
    rcd=rcd+nf90_put_att(nc_id,sfc_frc_bln_id,"long_name","Surface area fraction of soil blend")
    rcd=rcd+nf90_put_att(nc_id,sfc_typ_id,"long_name","Surface type code")
    rcd=rcd+nf90_put_att(nc_id,slat_id,"long_name","Staggered latitude for FV grid")
    rcd=rcd+nf90_put_att(nc_id,slon_id,"long_name","Staggered longitude for FV grid")
    rcd=rcd+nf90_put_att(nc_id,smp_sat_id,"long_name","Saturated soil matric potential (sand-dependent)")
    rcd=rcd+nf90_put_att(nc_id,smp_sfc_id,"long_name","Soil matric potential")
    rcd=rcd+nf90_put_att(nc_id,src_flg_id,"long_name","Source flag")
    rcd=rcd+nf90_put_att(nc_id,src_frq_id,"long_name","Source frequency")
    rcd=rcd+nf90_put_att(nc_id,src_odxc_id,"long_name","Extinction optical depth, source regions only")
    rcd=rcd+nf90_put_att(nc_id,src_str_id,"long_name","Source strength")
    rcd=rcd+nf90_put_att(nc_id,src_str_mdl_id,"long_name","Source strength, CCM")
    rcd=rcd+nf90_put_att(nc_id,src_thr_id,"long_name","Source strength threshold")
    rcd=rcd+nf90_put_att(nc_id,time_id,"long_name","Day of Year")
    rcd=rcd+nf90_put_att(nc_id,vai_clm_id,"long_name","Vegetation area index, one-sided")
    rcd=rcd+nf90_put_att(nc_id,vai_lsm_id,"long_name","Vegetation area index, one-sided")
    rcd=rcd+nf90_put_att(nc_id,vai_ttl_clm_id,"long_name","Total vegetation area index, one-sided")
    rcd=rcd+nf90_put_att(nc_id,vwc_dry_id,"long_name","Dry volumetric water content (no E-T)")
    rcd=rcd+nf90_put_att(nc_id,vwc_opt_id,"long_name","E-T optimal volumetric water content")
    rcd=rcd+nf90_put_att(nc_id,vwc_rel_id,"long_name","Water content relative to saturation")
    rcd=rcd+nf90_put_att(nc_id,vwc_sat_id,"long_name","Saturated volumetric water content (sand-dependent)")
    rcd=rcd+nf90_put_att(nc_id,vwc_sfc_id,"long_name","Volumetric water content")
    rcd=rcd+nf90_put_att(nc_id,wtl_frc_id,"long_name","Wetland Fraction")
    ! Add units
    rcd=rcd+nf90_put_att(nc_id,aod_nmd_frc_id,"units","fraction")
    rcd=rcd+nf90_put_att(nc_id,area_id,"units","meter2")
    rcd=rcd+nf90_put_att(nc_id,bsn_enm_id,"units","enumerated")
    rcd=rcd+nf90_put_att(nc_id,bsn_sz_id,"units","meter2")
    rcd=rcd+nf90_put_att(nc_id,csn_cyc_id,"units","fraction")
    rcd=rcd+nf90_put_att(nc_id,csn_cyc_mdl_id,"units","fraction")
    rcd=rcd+nf90_put_att(nc_id,flw_acm_fct_id,"units","fraction")
    rcd=rcd+nf90_put_att(nc_id,fsh_fct_id,"units","unknown")
    rcd=rcd+nf90_put_att(nc_id,grd_mrd_lng_id,"units","meter")
    rcd=rcd+nf90_put_att(nc_id,grd_rds_rth_id,"units","meter")
    rcd=rcd+nf90_put_att(nc_id,grd_znl_lng_id,"units","meter")
    rcd=rcd+nf90_put_att(nc_id,gw_id,"units","fraction")
    rcd=rcd+nf90_put_att(nc_id,hgt_sfc_gpm_id,"units","geopotentialmeter")
    rcd=rcd+nf90_put_att(nc_id,hgt_sfc_id,"units","meter")
    rcd=rcd+nf90_put_att(nc_id,hgt_sfc_std_dvn_id,"units","meter")
    rcd=rcd+nf90_put_att(nc_id,hgt_zpd_lsm_id,"units","meter")
    rcd=rcd+nf90_put_att(nc_id,lai_clm_id,"units","meter2 meter-2")
    rcd=rcd+nf90_put_att(nc_id,lai_lsm_id,"units","meter2 meter-2")
    rcd=rcd+nf90_put_att(nc_id,lai_ttl_clm_id,"units","meter2 meter-2")
    rcd=rcd+nf90_put_att(nc_id,lak_frc_id,"units","fraction")
    rcd=rcd+nf90_put_att(nc_id,lat_fv_id,"units","degrees_north")
    rcd=rcd+nf90_put_att(nc_id,lat_grd_id,"units","degrees_north")
    rcd=rcd+nf90_put_att(nc_id,lat_id,"units","degrees_north")
    rcd=rcd+nf90_put_att(nc_id,lat_sz_id,"units","degree")
    rcd=rcd+nf90_put_att(nc_id,lnd_frc_clm_id,"units","fraction")
    rcd=rcd+nf90_put_att(nc_id,lnd_frc_dry_id,"units","fraction")
    rcd=rcd+nf90_put_att(nc_id,lnd_frc_id,"units","fraction")
    rcd=rcd+nf90_put_att(nc_id,lnd_msk_id,"units","index")
    rcd=rcd+nf90_put_att(nc_id,lon_grd_id,"units","degrees_east")
    rcd=rcd+nf90_put_att(nc_id,lon_id,"units","degrees_east")
    rcd=rcd+nf90_put_att(nc_id,lon_sz_id,"units","degree")
    rcd=rcd+nf90_put_att(nc_id,mbl_bsn_fct_id,"units","fraction")
    rcd=rcd+nf90_put_att(nc_id,msk_rgn_id,"units","bool")
    rcd=rcd+nf90_put_att(nc_id,mss_frc_Al_id,"units","kilogram kilogram-1")
    rcd=rcd+nf90_put_att(nc_id,mss_frc_C_org_id,"units","kilogram kilogram-1")
    rcd=rcd+nf90_put_att(nc_id,mss_frc_CaCO3_id,"units","kilogram kilogram-1")
    rcd=rcd+nf90_put_att(nc_id,mss_frc_K_id,"units","kilogram kilogram-1")
    rcd=rcd+nf90_put_att(nc_id,mss_frc_N_id,"units","kilogram kilogram-1")
    rcd=rcd+nf90_put_att(nc_id,mss_frc_Na_id,"units","kilogram kilogram-1")
    rcd=rcd+nf90_put_att(nc_id,mss_frc_P2O5_id,"units","kilogram kilogram-1")
    rcd=rcd+nf90_put_att(nc_id,mss_frc_cly_id,"units","kilogram kilogram-1")
    rcd=rcd+nf90_put_att(nc_id,mss_frc_slt_id,"units","kilogram kilogram-1")
    rcd=rcd+nf90_put_att(nc_id,mss_frc_snd_id,"units","kilogram kilogram-1")
    rcd=rcd+nf90_put_att(nc_id,mth_max_id,"units","month")
    rcd=rcd+nf90_put_att(nc_id,odxc_tms_id,"units","fraction")
    rcd=rcd+nf90_put_att(nc_id,oro_id,"units","fraction")
    rcd=rcd+nf90_put_att(nc_id,pft_clm_id,"units","enumerated")
    rcd=rcd+nf90_put_att(nc_id,pft_frc_clm_id,"units","fraction")
    rcd=rcd+nf90_put_att(nc_id,rfl_sfc_lnd_nrm_mds_lnr_id,"units","fraction")
    rcd=rcd+nf90_put_att(nc_id,rfl_sfc_lnd_nrm_mds_sqr_id,"units","fraction")
    rcd=rcd+nf90_put_att(nc_id,rgh_mmn_id,"units","meter")
    rcd=rcd+nf90_put_att(nc_id,sfc_acm_fct_id,"units","fraction")
    rcd=rcd+nf90_put_att(nc_id,sfc_frc_bln_id,"units","meter2 meter-2")
    rcd=rcd+nf90_put_att(nc_id,sfc_typ_id,"units","index")
    rcd=rcd+nf90_put_att(nc_id,slat_id,"units","degrees_north")
    rcd=rcd+nf90_put_att(nc_id,slon_id,"units","degrees_east")
    rcd=rcd+nf90_put_att(nc_id,smp_sat_id,"units","millimeter H2O")
    rcd=rcd+nf90_put_att(nc_id,smp_sfc_id,"units","millimeter H2O")
    rcd=rcd+nf90_put_att(nc_id,src_flg_id,"units","bool")
    rcd=rcd+nf90_put_att(nc_id,src_frq_id,"units","fraction")
    rcd=rcd+nf90_put_att(nc_id,src_odxc_id,"units","fraction")
    rcd=rcd+nf90_put_att(nc_id,src_str_id,"units","fraction")
    rcd=rcd+nf90_put_att(nc_id,src_str_mdl_id,"units","fraction")
    rcd=rcd+nf90_put_att(nc_id,src_thr_id,"units","fraction")
    rcd=rcd+nf90_put_att(nc_id,time_id,"units","day")
    rcd=rcd+nf90_put_att(nc_id,vai_clm_id,"units","meter2 meter-2")
    rcd=rcd+nf90_put_att(nc_id,vai_lsm_id,"units","meter2 meter-2")
    rcd=rcd+nf90_put_att(nc_id,vai_ttl_clm_id,"units","meter2 meter-2")
    rcd=rcd+nf90_put_att(nc_id,vwc_dry_id,"units","meter3 meter-3")
    rcd=rcd+nf90_put_att(nc_id,vwc_opt_id,"units","meter3 meter-3")
    rcd=rcd+nf90_put_att(nc_id,vwc_rel_id,"units","fraction")
    rcd=rcd+nf90_put_att(nc_id,vwc_sat_id,"units","meter3 meter-3")
    rcd=rcd+nf90_put_att(nc_id,vwc_sfc_id,"units","meter3 meter-3")
    rcd=rcd+nf90_put_att(nc_id,wtl_frc_id,"units","fraction")

    ! Add missing_value and _FillValue to selected fields
    rcd=rcd+nf90_put_att(nc_id,mss_frc_Al_id,"missing_value",mss_val)
    rcd=rcd+nf90_put_att(nc_id,mss_frc_C_org_id,"missing_value",mss_val)
    rcd=rcd+nf90_put_att(nc_id,mss_frc_CaCO3_id,"missing_value",mss_val)
    rcd=rcd+nf90_put_att(nc_id,mss_frc_K_id,"missing_value",mss_val)
    rcd=rcd+nf90_put_att(nc_id,mss_frc_N_id,"missing_value",mss_val)
    rcd=rcd+nf90_put_att(nc_id,mss_frc_Na_id,"missing_value",mss_val)
    rcd=rcd+nf90_put_att(nc_id,mss_frc_P2O5_id,"missing_value",mss_val)
    rcd=rcd+nf90_put_att(nc_id,mss_frc_cly_id,"missing_value",mss_val)
    rcd=rcd+nf90_put_att(nc_id,mss_frc_slt_id,"missing_value",mss_val)
    rcd=rcd+nf90_put_att(nc_id,mss_frc_snd_id,"missing_value",mss_val)

    rcd=rcd+nf90_put_att(nc_id,mss_frc_Al_id,"_FillValue",mss_val)
    rcd=rcd+nf90_put_att(nc_id,mss_frc_C_org_id,"_FillValue",mss_val)
    rcd=rcd+nf90_put_att(nc_id,mss_frc_CaCO3_id,"_FillValue",mss_val)
    rcd=rcd+nf90_put_att(nc_id,mss_frc_K_id,"_FillValue",mss_val)
    rcd=rcd+nf90_put_att(nc_id,mss_frc_N_id,"_FillValue",mss_val)
    rcd=rcd+nf90_put_att(nc_id,mss_frc_Na_id,"_FillValue",mss_val)
    rcd=rcd+nf90_put_att(nc_id,mss_frc_P2O5_id,"_FillValue",mss_val)
    rcd=rcd+nf90_put_att(nc_id,mss_frc_cly_id,"_FillValue",mss_val)
    rcd=rcd+nf90_put_att(nc_id,mss_frc_slt_id,"_FillValue",mss_val)
    rcd=rcd+nf90_put_att(nc_id,mss_frc_snd_id,"_FillValue",mss_val)
    ! Other field-specific attributes
    rcd=rcd+nf90_put_att(nc_id,sfc_frc_bln_id,"key", &
         "1 = ASS = Aluminosilicated silt, 2 = FS Fine sand, 3 = SS = Salts, 4 = CS = Coarse sand")
    ! Now that all dimensions, variables, and attributes have been defined, make call to end define mode
    rcd=rcd+nf90_enddef(nc_id)
    ! Write data
    rcd=rcd+nf90_put_var(nc_id,time_id,time_out,start=(/1/),count=(/time_out_nbr/)) ! Write record coordinate first to set record dimension size for rest of record variables
    rcd=rcd+nf90_put_var(nc_id,aod_nmd_frc_id,aod_nmd_frc)
    rcd=rcd+nf90_put_var(nc_id,area_id,area)
    rcd=rcd+nf90_put_var(nc_id,bsn_enm_id,bsn_enm)
    rcd=rcd+nf90_put_var(nc_id,bsn_sz_id,bsn_sz)
    rcd=rcd+nf90_put_var(nc_id,csn_cyc_id,csn_cyc)
    rcd=rcd+nf90_put_var(nc_id,csn_cyc_mdl_id,csn_cyc_mdl)
    rcd=rcd+nf90_put_var(nc_id,flw_acm_fct_id,flw_acm_fct)
    rcd=rcd+nf90_put_var(nc_id,fsh_fct_id,fsh_fct)
    rcd=rcd+nf90_put_var(nc_id,grd_mrd_lng_id,grd_mrd_lng)
    rcd=rcd+nf90_put_var(nc_id,grd_rds_rth_id,grd_rds_rth)
    rcd=rcd+nf90_put_var(nc_id,grd_znl_lng_id,grd_znl_lng)
    rcd=rcd+nf90_put_var(nc_id,gw_id,gw)
    rcd=rcd+nf90_put_var(nc_id,hgt_sfc_gpm_id,hgt_sfc_gpm)
    rcd=rcd+nf90_put_var(nc_id,hgt_sfc_id,hgt_sfc)
    rcd=rcd+nf90_put_var(nc_id,hgt_sfc_std_dvn_id,hgt_sfc_std_dvn)
    rcd=rcd+nf90_put_var(nc_id,hgt_zpd_lsm_id,hgt_zpd_lsm)
    rcd=rcd+nf90_put_var(nc_id,lai_clm_id,lai_clm)
    rcd=rcd+nf90_put_var(nc_id,lai_lsm_id,lai_lsm)
    rcd=rcd+nf90_put_var(nc_id,lai_ttl_clm_id,lai_ttl_clm)
    rcd=rcd+nf90_put_var(nc_id,lak_frc_id,lak_frc)
    rcd=rcd+nf90_put_var(nc_id,lat_fv_id,lat_fv)
    rcd=rcd+nf90_put_var(nc_id,lat_grd_id,lat_grd)
    rcd=rcd+nf90_put_var(nc_id,lat_id,lat)
    rcd=rcd+nf90_put_var(nc_id,lat_sz_id,lat_sz)
    rcd=rcd+nf90_put_var(nc_id,lnd_frc_clm_id,lnd_frc_clm)
    rcd=rcd+nf90_put_var(nc_id,lnd_frc_dry_id,lnd_frc_dry)
    rcd=rcd+nf90_put_var(nc_id,lnd_frc_id,lnd_frc)
    rcd=rcd+nf90_put_var(nc_id,lnd_msk_id,lnd_msk)
    rcd=rcd+nf90_put_var(nc_id,lon_grd_id,lon_grd)
    rcd=rcd+nf90_put_var(nc_id,lon_id,lon)
    rcd=rcd+nf90_put_var(nc_id,lon_sz_id,lon_sz)
    rcd=rcd+nf90_put_var(nc_id,mbl_bsn_fct_id,mbl_bsn_fct)
    rcd=rcd+nf90_put_var(nc_id,msk_rgn_id,msk_rgn)
    rcd=rcd+nf90_put_var(nc_id,mss_frc_Al_id,mss_frc_Al)
    rcd=rcd+nf90_put_var(nc_id,mss_frc_C_org_id,mss_frc_C_org)
    rcd=rcd+nf90_put_var(nc_id,mss_frc_CaCO3_id,mss_frc_CaCO3)
    rcd=rcd+nf90_put_var(nc_id,mss_frc_K_id,mss_frc_K)
    rcd=rcd+nf90_put_var(nc_id,mss_frc_N_id,mss_frc_N)
    rcd=rcd+nf90_put_var(nc_id,mss_frc_Na_id,mss_frc_Na)
    rcd=rcd+nf90_put_var(nc_id,mss_frc_P2O5_id,mss_frc_P2O5)
    rcd=rcd+nf90_put_var(nc_id,mss_frc_cly_id,mss_frc_cly)
    rcd=rcd+nf90_put_var(nc_id,mss_frc_slt_id,mss_frc_slt)
    rcd=rcd+nf90_put_var(nc_id,mss_frc_snd_id,mss_frc_snd)
    rcd=rcd+nf90_put_var(nc_id,mth_max_id,mth_max)
    rcd=rcd+nf90_put_var(nc_id,odxc_tms_id,odxc_tms)
    rcd=rcd+nf90_put_var(nc_id,oro_id,oro)
    rcd=rcd+nf90_put_var(nc_id,pft_clm_id,pft_clm)
    rcd=rcd+nf90_put_var(nc_id,pft_frc_clm_id,pft_frc_clm)
    rcd=rcd+nf90_put_var(nc_id,rfl_sfc_lnd_nrm_mds_lnr_id,rfl_sfc_lnd_nrm_mds_lnr)
    rcd=rcd+nf90_put_var(nc_id,rfl_sfc_lnd_nrm_mds_sqr_id,rfl_sfc_lnd_nrm_mds_sqr)
    rcd=rcd+nf90_put_var(nc_id,rgh_mmn_id,rgh_mmn)
    rcd=rcd+nf90_put_var(nc_id,sfc_acm_fct_id,sfc_acm_fct)
    rcd=rcd+nf90_put_var(nc_id,sfc_frc_bln_id,sfc_frc_bln)
    rcd=rcd+nf90_put_var(nc_id,sfc_typ_id,sfc_typ)
    rcd=rcd+nf90_put_var(nc_id,slat_id,slat)
    rcd=rcd+nf90_put_var(nc_id,slon_id,slon)
    rcd=rcd+nf90_put_var(nc_id,smp_sat_id,smp_sat)
    rcd=rcd+nf90_put_var(nc_id,smp_sfc_id,smp_sfc)
    rcd=rcd+nf90_put_var(nc_id,src_flg_id,src_flg)
    rcd=rcd+nf90_put_var(nc_id,src_frq_id,src_frq)
    rcd=rcd+nf90_put_var(nc_id,src_odxc_id,src_odxc)
    rcd=rcd+nf90_put_var(nc_id,src_str_id,src_str)
    rcd=rcd+nf90_put_var(nc_id,src_str_mdl_id,src_str_mdl)
    rcd=rcd+nf90_put_var(nc_id,src_thr_id,src_thr)
    rcd=rcd+nf90_put_var(nc_id,vai_clm_id,vai_clm)
    rcd=rcd+nf90_put_var(nc_id,vai_lsm_id,vai_lsm)
    rcd=rcd+nf90_put_var(nc_id,vai_ttl_clm_id,vai_ttl_clm)
    rcd=rcd+nf90_put_var(nc_id,vwc_dry_id,vwc_dry)
    rcd=rcd+nf90_put_var(nc_id,vwc_opt_id,vwc_opt)
    rcd=rcd+nf90_put_var(nc_id,vwc_rel_id,vwc_rel)
    rcd=rcd+nf90_put_var(nc_id,vwc_sat_id,vwc_sat)
    rcd=rcd+nf90_put_var(nc_id,vwc_sfc_id,vwc_sfc)
    rcd=rcd+nf90_put_var(nc_id,wtl_frc_id,wtl_frc)

    ! Close output file
    rcd=nf90_wrp_close(nc_id,fl_out,"Wrote results to",sbr_nm=sbr_nm)
    
    if (dbg_lvl >= dbg_sbr) write(6,"(a)") "Exiting bds_out()"
    return 
  end subroutine bds_out                       ! end bds_out()
  
end module bds_drv ! [mdl] BDS driver, processor, output routines
