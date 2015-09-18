! $Id$ -*-f90-*-

! Purpose: Soil texture and hydrology

! Usage: 
! use soi_txt ! [mdl] Soil texture and hydrology

module soi_txt ! [mdl] Soil texture and hydrology
  implicit none
  private ! [stt] Symbols are private unless individually qualified as public
  public::soi_hyd_prc ! [sbr] Hydrologic properties of soil texture  
  public::soi_txt_get ! [sbr] Create soil texture field

contains
  
  subroutine soi_txt_get( & ! 
       fl_in, & ! I
       lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
       lat_out_grd,lat_out_nbr,lon_out_grd,lon_out_nbr, & ! I
       ovr_nbr_max,area_out,lnd_msk_out,sfc_typ_LSM_out,flg_soi_nn, & ! I
       mss_frc_cly_out,mss_frc_slt_out,mss_frc_snd_out) ! O
    ! Purpose: Retrieve soil texture and rebin to requested grid
    ! soi_txt_get() is called by lnd_bnd_cnv_rbn() in bds_drv.F90
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use map_grd ! [mdl] Map grids and regridding
    use netcdf ! [mdl] netCDF interface
    use nf90_utl ! [mdl] netCDF utilities
    use sng_mdl,only:ftn_strlen,ftn_strstr ! [mdl] String manipulation
    use utl_mdl,only:mnt_dcr_chk,mnt_ncr_chk ! [mdl] Utility functions (date_time_get,mnt_chk...)
    implicit none
    ! Parameters
    character(len=*),parameter::sbr_nm="soi_txt_get" ! [sng] Subroutine name
    integer,parameter::fl_in_unit=73 ! Unit for reading input data
    integer,parameter::soil_txt_LSM_max=4 ! Maximum value of LSM soil textures
    integer,parameter::soil_txt_LSM_nbr=5 ! [nbr] Number of LSM soil textures
    integer,parameter::unit_dgn=74 ! Unit for writing general Map diagnoses
    real,parameter::mss_frc_snd_loam=0.43 ! [frc] Soil texture sand for loam (Bon96)
    real,parameter::mss_frc_slt_loam=0.39 ! [frc] Soil texture silt for loam (Bon96)
    real,parameter::mss_frc_cly_loam=0.18 ! [frc] Soil texture clay for loam (Bon96)
    ! Commons
    ! Input
    character(len=*),intent(in)::fl_in        ! [sng] Input file
    integer,intent(in)::lat_in_nbr        ! [nbr] Number of latitudes
    integer,intent(in)::lat_out_nbr       ! [nbr] Number of latitudes
    integer,intent(in)::lon_in_nbr        ! [nbr] Number of longitudes
    integer,intent(in)::lon_out_nbr       ! [nbr] Number of longitudes
    integer,intent(in)::ovr_nbr_max       ! [nbr] Maximum number of input cells which overlap any output cell
    integer,intent(in)::lnd_msk_out(lon_out_nbr,lat_out_nbr) ! [msk] Land mask (integer 0 or 1)
    integer,intent(in)::sfc_typ_LSM_out(lon_out_nbr,lat_out_nbr) ! [enm] Surface type 
    real,intent(in)::area_out(lon_out_nbr,lat_out_nbr) ! [m2] Area of gridcells
    real,intent(in)::lat_in_grd(lat_in_nbr+1) ! [dgr] Interface latitudes
    real,intent(in)::lat_out_grd(lat_out_nbr+1) ! [dgr] Interface latitudes
    real,intent(in)::lon_in_grd(lon_in_nbr+1) ! [dgr] Interface longitudes
    real,intent(in)::lon_out_grd(lon_out_nbr+1) ! [dgr] Interface longitudes
    logical,intent(in)::flg_soi_nn
    ! Output
    real,intent(out)::mss_frc_cly_out(lon_out_nbr,lat_out_nbr) ! [frc] Soil texture clay
    real,intent(out)::mss_frc_slt_out(lon_out_nbr,lat_out_nbr) ! [frc] Soil texture silt
    real,intent(out)::mss_frc_snd_out(lon_out_nbr,lat_out_nbr) ! [frc] Soil texture sand
    ! Local
    character(80)::fl_dgn_soi_txt ! Diagnostics file
    integer err_nbr           ! [nbr] Number of errors processing surface types
    integer idx               ! [idx] Counting index
    integer lat_in_idx        ! [idx] Counting index for lat
    integer lat_in_idx_err    ! Latitude of gridcell which caused error
    integer lat_out_idx       ! [idx] Counting index for lat
    integer loam_nbr          ! [nbr] Number of points set to loam texture
    integer lon_in_idx        ! [idx] Counting index for lon
    integer lon_in_idx_err    ! Longitude of gridcell which caused error
    integer lon_out_idx       ! [idx] Counting index for lon
    integer mss_val_ntg       ! [enm] Code for missing/invalid data
    integer ovr_idx           ! [idx] Counting index
    integer ovr_lat_idx(lon_out_nbr,lat_out_nbr,ovr_nbr_max) ! [idx] Map into input grid of latitude indices of overlap cells
    integer ovr_lon_idx(lon_out_nbr,lat_out_nbr,ovr_nbr_max) ! [idx] Map into input grid of longitude indices of overlap cells
    integer ovr_nbr(lon_out_nbr,lat_out_nbr) ! [nbr] Number of input gridcells which overlap each output gridcell
    integer ovr_nbr_crr       ! [nbr] Current number of overlapping gridcells
    integer rcd               ! [enm] Return success code
    integer rnk_1st_idx       ! [idx] Index of most overlapping cell
    integer lnd_msk_in(lon_in_nbr,lat_in_nbr) ! [msk] Land mask 
    logical mnt_ncr           ! [flg] Monotonic and increasing flag
    real area_in(lon_in_nbr,lat_in_nbr) ! [m2] Area of gridcells
    real eps_rlt              ! [frc] Relative error allowed in frc_ttl
    real frc_ttl              ! [frc] Total fraction of gridcell accounted for by all soil textures
    real ovr_wgt(lon_out_nbr,lat_out_nbr,ovr_nbr_max) ! [frc] Weight of overlapping input gridcells onto each output gridcell
    real ovr_wgt_crr          ! [frc] Overlap weight of current gridcell
    real mss_frc_cly_in(lon_in_nbr,lat_in_nbr) ! [frc] Soil texture clay
    real mss_frc_slt_in(lon_in_nbr,lat_in_nbr) ! [frc] Soil texture silt
    real mss_frc_snd_in(lon_in_nbr,lat_in_nbr) ! [frc] Soil texture sand
    character(35) txt_sng(0:soil_txt_LSM_max) ! [sng] String name for each soil type

    ! Local, MF
    integer lat_idx2         ! [idx] Counting index for nearest neighbor loop
    integer lon_idx2         ! [idx] Counting index for nearest neighbor loop
    integer lat_idx_min_dst  ! [idx] Lat index of nearest valid land gridcell for NN option
    integer lon_idx_min_dst  ! [idx] Lon index of nearest valid land gridcell for NN option
    real distance_min        ! [km] Distance of closest valid gridcell for NN option
    real lat1_tmp            ! [nbr] Latitude1 used in NN distance calculation
    real lat2_tmp            ! [nbr] Latitude2 used in NN distance calculation
    real lon1_tmp            ! [nbr] Longitude1 used in NN distance calculation
    real lon2_tmp            ! [nbr] Longitude2 used in NN distance calculation
    real distance_tmp        ! [km] Distance between current 2 gridcells being evaluated in NN loop
    real,parameter::degtorad = 57.2958  ! [nbr] conversion factor for degrees to radians
    integer temp_idx         ! [idx] Counting index for soi_txt text output file
    
    character(80)::fl_soi_txt_out="/data/mflanner/data/map_IGBP.txt" ! [fl] IGBP-style text output file (option)


    save txt_sng
    data txt_sng( 0)/"no soil: ocean, glacier            "/
    data txt_sng( 1)/"clays                              "/
    data txt_sng( 2)/"sands                              "/
    data txt_sng( 3)/"loams                              "/
    data txt_sng( 4)/"silts                              "/
    ! Variables only needed for IBIS dataset
    integer cly_pct_ntg       ! [%] Clay percentage
    integer lat_dmn_id        ! [id] Dimension ID for lat
    integer lat_nbr           ! [nbr] Number of latitudes in input file
    integer lon_dmn_id        ! [id] Dimension ID for lon
    integer lon_nbr           ! [nbr] Number of longitudes in input file
    integer mss_val_id        ! [id] Attribute ID
    integer nc_id             ! [id] File handle
    integer snd_pct_ntg       ! [%] Sand percentage
    integer slt_pct_ntg       ! [%] Silt percentage
    integer txt_IBIS_dmn_nbr  ! [nbr] Number of dimensions in disk txt_IBIS field
    integer lat_in_id         ! [id] Variable ID
    integer txt_IBIS_id       ! [id] Variable ID
    integer txt_pct_ntg       ! [%] Encoded texture percentages
    logical mss_flg           ! Variable has missing_value attribute
    logical flg_IBIS          ! Input file is IBIS 
    logical flg_Webb          ! Input file is Webb
    logical lat_mnt_dcr       ! IBIS latitude monotonically decreases
    real mss_val              ! Missing value
    real lat_in(lat_in_nbr) ! Midpoint latitudes
    real txt_IBIS_in(lon_in_nbr,lat_in_nbr) ! [frc] Soil texture 
    ! Debugging
    integer lon_idx_dbg       ! Longitude for verbose output
    integer lat_idx_dbg       ! Longitude for verbose output
    
    ! Main Code
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Entering soi_txt_get()"
    
    ! Check for monotonically increasing grids
    mnt_ncr=mnt_ncr_chk(lon_out_grd,lon_out_nbr+1)
    if (.not.mnt_ncr) stop "lon_out_grd not monotonically increasing in soi_txt_get()"
    mnt_ncr=mnt_ncr_chk(lat_out_grd,lat_out_nbr+1)
    if (.not.mnt_ncr) stop "lat_out_grd not monotonically increasing in soi_txt_get()"
    
    ! Initialize scalars
    eps_rlt=1.0e-5            ! [frc] Relative error allowed in frc_ttl
    flg_IBIS=.false.
    loam_nbr=0                ! [nbr] Number of points set to loam texture
    mss_val_ntg=-99999        ! 
    lon_idx_dbg=90
    lat_idx_dbg=50
    
    mss_frc_cly_out=0.0       ! [frc]
    mss_frc_slt_out=0.0       ! [frc]
    mss_frc_snd_out=0.0       ! [frc]
    
    ! Compute gridcell area 
    call map_area_get(lat_in_grd,lat_in_nbr,lon_in_grd,lon_in_nbr, & ! I
         area_in) ! O
    
    ! Determine whether input data is IBIS netCDF format or Webb text format
    if (ftn_strstr(fl_in,"IBIS") >= 1) then
       flg_IBIS=.true.
    else                      ! not IBIS
       flg_IBIS=.false.
    endif                     ! endif
    flg_Webb=.not.flg_IBIS
    if (flg_IBIS) then
       fl_dgn_soi_txt="/tmp/zender/map/soi_txt_IGBP_dgn.1x1"//char(0) ! Diagnostics file
    else                      ! not IGBP
       fl_dgn_soi_txt="/tmp/zender/map/soi_txt_webb_dgn.1x1"//char(0) ! Diagnostics file
    endif                     ! not IGBP
    
    if (flg_Webb) then
       ! Webb"s data is stored as a three column formatted text file
       ! An earlier version of these data is the texture dataset of Zobler (1986?), which is part of NCAR SCD DSS DS770.0
       ! The three soil percentages for each (latitude,longitude) point are written on one row
       ! The dataset resolution is 1 x 1 degree and the data are written from the South pole to the North pole, East to West, beginning at the date line
       ! Read strategy is read the entire file in a looped read statement
       open (fl_in_unit,file=fl_in,status="old",iostat=rcd)
       if (rcd /= 0) write (6,"(4a,i4)") prg_nm(1:ftn_strlen(prg_nm)), & 
            ": ERROR soi_txt_get() unable to open ",fl_in(1:ftn_strlen(fl_in))
       do lat_in_idx=1,lat_in_nbr ! NB: Outer loop over lat
          do lon_in_idx=1,lon_in_nbr
             read (fl_in_unit,"(3f10.4)") &
                  mss_frc_snd_in(lon_in_idx,lat_in_idx), &
                  mss_frc_slt_in(lon_in_idx,lat_in_idx), &
                  mss_frc_cly_in(lon_in_idx,lat_in_idx)
             ! Convert input from percent to fraction
             mss_frc_cly_in(lon_in_idx,lat_in_idx)=0.01*mss_frc_cly_in(lon_in_idx,lat_in_idx) ! [%] -> [frc]
             mss_frc_slt_in(lon_in_idx,lat_in_idx)=0.01*mss_frc_slt_in(lon_in_idx,lat_in_idx) ! [%] -> [frc]
             mss_frc_snd_in(lon_in_idx,lat_in_idx)=0.01*mss_frc_snd_in(lon_in_idx,lat_in_idx) ! [%] -> [frc]
          end do ! end loop over lon
       end do ! end loop over lat
       close (fl_in_unit)
    endif                     ! endif Webb soil texture
    
    if (flg_IBIS) then
       ! This IBIS soil dataset is a rebinned subset of the soil dataset assembled 
       ! by the Global Soil Data Task of the Data and Information System (DIS) of the 
       ! International Geosphere-Biosphere Programme (IGBP).
       ! A. J. Carter and R. J. Scholes of Environmentek CSIR South Africa created the
       ! original IGBP dataset and the SoilData program.
       ! Jon Foley and Veronica Fisher of U. Wisconsin generated the multilayer IBIS dataset from the IGBP dataset
       ! I used NCO to extract only the necessary fields into this single layer dataset
       ! The dataset resolution is 1 x 1 degree and the data are written from the North pole to the South pole, East to West, beginning at the date line
       ! One caveat in using the IBIS data is that it may have a different land/sea mask than LSM (Webb) data
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
       if (lat_nbr /= lat_in_nbr) stop "lat_nbr /= lat_in_nbr in soi_txt_get()"
       if (lon_nbr /= lon_in_nbr) stop "lon_nbr /= lon_in_nbr in soi_txt_get()"
       ! Get variable IDs
       rcd=nf90_wrp_inq_varid(nc_id,"tex",txt_IBIS_id)
       rcd=nf90_wrp_inq_varid(nc_id,"lat",lat_in_id)
       ! Get number of dimensions
       rcd=rcd+nf90_inquire_variable(nc_id,txt_IBIS_id,ndims=txt_IBIS_dmn_nbr)
       if (txt_IBIS_dmn_nbr /= 2) then
          write (6,"(2a,i1,a)") prg_nm(1:ftn_strlen(prg_nm)), & 
               ": ERROR soi_txt_get() reports tex has ",txt_IBIS_dmn_nbr," dimensions"
          stop
       endif                  ! endif err
       rcd=rcd+nf90_inquire_attribute(nc_id,txt_IBIS_id,"missing_value",attnum=mss_val_id)
       if (rcd==nf90_noerr) then
          mss_flg=.true.
       else
          mss_flg=.false.
          write (6,"(3a)") prg_nm(1:ftn_strlen(prg_nm)), & 
               ": WARNING soi_txt_get() reports no missing_value in ",fl_in(1:ftn_strlen(fl_in))
          rcd=nf90_noerr
       endif                  ! endif
       rcd=rcd+nf90_get_att(nc_id,txt_IBIS_id,"missing_value",mss_val)
       if (mss_val /= -9999.0) write (6,"(2a,f12.6)") prg_nm(1:ftn_strlen(prg_nm)), &
            ": WARNING soi_txt_get() reports mss_val = ",mss_val
       ! Get data
       rcd=nf90_wrp(nf90_get_var(nc_id,txt_IBIS_id,txt_IBIS_in),"get_var txt_IBIS_in")
       rcd=nf90_wrp(nf90_get_var(nc_id,lat_in_id,lat_in),"get_var lat_in")
       rcd=nf90_wrp_close(nc_id,fl_in,'',sbr_nm=sbr_nm)
       
       lat_mnt_dcr=mnt_dcr_chk(lat_in,lat_in_nbr)
       if (.not.lat_mnt_dcr) stop "ERROR: IBIS lat_in increases in soi_txt_get()"
       do lat_in_idx=1,lat_in_nbr
          do lon_in_idx=1,lon_in_nbr
             if (lat_mnt_dcr) then 
                ! Remap IBIS grid so latitude monotonically increases from South Pole
                lat_out_idx=lat_in_nbr-lat_in_idx+1
             else
                lat_out_idx=lat_in_idx
             endif            ! endif
             ! IBIS files store all 3 soil types in one float using encoding defined in global attribute "history":
             ! "V Fisher, 04/16/98: 2550 = 25%sand, 50%clay, 25%silt"
             if (mss_flg.and.txt_IBIS_in(lon_in_idx,lat_in_idx)==mss_val) then 
                mss_frc_cly_in(lon_in_idx,lat_out_idx)=0.0 ! [frc]
                mss_frc_slt_in(lon_in_idx,lat_out_idx)=0.0 ! [frc]
                mss_frc_snd_in(lon_in_idx,lat_out_idx)=0.0 ! [frc]
             else
                txt_pct_ntg=nint(txt_IBIS_in(lon_in_idx,lat_in_idx))
                cly_pct_ntg=mod(txt_pct_ntg,100) ! [%]
                snd_pct_ntg=nint((txt_pct_ntg-cly_pct_ntg)/100.0) ! [%]
                slt_pct_ntg=100-(cly_pct_ntg+snd_pct_ntg) ! [%]
                mss_frc_cly_in(lon_in_idx,lat_out_idx)=0.01*cly_pct_ntg ! [%] -> [frc]
                mss_frc_slt_in(lon_in_idx,lat_out_idx)=0.01*slt_pct_ntg ! [%] -> [frc]
                mss_frc_snd_in(lon_in_idx,lat_out_idx)=0.01*snd_pct_ntg ! [%] -> [frc]
                if (dbg_lvl==dbg_old) then
                   write (6,"(a,i5,a,3(i2,a))") "code = ",txt_pct_ntg,"%, (clay, sand, silt) = (", &
                        cly_pct_ntg,"%, ",snd_pct_ntg,"%, ",slt_pct_ntg,"% )"
                endif         ! endif dbg
             endif            ! endif
          end do ! end loop over lon
       end do ! end loop over lat
    endif                     ! endif IBIS soil texture
    write (6,"(a,1x,a)") "Read soil texture data from",fl_in(1:ftn_strlen(fl_in))
    
    if (dbg_lvl==dbg_crr) then
       ! Write soil texture data in CSM/CCM/LSM format (three values per line)
       ! These files should be very similar to LSM input files:
       ! diff -c -w /tmp/zender/map/soi_txt_IGBP_dgn.1x1 /fs/cgd/csm/input/lnd/soi_txt_webb.1x1 | m
       ! diff -c -w /tmp/zender/map/soi_txt_webb_dgn.1x1 /fs/cgd/csm/input/lnd/soi_txt_webb.1x1 | m
       open (unit_dgn,file=fl_dgn_soi_txt,status="unknown",iostat=rcd)
       if (rcd /= 0) write (6,"(4a,i4)") prg_nm(1:ftn_strlen(prg_nm)), & 
            ": ERROR soi_txt_get() unable to open ",fl_dgn_soi_txt(1:ftn_strlen(fl_dgn_soi_txt))
       do lat_in_idx=1,lat_in_nbr ! NB: Outer loop over lat
          do lon_in_idx=1,lon_in_nbr
             write (unit_dgn,"(3(f10.4))") &
                  mss_frc_snd_in(lon_in_idx,lat_in_idx)*100.0, & ! [frc] -> [pct]
                  mss_frc_slt_in(lon_in_idx,lat_in_idx)*100.0, & ! [frc] -> [pct]
                  mss_frc_cly_in(lon_in_idx,lat_in_idx)*100.0 ! [frc] -> [pct]
          end do ! end loop over lon
       end do ! end loop over lat
       close (unit_dgn)
       write (6,"(a,1x,a)") "Wrote soil texture in LSM format to",fl_dgn_soi_txt(1:ftn_strlen(fl_dgn_soi_txt))
    endif                     ! endif dbg
    
    ! Set land mask
    do lat_in_idx=1,lat_in_nbr
       do lon_in_idx=1,lon_in_nbr
          if ( mss_frc_snd_in(lon_in_idx,lat_in_idx)==0.0.and. &
               mss_frc_slt_in(lon_in_idx,lat_in_idx)==0.0.and. &
               mss_frc_cly_in(lon_in_idx,lat_in_idx)==0.0) then
             lnd_msk_in(lon_in_idx,lat_in_idx)=0
          else
             lnd_msk_in(lon_in_idx,lat_in_idx)=1
          endif               ! endif
       end do ! end loop over lon
    end do ! end loop over lat
    
    ! Sanity check
    err_nbr=0
    do lat_in_idx=1,lat_in_nbr
       do lon_in_idx=1,lon_in_nbr
          frc_ttl=mss_frc_snd_in(lon_in_idx,lat_in_idx)+ &
               mss_frc_slt_in(lon_in_idx,lat_in_idx)+ &
               mss_frc_cly_in(lon_in_idx,lat_in_idx)
          if ( &
               (mss_frc_snd_in(lon_in_idx,lat_in_idx) < 0.0.or.mss_frc_snd_in(lon_in_idx,lat_in_idx) > 1.0).or. &
               (mss_frc_slt_in(lon_in_idx,lat_in_idx) < 0.0.or.mss_frc_slt_in(lon_in_idx,lat_in_idx) > 1.0).or. &
               (mss_frc_cly_in(lon_in_idx,lat_in_idx) < 0.0.or.mss_frc_cly_in(lon_in_idx,lat_in_idx) > 1.0).or. &
               (lnd_msk_in(lon_in_idx,lat_in_idx)==1.and.abs(frc_ttl-1.0) > eps_rlt)) then
             err_nbr=err_nbr+1
             lon_in_idx_err=lon_in_idx
             lat_in_idx_err=lat_in_idx
          end if              ! endif err
       end do ! end loop over lon
    end do ! end loop over lat
    if (err_nbr > 0) then
       frc_ttl=mss_frc_snd_in(lon_in_idx_err,lat_in_idx_err)+ &
            mss_frc_slt_in(lon_in_idx_err,lat_in_idx_err)+ &
            mss_frc_cly_in(lon_in_idx_err,lat_in_idx_err)
       write (6,"(a,i3,a)") "ERROR soi_txt_get() reports ",err_nbr," soil texture errors during input"
       write (6,"(a)") "Soil texture properties:"
       write (6,"(4(a,7x))") "Sand","Silt","Clay","Total"
       write (6,"(4(f9.6,1x))") &
            mss_frc_snd_in(lon_in_idx_err,lat_in_idx_err), &
            mss_frc_slt_in(lon_in_idx_err,lat_in_idx_err), &
            mss_frc_cly_in(lon_in_idx_err,lat_in_idx_err), &
            frc_ttl
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
    
    ! Process each cell on output grid
    do lat_out_idx=1,lat_out_nbr
       do lon_out_idx=1,lon_out_nbr
          ovr_nbr_crr=ovr_nbr(lon_out_idx,lat_out_idx)
          rnk_1st_idx=mss_val_ntg
          ! Find input cell with largest overlap
          do ovr_idx=1,ovr_nbr_crr
             ovr_wgt_crr=ovr_wgt(lon_out_idx,lat_out_idx,ovr_idx)
             if (ovr_wgt_crr > 0.0) then
                if (rnk_1st_idx==mss_val_ntg) then
                   rnk_1st_idx=ovr_idx ! First time array is > 0.0
                else
                   if (ovr_wgt_crr > ovr_wgt(lon_out_idx,lat_out_idx,rnk_1st_idx)) rnk_1st_idx=ovr_idx
                endif         ! endif
             endif            ! endif
          end do ! end loop over dat
          if (rnk_1st_idx==mss_val_ntg) then
             write (6,"(a,i3,a,f9.5)") "soi_txt_get(): ERROR occured at lon_out_grd(",lon_out_idx,") = ",lon_out_grd(lon_out_idx)
             write (6,"(a,i3,a,f9.5)") "soi_txt_get(): ERROR occured at lat_out_grd(",lat_out_idx,") = ",lat_out_grd(lat_out_idx)
             write (6,"(a,i3)") "soi_txt_get(): # overlapping input gridcells = ",ovr_nbr_crr
             write (6,"(a)") "lon_in_idx lat_in_idx mss_frc_snd mss_frc_slt mss_frc_cly ovr_wgt" 
             do ovr_idx=1,ovr_nbr_crr ! Overlap cell index
                lon_in_idx=ovr_lon_idx(lon_out_idx,lat_out_idx,ovr_idx)
                lat_in_idx=ovr_lat_idx(lon_out_idx,lat_out_idx,ovr_idx)
                write (6,"(i3,1x,i3,1x,i2,1x,f9.6)") lon_in_idx,lat_in_idx, &
                     mss_frc_snd_out(lon_out_idx,lat_out_idx), &
                     mss_frc_slt_out(lon_out_idx,lat_out_idx), &
                     mss_frc_cly_out(lon_out_idx,lat_out_idx), &
                     ovr_wgt(lon_out_idx,lat_out_idx,ovr_idx)
             end do ! end loop over overlapping cells
             stop
          endif               ! endif err
          
          ! GBB recommends against performing a weighted average of soil textures
          ! Assign output cell sand, silt, clay fractions of input cell with largest overlap
          lon_in_idx=ovr_lon_idx(lon_out_idx,lat_out_idx,rnk_1st_idx)
          lat_in_idx=ovr_lat_idx(lon_out_idx,lat_out_idx,rnk_1st_idx)
          ! Warn if datasets do not agree on land mask of cell with largest overlap
          if (dbg_lvl==dbg_old) then
             if (lnd_msk_out(lon_out_idx,lat_out_idx)==1.neqv.lnd_msk_in(lon_in_idx,lat_in_idx)==1) then
                write (6,"(2a)") prg_nm(1:ftn_strlen(prg_nm)), & 
                     ": WARNING soi_txt_get() reports land masks are not equivalent"
                write (6,"(a,i1)") "lnd_msk_out = ",lnd_msk_out(lon_out_idx,lat_out_idx)
                write (6,"(a,i1)") "lnd_msk_in  = ",lnd_msk_in(lon_in_idx,lat_in_idx)
                write (6,"(a)") "Output cell edge locations:"
                write (6,"(a)") "(idx)   Lat Sth (idx)   Lat Nrt (idx)   Lon Wst (idx)   Lon Est"
                write (6,"(4(a1,i3,a1,1x,f9.4,1x))") &
                     "(",lat_out_idx,")",lat_out_grd(lat_out_idx),"(",lat_out_idx+1,")",lat_out_grd(lat_out_idx+1), &
                     "(",lon_out_idx,")",lon_out_grd(lon_out_idx),"(",lon_out_idx+1,")",lon_out_grd(lon_out_idx+1)
             endif ! endif landmask error
          endif ! endif dbg
          ! Use texture data over land only
          if (lnd_msk_out(lon_out_idx,lat_out_idx)==1) then
             !        if (lnd_msk_in(lon_in_idx,lat_in_idx)==1) then
             mss_frc_snd_out(lon_out_idx,lat_out_idx)=mss_frc_snd_in(lon_in_idx,lat_in_idx) ! [frc]
             mss_frc_slt_out(lon_out_idx,lat_out_idx)=mss_frc_slt_in(lon_in_idx,lat_in_idx) ! [frc]
             mss_frc_cly_out(lon_out_idx,lat_out_idx)=mss_frc_cly_in(lon_in_idx,lat_in_idx) ! [frc]
          endif               ! endif 
          ! Process texture data in accordance with LSM algorithms
          ! NB: Following code from LSM:soiltex.F
          if (sfc_typ_LSM_out(lon_out_idx,lat_out_idx)==0.or.sfc_typ_LSM_out(lon_out_idx,lat_out_idx)==1) then
             ! Ocean and glaciers have no soil
             mss_frc_snd_out(lon_out_idx,lat_out_idx)=0.0 ! [frc]
             mss_frc_slt_out(lon_out_idx,lat_out_idx)=0.0 ! [frc]
             mss_frc_cly_out(lon_out_idx,lat_out_idx)=0.0 ! [frc]


             ! If LSM = land, but input data has no soil: 
             !   If flg_soi_nn==FALSE (original method), set soil to LOAM
             !   If flg_soi_nn==TRUE (MF added) then set undefined land soil type 
             !   to land type of nearest (originally-defined) land gridcell
             ! NB: This block is key for making sure IBIS has soil data on LSM grid
          else if ( &          ! not Ocean or Ice,
               (mss_frc_snd_out(lon_out_idx,lat_out_idx)==0.0).and. &
               (mss_frc_slt_out(lon_out_idx,lat_out_idx)==0.0).and. &
               (mss_frc_cly_out(lon_out_idx,lat_out_idx)==0.0)) then
             
             loam_nbr=loam_nbr+1
             if(.not.flg_soi_nn) then
                
                mss_frc_snd_out(lon_out_idx,lat_out_idx)=mss_frc_snd_loam ! [frc] Loam
                mss_frc_slt_out(lon_out_idx,lat_out_idx)=mss_frc_slt_loam ! [frc] Loam 
                mss_frc_cly_out(lon_out_idx,lat_out_idx)=mss_frc_cly_loam ! [frc] Loam
                
                if (dbg_lvl==dbg_old) then
                   write (6,"(a,i2)") "soi_txt_get(): Loam texture to surface type ",sfc_typ_LSM_out(lon_out_idx,lat_out_idx)
                   write (6,"(a)") "Output cell edge locations:"
                   write (6,"(a)") "(idx)   Lat Sth (idx)   Lat Nrt (idx)   Lon Wst (idx)   Lon Est"
                   write (6,"(4(a1,i3,a1,1x,f9.4,1x))") &
                        "(",lat_out_idx,")",lat_out_grd(lat_out_idx),"(",lat_out_idx+1,")",lat_out_grd(lat_out_idx+1), &
                        "(",lon_out_idx,")",lon_out_grd(lon_out_idx),"(",lon_out_idx+1,")",lon_out_grd(lon_out_idx+1)
                endif            ! endif dbg
                
                
             else
                ! Nearest neighbor
                distance_min = 100000    ! [km] Initialized to high value
                lat1_tmp = lat_in_grd(lat_in_idx) / degtorad  ! [nbr] radians
                lon1_tmp = lon_in_grd(lon_in_idx) / degtorad  ! [nbr] radians
                
                ! Find nearest neighbor
                do lat_idx2=1,lat_in_nbr
                   do lon_idx2=1,lon_in_nbr
                      ! If neighboring cell on input grid has defined soil texture
                      if ( (mss_frc_snd_in(lon_idx2,lat_idx2).gt.0.0).or. &
                           (mss_frc_slt_in(lon_idx2,lat_idx2).gt.0.0).or. &
                           (mss_frc_cly_in(lon_idx2,lat_idx2).gt.0.0) )   then
                         
                         lat2_tmp = lat_in_grd(lat_idx2) / degtorad ! [nbr] radians
                         lon2_tmp = lon_in_grd(lon_idx2) / degtorad ! [nbr] radians
            
                         ! Great Circle Distance Formula (kilometers)
                         distance_tmp=6378.7*acos((sin(lat1_tmp)*sin(lat2_tmp))+ &
                              (cos(lat1_tmp)*cos(lat2_tmp)*cos(abs(lon2_tmp-lon1_tmp)))) ! [km]
                         if (distance_tmp.lt.distance_min) then
                            distance_min = distance_tmp
                            lat_idx_min_dst = lat_idx2 
                            lon_idx_min_dst = lon_idx2
                         endif        ! endif
                      endif           ! endif defined soil texture
                   enddo ! end loop over lon
                enddo ! end loop over lat
        
                mss_frc_snd_out(lon_out_idx,lat_out_idx)=mss_frc_snd_in(lon_idx_min_dst,lat_idx_min_dst)
                mss_frc_slt_out(lon_out_idx,lat_out_idx)=mss_frc_slt_in(lon_idx_min_dst,lat_idx_min_dst)
                mss_frc_cly_out(lon_out_idx,lat_out_idx)=mss_frc_cly_in(lon_idx_min_dst,lat_idx_min_dst)
                
             endif            ! endif not nearest neighbor algorithm
          endif               ! endifnot Ocean or Ice     
          
          if (dbg_lvl==dbg_old.and.lat_out_idx==lat_idx_dbg.and.lon_out_idx==lon_idx_dbg) then
             write (6,"(a)") "DEBUG: soi_txt_get()"
             write (6,"(a,f5.2)") "mss_frc_cly_out = ",mss_frc_cly_out(lon_out_idx,lat_out_idx)
             write (6,"(a)") "(idx)   Lat Sth (idx)   Lat Nrt (idx)   Lon Wst (idx)   Lon Est"
             write (6,"(4(a1,i3,a1,1x,f9.4,1x))") &
                  "(",lat_out_idx,")",lat_out_grd(lat_out_idx),"(",lat_out_idx+1,")",lat_out_grd(lat_out_idx+1), &
                  "(",lon_out_idx,")",lon_out_grd(lon_out_idx),"(",lon_out_idx+1,")",lon_out_grd(lon_out_idx+1)
          endif               ! endif dbg

          ! Sanity checks
          frc_ttl=mss_frc_snd_out(lon_out_idx,lat_out_idx)+ & ! [frc]
               mss_frc_slt_out(lon_out_idx,lat_out_idx)+ &
               mss_frc_cly_out(lon_out_idx,lat_out_idx)
                       
          if (sfc_typ_LSM_out(lon_out_idx,lat_out_idx) > 1) then
             if (abs(frc_ttl-1.0) >= eps_rlt) then
                write (6,"(a,f15.8)") "soi_txt_get(): ERROR Output textures sum to ",frc_ttl
                write (6,"(a)") "Output cell edge locations:"
                write (6,"(a)") "(idx)   Lat Sth (idx)   Lat Nrt (idx)   Lon Wst (idx)   Lon Est"
                write (6,"(4(a1,i3,a1,1x,f9.4,1x))") &
                     "(",lat_out_idx,")",lat_out_grd(lat_out_idx),"(",lat_out_idx+1,")",lat_out_grd(lat_out_idx+1), &
                     "(",lon_out_idx,")",lon_out_grd(lon_out_idx),"(",lon_out_idx+1,")",lon_out_grd(lon_out_idx+1)
                stop
             endif            ! endif
          endif               ! endif dry land

       end do ! end loop over lon
    end do ! end loop over lat

    if (loam_nbr > 0) then
       if (flg_soi_nn) then
          write (6,"(a,i4,a)") "DEBUG: soi_txt_get() reports ",loam_nbr," points set to nearest neighbor"
       else
          write (6,"(a,i4,a)") "DEBUG: soi_txt_get() reports ",loam_nbr," points set to loam"
       endif
    endif         ! endif


    ! Write output file in text format?  Requires 2x2 output resolution
    if (.false.) then
       ! Write soil texture data in CSM/CCM/LSM format (three values per line)
       ! These files should be very similar to LSM input files:
       ! Assumes input data is on 2x2 grid. Writes to 1x1 grid
       ! diff -c -w /tmp/zender/map/soi_txt_IGBP_dgn.1x1 /fs/cgd/csm/input/lnd/soi_txt_webb.1x1 | m
       ! diff -c -w /tmp/zender/map/soi_txt_webb_dgn.1x1 /fs/cgd/csm/input/lnd/soi_txt_webb.1x1 | m
       open (unit_dgn,file=fl_soi_txt_out,status="unknown",iostat=rcd)
       if (rcd /= 0) write (6,"(4a,i4)") prg_nm(1:ftn_strlen(prg_nm)), & 
            ": ERROR soi_txt_get() unable to open ",fl_soi_txt_out(1:ftn_strlen(fl_dgn_soi_txt))
       
       do lat_out_idx=1,lat_out_nbr ! NB: Outer loop over lat
          do temp_idx=1,2
             do lon_out_idx=1,lon_out_nbr
                write (unit_dgn,"(3(f10.4))") &
                     mss_frc_snd_out(lon_out_idx,lat_out_idx)*100.0, & ! [frc] -> [pct]
                     mss_frc_slt_out(lon_out_idx,lat_out_idx)*100.0, & ! [frc] -> [pct]
                     mss_frc_cly_out(lon_out_idx,lat_out_idx)*100.0 ! [frc] -> [pct]
                write (unit_dgn,"(3(f10.4))") &
                     mss_frc_snd_out(lon_out_idx,lat_out_idx)*100.0, & ! [frc] -> [pct]
                     mss_frc_slt_out(lon_out_idx,lat_out_idx)*100.0, & ! [frc] -> [pct]
                     mss_frc_cly_out(lon_out_idx,lat_out_idx)*100.0 ! [frc] -> [pct]
             end do ! end loop over lon
          end do ! end loop tmp_idx
       end do ! end loop over lat
       close (unit_dgn)
       write (6,"(a,1x,a)") "Wrote soil texture in LSM format to",fl_soi_txt_out(1:ftn_strlen(fl_dgn_soi_txt))
    endif                     ! endif write output file



    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Exiting soi_txt_get()"
    return 
  end subroutine soi_txt_get
  
  subroutine soi_hyd_prc( & ! [sbr] Hydrologic properties of soil texture
       lat,lat_nbr,lon,lon_nbr,time_nbr, & ! I
       mss_frc_cly,mss_frc_slt,mss_frc_snd, & ! I
       vwc_sfc, & ! I
       smp_sat,smp_sfc,vwc_dry,vwc_opt,vwc_rel,vwc_sat) ! O
    ! Purpose: Compute and return hydrologic properties
    ! soi_hyd_prc() is called by bds_prc()
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use map_grd ! [mdl] Map grids and regridding
    implicit none
    ! Parameters
    ! Commons
    ! Input
    integer lat_nbr ! [nbr] Number of latitudes
    integer lon_nbr ! [nbr] Number of longitudes
    integer time_nbr ! [nbr] Number of times
    real lat(lat_nbr) ! [dgr] Centered latitudes
    real lon(lon_nbr) ! [dgr] Centered longitudes
    real mss_frc_cly(lon_nbr,lat_nbr) ! [frc] Soil texture clay
    real mss_frc_slt(lon_nbr,lat_nbr) ! [frc] Soil texture silt
    real mss_frc_snd(lon_nbr,lat_nbr) ! [frc] Soil texture sand
    real vwc_sfc(lon_nbr,lat_nbr,time_nbr) ! [m3 m-3] Volumetric water content
    ! Output
    real smp_sat(lon_nbr,lat_nbr) ! [mm H2O] Saturated soil matric potential (sand-dependent)
    real smp_sfc(lon_nbr,lat_nbr,time_nbr) ! [mm H2O] Soil matric potential
    real vwc_dry(lon_nbr,lat_nbr) ! [m3 m-3] Dry volumetric water content (no E-T)
    real vwc_opt(lon_nbr,lat_nbr) ! [m3 m-3] E-T optimal volumetric water content 
    real vwc_rel(lon_nbr,lat_nbr,time_nbr) ! [frc] Water content relative to saturation
    real vwc_sat(lon_nbr,lat_nbr) ! [m3 m-3] Saturated volumetric water content (sand-dependent)
    ! Local
    integer lat_idx           ! [idx] Counting index
    integer lon_idx           ! [idx] Counting index
    integer time_idx          ! [idx] Counting index
    real smp_xpn_b(lon_nbr,lat_nbr) ! [frc] Exponent "b" for smp (clay-dependent)

    ! Main Code
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Entering soi_hyd_prc()"
    ! Initialize arrays
    smp_sat=0.0               ! [m3 m-3]
    vwc_dry=0.0               ! [m3 m-3]
    vwc_opt=0.0               ! [m3 m-3]
    vwc_sat=0.0               ! [m3 m-3]
    smp_sfc=0.0               ! [mm H2O]
    vwc_rel=0.0               ! [frc]
    ! Sanity checks
    if (time_nbr /= 12) stop "ERROR: time_nbr /= 12 in soi_hyd_prc()"
    
    ! Initialize time-invariant soil hydraulic properties
    ! See GBB96 p. 98, implemented in CCM:lsm/lsmtci()
    do lat_idx=1,lat_nbr
       do lon_idx=1,lon_nbr
          smp_xpn_b(lon_idx,lat_idx)= & ! [frc] Exponent "b" for smp (clay-dependent)
               2.91+0.159*mss_frc_cly(lon_idx,lat_idx)*100.0
          ! NB: Adopt convention that matric potential is positive definite
          smp_sat(lon_idx,lat_idx)= & ! [mm H2O] Saturated soil matric potential (sand-dependent)
               10.0*(10.0**(1.88-0.0131*mss_frc_snd(lon_idx,lat_idx)*100.0))
          vwc_sat(lon_idx,lat_idx)= & ! [m3 m-3] Saturated volumetric water content (sand-dependent)
               0.489-0.00126*mss_frc_snd(lon_idx,lat_idx)*100.0
          vwc_dry(lon_idx,lat_idx)= & ! [m3 m-3] Dry volumetric water content (no E-T)
               vwc_sat(lon_idx,lat_idx)* &
               (316230.0/smp_sat(lon_idx,lat_idx))**(-1.0/smp_xpn_b(lon_idx,lat_idx))
          vwc_opt(lon_idx,lat_idx)= & ! [m3 m-3] E-T optimal volumetric water content
               vwc_sat(lon_idx,lat_idx)* &
               (158490.0/smp_sat(lon_idx,lat_idx))**(-1.0/smp_xpn_b(lon_idx,lat_idx))
       end do ! end loop over lon
    end do ! end loop over lat
    
    ! Construct time-varying hydraulic properties
    ! See GBB96 p. 97
    do time_idx=1,time_nbr
       do lat_idx=1,lat_nbr
          do lon_idx=1,lon_nbr
             vwc_rel(lon_idx,lat_idx,time_idx)= & ! [frc] Water content relative to saturation
                  max( &
                  vwc_sfc(lon_idx,lat_idx,time_idx)/vwc_sat(lon_idx,lat_idx), &
                  0.05)
             ! NB: Adopt convention that matric potential is positive definite
             smp_sfc(lon_idx,lat_idx,time_idx)= & ! [mm H2O] Soil matric potential
                  smp_sat(lon_idx,lat_idx)* &
                  (vwc_rel(lon_idx,lat_idx,time_idx)**(-smp_xpn_b(lon_idx,lat_idx)))
          end do ! end loop over lon
       end do ! end loop over lat
    end do ! end loop over time
    
    if (dbg_lvl >= dbg_sbr) write (6,"(a)") "Exiting soi_hyd_prc()"
    return 
  end subroutine soi_hyd_prc
  
end module soi_txt ! [mdl] Soil texture and hydrology
