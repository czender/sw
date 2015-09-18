module cape_mdl
  !use shr_kind_mod,only:r8=>shr_kind_r8 ! [mdl] Precision r8, i8, ...
  use tdy_mdl,only:svp_H2O_lqd_prk78,svp_H2O_ice_prk78 ! [mdl] Thermodynamics
  ! [mdl] Fundamental and derived physical constants
  use phys_cst_mdl,only:gas_cst_dry_air,tpt_frz_pnt,spc_heat_dry_air, &
       gas_cst_H2O,cp_vpr_rcp_cp_dry_m1,eps_H2O,eps_H2O_rcp_m1, &
       spc_heat_H2O_vpr,kappa_dry_air
  
  implicit none

contains

subroutine pcl_tdy_pfl(alt, &
     grv, &
     spc_heat_mst_air, &
     tpt_nvr, &
     tpt_sfc, &
     CAPE, &
     CINE, &     
     lev_LCL, &
     lev_LFC, &
     lev_LNB, &
     lev_nbr, &
     nrg_dry, &
     nrg_mst, &
     pcl_DLR, &
     pcl_MLR, &
     pcl_PLR, &
     prs, &
     prs_sfc, &
     r_H2O, &
     RH, &
     tpt_dwp, &
     tpt_pcl, &
     tpt_pcl_sfc)
  
  implicit none

  ! Input
  real,dimension(:),intent(in)::alt              ! [m] Midlayer Height of Measurement
  real,dimension(:),intent(in)::grv              ! [m s-2] Midlayer Acceleration Due to Gravity
  integer          ,intent(in)::lev_nbr          ! [nbr] Number of levels in the sounding file 
  real,dimension(:),intent(in)::prs              ! [Pa] Midlayer Pressure
  real,dimension(:),intent(in)::r_H2O            ! [kg kg-1] Dry mass water vapor mixing ratio
  real,dimension(:),intent(in)::RH               ! [frc] Midlayer Relative Humidity
  real,dimension(:),intent(in)::spc_heat_mst_air ! [J kg-1 K-1] Midlayer Specific Heat of Moist Air
  real,dimension(:),intent(in)::tpt_nvr          ! [K] Midlayer Environment Air Temperature
  real             ,intent(in)::tpt_sfc          ! [K] Surface Temperature
  real             ,intent(in)::prs_sfc          ! [Pa] Surface Pressure
  
  ! Output
  real,intent(out)::CAPE             ! [J kg-1] Convective Available Potential Energy
  real,intent(out)::CINE             ! [J kg-1] Convective Inhibition Energy
  real,intent(out)::lev_LFC          ! [m] Level of Free Convection
  real,intent(out)::lev_LCL          ! [m] Lifted Condensation Level
  real,intent(out)::lev_LNB          ! [m] Level of Neutral Buoyancy
  real,intent(out)::nrg_dry(lev_nbr) ! [J kg-1] Dry static energy
  real,intent(out)::nrg_mst(lev_nbr) ! [J kg-1] Moist static energy
  real,intent(out)::pcl_DLR(lev_nbr) ! [K m-1] Dry Adiabatic Lapse Rate
  real,intent(out)::pcl_MLR(lev_nbr) ! [K m-1] Moist (Unsaturated) Adiabatic Lapse Rate
  real,intent(out)::pcl_PLR(lev_nbr) ! [K m-1] Pseudo Adiabatic Lapse Rate
  real,intent(out)::tpt_dwp(lev_nbr) ! [K] Dewpoint Temperature
  real,intent(out)::tpt_pcl(lev_nbr) ! [K] Air Parcel Temperature
  real,intent(out)::tpt_pcl_sfc      ! [K] Air Parcel Surface Temperature

  ! Static Variables
  !real lev_100      ! [m] Height at top of lowest 100mb layer CEWI
  !real H2O_100      ! [Pa] Vapor Pressure at top of lowest 100mb layer
  !real H2O_b        ! [m] Intercept of vapor pressure
  !real H2O_slope    ! [m Pa-1] Slope of vapor pressure
  real ppr_H2O_sfc  ! [Pa] Surface Air Vapor Pressure
  !real prs_100      ! [Pa] Pressure at top of lowest 100mb=10000Pa layer
  real r_H2O_sfc    ! [kg kg-1] Dry mass water vapor mixing ratio of lowest 10000Pa layer
  real svp_H2O_pcl  ! [Pa] Air Parcel Saturation Vapor Pressure
  !real theta_nxt_lyr ! [K] Potential temperature of next layer
  !real theta_prv_lyr ! [K] Potential temperature of previous layer
  real tpt_LFC      ! [K] Temperature of parcel and environment at LFC
  real tpt_LNB      ! [K] Temperature at Level of Neutral Buoyancy
  !real tpt_100      ! [K] Temperature at top of lowest 100mb layer
  !real tpt_b        ! [m] Intercept of parcel temperature
  !real tpt_slope    ! [m K-1] Slope of environment temperature
  !real wgt_hgt

  !Parameters
  real,parameter::hgt_sfc=0. ! [m] Surface Height
    
  integer idx ! [idx] Counting Index For Midlayer Levels
  logical LCL ! Lifted Condensation Level Reached
  logical LFC ! Level of Free Convection Reached
  logical LNB ! Level of Neutral Buoyancy Reached
  
  ! Initialize variables
  pcl_DLR(:)=0.! SCAPPS 2/21/07
  pcl_MLR(:)=0.! SCAPPS 2/21/07
  pcl_PLR(:)=0.! SCAPPS 2/21/07
  tpt_pcl(:)=0.! SCAPPS 2/21/07
  tpt_dwp(:)=0.! SCAPPS 2/21/07
  nrg_dry(:)=0.! SCAPPS 2/21/07
  nrg_mst(:)=0.! SCAPPS 2/21/07
  
  !lev_100=0.      ! [m] Height at top of lowest 100mb layer CEWI
  LFC    =.false. ! [flg] Level of Free Convection Reached 
  LCL    =.false. ! [flg] Lifted Condensation Level Reached
  LNB    =.false. ! [flg] Level of Neutral Buoyancy Reached
  CINE   =0.      ! SCAPPS 2/21/07
  CAPE   =0.      ! SCAPPS 2/21/07
  
  do idx=lev_nbr,1,-1
     pcl_DLR(idx)=(grv(idx)/spc_heat_dry_air)      ! [K m-1] Midlayer Dry adiabatic lapse rate   
     pcl_MLR(idx)=(grv(idx)/spc_heat_mst_air(idx)) ! [K m-1] Midlayer Moist (Unsaturated) adiabatic lapse rate
     tpt_dwp(idx)=tpt_nvr(idx)/(1.0+(-gas_cst_H2O*tpt_nvr(idx)* &
          log(RH(idx))/spc_ltn_heat(tpt_nvr(idx)))) ! [K] Tsonis P. 100
     if (tpt_dwp(idx)>tpt_nvr(idx)) then
        tpt_dwp(idx)=tpt_nvr(idx)
     end if ! end if tpt_dwp(idx)>tpt_nvr(idx)
  enddo ! end loop over lev_nbr
  

  !Lift a parcel with weighted mean values of temp and moisture in lowest 100mb layer
  !The parcels above 1000mb are warmed adiabatically as they descend to 1000mb  
  !prs_100=prs_sfc-10000. ! [Pa] Pressure at top of 10000Pa layer
  !do idx=lev_nbr,1,-1
  !   if (idx==lev_nbr)then
  !      if (prs(idx)==prs_100) then
  !         lev_100=alt(idx) ! [m] 10000Pa level
  !      elseif (prs(idx)<prs_100 .AND. prs_sfc>prs_100) then
  !         lev_100=hgt_sfc+(scl_hgt(idx)*log(prs_sfc/prs_100))
  !      end if ! end if prs(idx)==prs_100
  !   else
  !      if (prs(idx)==prs_100) then
  !         lev_100=alt(idx) ! [m] 10000Pa level
  !      elseif (prs(idx)<prs_100 .AND. prs(idx+1)>prs_100) then
  !         lev_100=alt(idx+1)+(scl_hgt(idx+1)*log(prs(idx+1)/prs_100))
  !      end if ! end if prs(idx)==prs_100
  !   end if
  !enddo ! end loop over lev_nbr

  !do idx=lev_nbr,1,-1
  !   if (idx==lev_nbr)then
  !      wgt_hgt=alt(idx)/lev_100
        !theta_prv_lyr=tpt_sfc*((prs_sfc/100.)/(prs(idx)/100.))**kappa_dry_air
        !theta_nxt_lyr=tpt_nvr(idx)*((prs_sfc/100.)/(prs(idx)/100.))**kappa_dry_air
        !tpt_pcl_sfc=0.5*(theta_prv_lyr+theta_nxt_lyr)*wgt_hgt
  !      tpt_pcl_sfc=0.5*(tpt_sfc+tpt_ptn(idx))*wgt_hgt
  !      r_H2O_sfc=r_H2O(idx)*wgt_hgt
  !   elseif (alt(idx)<=lev_100) then
  !      wgt_hgt=(alt(idx)-alt(idx+1))/lev_100
        !theta_prv_lyr=tpt_nvr(idx+1)*((prs_sfc/100.)/(prs(idx+1)/100.))**kappa_dry_air
        !theta_nxt_lyr=tpt_nvr(idx)*((prs_sfc/100.)/(prs(idx)/100.))**kappa_dry_air
        !tpt_pcl_sfc=tpt_pcl_sfc+(0.5*(theta_prv_lyr+theta_nxt_lyr)*wgt_hgt)
  !      tpt_pcl_sfc=tpt_pcl_sfc+(0.5*(tpt_ptn(idx+1)+tpt_ptn(idx))*wgt_hgt)
  !      r_H2O_sfc=r_H2O_sfc+(0.5*(r_H2O(idx+1)+r_H2O(idx))*wgt_hgt)
  !   elseif (alt(idx)>lev_100 .AND. alt(idx+1)<lev_100) then
        !Calculate Slope
  !      tpt_slope=(alt(idx)-alt(idx+1))/(tpt_nvr(idx)-tpt_nvr(idx+1))
        !Calculate y-intercept
  !      tpt_b=alt(idx+1)-(tpt_slope*tpt_nvr(idx+1))
        !Calculate temp and moisture at 100mb level
  !      tpt_100=(lev_100-tpt_b)/tpt_slope
        !Calculate Slope
  !      H2O_slope=(alt(idx)-alt(idx+1))/(r_H2O(idx)-r_H2O(idx+1))
        !Calculate y-intercept
  !      H2O_b =alt(idx+1)-(H2O_slope*r_H2O(idx+1))
        !Calculate temp and moisture at 100mb level
  !      H2O_100=(lev_100-H2O_b)/H2O_slope
  !      wgt_hgt=(lev_100-alt(idx+1))/lev_100
        !theta_prv_lyr=tpt_nvr(idx+1)*((prs_sfc/100.)/(prs(idx+1)/100.))**kappa_dry_air
        !theta_nxt_lyr=tpt_100*((prs_sfc/100.)/(prs_100/100.))**kappa_dry_air
        !tpt_pcl_sfc=tpt_pcl_sfc+(0.5*(theta_prv_lyr+theta_nxt_lyr)*wgt_hgt)
  !      tpt_pcl_sfc=tpt_pcl_sfc+(0.5*(tpt_ptn(idx+1)+tpt_ptn(idx))*wgt_hgt)
  !      r_H2O_sfc=r_H2O_sfc+(0.5*(r_H2O(idx+1)+H2O_100)*wgt_hgt)
  !   end if ! End if idx==lev_nbr
  !enddo ! end loop over lev_nbr
  
  tpt_pcl_sfc=tpt_sfc
  r_H2O_sfc=r_H2O(lev_nbr)

  !Calculate water vapor pressure at surface
  ppr_H2O_sfc=prs_sfc*r_H2O_sfc/(eps_H2O+r_H2O_sfc) ! [Pa] Water vapor pressure of lowest 10000Pa layer
  
  ! Determine Air Parcel Properties At Surface
  svp_H2O_pcl=svp_H2O_get(tpt_pcl_sfc) ! [Pa] Saturation Vapor Pressure
  
  if (ppr_H2O_sfc>=svp_H2O_pcl) then
     lev_LCL=0.0 ! [m] Lifted Condensation Level
     LCL=.TRUE.  ! [flg] Lifted Condensation Level Reached
     
     if(tpt_pcl_sfc>tpt_sfc)then
        LFC=.TRUE.  ![flg] Level of free convection
     end if
  end if ! End if ppr_H2O_sfc>=svp_H2O_pcl

  do idx=lev_nbr,1,-1
     if (idx==lev_nbr) then
        ! Lift Parcel To First MidLayer From Surface
        call lift_pcl(CAPE, &
             CINE, &
             LCL, &
             LFC, &
             LNB, &
             grv(idx), &
             grv(idx), &
             alt(idx), &
             hgt_sfc, &
             lev_LCL, &
             lev_LFC, &
             lev_LNB, &
             nrg_dry, &
             nrg_mst, &
             ppr_H2O_sfc, &
             pcl_MLR(idx), &
             pcl_PLR, &
             prs_sfc, &
             prs(idx), &
             r_H2O_sfc, &
             tpt_LFC, &
             tpt_LNB, &
             tpt_nvr(idx), &
             tpt_sfc, &
             tpt_pcl(idx), &
             tpt_pcl_sfc, &
             idx)
     else
        call lift_pcl(CAPE, &
             CINE, &       
             LCL, &
             LFC, &
             LNB, &
             grv(idx), &
             grv(idx+1), &
             alt(idx), &
             alt(idx+1), &
             lev_LCL, &
             lev_LFC, &
             lev_LNB, &
             nrg_dry, &
             nrg_mst, &
             ppr_H2O_sfc, &
             pcl_MLR(idx), &
             pcl_PLR, &
             prs(idx+1), &
             prs(idx), &
             r_H2O_sfc, &
             tpt_LFC, &
             tpt_LNB, &
             tpt_nvr(idx), &
             tpt_nvr(idx+1), &
             tpt_pcl(idx), &
             tpt_pcl(idx+1), &
             idx)

     end if ! end if idx==lev_nbr
  enddo ! end loop over lev_nbr

  CINE=ABS(CINE)
  if (CAPE<=0.0) then
     CAPE=-1.0
     lev_LFC=-1.0
     lev_LNB=-1.0
  end if ! end if CAPE<=0.0     
end subroutine pcl_tdy_pfl

subroutine lift_pcl(CAPE, &
     CINE, &
     LCL,&
     LFC, &
     LNB, &
     grv_nxt_lyr, &
     grv_prv_lyr, &
     hgt_nxt_lyr, &
     hgt_prv_lyr, &
     lev_LCL, &
     lev_LFC, &
     lev_LNB, &
     nrg_dry, &
     nrg_mst, &
     ppr_H2O_sfc, &
     pcl_MLR, &
     pcl_PLR, &
     prs_prv_lyr, &
     prs_nxt_lyr, &
     r_H2O_sfc, &
     tpt_LFC, &
     tpt_LNB, &
     tpt_nvr_nxt_lyr, &
     tpt_nvr_prv_lyr, &
     tpt_pcl_nxt_lyr, &
     tpt_pcl_prv_lyr, &
     idx)
   
  implicit none
  
  ! In Variables
  real,intent(in)::grv_nxt_lyr     ! [m s-2] Next Midlayer Acceleration Due to Gravity
  real,intent(in)::grv_prv_lyr     ! [m s-2] Previous Midlayer Acceleration Due to Gravity
  integer,intent(in)::idx          ! [idx] Counting Index For Midlayer Levels
  real,intent(in)::hgt_nxt_lyr     ! [m] Height of Next Midlayer Measurement
  real,intent(in)::hgt_prv_lyr     ! [m] Height of Previous Midlayer Measurement
  real,intent(in)::pcl_MLR         ! [K m-1] Moist Adiabatic Lapse Rate
  real,intent(in)::ppr_H2O_sfc     ! [Pa] Lowest 10000Pa layer H2O Partial Pressure
  real,intent(in)::prs_prv_lyr     ! [Pa] Midlayer Pressure
  real,intent(in)::prs_nxt_lyr     ! [Pa] Next Midlayer Pressure
  real,intent(in)::r_H2O_sfc       ! [kg kg-1] Dry mass water vapor mixing ratio at lowest 10000Pa
  real,intent(in)::tpt_nvr_nxt_lyr ! [K] Environment Temperature at Next Midlayer
  real,intent(in)::tpt_nvr_prv_lyr ! [K] Environment Temperature at Previous Midlayer
  real,intent(in)::tpt_pcl_prv_lyr ! [K] Previous Layer Parcel Temperature

  ! Inout Variables
  logical,intent(inout)::LCL  ! Lifted Condensation Level Reached
  logical,intent(inout)::LFC  ! Level of Free Convection Reached
  real,intent(inout)::lev_LCL ! [m] Height of Lifted Condensation Level
  real,intent(inout)::lev_LFC ! Level of Free Convection Height
  real,dimension(:),intent(inout)::pcl_PLR ! [K m-1] Pseudo-Adiabatic Lapse Rate
  real,intent(inout)::tpt_LFC ! [K] Temperature at Level of Free Convection

  ! Out Variables
  real,intent(inout)::CAPE            ! [J kg-1] Convective Available Potential Energy
  real,intent(inout)::CINE            ! [J kg-1] Convective Inhibition Energy 
  logical,intent(inout)::LNB          ! Level of Neutral Buoyancy Reached
  real,intent(out)::tpt_LNB           ! [K] Temperature at Level of Neutral Buoyancy
  real,intent(out)::lev_LNB           ! [m] Level of Neutral Buoyancy Height
  real,dimension(:),intent(inout)::nrg_dry ! [J kg-1] Dry static energy
  real,dimension(:),intent(inout)::nrg_mst ! [J kg-1] Moist static energy
  real,intent(out)::tpt_pcl_nxt_lyr   ! [K] Parcel Temperature at Next Midlayer Height
  
  ! Static Variables
  real CAPE_dltz ! [m] Thickness of LFC column
  real CINE_dltz ! [m] Thickness of LFC column
  real CAPE_fnc_val_prv
  real CAPE_fnc_val_nxt
  real CINE_fnc_val_prv
  real CINE_fnc_val_nxt
  real hgt             ! [m] Height between previous and next midlayer altitudes
  real PLR             ! [K m-1] Pseudo-adiabatic lapse rate
  real b_int           ! [m] Intercept
  real prs_LCL         ! [Pa] LCL pressure
  real slope           ! [m Pa-1] Slope dlt_hgt/dlt_prs
  real spc_heat_mst_air! [J kg-1 K-1] Specific heat of moist air
  real svp_LCL         ! [Pa] Saturation Vapor Pressure at LCL
  real svp_H2O_prv_lyr ! [Pa] Saturation Vapor Pressure for H2O previous layer
  real svp_H2O_nxt_lyr ! [Pa] Saturation Vapor Pressure for H2O next layer
  real tpt_LCL         ! [K] LCL temperature
  real tpt_nvr_LCL     ! [K] environment temperature at LCL
  
  hgt=hgt_nxt_lyr-hgt_prv_lyr ! [m] Layer Thickness
  
  svp_H2O_prv_lyr=svp_H2O_get(tpt_pcl_prv_lyr)      ! [Pa] Saturation Vapor Pressure
  
  if (LCL) then
     call lift_pcl_PLR(tpt_pcl_prv_lyr,prs_prv_lyr,grv_prv_lyr,hgt,PLR,tpt_pcl_nxt_lyr)
     pcl_PLR(idx)=PLR
     !Calculate Dry and Moist Static Energy
     nrg_dry(idx)=(spc_heat_dry_air*tpt_pcl_nxt_lyr)+(grv_nxt_lyr*hgt_nxt_lyr) ! [J kg-1] Dry static energy
     
     ! 0.87 [frc] for computing moist specific heatIrG81 pp. 77
     !cp_vpr_rcp_cp_dry_m1=spc_heat_H2O_vpr/spc_heat_dry_air-1.0 
     
     spc_heat_mst_air=spc_heat_dry_air*(1.0+cp_vpr_rcp_cp_dry_m1* &
          (pcl_H2O_svm_get(tpt_pcl_nxt_lyr,prs_nxt_lyr)))
     nrg_mst(idx)=(spc_heat_mst_air*tpt_pcl_nxt_lyr)+(grv_nxt_lyr*hgt_nxt_lyr)+ &
          ((pcl_H2O_svm_get(tpt_pcl_nxt_lyr,prs_nxt_lyr))* &
          spc_ltn_heat(tpt_pcl_nxt_lyr)) ! [J kg-1] Moist static energy
  else
     ! Use pcl_MLR To Lift Parcel To Next Level
     tpt_pcl_nxt_lyr=tpt_pcl_prv_lyr-(pcl_MLR*hgt) ! [K] Lowest Level Unsaturated Parcel Temperature
     
     pcl_PLR(idx)=0.0

     !Calculate Dry and Moist Static Energy
     nrg_dry(idx)=(spc_heat_dry_air*tpt_pcl_nxt_lyr)+(grv_nxt_lyr*hgt_nxt_lyr) ! [J kg-1] Dry static energy
        
     !cp_vpr_rcp_cp_dry_m1=spc_heat_H2O_vpr/spc_heat_dry_air-1.0 ! 0.87 [frc] for computing moist specific heatIrG81 pp. 77
     spc_heat_mst_air=spc_heat_dry_air*(1.0+(cp_vpr_rcp_cp_dry_m1*r_H2O_sfc))

     nrg_mst(idx)=(spc_heat_mst_air*tpt_pcl_nxt_lyr)+(grv_nxt_lyr*hgt_nxt_lyr)+ &
          ((pcl_H2O_svm_get(tpt_pcl_nxt_lyr,prs_nxt_lyr))*spc_ltn_heat(tpt_pcl_nxt_lyr)) ! [J kg-1] Moist static energy

        svp_H2O_nxt_lyr=svp_H2O_get(tpt_pcl_nxt_lyr)      ! [Pa] Saturation Vapor Pressure

     if (ppr_H2O_sfc==svp_H2O_nxt_lyr) then
        lev_LCL=hgt_nxt_lyr
        tpt_LCL=tpt_pcl_nxt_lyr
        LCL=.TRUE.
        
        if(tpt_LCL>tpt_nvr_nxt_lyr)then
           LFC=.TRUE.  ! [flg] Level of free convection
        end if
     elseif (ppr_H2O_sfc>svp_H2O_nxt_lyr) then
        LCL=.TRUE.
        
        call point_of_int(hgt,ppr_H2O_sfc,svp_H2O_prv_lyr,hgt_prv_lyr, &
             ppr_H2O_sfc,svp_H2O_nxt_lyr,lev_LCL,svp_LCL)
        
        tpt_LCL=tpt_pcl_prv_lyr-(pcl_MLR*(lev_LCL-hgt_prv_lyr))
        slope=hgt/(tpt_nvr_nxt_lyr-tpt_nvr_prv_lyr)
        b_int=hgt_nxt_lyr-(slope*tpt_nvr_nxt_lyr)
        tpt_nvr_LCL=(lev_LCL-b_int)/slope
        
        if(tpt_LCL>tpt_nvr_LCL)then
           LFC=.TRUE.  ! [flg] Level of free convection
        end if
        
        ! Recompute pcl tpt based on new LCL hgt
        slope=hgt/(prs_nxt_lyr-prs_prv_lyr)
        b_int=hgt_nxt_lyr-(slope*prs_nxt_lyr)
        prs_LCL=(lev_LCL-b_int)/slope
        
        call lift_pcl_PLR(tpt_LCL,prs_LCL,grv_prv_lyr,hgt_nxt_lyr-lev_LCL,PLR,tpt_pcl_nxt_lyr)
        pcl_PLR(idx)=PLR
        !Calculate Dry and Moist Static Energy
        nrg_dry(idx)=(spc_heat_dry_air*tpt_pcl_nxt_lyr)+(grv_nxt_lyr*hgt_nxt_lyr) ! [J kg-1] Dry static energy
        
        !cp_vpr_rcp_cp_dry_m1=spc_heat_H2O_vpr/spc_heat_dry_air-1.0 ! 0.87 [frc] for computing moist specific heatIrG81 pp. 77
        spc_heat_mst_air=spc_heat_dry_air*(1.0+cp_vpr_rcp_cp_dry_m1* &
             (pcl_H2O_svm_get(tpt_pcl_nxt_lyr,prs_nxt_lyr)))
        nrg_mst(idx)=(spc_heat_mst_air*tpt_pcl_nxt_lyr)+(grv_nxt_lyr*hgt_nxt_lyr)+ &
             ((pcl_H2O_svm_get(tpt_pcl_nxt_lyr,prs_nxt_lyr))* &
             spc_ltn_heat(tpt_pcl_nxt_lyr)) ! [J kg-1] Moist static energy
     end if ! End if ppr_H2O_sfc>=svp_H2O_nxt_lyr
  end if ! end if LCL
  
  ! LFC and LNB Check
  if (.NOT. LFC .AND. tpt_pcl_prv_lyr<=tpt_nvr_prv_lyr .AND. tpt_pcl_nxt_lyr<=tpt_nvr_nxt_lyr) then
     CINE_dltz=hgt
     CINE_fnc_val_prv=((tpt_pcl_prv_lyr-tpt_nvr_prv_lyr)/tpt_nvr_prv_lyr)*grv_prv_lyr
     CINE_fnc_val_nxt=((tpt_pcl_nxt_lyr-tpt_nvr_nxt_lyr)/tpt_nvr_nxt_lyr)*grv_nxt_lyr
     CINE=CINE+trap_rule(CINE_fnc_val_nxt,CINE_fnc_val_prv,CINE_dltz)
     
  else if (.NOT. LFC .AND. tpt_pcl_prv_lyr<=tpt_nvr_prv_lyr .AND. tpt_pcl_nxt_lyr>tpt_nvr_nxt_lyr) then
     LFC=.TRUE.
     call point_of_int(hgt,tpt_pcl_prv_lyr,tpt_nvr_prv_lyr,hgt_prv_lyr,tpt_pcl_nxt_lyr, &
          tpt_nvr_nxt_lyr,lev_LFC,tpt_LFC)
     
     CINE_dltz=lev_LFC-hgt_prv_lyr
     CINE_fnc_val_prv=((tpt_pcl_prv_lyr-tpt_nvr_prv_lyr)/tpt_nvr_prv_lyr)*grv_prv_lyr
     CINE_fnc_val_nxt=0.0
     CINE=CINE+trap_rule(CINE_fnc_val_nxt,CINE_fnc_val_prv,CINE_dltz)
     
     CAPE_dltz=hgt_nxt_lyr-lev_LFC
     CAPE_fnc_val_prv=0.0
     CAPE_fnc_val_nxt=((tpt_pcl_nxt_lyr-tpt_nvr_nxt_lyr)/tpt_nvr_nxt_lyr)*grv_nxt_lyr
     CAPE=CAPE+trap_rule(CAPE_fnc_val_nxt,CAPE_fnc_val_prv,CAPE_dltz)
  else if (LFC .AND. .NOT. LNB .AND. tpt_pcl_nxt_lyr>tpt_nvr_nxt_lyr ) then
     CAPE_dltz=hgt
     CAPE_fnc_val_prv=((tpt_pcl_prv_lyr-tpt_nvr_prv_lyr)/tpt_nvr_prv_lyr)*grv_prv_lyr
     CAPE_fnc_val_nxt=((tpt_pcl_nxt_lyr-tpt_nvr_nxt_lyr)/tpt_nvr_nxt_lyr)*grv_nxt_lyr
     CAPE=CAPE+trap_rule(CAPE_fnc_val_nxt,CAPE_fnc_val_prv,CAPE_dltz)
  else if (LFC .AND. .NOT. LNB .AND. tpt_pcl_prv_lyr>=tpt_nvr_prv_lyr .AND. &
       tpt_pcl_nxt_lyr<tpt_nvr_nxt_lyr) then
     LNB=.TRUE.

     call point_of_int(hgt,tpt_pcl_prv_lyr,tpt_nvr_prv_lyr,hgt_prv_lyr,tpt_pcl_nxt_lyr, &
          tpt_nvr_nxt_lyr,lev_LNB,tpt_LNB)
 
     CAPE_dltz=lev_LNB-hgt_prv_lyr
     CAPE_fnc_val_prv=((tpt_pcl_prv_lyr-tpt_nvr_prv_lyr)/tpt_nvr_prv_lyr)*grv_prv_lyr
     CAPE_fnc_val_nxt=0.0
     CAPE=CAPE+trap_rule(CAPE_fnc_val_nxt,CAPE_fnc_val_prv,CAPE_dltz)
  end if ! End if .NOT. LFC .AND. tpt_pcl_prv_lyr==tpt_nvr_prv_lyr

end subroutine lift_pcl

subroutine lift_pcl_PLR (tpt_pcl_prv_lyr,prs_prv,grv,dlt_hgt,PLR,tpt_pcl_nxt_lyr)
  implicit none
  
  real,intent(in)::dlt_hgt
  real,intent(in)::grv
  real,intent(in)::tpt_pcl_prv_lyr
  real,intent(in)::prs_prv
  real,intent(out)::PLR
  real,intent(out)::tpt_pcl_nxt_lyr
  
  real drsdt_H2O        ! [K-1] Change in Mixing Ratio per Change in Temperature
  real svp_H2O_prv_lyr  ! [Pa] Saturation vapor pressure at previous layer
  real pcl_DLR          ! [K m-1] Midlayer dry adiabatic lapse rate
  real pcl_H2O_svm      ! [kg kg-1] Saturation mixing ratio of parcel
  real tpt_pcl_vrt      ! [K] Virtual temperature of the parcel
      
  pcl_H2O_svm=pcl_H2O_svm_get(tpt_pcl_prv_lyr,prs_prv) ! [kg kg-1] Saturation mixing ratio of parcel
  
  tpt_pcl_vrt=(1.0+eps_H2O_rcp_m1*pcl_H2O_svm)*tpt_pcl_prv_lyr ! [K] Virtual temperature of parcel
  
  svp_H2O_prv_lyr=svp_H2O_get(tpt_pcl_vrt) ! [Pa] Saturation Vapor Pressure
  
  pcl_DLR=grv/spc_heat_dry_air ! [K m-1] Dry lapse rate
  
  ! Pseudo Adiabatic Lapse Rate At Each Level
  drsdt_H2O=dsvpdt_H2O(tpt_pcl_vrt)*(gas_cst_dry_air/(gas_cst_H2O*(prs_prv-svp_H2O_prv_lyr))) ! [K-1]
  PLR=pcl_DLR/(1.0+((spc_ltn_heat(tpt_pcl_vrt)/spc_heat_dry_air)*drsdt_H2O)) ! [K m-1] Midlayer Pseudo-Adiabatic Lapse Rate
  
  tpt_pcl_nxt_lyr=tpt_pcl_prv_lyr-(PLR*dlt_hgt) ! [K] Parcel temperature after lifting

end subroutine lift_pcl_PLR

real function pcl_H2O_svm_get(tpt,pcl_prs)
  implicit none
  
  real,intent(in)::tpt     ! [K] Temperature
  real,intent(in)::pcl_prs ! [Pa] Pressure
  real svp_H2O_prv_lyr     ! [Pa] Saturation vapor pressure at previous midlayer

  svp_H2O_prv_lyr=svp_H2O_get(tpt) ! [Pa] Saturation Vapor Pressure

  pcl_H2O_svm_get=(eps_H2O*(svp_H2O_prv_lyr/(pcl_prs-svp_H2O_prv_lyr))) ! [kg kg-1] Saturation water vapor mixing ratio
  
  return
end function pcl_H2O_svm_get

subroutine point_of_int(dlt_hgt,x1_prv_lyr,x2_prv_lyr,hgt_prv_lyr,x1_nxt_lyr,x2_nxt_lyr,y_int_pnt,x_int_pnt)
  implicit none
  
  ! In Variables
  real,intent(in)::dlt_hgt ! [m] Change in height
  real,intent(in)::hgt_prv_lyr ! [m] Height intercept or vertical axis
  real,intent(in)::x1_nxt_lyr ! Next layer value for variable 1
  real,intent(in)::x2_nxt_lyr ! Next layer value for variable 2
  real,intent(in)::x1_prv_lyr ! Previous layer value for variable 1
  real,intent(in)::x2_prv_lyr ! Previous layer value for variable 2
  
  ! Out variables
  real,intent(out)::x_int_pnt ! x-value of intersection
  real,intent(out)::y_int_pnt ! y-value of intersection
  
  ! Static variables
  real b_var_1 ! Intercept of variable 1
  real b_var_2 ! Intercept of variable 2
  real slope_var_1 ! Slope of variable 1
  real slope_var_2 ! Slope of variable 2
  
  if (x1_nxt_lyr==x1_prv_lyr) then
     slope_var_2=dlt_hgt/(x2_nxt_lyr-x2_prv_lyr)
     b_var_2=hgt_prv_lyr-(slope_var_2*x2_prv_lyr)
     x_int_pnt=x1_nxt_lyr
     y_int_pnt=b_var_2+(slope_var_2*x_int_pnt)
  else
     ! Compute Slope
     slope_var_1=dlt_hgt/(x1_nxt_lyr-x1_prv_lyr)
     slope_var_2=dlt_hgt/(x2_nxt_lyr-x2_prv_lyr)
     
     ! Calculate b=y-mx
     b_var_1=hgt_prv_lyr-(slope_var_1*x1_prv_lyr)
     b_var_2=hgt_prv_lyr-(slope_var_2*x2_prv_lyr)
     
     ! Calculate intersect points
     x_int_pnt=(b_var_1-b_var_2)/(slope_var_2-slope_var_1)
     y_int_pnt=b_var_1+(slope_var_1*x_int_pnt)
  end if ! end if x1_nxt_lyr==x1_prv_lyr
end subroutine point_of_int

real function dsvpdt_H2O(tpt) 
  ! [Pa K-1] Derivative of saturation vapor pressure over planar ice or liquid water
  ! Purpose: "Fast scalar" derivative of saturation vapor pressure over planar ice water
  ! Temperature conversion, and included files
  ! Input must be temperature in Celsius in range -50 C < T < 50 C
  ! Output is derivative of saturation vapor pressure over planar ice water in Pa K-1
  ! Taken from Lowe and Ficke (1974) as reported in PrK78 p. 625 */
  implicit none

  real,intent(in)::tpt ! [K]
  real tpt_cls ! [C]
  real,dimension(7)::cff_lqd(0:6)= &
       (/4.438099984e-01,2.857002636e-02,7.938054040e-04,1.215215065e-05, &
       1.036561403e-07,3.532421810e-10,-7.090244804e-13/)
  real,dimension(7)::cff_ice(0:6)= &
       (/5.030305237e-01,3.773255020e-02,1.267995369e-03,2.477563108e-05, &
       3.005693132e-07,2.158542548e-09,7.131097725e-12/)
  
  tpt_cls = tpt-tpt_frz_pnt ! [K] --> [C]
  
  if (tpt_cls<0.0) then
     dsvpdt_H2O=cff_ice(0)+tpt_cls*(cff_ice(1)+tpt_cls*(cff_ice(2)+tpt_cls*(cff_ice(3)+ &
          tpt_cls*(cff_ice(4)+tpt_cls*(cff_ice(5)+cff_ice(6)*tpt_cls))))) ! [mb K-1]
  else
     dsvpdt_H2O=cff_lqd(0)+tpt_cls*(cff_lqd(1)+tpt_cls*(cff_lqd(2)+tpt_cls*(cff_lqd(3)+ &
          tpt_cls*(cff_lqd(4)+tpt_cls*(cff_lqd(5)+cff_lqd(6)*tpt_cls))))) ! [mb K-1]
  end if  ! End if tpt_cls<0.0
  dsvpdt_H2O=dsvpdt_H2O*100.0 ![mb K-1] -> [Pa K-1]
  return
end function dsvpdt_H2O

real function trap_rule(fnc_val_nxt,fnc_val_prv,dlt_z)
  ! Performs a trapezoidal numerical method to integrate
  implicit none

  real,intent(in)::fnc_val_nxt ! [K] Next Level Calc
  real,intent(in)::fnc_val_prv ! [m s-2] Previous Level Calc
  real,intent(in)::dlt_z ! [K] Temperature of Next Level
 
  trap_rule=(0.5*dlt_z*(fnc_val_prv+fnc_val_nxt))
  
  return
end function trap_rule

real function spc_ltn_heat(tpt) ! [J kg-1] Latent heat
  ! This function returns the latent heat for water
  !The latent heat is a function of temperature t [K]. 
  !the formulas are polynomial approximations to the values 
  !in table 92, p. 343 of the smithsonian meteorological tables, 
  !sixth revised edition, 1963 by roland list. The approximations 
  !were developed by eric smith at colorado state university.
  implicit none
  
  real,intent(in)::tpt
  !Evaporation/Condensation
  real,parameter::a0=3337118.5
  real,parameter::a1=-3642.8583
  real,parameter::a2=2.1263947
  !Melting/Freezing
  real,parameter::b0=-1161004.0
  real,parameter::b1=9002.2648
  real,parameter::b2=-12.931292
  !Sublimation/Deposition
  real,parameter::c0=2632536.8
  real,parameter::c1=1726.9659
  real,parameter::c2=-3.6248111
  
  if(tpt>273.15)then
     spc_ltn_heat=a0+a1*tpt+a2*tpt*tpt
  elseif(tpt<=273.15 .AND. tpt>271.0)then
     spc_ltn_heat=c0+c1*tpt+c2*tpt*tpt
     spc_ltn_heat=spc_ltn_heat+(b0+b1*tpt+b2*tpt*tpt)
  else
     spc_ltn_heat=c0+c1*tpt+c2*tpt*tpt
  end if ! end if tpt>273.15
   
  return
end function spc_ltn_heat

real function tpt_vrt_get(tpt,r_H2O)
  implicit none
  
  real,intent(in)::tpt ! [K] Temperature
  real,intent(in)::r_H2O ! [kg kg-1] water vapor mixing ratio
  
  tpt_vrt_get=(1.0+eps_H2O_rcp_m1*r_H2O)*tpt
  return
end function tpt_vrt_get

real function svp_H2O_get(tpt)
  implicit none
  
  real,intent(in)::tpt ! [K] Temperature
  
  if (tpt<tpt_frz_pnt) then
     svp_H2O_get=svp_H2O_ice_PrK78(tpt) ! [Pa] Saturation Vapor Pressure
  else
     svp_H2O_get=svp_H2O_lqd_PrK78(tpt) ! [Pa] Saturation Vapor Pressure
  end if !End if (tpt<tpt_frz_pnt)
  
  return
end function svp_H2O_get

end module cape_mdl
