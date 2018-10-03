! $Id$

! Purpose: Dimers and collision complexes

! Copyright (C) 1994--2018 Charlie Zender
! This software is distributed under the terms of the GNU General Public License
! See http://www.gnu.org/copyleft/gpl.html for full license text

! These routines are heavily used by clm

! Usage:
! use dmr_mdl ! [mdl] Dimers, collision complexes

module dmr_mdl ! [mdl] Dimers, collision complexes
  implicit none
  public ! [stt] Symbols are public unless individually qualified as private
  
contains 
  
  subroutine cnc_H2OH2O_Chy97(lvl_nbr,alt_dlt,tpt,RH_lqd,mpl_H2O,cnc_H2OH2O)
    ! Purpose:
    ! Compute H2O dimer concentrations from Petr Chylek (personal communication, 1997) formulations
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use phys_cst_mdl,only:tpt_frz_pnt,Avagadro,mmw_H2O ! [mdl] Fundamental and derived physical constants
    implicit none
    ! Parameters
    integer,parameter::lvl_nbr_max=200
    real,parameter::a=1.0e19 ! Scale of xsx formula
    real,parameter::b=5.62e-2 ! Linear temperature dependence coefficient in xsx formula
    real,parameter::c=1.75e-3 ! Quadratic temperature dependence coefficient in xsx formula
    ! Input Arguments
    integer lvl_nbr           ! dimension size
    real RH_lqd(lvl_nbr)      ! fraction
    real mpl_H2O(lvl_nbr)     ! [kg m-2]
    real tpt(lvl_nbr)         ! [K]
    real alt_dlt(lvl_nbr)     ! m
    ! Input/Output Arguments
    ! Output Arguments
    real cnc_H2OH2O(lvl_nbr) ! mlc/m3
    ! Local workspace
    integer idx
    real cnc_H2O(lvl_nbr_max) !  mlc/m3
    real k_H2O_H2O(lvl_nbr_max) ! equilibrium rate constant for H2O + H2O <-> H2O - H2O
    real tpt_cls              ! [C]
    real mpl_H2O_gxcm2        ! g/cm2
    real npl_H2OH2O(lvl_nbr_max) !  mlc m-2
    real npl_H2OH2O_cm2       ! mlc/cm2
    ! Main code
    
    do idx=1,lvl_nbr
       ! Parameterization of Chy97 was derived for -15.0 C < T
       tpt_cls=tpt(idx)-tpt_frz_pnt
       ! Anything colder than -15.0 C is treated as -15.0 C
       if (tpt_cls < -15.0) tpt_cls=-15.0
       mpl_H2O_gxcm2=mpl_H2O(idx)*1000.0/10000.0 ! [kg m-2] -> g/cm2
       npl_H2OH2O_cm2=a*(1.0+b*tpt_cls+c*tpt_cls*tpt_cls)*mpl_H2O_gxcm2*RH_lqd(idx) ! mlc/cm2
       npl_H2OH2O(idx)=npl_H2OH2O_cm2*1.0e4 ! mlc/cm2 -> mlc m-2
       cnc_H2OH2O(idx)=npl_H2OH2O(idx)/alt_dlt(idx) ! mlc/m3
       ! Diagnostic output
       cnc_H2O(idx)=mpl_H2O(idx)*Avagadro/(mmw_H2O*alt_dlt(idx)) ! [kg m-2] -> mlc/m3
       k_H2O_H2O(idx)=cnc_H2O(idx)*(cnc_H2O(idx)/cnc_H2OH2O(idx)) ! [m3 mlc-1]
    enddo                     ! end loop over lvl
    
    if (dbg_lvl == 5) then
       write (6,'(a,i3)') 'Method of Chy97:'
       write (6,'(a3,1x,6(a11,1x))') 'idx','t','mpl_H2O','RH_lqd','cnc_H2O','cnc_H2OH2O','k_H2O_H2O'
       write (6,'(a3,1x,6(a11,1x))') '','K','[kg m-2]','%','mlc/m3','mlc/m3','m3/mlc'
       do idx=1,lvl_nbr
          write (6,'(i3,1x,3(f11.7,1x),3(es11.4,1x))')  &
               idx,tpt(idx),mpl_H2O(idx),RH_lqd(idx)*100.0,cnc_H2O(idx),cnc_H2OH2O(idx),k_H2O_H2O(idx)
       enddo                  ! end loop over lvl
    endif                     ! end if dbg
    
    return
  end subroutine cnc_H2OH2O_Chy97                       ! end cnc_H2OH2O_Chy97()
  
  subroutine cnc_H2OH2O_CFT99(lvl_nbr,alt_dlt,tpt,RH_lqd,mpl_mst_air,mpl_H2O,cnc_H2OH2O) 
    ! Purpose:
    ! Compute H2O dimer concentrations from method of CFT99.
    use phys_cst_mdl,only:mmw_H2OH2O,mmw_H2O,Avagadro ! [mdl] Fundamental and derived physical constants
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    implicit none
    ! Parameters
    integer,parameter::lvl_nbr_max=200
    integer,parameter::pzn_nbr=1 ! Which CFT99 pzn to use?
    ! Input Arguments
    integer lvl_nbr           ! dimension size
    real RH_lqd(lvl_nbr)      ! fraction
    real mpl_H2O(lvl_nbr)     ! [kg m-2]
    real mpl_mst_air(lvl_nbr) ! [kg m-2]
    real tpt(lvl_nbr)         ! [K]
    real alt_dlt(lvl_nbr)     ! m
    ! Input/Output Arguments
    ! Output Arguments
    real cnc_H2OH2O(lvl_nbr)  ! mlc/m3
    ! Local workspace
    integer idx
    real cnc_H2O(lvl_nbr_max) !  mlc m-3
    real mpl_H2OH2O(lvl_nbr_max) ! kg m-2
    real q_H2OH2O(lvl_nbr_max) ! kg kg-1
    real k_H2O_H2O(lvl_nbr_max) ! equilibrium rate constant for H2O + H2O <-> H2O - H2O
    real a_pzn_CFT99          ! Parameter "a" in CFT99 eqn 2 for mmr of H2OH2O
    real b_pzn_CFT99          ! Parameter "b" in CFT99 eqn 2 for mmr of H2OH2O
    real mmr_H2OH2O_rcp_mmr_H2O ! Quantity "C" predicted by CFT99 eqn 2 pzn
    ! Commons
    ! Main code
    
    ! Set parameters for CFT99 eqn 2 for mass of H2OH2O relative to mass of H2O
    ! CFT99 parameterizations range from 1--4, with pzn 1 predicting fewest dimers, pzn 2 being the best guess, and pzn 4 predicting most dimers
    if (pzn_nbr == 1) then
       a_pzn_CFT99=144.0      ! Parameter "a" for mmr of H2OH2O from CFT99 eqn 2, pzn 1 (Munoz-Caro)
       b_pzn_CFT99=3535.0     ! Parameter "b" for mmr of H2OH2O from CFT99 eqn 2, pzn 1 (Munoz-Caro)
    else if (pzn_nbr == 2) then
       a_pzn_CFT99=811.0      ! Parameter "a" for mmr of H2OH2O from CFT99 eqn 2, pzn 2 (Experimental)
       b_pzn_CFT99=3800.0     ! Parameter "b" for mmr of H2OH2O from CFT99 eqn 2, pzn 2 (Experimental)
    else if (pzn_nbr == 3) then
       a_pzn_CFT99=337.0      ! Parameter "a" for mmr of H2OH2O from CFT99 eqn 2, pzn 3 (Slanina, ChG97)
       b_pzn_CFT99=3418.0     ! Parameter "b" for mmr of H2OH2O from CFT99 eqn 2, pzn 3 (Slanina, ChG97)
    else if (pzn_nbr == 4) then
       a_pzn_CFT99=395.       ! Parameter "a" for mmr of H2OH2O from CFT99 eqn 2, pzn 4 (Munoz-Caro scaled)
       b_pzn_CFT99=3422.      ! Parameter "b" for mmr of H2OH2O from CFT99 eqn 2, pzn 4 (Munoz-Caro scaled)
    endif                     ! endif pzn
    
    do idx=1,lvl_nbr
       mmr_H2OH2O_rcp_mmr_H2O=a_pzn_CFT99*exp(-b_pzn_CFT99/tpt(idx)) ! fraction
       ! NB: Multiplying by RH_lqd corrects CFT99 pzn for unsaturated atmospheres
       mpl_H2OH2O(idx)=RH_lqd(idx)*mpl_H2O(idx)*mmr_H2OH2O_rcp_mmr_H2O ! kg m-2
       q_H2OH2O(idx)=mpl_H2OH2O(idx)/mpl_mst_air(idx) ! kg kg-1
       q_H2OH2O(idx)=q_H2OH2O(idx)+0.0 ! kg kg-1 CEWI
       cnc_H2OH2O(idx)=Avagadro*mpl_H2OH2O(idx)/(alt_dlt(idx)*mmw_H2OH2O) ! mlc m-3
       ! Diagnostic output
       cnc_H2O(idx)=mpl_H2O(idx)*Avagadro/(mmw_H2O*alt_dlt(idx)) ! [kg m-2] -> mlc/m3
       ! 20000913: RH_lqd = 0.0 causes SIGFPE since cnc_H2OH2O will also be 0.0
       if (RH_lqd(idx) == 0.0) then
          k_H2O_H2O(idx)=0.0 ! [m3 mlc-1]
       else
          k_H2O_H2O(idx)=cnc_H2O(idx)*(cnc_H2O(idx)/cnc_H2OH2O(idx)) ! [m3 mlc-1]
       endif                  ! endif
    enddo                     ! end loop over lvl
    
    if (dbg_lvl == 5) then
       write (6,'(a,i3)') 'Method of CFT99:'
       write (6,'(a3,1x,6(a11,1x))') 'idx','t','mpl_H2O','RH_lqd','cnc_H2O','cnc_H2OH2O','k_H2O_H2O'
       write (6,'(a3,1x,6(a11,1x))') '','K','[kg m-2]','%','mlc/m3','mlc/m3','m3/mlc'
       do idx=1,lvl_nbr
          write (6,'(i3,1x,3(f11.7,1x),3(es11.4,1x))')  &
               idx,tpt(idx),mpl_H2O(idx),RH_lqd(idx)*100.0,cnc_H2O(idx),cnc_H2OH2O(idx),k_H2O_H2O(idx)
       enddo                  ! end loop over lvl
    endif                     ! end if dbg
    
    return
  end subroutine cnc_H2OH2O_CFT99                       ! end cnc_H2OH2O_CFT99()
  
  subroutine cnc_H2OH2O_Mun97(lvl_nbr,alt_dlt,q_H2O,prs,tpt,RH_lqd,mpl_mst_air,cnc_H2OH2O)
    ! Purpose:
    ! Compute H2O dimer concentrations from Mun97 formulations
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use phys_cst_mdl,only:Avagadro,mmw_H2O,mmw_H2OH2O,mmw_dry_air,prs_STP,tpt_frz_pnt ! [mdl] Fundamental and derived physical constants
    implicit none
    ! Parameters
    integer,parameter::lvl_nbr_max=200
    ! Input Arguments
    integer lvl_nbr           ! dimension size
    real RH_lqd(lvl_nbr)      ! fraction
    real mpl_mst_air(lvl_nbr) ! [kg m-2]
    real q_H2O(lvl_nbr)       ! kg/kg
    real prs(lvl_nbr)         ! Pa
    real tpt(lvl_nbr)         ! [K]
    real alt_dlt(lvl_nbr)     ! m
    ! Input/Output Arguments
    ! Output Arguments
    real cnc_H2OH2O(lvl_nbr) ! mlc/m3
    ! Local workspace
    integer idx
    real atm_H2O              ! atm
    real atm_H2OH2O           ! atm
    real cnc_H2O(lvl_nbr_max) ! mlc/m3
    real k_H2O_H2O(lvl_nbr_max) ! equilibrium rate constant for H2O + H2O <-> H2O - H2O (Mun97 p. 4133)
    real mpl_H2O(lvl_nbr_max) ! [kg m-2]
    real nrm_Mun97_CFB79      ! Normalize model parameterization of Mun97 to observations of CFB79
    real q_H2OH2O(lvl_nbr_max) ! kg H2O/kg air
    real tpt_tmp              ! [K]
    real vmr_H2O(lvl_nbr_max) ! mlc H2O/mlc air
    real vmr_H2OH2O(lvl_nbr_max) ! mlc H2O/mlc air
    ! Main code
    
    ! Initialize some constants
    ! Ratio of CFB79 experimental results at 373 K to Mun97 model results (Mun97 p. 4133)
    ! Experiment observed more dimers (higher K_H2O_H2O) than model
    nrm_Mun97_CFB79=0.0110/0.0034 ! 3.235 atm-1/atm-1 
    
    do idx=1,lvl_nbr
       ! Parameterization of Mun97 was derived for 273 < T < 373 K
       tpt_tmp=tpt(idx)
       ! Anything colder than 0 C is treated as 0 C
       if (tpt(idx) < tpt_frz_pnt) tpt_tmp=tpt_frz_pnt
       k_H2O_H2O(idx)=exp(-10.27+1595.29/tpt_tmp+16201558.49/tpt_tmp**3) ! atm-1 (Mun97 p. 4133)
       k_H2O_H2O(idx)=k_H2O_H2O(idx)*nrm_Mun97_CFB79 ! atm-1
       !     k_H2O_H2O(idx)=k_H2O_H2O(idx)* ! atm-1 -> m3 mlc-1
       vmr_H2O(idx)=q_H2O(idx)*mmw_dry_air/mmw_H2O
       atm_H2O=vmr_H2O(idx)*prs(idx)/prs_STP ! atm
       atm_H2OH2O=atm_H2O*(atm_H2O*k_H2O_H2O(idx)) ! atm
       vmr_H2OH2O(idx)=atm_H2OH2O*prs_STP/prs(idx)
       q_H2OH2O(idx)=vmr_H2OH2O(idx)*mmw_H2OH2O/mmw_dry_air
       ! Conversion from q -> cnc worked out on CSZ V p. 62
       cnc_H2OH2O(idx)=Avagadro*q_H2OH2O(idx)*mpl_mst_air(idx)/(alt_dlt(idx)*mmw_H2OH2O)
       ! Diagnostic output
       cnc_H2O(idx)=Avagadro*q_H2O(idx)*mpl_mst_air(idx)/(alt_dlt(idx)*mmw_H2O)
       mpl_H2O(idx)=cnc_H2O(idx)*alt_dlt(idx)*mmw_H2O/Avagadro
       k_H2O_H2O(idx)=cnc_H2O(idx)*(cnc_H2O(idx)/cnc_H2OH2O(idx)) ! [m3 mlc-1]
    enddo
    
    if (dbg_lvl == 5) then
       write (6,'(a,i3)') 'Method of Mun97:'
       write (6,'(a3,1x,6(a11,1x))') 'idx','t','mpl_H2O','RH_lqd','cnc_H2O','cnc_H2OH2O','k_H2O_H2O'
       write (6,'(a3,1x,6(a11,1x))') '','K','[kg m-2]','%','mlc/m3','mlc/m3','m3/mlc'
       do idx=1,lvl_nbr
          write (6,'(i3,1x,3(f11.7,1x),3(es11.4,1x))')  &
               idx,tpt(idx),mpl_H2O(idx),RH_lqd(idx)*100.0,cnc_H2O(idx),cnc_H2OH2O(idx),k_H2O_H2O(idx)
       enddo                  ! end loop over lvl
    endif                     ! end if dbg
    
    return
  end subroutine cnc_H2OH2O_Mun97                       ! end cnc_H2OH2O_Mun97()
  
  subroutine cnc_H2OH2O_Zen97( &
       RH_lqd, &
       cnc_H2OH2O, &
       cnc_mst_air, &
       alt_dlt, &
       lvl_nbr, &
       mpl_CWP, &
       mpl_mst_air, &
       prs, &
       q_H2O, &
       tpt &
       )
    ! Purpose:
    ! Compute H2O dimer concentrations from Zen97 formulations
    ! Shawn Kathman at PNL <shawn.kathmann@pnl.gov> is interested in 
    ! working on ab initio water structure studies and can compute
    ! water N-mer abundances given environmental conditions.
    use sng_mdl,only:ftn_strlen   ! [mdl] String manipulation
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use phys_cst_mdl,only:Joules_per_calorie,Avagadro,Boltzmann,mmw_H2O,mmw_H2OH2O,eps_H2O ! [mdl] Fundamental and derived physical constants
    implicit none
    ! Parameters
    integer,parameter::lvl_nbr_max=200
    ! Input Arguments
    integer lvl_nbr           ! dimension size
    real RH_lqd(lvl_nbr)      ! fraction
    real cnc_mst_air(lvl_nbr) ! mlc/m3
    real mpl_CWP(lvl_nbr)     ! [kg m-2]
    real mpl_mst_air(lvl_nbr) ! [kg m-2]
    real q_H2O(lvl_nbr)       ! kg/kg
    real prs(lvl_nbr)         ! Pa
    real tpt(lvl_nbr)         ! [K]
    real alt_dlt(lvl_nbr)     ! m
    ! Input/Output Arguments
    real cnc_H2OH2O(lvl_nbr) ! mlc/m3
    ! Output Arguments
    ! Local workspace
    real(selected_real_kind(p=12))::pi
    real(selected_real_kind(p=12))::wgt_Bzm_dmr(lvl_nbr_max) ! Boltzmann weight of dimer evaporation
    real(selected_real_kind(p=12))::wgt_Bzm_mnm(lvl_nbr_max) ! Boltzmann weight of monomer evaporation
    real(selected_real_kind(p=12))::wgt_Bzm_nrm(lvl_nbr_max) ! Normalized Boltzmann weight of evaporation
    integer idx
    real CWC(lvl_nbr_max)     ! Condensed water content
    real MFP_H2OH2O(lvl_nbr_max) ! Mean free path (between collisions) of H2OH2O
    real area_spc             ! Specific surface area (area per unit mass) of hydrometeor distribution 
    real cnc_H2O(lvl_nbr_max) ! Concentration of H2O vapor
    real cnc_H2OH2O_CSZ(lvl_nbr_max) ! Equilibrium concentration of H2O dimers from condensate source
    real cnd_cff_H2O          ! Condensation coefficient for H2O vapor -> liquid (fraction of impinging H2O molecules sticking to surface)
    real dmr_H20_dmt          ! Diameter of H2O dimers
    real evp_fsh_H2OH2O(lvl_nbr_max) ! Evaporation efficiency of H2O dimers relative to monomers
    real flx_H2O_nvr_sfc(lvl_nbr_max) ! Flux of vapor from environment to unit hydrometeor surface
    real flx_H2O_sfc_nvr(lvl_nbr_max) ! Flux of vapor from unit hydrometeor surface to environment
    real flx_H2OH2O_sfc_nvr(lvl_nbr_max) ! Flux of H2O dimers from unit hydrometeor surface to environment
    real frq_cll_H2OH2O(lvl_nbr_max) ! Collision rate of an H2O dimer with other molecules
    real heat_ltn_vpr         ! Latent heat of evaporation
    real nrg_bnd_H2OH2O       ! Binding energy of H2O dimers
    real nrg_bnd_H2OH2O_kcm   ! Binding energy of H2O dimers in kcal/mol
    real nrg_evp_dmr          ! Energy required to evaporate a dimer from hydrometeor surface
    real nrg_evp_mnm          ! Energy required to evaporate a monomer from hydrometeor surface
    real nrg_knt_avg(lvl_nbr_max) ! Mean kinetic (translational) energy of molecules
    real nrg_rat(lvl_nbr_max) ! Ratio of dimer binding energy to kinetic energy
    real pp_H2O(lvl_nbr_max)  ! Partial pressure of H2O
    real r_H2O(lvl_nbr_max)   ! Dry-mass mixing ratio (r) of H2O 
    ! real rds_fct=10.0e-6     ! [m] Effective radius of hydrometeor distribution
    real sfc_area(lvl_nbr_max) ! Surface area of hydrometeors per unit volume
    !real sigma=1.4                ! Standard deviation of hydrometeor distribution
    !real spd_MB_RMS_H2O(lvl_nbr_max) ! RMS speed of H2O in Maxwellian velocity distribution
    real spd_MB_RMS_H2OH2O(lvl_nbr_max) ! RMS speed of H2O dimers in Maxwellian velocity distribution
    real spd_MB_avg_H2O(lvl_nbr_max) ! Mean speed of H2O in Maxwellian velocity distribution
    real spd_MB_avg_H2OH2O(lvl_nbr_max) ! Mean speed of H2O dimers in Maxwellian velocity distribution
    real src_H2OH2O(lvl_nbr_max) ! Source (production rate per unit volume) of H2O dimers
    real tau_cll_H2OH2O(lvl_nbr_max) ! Mean time between collisions for an H2O dimer
    real time_H2OH2O(lvl_nbr_max) ! Mean lifetime of H2O dimers
    ! Main code
    
    ! Initialize some constants
    pi=4.D0*atan(1.D0)
    area_spc=77.5             ! m2/kg
    heat_ltn_vpr=2.485e+06    ! J/kg (at 0 C) NB: generalize this as L_v = f(T)
    
    ! Use standard sticking efficiency
    cnd_cff_H2O=0.035         ! fraction PrK78 p. 133 eqn 5-53
    ! Experimental uncertainty in binding energy is +/- 0.7 kcal/mol, or ~15% CFB79 p. 2710
    nrg_bnd_H2OH2O_kcm=5.44  ! kcal/mol, Mun97 p. 4131, DaC97 p. 8152, CFB79 p. 2710
    nrg_bnd_H2OH2O=nrg_bnd_H2OH2O_kcm*1000.0*Joules_per_calorie/Avagadro ! kcal/mol -> J/mlc
    
    ! Assume energy required to liberate a monomer is (1/Avagadro)*(latent heat of vaporization of one mole)
    nrg_evp_mnm=heat_ltn_vpr*mmw_H2O/Avagadro ! J
    ! Assume energy required to liberate a dimer is twice as much as a monomer
    nrg_evp_dmr=2.0*nrg_evp_mnm ! J
    
    ! Compute necessary thermodynamic fields
    do idx=1,lvl_nbr
       CWC(idx)=mpl_CWP(idx)/alt_dlt(idx) ! kg/m3
       r_H2O(idx)=q_H2O(idx)/(1.0-q_H2O(idx))
       pp_H2O(idx)=prs(idx)*r_H2O(idx)/(eps_H2O+r_H2O(idx))
       pp_H2O(idx)=0.0+pp_H2O(idx) ! CEWI
    enddo                     ! end loop over lvl
    
    ! Assume interaction cross section of H2O dimers is geometric diameter for MFP purposes
    dmr_H20_dmt=3.0e-10        ! m Dac97 p. 8152 computes most favorable O-O distance is 2.87 A
    
    ! Compute necessary kinetic fields
    do idx=1,lvl_nbr
       spd_MB_avg_H2O(idx)=sqrt(8.0*Boltzmann*tpt(idx)*Avagadro/(pi*mmw_H2O)) ! m/s PrK78 p. 133 eqn. 5-50
       spd_MB_avg_H2OH2O(idx)=sqrt(8.0*Boltzmann*tpt(idx)*Avagadro/(pi*mmw_H2OH2O)) ! m/s PrK78 p. 133 eqn. 5-50
       !spd_MB_RMS_H2O(idx)=sqrt(3.0*Boltzmann*tpt(idx)*Avagadro/mmw_H2O) ! m/s KiK80 p. 393
       spd_MB_RMS_H2OH2O(idx)=sqrt(3.0*Boltzmann*tpt(idx)*Avagadro/mmw_H2OH2O) ! m/s KiK80 p. 393
       MFP_H2OH2O(idx)=1.0/(cnc_mst_air(idx)*pi*dmr_H20_dmt*dmr_H20_dmt) ! m KiK80 p. 395 eqn. 13
       nrg_knt_avg(idx)=1.5*Boltzmann*tpt(idx) ! J
       nrg_rat(idx)=nrg_bnd_H2OH2O/nrg_knt_avg(idx)
       frq_cll_H2OH2O(idx)=spd_MB_RMS_H2OH2O(idx)/MFP_H2OH2O(idx) ! s-1 KiK80 p. 397 eqn 16b
       tau_cll_H2OH2O(idx)=1.0/frq_cll_H2OH2O(idx) ! s
       ! Dimers survive until they collide with more energy than dimer binding energy 
       ! Assume each dimer survives a # of collisions proportional to ratio of binding energy to kinetic energy
       !     time_H2OH2O(idx)=nrg_rat(idx)*tau_cll_H2OH2O(idx) ! s
       time_H2OH2O(idx)=60.0 ! s
    enddo                     ! end loop over lvl
    
    ! Compute parameters controlling surface dimer flux
    do idx=1,lvl_nbr
       ! Assume dimer evaporation efficiency is Boltzmann weighted 
       wgt_Bzm_mnm(idx)=exp(-nrg_evp_mnm/(Boltzmann*tpt(idx)))
       wgt_Bzm_dmr(idx)=exp(-nrg_evp_dmr/(Boltzmann*tpt(idx)))
       !     wgt_Bzm_nrm(idx)=1./(exp(1.)+exp(2.))
       wgt_Bzm_nrm(idx)=1./(wgt_Bzm_mnm(idx)+wgt_Bzm_dmr(idx))
       evp_fsh_H2OH2O(idx)=wgt_Bzm_nrm(idx)*wgt_Bzm_dmr(idx) ! fraction (mlc H2OH2O per mlc H2O evaporated)
       ! Override with constant evaporation efficiency
       evp_fsh_H2OH2O(idx)=0.001
       cnc_H2O(idx)=Avagadro*q_H2O(idx)*mpl_mst_air(idx)/(alt_dlt(idx)*mmw_H2O)
       flx_H2O_nvr_sfc(idx)=cnc_H2O(idx)*spd_MB_avg_H2O(idx)/4.0 ! mlc m-2 s-1 PrK78 p. 133 eqn. 5-49
       ! In steady state equilibrium, condensation equals evaporation
       flx_H2O_sfc_nvr(idx)=flx_H2O_nvr_sfc(idx)*cnd_cff_H2O ! mlc m-2 s-1
       flx_H2OH2O_sfc_nvr(idx)=flx_H2O_sfc_nvr(idx)*evp_fsh_H2OH2O(idx) ! mlc m-2 s-1
       
       ! Apply these rates to droplet surface area
       sfc_area(idx)=CWC(idx)*area_spc ! m2/m3
       src_H2OH2O(idx)=flx_H2OH2O_sfc_nvr(idx)*sfc_area(idx) ! mlc m-3 s-1
       cnc_H2OH2O_CSZ(idx)=src_H2OH2O(idx)*time_H2OH2O(idx) ! mlc m-3
    enddo                     ! end loop over lvl
    
    if (dbg_lvl == 5) then
       write (6,'(a,i3)') 'Method of Zen97:'
       write (6,'(a,i3)') 'ThermoKinetic Properties:'
       write (6,'(a3,1x,6(a11,1x))') 'idx','nrg_bnd_dmr','nrg_knt_avg','nrg_rat','spd_avg_H2O','spd_avg_dmr','spd_RMS_dmr'
       write (6,'(a3,1x,6(a11,1x))') '','J','J','','m/s','m/s','m/s'
       do idx=1,lvl_nbr
          write (6,'(i3,1x,6(es11.4,1x))')  &
               idx, &
               nrg_bnd_H2OH2O,nrg_knt_avg(idx),nrg_rat(idx), &
               spd_MB_avg_H2O(idx),spd_MB_avg_H2OH2O(idx),spd_MB_RMS_H2OH2O(idx)
       enddo                  ! end loop over lvl
    endif                     ! end if dbg
    
    if (dbg_lvl == 5) then
       write (6,'(a,i3)') 'Collision Statistics:'
       write (6,'(a3,1x,4(a11,1x))') 'idx','MFP_dmr','frq_cll_dmr','tau_dmr','time_dmr'
       write (6,'(a3,1x,4(a11,1x))') '','m','mlc/s','s','s'
       do idx=1,lvl_nbr
          write (6,'(i3,1x,4(es11.4,1x))')  &
               idx, &
               MFP_H2OH2O(idx),frq_cll_H2OH2O(idx),tau_cll_H2OH2O(idx), &
               time_H2OH2O(idx)
       enddo                  ! end loop over lvl
    endif                     ! end if dbg
    
    if (dbg_lvl == 5) then
       write (6,'(a,i3)') 'Surface fluxes:'
       write (6,'(a3,1x,6(a11,1x))') 'idx', &
            'evp_fsh','H2O_nvr_sfc','H2O_sfc_nvr', &
            'dmr_sfc_nvr','src_H2OH2O','cnc_H2OH2O_CSZ'
       write (6,'(a3,1x,6(a11,1x))') '','mlc/mlc','mlc/m2/s','mlc/m2/s','mlc/m2/s','mlc/m3/s','mlc/m3'
       do idx=1,lvl_nbr
          write (6,'(i3,1x,6(es11.4,1x))')  &
               idx, &
               evp_fsh_H2OH2O(idx),flx_H2O_nvr_sfc(idx),flx_H2O_sfc_nvr(idx), &
               flx_H2OH2O_sfc_nvr(idx),src_H2OH2O(idx),cnc_H2OH2O_CSZ(idx)
       enddo                  ! end loop over lvl
    endif                     ! end if dbg
    
    if (dbg_lvl == 5) then
       write (6,'(a,i3)') 'Summary:'
       write (6,'(a3,1x,6(a11,1x))') 'idx','t','mpl_CWP','RH_lqd','cnc_H2O','cnc_H2OH2O','cnc_CSZ'
       write (6,'(a3,1x,6(a11,1x))') '','K','g/m2','%','mlc/m3','mlc/m3','mlc/m3'
       do idx=1,lvl_nbr
          write (6,'(i3,1x,6(es11.4,1x))')  &
               idx, &
               tpt(idx),mpl_CWP(idx)*1000.0,RH_lqd(idx)*100.0, &
               cnc_H2O(idx),cnc_H2OH2O(idx),cnc_H2OH2O_CSZ(idx)
       enddo                  ! end loop over lvl
    endif                     ! end if dbg
    
    ! Combine condensate equilibrium concentration with the gas-phase equilibrium concentration of dimers
    if (dbg_lvl == 6.and..true.) then
       do idx=1,lvl_nbr
          cnc_H2OH2O(idx)=cnc_H2OH2O(idx)+cnc_H2OH2O_CSZ(idx) ! mlc m-3
       enddo                  ! end loop over lvl
    endif                     ! endif Zen97
    
    if (dbg_lvl == 7.and..true.) then
       write (6,'(a,a)') prg_nm(1:ftn_strlen(prg_nm)), &
            ': WARNING Artificially enhancing in-cloud dimer concentration by factor of 5'
       do idx=1,lvl_nbr
          if (CWC(idx) > 0.0) cnc_H2OH2O(idx)=5.0*cnc_H2OH2O(idx) ! mlc m-3
       enddo                  ! end loop over lvl
    endif                     ! endif Zen97
    
    return
  end subroutine cnc_H2OH2O_Zen97                       ! end cnc_H2OH2O_Zen97()
  
end module dmr_mdl ! [mdl] Dimers, collision complexes
