c     $Id$

      program mie2nc

c     Purpose: Turn BPB's Mie property files into netCDF files
c     Note: BPB's Mie files are superceded by mie program as of 20010615

c     Compilation:
c     cd $HOME/aca; make -W mie2nc.F OPTS=D mie2nc; cd -
c     cd $HOME/aca; make -W mie2nc.F mie2nc; cd -
c     cd $HOME/aca; make OPTS=D mie2nc; cd -

c     Usage:
c     For liquid water spheres (the default) you can just use mie2nc.
c     Other phases and substances requires explicitly setting # of bands
c     mie2nc -n 73 -i ${DATA}/aca/lqd_05.dat -o ${DATA}/aca/lqd_05.nc
c     mie2nc -n 74 -i ${DATA}/aca/lqd_07.dat -o ${DATA}/aca/lqd_07.nc
c     mie2nc -n 74 -i ${DATA}/aca/lqd_10.dat -o ${DATA}/aca/lqd_10.nc
c     mie2nc -n 194 -i ${DATA}/aca/ice_20.dat -o ${DATA}/aca/ice_20.nc
c     mie2nc -n 194 -i ${DATA}/aca/ice_05.dat -o ${DATA}/aca/ice_05.nc

c     Currently the code processes input ASCII data files that look like:

c     re=10micron  log normal sig=1.6  water droplet mie properties
c     wvlngth re      qe       qs      we       ge      sg sct       sg ext     psi(m2/g) 
c     0.225 10.0 2.049779 2.049763 0.999992 0.899230  3.3193E-06  3.3193E-06  1.5373E-01
c     0.250 10.0 2.059366 2.059375 0.999999 0.896990  3.3348E-06  3.3348E-06  1.5445E-01
c     0.275 10.0 2.059684 2.059706 0.999999 0.895371  3.3354E-06  3.3353E-06  1.5448E-01

      use flx_slr_mdl ! [mdl] Solar spectral fluxes
      use netcdf ! [mdl] netCDF interface
      use nf90_utl ! [mdl] netCDF utilities
      use sng_mdl ! [mdl] String manipulation
      use utl_mdl ! [mdl] Utility functions (date_time_get,mnt_chk...)
      implicit none
c     Parameters
      character*(*) CVS_Id
      integer fl_in_unit
      integer wvl_nbr_max
      parameter(fl_in_unit=73,
     $     wvl_nbr_max=200,
     $     CVS_Id='$Id$')
      integer,parameter::bnd_nbr_max=wvl_nbr_max

c     Commons
#include <parameter.com> /* Fundamental and defined physical constants */ 
#include <dbg.com> /* Debugging constants */
c     Input Arguments
c     Input/Output Arguments
c     Output Arguments
c     Local workspace
      character argv*80
      character cmd_ln*200
      character lbl_sng*80
      character fl_in*80
      character fl_out*80
      character*26::lcl_date_time ! Time formatted as Day Mth DD HH:MM:SS TZ YYYY
      character prg_ID*200
      character CVS_Date*28
      character CVS_Revision*16

      integer arg
      integer exit_status       ! program exit status
      integer arg_nbr
      integer rcd               ! return success code

      real float_foo

      integer wvl_dmn_id        ! dimension ID for bands
      integer grd_dmn_id        ! dimension ID for grid
      integer wvl_idx           ! counting index
      integer nc_id             ! file handle
      integer wvl_nbr           ! dimension size

      integer abs_cff_mss_id
      integer abs_fsh_id
      integer abs_xsx_id
      integer asm_prm_id
      integer bandwidth_id
      integer bnd_id
      integer ext_cff_mss_id
      integer ext_fsh_id
      integer ext_xsx_id
      integer frc_slr_flx_LaN68_id
      integer frc_slr_flx_ThD71_id
      integer rds_ffc_id
      integer sca_cff_mss_id
      integer sca_fsh_id
      integer sca_xsx_id
      integer gsd_anl_id
      integer ss_alb_id
      integer ss_co_alb_id
      integer wvl_id            ! coordinate ID
      integer wvl_ctr_id
      integer wvl_grd_id
      integer wvl_max_id
      integer wvl_min_id
      integer wvl_dlt_id

      real abs_cff_mss(wvl_nbr_max)
      real abs_fsh(wvl_nbr_max)
      real abs_xsx(wvl_nbr_max)
      real asm_prm(wvl_nbr_max)
      real bandwidth
      real bnd(wvl_nbr_max)
      real wvl(wvl_nbr_max)     ! coordinate variable
      real ext_cff_mss(wvl_nbr_max)
      real ext_fsh(wvl_nbr_max)
      real ext_xsx(wvl_nbr_max)
      real frc_slr_flx_LaN68(wvl_nbr_max)
      real frc_slr_flx_ThD71(wvl_nbr_max)
      real rds_ffc
      real sca_cff_mss(wvl_nbr_max)
      real sca_fsh(wvl_nbr_max)
      real sca_xsx(wvl_nbr_max)
      real gsd_anl
      real ss_alb(wvl_nbr_max)
      real ss_co_alb(wvl_nbr_max)
      real wvl_ctr(wvl_nbr_max)
      real wvl_grd(wvl_nbr_max+1)
      real wvl_max(wvl_nbr_max)
      real wvl_min(wvl_nbr_max)
      real wvl_dlt(wvl_nbr_max)

c     Local variables
      real mss_nrm
      real dns_hydrometeor
      real vlm_nrm

c     Main code

c     Initialize default values
      dbg_lvl=0
      dns_hydrometeor=1000.0 ! [kg m-3]
      exit_status=0
      float_foo=0
      fl_in='/data/zender/aca/lqd_10.dat'
      wvl_nbr=74
      CVS_Date='$Date$'
      CVS_Revision='$Revision$'
      fl_out='/data/zender/aca/lqd_10.nc'
      rcd=nf90_noerr              ! nf90_noerr == 0

c     Retrieve command line arguments
      call date_time_get(lcl_date_time)
      call ftn_cmd_ln_sng(cmd_ln)
      call ftn_prg_ID_mk(CVS_Id,CVS_Revision,CVS_Date,prg_ID)
      write (6,'(a)') prg_ID(1:ftn_strlen(prg_ID))
      arg_nbr=command_argument_count()
      do arg=1,arg_nbr
         call getarg(arg,argv)
         if (argv(1:2) == '-D') then
            call getarg(arg+1,argv)
            read(argv,'(i4)') dbg_lvl
         endif
         if (argv(1:2) == '-f') then
            call getarg(arg+1,argv)
            read(argv,'(f8.3)') float_foo
         endif
         if (argv(1:2) == '-i') then
            call getarg(arg+1,argv)
            read(argv,'(a)') fl_in
         endif
         if (argv(1:2) == '-n') then
            call getarg(arg+1,argv)
            read(argv,'(i4)') wvl_nbr
         endif
         if (argv(1:2) == '-o') then
            call getarg(arg+1,argv)
            read(argv,'(a)') fl_out
         endif
         if (argv(1:2) == '-v') then
            write(6,'(a)') CVS_Id
            goto 1000
         endif
      end do

c     Compute any quantities that might depend on command line input
      call ftn_strnul(fl_in)
      call ftn_strnul(fl_out)

      open (fl_in_unit,file=fl_in,status='old',iostat=rcd)

      read (fl_in_unit,'(a4,f2.0,a23,f3.1)') lbl_sng,rds_ffc,lbl_sng,gsd_anl
      read (fl_in_unit,'(a80)') lbl_sng
      do wvl_idx=1,wvl_nbr
         read (fl_in_unit,*) 
     $        wvl_ctr(wvl_idx),
     $        rds_ffc,
     $        ext_fsh(wvl_idx),
     $        sca_fsh(wvl_idx),
     $        ss_alb(wvl_idx),
     $        asm_prm(wvl_idx),
     $        sca_xsx(wvl_idx),
     $        ext_xsx(wvl_idx),
     $        ext_cff_mss(wvl_idx)

      enddo

      close (fl_in_unit)
      write (6,'(a20,1x,a)') 'Read input data from',fl_in(1:ftn_strlen(fl_in))

c     Convert the input data to SI units where necessary.
      do wvl_idx=1,wvl_nbr
         wvl_ctr(wvl_idx)=wvl_ctr(wvl_idx)*1.0e-6 ! um -> m
         sca_xsx(wvl_idx)=sca_xsx(wvl_idx)*1.0e-4 ! cm2 -> m2
         ext_xsx(wvl_idx)=ext_xsx(wvl_idx)*1.0e-4 ! cm2 -> m2
         ext_cff_mss(wvl_idx)=ext_cff_mss(wvl_idx)*1000.0 ! [m2 g-1] -> m2 kg-1
      enddo
      gsd_anl=gsd_anl*1.0e-6
      rds_ffc=rds_ffc*1.0e-6

c     Compute diagnostic variables  
      bandwidth=wvl(2)-wvl(1)
      do wvl_idx=1,wvl_nbr
         wvl(wvl_idx)=wvl_ctr(wvl_idx)
         bnd(wvl_idx)=wvl_ctr(wvl_idx)
      enddo
      do wvl_idx=2,wvl_nbr
         wvl_min(wvl_idx)=0.5*(wvl_ctr(wvl_idx-1)+wvl_ctr(wvl_idx))
      enddo
      wvl_min(1)=wvl_min(2)-(wvl_min(3)-wvl_min(2))
      do wvl_idx=1,wvl_nbr-1
         wvl_max(wvl_idx)=0.5*(wvl_ctr(wvl_idx)+wvl_ctr(wvl_idx+1))
      enddo
      wvl_max(wvl_nbr)=wvl_max(wvl_nbr-1)+
     $     (wvl_max(wvl_nbr-1)-wvl_max(wvl_nbr-2))
      do wvl_idx=1,wvl_nbr
         wvl_dlt(wvl_idx)=wvl_max(wvl_idx)-wvl_min(wvl_idx)
         wvl_grd(wvl_idx)=wvl_min(wvl_idx)
      enddo
      wvl_grd(wvl_nbr+1)=wvl_max(wvl_nbr)

c     A word to posterity about the input BPB provides: 
c     The cross-sections are per unit concentration. 
c     The cross-sections were computed by
c     integrating the Mie optical efficiencies times the hydrometeor
c     cross sectional areas times the distribution over all radii.
c     These integrals are proportional to the second moment of the
c     radius.
c     To convert the cross sections to optical depth,
c     multiply them by the absolute concentration N_tot of hydrometeors 
c     (#/m^3) and then by the path length (m). N_tot is the total
c     number of hydrometeors of all sizes in a given volume, a number
c     most likely to be obtained from observations. However, the
c     number concentration of hydrometeors isn't a very useful 
c     quantity to work with in GCMs, so we convert to mass 
c     absorption/extinction/scattering coefficients instead.
c     These mass coefficients are proportional to the third moment 
c     of the radius and so contain information not in the cross 
c     sections.
c     BPB supplies the mass extinction coefficient. To convert this
c     to optical depth, simply multiply by the total path of 
c     hydrometeors (kg/m^2) in the homogeneous layer. I extract
c     and save the mass absorption/scattering coefficients from
c     the other input variables.

c     mss_nrm is the mass of H2O per unit concentration,
c     For instance, if the number density ( = concentration) of
c     the hydrometeor distribution were 100 per cubic meter, then
c     the H2O mass would be 100*mss_nrm.

c     vlm_nrm is analogous to mss_nrm, i.e., it is the volume of
c     H2O per unit concentration. vlm_nrm is originally obtained by
c     computing the integral of the volume times the normalized
c     distribution over all radii, i.e., vlm_nrm is 4*pi/3 times
c     the third moment of the hydrometeor radius.

c     dns_hydrometeor is the H2O density BPB used to convert
c     H2O volume to H2O mass.

      mss_nrm=ext_xsx(1)/ext_cff_mss(1)
      vlm_nrm=mss_nrm/dns_hydrometeor

      do wvl_idx=1,wvl_nbr
         abs_fsh(wvl_idx)=ext_fsh(wvl_idx)-sca_fsh(wvl_idx)
         abs_xsx(wvl_idx)=ext_xsx(wvl_idx)-sca_xsx(wvl_idx)
         ss_co_alb(wvl_idx)=1.0-ss_alb(wvl_idx)
      enddo
      do wvl_idx=1,wvl_nbr
         if (abs_fsh(wvl_idx) < 0.0) then
            abs_fsh(wvl_idx)=0.0
         endif
         if (abs_xsx(wvl_idx) < 0.0) then
            abs_xsx(wvl_idx)=0.0
         endif
         if (ss_co_alb(wvl_idx) < 0.0) then
            ss_co_alb(wvl_idx)=0.0
         endif
      enddo
      do wvl_idx=1,wvl_nbr
         abs_cff_mss(wvl_idx)=abs_xsx(wvl_idx)/mss_nrm
         sca_cff_mss(wvl_idx)=sca_xsx(wvl_idx)/mss_nrm
         frc_slr_flx_ThD71(wvl_idx)= ! NB: arguments in microns
     $        slfftd(1.0e6*wvl_min(wvl_idx),1.0e6*wvl_max(wvl_idx))
         frc_slr_flx_LaN68(wvl_idx)= ! NB: arguments in microns
     $        slffln(1.0e6*wvl_min(wvl_idx),1.0e6*wvl_max(wvl_idx))
      enddo

c     Begin NetCDF output routines
      rcd=rcd+nf90_create(fl_out,nf90_clobber,nc_id)
c     Define dimension IDs
      rcd=rcd+nf90_def_dim(nc_id,'wvl',wvl_nbr,wvl_dmn_id)
      if (rcd /= nf90_noerr) call nf90_err_exit(rcd,fl_out)
      rcd=rcd+nf90_def_dim(nc_id,'grd',wvl_nbr+1,grd_dmn_id)
      if (rcd /= nf90_noerr) call nf90_err_exit(rcd,fl_out)
c     Variable definitions
      rcd=rcd+nf90_def_var(nc_id,'abs_cff_mss',nf90_float,wvl_dmn_id,abs_cff_mss_id)
      rcd=rcd+nf90_def_var(nc_id,'abs_fsh',nf90_float,wvl_dmn_id,abs_fsh_id)
      rcd=rcd+nf90_def_var(nc_id,'abs_xsx',nf90_float,wvl_dmn_id,abs_xsx_id)
      rcd=rcd+nf90_def_var(nc_id,'asm_prm',nf90_float,wvl_dmn_id,asm_prm_id)
      rcd=rcd+nf90_def_var(nc_id,'bandwidth',nf90_float,bandwidth_id)
      rcd=rcd+nf90_def_var(nc_id,'bnd',nf90_float,wvl_dmn_id,bnd_id)
      rcd=rcd+nf90_def_var(nc_id,'ext_cff_mss',nf90_float,wvl_dmn_id,ext_cff_mss_id)
      rcd=rcd+nf90_def_var(nc_id,'ext_fsh',nf90_float,wvl_dmn_id,ext_fsh_id)
      rcd=rcd+nf90_def_var(nc_id,'ext_xsx',nf90_float,wvl_dmn_id,ext_xsx_id)
      rcd=rcd+nf90_def_var(nc_id,'frc_slr_flx_LaN68',nf90_float,wvl_dmn_id,frc_slr_flx_LaN68_id)
      rcd=rcd+nf90_def_var(nc_id,'frc_slr_flx_ThD71',nf90_float,wvl_dmn_id,frc_slr_flx_ThD71_id)
      rcd=rcd+nf90_def_var(nc_id,'rds_ffc',nf90_float,rds_ffc_id)
      rcd=rcd+nf90_def_var(nc_id,'sca_cff_mss',nf90_float,wvl_dmn_id,sca_cff_mss_id)
      rcd=rcd+nf90_def_var(nc_id,'sca_fsh',nf90_float,wvl_dmn_id,sca_fsh_id)
      rcd=rcd+nf90_def_var(nc_id,'sca_xsx',nf90_float,wvl_dmn_id,sca_xsx_id)
      rcd=rcd+nf90_def_var(nc_id,'gsd_anl',nf90_float,gsd_anl_id)
      rcd=rcd+nf90_def_var(nc_id,'ss_alb',nf90_float,wvl_dmn_id,ss_alb_id)
      rcd=rcd+nf90_def_var(nc_id,'ss_co_alb',nf90_float,wvl_dmn_id,ss_co_alb_id)
      rcd=rcd+nf90_def_var(nc_id,'wvl',nf90_float,wvl_dmn_id,wvl_id)
      rcd=rcd+nf90_def_var(nc_id,'wvl_ctr',nf90_float,wvl_dmn_id,wvl_ctr_id)
      rcd=rcd+nf90_def_var(nc_id,'wvl_grd',nf90_float,grd_dmn_id,wvl_grd_id)
      rcd=rcd+nf90_def_var(nc_id,'wvl_max',nf90_float,wvl_dmn_id,wvl_max_id)
      rcd=rcd+nf90_def_var(nc_id,'wvl_min',nf90_float,wvl_dmn_id,wvl_min_id)
      rcd=rcd+nf90_def_var(nc_id,'wvl_dlt',nf90_float,wvl_dmn_id,wvl_dlt_id)

c     Add global attibutes
      rcd=rcd+nf90_put_att(nc_id,nf90_global,'CVS_Id',CVS_Id)
      rcd=rcd+nf90_put_att(nc_id,nf90_global,'creation_date',lcl_date_time)
      rcd=rcd+nf90_put_att(nc_id,nf90_global,'prg_ID',prg_ID(1:ftn_strlen(prg_ID)))
      rcd=rcd+nf90_put_att(nc_id,nf90_global,'cmd_ln',cmd_ln(1:ftn_strlen(cmd_ln)))

c     Add english text descriptions
      rcd=rcd+nf90_put_att(nc_id,abs_cff_mss_id,'long_name','Mass absorption coefficient of distribution')
      rcd=rcd+nf90_put_att(nc_id,abs_fsh_id,'long_name','Absorption efficiency')
      rcd=rcd+nf90_put_att(nc_id,abs_xsx_id,'long_name','Absorption cross section per unit concentration')
      rcd=rcd+nf90_put_att(nc_id,asm_prm_id,'long_name','Asymmetry parameter')
      rcd=rcd+nf90_put_att(nc_id,bandwidth_id,'long_name','Uniform width of bands')
      rcd=rcd+nf90_put_att(nc_id,bnd_id,'long_name','Band center')
      rcd=rcd+nf90_put_att(nc_id,ext_cff_mss_id,'long_name','Mass extinction coefficient of distribution')
      rcd=rcd+nf90_put_att(nc_id,ext_fsh_id,'long_name','Extinction efficiency')
      rcd=rcd+nf90_put_att(nc_id,ext_xsx_id,'long_name','Extinction cross section per unit concentration')
      rcd=rcd+nf90_put_att(nc_id,frc_slr_flx_LaN68_id,'long_name','Fraction of solar flux: Labs & Neckel 1968')
      rcd=rcd+nf90_put_att(nc_id,frc_slr_flx_ThD71_id,'long_name','Fraction of solar fluxL Thekeakara & Drummond 1971')
      rcd=rcd+nf90_put_att(nc_id,rds_ffc_id,'long_name','Effective radius')
      rcd=rcd+nf90_put_att(nc_id,sca_cff_mss_id,'long_name','Mass scattering coefficient of distribution')
      rcd=rcd+nf90_put_att(nc_id,sca_fsh_id,'long_name','Scattering efficiency')
      rcd=rcd+nf90_put_att(nc_id,sca_xsx_id,'long_name','Scattering cross section per unit concentration')
      rcd=rcd+nf90_put_att(nc_id,gsd_anl_id,'long_name','Geometric standard deviation of particle distribution')
      rcd=rcd+nf90_put_att(nc_id,ss_alb_id,'long_name','Single scattering albedo')
      rcd=rcd+nf90_put_att(nc_id,ss_co_alb_id,'long_name','Single scattering co-albedo')
      rcd=rcd+nf90_put_att(nc_id,wvl_id,'long_name','Wavelength at band center')
      rcd=rcd+nf90_put_att(nc_id,wvl_ctr_id,'long_name','Band center wavelength')
      rcd=rcd+nf90_put_att(nc_id,wvl_grd_id,'long_name','Wavelength grid')
      rcd=rcd+nf90_put_att(nc_id,wvl_max_id,'long_name','Band maximum wavelength')
      rcd=rcd+nf90_put_att(nc_id,wvl_min_id,'long_name','Band minimum wavelength')
      rcd=rcd+nf90_put_att(nc_id,wvl_dlt_id,'long_name','Bandwidth')

c     Add units
      rcd=rcd+nf90_put_att(nc_id,abs_cff_mss_id,'units','meter2 kilogram-1')
      rcd=rcd+nf90_put_att(nc_id,abs_fsh_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,abs_xsx_id,'units','meter2')
      rcd=rcd+nf90_put_att(nc_id,asm_prm_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,bandwidth_id,'units','meter')
      rcd=rcd+nf90_put_att(nc_id,bnd_id,'units','meter')
      rcd=rcd+nf90_put_att(nc_id,ext_cff_mss_id,'units','meter2 kilogram-1')
      rcd=rcd+nf90_put_att(nc_id,ext_fsh_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,ext_xsx_id,'units','meter2')
      rcd=rcd+nf90_put_att(nc_id,frc_slr_flx_LaN68_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,frc_slr_flx_ThD71_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,rds_ffc_id,'units','meter')
      rcd=rcd+nf90_put_att(nc_id,sca_cff_mss_id,'units','meter2 kilogram-1')
      rcd=rcd+nf90_put_att(nc_id,sca_fsh_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,sca_xsx_id,'units','meter2')
      rcd=rcd+nf90_put_att(nc_id,gsd_anl_id,'units','meter')
      rcd=rcd+nf90_put_att(nc_id,ss_alb_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,ss_co_alb_id,'units','fraction')
      rcd=rcd+nf90_put_att(nc_id,wvl_id,'units','meter')
      rcd=rcd+nf90_put_att(nc_id,wvl_ctr_id,'units','meter')
      rcd=rcd+nf90_put_att(nc_id,wvl_grd_id,'units','meter')
      rcd=rcd+nf90_put_att(nc_id,wvl_max_id,'units','meter')
      rcd=rcd+nf90_put_att(nc_id,wvl_min_id,'units','meter')
      rcd=rcd+nf90_put_att(nc_id,wvl_dlt_id,'units','meter')

c     Now all dimensions, variables, and attributes have been defined, make call to end define mode
      rcd=rcd+nf90_enddef(nc_id)

c     Write data
      rcd=rcd+nf90_put_var(nc_id,abs_cff_mss_id,abs_cff_mss)
      rcd=rcd+nf90_put_var(nc_id,abs_fsh_id,abs_fsh)
      rcd=rcd+nf90_put_var(nc_id,abs_xsx_id,abs_xsx)
      rcd=rcd+nf90_put_var(nc_id,asm_prm_id,asm_prm)
      rcd=rcd+nf90_put_var(nc_id,bandwidth_id,bandwidth)
      rcd=rcd+nf90_put_var(nc_id,bnd_id,bnd)
      rcd=rcd+nf90_put_var(nc_id,ext_cff_mss_id,ext_cff_mss)
      rcd=rcd+nf90_put_var(nc_id,ext_fsh_id,ext_fsh)
      rcd=rcd+nf90_put_var(nc_id,ext_xsx_id,ext_xsx)
      rcd=rcd+nf90_put_var(nc_id,frc_slr_flx_LaN68_id,frc_slr_flx_LaN68)
      rcd=rcd+nf90_put_var(nc_id,frc_slr_flx_ThD71_id,frc_slr_flx_ThD71)
      rcd=rcd+nf90_put_var(nc_id,rds_ffc_id,rds_ffc)
      rcd=rcd+nf90_put_var(nc_id,sca_cff_mss_id,sca_cff_mss)
      rcd=rcd+nf90_put_var(nc_id,sca_fsh_id,sca_fsh)
      rcd=rcd+nf90_put_var(nc_id,sca_xsx_id,sca_xsx)
      rcd=rcd+nf90_put_var(nc_id,gsd_anl_id,gsd_anl)
      rcd=rcd+nf90_put_var(nc_id,ss_alb_id,ss_alb)
      rcd=rcd+nf90_put_var(nc_id,ss_co_alb_id,ss_co_alb)
      rcd=rcd+nf90_put_var(nc_id,wvl_id,wvl)
      rcd=rcd+nf90_put_var(nc_id,wvl_ctr_id,wvl_ctr)
      rcd=rcd+nf90_put_var(nc_id,wvl_grd_id,wvl_grd)
      rcd=rcd+nf90_put_var(nc_id,wvl_max_id,wvl_max)
      rcd=rcd+nf90_put_var(nc_id,wvl_min_id,wvl_min)
      rcd=rcd+nf90_put_var(nc_id,wvl_dlt_id,wvl_dlt)

      rcd=rcd+nf90_close(nc_id)
      write (6,'(a16,1x,a)') 'Wrote results to',fl_out(1:ftn_strlen(fl_out))
      if (rcd /= nf90_noerr) write (6,'(a,a,i4,a)') prg_nm,': ERROR rcd = ',rcd,' on exit'

 1000 continue

      call exit(exit_status)
      end




















