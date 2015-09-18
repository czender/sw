c
c $Id$
c
c $Log: not supported by cvs2svn $
c Revision 1.1.1.1  1998-09-15 02:06:41  zender
c Imported sources
c
c Revision 5.4  1995/06/19  04:54:58  zender
c have made many changes. got rid of fp_in,out,err. added flat updraft
c switch. changed -S -s to both refer to ice. got rid of while loop that
c became infinite when there was no cloud. changed vertical netCDF coord.
c to altitude. fixed bug in appending date and rev. to cmdline.
c
c Revision 5.3  1995/05/24  23:11:17  zender
c updated and validated the cloud model for ansi compatability,
c removed berkeley calls, streamlined architecture dependent
c fortran --> c calling methodology. validated results against
c comps II runs.
c
c Revision 5.2  1995/05/23  00:43:16  zender
c synchronization checkin prior to temperature sensitivity study
c for SPCP stable cloud mods.
c
c Revision 5.1  1994/02/14  21:28:54  zender
c Added DMC94 homogeneous ammonium sulphate solution nucleation,
c added -D 79 SWSCF/SWCF diagnostic. About to change default shape
c to hexagonal ice crystals, but will save bullet routines.
c NetCDF routines are not keeping up with the changes....
c
c Revision 5.0  1993/08/30  01:55:35  zender
c this version was used to make the paper (preprint at least).
c about to start working on an animation.
c
c Revision 4.12  1993/07/21  23:26:28  zender
c added Liou_IR_fudge() as a default routine whenever LIOU is true.
c
c Revision 4.11  1993/06/19  22:05:25  zender
c added ebert&curry parameterization to LW computations, and
c made this the default for the region outside of the cloud grid.
c
c Revision 4.10  1993/06/11  04:35:31  zender
c Played w/ Makefiles, made all CRAY2 refs. into plain CRAY,
c and put all shavano stuff on peace. changed surface type to
c land for predictable albedos, changed tropics SAZ to 30 deg.,
c removed gaussians from initial distributions. changed plev = 113.
c slingo option doesn't work for some reason (-D 53)
c
c Revision 4.9  1993/06/09  00:59:14  zender
c prepared the -Z initial_distribution_type switch to handle everyone's
c different cloud. Made Knollenberg fit his observations, prepared
c optical properties for 50 sizes and removed all num_size == 40 specific
c stuff. changed tau_CCN to 1. s.
c
c Revision 4.8  1993/05/29  00:01:13  zender
c fixed the Liou_fudge() bug, fixed a wind_speed bug, added
c Dowling and Radke's crystal distribution function as default.
c
c Revision 4.7  1993/05/27  16:08:24  zender
c A synchronization check-in for all cloud programs.
c
c Revision 4.6  1993/05/27  14:41:42  zender
c intermediate version to allow for a bug search. -X option now either calls
c Liou_fudge() (default) or Liou_interp(), but neither works well anymore.
c even the Mie option crashes in Tropics mode!
c
c Revision 4.4  1993/04/22  23:47:27  zender
c big change was integrated Liou ice crystal data with -X,
c and rearranging the spectral intervals in CCM2 colmod code
c so that CCM2 bands 14,15 are now my bands 15,16, and CCM2 16
c went to 14.  hex crosssections are still too big compared to
c corresponding bullet/spheres, tinkering with interp. geometry.
c
c Revision 4.3  1993/04/16  01:12:21  zender
c incorporated some Liou data, changed integration of optical
c parameters et al. to depend on optical cross sections, not
c efficiencies.  LIOU switch doesn't work, bug in bilinear_interp().
c
c Revision 4.2  1993/03/08  00:57:47  zender
c this version contains working working data i/o in 3 formats:
c ASCII, CZ_BINARY (my vanilla binary format), and NETCDF.
c since netCDF is transparent across platforms (i.e., cray and sun),
c the next version will only allow netCDF i/o.
c
c Revision 4.0  1993/02/18  01:49:55  zender
c this is pretty much the working comps version.  next additions
c will be major:  Liou's habit params, real aggregation, eddy
c diffusion.
c
c Revision 3.4  1993/01/23  03:55:36  zender
c Not much different.  Added start from scratch growth -D72, tuned
c iographs to make beautiful contour plots. About to add AGGREGATION.
c
c Revision 3.3  1993/01/08  02:19:41  zender
c cloud top, base now determined by Concentration not by IWC.
c Improved ESS selection, added IR optical depth routine, SW
c properties alb., trans., abs.  Switched to interface fluxes.
c Fixed Cray errors.
c
c Revision 3.2  1993/01/05  17:42:16  zender
c Implemented my own in line LW flux divergence code to extract
c crystal heating from.  Crashes on the Cray. About to tune the
c apportion routine and cloudiness param (and Eulerian growth??)
c so they work with non-monotonic growth rates vs. crystal size.
c
c Revision 3.1  1992/11/25  17:24:43  zender
c implemented routine to automatically use/create a disk based
c Mie parameter file whenver there are 40 crystal sizes. about
c to eliminate trace gas absorption within cloud grid to get crystal
c flux divergences.
c
c Revision 3.0  1992/11/16  02:41:01  zender
c size/spectral-weighted crystal rad. heating rates are implemented.
c old mean intensity method deadwood is stripped out.  must now
c compensate flux divergences for intrinsic trace gas heating.
c planck_avg_abs_eff not needed, just retained for a diagnostic.
c
      subroutine colmod(dayyr,rlat,pmid2,t2,h2ommr2,o3mmr2,cldfrc2,
     &     clwp2,ps2,ts,tg,oro,rghns,sndpth,alvss,alvsw,alnis,alniw,
     &     frcts,idebug,error_flag,
     &     envpint,
     &     bccmil,tccmil,numlay,ncclev,
     &     cbi,cti,ncldlay,
     &     geff,aeff,taueff,crecz,swfdiv,lwfdiv,
     &     fupsw,fdsw,fuplw,fdlw,fupcbi,fdncti,qrscld,qrlcld,qrcld,
     &     lwcf,swcf,albedo,swscf,lwmac)
c     
c ccm2 column radiation model on suns
c
c-----------------------------------------------------------------------
c
c     !   Basic grid point resolution parameters
c
      integer plev,plon,plat,plonp2,plevp,pcnst,plond,platd,plevd
      integer nxpt,jintmx,i1,j1,j1m
c
      parameter(plev=  133,        ! number of vertical levels
     $          plevp=plev+1,     ! number of vertical interfaces
     $          plon=1,           ! number of longitudes (T42)
     $          plat=1,           ! number of latitudes (T42)
     $          plonp2=plon+2,    ! length necessary for fft
     $          pcnst=1,          ! number of constituents
     $          nxpt=1,           !
     $          jintmx=2,         !
     $          plond=plon + 1 + 2*nxpt, ! slt extended domain longitude
     $          platd=plat + 2*nxpt + 2*jintmx, ! slt extended domain lat.
     $          plevd=plev*(3+pcnst), ! fold plev,pcnst indices into one
     $          i1=1+nxpt,        ! model starting longitude (3-d)
     $          j1=1+nxpt+jintmx, ! model starting latitude (3-d)
     $          j1m=j1-1)         ! model starting offset (3-d)
c
c-----------------------------------------------------------------------
C
C model control time variables
C
      integer
     $      nrstrt  ,nstep   ,nstepr  ,nestep  ,nstop   ,
     $      mdbase  ,msbase  ,mdcur   ,mscur   ,
     $      mbdate  ,mbsec   ,mcdate  ,mcsec   ,
     $      nndbas  ,nnsbas  ,nnbdat  ,nnbsec
      real  calday
c
      common /comtim/
     $      calday  ,nrstrt  ,nstep   ,nstepr  ,nestep  ,nstop   ,
     $      mdbase  ,msbase  ,mdcur   ,mscur   ,
     $      mbdate  ,mbsec   ,mcdate  ,mcsec   ,
     $      nndbas  ,nnsbas  ,nnbdat  ,nnbsec
C
c-----------------------------------------------------------------------
c
c local arguments
c
      integer nrow             ! latitude row index 
      real clat,               ! Current latitude (radians)
     $     oro(plond),         ! land/ocean/sea ice flag
     $     sndpth(plond),      ! snow depth (liquid water equivalent)
     $     ts(plond),          ! surface air temperature
     $     tg(plond),          ! surface (skin) temperature
     $     ps(plond),          ! surface pressure
     $     pmid(plond,plev),   ! model level pressures
     $     pint(plond,plevp),  ! model interface pressures
     $     pmln(plond,plev),   ! natural log of pmid
     $     piln(plond,plevp),  ! natural log of pint
     $     t(plond,plev),      ! model level temperatures
     $     h2ommr(plond,plev), ! model level specific humidity
     $     cldfrc(plond,plevp),! fractional cloud cover
     $     effcld(plond,plevp),! effective fractional cloud cover
     $     clwp(plond,plev),   ! cloud liquid water path
     $     plol(plond,plevp),  ! o3 pressure weighted path lengths (cm)
     $     plos(plond,plevp)   ! o3 path lengths (cm)
c
c output solar 
c
      real solin(plond),      ! solar incident flux
     $     sabtp(plond),      ! total column absorbed solar flux
     $     frsa(plond),       ! surface absorbed solar flux
     $     clrst(plond),      ! clr sky total column abs solar flux
     $     clrss(plond),      ! clr sky surface abs solar flux
     $     qrs(plond,plev)    ! solar heating rate
c
c output longwave
c
      real firtp(plond),      ! net outgoing lw flx at model top
     $     frla(plond),       ! srf longwave cooling (up-dwn) flux
     $     clrlt(plond),      ! clr sky lw flx at model top
     $     clrls(plond),      ! clr sky lw flx at srf (up-dwn)
     $     qrl(plond,plev),   ! longwave cooling rate
     $     slwd(plond)        ! srf down longwave flux
c
c surface radiative heating
c
      real srfrad(plond)       ! srf radiative heat flux
c
c local workspace
c
      real coszrs(plond),      ! cosine solar zenith angle
     $     loctim(plond)       ! local time of solar computation
c------------------------------------------------------------------------
c
c     zender cirrus cloud hack workspace
c
      real alvss,! grid box alb for vs over strng zn srfs
     $     alvsw, ! grid box alb for vs over weak  zn srfs
     $     alnis, ! grid box alb for ni over strng zn srfs
     $     alniw, ! grid box alb for ni over weak  zn srfs
     $     frcts  ! fraction of area in grid box strng zn 
c
      real albvss, ! grid box alb for vs over strng zn srfs
     $     albvsw, ! grid box alb for vs over weak  zn srfs
     $     albnis, ! grid box alb for ni over strng zn srfs
     $     albniw, ! grid box alb for ni over weak  zn srfs
     $     frctst  ! fraction of area in grid box strng zn 
c
      common /crdalb/albvss(plond,plat),albvsw(plond,plat),
     $               albnis(plond,plat),albniw(plond,plat),
     $               frctst(plond,plat)
c
cz      equivalence (alvss,albvss(1,1)),(alvsw,albvsw(1,1))
cz      equivalence (alniw,albniw(1,1)),(alnis,albnis(1,1))
c
      real vegtyp, ! surface thermal type, based on veg type
     $     rghnss, ! aerodynamic roughness length
     $     evapf , ! constant surface evaporability
     $     vevapf, ! variable surface evaporability
     $     snwjan, ! snow cover (liq water equiv) for january
     $     snwjly  ! snow cover (liq water equiv) for july
c
      common /crdsrf/vegtyp(plond,plat),rghnss(plond,plat),
     $               evapf (plond,plat),vevapf(plond,plat),
     $               snwjan(plond,plat),snwjly(plond,plat)
c
      real dayyr(plond),       ! day of year
     $     rlat(plond),        ! latitude input
     $     o3mmr(plond,plev),  ! o3 mass mixing ratio
     $     rghns,              ! aerodynamic roughness length
     $     ps2,                ! ground pressure
     $     lwcf,               ! longwave cloud forcing
     $     swcf,               ! shortwave cloud forcing
     $     albedo,             ! planetary (cloudy) albedo
     $     lwmac(0:numlay+1)   ! longwave mass absorption coeff (m^2/kg)
c
c     Create space to hold the swap arrays for inverting indice order
      real
     $     pmid2(plev),   ! model level pressures
     $     t2(plev),      ! model level temperatures
     $     h2ommr2(plev), ! model level specific humidity
     $     o3mmr2(plev),  ! o3 mass mixing ratio
     $     cldfrc2(plevp),! fractional cloud cover
     $     clwp2(plev)    ! cloud liquid (or ice) water path
c
      integer nspint  ! number of spectral intervals across solar spectrum
      parameter ( nspint = 18 )
c
      integer j,      !index for index trickery
     &     kcz,       !index for index trickery
     &     idebug,    !debugging level from cloud model
     &     ns,        !index over spectral intervals
     &     bccmil,    !bottom ccm2 interface level
     &     tccmil,    !top ccm2 interface level
     &     numlay,    !num_layer in cloud model
     &     ncclev,    !num_ccm2_level (with overlaps already subtracted)
c     ncclev is the # of ccm2 levels NOT in the cloud grid, therefore
c     ncclev+numlay MUST equal plev.
     &     cbi,       !cloud bottom index
     &     cti,       !cloud top index
     &     ncldlay    !number of layers in cloud itself
c
      real
     &     envpint(0:numlay+1),   ! interface pressures in cloud grid
     &     geff(0:nspint+1,0:numlay+1),   !effective asymmetry factor
     &     aeff(0:nspint+1,0:numlay+1),   !effective singel scat albedo
     &     taueff(0:nspint+1,0:numlay+1), !effective extinction optical depth
     &     crecz(0:numlay+1), !effective radius of particles (meters)
     &     swfdiv(0:nspint+1,0:numlay+1),  !SW spectral flux divergence
     &     lwfdiv(0:numlay+1), !LW flux divergence
     &     swscf,              ! shortwave surface cloud forcing
     &     fupsw(0:numlay+1), !total shortwave upward directed flux
     &     fdsw(0:numlay+1), !total shortwave downward directed flux
     &     fuplw(0:numlay+1), !total longwave upward directed flux
     &     fdlw(0:numlay+1), !total longwave downward directed flux
     &     fupcbi, !longwave upward directed flux at cloud base interface
     &     fdncti, !longwave downward directed flux at cloud top interface
     &     qrscld(0:numlay+1), !total shortwave local heating rate
     &     qrlcld(0:numlay+1), !total longwave local heating rate
     &     qrcld(0:numlay+1) !total local heating rate
c
c     end of zender workspace alterations
c-----------------------------------------------------------------------
c
      if( ncclev+numlay.ne.plev ) then
         print *,'Parameter plev needs to be reset and re-compiled for',
     $        ' this geometry.'
         print *,'plev = ',plev,' != ncclev + numlay = ',ncclev,' + ',
     $        numlay,' = ',ncclev+numlay
         print *,'Exiting program pronto.'
c     set the error_flag for exiting
         error_flag=-1
         return
      endif
c
c     Assign the common block variables to their command-line 
c     counterparts
      albvss(1,1)=alvss
      albvsw(1,1)=alvsw
      albnis(1,1)=alnis
      albniw(1,1)=alniw
      frctst(1,1)=frcts
      rghnss(1,1)=rghns
      ps(1)=ps2
c
      if( idebug.eq.44 ) then
         print *,'before rearranging . . .'
         print *,'level   p(Pa)    t(k)   h2ommr (g/g) o3mmr ',
     +        '(g/g) cld cvr  cld i/lwp(kg/m2)'
         do 45 k=1,plev
            write(6,99) k   ,pmid2(k),t2(k),h2ommr2(k),o3mmr2(k)
     +           ,cldfrc2(k),clwp2(k)
 45      continue
      endif	

c     Convert from the SI units of the cloud model to ccm2rad units
c     and reverse the indexical ordering
      do 50 i=1,1
         do 55 k=1,plev
            j=plev-k+1
            pmid(i,j)=pmid2(k)/100.  !convert Pa -> mb
            t(i,j)=t2(k)
            h2ommr(i,j)=h2ommr2(k)
            o3mmr(i,j)=o3mmr2(k)
            cldfrc(i,j)=cldfrc2(k)
            clwp(i,j)=clwp2(k)*1000. !convert kg/m^2 -> g/m^2
 55      continue
 50   continue

      if( idebug.eq.44 ) then
         print *,'after rearranging . . .'
         print *,'level   p(mb)    t(k)   h2ommr (g/g) o3mmr ',
     +        '(g/g) cld cvr  cld i/lwp(g/m2)'
         do 70 i=1,1
            do 75 k=1,plev
               write(6,99) k   ,pmid(i,k),t(i,k),h2ommr(i,k),o3mmr(i,k)
     +              ,cldfrc(i,k),clwp(i,k)
 75         continue
 70      continue
      endif

      if( idebug.eq.44 ) then
         print *,'input parameters . . .'
         do 76 k=1,numlay
            print *,'level spint geff aeff taueff crecz 
     &           lwmac swfdiv lwfdiv'
            do 77 ns=1,nspint
               write(6,98) k,ns,geff(ns,k),aeff(ns,k),
     &              taueff(ns,k),crecz(k),
     &              lwmac(k),swfdiv(ns,k),lwfdiv(k)
 77         continue
 76      continue
      endif
c
c     reverse the indexical ordering in altitude of the arrays
cz         do 86 k=1,numlay
cz         crecz(k)=crecz2(numlay+1-k)
cz            do 87 ns=1,nspint
cz               geff(ns,k)=geff2(ns,numlay+1-k)
cz               aeff(ns,k)=aeff2(ns,numlay+1-k)
cz               taueff(ns,k)=taueff2(ns,numlay+1-k)
cz 87         continue
cz 86      continue
c
 98   format(2x,i3,1x,i3,9(1pe11.4,1x))
 99   format(2x,i3,1x,6(1pe11.4,1x))
c
c------------------------------------------------------------------------
      if( idebug.eq.45.or.idebug.eq.53 ) then
         write(6,*) ' .... Begin CCM2 Radiation ....'
      endif
c     
c......... set parameters in common blocks:
c
      call stparm
c
      if( idebug.eq.45.or.idebug.eq.53 ) then
         write(6,*) ' .... Read in profile data ....'
      endif
c
      call getdat(gravx  ,cpairx ,epsilox,stebolx,   nrow,
     $            clat   ,oro    ,sndpth ,ts     ,     tg,
     $            ps     ,pmid   ,pint   ,pmln   ,   piln,
     $            t      ,h2ommr ,plol   ,plos   , cldfrc,
     $            clwp   ,effcld ,
     &     dayyr,rlat,o3mmr,
     &     lwmac,crecz,
     &     bccmil,tccmil,numlay,ncclev,
     &     cbi,cti,ncldlay,
     &     idebug)
c
      call radini(gravx,cpairx,epsilox,stebolx)
c
c compute radiation
c
c-----------------------------------------------------------------------
      call radctl(nrow   ,clat   ,oro    ,sndpth ,ts     ,
     $            tg     ,ps     ,pmid   ,pint   ,pmln   ,
     $            piln   ,t      ,h2ommr ,cldfrc ,effcld ,
     $            clwp   ,plol   ,plos   ,solin  ,sabtp  ,
     $            frsa   ,clrst  ,clrss  ,qrs    ,firtp  ,
     $            frla   ,clrlt  ,clrls  ,qrl    ,srfrad ,
     &     fupsw,fdsw,fuplw,fdlw,fupcbi,fdncti,
     &     geff,aeff,taueff,crecz,swfdiv,lwfdiv,
     &     bccmil,tccmil,numlay,ncclev,
     &     cbi,cti,ncldlay,
     &     idebug)
c-----------------------------------------------------------------------
c
      call cmpsol(clat,coszrs,eccf,loctim)
c
c write out final results:
c
      if( idebug.eq.45.or.idebug.eq.53 ) then
         write(6,*) ' ----- solar and longwave results ----- '
         write(6,*) ' ------------------ '
c     
         pie  = 4.*atan(1.)
cz    rlat = 180.*clat/pie
c     
         do 100 i=1,1
c     
            write(6,*) ' calendar day of year      = ',calday
            write(6,*) ' earth-sun distance factor = ',eccf
cz    write(6,*) ' earth latitude            = ',rlat
            write(6,*) ' earth latitude            = ',rlat(i)
            write(6,*) ' cosine solar zenith angle = ',coszrs(i)
            write(6,*) ' local solar time          = ',loctim(i)
            sol = solin(i)
            write(6,*) ' solar insolation          = ',sol,' wm-2' 
            sab = sabtp(i)
            if( sol .le. 0. ) then
               alb = -999.
            else
               alb = (sol-sab)/sol   
            endif
            write(6,*) ' solar total albedo        = ',alb
            write(6,*) ' solar total absorbed      = ',sab,' wm-2' 
            frs = frsa(i)
            frsatm = sab - frs
            write(6,*) ' solar absorbed atmosphere = ',frsatm,' wm-2' 
            write(6,*) ' solar absorbed surface    = ',frs,' wm-2' 
            clt = clrst(i)
            write(6,*) ' solar clear toa absorbed  = ',clt,' wm-2' 
            cls = clrss(i)
            clatm = clt - cls
            write(6,*) ' solar clear atm absorbed  = ',clatm,' wm-2' 
            write(6,*) ' solar clear srf absorbed  = ',cls,' wm-2' 
            scf = sab - clt
            write(6,*) ' solar cloud forcing       = ',scf,' wm-2' 
            if( sol .gt. 0. ) then
               albc = (sol - clt) / sol
               write(6,*) ' solar clear sky albedo    = ',albc
            endif
            write(6,*) ' ------------------ '
c     
            flw = firtp(i)
            write(6,*) ' longwave net up at top    = ',flw,' wm-2' 
            fla = frla (i)
            write(6,*) ' longwave surface net up   = ',fla,' wm-2' 
            slwd(i) = srfrad(i) - frsa(i)
            fld = slwd (i)
            write(6,*) ' longwave down at surface  = ',fld,' wm-2' 
            clt = clrlt(i)
            write(6,*) ' longwave clear outgoing   = ',clt,' wm-2' 
            cls = clrls(i)
            write(6,*) ' longwave clear net srf    = ',cls,' wm-2' 
            cfl = clt - flw
            write(6,*) ' longwave cloud forcing    = ',cfl,' wm-2' 
c     
            write(6,*) ' ------------------ '
            cfn = scf + cfl
            write(6,*) ' net cloud forcing         = ',cfn,' wm-2' 
            write(6,*) ' ------------------------------ '
            write(6,*) ' --- heating rates in K/day --- '
            write(6,*) '   level   pr(mb)     qrs           qrl '
c     
            do 200 k=1,plev
               qrsday = qrs(i,k) * 86400.
               qrlday = qrl(i,k) * 86400.
               pmb    = pmid(i,k) / 100.
               write(6,199) k,pmb,qrsday,qrlday
 199           format(2x,i4,2x,f8.3,2x,2(f12.7,2x))
 200        continue
c     
 100     continue
c     
         write(6,*) ' --- done with Column Model Rad computation ---'
      endif !endif debug
cz    here's where the short and longwave heating rates are transferred
cz    to the cloud model arrays.  The units are degrees Kelvin per second.
cz    altitude loop runs down from top of cloud (so indices increase in F77)

      scf = sabtp(1) - clrst(1)
      cfl = clrlt(1) - firtp(1)
      lwcf=cfl
      swcf=scf
cz
cz find the surface cloud forcing
cz
      frs= frsa(1)
      cls= clrss(1)
      swscf = frs-cls
cz
cz find the albedo
cz
      sol = solin(1)
      sab = sabtp(1)
      if( sol .le. 0. ) then
         alb = -999.
      else
         alb = (sol-sab)/sol   
      endif
      albedo = alb
cz
      do 695 i=1,plon
         do 696 k=plev+1-bccmil-numlay,plev+1-bccmil-1
            kcz=plev-k-bccmil+1
            qrscld(kcz) = qrs(i,k)
            qrlcld(kcz) = qrl(i,k)
            qrcld(kcz) = qrs(i,k)+qrl(i,k)
 696     continue
         do 697 k=plev+1-bccmil-numlay,plev+1-bccmil
            kcz=plev-k-bccmil+1
            envpint(kcz) = pint(i,k)
 697  continue
 695  continue
c     
cz      stop
      return
      end

cdir$ nolist
      subroutine stparm
c     
c     set certain model parameters for use in the column radiation model
c
c anncyc must be set for use in ozone computation
c nstep must be 0, as well as irad=1 and iradae=1 so that longwave 
c absorptivity and emissivity computation will be done
c
c-----------------------------------------------------------------------
c
c     !   Model grid point parameters. Include basic parameter deck.
c
c
c     !   Basic grid point resolution parameters
c
      integer plev,plon,plat,plonp2,plevp,pcnst,plond,platd,plevd
      integer nxpt,jintmx,i1,j1,j1m
c
      parameter(plev=  133,        ! number of vertical levels
     $          plevp=plev+1,     ! number of vertical interfaces
     $          plon=1,           ! number of longitudes (T42)
     $          plat=1,           ! number of latitudes (T42)
     $          plonp2=plon+2,    ! length necessary for fft
     $          pcnst=1,          ! number of constituents
     $          nxpt=1,           !
     $          jintmx=2,         !
     $          plond=plon + 1 + 2*nxpt, ! slt extended domain longitude
     $          platd=plat + 2*nxpt + 2*jintmx, ! slt extended domain lat.
     $          plevd=plev*(3+pcnst), ! fold plev,pcnst indices into one
     $          i1=1+nxpt,        ! model starting longitude (3-d)
     $          j1=1+nxpt+jintmx, ! model starting latitude (3-d)
     $          j1m=j1-1)         ! model starting offset (3-d)
c
c
c     !   Auxiliary grid point resolution parameters
c
      integer
     $    plnlv,pln2lv,plndlv,plat2,
     $    plevmx,pscal,ptflen,pdflen,
     $    pbflnb,pbflna,pflenb,pflena,plnbuf,pflds,pftyp,prbd,
     $    plenhi,plenhc,plenbr,plenhr,
     $    pmulti,psingl,ptplen,ptapes,phtbf
c
      parameter(
     $    plnlv=plon*plev,      !  length of multilevel ht field
     $    pln2lv=plonp2*plev,   !  length of multilevel buffer field
     $    plndlv=plond*plev,       ! length of multilevel 3-d field
     $    plat2=plat/2)         !
C
c
      parameter(
     $    plevmx=4,                        ! number of subsurface levels
     $    pscal=4,                         ! number of passive scalars
     $    ptflen=((3+2*pcnst)*plev + (9+pcnst))*plond,! len time dep fld
     $    pdflen=(4*plev + 7)*plond)       ! length of diagnostic fields
      parameter(
     $    pmulti=22 + pcnst*5,            ! number of multilevel fields
c                                         ! primary history tape
     $    psingl=35+plevmx,               ! number of single level field
     $    ptplen=pmulti+psingl)           ! total fields to primary tape
      parameter(
     $    pbflnb=2*plndlv+ptflen+pdflen,   ! length of nprgtl buffer
     $    pbflna=(4+6*plev)*plond,        ! length of nprg
     $    phtbf=(pmulti*plev+psingl)*plon, ! length of history tape buf
     $    pflenb=((pbflnb+pbflna+phtbf)/512+1)*512,  !  padded length of
c                                                   !   + nprgm1 + hbuf
     $    pflena=(pbflna/512+1)*512,       ! padded length of nprg
     $    plnbuf=pflenb+pflena,           ! total buffer length
     $    pflds=14,            ! number of fields on initial dataset
     $    ptapes=1)            ! maximum number of history tapes allowed
c
c     !       history tape header parameters
c
      parameter (
     $   prbd=3,                           ! records before data on ht
     $   pftyp=42,                         ! units digit of ht format ty
     $   plenhi=37+3*ptplen,               ! length of integer header re
     $   plenhc=76+2*ptplen,               ! length of character header
     $   plenbr=3*(2*plev+1)+2*plat,       ! length of real header buffe
     $   plenhr=plenbr)                    ! length of real header recor
c-----------------------------------------------------------------------
C
C FLAGS ASSOCIATED WITH ANNUAL CYCLE
C
C IF ANNCYC=.TRUE., SEA SURFACE TEMPERATURES
C AND OZONE PATH LENGTHS GO THROUGH AN ANNUAL
C CYCLE DETERMINED BY SPLINE COEFFICIENTS FOR
C EACH GRID POINT
C
C IF ANNCYC=.FALSE.,SST'S AND O3 PATH LENGTHS
C CONSTANT IN TIME
C
C IF HYDRO =.TRUE., HYDROLOGIC CYCLE IN MODEL
C WILL BE ENABLED (REQUIRES ANNCYC TO BE .TRUE.)
C
C ITSST IS SEA SURFACE TEMPERATURE UPDATE FREQUENCY (ITERS)
C
      common/comanc/anncyc,hydro,itsst
      logical       anncyc,hydro
      integer itsst
c-----------------------------------------------------------------------
C
C model control time variables
C
      integer
     $      nrstrt  ,nstep   ,nstepr  ,nestep  ,nstop   ,
     $      mdbase  ,msbase  ,mdcur   ,mscur   ,
     $      mbdate  ,mbsec   ,mcdate  ,mcsec   ,
     $      nndbas  ,nnsbas  ,nnbdat  ,nnbsec
      real  calday
c
      common /comtim/
     $      calday  ,nrstrt  ,nstep   ,nstepr  ,nestep  ,nstop   ,
     $      mdbase  ,msbase  ,mdcur   ,mscur   ,
     $      mbdate  ,mbsec   ,mcdate  ,mcsec   ,
     $      nndbas  ,nnsbas  ,nnbdat  ,nnbsec
C
c-----------------------------------------------------------------------
c
c radiation control variables
c
c fradsw = .t. iff full shortwave computation
c fradlw = .t. iff full longwave computation
c
c irad = iteration frequency for radiation computation
c iradae = iteration frequency for absorptivity/
c emissivity computation
c
c
      integer iradae,irad,naclw,nacsw,fnlw,fnsw
      logical aeres
c
      common/crdctl/iradae, irad, naclw, nacsw, fnlw, fnsw, aeres
c
c-----------------------------------------------------------------------
C
CL            HISTORY TAPE HEADER, RECORD NUMBER 1 (INTEGER VALUES)
CL            PROVIDES DATA FILE DESCRIPTION FOR HISTORY TAPE OUTPUT
C
      COMMON /COMHDI/
     +   LENHDI(1), MFTYP,    MFILH,    MFILTH,     NRBD
     +  ,MAXSIZ,    NDAVU,    MSPHER,   MLON,       NLONW
     +  ,MOREC,     MLEV,     MTRM,     MTRN,       MTRK
     +  ,NFLDH,     NSTEPH,   NSTPRH,   NITSLF
     +  ,NDBASE,   NSBASE,    NDCUR,    NSCUR
     +  ,NBDATE,   NBSEC,     NCDATE,   NCSEC
     +  ,MDT,      MHISF,     MFSTRT,   LENHDC,     LENHDR
     +  ,MPSIG,    MPLAT,     MPWTS,    MPFLDS,     MPCFLD
C             ARRAY OF INTEGER FIELD LIST INFORMATION
     $  ,mflds(3,ptplen)
C
CL            HISTORY TAPE HEADER, RECORD NUMBER 2 (CHARACTER VALUES)
C
      COMMON /COMHDC/
     +   MCASE,    MCSTIT
     +  ,LNHSTC,   LDHSTC,   LTHSTC,   LSHSTC
     +  ,LNHSTP,   LDHSTP,   LTHSTP,   LSHSTP
     +  ,LNHSTF,   LDHSTF,   LTHSTF,   LSHSTF
     +  ,LNHSTI,   LDHSTI,   LTHSTI,   LSHSTI
     +  ,LNHSTA,   LDHSTA,   LTHSTA,   LSHSTA
C             ARRAY OF CHARACTER FIELD LIST INFORMATION
     $  ,mcflds(2,ptplen)
      CHARACTER*8 MCASE
     + ,LDHSTC,   LTHSTC,   LSHSTC
     + ,LDHSTP,   LTHSTP,   LSHSTP
     + ,LDHSTF,   LTHSTF,   LSHSTF
     + ,LDHSTI,   LTHSTI,   LSHSTI
     + ,LDHSTA,   LTHSTA,   LSHSTA
     + ,MCFLDS
      CHARACTER*80 MCSTIT
     + ,LNHSTC,   LNHSTP,   LNHSTF,   LNHSTI,   LNHSTA
C
CL            HISTORY TAPE HEADER, RECORD NUMBER 3 (REAL VALUES)
C
      COMMON /COMHDR/
C             BUFFER BUFHD CONTAINS SIGMA,LATITUDE AND GAUSSIAN WEIGHTS
     +   REALHD(1),BUFHD(PLENBR)
C
c     !  Store certain header values for all history tapes
c
      common /comhdt/
     $       nstepht(ptapes)
      integer nstepht
C
      common /comhdtc/
     $       ldhstct(ptapes), lthstct(ptapes), lshstct(ptapes)
C
      character*8 ldhstct, lthstct, lshstct
c-----------------------------------------------------------------------
c
c... from 'comanc'
c
      hydro  = .false.
      anncyc = .true.
c
c... from 'comtim'
c
      nstep  =  0
c
c... from 'crdctl'
c
      irad   = 1
      iradae = 1
c
      return
      end
c-----------------------------------------------------------------------
      subroutine cmpsol(clat,coszrs,eccf,loctim)
c
c compute solar zenith angle and other diagnositics from 
c the time of year and latitude:
c
c-----------------------------------------------------------------------
c
c     !   Basic grid point resolution parameters
c
      integer plev,plon,plat,plonp2,plevp,pcnst,plond,platd,plevd
      integer nxpt,jintmx,i1,j1,j1m
c
      parameter(plev=  133,        ! number of vertical levels
     $          plevp=plev+1,     ! number of vertical interfaces
     $          plon=1,           ! number of longitudes (T42)
     $          plat=1,           ! number of latitudes (T42)
     $          plonp2=plon+2,    ! length necessary for fft
     $          pcnst=1,          ! number of constituents
     $          nxpt=1,           !
     $          jintmx=2,         !
     $          plond=plon + 1 + 2*nxpt, ! slt extended domain longitude
     $          platd=plat + 2*nxpt + 2*jintmx, ! slt extended domain lat.
     $          plevd=plev*(3+pcnst), ! fold plev,pcnst indices into one
     $          i1=1+nxpt,        ! model starting longitude (3-d)
     $          j1=1+nxpt+jintmx, ! model starting latitude (3-d)
     $          j1m=j1-1)         ! model starting offset (3-d)
c
c-----------------------------------------------------------------------
c
c radiation constants
c
      real gravit,   ! gravitational acceleration
     $        rga,   ! 1 over gravit
     $      cpair,   ! heat capacity air at constant pressure
     $     epsilo,   ! ratio mmw h2o to mmw air
     $       sslp,   ! standard pressure
     $     stebol,   ! stephan boltzmann constant
     $     rgsslp,   ! 0.5 / (gravit*sslp)
     $     co2vmr,   ! co2 volume mixing ratio
     $      dpfo3,   ! Doppler factor for o3
     $     dpfco2,   ! Doppler factor for co2
     $     dayspy,   ! solar days in one year
     $        pie    ! pie
c
      common/crdcon/gravit,    rga,  cpair,  epsilo,   sslp,
     $              stebol, rgsslp, co2vmr,   dpfo3, dpfco2,
     $              dayspy,    pie
c
c-----------------------------------------------------------------------
C
C model control time variables
C
      integer
     $      nrstrt  ,nstep   ,nstepr  ,nestep  ,nstop   ,
     $      mdbase  ,msbase  ,mdcur   ,mscur   ,
     $      mbdate  ,mbsec   ,mcdate  ,mcsec   ,
     $      nndbas  ,nnsbas  ,nnbdat  ,nnbsec
      real  calday
c
      common /comtim/
     $      calday  ,nrstrt  ,nstep   ,nstepr  ,nestep  ,nstop   ,
     $      mdbase  ,msbase  ,mdcur   ,mscur   ,
     $      mbdate  ,mbsec   ,mcdate  ,mcsec   ,
     $      nndbas  ,nnsbas  ,nnbdat  ,nnbsec
C
c-----------------------------------------------------------------------
c
c input/output arguments
c
       real coszrs(plon),  ! cosine solar zenith angle
     $      loctim(plon)   ! local time in hours

c compute eccentricity factor (sun-earth distance factor)
c
      theta = 2.*pie*calday/dayspy
      eccf  = 1.000110 + .034221*cos(theta) + .001280*sin(theta) +
     $        .000719*cos(2.*theta) + .000077*sin(2.*theta)
c
c solar declination in radians:
c
      delta = .006918 - .399912*cos(theta) + .070257*sin(theta) -
     $        .006758*cos(2.*theta) + .000907*sin(2.*theta) -
     $        .002697*cos(3.*theta) + .001480*sin(3.*theta)
c
c compute local cosine solar zenith angle,
c
      sinc = sin(clat)
      sind = sin(delta)
      cosc = cos(clat)
      cosd = cos(delta)
c
c calday is the calender day for greenwich, including fraction
c of day; the fraction of the day represents a local time at
c greenwich; to adjust this to produce a true instantaneous time
c for other longitudes, we must correct for the local time change:
c
      do 10 i=1,plon
c
         phi       = calday + (real(i-1)/real(plon))
         cphase    = cos(2.*pie*phi)
         coszrs(i) = sinc*sind - cosc*cosd*cphase
         loctim(i) = 12.*acos(cphase)/pie
c
   10 continue
c
      return
      end
c-----------------------------------------------------------------------
      subroutine radoz2(jlat,pmidm1,pintm1,plol,plos)
      return
      end
c-----------------------------------------------------------------------
      subroutine radrda(nunit,albvss,plon,plat,conv)
      return
      end
c-----------------------------------------------------------------------
      integer function norsou(nrow)
      norsou = 0
      return
      end
c-----------------------------------------------------------------------
      subroutine getdat(gravx  ,cpairx ,epsilox,stebolx,   nrow,
     $                  clat   ,oro    ,sndpth ,ts     ,     tg,
     $                  ps     ,pmid   ,pint   ,pmln   ,   piln,
     $                  t      ,h2ommr ,plol   ,plos   , cldfrc,
     $                  clwp   , effcld,
     &     dayyr,rlat,o3mmr,
     &     lwmac,crecz,
     &     bccmil,tccmil,numlay,ncclev,
     &     cbi,cti,ncldlay,
     &     idebug)
c-----------------------------------------------------------------------
c
c interface routine for column model that both initializes
c certain constants and reads external data:
c 
c o3 mass mixing ratios are read in, but the model also requires the
c path lengths; they are computed here
c
c also, from the cloud input (fraction and liquid water path), the
c cloud longwave emissivity must be computed; this is done here
c
c This routine is heavily modified to now accept input from the 
c cirrus cloud model instead of a file.  The cirrus cloud model
c actually reads the file itself then inserts the cloud and passes
c the relevant arrays to getdat.
c
c-----------------------------------------------------------------------
      implicit none
c-----------------------------------------------------------------------
c
c     !   Basic grid point resolution parameters
c
      integer plev,plon,plat,plonp2,plevp,pcnst,plond,platd,plevd
      integer nxpt,jintmx,i1,j1,j1m
c
      parameter(plev=  133,        ! number of vertical levels
     $          plevp=plev+1,     ! number of vertical interfaces
     $          plon=1,           ! number of longitudes (T42)
     $          plat=1,           ! number of latitudes (T42)
     $          plonp2=plon+2,    ! length necessary for fft
     $          pcnst=1,          ! number of constituents
     $          nxpt=1,           !
     $          jintmx=2,         !
     $          plond=plon + 1 + 2*nxpt, ! slt extended domain longitude
     $          platd=plat + 2*nxpt + 2*jintmx, ! slt extended domain lat.
     $          plevd=plev*(3+pcnst), ! fold plev,pcnst indices into one
     $          i1=1+nxpt,        ! model starting longitude (3-d)
     $          j1=1+nxpt+jintmx, ! model starting latitude (3-d)
     $          j1m=j1-1)         ! model starting offset (3-d)
c
c-----------------------------------------------------------------------
C
C model control time variables
C
      integer
     $      nrstrt  ,nstep   ,nstepr  ,nestep  ,nstop   ,
     $      mdbase  ,msbase  ,mdcur   ,mscur   ,
     $      mbdate  ,mbsec   ,mcdate  ,mcsec   ,
     $      nndbas  ,nnsbas  ,nnbdat  ,nnbsec
      real  calday
c
      common /comtim/
     $      calday  ,nrstrt  ,nstep   ,nstepr  ,nestep  ,nstop   ,
     $      mdbase  ,msbase  ,mdcur   ,mscur   ,
     $      mbdate  ,mbsec   ,mcdate  ,mcsec   ,
     $      nndbas  ,nnsbas  ,nnbdat  ,nnbsec
C
c-----------------------------------------------------------------------
c
c surface albedo data
c
c vs = 0.2 - 0.7 micro-meters wavelength range
c ni = 0.7 - 5.0 micro-meters wavelength range
c
c s  = strong zenith angle dependent surfaces
c w  = weak   zenith angle dependent surfaces
c
c the albedos are computed for a model grid box by ascribing values to
c high resolution points from a vegetation dataset, then linearlly 
c averaging to the grid box value; ocean and land values are averaged
c together along coastlines; the fraction of every grid box that has
c strong zenith angle dependence is included also.
c
      real albvss, ! grid box alb for vs over strng zn srfs
     $     albvsw, ! grid box alb for vs over weak  zn srfs
     $     albnis, ! grid box alb for ni over strng zn srfs
     $     albniw, ! grid box alb for ni over weak  zn srfs
     $     frctst  ! fraction of area in grid box strng zn 
c
      common /crdalb/albvss(plond,plat),albvsw(plond,plat),
     $               albnis(plond,plat),albniw(plond,plat),
     $               frctst(plond,plat)
c
c
c surface boundary data
c
c vegtyp is used to specify the thermal properites of the surface,
c as well as determine the location of permanent land ice points;
c it is the dominant surface type within the model grid box based
c on a high resolution vegetation type dataset;
c it encoded in the following manner:
c
c   1        ocean
c   2        sea ice
c   3        permanent land ice
c   4        tropical evergreen forest
c   5        deciduous forest
c   6        grassland/tundra
c   7        desert
c
c rghnss is the aerodynamic roughness length for the grid box, computed
c by linear averaging of the values ascribed to high resolution 
c vegetation dataset values; ocean and land values are averaged together
c at coastlines.
c
c evapf is the ratio of actual to potential evaporation, and is computed
c from the high resolution vegetation dataset in a manner similar to the
c aerodynamic roughness.
c
c vevapf allows for variable snow cover, where the underlying
c evaporability factor is modified; see radalb.
c
c snwjan and snwjly are mean climatological snow depths (liquid water
c equivalent) used to compute the actual daily values of snow cover.
c
      real vegtyp, ! surface thermal type, based on veg type
     $     rghnss, ! aerodynamic roughness length
     $     evapf , ! constant surface evaporability
     $     vevapf, ! variable surface evaporability
     $     snwjan, ! snow cover (liq water equiv) for january
     $     snwjly  ! snow cover (liq water equiv) for july
c
      common /crdsrf/vegtyp(plond,plat),rghnss(plond,plat),
     $               evapf (plond,plat),vevapf(plond,plat),
     $               snwjan(plond,plat),snwjly(plond,plat)
c
c-----------------------------------------------------------------------
c
c radiation constants
c
      real gravit,   ! gravitational acceleration
     $        rga,   ! 1 over gravit
     $      cpair,   ! heat capacity air at constant pressure
     $     epsilo,   ! ratio mmw h2o to mmw air
     $       sslp,   ! standard pressure
     $     stebol,   ! stephan boltzmann constant
     $     rgsslp,   ! 0.5 / (gravit*sslp)
     $     co2vmr,   ! co2 volume mixing ratio
     $      dpfo3,   ! Doppler factor for o3
     $     dpfco2,   ! Doppler factor for co2
     $     dayspy,   ! solar days in one year
     $        pie    ! pie
c
      common/crdcon/gravit,    rga,  cpair,  epsilo,   sslp,
     $              stebol, rgsslp, co2vmr,   dpfo3, dpfco2,
     $              dayspy,    pie
c
c-----------------------------------------------------------------------
c
c output arguments
c
      real gravx,       ! gravitational acceleration (m/s**2)
     $     cpairx,      ! heat capacity dry air at constant prs (J/kg/K)
     $     epsilox,     ! ratio mean mol weight h2o to dry air
     $     stebolx      ! Sefan-Boltzmann constant (W/m**2/K**4)
c
      integer nrow      ! model latitude index
c
      real clat,               ! model latitude in radians
     $     oro(plond),         ! land surface flag
     $     sndpth(plond),      ! snow depth (liquid water equivalent)
     $     ts(plond),          ! surface (air)  temperature
     $     tg(plond),          ! surface (skin) temperature
     $     ps(plond),          ! model surface pressure field
     $     pmid(plond,plev),   ! pressure at model mid-levels 
     $     pint(plond,plevp),  ! pressure at model interfaces 
     $     pmln(plond,plev),   ! ln(pmid)
     $     piln(plond,plevp),  ! ln(pint)
     $     t(plond,plev),      ! atmospheric temperature
     $     h2ommr(plond,plev), ! moisture field
     $     plol(plond,plevp),  ! o3 pressure weighted path length
     $     plos(plond,plevp),  ! o3 path length
     $     cldfrc(plond,plevp),! cloud fraction
     $     clwp(plond,plev),   ! cloud liquid water path (g/m**2)
     $     effcld(plond,plevp) ! effective cloud fraction
c
c local workspace
c
      real dayyr(plond),       ! day of year
     $     rlat(plond),        ! latitude input
     $     o3mmr(plond,plev),  ! o3 mass mixing ratio
     $     emis(plond,plev),   ! cloud emissivity for longwave
     $     ptop,               ! top layer interface pressure
     $     pbot                ! bottom layer interface pressure
c
cz      integer lev(plev),       ! level input
      integer 
     $                i,       ! longitude index
     $                k        ! level  index
c
cz      character*80 label
c
      real     v0,  ! volume of a gas at stp (cm**3/mol)
     $         p0,  ! standard pressure (dynes/cm**2)
     $        amd,  ! effective molecular weight of dry air (g/mol)
     $        amo,  ! molecular weight of ozone (g/mol)
     $        cpl,  ! constant in ozone path length to mixing ratio
     $      cpwpl,  ! pressure weighted ozone path length constant
     $       vmmr   ! ozone volume mixing ratio
c
      data v0    /  22413.6   /
      data p0    /  1.01325e6 /
      data amd   /  28.9644   /
      data amo   /  48.0000   /
c
c------------------------------------------------------------------------
c
c     zender cirrus cloud hack workspace
c
      real
     $     lwmac(0:numlay+1), ! longwave mass absorption coeff (m^2/kg)
     &     crecz(0:numlay+1)  !effective radius of particles (meters)
c
      integer
     &     idebug,    !debugging level from cloud model
     &     bccmil,    !bottom ccm2 interface level
     &     tccmil,    !top ccm2 interface level
     &     numlay,    !num_layer in cloud model
     &     ncclev,    !num_ccm2_level (with overlaps already subtracted)
     &     cbi,       !cloud bottom index
     &     cti,       !cloud top index
     &     ncldlay    !number of layers in cloud itself
c
c     end of zender workspace alterations
c------------------------------------------------------------------------
c
c set fundamental constants (mks):
c
      gravx   =   9.80616
      cpairx  =   1.00464e3
      epsilox =   0.622
      stebolx =   5.67e-8
c
      nrow    =   1
c
c begin read of data:
c
      do 100 i=1,1
c
c        read(5,101)  label
c  101   format(a80)
c        write(6,*)   label
cc
c        read(5,101)  label
c        write(6,*)   label
cc
c        read(5,101)  label
c        write(6,*)   label
cc
c        read(5,*)    dayyr(i)
         if( idebug.eq.45.or.idebug.eq.53 ) then
            write(6,*) ' day of year (1..365)  = ',dayyr(i)
         endif
cc
c        read(5,*)    rlat(i)
         if( idebug.eq.45.or.idebug.eq.53 ) then
            write(6,*) ' latitude (-90 to +90) = ',rlat(i)
         endif
cc
c        read(5,101)  label
c        write(6,*)   label
cc
         if( idebug.eq.45.or.idebug.eq.53 ) then
            print *,'level   p(mb)    t(k)   h2ommr (g/g) o3mmr ',
     +           '(g/g) cld cvr  cld lwp(g/m2)'
         endif
         do 200 k=1,plev
c          read(5,*) lev(k),pmid(i,k),t(i,k),h2ommr(i,k),o3mmr(i,k)
c     +              ,cldfrc(i,k),clwp(i,k)
            if( idebug.eq.45.or.idebug.eq.53 ) then
               write(6,99) k   ,pmid(i,k),t(i,k),h2ommr(i,k),o3mmr(i,k)
     +              ,cldfrc(i,k),clwp(i,k)
            endif
 200     continue
 99      format(2x,i3,1x,6(1pe11.4,1x))
c
c..... currently, cld requires extra 'below surface' value:
c
        cldfrc(i,plevp) = 0.0
c
c        read(5,*)     ps(i)
        if( idebug.eq.45.or.idebug.eq.53 ) then
           write(6,*) ' srf pressure = ',ps(i)
        endif
c
c.......... convert pressures from mb to pascals and define 
c.......... interface pressures:
c
        ps(i) = ps(i) * 100.
        do 125 k=1,plev
c
          pmid(i,k) = pmid(i,k) * 100.
          pmln(i,k) = alog(pmid(i,k))
c
  125   continue
        do 150 k=1,plevp
c
          if( k .eq. 1 ) then
            pint(i,k) = pmid(i,k) / 2.0
          else if ( k .gt. 1 .and. k .le. plev ) then
            pint(i,k) = 0.5 * (pmid(i,k-1) + pmid(i,k))
          else if ( k .eq. plevp ) then
            pint(i,k) = ps(i)
          endif
          piln(i,k) = alog(pint(i,k))
c
  150   continue
c
c        read(5,*)       ts(i)
        if( idebug.eq.45.or.idebug.eq.53 ) then
           write(6,*) ' surface air temperature  = ',ts(i)
        endif
c
c        read(5,*)       tg(i)
        if( idebug.eq.45.or.idebug.eq.53 ) then
           write(6,*) ' surface skin temperature = ',tg(i)
        endif
c
c        read(5,*)       oro(i)
        if( idebug.eq.45.or.idebug.eq.53 ) then
           write(6,*) ' surface type (oro) flag  = ',oro(i)
        endif
c
c        read(5,*)       rghnss(i,1)
        if( idebug.eq.45.or.idebug.eq.53 ) then
           write(6,*) ' surface roughness        = ',rghnss(i,1)
        endif
c
c        read(5,*)       sndpth(i)
        if( idebug.eq.45.or.idebug.eq.53 ) then
           write(6,*) ' snow depth (liq wat)     = ',sndpth(i)
        endif
c
c        read(5,*)       albvss(i,1)
        if( idebug.eq.45.or.idebug.eq.53 ) then
           write(6,*) ' albvss     = ',albvss(i,1)
        endif
c
c        read(5,*)       albnis(i,1)
        if( idebug.eq.45.or.idebug.eq.53 ) then
           write(6,*) ' albnis     = ',albnis(i,1)
        endif
c
c        read(5,*)       albvsw(i,1)
        if( idebug.eq.45.or.idebug.eq.53 ) then
           write(6,*) ' albvsw     = ',albvsw(i,1)
        endif
c
c        read(5,*)       albniw(i,1)
        if( idebug.eq.45.or.idebug.eq.53 ) then
           write(6,*) ' albniw     = ',albniw(i,1)
        endif
c
c        read(5,*)       frctst(i,1)
        if( idebug.eq.45.or.idebug.eq.53 ) then
           write(6,*) ' frctst     = ',frctst(i,1)
        endif
c
  100 continue
c
      calday = dayyr(1)
      pie    = 4.*atan(1.)
      clat   = rlat(1)*(pie/180.)
c
c compute ozone path lengths from mixing ratio:
c
c constants for following sums:
c
      gravit = gravx * 100.
      cpl    = v0  / (amd * gravit)
      cpwpl  = 0.5 * v0 / (amd * gravit * p0)
      vmmr   = amd / amo
c
      do 225 i=1,plon
c
c set top level to space path lengths:
c
       pbot = pint(i,1) * 10.
c
       plos(i,1) = cpl * vmmr * o3mmr(i,1) * pbot
c
       plol(i,1) = cpwpl * vmmr * o3mmr(i,1) *
     +                              (pbot*pbot)
c
       ptop = 0.0
c
c set rest of level path lengths:
c
       do 250 k=2,plevp
c
          ptop = pint(i,k-1) * 10.
          pbot = pint(i,k)   * 10.
c
          plos(i,k) = plos(i,k-1) +
     +             (cpl * vmmr * o3mmr(i,k-1) * (pbot - ptop))
c
          plol(i,k) = plol(i,k-1) +
     +                    (cpwpl * vmmr * o3mmr(i,k-1) *
     +                    (pbot*pbot - ptop*ptop))
c
  250 continue
  225 continue
c
c compute effective cloud cover
c
      call cldems(clwp,emis,
     &     lwmac,crecz,
     &     bccmil,tccmil,numlay,ncclev,
     &     cbi,cti,ncldlay,
     &     idebug)
c
      do 300 k=1,plev
         do 400 i=1,plon
            effcld(i,k) = cldfrc(i,k)*emis(i,k)
  400    continue
  300 continue
c
c cloud cover at surface interface always zero
c
      do 500 i=1,plon
         effcld(i,plevp) = 0.
         cldfrc(i,plevp) = 0.
  500 continue
c
      if( idebug.eq.45.or.idebug.eq.53 ) then
         write(6,*) ' -----end of profile data input----- '
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine cldems(clwp,emis,
     &     lwmac,crecz,
     &     bccmil,tccmil,numlay,ncclev,
     &     cbi,cti,ncldlay,
     &     idebug)
c-----------------------------------------------------------------------
c
c compute cloud emissivity using cloud liquid water path (g/m**2)
c
c-----------------------------------------------------------------------
      implicit none
c-----------------------------------------------------------------------
c
c     !   Basic grid point resolution parameters
c
      integer plev,plon,plat,plonp2,plevp,pcnst,plond,platd,plevd
      integer nxpt,jintmx,i1,j1,j1m
c
      parameter(plev=  133,        ! number of vertical levels
     $          plevp=plev+1,     ! number of vertical interfaces
     $          plon=1,           ! number of longitudes (T42)
     $          plat=1,           ! number of latitudes (T42)
     $          plonp2=plon+2,    ! length necessary for fft
     $          pcnst=1,          ! number of constituents
     $          nxpt=1,           !
     $          jintmx=2,         !
     $          plond=plon + 1 + 2*nxpt, ! slt extended domain longitude
     $          platd=plat + 2*nxpt + 2*jintmx, ! slt extended domain lat.
     $          plevd=plev*(3+pcnst), ! fold plev,pcnst indices into one
     $          i1=1+nxpt,        ! model starting longitude (3-d)
     $          j1=1+nxpt+jintmx, ! model starting latitude (3-d)
     $          j1m=j1-1)         ! model starting offset (3-d)
c
c-----------------------------------------------------------------------
c
c input argument
c
      real clwp(plond,plev)    ! cloud liquid water path (g/m**2)
c
c output argument
c
      real emis(plond,plev)    ! cloud emissivity (fraction)
c
c local workspace
c
      integer          i,  ! longitude index
     $                 k   ! level index
c
      real kabs         ! longwave absorption coefficiant (m**2/g)
c------------------------------------------------------------------------
c
c     zender cirrus cloud hack workspace
c
      real
     $     lwmac(0:numlay+1),  ! longwave mass absorption coeff (m^2/kg)
     $     diffusivity         ! coeff of lwmac = 1.66 
c
      real
     $     cre(plev),          ! cloud effective radius (micro-meters)
     &     crecz(0:numlay+1),  ! effective radius of particles (meters)
     &     cldefr,      ! Universal cloud effective radius in micro-meters
     &     alphae,      ! LW band weighted fitting param. (E&C table 3)
     &     gammae,      ! LW band weighted fitting param. (E&C table 3)
     &     eclwmac      ! Ebert/Curry LW mass abs. coeff. (E&C eqn. 12)
c
      integer
     &     kcz,       !index for index trickery
     &     idebug,    !debugging level from cloud model
     &     bccmil,    !bottom ccm2 interface level
     &     tccmil,    !top ccm2 interface level
     &     numlay,    !num_layer in cloud model
     &     ncclev,    !num_ccm2_level (with overlaps already subtracted)
     &     cbi,       !cloud bottom index
     &     cti,       !cloud top index
     &     ncldlay    !number of layers in cloud itself
c
c     end of zender workspace alterations
c------------------------------------------------------------------------
c
      data kabs  / 0.1 /
      data diffusivity  / 1.66 /
cz
      data cldefr / 10.0 /       ! microns
      data alphae / .29132e-03 / ! meter^2 gram^-1
      data gammae / 1.059344 /   ! micron^-1 
cz
      do 40 k=1,plev
         do 30 i=1,plon
cz    
cz To implement the longwave Ebert and Curry parameterization is 
cz straightforward:  E&C give fitting parameters alpha and gamma
cz weighted for five longwave bands from 4 microns on up.  Since
cz CCM2 uses only one LW band, just find the weighted average,
cz i.e., Planckian mean, alpha and gamma.
cz
cz debug 47 is use the Ebert Curry scheme
cz debug 53 is use the Slingo scheme
cz
cz The effective radius is needed everywhere for the Ebert&Curry scheme
cz but only outside the cloud grid for the full microphysical scheme.
cz
            if(idebug.eq.53) then
cz     all levels of Slingo scheme  
               emis(i,k) = 1. - exp(-kabs*clwp(i,k))
            else if((k.ge.(plev+1-bccmil-numlay)).and.
     &              (k.le.(plev+1-bccmil-1))) then
cz    it is inside cloud grid and not slingo style so find the cz layer
               kcz=plev-k-bccmil+1
               if(idebug.eq.47) then
cz    employing Ebert & Curry scheme internal to cloud grid
c     remember to convert the meters to microns
                  cre(k) = crecz(kcz)*1.e6
                  eclwmac=alphae+gammae/cre(k)
                  emis(i,k) = 1. - exp(-eclwmac*diffusivity*clwp(i,k))
               else
cz     cloud grid levels of complete cirrus cloud Mie/microphysical scheme
cz     convert the lwmac from m^2/kg to m^2/g before exponentiating
cz     since clwp is already in g/m^2
cz3456789012345678901234567890123456789012345678901234567890123456789012
                  emis(i,k)=
     $                 1.-exp(-lwmac(kcz)*diffusivity*clwp(i,k)/1000.)
               endif
            else
cz     levels outside cloud grid in both full mie/microphysical scheme
cz     and Ebert-Curry scheme. set the effective radius to the universal
cz    effective radius
               cre(k) = cldefr
               eclwmac=alphae+gammae/cre(k)
               emis(i,k) = 1. - exp(-eclwmac*diffusivity*clwp(i,k))
            endif
cz
czcz    the newest twist is to remove the cloud completely from the
czcz    longwave computations in order to obtain clear sky LW fluxes
czcz    at the cloud boundaries. in that case must zero all emissivities.
cz            emis(i,k) = 0.
 30      continue
 40   continue
c     
c     done
c
      return
      end
      subroutine radini(gravx,cpairx,epsilox,stebolx)
c---------------------------------------------------------------------
c
c modified version of initialization for radiation scheme; done
c for the column radiation model
c
c note: radiation scheme uses cgs units
c
c---------------------------------------------------------------------
c.......      implicit none
c---------------------------------------------------------------------
c
c     !   Model grid point parameters. Include basic parameter deck.
c
c
c     !   Basic grid point resolution parameters
c
      integer plev,plon,plat,plonp2,plevp,pcnst,plond,platd,plevd
      integer nxpt,jintmx,i1,j1,j1m
c
      parameter(plev=  133,        ! number of vertical levels
     $          plevp=plev+1,     ! number of vertical interfaces
     $          plon=1,           ! number of longitudes (T42)
     $          plat=1,           ! number of latitudes (T42)
     $          plonp2=plon+2,    ! length necessary for fft
     $          pcnst=1,          ! number of constituents
     $          nxpt=1,           !
     $          jintmx=2,         !
     $          plond=plon + 1 + 2*nxpt, ! slt extended domain longitude
     $          platd=plat + 2*nxpt + 2*jintmx, ! slt extended domain lat.
     $          plevd=plev*(3+pcnst), ! fold plev,pcnst indices into one
     $          i1=1+nxpt,        ! model starting longitude (3-d)
     $          j1=1+nxpt+jintmx, ! model starting latitude (3-d)
     $          j1m=j1-1)         ! model starting offset (3-d)
c
c
c     !   Auxiliary grid point resolution parameters
c
      integer
     $    plnlv,pln2lv,plndlv,plat2,
     $    plevmx,pscal,ptflen,pdflen,
     $    pbflnb,pbflna,pflenb,pflena,plnbuf,pflds,pftyp,prbd,
     $    plenhi,plenhc,plenbr,plenhr,
     $    pmulti,psingl,ptplen,ptapes,phtbf
c
      parameter(
     $    plnlv=plon*plev,      !  length of multilevel ht field
     $    pln2lv=plonp2*plev,   !  length of multilevel buffer field
     $    plndlv=plond*plev,       ! length of multilevel 3-d field
     $    plat2=plat/2)         !
C
c
      parameter(
     $    plevmx=4,                        ! number of subsurface levels
     $    pscal=4,                         ! number of passive scalars
     $    ptflen=((3+2*pcnst)*plev + (9+pcnst))*plond,! len time dep fld
     $    pdflen=(4*plev + 7)*plond)       ! length of diagnostic fields
      parameter(
     $    pmulti=22 + pcnst*5,            ! number of multilevel fields
c                                         ! primary history tape
     $    psingl=35+plevmx,               ! number of single level field
     $    ptplen=pmulti+psingl)           ! total fields to primary tape
      parameter(
     $    pbflnb=2*plndlv+ptflen+pdflen,   ! length of nprgtl buffer
     $    pbflna=(4+6*plev)*plond,        ! length of nprg
     $    phtbf=(pmulti*plev+psingl)*plon, ! length of history tape buf
     $    pflenb=((pbflnb+pbflna+phtbf)/512+1)*512,  !  padded length of
c                                                   !   + nprgm1 + hbuf
     $    pflena=(pbflna/512+1)*512,       ! padded length of nprg
     $    plnbuf=pflenb+pflena,           ! total buffer length
     $    pflds=14,            ! number of fields on initial dataset
     $    ptapes=1)            ! maximum number of history tapes allowed
c
c     !       history tape header parameters
c
      parameter (
     $   prbd=3,                           ! records before data on ht
     $   pftyp=42,                         ! units digit of ht format ty
     $   plenhi=37+3*ptplen,               ! length of integer header re
     $   plenhc=76+2*ptplen,               ! length of character header
     $   plenbr=3*(2*plev+1)+2*plat,       ! length of real header buffe
     $   plenhr=plenbr)                    ! length of real header recor
c---------------------------------------------------------------------
C
C Parameters related to spectral domain
C
      integer
     $    ptrm,ptrn,ptrk,pcray,pmax,pmaxp,pnmax,pmmax,plev2,
     $    par0,par2,pspt,psp,pspl
C
C Truncation parameters
C
      parameter(
     $    ptrm=42,              !  m truncation parameter (T42)
     $    ptrn=42,              !  n truncation parameter
     $    ptrk=42,              !  k truncation parameter
     $    pcray=64,             !  machine word size
C
C Dimensions for spectral arrays
C
     $    pmax=ptrn+1,
     $    pmaxp=pmax+1,
     $    pnmax=ptrk+1,
     $    pmmax=ptrm+1,
     $    plev2=plev*plev)
C
C Intermediate parameters
C
      parameter(
     $    par0=ptrm+ptrn-ptrk,
     $    par2=par0*(par0+1)/2)
C
C Dimensions for Legendre arrays
C
      parameter(
     $    pspt=(ptrn+1)*pmmax-par2,
     $    psp=2*pspt,
     $    pspl=psp*plev)
C
c---------------------------------------------------------------------
C
C FLAGS ASSOCIATED WITH ANNUAL CYCLE
C
C IF ANNCYC=.TRUE., SEA SURFACE TEMPERATURES
C AND OZONE PATH LENGTHS GO THROUGH AN ANNUAL
C CYCLE DETERMINED BY SPLINE COEFFICIENTS FOR
C EACH GRID POINT
C
C IF ANNCYC=.FALSE.,SST'S AND O3 PATH LENGTHS
C CONSTANT IN TIME
C
C IF HYDRO =.TRUE., HYDROLOGIC CYCLE IN MODEL
C WILL BE ENABLED (REQUIRES ANNCYC TO BE .TRUE.)
C
C ITSST IS SEA SURFACE TEMPERATURE UPDATE FREQUENCY (ITERS)
C
      common/comanc/anncyc,hydro,itsst
      logical       anncyc,hydro
      integer itsst
c---------------------------------------------------------------------
C
CL            HISTORY TAPE HEADER, RECORD NUMBER 1 (INTEGER VALUES)
CL            PROVIDES DATA FILE DESCRIPTION FOR HISTORY TAPE OUTPUT
C
      COMMON /COMHDI/
     +   LENHDI(1), MFTYP,    MFILH,    MFILTH,     NRBD
     +  ,MAXSIZ,    NDAVU,    MSPHER,   MLON,       NLONW
     +  ,MOREC,     MLEV,     MTRM,     MTRN,       MTRK
     +  ,NFLDH,     NSTEPH,   NSTPRH,   NITSLF
     +  ,NDBASE,   NSBASE,    NDCUR,    NSCUR
     +  ,NBDATE,   NBSEC,     NCDATE,   NCSEC
     +  ,MDT,      MHISF,     MFSTRT,   LENHDC,     LENHDR
     +  ,MPSIG,    MPLAT,     MPWTS,    MPFLDS,     MPCFLD
C             ARRAY OF INTEGER FIELD LIST INFORMATION
     $  ,mflds(3,ptplen)
C
CL            HISTORY TAPE HEADER, RECORD NUMBER 2 (CHARACTER VALUES)
C
      COMMON /COMHDC/
     +   MCASE,    MCSTIT
     +  ,LNHSTC,   LDHSTC,   LTHSTC,   LSHSTC
     +  ,LNHSTP,   LDHSTP,   LTHSTP,   LSHSTP
     +  ,LNHSTF,   LDHSTF,   LTHSTF,   LSHSTF
     +  ,LNHSTI,   LDHSTI,   LTHSTI,   LSHSTI
     +  ,LNHSTA,   LDHSTA,   LTHSTA,   LSHSTA
C             ARRAY OF CHARACTER FIELD LIST INFORMATION
     $  ,mcflds(2,ptplen)
      CHARACTER*8 MCASE
     + ,LDHSTC,   LTHSTC,   LSHSTC
     + ,LDHSTP,   LTHSTP,   LSHSTP
     + ,LDHSTF,   LTHSTF,   LSHSTF
     + ,LDHSTI,   LTHSTI,   LSHSTI
     + ,LDHSTA,   LTHSTA,   LSHSTA
     + ,MCFLDS
      CHARACTER*80 MCSTIT
     + ,LNHSTC,   LNHSTP,   LNHSTF,   LNHSTI,   LNHSTA
C
CL            HISTORY TAPE HEADER, RECORD NUMBER 3 (REAL VALUES)
C
      COMMON /COMHDR/
C             BUFFER BUFHD CONTAINS SIGMA,LATITUDE AND GAUSSIAN WEIGHTS
     +   REALHD(1),BUFHD(PLENBR)
C
c     !  Store certain header values for all history tapes
c
      common /comhdt/
     $       nstepht(ptapes)
      integer nstepht
C
      common /comhdtc/
     $       ldhstct(ptapes), lthstct(ptapes), lshstct(ptapes)
C
      character*8 ldhstct, lthstct, lshstct
c---------------------------------------------------------------------
C
C Semi-implicit timestep constants
C
      common/comimp/
     $     t0(plev),  bm1(plev2*pnmax),  tau(plev2),  aq(plev2),
     $     dtime,     twodt,        eps
      real t0, bm1, tau, aq, dtime, twodt, eps
C
      real xtau(plev,plev)
      equivalence (xtau,tau)
C
C dimension bm1(nlev*nlev*nmax)
C
c---------------------------------------------------------------------
C
C LOGICAL UNIT NUMBERS
C
      common/comlun/nsds    ,nsre    ,nsre1   ,nra1    ,
     $              nrb1    ,nprint  ,nout    ,nread   ,
     $              ninit   ,nsdv    ,nozone  ,nsst    ,
     $              nalb    ,nabem   ,lutag(99)
C
      integer nsds,    !  restart dataset unit
     $        nsre,    !  primary regeneration dataset unit
     $        nsre1,   !  secondary regeneration dataset
     $        nra1,    !  a work file
     $        nrb1,    !  b work file
     $        nprint,  !  alternate print unit
     $        nout,    !  print unit
     $        nread,   !  input unit
     $        ninit,   !  initial dataset unit
     $        nsdv,    !  standard deviation dataset
     $        nozone,  !  ozone dataset
     $        nsst,    !  sst dataset
     $        nalb,    !  albedo dataset
     $        nabem    !  absorptivity/emissivity dataset
      logical lutag
c---------------------------------------------------------------------
C
C model control time variables
C
      integer
     $      nrstrt  ,nstep   ,nstepr  ,nestep  ,nstop   ,
     $      mdbase  ,msbase  ,mdcur   ,mscur   ,
     $      mbdate  ,mbsec   ,mcdate  ,mcsec   ,
     $      nndbas  ,nnsbas  ,nnbdat  ,nnbsec
      real  calday
c
      common /comtim/
     $      calday  ,nrstrt  ,nstep   ,nstepr  ,nestep  ,nstop   ,
     $      mdbase  ,msbase  ,mdcur   ,mscur   ,
     $      mbdate  ,mbsec   ,mcdate  ,mcsec   ,
     $      nndbas  ,nnsbas  ,nnbdat  ,nnbsec
C
c---------------------------------------------------------------------
c
c surface albedo data
c
c vs = 0.2 - 0.7 micro-meters wavelength range
c ni = 0.7 - 5.0 micro-meters wavelength range
c
c s  = strong zenith angle dependent surfaces
c w  = weak   zenith angle dependent surfaces
c
c the albedos are computed for a model grid box by ascribing values to
c high resolution points from a vegetation dataset, then linearlly 
c averaging to the grid box value; ocean and land values are averaged
c together along coastlines; the fraction of every grid box that has
c strong zenith angle dependence is included also.
c
      real albvss, ! grid box alb for vs over strng zn srfs
     $     albvsw, ! grid box alb for vs over weak  zn srfs
     $     albnis, ! grid box alb for ni over strng zn srfs
     $     albniw, ! grid box alb for ni over weak  zn srfs
     $     frctst  ! fraction of area in grid box strng zn 
c
      common /crdalb/albvss(plond,plat),albvsw(plond,plat),
     $               albnis(plond,plat),albniw(plond,plat),
     $               frctst(plond,plat)
c
c
c surface boundary data
c
c vegtyp is used to specify the thermal properites of the surface,
c as well as determine the location of permanent land ice points;
c it is the dominant surface type within the model grid box based
c on a high resolution vegetation type dataset;
c it encoded in the following manner:
c
c   1        ocean
c   2        sea ice
c   3        permanent land ice
c   4        tropical evergreen forest
c   5        deciduous forest
c   6        grassland/tundra
c   7        desert
c
c rghnss is the aerodynamic roughness length for the grid box, computed
c by linear averaging of the values ascribed to high resolution 
c vegetation dataset values; ocean and land values are averaged together
c at coastlines.
c
c evapf is the ratio of actual to potential evaporation, and is computed
c from the high resolution vegetation dataset in a manner similar to the
c aerodynamic roughness.
c
c vevapf allows for variable snow cover, where the underlying
c evaporability factor is modified; see radalb.
c
c snwjan and snwjly are mean climatological snow depths (liquid water
c equivalent) used to compute the actual daily values of snow cover.
c
      real vegtyp, ! surface thermal type, based on veg type
     $     rghnss, ! aerodynamic roughness length
     $     evapf , ! constant surface evaporability
     $     vevapf, ! variable surface evaporability
     $     snwjan, ! snow cover (liq water equiv) for january
     $     snwjly  ! snow cover (liq water equiv) for july
c
      common /crdsrf/vegtyp(plond,plat),rghnss(plond,plat),
     $               evapf (plond,plat),vevapf(plond,plat),
     $               snwjan(plond,plat),snwjly(plond,plat)
c
c---------------------------------------------------------------------
c
c water vapor narrow band constants for lw computations
c
      real realk,st,a1,a2,b1,b2,
     $     coefa,coefb,coefc,coefd,
     $     coefe,coeff,coefg,coefh,
     $     coefi,coefj,coefk,
     $     c1,c2,c3,c4,c5,c6,c7,
     $     c8 ,c9 ,c10,c11,c12,c13,c14,c15,c16,c17,
     $     c18,c19,c20,c21,c22,c23,c24,c25,c26,c27,
     $     c28,c29,c30,c31,
     $     fwcoef,fwc1,fwc2,fc1,cfa1
c
      common/crdcae/realk(2), st(2), a1(2), a2(2), b1(2), b2(2),
c
c constant coefficients for water vapor absorptivity and emissivi
c
     $              coefa(3,4),coefb(4,4),coefc(3,4),coefd(4,4),
     $              coefe(3,4),coeff(6,2),coefg(2,4),coefh(2,4),
     $              coefi(6,2),coefj(3,2),coefk(3,2),
     $              c1(4),c2(4),c3(4),c4(4),c5(4),c6(4),c7(4),
     $              c8 ,c9 ,c10,c11,c12,c13,c14,c15,c16,c17,
     $              c18,c19,c20,c21,c22,c23,c24,c25,c26,c27,
     $              c28,c29,c30,c31,
c
c farwing correction constants for narrow-band emissivity model
c introduce farwing correction to account for the
c deficiencies in narrow-band model used to derive the
c emissivity. tuned with arkings line calculations.
c
     $              fwcoef,fwc1,fwc2,fc1,cfa1
c
c---------------------------------------------------------------------
c
c radiation constants
c
      real gravit,   ! gravitational acceleration
     $        rga,   ! 1 over gravit
     $      cpair,   ! heat capacity air at constant pressure
     $     epsilo,   ! ratio mmw h2o to mmw air
     $       sslp,   ! standard pressure
     $     stebol,   ! stephan boltzmann constant
     $     rgsslp,   ! 0.5 / (gravit*sslp)
     $     co2vmr,   ! co2 volume mixing ratio
     $      dpfo3,   ! Doppler factor for o3
     $     dpfco2,   ! Doppler factor for co2
     $     dayspy,   ! solar days in one year
     $        pie    ! pie
c
      common/crdcon/gravit,    rga,  cpair,  epsilo,   sslp,
     $              stebol, rgsslp, co2vmr,   dpfo3, dpfco2,
     $              dayspy,    pie
c
c---------------------------------------------------------------------
c
c ozone mixing ratios, pressures and times
c
c
      integer noz,  ! number of ozone data levels
     $      nozp1,  ! number of ozone data levels plus one
     $        koz,  ! number of levels of zone data
     $       lmoz,  ! first month index for ozone data storage
     $       nmoz,  ! second month index for ozone data storage
     $   monthloz,  ! first month index for ozone data read in
     $   monthnoz   ! second month index for ozone data read in
c
      real  cplos,  ! conversion factor for ozone path lengths
     $      cplol,  ! conversion factor for ozone prs weighted path lengths
     $      ozday,  ! day of year for each month of ozone data
     $      ozmix,  ! time interpolated zonal mean lat/height o3 vol mx ratio
     $     ozmixm,  ! two months of ozone data that bracket current date
     $      plol0,  ! prs wghtd ozone path lengths over original prs lvls (cm)
     $      plos0,  ! ozone path lengths between original prs levels (cm)
     $        pin,  ! interface pressures for input ozone data
     $        poz   ! level pressures for input ozone data
c
      parameter (noz=40)
      parameter (nozp1=noz+1) 
c
      common /crdozp/  cplos,    cplol,       
     $                   koz,     lmoz,        nmoz,
     $              monthloz, monthnoz, ozday(0:13),
     $       ozmix(plat,noz),    ozmixm(plat,noz,2),
     $     plol0(plat,nozp1),     plos0(plat,nozp1),
     $            pin(nozp1),  poz(noz)
c
c---------------------------------------------------------------------
c
c input arguments
c
      real gravx,       ! gravitational acceleration (m/s**2)
     $     cpairx,      ! heat capacity dry air at constant prs (J/kg/K)
     $     epsilox,     ! ratio mean mol weight h2o to dry air
     $     stebolx      ! Sefan-Boltzmann constant (W/m**2/K**4)
c
c---------------------------------------------------------------------
c
c local workspace
c
      integer   iband   ! h2o band index
C
C /CRDCAE/
C H2O EMISSIVITY AND ABSORTIVITY COEFFICIENTS
C
      data coefd/7.03047e-01,-2.63501e-03,-1.57023e-06,0.0,
     $           5.29269e-01,-3.14754e-03, 4.39595e-06,0.0,
     $           7.88193e-02, 1.31290e-03, 4.25827e-06,-1.23982e-08,
     $           1.62744e-01, 2.22847e-03, 2.60102e-06,-4.30133e-08/
C
      data coefb/8.85675e+00,-3.51620e-02, 2.38653e-04,-1.71439e-06,
     $           5.73841e+00,-1.91919e-02, 1.65993e-04,-1.54665e-06,
     $           6.64034e+00, 1.56651e-02,-9.73357e-05, 0.0,
     $           7.09281e+00, 1.40056e-02,-1.15774e-04, 0.0/
C
      data coefe/3.93137e-02,-4.34341e-05,3.74545e-07,
     $           3.67785e-02,-3.10794e-05,2.94436e-07,
     $           7.42500e-02, 3.97397e-05,0.0,
     $           7.52859e-02, 4.18073e-05,0.0/
C
      data coefa/1.01400e+00,6.41695e-03,2.85787e-05,
     $           1.01320e+00,6.86400e-03,2.96961e-05,
     $           1.02920e+00,1.01680e-02,5.30226e-05,
     $           1.02743e+00,9.85113e-03,5.00233e-05/
C
      data coefc/9.90127e-01,1.22475e-03,4.90135e-06,
     $           9.89753e-01,1.97081e-03,3.42046e-06,
     $           9.75230e-01,1.03341e-03,0.0,
     $           9.77366e-01,8.60014e-04,0.0/
C
      data coeff/2.2037 e-01,1.39719e-03,-7.32011e-06,
     $          -1.40262e-08,2.13638e-10,-2.35955e-13,
     $           3.07431e-01,8.27225e-04,-1.30067e-05,
     $           3.49847e-08,2.07835e-10,-1.98937e-12/
C
      data coefg/9.04489e+00,-9.56499e-03,
     $           1.80898e+01,-1.91300e-02,
     $           8.72239e+00,-9.53359e-03,
     $           1.74448e+01,-1.90672e-02/
C
      data coefh/5.46557e+01,-7.30387e-02,
     $           1.09311e+02,-1.46077e-01,
     $           5.11479e+01,-6.82615e-02,
     $           1.02296e+02,-1.36523e-01/
C
      data coefi/3.31654e-01,-2.86103e-04,-7.87860e-06,
     $           5.88187e-08,-1.25340e-10,-1.37731e-12,
     $           3.14365e-01,-1.33872e-03,-2.15585e-06,
     $           6.07798e-08,-3.45612e-10,-9.34139e-15/
C
      data coefj/2.82096e-02,2.47836e-04,1.16904e-06,
     $           9.27379e-02,8.04454e-04,6.88844e-06/
C
      data coefk/2.48852e-01,2.09667e-03,2.60377e-06,
     $           1.03594e+00,6.58620e-03,4.04456e-06/
C
C NARROW BAND DATA FOR H2O
C 200CM DATA FOR 800-1000 CM-1 AND 1000-1200 CM-1.
C
      data realk/  0.18967069430426e-04, 0.70172244841851e-04   /
      data   st /  0.31930234492350e-03, 0.97907319939060e-03   /
      data   a1 /  0.28775403075736e-01, 0.23236701470511e-01   /
      data   a2 / -0.57966222388131e-04,-0.95105504388411e-04   /
      data   b1 /  0.29927771523756e-01, 0.21737073577293e-01   /
      data   b2 / -0.86322071248593e-04,-0.78543550629536e-04   /
c
c---------------------------------------------------------------------
c
c set general radiation constants; convert to cgs units where appropriate:
c
      gravit  =  100.*gravx  
      rga     =  1./gravit   
      cpair   =  1.e4*cpairx
      epsilo  =  epsilox
      sslp    =  1.013250e6
      stebol  =  1.e3*stebolx    
      rgsslp  =  0.5/(gravit*sslp)
      co2vmr  =  3.3e-4          
      dpfo3   =  2.5e-3      
      dpfco2  =  5.0e-3      
      dayspy  =  365.
      pie     =  4.*atan(1.)
c
c coefficients for h2o emissivity and absortivity.
c
      do 10 iband=1,4
         c1(iband) = coefe(3,iband)/coefe(2,iband)
         c2(iband) = coefb(3,iband)/coefb(2,iband)
         c3(iband) = coefb(4,iband)/coefb(3,iband)
         c4(iband) = coefd(3,iband)/coefd(2,iband)
         c5(iband) = coefd(4,iband)/coefd(3,iband)
         c6(iband) = coefa(3,iband)/coefa(2,iband)
         c7(iband) = coefc(3,iband)/coefc(2,iband)
   10 continue
      c8   = coeff(3,1)/coeff(2,1)
      c9   = coeff(3,2)/coeff(2,2)
      c10  = coeff(4,1)/coeff(3,1)
      c11  = coeff(4,2)/coeff(3,2)
      c12  = coeff(5,1)/coeff(4,1)
      c13  = coeff(5,2)/coeff(4,2)
      c14  = coeff(6,1)/coeff(5,1)
      c15  = coeff(6,2)/coeff(5,2)
      c16  = coefj(3,1)/coefj(2,1)
      c17  = coefk(3,1)/coefk(2,1)
      c18  = coefi(3,1)/coefi(2,1)
      c19  = coefi(3,2)/coefi(2,2)
      c20  = coefi(4,1)/coefi(3,1)
      c21  = coefi(4,2)/coefi(3,2)
      c22  = coefi(5,1)/coefi(4,1)
      c23  = coefi(5,2)/coefi(4,2)
      c24  = coefi(6,1)/coefi(5,1)
      c25  = coefi(6,2)/coefi(5,2)
      c26  = coefj(3,2)/coefj(2,2)
      c27  = coefk(3,2)/coefk(2,2)
      c28  = .5
      c29  = .002053
      c30  = .1
      c31  = 3.0e-5
      cfa1 = .61
c
c initialize further longwave constants referring to 
c far wing correction:
c
      fwcoef  = .1
      fwc1    = .30
      fwc2    = 4.5
      fc1     = 2.6
c
c...... remove surface data and ozone io statements:
c
c
c done
c
      return
      end
      subroutine radctl(nrow   ,clat   ,oro    ,sndpth ,ts     ,
     $                  tg     ,ps     ,pmid   ,pint   ,pmln   ,
     $                  piln   ,t      ,h2ommr ,cld    ,effcld ,
     $                  clwp   ,plol   ,plos   ,solin  ,sabtp  ,
     $                  frsa   ,clrst  ,clrss  ,qrs    ,firtp  ,
     $                  frla   ,clrlt  ,clrls  ,qrl    ,srfrad ,
     &     fupsw,fdsw,fuplw,fdlw,fupcbi,fdncti,
     &     geff,aeff,taueff,crecz,swfdiv,lwfdiv,
     &     bccmil,tccmil,numlay,ncclev,
     &     cbi,cti,ncldlay,
     &     idebug)
c-----------------------------------------------------------------------
c
c main computational entry for radiation computation
c
c computations done for one latitude line
c
c radiation uses cgs units, so conversions must be done from 
c model fields to radiation fields
c
c
c calling sequence:
c
c     .
c
c    radini          initializes radiation constants, ozone, 
c                    and surface data
c     .
c
c    radoz1          updates ozone and does time interpolation
c
c     .
c
c    phys            calls physics routines
c
c      radctl        interface for radiation scheme
c
c        radoz2      interpolates ozone paths to model interface pressures
c
c        radinp      converts units of model fields and computes ozone
c                    mixing ratio for solar scheme
c
c        radcsw      performs solar computation
c
c          radalb    computes surface albedos
c          radded    computes delta-Eddington solution 
c          radclr    computes diagnostic clear sky fluxes
c
c        radclw      performs longwave computation
c
c          radtpl    computes path quantities
c          radems    computes emissivity
c          radabs    computes absorptivity
c
c        radout      converts radiation fluxes to mks units; computes
c                    surface radiative flux for surface computations
c
c
c-----------------------------------------------------------------------
      implicit none
c-----------------------------------------------------------------------
c
c     !   Basic grid point resolution parameters
c
      integer plev,plon,plat,plonp2,plevp,pcnst,plond,platd,plevd
      integer nxpt,jintmx,i1,j1,j1m
c
      parameter(plev=  133,        ! number of vertical levels
     $          plevp=plev+1,     ! number of vertical interfaces
     $          plon=1,           ! number of longitudes (T42)
     $          plat=1,           ! number of latitudes (T42)
     $          plonp2=plon+2,    ! length necessary for fft
     $          pcnst=1,          ! number of constituents
     $          nxpt=1,           !
     $          jintmx=2,         !
     $          plond=plon + 1 + 2*nxpt, ! slt extended domain longitude
     $          platd=plat + 2*nxpt + 2*jintmx, ! slt extended domain lat.
     $          plevd=plev*(3+pcnst), ! fold plev,pcnst indices into one
     $          i1=1+nxpt,        ! model starting longitude (3-d)
     $          j1=1+nxpt+jintmx, ! model starting latitude (3-d)
     $          j1m=j1-1)         ! model starting offset (3-d)
c
c-----------------------------------------------------------------------
C
C FLAGS ASSOCIATED WITH ANNUAL CYCLE
C
C IF ANNCYC=.TRUE., SEA SURFACE TEMPERATURES
C AND OZONE PATH LENGTHS GO THROUGH AN ANNUAL
C CYCLE DETERMINED BY SPLINE COEFFICIENTS FOR
C EACH GRID POINT
C
C IF ANNCYC=.FALSE.,SST'S AND O3 PATH LENGTHS
C CONSTANT IN TIME
C
C IF HYDRO =.TRUE., HYDROLOGIC CYCLE IN MODEL
C WILL BE ENABLED (REQUIRES ANNCYC TO BE .TRUE.)
C
C ITSST IS SEA SURFACE TEMPERATURE UPDATE FREQUENCY (ITERS)
C
      common/comanc/anncyc,hydro,itsst
      logical       anncyc,hydro
      integer itsst
c-----------------------------------------------------------------------
C
C model control time variables
C
      integer
     $      nrstrt  ,nstep   ,nstepr  ,nestep  ,nstop   ,
     $      mdbase  ,msbase  ,mdcur   ,mscur   ,
     $      mbdate  ,mbsec   ,mcdate  ,mcsec   ,
     $      nndbas  ,nnsbas  ,nnbdat  ,nnbsec
      real  calday
c
      common /comtim/
     $      calday  ,nrstrt  ,nstep   ,nstepr  ,nestep  ,nstop   ,
     $      mdbase  ,msbase  ,mdcur   ,mscur   ,
     $      mbdate  ,mbsec   ,mcdate  ,mcsec   ,
     $      nndbas  ,nnsbas  ,nnbdat  ,nnbsec
C
c-----------------------------------------------------------------------
c
c radiation control variables
c
c fradsw = .t. iff full shortwave computation
c fradlw = .t. iff full longwave computation
c
c irad = iteration frequency for radiation computation
c iradae = iteration frequency for absorptivity/
c emissivity computation
c
c
      integer iradae,irad,naclw,nacsw,fnlw,fnsw
      logical aeres
c
      common/crdctl/iradae, irad, naclw, nacsw, fnlw, fnsw, aeres
c
c-----------------------------------------------------------------------
c
c input arguments
c
      integer nrow             ! latitude row index 
      real clat,               ! Current latitude (radians)
     $     oro(plond),         ! land/ocean/sea ice flag
     $     sndpth(plond),      ! snow depth (liquid water equivalent)
     $     ts(plond),          ! surface air temperature
     $     tg(plond),          ! surface (skin) temperature
     $     ps(plond),          ! surface pressure
     $     pmid(plond,plev),   ! model level pressures
     $     pint(plond,plevp),  ! model interface pressures
     $     pmln(plond,plev),   ! natural log of pmid
     $     piln(plond,plevp),  ! natural log of pint
     $     t(plond,plev),      ! model level temperatures
     $     h2ommr(plond,plev), ! model level specific humidity
     $     cld(plond,plevp),   ! fractional cloud cover
     $     effcld(plond,plevp),! effective fractional cloud cover
     $     clwp(plond,plev),   ! cloud liquid water path
     $     plol(plond,plevp),  ! o3 pressure weighted path lengths (cm)
     $     plos(plond,plevp)   ! o3 path lengths (cm)
c
c output arguments
c
c output solar 
c
      real solin(plond),      ! solar incident flux
     $     sabtp(plond),      ! total column absorbed solar flux
     $     frsa(plond),       ! surface absorbed solar flux
     $     clrst(plond),      ! clr sky total column abs solar flux
     $     clrss(plond),      ! clr sky surface abs solar flux
     $     qrs(plond,plev)    ! solar heating rate
c
c output longwave
c
      real firtp(plond),      ! net outgoing lw flx at model top
     $     frla(plond),       ! srf longwave cooling (up-dwn) flux
     $     clrlt(plond),      ! clr sky lw flx at model top
     $     clrls(plond),      ! clr sky lw flx at srf (up-dwn)
     $     qrl(plond,plev),   ! longwave cooling rate
     $     slwd(plond)        ! srf down longwave flux
c
c surface radiative heating
c
      real srfrad(plond)       ! srf radiative heat flux
c
c local workspace
c
      real pbr(plond,plev),    ! model mid-level pressures (dynes/cm2)
     $     pnm(plond,plevp),   ! model interface pressures (dynes/cm2)
     $     o3mmr(plond,plev),  ! ozone mass mixing ratio
     $     plco2(plond,plevp), ! prs weighted co2 path
     $     plh2o(plond,plevp), ! prs weighted h2o path
     $     tclrsf(plond,plevp),! total clear sky fraction, level to space
     $     coszrs(plond),      ! cosine solar zenith angle
     $     eccf                ! earth/sun distance factor
c
      integer  jlat,           ! north/south latitude index
     $         norsou
      external norsou,         ! sets north/south latitude index from nrow
     $         radinp,         ! computes latitude dependent radiation input
     $         radcsw,         ! computes solar radiation
     $         radclw,         ! computes longwave radiation
     $         radout          ! sets radiation output and converts units
c------------------------------------------------------------------------
c
c     zender cirrus cloud hack workspace
      integer nspint  ! number of spectral intervals across solar spectrum
      parameter ( nspint = 18 )
c
      integer 
     &     idebug,    !debugging level from cloud model
     &     bccmil,    !bottom ccm2 interface level
     &     tccmil,    !top ccm2 interface level
     &     numlay,    !num_layer in cloud model
     &     ncclev,    !num_ccm2_level (with overlaps already subtracted)
     &     cbi,       !cloud bottom index
     &     cti,       !cloud top index
     &     ncldlay    !number of layers in cloud itself
C
      real
     &     geff(0:nspint+1,0:numlay+1),   !effective asymmetry factor
     &     aeff(0:nspint+1,0:numlay+1),   !effective singel scat albedo
     &     taueff(0:nspint+1,0:numlay+1), !effective extinction optical depth
     &     swfdiv(0:nspint+1,0:numlay+1),  !SW spectral flux divergence
     &     lwfdiv(0:numlay+1), !LW flux divergence
     &     crecz(0:numlay+1), !effective radius of particles (meters)
     &     fupsw(0:numlay+1), !total shortwave upward directed flux
     &     fdsw(0:numlay+1), !total shortwave downward directed flux
     &     fuplw(0:numlay+1), !total longwave upward directed flux
     &     fdlw(0:numlay+1), !total longwave downward directed flux
     &     fupcbi, !longwave upward directed flux at cloud base interface
     &     fdncti !longwave downward directed flux at cloud top interface
c
c     end of zender workspace alterations
c--------------------------------------------------------------------------
c
c latitude index from north to south
c
      jlat = norsou(nrow) + 1
c
c set ozone
c
      if(anncyc) then
         if(nstep.eq.0 .or.
     $      (mod(nstep-1,iradae).eq.0 .and. nstep.ne.1)) then
            call radoz2(jlat,pint,plol,plos)
         end if
      else 
         if(nstep.eq.0) call radoz2(jlat,pint,plol,plos)
      end if
c
c set latitude dependent radiation input
c
      call radinp(clat    ,pmid  ,pint   ,h2ommr ,cld    ,
     $            plos    ,pbr   ,pnm    ,plco2  ,plh2o  ,
     $            tclrsf  ,eccf  ,coszrs ,o3mmr,
     &     bccmil,tccmil,numlay,ncclev,
     &     cbi,cti,ncldlay,
     &     idebug)
c
c solar radiation computation
c
      call radcsw(jlat    ,oro   ,sndpth ,pnm    ,h2ommr ,
     $            cld     ,clwp  ,o3mmr  ,eccf   ,coszrs ,
     $            solin   ,qrs   ,frsa   ,sabtp  ,clrss  ,
     $            clrst   ,
     &     fupsw,fdsw,
     &     geff,aeff,taueff,crecz,swfdiv,
     &     bccmil,tccmil,numlay,ncclev,
     &     cbi,cti,ncldlay,
     &     idebug)
c
c longwave radiation computation
c
      call radclw(nrow   ,tg     ,plol   ,plos   ,t      ,
     $            h2ommr ,pbr    ,pnm    ,pmln   ,piln   ,
     $            plco2  ,plh2o  ,effcld ,tclrsf ,qrl    ,
     $            frla   ,firtp  ,clrls  ,clrlt  ,slwd   ,
     &     fuplw,fdlw,fupcbi,fdncti,
     &     lwfdiv,
     &     bccmil,tccmil,numlay,ncclev,
     &     cbi,cti,ncldlay,
     &     idebug)
c
c set output radiation fields
c
      call radout(solin  ,sabtp  ,frsa   ,clrst  ,clrss  ,
     $            qrs    ,firtp  ,frla   ,clrlt  ,clrls  ,
     $            qrl    ,slwd   ,srfrad )
c
c done
c
      return
      end

      subroutine radinp(clat    ,pmid    ,pint    ,h2ommr ,cld   ,
     $                  plos    ,pmidrd  ,pintrd  ,plco2  ,plh2o ,
     $                  tclrsf  ,eccf    ,coszrs  ,o3mmr,
     &     bccmil,tccmil,numlay,ncclev,
     &     cbi,cti,ncldlay,
     &     idebug)
c---------------------------------------------------------------------
c
c sets latitude and time dependent arrays for input to radiation
c
c works on a longitude/height slice of data for each latitude band
c
c computes solar input (earth-sun distance factor and cosine solar
c zenith angle), converts model pressures to cgs, computes 
c path length arrays needed for the longwave radiation code,
c and computes ozone mixing ratio.
c
c---------------------------------------------------------------------
      implicit none
c---------------------------------------------------------------------
c
c     !   Basic grid point resolution parameters
c
      integer plev,plon,plat,plonp2,plevp,pcnst,plond,platd,plevd
      integer nxpt,jintmx,i1,j1,j1m
c
      parameter(plev=  133,        ! number of vertical levels
     $          plevp=plev+1,     ! number of vertical interfaces
     $          plon=1,           ! number of longitudes (T42)
     $          plat=1,           ! number of latitudes (T42)
     $          plonp2=plon+2,    ! length necessary for fft
     $          pcnst=1,          ! number of constituents
     $          nxpt=1,           !
     $          jintmx=2,         !
     $          plond=plon + 1 + 2*nxpt, ! slt extended domain longitude
     $          platd=plat + 2*nxpt + 2*jintmx, ! slt extended domain lat.
     $          plevd=plev*(3+pcnst), ! fold plev,pcnst indices into one
     $          i1=1+nxpt,        ! model starting longitude (3-d)
     $          j1=1+nxpt+jintmx, ! model starting latitude (3-d)
     $          j1m=j1-1)         ! model starting offset (3-d)
c
c---------------------------------------------------------------------
C
C model control time variables
C
      integer
     $      nrstrt  ,nstep   ,nstepr  ,nestep  ,nstop   ,
     $      mdbase  ,msbase  ,mdcur   ,mscur   ,
     $      mbdate  ,mbsec   ,mcdate  ,mcsec   ,
     $      nndbas  ,nnsbas  ,nnbdat  ,nnbsec
      real  calday
c
      common /comtim/
     $      calday  ,nrstrt  ,nstep   ,nstepr  ,nestep  ,nstop   ,
     $      mdbase  ,msbase  ,mdcur   ,mscur   ,
     $      mbdate  ,mbsec   ,mcdate  ,mcsec   ,
     $      nndbas  ,nnsbas  ,nnbdat  ,nnbsec
C
c---------------------------------------------------------------------
c
c radiation constants
c
      real gravit,   ! gravitational acceleration
     $        rga,   ! 1 over gravit
     $      cpair,   ! heat capacity air at constant pressure
     $     epsilo,   ! ratio mmw h2o to mmw air
     $       sslp,   ! standard pressure
     $     stebol,   ! stephan boltzmann constant
     $     rgsslp,   ! 0.5 / (gravit*sslp)
     $     co2vmr,   ! co2 volume mixing ratio
     $      dpfo3,   ! Doppler factor for o3
     $     dpfco2,   ! Doppler factor for co2
     $     dayspy,   ! solar days in one year
     $        pie    ! pie
c
      common/crdcon/gravit,    rga,  cpair,  epsilo,   sslp,
     $              stebol, rgsslp, co2vmr,   dpfo3, dpfco2,
     $              dayspy,    pie
c
c---------------------------------------------------------------------
c
c input arguments
c
      real clat,               ! current latitude (radians)
     $     pmid(plond,plev),   ! pressure at model mid-levels (Pascals)
     $     pint(plond,plevp),  ! pressure at model interfaces (pascals)
     $     h2ommr(plond,plev), ! model moisture field
     $     cld(plond,plevp),   ! fractional cloud cover
     $     plos(plond,plevp)   ! o3 path length (cm)
c
c output arguments
c
      real pmidrd(plond,plev), ! pressure at model mid-levels (dynes/cm*2)
     $     pintrd(plond,plevp),! pressure at model interfaces (dynes/cm*2)
     $     plco2(plond,plevp), ! vert. pth lngth of co2 (prs-weighted)
     $     plh2o(plond,plevp), ! vert. pth lngth of h2o vap. (prs-weighted)
     $     tclrsf(plond,plevp) ! product of clr-sky fractions from top
c                              ! of atmos. to lvl k.
      real eccf,               ! earth-sun distance factor
     $     coszrs(plond),      ! cosine solar zenith angle
     $     o3mmr(plond,plev)   ! ozone mass mixing ratio
c
c local variables
c
      integer   i,  ! longitude loop index
     $          k   ! vertical loop index
c
      real    phi,  ! greenwich calendar day + local time + longitude offset
     $      theta,  ! earth orbit seasonal angle in radians
     $      delta,  ! solar declination angle  in radians
     $       sinc,  ! sine   of latitude
     $       cosc,  ! cosine of latitude
     $       sind,  ! sine   of declination
     $       cosd,  ! cosine of declination
     $         v0,  ! volume of a gas at stp (cm**3/mol)
     $         p0,  ! standard pressure (dynes/cm**2)
     $        amd,  ! effective molecular weight of dry air (g/mol)
     $        amo,  ! molecular weight of ozone (g/mol)
     $      amco2,  ! molecular weight of co2   (g/mol)
     $        cpl,  ! constant in ozone path length to mixing ratio
     $      cpwpl,  ! constant in co2 mixing ratio to path length conversion
     $       vmmr   ! ozone volume mixing ratio
c
      data v0    /  22413.6   /
      data p0    /  1.01325e6 /
      data amd   /  28.9644   /
      data amo   /  48.0000   /
      data amco2 /  44.0000   /
c
c-----------------------------------------------------------------------
c
c------------------------------------------------------------------------
c     zender cirrus cloud hack workspace
c
      integer 
     &     idebug,    !debugging level from cloud model
c     &     kcz,       !index for index trickery
     &     bccmil,    !bottom ccm2 interface level
     &     tccmil,    !top ccm2 interface level
     &     numlay,    !num_layer in cloud model
     &     ncclev,    !num_ccm2_level (with overlaps already subtracted)
     &     cbi,       !cloud bottom index
     &     cti,       !cloud top index
     &     ncldlay    !number of layers in cloud itself
c
c     end of zender workspace alterations
c-----------------------------------------------------------------------
c
c compute solar distance factor and cosine solar zenith angle usi
c day value where a round day (such as 213.0) refers to 0z at
c greenwich longitude.
c
c use formulas from paltridge and platt 1976  p. 57, p. 62,63.
c
c compute eccentricity factor (sun-earth distance factor)
c
      theta = 2.*pie*calday/dayspy
      eccf  = 1.000110 + .034221*cos(theta) + .001280*sin(theta) +
     $        .000719*cos(2.*theta) + .000077*sin(2.*theta)
c
c solar declination in radians:
c
      delta = .006918 - .399912*cos(theta) + .070257*sin(theta) -
     $        .006758*cos(2.*theta) + .000907*sin(2.*theta) -
     $        .002697*cos(3.*theta) + .001480*sin(3.*theta)
c
c compute local cosine solar zenith angle, 
c
      sinc = sin(clat)
      sind = sin(delta)
      cosc = cos(clat)
      cosd = cos(delta)
c
c calday is the calender day for greenwich, including fraction
c of day; the fraction of the day represents a local time at
c greenwich; to adjust this to produce a true instantaneous time
c for other longitudes, we must correct for the local time change:
c
      do 10 i=1,plon
c
         phi       = calday + (real(i-1)/real(plon))
         coszrs(i) = sinc*sind - cosc*cosd*cos(2.*pie*phi)
c
   10 continue
c
c convert pressure from pascals to dynes/cm2
c
      do 30 k=1,plev
         do 20 i=1,plon
            pmidrd(i,k) = pmid(i,k)*10.0
            pintrd(i,k) = pint(i,k)*10.0
   20    continue
   30 continue
      do 35 i=1,plon
         pintrd(i,plevp) = pint(i,plevp)*10.0
   35 continue
c
c compute path quantities used in the longwave radiation:
c
      vmmr  = amco2 / amd
      cpwpl = vmmr * 0.5 / (gravit * p0)
      do 40 i=1,plon
         plh2o(i,1)  = rgsslp*h2ommr(i,1)*pintrd(i,1)*pintrd(i,1)
cz the top ccm2 layer absorption should probablly never be zeroed.
         plco2(i,1)  = co2vmr*cpwpl*pintrd(i,1)*pintrd(i,1)
         tclrsf(i,1) = 1.
  40  continue
      do 50 k=1,plev
         do 55 i=1,plon
            plh2o(i,k+1)  = plh2o(i,k) + rgsslp*
     $             (pintrd(i,k+1)**2 - pintrd(i,k)**2)*h2ommr(i,k)
            plco2(i,k+1)  = co2vmr*cpwpl*pintrd(i,k+1)**2
cz alert
cz check whether this layer is with the cloud, if so then zero
cz the trace gas absorption
cz since plco2 is an interface quantity, safer to zero within one
cz index of top and bottom
cz            kcz=plev-k-bccmil+1
cz            if(kcz.ge.(cbi-1).and.kcz.le.(cti+1))then 
cz               plco2(i,k+1)  = 1.e-12*co2vmr*cpwpl*pintrd(i,k+1)**2
cz            else
cz               plco2(i,k+1)  = co2vmr*cpwpl*pintrd(i,k+1)**2
cz            endif
            tclrsf(i,k+1) = tclrsf(i,k)*(1.-cld(i,k+1))
  55     continue
  50  continue
c
c compute ozone mixing ratio from path lengths: 
c
c constants for following sums:
c
      cpl   = v0  / (amd * gravit)
      vmmr  = amd / amo
c
      do 60 k=1,plev
       do 65 i=1,plon
c
          o3mmr(i,k) = (plos(i,k+1)-plos(i,k))
     $                /(cpl*vmmr*(pintrd(i,k+1)-pintrd(i,k)))
c
  65   continue
  60  continue
c
c done
c
      return
      end
      subroutine radout(solin  ,sabtp  ,frsa   ,clrst  ,clrss  ,
     $                  qrs    ,firtp  ,frla   ,clrlt  ,clrls  ,
     $                  qrl    ,slwd   ,srfrad )
c-----------------------------------------------------------------------
c
c copy radiation output quantities to model buffer
c
c change units of the radiative fluxes from cgs to mks
c
c compute the total radiative heat flux at the surface for
c the surface temperature computation
c
c-----------------------------------------------------------------------
      implicit none
c-----------------------------------------------------------------------
c
c     !   Basic grid point resolution parameters
c
      integer plev,plon,plat,plonp2,plevp,pcnst,plond,platd,plevd
      integer nxpt,jintmx,i1,j1,j1m
c
      parameter(plev=  133,        ! number of vertical levels
     $          plevp=plev+1,     ! number of vertical interfaces
     $          plon=1,           ! number of longitudes (T42)
     $          plat=1,           ! number of latitudes (T42)
     $          plonp2=plon+2,    ! length necessary for fft
     $          pcnst=1,          ! number of constituents
     $          nxpt=1,           !
     $          jintmx=2,         !
     $          plond=plon + 1 + 2*nxpt, ! slt extended domain longitude
     $          platd=plat + 2*nxpt + 2*jintmx, ! slt extended domain lat.
     $          plevd=plev*(3+pcnst), ! fold plev,pcnst indices into one
     $          i1=1+nxpt,        ! model starting longitude (3-d)
     $          j1=1+nxpt+jintmx, ! model starting latitude (3-d)
     $          j1m=j1-1)         ! model starting offset (3-d)
c
c-----------------------------------------------------------------------
c
c input/output arguments
c
      real solin(plond),     ! instantaneous incident solar
     $     sabtp(plond),     ! total column absorbed solar flux
     $     frsa(plond),      ! surface absorbed solar flux
     $     clrst(plond),     ! clear sky total column abs solar flux
     $     clrss(plond),     ! clear sky surface absorbed solar flux
     $     qrs(plond,plev),  ! solar heating rate
     $     firtp(plond),     ! net up flux top of model (up-dwn flx)
     $     frla(plond),      ! longwave cooling of surface (up-dwn flx)
     $     clrlt(plond),     ! clr sky net up flx top of model (up-dwn flx)
     $     clrls(plond),     ! clr sky lw cooling of srf (up-dwn flx)
     $     qrl(plond,plev),  ! longwave cooling rate
     $     slwd(plond),      ! surface longwave down flux
     $     srfrad(plond)     ! surface radiative heat flux (frsa+swld)
c 
c local variables
c
      integer       i    ! longitude index
c
      real cgsmks     ! conversion factor for fluxes from cgs to mks
c
      data cgsmks / 1.e-3 /
c-----------------------------------------------------------------------
CDIR$ IVDEP
c
c compute total radiative heating flux for the surface, 
c converting units from cgs to mks:
c
      do 10 i=1,plon
         srfrad(i) = (frsa(i) + slwd(i)) * cgsmks
   10 continue
c
c convert units from cgs to mks in solar fluxes:
c
      do 20 i=1,plon
         solin(i) = solin(i) * cgsmks
         sabtp(i) = sabtp(i) * cgsmks
         frsa(i)  = frsa(i)  * cgsmks
         clrst(i) = clrst(i) * cgsmks
         clrss(i) = clrss(i) * cgsmks
   20 continue
c
c convert units from cgs to mks in longwave fluxes:
c
      do 30 i=1,plon
         firtp(i) = firtp(i) * cgsmks
         frla(i)  = frla(i)  * cgsmks
         clrlt(i) = clrlt(i) * cgsmks
         clrls(i) = clrls(i) * cgsmks
   30 continue
c
c done
c
      return
      end
      subroutine radclw(nrow    ,tg      ,plol    ,plos    ,tnm     ,
     $                  qnm     ,pmid    ,pint    ,pmln    ,piln    ,
     $                  plco2   ,plh2o   ,cld     ,tclrsf  ,qrl     ,
     $                  frla    ,firtp   ,clrls   ,clrlt   ,slwd    ,
     &     fuplw,fdlw,fupcbi,fdncti,
     &     lwfdiv,
     &     bccmil,tccmil,numlay,ncclev,
     &     cbi,cti,ncldlay,
     &     idebug)
CFPP$ NOCONCUR R
c-----------------------------------------------------------------------
c
c compute longwave radiation heating rates and boundary fluxes
c
c uses broad band absorptivity/emissivity method to compute clear
c sky; assumes randomly overlapped clouds with variable cloud
c emissivity to include effects of clouds.
c
c computes clear sky absorptivity/emissivity at lower frequency
c (in general) than the model radiation frequency; uses previously
c computed and stored values for efficiency
c
c-----------------------------------------------------------------------
      implicit none
c-----------------------------------------------------------------------
c
c     !   Basic grid point resolution parameters
c
      integer plev,plon,plat,plonp2,plevp,pcnst,plond,platd,plevd
      integer nxpt,jintmx,i1,j1,j1m
c
      parameter(plev=  133,        ! number of vertical levels
     $          plevp=plev+1,     ! number of vertical interfaces
     $          plon=1,           ! number of longitudes (T42)
     $          plat=1,           ! number of latitudes (T42)
     $          plonp2=plon+2,    ! length necessary for fft
     $          pcnst=1,          ! number of constituents
     $          nxpt=1,           !
     $          jintmx=2,         !
     $          plond=plon + 1 + 2*nxpt, ! slt extended domain longitude
     $          platd=plat + 2*nxpt + 2*jintmx, ! slt extended domain lat.
     $          plevd=plev*(3+pcnst), ! fold plev,pcnst indices into one
     $          i1=1+nxpt,        ! model starting longitude (3-d)
     $          j1=1+nxpt+jintmx, ! model starting latitude (3-d)
     $          j1m=j1-1)         ! model starting offset (3-d)
c
c-----------------------------------------------------------------------
      integer plevp2,plevp3,plevp4
      parameter (plevp2=plev+2,plevp3=plev+3,plevp4=plev+4)
c-----------------------------------------------------------------------
C
C LOGICAL UNIT NUMBERS
C
      common/comlun/nsds    ,nsre    ,nsre1   ,nra1    ,
     $              nrb1    ,nprint  ,nout    ,nread   ,
     $              ninit   ,nsdv    ,nozone  ,nsst    ,
     $              nalb    ,nabem   ,lutag(99)
C
      integer nsds,    !  restart dataset unit
     $        nsre,    !  primary regeneration dataset unit
     $        nsre1,   !  secondary regeneration dataset
     $        nra1,    !  a work file
     $        nrb1,    !  b work file
     $        nprint,  !  alternate print unit
     $        nout,    !  print unit
     $        nread,   !  input unit
     $        ninit,   !  initial dataset unit
     $        nsdv,    !  standard deviation dataset
     $        nozone,  !  ozone dataset
     $        nsst,    !  sst dataset
     $        nalb,    !  albedo dataset
     $        nabem    !  absorptivity/emissivity dataset
      logical lutag
c-----------------------------------------------------------------------
C
C model control time variables
C
      integer
     $      nrstrt  ,nstep   ,nstepr  ,nestep  ,nstop   ,
     $      mdbase  ,msbase  ,mdcur   ,mscur   ,
     $      mbdate  ,mbsec   ,mcdate  ,mcsec   ,
     $      nndbas  ,nnsbas  ,nnbdat  ,nnbsec
      real  calday
c
      common /comtim/
     $      calday  ,nrstrt  ,nstep   ,nstepr  ,nestep  ,nstop   ,
     $      mdbase  ,msbase  ,mdcur   ,mscur   ,
     $      mbdate  ,mbsec   ,mcdate  ,mcsec   ,
     $      nndbas  ,nnsbas  ,nnbdat  ,nnbsec
C
c-----------------------------------------------------------------------
c
c radiation constants
c
      real gravit,   ! gravitational acceleration
     $        rga,   ! 1 over gravit
     $      cpair,   ! heat capacity air at constant pressure
     $     epsilo,   ! ratio mmw h2o to mmw air
     $       sslp,   ! standard pressure
     $     stebol,   ! stephan boltzmann constant
     $     rgsslp,   ! 0.5 / (gravit*sslp)
     $     co2vmr,   ! co2 volume mixing ratio
     $      dpfo3,   ! Doppler factor for o3
     $     dpfco2,   ! Doppler factor for co2
     $     dayspy,   ! solar days in one year
     $        pie    ! pie
c
      common/crdcon/gravit,    rga,  cpair,  epsilo,   sslp,
     $              stebol, rgsslp, co2vmr,   dpfo3, dpfco2,
     $              dayspy,    pie
c
c-----------------------------------------------------------------------
c
c radiation control variables
c
c fradsw = .t. iff full shortwave computation
c fradlw = .t. iff full longwave computation
c
c irad = iteration frequency for radiation computation
c iradae = iteration frequency for absorptivity/
c emissivity computation
c
c
      integer iradae,irad,naclw,nacsw,fnlw,fnsw
      logical aeres
c
      common/crdctl/iradae, irad, naclw, nacsw, fnlw, fnsw, aeres
c
c-----------------------------------------------------------------------
C
C     !  TEMPORARY ARRAYS FOR 512 WORD MULTIPLE READS/WRITES
C
      integer nabebf,lngbuf
      real    pabem
c
      parameter (nabebf=(plond*plevp*plevp +
     $                   plond*plev*4 +
     $                   plond*plevp)/512 + 1)
      parameter (lngbuf=512*nabebf)
      common/abemub/pabem(plat)       ! position array for random access
c-----------------------------------------------------------------------
c
c input arguments
c
      integer nrow                ! model latitude index
      real tg(plond),             ! ground (skin) temperature
     $     plol(plond,plevp),     ! o3 pressure wghted path length
     $     plos(plond,plevp)      ! o3 path length
c
c input arguments which are only passed to other routines
c
      real tnm(plond,plev),        ! level temperature
     $     qnm(plond,plev),        ! level moisture field
     $     pmid(plond,plev),       ! level pressure
     $     pint(plond,plevp),      ! model interface pressure
     $     pmln(plond,plev),       ! ln(pmid)
     $     piln(plond,plevp),      ! ln(pint)
     $     plco2(plond,plevp),     ! path length co2
     $     plh2o(plond,plevp)      ! path length h2o
c
c input/Output arguments
c
cz      real cld(plond,plevp),       ! cloud cover
      real cld(plond,plevp),       ! (effective) cloud cover (=cld frac x emis.)
     $      tclrsf(plond,plevp)    ! clear sky fraction
c
c output arguments
c
      real
     $     qrl(plond,plev),        ! longwave heating rate
     $     frla(plond),            ! surface cooling flux
     $     firtp(plond),           ! net outgoing flux
     $     clrls(plond),           ! clear sky surface cooing
     $     clrlt(plond),           ! net clear sky outgoing flux
     $     slwd(plond)             ! down longwave flux at surface
c
c local workspace
c
      integer     i,   ! longitude index
     $            k,   ! level index
     $           k1,   ! level index
     $           k2,   ! level index
     $           k3,   ! level index
     $           km,   ! level index
     $          kmm,   ! level index
     $          km1,   ! level index
     $          km2,   ! level index
     $          km3,   ! level index
     $          km4    ! level index
c
      integer  itop,   ! level index
     $         icld,   ! level index
     $          klo,   ! level index
     $          khi,   ! level index
     $          mp1,   ! level index
     $         mp12,   ! level index
     $          mm1,   ! level index
     $        nptsc    ! level index
c
      real      sum,   ! sum
     $         tmp1,   ! temporary 1
     $     absbt(plond)
c
      real co2em(plond,plevp),   ! layer co2 normalized plnck function drvtv
     $     co2eml(plond,plev),   ! intrfc co2 normalized plnck function drvtv
     $     delt(plond),          ! diff t**4 mid layer to top interface
     $     delt1(plond),         ! diff t**4 lower intrfc to mid layer
     $     bk1(plond),           ! absrptvty for vertical quadrature
     $     bk2(plond),           ! absrptvty for vertical quadrature
     $     ful(plond,plevp),          ! total upwards longwave flux
     $     fsul(plond,plevp),         ! clear sky upwards longwave flux
     $     fdl(plond,plevp),          ! total downwards longwave flux
     $     fsdl(plond,plevp),         ! clear sky downwards longwave flux
     $     fclb4(plond,plev),         ! sig t**4 for cloud bottom interface
     $     fclt4(plond,plev),         ! sig t**4 for cloud top interface
     $     s(plond,plevp,plevp)       ! flx integral sum
      real absnxt(plond,plev,4),      ! nearest layer absorptivities
     $     abstot(plond,plevp,plevp), ! non-adjacent layer absorptivites
     $     emstot(plond,plevp),       ! total emissivity
     $     tplnka(plond,plevp),       ! planck fnctn tmp
     $     s2c(plond,plevp),          ! h2o cont amount
     $     s2t(plond,plevp),          ! h2o cont tmp
     $     w(plond,plevp),            ! h2o path
     $     tplnke(plond)              ! planck fnctn tmp
      real h2otr(plond,plevp),        ! h2o trnmsn for o3 overlap
     $     co2t(plond,plevp),         ! prs wghted tmp path
     $     tint(plond,plevp),         ! interface tmp
     $     tint4(plond,plevp),        ! interface tmp**4
     $     tlayr(plond,plevp),        ! level tmp
     $     tlayr4(plond,plevp)        ! level tmp**4
      integer ipos(plond),            ! array for specified condition
     $        indxc(plond),           ! indices for cloud covered points
     $        klov(plond),            ! cloud lowest level index
     $        khiv(plond)             ! cloud highest level index
      real absems(lngbuf)
      equivalence (abstot,absems(1                                   )),
     $            (absnxt,absems(1 + plond*plevp*plevp               )),
     $            (emstot,absems(1 + plond*plevp*plevp + plond*plev*4))
c
      external radtpl,    ! compute path lengths
     $         radems,    ! h2o,co2,o3 emissivity
     $         radabs,    ! h2o,co2,o3 absorptivity
     $       writeric,    ! write for abs/ems
     $        readric     ! read  for abs/ems
c
c------------------------------------------------------------------------
c     zender cirrus cloud hack workspace
c
      real cgsmks     ! conversion factor for fluxes from cgs to mks
c
      data cgsmks / 1.e-3 /
c
      integer 
     &     kcz,       !index for index trickery
     &     idebug,    !debugging level from cloud model
     &     bccmil,    !bottom ccm2 interface level
     &     tccmil,    !top ccm2 interface level
     &     numlay,    !num_layer in cloud model
     &     ncclev,    !num_ccm2_level (with overlaps already subtracted)
     &     cbi,       !cloud bottom index
     &     cti,       !cloud top index
     &     ncldlay    !number of layers in cloud itself
C
      real
     &     lwfdiv(0:numlay+1), !LW flux divergence
     &     fuplw(0:numlay+1), !total longwave upward directed flux
     &     fdlw(0:numlay+1), !total longwave downward directed flux
     &     fupcbi, !longwave upward directed flux at cloud base interface
     &     fdncti !longwave downward directed flux at cloud top interface
c
c     end of zender workspace alterations
c-----------------------------------------------------------------------
c
c  initialize and recompute the tclrsf array
c
      do 40 k=1,plev
         do 30 i=1,plon
            fclb4(i,k) = 0.
            fclt4(i,k) = 0.
            tclrsf(i,k+1) = tclrsf(i,k)*(1. - cld(i,k+1))
   30    continue
   40 continue
c
c calculate some temperatures needed to derive absorptivity and 
c emissivity, as well as some h2o path lengths
c
      call radtpl(tnm    ,tg     ,qnm    ,pmid   ,pint   ,
     $            plh2o  ,tplnka ,s2c    ,s2t    ,w      ,
     $            tplnke ,tint   ,tint4  ,tlayr  ,tlayr4,
     $            pmln   ,piln)
c
c do emissivity and absorptivity calculations
c only if abs/ems computation
c
      if(nstep.eq.0 .or.
     $   (mod(nstep-1,iradae).eq.0 .and. nstep.ne.1)) then
c
c compute total emissivity:
c
         call radems(s2c    ,s2t    ,w     ,tplnke ,plh2o ,
     $               pint   ,plco2  ,tint  ,tint4  ,tlayr ,
     $               tlayr4 ,plol   ,plos  ,co2em  ,co2eml,
     $               co2t   ,h2otr  ,emstot)
c
c compute total absorptivity:
c
         call radabs(pmid   ,pint   ,co2em ,co2eml ,tplnka,
     $               s2c    ,s2t    ,w     ,h2otr  ,plco2 ,
     $               plh2o  ,co2t   ,tint  ,tlayr  ,plol  ,
     $               plos   ,pmln   ,piln  ,abstot ,absnxt)
c
CMIC$ GUARD 3
        call writeric(nabem,absems(1),lngbuf,nrow)
CMIC$ END GUARD 3
      else       
c
c retrieve total absorptivity and emissivity from
c last abs/ems computation
c
CMIC$ GUARD 3
        call readric(nabem,absems(1),lngbuf,nrow)
CMIC$ END GUARD 3
      endif
c
c compute fluxes and cooling rates.
c Initialize longitude index subset.
c
      do 100 i=1,plon
         ipos(i) = 1
  100 continue
c
c find the lowest and highest level cloud for each grid point
c
      do 120 i=1,plon
         klov(i) = 0
         do 110 k=1,plev
            if(cld(i,plevp2-k) .gt. 0.0)then
               klov(i) = k
               go to 120
            endif
  110    continue
  120 continue
c
      do 140 i=1,plon
         khiv(i) = klov(i)
         itop    = klov(i)
         if(itop.eq.0) itop = 1
         do 130 k=plev,itop,-1
            if(cld(i,plevp2-k) .gt. 0.0)then
               khiv(i) = k
               go to 140
            endif
  130    continue
  140 continue
c
      do 160 i=1,plon
         if(klov(i).ne.0)then
            do 150 k=klov(i),khiv(i)
               fclt4(i,plevp-k) = stebol*tint4(i,plevp2-k)
               fclb4(i,plevp-k) = stebol*tint4(i,plevp3-k)
  150       continue
         endif
  160 continue
c
c compute sums used in integrals (all longitude points)
c
c definition of bk1 & bk2 depends on finite differencing.
c for trapezoidal rule bk1=bk2. trapezoidal rule applied for
c nonadjacent layers only.
c
c delt=t**4 in layer above current sigma level km.
c delt1=t**4 in layer below current sigma level km.
c
      do 210 i=1,plon
         delt(i) = tint4(i,plev) - tlayr4(i,plevp)
         delt1(i) = tlayr4(i,plevp) - tint4(i,plevp)
         s(i,plevp,plevp) = stebol*(delt1(i)*absnxt(i,plev,1) +
     $                              delt (i)*absnxt(i,plev,4))
         s(i,plev,plevp)  = stebol*(delt (i)*absnxt(i,plev,2) +
     $                              delt1(i)*absnxt(i,plev,3))
  210 continue
      do 230 k=1,plev-1
         do 220 i=1,plon
            bk2(i) = (abstot(i,k,plev) + abstot(i,k,plevp))*0.5
            bk1(i) = bk2(i)
            s(i,k,plevp) = stebol*(bk2(i)*delt(i) + bk1(i)*delt1(i))
  220    continue
  230 continue
c
c all k, km>1
c
      do 300 km=plev,2,-1
         do 240 i=1,plon
            delt(i)  = tint4(i,km-1) - tlayr4(i,km)
            delt1(i) = tlayr4(i,km) - tint4(i,km)
  240    continue
         do 290 k=plevp,1,-1
            if (k.eq.km) then
               do 250 i=1,plon
                  bk2(i) = absnxt(i,km-1,4)
                  bk1(i) = absnxt(i,km-1,1)
  250          continue
            else if(k.eq.km-1) then
               do 260 i=1,plon
                  bk2(i) = absnxt(i,km-1,2)
                  bk1(i) = absnxt(i,km-1,3)
  260          continue
            else
               do 270 i=1,plon
                  bk2(i) = (abstot(i,k,km-1) + abstot(i,k,km))*0.5
                  bk1(i) = bk2(i)
  270          continue
            endif
            do 280 i=1,plon
               s(i,k,km) = s(i,k,km+1) + stebol*
     $                    (bk2(i)*delt(i) + bk1(i)*delt1(i))
  280       continue
  290    continue
  300 continue
c
c computation of clear sky fluxes
c always set first level of fsul
c
      do 340 i=1,plon
         fsul(i,plevp) = stebol*(tg(i)**4)
  340 continue
c
c all longitude points
c
      do 360 k=1,plev
         do 350 i=1,plon
            tmp1 = fsul(i,plevp) - stebol*tint4(i,plevp)
            fsul(i,k) = fsul(i,plevp) - abstot(i,k,plevp)*tmp1 +
     $                  s(i,k,k+1)
  350    continue
  360 continue
c
c downward clear sky fluxes
c store intermediate quantities in down flux
c
      do 380 k=1,plevp
         do 370 i=1,plon
            fsdl(i,k) = stebol*(tplnke(i)**4)*emstot(i,k)
  370    continue
  380 continue
c
c store the downward emission from level 1
c = total gas emission * sigma t**4.
c fsdl does not yet include all terms
c
      do 390 i=1,plon
         absbt(i) = fsdl(i,plevp)
  390 continue
c
c fsdl(i,plevp) assumes isothermal layer
c
      do 420 k=2,plev
         do 410 i=1,plon
            fsdl(i,k) = fsdl(i,k) - (s(i,k,2) - s(i,k,k+1))
  410    continue
  420 continue
      do 430 i=1,plon
         fsdl(i,plevp) = absbt(i) - s(i,plevp,2)
  430 continue
c
c computation of the upward fluxes
c
c first set fluxes to clear sky values, then modify for clouds
c
      do 450 k=1,plevp
         do 440 i=1,plon
            ful(i,k) = fsul(i,k)
            fdl(i,k) = fsdl(i,k)
  440    continue
  450 continue
c
c modifications for clouds
c
c further qualify longitude subset for computations by changing
c ipos from 1 to 2 where there are clouds
c (total cloud fraction <= 1.e-3 treated as clear)
c
      do 460 i=1,plon
         if((1.-tclrsf(i,plevp)) .gt. 1.e-3) then
            ipos(i) = 2
         endif
  460 continue
c
c generate longitude index list (indxc) for cloud modifications
c
      call wheneq(plon,ipos,1,2,indxc,nptsc)
c
c compute downflux at level 1 for cloudy sky
c
      do 470 icld=1,nptsc
         i = indxc(icld)
c
c first clear sky flux plus flux from cloud at level 1
c
         fdl(i,plevp) = fsdl(i,plevp)*tclrsf(i,plev)/
     $     tclrsf(i,plevp-khiv(i)) + fclb4(i,plev-1)*cld(i,plev)
  470 continue
c
c flux emitted by other layers
c
      do 490 icld=1,nptsc
         i = indxc(icld)
         do 480 km=3,khiv(i)
            km1 = plevp - km
            km2 = plevp2 - km
            km4 = plevp4 - km
            tmp1 = cld(i,km2)*tclrsf(i,plev)/tclrsf(i,km2)
            fdl(i,plevp) = fdl(i,plevp) +
     $           (fclb4(i,km1) - s(i,plevp,km4))*tmp1
  480    continue
  490 continue
c
c begin outer longitude loop for cloud modifications
c (operate on points with clouds only)
c
      do 560 icld=1,nptsc
         i = indxc(icld)
         klo = klov(i)
         khi = khiv(i)
         mp1 = khi + 1
         mp12 = plevp2 - mp1
         mm1 = khi - 1
         do 510 k=klo,mm1
            k1 = plevp - k
            k2 = plevp2 - k
            k3 = plevp3 - k
            sum = fsul(i,k2)*(tclrsf(i,plevp)/tclrsf(i,k1))
            do 500 km=klo,k
               km1 = plevp - km
               km2 = plevp2 - km
               km3 = plevp3 - km
               sum = sum + (fclt4(i,km1) + s(i,k2,k3) - s(i,k2,km3))*
     $                     cld(i,km2)*(tclrsf(i,km1)/tclrsf(i,k1))
  500       continue
            ful(i,k2) = sum
  510    continue
         do 530 k=khi,plevp
            k2 = plevp2 - k
            k3 = plevp3 - k
            sum = fsul(i,k2)*tclrsf(i,plevp)/tclrsf(i,mp12)
            do 520 km=klo,khi
               km1 = plevp - km
               km2 = plevp2 - km
               km3 = plevp3 - km
               sum = sum + (cld(i,km2)*tclrsf(i,km1)/tclrsf(i,mp12))*
     $                     (fclt4(i,km1) + s(i,k2,k3) - s(i,k2,km3))
  520       continue
            ful(i,k2) = sum
  530    continue
c
c computation of the downward fluxes
c
         do 550 k=2,mm1
            k1 = plevp - k
            k2 = plevp2 - k
            k3 = plevp3 - k
            sum = 0.
            kmm = max0(k+1,klo)
            do 540 km=kmm,khi
               km1 = plevp - km
               km2 = plevp2 - km
               km4 = plevp4 - km
               sum = sum + (cld(i,km2)*tclrsf(i,k1)/tclrsf(i,km2))*
     $                     (fclb4(i,km1) - s(i,k2,km4) + s(i,k2,k3))
  540       continue
            fdl(i,k2) = sum + fsdl(i,k2)*(tclrsf(i,k1)/tclrsf(i,mp12))
  550    continue
c
c end cloud modification longitude loop
c
  560 continue
c
c back to original longitude subset.
c put down longwave flux in local array.
c
      do 570 i=1,plon
         slwd(i) = fdl(i,plevp)
         frla(i) = ful(i,plevp) - fdl(i,plevp)
  570 continue
c
c all longitudes: store history tape quantities
c
      do 580 i=1,plon
c
c net flux
c
         frla(i) = ful(i,plevp) - fdl(i,plevp)
c
c clear sky flux at top of atmosphere
c
         clrlt(i) = fsul(i,1)
         clrls(i) = fsul(i,plevp) - fsdl(i,plevp)
c
c outgoing ir
c
         firtp(i) = ful(i,1) - fdl(i,1)
  580 continue
c
c computation of longwave heating (k per sec)
c
      do 600 k=1,plev
         do 590 i=1,plon
            qrl(i,k) = (ful(i,k) - fdl(i,k) - ful(i,k+1) + fdl(i,k+1))*
     $                 gravit/((pint(i,k) - pint(i,k+1))*cpair)
  590    continue
  600 continue
cz    here's where the mean intensity field (at layer center) is computed
cz    and the units are converted MKS before returning to cloud model
cz    altitude loop runs down from top of cloud (so indices increase in F77)
      do 686 i=1,plon
         do 687 k=plev+1-bccmil-numlay,plev+1-bccmil
            kcz=plev-k-bccmil+1
            fuplw(kcz) = cgsmks*ful(i,k)
            fdlw(kcz) = cgsmks*fdl(i,k)
cz    the following lines would return the mid-layer fluxes
cz         do 687 k=plev+1-bccmil-numlay,plev+1-bccmil-1
cz            meanjir(kcz) = .5*cgsmks*
cz     &           (ful(i,k)+fdl(i,k)+
cz     &           ful(i,k+1)+fdl(i,k+1))/
cz     &           (2.*pie)
cz            fuplw(kcz) = .5*cgsmks*
cz     &           (ful(i,k)+ful(i,k+1))
cz            fdlw(kcz) = .5*cgsmks*
cz     &           (fdl(i,k)+fdl(i,k+1))
cz    Note this next step defines lwfdiv(0) which is meaningless since
cz    lwfdiv is a layer-midpoint-indexed array.
            lwfdiv(kcz) = 
     &           cgsmks*(ful(i,k) - fdl(i,k) - ful(i,k+1) + fdl(i,k+1))
 687     continue
         fupcbi=cgsmks*ful(i,plev+1-bccmil-1-cbi+2)
         fdncti=cgsmks*fdl(i,plev+1-bccmil-1-cti+1)
 686  continue
c     
c     done
c     
      return
      end
      subroutine radtpl(tnm    ,tg     ,qnm    ,pbr    ,pnm    ,
     $                  plh2o  ,tplnka ,s2c    ,s2t    ,w      ,
     $                  tplnke ,tint   ,tint4  ,tlayr  ,tlayr4 ,
     $                  pmln   ,piln   )
c-----------------------------------------------------------------------
c
c compute temperatures and path lengths for longwave radiation 
c
c-----------------------------------------------------------------------
      implicit none
c-----------------------------------------------------------------------
c
c     !   Basic grid point resolution parameters
c
      integer plev,plon,plat,plonp2,plevp,pcnst,plond,platd,plevd
      integer nxpt,jintmx,i1,j1,j1m
c
      parameter(plev=  133,        ! number of vertical levels
     $          plevp=plev+1,     ! number of vertical interfaces
     $          plon=1,           ! number of longitudes (T42)
     $          plat=1,           ! number of latitudes (T42)
     $          plonp2=plon+2,    ! length necessary for fft
     $          pcnst=1,          ! number of constituents
     $          nxpt=1,           !
     $          jintmx=2,         !
     $          plond=plon + 1 + 2*nxpt, ! slt extended domain longitude
     $          platd=plat + 2*nxpt + 2*jintmx, ! slt extended domain lat.
     $          plevd=plev*(3+pcnst), ! fold plev,pcnst indices into one
     $          i1=1+nxpt,        ! model starting longitude (3-d)
     $          j1=1+nxpt+jintmx, ! model starting latitude (3-d)
     $          j1m=j1-1)         ! model starting offset (3-d)
c
c ----------------------------------------------------------------------
c
c radiation constants
c
      real gravit,   ! gravitational acceleration
     $        rga,   ! 1 over gravit
     $      cpair,   ! heat capacity air at constant pressure
     $     epsilo,   ! ratio mmw h2o to mmw air
     $       sslp,   ! standard pressure
     $     stebol,   ! stephan boltzmann constant
     $     rgsslp,   ! 0.5 / (gravit*sslp)
     $     co2vmr,   ! co2 volume mixing ratio
     $      dpfo3,   ! Doppler factor for o3
     $     dpfco2,   ! Doppler factor for co2
     $     dayspy,   ! solar days in one year
     $        pie    ! pie
c
      common/crdcon/gravit,    rga,  cpair,  epsilo,   sslp,
     $              stebol, rgsslp, co2vmr,   dpfo3, dpfco2,
     $              dayspy,    pie
c
c ----------------------------------------------------------------------
c
c input arguments
c
      real tnm(plond,plev),     ! model level temperatures
     $     tg(plond),           ! surface skin temperature
     $     qnm(plond,plev),     ! model level specific humidity
     $     pbr(plond,plev),     ! prsr at model mid-levels (dynes/cm2)
     $     pnm(plond,plevp),    ! prsr at model interfaces (dynes/cm2)
     $     plh2o(plond,plevp)   ! prs wghtd h2o path
c
c output arguments
c
      real tplnka(plond,plevp), ! level tmp from interface tmps
     $     s2c(plond,plevp),    ! h2o continuum path length
     $     s2t(plond,plevp),    ! h2o tmp and prs wghtd path length
     $     w(plond,plevp),      ! h2o prs wghtd path length
     $     tplnke(plond),       ! equal to tplnka
     $     tint(plond,plevp),   ! layer interface temperature
     $     tint4(plond,plevp),  ! tint to the 4th power
     $     tlayr(plond,plevp),  ! k-1 level temperature
     $     tlayr4(plond,plevp), ! tlayr to the 4th power
     $     pmln(plond,plev),    ! ln(pmidm1)
     $     piln(plond,plevp)    ! ln(pintm1)
c
c local variables
c
      integer   i,   ! longitude index
     $          k    ! level index

      real r296,     ! inverse standard temperature for h2o continuum
     $     repsil,   ! inverse ratio mol weight h2o to dry air
     $         dy,   ! thickness of layer for tmp interpolation
     $       dpnm,   ! pressure thickness of layer
     $     dpnmsq,   ! prs squared difference across layer
     $       rtnm    ! inverse level temperature
c
      r296   = 1./296.
      repsil = 1./epsilo
c
c set the top and bottom intermediate level temperatures,
c top level planck temperature and top layer temp**4.
c
c tint is lower interface temperature
c (not available for bottom layer, so use ground temperature)
c
      do 30 i=1,plon
         tint(i,plevp)  = tg(i)
         tint4(i,plevp) = tint(i,plevp)**4
         tplnka(i,1)    = tnm(i,1)
         tint(i,1)      = tplnka(i,1)
         tlayr4(i,1)    = tplnka(i,1)**4
         tint4(i,1)     = tlayr4(i,1)
   30 continue
c
c intermediate level temperatures are calculated based on the tem
c at the full level below less dy*delta t (between the full level
c
      do 50 k=2,plev
         do 40 i=1,plon
            dy         = (piln(i,k)-pmln(i,k)) / 
     $                   (pmln(i,k-1)-pmln(i,k))
            tint(i,k)  = tnm(i,k) - dy*(tnm(i,k)-tnm(i,k-1))
            tint4(i,k) = tint(i,k)**4
   40    continue
   50 continue
c
c now set the layer temp=full level temperatures and establish a
c planck temperature for absorption (tplnka) which is the average
c the intermediate level temperatures.  note that tplnka is not
c equal to the full level temperatures.
c
      do 70 k=2,plevp
         do 60 i=1,plon
            tlayr(i,k)  = tnm(i,k-1)
            tlayr4(i,k) = tlayr(i,k)**4
            tplnka(i,k) = .5*(tint(i,k) + tint(i,k-1))
   60    continue
   70 continue
c
c calculate tplank for emissivity calculation.
c assume isothermal tplnke i.e. all levels=ttop.
c
      do 80 i=1,plon
         tplnke(i)  = tplnka(i,1)
         tlayr(i,1) = tint(i,1)
   80 continue
c
c now compute h2o path fields:
c
      do 90 i=1,plon
         s2t(i,1) = plh2o(i,1) * tnm(i,1)
         w(i,1)   = (plh2o(i,1)*2.) / pnm(i,1)
         s2c(i,1) = plh2o(i,1) * qnm(i,1) * repsil
   90 continue
      do 110 k=1,plev
         do 100 i=1,plon
            dpnm       = pnm(i,k+1) - pnm(i,k)
            dpnmsq     = pnm(i,k+1)**2 - pnm(i,k)**2
            rtnm       = 1./tnm(i,k)
            s2t(i,k+1) = s2t(i,k) + rgsslp*dpnmsq*qnm(i,k)*tnm(i,k)
            w(i,k+1)   = w(i,k)   + rga*qnm(i,k)*dpnm
            s2c(i,k+1) = s2c(i,k) + rgsslp*dpnmsq*qnm(i,k)*
     $             exp(1800.*(rtnm - r296))*qnm(i,k)*repsil
  100    continue
  110 continue
c
c done
c
      return
      end
      subroutine radabs(pbr    ,pnm    ,co2em ,co2eml ,tplnka,
     $                  s2c    ,s2t    ,w     ,h2otr  ,plco2 ,
     $                  plh2o  ,co2t   ,tint  ,tlayr  ,plol  ,
     $                  plos   ,pmln   ,piln  ,abstot ,absnxt)
c-----------------------------------------------------------------------
c
c compute absorptivities for h2o, co2, and o3
c
c
c h2o  ....  uses nonisothermal emissivity for water vapor from
c            Ramanathan, V. and  P.Downey, 1986: A Nonisothermal
c            Emissivity and Absorptivity Formulation for Water Vapor
c            Jouranl of Geophysical Research, vol. 91., D8, pp 8649-8666
c
c
c co2  ....  uses absorptance parameterization of the 15 micro-meter
c            (500 - 800 cm-1) band system of Carbon Dioxide, from
c            Kiehl, J.T. and B.P.Briegleb, 1991: A New Parameterization
c            of the Absorptance Due to the 15 micro-meter Band System
c            of Carbon Dioxide Jouranl of Geophysical Research,
c            vol. 96., D5, pp 9013-9019
c
c o3   ....  uses absorptance parameterization of the 9.6 micro-meter
c            band system of ozone, from Ramanathan, V. and R.E.Dickinson,
c            1979: The Role of stratospheric ozone in the zonal and
c            seasonal radiative energy balance of the earth-troposphere
c            system. Journal of the Atmospheric Sciences, Vol. 36,
c            pp 1084-1104
c
c
c computes individual absorptivities for non-adjacent layers, accounting 
c for band overlap, and sums to obtain the total; then, computes the
c nearest layer contribution.
c
c
c-----------------------------------------------------------------------
      implicit none
c-----------------------------------------------------------------------
c
c     !   Basic grid point resolution parameters
c
      integer plev,plon,plat,plonp2,plevp,pcnst,plond,platd,plevd
      integer nxpt,jintmx,i1,j1,j1m
c
      parameter(plev=  133,        ! number of vertical levels
     $          plevp=plev+1,     ! number of vertical interfaces
     $          plon=1,           ! number of longitudes (T42)
     $          plat=1,           ! number of latitudes (T42)
     $          plonp2=plon+2,    ! length necessary for fft
     $          pcnst=1,          ! number of constituents
     $          nxpt=1,           !
     $          jintmx=2,         !
     $          plond=plon + 1 + 2*nxpt, ! slt extended domain longitude
     $          platd=plat + 2*nxpt + 2*jintmx, ! slt extended domain lat.
     $          plevd=plev*(3+pcnst), ! fold plev,pcnst indices into one
     $          i1=1+nxpt,        ! model starting longitude (3-d)
     $          j1=1+nxpt+jintmx, ! model starting latitude (3-d)
     $          j1m=j1-1)         ! model starting offset (3-d)
c
c-----------------------------------------------------------------------
c
c water vapor narrow band constants for lw computations
c
      real realk,st,a1,a2,b1,b2,
     $     coefa,coefb,coefc,coefd,
     $     coefe,coeff,coefg,coefh,
     $     coefi,coefj,coefk,
     $     c1,c2,c3,c4,c5,c6,c7,
     $     c8 ,c9 ,c10,c11,c12,c13,c14,c15,c16,c17,
     $     c18,c19,c20,c21,c22,c23,c24,c25,c26,c27,
     $     c28,c29,c30,c31,
     $     fwcoef,fwc1,fwc2,fc1,cfa1
c
      common/crdcae/realk(2), st(2), a1(2), a2(2), b1(2), b2(2),
c
c constant coefficients for water vapor absorptivity and emissivi
c
     $              coefa(3,4),coefb(4,4),coefc(3,4),coefd(4,4),
     $              coefe(3,4),coeff(6,2),coefg(2,4),coefh(2,4),
     $              coefi(6,2),coefj(3,2),coefk(3,2),
     $              c1(4),c2(4),c3(4),c4(4),c5(4),c6(4),c7(4),
     $              c8 ,c9 ,c10,c11,c12,c13,c14,c15,c16,c17,
     $              c18,c19,c20,c21,c22,c23,c24,c25,c26,c27,
     $              c28,c29,c30,c31,
c
c farwing correction constants for narrow-band emissivity model
c introduce farwing correction to account for the
c deficiencies in narrow-band model used to derive the
c emissivity. tuned with arkings line calculations.
c
     $              fwcoef,fwc1,fwc2,fc1,cfa1
c
c-----------------------------------------------------------------------
c
c radiation constants
c
      real gravit,   ! gravitational acceleration
     $        rga,   ! 1 over gravit
     $      cpair,   ! heat capacity air at constant pressure
     $     epsilo,   ! ratio mmw h2o to mmw air
     $       sslp,   ! standard pressure
     $     stebol,   ! stephan boltzmann constant
     $     rgsslp,   ! 0.5 / (gravit*sslp)
     $     co2vmr,   ! co2 volume mixing ratio
     $      dpfo3,   ! Doppler factor for o3
     $     dpfco2,   ! Doppler factor for co2
     $     dayspy,   ! solar days in one year
     $        pie    ! pie
c
      common/crdcon/gravit,    rga,  cpair,  epsilo,   sslp,
     $              stebol, rgsslp, co2vmr,   dpfo3, dpfco2,
     $              dayspy,    pie
c
c-----------------------------------------------------------------------
c
c input arguments
c
      real pbr(plond,plev),           ! prssr at mid-levels (dynes/cm2)
     $     pnm(plond,plevp),          ! prssr at interfaces (dynes/cm2)
     $     co2em(plond,plevp),        ! co2 emissivity function
     $     co2eml(plond,plev),        ! co2 emissivity function
     $     tplnka(plond,plevp),       ! planck fnctn level temperature
     $     s2c(plond,plevp),          ! h2o continuum path length
     $     s2t(plond,plevp),          ! h2o tmp and prs wghted path
     $     w(plond,plevp),            ! h2o prs wghted path
     $     h2otr(plond,plevp),        ! h2o trnsmssn fnctn for o3 overlap
     $     plco2(plond,plevp),        ! co2 prs wghted path length
     $     plh2o(plond,plevp),        ! h2o prs wfhted path length
     $     co2t(plond,plevp),         ! tmp and prs wghted path length
     $     tint(plond,plevp),         ! interface temperatures
     $     tlayr(plond,plevp),        ! k-1 level temperatures
     $     plol(plond,plevp),         ! ozone prs wghted path length
     $     plos(plond,plevp)          ! ozone path length
      real pmln(plond,plev),          ! ln(pmidm1)
     $     piln(plond,plevp)          ! ln(pintm1)
c
c output arguments
c
      real abstot(plond,plevp,plevp), ! total absorptivity
     $     absnxt(plond,plev,4)       ! total nearest layer absorptivity
c
c local workspace
c
      integer       i,   ! logitude index
     $              k,   ! level index     
     $             k1,   ! level index     
     $             k2,   ! level index     
     $             kn,   ! nearest level index
     $          iband    ! band  index
c
      real pnew(plond),     !
     $     trline(plond,2), !
     $     u(plond),        !
     $     tbar(plond,4),   !
     $     emm(plond,4),    !
     $     o3emm(plond,4),  !
     $     o3bndi,          ! ozone band parameter
     $     temh2o(plond,4), !
     $     k21,             ! exponential coefficient used to calculate
c                           ! rotation band transmissivity in the 650-80
c                           ! cm-1 region (tr1)
     $     k22,             ! exponential coefficient used to calculate
c                           ! rotation band transmissivity in the 500-65
c                           !  cm-1 region (tr2)
     $     uc1(plond)       !
      real to3h2o(plond),   ! h2o trnsmsn for overlap with o3
     $     pi,              ! for co2 absorptivity computation 
     $     sqti(plond),     ! "
     $     et,              ! " 
     $     et2,             ! " 
     $     et4,             ! " 
     $     omet,            ! " 
     $     f1co2,           ! " 
     $     f2co2(plond),    ! " 
     $     f3co2(plond),    ! " 
     $     t1co2(plond),    ! " 
     $     sqwp,            ! " 
     $     f1sqwp(plond)    ! " 
      real oneme,           ! " 
     $     alphat,          ! " 
     $     wco2,            ! " 
     $     posqt,           ! " 
     $     u7,              ! " 
     $     u8,              ! " 
     $     u9,              ! " 
     $     u13,             ! " 
     $     rbeta7,          ! " 
     $     rbeta8,          ! " 
     $     rbeta9,          ! " 
     $     rbeta13          ! " 
      real tpatha(plond),     ! for absorptivity computation
     $     a,                 ! "
     $     abso(plond,6),     ! absorptivity for various gases/bands
     $     dtp(plond),        ! "
     $     dtx(plond),        ! "
     $     dty(plond),        ! "
     $     dtz(plond),        ! "
     $     term1(plond,4),    ! "
     $     term2(plond,4)     ! "
      real term3(plond,4),    ! "
     $     term4(plond,4),    ! "
     $     term5(plond,4),    ! "
     $     term6(plond,plevp),! "
     $     term7(plond,2),    ! "
     $     term8(plond,2),    ! "
     $     term9(plond,plevp),! "
     $     tr1,               ! "
     $     tr10(plond),       ! "
     $     tr2                ! "
      real tr5,               ! "
     $     tr6,               ! "
     $     tr9(plond),        ! "
     $     uc(plond)          ! "
      real sqrtu(plond),      ! " 
     $     fwk(plond),        ! " 
     $     fwku(plond),       ! " 
     $     r2st(2),           ! " 
     $     dtyp15(plond),     ! " 
     $     dtyp15sq(plond),   ! " 
     $     to3co2(plond),     ! " 
     $     dpnm(plond),       ! " 
     $     pnmsq(plond,plevp),! " 
     $     dw(plond),         ! " 
     $     uinpl(plond,4),    ! " 
     $     winpl(plond,4),    ! " 
     $     zinpl(plond,4),    ! " 
     $     pinpl(plond,4),    ! " 
     $     dplh2o(plond)      ! " 
      real            r80257, ! " 
     $                  r293, ! " 
     $                  r250, ! " 
     $                 r3205, ! " 
     $                  r300, ! " 
     $                 rsslp, ! " 
     $                r2sslp  ! " 
c
      real  ds2c,     ! used in absorptivity computation
     $      a11,      ! "
     $      a31,      ! "
     $      a21,      ! "
     $      a22,      ! "
     $      a23,      ! "
     $      t1t4,     ! "
     $      t2t5,     ! "
     $      rsum,     ! "
     $      a41,      ! "
     $      a51,      ! "
     $      a61       ! "
c
      real  phi,      ! "
     $      psi,      ! "
     $      cf812,    ! "
     $      ubar,     ! "
     $      pbar,     ! "
     $      g4        ! "
c
      real  dplos,    ! used in absorptivity computation
     $      dplol,    ! "
     $      tlocal,   ! "
     $      beta,     ! "
     $      rphat,    ! "
     $      tcrfac,   ! "
     $      tmp1,     ! "
     $      u1,       ! "
     $      realnu,   ! "
     $      tmp2,     ! "
     $      u2,       ! "
     $      rsqti     ! "
c
      real  tpath,    ! "
     $      tmp3,     ! "
     $      rdpnmsq,  ! "
     $      rdpnm,    ! "
     $      p1,       ! "
     $      p2,       ! "
     $      dtym10,   ! "
     $      dplco2,   ! "
     $      corfac,   ! "
     $      g2,       ! "
     $      te,       ! "
     $      denom     ! "
c
c transmission terms for various spectral intervals:
c
      real    trab1(plond),  ! h2o     0 -  800 cm-1
     $        trab2(plond),  ! h2o   500 -  800 cm-1
     $        trab3(plond),  ! co2   500 -  800 cm-1
     $        trab4(plond),  ! h2o   800 - 1000 cm-1
     $        trab5(plond),  ! o3     9.6 micro-meter band
     $        trab6(plond),  ! h2o  1000 - 1200 cm-1
     $        trab7(plond)   ! h2o  1200 - 2200 cm-1
c
      real bndfct, ! band absorptance parameter for co2
     $     absbnd  ! proportional to co2 band absorptance
c
      real dbvtit(plond,plevp),      ! intrfc drvtv plnck fnctn for o3
     $     dbvtly(plond,plev)        ! level drvtv plnck fnctn for o3
c
      real dbvt,t     ! planck fnctn tmp derivative for o3
c
      dbvt(t)=(-2.8911366682e-4+(2.3771251896e-6+1.1305188929e-10*t)*t)/
     $  (1.0+(-6.1364820707e-3+1.5550319767e-5*t)*t)
c
c-----------------------------------------------------------------------
c
c initialize
c
      do 2 k1=1,plevp
         do 3 i=1,plon
            dbvtit(i,k1) = dbvt(tint(i,k1))
    3    continue
    2 continue
      do 4 k=1,plev
         do 5 i=1,plon
            dbvtly(i,k) = dbvt(tlayr(i,k+1))
    5    continue
    4 continue
c
      r80257  = 1./8.0257e-04
      r293    = 1./293.
      r250    = 1./250.
      r3205   = 1./.3205
      r300    = 1./300.
      rsslp   = 1./sslp
      r2sslp  = 1./(2.*sslp)
      r2st(1) = 1./(2.*st(1))
      r2st(2) = 1./(2.*st(2))
      bndfct  = 2.0*22.18/(sqrt(196.)*300.)
c
c non-adjacent layer absorptivity:
c
c abso(i,1)     0 -  800 cm-1   h2o rotation band
c abso(i,2)  1200 - 2200 cm-1   h2o vibration-rotation band
c abso(i,3)   800 - 1200 cm-1   h2o window
c abso(i,4)   500 -  800 cm-1   h2o rotation band overlap with co2
c abso(i,5)   ?  cm-1   o3  9.6 micro-meter band
c abso(i,6)   500 -  800 cm-1   co2 15  micro-meter band
c
      do 8 k=1,plevp
         do 9 i=1,plon
            pnmsq(i,k) = pnm(i,k)**2
            dtx(i) = tplnka(i,k) - 250.
            term6(i,k) = coeff(1,2) + coeff(2,2)*dtx(i)*
     $                   (1. +  c9*dtx(i)*(1. + c11*dtx(i)*
     $                   (1. + c13*dtx(i)*(1. + c15*dtx(i)))))
            term9(i,k) = coefi(1,2) + coefi(2,2)*dtx(i)*
     $                   (1. + c19*dtx(i)*(1. + c21*dtx(i)*
     $                   (1. + c23*dtx(i)*(1. + c25*dtx(i)))))
    9    continue
    8 continue
c
c begin non-nearest layer level loops
c
      do 180 k1=plevp,1,-1
         do 170 k2=plevp,1,-1
            do 10 i=1,plon
               dplh2o(i) = plh2o(i,k1) - plh2o(i,k2)
               u(i)      = abs(dplh2o(i))
               sqrtu(i)  = sqrt(u(i))
               ds2c      = abs(s2c(i,k1) - s2c(i,k2))
               dw(i)     = abs(w(i,k1) - w(i,k2))
               uc1(i)    = (ds2c + 1.7e-3*u(i))*(1. +  2.*ds2c)/
     $                                          (1. + 15.*ds2c)
               uc(i)     = ds2c + 2.e-3*u(i)
   10       continue
            if(k1.eq.k2) then
               do 20 i=1,plon
                  pnew(i)   = 0.5
                  tpatha(i) = s2t(i,k1) - s2t(i,k2)
   20          continue
            else                   ! k1.ne.k2
               do 30 i=1,plon
                  pnew(i)   = u(i)/dw(i)
                  tpatha(i) = (s2t(i,k1) - s2t(i,k2))/dplh2o(i)
   30          continue
            end if
            do 40 i=1,plon
               dtx(i)      = tplnka(i,k2) - 250.
               dty(i)      = tpatha(i)    - 250.
               dtyp15(i)   = dty(i) + 15.
               dtyp15sq(i) = dtyp15(i)**2
               dtz(i)      = dtx(i) - 50.
               dtp(i)      = dty(i) - 50.
   40       continue
            do 60 iband=2,4,2
               do 50 i=1,plon
                  term1(i,iband) = coefe(1,iband) + coefe(2,iband)*
     $                             dtx(i)*(1. + c1(iband)*dtx(i))
                  term2(i,iband) = coefb(1,iband) + coefb(2,iband)*
     $                             dtx(i)*(1. + c2(iband)*dtx(i)*
     $                                     (1. + c3(iband)*dtx(i)))
                  term3(i,iband) = coefd(1,iband) + coefd(2,iband)*
     $                             dtx(i)*(1. + c4(iband)*dtx(i)*
     $                                     (1. + c5(iband)*dtx(i)))
                  term4(i,iband) = coefa(1,iband) + coefa(2,iband)*
     $                             dty(i)*(1. + c6(iband)*dty(i))
                  term5(i,iband) = coefc(1,iband) + coefc(2,iband)*
     $                             dty(i)*(1. + c7(iband)*dty(i))
   50       continue
   60    continue
            do 70 i=1,plon
c
c abso(i,1)     0 -  800 cm-1   h2o rotation band
c
               a11 = 0.44 + 3.380e-4*dtz(i) - 1.520e-6*dtz(i)*dtz(i)
               a31 = 1.05 - 6.000e-3*dtp(i) + 3.000e-6*dtp(i)*dtp(i)
               a21 = 1.00 + 1.717e-3*dtz(i) - 1.133e-5*dtz(i)*dtz(i)
               a22 = 1.00 + 4.443e-3*dtp(i) + 2.750e-5*dtp(i)*dtp(i)
               a23 = 1.00 + 3.600*sqrtu(i)
               corfac  = a31*(a11 + ((2.*a21*a22)/a23))
               t1t4    = term1(i,2)*term4(i,2)
               t2t5    = term2(i,2)*term5(i,2)
               a       = t1t4 + t2t5/(1. + t2t5*sqrtu(i)*corfac)
               fwk(i)  = fwcoef + fwc1/(1. + fwc2*u(i))
               fwku(i) = fwk(i)*u(i)
               rsum    = exp(-a*(sqrtu(i) + fwku(i)))
               abso(i,1) = (1. - rsum)*term3(i,2)
               trab1(i)  = rsum
   70       continue
            do 72 i=1,plon
c
c abso(i,2)  1200 - 2200 cm-1   h2o vibration-rotation band
c
               a41    = 1.75 - 3.960e-03*dtz(i)
               a51    = 1.00 + 1.3*sqrtu(i)
               a61    = 1.00 + 1.250e-03*dtp(i) 
     $                + 6.250e-05*dtp(i)*dtp(i)
               corfac = .29*(1. + a41/a51)*a61
               t1t4   = term1(i,4)*term4(i,4)
               t2t5   = term2(i,4)*term5(i,4)
               a      = t1t4 + t2t5/(1. + t2t5*sqrtu(i)*corfac)
               rsum   = exp(-a*(sqrtu(i) + fwku(i)))
               abso(i,2) = (1. - rsum)*term3(i,4)
               trab7(i)  = rsum
   72       continue
c
c line transmission in 800-1000 and 1000-1200 cm-1 intervals
c
            do 90 k=1,2
               do 80 i=1,plon
                  phi   = exp(a1(k)*dtyp15(i) + a2(k)*dtyp15sq(i))
                  psi   = exp(b1(k)*dtyp15(i) + b2(k)*dtyp15sq(i))
                  ubar  = dw(i)*phi*1.66*r80257
                  pbar  = pnew(i)*(psi/phi)
                  cf812 = cfa1 + (1. - cfa1)/(1. + ubar*pbar*10.)
                  g2    = 1. + ubar*4.0*st(k)*cf812/pbar
                  g4    = realk(k)*pbar*r2st(k)*(sqrt(g2) - 1.)
                  trline(i,k) = exp(-g4)
   80          continue
   90       continue
            do 100 i=1,plon
               term7(i,1) = coefj(1,1) + coefj(2,1)*dty(i)*
     $                                   (1. + c16*dty(i))
               term8(i,1) = coefk(1,1) + coefk(2,1)*dty(i)*
     $                                   (1. + c17*dty(i))
               term7(i,2) = coefj(1,2) + coefj(2,2)*dty(i)*
     $                                   (1. + c26*dty(i))
               term8(i,2) = coefk(1,2) + coefk(2,2)*dty(i)*
     $                                   (1. + c27*dty(i))
  100       continue
            do 110 i=1,plon
c
c abso(i,3)   800 - 1200 cm-1   h2o window
c abso(i,4)   500 -  800 cm-1   h2o rotation band overlap with co2
c
               k21    = term7(i,1) + term8(i,1)/
     $             (1. + (c30 + c31*(dty(i)-10.)*(dty(i)-10.))*sqrtu(i))
               k22    = term7(i,2) + term8(i,2)/
     $             (1. + (c28 + c29*(dty(i)-10.))*sqrtu(i))
               tr1    = exp(-(k21*(sqrtu(i) + fc1*fwku(i))))
               tr2    = exp(-(k22*(sqrtu(i) + fc1*fwku(i))))
               tr5    = exp(-((coefh(1,3) + coefh(2,3)*dtx(i))*uc1(i)))
               tr6    = exp(-((coefh(1,4) + coefh(2,4)*dtx(i))*uc1(i)))
               tr9(i)   = tr1*tr5
               tr10(i)  = tr2*tr6
               trab2(i) = 0.65*tr9(i) + 0.35*tr10(i)
               trab4(i) = exp(-(coefg(1,3) + coefg(2,3)*dtx(i))*uc(i))
               trab6(i) = exp(-(coefg(1,4) + coefg(2,4)*dtx(i))*uc(i))
               abso(i,3) = term6(i,k2)*(1. - .5*trab4(i)*trline(i,2) -
     $                                       .5*trab6(i)*trline(i,1))
               abso(i,4) = term9(i,k2)*.5*(tr1 - tr9(i) + tr2 - tr10(i))
  110       continue
            if(k1.eq.k2) go to 170
            if(k2.lt.k1) then
               do 120 i=1,plon
                  to3h2o(i) = h2otr(i,k1)/h2otr(i,k2)
  120          continue
            else        
               do 130 i=1,plon
                  to3h2o(i) = h2otr(i,k2)/h2otr(i,k1)
  130          continue
            end if
            do 140 i=1,plon
c
c abso(i,5)   ?  cm-1   o3  9.6 micro-meter band
c
               dpnm(i)  = pnm(i,k1) - pnm(i,k2)
               to3co2(i)=(pnm(i,k1)*co2t(i,k1) - pnm(i,k2)*co2t(i,k2))/
     $                   dpnm(i)
               te       = (to3co2(i)*r293)**.7
               dplos    = plos(i,k1) - plos(i,k2)
               dplol    = plol(i,k1) - plol(i,k2)
               u1       = 18.29*abs(dplos)/te
               u2       = .5649*abs(dplos)/te
               rphat    = dplol/dplos
               tlocal   = tint(i,k2)
               tcrfac   = sqrt(tlocal*r250)*te
               beta     = r3205*(rphat + dpfo3*tcrfac)
               realnu   = te/beta
               tmp1     = u1/sqrt(4. + u1*(1. + realnu))
               tmp2     = u2/sqrt(4. + u2*(1. + realnu))
               o3bndi    = 74.*te*alog(1. + tmp1 + tmp2)
               abso(i,5) = o3bndi*to3h2o(i)*dbvtit(i,k2)
               trab5(i)  = 1.-(o3bndi/(1060-980.))
  140       continue
            do 150 i=1,plon
c
c abso(i,6)   500 -  800 cm-1   co2 15  micro-meter band
c
               sqwp      = sqrt(abs(plco2(i,k1) - plco2(i,k2)))
               et        = exp(-480./to3co2(i))
               sqti(i)   = sqrt(to3co2(i))
               rsqti     = 1./sqti(i)
               et2       = et*et
               et4       = et2*et2
               omet      = 1. - 1.5*et2
               f1co2     = 899.70*omet*
     $                   (1. + 1.94774*et + 4.73486*et2)*rsqti
               f1sqwp(i) = f1co2*sqwp
               t1co2(i)  = 1./(1. + (245.18*omet*sqwp*rsqti))
               oneme     = 1. - et2
               alphat    = oneme**3*rsqti
               pi        = abs(dpnm(i))
               wco2      =  2.5221*co2vmr*pi*rga
               u7        =  4.9411e4*alphat*et2*wco2
               u8        =  3.9744e4*alphat*et4*wco2
               u9        =  1.0447e5*alphat*et4*et2*wco2
               u13       = 2.8388e3*alphat*et4*wco2
               tpath     = to3co2(i)
               tlocal    = tint(i,k2)
               tcrfac    = sqrt(tlocal*r250*tpath*r300)
               posqt     = ((pnm(i,k2) + pnm(i,k1))*r2sslp +
     $                     dpfco2*tcrfac)*rsqti
               rbeta7    = 1./(5.3228*posqt)
               rbeta8    = 1./(10.6576*posqt)
               rbeta9    = rbeta7
               rbeta13   = rbeta9
               f2co2(i)  = (u7/sqrt(4. + u7*(1. + rbeta7))) +
     $                     (u8/sqrt(4. + u8*(1. + rbeta8))) +
     $                     (u9/sqrt(4. + u9*(1. + rbeta9)))
               f3co2(i)  = u13/sqrt(4. + u13*(1. + rbeta13))
  150       continue
            if(k2.ge.k1) then
               do 152 i=1,plon
                  sqti(i) = sqrt(tlayr(i,k2))
  152          continue
            end if
c
            do 160 i=1,plon
               tmp1      = alog(1. + f1sqwp(i))
               tmp2      = alog(1. + f2co2(i))
               tmp3      = alog(1. + f3co2(i))
               absbnd    = (tmp1 + 2.*t1co2(i)*tmp2 + 2.*tmp3)
     $                     *sqti(i)
               abso(i,6) = trab2(i)*co2em(i,k2)*absbnd
               trab3(i)  = 1. - bndfct*absbnd
  160       continue
            do 165 i=1,plon
c
c sum total absorptivity
c
               abstot(i,k1,k2) = abso(i,1) + abso(i,2) + abso(i,3) +
     $                           abso(i,4) + abso(i,5) + abso(i,6)
  165       continue
c
  170    continue
  180 continue
c
c end of non-nearest layer level loops
c
c
c non-adjacent layer absorptivity:
c
c abso(i,1)     0 -  800 cm-1   h2o rotation band
c abso(i,2)  1200 - 2200 cm-1   h2o vibration-rotation band
c abso(i,3)   800 - 1200 cm-1   h2o window
c abso(i,4)   500 -  800 cm-1   h2o rotation band overlap with co2
c abso(i,5)   ?  cm-1   o3  9.6 micro-meter band
c abso(i,6)   500 -  800 cm-1   co2 15  micro-meter band
c
c begin nearest layer level loop
c
      do 360 k2=plev,1,-1
         do 190 i=1,plon
            tbar(i,1)   = 0.5*(tint(i,k2+1) + tlayr(i,k2+1))
            emm(i,1)    = 0.5*(co2em(i,k2+1) + co2eml(i,k2))
            tbar(i,2)   = 0.5*(tlayr(i,k2+1) + tint(i,k2))
            emm(i,2)    = 0.5*(co2em(i,k2) + co2eml(i,k2))
            tbar(i,3)   = 0.5*(tbar(i,2) + tbar(i,1))
            emm(i,3)    = emm(i,1)
            tbar(i,4)   = tbar(i,3)
            emm(i,4)    = emm(i,2)
            o3emm(i,1)  = 0.5*(dbvtit(i,k2+1) + dbvtly(i,k2))
            o3emm(i,2)  = 0.5*(dbvtit(i,k2) + dbvtly(i,k2))
            o3emm(i,3)  = o3emm(i,1)
            o3emm(i,4)  = o3emm(i,2)
            temh2o(i,1) = tbar(i,1)
            temh2o(i,2) = tbar(i,2)
            temh2o(i,3) = tbar(i,1)
            temh2o(i,4) = tbar(i,2)
            dpnm(i)     = pnm(i,k2+1) - pnm(i,k2)
  190    continue
         do 205 i=1,plon
            rdpnmsq    = 1./(pnmsq(i,k2+1) - pnmsq(i,k2))
            rdpnm      = 1./dpnm(i)
            p1         = .5*(pbr(i,k2) + pnm(i,k2+1))
            p2         = .5*(pbr(i,k2) + pnm(i,k2  ))
            uinpl(i,1) =  (pnmsq(i,k2+1) - p1**2)*rdpnmsq
            uinpl(i,2) = -(pnmsq(i,k2  ) - p2**2)*rdpnmsq
            uinpl(i,3) = -(pnmsq(i,k2  ) - p1**2)*rdpnmsq
            uinpl(i,4) =  (pnmsq(i,k2+1) - p2**2)*rdpnmsq
            winpl(i,1) = (.5*( pnm(i,k2+1) - pbr(i,k2)))*rdpnm
            winpl(i,2) = (.5*(-pnm(i,k2  ) + pbr(i,k2)))*rdpnm
            winpl(i,3) = (.5*( pnm(i,k2+1) + pbr(i,k2)) - pnm(i,k2  ))*
     $                   rdpnm
            winpl(i,4) = (.5*(-pnm(i,k2  ) - pbr(i,k2)) + pnm(i,k2+1))*
     $                   rdpnm
            tmp1       = 1./(piln(i,k2+1) - piln(i,k2))
            tmp2       = piln(i,k2+1) - pmln(i,k2)
            tmp3       = piln(i,k2  ) - pmln(i,k2)
            zinpl(i,1) = (.5*tmp2          )*tmp1
            zinpl(i,2) = (        - .5*tmp3)*tmp1
            zinpl(i,3) = (.5*tmp2 -    tmp3)*tmp1
            zinpl(i,4) = (   tmp2 - .5*tmp3)*tmp1
            pinpl(i,1) = 0.5*(p1 + pnm(i,k2+1))
            pinpl(i,2) = 0.5*(p2 + pnm(i,k2  ))
            pinpl(i,3) = 0.5*(p1 + pnm(i,k2  ))
            pinpl(i,4) = 0.5*(p2 + pnm(i,k2+1))
  205    continue
         do  350 kn=1,4
            do 210 i=1,plon
               u(i)     = uinpl(i,kn)*abs(plh2o(i,k2) - plh2o(i,k2+1))
               sqrtu(i) = sqrt(u(i))
               dw(i)    = abs(w(i,k2) - w(i,k2+1))
               pnew(i)  = u(i)/(winpl(i,kn)*dw(i))
               ds2c     = abs(s2c(i,k2) - s2c(i,k2+1))
               uc1(i)   = uinpl(i,kn)*ds2c
               uc1(i)   = (uc1(i) + 1.7e-3*u(i))*(1. + 2.*uc1(i))/
     $                                         (1. + 15.*uc1(i))
               uc(i)    = uinpl(i,kn)*ds2c + 2.e-3*u(i)
  210       continue
            do 230 i=1,plon
               dtx(i)      = temh2o(i,kn) - 250.
               dty(i)      = tbar(i,kn) - 250.
               dtyp15(i)   = dty(i) + 15.
               dtyp15sq(i) = dtyp15(i)**2
               dtz(i)      = dtx(i) - 50.
               dtp(i)      = dty(i) - 50.
  230       continue
            do 270 iband=2,4,2
               do 260 i=1,plon
                  term1(i,iband) = coefe(1,iband) + coefe(2,iband)*
     $                             dtx(i)*(1. + c1(iband)*dtx(i))
                  term2(i,iband) = coefb(1,iband) + coefb(2,iband)*
     $                             dtx(i)*(1. + c2(iband)*dtx(i)*
     $                                     (1. + c3(iband)*dtx(i)))
                  term3(i,iband) = coefd(1,iband) + coefd(2,iband)*
     $                             dtx(i)*(1. + c4(iband)*dtx(i)*
     $                                     (1. + c5(iband)*dtx(i)))
                  term4(i,iband) = coefa(1,iband) + coefa(2,iband)*
     $                             dty(i)*(1. + c6(iband)*dty(i))
                  term5(i,iband) = coefc(1,iband) + coefc(2,iband)*
     $                             dty(i)*(1. + c7(iband)*dty(i))
  260          continue
  270       continue
            do 280 i=1,plon
c
c abso(i,1)     0 -  800 cm-1   h2o rotation band
c
               a11 = 0.44 + 3.380e-4*dtz(i) - 1.520e-6*dtz(i)*dtz(i)
               a31 = 1.05 - 6.000e-3*dtp(i) + 3.000e-6*dtp(i)*dtp(i)
               a21 = 1.00 + 1.717e-3*dtz(i) - 1.133e-5*dtz(i)*dtz(i)
               a22 = 1.00 + 4.443e-3*dtp(i) + 2.750e-5*dtp(i)*dtp(i)
               a23 = 1.00 + 3.600*sqrtu(i)
               corfac    = a31*(a11 + ((2.*a21*a22)/a23))
               t1t4      = term1(i,2)*term4(i,2)
               t2t5      = term2(i,2)*term5(i,2)
               a         = t1t4 + t2t5/(1. + t2t5*sqrtu(i)*corfac)
               fwk(i)    = fwcoef + fwc1/(1. + fwc2*u(i))
               fwku(i)   = fwk(i)*u(i)
               rsum      = exp(-a*(sqrtu(i) + fwku(i)))
               abso(i,1) = (1. - rsum)*term3(i,2)
               trab1(i)  = rsum
  280       continue
            do 282 i=1,plon
c
c abso(i,2)  1200 - 2200 cm-1   h2o vibration-rotation band
c
               a41       = 1.75 - 3.960e-03*dtz(i)
               a51       = 1.00 + 1.3*sqrtu(i)
               a61       = 1.00 + 1.250e-03*dtp(i) 
     $                   + 6.250e-05*dtp(i)*dtp(i)
               corfac    = .29*(1. + a41/a51)*a61
               t1t4      = term1(i,4)*term4(i,4)
               t2t5      = term2(i,4)*term5(i,4)
               a         = t1t4 + t2t5/(1. + t2t5*sqrtu(i)*corfac)
               rsum      = exp(-a*(sqrtu(i) + fwku(i)))
               abso(i,2) = (1. - rsum)*term3(i,4)
               trab7(i)  = rsum
  282       continue
c
c line transmission in 800-1000 and 1000-1200 cm-1 intervals
c
            do 300 k=1,2
               do 290 i=1,plon
                  phi   = exp(a1(k)*dtyp15(i) + a2(k)*dtyp15sq(i))
                  psi   = exp(b1(k)*dtyp15(i) + b2(k)*dtyp15sq(i))
                  ubar  = dw(i)*phi*winpl(i,kn)*1.66*r80257
                  pbar  = pnew(i)*(psi/phi)
                  cf812 = cfa1 + (1. - cfa1)/(1. + ubar*pbar*10.)
                  g2    = 1. + ubar*4.0*st(k)*cf812/pbar
                  g4    = realk(k)*pbar*r2st(k)*(sqrt(g2) - 1.)
                  trline(i,k) = exp(-g4)
  290          continue
  300       continue
            do 310 i=1,plon
               term7(i,1) = coefj(1,1) + coefj(2,1)*dty(i)*
     $                                   (1. + c16*dty(i))
               term8(i,1) = coefk(1,1) + coefk(2,1)*dty(i)*
     $                                   (1. + c17*dty(i))
               term7(i,2) = coefj(1,2) + coefj(2,2)*dty(i)*
     $                                   (1. + c26*dty(i))
               term8(i,2) = coefk(1,2) + coefk(2,2)*dty(i)*
     $                                   (1. + c27*dty(i))
  310       continue
            do 320 i=1,plon
c
c abso(i,3)   800 - 1200 cm-1   h2o window
c abso(i,4)   500 -  800 cm-1   h2o rotation band overlap with co2
c
               dtym10    = dty(i) - 10.
               denom     = 1. + (c30 + c31*dtym10*dtym10)*sqrtu(i)
               k21       = term7(i,1) + term8(i,1)/denom
               denom     = 1. + (c28 + c29*dtym10       )*sqrtu(i)
               k22       = term7(i,2) + term8(i,2)/denom
               term9(i,2) = coefi(1,2) + coefi(2,2)*dtx(i)*
     $                     (1. + c19*dtx(i)*(1. + c21*dtx(i)*
     $                     (1. + c23*dtx(i)*(1. + c25*dtx(i)))))
               tr1     = exp(-(k21*(sqrtu(i) + fc1*fwku(i))))
               tr2     = exp(-(k22*(sqrtu(i) + fc1*fwku(i))))
               tr5     = exp(-((coefh(1,3) + coefh(2,3)*dtx(i))*uc1(i)))
               tr6     = exp(-((coefh(1,4) + coefh(2,4)*dtx(i))*uc1(i)))
               tr9(i)  = tr1*tr5
               tr10(i) = tr2*tr6
               trab2(i)= 0.65*tr9(i) + 0.35*tr10(i)
               trab4(i)= exp(-(coefg(1,3) + coefg(2,3)*dtx(i))*uc(i))
               trab6(i)= exp(-(coefg(1,4) + coefg(2,4)*dtx(i))*uc(i))
               term6(i,2) = coeff(1,2) + coeff(2,2)*dtx(i)*
     $                     (1. + c9*dtx(i)*(1. + c11*dtx(i)*
     $                     (1. + c13*dtx(i)*(1. + c15*dtx(i)))))
               abso(i,3)  = term6(i,2)*(1. - .5*trab4(i)*trline(i,2) -
     $                                       .5*trab6(i)*trline(i,1))
               abso(i,4)  = term9(i,2)*.5*(tr1 - tr9(i) + tr2 - tr10(i))
  320       continue
            do 330 i=1,plon
c
c abso(i,5)   ?  cm-1   o3  9.6 micro-meter band
c
               te        = (tbar(i,kn)*r293)**.7
               dplos     = abs(plos(i,k2+1) - plos(i,k2))
               u1        = zinpl(i,kn)*18.29*dplos/te
               u2        = zinpl(i,kn)*.5649*dplos/te
               tlocal    = tbar(i,kn)
               tcrfac    = sqrt(tlocal*r250)*te
               beta      = r3205*(pinpl(i,kn)*rsslp + dpfo3*tcrfac)
               realnu    = te/beta
               tmp1      = u1/sqrt(4. + u1*(1. + realnu))
               tmp2      = u2/sqrt(4. + u2*(1. + realnu))
               o3bndi    = 74.*te*alog(1. + tmp1 + tmp2)
               abso(i,5) = o3bndi*o3emm(i,kn)
     $                     *(h2otr(i,k2+1)/h2otr(i,k2))
               trab5(i)  = 1.-(o3bndi/(1060-980.))
  330       continue
            do 340 i=1,plon
c
c abso(i,6)   500 -  800 cm-1   co2 15  micro-meter band
c
               dplco2   = plco2(i,k2+1) - plco2(i,k2)
               sqwp     = sqrt(uinpl(i,kn)*dplco2)
               et       = exp(-480./tbar(i,kn))
               sqti(i)  = sqrt(tbar(i,kn))
               rsqti    = 1./sqti(i)
               et2      = et*et
               et4      = et2*et2
               omet     = (1. - 1.5*et2)
               f1co2    = 899.70*omet*
     $                    (1. + 1.94774*et + 4.73486*et2)*rsqti
               f1sqwp(i)= f1co2*sqwp
               t1co2(i) = 1./(1. + (245.18*omet*sqwp*rsqti))
               oneme    = 1. - et2
               alphat   = oneme**3*rsqti
               pi       = abs(dpnm(i))*winpl(i,kn)
               wco2     = 2.5221*co2vmr*pi*rga
               u7       = 4.9411e4*alphat*et2*wco2
               u8       = 3.9744e4*alphat*et4*wco2
               u9       = 1.0447e5*alphat*et4*et2*wco2
               u13      = 2.8388e3*alphat*et4*wco2
               tpath    = tbar(i,kn)
               tlocal   = tbar(i,kn)
               tcrfac   = sqrt((tlocal*r250)*(tpath*r300))
               posqt    = (pinpl(i,kn)*rsslp + dpfco2*tcrfac)*rsqti
               rbeta7   = 1./(5.3228*posqt)
               rbeta8   = 1./(10.6576*posqt)
               rbeta9   = rbeta7
               rbeta13  = rbeta9
               f2co2(i) = u7/sqrt(4. + u7*(1. + rbeta7)) +
     $                    u8/sqrt(4. + u8*(1. + rbeta8)) +
     $                    u9/sqrt(4. + u9*(1. + rbeta9))
               f3co2(i) = u13/sqrt(4. + u13*(1. + rbeta13))
               tmp1     = alog(1. + f1sqwp(i))
               tmp2     = alog(1. + f2co2(i))
               tmp3     = alog(1. + f3co2(i))
               absbnd   = (tmp1 + 2.*t1co2(i)*tmp2 + 2.*tmp3)
     $                    *sqti(i)
               abso(i,6)= trab2(i)*emm(i,kn)*absbnd
               trab3(i) = 1. - bndfct*absbnd
  340       continue
c
c compute total next layer absorptivity:
c
         do 370 i=1,plon
            absnxt(i,k2,kn) = abso(i,1) + abso(i,2) + abso(i,3) +
     $                        abso(i,4) + abso(i,5) + abso(i,6)
  370    continue
c
  350    continue
  360 continue
c
c end of nearest layer level loop
c
c
c done
c
      return
      end
      subroutine radems(s2c    ,s2t    ,w     ,tplnke ,plh2o ,
     $                  pnm    ,plco2  ,tint  ,tint4  ,tlayr ,
     $                  tlayr4 ,plol   ,plos  ,co2em  ,co2eml,
     $                  co2t   ,h2otr  ,emstot)
c-----------------------------------------------------------------------
c
c compute emissivity for h2o, co2, o3
c
c 
c h2o  ....  uses nonisothermal emissivity for water vapor from
c            Ramanathan, V. and  P.Downey, 1986: A Nonisothermal
c            Emissivity and Absorptivity Formulation for Water Vapor
c            Jouranl of Geophysical Research, vol. 91., D8, pp 8649-8666
c
c
c co2  ....  uses absorptance parameterization of the 15 micro-meter
c            (500 - 800 cm-1) band system of Carbon Dioxide, from
c            Kiehl, J.T. and B.P.Briegleb, 1991: A New Parameterization
c            of the Absorptance Due to the 15 micro-meter Band System
c            of Carbon Dioxide Jouranl of Geophysical Research, 
c            vol. 96., D5, pp 9013-9019
c
c o3   ....  uses absorptance parameterization of the 9.6 micro-meter 
c            band system of ozone, from Ramanathan, V. and R.E.Dickinson,
c            1979: The Role of stratospheric ozone in the zonal and 
c            seasonal radiative energy balance of the earth-troposphere
c            system. Journal of the Atmospheric Sciences, Vol. 36, 
c            pp 1084-1104
c
c 
c computes individual emissivities, accounting for band overlap, and
c sums to obtain the total.
c
c-----------------------------------------------------------------------
      implicit none
c-----------------------------------------------------------------------
c
c     !   Basic grid point resolution parameters
c
      integer plev,plon,plat,plonp2,plevp,pcnst,plond,platd,plevd
      integer nxpt,jintmx,i1,j1,j1m
c
      parameter(plev=  133,        ! number of vertical levels
     $          plevp=plev+1,     ! number of vertical interfaces
     $          plon=1,           ! number of longitudes (T42)
     $          plat=1,           ! number of latitudes (T42)
     $          plonp2=plon+2,    ! length necessary for fft
     $          pcnst=1,          ! number of constituents
     $          nxpt=1,           !
     $          jintmx=2,         !
     $          plond=plon + 1 + 2*nxpt, ! slt extended domain longitude
     $          platd=plat + 2*nxpt + 2*jintmx, ! slt extended domain lat.
     $          plevd=plev*(3+pcnst), ! fold plev,pcnst indices into one
     $          i1=1+nxpt,        ! model starting longitude (3-d)
     $          j1=1+nxpt+jintmx, ! model starting latitude (3-d)
     $          j1m=j1-1)         ! model starting offset (3-d)
c
c-----------------------------------------------------------------------
c
c water vapor narrow band constants for lw computations
c
      real realk,st,a1,a2,b1,b2,
     $     coefa,coefb,coefc,coefd,
     $     coefe,coeff,coefg,coefh,
     $     coefi,coefj,coefk,
     $     c1,c2,c3,c4,c5,c6,c7,
     $     c8 ,c9 ,c10,c11,c12,c13,c14,c15,c16,c17,
     $     c18,c19,c20,c21,c22,c23,c24,c25,c26,c27,
     $     c28,c29,c30,c31,
     $     fwcoef,fwc1,fwc2,fc1,cfa1
c
      common/crdcae/realk(2), st(2), a1(2), a2(2), b1(2), b2(2),
c
c constant coefficients for water vapor absorptivity and emissivi
c
     $              coefa(3,4),coefb(4,4),coefc(3,4),coefd(4,4),
     $              coefe(3,4),coeff(6,2),coefg(2,4),coefh(2,4),
     $              coefi(6,2),coefj(3,2),coefk(3,2),
     $              c1(4),c2(4),c3(4),c4(4),c5(4),c6(4),c7(4),
     $              c8 ,c9 ,c10,c11,c12,c13,c14,c15,c16,c17,
     $              c18,c19,c20,c21,c22,c23,c24,c25,c26,c27,
     $              c28,c29,c30,c31,
c
c farwing correction constants for narrow-band emissivity model
c introduce farwing correction to account for the
c deficiencies in narrow-band model used to derive the
c emissivity. tuned with arkings line calculations.
c
     $              fwcoef,fwc1,fwc2,fc1,cfa1
c
c-----------------------------------------------------------------------
c
c radiation constants
c
      real gravit,   ! gravitational acceleration
     $        rga,   ! 1 over gravit
     $      cpair,   ! heat capacity air at constant pressure
     $     epsilo,   ! ratio mmw h2o to mmw air
     $       sslp,   ! standard pressure
     $     stebol,   ! stephan boltzmann constant
     $     rgsslp,   ! 0.5 / (gravit*sslp)
     $     co2vmr,   ! co2 volume mixing ratio
     $      dpfo3,   ! Doppler factor for o3
     $     dpfco2,   ! Doppler factor for co2
     $     dayspy,   ! solar days in one year
     $        pie    ! pie
c
      common/crdcon/gravit,    rga,  cpair,  epsilo,   sslp,
     $              stebol, rgsslp, co2vmr,   dpfo3, dpfco2,
     $              dayspy,    pie
c
c-----------------------------------------------------------------------
c
c input arguments
c
      real s2c(plond,plevp),    ! h2o continuum path length
     $     s2t(plond,plevp),    ! tmp and prs wghted h2o path length
     $     w(plond,plevp),      ! h2o path length
     $     tplnke(plond),       ! layer planck temperature
     $     plh2o(plond,plevp),  ! h2o prs wghted path length
     $     pnm(plond,plevp),    ! model interface pressure
     $     plco2(plond,plevp),  ! prs wghted path of co2
     $     tint(plond,plevp),   ! model interface temperatures
     $     tint4(plond,plevp),  ! tint to the 4th power
     $     tlayr(plond,plevp),  ! k-1 model layer temperature
     $     tlayr4(plond,plevp), ! tlayr to the 4th power
     $     plol(plond,plevp),   ! pressure wghtd ozone path
     $     plos(plond,plevp)    ! ozone path
c
c output arguments
c
      real emstot(plond,plevp),  ! total emissivity
     $     co2em(plond,plevp),   ! layer co2 normalized plnck function drvtv
     $     co2eml(plond,plev),   ! intrfc co2 normalized plnck function drvtv
     $     co2t(plond,plevp),    ! tmp and prs weighted path length
     $     h2otr(plond,plevp)    ! h2o transmission over o3 band
c
c local workspace for h2o:
c
      integer   i,            ! longitude index
     $          k,            ! level index]
     $         k1,            ! level index
     $      iband             ! h2o band index
c
      real h2oems(plond,plevp),! h2o emissivity
     $     tpathe(plond),     ! used to compute h2o emissivity
     $     a(plond),          ! "
     $     corfac(plond),     ! "
     $     dtp(plond),        ! "
     $     dtx(plond),        ! "
     $     dty(plond),        ! "
     $     dtz(plond),        ! "
     $     emis(plond,4),     ! "
     $     rsum(plond),       ! "
     $     term1(plond,4),    ! "
     $     term2(plond,4)     ! "
      real term3(plond,4),    ! "
     $     term4(plond,4),    ! "
     $     term5(plond,4),    ! "
     $     term6(plond,2),    ! "
     $     term7(plond,2),    ! "
     $     term8(plond,2),    ! "
     $     term9(plond,2),    ! "
     $     tr1(plond),        ! "
     $     tr2(plond),        ! "
     $     tr3(plond)         ! "
      real tr4(plond),        ! "
     $     tr7(plond),        ! "
     $     tr8(plond),        ! "
     $     uc(plond),         ! "
     $     pnew(plond),       ! "
     $     trline(plond,2),   ! "
     $     k21(plond),        ! "
     $     k22(plond),        ! "
     $     u(plond),          ! "
     $     uc1(plond),        ! "
     $        r80257 
      real      a11,          ! used to compute h2o emissivity
     $          a31,          ! "
     $          a21,          ! "
     $          a22,          ! "
     $          a23,          ! "
     $          t1t4,         ! "
     $          t2t5,         ! "
     $          fwk,          ! "
     $          a41,          ! "
     $          a51,          ! "
     $          a61,          ! "
     $          phi,          ! "
     $          psi,          ! "
     $          ubar,         ! "
     $          g1,           ! "
     $          pbar,         ! "
     $          g3,           ! "
     $          g2,           ! "
     $          g4,           ! "
     $       cf812            ! "
      real troco2(plond,plevp)  ! h2o overlap factor for co2 absorption
c
c local workspace for co2:
c
      real co2ems(plond,plevp), ! co2 emissivity
     $     co2plk(plond),    ! used to compute co2 emissivity
     $        sum(plond),    ! "
     $               t1i,    ! "
     $              sqti,    ! "
     $                pi,    ! "
     $                et,    ! "
     $               et2,    ! "
     $               et4,    ! "
     $              omet,    ! "
     $                ex     ! "
      real         f1co2,    ! "
     $             f2co2,    ! "
     $             f3co2,    ! "
     $             t1co2,    ! "
     $              sqwp,    ! "
     $            f1sqwp,    ! "
     $             oneme,    ! "
     $            alphat,    ! "
     $              wco2,    ! "
     $             posqt,    ! "
     $            rbeta7,    ! "
     $            rbeta8,    ! "
     $            rbeta9,    ! "
     $           rbeta13     ! "
      real         tpath,    ! "
     $              tmp1,    ! "
     $              tmp2,    ! "
     $              tmp3,    ! "
     $            tlayr5,    ! "
     $             rsqti,    ! "
     $            exm1sq     ! "
c
      real u7,    ! absorber amount for various co2 band systems
     $     u8,    !    "
     $     u9,    !    "
     $     u13    !    "
c
      real r250,  ! inverse 250K
     $     r300,  ! inverse 300K
     $    rsslp   ! inverse standard sea-level pressure
c
c local workspace for o3:
c
      real o3ems(plond,plevp),  ! ozone emissivity 
     $     dbvtt(plond),   ! tmp drvtv of planck fctn for tplnke
     $               te,   ! temperature factor
     $               u1,   ! path length factor
     $               u2,   ! path length factor
     $             phat,   ! effecitive path length pressure
     $           tlocal,   ! local planck function temperature
     $           tcrfac,   ! scaled temperature factor
     $             beta,   ! absorption function factor with voigt effect
     $           realnu,   ! absorption function factor
     $           o3bndi    ! band absorption factor
c
c transmission terms for various spectral intervals:
c
      real    trem1(plond),  ! h2o     0 -  800 cm-1
     $        trem2(plond),  ! h2o   500 -  800 cm-1
     $        trem3(plond),  ! co2   500 -  800 cm-1
     $        trem4(plond),  ! h2o   800 - 1000 cm-1
     $        trem5(plond),  ! o3     9.6 micro-meter band
     $        trem6(plond),  ! h2o  1000 - 1200 cm-1
     $        trem7(plond)   ! h2o  1200 - 2200 cm-1
c
      real bndfct, ! band absorptance parameter for co2
     $     absbnd  ! proportional to co2 band absorptance
c
c derivative of planck function at 9.6 micro-meter wavelength, and
c an absorption function factor:
c
      real dbvt,fo3,t,ux,vx
c
      dbvt(t)=(-2.8911366682e-4+(2.3771251896e-6+1.1305188929e-10*t)*t)/
     $  (1.0+(-6.1364820707e-3+1.5550319767e-5*t)*t)
c
      fo3(ux,vx)=ux/sqrt(4.+ux*(1.+vx))
c
c-----------------------------------------------------------------------
c
c initialize
c
      r80257  = 1./8.0257e-04
c
      r250  = 1./250.
      r300  = 1./300.
      rsslp = 1./sslp
c
c planck function for co2
c
      do 110 i=1,plon
         ex        = exp(960./tplnke(i))
         co2plk(i) = 5.e8/((tplnke(i)**4)*(ex - 1.))
         co2t(i,1) = tplnke(i)
         sum(i)    = co2t(i,1)*pnm(i,1)
  110 continue
      k = 1
      do 140 k1=plevp,2,-1
         k = k + 1
         do 130 i=1,plon
            sum(i)    = sum(i) + tlayr(i,k)*(pnm(i,k)-pnm(i,k-1))
            ex        = exp(960./tlayr(i,k1))
            tlayr5    = tlayr(i,k1)*tlayr4(i,k1)
            co2eml(i,k1-1) 
     $                = 1.2e11*ex/(tlayr5*(ex - 1.)**2)
            co2t(i,k) = sum(i)/pnm(i,k)
  130    continue
  140 continue
      bndfct = 2.0*22.18/(sqrt(196.)*300.)
c
c initialize planck function derivative for o3 
c
      do 210 i=1,plon
         dbvtt(i) = dbvt(tplnke(i))
  210 continue
c
c begin interface loop
c
      do 100 k1=1,plevp
c-------------------------------------------------------------------------
c
c h2o emissivity
c
c emis(i,1)     0 -  800 cm-1   rotation band 
c emis(i,2)  1200 - 2200 cm-1   vibration-rotation band 
c emis(i,3)   800 - 1200 cm-1   window
c emis(i,4)   500 -  800 cm-1   rotation band overlap with co2
c
c for the p type continuum
c
         do 10 i=1,plon
            uc(i)     = s2c(i,k1) + 2.e-3*plh2o(i,k1)
            u(i)      = plh2o(i,k1)
            pnew(i)   = u(i)/w(i,k1)
c
c apply scaling factor for 500-800 continuum
c
            uc1(i)    = (s2c(i,k1) + 1.7e-3*plh2o(i,k1))*
     $                 (1. + 2.*s2c(i,k1))/(1. + 15.*s2c(i,k1))
            tpathe(i) = s2t(i,k1)/plh2o(i,k1)
   10    continue
         do 20 i=1,plon
            dtx(i) = tplnke(i) - 250.
            dty(i) = tpathe(i) - 250.
            dtz(i) = dtx(i) - 50.
            dtp(i) = dty(i) - 50.
   20    continue
         do 40 iband=1,3,2
            do 30 i=1,plon
               term1(i,iband) = coefe(1,iband) + coefe(2,iband)*
     $                          dtx(i)*(1. + c1(iband)*dtx(i))
               term2(i,iband) = coefb(1,iband) + coefb(2,iband)*
     $                          dtx(i)*(1. + c2(iband)*dtx(i)*
     $                                  (1. + c3(iband)*dtx(i)))
               term3(i,iband) = coefd(1,iband) + coefd(2,iband)*
     $                          dtx(i)*(1. +  c4(iband)*dtx(i)*
     $                                  (1. + c5(iband)*dtx(i)))
               term4(i,iband) = coefa(1,iband) + coefa(2,iband)*
     $                          dty(i)*(1. + c6(iband)*dty(i))
               term5(i,iband) = coefc(1,iband) + coefc(2,iband)*
     $                          dty(i)*(1. + c7(iband)*dty(i))
   30       continue
   40    continue
         do 50 i=1,plon
c
c emis(i,1)     0 -  800 cm-1   rotation band 
c
            a11  = .37 - 3.33e-5*dtz(i) + 3.33e-6*dtz(i)*dtz(i)
            a31  = 1.07 - 1.00e-3*dtp(i) + 1.475e-5*dtp(i)*dtp(i)
            a21  = 1.3870 + 3.80e-3*dtz(i) - 7.8e-6*dtz(i)*dtz(i)
            a22  = 1.0 - 1.21e-3*dtp(i) - 5.33e-6*dtp(i)*dtp(i)
            a23  = 0.9 + 2.62*sqrt(u(i))
            corfac(i) = a31*(a11 + ((a21*a22)/a23))
            t1t4 = term1(i,1)*term4(i,1)
            t2t5 = term2(i,1)*term5(i,1)
            a(i) = t1t4 + t2t5/(1. + t2t5*sqrt(u(i))*corfac(i))
            fwk  = fwcoef + fwc1/(1. + fwc2*u(i))
            rsum(i)   = exp(-a(i)*(sqrt(u(i)) + fwk*u(i)))
            emis(i,1) = (1. - rsum(i))*term3(i,1)
            trem1(i)  = rsum(i)
c
c emis(i,2)  1200 - 2200 cm-1   vibration-rotation band 
c
            a41      = 1.75 - 3.96e-3*dtz(i)
            a51      = 1.00 + 1.3*sqrt(u(i))
            a61      = 1.00 + 1.25e-3*dtp(i) + 6.25e-5*dtp(i)*dtp(i)
            corfac(i)= .3*(1. + (a41)/(a51))*a61
            t1t4     = term1(i,3)*term4(i,3)
            t2t5     = term2(i,3)*term5(i,3)
            a(i)     = t1t4 + t2t5/(1. + t2t5*sqrt(u(i))*corfac(i))
            fwk      = fwcoef + fwc1/(1. + fwc2*u(i))
            rsum(i)  = exp(-a(i)*(sqrt(u(i)) + fwk*u(i)))
            emis(i,2)= (1. - rsum(i))*term3(i,3)
            trem7(i) = rsum(i)
   50    continue
c
c line transmission in 800-1000 and 1000-1200 cm-1 intervals
c
         do 70 k=1,2
            do 60 i=1,plon
               phi  = a1(k)*(dty(i) + 15.) + a2(k)*(dty(i) + 15.)**2
               psi  = b1(k)*(dty(i) + 15.) + b2(k)*(dty(i) + 15.)**2
               phi  = exp(phi)
               psi  = exp(psi)
               ubar = w(i,k1)*phi
               ubar = (ubar*1.66)*r80257
               pbar = pnew(i)*(psi/phi)
               cf812 = cfa1 + ((1.-cfa1)/(1. + ubar*pbar*10.))
               g1   = (realk(k)*pbar)/(2.*st(k))
               g2   = 1. + (ubar*4.0*st(k)*cf812)/pbar
               g3   = sqrt(g2) - 1.
               g4   = g1*g3
               trline(i,k) = exp(-g4)
   60       continue
   70    continue
c
         do 80 i=1,plon
            term7(i,1) = coefj(1,1) + coefj(2,1)*dty(i)*(1.+c16*dty(i))
            term8(i,1) = coefk(1,1) + coefk(2,1)*dty(i)*(1.+c17*dty(i))
            term7(i,2) = coefj(1,2) + coefj(2,2)*dty(i)*(1.+c26*dty(i))
            term8(i,2) = coefk(1,2) + coefk(2,2)*dty(i)*(1.+c27*dty(i))
   80    continue
c
         do 90 i=1,plon
c
c emis(i,3)   800 - 1200 cm-1   window
c
            term6(i,1) = coeff(1,1) + coeff(2,1)*dtx(i)*
     $                  (1. +  c8*dtx(i)*(1. + c10*dtx(i)*
     $                  (1. + c12*dtx(i)*(1. + c14*dtx(i)))))
c
            trem4(i)  = exp(-(coefg(1,1)+coefg(2,1)*dtx(i))*uc(i))
     $                  *trline(i,2)
            trem6(i)  = exp(-(coefg(1,2)+coefg(2,2)*dtx(i))*uc(i))
     $                  *trline(i,1)        
c
            emis(i,3) = term6(i,1)*(1. - .5*trem4(i) -.5*trem6(i))
c
c emis(i,4)   500 -  800 cm-1   rotation band overlap with co2
c
            k21(i) = term7(i,1) + term8(i,1)/
     $           (1. + (c30 + c31*(dty(i)-10.)*(dty(i)-10.))*sqrt(u(i)))
            k22(i) = term7(i,2) + term8(i,2)/
     $           (1. + (c28 + c29*(dty(i)-10.))*sqrt(u(i)))
            term9(i,1) = coefi(1,1) + coefi(2,1)*dtx(i)*
     $                  (1. + c18*dtx(i)*(1. + c20*dtx(i)*
     $                   (1. + c22*dtx(i)*(1. + c24*dtx(i)))))
            fwk    = fwcoef + fwc1/(1.+fwc2*u(i))
            tr1(i) = exp(-(k21(i)*(sqrt(u(i)) + fc1*fwk*u(i))))
            tr2(i) = exp(-(k22(i)*(sqrt(u(i)) + fc1*fwk*u(i))))
            tr3(i) = exp(-((coefh(1,1) + coefh(2,1)*dtx(i))*uc1(i)))
            tr4(i) = exp(-((coefh(1,2) + coefh(2,2)*dtx(i))*uc1(i)))
            tr7(i) = tr1(i)*tr3(i)
            tr8(i) = tr2(i)*tr4(i)
            emis(i,4) = term9(i,1)*.5*(tr1(i)-tr7(i) + tr2(i)-tr8(i))
            h2oems(i,k1) = emis(i,1)+emis(i,2)+emis(i,3)+emis(i,4)
            troco2(i,k1) = 0.65*tr7(i) + 0.35*tr8(i)
            trem2(i)     = troco2(i,k1) 
   90    continue
c-----------------------------------------------------------------------
c
c co2 emissivity for 500-800 cm-1
c
         do 150 i=1,plon
c
            t1i    = exp(-480./co2t(i,k1))
            sqti   = sqrt(co2t(i,k1))
            rsqti  = 1./sqti
            et     = t1i
            et2    = et*et
            et4    = et2*et2
            omet   = 1. - 1.5*et2
            f1co2  = 899.70*omet*(1. + 1.94774*et + 4.73486*et2)*rsqti
            sqwp   = sqrt(plco2(i,k1))
            f1sqwp = f1co2*sqwp
            t1co2  = 1./(1. + 245.18*omet*sqwp*rsqti)
            oneme  = 1. - et2
            alphat = oneme**3*rsqti
            wco2   = 2.5221*co2vmr*pnm(i,k1)*rga
            u7     = 4.9411e4*alphat*et2*wco2
            u8     = 3.9744e4*alphat*et4*wco2
            u9     = 1.0447e5*alphat*et4*et2*wco2
            u13    = 2.8388e3*alphat*et4*wco2
c
            tpath  = co2t(i,k1)
            tlocal = tplnke(i)
            tcrfac = sqrt((tlocal*r250)*(tpath*r300))
            pi     = pnm(i,k1)*rsslp + 2.*dpfco2*tcrfac
            posqt  = pi/(2.*sqti)
            rbeta7 =  1./( 5.3288*posqt)
            rbeta8 = 1./ (10.6576*posqt)
            rbeta9 = rbeta7
            rbeta13= rbeta9
            f2co2  = (u7/sqrt(4. + u7*(1. + rbeta7))) +
     $               (u8/sqrt(4. + u8*(1. + rbeta8))) +
     $               (u9/sqrt(4. + u9*(1. + rbeta9)))
            f3co2  = u13/sqrt(4. + u13*(1. + rbeta13))
            tmp1   = alog(1. + f1sqwp)
            tmp2   = alog(1. +  f2co2)
            tmp3   = alog(1. +  f3co2)
            absbnd = (tmp1 + 2.*t1co2*tmp2 + 2.*tmp3)*sqti
            co2ems(i,k1)  = troco2(i,k1)*absbnd*co2plk(i)
            ex     = exp(960./tint(i,k1))
            exm1sq = (ex - 1.)**2
            co2em(i,k1) 
     $             = 1.2e11*ex/(tint(i,k1)*tint4(i,k1)*exm1sq)
            trem3(i) 
     $             = 1. - bndfct*absbnd
c
  150    continue
c-----------------------------------------------------------------------
c
c o3 emissivity
c
         do 220 i=1,plon
            h2otr(i,k1) = exp(-12.*s2c(i,k1))
            te          = (co2t(i,k1)/293.)**.7
            u1          = 18.29*plos(i,k1)/te
            u2          = .5649*plos(i,k1)/te
            phat        = plos(i,k1)/plol(i,k1)
            tlocal      = tplnke(i)
            tcrfac      = sqrt(tlocal/250.)*te
            beta        = (1./.3205)*((1./phat) + (dpfo3*tcrfac))
            realnu      = (1./beta)*te
            o3bndi      = 74.*te*(tplnke(i)/375.)*
     $         alog(1. + fo3(u1,realnu) + fo3(u2,realnu))
            o3ems(i,k1) = dbvtt(i)*h2otr(i,k1)*o3bndi
            trem5(i)    = 1.-(o3bndi/(1060-980.))
  220    continue
c
  100 continue
c
c end of interface loop
c
c-----------------------------------------------------------------------
c
c compute total emissivity:
c
      do 300 k1=1,plevp
         do 310 i=1,plon
            emstot(i,k1) = h2oems(i,k1)
     $                   + co2ems(i,k1) + o3ems(i,k1)
  310    continue
  300 continue
c
c done
c
      return
      end

      subroutine radalb(jlat  ,oro   ,sndpth ,coszrs,
     $                  albs  ,albl  ,albsd  ,albld )
c---------------------------------------------------------------------
c
c compute surface albedos
c
c computes surface albedos for direct/diffuse incident radiation for
c two spectral intervals:
c   s = 0.2-0.7 micro-meters
c   l = 0.7-5.0 micro-meters
c
c uses knowledge of surface type to specify albedo, as follows:
c
c
c ocean           uses solar zenith angle to compute albedo for direct
c                 radiation; diffuse radiation values constant; albedo
c                 independent of spectral interval and other physical
c                 factors such as ocean surface wind speed.
c  
c land without    albedos specified by two dimensional surface albedo 
c     snow        fields, which distinguish surfaces with strong solar
c                 zenith angle dependence from those with weaker solar
c                 zenith angle dependence; albedo independent of surface
c                 moisture or other physical factors.
c
c land with snow  snow depth (liquid water equivalent) used, along with
c                 aerodynamic roughness to define a horizontal fraction
c                 of land surface covered with snow; snow albedos are
c                 computed as functions of solar zenith angle; these snow
c                 albedos are then weighted by the horizontal fraction 
c                 of coverage with the underlying surface albedos 
c                 computed above to produce total grid mean albedo.
c
c land with ice   surface albedos specified as functions of spectral
c                 interval; combined with overlying snow in a similar 
c                 manner to the case of land with snow.
c
c ocean with      surface albedos specified; combined with overlying snow
c  sea  ice       in a similar manner to the case of land with snow.
c
c
c note, the code collects together surfaces of the same type for various
c computations in order to vectorize longitude loops
c
c for more information, see Briegleb, B.P., P.Minnis, V.Ramanathan,
c E.Harrison, 1986: Comparison of Regional Clear-Sky Albedos Inferred from
c Satellite Observations and Model Computations. JCAM Vol. 25, pp 214-226
c
c the details of the land surface albedo arrays can be found in the
c common block description below.
c
c---------------------------------------------------------------------
      implicit none
c---------------------------------------------------------------------
c
c     !   Basic grid point resolution parameters
c
      integer plev,plon,plat,plonp2,plevp,pcnst,plond,platd,plevd
      integer nxpt,jintmx,i1,j1,j1m
c
      parameter(plev=  133,        ! number of vertical levels
     $          plevp=plev+1,     ! number of vertical interfaces
     $          plon=1,           ! number of longitudes (T42)
     $          plat=1,           ! number of latitudes (T42)
     $          plonp2=plon+2,    ! length necessary for fft
     $          pcnst=1,          ! number of constituents
     $          nxpt=1,           !
     $          jintmx=2,         !
     $          plond=plon + 1 + 2*nxpt, ! slt extended domain longitude
     $          platd=plat + 2*nxpt + 2*jintmx, ! slt extended domain lat.
     $          plevd=plev*(3+pcnst), ! fold plev,pcnst indices into one
     $          i1=1+nxpt,        ! model starting longitude (3-d)
     $          j1=1+nxpt+jintmx, ! model starting latitude (3-d)
     $          j1m=j1-1)         ! model starting offset (3-d)
c
c --------------------------------------------------------------------
c
c surface albedo data
c
c vs = 0.2 - 0.7 micro-meters wavelength range
c ni = 0.7 - 5.0 micro-meters wavelength range
c
c s  = strong zenith angle dependent surfaces
c w  = weak   zenith angle dependent surfaces
c
c the albedos are computed for a model grid box by ascribing values to
c high resolution points from a vegetation dataset, then linearlly 
c averaging to the grid box value; ocean and land values are averaged
c together along coastlines; the fraction of every grid box that has
c strong zenith angle dependence is included also.
c
      real albvss, ! grid box alb for vs over strng zn srfs
     $     albvsw, ! grid box alb for vs over weak  zn srfs
     $     albnis, ! grid box alb for ni over strng zn srfs
     $     albniw, ! grid box alb for ni over weak  zn srfs
     $     frctst  ! fraction of area in grid box strng zn 
c
      common /crdalb/albvss(plond,plat),albvsw(plond,plat),
     $               albnis(plond,plat),albniw(plond,plat),
     $               frctst(plond,plat)
c
c
c surface boundary data
c
c vegtyp is used to specify the thermal properites of the surface,
c as well as determine the location of permanent land ice points;
c it is the dominant surface type within the model grid box based
c on a high resolution vegetation type dataset;
c it encoded in the following manner:
c
c   1        ocean
c   2        sea ice
c   3        permanent land ice
c   4        tropical evergreen forest
c   5        deciduous forest
c   6        grassland/tundra
c   7        desert
c
c rghnss is the aerodynamic roughness length for the grid box, computed
c by linear averaging of the values ascribed to high resolution 
c vegetation dataset values; ocean and land values are averaged together
c at coastlines.
c
c evapf is the ratio of actual to potential evaporation, and is computed
c from the high resolution vegetation dataset in a manner similar to the
c aerodynamic roughness.
c
c vevapf allows for variable snow cover, where the underlying
c evaporability factor is modified; see radalb.
c
c snwjan and snwjly are mean climatological snow depths (liquid water
c equivalent) used to compute the actual daily values of snow cover.
c
      real vegtyp, ! surface thermal type, based on veg type
     $     rghnss, ! aerodynamic roughness length
     $     evapf , ! constant surface evaporability
     $     vevapf, ! variable surface evaporability
     $     snwjan, ! snow cover (liq water equiv) for january
     $     snwjly  ! snow cover (liq water equiv) for july
c
      common /crdsrf/vegtyp(plond,plat),rghnss(plond,plat),
     $               evapf (plond,plat),vevapf(plond,plat),
     $               snwjan(plond,plat),snwjly(plond,plat)
c
c --------------------------------------------------------------------
c
c input arguments
c
      integer jlat         ! latitude index for two dimensional data arrays
      real oro(plond),     ! surface type flag (ocean, land, sea ice)
     $     sndpth(plond),  ! snow depth (liquid water equivalent)
     $     coszrs(plond)   ! cosine solar zenith angle
c
c output arguments
c
      real albs(plond),    ! srf alb for direct rad   0.2-0.7 micro-meters
     $     albl(plond),    ! srf alb for direct rad   0.7-5.0 micro-meters
     $     albsd(plond),   ! srf alb for diffuse rad  0.2-0.7 micro-meters
     $     albld(plond)    ! srf alb for diffuse rad  0.7-5.0 micro-meters
c
      external  wheneq     ! when equal function; gives indices for condition
c
c local workspace
c
      integer           i, ! longitude index
     $                 ii, ! lngtd indx for pnts satisfying certain conditons
     $        ipos(plond), ! flag for computation points
     $        indx(plond), ! indices for computation points
     $               npts  ! number of computation points
c
      logical ilf(plond)   ! logical for specifying certain types of surface
c
      real ocean(plond),   ! ocean flag   (=1 for ocean,   0 otherwise)
     $      sice(plond),   ! sea ice flag (=1 for sea ice, 0 otherwise)
     $               rs,   ! empirical factor for strng znth angl dependence
     $               rw,   ! empirical factor for weak  znth angl dependence
     $    frsnow(plond),   ! horizontal fraction of snow cover
     $           snwhgt,   ! physical snow height
     $           rghsnw    ! roughness for horizontal snow cover fraction
c
      real salbs(plond),   ! snow alb for direct rad  0.2-0.7 micro-meters
     $     salbl(plond),   ! snow alb for direct rad  0.7-5.0 micro-meters
     $     salbsd(plond),  ! snow alb for diffuse rad  0.2-0.7 micro-meters
     $     salbld(plond)   ! snow alb for diffuse rad  0.7-5.0 micro-meters
c
c albedos for snow, land ice, and sea ice:
c
      real snws,    ! snow albedo for 0.2-0.7 micro-meters
     $     snwl,    ! snow albedo for 0.7-5.0 micro-meters
     $     sices,   ! sea ice albedo for 0.2-0.7 micro-meters
     $     sicel    ! sea ice albedo for 0.7-5.0 micro-meters
c
      data snws    / .95  /
      data snwl    / .70  /
      data sices   / .70  /
      data sicel   / .50  /
c
c-----------------------------------------------------------------------
c
c set flags for ocean (with or without sea ice) and 
c ocean with sea ice points:
c
CDIR$ IVDEP
      do 10 i=1,plon
         ocean(i) = 0.
         sice(i)  = 0.
         if(int(oro(i)+.45).eq.0) ocean(i) = 1.
         if(int(oro(i)+.45).eq.2) sice(i)  = 1.
   10 continue
c
c initialize all surface albedos to zero
c
      do 20 i=1,plon
         albs(i)  = 0.
         albl(i)  = 0.
         albsd(i) = 0.
         albld(i) = 0.
   20 continue
c
c land surfaces 
c
      do 50 i=1,plon
c
         if(ocean(i).ne.1.0 .and. sice(i).eq.0.0 .and.
     $      coszrs(i).gt.0.0) then
              ipos(i) = 1
         else
              ipos(i) = 0
         endif
c
   50 continue
      call wheneq(plon,ipos,1,1,indx,npts)
CDIR$ IVDEP
      do 60 ii=1,npts
         i = indx(ii)
c
c use empirical factors to adjust surface albedos for zenith angle
c effects, distinquishing between strong and weakly dependent surfaces:
c
         rs = 1.4/(1. + .8*coszrs(i))
         rw = 1.1/(1. + .2*coszrs(i))
         albs(i)  = albvss(i,jlat)*frctst(i,jlat)*rs +
     $              albvsw(i,jlat)*(1. - frctst(i,jlat))*rw
         albl(i)  = albnis(i,jlat)*frctst(i,jlat)*rs +
     $              albniw(i,jlat)*(1. - frctst(i,jlat))*rw
         albsd(i) = albvss(i,jlat)*frctst(i,jlat) +
     $              albvsw(i,jlat)*(1. - frctst(i,jlat))
         albld(i) = albnis(i,jlat)*frctst(i,jlat) +
     $              albniw(i,jlat)*(1. - frctst(i,jlat))
   60 continue
c
c surfaces covered with sea ice:
c
      do 70 i=1,plon
         if( (sice(i).eq.1.0 .and. coszrs(i).gt.0.0) ) then
           ipos(i) = 1
         else
           ipos(i) = 0
         endif
   70 continue
c
      call wheneq(plon,ipos,1,1,indx,npts)
CDIR$ IVDEP
      do 80 ii=1,npts
         i = indx (ii)
         albs(i)  = sices
         albl(i)  = sicel
         albsd(i) = albs(i)
         albld(i) = albl(i)
   80 continue
c
c points with snowcover:
c
      do 90 i=1,plon
         ilf(i) = sndpth(i) .gt. 0. .and. coszrs(i) .gt. 0. 
         if( sndpth(i).gt.0. .and. coszrs(i).gt.0. ) then
           ipos(i) = 1
         else
           ipos(i) = 0
         endif
   90 continue
c
      call wheneq(plon,ipos,1,1,indx,npts)
      do 100 ii=1,npts
         i = indx(ii)
         salbsd(i) = snws
         salbld(i) = snwl
  100 continue
c
c direct snow albedos: distinguish between 2 zenith angle regimes
c
      do 110 i=1,plon
         ipos(i) = 0
         if(ilf(i)) ipos(i) = 2
         if(ilf(i) .and. coszrs(i).lt.0.5) ipos(i) = 1
  110 continue
c
c zenith angle regime 1 ( coszrs < 0.5 ).
c set direct snow albedos (limit to 0.98 max)
c
      call wheneq(plon,ipos,1,1,indx,npts)
      do 120 ii=1,npts
         i        = indx(ii)
         salbs(i) = amin1(0.98,salbsd(i) + (1. - salbsd(i))*0.5*
     $                    ((3./(1. + 4.*coszrs(i))) - 1.))
         salbl(i) = amin1(0.98,salbld(i) + (1. - salbld(i))*0.5*
     $                    ((3./(1. + 4.*coszrs(i))) - 1.))
  120 continue
c
c zenith angle regime 2 ( coszrs >= 0.5 )
c
      call wheneq(plon,ipos,1,2,indx,npts)
      do 130 ii=1,npts
         i        = indx(ii)
         salbs(i) = snws
         salbl(i) = snwl
  130 continue
c
c compute fraction of horizonal surface covered with snow:
c
      do 140 i=1,plon
         if(sndpth(i) .le. 0.) then
            snwhgt    = 0.0
            rghsnw    = 0.0
            frsnow(i) = 0.0
         else
            snwhgt    = 20. * sndpth(i)
            rghsnw    = amax1(rghnss(i,jlat),0.25)
            frsnow(i) = snwhgt/(rghsnw + snwhgt)
         endif
c
  140 continue
c
c points with snowcover:
c
      do 145 i=1,plon
         if( sndpth(i).gt.0. .and. coszrs(i).gt.0. ) then
           ipos(i) = 1
         else
           ipos(i) = 0
         endif
  145 continue
      call wheneq(plon,ipos,1,1,indx,npts)
c
c compute both diffuse and direct total albedos 
c
CDIR$ IVDEP
      do 150 ii=1,npts
         i = indx(ii)
         albs(i)  = albs(i) *(1.-frsnow(i)) + salbs(i) *frsnow(i)
         albl(i)  = albl(i) *(1.-frsnow(i)) + salbl(i) *frsnow(i)
         albsd(i) = albsd(i)*(1.-frsnow(i)) + salbsd(i)*frsnow(i)
         albld(i) = albld(i)*(1.-frsnow(i)) + salbld(i)*frsnow(i)
  150 continue
c
c local points over ice-free ocean with incident solar radiation:
c
      do 160 i=1,plon
         ipos(i)=0
         if(ocean(i).eq.1. .and. sice(i).eq.0. .and. coszrs(i).ge.0.)
     $     ipos(i) = 1
  160 continue
      call wheneq(plon,ipos,1,1,indx,npts)
c
c ocean albedos function of solar zenith angle only:
c
CDIR$ IVDEP
      do 170 ii=1,npts
         i = indx(ii)
         albl(i)  = (.026/(coszrs(i)**1.7 + .065)) +
     $           (.15*(coszrs(i)-.10)*(coszrs(i)-.50)*(coszrs(i)-1.))
         albs(i)  = albl(i)
         albld(i) = 0.06
         albsd(i) = 0.06
  170 continue
c
c done
c
c      do 1717 i=1,plon
c      write(6,1718) i,albs(i),albsd(i),albl(i),albld(i)
c1718  format(' i albs albsd albl albld = ',i2,1x,4(f7.5,1x))
c1717 continue
c
c
      return
      end
      subroutine radcsw(jlat    ,oro     ,sndpth  ,pint    ,h2ommr  ,
     $                  cld     ,clwp    ,o3mmr   ,eccf    ,coszrs  ,
     $                  solin   ,qrs     ,frsa    ,sabtp   ,clrss   ,
     $                  clrst   ,
     &     fupsw,fdsw,
     &     geff,aeff,taueff,crecz,swfdiv,
     &     bccmil,tccmil,numlay,ncclev,
     &     cbi,cti,ncldlay,
     &     idebug)
c-----------------------------------------------------------------------
c solar radiation code
c
c computes incident solar flux, solar heating rate, surface absorbed
c solar flux, and total column absorbed solar flux
c
c uses the delta-eddington method
c
c divides solar spectrum into 18 intervals from 0.2-5.0 micro-meters.
c solar flux fractions specified for each interval. allows for seasonally
c and diurnally varying solar input.  includes molecular, cloud, and
c surface scattering, along with h2o,o3,co2,o2,cloud, and surface 
c absorption. computes delta-eddington reflections and transmissions
c assuming homogeneously mixed layers.  computes surface albedos by 
c invoking radalb. adds the layers assuming scattering between layers to 
c be isotropic, and distinguishes direct solar beam from scattered 
c radiation down to the surface.
c
c longitude loops are broken into 1 or 2 sections, so that only daylight
c (i.e. coszrs > 0) computations are done.
c
c note that an extrae layer above the model top layer is added.
c
c cgs units are used.
c
c special diagnostic calculation of the clear sky surface and total column
c absorbed flux is also done; this calculation does not effect the rest
c of the model, but is included for cloud forcing diagnostics.
c
c for more details, see  briegleb, b.p.   delta-eddington approximation
c for solar radiation in the ncar community climate model, submitted
c to journal of geophysical research
c
c-----------------------------------------------------------------------
      implicit none
c-----------------------------------------------------------------------
c
c     !   Basic grid point resolution parameters
c
      integer plev,plon,plat,plonp2,plevp,pcnst,plond,platd,plevd
      integer nxpt,jintmx,i1,j1,j1m
c
      parameter(plev=  133,        ! number of vertical levels
     $          plevp=plev+1,     ! number of vertical interfaces
     $          plon=1,           ! number of longitudes (T42)
     $          plat=1,           ! number of latitudes (T42)
     $          plonp2=plon+2,    ! length necessary for fft
     $          pcnst=1,          ! number of constituents
     $          nxpt=1,           !
     $          jintmx=2,         !
     $          plond=plon + 1 + 2*nxpt, ! slt extended domain longitude
     $          platd=plat + 2*nxpt + 2*jintmx, ! slt extended domain lat.
     $          plevd=plev*(3+pcnst), ! fold plev,pcnst indices into one
     $          i1=1+nxpt,        ! model starting longitude (3-d)
     $          j1=1+nxpt+jintmx, ! model starting latitude (3-d)
     $          j1m=j1-1)         ! model starting offset (3-d)
c
c-----------------------------------------------------------------------
c
c radiation constants
c
      real gravit,   ! gravitational acceleration
     $        rga,   ! 1 over gravit
     $      cpair,   ! heat capacity air at constant pressure
     $     epsilo,   ! ratio mmw h2o to mmw air
     $       sslp,   ! standard pressure
     $     stebol,   ! stephan boltzmann constant
     $     rgsslp,   ! 0.5 / (gravit*sslp)
     $     co2vmr,   ! co2 volume mixing ratio
     $      dpfo3,   ! Doppler factor for o3
     $     dpfco2,   ! Doppler factor for co2
     $     dayspy,   ! solar days in one year
     $        pie    ! pie
c
      common/crdcon/gravit,    rga,  cpair,  epsilo,   sslp,
     $              stebol, rgsslp, co2vmr,   dpfo3, dpfco2,
     $              dayspy,    pie
c
c-----------------------------------------------------------------------
c
c input arguments
c
      integer jlat               ! latitude index used by 'radalb'
      real oro(plond),           ! land/ocean/seaice flag used by 'radalb'
     $     sndpth(plond),        ! snow depth used by 'radalb'
     $     pint(plond,plevp),    ! interface pressure 
     $     h2ommr(plond,plevp),  ! specific humidity (h2o mass mixing ratio)
     $     cld(plond,plevp),     ! fractional cloud cover
     $     clwp(plond,plev),     ! layer liquid water path
     $     o3mmr(plond,plev),    ! ozone mass mixing ratio
     $     eccf,                 ! ecctrcty fctr (earth-sun inv dstnc sqrd)
     $     coszrs(plond)         ! cosine solar zenith angle
c
c output arguments
c
      real solin(plond),     ! incident solar flux
     $     qrs(plond,plev),  ! solar heating rate
     $     frsa(plond),      ! surface absorbed solar flux
     $     sabtp(plond),     ! total column absorbed solar flux
     $     clrss(plond),     ! clear sky surface absorbed solar flux
     $     clrst(plond)      ! clear sky total column absorbed solar flux
c
c externals
c
      integer   isrchfgt,    ! search for first array element > 0
     $          isrchfle     ! search for first array element < 0
      external    radalb,    ! computes surface albedos
     $            radded,    ! computes delta-eddington solution
     $            radclr,    ! computes clear sky delta-edd solution
     $          isrchfgt,    ! search for first array element > 0
     $          isrchfle     ! search for first array element < 0
c
c local arrays and constants
c
      integer       ns,  ! spectral loop index
     $               i,  ! longitude loop index
     $               k,  ! level loop index
     $               n,  ! loop index for daylight
     $           nloop,  ! number of daylight loops
     $             isn,  ! starting daylight longitude index
     $             ien,  ! ending   daylight longitude index
     $           is(2),  ! daytime start indices
     $           ie(2),  ! daytime end indices
     $          indxsl   ! index for cloud particle properties
c
      real scon    ! solar constant
      data scon / 1.370e6 /
c
c a. slingo's data for cloud particle radiative properties
c (from 'a gcm parameterization for the shortwave properties of
c water clouds' jas vol. 46 may 1989  pp 1419-1427)
c
      real abar(4), ! a coefficient for extinction optical depth
     $     bbar(4), ! b coefficiant for extinction optical depth
     $     cbar(4), ! c coefficiant for single particle scat albedo
     $     dbar(4), ! d coefficiant for single particle scat albedo
     $     ebar(4), ! e coefficiant for asymmetry parameter
     $     fbar(4)  ! f coefficiant for asymmetry parameter
c
      data abar/ 2.817e-02, 2.682e-02,2.264e-02,1.281e-02/
      data bbar/ 1.305    , 1.346    ,1.454    ,1.641    /
      data cbar/-5.62e-08 ,-6.94e-06 ,4.64e-04 ,0.201    /
      data dbar/ 1.63e-07 , 2.35e-05 ,1.24e-03 ,7.56e-03 /
      data ebar/ 0.829    , 0.794    ,0.754    ,0.826    /
      data fbar/ 2.482e-03, 4.226e-03,6.560e-03,4.353e-03/
c
      real abari,   ! a coefficiant for current spectral interval
     $     bbari,   ! b coefficiant for current spectral interval
     $     cbari,   ! c coefficiant for current spectral interval
     $     dbari,   ! d coefficiant for current spectral interval
     $     ebari,   ! e coefficiant for current spectral interval
     $     fbari    ! f coefficiant for current spectral interval
c
c caution... a. slingo recommends no less than 4.0 micro-meters
c nor greater than 20 micro-meters
c
      real cldefr      !  universal cloud effective radius in micro-meters
      data cldefr / 10.0 /
c
      real delta  ! pressure (atmospheres) for stratospheric h2o limit:
      data delta  /  1.70e-3 /
c
      real o2mmr  ! o2 mass mixing ratio:
      data o2mmr / .23143 /
c
c co2 info:
c
      real mmwair,   ! mean molecular weight of air
     $     mmwco2,   ! mean molecular weight of co2
     $     co2mmr    ! co2 mass mixing ratio
      data mmwair / 28.9644 /
      data mmwco2 / 44.0000 /
c
      real albdir(plond),    ! current spc intrvl srf alb to direct rad
     $     albdif(plond)     ! current spc intrvl srf alb to diffuse rad
c
      real albs(plond),      ! 0.2-0.7 micro-meter srfc albedo to direct rad
     $     albl(plond),      ! 0.7-5.0 micro-meter srfc albedo to direct rad
     $     albsd(plond),     ! 0.2-0.7 micro-meter srfc albedo to diffuse rad
     $     albld(plond)      ! 0.7-5.0 micro-meter srfc albedo to diffuse rad
c
      integer nspint  ! number of spectral intervals across solar spectrum
      parameter ( nspint = 18 )
c
c next series depends on spectral interval
c
      real frcsol(nspint),      ! fraction of solar flux in each spec int
     $     wavmin(nspint),      ! min wavelength (micro-meters) of interval
     $     wavmax(nspint),      ! max wavelength (micro-meters) of interval
     $     raytau(nspint),      ! rayleigh scattering optical depth
     $     abh2o(nspint),       ! absorption coefficiant for h2o (cm2/g)
     $     abo3 (nspint),       ! absorption coefficiant for o3    "
     $     abco2(nspint),       ! absorption coefficiant for co2   "
     $     abo2 (nspint),       ! absorption coefficiant for o2    "
     $     ph2o(nspint),        ! weight of h2o in spectral interval
     $     po3 (nspint),        ! weight of o3  in spectral interval
     $     pco2(nspint),        ! weight of co2 in spectral interval
     $     po2 (nspint)         ! weight of o2  in spectral interval
c
c
c CAUTION! CAUTION! CAUTION! CAUTION! CAUTION! CAUTION! CAUTION! CAUTION!
c
c The most radical change i made to this code was switching briegleb's
c bands 14 and 15 (last two H2O bands) to be my bands 15 and 16. Therefore
c his CO2 band 16 is my band 14. all the data in the spectral tables has
c been rearranged accordingly, and this makes taking my spectral data
c and plugging it into unmodified CCM2 code a big mistake!!.
c
c
      data frcsol / .001488, .001389, .001290, .001686, .002877,
     $              .003869, .026336, .426131, .526861, .526861,
     $              .526861, .526861, .526861,
     $              .006239, .526861, .526861, .001834, .001834/
c
      data wavmin / .200,  .245,  .265,  .275,  .285,
     $              .295,  .305,  .350,  .700,  .701,
     $              .701,  .701,  .701, 
     $             2.630,  .702,  .702, 4.160, 4.160/
c
      data wavmax / .245,  .265,  .275,  .285,  .295,
     $              .305,  .350,  .700, 5.000, 5.000,
     $             5.000, 5.000, 5.000, 
     $             2.860, 5.000, 5.000, 4.550, 4.550/
c
      data raytau / 4.020, 2.180, 1.700, 1.450, 1.250,
     $              1.085, 0.730, 0.135, 0.020, .0001,
     $              .0001, .0001, .0001, 
     $              .0001, .0001, .0001, .0001, .0001/
c
c absorption coefficiants
c
      data abh2o /    .000,     .000,    .000,    .000,    .000,
     $                .000,     .000,    .000,    .002,    .035,
     $                .377,    1.950,   9.400,  
     $                .000,   44.600, 190.000,    .000,    .000/
c
      data abo3  /
     $ 5.370e+04, 13.080e+04,  9.292e+04, 4.530e+04, 1.616e+04,
     $ 4.441e+03,  1.775e+02,  2.101e+01,      .000,      .000,
     $  .000    ,   .000    ,   .000    ,     
     $  .000    ,    .000,      .000,      .000    ,   .000    /
c
      data abco2  /    .000,     .000,    .000,    .000,    .000,
     $                 .000,     .000,    .000,    .000,    .000,
     $                 .000,     .000,    .000,  
     $                 .094,     .000,    .000,  .196,   1.963/
c
      data abo2  /    .000,     .000,    .000,    .000,    .000,
     $                .000,     .000,1.11e-05,6.69e-05,    .000,
     $                .000,     .000,    .000,  
     $                .000,    .000,    .000,   .000,    .000/
c
c spectral interval weights
c
      data ph2o  /    .000,     .000,    .000,    .000,    .000,
     $                .000,     .000,    .000,    .505,    .210,
     $                .120,     .070,    .048, 
     $                .000,     .029,    .018,  .000,    .000/
      data po3   /   1.000,    1.000,   1.000,   1.000,   1.000,
     $               1.000,    1.000,   1.000,    .000,    .000,
     $                .000,     .000,    .000,  
     $                .000,     .000,    .000,   .000,    .000/
      data pco2  /    .000,     .000,    .000,    .000,    .000,
     $                .000,     .000,    .000,    .000,    .000,
     $                .000,     .000,    .000, 
     $               1.000,     .000,    .000,  .640,    .360/
      data po2   /    .000,     .000,    .000,    .000,    .000,
     $                .000,     .000,   1.000,   1.000,    .000,
     $                .000,     .000,    .000,  
     $                .000,     .000,    .000,   .000,    .000/
c
c diagnostic and accumulation arrays; note that sfltot, fswup, and
c fswdn are not used in the computation, but are retained for future use.
c
      real solflx(plond),         ! solar flux in current interval
     $     sfltot(plond),         ! spectrally summed total solar flux
     $     totfld(plond,0:plev),  ! spectrally summed flux divergence
     $     fswup(plond,0:plevp),  ! spectrally summed up flux
     $     fswdn(plond,0:plevp)   ! spectrally summed down flux
c
c cloud radiative property arrays
c
      real tauexc(plond,0:plev),     ! cloud extinction optical depth
     $         wc(plond,0:plev),     ! cloud single scattering albedo
     $         gc(plond,0:plev),     ! cloud assymetry parameter
     $         fc(plond,0:plev),     ! cloud forward scattered fraction
     $         cre(plev),            ! cloud effective radius (micro-meters)
     $         rcre(plev)            ! inverse of cloud effective radius
c
c various arrays and other constants:
c
      real pflx(plond,0:plevp), ! interface pressure, including extrae layer
     $     zenfac(plond),       ! square root of cosine solar zenith angle
     $     sqrco2,              ! square root of the co2 mass mixing ratio
     $     tmp1,                ! temporary constant array
     $     tmp2,                ! temporary constant array
     $     tmp3,                ! temporary constant array
     $     pdel,                ! pressure difference across layer
     $     path,                ! mass path of layer
     $     ptop,                ! lower interface pressure of extrae layer
     $     ptho2,               ! used to compute mass path of o2
     $     ptho3,               ! used to compute mass path of o3
     $     pthco2,              ! used to compute mass path of co2
     $     pthh2o,              ! used to compute mass path of h2o
     $     h2ostr,              ! inverse square root h2o mass mixing ratio
     $     wavmid,              ! spectral interval middle wavelength
     $     trayoslp             ! rayleigh opt dpth over standard pressure
c
      real rdenom,              ! multiple scattering term
     $     psf,                 ! fraction of solar flux in spectral interval
     $     gocp                 ! grav acc over heat capacity constant prs
c
c layer absorber amounts; note that 0 refers to the extrae layer 
c added above the top model layer
c
      real uh2o(plond,0:plev), ! layer absorber amount of h2o
     $      uo3(plond,0:plev), ! layer absorber amount of  o3
     $     uco2(plond,0:plev), ! layer absorber amount of co2
     $      uo2(plond,0:plev)  ! layer absorber amount of  o2
c
c total column absorber amounts:
c
      real uth2o(plond),  ! total column  absorber amount of  h2o
     $     uto3(plond),   ! total column  absorber amount of  o3
     $     utco2(plond),  ! total column  absorber amount of  co2
     $     uto2(plond)    ! total column  absorber amount of  o2
c
c these arrays are defined for plev model layers; 0 refers to the
c extrae layer on top:
c
      real rdir(plond,0:plev),   ! layer reflectivity to direct radiation
     $     rdif(plond,0:plev),   ! layer reflectivity to diffuse radiation
     $     tdir(plond,0:plev),   ! layer transmission to direct radiation
     $     tdif(plond,0:plev),   ! layer transmission to diffuse radiation
     $     explay(plond,0:plev), ! solar beam exp transmission for layer
     $     flxdiv(plond,0:plev)  ! flux divergence for layer
c
c these arrays are defined at model interfaces; 0 is the top of the
c extrae layer above the model top; plevp is the earth surface:
c
      real rupdir(plond,0:plevp), ! ref to dir radiation for layers below
     $     rupdif(plond,0:plevp), ! ref to dif radiation for layers below
     $     rdndif(plond,0:plevp), ! ref to dif radiation for layers above
     $     exptdn(plond,0:plevp), ! solar beam exp down transmission from top
     $     tottrn(plond,0:plevp), ! total transmission for layers above
     $     fluxup(plond,0:plevp), ! up   flux at model interface
     $     fluxdn(plond,0:plevp)  ! down flux at model interface
c
c------------------------------------------------------------------------
c     zender cirrus cloud hack workspace
c
c CAUTION! CAUTION! CAUTION! CAUTION! CAUTION! CAUTION! CAUTION! CAUTION!
c
c The most radical change i made to this code was switching briegleb's
c bands 14 and 15 (last two H2O bands) to be my bands 15 and 16. Therefore
c his CO2 band 16 is my band 14. all the data in the spectral tables has
c been rearranged accordingly, and this makes taking my spectral data
c and plugging it into unmodified CCM2 code a big mistake!!.
c
c Here is a list of the band # and the corresponding Slingo and Ebert-Curry 
c categories:
c  my band  slingo   E&C   CCM2 band     
c     1       1       1        1  
c     2       1       1        2
c     3       1       1        3 
c     4       1       1        4 
c     5       1       1        5 
c     6       1       1        6 
c     7       1       1        7 
c     8       1       1        8  
c     9       2       2        9 
c     10      3       3        10 
c     11      3       3        11 
c     12      3       4        12 
c     13      3       4        13 
c     14      4       5        16 
c     15      4       5        14 
c     16      4       5        15 
c     17      4       5        17 
c     18      4       5        18 
c     (clearly i'm equivalencing slingo bin 4 with E&C bin 5)      
c     
      real cgsmks     ! conversion factor for fluxes from cgs to mks
c
      data cgsmks / 1.e-3 /
c
      integer 
     &     idebug,    !debugging level from cloud model
     &     indxec,    ! index for cloud Ebert & Curry particle properties
     &     kcz,       !index for index trickery
     &     bccmil,    !bottom ccm2 interface level
     &     tccmil,    !top ccm2 interface level
     &     numlay,    !num_layer in cloud model
     &     ncclev,    !num_ccm2_level (with overlaps already subtracted)
     &     cbi,       !cloud bottom index
     &     cti,       !cloud top index
     &     ncldlay    !number of layers in cloud itself
C
      real
     &     geff(0:nspint+1,0:numlay+1),   !effective asymmetry factor
     &     aeff(0:nspint+1,0:numlay+1),   !effective singel scat albedo
     &     taueff(0:nspint+1,0:numlay+1), !effective extinction optical depth
     &     swfdiv(0:nspint+1,0:numlay+1),  !SW spectral flux divergence
     &     crecz(0:numlay+1), !effective radius of particles (meters)
     &     fupsw(0:numlay+1), !total shortwave upward directed flux
     &     fdsw(0:numlay+1) !total shortwave downward directed flux
c
c ebert & curry's data for cloud particle radiative properties
c (from 'a parameterization of ice cloud optical properties for
c climate models' jgr vol. 97  no. D4 march 20, 1992  pp 3831-3836)
c the 'e' suffix stands for "Ebert" and the parameters are 
c otherwise analogous to their liquid water slingo counterparts
c
      real 
     $     tmp1ec,                ! temporary constant array
     $     tmp2ec,                ! temporary constant array
     $     tmp3ec                 ! temporary constant array
c
      real ecwavmin(nspint) ! min wavelength (micro-meters) of interval
c
      data ecwavmin / .200,  .245,  .265,  .275,  .285,
     $              .295,  .305,  .350,  .700,  .701,
     $              .701,  .702,  .702,
     $             2.630,   .703,  .703, 4.160, 4.160/
c
      real abare(5), ! a coefficient for extinction optical depth
     $     bbare(5), ! b coefficiant for extinction optical depth
     $     cbare(5), ! c coefficiant for single particle scat albedo
     $     dbare(5), ! d coefficiant for single particle scat albedo
     $     ebare(5), ! e coefficiant for asymmetry parameter
     $     fbare(5)  ! f coefficiant for asymmetry parameter
c
      data abare/ 3.448e-03, 3.448e-03, 3.448e-03, 3.448e-03, 3.448e-03/
      data bbare/ 2.431    , 2.431    , 2.431    , 2.431    , 2.431    /
      data cbare/ .00001   , .00011   , .01240   , .03779   , .46658   /
      data dbare/ .0000    , 1.405e-05, 6.867e-04, 1.284e-03, 2.050e-05/
      data ebare/ .7661    , .7730    , .7865    , .8172    , .9595    /
      data fbare/ 5.851e-04, 5.665e-04, 7.204e-04, 7.463e-04, 1.076e-04/
c
      real abarei,   ! a coefficiant for current spectral interval
     $     bbarei,   ! b coefficiant for current spectral interval
     $     cbarei,   ! c coefficiant for current spectral interval
     $     dbarei,   ! d coefficiant for current spectral interval
     $     ebarei,   ! e coefficiant for current spectral interval
     $     fbarei    ! f coefficiant for current spectral interval
c
c caution... a. slingo recommends no less than 4.0 micro-meters
c nor greater than 20 micro-meters
c
c     end of zender workspace alterations
c-----------------------------------------------------------------------
c
c compute surface albedos:
c
      call radalb(jlat    ,oro     ,sndpth  ,coszrs  ,
     $            albs    ,albl    ,albsd   ,albld   )
c
c set cloud effective radius in micrometers:
c
cz      do 8 k=1,plev
cz         cre(k)  = cldefr
cz         rcre(k) = 1./cre(k)
cz    8 continue
c recall the effective radius isn't used unless:
c
c debug 47 is use the Ebert-Curry scheme
c debug 53 is use the Slingo scheme
c 
      if((idebug.eq.47).or.(idebug.eq.53)) then
         do 11 k=1,plev
            if((k.ge.(plev+1-bccmil-numlay)).and.
     &           (k.le.(plev+1-bccmil-1))) then
c     cloud grid levels of Ebert and Curry and Slingo schemes
               kcz=plev-k-bccmil+1
c     remember to convert the meters to microns
               cre(k) = crecz(kcz)*1.e6
            else
c     levels outside cloud grid in Ebert-Curry and Slingo schemes
               cre(k) = cldefr
            endif
            rcre(k) = 1./cre(k)
 11      continue
      endif
c     
c     initialize output fields:
c     
      do 13 i=1, plon
         sabtp(i) = 0.0
         frsa(i)  = 0.0
         solin(i) = 0.0
         clrss(i) = 0.0
         clrst(i) = 0.0
 13   continue
      do 15 k=1, plev
         do 20 i=1, plon
            qrs(i,k) = 0.0
   20    continue
   15 continue
c
c compute starting, ending daytime loop indices:
c
      nloop = 0
      is(1) = isrchfgt(plon,coszrs,1,0.0)
c
c if night everywhere, return:
c
      if(is(1).gt.plon) then
         print *,"WARNING: It's nightie night time in the SW model!"
         return            
      endif
      ie(1) = isrchfle(plon-is(1),coszrs(is(1)+1),1,0.0) + is(1) - 1
      nloop = 1
c
c possibly 2 daytime loops needed:
c
      if(ie(1).ne.plon) then              
         is(2) = isrchfgt(plon-ie(1),coszrs(ie(1)+1),1,0.0) + ie(1)
         if(is(2).le.plon) then
            nloop = 2
            ie(2) = plon
         end if
      end if
c
c define solar incident radiation and interface pressures:
c
      do 25 n=1,nloop
        do 26 i=is(n),ie(n)
           solin(i) = scon*eccf*coszrs(i)
           pflx(i,0) = 0.
   26   continue
   25 continue
      do 33 k=1,plevp
         do 34 n=1,nloop
            do 36 i=is(n),ie(n)
               pflx(i,k) = pint(i,k)
   36       continue
   34    continue
   33 continue
c
c compute optical paths:
c
      tmp1   = 0.5/(gravit*sslp)
      co2mmr = co2vmr*(mmwco2/mmwair)
      sqrco2 = sqrt(co2mmr)
      do 38 n=1,nloop
         do 40 i=is(n),ie(n)
            ptop      = pflx(i,1)
            ptho2     = o2mmr * ptop / gravit
            ptho3     = o3mmr(i,1) * ptop / gravit
            pthco2    = sqrco2 * (ptop / gravit)
c
cz the next line is why you can't set h2ommr to exactly zero
            h2ostr    = sqrt( 1. / h2ommr(i,1) )
            zenfac(i) = sqrt(coszrs(i))
            pthh2o    = ptop**2*tmp1 +
     $                (ptop*rga)*(h2ostr*zenfac(i)*delta)
c
            uh2o(i,0) = h2ommr(i,1)*pthh2o
            uco2(i,0) = zenfac(i)*pthco2
            uo2 (i,0) = zenfac(i)*ptho2
            uo3 (i,0) = ptho3
   40    continue
   38 continue
c
      tmp2 = delta/gravit
      do 43 k=1,plev
        do 44 n=1,nloop
          do 46 i=is(n),ie(n)
            pdel   = pflx(i,k+1) - pflx(i,k)
            path   = pdel / gravit
            ptho2  = o2mmr * path
            ptho3  = o3mmr(i,k) * path
            pthco2 = sqrco2 * path
c
            h2ostr = sqrt(1.0/h2ommr(i,k))
            pthh2o = (pflx(i,k+1)**2 - pflx(i,k)**2)*tmp1
     $          + pdel*h2ostr*zenfac(i)*tmp2
c
            uh2o(i,k) = h2ommr(i,k)*pthh2o
            uco2(i,k) = zenfac(i)*pthco2
            uo2 (i,k) = zenfac(i)*ptho2
            uo3 (i,k) = ptho3
   46     continue
   44   continue
   43 continue
c
c compute column absorber amounts for the clear sky computation:
c
      do 47 n=1,nloop
         do 48 i=is(n),ie(n)
            uth2o(i) = 0.0
            uto3(i)  = 0.0
            utco2(i) = 0.0
            uto2(i)  = 0.0
            do 49 k=0,plev
               uth2o(i) = uth2o(i) + uh2o(i,k)
               uto3(i)  = uto3(i)  + uo3(i,k)
               utco2(i) = utco2(i) + uco2(i,k)
               uto2(i)  = uto2(i)  + uo2(i,k)
   49       continue
   48    continue
   47 continue
c
c initialize spectrally integrated totals:
c
      do 50 k=0,plev
         do 52 i=1,plon
            totfld(i,k) = 0.0
            fswup (i,k) = 0.0
            fswdn (i,k) = 0.0
   52    continue
   50 continue
      do 58 i=1,plon
         sfltot(i)       = 0.0
         fswup (i,plevp) = 0.0
         fswdn (i,plevp) = 0.0
   58 continue
c
c set cloud properties for top (0) layer; so long as tauexc is zero, 
c there is no cloud above top of model; the other cloud properties
c are arbitrary:
c
      do 206 n=1,nloop
         do 208 i=is(n),ie(n)
            tauexc(i,0) = 0.
            wc(i,0)     = 0.999999
            gc(i,0)     = 0.85
            fc(i,0)     = 0.725
  208    continue
  206 continue
c
c begin spectral loop
c
      do 100 ns=1,nspint
c
c set index for cloud particle properties based on the wavelength,
c according to a. slingo (1989) equations 1-3:
c use index 1 (0.25 to 0.69 micrometers) for visible
c use index 2 (0.69 - 1.19 micrometers) for near-infrared
c use index 3 (1.19 to 2.38 micrometers) for near-infrared
c use index 4 (2.38 to 4.00 micrometers) for near-infrared
c
c for Ebert & Curry cirrus param 
c use index 1 (0.25 to 0.7 micrometers) for visible
c use index 2 (0.7 - 1.3 micrometers) for near-infrared
c use index 3 (1.3 to 1.9 micrometers) for near-infrared
c use index 4 (1.9 to 2.5 micrometers) for near-infrared
c use index 5 (2.5 to 3.5 micrometers) for near-infrared
c and for now just extend index 5 params out to 4 microns
c
c note that the minimum wavelength is encoded (with .001, .002, .003)
c in order to specify the index appropriate for the near-infrared
c cloud absorption properties
c
        if(wavmax(ns) .le. 0.7) then
           indxsl = 1
        else if(wavmin(ns) .eq. 0.700) then
           indxsl = 2
        else if(wavmin(ns) .eq. 0.701) then
           indxsl = 3
        else if(wavmin(ns) .eq. 0.702 .or. wavmin(ns) .gt. 2.38) then
           indxsl = 4
        endif
c
        if(wavmax(ns) .le. 0.7) then
           indxec = 1
        else if(ecwavmin(ns) .eq. 0.700) then
           indxec = 2
        else if(ecwavmin(ns) .eq. 0.701) then
           indxec = 3
        else if(ecwavmin(ns) .eq. 0.702) then
           indxec = 4
        else if(ecwavmin(ns) .eq. 0.703 .or. ecwavmin(ns) .gt.2.5) then
           indxec = 5
        endif
c
c set cloud extinction optical depth, single scatter albedo,
c asymmetry parameter, and forward scattered fraction:
c 
        abari = abar(indxsl)
        bbari = bbar(indxsl)
        cbari = cbar(indxsl)
        dbari = dbar(indxsl)
        ebari = ebar(indxsl)
        fbari = fbar(indxsl)
c
        abarei = abare(indxec)
        bbarei = bbare(indxec)
        cbarei = cbare(indxec)
        dbarei = dbare(indxec)
        ebarei = ebare(indxec)
        fbarei = fbare(indxec)
c
        do 228 k=1,plev
           tmp1 = abari + bbari*rcre(k)
           tmp2 = 1. - cbari - dbari*cre(k)
           tmp3 = fbari*cre(k)
           tmp1ec = abarei + bbarei*rcre(k)
           tmp2ec = 1. - cbarei - dbarei*cre(k)
           tmp3ec = fbarei*cre(k)
           do 229 n=1,nloop
              do 230 i=is(n),ie(n)
c
c do not let single scatter albedo be 1; delta-eddington solution
c for non-conservative case:
c
c debug 47 is use the Ebert Curry scheme
c debug 53 is use the Slingo scheme
c
                 if(idebug.eq.53) then
c     all levels of Slingo scheme  
                    tauexc(i,k) = clwp(i,k)*tmp1
                    wc(i,k)     = amin1(tmp2,.999999)
                    gc(i,k)     = ebari + tmp3
                 else if((k.ge.(plev+1-bccmil-numlay)).and.
     &                   (k.le.(plev+1-bccmil-1)).and.
     &                   idebug.ne.47) then
c     cloud grid levels of complete cirrus cloud Mie/microphysical scheme
                    kcz=plev-k-bccmil+1
                    tauexc(i,k) = taueff(ns,kcz)
                    wc(i,k) = amin1(aeff(ns,kcz),.999999)
                    gc(i,k) = geff(ns,kcz)
                 else
c     levels outside cloud grid in full cirrus cloud scheme
c     and entire grid with Ebert-Curry scheme
                    tauexc(i,k) = clwp(i,k)*tmp1ec
                    wc(i,k)     = amin1(tmp2ec,.999999)
                    gc(i,k)     = ebarei + tmp3ec
                 endif
c     
c cloud fraction incorporated into cloud extinction optical depth
c 
                 tauexc(i,k) = tauexc(i,k)*cld(i,k)*sqrt(cld(i,k))
                 fc(i,k)     = gc(i,k)*gc(i,k)
  230          continue
  229      continue
  228   continue

c
c set reflectivities for surface based on mid-point wavelength
c
        wavmid = 0.5*(wavmin(ns) + wavmax(ns))
c
c wavelength less  than 0.7 micro-meter
c
        if(wavmid .lt. 0.7 ) then
           do 412 n=1,nloop
              do 414 i=is(n),ie(n)
                 albdir(i)  = albs(i)
                 albdif(i)  = albsd(i)
  414         continue
  412      continue
c
c wavelength greater than 0.7 micro-meter
c
        else   
           do 422 n=1,nloop
              do 424 i=is(n),ie(n)
                 albdir(i)  = albl(i)
                 albdif(i)  = albld(i)
  424         continue
  422      continue
        endif
        trayoslp = raytau(ns)/sslp
c
c layer input properties now completely specified; compute the 
c delta-eddington solution reflectivities and transmissivities 
c for each layer, starting from the top and working downwards:
c
        call radded(coszrs   ,trayoslp,pflx    ,abh2o(ns),abo3(ns),
     $              abco2(ns),abo2(ns),uh2o    ,uo3      ,uco2    ,
     $              uo2      ,tauexc  ,wc      ,gc       ,fc      ,
     $              nloop    ,is      ,ie      ,rdir     ,rdif    ,
     $              tdir     ,tdif    ,explay  ,exptdn   ,rdndif  ,
     $              tottrn,
     &     bccmil,tccmil,numlay,ncclev,
     &     cbi,cti,ncldlay,
     &     idebug)
c
c compute reflectivity to direct and diffuse radiation for layers below 
c by adding succesive layers starting from the surface and working
c upwards:
c
        do 560 n=1,nloop
           do 570 i=is(n),ie(n)
              rupdir(i,plevp) = albdir(i)
              rupdif(i,plevp) = albdif(i)
  570      continue
  560   continue
        do 580 k=plev,0,-1
           do 582 n=1,nloop
              do 584 i=is(n),ie(n)
                 rdenom = 1./( 1. - rdif(i,k)*rupdif(i,k+1))
                 rupdir(i,k) = rdir(i,k) + tdif(i,k)*
     $                 (rupdir(i,k+1)*explay(i,k) +
     $                  rupdif(i,k+1)*(tdir(i,k)-explay(i,k)))*rdenom
                 rupdif(i,k) = rdif(i,k) +
     $                        rupdif(i,k+1)*tdif(i,k)**2*rdenom
  584         continue
  582      continue
  580   continue
c
c compute up and down fluxes for each interface, using the added
c atmospheric layer properties at each interface:
c
        do 600 k=0,plevp
           do 610 n=1,nloop
              do 620 i=is(n),ie(n)
                 rdenom = 1./(1. - rdndif(i,k)*rupdif(i,k))
                 fluxup(i,k) = (exptdn(i,k)*rupdir(i,k) +
     $                  (tottrn(i,k)-exptdn(i,k))*rupdif(i,k))*rdenom
                 fluxdn(i,k)=exptdn(i,k) + (tottrn(i,k) - exptdn(i,k) +
     $                  exptdn(i,k)*rupdir(i,k)*rdndif(i,k))*rdenom
  620         continue
  610      continue
  600   continue
c
c compute flux divergence in each layer using the interface 
c up and down fluxes:
c
        do 630 k=0,plev
           do 640 n=1,nloop
              do 650 i=is(n),ie(n)
                 flxdiv(i,k) = (fluxdn(i,k)-fluxdn(i,k+1)) +
     $                         (fluxup(i,k+1)-fluxup(i,k))
  650         continue
  640      continue
  630   continue
c
c monochromatic computation  completed; accumulate in totals;
c adjust fraction within spectral interval to allow for the
c possibility of sub-divisions within a particular interval:
c
        psf = 1.0
        if(ph2o(ns).ne.0.) psf = psf*ph2o(ns)
        if(pco2(ns).ne.0.) psf = psf*pco2(ns)
        if(po2 (ns).ne.0.) psf = psf*po2 (ns)
        do 656 n=1,nloop
           do 658 i=is(n),ie(n)
              solflx(i)  = solin(i)*frcsol(ns)*psf
              sabtp(i)   = sabtp(i) + solflx(i)*
     $                     (fluxdn(i,0) - fluxup(i,0))
              frsa(i)    = frsa(i)  + solflx(i)*
     $                     (fluxdn(i,plevp) - fluxup(i,plevp))
              sfltot(i)  = sfltot(i) + solflx(i)
              fswup(i,0) = fswup(i,0) + solflx(i)*fluxup(i,0)
              fswdn(i,0) = fswdn(i,0) + solflx(i)*fluxdn(i,0)
  658      continue
  656   continue
        do 680 k=0,plev
           do 682 n=1,nloop
              do 684 i=is(n),ie(n)
                 totfld(i,k)  = totfld(i,k)  + solflx(i)*flxdiv(i,k)
                 fswup(i,k+1) = fswup(i,k+1) + solflx(i)*fluxup(i,k+1)
                 fswdn(i,k+1) = fswdn(i,k+1) + solflx(i)*fluxdn(i,k+1)
  684         continue
  682      continue
  680   continue
c
cz    here's where the mean intensity field (at layer center) is computed
cz    and the units are converted MKS before returning to cloud model
cz    altitude loop runs down from top of cloud (so indices increase in F77)
        do 685 n=1,nloop
           do 686 i=is(n),ie(n)
              do 687 k=plev+1-bccmil-numlay,plev+1-bccmil-1
                 kcz=plev-k-bccmil+1
cz                 meanj(ns,kcz) = .5*
cz     &                cgsmks*(solflx(i)*
cz     &                (fluxup(i,k)+fluxdn(i,k)+
cz     &                fluxup(i,k+1)+fluxdn(i,k+1)))/
cz     &                (2.*pie)
cz    the -'ve sign is to switch swfdiv to normal coords, matching lwfdiv
                 swfdiv(ns,kcz) = -cgsmks*solflx(i)*flxdiv(i,k)
 687          continue
 686       continue
 685    continue
c
c
c following code is the diagnostic clear sky computation:
c
c compute delta-eddington solution reflectivities and transmissivities 
c for the entire column; note, for convenience, we use the same 
c reflectivity and transmissivity arrays as for the full calculation
c above, where 0 for layer quantities refers to the entire atmospheric
c column, and where 0 for interface quantities refers to top of atmos-
c phere, while 1 refers to the surface:
c
        call radclr(coszrs   ,trayoslp,pflx    ,abh2o(ns),abo3(ns),
     $              abco2(ns),abo2(ns),uth2o   ,uto3     ,utco2    ,
     $              uto2     ,nloop   ,is      ,ie       ,rdir     ,
     $              rdif     ,tdir    ,tdif    ,explay   ,exptdn   ,
     $              rdndif   ,tottrn  )
c
c compute reflectivity to direct and diffuse radiation for 
c entire column; 0,1 on layer quantities refers to two effective layers
c overlying surface; 0 on interface quantities refers to top of column; 
c 2 on interface quantities refers to the surface:
c
        do 705 n=1,nloop
           do 710 i=is(n),ie(n)
              rupdir(i,2) = albdir(i)
              rupdif(i,2) = albdif(i)
  710      continue
  705   continue
c
        do 715 k=1,0,-1
           do 720 n=1,nloop
              do 725 i=is(n),ie(n)
                 rdenom = 1./( 1. - rdif(i,k)*rupdif(i,k+1))
                 rupdir(i,k) = rdir(i,k) + tdif(i,k)*
     $                 (rupdir(i,k+1)*explay(i,k) +
     $                  rupdif(i,k+1)*(tdir(i,k)-explay(i,k)))*rdenom
                 rupdif(i,k) = rdif(i,k) +
     $                        rupdif(i,k+1)*tdif(i,k)**2*rdenom
  725         continue
  720      continue
  715   continue
c
c compute up and down fluxes for each interface, using the added
c atmospheric layer properties at each interface:
c
        do 730 k=0,2
           do 735 n=1,nloop
              do 740 i=is(n),ie(n)
                 rdenom = 1./(1. - rdndif(i,k)*rupdif(i,k))
                 fluxup(i,k) = (exptdn(i,k)*rupdir(i,k) +
     $                  (tottrn(i,k)-exptdn(i,k))*rupdif(i,k))*rdenom
                 fluxdn(i,k)=exptdn(i,k) + (tottrn(i,k) - exptdn(i,k) +
     $                  exptdn(i,k)*rupdir(i,k)*rdndif(i,k))*rdenom
  740         continue
  735      continue
  730   continue
c
        do 745 n=1,nloop
           do 750 i=is(n),ie(n)

              clrst(i)   = clrst(i) + solflx(i)*
     $                     (fluxdn(i,0) - fluxup(i,0))
              clrss(i)   = clrss(i) + solflx(i)*
     $                     (fluxdn(i,2) - fluxup(i,2))
cz Subtract the SW clear sky surface downwelling flux from the 
cz surface SW cloud forcing. 
cz
cz              swscf      = swscf-solflx(i)*fluxdn(i,2)
  750      continue
  745   continue
c
c end of clear sky computation
c
c
c end of spectral interval loop
c
  100 continue
cz    here's where the net up and down short-wave fluxes are stored
cz    in the cloud model arrays and the units are converted to MKS 
cz    before returning to cloud model. REMEMBER: this retrieves the
cz    interface fluxes, so cz index zero is true cloud base.
cz    altitude loop runs down from top of cloud (so indices increase in F77)
      do 691 i=1,plon
         do 692 k=plev+1-bccmil-numlay,plev+1-bccmil
            kcz=plev-k-bccmil+1
            fupsw(kcz) = cgsmks*fswup(i,k)
            fdsw(kcz) = cgsmks*fswdn(i,k)
cz    the following lines would return the mid-layer fluxes
cz         do 692 k=plev+1-bccmil-numlay,plev+1-bccmil-1
cz            fupsw(kcz) = .5*cgsmks*
cz     &           (fswup(i,k)+fswup(i,k+1))
cz            fdsw(kcz) = .5*cgsmks*
cz     &           (fswdn(i,k)+fswdn(i,k+1))
 692     continue
 691  continue
cz
cz Total up the surface cloud forcing and convert it to MKS.
cz swscf=F_SW_down_surf_cloudy-F_SW_down_surf_clear 
cz The clear sky terms have been subtracted within the spectral loop,
cz so now just add on the total all sky term.
cz      swscf=swscf+fswdn(1,plevp)
cz      swscf=swscf*cgsmks
cz
c
c compute solar heating rate (k/s)
c
      gocp = gravit/cpair
      do 800 k=1,plev
         do 802 n=1,nloop
           isn = is(n)
           ien = ie(n)
           do 804 i=isn,ien
              qrs(i,k) = -gocp*totfld(i,k)/(pint(i,k) - pint(i,k+1))
  804      continue
  802    continue
  800 continue
c
c done
c
      return
      end
      subroutine radded(coszrs  ,trayoslp,pflx    ,abh2o   ,abo3    ,
     $                  abco2   ,abo2    ,uh2o    ,uo3     ,uco2    ,
     $                  uo2     ,tauexc  ,wc      ,gc      ,fc      ,
     $                  nloop   ,is      ,ie      ,rdir    ,rdif    ,
     $                  tdir    ,tdif    ,explay  ,exptdn  ,rdndif  ,
     $                  tottrn,  
     &     bccmil,tccmil,numlay,ncclev,
     &     cbi,cti,ncldlay,
     &     idebug)
c-----------------------------------------------------------------------
c
c delta-Eddington solution
c
c computes layer reflectivities and transmissivities, from the
c top down to the surface using the delta-Eddington solutions for
c each layer
c
c if total transmission to the interface above a particular
c layer is less than trmin, then no further delta-Eddington solutions
c are evaluated for layers below
c
c for more details, see  Briegleb, B.P.   Delta-Eddington approximation
c for solar radiation in the NCAR community climate model, submitted
c to Journal of Geophysical Research
c
c-----------------------------------------------------------------------
      implicit none
c-----------------------------------------------------------------------
c
c     !   Basic grid point resolution parameters
c
      integer plev,plon,plat,plonp2,plevp,pcnst,plond,platd,plevd
      integer nxpt,jintmx,i1,j1,j1m
c
      parameter(plev=  133,        ! number of vertical levels
     $          plevp=plev+1,     ! number of vertical interfaces
     $          plon=1,           ! number of longitudes (T42)
     $          plat=1,           ! number of latitudes (T42)
     $          plonp2=plon+2,    ! length necessary for fft
     $          pcnst=1,          ! number of constituents
     $          nxpt=1,           !
     $          jintmx=2,         !
     $          plond=plon + 1 + 2*nxpt, ! slt extended domain longitude
     $          platd=plat + 2*nxpt + 2*jintmx, ! slt extended domain lat.
     $          plevd=plev*(3+pcnst), ! fold plev,pcnst indices into one
     $          i1=1+nxpt,        ! model starting longitude (3-d)
     $          j1=1+nxpt+jintmx, ! model starting latitude (3-d)
     $          j1m=j1-1)         ! model starting offset (3-d)
c
c-----------------------------------------------------------------------
c
c input arguments
c
      real coszrs(plond),         ! cosine zenith angle
     $     trayoslp,              ! tray/sslp
     $     pflx(plond,0:plevp),   ! interface pressure
     $     abh2o,                 ! absorption coefficiant for h2o
     $     abo3 ,                 ! absorption coefficiant for o3
     $     abco2,                 ! absorption coefficiant for co2
     $     abo2 ,                 ! absorption coefficiant for o2
     $     uh2o(plond,0:plev),    ! layer absorber amount of h2o
     $      uo3(plond,0:plev),    ! layer absorber amount of  o3
     $     uco2(plond,0:plev),    ! layer absorber amount of co2
     $      uo2(plond,0:plev)     ! layer absorber amount of  o2
      real tauexc(plond,0:plev),  ! cloud extinction optical depth
     $         wc(plond,0:plev),  ! cloud single scattering albedo
     $         gc(plond,0:plev),  ! cloud assymetry parameter
     $         fc(plond,0:plev)   ! cloud forward scattered fraction
      integer nloop,              ! number of loops (1 or 2)
     $        is(2),              ! starting index for 1 or 2 loops
     $        ie(2)               ! ending index for 1 or 2 loops
c
c input/output arguments
c
c following variables are defined for each layer; 0 refers to extrae
c layer above top of model:
c
      real rdir(plond,0:plev),   ! layer reflectivity to direct radiation
     $     rdif(plond,0:plev),   ! layer reflectivity to diffuse radiation
     $     tdir(plond,0:plev),   ! layer transmission to direct radiation
     $     tdif(plond,0:plev),   ! layer transmission to diffuse radiation
     $     explay(plond,0:plev)  ! solar beam exp transmission for layer
c
c (note that the following variables are defined on interfaces, with
c the index k referring to the top interface of the kth layer:
c exptdn,rdndif,tottrn; for example, tottrn(k=5) refers to the total
c transmission to the top interface of the 5th layer; plevp refers to the
c earth surface
c    
      real rdndif(plond,0:plevp), ! added dif ref for layers above
     $     exptdn(plond,0:plevp), ! solar beam exp down transmission from top
     $     tottrn(plond,0:plevp)  ! total transmission for layers above
c
      external  resetr,     ! resets array elements to zero
     $          whenfgt     ! collect indices for greater than condition
c
c Local workspace
c
      integer          i,   ! longitude index
     $                 k,   ! level index
     $                nn,   ! index of longitude loops (max=nloop)
     $                ii,   ! longitude index
     $              nval,   ! number of longitude values satisfying criteria
     $      index(plond)    ! array of longitude indices
c
      real taugab(plond),   ! layer total gas absorption optical depth
     $     tauray(plond),   ! layer rayleigh optical depth
     $     taucsc       ,   ! layer cloud scattering optical depth
     $     tautot       ,   ! total layer optical depth
     $       wtot       ,   ! total layer single scatter albedo
     $       gtot       ,   ! total layer asymmetry parameter
     $       ftot           ! total layer forward scatter fraction
c
c minimum total transmission below which no layer computation are done:
c
      real trmin,    ! transmission cutoff
     $      wray,    ! rayleigh single scatter albedo
     $      gray,    ! rayleigh asymetry parameter
     $      fray     ! rayleigh forward scattered fraction
c
      data trmin  /  1.e-3    /
      data  wray  /  0.999999 /
      data  gray  /  0.0      /
      data  fray  /  0.1      /
c
      real wtau,    !  rayleigh layer scattering optical depth
     $       wt,    !  layer total single scattering albedo
     $       ts,    !  layer scaled extinction optical depth
     $       ws,    !  layer scaled single scattering albedo
     $       gs     !  layer scaled asymmetry parameter
c
      real rdenom,    !  mulitiple scattering term
     $    rdirexp,    !  layer direct ref times exp transmission
     $    tdnmexp     !  total transmission minus exp transmission
c
c------------------------------------------------------------------------
c     zender cirrus cloud hack workspace
c
      integer 
     &     idebug,    !debugging level from cloud model
     &     kcz,       !index for index trickery
     &     bccmil,    !bottom ccm2 interface level
     &     tccmil,    !top ccm2 interface level
     &     numlay,    !num_layer in cloud model
     &     ncclev,    !num_ccm2_level (with overlaps already subtracted)
     &     cbi,       !cloud bottom index
     &     cti,       !cloud top index
     &     ncldlay    !number of layers in cloud itself
c
c     end of zender workspace alterations
c-----------------------------------------------------------------------
c
c statement functions for delta-Eddington solution for each layer
c
      real alpha,gamma,el,taus,omgs,asys,u,n,lm,ne
      real w,uu,g,e,f,t,et
c
c intermediate terms for delta-Eddington solution
c
      real alp,gam,ue,arg,extins,amg,apg
c
      alpha(w,uu,g,e) = .75*w*uu*((1. + g*(1-w))/(1. - e*e*uu*uu))
      gamma(w,uu,g,e) = .50*w*((3.*g*(1.-w)*uu*uu + 1.)/(1.-e*e*uu*uu))
      el(w,g)         = sqrt(3.*(1-w)*(1. - w*g))
      taus(w,f,t)     = (1. - w*f)*t
      omgs(w,f)       = (1. - f)*w/(1. - w*f)
      asys(g,f)       = (g - f)/(1. - f)
      u(w,g,e)        = 1.5*(1. - w*g)/e
      n(uu,et)        = ((uu+1.)*(uu+1.)/et ) - ((uu-1.)*(uu-1.)*et)
c
c-----------------------------------------------------------------------
c
c initialize all total transmission values to 0, so that nighttime values
c from previous computations are not used:
c
      call resetr(tottrn,plond*plevp,0.)
c
c compute total direct beam transmission, total transmission, and
c reflectivity for diffuse radiation (from below) for all layers
c above each interface by starting from the top and adding layers
c down:
c
c for the extrae layer above model top:
c
      do 10 nn=1,nloop
         do 20 i=is(nn),ie(nn)
c
            tauray(i) = trayoslp*(pflx(i,1)-pflx(i,0))
cz the top ccm2 layer absorption should probablly never be zeroed.
            taugab(i) = abh2o*uh2o(i,0) + abo3*uo3(i,0) +
     $                  abco2*uco2(i,0) + abo2*uo2(i,0)
c
            tautot  = tauexc(i,0) + tauray(i) + taugab(i)
            taucsc  = tauexc(i,0)*wc(i,0)
            wtau    = wray*tauray(i)
            wt      = wtau + taucsc
            wtot = wt/tautot
            gtot = (wtau*gray + gc(i,0)*taucsc)/wt
            ftot = (wtau*fray + fc(i,0)*taucsc)/wt
c
            ts   = taus(wtot,ftot,tautot)
            ws   = omgs(wtot,ftot)
            gs   = asys(gtot,ftot)
            lm   = el(ws,gs)
            alp  = alpha(ws,coszrs(i),gs,lm)
            gam  = gamma(ws,coszrs(i),gs,lm)
            ue   = u(ws,gs,lm)
c
c limit argument of exponential to 25, incase lm*ts very large:
c
            arg  = amin1(lm*ts,25.)
            extins = exp(-arg)
            ne = n(ue,extins)
c
            rdif(i,0) = (ue+1.)*(ue-1.)*(1./extins - extins)/ne
            tdif(i,0) = 4.*ue/ne
c
c limit argument of exponential to 25, incase coszrs is very small:
c
            arg       = amin1(ts/coszrs(i),25.)
            explay(i,0) = exp(-arg)
c
            apg = alp + gam
            amg = alp - gam
            rdir(i,0) = amg*(tdif(i,0)*explay(i,0) - 1.) + apg*rdif(i,0)
            tdir(i,0) = apg*tdif(i,0) +
     $                  (amg*rdif(i,0) - (apg-1.))*explay(i,0)
c
c under rare conditions, reflectivies and transmissivities can be 
c negative; zero out any negative values
c
            rdir(i,0) = amax1(rdir(i,0),0.0)
            tdir(i,0) = amax1(tdir(i,0),0.0)
            rdif(i,0) = amax1(rdif(i,0),0.0)
            tdif(i,0) = amax1(tdif(i,0),0.0)
c
c initialize top interface of extrae layer:
c
            exptdn(i,0) =   1.0
            rdndif(i,0) =   0.0
            tottrn(i,0) =   1.0
c
            rdndif(i,1) = rdif(i,0)
            tottrn(i,1) = tdir(i,0)
c
   20    continue
   10 continue
c
c now, continue down one layer at a time; if the total transmission
c to the interface just above a given layer is less than trmin, then
c no delta-Eddington computation for that layer is done:
c
      do 30 k=1,plev
c
c initialize current layer properties to zero; only if total transmission
c to the top interface of the current layer exceeds the minimum, will 
c these values be computed below:
c
         do 40 nn=1,nloop
            do 50 i=is(nn),ie(nn)
c
               rdir(i,k)   =  0.0
               rdif(i,k)   =  0.0
               tdir(i,k)   =  0.0
               tdif(i,k)   =  0.0
               explay(i,k) =  0.0
c
c calculates the solar beam transmission, total transmission, and
c reflectivity for diffuse radiation from below at the top of the
c current layer:
c
               exptdn(i,k) = exptdn(i,k-1)*explay(i,k-1)
               rdenom      = 1./(1. - rdif(i,k-1)*rdndif(i,k-1))
               rdirexp     = rdir(i,k-1)*exptdn(i,k-1)
               tdnmexp     = tottrn(i,k-1) - exptdn(i,k-1)
               tottrn(i,k) = exptdn(i,k-1)*tdir(i,k-1) + tdif(i,k-1)*
     $                      (tdnmexp + rdndif(i,k-1)*rdirexp)*rdenom
               rdndif(i,k) = rdif(i,k-1)  +
     $                (rdndif(i,k-1)*tdif(i,k-1))*(tdif(i,k-1)*rdenom)
c
   50       continue
   40    continue
c
c compute next layer delta-Eddington solution only if total transmission
c of radiation to the interface just above the layer exceeds trmin.
c
         call whenfgt(plon,tottrn(1,k),1,trmin,index,nval)
         if(nval.gt.0) then
CDIR$ IVDEP
            do 60 ii=1,nval
               i=index(ii)
c
               tauray(i) = trayoslp*(pflx(i,k+1)-pflx(i,k))
cz alert
cz               taugab(i) = abh2o*uh2o(i,k) + abo3*uo3(i,k) +
cz     $                     abco2*uco2(i,k) + abo2*uo2(i,k)
cz check whether this layer is with the cloud, if so then zero
cz the trace gas absorption
               kcz=plev-k-bccmil+1
               if(kcz.ge.cbi.and.kcz.le.cti)then 
                  taugab(i) = 0.
               else
                  taugab(i) = abh2o*uh2o(i,k) + abo3*uo3(i,k) +
     $                 abco2*uco2(i,k) + abo2*uo2(i,k)
               endif
c
               tautot = tauexc(i,k) + tauray(i) + taugab(i)
               taucsc    = tauexc(i,k)*wc(i,k)
               wtau      = wray*tauray(i)
               wt        = wtau + taucsc
               wtot   = wt/tautot
               gtot   = (wtau*gray + gc(i,k)*taucsc)/wt
               ftot   = (wtau*fray + fc(i,k)*taucsc)/wt
c
               ts   = taus(wtot,ftot,tautot)
               ws   = omgs(wtot,ftot)
               gs   = asys(gtot,ftot)
               lm   = el(ws,gs)
               alp  = alpha(ws,coszrs(i),gs,lm)
               gam  = gamma(ws,coszrs(i),gs,lm)
               ue   = u(ws,gs,lm)
c
c limit argument of exponential to 25, incase lm very large:
c
               arg  = amin1(lm*ts,25.)
               extins = exp(-arg)
               ne = n(ue,extins)
c
               rdif(i,k) = (ue+1.)*(ue-1.)*(1./extins - extins)/ne
               tdif(i,k)   =   4.*ue/ne
c
c limit argument of exponential to 25, incase coszrs is very small:
c
               arg       = amin1(ts/coszrs(i),25.)
               explay(i,k) = exp(-arg)
c
               apg = alp + gam
               amg = alp - gam
               rdir(i,k) = amg*(tdif(i,k)*explay(i,k) - 1.) +
     $                     apg*rdif(i,k)
               tdir(i,k) = apg*tdif(i,k) +
     $                     (amg*rdif(i,k) - (apg-1.))*explay(i,k)
c
c under rare conditions, reflectivies and transmissivities can be 
c negative; zero out any negative values
c
               rdir(i,k) = amax1(rdir(i,k),0.0)
               tdir(i,k) = amax1(tdir(i,k),0.0)
               rdif(i,k) = amax1(rdif(i,k),0.0)
               tdif(i,k) = amax1(tdif(i,k),0.0)
   60       continue
         endif
c
   30 continue
c
c compute total direct beam transmission, total transmission, and
c reflectivity for diffuse radiation (from below) for all layers
c above the surface:
c
      k = plevp
      do 70 nn=1,nloop
         do 80 i=is(nn),ie(nn)
            exptdn(i,k) = exptdn(i,k-1)*explay(i,k-1)
            rdenom = 1./(1. - rdif(i,k-1)*rdndif(i,k-1))
            rdirexp = rdir(i,k-1)*exptdn(i,k-1)
            tdnmexp = tottrn(i,k-1) - exptdn(i,k-1)
            tottrn(i,k) = exptdn(i,k-1)*tdir(i,k-1) + tdif(i,k-1)*
     $                   (tdnmexp + rdndif(i,k-1)*rdirexp)*rdenom
            rdndif(i,k) = rdif(i,k-1)  +
     $               (rdndif(i,k-1)*tdif(i,k-1))*(tdif(i,k-1)*rdenom)
   80    continue
   70 continue
c
c done
c 
      return
      end
      subroutine radclr(coszrs  ,trayoslp,pflx    ,abh2o   ,abo3    ,
     $                  abco2   ,abo2    ,uth2o   ,uto3    ,utco2   ,
     $                  uto2    ,nloop   ,is      ,ie      ,rdir    ,
     $                  rdif    ,tdir    ,tdif    ,explay  ,exptdn  ,
     $                  rdndif  ,tottrn  )
c-----------------------------------------------------------------------
c
c delta-Eddington solution for special clear sky computation
c
c computes total reflectivities and transmissivities for two atmospheric
c layers: an overlying purely ozone absorbing layer, and the rest of the
c column below.
c
c for more details, see  Briegleb, B.P.   Delta-Eddington approximation
c for solar radiation in the NCAR community climate model, submitted
c to Journal of Geophysical Research
c
c-----------------------------------------------------------------------
      implicit none 
c-----------------------------------------------------------------------
c
c     !   Basic grid point resolution parameters
c
      integer plev,plon,plat,plonp2,plevp,pcnst,plond,platd,plevd
      integer nxpt,jintmx,i1,j1,j1m
c
      parameter(plev=  133,        ! number of vertical levels
     $          plevp=plev+1,     ! number of vertical interfaces
     $          plon=1,           ! number of longitudes (T42)
     $          plat=1,           ! number of latitudes (T42)
     $          plonp2=plon+2,    ! length necessary for fft
     $          pcnst=1,          ! number of constituents
     $          nxpt=1,           !
     $          jintmx=2,         !
     $          plond=plon + 1 + 2*nxpt, ! slt extended domain longitude
     $          platd=plat + 2*nxpt + 2*jintmx, ! slt extended domain lat.
     $          plevd=plev*(3+pcnst), ! fold plev,pcnst indices into one
     $          i1=1+nxpt,        ! model starting longitude (3-d)
     $          j1=1+nxpt+jintmx, ! model starting latitude (3-d)
     $          j1m=j1-1)         ! model starting offset (3-d)
c
c-----------------------------------------------------------------------
c
c input arguments
c
      real coszrs(plond),         ! cosine zenith angle
     $     trayoslp,              ! tray/sslp
     $     pflx(plond,0:plevp),   ! interface pressure
     $     abh2o,                 ! absorption coefficiant for h2o
     $     abo3 ,                 ! absorption coefficiant for o3
     $     abco2,                 ! absorption coefficiant for co2
     $     abo2 ,                 ! absorption coefficiant for o2
     $     uth2o(plond),          ! total column absorber amount of h2o
     $     uto3(plond),           ! total column absorber amount of  o3
     $     utco2(plond),          ! total column absorber amount of co2
     $     uto2(plond)            ! total column absorber amount of  o2
      integer nloop,              ! number of loops (1 or 2)
     $        is(2),              ! starting index for 1 or 2 loops
     $        ie(2)               ! ending index for 1 or 2 loops
c
c input/output arguments
c
c following variables are defined for each layer; note, we use layer 0 to
c refer to the entire atmospheric column:
c
      real rdir(plond,0:plev),   ! layer reflectivity to direct radiation
     $     rdif(plond,0:plev),   ! layer refflectivity to diffuse radiation
     $     tdir(plond,0:plev),   ! layer transmission to direct radiation
     $     tdif(plond,0:plev),   ! layer transmission to diffuse radiation
     $     explay(plond,0:plev)  ! solar beam exp transmission for layer
c
c (note that the following variables are defined on interfaces, with
c the index k referring to the top interface of the kth layer:
c exptdn,rdndif,tottrn; for example, tottrn(k=5) refers to the total
c transmission to the top interface of the 5th layer.
c
      real exptdn(plond,0:plevp), ! solar beam exp down transmission from top
     $     rdndif(plond,0:plevp), ! added dif ref for layers above
     $     tottrn(plond,0:plevp)  ! total transmission for layers above
c
      external  resetr,     ! resets array elements to zero
     $          whenfgt     ! collect indices for greater than condition
c
c Local workspace
c
      integer          i,   ! longitude index
     $                 k,   ! level index
     $                nn,   ! index of longitude loops (max=nloop)
     $                ii,   ! longitude index
     $              nval,   ! number of longitude values satisfying criteria
     $      index(plond)    ! array of longitude indices
c
      real taugab(plond),   ! total column gas absorption optical depth
     $     tauray(plond),   ! column rayleigh optical depth
     $     tautot       ,   ! total column optical depth
     $       wtot       ,   ! total column single scatter albedo
     $       gtot       ,   ! total column asymmetry parameter
     $       ftot           ! total column forward scatter fraction
c
c minimum total transmission below which no layer computation are done:
c
      real  trmin,   ! minimum total transmission allowed
     $      wray,    ! rayleigh single scatter albedo
     $      gray,    ! rayleigh asymetry parameter
     $      fray     ! rayleigh forward scattered fraction
c
      data trmin  /  1.e-3    /
      data  wray  /  0.999999 /
      data  gray  /  0.0      /
      data  fray  /  0.1      /
c
      real   ts,    !  column scaled extinction optical depth
     $       ws,    !  column scaled single scattering albedo
     $       gs     !  column scaled asymmetry parameter
c
      real rdenom,    !  mulitiple scattering term
     $    rdirexp,    !  layer direct ref times exp transmission
     $    tdnmexp     !  total transmission minus exp transmission
c
c statement functions for delta-Eddington solution for entire column:
c
      real alpha,gamma,el,taus,omgs,asys,u,n,lm,ne
      real w,uu,g,e,f,t,et
c
c intermediate terms for delta-Eddington solution
c
      real alp,gam,ue,arg,extins,amg,apg
c
      alpha(w,uu,g,e) = .75*w*uu*((1. + g*(1-w))/(1. - e*e*uu*uu))
      gamma(w,uu,g,e) = .50*w*((3.*g*(1.-w)*uu*uu + 1.)/(1.-e*e*uu*uu))
      el(w,g)         = sqrt(3.*(1-w)*(1. - w*g))
      taus(w,f,t)     = (1. - w*f)*t
      omgs(w,f)       = (1. - f)*w/(1. - w*f)
      asys(g,f)       = (g - f)/(1. - f)
      u(w,g,e)        = 1.5*(1. - w*g)/e
      n(uu,et)        = ((uu+1.)*(uu+1.)/et ) - ((uu-1.)*(uu-1.)*et)
c
c-----------------------------------------------------------------------
c
c initialize all total transmission values to 0, so that nighttime values
c from previous computations are not used:
c
      call resetr(tottrn,plond*2,0.)
c
c compute total direct beam transmission, total transmission, and
c reflectivity for diffuse radiation (from below) for all layers
c above each interface by starting from the top and adding layers
c down:
c
c the top layer is assumed to be a purely absorbing ozone layer, and
c that the mean diffusivity for diffuse transmission is 1.66:
c
      do 10 nn=1,nloop
         do 20 i=is(nn),ie(nn)
c
            taugab(i) = abo3*uto3(i) 
c
c limit argument of exponential to 25, incase coszrs is very small:
c
            arg         = amin1(taugab(i)/coszrs(i),25.)
            explay(i,0) = exp(-arg)
            tdir(i,0)   = explay(i,0)
c
c same limit for diffuse transmission:
c
            arg         = amin1(1.66*taugab(i),25.)
            tdif(i,0)   = exp(-arg)
c
            rdir(i,0)   = 0.0
            rdif(i,0)   = 0.0
c
c initialize top interface of extrae layer:
c
            exptdn(i,0) =   1.0
            rdndif(i,0) =   0.0
            tottrn(i,0) =   1.0
c
            rdndif(i,1) = rdif(i,0)
            tottrn(i,1) = tdir(i,0)
c
   20    continue
   10 continue
c
c now, complete the rest of the column; if the total transmission
c through the top ozone layer is less than trmin, then
c no delta-Eddington computation for the underlying column is done:
c
      do 30 k=1,1
c
c initialize current layer properties to zero; only if total transmission
c to the top interface of the current layer exceeds the minimum, will 
c these values be computed below:
c
         do 40 nn=1,nloop
            do 50 i=is(nn),ie(nn)
c
               rdir(i,k)   =  0.0
               rdif(i,k)   =  0.0
               tdir(i,k)   =  0.0
               tdif(i,k)   =  0.0
               explay(i,k) =  0.0
c
c calculates the solar beam transmission, total transmission, and
c reflectivity for diffuse radiation from below at the top of the
c current layer:
c
               exptdn(i,k) = exptdn(i,k-1)*explay(i,k-1)
               rdenom      = 1./(1. - rdif(i,k-1)*rdndif(i,k-1))
               rdirexp     = rdir(i,k-1)*exptdn(i,k-1)
               tdnmexp     = tottrn(i,k-1) - exptdn(i,k-1)
               tottrn(i,k) = exptdn(i,k-1)*tdir(i,k-1) + tdif(i,k-1)*
     $                      (tdnmexp + rdndif(i,k-1)*rdirexp)*rdenom
               rdndif(i,k) = rdif(i,k-1)  +
     $                (rdndif(i,k-1)*tdif(i,k-1))*(tdif(i,k-1)*rdenom)
c
   50       continue
   40    continue
c
c compute next layer delta-Eddington solution only if total transmission
c of radiation to the interface just above the layer exceeds trmin.
c
         call whenfgt(plon,tottrn(1,k),1,trmin,index,nval)
         if(nval.gt.0) then
CDIR$ IVDEP
            do 60 ii=1,nval
               i=index(ii)
c
c remember, no ozone absorption in this layer:
c
               tauray(i) = trayoslp*pflx(i,plevp)
               taugab(i) = abh2o*uth2o(i) +
     $                     abco2*utco2(i) + abo2*uto2(i)
c
               tautot    = tauray(i) + taugab(i)
               wtot      = (wray*tauray(i))/tautot
               gtot      = gray
               ftot      = fray
c
               ts        = taus(wtot,ftot,tautot)
               ws        = omgs(wtot,ftot)
               gs        = asys(gtot,ftot)
               lm        = el(ws,gs)
               alp       = alpha(ws,coszrs(i),gs,lm)
               gam       = gamma(ws,coszrs(i),gs,lm)
               ue        = u(ws,gs,lm)
c
c limit argument of exponential to 25, incase lm very large:
c
               arg       = amin1(lm*ts,25.)
               extins    = exp(-arg)
               ne        = n(ue,extins)
c
               rdif(i,k) = (ue+1.)*(ue-1.)*(1./extins - extins)/ne
               tdif(i,k) =   4.*ue/ne
c
c limit argument of exponential to 25, incase coszrs is very small:
c
               arg       = amin1(ts/coszrs(i),25.)
               explay(i,k) = exp(-arg)
c
               apg       = alp + gam
               amg       = alp - gam
               rdir(i,k) = amg*(tdif(i,k)*explay(i,k) - 1.) +
     $                     apg*rdif(i,k)
               tdir(i,k) = apg*tdif(i,k) +
     $                     (amg*rdif(i,k) - (apg-1.))*explay(i,k)
c
c under rare conditions, reflectivies and transmissivities can be 
c negative; zero out any negative values
c
               rdir(i,k) = amax1(rdir(i,k),0.0)
               tdir(i,k) = amax1(tdir(i,k),0.0)
               rdif(i,k) = amax1(rdif(i,k),0.0)
               tdif(i,k) = amax1(tdif(i,k),0.0)
   60       continue
         endif
c
   30 continue
c
c compute total direct beam transmission, total transmission, and
c reflectivity for diffuse radiation (from below) for both layers
c above the surface:
c
      k = 2
      do 70 nn=1,nloop
         do 80 i=is(nn),ie(nn)
            exptdn(i,k) = exptdn(i,k-1)*explay(i,k-1)
            rdenom = 1./(1. - rdif(i,k-1)*rdndif(i,k-1))
            rdirexp = rdir(i,k-1)*exptdn(i,k-1)
            tdnmexp = tottrn(i,k-1) - exptdn(i,k-1)
            tottrn(i,k) = exptdn(i,k-1)*tdir(i,k-1) + tdif(i,k-1)*
     $                   (tdnmexp + rdndif(i,k-1)*rdirexp)*rdenom
            rdndif(i,k) = rdif(i,k-1)  +
     $               (rdndif(i,k-1)*tdif(i,k-1))*(tdif(i,k-1)*rdenom)
   80    continue
   70 continue
c
c done
c
      return
      end
      subroutine wheneq(plon,ipos,i1,i2,indx,npts)
c....... dummy routine for cray routine
      integer plon,ipos(plon),indx(plon)
       npts = 0
       do 10 i=1,plon
        if( ipos(i) .eq. i2 ) then
          npts       = npts + 1
          indx(npts) = i
        endif
   10  continue
      return
      end
      subroutine whenne(plon,ipos,i1,i2,indx,npts)
c....... dummy routine for cray routine
      integer plon,ipos(plon),indx(plon)
       npts = 0
       do 10 i=1,plon
        if( ipos(i) .ne. i2 ) then
          npts       = npts + 1
          indx(npts) = i
        endif
   10  continue
      return
      end
      subroutine whenfgt(nmbr,tin,len,tmin,index,nval)
c
c....... dummy routine for cray routine
c
      dimension tin(nmbr),index(nmbr)
      nval     = 0
      index(1) = 0
      do 100 n=1,nmbr
        if( tin(n) .gt. tmin ) then
          nval        = nval + 1
          index(nval) = n
        endif
  100 continue
      return
      end
      subroutine whenflt(nmbr,tin,len,tmin,index,nval)
c
c....... dummy routine for cray routine
c
      dimension tin(nmbr),index(nmbr)
      nval     = 0
      index(1) = 0
      do 100 n=1,nmbr
        if( tin(n) .lt. tmin ) then
          nval        = nval + 1
          index(nval) = n
        endif
  100 continue
      return
      end
      integer function isrchfgt(number,array,indx,value)
c
c search for first array element in array (up to number
c of elements), starting with index indx, that is greater 
c than value; if not found, set to number+1; otherwise, index
c
      real array(number)
c
      isrchfgt = 0
      do 100 n=1,number
        if( array(indx+n-1) .gt. value ) then
          isrchfgt = indx+n-1
          goto 101
        endif
  100 continue
  101 if( isrchfgt .eq. 0 ) then
       isrchfgt = number + 1
      endif
      return
      end
      integer function isrchfle(number,array,indx,value)
c
c search for first array element in array (up to number
c of elements), starting with index indx, that is greater 
c than value; if not found, set to zero; otherwise, index
c
      real array(number)
c
      isrchfle = 0
      do 100 n=1,number
        if( array(indx+n-1) .gt. value ) then
          isrchfle = indx+n-1
          goto 101
        endif
  100 continue
  101 if( isrchfle .eq. 0 ) then
       isrchfle = number + 1
      endif
      return
      end
      subroutine resetr(a,length,val)
      dimension a(1)
       do 10 n=1,length
         a(n) = val
   10  continue
      return
      end
      subroutine writeric(nabem,absems,lngbuf,nrow)
c
c........ dummy routine for write of ab/em data
c
      return 
      end
      subroutine readric(nabem,absems,lngbuf,nrow)
c
c........ dummy routine for read of ab/em data
c
      return 
      end
