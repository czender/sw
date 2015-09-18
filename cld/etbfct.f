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
c Revision 1.1  1993/04/28  00:02:50  zender
c Initial revision
c
c
      subroutine etbfct(rhoo,rhon,n,lbc,rbc)
c
c
c     **************************************************************************
c
c     etbfct(rhoo,rhon,n,lbc,rbc)
c     originator-j.p. boris        oct.1975
c
c     this version has been modified to compile on the cyber 180-855
c     and uses the cyber sign and abs intrinsic functions in loop 21
c     rather than mask1,-2 and -3 as in the original code.
c     grt (august 28, 1985)
c
c     loop optimization for cyber started october 1, 1985. (grt)
c
c
c     description:  this routine solves generalized continuity equations
c     of the form   drho/dt=-div(rho*v)+sources   in either 
c     cartesian, cylindrical or spherical geometry.  the finite-difference
c     grid can be eulerian, sliding rezone, or lagrangian and can be
c     arbitrarily spaced.  the algorithm used is a low-phase-error fct
c     algorithm, vectorized and optimized for speed.  a complete
c     description appears in nrl memo report #3237, "flux-corrected
c     transport modules for solving generalized continuity equations."
c
c     arguments:  in this routine the right boundary at radr is half a
c     cell beyond the last grid point n at radn(n), and the left boundary
c     at radl is half a cell before the first grid point at radn(1).
c     rhoo     real array(n)   grid point densities at start of step   i
c     rhon     real array(n)   grid point densities at end of step     o
c     n        integer         number of interior grid points          i
c     lbc      real            left boundary condition factor          i
c     rbc      real            right boundary condition factor         i
c
c     language and limitations:  the subroutine etbfct is a multiple-
c     entry fortran routine is single precision (60 bits cyber).  the 
c     parameter statement is used to set symbolically the internal array
c     dimensions.  underflows are possible when the function being
c     transported has many zeroes.  the calculations generally misconserve
c     by one or two bits per cycle.  the relative phase and amplitude 
c     errors (for smooth functions) are typically several percent for 
c     characteristic lengths of 1-2 cells (wavelengths of order 10 cells).
c
c     entry points:  ogride, ngride, veloce, source, consre.  see # 11
c     of the detailed documention (or the listing below) for the
c     explanation and use of the arguments to these other entries.
c
c     no auxiliary or library routines are called by etbfct.
c
c     modifications: 
c              geometry additions in ngride to include the
c              cylindrical geometry for the cloud model.
c                      alpha=  4       inner cylinder  (arad)
c                              5       outer region    (brad)
c              the areas and volumes are computed with the
c              cloud areas in common /varcon/, i.e., areaa
c              and areab, calculated in subroutine readin.
c              (november 14, 1985...grt)
c
c     ********************************************************************
c
c
c
c
      real     rhoo(n), rhon(n)
      real     u(n), radn(n), c(n), d(n), rho(n)
      integer  alpha
c
      parameter (npt=202)
      logical  lsourc
      real     rbc, lbc
      real     source(npt), scrh(npt), rhot(npt), diff(npt) 
      real     adugth(npt), flxh(npt), nulh(npt), mulh(npt) 
      real     lnrhot(npt), fsgn(npt), fabs(npt), epsh(npt) 
      real     lorhot(npt), terp(npt), term(npt), adudth(npt)
      real     lo(npt), ln(npt), lh(npt), rlo(npt), rln(npt)
      real     rnh(npt), roh(npt), rlh(npt), ah(npt)
c
c
      equivalence (epsh(1), scrh(1))
      equivalence  (flxh(1), scrh(1)), (fsgn(1), rhot(1))
      equivalence  (fabs(1), scrh(1)), (terp(1), source(1)) 
      equivalence  (term(1), source(1)) 
c
c
c     mask1,2,3 are hexadecimal for 32 bit machine
c     data mask1,mask2,mask3 /z"80000000",0.998,z"7ffffffff"/
      data pi/3.141592654/,ftpi/4.188790207/
      data roh/npt*1.0/,source/npt*0.0/, lsourc/.false./
	data tlim/1.0e+10/
c
c
c  calculate the diffusive and convective fluxes
      np=n+1
      do 11 i=2,n
      flxh(i)=0.5*adudth(i)*(rhoo(i)+rhoo(i-1))
 11   diff(i)=nulh(i)*(rhoo(i)-rhoo(i-1))
      rhol=rhoo(1)*lbc
      rhor=rhoo(n)*rbc
      diff(1)=nulh(1)*(rhoo(1)-rhol)
      diff(np)=nulh(np)*(rhor-rhoo(n))
      flxh(1)=0.5*adudth(1)*(rhoo(1)+rhol)
      flxh(np)=0.5*adudth(np)*(rhor+rhoo(n))
c
c  calculate lambdao*rhot, the transported mass elements.
      do 12 i=1,n
        	lorhot(i)=lo(i)*rhoo(i)-flxh(i+1)+flxh(i)
		if(abs(lorhot(i)).lt.1.0e-15) lorhot(i)=0.00 
12	continue
c
c  add in the source terms as appropriate.
      if(.not. lsourc) go to 14
      do 13 i=1,n
 13   lorhot(i)=lorhot(i)+source(i)
c
c  calculate the phoenical antidiffusive fluxes here.
 14   rhot(1)=lorhot(1)*rlo(1)
      do 17 i=2,n
      rhot(i)=lorhot(i)*rlo(i)
 17   flxh(i)=mulh(i)*(rhot(i)-rhot(i-1))
      flxh(1)=mulh(1)*(rhot(1)-lbc*rhot(1))
      flxh(np)=mulh(np)*(rbc*rhot(n)-rhot(n))
c
c  diffuse the solution rhot using old fluxes.
      lnrhot(1)=lorhot(1)+diff(2)-diff(1)
      rhot(1)=lnrhot(1)*rln(1)
      do 20 i=2,n
      lnrhot(i)=lorhot(i)+diff(i+1)-diff(i)
c
c  calculate the transported/diffused density and grid differences.
      rhot(i)=lnrhot(i)*rln(i)
 20   diff(i)=rhot(i)-rhot(i-1)
      diff(1)=rhot(1)-lbc*rhot(1)
      diff(np)=rbc*rhot(n)-rhot(n)
c
c  calculate the sign and magnitude of the antidiffusive flux.
      do 21 i=1,np
      fabs(i)=abs(flxh(i))
 21   fsgn(i)=sign(1.0,diff(i))
c
c  calculate the flux-limiting changes on the right and the left
c
c  terp(n) and term(2) are set to 1.0e+75
c  to force the antidiffusion to remove most
c  of the diffusion at the boundaries.
c  (grt 12-85)
      do 23 i=1,n-1 
 23   terp(i)=fsgn(i)*ln(i)*diff(i+1)
c     terp(n)=1.0e+75
c     terp(np)=1.0e+75
      terp(n)=tlim
      terp(np)=tlim
      do 24 i=1,np
 24   fabs(i)=amin1(terp(i),fabs(i))
      do 25 i=3,np
 25   term(i)=fsgn(i)*ln(i-1)*diff(i-1) 
      term(1)=tlim
      term(2)=tlim
c
c  correct the fluxes completely now.
      do 35 i=1,np
 35   diff(i)=amin1(fabs(i),term(i))
      do 36 i=1,np
 36   flxh(i)=amax1(0.0, diff(i))
      do 37 i=1,np
 37   flxh(i)=fsgn(i)*flxh(i) 
c
c
c
c  calculate the new flux-corrected densities.
      do 41 i=1,n
      source(i)=0.0 
 41   rhon(i)=rln(i)*(lnrhot(i)-flxh(i+1)+flxh(i))
      lsourc=.false.
      return
c
c
c     --------------------------------------------------------------------------
c
      entry veloce(u, n, ul, ur, dt)
c
c
c     **************************************************************************
c
c     veloce (u,n,ul,ur,dt)
c     description:  this entry calculates all velocity-dependant coeffs.
c
c     arguments: 
c     u        real array(n)   flow velocity at the grid point         i
c     n        integer         number of interior grid points          i
c     ul       real            velocity of flow at left boundary       i
c     ur       real            velocity of flow at right boundary      i
c     dt       real            stepsize for the time integration       i
c
c     **************************************************************************
c
c
c  calculate the interface area x velocity differential x dt.
      np=n+1
      dth=0.5*dt
      do 101 i=2,n
 101  adudth(i)=ah(i)*dth*(u(i)+u(i-1))-adugth(i) 
      adudth(1)=ah(1)*dt*ul-adugth(1)
      adudth(np)=ah(np)*dt*ur-adugth(np)
c
c  calculate the half-cell epsilon (v*dt/dx)
      do 104 i=1,np 
      epsh(i)=adudth(i)*rlh(i)
c
c  next calculate the diffusion and antidiffusion coefficients.
c  variation with epsilon means fourth-order accurate phases.
      nulh(i)=0.1666666667+0.3333333333*epsh(i)*epsh(i)
      mulh(i)=0.25-0.5*nulh(i)
      nulh(i)=lh(i)*nulh(i)
 104  mulh(i)=lh(i)*mulh(i)
      return
c
c
c     --------------------------------------------------------------------------
c
      entry ngride(radn, n, radl, radr, alpha)
c
c
c     **************************************************************************
c
c     ngride(raqdn, n, radl, radr, alpha)
c     description:  this entry sets new geometry variables and coeffs.
c
c     arguments: 
c     radn     real array(n)   new grid point positions                i
c     n        integer         number of interior grid points          i
c     radl     real            position of the left boundary           i
c     radr     real            position of the right boundary          i
c     alpha    integer         =1  for cartesian geometry              i
c                              =2  for cylindrical geometry            i
c                              =3  for spherical geometry              i
c                              =4  for inner cloud region              i
c                              =5  for outer cloud region              i
c     **************************************************************************
c
c
c
c  calculate the new half-cell positions and grid changes.
      np=n+1
      do 202 i=2,n
 202  rnh(i)=0.5*(radn(i)+radn(i-1))
      rnh(1)=radl
      rnh(np)=radr
c
c  calculate the three coordinate systems
      go to (203, 206, 209), alpha
c
c  cartesian coordinates.
 203  do 204 i=1,np 
 204  ah(i)=1.0
      do 205 i=1,n
 205  ln(i)=rnh(i+1)-rnh(i)
      go to 250
c
c  cylindrical coordinates.
 206  do 207 i=1,np 
      diff(i)=rnh(i)*rnh(i)
 207  ah(i)=pi*(roh(i)+rnh(i))
      do 208 i=1,n
 208  ln(i)=pi*(diff(i+1)-diff(i))
      go to 250
c
c  spherical coordinates
 209  do 210 i=1,np 
      diff(i)=rnh(i)*rnh(i)*rnh(i)
 210  scrh(i)=(roh(i)+rnh(i))*roh(i)
      do 211 i=1,np 
 211  ah(i)=ftpi*(scrh(i)+rnh(i)*rnh(i))
      do 212 i=1,n
 212  ln(i)=ftpi*(diff(i+1)-diff(i))
      go to 250
c  now the geometric variables which are system independant.
 250  do 255 i=2,n
 255  lh(i)=0.5*(ln(i)+ln(i-1))
      lh(1)=ln(1)
      lh(np)=ln(n)
      do 260 i=1,n
 260  rln(i)=1.0/ln(i)
      do 265 i=1,np 
 265  adugth(i)=ah(i)*(rnh(i)-roh(i))
      do 270 i=2,n
 270  rlh(i)=0.5*(rln(i)+rln(i-1))
      rlh(1)=rln(1) 
      rlh(np)=rln(n)
      return
c
c
c     --------------------------------------------------------------------------
c
      entry sources(n, dt, modes, c, d, dl, dr)
c
c     **************************************************************************
c
c     sources(n, dt, modes, c, d, dl, dr)
c     description:  this entry accumulates different source terms.
c
c     arguments: 
c     n        integer         number of interior grid points          i
c     dt       real            stepsize for the time integration       i
c     modes    integer         =1  computes +div(d)                    i
c                              =2  computes +c*grad(d)                 i
c                              =3  adds +d to the sources              i
c     c        real array(n)   array of source variables at grid pts   i
c     d        real array(n)   array of source variables at grid pts   i
c     dl       real            left boundary value of d                i
c     dr       real            right boundary value of d               i
c
c     **************************************************************************
c
c
c
      np=n+1
      dth=0.5*dt
      dtq=0.25*dt
      go to (310, 320, 330), modes
c
c  +div(d) is computed conservatively and added to the sources.
 310  do 311 i=2,n
 311  scrh(i)=dth*ah(i)*(d(i)+d(i-1))
      scrh(1)=dt*ah(1)*dl
      scrh(np)=dt*ah(np)*dr
      do 312 i=n,1,-1
 312  source(i)=source(i)+scrh(i+1)-scrh(i)
      lsourc=.true. 
      return
c
c  +c*grad(d) is computed efficiently and added to the sources.
 320  do 321 i=2,n
 321  scrh(i)=dtq*(d(i)+d(i-1))
      scrh(1)=dth*dl
      scrh(np)=dth*dr
      do 323 i=n,1,-1
      diff(i)=scrh(i+1)-scrh(i)
 323  source(i)=source(i)+c(i)*(ah(i+1)+ah(i))*diff(i)
      lsourc=.true. 
      return
c
c  +d is added to the sources in an explicit formulation.
 330  do 331 i=1,n
 331  source(i)=source(i)+dt*lo(i)*d(i) 
      lsourc=.true. 
      return
c
c
c     --------------------------------------------------------------------------
c
      entry ogride(n)
c
c     **************************************************************************
c
c     ogride(n)
c     description:  this entry copies old grid and geometry variables
c
c     arguments: 
c     n        integer         number of interior grid points          i
c
c     **************************************************************************
c
c  copy the previously new grid values to be used as the old grid
      np=n+1
      do 401 i=1,n
      lo(i)=ln(i)
 401  rlo(i)=rln(i) 
      do 402 i=1,np 
 402  roh(i)=rnh(i) 
      return
c
c
c     --------------------------------------------------------------------------
c
      entry consre(rho, n, csum)
c
c     **************************************************************************
c
c     consre(rho, n, csum)
c     description:  this entry computes the ostensibly conserved sum.
c
c     arguments: 
c     rho      real array(n)   grid point values for conservation sum  i
c     n        integer         number of interior grid points          i
c     csum     real            value of the conservation sum of rho    o
c
c     **************************************************************************
c
c
c
c  compute the ostensibly conserved total mass (beware your b.c.)
      csum=0.0
      do 501 i=1,n
 501  csum=csum+ln(i)*rho(i)
      return
c
c
      end 
