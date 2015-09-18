      program fct
c
      parameter(n=200)     
c
      dimension rhoo(n),rhon(n)
      dimension u(n)
      dimension radn(n)
c
      real dt
      real lbc,rbc
      real radl,radr
c
      integer ul,ur
c
c

      cloud_base=10000.
      cloud_thick=1000.
      bin_size=20.e-6	
      crystal_length=50.e-6
      crystal_mass=mass_of_length(crystal_length)
      conc_LBC_add=0.
      conc_RBC_add=0.
      conc_LBC_mult=0.
      conc_RBC_mult=0.
      alpha = .2
      bpd = 10.
      num_layer = 50
      nsteps = 60
      lbc=0.
      rbc=0.
c
c     set up the grid variables
      call ngride(raqdn, n, radl, radr, alpha)
c     compute the velocity variables
      call veloce(u, n, ul, ur, dt)
c     find out what the conserved sum is
      call consre(rho,n,csum)
c     call the advection scheme
      call etbfct(rhoo,rhon,n,lbc,rbc)
c     copy the new values over the old values
      call ogride(n)
c
      end
c

