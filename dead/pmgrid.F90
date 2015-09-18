! $Id$ -*-f90-*-

! Purpose: Supply spatial resolution parameters as Fortran variables
! Pre-processor defined tokens are stored in Fortran variables here

! Usage: 
! use pmgrid ! [mdl] Spatial resolution parameters

#include <params.h> /* Preprocessor tokens */

module pmgrid ! [mdl] Spatial resolution parameters
  implicit none
  save ! [stt] Changes to common variables are sticky
  
  ! Parameters
#ifdef BXM
  integer,parameter::pcnst= PCNST ! [nbr] number of advected constituents
  integer,parameter::plon= PLON ! [nbr] number of longitudes
  integer,parameter::plat= PLAT ! [nbr] number of latitudes
  integer,parameter::plev= PLEV ! [nbr] number of vertical levels
  integer,parameter::plevp= plev + 1 ! [nbr] plev + 1
  integer,parameter::plond= plon ! [nbr] slt extended domain longitude
#endif /* endif BXM */
  
#if (defined CCM) && (!defined BXM)
  integer,parameter::plon   = PLON
  integer,parameter::plev   = PLEV
  integer,parameter::plat   = PLAT
  integer,parameter::pcnst  = PCNST
  integer,parameter::pnats  = PNATS
  integer,parameter::plevmx = 4
  integer,parameter::plevp  = plev + 1
  integer,parameter::nxpt   = 1
  integer,parameter::jintmx = 1
  integer,parameter::plond  = plon + 1 + 2*nxpt
  integer,parameter::platd  = plat + 2*nxpt + 2*jintmx
  integer,parameter::p3d    = 3 + pcnst + pnats
  integer,parameter::plevd  = plev*p3d
  integer,parameter::i1     = 1 + nxpt
  integer,parameter::j1     = 1 + nxpt + jintmx
  integer,parameter::numbnd = nxpt + jintmx

#if ( defined SPMD )
  ! Re-do SPMD
#endif

#endif /* endif CCM */
  
#if (!defined CCM) && (!defined BXM)
  
#ifndef PARAMS_H
#include <params.h>
#endif
  
  ! Set kind of real variables to have at least 12 digits of precision.
  ! integer, parameter :: REALKIND = selected_real_kind( p = 12 )
  ! integer, parameter :: REALKIND = selected_real_kind( p = RPREC )
  
  ! Basic grid point resolution parameters
  
  integer,parameter::plon = PLON        ! number of longitudes
  integer,parameter::plat = PLAT        ! number of latitudes
  integer,parameter::plev = PLEV        ! number of vertical levels
  integer,parameter::pcnst = PCNST      ! number of advected constituents
  integer,parameter::pnats = PNATS      ! number of non-advected trace species
  integer,parameter::plevp = plev + 1
  integer,parameter::plevd = 2*plev
  integer,parameter::nxpt = 1           ! no. of pts outside active domain of interpolant
  integer,parameter::jintmx = 1         ! number of extra latitudes in polar region
  integer,parameter::plond = plon + 1 + 2*nxpt
  integer,parameter::platd = plat + 2*nxpt + 2*jintmx
  integer,parameter::i1 = nxpt + 1      ! model starting longitude (3-d)
  integer,parameter::j1 = jintmx + nxpt + 1 ! model starting latitude (3-d)
  integer,parameter::j1m = j1 - 1       ! model starting offset (3-d)
  integer,parameter::padv = 2           ! number of miscellaneous advected fields
  integer,parameter::mxdynflds = 42     ! maximum number of dynamics input fields

#ifdef DST
  integer,parameter::mxoutflds = 9*pcnst + pnats + 140 ! maximum number of history output fields
#else 
  integer,parameter::mxoutflds = 5*pcnst + pnats + 70 ! maximum number of history output fields
#endif /* not DST */
  
  logical,parameter::masterproc=.true.  ! for CCM compatibility.  true => shared memory code
  
#endif /* endif MATCH */
  
end module pmgrid
