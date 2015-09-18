! Program tests whether writing to unit 0 is supported in Fortran 90 (it should be).
! Result should be that the numeral "73" is written to unit 0 (normally stderr).
! Sun f90 causes core dump on Solaris 2.6
! SGI f90 works fine
! Cray f90 works fine
! Please contact Charlie Zender <zender@ncar.ucar.edu> (303) 497-1612 if you fix this.
! Usage:
! f90 -o bug bug.f90; ./bug
program bug
  write (6,'( &
       &i4.4,a1, & 
       &i2.2,a1 &
       &)') &
       2002,'a',98,'b'
end Program bug
