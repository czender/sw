program hello
  ! g95 hello.f90
  integer idx
  character(len=*),parameter::sng="foobar"
  idx=73
  write(6,*) "Hello, World! from hello.f90 program"
  write (6,'(a,i2)') 'DEBUG1 '//sng//' DEBUG2',idx
end program hello
