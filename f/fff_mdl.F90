! $Id$

! Purpose: Fortran features and regression tests

! Copyright (C) 1994--present Charlie Zender
! This software is distributed under the terms of the GNU General Public License
! See http://www.gnu.org/copyleft/gpl.html for full license text

! Usage:
! use fff_mdl ! [mdl] Fortran features and regression tests

module fff_mdl ! [mdl] Fortran features and regression tests
  implicit none
  private ! [stt] Symbols are private unless individually qualified as public
  public::date_time_tst ! [sbr] Test date and time routines
  public::gsl_tst ! [sbr] Test scientific library functions
  public::ieee_tst ! [sbr] Test IEEE error trapping
  public::omp_tst ! [sbr] Test OpenMP routines
  public::prc_tst ! [sbr] Test precision handling
  public::prn_tst ! [sbr] Test formatted print statements
  public::ptr_tst ! [sbr] Test Fortran pointers
  public::rbn_tst ! [sbr] Test rbn_vec() subroutine and subsidiaries
  public::sbr_foo ! [sbr] Foo
  public::syn_tst ! [sbr] Test Fortran syntax
  
contains 
  
  subroutine date_time_tst( & ! [sbr] Test date and time routines
       lcl_date_time) ! O [sng] Time formatted as Day Mth DD HH:MM:SS TZ YYYY
    ! Purpose: Test date and time routines
    use shr_kind_mod,only:r8=>shr_kind_r8,DBLKIND ! Precision r8, i8, ...
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    implicit none
    ! Parameters
    ! Commons
    ! Input
    ! Output 
    character(26),intent(out)::lcl_date_time ! Time formatted as Day Mth DD HH:MM:SS TZ YYYY
    ! Local workspace
    character(8) date
    character(10) time
    character(5) zone
    character(4) tz_sng ! [sng] Current time representation, e.g., GMT
    character(32) gmt_sng ! [sng] GMT formatted as Day Mth DD HH:MM:SS TZ YYYY
    integer time_val(8)
    real(DBLKIND) time_gmt ! [s] Seconds between 1969 and GMT of simulation
    ! Main code
    if (dbg_lvl >= dbg_sbr) write(6,"(a)") "Entering date_time_tst()"
    
    call date_and_time(date,time,zone,time_val) ! Fortran90 DATE_AND_TIME intrinsic
    write (6,"(a)") "Fortran90 intrinsic date_and_time() routine:"
    write (6,"(6a)") "Current date = ",date,", time = ",time,", zone = ",zone
    write (6,"(8(a,i4))") &
         "year = ",time_val(1), &
         ", month = ",time_val(2), &
         ", day = ",time_val(3), &
         ", dlt UTC = ",time_val(4), &
         ", hour = ",time_val(5), &
         ", minute = ",time_val(6), &
         ", second = ",time_val(7), &
         ", millisecond = ",time_val(8)
    write (lcl_date_time,"(i4.4,5(a1,i2.2),a1,i3.3)") &
         time_val(1),"/",time_val(2),"/",time_val(3)," ", &
         time_val(5),":",time_val(6),":",time_val(7),".",time_val(8)
    write (6,"(2a)") "Assembled into lcl_date_time = ",lcl_date_time
    
    tz_sng="GMT" ! [sng] Current time representation, e.g., GMT
    ! Answer for time_gmt=815079107.6 should be 1995/10/30 DOY 303 GMT 18:51:47.600000
    time_gmt=815079107.6_DBLKIND ! [s] Seconds between 1969 and GMT of simulation
    ! Create standard date-time string for given decimal time and time-zone
    call unix2gmt_sng(time_gmt,gmt_sng,tz_sng)
    write (6,"(a)") "Testing Fortran interface to C time routines: unix2gmt_sng()"
    write (6,"(a,es20.12,2a)") "time_gmt = ",time_gmt,", gmt_sng = ",gmt_sng
    write (6,"(a,f15.3)") "Compare to results of pure C: time_unix2gmt ",time_gmt
    
    if (dbg_lvl >= dbg_sbr) write(6,"(a)") "Exiting date_time_tst()"
    return
  end subroutine date_time_tst ! end date_time_tst()
  
  subroutine gsl_tst( & ! [sbr] Test scientific library functions
       cmp_prc_foo & ! [frc] Computational precision
       )
    ! Purpose: Test scientific library functions
    use shr_kind_mod,only:r8=>shr_kind_r8 ! Precision r8, i8, ...
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use erf_mdl,only:erf ! [mdl] Error functions erf(), erfc(), erfcx()
    use gmm_mdl,only:gamma ! [mdl] Gamma function gamma()
    implicit none
    ! Parameters
    ! Commons
    ! Input
    real(r8),intent(in)::cmp_prc_foo ! [frc] Computational precision
    ! Output 
    ! Local workspace
    real(r8) pi ! [frc] 3
    ! Externals
#if 0
    ! Deprecated code to access functions from external libraries
    real(r8),external::erf ! Error function (libspecfun)
    real(r8),external::gamma ! Gamma function (libspecfun)
    real(r8),external::gammq ! Incomplete gamma function (librecipes_f)
#endif /* !0 */
    ! Main code
    if (dbg_lvl >= dbg_sbr) write(6,"(a)") "Entering gsl_tst()"
    write (6,"(a)") "Testing Scientific library routines in gsl_tst()..."
    ! Initialize defaults
    pi=4.0_r8*atan(1.0_r8)        ! [frc] 3
    
    ! Test ambiguous-looking math syntax
    write (6,"(a)") "Testing ambiguous-looking math syntax:"
    ! csz 20061125 following line causes ICE on gfortran 4.1.2:
    !    write (6,"(a,es15.8)") "2.0_r8**(-0.0_r8) = ",2.0_r8**(-0.0_r8)
    write (6,"(a,es15.8)") "cmp_prc_foo = ",cmp_prc_foo
    write (6,"(a,es15.8,a,es15.8)") "1.0_r8/(1.0_r8*2.0_r8)*",cmp_prc_foo," = ",1.0_r8/(1.0_r8*2.0_r8)*cmp_prc_foo
    
    ! Test error function
    write (6,"(a)") "Testing error function erf():"
    write (6,"(3(a,es15.8))") "erf(",0.0_r8,") = ",erf(0.0_r8)," =?= ",0.0_r8
    write (6,"(3(a,es15.8))") "erf(",1.0_r8,") = ",erf(1.0_r8)," =?= ",8.42700793E-01_r8
    write (6,"(2(a,es15.8))") "erf(",cmp_prc_foo,") = ",erf(cmp_prc_foo)
    if (abs(0.8427_r8-erf(1.0_r8))/0.8427_r8 > 0.001_r8) error stop "gsl_tst() reports Error function error"
    if (erf(0.0_r8) /= 0.0_r8) error stop "gsl_tst() reports Error function error"
    
    ! Test gamma function
    write (6,"(a)") "Testing gamma function gamma():"
    write (6,"(3(a,es15.8))") "gamma(",0.5_r8,") = ",gamma(0.5_r8)," =?= ",sqrt(pi)
    write (6,"(2(a,es15.8))") "gamma(",cmp_prc_foo,") = ",gamma(cmp_prc_foo)
    if (abs(sqrt(pi)-gamma(0.5_r8))/sqrt(pi) > 1.0e-5_r8) &
         error stop "gsl_tst() reports Gamma function error"
    
#if 0
    ! Incomplete gamma function gammq() is only available in librecipes_f
    ! Deprecate to remove dependence of fff on librecipes_f
    ! fxm: Add incomplete gamma function gammq() to gmm_mdl
    ! Test incomplete gamma function
    write (6,"(a)") "Testing incomplete gamma function gammq():"
    write (6,"(3(a,es15.8))") "gammq(",1.0_r8,",",cmp_prc_foo,") = ",gammq(1.0_r8,cmp_prc_foo)
    write (6,"(3(a,es15.8))") "gammq(",1.0_r8,",",cmp_prc_foo,") = ",gammq(1.0_r8,cmp_prc_foo)
    if (gammq(73.0_r8,0.0_r8) /= 1.0_r8) error stop "Incomplete gamma function error"
    if (abs(exp(-1.0_r8)-1.0_r8*gammq(1.0_r8,1.0_r8))/exp(-1.0_r8) > 1.0e-5_r8) &
         error stop "gsl_tst() reports Incomplete gamma function error"
#endif /* !0 */
    
    if (dbg_lvl >= dbg_sbr) write(6,"(a)") "Exiting gsl_tst()"
    return
  end subroutine gsl_tst ! end gsl_tst()
  
  subroutine ieee_tst( & ! [sbr] Test IEEE error trapping
       cmp_prc_foo &            ! [frc] Computational precision
       )
    ! Purpose: Test IEEE error trapping
    use shr_kind_mod,only:r8=>shr_kind_r8 ! Precision r8, i8, ...
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    implicit none
    ! Parameters
    ! Commons
    ! Input
    real(r8),intent(in)::cmp_prc_foo ! [frc] Computational precision
    ! Output 
    ! Local workspace
    real(r8) one_over_cmp_prc_foo     ! [frc] 1.0/cmp_prc_foo
    ! Main code 
    if (dbg_lvl >= dbg_sbr) write(6,"(a)") "Entering ieee_tst()"
    write (6,"(a)") "Testing IEEE error trapping..."
    ! Initialize defaults
    ! Generate SIGFPE
    write (6,"(a,es9.2,a)") "About to attempt to divide by ",cmp_prc_foo," ..."
    one_over_cmp_prc_foo=1.0_r8/cmp_prc_foo ! [frc] 1.0_r8/cmp_prc_foo
    write (6,"(2(a,es9.2))") "1.0_r8/",cmp_prc_foo," = ",one_over_cmp_prc_foo
    
    if (dbg_lvl >= dbg_sbr) write(6,"(a)") "Exiting ieee_tst()"
    return
  end subroutine ieee_tst ! end ieee_tst()
  
  subroutine omp_tst( & ! [sbr] Test OpenMP routines
       cmp_prc_foo, & ! [frc] Computational precision
       int_foo & ! [nbr] Integer
       )
    ! Purpose: Test OpenMP routines
    ! OpenMP API: 
    ! http://www.openmp.org/specs
    ! file:/usr/local/pgi/doc/pgiws_ug/pgi31u.htm
    ! Sample OpenMP programs:
    ! Fortran: babyblue.ucar.edu:/usr/local/examples/openmp
    ! C: utefe.ucar.edu:~rosinski/timing
    use shr_kind_mod,only:r8=>shr_kind_r8 ! Precision r8, i8, ...
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use sng_mdl,only:ftn_strlen ! [mdl] String manipulation
    implicit none
    ! Parameters
    ! Commons
    ! Input
    integer int_foo           ! [nbr] Integer
    real(r8),intent(in)::cmp_prc_foo  ! [frc] Computational precision
    ! Output 
    ! Allocatable variables
    integer,dimension(:),allocatable::thr_idx ! [idx] OpenMP thread index
    ! Local workspace
    integer idx ! [idx] Counting index
    integer rcd ! [rcd] Return success code
#ifdef _OPENMP
    integer thr_nbr ! [nbr] OpenMP number of threads
    integer omp_get_num_threads ! [nbr] OpenMP number of threads
    integer omp_get_thread_num ! [idx] OpenMP thread index
#endif /* !_OPENMP */
    ! Main code
    if (dbg_lvl >= dbg_sbr) write(6,"(a)") "Entering omp_tst()"
    
    rcd=0 ! [rcd] Return success code
    write (6,"(a)") "Testing OpenMP implementation in omp_tst()..."
    !$omp parallel default(none) private(idx) shared(int_foo,cmp_prc_foo,rcd,thr_nbr,thr_idx,prg_nm)
#ifdef _OPENMP /* OpenMP-compliant compilers define _OPENMP=YYYYMM = year and month of OpenMP specification */
    !$omp single
    thr_nbr=omp_get_num_threads() ! [nbr] OpenMP number of threads
    ! Allocate space for dynamic arrays
    allocate(thr_idx(thr_nbr),stat=rcd) ! [idx] OpenMP thread index
    ! Initialize indices for all threads
    thr_idx(:)=-1
    write (6,"(2a,i3,a)") prg_nm(1:ftn_strlen(prg_nm)),": INFO OpenMP threading with ",thr_nbr," threads"
    !$omp end single
    thr_idx(omp_get_thread_num())=omp_get_thread_num()
#else /* not _OPENMP */
    write (6,"(2a)") prg_nm(1:ftn_strlen(prg_nm)),": INFO OpenMP was not enabled during compilation"
    write (6,"(2a)") prg_nm(1:ftn_strlen(prg_nm)),": INFO Not attempting OpenMP threading"
#endif /* not _OPENMP */
    !$omp do
    do idx=0,int_foo
       call syn_tst(cmp_prc_foo)
    enddo                     ! end loop over idx
    !$omp end do
#ifdef _OPENMP /* OpenMP-compliant compilers define _OPENMP=YYYYMM = year and month of OpenMP specification */
    !$omp single
    write (6,"(a)") "Values != -1 prove that parallel loop was executed:"
    write (6,"(a3,1x,a3)") "thr","idx"
    do idx=0,thr_nbr-1
       write (6,"(i3,1x,i3)") idx,thr_idx(idx)
    enddo                     ! end loop over idx
    !$omp end single
#endif /* not _OPENMP */
    !$omp end parallel
    !$omp parallel default(none) shared(prg_nm)
    !$omp single
    !$write (6,"(a,a36,i3,a8)") prg_nm(1:ftn_strlen(prg_nm)),": INFO OpenMP multi-threading using ",omp_get_num_threads()," threads"
    !$omp end single
    !$omp end parallel
    
    ! De-allocate dynamic variables
    if (allocated(thr_idx)) deallocate(thr_idx,stat=rcd) ! [idx] OpenMP thread index
    if(rcd /= 0) error stop "deallocate() failed for thr_idx"
    
    if (dbg_lvl >= dbg_sbr) write(6,"(a)") "Exiting omp_tst()"
    return
  end subroutine omp_tst ! end omp_tst()
  
  subroutine prc_tst( & ! [sbr] Test precision handling
       cmp_prc_foo, & ! [frc] Computational precision
       dbl_foo, & ! [frc] Double
       flt_foo, & ! [frc] Float
       int_foo & ! [nbr] Integer
       )
    ! Purpose: Test precision handling
    ! prc_tst does not use "precision" module but rather tests precision handling
    ! prc_tst names precisions with the sensible convention that NCAR should adopt:
    ! r8 is 8 byte double precision
    ! r4 is 4 byte single precision
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    implicit none
    integer,parameter::r4=selected_real_kind(p=6) ! r4: 4B (C float) default, 8B (C double) possible
    integer,parameter::r8=selected_real_kind(p=12) ! r8: 8B (C double) default, 4B (C float) possible
    integer,parameter::i8=selected_int_kind(13) ! i8: 8B (C long long) default, 4B (C int) possible
    integer,parameter::DBLKIND=selected_real_kind(p=12) ! DBLKIND: Fixed 8B (double precision) always
    integer,parameter::r_ntv=kind(1.0) ! Native real kind
    ! Parameters
    ! Commons
    ! Input
#ifdef PRC_DBL
    real(r8),intent(in)::cmp_prc_foo ! [frc] Computational precision
#endif /* !PRC_DBL */
#ifdef PRC_FLT
    real(r4),intent(in)::cmp_prc_foo ! [frc] Computational precision
#endif /* !PRC_FLT */
#ifdef PRC_NTV
    real(r_ntv),intent(in)::cmp_prc_foo ! [frc] Computational precision
#endif /* !PRC_NTV */
    real(DBLKIND),intent(in)::dbl_foo ! [frc] Double
    real(r4),intent(in)::flt_foo ! [frc] Float
    integer,intent(in)::int_foo ! [nbr] Integer
    ! Output 
    ! Local workspace
    real(DBLKIND)::dblkind_var=0.0_DBLKIND ! [frc] DBLKIND precision real variable
    real(r8)::r8_var=0.0_r8 ! [frc] r8 precision real variable
    real(r4)::r4_var=0.0_r4 ! [frc] r4 precision real variable
    real::r_var=0.0 ! [frc] Default real variable
    integer(kind=1)::i1_var=0 ! [nbr] i1 kind integer variable
    integer(kind=2)::i2_var=0 ! [nbr] i2 kind integer variable
    integer(kind=4)::i4_var=0 ! [nbr] i4 kind integer variable
    integer(i8)::i8_var=0 ! [nbr] i8 kind integer variable
    integer::i_var=0 ! [nbr] Default integer variable
    ! Main code
    if (dbg_lvl >= dbg_sbr) write(6,"(a)") "Entering prc_tst()"
    write (6,"(a)") "Testing precision and range of variables..."
    r4_var=flt_foo+r4_var ! [frc] Float CEWI
    ! Initialize defaults
    write (6,230) &
         "Reals, format = Variable_name(kind_specification):", &
         "r_var: precision()=",precision(r_var),", range()=",range(r_var),", kind()=",kind(r_var), &
         "r4_var(r4): precision()=",precision(r4_var),", range()=",range(r4_var),", kind()=",kind(r4_var), &
         "r8_var(r8): precision()=",precision(r8_var),", range()=",range(r8_var),", kind()=",kind(r8_var), &
         "dblkind_var(DBLKIND): precision()=",precision(dblkind_var),", range()=",range(dblkind_var), & 
         ", kind()=",kind(dblkind_var), &
#ifdef PRC_DBL
         "cmp_prc_foo(r8): precision()=",precision(cmp_prc_foo),", range()=",range(cmp_prc_foo),", kind()=",kind(cmp_prc_foo), &
#endif /* !PRC_DBL */
#ifdef PRC_FLT
         "cmp_prc_foo(r4): precision()=",precision(cmp_prc_foo),", range()=",range(cmp_prc_foo),", kind()=",kind(cmp_prc_foo), &
#endif /* !PRC_FLT */
         "dbl_foo(DBLKIND): precision()=",precision(dbl_foo),", range()=",range(dbl_foo),", kind()=",kind(dbl_foo), &
         "Integers, format=Variable_name(kind_specification):", &
         "i_var: range()=",range(i_var),", kind()=",kind(i_var), &
         "i1_var(kind=1): range()=",range(i1_var),", kind()=",kind(i1_var), &
         "i2_var(kind=2): range()=",range(i2_var),", kind()=",kind(i2_var), &
         "i4_var(kind=4): range()=",range(i4_var),", kind()=",kind(i4_var), &
         "i8_var(kind=8): range()=",range(i8_var),", kind()=",kind(i8_var), &
         "int_foo: range()=",range(int_foo),", kind()=",kind(int_foo)
230 format( &
         a,/, &
         6(a,i2,a,i3,a,i2,/), &
         a,/, &
         2(a,i3,a,i2,/) &
         )
    
    if (dbg_lvl >= dbg_sbr) write(6,"(a)") "Exiting prc_tst()"
    return
  end subroutine prc_tst ! end prc_tst()
  
  subroutine prn_tst( & ! [sbr] Test formatted print statements
       cmp_prc_foo &            ! [frc] Computational precision
       )
    ! Purpose: Test formatted print statements
    use shr_kind_mod,only:r8=>shr_kind_r8 ! Precision r8, i8, ...
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    implicit none
    ! Parameters
    ! Commons
    ! Input
    real(r8),intent(in)::cmp_prc_foo ! [frc] Computational precision
    ! Output 
    ! Local workspace
    ! Main code
    if (dbg_lvl >= dbg_sbr) write(6,"(a)") "Entering prn_tst()"
    write (6,"(a)") "Testing formatted print statements..."
    ! Initialize defaults
    write (6,*) "write(6,*) produces ",cmp_prc_foo
    write (6,230) &
         "e8.1 = ",cmp_prc_foo, &
         "es8.1 = ",cmp_prc_foo, &
         "en8.1 = ",cmp_prc_foo, &
         "f9.5 = ",cmp_prc_foo, &
         "f8.5 = ",cmp_prc_foo, &
         "g8.1 = ",cmp_prc_foo
230 format( &
         a15,e8.1,/, &
         a15,es8.1,/, &
         a15,en8.1,/, &
         a15,f9.5,/, &
         a15,f8.5,/, &
         a15,g8.1,/ &
         )
    
    if (dbg_lvl >= dbg_sbr) write(6,"(a)") "Exiting prn_tst()"
    return
  end subroutine prn_tst ! end prn_tst()
  
  subroutine ptr_tst( & ! [sbr] Test Fortran pointers
       one_dmn_var, & ! [frc] 
       one_dmn_var_ptr & ! [ptr] 
       )
    ! Purpose: Test Fortran pointers
    use shr_kind_mod,only:r8=>shr_kind_r8 ! Precision r8, i8, ...
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    implicit none
    ! Parameters
    ! Commons
    ! Input
    real(r8),dimension(:),target,intent(in)::one_dmn_var
    real(r8),dimension(:),pointer::one_dmn_var_ptr
    ! Output 
    ! Local workspace
    ! Main code
    if (dbg_lvl >= dbg_sbr) write(6,"(a)") "Entering ptr_tst()"
    write (6,"(a)") "Testing Fortran pointers..."
    ! Associate pointer with variable
    one_dmn_var_ptr=>one_dmn_var
    if (dbg_lvl >= dbg_sbr) write(6,"(a)") "Exiting ptr_tst()"
    return
  end subroutine ptr_tst ! end ptr_tst()
  
  subroutine rbn_tst( & ! [sbr] Test rbn_vec() subroutine and subsidiaries
       xtr_typ_LHS, & ! I [msk] Extrapolation flags for LHS
       xtr_typ_RHS) ! I [msk] Extrapolation flags for RHS
    ! Purpose: Test rbn_vec() subroutine and subsidiaries
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use shr_kind_mod,only:r8=>shr_kind_r8 ! Precision r8, i8, ...
    use vec_mdl ! [mdl] Vector manipulation, interpolation, rebinning
    implicit none
    ! Parameters
    integer,parameter::crd_in_nbr=5
    integer,parameter::grd_in_nbr=crd_in_nbr+1
    
    integer,parameter::crd_out_nbr=2
    integer,parameter::grd_out_nbr=crd_out_nbr+1
    
    ! Input 
    integer,intent(in)::xtr_typ_LHS ! I [msk] Extrapolation flags for LHS
    integer,intent(in)::xtr_typ_RHS ! I [msk] Extrapolation flags for RHS
    
    ! Local workspace
    integer idx
    
    real(r8) crd_in_max(crd_in_nbr)
    real(r8) crd_in_min(crd_in_nbr)
    real(r8) crd_out_max(crd_out_nbr)
    real(r8) crd_out_min(crd_out_nbr)
    real(r8) dat_in_idx(crd_in_nbr)
!    real(r8) dat_in_neg_one(crd_in_nbr)
!    real(r8) dat_in_one(crd_in_nbr)
    real(r8) dat_out(crd_out_nbr)
    real(r8) grd_in(grd_in_nbr)
    real(r8) grd_out(grd_out_nbr)
    
    ! Main code
    if (dbg_lvl >= dbg_sbr) write(6,"(a)") "Entering rbn_tst()"
    do idx=1,crd_in_nbr
       crd_in_min(idx)=idx-1.0_r8
       grd_in(idx)=crd_in_min(idx)
!       dat_in_one(idx)=1.0_r8
       dat_in_idx(idx)=idx
!       dat_in_neg_one(idx)=-1.0_r8
    end do                    ! end loop over crd_in
    do idx=1,crd_in_nbr-1
       crd_in_max(idx)=crd_in_min(idx+1)
    end do                    ! end loop over crd_in
    crd_in_max(crd_in_nbr)=crd_in_nbr
    grd_in(crd_in_nbr+1)=crd_in_max(crd_in_nbr)
    
    do idx=1,crd_out_nbr
       crd_out_min(idx)=idx-0.5_r8
       grd_out(idx)=crd_out_min(idx)
    end do                    ! end loop over crd_out
    do idx=1,crd_out_nbr-1
       crd_out_max(idx)=crd_out_min(idx+1)
    end do                    ! end loop over crd_in
    crd_out_max(crd_out_nbr)=crd_out_min(crd_out_nbr)+1.0_r8
    grd_out(crd_out_nbr+1)=crd_out_max(crd_out_nbr)
    
    call rbn_vec(crd_in_nbr,grd_in,dat_in_idx, &
         crd_out_nbr,grd_out,dat_out, &
         xtr_typ_LHS,xtr_typ_RHS)
    
    if (dbg_lvl >= dbg_sbr) write(6,"(a)") "Exiting rbn_tst()"
    return
  end subroutine rbn_tst ! end rbn_tst()
  
  subroutine sbr_foo( & ! [sbr] Foo
       arr, & ! []
       nbr_lmn) ! []
    ! Purpose: Foo
    use shr_kind_mod,only:r8=>shr_kind_r8 ! Precision r8, i8, ...
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    implicit none
    ! Parameters
    ! Commons
    ! Input
    integer nbr_lmn
    real(r8) arr(nbr_lmn)
    ! Input/Output
    ! Output
    ! Local workspace
    integer lmn
    ! Main code
    if (dbg_lvl >= dbg_sbr) write(6,"(a)") "Entering sbr_foo()"
    write (6,"(a)") "From subroutine"
    do lmn=1,nbr_lmn
       write (6,"(i4,1(2x,es14.7))") lmn,arr(lmn)
    end do
    if (dbg_lvl >= dbg_sbr) write(6,"(a)") "Exiting sbr_foo()"
    return
  end subroutine sbr_foo
  
  subroutine syn_tst( & ! [sbr] Test Fortran syntax
       cmp_prc_foo &           ! [frc] Computational precision
       )
    ! Purpose: Test Fortran syntax
    use shr_kind_mod,only:r8=>shr_kind_r8 ! Precision r8, i8, ...
    use dbg_mdl ! [mdl] Debugging constants, prg_nm, dbg_lvl
    use sng_mdl,only:ftn_strlen ! [mdl] String manipulation
    implicit none
    ! Parameters
    ! Commons
    ! Input
    real(r8),intent(in)::cmp_prc_foo  ! [frc] Computational precision
    ! Output 
    ! Local workspace
    real(r8) cmp_prc_foo_tmp          ! [frc] Temporary variable
    ! Main code
    cmp_prc_foo_tmp=cmp_prc_foo       ! [frc] Temporary variable CEWI
    cmp_prc_foo_tmp=0.0_r8+cmp_prc_foo_tmp       ! [frc] Temporary variable CEWI
    if (dbg_lvl >= dbg_sbr) write(6,"(a)") "Entering syn_tst()"
    write (6,"(a)") "Testing Fortran syntax..."
    ! Initialize defaults
    ! if construct
    if (.true.) write (6,"(a,a)") prg_nm(1:ftn_strlen(prg_nm)),": if construct works as expected"
    ! if-then construct
    if (.true.) then
       write (6,"(a,a)") prg_nm(1:ftn_strlen(prg_nm)),": if-then construct works as expected"
    endif                     ! endif
    ! if-then-else construct
    if (.false.) then
       write (6,"(a,a)") prg_nm(1:ftn_strlen(prg_nm)),": ERROR: if-then-else construct is broken"
       error stop
    else
       write (6,"(a,a)") prg_nm(1:ftn_strlen(prg_nm)),": if-then-else construct works as expected"
    endif                     ! endif
    
    if (dbg_lvl >= dbg_sbr) write(6,"(a)") "Exiting syn_tst()"
    return
  end subroutine syn_tst ! end syn_tst()
  
end module fff_mdl ! [mdl] Fortran features and regression tests
