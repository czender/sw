// $Id$ 

// Purpose: OpenMP (OMP) shared memory parallelism (SMP) utilities

/* Copyright (C) 1997--present Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

#include <omp.hh> // OpenMP utilities

void
cxx11_hello() // [fnc] Test CXX-11 multi-threading environment
{
  // Local
  const std::string sbr_nm("cxx11_hello"); // [sng] Name of subroutine

  std::cout << sbr_nm << " says Hello from thread " << std::this_thread::get_id() << std::endl;
} // end cxx11_hello()

int // O [nbr] Thread number in parallel regions
cxx11_ini // [fnc] Set up CXX-11 multi-threading environment
(const int thr_nbr) // I [nbr] User-requested thread number
{
  // Output
  int thr_nbr_act=0; // O [nbr] Thread number in parallel regions
  // Local
  const std::string sbr_nm("cxx11_ini"); // [sng] Name of subroutine

  std::cout << "Hello from thread " << std::this_thread::get_id() << std::endl;

  return thr_nbr_act; // O [nbr] Thread number in parallel regions
} // end cxx11_ini()

int // O [nbr] Thread number in parallel regions
openmp_ini // [fnc] Set up OpenMP multi-threading environment
(const int thr_nbr) // I [nbr] User-requested thread number
{
  /* Purpose: Initialize OpenMP multi-threading environment
     Honor user-requested thread number, balance against known code efficiency,
     print diagnostics */
  // Output
  int thr_nbr_act=0; // O [nbr] Thread number in parallel regions
  // Local
  const std::string sbr_nm("openmp_ini"); // [sng] Name of subroutine

  // using namespace std; // std is default namespace

#ifdef _OPENMP
  /* System allocates OMP_NUM_THREADS if possible
     ncwa is I/O bottlenecked beyond about thr_nbr_max_fsh=4 threads
     If OMP_NUM_THREADS > 4 then NCO will not be using threads efficiently
     Strategy: Determine maximum number of threads system will allocate (thr_nbr_max)
     Reduce maximum number of threads available to system to thr_nbr_max_fsh
     Play nice: Set dynamic threading so that system can make efficiency decisions
     When dynamic threads are set, then system will never allocate more than thr_nbr_max_fsh
  */
  const int prc_nbr_max(omp_get_num_procs()); // [nbr] Maximum number of processors available
  const int thr_nbr_max(omp_get_max_threads()); // [nbr] Maximum number of threads system allows
  const std::string nvr_OMP_NUM_THREADS((std::getenv("OMP_NUM_THREADS")) ? std::getenv("OMP_NUM_THREADS") : ""); // [sng] Environment variable OMP_NUM_THREADS

  bool USR_SPC_THR_RQS=false;

  int dyn_thr(1); // [flg] Allow system to dynamically set number of threads 
  int ntg_OMP_NUM_THREADS=CEWI_int; // [nbr] OMP_NUM_THREADS environment variable
  int thr_nbr_max_fsh=32; // [nbr] Maximum number of threads program can use efficiently 
  int thr_nbr_rqs=CEWI_int; // [nbr] Number of threads to request 

  if(thr_nbr < 0) err_prn(prg_nm_get(),sbr_nm,"User-requested thread number = "+nbr2sng(thr_nbr)+" is less than zero");

  if(thr_nbr > 0) USR_SPC_THR_RQS=true;

  if(dbg_lvl_get() > 2){
    if(nvr_OMP_NUM_THREADS != "") ntg_OMP_NUM_THREADS=static_cast<int>(std::strtol(nvr_OMP_NUM_THREADS.c_str(),(char **)NULL,10));
    std::cout << prg_nm_get() << ": INFO Environment variable OMP_NUM_THREADS ";
    if(ntg_OMP_NUM_THREADS > 0) std::cout << "= " << ntg_OMP_NUM_THREADS; else std::cout << "does not exist";
    std::cout << std::endl;
    std::cout << prg_nm_get() << ": INFO Maximum number of threads system allows is " << thr_nbr_max << std::endl;
    std::cout << prg_nm_get() << ": INFO Number of processors available is " << prc_nbr_max << std::endl;
  } // endif dbg

  if(USR_SPC_THR_RQS){
    // Honor user-specified thread request...
    thr_nbr_rqs=thr_nbr; // [nbr] Number of threads to request
    // ...if possible...
    if(dbg_lvl_get() > 2) std::cout << prg_nm_get() << ": INFO User requested " << thr_nbr << " threads" << std::endl;
    if(thr_nbr > thr_nbr_max){
      std::cout << prg_nm_get() << ": WARNING Reducing user-requested thread number = " << thr_nbr << " to maximum thread number allowed = " << thr_nbr_max << std::endl;
      thr_nbr_rqs=thr_nbr_max; // [nbr] Number of threads to request 
    } // endif
  }else{
    // Automatic thread allocation algorithm 

    // Request maximum number of threads permitted 
    thr_nbr_rqs=thr_nbr_max; // [nbr] Number of threads to request 

    // Disable threading on per-program basis to play nicely with others 
    /* ncrcat is extremely I/O intensive 
       Maximum efficiency when one thread reads from input file while other writes to output file */
    //    if(prg_nm_get() == "ccc") thr_nbr_max_fsh=2;
    if(!std::strcmp(prg_nm_get(),"ccc")) thr_nbr_max_fsh=2;
    
    // Play nice with others 
    (void)omp_set_dynamic(dyn_thr); // [flg] Allow system to dynamically set number of threads
    if(dbg_lvl_get() > 0) std::cout << prg_nm_get() << ": INFO " << (dyn_thr ? "Allowing" : "Not allowing") << " OS to utilize dynamic threading" << std::endl;
    dyn_thr=omp_get_dynamic(); // [flg] Allow system to dynamically set number of threads
    if(dbg_lvl_get() > 0) std::cout << prg_nm_get() << ": INFO System will " << (dyn_thr ? "" : "not ") << "utilize dynamic threading" << std::endl;

    // Apply program/system limitations
    if(thr_nbr_max > thr_nbr_max_fsh){
      if(dbg_lvl_get() > 0) std::cout << prg_nm_get() << ": INFO Reducing default thread number from " << thr_nbr_max << " to " << thr_nbr_max_fsh << ", a \"play-nice\" number set in openmp_ini()" << std::endl;
      thr_nbr_rqs=thr_nbr_max_fsh; // [nbr] Number of threads to request
    } // endif
  } // endif

  // Set thread number 
  if(omp_in_parallel()){
    err_prn(prg_nm_get(),sbr_nm,"Attempted to set thread number from within parallel region");
  }else{
    omp_set_num_threads(thr_nbr_rqs); 
    if(dbg_lvl_get() > 0) std::cout << prg_nm_get() << ": INFO " << sbr_nm << "() requested " << thr_nbr_rqs << " threads from system" << std::endl;
  } // end error

#pragma omp parallel default(none) shared(std::cout,thr_nbr_act)
  { // begin OpenMP parallel
#pragma omp single nowait
    { // begin OpenMP single
      thr_nbr_act=omp_get_num_threads(); // [nbr] Thread number in parallel regions
      if(dbg_lvl_get() > 0) std::cout << prg_nm_get() << ": INFO Parallel regions spawn teams of " << thr_nbr_act << " threads" << std::endl;
    } // end OpenMP single
  } // end OpenMP parallel
#else // !_OPENMP
  thr_nbr_act+=0*thr_nbr; // CEWI
  if(dbg_lvl_get() > 0) std::cout << prg_nm_get() << ": INFO Not attempting OpenMP threading" << std::endl;
#endif // !_OPENMP
  
  return thr_nbr_act; // O [nbr] Thread number in parallel regions
} // end openmp_ini()
