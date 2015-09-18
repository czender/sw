// $Id$ 

// Purpose: Vector utilities for C++ programs

/* Copyright (C) 1997--2014 Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

/* fxm: vec routines may segfault if compiled with bounds checking/electric fence
   Reason appears to be that ntp_vec() and rbn_vec() _appear_ to access illegal memory
   _appear_ because I checked these routines so many times, I am sure there is no bug
   Three such potential segfaults are marked in the code */

/* ntp_vec() and rbn_vec() contain some while() conditions that depend on previous conditions being true 
   Compiler optimization level determines whether subsequent conditions are tested when one is false
   If previous condition is false, optimized code never checks subsequent conditions, so everything works
   If previous condition is false, unoptimized code checks subsequent conditions which may be illegal */

#include <vec.hh> // Vector functions ntp_vec(), rbn_vec()

void vec_tst(Xtr_cls xtr_LHS,Xtr_cls xtr_RHS)
{
  // Purpose: Test ntp_vec() and rbn_vec()
  int rcd(0); // [enm] Return success code
  const long crd_in_nbr(5);
  const long grd_in_nbr=crd_in_nbr+1;
  const long crd_out_nbr(7);
  const long grd_out_nbr=crd_out_nbr+1;

  long idx;
  long out_nbr;

  prc_cmp *crd_in=new prc_cmp[crd_in_nbr];
  prc_cmp *crd_in_max=new prc_cmp[crd_in_nbr];
  prc_cmp *crd_in_min=new prc_cmp[crd_in_nbr];
  prc_cmp *crd_out=new prc_cmp[crd_out_nbr];
  prc_cmp *crd_out_max=new prc_cmp[crd_out_nbr];
  prc_cmp *crd_out_min=new prc_cmp[crd_out_nbr];
  prc_cmp *dat_in_idx=new prc_cmp[crd_in_nbr];
  prc_cmp *dat_in_neg_one=new prc_cmp[crd_in_nbr];
  prc_cmp *dat_in_one=new prc_cmp[crd_in_nbr];
  prc_cmp *dat_out_crd_out=new prc_cmp[crd_out_nbr];
  prc_cmp *grd_in=new prc_cmp[grd_in_nbr];
  prc_cmp *grd_out=new prc_cmp[grd_out_nbr];
  prc_cmp *dat_out_grd_in=new prc_cmp[grd_in_nbr]; // Tests output data on input grd
  prc_cmp *dat_out_crd_in=new prc_cmp[crd_in_nbr]; // Tests output data on input crd

  // Main code
  std::string sbr_nm("vec_tst");
  unsigned short dbg_lvl=dbg_lvl_get(); // Debugging level
  if(dbg_lvl >= dbg_sbr) dbg_prn(sbr_nm,"Entering...");

  // Define input grid and data
  for(idx=0;idx<crd_in_nbr;idx++){
    crd_in_min[idx]=idx;
    grd_in[idx]=crd_in_min[idx];
    dat_in_one[idx]=1.0;
    dat_in_idx[idx]=idx;
    dat_in_neg_one[idx]=-1.0;
  } // end loop over crd_in
  for(idx=0;idx<crd_in_nbr-1;idx++){ // NB: Loop ends at crd_in_nbr-2
    crd_in_max[idx]=crd_in_min[idx+1];
  } // end loop over crd_in
  if (crd_in_nbr > 1) crd_in_max[crd_in_nbr-1]=2*crd_in_max[crd_in_nbr-2]-crd_in_min[crd_in_nbr-2]; else crd_in_max[0]=crd_in_min[0]+1.0;
  grd_in[crd_in_nbr]=crd_in_max[crd_in_nbr-1];
  for(idx=0;idx<crd_in_nbr;idx++){
    crd_in[idx]=0.5*(crd_in_min[idx]+crd_in_max[idx]);
  } // end loop over crd_in

  for(idx=0;idx<crd_out_nbr;idx++){
    crd_out_min[idx]=idx-1.5;
    grd_out[idx]=crd_out_min[idx];
  } // end loop over crd_out
  for(idx=0;idx<crd_out_nbr-1;idx++){
    crd_out_max[idx]=crd_out_min[idx+1];
  } // end loop over crd_out
  if (crd_out_nbr > 1) crd_out_max[crd_out_nbr-1]=2*crd_out_max[crd_out_nbr-2]-crd_out_min[crd_out_nbr-2]; else crd_out_max[0]=crd_out_min[0]+1.0;
  grd_out[crd_out_nbr]=crd_out_max[crd_out_nbr-1];
  for(idx=0;idx<crd_out_nbr;idx++){
    crd_out[idx]=0.5*(crd_out_min[idx]+crd_out_max[idx]);
  } // end loop over crd_out
      
  if(dbg_lvl == dbg_crr){
    // Print extrapolation flags
    err_prn(sbr_nm,"xtr_LHS = "+xtr_LHS.xtr_sng);
    err_prn(sbr_nm,"xtr_RHS = "+xtr_RHS.xtr_sng);
  } // endif dbg

  // Pick input array to test
  prc_cmp *dat_in(dat_in_idx);
  // Pick output array to hold answer
  prc_cmp *dat_out(dat_out_crd_out);
  out_nbr=crd_out_nbr;
  // prc_cmp *dat_out(dat_out_crd_in);
  //  out_nbr=crd_in_nbr;
  // prc_cmp *dat_out(dat_out_grd_in);
  // out_nbr=grd_in_nbr;

  if(dbg_lvl == dbg_crr){
    // Print input coordinates
    err_prn(sbr_nm,"Input crd_in");
    std::cerr << "idx\t" << "crd_in\t" << "crd_min\t" << "crd_max\t" << "dat_in\t" << std::endl;
    for(idx=0;idx<crd_in_nbr;idx++) std::cerr << idx << "\t" << crd_in[idx] << "\t" << crd_in_min[idx] << "\t" << crd_in_max[idx] << "\t" << dat_in[idx] << std::endl;
  } // endif dbg

  // Interpolate identity
  // rcd+=ntp_vec(crd_in_nbr,crd_in,dat_in,out_nbr,crd_in,dat_out,xtr_LHS,xtr_RHS);
  // Interpolate from coordinate to own grid
  // rcd+=ntp_vec(crd_in_nbr,crd_in,dat_in,out_nbr,grd_in,dat_out,xtr_LHS,xtr_RHS);
  // Interpolate generic arrays
  rcd+=ntp_vec(crd_in_nbr,crd_in,dat_in,out_nbr,crd_out,dat_out,xtr_LHS,xtr_RHS);
  // Rebin identity
  // rbn_vec(crd_in_nbr,grd_in,dat_in,out_nbr,grd_in,dat_out,xtr_LHS,xtr_RHS);
  // Rebin generic arrays
  // rbn_vec(crd_in_nbr,grd_in,dat_in,out_nbr,grd_out,dat_out,xtr_LHS,xtr_RHS);

  //  long out_nbr=sizeof(dat_out)/sizeof(dat_out[0]);
  if(dbg_lvl == dbg_crr){
    // Print results
    err_prn(sbr_nm,"Coordinate output");
    std::cerr << "idx\t" << "crd_out\t" << "crd_min\t" << "crd_max\t" << "dat_out\t" << std::endl;
    for(idx=0;idx<out_nbr;idx++) std::cerr << idx << "\t" << crd_out[idx] << "\t" << crd_out_min[idx] << "\t" << crd_out_max[idx] << "\t" << dat_out[idx] << std::endl;
  } // endif dbg

  // Free dynamic memory
  delete []crd_in;
  delete []crd_in_max;
  delete []crd_in_min;
  delete []crd_out;
  delete []crd_out_max;
  delete []crd_out_min;
  delete []dat_in_idx;
  delete []dat_in_neg_one;
  delete []dat_in_one;
  delete []dat_out_crd_out;
  delete []dat_out_grd_in;
  delete []dat_out_crd_in;
  delete []grd_in;
  delete []grd_out;

  if(dbg_lvl >= dbg_sbr) dbg_prn(sbr_nm,"Exiting...");
} // end vec_tst()

int // [enm] Return success code
rbn_vec // [fnc] Rebin a vector from one grid onto another grid
(const long in_nbr, // I [nbr] Input grid size
 prc_cmp *grd_in, // I [crd] Input grid
 prc_cmp *dat_in, // I [frc] Input data
 const long out_nbr, // I [nbr] Output grid size
 prc_cmp *grd_out, // I [crd] Output grid
 prc_cmp *dat_out, // O [frc] Output (re-binned) data
 const Xtr_cls &xtr_LHS, // I [enm] LHS extrapolation flags
 const Xtr_cls &xtr_RHS) // I [enm] RHS extrapolation flags
{
  /* Purpose: Rebin a vector from one grid onto another grid
     Input and output coordinate vectors should be monotonic 
     Input and output vectors are specified with grid arrays definining interface coordinates
     
     grd_in defines the coordinate grid associated with the values dat_in
     grd_out defines a new (user) coordinate grid onto which dat_in is interpolated
     grd_in and grd_out may be reversed (twice) and so cannot be declared intent::in or const
     
     Output (interpolated) values are set so that output data have the same integral properties as input data
     In plain English this means the interpolated values are average values in some neighborhood of the specified grid point
     Use this option to interpolate noisy input data, e.g., solar spectra or NO2 absorption cross sections, to more regular grids */
  
  // Output
  int rcd(0); // [enm] Return success code
  // Local Workspace
  long *brk_lft_idx=new long[out_nbr+1];
  long *brk_rgt_idx=new long[out_nbr+1];
  long grd_out_idx;
  long in_idx;
  long out_idx;
  long out_idxp1;
  long trp_idx;
  long trp_srt_idx;
  long trp_nbr;
  
  prc_cmp *crd_in=new prc_cmp[in_nbr]; // Input coordinates (midpoints of grid)
  prc_cmp *crd_out=new prc_cmp[out_nbr]; // Output coordinates (midpoints of grid)
  prc_cmp *dat_ntp_in=new prc_cmp[in_nbr+1]; // Linearly interpolated data on grd_in grid
  prc_cmp *dat_ntp_out=new prc_cmp[out_nbr+1]; // Linearly interpolated data on grd_out grid
  prc_cmp *trp_area=new prc_cmp[in_nbr]; // Area of each trapezoid defined by input data
  
  prc_cmp crd_dlt_lft;
  prc_cmp crd_dlt_rgt;
  prc_cmp crd_dlt_ttl;
  prc_cmp ovl_avg; // Overlap average
  prc_cmp ovl_frc; // Overlap fraction
  prc_cmp trp_area_lft;
  prc_cmp trp_area_rgt;
  prc_cmp trp_area_ttl;
  
  // Main Code
  std::string sbr_nm("rbn_vec");
  unsigned short dbg_lvl=dbg_lvl_get();
  if(dbg_lvl >= dbg_sbr) dbg_prn(sbr_nm,"Entering...");
  
  // Vet input
  if(in_nbr < 1) err_prn(sbr_nm,"in_nbr < 1");
  if(out_nbr < 1) err_prn(sbr_nm,"out_nbr < 1");
  
  // Initialize default values
  const long grd_in_nbr(in_nbr+1);
  const long grd_out_nbr(out_nbr+1);
  
  // Check for monotonicity
  bool mnt=mnt_chk(grd_in,grd_in_nbr);
  if(!mnt) err_prn(sbr_nm,"grd_in not monotonic");
  mnt=mnt_chk(grd_out,grd_out_nbr);
  if(!mnt) err_prn(sbr_nm,"grd_out not monotonic");
  
  // Find direction of monotonicity
  bool in_ncr=mnt_ncr(grd_in,grd_in_nbr);
  bool in_dcr=!in_ncr;
  bool out_ncr=mnt_ncr(grd_out,grd_out_nbr);
  bool out_dcr=!out_ncr;
  
  // Convert input and output arrays to monotonic increasing arrays
  if(in_dcr){
    if(dbg_lvl >= dbg_sbr) dbg_prn("Reversing input grid"+sbr_nm);
    rvr_vec(grd_in,grd_in_nbr);
    rvr_vec(dat_in,in_nbr);
  } // endif in_dcr
  if(out_dcr){
    if(dbg_lvl >= dbg_sbr) dbg_prn("Reversing output grid"+sbr_nm);
    rvr_vec(grd_out,grd_out_nbr);
  } // endif out_dcr
  
  // Input and output coordinates both increase monotonically
  if(xtr_LHS.xtr_err && xtr_RHS.xtr_err){
    // Vet input more stringently
    if(grd_out[grd_out_nbr-1] < grd_in[0]){
      std::cerr << "ERROR: rbn_vec() has grd_out[grd_out_nbr-1 = " << grd_out_nbr-1 << "] = " << grd_out[grd_out_nbr-1] << " < grd_in[0] = " << grd_in[0] << std::endl;
      err_prn("Exiting");
    } // endif
    if(grd_out[0] > grd_in[grd_in_nbr-1]){
      std::cerr << "ERROR: rbn_vec() has grd_out[0] = " << grd_out[0] << " > grd_in[grd_in_nbr-1 = " << grd_in_nbr-1 << "] = " << grd_in[grd_in_nbr-1] << std::endl;
      err_prn("Exiting");
    } // endif
  } // endif xtr_err
  
  // Initialize coordinate arrays
  for(in_idx=0;in_idx<in_nbr;in_idx++){
    crd_in[in_idx]=0.5*(grd_in[in_idx]+grd_in[in_idx+1]);
  } // end loop over input coordinates
  for(out_idx=0;out_idx<out_nbr;out_idx++){
    crd_out[out_idx]=0.5*(grd_out[out_idx]+grd_out[out_idx+1]);
  } // end loop over output coordinates
  
  // Assemble LHS bracketing indices for each output gridpoint
  for(grd_out_idx=0;grd_out_idx<grd_out_nbr;grd_out_idx++){
    brk_lft_idx[grd_out_idx]=0;
    // Initialize bracketing index with previous value if possible
    if(grd_out_idx > 0){
      if(brk_lft_idx[grd_out_idx-1] >= 0 && brk_lft_idx[grd_out_idx-1] < grd_in_nbr)
	brk_lft_idx[grd_out_idx]=brk_lft_idx[grd_out_idx-1];
    } // endif
    // Order of conditions is important since second condition is illegal if brk_lft_idx >= grd_in_nbr
    while((brk_lft_idx[grd_out_idx] < grd_in_nbr) &&
	  (grd_in[brk_lft_idx[grd_out_idx]] <= grd_out[grd_out_idx])){
      brk_lft_idx[grd_out_idx]++;
    } // end while
    // Decrement bracketing index since last loop iteration overshoots target
    brk_lft_idx[grd_out_idx]--;

    /*    // If next input gridpoint equals current output gridpoint put
    if(brk_lft_idx[grd_out_idx] < grd_in_nbr-1)
      if(grd_in[brk_lft_idx[grd_out_idx]+1] == grd_out[grd_out_idx])
      brk_lft_idx[grd_out_idx]++; */

    /*    // fxm: 20000819 potential segfault here when brk_lft_idx[grd_out_idx] == grd_in_nbr
    if(grd_in[brk_lft_idx[grd_out_idx]] != grd_out[grd_out_idx]){
      // Decrement bracketing index unless while loop was broken by an equality
      brk_lft_idx[grd_out_idx]=brk_lft_idx[grd_out_idx]-1;
      } // endif */

  } // end loop over output coordinates
  
  // Assemble RHS bracketing indices for each output gridpoint
  for(grd_out_idx=grd_out_nbr-1;grd_out_idx>=0;grd_out_idx--){
    brk_rgt_idx[grd_out_idx]=grd_in_nbr-1;
    // Initialize bracketing index with previous value if possible
    if(grd_out_idx < grd_out_nbr-1){
      if(brk_rgt_idx[grd_out_idx+1] < grd_in_nbr && brk_rgt_idx[grd_out_idx+1] >= 0) 
	brk_rgt_idx[grd_out_idx]=brk_rgt_idx[grd_out_idx+1];
    } // endif
    // Order of conditions is important since second condition is illegal if brk_rgt_idx < 0
    while((brk_rgt_idx[grd_out_idx] >= 0) &&
	  (grd_in[brk_rgt_idx[grd_out_idx]] >= grd_out[grd_out_idx])){
      brk_rgt_idx[grd_out_idx]--;
    } // end while
    // Increment bracketing index since last loop iteration overshoots target
    brk_rgt_idx[grd_out_idx]++;

    /*    // fxm: 20000819 potential segfault here when brk_rgt_idx[grd_out_idx] == -1
    if(grd_in[brk_rgt_idx[grd_out_idx]] != grd_out[grd_out_idx]){
      // Increment bracketing index unless while loop was broken by an equality
      brk_rgt_idx[grd_out_idx]++;
      } // endif */
  } // end loop over output coordinates
  
  // Vet bracketing indices
  for(grd_out_idx=0;grd_out_idx<grd_out_nbr;grd_out_idx++){
    if(brk_lft_idx[grd_out_idx] > brk_rgt_idx[grd_out_idx]){
      err_prn(sbr_nm,"brk_lft_idx[grd_out_idx] > brk_rgt_idx[grd_out_idx]");
    } // endif
  } // end loop over output coordinates
  
  // Precompute all linearly interpolated values for use below
  Xtr_cls xtr_fll_lnr("xtr_fll_lnr");
  // Input values interpolated to input grid (Recall values are input on coordinates, not grid)
  rcd+=ntp_vec(in_nbr,crd_in,dat_in,grd_in_nbr,grd_in,dat_ntp_in,xtr_fll_lnr,xtr_fll_lnr);
  // Input values interpolated to output grid
  rcd+=ntp_vec(in_nbr,crd_in,dat_in,grd_out_nbr,grd_out,dat_ntp_out,xtr_fll_lnr,xtr_fll_lnr);
  // Precompute area of all trapezoids defined by input grid for use below
  trp_area_vec(in_nbr,grd_in,dat_ntp_in,trp_area);
  
  // Assign an output value to each output region
  for(out_idx=0;out_idx<out_nbr;out_idx++){
    out_idxp1=out_idx+1;
    // Traverse truth table by examining all possible values of brk_lft_idx[out_idx]
    if(brk_lft_idx[out_idx] < -1){
      err_prn(sbr_nm,"brk_lft_idx[out_idx] < -1");
    }else if(brk_lft_idx[out_idx] == -1){
      // LHS extrapolation required
      if(xtr_LHS.xtr_vrb){
	std::cerr << "WARNING: rbn_vec() output bin value centered at crd_out[" << out_idx << "] between interface coordinates grd_out[" << out_idx << "] = " << grd_out[out_idx] << " and grd_out[" << out_idxp1 << "] = " << grd_out[out_idxp1] << " requires LHS extrapolation beyond leftmost valid interface coordinates at grd_in[0] = " << grd_in[0] << " and grd_in[1] = " << grd_in[1] << ". Nearest valid datum is dat_in[0] = " << dat_in[0];
      } // endif xtr_vrb
      if(brk_lft_idx[out_idxp1] == -1){
	// Fully degenerate case: current output bin has no overlap with input bins
	if(xtr_LHS.xtr_vrb) std::cerr << "Full LHS extrpolation required. Distances of bin endpoints from valid data are " << grd_in[0]-grd_out[out_idx] << " and " << grd_in[0]-grd_out[out_idxp1] << std::endl;
	// Extrapolation options are presented in decreasing order of preference
	if(!xtr_LHS.xtr_fll){
	  err_prn(sbr_nm,"Full LHS extrapolation required but not permitted");
	}else if(xtr_LHS.xtr_fll_nil){
	  dat_out[out_idx]=0.0;
	}else if(xtr_LHS.xtr_fll_ngh){
	  dat_out[out_idx]=dat_in[0];
	}else if(xtr_LHS.xtr_fll_lnr){
	  err_prn(sbr_nm,"xtr_fll_lnr not implemented yet");
	}else{
	  err_prn(sbr_nm,"Unknown xtr_typ_LHS in fll branch");
	} // endif xtr_typ_LHS
      }else if(brk_lft_idx[out_idxp1] > -1){ // brk_lft_idx[out_idxp1] == -1
	// Half degenerate case: current output bin has partial overlap with input bins
	ovl_frc=(grd_out[out_idxp1]-grd_in[0])/(grd_out[out_idxp1]-grd_out[out_idx]);
	if(xtr_LHS.xtr_vrb) std::cerr << "Partial LHS extrapolation required. Overlap fraction is " << ovl_frc << std::endl;
	// Extrapolation options are presented in decreasing order of preference
	if(!xtr_LHS.xtr_prt){
	  err_prn(sbr_nm,"Partial LHS extrapolation required but not permitted");
	}else if(xtr_LHS.xtr_prt_frc){
	  // Compute average of overlap input bins
	  crd_dlt_rgt=grd_out[out_idxp1]-grd_in[brk_lft_idx[out_idxp1]];
	  trp_area_rgt=0.5*crd_dlt_rgt*(dat_ntp_out[out_idxp1]+dat_ntp_in[brk_lft_idx[out_idxp1]]);
	  crd_dlt_ttl=crd_dlt_rgt;
	  trp_area_ttl=trp_area_rgt;
	  trp_nbr=brk_lft_idx[out_idxp1];
	  for(trp_idx=0;trp_idx<trp_nbr;trp_idx++){
	    in_idx=trp_idx;
	    trp_area_ttl+=trp_area[in_idx];
	    crd_dlt_ttl+=grd_in[in_idx+1]-grd_in[in_idx];
	  } // end loop over inner trapezoids
	  if(crd_dlt_ttl <= 0.0) err_prn(sbr_nm,"crd_dlt_ttl <= 0.0 in xtr_LHS branch");
	  ovl_avg=trp_area_ttl/crd_dlt_ttl ;
	  if(xtr_LHS.xtr_prt_wgt) dat_out[out_idx]=ovl_frc*ovl_avg; else dat_out[out_idx]=ovl_avg;
	}else if(xtr_LHS.xtr_prt_lnr){
	  err_prn(sbr_nm,"xtr_prt_lnr not implemented yet");
	}else if(xtr_LHS.xtr_prt_ngh){
	  dat_out[out_idx]=dat_in[0];
	}else if(xtr_LHS.xtr_prt_nil){
	  dat_out[out_idx]=0.0;
	}else{
	  err_prn(sbr_nm,"Unknown xtr_typ_LHS in prt branch");
	} // endif xtr_typ_LHS
      }else{ // brk_lft_idx[out_idxp1] == -1
	err_prn(sbr_nm,"Unforeseen brk_lft_idx in xtr_typ_LHS");
      } // brk_lft_idx[out_idxp1] == -1
      if(xtr_LHS.xtr_vrb) std::cerr << xtr_LHS.xtr_sng << " yields dat_out[" << out_idx << "] = " << dat_out[out_idx] << std::endl;
    }else if((brk_lft_idx[out_idx] < grd_in_nbr) && 
	     (brk_rgt_idx[out_idxp1] < grd_in_nbr)){ // brk_lft_idx[out_idx] == -1
      // Normal case: 99% of execution time for large arrays should be in this block
      // Current and next output gridpoints are both interpolable
      // Find out how many input crds lay between current and next output gridpoints
      trp_nbr=brk_lft_idx[out_idxp1]-brk_rgt_idx[out_idx];
      if(trp_nbr < 0){
	// Case A: Current and next bracketing indices are the same
	// This means the output crd resolution is finer than the input crd resolution
	// Use trapezoidal are current point (this is the same as using just a LHS trapezoid)
	dat_out[out_idx]=0.5*(dat_ntp_out[out_idx]+dat_ntp_out[out_idxp1]);
      }else if(trp_nbr <= in_nbr){
	// Case B: Current and previous bracketing indices differ by one
	// This means the output crd resolution roughly equals the input crd resolution
	// This case can be handled by the following code since the number of inner trapezoids
	// is 0 and the loop will not be executed. 
	// The input/output area will just be the sum of the LHS and RHS trapezoids
	
	// Case C: Current and previous bracketing indices differ by many.
	// This means the output crd resolution is coarser than the input crd resolution.
	// First, figure out area under the RHS and LHS bookend trapezoids:
	crd_dlt_lft=grd_in[brk_rgt_idx[out_idx]]-grd_out[out_idx];
	trp_area_lft=0.5*crd_dlt_lft*(dat_ntp_out[out_idx]+dat_ntp_in[brk_rgt_idx[out_idx]]);
	crd_dlt_rgt=grd_out[out_idxp1]-grd_in[brk_lft_idx[out_idxp1]];
	trp_area_rgt=0.5*crd_dlt_rgt*(dat_ntp_out[out_idxp1]+dat_ntp_in[brk_lft_idx[out_idxp1]]);
	// NB: We use trapezoids to integrate. But we are unable to use the classical 
	// Trapezoidal rule because the abscissas are, in general, unevenly spaced.
	// Therefore we use the trapezoidal rule generalized to unevenly spaced abscissas.
	// If the input grid equals the output grid, there is one inner trapezoid.
	
	// Trapezoids are indexed with the C (0-based) convention where
	// Trapezoid 0 is bounded by the abscissae grd_in[0] and grd_in[1] and
	// Trapezoid N-1 is bounded by the abscissae grd_in[N-1] and grd_in[N].
	crd_dlt_ttl=crd_dlt_lft+crd_dlt_rgt;
	trp_area_ttl=trp_area_lft+trp_area_rgt;
	if(brk_lft_idx[out_idx] != brk_rgt_idx[out_idx]) trp_srt_idx=brk_lft_idx[out_idx]; else trp_srt_idx=brk_lft_idx[out_idx]-1;
	for(trp_idx=0;trp_idx<trp_nbr;trp_idx++){
	  in_idx=trp_srt_idx+trp_idx+1; // NB: trp_idx+1 is correct
	  trp_area_ttl+=trp_area[in_idx];
	  crd_dlt_ttl+=grd_in[in_idx+1]-grd_in[in_idx];
	} // end loop over inner trapezoids
	// This relationship results from requiring output rectangle to be equal in area to 
	// total of input trapezoids (between current and next gridpoints).
	if(crd_dlt_ttl <= 0.0) err_prn(sbr_nm,"crd_dlt_ttl <= 0.0 in main branch"); else dat_out[out_idx]=trp_area_ttl/crd_dlt_ttl; 
      }else if(trp_nbr > in_nbr){
	std::cerr << "brk_lft_idx[" << out_idx << "] = " << brk_lft_idx[out_idx] << ", brk_rgt_idx[" << out_idx << "] = " << brk_rgt_idx[out_idx] << std::endl;
	std::cerr << "brk_lft_idx[" << out_idxp1 << "] = " << brk_lft_idx[out_idxp1] << ", brk_rgt_idx[" << out_idxp1 << "] = " << brk_rgt_idx[out_idxp1] << std::endl;
	std::cerr << "trp_nbr = " << trp_nbr << std::endl;
	err_prn(sbr_nm,"Unforeseen trp_nbr"); 
      } // endif trp_nbr == 0
    }else if(brk_rgt_idx[out_idxp1] == grd_in_nbr){ // brk_lft_idx[out_idx] == -1
      // RHS extrapolation required
      if(xtr_RHS.xtr_vrb) std::cerr << "WARNING: rbn_vec() output bin value centered at crd_out[" << out_idx << "] between interface coordinates grd_out[" << out_idx << "] = " << grd_out[out_idx] << " and grd_out[" << out_idxp1 << "] = " << grd_out[out_idxp1] << " requires RHS extrapolation beyond rightmost valid interface coordinates at grd_in[" << grd_in_nbr-2 << "] = " << grd_in[grd_in_nbr-2] << " and grd_in[" << grd_in_nbr-1 << "] = " << grd_in[grd_in_nbr-1] << ". Nearest valid datum is dat_in[" << in_nbr-1 << "] = " << dat_in[in_nbr-1] << std::endl;
      if(brk_rgt_idx[out_idx] == grd_in_nbr){
	// Fully degenerate case: current output bin has no overlap with input bins
	if(xtr_RHS.xtr_vrb) std::cerr << "Full RHS extrpolation required. Distances of bin endpoints from valid data are " << grd_out[out_idx]-grd_in[grd_in_nbr-1] << " and " << grd_out[out_idxp1]-grd_in[grd_in_nbr-1] << std::endl;
	// Extrapolation options are presented in decreasing order of preference
	if(!xtr_RHS.xtr_fll){
	  err_prn(sbr_nm,"Full RHS extrapolation required but not permitted");
	}else if(xtr_RHS.xtr_fll_nil){
	  dat_out[out_idx]=0.0;
	}else if(xtr_RHS.xtr_fll_ngh){
	  dat_out[out_idx]=dat_in[in_nbr-1];
	}else if(xtr_RHS.xtr_fll_lnr){
	  err_prn(sbr_nm,"xtr_fll_lnr not implemented yet");
	}else{
	  err_prn(sbr_nm,"Unknown xtr_typ_RHS");
	} // endif xtr_typ_RHS
      }else{ // brk_rgt_idx[out_idx] == grd_in_nbr
	// Half degenerate case: current output bin has partial overlap with input bins
	ovl_frc=(grd_in[grd_in_nbr-1]-grd_out[out_idx])/(grd_out[out_idxp1]-grd_out[out_idx]);
	if(xtr_RHS.xtr_vrb) std::cerr << "Partial RHS extrapolation required. Overlap fraction is " << ovl_frc << std::endl;
	// Extrapolation options are presented in decreasing order of preference
	if(!xtr_RHS.xtr_prt){
	  err_prn(sbr_nm,"Partial RHS extrapolation required but not permitted");
	}else if(xtr_RHS.xtr_prt_frc){
	  // Compute average of overlap input bins
	  crd_dlt_lft=grd_in[brk_rgt_idx[out_idx]]-grd_out[out_idx];
	  trp_area_lft=0.5*crd_dlt_lft*(dat_ntp_out[out_idx]+dat_ntp_in[brk_rgt_idx[out_idx]]);
	  crd_dlt_ttl=crd_dlt_lft;
	  trp_area_ttl=trp_area_lft;
	  trp_nbr=grd_in_nbr-1-brk_rgt_idx[out_idx]; // NB: grd_in_nbr-1 is correct
	  for(trp_idx=0;trp_idx<trp_nbr;trp_idx++){
	    in_idx=brk_lft_idx[out_idx]+trp_idx+1; // NB: trp_idx+1 is correct
	    trp_area_ttl+=trp_area[in_idx];
	    crd_dlt_ttl+=grd_in[in_idx+1]-grd_in[in_idx];
	  } // end loop over inner trapezoids
	  if(crd_dlt_ttl <= 0.0) err_prn(sbr_nm,"crd_dlt_ttl <= 0.0 in xtr_RHS branch");
	  ovl_avg=trp_area_ttl/crd_dlt_ttl ;
	  if(xtr_RHS.xtr_prt_wgt) dat_out[out_idx]=ovl_frc*ovl_avg; else dat_out[out_idx]=ovl_avg;
	}else if(xtr_RHS.xtr_prt_lnr){
	  err_prn(sbr_nm,"xtr_prt_lnr not implemented yet");
	}else if(xtr_RHS.xtr_prt_ngh){
	  dat_out[out_idx]=dat_in[in_nbr-1];
	}else if(xtr_RHS.xtr_prt_nil){
	  dat_out[out_idx]=0.0;
	}else{
	  err_prn(sbr_nm,"Unknown xtr_typ_RHS in prt branch");
	} // endif xtr_typ_RHS
      } // brk_rgt_idx[out_idx] == grd_in_nbr
      if(xtr_RHS.xtr_vrb) std::cerr << xtr_RHS.xtr_sng << " yields dat_out[" << out_idx << "] = " << dat_out[out_idx] << std::endl;
    }else{
      err_prn(sbr_nm,"Incomplete truth table");
      std::cerr << "brk_lft_idx[" << out_idx << "] = " << brk_lft_idx[out_idx] << ", brk_rgt_idx[" << out_idx << "] = " << brk_rgt_idx[out_idx] << std::endl;
      err_prn(sbr_nm,"Unforeseen bracketing values"); 
    } // endelse RHS extrapolate
  } // end loop over output coordinates
  
  if(dbg_lvl == dbg_crr){
    // Print input coordinates
    err_prn(sbr_nm,"Input crd_in");
    std::cerr << "idx\t" << "crd_in\t" << "crd_min\t" << "crd_max\t" << "dat_in\t" << std::endl;
    for(long idx=0;idx<in_nbr;idx++) std::cerr << idx << "\t" << crd_in[idx] << "\t" << grd_in[idx] << "\t" << grd_in[idx+1] << "\t" << dat_in[idx] << std::endl;
  } // endif dbg
  
  if(dbg_lvl == dbg_crr){
    // Print input grid
    err_prn(sbr_nm,"Input grd_in");
    std::cerr << "idx\t" << "grd_in\t" << "dat_ntp_in\t" << std::endl;
    for(long idx=0;idx<grd_in_nbr;idx++) std::cerr << idx << "\t" << grd_in[idx] << "\t" << dat_ntp_in[idx] << std::endl;
  } // endif dbg
  
  if(dbg_lvl == dbg_crr){
    // Print bracketing indices
    err_prn(sbr_nm,"Grid output");
    std::cerr << "idx\t" << "grd_out\t" << "brk_lft\t" << "brk_rgt\t" << "dat_ntp_out\t" << std::endl;
    for(long idx=0;idx<grd_out_nbr;idx++) std::cerr << idx << "\t" << grd_out[idx] << "\t" << brk_lft_idx[idx] << "\t" << brk_rgt_idx[idx] << "\t" << dat_ntp_out[idx] << std::endl;
  } // endif dbg
  
  if(dbg_lvl == dbg_crr){
    // Print trapezoid info
    err_prn(sbr_nm,"Trapezoid output");
    std::cerr << "idx\t" << "grd_out\t" << "outp1\t" << "trp_nbr\t" << "trp_idx\t" << std::endl;
    for(out_idx=0;out_idx<out_nbr;out_idx++){
      trp_nbr=brk_lft_idx[out_idx+1]-brk_rgt_idx[out_idx];
      std::cerr << out_idx << "\t" << grd_out[out_idx] << "\t" << grd_out[out_idx+1] << "\t" << trp_nbr << "\t";
      if(brk_lft_idx[out_idx] != brk_rgt_idx[out_idx]) trp_srt_idx=brk_lft_idx[out_idx]; else trp_srt_idx=brk_lft_idx[out_idx]-1;
      for(trp_idx=0;trp_idx<trp_nbr-1;trp_nbr++){ // NB: Loop ends early
	in_idx=trp_srt_idx+trp_idx+1; // NB: trp_idx+1 is correct
	std::cerr << in_idx << ","; 
      } // end loop over trp
      if (trp_nbr > 0) std::cerr << trp_srt_idx+trp_nbr << std::endl; else std::cerr << "none" << std::endl;
    } // end loop over crd_out
  } // endif dbg
  
  if(dbg_lvl == dbg_crr){
    // Print results of rebinning
    err_prn(sbr_nm,"Coordinate output");
    std::cerr << "idx\t" << "crd_out\t" << "crd_min\t" << "crd_max\t" << "dat_out\t" << std::endl;
    for(long idx=0;idx<out_nbr;idx++) std::cerr << idx << "\t" << crd_out[idx] << "\t" << grd_out[idx] << "\t" << grd_out[idx+1] << "\t" << dat_out[idx] << std::endl;
  } // endif dbg
  
  // Convert input and output arrays to original direction of monotonicity
  if(in_dcr){
    if(dbg_lvl >= dbg_sbr) err_prn(sbr_nm,"Un-reversing input grid");
    rvr_vec(grd_in,grd_in_nbr);
    rvr_vec(dat_in,in_nbr);
  } // endif in_dcr
  if(out_dcr){
    if(dbg_lvl >= dbg_sbr) err_prn(sbr_nm,"Un-reversing output grid");
    rvr_vec(grd_out,grd_out_nbr);
    rvr_vec(dat_out,out_nbr);
  } // endif in_dcr
  
  if(dbg_lvl >= dbg_sbr) dbg_prn(sbr_nm,"Exiting...");
  return rcd;
} // end rbn_vec()

void
trp_area_vec
(const long in_nbr, // I [nbr] Input coordinate size
 const prc_cmp *grd_in, // I [crd] Input grid
 const prc_cmp *dat_in, // I [frc] Input data
 prc_cmp *trp_area) // O [frc] Area of each trapezoid defined by input data
{
  /* Purpose: Compute area of each trapezoid in input array of trapezoids
     N trapezoids are centered within N+1 gridpoints
     Trapezoids are indexed with the C (0-based) convention where
     Trapezoid 0 is bounded by the abscissae grd_in[0] and grd_in[1]
     Trapezoid 1 is bounded by the abscissae grd_in[1] and grd_in[2]
     Trapezoid k is bounded by the abscissae grd_in[k-1] and grd_in[k]
     Trapezoid N-2 is bounded by the abscissae grd_in[N-2] and grd_in[N-1]
     Trapezoid N-1 is bounded by the abscissae grd_in[N-1] and grd_in[N] */
  // Main Code
  for(long in_idx=0;in_idx<in_nbr;in_idx++){
    trp_area[in_idx]=0.5*(grd_in[in_idx+1]-grd_in[in_idx])*(dat_in[in_idx+1]+dat_in[in_idx]);
  } // end loop over trapezoids
} // end trp_area_vec()
