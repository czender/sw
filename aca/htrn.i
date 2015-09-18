// Yorick commands for examining HITRAN data
// #include "/home/zender/aca/aca.i"
// htrn,"swnb_H2O"
// htrn,"swnb_H2O_all"
// htrn,"swnb_OH"
// htrn,"swnb_OH_all"
// htrn,"swnb_O3_all"
// htrn,"swnb_N2O"
// htrn,"swnb_NO"
// htrn,"swnb_N2_all"
// htrn,"swnb_NO2"
// htrn,"swnb_SO2"
// htrn,"swnb_H2O2"

// #include "style.i"

func htrn(mlc_sng)
{
  fma;
  window,style="vgbox.gs";

  fl_in_dir="/data/zender/aca";
  if(is_void(mlc_sng)) fl_in_stb="swnb_H20.nc"; else fl_in_stb=mlc_sng;
  fl_in=fl_in_dir+"/"+fl_in_stb+".nc";

  // Get netCDF data
  nc_in=openb(fl_in);
  wvl=nc_in.wvl_ctr;
  // dat=nc_in.abs_xsx_O3_cold
  dat=nc_in.S_d_abs_cff_mss;
  // dat=nc_in.S_p_abs_cff_mss;
  close,nc_in;
  
  list=where(dat > 0);
  nbr_abs_bnd=numberof(list);
  dbg=0
    if(dbg){
      if(nbr_abs_bnd){
	for(idx=1;idx<=nbr_abs_bnd;idx++){
	  idx_org=list(idx);
	  //    write(format="wvl(%d) dat(%g) = %g\n",list,wvl,dat);
	  write(format="wvl(%d) dat(%g) = %g\n",idx_org,wvl(idx_org),dat(idx_org));
	  //    print,"wvl(",idx_org,") = ",wvl(idx_org)*1.e6," microns, dat = ",dat(idx_org),"\n";
	} // end loop over idx
      }else{
	print,"\n";
      } // endif nbr_abs_bnd
    } // endif dbg
  
  limits,0,5,,"e";
  logxy,0,1;
  plg,dat,wvl*1.e6,marks=0;
  // plg,nsl*1.e-9,wvl*1.e6,legend="Surface Insolation",marks=0
} // end htrn()




