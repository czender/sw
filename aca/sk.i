// Yorick commands for computing Stefan Kinne's diagnostics
// #include "/home/zender/aca/sk.i"
// sk,"sk_01"

// #include "style.i"

func sk(fl)
{
  fma;
  window,style="vgbox.gs";

  fl_in_dir="/data/zender/tmp";
  if(is_void(fl)) fl_in_stb="sk_01"; else fl_in_stb=fl;
  fl_in=fl_in_dir+"/"+fl_in_stb+".nc";

  // Get netCDF data
  nc_in=openb(fl_in);
  wvl=nc_in.wvl_ctr;
  wvl_sz=nc_in.wvl_sz;
  flx_spc_dwn_sfc=nc_in.flx_spc_dwn_sfc;
  flx_bb_dwn_sfc=nc_in.flx_bb_dwn_sfc;
  close,nc_in;
  
  vsb=where(wvl < .7e-6);
  nir=where(wvl > .7e-6);
  flx_vsb_dwn_sfc=sum(flx_spc_dwn_sfc(vsb)*wvl_sz(vsb))
  flx_nir_dwn_sfc=sum(flx_spc_dwn_sfc(nir)*wvl_sz(nir))
  write(format="flx_vsb_dwn_sfc = %g W/m2, flx_nir_dwn_sfc = %g W/m2\n",flx_vsb_dwn_sfc,flx_nir_dwn_sfc);
  write(format="flx_bb_dwn_sfc = %g, %g W/m2\n",flx_bb_dwn_sfc,flx_vsb_dwn_sfc+flx_nir_dwn_sfc);
  
} // end sk()




