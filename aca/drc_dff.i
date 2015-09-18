// Yorick commands for computing direct-diffuse flux diagnostics
// #include "/home/zender/aca/drc_dff.i"
// drc_dff,"951011_1200_arese_mdl_clr_aer"
// drc_dff,"951015_1200_arese_mdl_clr_aer"
// drc_dff,"951030_1200_arese_mdl_clr_aer"

// #include "style.i"

func drc_dff(fl)
{
  fma;
  window,style="vgbox.gs";

  fl_in_dir="/data/zender/arese/mdl";
  if(is_void(fl)) fl_in_stb="951015_1200_arese_mdl_clr_aer"; else fl_in_stb=fl;
  fl_in=fl_in_dir+"/"+fl_in_stb+".nc";

  // Get netCDF data
  nc_in=openb(fl_in);
  wvl=nc_in.wvl_ctr;
  wvl_sz=nc_in.wvl_sz;
  flx_spc_dwn_sfc=nc_in.flx_spc_dwn_sfc;
  flx_bb_dwn_sfc=nc_in.flx_bb_dwn_sfc;
  close,nc_in;
  
  // Get spectral response function data
  fl_in="/data/zender/aca/nst_TSBR.nc";
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
  flx_nst_dwn_sfc=sum(flx_spc_dwn_sfc(nir)*wvl_sz(nir))
  write(format="filtered flx_vsb_dwn_sfc = %g W/m2, flx_nir_dwn_sfc = %g W/m2\n",flx_vsb_dwn_sfc,flx_nir_dwn_sfc);
  write(format="unfiltered flx_bb_dwn_sfc = %g W/m2\n",flx_bb_dwn_sfc);
  write(format="filtered flx_nst_dwn_sfc = %g W/m2\n",flx_vsb_dwn_sfc+flx_nir_dwn_sfc);
  
} // end drc_dff()




