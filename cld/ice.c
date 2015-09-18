/* $Author: zender $
 * $Date$
 * $Id$
 * $Revision$
 * $Locker:  $
 * $RCSfile: ice.c,v $
 * $Source: /home/zender/cvs/cld/ice.c,v $
 * $Id$
 * $State: Exp $
 * */

/* Purpose: Generation and manipulation of radiative parameters and 
   heating rates. */ 

/* $Log: not supported by cvs2svn $
/* Revision 1.1.1.1  1998-09-15 02:06:42  zender
/* Imported sources
/*
 * Revision 5.4  1995/06/19  04:54:58  zender
 * have made many changes. got rid of fp_in,out,err. added flat updraft
 * switch. changed -S -s to both refer to ice. got rid of while loop that
 * became infinite when there was no cloud. changed vertical netCDF coord.
 * to altitude. fixed bug in appending date and rev. to cmdline.
 *
 * Revision 5.3  1995/05/24  23:11:17  zender
 * updated and validated the cloud model for ansi compatability,
 * removed berkeley calls, streamlined architecture dependent
 * fortran --> c calling methodology. validated results against
 * comps II runs.
 *
 * Revision 5.2  1995/05/23  00:43:16  zender
 * synchronization checkin prior to temperature sensitivity study
 * for SPCP stable cloud mods.
 *
 * Revision 5.1  1994/02/14  21:28:54  zender
 * Added DMC94 homogeneous ammonium sulphate solution nucleation,
 * added -D 79 SWSCF/SWCF diagnostic. About to change default shape
 * to hexagonal ice crystals, but will save bullet routines.
 * NetCDF routines are not keeping up with the changes....
 *
 * Revision 5.0  1993/08/30  01:55:35  zender
 * this version was used to make the paper (preprint at least).
 * about to start working on an animation.
 *
 * Revision 4.12  1993/07/21  23:26:28  zender
 * added Liou_IR_fudge() as a default routine whenever LIOU is true.
 *
 * Revision 4.11  1993/06/19  22:05:25  zender
 * added ebert&curry parameterization to LW computations, and
 * made this the default for the region outside of the cloud grid.
 *
 * Revision 4.10  1993/06/11  04:35:31  zender
 * Played w/ Makefiles, made all CRAY2 refs. into plain CRAY,
 * and put all shavano stuff on peace. changed surface type to
 * land for predictable albedos, changed tropics SAZ to 30 deg.,
 * removed gaussians from initial distributions. changed plev = 113.
 * slingo option doesn't work for some reason (-D 53)
 *
 * Revision 4.9  1993/06/09  00:59:14  zender
 * prepared the -Z initial_distribution_type switch to handle everyone's
 * different cloud. Made Knollenberg fit his observations, prepared
 * optical properties for 50 sizes and removed all num_size == 40 specific
 * stuff. changed tau_CCN to 1. s.
 *
 * Revision 4.8  1993/05/29  00:01:13  zender
 * fixed the Liou_fudge() bug, fixed a wind_speed bug, added
 * Dowling and Radke's crystal distribution function as default.
 *
 * Revision 4.7  1993/05/27  16:08:24  zender
 * A synchronization check-in for all cloud programs.
 *
 * Revision 4.6  1993/05/27  14:41:42  zender
 * intermediate version to allow for a bug search. -X option now either calls
 * Liou_fudge() (default) or Liou_interp(), but neither works well anymore.
 * even the Mie option crashes in Tropics mode!
 *
 * Revision 4.5  1993/04/30  01:54:34  zender
 * added -M (enviromental profile type) and polished up -X (radiative
 * crystal type) switches. also ANSI'd some more routines. next step
 * is to get rid of bullet habit completely, and grow things as columns.
 *
 * Revision 4.4  1993/04/22  23:47:27  zender
 * big change was integrated Liou ice crystal data with -X,
 * and rearranging the spectral intervals in CCM2 colmod code
 * so that CCM2 bands 14,15 are now my bands 15,16, and CCM2 16
 * went to 14.  hex crosssections are still too big compared to
 * corresponding bullet/spheres, tinkering with interp. geometry.
 *
 * Revision 4.3  1993/04/16  01:12:21  zender
 * incorporated some Liou data, changed integration of optical
 * parameters et al. to depend on optical cross sections, not
 * efficiencies.  LIOU switch doesn't work, bug in bilinear_interp().
 *
 * Revision 4.2  1993/03/08  00:57:47  zender
 * this version contains working working data i/o in 3 formats:
 * ASCII, CZ_BINARY (my vanilla binary format), and NETCDF.
 * since netCDF is transparent across platforms (i.e., cray and sun),
 * the next version will only allow netCDF i/o.
 *
 * Revision 4.1  1993/02/27  00:05:18  zender
 * incorporated ascii writing automatically on the CRAY in main,
 * merged ioinput and iooutput into iographs. added albedo vs.
 * emissivity graph. alphabetized output. might work on netCDF now.
 *
 */

#include <stdio.h>              /* stderr, FILE, NULL, etc. */
#include <math.h>               /* M_PI, and the intrinsic functions */ 

/* my header files */
#include "defs.h"               /* YES, NO, FILESIZE, and more */
#include "globals.h"            /* externs announced in main() */

#define SIZE_MIE_DISKIO 40

void Bohren_Mie(wavelength,num_band,
		ice_real_idx,ice_imag_idx,
		equiv_radius,num_size,
		size_param,
		absorption_eff,
		scattering_eff,
		extinction_eff,
		absorption_x_sec,
		scattering_x_sec,
		extinction_x_sec,
		equiv_rad_squared,
		asymmetry_param)
     float *wavelength,*ice_real_idx,*ice_imag_idx;
     float *equiv_radius;
     float **absorption_x_sec;
     float **scattering_x_sec;
     float **extinction_x_sec;
     float *equiv_rad_squared;
     float **absorption_eff,**extinction_eff,**scattering_eff,**size_param;
     float **asymmetry_param;
     int num_size,num_band;
{
  void FORTRAN_callbh();
  
  Boolean NO_MIE_FILE_FLAG;

  float equiv_radius_microns;
  float float_foo;
  float pi_r_squared;
  float wavelength_microns;
  
  int size;
  int num_angles;
  int band;
  
  FILE *fp_mie;
  char mie_file[FILESIZE];

  NO_MIE_FILE_FLAG=False;

  if(num_band > 19){
    (void)sprintf(mie_file,"mie_%2i_LW.dat",num_size);
  }else{
    (void)sprintf(mie_file,"mie_%2i_SW.dat",num_size);
  } /* end else */

  /* See if the parameters are on disk */ 
  if((fp_mie = fopen( mie_file, "r")) == NULL){
    (void)fprintf(stdout,"\nNo existing %s. Creating one . . .\n",mie_file);
      NO_MIE_FILE_FLAG = True;
  } /* end if */

  /* If they are then use them and leave the routine */
  if(!NO_MIE_FILE_FLAG){
    if(
       (fread((char *)size_param[0],
	      (num_band+2)*(num_size+2)*sizeof(float),1,fp_mie) != 1) ||
       (fread((char *)absorption_eff[0],
	      (num_band+2)*(num_size+2)*sizeof(float),1,fp_mie) != 1) ||
       (fread((char *)scattering_eff[0],
	      (num_band+2)*(num_size+2)*sizeof(float),1,fp_mie) != 1) ||
       (fread((char *)extinction_eff[0],
	      (num_band+2)*(num_size+2)*sizeof(float),1,fp_mie) != 1) ||
       (fread((char *)asymmetry_param[0],
	      (num_band+2)*(num_size+2)*sizeof(float),1,fp_mie) != 1) ||
       False){
      (void)fprintf(stdout,"Unable to read Mie file in Bohren_Mie()\n");
      exit(1);
    } /* endif read */ 
    (void)fclose(fp_mie);
    (void)fprintf(stdout,"Read Mie data from %s\n",mie_file);
  } /* endif */ 

  /* If they aren't on disk then compute them */
  if(NO_MIE_FILE_FLAG){
    /* number of angles between 0 and 90 degrees at which scattering matrix
       is calculated. 2 is the minimum non-crashable number */ 
    num_angles=2;
    
    for(band=1;band<=num_band;band++){
      wavelength_microns=wavelength[band]*MICRONS_PER_METER;
      for(size=1;size<=num_size;size++){
	equiv_radius_microns=equiv_radius[size]*MICRONS_PER_METER;
	
	/* The wavelength and radius must be in the same units so why not
	   scale them to microns in case it helps with rounding errors? */ 
	FORTRAN_callbh
	    (&num_angles,
	     &wavelength_microns,
	     &equiv_radius_microns,
	     ice_real_idx+band,
	     ice_imag_idx+band,
	     size_param[band]+size,
	     scattering_eff[band]+size,
	     extinction_eff[band]+size,
	     asymmetry_param[band]+size,
	     &debug);
	extinction_eff[band][size]=
	  (extinction_eff[band][size] > scattering_eff[band][size]) ? 
	    extinction_eff[band][size] : scattering_eff[band][size];
      absorption_eff[band][size]=
	extinction_eff[band][size]-scattering_eff[band][size];
      } /* end loop over sizes */
    } /* end loop over bands */
  } /* end if */
  
  /* find the optical cross-sections */ 
  for(size=1;size<=num_size;size++){
    pi_r_squared=M_PI*equiv_rad_squared[size];
    for(band=1;band<=num_band;band++){
      absorption_x_sec[band][size]=absorption_eff[band][size]*pi_r_squared;
      scattering_x_sec[band][size]=scattering_eff[band][size]*pi_r_squared;
      extinction_x_sec[band][size]=extinction_eff[band][size]*pi_r_squared;
    } /* end loop over bands */
  } /* end loop over sizes */

  /* If there was no disk file saving them then
     preserve the computationally expensive Mie parameters for next time */ 
  if(NO_MIE_FILE_FLAG){
    if((fp_mie = fopen( mie_file, "w")) == NULL){
      (void)fprintf(stdout,"\nError in opening mie file %s\n",mie_file);
      exit(1);
    } /* end if */
    if(
       (fwrite((char *)size_param[0],
	       (num_band+2)*(num_size+2)*sizeof(float),1,fp_mie) != 1) ||
       (fwrite((char *)absorption_eff[0],
	       (num_band+2)*(num_size+2)*sizeof(float),1,fp_mie) != 1) ||
       (fwrite((char *)scattering_eff[0],
	       (num_band+2)*(num_size+2)*sizeof(float),1,fp_mie) != 1) ||
       (fwrite((char *)extinction_eff[0],
	       (num_band+2)*(num_size+2)*sizeof(float),1,fp_mie) != 1) ||
       (fwrite((char *)asymmetry_param[0],
	       (num_band+2)*(num_size+2)*sizeof(float),1,fp_mie) != 1) ||
       False){
      (void)fprintf(stdout,"Unable to write Mie file in main()\n");
      exit(1);
    } /* endif fwrite */ 
    (void)fclose(fp_mie);
    (void)fprintf(stdout,"Wrote Mie data to %s\n",mie_file);
  } /* end if */
  
  if(debug == 40){
    for(band=1;band<=num_band;band++){
      (void)fprintf(stdout,"\nwband\tlambda\tsbin\tequiv_r\tnreal\tnimag\tX\tQext\tQsca\tQabs\tg\tomega\n");
      for(size=1;size<=num_size;size++){
	(void)fprintf(stdout,"%i\t%.3f\t%i\t%.2f\t%.3f\t%7.1e\t%.1f\t%.4f\t%.4f\t%.4f\t%.4f\t%.6f\n",
		      band,
		      wavelength[band]*MICRONS_PER_METER,
		      size,
		      equiv_radius[size]*MICRONS_PER_METER,
		      ice_real_idx[band],
		      ice_imag_idx[band],
		      size_param[band][size],
		      extinction_eff[band][size],
		      scattering_eff[band][size],
		      absorption_eff[band][size],
		      asymmetry_param[band][size],
		      scattering_eff[band][size]/extinction_eff[band][size]
		      );
      } /* end loop over sizes */
    } /* end loop over bands */
  } /* end debug */

} /* end Bohren_Mie */ 

void apportion_crystal_heating
  (crystal_heat_SW,
   crystal_heat_LW,
   crystal_heat_net,
   concentration,
   total_layer_conc,
   SW_spec_flux_div,
   LW_flux_div,
   absorption_x_sec,
   planck_avg_abs_x_sec,
   num_band,
   num_layer,
   num_size,
   cloud_base_idx,
   cloud_top_idx,
   num_cloudy_layers)
float **absorption_x_sec,**concentration;
float **crystal_heat_LW,**crystal_heat_SW,**crystal_heat_net;
float **planck_avg_abs_x_sec,**SW_spec_flux_div;
float *LW_flux_div,*total_layer_conc;
int num_band,num_layer,num_size;
int cloud_base_idx,cloud_top_idx,num_cloudy_layers;
{
  float dot_product_VEC();

  float layer_SW_xsection;
  float layer_LW_xsection;

  int band;
  int layer;
  int size;

  /* must make sure to zero all variables which are changed by increments */
  (void)memset((char *)crystal_heat_SW[0],'\0',(num_layer+2)*(num_size+2)*sizeof(float));
  (void)memset((char *)crystal_heat_LW[0],'\0',(num_layer+2)*(num_size+2)*sizeof(float));

/*  for(layer=1;layer<=num_layer;layer++){*/
  for(layer=cloud_base_idx;layer<=cloud_top_idx;layer++){
    if(total_layer_conc[layer] > 0.){
      for(band=1;band<=num_band;band++){
	/* sometimes the absorption efficiency can be zero for all sizes
	   in a band, or only zero for those with non-zero concentrations,
	   i.e. band 4,5 */ 
 	if((layer_SW_xsection=
	    dot_product_VEC(absorption_x_sec[band]+1,concentration[layer]+1,
			    num_size)) > 0.){
/*	  layer_SW_xsection=0.;*/
	  /* First find the "denominator" to normalize the spectral heating */ 
/*	  for(size=1;size<=num_size;size++){*/
/*	    layer_SW_xsection+=*/
/*	      concentration[layer][size]*absorption_x_sec[band][size];*/
/*	  } *//* end loop over sizes */
	  
	  /* Ready to accumulate the SW crystal heating integral */ 
	  for(size=1;size<=num_size;size++){
	    if(concentration[layer][size] > 0.){
	      /* Recall that crystal heating
		 corresponds to a negative convergence (i.e. a divergence) of the 
		 net radiative flux, so multiply by -1 to get signs correct. */ 
	      crystal_heat_SW[layer][size]+=
		-SW_spec_flux_div[layer][band]*absorption_x_sec[band][size]/
		  layer_SW_xsection;
	    } /* end if conc of this size == 0 */
	  } /* end loop over sizes */
	} /* end if dot_product absorption_x_sec * conc > 0. */
      } /* end loop over bands */
      /* Do the analogous computation for the single LW band */
      layer_LW_xsection=0.;
      for(size=1;size<=num_size;size++){
	layer_LW_xsection+=
	  concentration[layer][size]*planck_avg_abs_x_sec[layer][size];
      } /* end loop over sizes */
      
      for(size=1;size<=num_size;size++){
	if(concentration[layer][size] > 0.){
	  crystal_heat_LW[layer][size]+=
	    -LW_flux_div[layer]*planck_avg_abs_x_sec[layer][size]/
	      layer_LW_xsection;
	} /* end if conc of this size == 0 */
      } /* end loop over sizes */
    } /* end if there are crystals in layer */ 
  } /* end loop over layers */
  
  for(layer=1;layer<=num_layer;layer++){
    for(size=1;size<=num_size;size++){
      crystal_heat_net[layer][size]=
	crystal_heat_SW[layer][size]+crystal_heat_LW[layer][size];
    } /* end loop over sizes */
  } /* end loop over layers */

} /* end apportion crystal heating */ 

void integrate_LW_mass_abs_coeff
  (LW_absorption_x_sec,
   LW_spec_band_weight,
   LW_spec_vol_abs_coeff,
   LW_mass_abs_coeff,
   concentration,
   IWC,
   num_LW_band,
   num_layer,
   num_size)
float **LW_absorption_x_sec,**LW_spec_band_weight,**LW_spec_vol_abs_coeff;
float **concentration;
float *LW_mass_abs_coeff;
float *IWC;
int num_LW_band,num_layer,num_size;
{
  int band;
  int layer;
  int size;

  float float_foo;
  float float_foo2;
  float float_foo3;
  
  /* must make sure to zero all variables which are changed by increments */
  (void)memset((char *)LW_spec_vol_abs_coeff[0],'\0',(num_layer+2)*(num_LW_band+2)*sizeof(float));
  (void)memset((char *)LW_mass_abs_coeff,'\0',(num_layer+2)*sizeof(float));

  /* Weigh the grid specific optical parameters by the current 
     distribution of particles */ 
  for(layer=1;layer<=num_layer;layer++){
    for(band=1;band<=num_LW_band;band++){
      for(size=1;size<=num_size;size++){
	LW_spec_vol_abs_coeff[layer][band]+=
	  concentration[layer][size]*LW_absorption_x_sec[band][size];
      } /* end loop over sizes */
      /* find the Planckian emission band-weighted layer average mass
       absorption coefficient */ 
      LW_mass_abs_coeff[layer]+=LW_spec_band_weight[layer][band]*
	LW_spec_vol_abs_coeff[layer][band];
    } /* end loop over bands */
    /* Turn the volume absorption coefficient sigma [m^-1] into 
       the mass absorption coefficient by dividing by the 
       ice water content.  Now kappa, the LW_mass_abs_coeff[layer], 
       will be in m^2/kg */
    if(IWC[layer] > EMINUSTWENTY){
      LW_mass_abs_coeff[layer]/=IWC[layer];
    }else{
      LW_mass_abs_coeff[layer]=0.;
    } /* end else */ 
  } /* end loop over layers */

} /* end integrate_LW_mass_abs_coeff */

void init_IR_params
  (float **LW_planck_radiance,
   float *IR_bandwidth,
   float *IR_wavelength,
   float **LW_spec_band_weight,
   float **LW_absorption_x_sec,
   float **planck_avg_abs_x_sec,
   float *env_temperature,
   int num_LW_band,
   int num_layer,
   int num_size)
{
  float total_VEC();

  float Planck_function();
  float dot_product_VEC();

  float total_intensity;
  float sigma_T_fourth;
  
  int band;
  int layer;
  int size;

  for(layer=1;layer<=num_layer;layer++){
    for(band=1;band<=num_LW_band;band++){
      LW_planck_radiance[layer][band]=
	Planck_function(env_temperature[layer],IR_wavelength[band]);
    } /* end loop over bands */
  } /* end loop over layers */

  for(layer=1;layer<=num_layer;layer++){
/*    sigma_T_fourth=stefan_boltzmann_constant*pow(env_temperature[layer],4.);*/
    total_intensity=dot_product_VEC
      (LW_planck_radiance[layer]+1,IR_bandwidth+1,num_LW_band);
    for(band=1;band<=num_LW_band;band++){
      LW_spec_band_weight[layer][band]=LW_planck_radiance[layer][band]*
	IR_bandwidth[band]/total_intensity;
/*      IR_bandwidth[band]/(sigma_T_fourth/M_PI);*/
    } /* end loop over bands */
  } /* end loop over layers */

  if(debug == 52){
    (void)fprintf(stdout,"zbin\twbin\tlambda\tdlambda\tTenv\tweight\tB(T,w)\n");
    (void)fprintf(stdout,"\t\tmicrons\tmicrons\tK\t\tW/m2/u/ster\n");
    for(layer=1;layer<=num_layer;layer++){
      for(band=1;band<=num_LW_band;band++){
	(void)fprintf(stdout,"%i\t%i\t%.2f\t%.2f\t%.2f\t%.4f\t%.4f\n",
		      layer,band,
		      IR_wavelength[band]*MICRONS_PER_METER,
		      IR_bandwidth[band]*MICRONS_PER_METER,
		      env_temperature[layer],
		      LW_spec_band_weight[layer][band],
		      LW_planck_radiance[layer][band]*METERS_PER_MICRON);
      } /* end loop over bands */
      (void)fprintf(stdout,"Total of %i LW band weights = %g\n\n",num_LW_band,
		    total_VEC(LW_spec_band_weight[layer]+1,num_LW_band));
    } /* end loop over layers */;
  } /* end debug */

  /* must make sure to zero all variables which are changed by increments */
  (void)memset((char *)planck_avg_abs_x_sec[0],'\0',(num_layer+2)*(num_size+2)*sizeof(float));

  if(debug == 52){
    (void)fprintf(stdout,"zbin\tTenv\tsum(B)\tsigT4/pi\n");
    (void)fprintf(stdout,"\tK\tW/m^2\tW/m^2\n");
  } /* end debug */

  for(layer=1;layer<=num_layer;layer++){
    sigma_T_fourth=stefan_boltzmann_constant*pow(env_temperature[layer],4.);
    for(size=1;size<=num_size;size++){
      total_intensity=0.;
      for(band=1;band<=num_LW_band;band++){
	planck_avg_abs_x_sec[layer][size]+=
	  LW_absorption_x_sec[band][size]*LW_planck_radiance[layer][band]*
	    IR_bandwidth[band];
	total_intensity+=LW_planck_radiance[layer][band]*IR_bandwidth[band];
      } /* end loop over bands */
      /* this should be a good approximation with 4u < lambda < 84u */ 
/*      planck_avg_abs_x_sec[layer][size]/=sigma_T_fourth;*/
      planck_avg_abs_x_sec[layer][size]/=total_intensity;
    } /* end loop over sizes */
    if(debug == 52){
      (void)fprintf(stdout,"%i\t%.2f\t%g\t%g\n",
		    layer,
		    env_temperature[layer],
		    total_intensity,
		    sigma_T_fourth/M_PI);
    } /* end debug */
  } /* end loop over layers */
} /* end init_IR_params */

void integrate_SW_optical_properties
  (
   float **asymmetry_param,
   float **concentration,
   float dz,
   float **eff_asymmetry_param,
   float **eff_ss_albedo,
   float **eff_tau_extinction,
   float **eff_tau_scatter,
   float **extinction_x_sec,
   int num_band,
   int num_layer,
   int num_size,
   float **scattering_x_sec
   )
{
  int band;
  int layer;
  int size;
  
  static float Mie_Liou_g_correction_factor[NUM_CCM2_SPECTRAL_BANDS+2] = {
    /* the number by which the asymmetry factor g for a cirrus distribution
       of spherical particles must be multiplied to agree with Liou's g for 
       the Cs cirrus size distribution of hexagonal ice crystals */ 
    0.,
    .8869, /* CCM2 bands 1-8, corresponding to Liou's .2--.7 micron band */ 
    .8869,
    .8869,
    .8869,
    .8869,
    .8869,
    .8869,
    .8869,
    .8994, /* CCM2 bands 9-11, corresponding to Liou's .7--1.3 micron band */ 
    .8994,
    .8994,
    .9182, /* CCM2 bands 12-14, corresponding to Liou's 1.3--2.5 micron band */
    .9182, 
    .9182, 
    .9878, /* CCM2 bands 15-18, corresponding to Liou's 2.5--4.0 micron band */
    .9878,
    .9878,
    .9878,
    0.};

  static float Mie_Liou_omega_correction_factor[NUM_CCM2_SPECTRAL_BANDS+2] = {
    /* the number by which the asymmetry factor g for a cirrus distribution
       of spherical particles must be multiplied to agree with Liou's g for 
       the Cs cirrus size distribution of hexagonal ice crystals */ 
    0.,
    1.000002, /* CCM2 bands 1-8, corresponding to Liou's .2--.7 micron band */ 
    1.000002,
    1.000002,
    1.000002,
    1.000002,
    1.000002,
    1.000002,
    1.000002,
    1.000159, /* CCM2 bands 9-11, corresponding to Liou's .7--1.3 micron band */ 
    1.000159,
    1.000159,
    1.00926, /* CCM2 bands 12-14, corresponding to Liou's 1.3--2.5 micron band */
    1.00926, 
    1.00926, 
    1.036723, /* CCM2 bands 15-18, corresponding to Liou's 2.5--4.0 micron band */
    1.036723,
    1.036723,
    1.036723,
    0.};

  /* Make sure that all quantities which are "incremented" are zeroed 1st */ 
  (void)memset((char *)eff_asymmetry_param[0],'\0',(num_layer+2)*(num_band+2)*sizeof(float));
  (void)memset((char *)eff_tau_extinction[0],'\0',(num_layer+2)*(num_band+2)*sizeof(float));
  (void)memset((char *)eff_tau_scatter[0],'\0',(num_layer+2)*(num_band+2)*sizeof(float));
      
  /* Weigh the grid specific optical parameters by the current 
     distribution of particles */ 
  for(layer=1;layer<=num_layer;layer++){
    for(band=1;band<=num_band;band++){
      for(size=1;size<=num_size;size++){
	eff_asymmetry_param[layer][band]+=concentration[layer][size]*
	  asymmetry_param[band][size]*scattering_x_sec[band][size];
	eff_tau_extinction[layer][band]+=
	  concentration[layer][size]*extinction_x_sec[band][size];
	eff_tau_scatter[layer][band]+=
	  concentration[layer][size]*scattering_x_sec[band][size];
      } /* end loop over sizes */
      
      /* Apply the trapezoidal rule correction to the end values */ 
      eff_asymmetry_param[layer][band]-=.5*
	(concentration[layer][1]*
	 asymmetry_param[band][1]*scattering_x_sec[band][1]+
	 concentration[layer][num_size]*
	 asymmetry_param[band][num_size]*scattering_x_sec[band][num_size]);
      
      eff_tau_extinction[layer][band]-=.5*
	(concentration[layer][1]*extinction_x_sec[band][1]+
	 concentration[layer][num_size]*extinction_x_sec[band][num_size]);
      
      eff_tau_scatter[layer][band]-=
	.5*(concentration[layer][1]*scattering_x_sec[band][1]+
	    concentration[layer][num_size]*scattering_x_sec[band][num_size]);
      
      /* NB:  it might be wise to insert a check verifying that the
	 eff_ss_albedo is less than or equal to one right here */ 
      eff_tau_extinction[layer][band]=
	(eff_tau_extinction[layer][band] > eff_tau_scatter[layer][band]) ?
	  eff_tau_extinction[layer][band] : eff_tau_scatter[layer][band];

      if(eff_tau_extinction[layer][band] > EMINUSFIFTEEN){
	eff_ss_albedo[layer][band]=
	  eff_tau_scatter[layer][band]/eff_tau_extinction[layer][band];
	if(debug == 74){
	  /* Apply the K.-N. Liou single scatter parameter correction factor, i.e.,
	     turn the Mie spherical omega into a hexagonal ice crystal omega */ 
	  eff_ss_albedo[layer][band]*=Mie_Liou_omega_correction_factor[band];
	} /* end debug */
	if(eff_ss_albedo[layer][band] > 1.){
	  eff_ss_albedo[layer][band]=1.;
	} /* end if */
      }else{
	eff_ss_albedo[layer][band]=BOGUS_FLOAT;
	eff_tau_extinction[layer][band]=0.;
      } /* end if */

      if(eff_tau_scatter[layer][band] > EMINUSFIFTEEN){
	eff_asymmetry_param[layer][band]/=eff_tau_scatter[layer][band];
	if(debug == 74){
	  /* Apply the K.-N. Liou asymmetry parameter correction factor, i.e.,
	     turn the Mie spherical g into a hexagonal ice crystal g */ 
	  eff_asymmetry_param[layer][band]*=Mie_Liou_g_correction_factor[band];
	} /* end debug */
      }else{
	eff_asymmetry_param[layer][band]=BOGUS_FLOAT;
	eff_tau_scatter[layer][band]=0.;
      } /* end if */

      /* Now multiply the extinctions per meter by the path lengths */
      eff_tau_extinction[layer][band]*=dz;
      eff_tau_scatter[layer][band]*=dz;
      
    } /* end loop over bands */
  } /* end loop over layers */
} /* end integrate_Liou_optical_properties() routine */ 

void integrate_LW_abs_opt_depth
  (LW_abs_opt_depth,
   concentration,
   dz,
   planck_avg_abs_x_sec,
   num_layer,
   num_size)
float **concentration,**planck_avg_abs_x_sec;
float *LW_abs_opt_depth;
float dz;
int num_layer,num_size;
{
  int layer;
  int size;
  
  /* Make sure that all quantities which are "incremented" are zeroed 1st */ 
  (void)memset((char *)LW_abs_opt_depth,'\0',(num_layer+2)*sizeof(float));
      
  /* Weigh the grid specific optical parameters by the current 
     distribution of particles */ 
  for(layer=1;layer<=num_layer;layer++){
    for(size=1;size<=num_size;size++){
      LW_abs_opt_depth[layer]+=
	concentration[layer][size]*planck_avg_abs_x_sec[layer][size];
    } /* end loop over sizes */
    
    /* Apply the trapezoidal rule correction to the end values */ 
    LW_abs_opt_depth[layer]-=
      .5*(concentration[layer][1]*planck_avg_abs_x_sec[layer][1]+
	  concentration[layer][num_size]*planck_avg_abs_x_sec[layer][num_size]);
    
    /* Now multiply the extinctions per meter by the path lengths */
    LW_abs_opt_depth[layer]*=dz;
  } /* end loop over layers */
} /* end integrate_LW_abs_opt_depth() */ 

void Warren_SW_ice_init(wavelength,ice_real_idx,ice_imag_idx)
     float *wavelength,*ice_real_idx,*ice_imag_idx;
{
  /* initializes SW indices of refraction required by the mie scattering
     program */ 
  
  void FORTRAN_refice();

  float foo_temperature = 230.;
  float wavelength_microns;

  int band;
  
  static float ccm2_spectral_band[NUM_CCM2_SPECTRAL_BANDS+2][3] = {
    /* Note wavelength data is in microns, format is
       [beginning wavelength, middle wavelength, end wavelength] */ 
    0.,0.,0.,
    .200,.2225,.245, 
    .245,.255,.265,
    .265,.270,.275,
    .275,.280,.285,
    .285,.290,.295, /* band 5 */ 
    .295,.300,.305,
    .305,.3275,.350,
    .350,.55,.700,
    .700,1.0,5.000,
    .701,1.3,5.000, /* band 10 */ 
    .701,1.6,5.000,
    .702,2.0,5.000,
    .702,2.5,5.000,
    2.630,2.745,2.860,
    .703,3.0,5.000, /* band 15 */ 
    .703,3.5,5.000,
    4.160,4.25,4.550,
    4.160,4.4,4.550,
    /* These are the spectral intervals used until 4/22/93, when 
       implementation of a wavelength weighted interpolation routine,
       Liou_interp(), made it necessary that centers of wavebands be
       montonically increasing, or at least non-decreasing.  Here
       band 16 was causing problems.
       .700,1.0,5.000,
       .700,1.6,5.000,
       .700,1.6,5.000,
       .700,2.2,5.000,
       .700,2.2,5.000,
       .700,3.25,5.000,
       .700,3.25,5.000,
       2.630,2.745,2.860,
       4.160,4.355,4.550,
       4.160,4.355,4.550, */ 
    0.,0.,0.,
    };
  
  for(band=1;band<=NUM_CCM2_SPECTRAL_BANDS;band++){
    wavelength_microns=ccm2_spectral_band[band][1];
    wavelength[band]=METERS_PER_MICRON*wavelength_microns;
    /* get the indices of refraction */ 
    FORTRAN_refice(&wavelength_microns,&foo_temperature,
	    ice_real_idx+band,ice_imag_idx+band);
  } /* end loop over bands */
} /* end Warren_SW_ice_init */ 

void Warren_LW_ice_init(IR_bandwidth,IR_wavelength,
			IR_ice_real_idx,IR_ice_imag_idx)
     float *IR_bandwidth;
     float *IR_wavelength;
     float *IR_ice_imag_idx;
     float *IR_ice_real_idx;
{
  /* initializes LW indices of refraction required by the mie scattering
     program */ 
  
  void FORTRAN_refice();

  float foo_temperature = 230.;
  float IR_wavelength_microns;

  int band;
  
  static float IR_spectral_band[NUM_IR_SPECTRAL_BANDS+2][3]={
    /* Note wavelength is in microns, format is
       [beginning wavelength, middle wavelength, end wavelength] */ 
    0.,0.,0., /* band [0] */ 
    4.,4.5,5.,
    5.,5.5,6.,
    6.,6.5,7.,
    7.,7.5,8.,
    8.,8.5,9., /* band [5][*] */ 
    9.,9.375,9.75,
    9.75,10.,10.25,
    10.25,10.5,10.75,
    10.75,11.,11.25,
    11.25,11.5,11.75, /* band [10][*] */ 
    11.75,12.,12.25,
    12.25,12.625,13.,
    /* Old spectral bands until 4/22/93
    9.,9.5,10.,
    10.,10.5,11.,
    11.,11.5,12.,
    12.,12.5,13., */ 
    13.,14.,15.,
    15.,16.,17.,
    17.,18.,19., /* band [15][*] */ 
    19.,21.,23.,
    23.,25.,27.,
    27.,29.,31.,
    31.,33.,35.,
    35.,37.,39., /* band [20][*] */ 
    39.,41.,43.,
    43.,45.,47.,
    47.,49.,51.,
    51.,53.,55.,
    55.,57.,59., /* band [25][*] */ 
    59.,61.,63.,
    63.,65.,67.,
    67.,69.,71.,
    71.,73.,75.,
    75.,77.,79., /* band [30][*] */ 
    79.,81.,83.,
    0.,0.,0. /* band [num_LW_band+1] */
    };
  
  for(band=1;band<=NUM_IR_SPECTRAL_BANDS;band++){
    IR_bandwidth[band]=METERS_PER_MICRON*
      (IR_spectral_band[band][2]-IR_spectral_band[band][0]);
    IR_wavelength_microns=IR_spectral_band[band][1];
    IR_wavelength[band]=METERS_PER_MICRON*IR_wavelength_microns;

    /* get the indices of refraction */ 
    FORTRAN_refice(&IR_wavelength_microns,&foo_temperature,
	    IR_ice_real_idx+band,IR_ice_imag_idx+band);
  } /* end loop over bands */
} /* end Warren_LW_ice_init */ 
