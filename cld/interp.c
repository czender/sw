/* $Author: zender $
 * $Date$
 * $Id$
 * $Revision$
 * $Locker:  $
 * $RCSfile: interp.c,v $
 * $Source: /home/zender/cvs/cld/interp.c,v $
 * $Id$
 * $State: Exp $
 * */

/* Purpose: Initialization and interpolation routines for hardcoded data. */ 

/* $Log: not supported by cvs2svn $
/* Revision 1.1.1.1  1998-09-15 02:06:43  zender
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
 * Revision 4.0  1993/02/18  01:49:55  zender
 * this is pretty much the working comps version.  next additions
 * will be major:  Liou's habit params, real aggregation, eddy
 * diffusion.
 *
 * Revision 3.4  1993/01/23  03:55:36  zender
 * Not much different.  Added start from scratch growth -D72, tuned
 * iographs to make beautiful contour plots. About to add AGGREGATION.
 *
 * Revision 3.3  1993/01/08  02:19:41  zender
 * cloud top, base now determined by Concentration not by IWC.
 * Improved ESS selection, added IR optical depth routine, SW
 * properties alb., trans., abs.  Switched to interface fluxes.
 * Fixed Cray errors.
 *
 * Revision 3.2  1993/01/05  17:42:16  zender
 * Implemented my own in line LW flux divergence code to extract
 * crystal heating from.  Crashes on the Cray. About to tune the
 * apportion routine and cloudiness param (and Eulerian growth?)
 * so they work with non-monotonic growth rates vs. crystal size.
 *
 * Revision 3.1  1992/11/25  17:24:43  zender
 * implemented routine to automatically use/create a disk based
 * Mie parameter file whenver there are 40 crystal sizes. about
 * to eliminate trace gas absorption with cloud grid to get crystal
 * flux divergences.
 *
 * Revision 3.0  1992/11/16  02:41:01  zender
 * size/spectral-weighted crystal rad. heating rates are implemented.
 * old mean intensity method deadwood is stripped out.  must now
 * compensate flux divergences for intrinsic trace gas heating.
 * planck_avg_abs_eff not needed, just retained for a diagnostic.
 *
 * */

#include <stdio.h>              /* stderr, FILE, NULL, etc. */
#include <math.h>

/* my header files */
#include "defs.h"               /* YES, NO, FILESIZE, and more */
#include "globals.h"            /* externs announced in main() */

void init_ccm2_input
  (atmospheric_profile_type,
   ccm2_cloud_cover,
   ccm2_env_pressure,
   ccm2_env_temperature,
   ccm2_ice_path,
   ccm2_mmr_O3,
   ccm2_mmr_vapor,
   env_pressure,
   mmr_O3,
   
   ASNIR,
   ASVIS,
   AWNIR,
   AWVIS,
   day_of_year,
   frac_strng_zen_ang_srf,
   ground_pressure,
   latitude,
   skin_temp,
   snow_cover,
   surf_air_temp,
   surf_roughness,
   surf_type_flag,
   
   num_ccm2_level,
   num_layer,
   bottom_ccm2_interface_level,
   top_ccm2_interface_level)
float *ccm2_cloud_cover,*ccm2_env_pressure,*ccm2_env_temperature;
float *ccm2_ice_path,*ccm2_mmr_O3,*ccm2_mmr_vapor;
float *env_pressure,*mmr_O3;

float *ASNIR;
float *ASVIS;
float *AWNIR;
float *AWVIS;
float *day_of_year;
float *frac_strng_zen_ang_srf;
float *ground_pressure;
float *latitude;
float *skin_temp;
float *snow_cover;
float *surf_air_temp;
float *surf_roughness;
float *surf_type_flag;

int *bottom_ccm2_interface_level,*top_ccm2_interface_level;
int *num_ccm2_level;
int atmospheric_profile_type;
int num_layer;
{
  /* this routine shows reads in the ccm2 colmod style input file,
   finds out how much storage the external + cloud grid size arrays
   will take, allocates the arrays, and uses an interpolation 
   routine spline() to initialize the cirrus grid mmr_O3 array. 
   Returns all the ccm2 arrays with space set aside for the insertion
   of the cloud grid arrays. */ 

  void FORTRAN_spline();

  char ccm2_in_file[FILESIZE];
  char in_buf[FILESIZE];

  FILE *fp_ccm2_in;

  int ccm2_level;
  int int_foo;

  /* Initialize the parameters */ 
  if(atmospheric_profile_type == MID_LATITUDE_SUMMER){
    (void)strcpy(ccm2_in_file,"mls35_cld_crm.txt"); 
  }else if(atmospheric_profile_type == TROPICS){
    (void)strcpy(ccm2_in_file,"trp35_cld_crm.txt"); 
  } /* end else */

  if( (fp_ccm2_in = fopen( ccm2_in_file, "r")) == NULL) {
    (void)fprintf(stdout,"\nError in opening ccm2 file %s\n",ccm2_in_file);
    exit(1);
  } /* end if */
  
  (void)fscanf(fp_ccm2_in,"%s\n",in_buf);
  if(debug == 42)(void)fprintf(stdout,"%s\n",in_buf);
  (void)fscanf(fp_ccm2_in," %[mls]",in_buf);
  if(debug == 42)(void)fprintf(stdout,"%s",in_buf);
  (void)fscanf(fp_ccm2_in," %i",num_ccm2_level);
  if(debug == 42)(void)fprintf(stdout," %i ",*num_ccm2_level);
  (void)fscanf(fp_ccm2_in," %50c",in_buf);
/*  (void)fscanf(fp_ccm2_in," %[level ccm2]",in_buf);*/
  if(debug == 42)(void)fprintf(stdout,"%s\n",in_buf);
  (void)fscanf(fp_ccm2_in,"%s\n",in_buf);
  if(debug == 42)(void)fprintf(stdout,"%s\n",in_buf);
  (void)fscanf(fp_ccm2_in,"%f\n",day_of_year);
  if(debug == 42)(void)fprintf(stdout,"%f\t",*day_of_year);
  (void)fscanf(fp_ccm2_in,"%[day of year]",in_buf);
  if(debug == 42)(void)fprintf(stdout,"%s\n",in_buf);
  (void)fscanf(fp_ccm2_in,"%f\n",latitude);
  if(debug == 42)(void)fprintf(stdout,"%f\t",*latitude);
  (void)fscanf(fp_ccm2_in,"%s\n",in_buf);
  if(debug == 42)(void)fprintf(stdout,"%s\n",in_buf);
  (void)fscanf(fp_ccm2_in," %[level p(mb) t(k) h2ommr(g/g) o3mmr(g/g) cldcvr cldlwp(g/m2)] ",in_buf);
  if(debug == 42)(void)fprintf(stdout,"%s\n",in_buf);

  /* now for the real data stream */ 
  ccm2_level=*num_ccm2_level;
  while(ccm2_level > 1){
    (void)fscanf(fp_ccm2_in,"%i",&ccm2_level);
    (void)fscanf(fp_ccm2_in,"%f %f %f %f %f %f",
		 ccm2_env_pressure+ccm2_level,
		 ccm2_env_temperature+ccm2_level,
		 ccm2_mmr_vapor+ccm2_level,
		 ccm2_mmr_O3+ccm2_level,
		 ccm2_cloud_cover+ccm2_level,
		 ccm2_ice_path+ccm2_level);

    /* Convert the ccm2 pressures to Pascals while they are in my model,
     just for consistency.  They will all be reconverted to mb before 
     ccm2rad is called. */ 
    ccm2_env_pressure[ccm2_level]*=PASCALS_PER_MB;

    /* Convert the ccm2 cloud IWP's to kg/m^2 while they are in my model,
     just for consistency.  They will all be reconverted to g/m^2 before 
     ccm2 is called. */ 
    ccm2_ice_path[ccm2_level]*=KILOGRAMS_PER_GRAM;

  } /* end loop over ccm2_level */

  if(debug == 42){
    for(ccm2_level=*num_ccm2_level;ccm2_level>=1;ccm2_level--){
      (void)fprintf(stdout," %i %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n",
		    ccm2_level,
		    ccm2_env_pressure[ccm2_level],
		    ccm2_env_temperature[ccm2_level],
		    ccm2_mmr_vapor[ccm2_level],
		    ccm2_mmr_O3[ccm2_level],
		    ccm2_cloud_cover[ccm2_level],
		    ccm2_ice_path[ccm2_level]);
    } /* end loop over ccm2_level */
  } /* end debug */
  
  (void)fscanf(fp_ccm2_in,"%f\n",ground_pressure);
  if(debug == 42)(void)fprintf(stdout,"%f\n",*ground_pressure);
  (void)fscanf(fp_ccm2_in,"%f %[surface air temperature]",surf_air_temp,in_buf);
  if(debug == 42)(void)fprintf(stdout,"\t%f %s\n",*surf_air_temp,in_buf);
  (void)fscanf(fp_ccm2_in,"%f %[ground (skin) temperature]",skin_temp,in_buf);
  if(debug == 42)(void)fprintf(stdout,"\t%f %s\n",*skin_temp,in_buf);
  (void)fscanf(fp_ccm2_in,"%f %[surface type flag]",surf_type_flag,in_buf);
  if(debug == 42)(void)fprintf(stdout,"\t%f %s\n",*surf_type_flag,in_buf);
  (void)fscanf(fp_ccm2_in,"%f %[surface roughness]\n",surf_roughness,in_buf);
  if(debug == 42)(void)fprintf(stdout,"\t%f %s\n",*surf_roughness,in_buf);
  (void)fscanf(fp_ccm2_in,"%f %[SNOW COVER]",snow_cover,in_buf);
  if(debug == 42)(void)fprintf(stdout,"\t%f %s\n",*snow_cover,in_buf);
  (void)fscanf(fp_ccm2_in,"%f %[ASVIS]",ASVIS,in_buf);
  if(debug == 42)(void)fprintf(stdout,"\t%f %s\n",*ASVIS,in_buf);
  (void)fscanf(fp_ccm2_in,"%f %[ASNIR]",ASNIR,in_buf);
  if(debug == 42)(void)fprintf(stdout,"\t%f %s\n",*ASNIR,in_buf);
  (void)fscanf(fp_ccm2_in,"%f %[AWVIS]",AWVIS,in_buf);
  if(debug == 42)(void)fprintf(stdout,"\t%f %s\n",*AWVIS,in_buf);
  (void)fscanf(fp_ccm2_in,"%f %[AWNIR]",AWNIR,in_buf);
  if(debug == 42)(void)fprintf(stdout,"\t%f %s\n",*AWNIR,in_buf);
  (void)fscanf(fp_ccm2_in,"%f %[fraction strng zen ang srf]",frac_strng_zen_ang_srf,
	       in_buf);
  if(debug == 42)(void)fprintf(stdout,"\t%f %s\n",*frac_strng_zen_ang_srf,in_buf);
  (void)fclose(fp_ccm2_in);
  (void)fprintf(stdout,"Read environmental profile data from %s\n",ccm2_in_file);

  /* Now must splice my environmental data from within the cloud layers
     into the ccm2 environmental profile just read in. First find the two
     pressure levels which circumscribe the pressure in the computational
     grid of my cloud. */
  
  for(ccm2_level=1;ccm2_level<=*num_ccm2_level;ccm2_level++){
    if(ccm2_env_pressure[ccm2_level] < env_pressure[1])
      break;
  } /* end loop over ccm2_level */
  /* this indexes the last level beneath the cloud grid in a system 
   where index 1 is the nearest to the ground and e.g., 35 is TOA */ 
  *bottom_ccm2_interface_level=ccm2_level-1;

  for(ccm2_level=1;ccm2_level<=*num_ccm2_level;ccm2_level++){
    if(ccm2_env_pressure[ccm2_level] < env_pressure[num_layer])
      break;
  } /* end loop over ccm2_level */
  /* this indexes the first level above the cloud grid in a system 
   where index 1 is the nearest to the ground and e.g., 35 is TOA */ 
  *top_ccm2_interface_level=ccm2_level;

  /* Given the environmental ozone profile, interpolate a suitable profile
   over the cloud grid */
  int_foo=*num_ccm2_level-1;
  /* set debug_value = 2 for spline debugging */ 
  FORTRAN_spline(ccm2_env_pressure+1,ccm2_mmr_O3+1,&int_foo,env_pressure+1,
	 mmr_O3+1,&num_layer,&debug_value);

  /* First move the what's above the cloud to make space for the cloud */ 
  (void)memmove(
	ccm2_cloud_cover+*bottom_ccm2_interface_level+num_layer+1,ccm2_cloud_cover+*top_ccm2_interface_level,
	(*num_ccm2_level-*top_ccm2_interface_level+1)*sizeof(float));
  (void)memmove(
	ccm2_env_pressure+*bottom_ccm2_interface_level+num_layer+1,ccm2_env_pressure+*top_ccm2_interface_level,
	(*num_ccm2_level-*top_ccm2_interface_level+1)*sizeof(float));
  (void)memmove(
	ccm2_env_temperature+*bottom_ccm2_interface_level+num_layer+1,ccm2_env_temperature+*top_ccm2_interface_level,
	(*num_ccm2_level-*top_ccm2_interface_level+1)*sizeof(float));
  (void)memmove(
	ccm2_ice_path+*bottom_ccm2_interface_level+num_layer+1,ccm2_ice_path+*top_ccm2_interface_level,
	(*num_ccm2_level-*top_ccm2_interface_level+1)*sizeof(float));
  (void)memmove(
	ccm2_mmr_O3+*bottom_ccm2_interface_level+num_layer+1,ccm2_mmr_O3+*top_ccm2_interface_level,
	(*num_ccm2_level-*top_ccm2_interface_level+1)*sizeof(float));
  (void)memmove(
	ccm2_mmr_vapor+*bottom_ccm2_interface_level+num_layer+1,ccm2_mmr_vapor+*top_ccm2_interface_level,
	(*num_ccm2_level-*top_ccm2_interface_level+1)*sizeof(float));

} /* end ccm2input */ 

void write_ccm2_output
  (ccm2_cloud_cover,
   ccm2_env_pressure,
   ccm2_env_temperature,
   ccm2_ice_path,
   ccm2_mmr_O3,
   ccm2_mmr_vapor,
   
   ASNIR,
   ASVIS,
   AWNIR,
   AWVIS,
   day_of_year,
   frac_strng_zen_ang_srf,
   ground_pressure,
   latitude,
   skin_temp,
   snow_cover,
   surf_air_temp,
   surf_roughness,
   surf_type_flag,
   
   num_ccm2_level,
   num_layer
   )
float *ccm2_cloud_cover,*ccm2_env_pressure,*ccm2_env_temperature;
float *ccm2_ice_path,*ccm2_mmr_O3,*ccm2_mmr_vapor;

float ASNIR;
float ASVIS;
float AWNIR;
float AWVIS;
float day_of_year;
float frac_strng_zen_ang_srf;
float ground_pressure;
float latitude;
float skin_temp;
float snow_cover;
float surf_air_temp;
float surf_roughness;
float surf_type_flag;

int num_ccm2_level;
int num_layer;
{
  /* this routine writes out the cloud variables in a style compatible with
     ccm2 colmod style input files */ 
  
  char ccm2_out_file[FILESIZE];
  
  FILE *fp_ccm2_out;
  
  int layer;
  int num_total_layers;

  /* Initialize variables */ 
  (void)strcpy(ccm2_out_file,"cld_crm.txt"); 
  num_total_layers=num_layer+num_ccm2_level;

  if( (fp_ccm2_out = fopen( ccm2_out_file, "w")) == NULL) {
    (void)fprintf(stdout,"\nError in opening ccm2 file %s\n",ccm2_out_file);
    exit(1);
  } /* end if */

  (void)fprintf(fp_ccm2_out,"--------------------------------------\n");
  (void)fprintf(fp_ccm2_out," mls %i level ccm2\n",num_total_layers);
  (void)fprintf(fp_ccm2_out,"--------------------------------------\n");
  (void)fprintf(fp_ccm2_out,"%f\tday of year\n",day_of_year);
  (void)fprintf(fp_ccm2_out,"%f\tlatitude\n",latitude);
  (void)fprintf(fp_ccm2_out,"level p(mb) t(k) h2ommr(g/g) o3mmr(g/g) cldcvr cldlwp(g/m2)\n");

  /* now for the real data stream */ 
  for(layer=num_total_layers;layer>=1;layer--){
    (void)fprintf(fp_ccm2_out," %i %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n",
		  layer,ccm2_env_pressure[layer]*MB_PER_PASCAL,
		  ccm2_env_temperature[layer],
		  ccm2_mmr_vapor[layer],
		  ccm2_mmr_O3[layer],
		  ccm2_cloud_cover[layer],
		  ccm2_ice_path[layer]*GRAMS_PER_KILOGRAM);
  } /* end loop over layers */

  (void)fprintf(fp_ccm2_out,"%f\n",ground_pressure);
  (void)fprintf(fp_ccm2_out,"%f\tsurface air temperature\n",surf_air_temp);
  (void)fprintf(fp_ccm2_out,"%f\tground (skin) temperature\n",skin_temp);
  (void)fprintf(fp_ccm2_out,"%f\tsurface type flag\n",surf_type_flag);
  (void)fprintf(fp_ccm2_out,"%f\tsurface roughness\n",surf_roughness);
  (void)fprintf(fp_ccm2_out,"%f\tSNOW COVER\n",snow_cover);
  (void)fprintf(fp_ccm2_out,"%f\tASVIS\n",ASVIS);
  (void)fprintf(fp_ccm2_out,"%f\tASNIR\n",ASNIR);
  (void)fprintf(fp_ccm2_out,"%f\tAWVIS\n",AWVIS);
  (void)fprintf(fp_ccm2_out,"%f\tAWNIR\n",AWNIR);
  (void)fprintf(fp_ccm2_out,"%f\tfrac strng zen ang srf\n",frac_strng_zen_ang_srf);
  (void)fflush(fp_ccm2_out);
  (void)fclose(fp_ccm2_out);

  (void)fprintf(stdout,"Wrote %s file\n",ccm2_out_file);

} /* end ccm2output */ 

void init_temp_pres
  (float *altitude,
   float *env_pressure,
   float *env_temperature,
   int num_layer,
   int atmospheric_profile_type)
{
  void FORTRAN_spline
    (float */* std_atm_alt+1 */,
     float */* std_atm_temp+1 */,
     int */* &int_foo */,
     float */* altitude+1 */,
     float */* env_temperature+1 */,
     int */* &num_layer */,
     int */* &debug_value */);
  
  float *std_atm_alt;
  float *std_atm_pres;
  float *std_atm_temp;

  int int_foo;
  int layer;
  int num_std_atm_level;

  float mls_altitude_pressure_temperature[NUM_MLS_STD_ATM_LEVELS+2][3]={
    /* 45 degrees N July profile */ 
    /* format is altitude (m), pressure (Pa), temperature (K) */ 
    0.,0.,0., /* element [0] */ 

    4000.,6.280e4,273.57,
    4250.,6.087e4,272.04,
    4500.,5.899e4,270.51,
    4750.,5.715e4,268.98,

    5000.,5.536e4,267.45,
    5250.,5.362e4,265.92,
    5500.,5.192e4,264.39,
    5750.,5.027e4,262.86,

    6000.,4.866e4,261.33,
    6250.,4.709e4,259.70,
    6500.,4.557e4,258.07,
    6750.,4.408e4,256.44,
    
    7000.,4.264e4,254.81,
    7250.,4.123e4,253.17,
    7500.,3.986e4,251.54,
    7750.,3.853e4,249.91,

    8000.,3.724e4,248.28,
    8250.,3.598e4,246.65,
    8500.,3.475e4,245.03,
    8750.,3.356e4,243.40,

    9000.,3.240e4,241.77,
    9250.,3.128e4,240.15,
    9500.,3.018e4,238.52,
    9750.,2.912e4,236.90,

    10000.,2.809e4,235.27,
    10250.,2.709e4,233.65,
    10500.,2.612e4,232.02,
    10750.,2.517e4,230.40,

    11000.,2.426e4,228.77,
    11500.,2.250e4,225.53,

    12000.,2.086e4,222.30,
    12500.,1.931e4,219.06,

    13000.,1.786e4,215.82,
    13500.,1.650e4,215.65,
    
    14000.,1.525e4,215.65,
    14500.,1.409e4,215.65,

    15000.,1.525e4,215.65,
    15500.,1.409e4,215.65,
    
    0.,0.,0. /* element [NUM_MLS_STD_ATM_LEVELS+1] */
    };

  float tropical_altitude_pressure_temperature[NUM_TROPICS_STD_ATM_LEVELS+2][3]={
    /* Tropical ICRCCM profile */ 
    /* format is altitude (m), pressure (Pa), temperature (K) */ 
    0.,0.,0., /* element [0] */ 

    56.50,100650.,299.66,
    201.35,99000.,298.79,
    379.89,97000.,297.72,
    561.81,95000.,296.64,
    747.26,93000.,295.53, /* element [5] */

    936.39,91000.,294.41,
    1129.35,89000.,293.27,
    1326.31,87000.,292.12,
    1527.46,85000.,290.96,
    1733.00,83000.,289.82, /* element [10] */

    1943.16,81000.,288.71,
    2158.19,79000.,287.66,
    2378.35,77000.,286.65,
    2603.92,75000.,285.64,
    2835.16,73000.,284.57, /* element [15] */

    3072.34,71000.,283.37,
    3315.73,69000.,282.01,
    3565.62,67000.,280.50,
    3822.38,65000.,278.89,
    4086.39,63000.,277.20, /* element [20] */

    4358.10,61000.,275.46,
    4637.98,59000.,273.66,
    4926.59,57000.,271.82,
    5224.52,55000.,269.92,
    5532.42,53000.,267.96, /* element [25] */

    5851.03,51000.,265.94,
    6181.16,49000.,263.86,
    6523.75,47000.,261.71,
    6879.82,45000.,259.48,
    7250.54,43000.,257.17, /* element [30] */

    7637.25,41000.,254.78,
    8041.48,39000.,252.28,
    8464.99,37000.,249.68,
    8909.84,35000.,246.97,
    9378.42,33000.,244.13, /* element [35] */

    9873.57,31000.,241.14,
    10398.69,29000.,238.00,
    10957.88,27000.,234.68,
    11556.15,25000.,231.15,
    12199.74,23000.,227.39, /* element [40] */

    12896.53,21000.,223.35,
    13656.73,19000.,218.99,
    14493.83,17000.,214.24,
    15426.27,15000.,209.02,
    16480.14,13000.,203.20, /* element [45] */

    17695.00,11000.,196.89,
    18922.04,9288.,195.75,
    20043.52,7967.,199.18,
    21173.16,6833.,202.81,
    22311.14,5861.,206.51, /* element [50] */

    23457.50,5027.,210.25,
    24611.36,4311.,213.60,
    25770.57,3698.,215.99,
    28108.47,2720.,220.34,
    29282.41,2333.,222.55, /* element [55] */

    30461.69,2001.,224.78,
    31645.44,1716.,227.03,
    32834.84,1472.,229.30,
    34029.24,1263.,231.60,
    35228.23,1083.,233.92, /* element [60] */

    36433.32,929.,236.27,
    37643.45,797.,238.64,
    38858.98,683.,241.04,
    40079.98,586.,243.46,
    41306.43,503.,245.90, /* element [65] */

    42538.41,431.,248.37,
    43775.97,370.,250.87,
    45019.20,317.,253.39,
    46268.10,272.,255.94,
    47522.85,233.,258.52, /* element [70] */

    48783.75,200.,261.10,
    50049.61,172.,263.59,
    51320.52,147.,265.57,
    52594.36,126.,266.96,
    53870.33,108.,268.20, /* element [75] */

    55149.36,93.,269.37,
    56430.21,80.,270.06,
    57708.64,68.,268.90,
    58979.65,59.,265.97,
    60243.62,50.,262.64, /* element [80] */

    61501.37,43.,259.30,
    62749.45,37.,256.00,
    63990.67,32.,252.75,
    65227.06,27.,249.53,
    66454.28,23.,246.35, /* element [85] */

    67674.69,20.,243.21,
    68889.26,17.,240.11,
    70094.78,15.,237.05,
    71291.63,13.,234.03,
    72486.64,11.,231.05, /* element [90] */

    73672.12,9.,228.10,

    0.,0.,0. /* element [NUM_TROPICS_STD_ATM_LEVELS+1] */
    };
  
  if(atmospheric_profile_type == MID_LATITUDE_SUMMER){
    num_std_atm_level=NUM_MLS_STD_ATM_LEVELS;
  }else if(atmospheric_profile_type == TROPICS){
    num_std_atm_level=NUM_TROPICS_STD_ATM_LEVELS;
  } /* end else */

  if(
     ((std_atm_alt=(float *)malloc((num_std_atm_level+2)*sizeof(float))) == NULL ) ||
     ((std_atm_pres=(float *)malloc((num_std_atm_level+2)*sizeof(float))) == NULL ) ||
     ((std_atm_temp=(float *)malloc((num_std_atm_level+2)*sizeof(float))) == NULL ) ||
     False ){
    (void)fprintf(stdout,"Unable to allocate array in interp()\n");
    exit(1);
  } /* end if */

  if(atmospheric_profile_type == MID_LATITUDE_SUMMER){
    for(layer=1;layer<=num_std_atm_level;layer++){
      std_atm_alt[layer]=mls_altitude_pressure_temperature[layer][0];
      std_atm_pres[layer]=mls_altitude_pressure_temperature[layer][1];
      std_atm_temp[layer]=mls_altitude_pressure_temperature[layer][2];
    } /* end loop over layers */
  }else if(atmospheric_profile_type == TROPICS){
    for(layer=1;layer<=num_std_atm_level;layer++){
      std_atm_alt[layer]=tropical_altitude_pressure_temperature[layer][0];
      std_atm_pres[layer]=tropical_altitude_pressure_temperature[layer][1];
      std_atm_temp[layer]=tropical_altitude_pressure_temperature[layer][2];
    } /* end loop over layers */
  } /* end else */

  if(debug == 57){
    (void)fprintf(stdout,"Initial Thermodynamic Grid:\n");
    (void)fprintf(stdout,"layer\talt.\tTemp\tpres\n");
    (void)fprintf(stdout,"\t km\t K\t mbar\n");
    for(layer=1;layer<=num_std_atm_level;layer++){
      (void)fprintf(stdout,"%i\t%.3f\t%.2f\t%.1f\n",layer,std_atm_alt[layer]*KILOMETERS_PER_METER,std_atm_temp[layer],std_atm_pres[layer]*MB_PER_PASCAL);
    } /* end loop over layers */
  } /* end debug */
  
  int_foo=num_std_atm_level-1;
  /* set debug_value = 2 for spline debugging */ 
  
  FORTRAN_spline(std_atm_alt+1,std_atm_temp+1,&int_foo,altitude+1,
	  env_temperature+1,&num_layer,&debug_value);
  FORTRAN_spline(std_atm_alt+1,std_atm_pres+1,&int_foo,altitude+1,
	  env_pressure+1,&num_layer,&debug_value);
  
  if(debug == 57){
    (void)fprintf(stdout,"\n");
    (void)fprintf(stdout,"Spline fitting results spline.f:\n");
    (void)fprintf(stdout,"layer\talt.\tTemp\tpres\n");
    (void)fprintf(stdout,"\t km\t K\t mbar\n");
    for(layer=1;layer<=num_layer;layer++){
      (void)fprintf(stdout,"%i\t%.3f\t%.2f\t%.1f\n",layer,altitude[layer]*KILOMETERS_PER_METER,env_temperature[layer],env_pressure[layer]*MB_PER_PASCAL);
    } /* end loop over layers */
  } /* end debug */
  
  free(std_atm_alt);
  free(std_atm_temp);
  free(std_atm_pres);
  
} /* end init_temp_pres() */

#define NUM_LW_HEX_BANDS 12
void Liou_IR_fudge
  (int num_LW_band,
   int num_size,
   float **LW_absorption_x_sec, /* returned scaled absorption */ 
   float **LW_scattering_x_sec, /* returned scaled scattering */ 
   float **LW_extinction_x_sec, /* returned scaled extinction */ 
   float **LW_asymmetry_param) /* returned scaled asymmetry */ 
{
  int band;
  int size;

  int CCM2_IR_band2hex_band[NUM_IR_SPECTRAL_BANDS+2]={
    0, /* band [0] */ 
    1,
    2,
    3,
    4,
    5, /* band [5] */ 
    6,
    6,
    7,
    7,
    7, /* band [10] */ 
    7,
    8,
    8,
    9,
    9, /* band [15] */ 
    10,
    10,
    11,
    11,
    12, /* band [20] */ 
    12,
    12,
    12,
    12,
    12, /* band [25] */ 
    12,
    12,
    12,
    12,
    12, /* band [30] */ 
    12,
    0 /* element [NUM_IR_SPECTRAL_BANDS+1] */ 
    };

  float hex_data_wavelength_low_microns[NUM_LW_HEX_BANDS+2]={
    /* these are the band centers in microns */ 
    0., /* element [0] */ 
    4.0,5.26,5.88,7.14,8.0,
    9.09,10.2,12.5,14.9,18.5,
    25.0,35.7,
    0. /* element [NUM_LW_HEX_BANDS+1] */ 
    };
  float hex_data_wavelength_microns[NUM_LW_HEX_BANDS+2]={
    /* these are the band centers in microns */ 
    0., /* element [0] */ 
    4.9,5.6,6.5,7.6,8.5,
    9.6,11.3,13.7,16.6,21.5,
    30.0,70.0,
    0. /* element [NUM_LW_HEX_BANDS+1] */ 
    };
  float hex_data_wavelength_high_microns[NUM_LW_HEX_BANDS+2]={
    /* these are the band centers in microns */ 
    0., /* element [0] */ 
    5.26,5.88,7.14,8.0,9.09,
    10.2,12.5,14.9,18.5,25.0,
    35.7,200.,
    0. /* element [NUM_LW_HEX_BANDS+1] */ 
    };

  float hex_LW_ext_coeff[NUM_LW_HEX_BANDS+2]={
    /* extinction coefficient in units of km^-1 */ 
    0., /* element [0] */ 
    .4075,.4366,.5123,.4923,.4587,
    .4158,.4559,.4596,.4433,.4286,
    .4258,.4184,
    0. /* element [NUM_LW_HEX_BANDS+1] */ 
    };
  float hex_LW_ss_albedo[NUM_LW_HEX_BANDS+2]={
    0., /* element [0] */ 
    .4285,.6586,.6962,.5555,.5161,
    .4937,.6133,.6191,.5804,.5563,
    .5915,.6506,
    0. /* element [NUM_LW_HEX_BANDS+1] */ 
    };
  float hex_LW_asym_param[NUM_LW_HEX_BANDS+2]={
    0., /* element [0] */ 
    .7499,.8668,.8089,.8722,.9044,
    .9480,.9398,.9215,.9279,.9372,
    .9278,.9002,
    0. /* element [NUM_LW_HEX_BANDS+1] */ 
    };

  float sphere_LW_ext_coeff[NUM_LW_HEX_BANDS+2]={
    /* extinction coefficient in units of km^-1 */ 
    0., /* element [0] */ 
    .4383,.4760,.5208,.4945,.4599,
    .4256,.4693,.4675,.4473,.4361,
    .4298,.4290,
    0. /* element [NUM_LW_HEX_BANDS+1] */ 
    };
  float sphere_LW_ss_albedo[NUM_LW_HEX_BANDS+2]={
    0., /* element [0] */ 
    .4296,.6516,.6716,.5327,.4995,
    .4827,.5834,.5838,.5389,.5180,
    .5509,.6219,
    0. /* element [NUM_LW_HEX_BANDS+1] */ 
    };
  float sphere_LW_asym_param[NUM_LW_HEX_BANDS+2]={
    0., /* element [0] */ 
    .7590,.8694,.8083,.8813,.9094,
    .9519,.9428,.9268,.9347,.9467,
    .9355,.9123,
    0. /* element [NUM_LW_HEX_BANDS+1] */ 
    };
    
  for(band=1;band<=num_LW_band;band++){
    for(size=1;size<=num_size;size++){
      LW_absorption_x_sec[band][size]*=
	(hex_LW_ext_coeff[CCM2_IR_band2hex_band[band]]*
	 (1.-hex_LW_ss_albedo[CCM2_IR_band2hex_band[band]]))/
	   (sphere_LW_ext_coeff[CCM2_IR_band2hex_band[band]]*
	    (1.-sphere_LW_ss_albedo[CCM2_IR_band2hex_band[band]]));
      LW_scattering_x_sec[band][size]*=
	(hex_LW_ext_coeff[CCM2_IR_band2hex_band[band]]*
	 hex_LW_ss_albedo[CCM2_IR_band2hex_band[band]])/
	   (sphere_LW_ext_coeff[CCM2_IR_band2hex_band[band]]*
	    sphere_LW_ss_albedo[CCM2_IR_band2hex_band[band]]);
      LW_extinction_x_sec[band][size]*=
	hex_LW_ext_coeff[CCM2_IR_band2hex_band[band]]/
	  sphere_LW_ext_coeff[CCM2_IR_band2hex_band[band]];
      LW_asymmetry_param[band][size]*=
	hex_LW_asym_param[CCM2_IR_band2hex_band[band]]/
	  sphere_LW_asym_param[CCM2_IR_band2hex_band[band]];
    } /* end loop over sizes */
  } /* end loop over bands */
  
  /* Ensure the scattering is never greater than the extinction and 
     compute the absorption cross-section and make sure the asymmetry
     parameter is always between plus and minus one. */ 
  for(band=1;band<=num_LW_band;band++){
    for(size=1;size<=num_size;size++){
      LW_extinction_x_sec[band][size]=
	(LW_extinction_x_sec[band][size] > LW_scattering_x_sec[band][size]) ?
	  LW_extinction_x_sec[band][size] : LW_scattering_x_sec[band][size];
      LW_absorption_x_sec[band][size]=
	LW_extinction_x_sec[band][size]-LW_scattering_x_sec[band][size];
      LW_asymmetry_param[band][size]=
	(LW_asymmetry_param[band][size] <= 1.) ? 
	  LW_asymmetry_param[band][size] : 1.;
      LW_asymmetry_param[band][size]=
	(LW_asymmetry_param[band][size] >= -1.) ? 
	  LW_asymmetry_param[band][size] : -1.;
    } /* end loop over sizes */
  } /* end loop over bands */
  
} /* end Liou_IR_fudge */ 

#define NUM_HEX_BANDS 5
#define NUM_HEX_SIZES 5
void Liou_fudge
  (float *wavelength,int num_band,
   float *hex_column_length, int num_size,
   float *hex_equiv_radius,
   float *hex_equiv_rad_squared, 
   float **absorption_x_sec, /* returned scaled absorption */ 
   float **scattering_x_sec, /* returned scaled scattering */ 
   float **extinction_x_sec, /* returned scaled extinction */ 
   float **asymmetry_param) /* returned scaled asymmetry */ 
/* NB: the last four optical properties must be supplie as equal area
   sphere properties and will be returned as hexagonal properties. */ 
{
  void FORTRAN_callbh();
  void FORTRAN_refice();

  Boolean BAND_AVERAGED_INDICES=False;

  float ice_real_idx[NUM_HEX_BANDS+2];
  float ice_imag_idx[NUM_HEX_BANDS+2];

  float hex_ext_x_sec[NUM_HEX_BANDS+2][NUM_HEX_SIZES+2];
  float hex_scat_x_sec[NUM_HEX_BANDS+2][NUM_HEX_SIZES+2];
  float hex_abs_x_sec[NUM_HEX_BANDS+2][NUM_HEX_SIZES+2];

  float sphere_abs_eff[NUM_HEX_BANDS+2][NUM_HEX_SIZES+2];
  float sphere_scat_eff[NUM_HEX_BANDS+2][NUM_HEX_SIZES+2];
  float sphere_ext_eff[NUM_HEX_BANDS+2][NUM_HEX_SIZES+2];
  float sphere_size_param[NUM_HEX_BANDS+2][NUM_HEX_SIZES+2];

  float sphere_ext_x_sec[NUM_HEX_BANDS+2][NUM_HEX_SIZES+2];
  float sphere_scat_x_sec[NUM_HEX_BANDS+2][NUM_HEX_SIZES+2];
  float sphere_abs_x_sec[NUM_HEX_BANDS+2][NUM_HEX_SIZES+2];
  float sphere_asym_param[NUM_HEX_BANDS+2][NUM_HEX_SIZES+2];

  float sphere2hex_abs_factor[NUM_HEX_BANDS+2][NUM_HEX_SIZES+2];
  float sphere2hex_asym_factor[NUM_HEX_BANDS+2][NUM_HEX_SIZES+2];
  float sphere2hex_ext_factor[NUM_HEX_BANDS+2][NUM_HEX_SIZES+2];
  float sphere2hex_scat_factor[NUM_HEX_BANDS+2][NUM_HEX_SIZES+2];

  float hex_data_aspect_ratio[NUM_HEX_SIZES+2];
  float hex_data_crystal_mass[NUM_HEX_SIZES+2];
  float hex_data_diam[NUM_HEX_SIZES+2];
  float hex_data_equiv_radius[NUM_HEX_SIZES+2];
  float hex_data_equiv_rad_squared[NUM_HEX_SIZES+2];
  float hex_data_length[NUM_HEX_SIZES+2];
  float hex_data_surface_area[NUM_HEX_SIZES+2];
  float hex_data_volume[NUM_HEX_SIZES+2];
  float hex_data_wavelength[NUM_HEX_BANDS+2];
  float hex_data_x_sec_area[NUM_HEX_SIZES+2];

  float basal_plane_area;
  float basal_plane_radius;
  float float_foo;
  float float_foo2;
  float float_foo3;
  float foo_temperature = 230.; 
  float hex_data_equiv_rad_microns;
  float pi_r_squared;

  #define NUM_INDEX_REFRACTION_QUAD_POINTS 20
  float delta_quad_wavelength_microns;
  float imag_idx;
  float quad_wavelength_microns;
  float real_idx;
  float total_imag_idx;
  float total_real_idx;

  int band;
  int hex_size;
  int layer;
  int num_angles = 2;
  int quad_point;
  int size;

  int CCM2_band2hex_band[NUM_CCM2_SPECTRAL_BANDS+2]={
    0, /* element [0] */ 
    1,1,1,1,1,1,1,1, /* bands 1--8 */ 
    2,2, /* bands 9--10 */ 
    3, /* band 11 */ 
    4,4,4, /* bands 12-14 */ 
    5,5,5,5, /* bands 15--18 */ 
    0 /* element [NUM_CCM2_SPECTRAL_BANDS+1] */ 
    };
  float hex_data_wavelength_low_microns[NUM_HEX_BANDS+2]={
    /* these are the band centers in microns */ 
    0., /* element [0] */ 
    .4,.7,1.3,1.9,2.5,
    0. /* element [NUM_HEX_BANDS+1] */ 
    };
  float hex_data_wavelength_microns[NUM_HEX_BANDS+2]={
    /* these are the band centers in microns */ 
    0., /* element [0] */ 
    .55,
    1.0,
    1.6,
    2.2,
    3.0,
    0. /* element [NUM_HEX_BANDS+1] */ 
    };
  float hex_data_wavelength_high_microns[NUM_HEX_BANDS+2]={
    /* these are the band centers in microns */ 
    0., /* element [0] */ 
    .7,1.3,1.9,2.5,3.5,
    0. /* element [NUM_HEX_BANDS+1] */ 
    };
  float ice_data_real_idx[NUM_HEX_BANDS+2]={
    /* these are the band centers in microns */ 
    0., /* element [0] */ 
    1.311,1.302,1.290,1.263,1.242,
    0. /* element [NUM_HEX_BANDS+1] */ 
    };
  float ice_data_imag_idx[NUM_HEX_BANDS+2]={
    /* these are the band centers in microns */ 
    0., /* element [0] */ 
    3.110e-9,1.931e-6,2.128e-4,7.997e-4,1.424e-1,
    0. /* element [NUM_HEX_BANDS+1] */ 
    };
  float hex_data_length_microns[NUM_HEX_SIZES+2]={
    0., /* element [0] */ 
    20.,
    50.,
    120.,
    300.,
    750.,
    0. /* element [NUM_HEX_SIZES+1] */ 
    };
  float hex_data_diam_microns[NUM_HEX_SIZES+2]={
    0., /* element [0] */ 
    20.,
    40.,
    60.,
    100.,
    160.,
    0. /* element [NUM_HEX_SIZES+1] */ 
    };
  float hex_ext_x_sec_CGS[NUM_HEX_BANDS+2][NUM_HEX_SIZES+2]={
    /* extinction cross section in units of cm^2 */ 
    0.,0.,0.,0.,0.,0.,0., /* elements [0][*] */ 
    0.,8.5981e-6,4.0392e-5,1.3138e-4,5.1495e-4,1.9663e-3,0.,
    0.,8.5981e-6,4.0392e-5,1.3138e-4,5.1495e-4,1.9663e-3,0.,
    0.,8.5981e-6,4.0392e-5,1.3138e-4,5.1495e-4,1.9663e-3,0.,
    0.,8.5981e-6,4.0392e-5,1.3138e-4,5.1495e-4,1.9663e-3,0.,
    0.,8.5981e-6,4.0392e-5,1.3138e-4,5.1495e-4,1.9663e-3,0.,
    0.,0.,0.,0.,0.,0.,0. /* elements [NUM_HEX_BANDS+1][*] */ 
    };
  float hex_ss_albedo[NUM_HEX_BANDS+2][NUM_HEX_SIZES+2]={
    0.,0.,0.,0.,0.,0.,0., /* elements [0][*] */ 
    0.,0.999999278,0.999998462,0.999997444,0.999995471,0.999992424,0.,
    0.,0.99975683 ,0.999482194,0.999140117,0.998478633,0.997460032,0.,
    0.,0.983894634,0.966534453,0.946369974,0.911990537,0.870501324,0.,
    0.,0.959152259,0.918214765,0.875399152,0.812903237,0.749717993,0.,
    0.,0.53411942 ,0.53138714 ,0.53089312 ,0.53070594 ,0.53064994 ,0.,
    0.,0.,0.,0.,0.,0.,0. /* elements [NUM_HEX_BANDS+1][*] */ 
    };
  float hex_asym_param[NUM_HEX_BANDS+2][NUM_HEX_SIZES+2]={
    0.,0.,0.,0.,0.,0.,0., /* elements [0][*] */ 
    0.,0.770403,0.777955,0.815471,0.842867,0.859161,0.,
    0.,0.776228,0.784654,0.821623,0.849228,0.865001,0.,
    0.,0.789597,0.805137,0.846257,0.880174,0.903669,0.,
    0.,0.818510,0.841624,0.882939,0.916478,0.937994,0.,
    0.,0.956986,0.967520,0.971743,0.973668,0.974627,0.,
    0.,0.,0.,0.,0.,0.,0. /* elements [NUM_HEX_BANDS+1][*] */ 
    };

  /* Convert all units to MKS */ 
  for(band=1;band<=NUM_HEX_BANDS;band++){
    hex_data_wavelength[band]=hex_data_wavelength_microns[band]*METERS_PER_MICRON;
    for(size=1;size<=NUM_HEX_SIZES;size++){
      hex_ext_x_sec[band][size]=hex_ext_x_sec_CGS[band][size]*
	METERS_PER_CENTIMETER*METERS_PER_CENTIMETER;
      hex_scat_x_sec[band][size]=
	hex_ss_albedo[band][size]*hex_ext_x_sec[band][size];
      hex_abs_x_sec[band][size]=
	hex_ext_x_sec[band][size]-hex_scat_x_sec[band][size];
    } /* end loop over sizes */
  } /* end loop over bands */

  for(size=1;size<=NUM_HEX_SIZES;size++){
    hex_data_length[size]=hex_data_length_microns[size]*METERS_PER_MICRON;
    hex_data_diam[size]=hex_data_diam_microns[size]*METERS_PER_MICRON;
  } /* end loop over sizes */

  /* Figure out what the geometric cross-sectional area, volume, and mass
     of each of the hexagonal columns in the Liou data is */ 
  for(size=1;size<=NUM_HEX_SIZES;size++){
    basal_plane_radius=hex_data_diam[size]/2.;
    float_foo=3.*sqrt(3.)/2.;
    basal_plane_area=float_foo*basal_plane_radius*basal_plane_radius;
    hex_data_surface_area[size]=2.*basal_plane_area+
      6.*basal_plane_radius*hex_data_length[size];
    hex_data_equiv_rad_squared[size]=.25*hex_data_surface_area[size]/M_PI;
    hex_data_equiv_radius[size]=sqrt(hex_data_equiv_rad_squared[size]);

    /* note that the following is ebert and curry's approximation, 
     takano and liou do this correctly */ 
    hex_data_x_sec_area[size]=hex_data_surface_area[size]/4.;
    hex_data_volume[size]=basal_plane_area*hex_data_length[size];
    hex_data_crystal_mass[size]=density_of_ice*hex_data_volume[size];
    hex_data_aspect_ratio[size]=hex_data_length[size]/hex_data_diam[size];
  } /* end loop over sizes */

  if(debug == 75){
    (void)fprintf(stdout,"sizebin\tLength\tDiam\tAspect\tr_s\tr_s^2\tSurf A\tX-sec A\tMass\n");
    (void)fprintf(stdout,"\tmicrons\tmicrons\t\tmicrons\tmu^2\tmu^2\tmu^2\tng\n");
    for(size=1;size<=NUM_HEX_SIZES;size++){
      (void)fprintf
	(stdout,"%i\t%.1f\t%.1f\t%.1f\t%.1f\t%g\t%.0f\t%.0f\t%.0f\n",
	 size,
	 hex_data_length[size]*MICRONS_PER_METER,
	 hex_data_diam[size]*MICRONS_PER_METER,
	 hex_data_aspect_ratio[size],
	 hex_data_equiv_radius[size]*MICRONS_PER_METER,
	 hex_data_equiv_rad_squared[size]*SQUARE_MICRONS_PER_SQUARE_METER,
	 hex_data_surface_area[size]*SQUARE_MICRONS_PER_SQUARE_METER,
	 hex_data_x_sec_area[size]*SQUARE_MICRONS_PER_SQUARE_METER,
	 hex_data_crystal_mass[size]*NANOGRAMS_PER_KILOGRAM);
    } /* end loop over sizes */
  } /* end debug */

  if(BAND_AVERAGED_INDICES){
    /* for each band find a mean index of refraction by averaging over 
       NUM_INDEX_REFRACTION_QUAD_POINTS wavelengths in the band */ 
    for(band=1;band<=NUM_HEX_BANDS;band++){
      delta_quad_wavelength_microns=
	(hex_data_wavelength_high_microns[band]-
	 hex_data_wavelength_low_microns[band])/
	   ((float)NUM_INDEX_REFRACTION_QUAD_POINTS-1.);
      /* Reset the accumulating variables */ 
      total_real_idx=0.;
      total_imag_idx=0.;
      quad_wavelength_microns=hex_data_wavelength_low_microns[band];
      for(quad_point=1;quad_point<=NUM_INDEX_REFRACTION_QUAD_POINTS;quad_point++){
	FORTRAN_refice(&quad_wavelength_microns,&foo_temperature,
		&real_idx,&imag_idx);
	total_real_idx+=real_idx;
	total_imag_idx+=imag_idx;
	if(debug == 78){
	  (void)fprintf(stdout,"%i\t%10.3f\t%10.3f\t%e\t%10.3f\t%e\n",
			quad_point,
			quad_wavelength_microns,
			real_idx,
			imag_idx,
			total_real_idx/quad_point,
			total_imag_idx/quad_point);
	} /* end debug */
	quad_wavelength_microns+=delta_quad_wavelength_microns;
	
      } /* end loop over quad_points */
      ice_real_idx[band]=total_real_idx/NUM_INDEX_REFRACTION_QUAD_POINTS;
      ice_imag_idx[band]=total_imag_idx/NUM_INDEX_REFRACTION_QUAD_POINTS;
    } /* end loop over bands */
  } /* end if BAND_AVERAGED_INDICES */
  
  if(!BAND_AVERAGED_INDICES){
    /* CAUTION! CAUTION! CAUTION! CAUTION! CAUTION! CAUTION! CAUTION! */ 
    /* This moots the band-averaging of the indices of refraction above */ 
    for(band=1;band<=NUM_HEX_BANDS;band++){
      ice_real_idx[band]=ice_data_real_idx[band];
      ice_imag_idx[band]=ice_data_imag_idx[band];
    } /* end loop over bands */
  } /* end if !BAND_AVERAGED_INDICES */
  
  /* find out what the optical properties for the corresponding 
     equivalent radius spheres are: if the computationally expensive
     parameters have been saved to disk, retrieve them; if not,
     then compute and store them for next time */ 

  for(band=1;band<=NUM_HEX_BANDS;band++){
    for(size=1;size<=NUM_HEX_SIZES;size++){
      hex_data_equiv_rad_microns=hex_data_equiv_radius[size]*MICRONS_PER_METER;
      
      /* The wavelength and radius must be in the same units so why not
	 scale them to microns in case it helps with rounding errors? */ 
      FORTRAN_callbh
	(&num_angles,
	 hex_data_wavelength_microns+band,
	 &hex_data_equiv_rad_microns,
	 ice_real_idx+band,
	 ice_imag_idx+band,
	 sphere_size_param[band]+size,
	 sphere_scat_eff[band]+size,
	 sphere_ext_eff[band]+size,
	 sphere_asym_param[band]+size,
	 &debug);
      
      sphere_ext_eff[band][size]=
	(sphere_ext_eff[band][size] > sphere_scat_eff[band][size]) ? 
	  sphere_ext_eff[band][size] : sphere_scat_eff[band][size];
      sphere_abs_eff[band][size]=
	sphere_ext_eff[band][size]-sphere_scat_eff[band][size];
    } /* end loop over sizes */
  } /* end loop over bands */
  
  if(debug == 76){
    for(band=1;band<=NUM_HEX_BANDS;band++){
      (void)fprintf(stdout,"\nwband\tlambda\tsbin\tequiv_r\tnreal\tnimag\tX\tQext\tQsca\tQabs\tg\tomega\n");
      for(size=1;size<=NUM_HEX_SIZES;size++){
	(void)fprintf
	  (stdout,"%i\t%.3f\t%i\t%.2f\t%.3f\t%7.1e\t%.1f\t%.4f\t%.4f\t%.4f\t%.4f\t%.6f\n",
	   band,
	   hex_data_wavelength[band]*MICRONS_PER_METER,
	   size,
	   hex_data_equiv_radius[size]*MICRONS_PER_METER,
	   ice_real_idx[band],
	   ice_imag_idx[band],
	   sphere_size_param[band][size],
	   sphere_ext_eff[band][size],
	   sphere_scat_eff[band][size],
	   sphere_abs_eff[band][size],
	   sphere_asym_param[band][size],
	   sphere_scat_eff[band][size]/sphere_ext_eff[band][size]);
      } /* end loop over sizes */
    } /* end loop over bands */
  } /* end debug */

  /* find the optical cross-sections */ 
  for(size=1;size<=NUM_HEX_SIZES;size++){
    pi_r_squared=M_PI*hex_data_equiv_rad_squared[size];
    for(band=1;band<=NUM_HEX_BANDS;band++){
      sphere_abs_x_sec[band][size]=sphere_abs_eff[band][size]*pi_r_squared;
      sphere_scat_x_sec[band][size]=sphere_scat_eff[band][size]*pi_r_squared;
      sphere_ext_x_sec[band][size]=sphere_ext_eff[band][size]*pi_r_squared;
    } /* end loop over bands */
  } /* end loop over sizes */
  
  /* CAUTION! CAUTION! CAUTION! CAUTION! CAUTION! CAUTION! CAUTION! */ 
  sphere_abs_x_sec[1][1]*=891.8e-12/sphere_ext_x_sec[1][1];
  sphere_scat_x_sec[1][1]*=891.8e-12/sphere_ext_x_sec[1][1];
  sphere_ext_x_sec[1][1]*=891.8e-12/sphere_ext_x_sec[1][1];

  if(debug == 76){
    (void)fprintf(stdout,"X-sections and individual optical properties of crystals:\n");
    for(band=1;band<=NUM_HEX_BANDS;band++){
      (void)fprintf
	(stdout,"%10s%10s%10s%10s%12s%12s%12s%12s%13s%13s%10s%10s%10s%10s\n",
	 "wavebin","wavelength","sizebin","hex length","Hex ext-xs",
	 "Sph ext-xs","Hex sc-xs","Sph sc-xs","Hex abs-xs","Sph abs-xs",
	 "Hex co-alb","Sph co-alb","Hex asym","Sph asym");
      (void)fprintf
	(stdout,"%10s%10s%10s%10s%12s%12s%12s%12s%13s%13s%10s%10s%10s%10s\n",
	 "index","microns","index","microns","micron^2",
	 "micron^2","micron^2","micron^2","micron^2","micron^2",
	 "","","","");
      for(size=1;size<=NUM_HEX_SIZES;size++){
	(void)fprintf
	  (stdout,"%10i%10.3f%10i%10.3f%12.3f%12.3f%12.3f%12.3f%13.5f%13.5f%10.7f%10.7f%10.3f%10.3f\n",
	   band,
	   hex_data_wavelength[band]*MICRONS_PER_METER,
	   size,
	   hex_data_length[size]*MICRONS_PER_METER,
	   hex_ext_x_sec[band][size]*SQUARE_MICRONS_PER_SQUARE_METER,
	   sphere_ext_x_sec[band][size]*SQUARE_MICRONS_PER_SQUARE_METER,
	   hex_scat_x_sec[band][size]*SQUARE_MICRONS_PER_SQUARE_METER,
	   sphere_scat_x_sec[band][size]*SQUARE_MICRONS_PER_SQUARE_METER,
	   hex_abs_x_sec[band][size]*SQUARE_MICRONS_PER_SQUARE_METER,
	   sphere_abs_x_sec[band][size]*SQUARE_MICRONS_PER_SQUARE_METER,
	   1.-hex_scat_x_sec[band][size]/hex_ext_x_sec[band][size],
	   1.-sphere_scat_x_sec[band][size]/sphere_ext_x_sec[band][size],
	   hex_asym_param[band][size],
	   sphere_asym_param[band][size]);
      } /* end loop over sizes */
    } /* end loop over bands */
  } /* end debug */
  
  /* create the matrices of scaling factors to go from Mie parameters to
     hexagonal parameters by scaling */ 
  for(band=1;band<=NUM_HEX_BANDS;band++){
    for(size=1;size<=NUM_HEX_SIZES;size++){
      sphere2hex_scat_factor[band][size]=
	hex_scat_x_sec[band][size]/sphere_scat_x_sec[band][size];
      sphere2hex_ext_factor[band][size]=
	hex_ext_x_sec[band][size]/sphere_ext_x_sec[band][size];
      sphere2hex_abs_factor[band][size]=
	(sphere_abs_x_sec[band][size] > 0.) ? 
	  hex_abs_x_sec[band][size]/sphere_abs_x_sec[band][size] : 1.;
      sphere2hex_asym_factor[band][size]=
	hex_asym_param[band][size]/sphere_asym_param[band][size];
      /* Make sure the scaling factors are physically realistic */ 
      sphere2hex_ext_factor[band][size]=
	(sphere2hex_ext_factor[band][size] < 1.) ?
	  sphere2hex_ext_factor[band][size] : 1.;
      sphere2hex_abs_factor[band][size]=
	(sphere2hex_abs_factor[band][size] < 1.) ?
	  sphere2hex_abs_factor[band][size] : 1.;
    } /* end loop over sizes */
  } /* end loop over bands */
  
  if(debug == 76){
    (void)fprintf(stdout,"Scaling matrices: sphere2hex\n");
    for(band=1;band<=NUM_HEX_BANDS;band++){
      (void)fprintf
	(stdout,"%10s%10s%10s%10s%10s%10s%10s%10s\n",
	 "wavebin","wavelength","sizebin","hex length","s2h ext",
	 "s2h scat","s2h abs","s2h asym");
      (void)fprintf
	(stdout,"%10s%10s%10s%10s%10s%10s%10s\n",
	 "index","microns","index","microns","","","");
      for(size=1;size<=NUM_HEX_SIZES;size++){
	(void)fprintf
	  (stdout,"%10i%10.3f%10i%10.3f%10.6f%10.6f%10.6f%10.6f\n",
	   band,
	   hex_data_wavelength[band]*MICRONS_PER_METER,
	   size,
	   hex_data_length[size]*MICRONS_PER_METER,
	   sphere2hex_ext_factor[band][size],
	   sphere2hex_scat_factor[band][size],
	   sphere2hex_abs_factor[band][size],
	   sphere2hex_asym_factor[band][size]);
      } /* end loop over sizes */
    } /* end loop over bands */
  } /* end debug */
  
  /* multiply the Mie parameters by the scaling parameters to get the
     hex parameters */ 
  hex_size=1;
  for(band=1;band<=num_band;band++){
    for(size=1;size<=num_size;size++){
      if(hex_size < NUM_HEX_SIZES){
	while(hex_column_length[size] > hex_data_length[hex_size]*1.6){
	  hex_size++;
	} /* end while */ 
      } /* end if */
      float_foo=extinction_x_sec[band][size]*
	sphere2hex_ext_factor[CCM2_band2hex_band[band]][hex_size];
      float_foo2=scattering_x_sec[band][size]*
	sphere2hex_scat_factor[CCM2_band2hex_band[band]][hex_size];
      float_foo3=absorption_x_sec[band][size]*
	sphere2hex_abs_factor[CCM2_band2hex_band[band]][hex_size];
      if(float_foo-float_foo2 < 0. && float_foo3 > 0.){
	absorption_x_sec[band][size]*=
	  sphere2hex_abs_factor[CCM2_band2hex_band[band]][hex_size];
	scattering_x_sec[band][size]*=
	  sphere2hex_scat_factor[CCM2_band2hex_band[band]][hex_size];
	extinction_x_sec[band][size]=
	  scattering_x_sec[band][size]+absorption_x_sec[band][size];
      }else{
	extinction_x_sec[band][size]*=
	  sphere2hex_ext_factor[CCM2_band2hex_band[band]][hex_size];
	scattering_x_sec[band][size]*=
	  sphere2hex_scat_factor[CCM2_band2hex_band[band]][hex_size];
	absorption_x_sec[band][size]=
	  extinction_x_sec[band][size]-scattering_x_sec[band][size];
      } /* end else */
      asymmetry_param[band][size]*=
	sphere2hex_asym_factor[CCM2_band2hex_band[band]][hex_size];
    } /* end loop over sizes */
  } /* end loop over bands */
  
  /* Ensure the scattering is never greater than the extinction and 
     compute the absorption cross-section and make sure the asymmetry
     parameter is always between plus and minus one. */ 
  for(band=1;band<=num_band;band++){
    for(size=1;size<=num_size;size++){
      extinction_x_sec[band][size]=
	(extinction_x_sec[band][size] > scattering_x_sec[band][size]) ?
	  extinction_x_sec[band][size] : scattering_x_sec[band][size];
      absorption_x_sec[band][size]=
	extinction_x_sec[band][size]-scattering_x_sec[band][size];
      asymmetry_param[band][size]=
	(asymmetry_param[band][size] <= 1.) ? asymmetry_param[band][size] : 1.;
      asymmetry_param[band][size]=
	(asymmetry_param[band][size] >= -1.) ? asymmetry_param[band][size] : -1.;
    } /* end loop over sizes */
  } /* end loop over bands */
  
} /* end Liou_fudge() */ 

void Liou_interp
  (float *wavelength,int num_band,
   float *crystal_length, int num_size,
   float *equiv_rad_squared, 
   float **extinction_x_sec, /* returned interpolated extinction */ 
   float **scattering_x_sec, /* returned interpolated scattering */ 
   float **absorption_x_sec, /* returned interpolated absorption */ 
   float **asymmetry_param) /* returned interpolated asymmetry */ 
{
  void bilinear_interp
    (float *,int,
     float *,int,
     float *,int,
     float *,int,
     float [NUM_HEX_BANDS+2][NUM_HEX_SIZES+2],
     float [NUM_HEX_BANDS+2][NUM_HEX_SIZES+2],
     float [NUM_HEX_BANDS+2][NUM_HEX_SIZES+2],
     float **,
     float **,
     float **);
  
  float basal_plane_area;
  float basal_plane_radius;
  float float_foo;

  float hex_ext_x_sec[NUM_HEX_BANDS+2][NUM_HEX_SIZES+2];
  float hex_scat_x_sec[NUM_HEX_BANDS+2][NUM_HEX_SIZES+2];

  float hex_aspect_ratio[NUM_HEX_SIZES+2];
  float hex_crystal_mass[NUM_HEX_SIZES+2];
  float hex_diam[NUM_HEX_SIZES+2];
  float hex_equiv_radius[NUM_HEX_SIZES+2];
  float hex_equiv_rad_squared[NUM_HEX_SIZES+2];
  float hex_length[NUM_HEX_SIZES+2];
  float hex_surface_area[NUM_HEX_SIZES+2];
  float hex_volume[NUM_HEX_SIZES+2];
  float hex_wavelength[NUM_HEX_BANDS+2];
  float hex_x_sec_area[NUM_HEX_SIZES+2];

  int band;
  int int_foo;
  int layer;
  int size;

  float hex_wavelength_microns[NUM_HEX_BANDS+2]={
    /* these are the band centers in microns */ 
    0., /* element [0] */ 
    .55,
    1.0,
    1.6,
    2.2,
    3.0,
    0. /* element [NUM_HEX_BANDS+1] */ 
    };
  float hex_length_microns[NUM_HEX_SIZES+2]={
    0., /* element [0] */ 
    20.,
    50.,
    120.,
    300.,
    750.,
    0. /* element [NUM_HEX_SIZES+1] */ 
    };
  float hex_diam_microns[NUM_HEX_SIZES+2]={
    0., /* element [0] */ 
    20.,
    40.,
    60.,
    100.,
    160.,
    0. /* element [NUM_HEX_SIZES+1] */ 
    };
  float hex_ext_x_sec_CGS[NUM_HEX_BANDS+2][NUM_HEX_SIZES+2]={
    /* extinction cross section in units of cm^2 */ 
    0.,0.,0.,0.,0.,0.,0., /* elements [0][*] */ 
    0.,8.5981e-6,4.0392e-5,1.3138e-4,5.1495e-4,1.9663e-3,0.,
    0.,8.5981e-6,4.0392e-5,1.3138e-4,5.1495e-4,1.9663e-3,0.,
    0.,8.5981e-6,4.0392e-5,1.3138e-4,5.1495e-4,1.9663e-3,0.,
    0.,8.5981e-6,4.0392e-5,1.3138e-4,5.1495e-4,1.9663e-3,0.,
    0.,8.5981e-6,4.0392e-5,1.3138e-4,5.1495e-4,1.9663e-3,0.,
    0.,0.,0.,0.,0.,0.,0. /* elements [NUM_HEX_BANDS+1][*] */ 
    };
  float hex_ss_albedo[NUM_HEX_BANDS+2][NUM_HEX_SIZES+2]={
    0.,0.,0.,0.,0.,0.,0., /* elements [0][*] */ 
    0.,0.999999278,0.999998462,0.999997444,0.999995471,0.999992424,0.,
    0.,0.99975683 ,0.999482194,0.999140117,0.998478633,0.997460032,0.,
    0.,0.983894634,0.966534453,0.946369974,0.911990537,0.870501324,0.,
    0.,0.959152259,0.918214765,0.875399152,0.812903237,0.749717993,0.,
    0.,0.53411942 ,0.53138714 ,0.53089312 ,0.53070594 ,0.53064994 ,0.,
    0.,0.,0.,0.,0.,0.,0. /* elements [NUM_HEX_BANDS+1][*] */ 
    };
  float hex_asym_param[NUM_HEX_BANDS+2][NUM_HEX_SIZES+2]={
    0.,0.,0.,0.,0.,0.,0., /* elements [0][*] */ 
    0.,0.770403,0.777955,0.815471,0.842867,0.859161,0.,
    0.,0.776228,0.784654,0.821623,0.849228,0.865001,0.,
    0.,0.789597,0.805137,0.846257,0.880174,0.903669,0.,
    0.,0.818510,0.841624,0.882939,0.916478,0.937994,0.,
    0.,0.956986,0.967520,0.971743,0.973668,0.974627,0.,
    0.,0.,0.,0.,0.,0.,0. /* elements [NUM_HEX_BANDS+1][*] */ 
    };

  /* Convert all units to MKS */ 
  for(band=1;band<=NUM_HEX_BANDS;band++){
    hex_wavelength[band]=hex_wavelength_microns[band]*METERS_PER_MICRON;
    for(size=1;size<=NUM_HEX_SIZES;size++){
      hex_ext_x_sec[band][size]=hex_ext_x_sec_CGS[band][size]*
	METERS_PER_CENTIMETER*METERS_PER_CENTIMETER;
      hex_scat_x_sec[band][size]=
	hex_ss_albedo[band][size]*hex_ext_x_sec[band][size];
    } /* end loop over sizes */
  } /* end loop over bands */
  for(size=1;size<=NUM_HEX_SIZES;size++){
    hex_length[size]=hex_length_microns[size]*METERS_PER_MICRON;
    hex_diam[size]=hex_diam_microns[size]*METERS_PER_MICRON;
  } /* end loop over sizes */

  /* Figure out what the geometric cross-sectional area, volume, and mass
     of each of the hexagonal columns is */ 
  for(size=1;size<=NUM_HEX_SIZES;size++){
    basal_plane_radius=hex_diam[size]/2.;
    float_foo=3.*sqrt(3.)/2.;
    basal_plane_area=float_foo*basal_plane_radius*basal_plane_radius;
    hex_surface_area[size]=2.*basal_plane_area+
      6.*basal_plane_radius*hex_length[size];
    hex_equiv_rad_squared[size]=.25*hex_surface_area[size]/M_PI;
    hex_equiv_radius[size]=sqrt(hex_equiv_rad_squared[size]);

    /* note that the following is ebert and curry's approximation, 
     takano and liou do this correctly */ 
    hex_x_sec_area[size]=hex_surface_area[size]/4.;
    hex_volume[size]=basal_plane_area*hex_length[size];
    hex_crystal_mass[size]=density_of_ice*hex_volume[size];
    hex_aspect_ratio[size]=hex_length[size]/hex_diam[size];
  } /* end loop over sizes */

  if(debug == 75){
    (void)fprintf(stdout,"sizebin\tLength\tDiam\tAspect\tr_e\tr_e^2\tSurf A\tX-sec A\tMass\n");
    (void)fprintf(stdout,"\tmicrons\tmicrons\t\tmicrons\tmu^2\tmu^2\tmu^2\tng\n");
    for(size=1;size<=NUM_HEX_SIZES;size++){
      (void)fprintf(stdout,"%i\t%.1f\t%.1f\t%.1f\t%.1f\t%g\t%.0f\t%.0f\t%.0f\n",
		    size,
		    hex_length[size]*MICRONS_PER_METER,
		    hex_diam[size]*MICRONS_PER_METER,
		    hex_aspect_ratio[size],
		    hex_equiv_radius[size]*MICRONS_PER_METER,
		    hex_equiv_rad_squared[size]*SQUARE_MICRONS_PER_SQUARE_METER,
		    hex_surface_area[size]*SQUARE_MICRONS_PER_SQUARE_METER,
		    hex_x_sec_area[size]*SQUARE_MICRONS_PER_SQUARE_METER,
		    hex_crystal_mass[size]*NANOGRAMS_PER_KILOGRAM);
    } /* end loop over sizes */
  } /* end debug */

  /* Interpolate Liou's hexagonal column properties to my size/wavelength grid.
     Note that interpolation based on the surface area geometry seems to be much 
     more accurate than interpolation based on the length geometry.
     Remember the Liou geometry must match the interpolated geometry,
     i.e., both crystal lengths or areas...  */
  bilinear_interp
    (wavelength,num_band, /* wavelengths at interp points */ 
/*     crystal_length,num_size, */
     equiv_rad_squared,num_size, 
     hex_wavelength,NUM_HEX_BANDS, /* Liou wavelength data */ 
/*     hex_length,NUM_HEX_SIZES, */
     hex_equiv_rad_squared,NUM_HEX_SIZES, 
     hex_ext_x_sec, /* Liou extinction data */ 
     hex_scat_x_sec, /* Liou scattering data */ 
     hex_asym_param, /* Liou asymmetry parameter data */ 
     extinction_x_sec, /* interpolated extinction */ 
     scattering_x_sec, /* interpolated scattering */ 
     asymmetry_param); /* interpolated asymmetry */ 

  /* Ensure the scattering is never greater than the extinction and 
     compute the absorption cross-section */ 
  for(band=1;band<=num_band;band++){
    for(size=1;size<=num_size;size++){
      extinction_x_sec[band][size]=
	(extinction_x_sec[band][size] > scattering_x_sec[band][size]) ?
	  extinction_x_sec[band][size] : scattering_x_sec[band][size];
      absorption_x_sec[band][size]=
	extinction_x_sec[band][size]-scattering_x_sec[band][size];
    } /* end loop over sizes */
  } /* end loop over bands */

  /* ensure that the smallest absorption x-section grows at least quadratically 
     with size */ 
  for(band=1;band<=num_band;band++){
    for(size=10;size>=1;size--){
      if(absorption_x_sec[band][size] >= absorption_x_sec[band][size+1]){
	absorption_x_sec[band][size]=absorption_x_sec[band][size+1]*
	  equiv_rad_squared[size+1]/equiv_rad_squared[size+2];
	extinction_x_sec[band][size]=extinction_x_sec[band][size+1]*
	  equiv_rad_squared[size+1]/equiv_rad_squared[size+2];
	scattering_x_sec[band][size]=
	  extinction_x_sec[band][size]-absorption_x_sec[band][size];
      } /* end if */
    } /* end loop over sizes */
  } /* end loop over bands */

} /* end Liou_interp() */ 

void bilinear_interp
  (float *x, /* i.e., wavelength at interpolation = x */ 
   int num_x_idx, /* size of x array */ 
   float *y, /* i.e., size at interpolation = y */
   int num_y_idx, /* size of y array */ 
   float *abscissa_grid, /* Liou abscissa vector */ 
   int num_abscissa, /* length of abscissa vector */ 
   float *ordinate_grid, /* Liou ordinate vector */ 
   int num_ordinate, /* length of ordinate vector */ 
   /* Liou extinction data */ 
   float hex_ext_x_sec[NUM_HEX_BANDS+2][NUM_HEX_SIZES+2], 
   /* Liou scattering data */ 
   float hex_scat_x_sec[NUM_HEX_BANDS+2][NUM_HEX_SIZES+2], 
   /* Liou asymmetry parameter data */ 
   float hex_asym_param[NUM_HEX_BANDS+2][NUM_HEX_SIZES+2], 
   float **extinction_x_sec, /* interpolated extinction */ 
   float **scattering_x_sec, /* interpolated scattering */ 
   float **asymmetry_param) /* interpolated asymmetry */ 
{
  /* Reference: Numerical Recipes in C, Ch. 3.6 p. 104 */ 
  
  float max_VEC(float *,int);
  float min_VEC(float *,int);

  int find_smaller_neighbor(float,float *,int);

  float *x_fraction;
  float *y_fraction;

  int *left_neighbor;
  int *right_neighbor;
  int *lower_neighbor;
  int *upper_neighbor;

  float float_foo;

  int abscissa;
  int ordinate;
  int x_idx;
  int y_idx;

  if(
     ((left_neighbor=(int *)malloc((num_x_idx+2)*sizeof(int))) == NULL ) ||
     ((right_neighbor=(int *)malloc((num_x_idx+2)*sizeof(int))) == NULL ) ||
     ((lower_neighbor=(int *)malloc((num_y_idx+2)*sizeof(int))) == NULL ) ||
     ((upper_neighbor=(int *)malloc((num_y_idx+2)*sizeof(int))) == NULL ) ||
     False ){
    (void)fprintf(stdout,"Unable to allocate array in interp\n");
    exit(1);
  } /* end if */

  if(
     ((x_fraction=(float *)malloc((num_x_idx+2)*sizeof(float))) == NULL ) ||
     ((y_fraction=(float *)malloc((num_y_idx+2)*sizeof(float))) == NULL ) ||
     False ){
    (void)fprintf(stdout,"Unable to allocate array in interp\n");
    exit(1);
  } /* end if */

  /* Before performing the bilinear interpolation it is necessary to put boundary
     values of Liou's data that will "circumscribe" the any possible interpolation
     point, so find the minimum and maximum interpolation x and y, and place a
     "real" data point at that spot in the data grid */ 
  abscissa_grid[0]=min_VEC(x+1,num_x_idx);
  abscissa_grid[num_abscissa+1]=max_VEC(x+1,num_x_idx);
  ordinate_grid[0]=min_VEC(y+1,num_y_idx);
  ordinate_grid[num_ordinate+1]=max_VEC(y+1,num_y_idx);

  /* Must also fill in the corners of the input data grids so the bilinear
   interpolation will work transparently at the new boundaries. 
   Use the appropriate first order interpolation approximation that
   z = z1 + (x-x1)*((z2 -z1)/(x2-x1)) OR z = z2 + (x-x2)*((z2 -z1)/(x2-x1)) 
   z = z1 + (y-y1)*((z2 -z1)/(y2-y1)) OR z = z2 + (y-y2)*((z2 -z1)/(y2-y1)) 
   to fill in the needed boundaries, this applies to x-sections only.

   The asymmetry parameter cannot be expected to scale accurately with 
   crude resolution wavelength/size resolution, so just set boundaries equal
   to nearest neighbors. 

   Also, extrapolating the extinction and absorption to the smallest sizes
   leads to unphysical results in the form of scattering being larger than
   extinction, so just extrapolate the smallest scattering and extinction 
   by the same amount! */ 

  /* first go through an abscissa and fill in the side edges of the window pane,
     the ordinate loop will take care of the top and bottom edges including all
     four corners */ 
  for(abscissa=1;abscissa<=num_abscissa;abscissa++){
    hex_ext_x_sec[abscissa][0]=
      hex_ext_x_sec[abscissa][1]+
	(ordinate_grid[0]-ordinate_grid[1])*
	  (hex_ext_x_sec[abscissa][2]-hex_ext_x_sec[abscissa][1])/
	    (ordinate_grid[2]-ordinate_grid[1]);
    hex_ext_x_sec[abscissa][num_ordinate+1]=
      hex_ext_x_sec[abscissa][num_ordinate]+
	(ordinate_grid[num_ordinate+1]-ordinate_grid[num_ordinate])*
	  (hex_ext_x_sec[abscissa][num_ordinate]-hex_ext_x_sec[abscissa][num_ordinate-1])/
	    (ordinate_grid[num_ordinate]-ordinate_grid[num_ordinate-1]);

    hex_scat_x_sec[abscissa][0]=
      hex_scat_x_sec[abscissa][1]+
	(ordinate_grid[0]-ordinate_grid[1])*
	  (hex_scat_x_sec[abscissa][2]-hex_scat_x_sec[abscissa][1])/
	    (ordinate_grid[2]-ordinate_grid[1]);
    hex_scat_x_sec[abscissa][num_ordinate+1]=
      hex_scat_x_sec[abscissa][num_ordinate]+
	(ordinate_grid[num_ordinate+1]-ordinate_grid[num_ordinate])*
	  (hex_scat_x_sec[abscissa][num_ordinate]-hex_scat_x_sec[abscissa][num_ordinate-1])/
	    (ordinate_grid[num_ordinate]-ordinate_grid[num_ordinate-1]);

    hex_asym_param[abscissa][0]=hex_asym_param[abscissa][1];
    hex_asym_param[abscissa][num_ordinate+1]=
      hex_asym_param[abscissa][num_ordinate];

    /* This kludge is necessary to preserve absorption in the smallest
       particles in the cirrus grid which might otherwise get extrapolated
       out due to Liou's coarse data mesh.  It also fixes the non-physical
       absorption which would occur in the UV. */ 
    float_foo=hex_scat_x_sec[abscissa][1]-hex_scat_x_sec[abscissa][0];
    hex_ext_x_sec[abscissa][0]=hex_ext_x_sec[abscissa][1]-float_foo;

    /* Ensure the extinction and scattering x-sections are never negative */ 
    hex_ext_x_sec[abscissa][0]=(hex_ext_x_sec[abscissa][0] > 0.) ?
      hex_ext_x_sec[abscissa][0] : 
	hex_ext_x_sec[abscissa][1]*(ordinate_grid[0]/ordinate_grid[1]);
    hex_scat_x_sec[abscissa][0]=(hex_scat_x_sec[abscissa][0] > 0.) ?
      hex_scat_x_sec[abscissa][0] : 
	hex_scat_x_sec[abscissa][1]*(ordinate_grid[0]/ordinate_grid[1]);
    
  } /* end loop over abscissas */

  /* the ordinate loop will be used to set the four corners
     of the array as well, so index from 0 to num_ordinate+1 */ 
  for(ordinate=0;ordinate<=num_ordinate+1;ordinate++){
    hex_ext_x_sec[0][ordinate]=
      hex_ext_x_sec[1][ordinate]+
	(abscissa_grid[0]-abscissa_grid[1])*
	  (hex_ext_x_sec[2][ordinate]-hex_ext_x_sec[1][ordinate])/
	    (abscissa_grid[2]-abscissa_grid[1]);
    hex_ext_x_sec[num_abscissa+1][ordinate]=
      hex_ext_x_sec[num_abscissa][ordinate]+
	(abscissa_grid[num_abscissa+1]-abscissa_grid[num_abscissa])*
	  (hex_ext_x_sec[num_abscissa][ordinate]-hex_ext_x_sec[num_abscissa-1][ordinate])/
	    (abscissa_grid[num_abscissa]-abscissa_grid[num_abscissa-1]);

    hex_scat_x_sec[0][ordinate]=
      hex_scat_x_sec[1][ordinate]+
	(abscissa_grid[0]-abscissa_grid[1])*
	  (hex_scat_x_sec[2][ordinate]-hex_scat_x_sec[1][ordinate])/
	    (abscissa_grid[2]-abscissa_grid[1]);
    hex_scat_x_sec[num_abscissa+1][ordinate]=
      hex_scat_x_sec[num_abscissa][ordinate]+
	(abscissa_grid[num_abscissa+1]-abscissa_grid[num_abscissa])*
	  (hex_scat_x_sec[num_abscissa][ordinate]-hex_scat_x_sec[num_abscissa-1][ordinate])/
	    (abscissa_grid[num_abscissa]-abscissa_grid[num_abscissa-1]);

    hex_asym_param[0][ordinate]=hex_asym_param[1][ordinate];
    hex_asym_param[num_abscissa+1][ordinate]=
      hex_asym_param[num_abscissa][ordinate];

    /* Ensure the extinction and scattering x-sections are never negative */ 
    hex_ext_x_sec[0][ordinate]=(hex_ext_x_sec[0][ordinate] > 0.) ?
      hex_ext_x_sec[0][ordinate] : 
	hex_ext_x_sec[1][ordinate]*(abscissa_grid[0]/abscissa_grid[1]);
    hex_scat_x_sec[0][ordinate]=(hex_scat_x_sec[0][ordinate] > 0.) ?
      hex_scat_x_sec[0][ordinate] : 
	hex_scat_x_sec[1][ordinate]*(abscissa_grid[0]/abscissa_grid[1]);

    /* because scattering xsections are falling rapidly at 4 microns as
       absorption picks up, this crude interpolation will lead to negative
       (bogus) results, and make sure the kludge has scat. xsection
       decreasing with increasing wavelength */
    hex_scat_x_sec[num_abscissa+1][ordinate]=
      (hex_scat_x_sec[num_abscissa+1][ordinate] > 0.) ?
	hex_scat_x_sec[num_abscissa+1][ordinate] : 
	hex_scat_x_sec[num_abscissa][ordinate]*
	  (abscissa_grid[num_abscissa]/abscissa_grid[num_abscissa+1]);

  } /* end loop over ordinates */

  /* Search to find the grid square in which the point (x1,x2) falls.
     The first and last interpolation points might well fall exactly
     on the boundaries of the data grid, so be sure to use 
     "less than or equals" and "greater than or equals" in the
     comparisons. */
  for(x_idx=1;x_idx<=num_x_idx;x_idx++){
    left_neighbor[x_idx]=find_smaller_neighbor
      (x[x_idx],abscissa_grid,num_abscissa);
    if(left_neighbor[x_idx] == num_abscissa+1){
      left_neighbor[x_idx] = num_abscissa;
    } /* end if */
    right_neighbor[x_idx]=left_neighbor[x_idx]+1;
    x_fraction[x_idx]=
      (x[x_idx]-abscissa_grid[left_neighbor[x_idx]])/
      (abscissa_grid[right_neighbor[x_idx]]-abscissa_grid[left_neighbor[x_idx]]);
  } /* end loop over x_idx */
  
  for(y_idx=1;y_idx<=num_y_idx;y_idx++){
    lower_neighbor[y_idx]=find_smaller_neighbor
      (y[y_idx],ordinate_grid,num_ordinate);
    if(lower_neighbor[y_idx] == num_ordinate+1){
      lower_neighbor[y_idx] = num_ordinate;
    } /* end if */
    upper_neighbor[y_idx]=lower_neighbor[y_idx]+1;
    y_fraction[y_idx]=
      (y[y_idx]-ordinate_grid[lower_neighbor[y_idx]])/
      (ordinate_grid[upper_neighbor[y_idx]]-ordinate_grid[lower_neighbor[y_idx]]);
  } /* end loop over y_idx */

  if(debug == 75){
    (void)fprintf(stdout,"x-absc\tabs-grid\n");
    for(abscissa=0;abscissa<=num_abscissa+1;abscissa++){
      (void)fprintf(stdout,"%i\t%g\n",abscissa,abscissa_grid[abscissa]);
    } /* end loop over abscissas */

    (void)fprintf(stdout,"y-ord\tord-grid\n");
    for(ordinate=0;ordinate<=num_ordinate+1;ordinate++){
      (void)fprintf(stdout,"%i\t%g\n",ordinate,ordinate_grid[ordinate]);
    } /* end loop over ordinates */

    for(abscissa=0;abscissa<=num_abscissa+1;abscissa++){
      for(ordinate=0;ordinate<=num_ordinate+1;ordinate++){
	(void)fprintf(stdout,"hex_ext_x_sec[%i][%i] = %g cm^2, hex_scat_x_sec[%i][%i] = %g cm^2\n",
		      abscissa,
		      ordinate,
		      hex_ext_x_sec[abscissa][ordinate]*SQUARE_CMS_PER_SQUARE_METER,
		      abscissa,
		      ordinate,
		      hex_scat_x_sec[abscissa][ordinate]*SQUARE_CMS_PER_SQUARE_METER);
      } /* end loop over ordinates */
    } /* end loop over abscissas */

    (void)fprintf(stdout,"x_indx\tx[x_i]\tl_idx\tabs[li]\tr_idx\tabs[ri]\txfrac\n");
    for(x_idx=1;x_idx<=num_x_idx;x_idx++){
      (void)fprintf(stdout,"%i\t%g\t%i\t%g\t%i\t%g\t%g\n",
		    x_idx,
		    x[x_idx],
		    left_neighbor[x_idx],
		    abscissa_grid[left_neighbor[x_idx]],
		    right_neighbor[x_idx],
		    abscissa_grid[right_neighbor[x_idx]],
		    x_fraction[x_idx]);
    } /* end loop over x_idx */

    (void)fprintf(stdout,"y_indx\ty[y_i]\tl_idx\tord[li]\tu_idx\tord[ui]\tyfrac\n");
    for(y_idx=1;y_idx<=num_y_idx;y_idx++){
      (void)fprintf(stdout,"%i\t%g\t%i\t%g\t%i\t%g\t%g\n",
		    y_idx,
		    y[y_idx],
		    lower_neighbor[y_idx],
		    ordinate_grid[lower_neighbor[y_idx]],
		    upper_neighbor[y_idx],
		    ordinate_grid[upper_neighbor[y_idx]],
		    y_fraction[y_idx]);
    } /* end loop over y_idx */

  } /* end debug */
  
  /* Now that the interpolation ratios have been determined, fill
   in the interpolated arrays for each field desired */ 
  for(x_idx=1;x_idx<=num_x_idx;x_idx++){
    for(y_idx=1;y_idx<=num_y_idx;y_idx++){
      extinction_x_sec[x_idx][y_idx]=
	(1.-x_fraction[x_idx])*(1.-y_fraction[y_idx])*
	  hex_ext_x_sec[left_neighbor[x_idx]][lower_neighbor[y_idx]];
      extinction_x_sec[x_idx][y_idx]+=
	x_fraction[x_idx]*(1.-y_fraction[y_idx])*
	  hex_ext_x_sec[right_neighbor[x_idx]][lower_neighbor[y_idx]];
      extinction_x_sec[x_idx][y_idx]+=
	x_fraction[x_idx]*y_fraction[y_idx]*
	  hex_ext_x_sec[right_neighbor[x_idx]][upper_neighbor[y_idx]];
      extinction_x_sec[x_idx][y_idx]+=
	(1.-x_fraction[x_idx])*y_fraction[y_idx]*
	  hex_ext_x_sec[left_neighbor[x_idx]][upper_neighbor[y_idx]];
      extinction_x_sec[x_idx][y_idx]=
	(extinction_x_sec[x_idx][y_idx] > 0.) ?
	  extinction_x_sec[x_idx][y_idx] : 0.;

      scattering_x_sec[x_idx][y_idx]=
	(1.-x_fraction[x_idx])*(1.-y_fraction[y_idx])*
	  hex_scat_x_sec[left_neighbor[x_idx]][lower_neighbor[y_idx]];
      scattering_x_sec[x_idx][y_idx]+=
	x_fraction[x_idx]*(1.-y_fraction[y_idx])*
	  hex_scat_x_sec[right_neighbor[x_idx]][lower_neighbor[y_idx]];
      scattering_x_sec[x_idx][y_idx]+=
	x_fraction[x_idx]*y_fraction[y_idx]*
	  hex_scat_x_sec[right_neighbor[x_idx]][upper_neighbor[y_idx]];
      scattering_x_sec[x_idx][y_idx]+=
	(1.-x_fraction[x_idx])*y_fraction[y_idx]*
	  hex_scat_x_sec[left_neighbor[x_idx]][upper_neighbor[y_idx]];
      scattering_x_sec[x_idx][y_idx]=
	(scattering_x_sec[x_idx][y_idx] > 0.) ?
	  scattering_x_sec[x_idx][y_idx] : 0.;

      asymmetry_param[x_idx][y_idx]=
	(1.-x_fraction[x_idx])*(1.-y_fraction[y_idx])*
	  hex_asym_param[left_neighbor[x_idx]][lower_neighbor[y_idx]];
      asymmetry_param[x_idx][y_idx]+=
	x_fraction[x_idx]*(1.-y_fraction[y_idx])*
	  hex_asym_param[right_neighbor[x_idx]][lower_neighbor[y_idx]];
      asymmetry_param[x_idx][y_idx]+=
	x_fraction[x_idx]*y_fraction[y_idx]*
	  hex_asym_param[right_neighbor[x_idx]][upper_neighbor[y_idx]];
      asymmetry_param[x_idx][y_idx]+=
	(1.-x_fraction[x_idx])*y_fraction[y_idx]*
	  hex_asym_param[left_neighbor[x_idx]][upper_neighbor[y_idx]];
      asymmetry_param[x_idx][y_idx]=
	(asymmetry_param[x_idx][y_idx] > -1.) ?
	  asymmetry_param[x_idx][y_idx] : -1.;
      asymmetry_param[x_idx][y_idx]=
	(asymmetry_param[x_idx][y_idx] < 1.) ?
	  asymmetry_param[x_idx][y_idx] : 1.;

    } /* end loop over x_idx */
  } /* end loop over y_idx */

  /* Get rid of the excess memory */ 
  free(left_neighbor);
  free(right_neighbor);
  free(lower_neighbor);
  free(upper_neighbor);
  free(x_fraction);
  free(y_fraction);
} /* end bilinear_interp() */ 

int find_smaller_neighbor
  (float abscissa,
   float *abscissa_grid,
   int num_abscissa)
{
  /* returns the index i into the abscissa_grid array which points
     to the greatest abscissa_grid less than or equal the to abscissa.
     i.e. abscissa_grid[i] = GLB(abscissa). */

  /* NB: assumes abscissa_grid is indexed [0..num_abscissa+1] */ 
  
  int index;

  index=0;
  while((abscissa_grid[index] <= abscissa) && (index <= num_abscissa+1)){
    index++;
  } /* end while */ 
  return --index;
} /* end find_smaller_neighbor() */ 

