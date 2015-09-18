static char rcs_Id[] = "$Id$\n";
static char rcs_Revision[] = "$Revision$";

/* $Author: zender $
 * $Date$ 
 * $Locker:  
 * $RCSfile: main.c,v $
 * $Source: /home/zender/cvs/cld/main.c,v $
 * $Id$
 * $State: Exp $
 */

/* Purpose: Thesis project modelling of cirrus cloud. */ 

/* Example usage:
   the 6-12km altitude scale
   cld -R -n 0 -x 35 -l 80 -G 6000 -g 6000 -D 0 > & ! foo ; less foo
   RAMASWAMY initialization:
   cld -F -U -Z 3 -N 3.3e4 -n 1 -x 3 -l 80 -m 2.00e-4 -b 5 -D 64 > & ! foo2 ; less foo2
   ZHANG initialization:
   cld -F -U -Z 4 -n 1 -x 40 -l 80 -w .05 -D 64 > & ! foo2 ; less foo2
   TROPICS initialization:
   cld -M 1 -H 15000 -h 1000 -g 4000 -G 13000 -X 

########################################################################
# Perform the new standard cloud integration
########################################################################
cld -D 55 -E -e err -G 10000 -g 8000 -H 15000 -h 2000 -k 3. -l 103 -M 1 -m 3.e-6 -N 1.e6 -n 1200 -o nc -p 20 -r 5 -S 1.2 -s .4 -X -x 50

########################################################################
# Crashes flamenco
########################################################################
cld -D 55 -e err -G 10000 -g 8000 -H 15000 -h 2000 -k 3. -l 127 -M 1 -m 3.e-6 -N 1.e6 -n 1200 -o nc -p 20 -r 5 -S 1.2 -s .4 -X -x 50

########################################################################
# Modifications to standard: truncated distribution (-x 40 -m 10.e-6) or
# (-x 35 -m 20.e-6)
########################################################################
cld -D 55 -E -e err -G 10000 -g 8000 -H 15000 -h 2000 -k 3. -l 103 -M 1 -m 20.e-6 -N 1.e6 -n 1200 -o nc -p 20 -r 5 -S 1.2 -s .4 -X -x 35

########################################################################
# Modifications to standard: radiative spheres (no -X)
########################################################################
#cld -D 55 -E -e err -G 10000 -g 8000 -H 15000 -h 2000 -k 3. -l 103 -M 1 -m 3.e-6 -N 1.e6 -n 1200 -o nc -p 20 -r 5 -S 1.2 -s .4 -x 50

########################################################################
# Modifications to standard: radiative spheres, truncated 
########################################################################
#cld -D 55 -E -e err -G 10000 -g 8000 -H 15000 -h 2000 -k 3. -l 103 -M 1 -m 20.e-6 -N 1.e6 -n 1200 -o nc -p 20 -r 5 -S 1.2 -s .4 -x 35

########################################################################
# Modifications to standard: 4 hour integration (-k 3. -n 4800)
########################################################################
#cld -D 55 -E -e err -G 10000 -g 8000 -H 15000 -h 2000 -k 3. -l 103 -M 1 -m 3.e-6 -N 1.e6 -n 4800 -o nc -p 20 -r 5 -S 1.2 -s .4 -X -x 50

########################################################################
# Modifications to standard: fast updraft (-w .20)
########################################################################
#cld -D 55 -E -e err -G 10000 -g 8000 -H 15000 -h 2000 -k 3. -l 103 -M 1 -m 3.e-6 -N 1.e6 -n 1200 -o nc -p 20 -r 5 -S 1.2 -s .4 -X -x 50 -w .20
   */

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
 * Revision 4.1  1993/02/27  00:05:18  zender
 * incorporated ascii writing automatically on the CRAY in main,
 * merged ioinput and iooutput into iographs. added albedo vs.
 * emissivity graph. alphabetized output. might work on netCDF now.
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
 *  */

/* standard header files */
#include <stdio.h>              /* stderr, FILE, NULL, etc. */
#include <time.h>               /* machine time */
#include <string.h>             /* strcmp. . . */
#include <math.h>               /* sin cos cos sin 3.14159 */
#include <stdlib.h>             /* atof, atoi, malloc, getopt */ 
#include <unistd.h>             /* all sorts of POSIX stuff */ 
#include <netcdf.h>             /* netCDF routines */
#if ( defined Solaris )
#include "/opt/SUNWspro/SC3.0/include/cc/sunmath.h"            /* IEEE signal handling */
#include <siginfo.h>            /* IEEE stuff, je ne sais pas */ 
#include <ucontext.h>           /* IEEE stuff, je ne sais pas */ 
/* #include <ieeefp.h>*/
#endif
#if ( defined SunOS )
#include <signal.h>            /* IEEE stuff, je ne sais pas */ 
#endif

/* my header files */
#include "defs.h"               /* YES, NO, FILESIZE, and more */

#define MIN_CONC_FOR_INIT 10.         /* lower concentrations will be zeroed */ 
#define RAMASWAMY_TEMPERATURE 228.    /* for validation initializations */ 
#define RAMASWAMY_PRESSURE 25000.     /* for validation initializations */ 
#define LW_DIFFUSIVITY 1.66           /* Should be between 1.65 and 2. */ 
  /* Fewer than 60 sizes results in underestimating the concentration,
     i.e, the integral of the distribution will not converge to
     CCN_tot_conc */
#define NUM_CCN_SIZES 60  

/* Global variables declared here */ 
int debug=0; /* Option D */
int debug_value = 0; /* Option d */

main(int argc,char **argv)
{
  /* NB:  Everything labelled divergence is really a convergence */ 

#if ( defined SunOS )
  void ieee_exception_handler
    (int,
     int,
     struct sigcontext *,
     char *);
#endif
#if ( defined Solaris )
  void ieee_exception_handler
    (int,
     siginfo_t *, 
     ucontext_t *);
#endif

  /* "Prototype" the functions first, then the variables */
  float Aspect_ratio
    (float /* hex_crystal_length */,
     float */* daspect_ratiodL */); 
  float Capacitance
    (float /* length */,
     float /* diameter */);
  float Effective_Radius();
  float Eqm_Vap_Ice(float /* temperature */ );
  float Eqm_Vap_Liquid(float /* temperature */ );
  float Fall_Speed_Complex
    (float /* crystal_length */,
     float /* env_pressure */,
     float /* env_temperature */, 
     float /* env_density */, 
     int /* ice_crystal_habit */);
  float Fall_Speed_Simple();
  float IWC_of_temp();
  float PDF_of_length(float /* length */,int /* initial_distribution_type */,
		      float /* cloud_scale_heights_above_center */);
  float Thermal_Conductivity();
  float Vapor_Diffusivity();
  float conc_left_redistribution();
  float conc_right_redistribution();
  float dot_product_VEC();
  float ice_density_of_length();
  float bullet_length_of_mass(float /* mass */ );
  float bullet_length_to_equiv_rad(float /* length */ );
  float homogeneous_nucleation_DMC94
    (float * /* CCN_conc */,
     float * /* CCN_diameter */,
     float /* temperature */,
     int /* num_CCN_size */, 
     float /* saturation_liquid */);
  float mass_of_bullet_length(float /* length */ );
  float max_MAT();
  float max_abs_MAT();
  float max_abs_MAT_offset();
  float max_VEC();
  float max_abs_VEC();
  float min_MAT();
  float min_VEC();
  float heterogeneous_nucleation_MDC92();
  float size_dist_of_IWC(float /* length */,float /* temperature */);
  float total_VEC(float *,int);
  
  int Cloud_Base_idx();
  int Cloud_Top_idx();
  int find_nearest_mass_bin();
  
  void FORTRAN_callbh();
  void FORTRAN_colmod();
  void FORTRAN_consre();
  void FORTRAN_etbfct();
  void FORTRAN_ngride();
  void FORTRAN_ogride();
  void FORTRAN_veloce();

  void Bohren_Mie();
  void Exit_gracefully(void);
  void Kinetic_Ventilation();
  void Liou_IR_fudge
    (int /* num_LW_band */,
     int /* num_size */,
     float **/* LW_absorption_x_sec */, /* returned scaled absorption */ 
     float **/* LW_scattering_x_sec */, /* returned scaled scattering */ 
     float **/* LW_extinction_x_sec */, /* returned scaled extinction */ 
     float **/* LW_asymmetry_param */); /* returned scaled asymmetry */ 
  void Liou_fudge(float *,int,
		  float *,int,
		  float *,
		  float *,
		  float **,
		  float **,
		  float **,
		  float **);
  void Liou_interp(float *,int,
		  float *,int,
		  float *,
		  float **,
		  float **,
		  float **,
		  float **);
  void Warren_SW_ice_init();
  void Warren_LW_ice_init();
  void apportion_crystal_heating();
  void hex_column_dimensions_of_mass
    (float,
     float *,
     float *,
     float *,
     float *);
  void init_CCN_conc
    (float */* CCN_mass */,
     float */* CCN_diameter */,
     float */* CCN_conc */,
     float /* CCN_tot_conc */,
     float /* CCN_scaling_diameter */,
     int /* num_CCN_size */);
  void init_ccm2_input();
  void init_temp_pres
    (float */* altitude */,
     float */* env_pressure */,
     float */* env_temperature */,
     int /* num_layer */,
     int /* atmospheric_profile_type */);
  void init_IR_params
    (float **/* LW_planck_radiance */, 
     float */* IR_bandwidth */,
     float */* IR_wavelength */,
     float **/* LW_spec_band_weight */,
     float **/* LW_absorption_x_sec */,
     float **/* planck_avg_abs_x_sec */,
     float */* env_temperature */,
     int /* num_LW_band */,
     int /* num_layer */,
     int /* num_size */);
  void integrate_LW_mass_abs_coeff();
  void integrate_SW_optical_properties
    (float **,
     float **,
     float,
     float **,
     float **,
     float **,
     float **,
     float **,
     int,
     int,
     int,
     float **);
  void integrate_LW_abs_opt_depth();
  void max_FLOATVEC_idxed();
  void max_INTVEC_idxed();
  void print_debug(void);
  void print_usage(char *);
  void refraction();
  void string_cmdline(int,char *[],char *,int);
  void transpose_conc();
  void transpose_trans();
  void write_ccm2_output();
  
  Boolean AGGREGATION;
  Boolean ERROR_FLAG;
  Boolean EULER;
  Boolean EVEN_STEP;
  Boolean HETEROGENEOUS_NUCLEATION;
  Boolean HOMOGENEOUS_NUCLEATION;
  Boolean ICE_ADVECT;
  Boolean LAGRANGE;
  Boolean LIOU;
  Boolean MOVIE_NETCDF_OUTPUT;
  Boolean NETCDF_OUTPUT;
  Boolean PAUSE;
  Boolean RADIATION;
  Boolean RAD_FEEDBACK;
  Boolean UPDRAFT_PROFILE_IS_FLAT;
  Boolean VAPOR_ADVECT;
  Boolean VERBOSE;
  Boolean WEED_NEGATIVE_VALUES;
  
  char *char_ptr_foo;
  char *char_ptr_foo2;
  char *option_string;
  char *time_buf_start;

  char *atmospheric_profile_string_array[]={
    "MID_LATITUDE_SUMMER","TROPICS"}; 
  char *initial_distribution_string_array[]={
    "HEYMSFIELD","DOWLING_RADKE","KNOLLENBERG","RAMASWAMY","ZHANG"};
  char *ice_crystal_habit_string_array[]={
  "HEXAGONAL_COLUMNS","HEXAGONAL_BULLETS","SPHERES"};

  char cmdline[CMDLINE_SIZE];
  char err_file[FILESIZE];
  char in_file[FILESIZE];
  char out_file[FILESIZE];
  
  extern char *optarg;
  extern int optind;
  
  float **IWC_snapshot;
  float **LW_absorption_eff;
  float **LW_asymmetry_param;
  float **LW_extinction_eff;
  float **LW_planck_radiance;
  float **LW_scattering_eff;
  float **LW_size_param;
  float **LW_spec_band_weight;
  float **LW_spec_vol_abs_coeff;
  float **SW_spec_flux_div;
  float **absorption_eff;
  float **absorption_x_sec;
  float **extinction_x_sec;
  float **scattering_x_sec;
  float **LW_absorption_x_sec;
  float **LW_extinction_x_sec;
  float **LW_scattering_x_sec;
  float **asymmetry_param;
  float **conc_transpose;
  float **conc_transpose_new;
  float **concentration;
  float **concentration_new;
  float **concentration_swap;
  float **conductivity_correction;
  float **crystal_temperature;
  float **crystal_heat_net;
  float **crystal_heat_SW;
  float **crystal_heat_LW;
  float **daltitude_dtime;
  float **diffusivity_correction;
  float **distribution;
  float **distribution_new;
  float **dmass_dtime;
  float **eff_asymmetry_param;
  float **eff_ss_albedo;
  float **eff_tau_extinction;
  float **eff_tau_scatter;
  float **extinction_eff;
  float **fall_speed;
  float **ice_density;
  float **lagrange_mass;
  float **number_snapshot;
  float **planck_avg_abs_x_sec;
  float **scattering_eff;
  float **size_param;
  
  float * IWC_snapshot_ptr;
  float * LW_absorption_eff_ptr;
  float * LW_asymmetry_param_ptr;
  float * LW_extinction_eff_ptr;
  float * LW_planck_radiance_ptr;
  float * LW_scattering_eff_ptr;
  float * LW_size_param_ptr;
  float * LW_spec_band_weight_ptr;
  float * LW_spec_vol_abs_coeff_ptr;
  float * SW_spec_flux_div_ptr;
  float * absorption_eff_ptr;
  float * absorption_x_sec_ptr;
  float * extinction_x_sec_ptr;
  float * scattering_x_sec_ptr;
  float * LW_absorption_x_sec_ptr;
  float * LW_extinction_x_sec_ptr;
  float * LW_scattering_x_sec_ptr;
  float * asymmetry_param_ptr;
  float * conc_transpose_new_ptr;
  float * conc_transpose_ptr;
  float * concentration_new_ptr;
  float * concentration_ptr;
  float * conductivity_correction_ptr;
  float * crystal_temperature_ptr;
  float * crystal_heat_net_ptr;
  float * crystal_heat_SW_ptr;
  float * crystal_heat_LW_ptr;
  float * daltitude_dtime_ptr;
  float * diffusivity_correction_ptr;
  float * distribution_new_ptr;
  float * distribution_ptr;
  float * dmass_dtime_ptr;
  float * eff_asymmetry_param_ptr;
  float * eff_ss_albedo_ptr;
  float * eff_tau_extinction_ptr;
  float * eff_tau_scatter_ptr;
  float * extinction_eff_ptr;
  float * fall_speed_ptr;
  float * ice_density_ptr;
  float * lagrange_mass_ptr;
  float * number_snapshot_ptr;
  float * planck_avg_abs_x_sec_ptr;
  float * scattering_eff_ptr;
  float * size_param_ptr;

  float *crystal_diameter;
  float *aspect_ratio;
  float *saturation_liquid;
  float *CCN_activated;
  float *CCN_conc;
  float *CCN_mass;
  float *CCN_diameter;
  float *hex_column_diameter;
  float *hex_column_length;
  float *hex_equiv_radius;
  float *hex_equiv_rad_squared;
  float *transmissivity;
  float *emissivity;
  float *blackbody_flux;
  float *LW_flux_div;
  float *SW_flux_div;
  float *IR_bandwidth;
  float *IR_wavelength;
  float *IR_ice_imag_idx;
  float *IR_ice_real_idx;
  float *LW_cloud_forcing;
  float *LW_mass_abs_coeff;
  float *cloud_optical_depth;
  float *albedo_of_time;
  float *emissivity_of_time;
  float *IWP_of_time;
  float *SW_cloud_forcing;
  float *altitude;
  float *capacitance;
  float *ccm2_cloud_cover;
  float *ccm2_env_pressure;
  float *ccm2_env_temperature;
  float *ccm2_ice_path;
  float *ccm2_mmr_O3;
  float *ccm2_mmr_vapor;
  float *cloud_base_of_time;
  float *cloud_cover;
  float *cloud_effective_radius;
  float *cloud_top_of_time;
  float *crystal_length;
  float *crystal_mass;
  float *delta_length;
  float *delta_mass;
  float *denv_temp_dt;
  float *dew_point_temp;
  float *dvapor_density_dt;
  float *effective_radius;
  float *env_density;
  float *env_pressure;
  float *env_pressure_int;
  float *env_temperature;
  float *eqm_vap_ice;
  float *eqm_vap_liquid;
  float *equiv_rad_cubed;
  float *equiv_rad_squared;
  float *equiv_radius;
  float *flux_up_SW;
  float *flux_down_SW;
  float *flux_up_LW;
  float *flux_down_LW;
  float *flux_up_LW_int;
  float *flux_down_LW_int;
  float *frost_point_temp;
  float *heating_rate_LW;
  float *heating_rate_SW;
  float *heating_rate_net;
  float *ice_imag_idx;
  float *IWP;
  float *IWC;
  float *ice_real_idx;
  float *latent_heating;
  float *max_length_of_time;
  float *mmr_O3;
  float *mmr_vapor;
  float *orig_env_temp;
  float *orig_mmr_vapor;
  float *potential_temp;
  float *pp_vapor;
  float *prism_radius;
  float *saturation_ice_of_time;
  float *saturation_ice;
  float *advective_heating;
  float *surface_area_of_time;
  float *thermal_conductivity;
  float *time_array;
  float *time_snapshot;
  float *total_layer_conc;
  float *vapor_density;
  float *vapor_density_new;
  float *vapor_density_swap;
  float *vapor_diffusivity;
  float *vapor_path;
  float *vapor_path_of_time;
  float *optical_depth;
  float *LW_abs_opt_depth;
  float *water_path_of_time;
  float *wavelength;
  float *wind_speed;
  float *wind_speed_new;
  float *wind_speed_swap;
  
  float ASNIR;
  float ASVIS;
  float AWNIR;
  float AWVIS;
  float CCN_scaling_diameter;
  float CCN_tot_conc;
  float IWC_init;
  float INC_init;
  float S_ice_ambient;
  float S_ice_cloud;
  float S_liquid_ambient;
  float SW_surf_cloud_forcing;
  float alpha_growth;
  float alpha_transport;
  float altitude_LBC;
  float altitude_RBC;
  float cloud_base;
  float cloud_thick;
  float conc_LBC;
  float conc_RBC;
  float curvature;
  float daltitude_dtime_LBC;
  float daltitude_dtime_RBC;
  float day_of_year;
  float denom_1;
  float denom_2;
  float denom_3;
  float dist_LBC;
  float dist_RBC;
  float dmass_dtime_LBC;
  float dmass_dtime_RBC;
  float dt;
  float dt_growth;
  float dt_net_fall_speed;
  float dt_wind_speed;
  float dz;
  float float_foo;
  float float_foo2;
  float float_foo3;
  float float_foo_input;
  float flux_up_LW_cloud_base_int;
  float flux_down_LW_cloud_top_int;
  float frac_strng_zen_ang_srf;
  float grid_base;
  float grid_thick;
  float ground_pressure;
  float latitude;
  float mass_LBC;
  float mass_RBC;
  float numerator_1;
  float numerator_2;
  float skin_temp;
  float smallest_mass;
  float smallest_length;
  float snow_cover;
  float solute_deliquescence_point;
  float surf_air_temp;
  float surf_roughness;
  float surf_type_flag;
  float tau_CCN_timescale;
  float vapor_LBC;
  float vapor_RBC;
  float wind_speed_max;
  float wind_speed_LBC;
  float wind_speed_RBC;
  
  int *Euler_stable_size;
  
  int agg_step;
  int bottom_ccm2_interface_level;
  int band;
  int bins_per_doubling;
  int cloud_base_idx;
  int cloud_top_idx;
  int etbfct_alpha;
  int atmospheric_profile_type;
  int frame;
  int ice_crystal_habit;
  int initial_distribution_type;
  int int_foo;
  int int_foo2;
  int layer;
  int new_size;
  int num_CCN_size;
  int num_Euler_stable_size;
  int num_LW_band;
  int num_band;
  int num_ccm2_level;
  int num_cloudy_layers;
  int num_frame;
  int num_layer;
  int num_size;
  int num_step;
  int opt;
  int plot_step;
  int rad_step;
  int size;
  int step; 
  int time_step;
  int top_ccm2_interface_level;
  
  time_t clock;
  
  int smallest_rimee_bin;
  int largest_rimee_bin;
  int smallest_rimer_bin;
  int largest_rimer_bin;
  float collision_kernel;
  float pi_over_four;
  float num_collisions;
  float num_collections;
  float aggregate_mass;
  float sticking_efficiency;

  /* netCDF declarations */ 
  int cdfid;
  
  /* containers for scalar attributes */
  char char_val;

  /* dimension ID's (which are ints) and array to hold them */
  int altitude_dim;
  int num_layerp2_dim, num_bandp2_dim, num_framep1_dim, num_sizep2_dim, num_stepp1_dim, num_CCN_sizep2_dim;
  int dims[MAX_NUM_NETCDF_DIMS];

  /* dimension sizes (which are longs) and array to hold them */ 
  long count[MAX_NUM_NETCDF_DIMS];
  long num_layerp2,num_bandp2,num_framep1,num_sizep2,num_stepp1,num_CCN_sizep2;

  /* indicial starting offsets (which are longs) */
  long start[MAX_NUM_NETCDF_DIMS];
  
  /* scalar variable ids */
  int dt_id, dz_id, num_CCN_size_id, num_band_id, num_ccm2_level_id, num_cloudy_layers_id, num_frame_id, num_layer_id, num_size_id, num_step_id, agg_step_id, rad_step_id, plot_step_id, PAUSE_id;

  /* two dimensional variable ids */
  int IWC_snapshot_id, number_snapshot_id, concentration_id, distribution_id, crystal_heat_SW_id, crystal_heat_LW_id, eff_asymmetry_param_id, eff_ss_albedo_id, eff_tau_extinction_id, daltitude_dtime_id;

  /* one dimensional variable ids */
  int saturation_liquid_id, aspect_ratio_id, crystal_diameter_id, CCN_activated_id, CCN_conc_id, CCN_diameter_id, CCN_mass_id, IWC_id, IWP_id, IWP_of_time_id, LW_cloud_forcing_id, SW_cloud_forcing_id, advective_heating_id, albedo_of_time_id, altitude_id, cloud_base_of_time_id, cloud_effective_radius_id, cloud_top_of_time_id, crystal_length_id, crystal_mass_id, delta_length_id, delta_mass_id, effective_radius_id, emissivity_of_time_id, env_density_id, env_pressure_id, env_temperature_id, flux_down_LW_id, flux_down_SW_id, flux_up_LW_id, flux_up_SW_id, heating_rate_LW_id, heating_rate_SW_id, heating_rate_net_id, latent_heating_id, max_length_of_time_id, mmr_vapor_id, cloud_optical_depth_id, orig_env_temp_id, orig_mmr_vapor_id, potential_temp_id, pp_vapor_id, saturation_ice_of_time_id, saturation_ice_id, surface_area_of_time_id, time_array_id, time_snapshot_id, vapor_density_id, vapor_path_id, vapor_path_of_time_id, optical_depth_id, water_path_of_time_id, wind_speed_id;
  
  /* Movie-specific id's */ 
  /* three dimensional variable ids */
  int movie_concentration_id;
  /* two dimensional variable ids */
  int movie_effective_radius_id;
  /* one dimensional variable ids */
  int movie_cloud_effective_radius_id;
    
  /* Aggregation parameters */ 
  smallest_rimee_bin = 10;
  largest_rimee_bin = 35;
  smallest_rimer_bin = 10;
  largest_rimer_bin = 35;
  sticking_efficiency = 1.;
  pi_over_four = M_PI/4.;

  /* set defaults */
  AGGREGATION = False; /* Option K */ 
  CCN_scaling_diameter = .075e-6;
  CCN_tot_conc = 200.e6;
  ERROR_FLAG = False; /* Option B */
  EULER = True; /* Option B */
  EVEN_STEP = False;
  HETEROGENEOUS_NUCLEATION = True; /* Option F */
  HOMOGENEOUS_NUCLEATION = True; /* Option L */
  ICE_ADVECT = True; /* Option C */
  INC_init = 1.e6; /* m^-3 */ /* Option N */
  IWC_init = .00001; /* kg/m^3 (10 mg/m^3) */ 
  LAGRANGE = True; /* Option A */
  LIOU = False; 
  MOVIE_NETCDF_OUTPUT = False; /* Option z */ 
  NETCDF_OUTPUT = True;
  PAUSE = True; /* ready to roll */ /* Option P */
  RADIATION = True; /* Option R */ 
  RAD_FEEDBACK = False; /* Option Q */ 
  S_ice_cloud = 1.2; /* ice sat. ratio in cloud */ /* Option S */
  S_liquid_ambient = .4; /* liquid sat. ratio beneath and above cloud */ /* Option s */
  UPDRAFT_PROFILE_IS_FLAT = False; /* Option u */
  VAPOR_ADVECT = True; /* Option U */
  VERBOSE = True; /* print out WARNINGS? */ /* Option V */
  WEED_NEGATIVE_VALUES = True; /* Option W */
  agg_step = 10; /* time steps per AGGREGATION call */ /* Option ? */
  alpha_growth = .4; /* Courant safety factor dt <= alpha*dz/max(w) */ /* Option a */
  alpha_transport = .2; /* Courant safety factor dt <= alpha*dz/max(w) */ /* Option j */
  bins_per_doubling = 2; /* resolution of mass mesh */ /* Option b */
  cloud_base = 1.e4; /* m */ /* Option H */
  cloud_thick = 1000.; /* m */ /* Option h */
  conc_LBC = 1.; /* boundary on crystal concentration due to fall speed advection */
  conc_RBC = 1.; /* boundary on crystal concentration due to fall speed advection */
  dist_LBC = 1.; /* boundary on crystal distribution due to growth advection */
  dist_RBC = 1.; /* boundary on crystal distribution due to growth advection */
  dt = 100.; /* dt must be initialized to arbitrary for Courant scheme to work */ 
  etbfct_alpha = 1.; /* Cartesian geometry */
  atmospheric_profile_type = TROPICS; /* Option M */ 
  float_foo_input = 0.; /* Option f */ 
  grid_base = 8000.; /* m */ /* Option G */
  grid_thick = 4000.; /* m */ /* Option g  */
  ice_crystal_habit = HEXAGONAL_COLUMNS;
  initial_distribution_type = KNOLLENBERG; /* Option Z */ 
  num_CCN_size = NUM_CCN_SIZES;
  num_LW_band = NUM_IR_SPECTRAL_BANDS;
  num_band = NUM_CCM2_SPECTRAL_BANDS;
  num_ccm2_level = 35; /* Option c */ 
  num_layer = 100; /* Option l */
  num_size = 40; /* Option x */
  num_step = 1;  /* Option n */
  plot_step = 20; /* time steps per graph */ /* Option p */
  rad_step = 1; /* time steps per RADIATION call */ /* Option r */
  smallest_length = 10.e-6; /* length in m of smallest mass bin */ /* Option m */
  solute_deliquescence_point = NH42SO4_deliquescence_point;
  tau_CCN_timescale = 1.; /* s */ /* Option t */
  time_step = 0; /* for initialization of time arrays */ 
  vapor_LBC = 1.; /* boundaries on vapor density due to wind advection */
  vapor_RBC = 1.;
  wind_speed_max = .05; /* m/s wind speed at initial base of grid */ /* Option w */
  
  (void)strcpy(in_file,"stdin"); /* Option i */
  (void)strcpy(out_file,"cld.nc"); /* Option o */
  (void)strcpy(err_file,"stderr"); /* Option e */
  
  /* parse command line arguments */
  option_string="Aa:Bb:Cc:D:d:e:Ff:G:g:H:h:i:J:j:Kk:L:l:M:m:N:n:o:P:p:QqRr:S:s:t:UuVvWw:Xx:Z:z";
  while((opt = getopt(argc,argv,option_string)) != EOF){
    switch(opt){
    case 'A':
      /* Toggle Lagrangian growth phase. Default is ON */
      LAGRANGE = !LAGRANGE;
      break;
    case 'a':
      /* The Courant condition safety factor dt <= alpha*dz/max(w) Default is .4 */
      alpha_growth = atof(optarg);
      break;
    case 'B':
      /* Toggle Eulerian growth phase. Default is ON */
      EULER = !EULER;
      break;
    case 'b':
      /* resolution of mass mesh. Default is 2 */
      bins_per_doubling = atoi(optarg);
      break;
    case 'C':
      /* Toggle velocity advection phase. Default is ON */
      ICE_ADVECT = !ICE_ADVECT;
      break;
    case 'c':
      /* number of externally specified levels for ccm2 rad code. Default is 35 */
      num_ccm2_level = atoi(optarg);
      break;
    case 'D':
      /* The debugging level.  Default is 0. */
      debug = atoi(optarg);
      break;
    case 'd':
      /* the spare int input for debugging. Default is 0 */
      debug_value = atoi(optarg);
      break;
    case 'e':
      /* get the error file name. Default is stderr */
      (void)strcpy(err_file,optarg);
      break;
    case 'F':
      /* Toggle the hetero. nuc. pzn. Default is True */ 
      HETEROGENEOUS_NUCLEATION = !HETEROGENEOUS_NUCLEATION;
      break;
    case 'f':
      /* Set the generic tuning parameter.  Default is 0. */
      float_foo_input = atof(optarg);
/*      solute_deliquescence_point = float_foo_input;*/
      S_liquid_ambient = float_foo_input;
      break;
    case 'G':
      /* Default is 8000. m */ 
      grid_base = atof(optarg);
      break;
    case 'g':
      /* Default is 4000. m */ 
      grid_thick = atof(optarg);
      break;
    case 'H':
      /* Default is 1.e4 m */ 
      cloud_base = atof(optarg);
      break;
    case 'h':
      /* Default is 1000. m */ 
      cloud_thick = atof(optarg);
      break;
    case 'i':
      /* get the input file name. Default is cloud.nc */
      (void)strcpy(in_file,optarg);
      break;
    case 'J':
      /* Toggle ice crystal habit. Default is HEXAGONAL_COLUMNS */
      ice_crystal_habit = atoi(optarg); 
      break;
    case 'j':
      /* The Courant condition safety factor dt <= alpha*dz/max(w) Default is .2 */
      alpha_transport = atof(optarg);
      break;
    case 'K':
      /* Toggle crystal aggregation (collisions). Default is ON */
      AGGREGATION=!AGGREGATION;
      break;
    case 'k':
      /* timestep. Courant transport condition overrides this. Default is 100. s */ 
      dt = atof(optarg);
      break;
    case 'L':
      /* Toggle HOMOGENEOUS_NUCLEATION. Default is True */
      HOMOGENEOUS_NUCLEATION = !HOMOGENEOUS_NUCLEATION;
      break;
    case 'l':
      /* Default is 100 */
      num_layer = atoi(optarg); 
      break;
    case 'M':
      /* Default is TROPICS */
      atmospheric_profile_type = atoi(optarg); 
      break;
    case 'm':
      /* Default is 10 microns */ 
      smallest_length = atof(optarg);
      break;
    case 'N':
      /* Default is 1.e6 m^-3 */ 
      INC_init = atof(optarg);
      break;
    case 'n':
      /* Default is 1 */
      num_step = atoi(optarg);
      break;
    case 'o':
      /* get the output file name. Default is stdout */
      (void)strcpy(out_file,optarg);
      break;
    case 'P':
      /* Wait in between graphs? Default is True */ 
      PAUSE=!PAUSE;
      break;
    case 'p':
      /* Default is 20 time steps per snapshot */ 
      plot_step = atoi(optarg);
      break;
    case 'Q':
      /* toggle flux convergence feedback into env T. Default is True */ 
      RAD_FEEDBACK=!RAD_FEEDBACK;
      break;
    case 'q':
      /* print out the debug options */
      print_debug();
      exit(1);
      break;
    case 'R':
      /* toggle running mie and ccm2 radiation routines. Default is True */ 
      RADIATION=!RADIATION;
      break;
    case 'r':
      /* Default is update radiation values every step */ 
      rad_step = atoi(optarg);
      break;
    case 'S':
      /* Default is 1.2 ice sat ratio in cloud */ 
      S_ice_cloud = atof(optarg);
      break;
    case 's':
      /* Default is .9 ice sat ratio beneath cloud */ 
      S_liquid_ambient = atof(optarg);
      break;
    case 't':
      /* Default is 1. s */ 
      tau_CCN_timescale = atof(optarg);
      break;
    case 'U':
      /* Toggle wind advection phase. Default is ON */
      VAPOR_ADVECT = False;
      break;
    case 'u':
      /* Toggle updraft profile shape. Default is OFF */
      UPDRAFT_PROFILE_IS_FLAT = !UPDRAFT_PROFILE_IS_FLAT;
      break;
    case 'v':
      /* print the RCS program info */
      (void)fprintf(stderr,rcs_Id);
      (void)fprintf(stderr,rcs_Revision);
      (void)fprintf(stderr,"$Author: zender $\n");
      (void)fprintf(stderr,"$Date$\n");
      (void)fprintf(stderr,"$Locker:  $\n");
      (void)fprintf(stderr,"$RCSfile: main.c,v $\n");
      (void)fprintf(stderr,"$Source: /home/zender/cvs/cld/main.c,v $\n");
      (void)fprintf(stderr,"$Id$\n");
      (void)fprintf(stderr,"$State: Exp $\n");
      exit(0);
      break;
    case 'V':
      /* toggle verbose printing out of WARNINGS. Default is True */ 
      VERBOSE=!VERBOSE;
      break;
    case 'w':
      /* m/s wind speed at initial base of cloud */ /* Option w */
      wind_speed_max = atof(optarg);
      break;
    case 'W':
      /* Weed out negative data values from ETBFCT? Default is True */ 
      WEED_NEGATIVE_VALUES=!WEED_NEGATIVE_VALUES;
      break;
    case 'X':
      /* Toggle Liou hex crystal treatment. Default is OFF */
      LIOU = !LIOU;
      break;
    case 'x':
      /* Default is 40 */
      num_size = atoi(optarg);
      break;
    case 'Z':
      /* Toggle initial ice crystal distribution. Default is KNOLLENBERG */
      initial_distribution_type = atoi(optarg); 
      break;
    case 'z':
      /* Turn on movie mode, dump data ever # timesteps */
      MOVIE_NETCDF_OUTPUT = !MOVIE_NETCDF_OUTPUT;
      break;
    case '?':
      /* print proper usage */
      (void)print_usage(option_string);
      exit(1);
    } /* end switch */
  } /* end while loop */
  
  /* Attempt to set up IEEE signal handling */ 
#if ( defined sun )
  if(ieee_handler("set","common",(sigfpe_handler_type)ieee_exception_handler)!=0){
    (void)fprintf(stdout,"IEEE trapping not supported here.\n");
  } /* end if */
#endif

  /* start the clock and save the command line */  
  string_cmdline(argc,argv,cmdline,CMDLINE_SIZE);
  (void)fprintf(stdout,"Command Line: %s\n",cmdline);
  clock=time((time_t *)NULL);
  time_buf_start=ctime(&clock);
  (void)fprintf(stdout,"\tstart = %s",time_buf_start);
  strcat(cmdline," Version ");
  char_ptr_foo2=strchr(rcs_Revision,':');
  strncat(cmdline,char_ptr_foo2+2,4);
  strcat(cmdline," ");
  strcat(cmdline,time_buf_start);
  
  /* any flexible sizes to be defined before the malloc'ing? */ 
  num_frame = num_step/plot_step;

  /* NB: CRAY doesn't like these huge if(malloc...) statements, they
   exhaust the heap space, so break them up into smaller statements, 30
   malloc's is ok, but 60 is too many */ 
  if(
     ((saturation_liquid=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((aspect_ratio=(float *)malloc((num_size+2)*sizeof(float))) == NULL ) ||
     ((crystal_diameter=(float *)malloc((num_size+2)*sizeof(float))) == NULL ) ||
     ((CCN_activated=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((CCN_conc=(float *)malloc((num_CCN_size+2)*sizeof(float))) == NULL ) ||
     ((CCN_diameter=(float *)malloc((num_CCN_size+2)*sizeof(float))) == NULL ) ||
     ((CCN_mass=(float *)malloc((num_CCN_size+2)*sizeof(float))) == NULL ) ||
     ((LW_flux_div=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((SW_flux_div=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((IR_wavelength=(float *)malloc((num_LW_band+2)*sizeof(float))) == NULL ) ||
     ((IR_bandwidth=(float *)malloc((num_LW_band+2)*sizeof(float))) == NULL ) ||
     ((IR_ice_imag_idx=(float *)malloc((num_LW_band+2)*sizeof(float))) == NULL ) ||
     ((IR_ice_real_idx=(float *)malloc((num_LW_band+2)*sizeof(float))) == NULL ) ||
     ((SW_cloud_forcing=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((LW_cloud_forcing=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((LW_mass_abs_coeff=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((cloud_optical_depth=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((albedo_of_time=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((emissivity_of_time=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((IWP_of_time=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((altitude=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((blackbody_flux=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((capacitance=(float *)malloc((num_size+2)*sizeof(float))) == NULL ) ||
     ((cloud_base_of_time=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((cloud_cover=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((cloud_effective_radius=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((cloud_top_of_time=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((crystal_length=(float *)malloc((num_size+2)*sizeof(float))) == NULL ) ||
     ((crystal_mass=(float *)malloc((num_size+2)*sizeof(float))) == NULL ) ||
     ((delta_length=(float *)malloc((num_size+2)*sizeof(float))) == NULL ) ||
     ((delta_mass=(float *)malloc((num_size+2)*sizeof(float))) == NULL ) ||
     ((denv_temp_dt=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((dew_point_temp=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((dvapor_density_dt=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((effective_radius=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((emissivity=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     False ){
    (void)fprintf(stdout,"Unable to allocate array in main\n");
    exit(1);
  } /* end if */

  if(
     ((env_density=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((env_pressure=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((env_pressure_int=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((env_temperature=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((eqm_vap_ice=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((eqm_vap_liquid=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((equiv_rad_cubed=(float *)malloc((num_size+2)*sizeof(float))) == NULL ) ||
     ((equiv_rad_squared=(float *)malloc((num_size+2)*sizeof(float))) == NULL ) ||
     ((equiv_radius=(float *)malloc((num_size+2)*sizeof(float))) == NULL ) ||
     ((flux_up_SW=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((flux_down_SW=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((flux_up_LW=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((flux_down_LW=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((flux_up_LW_int=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((flux_down_LW_int=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((frost_point_temp=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((hex_column_diameter=(float *)malloc((num_size+2)*sizeof(float))) == NULL ) ||
     ((hex_column_length=(float *)malloc((num_size+2)*sizeof(float))) == NULL ) ||
     ((hex_equiv_radius=(float *)malloc((num_size+2)*sizeof(float))) == NULL ) ||
     ((hex_equiv_rad_squared=(float *)malloc((num_size+2)*sizeof(float))) == NULL ) ||
     ((heating_rate_LW=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((heating_rate_SW=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((heating_rate_net=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((ice_imag_idx=(float *)malloc((num_band+2)*sizeof(float))) == NULL ) ||
     ((IWP=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((IWC=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((ice_real_idx=(float *)malloc((num_band+2)*sizeof(float))) == NULL ) ||
     ((latent_heating=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     False ){
    (void)fprintf(stdout,"Unable to allocate array in main\n");
    exit(1);
  } /* end if */

  if(
     ((max_length_of_time=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((mmr_O3=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((mmr_vapor=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((orig_env_temp=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((orig_mmr_vapor=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((potential_temp=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((pp_vapor=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((prism_radius=(float *)malloc((num_size+2)*sizeof(float))) == NULL ) ||
     ((saturation_ice_of_time=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((saturation_ice=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((advective_heating=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((surface_area_of_time=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((thermal_conductivity=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((time_array=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((time_snapshot=(float *)malloc((num_frame+1)*sizeof(float))) == NULL ) ||
     ((total_layer_conc=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((transmissivity=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((vapor_density=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((vapor_density_new=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((vapor_diffusivity=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((vapor_path=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((vapor_path_of_time=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((optical_depth=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((LW_abs_opt_depth=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((water_path_of_time=(float *)malloc((num_step+1)*sizeof(float))) == NULL ) ||
     ((wavelength=(float *)malloc((num_band+2)*sizeof(float))) == NULL ) ||
     ((wind_speed=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     ((wind_speed_new=(float *)malloc((num_layer+2)*sizeof(float))) == NULL ) ||
     False ){
    (void)fprintf(stdout,"Unable to allocate array in main\n");
    exit(1);
  } /* end if */
  
  if(
     ((ccm2_cloud_cover=
       (float *)malloc((num_ccm2_level+num_layer+2)*sizeof(float))) == NULL ) ||
     ((ccm2_env_pressure=
       (float *)malloc((num_ccm2_level+num_layer+2)*sizeof(float))) == NULL ) ||
     ((ccm2_env_temperature=
       (float *)malloc((num_ccm2_level+num_layer+2)*sizeof(float))) == NULL ) ||
     ((ccm2_ice_path=
       (float *)malloc((num_ccm2_level+num_layer+2)*sizeof(float))) == NULL ) ||
     ((ccm2_mmr_O3=
       (float *)malloc((num_ccm2_level+num_layer+2)*sizeof(float))) == NULL ) ||
     ((ccm2_mmr_vapor=
       (float *)malloc((num_ccm2_level+num_layer+2)*sizeof(float))) == NULL ) ||
     False ){
    (void)fprintf(stdout,"Unable to allocate array in main\n");
    exit(1);
  } /* end if */
  
  if(
     ((Euler_stable_size=(int *)malloc((num_layer+2)*sizeof(int))) == NULL ) ||
     False ){
    (void)fprintf(stdout,"Unable to allocate array in main\n");
    exit(1);
  } /* end if */
  
  if(
     ((LW_asymmetry_param_ptr=
       (float *)malloc((num_LW_band+2)*(num_size+2)*sizeof(float))) == NULL ) ||
     ((LW_size_param_ptr=
       (float *)malloc((num_LW_band+2)*(num_size+2)*sizeof(float))) == NULL ) ||
     ((LW_scattering_eff_ptr=
       (float *)malloc((num_LW_band+2)*(num_size+2)*sizeof(float))) == NULL ) ||
     ((LW_extinction_eff_ptr=
       (float *)malloc((num_LW_band+2)*(num_size+2)*sizeof(float))) == NULL ) ||
     ((LW_planck_radiance_ptr=
       (float *)malloc((num_layer+2)*(num_LW_band+2)*sizeof(float))) == NULL ) ||
     ((IWC_snapshot_ptr=
       (float *)malloc((num_frame+1)*(num_layer+2)*sizeof(float))) == NULL ) ||
     ((number_snapshot_ptr=
       (float *)malloc((num_frame+1)*(num_layer+2)*sizeof(float))) == NULL ) ||
     ((asymmetry_param_ptr=
       (float *)malloc((num_band+2)*(num_size+2)*sizeof(float))) == NULL ) ||
     ((size_param_ptr=
       (float *)malloc((num_band+2)*(num_size+2)*sizeof(float))) == NULL ) ||
     ((crystal_temperature_ptr=
       (float *)malloc((num_layer+2)*(num_size+2)*sizeof(float))) == NULL ) ||
     ((crystal_heat_net_ptr=
       (float *)malloc((num_layer+2)*(num_size+2)*sizeof(float))) == NULL ) ||
     ((crystal_heat_SW_ptr=
       (float *)malloc((num_layer+2)*(num_size+2)*sizeof(float))) == NULL ) ||
     ((crystal_heat_LW_ptr=
       (float *)malloc((num_layer+2)*(num_size+2)*sizeof(float))) == NULL ) ||
     ((absorption_eff_ptr=
       (float *)malloc((num_band+2)*(num_size+2)*sizeof(float))) == NULL ) ||
     ((absorption_x_sec_ptr=
       (float *)malloc((num_band+2)*(num_size+2)*sizeof(float))) == NULL ) ||
     ((extinction_x_sec_ptr=
       (float *)malloc((num_band+2)*(num_size+2)*sizeof(float))) == NULL ) ||
     ((scattering_x_sec_ptr=
       (float *)malloc((num_band+2)*(num_size+2)*sizeof(float))) == NULL ) ||
     ((scattering_eff_ptr=
       (float *)malloc((num_band+2)*(num_size+2)*sizeof(float))) == NULL ) ||
     ((extinction_eff_ptr=
       (float *)malloc((num_band+2)*(num_size+2)*sizeof(float))) == NULL ) ||
     ((concentration_ptr=
       (float *)malloc((num_layer+2)*(num_size+2)*sizeof(float))) == NULL ) ||
     ((concentration_new_ptr=
       (float *)malloc((num_layer+2)*(num_size+2)*sizeof(float))) == NULL ) ||
     ((conc_transpose_ptr=
       (float *)malloc((num_size+2)*(num_layer+2)*sizeof(float))) == NULL ) ||
     ((conc_transpose_new_ptr=
       (float *)malloc((num_size+2)*(num_layer+2)*sizeof(float))) == NULL ) ||
     ((daltitude_dtime_ptr=
       (float *)malloc((num_size+2)*(num_layer+2)*sizeof(float))) == NULL ) ||
     ((distribution_ptr=
       (float *)malloc((num_layer+2)*(num_size+2)*sizeof(float))) == NULL ) ||
     ((distribution_new_ptr=
       (float *)malloc((num_layer+2)*(num_size+2)*sizeof(float))) == NULL ) ||
     False ){
    (void)fprintf(stdout,"Unable to allocate array in main\n");
    exit(1);
  } /* end if */

  if(
     ((SW_spec_flux_div_ptr=
       (float *)malloc((num_layer+2)*(num_band+2)*sizeof(float))) == NULL ) ||
     ((LW_spec_vol_abs_coeff_ptr=
       (float *)malloc((num_layer+2)*(num_LW_band+2)*sizeof(float))) == NULL ) ||
     ((LW_absorption_eff_ptr=
       (float *)malloc((num_LW_band+2)*(num_size+2)*sizeof(float))) == NULL ) ||
     ((LW_absorption_x_sec_ptr=
       (float *)malloc((num_LW_band+2)*(num_size+2)*sizeof(float))) == NULL ) ||
     ((LW_extinction_x_sec_ptr=
       (float *)malloc((num_LW_band+2)*(num_size+2)*sizeof(float))) == NULL ) ||
     ((LW_scattering_x_sec_ptr=
       (float *)malloc((num_LW_band+2)*(num_size+2)*sizeof(float))) == NULL ) ||
     ((LW_spec_band_weight_ptr=
       (float *)malloc((num_layer+2)*(num_LW_band+2)*sizeof(float))) == NULL ) ||
     ((dmass_dtime_ptr=
       (float *)malloc((num_layer+2)*(num_size+2)*sizeof(float))) == NULL ) ||
     ((eff_asymmetry_param_ptr=
       (float *)malloc((num_layer+2)*(num_band+2)*sizeof(float))) == NULL ) ||
     ((eff_ss_albedo_ptr=
       (float *)malloc((num_layer+2)*(num_band+2)*sizeof(float))) == NULL ) ||
     ((eff_tau_extinction_ptr=
       (float *)malloc((num_layer+2)*(num_band+2)*sizeof(float))) == NULL ) ||
     ((eff_tau_scatter_ptr=
       (float *)malloc((num_layer+2)*(num_band+2)*sizeof(float))) == NULL ) ||
     ((fall_speed_ptr=
       (float *)malloc((num_size+2)*(num_layer+2)*sizeof(float))) == NULL ) ||
     ((ice_density_ptr=
       (float *)malloc((num_layer+2)*(num_size+2)*sizeof(float))) == NULL ) ||
     ((diffusivity_correction_ptr=
       (float *)malloc((num_layer+2)*(num_size+2)*sizeof(float))) == NULL ) ||
     ((conductivity_correction_ptr=
       (float *)malloc((num_layer+2)*(num_size+2)*sizeof(float))) == NULL ) ||
     ((lagrange_mass_ptr=
       (float *)malloc((num_layer+2)*(num_size+2)*sizeof(float))) == NULL ) ||
     ((planck_avg_abs_x_sec_ptr=
       (float *)malloc((num_layer+2)*(num_size+2)*sizeof(float))) == NULL ) ||
     False ){
    (void)fprintf(stdout,"Unable to allocate array in main\n");
    exit(1);
  } /* end if */
  
  if(
     ((SW_spec_flux_div=
       (float **)malloc((num_layer+2)*sizeof(float *))) == NULL ) ||
     ((LW_asymmetry_param=
       (float **)malloc((num_LW_band+2)*sizeof(float *))) == NULL ) ||
     ((LW_size_param=
       (float **)malloc((num_LW_band+2)*sizeof(float *))) == NULL ) ||
     ((LW_scattering_eff=
       (float **)malloc((num_LW_band+2)*sizeof(float *))) == NULL ) ||
     ((LW_extinction_eff=
       (float **)malloc((num_LW_band+2)*sizeof(float *))) == NULL ) ||
     ((LW_planck_radiance=
       (float **)malloc((num_layer+2)*sizeof(float *))) == NULL ) ||
     ((IWC_snapshot=
       (float **)malloc((num_frame+1)*sizeof(float *))) == NULL ) ||
     ((number_snapshot=
       (float **)malloc((num_frame+1)*sizeof(float *))) == NULL ) ||
     ((asymmetry_param=
       (float **)malloc((num_band+2)*sizeof(float *))) == NULL ) ||
     ((size_param=
       (float **)malloc((num_band+2)*sizeof(float *))) == NULL ) ||
     ((crystal_temperature=
       (float **)malloc((num_layer+2)*sizeof(float *))) == NULL ) ||
     ((crystal_heat_net=
       (float **)malloc((num_layer+2)*sizeof(float *))) == NULL ) ||
     ((crystal_heat_SW=
       (float **)malloc((num_layer+2)*sizeof(float *))) == NULL ) ||
     ((crystal_heat_LW=
       (float **)malloc((num_layer+2)*sizeof(float *))) == NULL ) ||
     ((absorption_eff=
       (float **)malloc((num_band+2)*sizeof(float *))) == NULL ) ||
     ((absorption_x_sec=
       (float **)malloc((num_band+2)*sizeof(float *))) == NULL ) ||
     ((extinction_x_sec=
       (float **)malloc((num_band+2)*sizeof(float *))) == NULL ) ||
     ((scattering_x_sec=
       (float **)malloc((num_band+2)*sizeof(float *))) == NULL ) ||
     ((scattering_eff=
       (float **)malloc((num_band+2)*sizeof(float *))) == NULL ) ||
     ((extinction_eff=
       (float **)malloc((num_band+2)*sizeof(float *))) == NULL ) ||
     ((concentration=
       (float **)malloc((num_layer+2)*sizeof(float *))) == NULL ) ||
     ((concentration_new=
       (float **)malloc((num_layer+2)*sizeof(float *))) == NULL ) ||
     ((conc_transpose=
       (float **)malloc((num_size+2)*sizeof(float *))) == NULL ) ||
     ((conc_transpose_new=
       (float **)malloc((num_size+2)*sizeof(float *))) == NULL ) ||
     ((daltitude_dtime=
       (float **)malloc((num_size+2)*sizeof(float *))) == NULL ) ||
     ((distribution=
       (float **)malloc((num_layer+2)*sizeof(float *))) == NULL ) ||
     ((distribution_new=
       (float **)malloc((num_layer+2)*sizeof(float *))) == NULL ) ||
     ((dmass_dtime=
       (float **)malloc((num_layer+2)*sizeof(float *))) == NULL ) ||
     ((eff_asymmetry_param=
       (float **)malloc((num_layer+2)*sizeof(float *))) == NULL ) ||
     ((eff_ss_albedo=
       (float **)malloc((num_layer+2)*sizeof(float *))) == NULL ) ||
     ((eff_tau_extinction=
       (float **)malloc((num_layer+2)*sizeof(float *))) == NULL ) ||
     ((eff_tau_scatter=
       (float **)malloc((num_layer+2)*sizeof(float *))) == NULL ) ||
     ((fall_speed=
       (float **)malloc((num_size+2)*sizeof(float *))) == NULL ) ||
     ((ice_density=
       (float **)malloc((num_layer+2)*sizeof(float *))) == NULL ) ||
     ((diffusivity_correction=
       (float **)malloc((num_layer+2)*sizeof(float *))) == NULL ) ||
     ((conductivity_correction=
       (float **)malloc((num_layer+2)*sizeof(float *))) == NULL ) ||
     ((lagrange_mass=
       (float **)malloc((num_layer+2)*sizeof(float *))) == NULL ) ||
     ((LW_spec_vol_abs_coeff=
       (float **)malloc((num_layer+2)*sizeof(float *))) == NULL ) ||
     ((LW_absorption_eff=
       (float **)malloc((num_LW_band+2)*sizeof(float *))) == NULL ) ||
     ((LW_absorption_x_sec=
       (float **)malloc((num_LW_band+2)*sizeof(float *))) == NULL ) ||
     ((LW_extinction_x_sec=
       (float **)malloc((num_LW_band+2)*sizeof(float *))) == NULL ) ||
     ((LW_scattering_x_sec=
       (float **)malloc((num_LW_band+2)*sizeof(float *))) == NULL ) ||
     ((LW_spec_band_weight=
       (float **)malloc((num_layer+2)*sizeof(float *))) == NULL ) ||
     ((planck_avg_abs_x_sec=
       (float **)malloc((num_layer+2)*sizeof(float *))) == NULL ) ||
     False ){
    (void)fprintf(stdout,"Unable to allocate array in main\n");
    exit(1);
  } /* end if */
  
  /* REMEMBER: these matrices are indexed as [layer][size] */
  for(layer=0;layer<num_layer+2;layer++){
    concentration[layer]=concentration_ptr+(num_size+2)*layer;    
    concentration_new[layer]=concentration_new_ptr+(num_size+2)*layer;    
    crystal_temperature[layer]=crystal_temperature_ptr+(num_size+2)*layer;    
    crystal_heat_net[layer]=crystal_heat_net_ptr+(num_size+2)*layer;    
    crystal_heat_SW[layer]=crystal_heat_SW_ptr+(num_size+2)*layer;    
    crystal_heat_LW[layer]=crystal_heat_LW_ptr+(num_size+2)*layer;    
    distribution[layer]=distribution_ptr+(num_size+2)*layer;    
    distribution_new[layer]=distribution_new_ptr+(num_size+2)*layer;    
    dmass_dtime[layer]=dmass_dtime_ptr+(num_size+2)*layer;
    ice_density[layer]=ice_density_ptr+(num_size+2)*layer;
    diffusivity_correction[layer]=diffusivity_correction_ptr+(num_size+2)*layer;
    conductivity_correction[layer]=conductivity_correction_ptr+(num_size+2)*layer;
    lagrange_mass[layer]=lagrange_mass_ptr+(num_size+2)*layer;
    planck_avg_abs_x_sec[layer]=planck_avg_abs_x_sec_ptr+(num_size+2)*layer;
  } /* end loop over layers */
  
  /* REMEMBER: these matrices are indexed as [layer][band] */
  for(layer=0;layer<num_layer+2;layer++){
    SW_spec_flux_div[layer]=SW_spec_flux_div_ptr+(num_band+2)*layer;
    eff_asymmetry_param[layer]=eff_asymmetry_param_ptr+(num_band+2)*layer;    
    eff_ss_albedo[layer]=eff_ss_albedo_ptr+(num_band+2)*layer;    
    eff_tau_extinction[layer]=eff_tau_extinction_ptr+(num_band+2)*layer;    
    eff_tau_scatter[layer]=eff_tau_scatter_ptr+(num_band+2)*layer;    
  } /* end loop over layers */
  
  /* REMEMBER: these are the only LW matrices indexed as [layer][band] */
  for(layer=0;layer<num_layer+2;layer++){
    LW_spec_vol_abs_coeff[layer]=LW_spec_vol_abs_coeff_ptr+(num_LW_band+2)*layer;
    LW_spec_band_weight[layer]=LW_spec_band_weight_ptr+(num_LW_band+2)*layer;
    LW_planck_radiance[layer]=LW_planck_radiance_ptr+(num_LW_band+2)*layer;    
  } /* end loop over layers */

  /* REMEMBER: these are the only LW matrices indexed as [band][size] */
  for(band=0;band<num_LW_band+2;band++){
    LW_absorption_x_sec[band]=LW_absorption_x_sec_ptr+(num_size+2)*band;    
    LW_extinction_x_sec[band]=LW_extinction_x_sec_ptr+(num_size+2)*band;    
    LW_scattering_x_sec[band]=LW_scattering_x_sec_ptr+(num_size+2)*band;    
    LW_absorption_eff[band]=LW_absorption_eff_ptr+(num_size+2)*band;
    LW_asymmetry_param[band]=LW_asymmetry_param_ptr+(num_size+2)*band;    
    LW_size_param[band]=LW_size_param_ptr+(num_size+2)*band;    
    LW_scattering_eff[band]=LW_scattering_eff_ptr+(num_size+2)*band;    
    LW_extinction_eff[band]=LW_extinction_eff_ptr+(num_size+2)*band;    
  } /* end loop over bands */
  
  /* REMEMBER: these are the only matrices indexed as [size][layer] */
  for(size=0;size<num_size+2;size++){
    conc_transpose[size]=conc_transpose_ptr+(num_layer+2)*size;    
    conc_transpose_new[size]=conc_transpose_new_ptr+(num_layer+2)*size;    
    daltitude_dtime[size]=daltitude_dtime_ptr+(num_layer+2)*size;
    fall_speed[size]=fall_speed_ptr+(num_layer+2)*size;
  } /* end loop over sizes */
  
  /* REMEMBER: these are the only SW matrices indexed as [band][size] */
  for(band=0;band<num_band+2;band++){
    asymmetry_param[band]=asymmetry_param_ptr+(num_size+2)*band;    
    size_param[band]=size_param_ptr+(num_size+2)*band;    
    absorption_eff[band]=absorption_eff_ptr+(num_size+2)*band;    
    absorption_x_sec[band]=absorption_x_sec_ptr+(num_size+2)*band;    
    extinction_x_sec[band]=extinction_x_sec_ptr+(num_size+2)*band;    
    scattering_x_sec[band]=scattering_x_sec_ptr+(num_size+2)*band;    
    scattering_eff[band]=scattering_eff_ptr+(num_size+2)*band;    
    extinction_eff[band]=extinction_eff_ptr+(num_size+2)*band;    
  } /* end loop over bands */
  
  /* REMEMBER: these are the only matrices indexed as [frame][layer] */
  for(frame=0;frame<num_frame+1;frame++){
    IWC_snapshot[frame]=IWC_snapshot_ptr+(num_layer+2)*frame;    
    number_snapshot[frame]=number_snapshot_ptr+(num_layer+2)*frame;    
  } /* end loop over frames */
  
  /* Bootstrap variables */
  if(True){

    if(ice_crystal_habit == HEXAGONAL_COLUMNS){
      crystal_diameter[1]=smallest_length/Aspect_ratio(smallest_length,&float_foo);
      smallest_mass=3.*sqrt(3.)*crystal_diameter[1]*crystal_diameter[1]*
      smallest_length*density_of_ice/8.;
    }else if(ice_crystal_habit == HEXAGONAL_BULLETS){
      smallest_mass=mass_of_bullet_length(smallest_length);
    }else if(ice_crystal_habit == SPHERES){
      smallest_mass=M_PI*pow(smallest_length,3.)*density_of_ice/6.;
    } /* end else */ 
    for(size=1;size<=num_size;size++){
      crystal_mass[size]=smallest_mass*
	pow(2.,((size-1.)/(float)bins_per_doubling));
    } /* end loop over sizes */
    float_foo=pow(2.,(1./bins_per_doubling));
    mass_LBC=smallest_mass*(1.+1./float_foo)/2.;
    mass_RBC=crystal_mass[num_size]*(1.+float_foo)/2.;
    
    for(size=2;size<=num_size-1;size++){
      delta_mass[size]=.5*(crystal_mass[size+1]-crystal_mass[size-1]);
    } /* end loop over sizes */
    delta_mass[1]=.5*(crystal_mass[1+1]+crystal_mass[1])-mass_LBC;
    delta_mass[num_size]=mass_RBC-
      .5*(crystal_mass[num_size]+crystal_mass[num_size-1]);
    
    /* Recall mass_LBC and mass_RBC are really interface mass at the
       boundaries of the grid, they're just stored here for convenience. */ 
    crystal_mass[0]=mass_LBC;
    crystal_mass[num_size+1]=mass_RBC;
    for(size=0;size<=num_size+1;size++){
      hex_column_dimensions_of_mass
	(crystal_mass[size],
	 hex_column_length+size,
	 hex_column_diameter+size,
	 hex_equiv_radius+size,
	 hex_equiv_rad_squared+size);
      if(ice_crystal_habit == HEXAGONAL_COLUMNS){
	crystal_length[size]=hex_column_length[size];
	crystal_diameter[size]=hex_column_diameter[size];
	aspect_ratio[size]=Aspect_ratio(crystal_length[size],&float_foo);
	equiv_radius[size]=hex_equiv_radius[size];
	equiv_rad_squared[size]=hex_equiv_rad_squared[size];
	prism_radius[size]=sqrt(3.)*hex_column_diameter[size]/4.;
      }else if(ice_crystal_habit == HEXAGONAL_BULLETS){
	crystal_length[size]=bullet_length_of_mass(crystal_mass[size]);
	crystal_diameter[size]=
	  2.*prism_radius_of_bullet_length(crystal_length[size]);
	aspect_ratio[size]=crystal_length[size]/crystal_diameter[size];
	equiv_radius[size]=bullet_length_to_equiv_rad(crystal_length[size]);
	equiv_rad_squared[size]=equiv_radius[size]*equiv_radius[size];
	prism_radius[size]=crystal_diameter[size]/2.;
      }else if(ice_crystal_habit == SPHERES){
	float_foo=6.*crystal_mass[size]/(density_of_ice*M_PI);
	crystal_length[size]=pow(float_foo,1./3.);
	crystal_diameter[size]=crystal_length[size];
	aspect_ratio[size]=1.;
	equiv_radius[size]=crystal_length[size]/2.;
	equiv_rad_squared[size]=equiv_radius[size]*equiv_radius[size];
	prism_radius[size]=equiv_radius[size];
      } /* end else */ 
      equiv_rad_cubed[size]=equiv_rad_squared[size]*equiv_radius[size];
      capacitance[size]=Capacitance(crystal_length[size],crystal_diameter[size]);
    } /* end loop over sizes */

    for(size=2;size<=num_size-1;size++){
      delta_length[size]=.5*(crystal_length[size+1]-crystal_length[size-1]);
    } /* end loop over sizes */
    delta_length[1]=.5*(crystal_length[1+1]+crystal_length[1])-crystal_length[0];
    delta_length[num_size]=crystal_length[num_size+1]-
      .5*(crystal_length[num_size]+crystal_length[num_size-1]);

    dz=(float)grid_thick/num_layer;
    for(layer=1;layer<=num_layer;layer++){
      altitude[layer]=grid_base+(layer-1)*dz+dz/2.; /* m */
    } /* end loop over layers */
    altitude_LBC=altitude[0]=grid_base;
    altitude_RBC=altitude[num_layer+1]=grid_base+grid_thick;
    
    init_temp_pres
      (altitude,env_pressure,env_temperature,num_layer,atmospheric_profile_type);
    if(min_VEC(env_temperature+1,num_layer) < -213.){
      (void)fprintf(stdout,"WARNING: size_dist_of_IWC() forced to kludge data for low temperature < -60.C\n");
    } /* endif */ 

    /* Set the saturation vapor pressures based on the layer temperature */ 
    for(layer=1;layer<=num_layer;layer++){
      eqm_vap_ice[layer]=Eqm_Vap_Ice(env_temperature[layer]);
      eqm_vap_liquid[layer]=Eqm_Vap_Liquid(env_temperature[layer]);
    } /* end loop over layers */

#undef NOT_DEFINED
#ifdef NOT_DEFINED
    /* This is the algorithm used in my Comps II */ 
    for(layer=1;layer<=num_layer;layer++){
      /* if the layer falls within the cloud then initialize it with a higher
	 water vapor content */
      if((altitude[layer] >= cloud_base) && 
	 (altitude[layer] <= cloud_base+cloud_thick)){
	/* CAUTION:  the Cray standard math lib pow() function doesn't
	   like raising a negative number to a fractional power! Will crash! */ 
	float_foo=(cloud_base+cloud_thick/2.-altitude[layer])*
	  (cloud_base+cloud_thick/2.-altitude[layer]);
	pp_vapor[layer]=S_liquid_ambient*eqm_vap_liquid[layer];
	S_ice_ambient=pp_vapor[layer]/eqm_vap_ice[layer];
	saturation_ice[layer]=S_ice_ambient+(S_ice_cloud-S_ice_ambient)*
	  exp(-10.*float_foo/(cloud_thick*cloud_thick)); 
      }else{
	saturation_liquid[layer]=S_liquid_ambient;
	pp_vapor[layer]=saturation_liquid[layer]*eqm_vap_liquid[layer];
	saturation_ice[layer]=pp_vapor[layer]/eqm_vap_ice[layer];
      } /* end else */ 
    } /* end loop over layers */
#endif
    
#undef NOT_DEFINED
#ifdef NOT_DEFINED
    /* This is the algorithm used for the stable ice pzns in ANV */ 
    /* This block uses unsmoothed S_ice_cloud in the cloud, and S_liquid_ambient outside the cloud. */ 
    for(layer=1;layer<=num_layer;layer++){
      /* if the layer falls within the cloud then initialize it with a higher
	 water vapor content */
      if((altitude[layer] >= cloud_base) && 
	 (altitude[layer] <= cloud_base+cloud_thick)){
	saturation_ice[layer]=S_ice_cloud;
	pp_vapor[layer]=S_ice_cloud*eqm_vap_ice[layer];
	saturation_liquid[layer]=pp_vapor[layer]/eqm_vap_liquid[layer];
      }else{
	saturation_liquid[layer]=S_liquid_ambient;
	pp_vapor[layer]=saturation_liquid[layer]*eqm_vap_liquid[layer];
	saturation_ice[layer]=pp_vapor[layer]/eqm_vap_ice[layer];
      } /* end else */ 
    } /* end loop over layers */
#endif
    
/*#undef NOT_DEFINED*/
/*#ifdef NOT_DEFINED*/
    /* This block defines the values that used to be written by the above block.
       This block uses S_ice_cloud for the entire column profile. */ 
    for(layer=1;layer<=num_layer;layer++){
	saturation_ice[layer]=S_ice_cloud;
	pp_vapor[layer]=S_ice_cloud*eqm_vap_ice[layer];
	saturation_liquid[layer]=pp_vapor[layer]/eqm_vap_liquid[layer];
    } /* end loop over layers */
/*#endif*/
    
    init_CCN_conc
      (CCN_mass,
       CCN_diameter,
       CCN_conc,
       CCN_tot_conc,
       CCN_scaling_diameter,
       num_CCN_size);

    for(layer=1;layer<=num_layer;layer++){
      /* only put particles in the cloud, not the ambient air */
      if((altitude[layer] >= cloud_base) && 
	 (altitude[layer] <= cloud_base+cloud_thick)){
	/* This is the distance to the center of the cloud in meters */ 
	float_foo2=altitude[layer]-(cloud_base+cloud_thick/2.);
	/* This is the distance to center as a fraction of the radius */ 
	float_foo3=2.*float_foo2/cloud_thick;
	for(size=1;size<=num_size;size++){
	  /* if there are "too few" particles, then initialize to zero, this
	     helps make sure that the big particles have some growing room */
	  
	  /* CAUTION:  the Cray standard math lib pow() function doesn't
	     like raising a negative number to a fractional power! Will crash! */ 
	  if(initial_distribution_type == HEYMSFIELD){
	    concentration[layer][size]=delta_length[size]*
	      size_dist_of_IWC(crystal_length[size],env_temperature[layer]);
/*		*exp(-10.*float_foo2*float_foo2/(cloud_thick*cloud_thick));*/
	  }else if(initial_distribution_type == DOWLING_RADKE){
	    concentration[layer][size]=INC_init*delta_length[size]*
	      PDF_of_length
		(crystal_length[size],initial_distribution_type,float_foo3);
/*		  *exp(-10.*float_foo2*float_foo2/(cloud_thick*cloud_thick));*/
	  }else if(initial_distribution_type == KNOLLENBERG){
	    concentration[layer][size]=INC_init*delta_length[size]*
	      PDF_of_length
		(crystal_length[size],initial_distribution_type,float_foo3);
	    if(crystal_length[size] < 30.e-6 && float_foo3 > 0.){
	      concentration[layer][size]*=
		(30.-MICRONS_PER_METER*crystal_length[size])*.05*float_foo3+1.;
	    } /* end if */
	  } /* end else */
	  if(concentration[layer][size] < MIN_CONC_FOR_INIT ||
	     crystal_length[size] > 600.e-6 || size > num_size-5){
	    concentration[layer][size]=0.;
	  } /* end if */
	  distribution[layer][size]=concentration[layer][size]/delta_mass[size];
	} /* end loop over sizes */
      }else{  /* make sure the non cloud-bins are initialized to zero */
	(void)memset((char *)concentration[layer],'\0',(num_size+2)*sizeof(float));
	(void)memset((char *)distribution[layer],'\0',(num_size+2)*sizeof(float)); 
      } /* end else */
    } /* end loop over layers */
    
    if(initial_distribution_type == RAMASWAMY){
      for(layer=1;layer<=num_layer;layer++){
	/* the radiation model needs a gradient in pressures between layers
	 because it does finite differencing over a pressure grid */ 
	env_pressure[layer]=RAMASWAMY_PRESSURE-
	  .01*layer*RAMASWAMY_PRESSURE/num_layer;
	env_temperature[layer]=RAMASWAMY_TEMPERATURE;
	if((altitude[layer] >= 10000.) && (altitude[layer] <= 11000.)){
	  saturation_ice[layer]=S_ice_cloud;
	}else{
	  saturation_ice[layer]=S_ice_ambient;
	} /* end else */ 
	if((altitude[layer] >= 10000.) && (altitude[layer] <= 11000.)){
	  for(size=1;size<=num_size;size++){
	    concentration[layer][size]=INC_init;
	    distribution[layer][size]=concentration[layer][size]/delta_mass[size];
	  } /* end loop over sizes */
	  concentration[layer][num_size]=0.;
	  concentration[layer][num_size-1]=0.;
	  distribution[layer][num_size]=
	    concentration[layer][num_size]/delta_mass[num_size];
	  distribution[layer][num_size-1]=
	    concentration[layer][num_size-1]/delta_mass[num_size-1];
	}else{
	  (void)memset((char *)concentration[layer],'\0',(num_size+2)*sizeof(float));
	  (void)memset((char *)distribution[layer],'\0',(num_size+2)*sizeof(float)); 
	} /* end else */ 
      } /* end loop over layers */

      for(size=1;size<=num_size;size++){
	capacitance[size]=equiv_radius[size];
      } /* end loop over sizes */
    } /* end if initial_distribution_type == RAMASWAMY */

    /* Use Zhang's inititialzation scheme?  */ 
    if(initial_distribution_type == ZHANG){
      for(layer=1;layer<=num_layer;layer++){
	if(altitude[layer] > 10000. && altitude[layer] < 11000.){
	  saturation_ice[layer]=S_ice_cloud;
	  for(size=1;size<=num_size;size++){
	    concentration[layer][size]=
	      ((float_foo=delta_length[size]*
		size_dist_of_IWC(crystal_length[size],env_temperature[layer]))
/*		size_dist_of_IWC(crystal_length[size],float_foo_input))*/
	       > MIN_CONC_FOR_INIT) ? float_foo/3. : 0.; 
	    distribution[layer][size]=concentration[layer][size]/delta_mass[size];
	  } /* end loop over sizes */
	}else{
	  saturation_ice[layer]=S_ice_ambient;
	  (void)memset((char *)concentration[layer],'\0',(num_size+2)*sizeof(float));
	  (void)memset((char *)distribution[layer],'\0',(num_size+2)*sizeof(float)); 
	} /* end else */
      } /* end loop over layers */

      for(size=1;size<=num_size;size++){
	capacitance[size]=equiv_radius[size];
      } /* end loop over sizes */
    } /* end if initial_distribution_type == ZHANG */

    /* Let the cloud grow from scratch */ 
    if(debug == 71){
      (void)memset((char *)concentration_ptr,'\0',(num_layer+2)*(num_size+2)*sizeof(float));
      (void)memset((char *)distribution_ptr,'\0',(num_layer+2)*(num_size+2)*sizeof(float)); 
    } /* end debug */

    /* Impose some conditions beyond the computational domain so that contour
       plots will end smoothly */
/*    bcopy(concentration[1],concentration[0],*/
/*	  (num_size+2)*sizeof(float));*/
    (void)memmove(concentration[0],concentration[1],
	  (num_size+2)*sizeof(float));
    (void)memmove(concentration[num_layer+1],concentration[num_layer],
	  (num_size+2)*sizeof(float));
    
    if (debug == 11){
      for(layer=1;layer<=num_layer;layer++){
	for(size=1;size<=num_size;size++){
	  (void)fprintf
	    (stdout,"initial concentration[%i][%i] = %15.10f per liter\n",
	     layer,size,concentration[layer][size]/LITERS_PER_CUBIC_METER);
	} /* end loop over sizes */
      } /* end loop over layers */
    } /* end debug */
    
    for(layer=1;layer<=num_layer;layer++){
      eqm_vap_ice[layer]=Eqm_Vap_Ice(env_temperature[layer]);
      eqm_vap_liquid[layer]=Eqm_Vap_Liquid(env_temperature[layer]);
      pp_vapor[layer]=saturation_ice[layer]*eqm_vap_ice[layer];
      saturation_liquid[layer]=pp_vapor[layer]/eqm_vap_liquid[layer];
      mmr_vapor[layer]=pp_vapor[layer]*epsilon_vapor/
	(env_pressure[layer]-pp_vapor[layer]); /* kg/kg */
      env_density[layer]=env_pressure[layer]/
	(env_temperature[layer]*gas_const_dry_air); 
      vapor_density[layer]=mmr_vapor[layer]*env_density[layer];
      IWC[layer]=dot_product_VEC(concentration[layer]+1,crystal_mass+1,num_size); 
      IWP[layer]=dz*IWC[layer];
      total_layer_conc[layer]=total_VEC(concentration[layer]+1,num_size);
      vapor_path[layer]=dz*env_density[layer]*mmr_vapor[layer];
      
      effective_radius[layer]=
	(max_VEC(concentration[layer]+1,num_size) > EMINUSNINE) ?
	  Effective_Radius(concentration[layer],equiv_rad_squared,
			   equiv_rad_cubed,num_size) : ONE_MICRON;
    } /* end loop over layers */
    
    for(layer=1;layer<=num_layer;layer++){
      for(size=1;size<=num_size;size++){
	fall_speed[size][layer]=
/*	    Fall_Speed_Simple(crystal_length[size],env_pressure[layer]);*/
	  Fall_Speed_Complex(crystal_length[size],
			     env_pressure[layer],
			     env_temperature[layer],
			     env_density[layer],
			     ice_crystal_habit);
      } /* end loop over sizes */
    } /* end loop over layers */
    
    /* compute the combined fall speed */ 
    (void)memmove(daltitude_dtime[0],fall_speed[0],
	  (num_size+2)*(num_layer+2)*sizeof(float));

    /* zero the wind at the boundaries */ 
    (void)memset((char *)wind_speed,'\0',(num_layer+2)*sizeof(float));
    if(UPDRAFT_PROFILE_IS_FLAT){
      for(layer=1;layer<=num_layer;layer++){
	wind_speed[layer]=wind_speed_max;
      } /* end loop over layers */
    }else{
      /* This is the algorithm used in my comps II */ 
      for(layer=1;layer<=num_layer;layer++){
	if((altitude[layer] >= altitude[1]+4.*dz) && 
	   (altitude[layer] <= altitude[num_layer]-4.*dz)){
	  float_foo=fabs(grid_base+grid_thick/2.-altitude[layer])/
	    (grid_thick/2.-4.*dz);
	  wind_speed[layer]=wind_speed_max*(1.-float_foo*float_foo);
	} /* end if */
      } /* end loop over layers */
    }/* endif */

    if(VAPOR_ADVECT){
      for(layer=1;layer<=num_layer;layer++){
	for(size=1;size<=num_size;size++){
	  daltitude_dtime[size][layer]+=wind_speed[layer];
	} /* end loop over sizes */
      } /* end loop over layers */
    }else{
      for(layer=1;layer<=num_layer;layer++){
	wind_speed[layer]=0.;
      } /* end loop over layers */
    }/* endif */
    
    /* Bootstrap variables which are derivative of user specified
       variables or otherwise unassignable initially */     
    for(layer=1;layer<=num_layer;layer++){
      vapor_diffusivity[layer]=Vapor_Diffusivity
	(env_temperature[layer],env_pressure[layer]);
      thermal_conductivity[layer]=
	Thermal_Conductivity(env_temperature[layer]);
      for(size=1;size<=num_size;size++){
	ice_density[layer][size]=ice_density_of_length(/*crystal_length[size]*/); 
      } /* end loop over sizes */
    } /* end loop over layers */
    
    Kinetic_Ventilation(conductivity_correction,
			diffusivity_correction,
			thermal_conductivity,vapor_diffusivity,
			daltitude_dtime,
			env_temperature,env_density,
			crystal_length,prism_radius,equiv_rad_squared,
			capacitance,
			num_layer,num_size);
    
    /* override time step choice if the Courant criterium for stability 
       due to the vertical advection will be violated */
    dt_net_fall_speed=alpha_transport*dz/fabs
      (max_abs_MAT_offset(daltitude_dtime,num_size,num_layer));
    dt = (dt_net_fall_speed < dt) ? dt_net_fall_speed : dt;
    if(wind_speed_max != 0.){
      if(VAPOR_ADVECT){
	dt_wind_speed=alpha_transport*dz/fabs(max_abs_VEC(wind_speed+1,num_layer));
	dt = (dt_wind_speed < dt) ? dt_wind_speed : dt;
      }else{
	dt_wind_speed = 0.;
      } /* end else */
    }else{
      dt_wind_speed = 0.;
    } /* end else */

    /* Initialize miscellaneous things here, regardless of RADIATION */
    (void)memset((char *)CCN_activated,'\0',(num_layer+2)*sizeof(float));
    (void)memset((char *)LW_flux_div,'\0',(num_layer+2)*sizeof(float));
    (void)memset((char *)SW_flux_div,'\0',(num_layer+2)*sizeof(float));
    (void)memset((char *)LW_cloud_forcing,'\0',(num_step+1)*sizeof(float));
    (void)memset((char *)LW_mass_abs_coeff,'\0',(num_layer+2)*sizeof(float));
    (void)memset((char *)SW_cloud_forcing,'\0',(num_step+1)*sizeof(float));
    (void)memset((char *)advective_heating,'\0',(num_layer+2)*sizeof(float));
    (void)memset((char *)dew_point_temp,'\0',(num_layer+2)*sizeof(float));
    (void)memset((char *)flux_down_LW,'\0',(num_layer+2)*sizeof(float));
    (void)memset((char *)flux_down_LW_int,'\0',(num_layer+2)*sizeof(float));
    (void)memset((char *)flux_down_SW,'\0',(num_layer+2)*sizeof(float));
    (void)memset((char *)flux_up_LW,'\0',(num_layer+2)*sizeof(float));
    (void)memset((char *)flux_up_LW_int,'\0',(num_layer+2)*sizeof(float));
    (void)memset((char *)flux_up_SW,'\0',(num_layer+2)*sizeof(float));
    (void)memset((char *)frost_point_temp,'\0',(num_layer+2)*sizeof(float));
    (void)memset((char *)heating_rate_LW,'\0',(num_layer+2)*sizeof(float));
    (void)memset((char *)heating_rate_SW,'\0',(num_layer+2)*sizeof(float));
    (void)memset((char *)heating_rate_net,'\0',(num_layer+2)*sizeof(float));
    (void)memset((char *)latent_heating,'\0',(num_layer+2)*sizeof(float));
    (void)memset((char *)cloud_optical_depth,'\0',(num_step+1)*sizeof(float));
    (void)memset((char *)optical_depth,'\0',(num_layer+2)*sizeof(float));
    (void)memset((char *)LW_abs_opt_depth,'\0',(num_layer+2)*sizeof(float));

    /* Initialize the Radiative and CCM2 arrays here */ 
    if(RADIATION){
      for(layer=1;layer<=num_layer;layer++){
	cloud_cover[layer]=1.;
	for(size=1;size<=num_size;size++){
	  crystal_temperature[layer][size]=env_temperature[layer];
	} /* end loop over sizes */
      } /* end loop over layers */
      
      Warren_SW_ice_init(wavelength,ice_real_idx,ice_imag_idx);
      
      /* get the SW optical efficiencies */ 
      Bohren_Mie(wavelength,num_band,
		 ice_real_idx,ice_imag_idx,
		 hex_equiv_radius,num_size,
		 size_param,
		 absorption_eff,
		 scattering_eff,
		 extinction_eff,
		 absorption_x_sec,
		 scattering_x_sec,
		 extinction_x_sec,
		 hex_equiv_rad_squared,
		 asymmetry_param);
      
      /* Decide whether to use the spherical or hex-crystal properties */ 
      if(LIOU){
	if(False){
	  Liou_interp(wavelength,num_band,
		      hex_column_length,num_size,
		      hex_equiv_rad_squared,
		      extinction_x_sec,
		      scattering_x_sec,
		      absorption_x_sec,
		      asymmetry_param);
	}else{
	  Liou_fudge(wavelength,num_band,
		     hex_column_length,num_size,
		     hex_equiv_radius,
		     hex_equiv_rad_squared,
		     absorption_x_sec,
		     scattering_x_sec,
		     extinction_x_sec,
		     asymmetry_param);
	} /* end else */
      } /* end if */ 

      /* Free all the unneeded memory */ 
      free(absorption_eff_ptr);
      free(extinction_eff_ptr);
      free(scattering_eff_ptr);
      free(size_param_ptr);

      Warren_LW_ice_init(IR_bandwidth,IR_wavelength,
			 IR_ice_real_idx,IR_ice_imag_idx);
      
      /* get the LW optical efficiencies */ 
      Bohren_Mie(IR_wavelength,num_LW_band,
		 IR_ice_real_idx,IR_ice_imag_idx,
		 hex_equiv_radius,num_size,
		 LW_size_param,
		 LW_absorption_eff,
		 LW_scattering_eff,
		 LW_extinction_eff,
		 LW_absorption_x_sec,
		 LW_scattering_x_sec,
		 LW_extinction_x_sec,
		 hex_equiv_rad_squared,
		 LW_asymmetry_param);
  
      if(True){
	/* NB: this routine only changes the cross-sections, the 
	   efficiencies should also be changed if they are used again. */ 
	Liou_IR_fudge
	  (num_LW_band,
	   num_size,
	   LW_absorption_x_sec, /* returned scaled absorption */ 
	   LW_scattering_x_sec, /* returned scaled scattering */ 
	   LW_extinction_x_sec, /* returned scaled extinction */ 
	   LW_asymmetry_param); /* returned scaled asymmetry */ 
      } /* end if */

      /* Integrate the crystal cooling rate over the longwave spectrum.
	 This routine is where planck_avg_abs_x_sec is computed */
      init_IR_params
	(LW_planck_radiance,
	 IR_bandwidth,
	 IR_wavelength,
	 LW_spec_band_weight,
	 LW_absorption_x_sec,
	 planck_avg_abs_x_sec,
	 env_temperature,
	 num_LW_band,
	 num_layer,
	 num_size);
      
      /* Free all the unneeded memory */ 
      free(LW_absorption_eff_ptr);
      free(LW_asymmetry_param_ptr);
      free(LW_extinction_eff_ptr);
      free(LW_extinction_x_sec_ptr);
      free(LW_planck_radiance_ptr);
      free(LW_scattering_eff_ptr);
      free(LW_scattering_x_sec_ptr);
      free(LW_size_param_ptr);
      
      if(debug == 62){
	for(layer=1;layer<=num_layer;layer++){
	  /* Set the LWmac to its CCM2 default in m^2/kg, it will 
	     be converted to m^2/g in the ccm2rad code itself */ 
	  LW_mass_abs_coeff[layer]=60.; 
	} /* end loop over layers */
      } /* end debug */

      init_ccm2_input
	(atmospheric_profile_type,

	 ccm2_cloud_cover,
	 ccm2_env_pressure,
	 ccm2_env_temperature,
	 ccm2_ice_path,
	 ccm2_mmr_O3,
	 ccm2_mmr_vapor,
	 env_pressure,
	 mmr_O3,
	 
	 &ASNIR,
	 &ASVIS,
	 &AWNIR,
	 &AWVIS,
	 &day_of_year,
	 &frac_strng_zen_ang_srf,
	 &ground_pressure,
	 &latitude,
	 &skin_temp,
	 &snow_cover,
	 &surf_air_temp,
	 &surf_roughness,
	 &surf_type_flag,
	 
	 &num_ccm2_level,
	 num_layer,
	 &bottom_ccm2_interface_level,
	 &top_ccm2_interface_level);
      
      /* NOTE: I don't like the way the ccm2 arrays are not dynamically 
	 allocated to the exactly correct size, now num_ccm2_level is changed
	 because of discarded levels, so each ccm2 array
	 has a little extra padding, it's too much like a FORTRAN parameter 
	 dimension procedure, but it should work */ 
      num_ccm2_level-=top_ccm2_interface_level-bottom_ccm2_interface_level-1;
      
      /* Now move the cloud into the space just vacated */ 
      (void)memmove(ccm2_cloud_cover+bottom_ccm2_interface_level+1,cloud_cover+1,
	    (num_layer+0)*sizeof(float));
    }else{
      (void)memset((char *)cloud_cover,'\0',(num_layer+2)*sizeof(float));
    } /* end else RADIATION*/
    
    /* Save some of the original arrays for plotting purposes */
    (void)memmove(orig_env_temp,env_temperature,(num_layer+2)*sizeof(float));
    (void)memmove(orig_mmr_vapor,mmr_vapor,(num_layer+2)*sizeof(float));
    
    /* Initialize time series arrays */
    cloud_base_idx=Cloud_Base_idx(IWC,num_layer,time_step,total_layer_conc);
    cloud_top_idx=Cloud_Top_idx(IWC,num_layer,time_step,total_layer_conc);
    num_cloudy_layers=cloud_top_idx-cloud_base_idx+1;

    cloud_base_of_time[time_step]=altitude[cloud_base_idx];
    cloud_top_of_time[time_step]=altitude[cloud_top_idx];

    cloud_effective_radius[time_step]=0.;
    for(layer=cloud_base_idx;layer<=cloud_top_idx;layer++){
      cloud_effective_radius[time_step]+=effective_radius[layer]*
	total_layer_conc[layer];
    } /* end loop over layers */
    cloud_effective_radius[time_step]/=((float_foo=total_VEC(total_layer_conc+cloud_base_idx,num_cloudy_layers)) > EMINUSTWELVE) ? float_foo : 1.;

    saturation_ice_of_time[time_step]=max_VEC(saturation_ice+1,num_layer);
    
    max_length_of_time[time_step]=0.;
    for(size=num_size;size>=1;size--){
      for(layer=1;layer<=num_layer;layer++){
	if(concentration[layer][size] >= EMINUSTHREE)
	  max_length_of_time[time_step]=crystal_length[size];
      } /* end loop over layers */
      if(max_length_of_time[time_step] > 0.) break;
    } /* end loop over sizes */
    
    surface_area_of_time[time_step]=0.;
    for(layer=1;layer<=num_layer;layer++){
      surface_area_of_time[time_step]+=4.*M_PI*
	dot_product_VEC(equiv_rad_squared+1,concentration[layer]+1,num_size);
    } /* end loop over layers */
    
    IWP_of_time[time_step]=dz*total_VEC(IWC+1,num_layer);
    vapor_path_of_time[time_step]=total_VEC(vapor_path+1,num_layer);
    water_path_of_time[time_step]=IWP_of_time[time_step]+
      vapor_path_of_time[time_step];
    time_array[time_step]=0.;
    
    /* Initialize the snapshots */ 
    frame=0;
    time_snapshot[frame]=time_array[time_step];
    for(layer=1;layer<=num_layer;layer++){
      IWC_snapshot[frame][layer]=IWC[layer];
      number_snapshot[frame][layer]=total_layer_conc[layer];
    } /* end loop over layers */

  } /* endif bootstrapping */
  
  /* prepare the output file */ 
  if(NETCDF_OUTPUT){

    /* Write the data to a netCDF file */ 
    /* enter define mode while opening the file */
    cdfid=nccreate(out_file, NC_CLOBBER);
    
    /* define sizes of dimensions */
    num_layerp2=(long)(num_layer+2);
    num_bandp2=(long)(num_band+2);
    num_framep1=(long)(num_frame+1);
    num_sizep2=(long)(num_size+2);
    num_CCN_sizep2=(long)(num_CCN_size+2);
    num_stepp1=(long)(num_step+1);

    altitude_dim=ncdimdef(cdfid,"altitude",num_layer+2);
/*    band_dim=ncdimdef(cdfid,"band", num_band+2);*/
/*    frame_dim=ncdimdef(cdfid,"frame",num_frame+1);*/
/*    size_dim=ncdimdef(cdfid,"size",num_size+2);*/
/*    CCN_size_dim=ncdimdef(cdfid,"CCN_size",num_CCN_size+2);*/
/*    time_dim=ncdimdef(cdfid,"time",num_step+1);*/

    num_layerp2_dim=ncdimdef(cdfid,"num_layerp2",num_layerp2);
    num_bandp2_dim=ncdimdef(cdfid,"num_bandp2", num_bandp2);
    num_framep1_dim=ncdimdef(cdfid,"num_framep1",num_framep1);
    num_sizep2_dim=ncdimdef(cdfid,"num_sizep2",num_sizep2);
    num_CCN_sizep2_dim=ncdimdef(cdfid,"num_CCN_sizep2",num_CCN_sizep2);
    num_stepp1_dim=ncdimdef(cdfid,"num_stepp1",num_stepp1);
    char_val='\000';

    /* define scalar variables */
    dt_id=ncvardef(cdfid,"dt",NC_FLOAT,0,dims);
    dz_id=ncvardef(cdfid,"dz",NC_FLOAT,0,dims);
    num_CCN_size_id=ncvardef(cdfid,"num_CCN_size",NC_LONG,0,dims);
    num_band_id=ncvardef(cdfid,"num_band",NC_LONG,0,dims);
    num_ccm2_level_id=ncvardef(cdfid,"num_ccm2_level",NC_LONG,0,dims);
    num_cloudy_layers_id=ncvardef(cdfid,"num_cloudy_layers",NC_LONG,0,dims);
    num_frame_id=ncvardef(cdfid,"num_frame",NC_LONG,0,dims);
    num_layer_id=ncvardef(cdfid,"num_layer",NC_LONG,0,dims);
    num_size_id=ncvardef(cdfid,"num_size",NC_LONG,0,dims);
    num_step_id=ncvardef(cdfid,"num_step",NC_LONG,0,dims);
    agg_step_id=ncvardef(cdfid,"agg_step",NC_LONG,0,dims);
    rad_step_id=ncvardef(cdfid,"rad_step",NC_LONG,0,dims);
    plot_step_id=ncvardef(cdfid,"plot_step",NC_LONG,0,dims);
    PAUSE_id=ncvardef(cdfid,"PAUSE",NC_LONG,0,dims);
    
    /* define two dimensional variables */
    dims[0]=num_framep1_dim;
    dims[1]=altitude_dim;
    IWC_snapshot_id=ncvardef(cdfid,"IWC_snapshot",NC_FLOAT,2,dims);
    dims[0]=num_framep1_dim;
    dims[1]=altitude_dim;
    number_snapshot_id=ncvardef(cdfid,"number_snapshot",NC_FLOAT,2,dims);
    dims[0]=altitude_dim;
    dims[1]=num_sizep2_dim;
    concentration_id=ncvardef(cdfid,"concentration",NC_FLOAT,2,dims);
    dims[0]=altitude_dim;
    dims[1]=num_sizep2_dim;
    distribution_id=ncvardef(cdfid,"distribution",NC_FLOAT,2,dims);
    dims[0]=altitude_dim;
    dims[1]=num_sizep2_dim;
    crystal_heat_SW_id=ncvardef(cdfid,"crystal_heat_SW",NC_FLOAT,2,dims);
    dims[0]=altitude_dim;
    dims[1]=num_sizep2_dim;
    crystal_heat_LW_id=ncvardef(cdfid,"crystal_heat_LW",NC_FLOAT,2,dims);
    dims[0]=altitude_dim;
    dims[1]=num_bandp2_dim;
    eff_asymmetry_param_id=ncvardef(cdfid,"eff_asymmetry_param",NC_FLOAT,2,dims);
    dims[0]=altitude_dim;
    dims[1]=num_bandp2_dim;
    eff_ss_albedo_id=ncvardef(cdfid,"eff_ss_albedo",NC_FLOAT,2,dims);
    dims[0]=altitude_dim;
    dims[1]=num_bandp2_dim;
    eff_tau_extinction_id=ncvardef(cdfid,"eff_tau_extinction",NC_FLOAT,2,dims);
    dims[0]=num_sizep2_dim;
    dims[1]=altitude_dim;
    daltitude_dtime_id=ncvardef(cdfid,"daltitude_dtime",NC_FLOAT,2,dims);
    
    /* define one dimensional variables */

    CCN_activated_id=ncvardef(cdfid,"CCN_activated",NC_FLOAT,1,&altitude_dim);
    CCN_conc_id=ncvardef(cdfid,"CCN_conc",NC_FLOAT,1,&num_CCN_sizep2_dim);
    CCN_diameter_id=ncvardef(cdfid,"CCN_diameter",NC_FLOAT,1,&num_CCN_sizep2_dim);
    CCN_mass_id=ncvardef(cdfid,"CCN_mass",NC_FLOAT,1,&num_CCN_sizep2_dim);
    aspect_ratio_id=ncvardef(cdfid,"aspect_ratio",NC_FLOAT,1,&num_sizep2_dim);
    crystal_diameter_id=ncvardef(cdfid,"crystal_diameter",NC_FLOAT,1,&num_sizep2_dim);
    saturation_liquid_id=ncvardef(cdfid,"saturation_liquid",NC_FLOAT,1,&altitude_dim);
    IWC_id=ncvardef(cdfid,"IWC",NC_FLOAT,1,&altitude_dim);
    IWP_id=ncvardef(cdfid,"IWP",NC_FLOAT,1,&altitude_dim);
    IWP_of_time_id=ncvardef(cdfid,"IWP_of_time",NC_FLOAT,1,&num_stepp1_dim);
    LW_cloud_forcing_id=ncvardef(cdfid,"LW_cloud_forcing",NC_FLOAT,1,&num_stepp1_dim);
    SW_cloud_forcing_id=ncvardef(cdfid,"SW_cloud_forcing",NC_FLOAT,1,&num_stepp1_dim);
    advective_heating_id=ncvardef(cdfid,"advective_heating",NC_FLOAT,1,&altitude_dim);
    albedo_of_time_id=ncvardef(cdfid,"albedo_of_time",NC_FLOAT,1,&num_stepp1_dim);
    altitude_id=ncvardef(cdfid,"altitude",NC_FLOAT,1,&altitude_dim);
    cloud_base_of_time_id=ncvardef(cdfid,"cloud_base_of_time",NC_FLOAT,1,&num_stepp1_dim);
    cloud_effective_radius_id=ncvardef(cdfid,"cloud_effective_radius",NC_FLOAT,1,&num_stepp1_dim);
    cloud_top_of_time_id=ncvardef(cdfid,"cloud_top_of_time",NC_FLOAT,1,&num_stepp1_dim);
    crystal_length_id=ncvardef(cdfid,"crystal_length",NC_FLOAT,1,&num_sizep2_dim);
    crystal_mass_id=ncvardef(cdfid,"crystal_mass",NC_FLOAT,1,&num_sizep2_dim);
    delta_length_id=ncvardef(cdfid,"delta_length",NC_FLOAT,1,&num_sizep2_dim);
    delta_mass_id=ncvardef(cdfid,"delta_mass",NC_FLOAT,1,&num_sizep2_dim);
    effective_radius_id=ncvardef(cdfid,"effective_radius",NC_FLOAT,1,&altitude_dim);
    emissivity_of_time_id=ncvardef(cdfid,"emissivity_of_time",NC_FLOAT,1,&num_stepp1_dim);
    env_density_id=ncvardef(cdfid,"env_density",NC_FLOAT,1,&altitude_dim);
    env_pressure_id=ncvardef(cdfid,"env_pressure",NC_FLOAT,1,&altitude_dim);
    env_temperature_id=ncvardef(cdfid,"env_temperature",NC_FLOAT,1,&altitude_dim);
    flux_down_LW_id=ncvardef(cdfid,"flux_down_LW",NC_FLOAT,1,&altitude_dim);
    flux_down_SW_id=ncvardef(cdfid,"flux_down_SW",NC_FLOAT,1,&altitude_dim);
    flux_up_LW_id=ncvardef(cdfid,"flux_up_LW",NC_FLOAT,1,&altitude_dim);
    flux_up_SW_id=ncvardef(cdfid,"flux_up_SW",NC_FLOAT,1,&altitude_dim);
    heating_rate_LW_id=ncvardef(cdfid,"heating_rate_LW",NC_FLOAT,1,&altitude_dim);
    heating_rate_SW_id=ncvardef(cdfid,"heating_rate_SW",NC_FLOAT,1,&altitude_dim);
    heating_rate_net_id=ncvardef(cdfid,"heating_rate_net",NC_FLOAT,1,&altitude_dim);
    latent_heating_id=ncvardef(cdfid,"latent_heating",NC_FLOAT,1,&altitude_dim);
    max_length_of_time_id=ncvardef(cdfid,"max_length_of_time",NC_FLOAT,1,&num_stepp1_dim);
    mmr_vapor_id=ncvardef(cdfid,"mmr_vapor",NC_FLOAT,1,&altitude_dim);
    cloud_optical_depth_id=ncvardef(cdfid,"cloud_optical_depth",NC_FLOAT,1,&num_stepp1_dim);
    orig_env_temp_id=ncvardef(cdfid,"orig_env_temp",NC_FLOAT,1,&altitude_dim);
    orig_mmr_vapor_id=ncvardef(cdfid,"orig_mmr_vapor",NC_FLOAT,1,&altitude_dim);
    potential_temp_id=ncvardef(cdfid,"potential_temp",NC_FLOAT,1,&altitude_dim);
    pp_vapor_id=ncvardef(cdfid,"pp_vapor",NC_FLOAT,1,&altitude_dim);
    saturation_ice_of_time_id=ncvardef(cdfid,"saturation_ice_of_time",NC_FLOAT,1,&num_stepp1_dim);
    saturation_ice_id=ncvardef(cdfid,"saturation_ice",NC_FLOAT,1,&altitude_dim);
    surface_area_of_time_id=ncvardef(cdfid,"surface_area_of_time",NC_FLOAT,1,&num_stepp1_dim);
    time_array_id=ncvardef(cdfid,"time_array",NC_FLOAT,1,&num_stepp1_dim);
    time_snapshot_id=ncvardef(cdfid,"time_snapshot",NC_FLOAT,1,&num_framep1_dim);
    vapor_density_id=ncvardef(cdfid,"vapor_density",NC_FLOAT,1,&altitude_dim);
    vapor_path_id=ncvardef(cdfid,"vapor_path",NC_FLOAT,1,&altitude_dim);
    vapor_path_of_time_id=ncvardef(cdfid,"vapor_path_of_time",NC_FLOAT,1,&num_stepp1_dim);
    optical_depth_id=ncvardef(cdfid,"optical_depth",NC_FLOAT,1,&altitude_dim);
    water_path_of_time_id=ncvardef(cdfid,"water_path_of_time",NC_FLOAT,1,&num_stepp1_dim);
    wind_speed_id=ncvardef(cdfid,"wind_speed",NC_FLOAT,1,&altitude_dim);
    
    /* assign attributes */
    ncattput(cdfid,CCN_activated_id,"units",NC_CHAR,7,(void *)"meter-3");
    ncattput(cdfid,CCN_conc_id,"units",NC_CHAR,7,(void *)"meter-3");
    ncattput(cdfid,CCN_diameter_id,"units",NC_CHAR,5,(void *)"meter");
    ncattput(cdfid,CCN_mass_id,"units",NC_CHAR,8,(void *)"kilogram");
    ncattput(cdfid,aspect_ratio_id,"units",NC_CHAR,5,(void *)"ratio");
    ncattput(cdfid,crystal_diameter_id,"units",NC_CHAR,5,(void *)"meter");
    ncattput(cdfid,saturation_liquid_id,"units",NC_CHAR,5,(void *) "ratio");
    ncattput(cdfid,IWC_snapshot_id,"units",NC_CHAR,16,(void *)"kilogram meter-3");
    ncattput(cdfid,number_snapshot_id,"units",NC_CHAR,7,(void *)"meter-3");
    ncattput(cdfid,concentration_id,"units",NC_CHAR,7,(void *)"meter-3");
    ncattput(cdfid,distribution_id,"units",NC_CHAR,15,(void *)"meter-3 kilogram-1");
    ncattput(cdfid,crystal_heat_SW_id,"units",NC_CHAR,4,(void *)"watt");
    ncattput(cdfid,crystal_heat_LW_id,"units",NC_CHAR,4,(void *)"watt");
    ncattput(cdfid,eff_asymmetry_param_id,"units",NC_CHAR,0,(void *) &char_val);
    ncattput(cdfid,eff_ss_albedo_id,"units",NC_CHAR,0,(void *) &char_val);
    ncattput(cdfid,eff_tau_extinction_id,"units",NC_CHAR,0,(void *) &char_val);
    ncattput(cdfid,IWC_id,"units",NC_CHAR,16,(void *)"kilogram meter-3");
    ncattput(cdfid,IWP_id,"units",NC_CHAR,16,(void *)"kilogram meter-2");
    ncattput(cdfid,IWP_of_time_id,"units",NC_CHAR,16,(void *)"kilogram meter-2");
    ncattput(cdfid,LW_cloud_forcing_id,"units",NC_CHAR,0,(void *) &char_val);
    ncattput(cdfid,SW_cloud_forcing_id,"units",NC_CHAR,0,(void *) &char_val);
    ncattput(cdfid,advective_heating_id,"units",NC_CHAR,15,(void *)"kelvin second-1");
    ncattput(cdfid,albedo_of_time_id,"units",NC_CHAR,0,(void *) &char_val);
    ncattput(cdfid,altitude_id,"units",NC_CHAR,5,(void *)"meter");
    ncattput(cdfid,cloud_base_of_time_id,"units",NC_CHAR,0,(void *) &char_val);
    ncattput(cdfid,cloud_effective_radius_id,"units",NC_CHAR,0,(void *) &char_val);
    ncattput(cdfid,cloud_top_of_time_id,"units",NC_CHAR,0,(void *) &char_val);
    ncattput(cdfid,crystal_length_id,"units",NC_CHAR,5,(void *)"meter");
    ncattput(cdfid,crystal_mass_id,"units",NC_CHAR,8,(void *)"kilogram");
    ncattput(cdfid,delta_length_id,"units",NC_CHAR,5,(void *)"meter");
    ncattput(cdfid,delta_mass_id,"units",NC_CHAR,8,(void *)"kilogram");
    ncattput(cdfid,effective_radius_id,"units",NC_CHAR,5,(void *)"meter");
    ncattput(cdfid,emissivity_of_time_id,"units",NC_CHAR,0,(void *) &char_val);
    ncattput(cdfid,env_density_id,"units",NC_CHAR,16,(void *)"kilogram meter-3");
    ncattput(cdfid,env_pressure_id,"units",NC_CHAR,6,(void *)"pascal");
    ncattput(cdfid,env_temperature_id,"units",NC_CHAR,6,(void *)"kelvin");
    ncattput(cdfid,flux_down_LW_id,"units",NC_CHAR,12,(void *)"watt meter-2");
    ncattput(cdfid,flux_down_SW_id,"units",NC_CHAR,12,(void *)"watt meter-2");
    ncattput(cdfid,flux_up_LW_id,"units",NC_CHAR,12,(void *)"watt meter-2");
    ncattput(cdfid,flux_up_SW_id,"units",NC_CHAR,12,(void *)"watt meter-2");
    ncattput(cdfid,heating_rate_LW_id,"units",NC_CHAR,15,(void *)"kelvin second-1");
    ncattput(cdfid,heating_rate_SW_id,"units",NC_CHAR,15,(void *)"kelvin second-1");
    ncattput(cdfid,heating_rate_net_id,"units",NC_CHAR,15,(void *)"kelvin second-1");
    ncattput(cdfid,latent_heating_id,"units",NC_CHAR,15,(void *)"kelvin second-1");
    ncattput(cdfid,max_length_of_time_id,"units",NC_CHAR,5,(void *)"meter");
    ncattput(cdfid,mmr_vapor_id,"units",NC_CHAR,19,(void *)"kilogram kilogram-1");
    ncattput(cdfid,cloud_optical_depth_id,"units",NC_CHAR,0,(void *) &char_val);
    ncattput(cdfid,orig_env_temp_id,"units",NC_CHAR,6,(void *)"kelvin");
    ncattput(cdfid,orig_mmr_vapor_id,"units",NC_CHAR,19,(void *)"kilogram kilogram-1");
    ncattput(cdfid,potential_temp_id,"units",NC_CHAR,6,(void *)"kelvin");
    ncattput(cdfid,pp_vapor_id,"units",NC_CHAR,6,(void *)"pascal");
    ncattput(cdfid,saturation_ice_of_time_id,"units",NC_CHAR,0,(void *) &char_val);
    ncattput(cdfid,saturation_ice_id,"units",NC_CHAR,0,(void *) &char_val);
    ncattput(cdfid,surface_area_of_time_id,"units",NC_CHAR,6,(void *)"meter2");
    ncattput(cdfid,time_array_id,"units",NC_CHAR,6,(void *)"second");
    ncattput(cdfid,time_snapshot_id,"units",NC_CHAR,6,(void *)"second");
    ncattput(cdfid,vapor_density_id,"units",NC_CHAR,16,(void *)"kilogram meter-3");
    ncattput(cdfid,vapor_path_id,"units",NC_CHAR,16,(void *)"kilogram meter-2");
    ncattput(cdfid,vapor_path_of_time_id,"units",NC_CHAR,16,(void *)"kilogram meter-2");
    ncattput(cdfid,optical_depth_id,"units",NC_CHAR,0,(void *) &char_val);
    ncattput(cdfid,water_path_of_time_id,"units",NC_CHAR,16,(void *)"kilogram meter-2");
    ncattput(cdfid,wind_speed_id,"units",NC_CHAR,14,(void *)"meter second-1");
    
    ncattput(cdfid,NC_GLOBAL,"cmdline",NC_CHAR,strlen(cmdline),(void *)cmdline);

    /* leave define mode */
    ncendef (cdfid); 
    
  } /* endif(NETCDF_OUTPUT) */ 

  if(MOVIE_NETCDF_OUTPUT){
    /* re-enter define mode */
    ncredef (cdfid); 

    /* define three dimensional variables */
    dims[0]=num_framep1_dim;
    dims[1]=altitude_dim;
    dims[2]=num_sizep2_dim;
    movie_concentration_id=ncvardef(cdfid,"movie_concentration",NC_FLOAT,3,dims);

    /* define two dimensional variables */
    dims[0]=num_framep1_dim;
    dims[1]=altitude_dim;
    movie_effective_radius_id=ncvardef(cdfid,"movie_effective_radius",NC_FLOAT,2,dims);

    /* define one dimensional variables */
    movie_cloud_effective_radius_id=ncvardef(cdfid,"movie_cloud_effective_radius",NC_FLOAT,1,&num_framep1_dim);

    /* assign attributes */
    ncattput(cdfid,movie_concentration_id,"units",NC_CHAR,7,(void *)"meter-3");
    ncattput(cdfid,movie_effective_radius_id,"units",NC_CHAR,5,(void *)"meter");
    ncattput(cdfid,movie_cloud_effective_radius_id,"units",NC_CHAR,5,(void *)"meter");

    /* leave define mode */
    ncendef (cdfid); 

    /* output the initial conditions frame */ 
    /* append a scalar to a one-dimensional array */
    start[0]=(long)frame;
    ncvarput1(cdfid,movie_cloud_effective_radius_id,start,(void *)(cloud_effective_radius+frame));
    
    /* append a one-dimensional array to a two-dimensional array */
    start[0]=(long)frame;
    start[1]=0L;
    
    count[0]=1L;
    count[1]=num_layerp2;
    ncvarput(cdfid,movie_effective_radius_id,start,count,(void *)effective_radius);
    
    /* append a two dimensional array to a three-dimensional array */
    start[0]=(long)frame;
    start[1]=0L;
    start[2]=0L;
    
    count[0]=1L;
    count[1]=num_layerp2;
    count[2]=num_sizep2;
    ncvarput(cdfid,movie_concentration_id,start,count,(void *)concentration_ptr);
  } /* endif(MOVIE_NETCDF_OUTPUT) */ 

  /* Now must compute the simultaneous solution of the energy-balance
     and of the mass-continuity equations */
  (void)fprintf(stdout,"\nCirrus Cloud Evolution Model %s\n",rcs_Revision);
  (void)fprintf(stdout,"Debugging level: %i\n",debug);
  (void)fprintf(stdout,"%s, %s, %s\n",
		atmospheric_profile_string_array[atmospheric_profile_type],
		initial_distribution_string_array[initial_distribution_type],
		ice_crystal_habit_string_array[ice_crystal_habit]);
  if(!AGGREGATION)
    (void)fprintf(stdout,"Not performing aggregation phase.\n");
  if(!EULER)
    (void)fprintf(stdout,"Not performing Eulerian growth phase.\n");
  if(!ICE_ADVECT)
    (void)fprintf(stdout,"Not performing net fall velocity advection phase.\n");
  if(!LAGRANGE)
    (void)fprintf(stdout,"Not performing Lagrangian growth phase.\n");
  if(!VAPOR_ADVECT)
    (void)fprintf(stdout,"Not performing water vapor wind advection phase.\n");
  if(UPDRAFT_PROFILE_IS_FLAT)
    (void)fprintf(stdout,"Not using parabolic updraft profile, using constant profile instead.\n");
  if(!HETEROGENEOUS_NUCLEATION)
    (void)fprintf(stdout,"Not performing MDC92-style heterogeneous nucleation.\n");
  if(!HOMOGENEOUS_NUCLEATION)
    (void)fprintf(stdout,"Not performing DMC94-style homogeneous nucleation.\n");
  if(!RADIATION)
    (void)fprintf(stdout,"Not performing Radiation calculations.\n");
  if(!RAD_FEEDBACK)
    (void)fprintf(stdout,"Not feeding back flux convergences into env. temp.\n");
  if(!LIOU)
    (void)fprintf(stdout,"Not radiatively treating crystals as hex columns.\n");
  if(!VERBOSE)
    (void)fprintf(stdout,"Not printing out WARNING messages.\n");
  if(!WEED_NEGATIVE_VALUES)
    (void)fprintf(stdout,"Not weeding out negative values after ETBFCT.\n");
  (void)fprintf(stdout,"num_size = %i and num_layer = %i\n",num_size,num_layer);
  (void)fprintf(stdout,"dz = %.1f meters\n",dz);
  (void)fprintf(stdout,"Integration period = %i seconds = %.2f minutes\n",
		(int)(num_step*dt),num_step*dt/60.);
  (void)fprintf(stdout,"%i frames will be saved (every %i steps = %.2f minutes).\n",
		num_frame+1,plot_step,plot_step*dt/60.);
  if(RADIATION)
    (void)fprintf(stdout,"Radiation routine updates every %i steps = %.2f minutes.\n",rad_step,rad_step*dt/60.);
  if(AGGREGATION)
    (void)fprintf(stdout,"Aggregation routine updates every %i steps = %.2f minutes.\n",agg_step,agg_step*dt/60.);
  (void)fprintf(stdout,"wind speed max = %.1f cm/s, ",wind_speed_max*CENTIMETERS_PER_METER);
  (void)fprintf(stdout,"maximum net fall speed = %g m/s\n",max_abs_MAT_offset(daltitude_dtime,num_size,num_layer));
  (void)fprintf(stdout,"Initial stepsize = %g seconds, net fall speed Courant = %g s, wind speed Courant = %g s\n",dt,dt_net_fall_speed,dt_wind_speed);
  (void)fprintf(stdout,"printing one \".\" per time step:\n");
  
  for(time_step=1;time_step<=num_step;time_step++){ /* begin outer time loop */
    
    if(debug == 19){
      (void)fprintf(stdout,"Begin time_step = %i. Concentration of crystals(#/m^3): size bin (across) vs layer (down)\n",time_step);
      for(size=1;size<=num_size;size++){
	(void)fprintf(stdout,"\t%i",size);
      } /* end loop over sizes */
      for(layer=1;layer<=num_layer;layer++){
	(void)fprintf(stdout,"\n%i\t",layer);
	for(size=1;size<=num_size;size++){
	  (void)fprintf(stdout,"%g\t",concentration[layer][size]);
	} /* end loop over sizes */
      } /* end loop over layers */
      (void)fprintf(stdout,"\n");
    } /* end debug */
    
    if(RADIATION && ((time_step-1)%rad_step == 0)){
      /* Now move the cloud into the space just vacated */ 
      (void)memmove(
	    ccm2_env_pressure+bottom_ccm2_interface_level+1,env_pressure+1,
	    (num_layer+0)*sizeof(float));
      (void)memmove(
	    ccm2_env_temperature+bottom_ccm2_interface_level+1,env_temperature+1,
	    (num_layer+0)*sizeof(float));
      (void)memmove(
	    ccm2_ice_path+bottom_ccm2_interface_level+1,IWP+1,
	    (num_layer+0)*sizeof(float));
      (void)memmove(
	    ccm2_mmr_vapor+bottom_ccm2_interface_level+1,mmr_vapor+1,
	    (num_layer+0)*sizeof(float));
      (void)memmove(
	    ccm2_mmr_O3+bottom_ccm2_interface_level+1,mmr_O3+1,
	    (num_layer+0)*sizeof(float));

      /* Now eliminate trace gases between cloud top and cloud bottom */
      /* commented out pending further review */ 
/*      for(layer=1;layer<=num_cloudy_layers;layer++){*/
/*	ccm2_mmr_vapor[bottom_ccm2_interface_level+cloud_base_idx+layer-1]=*/
/*	  EMINUSFIFTEEN;*/
/*	ccm2_mmr_O3[bottom_ccm2_interface_level+cloud_base_idx+layer-1]=*/
/*	  EMINUSFIFTEEN;*/
/*      }*/ /* end loop over layers */
/*      (void)memset((char *)ccm2_mmr_O3+bottom_ccm2_interface_level+cloud_base_idx,'\0',*/
/*	    (num_cloudy_layers+0)*sizeof(float));*/

      integrate_SW_optical_properties
	(asymmetry_param,
	 concentration,
	 dz,
	 eff_asymmetry_param,
	 eff_ss_albedo,
	 eff_tau_extinction,
	 eff_tau_scatter,
	 extinction_x_sec,
	 num_band,
	 num_layer,
	 num_size,
	 scattering_x_sec);
      
      if(debug == 43){
	(void)fprintf(stdout,"cloud_base_idx = %i = %.2f, cloud_top_idx = %i = %.2f\n",
		      cloud_base_idx,altitude[cloud_base_idx],
		      cloud_top_idx,altitude[cloud_top_idx]);
	for(layer=1;layer<=num_layer;layer++){
	  (void)fprintf(stdout,"\nzbin\talt(km)\twband\tlambda\tTau,s,e\tTau,e,e\tw,eff\tg,eff\n");
	  for(band=1;band<=num_band;band++){
	    (void)fprintf(stdout,"%i\t%.2f\t%i\t%.3f\t%7.1e\t%7.1e\t%.3f\t%.3f\n",
			  layer,
			  altitude[layer]*KILOMETERS_PER_METER,
			  band,
			  wavelength[band]*MICRONS_PER_METER,
			  eff_tau_scatter[layer][band],
			  eff_tau_extinction[layer][band],
			  eff_ss_albedo[layer][band],
			  eff_asymmetry_param[layer][band]);
	  } /* end loop over bands */
	} /* end loop over layers */
      } /* end debug */
      
      if(debug != 62){
	integrate_LW_mass_abs_coeff
	  (LW_absorption_x_sec,
	   LW_spec_band_weight,
	   LW_spec_vol_abs_coeff,
	   LW_mass_abs_coeff,
	   concentration,
	   IWC,
	   num_LW_band,
	   num_layer,
	   num_size);
      } /* end debug */

      if(debug == 44){
	(void)fprintf(stdout,"cloud_base_idx = %i = %.2f, cloud_top_idx = %i = %.2f, num_cloudy_layers = %i\n",
		      cloud_base_idx,altitude[cloud_base_idx],
		      cloud_top_idx,altitude[cloud_top_idx],
		      num_cloudy_layers);
	(void)fprintf(stdout,"bccmil = %i, tccmil = %i, nccmlev = %i\n",
		      bottom_ccm2_interface_level,top_ccm2_interface_level,
		      num_ccm2_level);
      } /* end debug */
      

      FORTRAN_colmod
	(&day_of_year,  /* dayyr */
	 &latitude,  /* rlat */
	 ccm2_env_pressure+1,  /* pmid */
	 ccm2_env_temperature+1,  /* t */
	 ccm2_mmr_vapor+1,  /* h2ommr */
	 ccm2_mmr_O3+1,  /* o3mmr */
	 ccm2_cloud_cover+1,  /* cldfrc */
	 ccm2_ice_path+1,  /* clwp */
	 &ground_pressure,  /* ps */
	 &surf_air_temp,  /* ts */
	 &skin_temp,  /* tg */
	 &surf_type_flag,  /* oro */
	 &surf_roughness,  /* rghnss */
	 &snow_cover,  /* sndpth */
	 &ASVIS,  /* albvss */
	 &AWVIS,  /* albvsw */
	 &ASNIR,  /* albnis */
	 &AWNIR,  /* albniw */
	 &frac_strng_zen_ang_srf,  /* frctst */ 
	 /* the above are all needed to drive the old colmod, 
	    everything below is used to implement my cloud 
	    parameterization */ 
	 &debug,  /* idebug */ 
	 &ERROR_FLAG,
	 env_pressure_int, /* envpint */ 
	 &bottom_ccm2_interface_level,  /* bccmil */ 
	 &top_ccm2_interface_level,  /* tccmil */ 
	 &num_layer,  /* numlay */ 
	 &num_ccm2_level,  /* ncclev */ 
	 &cloud_base_idx, /* cbi */ 
	 &cloud_top_idx, /* cti */ 
	 &num_cloudy_layers, /* ncldlay */ 
	 eff_asymmetry_param_ptr,  /* geff */ 
	 eff_ss_albedo_ptr,  /* aeff */ 
	 eff_tau_extinction_ptr,  /* taueff */ 
	 effective_radius, /* crecz */ 
	 SW_spec_flux_div_ptr,  /* swfdiv */ 
	 LW_flux_div,  /* lwfdiv */ 
	 flux_up_SW, /* fupsw */ 
	 flux_down_SW, /* fdsw */ 
	 flux_up_LW, /* fuplw */ 
	 flux_down_LW, /* fdlw */ 
	 &flux_up_LW_cloud_base_int, /* fupcbi */ 
	 &flux_down_LW_cloud_top_int, /* fdncti */ 
	 heating_rate_SW, /* qrscld */ 
	 heating_rate_LW, /* qrlcld */ 
	 heating_rate_net, /* qrcld */ 
	 LW_cloud_forcing+time_step-1, /* lwcf */ 
	 SW_cloud_forcing+time_step-1, /* swcf */ 
	 albedo_of_time+time_step-1, /* albedo */ 
	 &SW_surf_cloud_forcing, /* swscf */ 
	 LW_mass_abs_coeff); /* lwmac */ 
      
      if(ERROR_FLAG){
	(void)fprintf(stdout,"Unrecoverable error in column model, probably need to reset PLEV parameter and recompile. Exiting....\n");
	Exit_gracefully();
	exit(1);
      } /* end if */
      
      for(layer=cloud_base_idx;layer<=cloud_top_idx;layer++){
	float_foo=LW_DIFFUSIVITY*LW_mass_abs_coeff[layer]*IWP[layer];
	transmissivity[layer]=exp(-float_foo);
	emissivity[layer]=1.-transmissivity[layer];
	blackbody_flux[layer]=
	  stefan_boltzmann_constant*pow(env_temperature[layer],4.);
      } /* end loop over layers */

      /* Bootstrap the LW interface fluxes using the interface fluxes
       at top and bottom */ 
      flux_up_LW_int[cloud_base_idx-1]=flux_up_LW_cloud_base_int;
      flux_down_LW_int[cloud_top_idx]=flux_down_LW_cloud_top_int;
      for(layer=cloud_base_idx;layer<=cloud_top_idx;layer++){
	flux_up_LW_int[layer]=emissivity[layer]*blackbody_flux[layer]+
	  transmissivity[layer]*flux_up_LW_int[layer-1];
      } /* end loop over layers */
      for(layer=cloud_top_idx;layer>=cloud_base_idx;layer--){
	flux_down_LW_int[layer-1]=emissivity[layer]*blackbody_flux[layer]+
	  transmissivity[layer]*flux_down_LW_int[layer];
      } /* end loop over layers */

      /* Recall the quantities returned from CCM2RAD are still flux
	 differences, now take their gradients.  Remember that although
	 they are labelled as divergences (of the radiation field) they
	 are really convergences (on the ice crystals), i.e., a positive
	 flux divergence leads to a crystal heating. */ 
      for(layer=1;layer<=num_layer;layer++){
	for(band=1;band<=num_band;band++){
	  SW_spec_flux_div[layer][band]/=dz;
	} /* end loop over bands */
	LW_flux_div[layer]/=dz;
	SW_flux_div[layer]=total_VEC(SW_spec_flux_div[layer]+1,num_band);
      } /* end loop over layers */

      /* Use the layer interface fluxes to find the layer fluxes and 
	 flux divergences.  These intra-cloud divergences replace the
	 heating calculated by CCM2RAD within the cloud. */ 
      for(layer=cloud_base_idx;layer<=cloud_top_idx;layer++){
	LW_flux_div[layer]=
	  (flux_up_LW_int[layer]-flux_down_LW_int[layer]-
	   flux_up_LW_int[layer-1]+flux_down_LW_int[layer-1])/dz;
      } /* end loop over layers */

      /* The following lines correct for the artificially large flux
	 divergences at the interface between the high resolution cloud grid
	 and the low resolution ccm2 model */ 
      LW_flux_div[1]=LW_flux_div[2];
      LW_flux_div[num_layer]=LW_flux_div[num_layer-1];
      SW_flux_div[1]=SW_flux_div[2];
      SW_flux_div[num_layer]=SW_flux_div[num_layer-1];
      for(band=1;band<=num_band;band++){
	SW_spec_flux_div[1][band]=SW_spec_flux_div[2][band];
	SW_spec_flux_div[num_layer][band]=SW_spec_flux_div[num_layer-1][band];
      } /* end loop over bands */

      /* Compute the bulk Radiative parameters */ 
      float_foo3=0.;
      for(layer=cloud_base_idx;layer<=cloud_top_idx;layer++){
	/* accumulate the sigma's: this counts all crystals in the cloud */ 
	float_foo3+=LW_mass_abs_coeff[layer]*IWC[layer];
      } /* end loop over layers */;
      /* normalize the total sigma to represent crystals in a cubic meter */ 
      float_foo3/=num_cloudy_layers;
      float_foo=LW_DIFFUSIVITY*float_foo3*dz*num_cloudy_layers;
      float_foo2=1.-exp(-float_foo);
      emissivity_of_time[time_step-1]=float_foo2;
/*      albedo_of_time[time_step-1]=*/
/*	flux_up_SW[cloud_top_idx]/flux_down_SW[cloud_top_idx];*/

      /* Apportion the flux divergence onto crystals weighted by their 
	 concentration and absorption efficiency. */
      apportion_crystal_heating
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
	 num_cloudy_layers);

      if(debug_value == 1){
	for(layer=1;layer<=num_layer;layer++){
	  if(total_layer_conc[layer] > 0.){
	    float_foo=-SW_flux_div[layer]/total_layer_conc[layer];
	    float_foo2=-LW_flux_div[layer]/total_layer_conc[layer];
	    float_foo3=float_foo+float_foo2;
	  }else{
	    float_foo=0.;
	    float_foo2=0.;
	    float_foo3=0.;
	  } /* end else */ 
	  for(size=1;size<=num_size;size++){
	    crystal_heat_SW[layer][size]=float_foo;
	    crystal_heat_LW[layer][size]=float_foo2;
	    crystal_heat_net[layer][size]=float_foo3;
	  } /* end loop over sizes */
	} /* end loop over layers */
      } /* end debug */

      if(debug == 65){
	(void)fprintf(stdout,"zbin\talt\tsbin\tlength\tConc\tRSW q\tRLW q\tRTot q\tZSW q\tZLW q\tZNet q\n");
	(void)fprintf(stdout,"\tkm\t\tmicrons\t#/l\tuW\tuW\tuW\tuW\tuW\tuW\n");
	for(layer=cloud_base_idx;layer<=cloud_top_idx;layer++){
	  if(total_layer_conc[layer] > 0.){
	    float_foo=-SW_flux_div[layer]/total_layer_conc[layer];
	    float_foo2=-LW_flux_div[layer]/total_layer_conc[layer];
	    float_foo3=float_foo+float_foo2;
	  }else{
	    float_foo=0.;
	    float_foo2=0.;
	    float_foo3=0.;
	  } /* end else */ 
	  for(size=1;size<=num_size;size++){
	    (void)fprintf(stdout,"%i\t%.2f\t%i\t%.2f\t%.1f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",
			  layer,
			  altitude[layer]*KILOMETERS_PER_METER,
			  size,
			  crystal_length[size]*MICRONS_PER_METER,
			  concentration[layer][size]*CUBIC_METERS_PER_LITER,
			  float_foo*MICROWATTS_PER_WATT,
			  float_foo2*MICROWATTS_PER_WATT,
			  float_foo3*MICROWATTS_PER_WATT,
			  crystal_heat_SW[layer][size]*MICROWATTS_PER_WATT,
			  crystal_heat_LW[layer][size]*MICROWATTS_PER_WATT,
			  crystal_heat_net[layer][size]*MICROWATTS_PER_WATT);
	  } /* end loop over sizes */
	} /* end loop over layers */
      } /* end debug */
      
      if(debug == 67){
	(void)fprintf(stdout,"zbin\talt\tN*ZSW q\tdFSW/dz\tN*ZLW q\tdFLW/dz\n");
	(void)fprintf(stdout,"\tkm\tmW/m^3\tmW/m^3\tmW/m^3\tmW/m^3\n");
	for(layer=1;layer<=num_layer;layer++){
	  float_foo=dot_product_VEC(concentration[layer]+1,
				    crystal_heat_SW[layer]+1,num_size);
	  float_foo2=dot_product_VEC(concentration[layer]+1,
				     crystal_heat_LW[layer]+1,num_size);
	  (void)fprintf(stdout,"%i\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",
			layer,
			altitude[layer]*KILOMETERS_PER_METER,
			float_foo*MILLIWATTS_PER_WATT,
			SW_flux_div[layer]*MILLIWATTS_PER_WATT,
			float_foo2*MILLIWATTS_PER_WATT,
			LW_flux_div[layer]*MILLIWATTS_PER_WATT);
	} /* end loop over layers */
      } /* end debug */
      
    } /* end if RADIATION */

    if(debug == 63){
      (void)fprintf(stdout,"zbin\talt\tsbin\tlength\tConc\talpha\tbeta\tLatent\tSW Rad\tLW Rad\tNet Rad\tLat+Rad\n");
      (void)fprintf(stdout,"\tkm\t\tmicrons\t#/l\t\t\tuW\tuW\tuW\tuW\tuW\n");
    } /* end debug */

    /* Now compute layer dependent vectors in the main loop */
    for(layer=1;layer<=num_layer;layer++){
      
      /* Now compute layer,size dependent vectors in this layer loop */
      for(size=1;size<=num_size;size++){
	
	denom_1=env_temperature[layer]*gas_const_vapor/
	  (eqm_vap_ice[layer]*vapor_diffusivity[layer]*
	   diffusivity_correction[layer][size]);
	curvature=2.*surf_free_energy_ice/
	  (ice_density[layer][size]*gas_const_vapor*
  	   prism_radius[size]*env_temperature[layer]);
	numerator_1=4.*M_PI*capacitance[size]*
	  (saturation_ice[layer]-1.-curvature);
	denom_3=gas_const_vapor*env_temperature[layer]*env_temperature[layer]*
	  thermal_conductivity[layer]*conductivity_correction[layer][size];
	denom_2=latent_heat_sub*latent_heat_sub*(curvature+1.)/denom_3;
	
	if(RADIATION){
	  numerator_2=latent_heat_sub*(curvature+1.)*
	    crystal_heat_net[layer][size]/denom_3;
	}else{
	  numerator_2=0.;
	} /* end else */
	
	/* Compute the rate of change of mass per particle */
	dmass_dtime[layer][size]=(numerator_1-numerator_2)/(denom_1+denom_2);
	lagrange_mass[layer][size]=crystal_mass[size]+dmass_dtime[layer][size]*dt;
	
	if(debug == 63){
	  if(RADIATION){
	    float_foo=numerator_1/numerator_2;
	    float_foo2=numerator_1*latent_heat_sub/
	      (crystal_heat_net[layer][size]*denom_1);
	  }else{
	    float_foo=0.;
	    float_foo2=0.;
	  } /* end else */

	  (void)fprintf(stdout,"%i\t%.2f\t%i\t%.2f\t%.1f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",
			layer,
			altitude[layer]*KILOMETERS_PER_METER,
			size,
			crystal_length[size]*MICRONS_PER_METER,
			concentration[layer][size]*CUBIC_METERS_PER_LITER,
			float_foo,
			float_foo2,
			latent_heat_sub*dmass_dtime[layer][size]*
			MICROWATTS_PER_WATT,
			crystal_heat_SW[layer][size]*MICROWATTS_PER_WATT,
			crystal_heat_LW[layer][size]*MICROWATTS_PER_WATT,
			crystal_heat_net[layer][size]*MICROWATTS_PER_WATT,
			(latent_heat_sub*dmass_dtime[layer][size]+
			 crystal_heat_net[layer][size])*MICROWATTS_PER_WATT);
      } /* end debug */
      
    } /* end loop over sizes */
  } /* end loop over layers */
  
  /* now must update the concentrations, which can change in both
     layer space and in size space.  I arbitrarily choose the following
     order of processes:
       
       1. Lagrangian growth
       2. Eulerian growth
       3. Activation of new ice nuclei
       4. Mixing due to net fall speed
       5. Aggregation of smaller crystals into larger ones
       6. Advection of water vapor */ 
    
    /* must make sure to zero all variables which are changed by increments */
    (void)memset((char *)CCN_activated,'\0',(num_layer+2)*sizeof(float));
    (void)memset((char *)concentration_new[0],'\0',(num_layer+2)*(num_size+2)*sizeof(float));
    (void)memset((char *)denv_temp_dt,'\0',(num_layer+2)*sizeof(float));
    (void)memset((char *)dvapor_density_dt,'\0',(num_layer+2)*sizeof(float));
    
    /* Find which part of the mass grid should be computed by Lagrangian
       techniques, and which by Eulerian advection.
       Euler_stable_size[layer] is theoretically the first bin size stable
       in an Eulerian advection scheme, one past the last in the Lagrangian. */
    for(layer=1;layer<=num_layer;layer++){
      Euler_stable_size[layer]=1;
      dt_growth=alpha_growth*delta_mass[1]/fabs(dmass_dtime[layer][1]);
      while(dt_growth < dt){
	Euler_stable_size[layer]++;
	dt_growth=alpha_growth*delta_mass[Euler_stable_size[layer]]/
	  fabs(dmass_dtime[layer][Euler_stable_size[layer]]);
      } /* end while */

      /* In case the size/spectral dependence of the crystal heating has
       placed some Lagrangian bins after the first Eulerian bin, correct
       the problem.  Testing shows this only happens when ESS is small. */ 
      if(Euler_stable_size[layer] <= 5){
	int_foo=Euler_stable_size[layer]-1;
	for(size=Euler_stable_size[layer];size<=5;size++){
	  dt_growth=alpha_growth*delta_mass[size]/fabs(dmass_dtime[layer][size]);
	  if(dt_growth < dt){
	    int_foo=size;
	  } /* end if */
	} /* end loop over sizes */
	Euler_stable_size[layer]=int_foo+1;
      } /* end if */

      if(debug == 31){
/*	if(Euler_stable_size[layer] == 1){*/
/*	  (void)fprintf(stdout,"ESS = 1 in layer = %i and timestep = %i\n",*/
/*			layer,time_step); */
/*	} *//* endif */
	for(size=Euler_stable_size[layer];size<=num_size;size++){
	  dt_growth=alpha_growth*delta_mass[size]/fabs(dmass_dtime[layer][size]);
	  if(dt_growth < dt){
	    (void)fprintf(stdout,"Size %i is Euler unstable but ESS = %i.  Layer = %i and time_step = %i\n",size,Euler_stable_size[layer],layer,time_step);
	  } /* end if */
	} /* end loop over sizes */
      } /* end debug */

      if(Euler_stable_size[layer] >= num_size)
	(void)fprintf(stdout,"WARNING:  ESS = %i but num_size = %i in layer = %i and timestep = %i\n",Euler_stable_size[layer],num_size,layer,time_step); 
    } /* end loop over layers */
    
    if(debug == 31){
      max_INTVEC_idxed(Euler_stable_size+cloud_base_idx,num_cloudy_layers,
			 &int_foo,&int_foo2);
      (void)fprintf(stdout,"Time = %i, max(ESS) in cloud = %i @ layer %i (cbi = %i, cti = %i)\n",time_step,int_foo,int_foo2+cloud_base_idx,
		    cloud_base_idx,cloud_top_idx); 
    } /* end debug */

    if(EVEN_STEP){
      /* flip-flop on the ESS to ensure particles communicate between grids */
      for(layer=1;layer<=num_layer;layer++){
	Euler_stable_size[layer]++;
      } /* end loop over layers */
    } /* end if */
    EVEN_STEP = !EVEN_STEP;
    
    if(debug == 23){
      (void)fprintf(stdout,"\nEuler Stable size (index): level (across)\n");
      for(layer=1;layer<=num_layer;layer++){
	(void)fprintf(stdout,"%i\t",Euler_stable_size[layer]);
      } /* end loop over layers */
    } /* end debug */
    if(debug == 24){
      (void)fprintf(stdout,"\nEuler Stable size (microns): level (across)\n");
      for(layer=1;layer<=num_layer;layer++){
	(void)fprintf(stdout,"%.1f\t",
		      crystal_length[Euler_stable_size[layer]]*1.e6);
      } /* end loop over layers */
    } /* end debug */
    
    if(LAGRANGE){
      for(layer=1;layer<=num_layer;layer++){
	for(size=1;size<=Euler_stable_size[layer]-1;size++){
	  if(concentration[layer][size] > 0.){
	    new_size=find_nearest_mass_bin(lagrange_mass[layer][size],
					   crystal_mass,num_size);
	    if(new_size == SIZE_TOO_SMALL){
	      /* perform the vapor sublimation */
	      /* doing nothing here has the effect of sublimating all the
		 particles from this bin, since none are placed in the 
		 conc_new array */
	      if(debug == 30)
		/* this message can't possibly tell the right dv/dt */
		(void)fprintf(stdout,"sublimating %.2f ng vapor back to field at time_step = %i and layer = %i and size = %i\n",dvapor_density_dt[layer]*dt*1.e12,time_step,layer,size);
	    }else if(new_size == SIZE_TOO_BIG){
	      if(VERBOSE)
		(void)fprintf(stdout,"WARNING: largest particle size needs to grow at time_step = %i and layer = %i and size = %i\n",time_step,layer,size);
	    }else{
	      /* perform the mass-grid adjustment */
	      float_foo=conc_left_redistribution
		(concentration[layer][size],lagrange_mass[layer][size],
		 crystal_mass[new_size],crystal_mass[new_size+1]);
	      concentration_new[layer][new_size]+=float_foo;
	      concentration_new[layer][new_size+1]+=
		concentration[layer][size]-float_foo; 
	    } /* end if */
	  } /* end if */
	} /* end loop over sizes */
      } /* end loop over layers */
      
      /* At this point the old concentrations are from the last time step and
	 the new concentrations are from this time step with Lagrangian growth,
	 but no Eulerian growth or velocity advection. Can't swap them until
	 both growth processes are complete or else some particles would get
	 a chance to grow twice. */
    } /* end if LAGRANGE */
    
    if(debug == 28){
      for(layer=1;layer<=num_layer;layer++){
	IWC[layer]=dot_product_VEC(concentration[layer]+1,crystal_mass+1,num_size);
	(void)fprintf(stdout,"post-Lagrangian pre-Eulerian growth IWC[%i] = %g\n",layer,IWC[layer]); 
      } /* end loop over layers */
    } /* end debug */
    
    if(EULER){
      for(layer=1;layer<=num_layer;layer++){
	for(size=Euler_stable_size[layer];size<=num_size;size++){
	  distribution[layer][size]=concentration[layer][size]/delta_mass[size];
	} /* end loop over sizes */
      } /* end loop over layers */
      
      for(layer=1;layer<=num_layer;layer++){
	/* A word about these boundary conditions:  
	   always having fixed LBC means no possible concentration
	   fluxes from ETBFCT will be encountered. */
	
	dmass_dtime_LBC=0.;
	dmass_dtime_RBC=0.;
	mass_LBC=.5*(crystal_mass[Euler_stable_size[layer]-1]+
		     crystal_mass[Euler_stable_size[layer]]);
	
	num_Euler_stable_size=num_size-Euler_stable_size[layer]+1;

	FORTRAN_ngride(crystal_mass+Euler_stable_size[layer],&num_Euler_stable_size,
		&mass_LBC,&mass_RBC,&etbfct_alpha);
	FORTRAN_ogride(&num_Euler_stable_size);
	FORTRAN_ngride(crystal_mass+Euler_stable_size[layer],&num_Euler_stable_size,
		&mass_LBC,&mass_RBC,&etbfct_alpha);
	FORTRAN_veloce(dmass_dtime[layer]+Euler_stable_size[layer],
		&num_Euler_stable_size,&dmass_dtime_LBC,&dmass_dtime_RBC,&dt);
	/*	FORTRAN_consre(distribution[layer]+Euler_stable_size[layer],
		&num_Euler_stable_size,&float_foo_1); */
	FORTRAN_etbfct(distribution[layer]+Euler_stable_size[layer],
		distribution_new[layer]+Euler_stable_size[layer],
		&num_Euler_stable_size,&dist_LBC,&dist_RBC);

      } /* end loop over layers */
      
      if(WEED_NEGATIVE_VALUES){
	for(layer=1;layer<=num_layer;layer++){
	  for(size=1;size<=num_size;size++){
	    if(distribution_new[layer][size] < 0.)
	      distribution_new[layer][size]=0.;
	  } /* end loop over sizes */
	} /* end loop over layers */
      }/* end WEED */
      
      /* Add the Eulerian advected masses to the Lagrangian ones already 
	 stored in the new concentration array.
	 Find the change in these thermodynamic properties due to the
	 Eulerian mass growth */
      for(layer=1;layer<=num_layer;layer++){
	for(size=Euler_stable_size[layer];size<=num_size;size++){
	  concentration_new[layer][size]+=distribution_new[layer][size]*
	    delta_mass[size];
	} /* end loop over sizes */
      } /* end loop over layers */
    } /* end if EULER */
    
    /* Seed the cloud with active deposition/condensation ice-nuclei */
    if(HETEROGENEOUS_NUCLEATION){
      for(layer=1;layer<=num_layer;layer++){
	if(env_temperature[layer] >= 253.15 && saturation_ice[layer] >= 1.){
	  concentration_new[layer][1]+=
	    heterogeneous_nucleation_MDC92
	      (saturation_ice[layer],dt,tau_CCN_timescale);
	} /* end if */ 
      } /* end loop over layers */
    } /* end if HETEROGENEOUS_NUCLEATION */
    
    /* Activate some deliquesent sulphate aerosol droplets */
    if(HOMOGENEOUS_NUCLEATION){
      for(layer=1;layer<=num_layer;layer++){
	if(saturation_liquid[layer] >= solute_deliquescence_point){
	  CCN_activated[layer]=homogeneous_nucleation_DMC94
	    (CCN_conc,
	     CCN_diameter,
	     env_temperature[layer],
	     num_CCN_size,
	     saturation_liquid[layer]);
	  concentration_new[layer][1]+=CCN_activated[layer];
	} /* end if */ 
      } /* end loop over layers */
    } /* end if HOMOGENEOUS_NUCLEATION */
    
    if(EULER || LAGRANGE || HETEROGENEOUS_NUCLEATION || HOMOGENEOUS_NUCLEATION 
       || RADIATION){
      /* Compute the change in the vapor density and temperature of the ambient air
	 in each layer, following notation of Zhang '89 p.310 and method of
	 Heymsfield '75 p. 821 and czp 64-65. Make sure this is done before
	 the vapor is advected. */
      
      /* the change in temperature depends only on the local change in the
	 vapor concentration in the field so update these two variables
	 before beginning the spatial advection phases */
      
      env_temperature[0]=env_temperature[1]-
	(env_temperature[2]-env_temperature[1]);
      env_temperature[num_layer+1]=env_temperature[num_layer]+
	(env_temperature[num_layer]-env_temperature[num_layer-1]);
      
      for(layer=1;layer<=num_layer;layer++){
	for(size=1;size<=num_size;size++){
	  dvapor_density_dt[layer]+=-crystal_mass[size]*
	    (concentration_new[layer][size]-concentration[layer][size])/dt;
	} /* end loop over sizes */
	
	vapor_density[layer]+=dt*dvapor_density_dt[layer];
	latent_heating[layer]=-latent_heat_sub*dvapor_density_dt[layer]/
	  (spec_heat_air*env_density[layer]);
	/* recall that Gamma is defined as minus the temperature gradient */ 
	advective_heating[layer]=
	  -wind_speed[layer]*(adiabatic_lapse+
	     (.5*(env_temperature[layer+1]-env_temperature[layer-1])/dz));
	denv_temp_dt[layer]=latent_heating[layer]+advective_heating[layer];
      } /* end loop over layers */
      
      if(RAD_FEEDBACK){
	for(layer=1;layer<=num_layer;layer++){
	  /* add in the the non-convection compensated radiative heating rate */ 
	  denv_temp_dt[layer]+=heating_rate_net[layer];
	} /* end loop over layers */
      } /* endif RAD_FEEDBACK */ 
      
      if(debug == 20){
	for(layer=1;layer<=num_layer;layer++){
	  (void)fprintf(stdout,"denv_temp_dt[%i] = %g\n",
			layer,denv_temp_dt[layer]); 
	} /* end loop over layers */
      } /* end debug */
      
      for(layer=1;layer<=num_layer;layer++){
	/* NB: were the fixed density profile/exponential wind scheme to 
	 be used, here is where denv_temp_dt would be divided by (1.-kappa) */ 
	env_temperature[layer]+=dt*denv_temp_dt[layer];
      } /* end loop over layers */
      
      /* At this point the old concentrations are from this time step and
	 have undergone both Lagrangian growth, Eulerian growth, and Nucleation
	 but no velocity advection. Now we'll swap them and perform the 
	 velocity advection phase. */
      concentration_swap=concentration_new;
      concentration_new=concentration;
      concentration=concentration_swap; 
    } /* end if EULER || LAGRANGE || HETEROGENEOUS_NUCLEATION || HOMOGENEOUS_NUCLEATION || RADIATION */

    if(debug == 28){
      for(layer=1;layer<=num_layer;layer++){
	IWC[layer]=dot_product_VEC(concentration[layer]+1,crystal_mass+1,num_size);
	(void)fprintf(stdout,"post-Lagrangian post-Eulerian growth IWC[%i] = %g\n",layer,IWC[layer]); 
      } /* end loop over layers */
    }
    
    if(AGGREGATION && (time_step%agg_step == 0)){
      if(VERBOSE)
	(void)fprintf(stdout,"\n\ntime = %i = %g seconds\n",
		      time_step,time_array[time_step-1]+dt);
      for(layer=cloud_base_idx;layer<=cloud_top_idx;layer++){
	if(VERBOSE)
	  (void)fprintf(stdout,"\n\nlayer = %i = %.2f km\n",
			layer,altitude[layer]*KILOMETERS_PER_METER);
	for(size=smallest_rimer_bin;size<=largest_rimer_bin;size++){
	  if(concentration[layer][size] > 1.){
	    if(VERBOSE)
	      (void)fprintf(stdout,"\nrimer size = %i, conc = %g\n",
			    size,concentration[layer][size]);
	    for(int_foo=smallest_rimee_bin;int_foo<=largest_rimee_bin;int_foo++){
	      if(concentration[layer][int_foo] > 1.){
		collision_kernel=pi_over_four*
		  (crystal_length[size]+crystal_length[int_foo])*
		    (crystal_length[size]+crystal_length[int_foo])*
		      fabs(fall_speed[size][layer]-fall_speed[int_foo][layer]);
		num_collisions=
		  concentration[layer][int_foo]*collision_kernel*dt*agg_step;
		/* prevent overcounting by multiplying by one half */ 
		num_collections=.5*num_collisions*sticking_efficiency;
		aggregate_mass=
		  crystal_mass[size]+num_collections*crystal_mass[int_foo];
		new_size=
		  find_nearest_mass_bin(aggregate_mass,crystal_mass,num_size);
		if(new_size != size){
		  (void)fprintf(stdout,"new_size != size in AGGREGATE\n");
		} /* endif */ 
		if(new_size == SIZE_TOO_SMALL){
		  (void)fprintf(stdout,"new_size too small in AGGREGATE\n");
		}else if(new_size == SIZE_TOO_BIG){
		  (void)fprintf(stdout,"new_size too big in AGGREGATE\n");
		}else{
		  /* NB:  algorithm works correctly only when new_size == size */ 
		  float_foo=conc_left_redistribution
		    (concentration[layer][size],aggregate_mass,
		     crystal_mass[new_size],crystal_mass[new_size+1]);

		  /* only carry through if it matters */ 
		  if(concentration[layer][size]-float_foo < EMINUSTHREE){
		    continue; /* go on the next rimee size */ 
		  } /* end if */

		  if(VERBOSE){
		    (void)fprintf(stdout,"\nrimee size = %i, conc = %g\n",
				  int_foo,concentration[layer][int_foo]);
		    (void)fprintf(stdout,"coll kernel: %g, num_collisions: %g, num_collections: %g\n",collision_kernel,num_collisions,num_collections);
		    (void)fprintf(stdout,"rimer mass: %g, rimee mass: %g, aggretate mass: %g\n",crystal_mass[size],crystal_mass[int_foo],aggregate_mass);
		    (void)fprintf(stdout,"removing %e rimees m^-3 from rimee bin %i\n",concentration[layer][size]*num_collections,int_foo);
		    (void)fprintf(stdout,"adding %e aggregates m^-3 to rimed bin %i\n",concentration[layer][size]-float_foo,new_size+1);
		    (void)fprintf(stdout,"replacing %e rimers by %e crystals m^-3 in rimer bin %i\n",concentration[layer][size],float_foo,new_size);
		  } /* end if */

		  /* keep track of the total mass flow */ 
		  float_foo2=-concentration[layer][size]*num_collections*crystal_mass[int_foo];
		  /* adjust the rimee the mass-bin */
		  concentration[layer][int_foo]-=
		    concentration[layer][size]*num_collections;

		  /* keep track of the total mass flow */ 
		  float_foo2+=(concentration[layer][size]-float_foo)*crystal_mass[new_size+1]; 
		  /* apportion the aggregates into the larger bin */ 
		  concentration[layer][new_size+1]+=
		    concentration[layer][size]-float_foo; 
		  
		  /* keep track of the total mass flow */ 
		  float_foo2-=(concentration[layer][size]-float_foo)*crystal_mass[new_size]; 
		  /* remaining aggregates stay in the the rimer bin */ 
		  concentration[layer][new_size]=float_foo;
		  
		  if(VERBOSE)
		    (void)fprintf(stdout,"net mass change = %e\n",float_foo2);
		} /* end if */
	      } /* end if conc rimees > 1 */
	    } /* end loop over int_foos */
	  } /* end if conc rimers > 1 */
	} /* end loop over sizes */
      } /* end loop over layers */
    } /* end if */
    
    if(ICE_ADVECT){
      /* must fill the transposed concentration array once per time step */
      transpose_conc(concentration,conc_transpose,num_layer,num_size);
      
      for(size=1;size<=num_size;size++){
	daltitude_dtime_LBC=daltitude_dtime[size][1];
	daltitude_dtime_RBC=daltitude_dtime[size][num_layer]; 
	if(debug_value == 3){
	  daltitude_dtime_LBC=0.;
	  daltitude_dtime_RBC=0.;
	}

	FORTRAN_ngride(altitude+1,&num_layer,&altitude_LBC,&altitude_RBC,&etbfct_alpha);
	FORTRAN_ogride(&num_layer);
	FORTRAN_ngride(altitude+1,&num_layer,&altitude_LBC,&altitude_RBC,&etbfct_alpha);
	FORTRAN_veloce(daltitude_dtime[size]+1,&num_layer,&daltitude_dtime_LBC,
		&daltitude_dtime_RBC,&dt);
	/* FORTRAN_consre(conc_transpose[size]+1,&num_layer,&float_foo); */
	FORTRAN_etbfct(conc_transpose[size]+1,conc_transpose_new[size]+1,
		&num_layer,&conc_LBC,&conc_RBC);
      } /* end loop over sizes */

      /* must now place the transposed new concentration into conc_new */
      transpose_trans(conc_transpose_new,concentration_new,num_size,num_layer);
      
      /* At this point the old concentrations are from the this time step and
	 have undergone four mixing/growth processes. Only the water vapor
	 advection remains to be done.  Now we'll swap them to prepare for
	 the next time step. */
      concentration_swap=concentration_new;
      concentration_new=concentration;
      concentration=concentration_swap; 
      
      if(WEED_NEGATIVE_VALUES){
	for(layer=1;layer<=num_layer;layer++){
	  for(size=1;size<=num_size;size++){
	    if(concentration[layer][size] < 0.)
	      concentration[layer][size]=0.;
	  } /* end loop over sizes */
	} /* end loop over layers */
      }/* end WEED */
      
    } /* end if ICE_ADVECT */
    
    if(VAPOR_ADVECT){
      /* Setting the wind speed b.c.'s = 0 would conserve the initial total
	 water path, except possibly for ice fallout. Setting them equal to
	 their nearest neighbors ensures 'free' transport of vapor from below */
      wind_speed_LBC=wind_speed[1];
      wind_speed_RBC=wind_speed[num_layer]; 
      if(debug_value == 3){
	wind_speed_LBC=0.;
	wind_speed_RBC=0.;
      }

      FORTRAN_ngride(altitude+1,&num_layer,&altitude_LBC,&altitude_RBC,&etbfct_alpha);
      FORTRAN_ogride(&num_layer);
      FORTRAN_ngride(altitude+1,&num_layer,&altitude_LBC,&altitude_RBC,&etbfct_alpha);
      FORTRAN_veloce(wind_speed+1,&num_layer,&wind_speed_LBC,
	      &wind_speed_RBC,&dt);
      if(debug == 37){
	FORTRAN_consre(vapor_density+1,&num_layer,&float_foo); 
	(void)fprintf(stdout,"Vapor Path before time_step = %i is %g g/m^2\n",
		      time_step,float_foo*1000.); 
      }
      FORTRAN_etbfct(vapor_density+1,vapor_density_new+1,
	      &num_layer,&vapor_LBC,&vapor_RBC);
      if(debug == 37){
	FORTRAN_consre(vapor_density_new+1,&num_layer,&float_foo); 
	(void)fprintf(stdout,"Vapor Path after time_step = %i is %g g/m^2\n",
		      time_step,float_foo*1000.); 
      }

      /* At this point the old vapor densities are from this time step and
	 have undergone all the 3 growth processes. The new vapor densities
	 have also been advected by the large scale wind field.
	 Now we'll swap them to prepare for the next time step. */
      vapor_density_swap=vapor_density_new;
      vapor_density_new=vapor_density;
      vapor_density=vapor_density_swap; 
      
      if(WEED_NEGATIVE_VALUES){
	for(layer=1;layer<=num_layer;layer++){
	  if(vapor_density[layer] <= 0.){
	    (void)fprintf
	      (stdout,"WARNING: vapor_density[%i] = %f at time_step = %i\n",
	       layer,vapor_density[layer],time_step);
	    /* vapor densities too close to zero will crash radcsw() in 
	       the H2O path length computations */ 
	    vapor_density[layer]=EMINUSFIFTEEN;
	  } /* endif */ 
	} /* end loop over layers */
      }/* end WEED */
    } /* end if VAPOR_ADVECT */
    
    if(debug_value == 5){
      if(VAPOR_ADVECT){
	wind_speed[0]=wind_speed[1];
	for(layer=1;layer<=num_layer;layer++){
	  wind_speed_new[layer]=wind_speed[layer]+dt*
	    (-wind_speed[layer]*(wind_speed[layer]-wind_speed[layer-1])/dz +
	     gravity*(1./orig_env_temp[layer])*
	     (env_temperature[layer]-orig_env_temp[layer]) -
	     gravity*IWC[layer]/env_density[layer]);
	} /* end loop over layers */
	
	/* interchange the wind_speed vectors */ 
	wind_speed_swap=wind_speed_new;
	wind_speed_new=wind_speed;
	wind_speed=wind_speed_swap; 
	
	/* adjust the combined fall speeds */ 
	for(layer=1;layer<=num_layer;layer++){
	  for(size=1;size<=num_size;size++){
	    daltitude_dtime[size][layer]=fall_speed[size][layer]+wind_speed[layer];
	  } /* end loop over sizes */
	} /* end loop over layers */
	
      }/* endif VAPOR_ADVECT */
    } /* end debug_value */
    
    /* Now all the prognostic variables have been adjusted for the current time 
       step, so the hard part is over. It remains to recalculate those diagnostic
       variables which are useful in calculations and visualization */
    for(layer=1;layer<=num_layer;layer++){
      eqm_vap_ice[layer]=Eqm_Vap_Ice(env_temperature[layer]);
      eqm_vap_liquid[layer]=Eqm_Vap_Liquid(env_temperature[layer]);
      mmr_vapor[layer]=vapor_density[layer]/env_density[layer];
      pp_vapor[layer]=env_pressure[layer]*mmr_vapor[layer]/
	(epsilon_vapor+mmr_vapor[layer]);
      saturation_ice[layer]=pp_vapor[layer]/eqm_vap_ice[layer];
      saturation_liquid[layer]=pp_vapor[layer]/eqm_vap_liquid[layer];
      
      IWC[layer]=dot_product_VEC(concentration[layer]+1,crystal_mass+1,num_size); 
      IWP[layer]=dz*IWC[layer];
      total_layer_conc[layer]=total_VEC(concentration[layer]+1,num_size);
      vapor_path[layer]=dz*env_density[layer]*mmr_vapor[layer];
      
      effective_radius[layer]=
	(max_VEC(concentration[layer]+1,num_size) > EMINUSNINE) ?
	  Effective_Radius(concentration[layer],equiv_rad_squared,
			   equiv_rad_cubed,num_size) : ONE_MICRON;
    } /* end loop over layers */
    
    /* Update the time series arrays */
    cloud_base_idx=Cloud_Base_idx(IWC,num_layer,time_step,total_layer_conc);
    cloud_top_idx=Cloud_Top_idx(IWC,num_layer,time_step,total_layer_conc);
    num_cloudy_layers=cloud_top_idx-cloud_base_idx+1;

    cloud_base_of_time[time_step]=altitude[cloud_base_idx];
    cloud_top_of_time[time_step]=altitude[cloud_top_idx];

    cloud_effective_radius[time_step]=0.;
    for(layer=cloud_base_idx;layer<=cloud_top_idx;layer++){
      cloud_effective_radius[time_step]+=effective_radius[layer]*
	total_layer_conc[layer];
    } /* end loop over layers */
    cloud_effective_radius[time_step]/=((float_foo=total_VEC(total_layer_conc+cloud_base_idx,num_cloudy_layers)) > EMINUSTWELVE) ? float_foo : 1.;

    saturation_ice_of_time[time_step]=max_VEC(saturation_ice+1,num_layer);
    if(saturation_ice_of_time[time_step]/saturation_ice_of_time[time_step-1] > 1.03){
      (void)fprintf(stdout,"Hiccup Belch and Fart at time step = %i\n",time_step);
    } /* end if */
    
    surface_area_of_time[time_step]=0.;
    for(layer=1;layer<=num_layer;layer++){
      surface_area_of_time[time_step]+=4.*M_PI*
	dot_product_VEC(equiv_rad_squared+1,concentration[layer]+1,num_size);
    } /* end loop over layers */
    
    if(RADIATION){
      if((time_step-1)%rad_step == 0){
	cloud_optical_depth[time_step-1]=0.;
	for(layer=cloud_base_idx;layer<=cloud_top_idx;layer++){
	  cloud_optical_depth[time_step-1]+=eff_tau_extinction[layer][8];
	} /* end loop over layers */
      }else{
	albedo_of_time[time_step-1]=albedo_of_time[time_step-2];
	emissivity_of_time[time_step-1]=emissivity_of_time[time_step-2];
	SW_cloud_forcing[time_step-1]=SW_cloud_forcing[time_step-2];
	LW_cloud_forcing[time_step-1]=LW_cloud_forcing[time_step-2];
	cloud_optical_depth[time_step-1]=cloud_optical_depth[time_step-2];
      } /* end if */
    } /* end if */

    max_length_of_time[time_step]=0.;
    for(size=num_size;size>=1;size--){
      for(layer=1;layer<=num_layer;layer++){
	if(concentration[layer][size] >= EMINUSTHREE)
	  max_length_of_time[time_step]=crystal_length[size];
      } /* end loop over layers */
      if(max_length_of_time[time_step] > 0.) break;
    } /* end loop over sizes */
    
    IWP_of_time[time_step]=total_VEC(IWP+1,num_layer);
    vapor_path_of_time[time_step]=total_VEC(vapor_path+1,num_layer);
    water_path_of_time[time_step]=IWP_of_time[time_step]+
      vapor_path_of_time[time_step];
    time_array[time_step]=time_array[time_step-1]+dt;
    (void)fprintf(stdout,".");
    (void)fflush(stdout);
    
    /* time to update the snapshots? */ 
    if(time_step%plot_step == 0){
      frame++;
      time_snapshot[frame]=time_array[time_step];
      for(layer=1;layer<=num_layer;layer++){
	IWC_snapshot[frame][layer]=IWC[layer];
	number_snapshot[frame][layer]=total_layer_conc[layer];
      } /* end loop over layers */

      if(MOVIE_NETCDF_OUTPUT){
	/* append a scalar to a one-dimensional array */
	start[0]=(long)frame;
	ncvarput1(cdfid,movie_cloud_effective_radius_id,start,(void *)(cloud_effective_radius+frame));
	
	/* append a one-dimensional array to a two-dimensional array */
	start[0]=(long)frame;
	start[1]=0L;
	
	count[0]=1L;
	count[1]=num_layerp2;
	ncvarput(cdfid,movie_effective_radius_id,start,count,(void *)effective_radius);
	
	/* append a two dimensional array to a three-dimensional array */
	start[0]=(long)frame;
	start[1]=0L;
	start[2]=0L;
	
	count[0]=1L;
	count[1]=num_layerp2;
	count[2]=num_sizep2;
	ncvarput(cdfid,movie_concentration_id,start,count,(void *)concentration_ptr);
      } /* endif MOVIE_NETCDF_OUTPUT */ 
    } /* end if snapshot update */
    
    if(debug == 12){
      (void)fprintf(stdout,"time step = %i, concentration[%i][%i] = %g\n",
		    time_step,1,1,concentration[1][1]);
    } /* end debug */
    
    if(debug == 79){
      (void)fprintf(stdout,"%i TOA SWCF = %g, Surf SWCF = %g, Surf/TOA = %g\n",
		    time_step-1,
		    SW_cloud_forcing[time_step-1],
		    SW_surf_cloud_forcing,
		    SW_surf_cloud_forcing/SW_cloud_forcing[time_step-1]);
    } /* end debug */

  } /* end outer time loop */
  
  /* decrement time_step by one so it represents actual # of steps computed */
  time_step--;
  
  (void)fprintf(stdout,"\n");
  
  /* Finish up some calculations useful for plotting */ 
  for(layer=1;layer<=num_layer;layer++){
    for(size=1;size<=num_size;size++){
      distribution[layer][size]=concentration[layer][size]/delta_mass[size];
    } /* end loop over sizes */
  } /* end loop over layers */
  
  if(RADIATION){
    for(layer=num_layer;layer>=1;layer--){
      optical_depth[layer]=eff_tau_extinction[layer][8];
    } /* end loop over layers */
  } /* end if */

  (void)fprintf(stdout,"\nstep\ttime\ttime\tBase\tTop\tMax L\tReff\tSat_i\tIWP\tH2O\n");
  (void)fprintf(stdout," #\t sec\t min\t km \t km\tmicrons\tmicrons\t\tg/m^2\tg/m^2\n");
  for(step=0;step<=num_step;step++){
    (void)fprintf(stdout,"%i\t%.1f\t%.2f\t%.2f\t%.2f\t%.1f\t%.1f\t%.4f\t%.4f\t%.4f\n",
		  step,
		  time_array[step],
		  time_array[step]/60.,
		  cloud_base_of_time[step]*KILOMETERS_PER_METER,
		  cloud_top_of_time[step]*KILOMETERS_PER_METER,
		  max_length_of_time[step]*MICRONS_PER_METER,
		  cloud_effective_radius[step]*MICRONS_PER_METER,
		  saturation_ice_of_time[step],
		  IWP_of_time[step]*GRAMS_PER_KILOGRAM,
		  water_path_of_time[step]*GRAMS_PER_KILOGRAM);
  } /* end loop over time steps */
  
  if(RADIATION){
    (void)fprintf(stdout,"\nstep\ttime\ttime\tLWCF\tSWCF\tNetCF\tTau\tAlbedo\tEmissivity\n");
    (void)fprintf(stdout," #\t sec\t min\tW/m^2\tW/m^2\tW/m^2\t\t\t\n");
    for(step=0;step<=num_step-1;step++){
      (void)fprintf(stdout,"%i\t%.1f\t%.2f\t%.2f\t%.2f\t%.2f\t%.3f\t%.3f\t%.3f\n",
		    step,
		    time_array[step],
		    time_array[step]/60.,
		    LW_cloud_forcing[step],
		    SW_cloud_forcing[step],
		    LW_cloud_forcing[step]+SW_cloud_forcing[step],
		    cloud_optical_depth[step],
		    albedo_of_time[step],
		    emissivity_of_time[step]);
    } /* end loop over time steps */
  } /* end if */
  
  if(debug == 10){
    (void)fprintf(stdout,"\nzbin\talt\tTemp K\tTemp C\tS_ice\tS_liq\tP\tr(Vap)\tr(Ice)\n");
    (void)fprintf(stdout,"\tkm\tK\tC\t\t\tmb\tg/kg\tg/kg\n");
    for(layer=1;layer<=num_layer;layer++){
      (void)fprintf(stdout,"%i\t%.2f\t%.1f\t%.1f\t%.2f\t%.2f\t%.0f\t%.3f\t%.3f\n",
		    layer,
		    altitude[layer]/1000.,
		    env_temperature[layer],
		    env_temperature[layer]-water_freezing_point,
		    saturation_ice[layer],
		    saturation_liquid[layer],
		    env_pressure[layer]/100.,
		    mmr_vapor[layer]*1000.,
		    IWC[layer]*1000./env_density[layer]);
    } /* end loop over layers */
  } /* endif debug */
  
  if(debug == 25){
    for(size=1;size<=num_size;size++){
      (void)fprintf
	(stdout,"at end of run, concentration[%i][%i] = %g per liter\n",
	 (int)float_foo_input,size,
	 concentration[(int)float_foo_input][size]/LITERS_PER_CUBIC_METER);
    } /* end loop over sizes */
  } /* end debug */ 
  
  if(NETCDF_OUTPUT){
    /* store the scalars */
    start[0]=0L;

    ncvarput1(cdfid,dt_id,start,(void *)&dt);
    ncvarput1(cdfid,dz_id,start,(void *)&dz);
    ncvarput1(cdfid,num_CCN_size_id,start,(void *)&num_CCN_size);
    ncvarput1(cdfid,num_band_id,start,(void *)&num_band);
    ncvarput1(cdfid,num_ccm2_level_id,start,(void *)&num_ccm2_level);
    ncvarput1(cdfid,num_cloudy_layers_id,start,(void *)&num_cloudy_layers);
    ncvarput1(cdfid,num_frame_id,start,(void *)&num_frame);
    ncvarput1(cdfid,num_layer_id,start,(void *)&num_layer);
    ncvarput1(cdfid,num_size_id,start,(void *)&num_size);
    ncvarput1(cdfid,num_step_id,start,(void *)&num_step);
    ncvarput1(cdfid,agg_step_id,start,(void *)&agg_step);
    ncvarput1(cdfid,rad_step_id,start,(void *)&rad_step);
    ncvarput1(cdfid,plot_step_id,start,(void *)&plot_step);
    ncvarput1(cdfid,PAUSE_id,start,(void *)&PAUSE);

    /* store the one-dimensional arrays */
    start[0]=0L;

    ncvarput(cdfid,CCN_activated_id,start,&num_layerp2,(void *)CCN_activated);
    ncvarput(cdfid,CCN_conc_id,start,&num_CCN_sizep2,(void *)CCN_conc);
    ncvarput(cdfid,CCN_diameter_id,start,&num_CCN_sizep2,(void *)CCN_diameter);
    ncvarput(cdfid,CCN_mass_id,start,&num_CCN_sizep2,(void *)CCN_mass);
    ncvarput(cdfid,aspect_ratio_id,start,&num_sizep2,(void *)aspect_ratio);
    ncvarput(cdfid,crystal_diameter_id,start,&num_sizep2,(void *)crystal_diameter);
    ncvarput(cdfid,saturation_liquid_id,start,&num_layerp2,(void *)saturation_liquid);
    ncvarput(cdfid,IWC_id,start,&num_layerp2,(void *)IWC);
    ncvarput(cdfid,IWP_id,start,&num_layerp2,(void *)IWP);
    ncvarput(cdfid,IWP_of_time_id,start,&num_stepp1,(void *)IWP_of_time);
    ncvarput(cdfid,LW_cloud_forcing_id,start,&num_stepp1,(void *)LW_cloud_forcing);
    ncvarput(cdfid,SW_cloud_forcing_id,start,&num_stepp1,(void *)SW_cloud_forcing);
    ncvarput(cdfid,advective_heating_id,start,&num_layerp2,(void *)advective_heating);
    ncvarput(cdfid,albedo_of_time_id,start,&num_stepp1,(void *)albedo_of_time);
    ncvarput(cdfid,altitude_id,start,&num_layerp2,(void *)altitude);
    ncvarput(cdfid,cloud_base_of_time_id,start,&num_stepp1,(void *)cloud_base_of_time);
    ncvarput(cdfid,cloud_effective_radius_id,start,&num_stepp1,(void *)cloud_effective_radius);
    ncvarput(cdfid,cloud_top_of_time_id,start,&num_stepp1,(void *)cloud_top_of_time);
    ncvarput(cdfid,crystal_length_id,start,&num_sizep2,(void *)crystal_length);
    ncvarput(cdfid,crystal_mass_id,start,&num_sizep2,(void *)crystal_mass);
    ncvarput(cdfid,delta_length_id,start,&num_sizep2,(void *)delta_length);
    ncvarput(cdfid,delta_mass_id,start,&num_sizep2,(void *)delta_mass);
    ncvarput(cdfid,effective_radius_id,start,&num_layerp2,(void *)effective_radius);
    ncvarput(cdfid,emissivity_of_time_id,start,&num_stepp1,(void *)emissivity_of_time);
    ncvarput(cdfid,env_density_id,start,&num_layerp2,(void *)env_density);
    ncvarput(cdfid,env_pressure_id,start,&num_layerp2,(void *)env_pressure);
    ncvarput(cdfid,env_temperature_id,start,&num_layerp2,(void *)env_temperature);
    ncvarput(cdfid,flux_down_LW_id,start,&num_layerp2,(void *)flux_down_LW);
    ncvarput(cdfid,flux_down_SW_id,start,&num_layerp2,(void *)flux_down_SW);
    ncvarput(cdfid,flux_up_LW_id,start,&num_layerp2,(void *)flux_up_LW);
    ncvarput(cdfid,flux_up_SW_id,start,&num_layerp2,(void *)flux_up_SW);
    ncvarput(cdfid,heating_rate_LW_id,start,&num_layerp2,(void *)heating_rate_LW);
    ncvarput(cdfid,heating_rate_SW_id,start,&num_layerp2,(void *)heating_rate_SW);
    ncvarput(cdfid,heating_rate_net_id,start,&num_layerp2,(void *)heating_rate_net);
    ncvarput(cdfid,latent_heating_id,start,&num_layerp2,(void *)latent_heating);
    ncvarput(cdfid,max_length_of_time_id,start,&num_stepp1,(void *)max_length_of_time);
    ncvarput(cdfid,mmr_vapor_id,start,&num_layerp2,(void *)mmr_vapor);
    ncvarput(cdfid,cloud_optical_depth_id,start,&num_stepp1,(void *)cloud_optical_depth);
    ncvarput(cdfid,orig_env_temp_id,start,&num_layerp2,(void *)orig_env_temp);
    ncvarput(cdfid,orig_mmr_vapor_id,start,&num_layerp2,(void *)orig_mmr_vapor);
    ncvarput(cdfid,potential_temp_id,start,&num_layerp2,(void *)potential_temp);
    ncvarput(cdfid,pp_vapor_id,start,&num_layerp2,(void *)pp_vapor);
    ncvarput(cdfid,saturation_ice_of_time_id,start,&num_stepp1,(void *)saturation_ice_of_time);
    ncvarput(cdfid,saturation_ice_id,start,&num_layerp2,(void *)saturation_ice);
    ncvarput(cdfid,surface_area_of_time_id,start,&num_stepp1,(void *)surface_area_of_time);
    ncvarput(cdfid,time_array_id,start,&num_stepp1,(void *)time_array);
    ncvarput(cdfid,time_snapshot_id,start,&num_framep1,(void *)time_snapshot);
    ncvarput(cdfid,vapor_density_id,start,&num_layerp2,(void *)vapor_density);
    ncvarput(cdfid,vapor_path_id,start,&num_layerp2,(void *)vapor_path);
    ncvarput(cdfid,vapor_path_of_time_id,start,&num_stepp1,(void *)vapor_path_of_time);
    ncvarput(cdfid,optical_depth_id,start,&num_layerp2,(void *)optical_depth);
    ncvarput(cdfid,water_path_of_time_id,start,&num_stepp1,(void *)water_path_of_time);
    ncvarput(cdfid,wind_speed_id,start,&num_layerp2,(void *)wind_speed);

    /* store two-dimensional arrays */
    start[0]=0L;
    start[1]=0L;

    count[0]=num_framep1; 
    count[1]=num_layerp2;
    ncvarput(cdfid,IWC_snapshot_id,start,count,(void *)IWC_snapshot_ptr);
    count[0]=num_framep1; 
    count[1]=num_layerp2;
    ncvarput(cdfid,number_snapshot_id,start,count,(void *)number_snapshot_ptr);
    count[0]=num_layerp2;
    count[1]=num_sizep2;
    ncvarput(cdfid,concentration_id,start,count,(void *)concentration_ptr);
    count[0]=num_layerp2;
    count[1]=num_sizep2;
    ncvarput(cdfid,distribution_id,start,count,(void *)distribution_ptr);
    count[0]=num_layerp2;
    count[1]=num_sizep2;
    ncvarput(cdfid,crystal_heat_SW_id,start,count,(void *)crystal_heat_SW_ptr);
    count[0]=num_layerp2;
    count[1]=num_sizep2;
    ncvarput(cdfid,crystal_heat_LW_id,start,count,(void *)crystal_heat_LW_ptr);
    count[0]=num_layerp2;
    count[1]=num_bandp2;
    ncvarput(cdfid,eff_asymmetry_param_id,start,count,(void *)eff_asymmetry_param_ptr);
    count[0]=num_layerp2;
    count[1]=num_bandp2;
    ncvarput(cdfid,eff_ss_albedo_id,start,count,(void *)eff_ss_albedo_ptr);
    count[0]=num_layerp2;
    count[1]=num_bandp2;
    ncvarput(cdfid,eff_tau_extinction_id,start,count,(void *)eff_tau_extinction_ptr);
    count[0]=num_sizep2;
    count[1]=num_layerp2;
    ncvarput(cdfid,daltitude_dtime_id,start,count,(void *)daltitude_dtime_ptr);

    ncclose(cdfid);
    (void)fprintf(stdout,"Wrote netCDF cloud data to %s\n",out_file);
  } /* end if netCDF output */ 
  
  if((debug == 19) || (debug == 29)){
    (void)fprintf(stdout,"\nFinal Concentration of crystals(#/liter): size bin (across) vs layer (down)\n");
    for(size=1;size<=num_size;size++){
      (void)fprintf(stdout,"\t%i",size);
    } /* end loop over sizes */
    for(layer=1;layer<=num_layer;layer++){
      (void)fprintf(stdout,"\n%i\t",layer);
      for(size=1;size<=num_size;size++){
	(void)fprintf
	  (stdout,"%.3f\t",concentration[layer][size]/LITERS_PER_CUBIC_METER);
      } /* end loop over sizes */
    } /* end loop over layers */

    (void)fprintf(stdout,"\nFinal Distribution of crystals(#/liter/micron): size bin (across) vs layer (down)\n");
    for(size=1;size<=num_size;size++){
      (void)fprintf(stdout,"\t%i",size);
    } /* end loop over sizes */
    for(layer=1;layer<=num_layer;layer++){
      (void)fprintf(stdout,"\n%i\t",layer);
      for(size=1;size<=num_size;size++){
	(void)fprintf
	  (stdout,"%.3f\t",concentration[layer][size]/
	   (delta_length[size]*LITERS_PER_CUBIC_METER*MICRONS_PER_METER));
      } /* end loop over sizes */
    } /* end loop over layers */
  }/* end debug */
  
  if(debug == 13){
    for(layer=1;layer<=num_layer;layer++){
	(void)fprintf(stdout,"thermal conductivity[%i] = %g\n",
		      layer,thermal_conductivity[layer]); 
    } /* end loop over layers */
  } /* end debug */
  
  if(debug == 21){
    for(layer=1;layer<=num_layer;layer++){
      (void)fprintf(stdout,"vapor_diffusivity[%i] = %g\n",layer,
		    vapor_diffusivity[layer]); 
    } /* end loop over layers */
  } /* end debug */
  
  if(debug == 15){
    (void)fprintf(stdout,"\nDaltitude/Dtime (cm/s) of crystals: size bin (across) vs layer (down)\n");
    for(size=1;size<=num_size;size++){
      (void)fprintf(stdout,"\t%i",size);
    } /* end loop over sizes */
    for(layer=1;layer<=num_layer;layer++){
      (void)fprintf(stdout,"\n%i\t",layer);
      for(size=1;size<=num_size;size++){
	(void)fprintf(stdout,"%.1f\t",daltitude_dtime[size][layer]*100.);
      } /* end loop over sizes */
    } /* end loop over layers */
  } /* end debug */
  
  if(debug == 16){
    (void)fprintf(stdout,"\nEquiv radius (microns) of crystals: size bin (across)\n");
    for(size=1;size<=num_size;size++){
      (void)fprintf(stdout,"%.1f\t",equiv_radius[size]*1.e6);
    } /* end loop over sizes */
  } /* end debug */
  
  if(debug == 17){
    (void)fprintf(stdout,"\nPrism radius (microns) of crystals: size bin (across)\n");
    for(size=1;size<=num_size;size++){
      (void)fprintf(stdout,"%.1f\t",prism_radius[size]*1.e6);
    } /* end loop over sizes */
  } /* end debug */
  
  if(debug == 18){
    (void)fprintf(stdout,"\nDmass/Dtime (ng/s) of crystals: size bin (across) vs layer (down)\n");
    for(size=1;size<=num_size;size++){
      (void)fprintf(stdout,"\t%i",size);
    } /* end loop over sizes */
    for(layer=1;layer<=num_layer;layer++){
      (void)fprintf(stdout,"\n%i\t",layer);
      for(size=1;size<=num_size;size++){
	(void)fprintf(stdout,"%.2f\t",dmass_dtime[layer][size]*1.e12);
      } /* end loop over sizes */
    } /* end loop over layers */
  } /* end debug */
  
  if(debug == 22){
    (void)fprintf(stdout,"%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s\n",
		  "bin","Mass","dM","Length","dL",
		  "r_equiv","S area","Hex L","Hex D","aspect",
		  "Hex r_s","Hex area");
    (void)fprintf(stdout,"%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s\n",
		  "index","ng","ng","microns","microns",
		  "microns","micr^2","microns","microns","ratio",
		  "microns","micr^2");
    for(size=1;size<=num_size;size++){
      (void)fprintf(stdout,"%9i%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f\n",
		    size,
		    crystal_mass[size]*NANOGRAMS_PER_KILOGRAM,
		    delta_mass[size]*NANOGRAMS_PER_KILOGRAM,
		    crystal_length[size]*MICRONS_PER_METER,
		    delta_length[size]*MICRONS_PER_METER,
		    equiv_radius[size]*MICRONS_PER_METER,
		    equiv_rad_squared[size]*4.*M_PI*
		    MICRONS_PER_METER*MICRONS_PER_METER,
		    hex_column_length[size]*MICRONS_PER_METER,
		    hex_column_diameter[size]*MICRONS_PER_METER,
		    hex_column_length[size]/hex_column_diameter[size],
		    hex_equiv_radius[size]*MICRONS_PER_METER,
		    hex_equiv_rad_squared[size]*4.*M_PI*
		    MICRONS_PER_METER*MICRONS_PER_METER);
    } /* end loop over sizes */
  } /* end debug */
  
  if(debug == 32){
    for(layer=1;layer<=num_layer;layer++){
      (void)fprintf(stdout,"effective_radius[%i]= %.1f microns\n",
		    layer,effective_radius[layer]*1.e6);
    } /* end loop over layers */
  } /* end debug */
  
  if(debug == 34){
    if(HETEROGENEOUS_NUCLEATION){
      for(layer=1;layer<=num_layer;layer++){
	if(env_temperature[layer] >= 253.15 && saturation_ice[layer] >= 1.){
	  float_foo=heterogeneous_nucleation_MDC92
	    (saturation_ice[layer],dt,tau_CCN_timescale);
	  (void)fprintf(stdout,"heterogeneous nucleation: %g crystals/m^3 into layer = %i at time_step = %i\n",float_foo,layer,time_step);	
	} /* end if */ 
      } /* end loop over layers */
    } /* end if */
    if(HOMOGENEOUS_NUCLEATION){
      for(layer=1;layer<=num_layer;layer++){
	(void)fprintf(stdout,"CCN_activated[%i] = %g m^-3\n",
		      layer,CCN_activated[layer]);
      } /* end loop over layers */
    } /* end if */
  } /* endif debug */
  
  if(debug == 35){
    for(layer=1;layer<=num_layer;layer++){
      (void)fprintf(stdout,"wind_speed[%i] = %.1f cm/s\n",
		    layer,wind_speed[layer]*CENTIMETERS_PER_METER);
    } /* end loop over layers */
  } /* end debug */
  
  if(debug == 36){
    (void)fprintf(stdout,"\nFall Speed (cm/s) of crystals: size bin (across) vs layer (down)\n");
    for(size=1;size<=num_size;size++){
      (void)fprintf(stdout,"\t%i",size);
    } /* end loop over sizes */
    for(layer=1;layer<=num_layer;layer++){
      (void)fprintf(stdout,"\n%i\t",layer);
      for(size=1;size<=num_size;size++){
	(void)fprintf(stdout,"%.1f\t",fall_speed[size][layer]*CENTIMETERS_PER_METER);
      } /* end loop over sizes */
    } /* end loop over layers */
  } /* end debug */
  
  if(debug == 38){
    fprintf(stdout,"band\tlambda\tn_real\tn_imag\n");
    for(band=1;band<=num_band;band++){
      fprintf(stdout,"%i\t%.3f\t%.3f\t%g\n",
	      band,wavelength[band]*MICRONS_PER_METER,
	      ice_real_idx[band],ice_imag_idx[band]);
    } /* end loop over bands */
  } /* end debug */ 
  
  write_ccm2_output
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
     num_layer);
  
  if(debug == 48){
    for(layer=1;layer<=num_layer;layer++){
      for(size=1;size<=num_size;size++){
	
	/* Compute the crystal surface temperature */ 
	crystal_temperature[layer][size]=
	  env_temperature[layer]*
	    (1.+
	     (latent_heat_sub*dmass_dtime[layer][size]+
	      crystal_heat_net[layer][size])/
	     (4.*M_PI*env_temperature[layer]*capacitance[size]*
	      thermal_conductivity[layer]*
	      conductivity_correction[layer][size]));
	
	(void)fprintf(stdout,"crystal_temp[%i][%i] = %.2f, env_temp = %.2f\n",
		      layer,size,crystal_temperature[layer][size],
		      env_temperature[layer]);
      } /* end loop over sizes */
    } /* end loop over layers */
  } /* end debug */
  
  if(debug == 50){
    (void)fprintf(stdout,"zbin\talt\tsbin\tradius\tIRQabs\tcSWheat\tcLWcool\tq_rad\tls*dmdt\n");
    (void)fprintf(stdout,"\tkm\t\tmicrons\tm^2\tuW\tuW\tuW\tuW\n");
    for(layer=cloud_base_idx;layer<=cloud_top_idx;layer++){
      for(size=1;size<=num_size;size++){
	(void)fprintf(stdout,"%i\t%.2f\t%i\t%.2f\t%.4f\t%.2f\t%.2f\t%.2f\t%.2f\n",
		      layer,
		      altitude[layer]*KILOMETERS_PER_METER,
		      size,
		      hex_equiv_radius[size]*MICRONS_PER_METER,
		      planck_avg_abs_x_sec[layer][size],
		      crystal_heat_SW[layer][size]*MICROWATTS_PER_WATT,
		      crystal_heat_LW[layer][size]*MICROWATTS_PER_WATT,
		      crystal_heat_net[layer][size]*MICROWATTS_PER_WATT,
		      latent_heat_sub*dmass_dtime[layer][size]*MICROWATTS_PER_WATT);
      } /* end loop over sizes */
    } /* end loop over layers */
  } /* end debug */
  
  if(debug == 54){
    (void)fprintf(stdout,"zbin\talt\tsbin\tlength\tD(z)\tf_M\tK(z)\tf_T\n");
    (void)fprintf(stdout,"\tkm\t\tmicrons\tm2/s\t\tJ/(msC)\t\n");
    for(layer=1;layer<=num_layer;layer++){
      for(size=1;size<=num_size;size++){
	(void)fprintf(stdout,"%i\t%.2f\t%i\t%.2f\t%7.1e\t%.3f\t%.4f\t%.3f\n",
		      layer,
		      altitude[layer]*KILOMETERS_PER_METER,
		      size,
		      crystal_length[size]*MICRONS_PER_METER,
		      vapor_diffusivity[layer],
		      diffusivity_correction[layer][size],
		      thermal_conductivity[layer],
		      conductivity_correction[layer][size]);
      } /* end loop over sizes */
    } /* end loop over layers */
  } /* end debug */
  
  if(debug == 55){
    (void)fprintf(stdout,"IWC and INC and r_e and activated haze by level:\n");
    (void)fprintf(stdout,"%10s%10s%14s%10s%14s\n",
		  "layer","IWC","INC","r_e","Act. haze");
    (void)fprintf(stdout,"%10s%10s%14s%10s%14s\n",
		  "index","mg/m^3","#/liter","microns","#/liter");
    for(layer=1;layer<=num_layer;layer++){
      (void)fprintf
	(stdout,"%10i%10.3f%14.3f%10.3f%14.3f\n",
	 layer,
	 IWC[layer]*MILLIGRAMS_PER_KILOGRAM,
	 total_layer_conc[layer]/LITERS_PER_CUBIC_METER,
	 effective_radius[layer]*MICRONS_PER_METER,
	 CCN_activated[layer]/LITERS_PER_CUBIC_METER);
    } /* end loop over layers */
    (void)fprintf(stdout,"\nMean IWC = %.4f g/m^3\n",GRAMS_PER_KILOGRAM*
		  IWP_of_time[time_step]/(num_cloudy_layers*dz));
    (void)fprintf(stdout,"Mean INC = %.0f per liter\n",
		  total_VEC(total_layer_conc+cloud_base_idx,num_cloudy_layers)/ 
		  (num_cloudy_layers*LITERS_PER_CUBIC_METER));
    (void)fprintf(stdout,"Mean eff. radius = %.1f microns\n",
		  cloud_effective_radius[time_step]*MICRONS_PER_METER);
  } /* end debug */

  if(debug == 58 && RADIATION){
    (void)fprintf(stdout,"[layer][size] = [%i][%i] thermo/rad. = %f\n",
		  layer,size,numerator_1/numerator_2);
  } /* end debug */
  
  if(debug == 59){
    (void)fprintf(stdout,"zbin\talt\tSW up\tSW dn\tLW up\tLW dn\tdFSW/dz\tdFLW/dz\tQSW\tQLW\tQnet\n");
    (void)fprintf(stdout,"\tkm\tW/m^2\tW/m^2\tW/m^2\tW/m^2\tmW/m^3\tmW/m^3\tK/hour\tK/hour\tK/hour\n");
    for(layer=1;layer<=num_layer;layer++){
      (void)fprintf(stdout,"%i\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",
		    layer,
		    altitude[layer]*KILOMETERS_PER_METER,
		    flux_up_SW[layer],
		    flux_down_SW[layer],
		    flux_up_LW[layer],
		    flux_down_LW[layer],
		    SW_flux_div[layer]*MILLIWATTS_PER_WATT,
		    LW_flux_div[layer]*MILLIWATTS_PER_WATT,
		    heating_rate_SW[layer]*SECONDS_PER_HOUR,
		    heating_rate_LW[layer]*SECONDS_PER_HOUR,
		    heating_rate_net[layer]*SECONDS_PER_HOUR);
    } /* end loop over layers */
  } /* end debug */
  
  if(debug == 60){
    (void)fprintf(stdout,"zbin\talt\tLH\tAH\tRH\tNet Heat\n");
    (void)fprintf(stdout,"\tkm\tK/hour\tK/hour\tK/hour\tK/hour\n");
    for(layer=1;layer<=num_layer;layer++){
      (void)fprintf(stdout,"%i\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t\n",
		    layer,
		    altitude[layer]*KILOMETERS_PER_METER,
		    latent_heating[layer]*SECONDS_PER_HOUR,
		    advective_heating[layer]*SECONDS_PER_HOUR,
		    heating_rate_net[layer]*SECONDS_PER_HOUR,
		    (latent_heating[layer]+advective_heating[layer]
		     +heating_rate_net[layer])*SECONDS_PER_HOUR);
    } /* end loop over layers */
  } /* end debug */
  
  if(debug == 61){
    for(layer=1;layer<=num_layer;layer++){
      potential_temp[layer]=env_temperature[layer]*
	pow(101300./env_pressure[layer],kappa);
    } /* end loop over layers */
    potential_temp[0]=potential_temp[1]-
      (potential_temp[2]-potential_temp[1]);
    potential_temp[num_layer+1]=potential_temp[num_layer]+
      (potential_temp[num_layer]-potential_temp[num_layer-1]);

    (void)fprintf(stdout,"zbin\talt\tTemp\tdT/dz\tTheta\tdTh/dz\tAdv\tRad\tLat\n");
    (void)fprintf(stdout,"\tkm\tK\tK/km\tK\tK/km\tK/hour\tK/hour\tK/hour\n");
    for(layer=1;layer<=num_layer;layer++){
      float_foo=(.5*(env_temperature[layer+1]-env_temperature[layer-1])/dz);
      float_foo2=(.5*(potential_temp[layer+1]-potential_temp[layer-1])/dz);
      (void)fprintf(stdout,"%i\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",
		    layer,
		    altitude[layer]/1000.,
		    env_temperature[layer],
		    float_foo*METERS_PER_KILOMETER,
		    potential_temp[layer],
		    float_foo2*METERS_PER_KILOMETER,
		    advective_heating[layer]*SECONDS_PER_HOUR,
		    heating_rate_net[layer]*SECONDS_PER_HOUR,
		    latent_heating[layer]*SECONDS_PER_HOUR);
    } /* end loop over layers */
  } /* endif debug */
  
  if(debug == 64){
    if(RADIATION){
      for(layer=1;layer<=num_layer;layer++){
	for(size=1;size<=num_size;size++){
	  float_foo=4.*M_PI*capacitance[size]*(saturation_ice[layer]-1.)/
	    (latent_heat_sub*crystal_heat_net[layer][size]/
	     (gas_const_vapor*env_temperature[layer]*env_temperature[layer]*
	      thermal_conductivity[layer]*conductivity_correction[layer][size]));
	  float_foo2=4.*M_PI*capacitance[size]*
	    (saturation_ice[layer]-1.)*latent_heat_sub/
	      (crystal_heat_net[layer][size]*env_temperature[layer]*gas_const_vapor/
	       (eqm_vap_ice[layer]*vapor_diffusivity[layer]*
		diffusivity_correction[layer][size]));
	} /* end loop over sizes */
      } /* end loop over layers */
    }else{
      float_foo=0.;
      float_foo2=0.;
    } /* end else */
    
    (void)fprintf(stdout,"zbin\talt\tsbin\tlength\tConc\talpha\tbeta\tLatent\tSW Rad\tLW Rad\tTot Rad\tH=L+R\n");
    (void)fprintf(stdout,"\tkm\t\tmicrons\t#/l\t\t\tuW\tuW\tuW\tuW\tuW\n");
    for(layer=1;layer<=num_layer;layer++){
      for(size=1;size<=num_size;size++){
	/* Compute the crystal surface temperature */ 
	(void)fprintf(stdout,"%i\t%.2f\t%i\t%.2f\t%.1f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",
		      layer,
		      altitude[layer]*KILOMETERS_PER_METER,
		      size,
		      crystal_length[size]*MICRONS_PER_METER,
		      concentration[layer][size]*CUBIC_METERS_PER_LITER,
		      float_foo,
		      float_foo2,
		      latent_heat_sub*dmass_dtime[layer][size]*MICROWATTS_PER_WATT,
		      crystal_heat_SW[layer][size]*MICROWATTS_PER_WATT,
		      crystal_heat_LW[layer][size]*MICROWATTS_PER_WATT,
		      crystal_heat_net[layer][size]*MICROWATTS_PER_WATT,
		      (latent_heat_sub*dmass_dtime[layer][size]+
		      crystal_heat_net[layer][size])*MICROWATTS_PER_WATT);
      } /* end loop over sizes */
    } /* end loop over layers */
  } /* end debug */

  if(debug == 66){
    (void)fprintf(stdout,"zbin\talt\tband\tlambda\tngSWdiv\tngLWdiv\tdFSW/dz\tdFLW/dz\n");
    (void)fprintf(stdout,"\tkm\t\tmicrons\tmW/m^3\tmW/m^3\tmW/m^3\tmW/m^3\n");
    for(layer=1;layer<=num_layer;layer++){
      float_foo=0.;
      float_foo2=.5*(flux_up_SW[layer+1]-flux_down_SW[layer+1]-
		     flux_up_SW[layer-1]+flux_down_SW[layer-1])/dz;
      float_foo3=.5*(flux_up_LW[layer+1]-flux_down_LW[layer+1]-
		     flux_up_LW[layer-1]+flux_down_LW[layer-1])/dz;
      for(band=1;band<=num_band;band++){
	float_foo+=SW_spec_flux_div[layer][band];
	(void)fprintf(stdout,"%i\t%.2f\t%i\t%.2f\t%.2f\t%.2f\t%7.1e\t%.2f\n",
		      layer,
		      altitude[layer]*KILOMETERS_PER_METER,
		      band,
		      wavelength[band]*MICRONS_PER_METER,
		      float_foo2*MILLIWATTS_PER_WATT,
		      float_foo3*MILLIWATTS_PER_WATT,
		      SW_spec_flux_div[layer][band]*MILLIWATTS_PER_WATT,
		      LW_flux_div[layer]*MILLIWATTS_PER_WATT);
      } /* end loop over bands */
      (void)fprintf(stdout,"total SW flux divergence = %.2f mW/m^3\n\n",
		    float_foo*MILLIWATTS_PER_WATT);
    } /* end loop over layers */
  } /* end debug */
  
  if(debug == 69){
    (void)fprintf(stdout,"layer\talt\tpres\tintpres\temis.\ttrans.\tbbodyF\tFintup\tFupCCM2\tFintdn\tFdnCCM2\tdFdz\tdFdzCCM\n");
    (void)fprintf(stdout,"\tkm\tmb\tmb\t\t\tW/m^2\tW/m^2\tW/m^2\tW/m^2\tW/m^2\tmW/m^3\tmW/m^3\n");
    for(layer=0;layer<=num_layer;layer++){
      (void)fprintf(stdout,"%i\t%.2f\t%.1f\t%.1f\t%.4f\t%.4f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",
		    layer,
		    altitude[layer]*KILOMETERS_PER_METER,
		    env_pressure[layer]*MB_PER_PASCAL,
		    env_pressure_int[layer]*MB_PER_PASCAL,
		    emissivity[layer],
		    transmissivity[layer],
		    blackbody_flux[layer],
		    flux_up_LW_int[layer],
		    flux_up_LW[layer],
		    flux_down_LW_int[layer],
		    flux_down_LW[layer],
		    LW_flux_div[layer]*MILLIWATTS_PER_WATT,
		    -heating_rate_LW[layer]*MILLIWATTS_PER_WATT*
		    env_density[layer]*spec_heat_air);
    } /* end loop over layers */
  } /* end debug */
  
  if(debug == 70){
    integrate_LW_abs_opt_depth(
			       LW_abs_opt_depth,
			       concentration,
			       dz,
			       planck_avg_abs_x_sec,
			       num_layer,
			       num_size);
    
    (void)fprintf(stdout,"\nTimestep = %i = %.1f minutes\n",time_step,
		  time_array[time_step]/60.);
    (void)fprintf(stdout,"zbin\talt\tVTauExt\tLWTauAbs\n");

    float_foo2=0.;
    for(layer=1;layer<=num_layer;layer++){
      float_foo2+=LW_abs_opt_depth[layer];
      (void)fprintf(stdout,"%i\t%.2f\t%.2f\t%.2f\n",
		    layer,
		    altitude[layer]*KILOMETERS_PER_METER,
		    optical_depth[layer],
		    LW_abs_opt_depth[layer]);
    } /* end loop over layers */
    (void)fprintf(stdout,"Total visible (525 nm) extinction optical depth = %f\n",
		  cloud_optical_depth[time_step-1]);
    (void)fprintf(stdout,"Total Planckian Avg. LW absorption optical depth = %f\n",
		  float_foo2);
  } /* end debug */
  
  if(debug == 72){
    float_foo3=0.;
    for(layer=cloud_base_idx;layer<=cloud_top_idx;layer++){
      (void)fprintf(stdout,"zbin\tAlt\twbin\tlambda\tLWsbw\tLWsvac\tLWsmac\tTotLWmac\n");
      (void)fprintf(stdout,"\tkm\t\tmicrons\t\tkm^-1\tm^2/g\tm^2/g\n");
      for(band=1;band<=num_LW_band;band++){
	(void)fprintf(stdout,"%i\t%.2f\t%i\t%.2f\t%.4f\t%7.1e\t%.4f\t%8.2e\n",
		      layer,
		      altitude[layer]*KILOMETERS_PER_METER,
		      band,
		      IR_wavelength[band]*MICRONS_PER_METER,
		      LW_spec_band_weight[layer][band],
		      LW_spec_vol_abs_coeff[layer][band]*METERS_PER_KILOMETER,
		      (IWC[layer] > EMINUSTWENTY) ? 
		      LW_spec_vol_abs_coeff[layer][band]*
		      KILOGRAMS_PER_GRAM/IWC[layer] : 0.,
		      LW_mass_abs_coeff[layer]*KILOGRAMS_PER_GRAM);
      } /* end loop over bands */
      (void)fprintf(stdout,"Total of %i LW band weights = %g\n",num_LW_band,
		    total_VEC(LW_spec_band_weight[layer]+1,num_LW_band));
      float_foo=LW_DIFFUSIVITY*LW_mass_abs_coeff[layer]*IWP[layer];
      float_foo2=1.-exp(-float_foo);
      (void)fprintf(stdout,"kappa total = %8.1e m^2/g, IWP = %8.1e g/m^2\n",
		    LW_mass_abs_coeff[layer]*KILOGRAMS_PER_GRAM,
		    IWP[layer]*GRAMS_PER_KILOGRAM);
      (void)fprintf(stdout,"sigma total = %.4f km^-1\n",
		    LW_mass_abs_coeff[layer]*IWC[layer]*
		    METERS_PER_KILOMETER);
      (void)fprintf(stdout,"beta * kappa total * IWP = %.4f, emissivity = %.4f\n\n",
		    float_foo,float_foo2);
      /* accumulate the sigma's: this counts all crystals in the cloud */ 
      float_foo3+=LW_mass_abs_coeff[layer]*IWC[layer];
    } /* end loop over layers */;
    /* normalize the total sigma to represent crystals in a cubic meter */ 
    float_foo3/=num_cloudy_layers;
    float_foo=IWP_of_time[time_step]/(num_cloudy_layers*dz);
    (void)fprintf(stdout,"Bulk Parameters:\n");
    (void)fprintf(stdout,"IWP = %.4f g/m^2, Delta z = %.2f km\n",
		  IWP_of_time[time_step]*GRAMS_PER_KILOGRAM,
		  (num_cloudy_layers*dz)*KILOMETERS_PER_METER);
    (void)fprintf(stdout,"mean IWC = IWP/Delta z = %.4f mg/m^3\n",
		  float_foo*MILLIGRAMS_PER_KILOGRAM);
    (void)fprintf(stdout,"LW volume abs coeff sigma = %.4f km^-1\n",
		  float_foo3*METERS_PER_KILOMETER);
    (void)fprintf(stdout,"LW mass abs coeff kappa = %.4f m^2/g\n",
		  (float_foo3/float_foo)*KILOGRAMS_PER_GRAM);
    float_foo=LW_DIFFUSIVITY*float_foo3*dz*num_cloudy_layers;
    float_foo2=1.-exp(-float_foo);
    (void)fprintf(stdout,"emissivity epsilon = %.4f\n\n",float_foo2);
  } /* end debug */

  if(debug == 73){
    max_FLOATVEC_idxed(IWC+cloud_base_idx,num_cloudy_layers,
			 &float_foo,&int_foo2);
    int_foo=int_foo2+cloud_base_idx;
    (void)fprintf(stdout,"Thickness of cloud = %.2f km\n",
		  (altitude[cloud_top_idx]-altitude[cloud_base_idx])*
		  KILOMETERS_PER_METER);
    (void)fprintf
      (stdout,"max(IWC) in cloud = %.3f mg/m^3 @ layer %i (cbi = %i, cti = %i)\n",
       float_foo*MILLIGRAMS_PER_KILOGRAM,int_foo,cloud_base_idx,cloud_top_idx); 

    (void)fprintf(stdout,"Shortwave Albedo = %.3f\n",albedo_of_time[time_step-1]);
    (void)fprintf
      (stdout,"Shortwave Absorptance = %.3f\n",
       ((flux_up_SW[cloud_top_idx]-flux_down_SW[cloud_top_idx])-
	(flux_up_SW[cloud_base_idx-1]-flux_down_SW[cloud_base_idx-1]))/
       (flux_up_SW[cloud_top_idx]-flux_down_SW[cloud_top_idx]));
    (void)fprintf(stdout,"Shortwave Transmittance = %.3f\n",
		  flux_down_SW[cloud_base_idx-1]/flux_down_SW[cloud_top_idx]);
    (void)fprintf
      (stdout,"emissivity (method):\n\t %.3f (cz bulk cloud -D 72)\n\t %.3f (stephens 90)\n\t %.3f (stackhouse 91)\n",
       emissivity_of_time[time_step-1],
       flux_down_LW[cloud_base_idx]/blackbody_flux[int_foo],
       (flux_up_LW[cloud_base_idx-1]-flux_up_LW[cloud_top_idx])/
       (flux_up_LW[cloud_base_idx-1]-
	blackbody_flux[(int_foo)]));
  } /* end debug */
  
  if(debug == 76){
    (void)fprintf(stdout,"X-sections and individual optical properties of crystals:\n");
    for(band=1;band<=num_band;band++){
      (void)fprintf
	(stdout,"%10s%10s%10s%10s%12s%12s%13s%10s%10s\n",
	 "wavebin","wavelength","sizebin","hex length","ext-xs",
	 "scat-xs","abs-xs","co-albedo","asym param");
      (void)fprintf(stdout,"%10s%10s%10s%10s%12s%12s%13s%10s%10s\n",
		    "index","microns","index","microns","micron^2",
		    "micron^2","micron^2","","");
      for(size=1;size<=num_size;size++){
	(void)fprintf
	  (stdout,"%10i%10.3f%10i%10.3f%12.3f%12.3f%13.5f%10.6f%10.3f\n",
	   band,
	   wavelength[band]*MICRONS_PER_METER,
	   size,
	   hex_column_length[size]*MICRONS_PER_METER,
	   extinction_x_sec[band][size]*SQUARE_MICRONS_PER_SQUARE_METER,
	   scattering_x_sec[band][size]*SQUARE_MICRONS_PER_SQUARE_METER,
	   absorption_x_sec[band][size]*SQUARE_MICRONS_PER_SQUARE_METER,
	   1.-scattering_x_sec[band][size]/extinction_x_sec[band][size],
	   asymmetry_param[band][size]);
      } /* end loop over sizes */
    } /* end loop over bands */
  } /* end debug */
  
  /* Find out if there are any outstanding exceptions */ 
#if ( defined sun )
  ieee_retrospective(stdout);
#endif

  /* Clear all signal handling traps before exiting */ 
#if ( defined sun )
  ieee_handler("clear","common",(sigfpe_handler_type)ieee_exception_handler);
#endif

  Exit_gracefully();
} /* end of main() */

#if ( defined SunOS )
/* See the Sun Numerical Computation Guide p. 145 */ 
void ieee_exception_handler
(int signal,
 int code,
 struct sigcontext *scp, 
 char *address)
#endif
#if ( defined Solaris )
void ieee_exception_handler
(int signal,
 siginfo_t *sip, 
 ucontext_t *uap)
#endif
#if ( defined AIX )
void ieee_exception_handler()
#endif
#if ( defined IRIX )
void ieee_exception_handler()
#endif
#if ( ! defined UNICOS )
{
  /* See /usr/include/signal.h for SIGFPE codes.
     This example is taken from the Sun Numerical Computation Guide, p. 146 */ 
/*  (void)fprintf(stdout,"floating point exception %x at address %x\n",*/
/*                code,*/
/*                address);*/
} /* end ieee_exception_handler() */ 
#endif

void print_debug(void)
{
#define NUM_DEBUG_LEVELS 80
  
  int level;
  static char *debug_action[]={
    /* Debug options 0--10 */ 
    "Do nothing. The default debugging mode.",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "Environmental Profiles: prints k,z (km),T (K), T (C), S, P (mb), r_vap (g/kg), r_ice (g/kg)",
    /* Debug options 11--20 */ 
    "Initial concentration of entire distribtuion in #/liter",
    "Time evolution of concentration of one particular layer & size",
    "Thermal conductivity of each layer",
    "",
    "daltitude_dtime (cm/s) of each layer & size",
    "equiv_radius (microns) of each size",
    "prism_radius (microns) of each size",
    "dmass_dtime (ng/s) of each layer & size at end",
    "Concentration (#/m^3) of each layer & size & timestep",
    "d_env_temp_dt of each layer & timestep",
    /* Debug options 21--30 */ 
    "vapor_diffusivity of each layer at end",
    "Crystal dimensions: j, Mass (ng),dMass (ng), L (u), dL (u), r_s (u), S (u^2), Hex L (u), Hex D (u), Aspect, Hex r_s (u), Hex S (u^2)",
    "Euler_stable_size of each layer & timestep",
    "Length (microns) of Euler_stable_size of each layer & timestep",
    "Concentration (#/liter) of one particular layer & size at end",
    "No opt",
    "No opt",
    "IWC (kg/m^3) at each layer & timestep post-Lagrangian but pre-Eulerian growth",
    "Concentration (#/liter) of each layer & size at end",
    "Tells how much sublimate (ng) (if any) from each layer & size & timestep",
    /* Debug options 31--40 */ 
    "Maximum Euler_stable_size in cloud, cloud_base, and cloud_top at each timestep, ",
    "Effective radius of each layer at end",
    "No opt",
    "Nucleation (#/m^3) (if any) into each layer at end",
    "wind_speed (cm/s) of each layer at end",
    "fall_speed (cm/s) of each layer & size",
    "Vapor path conservation check before and after wind transport",
    "SW band, lambda (microns), ice real index, ice imag index",
    "Mie_Bohren_c internal debugging info",
    "Mie optical properties table: band, lambda (u), size, r_s (u), n_real, n_imag, X, Q_ext, Q_sca, Q_abs, g, omega",
    /* Debug options 41--50 */ 
    "Runs Bohren's test case Mie parameters",
    "CCM2 environmental input, e.g., from ccm2rad.trp35.dat",
    "Bulk layer optical properties of each layer & band & timestep: cbi, z (km), cti, z (km), band, lambda (u), eff_tau_sca, eff_tau_ext, eff_omega, eff_g",
    "Cloud grid info: cbi, cti, num_cloud_layers, z (km), bccmil, tccmil, nccmlev",
    "Outputs full colmod style info each timestep",
    "No opt",
    "Ebert & Curry treatment of SW & LW optical properties by bulk IWC and r_e",
    "Crystal surface temperature and environmental temperature at each layer & size at end",
    "No opt",
    "Crystal Power absorption at each layer & size at end: hex_r_e (u), planck_q_abs_avg, q_SW (uW), q_LW (uW), latent heat (uW)",
    /* Debug options 51--60 */ 
    "No opt",
    "IR wavelength and band diagnostics at each layer and band: lambda (u), dlambda (u), T (K), band weight, IR radiance, total intens., sigma T^4",
    "Slingo treatment of all bulk IWC as LWC and use r_e (option crashes)",
    "Diffusivities of each layer & size: z (km), L (u), D_M, f_M, K_T, f_T",
    "Bulk distribution properties of each layer at end: IWC (mg), N (#/l), r_e (u), mean IWC, mean INC, mean r_e",
    "No opt",
    "Initial cloud lookup environmental table at grid points and interpolated state variables at each layer: z (km), T (K), p (mb)",
    "Latent/Radiative heating of one particular size & layer at end",
    "Environmental fluxes and heating rates: z (km), FupSW (W/m^2), FdownSW, FupLW, FdownLW, SWFdiv (mW/m^3), LWFdiv, HSW (K/hr), HLW, Hnet",
    "Heat balance of each layer at end: z(km), latent (K/hr), advective, net radiative, their sum",
    /* Debug options 61--70 */ 
    "Stability Check: z (km), T (K), dT/dz, Theta (K), dTheta/dz, advective (K/hr), radiative net, latent",
    "Reset LWmac kappa to its original CCM2 value of 60 m^2/kg",
    "Latent/Radiative heating of each layer & size & timestep: z (km), L (u), N (#/l), alpha, beta, latent (uW/m^3), q_SW, q_LW, q_net, latent + net",
    "Latent/Radiative heating of each layer & size at end: z (km), L (u), N (#/l), alpha, beta, latent (uW/m^3), q_SW, q_LW, q_net, latent + net",
    "Comparison of q apportioned equally by crystal conc. vs. weighted by size & band at each layer & size & timestep: z (km), L (u), simple q_SW, q_LW, q_net, actual q_SW, q_LW, q_net",
    "Comparison of flux divergences computed by differencing returned total fluxes vs. summing their spectral totals at end: z (km), lambda (u), ngSWdiv, ngLWdiv, SWdiv, LSdiv, total SWdiv",
    "Check sum comparison of total radiative heating vs. computed flux divergences of each layer & timestep: z (km), total q_SW, SWdiv, total q_LW, LWdiv",
    "No opt",
    "LW validation of bulk (interface) method vs. CCM2 of each layer at end: z (km), p (mb), p_int (mb), epsilon, tranmissivity, sigma T^4, FupLW_int, FupLW, FdnLW_int, FdnLW, FdivLW, FdivLW_int",
    "Display visible extinction and LW absorption optical depths of each layer at end",
    /* Debug options 71--80 */ 
    "Concentration = 0 initially in each layer & size",
    "LW mass absorption validation at each layer & band at end: z (km), lambda (u), LW spec band weight, LW spec vol abs coeff (km^-1), LWsvac/IWC (m^2/g), LWmac (m^2/g). Then layer totals: weights, kappa, LWsvac, beta*kappa*IWP, epsilon. Then cloud-integrated totals: IWP, mean IWC, LWvac, LWmac, total epsilon",
    "Albedo-Emissivity diagnostics at end: thickness, max(IWC) and its height, albedo, absorptance, transmittance, emissivity by zender, stephens, stackhouse methods",
    "Apply Mie_Liou correction fudge factors in the SW",
    "Output from Liou_fudge hex column crystal dimensions of each size L (u), D (u), Aspect, r_s (u), r_s^2 (u^2), A  (u^2), x-Sec A (u^2), Mass (ng) and bilinear_interp matrices (when applicable)",
    "X-sections and individual optical properties of crystals (whether hexagonal or not",
    "Hex column dimension diagnostics from newton-raphson routine: initial mass (kg), volume (m^3), L (u), iteration, L (u), epsilon, final D (u), aspect",
    "Band-averaged Mie calcs. at each quad. point: lambda (u), n_real, n_imag, avg_n_real, avg_n_imag",
    "SWCF and surface SWCF diagnostics",
    ""
    };

  (void)fprintf(stderr,"There are about %i debugging options:\n\n",
		NUM_DEBUG_LEVELS-10);
  (void)fprintf(stderr,"Option\tAction\n");
  (void)fprintf(stderr,"%i: %s\n",0,debug_action[0]);
  for(level=10;level<=NUM_DEBUG_LEVELS;level++){
    (void)fprintf(stderr,"%i: %s\n",level,debug_action[level]);
  } /* end loop over levels */
} /* end of print_debug() */ 

void print_usage(char *option_string)
{
  (void)fprintf(stderr,
		"\nusage: cld [-options] where options are one or more of:\n\n");
  (void)fprintf(stderr,"%s\n\n",option_string);
  (void)fprintf(stderr,"-A LAGRANGE Toggle Lagrangian growth. Default is ON\n");
  (void)fprintf(stderr,"-B EULER Toggle Eulerian growth phase. Default is ON\n");
  (void)fprintf(stderr,"-C ICE_ADVECT Toggle velocity advection. Default is ON\n");
  (void)fprintf(stderr,"-D debug Set the debugging level.  Default is 0\n");
  (void)fprintf(stderr,"-F HETEROGENEOUS_NUCLEATION Toggle. Default is True\n");
  (void)fprintf(stderr,"-G grid_base Default is 8000. m\n");
  (void)fprintf(stderr,"-H cloud_base Default is 1.e4 m\n");
  (void)fprintf(stderr,"-J ice_crystal_habit. Default is HEXAGONAL_COLUMNS\n");
  (void)fprintf(stderr,"-K AGGREGATION Toggle aggregation. Default is ON\n");
  (void)fprintf(stderr,"-L HOMOGENEOUS_NUCLEATION Toggle. Default is True\n");
  (void)fprintf(stderr,"-M atmospheric_profile_type Default is TROPICS\n");
  (void)fprintf(stderr,"-N INC_init Default is 1.e6 m^-3\n");
  (void)fprintf(stderr,"-P PAUSE Wait in between graphs? Default is True\n");
  (void)fprintf(stderr,"-Q RAD_FEEDBACK toggle flux feedback. Default is True\n");
  (void)fprintf(stderr,"-R RADIATION toggle radiation routines. Default is True\n");
  (void)fprintf(stderr,"-S S_ice_cloud Default is 1.2 ice sat ratio in cloud\n");
  (void)fprintf(stderr,"-U VAPOR_ADVECT toggle. Default is True\n");
  (void)fprintf(stderr,"-V VERBOSE toggle printing of WARNINGS. Default is True\n");
  (void)fprintf(stderr,"-X LIOU toggle hex crystals. Default is False\n");
  (void)fprintf(stderr,"-W WEED_NEGATIVE_VALUES toggle. Default is True\n");
  (void)fprintf(stderr,"-Z initial_distribution_type option. Default is KNOLLENBERG\n");
  (void)fprintf(stderr,"-a alpha_growth Courant condition safety factor. Default is .4\n");
  (void)fprintf(stderr,"-b bins_per_doubling Default is 2\n");
  (void)fprintf(stderr,"-c num_ccm2_level Default is 35\n");
  (void)fprintf(stderr,"-d debug_value Default is 0\n");
  (void)fprintf(stderr,"-e err_file Get the error file name. Default is stderr\n");
  (void)fprintf(stderr,"-f float_foo_input Generic tuning param. Default is 0.\n");
  (void)fprintf(stderr,"-g grid_thick Default is 4000. m\n");
  (void)fprintf(stderr,"-h cloud_thick Default is 1000. m\n");
  (void)fprintf(stderr,"-i in_file get the input file name. Default is stdin\n");
  (void)fprintf(stderr,"-j alpha_transport Courant condition safety factor. Default is .2\n");
  (void)fprintf(stderr,"-k dt timestep (Courant overrides).  Default is 100. s\n");
  (void)fprintf(stderr,"-l num_layer Default is 100\n");
  (void)fprintf(stderr,"-m smallest_length Default is 10.e-6\n");
  (void)fprintf(stderr,"-n num_step Default is 1\n");
  (void)fprintf(stderr,"-o out_file Output file name. Default is cloud.nc\n");
  (void)fprintf(stderr,"-p plot_step Default is 20 time steps per snapshot\n");
  (void)fprintf(stderr,"-q print out the debugging actions\n");
  (void)fprintf(stderr,"-r rad_step radiation update interval. Default is 1\n");
  (void)fprintf(stderr,"-s S_liq_ambient Default is .4 liq sat ratio beneath cloud\n");
  (void)fprintf(stderr,"-t tau_CCN_timescale Default is 1. s\n");
  (void)fprintf(stderr,"-u UPDRAFT_PROFILE_IS_FLAT toggle. Default is False\n");
  (void)fprintf(stderr,"-v print the RCS program version\n");
  (void)fprintf(stderr,"-w wind_speed_max Default is .05 m/s\n");
  (void)fprintf(stderr,"-x num_size Default is 40\n");
  (void)fprintf(stderr,"-z MOVIE_NETCDF_OUTPUT Default is False\n");
  (void)fprintf(stderr,"\n");
} /* end print_usage() */ 

void Exit_gracefully(void)
{
  char *time_buf_finish;
  time_t clock;
  
  /* end the clock */  
  
  clock=time((time_t *)NULL);
  time_buf_finish=ctime(&clock);
  (void)fprintf(stdout,"\tfinish = %s\n",time_buf_finish);
  
  exit(0);
} /* end Exit_gracefully() */ 
