/* $Author: zender $
 * $Date$
 * $Id$
 * $Revision$
 * $Locker:  $
 * $RCSfile: utilities.c,v $
 * $Source: /home/zender/cvs/cld/utilities.c,v $
 * $Id$
 * $State: Exp $
 * */

/* Purpose: Some generic and some cloud-specific utility routines. */ 

/* $Log: not supported by cvs2svn $
/* Revision 1.1.1.1  1998-09-15 02:06:44  zender
/* Imported sources
/*
 * Revision 5.4  1995/06/19  04:54:58  zender
 * have made many changes. got rid of fp_in,out,err. added flat updraft
 * switch. changed -S -s to both refer to ice. got rid of while loop that
 * became infinite when there was no cloud. changed vertical netCDF coord.
 * to altitude. fixed bug in appending date and rev. to cmdline.
 *
 * Revision 5.3  1995/05/24  23:11:17  zender
 * updated and validated the clouds model for ansi compatability,
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
 * different clouds. Made Knollenberg fit his observations, prepared
 * optical properties for 50 sizes and removed all num_size == 40 specific
 * stuff. changed tau_CCN to 1. s.
 *
 * Revision 4.8  1993/05/29  00:01:13  zender
 * fixed the Liou_fudge() bug, fixed a wind_speed bug, added
 * Dowling and Radke's crystal distribution function as default.
 *
 * Revision 4.7  1993/05/27  16:08:24  zender
 * A synchronization check-in for all clouds programs.
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
 *  */

/* standard header files */
#include <stdio.h>              /* stderr, FILE, NULL, etc. */
#include <time.h>               /* machine time */
#include <math.h>               /* sin cos cos sin 3.14159 */

/* my header files */
#include "defs.h"               /* YES, NO, FILESIZE, and more */
#include "globals.h"            /* externs announced in main() */

/* determines threshold for cloud base and top */ 
#define MIN_IWC_FOR_CLOUD .1e-6 /* kg/m^3 */         
#define MIN_CONC_FOR_CLOUD 50. /* #/m^3 */         

void Kinetic_Ventilation(conductivity_correction,
			 diffusivity_correction,
			 thermal_conductivity,vapor_diffusivity,
			 daltitude_dtime,
			 env_temperature,env_density,
			 crystal_length,prism_radius,equiv_rad_squared,
			 capacitance,
			 num_layer,num_size)
     float **conductivity_correction;
     float **diffusivity_correction;
     float **daltitude_dtime;
     float *env_temperature,*env_density;
     float *crystal_length,*prism_radius,*equiv_rad_squared;
     float *thermal_conductivity;
     float *vapor_diffusivity;
     float *capacitance;
     int num_layer,num_size;
{
  /* Figures out the mechanical corrections to the vapor diffusivity of air 
     as per Hall and Pruppacher '76 czr54 */ 

  /* nomenclature conversion between the routine and Hall & Pruppacher:
     conductivity_ventilation_correction = f_sub_q
     diffusivity_ventilation_correction = f_sub_m
     conductivity_kinetic_correction = F_sub_alpha
     diffusivity_kinetic_correction = F_sub_beta */ 

  float Prandtl_number;
  float Reynolds_number;
  float Schmidt_number;
  float alpha;
  float beta; 
  float collision_speed;
  float conductivity_kinetic_correction;
  float conductivity_ventilation_correction;
  float diffusivity_kinetic_correction;
  float diffusivity_ventilation_correction;
  float dynamic_viscosity;
  float l_star;
  float l_sub_m_star;
  float l_sub_m_star_factor;
  float l_sub_q_star;
  float l_sub_q_star_factor;
  float number_product;
  float perimeter;
  float r_star;
  float surface_area;

  int layer;
  int size;

  /* alpha is the fraction of air molecules which achieve thermal
     equilibrium while striking an ice surface. */
  alpha = 1.; 

  /* beta is the ratio of the actual number of water molecules leaving 
     the ice surface during evaporation to the number predicted by kinetic
     theory on the basis of the ice surface's equilibrium vapor pressure. */ 
  beta = .1; 

  for(layer=1;layer<=num_layer;layer++){
    dynamic_viscosity=1.832e-4* /* formula for mu is in cgs = cm/g/s */ 
      pow(env_temperature[layer]/296.16,1.5)*
	(296.16+120.)/(env_temperature[layer]+120.);
    dynamic_viscosity*=.1; /* convert from cm/g/s --> m/kg/s */ 

    Prandtl_number=dynamic_viscosity*spec_heat_air/thermal_conductivity[layer];

    Schmidt_number=dynamic_viscosity/(env_density[layer]*vapor_diffusivity[layer]);

    collision_speed=sqrt(8.*gas_const_dry_air*env_temperature[layer]/M_PI);

    l_sub_m_star_factor=sqrt(2.*M_PI/(gas_const_vapor*env_temperature[layer]))*
      vapor_diffusivity[layer]/(2.*beta/(2.-beta));

    l_sub_q_star_factor=4./(env_density[layer]*collision_speed*alpha*
			    spec_heat_air);

    for(size=1;size<=num_size;size++){
      perimeter=2.*(prism_radius[size]+crystal_length[size]);
      surface_area=4.*M_PI*equiv_rad_squared[size];
      l_star=surface_area/perimeter;
      r_star=equiv_rad_squared[size]/capacitance[size];
      
      Reynolds_number=fabs(daltitude_dtime[size][layer])*
	l_star*env_density[layer]/dynamic_viscosity;
	
      number_product=pow(Schmidt_number,1./3.)*sqrt(Reynolds_number);
      if(number_product <= 1.){
	diffusivity_ventilation_correction=1.+.14*number_product*number_product;
      }else{
	diffusivity_ventilation_correction=.86+.28*number_product;
      } /* end else */

      number_product=pow(Prandtl_number,1./3.)*sqrt(Reynolds_number);
      if(number_product <= 1.){
	conductivity_ventilation_correction=1.+.14*number_product*number_product;
      }else{
	conductivity_ventilation_correction=.86+.28*number_product;
      } /* end else */

      l_sub_m_star=diffusivity_ventilation_correction*l_sub_m_star_factor;

      l_sub_q_star=thermal_conductivity[layer]*l_sub_q_star_factor*
	conductivity_ventilation_correction;
	
      /* Get the total kinetic correction factors */ 
      diffusivity_kinetic_correction=r_star/(r_star+l_sub_m_star);
      conductivity_kinetic_correction=r_star/(r_star+l_sub_q_star);

      /* Multiply the kinetic correction by the ventilation correction to
       get the total kinetic-ventilation correction */ 
      diffusivity_correction[layer][size]=diffusivity_kinetic_correction*
	diffusivity_ventilation_correction;
      conductivity_correction[layer][size]=conductivity_kinetic_correction*
	conductivity_ventilation_correction;
      
    } /* end loop over sizes */
  } /* end loop over layers */
} /* end Kinetic_Ventilation() */

float Planck_function(temperature,wavelength)
     float temperature,wavelength;
{
  /* This function returns the specific intensity of radiation emitted
   by a black body in SI/MKS i.e. J/m^2/s/m/ster */ 
  
  float denom;
  float exponent;
  float specific_intensity;

  /* beware of overflows */ 
  if((exponent=(hcfirst/(wavelength*boltzmann_constant*temperature))) > 88.){
    denom=1.e38;
  }else{
    denom=exp(exponent)-1.;
  } /* end else */ 

  specific_intensity=2.*hcsquared/(pow(wavelength,5.)*denom);

  return specific_intensity;
} /* end Planck_function */ 

float d2r(theta_degrees)
     float theta_degrees;
{
  return theta_degrees*M_PI/180.;
} /* end d2r */ 

float r2d(theta_radians)
     float theta_radians;
{
  return theta_radians*180./M_PI;
} /* end r2d */ 

float Effective_Radius(concentration,equiv_rad_squared,
		       equiv_rad_cubed,num_size)
     float *concentration,*equiv_rad_squared,*equiv_rad_cubed;
     int num_size;
{
  /* Input: concentration,equiv_rad_squared,equiv_rad_cubed,num_size
     Output: answer */
  /* Note that the argument arrays are assumed to be indexed [1..num_size] */

  float dot_product_VEC();
  float max_VEC();

  float effective_radius;
  float second_moment_of_dist;
  float third_moment_of_dist;

  /* Trapezoidal rule is used for integration */
  second_moment_of_dist = 
    dot_product_VEC(concentration+1,equiv_rad_squared+1,num_size);
  second_moment_of_dist-=.5*
    (concentration[1]*equiv_rad_squared[1]+
     concentration[num_size]*equiv_rad_squared[num_size]);

  third_moment_of_dist = 
    dot_product_VEC(concentration+1,equiv_rad_cubed+1,num_size);
  third_moment_of_dist-=.5*
    (concentration[1]*equiv_rad_cubed[1]+
     concentration[num_size]*equiv_rad_cubed[num_size]);

  if((second_moment_of_dist > 0.) && (third_moment_of_dist > 0.)){
    effective_radius=third_moment_of_dist/second_moment_of_dist;
  }else{
    (void)fprintf(stdout,"2nd, 3rd moments of dist are %g, %g\n",
		  second_moment_of_dist,third_moment_of_dist);
  } /* end else */
    
  return effective_radius; /* m */
}

float homogeneous_nucleation_DMC94
(float *CCN_conc,
 float *CCN_diameter,
 float temperature,
 int num_CCN_size,
 float saturation_liquid)
{
  /* Perform the homogeneous freezing of deliquesent but unactivated
     aerosol particles.  These are ammonium sulphate solution particles. 
     First assumption: all layers have the same distribution of dry sulphate 
     mass. */ 
  /* Nota Bene: the DMC pzns. are written for CGS units. This pzn appears
     in DMC94 p. 80. */ 

  float total_VEC(float *,int);

  static float *CCN_frac_freezing;
  static float *CCN_conc_freezing;

  float CCN_density_cgs;
  float CCN_diam_cgs;
  float a;
  float b;
  float c1;
  float c2;
  float float_foo;
  float solute_density;
  float t_celsius;
  float tot_CCN_freezing;

  int size;

  if(CCN_frac_freezing == NULL){
    if(
       ((CCN_conc_freezing=(float *)malloc((num_CCN_size+2)*sizeof(float))) == NULL ) ||
       ((CCN_frac_freezing=(float *)malloc((num_CCN_size+2)*sizeof(float))) == NULL ) ||
       False ){
      (void)fprintf(stdout,"Unable to allocate array in main\n");
      exit(1);
    } /* end if */
  } /* end if */

  solute_density=ammonium_sulfate_density;
  t_celsius=temperature-triple_point_of_water;
  CCN_density_cgs=solute_density*GRAMS_PER_KILOGRAM/CUBIC_CMS_PER_CUBIC_METER;

  c1 = -14.65 - (1.045*t_celsius);
  c2 = -492.35- (8.34*t_celsius) - (.0061*t_celsius*t_celsius);
  float_foo = c1 + c2*(1.-saturation_liquid);
  float_foo = exp(float_foo*M_LN10);
  a = M_PI*M_PI*CCN_density_cgs*CCN_density_cgs*float_foo/6.;
  b = 6.;

  for(size=1;size<=num_CCN_size;size++){
    CCN_diam_cgs = CCN_diameter[size]*CENTIMETERS_PER_METER;
    float_foo = -a*pow(CCN_diam_cgs,b);
    CCN_frac_freezing[size] = 1.-exp(float_foo);
    CCN_conc_freezing[size] = CCN_frac_freezing[size]*CCN_conc[size];
  } /* end loop over sizes */

  tot_CCN_freezing=total_VEC(CCN_conc_freezing+1,num_CCN_size);
  return tot_CCN_freezing;

} /* end homogeneous_nucleation_DMC94() */ 

float CCN_PDF_DMC94
(float CCN_diameter,
 float scaling_diameter)
{
  /* Implements the DMC94 gamma distribution. For an understanding of
     what CSU types mean by all their distribution terminology see the
     RAMS microphysics documentation, FTVC89 */ 
  float CCN_PDF;

  CCN_PDF=(1./scaling_diameter)*exp(-CCN_diameter/scaling_diameter);
  return CCN_PDF;
} /* end CCN_PDF_DMC94 */ 

void init_CCN_conc
(float *CCN_mass, /* returned mass grid */ 
 float *CCN_diameter, /* returned length grid */ 
 float *CCN_conc, /* returned concentration */ 
 float CCN_tot_conc, /* input total concentration */ 
 float CCN_scaling_diameter, /* input mean diameter */ 
 int num_CCN_size)
{
  /* Return the concentration and mass and length grids for the dry solute 
     mass distribution. */ 

  float *delta_diameter;

  float CCN_PDF;
  float CCN_diameter_LBC;
  float CCN_diameter_RBC;
  float CCN_mass_LBC;
  float CCN_mass_RBC;
  float max_mass;
  float min_mass;
  float solute_density;

  int size;
  
  if(
     ((delta_diameter=(float *)malloc((num_CCN_size+2)*sizeof(float))) == NULL ) ||
     False ){
    (void)fprintf(stdout,"Unable to allocate array in init_CCN_conc()\n");
    exit(1);
  } /* end if */

  /* Initialize the free parameters */ 
  max_mass = 1.e-15; /* kg */ 
  min_mass = max_mass/pow(2.,((num_CCN_size-1.)/2.));
  solute_density = ammonium_sulfate_density;
  
  for(size=1;size<=num_CCN_size;size++){
    CCN_mass[size]=min_mass*pow(2.,((size-1.)/2.)); 
  } /* end loop over sizes */ 
  CCN_mass_LBC=CCN_mass[1]*(1.+1./sqrt(2.))/2.;
  CCN_mass_RBC=CCN_mass[num_CCN_size]*(1.+sqrt(2.))/2.;
  
  for(size=1;size<=num_CCN_size;size++){
    CCN_diameter[size]=6.*CCN_mass[size]/(M_PI*solute_density);
    CCN_diameter[size]=pow(CCN_diameter[size],(1./3.));
  } /* end loop over sizes */ 
  CCN_diameter_LBC=6.*CCN_mass_LBC/(M_PI*solute_density);
  CCN_diameter_LBC=pow(CCN_diameter_LBC,(1./3.));
  CCN_diameter_RBC=6.*CCN_mass_RBC/(M_PI*solute_density);
  CCN_diameter_RBC=pow(CCN_diameter_RBC,(1./3.));
  
  for(size=2;size<=num_CCN_size-1;size++){
    delta_diameter[size]=.5*(CCN_diameter[size+1]-CCN_diameter[size-1]); 
  }
  delta_diameter[1]=.5*(CCN_diameter[1+1]+CCN_diameter[1])-CCN_diameter_LBC;
  delta_diameter[num_CCN_size]=CCN_diameter_RBC-
    .5*(CCN_diameter[num_CCN_size]+CCN_diameter[num_CCN_size-1]);
  
  for(size=1;size<=num_CCN_size;size++){
    /* Note that this routine returns the PDF for the _length_
       distribution, not the mass distribution, as is true for the
       pristine ice crystals. */
/*    CCN_PDF=CCN_PDF_DMC94(CCN_diameter[size],CCN_scaling_diameter);*/
    CCN_PDF=(1./CCN_scaling_diameter)*exp(-CCN_diameter[size]/CCN_scaling_diameter);
    CCN_conc[size]=CCN_PDF*CCN_tot_conc*delta_diameter[size];
  } /* end loop over sizes */ 

  free(delta_diameter);

} /* end init_CCN_conc() */

float heterogeneous_nucleation_MDC92(saturation_ice,dt,tau_CCN_timescale)
     float saturation_ice;
     float dt;
     float tau_CCN_timescale;
{
  /* the evolved version of Fletcher's 1962 dragon, published by 
     Meyers et al in JAM 1992 V31 p.708.  NB: Meyers disclaims any
     similarities to reality below -20 C */ 
  float num_new_particles;
  static float a = -.639;
  static float b = .1296;
  
  num_new_particles=exp(a+b*(100.*(saturation_ice-1.))); /* #/liter */
  /* e.g., 
     S=1.01 --> nuc.= .6 per liter * dt/tau
     S=1.05 --> nuc.= 1.0 per liter * dt/tau
     S=1.2 --> nuc.= 7 per liter * dt/tau
     S=1.5 --> nuc.= 344 per liter * dt/tau */ 

  num_new_particles/=tau_CCN_timescale;
  num_new_particles*=dt;

  return num_new_particles*LITERS_PER_CUBIC_METER; /* #/m^3 */
}

float Water_Mass(concentration,crystal_mass,
		 mmr_vapor,
		 env_density,
		 num_layer,num_size,dz)
     float **concentration;
     float *crystal_mass;
     float *mmr_vapor;
     float *env_density;
     float dz;
     int num_layer,num_size;
{
  /* INPUT/UNCHANGED: concentration,crystal_mass,
     CHANGED:  */
  
  int layer,size;
  float water_mass,mass_of_air_in_column;
  
  water_mass = 0.;
  for(layer=0;layer<num_layer;layer++){
    /* add in the water mass contained in the vapor field */
    mass_of_air_in_column=env_density[layer]*dz;
    water_mass+=mass_of_air_in_column*mmr_vapor[layer];
    for(size=0;size<num_size;size++){
      water_mass+=concentration[layer][size]*crystal_mass[size]*dz; /* kg/m^2 */
    } /* end loop over sizes */
  } /* end loop over layers */

  return water_mass;
}

int Cloud_Base_idx(IWC,num_layer,time_step,total_layer_conc)
     float *IWC;
     float *total_layer_conc;
     int num_layer;
     int time_step;
{
  int cloud_base_idx;
  int layer;

  cloud_base_idx=1;
  for(layer=1;layer<=num_layer;layer++){
    if(IWC[layer] > MIN_IWC_FOR_CLOUD){
/*    if(total_layer_conc[layer] > MIN_CONC_FOR_CLOUD){*/
      cloud_base_idx=layer;
      break;
    } /* end if */ 
  } /* end loop over layers */
  if(layer == num_layer+1){
    (void)fprintf(stdout,"No cloud base at timestep %i, assigning cbi = %i\n",
		  time_step,1);
    cloud_base_idx=1;
  } /* end if */

  return cloud_base_idx;
}

int Cloud_Top_idx(IWC,num_layer,time_step,total_layer_conc)
     float *IWC;
     float *total_layer_conc;
     int num_layer;
     int time_step;
{
  int cloud_top_idx;
  int layer;

  cloud_top_idx=num_layer;
  for(layer=num_layer;layer>=1;layer--){
    if(IWC[layer] > MIN_IWC_FOR_CLOUD){
/*    if(total_layer_conc[layer] > MIN_CONC_FOR_CLOUD){*/
      cloud_top_idx=layer;
      break;
    } /* end if */ 
  } /* end loop over layers */
  if(layer == 0){
    (void)fprintf(stdout,"No cloud top at timestep %i, assigning cti = %i\n",
		  time_step,num_layer);
    cloud_top_idx=num_layer;
  } /* end if */

  return cloud_top_idx;
}

void transpose_conc(concentration,conc_transpose,num_layer,num_size)
     float **concentration,**conc_transpose;
     int num_layer,num_size;
{
  /* Fill the conc_transpose matrix with the transpose of the 
     concentration matrix */

  int layer,size;

  for(layer=0;layer<=num_layer+1;layer++){
    for(size=0;size<=num_size+1;size++){
      conc_transpose[size][layer]=concentration[layer][size];
    } /* end loop over sizes */
  } /* end loop over layers */
}

void transpose_trans(conc_transpose_new,concentration_new,num_size,num_layer)
     float **conc_transpose_new,**concentration_new;
     int num_size,num_layer;
{
  /* Fill the concentration_new matrix with the transpose of the 
     conc_transpose_new matrix, i.e. transpose back */
  int layer,size;

  for(size=0;size<=num_size+1;size++){
    for(layer=0;layer<=num_layer+1;layer++){
      concentration_new[layer][size]=conc_transpose_new[size][layer];
    } /* end loop over sizes */
  } /* end loop over layers */
}

float Dlength_Dmass(crystal_length)
     float crystal_length;
{
  /* converts all the SI units to cgs then applies the derivative of
     ramaswamy's mass/length parameterization (p. 2290 eqn. 4) to
     find dl/dm then converts everything back to SI */
  
  float dlength_dmass;
  
  dlength_dmass=28.37974367*  /* in cm/gm */
    pow(100.*crystal_length /* converts m --> cm */,
	-1.572);
  dlength_dmass*=10.;  /* converts cm/gm --> m/kg */
  return dlength_dmass;
}

float max_VEC(vector,n_elements)
     float *vector;
     int n_elements;
{
  /* Input: vector,n_elements
     Output: answer */
  
  /* returns the maximum value in the vector */
  /* Note that the argument array is assumed to be indexed [0..n_elements-1] */
  
  float big_mo_fo;
  int i;
  
  big_mo_fo=vector[0];
  for(i=1;i<n_elements;i++){
    big_mo_fo = ( vector[i] > big_mo_fo ) ? vector[i] : big_mo_fo;
  }
  return big_mo_fo;
}

void max_INTVEC_idxed(vector,n_elements,big_mo_fo,index_mo_fo)
     int *vector;
     int n_elements;
     int *big_mo_fo;
     int *index_mo_fo;
{
  /* Input: vector,n_elements
     Output: answer */
  
  /* returns the maximum value in the integer vector, and the index into
   the vector at that maximum */
  /* Note that the argument array is assumed to be indexed [0..n_elements-1] */
  
  int i;
  
  *index_mo_fo=0;
  *big_mo_fo=vector[0];
  for(i=1;i<n_elements;i++){
    if(vector[i] > *big_mo_fo){
      *big_mo_fo = vector[i];
      *index_mo_fo = i;
    } /* end if */
  } /* endfor */ 
} /* end max_INTVEC_idxed */ 

void max_FLOATVEC_idxed(vector,n_elements,big_mo_fo,index_mo_fo)
     float *vector;
     int n_elements;
     float *big_mo_fo;
     int *index_mo_fo;
{
  /* Input: vector,n_elements
     Output: answer */
  
  /* returns the maximum value in the float vector, and the index into
   the vector at that maximum */
  /* Note that the argument array is assumed to be indexed [0..n_elements-1] */
  
  int i;
  
  *index_mo_fo=0;
  *big_mo_fo=vector[0];
  for(i=1;i<n_elements;i++){
    if(vector[i] > *big_mo_fo){
      *big_mo_fo = vector[i];
      *index_mo_fo = i;
    } /* end if */
  } /* endfor */ 
} /* end max_FLOATVEC_idxed */ 

float max_abs_VEC(vector,n_elements)
     float *vector;
     int n_elements;
{
  /* Input: vector,n_elements
     Output: answer */
  
  /* returns the value whose absolute magnitude is the greatest in the vector */
  /* Note that the argument array is assumed to be indexed [0..n_elements-1] */
  
  float big_mo_fo;
  int i;
  
  big_mo_fo=vector[0];
  for(i=1;i<n_elements;i++){
    big_mo_fo = ( fabs(vector[i]) > fabs(big_mo_fo) ) ? vector[i] : big_mo_fo;
  }
  return big_mo_fo;
}

float min_VEC(vector,n_elements)
     float *vector;
     int n_elements;
{
  /* Input: vector,n_elements
     Output: answer */
  
  /* returns the minimum value in the vector */
  /* Note that the argument array is assumed to be indexed [0..n_elements-1] */
  
  float tiny_tim;
  int i;
  
  tiny_tim=vector[0];
  for(i=1;i<n_elements;i++){
    tiny_tim = ( vector[i] < tiny_tim ) ? vector[i] : tiny_tim;
  }
  return tiny_tim;
}

float total_VEC(float *vector,int n_elements)
{
  /* Input: vector,n_elements
     Output: answer */
  
  /* returns the total of all the elements of the vector */
  /* Note that the argument arrays are assumed to be indexed [0..n_elements-1] */
  
  float total;
  int i;
  
  total=vector[0];
  for(i=1;i<n_elements;i++){
    total+=vector[i];
  }
  return total;
}

float dot_product_VEC(vector1,vector2,n_elements)
     float *vector1,*vector2;
     int n_elements;
{
  /* Input: vector1,vector2,n_elements
     Output: answer */
  
  /* returns the dot product of the two vectors pointed to by vectors 1 & 2 */
  /* Note that the argument arrays are assumed to be indexed [0..n_elements-1] */
  
  float total;
  int i;
  
  total=vector1[0]*vector2[0];
  for(i=1;i<n_elements;i++){
    total+=vector1[i]*vector2[i];
  }
  return total;
}

float max_MAT(matrix,num_rows,num_cols)
     float **matrix;
     int num_rows,num_cols;
{
  /* Input: matrix,num_rows,num_cols
     Output: answer */
  
  /* returns the maximum value in the matrix */
  
  float big_mo_fo;
  int i,j;
  
  big_mo_fo=matrix[0][0];
  for(i=0;i<num_rows;i++){
    for(j=0;j<num_cols;j++){
      big_mo_fo = ( matrix[i][j] > big_mo_fo ) ? matrix[i][j] : big_mo_fo;
    }
  }
  return big_mo_fo;
}

float max_abs_MAT(matrix,num_rows,num_cols)
     float **matrix;
     int num_rows,num_cols;
{
  /* Input: matrix,num_rows,num_cols
     Output: answer */
  
  /* returns the element with largest absolute value in the matrix */
  /* Note that the argument matrix is assumed to be indexed 
     [0..num_rows-1][0..num_cols-1] */
  
  float big_mo_fo;
  int i,j;
  
  big_mo_fo=matrix[0][0];
  for(i=0;i<num_rows;i++){
    for(j=0;j<num_cols;j++){
      big_mo_fo = ( fabs(matrix[i][j]) > fabs(big_mo_fo) ) ? 
	matrix[i][j] : big_mo_fo;
    }
  }
  return big_mo_fo;
}

float max_abs_MAT_offset(matrix,num_rows,num_cols)
     float **matrix;
     int num_rows,num_cols;
{
  /* Input: matrix,num_rows,num_cols
     Output: answer */
  
  /* returns the element with largest absolute value in the matrix */
  
  /* Note that the argument matrix is assumed to be indexed 
     [1..num_rows][1..num_cols], i.e.
     assumes the matrix was dimensioned as a contiguous block of size
     (num_rows+2)*(num_cols+2) but has a symmetrical picture frame of unwanted
     data around it, row_offset rows deep on top and bottom, and
     col_offset columns wide on the left and right */

  float big_mo_fo;
  int i,j;
  
  big_mo_fo=matrix[1][1];
  for(i=1;i<=num_rows;i++){
    for(j=1;j<=num_cols;j++){
      big_mo_fo = ( fabs(matrix[i][j]) > fabs(big_mo_fo) ) ? 
	matrix[i][j] : big_mo_fo;
    }
  }
  return big_mo_fo;
}

float min_MAT(matrix,num_rows,num_cols)
     float **matrix;
     int num_rows,num_cols;
{
  /* Input: matrix,num_rows,num_cols
     Output: answer */
  
  /* returns the minimum value in the matrix */
  /* Note that the argument matrix is assumed to be indexed 
     [0..num_rows-1][0..num_cols-1] */
  
  float tiny_tim;
  int i,j;
  
  tiny_tim=matrix[0][0];
  for(i=0;i<num_rows;i++){
    for(j=0;j<num_cols;j++){
      tiny_tim = ( matrix[i][j] < tiny_tim ) ? matrix[i][j] : tiny_tim;
    }
  }
  return tiny_tim;
}

float Eqm_Vap_Ice(float temperature)
{
  /* return the equilibrium water vapor pressure over bulk solid ice, in Pa,
     for the given temperature in Kelvin
     REMEMBER: value is returned in Pa NOT mb */
  
  float vapor_pressure;
  
  vapor_pressure=24.29-6148./temperature;
  vapor_pressure=100.*exp(vapor_pressure); /* mb -> Pa */
  
  return vapor_pressure;
} /* end Eqm_Vap_Ice() */ 

float Eqm_Vap_Liquid(float temperature)
{
  /* return the equilibrium water vapor pressure over bulk liquid water, in Pa, 
     for the given temperature in Kelvin
     REMEMBER: value is returned in Pa NOT mb */
  
  float vapor_pressure;
  
  vapor_pressure=21.6-5420./temperature;
  vapor_pressure=100.*exp(vapor_pressure); /* mb -> Pa */
  
  return vapor_pressure;
} /* end Eqm_Vap_Liquid() */ 

float prism_radius_of_bullet_length(float length)
{
  /* return the radius of a bullet according to the parameterization
     of Heymsfield '72 p. 1351
     This prism radius is the fancy way to compute the curvature according 
     to Pruppacher & Klett.
     The prism radius is pictured on PrK78 p. 123. It is the radius
     of the circle inscribed in the basal plane.
     REMEMBER: this only works up to lengths of 3000 microns (3 mm) */
  
  float length_mm,width_mm,prism_radius;
  
  length_mm=length*1000.		; /* m -> mm */
  width_mm=.25*pow(length_mm,.7856)	; /* mm */
  prism_radius=.5*width_mm/1000.	; /* mm -> m, width -> radius */
  
  return prism_radius;
} /* end prism_radius_of_bullet_length() */ 

float bullet_length_to_equiv_rad(float crystal_length)
{
  /* returns the radius of the sphere the same surface area as a bullet
     shaped ice-particle via the parameterization given in Ramaswamy and
     Detwiler '86 p. 2290 and also in Zhang '89.
     Equivalent radius is used not only in Mie calculations, but for
     the radius in denom_2

     REMEMBER: input and output in meters */
  
  float length_cm,surface_area,equiv_radius;
  
  length_cm=crystal_length*100.;
  surface_area=.0334*pow(length_cm,1.572)+.505*pow(length_cm,1.786);
  equiv_radius=sqrt(.25*surface_area/M_PI); /* cm */
  equiv_radius*=METERS_PER_CENTIMETER; /* cm --> m */
  
  return equiv_radius;
} /* end bullet_length_to_equiv_rad() */ 

float equiv_rad_to_bullet_length(float equiv_rad)
{
  /* returns the length of the ice crystal bullet corresponding to
     the given effective (equivalent surface area sphere) radius. 
     Newton-Raphson must be used because the equation is screwy. */

  float crystal_length;
  float dAdl;
  float tolerance;
  float fofl;
  float length_cm;
  float length_cm_new;
  float sphere_area_cm;
  float surface_area_cm;

  int count;

  sphere_area_cm=4.*M_PI*equiv_rad*equiv_rad*
    CENTIMETERS_PER_METER*CENTIMETERS_PER_METER;

  length_cm=2.*equiv_rad*CENTIMETERS_PER_METER;

  if(debug == 16){
    (void)fprintf(stdout,"input equiv_rad=%f\n",equiv_rad);
    (void)fprintf(stdout,"sphere_area_cm2=%f\n",sphere_area_cm);
  } /* end debug */

  count=1;
  tolerance=1.;
  while(fabs(tolerance) > EMINUSFOUR && count < 25){
    count++;
    surface_area_cm=.0334*pow(length_cm,1.572)+.505*pow(length_cm,1.786);
    tolerance=(surface_area_cm-sphere_area_cm)/sphere_area_cm;
    dAdl=.0525*pow(length_cm,.572)+.90193*pow(length_cm,.786);
    fofl=surface_area_cm-sphere_area_cm;
    length_cm-=(fofl/dAdl);
    if(debug == 16){
      (void)fprintf(stdout,"count = %i length_cm = %f epsilon = %f\n",
		    count,length_cm,tolerance);
    } /* end debug */
  } /* end while */ 

  crystal_length=length_cm*METERS_PER_CENTIMETER;
  return crystal_length;
} /* end equiv_rad_to_bullet_length() */ 

float equiv_rad_to_hex_column_length(float equiv_rad)
{
  /* returns the length of the ice crystal bullet corresponding to
     the given effective (equivalent surface area sphere) radius. 
     Newton-Raphson must be used because the equation is screwy. */

  float Aspect_ratio
    (float /* hex_crystal_length */,
     float */* daspect_ratiodL */); 

  float aspect_ratio;
  float aspect_ratio_squared;
  float crystal_length;
  float daspect_ratiodL;
  float dfdL;
  float fofL;
  float guess_length;
  float hex_diameter;
  float hex_length;
  float sphere_area;
  float surface_area;
  float tolerance;

  int count;

  sphere_area=4.*M_PI*equiv_rad*equiv_rad;

  guess_length=2.4*equiv_rad;

  if(debug == 80){
    (void)fprintf(stdout,"input equiv_rad = %g microns\n",
		  equiv_rad*MICRONS_PER_METER);
    (void)fprintf(stdout,"sphere_area = %g square microns\n",
		  sphere_area*MICRONS_PER_METER*MICRONS_PER_METER);
    (void)fprintf(stdout,"guess length = %g microns\n",
		  guess_length*MICRONS_PER_METER);
  } /* end debug */

  count=1;
  tolerance=1.;
  hex_length=guess_length;
  while(fabs(tolerance) > EMINUSFOUR && count < 25){
    count++;
    aspect_ratio=Aspect_ratio(hex_length,&daspect_ratiodL);
    hex_diameter=hex_length/aspect_ratio;
    surface_area=3.*sqrt(3.)*hex_diameter*hex_diameter/4.; /* basal faces */ 
    surface_area+=3.*hex_length*hex_diameter; /* sidewalls */ 
    aspect_ratio_squared=aspect_ratio*aspect_ratio;
    fofL=surface_area-sphere_area;
    tolerance=fofL/sphere_area;
    /* see NCARIIICZP#57 */ 
    dfdL=(sqrt(3.)/4.)*(aspect_ratio-daspect_ratiodL*hex_length)+
      aspect_ratio_squared-aspect_ratio*daspect_ratiodL/2.;
    dfdL*=6.*hex_length/(aspect_ratio*aspect_ratio_squared);
    hex_length-=(fofL/dfdL);
    if(debug == 80){
      (void)fprintf(stdout,"count = %i length  = %g microns, epsilon = %g\n",
		    count,hex_length*MICRONS_PER_METER,tolerance);
    } /* end debug */
  } /* end while */ 

  return hex_length;
} /* end equiv_rad_to_hex_column_length() */ 

float Fall_Speed_Simple(crystal_length,env_pressure)
     float crystal_length,env_pressure;
{
  /* returns the fall speed of the bullet shaped ice crystal
     via the parameterization given in Heymsfield '72, the highly
     parameterized form on p. 1356. Maybe this is the form used
     in both Ramaswamy's and Zhang's models. As Heymsfield notes
     this parameterization overestimates the fall speed by ~ 25%
     
     REMEMBER: input in meters, output in meters/s */
  
  float p_1000=1000.; /* mb */
  float p_env_mb;
  float length_mm,U_1000,U;
  
  length_mm=crystal_length*1000.;  /* m --> mm */
  p_env_mb=env_pressure/100.; /* Pa --> hPa = mb */
  
  U_1000=-1.2+169.65*length_mm-84.5*length_mm*length_mm; /* cm/s */
  U=U_1000*sqrt(p_1000/p_env_mb);
  U/=100.			; /* cm/s --> m/s */
  
  return -U; /* make sure the speed is negative, to reflect falling */
}

float Fall_Speed_Complex
(float crystal_length,
 float env_pressure,
 float env_temperature,
 float env_density,
 int ice_crystal_habit)
{
  /* returns the fall speed of the bullet shaped ice crystal
     via the parameterization given in Heymsfield '72, the highly
     complex form on p. 1353. Ignores "delta" term contribution
     until I can get it working right.
     
     REMEMBER: input in meters, output in meters/s */
  
  float B;
  float C;
  float F;
  float G;
  float H;
  float U;
/*  float W_over_L;*/
  float X;
  float delta;
  float length_mm;
  float mu;
  float p_env_mb;
  float rho_f;
  float rho_s;
  float t_celsius;

  G=-1.10114;
  B=1.05687;
  C=-.09244;
  H=.00535;

  length_mm=crystal_length*MILLIMETERS_PER_METER; /* m --> mm */ 
  p_env_mb=env_pressure*MB_PER_PASCAL;     /* Pa --> hPa = mb */ 
  t_celsius=env_temperature-water_freezing_point; /* K --> C */ 
  rho_s=density_of_ice*MKS_DENSITY_TO_CGS_DENSITY; /* gm/cm^3 */ 
  /*rho_s=.78*pow(length_mm,-.0038);*/   /* gm/cm^3 */ 
  rho_f=env_density*MKS_DENSITY_TO_CGS_DENSITY;    /* kg/m^3 --> g/cm^3 */ 

  /* this is supposedly the dynamic viscosity */ 
  mu=4.301e-2*pow(t_celsius+273.,2.5)/(p_env_mb*(t_celsius+273.+120.));

  if(ice_crystal_habit == HEXAGONAL_COLUMNS){
    if(length_mm <= .2){
      F=20.*mu/length_mm;
      X=1.59e-1*rho_s*pow(length_mm,3.)/(mu*mu*rho_f);
      /*    W_over_L=.5;*/
    }else{
      F=mu/(1.973e-2*pow(length_mm,.414));
      X=9.78e-3*rho_s*pow(length_mm,1.242)/(mu*mu*rho_f);
      /*    W_over_L=.1973*pow(length_mm,-.586);*/
    } /* end else */ 
  }else if(ice_crystal_habit == HEXAGONAL_BULLETS){
    if(length_mm <= .3){
      F=mu/(2.5e-2*pow(length_mm,.7856));
      X=2.34e-2*rho_s*pow(length_mm,2.304)/(mu*mu*rho_f);
      /*    W_over_L=.25*pow(length_mm,-.214);*/
    }else{
      F=mu/(1.85e-2*pow(length_mm,.532));
      X=8.77e-3*rho_s*pow(length_mm,1.475)/(mu*mu*rho_f);
      /*    W_over_L=.185*pow(length_mm,-.486);*/
    } /* end else */ 
  }else if(ice_crystal_habit == SPHERES){
    /* Change this to the real thing when you get a chance... */ 
    if(length_mm <= .3){
      F=mu/(2.5e-2*pow(length_mm,.7856));
      X=2.34e-2*rho_s*pow(length_mm,2.304)/(mu*mu*rho_f);
      /*    W_over_L=.25*pow(length_mm,-.214);*/
    }else{
      F=mu/(1.85e-2*pow(length_mm,.532));
      X=8.77e-3*rho_s*pow(length_mm,1.475)/(mu*mu*rho_f);
      /*    W_over_L=.185*pow(length_mm,-.486);*/
    } /* end else */ 
  } /* end else */ 
  
/*  delta=exp10(-.090186+1.0034*log10(X)-.10142*pow(log10(X),2.) + */
/*	      .0083*pow(log10(X),3.) - */
/*	      exp10(-1.10114+1.05687*log10(X) - .09244*pow(log10(X),2.) + */
/*		    .00535*pow(log10(X),3.))* */
/*	      (-2.56*sqrt(W_over_L)+1.81)*/

  U=log10(X);
  U=G+B*U+C*U*U+H*U*U*U;
  U=F*exp(U*log(10.));

  U/=CENTIMETERS_PER_METER;                      /* cm/s --> m/s */
  return -U; /* make sure the speed is negative, to reflect falling */
}

float Capacitance(float length,float diameter)
{
  /* returns the capacitance of the (presumably) bullet shaped ice crystal
     in meters.  Uses prolate spheroid approximation of Pruppacher & Klett, 
     p. 449 where a=length/2 and b=radius of sphere of equal surface area
     REMEMBER: this function must be compiled after length_to_equiv_rad()
     function changed to returned capacitance values recommended by 
     Heymsfield '75 p. 821 */
  
  /* Note that in most of my work "length" refers to the c-axis while
     "diameter" refers to the a-axis. */ 
  float aspect_ratio;
  float asphericity;
  float capacitance;
  float eccentricity;
  float semi_major_axis;
  float semi_minor_axis;

  aspect_ratio=length/diameter;
  if(aspect_ratio > 1.1){
    /* prolate spheroids */ 
    semi_major_axis=length/2.;
    semi_minor_axis=diameter/2.;
    asphericity=sqrt
      (semi_major_axis*semi_major_axis-semi_minor_axis*semi_minor_axis);
    capacitance=asphericity/log((semi_major_axis+asphericity)/semi_minor_axis);
  }else if(aspect_ratio < .9){
    /* oblate spheroids */ 
    semi_major_axis=diameter/2.;
    semi_minor_axis=length/2.;
    eccentricity=sqrt(1.-semi_minor_axis*semi_minor_axis/
		      (semi_major_axis*semi_major_axis));
    capacitance=semi_major_axis*eccentricity/asin(eccentricity);
  }else{
    /* spheres */ 
    capacitance=diameter/2.;
  } /* end else */

  /* see czp #45 and czp#98 */
  /*  capacitance=length/2.;*/
  return capacitance;
}

float mass_of_bullet_length(float length) 
{
  /* returns the mass of the (presumably) bullet shaped ice crystal
     in kg.  Uses parameterization by Ramaswamy & Detwiler '86 p. 2290
     REMEMBER: input length in meters */
  
  float length_cm,mass;
  
  length_cm=length*100.;
  mass=.0137*pow(length_cm,2.572); /* g */
  return mass/1000.; /* g --> kg */
} /* end mass_of_bullet_length() */ 

float bullet_length_of_mass(float mass)
{
  /* returns the length of the (presumably) bullet shaped ice crystal
     in m.  Uses parameterization by Ramaswamy & Detwiler '86 p. 2290
     REMEMBER: input mass in kg */
  
  float mass_g,length_cm;
  
  mass_g=mass*1000.;
  length_cm=pow(mass_g/.0137,1./2.572); /* cm */
  return length_cm/100.; /* m */
} /* end bullet_length_of_mass() */ 

float Aspect_ratio
  (float hex_crystal_length,
   float *daspect_ratiodL)
{
  /* Implements the AuV70 hexagonal column aspect ratio. The pzn.
     coefficients were derived in IDL by a  
     fit to Liou's implementation of Auer and Veal's data.
     Also returns the rate of
     change of aspect ratio with respect to length, because that's
     needed by the newton-raphson routine. */ 

  float aspect_ratio_coeffs[]={5.82175,-4.99677};
  
  float aspect_ratio;

/*  float aspect_ratio_coeffs[]={0.998160,5843.33,-861429.};*/
/*    aspect_ratio= */
/*      aspect_ratio_coeffs[0]+aspect_ratio_coeffs[1]*hex_crystal_length+*/
/*      aspect_ratio_coeffs[2]*hex_crystal_length*hex_crystal_length;*/
/*    daspect_ratiodL=aspect_ratio_coeffs[1]+2.*aspect_ratio_coeffs[2]*hex_length;*/

  aspect_ratio=aspect_ratio_coeffs[0]+aspect_ratio_coeffs[1]*
    exp(-hex_crystal_length/.00050);

  *daspect_ratiodL=aspect_ratio_coeffs[1]*(-hex_crystal_length/.00050)*
    exp(-hex_crystal_length/.00050);

  return aspect_ratio;
} /* end Aspect_ratio() */ 

void hex_column_dimensions_of_mass
  (float hex_crystal_mass,
   float *length,
   float *diameter,
   float *equiv_rad,
   float *equiv_rad_squared)
{
  /* returns the length of the hex column habit ice crystal in m.  
     Uses geometry found e.g., czp#112 where the aspect ratio is 
     determined by a separate function.
     REMEMBER: input mass in kg */

  float Aspect_ratio
    (float /* hex_crystal_length */,
     float */* daspect_ratiodL */); 

  float approx_volume;
  float aspect_ratio;
  float aspect_ratio_squared;
  float basal_plane_area;
  float basal_plane_radius;
  float daspect_ratiodL;
  float dfdL;
  float float_foo2;
  float float_foo3;
  float float_foo;
  float fofL;
  float guess_length;
  float hex_diameter;
  float hex_equiv_rad_squared;
  float hex_equiv_rad;
  float hex_length;
  float hex_surface_area;
  float hex_volume;
  float tolerance;

  int count;

  hex_volume=hex_crystal_mass/density_of_ice;
  float_foo=3.*sqrt(3.)/2.;
  float_foo2=float_foo/4.;
  float_foo3=4.*3.5*hex_volume/float_foo;
  guess_length=pow(float_foo3,1./3.);
  guess_length=(guess_length < 600.e-6) ? guess_length : 2.*guess_length;

  if(debug == 77){
    (void)fprintf(stdout,"input mass = %g\n",hex_crystal_mass);
    (void)fprintf(stdout,"input volume = %g\n",hex_volume);
    (void)fprintf(stdout,"guess length = %g microns\n",guess_length*MICRONS_PER_METER);
  } /* end debug */

  /* Initialize all the arguments before entering the loop */ 
  count=1;
  tolerance=1.;
  hex_length=guess_length;
  while(fabs(tolerance) > EMINUSFOUR && count < 25){
    count++;
    aspect_ratio=Aspect_ratio(hex_length,&daspect_ratiodL);
    hex_diameter=hex_length/aspect_ratio;
    aspect_ratio_squared=aspect_ratio*aspect_ratio;
    approx_volume=float_foo2*hex_length*hex_diameter*hex_diameter;
    fofL=approx_volume-hex_volume;
    tolerance=fofL/hex_volume;
    dfdL=float_foo2*(3.*hex_length*hex_length*aspect_ratio_squared-
		     2.*hex_length*hex_length*hex_length*
		     aspect_ratio*daspect_ratiodL)/
		       (aspect_ratio_squared*aspect_ratio_squared);
    hex_length-=(fofL/dfdL);
    if(debug == 77){
      (void)fprintf(stdout,"count = %i length  = %g microns, epsilon = %g\n",
		    count,hex_length*MICRONS_PER_METER,tolerance);
    } /* end debug */
  } /* end while */ 

  aspect_ratio=Aspect_ratio(hex_length,&daspect_ratiodL);
  hex_diameter=hex_length/aspect_ratio;

  basal_plane_radius=hex_diameter/2.;
  basal_plane_area=float_foo*basal_plane_radius*basal_plane_radius;
  hex_surface_area=2.*basal_plane_area+6.*basal_plane_radius*hex_length;
  hex_equiv_rad_squared=.25*hex_surface_area/M_PI;
  hex_equiv_rad=sqrt(hex_equiv_rad_squared);

  *length=hex_length;
  *diameter=hex_diameter;
  *equiv_rad=hex_equiv_rad;
  *equiv_rad_squared=hex_equiv_rad_squared;

  if(debug == 77){
    (void)fprintf(stdout,"length = %g microns\n",hex_length*MICRONS_PER_METER);
    (void)fprintf(stdout,"diameter = %g microns\n",hex_diameter*MICRONS_PER_METER);
    (void)fprintf(stdout,"aspect ratio = %g\n\n",aspect_ratio);
  } /* end debug */
} /* end hex_column_dimensions_of_mass() */ 

float IWC_of_temp(temperature)
     float temperature;
{
  /* return the Ice Water Content (IWC) of a cirrus cloud as parameterized
     by Liou '86 and used in Zhang et al '89
     for the given temperature in Kelvin
     REMEMBER: IWC is returned in kg/m^3 */
  
  float t_celsius,foo,IWC;
  
  t_celsius=temperature-water_freezing_point;
  foo=-.2443e-3*pow(fabs(t_celsius)-20.,2.455);
  IWC=-7.6+4.*exp(foo);
  IWC=exp(IWC); /* g/m^3 */
  return IWC/1000.; /* kg/m^3 */
}

float PDF_of_length(float length,int initial_distribution_type,
		    float cloud_scale_heights_above_center)
{
  /* return the (usually initial) probability size dist. of a cirrus cloud as 
     parameterized by Dowling and Radke (JAM v.29, 1990 pp.970--978) from
     an ensemble of many observations.
     the length is input in meters, and the number of particles per cubic
     meter air per meter of length is returned.  This number should then
     be multiplied by a size domain (in meters) to obtain an absolute 
     concentration in #/m^3. */

  float float_foo;
  float length_microns;
  float log_length_microns;
  float PDF;
  float denom_factor;
  float multiplier;
  float pow_factor;

  length_microns=length*MICRONS_PER_METER; /* m -> microns */
  if(initial_distribution_type == DOWLING_RADKE){
    float_foo=exp(-length_microns/500.);
    PDF=float_foo/(3.4225*(length_microns+10.));
  }else if(initial_distribution_type == KNOLLENBERG){
    multiplier=.5*cloud_scale_heights_above_center+.5;
    log_length_microns=log(length_microns);
    multiplier*=log_length_microns-1.;
    float_foo=length_microns/10.;
    float_foo+=cloud_scale_heights_above_center/3.;
    denom_factor=2000.-100.*cloud_scale_heights_above_center;
    pow_factor=2.5+cloud_scale_heights_above_center/4.;
    PDF=20.*(2.*exp(-float_foo)+1.)*
      (1./pow(length_microns,pow_factor))*exp(-length_microns/denom_factor);
    PDF+=multiplier*PDF;
  } /* end else */
  
  /* Convert from prob. per micron of length -> prob. per meter of length */
  return PDF*MICRONS_PER_METER; 
} /* end PDF_of_length() */ 

float size_dist_of_IWC(float length,float temperature)
{
  /* return the (usually initial) size dist. of a cirrus cloud as parameterized
     by Heymsfield and Platt '85 and used in Zhang et al '89
     the length is input in meters, and the number of particles per cubic
     meter air per meter of length is returned.  This number should then
     be multiplied by a size domain (in meters) to obtain an absolute 
     concentration in #/m^3.
     for the given temperature in Kelvin
     
     REMEMBER: Neither of the distributions is fully satisfactory. */
  
  float IWC_of_temp();

  float D_sub_0;
  float IWC;
  float IWC_grams;
  float T_celsius;
  float a1;
  float a2;
  float length_microns;
  float size_dist;

  int index;

  /* Note: these contain some corrections supplied by Heymsfield:
   specifically e.g. column H, row 8, column F row 8.
   The indexing is chosen here to match his curves p. 850 */   
  static float b1[9]={0.,-2.56,-2.51,-2.21,-2.29,-3.23,-3.15,-3.83,-2.85};
  static float b2[9]={0.,-3.74,-4.49,-3.94,-4.37,-3.23,-3.15,-3.83,-2.85};

  /* #/g/micron */
  static float Nsub100byIWC[9]={0.,5170,7000,7430,19800,7500,5600,3890,5580}; 
  static float Nsub1000byIWC[9]={0.,12.0,10.4,13.7,10.3,4.86,4.00,.86,8.06}; 
  
  length_microns=length*MICRONS_PER_METER; /* m -> microns */
  T_celsius=temperature-water_freezing_point;
  
  if(T_celsius > -20){
    (void)fprintf(stdout,"temperature too high in function size_dist_of_IWC\n");
    return 0.;
  } /* endif */ 
  
  D_sub_0 = 500.;
  if (T_celsius > -25){
    index = 1; 
    D_sub_0 = 800.;
  } else if (T_celsius > -30){
    index = 2;
    D_sub_0 = 800.;
  } else if (T_celsius > -35){
    index = 3;
    D_sub_0 = 600.;
  } else if (T_celsius > -40){
    index = 4;
    D_sub_0 = 350.;
  } else if (T_celsius > -45){
    index = 5;
  } else if (T_celsius > -50){
    index = 6;
  } else if (T_celsius > -55){
    index = 7;
  } else if (T_celsius > -60){
    index = 8;
  } else if (T_celsius < -60){
    index = 8;
  } /* end else */ 
  
  IWC=IWC_of_temp(temperature); /* kg/m^3 */
  IWC_grams=IWC*GRAMS_PER_KILOGRAM; /* g/m^3 */
  
  /* The following distribution is what Heymsfield spells out, but 
     it makes the distribution independent of IWC, which is unintuitive. */
  if(length_microns <= D_sub_0){
    a1=Nsub100byIWC[index]/pow(100.,b1[index]);
    size_dist=a1*IWC_grams*pow(length_microns,b1[index]); /* #/m^3/micron */ 
  }else{
    a2=Nsub1000byIWC[index]/pow(100.,b2[index]);
    size_dist=a2*IWC_grams*pow(length_microns,b2[index]); /* #/m^3/micron */ 
  } /* end else */
  
  /* The following distribution is a combination between Liou's 
     IWC parameterization and Heymsfield's distribution which retains
     some temperature dependence. */
  
/*  a1=Nsub100byIWC[index]/pow(100.,b1[index]);*/
/*  size_dist=a1*IWC_grams*pow(length_microns,b1[index]); *//* #/m^3/micron */
  
  /* Is it a1 TIMES IWC (Heymsfield), or a1 OF IWC (Zhang)? */
  
  return size_dist*MICRONS_PER_METER; /* microns -> meters */
} /* end size_dist_of_IWC() */ 

float ice_density_of_length(/* length */)
/*     float length;*/
{
  /* return the ice_density of a cirrus cloud crystal as parameterized
     by its length by Heymsfield '72 p. 1352-53 and used in Zhang et al '89
     for the length given in meters
     REMEMBER: ice_density is returned in kg/m^3 */
  
  /* ice_density=.81*length^(-.054)	;g/m^3
     return,ice_density*1.e6/1000.	;kg/m^3 */
  
  /* Hardwiring in .8 g/cm^3 seems most reasonable, see czp 44. */
  
  return density_of_ice;
}

float Thermal_Conductivity(temperature)
     float temperature;
{
  /* computes the thermal conductivity of moist air k assuming 
     k ~ k_sub_a as per Pruppacher and Klett p. 418.
     units of thermal conductivity are J/m/s/K
     REMEMBER to input temperature in Kelvin! */
  
  float thermal_conductivity,t_celsius;
  
  t_celsius=temperature-water_freezing_point;
  thermal_conductivity=(5.69+.017*t_celsius)*EMINUSFIVE; /* cal/cm/s/C */

  /* convert from cal/cm --> J/m */
  thermal_conductivity*=JOULES_PER_CALORIE/METERS_PER_CENTIMETER; 

  return thermal_conductivity;
}

float Vapor_Diffusivity(temperature,env_pressure)
     float temperature,env_pressure;
{
  /* computes the diffusivity of water vapor in air
     as per Hall and Pruppacher '76 czr54 p.1997 and
     Pruppacher and Klett p. 413-418
     units of vapor diffusivity are m^2/s
     REMEMBER to input temperature in Kelvin! */

  float env_p_mb,vapor_diffusivity;
  
  env_p_mb=env_pressure*MB_PER_PASCAL;
  vapor_diffusivity=.211*pow(temperature/water_freezing_point,1.94)*(1013.25/env_p_mb); 
  /* cm^2/s */

  return vapor_diffusivity*SQUARE_METERS_PER_SQUARE_CM; /* cm^2/s --> m^2/s */
}

float conc_left_redistribution(concentration,mass,mass_left,mass_right)
     float concentration,mass,mass_left,mass_right;
{
  /* returns n prime sub i according to the Lagrangian mass, concentration
   conserving method of Norville, PhD Thesis, p. 18 */
  
  float conc_left;

  conc_left=concentration*(mass-mass_right)/(mass_left-mass_right);
  return conc_left;
}

float conc_right_redistribution(concentration,mass,mass_left,mass_right)
     float concentration,mass,mass_left,mass_right;
{
  /* returns n prime sub i according to the Lagrangian mass, concentration
   conserving method of Norville, PhD Thesis, p. 18 */
  
  float conc_right;

  conc_right=concentration*(mass-mass_left)/(mass_right-mass_left);
  return conc_right;
}

int find_nearest_mass_bin(lagrange_mass,crystal_mass,num_size)
     float *crystal_mass;
     float lagrange_mass;
     int num_size;
{
  /* returns the index i into the crystal_mass array which points
     to the greatest mass less than or equal the to lagrange_mass.
     i.e. crystal_mass[i] = GLB(lagrange_mass).
     NB: returns 0 when particles have shrunk smaller than the
     smallest mass bin */
  
  int i;

  i=1;
  while((crystal_mass[i] <= lagrange_mass) && (i <= num_size)){
    i++;
  }
  if(i > num_size){
    /* i.e. there is no crystal_mass[i] greater than the Lagrange mass */
    i = SIZE_TOO_BIG+1;
  }else if(i == 1){
    /* i.e. there is no crystal_mass[i] smaller than the Lagrange mass */
    i = SIZE_TOO_SMALL+1;
  }
  return --i;
}

/* Module:	cmdparse.c (Command Parse)
 * Purpose:	Set options from command line
 * Subroutine:	parse_cmdline()			returns: int
 * Subroutine:	string_cmdline()		returns: void
 * Subroutine:	usage()				returns: int
 * Xlib calls:	none
 * Copyright:	1989 Smithsonian Astrophysical Observatory
 *		You may do anything you like with this file except remove
 *		this copyright.  The Smithsonian Astrophysical Observatory
 *		makes no representations about the suitability of this
 *		software for any purpose.  It is provided "as is" without
 *		express or implied warranty.
 * Modified:	{0} Michael VanHilst	initial version	       5 January 1989
 *              {1} MVH BSDonly strings.h compatability           19 Feb 1990
 *		{n} <who> -- <does what> -- <when>
 */

void string_cmdline(int argc,char *argv[],char *cmdline,int linemax)
{
  int i;

  if( argc <= 0 ) {
    cmdline[0] = '\0';
  } else {
    (void)strcpy(cmdline, argv[0]);
    for( i=1; i<argc; i++ ) {
      (void)strncat(cmdline, " ", linemax);
      (void)strncat(cmdline, argv[i], linemax);
    }
  }
}

