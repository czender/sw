/* $Id$
 *
 * $Log: not supported by cvs2svn $
 * Revision 1.1.1.1  1998-09-15 02:06:41  zender
 * Imported sources
 *
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
 * Revision 1.1  1993/04/22  23:44:46  zender
 * Initial revision
 *
 ***/

/* Max size of a file name */  
#define FILESIZE 80                     

/* Max size of a graph label */
#define LABEL_SIZE 100

/* Max size of the shell command */ 
#define CMDLINE_SIZE 300                

/* used for Lagrangian growth */    
#define SIZE_TOO_BIG 7373               

/* used for Lagrangian growth */    
#define SIZE_TOO_SMALL 3737             

/* used for arbitrary optical values */ 
#define BOGUS_FLOAT .737373737          

/* the maximum number of dimensions of data stored in any NetCDF variable */ 
#define MAX_NUM_NETCDF_DIMS 3

/* current # of intervals used to compute IR cooling by particles */
#define NUM_IR_SPECTRAL_BANDS 31

/* current # of intervals and pseudo-intervals in the CCM2 model */  
#define NUM_CCM2_SPECTRAL_BANDS 18      

/* number of points in the USSA atmosphere in init_temp_pres() */ 
#define NUM_MLS_STD_ATM_LEVELS 38           

/* number of points in the USSA atmosphere in init_temp_pres() */ 
#define NUM_TROPICS_STD_ATM_LEVELS 91

/* atmospheric_profile_type used in init_temp_pres() and main() */ 
enum atmospheric_profile_type{
  MID_LATITUDE_SUMMER,
  TROPICS};

/* form of initial ice crystal distribution function employed */ 
enum initial_distribution_type{
  HEYMSFIELD,
  DOWLING_RADKE,
  KNOLLENBERG,
  RAMASWAMY,
  ZHANG};

/* what shape will the crystals be for all thermodynamical and dynamical
 processes? */ 
enum ice_crystal_habit{
  HEXAGONAL_COLUMNS,
  HEXAGONAL_BULLETS,
  SPHERES};

#ifdef cray
#define FORTRAN_callbh CALLBH
#define FORTRAN_colmod COLMOD
#define FORTRAN_consre CONSRE
#define FORTRAN_etbfct ETBFCT
#define FORTRAN_ngride NGRIDE
#define FORTRAN_ogride OGRIDE
#define FORTRAN_refice REFICE
#define FORTRAN_spline SPLINE
#define FORTRAN_veloce VELOCE
#endif
#ifdef sun
#define FORTRAN_callbh callbh_
#define FORTRAN_colmod colmod_
#define FORTRAN_consre consre_
#define FORTRAN_etbfct etbfct_
#define FORTRAN_ngride ngride_
#define FORTRAN_ogride ogride_
#define FORTRAN_refice refice_
#define FORTRAN_spline spline_
#define FORTRAN_veloce veloce_
#endif
#ifdef sgi
#define FORTRAN_callbh callbh_
#define FORTRAN_colmod colmod_
#define FORTRAN_consre consre_
#define FORTRAN_etbfct etbfct_
#define FORTRAN_ngride ngride_
#define FORTRAN_ogride ogride_
#define FORTRAN_refice refice_
#define FORTRAN_spline spline_
#define FORTRAN_veloce veloce_
#endif
#ifdef ibm
#define FORTRAN_callbh callbh
#define FORTRAN_colmod colmod
#define FORTRAN_consre consre
#define FORTRAN_etbfct etbfct
#define FORTRAN_ngride ngride
#define FORTRAN_ogride ogride
#define FORTRAN_refice refice
#define FORTRAN_spline spline
#define FORTRAN_veloce veloce
#endif

#define EMINUSTHREE 1.e-3
#define EMINUSFOUR 1.e-4
#define EMINUSFIVE 1.e-5
#define EMINUSSIX 1.e-6
#define EMINUSNINE 1.e-9
#define EMINUSTEN 1.e-10
#define EMINUSTWELVE 1.e-12
#define EMINUSFIFTEEN 1.e-15
#define EMINUSTWENTY 1.e-20
#define EMINUSTHIRTY 1.e-30
#define EMINUSTHIRTYSIX 1.e-36

#define Boolean int
#define True 1
#define False 0
#define YES 1
#define NO 0

#if ( ! defined M_PI )
#define M_PI 3.14159265358979323846     /* for Sun compatability */ 
#endif

#define H2SO4_deliquescence_point 0.    /* liquid saturation ratio */ 
#define NH42SO4_deliquescence_point .82 /* liquid saturation ratio */ 
#define adiabatic_lapse .00975786       /* K/m = gravity/spec_heat_air */ 
#define ammonium_sulfate_density 1.769e3 /* kg/m^3 DMC94 p. 80 */ 
#define boltzmann_constant 1.38063e-23  /* J/K */ 
#define density_of_ice 800.             /* kg/m^3 = .8 g/cm^3 */ 
#define dry_lapse_rate .01		/* K/m = 10 K/km */
#define epsilon_vapor .622		/* molec wgt vapor/molec wgt dry air */
#define gas_const_dry_air 287.05	/* J/kg/K */
#define gas_const_vapor 461.51		/* J/kg/K */
#define gravity 9.80665                 /* m/s^2 */ 
#define hcfirst 1.986488377e-25         /* Jm conserve precision */ 
#define hcsquared 5.9553531e-17         /* Jm^2/s conserve precision */ 
#define kappa .2856                     /* gas_const_dry_air/spec_heat_air */ 
#define latent_heat_sub 2.834e6    	/* J/kg (latent heat of sublimation) */
#define moist_lapse_rate .0065		/* K/m = 6.5 K/km */
#define planck_constant 6.62620e-34     /* Js */  
#define spec_heat_air 1005.             /* J/kg/K */
#define speed_of_light 2.99793e8        /* m/s */ 
#define stefan_boltzmann_constant 5.67e-8 /* J/m^2/K^4/s */ 
#define surf_free_energy_ice .106	/* J/m^2 = N/m Pruppacher & Klett p. 121  */
#define triple_point_of_water 273.16    /* K */ 
#define water_freezing_point 273.15     /* K */  

#define CENTIMETERS_PER_METER 100.
#define CUBIC_CMS_PER_CUBIC_METER 1.e6
#define CUBIC_METERS_PER_LITER .001
#define JOULES_PER_CALORIE 4.1855       /* Cal = energy to heat 1g H20 1C @ 15C */ 
#define GRAMS_PER_KILOGRAM 1000.
#define KILOGRAMS_PER_GRAM .001
#define KILOGRAMS_PER_MILLIGRAM 1.e-6
#define KILOMETERS_PER_METER .001
#define LITERS_PER_CUBIC_METER 1000.
#define MB_PER_PASCAL .01
#define METERS_PER_CENTIMETER .01
#define METERS_PER_KILOMETER 1000.
#define METERS_PER_MICRON 1.e-6
#define MICRONS_PER_METER 1.e6
#define MICROJOULES_PER_JOULE 1.e6
#define MICROWATTS_PER_WATT 1.e6
#define MILLIGRAMS_PER_KILOGRAM 1.e6
#define MILLIJOULES_PER_JOULE 1.e3
#define MILLIMETERS_PER_METER 1.e3
#define MILLIWATTS_PER_WATT 1.e3
#define MINUTES_PER_SECOND .016666666
#define MKS_DENSITY_TO_CGS_DENSITY .001
#define NANOJOULES_PER_JOULE 1.e9
#define NANOWATTS_PER_WATT 1.e9
#define NANOGRAMS_PER_KILOGRAM 1.e12
#define ONE_MICRON 1.e-6
#define PASCALS_PER_MB 100.
#define SECONDS_PER_DAY 86400.
#define SECONDS_PER_HOUR 3600.
#define SQUARE_CMS_PER_SQUARE_METER 1.e4
#define SQUARE_METERS_PER_SQUARE_CM 1.e-4
#define SQUARE_MICRONS_PER_SQUARE_METER 1.e12

typedef struct {
  float real;
  float imag;
} complex;



