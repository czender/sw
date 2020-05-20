// $Id$ 

// Purpose: Land Surface Model definitions

/* Copyright (C) 1997--present Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

/* Note that most of these arrays have been manipulated for easier C++ access
   Multidimensional arrays are transposed to preserve identical array storage in memory
   All arrays with a plant type (pnt_typ) dimension have been padded by one to allow indexing by the pnt_typ map */

// Usage:
// #include <lsm.hh> // Land Surface Model NB: Live code, not a header

#ifndef LSM_HH // Contents have not yet been inserted in current source file  
#define LSM_HH

namespace lsm{ // [nms] LSM namespace
  
  // Soil parameters (CCM:lsm/soiconi.F)
  const long nbr_LSM_soi_typ(5); // dry land, glacier, deep lake, shallow lake, wetland
  const prc_cmp sfc_ems[nbr_LSM_soi_typ+1]= // [frc] Surface emissivity
  {1.0e36, // Zeroth element enables array addressing by soi_typ key
   0.96,0.97,0.97,0.97,0.96};
  const prc_cmp rgh_lng[nbr_LSM_soi_typ+1]= // [m] Roughness length momentum
  {1.0e36, // Zeroth element enables array addressing by soi_typ key
   0.05,0.05,0.001,0.001,0.05};
  
  // Dimension sizes
  const long pft_nbr_CLM(16); // [nbr] Number of CLM plant functional types
  const long sgs_nbr_CLM(4); // [nbr] Number of CLM sub-gridscale patches
  const long nbr_LSM_pnt_typ(14); // [nbr] Number of LSM plant types
  const long nbr_LSM_sfc_typ(29); // [nbr] Number of LSM surface types
  const long nbr_LSM_sgs_sfc(3); // [nbr] Number of LSM sub-gridscale patches
  const long nbr_mth_per_yr(12); // [nbr] Number of months per year
  
  // Vegetation parameters
  const prc_cmp cwpvt[nbr_LSM_pnt_typ+1]= // [frc] Canopy wind parameter
  {1.0e36, // Zeroth element enables array addressing by pnt_typ key
   3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0};
  const prc_cmp z0mvt[nbr_LSM_pnt_typ+1]= // [m] Momentum roughness length
  {1.0e36, // Zeroth element enables array addressing by pnt_typ key
   0.94,0.77,2.62,1.10,0.99,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.06,0.00};
  const prc_cmp zpdvt[nbr_LSM_pnt_typ+1]= // [m] Zero plane displacement height
  {1.0e36, // Zeroth element enables array addressing by pnt_typ key
   11.39,9.38,23.45,13.40,12.06,0.34,0.34,0.34,0.34, 0.34, 0.34, 0.34,0.34,0.00};
  const prc_cmp hvt[nbr_LSM_pnt_typ+1]= // [m] Height at top of canopy
  {1.0e36, // Zeroth element enables array addressing by pnt_typ key
   17.0,14.0,35.0,20.0,18.0,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.0};
  const prc_cmp dleaf[nbr_LSM_pnt_typ+1]= // [m] Characteristic leaf dimension
  {1.0e36, // Zeroth element enables array addressing by pnt_typ key
   0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.00};
  
  // [m2 m-2] Monthly leaf area index + stem area index, one-sided
  // These LSM data are mid-month values, i.e., valid on the 15th of the month
  // C storage tai[nbr_LSM_pnt_typ][nbr_mth_per_yr] same as Fortran storage tai(nbr_LSM_pnt_typ,nbr_mth_per_yr)
  const prc_cmp tai[nbr_LSM_pnt_typ+1][nbr_mth_per_yr]={
    {1.0e36,1.0e36,1.0e36,1.0e36,1.0e36,1.0e36,1.0e36,1.0e36,1.0e36,1.0e36,1.0e36,1.0e36}, // Zeroth row enables array addressing by pnt_typ key
    {4.5,4.7,5.0,5.1,5.3,5.5,5.3,5.3,5.2,4.9,4.6,4.5}, //  1 Needleleaf evergreen tree
    {0.3,0.3,0.3,1.0,1.6,2.4,4.3,2.9,2.0,1.3,0.8,0.5}, //  2 Needleleaf deciduous tree
    {5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0}, //  3 Broadleaf evergreen tree
    {0.4,0.4,0.7,1.6,3.5,5.1,5.4,4.8,3.8,1.7,0.6,0.4}, //  4 Broadleaf deciduous tree
    {1.2,1.0,0.9,0.8,0.8,1.0,2.0,3.7,3.2,2.7,1.9,1.2}, //  5 Tropical seasonal tree
    {0.7,0.8,0.9,1.0,1.5,3.4,4.3,3.8,1.8,1.0,0.9,0.8}, //  6 Cool grass (C3)
    {1.3,1.3,1.3,1.3,1.3,1.3,1.3,1.3,1.3,1.3,1.3,1.3}, //  7 Evergreen shrub
    {1.0,1.0,0.8,0.3,0.6,0.0,0.1,0.3,0.5,0.6,0.7,0.9}, //  8 Deciduous shrub
    {0.1,0.1,0.1,0.1,0.1,0.3,1.5,1.7,1.4,0.1,0.1,0.1}, //  9 Arctic deciduous shrub
    {0.7,0.8,0.9,1.0,1.5,3.4,4.3,3.8,1.8,1.0,0.9,0.8}, // 10 Arctic grass
    {0.0,0.0,0.0,0.0,1.0,2.0,3.0,3.0,1.5,0.0,0.0,0.0}, // 11 Crop
    {0.0,0.0,0.0,0.0,1.0,2.0,3.0,3.0,1.5,0.0,0.0,0.0}, // 12 Irrigated crop
    {0.7,0.8,0.9,1.0,1.5,3.4,4.3,3.8,1.8,1.0,0.9,0.8}, // 13 Warm grass (C4)
    {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}  // 14 Bare ground
  }; // end tai[]
  
  // [m2 m-2] Monthly leaf area index, one-sided
  // These LSM data are mid-month values, i.e., valid on the 15th of the month
  // C storage gai[nbr_LSM_pnt_typ][nbr_mth_per_yr] same as Fortran storage gai(nbr_LSM_pnt_typ,nbr_mth_per_yr)
  const prc_cmp gai[nbr_LSM_pnt_typ+1][nbr_mth_per_yr]={
    {1.0e36,1.0e36,1.0e36,1.0e36,1.0e36,1.0e36,1.0e36,1.0e36,1.0e36,1.0e36,1.0e36,1.0e36}, // Zeroth row enables array addressing by pnt_typ key
    {4.1,4.2,4.6,4.8,4.9,5.0,4.8,4.7,4.6,4.2,4.0,4.0}, //  1 Needleleaf evergreen tree
    {0.0,0.0,0.0,0.6,1.2,2.0,2.6,1.7,1.0,0.5,0.2,0.0}, //  2 Needleleaf deciduous tree
    {4.5,4.5,4.5,4.5,4.5,4.5,4.5,4.5,4.5,4.5,4.5,4.5}, //  3 Broadleaf evergreen tree
    {0.0,0.0,0.3,1.2,3.0,4.7,4.5,3.4,1.2,0.3,0.0,0.0}, //  4 Broadleaf deciduous tree
    {0.8,0.7,0.4,0.5,0.5,0.7,1.7,3.0,2.5,1.6,1.0,1.0}, //  5 Tropical seasonal tree
    {0.4,0.5,0.6,0.7,1.2,3.0,3.5,1.5,0.7,0.6,0.5,0.4}, //  6 Cool grass (C3)
    {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0}, //  7 Evergreen shrub
    {0.9,0.8,0.2,0.2,0.0,0.0,0.0,0.2,0.4,0.5,0.6,0.8}, //  8 Deciduous shrub
    {0.0,0.0,0.0,0.0,0.0,0.2,1.4,1.2,0.0,0.0,0.0,0.0}, //  9 Arctic deciduous shrub
    {0.4,0.5,0.6,0.7,1.2,3.0,3.5,1.5,0.7,0.6,0.5,0.4}, // 10 Arctic grass
    {0.0,0.0,0.0,0.0,1.0,2.0,3.0,3.0,1.5,0.0,0.0,0.0}, // 11 Crop
    {0.0,0.0,0.0,0.0,1.0,2.0,3.0,3.0,1.5,0.0,0.0,0.0}, // 12 Irrigated crop
    {0.4,0.5,0.6,0.7,1.2,3.0,3.5,1.5,0.7,0.6,0.5,0.4}, // 13 Warm grass (C4)
    {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}  // 14 Bare ground
  }; // end gai[]
  
  // [idx] LSM plant type (1..14 = nbr_LSM_pnt_typ)
  // C storage pnt_typ[sgs_sfc][sfc_typ] differs from Fortran storage pnt_typ(sfc_typ,sgs_sfc)
  const long pnt_typ[nbr_LSM_sgs_sfc][nbr_LSM_sfc_typ]=
  {{ 0, // pnt_typ[0][*]
     14,  14,   1,   2,   4,   1,   1, 
     4,   1,   3,   5,  13,   1,   2,
     11,  11,   6,  13,   9,   7,   8,
     8,  12,  11,  12,  11,   3,  14},
   {0, // pnt_typ[1][*]
    14,  14,  14,  14,  14,   4,  14,
    14,   4,  14,  14,   5,  10,  10,
    4,   4,  13,   6,  10,  14,  14,
    14,  14,  14,  14,  14,  14,  14},
   {0, // pnt_typ[2][*]
    14,  14,  14,  14,  14,  14,  14,
    14,  14,  14,  14,  14,  14,  14,
    1,   1,  14,  14,  14,  14,  14,
    14,  14,  14,  14,  14,  14,  14}}; // end pnt_typ[]
  
  // [frc] Weight of corresponding plant type (sums to 1.0)
  // C storage pnt_frc[sgs_sfc][sfc_typ] differs from Fortran storage pnt_frc(sfc_typ,sgs_sfc)
  const prc_cmp pnt_frc[nbr_LSM_sgs_sfc][nbr_LSM_sfc_typ]=
  {{0.00, // pnt_frc[0][*]
    1.00,1.00,0.75,0.50,0.75,0.37,0.75,
    0.75,0.37,0.95,0.75,0.70,0.25,0.25,
    0.40,0.40,0.60,0.60,0.30,0.80,0.80,
    0.10,0.85,0.85,0.85,0.85,0.80,1.00},
   {0.00, // pnt_frc[1][*]
    0.00,0.00,0.25,0.50,0.25,0.37,0.25,
    0.25,0.37,0.05,0.25,0.30,0.25,0.25,
    0.30,0.30,0.20,0.20,0.30,0.20,0.20,
    0.90,0.15,0.15,0.15,0.15,0.20,0.00},
   {0.00, // pnt_frc[2][*]
    0.00,0.00,0.00,0.00,0.00,0.26,0.00,
    0.00,0.26,0.00,0.00,0.00,0.50,0.50,
    0.30,0.30,0.20,0.20,0.40,0.00,0.00,
    0.00,0.00,0.00,0.00,0.00,0.00,0.00}}; // end pnt_frc[]
  
} // end LSM namespace lsm

#endif // LSM_HH
