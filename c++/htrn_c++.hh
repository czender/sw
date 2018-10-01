// $Id$

// Purpose: HITRAN definitions used by C++ programs

/* Copyright (C) 1997--2018 Charlie Zender
   This software is distributed under the terms of the GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

/* Usage: /Users/zender/sw/c++/htrn_c++.hh automatically generated by htrn on Mon Oct  1 10:42:00 2018
   Command: htrn */

// #include <htrn_c++.hh> // HITRAN line database definitions

#ifndef HTRN_CCC_HH // Contents have not yet been inserted in current source file
#define HTRN_CCC_HH

// C++ headers
#include <string> // Standard C++ string class

// Standard C headers

// Personal headers

// Forward declarations

// Namespaces

// Typedefs

namespace htrn{ // [nms] HITRAN namespace

  // HITRAN dimensions
  const int mlc_nbr_max_htrn(49); // [nbr] Number of gases in HITRAN database
  const int iso_nbr_max_htrn(125); // [nbr] Number of isotopomers in HITRAN database
  const int iso_per_mlc_nbr_max_htrn(11); // [nbr] Maximum number of isotopes of a molecule in HITRAN database

  // Defining mean molecular weight in SI units instead of g mol-1 keeps Avagadro's number the same and still gets rid of all those pesky factors of 1000.0
  const double mmw_H2O(1.8015261e-02); // [kg mol-1] mmw of H2O
  const double mmw_CO2(4.4009743e-02); // [kg mol-1] mmw of CO2
  const double mmw_O3(4.7997832e-02); // [kg mol-1] mmw of O3
  const double mmw_N2O(4.4012674e-02); // [kg mol-1] mmw of N2O
  const double mmw_CO(2.8010445e-02); // [kg mol-1] mmw of CO
  const double mmw_CH4(1.6043074e-02); // [kg mol-1] mmw of CH4
  const double mmw_O2(3.1998575e-02); // [kg mol-1] mmw of O2
  const double mmw_NO(3.0005630e-02); // [kg mol-1] mmw of NO
  const double mmw_SO2(6.4046674e-02); // [kg mol-1] mmw of SO2
  const double mmw_NO2(4.5992904e-02); // [kg mol-1] mmw of NO2
  const double mmw_NH3(1.7030201e-02); // [kg mol-1] mmw of NH3
  const double mmw_HNO3(6.2999296e-02); // [kg mol-1] mmw of HNO3
  const double mmw_OH(1.7006907e-02); // [kg mol-1] mmw of OH
  const double mmw_HF(2.0006386e-02); // [kg mol-1] mmw of HF
  const double mmw_HCl(3.6460710e-02); // [kg mol-1] mmw of HCl
  const double mmw_HBr(8.0911588e-02); // [kg mol-1] mmw of HBr
  const double mmw_HI(1.2791245e-01); // [kg mol-1] mmw of HI
  const double mmw_ClO(5.1447643e-02); // [kg mol-1] mmw of ClO
  const double mmw_OCS(6.0071832e-02); // [kg mol-1] mmw of OCS
  const double mmw_H2CO(3.0025657e-02); // [kg mol-1] mmw of H2CO
  const double mmw_HOCl(5.2455469e-02); // [kg mol-1] mmw of HOCl
  const double mmw_N2(2.8013603e-02); // [kg mol-1] mmw of N2
  const double mmw_HCN(2.7025618e-02); // [kg mol-1] mmw of HCN
  const double mmw_CH3Cl(5.0476203e-02); // [kg mol-1] mmw of CH3Cl
  const double mmw_H2O2(3.4005480e-02); // [kg mol-1] mmw of H2O2
  const double mmw_C2H2(2.6037999e-02); // [kg mol-1] mmw of C2H2
  const double mmw_C2H6(3.0069000e-02); // [kg mol-1] mmw of C2H6
  const double mmw_PH3(3.3997238e-02); // [kg mol-1] mmw of PH3
  const double mmw_COF2(6.6002869e-02); // [kg mol-1] mmw of COF2
  const double mmw_SF6(1.4596249e-01); // [kg mol-1] mmw of SF6
  const double mmw_H2S(3.4079355e-02); // [kg mol-1] mmw of H2S
  const double mmw_HCOOH(4.6005480e-02); // [kg mol-1] mmw of HCOOH
  const double mmw_HO2(3.2997655e-02); // [kg mol-1] mmw of HO2
  const double mmw_O(1.5994915e-02); // [kg mol-1] mmw of O
  const double mmw_ClONO2(9.7440548e-02); // [kg mol-1] mmw of ClONO2
  const double mmw_NO_ion(2.9997989e-02); // [kg mol-1] mmw of NO+
  const double mmw_HOBr(9.6906346e-02); // [kg mol-1] mmw of HOBr
  const double mmw_C2H4(2.8053350e-02); // [kg mol-1] mmw of C2H4
  const double mmw_CH3OH(3.2026215e-02); // [kg mol-1] mmw of CH3OH
  const double mmw_CH3Br(9.4927081e-02); // [kg mol-1] mmw of CH3Br
  const double mmw_CH3CN(4.1026549e-02); // [kg mol-1] mmw of CH3CN
  const double mmw_CF4(8.7993616e-02); // [kg mol-1] mmw of CF4
  const double mmw_C4H2(5.0015650e-02); // [kg mol-1] mmw of C4H2
  const double mmw_HC3N(5.1010899e-02); // [kg mol-1] mmw of HC3N
  const double mmw_H2(2.0159634e-03); // [kg mol-1] mmw of H2
  const double mmw_CS(4.4072299e-02); // [kg mol-1] mmw of CS
  const double mmw_SO3(7.9956820e-02); // [kg mol-1] mmw of SO3
  const double mmw_C2N2(5.2006148e-02); // [kg mol-1] mmw of C2N2
  const double mmw_COCl2(9.8711621e-02); // [kg mol-1] mmw of COCl2
  const double mmw_1H2_16O(1.8010565e-02); // [kg mol-1] mmw of H2O isotope 1, 1H2_16O
  const double mmw_1H2_18O(2.0014811e-02); // [kg mol-1] mmw of H2O isotope 2, 1H2_18O
  const double mmw_1H2_17O(1.9014780e-02); // [kg mol-1] mmw of H2O isotope 3, 1H2_17O
  const double mmw_1H_2H_16O(1.9016740e-02); // [kg mol-1] mmw of H2O isotope 4, 1H_2H_16O
  const double mmw_1H_2H_18O(2.1020985e-02); // [kg mol-1] mmw of H2O isotope 5, 1H_2H_18O
  const double mmw_1H_2H_17O(2.0020956e-02); // [kg mol-1] mmw of H2O isotope 6, 1H_2H_17O
  const double mmw_2H_2H_16O(2.0022915e-02); // [kg mol-1] mmw of H2O isotope 7, 2H_2H_16O
  const double mmw_12C_16O2(4.3989830e-02); // [kg mol-1] mmw of CO2 isotope 1, 12C_16O2
  const double mmw_13C_16O2(4.4993185e-02); // [kg mol-1] mmw of CO2 isotope 2, 13C_16O2
  const double mmw_16O_12C_18O(4.5994076e-02); // [kg mol-1] mmw of CO2 isotope 3, 16O_12C_18O
  const double mmw_16O_12C_17O(4.4994045e-02); // [kg mol-1] mmw of CO2 isotope 4, 16O_12C_17O
  const double mmw_16O_13C_18O(4.6997431e-02); // [kg mol-1] mmw of CO2 isotope 5, 16O_13C_18O
  const double mmw_16O_13C_17O(4.5997400e-02); // [kg mol-1] mmw of CO2 isotope 6, 16O_13C_17O
  const double mmw_12C_18O2(4.7998322e-02); // [kg mol-1] mmw of CO2 isotope 7, 12C_18O2
  const double mmw_17O_12C_18O(4.6998291e-02); // [kg mol-1] mmw of CO2 isotope 8, 17O_12C_18O
  const double mmw_12C_17O2(4.5998262e-02); // [kg mol-1] mmw of CO2 isotope 9, 12C_17O2
  const double mmw_13C_18O2(4.9001675e-02); // [kg mol-1] mmw of CO2 isotope 10, 13C_18O2
  const double mmw_18O_13C_17O(4.8001646e-02); // [kg mol-1] mmw of CO2 isotope 11, 18O_13C_17O
  const double mmw_16O3(4.7984745e-02); // [kg mol-1] mmw of O3 isotope 1, 16O3
  const double mmw_16O_16O_18O(4.9988991e-02); // [kg mol-1] mmw of O3 isotope 2, 16O_16O_18O
  const double mmw_16O_18O_16O(4.9988991e-02); // [kg mol-1] mmw of O3 isotope 3, 16O_18O_16O
  const double mmw_16O_16O_17O(4.8988960e-02); // [kg mol-1] mmw of O3 isotope 4, 16O_16O_17O
  const double mmw_16O_17O_16O(4.8988960e-02); // [kg mol-1] mmw of O3 isotope 5, 16O_17O_16O
  const double mmw_14N2_16O(4.4001062e-02); // [kg mol-1] mmw of N2O isotope 1, 14N2_16O
  const double mmw_14N_15N_16O(4.4998096e-02); // [kg mol-1] mmw of N2O isotope 2, 14N_15N_16O
  const double mmw_15N_14N_16O(4.4998096e-02); // [kg mol-1] mmw of N2O isotope 3, 15N_14N_16O
  const double mmw_14N2_18O(4.6005308e-02); // [kg mol-1] mmw of N2O isotope 4, 14N2_18O
  const double mmw_14N2_17O(4.5005278e-02); // [kg mol-1] mmw of N2O isotope 5, 14N2_17O
  const double mmw_12C_16O(2.7994915e-02); // [kg mol-1] mmw of CO isotope 1, 12C_16O
  const double mmw_13C_16O(2.8998270e-02); // [kg mol-1] mmw of CO isotope 2, 13C_16O
  const double mmw_12C_18O(2.9999161e-02); // [kg mol-1] mmw of CO isotope 3, 12C_18O
  const double mmw_12C_17O(2.8999130e-02); // [kg mol-1] mmw of CO isotope 4, 12C_17O
  const double mmw_13C_18O(3.1002516e-02); // [kg mol-1] mmw of CO isotope 5, 13C_18O
  const double mmw_13C_17O(3.0002485e-02); // [kg mol-1] mmw of CO isotope 6, 13C_17O
  const double mmw_12C_1H4(1.6031300e-02); // [kg mol-1] mmw of CH4 isotope 1, 12C_1H4
  const double mmw_13C_1H4(1.7034655e-02); // [kg mol-1] mmw of CH4 isotope 2, 13C_1H4
  const double mmw_12C_1H3_2H(1.7037475e-02); // [kg mol-1] mmw of CH4 isotope 3, 12C_1H3_2H
  const double mmw_13C_1H3_2H(1.8040830e-02); // [kg mol-1] mmw of CH4 isotope 4, 13C_1H3_2H
  const double mmw_16O2(3.1989830e-02); // [kg mol-1] mmw of O2 isotope 1, 16O2
  const double mmw_16O_18O(3.3994076e-02); // [kg mol-1] mmw of O2 isotope 2, 16O_18O
  const double mmw_16O_17O(3.2994045e-02); // [kg mol-1] mmw of O2 isotope 3, 16O_17O
  const double mmw_14N_16O(2.9997989e-02); // [kg mol-1] mmw of NO isotope 1, 14N_16O
  const double mmw_15N_16O(3.0995023e-02); // [kg mol-1] mmw of NO isotope 2, 15N_16O
  const double mmw_14N_18O(3.2002234e-02); // [kg mol-1] mmw of NO isotope 3, 14N_18O
  const double mmw_32S_16O2(6.3961901e-02); // [kg mol-1] mmw of SO2 isotope 1, 32S_16O2
  const double mmw_34S_16O2(6.5957695e-02); // [kg mol-1] mmw of SO2 isotope 2, 34S_16O2
  const double mmw_14N_16O2(4.5992904e-02); // [kg mol-1] mmw of NO2 isotope 1, 14N_16O2
  const double mmw_14N_1H3(1.7026549e-02); // [kg mol-1] mmw of NH3 isotope 1, 14N_1H3
  const double mmw_15N_1H3(1.8023583e-02); // [kg mol-1] mmw of NH3 isotope 2, 15N_1H3
  const double mmw_1H_14N_16O3(6.2995644e-02); // [kg mol-1] mmw of HNO3 isotope 1, 1H_14N_16O3
  const double mmw_1H_15N_16O3(6.3992680e-02); // [kg mol-1] mmw of HNO3 isotope 2, 1H_15N_16O3
  const double mmw_16O_1H(1.7002740e-02); // [kg mol-1] mmw of OH isotope 1, 16O_1H
  const double mmw_18O_1H(1.9006986e-02); // [kg mol-1] mmw of OH isotope 2, 18O_1H
  const double mmw_16O_2H(1.8008915e-02); // [kg mol-1] mmw of OH isotope 3, 16O_2H
  const double mmw_1H_19F(2.0006229e-02); // [kg mol-1] mmw of HF isotope 1, 1H_19F
  const double mmw_2H_19F(2.1012404e-02); // [kg mol-1] mmw of HF isotope 2, 2H_19F
  const double mmw_1H_35Cl(3.5976678e-02); // [kg mol-1] mmw of HCl isotope 1, 1H_35Cl
  const double mmw_1H_37Cl(3.7973729e-02); // [kg mol-1] mmw of HCl isotope 2, 1H_37Cl
  const double mmw_2H_35Cl(3.6982853e-02); // [kg mol-1] mmw of HCl isotope 3, 2H_35Cl
  const double mmw_2H_37Cl(3.8979904e-02); // [kg mol-1] mmw of HCl isotope 4, 2H_37Cl
  const double mmw_1H_79Br(7.9926160e-02); // [kg mol-1] mmw of HBr isotope 1, 1H_79Br
  const double mmw_1H_81Br(8.1924115e-02); // [kg mol-1] mmw of HBr isotope 2, 1H_81Br
  const double mmw_2H_79Br(8.0932336e-02); // [kg mol-1] mmw of HBr isotope 3, 2H_79Br
  const double mmw_2H_81Br(8.2930289e-02); // [kg mol-1] mmw of HBr isotope 4, 2H_81Br
  const double mmw_1H_127I(1.2791230e-01); // [kg mol-1] mmw of HI isotope 1, 1H_127I
  const double mmw_2H_127I(1.2891847e-01); // [kg mol-1] mmw of HI isotope 2, 2H_127I
  const double mmw_35Cl_16O(5.0963768e-02); // [kg mol-1] mmw of ClO isotope 1, 35Cl_16O
  const double mmw_37Cl_16O(5.2960819e-02); // [kg mol-1] mmw of ClO isotope 2, 37Cl_16O
  const double mmw_16O_12C_32S(5.9966986e-02); // [kg mol-1] mmw of OCS isotope 1, 16O_12C_32S
  const double mmw_16O_12C_34S(6.1962780e-02); // [kg mol-1] mmw of OCS isotope 2, 16O_12C_34S
  const double mmw_16O_13C_32S(6.0970341e-02); // [kg mol-1] mmw of OCS isotope 3, 16O_13C_32S
  const double mmw_16O_12C_33S(6.0966371e-02); // [kg mol-1] mmw of OCS isotope 4, 16O_12C_33S
  const double mmw_18O_12C_32S(6.1971231e-02); // [kg mol-1] mmw of OCS isotope 5, 18O_12C_32S
  const double mmw_1H2_12C_16O(3.0010565e-02); // [kg mol-1] mmw of H2CO isotope 1, 1H2_12C_16O
  const double mmw_1H2_13C_16O(3.1013920e-02); // [kg mol-1] mmw of H2CO isotope 2, 1H2_13C_16O
  const double mmw_1H2_12C_18O(3.2014811e-02); // [kg mol-1] mmw of H2CO isotope 3, 1H2_12C_18O
  const double mmw_1H_16O_35Cl(5.1971593e-02); // [kg mol-1] mmw of HOCl isotope 1, 1H_16O_35Cl
  const double mmw_1H_16O_37Cl(5.3968644e-02); // [kg mol-1] mmw of HOCl isotope 2, 1H_16O_37Cl
  const double mmw_14N2(2.8006148e-02); // [kg mol-1] mmw of N2 isotope 1, 14N2
  const double mmw_14N_15N(2.9003182e-02); // [kg mol-1] mmw of N2 isotope 2, 14N_15N
  const double mmw_1H_12C_14N(2.7010899e-02); // [kg mol-1] mmw of HCN isotope 1, 1H_12C_14N
  const double mmw_1H_13C_14N(2.8014254e-02); // [kg mol-1] mmw of HCN isotope 2, 1H_13C_14N
  const double mmw_1H_12C_15N(2.8007933e-02); // [kg mol-1] mmw of HCN isotope 3, 1H_12C_15N
  const double mmw_12C_1H_35Cl(4.9992328e-02); // [kg mol-1] mmw of CH3Cl isotope 1, 12C_1H_35Cl
  const double mmw_12C_1H_37Cl(5.1989379e-02); // [kg mol-1] mmw of CH3Cl isotope 2, 12C_1H_37Cl
  const double mmw_1H2_16O2(3.4005480e-02); // [kg mol-1] mmw of H2O2 isotope 1, 1H2_16O2
  const double mmw_12C2_1H2(2.6015650e-02); // [kg mol-1] mmw of C2H2 isotope 1, 12C2_1H2
  const double mmw_1H_12C_13C_1H(2.7019005e-02); // [kg mol-1] mmw of C2H2 isotope 2, 1H_12C_13C_1H
  const double mmw_1H_12C_12C_2H(2.7021825e-02); // [kg mol-1] mmw of C2H2 isotope 3, 1H_12C_12C_2H
  const double mmw_12C2_1H6(3.0046950e-02); // [kg mol-1] mmw of C2H6 isotope 1, 12C2_1H6
  const double mmw_12C_13C_1H6(3.1050305e-02); // [kg mol-1] mmw of C2H6 isotope 2, 12C_13C_1H6
  const double mmw_31P_1H3(3.3997238e-02); // [kg mol-1] mmw of PH3 isotope 1, 31P_1H3
  const double mmw_12C_16O_19F2(6.5991722e-02); // [kg mol-1] mmw of COF2 isotope 1, 12C_16O_19F2
  const double mmw_13C_16O_19F2(6.6995083e-02); // [kg mol-1] mmw of COF2 isotope 2, 13C_16O_19F2
  const double mmw_32S_19F6(1.4596249e-01); // [kg mol-1] mmw of SF6 isotope 1, 32S_19F6
  const double mmw_1H2_32S(3.3987721e-02); // [kg mol-1] mmw of H2S isotope 1, 1H2_32S
  const double mmw_1H_34S_1H(3.5983515e-02); // [kg mol-1] mmw of H2S isotope 2, 1H_34S_1H
  const double mmw_1H_33S_1H(3.4987105e-02); // [kg mol-1] mmw of H2S isotope 3, 1H_33S_1H
  const double mmw_1H2_12C_16O2(4.6005480e-02); // [kg mol-1] mmw of HCOOH isotope 1, 1H2_12C_16O2
  const double mmw_1H_16O2(3.2997655e-02); // [kg mol-1] mmw of HO2 isotope 1, 1H_16O2
  const double mmw_16O(1.5994915e-02); // [kg mol-1] mmw of O isotope 1, 16O
  const double mmw_35Cl_16O_14N_16O2(9.6956672e-02); // [kg mol-1] mmw of ClONO2 isotope 1, 35Cl_16O_14N_16O2
  const double mmw_37Cl_16O_14N_16O2(9.8953723e-02); // [kg mol-1] mmw of ClONO2 isotope 2, 37Cl_16O_14N_16O2
  const double mmw_14N_16O_ion(2.9997989e-02); // [kg mol-1] mmw of NO+ isotope 1, 14N_16O+
  const double mmw_1H_16O_79Br(9.5921076e-02); // [kg mol-1] mmw of HOBr isotope 1, 1H_16O_79Br
  const double mmw_1H_16O_81Br(9.7919027e-02); // [kg mol-1] mmw of HOBr isotope 2, 1H_16O_81Br
  const double mmw_12C2_1H4(2.8031300e-02); // [kg mol-1] mmw of C2H4 isotope 1, 12C2_1H4
  const double mmw_12C_13C_1H4(2.9034655e-02); // [kg mol-1] mmw of C2H4 isotope 2, 12C_13C_1H4
  const double mmw_12C2_1H3_16O_1H(3.2026215e-02); // [kg mol-1] mmw of CH3OH isotope 1, 12C2_1H3_16O_1H
  const double mmw_12C_1H3_79Br(9.3941811e-02); // [kg mol-1] mmw of CH3Br isotope 1, 12C_1H3_79Br
  const double mmw_12C_1H3_81Br(9.5939764e-02); // [kg mol-1] mmw of CH3Br isotope 2, 12C_1H3_81Br
  const double mmw_12C2_1H3_14N(4.1026549e-02); // [kg mol-1] mmw of CH3CN isotope 1, 12C2_1H3_14N
  const double mmw_12C_F4(8.7993616e-02); // [kg mol-1] mmw of CF4 isotope 1, 12C_F4
  const double mmw_12C4_1H2(5.0015650e-02); // [kg mol-1] mmw of C4H2 isotope 1, 12C4_1H2
  const double mmw_1H_12C3_14N(5.1010899e-02); // [kg mol-1] mmw of HC3N isotope 1, 1H_12C3_14N
  const double mmw_1H(2.0156500e-03); // [kg mol-1] mmw of H2 isotope 1, 1H
  const double mmw_2H(3.0218250e-03); // [kg mol-1] mmw of H2 isotope 2, 2H
  const double mmw_12C_32S(4.3971036e-02); // [kg mol-1] mmw of CS isotope 1, 12C_32S
  const double mmw_13C_32S(4.5966787e-02); // [kg mol-1] mmw of CS isotope 2, 13C_32S
  const double mmw_12C_34S(4.4974368e-02); // [kg mol-1] mmw of CS isotope 3, 12C_34S
  const double mmw_13C_34S(4.4970399e-02); // [kg mol-1] mmw of CS isotope 4, 13C_34S
  const double mmw_32S_16O3(7.9956820e-02); // [kg mol-1] mmw of SO3 isotope 1, 32S_16O3
  const double mmw_12C2_14N2(5.2006148e-02); // [kg mol-1] mmw of C2N2 isotope 1, 12C2_14N2
  const double mmw_12C_16O_35Cl2(9.7932620e-02); // [kg mol-1] mmw of COCl2 isotope 1, 12C_16O_35Cl2
  const double mmw_12C_16O_35Cl_37Cl(9.9929670e-02); // [kg mol-1] mmw of COCl2 isotope 2, 12C_16O_35Cl_37Cl

  // Integer Fortran (1-based) indices for all HITRAN molecules
  const int idx_H2O(1); // [enm] HITRAN molecule number for H2O
  const int idx_CO2(2); // [enm] HITRAN molecule number for CO2
  const int idx_O3(3); // [enm] HITRAN molecule number for O3
  const int idx_N2O(4); // [enm] HITRAN molecule number for N2O
  const int idx_CO(5); // [enm] HITRAN molecule number for CO
  const int idx_CH4(6); // [enm] HITRAN molecule number for CH4
  const int idx_O2(7); // [enm] HITRAN molecule number for O2
  const int idx_NO(8); // [enm] HITRAN molecule number for NO
  const int idx_SO2(9); // [enm] HITRAN molecule number for SO2
  const int idx_NO2(10); // [enm] HITRAN molecule number for NO2
  const int idx_NH3(11); // [enm] HITRAN molecule number for NH3
  const int idx_HNO3(12); // [enm] HITRAN molecule number for HNO3
  const int idx_OH(13); // [enm] HITRAN molecule number for OH
  const int idx_HF(14); // [enm] HITRAN molecule number for HF
  const int idx_HCl(15); // [enm] HITRAN molecule number for HCl
  const int idx_HBr(16); // [enm] HITRAN molecule number for HBr
  const int idx_HI(17); // [enm] HITRAN molecule number for HI
  const int idx_ClO(18); // [enm] HITRAN molecule number for ClO
  const int idx_OCS(19); // [enm] HITRAN molecule number for OCS
  const int idx_H2CO(20); // [enm] HITRAN molecule number for H2CO
  const int idx_HOCl(21); // [enm] HITRAN molecule number for HOCl
  const int idx_N2(22); // [enm] HITRAN molecule number for N2
  const int idx_HCN(23); // [enm] HITRAN molecule number for HCN
  const int idx_CH3Cl(24); // [enm] HITRAN molecule number for CH3Cl
  const int idx_H2O2(25); // [enm] HITRAN molecule number for H2O2
  const int idx_C2H2(26); // [enm] HITRAN molecule number for C2H2
  const int idx_C2H6(27); // [enm] HITRAN molecule number for C2H6
  const int idx_PH3(28); // [enm] HITRAN molecule number for PH3
  const int idx_COF2(29); // [enm] HITRAN molecule number for COF2
  const int idx_SF6(30); // [enm] HITRAN molecule number for SF6
  const int idx_H2S(31); // [enm] HITRAN molecule number for H2S
  const int idx_HCOOH(32); // [enm] HITRAN molecule number for HCOOH
  const int idx_HO2(33); // [enm] HITRAN molecule number for HO2
  const int idx_O(34); // [enm] HITRAN molecule number for O
  const int idx_ClONO2(35); // [enm] HITRAN molecule number for ClONO2
  const int idx_NO_ion(36); // [enm] HITRAN molecule number for NO+
  const int idx_HOBr(37); // [enm] HITRAN molecule number for HOBr
  const int idx_C2H4(38); // [enm] HITRAN molecule number for C2H4
  const int idx_CH3OH(39); // [enm] HITRAN molecule number for CH3OH
  const int idx_CH3Br(40); // [enm] HITRAN molecule number for CH3Br
  const int idx_CH3CN(41); // [enm] HITRAN molecule number for CH3CN
  const int idx_CF4(42); // [enm] HITRAN molecule number for CF4
  const int idx_C4H2(43); // [enm] HITRAN molecule number for C4H2
  const int idx_HC3N(44); // [enm] HITRAN molecule number for HC3N
  const int idx_H2(45); // [enm] HITRAN molecule number for H2
  const int idx_CS(46); // [enm] HITRAN molecule number for CS
  const int idx_SO3(47); // [enm] HITRAN molecule number for SO3
  const int idx_C2N2(48); // [enm] HITRAN molecule number for C2N2
  const int idx_COCl2(49); // [enm] HITRAN molecule number for COCl2

  // Integer Fortran (1-based) indices for all HITRAN isotopomers
  const int idx_1H2_16O(1); // [enm] HITRAN isotopomer number for H2O isotope 1, 1H2_16O
  const int idx_1H2_18O(2); // [enm] HITRAN isotopomer number for H2O isotope 2, 1H2_18O
  const int idx_1H2_17O(3); // [enm] HITRAN isotopomer number for H2O isotope 3, 1H2_17O
  const int idx_1H_2H_16O(4); // [enm] HITRAN isotopomer number for H2O isotope 4, 1H_2H_16O
  const int idx_1H_2H_18O(5); // [enm] HITRAN isotopomer number for H2O isotope 5, 1H_2H_18O
  const int idx_1H_2H_17O(6); // [enm] HITRAN isotopomer number for H2O isotope 6, 1H_2H_17O
  const int idx_2H_2H_16O(7); // [enm] HITRAN isotopomer number for H2O isotope 7, 2H_2H_16O
  const int idx_12C_16O2(8); // [enm] HITRAN isotopomer number for CO2 isotope 1, 12C_16O2
  const int idx_13C_16O2(9); // [enm] HITRAN isotopomer number for CO2 isotope 2, 13C_16O2
  const int idx_16O_12C_18O(10); // [enm] HITRAN isotopomer number for CO2 isotope 3, 16O_12C_18O
  const int idx_16O_12C_17O(11); // [enm] HITRAN isotopomer number for CO2 isotope 4, 16O_12C_17O
  const int idx_16O_13C_18O(12); // [enm] HITRAN isotopomer number for CO2 isotope 5, 16O_13C_18O
  const int idx_16O_13C_17O(13); // [enm] HITRAN isotopomer number for CO2 isotope 6, 16O_13C_17O
  const int idx_12C_18O2(14); // [enm] HITRAN isotopomer number for CO2 isotope 7, 12C_18O2
  const int idx_17O_12C_18O(15); // [enm] HITRAN isotopomer number for CO2 isotope 8, 17O_12C_18O
  const int idx_12C_17O2(16); // [enm] HITRAN isotopomer number for CO2 isotope 9, 12C_17O2
  const int idx_13C_18O2(17); // [enm] HITRAN isotopomer number for CO2 isotope 10, 13C_18O2
  const int idx_18O_13C_17O(18); // [enm] HITRAN isotopomer number for CO2 isotope 11, 18O_13C_17O
  const int idx_16O3(19); // [enm] HITRAN isotopomer number for O3 isotope 1, 16O3
  const int idx_16O_16O_18O(20); // [enm] HITRAN isotopomer number for O3 isotope 2, 16O_16O_18O
  const int idx_16O_18O_16O(21); // [enm] HITRAN isotopomer number for O3 isotope 3, 16O_18O_16O
  const int idx_16O_16O_17O(22); // [enm] HITRAN isotopomer number for O3 isotope 4, 16O_16O_17O
  const int idx_16O_17O_16O(23); // [enm] HITRAN isotopomer number for O3 isotope 5, 16O_17O_16O
  const int idx_14N2_16O(24); // [enm] HITRAN isotopomer number for N2O isotope 1, 14N2_16O
  const int idx_14N_15N_16O(25); // [enm] HITRAN isotopomer number for N2O isotope 2, 14N_15N_16O
  const int idx_15N_14N_16O(26); // [enm] HITRAN isotopomer number for N2O isotope 3, 15N_14N_16O
  const int idx_14N2_18O(27); // [enm] HITRAN isotopomer number for N2O isotope 4, 14N2_18O
  const int idx_14N2_17O(28); // [enm] HITRAN isotopomer number for N2O isotope 5, 14N2_17O
  const int idx_12C_16O(29); // [enm] HITRAN isotopomer number for CO isotope 1, 12C_16O
  const int idx_13C_16O(30); // [enm] HITRAN isotopomer number for CO isotope 2, 13C_16O
  const int idx_12C_18O(31); // [enm] HITRAN isotopomer number for CO isotope 3, 12C_18O
  const int idx_12C_17O(32); // [enm] HITRAN isotopomer number for CO isotope 4, 12C_17O
  const int idx_13C_18O(33); // [enm] HITRAN isotopomer number for CO isotope 5, 13C_18O
  const int idx_13C_17O(34); // [enm] HITRAN isotopomer number for CO isotope 6, 13C_17O
  const int idx_12C_1H4(35); // [enm] HITRAN isotopomer number for CH4 isotope 1, 12C_1H4
  const int idx_13C_1H4(36); // [enm] HITRAN isotopomer number for CH4 isotope 2, 13C_1H4
  const int idx_12C_1H3_2H(37); // [enm] HITRAN isotopomer number for CH4 isotope 3, 12C_1H3_2H
  const int idx_13C_1H3_2H(38); // [enm] HITRAN isotopomer number for CH4 isotope 4, 13C_1H3_2H
  const int idx_16O2(39); // [enm] HITRAN isotopomer number for O2 isotope 1, 16O2
  const int idx_16O_18O(40); // [enm] HITRAN isotopomer number for O2 isotope 2, 16O_18O
  const int idx_16O_17O(41); // [enm] HITRAN isotopomer number for O2 isotope 3, 16O_17O
  const int idx_14N_16O(42); // [enm] HITRAN isotopomer number for NO isotope 1, 14N_16O
  const int idx_15N_16O(43); // [enm] HITRAN isotopomer number for NO isotope 2, 15N_16O
  const int idx_14N_18O(44); // [enm] HITRAN isotopomer number for NO isotope 3, 14N_18O
  const int idx_32S_16O2(45); // [enm] HITRAN isotopomer number for SO2 isotope 1, 32S_16O2
  const int idx_34S_16O2(46); // [enm] HITRAN isotopomer number for SO2 isotope 2, 34S_16O2
  const int idx_14N_16O2(47); // [enm] HITRAN isotopomer number for NO2 isotope 1, 14N_16O2
  const int idx_14N_1H3(48); // [enm] HITRAN isotopomer number for NH3 isotope 1, 14N_1H3
  const int idx_15N_1H3(49); // [enm] HITRAN isotopomer number for NH3 isotope 2, 15N_1H3
  const int idx_1H_14N_16O3(50); // [enm] HITRAN isotopomer number for HNO3 isotope 1, 1H_14N_16O3
  const int idx_1H_15N_16O3(51); // [enm] HITRAN isotopomer number for HNO3 isotope 2, 1H_15N_16O3
  const int idx_16O_1H(52); // [enm] HITRAN isotopomer number for OH isotope 1, 16O_1H
  const int idx_18O_1H(53); // [enm] HITRAN isotopomer number for OH isotope 2, 18O_1H
  const int idx_16O_2H(54); // [enm] HITRAN isotopomer number for OH isotope 3, 16O_2H
  const int idx_1H_19F(55); // [enm] HITRAN isotopomer number for HF isotope 1, 1H_19F
  const int idx_2H_19F(56); // [enm] HITRAN isotopomer number for HF isotope 2, 2H_19F
  const int idx_1H_35Cl(57); // [enm] HITRAN isotopomer number for HCl isotope 1, 1H_35Cl
  const int idx_1H_37Cl(58); // [enm] HITRAN isotopomer number for HCl isotope 2, 1H_37Cl
  const int idx_2H_35Cl(59); // [enm] HITRAN isotopomer number for HCl isotope 3, 2H_35Cl
  const int idx_2H_37Cl(60); // [enm] HITRAN isotopomer number for HCl isotope 4, 2H_37Cl
  const int idx_1H_79Br(61); // [enm] HITRAN isotopomer number for HBr isotope 1, 1H_79Br
  const int idx_1H_81Br(62); // [enm] HITRAN isotopomer number for HBr isotope 2, 1H_81Br
  const int idx_2H_79Br(63); // [enm] HITRAN isotopomer number for HBr isotope 3, 2H_79Br
  const int idx_2H_81Br(64); // [enm] HITRAN isotopomer number for HBr isotope 4, 2H_81Br
  const int idx_1H_127I(65); // [enm] HITRAN isotopomer number for HI isotope 1, 1H_127I
  const int idx_2H_127I(66); // [enm] HITRAN isotopomer number for HI isotope 2, 2H_127I
  const int idx_35Cl_16O(67); // [enm] HITRAN isotopomer number for ClO isotope 1, 35Cl_16O
  const int idx_37Cl_16O(68); // [enm] HITRAN isotopomer number for ClO isotope 2, 37Cl_16O
  const int idx_16O_12C_32S(69); // [enm] HITRAN isotopomer number for OCS isotope 1, 16O_12C_32S
  const int idx_16O_12C_34S(70); // [enm] HITRAN isotopomer number for OCS isotope 2, 16O_12C_34S
  const int idx_16O_13C_32S(71); // [enm] HITRAN isotopomer number for OCS isotope 3, 16O_13C_32S
  const int idx_16O_12C_33S(72); // [enm] HITRAN isotopomer number for OCS isotope 4, 16O_12C_33S
  const int idx_18O_12C_32S(73); // [enm] HITRAN isotopomer number for OCS isotope 5, 18O_12C_32S
  const int idx_1H2_12C_16O(74); // [enm] HITRAN isotopomer number for H2CO isotope 1, 1H2_12C_16O
  const int idx_1H2_13C_16O(75); // [enm] HITRAN isotopomer number for H2CO isotope 2, 1H2_13C_16O
  const int idx_1H2_12C_18O(76); // [enm] HITRAN isotopomer number for H2CO isotope 3, 1H2_12C_18O
  const int idx_1H_16O_35Cl(77); // [enm] HITRAN isotopomer number for HOCl isotope 1, 1H_16O_35Cl
  const int idx_1H_16O_37Cl(78); // [enm] HITRAN isotopomer number for HOCl isotope 2, 1H_16O_37Cl
  const int idx_14N2(79); // [enm] HITRAN isotopomer number for N2 isotope 1, 14N2
  const int idx_14N_15N(80); // [enm] HITRAN isotopomer number for N2 isotope 2, 14N_15N
  const int idx_1H_12C_14N(81); // [enm] HITRAN isotopomer number for HCN isotope 1, 1H_12C_14N
  const int idx_1H_13C_14N(82); // [enm] HITRAN isotopomer number for HCN isotope 2, 1H_13C_14N
  const int idx_1H_12C_15N(83); // [enm] HITRAN isotopomer number for HCN isotope 3, 1H_12C_15N
  const int idx_12C_1H_35Cl(84); // [enm] HITRAN isotopomer number for CH3Cl isotope 1, 12C_1H_35Cl
  const int idx_12C_1H_37Cl(85); // [enm] HITRAN isotopomer number for CH3Cl isotope 2, 12C_1H_37Cl
  const int idx_1H2_16O2(86); // [enm] HITRAN isotopomer number for H2O2 isotope 1, 1H2_16O2
  const int idx_12C2_1H2(87); // [enm] HITRAN isotopomer number for C2H2 isotope 1, 12C2_1H2
  const int idx_1H_12C_13C_1H(88); // [enm] HITRAN isotopomer number for C2H2 isotope 2, 1H_12C_13C_1H
  const int idx_1H_12C_12C_2H(89); // [enm] HITRAN isotopomer number for C2H2 isotope 3, 1H_12C_12C_2H
  const int idx_12C2_1H6(90); // [enm] HITRAN isotopomer number for C2H6 isotope 1, 12C2_1H6
  const int idx_12C_13C_1H6(91); // [enm] HITRAN isotopomer number for C2H6 isotope 2, 12C_13C_1H6
  const int idx_31P_1H3(92); // [enm] HITRAN isotopomer number for PH3 isotope 1, 31P_1H3
  const int idx_12C_16O_19F2(93); // [enm] HITRAN isotopomer number for COF2 isotope 1, 12C_16O_19F2
  const int idx_13C_16O_19F2(94); // [enm] HITRAN isotopomer number for COF2 isotope 2, 13C_16O_19F2
  const int idx_32S_19F6(95); // [enm] HITRAN isotopomer number for SF6 isotope 1, 32S_19F6
  const int idx_1H2_32S(96); // [enm] HITRAN isotopomer number for H2S isotope 1, 1H2_32S
  const int idx_1H_34S_1H(97); // [enm] HITRAN isotopomer number for H2S isotope 2, 1H_34S_1H
  const int idx_1H_33S_1H(98); // [enm] HITRAN isotopomer number for H2S isotope 3, 1H_33S_1H
  const int idx_1H2_12C_16O2(99); // [enm] HITRAN isotopomer number for HCOOH isotope 1, 1H2_12C_16O2
  const int idx_1H_16O2(100); // [enm] HITRAN isotopomer number for HO2 isotope 1, 1H_16O2
  const int idx_16O(101); // [enm] HITRAN isotopomer number for O isotope 1, 16O
  const int idx_35Cl_16O_14N_16O2(102); // [enm] HITRAN isotopomer number for ClONO2 isotope 1, 35Cl_16O_14N_16O2
  const int idx_37Cl_16O_14N_16O2(103); // [enm] HITRAN isotopomer number for ClONO2 isotope 2, 37Cl_16O_14N_16O2
  const int idx_14N_16O_ion(104); // [enm] HITRAN isotopomer number for NO+ isotope 1, 14N_16O+
  const int idx_1H_16O_79Br(105); // [enm] HITRAN isotopomer number for HOBr isotope 1, 1H_16O_79Br
  const int idx_1H_16O_81Br(106); // [enm] HITRAN isotopomer number for HOBr isotope 2, 1H_16O_81Br
  const int idx_12C2_1H4(107); // [enm] HITRAN isotopomer number for C2H4 isotope 1, 12C2_1H4
  const int idx_12C_13C_1H4(108); // [enm] HITRAN isotopomer number for C2H4 isotope 2, 12C_13C_1H4
  const int idx_12C2_1H3_16O_1H(109); // [enm] HITRAN isotopomer number for CH3OH isotope 1, 12C2_1H3_16O_1H
  const int idx_12C_1H3_79Br(110); // [enm] HITRAN isotopomer number for CH3Br isotope 1, 12C_1H3_79Br
  const int idx_12C_1H3_81Br(111); // [enm] HITRAN isotopomer number for CH3Br isotope 2, 12C_1H3_81Br
  const int idx_12C2_1H3_14N(112); // [enm] HITRAN isotopomer number for CH3CN isotope 1, 12C2_1H3_14N
  const int idx_12C_F4(113); // [enm] HITRAN isotopomer number for CF4 isotope 1, 12C_F4
  const int idx_12C4_1H2(114); // [enm] HITRAN isotopomer number for C4H2 isotope 1, 12C4_1H2
  const int idx_1H_12C3_14N(115); // [enm] HITRAN isotopomer number for HC3N isotope 1, 1H_12C3_14N
  const int idx_1H(116); // [enm] HITRAN isotopomer number for H2 isotope 1, 1H
  const int idx_2H(117); // [enm] HITRAN isotopomer number for H2 isotope 2, 2H
  const int idx_12C_32S(118); // [enm] HITRAN isotopomer number for CS isotope 1, 12C_32S
  const int idx_13C_32S(119); // [enm] HITRAN isotopomer number for CS isotope 2, 13C_32S
  const int idx_12C_34S(120); // [enm] HITRAN isotopomer number for CS isotope 3, 12C_34S
  const int idx_13C_34S(121); // [enm] HITRAN isotopomer number for CS isotope 4, 13C_34S
  const int idx_32S_16O3(122); // [enm] HITRAN isotopomer number for SO3 isotope 1, 32S_16O3
  const int idx_12C2_14N2(123); // [enm] HITRAN isotopomer number for C2N2 isotope 1, 12C2_14N2
  const int idx_12C_16O_35Cl2(124); // [enm] HITRAN isotopomer number for COCl2 isotope 1, 12C_16O_35Cl2
  const int idx_12C_16O_35Cl_37Cl(125); // [enm] HITRAN isotopomer number for COCl2 isotope 2, 12C_16O_35Cl_37Cl
} // end HITRAN namespace htrn

// Prototype global functions with C++ linkages
int // O [enm] Return success code
mmw_mlc_get // [fnc] Mean molecular weight of HITRAN molecules
(double *mmw_mlc) // O [kg mol-1] Mean molecular weight of HITRAN molecules
; // end mmw_mlc_get_fnc() prototype

int // O [enm] Return success code
mmw_iso_get // [fnc] Mean molecular weight of HITRAN isotopomers
(double *mmw_iso) // O [kg mol-1] Mean molecular weight of HITRAN isotopomers
; // end mmw_iso_get_fnc() prototype

int // O [enm] Return success code
mlc_sng_get // [fnc] HITRAN molecule names
(std::string *mlc_sng) // O [sng] HITRAN molecule names
; // end mlc_sng_get_fnc() prototype

int // O [enm] Return success code
iso_sng_get // [fnc] HITRAN isotopomer names
(std::string *iso_sng) // O [sng] HITRAN isotopomer names
; // end iso_sng_get_fnc() prototype

int // O [enm] Return success code
rtl_fnc_tpt_xpn_get // [fnc] Exponent defining temperature dependence of rotational partition function
(double *xpn) // O [frc] Exponent defining temperature dependence of rotational partition function
; // end rtl_fnc_tpt_xpn_get_fnc() prototype

int // O [enm] Return success code
iso_idx_map_get // [fnc] Map [mlc_id,iso_id]->istpmr_id
(short map[htrn::mlc_nbr_max_htrn+1][htrn::iso_per_mlc_nbr_max_htrn+1]) // O [map] Map [mlc_id,iso_id]->istpmr_id
; // end iso_idx_map_fnc() prototype

// Define inline'd functions in header so source is visible to calling files

#endif // HTRN_CCC_HH
