True=1
False=0
;
; fields: structure to hold predefined fields informations.
;
fld_nbr_in_database=22
foo={fld_sct, $
	index:0, $
	name:'', $
	english:'', $
	long_nm:'', $
	unit:'', $
	scale:0.0, $
	symbol:'', $
	abbrev:'', $
	bin_sz:0.0}
fld_lst=replicate({fld_sct},fld_nbr_in_database)
;
fld_lst(0)={fld_sct,0,'dom_avg_precip','dom_avg_precip: The precipitation averaged over the entire domain','!5Domain-Averaged Precipitation','!5mm day!E-1!N',8.64e7,'!5!8P!5!IT!N','!5Precip',1.}
;
fld_lst(1)={fld_sct,0,'dom_avg_CWP_partial','dom_avg_CWP_partial: The condensed water path above a given level','!5!8CWP!5','!5g m!E-2!N',1.0e3,'!5!8CWP!5','!5!8CWP!5',25.}
;
fld_lst(2)={fld_sct,0,'dom_avg_IWP_partial','dom_avg_IWP_partial: The ice water path above a given level','!5!8IWP!5','!5g m!E-2!N',1.0e3,'!5!8IWP!5','!5!8IWP!5',25.}
;
fld_lst(3)={fld_sct,0,'conv_avg_mass_flux','conv_avg_mass_flux: The sum of the convective mass fluxes at a given level averaged over the entire domain','!5!8M!I!5c!N','!5g m!E-2!N s!E-1!N',1.0e3,'!5!8M!I!5c!N','!5!8M!I!5c!N',5.}
;
fld_lst(4)={fld_sct,0,'cloud_frac','cloud_frac: The fraction of grid cells at a given level that met the cloudiness criteria','!5Cloud Fraction','!5%',1.0e2,'!8A!I!5c!N','!5Cloud Frac.',5.}
;
fld_lst(5)={fld_sct,0,'dom_avg_q_sat_ice','dom_avg_q_sat_ice: The domain average saturated vapor specific humidity with respect to an ice water surface','!5Ice Sat. Spec. Hum. !8Q!5!Ii!N','!5g/kg',1.0e3,'!8Q!5!Ii!N','!5Sat. Spec. Hum. !8Q!5!Ii!N',.1}
;
fld_lst(6)={fld_sct,0,'dom_avg_q_sat_liq','dom_avg_q_sat_liq: The domain average saturated vapor specific humidity with respect to a liquid water surface','!5g/kg','!5Liquid Saturated Specific Humidity',1.0e3,'!8Q!5!Il!N','!5Sat. Spec. Hum. !8Q!5!Il!N',.1}
;
fld_lst(7)={fld_sct,0,'dom_avg_q_vapor','dom_avg_q_vapor: ','!5Specific Humidity !8q!5','!5g/kg',1.0e3,'!8q!5','!5Spec. Humidity',.1}
;
;fld_lst(8)={fld_sct,0,'dIWPdt','dIWPdt: time derivative of partial ice water path','!9D!8!It!5!N!8IWP!5','!5g m!E-2!N s!E-1!N',1.0e3,'!9D!8!It!5!N!8IWP!5','!9D!8!It!5!N!8IWP!5',1.}
fld_lst(8)={fld_sct,0,'dIWPdt','dIWPdt: time derivative of partial ice water path','!9D!8!It!5!N!8IWP!5','!5g m!E-2!N (20 min)!E-1!N',1200*1.0e3,'!9D!8!It!5!N!8IWP!5','!9D!8!It!5!N!8IWP!5',1.}
;
fld_lst(9)={fld_sct,0,'dCWPdt','dCWPdt: time derivative of partial condensed water path','!9D!8!It!5!N!8CWP!5','!5g m!E-2!N s!E-1!N',1.0e3,'!9D!8!It!5!N!8CWP!5','!9D!8!It!5!N!8CWP!5',1.}
;
fld_lst(10)={fld_sct,0,'conv','rhs: convective source term','!5Convective Source','!5g kg!E-1!N s!E-1!N',1.0e0,'!5!8M!5!Ic!N','!5Conv. Source',1.}
;
fld_lst(11)={fld_sct,0,'precip','rhs: precipitation sink term','!5Precipitation Sink','!5g kg!E-1!N s!E-1!N',1.0e0,'!8IWP!5','!5Precipitation',1.}
;
fld_lst(12)={fld_sct,0,'evap','rhs: evaporation sink term','!5Evaporation Sink','!5g kg!E-1!N s!E-1!N',1.0e0,'!5(!8q!5 - !8q!5!Ii!N)!8P!5!E1/2!N','!5Evaporation',1.}
;
fld_lst(13)={fld_sct,0,'cond','rhs: condensation source term','!5Condensation Source','!5%',1.0e2,'!8q!5/!8q!5!Ii!N','!5Condensation',1.}
;
fld_lst(14)={fld_sct,0,'conv_precip','rhs: convective source + precipitation sink terms','!5Convective Source - Precipitation Sink','!5g m!E-2!N (20 min)!E-1!N',1200*1.0e3,'!7a!8M!5!Ic!N - !7b!8IWP!5','!5Conv. - Precip.',1.}
;
fld_lst(15)={fld_sct,0,'conv_cond_precip','rhs: convective source + condensational source - precipitation sink','!5Conv. Cond. Precip. RHS','!5g m!E-2!N (20 min)!E-1!N',1200*1.0e3,'!5RHS','!5RHS',1.}
;
fld_lst(16)={fld_sct,0,'dIWCdt','dIWCdt: time derivative of ice water content','!9D!8!It!5!N!8IWC!5','!5g m!E-3!N (20 min)!E-1!N',1200*1.0e3,'!9D!8!It!5!N!8IWC!5','!9D!8!It!5!N!8IWC!5',1.}
;
fld_lst(17)={fld_sct,0,'dom_avg_IWC','dom_avg_IWC: The local ice water content','!5!8IWC!5','!5g m!E-3!N',1.0e3,'!5!8IWC!5','!5!8IWC!5',25.}
;
fld_lst(18)={fld_sct,0,'dRHdt','dRHdt: time derivative of relative humidity w/r/t ice','!9D!8!It!5!N!8RH!5!Ii!N','!5% (20 min)!E-1!N',1200*1.0e2,'!9D!8!It!5!N!8RH!5!Ii!N','!9D!8!It!5!N!8RH!5!Ii!N',1.}
;
fld_lst(19)={fld_sct,0,'dom_avg_RH_ice','dom_avg_RH_ice: relative humidity w/r/t ice','!5!8RH!5!Ii!N','!5%',1.0e2,'!5!8RH!5!Ii!N','!5!8RH!5!Ii!N',25.}
;
fld_lst(20)={fld_sct,0,'dqdt','dqdt: time derivative of specific humidity','!9D!8!It!5!N!8q!5','!5g kg!E-1!N (20 min)!E-1!N',1200*1.0e3,'!9D!8!It!5!N!8q!5','!9D!8!It!5!N!8q!5',1.}
;
fld_lst(21)={fld_sct,0,'conv_frac','conv_frac: The fraction of grid cells at a given level that met the convective criteria','!5Convective Fraction','!5%',1.0e2,'!8F!I!5c!N','!5Conv Frac.',.1}
;
;fld_lst(0)={fld_sct,0,'foo','foo: ','!5','!5',1.0e0,'!5','!5',1.}
;
;fld_lst=fld_lst(sort(fld_lst.name))
fld_lst.index=indgen(fld_nbr_in_database)
;
