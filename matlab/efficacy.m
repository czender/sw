% Report change in temperature, temperature trend, efficacy, and
% range of uncertainty for control and experiment annual-mean,
% global-mean timeseries using the caseids 'xpt' and 'ctl' below:


clear;

% =1: Analyze experiment only
% =2: Analyze experiment and control
flg_cfg = 2;

% TOA forcing/ surface forcing ratio:
Ft_ovr_Fs = 0.91;

% CO2 climate sensitivity (delta T / adjusted forcing)
% CAM3/SOM sensitivity from Kiehl et al., 2006
sens_co2 = 0.69; 

% case id's of experiment and control:
xpt='cspdfdb07';
ctl='cspdfdb08';

fl1=strcat('/data/zender/anl_',xpt,'/',xpt,'_ts_ANN_xyt.nc');
fl2=strcat('/data/zender/anl_',ctl,'/',ctl,'_ts_ANN_xyt.nc');

% load relevant data:
nc1=netcdf(fl1,'nowrite');
txpt=nc1{'TREFHT'}(:);
frc=nc1{'SOTFRC_SFC'}(:);
frct=nc1{'SOTAFRC_TOP'}(:);
frc_lnd=nc1{'SOTFRC_SFCL'}(:);
frc_ice=nc1{'SOTFRC_SFCI'}(:);
nc1=close(nc1);


if (flg_cfg > 1)
  nc1=netcdf(fl2,'nowrite');
  tctl=nc1{'TREFHT'}(:);
  nc1=close(nc1);
end;
  

%%%%%  1. Temperature trends (95% confidence interval of slope):
design(1:length(txpt),1)=1;
design(1:length(txpt),2)=[1:length(txpt)];

[B,BINT]=regress(txpt,design);
tslop_xpt = B(2);
tslop_rng_xpt = BINT(2,:)

if (flg_cfg > 1)
  design(1:length(tctl),1)=1;
  design(1:length(tctl),2)=[1:length(tctl)];

  [B,BINT]=regress(tctl,design);
  tslop_ctl = B(2);
  tslop_rng_ctl = BINT(2,:)
end;
  

%%%%% 2. Difference in temperature and standard error:

if (flg_cfg > 1)
  n1=length(txpt);
  n2=length(tctl);
  s1=std(txpt);
  s2=std(tctl);
  
  % standard error (definition 1): (Prather):
  %se1=sqrt( s1^2 + s2^2 );

  % standard error (definition 2): unpaired pools, equal variance: 
  sp=sqrt( ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) );
  se_t=sp*sqrt((1/n1) + (1/n2));

  % standard error (definition 3): unpaired pools, unequal variance:
  %se_t=sqrt( (s1^2/n1) + (s2^2/n2));
  
  dT = mean(txpt)-mean(tctl)
  dTpm = se_t
end;
  


%%%%% 3. Forcing and standard error:
Fs = mean(frc)
Fspm = std(frc)
Fspm_rat = Fspm / Fs

Ft = mean(frct)
Ftpm = std(frct)
Ftpm_rat = Ftpm / Ft

frc_lnd = mean(frc_lnd)/(mean(frc_lnd)+mean(frc_ice))


%%%%% 4. Efficacy, and standard error, derived from temperature and
%%%%%    forcing standard errors in quadrature: 
if (flg_cfg > 1)
  err_frch = abs(1/Fs - 1/(Fs-Fspm));
  err_frcl = abs(1/Fs - 1/(Fs+Fspm));

  err_effh = sqrt( (dTpm/dT)^2 + (err_frch/(1/Fs))^2 );
  err_effl = sqrt( (dTpm/dT)^2 + (err_frcl/(1/Fs))^2 );

  effc = dT/(Fs*Ft_ovr_Fs) / sens_co2
  effh = effc + effc*err_effh;
  effl = effc - effc*err_effl;
  effph = effh-effc
  effpm = effl-effc

  % obsolete: efficacy error derived only from temperature change uncertainty
  effc2 = ((mean(txpt)-mean(tctl))/(mean(frc)*Ft_ovr_Fs)) / sens_co2;
  effh2 = ((mean(txpt)-mean(tctl)+se_t)/(mean(frc)*Ft_ovr_Fs)) / sens_co2;
  effl2 = ((mean(txpt)-mean(tctl)-se_t)/(mean(frc)*Ft_ovr_Fs)) / sens_co2;
  effpm2 = effh2-effc2;
end;




