% Computes band-averaged Mie parameters (mass extinction coefficient,
% single-scatter albedo, and asymetry parameter) by weighting
% according to surface-incident solar flux and, for ss_alb,
% semi-infinite albedo (from Wiscombe and Warren, 1980). The
% semi-infinite albedo depends strongly (and non-linearly) on the
% single-scattering albedo. Hence, weighting ss_alb also by this
% albedo predicts broadband ss_alb that predicts albedo quite well.
% Typical errors for 5-band approximation are ~5% albedo, ~10% sub-2cm
% absorption. See flx_abs_comparison.txt for detailed info.

% REQUIRED INPUT: 470-band Mie parameter file, computed at 10-nm
% resolution from 0.3-5.0 microns.

% Usage: 
% scp ~/matlab/mie_brd_bnd_cnv.m biogenic.ess.uci.edu:matlab
% scp ~/matlab/mie_brd_bnd_cnv.m esmf.ess.uci.edu:matlab
% matlab -nojvm
% mie_brd_bnd_cnv
% cp ${DATA}/dst/snicar/aer_dst_bln*.nc /ptmp/zender/inputdata_cam/snicar
% rsync 'esmf.ess.uci.edu:${DATA}/dst/snicar' /ptmp/zender/inputdata_cam

clear;

% flag: 
% =1: weight by clear-sky surface-incident flux
% =2: weight by cloudy-sky surface-incident flux
flg_spc = 1;


% flag:
% =1: only weight ice parameters
% =2: only weight soot parameters (no ice)
% =3: only weight dust parameters (no ice)
flg_aer = 3;
%fl_in_aer='/data/zender/dst/mie/aer_SiO2_Fe2O3_doccd_03_snicar.nc';
%fl_out_aer='/data/zender/dst/snicar/aer_SiO2_Fe2O3_doccd_03_snicar.nc';
fl_in_aer='/data/zender/dst/mie/aer_dst_bln_20060904_04_snicar.nc';
fl_out_aer='/data/zender/dst/snicar/aer_dst_bln_20060904_04_snicar.nc';


% flag:
% =0: uncoated. use ext_cff_mss
% =1: coated: use ext_cff_mss_cor
flg_cot = 0;


wrkdir='/data/mflanner/mie/mie6/';
write_dir='/data/mflanner/mie_clm/mie_clm17/';
fl_stb1='mie_ice_';
fl_stb2='.0.nc';


% Starting and ending grain radii (microns)
r1=30;
r2=1000;


% wavelength beginning and ending indecies for each band:

%bnd_min_idx = [1 41 91];
%bnd_max_idx = [40 90 470];

%bnd_min_idx = [1 41 81 111];
%bnd_max_idx = [40 80 110 470];

%bnd_min_idx = [1 41 81 111 151];
%bnd_max_idx = [40 80 110 150 470];

%bnd_min_idx = [1 41 81 111 151];
%bnd_max_idx = [40 80 110 150 470];

%bnd_min_idx = [1 41 71 91 131];
%bnd_max_idx = [40 70 90 130 470];

% BANDS:
% 1) 0.3-0.7 um
% 2) 0.7-1.0 um
% 3) 1.0-1.2 um
% 4) 1.2-1.5 um
% 5) 1.5-5.0 um
bnd_min_idx = [1 41 71 91 121];
bnd_max_idx = [40 70 90 120 470];

%bnd_min_idx = [1 41 71 121];
%bnd_max_idx = [40 70 120 470];



%wvl(1)=0.5E-6;
%wvl(2)=2.85E-6;

%wvl(1)=0.5E-6;
%wvl(2)=0.95E-6;
%wvl(3)=3.1E-6;

%wvl(1)=0.5E-6;
%wvl(2)=0.9E-6;
%wvl(3)=1.25E-6;
%wvl(4)=3.2E-6;

%wvl(1)=0.5E-6;
%wvl(2)=0.9E-6;
%wvl(3)=1.25E-6;
%wvl(4)=1.6E-6;
%wvl(5)=3.4E-6;

wvl(1)=0.5E-6;
wvl(2)=0.85E-6;
wvl(3)=1.1E-6;
wvl(4)=1.35E-6;
wvl(5)=3.25E-6;

%wvl(1)=0.5E-6;
%wvl(2)=0.85E-6;
%wvl(3)=1.25E-6;
%wvl(4)=3.25E-6;


date_str=date;

if (flg_aer == 1)
  for i=r1:r2
    i
    % Read snow Mie parameters
    if ((i<100) & (i>=30))
      s1=int2str(0);
      s2=int2str(i);
      fl=strcat(fl_stb1,s1,s1,s2,fl_stb2);
      fl_in=strcat(wrkdir,fl);
    
    elseif ((i>=100) & (i<1000))
      s1=int2str(0);
      s2=int2str(i);
      fl=strcat(fl_stb1,s1,s2,fl_stb2);
      fl_in=strcat(wrkdir,fl);
    elseif (i==1000)
      s2=int2str(i);
      fl=strcat(fl_stb1,s2,fl_stb2);
      fl_in=strcat(wrkdir,fl);
    end;
    
    
    omega_snw=       ncread(fl_in,'ss_alb');
    ext_cff_mss_snw= ncread(fl_in,'ext_cff_mss');
    g_snw=           ncread(fl_in,'asm_prm');
    %  flx_slr_frc=     ncread(fl_in,'flx_slr_frc');
  
    if (flg_spc == 1)
      load mlw_sfc_flx_frc_clr.txt;
      flx_slr_frc = mlw_sfc_flx_frc_clr;
    elseif (flg_spc == 2)
      load mlw_sfc_flx_frc_cld.txt;
      flx_slr_frc = mlw_sfc_flx_frc_cld;
    end; 
    
    % get semi-infinite albedo from Wiscombe and Warren, 1980
    alb_inf_dfs=     alb_ww(i);
  

    for j=1:length(bnd_max_idx)
      % single-scatter albedo:
      % ORIGINAL METHOD: weight by solar flux:
      % omega(j) = sum(flx_slr_frc(bnd_min_idx(j):bnd_max_idx(j)).*omega_snw(bnd_min_idx(j):bnd_max_idx(j)))/sum(flx_slr_frc(bnd_min_idx(j):bnd_max_idx(j)));
      
      % NEW METHOD: weight by solar flux and semi-infinite albedo
      omega(j) = sum(flx_slr_frc(bnd_min_idx(j):bnd_max_idx(j)).*alb_inf_dfs(bnd_min_idx(j):bnd_max_idx(j)).*omega_snw(bnd_min_idx(j):bnd_max_idx(j)))/sum(flx_slr_frc(bnd_min_idx(j):bnd_max_idx(j)).*alb_inf_dfs(bnd_min_idx(j):bnd_max_idx(j)));
      
      % asymmetry parameter: weight by solar flux
      asm_prm(j) = sum(flx_slr_frc(bnd_min_idx(j):bnd_max_idx(j)).*g_snw(bnd_min_idx(j):bnd_max_idx(j)))/sum(flx_slr_frc(bnd_min_idx(j):bnd_max_idx(j)));
      
      % mass extinction coefficient: weight by solar flux
      ext_cff_mss(j) = sum(flx_slr_frc(bnd_min_idx(j):bnd_max_idx(j)).*ext_cff_mss_snw(bnd_min_idx(j):bnd_max_idx(j)))/sum(flx_slr_frc(bnd_min_idx(j):bnd_max_idx(j)));
      
    end
    
   
  
  
    % --------Write to output file---------


    fl_out = strcat(write_dir,fl);
  
    % Create NetCDF file.
    nc = netcdf(fl_out, 'clobber');              

    % Global attributes.
    nc.description = 'Broadband Mie parameters for ice spheres. Derived from 470 spectral bands, weighted by clear-sky mid-latitude winter surface incident flux. Single-scattering albedo also weighted by semi-infinite albedo from Wiscombe and Warren, 1980.';
    nc.author = 'Mark Flanner';
    nc.date = date_str;
  
    % Define dimensions:
    nc('wvl') = length(wvl);
    
    % Define variables.
    nc{'wvl'} = 'wvl';
    nc{'ss_alb'} = 'wvl';
    nc{'ext_cff_mss'} = 'wvl';
    nc{'asm_prm'} = 'wvl';

    
    % Define attributes
    nc{'wvl'}.units = 'meters';
    nc{'ss_alb'}.units = 'Fraction';      
    nc{'ext_cff_mss'}.units = 'meter2 kilogram-1';
    nc{'asm_prm'}.units = 'unitless (cos theta)';
    
    % Store the data
    nc{'wvl'}(:) = wvl;
    nc{'ss_alb'}(:) = omega;
    nc{'ext_cff_mss'}(:) = ext_cff_mss;
    nc{'asm_prm'}(:) = asm_prm;
    
    nc = close(nc);                                      % Close the file.
    
  end

  
  
%%%%%%%  SOOT  %%%%%%%%
elseif (flg_aer==2)
  
  
  %wrkdir='/data/mflanner/mie/';

  %fl='mie_sot_0000.1_um_vma.nc';
  %fl='mie_sot_0000.0118_um_nma.nc';
  %fl='mie_sot_OPAC.nc';

  %fl='miecot_slfsot_ChC90_dns_1568.nc';
  %fl_in=strcat(wrkdir,fl);
  
  fl_in = fl_in_aer;
  

  % NOTE: USE "ext_cff_mss_cor" FOR COATED SPHERES!
  omega_sot=       ncread(fl_in,'ss_alb');
  if (flg_cot==0)
    ext_cff_mss_dst= ncread(fl_in,'ext_cff_mss');
  elseif(flg_cot==1)
    ext_cff_mss_dst= ncread(fl_in,'ext_cff_mss_cor');
  end;
  g_sot=           ncread(fl_in,'asm_prm');
  %flx_slr_frc=     ncread(fl_in,'flx_slr_frc');
  

  % load weighted clear/cloudy sky incident flux:
  load mlw_sfc_flx_frc_cmb.txt;
  flx_slr_frc = mlw_sfc_flx_frc_cmb;
  

  % weight all Mie parameters by surface incident flux:
  for i=1:length(bnd_max_idx)
    omega(i) = sum(flx_slr_frc(bnd_min_idx(i):bnd_max_idx(i)).*omega_sot(bnd_min_idx(i):bnd_max_idx(i)))/sum(flx_slr_frc(bnd_min_idx(i):bnd_max_idx(i)));
    
    asm_prm(i) = sum(flx_slr_frc(bnd_min_idx(i):bnd_max_idx(i)).*g_sot(bnd_min_idx(i):bnd_max_idx(i)))/sum(flx_slr_frc(bnd_min_idx(i):bnd_max_idx(i)));
      
    ext_cff_mss(i) = sum(flx_slr_frc(bnd_min_idx(i):bnd_max_idx(i)).*ext_cff_mss_sot(bnd_min_idx(i):bnd_max_idx(i)))/sum(flx_slr_frc(bnd_min_idx(i):bnd_max_idx(i)));
    
  end
   
  
  
  % --------Write to output file---------

  %write_dir='/data/mflanner/mie_clm2/';
  %fl_out = strcat(write_dir,fl);
  
  fl_out=fl_out_aer;
 
  % Create NetCDF file.
  nc = netcdf(fl_out, 'clobber');              

  % Global attributes.
  nc.description = 'Broadband Mie parameters for black carbon. Derived from 470 spectral bands, weighted by 70% clear-sky/30% cloudy-sky mid-latitude winter surface incident flux.';
  nc.author = 'Mark Flanner';
  nc.date = date_str;
  
  % Define dimensions:
  nc('wvl') = length(wvl);

  % Define variables.
  nc{'wvl'} = 'wvl';
  nc{'ss_alb'} = 'wvl';
  nc{'ext_cff_mss'} = 'wvl';
  nc{'asm_prm'} = 'wvl';


  % Define attributes
  nc{'wvl'}.units = 'meters';
  nc{'ss_alb'}.units = 'fraction';      
  nc{'ext_cff_mss'}.units = 'meter2 kilogram-1';
  nc{'asm_prm'}.units = 'unitless (cos theta)';

  % Store the data
  
  nc{'wvl'}(:) = wvl;
  nc{'ss_alb'}(:) = omega;
  nc{'ext_cff_mss'}(:) = ext_cff_mss;
  nc{'asm_prm'}(:) = asm_prm;

  nc = close(nc);                 
  

  
%%%%%%%  DUST  %%%%%%%%%
elseif (flg_aer==3)

  %wrkdir='/data/mflanner/mie/';
  %fl='miecot_slfsot_ChC90_dns_1568.nc';
  %fl_in=strcat(wrkdir,fl);

  fl_in=fl_in_aer;
  

  % NOTE: USE "ext_cff_mss_cor" FOR COATED SPHERES!
  omega_dst=       ncread(fl_in,'ss_alb');
  if (flg_cot==0)
    ext_cff_mss_dst= ncread(fl_in,'ext_cff_mss');
  elseif(flg_cot==1)
    ext_cff_mss_dst= ncread(fl_in,'ext_cff_mss_cor');
  end;
  g_dst=           ncread(fl_in,'asm_prm');
  %flx_slr_frc=     ncread(fl_in,'flx_slr_frc');


  % load weighted clear/cloudy sky incident flux:
  load mlw_sfc_flx_frc_cmb.txt;
  flx_slr_frc = mlw_sfc_flx_frc_cmb;


  % weight all Mie parameters by surface incident flux:
  for i=1:length(bnd_max_idx)
    omega(i) = sum(flx_slr_frc(bnd_min_idx(i):bnd_max_idx(i)).*omega_dst(bnd_min_idx(i):bnd_max_idx(i)))/sum(flx_slr_frc(bnd_min_idx(i):bnd_max_idx(i)));
  
    asm_prm(i) = sum(flx_slr_frc(bnd_min_idx(i):bnd_max_idx(i)).*g_dst(bnd_min_idx(i):bnd_max_idx(i)))/sum(flx_slr_frc(bnd_min_idx(i):bnd_max_idx(i)));
    
    ext_cff_mss(i) = sum(flx_slr_frc(bnd_min_idx(i):bnd_max_idx(i)).*ext_cff_mss_dst(bnd_min_idx(i):bnd_max_idx(i)))/sum(flx_slr_frc(bnd_min_idx(i):bnd_max_idx(i)));
  end
   
  
  
  % --------Write to output file---------

  %write_dir='/data/mflanner/mie_clm2/';
  %fl_out = strcat(write_dir,fl);
 
  fl_out=fl_out_aer;
 
  % Create NetCDF file.
  nc = netcdf(fl_out, 'clobber');              

  % Global attributes.
  nc.description = 'Broadband Mie parameters for dust. Derived from 470 spectral bands, weighted by 70% clear-sky/30% cloudy-sky mid-latitude winter surface incident flux.';
  nc.author = 'Mark Flanner';
  nc.date = date_str;
  
  % Define dimensions:
  nc('wvl') = length(wvl);

  % Define variables.
  nc{'wvl'} = 'wvl';
  nc{'ss_alb'} = 'wvl';
  nc{'ext_cff_mss'} = 'wvl';
  nc{'asm_prm'} = 'wvl';


  % Define attributes
  nc{'wvl'}.units = 'meters';
  nc{'ss_alb'}.units = 'fraction';      
  nc{'ext_cff_mss'}.units = 'meter2 kilogram-1';
  nc{'asm_prm'}.units = 'unitless (cos theta)';

  % Store the data
  
  nc{'wvl'}(:) = wvl;
  nc{'ss_alb'}(:) = omega;
  nc{'ext_cff_mss'}(:) = ext_cff_mss;
  nc{'asm_prm'}(:) = asm_prm;

  nc = close(nc);                 
  
  
  
end;
