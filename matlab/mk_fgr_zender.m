% dir         input directory:
% file        input netcdf file:
% fld         data field name
% ttl         figure title
% file_out    output eps file name and path 
% flg_anm     1=anamoly plot, 0=absolute plot
% flg_hl=1    1=set minimum and maximum plot values, 0=auto calculate (has bug)
% min_set     minimum plot value if flg_hl=1
% max_set     maximum plot value if flg_hl=1
% flg_wht     1=minimum value is white, 2=maximum value is white, 3=no white
% lat_bnd     1=plot from latitude -90 to 90, 2=plot from 0 to 90 
% flg_clrbar  1=plot colorbar, 0=no colorbar
% str_frmt    numerical string format
% flg_clr     1=color plot, 0=black and white (in development)
% flg_std     1=computes Students t-test
%             (requires the following two input files and ignores 'file')
%             0=do not compute student's t-test
% file_std1   monthly, zonal-mean timeseries 1
% file_std2   monthly, zonal-mean timeseries 2
% alpha       two-tailed significance level


% figure 1
figure;
dir='/data/mflanner/zender/';
file='cspdfdb07_clm_0112_x.nc';
fld='SOTFRC_SFC';
ttl='Present BC Snow Forcing (W m^-^2)';
file_out='/home/zender/ppr_ZFM07/fgr/fgr_cspdfdb07_x_sotfrcsfc_clr.eps'; 
flg_anm=0;
flg_hl=1;
min_set=0;
max_set=1.5;
flg_wht=1;
lat_bnd=2;
flg_clrbar=1;
str_frmt='%3.2f';
flg_clr=1;

temp = plot_contour_2d(dir,file,fld,ttl,file_out,flg_anm,flg_hl,min_set,max_set,flg_wht,lat_bnd,flg_clrbar,str_frmt,flg_clr);


% figure 2
figure;
dir='/data/mflanner/zender/';
file='cspdfdb07_clm_0112_x.nc';
fld='DSTFRC_SFC';
ttl='Present Dust Snow Forcing (W m^-^2)';
file_out='/home/zender/ppr_ZFM07/fgr/fgr_cspdfdb07_x_dstfrcsfc_clr.eps'; 
flg_anm=0;
flg_hl=1;
min_set=0;
max_set=1.5;
flg_wht=1;
lat_bnd=2;
flg_clrbar=1;
str_frmt='%3.2f';
flg_clr=1;

temp = plot_contour_2d(dir,file,fld,ttl,file_out,flg_anm,flg_hl,min_set,max_set,flg_wht,lat_bnd,flg_clrbar,str_frmt,flg_clr);


% figure 3
figure;
dir='/data/mflanner/zender/';
file='cspdfdb07_clm_0112_x.nc';
fld='AERFRC_SFC';
ttl='Present BC+Dust Snow Forcing (W m^-^2)';
file_out='/home/zender/ppr_ZFM07/fgr/fgr_cspdfdb07_x_aerfrcsfc_clr.eps'; 
flg_anm=0;
flg_hl=1;
min_set=0;
max_set=1.5;
flg_wht=1;
lat_bnd=2;
flg_clrbar=1;
str_frmt='%3.2f';
flg_clr=1;

temp = plot_contour_2d(dir,file,fld,ttl,file_out,flg_anm,flg_hl,min_set,max_set,flg_wht,lat_bnd,flg_clrbar,str_frmt,flg_clr);



% figure 4
figure;
dir='/data/mflanner/zender/';
file='cspdlgm03_clm_0112_x.nc';
fld='SOTFRC_SFC';
ttl='LGM BC Snow Forcing (W m^-^2)';
file_out='/home/zender/ppr_ZFM07/fgr/fgr_cspdlgm03_x_sotfrcsfc_clr.eps'; 
flg_anm=0;
flg_hl=1;
min_set=0;
max_set=1.0;
flg_wht=1;
lat_bnd=2;
flg_clrbar=1;
str_frmt='%3.1f';
flg_clr=1;

temp = plot_contour_2d(dir,file,fld,ttl,file_out,flg_anm,flg_hl,min_set,max_set,flg_wht,lat_bnd,flg_clrbar,str_frmt,flg_clr);


% figure 5
figure;
dir='/data/mflanner/zender/';
file='cspdlgm03_clm_0112_x.nc';
fld='DSTFRC_SFC';
ttl='LGM Dust Snow Forcing (W m^-^2)';
file_out='/home/zender/ppr_ZFM07/fgr/fgr_cspdlgm03_x_dstfrcsfc_clr.eps'; 
flg_anm=0;
flg_hl=1;
min_set=0;
max_set=5.0;
flg_wht=1;
lat_bnd=2;
flg_clrbar=1;
str_frmt='%3.1f';
flg_clr=1;

temp = plot_contour_2d(dir,file,fld,ttl,file_out,flg_anm,flg_hl,min_set,max_set,flg_wht,lat_bnd,flg_clrbar,str_frmt,flg_clr);


% figure 6
figure;
dir='/data/mflanner/zender/';
file='cspdlgm03_clm_0112_x.nc';
fld='AERFRC_SFC';
ttl='LGM BC+Dust Snow Forcing (W m^-^2)';
file_out='/home/zender/ppr_ZFM07/fgr/fgr_cspdlgm03_x_aerfrcsfc_clr.eps'; 
flg_anm=0;
flg_hl=1;
min_set=0;
max_set=5.0;
flg_wht=1;
lat_bnd=2;
flg_clrbar=1;
str_frmt='%3.1f';
flg_clr=1;

temp = plot_contour_2d(dir,file,fld,ttl,file_out,flg_anm,flg_hl,min_set,max_set,flg_wht,lat_bnd,flg_clrbar,str_frmt,flg_clr);


% figure 7
figure;
dir='/data/mflanner/zender/';
file='foo';
fld='TREFHT';
ttl='Pres. BC+Dust/Snow Temp. Change (\circ C)';
file_out='/home/zender/ppr_ZFM07/fgr/fgr_cspdfdb07-cspdfdb08_x_trefht.eps'; 
flg_anm=1;
flg_hl=1;
min_set=-4;
max_set=4;
flg_wht=1;
lat_bnd=2;
flg_clrbar=1;
str_frmt='%3.1f';
flg_clr=1;

flg_std=1;
file_std1='/data/mflanner/zender/cspdfdb07_ts_x.nc';
file_std2='/data/mflanner/zender/cspdfdb08_ts_x.nc';
alpha=0.05;

temp = plot_contour_2d(dir,file,fld,ttl,file_out,flg_anm,flg_hl,min_set,max_set,flg_wht,lat_bnd,flg_clrbar,str_frmt,flg_clr,flg_std,file_std1,file_std2,alpha);


% figure 8
figure;
dir='/data/mflanner/zender/';
file='foo';
fld='QMELT2';
ttl='Pres. BC+Dust/Snow QMELT Change (mm day^-^1)';
file_out='/home/zender/ppr_ZFM07/fgr/fgr_cspdfdb07-cspdfdb08_x_qmelt.eps'; 
flg_anm=1;
flg_hl=1;
min_set=-1;
max_set=1;
flg_wht=1;
lat_bnd=2;
flg_clrbar=1;
str_frmt='%3.2f';
flg_clr=1;

flg_std=1;
file_std1='/data/mflanner/zender/cspdfdb07_ts_x.nc';
file_std2='/data/mflanner/zender/cspdfdb08_ts_x.nc';
alpha=0.05;

temp = plot_contour_2d(dir,file,fld,ttl,file_out,flg_anm,flg_hl,min_set,max_set,flg_wht,lat_bnd,flg_clrbar,str_frmt,flg_clr,flg_std,file_std1,file_std2,alpha);


% figure 9
figure;
dir='/data/mflanner/zender/';
file='foo';
fld='ALBS';
ttl='Pres. BC+Dust/Snow ALBS Change';
file_out='/home/zender/ppr_ZFM07/fgr/fgr_cspdfdb07-cspdfdb08_x_albs.eps'; 
flg_anm=1;
flg_hl=1;
min_set=-0.2;
max_set=0.2;
flg_wht=1;
lat_bnd=2;
flg_clrbar=1;
str_frmt='%3.2f';
flg_clr=1;

flg_std=1;
file_std1='/data/mflanner/zender/cspdfdb07_ts_x.nc';
file_std2='/data/mflanner/zender/cspdfdb08_ts_x.nc';
alpha=0.05;

temp = plot_contour_2d(dir,file,fld,ttl,file_out,flg_anm,flg_hl,min_set,max_set,flg_wht,lat_bnd,flg_clrbar,str_frmt,flg_clr,flg_std,file_std1,file_std2,alpha);



% figure 10
figure;
dir='/data/mflanner/zender/';
file='foo';
fld='TREFHT';
ttl='LGM BC+Dust/Snow Temp. Change (\circ C)';
file_out='/home/zender/ppr_ZFM07/fgr/fgr_cspdlgm03-cspdlgm04_x_trefht.eps'; 
flg_anm=1;
flg_hl=1;
min_set=-4;
max_set=4;
flg_wht=1;
lat_bnd=2;
flg_clrbar=1;
str_frmt='%3.1f';
flg_clr=1;

flg_std=1;
file_std1='/data/mflanner/zender/cspdlgm03_ts_x.nc';
file_std2='/data/mflanner/zender/cspdlgm04_ts_x.nc';
alpha=0.05;

temp = plot_contour_2d(dir,file,fld,ttl,file_out,flg_anm,flg_hl,min_set,max_set,flg_wht,lat_bnd,flg_clrbar,str_frmt,flg_clr,flg_std,file_std1,file_std2,alpha);



% figure 11
figure;
dir='/data/mflanner/zender/';
file='foo';
fld='QMELT2';
ttl='LGM BC+Dust/Snow QMELT Change (mm day^-^1)';
file_out='/home/zender/ppr_ZFM07/fgr/fgr_cspdlgm03-cspdlgm04_x_qmelt.eps'; 
flg_anm=1;
flg_hl=1;
min_set=-1;
max_set=1;
flg_wht=1;
lat_bnd=2;
flg_clrbar=1;
str_frmt='%3.2f';
flg_clr=1;

flg_std=1;
file_std1='/data/mflanner/zender/cspdlgm03_ts_x.nc';
file_std2='/data/mflanner/zender/cspdlgm04_ts_x.nc';
alpha=0.05;

temp = plot_contour_2d(dir,file,fld,ttl,file_out,flg_anm,flg_hl,min_set,max_set,flg_wht,lat_bnd,flg_clrbar,str_frmt,flg_clr,flg_std,file_std1,file_std2,alpha);



% figure 12
figure;
dir='/data/mflanner/zender/';
file='foo';
fld='ALBS';
ttl='LGM BC+Dust/Snow ALBS Change';
file_out='/home/zender/ppr_ZFM07/fgr/fgr_cspdlgm03-cspdlgm04_x_albs.eps'; 
flg_anm=1;
flg_hl=1;
min_set=-0.2;
max_set=0.2;
flg_wht=1;
lat_bnd=2;
flg_clrbar=1;
str_frmt='%3.2f';
flg_clr=1;

flg_std=1;
file_std1='/data/mflanner/zender/cspdlgm03_ts_x.nc';
file_std2='/data/mflanner/zender/cspdlgm04_ts_x.nc';
alpha=0.05;

temp = plot_contour_2d(dir,file,fld,ttl,file_out,flg_anm,flg_hl,min_set,max_set,flg_wht,lat_bnd,flg_clrbar,str_frmt,flg_clr,flg_std,file_std1,file_std2,alpha);

