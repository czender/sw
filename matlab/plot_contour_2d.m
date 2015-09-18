
function min_min = plot_contour_2d(dir,file,fld,ttl,file_out,flg_anm,flg_hl,min_set,max_set,flg_wht,lat_bnd,flg_clrbar,str_frmt,flg_clr,flg_std,file_std1,file_std2,alpha)


% set all input parameters manually:
if (1==0)
  clear;
  
  dir='/data/mflanner/dif_camsombc_1998cnt5_camsombc_1998atm2/';
  file='camsombc_1998cnt5-camsombc_1998atm2_clm_0112_x.nc';
  fld='QMELT';

  %dir='/data/mflanner/anl_camsombc_1998cnt5/';
  %file='camsombc_1998cnt5_clm_0112_x.nc';
  %fld='SOTFRC_SFC';

  ttl=strcat('1998 Mid BC Snow - Control, 2m T (\circ C)');
  file_out='/home/mflanner/ppr_FZR06/foo.eps';
  
  % flg_anm: 0= normal (rainbow) plot; 1= anamoly plot (blues-reds)
  flg_anm=1;
  
  % flg_hl: 1= manually set high and low, 0= auto
  flg_hl=1;
  min_set=-2;
  max_set=2;
  
  % flg_wht: set white color for: 1=min; 2=max; or 3=none
  flg_wht=3;
  
  % lat_bnd: plot boundaries: 1= -90->90; 2= 0->90
  lat_bnd=2;
  
  % flg_clrbar: 0= don't plot colorbar; 1=plot
  flg_clrbar=1;
  str_frmt='%3.2f'; % for tick labels on colorbar

  flg_std=1;
  file_std1='/data/mflanner/anl_camsombc_1998cnt5/camsombc_1998cnt5_ts_x.nc';
  file_std2='/data/mflanner/anl_camsombc_1998atm2/camsombc_1998atm2_ts_x.nc';
  alpha=0.05;

end;
  
if (nargin < 15)
  flg_std=0;
end;

% retrieve data:
if (flg_std==0)
  fl=strcat(dir,file);

  nc = netcdf(fl,'nowrite');
  data_dif = nc{fld}(:,:,:);
  lat = nc{'lat'}(:);
  nc=close(nc);
  
  data_dif = squeeze(data_dif);
  data_dif = data_dif';
  data_dif(data_dif>1e30)=0;
  
elseif (flg_std==1)
  nc = netcdf(file_std1,'nowrite');
  data_std1 = nc{fld}(:,:,:);
  lat = nc{'lat'}(:);
  nc=close(nc);
  
  nc = netcdf(file_std2,'nowrite');
  data_std2 = nc{fld}(:,:,:);
  lat = nc{'lat'}(:);
  nc=close(nc);
  
  if (length(data_std1) ~= length(data_std2))
    error('student T-test arrays not of equal length');
  end;
  
  % reshape input std arrays to [year,month,lat]
  for i=1:length(data_std1)
    j=mod(i,12);
    if (j==0)
      j=12;
    end;
    k=ceil(i/12);
    
    data_std1a(k,j,:)=data_std1(i,:);
    data_std2a(k,j,:)=data_std2(i,:);
  end
  
  % perform the unpaired t-test, storing ones and zeroes in matrix std.
  std=ttest2(data_std1a,data_std2a,alpha); 
  std=squeeze(std);
  std=std';
  
  % calcualte the difference, for plotting;
  data_dif=mean(data_std1a,1)-mean(data_std2a,1);
  data_dif = squeeze(data_dif);
  data_dif = data_dif';
  data_dif(data_dif>1e30)=0;
     
end;

% set x dimension
time=[1:12];

% set min and max values for plotting
if (flg_hl==0)
  min1 = min(min(data_dif));
  max1 = max(max(data_dif));

  if (flg_anm==0)
    min_min = min1;
    max_max = max1;
  else
    t1=abs(min1);
    if (t1 > max1)
      max_max = abs(min1);
      min_min = min1;
    else
      max_max = max1;
      min_min = -max1;
    end;
  end;
  
else
  min_min = min_set;
  max_max = max_set;
end;

% set contour levels:
% set first contour level REALLY low, so that all data is
% 'filled' in. (Shading is only done ABOVE each contour value)
if (flg_anm==0)
  clevs_int = (max_max-min_min) / 10;
  v=[min_min:clevs_int:max_max];
  v(1) = min_min-30*clevs_int;
elseif (flg_anm==1)
  clevs_int = (max_max-min_min) / 9;
  v=[min_min-clevs_int:clevs_int:max_max];
  v(1) = min_min-30*clevs_int;
 
  if (1==0)
    v(1)=min_min-clevs_int;
    v(2)=v(1)+clevs_int;
    v(3)=v(2)+clevs_int;
    v(4)=v(3)+clevs_int;
    v(5)=v(4)+clevs_int;
    v(6)=v(5)+clevs_int;
    v(7)=v(6)+clevs_int;
    v(8)=v(7)+clevs_int;
    v(9)=v(8)+clevs_int;
    v(10)=v(9)+clevs_int;
    v(11)=v(10)+clevs_int;
  end;
end;

% plot contours
[temp hp] = contourf(time,lat,data_dif,v);

% set colormap
define_colors;

if (flg_anm==0)
  if (flg_wht==1)
    colormap(clr_rnb1b);
  elseif(flg_wht==2)
    colormap(clr_rnb2);
  elseif(flg_wht==3)
    %colormap(clr_rnb3);
    colormap(clr_rnb3);
  end;
elseif (flg_anm==1)
    colormap(clr_anm2);
end;

% set bounds of colormap to match contours:
% lower bound set from min_min, since v(1) is really low.
if (flg_anm==0)
  %caxis([v(1) v(end)+clevs_int-clevs_int/1000]);
  caxis([min_min v(end)+clevs_int-clevs_int/1000]);
  caxis('manual');
elseif(flg_anm==1)
  %caxis([v(1) v(end)+clevs_int-clevs_int/1000]);
  caxis([min_min-clevs_int v(end)+clevs_int-clevs_int/1000]);
  caxis('manual');
end;

% fine-tune the plot:

if (lat_bnd==1)
  axis([1 12 lat(1) lat(end)]);
elseif (lat_bnd==2)
  axis([1 12 0 lat(end)]);
end;

set(gca,'XTick',1:1:12);
set(gca,'XTickLabel',{'J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D'},'fontsize',20);

if (lat_bnd==1)
  set(gca,'YTick',-80:20:80,'fontsize',20);
elseif(lat_bnd==2)
  set(gca,'YTick',0:10:80,'fontsize',20);
end;

if (flg_std==1)
  k=1;
  [row,clm]=find(std==1);
  for i=1:length(row)
    hold on;
    hstd(k)=plot(clm(i),lat(row(i)),'kx');
    k=k+1;
  end
  set(hstd(:),'MarkerSize',20);
end;


xlabel('Month','fontsize',22);
ylabel('Latitude ({\circ}N)','fontsize',22);
title(ttl,'fontsize',20);
grid on;

% make colorbar:
if (flg_clrbar==1)
  hcb = colorbar;
  
  cblim=get(hcb,'YLim');
  tickint=(cblim(2)-cblim(1))/10;
  ticks=[cblim(1)+tickint/2:tickint:cblim(2)-tickint/2];
  
  for i=1:length(ticks)
    temp = sprintf(str_frmt,ticks(i));
    ticks2(i,1:length(temp)) = temp;
  end
  
  set(hcb,'YTick',ticks);
  set(hcb,'YTickLabel',ticks2);
  set(hcb,'YTickMode','manual');
  set(hcb,'fontsize',18);
end;
 
% save as eps file:
saveas(gcf,file_out,'epsc');
