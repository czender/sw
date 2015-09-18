function [data] = ncread(file, fld_nm);

nc = netcdf(file,'nowrite');
data = nc{fld_nm}(:);
nc = close(nc);
