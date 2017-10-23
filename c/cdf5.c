// Purpose: Test behavior of large CDF5 files
// 20170821: Original test by Wei-keng Liou 
// 20170909: Rewritten to expose CDF5 bug by Charlie Zender
// 20170909: A modified version of this has been added to netCDF 

// gcc -std=c99 -I/opt/local/include -o ~/bin/cdf5 ~/sw/c/cdf5.c -L/opt/local/lib -lnetcdf -lhdf5_hl -lhdf5 -lcurl
// gcc -std=c99 -I/usr/local/include -o ~/bin/cdf5 ~/sw/c/cdf5.c -L/usr/local/lib -lnetcdf -lhdf5_hl -lhdf5 -lcurl

#include <stdio.h>
#include <stdlib.h>
#include <netcdf.h>

#define BYT_PER_GiB 1073741824LL

// Test 8 GiB per variable, result: netCDF 4.4.x CDF5 is buggy
//#define DIM (BYT_PER_GiB*8/8)
// Test 5 GiB per variable, result: netCDF 4.4.x CDF5 is buggy
#define DIM (BYT_PER_GiB*5/8)
// Test 3 GiB per variable, result: netCDF 4.4.x CDF5 is fine
//#define DIM (BYT_PER_GiB*3/8)
#define FILENAME "cdf5.nc"

#define ERR {if(err!=NC_NOERR){printf("Error at line %d in %s: %s\n", __LINE__,__FILE__, nc_strerror(err));nerrs++;}}

int main(int argc, char *argv[])
{
  int err, nerrs=0, ncid, dimid, varid[2];
  long long *buf1,*buf2;
  long long avg1,avg2;
  avg1=avg2=0LL;
  
  err = nc_create(FILENAME, NC_CLOBBER|NC_64BIT_DATA, &ncid); ERR;
  //err = nc_create(FILENAME, NC_CLOBBER|NC_NETCDF4, &ncid); ERR;
  err = nc_def_dim(ncid, "dim", DIM, &dimid); ERR;
  // csz 20170830 test record dimension
  // err = nc_def_dim(ncid, "dim", NC_UNLIMITED, &dimid); ERR;
  err = nc_def_var(ncid, "var1", NC_INT64, 1, &dimid, &varid[0]); ERR;
  err = nc_def_var(ncid, "var2", NC_INT64, 1, &dimid, &varid[1]); ERR;
  err = nc_set_fill(ncid, NC_NOFILL, NULL); ERR;
  err = nc_enddef(ncid); ERR;

  buf1 = (long long *) malloc(DIM * sizeof(long long));
  buf2 = (long long *) malloc(DIM * sizeof(long long));
  // 20170831 Write index-dependent values into array so truncation does not yield false-negative answer
  for (long long i=0; i<DIM; i++) buf1[i] = buf2[i] = i;
  err = nc_put_var_longlong(ncid, varid[0], buf1); ERR;
  err = nc_put_var_longlong(ncid, varid[1], buf2); ERR;

  err = nc_close(ncid); ERR;

  err = nc_open(FILENAME, NC_NOWRITE, &ncid); ERR;
  err = nc_inq_varid(ncid, "var1", &varid[0]); ERR;
  err = nc_inq_varid(ncid, "var2", &varid[1]); ERR;
  for (long long i=0; i<DIM; i++) buf1[i] = buf2[i] = 0LL;;
  err = nc_get_var_longlong(ncid, varid[0], buf1); ERR;
  err = nc_get_var_longlong(ncid, varid[1], buf2); ERR;
  err = nc_close(ncid); ERR;

  for (long long i=0; i<DIM; i++) avg1 += buf1[i];
  for (long long i=0; i<DIM; i++) avg2 += buf2[i];
  printf("total1 = %lld, expected1 = %lld\n",avg1,(DIM-1LL)*DIM/2LL);
  printf("total2 = %lld, expected2 = %lld\n",avg2,(DIM-1LL)*DIM/2LL);
  printf("avg1 = %f\n",avg1*1.0/DIM);
  printf("avg2 = %f\n",avg2*1.0/DIM);
  free(buf1);
  free(buf2);

  return (nerrs > 0);
}
