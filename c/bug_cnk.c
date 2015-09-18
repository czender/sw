/* $Id$ */

/* Purpose: Demonstrate netCDF4 chunking behavior */

/* Usage:
   cd ~/c;./bug_cnk
   cd ~/c;gcc -std=c99 -o bug_cnk bug_cnk.c -L/usr/local/lib -lnetcdf -lhdf5_hl -lhdf5 -lcurl */

#include <stdio.h>
#include <netcdf.h> /* netCDF definitions and C library */

  /* Glue code */
#define FILE_NAME "./bug_cnk.nc"
  int ncid; // [id] netCDF file ID

#define ERR do { \
fflush(stdout); /* Make sure our stdout is synced with stderr. */ \
err++; \
fprintf(stderr, "Sorry! Unexpected result, %s, line: %d\n", \
        __FILE__, __LINE__);                                \
} while (0)
int err=0; // global

int main(){
  /* Code from netCDF4 C Users Manual p. 90 */
#define NDIMS6 1
#define DIM6_NAME "D5"
#define VAR_NAME6 "V5"
#define DIM6_LEN 100
  int dimids[NDIMS6], dimids_in[NDIMS6];
  int varid;
  int ndims, nvars, natts, unlimdimid;
  nc_type xtype_in;
  char name_in[NC_MAX_NAME + 1];
  int data[DIM6_LEN], data_in[DIM6_LEN];
  size_t chunksize_in[NDIMS6];
  int storage_in;
  int i, d;
  for (i = 0; i < DIM6_LEN; i++)
    data[i] = i;
  /* Create a netcdf-4 file with one dim and one var. */
  if (nc_create(FILE_NAME, NC_NETCDF4, &ncid)) ERR;
  if (nc_def_dim(ncid, DIM6_NAME, NC_UNLIMITED, &dimids[0])) ERR;
  // if (nc_def_dim(ncid, DIM6_NAME, DIM6_LEN, &dimids[0])) ERR;
  if (dimids[0] != 0) ERR;
  if (nc_def_var(ncid, VAR_NAME6, NC_INT, NDIMS6, dimids, &varid)) ERR;
  //if (nc_def_var_chunking(ncid, varid, NC_CONTIGUOUS, NULL)) ERR;
  if (nc_put_var_int(ncid, varid, data)) ERR;
  /* Check stuff. */
  if (nc_inq_var_chunking(ncid, 0, &storage_in, chunksize_in)) ERR;
   if (storage_in != NC_CONTIGUOUS) ERR;
}
