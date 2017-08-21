// Purpose: Test behavior of large CDF5 files

// gcc -std=c99 -I/opt/local/include -o ~/bin/cdf5 ~/sw/c/cdf5.c -L/opt/local/lib -lnetcdf -lhdf5_hl -lhdf5 -lcurl

#include <stdio.h>
#include <stdlib.h>
#include <netcdf.h>

#define DIM 1073741824

#define ERR {if(err!=NC_NOERR){printf("Error at line %d in %s: %s\n", __LINE__,__FILE__, nc_strerror(err));nerrs++;}}

int main(int argc, char *argv[])
{
    int err, nerrs=0, ncid, dimid, varid[3], int_buf1, int_buf2;
    size_t i;
    double *buf, avg=0.0;

    err = nc_create("test.nc", NC_CLOBBER|NC_CDF5, &ncid); ERR
    err = nc_def_dim(ncid, "dim", DIM, &dimid); ERR
    err = nc_def_var(ncid, "var", NC_DOUBLE, 1, &dimid, &varid[0]); ERR
    err = nc_def_var(ncid, "one", NC_INT, 0, NULL, &varid[1]); ERR
    err = nc_def_var(ncid, "two", NC_INT, 0, NULL, &varid[2]); ERR
    err = nc_set_fill(ncid, NC_NOFILL, NULL); ERR
    err = nc_enddef(ncid); ERR

    buf = (double*) malloc(DIM * sizeof(double));
    for (i=0; i<DIM; i++) buf[i] = 1.0;

    err = nc_put_var_double(ncid, varid[0], buf); ERR
    int_buf1 = 1;
    err = nc_put_var_int(ncid, varid[1], &int_buf1); ERR
    int_buf2 = 2;
    err = nc_put_var_int(ncid, varid[2], &int_buf2); ERR
    err = nc_close(ncid); ERR

    err = nc_open("cdf5.nc", NC_NOWRITE, &ncid); ERR
    err = nc_inq_varid(ncid, "var", &varid[0]); ERR
    err = nc_inq_varid(ncid, "one", &varid[1]); ERR
    err = nc_inq_varid(ncid, "two", &varid[2]); ERR
    for (i=0; i<DIM; i++) buf[i] = 0.0;
    err = nc_get_var_double(ncid, varid[0], buf); ERR
    int_buf1 = int_buf2 = 0;
    err = nc_get_var_int(ncid, varid[1], &int_buf1); ERR
    err = nc_get_var_int(ncid, varid[2], &int_buf2); ERR
    err = nc_close(ncid); ERR

    printf("get var one = %d\n",int_buf1);
    printf("get var two = %d\n",int_buf2);
    for (i=0; i<DIM; i++) avg += buf[i];
    avg /= DIM;
    printf("avg = %f\n",avg);
    free(buf);

    return (nerrs > 0);
}
