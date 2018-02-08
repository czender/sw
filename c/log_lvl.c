// Purpose: Test nc_set_log_level()

// gcc -std=c99 -I/opt/local/include -o ~/bin/log_lvl ~/sw/c/log_lvl.c -L/opt/local/lib -lnetcdf -lhdf5_hl -lhdf5 -lcurl
// gcc -std=c99 -I/usr/local/include -o ~/bin/log_lvl ~/sw/c/log_lvl.c -L/usr/local/lib -lnetcdf -lhdf5_hl -lhdf5 -lcurl

#include <stdio.h>
#include <stdlib.h>
#include <netcdf.h>

#define ERR {if(err!=NC_NOERR){printf("Error at line %d in %s: %s\n", __LINE__,__FILE__, nc_strerror(err));nerrs++;}}

int main(int argc, char *argv[])
{
  int ncid;
  nc_set_log_level(3);
  nc_create("file.nc", NC_NETCDF4, &ncid);
  nc_close(ncid);
}

