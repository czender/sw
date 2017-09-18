// Purpose: Test behavior of large CDF5 files
// 20170821: Original test by Wei-keng Liou 
// 20170909: Rewritten to expose CDF5 bug by Charlie Zender
// 20170918: Rewritten with PnetCDF to further test CDF5 bug

// scp ~/sw/c/pcdf5.c skyglow.ess.uci.edu:sw/c

// mpicc -std=c99 -I/usr/local/parallel/include -o ~/bin/pcdf5 ~/sw/c/pcdf5.c -L/usr/local/parallel/lib -L/usr/local/parallel/lib -L/usr/lib64/hdf -lnetcdf -lpnetcdf -ljpeg -lmfhdf -ldf -lhdf5_hl -lhdf5 -ldl -lm -lz -lcurl -ljpeg # Skyglow
  
#include <stdio.h>
#include <stdlib.h>
#include <netcdf.h>
#include <mpi.h>
#include <netcdf.h>
#include <netcdf_par.h>

#define BYT_PER_GiB 1073741824LL
//Test 5 GiB per variable, result: netCDF 4.4.x CDF5 is buggy
#define DIM (BYT_PER_GiB*5/8)

#define FILENAME "pcdf5.nc"

static void 
report(int rcd, char* file, int line){
  fflush(stdout); /* Make sure our stdout is synced with stderr.*/
  if(rcd != 0) fprintf(stderr, "Sorry! Unexpected result, %s, line: %d: status=%d\n",__FILE__,__LINE__,rcd);
}
#define Error(rcd) {if(rcd!=NC_NOERR){report(rcd,__FILE__,__LINE__);nerrs++;}}

int write(int ncid, int parallel)
{
  int err, nerrs=0, cmode, varid[2], dimid, rcd;
  size_t start[10], count[10];
  long long *buf1,*buf2;
  
  rcd=nc_def_dim(ncid, "dim", DIM, &dimid); Error(rcd);
  rcd=nc_def_var(ncid, "var1",NC_INT64, 1, &dimid, &varid[0]); Error(rcd);
  rcd=nc_def_var(ncid, "var2",NC_INT64, 1, &dimid, &varid[1]); Error(rcd);
  rcd=nc_set_fill(ncid, NC_NOFILL, NULL); Error(rcd);
  rcd=nc_enddef(ncid); Error(rcd);

  if(parallel){
    /* NC_INDEPENDENT is default */
    rcd=nc_var_par_access(ncid, varid[0], NC_INDEPENDENT); Error(rcd);
    rcd=nc_var_par_access(ncid, varid[1], NC_INDEPENDENT); Error(rcd);
  } /* !parallel */

  buf1 = (long long *) malloc(DIM * sizeof(long long));
  buf2 = (long long *) malloc(DIM * sizeof(long long));
  // 20170831 Write index-dependent values into array so truncation does not yield false-negative answer
  for (long long i=0; i<DIM; i++) buf1[i] = buf2[i] = i;

  if(parallel){
    /* Parallel write in 1 GiB chunks */
    int nbr_gib;
    nbr_gib=DIM*sizeof(long long)/BYT_PER_GiB;
    printf("nbr_gib = %d\n",nbr_gib);
    for(int i=0;i<nbr_gib;i++){
      start[i]=i*DIM*8/nbr_gib;
      count[0]=DIM*8/nbr_gib;
      rcd=nc_put_vara_longlong(ncid, varid[0], &start[i], &count[i], buf1); Error(rcd);
      rcd=nc_put_vara_longlong(ncid, varid[1], &start[i], &count[i], buf2); Error(rcd);
    } /* !i */
  }else{
    err = nc_put_var_longlong(ncid, varid[0], buf1); Error(rcd);
    err = nc_put_var_longlong(ncid, varid[1], buf2); Error(rcd);
  }

  free(buf1);
  free(buf2);

  return NC_NOERR;
} /* !write() */

int read(int ncid, int parallel)
{
  int err, nerrs=0, cmode, varid[2], dimid, rcd;
  size_t start[10], count[10];
  long long *buf1,*buf2;
  long long avg1,avg2;
  avg1=avg2=0LL;

  buf1 = (long long *) malloc(DIM * sizeof(long long));
  buf2 = (long long *) malloc(DIM * sizeof(long long));
  for (long long i=0; i<DIM; i++) buf1[i] = buf2[i] = 0LL;

  err = nc_inq_varid(ncid, "var1", &varid[0]); Error(rcd);
  err = nc_inq_varid(ncid, "var2", &varid[1]); Error(rcd);

  if(parallel){
    /* Parallel read in 1 GiB chunks */
    int nbr_gib;
    nbr_gib=DIM*sizeof(long long)/BYT_PER_GiB;
    printf("nbr_gib = %d\n",nbr_gib);
    for(int i=0;i<nbr_gib;i++){
      start[i]=i*DIM*8/nbr_gib;
      count[0]=DIM*8/nbr_gib;
      rcd=nc_get_vara_longlong(ncid, varid[0], &start[i], &count[i], buf1); Error(rcd);
      rcd=nc_get_vara_longlong(ncid, varid[1], &start[i], &count[i], buf2); Error(rcd);
    } /* !i */
  }else{
    err = nc_get_var_longlong(ncid, varid[0], buf1); Error(rcd);
    err = nc_get_var_longlong(ncid, varid[1], buf2); Error(rcd);
  }

  for (long long i=0; i<DIM; i++) avg1 += buf1[i];
  for (long long i=0; i<DIM; i++) avg2 += buf2[i];
  printf("total1 = %lld, expected1 = %lld\n",avg1,(DIM-1LL)*DIM/2LL);
  printf("total2 = %lld, expected2 = %lld\n",avg2,(DIM-1LL)*DIM/2LL);
  printf("avg1 = %f\n",avg1*1.0/DIM);
  printf("avg2 = %f\n",avg2*1.0/DIM);

  free(buf1);
  free(buf2);

  return NC_NOERR;
} /* !read() */

int main(int argc, char* argv[])
{
  int rank, nprocs, nerrs=0, ncid, cmode, rcd;
  
  MPI_Comm comm=MPI_COMM_SELF;
  MPI_Info info=MPI_INFO_NULL;
  
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(nprocs > 1 && rank == 0) printf("This test program is intended to run on ONE process\n");
  if(rank > 0) goto fn_exit;
  
#ifdef DISABLE_PNETCDF_ALIGNMENT
  MPI_Info_create(&info);
  MPI_Info_set(info, "nc_header_align_size", "1");
  MPI_Info_set(info, "nc_var_align_size",    "1");
#endif
  
  /* PnetCDF->CDF5 */
  printf("\nWrite using PnetCDF; Read using CDF5\n");
  cmode = NC_PNETCDF | NC_CLOBBER;
  rcd=nc_create_par(FILENAME, cmode, comm, info, &ncid); Error(rcd);
  rcd=write(ncid,1); Error(rcd);
  rcd=nc_close(ncid); Error(rcd);
  
  cmode = NC_CDF5 | NC_NOCLOBBER;
  rcd=nc_open(FILENAME, cmode, &ncid); Error(rcd);
  rcd=read(ncid,0); Error(rcd);
  rcd=nc_close(ncid); Error(rcd);
  
  /* CDF5->PnetCDF */
  printf("\nWrite using CDF5; Read using PnetCDF\n");
  cmode = NC_CDF5 | NC_CLOBBER;
  rcd=nc_create(FILENAME, cmode, &ncid); Error(rcd);
  rcd=write(ncid,0); Error(rcd);
  rcd=nc_close(ncid); Error(rcd);
  
  cmode = NC_PNETCDF | NC_NOCLOBBER;
  rcd=nc_open_par(FILENAME, cmode, comm, info, &ncid); Error(rcd);
  rcd=read(ncid,1); Error(rcd);
  rcd=nc_close(ncid); Error(rcd);
  
  if (info != MPI_INFO_NULL) MPI_Info_free(&info);
  
 fn_exit:
  MPI_Finalize();
  return 0;
} /* !main() */

