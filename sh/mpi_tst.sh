#/bin/sh

# Purpose: Compile and mpirun simple MPI program to check for sane MPI environment

# Usage:
# cd ~/sh;chmod a+x mpi_tst.sh;./mpi_tst.sh;cd -
# ./mpi_tst.sh 3
# export MPI_PRC=4; ./mpi_tst.sh 

# Queue usage:
# cd ~/sh;bsub mpi_tst.sh
# cd ~/sh;qsub mpi_tst.sh

# Distribution:
# scp ~/sh/mpi_tst.sh greenplanet.ps.uci.edu:sh

# Source: Daniel Wang and Charlie Zender

# LSF: 70 8-way nodes available on blueice, 12 on GreenPlanet
#BSUB -a mpichp4             # fxm
#BSUB -J mpi_tst             # Job name
#BSUB -x                     # Exclusive use of node (not_shared)
#BSUB -n 8                   # Total tasks and threads (processors) needed
#BSUB -N                     # fxm
#BSUB -B                     # Send mail at dispatch and intiation times
#BSUB -R "span[ptile=8]"     # MPI tasks per node (blueice=8,greenplanet=8)
#BSUB -u zender@uci.edu      # Email notification
#BSUB -o lsf_%J.out          # Output filename
#BSUB -e lsf_%J.err          # Error filename
#BSUB -q zender              # Queue name
#BSUB -W 6:00                # 6 hour wallclock limit (required)

# Torque (e.g., greenplanet.ps.uci.edu, ipcc.ess.uci.edu, mpc.uci.edu)
## Job Name
#PBS -N mpi_tst
## Export all environment variables to job
#PBS -V
## Combine stdout with stderr
#PBS -j oe
## Send mail to this email address
#PBS -M zender@uci.edu
## Notify user via email at end or if aborted
#PBS -m ea
## PBS output file
#PBS -o pbs_mpi_tst.txt
## Queue name
#PBS -q zender
## Number of nodes:mpi_processes_per_node (see below)
#PBS -l nodes=1:ppn=8
##PBS -l nodes=compute-0-10+compute-0-11

# Grid Engine (aka SGE) (e.g., ipcc.ess.uci.edu, pbs.ess.uci.edu)
## Job Name
#$ -N mpi_tst
## Shell
#$ -S /bin/csh
## Exports all environment variables to the job
#$ -V
## Combine stdout with stderr
#$ -j yes
## Notify user via email at end or if aborted
#$ -m ea
## SGE output file
#$ -o sge_mpi_tst.txt
## Parallel environment (NB: This sets NSLOTS = total # of MPI processes on all nodes)
## Daniel created "parallel environments" (PEs) that allocate one node per MPI process: mpich-1, mpi-1, lam-1
#$ -pe openmpi 2

# Let environment override default
if [ -z "${MPI_PRC}" ]; then
    MPI_PRC='2' # [nbr] Number of MPI processes to spawn
fi # endif

# Let command line override environment
if [ -n "${1}" ]; then
    MPI_PRC=${1} # [nbr] Number of MPI processes to spawn
fi # !$1

TMP_TPL=foo.mpi.${USER}
TMP_BIN=${TMP_TPL}
TMP_SRC=${TMP_TPL}.c
TMP_OBJ=${TMP_TPL}.o

# Note: On pbs.ess.uci.edu, ${HOME} is only directory that:
# 1) ...is exported to compute nodes
# 2) ...is writable by us
# Write to current working directory

echo "$0: Testing whether current 'mpicc' and 'mpirun' work..."

# Cleanup residue from last invocation
rm -f ${TMP_BIN} ${TMP_SRC} ${TMP_OBJ}

# (Try to) Compile and run minimal MPI program
# mpicc does not work via stdin so make temporary source file first
#mpicc -x c -o $TMP_FILENAME -
cat > ${TMP_SRC} <<EOF
#include <stdio.h>
#include <mpi.h>
int main(int argc, char** argv){
#ifdef HOST_NAME_MAX
# define NCO_HOST_NAME_MAX HOST_NAME_MAX
#else
# define NCO_HOST_NAME_MAX 256
#endif
  char mpi_hst_nm[NCO_HOST_NAME_MAX];
  const int mpi_rnk_mgr=0;
  int mpi_hst_nm_lng;
  int mpi_prc_nbr=-1;
  int mpi_prc_rnk;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&mpi_prc_nbr);
  MPI_Comm_rank(MPI_COMM_WORLD,&mpi_prc_rnk);
  MPI_Get_processor_name(mpi_hst_nm,&mpi_hst_nm_lng);
  if(mpi_prc_rnk == mpi_rnk_mgr) (void)fprintf(stdout,"MPI: Manager reports total number of MPI processes is %d\n",mpi_prc_nbr);
  (void)fprintf(stdout,"MPI: MPI process ID mpi_prc_id = %d executing on host = %s\n",mpi_prc_rnk,mpi_hst_nm); 
  MPI_Finalize();
  return 0;
}
EOF
echo "$0: 'mpicc' resolves to this executable: `which mpicc`"
CMP_CMD="mpicc -o ${TMP_BIN} ${TMP_SRC}"
echo "$0: Compiling with '${CMP_CMD}' ..."
eval ${CMP_CMD}
if [ "$?" -ne '0' ]; then
    echo "$0: Error compiling ${TMP_SRC} into ${TMP_BIN}"
    exit 1
fi

# GreenPlanet MPI environments and results on interactive node:
# mpich_pgi_1.2.7p1 : 
# mvapich_pgi-1.1.0 : fails with Need to provide machinefile or a list of machines. Without hostfile option, hostnames must be specified on command line.
# mvapich2_pgi-1.2p1: fails with mpiexec_greenplanet.ps.uci.edu: cannot connect to local mpd (/tmp/mpd2.console_zender); (command works if mpd started by hand on local node, but how to get mpd to automatically start?)
# openmpi_pgi-1.3.2 : Works!

# GreenPlanet MPI environments and results via qsub to zender queue:
# mpich_pgi_1.2.7p1 : 
# mvapich_pgi-1.1.0 : Works!
# mvapich2_pgi-1.2p1: fails with mpiexec_compute-1-24.local: cannot connect to local mpd
# openmpi_pgi-1.3.2 : fails with foo.mpi.zender: Symbol `ompi_mpi_comm_world' has different size in shared object, consider re-linking

echo "$0: 'mpirun' resolves to this executable: `which mpirun`"
if [ -z ${PBS_NODEFILE} ]; then
    MPI_CMD="mpirun -np ${MPI_PRC} ${TMP_BIN}"
else
    MPI_CMD="mpirun -np ${MPI_PRC} -machinefile ${PBS_NODEFILE} ${TMP_BIN}"
fi # endif
#MPI_CMD="mpirun -np ${MPI_PRC} ${TMP_BIN}"
#MPI_CMD="mpirun -mca btl tcp,self -np ${MPI_PRC} ${TMP_BIN}"
#MPI_CMD="mpirun -mca btl ib -np ${MPI_PRC} ${TMP_BIN}"
echo "$0: Running simple MPI program with '${MPI_CMD}' ..."
eval ${MPI_CMD}
if [ "$?" -ne '0' ]; then
    echo "$0: Error mpirun'ing ${TMP_BIN}"
    exit 1
fi

# Cleanup
echo "$0: Cleaning up..."
rm -f ${TMP_BIN} ${TMP_SRC} ${TMP_OBJ}
