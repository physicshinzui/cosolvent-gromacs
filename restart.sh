#!/bin/bash
#PBS -l select=1:ncpus=16:mpiprocs=2:ompthreads=2:ngpus=1
#PBS -l walltime=5:00:00

#===========Preamble=============
if [ ! -z "${PBS_O_WORKDIR}" ]; then
  cd "${PBS_O_WORKDIR}"
fi

. /apl/hpc-x/2.11/hpcx-rebuild-gcc11.sh
hpcx_load
. /apl/gromacs/2022.4-CUDA/bin/GMXRC.bash
export LD_LIBRARY_PATH="/apl/cuda/12.0/lib64:${LD_LIBRARY_PATH}:/apl/pbs/22.05.11/lib"
export OMPI_MCA_btl=^openib

OMP_NUM_THREADS=16
NUM_MPI=1
module -s load gromacs/2022.4-CUDA
#==============================

set -eu
tpr=npt_prod_$id.tpr
cpt=npt_prod_$id.cpt
# jsub -v id=1 restart.sh

#==== If the simulation described by tpr file has completed and should be extended, then
gmx convert-tpr -s $tpr -until 50000 -o ${tpr}
gmx mdrun -deffnm npt_prod_$id -s $tpr -cpi $cpt

#cf. https://manual.gromacs.org/current/user-guide/managing-simulations.html
