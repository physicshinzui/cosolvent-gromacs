#!/bin/bash
#PBS -l select=1:ncpus=16:mpiprocs=2:ompthreads=2:ngpus=1
#PBS -l walltime=00:30:00

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
export GMXLIB=$HOME/noble_gas/top
#==============================

set -eu
id=xe
GMX="gmx"

#================Equilibriation step==============
echo "NVT equilibration runs are running..."
cat ./templates/nvt_eq.mdp | sed -e "s!#{RAND}!${RANDOM}!g" > nvt_eq_${id}.mdp
${GMX} grompp -f nvt_eq_${id}.mdp \
              -c em2.gro \
              -r em2.gro \
              -p topol.top \
              -po mdout_nvt_eq.mdp \
              -o nvt_eq_${id}.tpr
${GMX} mdrun -deffnm nvt_eq_${id} -ntomp $OMP_NUM_THREADS -ntmpi $NUM_MPI

echo "NPT equilibration runs are running..."
cp ./templates/npt_eq.mdp npt_eq_${id}.mdp
${GMX} grompp -f npt_eq_${id}.mdp \
              -c nvt_eq_${id}.gro \
              -r nvt_eq_${id}.gro \
              -p topol.top  \
              -po mdout_npt_eq.mdp \
              -o npt_eq_${id}.tpr -maxwarn 1 # if Gromacs says that Berendsen barostat does not guarantee canonical distribution, you should add this option. 
${GMX} mdrun -deffnm npt_eq_${id} -ntomp $OMP_NUM_THREADS -ntmpi $NUM_MPI
#==============================================

#===========Production run=====================
echo "NPT runs are running..."
cp ./templates/npt_prod.mdp npt_prod_${id}.mdp
$GMX grompp -f npt_prod_${id}.mdp  \
            -c npt_eq_${id}.gro    \
            -t npt_eq_${id}.cpt    \
            -p topol.top           \
            -po mdout_npt_prod.mdp \
            -o npt_prod_${id}.tpr
# - Starting coordinates can be read from trajectory with -t
#   - Only if this information is absent will the coordinates in the -c file be used.

$GMX mdrun -deffnm npt_prod_${id} -nsteps 500 -ntomp $OMP_NUM_THREADS -ntmpi $NUM_MPI
#==============================================
