#!/bin/bash
set -eu
GMX=gmx
${GMX} pdb2gmx -f ../xebox.gro -o xe.gro -water none
${GMX} editconf -f xe.gro \
              -o xe.gro    \
              -d 0.0 \
              -bt cubic

${GMX} grompp -f em.mdp \
              -c xe.gro \
              -r xe.gro \
              -p topol.top  \
              -po mdout_em.mdp \
              -o em.tpr -maxwarn 1

${GMX} mdrun -deffnm em 

${GMX} grompp -f nvt_prod.mdp \
              -c em.gro \
              -r em.gro \
              -p topol.top  \
              -po mdout_nvt.mdp \
              -o nvt.tpr -maxwarn 1
${GMX} mdrun -deffnm nvt  -nsteps 100000

${GMX} grompp -f npt_prod.mdp \
              -c nvt.gro \
              -r nvt.gro \
              -p topol.top  \
              -po mdout_npt.mdp \
              -o npt.tpr -maxwarn 1
${GMX} mdrun -deffnm npt  -nsteps 100000
