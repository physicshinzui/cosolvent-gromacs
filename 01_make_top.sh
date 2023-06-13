#!/bin/bash

#TODO: I wanna calculate Molar of noble gas automatially. 

set -eu
cat << EOF
Usage: 
    bash $0 [PDB file] 

Args:
    PDB file: structure.pdb

Returns: 
    - system.top
EOF

PDB=$1
FF=ff14SB 
tp=xe
nmol=10

# Create Amber topology and coordinate files
sed -e "s!#{INPUT}!${PDB}!g" -e "s!#{FF}!${FF}!g" templates/template_tleap.in > tleap.in
tleap -f tleap.in

# Convert Amber topology and coordinate files to Gromacs ones.
python amber2gmx.py system.prmtop system.inpcrd

# NOTE:  ==============================================================================
# Now we get a topology file and coordinate file of a protein, 
# but the protein system is solvated with cubic box, which may not be better than dodecahedron box in terms of computational cost.
# So, I want to use dodecahedron box rather than cubic one, but unfortunately, dodecahedron box is not implemented in tleap. 
# How do we use dodecahedron box?
#     -> Remove water and ion-related lines from the topology and coordinate files and solvate the protein again using Gromacs functionality! 
# ==============================================================================

# Remove non polymer
echo Protein | gmx trjconv -f system.gro -s system.gro -o protein.gro

# Delete the last three lines for Na+, Cl-, and WAT
# NOTE: the topological parameters still remains. The following lines just remove the number of Na+, Cl-, and WAT. 
sed -i.bak1 '$d' system.top
sed -i.bak2 '$d' system.top
sed -i.bak3 '$d' system.top

# Rename atom types for Gromacs naming convention (s/Amber/Gromacs/g)
sed -i.bak -e "s/WAT/SOL/g"  system.top
sed -i.bak '/SOL/s/H1/HW1/g' system.top
sed -i.bak '/SOL/s/H2/HW2/g' system.top
sed -i.bak '/SOL/s/ O / OW/g' system.top
# Note: Afterwards, Na+, Cl-, and WAT(=SOL) will be added by gmx commands

# Create a system with a box 
gmx editconf -f protein.gro -o newbox.gro -bt dodecahedron -d 1.0

# Solvate 
gmx solvate -cp newbox.gro -cs spc216.gro -o mol_solv.gro -p system.top

#==============insert noble gas atoms and modify a system topology file========
## Insert noble gas atoms with which SOL is replaced. 
gmx insert-molecules -f mol_solv.gro -ci ./box/noble/${tp}.pdb -nmol $nmol -replace SOL -o mol_solv.gro

## Insert #include "noble_LJTS.itp" before the first [ moleculetype ] in system.top
## beacuse "*.itp" has to be placed before any moleculetype section.
cp ../itp/noble_LJTS.itp noble_LJTS.itp
sed -i '0,/\[ moleculetype \]/s/\[ moleculetype \]/#include "noble_LJTS.itp"\n\n&/' system.top

## Place the new number of water molecules
new_nwat=`grep OW mol_solv.gro | wc -l`
sed -i "/^\[ molecules \]/,/^$/s/SOL\s\+[0-9]\+$/SOL          $new_nwat/" system.top

## Place the nubmer of noble gas atoms
echo "${tp^^}          $nmol" >> system.top 
#============================================

# Add ions; NOTE: Neutralisation is performed only.
gmx grompp -f templates/ions.mdp -c mol_solv.gro -p system.top -o ions.tpr #-maxwarn 1
echo SOL | gmx genion -s ions.tpr -o mol_solv_ions.gro -p system.top -pname Na+ -nname Cl- -neutral

### Create index file and posres.itp for protein chains
echo "Backbone" | gmx genrestr -f mol_solv_ions.gro -o posre.itp
python scripts/insert_posre_itp.py system.top tmp
mv tmp system.top

echo "Energy minimisation 1 ..."
gmx grompp -f ./templates/em1.mdp \
              -c mol_solv_ions.gro \
              -r mol_solv_ions.gro \
              -p system.top \
              -po mdout_em1.mdp \
              -o em1.tpr -maxwarn 1
gmx mdrun -deffnm em1

echo "Energy minimisation 2 ..."
gmx grompp -f ./templates/em2.mdp \
              -c em1.gro \
              -p system.top \
              -po mdout_em2.mdp \
              -o em2.tpr -maxwarn 1
gmx mdrun -deffnm em2

rm \#*
