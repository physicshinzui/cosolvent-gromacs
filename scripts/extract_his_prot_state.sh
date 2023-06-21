#!/bin/bash 
pdb=$1
gmx pdb2gmx -f $1 -ff amber03 -water tip3p -o conf.gro &> pdb2gmx.log
grep "Will use" pdb2gmx.log | awk '{print $3 "    " $6}' > his_states.txt
cat his_states.txt | sed -e "s/HISD/HID/g" -e "s/HISE/HIE/g" -e "s/HISH/HIP/g" > tmp 
mv tmp his_states.txt
rm topol_*.itp posre_*.itp topol.top conf.gro
