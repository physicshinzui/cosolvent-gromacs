#!/bin/bash 

## For Ne
gmx editconf -f noble/ne.pdb -o box.gro -box 1.0 1.0 1.0
gmx genconf -f box.gro -nbox 10 -o noble/nebox.gro

## For Ar
gmx editconf -f noble/ar.pdb -o box.gro -box 1.0 1.0 1.0 
gmx genconf -f box.gro -nbox 10 -o noble/arbox.gro

## For Kr
gmx editconf -f noble/kr.pdb -o box.gro -box 1.0 1.0 1.0
gmx genconf -f box.gro -nbox 10 -o noble/krbox.gro 

## For Xe
gmx editconf -f noble/xe.pdb -o box.gro -box 1.0 1.0 1.0
gmx genconf -f box.gro -nbox 10 -o noble/xebox.gro
# -rot: If you want to make a box consisting of organic molecules, you should add this option, which enable the generation of randomly orient conformations.

## For Benzene
#gmx editconf -f ben.pdb -o box.gro -box 1.0 1.0 1.0
#gmx genconf -f box.gro -nbox 10 -o benbox.gro -rot -seed 1122
# -rot: If you want to make a box consisting of organic molecules, you should add this option, which enable the generation of randomly orient conformations.

