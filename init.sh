#!/bin/bash 

ls tleap.in system* *.gro *.tpr *# *.edr *.trr *.xtc *.top *.itp *.mdp *.log *.cpt
read -p "Are you really sure?[Enter]"
rm -v tleap.in system* *.gro *.tpr *# *.edr *.trr *.xtc *.top *.itp *.mdp *.log *.cpt
