#!/bin/bash 

ls *.gro *.tpr *# *.edr *.trr *.xtc *.top *.itp *.mdp *.log *.cpt
read -p "Are you really sure?[Enter]"
read -p "Are you really sure?[Enter]"
mkdir tmp ; cp -rp *.gro *.tpr *# *.edr *.trr *.xtc *.top *.itp *.mdp *.log *.cpt tmp

rm -v *.gro *.tpr *# *.edr *.trr *.xtc *.top *.itp *.mdp *.log *.cpt
