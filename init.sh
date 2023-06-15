#!/bin/bash 

ls p4a* tleap.in system* *.gro *.tpr *# *.edr *.trr *.xtc *.top *.itp *.mdp *.log *.cpt
read -p "Are you really sure?[Enter]"
rm -v p4a* tleap.in system* *.gro *.tpr *# *.edr *.trr *.xtc *.top *.itp *.mdp *.log *.cpt
