from Bio.PDB import PDBParser
import sys
import os
import subprocess

"""
Input: 
   structure file including hydrogens.
   NOTE: this file should be an inital structure of minimisation or md simulation.

Output: 
   index.ndx : chain indices are added.
"""

pdb = sys.argv[1]
parser = PDBParser()
structure = parser.get_structure("name", pdb)

# Make atom selections for each chain 
atom_ranges = []
counter = 0
nchains = 1
for model in structure:
    for i, chain in enumerate(model):
        first_atom_id = counter + 1
#        print("#Chain: ", chain)
        for atom in chain.get_atoms():
#            print(counter)
            counter += 1
        atom_ranges.append(f"a {first_atom_id}-{counter}") #for gromacs make_ndx command
        nchains += 1

print("atom ranges for each chain", atom_ranges)

# generate an index.ndx with the atom selections for each chain
selections = " & \"Backbone\"\n".join(atom_ranges) + " & \"Backbone\"\n"
print(selections)
command = f"""
gmx make_ndx -f {pdb} <<EOF
{selections}
q
EOF
"""
os.system(command)

# Rename the selections in index.ndx
with open("index.ndx", 'r') as fin, open("rename_index.ndx", 'w') as fout:
    content = fin.read()
    for ichain, r in enumerate(atom_ranges):
        content = content.replace(r.replace(' ', '_'), f"chain{ichain}")
    fout.write(content)
os.system("mv rename_index.ndx index.ndx")
