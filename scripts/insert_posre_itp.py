
import sys
import pdbfixer 
import shutil

def count_nchains(pdbfilename):
    return len(list(pdbfixer.PDBFixer(filename=pdbfilename).topology.chains()))
    
top = sys.argv[1]
pdb = sys.argv[2]
# PDB is used to count the number of chains

# Insert ifdef lines at the end of each moleculetype line.
nchains = count_nchains(pdb)
print(nchains)
with open(top, 'r') as infile, open('tmp', 'w') as outfile:
    ichain = 0
    for line in infile:
        if '[ moleculetype' in line:
            if ichain == 0:
                ichain += 1
            elif ichain <= nchains:
                outfile.write(f"#ifdef POSRES\n#include \"posre{ichain - 1}.itp\"\n#endif\n\n")
                ichain += 1
        outfile.write(line)
shutil.move('tmp', top)
