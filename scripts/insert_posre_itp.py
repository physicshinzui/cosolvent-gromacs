
import sys
import pdbfixer 
import shutil

def count_nchains(pdbfilename):
    return len(list(pdbfixer.PDBFixer(filename=pdbfilename).topology.chains()))
    
top = sys.argv[1]
pdb = sys.argv[2]

# Insert ifdef lines before the 2nd moleculetype line.
moleculetype_target_counter = count_nchains(pdb) + 1
with open(top, 'r') as infile, open('tmp', 'w') as outfile:
    count = 0
    for line in infile:
        if '[ moleculetype' in line:
            count += 1
            if count == moleculetype_target_counter:
                outfile.write('#ifdef POSRES\n#include "posre.itp"\n#endif\n\n')
        outfile.write(line)
shutil.move('tmp', top)
