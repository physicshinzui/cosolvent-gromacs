
import sys
input_file = sys.argv[1]
output_file = sys.argv[2]

# Insert ifdef lines before the 2nd moleculetype line.
moleculetype_target_counter = 2
with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    count = 0
    for line in infile:
        if '[ moleculetype' in line:
            count += 1
            if count == moleculetype_target_counter:
                outfile.write('#ifdef POSRES\n#include "posre.itp"\n#endif\n\n')
        outfile.write(line)
