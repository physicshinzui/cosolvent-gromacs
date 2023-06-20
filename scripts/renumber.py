import sys

def update_posre(filename, f_renumbered):
    with open(filename, 'r') as fin, open(f_renumbered, 'w') as fout:
        first_atom_id = None
        for i, line in enumerate(fin):
            skipped_condition = line.strip().startswith((';', '[')) or not line.strip()
            print(i, "is skipped?:", skipped_condition)
            if skipped_condition:
                print("line: ", line.strip())
                fout.write(line)
                continue

            # For DBG
            #if i == 10: 
            ##    print(line)
            #    break

            if not first_atom_id:
                first_atom_id = line.strip().split()[0]
                ref_val = int(first_atom_id) - 1
                print("First atom id was obtained!: ", first_atom_id)
                print("    Reference value was set to ", ref_val)
            
            sp_line = line.strip().split()
            sp_line[0] = str(int(sp_line[0]) - ref_val)
            print("DBG: before ", line.strip().split())
            print("DBG: after  ", sp_line)
            new_line = '       '.join(sp_line)
            fout.write(new_line+'\n')
             
filename = sys.argv[1] # posre.itp
f_renumbered = sys.argv[2] # renamed_posre.itp
#filename = 'posre1.itp'
#filename = 'posre0.itp'
update_posre(filename, f_renumbered)
