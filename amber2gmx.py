import parmed as pmd
import sys

usage=f"""
Usage 
    python {sys.argv[0]} [top] [crd]
"""
print(usage)
top, crd = sys.argv[1], sys.argv[2]
# convert AMBER topology to GROMACS, CHARMM formats
amber = pmd.load_file(top, crd)

# Save a GROMACS topology and GRO files
out = top.split(".")[0]
amber.save(f'{out}.top')
amber.save(f'{out}.gro')
