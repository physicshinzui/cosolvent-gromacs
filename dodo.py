from glob import glob 
import os, sys 
import shutil
from pathlib import Path
import subprocess

path_to_pdbs =glob("/home/users/dva/noble_gas/data/DUD-E-mine/processed_pdbs/*.pdb")
print(path_to_pdbs)
finished_pdbids = glob("????")
print("Finished pdbs", finished_pdbids)
for i, path_to_pdb in enumerate(path_to_pdbs):
    if i > 2: sys.exit()
    pdbid = path_to_pdb.split('/')[-1].split('_')[0]
    if pdbid in finished_pdbids:
        print(f"Skipped: PDBID {pdbid} has been finished.")
        continue
    print(f"PDBID: {pdbid}")
    print(f"Path to: {path_to_pdb}")
    
    # Make a directory for a target protein
    subprocess.run(['bash', '../new.sh', pdbid])
    
    # Go to the dir and create a system
    shutil.copy(path_to_pdb, "./"+pdbid)
    os.chdir(pdbid)

    # Create a system 
    pdbfilename = path_to_pdb.split('/')[-1]
    subprocess.run(['bash', '01_create_system.sh', pdbfilename])
    
    # Submit runs
    subprocess.run(['bash', 'multi_runs.sh'])

    os.chdir('..')
