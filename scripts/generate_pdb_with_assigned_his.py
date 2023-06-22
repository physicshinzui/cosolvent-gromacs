#!/usr/bin/env python3
import sys
import os
import re

def read_his_state(filename):
    fin = open(filename, "r")

    dict_resi_w_hisname = {}
    for line in fin:
        if line[0] == "#": 
            continue
        hisname_amber = line.split()[2]
        resi = int(line.split()[3])
        if resi in dict_resi_w_hisname.keys():
            sys.exit(f"{resi} has already been included in keys.{dict_resi_w_hisname.keys()}")
        dict_resi_w_hisname[resi] = hisname_amber
    return dict_resi_w_hisname

def write_pdb_with_assigned_his(pdb, d):
    """
    Write a PDB by replacing 'HIS' with Amber's notation for histidine residue.

    Args: 
        pdb: PDB file
        d: a dictionary having keys are residue numbers and values are histidine's 3-letter name of Amber, e.g., HID, HIE, HIP

    Return: 
        None
    """
    fin = open(pdb, "r")
    fout = open("assigned_his.pdb","w")
    for line in fin:
        if line[0:6] != "ATOM  " and line[0:6] != "HETATM": 
            continue
        resi = int(line[22:27].strip())

        if resi in d.keys():
            fout.write(line.replace("HIS",d[resi]))
        else:
            fout.write(line)
    fin.close()
    fout.close()
    return None

def grep_his_assignment(pdb2gmxlog):
    fout = open("assigned_his_states.dat", "w")
    fout.write("#Gromacs -> Amber\n")
    with open(pdb2gmxlog, 'r') as fin:
        for line in fin:
            m = re.search("Will use HIS", line)

            if m: 
                his_type = line.split()[2]
                his_resi = line.split()[5]
                if his_type == "HISE": 
                    his_type_amber = "HIE"
                elif his_type == "HISD": 
                    his_type_amber = "HID"
                elif his_type == "HISH": 
                    his_type_amber = "HIP"
                else:
                    sys.exit(f"His type {his_type} seems invalid!")

                fout.write(f"{his_type} -> {his_type_amber} {his_resi}\n")

    return None

def main():
    help=f"""
    $0 is a script to get histidine states estimated by gmx pdb2gmx command. 
    , and give a list of pairs between Gromacs and Amber histidine 3-letters

    - ε-histidine is the default "HIE"
    - If the residue is called "HID" in the PDB file, the resulting residue for AMBER will become δ-histidine, 
    - "HIP" will yield the protonated form.

    # Notation list. Gro:Amber
    HISE:HIE
    HISD:HID
    HISH:HIP
    """
    print(help)
    PDB=sys.argv[1]
    GMX='gmx'
    os.system(f"echo \"1\n1\" | {GMX} pdb2gmx -f {PDB} -o pdb2gmx.pdb &> pdb2gmx.log")
    os.system(f"rm topol* posre* pdb2gmx.pdb")
    grep_his_assignment('pdb2gmx.log')
    d = read_his_state("assigned_his_states.dat")
    
    print("---Renamed Histidines with Amber-notation---")
    for key, val in d.items():
        print(key, val)

    write_pdb_with_assigned_his(PDB, d)
    
    print("\nOUTPUT: assigned_his.pdb")

if __name__ == "__main__":
    main()
