import subprocess
import re
import sys
import os
"""
    Args:
        input gro
        concentration

    Return:
        number of noble gas atoms
"""

input_gro = sys.argv[1]
concentration_molar = float(sys.argv[2])

# Run gmx command and get output
output = subprocess.check_output(f"gmx editconf -f {input_gro} -o tmp.gro", shell=True)

# Extract volume from output
volume = float(re.search(r"Volume: (\d+\.\d+)", output.decode()).group(1))
#print("Volume nm^3 = ", volume)

# Calculate number of atoms for 0.1 M concentration (in nm^-3)
avogadro_number = 6.02214076e23 # per mol
volume_litres = volume * 1e-24  # convert nm^3 to L

number_of_atoms = int(concentration_molar * avogadro_number * volume_litres)

#print("#atoms to be inserted: ", number_of_atoms)
print(number_of_atoms)
os.remove("tmp.gro")
