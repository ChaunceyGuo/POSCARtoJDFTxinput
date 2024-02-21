#A python3 script for conversion of POSCAR file from VASP to input.in for JDFTx computations. 
#Writen by ChatGPT and Xiangyu Guo
#---------------------------------------------------------
##generate common.in file for jdftx computations
with open('common.in', 'w') as f:
    f.write('#---Pseudopotentials ---\n')
    f.write('ion-species GBRV/$ID_pbe.uspp\n')
    f.write('elec-cutoff 20 #Ecutsforpsiandrho\n')
    f.write('#---Geometry---\n')
    f.write('coulomb-interaction Slab 001 #Make z nonperiodic\n')
    f.write('coulomb-truncation-embed 0.5 0.5 0.5 #Specify center\n')
    f.write('cache-projectors no\n')
    f.write('#coordination\n')
    f.write('lattice Triclinic XXX\n')
    f.write('latt-scale 1.8897260 1.8897260 1.8897260\n')
    f.write('coords-type Lattice\n')
    f.write('#---Electronic ---\n')
    f.write('kpoint-folding 3 3 1 #Gamma-centeredk-mesh\n')
    f.write('elec-smearing Gauss 0.0073499\n')
    f.write('coulomb-truncation-ion-margin 1\n')
    f.write('#target-mu -0.2447504 #U=2.0V vs. SHE #Fixechempotential\n')
    f.write('#target-mu -0.1712518 #U=2.0V vs. SHE #Fixechempotential\n')
    f.write('#For the CANDLE solvation model, mu_SHE = -4.66 eV\n')
    f.write('#muSHE = 0.1712518 U=0 V vs. SHE\n')
    f.write('core-overlap-check None\n')
    f.write('#---Fluid---\n')
    f.write('#elec-initial-charge 0.5\n')
    f.write('symmetry-threshold 0.0005\n')
    f.write('spintype z-spin\n')
    f.write('elec-initial-magnetization 10 no\n')
    f.write('fluid LinearPCM #Classof solvation model\n')
    f.write('pcm-variant CANDLE #Specificmodelwithinclass\n')
    f.write('fluid-solvent H2O #Aqueous electrolyte\n')
    f.write('fluid-cation Na+ 1. #1 mol/L Na+ cation\n')
    f.write('fluid-anion F- 1. #1 mol/L F- anion\n')
    f.write('#---opt---\n')
    f.write('electronic-minimize nIterations 500\n')
    f.write('ionic-minimize  nIterations 500\n')
    f.write('#---Outputs ---\n')
    f.write('dump-name NbTe2.$VAR      #This will overwrite outputs from successive runs\n')
    f.write('initial-state NbTe2.$VAR  #This will initialize from the preceding calculation\n')
    f.write('dump End State\n')
#---------------------------------------------------------
#read lattice informtion from POSCAR file
import numpy as np

# Load the POSCAR file
with open('POSCAR', 'r') as f:
    lines = f.readlines()

# Extract the lattice vectors information
lat_vecs = np.array([list(map(float, line.split())) for line in lines[2:5]])

# Calculate the length and angles of the vectors
a, b, c = np.sqrt(np.sum(lat_vecs**2, axis=1))
alpha = np.arccos(np.dot(lat_vecs[1], lat_vecs[2]) / (b*c)) * 180 / np.pi
beta = np.arccos(np.dot(lat_vecs[0], lat_vecs[2]) / (a*c)) * 180 / np.pi
gamma = np.arccos(np.dot(lat_vecs[0], lat_vecs[1]) / (a*b)) * 180 / np.pi

# Format the XXX string with the extracted parameters
xxx_str = "{:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f}".format(a, b, c, alpha, beta, gamma)

# Replace the XXX in the script
with open('common.in', 'r') as f:
    input_str = f.read()

input_str = input_str.replace("lattice Triclinic XXX", "lattice Triclinic " + xxx_str)

with open('common.in', 'w') as f:
    f.write(input_str)

#---------------------------------------------------------
import numpy as np
import re

def parse_contcar(POSCAR):
    with open(POSCAR, 'r') as f:
        lines = f.readlines()

    # Get the chemical symbols from the first line of the file
    symbols = re.findall('[A-Za-z]+', lines[5])

    # Get the number of each type of atom from the second line of the file
    num_atoms = list(map(int, lines[6].split()))

    # Get the coordinates of each atom from the rest of the file, considering Selective dynamics
    coords = []
    for i in range(8, len(lines)):
        if lines[i].strip() == '' or lines[i].strip().startswith('Direct'):
            continue
        parts = lines[i].split()
        coords.append(list(map(float, parts[:3])))

    # Generate the output file
    with open('out.data', 'w') as f:
        for i, sym in enumerate(symbols):
            for j in range(num_atoms[i]):
                ion = "ion " + sym 
                direct_coords = " ".join([str(x) for x in coords.pop(0)])
                f.write(ion + " " + direct_coords + " " + "0" + "\n")

if __name__ == '__main__':
    parse_contcar('POSCAR')

#---------------------------------------------------------
# Read the data from out.data
with open('out.data', 'r') as f:
    out_data = f.read()

# Open the input script
with open('common.in', 'r') as f:
    input_lines = f.readlines()

# Find the line numbers for "coords-type Lattice" and "#---Electronic ---"
start_line = None
end_line = None
for i, line in enumerate(input_lines):
    if "coords-type Lattice" in line:
        start_line = i
    if "#---Electronic ---" in line:
        end_line = i
        break

# Insert the out.data contents between the start and end lines
if start_line is not None and end_line is not None:
    output_lines = input_lines[:start_line+1] + [out_data] + input_lines[end_line:]
else:
    output_lines = input_lines

# Write the modified script to a new file
with open('input.in', 'w') as f:
    f.writelines(output_lines)
