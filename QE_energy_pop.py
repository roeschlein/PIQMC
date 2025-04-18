import numpy as np
import os

line_index = []
on_index = []
energy = []
occupancy = []
eigenvalues = []

# Read the file and find indices
with open('./' + os.environ['outfile'], 'r') as readfile:
    for num, line in enumerate(readfile):
        if 'End of self-consistent calculation' in line:
            line_index.append(num)
        if 'occupation numbers' in line:
            on_index.append(num)

start_ind = int(max(line_index))
on_ind = int(max(on_index))

# Extract energy and occupation numbers
with open('./' + os.environ['outfile'], 'r') as orbitals:
    index = [s.strip() for s in orbitals.readlines()]
    
    energy = [float(j) for i in range(start_ind + 4, on_ind - 1) for j in index[i].split()]
    occupancy = [float(j) for i in range(on_ind + 1, 2 * on_ind - start_ind - 4) for j in index[i].split()]

# Create eigenvalues list
eigenvalues = [[i + 1, e, occupancy[i]] for i, e in enumerate(energy)]

# Write energy values
with open('./eigenvalues', 'w') as eigen:
    for i in energy:
        eigen.write(f"{i:.4f}\n")

# Write formatted energy-population data
with open('./energy_pop', 'w') as eigen:
    for i, e, o in eigenvalues:
        eigen.write(f"{i:8}{e:12.4f}{o:12.5f}\n")

print("Files 'eigenvalues' and 'energy_pop' have been successfully written.")

