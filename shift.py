import numpy as np

starting_index = int(input("Index of lowest KS eigenvalue in active window: "))
ending_index = int(input("Index of highest KS eigenvalue in active window: "))

# Load data (assuming columns: index, eigenvalue, occupation)
data = np.loadtxt("energy_pop")

eigenvalues = data[:, 1]  # Second column: eigenvalues
occupations = data[:, 2]  # Third column: occupations

# Find HOMO (highest occupied molecular orbital)
homo = max(eigenvalues[occupations > 0])

# Find LUMO (lowest unoccupied molecular orbital)
lumo = min(eigenvalues[occupations == 0])

# Compute mid-gap energy
mid_gap = (homo + lumo) / 2

# Shift eigenvalues
data[:, 1] -= mid_gap

# Save only the active window eps
np.savetxt("energy_pop_shifted.dat", data[starting_index-1:ending_index], fmt="%.6f")

