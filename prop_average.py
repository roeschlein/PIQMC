import glob
import numpy as np

def read_valid_propagator(filepath):
    with open(filepath, 'r') as f:
        header = f.readline().strip()
        if not header.endswith("n_configs = 10"):
            print(f"Skipping {filepath}: n_configs != 10")
            return None
        data = []
        for line in f:
            tokens = line.strip().split()
            if len(tokens) != 4:
                continue
            try:
                tau_idx = int(tokens[0])
                orb_idx = int(tokens[1])
                real = float(tokens[2])
                imag = float(tokens[3])
                if np.isnan(real) or np.isnan(imag):
                    print(f"Skipping {filepath}: contains NaNs")
                    return None
                data.append([tau_idx, orb_idx, real, imag])
            except ValueError:
                print(f"Skipping {filepath}: parsing error")
                return None
        return np.array(data)

# Gather all valid files
all_files = sorted(glob.glob("propagator_*.dat"))
valid_data = []

for filename in all_files:
    data = read_valid_propagator(filename)
    if data is not None:
        valid_data.append(data)

if not valid_data:
    print("No valid propagator files found.")
    exit()

# Ensure all valid arrays are the same shape
shapes = [d.shape for d in valid_data]
if len(set(shapes)) != 1:
    raise ValueError("Valid propagator files do not have consistent shapes.")

# Average the data
avg_data = sum(valid_data) / len(valid_data)

# Write to final output
with open("propagator_final.dat", "w") as f:
    f.write(f"# Averaged over {len(valid_data)} files, n_configs = 10\n")
    for row in avg_data:
        tau, orb, real, imag = int(row[0]), int(row[1]), row[2], row[3]
        f.write(f"{tau} {orb} {real:.10e} {imag:.10e}\n")

print(f"Done. Wrote averaged data to propagator_final.dat using {len(valid_data)} valid files.")

