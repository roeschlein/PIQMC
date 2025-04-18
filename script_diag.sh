#!/bin/bash

read num_bose < matsubara_params

for k_boson in $(seq 1 ${num_bose})
do
	cat <<EOF> diag_${k_boson}.py
import numpy as np
import struct

def read_trim_xyz(filename="trim_xyz"):
    with open(filename, 'r') as f:
        x1, x2 = map(int, f.readline().split())
        y1, y2 = map(int, f.readline().split())
        z1, z2 = map(int, f.readline().split())
    return x1, x2, y1, y2, z1, z2

def read_sparse_bin(filename):
    with open(filename, 'rb') as f:
        # Read number of non-zero elements
        sparse_index = struct.unpack('q', f.read(8))[0]

        # Read values
        values = np.fromfile(f, dtype=np.float64, count=sparse_index)
        rows = np.fromfile(f, dtype=np.int64, count=sparse_index)
        cols = np.fromfile(f, dtype=np.int64, count=sparse_index)

    return sparse_index, values, rows, cols

def build_matrix(N, values, rows, cols, symmetric=True):
    M = np.zeros((N, N), dtype=np.float64)
    for i in range(len(values)):
        r = rows[i] - 1  # Fortran to Python index
        c = cols[i] - 1
        if 0 <= r < N and 0 <= c < N:
            M[r, c] = values[i]
            if symmetric and r != c:
                M[c, r] = values[i]
        else:
            print(f"Out of bounds: ({r+1}, {c+1})")  # Print 1-based index
    return M

def main():
    # Read domain bounds
    x1, x2, y1, y2, z1, z2 = read_trim_xyz()
    x_size = x2 - x1 + 1
    y_size = y2 - y1 + 1
    z_size = z2 - z1 + 1
    N = x_size * y_size * z_size

    # Read and build W
    _, W_values, W_rows, W_cols = read_sparse_bin("W_${k_boson}.bin")
    W = build_matrix(N, W_values, W_rows, W_cols, symmetric=True)

    # Read and build L
    _, L_values, L_rows, L_cols = read_sparse_bin("laplacian.bin")
    L = build_matrix(N, L_values, L_rows, L_cols, symmetric=False)

    # Build S_k
    S_k = -W - L

    print(f"S_k symmetry error: {np.linalg.norm(S_k - S_k.T) / np.linalg.norm(S_k)}")
    S_k = 0.5 * (S_k + S_k.T)

    # Diagonalize S_k
    eigvals, eigvecs = np.linalg.eigh(S_k)

    # Save eigenvalues
    np.savetxt("eigv_t${k_boson}", eigvals)

    # Save eigenvectors
    eigvecs.tofile("eigvec_t${k_boson}")

if __name__ == "__main__":
    main()
EOF
done

