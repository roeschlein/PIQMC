#!/bin/bash


for A0_int in $(seq 1 100)
do
	cat <<EOF> prop_${A0_int}.py
import numpy as np
import scipy.linalg

def read_delta_xyz():
    with open('../../delta_xyz', 'r') as f:
        x, y, z = map(float, f.readline().split())
        N1, N2, N3 = map(int, f.readline().split())
    return x, y, z, N1, N2, N3

def read_trim_xyz():
    with open('../../trim_xyz', 'r') as f:
        x1, x2 = map(int, f.readline().split())
        y1, y2 = map(int, f.readline().split())
        z1, z2 = map(int, f.readline().split())
    return x1, x2, y1, y2, z1, z2

def read_A0_tau(index, Nt, phi_x, phi_y, phi_z):
    phi_size = phi_x * phi_y * phi_z
    try:
        A0_data = np.loadtxt(f'../A0/A0_{index+1}.dat')
    except (OSError, FileNotFoundError):
        print(f"Warning: A0_{index+1}.dat not found. Skipping.")
        return None

    if A0_data.ndim != 1:
        print(f"Warning: A0_{index+1}.dat is malformed. Expected a single column. Skipping.")
        return None

    if np.isnan(A0_data).any():
        print(f"Warning: A0_{index+1}.dat contains NaNs. Skipping.")
        return None

    if len(A0_data) != Nt * phi_size:
        print(f"Warning: A0_{index+1}.dat has unexpected size. Expected {Nt * phi_size}, got {len(A0_data)}. Skipping.")
        return None

    A0_flat = A0_data.reshape((Nt, phi_size))
    A0_reshaped = np.zeros((Nt, phi_x, phi_y, phi_z), dtype=np.float64)
    for x in range(phi_x):
        for y in range(phi_y):
            for z in range(phi_z):
                i = x * phi_y * phi_z + y * phi_z + z
                A0_reshaped[:, x, y, z] = A0_flat[:, i]
    return A0_reshaped

def read_active_window():
    with open('../../active_window', 'r') as f:
        f.readline()
        KS_low = int(f.readline().split()[0])
        KS_high = int(f.readline().split()[0])
    return KS_high - KS_low + 1

def load_matsubara_params():
    with open('../../matsubara_params', 'r') as f:
        f.readline()
        f.readline()
        T = float(f.readline().split()[0])
    return T

def read_phi_bin(N1, N2, N3, active_window_size):
    total_elements = active_window_size * N1 * N2 * N3
    phi_flat = np.fromfile('../../phi.bin', dtype=np.complex128, count=total_elements)
    return phi_flat.reshape((active_window_size, N1, N2, N3), order='F')

def construct_M(Aij, delta, active_window, eps):
    Nt, dim, _ = Aij.shape
    exp_minus_iAij = np.array([scipy.linalg.expm(-1j * delta * Aij[t]) for t in range(Nt)])
    M = np.zeros((Nt * active_window, Nt * active_window), dtype=np.complex128)
    for i in range(active_window):
        for j in range(active_window):
            for tau1 in range(Nt):
                for tau2 in range(Nt):
                    idx_1 = tau1 * active_window + i
                    idx_2 = tau2 * active_window + j
                    if i == j:
                        if tau1 == tau2:
                            M[idx_1, idx_2] = delta * eps[i] - exp_minus_iAij[tau1, i, i]
                        elif tau1 + 1 == tau2:
                            M[idx_1, idx_2] = 1
                        elif tau1 == Nt - 1 and tau2 == 0:
                            M[idx_1, idx_2] = -1
                    elif tau1 == tau2:
                        M[idx_1, idx_2] = -exp_minus_iAij[tau1, i, j]
    return M

def compute_correlator(M, Nt, active_window):
    invM = np.linalg.inv(M)
    corr = np.zeros((Nt, active_window), dtype=np.complex128)
    for e in range(active_window):
        for tau in range(Nt):
            val = 0.0
            for t0 in range(Nt):
                t = (t0 + tau) % Nt
                sign = 1 if (t0 + tau) < Nt else -1
                col = t * active_window + e
                row = t0 * active_window + e
                val += sign * invM[row, col]
            corr[tau, e] = val / Nt
    return corr

def running_average(avg, N_avg, new_val):
    N_avg += 1
    a = 1 / N_avg
    b = 1 - a
    avg = a * new_val + b * avg
    return avg, N_avg

if __name__ == "__main__":
    T = load_matsubara_params()
    x1, x2, y1, y2, z1, z2 = read_trim_xyz()
    phi_x, phi_y, phi_z = x2 - x1 + 1, y2 - y1 + 1, z2 - z1 + 1
    beta = 1.0 / T
    Nt = 120
    delta = beta / Nt
    active_window = read_active_window()
    x, y, z, N1, N2, N3 = read_delta_xyz()
    phi = read_phi_bin(N1, N2, N3, active_window)
    phi_trim = phi[:, x1-1:x2, y1-1:y2, z1-1:z2]
    eps = np.loadtxt("../../energy_pop_shifted.dat")[:, 1]

    N_avg = 0
    Monte_Carlo_avg = np.zeros((Nt, active_window), dtype=np.complex128)

    for index in range((${A0_int}-1)*10,${A0_int}*10):
        print(fr"Attempting A0_{index+1}...")
        A0 = read_A0_tau(index, Nt=Nt, phi_x=phi_x, phi_y=phi_y, phi_z=phi_z)
        if A0 is None:
            continue
        Aij = (x/N1 * y/N2 * z/N3) * np.einsum('ixyz,txyz,jxyz->tij', np.conj(phi_trim), A0, phi_trim)
        M = construct_M(Aij, delta, active_window, eps)
        corr = compute_correlator(M, Nt, active_window)
        Monte_Carlo_avg, N_avg = running_average(Monte_Carlo_avg, N_avg, corr)
        with open("propagator_${A0_int}.dat", 'w') as f_out:
            f_out.write(f"tau_index orbital_index real_part imag_part n_configs = {N_avg}\n")
            for tau in range(Nt):
                for e in range(active_window):
                    re = np.real(Monte_Carlo_avg[tau, e])
                    im = np.imag(Monte_Carlo_avg[tau, e])
                    f_out.write(f"{tau} {e} {re:.10e} {im:.10e}\n")
EOF
done

