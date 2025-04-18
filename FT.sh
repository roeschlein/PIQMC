#!/bin/bash

# Read grid and box sizes from delta_xyz
read box_size_x box_size_y box_size_z < delta_xyz
read N1 N2 N3 < <(sed -n '2p' delta_xyz)
read HO_index < <(sed -n '1p' active_window)
read starting_index < <(sed -n '2p' active_window)
read ending_index < <(sed -n '3p' active_window)

AW_size=$((ending_index - starting_index + 1))

echo "Enter the number of orbitals per job: "
read KS_per_job

echo "Unit cell dimensions were read from delta_xyz as ${box_size_x}, ${box_size_y}, ${box_size_z}"
echo "Chosen FFT real grid dimensions were chosen to be ${N1}, ${N2}, ${N3}"

echo "Active window parameters: AW_size=${AW_size}, HO_index=${HO_index}, starting_index=${starting_index}, ending_index=${ending_index}"


let starting_index=starting_index-1
final_KS=$ending_index
let ending_index=ending_index-2

active_window=(${ending_index}-${starting_index})
num_jobs=$((active_window/KS_per_job))
let num_jobs=num_jobs-1

for bash_index in $(seq 0 ${num_jobs})
do
	start_index=$((bash_index * ${KS_per_job} + starting_index + 1))
	end_index=$(((bash_index + 1) * ${KS_per_job} + starting_index))
    cat <<EOF > FT_${start_index}_${end_index}.py
import numpy as np
import finufft
import multiprocessing as mp
import os
import warnings

A_to_bohr = 1.8897261254578281


warnings.simplefilter("ignore", UserWarning)

# Define real-space grid (uniform)
N1, N2, N3 = ${N1}, ${N2}, ${N3}  # Grid size
box_size_x, box_size_y, box_size_z = ${box_size_x}, ${box_size_y}, ${box_size_z}  # Box sizes

# Function to load data in parallel
def load_data(filename):
    if os.path.exists(filename):
        return np.loadtxt(filename)
    else:
        raise FileNotFoundError(f"{filename} not found.")
        
def FT_evc(index_number):
    with mp.Pool(2) as pool:
        g_vectors, evc_data = pool.map(load_data, ["gkvectors", "evc" + str(index_number)])

    evc_complex = evc_data[:, 0] + 1j * evc_data[:, 1]
    g_vectors = np.ascontiguousarray(g_vectors)
    evc_complex = np.ascontiguousarray(evc_complex)

    #delta_xyz is given in A, we won't use bohr radius units after the FFT
    delta_x = box_size_x / N1
    delta_y = box_size_y / N2
    delta_z = box_size_z / N3

    g_vectors[:, 0] *= (box_size_x * A_to_bohr) / N1
    g_vectors[:, 1] *= (box_size_y * A_to_bohr) / N2
    g_vectors[:, 2] *= (box_size_z * A_to_bohr) / N3

    psi_r = finufft.nufft3d1(
        g_vectors[:, 0], g_vectors[:, 1], g_vectors[:, 2],
        evc_complex,
        (N1, N2, N3),
    )

    normalization_factor = np.sqrt(box_size_x * box_size_y * box_size_z)  
    psi_r /= normalization_factor
    psi_r.flatten().tofile("evcr" + str(index_number) + ".bin")
    inner_product_real = np.sum(np.conj(psi_r) * psi_r) * (delta_x * delta_y * delta_z)
    print(f"<phi_{index_number}|phi_{index_number}>= {inner_product_real}")


for i in range(int($bash_index*$KS_per_job+$starting_index), int($bash_index*$KS_per_job + $KS_per_job+$starting_index)):
    index = i + 1
    FT_evc(index)
EOF

final_label=$(((num_jobs+1) * KS_per_job + starting_index + 1))
let end_index=end_index+KS_per_job

cat <<EOF > FT_${final_label}_${final_KS}.py
import numpy as np
import finufft
import multiprocessing as mp
import os
import warnings

A_to_bohr = 1.8897261254578281

warnings.simplefilter("ignore", UserWarning)

# Define real-space grid (uniform)
N1, N2, N3 = ${N1}, ${N2}, ${N3}  # Grid size
box_size_x, box_size_y, box_size_z = ${box_size_x}, ${box_size_y}, ${box_size_z}  # Box sizes

def load_data(filename):
    if os.path.exists(filename):
        return np.loadtxt(filename)
    else:
        raise FileNotFoundError(f"{filename} not found.")

def FT_evc(index_number):
    with mp.Pool(2) as pool:
        g_vectors, evc_data = pool.map(load_data, ["gkvectors", "evc" + str(index_number)])

    evc_complex = evc_data[:, 0] + 1j * evc_data[:, 1]
    g_vectors = np.ascontiguousarray(g_vectors)
    evc_complex = np.ascontiguousarray(evc_complex)

    #delta_xyz is given in A, we won't use bohr radius units after the FFT
    delta_x = box_size_x / N1
    delta_y = box_size_y / N2
    delta_z = box_size_z / N3

    g_vectors[:, 0] *= (box_size_x * A_to_bohr) / N1
    g_vectors[:, 1] *= (box_size_y * A_to_bohr) / N2
    g_vectors[:, 2] *= (box_size_z * A_to_bohr) / N3

    psi_r = finufft.nufft3d1(
        g_vectors[:, 0], g_vectors[:, 1], g_vectors[:, 2],
        evc_complex,
        (N1, N2, N3),
    )

    normalization_factor = np.sqrt(box_size_x * box_size_y * box_size_z)
    psi_r /= normalization_factor
    psi_r.flatten().tofile("evcr" + str(index_number) + ".bin")
    inner_product_real = np.sum(np.conj(psi_r) * psi_r) * (delta_x * delta_y * delta_z)
    print(f"<phi_{index_number}|phi_{index_number}>= {inner_product_real}")
"""
grid_points = np.array(np.meshgrid(
    np.linspace(0, box_size_x, N1, endpoint=False),
    np.linspace(0, box_size_y, N2, endpoint=False),
    np.linspace(0, box_size_z, N3, endpoint=False), indexing='ij'
)).reshape(3, -1).T

np.savetxt("GRID", grid_points)
"""
for i in range(int($final_label-1), int($final_KS)):
    index = i + 1
    FT_evc(index)
EOF
done

# Read the number of processors per node
echo "Enter number of processors per node: "
read xstep
shopt -s nullglob  # Prevent errors when no files match
files=(FT_*)       # Store matching files in an array
njobs=${#files[@]} # Get the count of matched files

echo "There were $njobs jobs produced."
let njobs=njobs-2
let starting_index=starting_index+1
let ending_index=ending_index+2

let norb=KS_per_job


# Start the loop from 0 to bash_index_end, and create PBS files accordingly
a=0
b=${njobs}  # Set the end of the bash index loop

while [ $a -le $b ]
do
    # Compute the start and end indices for the Python script filenames
    ((xc1=a*${norb}+${starting_index}))
    ((xc2=(a+1)*${norb}+${starting_index}-1))

    # Construct the Python file name
    python_file="FT_${xc1}_${xc2}.py"

    # Generate the corresponding .pbs file for the Python script
    pbs_file="FT_${xc1}_${xc2}.pbs"

    printf "#!/bin/bash
#PBS -q default
#PBS -N FT_${xc1}_${xc2}
#PBS -j oe
#PBS -l nodes=1:ppn=${xstep}
#PBS -l mem=2gb
#PBS -l walltime=00:10:00
#PBS -W group_list=x-ccast-prj-krjevski

module load openmpi

cd \$PBS_O_WORKDIR

python3 ${python_file}
" > ${pbs_file}

    # Submit the job to the queue using qsub
    qsub ${pbs_file}

    # Increment the index
    ((a=a+1))  

done

((xc1_last=(${njobs}+1)*${norb}+${starting_index}))
last_python_file="FT_${xc1_last}_${ending_index}.py"
last_pbs_file="FT_${xc1_last}_${ending_index}.pbs"
    printf "#!/bin/bash
#PBS -q default
#PBS -N FT_${xc1_last}_${ending_index}
#PBS -j oe
#PBS -l nodes=1:ppn=${xstep}
#PBS -l mem=2gb
#PBS -l walltime=00:10:00
#PBS -W group_list=x-ccast-prj-krjevski

module load openmpi

cd \$PBS_O_WORKDIR

python3 ${last_python_file}
" > ${last_pbs_file}

    # Submit the job to the queue using qsub
    qsub ${last_pbs_file}

