#!/bin/bash

shopt -s nullglob  # Prevent errors when no files match
files=(diag_*.f90)       # Store matching files in an array
njobs=${#files[@]} # Get the count of matched files
echo "${njobs} bosonic matsubara frequencies found"

a=1
b=${njobs}

while [ $a -le $b ]
do
	file_name="diag_${a}.f90"
	pbs_file="diag_${a}.pbs"

	printf "#!/bin/bash
#PBS -q default
#PBS -N diag_${a}
#PBS -j oe
#PBS -l nodes=1:ppn=16
#PBS -l mem=10gb
#PBS -l walltime=08:00:00
#PBS -W group_list=x-ccast-prj-krjevski

cd \$PBS_O_WORKDIR
module load openmpi
module load intel-parallel-studio
module load lapack/gcc/64
module load openblas
ifort diag_${a}.f90 -mkl -O3 -o diag_${a}
export OMP_NUM_THREADS=16        # For MKL's multithreading
export MKL_NUM_THREADS=16        # Recommended to set both
export MKL_DYNAMIC=false
./diag_${a}
" > ${pbs_file}

	qsub ${pbs_file}
	((a=a+1))
done

