#!/bin/bash


# Read the number of processors per node
echo "Enter number of processors per node: "
read xstep
shopt -s nullglob  # Prevent errors when no files match
files=(W_*)       # Store matching files in an array
njobs=${#files[@]} # Get the count of matched files

a=1
b=${njobs}

while [ $a -le $b ]
do
	file_name="W_${a}.f90"
	pbs_file="W_${a}.pbs"

	printf "#!/bin/bash
#PBS -q default
#PBS -N W_${a}
#PBS -j oe
#PBS -l nodes=1:ppn=${xstep}
#PBS -l mem=50gb
#PBS -l walltime=00:80:00
#PBS -W group_list=x-ccast-prj-krjevski


cd \$PBS_O_WORKDIR
module load openmpi
module load intel-parallel-studio
ifort -qopenmp -integer-size 64 W_${a}.f90 -o W_${a}
chmod +x W_${a}
export OMP_NUM_THREADS=${xstep}


./W_${a}
" > ${pbs_file}
	
	qsub ${pbs_file}
	((a=a+1))
done


