#!/bin/bash

shopt -s nullglob  # Prevent errors when no files match
files=(diag_*)       # Store matching files in an array
njobs=${#files[@]} # Get the count of matched files
echo "${njobs} bosonic matsubara frequencies found"

a=1
b=${njobs}

while [ $a -le $b ]
do
	file_name="diag_${a}.py"
	pbs_file="diag_${a}.pbs"

	printf "#!/bin/bash
#PBS -q default
#PBS -N diag_${a}
#PBS -j oe
#PBS -l nodes=1:ppn=1
#PBS -l mem=10gb
#PBS -l walltime=00:40:00

#PBS -W group_list=x-ccast-prj-krjevski

####module load openmpi

cd  \$PBS_O_WORKDIR
module load anaconda3
python3 diag_${a}.py
" > ${pbs_file}

	qsub ${pbs_file}
	((a=a+1))
done

