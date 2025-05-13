#!/bin/bash

shopt -s nullglob  # Prevent errors when no files match
files=(prop_*.py)       # Store matching files in an array
njobs=${#files[@]} # Get the count of matched files
echo "Generating propagator scripts for <a(tau)a^\dag(0)>"

a=1
b=${njobs}

while [ $a -le $b ]
do
	file_name="prop_${a}.py"
	pbs_file="prop_${a}.pbs"

	printf "#!/bin/bash
#PBS -q default
#PBS -N prop_${a}
#PBS -j oe
#PBS -l nodes=1:ppn=1
#PBS -l mem=10gb
#PBS -l walltime=00:80:00
#PBS -W group_list=x-ccast-prj-krjevski

cd \$PBS_O_WORKDIR
module load anaconda3
python3 prop_${a}.py 
" > ${pbs_file}

	qsub ${pbs_file}
	((a=a+1))
done

