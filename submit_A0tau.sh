#!/bin/bash

shopt -s nullglob  # Prevent errors when no files match
files=(A0_config_gen_*.f90)       # Store matching files in an array
njobs=${#files[@]} # Get the count of matched files
echo "Generating 1000 gauge field configurations"

a=1
b=${njobs}

while [ $a -le $b ]
do
	file_name="A0_config_gen_${a}.f90"
	pbs_file="A0_config_gen_${a}.pbs"

	printf "#!/bin/bash
#PBS -q default
#PBS -N A0_tau_${a}
#PBS -j oe
#PBS -l nodes=1:ppn=1
#PBS -l mem=10gb
#PBS -l walltime=60:00:00
#PBS -W group_list=x-ccast-prj-krjevski

cd \$PBS_O_WORKDIR
gfortran A0_config_gen_${a}.f90 -o A0_config_gen_${a}
./A0_config_gen_${a}
" > ${pbs_file}

	qsub ${pbs_file}
	((a=a+1))
done

