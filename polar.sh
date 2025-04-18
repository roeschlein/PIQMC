#! /bin/bash

module load anaconda3 
python3 matsubara.py
gfortran M.f90 -o x
./x
gfortran Gm_ij.f90 -o x
./x
./script_Wm_xy.sh
./submit_Wm_xy.sh
