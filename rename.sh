#!/bin/bash

# Get sorted list of files based on the number in their names
files=( $(ls A0_*.dat | sed -E 's/A0_([0-9]+)\.dat/\1 \0/' | sort -n | awk '{print $2}') )

# Rename them sequentially to A0_1.dat, A0_2.dat, ...
count=1
for file in "${files[@]}"; do
    new_name="A0_${count}.dat"
    if [[ "$file" != "$new_name" ]]; then
        mv "$file" "$new_name"
    fi
    ((count++))
done

