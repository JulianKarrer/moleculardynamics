#!/bin/bash

output_file="hamiltonian_par.csv"

echo COMPILING PROGRAM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cd builddir
meson compile

# reset csv file and write header
rm -f $output_file
echo "nb_processes,n,e_total,e_pot,e_kin" >> $output_file

for i in {1,2,4,8,16,27}
do
  echo RUNNING MPI SIMULATION WITH $i PROCESSES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Use the time command to measure execution time
  { time -p mpirun --oversubscribe -n $i ./milestones/08/milestone08 | grep -E '^[0-9]' >> $output_file; } 2> time_tmp.txt

  # Extract and print the timing information
  real_time=$(grep "real" time_tmp.txt | awk '{print $2}')
  user_time=$(grep "user" time_tmp.txt | awk '{print $2}')
  sys_time=$(grep "sys" time_tmp.txt | awk '{print $2}')

  echo "Timing for $i processes: real $real_time s, user $user_time s, sys $sys_time s"
done

#cleanup
rm -f time_tmp.txt

# execute analysis script
echo RUNNING ANALYSIS SCRIPT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cd ../analysis
python3 analysis.py