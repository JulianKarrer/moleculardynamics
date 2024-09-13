# Molecular Dynamics in C++

Repository for coursework concerning the summer 2024 course [Molecular Dynamics in C++](https://pastewka.github.io/MolecularDynamics/) at the University of Freiburg, forked from [the MD Meson Template](https://github.com/imtek-simulation/meson-skeleton/).


- `/src` and `/tests` contain the library source code and tests in C++ using Eigen and gtest
- `/rust` contains a reimplementation of the program in Rust for comparison of runtimes.
- `/milestones` contain the source code for the binaries used to generate the data in the report
- `/analysis/analysis.py` is a standalone script that produces all figures used in the report from the `csv` files in the same folder
- `/initial_configurations` contain some `xyz` files that should be copied over to `/builddir` when setting up the project to ensure all milestones run smoothly (`lattice.xyz, cluster923.xyz` etc.)
- `/mackay_gold_clusters` contains more `xyz` files that were generated using `mackay_generate.sh` in the root directory


# Reproducing figures
Just run `/analysis/analysis.py`.

The corresponding `csv` files were copied over from `builddir` and can themselves be reproduced by running the code as described below for each imlestone.

# Running the Code
- Setup: ```
  cd <your repository>
  meson setup builddir --buildtype=release
  cd builddir
  meson compile
  meson test
  ```
- Milestone 1: Run  `builddir/milestones/01/milestone01` to check if the build system works
- Milestone 2: Run  `builddir/milestones/02/milestone02` for a test of energy conservation using Velocity Verlet 
- Milestone 3: Run `builddir/milestones/03/milestone03` and copy `lj54.xyz` to `builddir` to run a Lennard-Jones direct summation simulation, with results in `builddir/traj.xyz`
- Milestone 4 & 5: `builddir/milestones/05/milestone05` uses the parameters at the top of `milestones/05/main.cpp` to first equilibrate a lattice to a target temperature, writing results to `builddir/lattice.xyz`, then creates a `builddir/thermostat_temps.csv` file for plotting how temperature evolves over time for different relaxation times tau (not included in the report)
- Milestone 6: `builddir/milestones/06/milestone06` uses the simulation parameters at the top of `milestones/06/main.cpp` to produce timings of simulation runs with more and more atoms using a direct summation and neighbour list approach, saving results to `builddir/runtimes.csv`
- Milestone 7: `builddir/milestones/07/milestone07` requires the `mackay_gold_clusters` and copies of `cluster_923.xyz, 923_500K.xyz` copied to `builddir`, can be run using the parameters:
  -  `heat` to produce `923_500K.xyz` from `cluster_923.xyz`
  -  `timestep` to produce `builddir/dt_eam.csv`, finding an optimal time step size for the Ducastelle EAM potential with Cleri & Rosato parameterization
  - `timestepljds` to similarly produce `builddir/dt_ljds.csv` for finding an optimal time step size for the Lennard-Jones potential with direct summation

  after any of these optional functionalities, `gold_t_e.csv` is produced, using the parameters in `milestones/07/main.cpp` lines `188` onwards to heat up gold clusters of varying sizes slowly and record total and potential energies, temperature etc.
- Milestone 8: `builddir/milestones/08/milestone08` requires `cluster_923.xyz` and instead of being invoked manually can be used with the script `ms8_mpi_e_const.sh` to produce `hamiltonian_par.csv`, which demonstrates energy conservation for MPI parallelization
- Milestone 9: `builddir/milestones/09/milestone09` uses the parameters at the top of `milestones/09/main.cpp` to simulate a gold nanowire being pulled. Especially the strain rate, temperature and the filename at the top were adjusted for the report, producing correspondingly named `csv` files as a result. The commented out values of `PLOT_EVERY, STRETCH_BY_PERCENT, SCALE_AFTER_N` were used along with `long_whisker.xyz` to produce the MPI scaling tests in the number of processes in the report, performed using `ms9_scaling_test_mpi.sh` 

# Video of Gold Nanowire being pulled:
<a href="https://www.youtube.com/watch?v=aHFvy7gYslU" target="_blank"><img src="https://raw.githubusercontent.com/JulianKarrer/moleculardynamics/main/thumbnail.jpg"/></a>

