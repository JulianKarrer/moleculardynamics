#include "atoms.h"
#include "mpi.h"
#include "verlet.h"
#include "xyz.h"
#include <domain.h>
#include <ducastelle.h>
#include <iomanip>
#include <iostream>
#include <neighbors.h>
#include <thermostat.h>

// settings
const double CUTOFF_RADIUS{10.};
const double DT{0.1};
const size_t TIMESTEPS{5000};
const size_t PLOT_EVERY{10};

/// @brief Specify one or more arguments to run a test on what timestep size
/// is appropriate for this potential. Otherwise, run the simulation and
/// collect data.
int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    // create csv file for measurement results
    std::ofstream hamiltonian_par("hamiltonian_par.csv");
    hamiltonian_par << "n,e_total,e_pot,e_kin" << std::endl;

    // initialize gold cluster from file
    auto [names, positions,
          velocities]{read_xyz_with_velocities("cluster_923.xyz")};
    Atoms atoms{Atoms(names, positions)};

    // create neighbour list and initialize forces
    NeighborList neighbour_list;
    neighbour_list.update(atoms, CUTOFF_RADIUS);
    (void)ducastelle(atoms, neighbour_list);
    double current_temp{0.};

    // initialize a domain for periodic boundaries and decomposition for
    // parallelization
    Vec3_t volume{atoms.positions.rowwise().maxCoeff() -
                  atoms.positions.rowwise().minCoeff()};
    Domain domain(MPI_COMM_WORLD, volume, {2, 2, 2}, {1, 1, 1});
    // switch to decomposed state
    domain.enable(atoms);

    // first, heat up the system gently to some initial temperature so we don't
    // have to start from zero
    double e_pot{0.0};
    for (size_t i{0}; i < TIMESTEPS; i++) {
        // perform timestep
        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, DT,
                     atoms.masses);
        domain.update_ghosts(
            atoms,
            CUTOFF_RADIUS * 2.); // border radius is twice the cutoff radius
        neighbour_list.update(atoms, CUTOFF_RADIUS);
        double e_pot{ducastelle(atoms, neighbour_list)};
        verlet_step2(atoms.velocities, atoms.forces, DT, atoms.masses);
        domain.exchange_atoms(atoms);

        if (i % PLOT_EVERY == 0) {
        }
    }

    MPI_Finalize();
    return 0;
}