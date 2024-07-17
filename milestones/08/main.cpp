#include "atoms.h"
#include "mpi.h"
#include "mpi_support.h"
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
const double DT{1};
const double TEMPERATURE{500.};
const size_t TIMESTEPS{2000};
const size_t PLOT_EVERY{10};

/// @brief Communicate the current local potential and kinetic energies amongst
/// all MPI processes and let the rank 0 process write the global kinetic and
/// potential energies to stdout, where they may be piped into a csv.
/// @param atoms
/// @param e_pot_local
/// @param domain
/// @param m the mass of all atoms (we assume they are of the same kind)
/// @param nb_processes number of MPI processes
/// @param nb_global global number of atoms
void print_results(Atoms &atoms, const double e_pot_local, Domain &domain,
                   const double m, const int nb_processes,
                   const size_t nb_global) {
    // sum up energies using allreduce operations
    double e_kin_local{atoms.kinetic_energy(domain.nb_local(), m)};
    double e_pot_global{
        MPI::allreduce(e_pot_local, MPI_SUM, domain.communicator())};
    double e_kin_global{
        MPI::allreduce(e_kin_local, MPI_SUM, domain.communicator())};
    // assert that global values are less than local ones for the
    // negative potentials and global>local for the positive kinetic
    assert(e_kin_local <= e_kin_global);
    assert(e_pot_global <= e_pot_local);
    // energy only the first rank process may print the results
    if (domain.rank() == 0) {
        // write result to stdout
        std::cout << nb_processes << "," << nb_global << "," << std::fixed
                  << std::setprecision(15) << e_kin_global + e_pot_global << ","
                  << e_pot_global << "," << e_kin_global << std::endl;
    }
}

/// @brief Run simulations steps on the Mackay 923 Au initial configuration
/// using MPI to parallelize the computations, plotting the time evolution of
/// the Hamiltonian to confirm energy conservation or at least equal behaviour
/// between versions with different numbers of MPI processes.
int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    // initialize gold cluster from file
    auto [names, positions,
          velocities]{read_xyz_with_velocities("cluster_923.xyz")};
    Atoms atoms{Atoms(names, positions)};

    // fix the mass since all the atoms are gold and masses are not restored
    // when resizing between decomposed and unified domain
    double m{ELEM_NAME_TO_MASS.at("Au")};

    // create neighbour list and initialize forces
    NeighborList neighbour_list;
    neighbour_list.update(atoms, CUTOFF_RADIUS);
    (void)ducastelle(atoms, neighbour_list);

    // find the number of MPI processes
    int nb_processes{0};
    MPI_Comm_size(MPI_COMM_WORLD, &nb_processes);
    // initialize a domain for periodic boundaries and decomposition for
    // parallelization
    Vec3_t volume{atoms.positions.rowwise().maxCoeff() -
                  atoms.positions.rowwise().minCoeff()};
    // when in doubt and cbrt(nb_processes) is not an int, use more processes
    // for the x-axis
    int lower{(int)floor(cbrt((double)nb_processes))};
    int higher = nb_processes / (lower * lower);
    assert(higher * lower * lower == nb_processes);
    // use periodic boundaries
    Domain domain(MPI_COMM_WORLD, volume, {higher, lower, lower}, {0, 0, 0});

    // switch to decomposed state
    size_t nb_global{atoms.nb_atoms()};
    domain.enable(atoms);

    double e_pot_local{0.0};

    // initial force calculation
    domain.update_ghosts(atoms,
                         CUTOFF_RADIUS *
                             2.); // border radius is twice the cutoff radius
    neighbour_list.update(atoms, CUTOFF_RADIUS);
    (void)ducastelle(atoms, neighbour_list);
    // main simulation loop
    for (size_t i{0}; i < TIMESTEPS; i++) {
        // update positions and velocities based on forces
        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, DT, m);
        // communicate new ghost particles and update forces
        domain.exchange_atoms(atoms);
        domain.update_ghosts(atoms, CUTOFF_RADIUS * 2.);
        neighbour_list.update(atoms, CUTOFF_RADIUS);
        e_pot_local =
            ducastelle(atoms, neighbour_list, CUTOFF_RADIUS, domain.nb_local());
        // update masses using the newly computed forces
        verlet_step2(atoms.velocities, atoms.forces, DT, m);

        // exchange data on potential and kinetic energies every so often
        if (i % PLOT_EVERY == 0) {
            print_results(atoms, e_pot_local, domain, m, nb_processes,
                          nb_global);
        }
    }

    MPI_Finalize();
    return 0;
}