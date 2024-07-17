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
const double TEMPERATURE{100.};
const size_t TIMESTEPS{10000};
const size_t EQUILIBRATE_FOR{5000};
const size_t SCALE_EVERY{50};
const size_t PLOT_EVERY{TIMESTEPS / 100};
const double SPACING{4.079 / sqrt(2)};
const double BORDER{5 * SPACING};

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    // initialize gold cluster from file
    auto [names, positions,
          velocities]{read_xyz_with_velocities("whisker_small.xyz")};
    Atoms atoms{Atoms(names, positions)};
    // miniscule random velocities to get kinetic energy in all modes
    atoms.velocities.setRandom();
    atoms.velocities *= 1e-8;
    double m{ELEM_NAME_TO_MASS.at("Au")};
    std::ofstream traj("whisker_traj.xyz");

    // create neighbour list and initialize forces
    NeighborList neighbour_list;
    neighbour_list.update(atoms, CUTOFF_RADIUS);
    (void)ducastelle(atoms, neighbour_list);

    // find the number of MPI processes
    int nb_processes{0};
    MPI_Comm_size(MPI_COMM_WORLD, &nb_processes);

    // move the atoms such that the lowest z-value is zero and the volume
    // specified later surrounds the whisker in x,y directions
    Vec3_t min_pos{atoms.positions.rowwise().minCoeff()};
    atoms.positions.row(0) -= min_pos(0);
    atoms.positions.row(1) -= min_pos(1);
    atoms.positions.row(2) -= min_pos(2);
    min_pos = atoms.positions.rowwise().minCoeff();
    Vec3_t max_pos{atoms.positions.rowwise().maxCoeff()};
    Vec3_t volume{max_pos - min_pos};
    atoms.positions.row(0) += BORDER;
    atoms.positions.row(1) += BORDER;
    // increase the volume in x and y direction as much as desired so atoms
    // don't go out of bounds, only increase in z-direction by some small
    // epsilon so atoms are not lost when the domain is activated.
    volume(0) += 2.0 * BORDER;
    volume(1) += 2.0 * BORDER;
    volume(2) += SPACING / 2.0;

    // here, only decompose along the z-axis,  where the whisker extends and use
    // periodic boundary in z-direction
    Domain domain(MPI_COMM_WORLD, volume, {1, 1, nb_processes}, {0, 0, 1});

    if (domain.rank() == 0) {
        // write csv header to stdout
        std::cout << "temp,rate,v_x,v_y,v_z,sig_zz,e_pot" << std::endl;

        // std::cout << std::endl << volume << std::endl;
        // std::cout << std::endl << min_pos << std::endl;
        // std::cout << std::endl << max_pos << std::endl;
        // write_xyz(traj, atoms);
    }

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
    for (size_t i{0}; i < TIMESTEPS + EQUILIBRATE_FOR; i++) {
        if (i > EQUILIBRATE_FOR) {
            // continuously rescale the domain to stretch it in z-direction
            volume(2) += SPACING / (double)SCALE_EVERY;
            domain.scale(atoms, volume);
        }

        // update positions and velocities based on forces
        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, DT, m);
        // communicate new ghost particles and update forces
        domain.exchange_atoms(atoms);
        domain.update_ghosts(atoms, CUTOFF_RADIUS * 2.);
        neighbour_list.update(atoms, CUTOFF_RADIUS);
        auto [e_pot_l, stress_l] = ducastelle_stress(
            atoms, neighbour_list, domain.nb_local(), CUTOFF_RADIUS);
        // update masses using the newly computed forces
        verlet_step2(atoms.velocities, atoms.forces, DT, m);
        berendsen_thermostat_decomposed(atoms, TEMPERATURE, DT, 100 * DT,
                                        domain.nb_local(), m);

        if (i > EQUILIBRATE_FOR && i % PLOT_EVERY == 0) {
            double stress_g{
                MPI::allreduce(stress_l, MPI_SUM, domain.communicator()) * 1.0 /
                volume.prod()};
            double e_pot_global{
                MPI::allreduce(e_pot_l, MPI_SUM, domain.communicator())};

            domain.disable(atoms);
            if (domain.rank() == 0) {
                std::cout << std::fixed << std::setprecision(15)
                          << temperature_cur(atoms) << ","
                          << SPACING / (double)SCALE_EVERY << "," << volume(0)
                          << "," << volume(1) << "," << volume(2) << ","
                          << stress_g << "," << e_pot_global << std::endl;
                write_xyz(traj, atoms);
            }
            domain.enable(atoms);
            domain.update_ghosts(atoms, CUTOFF_RADIUS * 2.);
            neighbour_list.update(atoms, CUTOFF_RADIUS);
            (void)ducastelle(atoms, neighbour_list);
        }
    }

    if (domain.rank() == 0)
        std::cout << std::endl;
    MPI_Finalize();
    return 0;
}