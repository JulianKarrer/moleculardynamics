#include <gtest/gtest.h>
#include <lj_direct_summation.h>
#include <thermostat.h>
#include <verlet.h>

/// @brief Implements the test of the thermostat for the different initial
/// conditions defined in the `TEST`s in the following
void berendsen_test(Atoms &atoms, double t_init, double t_wish, double dt,
                    double assert_tolerance) {
    // Perform 10k simulation steps to make sure the system should have
    // calibrated
    for (uint step = 0; step < 10000; step++) {
        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, dt,
                     atoms.masses);
        (void)lj_direct_summation(atoms, 1., 1.);
        verlet_step2(atoms.velocities, atoms.forces, dt, atoms.masses);
        berendsen_thermostat(atoms, t_wish, dt, 100 * dt);
    }
    double t_end{temperature_cur(atoms)};

    // Print results for debugging and to catch errors
    std::cout << "t_init: " << t_init << "t_wish: " << t_wish
              << "t_end: " << t_end << "\n";

    // Assert that the system has equilibrated.
    ASSERT_NEAR(t_end, t_wish, assert_tolerance);
}

/// @brief Tests if the Berendsen thermostat equlibrates a system of 10 randomly
/// positioned atoms in [-1;1] with random unit initial velocities to twice the
/// initial temperature
TEST(BerendsenThermosTest, Equilibrate) {
    // intialize state of the system
    int nb_atoms{10};
    Positions_t positions(3, nb_atoms);
    positions.setRandom();
    Velocities_t velocities(3, nb_atoms);
    velocities.setRandom();
    velocities.colwise().normalize();
    Atoms atoms{Atoms(positions, velocities)};

    // Set target temperature and time step size
    double t_init{temperature_cur(atoms)};
    double t_wish{t_init * 2};
    double dt{0.01};

    // set low tolerance
    berendsen_test(atoms, t_init, t_wish, dt, 1e-4);
}

/// @brief Check if the above test also works from zero temperature to a target
/// temperature
TEST(BerendsenThermosTest, FromZero) {
    int nb_atoms{10};
    Positions_t positions(3, nb_atoms);
    positions.setRandom();
    Atoms atoms{Atoms(positions)};

    // check if the initial velocity is, in fact, zero
    double t_init{temperature_cur(atoms)};
    ASSERT_FLOAT_EQ(t_init, 0.0);
    // set the target temperature to 1000, otherwise perform the same test
    double t_wish{1000};
    double dt{0.01};

    // higher tolerance since initial conditions are more tricky
    berendsen_test(atoms, t_init, t_wish, dt, 0.1);
}

/// @brief Check if the thermostat can cool the system to zero K. This should
/// result in standstill. The initially densly packed atoms with random
/// unit velocities should also expand, meaning the bounding volume increases.
TEST(BerendsenThermosTest, ToZero) {
    int nb_atoms{10};
    Positions_t positions(3, nb_atoms);
    positions.setRandom();
    positions *= 0.5;
    Velocities_t velocities(3, nb_atoms);
    velocities.setRandom();
    velocities.colwise().normalize();
    Atoms atoms{Atoms(positions, velocities)};

    double initial_volume{(atoms.positions.rowwise().maxCoeff() -
                           atoms.positions.rowwise().minCoeff())
                              .prod()};
    double t_init{temperature_cur(atoms)};
    // set the target temperature to 0, otherwise perform the same test
    double t_wish{0};
    double dt{0.01};

    // very low tolerance, since this state should be easily reached
    berendsen_test(atoms, t_init, t_wish, dt, 1e-10);
    double final_volume{(atoms.positions.rowwise().maxCoeff() -
                         atoms.positions.rowwise().minCoeff())
                            .prod()};
    // assert that the bounding volume has expanded.
    std::cout << "initial volume: " << initial_volume
              << "final volume: " << final_volume << "\n";
    ASSERT_GT(final_volume, initial_volume);
}
