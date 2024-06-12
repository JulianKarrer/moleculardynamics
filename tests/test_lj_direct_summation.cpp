/*
 * Copyright 2021 Lars Pastewka
 *
 * ### MIT license
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <gtest/gtest.h>

#include "lj_direct_summation.h"
#include <verlet.h>

TEST(LJDirectSummationTest, Forces) {
    constexpr int nb_atoms = 10;
    constexpr double epsilon =
        0.7; // choose different to 1 to pick up missing factors
    constexpr double sigma = 0.3;
    constexpr double delta = 0.0001; // difference used for numerical (finite
                                     // difference) computation of forces

    Atoms atoms(nb_atoms);
    atoms.positions.setRandom(); // random numbers between -1 and 1

    // compute and store energy of the indisturbed configuration
    // cast summation result to void so the unused warning is surpressed
    (void)lj_direct_summation(atoms, epsilon, sigma);
    Forces_t forces0{atoms.forces};

    // loop over all atoms and compute forces from a finite differences
    // approximation
    for (int i{0}; i < nb_atoms; ++i) {
        // loop over all Cartesian directions
        for (int j{0}; j < 3; ++j) {
            // move atom to the right
            atoms.positions(j, i) += delta;
            double eplus{lj_direct_summation(atoms, epsilon, sigma)};
            // move atom to the left
            atoms.positions(j, i) -= 2 * delta;
            double eminus{lj_direct_summation(atoms, epsilon, sigma)};
            // move atom back to original position
            atoms.positions(j, i) += delta;

            // finite-differences forces
            double fd_force{-(eplus - eminus) / (2 * delta)};

            // check whether finite-difference and analytic forces agree
            if (abs(forces0(j, i)) > 1e-10) {
                EXPECT_NEAR(abs(fd_force - forces0(j, i)) / forces0(j, i), 0,
                            1e-5);
            } else {
                EXPECT_NEAR(fd_force, forces0(j, i), 1e-10);
            }
        }
    }
}

// Test if the potential has the minimum at the right distance and energy
// by setting 2 particles at r_min and checking if
// 1) they remain stationary
// 2) the potential energy at the minimum is -epsilon
TEST(LJDirectSummationTest, CorrectMinimumAndEquilibriumDistance) {
    constexpr int nb_atoms = 2;
    constexpr double epsilon = 0.7;
    constexpr double sigma = 0.3;
    Atoms atoms(nb_atoms);
    double r_min = pow(2, 1. / 6.) * sigma;
    atoms.positions(0, 0) = -r_min / 2.;
    atoms.positions(0, 1) = r_min / 2.;
    double dt{0.001};

    (void)lj_direct_summation(atoms, epsilon, sigma);
    for (int i{0}; i < 5000; ++i) {
        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, dt,
                     atoms.masses);
        double potential{lj_direct_summation(atoms, epsilon, sigma)};
        verlet_step2(atoms.velocities, atoms.forces, dt, atoms.masses);
        // assert that particles at r_min distance stay at that distance
        Vec3_t distance{atoms.positions.col(0) - atoms.positions.col(1)};
        EXPECT_NEAR(distance.norm(), r_min, 1e-6);
        // the potential energy at r_min shoould be -epsilon
        EXPECT_NEAR(potential, -epsilon, 1e-6);
    }
}

// Test if momentum is conserved
TEST(LJDirectSummationTest, MomentumConserved) {
    constexpr int nb_atoms{10};
    constexpr double epsilon{0.7};
    constexpr double sigma{0.3};
    constexpr double v_init{0.1};
    constexpr double dt{0.001};
    Atoms atoms(nb_atoms);
    atoms.positions.setRandom();

    // calculate the momentum of the randomly initialized atoms
    Vec3_t momentum_init{Vec3_t(0., 0., 0.)};
    for (int i{0}; i < nb_atoms; ++i) {
        Vec3_t v{Vec3_t(0., 0., 0.)};
        v.setRandom();
        v.normalize();
        ((Vec3_t)atoms.velocities.col(i)) = v * v_init;
        momentum_init += ((Vec3_t)atoms.velocities.col(i)) * atoms.masses(i);
    }

    (void)lj_direct_summation(atoms, epsilon, sigma);
    // execute 10_000 simulation steps
    for (int i{0}; i < 10000; ++i) {
        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, dt,
                     atoms.masses);
        (void)lj_direct_summation(atoms, epsilon, sigma);
        verlet_step2(atoms.velocities, atoms.forces, dt, atoms.masses);
    }

    // calculate the momentum again and compare it to the initial value
    Vec3_t momentum_after{Vec3_t(0., 0., 0.)};
    for (int i{0}; i < nb_atoms; ++i) {
        momentum_after += ((Vec3_t)atoms.velocities.col(i)) * atoms.masses(i);
    }
    EXPECT_NEAR(momentum_after.norm(), momentum_init.norm(), 1e-6);
}