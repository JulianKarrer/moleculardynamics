#include "lj.h"
#include "lj_direct_summation.h"
#include "neighbors.h"

/// @brief Modified `lj_potential` where the energy is shifted up such that V=0
/// at the cutoff distance
/// @param r
/// @param epsilon
/// @param sigma
/// @return
double lj_cutoff_potential(double r, double epsilon, double sigma,
                           double cutoff) {

    return r > cutoff ? 0.0
                      : (lj_potential(r, epsilon, sigma) +
                         lj_potential(r, epsilon, cutoff));
}

/// @brief compute forces derived from the truncated and shifted Lennard-Jones
/// potential using an O(n) neighbor search and store them in the `forces` field
/// of the `Atoms` class, overwriting existing entries.
/// @param atoms the atoms that interact. The field `forces` is written to.
/// @param epsilon the binding energy of the Lennard-Jones potential
/// (commonly given in electron volts)
/// @param sigma the distance for which the potential is zero
/// (commonly given in Ångström)
/// @param cutoff the cutoff distance for the potential, given in the same units
/// as `sigma`
/// @return the total potential energy of the system.
double ljts(Atoms &atoms, double epsilon, double sigma, double cutoff) {
    double energy_sum{0.};
    atoms.forces.setZero();
    // initialize neighbour list
    NeighborList neighbor_list;
    neighbor_list.update(atoms, cutoff);
    for (auto [i, j] : neighbor_list) {
        if (i < j) {
            Vec3_t r_ij{atoms.positions.col(i) - atoms.positions.col(j)};
            double r_ij_mag{r_ij.norm()};
            Vec3_t r_hat_ij{r_ij.normalized()};
            // f_ij = - \nabla E_{pot} = \sum_j V'(|r_ij|) r_hat_ij
            Vec3_t f_ij{r_hat_ij *
                        lj_potential_derivative(r_ij_mag, epsilon, sigma)};
            // since we loop over pairs, deposit f_ij at i and -f_ij at j
            Vec3_t force_i{atoms.forces.col(i)};
            Vec3_t force_j{atoms.forces.col(j)};
            atoms.forces.col(i) = force_i - f_ij;
            atoms.forces.col(j) = force_j + f_ij;
            // add the pairwise potential to the accumulated total potential
            // energy
            energy_sum += lj_cutoff_potential(r_ij_mag, epsilon, sigma, cutoff);
        }
    }
    // return the total potential energy of the system
    return energy_sum;
}