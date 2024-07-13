#include "lj_direct_summation.h"

/// @brief compute the Lennard-Jones potential as a function of atom distance
/// `r`
/// @param r distance between the two atoms
/// @param epsilon the binding energy of the Lennard-Jones potential
/// (commonly given in electron volts)
/// @param sigma the distance for which the potential has its minimum
/// (commonly given in Ångström)
/// @return the Lennard-Jones potential at the distance `r`
double lj_potential(double r, double epsilon, double sigma) {
    double sig_over_r{sigma / r};
    double s_2{sig_over_r * sig_over_r};
    double s_4{s_2 * s_2};
    double s_6{s_4 * s_2};
    double s_8{s_4 * s_4};
    double s_12{s_8 * s_4};
    // implement the potential as shown in the lecture
    return 4. * epsilon * (s_12 - s_6);
}

/// @brief the derivative of the Lennard-Jones potential with respect to the
/// distance r
/// @param r distance between the two atoms
/// @param epsilon the binding energy of the Lennard-Jones potential
/// (commonly given in electron volts)
/// @param sigma the distance for which the potential is zero
/// (commonly given in Ångström)
/// @return the derivative of the Lennard-Jones potential
double lj_potential_derivative(double r, double epsilon, double sigma) {
    double sig_over_r = sigma / r;
    double s_2{sig_over_r * sig_over_r};
    double s_4{s_2 * s_2};
    double s_6{s_4 * s_2};
    double s_8{s_4 * s_4};
    double s_12{s_8 * s_4};
    // implement the derivative of the potential (note the additional /r terms
    // and multipliers)
    return epsilon * (-48. * s_12 / r + 24. * s_6 / r);
}

/// @brief compute forces derived from the Lennard-Jones potential through a
/// O(n^2) summation over unique pairs in (~ 0.5n^2 o(atoms.nb_atoms() - 1)
/// @param sigma the distance for which the potential is zero
/// (commonly given in Ångström)
/// @return the total potential energy of the system.
double lj_direct_summation(Atoms &atoms, double epsilon, double sigma) {
    // reserve a variable to store the total energy accumulated in the loop
    double energy_sum{0.};
    atoms.forces.setZero();
    // loop over unique pairs by asserting i<j in the reversed inner loop
    for (size_t i{0}; i < atoms.nb_atoms(); ++i) {
        for (size_t j{i + 1}; j < atoms.nb_atoms(); ++j) {
            Vec3_t r_ij{atoms.positions.col(i) - atoms.positions.col(j)};
            double r_ij_mag{r_ij.norm()};
            Vec3_t r_hat_ij{r_ij / r_ij_mag};
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
            energy_sum += lj_potential(r_ij_mag, epsilon, sigma);
        }
    }
    // return the total potential energy of the system
    return energy_sum;
}
