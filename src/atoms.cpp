#include "atoms.h"
#include <algorithm>
#include <iostream>

/// @brief Initialize a set of `n` atoms. Masses are 1, all vectors are 0.
/// @param n number of atoms
Atoms::Atoms(const size_t n)
    : positions{3, n},
      velocities{3, n},
      forces{3, n},
      masses{n} {
    positions.setZero();
    velocities.setZero();
    forces.setZero();
    masses.setOnes();
};

/// @brief Initialize a set of `n` atoms with masses of `1` and zero velocities
/// on a body-centered cubic lattice, given a spacing between grid points.
/// The bounding volume of the lattice is optimized to be as close to cube
/// shaped as possible and centred around the origin.
///
/// https://www5.in.tum.de/lehre/vorlesungen/sci_compII/ss06/lectures/01_molecular_dynamics.pdf
/// (slide 18)
/// @param n number of atoms
/// @param spacing 0.5*spacing between atoms in the lattice
Atoms::Atoms(const size_t n, double spacing)
    : positions{3, n},
      velocities{3, n},
      forces{3, n},
      masses{n} {
    positions.setZero();
    velocities.setZero();
    forces.setZero();
    masses.setOnes();

    // decompose uint N into uints x,y,z such that x*y*z=N and
    // max(x,y,z)-min(x,y,z) is minimal (~ find the must cube-y cuboid of
    // integer sizes and volume N).
    int root{(int)round(cbrt(n))};
    // the spacing around the cube root around to search for decompositions
    // the search complexity is O((2.*searchRadiusMax)^2) = O(1)
    int searchRadiusMax{100};
    // remember the best guess so far, start from a valid guess
    size_t best_x, best_y = 1;
    size_t best_z{n};
    size_t best_cubeness{n - 1};
    // any factor of the decomposition must be in [1; n]
    size_t lower{(size_t)std::max(root - searchRadiusMax, 1)};
    size_t upper{(size_t)std::min(root + searchRadiusMax, (int)n)};
    for (size_t x{lower}; x <= upper; x++) {
        for (size_t y{lower}; y <= upper; y++) {
            size_t z{n / (x * y)};
            // the goodness of the decomposition is the distance between the
            // maximum and minimum factor
            size_t cubeness{std::max(x, std::max(y, z)) -
                            std::min(x, std::min(y, z))};
            if (n == x * y * z && cubeness < best_cubeness) {
                best_x = x;
                best_y = y, best_z = z, best_cubeness = cubeness;
            }
        }
    }

    std::cout << "decomposition: " << best_x << "x " << best_y << "y " << best_z
              << "z " << "\n";
    // having found an XYZ decomposition, initialize the lattice
    size_t i{0};
    for (size_t x{0}; x < best_x; x++) {
        for (size_t y{0}; y < best_y; y++) {
            for (size_t z{0}; z < best_z; z++) {
                assert(i < n);
                double delta = spacing; // * M_SQRT2;
                // offset xy plane by 0.5 in each z-step to form a
                double offset{z % 2 == 0 ? delta * 0.25 : -delta * 0.25};
                // place the particles
                positions(0, i) = x * delta - best_x * delta * 0.5 + offset;
                positions(1, i) = y * delta - best_y * delta * 0.5 + offset;
                positions(2, i) = z * delta - best_z * delta * 0.5;
                i++;
            }
        }
    }
    std::cout << "grid ok\n";
};

/// @brief Initialize a set of atoms at positions `p`. Masses are 1, all other
/// vectors are 0.
/// @param p initial positions of the atoms
Atoms::Atoms(const Positions_t &p)
    : positions{p},
      velocities{3, p.cols()},
      forces{3, p.cols()},
      masses{p.cols()} {
    velocities.setZero();
    forces.setZero();
    masses.setOnes();
};

/// @brief Initialize a set of atoms at positions `p` with velocities `v`.
/// Masses are 1, all other vectors are 0.
/// @param p initial positions of the atoms
/// @param v initial velocities of the atoms
Atoms::Atoms(const Positions_t &p, const Velocities_t &v)
    : positions{p},
      velocities{v},
      forces{3, p.cols()},
      masses{p.cols()} {
    assert(p.cols() == v.cols());
    forces.setZero();
    masses.setOnes();
};

/// @brief Initialize a set of atoms at positions `p` with velocities `v` with
/// given masses `m`. All other vectors are 0.
/// @param p initial positions of the atoms
/// @param v initial velocities of the atoms
/// @param m initial masses of the atoms
Atoms::Atoms(const Positions_t &p, const Velocities_t &v, const Masses_t &m)
    : positions{p},
      velocities{v},
      forces{3, p.cols()},
      masses{m} {
    assert(p.cols() == v.cols());
    forces.setZero();
};

/// @brief Query the number of atoms in the system.
/// @return the number of atoms
size_t Atoms::nb_atoms() const {
    return positions.cols();
};

/// @brief Compute the total kinetic energy of all atoms $E_kin = \\sum_i m_i
/// v_i\\cdot v_i$
/// @return the total kinetic energy of the system
double Atoms::kinetic_energy() {
    // use `colwise()` and `transpose()` to perform the sum using only Eigen
    // functions
    return ((Eigen::ArrayXd)(velocities.colwise().squaredNorm().transpose() *
                             masses * 0.5))
        .sum();
}