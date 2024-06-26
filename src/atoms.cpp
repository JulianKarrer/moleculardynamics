#include "atoms.h"
#include <algorithm>
#include <iostream>

/// @brief Initialize positions on an oblique lattice, given a spacing between
/// grid points. The bounding volume of the lattice is optimized to be close to
/// cube shaped and centred around the origin.
///
/// https://www5.in.tum.de/lehre/vorlesungen/sci_compII/ss06/lectures/01_molecular_dynamics.pdf
/// (slide 18)
/// @param n number of atoms
/// @param spacing 0.5*spacing between atoms in the lattice
void initialize_lattice(Positions_t &positions, size_t n, double spacing) {
    size_t cube_length{(size_t)std::ceil(std::cbrt(n))};
    size_t i{0};
    for (size_t x{0}; x < cube_length + 1; x++) {
        for (size_t y{0}; y < cube_length + 1; y++) {
            for (size_t z{0}; z < cube_length + 1; z++) {
                if (i >= n) {
                    return;
                }
                double delta = spacing; // * M_SQRT2;
                // offset xy plane by 0.5 in each z-step to form a
                double offset{z % 2 == 0 ? delta * 0.25 : -delta * 0.25};
                // place the particles
                positions(0, i) =
                    x * delta - cube_length * delta * 0.5 + offset;
                positions(1, i) =
                    y * delta - cube_length * delta * 0.5 + offset;
                positions(2, i) =
                    z * delta / M_SQRT2 - cube_length * delta / M_SQRT2 * 0.5;
                i++;
            }
        }
    }
}

/// @brief Initialize a set of `n` atoms. Masses are 1, all vectors are 0.
/// @param n number of atoms
Atoms::Atoms(const size_t n)
    : positions{3, n},
      velocities{3, n},
      forces{3, n},
      masses{n},
      names(n, "H") {
    positions.setZero();
    velocities.setZero();
    forces.setZero();
    masses.setOnes();
};

/// @brief Initialize a set of `n` atoms with masses of `1` and zero velocities
/// on a lattice, given a spacing between grid points.
/// The bounding volume of the lattice is optimized to be as close to cube
/// shaped as possible and centred around the origin.
/// @param n number of atoms
/// @param spacing 0.5*spacing between atoms in the lattice
Atoms::Atoms(const size_t n, double spacing)
    : positions{3, n},
      velocities{3, n},
      forces{3, n},
      masses{n},
      names(n, "H") {
    positions.setZero();
    velocities.setZero();
    forces.setZero();
    masses.setOnes();
    initialize_lattice(positions, n, spacing);
};

/// @brief Initialize a set of atoms at positions `p`. Masses are 1, all other
/// vectors are 0.
/// @param p initial positions of the atoms
Atoms::Atoms(const Positions_t &p)
    : positions{p},
      velocities{3, p.cols()},
      forces{3, p.cols()},
      masses{p.cols()},
      names(p.cols(), "H") {
    velocities.setZero();
    forces.setZero();
    masses.setOnes();
};

/// @brief Initialize a set of atoms at positions `p` with `names` that
/// correspond to CamelCase names of elements and determine the mass of the
/// particle in units of `u`. Initial velocities are zero.
/// @param p initial positions of the atoms
Atoms::Atoms(const Names_t names, Positions_t &p)
    : positions{p},
      velocities{3, p.cols()},
      forces{3, p.cols()},
      masses{p.cols()},
      names{names} {
    velocities.setZero();
    forces.setZero();
    for (auto i{0}; i < p.cols(); i++) {
        masses(i) = ELEM_NAME_TO_MASS.at(names[i]);
    }
};

/// @brief Initialize a set of atoms at positions `p` with velocities `v`.
/// Masses are 1, all other vectors are 0.
/// @param p initial positions of the atoms
/// @param v initial velocities of the atoms
Atoms::Atoms(const Positions_t &p, const Velocities_t &v)
    : positions{p},
      velocities{v},
      forces{3, p.cols()},
      masses{p.cols()},
      names(p.cols(), "H") {
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
      masses{m},
      names(p.cols(), "H") {
    assert(p.cols() == v.cols());
    forces.setZero();
};

/// @brief Initialize a set of atoms at positions `p` with velocities `v` with
/// given names containing their atomic symbols. All other vectors are 0.
/// @param p initial positions of the atoms
/// @param v initial velocities of the atoms
/// @param names names of the atoms, ie. "H", "Au", ...
Atoms::Atoms(const Positions_t &p, const Velocities_t &v, const Names_t names)
    : positions{p},
      velocities{v},
      forces{3, p.cols()},
      masses{p.cols()},
      names(names) {
    assert(p.cols() == v.cols());
    forces.setZero();
    for (auto i{0}; i < p.cols(); i++) {
        masses(i) = ELEM_NAME_TO_MASS.at(names[i]);
    }
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