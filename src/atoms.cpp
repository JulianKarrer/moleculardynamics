#include "atoms.h"
#include <algorithm>

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