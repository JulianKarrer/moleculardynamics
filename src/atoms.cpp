#include "atoms.h"

Atoms::Atoms(const Positions_t &p)
    : positions{p},
      velocities{3, p.cols()},
      forces{3, p.cols()} {
    velocities.setZero();
    forces.setZero();
};

Atoms::Atoms(const Positions_t &p, const Velocities_t &v)
    : positions{p},
      velocities{v},
      forces{3, p.cols()} {
    assert(p.cols() == v.cols());
    forces.setZero();
};

Atoms::Atoms(const size_t n) : positions{3, n}, velocities{3, n}, forces{3, n} {
    positions.setZero();
    velocities.setZero();
    forces.setZero();
};

size_t Atoms::nb_atoms() const {
    return positions.cols();
};

double Atoms::kinetic_energy(const double mass) {
    double sum{0.};
    for (size_t i{0}; i < nb_atoms(); ++i) {
        Vec3_t v_i{velocities.col(i)};
        sum += 0.5 * mass * v_i.dot(v_i);
    }
    return sum;
}