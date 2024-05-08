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

size_t Atoms::nb_atoms() const {
    return positions.cols();
};