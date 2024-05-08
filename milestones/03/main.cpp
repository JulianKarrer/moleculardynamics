#include "verlet.h"
#include <cmath>
#include <iostream>

class Atoms {
  public:
    Positions_t positions;
    Velocities_t velocities;
    Forces_t forces;

    Atoms(const Positions_t &p)
        : positions{p},
          velocities{3, p.cols()},
          forces{3, p.cols()} {
        velocities.setZero();
        forces.setZero();
    }

    Atoms(const Positions_t &p, const Velocities_t &v)
        : positions{p},
          velocities{v},
          forces{3, p.cols()} {
        assert(p.cols() == v.cols());
        forces.setZero();
    }

    size_t nb_atoms() const {
        return positions.cols();
    }
};

int main(int argc, char *argv[]) {
    std::cout << "TEST: milestone 3 compiles!\n";
    return 0;
}
