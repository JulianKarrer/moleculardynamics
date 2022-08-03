#include "verlet.h"
#include <Eigen/Dense>
#include <iostream>

#ifdef USE_MPI
#include <mpi.h>
#endif

using Eigen::MatrixXd;

int main(int argc, char *argv[]) {
    int rank = 0, size = 1;

    // Below is some MPI code, try compiling with `cmake -DUSE_MPI=ON ..`
#ifdef USE_MPI
    MPI_Init(&argc, &argv);

    // Retrieve process infos
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

    std::cout << "Hello I am rank " << rank << " of " << size << "\n";

    // Below is some Eigen code
    MatrixXd m(2, 2);
    m(0, 0) = 3;
    m(1, 0) = 2.5;
    m(0, 1) = -1;
    m(1, 1) = m(1, 0) + m(0, 1);

    // C++11 style loop over elements of an Eigen matrix
    for (auto&& value : m.reshaped())
        value += 1.;

    // Output matrix values
    if (rank == 0)
        std::cout << m << std::endl;

    // Below is some code from our MD library
    double x, y, z, vx, vy, vz;
    scalar_verlet_step1(x, y, z, vx, vy, vz, 0, 0, 0, 1);
    scalar_verlet_step2(vx, vy, vz, 0, 0, 0, 1);

#ifdef USE_MPI
    MPI_Finalize();
#endif

    return 0;
}
