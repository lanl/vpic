#ifndef VPIC_INTERNAL_DIAGNOSTICS_H
#define VPIC_INTERNAL_DIAGNOSTICS_H

#include <mpi.h>
#include <algorithm>
#include <random>

class FULL_GATHER;

// TODO: the MPI stuff could be moved into the message passing layer of VPIC
class internal_diagnostics_t
{

    public:
        /**
         * @brief WARNING: This likely implies MPI comms
         */
        template <class T> void detect_load_imbalance(int rank, int nproc, double in_time, int np_total)
        {
            detect_load_imbalance<FULL_GATHER>(rank, nproc, in_time, np_total);
        }
};


// inline prevents multiply defined variants, but this should be moved out of
// the header and into a sensible translation unit
template<> inline void
internal_diagnostics_t::detect_load_imbalance<FULL_GATHER>(int rank, int nproc, double in_time, int np_total)
{
    /*
    // Grab a random number in some range
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 eng(rd()); // seed the generator
    std::uniform_int_distribution<> distr(25, 63); // define the range
    float in_time = distr(eng);
    std::cout << "Rank " << rank << " picked number " << in_time << std::endl;
    */

    // TODO: use c++ memory allocation

    const int root = 0;
    const int num_vals = 2;

    using pair_t = std::pair<double, double>;
    pair_t send = { in_time, np_total };
    std::vector< pair_t > recv(nproc);

    MPI_Gather(
            &send,
            num_vals,
            MPI_DOUBLE,
            recv.data(),
            num_vals,
            MPI_DOUBLE,
            root,
            MPI_COMM_WORLD
            );

    if ( rank == root )
    {
        std::cout << "My time " << in_time << " np " << np_total << std::endl;
        const auto p = std::minmax_element(recv.begin(), recv.end());
        auto min = p.first->first;
        auto max = p.second->first;

        std::cout <<
            "Max :" << max <<
            " Min: " << min <<
            " Ratio: " << min/max <<
            std::endl;
    }
}

#endif // end guard
