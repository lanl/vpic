#ifndef VPIC_INTERNAL_DIAGNOSTICS_H
#define VPIC_INTERNAL_DIAGNOSTICS_H

#include <mpi.h>
#include <algorithm>
#include <random>

class FULL_GATHER;
class NONE;

namespace Experimental {
    // TODO: convert this to a policy

    // TODO: the MPI stuff could be moved into the message passing layer of VPIC
    class internal_diagnostics_t
    {

        public:
            /**
             * @brief WARNING: This likely implies MPI comms
             */
            template <class T = NONE> void detect_load_imbalance(int rank, int nproc, double in_time, int np_total)
            {
                // Let it be a no-op
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

            // TODO: we probably have no guarantee that pair is locally in memory
            using pair_t = std::array<double, 2>;
            pair_t send = { in_time, (double)np_total };
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
                // This gives an array of array[2]
                const auto p = std::minmax_element(recv.begin(), recv.end());
                auto min_t = p.first[0][0];
                auto max_t = p.second[0][0];

                // This is the p for the proc which took most time
                auto min_p = p.first[0][1];
                auto max_p = p.second[0][1];

                std::cout <<
                    "Time " <<
                    " Max :" << max_t <<
                    " Min: " << min_t <<
                    " Ratio: " << min_t/max_t <<
                    std::endl;
                std::cout <<
                    "Np for: " <<
                    " Max_t :" << max_p <<
                    ", Min_t " << min_p <<
                    ", Ratio: " << min_p/max_p <<
                    std::endl;
            }
        }

} // end namespace expirmental

#endif // end guard
