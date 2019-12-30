#define CATCH_CONFIG_MAIN  // This tells Catch to provide a =
#include "catch.hpp"

#include <chrono>

#include "src/species_advance/species_advance.h"
#include "src/vpic/vpic.h"
#include "test/integrated/particle_push/advance_p.h"


class vpic_test_simulation : public vpic_simulation
{
    public:
        void test()
        {
            double L  = 1e2;
#ifndef NUM_PARTICLES
#define NUM_PARTICLES 1048576 //32768
#endif
            int npart = NUM_PARTICLES; // This is independent of the number of cells
            int nstep = 100;

            define_units( 1, 1 );
            define_timestep( 1 );

            define_periodic_grid( 0, 0, 0,   // Grid low corner
                                  L, L, L,   // Grid high corner
                                  8, 8, 8,   // Grid resolution
                                  1, 1, 1 ); // Processor configuration

            define_material( "vacuum", 1.0, 1.0, 0.0 );
            define_field_array();

            field(1,1,1).ex  = 1;
            field(1,2,1).ex  = 1;
            field(1,1,2).ex  = 1;
            field(1,2,2).ex  = 1;

            field(1,1,1).ey  = 2;
            field(1,1,2).ey  = 2;
            field(2,1,1).ey  = 2;
            field(2,1,2).ey  = 2;

            field(1,1,1).ez  = 3;
            field(2,1,1).ez  = 3;
            field(1,2,1).ez  = 3;
            field(2,2,1).ez  = 3;

            species_t * sp =
                define_species( "test_species", 1., 1., npart, npart, 0, 0 );

            for (int i = 0; i < npart; i++)
            {
                float x = uniform( rng(0), 0, L);
                float y = uniform( rng(0), 0, L);
                float z = uniform( rng(0), 0, L);

                inject_particle( sp , x, y, z, 0., 0., 0., 1., 0., 0);
            }

            load_interpolator_array( interpolator_array, field_array );

#ifdef ENABLE_SORT
            std::cout << "Sort is enabled." << std::endl;
            sort_p( sp );
#endif
            auto start = std::chrono::system_clock::now();

            for( int n=0; n<nstep; n++)
            {
                uncenter_p( sp, interpolator_array );
                center_p( sp, interpolator_array );
            }
            auto end = std::chrono::system_clock::now();

            std::chrono::duration<double> elapsed_seconds = end-start;
            auto time_per_p = elapsed_seconds.count() / (sp->np*nstep);
            std::cout << "Elapsed time: " << elapsed_seconds.count() << "" << std::endl;
            std::cout << "Particles: " << sp->np << ", Steps: " << nstep << std::endl;
            std::cout << "Time per particle step: " << time_per_p << ", Particle/s: " << 1.0/time_per_p << std::endl;

        }
};

TEST_CASE( "repeatedly center and uncenter particles in a timed loop", "[uncenter]" )
{

#ifndef NUM_THREADS
#define NUM_THREADS 1
#endif

    int pargc = 1;

    if (NUM_THREADS > 1) { pargc = 3; }

    char str[] = "bin/vpic";
    char **pargv = (char **) malloc( (pargc+1) *sizeof(char **));
    pargv[0] = str;

    if (NUM_THREADS > 1)
    {
        char str2[] = "--tpp";
        char str3[] = "" EXPAND_AND_STRINGIFY(NUM_THREADS);
        pargv[1] = str2;
        pargv[2] = str3;
    }

    boot_services( &pargc, &pargv );

    vpic_test_simulation* simulation = new vpic_test_simulation();
    simulation->test();

    delete simulation;
    if( world_rank==0 ) log_printf( "normal exit\n" );

    halt_mp();
}
