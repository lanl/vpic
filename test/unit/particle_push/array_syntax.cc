#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()
#include "catch.hpp"

#include "src/species_advance/species_advance.h"
#include "src/vpic/vpic.h"
#include "test/integrated/particle_push/advance_p.h"

void vpic_simulation::user_diagnostics() {}

void
vpic_simulation::user_initialization( int num_cmdline_arguments,
                                      char ** cmdline_argument )
{
  double L  = 1e2;
  int npart = 127;
  int nstep = 100;

  define_units( 1, 1 );
  define_timestep( 1 );
  define_periodic_grid( 0, 0, 0,   // Grid low corner
                        L, L, L,   // Grid high corner
                        1, 1, 1,   // Grid resolution
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

  species_t* sp2 =
    define_species( "test_species2", 1., 1., npart, npart, 0, 0 );

  for (int i = 0; i < npart; i++)
  {
      float x = uniform( rng(0), 0, L);
      float y = uniform( rng(0), 0, L);
      float z = uniform( rng(0), 0, L);

      // Put two sets of particle in the exact same space
      inject_particle( sp2, x, y, z, 0., 0., 0., 1., 0., 0);
      inject_particle( sp , x, y, z, 0., 0., 0., 1., 0., 0);
  }

  // Create a second accumulator_array
  accumulator_array_t* accumulator_array2 = new_accumulator_array( grid );

  clear_accumulator_array(accumulator_array);
  clear_accumulator_array(accumulator_array2);

  // Hack into vpic internals
  int failed = 0;
  load_interpolator_array( interpolator_array, field_array );
  for( int n=0; n<nstep; n++ ) {

    advance_p( sp, accumulator_array, interpolator_array );
    advance_p2( sp2, accumulator_array2, interpolator_array );

    // This is how many pipelines there are inside the array
    for (int n = 0; n < accumulator_array->n_pipeline+1; n++)
    {
        accumulator_t* a = accumulator_array->a + (n * accumulator_array2->stride);
        accumulator_t* a2 = accumulator_array2->a + (n * accumulator_array2->stride);
        for (int i = 0; i < grid->nv; i++)
        {
            if (
                    (a[i].jx[0] != a2[i].jx[0]) ||
                    (a[i].jx[1] != a2[i].jx[1]) ||
                    (a[i].jx[2] != a2[i].jx[2]) ||
                    (a[i].jx[3] != a2[i].jx[3]) ||
                    (a[i].jy[0] != a2[i].jy[0]) ||
                    (a[i].jy[1] != a2[i].jy[1]) ||
                    (a[i].jy[2] != a2[i].jy[2]) ||
                    (a[i].jy[3] != a2[i].jy[3]) ||
                    (a[i].jz[0] != a2[i].jz[0]) ||
                    (a[i].jz[1] != a2[i].jz[1]) ||
                    (a[i].jz[2] != a2[i].jz[2]) ||
                    (a[i].jz[3] != a2[i].jz[3])
            )
            {
                std::cout << " Failed at " << i << std::endl;
                failed++;
            }
        }
      if( failed )
      {  std::cout << "FAIL" << std::endl;
      }
      REQUIRE_FALSE(failed);
    }

  }

  std::cout << "pass" << std::endl;
}

TEST_CASE( "vectors can be sized and resized", "[vector]" )
{

    int pargc = 0;
    char str[] = "bin/vpic";
    char **pargv = (char **) malloc(sizeof(char **));
    pargv[0] = str;
    boot_services( &pargc, &pargv );

    SECTION( "resizing bigger changes size and capacity" )
    {
        int num_particles = 64;

        vpic_simulation* simulation = new vpic_simulation;
        simulation->initialize( pargc, pargv );

        simulation->finalize();
        delete simulation;
        if( world_rank==0 ) log_printf( "normal exit\n" );

        halt_mp();
    }
}
