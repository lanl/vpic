#define IN_particle_pipeline
#include <particle_pipelines.h>
#include <libspe2.h>

extern spe_program_handle_t advance_p_pipeline_spu;

int
main( int argc,
      char ** argv ) {
  int NP, NM, NS;

  DECLARE_ALIGNED_ARRAY( advance_p_pipeline_args_t, 128, args, 1 );
  DECLARE_ALIGNED_ARRAY( particle_mover_seg_t, 128, seg, MAX_PIPELINE+1 );

  if( argc>1 ) NP = atoi(argv[1]);
  else         NP = 16384;

  NM = NP;

  if( argc>2 ) NS = atoi(argv[2]);
  else         NS = 8;

  MESSAGE(( "Booting" ));
  
  spu.boot( NS, 0 );

  MESSAGE(( "Setting args" ));

# define NX 1
# define NY 1
# define NZ 1
# define NV (NX+2)*(NY+2)*(NZ+2)

  MESSAGE(( "NP = %i, NM = %i, NS = %i, MESH = (%i,%i,%i)", NP, NM, NS, NX, NY, NZ ));
  MESSAGE(( "particle_t       = %i", (int)sizeof(particle_t)       ));
  MESSAGE(( "particle_mover_t = %i", (int)sizeof(particle_mover_t) ));
  MESSAGE(( "accumulator_t    = %i", (int)sizeof(accumulator_t)    ));
  MESSAGE(( "interpolator_t   = %i", (int)sizeof(interpolator_t)   ));

  args->p0       = malloc_aligned( NP*sizeof(particle_t),       128 );
  args->pm       = malloc_aligned( NM*sizeof(particle_mover_t), 128 );
  args->a0       = malloc_aligned( NV*sizeof(accumulator_t),    128 );
  args->f0       = malloc_aligned( NV*sizeof(interpolator_t),   128 );
  args->seg      = seg;
  args->neighbor = malloc_aligned( NV*6*sizeof(int64_t),        128 );
  args->rangel   = 0;
  args->rangeh   = NV-1;
  args->qdt_2mc  = 0;
  args->cdt_dx   = 0;
  args->cdt_dy   = 0;
  args->cdt_dz   = 0;
  args->np       = NP;
  args->max_nm   = NM;
  args->nx       = NX;
  args->ny       = NY;
  args->nz       = NZ;

  memset( args->p0, 0, NP*sizeof(args->p0[0]) );
  memset( args->pm, 0, NM*sizeof(args->pm[0]) );
  memset( args->a0, 0, (1+spu.n_pipeline)*NV*sizeof(args->a0[0])  );
  memset( (interpolator_t *)args->f0, 0, NV*sizeof(args->f0[0]) );
  memset( (int64_t *)args->neighbor, 0, NV*6*sizeof(int64_t) );

  MESSAGE(( "Dispatching" ));

  spu.dispatch( SPU_PIPELINE(advance_p_pipeline_spu), args, 0 );

  MESSAGE(( "Waiting" ));

  spu.wait();

  MESSAGE(( "Halting" ));

  spu.halt();

  MESSAGE(( "Finished" ));

  return 0;
}

