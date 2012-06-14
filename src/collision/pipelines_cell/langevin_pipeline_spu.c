#define IN_collision
#define HAS_SPU_PIPELINE
#include <collision_private.h>

#define NEED_frandn
#include <rng_spe.h>

#include <spu_mfcio.h>

// DMA tag usage:
// Tags  0:15 - input stream
// Tags 16:31 - output stream

// Define the buffered I/O streams

#define IN_TAG(b)  (b)
#define OUT_TAG(b) (16+(b))
#define N_BUF 16   /* Must be a power of two */

// FIXME: MOVE TO A SHARED LOCATION
#define OFFSET_OF(t,f) ((size_t)(&(((t *)(NULL))->f)))

void
_SPUEAR_langevin_pipeline_spu( langevin_pipeline_args_t * RESTRICT l,
                               int pipeline_rank,
                               int n_pipeline ) {
  double n_target = (double)l->np / (double)n_pipeline;
  int i = (int)( 0.5 + n_target*(double) pipeline_rank    );
  int n = (int)( 0.5 + n_target*(double)(pipeline_rank+1) ) - i;

  MEM_PTR( particle_t, 128 ) pi_ux =
    l->p + i*sizeof(particle_t) + OFFSET_OF(particle_t,ux);
  rng_t * RESTRICT rng = mfc_get_rng( l->rng[pipeline_rank] );
  const float decay = l->decay;
  const float drive = l->drive;
  float dx, dy, dz;
  int b;

  DECLARE_ALIGNED_ARRAY( float, 16, u_in,  4*N_BUF );
  DECLARE_ALIGNED_ARRAY( float, 16, u_out, 4*N_BUF );
  
  // Fill up the pipeline

  for( b=0; b<N_BUF; b++ )
    if( LIKELY(b<n) ) mfc_get( u_in+4*b, pi_ux + b*sizeof(particle_t),
                               4*sizeof(float), IN_TAG(b), 0, 0 );

  b = 0;
  for( ; n; n-- ) {

    // Compute the random term for this particle

    dx = drive*frandn(rng);
    dy = drive*frandn(rng);
    dz = drive*frandn(rng);

    // Finishing reading in this particle, update it, and start
    // writing it back out

    mfc_write_tag_mask( (1<<IN_TAG(b)) | (1<<OUT_TAG(b)) );
    mfc_read_tag_status_all();

    // FIXME: VECTORIZE 
    u_out[4*b+0] = decay*u_in[4*b+0] + dx;
    u_out[4*b+1] = decay*u_in[4*b+1] + dy;
    u_out[4*b+2] = decay*u_in[4*b+2] + dz;
    u_out[4*b+3] =       u_in[4*b+3];

    mfc_put( u_out+4*b, pi_ux, 4*sizeof(float), OUT_TAG(b), 0, 0 );
    if( LIKELY( N_BUF<n ) ) { 
      mfc_get( u_in+4*b, pi_ux + N_BUF*sizeof(particle_t), 4*sizeof(float),
               IN_TAG(b), 0, 0  );
    } 

    // Advance to the next particle

    pi_ux += sizeof(particle_t); b++; if( b==N_BUF ) b = 0;
  }

  mfc_put_rng( rng, l->rng[pipeline_rank] );
}

