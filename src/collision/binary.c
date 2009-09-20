#define IN_collision
/*#define HAS_SPU_PIPLINE*/
/*#define HAS_V4_PIPELINE*/
#include "collision_private.h"

/* FIXME: ADD SAMPLE TO UNARY */

/* Private interface *********************************************************/

typedef struct binary_collision_model {
  char * name;
  binary_rate_constant_func_t rate_constant;
  binary_collision_func_t collision;
  void * params;
  species_t  * spi;
  species_t  * spj;
  rng_pool_t * rp;
  double sample;
  int interval;
  int n_large_pr[ MAX_PIPELINE ];
} binary_collision_model_t;

void
binary_pipeline( binary_collision_model_t * RESTRICT cm,
                 int pipeline_rank,
                 int n_pipeline ) {
  if( pipeline_rank==n_pipeline ) return; /* No host straggler cleanup */

  binary_rate_constant_func_t rate_constant = cm->rate_constant;
  binary_collision_func_t     collision     = cm->collision;

  /**/  void       * RESTRICT params        = cm->params;
  /**/  species_t  * RESTRICT spi           = cm->spi;
  /**/  species_t  * RESTRICT spj           = cm->spj;
  /**/  rng_t      * RESTRICT rng           = cm->rp->rng[ pipeline_rank ];

  /**/  particle_t * RESTRICT spi_p         = spi->p;
  const int        * RESTRICT spi_partition = spi->partition;
  const grid_t     * RESTRICT g             = spi->g;

  /**/  particle_t * RESTRICT spj_p         = spj->p;
  const int        * RESTRICT spj_partition = spj->partition;

  const double sample        = (spi_p==spj_p ? 0.5 : 1)*cm->sample;
  const float  dtinterval_dV = ( g->dt * (float)cm->interval ) / g->dV;

  float pr_norm, pr_coll, wk, wl, w_max, w_min;
  int v, v1, k, k0, nk, rk, l, l0, nl, rl, np, nc, type, n_large_pr = 0;

  /* Stripe the (mostly non-ghost) voxels over threads for load balance */

  v  = VOXEL( 0,0,0,             g->nx,g->ny,g->nz ) + pipeline_rank;
  v1 = VOXEL( g->nx,g->ny,g->nz, g->nx,g->ny,g->nz ) + 1;
  for( ; v<v1; v+=n_pipeline ) {
 
    /* Find the species i computational particles, k, and the species j
       computational particles, l, in this voxel, determine the number
       of computational particle pairs, np and the number of candidate
       pairs, nc, to test for collisions within this voxel.  Also,
       split the range of the fastest integer rng into intervals
       suitable for sampling pairs. */

    k0 = spi_partition[v  ];
    nk = spi_partition[v+1] - k0;
    if( !nk ) continue; /* Nothing to do */
    rk = UINT_MAX / (unsigned)nk;

    if( spi==spj ) {

      /* For intraspecies collisions:
           np = nk(nk+1)/2
         and:
           nc = round( sample nk / 2 )
         such that, for sample==1, on average every particle is tested
         for collision once.  Note that the below pair sampling method
         allows for the possibility of a computational particle
         colliding with itself.  This isn't so silly when considering
         that a computational particle represents many physical
         particles.  At the same time, most microscopic physics
         processes will give a zero collision rate constant for such as
         the implied colliding physical particles are comoving. */

      l0 = k0;
      nl = nk;
      rl = rk;
      np = nk*(nk+1) >> 1;
      nc = (int)( 0.5 + sample*(double)nk );

    } else {

      /* For interspecies collisions:
           np = nk nl
         and:
           nc = round( sample max( nk, nl ) )
         such that, for sample==1, on average every particle is thest
         for collision at least once. */

      l0 = spj_partition[v  ];
      nl = spj_partition[v+1] - l0;
      if( !nl ) continue; /* Nothing to do */
      rl = UINT_MAX / (unsigned)nl;
      np = nk*nl;
      nc = (int)( 0.5 + sample*(double)(nk>nl ? nk : nl) );

    }

    /* Determine the collision rate to probability normalization:
         pr_norm = ( dt interval np ) / ( dV nc ) */

    pr_norm = dtinterval_dV*((float)np / (float)nc);

    /* For each candidate pair */

    for( ; nc; nc-- ) {

      /* Pick a pair of computational particles uniformly at random
         from all pairs of particles in the voxel.  Note that the
         while test virtually always fails (this manner of
         splitting up the nk, nl guarantees a uniform prob 
         of getting k on 0:nk-1 and l on 0:nl-1 and uses the
         preferred high order randgen bits). */
  
      do { k = (int)(uirand(rng)/rk); } while( k==nk ); k += k0;
      do { l = (int)(uirand(rng)/rl); } while( l==nl ); l += l0;
     
      /* Compute the probability that a physical particle in the
         species whose candidate computational particle has the least
         weight (and comoving with this computational particle) will
         collide off a beam of physical particles of the density and
         momentum of the computational particle in the other species.
         If this probability is bigger than one, make a note for
         diagnostic use. */

      wk = spi_p[k].w;
      wl = spj_p[l].w;
      w_max = (wk>wl) ? wk : wl;
      pr_coll = w_max * pr_norm *
        rate_constant( params, spi, spj, &spi_p[k], &spj_p[l] );
      if( pr_coll>1 ) n_large_pr++;

      /* Yes, >= so that 0 rate constants guarantee no collision and
         yes, _c0, so that 1 probabilities guarantee a collision */
      if( frand_c0(rng)>=pr_coll ) continue; /* Didn't collide */

      /* k and l had a collision.  Determine which computational
         particles should be updated by the collision process.
         We should always update the particle of least weight.
         The other particle should be updated with probability of
         w_min / w_max, such that, on average, detailed balance is
         preserved. */
    
      w_min = (wk>wl) ? wl : wk;
      type = 1; if( wl==w_min ) type++;
      if( w_max==w_min || w_max*frand_c0(rng)<w_min ) type = 3;
      collision( params, spi, spj, &spi_p[k], &spj_p[l], rng, type );

    } 
  }

  cm->n_large_pr[pipeline_rank] = n_large_pr;
}

void
apply_binary_collision_model( binary_collision_model_t * cm ) {
  int p, n_large_pr = 0;
  if( cm->interval<1 || (cm->spi->g->step % cm->interval) ) return;
  if( cm->spi->last_sorted!=cm->spi->g->step ) sort_p( cm->spi );
  if( cm->spj->last_sorted!=cm->spi->g->step ) sort_p( cm->spj );
  EXEC_PIPELINES( binary, cm, 0 );
  WAIT_PIPELINES();
  for( p=0; p<N_PIPELINE; p++ ) n_large_pr += cm->n_large_pr[p];
  if( n_large_pr )
    WARNING(( "%i computational particle pairs between species \"%s\" and "
              "species \"%s\" encountered a large collision probability in "
              "collision model \"%s\".  The collision rate for such pairs "
              "will be lower than it should be physically.  Consider lowering "
              "the collision operator interval, increasing the sampling or "
              "reducing the timestep.",
              n_large_pr, cm->spi->name, cm->spj->name, cm->name ));
}

void
checkpt_binary_collision_model( const collision_op_t * cop ) {
  const binary_collision_model_t * cm =
    (const binary_collision_model_t *)cop->params;
  CHECKPT( cm, 1 );
  CHECKPT_STR( cm->name );
  CHECKPT_SYM( cm->rate_constant );
  CHECKPT_SYM( cm->collision );
  CHECKPT_PTR( cm->params );
  CHECKPT_PTR( cm->spi );
  CHECKPT_PTR( cm->spj );
  CHECKPT_PTR( cm->rp );
  checkpt_collision_op_internal( cop );
}

collision_op_t *
restore_binary_collision_model( void ) {
  binary_collision_model_t * cm;
  RESTORE( cm );
  RESTORE_STR( cm->name );
  RESTORE_SYM( cm->rate_constant );
  RESTORE_SYM( cm->collision );
  RESTORE_PTR( cm->params );
  RESTORE_PTR( cm->spi );
  RESTORE_PTR( cm->spj );
  RESTORE_PTR( cm->rp );
  return restore_collision_op_internal( cm );
}

void
delete_binary_collision_model( collision_op_t * cop ) {
  binary_collision_model_t * cm = (binary_collision_model_t *)cop->params;
  FREE( cm->name );
  FREE( cm );
  delete_collision_op_internal( cop );
}

/* Public interface **********************************************************/

collision_op_t *
binary_collision_model( const char       * RESTRICT name,
                        binary_rate_constant_func_t rate_constant,
                        binary_collision_func_t     collision,
                        /**/  void       * RESTRICT params,
                        /**/  species_t  * RESTRICT spi,
                        /**/  species_t  * RESTRICT spj,
                        /**/  rng_pool_t * RESTRICT rp,
                        double                      sample,
                        int                         interval ) {
  binary_collision_model_t * cm;
  size_t len = name ? strlen(name) : 0;

  if( !rate_constant || !collision || !spi || !spj || spi->g!=spj->g ||
      !rp || rp->n_rng<N_PIPELINE ) ERROR(( "Bad args" ));
  if( len==0 ) ERROR(( "Cannot specify a nameless collision model" ));
  if( params && !object_id( params ) )
    ERROR(( "collision model parameters must be checkpoint registered" ));

  MALLOC( cm, 1 );
  MALLOC( cm->name, len+1 );
  strcpy( cm->name, name ); 
  cm->rate_constant = rate_constant;
  cm->collision     = collision;
  cm->params        = params;
  cm->spi           = spi;
  cm->spj           = spj;
  cm->rp            = rp;
  cm->sample        = sample;
  cm->interval      = interval;
  return new_collision_op_internal( cm,
                                    (collision_op_func_t)apply_binary_collision_model,
                                    delete_binary_collision_model,
                                    (checkpt_func_t)checkpt_binary_collision_model,
                                    (restore_func_t)restore_binary_collision_model,
                                    NULL );
}

