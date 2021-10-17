#define IN_collision

/* #define HAS_V4_PIPELINE */

#include "collision_pipeline.h"

#include "../chemical.h"

#include "../../util/pipelines/pipelines_exec.h"

unsigned long long
multichoose(const int n, const int k)
{
  // multichoose(n,k) = choose(n+k-1,k) = prod_i=1^k (n+i-1)/i
  unsigned long long y=1, i=1, j=n;
  for( ; i <= k ; ++i, ++j) {
    if      ( j%i == 0 ) y *= j/i;
    else if ( y%i == 0 ) y  = (y/i)*j;
    else                 y  = (y*j)/i;
  }
  return y;
}

int
compare_particle_mover(const void *a, const void *b)
{
  return ((particle_mover_t*)a)->i - ((particle_mover_t*)b)->i ;
}

void
chemical_pipeline_scalar( chemical_collision_model_t * RESTRICT cm,
                          int pipeline_rank,
                          int n_pipeline )
{
  if( pipeline_rank==n_pipeline ) return; /* No host straggler cleanup */

  chemical_rate_constant_func_t rate_constant = cm->rate_constant;
  chemical_collision_func_t     collision     = cm->collision;

  /**/  void        * RESTRICT params     = cm->params;
  /**/  species_t  ** RESTRICT reactants  = cm->reactants;
  /**/  species_t  ** RESTRICT products   = cm->products;
  const int         * RESTRICT consumable = cm->consumable;
  /**/  rng_t       * RESTRICT rng        = cm->rp->rng[ pipeline_rank ];
  const grid_t      * RESTRICT g          = cm->reactants[0]->g;

  const int    N             = cm->n_reactants;
  const int    M             = cm->n_products;
  const double sample        = cm->sample;
  const float  dtinterval_dV = ( g->dt * (float)cm->interval ) / g->dV;

  species_t * sp;
  float pr_norm, pr_coll, w_prod, w_min, pr_max = 0;
  int n, m, k, v, v1, nc, type, n_large_pr = 0, n_tested = 0, unaccumulated = 0;
  unsigned long long np;

  // Reactant variables
  //  i     : Particle index for reactant i
  //  i0    : Start of particles for reactant i in this voxel
  //  ni    : Number of particles of reactant i in this voxel
  //  ri    : Scale from random integer to particle index
  //  im    : Index of next available reactant mover
  //  nm    : Number of available reactant movers
  //  wi    : Statistical weight of reactant particle i
  //  mass  : Species mass of reactant i
  //  mpart : Temporary variable for particle mass (weight * species mass)
  //  mtot  : Temporary variable for sum(mpart)
  //  com   : center of mass coordinates
  //  n_consumed : Number of reactant i consumed
  //  reactant_particles : Pointers to the reactant particles.
  int i[N], i0[N], ni[N], ri[N], im[N], nm[N], n_modified[N];
  float wi[N], mass[N], com[3], mpart, mtot;
  particle_t *reactant_particles[N];

  // Product variables
  //  j  : Index of the next available particle for product j
  //  j0 : Index of first available particle for product j
  //  nj : How much free space is available for product j
  //  product_particles : Pointers to the product particles
  int j[M], j0[M], nj[M], n_produced=0;
  particle_t *product_particles[M];

  // Determine the free space available for this pipeline to inject product
  // particles.
  FOR_PRODUCTS {
    DISTRIBUTE(sp->max_np-sp->np, 16, pipeline_rank, n_pipeline, j0[m], nj[m]);
    j0[m] += sp->np;
    j[m]   = j0[m];
  }

  // Determine which particle movers belong to this pipeline. We are using
  // particle movers to temporarially hold deleted particles since they are
  // already allocated and sized appropriately.
  FOR_REACTANTS {
    DISTRIBUTE(sp->max_nm, 8, pipeline_rank, n_pipeline, im[n], nm[n]);
    n_modified[n] = 0;
    mass[n] = sp->m;
  }

  // Stripe the (mostly non-ghost) voxels over threads for load balance.
  v  = VOXEL( 0,0,0,             g->nx,g->ny,g->nz ) + pipeline_rank;
  v1 = VOXEL( g->nx,g->ny,g->nz, g->nx,g->ny,g->nz ) + 1;
  for( ; v<v1; v+=n_pipeline ) {

    /* Find the species i computational particles, k, and the species j
       computational particles, l, in this voxel, determine the number
       of computational particle pairs, np and the number of candidate
       pairs, nc, to test for collisions within this voxel.  Also,
       split the range of the fastest integer rng into intervals
       suitable for sampling pairs. */

    np = 1;
    CLEAR(ni, N);
    FOR_REACTANTS if( !ni[n] ) {

        // Compute number of particles for each unique reactant species.
        i0[n] = sp->partition[v  ];
        ni[n] = sp->partition[v+1] - i0[n];
        if( !ni[n] ) goto next_voxel;
        ri[n] = UINT_MAX / (unsigned)ni[n];

        // Check to see if this species is repeated in the reactants, and
        // if so, count the number of repetitions and update arrays.
        for( k=1, m=n+1 ; m<N ; ++m )
          if( sp == reactants[m] ){
            ++k;
            i0[m] = i0[n];
            ni[m] = ni[n];
            ri[m] = ri[n];
          }

        // Number of combinations = prod (ni multichoose k)
        np *= multichoose(ni[n],k);

      }

    // Determine how many combinations to sample and probability normalization.
    nc        = (int)(1 + sample*(double)np);
    pr_norm   = dtinterval_dV * ((float)np / (float)nc);
    n_tested += nc;

    // Sample the collision integral.
    for( ; nc; nc-- ) {

      // Pick computational particles uniformly at random from all  particles in
      // the voxel.  Note that the while test virtually always fails (this
      // manner of splitting up the nk, nl guarantees a uniform prob of getting
      // k on 0:nk-1 and l on 0:nl-1 and uses the preferred high order randgen
      // bits).

      w_prod = 1;
      FOR_REACTANTS {
        do { i[n] = (int)(uirand(rng)/ri[n]); } while( i[n]==ni[n] );
        i[n] += i0[n];
        reactant_particles[n] = sp->p + i[n];

        // Compute prod(w) and find min(w)
        wi[n]   = reactant_particles[n]->w;
        w_prod *= wi[n];
        if( n == 0 || wi[n] < w_min ) w_min = wi[n];
      }

      // There is a possibility that if a particle is picked for collisions
      // that it's weight could be <= 0. In this case, we should really just
      // pick another set of particles in order to maintain the correct
      // probabilities, however this could lead to an infinite loop if _all_
      // particles are consumed. The correct way to do this would be to
      // recompute probabilities after each collision, taking into account
      // particle consumption, but this is a lot of overhead. For now, we
      // will just skip this collision and take the loss in accuracy.
      if( w_min <= 0 ) continue;


      /* Compute the probability that a physical particle in the
         species whose candidate computational particle has the least
         weight (and comoving with this computational particle) will
         collide off a beam of physical particles of the density and
         momentum of the computational particle in the other species.
         If this probability is bigger than one, make a note for
         diagnostic use. */

      pr_coll = w_prod / w_min * pr_norm *
                rate_constant( params, reactants, reactant_particles );

      if( pr_coll >      1 ) n_large_pr++;
      if( pr_coll > pr_max ) pr_max = pr_coll;

      /* Yes, >= so that 0 rate constants guarantee no collision and
         yes, _c0, so that 1 probabilities guarantee a collision */
      if( frand_c0(rng)>=pr_coll ) continue; /* Didn't collide */

      // Reactants had a collision.  Determine which computational particles
      // should be updated by the collision process. We should always update the
      // particle of least weight. The other particles should be updated with
      // probability of w_min / w[i] such that, on average, detailed balance is
      // preserved.

      type = 0;
      com[0] = com[1] = com[2] = 0;
      FOR_REACTANTS{
        if( wi[n] == w_min || wi[n]*frand_c0(rng)<w_min) type |= (1 << n);
        mpart   = wi[n] * mass[n];
        com[0] += mpart * reactant_particles[n]->dx;
        com[1] += mpart * reactant_particles[n]->dy;
        com[2] += mpart * reactant_particles[n]->dz;
        mtot   += mpart;
      }
      com[0] /= mtot;
      com[1] /= mtot;
      com[2] /= mtot;

      // Now create the product particles. They should all be created with
      // a statistical weight of w_min and located at the center of mass.
      // NOTE: If cells are concave, then COM may be outside of the cell!

      FOR_PRODUCTS {
        product_particles[m]     = sp->p + j[m];
        product_particles[m]->i  = v;
        product_particles[m]->dx = com[0];
        product_particles[m]->dy = com[1];
        product_particles[m]->dz = com[2];
        product_particles[m]->w  = w_min;
        ++j[m];
      }

      // Do the collision. The collision function should set u{x,y,z} and d{x,y,z}
      // for the product particles, and may also modify u{x,y,z} and w for the
      // reactant particles.
      collision( params, reactants, reactant_particles,
                 products, product_particles, rng, type );

      // Check for reactant consumption which is flagged by a change in particle
      // statistical weight post-collison. These modified particles (along with
      // product particles) need to undergo rhob accumulation, however doing this
      // from within the pipeline is significant effort since pipelines would want
      // to write to the same location. To do this, we would either need to use
      // blocking constructs or have each pipeline allocate a field array and then
      // reduce across pipelines. It is much easier to offload this work to the host,
      // although this raises an issue if more particles are modified than exist
      // particle movers. This should not be a common event, and can be avoided by
      // either allocating more movers or decreasing the collisional timestep.
      FOR_REACTANTS if( consumable[n] ) {
          reactant_particles[n]->w -= w_min;
          if( nm[n] <= n_modified[n] ) {
            ++unaccumulated;
          }
          else {
            sp->pm[im[n]].i     = i[n];
            sp->pm[im[n]].dispx = w_min;
            ++im[n];
            ++n_modified[n];
          }
        }

      // Check that we have enough product storage left to continue to process
      // collisions. The use of goto is not ideal, but cleaner than 3x break.
      FOR_PRODUCTS if( j[m] - j0[m] >= nj[m] ) goto final;

    }

    next_voxel:;

  }

  final:

    // Save pipeline-produced information to the model structure.
    FOR_REACTANTS {
      cm->n_modified[N*pipeline_rank + n] = n_modified[n];
      cm->reactant_movers[N*pipeline_rank + n] = sp->pm + im[n] - n_modified[n];
    }

    cm->n_produced[pipeline_rank] = j[0]-j0[0];
    FOR_PRODUCTS cm->product_particles[M*pipeline_rank + m] = sp->p + j0[m];

    cm->unaccumulated[pipeline_rank] = unaccumulated;
    cm->n_large_pr[pipeline_rank]    = n_large_pr;
    cm->pr_max[pipeline_rank]        = pr_max;
    cm->n_tested[pipeline_rank]      = n_tested;
}

void
apply_chemical_collision_model_pipeline( chemical_collision_model_t * cm )
{
  species_t *sp;
  double pr_max = 0;
  int n, m, p, np, nm, n_large_pr = 0, n_tested = 0, unaccumulated = 0;
  particle_t *p0;
  particle_mover_t *m0;
  float dw, w0;
  /**/  species_t  ** RESTRICT reactants = cm->reactants;
  /**/  species_t  ** RESTRICT products  = cm->products;
  const grid_t      * RESTRICT g         = reactants[0]->g;
  /**/  field_t     * RESTRICT f         = cm->fa->f;
  const int N = cm->n_reactants;
  const int M = cm->n_products;

  if( cm->interval<1 || (g->step % cm->interval) ) return;
  FOR_REACTANTS if( sp->last_sorted != g->step ) sort_p( sp );
  EXEC_PIPELINES( chemical, cm, 0 );
  WAIT_PIPELINES();

  // Combine elements from all pipelines and do some reductions.
  for( p=0; p<N_PIPELINE; p++ ) {

    // Record excessive probabilities.
    n_large_pr    += cm->n_large_pr[p];
    n_tested      += cm->n_tested[p];
    unaccumulated += cm->unaccumulated[p];
    if( cm->pr_max[p] > pr_max ) pr_max = cm->pr_max[p];

    // Eliminate holes in product particle arrays and accumulate bound charge.
    FOR_PRODUCTS {
      np = cm->n_produced[p];
      p0 = cm->product_particles[p*M+m];
      if( p0 != sp->p + sp->np ) MOVE(sp->p + sp->np, p0, np);
      for( ; np ; --np, ++(sp->np) ) accumulate_rhob(f, sp->p+sp->np, g, -sp->q);
    }

    // Eliminate holes in reactant particle movers.
    FOR_REACTANTS {
      nm = cm->n_modified[p*N+n];
      m0 = cm->reactant_movers[p*N+n];
      if( m0 != sp->pm + sp->nm ) MOVE(sp->pm + sp->nm, m0, nm);
      sp->nm += nm;
    }

  }

  // Since we striped voxels for better load balancing, the reactant particle
  // movers are out of order. We need to sort them before we can backfill. This
  // could be optimized significantly, but we expect len(pm) << len(p) so
  // complex sorting functions are more effort than they are worth here.
  FOR_REACTANTS {
    nm = 0; // Number of particles deleted from this species
    qsort(sp->pm, sp->nm, sizeof(particle_mover_t), compare_particle_mover);
    for( ; sp->nm ; --(sp->nm) ){
      p0 = sp->p + sp->pm[sp->nm-1].i;

      // Fractional accumulation of rhob. Charge = dw q
      dw = sp->pm[sp->nm-1].dispx;
      w0 = p0->w;
      p0->w = dw;
      accumulate_rhob(f, p0, g, sp->q);
      p0->w = w0;

      // Check for complete consumption and backfill.
      if(w0 <= 0) {
        *p0 = sp->p[ --(sp->np) ];
        ++nm;
      }

    }
    if( nm ) sort_p(sp);
  }

  if( n_large_pr ) {
    WARNING(( "%i/%i computational particle groupings encountered a large "
              "collision probability in collision model \"%s\".  The "
              "collision rate for such groupings will be lower than it should "
              "be physically.  Consider lowering the collision operator "
              "interval, increasing the sampling or reducing the timestep.",
              n_large_pr, n_tested, cm->name ));
    if( cm->dynamic_sampling ) cm->sample *= 2*pr_max;
  }

  if( unaccumulated )
    WARNING(( "%i reactant particles were consumed without proper rhob "
              "accumulation in collision model \"%s\". Local charge "
              "conservation may be broken.", unaccumulated, cm->name ));

}
