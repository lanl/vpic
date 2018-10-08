#define IN_collision

#include "../util/pipelines/pipelines_exec.h"

/* Private interface *********************************************************/

void
langevin_pipeline_scalar( langevin_pipeline_args_t * RESTRICT args,
                          int pipeline_rank,
                          int n_pipeline )
{
  if ( pipeline_rank == n_pipeline )
  {
    return; /* No host straggler cleanup */
  }

  particle_t * RESTRICT p     = args->p;
  rng_t      * RESTRICT rng   = args->rng[ pipeline_rank ];
  float                 decay = args->decay;
  float                 drive = args->drive;

  double n_target = (double)args->np / (double)n_pipeline;

  /**/  int i  = (int)( 0.5 + n_target * (double)  pipeline_rank    );
  const int i1 = (int)( 0.5 + n_target * (double) (pipeline_rank+1) );

  for( ; i < i1; i++ )
  {
    p[i].ux = decay * p[i].ux + drive * frandn(rng);
    p[i].uy = decay * p[i].uy + drive * frandn(rng);
    p[i].uz = decay * p[i].uz + drive * frandn(rng);
  }
}

#if defined(V4_ACCELERATION) && defined(HAS_V4_PIPELINE)

#error "V4 pipeline not implemented"

#endif

void
apply_langevin_pipeline( langevin_t * l )
{
  if ( l->interval < 1                  ||
       ( l->sp->g->step % l->interval ) )
  {
    return;
  }

  /* Decay and drive have a fun derivation.  We want to integrate the
     stochastic equation:
       du = -nu u dt + sqrt( 2 kT / mc ) dW
     For small dt, this is:
       u_1 = u_0 ( 1- nu dt ) + RANDN( 2 kT nu dt / mc )
     where RANDN( var ) is a normal random number with _variance_ var.
     Let:
       a = nu dt
       b = 2 kT nu dt / mc
     Then:
       u_1 = (1-a) u_0 + RANDN( b )
     We can get more accurate by making N substeps of length dt/N.
     Then:
       u_{n+1} = (1-a/N) u_n + RANDN( b/N )
     such that:
       u_N = (1-a/N)^N u_0 + sum_{n=0:N-1} (1-a/N)^n RANDN( b/N )
     Noting that sum of N normal random numbers is a normal random
     number whose variance is the sum of the N variances, we have:
       u_N = (1-a/N)^N u_0 + RANDN( sum_{n=0:N-1} (1-a/N)^{2n} b/N )
     Analytically summing the variances yields:
       u_N = (1-a/N)^N u_0 + RANDN( [1-(1-a/N)^{2N} / ( 1-(1-a/N)^2 )] b/N )
     In the continuum limit (N goes to infinity):
       u_N = decay u_0 + RANDN( drive^2 )
     or:
       u_N = decay u_0 + drive RANDN( 1 )
     where:
       decay   = lim (1-a/N)^N = exp(-a)
       drive^2 = lim variance sum
               = [(1-exp(-2a) b] / [(1 - 1 + 2a/N) N]
               = ( 1-exp(-2a) ) b / (2a)
     subtituting a and b into decay and drive yields:
       decay   = exp(-nu dt)
       drive   = sqrt( (1-exp(-2 nu dt)) kT / mc )
     In the limit nu dt small:
       decay -> 1 - nu dt
       drive -> sqrt( 2 nu dt kT / mc )
     reproducing the infinitesimal stochastic differential equation.
     In the limit nu dt large:
       decay -> 0
       drive -> sqrt( kT / mc )
     which is equivalent to resampling the momentum with the
     desired temperature. */

  float nudt = l->nu * (float) l->interval * l->sp->g->dt;

  DECLARE_ALIGNED_ARRAY( langevin_pipeline_args_t, 128, args, 1 );

  args->p     = l->sp->p;

  COPY( args->rng, l->rp->rng, N_PIPELINE );

  args->decay = exp( -nudt );
  args->drive = sqrt( ( -expm1( -2 * nudt ) * l->kT ) / ( l->sp->m * l->sp->g->cvac ) );
  args->np    = l->sp->np;

  EXEC_PIPELINES( langevin, args, 0 );

  WAIT_PIPELINES();
}
