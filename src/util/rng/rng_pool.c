#include <rng.h>
#include <checkpt.h>

/* Private API ***************************************************************/

void
checkpt_rng_pool ( const rng_pool_t * rp ) {
  int n;
  CHECKPT( rp, 1 );
  CHECKPT( rp->rng, rp->n_rng );
  for( n=0; n<rp->n_rng; n++ ) CHECKPT_PTR( rp->rng[n] );
}

rng_pool_t *
restore_rng_pool( void ) {
  rng_pool_t * rp;
  int n;
  RESTORE( rp );
  RESTORE( rp->rng );
  for( n=0; n<rp->n_rng; n++ ) RESTORE_PTR( rp->rng[n] );
  return rp;
}

/* Public API ****************************************************************/

rng_pool_t *
new_rng_pool( int n_rng,
              int seed,
              int sync ) {
  rng_pool_t * rp;
  int n;
  if( n_rng<1 ) ERROR(( "Bad args" ));
  MALLOC( rp, 1 );
  MALLOC( rp->rng, n_rng );
  for( n=0; n<n_rng; n++ ) rp->rng[n] = new_rng( 0 );
  rp->n_rng = n_rng;
  seed_rng_pool( rp, seed, sync );
  REGISTER_OBJECT( rp, checkpt_rng_pool, restore_rng_pool, NULL );
  return rp;
}

void
delete_rng_pool( rng_pool_t * rp ) {
  int n;
  if( !rp ) return;
  UNREGISTER_OBJECT( rp );
  for( n=0; n<rp->n_rng; n++ ) delete_rng( rp->rng[n] );
  FREE( rp->rng );
  FREE( rp );
}

rng_pool_t *
seed_rng_pool( rng_pool_t * RESTRICT rp,
               int seed,
               int sync ) {
  int n;
  if( !rp ) ERROR(( "Bad args" ));
  seed = (sync ? world_size : world_rank) + (world_size+1)*rp->n_rng*seed;
  for( n=0; n<rp->n_rng; n++ ) seed_rng( rp->rng[n], seed + (world_size+1)*n );
  return rp;
}

