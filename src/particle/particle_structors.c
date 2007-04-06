/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Original version
 *
 */

#include <particle.h>

particle_t * ALIGNED new_particle_array( int np ) {
  particle_t * ALIGNED p;

  if( np<1 ) ERROR(("Bad np"));

  p = (particle_t * ALIGNED)
    malloc_aligned( np*sizeof(particle_t), preferred_alignment );
  if( p==NULL ) ERROR(("Failed to allocate particle array."));
  memset( p, 0, np*sizeof(particle_t) );

  return p;
}

void delete_particle_array( particle_t ** ALIGNED p ) {
  if( p==NULL ) return;
  if( *p!=NULL ) free_aligned(*p);
  *p = NULL;
}

particle_mover_t * ALIGNED new_particle_mover( int nm ) {
  particle_mover_t * ALIGNED pm;

  if( nm<1 ) ERROR(("Bad nm"));

  pm = (particle_mover_t * ALIGNED)
    malloc_aligned( nm*sizeof(particle_mover_t), preferred_alignment );
  if( pm==NULL ) ERROR(("Failed to allocate particle mover."));
  memset( pm, 0, nm*sizeof(particle_mover_t) );

  return pm;
}

void delete_particle_mover( particle_mover_t ** ALIGNED pm ) {
  if( pm==NULL ) return;
  if( *pm!=NULL ) free_aligned(*pm);
  *pm = NULL;
}
