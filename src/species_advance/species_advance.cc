/*
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Original version
 *
 */

#include "species_advance.h"

/* Private interface *********************************************************/

void
checkpt_species( const species_t * sp )
{
  CHECKPT( sp, 1 );
  CHECKPT_STR( sp->name );
  checkpt_data( sp->p,
                sp->np     * sizeof(particle_t),
                sp->max_np * sizeof(particle_t), 1, 1, 128 );
  checkpt_data( sp->pm,
                sp->nm     * sizeof(particle_mover_t),
                sp->max_nm * sizeof(particle_mover_t), 1, 1, 128 );
  CHECKPT_ALIGNED( sp->partition, sp->g->nv+1, 128 );
  CHECKPT_PTR( sp->g );
  CHECKPT_PTR( sp->next );

  #ifdef VPIC_GLOBAL_PARTICLE_ID
  // RFB: Add checkpointing of particle list
  if(sp->has_ids) {
    checkpt_data( sp->p_id,
                  sp->np     * sizeof(size_t),
                  sp->max_np * sizeof(size_t), 1, 1, 128 );
  }
  #endif

}

species_t *
restore_species( void )
{
  species_t * sp;
  RESTORE( sp );
  RESTORE_STR( sp->name );
  sp->p  = (particle_t *)        restore_data();
  sp->pm = (particle_mover_t *)  restore_data();
  RESTORE_ALIGNED( sp->partition );
  RESTORE_PTR( sp->g );
  RESTORE_PTR( sp->next );
  #ifdef VPIC_GLOBAL_PARTICLE_ID
  if(sp->has_ids) {
    sp->p_id  = (size_t*)        restore_data();
  } else {
    sp->p_id  = (size_t*)        nullptr;
  }
  #endif
  return sp;
}

void
delete_species( species_t * sp )
{
  UNREGISTER_OBJECT( sp );
  FREE_ALIGNED( sp->partition );
  FREE_ALIGNED( sp->pm );
  FREE_ALIGNED( sp->p );
  #ifdef VPIC_GLOBAL_PARTICLE_ID
  if(sp->has_ids) {
    FREE_ALIGNED( sp->p_id );
  }
  #endif
  FREE( sp->name );
  FREE( sp );
}

#ifdef VPIC_GLOBAL_PARTICLE_ID
  std::string make_tracer_name_unique(const std::string prefix, species_t* sp_list) {
    // check if a species with that name exists in that species_list
    species_t * sp = find_species_name(prefix.c_str(), sp_list);
    if(!sp) { // No species with that name found
      return prefix; // We can just use that name
    } else {
     // We have to append things to try and make the name unique
     int postfix = 2;
     std::string name = prefix;
     while(sp) {
       name = prefix + std::to_string(postfix);
       sp = find_species_name(name.c_str(), sp_list);
     }
     return name;
    }
  }
#endif

/* Public interface **********************************************************/

int
num_species( const species_t * sp_list )
{
  return sp_list ? sp_list->id+1 : 0;
}

void
delete_species_list( species_t * sp_list )
{
  species_t * sp;
  while( sp_list ) {
    sp = sp_list;
    sp_list = sp_list->next;
    delete_species( sp );
  }
}

species_t *
find_species_id( species_id id,
                 species_t * sp_list )
{
  species_t * sp;
  LIST_FIND_FIRST( sp, sp_list, sp->id==id );
  return sp;
}

species_t *
find_species_name( const char * name,
                   species_t * sp_list )
{
  species_t * sp;
  if( !name ) return NULL;
  LIST_FIND_FIRST( sp, sp_list, strcmp( sp->name, name )==0 );
  return sp;
}

species_t *
append_species( species_t * sp,
                species_t ** sp_list )
{
  if( !sp || !sp_list ) ERROR(( "Bad args" ));
  if( sp->next ) ERROR(( "Species \"%s\" already in a list", sp->name ));
  if( find_species_name( sp->name, *sp_list ) )
    ERROR(( "There is already a species in the list named \"%s\"", sp->name ));
  if( (*sp_list) && sp->g!=(*sp_list)->g )
    ERROR(( "Species \"%s\" uses a different grid from this list", sp->name ));
  sp->id   = num_species( *sp_list );
  sp->next = *sp_list;
  *sp_list = sp;
  return sp;
}

species_t *
species( const char * name,
         float q,
         float m,
         size_t max_local_np,
         size_t max_local_nm,
         int sort_interval,
         int sort_out_of_place,
         grid_t * g
         )
{
  species_t * sp;
  int len = name ? strlen(name) : 0;

  if( !len ) ERROR(( "Cannot create a nameless species" ));
  if( !g ) ERROR(( "NULL grid" ));
  if( g->nv == 0) ERROR(( "Allocate grid before defining species." ));
  if( max_local_np<1 ) max_local_np = 1;
  if( max_local_nm<1 ) max_local_nm = 1;

  MALLOC( sp, 1 );
  CLEAR( sp, 1 );

  MALLOC( sp->name, len+1 );
  strcpy( sp->name, name );

  sp->q = q;
  sp->m = m;

  MALLOC_ALIGNED( sp->p, max_local_np, 128 );
  #ifdef VPIC_GLOBAL_PARTICLE_ID
  sp->has_ids = 0;    // By default a newly created species will not have IDs.
  sp->p_id = nullptr; // if we ever dereference this nullptr, I screwed up
  //  MALLOC_ALIGNED( sp->p_id, max_local_np, 128 );
  #endif
  sp->max_np = max_local_np;

  MALLOC_ALIGNED( sp->pm, max_local_nm, 128 );
  sp->max_nm = max_local_nm;

  sp->last_sorted       = INT64_MIN;
  sp->sort_interval     = sort_interval;
  sp->sort_out_of_place = sort_out_of_place;
  MALLOC_ALIGNED( sp->partition, g->nv+1, 128 );

  sp->g = g;

  /* id, next are set by append species */

  REGISTER_OBJECT( sp, checkpt_species, restore_species, NULL );
  return sp;
}

species_t * tracerspecies_by_percentage(const species_t* parentspecies,
                                        const float percentage,
                                        std::string name,
                                        species_t* sp_list,
                                        grid_t* grid) {

  // REVIEW change the provided name if need be and surprise the user, or fail loudly?
  //std::string name = make_tracer_name_unique(tracername, sp_list);
  if(find_species_name(name.c_str(), sp_list)) {
    ERROR(( "Species with name %d already exists", name.c_str() ));
  }
  const float q = parentspecies->q;
  const float m = parentspecies->m;
  const size_t max_local_np = ceil(percentage/100.0 * parentspecies->max_np) + 1;
  const size_t max_local_nm = ceil(percentage/100.0 * parentspecies->max_nm) + 1;
  const int sort_interval = parentspecies->sort_interval;
  const int sort_out_of_place = parentspecies->sort_out_of_place;

  species_t* tracerspecies = species(name.c_str(), q, m, max_local_np, max_local_nm, sort_interval, sort_out_of_place, grid);
  if(!tracerspecies) ERROR(( "Creation of tracerspecies failed" ));
  // Grab into the species and make it have IDs
  tracerspecies->has_ids = 1;
  MALLOC_ALIGNED( tracerspecies->p_id, max_local_np, 128 );

  // Move the desired percentage of particles from the parent species to the tracer species
  // REVIEW: Should that be copy instead of move?
}
