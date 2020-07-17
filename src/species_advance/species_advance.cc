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
  #ifdef VPIC_PARTICLE_ANNOTATION
  if(sp->has_annotation) {
    checkpt_data( sp->p_annotation,
                  sp->np     * sp->has_annotation * sizeof(float),
                  sp->max_np * sp->has_annotation * sizeof(float), 1, 1, 128);
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
  #ifdef VPIC_PARTICLE_ANNOTATION
  if(sp->has_annotation) {
    sp->p_annotation  = (float*) restore_data();
  } else {
    sp->p_annotation  = (float*) nullptr;
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
  #ifdef VPIC_PARTICLE_ANNOTATION
  if(sp->has_annotation) {
    FREE_ALIGNED( sp->p_annotation );
  }
  #endif
  FREE( sp->name );
  FREE( sp );
}

/**
 * @brief Modify a species name to make sure if does not exist yet in a species list
 *
 * @param prefix This is the name suggestion to start with. The final name will start with this string
 * @param sp_list The species list to check for name collisions
 *
 * @return A name that starts with prefix and has a suitable (possibly empty) suffix append to make sure it is unqiue with the species list
 */
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
     postfix++;
  }
  return name;
 }
}

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
  #endif
  #ifdef VPIC_PARTICLE_ANNOTATION
  sp->has_annotation = 0;    // By default a newly created species will not have annotation buffers
  sp->p_annotation = nullptr; // if we ever dereference this nullptr, I screwed up
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

/**
 * @brief Create a separate species by copying/moving some particles from a parent species
 *
 * @note These tracerspecies functions might still get renamed before we hand them to users
 *
 * @param parentspecies The species from which we source particles. When using "move" it will be modified.
 * @param skip Grab on particle, skip the next skip-1. Can be non-integer.
 * @param copyormove An enum that can be Tracertype::Copy in which case the particle stays in the parrent species and a zero-weight copy is added to the new species or Tracertype::Move in which case the particle is removed from the parent species (and keeps it's statistical weight)
 * @param name The name for the newly created species
 * @param sp_list The list of species we intend to add the newly created species to. Allow to check that param nbame will not clash with any existing species in that list
 * @param grid The global simulation grid. A reference to it will be stored inside the newly created tracer species
 *
 * @return The newly created tracer species
 */
species_t * tracerspecies_by_skip(species_t* parentspecies,
                                  const float skip,
                                  const Tracertype copyormove,
                                  std::string name,
                                  species_t* sp_list,
                                  grid_t* grid) {

  if(find_species_name(name.c_str(), sp_list)) {
    ERROR(( "Species with name %d already exists", name.c_str() ));
  }
  const float q = parentspecies->q;
  const float m = parentspecies->m;
  const size_t max_local_np = ceil(parentspecies->max_np/skip) + 1;
  const size_t max_local_nm = ceil(parentspecies->max_nm/skip) + 1;
  const int sort_interval = parentspecies->sort_interval;
  const int sort_out_of_place = parentspecies->sort_out_of_place;

  species_t* tracerspecies = species(name.c_str(), q, m, max_local_np, max_local_nm, sort_interval, sort_out_of_place, grid);
  if(!tracerspecies) ERROR(( "Creation of tracerspecies failed" ));

  // If we do compile without global_particle_IDs the resulting species will
  // not actually be a good tracer species. But this function might be useful
  // to peel of a fration of particles into a new species for other uses.
  #ifdef VPIC_GLOBAL_PARTICLE_ID
    // Grab into the species and make it have IDs
    tracerspecies->has_ids = 1;
    MALLOC_ALIGNED( tracerspecies->p_id, max_local_np, 128 );
  #endif

  // If the parentspecies has annotations we should have themon the tracers as well
  #ifdef VPIC_PARTICLE_ANNOTATION
  if(parentspecies->has_annotation){
    tracerspecies->allocate_annotation_buffer(parentspecies->has_annotation);
  }
  #endif

  // Select the desired fraction of particles from the parent species and add to the tracer species
  int step = 0;
  const int np = parentspecies->np;
  for(int i = np-1; i>= 0; i--) {
    if(np-i >= (step+1) * skip) {
      // Copy that particle over
      tracerspecies->p[step] = parentspecies->p[i];
      tracerspecies->np++;
      #ifdef VPIC_GLOBAL_PARTICLE_ID
        // Create an ID
        tracerspecies->p_id[step] = tracerspecies->generate_particle_id(step, tracerspecies->max_np);
      #endif

      #ifdef VPIC_PARTICLE_ANNOTATION
      // Copy over annotation
      if(parentspecies->has_annotation){
        for(int j=0; j<parentspecies->has_annotation; j++) {
          const float v = parentspecies->get_annotation(i,j);
          tracerspecies->set_annotation(step, j, v);
        }
      }
      #endif

      if(copyormove == Tracertype::Move) {
        // Remove from parent species
        parentspecies->p[i] = parentspecies->p[parentspecies->np-1];
        parentspecies->np--;
      } else if (copyormove == Tracertype::Copy) {
        // Copied tracers should have zero statistical weight
        tracerspecies->p[step].w = 0.;
      } else {
        ERROR(( "Invalid enum value for copyormove" ));
      }
      // Increment step
      step++;
    }
  }

  return tracerspecies;
}

/**
 * @brief Create a separate species by copying/moving some particles from a parent species
 *
 * @note These tracerspecies functions might still get renamed before we hand them to users
 * @note A predicate(particle_t -> bool) might map badly to GPU. maybe offer predicate(dx,dy.dz,i,ux,uy,uz,w,id->bool there)
 *
 * @param parentspecies The species from which we source particles. When using "move" it will be modified.
 * @param f The user supplied predicate function that takes a particle_t and returns true if the particle should be used in the tracer species
 * @param copyormove An enum that can be Tracertype::Copy in which case the particle stays in the parrent species and a zero-weight copy is added to the new species or Tracertype::Move in which case the particle is removed from the parent species (and keeps it's statistical weight)
 * @param name The name for the newly created species
 * @param sp_list The list of species we intend to add the newly created species to. Allow to check that param nbame will not clash with any existing species in that list
 * @param grid The global simulation grid. A reference to it will be stored inside the newly created tracer species
 *
 * @return The newly created tracer species
 */
species_t * tracerspecies_by_predicate(species_t* parentspecies,
                                       std::function <bool (particle_t)> f,
                                       const Tracertype copyormove,
                                       std::string name,
                                       species_t* sp_list,
                                       grid_t* grid) {

  if(find_species_name(name.c_str(), sp_list)) {
    ERROR(( "Species with name %d already exists", name.c_str() ));
  }
  const float q = parentspecies->q;
  const float m = parentspecies->m;
  const size_t count_true = std::count_if( parentspecies->p, parentspecies->p + parentspecies->np, f);
  const size_t max_local_np = ceil(parentspecies->max_np * count_true/float(parentspecies->np)) + 1;
  const size_t max_local_nm = ceil(parentspecies->max_nm * count_true/float(parentspecies->np)) + 1;
  const int sort_interval = parentspecies->sort_interval;
  const int sort_out_of_place = parentspecies->sort_out_of_place;

  species_t* tracerspecies = species(name.c_str(), q, m, max_local_np, max_local_nm, sort_interval, sort_out_of_place, grid);
  if(!tracerspecies) ERROR(( "Creation of tracerspecies failed" ));

  // If we do compile without global_particle_IDs the resulting species will
  // not actually be a good tracer species. But this function might be useful
  // to peel of a fration of particles into a new species for other uses.
  #ifdef VPIC_GLOBAL_PARTICLE_ID
    // Grab into the species and make it have IDs
    tracerspecies->has_ids = 1;
    MALLOC_ALIGNED( tracerspecies->p_id, max_local_np, 128 );
  #endif

  // If the parentspecies has annotations we should have them on the tracers as well
  #ifdef VPIC_PARTICLE_ANNOTATION
  if(parentspecies->has_annotation){
    tracerspecies->allocate_annotation_buffer(parentspecies->has_annotation);
  }
  #endif

  // Select the desired fraction of particles from the parent species and add to the tracer species
  int step = 0;
  for(int i = 0; i < parentspecies->np; i++) {
    if(f(parentspecies->p[i])) { // This particle was picked by the user provided predicate
      // Copy that particle over
      tracerspecies->p[step] = parentspecies->p[i];
      tracerspecies->np++;
      #ifdef VPIC_GLOBAL_PARTICLE_ID
        // Create an ID
        tracerspecies->p_id[step] = tracerspecies->generate_particle_id(step, tracerspecies->max_np);
      #endif

      #ifdef VPIC_PARTICLE_ANNOTATION
      if(parentspecies->has_annotation){
        for(int j=0; j<parentspecies->has_annotation; j++) {
          const float v = parentspecies->get_annotation(i,j);
          tracerspecies->set_annotation(step, j, v);
        }
      }
      #endif

      if(copyormove == Tracertype::Move) {
        // Remove from parent species
        parentspecies->p[i] = parentspecies->p[parentspecies->np-1];
        parentspecies->np--;
        // Reduce i by one to also check the particle that was copied in from the end of the array
        i--;
      } else if (copyormove == Tracertype::Copy) {
        // Copied tracers should have zero statistical weight
        tracerspecies->p[step].w = 0.;
      } else {
        ERROR(( "Invalid enum value for copyormove" ));
      }
      // Increment step
      step++;
    }
  }
  return tracerspecies;
}
