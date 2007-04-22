/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Original version
 *
 */

#include <species.h>

// Note: new_species added the created species to the head of the
// species list. Further the species ids are simply incremented from
// the previous head of the list. The first species is numbered 0. As
// a result, the total number of species in the species list is the id
// of the species at the head of the list plus one.

int num_species( const species_t *sp_list ) {
  if( sp_list==NULL ) return 0;
  return sp_list->id+1;
}

species_id
new_species( const char *name,
             float q_m,
             int max_local_np,
             int max_local_nm,
             int sort_interval,
             species_t **sp_list ) {
  species_t *sp;
  int len;

  if( sp_list==NULL ) ERROR(("Invalid species list."));
  // Note: strlen does not include terminating NULL
  len = (name==NULL) ? 0 : strlen(name);
  if( len<=0 ) ERROR(("Cannot create a nameless species."));
  if( find_species_name(name,*sp_list)!=NULL ) 
    ERROR(("There is already a species named \"%s\".",name));
  if( max_local_np<1 )
    ERROR(("Invalid max_local_np requested for species \"%s\".",name));
  if( max_local_nm<1 )
    ERROR(("Invalid max_local_nm requested for species \"%s\".",name));
  
  // Note: Since a sp->name is declared as a 1-element char array, the
  // terminating NULL is included in sizeof(species_t)
  sp = (species_t *)malloc(sizeof(species_t)+len);
  if( sp==NULL ) ERROR(("Unable to allocate species \"%s\".", name));

  sp->id = num_species(*sp_list);
  sp->np = 0;
  sp->max_np = max_local_np;
  sp->p = new_particle_array( max_local_np );
  if( sp->p==NULL )
    ERROR(("Could not allocate particle array for species \"%s\".",name));
  sp->nm = 0;
  sp->max_nm = max_local_nm;
  sp->pm = new_particle_mover( max_local_nm );
  if( sp->pm==NULL )
    ERROR(("Could not allocate mover array for species \"%s\".",name));
  sp->q_m = q_m;
  sp->sort_interval = sort_interval;
  sp->copy = NULL;
  sp->next = *sp_list;
  strcpy( sp->name, name );
  
  *sp_list = sp;
  return sp->id;
}

void delete_species_list( species_t **sp_list ) {
  species_t *sp;

  if( sp_list==NULL ) return;
  while( *sp_list!=NULL ) {
    sp = *sp_list;
    delete_particle_array(&sp->p);
    delete_particle_mover(&sp->pm);
    free_aligned(sp->copy); 
    *sp_list = sp->next;
    free(sp);
  }
}

species_t *find_species_id( species_id id, species_t *sp_list ) {
  species_t *sp;
  LIST_FIND_FIRST(sp,sp_list,sp->id==id);
  return sp;
}

species_t *find_species_name( const char *name, species_t *sp_list ) {
  species_t *sp;
  if( name==NULL ) return NULL;
  LIST_FIND_FIRST(sp,sp_list,strcmp(sp->name,name)==0);
  return sp;
}

/*****************************************************************************/

species_t **new_species_lookup( species_t *sp_list ) {
  species_t *sp;
  species_t **sp_lookup;

  if( sp_list==NULL ) ERROR(("Empty species list"));
  sp_lookup = (species_t **)malloc(num_species(sp_list)*sizeof(species_t **));
  if( sp_lookup==NULL ) ERROR(("Could not allocate species lookup"));
  LIST_FOR_EACH(sp,sp_list) sp_lookup[sp->id] = sp;

  return sp_lookup;
}

void delete_species_lookup( species_t ***sp_lookup ) {
  if( sp_lookup==NULL ) return;
  if( *sp_lookup==NULL ) return;
  free(*sp_lookup);
  *sp_lookup = NULL;
}

