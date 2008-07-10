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

int
num_species( const species_t *sp_list ) {
  const species_t * sp;
  int n = 0;
  LIST_FOR_EACH( sp, sp_list ) n++;
  return n;
}

species_t *
new_species( const char *name,
             float q_m,
             int max_local_np,
             int max_local_nm,
             int sort_interval,
             int sort_out_of_place,
             species_t **sp_list ) {
  char * buf;
  species_t * sp;
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
  MALLOC( buf, sizeof(sp[0])+len ); sp = (species_t *)buf;
  sp->id = num_species(*sp_list);
  sp->np = 0;
  sp->max_np = max_local_np;
  MALLOC_ALIGNED( sp->p, max_local_np, 128 );
  sp->nm = 0;
  sp->max_nm = max_local_nm;
  MALLOC_ALIGNED( sp->pm, max_local_nm, 128 );
  sp->q_m = q_m;
  sp->sort_interval = sort_interval;
  sp->sort_out_of_place = sort_out_of_place;
  sp->partition = NULL; // FIXME: ALLOCATE THIS HERE! SP SHOULD KNOW ABOUT GRID
  strcpy( sp->name, name );
  
  sp->next = *sp_list;
  *sp_list = sp;
  return sp;
}

void
delete_species_list( species_t **sp_list ) {
  species_t *sp;

  if( sp_list==NULL ) return;
  while( *sp_list!=NULL ) {
    sp = *sp_list;
    FREE_ALIGNED( sp->partition );
    FREE_ALIGNED( sp->pm );
    FREE_ALIGNED( sp->p );
    *sp_list = sp->next;
    FREE( sp );
  }
}

species_t *
find_species_id( species_id id,
                 species_t *sp_list ) {
  species_t *sp;
  LIST_FIND_FIRST(sp,sp_list,sp->id==id);
  return sp;
}

species_t *
find_species_name( const char *name,
                   species_t *sp_list ) {
  species_t *sp;
  if( name==NULL ) return NULL;
  LIST_FIND_FIRST(sp,sp_list,strcmp(sp->name,name)==0);
  return sp;
}
