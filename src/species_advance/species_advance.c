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

/* Though these functions are not part of the public API, they must
   not be declared static. */

void
checkpt_species( const species_t * sp ) {
  CHECKPT( (const char *)sp, sizeof(*sp)+strlen(sp->name) );
  checkpt_data( sp->p,
                sp->np*sizeof(particle_t),
                sp->max_np*sizeof(particle_t), 1, 1, 128 );
  checkpt_data( sp->pm,
                sp->nm*sizeof(particle_mover_t),
                sp->max_nm*sizeof(particle_mover_t), 1, 1, 128 );
  CHECKPT_ALIGNED( sp->partition,
                   (sp->g->nx+2)*(sp->g->ny+2)*(sp->g->nz+2)+1, 128 );
  CHECKPT_PTR( sp->g );
  CHECKPT_PTR( sp->next );
}

species_t *
restore_species( void ) {
  species_t * sp;
  RESTORE( sp );
  sp->p  = (particle_t *)      restore_data();
  sp->pm = (particle_mover_t *)restore_data();
  RESTORE_ALIGNED( sp->partition );
  RESTORE_PTR( sp->g );
  RESTORE_PTR( sp->next );
  return sp;
}

int
num_species( const species_t * sp_list ) {
  return sp_list ? sp_list->id+1 : 0;
}

species_t *
new_species( const char * name,
             float q,
             float m,
             int max_local_np,
             int max_local_nm,
             int sort_interval,
             int sort_out_of_place,
             grid_t * g,
             species_t ** sp_list ) {
  species_t * sp;
  char * buf;
  int len = (!name ? 0 : strlen(name)); // len does not incl terminal '\0'

  if( len<=0 ) ERROR(( "Cannot create a nameless species" ));
  if( !sp_list ) ERROR(( "NULL species list" ));
  if( find_species_name( name, *sp_list ) ) 
    ERROR(( "There is already a species named \"%s\".", name ));

  if( max_local_np<1 ) ERROR(( "Bad max_local_np (%i) requested for \"%s\".",
                               max_local_np, name ));

  if( max_local_nm<1 ) ERROR(( "Bad max_local_nm (%i) requested for \"%s\".",
                               max_local_nm, name ));
  if( !g ) ERROR(( "NULL grid" ));
  if( (*sp_list) && (*sp_list)->g!=g )
    ERROR(( "Species lists should all use the same grid" ));


  MALLOC( buf, sizeof(*sp) + len ); // sizeof(*sp) includes terminal '\0'
  sp = (species_t *)buf;

  sp->id = num_species(*sp_list);

  sp->q = q;
  sp->m = m;

  sp->np = 0;
  sp->max_np = max_local_np;
  MALLOC_ALIGNED( sp->p, max_local_np, 128 );

  sp->nm = 0;
  sp->max_nm = max_local_nm;
  MALLOC_ALIGNED( sp->pm, max_local_nm, 128 );

  sp->sort_interval = sort_interval;
  sp->sort_out_of_place = sort_out_of_place;
  MALLOC_ALIGNED( sp->partition, (g->nx+2)*(g->ny+2)*(g->nz+2)+1, 128 );

  sp->g = g;   
  sp->next = *sp_list;
  strcpy( sp->name, name );

  REGISTER_OBJECT( sp, checkpt_species, restore_species, NULL );

  *sp_list = sp;
  return sp;
}

static void
delete_species( species_t * sp ) {
  if( sp==NULL ) return;
  UNREGISTER_OBJECT( sp );
  FREE_ALIGNED( sp->partition );
  FREE_ALIGNED( sp->pm );
  FREE_ALIGNED( sp->p );
  FREE( sp );
}

void
delete_species_list( species_t ** sp_list ) {
  species_t * sp;
  if( sp_list==NULL ) return;
  while( *sp_list!=NULL ) {
    sp = *sp_list;
    *sp_list = sp->next;
    delete_species( sp );
  }
}

species_t *
find_species_id( species_id id,
                 species_t * sp_list ) {
  species_t * sp;
  LIST_FIND_FIRST( sp, sp_list, sp->id==id );
  return sp;
}

species_t *
find_species_name( const char * name,
                   species_t * sp_list ) {
  species_t * sp;
  if( name==NULL ) return NULL;
  LIST_FIND_FIRST( sp, sp_list, strcmp( sp->name, name )==0 );
  return sp;
}
