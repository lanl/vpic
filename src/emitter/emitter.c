#include "emitter.h"

emitter_t *
new_emitter( const char * name,
             species_t * sp,
             emission_model_t emission_model,
             int max_component,
             emitter_t **e_list ) {
  char * buf;
  emitter_t * e;
  int len;

  if( e_list==NULL ) ERROR(("Invalid emitter list."));
  len = (name==NULL) ? 0 : strlen(name);
  if( len<=0 ) ERROR(("Cannot create a nameless emitter."));
  if( find_emitter_name(name,*e_list)!=NULL )
    ERROR(("There is already a emitter named \"%s\".",name)); 
  if( max_component<1 ) {

    // BJA: commented ERROR() out for the following reason: If an
    // emitter is localized and does not reside on a given processor,
    // then the number of faces will be 0 on processors which don't
    // contain the emitter.  The proper behavior should be to silently
    // return NULL rather than throw an error condition.

    // FIXME: KJB-BOTH BEHAVIORS ARE UNCOOL.  I'LL FIX THIS LATER.

    // ERROR(("Invalid max_component requested for emitter \"%s\".",name));
    return NULL;
  }

  // Create the emitter

  // sizeof(e[0]) includes the termininating null of name
  MALLOC( buf, sizeof(e[0])+len ); e = (emitter_t *)buf;

  // Initialize the emitter

  e->sp             = sp;
  e->emission_model = emission_model;
  CLEAR( e->model_parameters, MAX_EMISSION_MODEL_SIZE );
  MALLOC_ALIGNED( e->component, max_component, 16 );
  e->n_component   = 0;
  e->max_component = max_component;

  e->next = *e_list;
  strcpy( e->name, name );

  *e_list = e;
  return e;
}

void
delete_emitter_list( emitter_t **e ) {
  emitter_t * next;

  if( e==NULL ) return;
  while( *e!=NULL ) {
    next = (*e)->next;
    FREE_ALIGNED( (*e)->component );
    FREE( (*e) );
    (*e) = next;
  }
}

emitter_t *
find_emitter_name( const char *name, emitter_t *e_list ) {
  emitter_t *e;
  if( name==NULL ) return NULL;
  LIST_FIND_FIRST(e,e_list,strcmp(e->name,name)==0);
  return e;
}

