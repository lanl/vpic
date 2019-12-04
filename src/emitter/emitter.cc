#define IN_emitter
#include "emitter_private.h"

/* Private interface **********************************************************/

void
checkpt_emitter_internal( const emitter_t * RESTRICT e ) {
  CHECKPT( e, 1 );
  CHECKPT_SYM( e->emit );
  CHECKPT_SYM( e->delete_e );
  CHECKPT_ALIGNED( e->component, e->n_component, 128 );
  CHECKPT_PTR( e->next );
}

emitter_t *
restore_emitter_internal( void * params ) {
  emitter_t * e;
  RESTORE( e );
  e->params = params;
  RESTORE_SYM( e->emit );
  RESTORE_SYM( e->delete_e );
  RESTORE_ALIGNED( e->component );
  RESTORE_PTR( e->next );
  return e;
}

emitter_t *
new_emitter_internal( void * params,
                      emit_func_t emit,
                      delete_emitter_func_t delete_e,
                      checkpt_func_t checkpt,
                      restore_func_t restore,
                      reanimate_func_t reanimate ) {
  emitter_t * e;
  MALLOC( e, 1 );
  CLEAR( e, 1 );
  e->params   = params;
  e->emit     = emit;
  e->delete_e = delete_e;
  /* next set by append_emitter */
  REGISTER_OBJECT( e, checkpt, restore, reanimate );
  return e;
}

void
delete_emitter_internal( emitter_t * e ) {
  UNREGISTER_OBJECT( e );
  FREE_ALIGNED( e->component );
  FREE( e );
}

/* Public interface ***********************************************************/

int
num_emitter( const emitter_t * RESTRICT e_list ) {
  const emitter_t * RESTRICT e;
  int n = 0;
  LIST_FOR_EACH( e, e_list ) n++;
  return n;
}

void
apply_emitter_list( emitter_t * RESTRICT e_list ) {
  emitter_t * e;
  LIST_FOR_EACH( e, e_list )
    e->emit( e->params, e->component, e->n_component );
}

void
delete_emitter_list( emitter_t * e_list ) {
  emitter_t * e;
  while( e_list ) {
    e = e_list;
    e_list = e_list->next;
    e->delete_e( e );
  }
}

emitter_t *
append_emitter( emitter_t * e,
                emitter_t ** e_list ) {
  emitter_t * ee;
  if( !e || !e_list ) ERROR(( "Bad args" ));
  LIST_FOR_EACH( ee, *e_list ) if( ee==e ) return e;
  if( e->next ) ERROR(( "Emitter already in a list" ));
  e->next = *e_list;
  *e_list = e;
  return e;
}

int32_t * ALIGNED(128)
size_emitter( emitter_t * RESTRICT e,
              int n_component ) {
  if( !e || e->n_component || n_component<0 ) ERROR(( "Bad args" ));
  MALLOC_ALIGNED( e->component, n_component, 128 );
  e->n_component = n_component;
  return e->component;
}

