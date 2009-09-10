#include "emitter.h"

int
num_emitter( const emitter_t * RESTRICT e_list ) {
  const emitter_t * RESTRICT e;
  int n = 0;
  LIST_FOR_EACH( e, e_list ) n++;
  return n;
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

void
size_emitter( emitter_t * RESTRICT e,
              int n_component ) {
  if( !e || e->n_component || n_component<0 ) ERROR(( "Bad args" ));
  MALLOC_ALIGNED( e->component, n_component, 128 );
  e->n_component = n_component;
}

/******************************************************************************/

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
  RESTORE_SYM( e->emit );
  RESTORE_SYM( e->delete_e );
  RESTORE_ALIGNED( e->component );
  RESTORE_PTR( e->next );
  e->params = params;
  return e;
}

emitter_t *
new_emitter_internal( emitter_t ** e_list,
                      emit_func_t emit,
                      delete_emitter_func_t delete_e,
                      void * params,
                      checkpt_func_t checkpt,
                      restore_func_t restore,
                      reanimate_func_t reanimate ) {
  emitter_t * e;

  if( !e_list || !emit || !delete_e ) ERROR(( "Bad args" ));

  MALLOC( e, 1 );
  CLEAR( e, 1 );

  e->emit     = emit;
  e->delete_e = delete_e;
  e->params   = params;
  e->next     = NULL;     /* Set by define_emitter */

  REGISTER_OBJECT( e, checkpt, restore, reanimate );
  return e;
}

void
delete_emitter_internal( emitter_t * e ) {
  if( !e ) return;
  UNREGISTER_OBJECT( e );
  FREE_ALIGNED( e->component );
  FREE( e );
}

