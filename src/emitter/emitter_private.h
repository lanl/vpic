#ifndef _emitter_private_h_
#define _emitter_private_h_

#ifndef IN_emitter
#error "Do not include emitter_private.h; use emitter.h"
#endif

#include <emitter.h>

typedef void
(*emit_func_t)( /**/  void * RESTRICT              params,
                const int  * RESTRICT ALIGNED(128) component,
                int                                n_component );

typedef void
(*delete_emitter_func_t)( emitter_t * RESTRICT e );

struct emitter {
  void * params;
  emit_func_t emit;
  delete_emitter_func_t delete_e;
  int * ALIGNED(128) component;
  int n_component;
  emitter_t * next;
};

BEGIN_C_DECLS

// In emitter.c

void
checkpt_emitter_internal( const emitter_t * e );

emitter_t *
restore_emitter_internal( void * params );

emitter_t *
new_emitter_internal( void * params,
                      emit_func_t emit,
                      delete_emitter_func_t delete_e,
                      checkpt_func_t checkpt,
                      restore_func_t restore,
                      reanimate_func_t reanimate );

void
delete_emitter_internal( emitter_t * e );

END_C_DECLS

#endif // _emitter_private_h_
