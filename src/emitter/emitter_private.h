#ifndef _emitter_private_h_
#define _emitter_private_h_

#ifndef IN_emitter
#error "Do not include emitter_private.h; use emitter.h"
#endif

#include "emitter.h"

BEGIN_C_DECLS

// In structors.c

void
checkpt_emitter_internal( const emitter_t * e );

emitter_t *
restore_emitter_internal( void );

emitter_t *
new_emitter_internal( emitter_t ** e_list,
                      emit_func_t emit,
                      void * params,
                      delete_emitter_func_t delete_e,
                      checkpt_func_t checkpt,
                      restore_func_t restore,
                      reanimate_func_t reanimate );

void
delete_emitter_internal( emitter_t * e );

END_C_DECLS

#endif // _emitter_private_h_
