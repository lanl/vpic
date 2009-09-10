#ifndef _field_advance_private_h_
#define _field_advance_private_h_

#ifndef IN_field_advance
#error "Do not include field_advance_private.h; use field_advance.h"
#endif

#include "../field_advance.h"

BEGIN_C_DECLS

/* Checkpoint the symbols for the given field advance kernel */

void
checkpt_field_advance_kernels( const field_advance_kernels_t * kernel );

/* Checkpoint the symbols for the given field advance kernel */

void
restore_field_advance_kernels( field_advance_kernels_t * kernel );

END_C_DECLS

#endif // _field_advance_private_h_
