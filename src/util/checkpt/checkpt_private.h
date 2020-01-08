#ifndef _checkpt_private_h_
#define _checkpt_private_h_

#ifndef IN_checkpt
#error "Do not include checkpt_private.h; use checkpt.h"
#endif

#include "checkpt.h"

struct checkpt;
typedef struct checkpt checkpt_t;

BEGIN_C_DECLS

checkpt_t* checkpt_open_rdonly( const char* name );

checkpt_t* checkpt_open_wronly( const char* name );

void checkpt_close( checkpt_t* checkpt );

void checkpt_read( checkpt_t* checkpt, void* data, size_t sz );

void checkpt_write( checkpt_t* checkpt, const void* data, size_t sz );

END_C_DECLS

#endif /* _checkpt_private_h_ */
