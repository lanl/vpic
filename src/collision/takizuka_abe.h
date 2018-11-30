#ifndef _takizuka_abe_h_
#define _takizuka_abe_h_

#include "collision_private.h"

typedef struct takizuka_abe {
  char * name;
  species_t  * spi;
  species_t  * spj;
  rng_pool_t * rp;
  int interval;
  double nu0;
} takizuka_abe_t;

void
apply_takizuka_abe_pipeline( takizuka_abe_t * l );

#endif /* _takizuka_abe_h_ */
