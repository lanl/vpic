#ifndef _langevin_h_
#define _langevin_h_

typedef struct langevin
{
  species_t  * sp;
  rng_pool_t * rp;
  float kT;
  float nu;
  int interval;
} langevin_t;

#endif /* _langevin_h_ */
