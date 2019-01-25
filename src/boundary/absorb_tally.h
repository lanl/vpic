#ifndef _absorb_tally_h_
#define _absorb_tally_h_

typedef struct absorb_tally {
  /**/  species_t     * sp_list;
  const field_array_t * fa;
  /**/  int           * tally;
} absorb_tally_t;


int
interact_absorb_tally( absorb_tally_t      * RESTRICT at,
                       species_t           * RESTRICT sp,
                       particle_t          * RESTRICT p,
                       particle_mover_t    * RESTRICT pm,
                       particle_injector_t * RESTRICT pi,
                       int                            max_pi,
                       int                            face );

#endif /* _absorb_tally_h_ */
