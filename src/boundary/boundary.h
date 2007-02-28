#ifndef _boundary_handler_h_
#define _boundary_handler_h_

#include <common.h>
#include <mtrand.h>
#include <grid.h>
#include <field.h>
#include <particle.h>
#include <species.h>
#include <mtrand.h>

BEGIN_C_DECLS

/* in maxwellian_reflux.c */

typedef struct _maxwellian_reflux_t {
  int nspec;                    /* number of species interacting with boundary */ 
  species_id id[32];            /* array of species ids */ 
  double ut_perp[32];           /* array of perp thermal velocities. */ 
  double ut_para[32];           /* array of parallel thermal velocities */ 
} maxwellian_reflux_t; 

void maxwellian_reflux( void * params, particle_t *r, particle_mover_t *pm, 
                        field_t *f, accumulator_t *a, const grid_t *g,
                        species_t *s, particle_injector_t **ppi, mt_handle rng, int face ); 

/* In absorb_tally.c */ 

typedef struct _absorb_tally_t {
  int nspec;                    /* number of species interacting with boundary */ 
  species_id id[32];            /* array of species ids */ 
  int nabs[32];                 /* array of numbers of particles absorbed. */ 
} absorb_tally_t; 

void absorb_tally( void * params, particle_t *r, particle_mover_t *pm, 
                   field_t *f, accumulator_t *a, const grid_t *g,
                   species_t *s, particle_injector_t **ppi, mt_handle rng, int face ); 

/* In absorb_tally.c */ 

typedef struct _link_boundary_t {
  char fbase[256];  /* base of file name to contain link info */
  double n_out;     /* number of writes so far on this node (double to accomodate long long runs) */
} link_boundary_t;

void link_boundary( void * params, particle_t *r, particle_mover_t *pm, 
                    field_t *f, accumulator_t *a, const grid_t *g,
                    species_t *s, particle_injector_t **ppi, mt_handle rng, int face ); 

END_C_DECLS

#endif 
