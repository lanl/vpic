/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Original version (data structures based on earlier
 *                    V4PIC versions)
 *
 */

#ifndef _species_h_
#define _species_h_

#include <common.h>
#include <particle.h> /* For particle_t, particle_mover_t */
#include <mtrand.h>   /* For mt_handle */ 

enum species_enums {
  invalid_species_id = -1
};

typedef int species_id;

typedef struct _species_t {
  species_id id;           /* Unique identifier for a species */
  int np, max_np;          /* Current and maximum local number of particles */
  particle_t * ALIGNED p;  /* Array of particles for the species */
  int nm, max_nm;          /* Number of movers in use */
  particle_mover_t *pm;    /* Array of particles for the mover list */
  float q_m;               /* Species charge to mass ratio */
  int sort_interval;       /* How often to sort the species */
  int * ALIGNED copy;      /* Static array of length (nx+2)*(ny+2)*(nz+2).  Each 
                              value corresponds to the associated particle array 
                              index of the first particle in the cell.  Array 
                              is allocated and values computed in sort_p.  
                              Purpose is for implementing collision models */ 
  struct _species_t *next; /* Next species in the list */
  char name[1];            /* Name is resized on allocation */
} species_t;

BEGIN_C_DECLS

/* In species.c */
extern int num_species( const species_t *sp_list );
extern species_id new_species( const char *name,
                               float q_m,
                               int max_local_np,
                               int max_local_ng,
                               int sort_interval,
                               species_t **sp_list );
extern void delete_species_list( species_t **sp_list );
extern species_t *find_species_id( species_id id, species_t *sp_list );
extern species_t *find_species_name( const char *name, species_t *sp_list );

extern species_t **new_species_lookup( species_t *sp_list );
extern void delete_species_lookup( species_t ***sp_array );

END_C_DECLS

#endif

