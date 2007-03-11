/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Original version (data structures were revised from
 *                    earlier V4PIC versions)
 *
 */

#ifndef _particle_h_
#define _particle_h_

#include <common.h>
#include <grid.h> 
#include <field.h> 
#include <mtrand.h> /* UGLY: For mt_handle in boundary_p */

typedef struct _particle_t {
  float dx, dy, dz; /* Particle position in cell coordinates (on [-1,1]) */
  INT32_TYPE i;     /* Index of cell containing the particle */
  float ux, uy, uz; /* Particle normalized momentum */
  float q;          /* Particle charge */
} particle_t;

/* WARNING: FUNCTIONS THAT USE A PARTICLE_MOVER ASSUME THAT
   EVERYBODY WHO USES THAT PARTICLE MOVER WILL BE GIVEN A
   POINTER TO PARTICLE 0 */

typedef struct _particle_mover_t {
  float dispx, dispy, dispz; /* Displacement of particle */
  INT32_TYPE i;           /* Index of the particle to move */
} particle_mover_t;

typedef struct _particle_injector_t {
  float dx, dy, dz; INT32_TYPE i;
  float ux, uy, uz, q;
  float dispx, dispy, dispz;
} particle_injector_t;

BEGIN_C_DECLS
struct _species_t; 

/* In boundary_p.c */
extern int boundary_p( particle_mover_t * RESTRICT ALIGNED pm,
                       int nm,
                       int max_nm,
                       particle_t * RESTRICT ALIGNED p,
                       int np,
                       int max_np,
                       field_t * RESTRICT ALIGNED f,
                       accumulator_t * RESTRICT ALIGNED a,
                       const grid_t * RESTRICT g,
		       struct _species_t * RESTRICT sp,
                       mt_handle rng );

/* In hydro_p.c */
extern void accumulate_hydro_p( hydro_t * RESTRICT ALIGNED h0,
                                const particle_t * RESTRICT ALIGNED p,
                                int np,
                                float q_m,
                                const interpolator_t * RESTRICT ALIGNED f0,
                                const grid_t * RESTRICT g );

/* In move_p.c */
extern int inject_p( particle_t * RESTRICT ALIGNED p,        /* Array to inject into */
                     int i,                                  /* Where to inject */
                     particle_mover_t * RESTRICT ALIGNED pm, /* Free mover */
                     field_t * RESTRICT ALIGNED f,
                     accumulator_t * RESTRICT ALIGNED a,
                     const particle_injector_t * RESTRICT pi,
                     const grid_t * RESTRICT g );

/* In sort_p.c */
extern void sort_p( struct _species_t * RESTRICT sp, 
                    const grid_t * RESTRICT g );

/* In rho_p.c */
extern void accumulate_rho_p( field_t * RESTRICT ALIGNED f,
                              const particle_t * RESTRICT ALIGNED p,
                              int np,
                              const grid_t * RESTRICT g );

/* In pstructors.c */
extern particle_t * ALIGNED new_particle_array( int np );
extern particle_mover_t * ALIGNED new_particle_mover( int nm );
extern void delete_particle_array( particle_t ** ALIGNED p );
extern void delete_particle_mover( particle_mover_t ** ALIGNED pm );

/* In advance_p.cpp */
extern int advance_p( particle_t * RESTRICT ALIGNED p,
                      int n,
                      const float q_m,
                      particle_mover_t * RESTRICT ALIGNED pm,
                      int nm,       
                      accumulator_t * RESTRICT ALIGNED a,
                      const interpolator_t * RESTRICT ALIGNED f,
                      const grid_t * RESTRICT g );

/* In center_p.cpp */
extern void center_p( particle_t * RESTRICT ALIGNED p,
                      int np,
                      const float q_m,
                      const interpolator_t * RESTRICT ALIGNED f,
                      const grid_t * RESTRICT g );

/* In energy.cpp */
extern double energy_p( const particle_t * RESTRICT ALIGNED p,
                        int np,
                        float q_m,
                        const interpolator_t * RESTRICT ALIGNED f,
                        const grid_t * RESTRICT g );

/* In uncenter_p.cpp */
extern void uncenter_p( particle_t * RESTRICT ALIGNED p,
                        int np,
                        const float q_m,
                        const interpolator_t * RESTRICT ALIGNED f,
                        const grid_t * RESTRICT g );

/* INTERNAL USE ONLY FUNCTIONS */

/* In move_p.c */
extern int move_p( particle_t * RESTRICT ALIGNED p, 
                   particle_mover_t * RESTRICT ALIGNED m,
		   accumulator_t * RESTRICT ALIGNED a,
		   const grid_t * RESTRICT g );
extern int remove_p( particle_t * RESTRICT ALIGNED r,
                     particle_t * RESTRICT ALIGNED p,
                     int np,
                     field_t * RESTRICT ALIGNED f,
                     const grid_t * RESTRICT g );

END_C_DECLS

#endif
