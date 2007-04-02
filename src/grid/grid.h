/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Original version
 *
 */

#ifndef _grid_h_
#define _grid_h_

#include <common.h>
#include <mp.h> /* For mp_handle */
#include <mtrand.h> /* For mt_handle...SEE NOTE BELOW ABOUT BABIES */


/* Define a "pointer to boundary handler function" type. */

/* THESE PREDECLARATIONS ARE UGLY AND MAKE BABY JESUS CRY */

struct particle;
struct particle_mover;
struct field;
struct accumulator;
struct grid;
struct species;
struct particle_injector;

typedef void (*boundary_handler_t)( void                  * params,
                                    struct particle       * r,
                                    struct particle_mover * pm,       
                                    struct field          * f,
                                    struct accumulator    * a,
                                    const struct grid     * g,
                                    struct species        * s, 
                                    struct particle_injector ** ppi, 
                                    mt_handle     rng, 
                                    int           face );

enum boundary_handler_enums {
  INVALID_BOUNDARY       = 0xBADF00D,
  MAX_BOUNDARY_DATA_SIZE = 1024 /* Sized to hold data characterizing a boundary model. */
};

/* FIXME: MODEL_PARAMETER DOES NOT HAVE GUARANTEED NICE ALIGNMENT */

typedef struct boundary {
 boundary_handler_t handler;
 char params[MAX_BOUNDARY_DATA_SIZE]; 
} boundary_t;

#define BOUNDARY(i,j,k) INDEX_FORTRAN_3(i,j,k,-1,1,-1,1,-1,1)

enum grid_enums {

  /* Phase 2 boundary conditions */
  anti_symmetric_fields = -1, /* E_tang = 0 */
  pec_fields            = -1,
  metal_fields          = -1,
  symmetric_fields      = -2, /* B_tang = 0, B_norm = 0 */
  pmc_fields            = -3, /* B_tang = 0, B_norm floats */
  absorb_fields         = -4, /* Gamma = 0 */

  /* Phase 3 boundary conditions */
  reflect_particles = -1, /* Cell boundary should reflect particles */
  absorb_particles  = -2  /* Cell boundary should absorb particles */

  /* Symmetry in the field boundary conditions refers to image charge sign
     Anti-symmetric -> Image charges are opposite signed (ideal metal)
                       Boundary rho/j are accumulated over partial cell+image
     Symmetric      -> Image charges are same signed (symmetry plane or pmc)
                       Boundary rho/j are accumulated over partial cell+image
     Absorbing      -> No image charges
                       Boundary rho/j are accumulated over partial cell only

     rho     -> Anti-symmetric      | rho     -> Symmetric
     jf_tang -> Anti-symmetric      | jf_tang -> Symmetric
     E_tang  -> Anti-symmetric      | E_tang  -> Symmetric
     B_norm  -> Anti-symmetric + DC | B_norm  -> Symmetric      (see note)
     B_tang  -> Symmetric           | B_tang  -> Anti-symmetric
     E_norm  -> Symmetric           | E_norm  -> Anti-symmetric (see note)
     div B   -> Symmetric           | div B   -> Anti-symmetric
     
     Note: B_norm is tricky. For a symmetry plane, B_norm on the boundary must
     be zero as there are no magnetic charges (a non-zero B_norm would imply
     an infinitesimal layer of magnetic charge). However, if a symmetric
     boundary is interpreted as a perfect magnetic conductor, B_norm could be
     present due to magnetic conduction surface charges. Even though there are
     no bulk volumetric magnetic charges to induce a surface magnetic charge,
     I think that radiation/waveguide modes/etc could (the total surface
     magnetic charge in the simulation would be zero though). As a result,
     symmetric and pmc boundary conditions are treated separately. Symmetric
     and pmc boundaries are identical except the symmetric boundaries
     explicitly zero boundary B_norm. Note: anti-symmetric and pec boundary
     conditions would have the same issue if norm E was located directly on
     the boundary. However, it is not so this problem does not arise.

     Note: Absorbing boundary conditions make no effort to clean divergence
     errors on them. They assume that the ghost div b is zero and force the
     surface div e on them to be zero. This means ghost norm e can be set to
     any value on absorbing boundaries. */ 
};

typedef struct grid {
  mp_handle mp;           /* Communications handle */
  float dt, cvac, eps0;   /* System of units */
  float damp;             /* Radiation damping parameter */

  /* Phase 2 grid data structures */
  float x0, y0, z0;       /* Corner of cell 1,1,1 */
  float dx, dy, dz;       /* Cell dimensions */
  int nx, ny, nz;         /* Number of cells in domain */
  int bc[27];             /* (-1:1,-1:1,-1:1) FORTRAN indexed array of
			     boundary conditions to apply at domain edge
                             0 ... nproc-1 ... comm boundary condition
                             <0 ... locally applied boundary condition */

  /* Phase 3 grid data structures */
  /* NOTE: LOCAL_CELL_ID LIMITS NUMBER OF CELLS TO 2^31 (INCLUDING GHOSTS) PER
    NODE. CELL ADJACENCY INDEXING FURTHER LIMITS TO 2^31 / 6. THE EMITTER
    COMPONENT ID STRATEGY FURTHER LIMITS TO 2^27 PER NODE. THE LIMIT
    IS 2^64 OVER ALL NODES THOUGH. */

  INT64_TYPE * ALIGNED range;    /* (0:nproc) indexed array giving range of global
                                    indexes of cells owned by each processor.
                                    Replicated on each processor.
                                    (range[rank]:range[rank+1]-1) are global cells
                                    owned by processor "rank". Note:
                                    range[rank+1]-range[rank] <~ 2^31 / 6 */
  INT64_TYPE * ALIGNED neighbor; /* (0:5,0:local_num_cells-1) FORTRAN indexed array
                                    neighbor(0:5,lidx) are the global indexes of
                                    neighboring cells of the cell with local index
                                    "lidx". Negative if neighbor is a boundary
                                    condition. */
  INT64_TYPE rangel, rangeh;   /* Redundant for move_p performance reasons:
                                  rangel = range[rank]
                                  rangeh = range[rank+1]-1.
                                  Note: rangeh-rangel <~ 2^31 / 6 */

  /* Enable user-defined boundary handlers */
  int nb;                 /* Number of custom boundary conditions */
  boundary_t * boundary;  /* Head of array of boundary_t. */
} grid_t;

BEGIN_C_DECLS

/* In boundary_handler.c */

/* Note that boundarys _copy_ the state given by ip into their own storage.
   Thus, to change the state of a given boundary, one has to look up 
   the boundary and change it there.
   Example usage:

   refux_boundary_params_t bc;

   bc.species[0] = ion;
   bc.ux[0]      = sqrt(k*T/(m*c^2))
   ...
   reflux_boundary = add_boundary( g, reflux_handler, &bc ); */

#define add_boundary(g,bh,ip) \
  IUO_add_boundary((g),(bh),(ip),sizeof(*(ip)))
int
IUO_add_boundary( grid_t *g,
		  boundary_handler_t bh,
		  const void * initial_params,
		  int sz );

/* In structors.c */

grid_t *
new_grid(void);

void
delete_grid( grid_t **g );

/* In ops.c */

void
size_grid( grid_t * g, int lnx, int lny, int lnz );

void
join_grid( grid_t * g, int bound, int rank );

void
set_fbc( grid_t *g, int bound, int fbc );

void
set_pbc( grid_t *g, int bound, int pbc );

/* In partition.c */

void
partition_periodic_box( grid_t *g,
			double glx, double gly, double glz,
			int gnx, int gny, int gnz,
			int gpx, int gpy, int gpz );

void
partition_absorbing_box( grid_t *g,
			 double glx, double gly, double glz,
			 int gnx, int gny, int gnz,
			 int gpx, int gpy, int gpz,
			 int pbc );

void
partition_metal_box( grid_t *g,
		     double glx, double gly, double glz,
		     int gnx, int gny, int gnz,
		     int gpx, int gpy, int gpz );

END_C_DECLS

#endif
