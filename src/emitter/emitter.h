// FIXME: WRITE SIMPLE EMITTERS FOR THINGS LIKE CONSTANT CURRENT PARTICLE
// BEAMS

// FIXME: COULD ADJUST API OF THESE TO DO MULTIPLE SPECIES OF EMISSION
// FROM THE SAME EMITTER GEOMETRY (E.G. EMITTER_GEOMETRY OBJECT).

#ifndef _emitter_h_
#define _emitter_h_

#include "../species_advance/standard/spa.h"

struct emitter;
typedef struct emitter emitter_t;

// Every local cell has 27 components associated with it (6 faces, 12
// edges, 8 corners and cell body). All components in a local grid
// simulation can be uniquely enumerated by 32*local_cell_id +
// component_type. 32 is used instead of 27 for the following reasons:
// - Trivial to compute the component id by bit ops
// - Trivial to compute cell and component by bit ops
// - Regions 27-31 can be used for future expansion.
// The face, edge and corner, body components types are enumerated by
// a (-1:1,-1:1,-1:1) FORTRAN style indexing calculation. Note that
// this allows distinctions like which side of a cell a face is on.

#define COMPONENT_ID( local_cell, component_type ) \
  (((local_cell)<<5) | (component_type))
#define EXTRACT_LOCAL_CELL( component_id )     ((component_id)>>5)
#define EXTRACT_COMPONENT_TYPE( component_id ) ((component_id)&31)

BEGIN_C_DECLS

// In emitter.c

int
num_emitter( const emitter_t * e_list );

void
apply_emitter_list( emitter_t * e_list );

void
delete_emitter_list( emitter_t * e_list );

// Note that this append is hacked to silently return if the given
// emitter is already part of the list.  This allows the emitter
// initialization in vpic.hxx / deck_wrappers.cxx to get around
// some limitations of strict C++. 

emitter_t *
append_emitter( emitter_t * e,
                emitter_t ** e_list );

// Each emitter must be sized once and only once.  Returns
// the buffer were the emitter components should be stored.

int32_t * ALIGNED(128)
size_emitter( emitter_t * e,
              int n_component );

// In child-langmuir.c

#define CHILD_LANGMUIR sqrt(32./81.)
#define CCUBE          sqrt(1./6.)
#define IVORY          sqrt(1./6.)

emitter_t *
child_langmuir( /**/  species_t            * RESTRICT sp,  // Species to emit
                const interpolator_array_t * RESTRICT ia,  // For field interpolation
                /**/  field_array_t        * RESTRICT fa,  // For rhob accum (inject)
                /**/  accumulator_array_t  * RESTRICT aa,  // For Jf accum (aging)
                /**/  mt_rng_t             **         rng, // Random number source
                int   n_emit_per_face, // Particles to emit per face per step
                float ut_perp,         // Perpendicular normalized thermal momentum
                float ut_para,         // Parallel normalized thermal momentum
                float thresh_e_norm,   // Only emit if E_norm>thresh_e_norm
                float norm );          // Child-langmuir normalization

END_C_DECLS

#endif // _emitter_h_
