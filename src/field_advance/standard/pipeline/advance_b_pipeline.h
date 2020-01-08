#ifndef _advance_b_pipeline_h_
#define _advance_b_pipeline_h_

#ifndef IN_advance_b_pipeline
#error "Only include advance_b_pipeline.h in advance_b_pipeline source files."
#endif

#include "../../field_advance.h"

typedef struct pipeline_args
{
    field_t *ALIGNED( 128 ) f;
    const grid_t *g;
    float frac;
} pipeline_args_t;

#define DECLARE_STENCIL()                                                      \
    field_t *ALIGNED( 128 ) f = args->f;                                       \
    const grid_t *g = args->g;                                                 \
                                                                               \
    const int nx = g->nx;                                                      \
    const int ny = g->ny;                                                      \
    const int nz = g->nz;                                                      \
                                                                               \
    const float frac = args->frac;                                             \
    const float px = ( nx > 1 ) ? frac * g->cvac * g->dt * g->rdx : 0;         \
    const float py = ( ny > 1 ) ? frac * g->cvac * g->dt * g->rdy : 0;         \
    const float pz = ( nz > 1 ) ? frac * g->cvac * g->dt * g->rdz : 0;         \
                                                                               \
    field_t *ALIGNED( 16 ) f0;                                                 \
    field_t *ALIGNED( 16 ) fx, *ALIGNED( 16 ) fy, *ALIGNED( 16 ) fz;           \
    int x, y, z

#define f( x, y, z ) f[VOXEL( x, y, z, nx, ny, nz )]

#define INIT_STENCIL()                                                         \
    f0 = &f( x, y, z );                                                        \
    fx = &f( x + 1, y, z );                                                    \
    fy = &f( x, y + 1, z );                                                    \
    fz = &f( x, y, z + 1 )

#define NEXT_STENCIL()                                                         \
    f0++;                                                                      \
    fx++;                                                                      \
    fy++;                                                                      \
    fz++;                                                                      \
    x++;                                                                       \
    if ( x > nx )                                                              \
    {                                                                          \
        y++;                                                                   \
        x = 1;                                                                 \
        if ( y > ny )                                                          \
            z++;                                                               \
        if ( y > ny )                                                          \
            y = 1;                                                             \
        INIT_STENCIL();                                                        \
    }

// WTF!  Under -ffast-math, gcc-4.1.1 thinks it is okay to treat the
// below as
//   f0->cbx = ( f0->cbx + py*( blah ) ) - pz*( blah )
// even with explicit parenthesis are in there!  Oh my ...
// -fno-unsafe-math-optimizations must be used

#define UPDATE_CBX()                                                           \
    f0->cbx -= ( py * ( fy->ez - f0->ez ) - pz * ( fz->ey - f0->ey ) )
#define UPDATE_CBY()                                                           \
    f0->cby -= ( pz * ( fz->ex - f0->ex ) - px * ( fx->ez - f0->ez ) )
#define UPDATE_CBZ()                                                           \
    f0->cbz -= ( px * ( fx->ey - f0->ey ) - py * ( fy->ex - f0->ex ) )

void advance_b_pipeline_scalar( pipeline_args_t *args, int pipeline_rank,
                                int n_pipeline );

void advance_b_pipeline_v4( pipeline_args_t *args, int pipeline_rank,
                            int n_pipeline );

void advance_b_pipeline_v8( pipeline_args_t *args, int pipeline_rank,
                            int n_pipeline );

void advance_b_pipeline_v16( pipeline_args_t *args, int pipeline_rank,
                             int n_pipeline );

#endif // _advance_b_pipeline_h_
