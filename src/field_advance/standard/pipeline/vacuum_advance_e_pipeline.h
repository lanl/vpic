#ifndef _vacuum_advance_e_pipeline_h_
#define _vacuum_advance_e_pipeline_h_

#ifndef IN_vacuum_advance_e_pipeline
#error                                                                         \
    "Only include vacuum_advance_e_pipeline.h in vacuum_advance_e_pipeline source files."
#endif

#include "../sfa_private.h"

typedef struct pipeline_args
{
    field_t* ALIGNED( 128 ) f;
    const sfa_params_t* p;
    const grid_t* g;
} pipeline_args_t;

#define DECLARE_STENCIL()                                                      \
    field_t* ALIGNED( 128 ) f = args->f;                                       \
    const material_coefficient_t* ALIGNED( 128 ) m = args->p->mc;              \
    const grid_t* g = args->g;                                                 \
    const int nx = g->nx, ny = g->ny, nz = g->nz;                              \
                                                                               \
    const float decayx = m->decayx, drivex = m->drivex;                        \
    const float decayy = m->decayy, drivey = m->drivey;                        \
    const float decayz = m->decayz, drivez = m->drivez;                        \
    const float damp = args->p->damp;                                          \
    const float px_muz =                                                       \
        ( ( nx > 1 ) ? ( 1 + damp ) * g->cvac * g->dt * g->rdx : 0 ) *         \
        m->rmuz;                                                               \
    const float px_muy =                                                       \
        ( ( nx > 1 ) ? ( 1 + damp ) * g->cvac * g->dt * g->rdx : 0 ) *         \
        m->rmuy;                                                               \
    const float py_mux =                                                       \
        ( ( ny > 1 ) ? ( 1 + damp ) * g->cvac * g->dt * g->rdy : 0 ) *         \
        m->rmux;                                                               \
    const float py_muz =                                                       \
        ( ( ny > 1 ) ? ( 1 + damp ) * g->cvac * g->dt * g->rdy : 0 ) *         \
        m->rmuz;                                                               \
    const float pz_muy =                                                       \
        ( ( nz > 1 ) ? ( 1 + damp ) * g->cvac * g->dt * g->rdz : 0 ) *         \
        m->rmuy;                                                               \
    const float pz_mux =                                                       \
        ( ( nz > 1 ) ? ( 1 + damp ) * g->cvac * g->dt * g->rdz : 0 ) *         \
        m->rmux;                                                               \
    const float cj = g->dt / g->eps0;                                          \
                                                                               \
    field_t* ALIGNED( 16 ) f0;                                                 \
    field_t *ALIGNED( 16 ) fx, *ALIGNED( 16 ) fy, *ALIGNED( 16 ) fz;           \
    int x, y, z

#define f( x, y, z ) f[VOXEL( x, y, z, nx, ny, nz )]

#define INIT_STENCIL()                                                         \
    f0 = &f( x, y, z );                                                        \
    fx = &f( x - 1, y, z );                                                    \
    fy = &f( x, y - 1, z );                                                    \
    fz = &f( x, y, z - 1 )

#define NEXT_STENCIL()                                                         \
    f0++;                                                                      \
    fx++;                                                                      \
    fy++;                                                                      \
    fz++;                                                                      \
    x++;                                                                       \
    if ( x > nx )                                                              \
    {                                                                          \
        y++;                                                                   \
        x = 2;                                                                 \
        if ( y > ny )                                                          \
            z++;                                                               \
        if ( y > ny )                                                          \
            y = 2;                                                             \
        INIT_STENCIL();                                                        \
    }

#define UPDATE_EX()                                                            \
    f0->tcax =                                                                 \
        ( py_muz * ( f0->cbz - fy->cbz ) - pz_muy * ( f0->cby - fz->cby ) ) -  \
        damp * f0->tcax;                                                       \
    f0->ex = decayx * f0->ex + drivex * ( f0->tcax - cj * f0->jfx )

#define UPDATE_EY()                                                            \
    f0->tcay =                                                                 \
        ( pz_mux * ( f0->cbx - fz->cbx ) - px_muz * ( f0->cbz - fx->cbz ) ) -  \
        damp * f0->tcay;                                                       \
    f0->ey = decayy * f0->ey + drivey * ( f0->tcay - cj * f0->jfy )

#define UPDATE_EZ()                                                            \
    f0->tcaz =                                                                 \
        ( px_muy * ( f0->cby - fx->cby ) - py_mux * ( f0->cbx - fy->cbx ) ) -  \
        damp * f0->tcaz;                                                       \
    f0->ez = decayz * f0->ez + drivez * ( f0->tcaz - cj * f0->jfz )

void vacuum_advance_e_pipeline_scalar( pipeline_args_t* args, int pipeline_rank,
                                       int n_pipeline );

void vacuum_advance_e_pipeline_v4( pipeline_args_t* args, int pipeline_rank,
                                   int n_pipeline );

void vacuum_advance_e_pipeline_v8( pipeline_args_t* args, int pipeline_rank,
                                   int n_pipeline );

void vacuum_advance_e_pipeline_v16( pipeline_args_t* args, int pipeline_rank,
                                    int n_pipeline );

#endif // _vacuum_advance_e_pipeline_h_
