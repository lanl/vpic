#ifndef _energy_f_pipeline_h_
#define _energy_f_pipeline_h_

#ifndef IN_energy_f_pipeline
#error "Only include energy_f_pipeline.h in energy_f_pipeline source files."
#endif

#include "../sfa_private.h"

typedef struct pipeline_args
{
    const field_t *ALIGNED( 128 ) f;
    const sfa_params_t *p;
    const grid_t *g;
    double en[MAX_PIPELINE + 1][6];
} pipeline_args_t;

#define DECLARE_STENCIL()                                                      \
    const field_t *ALIGNED( 128 ) f = args->f;                                 \
    const material_coefficient_t *ALIGNED( 128 ) m = args->p->mc;              \
    const grid_t *g = args->g;                                                 \
    const int nx = g->nx, ny = g->ny, nz = g->nz;                              \
                                                                               \
    const field_t *ALIGNED( 16 ) f0;                                           \
    const field_t *ALIGNED( 16 ) fx, *ALIGNED( 16 ) fy, *ALIGNED( 16 ) fz;     \
    const field_t *ALIGNED( 16 ) fyz, *ALIGNED( 16 ) fzx, *ALIGNED( 16 ) fxy;  \
    double en_ex = 0, en_ey = 0, en_ez = 0, en_bx = 0, en_by = 0, en_bz = 0;   \
    int x, y, z

#define f( x, y, z ) f[VOXEL( x, y, z, nx, ny, nz )]

#define INIT_STENCIL()                                                         \
    f0 = &f( x, y, z );                                                        \
    fx = &f( x + 1, y, z );                                                    \
    fy = &f( x, y + 1, z );                                                    \
    fz = &f( x, y, z + 1 );                                                    \
    fyz = &f( x, y + 1, z + 1 );                                               \
    fzx = &f( x + 1, y, z + 1 );                                               \
    fxy = &f( x + 1, y + 1, z )

#define NEXT_STENCIL()                                                         \
    f0++;                                                                      \
    fx++;                                                                      \
    fy++;                                                                      \
    fz++;                                                                      \
    fyz++;                                                                     \
    fzx++;                                                                     \
    fxy++;                                                                     \
    x++;                                                                       \
    if ( x > nx )                                                              \
    {                                                                          \
        /**/ y++;                                                                \
        x = 1;                                                                 \
        if ( y > ny )                                                          \
            z++;                                                               \
        if ( y > ny )                                                          \
            y = 1;                                                             \
        INIT_STENCIL();                                                        \
    }

#define REDUCE_EN()                                                            \
    en_ex += 0.25 * ( m[f0->ematx].epsx * f0->ex * f0->ex +                    \
                      m[fy->ematx].epsx * fy->ex * fy->ex +                    \
                      m[fz->ematx].epsx * fz->ex * fz->ex +                    \
                      m[fyz->ematx].epsx * fyz->ex * fyz->ex );                \
    en_ey += 0.25 * ( m[f0->ematy].epsy * f0->ey * f0->ey +                    \
                      m[fz->ematy].epsy * fz->ey * fz->ey +                    \
                      m[fx->ematy].epsy * fx->ey * fx->ey +                    \
                      m[fzx->ematy].epsy * fzx->ey * fzx->ey );                \
    en_ez += 0.25 * ( m[f0->ematz].epsz * f0->ez * f0->ez +                    \
                      m[fx->ematz].epsz * fx->ez * fx->ez +                    \
                      m[fy->ematz].epsz * fy->ez * fy->ez +                    \
                      m[fxy->ematz].epsz * fxy->ez * fxy->ez );                \
    en_bx += 0.5 * ( m[f0->fmatx].rmux * f0->cbx * f0->cbx +                   \
                     m[fx->fmatx].rmux * fx->cbx * fx->cbx );                  \
    en_by += 0.5 * ( m[f0->fmaty].rmuy * f0->cby * f0->cby +                   \
                     m[fy->fmaty].rmuy * fy->cby * fy->cby );                  \
    en_bz += 0.5 * ( m[f0->fmatz].rmuz * f0->cbz * f0->cbz +                   \
                     m[fz->fmatz].rmuz * fz->cbz * fz->cbz )

void energy_f_pipeline_scalar( pipeline_args_t *args, int pipeline_rank,
                               int n_pipeline );

#endif // _energy_f_pipeline_h_
