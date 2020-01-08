#ifndef _compute_curl_b_pipeline_h_
#define _compute_curl_b_pipeline_h_

#ifndef IN_compute_curl_b_pipeline
#error "Only include compute_curl_b_pipeline.h in compute_curl_b_pipeline source files."
#endif

#include "../sfa_private.h"

typedef struct pipeline_args
{
  field_t            * ALIGNED(128) f;
  const sfa_params_t *              p;
  const grid_t       *              g;
} pipeline_args_t;

#define DECLARE_STENCIL()                                        \
        field_t                * ALIGNED(128) f = args->f;       \
  const material_coefficient_t * ALIGNED(128) m = args->p->mc;   \
  const grid_t                 *              g = args->g;       \
  const int nx = g->nx, ny = g->ny, nz = g->nz;                  \
                                                                 \
  const float px = (nx>1) ? g->cvac*g->dt*g->rdx : 0;            \
  const float py = (ny>1) ? g->cvac*g->dt*g->rdy : 0;            \
  const float pz = (nz>1) ? g->cvac*g->dt*g->rdz : 0;            \
                                                                 \
  field_t * ALIGNED(16) f0;                                      \
  field_t * ALIGNED(16) fx, * ALIGNED(16) fy, * ALIGNED(16) fz;  \
  int x, y, z

#define f(x,y,z) f[ VOXEL( x, y, z, nx, ny, nz ) ]

#define INIT_STENCIL()      \
  f0 = &f( x,   y,   z   ); \
  fx = &f( x-1, y,   z   ); \
  fy = &f( x,   y-1, z   ); \
  fz = &f( x,   y,   z-1 )

#define NEXT_STENCIL()                        \
  f0++; fx++; fy++; fz++; x++;                \
  if ( x > nx )                               \
  {                                           \
                  y++;               x = 2;   \
    if ( y > ny ) z++; if ( y > ny ) y = 2;   \
    INIT_STENCIL();                           \
  }

#define UPDATE_EX()                                     \
  f0->tcax = ( py * ( f0->cbz * m[f0->fmatz].rmuz -     \
		      fy->cbz * m[fy->fmatz].rmuz ) -   \
               pz * ( f0->cby * m[f0->fmaty].rmuy -     \
		      fz->cby * m[fz->fmaty].rmuy ) )

#define UPDATE_EY()                                     \
  f0->tcay = ( pz * ( f0->cbx * m[f0->fmatx].rmux -     \
		      fz->cbx * m[fz->fmatx].rmux ) -   \
               px * ( f0->cbz * m[f0->fmatz].rmuz -     \
		      fx->cbz * m[fx->fmatz].rmuz ) )

#define UPDATE_EZ()                                     \
  f0->tcaz = ( px * ( f0->cby * m[f0->fmaty].rmuy -     \
		      fx->cby * m[fx->fmaty].rmuy ) -   \
               py * ( f0->cbx * m[f0->fmatx].rmux -     \
		      fy->cbx * m[fy->fmatx].rmux ) )

void
compute_curl_b_pipeline_scalar( pipeline_args_t * args,
                                int pipeline_rank,
                                int n_pipeline );

void
compute_curl_b_pipeline_v4( pipeline_args_t * args,
                            int pipeline_rank,
                            int n_pipeline );

void
compute_curl_b_pipeline_v8( pipeline_args_t * args,
                            int pipeline_rank,
                            int n_pipeline );

void
compute_curl_b_pipeline_v16( pipeline_args_t * args,
                             int pipeline_rank,
                             int n_pipeline );

#endif // _compute_curl_b_pipeline_h_
