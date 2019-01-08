#ifndef _vacuum_clean_div_e_pipeline_h_
#define _vacuum_clean_div_e_pipeline_h_

#ifndef IN_vacuum_clean_div_e_pipeline
#error "Only include vacuum_clean_div_e_pipeline.h in vacuum_clean_div_e_pipeline source files."
#endif

#include "../sfa_private.h"

typedef struct pipeline_args
{
  field_t            * ALIGNED(128) f;
  const sfa_params_t *              p;
  const grid_t       *              g;
} pipeline_args_t;

#define DECLARE_STENCIL()                                                \
  field_t                      * ALIGNED(128) f = args->f;               \
  const material_coefficient_t * ALIGNED(128) m = args->p->mc;           \
  const grid_t                 *              g = args->g;               \
  const int nx = g->nx, ny = g->ny, nz = g->nz;                          \
                                                                         \
  const float _rdx = (nx>1) ? g->rdx : 0;                                \
  const float _rdy = (ny>1) ? g->rdy : 0;                                \
  const float _rdz = (nz>1) ? g->rdz : 0;                                \
  const float alphadt = 0.3888889/( _rdx*_rdx + _rdy*_rdy + _rdz*_rdz ); \
  const float px   = (alphadt*_rdx)*m->drivex;                           \
  const float py   = (alphadt*_rdy)*m->drivey;                           \
  const float pz   = (alphadt*_rdz)*m->drivez;                           \
                                                                         \
  field_t * ALIGNED(16) f0;                                              \
  field_t * ALIGNED(16) fx, * ALIGNED(16) fy, * ALIGNED(16) fz;          \
  int x, y, z
                     
#define f(x,y,z) f[ VOXEL(x,y,z,nx,ny,nz) ]

#define INIT_STENCIL()  \
  f0 = &f(x,  y,  z  ); \
  fx = &f(x+1,y,  z  ); \
  fy = &f(x,  y+1,z  ); \
  fz = &f(x,  y,  z+1)

#define NEXT_STENCIL()                \
  f0++; fx++; fy++; fz++; x++;        \
  if( x>nx ) {                        \
    /**/       y++;            x = 1; \
    if( y>ny ) z++; if( y>ny ) y = 1; \
    INIT_STENCIL();                   \
  }

#define MARDER_EX() f0->ex += px*(fx->div_e_err-f0->div_e_err)
#define MARDER_EY() f0->ey += py*(fy->div_e_err-f0->div_e_err)
#define MARDER_EZ() f0->ez += pz*(fz->div_e_err-f0->div_e_err)

static void
vacuum_clean_div_e_pipeline_scalar( pipeline_args_t * args,
                                    int pipeline_rank,
                                    int n_pipeline );

#endif // _vacuum_clean_div_e_pipeline_h_
