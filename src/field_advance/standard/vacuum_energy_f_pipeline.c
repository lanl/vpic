// FIXME: USE THE DISCRETIZED VARIATIONAL DEFINITION OF ENERGY

#define IN_sfa

#include "sfa_private.h"

#include "../../util/pipelines/pipelines_exec.h"

typedef struct pipeline_args
{
  const field_t      * ALIGNED(128) f;
  const sfa_params_t *              p;
  const grid_t       *              g;
  double en[MAX_PIPELINE+1][6];
} pipeline_args_t;

#define DECLARE_STENCIL()                                                  \
  const field_t                * ALIGNED(128) f = args->f;                 \
  const material_coefficient_t * ALIGNED(128) m = args->p->mc;             \
  const grid_t                 *              g = args->g;                 \
  const int nx = g->nx, ny = g->ny, nz = g->nz;                            \
                                                                           \
  const float qepsx = 0.25*m->epsx;                                        \
  const float qepsy = 0.25*m->epsy;                                        \
  const float qepsz = 0.25*m->epsz;                                        \
  const float hrmux = 0.50*m->rmux; /* was previously 0.25 in master */	   \
  const float hrmuy = 0.50*m->rmuy; /* was previously 0.25 in master */	   \
  const float hrmuz = 0.50*m->rmuz; /* was previously 0.25 in master */	   \
                                                                           \
  const field_t * ALIGNED(16) f0;                                          \
  const field_t * ALIGNED(16) fx,  * ALIGNED(16) fy,  * ALIGNED(16) fz;    \
  const field_t * ALIGNED(16) fyz, * ALIGNED(16) fzx, * ALIGNED(16) fxy;   \
  double en_ex = 0, en_ey = 0, en_ez = 0, en_bx = 0, en_by = 0, en_bz = 0; \
  int x, y, z

#define f(x,y,z) f[ VOXEL(x,y,z, nx,ny,nz) ]

#define INIT_STENCIL()   \
  f0  = &f(x,  y,  z  ); \
  fx  = &f(x+1,y,  z  ); \
  fy  = &f(x,  y+1,z  ); \
  fz  = &f(x,  y,  z+1); \
  fyz = &f(x,  y+1,z+1); \
  fzx = &f(x+1,y,  z+1); \
  fxy = &f(x+1,y+1,z  )

#define NEXT_STENCIL()                              \
  f0++; fx++; fy++; fz++; fyz++; fzx++; fxy++; x++; \
  if( x>nx ) {                                      \
    /**/       y++;            x = 1;               \
    if( y>ny ) z++; if( y>ny ) y = 1;               \
    INIT_STENCIL();                                 \
  }

#define REDUCE_EN()                     \
  en_ex += qepsx*(  f0->ex * f0->ex +   \
                    fy->ex * fy->ex +   \
                    fz->ex * fz->ex +   \
                   fyz->ex *fyz->ex );  \
  en_ey += qepsy*(  f0->ey * f0->ey +   \
                    fz->ey * fz->ey +   \
                    fx->ey * fx->ey +   \
                   fzx->ey *fzx->ey );  \
  en_ez += qepsz*(  f0->ez * f0->ez +   \
                    fx->ez * fx->ez +   \
                    fy->ez * fy->ez +   \
                   fxy->ez *fxy->ez );  \
  en_bx += hrmux*(  f0->cbx* f0->cbx +  \
                    fx->cbx* fx->cbx ); \
  en_by += hrmuy*(  f0->cby* f0->cby +  \
                    fy->cby* fy->cby ); \
  en_bz += hrmuz*(  f0->cbz* f0->cbz +  \
                    fz->cbz* fz->cbz )
 
void
vacuum_energy_f_pipeline_scalar( pipeline_args_t * args,
                                 int pipeline_rank,
                                 int n_pipeline )
{
  DECLARE_STENCIL();
  
  int n_voxel;

  DISTRIBUTE_VOXELS( 1,nx, 1,ny, 1,nz, 16,
                     pipeline_rank, n_pipeline,
                     x, y, z, n_voxel );
  
  INIT_STENCIL();
  for( ; n_voxel; n_voxel-- )
  {
    REDUCE_EN();

    NEXT_STENCIL();
  }

  args->en[pipeline_rank][0] = en_ex;
  args->en[pipeline_rank][1] = en_ey;
  args->en[pipeline_rank][2] = en_ez;
  args->en[pipeline_rank][3] = en_bx;
  args->en[pipeline_rank][4] = en_by;
  args->en[pipeline_rank][5] = en_bz;
}

#if defined(V4_ACCELERATION) && defined(HAS_V4_PIPELINE)

#error "Not implemented"

#endif

void
vacuum_energy_f_pipeline( double * global,
                          const field_array_t * RESTRICT fa )
{
  if ( !global || !fa )
  {
    ERROR( ( "Bad args" ) );
  }

  // Have each pipeline and the host handle a portion of the
  // local voxels
  
  pipeline_args_t args[1];

  args->f = fa->f;
  args->p = (sfa_params_t *)fa->params;
  args->g = fa->g;

  EXEC_PIPELINES( vacuum_energy_f, args, 0 );

  WAIT_PIPELINES();

  // Reduce results from each pipelines
  
  int p;

  for( p=1; p<=N_PIPELINE; p++ )
  {
    args->en[0][0] += args->en[p][0]; args->en[0][1] += args->en[p][1];
    args->en[0][2] += args->en[p][2]; args->en[0][3] += args->en[p][3];
    args->en[0][4] += args->en[p][4]; args->en[0][5] += args->en[p][5];
  }
    
  // Convert to physical units and reduce results between nodes
  
  double v0 = 0.5*fa->g->eps0*fa->g->dV;

  args->en[0][0] *= v0; args->en[0][1] *= v0;
  args->en[0][2] *= v0; args->en[0][3] *= v0;
  args->en[0][4] *= v0; args->en[0][5] *= v0;

  // Reduce results between nodes

  mp_allsum_d( args->en[0], global, 6 );
}
