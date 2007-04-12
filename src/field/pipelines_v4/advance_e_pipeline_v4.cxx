#include <field_pipelines.h>

#ifdef V4_ACCELERATION
using namespace v4;

// Note: This is similar to compute_curl_b_pipeline

#if 0 // Original non-pipelined non-vectorized version 
  for( z=2; z<=nz; z++ ) {
    for( y=2; y<=ny; y++ ) {
      f0 = &f(2,y,  z);
      fx = &f(1,y,  z);
      fy = &f(2,y-1,z);
      fz = &f(2,y,  z-1);
      for( x=2; x<=nx; x++ ) {

        f0_cbx_rmux = f0->cbx*m[f0->fmatx].rmux;
        f0_cby_rmuy = f0->cby*m[f0->fmaty].rmuy;
        f0_cbz_rmuz = f0->cbz*m[f0->fmatz].rmuz;

        f0->tcax = ( py*( f0_cbz_rmuz - fy->cbz*m[fy->fmatz].rmuz ) -   
                     pz*( f0_cby_rmuy - fz->cby*m[fz->fmaty].rmuy ) ) -
                   damp*f0->tcax; 

        f0->tcay = ( pz*( f0_cbx_rmux - fz->cbx*m[fz->fmatx].rmux ) - 
                     px*( f0_cbz_rmuz - fx->cbz*m[fx->fmatz].rmuz ) ) -
                   damp*f0->tcay;

        f0->tcaz = ( px*( f0_cby_rmuy - fx->cby*m[fx->fmaty].rmuy ) - 
                     py*( f0_cbx_rmux - fy->cbx*m[fy->fmatx].rmux ) ) -
                   damp*f0->tcaz; 

        f0->ex = m[f0->ematx].decayx*f0->ex +
                 m[f0->ematx].drivex*( f0->tcax - cj*f0->jfx );

        f0->ey = m[f0->ematy].decayy*f0->ey +
                 m[f0->ematy].drivey*( f0->tcay - cj*f0->jfy );

        f0->ez = m[f0->ematz].decayz*f0->ez +
                 m[f0->ematz].drivez*( f0->tcaz - cj*f0->jfz );

	f0++; fx++; fy++; fz++;

      }
    }
  }
#endif

#define f(x,y,z) f[INDEX_FORTRAN_3(x,y,z,0,nx+1,0,ny+1,0,nz+1)]

void
advance_e_pipeline_v4( advance_e_pipeline_args_t * args,
                       int pipeline_rank,
                       int n_pipeline ) {
  field_t                      * ALIGNED f = args->f;
  const material_coefficient_t * ALIGNED m = args->m;
  const grid_t                 *         g = args->g;

  int x, y, z, n_voxel;
  field_t *f0, *fx, *fy, *fz;

  const int nx = g->nx;
  const int ny = g->ny;
  const int nz = g->nz;

  const float damp = g->damp;
  const float px = (nx>1) ? (1+damp)*g->cvac*g->dt/g->dx : 0;
  const float py = (ny>1) ? (1+damp)*g->cvac*g->dt/g->dy : 0;
  const float pz = (nz>1) ? (1+damp)*g->cvac*g->dt/g->dz : 0;
  const float cj = g->dt/g->eps0;

  float f00_cbx_rmux, f00_cby_rmuy, f00_cbz_rmuz;

  const v4float vdamp(damp);
  const v4float vpx(px);
  const v4float vpy(py);
  const v4float vpz(pz);
  const v4float vcj(cj);

  v4float f0_ex,   f0_ey,   f0_ez;
  v4float f0_cbx,  f0_cby,  f0_cbz;
  v4float f0_tcax, f0_tcay, f0_tcaz;
  v4float f0_jfx,  f0_jfy,  f0_jfz;
  v4float          fx_cby,  fx_cbz;
  v4float fy_cbx,           fy_cbz;
  v4float fz_cbx,  fz_cby;
  v4float m_f0_rmux, m_f0_rmuy, m_f0_rmuz;
  v4float            m_fx_rmuy, m_fx_rmuz;
  v4float m_fy_rmux,            m_fy_rmuz;
  v4float m_fz_rmux, m_fz_rmuy;
  v4float m_f0_decayx, m_f0_drivex;
  v4float m_f0_decayy, m_f0_drivey;
  v4float m_f0_decayz, m_f0_drivez;

  v4float f0_cbx_rmux, f0_cby_rmuy, f0_cbz_rmuz;

  field_t * ALIGNED f00, * ALIGNED f01, * ALIGNED f02, * ALIGNED f03; // Voxel quad
  field_t * ALIGNED fx0, * ALIGNED fx1, * ALIGNED fx2, * ALIGNED fx3; // Voxel quad +x neighbors
  field_t * ALIGNED fy0, * ALIGNED fy1, * ALIGNED fy2, * ALIGNED fy3; // Voxel quad +y neighbors
  field_t * ALIGNED fz0, * ALIGNED fz1, * ALIGNED fz2, * ALIGNED fz3; // Voxel quad +z neighbors

  // Process the voxels assigned to this pipeline 
  
  n_voxel = distribute_voxels_v4( 2,nx, 2,ny, 2,nz,
                                  pipeline_rank, n_pipeline,
                                  &x, &y, &z );

  // Process the bulk of the voxels 4 at a time
                               
# define LOAD_PTRS()    \
  f0 = &f(x,  y,  z  ); \
  fx = &f(x-1,y,  z  ); \
  fy = &f(x,  y-1,z  ); \
  fz = &f(x,  y,  z-1)

# define ADVANCE_PTRS(n) \
  f0##n = f0++;          \
  fx##n = fx++;          \
  fy##n = fy++;          \
  fz##n = fz++;          \
  x++;                   \
  if( x>nx ) {           \
    x=2, y++;            \
    if( y>ny ) y=2, z++; \
    LOAD_PTRS();         \
  }

  LOAD_PTRS();
  ADVANCE_PTRS(0); ADVANCE_PTRS(1); ADVANCE_PTRS(2); ADVANCE_PTRS(3);

  for( ; n_voxel>3; n_voxel-=4 ) {

    load_4x4_tr( &f00->ex,    &f01->ex,    &f02->ex,    &f03->ex,    f0_ex,    f0_ey,    f0_ez,    f0_cbx );
    load_4x2_tr( &f00->cby,   &f01->cby,   &f02->cby,   &f03->cby,   f0_cby,   f0_cbz                     );
    load_4x4_tr( &f00->tcaz,  &f01->tcaz,  &f02->tcaz,  &f03->tcaz,  f0_tcaz,  f0_jfx,   f0_jfy,   f0_jfz );

    load_4x2_tr( &fx0->cby,   &fx1->cby,   &fx2->cby,   &fx3->cby,   fx_cby,   fx_cbz                     );

    load_4x1_tr( &fy0->cbx,   &fy1->cbx,   &fy2->cbx,   &fy3->cbx,                                 fy_cbx );
    load_4x1_tr( &fy0->cbz,   &fy1->cbz,   &fy2->cbz,   &fy3->cbz,             fy_cbz                     );

    load_4x1_tr( &fz0->cbx,   &fz1->cbx,   &fz2->cbx,   &fz3->cbx,                                 fy_cbx );
    load_4x1_tr( &fz0->cby,   &fz1->cby,   &fz2->cby,   &fz3->cby,   fz_cby                               );

#   define LOAD_RMU( V, D ) load_4x1_tr( &m[f##V##0->fmat##D].rmu##D, \
                                         &m[f##V##1->fmat##D].rmu##D, \
                                         &m[f##V##2->fmat##D].rmu##D, \
                                         &m[f##V##3->fmat##D].rmu##D, \
                                         m_f##V##_rmu##D )

    LOAD_RMU(0,x); LOAD_RMU(0,y); LOAD_RMU(0,z);
    /**/           LOAD_RMU(x,y); LOAD_RMU(x,z);
    LOAD_RMU(y,x);                LOAD_RMU(y,z);
    LOAD_RMU(z,x); LOAD_RMU(z,y);
    
    load_4x2_tr( &m[f00->ematx].decayx, &m[f01->ematx].decayx,
                 &m[f02->ematx].decayx, &m[f03->ematx].decayx,
                 m_f0_decayx, m_f0_drivex );
    load_4x2_tr( &m[f00->ematy].decayy, &m[f01->ematy].decayy,
                 &m[f02->ematy].decayy, &m[f03->ematy].decayy,
                 m_f0_decayy, m_f0_drivey );
    load_4x2_tr( &m[f00->ematz].decayz, &m[f01->ematz].decayz,
                 &m[f02->ematz].decayz, &m[f03->ematz].decayz,
                 m_f0_decayz, m_f0_drivez );

    f0_cbx_rmux = f0_cbx * m_f0_rmux;
    f0_cby_rmuy = f0_cby * m_f0_rmuy;
    f0_cbz_rmuz = f0_cbz * m_f0_rmuz;

    f0_tcax = fnms( vdamp,f0_tcax,
                    fms( vpy,fnms( fy_cbz,m_fy_rmuz, f0_cbz_rmuz ),
                         vpz*fnms( fz_cby,m_fz_rmuy, f0_cby_rmuy ) ) );

    f0_tcay = fnms( vdamp,f0_tcay,
                    fms( vpz,fnms( fz_cbx,m_fz_rmux, f0_cbx_rmux ),
                         vpx*fnms( fx_cbz,m_fx_rmuz, f0_cbz_rmuz ) ) );

    f0_tcaz = fnms( vdamp,f0_tcaz,
                    fms( vpx,fnms( fx_cby,m_fx_rmuy, f0_cby_rmuy ),
                         vpy*fnms( fy_cbx,m_fy_rmux, f0_cbx_rmux ) ) );

    f0_ex   = fms( m_f0_decayx,f0_ex,
                   m_f0_drivex*fms( vcj,f0_jfx, f0_tcax ) );

    f0_ey   = fms( m_f0_decayy,f0_ey,
                   m_f0_drivey*fms( vcj,f0_jfy, f0_tcay ) );

    f0_ez   = fms( m_f0_decayz,f0_ez,
                   m_f0_drivez*fms( vcj,f0_jfz, f0_tcaz ) );

    store_4x4_tr( f0_ex, f0_ey, f0_ez,   f0_cbx,  &f00->ex,    &f01->ex,    &f02->ex,    &f03->ex   ); // Could use store_4x3_tr
    store_4x2_tr(               f0_tcax, f0_tcay, &f00->tcax,  &f01->tcax,  &f02->tcax,  &f03->tcax );
    store_4x1_tr( f0_tcaz,                        &f00->tcaz,  &f01->tcaz,  &f02->tcaz,  &f03->tcaz );

    ADVANCE_PTRS(0); ADVANCE_PTRS(1); ADVANCE_PTRS(2); ADVANCE_PTRS(3);

  }

  // Process straggler voxels with scalar processing
  // Note that the pointers needed for straggler voxel processing were already
  // loaded implicitly above.

  // FIXME: CONSIDER A MODEL WHERE THE HOST PROCESSES STRAGGLERS.

  for( ; n_voxel; n_voxel-- ) {

    f00_cbx_rmux = f00->cbx*m[f00->fmatx].rmux;
    f00_cby_rmuy = f00->cby*m[f00->fmaty].rmuy;
    f00_cbz_rmuz = f00->cbz*m[f00->fmatz].rmuz;

    f00->tcax = ( py*( f00_cbz_rmuz - fy0->cbz*m[fy0->fmatz].rmuz ) -
                  pz*( f00_cby_rmuy - fz0->cby*m[fz0->fmaty].rmuy ) ) -
                damp*f00->tcax; 

    f00->tcay = ( pz*( f00_cbx_rmux - fz0->cbx*m[fz0->fmatx].rmux ) -
                  px*( f00_cbz_rmuz - fx0->cbz*m[fx0->fmatz].rmuz ) ) -
                damp*f00->tcay;

    f00->tcaz = ( px*( f00_cby_rmuy - fx0->cby*m[fx0->fmaty].rmuy ) -
                  py*( f00_cbx_rmux - fy0->cbx*m[fy0->fmatx].rmux ) ) -
                damp*f00->tcaz; 

    f00->ex = m[f00->ematx].decayx*f00->ex +
              m[f00->ematx].drivex*( f00->tcax - cj*f00->jfx );

    f00->ey = m[f00->ematy].decayy*f00->ey +
              m[f00->ematy].drivey*( f00->tcay - cj*f00->jfy );

    f00->ez = m[f00->ematz].decayz*f00->ez +
              m[f00->ematz].drivez*( f00->tcaz - cj*f00->jfz );

    f00 = f01; fx0 = fx1; fy0 = fy1; fz0 = fz1;
    f01 = f02; fx1 = fx2; fy1 = fy2; fz1 = fz2;

  }

}

#endif
