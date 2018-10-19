#define IN_sfa

#include "sfa_private.h"

#if defined(V4_ACCELERATION)

using namespace v4;

void
advance_e_pipeline_v4( pipeline_args_t * args,
                       int pipeline_rank,
                       int n_pipeline )
{
  DECLARE_STENCIL();

  int n_voxel;

  DISTRIBUTE_VOXELS( 2,nx, 2,ny, 2,nz, 16,
                     pipeline_rank, n_pipeline,
                     x, y, z, n_voxel );

  const v4float vdamp( damp );
  const v4float vpx( px );
  const v4float vpy( py );
  const v4float vpz( pz );
  const v4float vcj( cj );

  v4float save0, save1, dummy;

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

  field_t * ALIGNED(16) f00, * ALIGNED(16) f01, * ALIGNED(16) f02, * ALIGNED(16) f03; // Voxel quad

  field_t * ALIGNED(16) fx0, * ALIGNED(16) fx1, * ALIGNED(16) fx2, * ALIGNED(16) fx3; // Voxel quad +x neighbors

  field_t * ALIGNED(16) fy0, * ALIGNED(16) fy1, * ALIGNED(16) fy2, * ALIGNED(16) fy3; // Voxel quad +y neighbors

  field_t * ALIGNED(16) fz0, * ALIGNED(16) fz1, * ALIGNED(16) fz2, * ALIGNED(16) fz3; // Voxel quad +z neighbors

  // Process the bulk of the voxels 4 at a time
                               
  INIT_STENCIL();

  for( ; n_voxel > 3; n_voxel -= 4 )
  {
    f00 = f0; fx0 = fx; fy0 = fy; fz0 = fz; NEXT_STENCIL();
    f01 = f0; fx1 = fx; fy1 = fy; fz1 = fz; NEXT_STENCIL();
    f02 = f0; fx2 = fx; fy2 = fy; fz2 = fz; NEXT_STENCIL();
    f03 = f0; fx3 = fx; fy3 = fy; fz3 = fz; NEXT_STENCIL();

    //------------------------------------------------------------------------//
    // Load field data.
    //------------------------------------------------------------------------//

    load_4x4_tr( &f00->ex, &f01->ex, &f02->ex, &f03->ex,
                 f0_ex, f0_ey, f0_ez, save0 );

    load_4x3_tr( &f00->cbx, &f01->cbx, &f02->cbx, &f03->cbx,
                 f0_cbx, f0_cby, f0_cbz );

    load_4x4_tr( &f00->tcax, &f01->tcax, &f02->tcax, &f03->tcax,
                 f0_tcax, f0_tcay, f0_tcaz, save1 );

    load_4x3_tr( &f00->jfx, &f01->jfx, &f02->jfx, &f03->jfx,
                 f0_jfx, f0_jfy, f0_jfz );

    load_4x3_tr( &fx0->cbx, &fx1->cbx, &fx2->cbx, &fx3->cbx,
                 dummy, fx_cby, fx_cbz );

    load_4x3_tr( &fy0->cbx, &fy1->cbx, &fy2->cbx, &fy3->cbx,
                 fy_cbx, dummy, fy_cbz );

    load_4x2_tr( &fz0->cbx, &fz1->cbx, &fz2->cbx, &fz3->cbx,
                 fz_cbx, fz_cby );

#   define LOAD_RMU(V,D) m_f##V##_rmu##D=v4float( m[f##V##0->fmat##D].rmu##D, \
                                                  m[f##V##1->fmat##D].rmu##D, \
                                                  m[f##V##2->fmat##D].rmu##D, \
                                                  m[f##V##3->fmat##D].rmu##D )

    LOAD_RMU(0,x); LOAD_RMU(0,y); LOAD_RMU(0,z);
                   LOAD_RMU(x,y); LOAD_RMU(x,z);
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

#   undef LOAD_RMU

    f0_cbx_rmux = f0_cbx * m_f0_rmux;
    f0_cby_rmuy = f0_cby * m_f0_rmuy;
    f0_cbz_rmuz = f0_cbz * m_f0_rmuz;

    f0_tcax = fnms( vdamp,
                    f0_tcax,
                    fms( vpy,
                         fnms( fy_cbz, m_fy_rmuz, f0_cbz_rmuz ),
                         vpz * fnms( fz_cby, m_fz_rmuy, f0_cby_rmuy ) ) );

    f0_tcay = fnms( vdamp,
                    f0_tcay,
                    fms( vpz,
                         fnms( fz_cbx, m_fz_rmux, f0_cbx_rmux ),
                         vpx * fnms( fx_cbz, m_fx_rmuz, f0_cbz_rmuz ) ) );

    f0_tcaz = fnms( vdamp,
                    f0_tcaz,
                    fms( vpx,
                         fnms( fx_cby, m_fx_rmuy, f0_cby_rmuy ),
                         vpy * fnms( fy_cbx, m_fy_rmux, f0_cbx_rmux ) ) );

    f0_ex = fma( m_f0_decayx, f0_ex, m_f0_drivex * fnms( vcj, f0_jfx, f0_tcax ) );
    f0_ey = fma( m_f0_decayy, f0_ey, m_f0_drivey * fnms( vcj, f0_jfy, f0_tcay ) );
    f0_ez = fma( m_f0_decayz, f0_ez, m_f0_drivez * fnms( vcj, f0_jfz, f0_tcaz ) );

    //------------------------------------------------------------------------//
    // Note:
    //------------------------------------------------------------------------//
    // Unlike load_4x3 versus load_4x4, store_4x4 is much more efficient than
    // store_4x3.
    //------------------------------------------------------------------------//

    store_4x4_tr( f0_ex, f0_ey, f0_ez, save0,
                  &f00->ex, &f01->ex, &f02->ex, &f03->ex );

    store_4x4_tr( f0_tcax, f0_tcay, f0_tcaz, save1,
                  &f00->tcax, &f01->tcax, &f02->tcax, &f03->tcax );
  }
}

#else

void
advance_e_pipeline_v4( pipeline_args_t * args,
                       int pipeline_rank,
                       int n_pipeline )
{
  // No v4 implementation.
  ERROR( ( "No advance_e_pipeline_v4 implementation." ) );
}

#endif
