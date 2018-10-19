#define IN_sfa

#include "compute_curl_b_pipeline.h"
#include "sfa_private.h"

#if defined(V16_ACCELERATION)

using namespace v16;

void
compute_curl_b_pipeline_v16( pipeline_args_t * args,
                             int pipeline_rank,
                             int n_pipeline )
{
  DECLARE_STENCIL();

  int n_voxel;

  DISTRIBUTE_VOXELS( 2,nx, 2,ny, 2,nz, 16,
                     pipeline_rank, n_pipeline,
                     x, y, z, n_voxel );

  const v16float vpx( px );
  const v16float vpy( py );
  const v16float vpz( pz );

  v16float save1, dummy;

  v16float f0_cbx,  f0_cby,  f0_cbz;
  v16float f0_tcax, f0_tcay, f0_tcaz;
  v16float          fx_cby,  fx_cbz;
  v16float fy_cbx,           fy_cbz;
  v16float fz_cbx,  fz_cby;
  v16float m_f0_rmux, m_f0_rmuy, m_f0_rmuz;
  v16float            m_fx_rmuy, m_fx_rmuz;
  v16float m_fy_rmux,            m_fy_rmuz;
  v16float m_fz_rmux, m_fz_rmuy;

  v16float f0_cbx_rmux, f0_cby_rmuy, f0_cbz_rmuz;

  field_t * ALIGNED(16) f000, * ALIGNED(16) f001, * ALIGNED(16) f002, * ALIGNED(16) f003; // Voxel block
  field_t * ALIGNED(16) f004, * ALIGNED(16) f005, * ALIGNED(16) f006, * ALIGNED(16) f007; // Voxel block
  field_t * ALIGNED(16) f008, * ALIGNED(16) f009, * ALIGNED(16) f010, * ALIGNED(16) f011; // Voxel block
  field_t * ALIGNED(16) f012, * ALIGNED(16) f013, * ALIGNED(16) f014, * ALIGNED(16) f015; // Voxel block

  field_t * ALIGNED(16) fx00, * ALIGNED(16) fx01, * ALIGNED(16) fx02, * ALIGNED(16) fx03; // Voxel block +x neighbors
  field_t * ALIGNED(16) fx04, * ALIGNED(16) fx05, * ALIGNED(16) fx06, * ALIGNED(16) fx07; // Voxel block +x neighbors
  field_t * ALIGNED(16) fx08, * ALIGNED(16) fx09, * ALIGNED(16) fx10, * ALIGNED(16) fx11; // Voxel block +x neighbors
  field_t * ALIGNED(16) fx12, * ALIGNED(16) fx13, * ALIGNED(16) fx14, * ALIGNED(16) fx15; // Voxel block +x neighbors

  field_t * ALIGNED(16) fy00, * ALIGNED(16) fy01, * ALIGNED(16) fy02, * ALIGNED(16) fy03; // Voxel block +y neighbors
  field_t * ALIGNED(16) fy04, * ALIGNED(16) fy05, * ALIGNED(16) fy06, * ALIGNED(16) fy07; // Voxel block +y neighbors
  field_t * ALIGNED(16) fy08, * ALIGNED(16) fy09, * ALIGNED(16) fy10, * ALIGNED(16) fy11; // Voxel block +y neighbors
  field_t * ALIGNED(16) fy12, * ALIGNED(16) fy13, * ALIGNED(16) fy14, * ALIGNED(16) fy15; // Voxel block +y neighbors

  field_t * ALIGNED(16) fz00, * ALIGNED(16) fz01, * ALIGNED(16) fz02, * ALIGNED(16) fz03; // Voxel block +z neighbors
  field_t * ALIGNED(16) fz04, * ALIGNED(16) fz05, * ALIGNED(16) fz06, * ALIGNED(16) fz07; // Voxel block +z neighbors
  field_t * ALIGNED(16) fz08, * ALIGNED(16) fz09, * ALIGNED(16) fz10, * ALIGNED(16) fz11; // Voxel block +z neighbors
  field_t * ALIGNED(16) fz12, * ALIGNED(16) fz13, * ALIGNED(16) fz14, * ALIGNED(16) fz15; // Voxel block +z neighbors

  // Process the bulk of the voxels 16 at a time

  INIT_STENCIL();

  for( ; n_voxel > 15; n_voxel -= 16 )
  {
    f000 = f0; fx00 = fx; fy00 = fy; fz00 = fz; NEXT_STENCIL();
    f001 = f0; fx01 = fx; fy01 = fy; fz01 = fz; NEXT_STENCIL();
    f002 = f0; fx02 = fx; fy02 = fy; fz02 = fz; NEXT_STENCIL();
    f003 = f0; fx03 = fx; fy03 = fy; fz03 = fz; NEXT_STENCIL();
    f004 = f0; fx04 = fx; fy04 = fy; fz04 = fz; NEXT_STENCIL();
    f005 = f0; fx05 = fx; fy05 = fy; fz05 = fz; NEXT_STENCIL();
    f006 = f0; fx06 = fx; fy06 = fy; fz06 = fz; NEXT_STENCIL();
    f007 = f0; fx07 = fx; fy07 = fy; fz07 = fz; NEXT_STENCIL();
    f008 = f0; fx08 = fx; fy08 = fy; fz08 = fz; NEXT_STENCIL();
    f009 = f0; fx09 = fx; fy09 = fy; fz09 = fz; NEXT_STENCIL();
    f010 = f0; fx10 = fx; fy10 = fy; fz10 = fz; NEXT_STENCIL();
    f011 = f0; fx11 = fx; fy11 = fy; fz11 = fz; NEXT_STENCIL();
    f012 = f0; fx12 = fx; fy12 = fy; fz12 = fz; NEXT_STENCIL();
    f013 = f0; fx13 = fx; fy13 = fy; fz13 = fz; NEXT_STENCIL();
    f014 = f0; fx14 = fx; fy14 = fy; fz14 = fz; NEXT_STENCIL();
    f015 = f0; fx15 = fx; fy15 = fy; fz15 = fz; NEXT_STENCIL();

    //------------------------------------------------------------------------//
    // Load field data.
    //------------------------------------------------------------------------//

    load_16x3_tr( &f000->cbx, &f001->cbx, &f002->cbx, &f003->cbx,
                  &f004->cbx, &f005->cbx, &f006->cbx, &f007->cbx,
                  &f008->cbx, &f009->cbx, &f010->cbx, &f011->cbx,
                  &f012->cbx, &f013->cbx, &f014->cbx, &f015->cbx,
                  f0_cbx, f0_cby, f0_cbz );

    load_16x4_tr( &f000->tcax, &f001->tcax, &f002->tcax, &f003->tcax,
                  &f004->tcax, &f005->tcax, &f006->tcax, &f007->tcax,
                  &f008->tcax, &f009->tcax, &f010->tcax, &f011->tcax,
                  &f012->tcax, &f013->tcax, &f014->tcax, &f015->tcax,
                  f0_tcax, f0_tcay, f0_tcaz, save1 );

    load_16x3_tr( &fx00->cbx, &fx01->cbx, &fx02->cbx, &fx03->cbx,
                  &fx04->cbx, &fx05->cbx, &fx06->cbx, &fx07->cbx,
                  &fx08->cbx, &fx09->cbx, &fx10->cbx, &fx11->cbx,
                  &fx12->cbx, &fx13->cbx, &fx14->cbx, &fx15->cbx,
                  dummy, fx_cby, fx_cbz );

    load_16x3_tr( &fy00->cbx, &fy01->cbx, &fy02->cbx, &fy03->cbx,
                  &fy04->cbx, &fy05->cbx, &fy06->cbx, &fy07->cbx,
                  &fy08->cbx, &fy09->cbx, &fy10->cbx, &fy11->cbx,
                  &fy12->cbx, &fy13->cbx, &fy14->cbx, &fy15->cbx,
                  fy_cbx, dummy, fy_cbz );

    load_16x2_tr( &fz00->cbx, &fz01->cbx, &fz02->cbx, &fz03->cbx,
                  &fz04->cbx, &fz05->cbx, &fz06->cbx, &fz07->cbx,
                  &fz08->cbx, &fz09->cbx, &fz10->cbx, &fz11->cbx,
                  &fz12->cbx, &fz13->cbx, &fz14->cbx, &fz15->cbx,
                  fz_cbx, fz_cby );

#   define LOAD_RMU(V,D) m_f##V##_rmu##D=v16float( m[f##V##00->fmat##D].rmu##D, \
                                                   m[f##V##01->fmat##D].rmu##D, \
                                                   m[f##V##02->fmat##D].rmu##D, \
                                                   m[f##V##03->fmat##D].rmu##D, \
                                                   m[f##V##04->fmat##D].rmu##D, \
                                                   m[f##V##05->fmat##D].rmu##D, \
                                                   m[f##V##06->fmat##D].rmu##D, \
                                                   m[f##V##07->fmat##D].rmu##D, \
                                                   m[f##V##08->fmat##D].rmu##D, \
                                                   m[f##V##09->fmat##D].rmu##D, \
                                                   m[f##V##10->fmat##D].rmu##D, \
                                                   m[f##V##11->fmat##D].rmu##D, \
                                                   m[f##V##12->fmat##D].rmu##D, \
                                                   m[f##V##13->fmat##D].rmu##D, \
                                                   m[f##V##14->fmat##D].rmu##D, \
                                                   m[f##V##15->fmat##D].rmu##D )

    LOAD_RMU(0,x); LOAD_RMU(0,y); LOAD_RMU(0,z);
                   LOAD_RMU(x,y); LOAD_RMU(x,z);
    LOAD_RMU(y,x);                LOAD_RMU(y,z);
    LOAD_RMU(z,x); LOAD_RMU(z,y);

#   undef LOAD_RMU

    f0_cbx_rmux = f0_cbx * m_f0_rmux;
    f0_cby_rmuy = f0_cby * m_f0_rmuy;
    f0_cbz_rmuz = f0_cbz * m_f0_rmuz;

    f0_tcax = fms( vpy,
                   fnms( fy_cbz, m_fy_rmuz, f0_cbz_rmuz ),
                   vpz * fnms( fz_cby, m_fz_rmuy, f0_cby_rmuy ) );

    f0_tcay = fms( vpz,
                   fnms( fz_cbx, m_fz_rmux, f0_cbx_rmux ),
                   vpx * fnms( fx_cbz, m_fx_rmuz, f0_cbz_rmuz ) );

    f0_tcaz = fms( vpx,
                   fnms( fx_cby, m_fx_rmuy, f0_cby_rmuy ),
                   vpy * fnms( fy_cbx, m_fy_rmux, f0_cbx_rmux ) );

    //------------------------------------------------------------------------//
    // Note:
    //------------------------------------------------------------------------//
    // Unlike load_16x3 versus load_16x4, store_16x4 is much more efficient
    // than store_16x3.
    //------------------------------------------------------------------------//

    store_16x4_tr( f0_tcax, f0_tcay, f0_tcaz, save1,
                   &f000->tcax, &f001->tcax, &f002->tcax, &f003->tcax,
                   &f004->tcax, &f005->tcax, &f006->tcax, &f007->tcax,
                   &f008->tcax, &f009->tcax, &f010->tcax, &f011->tcax,
                   &f012->tcax, &f013->tcax, &f014->tcax, &f015->tcax );
  }
}

#else

void
compute_curl_b_pipeline_v16( pipeline_args_t * args,
                             int pipeline_rank,
                             int n_pipeline )
{
  // No v16 implementation.
  ERROR( ( "No compute_curl_b_pipeline_v16 implementation." ) );
}

#endif
