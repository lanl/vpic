#define IN_sfa
#define IN_advance_b_pipeline

#include "advance_b_pipeline.h"
#include "sfa_private.h"

#if defined(V16_ACCELERATION)

using namespace v16;

void
advance_b_pipeline_v16( pipeline_args_t * args,
                        int pipeline_rank,
                        int n_pipeline )
{
  DECLARE_STENCIL();

  int n_voxel;

  DISTRIBUTE_VOXELS( 1,nx, 1,ny, 1,nz, 16,
                     pipeline_rank, n_pipeline,
                     x, y, z, n_voxel );

  const v16float vpx( px );
  const v16float vpy( py );
  const v16float vpz( pz );

  v16float f0_ex,  f0_ey,  f0_ez;  // Voxel block electric fields
  v16float f0_cbx, f0_cby, f0_cbz; // Voxel block magnetic fields
  v16float fx_ey, fx_ez;           // Voxel block +x neighbor fields
  v16float fy_ez, fy_ex;           // Voxel block +y neighbor fields
  v16float fz_ex, fz_ey;           // Voxel block +z neighbor fields
  v16float dummy;

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

    load_16x3_tr( &f000->ex, &f001->ex, &f002->ex, &f003->ex,
                  &f004->ex, &f005->ex, &f006->ex, &f007->ex,
                  &f008->ex, &f009->ex, &f010->ex, &f011->ex,
                  &f012->ex, &f013->ex, &f014->ex, &f015->ex,
                  f0_ex, f0_ey, f0_ez );

    load_16x3_tr( &f000->cbx, &f001->cbx, &f002->cbx, &f003->cbx,
                  &f004->cbx, &f005->cbx, &f006->cbx, &f007->cbx,
                  &f008->cbx, &f009->cbx, &f010->cbx, &f011->cbx,
                  &f012->cbx, &f013->cbx, &f014->cbx, &f015->cbx,
                  f0_cbx, f0_cby, f0_cbz );

    load_16x3_tr( &fx00->ex, &fx01->ex, &fx02->ex, &fx03->ex,
                  &fx04->ex, &fx05->ex, &fx06->ex, &fx07->ex,
                  &fx08->ex, &fx09->ex, &fx10->ex, &fx11->ex,
                  &fx12->ex, &fx13->ex, &fx14->ex, &fx15->ex,
                  dummy, fx_ey, fx_ez );

    load_16x3_tr( &fy00->ex, &fy01->ex, &fy02->ex, &fy03->ex,
                  &fy04->ex, &fy05->ex, &fy06->ex, &fy07->ex,
                  &fy08->ex, &fy09->ex, &fy10->ex, &fy11->ex,
                  &fy12->ex, &fy13->ex, &fy14->ex, &fy15->ex,
                  fy_ex, dummy, fy_ez );

    load_16x2_tr( &fz00->ex, &fz01->ex, &fz02->ex, &fz03->ex,
                  &fz04->ex, &fz05->ex, &fz06->ex, &fz07->ex,
                  &fz08->ex, &fz09->ex, &fz10->ex, &fz11->ex,
                  &fz12->ex, &fz13->ex, &fz14->ex, &fz15->ex,
                  fz_ex, fz_ey );

    f0_cbx += fnms( vpy, ( fy_ez - f0_ez ), vpz*( fz_ey - f0_ey ) );
    f0_cby += fnms( vpz, ( fz_ex - f0_ex ), vpx*( fx_ez - f0_ez ) );
    f0_cbz += fnms( vpx, ( fx_ey - f0_ey ), vpy*( fy_ex - f0_ex ) );

    store_16x3_tr( f0_cbx, f0_cby, f0_cbz,
                   &f000->cbx, &f001->cbx, &f002->cbx, &f003->cbx,
                   &f004->cbx, &f005->cbx, &f006->cbx, &f007->cbx,
                   &f008->cbx, &f009->cbx, &f010->cbx, &f011->cbx,
                   &f012->cbx, &f013->cbx, &f014->cbx, &f015->cbx );
  }
}

#else

void
advance_b_pipeline_v16( pipeline_args_t * args,
                        int pipeline_rank,
                        int n_pipeline )
{
  // No v16 implementation.
  ERROR( ( "No advance_b_pipeline_v16 implementation." ) );
}

#endif
