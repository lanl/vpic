#define IN_sfa
#define IN_clean_div_b_pipeline

#include "clean_div_b_pipeline.h"

#include "../sfa_private.h"

#if defined(V16_ACCELERATION)

using namespace v16;

void
clean_div_b_pipeline_v16( pipeline_args_t * args,
                          int pipeline_rank,
                          int n_pipeline )
{
  field_t      * ALIGNED(128) f = args->f;
  const grid_t *              g = args->g;

  field_t * ALIGNED(16) f0;
  field_t * ALIGNED(16) fx, * ALIGNED(16) fy, * ALIGNED(16) fz;
  int x, y, z, n_voxel;
  
  const int nx = g->nx;
  const int ny = g->ny;
  const int nz = g->nz;

  float px, py, pz, alphadt;

  px = ( nx > 1 ) ? g->rdx : 0;
  py = ( ny > 1 ) ? g->rdy : 0;
  pz = ( nz > 1 ) ? g->rdz : 0;

  alphadt = 0.3888889/( px*px + py*py + pz*pz );

  px *= alphadt;
  py *= alphadt;
  pz *= alphadt;

  const v16float vpx(px);
  const v16float vpy(py);
  const v16float vpz(pz);

  v16float f0_cbx, f0_cby, f0_cbz; // Voxel block magnetic fields
  v16float f0_div_b_err;           // Voxel block div b errs
  v16float fx_div_b_err;           // Voxel block -x neighbor div b err
  v16float fy_div_b_err;           // Voxel block -y neighbor div b err
  v16float fz_div_b_err;           // Voxel block -z neighbor div b err

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

  // Process voxels assigned to this pipeline 
  
  DISTRIBUTE_VOXELS( 2,nx, 2,ny, 2,nz, 16,
                     pipeline_rank, n_pipeline,
                     x, y, z, n_voxel );

  // Process bulk of voxels 16 at a time

# define LOAD_STENCIL()     \
  f0 = &f( x,   y,   z   ); \
  fx = &f( x-1, y,   z   ); \
  fy = &f( x,   y-1, z   ); \
  fz = &f( x,   y,   z-1 )

# define NEXT_STENCIL(n)      \
  f0##n = f0++;               \
  fx##n = fx++;               \
  fy##n = fy++;               \
  fz##n = fz++;               \
  x++;                        \
  if ( x > nx )               \
  {			      \
                  x = 2, y++; \
    if ( y > ny ) y = 2, z++; \
    LOAD_STENCIL();           \
  }

  LOAD_STENCIL();

  for( ; n_voxel > 15; n_voxel -= 16 )
  {
    NEXT_STENCIL(00);
    NEXT_STENCIL(01);
    NEXT_STENCIL(02);
    NEXT_STENCIL(03);
    NEXT_STENCIL(04);
    NEXT_STENCIL(05);
    NEXT_STENCIL(06);
    NEXT_STENCIL(07);
    NEXT_STENCIL(08);
    NEXT_STENCIL(09);
    NEXT_STENCIL(10);
    NEXT_STENCIL(11);
    NEXT_STENCIL(12);
    NEXT_STENCIL(13);
    NEXT_STENCIL(14);
    NEXT_STENCIL(15);

    load_16x4_tr( &f000->cbx, &f001->cbx, &f002->cbx, &f003->cbx,
                  &f004->cbx, &f005->cbx, &f006->cbx, &f007->cbx,
                  &f008->cbx, &f009->cbx, &f010->cbx, &f011->cbx,
                  &f012->cbx, &f013->cbx, &f014->cbx, &f015->cbx,
                  f0_cbx, f0_cby, f0_cbz, f0_div_b_err );

    fx_div_b_err = v16float( fx00->div_b_err, fx01->div_b_err, fx02->div_b_err, fx03->div_b_err,
                             fx04->div_b_err, fx05->div_b_err, fx06->div_b_err, fx07->div_b_err,
                             fx08->div_b_err, fx09->div_b_err, fx10->div_b_err, fx11->div_b_err,
                             fx12->div_b_err, fx13->div_b_err, fx14->div_b_err, fx15->div_b_err );

    fy_div_b_err = v16float( fy00->div_b_err, fy01->div_b_err, fy02->div_b_err, fy03->div_b_err,
                             fy04->div_b_err, fy05->div_b_err, fy06->div_b_err, fy07->div_b_err,
                             fy08->div_b_err, fy09->div_b_err, fy10->div_b_err, fy11->div_b_err,
                             fy12->div_b_err, fy13->div_b_err, fy14->div_b_err, fy15->div_b_err );

    fz_div_b_err = v16float( fz00->div_b_err, fz01->div_b_err, fz02->div_b_err, fz03->div_b_err,
                             fz04->div_b_err, fz05->div_b_err, fz06->div_b_err, fz07->div_b_err,
                             fz08->div_b_err, fz09->div_b_err, fz10->div_b_err, fz11->div_b_err,
                             fz12->div_b_err, fz13->div_b_err, fz14->div_b_err, fz15->div_b_err );

    f0_cbx = fma( f0_div_b_err - fx_div_b_err, px, f0_cbx );
    f0_cby = fma( f0_div_b_err - fy_div_b_err, py, f0_cby );
    f0_cbz = fma( f0_div_b_err - fz_div_b_err, pz, f0_cbz );

    store_16x4_tr( f0_cbx, f0_cby, f0_cbz, f0_div_b_err,
                   &f000->cbx, &f001->cbx, &f002->cbx, &f003->cbx,
                   &f004->cbx, &f005->cbx, &f006->cbx, &f007->cbx,
                   &f008->cbx, &f009->cbx, &f010->cbx, &f011->cbx,
                   &f012->cbx, &f013->cbx, &f014->cbx, &f015->cbx );
  }

# undef NEXT_STENCIL
# undef LOAD_STENCIL
}

#else

void
clean_div_b_pipeline_v16( pipeline_args_t * args,
                          int pipeline_rank,
                          int n_pipeline )
{
  // No v16 implementation.
  ERROR( ( "No clean_div_b_pipeline_v16 implementation." ) );
}

#endif
