/*
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Original version
 *
 */

#include "sf_interface.h"

static int
ha_n_pipeline( void )
{
#if defined(VPIC_USE_PTHREADS)                         // Pthreads case.
  int                          n = serial.n_pipeline;
  if ( n < thread.n_pipeline ) n = thread.n_pipeline;

#elif defined(VPIC_USE_OPENMP)                         // OpenMP case.
  int                          n = omp_helper.n_pipeline;

#else                                                  // Error case.
  #error "VPIC_USE_OPENMP or VPIC_USE_PTHREADS must be specified"

#endif

  return n; // max( {serial,thread}.n_pipeline )
}

// Though the checkpt / restore functions are not part of the public
// API, they must not be declared as static.

void
checkpt_hydro_array( const hydro_array_t * ha )
{
  CHECKPT( ha, 1 );

  CHECKPT_ALIGNED( ha->h,
                   (size_t) ( ha->n_pipeline + 1 ) * (size_t) ha->stride,
                   128 );

  CHECKPT_PTR( ha->g );
}

hydro_array_t *
restore_hydro_array( void )
{
  hydro_array_t * ha;

  RESTORE( ha );

  RESTORE_ALIGNED( ha->h );

  RESTORE_PTR( ha->g );

  if ( ha->n_pipeline != ha_n_pipeline() )
  {
    ERROR( ( "Number of pipelines restored is not the same as the number of "
             "pipelines checkpointed.  Did you change the number of threads "
             "per process between checkpt and restore?" ) );
  }

  return ha;
}

hydro_array_t *
new_hydro_array( grid_t * g )
{
  hydro_array_t * ha;

  if ( ! g )
  {
    ERROR( ( "NULL grid." ) );
  }

  MALLOC( ha, 1 );

  ha->n_pipeline = ha_n_pipeline();
  ha->stride     = POW2_CEIL( g->nv, 2 );
  ha->g          = g;

  MALLOC_ALIGNED( ha->h,
                  (size_t) ( ha->n_pipeline + 1 ) * (size_t) ha->stride,
                  128 );

  CLEAR( ha->h,
         (size_t) ( ha->n_pipeline + 1 ) * (size_t) ha->stride );

  REGISTER_OBJECT( ha,
                   checkpt_hydro_array,
                   restore_hydro_array,
                   NULL );

  return ha;
}

void
delete_hydro_array( hydro_array_t * ha )
{
  if ( ! ha )
  {
    return;
  }

  UNREGISTER_OBJECT( ha );

  FREE_ALIGNED( ha->h );

  FREE( ha );
}

#define hydro( x, y, z ) h0[ VOXEL( x, y, z, nx, ny, nz ) ]

// Generic looping.

#define XYZ_LOOP( xl, xh, yl, yh, zl, zh ) \
  for( z = zl; z <= zh; z++ )              \
    for( y = yl; y <= yh; y++ )            \
      for( x = xl; x <= xh; x++ )

// x_NODE_LOOP => Loop over all non-ghost nodes at plane x.

#define x_NODE_LOOP( x ) XYZ_LOOP( x, x, 1, ny+1, 1, nz+1 )

#define y_NODE_LOOP( y ) XYZ_LOOP( 1, nx+1, y, y, 1, nz+1 )

#define z_NODE_LOOP( z ) XYZ_LOOP( 1, nx+1, 1, ny+1, z, z )

void
synchronize_hydro_array( hydro_array_t * ha )
{
  int size, face, bc, x, y, z, nx, ny, nz;

  float *p, lw, rw;

  hydro_t *h0, *h;

  grid_t *g;

  if ( ! ha )
  {
    ERROR( ( "NULL hydro array." ) );
  }

  // First reduce the pipelines.

  reduce_hydro_array( ha );

  // Now begin to synchronize.

  h0 = ha->h;
  g  = ha->g;

  nx = g->nx;
  ny = g->ny;
  nz = g->nz;

  // Note: synchronize_hydro assumes that hydro has not been adjusted
  // at the local domain boundary. Because hydro fields are purely
  // diagnostic, correct the hydro along local boundaries to account
  // for accumulations over partial cell volumes.

  #define ADJUST_HYDRO( i, j, k, X, Y, Z )      \
  do                                            \
  {                                             \
    bc = g->bc[ BOUNDARY( i, j, k ) ];          \
    if ( bc < 0 || bc >= world_size )           \
    {                                           \
      face = ( i + j + k ) < 0 ? 1 : n##X + 1;  \
      X##_NODE_LOOP( face )                     \
      {                                         \
        h       = &hydro( x, y, z );            \
        h->jx  *= 2;                            \
        h->jy  *= 2;                            \
        h->jz  *= 2;                            \
        h->rho *= 2;                            \
        h->px  *= 2;                            \
        h->py  *= 2;                            \
        h->pz  *= 2;                            \
        h->ke  *= 2;                            \
        h->txx *= 2;                            \
        h->tyy *= 2;                            \
        h->tzz *= 2;                            \
        h->tyz *= 2;                            \
        h->tzx *= 2;                            \
        h->txy *= 2;                            \
      }                                         \
    }                                           \
  } while( 0 )

  ADJUST_HYDRO( (-1),   0,    0,  x, y, z );
  ADJUST_HYDRO(   0,  (-1),   0,  y, z, x );
  ADJUST_HYDRO(   0,    0,  (-1), z, x, y );
  ADJUST_HYDRO(   1,    0,    0,  x, y, z );
  ADJUST_HYDRO(   0,    1,    0,  y, z, x );
  ADJUST_HYDRO(   0,    0,    1,  z, x, y );

  #undef ADJUST_HYDRO

  #define BEGIN_RECV( i, j, k, X, Y, Z ) \
  begin_recv_port( i, j, k, ( 1 + 14 * ( n##Y + 1 ) * ( n##Z + 1 ) ) * sizeof(float), g )

  #define BEGIN_SEND( i, j, k, X, Y, Z )                     \
  BEGIN_PRIMITIVE                                            \
  {                                                          \
    size = ( 1 + 14 * (n##Y+1) * (n##Z+1) ) * sizeof(float); \
    p = (float *) size_send_port( i, j, k, size, g );        \
    if ( p )                                                 \
    {                                                        \
      ( *(p++) ) = g->d##X;                                  \
      face       = ( i + j + k ) < 0 ? 1 : n##X+1;           \
      X##_NODE_LOOP( face )                                  \
      {                                                      \
        h          = &hydro( x, y, z );                      \
        ( *(p++) ) = h->jx;                                  \
        ( *(p++) ) = h->jy;                                  \
        ( *(p++) ) = h->jz;                                  \
        ( *(p++) ) = h->rho;                                 \
        ( *(p++) ) = h->px;                                  \
        ( *(p++) ) = h->py;                                  \
        ( *(p++) ) = h->pz;                                  \
        ( *(p++) ) = h->ke;                                  \
        ( *(p++) ) = h->txx;                                 \
        ( *(p++) ) = h->tyy;                                 \
        ( *(p++) ) = h->tzz;                                 \
        ( *(p++) ) = h->tyz;                                 \
        ( *(p++) ) = h->tzx;                                 \
        ( *(p++) ) = h->txy;                                 \
      }                                                      \
      begin_send_port( i, j, k, size, g );                   \
    }                                                        \
  } END_PRIMITIVE

  #define END_RECV( i, j, k, X, Y, Z )                                 \
  BEGIN_PRIMITIVE                                                      \
  {                                                                    \
    p = (float *) end_recv_port( i, j, k, g );                         \
    if ( p )                                                           \
    {                                                                  \
      rw    = ( *(p++) );                     /* Remote g->d##X */     \
      lw    = rw + g->d##X;                                            \
      rw   /= lw;                                                      \
      lw    = g->d##X / lw;                                            \
      lw   += lw;                                                      \
      rw   += rw;                                                      \
      face  = ( i + j + k ) < 0 ? n##X+1 : 1; /* Twice weighted sum */ \
      X##_NODE_LOOP( face )                                            \
      {                                                                \
        h      = &hydro( x, y, z );                                    \
        h->jx  = lw * h->jx  + rw * ( *(p++) );                        \
        h->jy  = lw * h->jy  + rw * ( *(p++) );                        \
        h->jz  = lw * h->jz  + rw * ( *(p++) );                        \
        h->rho = lw * h->rho + rw * ( *(p++) );                        \
        h->px  = lw * h->px  + rw * ( *(p++) );                        \
        h->py  = lw * h->py  + rw * ( *(p++) );                        \
        h->pz  = lw * h->pz  + rw * ( *(p++) );                        \
        h->ke  = lw * h->ke  + rw * ( *(p++) );                        \
        h->txx = lw * h->txx + rw * ( *(p++) );                        \
        h->tyy = lw * h->tyy + rw * ( *(p++) );                        \
        h->tzz = lw * h->tzz + rw * ( *(p++) );                        \
        h->tyz = lw * h->tyz + rw * ( *(p++) );                        \
        h->tzx = lw * h->tzx + rw * ( *(p++) );                        \
        h->txy = lw * h->txy + rw * ( *(p++) );                        \
      }                                                                \
    }                                                                  \
  } END_PRIMITIVE

  #define END_SEND( i, j, k, X, Y, Z ) end_send_port( i, j, k, g )

  // Exchange x-faces.

  BEGIN_SEND( (-1), 0, 0, x, y, z );
  BEGIN_SEND(   1,  0, 0, x, y, z );
  BEGIN_RECV( (-1), 0, 0, x, y, z );
  BEGIN_RECV(   1,  0, 0, x, y, z );

  END_RECV( (-1), 0, 0, x, y, z );
  END_RECV(   1,  0, 0, x, y, z );
  END_SEND( (-1), 0, 0, x, y, z );
  END_SEND(   1,  0, 0, x, y, z );

  // Exchange y-faces.

  BEGIN_SEND( 0, (-1), 0, y, z, x );
  BEGIN_SEND( 0,   1,  0, y, z, x );
  BEGIN_RECV( 0, (-1), 0, y, z, x );
  BEGIN_RECV( 0,   1,  0, y, z, x );

  END_RECV( 0, (-1), 0, y, z, x );
  END_RECV( 0,   1,  0, y, z, x );
  END_SEND( 0, (-1), 0, y, z, x );
  END_SEND( 0,   1,  0, y, z, x );

  // Exchange z-faces.

  BEGIN_SEND( 0, 0, (-1), z, x, y );
  BEGIN_SEND( 0, 0,   1,  z, x, y );
  BEGIN_RECV( 0, 0, (-1), z, x, y );
  BEGIN_RECV( 0, 0,   1,  z, x, y );

  END_RECV( 0, 0, (-1), z, x, y );
  END_RECV( 0, 0,   1,  z, x, y );
  END_SEND( 0, 0, (-1), z, x, y );
  END_SEND( 0, 0,   1,  z, x, y );

  #undef BEGIN_RECV
  #undef BEGIN_SEND
  #undef END_RECV
  #undef END_SEND
}
