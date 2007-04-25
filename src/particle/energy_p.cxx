#include <particle.h>

#ifndef V4_ACCELERATION
#define ENERGY_P_PIPELINE (pipeline_func_t)energy_p_pipeline
#else
#define ENERGY_P_PIPELINE (pipeline_func_t)energy_p_pipeline_v4
#endif

typedef struct energy_p_pipeline_args {
  MEM_PTR( const particle_t,     128 ) p0;   // Particle array
  MEM_PTR( const interpolator_t, 128 ) f0;   // Interpolator array
  MEM_PTR( double,               128 ) en;   // Return values
  float                                q_m;  // Charge to mass ratio
  float                                cvac; // Speed of light
  float                                dt;   // Timestep
  int                                  np;   // Number of particles
# ifdef USE_CELL_SPUS // Align to 16-bytes
  char _pad[PAD(3*SIZEOF_MEM_PTR+3*sizeof(float)+sizeof(int),16)];
# endif
} energy_p_pipeline_args_t;

static void
energy_p_pipeline( energy_p_pipeline_args_t * args,
                   int pipeline_rank,
                   int n_pipeline ) {
  const interpolator_t * ALIGNED(128) f0  = args->f0;

  const particle_t     * ALIGNED(32)  p;
  const interpolator_t * ALIGNED(16)  f;

  const float qdt_2mc = 0.5*args->q_m*args->dt/args->cvac;
  const float one     = 1.;

  float dx, dy, dz;
  float v0, v1, v2;

  double en = 0;

  int first, n;

  // Determine which particles this pipeline processes

  DISTRIBUTE( args->np, 4, pipeline_rank, n_pipeline, first, n );
  p = args->p0 + first;

  // Process particles quads for this pipeline

  for(;n;n--,p++) {
    dx  = p->dx;
    dy  = p->dy;
    dz  = p->dz;
    f   = f0 + p->i;
    v0  = p->ux + qdt_2mc*(    ( f->ex    + dy*f->dexdy    ) +
                            dz*( f->dexdz + dy*f->d2exdydz ) );
    v1  = p->uy + qdt_2mc*(    ( f->ey    + dz*f->deydz    ) +
                            dx*( f->deydx + dz*f->d2eydzdx ) );
    v2  = p->uz + qdt_2mc*(    ( f->ez    + dx*f->dezdx    ) +
                            dy*( f->dezdy + dx*f->d2ezdxdy ) );
    v0  = v0*v0 + v1*v1 + v2*v2;
    v0 /= (float)sqrt(one+v0)+one;
    en += (double)v0*(double)p->q;
  }

  args->en[pipeline_rank] = en;
}

#ifdef V4_ACCELERATION

using namespace v4;

static void
energy_p_pipeline_v4( energy_p_pipeline_args_t * args,
                      int pipeline_rank,
                      int n_pipeline ) {
  const interpolator_t * ALIGNED(128) f0 = args->f0;

  const particle_t     * ALIGNED(128) p;
  const float          * ALIGNED(16)  vp0;
  const float          * ALIGNED(16)  vp1;
  const float          * ALIGNED(16)  vp2;
  const float          * ALIGNED(16)  vp3;

  const v4float qdt_2mc(0.5*args->q_m*args->dt/args->cvac);
  const v4float one(1.);

  v4float dx, dy, dz;
  v4float ex, ey, ez;
  v4float v0, v1, v2, q;
  v4int ii;

  double en0 = 0, en1 = 0, en2 = 0, en3 = 0;

  int first, nq;

  // Determine which particle quads this pipeline processes

  DISTRIBUTE( args->np, 4, pipeline_rank, n_pipeline, first, nq );
  p = args->p0 + first;
  nq >>= 2;

  // Process the particle quads for this pipeline

  for( ; nq; nq--, p+=4 ) {
    load_4x4_tr(&p[0].dx,&p[1].dx,&p[2].dx,&p[3].dx,dx,dy,dz,ii);

    // Interpolate fields
    vp0 = (float *)(f0 + ii(0));
    vp1 = (float *)(f0 + ii(1));
    vp2 = (float *)(f0 + ii(2));
    vp3 = (float *)(f0 + ii(3));
    load_4x4_tr(vp0,  vp1,  vp2,  vp3,  ex,v0,v1,v2); ex = fma( fma( dy, v2, v1 ), dz, fma( dy, v0, ex ) );
    load_4x4_tr(vp0+4,vp1+4,vp2+4,vp3+4,ey,v0,v1,v2); ey = fma( fma( dz, v2, v1 ), dx, fma( dz, v0, ey ) );
    load_4x4_tr(vp0+8,vp1+8,vp2+8,vp3+8,ez,v0,v1,v2); ez = fma( fma( dx, v2, v1 ), dy, fma( dx, v0, ez ) );

    // Update momentum to half step
    // (note Boris rotation does not change energy so it is unnecessary)
    load_4x4_tr(&p[0].ux,&p[1].ux,&p[2].ux,&p[3].ux,v0,v1,v2,q);
    v0 = fma( ex, qdt_2mc, v0 );
    v1 = fma( ey, qdt_2mc, v1 );
    v2 = fma( ez, qdt_2mc, v2 );

    // Accumulate energy
    v0  = fma( v0,v0, fma( v1,v1, v2*v2 ) );
    v0 /= sqrt(one+v0)+one; 
    en0 += (double)v0(0)*(double)q(0);
    en1 += (double)v0(1)*(double)q(1);
    en2 += (double)v0(2)*(double)q(2);
    en3 += (double)v0(3)*(double)q(3);
  }

  args->en[pipeline_rank] = en0 + en1 + en2 + en3;
}

#endif

double
energy_p( const particle_t     * ALIGNED(128) p0,
          const int                           np,
          const float                         q_m,
          const interpolator_t * ALIGNED(128) f0,
          const grid_t         *              g ) {
  DECLARE_ALIGNED_ARRAY( energy_p_pipeline_args_t, 128, args, 1 );
  DECLARE_ALIGNED_ARRAY( double, 128, en, MAX_PIPELINE+1 );
  double local, global;
  int rank;

  // FIXME: p0 NULL checking
  if( np<0     ) ERROR(("Bad number of particles"));
  if( f0==NULL ) ERROR(("Bad interpolator"));
  if( g==NULL  ) ERROR(("Bad grid"));

  // Have the pipelines do the bulk of particles in quads and have the
  // host do the final incomplete quad.

  args->p0   = p0;
  args->f0   = f0;
  args->en   = en;
  args->q_m  = q_m;
  args->cvac = g->cvac;
  args->dt   = g->dt;
  args->np   = np;

# ifdef USE_CELL_SPUS

  spu.dispatch( SPU_PIPELINE( energy_p_pipeline_spu ), args, 0 );
  energy_p_pipeline( args, spu.n_pipeline, spu.n_pipeline );
  spu.wait();

# else

  PSTYLE.dispatch( ENERGY_P_PIPELINE, args, 0 );
  energy_p_pipeline( args, PSTYLE.n_pipeline, PSTYLE.n_pipeline );
  PSTYLE.wait();

# endif

  local = 0;
  for( rank=0; rank<=PSTYLE.n_pipeline; rank++ )
    local += en[rank];
  mp_allsum_d( &local, &global, 1, g->mp );
  return (double)g->cvac*(double)g->cvac*global/(double)q_m;
}
