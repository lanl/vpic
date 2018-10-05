#ifndef _pipelines_openmp_h_
#define _pipelines_openmp_h_

#ifndef THREAD_REROUTE
#error "Do not include pipelines_openmp.h directly; use pipelines.h"
#endif

#include "../util_base.h"

#include <omp.h>

enum { MAX_PIPELINE = 272 };

// A pipeline function takes a pointer to arguments for the pipeline
// and a integer which gives the rank of the pipeline and the total
// number of pipelines dispatched.

// typedef void
// (*pipeline_func_t)( void * args,
//                     int pipeline_rank,
//                     int n_pipeline );

//----------------------------------------------------------------------------//
// Note from Evan Peters:
// _Pragma is supported in GCC 3.0+ and some Intel compilers.
//----------------------------------------------------------------------------//
// Note from Dave Nystrom:
// Need to figure out how portable this is.
//----------------------------------------------------------------------------//
// Note from Dave Nystrom:
// Review this implementation and see if it makes sense to try and package it
// and implement it more like the pthreads implementation.
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
// Generic macros that are used for all cases of vector acceleration as well
// as the standard case that does not use vector acceleration.
//----------------------------------------------------------------------------//

#define TOSTRING( a ) #a       //convert pragma directives to string

#define N_PIPELINE omp_helper.n_pipeline

#define WAIT_PIPELINES() _Pragma( TOSTRING( omp barrier ) )

#define PAD_STRUCT( sz )

//----------------------------------------------------------------------------//
// Macro defines to support v16 simd vector acceleration.  Uses thread
// dispatcher on the v16 pipeline and the caller does straggler cleanup with
// the scalar pipeline.
//----------------------------------------------------------------------------//

#if defined(V16_ACCELERATION) && defined(HAS_V16_PIPELINE)

# define EXEC_PIPELINES(name, args, str)                                   \
  _Pragma( TOSTRING( omp parallel num_threads(N_PIPELINE) shared(args) ) ) \
  {                                                                        \
    _Pragma( TOSTRING( omp for ) )                                         \
    for( int id = 0; id < N_PIPELINE; id++ )                               \
    {                                                                      \
      name##_pipeline_v16( args+id*sizeof(*args)*str, id, N_PIPELINE );    \
    }                                                                      \
  }                                                                        \
  name##_pipeline_scalar( args+str*N_PIPELINE, N_PIPELINE, N_PIPELINE );

# define PROTOTYPE_PIPELINE( name, args_t )  \
  void                                       \
  name##_pipeline_v16( args_t *args,         \
                       int pipeline_rank,    \
                       int n_pipeline );     \
                                             \
  void                                       \
  name##_pipeline_scalar( args_t *args,      \
                          int pipeline_rank, \
                          int n_pipeline )

//----------------------------------------------------------------------------//
// Macro defines to support v8 simd vector acceleration.  Uses thread
// dispatcher on the v8 pipeline and the caller does straggler cleanup with
// the scalar pipeline.
//----------------------------------------------------------------------------//

#elif defined(V8_ACCELERATION) && defined(HAS_V8_PIPELINE)

# define EXEC_PIPELINES(name, args, str)                                   \
  _Pragma( TOSTRING( omp parallel num_threads(N_PIPELINE) shared(args) ) ) \
  {                                                                        \
    _Pragma( TOSTRING( omp for ) )                                         \
    for( int id = 0; id < N_PIPELINE; id++ )                               \
    {                                                                      \
      name##_pipeline_v8( args+id*sizeof(*args)*str, id, N_PIPELINE );     \
    }                                                                      \
  }                                                                        \
  name##_pipeline_scalar( args+str*N_PIPELINE, N_PIPELINE, N_PIPELINE );

# define PROTOTYPE_PIPELINE( name, args_t )  \
  void                                       \
  name##_pipeline_v8( args_t *args,          \
                      int pipeline_rank,     \
                      int n_pipeline );      \
                                             \
  void                                       \
  name##_pipeline_scalar( args_t *args,      \
                          int pipeline_rank, \
                          int n_pipeline )

//----------------------------------------------------------------------------//
// Macro defines to support v4 simd vector acceleration.  Uses thread
// dispatcher on the v4 pipeline and the caller does straggler cleanup with
// the scalar pipeline.
//----------------------------------------------------------------------------//

#elif defined(V4_ACCELERATION) && defined(HAS_V4_PIPELINE)

# define EXEC_PIPELINES(name, args, str)                                   \
  _Pragma( TOSTRING( omp parallel num_threads(N_PIPELINE) shared(args) ) ) \
  {                                                                        \
    _Pragma( TOSTRING( omp for ) )                                         \
    for( int id = 0; id < N_PIPELINE; id++ )                               \
    {                                                                      \
      name##_pipeline_v4( args+id*sizeof(*args)*str, id, N_PIPELINE );     \
    }                                                                      \
  }                                                                        \
  name##_pipeline_scalar( args+str*N_PIPELINE, N_PIPELINE, N_PIPELINE );

# define PROTOTYPE_PIPELINE( name, args_t )  \
  void                                       \
  name##_pipeline_v4( args_t *args,          \
                      int pipeline_rank,     \
                      int n_pipeline );      \
                                             \
  void                                       \
  name##_pipeline_scalar( args_t *args,      \
                          int pipeline_rank, \
                          int n_pipeline )

//----------------------------------------------------------------------------//
// Macro defines to support the standard implementation which does not use
// explicit simd vectorization.  Uses thread dispatcher on the scalar pipeline
// and the caller does straggler cleanup with the scalar pipeline.
//----------------------------------------------------------------------------//

#else

# define EXEC_PIPELINES(name, args, str)                                   \
  _Pragma( TOSTRING( omp parallel num_threads(N_PIPELINE) shared(args) ) ) \
  {                                                                        \
    _Pragma( TOSTRING( omp for ) )                                         \
    for( int id = 0; id < N_PIPELINE; id++ )                               \
    {                                                                      \
      name##_pipeline_scalar( args+id*sizeof(*args)*str, id, N_PIPELINE ); \
    }                                                                      \
  }                                                                        \
  name##_pipeline_scalar( args+str*N_PIPELINE, N_PIPELINE, N_PIPELINE );

# define PROTOTYPE_PIPELINE( name, args_t )  \
  void                                       \
  name##_pipeline_scalar( args_t *args,      \
                          int pipeline_rank, \
                          int n_pipeline )

#endif

//----------------------------------------------------------------------------//
// A container object to mimic 'serial' and 'thread' objects.  Currently just
// useful for avoiding global var 'n_pipeline'.
//----------------------------------------------------------------------------//

typedef struct omp_container
{
  int n_pipeline;
  int dispatch_to_host;

  //const char * f_dump;
  //const char * e_dump;

  //boot gets the number of pipelines from the cmd line
  //and passes it on to the EXEC_PIPELINS macro eventually
  void
  (*boot)( int * pargc,
	   char *** pargv );
} omp_container_t;

BEGIN_C_DECLS

extern omp_container_t omp_helper;

END_C_DECLS

#endif // _pipelines_openmp_h_ 
