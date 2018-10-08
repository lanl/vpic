#ifndef _pipelines_exec_omp_h_
#define _pipelines_exec_omp_h_

#ifndef THREAD_REROUTE
#error "Do not include pipelines_exec_omp.h directly; use pipelines_exec.h"
#endif

//----------------------------------------------------------------------------//
// Generic macros that are used for all cases of vector acceleration as well
// as the standard case that does not use vector acceleration.
//----------------------------------------------------------------------------//

#define TOSTRING( a ) #a       //convert pragma directives to string

#define WAIT_PIPELINES() _Pragma( TOSTRING( omp barrier ) )

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

//# define PROTOTYPE_PIPELINE( name, args_t )	\
//  void					\
//  name##_pipeline_v16( args_t *args,		\
//                       int pipeline_rank,	\
//                       int n_pipeline );	\
//						\
//  void					\
//  name##_pipeline_scalar( args_t *args,	\
//                          int pipeline_rank,	\
//                          int n_pipeline )

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

//# define PROTOTYPE_PIPELINE( name, args_t )	\
//  void					\
//  name##_pipeline_v8( args_t *args,		\
//                      int pipeline_rank,	\
//                      int n_pipeline );	\
//						\
//  void					\
//  name##_pipeline_scalar( args_t *args,	\
//                          int pipeline_rank,	\
//                          int n_pipeline )

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

//# define PROTOTYPE_PIPELINE( name, args_t )	\
//  void					\
//  name##_pipeline_v4( args_t *args,		\
//                      int pipeline_rank,	\
//                      int n_pipeline );	\
//						\
//  void					\
//  name##_pipeline_scalar( args_t *args,	\
//                          int pipeline_rank,	\
//                          int n_pipeline )

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

//# define PROTOTYPE_PIPELINE( name, args_t )	\
//  void					\
//  name##_pipeline_scalar( args_t *args,	\
//                          int pipeline_rank,	\
//                          int n_pipeline )

#endif

#endif // _pipelines_exec_omp_h_ 
