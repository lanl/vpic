#ifndef _pipelines_exec_pth_h_
#define _pipelines_exec_pth_h_

#ifndef THREAD_REROUTE
#error "Do not include pipelines_exec_pth.h directly; use pipelines_exec.h."
#endif

//----------------------------------------------------------------------------//
// Generic macros that are used for all cases of vector acceleration as well
// as the standard case that does not use vector acceleration.
//----------------------------------------------------------------------------//

# define WAIT_PIPELINES() thread.wait()

//----------------------------------------------------------------------------//
// Macro defines to support v16 simd vector acceleration.  Uses thread
// dispatcher on the v16 pipeline and the caller does straggler cleanup with
// the scalar pipeline.
//----------------------------------------------------------------------------//

#if defined(V16_ACCELERATION) && defined(HAS_V16_PIPELINE)

# define EXEC_PIPELINES(name,args,str)                           \
  thread.dispatch( (pipeline_func_t)name##_pipeline_v16,         \
                   args, sizeof(*args), str );                   \
  name##_pipeline_scalar( args+str*N_PIPELINE, N_PIPELINE, N_PIPELINE )

//----------------------------------------------------------------------------//
// Macro defines to support v8 simd vector acceleration.  Uses thread
// dispatcher on the v8 pipeline and the caller does straggler cleanup with
// the scalar pipeline.
//----------------------------------------------------------------------------//

#elif defined(V8_ACCELERATION) && defined(HAS_V8_PIPELINE)

# define EXEC_PIPELINES(name,args,str)                           \
  thread.dispatch( (pipeline_func_t)name##_pipeline_v8,          \
                   args, sizeof(*args), str );                   \
  name##_pipeline_scalar( args+str*N_PIPELINE, N_PIPELINE, N_PIPELINE )

//----------------------------------------------------------------------------//
// Macro defines to support v4 simd vector acceleration.  Uses thread
// dispatcher on the v4 pipeline and the caller does straggler cleanup with
// the scalar pipeline.
//----------------------------------------------------------------------------//

#elif defined(V4_ACCELERATION) && defined(HAS_V4_PIPELINE)

# define EXEC_PIPELINES(name,args,str)                           \
  thread.dispatch( (pipeline_func_t)name##_pipeline_v4,          \
                   args, sizeof(*args), str );                   \
  name##_pipeline_scalar( args+str*N_PIPELINE, N_PIPELINE, N_PIPELINE )

//----------------------------------------------------------------------------//
// Macro defines to support the standard implementation which does not use
// explicit simd vectorization.  Uses thread dispatcher on the scalar pipeline
// and the caller does straggler cleanup with the scalar pipeline.
//----------------------------------------------------------------------------//

#else

# define EXEC_PIPELINES(name,args,str)                           \
  thread.dispatch( (pipeline_func_t)name##_pipeline_scalar,      \
                   args, sizeof(*args), str );                   \
  name##_pipeline_scalar( args+str*N_PIPELINE, N_PIPELINE, N_PIPELINE )

#endif

#endif // _pipelines_exec_pth_h_ 
