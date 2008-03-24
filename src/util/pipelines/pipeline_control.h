#ifndef pipeline_control_h
#define pipeline_control_h

#define FOR_SPU ( defined(CELL_SPU_BUILD) || \
	( defined(CELL_PPU_BUILD)    &&   \
	defined(USE_CELL_SPUS)    &&   \
	defined(HAS_SPU_PIPELINE) ) )

////////////////////////////////////////////////////////////////////////////////
// Pipeline Functions
////////////////////////////////////////////////////////////////////////////////

// use SPU dispatcher on the SPU pipeline
#if FOR_SPU

	#include <spe_events.h>

	#if defined (CELL_PPU_BUILD)

		#define INIT_PIPELINES(name,args,sz_args) \
			spu.dispatch(SPU_PIPELINE(name##_pipeline_spu), args, sz_args)

		#define EXEC_PIPELINES(name,args,sz_args) \
            name##_pipeline_signal(); \
			name##_pipeline( args, spu.n_pipeline, spu.n_pipeline )

		#define WAIT_PIPELINES() \
			spu.sync(UPDATE_COMPLETE)

		#define FINALIZE_PIPELINES() \
			spu.signal(END_EVENT_LOOP); \
			spu.wait()
		
		#define N_PIPELINE spu.n_pipeline

	#else

		// SPUs cannot dispatch pipelines

	#endif // CELL_PPU_BUILD

// Use thread dispatcher on the v4 pipeline
#elif defined(V4_ACCELERATION) && defined(HAS_V4_PIPELINE)

	#define INIT_PIPELINES(name,args,sz_args)

	#define EXEC_PIPELINES(name,args,sz_args) \
		thread.dispatch( (pipeline_func_t)name##_pipeline_v4, args, sz_args ); \
		name##_pipeline( args, thread.n_pipeline, thread.n_pipeline )

	#define WAIT_PIPELINES() \
		thread.wait()

	#define FINALIZE_PIPELINES()
	
	#define N_PIPELINE thread.n_pipeline

// Use thread dispatcher on the scalar pipeline
#else 

	#define INIT_PIPELINES(name,args,sz_args)

	#define EXEC_PIPELINES(name,args,sz_args) \
		thread.dispatch( (pipeline_func_t)name##_pipeline_v4, args, sz_args ); \
		name##_pipeline( args, thread.n_pipeline, thread.n_pipeline )

	#define WAIT_PIPELINES() \
		thread.wait()

	#define FINALIZE_PIPELINES()
	
	#define N_PIPELINE thread.n_pipeline

#endif // PIPELINE type

#endif // pipeline_control_h
