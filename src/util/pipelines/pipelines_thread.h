#ifndef _pipelines_thread_h_
#define _pipelines_thread_h_

void checkpt_thread( const pipeline_dispatcher_t * _thread );
pipeline_dispatcher_t* restore_thread( void );

#endif // _pipelines_thread_h_
