#include <pipeline.h> /* For common.h, datatypes and prototypes */
 
void
_dispatch_pipelines( pipeline_func_t pipeline,
                     void * args,
                     int size_args,
                     pipeline_request_t * request ) {
  int p;
  if( pipeline==NULL ) ERROR(( "Bad pipeline" ));
  
  for( p=0; p<n_pipeline; p++ ) {
    if( request!=NULL ) request->complete[p] = 0;
    pipeline( ((char *)args) + p*size_args, p );
    if( request!=NULL ) request->complete[p] = 1;
  }
}

void
wait_for_pipelines( pipeline_request_t * request ) {
  int p;
  if( request!=NULL ) for( p=0; p<n_pipeline; p+=request->complete[p] );
}

