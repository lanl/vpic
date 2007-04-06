#include <pipeline.h> /* For common.h, datatypes and prototypes */

static int busy = 0;
 
void
_dispatch_pipelines( pipeline_func_t pipeline,
                     void * args,
                     int size_args ) {
  int p;
  if( pipeline==NULL ) ERROR(( "Bad pipeline" ));

  if( busy ) ERROR(( "Pipelines are busy!" ));
  busy = 1;
  
  for( p=0; p<n_pipeline; p++ )
    pipeline( ((char *)args) + p*size_args, p );
}

void
wait_for_pipelines( void ) {
  if( !busy ) ERROR(( "Pipelines are not busy!" ));
  busy = 0;
}

