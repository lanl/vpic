#include <pipelines.h> /* For common.h, datatypes and prototypes */

static int busy = 0;

static void
serial_boot( int n_pipeline_requested ) {
  if( serial._n_pipeline!=0 ) ERROR(( "Pipelines already booted!" ));
  if( n_pipeline_requested < 1 || n_pipeline_requested > MAX_PIPELINE )
    ERROR(( "Invalid number of pipelines requested" ));

  serial._n_pipeline = n_pipeline_requested;
  busy = 0;
}

static void
serial_dispatch( pipeline_func_t pipeline,
                 void * args,
                 int size_args ) {
  int p;

  if( serial._n_pipeline==0 ) ERROR(( "Boot pipelines first!" ));
  if( busy ) ERROR(( "Pipelines are busy!" ));
  if( pipeline==NULL ) ERROR(( "Bad pipeline" ));

  busy = 1;
  for( p=0; p<serial._n_pipeline; p++ )
    pipeline( ((char *)args) + p*size_args, p );

}

static void
serial_wait( void ) {
  if( serial._n_pipeline==0 ) ERROR(( "Boot pipelines first!" ));
  if( !busy ) ERROR(( "Pipelines are not busy!" ));

  busy = 0;
}

static void
serial_halt( void ) {
  if( serial._n_pipeline==0 ) ERROR(( "Boot pipelines first!" ));
  if( busy ) ERROR(( "Pipelines are busy!" ));

  serial._n_pipeline = 0;
  busy = 0;
}

pipeline_dispatcher_t serial = {
  0,               /* n_pipeline */
  serial_boot,     /* boot */
  serial_halt,     /* halt */
  serial_dispatch, /* dispatch */
  serial_wait      /* wait */
};

