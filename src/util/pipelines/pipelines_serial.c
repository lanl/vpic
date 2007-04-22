#include <pipelines.h> // For util_base.h, datatypes and prototypes

static int busy = 0;
static int in_pipeline = 0;

static void
serial_boot( int n_pipeline_requested ) {
  if( serial.n_pipeline!=0 ) ERROR(( "Pipelines already booted!" ));
  if( n_pipeline_requested < 1 || n_pipeline_requested > MAX_PIPELINE )
    ERROR(( "Invalid number of pipelines requested" ));

  serial.n_pipeline = n_pipeline_requested;
  busy = 0;
}

static void
serial_dispatch( pipeline_func_t pipeline,
                 void * args,
                 int size_args ) {
  int p;

  if( serial.n_pipeline==0 ) ERROR(( "Boot pipelines first!" ));
  if( in_pipeline ) ERROR(( "Only the host can call this!" ));
  if( busy ) ERROR(( "Pipelines are busy!" ));
  if( pipeline==NULL ) ERROR(( "Bad pipeline" ));

  busy = 1;
  in_pipeline = 1;
  for( p=0; p<serial.n_pipeline; p++ )
    pipeline( ((char *)args) + p*size_args, p, serial.n_pipeline );
  in_pipeline = 0;
}

static void
serial_wait( void ) {
  if( serial.n_pipeline==0 ) ERROR(( "Boot pipelines first!" ));
  if( in_pipeline ) ERROR(( "Only the host can call this!" ));
  if( !busy ) ERROR(( "Pipelines are not busy!" ));
  busy = 0;
}

static void
serial_halt( void ) {
  if( serial.n_pipeline==0 ) ERROR(( "Boot pipelines first!" ));
  if( in_pipeline ) ERROR(( "Only the host can call this!" ));
  if( busy ) ERROR(( "Pipelines are busy!" ));

  serial.n_pipeline = 0;
  busy = 0;
  in_pipeline = 0;
}

pipeline_dispatcher_t serial = {
  0,               // n_pipeline
  serial_boot,     // boot
  serial_halt,     // halt
  serial_dispatch, // dispatch
  serial_wait      // wait
};

