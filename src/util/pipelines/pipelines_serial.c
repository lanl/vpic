#include <pipelines.h> // For util_base.h, datatypes and prototypes

static int Busy = 0;

static void
serial_boot( int n_pipeline,
             int dispatch_to_host ) {
  if( serial.n_pipeline!=0 ) ERROR(( "Serial dispatcher already booted!" ));
  if( n_pipeline<1 || n_pipeline>MAX_PIPELINE )
    ERROR(( "Invalid number of pipelines requested" ));
  serial.n_pipeline = n_pipeline;
  Busy = 0;
}

static void
serial_dispatch( pipeline_func_t func,
                 void * args,
                 int size_args ) {
  int id;

  if( serial.n_pipeline==0 ) ERROR(( "Boot serial dispatcher first!" ));
  if( Busy ) ERROR(( "Pipelines are busy!" ));
  Busy = 1;
  for( id=0; id<serial.n_pipeline; id++ )
    if( func ) func( ((char *)args) + id*size_args, id, serial.n_pipeline );
}

static void
serial_wait( void ) {
  if( serial.n_pipeline==0 ) ERROR(( "Boot serial dispatcher first!" ));
  if( !Busy ) ERROR(( "Pipelines are not busy!" ));
  Busy = 0;
}

static void
serial_halt( void ) {
  if( serial.n_pipeline==0 ) ERROR(( "Boot serial dispatcher first!" ));
  if( Busy ) ERROR(( "Pipelines are busy!" ));
  serial.n_pipeline = 0;
}

pipeline_dispatcher_t serial = {
  0,               // n_pipeline
  serial_boot,     // boot
  serial_halt,     // halt
  serial_dispatch, // dispatch
  serial_wait      // wait
};

