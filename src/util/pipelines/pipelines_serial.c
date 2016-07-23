#include "pipelines.h" // For util_base.h, datatypes and prototypes

static int Busy = 0;

/*****************************************************************************/

#include "../checkpt/checkpt.h"

void
checkpt_serial( const pipeline_dispatcher_t * _serial ) {
  CHECKPT_VAL( int, serial.n_pipeline );
}

pipeline_dispatcher_t *
restore_serial( void ) {
  int n_pipeline;
  RESTORE_VAL( int, n_pipeline );
  if( serial.n_pipeline!=n_pipeline )
    ERROR(( "--serial.n_pipeline changed between checkpt (%i) and "
            "restore (%i)", serial.n_pipeline, n_pipeline ));
  return &serial;
}

/*****************************************************************************/

static void
serial_boot( int * pargc,
             char *** pargv ) {
  if( serial.n_pipeline!=0 ) ERROR(( "Serial dispatcher already booted!" ));
  serial.n_pipeline = strip_cmdline_int(pargc,pargv,"--serial.n_pipeline",1);
  if( serial.n_pipeline<1 || serial.n_pipeline>MAX_PIPELINE )
    ERROR(( "Invalid number of pipelines requested (%i)", serial.n_pipeline ));
  REGISTER_OBJECT( &serial, checkpt_serial, restore_serial, NULL );
  Busy = 0;
}

static void
serial_dispatch( pipeline_func_t func,
                 void * args,
                 int sz,
                 int str ) {
  int id;

  if( serial.n_pipeline==0 ) ERROR(( "Boot serial dispatcher first!" ));
  if( Busy ) ERROR(( "Pipelines are busy!" ));
  Busy = 1;
  for( id=0; id<serial.n_pipeline; id++ )
    if( func ) func( ((char *)args) + id*sz*str, id, serial.n_pipeline );
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
  UNREGISTER_OBJECT( &serial );
  serial.n_pipeline = 0;
}

pipeline_dispatcher_t serial = {
  0,               // n_pipeline
  serial_boot,     // boot
  serial_halt,     // halt
  serial_dispatch, // dispatch
  serial_wait      // wait
};
