#include <libspe2.h>
#include <stdio.h>

extern spe_program_handle_t test_v4;

int
main( int argc,
      char **argv ) {
  spe_context_ptr_t context;
  unsigned int entry;

  printf( "Creating SPU context\n" );
  context = spe_context_create( 0, NULL );

  printf( "Loading SPU program\n" );
  spe_program_load( context, &test_v4 );

  printf( "Running SPU program\n" );
  entry = SPE_DEFAULT_ENTRY;
  spe_context_run( context, &entry, 0, NULL, NULL, NULL );

  printf( "Destroying SPU context\n" );
  spe_context_destroy( context );
  
  return 0;
}

