#define IN_checkpt
#include "checkpt_private.h"

#include <stdio.h>

checkpt_t *
checkpt_open_rdonly( const char * name ) {
  FILE * file;

  /* Check input args */

  if( !name ) ERROR(( "NULL name" ));

  /* Open the file for reading */

  file = fopen( name, "rb" );
  if( !file ) ERROR(( "Unable to open \"%s\" for checkpt read", name ));
  return (checkpt_t *)file;
}

checkpt_t *
checkpt_open_wronly( const char * name ) {
  FILE * file;

  /* Check input args */

  if( !name ) ERROR(( "NULL name" ));

  /* Open the file for writing */

  file = fopen( name, "wb" );
  if( !file ) ERROR(( "Unable to open \"%s\" for checkpt write", name ));
  return (checkpt_t *)file;
}

void
checkpt_close( checkpt_t * checkpt ) {

  /* Check input args */

  if( !checkpt ) return;

  /* Close the checpt file */

  fclose( (FILE *)checkpt );
}

void
checkpt_read( checkpt_t * checkpt,
              void * data,
              size_t sz ) {

  /* If an empty request, do nothing */

  if( !sz ) return;

  /* Check input args */

  if( !checkpt || !data ) ERROR(( "Invalid checkpt_read request" ));

  /* Read sz bytes from the checkpt file */

  if( fread( data, sz, 1, (FILE *)checkpt )!=1 )
    ERROR(( "Read error or unexpected end of checkpt file (%lu bytes)",
            (unsigned long)sz ));
}

void
checkpt_write( checkpt_t * checkpt,
               const void * data,
               size_t sz ) {

  /* If an empty request, return */

  if( !sz ) return;

  /* Check input args */

  if( !checkpt || !data ) ERROR(( "Invalid checkpt_write request" ));

  /* Write sz to the checkpt file */

  if( fwrite( data, sz, 1, (FILE *)checkpt )!=1 )
    ERROR(( "Write error or unexpected end of checkpt file (%lu bytes)",
            (unsigned long)sz ));
}
