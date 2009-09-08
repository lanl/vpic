/*#define VERBOSE_CHECKPOINTING*/
#ifndef NO_REVERSE_SYMBOL_TABLE_LOOKUP_SUPPORT
#define _GNU_SOURCE
#include <dlfcn.h>
#endif

#define IN_checkpt
#include "checkpt_private.h"

/* Boolean flag indicating whether or not checkpoint is booted. */

static int booted = 0;

/* If NULL, indicates that we are not in the middle of a checkpt or a
   restore.  Otherwise, it gives the handle of the stream used for I/O
   operations. */

static checkpt_t * checkpt = NULL;
static checkpt_t * restore = NULL;

/* The registry is a list of objects that need to checkpointed (in the
   order they should be checkpointed).  The registry gives each object
   a unique identifier that is invariant across a checkpt/restore and
   includes the details how to checkpt and restore each object. */

typedef struct registry {
  void * obj;
  checkpt_func_t checkpt_func;
  restore_func_t restore_func;
  reanimate_func_t reanimate_func;
  size_t id;
  struct registry * next;
} registry_t;

static registry_t * registry = NULL;

/* Counter used to dole out unique ids.  Zero ids are used to indicate
   an error condition. */

static size_t next_id = 1;

#ifdef VERBOSE_CHECKPOINTING

static void
dump_node( registry_t * node ) {
  Dl_info dli[1];
  if( node->checkpt_func ) dladdr( (void *)(size_t)node->checkpt_func, dli );
  else                     dli->dli_sname = NULL, dli->dli_fname = NULL;
  log_printf( "Object %lu:\n"
              "\taddr=%lx\n"
              "\tcheckpt_func=%lx (%s from %s)\n",
              (unsigned long)node->id,
              (unsigned long)node->obj,
              (unsigned long)node->checkpt_func,
              dli->dli_sname==NULL ? "(null)" : dli->dli_sname,
              dli->dli_fname==NULL ? "(null)" : dli->dli_fname );
  if( node->restore_func ) dladdr( (void *)(size_t)node->restore_func, dli );
  else                     dli->dli_sname = NULL, dli->dli_fname = NULL;
  log_printf( "\trestore_func=%lx (%s from %s)\n",
              (unsigned long)node->restore_func,
              dli->dli_sname==NULL ? "(null)" : dli->dli_sname,
              dli->dli_fname==NULL ? "(null)" : dli->dli_fname );
  if( node->reanimate_func ) dladdr( (void *)(size_t)node->reanimate_func,dli);
  else                       dli->dli_sname = NULL, dli->dli_fname = NULL;
  log_printf( "\treanimate_func=%lx (%s from %s)\n",
              (unsigned long)node->reanimate_func,
              dli->dli_sname==NULL ? "(null)" : dli->dli_sname,
              dli->dli_fname==NULL ? "(null)" : dli->dli_fname );
}

static void
dump_registry( void ) {
  registry_t * node;
  for( node=registry; node; node=node->next ) dump_node( node );
}

#else

static void
dump_node( registry_t * node ) {
}

static void
dump_registry( void ) {
}

#endif

void
boot_checkpt( int * pargc,
              char *** pargv ) {

  /* Check input args */

  if( booted ) ERROR(( "checkpt service already booted." ));

  /* Initialize the checkpt service state */

  registry = NULL;
  checkpt  = NULL;
  restore  = NULL;
  next_id  = 1;

  /* Mark the service as booted */

  booted = 1;
}

void
halt_checkpt( void ) {

  /* Check input args */

  if( !booted  ) ERROR(( "checkpt service not booted." ));
  if( registry ) {
    dump_registry();
    ERROR(( "halt called with some objects still registered" ));
  }
  if( checkpt  ) ERROR(( "currently writing a checkpt" ));
  if( restore  ) ERROR(( "currently reading a checkpt" ));

  /* Mark the service as halted */

  booted = 0;
}

size_t
object_id( const void * obj ) {
  registry_t * node;

  /* Check input args */

  if( !booted ) ERROR(( "checkpt service not booted" ));

  /* Search the registry for the object */

  for( node=registry; node; node=node->next )
    if( node->obj==obj ) break;

  /* If found, return the id, otherwise, return 0 */

  return node ? node->id : 0;
}

void *
object_ptr( size_t id ) {
  registry_t * node;

  /* Check input args */

  if( !booted ) ERROR(( "checkpt service not booted" ));

  /* Search the registry for the id */

  for( node=registry; node; node=node->next )
    if( node->id==id ) break;

  /* If found, return the id, otherwise, return NULL */

  return node ? node->obj : NULL;
}

void
register_object( void * obj,
                 checkpt_func_t checkpt_func,
                 restore_func_t restore_func,
                 reanimate_func_t reanimate_func ) {
  registry_t * node, * prev;

  /* Check input args */

  if( !booted ) ERROR(( "checkpt service not booted" ));
  if( checkpt ) ERROR(( "currently writing a checkpt" ));
  if( restore ) ERROR(( "currently reading a checkpt" ));

  /* Check that obj is valid and that obj isn't already registered.
     At the same time, find the last entry in the registry. */

  if( !obj ) ERROR(( "NULL obj" ));
  for( prev=NULL, node=registry; node; prev=node, node=node->next )
    if( node->obj==obj ) break;
  if( node ) ERROR(( "object already registered" ));
  
  /* Create the registry entry for this object */

  MALLOC( node, 1 );
  node->id             = next_id++;
  node->obj            = obj;
  node->checkpt_func   = checkpt_func;
  node->restore_func   = restore_func;
  node->reanimate_func = reanimate_func;
  node->next           = NULL;

  /* And append it to the end of the registry */

  if( prev ) prev->next = node;
  if( !registry ) registry = node;
}

void
unregister_object( void * obj ) {
  registry_t * node, * prev;

  /* Check input args */

  if( !booted ) ERROR(( "checkpt service not booted" ));
  if( checkpt ) ERROR(( "currently writing a checkpt" ));
  if( restore ) ERROR(( "currently reading a checkpt" ));

  /* Find the entry for this object in the register and the previous
     entry.  If the object is not registered, return an error */

  for( prev=NULL, node=registry; node; prev=node, node=node->next )
    if( node->obj==obj ) break;
  if( !node ) ERROR(( "object was not registered" ));

  /* Remove the node from the registry */

  if( node==registry ) registry = node->next;
  if( prev ) prev->next = node->next;
  node->next = NULL;

  /* And delete it */

  FREE( node );
}

void
checkpt_objects( const char * name ) {
  registry_t * node;

  /* Check input args */

  if( !booted ) ERROR(( "checkpt not booted" ));
  if( checkpt ) ERROR(( "currently writing a checkpt" ));
  if( restore ) ERROR(( "currently reading a checkpt" ));

  /* Open the checkpt serialization stream */

  checkpt = checkpt_open_wronly( name );
  CHECKPT_VAL( size_t, next_id );

  /* Checkpoint the objects */

  for( node=registry; node; node=node->next ) {
    dump_node( node );
    CHECKPT_VAL( size_t, 0x600DF00D );
    checkpt_raw( node, sizeof(*node) );
    checkpt_sym( (void *)(size_t)node->checkpt_func   );
    checkpt_sym( (void *)(size_t)node->restore_func   );
    checkpt_sym( (void *)(size_t)node->reanimate_func );
    if( node->checkpt_func ) node->checkpt_func( node->obj );
  }

  /* Mark that there are no more objects in the stream, close the
     serialization stream and indicate that we are no longer writing a
     checkpt */
  
  CHECKPT_VAL( size_t, 0xBADF00D );
  checkpt_close( checkpt );
  checkpt = NULL;
}

void
restore_objects( const char * name ) {
  registry_t * node, * prev;
  size_t prefix;

  /* Check input args */

  if( !booted ) ERROR(( "checkpt not booted" ));
  if( checkpt ) ERROR(( "currently writing a checkpt" ));
  if( restore ) ERROR(( "currently reading a checkpt" ));

  /* Delete all objects in the in favor of the checkpointed objects */

  node = registry;
  while( node ) {
    prev = node;
    node = node->next;
    FREE( prev );
  }
  registry = NULL;
  next_id = 0;

  /* Open the checkpt deserialization stream */

  restore = checkpt_open_rdonly( name );
  RESTORE_VAL( size_t, next_id );

  /* Restore the objects */

  prev = NULL;
  for(;;) {
    RESTORE_VAL( size_t, prefix );
    if( prefix== 0xBADF00D ) break;
    if( prefix!=0x600DF00D )
      ERROR(( "Malformed checkpt (expected an object header)" ));
    MALLOC( node, 1 );
    restore_raw( node, sizeof(*node) );
    node->checkpt_func   = (checkpt_func_t)  (size_t)restore_sym();
    node->restore_func   = (restore_func_t)  (size_t)restore_sym();
    node->reanimate_func = (reanimate_func_t)(size_t)restore_sym();
    node->next = NULL;
    if( !registry ) registry = node;
    if( prev ) prev->next = node;
    prev = node;
    dump_node( node );
    if( node->restore_func ) node->obj = node->restore_func();
  }

  /* Close the checkpt deserialization stream and indicate that we are
     no longer reading a checkpt */

  checkpt_close( restore );
  restore = NULL;
}

void
reanimate_objects( void ) {
  registry_t * node;

  /* Check input args */

  if( !booted ) ERROR(( "checkpt not booted" ));
  if( checkpt ) ERROR(( "currently writing a checkpt" ));
  if( restore ) ERROR(( "currently reading a checkpt" ));

  /* Call each objects reanimate function */

  for( node=registry; node; node=node->next ) {
    dump_node( node );
    if( node->reanimate_func ) node->reanimate_func( node->obj );
  }
}

/* Primitive checkpt helpers */

void
checkpt_raw( const void * data,
             size_t n_byte ) {

  /* Check input args */

  if( !checkpt ) ERROR(( "not writing a checkpt" ));
  if( !data && n_byte ) ERROR(( "NULL data" ));

  /* Write data to the serialization stream */

  if( n_byte ) checkpt_write( checkpt, data, n_byte );
}

void
restore_raw( void * data,
             size_t n_byte ) {

  /* Check input args */

  if( !restore ) ERROR(( "not reading a checkpt" ));
  if( !data && n_byte ) ERROR(( "NULL data" ));

  /* Read data from the deserialization stream */

  if( n_byte ) checkpt_read( restore, data, n_byte );
}

/* Composiite checkpt helpers */

void
checkpt_data( const void * _data,
              size_t sz_ele,
              size_t str_ele,
              size_t n_ele,
              size_t max_ele,
              size_t align ) {
  const char * data = (const char *)_data;
  size_t n;

  /* Check input args */

  if( !data && n_ele*sz_ele>0 ) ERROR(( "NULL data" ));
  if( n_ele>max_ele || sz_ele>str_ele ) ERROR(( "bad data layout" ));

  /* Write the data header */

  CHECKPT_VAL( size_t, 0xDA7A );
  CHECKPT_VAL( size_t, sz_ele ); CHECKPT_VAL( size_t, str_ele );
  CHECKPT_VAL( size_t, n_ele  ); CHECKPT_VAL( size_t, max_ele );
  CHECKPT_VAL( size_t, align  );

  /* Write out the individual elements */

  for( n=0; n<n_ele; n++ ) checkpt_raw( data+n*str_ele, sz_ele );
}

void *
restore_data( void ) {
  char * data;
  size_t n, sz_ele, str_ele, n_ele, max_ele, align;

  /* Read the data header */

  RESTORE_VAL( size_t, n );
  if( n!=0xDA7A ) ERROR(( "malformed checkpt (expected a data header)" ));
  RESTORE_VAL( size_t, sz_ele ); RESTORE_VAL( size_t, str_ele );
  RESTORE_VAL( size_t, n_ele  ); RESTORE_VAL( size_t, max_ele );
  RESTORE_VAL( size_t, align  );
  if( n_ele>max_ele || sz_ele>str_ele )
    ERROR(( "malformed checkpt (invalid data layout)" ));

  /* Allocate the data according to the header */

  if( align==0 ) MALLOC(         data, max_ele*str_ele        );
  else           MALLOC_ALIGNED( data, max_ele*str_ele, align );

  /* And read in the checkpointed elements */

  for( n=0; n<n_ele; n++ ) restore_raw( data+n*str_ele, sz_ele );
  return data;
}

void
checkpt_str( const char * str ) {
  
  /* If str is NULL, write a NULL string header to the stream.
     Otherwise, write a non-NULL string header and the string itself
     (not including the terminating '\0'). */

  if( !str ) CHECKPT_VAL( size_t, 0x2c11577 );
  else {
    size_t len = strlen(str);
    CHECKPT_VAL( size_t, 0x577 );
    CHECKPT_VAL( size_t, len   );
    checkpt_raw( str, len );
  }
}

char *
restore_str( void ) {
  char * str;
  size_t len;

  /* Read the string header and decide what kind of string we
     received */

  RESTORE_VAL( size_t, len );
  if( len==0x2c11577 ) return NULL; /* NULL string received */
  if( len!=0x577 )
    ERROR(( "malformed checkpt (expected a non-NULL string header)" ));

  /* We received a non-NULL string.  Get its length, allocate it,
     read it (appending the terminating '\0') and check that it
     restored okay. */

  RESTORE_VAL( size_t, len );
  MALLOC( str, len+1 );
  restore_raw( str, len );
  str[len] = '\0';
  if( strlen(str)!=len )
    ERROR(( "malformed checkpt (malformed string data)" ));

  return str;
}

void
checkpt_fptr( const void * ptr ) {

  /* If ptr is NULL, write a NULL pointer header to the stream.
     Otherwise, write a non-NULL pointer header and the unique
     identifier for the ptr itself. */

  if( !ptr ) CHECKPT_VAL( size_t, 0x2c11977 );
  else {
    size_t id = object_id( ptr );
    if( !id )
      ERROR(( "cannot checkpoint a pointer to an unregistered object" ));
    CHECKPT_VAL( size_t, 0x977 );
    CHECKPT_VAL( size_t, id    );
  }
}

void *
restore_fptr( void ) {
  size_t id;

  /* Read the pointer header and decide what kind of pointer we
     received */

  RESTORE_VAL( size_t, id );
  if( id==0x2c11977 ) return NULL; /* NULL pointer received */
  if( id!=0x977 )
    ERROR(( "malformed checkpt (expected a non-NULL pointer header)" ));

  /* We received a non-NULL pointer.  Get its unique identifier and
     determine what it points to locally. */

  RESTORE_VAL( size_t, id );
  if( !id ) ERROR(( "malformed checkpt (expected an object identifier)" ));

  return (void *)id;
}

void *
reanimate_fptr( void * ptr ) {
  if( ptr==NULL ) return NULL;
  ptr = object_ptr( (size_t)ptr );
  if( !ptr )
    ERROR(( "Unable to reanimate a pointer.  Mostly likely, either the "
            "checkpt is malformed or a restore function is trying to restore "
            "a pointer to an object that hasn't been restored yet (the second "
            "case shouldn't happen if, as is typically the case, the "
            "checkpointed object dependency network is a DAG ... use "
            "restore_fptr in the restore and reanimate_fptr functions in the "
            "reanimate functions to restore pointers in a non-DAG dependency "
            "network if not doing so currently)." ));
  return ptr;
}

void
checkpt_ptr( const void * ptr ) {
  checkpt_fptr( ptr );
}

void *
restore_ptr( void ) {
  return reanimate_fptr( restore_fptr() );
}

#ifdef NO_REVERSE_SYMBOL_TABLE_LOOKUP_SUPPORT

void
checkpt_sym( const void * saddr ) {

  /* If symbol address is NULL, checkpoint a NULL symbol header.
     Otherwise, write a address symbol header and the symbol
     address. */
  
  if( saddr==NULL ) CHECKPT_VAL( size_t, 0x2C11513B );
  else {
    static int first_time = 1;
    if( first_time ) {
      if( !world_rank ) 
        WARNING(( "Checkpointing was compiled without reverse symbol table "
                  "lookup support and has been asked to checkpoint at least "
                  "one symbol.  Checkpointing will likely only work correctly "
                  "if checkpointed symbols are statically linked in the "
                  "application and the exact same application binary is used "
                  "for both the checkpoint and restore process." ));
      first_time = 0;
    }

    CHECKPT_VAL( size_t, 0xADD7513B  );
    CHECKPT_VAL( const void *, saddr );    
  }
}

void *
restore_sym( void ) {
  static int first_time = 1;
  size_t type;
  void * saddr;

  /* Read the symbol header and determine the symbol type. */

  RESTORE_VAL( size_t, type );
  if( type==0x2C11513B ) return NULL; /* NULL symbol header */
  if( type!=0xADD7513B ) 
    ERROR(( "Malformed checkpt (expected symbol header)" ));

  /* Read the symbol address */

  if( first_time ) {
    if( !world_rank ) 
      WARNING(( "Checkpointing was compiled without reverse symbol table "
                "lookup support and has been asked to restore at least one "
                "symbol.  Checkpointing will likely only work correctly if "
                "checkpointed symbols are statically linked in the "
                "application and the exact same application binary is used "
                "for both the checkpoint and restore process." ));
    first_time = 0;
  }

  RESTORE_VAL( void *, saddr );
  return saddr;
}

#else

/* Note that this function is only used to look for non-NULL symbols */

static void *
find_saddr( const char * sname,
            const char * fname ) {
  void * saddr, * fbase;
  const char * err;

  /* If we aren't given a symbol name, return not found */

  if( sname==NULL ) return NULL;
  
  /* Note that the default library search path (RTLD_DEFAULT) includes
     the application itself.  The application should be linked with
     with the global symbol table exported (e.g. "-rdynamic" under
     gcc) for this to work. */

  dlerror();
  saddr = dlsym( RTLD_DEFAULT, sname );
  err = dlerror();
  if( !err ) return saddr;

  /* If we didn't find it in the library search path, try opening the
     library given by the explicit fname parameter (if provided).

     Note: Currently, this code path is disabled as fname returned by
     dladdr seems not be to useful (it typically is the application's
     name, which dlopen cannot actually open) and there are open
     issues of how to call dlopen (RTLD_LAZY, RTLD_NOW, RTLD_GLOBAL)
     and what and how to call dlclose. */

  if( 0 && fname ) {
    dlerror();
    fbase = dlopen( fname, RTLD_LAZY );
    err = dlerror();
    if( err ) WARNING(( "dlopen(%s,RTLD_LAZY): %s", fname, err ));
    else {
      dlerror();
      saddr = dlsym( fbase, sname );
      err = dlerror();
      if( !err ) return saddr;
    }
  }

  return NULL;
}

void
checkpt_sym( const void * saddr ) {
  Dl_info dli[1];
  const char * err;

  /* If the symbol is NULL, checkpoint a NULL symbol header */
  
  if( saddr==NULL ) {
    CHECKPT_VAL( size_t, 0x2C11513B );
    return;
  }

  /* Do a reverse symbol lookup */

  dlerror();
  dladdr( saddr, dli );
  err = dlerror();
  if( err ) ERROR(( "dladdr: %s", err ));

  if( saddr==find_saddr( dli->dli_sname, dli->dli_fname ) ) {

    /* If we can resolve the symbol using the strings given to us by
       the reverse lookup, checkpoint the symbol strings (presumably,
       these strings are invariant across a checkpt/restore and thus
       we should be able to do this resolution again on the
       restore). */

    CHECKPT_VAL( size_t, 0x577513B );
    checkpt_str( dli->dli_sname );
    checkpt_str( dli->dli_fname );

  } else {

    /* Otherwise, just checkpoint the raw symbol address (if we can't
       resolve the address with the strings now, there is virtually no
       chance we will be able to resolve it with the strings).  The
       raw symbol address is only invariant under very limited
       circumstances (e.g. symbol is in a statically linked section of
       the code and the exact same binary used for writing the
       checkpoint is used to read the checkpoint).  So, we give a
       stern warning to the user that bad things are possible when
       restoring this symbol. */

    WARNING(( "Unable to find a safely writable symbol that corresponds to "
              "address %lx (the closest match was \"%s\" from \"%s\").  "
              "Writing out the raw address instead and keeping my fingers "
              "crossed.\n", (unsigned long)saddr,
              dli->dli_sname==NULL ? "(null)" : dli->dli_sname,
              dli->dli_fname==NULL ? "(null)" : dli->dli_fname ));

    CHECKPT_VAL( size_t, 0xADD7513B  );
    CHECKPT_VAL( const void *, saddr );
    
  }
}

void *
restore_sym( void ) {
  size_t type;
  void * saddr;
  char * sname, * fname;
  Dl_info dli[1];
  const char * err;

  /* Read the symbol header */

  RESTORE_VAL( size_t, type );

  switch( type ) {

  default:         /* Invalid header */
    ERROR(( "Malformed checkpt (expected symbol header)" ));
    saddr = NULL;
    break;

  case 0x2C11513B: /* NULL symbol */
    saddr = NULL;
    break;

  case  0x577513B: /* string symbol */

    /* Restore the symbol strings */

    sname = restore_str();
    fname = restore_str();

    /* Find the symbol's address */

    saddr = find_saddr( sname, fname );   
    if( !saddr )
      ERROR(( "Unable to resolve the symbol \"%s\" (from \"%s\") in the "
              "current symbol table.  When this symbol was written, it was "
              "resolvable in old symbol table.  Thus, make sure you have "
              "not changed the application to eliminate this symbol or give "
              "it a slightly different name, that your paths to various "
              "libraries are set correctly, that the same library versions "
              "are being used and so forth to allow this symbol to be "
              "resolved correctly.", sname, fname ));

    /* Free the strings */

    FREE( sname );
    FREE( fname );
    break;

  case 0xADD7513B: /* address symbol */

    /* Read the old symbol's address. */

    RESTORE_VAL( void *, saddr );

    /* Check to see if the symbol address makes any sense currently. */
    
    dlerror();
    dladdr( saddr, dli );
    err = dlerror();
    if( err ) ERROR(( "dladdr: %s", err ));
    if( saddr==find_saddr( dli->dli_sname, dli->dli_fname ) )
      WARNING(( "Read a symbol that was saved unsafely as the symbol's "
                "address (%lx) in the old symbol table.  This address "
                "resolves to symbol \"%s\" (from \"%s\") in the current "
                "symbol table.  Keeping my fingers crossed that this "
                "corresponds to the saved symbol.", (unsigned long)saddr,
                dli->dli_sname, dli->dli_fname ));
    else
      WARNING(( "Read a symbol that was saved unsafely as the symbol's "
                "address (%lx) in the old symbol table.  This address was not "
                "found in the current symbol table; the closest match is "
                "symbol \"%s\" (from \"%s\") at address %lx.  Keeping my "
                "fingers crossed that the old address still means "
                "something.\n", (unsigned long)saddr, dli->dli_sname,
                dli->dli_fname, (unsigned long)dli->dli_saddr ));
    break;

  }

  return saddr;
}

#endif

void
_cxx_illegal_ptr_copy( void * lv_ref,
                       const void * rv ) {
  if( lv_ref ) (*((const void **)lv_ref)) = rv;
}

