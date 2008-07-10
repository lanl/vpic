/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Heavily revised and extended from earlier V4PIC versions
 *
 */

#include "vpic.hxx"
//#include <stdio.h>  // For fopen, fclose, fwrite, fprintf
#include <FileIO.hxx>

// FIXME: NEW FIELDS IN THE GRID READ/WRITE WAS HACKED UP TO BE BACKWARD
// COMPATIBLE WITH EXISTING EXTERNAL 3RD PARTY VISUALIZATION SOFTWARE.
// IN THE LONG RUN, THIS EXTERNAL SOFTWARE WILL NEED TO BE UPDATED.

/*****************************************************************************
 * ASCII dump IO
 *****************************************************************************/

void
vpic_simulation::dump_energies( const char *fname,
                                int append ) {
  double en_f[6], en_p;
  int rank = mp_rank(grid->mp);
  species_t *sp;
  FileIO fileIO;
  FileIOStatus status(fail);
  
  if( fname==NULL ) ERROR(("Invalid file name"));

  if( rank==0 ) {
    status = fileIO.open(fname, append ? io_write_append : io_write);
    if( status==fail ) ERROR(( "Could not open \"%s\".", fname ));
    else { 
      if( append==0 ) {
        fileIO.print( "%% Layout\n%% step ex ey ez bx by bz" );
        LIST_FOR_EACH(sp,species_list)
          fileIO.print( " \"%s\"", sp->name );
        fileIO.print( "\n" );
        fileIO.print( "%% timestep = %e\n", grid->dt );
      }
      fileIO.print( "%i", step );
    }
  }
  
  field_advance->method->energy_f( en_f, field_advance->f, field_advance->m, field_advance->g );
  if( rank==0 && status!=fail )
    fileIO.print( " %e %e %e %e %e %e",
                  en_f[0], en_f[1], en_f[2],
                  en_f[3], en_f[4], en_f[5] );
  
  LIST_FOR_EACH(sp,species_list) {
    en_p = energy_p( sp->p, sp->np, sp->q_m, interpolator, grid );
    if( rank==0 && status!=fail ) fileIO.print( " %e", en_p );
  }
  
  if( rank==0 && status!=fail ) {
    fileIO.print( "\n" );
    fileIO.close();
  }
}

// Note: dump_species/materials assume that names do not contain any \n!

void
vpic_simulation::dump_species( const char *fname ) {
  species_t *sp;
  FileIO fileIO;

  if( mp_rank(grid->mp)!=0 ) return;

  if( fname==NULL ) ERROR(( "Invalid file name" ));
  
  MESSAGE(("Dumping species to \"%s\"",fname));
  
  FileIOStatus status = fileIO.open(fname, io_write);
  
  if( status==fail ) ERROR(( "Could not open \"%s\".", fname ));
  
  LIST_FOR_EACH(sp,species_list)
    fileIO.print( "%s\n%i\n%e\n", sp->name, sp->id, sp->q_m );
  fileIO.close();
}

void
vpic_simulation::dump_materials( const char *fname ) {
  FileIO fileIO;
  material_t *m;

  if( mp_rank(grid->mp)==0 ) {
    if( fname==NULL ) ERROR(( "Invalid file name" ));
    MESSAGE(("Dumping materials to \"%s\"",fname));
    FileIOStatus status = fileIO.open(fname, io_write);
    if( status==fail ) ERROR(( "Could not open \"%s\".", fname ));
    LIST_FOR_EACH(m,material_list)
      fileIO.print( "%s\n%i\n%e %e %e\n%e %e %e\n%e %e %e\n",
               m->name, m->id,
               m->epsx,   m->epsy,   m->epsz,
               m->mux,    m->muy,    m->muz,
               m->sigmax, m->sigmay, m->sigmaz );
    fileIO.close();
  }
}

/*****************************************************************************
 * Binary dump IO
 *****************************************************************************/

/*
enum dump_types {
  grid_dump = 0,
  field_dump = 1,
  hydro_dump = 2,
  particle_dump = 3,
  restart_dump = 4
};
*/

namespace dump_type {
  const int grid_dump = 0;
  const int field_dump = 1;
  const int hydro_dump = 2;
  const int particle_dump = 3;
  const int restart_dump = 4;
} // namespace

#define WRITE_HEADER_V0(dump_type,sp_id,q_m,fileIO) BEGIN_PRIMITIVE {      \
  /* Binary compatibility information */                                 \
  WRITE( char,      CHAR_BIT,               fileIO );                      \
  WRITE( char,      sizeof(short int),      fileIO );                      \
  WRITE( char,      sizeof(int),            fileIO );                      \
  WRITE( char,      sizeof(float),          fileIO );                      \
  WRITE( char,      sizeof(double),         fileIO );                      \
  WRITE( short int, 0xcafe,                 fileIO );                      \
  WRITE( int,       0xdeadbeef,             fileIO );                      \
  WRITE( float,     1.0,                    fileIO );                      \
  WRITE( double,    1.0,                    fileIO );                      \
  /* Dump type and header format version */                              \
  WRITE( int,       0 /* Version */,        fileIO );                      \
  WRITE( int,       dump_type,              fileIO );                      \
  /* High level information */                                           \
  WRITE( int,       step,                   fileIO );                      \
  WRITE( int,       grid->nx,               fileIO );                      \
  WRITE( int,       grid->ny,               fileIO );                      \
  WRITE( int,       grid->nz,               fileIO );                      \
  WRITE( float,     grid->dt,               fileIO );                      \
  WRITE( float,     grid->dx,               fileIO );                      \
  WRITE( float,     grid->dy,               fileIO );                      \
  WRITE( float,     grid->dz,               fileIO );                      \
  WRITE( float,     grid->x0,               fileIO );                      \
  WRITE( float,     grid->y0,               fileIO );                      \
  WRITE( float,     grid->z0,               fileIO );                      \
  WRITE( float,     grid->cvac,             fileIO );                      \
  WRITE( float,     grid->eps0,             fileIO );                      \
  WRITE( float,     grid->damp,             fileIO );                      \
  WRITE( int,       mp_rank(grid->mp),      fileIO );                      \
  WRITE( int,       mp_nproc(grid->mp),     fileIO );                      \
  /* Species parameters */                                               \
  WRITE( int,       sp_id,                  fileIO );                      \
  WRITE( float,     q_m,                    fileIO );                      \
} END_PRIMITIVE

// Note dim _MUST_ be a pointer to an int

#define WRITE_ARRAY_HEADER(p,ndim,dim,fileIO) BEGIN_PRIMITIVE {    \
  WRITE( int, sizeof(p[0]), fileIO );                              \
  WRITE( int, ndim,         fileIO );                              \
  fileIO.write( dim, ndim );                        \
} END_PRIMITIVE

// The WRITE macro copies the output "value" into a temporary variable
// of the requested output "type" so that the write to the "file"
// occurs from a known binary data type. For example, if grid.dx were
// changed from a float to a double, routines which say
// WRITE(float,grid.dx,out) will still work fine without modification
// (you will lose some precision in the output file obviously). Note:
// No effort is made to convert the raw binary data from the root
// node's native format. Usually, integer data will be little endian
// 32-bit values and floating data with be little endian 32-bit IEEE
// single precision write copies. However, specialty types could be
// created so that the type cast __WRITE_tmp = (type)(value)
// automatically does the underlying conversion in C++

#define WRITE(type,value,fileIO) BEGIN_PRIMITIVE { \
  type __WRITE_tmp = (type)(value);              \
  fileIO.write( &__WRITE_tmp, 1 ); \
} END_PRIMITIVE

// Note: strlen does not include the terminating NULL
#define WRITE_STRING(string,fileIO) BEGIN_PRIMITIVE {       \
  int __WRITE_STRING_len = 0;                             \
  if( string!=NULL ) __WRITE_STRING_len = strlen(string); \
  fileIO.write( &__WRITE_STRING_len, 1 );    \
  if( __WRITE_STRING_len>0 )                              \
    fileIO.write( string, __WRITE_STRING_len );        \
} END_PRIMITIVE

#define READ(type,value,file) BEGIN_PRIMITIVE { \
  type __READ_tmp;                              \
  fread( &__READ_tmp, sizeof(type), 1, file );  \
  (value) = __READ_tmp;                         \
} END_PRIMITIVE

void
vpic_simulation::dump_grid( const char *fbase ) {
  char fname[256];
  FileIO fileIO;
  int dim[4];

  if( fbase==NULL ) ERROR(( "Invalid filename" ));

  if( mp_rank(grid->mp)==0 ) MESSAGE(("Dumping grid to \"%s\"",fbase));

  sprintf( fname, "%s.%i", fbase, mp_rank(grid->mp) );
  FileIOStatus status = fileIO.open(fname, io_write);
  if( status==fail ) ERROR(( "Could not open \"%s\".", fname ));

  WRITE_HEADER_V0( dump_type::grid_dump, invalid_species_id, 0, fileIO );

  dim[0] = 3;
  dim[1] = 3;
  dim[2] = 3;
  WRITE_ARRAY_HEADER( grid->bc, 3, dim, fileIO );
  fileIO.write( grid->bc, dim[0]*dim[1]*dim[2] );

  dim[0] = mp_nproc(grid->mp)+1;
  WRITE_ARRAY_HEADER( grid->range, 1, dim, fileIO );
  fileIO.write( grid->range, dim[0] );

  dim[0] = 6;
  dim[1] = grid->nx+2;
  dim[2] = grid->ny+2;
  dim[3] = grid->nz+2;
  WRITE_ARRAY_HEADER( grid->neighbor, 4, dim, fileIO );
  fileIO.write( grid->neighbor, dim[0]*dim[1]*dim[2]*dim[3] );

  fileIO.close();
}

void
vpic_simulation::dump_fields( const char *fbase, int ftag ) {
  char fname[256];
  FileIO fileIO;
  int dim[3];

  if( fbase==NULL ) ERROR(( "Invalid filename" ));

  if( mp_rank(grid->mp)==0 ) MESSAGE(("Dumping fields to \"%s\"",fbase));

  if( ftag ) sprintf( fname, "%s.%i.%i", fbase, step, mp_rank(grid->mp) );
  else       sprintf( fname, "%s.%i", fbase, mp_rank(grid->mp) );

  FileIOStatus status = fileIO.open(fname, io_write);
  if( status==fail ) ERROR(( "Could not open \"%s\".", fname ));

  WRITE_HEADER_V0( dump_type::field_dump, invalid_species_id, 0, fileIO );

  dim[0] = grid->nx+2;
  dim[1] = grid->ny+2;
  dim[2] = grid->nz+2;
  WRITE_ARRAY_HEADER( field_advance->f, 3, dim, fileIO );
  fileIO.write( field_advance->f, dim[0]*dim[1]*dim[2] );
  fileIO.close();
}

void
vpic_simulation::dump_hydro( const char *sp_name,
                             const char *fbase,
                             int ftag ) {
  species_t *sp;
  char fname[256];
  FileIO fileIO;
  int dim[3];

  sp = find_species_name(sp_name,species_list);
  if( sp==NULL ) ERROR(( "Invalid species \"%s\".", sp_name ));

  clear_hydro( hydro, grid );
  accumulate_hydro_p( hydro, sp->p, sp->np, sp->q_m, interpolator, grid );
  synchronize_hydro( hydro, grid );

  if( fbase==NULL ) ERROR(( "Invalid filename" ));

  if( mp_rank(grid->mp)==0 )
    MESSAGE(("Dumping \"%s\" hydro fields to \"%s\"",sp->name,fbase));

  if( ftag ) sprintf( fname, "%s.%i.%i", fbase, step, mp_rank(grid->mp) );
  else       sprintf( fname, "%s.%i", fbase, mp_rank(grid->mp) );
  FileIOStatus status = fileIO.open(fname, io_write);
  if( status==fail) ERROR(( "Could not open \"%s\".", fname ));
  
  WRITE_HEADER_V0(dump_type::hydro_dump,sp->id,sp->q_m,fileIO);

  dim[0] = grid->nx+2;
  dim[1] = grid->ny+2;
  dim[2] = grid->nz+2;
  WRITE_ARRAY_HEADER( hydro, 3, dim, fileIO );
  fileIO.write( hydro, dim[0]*dim[1]*dim[2] );
  fileIO.close();
}

void
vpic_simulation::dump_particles( const char *sp_name,
                                 const char *fbase,
                                 int ftag ) {
  species_t *sp;
  char fname[256];
  FileIO fileIO;
  int dim[1], buf_start, n_buf;
  static particle_t * ALIGNED(128) p_buf=NULL;
# define PBUF_SIZE 32768 // 1MB of particles
  
  sp = find_species_name(sp_name,species_list);
  if( sp==NULL ) ERROR(( "Invalid species name \"%s\".", sp_name ));

  if( fbase==NULL ) ERROR(( "Invalid filename" ));

  if( !p_buf ) MALLOC_ALIGNED( p_buf, PBUF_SIZE, 128 );

  if( mp_rank(grid->mp)==0 )
    MESSAGE(("Dumping \"%s\" particles to \"%s\"",sp->name,fbase));

  if( ftag ) sprintf( fname, "%s.%i.%i", fbase, step, mp_rank(grid->mp) );
  else       sprintf( fname, "%s.%i", fbase, mp_rank(grid->mp) );
  FileIOStatus status = fileIO.open(fname, io_write);
  if( status==fail ) ERROR(( "Could not open \"%s\".", fname ));
  
  WRITE_HEADER_V0( dump_type::particle_dump, sp->id, sp->q_m, fileIO );

  dim[0] = sp->np;
  WRITE_ARRAY_HEADER( p_buf, 1, dim, fileIO );

  // Copy a PBUF_SIZE hunk of the particle list into the particle
  // buffer, timecenter it and write it out. This is done this way to
  // guarantee the particle list unchanged while not requiring too
  // much memory.

  // FIXME: WITH A PIPELINED CENTER_P, PBUF NOMINALLY SHOULD BE QUITE
  // LARGE.

  for( buf_start=0; buf_start<sp->np; buf_start += PBUF_SIZE ) {
    n_buf = PBUF_SIZE;
    if( buf_start+n_buf > sp->np ) n_buf = sp->np - buf_start;
    COPY( p_buf, &sp->p[buf_start], n_buf );
    center_p( p_buf, n_buf, sp->q_m, interpolator, grid );
    fileIO.write( p_buf, n_buf );
  }

  fileIO.close();
}

// FIXME: dump_restart and restart would be much cleaner if the
// materials, fields and species modules had their own binary IO
// routines

// FIXME: put emission model and custom boundary conditions into the
// restart files

void
vpic_simulation::dump_restart( const char *fbase,
                               int ftag ) {
  char fname[256];
  FileIO fileIO;
  int dim[4];
  char *state;
  material_t *m;
  species_t *sp;

  if( fbase==NULL ) ERROR(( "Invalid filename" ));

  if( mp_rank(grid->mp)==0 ) MESSAGE(("Dumping restart to \"%s\"",fbase));

  if( ftag ) sprintf( fname, "%s.%i.%i", fbase, step, mp_rank(grid->mp) );
  else       sprintf( fname, "%s.%i", fbase, mp_rank(grid->mp) );
  FileIOStatus status = fileIO.open(fname, io_write);
  if( status==fail ) ERROR(( "Could not open \"%s\".", fname ));

  WRITE_HEADER_V0( dump_type::restart_dump, invalid_species_id, 0, fileIO );

  /************************************************
   * Restart saved directly initialized variables *
   ************************************************/

  // variables not already saved above

  WRITE( int, num_step,             fileIO );
  WRITE( int, status_interval,      fileIO );
  WRITE( int, clean_div_e_interval, fileIO );
  WRITE( int, clean_div_b_interval, fileIO );
  WRITE( int, sync_shared_interval, fileIO );
  WRITE( double, quota,             fileIO );
  WRITE( int, restart_interval,     fileIO ); 
  WRITE( int, hydro_interval,       fileIO ); 
  WRITE( int, field_interval,       fileIO ); 
  WRITE( int, particle_interval,    fileIO ); 

  /**********************************************
   * Restart saved helper initialized variables *
   **********************************************/

  // random number generator

  dim[0] = get_mt_rng_size(rng);
  MALLOC( state, dim[0] );
  get_mt_rng_state( rng, state );
  WRITE_ARRAY_HEADER( state, 1, dim, fileIO );
  fileIO.write( state, dim[0] );
  FREE( state );

  // materials

  WRITE( int, num_materials(material_list), fileIO );
  LIST_FOR_EACH(m,material_list) {
    WRITE_STRING( m->name, fileIO );
    WRITE( short int, m->id,     fileIO );
    WRITE( float,     m->epsx,   fileIO );
    WRITE( float,     m->epsy,   fileIO );
    WRITE( float,     m->epsz,   fileIO );
    WRITE( float,     m->mux,    fileIO );
    WRITE( float,     m->muy,    fileIO );
    WRITE( float,     m->muz,    fileIO );
    WRITE( float,     m->sigmax, fileIO );
    WRITE( float,     m->sigmay, fileIO );
    WRITE( float,     m->sigmaz, fileIO );
    WRITE( float,     m->zetax,  fileIO );
    WRITE( float,     m->zetay,  fileIO );
    WRITE( float,     m->zetaz,  fileIO );
  }

  // grid variables not already saved above

  WRITE( float,     grid->rdx, fileIO );
  WRITE( float,     grid->rdy, fileIO );
  WRITE( float,     grid->rdz, fileIO );
  WRITE( float,     grid->x1,  fileIO );
  WRITE( float,     grid->y1,  fileIO );
  WRITE( float,     grid->z1,  fileIO );

  dim[0] = 3;
  dim[1] = 3;
  dim[2] = 3;
  WRITE_ARRAY_HEADER( grid->bc, 3, dim, fileIO );
  fileIO.write( grid->bc, dim[0]*dim[1]*dim[2] );

  dim[0] = mp_nproc(grid->mp)+1;
  WRITE_ARRAY_HEADER( grid->range, 1, dim, fileIO );
  fileIO.write( grid->range, dim[0] );

  dim[0] = 6;
  dim[1] = grid->nx+2;
  dim[2] = grid->ny+2;
  dim[3] = grid->nz+2;
  WRITE_ARRAY_HEADER( grid->neighbor, 4, dim, fileIO );
  fileIO.write( grid->neighbor, dim[0]*dim[1]*dim[2]*dim[3] );

  dim[0] = grid->nx+2;
  dim[1] = grid->ny+2;
  dim[2] = grid->nz+2;
  WRITE_ARRAY_HEADER( grid->sfc, 3, dim, fileIO );
  fileIO.write( grid->sfc, dim[0]*dim[1]*dim[2] );

  // fields

  // FIXME: THIS IS NOT COMPLETELY ROBUST.  IT IDEALLY WOULD BE BASED
  // ON DYNAMIC LINKING (e.g. dlopen, ...) IN THE TRULY GENERAL CASE
  // (AND BUILDING VPIC WITH -rdynamic and -ldl AND POSSIBLY SOME
  // GNU DYNAMIC LINK EXTENSIONS IN THE COMMON CASE).  HOWEVER, SINCE
  // THERE ARE ENOUGH PORTABILITY CONCERNS ABOUT THIS AND THAT IT IS
  // THOGHT THERE ARE NO LIKELY VPIC USAGE SCENARIOS IN WHICH THIS
  // MIGHT FAIL ... WE KEEP OUR FINGERS CROSSED.  IN THE FUTURE,
  // THIS, THE FUNCTION DECL AND SOME OF THE BUILD SYSTEMS / MACHINES
  // SHOULD BE UPDATED TO DO THIS ROBUSTLY.  IN THE MEANTIME, THIS
  // IS THOUGHT TO BE SAFE SO LONG AS NON OF THE METHODS REFERRED
  // TO ARE IN A DYNAMICALLY LINKED LIBRARY (THIS IS TRUE NOW
  // AND EXPECTED TO BE TRUE FOR LIKELY USAGE SCENARIOS).

  WRITE( field_advance_methods_t, field_advance->method[0], fileIO );

  dim[0] = grid->nx+2;
  dim[1] = grid->ny+2;
  dim[2] = grid->nz+2;
  WRITE_ARRAY_HEADER( field_advance->f, 3, dim, fileIO );
  fileIO.write( field_advance->f, dim[0]*dim[1]*dim[2] );

  // species

  WRITE( int, num_species(species_list), fileIO );
  LIST_FOR_EACH(sp,species_list) {
    WRITE_STRING( sp->name,              fileIO );
    WRITE( int,   sp->id,                fileIO );
    WRITE( int,   sp->max_np,            fileIO );
    WRITE( int,   sp->max_nm,            fileIO );
    WRITE( float, sp->q_m,               fileIO );
    WRITE( int,   sp->sort_interval,     fileIO );
    WRITE( int,   sp->sort_out_of_place, fileIO );
    dim[0] = sp->np;
    WRITE_ARRAY_HEADER( sp->p, 1, dim, fileIO );
    if ( dim[0]>0 ) fileIO.write( sp->p, dim[0] ); 
  }

  // custom boundary condition data

  WRITE( int,        grid->nb,       fileIO );
  if( grid->nb>0 ) fileIO.write( grid->boundary, grid->nb );

  /****************************************
   * Restart-saved internal use variables *
   ****************************************/

  WRITE( double, p_time, fileIO );
  WRITE( double, g_time, fileIO );
  WRITE( double, u_time, fileIO );
  WRITE( double, f_time, fileIO );

  /****************************************
   * Restart-saved user defined variables *
   ****************************************/

  dim[0] = USER_GLOBAL_SIZE;
  WRITE_ARRAY_HEADER( user_global, 1, dim, fileIO );
  fileIO.write( user_global, dim[0] );

  fileIO.close();

  // Synchronize here in case to prevent ranks who immediately
  // terminate after a restart dump from disrupting ranks still in the
  // process of dumping

  // BJA - New modality of job termination involves putting a barrier
  // in the termination section of the input deck and removing it
  // here.  Aside from requiring fewer barriers, this is needed to
  // make turnstiles work.

  // mp_barrier( grid->mp );


}

void
vpic_simulation::restart( const char *fbase ) {
  int rank, nproc, size, ndim, dim[4], n;
  char fname[256], c, *state;
  FILE *handle = NULL;
  short int s;
  float f;
  double d;
  material_t *m, *last_m;
  species_t *sp, *last_sp;

# define ABORT(cond) if( cond ) ERROR(( #cond ))

  // Create an empty grid (creates the communicator too)
  grid  = new_grid();
  rank  = mp_rank(  grid->mp );
  nproc = mp_nproc( grid->mp );

  // Create a random number generator seeded with the rank. This will
  // be reseeded below.
  rng = new_mt_rng( rank );

  // Open the restart dump
  ABORT(fbase==NULL);
  if( rank==0 ) MESSAGE(("Restarting from \"%s\"",fbase));
  sprintf( fname, "%s.%i", fbase, rank );
  handle = fopen( fname, "rb" ); ABORT(handle==NULL);

  // Machine compatibility information
  READ(char,     c,handle);  ABORT(c!=CHAR_BIT         );
  READ(char,     c,handle);  ABORT(c!=2                );
  READ(char,     c,handle);  ABORT(c!=4                );
  READ(char,     c,handle);  ABORT(c!=4                );
  READ(char,     c,handle);  ABORT(c!=8                );
  READ(short int,s,handle);  ABORT(s!=(short int)0xcafe);
  READ(int,      n,handle);  ABORT(n!=(int)0xdeadbeef  );
  READ(float,    f,handle);  ABORT(f!=1.0              );
  READ(double,   d,handle);  ABORT(d!=1.0              );

  // Dump type and header version
  READ(int,n,handle); ABORT(n!=0           );
  READ(int,n,handle); ABORT(n!=dump_type::restart_dump);

  // High level information
  READ(int,  step,      handle); ABORT(step<0       );
  READ(int,  grid->nx,  handle); ABORT(grid->nx<1   );
  READ(int,  grid->ny,  handle); ABORT(grid->ny<1   );
  READ(int,  grid->nz,  handle); ABORT(grid->nz<1   );
  READ(float,grid->dt,  handle); ABORT(grid->dt<=0  );
  READ(float,grid->dx,  handle); ABORT(grid->dx<=0  );
  READ(float,grid->dy,  handle); ABORT(grid->dy<=0  );
  READ(float,grid->dz,  handle); ABORT(grid->dz<=0  );
  READ(float,grid->x0,  handle);
  READ(float,grid->y0,  handle);
  READ(float,grid->z0,  handle);
  READ(float,grid->cvac,handle); ABORT(grid->cvac<=0);
  READ(float,grid->eps0,handle); ABORT(grid->eps0<=0);
  READ(float,grid->damp,handle); ABORT(grid->damp<0 );
  READ(int,  n,         handle); ABORT(n!=rank      );
  READ(int,  n,         handle); ABORT(n!=nproc     );

  // Species parameters
  READ(int,  n,handle); ABORT(n!=invalid_species_id);
  READ(float,f,handle); ABORT(f!=0                 );
  
  /************************************************
   * Restart saved directly initialized variables *
   ************************************************/

  // variables not already restored above

  READ(int,    num_step,             handle); 
  READ(int,    status_interval,      handle);
  READ(int,    clean_div_e_interval, handle);
  READ(int,    clean_div_b_interval, handle);
  READ(int,    sync_shared_interval, handle);
  READ(double, quota,                handle); 
  READ(int,    restart_interval,     handle);
  READ(int,    hydro_interval,       handle);
  READ(int,    field_interval,       handle);
  READ(int,    particle_interval,    handle);

  /**********************************************
   * Restart saved helper initialized variables *
   **********************************************/

  // random number generator

  READ(int,size,  handle);           ABORT(size!=sizeof(char));
  READ(int,ndim,  handle);           ABORT(ndim!=1           );
  READ(int,dim[0],handle);           ABORT(dim[0]<0          );

  MALLOC( state, dim[0] );
  fread(state,size,dim[0],handle);
  set_mt_rng_state( rng, state, dim[0] );
  FREE( state );

  // materials ... material_list must be put together in the same
  // order it was originally in

  READ(int,n,handle); ABORT(n<1);
  for(last_m=NULL;n;n--) {
    char * buf;
    // Note: sizeof(m[0]) includes terminating '\0'
    READ(int,size,handle);                       ABORT(size<1 );
    MALLOC( buf, sizeof(m[0])+size ); m = (material_t *)buf;
    fread(m->name,1,size,handle);
    m->name[size]='\0';

    READ(short int,m->id,     handle); ABORT(m->id!=n-1);
    READ(float,    m->epsx,   handle);
    READ(float,    m->epsy,   handle);
    READ(float,    m->epsz,   handle);
    READ(float,    m->mux,    handle);
    READ(float,    m->muy,    handle);
    READ(float,    m->muz,    handle);
    READ(float,    m->sigmax, handle);
    READ(float,    m->sigmay, handle);
    READ(float,    m->sigmaz, handle);
    READ(float,    m->zetax,  handle);
    READ(float,    m->zetay,  handle);
    READ(float,    m->zetaz,  handle);

    m->next = NULL;
    if( last_m==NULL ) {
      material_list = m;
      last_m = m;
    } else {
      last_m->next = m;
      last_m = m;
    }
  }

  // grid variables not already saved above

  READ(float,grid->rdx, handle);
  READ(float,grid->rdy, handle);
  READ(float,grid->rdz, handle);
  READ(float,grid->x1,  handle);
  READ(float,grid->y1,  handle);
  READ(float,grid->z1,  handle);

  READ(int,size,  handle); ABORT(size!=sizeof(grid->bc[0]));
  READ(int,ndim,  handle); ABORT(ndim!=3          );
  READ(int,dim[0],handle); ABORT(dim[0]!=3        );
  READ(int,dim[1],handle); ABORT(dim[1]!=3        );
  READ(int,dim[2],handle); ABORT(dim[2]!=3        );
  fread( grid->bc, size, dim[0]*dim[1]*dim[2], handle );

  READ(int,size,handle);   ABORT(size!=sizeof(grid->range[0]));
  READ(int,ndim,handle);   ABORT(ndim!=1          );
  READ(int,dim[0],handle); ABORT(dim[0]!=nproc+1  );
  MALLOC_ALIGNED( grid->range, dim[0], 16 );
  fread( grid->range, size, dim[0], handle );

  READ(int,size,  handle); ABORT(size!=sizeof(grid->neighbor[0]) );
  READ(int,ndim,  handle); ABORT(ndim!=4           );
  READ(int,dim[0],handle); ABORT(dim[0]!=6         );
  READ(int,dim[1],handle); ABORT(dim[1]!=grid->nx+2);
  READ(int,dim[2],handle); ABORT(dim[2]!=grid->ny+2);
  READ(int,dim[3],handle); ABORT(dim[3]!=grid->nz+2);
  MALLOC_ALIGNED( grid->neighbor, dim[0]*dim[1]*dim[2]*dim[3], 128 );
  fread( grid->neighbor, size, dim[0]*dim[1]*dim[2]*dim[3], handle );

  READ(int,size,  handle); ABORT(size!=sizeof(grid->sfc[0]) );
  READ(int,ndim,  handle); ABORT(ndim!=3           );
  READ(int,dim[0],handle); ABORT(dim[0]!=grid->nx+2);
  READ(int,dim[1],handle); ABORT(dim[1]!=grid->ny+2);
  READ(int,dim[2],handle); ABORT(dim[2]!=grid->nz+2);
  MALLOC_ALIGNED( grid->sfc, dim[0]*dim[1]*dim[2], 128 );
  fread( grid->sfc, size, dim[0]*dim[1]*dim[2], handle );

  grid->rangel = grid->range[rank];
  grid->rangeh = grid->range[rank+1]-1;

  // fields

  // FIXME: SEE NOTE IN ABOVE ABOUT THIS!
  field_advance_methods_t fam[1];
  READ(field_advance_methods_t,fam[0],handle);
  field_advance = new_field_advance( grid, material_list, fam );
  field = field_advance->f; // FIXME: Temporary hack

  READ(int,size,  handle); ABORT(size!=sizeof(field[0]));
  READ(int,ndim,  handle); ABORT(ndim!=3              );
  READ(int,dim[0],handle); ABORT(dim[0]!=grid->nx+2   );
  READ(int,dim[1],handle); ABORT(dim[1]!=grid->ny+2   );
  READ(int,dim[2],handle); ABORT(dim[2]!=grid->nz+2   );
  fread( field_advance->f, size, dim[0]*dim[1]*dim[2], handle );

  // species ... species_list must be put together in the same order
  // it was originally in

  // FIXME: WHY IS NEW_SPECIES NOT CALLED?  (PROBABLY THE ABOVE
  // CONCERN ABOUT ORDERING ... RETOOL THIS TO CREATE THE
  // LIST AND THEN FLIP INTO CORRECT ORDERING?

  READ(int,n,handle); ABORT(n<1);
  for(last_sp=NULL;n;n--) {
    char * buf;
    // Note: sizeof(sp[0]) includes terminating '\0'
    READ(int,size,handle);                        ABORT(size<1 );
    MALLOC( buf, sizeof(sp[0])+size ); sp = (species_t *)buf;
    fread(sp->name,1,size,handle);
    sp->name[size]='\0';

    READ( int,   sp->id,                handle ); ABORT(sp->id!=n-1 );
    READ( int,   sp->max_np,            handle ); ABORT(sp->max_np<1);
    READ( int,   sp->max_nm,            handle ); ABORT(sp->max_nm<1);
    READ( float, sp->q_m,               handle );
    READ( int,   sp->sort_interval,     handle );
    READ( int,   sp->sort_out_of_place, handle );

    MALLOC_ALIGNED( sp->p, sp->max_np, 128 );
    READ(int,size,  handle); ABORT(size!=sizeof(sp->p[0])       );
    READ(int,ndim,  handle); ABORT(ndim!=1                      );
    READ(int,dim[0],handle); ABORT(dim[0]<0 || dim[0]>sp->max_np);
    if( dim[0]>0 ) fread( sp->p, size, dim[0], handle );
    sp->np = dim[0];

    MALLOC_ALIGNED( sp->pm, sp->max_nm, 128 );
    sp->nm = 0;

    sp->next = NULL;
    if( last_sp==NULL ) {
      species_list = sp;
      last_sp = sp;
    } else {
      last_sp->next = sp;
      last_sp = sp;
    }
    
    sp->partition=NULL;  /* FIXME: SHOULD REALLOCATE THIS */
  }

  // custom boundary condition data
  
  READ( int, grid->nb, handle ); ABORT( grid->nb<0 );
  if( grid->nb==0 ) grid->boundary = NULL;
  else {
    MALLOC( grid->boundary, grid->nb );
    fread( grid->boundary,  sizeof(grid->boundary[0]), grid->nb, handle );
  }

  /****************************************
   * Restart-saved internal use variables *
   ****************************************/

  READ( double, p_time, handle ); ABORT(p_time<0);
  READ( double, g_time, handle ); ABORT(g_time<0);
  READ( double, u_time, handle ); ABORT(u_time<0);
  READ( double, f_time, handle ); ABORT(f_time<0);

  /****************************************
   * Restart-saved user defined variables *
   ****************************************/

  READ(int,size, handle); ABORT(size!=sizeof(user_global[0]));
  READ(int,ndim, handle); ABORT(ndim!=1                     );
  READ(int,dim[0],handle); ABORT(dim[0]!=USER_GLOBAL_SIZE   );
  fread( user_global, size, dim[0], handle );

  fclose(handle);
  handle = NULL;

  /**********************************
   * Recreate other data structures *
   **********************************/

  interpolator = new_interpolator(grid);
  load_interpolator( interpolator, field_advance->f, field_advance->g );
  accumulator = new_accumulators(grid);
  hydro = new_hydro(grid);

  mp_size_recv_buffer(BOUNDARY(-1, 0, 0),(grid->ny+1)*(grid->nz+1)*sizeof(hydro[0]),grid->mp);
  mp_size_recv_buffer(BOUNDARY( 1, 0, 0),(grid->ny+1)*(grid->nz+1)*sizeof(hydro[0]),grid->mp);
  mp_size_recv_buffer(BOUNDARY( 0,-1, 0),(grid->nz+1)*(grid->nx+1)*sizeof(hydro[0]),grid->mp);
  mp_size_recv_buffer(BOUNDARY( 0, 1, 0),(grid->nz+1)*(grid->nx+1)*sizeof(hydro[0]),grid->mp);
  mp_size_recv_buffer(BOUNDARY( 0, 0,-1),(grid->nx+1)*(grid->ny+1)*sizeof(hydro[0]),grid->mp);
  mp_size_recv_buffer(BOUNDARY( 0, 0, 1),(grid->nx+1)*(grid->ny+1)*sizeof(hydro[0]),grid->mp);
  mp_size_send_buffer(BOUNDARY(-1, 0, 0),(grid->ny+1)*(grid->nz+1)*sizeof(hydro[0]),grid->mp);
  mp_size_send_buffer(BOUNDARY( 1, 0, 0),(grid->ny+1)*(grid->nz+1)*sizeof(hydro[0]),grid->mp);
  mp_size_send_buffer(BOUNDARY( 0,-1, 0),(grid->nz+1)*(grid->nx+1)*sizeof(hydro[0]),grid->mp);
  mp_size_send_buffer(BOUNDARY( 0, 1, 0),(grid->nz+1)*(grid->nx+1)*sizeof(hydro[0]),grid->mp);
  mp_size_send_buffer(BOUNDARY( 0, 0,-1),(grid->nx+1)*(grid->ny+1)*sizeof(hydro[0]),grid->mp);
  mp_size_send_buffer(BOUNDARY( 0, 0, 1),(grid->nx+1)*(grid->ny+1)*sizeof(hydro[0]),grid->mp);
  mp_barrier(grid->mp);
}

// Add capability to modify certain fields "on the fly" so that one
// can, e.g., extend a run, change a quota, or modify a dump interval
// without having to rerun from the start.
//
// File is specified in arg 3 slot in command line inputs.  File is in
// ASCII format with each field in the form: field val [newline].
//
// Allowable values of field variables are: num_steps, quota,
// restart_interval, hydro_interval, field_interval, particle_interval
//
// [x]_interval sets interval value for dump type [x].  Set interval
// to zero to turn off dump type.

# define SETIVAR( V, A, S )                                                     \
  {                                                                             \
    V=(A);                                                                      \
    if ( mp_rank(grid->mp)==0 ) MESSAGE(("* Modifying %s to value %d", S, A));  \
  }
# define SETDVAR( V, A, S )                                                     \
  {                                                                             \
    V=(A);                                                                      \
    if ( mp_rank(grid->mp)==0 ) MESSAGE(("* Modifying %s to value %le", S, A));  \
  }
# define ITEST( V, N, A ) if ( sscanf( line, N " %d", &iarg )==1 ) SETIVAR( V, A, N ); 
# define DTEST( V, N, A ) if ( sscanf( line, N " %le", &darg )==1 ) SETDVAR( V, A, N ); 

void
vpic_simulation::modify_runparams( const char *fname ) {
  FILE *handle=NULL; 
  char line[128];
  int iarg=0; 
  double darg=0; 

  // Open the modfile
  if( mp_rank(grid->mp)==0 ) MESSAGE(("Modifying run parameters from file \"%s\"", fname));
  handle = fopen( fname, "r" ); 
  if ( handle==NULL ) {
    ERROR(( "Modfile read failed %s(%i)[%i]: %s",
            __FILE__, __LINE__, mp_rank(grid->mp), "handle==NULL" ));
    if( handle!=NULL ) fclose(handle); handle=NULL;
    mp_abort(__LINE__,grid->mp);
  }
  // Parse modfile
  while ( fgets( line, 127, handle ) ) {
    DTEST( quota,             "quota",             darg ); 
    ITEST( num_step,          "num_step",          iarg ); 
    ITEST( restart_interval,  "restart_interval",  (iarg<0 ? 0 : iarg) );
    ITEST( hydro_interval,    "hydro_interval",    (iarg<0 ? 0 : iarg) );
    ITEST( field_interval,    "field_interval",    (iarg<0 ? 0 : iarg) );
    ITEST( particle_interval, "particle_interval", (iarg<0 ? 0 : iarg) );
  }
}

#undef SETIVAR
#undef SETDVAR
#undef ITEST
#undef DTEST
