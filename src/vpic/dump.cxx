/*
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Heavily revised and extended from earlier V4PIC versions
 *
 * snell - revised to add strided dumps, time history dumps, others  20080404
 */
 
#include "vpic.hxx"
#include "dumphead.h"
#include "dumpvars.h"
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
 
#define FILETYPE FILE*
// #define FILETYPE FileIO
FILETYPE fdgrd[NFILMX];
FILETYPE fdfld[NFILMX];
FILETYPE fdhyd[NFILMX];
FILETYPE fdpar[NFILMX];
FILETYPE fdheader;
FILETYPE fdcat;
 
FILETYPE fdhis[NFILHISMX];
char fbasehis[NFILHISMX][256];
 
namespace dump_type {
  const int grid_dump = 0;
  const int field_dump = 1;
  const int hydro_dump = 2;
  const int particle_dump = 3;
  const int restart_dump = 4;
  const int history_dump = 5;
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
 
FILE *fildesc;
 
#define F_WRITE_HEADER_V0(dump_type,sp_id,q_m,fildesc) BEGIN_PRIMITIVE { \
  /* Binary compatibility information */                     \
  F_WRITE( char,      CHAR_BIT,               fildesc );     \
  F_WRITE( char,      sizeof(short int),      fildesc );     \
  F_WRITE( char,      sizeof(int),            fildesc );     \
  F_WRITE( char,      sizeof(float),          fildesc );     \
  F_WRITE( char,      sizeof(double),         fildesc );     \
  F_WRITE( short int, 0xcafe,                 fildesc );     \
  F_WRITE( int,       0xdeadbeef,             fildesc );     \
  F_WRITE( float,     1.0,                    fildesc );     \
  F_WRITE( double,    1.0,                    fildesc );     \
  /* Dump type and header format version */                  \
  F_WRITE( int,       0 /* Version */,        fildesc );     \
  F_WRITE( int,       dump_type,              fildesc );     \
  /* High level information */                               \
  F_WRITE( int,       step,                   fildesc );     \
  F_WRITE( int,       imxstr-2,               fildesc );     \
  F_WRITE( int,       jmxstr-2,               fildesc );     \
  F_WRITE( int,       kmxstr-2,               fildesc );     \
  F_WRITE( float,     grid->dt,               fildesc );     \
  F_WRITE( float,     dxstr,                  fildesc );     \
  F_WRITE( float,     dystr,                  fildesc );     \
  F_WRITE( float,     dzstr,                  fildesc );     \
  F_WRITE( float,     grid->x0,               fildesc );     \
  F_WRITE( float,     grid->y0,               fildesc );     \
  F_WRITE( float,     grid->z0,               fildesc );     \
  F_WRITE( float,     grid->cvac,             fildesc );     \
  F_WRITE( float,     grid->eps0,             fildesc );     \
  F_WRITE( float,     grid->damp,             fildesc );     \
  F_WRITE( int,       mp_rank(grid->mp),      fildesc );     \
  F_WRITE( int,       mp_nproc(grid->mp),     fildesc );     \
  /* Species parameters */                                   \
  F_WRITE( int,       sp_id,                  fildesc );     \
  F_WRITE( float,     q_m,                    fildesc );     \
} END_PRIMITIVE
 
#define F_WRITE_HEADER_PAR(dump_type,sp_id,q_m,fildesc) BEGIN_PRIMITIVE { \
  /* Binary compatibility information */                     \
  F_WRITE( char,      CHAR_BIT,               fildesc );     \
  F_WRITE( char,      sizeof(short int),      fildesc );     \
  F_WRITE( char,      sizeof(int),            fildesc );     \
  F_WRITE( char,      sizeof(float),          fildesc );     \
  F_WRITE( char,      sizeof(double),         fildesc );     \
  F_WRITE( short int, 0xcafe,                 fildesc );     \
  F_WRITE( int,       0xdeadbeef,             fildesc );     \
  F_WRITE( float,     1.0,                    fildesc );     \
  F_WRITE( double,    1.0,                    fildesc );     \
  /* Dump type and header format version */                  \
  F_WRITE( int,       0 /* Version */,        fildesc );     \
  F_WRITE( int,       dump_type,              fildesc );     \
  /* High level information */                               \
  F_WRITE( int,       step,                   fildesc );     \
  F_WRITE( int,       grid->nx,               fildesc );     \
  F_WRITE( int,       grid->ny,               fildesc );     \
  F_WRITE( int,       grid->nz,               fildesc );     \
  F_WRITE( float,     grid->dt,               fildesc );     \
  F_WRITE( float,     grid->dx,               fildesc );     \
  F_WRITE( float,     grid->dy,               fildesc );     \
  F_WRITE( float,     grid->dz,               fildesc );     \
  F_WRITE( float,     grid->x0,               fildesc );     \
  F_WRITE( float,     grid->y0,               fildesc );     \
  F_WRITE( float,     grid->z0,               fildesc );     \
  F_WRITE( float,     grid->cvac,             fildesc );     \
  F_WRITE( float,     grid->eps0,             fildesc );     \
  F_WRITE( float,     grid->damp,             fildesc );     \
  F_WRITE( int,       mp_rank(grid->mp),      fildesc );     \
  F_WRITE( int,       mp_nproc(grid->mp),     fildesc );     \
  /* Species parameters */                                   \
  F_WRITE( int,       sp_id,                  fildesc );     \
  F_WRITE( float,     q_m,                    fildesc );     \
} END_PRIMITIVE
 
// Note dim _MUST_ be a pointer to an int
 
#define F_WRITE_ARRAY_HEADER(psiz,ndim,dim,fildesc) BEGIN_PRIMITIVE { \
  F_WRITE( int, psiz, fildesc );                           \
  F_WRITE( int, ndim,         fildesc );                           \
   fwrite( dim, 4, ndim, fildesc );              \
} END_PRIMITIVE
 
#define F_WRITE(type,value,fildesc) BEGIN_PRIMITIVE {  \
  type __F_WRITE_tmp = (type)(value);                  \
   fwrite( &__F_WRITE_tmp, sizeof(type), 1, fildesc ); \
} END_PRIMITIVE
 
#define F_READ(type,value,fildesc) BEGIN_PRIMITIVE { \
  type __F_READ_tmp;                                 \
  fread( &__F_READ_tmp, sizeof(type), 1, fildesc );  \
  (value) = __F_READ_tmp;                            \
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
  int siztmp;
  int ii;
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
 
  /* restart save additional user variables */
  WRITE( int, ndfld, fileIO );
  WRITE( int, ndhyd, fileIO );
  WRITE( int, ndpar, fileIO );
  WRITE( int, ndhis, fileIO );
  WRITE( int, ndgrd, fileIO );
  WRITE( int, head_option, fileIO );
  WRITE( int, istride, fileIO );
  WRITE( int, jstride, fileIO );
  WRITE( int, kstride, fileIO );
  WRITE( int, stride_option, fileIO );
  WRITE( int, pstride, fileIO );
  WRITE( int, nprobe, fileIO );
  WRITE( int, ifenergies, fileIO );
  WRITE( int, block_dump, fileIO );
  WRITE( int, stepdigit, fileIO );
  WRITE( int, rankdigit, fileIO );
  siztmp = NVARFLDMX;
  fileIO.write( iffldvar, siztmp);
  siztmp = NVARHYDMX;
  fileIO.write( ifhydvar, siztmp);
  siztmp = NVARPARMX;
  fileIO.write( ifparvar, siztmp);
  siztmp = 4*NVARHISMX;
  fileIO.write( ijkprobe,      siztmp);
  siztmp = 3*NVARHISMX;
  fileIO.write( xyzprobe,      siztmp);
 
  fileIO.close();
 
  // time history dump files are closed each time a restart dump is
  // written (except on cycle 0). this is not compulsory, as files
  // could be closed at end of run, but is done for reliability in
  // continuity of dump files on restart.
  if (nfilhis > 0 && step > 0) {
    for (ii=0; ii<nfilhis; ii++) {
      fclose(fdhis[ii]);
      strcpy(&fbasehis[ii][0]," ");
    }
    nfilhis = 0;
  }
 
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
  int siztmp;
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
 
  /* restart read additional user variables */
  READ(int,    ndfld, handle );
  READ(int,    ndhyd, handle );
  READ(int,    ndpar, handle );
  READ(int,    ndhis, handle );
  READ(int,    ndgrd, handle );
  READ(int,    head_option, handle );
  READ(int,    istride, handle );
  READ(int,    jstride, handle );
  READ(int,    kstride, handle );
  READ(int,    stride_option, handle );
  READ(int,    pstride, handle );
  READ(int,    nprobe, handle );
  READ(int,    ifenergies, handle );
  READ(int,    block_dump, handle );
  READ(int,    stepdigit, handle );
  READ(int,    rankdigit, handle );
  siztmp = NVARFLDMX;
  fread( iffldvar, 4, siztmp, handle);
  siztmp = NVARHYDMX;
  fread( ifhydvar, 4, siztmp, handle);
  siztmp = NVARPARMX;
  fread( ifparvar, 4, siztmp, handle);
  siztmp = 4*NVARHISMX;
  fread( ijkprobe, 4, siztmp, handle);
  siztmp = 3*NVARHISMX;
  fread( xyzprobe, 4, siztmp, handle);
 
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
// ndfld, ndhyd, ndpar, ndhis, ndgrd, head_option,
// istride, jstride, kstride, stride_option, pstride
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
    ITEST( ndfld, "ndfld", (iarg<0 ? 0 : iarg) );
    ITEST( ndhyd, "ndhyd", (iarg<0 ? 0 : iarg) );
    ITEST( ndpar, "ndpar", (iarg<0 ? 0 : iarg) );
    ITEST( ndhis, "ndhis", (iarg<0 ? 0 : iarg) );
    ITEST( ndgrd, "ndgrd", (iarg<0 ? 0 : iarg) );
    ITEST( head_option, "head_option", (iarg<0 ? 0 : iarg) );
    ITEST( istride, "istride", (iarg<1 ? 1 : iarg) );
    ITEST( jstride, "jstride", (iarg<1 ? 1 : iarg) );
    ITEST( kstride, "kstride", (iarg<1 ? 1 : iarg) );
    ITEST( stride_option, "stride_option", (iarg<1 ? 1 : iarg) );
    ITEST( pstride, "pstride", (iarg<1 ? 1 : iarg) );
    ITEST( stepdigit, "stepdigit", (iarg<0 ? 0 : iarg) );
    ITEST( rankdigit, "rankdigit", (iarg<0 ? 0 : iarg) );
  }
}
 
/**********************************************************************/
 
#define INT int
#define REAL float
#define DOUB double
#define LONG long int
#define LONGLONG long long int
 
#include <stddef.h>
// #include <stdio.h>   /* some systems may require stdio as well as FileIo */
#include <string.h>
#include <math.h>
 
#define LENCHR 132
#define NDIGSUFF 10
#define NBUFMX 1024
 
#define MAX0(A,B) ((A) > (B) ? (A) : (B))
#define MIN0(A,B) ((A) < (B) ? (A) : (B))
 
/* static array for real buffer for writing or reading files */
/* the buffer size can be changed but may need to be limited */
/* to reduce memory size especially on blade processors.     */
 
REAL buf[NBUFMX];
 
INT nbytwrd = 4;   /* for 4-byte operation */
 
INT idhead[2];
 
/**********************************************************************/
 
void vpic_simulation::vfld_dump( const char *fbase)
 
  /* fbase - character string base name for the dump file, */
  /*         will append cycle number and processor number */
//  #include <mpi.h>
 
{
  INT nv;
  INT nvartot, offset;
  INT ier;
  INT dim[3];
  INT ii, len;
  /* number of ghost cells at each boundary, currently set to 1 */
  INT iglo=1, igup=1, jglo=1, jgup=1, kglo=1, kgup=1;
 
//    int iermpi;
 
  // begin
 
//  synchronize_rhob( field, grid );
 
  if( fbase==NULL ) {
    ERROR(("Invalid filename"));
    return;
  }
 
  /* construct the field file name, assign the file */
  strcpy(fdmpnam,fbase);
  ncyc = step;
  nrnk = mp_rank(grid->mp);
  nprc = mp_nproc(grid->mp);
  fname_append(fdmpnam, &step, &stepdigit, &nrnk, &rankdigit);
  fdfld[0] = fopen(fdmpnam,"wb");
  if(fdfld[0] == NULL) {
    printf(" *** error - could not open file %s \n",fdmpnam);
    return;
  }
 
  /* open or open-append text file for catalog of dumps and times */
  strcpy(fcatnam,fbase);
  tim = dt * float(ncyc);
  len = strlen(fcatnam);
  strcpy(&fcatnam[len],".catalog");
  ii = -1;
  fname_append(fcatnam, &step, &ii, &nrnk, &rankdigit);
  if (step<=1) {
    fdcat = fopen(fcatnam,"w");
  }
  else {
    fdcat = fopen(fcatnam,"a");
  }
  fprintf(fdcat, "%s %d %e\n",&fdmpnam[0],step,(double)tim);
  fclose(fdcat);
 
  /* dump up to 16 field variables, as specified */
  nvarfld = 0;
  for (ii=0; ii<NVARFLDMX; ii++) {
    if (iffldvar[ii] >= 0) {
      strcpy(&namfld[nvarfld][0],&names_field[ii][0]);
      strcpy(&cenfld[nvarfld][0],&center_field[ii][0]);
      nvarfld++;
    }
  }
 
  ityphead = itypfld;
  if (block_dump >= 1) ityphead = ityphead+10;
  ityp = itypfld;
  if (block_dump >= 1) ityp = ityp+10;
  nvar = nvarfld;
  vhead_init();
  npar = 0;
  pstr = 0;
  nparstr = 0;
  ipspec = 0;
  qom = 0.0;
 
  if (step<=1) {
    /* construct the header file name, assign the file */
    strcpy(fheadnam,fbase);
    len = strlen(fheadnam);
    strcpy(&fheadnam[len],".header");
    ii = -1;
    fname_append(fheadnam, &step, &ii, &nrnk, &rankdigit);
    fdheader = fopen(fheadnam,"w");
    if(fdheader == NULL) {
      printf(" *** error - could not open file %s \n",fheadnam);
      return;
    }
    //  fwr_bin_head(fdheader, &ier);
    fwr_asc_head(fdheader, &ier);
 
    idhead[0]=105;
    idhead[1] = nvarfld;
    fprintf(fdheader, "%d %d\n", idhead[0], idhead[1]);
 
    for (nv=0; nv<nvarfld; nv++) {
      fprintf(fdheader, "%s\n", &namfld[nv][0]);
    }
 
    idhead[0]=107;
    idhead[1] = nvarfld;
    fprintf(fdheader, "%d %d\n", idhead[0], idhead[1]);
 
    for (nv=0; nv<nvarfld; nv++) {
      fprintf(fdheader, "%s\n", &cenfld[nv][0]);
    }
 
    fclose(fdheader);
  }
 
  if (head_option == 1) {
    fwr_bin_head(fdfld[0], &ier);
 
    idhead[0]=105;
    idhead[1] = nvarfld;
    fwrite(&idhead[0], nbytwrd, 2, fdfld[0]);
 
    for (nv=0; nv<nvarfld; nv++) {
      fwrite(&namfld[nv][0], 1, 8, fdfld[0]);
    }
 
    idhead[0]=107;
    idhead[1] = nvarfld;
    fwrite(&idhead[0], nbytwrd, 2, fdfld[0]);
 
    for (nv=0; nv<nvarfld; nv++) {
      fwrite(&cenfld[nv][0], 1, 8, fdfld[0]);
    }
  }
  else {
    ityp = dump_type::field_dump;
    if (block_dump >= 1) ityp = ityp+20;
    F_WRITE_HEADER_V0( ityp, invalid_species_id, 0, fdfld[0] );
    dim[0] = imxstr;
    dim[1] = jmxstr;
    dim[2] = kmxstr;
    F_WRITE_ARRAY_HEADER( nvar*nbytwrd, 3, dim, fdfld[0] );
    if (step <= 1 && mp_rank(grid->mp)==0) {
      printf("Number of field variables dumped: %d\n",nvarfld);
      printf("List of field variables dumped and centering:\n");
      for (nv=0; nv<nvarfld; nv++) {
        printf(" %s  %s\n",&namfld[nv][0],&cenfld[nv][0]);
      }
    }
  }
 
  stropt = stride_option;
  nvartot = NVARFLDMX;
  offset = ((int)sizeof(field[0]))/nbytwrd;
  if (block_dump <= 0) {
    fwr_stride(&field->ex, &nvartot, &iffldvar[0], &offset,
      &imx, &jmx, &kmx, &istr, &jstr, &kstr,
      &iglo, &igup, &jglo, &jgup, &kglo, &kgup, &stropt, fdfld[0], &ier);
  }
  else {
    fwr_stride_block(&field->ex, &nvartot, &iffldvar[0], &offset,
      &imx, &jmx, &kmx, &istr, &jstr, &kstr,
      &iglo, &igup, &jglo, &jgup, &kglo, &kgup, &stropt, fdfld[0], &ier);
  }
 
//
// example of writing variables sequentially in blockwise
// rather than pointwise form.
//
// stropt = stride_option;
// nvartot = 1;
// offset = ((int)sizeof(field[0]))/nbytwrd;
//
// fwr_stride(&field->ex, &nvartot, &iffldvar[0], &offset,
//   &imx, &jmx, &kmx, &istr, &jstr, &kstr,
//   &iglo, &igup, &jglo, &jgup, &kglo, &kgup, &stropt, fdfld[0], &ier);
// fwr_stride(&field->ey, &nvartot, &iffldvar[0], &offset,
//   &imx, &jmx, &kmx, &istr, &jstr, &kstr,
//   &iglo, &igup, &jglo, &jgup, &kglo, &kgup, &stropt, fdfld[0], &ier);
// fwr_stride(&field->ez, &nvartot, &iffldvar[0], &offset,
//   &imx, &jmx, &kmx, &istr, &jstr, &kstr,
//   &iglo, &igup, &jglo, &jgup, &kglo, &kgup, &stropt, fdfld[0], &ier);
// fwr_stride(&field->div_e_err, &nvartot, &iffldvar[0], &offset,
//   &imx, &jmx, &kmx, &istr, &jstr, &kstr,
//   &iglo, &igup, &jglo, &jgup, &kglo, &kgup, &stropt, fdfld[0], &ier);
// fwr_stride(&field->cbx, &nvartot, &iffldvar[0], &offset,
//   &imx, &jmx, &kmx, &istr, &jstr, &kstr,
//   &iglo, &igup, &jglo, &jgup, &kglo, &kgup, &stropt, fdfld[0], &ier);
// fwr_stride(&field->cby, &nvartot, &iffldvar[0], &offset,
//   &imx, &jmx, &kmx, &istr, &jstr, &kstr,
//   &iglo, &igup, &jglo, &jgup, &kglo, &kgup, &stropt, fdfld[0], &ier);
// fwr_stride(&field->cbz, &nvartot, &iffldvar[0], &offset,
//   &imx, &jmx, &kmx, &istr, &jstr, &kstr,
//   &iglo, &igup, &jglo, &jgup, &kglo, &kgup, &stropt, fdfld[0], &ier);
// fwr_stride(&field->div_b_err, &nvartot, &iffldvar[0], &offset,
//   &imx, &jmx, &kmx, &istr, &jstr, &kstr,
//   &iglo, &igup, &jglo, &jgup, &kglo, &kgup, &stropt, fdfld[0], &ier);
// fwr_stride(&field->tcax, &nvartot, &iffldvar[0], &offset,
//   &imx, &jmx, &kmx, &istr, &jstr, &kstr,
//   &iglo, &igup, &jglo, &jgup, &kglo, &kgup, &stropt, fdfld[0], &ier);
// fwr_stride(&field->tcay, &nvartot, &iffldvar[0], &offset,
//   &imx, &jmx, &kmx, &istr, &jstr, &kstr,
//   &iglo, &igup, &jglo, &jgup, &kglo, &kgup, &stropt, fdfld[0], &ier);
// fwr_stride(&field->tcaz, &nvartot, &iffldvar[0], &offset,
//   &imx, &jmx, &kmx, &istr, &jstr, &kstr,
//   &iglo, &igup, &jglo, &jgup, &kglo, &kgup, &stropt, fdfld[0], &ier);
// fwr_stride(&field->rhob, &nvartot, &iffldvar[0], &offset,
//   &imx, &jmx, &kmx, &istr, &jstr, &kstr,
//   &iglo, &igup, &jglo, &jgup, &kglo, &kgup, &stropt, fdfld[0], &ier);
// fwr_stride(&field->jfx, &nvartot, &iffldvar[0], &offset,
//   &imx, &jmx, &kmx, &istr, &jstr, &kstr,
//   &iglo, &igup, &jglo, &jgup, &kglo, &kgup, &stropt, fdfld[0], &ier);
// fwr_stride(&field->jfy, &nvartot, &iffldvar[0], &offset,
//   &imx, &jmx, &kmx, &istr, &jstr, &kstr,
//   &iglo, &igup, &jglo, &jgup, &kglo, &kgup, &stropt, fdfld[0], &ier);
// fwr_stride(&field->jfz, &nvartot, &iffldvar[0], &offset,
//   &imx, &jmx, &kmx, &istr, &jstr, &kstr,
//   &iglo, &igup, &jglo, &jgup, &kglo, &kgup, &stropt, fdfld[0], &ier);
// fwr_stride(&field->rhof, &nvartot, &iffldvar[0], &offset,
//   &imx, &jmx, &kmx, &istr, &jstr, &kstr,
//   &iglo, &igup, &jglo, &jgup, &kglo, &kgup, &stropt, fdfld[0], &ier);
 
  fclose(fdfld[0]);
}
 
/**********************************************************************/
 
void vpic_simulation::vhyd_dump( const char *sp_name, const char *fbase)
 
  /* sp_name - character string, particle species name     */
  /* fbase - character string base name for the dump file, */
  /*         will append cycle number and processor number */
 
{
  species_t *sp;
  INT nv;
  INT nvartot, offset;
  INT ier;
  INT dim[3];
  INT ii, len;
  /* number of ghost cells at each boundary, currently set to 1 */
  INT iglo=1, igup=1, jglo=1, jgup=1, kglo=1, kgup=1;
 
  // begin
  // FIXME: Potential deadlocks on error returns
 
  if( sp_name==NULL ) {
    ERROR(("Invalid species name"));
    return;
  }
 
  sp = find_species_name(sp_name,species_list);
  if( sp==NULL ) {
    ERROR(("Invalid species \"%s\".",sp_name));
    return;
  }
 
  clear_hydro( hydro, grid );
  accumulate_hydro_p( hydro, sp->p, sp->np, sp->q_m, interpolator, grid );
  synchronize_hydro( hydro, grid );
 
  if( fbase==NULL ) {
    ERROR(("Invalid filename"));
    return;
  }
 
  /* construct the hyd file name, assign the file */
  strcpy(fdmpnam,fbase);
  ncyc = step;
  nrnk = mp_rank(grid->mp);
  nprc = mp_nproc(grid->mp);
  fname_append(fdmpnam, &step, &stepdigit, &nrnk, &rankdigit);
  fdhyd[0] = fopen(fdmpnam,"wb");
  if(fdhyd[0] == NULL) {
    printf(" *** error - could not open file %s \n",fdmpnam);
    return;
  }
 
  /* open or open-append text file for catalog of dumps and times */
  strcpy(fcatnam,fbase);
  tim = dt * float(ncyc);
  len = strlen(fcatnam);
  strcpy(&fcatnam[len],".catalog");
  ii = -1;
  fname_append(fcatnam, &step, &ii, &nrnk, &rankdigit);
  if (step<=1) {
    fdcat = fopen(fcatnam,"w");
  }
  else {
    fdcat = fopen(fcatnam,"a");
  }
  fprintf(fdcat, "%s %d %e\n",&fdmpnam[0],step,(double)tim);
  fclose(fdcat);
 
  /* dump up to 14 hyd variables, as specified */
  nvarhyd = 0;
  for (ii=0; ii<NVARHYDMX; ii++) {
    if (ifhydvar[ii] >= 0) {
      strcpy(&namhyd[nvarhyd][0],&names_hydro[ii][0]);
      nvarhyd++;
    }
  }
 
 
  ityphead = ityphyd;
  if (block_dump >= 1) ityphead = ityphead+10;
  ityp = ityphyd;
  if (block_dump >= 1) ityp = ityp+10;
  nvar = nvarhyd;
  vhead_init();
  npar = 0;
  pstr = 0;
  nparstr = 0;
  ipspec = 0;
  qom = 0.0;
 
  if (step<=1) {
    /* construct the header file name, assign the file */
    strcpy(fheadnam,fbase);
    len = strlen(fheadnam);
    strcpy(&fheadnam[len],".header");
    ii = -1;
    fname_append(fheadnam, &step, &ii, &nrnk, &rankdigit);
    fdheader = fopen(fheadnam,"w");
    if(fdheader == NULL) {
      printf(" *** error - could not open file %s \n",fheadnam);
      return;
    }
    //  fwr_bin_head(fdheader, &ier);
    fwr_asc_head(fdheader, &ier);
 
    idhead[0]=105;
    idhead[1] = nvarhyd;
    fprintf(fdheader, "%d %d\n", idhead[0], idhead[1]);
 
    for (nv=0; nv<nvarhyd; nv++) {
      fprintf(fdheader, "%s\n", &namhyd[nv][0]);
    }
 
    fclose(fdheader);
  }
 
  if (head_option == 1) {
    fwr_bin_head(fdhyd[0], &ier);
 
    idhead[0]=105;
    idhead[1] = nvarhyd;
    fwrite(&idhead[0], nbytwrd, 2, fdhyd[0]);
 
    for (nv=0; nv<nvarhyd; nv++) {
      fwrite(&namhyd[nv][0], 1, 8, fdhyd[0]);
    }
  }
  else {
    ityp = dump_type::hydro_dump;
    if (block_dump >= 1) ityp = ityp+20;
    F_WRITE_HEADER_V0( ityp, sp->id, sp->q_m, fdhyd[0]);
    dim[0] = imxstr;
    dim[1] = jmxstr;
    dim[2] = kmxstr;
    F_WRITE_ARRAY_HEADER( nvar*nbytwrd, 3, dim, fdhyd[0] );
    if (step <= 1 && mp_rank(grid->mp)==0) {
      printf("Number of hydro variables dumped: %d\n",nvarhyd);
      printf("List of hydro variables dumped:\n");
      for (nv=0; nv<nvarhyd; nv++) {
        printf(" %s\n",&namhyd[nv][0]);
      }
    }
  }
 
  stropt = stride_option;
  nvartot = NVARHYDMX;
  offset = ((int)sizeof(hydro[0]))/nbytwrd;
  if (block_dump <= 0) {
    fwr_stride(&hydro->jx, &nvartot, &ifhydvar[0], &offset,
      &imx, &jmx, &kmx, &istr, &jstr, &kstr,
      &iglo, &igup, &jglo, &jgup, &kglo, &kgup, &stropt, fdhyd[0], &ier);
  }
  else {
    fwr_stride_block(&hydro->jx, &nvartot, &ifhydvar[0], &offset,
      &imx, &jmx, &kmx, &istr, &jstr, &kstr,
      &iglo, &igup, &jglo, &jgup, &kglo, &kgup, &stropt, fdhyd[0], &ier);
  }
 
//
// example of writing variables sequentially in blockwise
// rather than pointwise form.
//
// stropt = stride_option;
// nvartot = 1;
// offset = ((int)sizeof(hydro[0]))/nbytwrd;;
//
// fwr_stride(&hydro->jx, &nvartot, &ifhydvar[0], &offset,
//   &imx, &jmx, &kmx, &istr, &jstr, &kstr,
//   &iglo, &igup, &jglo, &jgup, &kglo, &kgup, &stropt, fdhyd[0], &ier);
// fwr_stride(&hydro->jy, &nvartot, &ifhydvar[0], &offset,
//   &imx, &jmx, &kmx, &istr, &jstr, &kstr,
//   &iglo, &igup, &jglo, &jgup, &kglo, &kgup, &stropt, fdhyd[0], &ier);
// fwr_stride(&hydro->jz, &nvartot, &ifhydvar[0], &offset,
//   &imx, &jmx, &kmx, &istr, &jstr, &kstr,
//   &iglo, &igup, &jglo, &jgup, &kglo, &kgup, &stropt, fdhyd[0], &ier);
// fwr_stride(&hydro->rho, &nvartot, &ifhydvar[0], &offset,
//   &imx, &jmx, &kmx, &istr, &jstr, &kstr,
//   &iglo, &igup, &jglo, &jgup, &kglo, &kgup, &stropt, fdhyd[0], &ier);
// fwr_stride(&hydro->px, &nvartot, &ifhydvar[0], &offset,
//   &imx, &jmx, &kmx, &istr, &jstr, &kstr,
//   &iglo, &igup, &jglo, &jgup, &kglo, &kgup, &stropt, fdhyd[0], &ier);
// fwr_stride(&hydro->py, &nvartot, &ifhydvar[0], &offset,
//   &imx, &jmx, &kmx, &istr, &jstr, &kstr,
//   &iglo, &igup, &jglo, &jgup, &kglo, &kgup, &stropt, fdhyd[0], &ier);
// fwr_stride(&hydro->pz, &nvartot, &ifhydvar[0], &offset,
//   &imx, &jmx, &kmx, &istr, &jstr, &kstr,
//   &iglo, &igup, &jglo, &jgup, &kglo, &kgup, &stropt, fdhyd[0], &ier);
// fwr_stride(&hydro->ke, &nvartot, &ifhydvar[0], &offset,
//   &imx, &jmx, &kmx, &istr, &jstr, &kstr,
//   &iglo, &igup, &jglo, &jgup, &kglo, &kgup, &stropt, fdhyd[0], &ier);
// fwr_stride(&hydro->txx, &nvartot, &ifhydvar[0], &offset,
//   &imx, &jmx, &kmx, &istr, &jstr, &kstr,
//   &iglo, &igup, &jglo, &jgup, &kglo, &kgup, &stropt, fdhyd[0], &ier);
// fwr_stride(&hydro->tyy, &nvartot, &ifhydvar[0], &offset,
//   &imx, &jmx, &kmx, &istr, &jstr, &kstr,
//   &iglo, &igup, &jglo, &jgup, &kglo, &kgup, &stropt, fdhyd[0], &ier);
// fwr_stride(&hydro->tzz, &nvartot, &ifhydvar[0], &offset,
//   &imx, &jmx, &kmx, &istr, &jstr, &kstr,
//   &iglo, &igup, &jglo, &jgup, &kglo, &kgup, &stropt, fdhyd[0], &ier);
// fwr_stride(&hydro->tyz, &nvartot, &ifhydvar[0], &offset,
//   &imx, &jmx, &kmx, &istr, &jstr, &kstr,
//   &iglo, &igup, &jglo, &jgup, &kglo, &kgup, &stropt, fdhyd[0], &ier);
// fwr_stride(&hydro->tzx, &nvartot, &ifhydvar[0], &offset,
//   &imx, &jmx, &kmx, &istr, &jstr, &kstr,
//   &iglo, &igup, &jglo, &jgup, &kglo, &kgup, &stropt, fdhyd[0], &ier);
// fwr_stride(&hydro->txy, &nvartot, &ifhydvar[0], &offset,
//   &imx, &jmx, &kmx, &istr, &jstr, &kstr,
//   &iglo, &igup, &jglo, &jgup, &kglo, &kgup, &stropt, fdhyd[0], &ier);
 
  fclose(fdhyd[0]);
}
 
/**********************************************************************/
 
void vpic_simulation::vpar_dump( const char *sp_name, const char *fbase)
 
  /* sp_name - character string, particle species name     */
  /* fbase - character string base name for the dump file, */
  /*         will append cycle number and processor number */
 
{
  species_t *sp;
/* snell */
// int buf_start, n_buf;
//  static particle_t * ALIGNED(128) p_new=NULL;
// # define PBUF_SIZE 32768 // 1MB of particles
/* snell */
 
 
  INT nv;
  INT npr, nvartot, nstr;
  INT ier;
  INT dim[2];
  INT ii, len;
 
  // begin
  // FIXME: Potential deadlocks on error returns
 
  if( sp_name==NULL ) {
    ERROR(("Invalid species name"));
    return;
  }
 
  sp = find_species_name(sp_name,species_list);
  if( sp==NULL ) {
    ERROR(("Invalid species name \"%s\".",sp_name));
    return;
  }
 
  if( fbase==NULL ) {
    ERROR(("Invalid filename"));
    return;
  }
  /* construct the par file name, assign the file */
  strcpy(fdmpnam,fbase);
  ncyc = step;
  nrnk = mp_rank(grid->mp);
  nprc = mp_nproc(grid->mp);
  fname_append(fdmpnam, &step, &stepdigit, &nrnk, &rankdigit);
  fdpar[0] = fopen(fdmpnam,"wb");
  if(fdpar[0] == NULL) {
    printf(" *** error - could not open file %s \n",fdmpnam);
    return;
  }
 
  /* open or open-append text file for catalog of dumps and times */
  strcpy(fcatnam,fbase);
  tim = dt * float(ncyc);
  len = strlen(fcatnam);
  strcpy(&fcatnam[len],".catalog");
  ii = -1;
  fname_append(fcatnam, &step, &ii, &nrnk, &rankdigit);
  if (step<=1) {
    fdcat = fopen(fcatnam,"w");
  }
  else {
    fdcat = fopen(fcatnam,"a");
  }
  fprintf(fdcat, "%s %d %e\n",&fdmpnam[0],step,(double)tim);
  fclose(fdcat);
 
  /* dump the 8 par variables; need to add selection option */
  nvarpar = 0;
  for (ii=0; ii<NVARPARMX; ii++) {
    if (ifparvar[ii] >= 0) {
      strcpy(&nampar[nvarpar][0],&names_particle[ii][0]);
      nvarpar++;
    }
  }
 
  ityphead = ityppar;
  ityp = ityppar;
  nvar = nvarpar;
  vhead_init();
  npar = sp->np;
  pstr = pstride;
  nparstr = npar/pstr;
  ipspec = sp->id;
  qom = sp->q_m;
 
  if (step<=1) {
    /* construct the header file name, assign the file */
    strcpy(fheadnam,fbase);
    len = strlen(fheadnam);
    strcpy(&fheadnam[len],".header");
    ii = -1;
    fname_append(fheadnam, &step, &ii, &nrnk, &rankdigit);
    fdheader = fopen(fheadnam,"w");
    if(fdheader == NULL) {
      printf(" *** error - could not open file %s \n",fheadnam);
      return;
    }
    //  fwr_bin_head(fdheader, &ier);
    fwr_asc_head(fdheader, &ier);
 
    idhead[0]=105;
    idhead[1] = nvarpar;
    fprintf(fdheader, "%d %d\n", idhead[0], idhead[1]);
 
    for (nv=0; nv<nvarpar; nv++) {
      fprintf(fdheader, "%s\n", &nampar[nv][0]);
    }
 
    fclose(fdheader);
  }
 
  if (head_option == 1) {
    fwr_bin_head(fdpar[0], &ier);
 
    idhead[0]=105;
    idhead[1] = nvarpar;
    fwrite(&idhead[0], nbytwrd, 2, fdpar[0]);
 
    for (nv=0; nv<nvarpar; nv++) {
      fwrite(&nampar[nv][0], 1, 8, fdpar[0]);
    }
  }
  else {
    F_WRITE_HEADER_PAR( dump_type::particle_dump, sp->id, sp->q_m, fdpar[0] );
    dim[0] = sp->np;
    F_WRITE_ARRAY_HEADER( nvar*nbytwrd, 1, dim, fdpar[0] );
    if (step <= 1 && mp_rank(grid->mp)==0) {
      printf("Number of particle variables dumped: %d\n",nvarpar);
      printf("List of particle variables dumped:\n");
      for (nv=0; nv<nvarpar; nv++) {
        printf(" %s\n",&nampar[nv][0]);
      }
    }
//  if ( !p_new ) {
//    p_new = (particle_t * ALIGNED(128))
//      malloc_aligned(PBUF_SIZE*sizeof(p_new[0]), 128);
//    if( p_new==NULL ) {
//      ERROR(("Failed to allocate particle buffer."));
//      return;
//    }
//    F_WRITE_ARRAY_HEADER( sizeof(p_new), 1, dim, fdpar[0] );
  }
 
  npr = sp->np;
  nvartot = NVARPARMX;
  nstr = pstride;
 
  fwr_par (&sp->p->dx, &npr, &nvartot, &ifparvar[0], &nstr, fdpar[0], &ier);
 
  fclose(fdpar[0]);
 
}
 
/**********************************************************************/
 
void vpic_simulation::vhis_dump( const char *sp_name, const char *fbase )
 
  /* sp_name - character string, particle species name     */
  /* fbase - character string base name for the dump file, */
  /*         will append cycle number and processor number */
 
{
  species_t *sp = 0;
  species_t *sp_local = 0;
 
  INT nv;
  INT nvartot;
  INT np, iprbtyp;
  INT ier;
  INT dim[3];
  INT ii, len;
 
  INT ifopenhis=0, nfil;
  INT iprb, jprb, kprb;
  INT ndxbuf;
  INT nglobal;
  INT iing;
  field_t *f0;
  hydro_t *h0;
  double en_f[6], en_p;
 
  // begin
  // FIXME: Potential deadlocks on error returns
 
  if( sp_name==NULL ) {
    ERROR(("Invalid species name"));
    return;
  }
  // FIXME currently must pass the species name "none" for problems
  // that do not have particle species present and that will not have
  // hydro dumps. no error checking is done. this may need a better
  // approach.
  if ( strcmp(sp_name,"none") != 0) {
    sp = find_species_name(sp_name,species_list);
    if( sp==NULL ) {
      ERROR(("Invalid species \"%s\".",sp_name));
      return;
    }
  }
 
  if ( strcmp(sp_name,"none") != 0) {
    clear_hydro( hydro, grid );
    accumulate_hydro_p( hydro, sp->p, sp->np, sp->q_m, interpolator, grid );
    synchronize_hydro( hydro, grid );
  }
 
//  synchronize_rhob( field, grid );
 
  if( fbase==NULL ) {
    ERROR(("Invalid filename"));
    return;
  }
 
  ifopenhis = 0;
  if (nfilhis > 0) {
    for (ii=0; ii<nfilhis; ii++) {
      if (strcmp(fbase,fbasehis[ii])==0) {
        ifopenhis = 1;
        nfil = ii+1;
        break;
      }
    }
  }
 
  if (ifopenhis == 0) {
    if (nfilhis >= NFILHISMX) {
      printf(" *** error - more than NFILHISMX time history dump files.\n");
      return;
    }
 
    nfilhis++;
    nfil = nfilhis;
    strcpy(&fbasehis[nfil-1][0],fbase);
    /* construct the his file name, assign the file */
    strcpy(fdmpnam,fbase);
    ncyc = step;
    nrnk = mp_rank(grid->mp);
    nprc = mp_nproc(grid->mp);
    fname_append(fdmpnam, &step, &stepdigit, &nrnk, &rankdigit);
    fdhis[nfil-1] = fopen(fdmpnam,"wb");
    if(fdhis[nfil-1] == NULL) {
      printf(" *** error - could not open file %s \n",fdmpnam);
      return;
    }
    printf("Opening time history dump file %s\n",fdmpnam);
 
    ifopenhis = 1;
 
    /* open or open-append text file for catalog of dumps and times */
    strcpy(fcatnam,fbase);
    tim = dt * float(ncyc);
    len = strlen(fcatnam);
    strcpy(&fcatnam[len],".catalog");
    ii = -1;
    fname_append(fcatnam, &step, &ii, &nrnk, &rankdigit);
    if (step<=1) {
      fdcat = fopen(fcatnam,"w");
    }
    else {
      fdcat = fopen(fcatnam,"a");
    }
    fprintf(fdcat, "%s %d %e\n",&fdmpnam[0],step,(double)tim);
    fclose(fdcat);
 
 
    /* set the names and nvarhis (number of history variables to dump */
    nvarhis = 0;
    nglobal = 0;
    strcpy(&namhis[0][0],"step    ");
    strcpy(&namhis[1][0],"time    ");
    strcpy(&namhis[2][0],"dump    ");
    nvarhis = nvarhis+3;
    nglobal = nglobal+3;
    if (ifenergies >= 1) {
      strcpy(&namhis[3][0],"enex    ");
      strcpy(&namhis[4][0],"eney    ");
      strcpy(&namhis[5][0],"enez    ");
      strcpy(&namhis[6][0],"enbx    ");
      strcpy(&namhis[7][0],"enby    ");
      strcpy(&namhis[8][0],"enbz    ");
      nvarhis = nvarhis+6;
      nglobal = nglobal+6;
      LIST_FOR_EACH(sp_local,species_list) {
        if (sp_local->id < 10) {sprintf(&namhis[nvarhis][0],"ensp%i   ",sp_local->id);}
        else if (sp_local->id < 100) {sprintf(&namhis[nvarhis][0],"ensp%i  ",sp_local->id);}
        else if (sp_local->id < 1000) {sprintf(&namhis[nvarhis][0],"ensp%i ",sp_local->id);}
        else {sprintf(&namhis[nvarhis][0],"ensp%i",sp_local->id);}
        nvarhis++;
        nglobal++;
      }
    }
    if ( (nprobe+nglobal) > NVARHISMX) {
      if (mp_rank(grid->mp)==0) {
        printf(" *** error - more than NVARHISMX time history globals plus probes.\n");
      }
    }
    for (ii=0; ii<nprobe; ii++) {
      iing = ii + nglobal;
      if      (ijkprobe[ii][0]==1) {vhis_nam(&namhis[iing][0],"ex",&iing);}
      else if (ijkprobe[ii][0]==2) {vhis_nam(&namhis[iing][0],"ey",&iing);}
      else if (ijkprobe[ii][0]==3) {vhis_nam(&namhis[iing][0],"ez",&iing);}
      else if (ijkprobe[ii][0]==4) {vhis_nam(&namhis[iing][0],"dive",&iing);}
      else if (ijkprobe[ii][0]==5) {vhis_nam(&namhis[iing][0],"cbx",&iing);}
      else if (ijkprobe[ii][0]==6) {vhis_nam(&namhis[iing][0],"cby",&iing);}
      else if (ijkprobe[ii][0]==7) {vhis_nam(&namhis[iing][0],"cbz",&iing);}
      else if (ijkprobe[ii][0]==8) {vhis_nam(&namhis[iing][0],"divb",&iing);}
      else if (ijkprobe[ii][0]==9) {vhis_nam(&namhis[iing][0],"tcax",&iing);}
      else if (ijkprobe[ii][0]==10) {vhis_nam(&namhis[iing][0],"tcay",&iing);}
      else if (ijkprobe[ii][0]==11) {vhis_nam(&namhis[iing][0],"tcaz",&iing);}
      else if (ijkprobe[ii][0]==12) {vhis_nam(&namhis[iing][0],"rhob",&iing);}
      else if (ijkprobe[ii][0]==13) {vhis_nam(&namhis[iing][0],"jfx",&iing);}
      else if (ijkprobe[ii][0]==14) {vhis_nam(&namhis[iing][0],"jfy",&iing);}
      else if (ijkprobe[ii][0]==15) {vhis_nam(&namhis[iing][0],"jfz",&iing);}
      else if (ijkprobe[ii][0]==16) {vhis_nam(&namhis[iing][0],"rhof",&iing);}
      else if (ijkprobe[ii][0]==17) {vhis_nam(&namhis[iing][0],"jx",&iing);}
      else if (ijkprobe[ii][0]==18) {vhis_nam(&namhis[iing][0],"jy",&iing);}
      else if (ijkprobe[ii][0]==19) {vhis_nam(&namhis[iing][0],"jz",&iing);}
      else if (ijkprobe[ii][0]==20) {vhis_nam(&namhis[iing][0],"rho",&iing);}
      else if (ijkprobe[ii][0]==21) {vhis_nam(&namhis[iing][0],"px",&iing);}
      else if (ijkprobe[ii][0]==22) {vhis_nam(&namhis[iing][0],"py",&iing);}
      else if (ijkprobe[ii][0]==23) {vhis_nam(&namhis[iing][0],"pz",&iing);}
      else if (ijkprobe[ii][0]==24) {vhis_nam(&namhis[iing][0],"ke",&iing);}
      else if (ijkprobe[ii][0]==25) {vhis_nam(&namhis[iing][0],"txx",&iing);}
      else if (ijkprobe[ii][0]==26) {vhis_nam(&namhis[iing][0],"tyy",&iing);}
      else if (ijkprobe[ii][0]==27) {vhis_nam(&namhis[iing][0],"tzz",&iing);}
      else if (ijkprobe[ii][0]==28) {vhis_nam(&namhis[iing][0],"tyz",&iing);}
      else if (ijkprobe[ii][0]==29) {vhis_nam(&namhis[iing][0],"tzx",&iing);}
      else if (ijkprobe[ii][0]==30) {vhis_nam(&namhis[iing][0],"txy",&iing);}
      nvarhis++;
    }
 
    ityphead = ityphis;
    ityp = ityphis;
    nvar = nvarhis;
    vhead_init();
    npar = 0;
    pstr = 0;
    nparstr = 0;
    ipspec = 0;
    if ( strcmp(sp_name,"none") != 0) {
      qom = sp->q_m;
    }
    else {
      qom = 0.0;
    }
 
    /* history files do not use the short header at present */
    if (head_option == 1 || head_option == 0) {
      fwr_bin_head(fdhis[nfil-1], &ier);
 
      idhead[0]=105;
      idhead[1] = nvarhis;
      fwrite(&idhead[0], nbytwrd, 2, fdhis[nfil-1]);
 
      for (nv=0; nv<nvarhis; nv++) {
        fwrite(&namhis[nv][0], 1, 8, fdhis[nfil-1]);
      }
    }
    else {
      F_WRITE_HEADER_V0( dump_type::history_dump, invalid_species_id, 0, fdhis[nfil-1] );
      dim[0] = imxstr;
      dim[1] = jmxstr;
      dim[2] = kmxstr;
      F_WRITE_ARRAY_HEADER( nvar*nbytwrd, 3, dim, fdhis[nfil-1] );
      if (step <= 1 && mp_rank(grid->mp)==0) {
        printf("Number of time history variables dumped: %d\n",nvarhis);
        printf("List of time history variables dumped:\n");
        for (nv=0; nv<nvarhis; nv++) {
          printf(" %s\n",&namhis[nv][0]);
        }
      }
    }
  }
 
  nvartot = NVARHISMX;
  ndxbuf = -1;   /* allow for c zero indexing */
 
  /* dump first three variables: step, time, and dump number */
  ndxbuf++;
  buf[ndxbuf] = (float)step;
  ndxbuf++;
  buf[ndxbuf] = (grid->dt)*((float)step);;
  ndxbuf++;
  buf[ndxbuf] = (float)(step/ndhis + 1);
 
  /* dump all global variables if requested (ifenergies>=1) */
  if (ifenergies >= 1) {
//    energy_f( en_f, field, material_coefficient, grid );
    field_advance->method->energy_f( en_f, field_advance->f, field_advance->m, field_advance->g );
    ndxbuf++;
    if (mp_rank(grid->mp)==0) buf[ndxbuf] = (float)en_f[0];
    else buf[ndxbuf] = 0.0;
    ndxbuf++;
    if (mp_rank(grid->mp)==0) buf[ndxbuf] = (float)en_f[1];
    else buf[ndxbuf] = 0.0;
    ndxbuf++;
    if (mp_rank(grid->mp)==0) buf[ndxbuf] = (float)en_f[2];
    else buf[ndxbuf] = 0.0;
    ndxbuf++;
    if (mp_rank(grid->mp)==0) buf[ndxbuf] = (float)en_f[3];
    else buf[ndxbuf] = 0.0;
    ndxbuf++;
    if (mp_rank(grid->mp)==0) buf[ndxbuf] = (float)en_f[4];
    else buf[ndxbuf] = 0.0;
    ndxbuf++;
    if (mp_rank(grid->mp)==0) buf[ndxbuf] = (float)en_f[5];
    else buf[ndxbuf] = 0.0;
    LIST_FOR_EACH(sp,species_list) {
      en_p = energy_p( sp->p, sp->np, sp->q_m, interpolator, grid );
      ndxbuf++;
      if (mp_rank(grid->mp)==0) buf[ndxbuf] = (float)en_p;
      else buf[ndxbuf] = 0.0;
    }
  }
 
  for (np=0; np<nprobe; np++) {
    ndxbuf++;
    iprbtyp = ijkprobe[np][0];
    if (ijkprobe[np][1]<=0 && ijkprobe[np][2]<=0 && ijkprobe[np][3]<=0) {
      iprbtyp = 0;
      buf[ndxbuf] = 0.0;
    }
    else {
      iprb = ijkprobe[np][1];
      jprb = ijkprobe[np][2];
      kprb = ijkprobe[np][3];
      ii = iprb + (jprb)*imx + (kprb)*imx*jmx;
      if (iprbtyp <= 16) {
        f0 = &field[ii];
        if      (iprbtyp ==  1) buf[ndxbuf] = f0->ex;
        else if (iprbtyp ==  2) buf[ndxbuf] = f0->ey;
        else if (iprbtyp ==  3) buf[ndxbuf] = f0->ez;
        else if (iprbtyp ==  4) buf[ndxbuf] = f0->div_e_err;
        else if (iprbtyp ==  5) buf[ndxbuf] = f0->cbx;
        else if (iprbtyp ==  6) buf[ndxbuf] = f0->cby;
        else if (iprbtyp ==  7) buf[ndxbuf] = f0->cbz;
        else if (iprbtyp ==  8) buf[ndxbuf] = f0->div_b_err;
        else if (iprbtyp ==  9) buf[ndxbuf] = f0->tcax;
        else if (iprbtyp == 10) buf[ndxbuf] = f0->tcay;
        else if (iprbtyp == 11) buf[ndxbuf] = f0->tcaz;
        else if (iprbtyp == 12) buf[ndxbuf] = f0->rhob;
        else if (iprbtyp == 13) buf[ndxbuf] = f0->jfx;
        else if (iprbtyp == 14) buf[ndxbuf] = f0->jfy;
        else if (iprbtyp == 15) buf[ndxbuf] = f0->jfz;
        else if (iprbtyp == 16) buf[ndxbuf] = f0->rhof;
      }
      else if (iprbtyp <= 30 && strcmp(sp_name,"none") != 0) {
        h0 = &hydro[ii];
        if      (iprbtyp == 17) buf[ndxbuf] = h0->jx;
        else if (iprbtyp == 18) buf[ndxbuf] = h0->jy;
        else if (iprbtyp == 19) buf[ndxbuf] = h0->jz;
        else if (iprbtyp == 20) buf[ndxbuf] = h0->rho;
        else if (iprbtyp == 21) buf[ndxbuf] = h0->px;
        else if (iprbtyp == 22) buf[ndxbuf] = h0->py;
        else if (iprbtyp == 23) buf[ndxbuf] = h0->pz;
        else if (iprbtyp == 24) buf[ndxbuf] = h0->ke;
        else if (iprbtyp == 25) buf[ndxbuf] = h0->txx;
        else if (iprbtyp == 26) buf[ndxbuf] = h0->tyy;
        else if (iprbtyp == 27) buf[ndxbuf] = h0->tzz;
        else if (iprbtyp == 28) buf[ndxbuf] = h0->tyz;
        else if (iprbtyp == 29) buf[ndxbuf] = h0->tzx;
        else if (iprbtyp == 30) buf[ndxbuf] = h0->txy;
      }
      else {printf(" *** error - unrecognized probe number %d\n",iprbtyp); }
    }
    if (ndxbuf >= NBUFMX-1) {
      fwr_buf(&ndxbuf, fdhis[nfil-1], &ier);
    }
  }
  /* write any remaining unwritten data to the file */
  if (ndxbuf >= 0) {
    fwr_buf(&ndxbuf, fdhis[nfil-1], &ier);
  }
 
}
 
/**********************************************************************/
 
void vpic_simulation::vgrd_dump( const char *fbase )
 
  /* fbase - character string base name for the dump file, */
  /*         will append cycle number and processor number */
 
 
/* this routine does not yet do anything substantive since we do not  */
/* know exactly what it is intended to do. the basic idea is that     */
/* the grid dump routine will save all the master static data such    */
/* as the geometry and materials at the beginning of the simulation   */
/* and perhaps at restarts. only the changing dynamic data will be    */
/* dumped to the other repeated files (field, hydro, particle, etc.)  */
 
{
  INT nv;
  INT nvartot, offset;
  INT ier;
  INT dim[3];
  INT ii, len;
  /* number of ghost cells at each boundary, currently set to 1 */
  INT iglo=1, igup=1, jglo=1, jgup=1, kglo=1, kgup=1;
 
  /* NOTE: grid variables currently start at word 17 in field struct */
  INT ifxyz=1, nvstrt=17;
 
  // begin
 
  if( fbase==NULL ) {
    ERROR(("Invalid filename"));
    return;
  }
 
  /* construct the grid dump file name, assign the file */
  strcpy(fdmpnam,fbase);
  ncyc = step;
  nrnk = mp_rank(grid->mp);
  nprc = mp_nproc(grid->mp);
  fname_append(fdmpnam, &step, &stepdigit, &nrnk, &rankdigit);
  fdgrd[0] = fopen(fdmpnam,"wb");
  if(fdgrd[0] == NULL) {
    printf(" *** error - could not open file %s \n",fdmpnam);
    return;
  }
 
  /* open or open-append text file for catalog of dumps and times */
  strcpy(fcatnam,fbase);
  tim = dt * float(ncyc);
  len = strlen(fcatnam);
  strcpy(&fcatnam[len],".catalog");
  ii = -1;
  fname_append(fcatnam, &step, &ii, &nrnk, &rankdigit);
  if (step<=1) {
    fdcat = fopen(fcatnam,"w");
  }
  else {
    fdcat = fopen(fcatnam,"a");
  }
  fprintf(fdcat, "%s %d %e\n",&fdmpnam[0],step,(double)tim);
  fclose(fdcat);
 
  /* dump the grid variables to a separate dump file */
  nvargrd = 0;
  if (ifxyz == 1) {
    strcpy(&namgrd[0][0],"xcoord  ");
    strcpy(&namgrd[1][0],"ycoord  ");
    strcpy(&namgrd[2][0],"zcoord  ");
    nvargrd = 3;
  }
  for (ii=0; ii<NVARGRDMX; ii++) {
    strcpy(&namgrd[nvargrd][0],&names_grid[ii][0]);
    nvargrd++;
  }
 
  ityphead = itypgrd;
  if (block_dump >= 1) ityphead = ityphead+10;
  ityp = itypgrd;
  if (block_dump >= 1) ityp = ityp+10;
  nvar = nvargrd;
  vhead_init();
  npar = 0;
  pstr = 0;
  nparstr = 0;
  ipspec = 0;
  qom = 0.0;
 
  if (step<=1) {
    /* construct the header file name, assign the file */
    strcpy(fheadnam,fbase);
    len = strlen(fheadnam);
    strcpy(&fheadnam[len],".header");
    ii = -1;
    fname_append(fheadnam, &step, &ii, &nrnk, &rankdigit);
    fdheader = fopen(fheadnam,"w");
    if(fdheader == NULL) {
      printf(" *** error - could not open file %s \n",fheadnam);
      return;
    }
    //  fwr_bin_head(fdheader, &ier);
    fwr_asc_head(fdheader, &ier);
 
    idhead[0]=105;
    idhead[1] = nvargrd;
    fprintf(fdheader, "%d %d\n", idhead[0], idhead[1]);
 
    for (nv=0; nv<nvargrd; nv++) {
      fprintf(fdheader, "%s\n", &namgrd[nv][0]);
    }
 
    fclose(fdheader);
  }
 
  if (head_option == 1) {
    fwr_bin_head(fdgrd[0], &ier);
 
    idhead[0]=105;
    idhead[1] = nvargrd;
    fwrite(&idhead[0], nbytwrd, 2, fdgrd[0]);
 
    for (nv=0; nv<nvargrd; nv++) {
      fwrite(&namgrd[nv][0], 1, 8, fdgrd[0]);
    }
  }
  else {
    ityp = dump_type::grid_dump;
    if (block_dump >= 1) ityp = ityp+20;
    F_WRITE_HEADER_V0( ityp, invalid_species_id, 0, fdgrd[0] );
    dim[0] = imxstr;
    dim[1] = jmxstr;
    dim[2] = kmxstr;
    F_WRITE_ARRAY_HEADER( nvar*nbytwrd, 3, dim, fdgrd[0] );
    if (step <= 1 && mp_rank(grid->mp)==0) {
      printf("Number of grid variables dumped: %d\n",nvargrd);
      printf("List of grid variables dumped:\n");
      for (nv=0; nv<nvargrd; nv++) {
        printf(" %s\n",&namgrd[nv][0]);
      }
    }
  }
 
  stropt = MIN0(2,stride_option);
  nvartot = nvargrd;
  offset = ((int)sizeof(field[0]))/nbytwrd;
 
  fwr_stride_grid(&field->ex, &nvartot, &ifxyz, &nvstrt, &offset,
    &imx, &jmx, &kmx, &istr, &jstr, &kstr,
    &iglo, &igup, &jglo, &jgup, &kglo, &kgup,
    &xmn, &ymn, &zmn, &dx, &dy, &dz,
    &stropt, fdgrd[0], &ier);
 
  fclose(fdgrd[0]);
}
 
/**********************************************************************/
 
void vpic_simulation::vhead_init( )
 
/* set the static header quantities that are universl to all */
/* types of dump files. the other quantities including       */
/* ityphead, ityp, nvar, and some particle quantities depend */
/* on the type of dump and must be set individually.         */
 
{
  /* number of ghost cells at each boundary, currently set to 1 */
  INT iglo=1, igup=1, jglo=1, jgup=1, kglo=1, kgup=1;
 
  static int ical=0;
 
  tim = dt * float(ncyc);
 
  ical++;
  if (ical>1) return;
 
  /* must set ityphead, ityp, and nvar for each individual dump */
 
  imx = grid->nx+2;
  jmx = grid->ny+2;
  kmx = grid->nz+2;
  istr = istride;
  jstr = jstride;
  kstr = kstride;
  imxstr = (imx-iglo-igup)/istr + iglo + igup;
  jmxstr = (jmx-jglo-jgup)/jstr + jglo + jgup;
  kmxstr = (kmx-kglo-kgup)/kstr + kglo + kgup;
  ndmp = 1;
  npar = 0;
  pstr = 0;
  nparstr = 0;
  ipspec = 0;
  dt = grid->dt;
  tim = dt * float(ncyc);
  cvac = grid->cvac;
  eps0 = grid->eps0;
  damp = grid->damp;
  qom = 0.;
  dx = grid->dx;
  dy = grid->dy;
  dz = grid->dz;
  xmn = grid->x0;
  xmx = xmn + dx * float(imx-iglo-igup);
  ymn = grid->y0;
  ymx = ymn + dy * float(jmx-jglo-jgup);
  zmn = grid->z0;
  zmx = zmn + dz * float(kmx-kglo-kgup);
  dxstr = dx*float(istr);
  dystr = dy*float(jstr);
  dzstr = dz*float(kstr);
}
 
/**********************************************************************/
 
void vpic_simulation::vhis_probe( const char *prbname,
      float xprb, float yprb, float zprb)
 
  /* prbname - character string, identifying name of probe (ex, etc.) */
  /* xprb, yprb, zprb - x,y,z physical coordinates of the probe       */
 
{
  int iprb, jprb, kprb;
  int ip, ipmx;
  int ioffdomain;
 
  // begin
 
  if (nprobe >= NVARHISMX) {
    if (mp_rank(grid->mp)==0) {
      printf(" *** error - more than NVARHISMX time history probes.\n");
      printf("             the probe %s is ignored.\n", prbname);
    }
    return;
  }
  if (strcmp(prbname,"field")==0 && (nprobe+6) >= NVARHISMX) {
    if (mp_rank(grid->mp)==0) {
      printf(" *** error - more than NVARHISMX time history probes.\n");
      printf("             the probe %s is ignored.\n", prbname);
    }
    return;
  }
  if (strcmp(prbname,"hydro")==0 && (nprobe+8) >= NVARHISMX) {
    if (mp_rank(grid->mp)==0) {
      printf(" *** error - more than NVARHISMX time history probes.\n");
      printf("             the probe %s is ignored.\n", prbname);
    }
    return;
  }
 
  if (grid->dx <= 0.0) {
    printf(" *** error - must define grid before history probes.\n");
    return;
  }
 
  ioffdomain = 0;
 
  /* if probe is not in this domain, turn off sentinels. it is */
  /* assumed that probe will be found in one and only one dom  */
  /* and sentinels will be turned on; needs to be verified.    */
 
  if (xprb < grid->x0 || xprb >= grid->x0 + grid->dx * float(grid->nx) ||
      yprb < grid->y0 || yprb >= grid->y0 + grid->dy * float(grid->ny) ||
      zprb < grid->z0 || zprb >= grid->z0 + grid->dz * float(grid->nz) ) {
    ioffdomain = 1;
    ipmx = 1;
    if (strcmp(prbname,"field")==0) ipmx=6;
    if (strcmp(prbname,"hydro")==0) ipmx=8;
    for (ip=1; ip<=ipmx; ip++) {
      ijkprobe[nprobe][1] = 0;
      ijkprobe[nprobe][2] = 0;
      ijkprobe[nprobe][3] = 0;
    }
  }
  iprb = 1 + (int)((xprb - grid->x0) / grid->dx);
  jprb = 1 + (int)((yprb - grid->y0) / grid->dy);
  kprb = 1 + (int)((zprb - grid->z0) / grid->dz);
  iprb = MIN0(MAX0(iprb,1),grid->nx );
  jprb = MIN0(MAX0(jprb,1),grid->ny );
  kprb = MIN0(MAX0(kprb,1),grid->nz );
 
  /* default set of six field probes for this cell */
  if (strcmp(prbname,"field")==0) {
    ipmx = 6;
    for (ip=1; ip<=ipmx; ip++) {
      if      (ip==1) {ijkprobe[nprobe][0] = 1;}
      else if (ip==2) {ijkprobe[nprobe][0] = 2;}
      else if (ip==3) {ijkprobe[nprobe][0] = 3;}
      else if (ip==4) {ijkprobe[nprobe][0] = 5;}
      else if (ip==5) {ijkprobe[nprobe][0] = 6;}
      else if (ip==6) {ijkprobe[nprobe][0] = 7;}
      if (ioffdomain==0) ijkprobe[nprobe][1] = iprb;
      if (ioffdomain==0) ijkprobe[nprobe][2] = jprb;
      if (ioffdomain==0) ijkprobe[nprobe][3] = kprb;
      nprobe++;
    }
    printf("Time history six field probes on domain %d, ",
     mp_rank(grid->mp));
    printf(" at %d %d %d\n",ijkprobe[nprobe-1][1],
     ijkprobe[nprobe-1][2],ijkprobe[nprobe-1][3]);
    return;
  }
 
  /* default set of eight hydro probes for this cell */
  if (strcmp(prbname,"hydro")==0) {
    ipmx = 8;
    for (ip=1; ip<=ipmx; ip++) {
      if      (ip==1) {ijkprobe[nprobe][0] = 17;}
      else if (ip==2) {ijkprobe[nprobe][0] = 18;}
      else if (ip==3) {ijkprobe[nprobe][0] = 19;}
      else if (ip==4) {ijkprobe[nprobe][0] = 20;}
      else if (ip==5) {ijkprobe[nprobe][0] = 21;}
      else if (ip==6) {ijkprobe[nprobe][0] = 22;}
      else if (ip==7) {ijkprobe[nprobe][0] = 23;}
      else if (ip==8) {ijkprobe[nprobe][0] = 24;}
      if (ioffdomain==0) ijkprobe[nprobe][1] = iprb;
      if (ioffdomain==0) ijkprobe[nprobe][2] = jprb;
      if (ioffdomain==0) ijkprobe[nprobe][3] = kprb;
      nprobe++;
    }
    printf("Time history eight hydro probes on domain %d, ",
     mp_rank(grid->mp));
    printf(" at %d %d %d\n",ijkprobe[nprobe-1][1],
     ijkprobe[nprobe-1][2],ijkprobe[nprobe-1][3]);
    return;
  }
 
  /* field variable probes, 1-16 */
  if      (strcmp(prbname,"ex")==0) {ijkprobe[nprobe][0] = 1;}
  else if (strcmp(prbname,"ey")==0) {ijkprobe[nprobe][0] = 2;}
  else if (strcmp(prbname,"ez")==0) {ijkprobe[nprobe][0] = 3;}
  else if (strcmp(prbname,"div_e_err")==0) {ijkprobe[nprobe][0] = 4;}
  else if (strcmp(prbname,"cbx")==0) {ijkprobe[nprobe][0] = 5;}
  else if (strcmp(prbname,"cby")==0) {ijkprobe[nprobe][0] = 6;}
  else if (strcmp(prbname,"cbz")==0) {ijkprobe[nprobe][0] = 7;}
  else if (strcmp(prbname,"div_b_err")==0) {ijkprobe[nprobe][0] = 8;}
  else if (strcmp(prbname,"tcax")==0) {ijkprobe[nprobe][0] = 9;}
  else if (strcmp(prbname,"tcay")==0) {ijkprobe[nprobe][0] = 10;}
  else if (strcmp(prbname,"tcaz")==0) {ijkprobe[nprobe][0] = 11;}
  else if (strcmp(prbname,"rhob")==0) {ijkprobe[nprobe][0] = 12;}
  else if (strcmp(prbname,"jfx")==0) {ijkprobe[nprobe][0] = 13;}
  else if (strcmp(prbname,"jfy")==0) {ijkprobe[nprobe][0] = 14;}
  else if (strcmp(prbname,"jfz")==0) {ijkprobe[nprobe][0] = 15;}
  else if (strcmp(prbname,"rhof")==0) {ijkprobe[nprobe][0] = 16;}
  /* hydro variable probes, 17-30 */
  else if (strcmp(prbname,"jx")==0) {ijkprobe[nprobe][0] = 17;}
  else if (strcmp(prbname,"jy")==0) {ijkprobe[nprobe][0] = 18;}
  else if (strcmp(prbname,"jz")==0) {ijkprobe[nprobe][0] = 19;}
  else if (strcmp(prbname,"rho")==0) {ijkprobe[nprobe][0] = 20;}
  else if (strcmp(prbname,"px")==0) {ijkprobe[nprobe][0] = 21;}
  else if (strcmp(prbname,"py")==0) {ijkprobe[nprobe][0] = 22;}
  else if (strcmp(prbname,"pz")==0) {ijkprobe[nprobe][0] = 23;}
  else if (strcmp(prbname,"ke")==0) {ijkprobe[nprobe][0] = 24;}
  else if (strcmp(prbname,"txx")==0) {ijkprobe[nprobe][0] = 25;}
  else if (strcmp(prbname,"tyy")==0) {ijkprobe[nprobe][0] = 26;}
  else if (strcmp(prbname,"tzz")==0) {ijkprobe[nprobe][0] = 27;}
  else if (strcmp(prbname,"tyz")==0) {ijkprobe[nprobe][0] = 28;}
  else if (strcmp(prbname,"tzx")==0) {ijkprobe[nprobe][0] = 29;}
  else if (strcmp(prbname,"txy")==0) {ijkprobe[nprobe][0] = 30;}
 
  if (ioffdomain==0) ijkprobe[nprobe][1] = iprb;
  if (ioffdomain==0) ijkprobe[nprobe][2] = jprb;
  if (ioffdomain==0) ijkprobe[nprobe][3] = kprb;
 
  /* check for illegal probe type, kill illegal probe */
  if (ijkprobe[nprobe][0] < 1 || ijkprobe[nprobe][0] > 30) {
    if (mp_rank(grid->mp)==0) {
      printf(" *** error - illegal probe type %d, %s\n",
       ijkprobe[nprobe][0], prbname);
      var_ls( );   /* print help list of variable names */
    }
    ijkprobe[nprobe][0] = 0;
    ijkprobe[nprobe][1] = 0;
    ijkprobe[nprobe][2] = 0;
    ijkprobe[nprobe][3] = 0;
  }
  else {
    printf("Time history probe %s ",prbname);
    printf(" on domain %d, ",mp_rank(grid->mp));
    printf(" at %d %d %d\n",ijkprobe[nprobe][1],
     ijkprobe[nprobe][2],ijkprobe[nprobe][3]);
    nprobe++;
  }
 
}
 
/**********************************************************************/
 
void vpic_simulation::vhis_probe( const int prbtyp,
      float xprb, float yprb, float zprb)
 
  /* prbtyp - integer type of probe, currently 1 to 30               */
  /* xprb, yprb, zprb - x,y,z physical coordinates of the probe       */
 
{
  int iprb, jprb, kprb;
  int ioffdomain;
 
  // begin
 
  if (nprobe >= NVARHISMX) {
    if (mp_rank(grid->mp)==0) {
      printf(" *** error - more than NVARHISMX time history probes,\n");
      printf("             the probe %d is ignored.\n", prbtyp);
    }
    return;
  }
 
  if (grid->dx <= 0.0) {
    printf(" *** error - must define grid before history probes.\n");
    return;
  }
 
  ioffdomain = 0;
 
  /* if the probe is not in this domain, turn off sentinels. it is */
  /* assumed that probe will be found in one and only one of the   */
  /* domains and will be turned on there; needs to be verified.    */
 
  if (xprb < grid->x0 || xprb >= grid->x0 + grid->dx * float(grid->nx) ||
      yprb < grid->y0 || yprb >= grid->y0 + grid->dy * float(grid->ny) ||
      zprb < grid->z0 || zprb >= grid->z0 + grid->dz * float(grid->nz) ) {
    ioffdomain = 1;
    ijkprobe[nprobe][1] = 0;
    ijkprobe[nprobe][2] = 0;
    ijkprobe[nprobe][3] = 0;
  }
  iprb = 1 + (int)((xprb - grid->x0) / grid->dx);
  jprb = 1 + (int)((yprb - grid->y0) / grid->dy);
  kprb = 1 + (int)((zprb - grid->z0) / grid->dz);
  iprb = MIN0(MAX0(iprb,1),grid->nx );
  jprb = MIN0(MAX0(jprb,1),grid->ny );
  kprb = MIN0(MAX0(kprb,1),grid->nz );
 
  ijkprobe[nprobe][0] = prbtyp;
 
  if (ioffdomain==0) ijkprobe[nprobe][1] = iprb;
  if (ioffdomain==0) ijkprobe[nprobe][2] = jprb;
  if (ioffdomain==0) ijkprobe[nprobe][3] = kprb;
 
  /* check for illegal probe type, kill illegal probe */
  if (ijkprobe[nprobe][0] < 1 || ijkprobe[nprobe][0] > 30) {
    if (mp_rank(grid->mp)==0) {
      printf(" *** error - illegal probe type %d\n",ijkprobe[nprobe][0]);
    }
    ijkprobe[nprobe][0] = 0;
    ijkprobe[nprobe][1] = 0;
    ijkprobe[nprobe][2] = 0;
    ijkprobe[nprobe][3] = 0;
  }
  else {
    printf("Time history probe %d ",ijkprobe[nprobe][0]);
    printf(" on domain %d, ",mp_rank(grid->mp));
    printf(" at %d %d %d\n",ijkprobe[nprobe][1],
     ijkprobe[nprobe][2],ijkprobe[nprobe][3]);
    nprobe++;
  }
 
}
 
/**********************************************************************/
 
void vpic_simulation::vhis_energies( const char *optname)
 
  /* optname - character string, "on" or "off", default is "on" */
 
/* turn on the output of global field and particle energies to the */
/* time history dump files. this requires synchronization and thus */
/* could slow the code significantly if frequent history dumps.    */
 
{
 
  // begin
 
  ifenergies = 1;
  if (strcmp(optname,"off")==0) ifenergies = 0;
 
}
 
/**********************************************************************/
 
void vpic_simulation::vhis_nam ( char *cnam, const char *cbase, const int *prbnum)
 
  /* cnam - character string, returned with constructed name */
  /* cbase - character string, 1 to 4 character base name    */
  /* prbnum - integer number to be appended to base name     */
 
/* construct a name for this history probe by appending a four-digit */
/* probe number to the base name, with zero fill inserted as needed. */
/* the base name must be 1 to 4 characters long, otherwise error.    */
 
{
  int ilen;
 
  // begin
 
  strcpy(&cnam[0],"        ");
  ilen = strlen(&cbase[0]);
  if (ilen>4) {
    printf(" *** error - time history variable name too long.\n");
    return;
  }
  strcpy(&cnam[0],&cbase[0]);
  sprintf(&cnam[ilen],"%04d",(*prbnum + 1));
  ilen=ilen+4;
  if      (ilen==4) strcpy(&cnam[ilen],"    ");
  else if (ilen==5) strcpy(&cnam[ilen],"   ");
  else if (ilen==6) strcpy(&cnam[ilen],"  ");
  else if (ilen==7) strcpy(&cnam[ilen]," ");
 
}
 
/**********************************************************************/
 
void vpic_simulation::var_off( const char *varname)
 
  /* varname - character string, name of variable to be turned off */
 
/* turn off dumping of a field, hydro, or particle variable varname */
 
{
 
  // begin
 
  if (mp_rank(grid->mp)==0) {
    printf("Turning off dump variable %s.\n",varname);
  }
 
  /* field variables, 16 */
  if      (strcmp(varname,"ex")==0) {iffldvar[0] = -1;}
  else if (strcmp(varname,"ey")==0) {iffldvar[1] = -1;}
  else if (strcmp(varname,"ez")==0) {iffldvar[2] = -1;}
  else if (strcmp(varname,"div_e_err")==0) {iffldvar[3] = -1;}
  else if (strcmp(varname,"cbx")==0) {iffldvar[4] = -1;}
  else if (strcmp(varname,"cby")==0) {iffldvar[5] = -1;}
  else if (strcmp(varname,"cbz")==0) {iffldvar[6] = -1;}
  else if (strcmp(varname,"div_b_err")==0) {iffldvar[7] = -1;}
  else if (strcmp(varname,"tcax")==0) {iffldvar[8] = -1;}
  else if (strcmp(varname,"tcay")==0) {iffldvar[9] = -1;}
  else if (strcmp(varname,"tcaz")==0) {iffldvar[10] = -1;}
  else if (strcmp(varname,"rhob")==0) {iffldvar[11] = -1;}
  else if (strcmp(varname,"jfx")==0) {iffldvar[12] = -1;}
  else if (strcmp(varname,"jfy")==0) {iffldvar[13] = -1;}
  else if (strcmp(varname,"jfz")==0) {iffldvar[14] = -1;}
  else if (strcmp(varname,"rhof")==0) {iffldvar[15] = -1;}
  /* hydro variables, 14 */
  else if (strcmp(varname,"jx")==0) {ifhydvar[0] = -1;}
  else if (strcmp(varname,"jy")==0) {ifhydvar[1] = -1;}
  else if (strcmp(varname,"jz")==0) {ifhydvar[2] = -1;}
  else if (strcmp(varname,"rho")==0) {ifhydvar[3] = -1;}
  else if (strcmp(varname,"px")==0) {ifhydvar[4] = -1;}
  else if (strcmp(varname,"py")==0) {ifhydvar[5] = -1;}
  else if (strcmp(varname,"pz")==0) {ifhydvar[6] = -1;}
  else if (strcmp(varname,"ke")==0) {ifhydvar[7] = -1;}
  else if (strcmp(varname,"txx")==0) {ifhydvar[8] = -1;}
  else if (strcmp(varname,"tyy")==0) {ifhydvar[9] = -1;}
  else if (strcmp(varname,"tzz")==0) {ifhydvar[10] = -1;}
  else if (strcmp(varname,"tyz")==0) {ifhydvar[11] = -1;}
  else if (strcmp(varname,"tzx")==0) {ifhydvar[12] = -1;}
  else if (strcmp(varname,"txy")==0) {ifhydvar[13] = -1;}
  /* particle variables, 8 */
  else if (strcmp(varname,"pardx")==0) {ifparvar[0] = -1;}
  else if (strcmp(varname,"pardy")==0) {ifparvar[1] = -1;}
  else if (strcmp(varname,"pardz")==0) {ifparvar[2] = -1;}
  else if (strcmp(varname,"parcell")==0) {ifparvar[3] = -1;}
  else if (strcmp(varname,"parux")==0) {ifparvar[4] = -1;}
  else if (strcmp(varname,"paruy")==0) {ifparvar[5] = -1;}
  else if (strcmp(varname,"paruz")==0) {ifparvar[6] = -1;}
  else if (strcmp(varname,"parq")==0) {ifparvar[7] = -1;}
  else {
    if (mp_rank(grid->mp)==0) {
      printf(" *** error - unknown variable name, %s\n",varname);
      var_ls( );   /* print help list of variable names */
    }
  }
 
}
 
/**********************************************************************/
 
void vpic_simulation::var_ls( )
 
/* print help list of available variable names (both field and hydro) */
 
{
 
  // begin
 
  if (mp_rank(grid->mp) != 0) return;
 
  printf("Legal variable names are:\n");
  printf("Field variables:\n");
  printf(" ex        ey        ez        div_e_err \n");
  printf(" cbx       cby       cbz       div_b_err \n");
  printf(" tcax      tcay      tcaz      rhob      \n");
  printf(" jfx       jfy       jfz       rhof      \n");
  printf("\n");
  printf("Hydro variables:\n");
  printf(" jx        jy        jz        rho       \n");
  printf(" px        py        pz        ke        \n");
  printf(" txx       tyy       tzz       tyz       \n");
  printf(" tzx       txy                           \n");
  printf("\n");
//   particle variables - not yet set up to turn off/on
// printf("Particle variables:\n");
// printf(" pardx     pardy     pardx     parcell   \n");
// printf(" parux     paruy     paruz     parq      \n");
// printf("\n");
 
}
 
/**********************************************************************/
 
void vpic_simulation::fname_append (
      char *fname, INT *step, INT *stepdigit, INT *nrnk, INT *rankdigit)
 
  /* fname - string, base of file name to be appended to */
  /* step - integer to be appended after first dot */
  /* stepdigit - integer number of digits to print for step */
  /* nrnk - integer to be appended after second dot */
  /* rankdigit - integer number of digits for nrnk */
 
/* if either stepdigit or rankdigit is a negative value, the */
/* digit string for that value is ignored (not appended).    */
 
/* construct a filename by appending two digit strings of specified   */
/* length, separated by dots; null-terminate and return the filename. */
 
{
  INT len, sdig, rdig, nc, ii, idiv, itmp;
  char chdig[11]="0123456789";
  char csuff[NDIGSUFF];
 
  // begin
 
  sdig = *stepdigit;
  rdig = *rankdigit;
  if (sdig == 0) sdig = 6;
  if (rdig == 0) rdig = 4;
  if (sdig < 0) sdig = 0;
  if (rdig < 0) rdig = 0;
 
  len = strlen(fname);
 
  /* append a dot and the cycle number to file name */
  if (sdig > 0) {
    for (ii=0; ii<NDIGSUFF; ii++) csuff[ii] = '\0';
    strcpy(&fname[len],".");
    len = len+1;
    for (nc=0; nc<sdig; nc++) {
      idiv = 1;
      for (ii=sdig-nc-1; ii>=1; ii--) idiv = 10*idiv;
      itmp = *step/idiv;
      itmp = itmp%10;
      csuff[nc] = chdig[itmp];
    }
    strcpy(&fname[len],&csuff[0]);
 
    len = len+sdig;
  }
 
  /* append a dot and the block number to file name */
  if (rdig > 0) {
    for (ii=0; ii<NDIGSUFF; ii++) csuff[ii] = '\0';
    strcpy(&fname[len],".");
    len = len+1;
    for (nc=0; nc<rdig; nc++) {
      idiv = 1;
      for (ii=rdig-nc-1; ii>=1; ii--) idiv = 10*idiv;
      itmp = *nrnk/idiv;
      itmp = itmp%10;
      csuff[nc] = chdig[itmp];
    }
    strcpy(&fname[len],&csuff[0]);
 
    len = len+rdig;
  }
  /* always null terminate the file name */
  fname[len] = '\0';
}
 
/**********************************************************************/
 
void vpic_simulation::nultrm (
      char *cstr)
 
  /* cstr - character string to be null terminated */
 
/* null terminate the string cstr */
 
{
  INT len;
 
  // begin
 
  len = strlen(cstr);
  strcpy(&cstr[len],"\0");
}
 
/**********************************************************************/
 
void vpic_simulation::nulset (
      char *cstr)
 
  /* cstr - character string to be initialize to null */
 
/* initialize the string cstr to null */
 
{ strcpy(&cstr[0],""); }
 
/**********************************************************************/
 
void vpic_simulation::fwr_ara (
      REAL *ara, INT *nara, FILETYPE fdunit, INT *ier)
 
  /* ara - real array to be written */
  /* nara - integer number of words to be written */
  /* fdunit - file unit descriptor, must be assigned */
  /* ier - integer error sentinel, 0 for success */
 
/* write a contiguous array using a buffer */
 
{
  INT nbat, nb, naralo, naraup, na, ndxbuf, nwr, nwritten;
 
  // begin
 
  *ier = 0;
 
  if (fdunit == NULL) {
    *ier++;
    printf(" *** error - file unit for output file not assigned.\n");
  }
  /* compute number of buffer batches needed, copy and output data */
  nbat = (*nara-1)/NBUFMX + 1;
  for (nb=1; nb<=nbat; nb++) {
    naralo = (nb-1)*NBUFMX;
    naraup = naralo + NBUFMX - 1;
    if (naraup > *nara-1) naraup = *nara-1;
    ndxbuf = -1;   /* allow for c zero indexing */
    for (na=naralo; na<=naraup; na++) {
      ndxbuf++;
      buf[ndxbuf] = ara[na];
    }
    nwr = nbytwrd*(ndxbuf+1);
    nwritten = fwrite(&buf[0], 1, nwr, fdunit);
    if (nwritten != nwr) {
      *ier++;
      printf(" *** error - fwrite to output file failed in fwr_ara.\n");
    }
  }
}
 
/**********************************************************************/
 
void vpic_simulation::fwr_stride (
      REAL *ara, INT *nvartot, INT *ifvar, INT *offset,
      INT *iara, INT *jara, INT *kara, INT *istr, INT *jstr, INT *kstr,
      INT *iglo, INT *igup, INT *jglo, INT *jgup, INT *kglo, INT *kgup,
      INT *stropt, FILETYPE fdunit, INT *ier)
 
  /* ara - real array to be processed and written */
  /* nvartot - integer nvartot number of variables to write */
  /* ifvar - integer array, flags to dump or not dump variable */
  /* offset - integer size in words of the ara struct that is passed */
  /*          to this function, essential to step through variables  */
  /* iara, jara, kara - integer dimensions of array */
  /* istr, jstr, kstr - integer strides in each dimension, */
  /* iglo, igup, jglo, jgup, kglo, kgup - integer number of ghost */
  /*                                      cells at each boundary  */
  /* stropt - integer option, 1 first index, 2 center index, 3 avg */
  /* fdunit - file unit descriptor, must be assigned */
  /* ier - integer error sentinel, 0 for success */
 
/* write a 1d, 2d, or 3d strided array using a buffer */
 
/* options are either first index stropt=1, centered   */
/* index stropt=2, or average over strides stropt=3.   */
/* other option may be added later, possibly including */
/* filtered or cosine averages. this module now works  */
/* only for strides that are evenly divisible into the */
/* respective dimensions. computing for non-divisible  */
/* strides or non-orthogonal data sets will involve    */
/* some complications.                                 */
 
{
  INT nnvar, iioff;
  INT ii, jj, kk, nv;
  INT iiara, jjara, kkara, iistr, jjstr, kkstr, iara_jara, ijk, ndxbuf;
  INT iav, jav, kav, idisp, jdisp, kdisp, iimod, jjmod, kkmod;
  INT iiglo, iigup, jjglo, jjgup, kkglo, kkgup;
  INT iigupmx, jjgupmx, kkgupmx;
  REAL vsum, rcipavg_fac;
int nwordswrit=0;
 
  // begin
 
  *ier = 0;
 
  if (fdunit == NULL) {
    *ier++;
    printf(" *** error - file unit for output file not assigned.\n");
  }
  /* test for illegal option for data output, do not write any data */
  if (*stropt < 0 || *stropt >3) {
    *ier++;
    printf(" *** error - illegal option for writing strided output.\n");
  }
 
  nnvar = *nvartot;
  if (nnvar <= 0) {
    *ier++;
    printf(" *** error - zero variables for dump output.\n");
  }
  iioff = *offset;
  iiara = *iara;
  jjara = *jara;
  kkara = *kara;
  iistr = *istr;
  jjstr = *jstr;
  kkstr = *kstr;
  if (iistr<1) iistr=1;
  if (jjstr<1) jjstr=1;
  if (kkstr<1) kkstr=1;
  iiglo = *iglo;
  iigup = *igup;
  jjglo = *jglo;
  jjgup = *jgup;
  kkglo = *kglo;
  kkgup = *kgup;
  iigupmx = iiara-iigup;
  jjgupmx = jjara-jjgup;
  kkgupmx = kkara-kkgup;
 
  if (((iiara-iiglo-iigup)/iistr)*iistr != (iiara-iiglo-iigup)) {
    *ier++;
    printf(" *** error - i stride not divisible into i dimension.\n");
  }
  if (((jjara-jjglo-jjgup)/jjstr)*jjstr != (jjara-jjglo-jjgup)) {
    *ier++;
    printf(" *** error - j stride not divisible into j dimension.\n");
  }
  if (((kkara-kkglo-kkgup)/kkstr)*kkstr != (kkara-kkglo-kkgup)) {
    *ier++;
    printf(" *** error - k stride not divisible into k dimension.\n");
  }
 
  iara_jara = iiara*jjara;   /* imax times jmax spacing to k level */
 
/* patch to shoehorn in even values into the dump arrays for debug */
// for (kk=1; kk<=kkara; kk++) {
//   for (jj=1; jj<=jjara; jj++) {
//     for (ii=1; ii<=iiara; ii++) {
//       ijk = iioff * ((ii) + (jj-1)*iiara +
//             (kk-1)*iara_jara -1);
//       for (nv=0; nv<nnvar; nv++) {
//         if (ifvar[nv] >= 0) {
//           ara[ijk+nv]=(float)(ii+jj+kk);
//         }
//       }
//     }
//   }
// }
 
  /* first index data point stropt=1, or centered index stropt=2 */
  if (*stropt <= 2) {
    if (*stropt == 2) {
      idisp = iistr/2;
      jdisp = jjstr/2;
      kdisp = kkstr/2;
    }
    else {
      idisp = 0;
      jdisp = 0;
      kdisp = 0;
    }
    ndxbuf = -1;   /* allow for c zero indexing */
    for (kk=1; kk<=kkara; kk++) {
      kkmod = (kk-kkglo-1)%kkstr;
      if (kk!=1 && kk!=kkara && kkmod!=0) continue;
      for (jj=1; jj<=jjara; jj++) {
        jjmod = (jj-jjglo-1)%jjstr;
        if (jj!=1 && jj!=jjara && jjmod!=0) continue;
        for (ii=1; ii<=iiara; ii++) {
          iimod = (ii-iiglo-1)%iistr;
          if (ii!=1 && ii!=iiara && iimod!=0) continue;
          if ( ( (ii<=iiglo || ii>iigupmx) && (jj<=jjglo || jj>jjgupmx || jjmod==0)
            && (kk<=kkglo || kk>kkgupmx || kkmod==0) )
            || ( (jj<=jjglo || jj>jjgupmx) && (ii<=iiglo || ii>iigupmx || iimod==0)
            && (kk<=kkglo || kk>kkgupmx || kkmod==0) )
            || ( (kk<=kkglo || kk>kkgupmx) && (ii<=iiglo || ii>iigupmx || iimod==0)
            && (jj<=jjglo || jj>jjgupmx || jjmod==0) ) ) {
            ijk = iioff * (ii + (jj-1)*iiara + (kk-1)*iara_jara - 1);
            for (nv=0; nv<nnvar; nv++) {
              if (ifvar[nv] >= 0) {
                ndxbuf++; nwordswrit++;
                buf[ndxbuf] = ara[ijk+nv];
                if (ndxbuf >= NBUFMX-1) {
                  fwr_buf(&ndxbuf, fdunit, ier);
                }
              }
            }
          }
          else if(iimod==0 && jjmod==0 && kkmod==0) {
            ijk = iioff * ((ii+idisp) + (jj+jdisp-1)*iiara +
                  (kk+kdisp-1)*iara_jara -1);
            for (nv=0; nv<nnvar; nv++) {
              if (ifvar[nv] >= 0) {
                ndxbuf++; nwordswrit++;
                buf[ndxbuf] = ara[ijk+nv];
                if (ndxbuf >= NBUFMX-1) {
                  fwr_buf(&ndxbuf, fdunit, ier);
                }
              }
            }
          }
        }
      }
    }
    /* write any remaining unwritten data to the file */
    if (ndxbuf >= 0) {
      fwr_buf(&ndxbuf, fdunit, ier);
    }
  }
 
  /* averaged data output stropt=3 */
  else if (*stropt == 3) {
    rcipavg_fac = (REAL)(iistr*jjstr*kkstr);
    rcipavg_fac = 1.0/rcipavg_fac;
    ndxbuf = -1;   /* allow for c zero indexing */
    for (kk=1; kk<=kkara; kk++) {
      kkmod = (kk-kkglo-1)%kkstr;
      if (kk!=1 && kk!=kkara && kkmod!=0) continue;
      for (jj=1; jj<=jjara; jj++) {
        jjmod = (jj-jjglo-1)%jjstr;
        if (jj!=1 && jj!=jjara && jjmod!=0) continue;
        for (ii=1; ii<=iiara; ii++) {
          iimod = (ii-iiglo-1)%iistr;
          if (ii!=1 && ii!=iiara && iimod!=0) continue;
          if ( ( (ii<=iiglo || ii>iigupmx) && (jj<=jjglo || jj>jjgupmx || jjmod==0)
            && (kk<=kkglo || kk>kkgupmx || kkmod==0) )
            || ( (jj<=jjglo || jj>jjgupmx) && (ii<=iiglo || ii>iigupmx || iimod==0)
            && (kk<=kkglo || kk>kkgupmx || kkmod==0) )
            || ( (kk<=kkglo || kk>kkgupmx) && (ii<=iiglo || ii>iigupmx || iimod==0)
            && (jj<=jjglo || jj>jjgupmx || jjmod==0) ) ) {
            ijk = iioff * (ii + (jj-1)*iiara + (kk-1)*iara_jara - 1);
            for (nv=0; nv<nnvar; nv++) {
              if (ifvar[nv] >= 0) {
                ndxbuf++;
                buf[ndxbuf] = ara[ijk+nv];
                if (ndxbuf >= NBUFMX-1) {
                  fwr_buf(&ndxbuf, fdunit, ier);
                }
              }
            }
          }
          else if(iimod==0 && jjmod==0 && kkmod==0) {
            /* sum over all included cells, compute average */
            for (nv=0; nv<nnvar; nv++) {
              if (ifvar[nv] >= 0) {
                vsum = 0.0;
                for (kav=kk; kav<kk+kkstr; kav++) {
                  for (jav=jj; jav<jj+jjstr; jav++) {
                    for (iav=ii; iav<ii+iistr; iav++) {
                      /* calculate the index into buf array, note the -1 */
                      ijk = iioff * (iav + (jav-1)*iiara + (kav-1)*iara_jara -1);
                      vsum = vsum + ara[ijk+nv];
                    }
                  }
                }
                /* multiply sum by average factor, put average in buffer */
                ndxbuf++;
                buf[ndxbuf] = vsum * rcipavg_fac;
                if (ndxbuf >= NBUFMX-1) {
                  fwr_buf(&ndxbuf, fdunit, ier);
                }
              }
            }
          }
        }
      }
    }
    /* write any remaining unwritten data to the file */
    if (ndxbuf >= 0) {
      fwr_buf(&ndxbuf, fdunit, ier);
    }
  }
 
/* printf(" WORDS WRITTEN TO STRIDED DUMP %d\n",nwordswrit); */
}
 
/**********************************************************************/
 
void vpic_simulation::fwr_stride_block (
      REAL *ara, INT *nvartot, INT *ifvar, INT *offset,
      INT *iara, INT *jara, INT *kara, INT *istr, INT *jstr, INT *kstr,
      INT *iglo, INT *igup, INT *jglo, INT *jgup, INT *kglo, INT *kgup,
      INT *stropt, FILETYPE fdunit, INT *ier)
 
  /* ara - real array to be processed and written */
  /* nvartot - integer nvartot number of variables to write */
  /* ifvar - integer array, flags to dump or not dump variable */
  /* offset - integer size in words of the ara struct that is passed */
  /*          to this function, essential to step through variables  */
  /* iara, jara, kara - integer dimensions of array */
  /* istr, jstr, kstr - integer strides in each dimension, */
  /* iglo, igup, jglo, jgup, kglo, kgup - integer number of ghost */
  /*                                      cells at each boundary  */
  /* stropt - integer option, 1 first index, 2 center index, 3 avg */
  /* fdunit - file unit descriptor, must be assigned */
  /* ier - integer error sentinel, 0 for success */
 
/* write a 1d, 2d, or 3d strided array using a buffer */
 
/* this routine writes data in block format with all  */
/* values for a given variable together, rather than  */
/* the standard point format with all variables for a */
/* cell together. block format is needed for some of  */
/* the analysis tools. the routine fwr_stride writes  */
/* data in the normal point format (see above).       */
 
/* options are either first index stropt=1, centered   */
/* index stropt=2, or average over strides stropt=3.   */
/* other option may be added later, possibly including */
/* filtered or cosine averages. this module now works  */
/* only for strides that are evenly divisible into the */
/* respective dimensions. computing for non-divisible  */
/* strides or non-orthogonal data sets will involve    */
/* some complications.                                 */
 
{
  INT nnvar, iioff;
  INT ii, jj, kk, nv;
  INT iiara, jjara, kkara, iistr, jjstr, kkstr, iara_jara, ijk, ndxbuf;
  INT iav, jav, kav, idisp, jdisp, kdisp, iimod, jjmod, kkmod;
  INT iiglo, iigup, jjglo, jjgup, kkglo, kkgup;
  INT iigupmx, jjgupmx, kkgupmx;
  REAL vsum, rcipavg_fac;
 
  // begin
 
  *ier = 0;
 
  if (fdunit == NULL) {
    *ier++;
    printf(" *** error - file unit for output file not assigned.\n");
  }
  /* test for illegal option for data output, do not write any data */
  if (*stropt < 0 || *stropt >3) {
    *ier++;
    printf(" *** error - illegal option for writing strided output.\n");
  }
 
  nnvar = *nvartot;
  if (nnvar <= 0) {
    *ier++;
    printf(" *** error - zero variables for dump output.\n");
  }
  iioff = *offset;
  iiara = *iara;
  jjara = *jara;
  kkara = *kara;
  iistr = *istr;
  jjstr = *jstr;
  kkstr = *kstr;
  if (iistr<1) iistr=1;
  if (jjstr<1) jjstr=1;
  if (kkstr<1) kkstr=1;
  iiglo = *iglo;
  iigup = *igup;
  jjglo = *jglo;
  jjgup = *jgup;
  kkglo = *kglo;
  kkgup = *kgup;
  iigupmx = iiara-iigup;
  jjgupmx = jjara-jjgup;
  kkgupmx = kkara-kkgup;
 
  if (((iiara-iiglo-iigup)/iistr)*iistr != (iiara-iiglo-iigup)) {
    *ier++;
    printf(" *** error - i stride not divisible into i dimension.\n");
  }
  if (((jjara-jjglo-jjgup)/jjstr)*jjstr != (jjara-jjglo-jjgup)) {
    *ier++;
    printf(" *** error - j stride not divisible into j dimension.\n");
  }
  if (((kkara-kkglo-kkgup)/kkstr)*kkstr != (kkara-kkglo-kkgup)) {
    *ier++;
    printf(" *** error - k stride not divisible into k dimension.\n");
  }
 
  iara_jara = iiara*jjara;   /* imax times jmax spacing to k level */
 
  /* first index data point stropt=1, or centered index stropt=2 */
  if (*stropt <= 2) {
    if (*stropt == 2) {
      idisp = iistr/2;
      jdisp = jjstr/2;
      kdisp = kkstr/2;
    }
    else {
      idisp = 0;
      jdisp = 0;
      kdisp = 0;
    }
    ndxbuf = -1;   /* allow for c zero indexing */
    for (nv=0; nv<nnvar; nv++) {
      if (ifvar[nv] >= 0) {
        for (kk=1; kk<=kkara; kk++) {
          kkmod = (kk-kkglo-1)%kkstr;
          if (kk!=1 && kk!=kkara && kkmod!=0) continue;
          for (jj=1; jj<=jjara; jj++) {
            jjmod = (jj-jjglo-1)%jjstr;
            if (jj!=1 && jj!=jjara && jjmod!=0) continue;
            for (ii=1; ii<=iiara; ii++) {
              iimod = (ii-iiglo-1)%iistr;
              if (ii!=1 && ii!=iiara && iimod!=0) continue;
              if ( ( (ii<=iiglo || ii>iigupmx) && (jj<=jjglo || jj>jjgupmx || jjmod==0)
                && (kk<=kkglo || kk>kkgupmx || kkmod==0) )
                || ( (jj<=jjglo || jj>jjgupmx) && (ii<=iiglo || ii>iigupmx || iimod==0)
                && (kk<=kkglo || kk>kkgupmx || kkmod==0) )
                || ( (kk<=kkglo || kk>kkgupmx) && (ii<=iiglo || ii>iigupmx || iimod==0)
                && (jj<=jjglo || jj>jjgupmx || jjmod==0) ) ) {
                ijk = iioff * (ii + (jj-1)*iiara + (kk-1)*iara_jara - 1);
                ndxbuf++;
                buf[ndxbuf] = ara[ijk+nv];
                if (ndxbuf >= NBUFMX-1) {
                  fwr_buf(&ndxbuf, fdunit, ier);
                }
              }
              else if(iimod==0 && jjmod==0 && kkmod==0) {
                ijk = iioff * ((ii+idisp) + (jj+jdisp-1)*iiara +
                      (kk+kdisp-1)*iara_jara -1);
                ndxbuf++;
                buf[ndxbuf] = ara[ijk+nv];
                if (ndxbuf >= NBUFMX-1) {
                  fwr_buf(&ndxbuf, fdunit, ier);
                }
              }
            }
          }
        }
        /* write any remaining unwritten data to the file */
        if (ndxbuf >= 0) {
          fwr_buf(&ndxbuf, fdunit, ier);
        }
      }
    }
  }
 
  /* averaged data output stropt=3 */
  else if (*stropt == 3) {
    rcipavg_fac = (REAL)(iistr*jjstr*kkstr);
    rcipavg_fac = 1.0/rcipavg_fac;
    ndxbuf = -1;   /* allow for c zero indexing */
    for (nv=0; nv<nnvar; nv++) {
      if (ifvar[nv] >= 0) {
        for (kk=1; kk<=kkara; kk++) {
          kkmod = (kk-kkglo-1)%kkstr;
          if (kk!=1 && kk!=kkara && kkmod!=0) continue;
          for (jj=1; jj<=jjara; jj++) {
            jjmod = (jj-jjglo-1)%jjstr;
            if (jj!=1 && jj!=jjara && jjmod!=0) continue;
            for (ii=1; ii<=iiara; ii++) {
              iimod = (ii-iiglo-1)%iistr;
              if (ii!=1 && ii!=iiara && iimod!=0) continue;
              if ( ( (ii<=iiglo || ii>iigupmx) && (jj<=jjglo || jj>jjgupmx || jjmod==0)
                && (kk<=kkglo || kk>kkgupmx || kkmod==0) )
                || ( (jj<=jjglo || jj>jjgupmx) && (ii<=iiglo || ii>iigupmx || iimod==0)
                && (kk<=kkglo || kk>kkgupmx || kkmod==0) )
                || ( (kk<=kkglo || kk>kkgupmx) && (ii<=iiglo || ii>iigupmx || iimod==0)
                && (jj<=jjglo || jj>jjgupmx || jjmod==0) ) ) {
                ijk = iioff * (ii + (jj-1)*iiara + (kk-1)*iara_jara - 1);
                ndxbuf++;
                buf[ndxbuf] = ara[ijk+nv];
                 if (ndxbuf >= NBUFMX-1) {
                  fwr_buf(&ndxbuf, fdunit, ier);
                }
              }
              else if(iimod==0 && jjmod==0 && kkmod==0) {
                /* sum over all included cells, compute average */
                vsum = 0.0;
                for (kav=kk; kav<kk+kkstr; kav++) {
                  for (jav=jj; jav<jj+jjstr; jav++) {
                    for (iav=ii; iav<ii+iistr; iav++) {
                      /* calculate the index into buf array, note the -1 */
                      ijk = iioff * (iav + (jav-1)*iiara + (kav-1)*iara_jara -1);
                      vsum = vsum + ara[ijk+nv];
                    }
                  }
                }
                /* multiply sum by average factor, put average in buffer */
                ndxbuf++;
                buf[ndxbuf] = vsum * rcipavg_fac;
                if (ndxbuf >= NBUFMX-1) {
                  fwr_buf(&ndxbuf, fdunit, ier);
                }
              }
            }
          }
        }
      }
    }
    /* write any remaining unwritten data to the file */
    if (ndxbuf >= 0) {
      fwr_buf(&ndxbuf, fdunit, ier);
    }
  }
 
}
 
/**********************************************************************/
 
void vpic_simulation::fwr_stride_grid (
      REAL *ara, INT *nvartot, INT *ifxyz, INT *nvstrt, INT *offset,
      INT *iara, INT *jara, INT *kara, INT *istr, INT *jstr, INT *kstr,
      INT *iglo, INT *igup, INT *jglo, INT *jgup, INT *kglo, INT *kgup,
      REAL *x0, REAL *y0, REAL *z0, REAL *dx, REAL *dy, REAL *dz,
      INT *stropt, FILETYPE fdunit, INT *ier)
 
  /* ara - array to be processed and written */
  /* nvartot - integer nvartot number of variables to write */
  /* ifxyz - integer option, 1 to compute and write x,y,z coordinates */
  /*         at the beginning of the grid data (first three variables)*/
  /* nvstrt - integer starting variable position of the grid data */
  /*          within the field struct; currently, 17th variable   */
  /* offset - integer size in words of the ara struct that is passed */
  /*          to this function, essential to step through variables  */
  /* iara, jara, kara - integer dimensions of array */
  /* istr, jstr, kstr - integer strides in each dimension, */
  /* iglo, igup, jglo, jgup, kglo, kgup - integer number of ghost */
  /*                                      cells at each boundary  */
  /* x0, y0, z0 - lower x,y,z coordinates of this domain */
  /* dx, dz, dz - delta x, delta y, delta z this domain */
  /*              (note - these are unstrided deltas)   */
  /* stropt - integer option, 1 first index, 2 center index */
  /* fdunit - file unit descriptor, must be assigned */
  /* ier - integer error sentinel, 0 for success */
 
/* write a 1d, 2d, or 3d strided grid array using a buffer */
 
/* the striding options stropt are: 1, first index;  */
/* or 2, centered index. no other options available. */
 
{
  INT nnvar, iioff;
  INT ii, jj, kk, nv, nnvstart;
  INT iiara, jjara, kkara, iistr, jjstr, kkstr, iara_jara, ijk, ndxbuf;
  INT idisp, jdisp, kdisp, iimod, jjmod, kkmod, iipos, jjpos, kkpos;
  INT iiglo, iigup, jjglo, jjgup, kkglo, kkgup;
  INT iigupmx, jjgupmx, kkgupmx;
 
  // begin
 
  *ier = 0;
 
  if (fdunit == NULL) {
    *ier++;
    printf(" *** error - file unit for output file not assigned.\n");
  }
  /* test for illegal option for data output, do not write any data */
  if (*stropt < 0 || *stropt > 2) {
    *ier++;
    printf(" *** error - illegal option for writing strided output.\n");
  }
 
  nnvar = *nvartot;
  if (nnvar <= 0) {
    *ier++;
    printf(" *** error - zero variables for dump output.\n");
  }
  nnvstart = *nvstrt;
  iioff = *offset;
  iiara = *iara;
  jjara = *jara;
  kkara = *kara;
  iistr = *istr;
  jjstr = *jstr;
  kkstr = *kstr;
  if (iistr<1) iistr=1;
  if (jjstr<1) jjstr=1;
  if (kkstr<1) kkstr=1;
  iiglo = *iglo;
  iigup = *igup;
  jjglo = *jglo;
  jjgup = *jgup;
  kkglo = *kglo;
  kkgup = *kgup;
  iigupmx = iiara-iigup;
  jjgupmx = jjara-jjgup;
  kkgupmx = kkara-kkgup;
 
  if (((iiara-iiglo-iigup)/iistr)*iistr != (iiara-iiglo-iigup)) {
    *ier++;
    printf(" *** error - i stride not divisible into i dimension.\n");
  }
  if (((jjara-jjglo-jjgup)/jjstr)*jjstr != (jjara-jjglo-jjgup)) {
    *ier++;
    printf(" *** error - j stride not divisible into j dimension.\n");
  }
  if (((kkara-kkglo-kkgup)/kkstr)*kkstr != (kkara-kkglo-kkgup)) {
    *ier++;
    printf(" *** error - k stride not divisible into k dimension.\n");
  }
 
  iara_jara = iiara*jjara;   /* imax times jmax spacing to k level */
 
  /* first index data point stropt=1, or centered index stropt=2 */
  if (*stropt <= 1) {
    idisp = 0;
    jdisp = 0;
    kdisp = 0;
  }
  else {
    idisp = iistr/2;
    jdisp = jjstr/2;
    kdisp = kkstr/2;
  }
 
  ndxbuf = -1;   /* allow for c zero indexing */
 
  for (kk=1; kk<=kkara; kk++) {
    kkmod = (kk-kkglo-1)%kkstr;
    if (kk!=1 && kk!=kkara && kkmod!=0) continue;
    for (jj=1; jj<=jjara; jj++) {
      jjmod = (jj-jjglo-1)%jjstr;
      if (jj!=1 && jj!=jjara && jjmod!=0) continue;
      for (ii=1; ii<=iiara; ii++) {
        iimod = (ii-iiglo-1)%iistr;
        if (ii!=1 && ii!=iiara && iimod!=0) continue;
        if ( ( (ii<=iiglo || ii>iigupmx) && (jj<=jjglo || jj>jjgupmx || jjmod==0)
          && (kk<=kkglo || kk>kkgupmx || kkmod==0) )
          || ( (jj<=jjglo || jj>jjgupmx) && (ii<=iiglo || ii>iigupmx || iimod==0)
          && (kk<=kkglo || kk>kkgupmx || kkmod==0) )
          || ( (kk<=kkglo || kk>kkgupmx) && (ii<=iiglo || ii>iigupmx || iimod==0)
          && (jj<=jjglo || jj>jjgupmx || jjmod==0) ) ) {
          ijk = iioff * (ii + (jj-1)*iiara + (kk-1)*iara_jara - 1);
          if (*ifxyz == 1) {
            iipos = MAX0((ii-iiglo-1),0);
            iipos = MIN0(iipos,iigupmx-2);
            jjpos = MAX0((jj-jjglo-1),0);
            jjpos = MIN0(jjpos,jjgupmx-2);
            kkpos = MAX0((kk-kkglo-1),0);
            kkpos = MIN0(kkpos,kkgupmx-2);
            if (ndxbuf >= NBUFMX-4) {
              fwr_buf(&ndxbuf, fdunit, ier);
            }
            ndxbuf++;
            buf[ndxbuf] = *x0 + iipos*(*dx);
            ndxbuf++;
            buf[ndxbuf] = *y0 + jjpos*(*dy);
            ndxbuf++;
            buf[ndxbuf] = *z0 + kkpos*(*dz);
          }
          for (nv=0; nv<nnvar; nv++) {
            ndxbuf++;
            memcpy(&buf[ndxbuf],&ara[ijk+nv+nnvstart],4);
            if (ndxbuf >= NBUFMX-1) {
              fwr_buf(&ndxbuf, fdunit, ier);
            }
          }
        }
        else if(iimod==0 && jjmod==0 && kkmod==0) {
          ijk = iioff * ((ii+idisp) + (jj+jdisp-1)*iiara +
                (kk+kdisp-1)*iara_jara -1);
          if (*ifxyz == 1) {
            iipos = MAX0((ii-iiglo-1),0);
            iipos = MIN0(iipos,iigupmx-2);
            jjpos = MAX0((jj-jjglo-1),0);
            jjpos = MIN0(jjpos,jjgupmx-2);
            kkpos = MAX0((kk-kkglo-1),0);
            kkpos = MIN0(kkpos,kkgupmx-2);
            if (ndxbuf >= NBUFMX-4) {
              fwr_buf(&ndxbuf, fdunit, ier);
            }
            ndxbuf++;
            buf[ndxbuf] = *x0 + iipos*(*dx);
            ndxbuf++;
            buf[ndxbuf] = *y0 + jjpos*(*dy);
            ndxbuf++;
            buf[ndxbuf] = *z0 + kkpos*(*dz);
          }
          for (nv=0; nv<nnvar; nv++) {
            ndxbuf++;
            memcpy(&buf[ndxbuf],&ara[ijk+nv+nnvstart],4);
            if (ndxbuf >= NBUFMX-1) {
              fwr_buf(&ndxbuf, fdunit, ier);
            }
          }
        }
      }
    }
  }
  /* write any remaining unwritten data to the file */
  if (ndxbuf >= 0) {
    fwr_buf(&ndxbuf, fdunit, ier);
  }
 
}
 
/**********************************************************************/
 
void vpic_simulation::fwr_par (
      REAL *ara, INT *nara, INT *nvartot, INT *ifvar, INT *nstr,
      FILETYPE fdunit, INT *ier)
 
  /* ara - real array to be processed and written */
  /* nara - integer length of array */
  /* nvartot - integer number of variables in the array */
  /* ifvar - integer array of flags to select variables to write, */
  /*          if 0 or positive, write; if negative, do not write  */
  /* nstr - integer stride to stride through array and write */
  /* fdunit - file unit descriptor, must be assigned */
  /* ier - integer error sentinel, 0 for success */
 
/* write a strided data set that is more than one variable in length  */
/* to output file. this function is intended for writing a particle   */
/* data structure to the output file but should work for any data set */
/* that consists of more than one variable. note that a patch had to  */
/* be made below because using ndxbuf=-1 as a sentinel for the end of */
/* data will not work with the strided data that has more than one    */
/* variable in each group (such as the particle data structure).      */
 
{
  INT nnara, nnvartot, nv, nnstr, ndxbuf, ivndx;
  INT ii;
 
  // begin
 
  *ier = 0;
 
  if (fdunit == NULL) {
    *ier++;
    printf(" *** error - file unit for output file not assigned.\n");
  }
 
  nnara = *nara;
  nnvartot = *nvartot;
  nnstr = *nstr;
  if (nnvartot <= 0) {
    *ier++;
    printf(" *** error - zero variables for particle dump output.\n");
  }
  if (nnstr <= 0) {
    *ier++;
    printf(" *** error - stride of zero for particle dump output.\n");
  }
  ndxbuf = -1;   /* allow for c zero indexing */
  for (ii=1; ii<=nnara; ii=ii+nnstr) {
    ivndx = (ii-1)*nnvartot;
    for (nv=0; nv<nnvartot; nv++) {
      if (ifvar[nv] >= 0) {
        ndxbuf++;
        memcpy(&buf[ndxbuf], &ara[ivndx], nbytwrd);
        if (ndxbuf >= NBUFMX-1) {
          fwr_buf(&ndxbuf, fdunit, ier);
        }
      }
      ivndx++;
    }
  }
  /* write any remaining unwritten data to the file */
  if (ndxbuf >= 0) {
    fwr_buf(&ndxbuf, fdunit, ier);
  }
}
 
/**********************************************************************/
 
void vpic_simulation::fwr_buf (
      INT *ndxbuf, FILETYPE fdunit, INT *ier)
 
  /* ndxbuf - integer number of words to write from buf */
  /* fdunit - file unit descriptor, must be open */
  /* ier - integer error sentinel, 0 for success */
 
/* write ndxbuf words in the buffer buf[] to output file fdunit */
 
/* ndxbuf is reset to -1 and returned to allow for c 0 indexing */
 
{
  INT nwr, nwritten;
 
  // begin
 
  nwr = nbytwrd*(*ndxbuf+1);
  nwritten = fwrite(&buf[0], 1, nwr, fdunit);
  if (nwritten != nwr) {
    *ier++;
    printf(" *** error - fwrite to output file failed in fwr_buf.\n");
  }
  *ndxbuf = -1;   /* allow for c zero indexing */
                  /* when resetting pointer    */
}
 
/**********************************************************************/
 
void vpic_simulation::fwr_bin_head (
      FILETYPE fdunit, INT *ier)
 
  /* fdunit - file unit descriptor to write, must be open */
  /* ier - integer error sentinel, 0 for success */
 
/* write the two-word type and length header, the 24-word integer */
/* header, and the 24-word float header at start of the dump file */
 
{
  INT nwr, nwritten;
 
  // begin
 
  *ier = 0;
 
  nwr = nbytwrd;
  nwritten = fwrite(&ityphead, 1, nwr, fdunit);
  if (nwritten != nwr) {
    *ier++;
    printf(" *** error - fwrite to output file failed in fwr_bin_head.\n");
  }
 
  idhead[0] = 101;
  idhead[1] = NWHEAD_INT;
  nwr = nbytwrd*2;
  nwritten = fwrite(&idhead[0], 1, nwr, fdunit);
  if (nwritten != nwr) {
    *ier++;
    printf(" *** error - fwrite to output file failed in fwr_bin_head.\n");
  }
  nwr = nbytwrd * NWHEAD_INT;
  nwritten = fwrite(&nrnk, 1, nwr, fdunit);
  if (nwritten != nwr) {
    *ier++;
    printf(" *** error - fwrite to output file failed in fwr_bin_head.\n");
  }
 
  idhead[0] = 103;
  idhead[1] = NWHEAD_REAL;
  nwr = nbytwrd*2;
  nwritten = fwrite(&idhead[0], 1, nwr, fdunit);
  if (nwritten != nwr) {
    *ier++;
    printf(" *** error - fwrite to output file failed in fwr_bin_head.\n");
  }
  nwr = nbytwrd * NWHEAD_REAL;
  nwritten = fwrite(&tim, 1, nwr, fdunit);
  if (nwritten != nwr) {
    *ier++;
    printf(" *** error - fwrite to output file failed in fwr_bin_head.\n");
  }
 
}
 
 
/**********************************************************************/
 
void vpic_simulation::fwr_asc_head (
      FILETYPE fdunit, INT *ier)
 
  /* fdunit - file unit descriptor to write, must be open */
  /* ier - integer error sentinel, 0 for success */
 
/* write the header to an ascii header file. this file contains   */
/* the one-word integer dump type specification, the two-word     */
/* integer type and length header, the integer data with one word */
/* per line, the two-word float type and length header, and the   */
/* float data with one word per line.                             */
 
{
  INT nint, nflt;
 
  // begin
 
  *ier = 0;
 
  fprintf(fdunit, "%d\n", ityphead);
 
  idhead[0] = 101;
  idhead[1] = NWHEAD_INT;
  fprintf(fdunit, "%d %d\n", idhead[0], idhead[1]);
  nint = NWHEAD_INT;
  prtint(fdunit, &nrnk, &nint);
 
  idhead[0] = 103;
  idhead[1] = NWHEAD_REAL;
  fprintf(fdunit, "%d %d\n", idhead[0], idhead[1]);
  nflt = NWHEAD_REAL;
  prtflt(fdunit, &tim, &nflt);
 
}
 
/**********************************************************************/
 
void vpic_simulation::prtint(FILE* fdasc, int *intara, int *nprt)
 
  /* fdasc - file descriptor of open ascii output file */
  /* fltara - int array to print to ascii file         */
  /* nprt - integer number of array members to print   */
 
/* print an integer array sequentially to an ascii file */
 
{
  int ii;
 
  for (ii=1; ii<=*nprt; ii++) {
    fprintf(fdasc, "%d\n", *intara++);
  }
}
 
/**********************************************************************/
 
void vpic_simulation::prtflt(FILE* fdasc, float *fltara, int *nprt)
 
  /* fdasc - file descriptor of open ascii output file */
  /* fltara - float array to print to ascii file       */
  /* nprt - integer number of array members to print   */
 
/* print a floating point array sequentially to an ascii file */
 
{
  int ii;
 
  for (ii=1; ii<=*nprt; ii++) {
    fprintf(fdasc, "%e\n", *fltara++);
  }
}
 
/**********************************************************************/
 
void vpic_simulation::prtstr(FILE* fdasc, char *strara, int *nprt)
 
  /* fdasc - file descriptor of open ascii output file */
  /* strara - character string array to print to file  */
  /* nprt - integer number of array members to print   */
 
/* print an array of character strings sequentially to an ascii file */
 
{
  int ii, jj, ic, slen;
 
  slen = strlen(&strara[0]) + 1;
  ic=0;
 
  for (ii=1; ii<=*nprt; ii++) {
    /* each string is printed one character at a time */
    for (jj=1; jj<=slen; jj++) {
      if (jj%slen == 0) {
        fprintf(fdasc, "\n");
      }
      else {
        fprintf(fdasc, "%c", strara[ic]);
      }
      ic++;
    }
  }
}
 
/**********************************************************************/
 
void vpic_simulation::define_glb( double xl,  double yl,  double zl,
  double xh,  double yh,  double zh,
  int gnx, int gny, int gnz, int gpx, int gpy, int gpz )
 
  /* xl - lower x boundary of entire spatial domain */
  /* yl - lower y boundary of entire spatial domain */
  /* zl - lower z boundary of entire spatial domain */
  /* xh - upper x boundary of entire spatial domain */
  /* yh - upper y boundary of entire spatial domain */
  /* zh - upper z boundary of entire spatial domain */
  /* gnx - number of x cells in entire spatial domain (no ghosts) */
  /* gny - number of y cells in entire spatial domain (no ghosts) */
  /* gnz - number of z cells in entire spatial domain (no ghosts) */
  /* gpx - number of processors requested in x dimension */
  /* gpy - number of processors requested in y dimension */
  /* gpz - number of processors requested in z dimension */
 
/* save the spatial limits, number of cells and number of processors */
/* for the entire spatial domain (not local processor). essential    */
/* for use in the strided postprocessor dumps. this function must be */
/* called in all the grid initialization procedures for all cases.   */
/* it has been included in the present geometry define routines in   */
/* vpic.hxx and must likewise be called for any corresponding future */
/* routines are added.                                               */
 
{
    xmnglb = xl;
    xmxglb = xh;
    ymnglb = yl;
    ymxglb = yh;
    zmnglb = zl;
    zmxglb = zh;
    icelglb = gnx;
    jcelglb = gny;
    kcelglb = gnz;
    itopo = gpx;
    jtopo = gpy;
    ktopo = gpz;
}
 
/**********************************************************************/
 
#ifdef BIT64
#define NBYTWRD 8
#define LONGLONGFLIP long long int
#define LONGFLIP long int
#else
#define NBYTWRD 4
#define LONGLONGFLIP long long int
#define LONGFLIP long int
#endif
 
/**********************************************************************/
 
/* vectorized byte flip functions, for arrays of 8-byte or 4-byte words */
 
void vpic_simulation::flipbyte ( LONGLONGFLIP *iara, LONGLONGFLIP *nara)
 
  /* iara - array of words to be byte flipped */
  /* nara - integer number of words to be flipped */
 
/* byte flip an array of 64-bit (8-byte) numbers, end-for-end.        */
/* this converts from compaq alpha format to standard or vice versa.  */
 
{
      LONGLONGFLIP ii, itmp;
      LONGLONGFLIP msk1=0377, msk2=0177400, msk3=077600000,
       msk4=037700000000, msk5=017740000000000,
       msk6=07760000000000000, msk7=03770000000000000000,
       msk8=01774000000000000000000;
 
      // begin
 
      for (ii=0; ii<*nara; ii++) {
        itmp = (iara[ii]>>56) & msk1;
        itmp = (itmp) | ((iara[ii]>>40) & msk2);
        itmp = (itmp) | ((iara[ii]>>24) & msk3);
        itmp = (itmp) | ((iara[ii]>>8)  & msk4);
        itmp = (itmp) | ((iara[ii]<<8)  & msk5);
        itmp = (itmp) | ((iara[ii]<<24) & msk6);
        itmp = (itmp) | ((iara[ii]<<40) & msk7);
        itmp = (itmp) | ((iara[ii]<<56) & msk8);
        iara[ii] = itmp;
      }
}
 
/**********************************************************************/
 
void vpic_simulation::flipbyte4 ( LONGFLIP *iara, LONGFLIP *nara)
 
  /* iara - array of words to be byte flipped */
  /* nara - integer number of words to be flipped */
 
/* byte flip an array of 32-bit (4-byte) numbers, end-for-end.        */
/* this converts from compaq alpha format to standard or vice versa.  */
 
{
      LONGFLIP ii, itmp;
      LONGFLIP msk1=0377, msk2=0177400, msk3=077600000, msk4=037700000000;
 
      // begin
 
      for (ii=0; ii<*nara; ii++) {
        itmp = (iara[ii]>>24) & msk1;
        itmp = (itmp) | ((iara[ii]>>8)  & msk2);
        itmp = (itmp) | ((iara[ii]<<8)  & msk3);
        itmp = (itmp) | ((iara[ii]<<24) & msk4);
        iara[ii] = itmp;
      }
}
/**********************************************************************/
 
#undef SETIVAR
#undef SETDVAR
#undef ITEST
#undef DTEST
