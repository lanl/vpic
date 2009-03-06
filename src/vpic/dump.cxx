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
#include "dumpmacros.h"
#include <FileIO.hxx>
#include <FileUtils.hxx>
#include <BitField.hxx>
#include <cassert>
 
#define VERBOSE 0

// FIXME: NEW FIELDS IN THE GRID READ/WRITE WAS HACKED UP TO BE BACKWARD
// COMPATIBLE WITH EXISTING EXTERNAL 3RD PARTY VISUALIZATION SOFTWARE.
// IN THE LONG RUN, THIS EXTERNAL SOFTWARE WILL NEED TO BE UPDATED.
 
int vpic_simulation::dump_mkdir(const char * dname) {
	return FileUtils::makeDirectory(dname);
} // dump_mkdir

int vpic_simulation::dump_cwd(char * dname, size_t size) {
	return FileUtils::getCurrentWorkingDirectory(dname, size);
} // dump_mkdir

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
    status = fileIO.open(fname, append ? io_append : io_write);
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
  const int history_dump = 5;
} // namespace
 
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
 
  /* IMPORTANT: these values are written in WRITE_HEADER_V0 */
  nxout = grid->nx;
  nyout = grid->ny;
  nzout = grid->nz;
  dxout = grid->dx;
  dyout = grid->dy;
  dzout = grid->dz;

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
 
  /* IMPORTANT: these values are written in WRITE_HEADER_V0 */
  nxout = grid->nx;
  nyout = grid->ny;
  nzout = grid->nz;
  dxout = grid->dx;
  dyout = grid->dy;
  dzout = grid->dz;

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
 
  /* IMPORTANT: these values are written in WRITE_HEADER_V0 */
  nxout = grid->nx;
  nyout = grid->ny;
  nzout = grid->nz;
  dxout = grid->dx;
  dyout = grid->dy;
  dzout = grid->dz;

  WRITE_HEADER_V0( dump_type::hydro_dump,sp->id,sp->q_m,fileIO);
 
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
 
  /* IMPORTANT: these values are written in WRITE_HEADER_V0 */
  nxout = grid->nx;
  nyout = grid->ny;
  nzout = grid->nz;
  dxout = grid->dx;
  dyout = grid->dy;
  dzout = grid->dz;

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
  FileIOUnswapped fileIO;
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
 
  /* IMPORTANT: these values are written in WRITE_HEADER_V0 */
  nxout = grid->nx;
  nyout = grid->ny;
  nzout = grid->nz;
  dxout = grid->dx;
  dyout = grid->dy;
  dzout = grid->dz;

  WRITE_HEADER_V0( dump_type::restart_dump, invalid_species_id,
  	0, fileIO );
 
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
  FileIOUnswapped fileIO;
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
  FileIOStatus status = fileIO.open( fname, io_read );
  if( status == fail) ABORT(true);
 
  // Machine compatibility information
  READ(char,     c,fileIO);  ABORT(c!=CHAR_BIT         );
  READ(char,     c,fileIO);  ABORT(c!=2                );
  READ(char,     c,fileIO);  ABORT(c!=4                );
  READ(char,     c,fileIO);  ABORT(c!=4                );
  READ(char,     c,fileIO);  ABORT(c!=8                );
  READ(short int,s,fileIO);  ABORT(s!=(short int)0xcafe);
  READ(int,      n,fileIO);  ABORT(n!=(int)0xdeadbeef  );
  READ(float,    f,fileIO);  ABORT(f!=1.0              );
  READ(double,   d,fileIO);  ABORT(d!=1.0              );
 
  // Dump type and header version
  READ(int,n,fileIO); ABORT(n!=0           );
  READ(int,n,fileIO); ABORT(n!=dump_type::restart_dump);
 
  // High level information
  READ(int,  step,      fileIO); ABORT(step<0       );
  READ(int,  grid->nx,  fileIO); ABORT(grid->nx<1   );
  READ(int,  grid->ny,  fileIO); ABORT(grid->ny<1   );
  READ(int,  grid->nz,  fileIO); ABORT(grid->nz<1   );
  READ(float,grid->dt,  fileIO); ABORT(grid->dt<=0  );
  READ(float,grid->dx,  fileIO); ABORT(grid->dx<=0  );
  READ(float,grid->dy,  fileIO); ABORT(grid->dy<=0  );
  READ(float,grid->dz,  fileIO); ABORT(grid->dz<=0  );
  READ(float,grid->x0,  fileIO);
  READ(float,grid->y0,  fileIO);
  READ(float,grid->z0,  fileIO);
  READ(float,grid->cvac,fileIO); ABORT(grid->cvac<=0);
  READ(float,grid->eps0,fileIO); ABORT(grid->eps0<=0);
  READ(float,grid->damp,fileIO); ABORT(grid->damp<0 );
  READ(int,  n,         fileIO); ABORT(n!=rank      );
  READ(int,  n,         fileIO); ABORT(n!=nproc     );
 
  // Species parameters
  READ(int,  n,fileIO); ABORT(n!=invalid_species_id);
  READ(float,f,fileIO); ABORT(f!=0                 );
 
  /************************************************
   * Restart saved directly initialized variables *
   ************************************************/
 
  // variables not already restored above
 
  READ(int,    num_step,             fileIO);
  READ(int,    status_interval,      fileIO);
  READ(int,    clean_div_e_interval, fileIO);
  READ(int,    clean_div_b_interval, fileIO);
  READ(int,    sync_shared_interval, fileIO);
  READ(double, quota,                fileIO);
  READ(int,    restart_interval,     fileIO);
  READ(int,    hydro_interval,       fileIO);
  READ(int,    field_interval,       fileIO);
  READ(int,    particle_interval,    fileIO);
 
  /**********************************************
   * Restart saved helper initialized variables *
   **********************************************/
 
  // random number generator
 
  READ(int,size,  fileIO);           ABORT(size!=sizeof(char));
  READ(int,ndim,  fileIO);           ABORT(ndim!=1           );
  READ(int,dim[0],fileIO);           ABORT(dim[0]<0          );
 
  MALLOC( state, dim[0] );
  fileIO.read(state, dim[0]);
  set_mt_rng_state( rng, state, dim[0] );
  FREE( state );
 
  // materials ... material_list must be put together in the same
  // order it was originally in
 
  READ(int,n,fileIO); ABORT(n<1);
  for(last_m=NULL;n;n--) {
    char * buf;
    // Note: sizeof(m[0]) includes terminating '\0'
    READ(int,size,fileIO);                       ABORT(size<1 );
    MALLOC( buf, sizeof(m[0])+size ); m = (material_t *)buf;
    fileIO.read(m->name, size);
    m->name[size]='\0';
 
    READ(short int,m->id,     fileIO); ABORT(m->id!=n-1);
    READ(float,    m->epsx,   fileIO);
    READ(float,    m->epsy,   fileIO);
    READ(float,    m->epsz,   fileIO);
    READ(float,    m->mux,    fileIO);
    READ(float,    m->muy,    fileIO);
    READ(float,    m->muz,    fileIO);
    READ(float,    m->sigmax, fileIO);
    READ(float,    m->sigmay, fileIO);
    READ(float,    m->sigmaz, fileIO);
    READ(float,    m->zetax,  fileIO);
    READ(float,    m->zetay,  fileIO);
    READ(float,    m->zetaz,  fileIO);
 
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
 
  READ(float,grid->rdx, fileIO);
  READ(float,grid->rdy, fileIO);
  READ(float,grid->rdz, fileIO);
  READ(float,grid->x1,  fileIO);
  READ(float,grid->y1,  fileIO);
  READ(float,grid->z1,  fileIO);
 
  READ(int,size,  fileIO); ABORT(size!=sizeof(grid->bc[0]));
  READ(int,ndim,  fileIO); ABORT(ndim!=3          );
  READ(int,dim[0],fileIO); ABORT(dim[0]!=3        );
  READ(int,dim[1],fileIO); ABORT(dim[1]!=3        );
  READ(int,dim[2],fileIO); ABORT(dim[2]!=3        );
  fileIO.read( grid->bc, dim[0]*dim[1]*dim[2]);
 
  READ(int,size,fileIO);   ABORT(size!=sizeof(grid->range[0]));
  READ(int,ndim,fileIO);   ABORT(ndim!=1          );
  READ(int,dim[0],fileIO); ABORT(dim[0]!=nproc+1  );
  MALLOC_ALIGNED( grid->range, dim[0], 16 );
  fileIO.read( grid->range, dim[0] );
 
  READ(int,size,  fileIO); ABORT(size!=sizeof(grid->neighbor[0]) );
  READ(int,ndim,  fileIO); ABORT(ndim!=4           );
  READ(int,dim[0],fileIO); ABORT(dim[0]!=6         );
  READ(int,dim[1],fileIO); ABORT(dim[1]!=grid->nx+2);
  READ(int,dim[2],fileIO); ABORT(dim[2]!=grid->ny+2);
  READ(int,dim[3],fileIO); ABORT(dim[3]!=grid->nz+2);
  MALLOC_ALIGNED( grid->neighbor, dim[0]*dim[1]*dim[2]*dim[3], 128 );
  fileIO.read( grid->neighbor, dim[0]*dim[1]*dim[2]*dim[3] );
 
  grid->rangel = grid->range[rank];
  grid->rangeh = grid->range[rank+1]-1;
 
  // fields
 
  // FIXME: SEE NOTE IN ABOVE ABOUT THIS!
  field_advance_methods_t fam[1];
  READ(field_advance_methods_t,fam[0],fileIO);
  field_advance = new_field_advance( grid, material_list, fam );
  field = field_advance->f; // FIXME: Temporary hack
 
  READ(int,size,  fileIO); ABORT(size!=sizeof(field[0]));
  READ(int,ndim,  fileIO); ABORT(ndim!=3              );
  READ(int,dim[0],fileIO); ABORT(dim[0]!=grid->nx+2   );
  READ(int,dim[1],fileIO); ABORT(dim[1]!=grid->ny+2   );
  READ(int,dim[2],fileIO); ABORT(dim[2]!=grid->nz+2   );

  MESSAGE(("Field dimensions: %d %d %d", dim[0], dim[1], dim[2]));
  fileIO.read( field_advance->f, dim[0]*dim[1]*dim[2] );
 
  // species ... species_list must be put together in the same order
  // it was originally in
 
  // FIXME: WHY IS NEW_SPECIES NOT CALLED?  (PROBABLY THE ABOVE
  // CONCERN ABOUT ORDERING ... RETOOL THIS TO CREATE THE
  // LIST AND THEN FLIP INTO CORRECT ORDERING?
 
  // FIXME: Crashes unnecessarily if we try restarting a
  // problem without any species.
  READ(int,n,fileIO); ABORT(n<1);
  for(last_sp=NULL;n;n--) {
    char * buf;
    // Note: sizeof(sp[0]) includes terminating '\0'
    READ(int,size,fileIO);                        ABORT(size<1 );
    MALLOC( buf, sizeof(sp[0])+size ); sp = (species_t *)buf;
	fileIO.read( sp->name, size);
    sp->name[size]='\0';
 
    READ( int,   sp->id,                fileIO ); ABORT(sp->id!=n-1 );
    READ( int,   sp->max_np,            fileIO ); ABORT(sp->max_np<1);
    READ( int,   sp->max_nm,            fileIO ); ABORT(sp->max_nm<1);
    READ( float, sp->q_m,               fileIO );
    READ( int,   sp->sort_interval,     fileIO );
    READ( int,   sp->sort_out_of_place, fileIO );
 
    MALLOC_ALIGNED( sp->p, sp->max_np, 128 );
    READ(int,size,  fileIO); ABORT(size!=sizeof(sp->p[0])       );
    READ(int,ndim,  fileIO); ABORT(ndim!=1                      );
    READ(int,dim[0],fileIO); ABORT(dim[0]<0 || dim[0]>sp->max_np);
    if( dim[0]>0 ) fileIO.read( sp->p, dim[0] );
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
 
  READ( int, grid->nb, fileIO ); ABORT( grid->nb<0 );
  if( grid->nb==0 ) grid->boundary = NULL;
  else {
    MALLOC( grid->boundary, grid->nb );
	fileIO.read( grid->boundary, grid->nb );
  }
 
  /****************************************
   * Restart-saved internal use variables *
   ****************************************/
 
  READ( double, p_time, fileIO ); ABORT(p_time<0);
  READ( double, g_time, fileIO ); ABORT(g_time<0);
  READ( double, u_time, fileIO ); ABORT(u_time<0);
  READ( double, f_time, fileIO ); ABORT(f_time<0);
 
  /****************************************
   * Restart-saved user defined variables *
   ****************************************/
 
  READ(int,size, fileIO); ABORT(size!=sizeof(user_global[0]));
  READ(int,ndim, fileIO); ABORT(ndim!=1                     );
  READ(int,dim[0],fileIO); ABORT(dim[0]!=USER_GLOBAL_SIZE   );
  fileIO.read( user_global, dim[0] );
 
  fileIO.close();
 
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

/*------------------------------------------------------------------------------
 * New dump logic
------------------------------------------------------------------------------*/

#include <iostream>

static FieldInfo fieldInfo[12] = {
	{ "Electric Field", "VECTOR", "3", "FLOATING_POINT", sizeof(float) },
	{ "Electric Field Divergence Error", "SCALAR", "1", "FLOATING_POINT",
		sizeof(float) },
	{ "Magnetic Field", "VECTOR", "3", "FLOATING_POINT", sizeof(float) },
	{ "Magnetic Field Divergence Error", "SCALAR", "1", "FLOATING_POINT",
		sizeof(float) },
	{ "TCA Field", "VECTOR", "3", "FLOATING_POINT", sizeof(float) },
	{ "Bound Charge Density", "SCALAR", "1", "FLOATING_POINT", sizeof(float) },
	{ "Free Current Field", "VECTOR", "3", "FLOATING_POINT", sizeof(float) },
	{ "Charge Density", "SCALAR", "1", "FLOATING_POINT", sizeof(float) },
	{ "Edge Material", "VECTOR", "3", "INTEGER", sizeof(material_id) },
	{ "Node Material", "SCALAR", "1", "INTEGER", sizeof(material_id) },
	{ "Face Material", "VECTOR", "3", "INTEGER", sizeof(material_id) },
	{ "Cell Material", "SCALAR", "1", "INTEGER", sizeof(material_id) }
}; // fieldInfo

static HydroInfo hydroInfo[5] = {
	{ "Current Density", "VECTOR", "3", "FLOATING_POINT", sizeof(float) },
	{ "Charge Density", "SCALAR", "1", "FLOATING_POINT", sizeof(float) },
	{ "Momentum Density", "VECTOR", "3", "FLOATING_POINT", sizeof(float) },
	{ "Kinetic Energy Density", "SCALAR", "1", "FLOATING_POINT",
		sizeof(float) },
	{ "Stress Tensor", "TENSOR", "6", "FLOATING_POINT", sizeof(float) }
	/*
	{ "STRESS_DIAGONAL", "VECTOR", "3", "FLOATING_POINT", sizeof(float) }
	{ "STRESS_OFFDIAGONAL", "VECTOR", "3", "FLOATING_POINT", sizeof(float) }
	*/
}; // hydroInfo

void vpic_simulation::create_field_list(char * strlist,
	DumpParameters & dumpParams) {

	strcpy(strlist, "");

	for(size_t i(0), pass(0); i<total_field_groups; i++) {
		if(dumpParams.output_vars.bitset(field_indeces[i])) {
			if(i>0 && pass) {
				strcat(strlist, ", ");
			}
			else {	
				pass = 1;
			} // if

			strcat(strlist, fieldInfo[i].name);

		} // if
	} // for
} // vpic_simulation::create_field_list

void vpic_simulation::create_hydro_list(char * strlist,
	DumpParameters & dumpParams) {

	strcpy(strlist, "");

	for(size_t i(0), pass(0); i<total_hydro_groups; i++) {
		if(dumpParams.output_vars.bitset(hydro_indeces[i])) {
			if(i>0 && pass) {
				strcat(strlist, ", ");
			}
			else {
				pass = 1;
			} // if
			strcat(strlist, hydroInfo[i].name);
		} // if
	} // for
} // vpic_simulation::create_field_list

void vpic_simulation::print_hashed_comment(FileIO & fileIO,
	const char * comment) {
		fileIO.print("################################################################################\n");
		fileIO.print("# %s\n", comment);
		fileIO.print("################################################################################\n");
} // vpic_simulation::print_hashed_comment

void vpic_simulation::global_header(const char * base,
	std::vector<DumpParameters *> dumpParams) {
	if(mp_rank(grid->mp) == 0) {
		/*
		 * Open the file for output
		 */
		char filename[256];
		sprintf(filename, "%s.vpc", base);

		FileIO fileIO;
		FileIOStatus status;

		status = fileIO.open(filename, io_write);
		if(status == fail) {
			ERROR(("Failed opening file: %s", filename));
		} // if

		print_hashed_comment(fileIO, "Header version information");
		fileIO.print("VPIC_HEADER_VERSION 1.0.0\n\n");

		print_hashed_comment(fileIO,
			"Header size for data file headers in bytes");
		fileIO.print("DATA_HEADER_SIZE 123\n\n");

		// Global grid inforation
		print_hashed_comment(fileIO, "Time step increment");
		fileIO.print("GRID_DELTA_T %f\n\n", grid->dt);

		print_hashed_comment(fileIO, "GRID_CVAC");
		fileIO.print("GRID_CVAC %f\n\n", grid->cvac);

		print_hashed_comment(fileIO, "GRID_EPS0");
		fileIO.print("GRID_EPS0 %f\n\n", grid->eps0);

		print_hashed_comment(fileIO, "Grid extents in the x-dimension");
		fileIO.print("GRID_EXTENTS_X %f %f\n\n", grid->x0, grid->x1);

		print_hashed_comment(fileIO, "Grid extents in the y-dimension");
		fileIO.print("GRID_EXTENTS_Y %f %f\n\n", grid->y0, grid->y1);

		print_hashed_comment(fileIO, "Grid extents in the z-dimension");
		fileIO.print("GRID_EXTENTS_Z %f %f\n\n", grid->z0, grid->z1);

		print_hashed_comment(fileIO, "Spatial step increment in x-dimension");
		fileIO.print("GRID_DELTA_X %f\n\n", grid->dx);

		print_hashed_comment(fileIO, "Spatial step increment in y-dimension");
		fileIO.print("GRID_DELTA_Y %f\n\n", grid->dy);

		print_hashed_comment(fileIO, "Spatial step increment in z-dimension");
		fileIO.print("GRID_DELTA_Z %f\n\n", grid->dz);

		print_hashed_comment(fileIO, "Domain partitions in x-dimension");
		fileIO.print("GRID_TOPOLOGY_X %d\n\n", px);

		print_hashed_comment(fileIO, "Domain partitions in y-dimension");
		fileIO.print("GRID_TOPOLOGY_Y %d\n\n", py);

		print_hashed_comment(fileIO, "Domain partitions in z-dimension");
		fileIO.print("GRID_TOPOLOGY_Z %d\n\n", pz);

		// Global data inforation
		assert(dumpParams.size() >= 2);

		print_hashed_comment(fileIO, "Field data information");
		fileIO.print("FIELD_DATA_DIRECTORY %s\n", dumpParams[0]->baseDir);
		fileIO.print("FIELD_DATA_BASE_FILENAME %s\n",
			dumpParams[0]->baseFileName);
		
		/*
		 * Create a variable list of field values to output.
		 */
		size_t numvars =
			std::min(dumpParams[0]->output_vars.bitsum(field_indeces,
			total_field_groups), total_field_groups);
		size_t * varlist = new size_t[numvars];
		for(size_t v(0), c(0); v<total_field_groups; v++) {
			if(dumpParams[0]->output_vars.bitset(field_indeces[v])) {
				varlist[c++] = v;
			} // if
		} // for

		// output variable list
		fileIO.print("FIELD_DATA_VARIABLES %d\n", numvars);

		for(size_t v(0); v<numvars; v++) {
			fileIO.print("\"%s\" %s %s %s %d\n", fieldInfo[varlist[v]].name,
			fieldInfo[varlist[v]].degree, fieldInfo[varlist[v]].elements,
			fieldInfo[varlist[v]].type, fieldInfo[varlist[v]].size);
		} // for

		fileIO.print("\n");

		delete varlist;
		varlist = NULL;

		/*
		 * Create a variable list for each species to output
		 */
		print_hashed_comment(fileIO, "Number of species with output data");
		fileIO.print("NUM_OUTPUT_SPECIES %d\n\n", dumpParams.size()-1);
		char species_comment[128];
		for(size_t i(1); i<dumpParams.size(); i++) {
			numvars =
				std::min(dumpParams[i]->output_vars.bitsum(hydro_indeces,
				total_hydro_groups), total_hydro_groups);

			sprintf(species_comment, "Species(%d) data information", (int)i);
			print_hashed_comment(fileIO, species_comment);
			fileIO.print("SPECIES_DATA_DIRECTORY %s\n",
				dumpParams[i]->baseDir);
			fileIO.print("SPECIES_DATA_BASE_FILENAME %s\n",
				dumpParams[i]->baseFileName);

			fileIO.print("HYDRO_DATA_VARIABLES %d\n", numvars);

			varlist = new size_t[numvars];
			for(size_t v(0), c(0); v<total_hydro_groups; v++) {
				if(dumpParams[i]->output_vars.bitset(hydro_indeces[v])) {
					varlist[c++] = v;
				} // if
			} // for

			for(size_t v(0); v<numvars; v++) {
				fileIO.print("\"%s\" %s %s %s %d\n", hydroInfo[varlist[v]].name,
				hydroInfo[varlist[v]].degree, hydroInfo[varlist[v]].elements,
				hydroInfo[varlist[v]].type, hydroInfo[varlist[v]].size);
			} // for


			delete[] varlist;
			varlist = NULL;

			if(i<dumpParams.size()-1) {
				fileIO.print("\n");
			} // if
		} // for

		fileIO.close();
	} // if
} // vpic_simulation::global_header

void vpic_simulation::field_dump(DumpParameters & dumpParams) {

	int32_t rank = mp_rank(grid->mp);

	/*
	 * Create directory for this time step
	 */
	char timeDir[256];
	sprintf(timeDir, "%s/T.%d", dumpParams.baseDir, step);
	dump_mkdir(timeDir);

	/*
	 * Open the file for output
	 */
 	char filename[256];
	sprintf(filename, "%s/T.%d/%s.%06d.%04d", dumpParams.baseDir, step,
		dumpParams.baseFileName, step, rank);

	FileIO fileIO;
	FileIOStatus status;

	status = fileIO.open(filename, io_write);
	if(status == fail) {
		if(rank == 0) {
			ERROR(("Failed opening file: %s", filename));
		} // if
	} // if

	// convenience
	const size_t istride(dumpParams.stride_x);
	const size_t jstride(dumpParams.stride_y);
	const size_t kstride(dumpParams.stride_z);

	/*
	 * Check stride values.
	 */
	if(remainder(grid->nx, istride) != 0) {
		if(rank == 0) {
			ERROR(("x stride must be an integer factor of nx"));
		} // if
	} // if
	if(remainder(grid->ny, jstride) != 0) {
		if(rank == 0) {
			ERROR(("y stride must be an integer factor of ny"));
		} // if
	} // if
	if(remainder(grid->nz, kstride) != 0) {
		if(rank == 0) {
			ERROR(("z stride must be an integer factor of nz"));
		} // if
	} // if

	int dim[3];

	/* IMPORTANT: these values are written in WRITE_HEADER_V0 */
	nxout = (grid->nx)/istride;
	nyout = (grid->ny)/jstride;
	nzout = (grid->nz)/kstride;
	dxout = (grid->dx)*istride;
	dyout = (grid->dy)*jstride;
	dzout = (grid->dz)*kstride;

	/* IMPORTANT: this depends on nxout, nyout, nzout */
	#define f(x,y,z) \
		f[INDEX_FORTRAN_3(x,y,z,0,nxout+1,0,nyout+1,0,nzout+1)]

	/*
	 * Banded output will write data as a single block-array as opposed to
	 * the Array-of-Structure format that is used for native storage.
	 *
	 * Additionally, the user can specify a stride pattern to reduce
	 * the resolution of the data that are output.  If a stride is
	 * specified for a particular dimension, VPIC will write the boundary
	 * plus every "stride" elements in that dimension.
	 */
	if(dumpParams.format == band) {

		WRITE_HEADER_V0(dump_type::field_dump, invalid_species_id, 0, fileIO);
	 
		dim[0] = nxout+2;
		dim[1] = nyout+2;
		dim[2] = nzout+2;

#if VERBOSE
		std::cerr << "nxout: " << nxout << std::endl;
		std::cerr << "nyout: " << nyout << std::endl;
		std::cerr << "nzout: " << nzout << std::endl;
		std::cerr << "nx: " << grid->nx << std::endl;
		std::cerr << "ny: " << grid->ny << std::endl;
		std::cerr << "nz: " << grid->nz << std::endl;
#endif

		WRITE_ARRAY_HEADER(field_advance->f, 3, dim, fileIO);

		/*
		 * Create a variable list of field values to output.
		 */
		size_t numvars = std::min(dumpParams.output_vars.bitsum(),
			total_field_variables);
		size_t * varlist = new size_t[numvars];

		for(size_t i(0), c(0); i<total_field_variables; i++) {
			if(dumpParams.output_vars.bitset(i)) { varlist[c++] = i; }
		} // for

		// more efficient for standard case
		if(istride == 1 && jstride == 1 && kstride == 1) {
			for(size_t v(0); v<numvars; v++) {
				for(size_t k(0); k<nzout+2; k++) {
					for(size_t j(0); j<nyout+2; j++) {
						for(size_t i(0); i<nxout+2; i++) {
							const uint32_t * fref =
								reinterpret_cast<uint32_t *>(
								&field_advance->f(i,j,k));
							fileIO.write(&fref[varlist[v]], 1);
#if VERBOSE
							if(rank==1 && v==0) {
								printf("%f ", field_advance->f(i,j,k).ex);
							} // if
#endif
						} // for
#if VERBOSE
						if(rank==1 && v==0) { printf("\nROW_BREAK\n"); }
#endif
					} // for
#if VERBOSE
					if(rank==1 && v==0) { printf("\nPLANE_BREAK\n"); }
#endif
				} // for
#if VERBOSE
				if(rank==1 && v==0) { printf("\nBLOCK_BREAK\n"); }
#endif
			} // for
		}
		else {
			for(size_t v(0); v<numvars; v++) {
				for(size_t k(0); k<nzout+2; k++) {
					const size_t koff = (k == 0) ? 0 : (k == nzout+1) ?
						grid->nz+1 : k*kstride-1;

					for(size_t j(0); j<nyout+2; j++) {
						const size_t joff = (j == 0) ? 0 : (j == nyout+1) ?
							grid->ny+1 : j*jstride-1;

						for(size_t i(0); i<nxout+2; i++) {
							const size_t ioff = (i == 0) ? 0 : (i == nxout+1) ?
								grid->nx+1 : i*istride-1;

#if VERBOSE
							if(rank == 0) {
								std::cout << "(" << ioff << " " <<
									joff << " " << koff << ")" << std::endl;
							} // if
#endif

							const uint32_t * fref =
								reinterpret_cast<uint32_t *>(
								&field_advance->f(ioff,joff,koff));
							fileIO.write(&fref[varlist[v]], 1);
						} // for
					} // for
				} // for
			} // for
		} // if

		delete[] varlist;
	}
	else { // band_interleave

		WRITE_HEADER_V0(dump_type::field_dump, invalid_species_id, 0, fileIO);
	 
		dim[0] = nxout+2;
		dim[1] = nyout+2;
		dim[2] = nzout+2;

		WRITE_ARRAY_HEADER(field_advance->f, 3, dim, fileIO);

		if(istride == 1 && jstride == 1 && kstride == 1) {
  			fileIO.write(field_advance->f, dim[0]*dim[1]*dim[2]);
		}
		else {
			for(size_t k(0); k<nzout+2; k++) {
				const size_t koff = (k == 0) ? 0 : (k == nzout+1) ?
					grid->nz+1 : k*kstride-1;

				for(size_t j(0); j<nyout+2; j++) {
					const size_t joff = (j == 0) ? 0 : (j == nyout+1) ?
						grid->ny+1 : j*jstride-1;

					for(size_t i(0); i<nxout+2; i++) {
						const size_t ioff = (i == 0) ? 0 : (i == nxout+1) ?
							grid->nx+1 : i*istride-1;

						fileIO.write(&field_advance->f(ioff,joff,koff), 1);
					} // for
				} // for
			} // for
		} // if
	} // if

	#undef f

	fileIO.close();
} // vpic_simulation::field_dump

void vpic_simulation::hydro_dump(const char * speciesname,
	DumpParameters & dumpParams) {

	/*
	 * Create directory for this time step
	 */
	char timeDir[256];
	sprintf(timeDir, "%s/T.%d", dumpParams.baseDir, step);
	dump_mkdir(timeDir);

	/*
	 * Open the file for output
	 */
 	char filename[256];
	sprintf(filename, "%s/T.%d/%s.%06d.%04d", dumpParams.baseDir, step,
		dumpParams.baseFileName, step, mp_rank(grid->mp));

	FileIO fileIO;
	FileIOStatus status;

	status = fileIO.open(filename, io_write);
	if(status == fail) {
		if(mp_rank(grid->mp) == 0) {
			ERROR(("Failed opening file: %s", filename));
		} // if
	} // if

	if(mp_rank(grid->mp) == 0) {
		MESSAGE(("Dumping %s species to %s", speciesname, filename));
	} // if

	species_t * sp = find_species_name(speciesname, species_list);
	if(sp == NULL) {
		ERROR(("Invalide species name: %s", speciesname));
	} // if

	clear_hydro(hydro, grid);
	accumulate_hydro_p(hydro, sp->p, sp->np, sp->q_m, interpolator, grid);
	synchronize_hydro(hydro, grid);

	// convenience
	const size_t istride(dumpParams.stride_x);
	const size_t jstride(dumpParams.stride_y);
	const size_t kstride(dumpParams.stride_z);

	/*
	 * Check stride values.
	 */
	if(remainder(grid->nx, istride) != 0) {
		if(mp_rank(grid->mp) == 0) {
			ERROR(("x stride must be an integer factor of nx"));
		} // if
	} // if
	if(remainder(grid->ny, jstride) != 0) {
		if(mp_rank(grid->mp) == 0) {
			ERROR(("y stride must be an integer factor of ny"));
		} // if
	} // if
	if(remainder(grid->nz, kstride) != 0) {
		if(mp_rank(grid->mp) == 0) {
			ERROR(("z stride must be an integer factor of nz"));
		} // if
	} // if

	int dim[3];

	/* IMPORTANT: these values are written in WRITE_HEADER_V0 */
	nxout = (grid->nx)/istride;
	nyout = (grid->ny)/jstride;
	nzout = (grid->nz)/kstride;
	dxout = (grid->dx)*istride;
	dyout = (grid->dy)*jstride;
	dzout = (grid->dz)*kstride;

	/* IMPORTANT: this depends on nxout, nyout, nzout */
	#define hydro(x,y,z) \
		hydro[INDEX_FORTRAN_3(x,y,z,0,nxout+1,0,nyout+1,0,nzout+1)]

	/*
	 * Banded output will write data as a single block-array as opposed to
	 * the Array-of-Structure format that is used for native storage.
	 *
	 * Additionally, the user can specify a stride pattern to reduce
	 * the resolution of the data that are output.  If a stride is
	 * specified for a particular dimension, VPIC will write the boundary
	 * plus every "stride" elements in that dimension.
	 */
	if(dumpParams.format == band) {

		WRITE_HEADER_V0(dump_type::hydro_dump, sp->id, sp->q_m, fileIO);
	 
		dim[0] = nxout+2;
		dim[1] = nyout+2;
		dim[2] = nzout+2;

		WRITE_ARRAY_HEADER(hydro, 3, dim, fileIO);

		/*
		 * Create a variable list of hydro values to output.
		 */
		size_t numvars = std::min(dumpParams.output_vars.bitsum(),
			total_hydro_variables);
		size_t * varlist = new size_t[numvars];
		for(size_t i(0), c(0); i<total_hydro_variables; i++) {
			if(dumpParams.output_vars.bitset(i)) { varlist[c++] = i;}
		} // for

		// More efficient for standard case
		if(istride == 1 && jstride == 1 && kstride == 1) {
			for(size_t v(0); v<numvars; v++) {
				for(size_t k(0); k<nzout+2; k++) {
					for(size_t j(0); j<nyout+2; j++) {
						for(size_t i(0); i<nxout+2; i++) {
							const uint32_t * href =
								reinterpret_cast<uint32_t *>(
								&hydro(i,j,k));
							fileIO.write(&href[varlist[v]], 1);
						} // for
					} // for
				} // for
			} // for
		}
		else {
			for(size_t v(0); v<numvars; v++) {
				for(size_t k(0); k<nzout+2; k++) {
					const size_t koff = (k == 0) ? 0 : (k == nzout+1) ?
						grid->nz+1 : k*kstride-1;

					for(size_t j(0); j<nyout+2; j++) {
						const size_t joff = (j == 0) ? 0 : (j == nyout+1) ?
							grid->ny+1 : j*jstride-1;

						for(size_t i(0); i<nxout+2; i++) {
							const size_t ioff = (i == 0) ? 0 : (i == nxout+1) ?
								grid->nx+1 : i*istride-1;

							const uint32_t * href =
								reinterpret_cast<uint32_t *>(
								&hydro(ioff,joff,koff));
							fileIO.write(&href[varlist[v]], 1);
						} // for
					} // for
				} // for
			} // for
		} // if

		delete[] varlist;
	}
	else { // band_interleave

		WRITE_HEADER_V0(dump_type::hydro_dump, sp->id, sp->q_m, fileIO);
	 
		dim[0] = nxout;
		dim[1] = nyout;
		dim[2] = nzout;

		WRITE_ARRAY_HEADER(hydro, 3, dim, fileIO);

		if(istride == 1 && jstride == 1 && kstride == 1) {
  			fileIO.write(hydro, dim[0]*dim[1]*dim[2]);
		}
		else {
			for(size_t k(0); k<nzout; k++) {
				const size_t koff = (k == 0) ? 0 : (k == nzout+1) ?
					grid->nz+1 : k*kstride-1;

				for(size_t j(0); j<nyout; j++) {
					const size_t joff = (j == 0) ? 0 : (j == nyout+1) ?
						grid->ny+1 : j*jstride-1;

					for(size_t i(0); i<nxout; i++) {
						const size_t ioff = (i == 0) ? 0 : (i == nxout+1) ?
							grid->nx+1 : i*istride-1;

						fileIO.write(&hydro(ioff,joff,koff), 1);
					} // for
				} // for
			} // for
		} // if
	} // if

	#undef hydro

	fileIO.close();
} // vpic_simulation::hydro_dump

#undef SETIVAR
#undef SETDVAR
#undef ITEST
#undef DTEST
