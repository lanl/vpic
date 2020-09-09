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

#include "vpic.h"
#include "../util/io/FileUtils.h"

#ifdef VPIC_ENABLE_HDF5
#include "hdf5.h" // from the lib
#endif

/* 1 means convert cell index to global */
#define OUTPUT_CONVERT_GLOBAL_ID 1

/* -1 means no ranks talk */
#define VERBOSE_rank -1

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

// Do not even define these functions when there are no IDs so that input decks that rely on them don't compile.
#ifdef VPIC_GLOBAL_PARTICLE_ID
int vpic_simulation::predicate_count(species_t* sp, std::function <bool (int)> f)
{
    if ((f == nullptr) || (!sp->has_ids)) {
      return sp->np;
    } else {
      return std::count_if( sp->p_id, sp->p_id + sp->np, f);
    }
}
int vpic_simulation::predicate_count(species_t* sp, std::function <bool (particle_t)> f)
{
    if (f != nullptr) return std::count_if( sp->p, sp->p + sp->np, f);
    else return sp->np;
}

void vpic_simulation::predicate_copy(species_t* sp_from, species_t* sp_to, std::function <bool (particle_t)> f)
{
    if (f != nullptr) {
        // std::copy_if( sp_from->p,    sp_from->p    + sp_from->np, sp_to->p   , f);
        // std::copy_if( sp_from->p_id, sp_from->p_id + sp_from->np, sp_to->p_id, f);
        // Manually loop over particles to do the 'cross copy' of p_id as well

        int next = 0; // track where we fill
        for (int i = 0; i < sp_from->np; i++)
        {
            if ( f(sp_from->p[i]) )
            {
                // copy i (inherently serial..)
                sp_to->p[next] = sp_from->p[i];
                if(sp_from->has_ids && sp_to->has_ids) {
                  sp_to->p_id[next] = sp_from->p_id[i];
                }
                #ifdef VPIC_PARTICLE_ANNOTATION
                if(sp_from->has_annotation && sp_to->has_annotation) {
                  for(int a = 0; a < sp_from->has_annotation && a < sp_to->has_annotation; a++) {
                    sp_to->p_annotation[next*sp_to->has_annotation + a] = sp_from->p_annotation[next*sp_from->has_annotation + a];
                  }
                }
                #endif
                next++;
            }

        }
        std::cout << "copied " << next << std::endl;
    }
}
void vpic_simulation::predicate_copy(species_t* sp_from, species_t* sp_to, std::function <bool (int)> f)
{
    if ((f != nullptr) && (sp_from->has_ids) )
    {
        //std::copy_if( sp->p_id, sp->p_id + sp->np, _sp.p, f);
        // Manually loop over particles to do the 'cross copy' from p_id->p

        int next = 0; // track where we fill
        for (int i = 0; i < sp_from->np; i++)
        {
            int this_id = sp_from->p_id[i];
            if ( f(this_id) )
            {
                // copy i (inherently serial..)
                sp_to->p[next] = sp_from->p[i];
                sp_to->p_id[next] = sp_from->p_id[i];
                #ifdef VPIC_PARTICLE_ANNOTATION
                if(sp_from->has_annotation && sp_to->has_annotation) {
                  for(int a = 0; a < sp_from->has_annotation && a < sp_to->has_annotation; a++) {
                    sp_to->p_annotation[next*sp_to->has_annotation + a] = sp_from->p_annotation[next*sp_from->has_annotation + a];
                  }
                }
                #endif
                next++;
            }

        }
        std::cout << "copied " << next << std::endl;
    }
}
#endif

void
vpic_simulation::dump_energies( const char *fname,
                                int append ) {
  double en_f[6], en_p;
  species_t *sp;
  FileIO fileIO;
  FileIOStatus status(fail);

  if( !fname ) ERROR(("Invalid file name"));

  if( rank()==0 ) {
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
      fileIO.print( "%li ", (long)step() );
    }
  }

  field_array->kernel->energy_f( en_f, field_array );
  if( rank()==0 && status!=fail )
    fileIO.print( "%e %e %e %e %e %e",
                  en_f[0], en_f[1], en_f[2],
                  en_f[3], en_f[4], en_f[5] );

  LIST_FOR_EACH(sp,species_list) {
    en_p = energy_p( sp, interpolator_array );
    if( rank()==0 && status!=fail ) fileIO.print( " %e", en_p );
  }

  if( rank()==0 && status!=fail ) {
    fileIO.print( "\n" );
    if( fileIO.close() ) ERROR(("File close failed on dump energies!!!"));
  }
}

// Note: dump_species/materials assume that names do not contain any \n!

void
vpic_simulation::dump_species( const char *fname ) {
  species_t *sp;
  FileIO fileIO;

  if( rank() ) return;
  if( !fname ) ERROR(( "Invalid file name" ));
  MESSAGE(( "Dumping species to \"%s\"", fname ));
  FileIOStatus status = fileIO.open(fname, io_write);
  if( status==fail ) ERROR(( "Could not open \"%s\".", fname ));
  LIST_FOR_EACH( sp, species_list )
    fileIO.print( "%s %i %e %e", sp->name, sp->id, sp->q, sp->m );
  if( fileIO.close() ) ERROR(( "File close failed on dump species!!!" ));
}

void
vpic_simulation::dump_materials( const char *fname ) {
  FileIO fileIO;
  material_t *m;
  if( rank() ) return;
  if( !fname ) ERROR(( "Invalid file name" ));
  MESSAGE(( "Dumping materials to \"%s\"", fname ));
  FileIOStatus status = fileIO.open(fname, io_write);
  if( status==fail ) ERROR(( "Could not open \"%s\"", fname ));
  LIST_FOR_EACH( m, material_list )
    fileIO.print( "%s\n%i\n%e %e %e\n%e %e %e\n%e %e %e\n",
                  m->name, m->id,
                  m->epsx,   m->epsy,   m->epsz,
                  m->mux,    m->muy,    m->muz,
                  m->sigmax, m->sigmay, m->sigmaz );
  if( fileIO.close() ) ERROR(( "File close failed on dump materials!!!" ));
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

void
vpic_simulation::dump_grid( const char *fbase ) {
  char fname[256];
  FileIO fileIO;
  int dim[4];

  if( !fbase ) ERROR(( "Invalid filename" ));
  if( rank()==0 ) MESSAGE(( "Dumping grid to \"%s\"", fbase ));

  sprintf( fname, "%s.%i", fbase, rank() );
  FileIOStatus status = fileIO.open(fname, io_write);
  if( status==fail ) ERROR(( "Could not open \"%s\".", fname ));

  /* IMPORTANT: these values are written in WRITE_HEADER_V0 */
  nxout = grid->nx;
  nyout = grid->ny;
  nzout = grid->nz;
  dxout = grid->dx;
  dyout = grid->dy;
  dzout = grid->dz;

  WRITE_HEADER_V0( dump_type::grid_dump, -1, 0, fileIO );

  dim[0] = 3;
  dim[1] = 3;
  dim[2] = 3;
  WRITE_ARRAY_HEADER( grid->bc, 3, dim, fileIO );
  fileIO.write( grid->bc, dim[0]*dim[1]*dim[2] );

  dim[0] = nproc()+1;
  WRITE_ARRAY_HEADER( grid->range, 1, dim, fileIO );
  fileIO.write( grid->range, dim[0] );

  dim[0] = 6;
  dim[1] = grid->nx+2;
  dim[2] = grid->ny+2;
  dim[3] = grid->nz+2;
  WRITE_ARRAY_HEADER( grid->neighbor, 4, dim, fileIO );
  fileIO.write( grid->neighbor, dim[0]*dim[1]*dim[2]*dim[3] );

  if( fileIO.close() ) ERROR(( "File close failed on dump grid!!!" ));
}

void
vpic_simulation::dump_fields( const char *fbase, int ftag ) {
  char fname[256];
  FileIO fileIO;
  int dim[3];

  if( !fbase ) ERROR(( "Invalid filename" ));

  if( rank()==0 ) MESSAGE(( "Dumping fields to \"%s\"", fbase ));

  if( ftag ) sprintf( fname, "%s.%li.%i", fbase, (long)step(), rank() );
  else       sprintf( fname, "%s.%i", fbase, rank() );

  FileIOStatus status = fileIO.open(fname, io_write);
  if( status==fail ) ERROR(( "Could not open \"%s\".", fname ));

  /* IMPORTANT: these values are written in WRITE_HEADER_V0 */
  nxout = grid->nx;
  nyout = grid->ny;
  nzout = grid->nz;
  dxout = grid->dx;
  dyout = grid->dy;
  dzout = grid->dz;

  WRITE_HEADER_V0( dump_type::field_dump, -1, 0, fileIO );

  dim[0] = grid->nx+2;
  dim[1] = grid->ny+2;
  dim[2] = grid->nz+2;
  WRITE_ARRAY_HEADER( field_array->f, 3, dim, fileIO );
  fileIO.write( field_array->f, dim[0]*dim[1]*dim[2] );
  if( fileIO.close() ) ERROR(( "File close failed on dump fields!!!" ));
}

void
vpic_simulation::dump_hydro( const char *sp_name,
                             const char *fbase,
                             int ftag ) {
  species_t *sp;
  char fname[256];
  FileIO fileIO;
  int dim[3];

  sp = find_species_name( sp_name, species_list );
  if( !sp ) ERROR(( "Invalid species \"%s\"", sp_name ));

  clear_hydro_array( hydro_array );
  accumulate_hydro_p( hydro_array, sp, interpolator_array );
  synchronize_hydro_array( hydro_array );

  if( !fbase ) ERROR(( "Invalid filename" ));

  if( rank()==0 )
    MESSAGE(("Dumping \"%s\" hydro fields to \"%s\"",sp->name,fbase));

  if( ftag ) sprintf( fname, "%s.%li.%i", fbase, (long)step(), rank() );
  else       sprintf( fname, "%s.%i", fbase, rank() );
  FileIOStatus status = fileIO.open(fname, io_write);
  if( status==fail) ERROR(( "Could not open \"%s\".", fname ));

  /* IMPORTANT: these values are written in WRITE_HEADER_V0 */
  nxout = grid->nx;
  nyout = grid->ny;
  nzout = grid->nz;
  dxout = grid->dx;
  dyout = grid->dy;
  dzout = grid->dz;

  WRITE_HEADER_V0( dump_type::hydro_dump,sp->id,sp->q/sp->m,fileIO);

  dim[0] = grid->nx+2;
  dim[1] = grid->ny+2;
  dim[2] = grid->nz+2;
  WRITE_ARRAY_HEADER( hydro_array->h, 3, dim, fileIO );
  fileIO.write( hydro_array->h, dim[0]*dim[1]*dim[2] );
  if( fileIO.close() ) ERROR(( "File close failed on dump hydro!!!" ));
}

void
vpic_simulation::dump_particles( const char *sp_name,
                                 const char *fbase,
                                 int ftag )
{
  species_t *sp;
  char fname[256];
  FileIO fileIO;
  int dim[2], buf_start;
  static particle_t * ALIGNED(128) p_buf = NULL;
# define PBUF_SIZE 32768 // 1MB of particles

  sp = find_species_name( sp_name, species_list );
  if( !sp ) ERROR(( "Invalid species name \"%s\".", sp_name ));

  if( !fbase ) ERROR(( "Invalid filename" ));

  if( !p_buf ) MALLOC_ALIGNED( p_buf, PBUF_SIZE, 128 );

  if( rank()==0 )
    MESSAGE(("Dumping \"%s\" particles to \"%s\"",sp->name,fbase));

  if( ftag ) sprintf( fname, "%s.%li.%i", fbase, (long)step(), rank() );
  else       sprintf( fname, "%s.%i", fbase, rank() );
  FileIOStatus status = fileIO.open(fname, io_write);
  if( status==fail ) ERROR(( "Could not open \"%s\"", fname ));

  /* IMPORTANT: these values are written in WRITE_HEADER_V0 */
  nxout = grid->nx;
  nyout = grid->ny;
  nzout = grid->nz;
  dxout = grid->dx;
  dyout = grid->dy;
  dzout = grid->dz;

  WRITE_HEADER_V0( dump_type::particle_dump, sp->id, sp->q/sp->m, fileIO );

  dim[0] = sp->np;
  WRITE_ARRAY_HEADER( p_buf, 1, dim, fileIO );

  // Copy a PBUF_SIZE hunk of the particle list into the particle
  // buffer, timecenter it and write it out. This is done this way to
  // guarantee the particle list unchanged while not requiring too
  // much memory.

  // FIXME: WITH A PIPELINED CENTER_P, PBUF NOMINALLY SHOULD BE QUITE
  // LARGE.

  particle_t * sp_p = sp->p;      sp->p      = p_buf;
  int sp_np         = sp->np;     sp->np     = 0;
  int sp_max_np     = sp->max_np; sp->max_np = PBUF_SIZE;
  for( buf_start=0; buf_start<sp_np; buf_start += PBUF_SIZE ) {
    sp->np = sp_np-buf_start; if( sp->np > PBUF_SIZE ) sp->np = PBUF_SIZE;
    COPY( sp->p, &sp_p[buf_start], sp->np );
    center_p( sp, interpolator_array );
    fileIO.write( sp->p, sp->np );
  }
  sp->p      = sp_p;
  sp->np     = sp_np;
  sp->max_np = sp_max_np;

  #ifdef VPIC_GLOBAL_PARTICLE_ID
  // append ID array at the end of the file
  if(sp->has_ids) {
    dim[0] = sp->np;
    WRITE_ARRAY_HEADER( sp->p_id, 1, dim, fileIO );
    // Maybe do this write in batched of PBUF_SIZE as well?
    fileIO.write(sp->p_id, sp->np);
  }
  #endif
  #ifdef VPIC_PARTICLE_ANNOTATION
  // append annotation buffer at the end of the file
  if(sp->has_annotation) {
    dim[0] = sp->np;
    dim[1] = sp->has_annotation;
    WRITE_ARRAY_HEADER( sp->p_annotation, 2, dim, fileIO );
    // Maybe do this write in batched of PBUF_SIZE as well?
    fileIO.write(sp->p_annotation, sp->np*sp->has_annotation);
  }
  #endif


  if( fileIO.close() ) ERROR(("File close failed on dump particles!!!"));
}

/*------------------------------------------------------------------------------
 * New dump logic
 *---------------------------------------------------------------------------*/

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

void
vpic_simulation::create_field_list( char * strlist,
                                    DumpParameters & dumpParams ) {
  strcpy(strlist, "");
  for(size_t i(0), pass(0); i<total_field_groups; i++)
    if(dumpParams.output_vars.bitset(field_indeces[i])) {
      if(i>0 && pass) strcat(strlist, ", ");
      else pass = 1;
      strcat(strlist, fieldInfo[i].name);
    }
}

void
vpic_simulation::create_hydro_list( char * strlist,
                                    DumpParameters & dumpParams ) {
  strcpy(strlist, "");
  for(size_t i(0), pass(0); i<total_hydro_groups; i++)
    if(dumpParams.output_vars.bitset(hydro_indeces[i])) {
      if(i>0 && pass) strcat(strlist, ", ");
      else pass = 1;
      strcat(strlist, hydroInfo[i].name);
    }
}

void
vpic_simulation::print_hashed_comment( FileIO & fileIO,
                                       const char * comment) {
  fileIO.print("################################################################################\n");
  fileIO.print("# %s\n", comment);
  fileIO.print("################################################################################\n");
}

void
vpic_simulation::global_header( const char * base,
                                std::vector<DumpParameters *> dumpParams ) {
  if( rank() ) return;

  // Open the file for output
  char filename[256];
  sprintf(filename, "%s.vpc", base);

  FileIO fileIO;
  FileIOStatus status;

  status = fileIO.open(filename, io_write);
  if(status == fail) ERROR(("Failed opening file: %s", filename));

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

  // Create a variable list of field values to output.
  size_t numvars = std::min(dumpParams[0]->output_vars.bitsum(field_indeces,
                                                              total_field_groups),
                            total_field_groups);
  size_t * varlist = new size_t[numvars];
  for(size_t v(0), c(0); v<total_field_groups; v++)
    if(dumpParams[0]->output_vars.bitset(field_indeces[v]))
      varlist[c++] = v;

  // output variable list
  fileIO.print("FIELD_DATA_VARIABLES %d\n", numvars);

  for(size_t v(0); v<numvars; v++)
    fileIO.print("\"%s\" %s %s %s %d\n", fieldInfo[varlist[v]].name,
                 fieldInfo[varlist[v]].degree, fieldInfo[varlist[v]].elements,
                 fieldInfo[varlist[v]].type, fieldInfo[varlist[v]].size);

  fileIO.print("\n");

  delete[] varlist;
  varlist = NULL;

  // Create a variable list for each species to output
  print_hashed_comment(fileIO, "Number of species with output data");
  fileIO.print("NUM_OUTPUT_SPECIES %d\n\n", dumpParams.size()-1);
  char species_comment[128];
  for(size_t i(1); i<dumpParams.size(); i++) {
    numvars = std::min(dumpParams[i]->output_vars.bitsum(hydro_indeces,
                                                         total_hydro_groups),
                       total_hydro_groups);

    sprintf(species_comment, "Species(%d) data information", (int)i);
    print_hashed_comment(fileIO, species_comment);
    fileIO.print("SPECIES_DATA_DIRECTORY %s\n",
                 dumpParams[i]->baseDir);
    fileIO.print("SPECIES_DATA_BASE_FILENAME %s\n",
                 dumpParams[i]->baseFileName);

    fileIO.print("HYDRO_DATA_VARIABLES %d\n", numvars);

    varlist = new size_t[numvars];
    for(size_t v(0), c(0); v<total_hydro_groups; v++)
      if(dumpParams[i]->output_vars.bitset(hydro_indeces[v]))
        varlist[c++] = v;

    for(size_t v(0); v<numvars; v++)
      fileIO.print("\"%s\" %s %s %s %d\n", hydroInfo[varlist[v]].name,
                   hydroInfo[varlist[v]].degree, hydroInfo[varlist[v]].elements,
                   hydroInfo[varlist[v]].type, hydroInfo[varlist[v]].size);


    delete[] varlist;
    varlist = NULL;

    if(i<dumpParams.size()-1) fileIO.print("\n");
  }


  if( fileIO.close() ) ERROR(( "File close failed on global header!!!" ));
}

void
vpic_simulation::field_dump( DumpParameters & dumpParams ) {

  // Create directory for this time step
  char timeDir[256];
  sprintf(timeDir, "%s/T.%ld", dumpParams.baseDir, (long)step());
  dump_mkdir(timeDir);

  // Open the file for output
  char filename[256];
  sprintf(filename, "%s/T.%ld/%s.%ld.%d", dumpParams.baseDir, (long)step(),
          dumpParams.baseFileName, (long)step(), rank());

  FileIO fileIO;
  FileIOStatus status;

  status = fileIO.open(filename, io_write);
  if( status==fail ) ERROR(( "Failed opening file: %s", filename ));

  // convenience
  const size_t istride(dumpParams.stride_x);
  const size_t jstride(dumpParams.stride_y);
  const size_t kstride(dumpParams.stride_z);

  // Check stride values.
  if(remainder(grid->nx, istride) != 0)
    ERROR(("x stride must be an integer factor of nx"));
  if(remainder(grid->ny, jstride) != 0)
    ERROR(("y stride must be an integer factor of ny"));
  if(remainder(grid->nz, kstride) != 0)
    ERROR(("z stride must be an integer factor of nz"));

  int dim[3];

  /* define to do C-style indexing */
# define f(x,y,z) f[ VOXEL(x,y,z, grid->nx,grid->ny,grid->nz) ]

  /* IMPORTANT: these values are written in WRITE_HEADER_V0 */
  nxout = (grid->nx)/istride;
  nyout = (grid->ny)/jstride;
  nzout = (grid->nz)/kstride;
  dxout = (grid->dx)*istride;
  dyout = (grid->dy)*jstride;
  dzout = (grid->dz)*kstride;

  /* Banded output will write data as a single block-array as opposed to
   * the Array-of-Structure format that is used for native storage.
   *
   * Additionally, the user can specify a stride pattern to reduce
   * the resolution of the data that are output.  If a stride is
   * specified for a particular dimension, VPIC will write the boundary
   * plus every "stride" elements in that dimension. */

  if(dumpParams.format == band) {

    WRITE_HEADER_V0(dump_type::field_dump, -1, 0, fileIO);

    dim[0] = nxout+2;
    dim[1] = nyout+2;
    dim[2] = nzout+2;

    if( rank()==VERBOSE_rank ) {
      std::cerr << "nxout: " << nxout << std::endl;
      std::cerr << "nyout: " << nyout << std::endl;
      std::cerr << "nzout: " << nzout << std::endl;
      std::cerr << "nx: " << grid->nx << std::endl;
      std::cerr << "ny: " << grid->ny << std::endl;
      std::cerr << "nz: " << grid->nz << std::endl;
    }

    WRITE_ARRAY_HEADER(field_array->f, 3, dim, fileIO);

    // Create a variable list of field values to output.
    size_t numvars = std::min(dumpParams.output_vars.bitsum(),
                              total_field_variables);
    size_t * varlist = new size_t[numvars];

    for(size_t i(0), c(0); i<total_field_variables; i++)
      if(dumpParams.output_vars.bitset(i)) varlist[c++] = i;

    if( rank()==VERBOSE_rank ) printf("\nBEGIN_OUTPUT\n");

    // more efficient for standard case
    if(istride == 1 && jstride == 1 && kstride == 1)
      for(size_t v(0); v<numvars; v++) {
      for(size_t k(0); k<nzout+2; k++) {
      for(size_t j(0); j<nyout+2; j++) {
      for(size_t i(0); i<nxout+2; i++) {
              const uint32_t * fref = reinterpret_cast<uint32_t *>(&field_array->f(i,j,k));
              fileIO.write(&fref[varlist[v]], 1);
              if(rank()==VERBOSE_rank) printf("%f ", field_array->f(i,j,k).ex);
              if(rank()==VERBOSE_rank) std::cout << "(" << i << " " << j << " " << k << ")" << std::endl;
      } if(rank()==VERBOSE_rank) std::cout << std::endl << "ROW_BREAK " << j << " " << k << std::endl;
      } if(rank()==VERBOSE_rank) std::cout << std::endl << "PLANE_BREAK " << k << std::endl;
      } if(rank()==VERBOSE_rank) std::cout << std::endl << "BLOCK_BREAK" << std::endl;
      }

    else

      for(size_t v(0); v<numvars; v++) {
      for(size_t k(0); k<nzout+2; k++) { const size_t koff = (k == 0) ? 0 : (k == nzout+1) ? grid->nz+1 : k*kstride-1;
      for(size_t j(0); j<nyout+2; j++) { const size_t joff = (j == 0) ? 0 : (j == nyout+1) ? grid->ny+1 : j*jstride-1;
      for(size_t i(0); i<nxout+2; i++) { const size_t ioff = (i == 0) ? 0 : (i == nxout+1) ? grid->nx+1 : i*istride-1;
              const uint32_t * fref = reinterpret_cast<uint32_t *>(&field_array->f(ioff,joff,koff));
              fileIO.write(&fref[varlist[v]], 1);
              if(rank()==VERBOSE_rank) printf("%f ", field_array->f(ioff,joff,koff).ex);
              if(rank()==VERBOSE_rank) std::cout << "(" << ioff << " " << joff << " " << koff << ")" << std::endl;
      } if(rank()==VERBOSE_rank) std::cout << std::endl << "ROW_BREAK " << joff << " " << koff << std::endl;
      } if(rank()==VERBOSE_rank) std::cout << std::endl << "PLANE_BREAK " << koff << std::endl;
      } if(rank()==VERBOSE_rank) std::cout << std::endl << "BLOCK_BREAK" << std::endl;
      }

    delete[] varlist;

  } else { // band_interleave

    WRITE_HEADER_V0(dump_type::field_dump, -1, 0, fileIO);

    dim[0] = nxout+2;
    dim[1] = nyout+2;
    dim[2] = nzout+2;

    WRITE_ARRAY_HEADER(field_array->f, 3, dim, fileIO);

    if(istride == 1 && jstride == 1 && kstride == 1)
      fileIO.write(field_array->f, dim[0]*dim[1]*dim[2]);
    else
      for(size_t k(0); k<nzout+2; k++) { const size_t koff = (k == 0) ? 0 : (k == nzout+1) ? grid->nz+1 : k*kstride-1;
      for(size_t j(0); j<nyout+2; j++) { const size_t joff = (j == 0) ? 0 : (j == nyout+1) ? grid->ny+1 : j*jstride-1;
      for(size_t i(0); i<nxout+2; i++) { const size_t ioff = (i == 0) ? 0 : (i == nxout+1) ? grid->nx+1 : i*istride-1;
            fileIO.write(&field_array->f(ioff,joff,koff), 1);
      }
      }
      }
  }

# undef f

  if( fileIO.close() ) ERROR(( "File close failed on field dump!!!" ));
}

void
vpic_simulation::hydro_dump( const char * speciesname,
                             DumpParameters & dumpParams ) {

  // Create directory for this time step
  char timeDir[256];
  sprintf(timeDir, "%s/T.%ld", dumpParams.baseDir, (long)step());
  dump_mkdir(timeDir);

  // Open the file for output
  char filename[256];
  sprintf( filename, "%s/T.%ld/%s.%ld.%d", dumpParams.baseDir, (long)step(),
           dumpParams.baseFileName, (long)step(), rank() );

  FileIO fileIO;
  FileIOStatus status;

  status = fileIO.open(filename, io_write);
  if(status == fail) ERROR(("Failed opening file: %s", filename));

  species_t * sp = find_species_name(speciesname, species_list);
  if( !sp ) ERROR(( "Invalid species name: %s", speciesname ));

  clear_hydro_array( hydro_array );
  accumulate_hydro_p( hydro_array, sp, interpolator_array );
  synchronize_hydro_array( hydro_array );

  // convenience
  const size_t istride(dumpParams.stride_x);
  const size_t jstride(dumpParams.stride_y);
  const size_t kstride(dumpParams.stride_z);

  // Check stride values.
  if(remainder(grid->nx, istride) != 0)
    ERROR(("x stride must be an integer factor of nx"));
  if(remainder(grid->ny, jstride) != 0)
    ERROR(("y stride must be an integer factor of ny"));
  if(remainder(grid->nz, kstride) != 0)
    ERROR(("z stride must be an integer factor of nz"));

  int dim[3];

  /* define to do C-style indexing */
# define hydro(x,y,z) hydro_array->h[VOXEL(x,y,z, grid->nx,grid->ny,grid->nz)]

  /* IMPORTANT: these values are written in WRITE_HEADER_V0 */
  nxout = (grid->nx)/istride;
  nyout = (grid->ny)/jstride;
  nzout = (grid->nz)/kstride;
  dxout = (grid->dx)*istride;
  dyout = (grid->dy)*jstride;
  dzout = (grid->dz)*kstride;

  /* Banded output will write data as a single block-array as opposed to
   * the Array-of-Structure format that is used for native storage.
   *
   * Additionally, the user can specify a stride pattern to reduce
   * the resolution of the data that are output.  If a stride is
   * specified for a particular dimension, VPIC will write the boundary
   * plus every "stride" elements in that dimension.
   */
  if(dumpParams.format == band) {

    WRITE_HEADER_V0(dump_type::hydro_dump, sp->id, sp->q/sp->m, fileIO);

    dim[0] = nxout+2;
    dim[1] = nyout+2;
    dim[2] = nzout+2;

    WRITE_ARRAY_HEADER(hydro_array->h, 3, dim, fileIO);

    /*
     * Create a variable list of hydro values to output.
     */
    size_t numvars = std::min(dumpParams.output_vars.bitsum(),
                              total_hydro_variables);
    size_t * varlist = new size_t[numvars];
    for(size_t i(0), c(0); i<total_hydro_variables; i++)
      if( dumpParams.output_vars.bitset(i) ) varlist[c++] = i;

    // More efficient for standard case
    if(istride == 1 && jstride == 1 && kstride == 1)

      for(size_t v(0); v<numvars; v++)
      for(size_t k(0); k<nzout+2; k++)
      for(size_t j(0); j<nyout+2; j++)
      for(size_t i(0); i<nxout+2; i++) {
              const uint32_t * href = reinterpret_cast<uint32_t *>(&hydro(i,j,k));
              fileIO.write(&href[varlist[v]], 1);
      }

    else

      for(size_t v(0); v<numvars; v++)
      for(size_t k(0); k<nzout+2; k++) { const size_t koff = (k == 0) ? 0 : (k == nzout+1) ? grid->nz+1 : k*kstride-1;
      for(size_t j(0); j<nyout+2; j++) { const size_t joff = (j == 0) ? 0 : (j == nyout+1) ? grid->ny+1 : j*jstride-1;
      for(size_t i(0); i<nxout+2; i++) { const size_t ioff = (i == 0) ? 0 : (i == nxout+1) ? grid->nx+1 : i*istride-1;
              const uint32_t * href = reinterpret_cast<uint32_t *>(&hydro(ioff,joff,koff));
              fileIO.write(&href[varlist[v]], 1);
      }
      }
      }

    delete[] varlist;

  } else { // band_interleave

    WRITE_HEADER_V0(dump_type::hydro_dump, sp->id, sp->q/sp->m, fileIO);

    dim[0] = nxout;
    dim[1] = nyout;
    dim[2] = nzout;

    WRITE_ARRAY_HEADER(hydro_array->h, 3, dim, fileIO);

    if(istride == 1 && jstride == 1 && kstride == 1)

      fileIO.write(hydro_array->h, dim[0]*dim[1]*dim[2]);

    else

      for(size_t k(0); k<nzout; k++) { const size_t koff = (k == 0) ? 0 : (k == nzout+1) ? grid->nz+1 : k*kstride-1;
      for(size_t j(0); j<nyout; j++) { const size_t joff = (j == 0) ? 0 : (j == nyout+1) ? grid->ny+1 : j*jstride-1;
      for(size_t i(0); i<nxout; i++) { const size_t ioff = (i == 0) ? 0 : (i == nxout+1) ? grid->nx+1 : i*istride-1;
            fileIO.write(&hydro(ioff,joff,koff), 1);
      }
      }
      }
  }

# undef hydro

  if( fileIO.close() ) ERROR(( "File close failed on hydro dump!!!" ));
}


void
vpic_simulation::init_buffered_particle_dump(const char * sp_name, const int N_timesteps, const double safety_factor) {
  species_t *sp = find_species_name( sp_name, species_list );
  if( !sp ) ERROR(( "Invalid species \"%s\"", sp_name ));

  // If any of those are still around from a previous init
  if(sp->output_buffer_dx) { FREE_ALIGNED(sp->output_buffer_dx); sp->output_buffer_dx = nullptr; }
  if(sp->output_buffer_dy) { FREE_ALIGNED(sp->output_buffer_dy); sp->output_buffer_dy = nullptr; }
  if(sp->output_buffer_dz) { FREE_ALIGNED(sp->output_buffer_dz); sp->output_buffer_dz = nullptr; }
  if(sp->output_buffer_i)  { FREE_ALIGNED(sp->output_buffer_i);  sp->output_buffer_i  = nullptr; }
  if(sp->output_buffer_ux) { FREE_ALIGNED(sp->output_buffer_ux); sp->output_buffer_ux = nullptr; }
  if(sp->output_buffer_uy) { FREE_ALIGNED(sp->output_buffer_uy); sp->output_buffer_uy = nullptr; }
  if(sp->output_buffer_uz) { FREE_ALIGNED(sp->output_buffer_uz); sp->output_buffer_uz = nullptr; }
  if(sp->output_buffer_w)  { FREE_ALIGNED(sp->output_buffer_w);  sp->output_buffer_w  = nullptr; }
  if(sp->output_buffer_id) { FREE_ALIGNED(sp->output_buffer_id); sp->output_buffer_id = nullptr; }
  if(sp->output_buffer_an) { FREE_ALIGNED(sp->output_buffer_an); sp->output_buffer_an = nullptr; }
  if(sp->output_buffer_ts) { FREE_ALIGNED(sp->output_buffer_ts); sp->output_buffer_ts = nullptr; }
  if(sp->buf_n_valid)      { FREE_ALIGNED(sp->buf_n_valid);      sp->buf_n_valid      = nullptr; }

  sp->buf_size = 0;

  // How many annotations per particles
  #ifdef VPIC_PARTICLE_ANNOTATION
    sp->buf_n_annotation = sp->has_annotation;
  #else
    sp->buf_n_annotation = 0;
  #endif

  // How many frames are we expected to store?
  sp->buf_n_frames = N_timesteps;

  // How many particles are we expected to store?
  sp->buf_n_particles = ceil(safety_factor * sp->max_np);

  // Buffer size
  sp->buf_size = sp->buf_n_frames * sp->buf_n_particles;

  // Allocate buffer
  fprintf(stderr, "<%d> Allocating buffers store up to %ld particles each in %ld timesteps\n", rank(), sp->buf_n_particles, sp->buf_n_frames);

  MALLOC_ALIGNED(sp->output_buffer_dx, sp->buf_size*sizeof(float), 128);
  MALLOC_ALIGNED(sp->output_buffer_dy, sp->buf_size*sizeof(float), 128);
  MALLOC_ALIGNED(sp->output_buffer_dz, sp->buf_size*sizeof(float), 128);
  MALLOC_ALIGNED(sp->output_buffer_i,  sp->buf_size*sizeof(int64_t), 128);
  MALLOC_ALIGNED(sp->output_buffer_ux, sp->buf_size*sizeof(float), 128);
  MALLOC_ALIGNED(sp->output_buffer_uy, sp->buf_size*sizeof(float), 128);
  MALLOC_ALIGNED(sp->output_buffer_uz, sp->buf_size*sizeof(float), 128);
  MALLOC_ALIGNED(sp->output_buffer_w,  sp->buf_size*sizeof(float), 128);
  #ifdef VPIC_GLOBAL_PARTICLE_ID
  if(sp->has_ids) {
    MALLOC_ALIGNED(sp->output_buffer_id, sp->buf_size*sizeof(size_t), 128);
  }
  #endif
  #ifdef VPIC_PARTICLE_ANNOTATION
  if(sp->has_annotation) {
    MALLOC_ALIGNED(sp->output_buffer_an, sp->buf_size*sp->buf_n_annotation*sizeof(float), 128);
  }
  #endif
  MALLOC_ALIGNED(sp->output_buffer_ts, sp->buf_size*sizeof(size_t), 128);
  MALLOC_ALIGNED(sp->buf_n_valid,      sp->buf_n_frames*sizeof(int64_t), 128);

  // Done
}

#ifdef OUTPUT_CONVERT_GLOBAL_ID
# define UNVOXEL(rank, ix, iy, iz, nx, ny, nz) BEGIN_PRIMITIVE {        \
int _ix, _iy, _iz;                                                  \
_ix  = (rank);        /* ix = ix+gpx*( iy+gpy*iz ) */       \
_iy  = _ix/int(nx);   /* iy = iy+gpy*iz */                  \
_ix -= _iy*int(nx);   /* ix = ix */                         \
_iz  = _iy/int(ny);   /* iz = iz */                         \
_iy -= _iz*int(ny);   /* iy = iy */                         \
(ix) = _ix;                                                         \
(iy) = _iy;                                                         \
(iz) = _iz;                                                         \
} END_PRIMITIVE
#endif

void
vpic_simulation::accumulate_buffered_particle_dump(const char * sp_name, const int frame) {
  species_t *sp = find_species_name( sp_name, species_list );
  if( !sp ) ERROR(( "Invalid species \"%s\"", sp_name ));

  if( !sp->buf_size ) ERROR(( "No buffer setup on species \"%s\"", sp_name));

  if(frame < 0 || frame >= sp->buf_n_frames) ERROR(( "Invalid frame %ld (has to be in (0;%ld( on species \"%s\"", frame, sp->buf_n_frames, sp_name ));

  if(sp->np > sp->buf_n_particles) WARNING(( "We have %d particles but buffer space for only %ld in species \"%s\". You will loose information. Chose a larger safety_factor next time.", sp->np, sp->buf_n_particles, sp_name));

  const int mpi_rank = rank();
  sp->buf_n_valid[frame] = 0;

  // Loop over particles
  for(int64_t n = 0; n < sp->np && n < sp->buf_n_particles; n++) {
    #ifdef OUTPUT_CONVERT_GLOBAL_ID
    // Convert Cell ID to global
    int64_t local_i = sp->p[n].i;

    int ix, iy, iz, rx, ry, rz;
    // Convert rank to local x/y/z
    UNVOXEL(mpi_rank, rx, ry, rz, grid->gpx, grid->gpy, grid->gpz);

    // Calculate local ix/iy/iz
    UNVOXEL(local_i, ix, iy, iz, grid->nx+2, grid->ny+2, grid->nz+2);

    // Account for the "first" ghost cell
    ix = ix - 1;
    iy = iy - 1;
    iz = iz - 1;

    // Convert ix/iy/iz to global
    int gix = ix + (grid->nx * rx);
    int giy = iy + (grid->ny * ry);
    int giz = iz + (grid->nz * rz);

    // calculate global grid sizes
    int gnx = grid->nx * grid->gpx;
    int gny = grid->ny * grid->gpy;
    int gnz = grid->nz * grid->gpz;

    // TODO: find a better way to account for the hard coded ghosts in VOXEL
    int64_t global_i = VOXEL(gix, giy, giz, gnx-2, gny-2, gnz-2);
    #endif

    const size_t index = frame*sp->buf_n_particles + n;

    // Particle Properties
    sp->output_buffer_dx[index] = sp->p[n].dx;
    sp->output_buffer_dy[index] = sp->p[n].dy;
    sp->output_buffer_dz[index] = sp->p[n].dz;

    #ifdef OUTPUT_CONVERT_GLOBAL_ID
    sp->output_buffer_i[index]  = global_i;
    #else
    sp->output_buffer_i[index]  = sp->p[n].i;
    #endif

    sp->output_buffer_ux[index] = sp->p[n].ux;
    sp->output_buffer_uy[index] = sp->p[n].uy;
    sp->output_buffer_uz[index] = sp->p[n].uz;
    sp->output_buffer_w[index]  = sp->p[n].w;

    #ifdef VPIC_GLOBAL_PARTICLE_ID
    if(sp->has_ids) {
      sp->output_buffer_id[index] = sp->p_id[n];
    }
    #endif

    // Annotations
    #ifdef VPIC_PARTICLE_ANNOTATION
    if(sp->has_annotation) {
      for(int a = 0; a < sp->has_annotation; a++) {
       const size_t out_index = sp->buf_n_particles * sp->buf_n_frames * a + n;
       const size_t in_index  = sp->has_annotation * n + a;
       sp->output_buffer_an[out_index] = sp->p_annotation[in_index];
      }
    }
    #endif

    sp->output_buffer_ts[index] = step();
    sp->buf_n_valid[frame]++;
  }
}

#undef UNVOXEL

void
vpic_simulation::write_buffered_particle_dump(const char * sp_name) {
  species_t *sp = find_species_name( sp_name, species_list );
  if( !sp ) ERROR(( "Invalid species \"%s\"", sp_name ));

#ifdef VPIC_ENABLE_HDF5
  char strbuf[4096];

  // MPI-Info object. necessary to create the property list for parallel file
  // open and to set IO hints
  MPI_Info info;
  MPI_Info_create(&info);
  MPI_Info_set(info, "romio_cb_read", "automatic");
  MPI_Info_set(info, "romio_cb_write", "automatic");
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, info);
  MPI_Info_free(&info);

  // Open file
  dump_mkdir("tracer");
  dump_mkdir("tracer/tracer1");
  sprintf(strbuf, "tracer/tracer1/T.%ld", step());
  dump_mkdir(strbuf);
  sprintf(strbuf, "tracer/tracer1/T.%ld/tracers.h5p", step());
  hid_t file = H5Fcreate(strbuf, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id); // FIXME: This limits us to a single species for now

  H5Pclose(plist_id);

  // Create group for the timestep
  sprintf(strbuf, "Step#%ld", step());
  hid_t group = H5Gcreate(file, strbuf, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hid_t subgroup = H5Gcreate(group, sp->name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // Find out which part of the global output falls to us
  int64_t total_entries, offset;
  int64_t local_entries = sp->buf_size;
  MPI_Allreduce(&local_entries, &total_entries, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Scan(&local_entries, &offset, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
  offset -= local_entries;

  plist_id = H5Pcreate(H5P_DATASET_XFER);
  if(total_entries > 0) {
    hid_t filespace, memspace, dset;
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    /* H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT); */
    filespace = H5Screate_simple(1, (hsize_t*) &total_entries, NULL);
    memspace =  H5Screate_simple(1, (hsize_t*) &local_entries, NULL);

    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, (hsize_t*)&offset, NULL, (hsize_t*)&local_entries, NULL);

    // dX
    dset = H5Dcreate(subgroup, "dX", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dset, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, sp->output_buffer_dx);
    H5Dclose(dset);
    // dY
    dset = H5Dcreate(subgroup, "dY", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dset, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, sp->output_buffer_dy);
    H5Dclose(dset);
    // dZ
    dset = H5Dcreate(subgroup, "dZ", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dset, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, sp->output_buffer_dz);
    H5Dclose(dset);
    // i
    dset = H5Dcreate(subgroup, "i", H5T_NATIVE_INT64, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dset, H5T_NATIVE_INT64, memspace, filespace, plist_id, sp->output_buffer_i);
    H5Dclose(dset);
    // Ux
    dset = H5Dcreate(subgroup, "Ux", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dset, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, sp->output_buffer_ux);
    H5Dclose(dset);
    // Uy
    dset = H5Dcreate(subgroup, "Uy", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dset, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, sp->output_buffer_uy);
    H5Dclose(dset);
    // Uz
    dset = H5Dcreate(subgroup, "Uz", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dset, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, sp->output_buffer_uz);
    H5Dclose(dset);
    // w
    dset = H5Dcreate(subgroup, "w", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dset, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, sp->output_buffer_w);
    H5Dclose(dset);
    // ID
    #ifdef VPIC_GLOBAL_PARTICLE_ID
    if(sp->has_ids) {
      dset = H5Dcreate(subgroup, "q", H5T_NATIVE_HSIZE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dwrite(dset, H5T_NATIVE_HSIZE, memspace, filespace, plist_id, sp->output_buffer_id);
      H5Dclose(dset);
    }
    #endif
    // dX
    #ifdef VPIC_PARTICLE_ANNOTATION
    if(sp->has_annotation) {
      for(int a = 0; a < sp->has_annotation; a++) {
        sprintf(strbuf, "Annotation_%d", a);
        float* data = &(sp->output_buffer_an[sp->buf_n_frames * sp->buf_n_particles * a]);
        dset = H5Dcreate(subgroup, strbuf, H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(dset, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, data);
        H5Dclose(dset);
      }
    }
    #endif
    // timestep
    dset = H5Dcreate(subgroup, "timestep", H5T_NATIVE_INT64, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dset, H5T_NATIVE_INT64, memspace, filespace, plist_id, sp->output_buffer_ts);
    H5Dclose(dset);

    H5Sclose(memspace);
    H5Sclose(filespace);
  }

  H5Gclose(subgroup);

  // Meta data
  hsize_t dcount[2], doffset[2], dset_dims[2];
  dcount[0] = sp->buf_n_frames;
  dcount[1] = 1;
  doffset[0] = 0;
  doffset[1] = rank();
  dset_dims[0] = dcount[0];
  dset_dims[1] = nproc();

  hid_t filespace = H5Screate_simple(2, dset_dims, NULL);
  hid_t memspace =  H5Screate_simple(2, dcount, NULL);
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, doffset, NULL, dcount, NULL);

  subgroup = H5Gcreate(group, "grid_meta_data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  sprintf(strbuf, "np_local_%s", sp->name);
  hid_t dset = H5Dcreate(subgroup, strbuf, H5T_NATIVE_INT64, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dset, H5T_NATIVE_INT64, memspace, filespace, plist_id, sp->buf_n_valid);
  H5Dclose(dset);

  H5Sclose(memspace);
  H5Sclose(filespace);
  H5Gclose(subgroup);

  // Cleanup
  H5Pclose(plist_id);
  H5Gclose(group);
  H5Fclose(file);
#else
  ERROR(("Only HDF5 format is current supported for buffered output of (annotated) particles\n"));
#endif
}

void
vpic_simulation::clear_buffered_particle_dump(const char * sp_name) {
  species_t *sp = find_species_name( sp_name, species_list );
  if( !sp ) ERROR(( "Invalid species \"%s\"", sp_name ));

  if(sp->buf_size) {
    memset(sp->output_buffer_dx, 0, sp->buf_size*sizeof(float));
    memset(sp->output_buffer_dy, 0, sp->buf_size*sizeof(float));
    memset(sp->output_buffer_dz, 0, sp->buf_size*sizeof(float));
    memset(sp->output_buffer_i,  0, sp->buf_size*sizeof(int64_t));
    memset(sp->output_buffer_ux, 0, sp->buf_size*sizeof(float));
    memset(sp->output_buffer_uy, 0, sp->buf_size*sizeof(float));
    memset(sp->output_buffer_uz, 0, sp->buf_size*sizeof(float));
    memset(sp->output_buffer_w,  0, sp->buf_size*sizeof(float));
    #ifdef VPIC_GLOBAL_PARTICLE_ID
    if(sp->has_ids) {
      memset(sp->output_buffer_id, 0, sp->buf_size*sizeof(size_t));
    }
    #endif
    #ifdef VPIC_PARTICLE_ANNOTATION
    if(sp->has_annotation) {
      memset(sp->output_buffer_an, 0, sp->buf_size*sp->buf_n_annotation*sizeof(float));
    }
    #endif

    // Set all timesteps to -1
    for(int64_t n = 0; n < sp->buf_n_frames * sp->buf_n_particles; n++) {
      sp->output_buffer_ts[n] = -1;
    }

    // No valid particles yet
    for(int n = 0; n < sp->buf_n_frames; n++) {
      sp->buf_n_valid[n] = 0;
    }
  }
}
