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

#include <cassert>

#include "vpic.h"
#include "dumpmacros.h"
#include "../util/io/FileUtils.h"

#ifdef VPIC_ENABLE_HDF5
#include "hdf5.h" // from the lib
#include "hdf5_header_info.h" // from vpic
#endif

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

// TODO: should this be an enum?
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

  WRITE_HEADER_V0( dump_type::grid_dump, -1, 0, step(), fileIO );

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
vpic_simulation::dump_fields( const char *fbase,
                              int ftag,
                              field_t *f )
{
  char fname[256];
  FileIO fileIO;
  int dim[3];

  if( !fbase ) ERROR(( "Invalid filename" ));

  if( rank()==0 ) MESSAGE(( "Dumping fields to \"%s\"", fbase ));

  if( ftag ) sprintf( fname, "%s.%li.%i", fbase, (long)step(), rank() );
  else       sprintf( fname, "%s.%i", fbase, rank() );

  FileIOStatus status = fileIO.open(fname, io_write);
  if( status==fail ) ERROR(( "Could not open \"%s\".", fname ));

  // default is to use VPIC native field data
  if ( f == NULL ) f = field_array->f;

  /* IMPORTANT: these values are written in WRITE_HEADER_V0 */
  nxout = grid->nx;
  nyout = grid->ny;
  nzout = grid->nz;
  dxout = grid->dx;
  dyout = grid->dy;
  dzout = grid->dz;

  WRITE_HEADER_V0( dump_type::field_dump, -1, 0, step(), fileIO );

  dim[0] = grid->nx+2;
  dim[1] = grid->ny+2;
  dim[2] = grid->nz+2;
  WRITE_ARRAY_HEADER( f, 3, dim, fileIO );
  fileIO.write( f, dim[0]*dim[1]*dim[2] );
  if( fileIO.close() ) ERROR(( "File close failed on dump fields." ));
}

void
vpic_simulation::dump_hydro( const char *sp_name,
                             const char *fbase,
                             int ftag,
                             hydro_t *h )
{
  species_t *sp;
  char fname[256];
  FileIO fileIO;
  int dim[3];

  sp = find_species_name( sp_name, species_list );
  if( !sp ) ERROR(( "Invalid species \"%s\"", sp_name ));

  // default behavior - do regular VPIC hydro dump

  if ( h == NULL)
  {
    h = hydro_array->h;
    clear_hydro_array( hydro_array );
    accumulate_hydro_p( hydro_array, sp, interpolator_array );
    synchronize_hydro_array( hydro_array );
  }

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

  WRITE_HEADER_V0( dump_type::hydro_dump, sp->id, sp->q / sp->m, step(), fileIO );

  dim[0] = grid->nx+2;
  dim[1] = grid->ny+2;
  dim[2] = grid->nz+2;
  WRITE_ARRAY_HEADER( h, 3, dim, fileIO );
  fileIO.write( h, dim[0]*dim[1]*dim[2] );
  if( fileIO.close() ) ERROR(( "File close failed on dump hydro." ));
}

#ifdef VPIC_ENABLE_HDF5
#define DUMP_DIR_FORMAT "./%s"

// TODO: rename or remove this
#define RANK_TO_INDEX2(rank, ix, iy, iz)                                      \
    BEGIN_PRIMITIVE                                                           \
    {                                                                         \
        int _ix, _iy, _iz;                                                    \
        _ix = (rank);                         /* ix = ix+gpx*( iy+gpy*iz ) */ \
        _iy = _ix / grid->gpx;  /* iy = iy+gpy*iz */            \
        _ix -= _iy * grid->gpx; /* ix = ix */                   \
        _iz = _iy / grid->gpy;  /* iz = iz */                   \
        _iy -= _iz * grid->gpy; /* iy = iy */                   \
        (ix) = _ix;                                                           \
        (iy) = _iy;                                                           \
        (iz) = _iz;                                                           \
    }                                                                         \
    END_PRIMITIVE

/* define to do C-style indexing */
#define hydro(x, y, z) hydro_array->h[VOXEL(x, y, z, grid->nx, grid->ny, grid->nz)]

// TODO: make function?
#define invert_field_xml_item(xml_file_name, speciesname_p, time_step, dims_4d, dims_3d, add_footer_flag)                                 \
  {                                                                                                                                       \
    FILE *fp;                                                                                                                             \
    fp = fopen(xml_file_name, "a");                                                                                                       \
    fprintf(fp, main_body_head, time_step);                                                                                               \
    if (field_dump_flag.enabledE())                                                                                                       \
      write_main_body_attribute(fp, main_body_attributeV, "E", dims_4d, dims_3d, speciesname_p, time_step, "ex", "ey", "ez");             \
    if (field_dump_flag.div_e_err)                                                                                                        \
      fprintf(fp, main_body_attributeS, "div_e_err", dims_3d, time_step, speciesname_p, time_step, time_step, "div_e_err");               \
    if (field_dump_flag.enabledCB())                                                                                                      \
      write_main_body_attribute(fp, main_body_attributeV, "B", dims_4d, dims_3d, speciesname_p, time_step, "cbx", "cby", "cbz");          \
    if (field_dump_flag.div_b_err)                                                                                                        \
      fprintf(fp, main_body_attributeS, "div_b_err", dims_3d, time_step, speciesname_p, time_step, time_step, "div_b_err");               \
    if (field_dump_flag.enabledTCA())                                                                                                     \
      write_main_body_attribute(fp, main_body_attributeV, "TCA", dims_4d, dims_3d, speciesname_p, time_step, "tcax", "tcay", "tcaz");     \
    if (field_dump_flag.rhob)                                                                                                             \
      fprintf(fp, main_body_attributeS, "rhob", dims_3d, time_step, speciesname_p, time_step, time_step, "rhob");                         \
    if (field_dump_flag.enabledJF())                                                                                                      \
      write_main_body_attribute(fp, main_body_attributeV, "JF", dims_4d, dims_3d, speciesname_p, time_step, "jfx", "jfy", "jfz");         \
    if (field_dump_flag.rhof)                                                                                                             \
      fprintf(fp, main_body_attributeS, "rhof", dims_3d, time_step, speciesname_p, time_step, time_step, "rhof");                         \
    if (field_dump_flag.enabledEMAT())                                                                                                    \
      write_main_body_attribute(fp, main_body_attributeV, "EMAT", dims_4d, dims_3d, speciesname_p, time_step, "ematx", "ematy", "ematz"); \
    if (field_dump_flag.nmat)                                                                                                             \
      fprintf(fp, main_body_attributeS, "nmat", dims_3d, time_step, speciesname_p, time_step, time_step, "nmat");                         \
    if (field_dump_flag.enabledFMAT())                                                                                                    \
      write_main_body_attribute(fp, main_body_attributeV, "FMAT", dims_4d, dims_3d, speciesname_p, time_step, "fmatx", "fmaty", "fmatz"); \
    if (field_dump_flag.cmat)                                                                                                             \
      fprintf(fp, main_body_attributeS, "cmat", dims_3d, time_step, speciesname_p, time_step, time_step, "cmat");                         \
    fprintf(fp, "%s", main_body_foot);                                                                                                          \
    if (add_footer_flag)                                                                                                                  \
      fputs(footer, fp);                                                                                                                  \
    fclose(fp);                                                                                                                           \
  }
void
vpic_simulation::dump_fields_hdf5( const char *fbase, int ftag )
{
    size_t step_for_viou = step();

    int mpi_size, mpi_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);


#ifdef DUMP_INFO_DEBUG
    printf("MPI rank = %d, size = %d \n", mpi_rank, mpi_size);
    printf("base dir for field: %s \n", global->fdParams.baseDir);
    printf("stride x y z  = (%ld, %ld, %ld)\n", global->fdParams.stride_x, global->fdParams.stride_y, global->fdParams.stride_z);
    printf("grid x, y z  = (%d, %d, %d) \n", grid->nx, grid->ny, grid->nz);
    printf("domain loc (x0, y0, z0) -> (x1, y1, z1) = (%f, %f, %f) -> (%f, %f, %f) \n", grid->x0, grid->y0, grid->z0, grid->x1, grid->y1, grid->z1);
    printf("global->topology_x, y, z =  %f, %f, %f \n ", global->topology_x, global->topology_y, global->topology_z);
    printf("grid -> sx, sy, sz =  (%d, %d, %d), nv=%d \n", grid->sx, grid->sy, grid->sz, grid->nv);
#endif

#define DUMP_FIELD_TO_HDF5(DSET_NAME, ATTRIBUTE_NAME, ELEMENT_TYPE)                                               \
    {                                                                                                             \
        dset_id = H5Dcreate(group_id, DSET_NAME, ELEMENT_TYPE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); \
        temp_buf_index = 0;                                                                                       \
        for (size_t i(1); i < grid->nx + 1; i++)                                                                  \
        {                                                                                                         \
            for (size_t j(1); j < grid->ny + 1; j++)                                                              \
            {                                                                                                     \
                for (size_t k(1); k < grid->nz + 1; k++)                                                          \
                {                                                                                                 \
                    temp_buf[temp_buf_index] = FIELD_ARRAY_NAME->fpp(i, j, k).ATTRIBUTE_NAME;                     \
                    temp_buf_index = temp_buf_index + 1;                                                          \
                }                                                                                                 \
            }                                                                                                     \
        }                                                                                                         \
        dataspace_id = H5Dget_space(dset_id);                                                                     \
        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, global_offset, NULL, global_count, NULL);               \
        H5Dwrite(dset_id, ELEMENT_TYPE, memspace, dataspace_id, plist_id, temp_buf);                              \
        H5Sclose(dataspace_id);                                                                                   \
        H5Dclose(dset_id);                                                                                        \
    }

    char fname[256];
    char field_scratch[128];
    char subfield_scratch[128];

    sprintf(field_scratch, DUMP_DIR_FORMAT, "field_hdf5");
    dump_mkdir(field_scratch);
    sprintf(subfield_scratch, "%s/T.%lld/", field_scratch, step_for_viou);
    dump_mkdir(subfield_scratch);

    sprintf(fname, "%s/%s_%lld.h5", subfield_scratch, "fields", step_for_viou);
    double el1 = uptime();
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
    hid_t file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    H5Pclose(plist_id);

    sprintf(fname, "Timestep_%lld", step_for_viou);
    hid_t group_id = H5Gcreate(file_id, fname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    el1 = uptime() - el1;
    //sim_log("TimeHDF5Open): " << el1 << " s"); //Easy to handle results for scripts
    double el2 = uptime();

    /*
// Create a variable list of field values to output.
size_t numvars = std::min(global->fdParams.output_vars.bitsum(), total_field_variables);
size_t * varlist = new size_t[numvars];

for(size_t i(0), c(0); i<total_field_variables; i++)
  if(global->fdParams.output_vars.bitset(i)) varlist[c++] = i;

printf("\nBEGIN_OUTPUT: numvars = %zd \n", numvars);*/

#define fpp(x, y, z) f[VOXEL(x, y, z, grid->nx, grid->ny, grid->nz)]
    /*
    typedef struct field {
    float ex,   ey,   ez,   div_e_err;     // Electric field and div E error
    float cbx,  cby,  cbz,  div_b_err;     // Magnetic field and div B error
    float tcax, tcay, tcaz, rhob;          // TCA fields and bound charge density
    float jfx,  jfy,  jfz,  rhof;          // Free current and charge density
    material_id ematx, ematy, ematz, nmat; // Material at edge centers and nodes
    material_id fmatx, fmaty, fmatz, cmat; // Material at face and cell centers
    } field_t;*/
    // Local voxel mesh resolution.  Voxels are
    // indexed FORTRAN style 0:nx+1,0:ny+1,0:nz+1
    // with voxels 1:nx,1:ny,1:nz being non-ghost
    // voxels.

    float *temp_buf = (float *)malloc(sizeof(float) * (grid->nx) * (grid->ny) * (grid->nz));
    hsize_t temp_buf_index;
    hid_t dset_id;
    //char  *field_var_name[] = {"ex","ey","ez","div_e_err","cbx","cby","cbz","div_b_err","tcax","tcay","tcaz","rhob","jfx","jfy","jfz","rhof"};
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    //Comment out for test only
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    //H5Sselect_hyperslab(filespace, H5S_SELECT_SET, (hsize_t *) &offset, NULL, (hsize_t *) &numparticles, NULL);

    //global->topology_x

    hsize_t field_global_size[3], field_local_size[3], global_offset[3], global_count[3];
    field_global_size[0] = (grid->nx * grid->gpx);
    field_global_size[1] = (grid->ny * grid->gpy);
    field_global_size[2] = (grid->nz * grid->gpz);

    field_local_size[0] = grid->nx;
    field_local_size[1] = grid->ny;
    field_local_size[2] = grid->nz;

    int gpx = grid->gpx;
    int gpy = grid->gpy;
    int gpz = grid->gpz;

    int mpi_rank_x, mpi_rank_y, mpi_rank_z;
    //RANK_TO_INDEX2(mpi_rank, mpi_rank_x, mpi_rank_y, mpi_rank_z);

    int _ix, _iy, _iz;
    _ix = (mpi_rank);
    _iy = _ix / grid->gpx;
    _ix -= _iy * grid->gpx;
    _iz = _iy / grid->gpy;
    _iy -= _iz * grid->gpy;
    int ix = _ix;
    int iy = _iy;
    int iz = _iz;

    mpi_rank_x = ix;
    mpi_rank_y = iy;
    mpi_rank_z = iz;

    global_offset[0] = (grid->nx) * mpi_rank_x;
    global_offset[1] = (grid->ny) * mpi_rank_y;
    global_offset[2] = (grid->nz) * mpi_rank_z;

    global_count[0] = (grid->nx);
    global_count[1] = (grid->ny);
    global_count[2] = (grid->nz);

#ifdef DUMP_INFO_DEBUG
    printf("global size   = " HSIZE_T ", " HSIZE_T ", " HSIZE_T "\n", field_global_size[0], field_global_size[1], field_global_size[2]);
    printf("global_offset = " HSIZE_T ", " HSIZE_T ", " HSIZE_T "\n", global_offset[0], global_offset[1], global_offset[2]);
    printf("global_count  = " HSIZE_T ", " HSIZE_T ", " HSIZE_T "\n", global_count[0], global_count[1], global_count[2]);
    printf("mpi-rank = %d, rank index = (%d, %d, %d) \n", mpi_rank, mpi_rank_x, mpi_rank_y, mpi_rank_z);
    fflush(stdout);
#endif

    hid_t filespace = H5Screate_simple(3, field_global_size, NULL);
    hid_t memspace = H5Screate_simple(3, field_local_size, NULL);
    hid_t dataspace_id;

    /*
    typedef struct field {
    float ex,   ey,   ez,   div_e_err;     // Electric field and div E error
    float cbx,  cby,  cbz,  div_b_err;     // Magnetic field and div B error
    float tcax, tcay, tcaz, rhob;          // TCA fields and bound charge density
    float jfx,  jfy,  jfz,  rhof;          // Free current and charge density
    material_id ematx, ematy, ematz, nmat; // Material at edge centers and nodes
    material_id fmatx, fmaty, fmatz, cmat; // Material at face and cell centers
    } field_t;*/

    if (field_dump_flag.ex)
        DUMP_FIELD_TO_HDF5("ex", ex, H5T_NATIVE_FLOAT);
    if (field_dump_flag.ey)
        DUMP_FIELD_TO_HDF5("ey", ey, H5T_NATIVE_FLOAT);
    if (field_dump_flag.ez)
        DUMP_FIELD_TO_HDF5("ez", ez, H5T_NATIVE_FLOAT);
    if (field_dump_flag.div_e_err)
        DUMP_FIELD_TO_HDF5("div_e_err", div_e_err, H5T_NATIVE_FLOAT);

    if (field_dump_flag.cbx)
        DUMP_FIELD_TO_HDF5("cbx", cbx, H5T_NATIVE_FLOAT);
    if (field_dump_flag.cby)
        DUMP_FIELD_TO_HDF5("cby", cby, H5T_NATIVE_FLOAT);
    if (field_dump_flag.cbz)
        DUMP_FIELD_TO_HDF5("cbz", cbz, H5T_NATIVE_FLOAT);
    if (field_dump_flag.div_b_err)
        DUMP_FIELD_TO_HDF5("div_b_err", div_b_err, H5T_NATIVE_FLOAT);

    if (field_dump_flag.tcax)
        DUMP_FIELD_TO_HDF5("tcax", tcax, H5T_NATIVE_FLOAT);
    if (field_dump_flag.tcay)
        DUMP_FIELD_TO_HDF5("tcay", tcay, H5T_NATIVE_FLOAT);
    if (field_dump_flag.tcaz)
        DUMP_FIELD_TO_HDF5("tcaz", tcaz, H5T_NATIVE_FLOAT);
    if (field_dump_flag.rhob)
        DUMP_FIELD_TO_HDF5("rhob", rhob, H5T_NATIVE_FLOAT);

    if (field_dump_flag.jfx)
        DUMP_FIELD_TO_HDF5("jfx", jfx, H5T_NATIVE_FLOAT);
    if (field_dump_flag.jfy)
        DUMP_FIELD_TO_HDF5("jfy", jfy, H5T_NATIVE_FLOAT);
    if (field_dump_flag.jfz)
        DUMP_FIELD_TO_HDF5("jfz", jfz, H5T_NATIVE_FLOAT);
    if (field_dump_flag.rhof)
        DUMP_FIELD_TO_HDF5("rhof", rhof, H5T_NATIVE_FLOAT);

    //H5T_NATIVE_SHORT  for material_id (typedef int16_t material_id)
    if (field_dump_flag.ematx)
        DUMP_FIELD_TO_HDF5("ematx", ematx, H5T_NATIVE_SHORT);
    if (field_dump_flag.ematy)
        DUMP_FIELD_TO_HDF5("ematy", ematy, H5T_NATIVE_SHORT);
    if (field_dump_flag.ematz)
        DUMP_FIELD_TO_HDF5("ematz", ematz, H5T_NATIVE_SHORT);
    if (field_dump_flag.nmat)
        DUMP_FIELD_TO_HDF5("nmat", nmat, H5T_NATIVE_SHORT);

    if (field_dump_flag.fmatx)
        DUMP_FIELD_TO_HDF5("fmatx", fmatx, H5T_NATIVE_SHORT);
    if (field_dump_flag.fmaty)
        DUMP_FIELD_TO_HDF5("fmaty", fmaty, H5T_NATIVE_SHORT);
    if (field_dump_flag.fmatz)
        DUMP_FIELD_TO_HDF5("fmatz", fmatz, H5T_NATIVE_SHORT);
    if (field_dump_flag.cmat)
        DUMP_FIELD_TO_HDF5("cmat", cmat, H5T_NATIVE_SHORT);

    el2 = uptime() - el2;
    //sim_log("TimeHDF5Write: " << el2 << " s");

    double el3 = uptime();

    //Write metadata (geo original and geo dx/dy/dz) for ArrayUDF
    float attr_data[2][3];
    attr_data[0][0] = grid->x0;
    attr_data[0][1] = grid->y0;
    attr_data[0][2] = grid->z0;
    attr_data[1][0] = grid->dx;
    attr_data[1][1] = grid->dy;
    attr_data[1][2] = grid->dz;
    hsize_t dims[2];
    dims[0] = 2;
    dims[1] = 3;
    hid_t va_geo_dataspace_id = H5Screate_simple(2, dims, NULL);
    hid_t va_geo_attribute_id = H5Acreate2(file_id, "VPIC-ArrayUDF-GEO", H5T_IEEE_F32BE, va_geo_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(va_geo_attribute_id, H5T_NATIVE_FLOAT, attr_data);
    H5Sclose(va_geo_dataspace_id);
    H5Aclose(va_geo_attribute_id);

    free(temp_buf);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Pclose(plist_id);
    H5Gclose(group_id);
    H5Fclose(file_id);

    el3 = uptime() - el3;
    //sim_log("TimeHDF5Close: " << el3 << " s");

    if (mpi_rank == 0)
    {
        char const *output_xml_file = "./field_hdf5/hdf5_field.xdmf";
        char dimensions_3d[128];
        sprintf(dimensions_3d, "%lld %lld %lld", field_global_size[0], field_global_size[1], field_global_size[2]);
        char dimensions_4d[128];
        sprintf(dimensions_4d, "%lld %lld %lld %d", field_global_size[0], field_global_size[1], field_global_size[2], 3);
        char orignal[128];
        sprintf(orignal, "%f %f %f", grid->x0, grid->y0, grid->z0);
        char dxdydz[128];
        sprintf(dxdydz, "%f %f %f", grid->dx, grid->dy, grid->dz);

        //int fields_interval = global->fields_interval;
        // TODO: make sure field interval is set
        int nframes = num_step / field_interval + 1;
        static int field_tframe = 0;

#ifdef DUMP_INFO_DEBUG
        printf("         meta file : %s \n", output_xml_file);
        printf(" array dims per var: %s \n", dimensions_3d);
        printf("array dims all vars: %s \n", dimensions_4d);
        printf("            orignal: %s \n", orignal);
        printf("             dxdydz: %s \n", dxdydz);
        printf("            nframes: %d \n", nframes);
        printf("    field_interval: %d \n", field_interval);
        printf("       current step: %lld \n", step_for_viou);
                printf("       current step: %lld \n", step_for_viou);

        //printf("    Simulation time: %f \n", grid->t0);
        printf("             tframe: %d \n", field_tframe);
#endif

        if (field_tframe >= 1)
        {
            if (field_tframe == (nframes - 1))
            {
                invert_field_xml_item(output_xml_file, "fields", step_for_viou, dimensions_4d, dimensions_3d, 1);
            }
            else
            {
                invert_field_xml_item(output_xml_file, "fields", step_for_viou, dimensions_4d, dimensions_3d, 0);
            }
        }
        else
        {
            create_file_with_header(output_xml_file, dimensions_3d, orignal, dxdydz, nframes, field_interval);
            if (field_tframe == (nframes - 1))
            {
                invert_field_xml_item(output_xml_file, "fields", step_for_viou, dimensions_4d, dimensions_3d, 1);
            }
            else
            {
                invert_field_xml_item(output_xml_file, "fields", step_for_viou, dimensions_4d, dimensions_3d, 0);
            }
        }
        field_tframe++;
    }
}
void vpic_simulation::dump_hydro_hdf5( const char *speciesname,
                             const char *fbase,
                             int ftag )
{
    size_t step_for_viou = step();

#define DUMP_HYDRO_TO_HDF5(DSET_NAME, ATTRIBUTE_NAME, ELEMENT_TYPE)                                               \
    {                                                                                                             \
        dset_id = H5Dcreate(group_id, DSET_NAME, ELEMENT_TYPE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); \
        temp_buf_index = 0;                                                                                       \
        for (size_t i(1); i < grid->nx + 1; i++)                                                                  \
        {                                                                                                         \
            for (size_t j(1); j < grid->ny + 1; j++)                                                              \
            {                                                                                                     \
                for (size_t k(1); k < grid->nz + 1; k++)                                                          \
                {                                                                                                 \
                    temp_buf[temp_buf_index] = hydro(i, j, k).ATTRIBUTE_NAME;                                     \
                    temp_buf_index = temp_buf_index + 1;                                                          \
                }                                                                                                 \
            }                                                                                                     \
        }                                                                                                         \
        dataspace_id = H5Dget_space(dset_id);                                                                     \
        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, global_offset, NULL, global_count, NULL);               \
        H5Dwrite(dset_id, ELEMENT_TYPE, memspace, dataspace_id, plist_id, temp_buf);                              \
        H5Sclose(dataspace_id);                                                                                   \
        H5Dclose(dset_id);                                                                                        \
    }
    //#define DUMP_INFO_DEBUG 1
    int mpi_size, mpi_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    species_t *sp = find_species_name(speciesname, species_list);
    if (!sp)
        ERROR(("Invalid species name: %s", speciesname));

#ifdef ENABLE_V407_SCIDAC
    clear_hydro( hydro, grid );
    accumulate_hydro_p( hydro, sp->p, sp->np, sp->q_m, interpolator, grid );
    synchronize_hydro( hydro, grid );
#else
    clear_hydro_array(hydro_array);
    accumulate_hydro_p(hydro_array, sp, interpolator_array);
    synchronize_hydro_array(hydro_array);
#endif
    /*#ifdef DUMP_INFO_DEBUG
printf("MPI rank = %d, size = %d \n", mpi_rank, mpi_size);
printf("base dir for field: %s \n", global->fdParams.baseDir);
printf("stride x y z  = (%ld, %ld, %ld)\n", global->fdParams.stride_x, global->fdParams.stride_y, global->fdParams.stride_z);
printf("grid x, y z  = (%d, %d, %d) \n", grid->nx, grid->ny, grid->nz);
printf("domain loc (x0, y0, z0) -> (x1, y1, z1) = (%f, %f, %f) -> (%f, %f, %f) \n", grid->x0, grid->y0, grid->z0, grid->x1, grid->y1, grid->z1);
printf("global->topology_x, y, z =  %f, %f, %f \n ", global->topology_x, global->topology_y, global->topology_z);
printf("grid -> sx, sy, sz =  (%d, %d, %d), nv=%d \n", grid->sx, grid->sy, grid->sz, grid->nv);
#endif*/

    char hname[256];
    char hydro_scratch[128];
    char subhydro_scratch[128];

    sprintf(hydro_scratch, "./%s", "hydro_hdf5");
    dump_mkdir(hydro_scratch);
    sprintf(subhydro_scratch, "%s/T.%lld/", hydro_scratch, step_for_viou);
    dump_mkdir(subhydro_scratch);

    sprintf(hname, "%s/hydro_%s_%lld.h5", subhydro_scratch, speciesname, step_for_viou);
    double el1 = uptime();
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
    hid_t file_id = H5Fcreate(hname, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    H5Pclose(plist_id);

    sprintf(hname, "Timestep_%lld", step_for_viou);
    hid_t group_id = H5Gcreate(file_id, hname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    el1 = uptime() - el1;
    //sim_log("TimeHDF5Open: " << el1 << " s"); //Easy to handle results for scripts
    double el2 = uptime();

    // Create a variable list of field values to output.
    //size_t numvars = std::min(global->fdParams.output_vars.bitsum(), total_field_variables);
    //size_t *varlist = new size_t[numvars];

    //for (size_t i(0), c(0); i < total_field_variables; i++)
    //    if (global->fdParams.output_vars.bitset(i))
    //        varlist[c++] = i;

    //printf("\nBEGIN_OUTPUT: numvars = %zd \n", numvars);


    //typedef struct hydro {
    //  float jx, jy, jz, rho; // Current and charge density => <q v_i f>, <q f>
    //  float px, py, pz, ke;  // Momentum and K.E. density  => <p_i f>, <m c^2 (gamma-1) f>
    //  float txx, tyy, tzz;   // Stress diagonal            => <p_i v_j f>, i==j
    //  float tyz, tzx, txy;   // Stress off-diagonal        => <p_i v_j f>, i!=j
    //  float _pad[2];         // 16-byte align
    //} hydro_t;

    //typedef struct hydro_array {
    //  hydro_t * ALIGNED(128) h;
    //  grid_t * g;
    //} hydro_array_t;

    float *temp_buf = (float *)malloc(sizeof(float) * (grid->nx) * (grid->ny) * (grid->nz));
    hsize_t temp_buf_index;
    hid_t dset_id;
    //char  *field_var_name[] = {"ex","ey","ez","div_e_err","cbx","cby","cbz","div_b_err","tcax","tcay","tcaz","rhob","jfx","jfy","jfz","rhof"};
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    //Comment out for test only
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    //H5Sselect_hyperslab(filespace, H5S_SELECT_SET, (hsize_t *) &offset, NULL, (hsize_t *) &numparticles, NULL);

    //global->topology_x

    hsize_t hydro_global_size[3], hydro_local_size[3], global_offset[3], global_count[3];
    hydro_global_size[0] = (grid->nx * grid->gpx);
    hydro_global_size[1] = (grid->ny * grid->gpy);
    hydro_global_size[2] = (grid->nz * grid->gpz);

    hydro_local_size[0] = grid->nx;
    hydro_local_size[1] = grid->ny;
    hydro_local_size[2] = grid->nz;

    int mpi_rank_x, mpi_rank_y, mpi_rank_z;
    RANK_TO_INDEX2(mpi_rank, mpi_rank_x, mpi_rank_y, mpi_rank_z);

    global_offset[0] = (grid->nx) * mpi_rank_x;
    global_offset[1] = (grid->ny) * mpi_rank_y;
    global_offset[2] = (grid->nz) * mpi_rank_z;

    global_count[0] = (grid->nx);
    global_count[1] = (grid->ny);
    global_count[2] = (grid->nz);

#ifdef DUMP_INFO_DEBUG
    printf("global size   = " HSIZE_T ", " HSIZE_T ", " HSIZE_T "\n", hydro_global_size[0], hydro_global_size[1], hydro_global_size[2]);
    printf("global_offset = " HSIZE_T ", " HSIZE_T ", " HSIZE_T "\n", global_offset[0], global_offset[1], global_offset[2]);
    printf("global_count  = " HSIZE_T ", " HSIZE_T ", " HSIZE_T "\n", global_count[0], global_count[1], global_count[2]);
    printf("mpi-rank = %d, rank index = (%d, %d, %d) \n", mpi_rank, mpi_rank_x, mpi_rank_y, mpi_rank_z);
    fflush(stdout);
#endif

    hid_t filespace = H5Screate_simple(3, hydro_global_size, NULL);
    hid_t memspace = H5Screate_simple(3, hydro_local_size, NULL);
    hid_t dataspace_id;

    //typedef struct hydro {
    //  float jx, jy, jz, rho; // Current and charge density => <q v_i f>, <q f>
    //  float px, py, pz, ke;  // Momentum and K.E. density  => <p_i f>, <m c^2 (gamma-1) f>
    //  float txx, tyy, tzz;   // Stress diagonal            => <p_i v_j f>, i==j
    //  float tyz, tzx, txy;   // Stress off-diagonal        => <p_i v_j f>, i!=j
    //  float _pad[2];         // 16-byte align
    //} hydro_t;

    if (hydro_dump_flag.jx)
        DUMP_HYDRO_TO_HDF5("jx", jx, H5T_NATIVE_FLOAT);
    if (hydro_dump_flag.jy)
        DUMP_HYDRO_TO_HDF5("jy", jy, H5T_NATIVE_FLOAT);
    if (hydro_dump_flag.jz)
        DUMP_HYDRO_TO_HDF5("jz", jz, H5T_NATIVE_FLOAT);
    if (hydro_dump_flag.rho)
        DUMP_HYDRO_TO_HDF5("rho", rho, H5T_NATIVE_FLOAT);

    if (hydro_dump_flag.px)
        DUMP_HYDRO_TO_HDF5("px", px, H5T_NATIVE_FLOAT);
    if (hydro_dump_flag.py)
        DUMP_HYDRO_TO_HDF5("py", py, H5T_NATIVE_FLOAT);
    if (hydro_dump_flag.pz)
        DUMP_HYDRO_TO_HDF5("pz", pz, H5T_NATIVE_FLOAT);
    if (hydro_dump_flag.ke)
        DUMP_HYDRO_TO_HDF5("ke", ke, H5T_NATIVE_FLOAT);

    if (hydro_dump_flag.txx)
        DUMP_HYDRO_TO_HDF5("txx", txx, H5T_NATIVE_FLOAT);
    if (hydro_dump_flag.tyy)
        DUMP_HYDRO_TO_HDF5("tyy", tyy, H5T_NATIVE_FLOAT);
    if (hydro_dump_flag.tzz)
        DUMP_HYDRO_TO_HDF5("tzz", tzz, H5T_NATIVE_FLOAT);

    if (hydro_dump_flag.tyz)
        DUMP_HYDRO_TO_HDF5("tyz", tyz, H5T_NATIVE_FLOAT);
    if (hydro_dump_flag.tzx)
        DUMP_HYDRO_TO_HDF5("tzx", tzx, H5T_NATIVE_FLOAT);
    if (hydro_dump_flag.txy)
        DUMP_HYDRO_TO_HDF5("txy", txy, H5T_NATIVE_FLOAT);

    el2 = uptime() - el2;
    //sim_log("TimeHDF5Write: " << el2 << " s");

    double el3 = uptime();

    //Write metadata (geo original and geo dx/dy/dz) for ArrayUDF
    float attr_data[2][3];
    attr_data[0][0] = grid->x0;
    attr_data[0][1] = grid->y0;
    attr_data[0][2] = grid->z0;
    attr_data[1][0] = grid->dx;
    attr_data[1][1] = grid->dy;
    attr_data[1][2] = grid->dz;
    hsize_t dims[2];
    dims[0] = 2;
    dims[1] = 3;
    hid_t va_geo_dataspace_id = H5Screate_simple(2, dims, NULL);
    hid_t va_geo_attribute_id = H5Acreate2(file_id, "VPIC-ArrayUDF-GEO", H5T_IEEE_F32BE, va_geo_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(va_geo_attribute_id, H5T_NATIVE_FLOAT, attr_data);
    H5Sclose(va_geo_dataspace_id);
    H5Aclose(va_geo_attribute_id);

    free(temp_buf);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Pclose(plist_id);
    H5Gclose(group_id);
    H5Fclose(file_id);

    el3 = uptime() - el3;
    //sim_log("TimeHDF5Close: " << el3 << " s");

    if (mpi_rank == 0)
    {
        char output_xml_file[128];
        sprintf(output_xml_file, "./%s/%s%s%s", "hydro_hdf5", "hydro-", speciesname, ".xdmf");
        char dimensions_3d[128];
        sprintf(dimensions_3d, "%lld %lld %lld", hydro_global_size[0], hydro_global_size[1], hydro_global_size[2]);
        char dimensions_4d[128];
        sprintf(dimensions_4d, "%lld %lld %lld %d", hydro_global_size[0], hydro_global_size[1], hydro_global_size[2], 3);
        char orignal[128];
        sprintf(orignal, "%f %f %f", grid->x0, grid->y0, grid->z0);
        char dxdydz[128];
        sprintf(dxdydz, "%f %f %f", grid->dx, grid->dy, grid->dz);

        int nframes = num_step / field_interval + 1;
        int fields_interval = field_interval;
        static int tframe = 0;

#ifdef DUMP_INFO_DEBUG
        printf("         meta file : %s \n", output_xml_file);
        printf(" array dims per var: %s \n", dimensions_3d);
        printf("array dims all vars: %s \n", dimensions_4d);
        printf("            orignal: %s \n", orignal);
        printf("             dxdydz: %s \n", dxdydz);
        printf("            nframes: %d \n", nframes);
        printf("    fields_interval: %d \n", fields_interval);
        printf("       current step: %lld \n", step_for_viou);
        printf("    Simulation time: %f \n", grid->t0);
        printf("             tframe: %d \n", tframe);
#endif

        char speciesname_new[128];
        sprintf(speciesname_new, "hydro_%s", speciesname);
        if (tframe >= 1)
        {
            if (tframe == (nframes - 1))
            {
                invert_hydro_xml_item(output_xml_file, speciesname_new, step_for_viou, dimensions_4d, dimensions_3d, 1);
            }
            else
            {
                invert_hydro_xml_item(output_xml_file, speciesname_new, step_for_viou, dimensions_4d, dimensions_3d, 0);
            }
        }
        else
        {
            create_file_with_header(output_xml_file, dimensions_3d, orignal, dxdydz, nframes, fields_interval);
            if (tframe == (nframes - 1))
            {
                invert_hydro_xml_item(output_xml_file, speciesname_new, step_for_viou, dimensions_4d, dimensions_3d, 1);
            }
            else
            {
                invert_hydro_xml_item(output_xml_file, speciesname_new, step_for_viou, dimensions_4d, dimensions_3d, 0);
            }
        }
        tframe++;
    }
}

// TODO": make the sp_name and speciesname varailbe naming consistent
void
vpic_simulation::dump_particles_hdf5( const char *sp_name,
                                 const char *fbase,
                                 int ftag )
{
    size_t step_for_viou = step();
    char fname[256];
    char group_name[256];
    char particle_scratch[128];
    char subparticle_scratch[128];

    int np_local;
    species_t *sp;

    float *Pf;
    int *Pi;

    // get the total number of particles. in this example, output only electrons
    sp = species_list;
    sprintf(particle_scratch, DUMP_DIR_FORMAT, "particle_hdf5");
    dump_mkdir(particle_scratch);
    sprintf(subparticle_scratch, "%s/T.%ld/", particle_scratch, step_for_viou);
    dump_mkdir(subparticle_scratch);

    // TODO: Allow the user to set this

    int stride_particle_dump = 1;
    while (sp)
    {
        np_local = (sp->np + stride_particle_dump - 1) / stride_particle_dump;

        // make a copy of the part of particle data to be dumped
        double ec1 = uptime();

        int sp_np = sp->np;
        int sp_max_np = sp->max_np;
        particle_t *ALIGNED(128) p_buf = NULL;
        if (!p_buf)
            MALLOC_ALIGNED(p_buf, np_local, 128);
        particle_t *sp_p = sp->p;
        sp->p = p_buf;
        sp->np = np_local;
        sp->max_np = np_local;

        for (long long iptl = 0, i = 0; iptl < sp_np; iptl += stride_particle_dump, ++i)
        {
            COPY(&sp->p[i], &sp_p[iptl], 1);
        }
    #ifdef ENABLE_V407_SCIDAC
        # define PBUF_SIZE 32768 // 1MB of particles
        for( int buf_start=0; buf_start<np_local; buf_start += PBUF_SIZE ) {
            int n_buf = PBUF_SIZE;
            if( buf_start+n_buf > np_local ) n_buf = np_local - buf_start;
                COPY( p_buf, &sp->p[buf_start], n_buf );
            center_p( p_buf, n_buf, sp->q_m, interpolator, grid );
        }
    #else
        center_p(sp, interpolator_array);
    #endif
        ec1 = uptime() - ec1;
        int mpi_rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

        //std::cout << "on mpi_rank: " << mpi_rank << ", time in copying particle data: " << ec1 << " s" << ", np_local = " << np_local << std::endl;
        //sim_log("on mpi_rank: " << mpi_rank << ", time in copying particle data: " << ec1 << " s" << ", np_local = " << np_local);

        Pf = (float *)sp->p;
        Pi = (int *)sp->p;

        // open HDF5 file in "particle/T.<step>/" subdirectory
        // filename: eparticle.h5p
        sprintf(fname, "%s/%s_%ld.h5", subparticle_scratch, sp->name, step_for_viou);
        sprintf(group_name, "/Timestep_%ld", step_for_viou);
        double el1 = uptime();

        hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
        hid_t file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
        hid_t group_id = H5Gcreate(file_id, group_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        H5Pclose(plist_id);

        long long total_particles, offset;
        long long numparticles = np_local;
        MPI_Allreduce(&numparticles, &total_particles, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        MPI_Scan(&numparticles, &offset, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        offset -= numparticles;

        hid_t filespace = H5Screate_simple(1, (hsize_t *)&total_particles, NULL);

        hsize_t memspace_count_temp = numparticles * 8;
        hid_t memspace = H5Screate_simple(1, &memspace_count_temp, NULL);
        plist_id = H5Pcreate(H5P_DATASET_XFER);

        //Comment out for test only
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, (hsize_t *)&offset, NULL, (hsize_t *)&numparticles, NULL);

        hsize_t memspace_start = 0, memspace_stride = 8, memspace_count = np_local;
        H5Sselect_hyperslab(memspace, H5S_SELECT_SET, &memspace_start, &memspace_stride, &memspace_count, NULL);

        el1 = uptime() - el1;
        //sim_log("Particle TimeHDF5Open): " << el1 << " s"); //Easy to handle results for scripts

        double el2 = uptime();

        hid_t dset_id = H5Dcreate(group_id, "dX", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        int ierr = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, Pf);
        H5Dclose(dset_id);
        //if (rank == 0) printf ("Written variable dX \n");

        dset_id = H5Dcreate(group_id, "dY", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        ierr = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, Pf + 1);
        H5Dclose(dset_id);
        //if (rank == 0) printf ("Written variable dY \n");

        dset_id = H5Dcreate(group_id, "dZ", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        ierr = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, Pf + 2);
        H5Dclose(dset_id);
        //if (rank == 0) printf ("Written variable dZ \n");

        dset_id = H5Dcreate(group_id, "i", H5T_NATIVE_INT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        ierr = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, Pi + 3);
        H5Dclose(dset_id);
        //if (rank == 0) printf ("Written variable i \n");

        dset_id = H5Dcreate(group_id, "Ux", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        ierr = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, Pf + 4);
        H5Dclose(dset_id);
        //if (rank == 0) printf ("Written variable  Ux \n");

        dset_id = H5Dcreate(group_id, "Uy", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        ierr = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, Pf + 5);
        H5Dclose(dset_id);
        //if (rank == 0) printf ("Written variable Uy \n");

        dset_id = H5Dcreate(group_id, "Uz", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        ierr = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, Pf + 6);
        H5Dclose(dset_id);
        //if (rank == 0) printf ("Written variable Uz \n");

        dset_id = H5Dcreate(group_id, "q", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        ierr = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, Pf + 7);
        H5Dclose(dset_id);
        //if (rank == 0) printf ("Written variable q \n");

        el2 = uptime() - el2;
        //sim_log("Particle TimeHDF5Write: " << el2 << " s");

        double el3 = uptime();
        H5Sclose(memspace);
        H5Sclose(filespace);
        H5Pclose(plist_id);
        H5Gclose(group_id);
        H5Fclose(file_id);
        el3 = uptime() - el3;
        //sim_log("Particle TimeHDF5Close: " << el3 << " s");

        sp->p = sp_p;
        sp->np = sp_np;
        sp->max_np = sp_max_np;
        FREE_ALIGNED(p_buf);

        // Write metadata if step() == 0
        char meta_fname[256];

        sprintf(meta_fname, "%s/grid_metadata_%s_%ld.h5", subparticle_scratch, sp->name, step_for_viou);

        double meta_el1 = uptime();

        hid_t meta_plist_id = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(meta_plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
        hid_t meta_file_id = H5Fcreate(meta_fname, H5F_ACC_TRUNC, H5P_DEFAULT, meta_plist_id);
        hid_t meta_group_id = H5Gcreate(meta_file_id, group_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Pclose(meta_plist_id);

        long long meta_total_particles, meta_offset;
        long long meta_numparticles = 1;
        MPI_Allreduce(&meta_numparticles, &meta_total_particles, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        MPI_Scan(&meta_numparticles, &meta_offset, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        meta_offset -= meta_numparticles;

        hid_t meta_filespace = H5Screate_simple(1, (hsize_t *)&meta_total_particles, NULL);
        hid_t meta_memspace = H5Screate_simple(1, (hsize_t *)&meta_numparticles, NULL);
        meta_plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(meta_plist_id, H5FD_MPIO_COLLECTIVE);
        H5Sselect_hyperslab(meta_filespace, H5S_SELECT_SET, (hsize_t *)&meta_offset, NULL, (hsize_t *)&meta_numparticles, NULL);
        meta_el1 = uptime() - meta_el1;
        //sim_log("Metafile TimeHDF5Open): " << meta_el1 << " s"); //Easy to handle results for scripts

        double meta_el2 = uptime();

        hid_t meta_dset_id = H5Dcreate(meta_group_id, "np_local", H5T_NATIVE_INT, meta_filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        ierr = H5Dwrite(meta_dset_id, H5T_NATIVE_INT, meta_memspace, meta_filespace, meta_plist_id, (int32_t *)&np_local);
        H5Dclose(meta_dset_id);
        //if (rank == 0) printf ("Written variable dX \n");

        meta_dset_id = H5Dcreate(meta_group_id, "nx", H5T_NATIVE_INT, meta_filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        ierr = H5Dwrite(meta_dset_id, H5T_NATIVE_INT, meta_memspace, meta_filespace, meta_plist_id, &grid->nx);
        H5Dclose(meta_dset_id);
        //if (rank == 0) printf ("Written variable dY \n");

        meta_dset_id = H5Dcreate(meta_group_id, "ny", H5T_NATIVE_INT, meta_filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        ierr = H5Dwrite(meta_dset_id, H5T_NATIVE_INT, meta_memspace, meta_filespace, meta_plist_id, &grid->ny);
        H5Dclose(meta_dset_id);
        //if (rank == 0) printf ("Written variable dZ \n");

        meta_dset_id = H5Dcreate(meta_group_id, "nz", H5T_NATIVE_INT, meta_filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        ierr = H5Dwrite(meta_dset_id, H5T_NATIVE_INT, meta_memspace, meta_filespace, meta_plist_id, &grid->nz);
        H5Dclose(meta_dset_id);
        //if (rank == 0) printf ("Written variable i \n");

        meta_dset_id = H5Dcreate(meta_group_id, "x0", H5T_NATIVE_FLOAT, meta_filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        ierr = H5Dwrite(meta_dset_id, H5T_NATIVE_FLOAT, meta_memspace, meta_filespace, meta_plist_id, &grid->x0);
        H5Dclose(meta_dset_id);

        meta_dset_id = H5Dcreate(meta_group_id, "y0", H5T_NATIVE_FLOAT, meta_filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        ierr = H5Dwrite(meta_dset_id, H5T_NATIVE_FLOAT, meta_memspace, meta_filespace, meta_plist_id, &grid->y0);
        H5Dclose(meta_dset_id);

        meta_dset_id = H5Dcreate(meta_group_id, "z0", H5T_NATIVE_FLOAT, meta_filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        ierr = H5Dwrite(meta_dset_id, H5T_NATIVE_FLOAT, meta_memspace, meta_filespace, meta_plist_id, &grid->z0);
        H5Dclose(meta_dset_id);

        meta_dset_id = H5Dcreate(meta_group_id, "dx", H5T_NATIVE_FLOAT, meta_filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        ierr = H5Dwrite(meta_dset_id, H5T_NATIVE_FLOAT, meta_memspace, meta_filespace, meta_plist_id, &grid->dx);
        H5Dclose(meta_dset_id);

        meta_dset_id = H5Dcreate(meta_group_id, "dy", H5T_NATIVE_FLOAT, meta_filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        ierr = H5Dwrite(meta_dset_id, H5T_NATIVE_FLOAT, meta_memspace, meta_filespace, meta_plist_id, &grid->dy);
        H5Dclose(meta_dset_id);

        meta_dset_id = H5Dcreate(meta_group_id, "dz", H5T_NATIVE_FLOAT, meta_filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        ierr = H5Dwrite(meta_dset_id, H5T_NATIVE_FLOAT, meta_memspace, meta_filespace, meta_plist_id, &grid->dz);
        H5Dclose(meta_dset_id);

        meta_el2 = uptime() - meta_el2;
        //sim_log("Metafile TimeHDF5Write: " << meta_el2 << " s");
        double meta_el3 = uptime();
        H5Sclose(meta_memspace);
        H5Sclose(meta_filespace);
        H5Pclose(meta_plist_id);
        H5Gclose(meta_group_id);
        H5Fclose(meta_file_id);
        meta_el3 = uptime() - meta_el3;
        //sim_log("Metafile TimeHDF5Close: " << meta_el3 << " s");

        sp = sp->next;
    }
}
#endif

void
vpic_simulation::dump_particles( const char *sp_name,
                                 const char *fbase,
                                 int ftag )
{
  species_t *sp;
  char fname[256];
  FileIO fileIO;
  int dim[1], buf_start;
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

  WRITE_HEADER_V0( dump_type::particle_dump, sp->id, sp->q / sp->m, step(), fileIO );

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
vpic_simulation::field_dump( DumpParameters & dumpParams,
			     field_t *f,
                             int64_t userStep )
{
  long dumpStep = ( userStep == -1 ) ? (long) step() : userStep;

  // Create directory for this time step
  char timeDir[256];

  sprintf( timeDir,
           "%s/T.%ld",
           dumpParams.baseDir,
           dumpStep );

  dump_mkdir( timeDir );

  // Open the file for output
  char filename[256];

  sprintf( filename,
           "%s/T.%ld/%s.%ld.%d",
           dumpParams.baseDir,
           dumpStep,
           dumpParams.baseFileName,
           dumpStep,
           rank() );

  FileIO fileIO;
  FileIOStatus status;

  status = fileIO.open(filename, io_write);
  if( status==fail ) ERROR(( "Failed opening file: %s", filename ));

  // default is to write field_array->f
  if ( f==NULL ) f = field_array->f;

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

  if ( dumpParams.format == band )
  {
    WRITE_HEADER_V0( dump_type::field_dump, -1, 0, dumpStep, fileIO );

    dim[0] = nxout+2;
    dim[1] = nyout+2;
    dim[2] = nzout+2;

    if ( rank() == VERBOSE_rank )
    {
      std::cerr << "nxout: " << nxout << std::endl;
      std::cerr << "nyout: " << nyout << std::endl;
      std::cerr << "nzout: " << nzout << std::endl;
      std::cerr << "nx: " << grid->nx << std::endl;
      std::cerr << "ny: " << grid->ny << std::endl;
      std::cerr << "nz: " << grid->nz << std::endl;
    }

    WRITE_ARRAY_HEADER(f, 3, dim, fileIO);

    // Create a variable list of field values to output.
    size_t numvars = std::min(dumpParams.output_vars.bitsum(),
                              total_field_variables);
    size_t * varlist = new size_t[numvars];

    for(size_t i(0), c(0); i<total_field_variables; i++)
      if(dumpParams.output_vars.bitset(i)) varlist[c++] = i;

    if( rank()==VERBOSE_rank ) printf("\nBEGIN_OUTPUT\n");

    // more efficient for standard case
    if ( istride == 1 &&
	 jstride == 1 &&
	 kstride == 1 )
    {
      for(size_t v(0); v<numvars; v++) {
      for(size_t k(0); k<nzout+2; k++) {
      for(size_t j(0); j<nyout+2; j++) {
      for(size_t i(0); i<nxout+2; i++) {
              const uint32_t * fref = reinterpret_cast<uint32_t *>(&f(i,j,k));
              fileIO.write(&fref[varlist[v]], 1);
              if(rank()==VERBOSE_rank) printf("%f ", f(i,j,k).ex);
              if(rank()==VERBOSE_rank) std::cout << "(" << i << " " << j << " " << k << ")" << std::endl;
      } if(rank()==VERBOSE_rank) std::cout << std::endl << "ROW_BREAK " << j << " " << k << std::endl;
      } if(rank()==VERBOSE_rank) std::cout << std::endl << "PLANE_BREAK " << k << std::endl;
      } if(rank()==VERBOSE_rank) std::cout << std::endl << "BLOCK_BREAK" << std::endl;
      }
    }

    else
    {
      for(size_t v(0); v<numvars; v++) {
      for(size_t k(0); k<nzout+2; k++) { const size_t koff = (k == 0) ? 0 : (k == nzout+1) ? grid->nz+1 : k*kstride-1;
      for(size_t j(0); j<nyout+2; j++) { const size_t joff = (j == 0) ? 0 : (j == nyout+1) ? grid->ny+1 : j*jstride-1;
      for(size_t i(0); i<nxout+2; i++) { const size_t ioff = (i == 0) ? 0 : (i == nxout+1) ? grid->nx+1 : i*istride-1;
              const uint32_t * fref = reinterpret_cast<uint32_t *>(&f(ioff,joff,koff));
              fileIO.write(&fref[varlist[v]], 1);
              if(rank()==VERBOSE_rank) printf("%f ", f(ioff,joff,koff).ex);
              if(rank()==VERBOSE_rank) std::cout << "(" << ioff << " " << joff << " " << koff << ")" << std::endl;
      } if(rank()==VERBOSE_rank) std::cout << std::endl << "ROW_BREAK " << joff << " " << koff << std::endl;
      } if(rank()==VERBOSE_rank) std::cout << std::endl << "PLANE_BREAK " << koff << std::endl;
      } if(rank()==VERBOSE_rank) std::cout << std::endl << "BLOCK_BREAK" << std::endl;
      }
    }

    delete[] varlist;

  }

  // band_interleave
  else
  {
    WRITE_HEADER_V0( dump_type::field_dump, -1, 0, dumpStep, fileIO );

    dim[0] = nxout+2;
    dim[1] = nyout+2;
    dim[2] = nzout+2;

    WRITE_ARRAY_HEADER(f, 3, dim, fileIO);

    if ( istride == 1 &&
	 jstride == 1 &&
	 kstride == 1 )
    {
      fileIO.write( f, dim[0] * dim[1] * dim[2] );
    }

    else
    {
      for(size_t k(0); k<nzout+2; k++) { const size_t koff = (k == 0) ? 0 : (k == nzout+1) ? grid->nz+1 : k*kstride-1;
      for(size_t j(0); j<nyout+2; j++) { const size_t joff = (j == 0) ? 0 : (j == nyout+1) ? grid->ny+1 : j*jstride-1;
      for(size_t i(0); i<nxout+2; i++) { const size_t ioff = (i == 0) ? 0 : (i == nxout+1) ? grid->nx+1 : i*istride-1;
            fileIO.write( &f( ioff, joff, koff ), 1 );
      }
      }
      }
    }
  }

# undef f

  if ( fileIO.close() ) ERROR(( "File close failed on field dump." ));
}

void
vpic_simulation::hydro_dump( const char * speciesname,
                             DumpParameters & dumpParams,
                             hydro_t *h,
                             int64_t userStep )
{
  long dumpStep = ( userStep == -1 ) ? (long) step() : userStep;

  // Create directory for this time step
  char timeDir[256];

  sprintf( timeDir,
           "%s/T.%ld",
           dumpParams.baseDir,
           dumpStep );

  dump_mkdir( timeDir );

  // Open the file for output
  char filename[256];

  sprintf( filename,
           "%s/T.%ld/%s.%ld.%d",
           dumpParams.baseDir,
           dumpStep,
           dumpParams.baseFileName,
           dumpStep,
           rank() );

  FileIO fileIO;
  FileIOStatus status;

  status = fileIO.open(filename, io_write);
  if(status == fail) ERROR(("Failed opening file: %s", filename));

  species_t * sp = find_species_name(speciesname, species_list);
  if( !sp ) ERROR(( "Invalid species name: %s", speciesname ));

  // default behavior is to accumulate hydro array and then write
  if ( h == NULL )
  {
    h = hydro_array->h;
    clear_hydro_array( hydro_array );
    accumulate_hydro_p( hydro_array, sp, interpolator_array );
    synchronize_hydro_array( hydro_array );
  }

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
# define hydro(x,y,z) h[VOXEL(x,y,z, grid->nx,grid->ny,grid->nz)]

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
  if ( dumpParams.format == band )
  {
    WRITE_HEADER_V0( dump_type::hydro_dump, sp->id, sp->q/sp->m, dumpStep, fileIO );

    dim[0] = nxout+2;
    dim[1] = nyout+2;
    dim[2] = nzout+2;

    WRITE_ARRAY_HEADER(h, 3, dim, fileIO);

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

  }

  // band_interleave
  else
  {
    WRITE_HEADER_V0( dump_type::hydro_dump, sp->id, sp->q/sp->m, dumpStep, fileIO );

    dim[0] = nxout;
    dim[1] = nyout;
    dim[2] = nzout;

    WRITE_ARRAY_HEADER(h, 3, dim, fileIO);

    if ( istride == 1 &&
	 jstride == 1 &&
	 kstride == 1 )
    {
      fileIO.write( h, dim[0] * dim[1] * dim[2] );
    }

    else
    {
      for(size_t k(0); k<nzout; k++) { const size_t koff = (k == 0) ? 0 : (k == nzout+1) ? grid->nz+1 : k*kstride-1;
      for(size_t j(0); j<nyout; j++) { const size_t joff = (j == 0) ? 0 : (j == nyout+1) ? grid->ny+1 : j*jstride-1;
      for(size_t i(0); i<nxout; i++) { const size_t ioff = (i == 0) ? 0 : (i == nxout+1) ? grid->nx+1 : i*istride-1;
            fileIO.write( &hydro( ioff, joff, koff ), 1 );
      }
      }
      }
    }
  }

# undef hydro

  if( fileIO.close() ) ERROR(( "File close failed on hydro dump!!!" ));
}
