#ifndef Dump_Strategy_h
#define Dump_Strategy_h

#include <unordered_map>
#include <vector>
#include <iostream>

//#define DUMP_INFO_DEBUG 1
//#define H5_ASYNC 1
#ifdef H5_ASYNC
#include "h5_vol_external_async_native.h"
#endif
//#define CHUNK_FLAG 1

//#define METADATA_COLL_WRITE 1
//#define TRUE 1

//#define HAS_FIELD_COMP 1
//#define HAS_HYDRO_COMP 1
//#define HAS_PARTICLE_COMP 1

//#define HAS_INDEPENDENT_IO 1

#include <mpi.h> // TODO: it would be good if this didn't have to know about MPI
#include <cassert>

// TODO: should I drop the ./src here?
#include "../util/io/FileIO.h"
#include "../util/util_base.h"
#include "../util/io/FileUtils.h"
#include "../field_advance/field_advance.h"
#include "../sf_interface/sf_interface.h"
#include "../species_advance/species_advance.h"

#include "dump.h"
#include "dumpmacros.h"

#ifdef VPIC_ENABLE_HDF5
#include "hdf5.h" // from the lib
#include "hdf5_header_info.h" // from vpic
//#define daos "daos:"
#define daos ""
#endif

#ifdef VPIC_ENABLE_OPENPMD
#include <openPMD/openPMD.hpp>
#endif

//#define N_FILE_N_PROCESS 1
//#define TEST_MPIIO 1

// TODO: delete this
#define _LOG_PREFIX \
  __FILE__ "(" EXPAND_AND_STRINGIFY(__LINE__) ")[" << rank << "]: "
#define io_log(x) do {                                \
    if( rank==0 ) {                                  \
      std::cerr << _LOG_PREFIX << x << std::endl;  \
      std::cerr.flush();                               \
    }                                                  \
  } while(0)


// Runtime inheritance is obviously not very "VPIC like", as we will [probably]
// incur a penalty for the vtable lookup, but given we're about to do IO this
// is very negligible.
class Dump_Strategy {
    public:
    int rank, nproc, num_step;

    Dump_Strategy(int _rank, int _nproc ) :
        rank(_rank),
        nproc(_nproc)
    { } // empty

    virtual ~Dump_Strategy() { };

    virtual void dump_fields(
        const char *fbase,
        int step,
        grid_t* grid,
        field_array_t* field_array,
        int ftag
    ) = 0;
    virtual void dump_hydro(
        const char *fbase,
        int step,
        hydro_array_t* hydro_array,
        species_t* sp,
        interpolator_array_t* interpolator_array,
        grid_t* grid,
        int ftag
    ) = 0;
    virtual void dump_particles(
        const char *fbase,
        species_t* sp,
        grid_t* grid,
        int step,
        interpolator_array_t* interpolator_array,
        int ftag
    ) = 0;
};

class BinaryDump : public Dump_Strategy {
    public:
        using Dump_Strategy::Dump_Strategy; // inherit constructor
        //BinaryDump(int _rank, int _nproc ) : Dump_Strategy(_rank, _nproc ){ } // empty

        // TODO: now we pass rank and step, ftag has odd semanticds
        void dump_fields(
                const char *fbase,
                int step,
                grid_t* grid,
                field_array_t* field_array,
                int ftag
        );
        void dump_hydro(
                const char *fbase,
                int step,
                hydro_array_t* hydro_array,
                species_t* sp,
                interpolator_array_t* interpolator_array,
                grid_t* grid,
                int ftag
        );
        void dump_particles(
                const char *fbase,
                species_t* sp,
                grid_t* grid,
                int step,
                interpolator_array_t* interpolator_array,
                int ftag
        );
};

#ifdef VPIC_ENABLE_HDF5

struct field_dump_flag_t
{
  bool ex = true, ey = true, ez = true, div_e_err = true;
  bool cbx = true, cby = true, cbz = true, div_b_err = true;
  bool tcax = true, tcay = true, tcaz = true, rhob = true;
  bool jfx = true, jfy = true, jfz = true, rhof = true;
  bool ematx = true, ematy = true, ematz = true, nmat = true;
  bool fmatx = true, fmaty = true, fmatz = true, cmat = true;
  void disableE()
  {
    ex = false, ey = false, ez = false, div_e_err = false;
  }

  void disableCB()
  {
    cbx = false, cby = false, cbz = false, div_b_err = false;
  }

  void disableTCA()
  {
    tcax = false, tcay = false, tcaz = false, rhob = false;
  }

  void disableJF()
  {
    jfx = false, jfy = false, jfz = false, rhof = false;
  }

  void disableEMAT()
  {
    ematx = false, ematy = false, ematz = false, nmat = false;
  }

  void disableFMAT()
  {
    fmatx = false, fmaty = false, fmatz = false, cmat = false;
  }

  void resetToDefaults()
  {
    ex = true, ey = true, ez = true, div_e_err = true;
    cbx = true, cby = true, cbz = true, div_b_err = true;
    tcax = true, tcay = true, tcaz = true, rhob = true;
    jfx = true, jfy = true, jfz = true, rhof = true;
    ematx = true, ematy = true, ematz = true, nmat = true;
    fmatx = true, fmaty = true, fmatz = true, cmat = true;
  }

  bool enabledE()
  {
    return ex && ey && ez;
  }

  bool enabledCB()
  {
    return cbx && cby && cbz;
  }

  bool enabledTCA()
  {
    return tcax && tcay && tcaz;
  }

  bool enabledJF()
  {
    return jfx && jfy && jfz;
  }

  bool enabledEMAT()
  {
    return ematx && ematy && ematz;
  }

  bool enabledFMAT()
  {
    return fmatx && fmaty && fmatz;
  }
};

struct hydro_dump_flag_t
{
  bool jx = true, jy = true, jz = true, rho = true;
  bool px = true, py = true, pz = true, ke = true;
  bool txx = true, tyy = true, tzz = true;
  bool tyz = true, tzx = true, txy = true;

  void disableJ()
  {
    jx = false, jy = false, jz = false, rho = false;
  }

  void disableP()
  {
    px = false, py = false, pz = false, ke = false;
  }

  void disableTD() //Stress diagonal
  {
    txx = false, tyy = false, tzz = false;
  }

  void disableTOD() //Stress off-diagonal
  {
    tyz = false, tzx = false, txy = false;
  }
  void resetToDefaults()
  {
    jx = true, jy = true, jz = true, rho = true;
    px = true, py = true, pz = true, ke = true;
    txx = true, tyy = true, tzz = true;
    tyz = true, tzx = true, txy = true;
  }

  bool enabledJ()
  {
    return jx && jy && jz;
  }

  bool enabledP()
  {
    return px && py && pz;
  }

  bool enabledTD()
  {
    return txx && tyy && tzz;
  }

  bool enabledTOD()
  {
    return tyz && tzx && txy;
  }
};
class HDF5Dump : public Dump_Strategy {
    std::unordered_map<species_id, size_t> tframe_map;
    public:
        using Dump_Strategy::Dump_Strategy; // inherit constructor

        // TODO: replace these with a common dump interface
        // Declare vars to use
        hydro_dump_flag_t hydro_dump_flag;
        field_dump_flag_t field_dump_flag;

#define DUMP_DIR_FORMAT "./%s"

// TODO: naming a macro so close to existing functions AND data is not a good
// define to do C-style indexing
#define _hydro(x, y, z) hydro_array->h[VOXEL(x, y, z, grid->nx, grid->ny, grid->nz)]

        // TODO: make function?
        void dump_fields(
            const char *fbase,
            int step,
            grid_t* grid,
            field_array_t* field_array,
            int ftag
        )
        {
            size_t step_for_viou = step;
            int file_exist_flag = 0;

            int mpi_size, mpi_rank;
            MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
            MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

#ifdef DUMP_INFO_DEBUG
            printf("MPI rank = %d, size = %d \n", mpi_rank, mpi_size);
            //printf("base dir for field: %s \n", fdParams.baseDir);
            //printf("stride x y z  = (%ld, %ld, %ld)\n", fdParams.stride_x, fdParams.stride_y, fdParams.stride_z);
            printf("grid x, y z  = (%d, %d, %d) \n", grid->nx, grid->ny, grid->nz);
            printf("domain loc (x0, y0, z0) -> (x1, y1, z1) = (%f, %f, %f) -> (%f, %f, %f) \n", grid->x0, grid->y0, grid->z0, grid->x1, grid->y1, grid->z1);
            //printf("global->topology_x, y, z =  %f, %f, %f \n ", global->topology_x, global->topology_y, global->topology_z);
            printf("grid -> sx, sy, sz =  (%d, %d, %d), nv=%d \n", grid->sx, grid->sy, grid->sz, grid->nv);
#endif
            char fname[256];
            char field_scratch[128];
            char subfield_scratch[128];

            sprintf(field_scratch, DUMP_DIR_FORMAT, "field_hdf5");
            FileUtils::makeDirectory(field_scratch);
            sprintf(subfield_scratch, "%s/T.%zu/", field_scratch, step_for_viou);
            FileUtils::makeDirectory(subfield_scratch);

            sprintf(fname, "%s%s/%s_%zu.h5", daos, subfield_scratch, "fields", step_for_viou);
            double el1 = uptime();

//    int file_exist(const char *filename)
//{
//    struct stat buffer;
//    return (stat(filename, &buffer) == 0);
//}

            struct stat buffer;
            if((stat(fname, &buffer) == 0)){
                file_exist_flag  = 1;
                if(!mpi_rank)
                    printf("Write original files /w HDF5! \n");
            }
            file_exist_flag = 0;

            hid_t plist_id;
            hid_t file_id;
            if(!file_exist_flag){
                plist_id = H5Pcreate(H5P_FILE_ACCESS);
                H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
                //H5Pset_alignment(plist_id, 4194304, 4194304);
                /*if(H5Pset_libver_bounds(plist_id, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST) < 0){
                    exit(-1);
                }*/

                #ifdef METADATA_COLL_WRITE
                if(!mpi_rank) printf("Enable collective metadata write !\n");
                H5Pset_coll_metadata_write(plist_id, TRUE);
                #endif
                file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
                H5Pclose(plist_id);
            }else{
                file_id = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT);
            }

            sprintf(fname, "Timestep_%zu", step_for_viou);
            hid_t group_id;
            if(!file_exist_flag){
                group_id = H5Gcreate(file_id, fname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            }else{
                group_id = H5Gopen(file_id, fname,  H5P_DEFAULT);
            }
            el1 = uptime() - el1;
            io_log("TimeHDF5Open): " << el1 << " s"); //Easy to handle results for scripts
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

            hsize_t temp_buf_index;
            hid_t dset_id;
            //char  *field_var_name[] = {"ex","ey","ez","div_e_err","cbx","cby","cbz","div_b_err","tcax","tcay","tcaz","rhob","jfx","jfy","jfz","rhof"};
            //Comment out for test only

            if(!file_exist_flag){
                plist_id = H5Pcreate(H5P_DATASET_XFER);
                //H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
                #ifdef HAS_INDEPENDENT_IO
                    if(!mpi_rank) printf("\n ###\n VPIC Independent I/O! \n ###\n");
                    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
                #else
                    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
                #endif
            }else{
                plist_id = H5P_DEFAULT;

            }

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

            // Convert rank to local decomposition
            int rx, ry, rz;
            UNVOXEL(mpi_rank, rx, ry, rz, grid->gpx, grid->gpy, grid->gpz);

            int mpi_rank_x, mpi_rank_y, mpi_rank_z;
            mpi_rank_x = rx;
            mpi_rank_y = ry;
            mpi_rank_z = rz;

            global_offset[0] = (grid->nx) * mpi_rank_x;
            global_offset[1] = (grid->ny) * mpi_rank_y;
            global_offset[2] = (grid->nz) * mpi_rank_z;

            global_count[0] = (grid->nx);
            global_count[1] = (grid->ny);
            global_count[2] = (grid->nz);

#ifdef DUMP_INFO_DEBUG
            if(mpi_rank < 4){
            printf("grid nx, ny nz  = (%d, %d, %d) \n", grid->nx, grid->ny, grid->nz);
            printf("mpi-rank = %d, rank index = (%d, %d, %d) \n", mpi_rank, mpi_rank_x, mpi_rank_y, mpi_rank_z);
            printf("global size   = %llu  %llu %llu \n", field_global_size[0], field_global_size[1], field_global_size[2]);
            printf("global_offset = %llu %llu %llu \n", global_offset[0], global_offset[1], global_offset[2]);
            printf("global_count  = %llu  %llu %llu \n", global_count[0], global_count[1], global_count[2]);
            fflush(stdout);
            }
#endif
            hid_t filespace;   //= H5Screate_simple(3, field_global_size, NULL);
            hid_t memspace; // = H5Screate_simple(3, field_local_size, NULL);
            if(!file_exist_flag){
                filespace = H5Screate_simple(3, field_global_size, NULL);
            }
            memspace = H5Screate_simple(3, field_local_size, NULL);

            hsize_t chunk_dims[3];
            chunk_dims[0] = 288; //grid->nx;  //8 x 8 x 8
            chunk_dims[1] = 24;  //grid->ny; //
            chunk_dims[2] = 24;  //grid->nz;

            hid_t dataspace_id;
            hid_t dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
#ifdef CHUNK_FLAG
            H5Pset_chunk(dcpl_id, 3, chunk_dims);
            if(!mpi_rank)   printf("Enable chunking !\n");
#endif

            /*
               typedef struct field {
               float ex,   ey,   ez,   div_e_err;     // Electric field and div E error
               float cbx,  cby,  cbz,  div_b_err;     // Magnetic field and div B error
               float tcax, tcay, tcaz, rhob;          // TCA fields and bound charge density
               float jfx,  jfy,  jfz,  rhof;          // Free current and charge density
               material_id ematx, ematy, ematz, nmat; // Material at edge centers and nodes
               material_id fmatx, fmaty, fmatz, cmat; // Material at face and cell centers
               } field_t;*/


#if HAS_FIELD_COMP

            if(!mpi_rank)
              printf("Using Field Compund type !\n");
            hid_t  field_comp_type_it = H5Tcreate (H5T_COMPOUND, sizeof(field_t));
            H5Tinsert(field_comp_type_it, "ex", HOFFSET(field_t, ex), H5T_NATIVE_FLOAT);
            H5Tinsert(field_comp_type_it, "ey", HOFFSET(field_t, ey), H5T_NATIVE_FLOAT);
            H5Tinsert(field_comp_type_it, "ez", HOFFSET(field_t, ez), H5T_NATIVE_FLOAT);
            H5Tinsert(field_comp_type_it, "div_e_err", HOFFSET(field_t, div_e_err), H5T_NATIVE_FLOAT);

            H5Tinsert(field_comp_type_it, "cbx", HOFFSET(field_t, cbx), H5T_NATIVE_FLOAT);
            H5Tinsert(field_comp_type_it, "cby", HOFFSET(field_t, cby), H5T_NATIVE_FLOAT);
            H5Tinsert(field_comp_type_it, "cbz", HOFFSET(field_t, cbz), H5T_NATIVE_FLOAT);
            H5Tinsert(field_comp_type_it, "div_b_err", HOFFSET(field_t, div_b_err), H5T_NATIVE_FLOAT);

            H5Tinsert(field_comp_type_it, "tcax", HOFFSET(field_t, tcax), H5T_NATIVE_FLOAT);
            H5Tinsert(field_comp_type_it, "tcay", HOFFSET(field_t, tcay), H5T_NATIVE_FLOAT);
            H5Tinsert(field_comp_type_it, "tcaz", HOFFSET(field_t, tcaz), H5T_NATIVE_FLOAT);
            H5Tinsert(field_comp_type_it, "rhob", HOFFSET(field_t, rhob), H5T_NATIVE_FLOAT);

            H5Tinsert(field_comp_type_it, "jfx", HOFFSET(field_t, jfx), H5T_NATIVE_FLOAT);
            H5Tinsert(field_comp_type_it, "jfy", HOFFSET(field_t, jfy), H5T_NATIVE_FLOAT);
            H5Tinsert(field_comp_type_it, "jfz", HOFFSET(field_t, jfz), H5T_NATIVE_FLOAT);
            H5Tinsert(field_comp_type_it, "rhof", HOFFSET(field_t, rhof), H5T_NATIVE_FLOAT);

            H5Tinsert(field_comp_type_it, "ematx", HOFFSET(field_t, ematx), H5T_NATIVE_SHORT);
            H5Tinsert(field_comp_type_it, "ematy", HOFFSET(field_t, ematy), H5T_NATIVE_SHORT);
            H5Tinsert(field_comp_type_it, "ematz", HOFFSET(field_t, ematz), H5T_NATIVE_SHORT);
            H5Tinsert(field_comp_type_it, "nmat", HOFFSET(field_t, nmat), H5T_NATIVE_SHORT);

            H5Tinsert(field_comp_type_it, "fmatx", HOFFSET(field_t, fmatx), H5T_NATIVE_SHORT);
            H5Tinsert(field_comp_type_it, "fmaty", HOFFSET(field_t, fmaty), H5T_NATIVE_SHORT);
            H5Tinsert(field_comp_type_it, "fmatz", HOFFSET(field_t, fmatz), H5T_NATIVE_SHORT);
            H5Tinsert(field_comp_type_it, "cmat", HOFFSET(field_t, cmat), H5T_NATIVE_SHORT);

            if(!file_exist_flag){
              dset_id = H5Dcreate(group_id, "FIELD_DATA_COMP_T", field_comp_type_it, filespace, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
            }else{
              dset_id = H5Dopen(group_id, "FIELD_DATA_COMP_T", H5P_DEFAULT);
            }
            dataspace_id = H5Dget_space(dset_id);
            H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, global_offset, NULL, global_count, NULL);

            field_t *temp_buf = (field_t *)malloc(sizeof(field_t) * (grid->nx) * (grid->ny) * (grid->nz));
            temp_buf_index = 0;
            for (size_t i(1); i < grid->nx + 1; i++)
              {
                for (size_t j(1); j < grid->ny + 1; j++)
                  {
                    for (size_t k(1); k < grid->nz + 1; k++)
                      {
                        temp_buf[temp_buf_index].ex = FIELD_ARRAY_NAME->fpp(i, j, k).ex;
                        temp_buf[temp_buf_index].ey = FIELD_ARRAY_NAME->fpp(i, j, k).ey;
                        temp_buf[temp_buf_index].ez = FIELD_ARRAY_NAME->fpp(i, j, k).ez;
                        temp_buf[temp_buf_index].div_e_err = FIELD_ARRAY_NAME->fpp(i, j, k).div_e_err;
                        temp_buf[temp_buf_index].cbx = FIELD_ARRAY_NAME->fpp(i, j, k).cbx;
                        temp_buf[temp_buf_index].cby = FIELD_ARRAY_NAME->fpp(i, j, k).cby;
                        temp_buf[temp_buf_index].cbz = FIELD_ARRAY_NAME->fpp(i, j, k).cbz;
                        temp_buf[temp_buf_index].div_b_err = FIELD_ARRAY_NAME->fpp(i, j, k).div_b_err;
                        temp_buf[temp_buf_index].tcax = FIELD_ARRAY_NAME->fpp(i, j, k).tcax;
                        temp_buf[temp_buf_index].tcay = FIELD_ARRAY_NAME->fpp(i, j, k).tcay;
                        temp_buf[temp_buf_index].tcaz = FIELD_ARRAY_NAME->fpp(i, j, k).tcaz;
                        temp_buf[temp_buf_index].rhob = FIELD_ARRAY_NAME->fpp(i, j, k).rhob;
                        temp_buf[temp_buf_index].jfx = FIELD_ARRAY_NAME->fpp(i, j, k).jfx;
                        temp_buf[temp_buf_index].jfy = FIELD_ARRAY_NAME->fpp(i, j, k).jfy;
                        temp_buf[temp_buf_index].jfz = FIELD_ARRAY_NAME->fpp(i, j, k).jfz;
                        temp_buf[temp_buf_index].rhof = FIELD_ARRAY_NAME->fpp(i, j, k).rhof;
                        temp_buf[temp_buf_index].ematx = FIELD_ARRAY_NAME->fpp(i, j, k).ematx;
                        temp_buf[temp_buf_index].ematy = FIELD_ARRAY_NAME->fpp(i, j, k).ematy;
                        temp_buf[temp_buf_index].ematz = FIELD_ARRAY_NAME->fpp(i, j, k).ematz;
                        temp_buf[temp_buf_index].nmat = FIELD_ARRAY_NAME->fpp(i, j, k).nmat;
                        temp_buf[temp_buf_index].fmatx = FIELD_ARRAY_NAME->fpp(i, j, k).fmatx;
                        temp_buf[temp_buf_index].fmaty = FIELD_ARRAY_NAME->fpp(i, j, k).fmaty;
                        temp_buf[temp_buf_index].fmatz = FIELD_ARRAY_NAME->fpp(i, j, k).fmatz;
                        temp_buf[temp_buf_index].cmat = FIELD_ARRAY_NAME->fpp(i, j, k).cmat;
                        temp_buf_index = temp_buf_index + 1;
                      }
                  }
              }

            H5Dwrite(dset_id, field_comp_type_it, memspace, dataspace_id, plist_id, temp_buf);
            H5Sclose(dataspace_id);
            H5Dclose(dset_id);
            H5Tclose(field_comp_type_it);
            free(temp_buf);
#else

#define DUMP_FIELD_TO_HDF5(DSET_NAME, ATTRIBUTE_NAME, ELEMENT_TYPE)     \
            {                                                           \
              if(!file_exist_flag){                                     \
                dset_id = H5Dcreate(group_id, DSET_NAME, ELEMENT_TYPE, filespace, H5P_DEFAULT, dcpl_id, H5P_DEFAULT); \
              }else{                                                    \
                dset_id = H5Dopen(group_id, DSET_NAME, H5P_DEFAULT);    \
                dataspace_id = H5Dget_space(dset_id);                   \
              }                                                         \
              temp_buf_index = 0;                                       \
              for (size_t i(1); i < grid->nx + 1; i++)                  \
                {                                                       \
                  for (size_t j(1); j < grid->ny + 1; j++)              \
                    {                                                   \
                      for (size_t k(1); k < grid->nz + 1; k++)          \
                        {                                               \
                          temp_buf[temp_buf_index] = FIELD_ARRAY_NAME->fpp(i, j, k).ATTRIBUTE_NAME; \
                          temp_buf_index = temp_buf_index + 1;          \
        }                                                               \
                    }                                                   \
                }                                                       \
              dataspace_id = H5Dget_space(dset_id);                     \
              H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, global_offset, NULL, global_count, NULL); \
              H5Dwrite(dset_id, ELEMENT_TYPE, memspace, dataspace_id, plist_id, temp_buf); \
              H5Sclose(dataspace_id);                                   \
              H5Dclose(dset_id);                                        \
        }

            float *temp_buf = (float *)malloc(sizeof(float) * (grid->nx) * (grid->ny) * (grid->nz));

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

            free(temp_buf);
#endif

            H5D_mpio_actual_io_mode_t actual_io_mode;
            H5Pget_mpio_actual_io_mode(plist_id,  &actual_io_mode);
            switch(actual_io_mode){
            case H5D_MPIO_NO_COLLECTIVE:
              io_log("H5Pget_mpio_actual_io_mode: H5D_MPIO_NO_COLLECTIVE: ");
              break;
            case H5D_MPIO_CHUNK_INDEPENDENT:
              io_log("H5Pget_mpio_actual_io_mode: H5D_MPIO_CHUNK_INDEPENDENT: ");
              break;
            case H5D_MPIO_CHUNK_COLLECTIVE:
              io_log("H5Pget_mpio_actual_io_mode: H5D_MPIO_CHUNK_COLLECTIVE: ");
              break;
            case H5D_MPIO_CHUNK_MIXED:
              io_log("H5Pget_mpio_actual_io_mode: H5D_MPIO_CHUNK_MIXED: ");
              break;
            case H5D_MPIO_CONTIGUOUS_COLLECTIVE:
              io_log("H5Pget_mpio_actual_io_mode: H5D_MPIO_CONTIGUOUS_COLLECTIVE: ");
              break;
            default :
              io_log("H5Pget_mpio_actual_io_mode: None returend: ");
              break;
            }

            H5D_mpio_actual_chunk_opt_mode_t actual_chunk_opt_mode;
            H5Pget_mpio_actual_chunk_opt_mode(plist_id, &actual_chunk_opt_mode);
             switch(actual_chunk_opt_mode){
                case H5D_MPIO_NO_CHUNK_OPTIMIZATION:
                    io_log("H5Pget_mpio_actual_chunk_opt_mode: H5D_MPIO_NO_CHUNK_OPTIMIZATION: ");
                    break;
                case H5D_MPIO_MULTI_CHUNK:
                    io_log("H5Pget_mpio_actual_chunk_opt_mode: H5D_MPIO_MULTI_CHUNK: ");
                    break;
              //  case H5D_MPIO_MULTI_CHUNK_NO_OPT:
              //      io_log("H5Pget_mpio_actual_chunk_opt_mode: H5D_MPIO_MULTI_CHUNK_NO_OPT: ");
               //     break;
                case H5D_MPIO_LINK_CHUNK:
                    io_log("H5Pget_mpio_actual_chunk_opt_mode: H5D_MPIO_LINK_CHUNK: ");
                    break;
                default :
                    io_log("H5Pget_mpio_actual_chunk_opt_mode: None returend: ");
                    break;
            }

            uint32_t local_no_collective_cause,  global_no_collective_cause;
            H5Pget_mpio_no_collective_cause(plist_id, &local_no_collective_cause, &global_no_collective_cause);

            switch(local_no_collective_cause){
                case H5D_MPIO_COLLECTIVE:
                    io_log("local_no_collective_cause: H5D_MPIO_COLLECTIVE: ");
                    break;
                case H5D_MPIO_SET_INDEPENDENT:
                    io_log("local_no_collective_cause: H5D_MPIO_SET_INDEPENDENT: ");
                    break;
                case H5D_MPIO_DATA_TRANSFORMS:
                    io_log("local_no_collective_cause: H5D_MPIO_DATA_TRANSFORMS: ");
                    break;
                //case H5D_MPIO_SET_MPIPOSIX:
                //    io_log("local_no_collective_cause: H5D_MPIO_SET_MPIPOSIX: ");
                //    break;
                case H5D_MPIO_NOT_SIMPLE_OR_SCALAR_DATASPACES:
                    io_log("local_no_collective_cause: H5D_MPIO_NOT_SIMPLE_OR_SCALAR_DATASPACES: ");
                    break;
                //case H5D_MPIO_POINT_SELECTIONS:
                //    io_log("local_no_collective_cause: H5D_MPIO_POINT_SELECTIONS: ");
                //    break;
                // case H5D_MPIO_FILTERS:
                //    io_log("local_no_collective_cause: H5D_MPIO_FILTERS: ");
                //    break;
                default :
                    io_log("local_no_collective_cause: None returend: ");
                    break;
            }


            switch(global_no_collective_cause){
                case H5D_MPIO_COLLECTIVE:
                    io_log("global_no_collective_cause: H5D_MPIO_COLLECTIVE: ");
                    break;
                case H5D_MPIO_SET_INDEPENDENT:
                    io_log("global_no_collective_cause: H5D_MPIO_SET_INDEPENDENT: ");
                    break;
                case H5D_MPIO_DATA_TRANSFORMS:
                    io_log("global_no_collective_cause: H5D_MPIO_DATA_TRANSFORMS: ");
                    break;
                //case H5D_MPIO_SET_MPIPOSIX:
                //    io_log("global_no_collective_cause: H5D_MPIO_SET_MPIPOSIX: ");
                //    break;
                case H5D_MPIO_NOT_SIMPLE_OR_SCALAR_DATASPACES:
                    io_log("global_no_collective_cause: H5D_MPIO_NOT_SIMPLE_OR_SCALAR_DATASPACES: ");
                    break;
                //case H5D_MPIO_POINT_SELECTIONS:
                //    io_log("global_no_collective_cause: H5D_MPIO_POINT_SELECTIONS: ");
                //    break;
                // case H5D_MPIO_FILTERS:
                 //   io_log("global_no_collective_cause: H5D_MPIO_FILTERS: ");
                 //   break;
                default :
                    io_log("global_no_collective_cause: None returend: ");
                    break;
            }

            el2 = uptime() - el2;
            io_log("TimeHDF5Write: " << el2 << " s");

            double el3 = uptime();

            /*
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
            if(!file_exist_flag){
            hid_t va_geo_dataspace_id = H5Screate_simple(2, dims, NULL);
            hid_t va_geo_attribute_id = H5Acreate2(file_id, "VPIC-ArrayUDF-GEO", H5T_IEEE_F32BE, va_geo_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
            H5Awrite(va_geo_attribute_id, H5T_NATIVE_FLOAT, attr_data);
            H5Sclose(va_geo_dataspace_id);
            H5Aclose(va_geo_attribute_id);
            }
            */
            if(!file_exist_flag)
                H5Sclose(filespace);
            H5Sclose(memspace);
            H5Pclose(plist_id);
            H5Gclose(group_id);
            H5Fclose(file_id);

            el3 = uptime() - el3;
            io_log("TimeHDF5Close: " << el3 << " s");

            if (mpi_rank == 0)
            {
                char const *output_xml_file = "./field_hdf5/hdf5_field.xdmf";
                char dimensions_3d[128];
                sprintf(dimensions_3d, "%llu %llu %llu", field_global_size[0], field_global_size[1], field_global_size[2]);
                char dimensions_4d[128];
                sprintf(dimensions_4d, "%llu %llu %llu %d", field_global_size[0], field_global_size[1], field_global_size[2], 3);
                char orignal[128];
                sprintf(orignal, "%f %f %f", grid->x0, grid->y0, grid->z0);
                char dxdydz[128];
                sprintf(dxdydz, "%f %f %f", grid->dx, grid->dy, grid->dz);


                // TODO: remove or let the user set
                int field_interval = 1;

                // TODO: remove this dependence on number of steps
                std::cout << "num_step " << num_step << std::endl;

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
                printf("       current step: %zd \n", step_for_viou);
                printf("       current step: %zd \n", step_for_viou);

                //printf("    Simulation time: %f \n", grid->t0);
                printf("             tframe: %d \n", field_tframe);
#endif

                // TODO: this footer dumping is more likely better done in a
                // destructor, rather than hoping a multiple division works out
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
        void dump_particles(
            const char *fbase,
            species_t* sp,
            grid_t* grid,
            int step,
            interpolator_array_t* interpolator_array,
            int ftag
        )
        {
            static int file_index = 0;
            file_index ++;
            int mpi_rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
            double dump_particles_uptime = uptime();
            time_t seconds = time(NULL);
           // printf("Atrank = %d, file_index = %d, dump_particles_uptime = %f, epoch_seconds = %ld  \n ", mpi_rank, file_index, dump_particles_uptime, seconds);


            size_t step_for_viou = step;
            char fname[256];
            char group_name[256];
            char particle_scratch[128];
            char subparticle_scratch[128];

            int np_local;

            float *Pf;
            int *Pi;

            // get the total number of particles. in this example, output only electrons
            //sp = species_list;
            sprintf(particle_scratch, DUMP_DIR_FORMAT, "particle_hdf5");
            FileUtils::makeDirectory(particle_scratch);
            sprintf(subparticle_scratch, "%s/T.%ld/", particle_scratch, step_for_viou);
            FileUtils::makeDirectory(subparticle_scratch);

            // TODO: Allow the user to set this
            int stride_particle_dump = 1;

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

            center_p(sp, interpolator_array);

            ec1 = uptime() - ec1;


            //if(!mpi_rank || mpi_rank == 2047 )
            //    std::cout << "on mpi_rank: " << mpi_rank << ", time in copying particle data: " << ec1 << " s" << ", np_local = " << np_local << std::endl;
	
            #ifndef N_FILE_N_PROCESS
	        int np_local_max, np_local_min;
	        MPI_Reduce(&np_local, &np_local_max, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
	        MPI_Reduce(&np_local, &np_local_min, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
            //io_log("on mpi_rank: " << mpi_rank << ", time in copying particle data: " << ec1 << " s" << ", np_local = " << np_local << ",np_local_max = " << np_local_max << ", local_min = "<< np_local_min);
            #endif

            Pf = (float *)sp->p;
            Pi = (int *)sp->p;

            // open HDF5 file in "particle/T.<step>/" subdirectory
            // filename: eparticle.h5p
            #ifndef N_FILE_N_PROCESS
            sprintf(fname, "%s%s/%s_%ld.h5", daos, subparticle_scratch, sp->name, step_for_viou);
            #else
            sprintf(fname, "%s%s/%s_%ld_p%d.h5", daos, subparticle_scratch, sp->name, step_for_viou, mpi_rank);
            #endif

            sprintf(group_name, "/Timestep_%ld", step_for_viou);
            double el1 = uptime();

            long long total_particles, offset;
            long long numparticles = np_local;

            #ifndef N_FILE_N_PROCESS
                MPI_Allreduce(&numparticles, &total_particles, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
                MPI_Scan(&numparticles, &offset, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
                offset -= numparticles;
            #else
                total_particles = np_local;
                offset = 0;
            #endif


            hid_t file_plist_id = H5Pcreate(H5P_FILE_ACCESS);

            #ifndef N_FILE_N_PROCESS
            H5Pset_fapl_mpio(file_plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
            #endif

            #ifdef H5_ASYNC
            if(!mpi_rank) printf("Enable async on particle data");

            assert(H5Pset_vol_async(file_plist_id));
            #endif

            #ifdef METADATA_COLL_WRITE
                if(!mpi_rank) printf("Enable collective metadata write !\n");
                H5Pset_coll_metadata_write(file_plist_id, TRUE);
            #endif

            hid_t file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, file_plist_id);
            //if(!mpi_rank )
            //    io_log("++Particle H5Fcreate) ");


            hid_t group_id = H5Gcreate(file_id, group_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            //if(!mpi_rank )
            //    io_log("++Particle H5Gcreate) ");

#ifdef  HAS_PARTICLE_COMP
        if(!mpi_rank)
            printf("Using Particle Compound type !\n");
        hid_t  particle_comp_type_it = H5Tcreate (H5T_COMPOUND, sizeof(particle_t));
        H5Tinsert(particle_comp_type_it, "dx", HOFFSET(particle_t, dx), H5T_NATIVE_FLOAT);
        H5Tinsert(particle_comp_type_it, "dy", HOFFSET(particle_t, dy), H5T_NATIVE_FLOAT);
        H5Tinsert(particle_comp_type_it, "dz", HOFFSET(particle_t, dz), H5T_NATIVE_FLOAT);

        H5Tinsert(particle_comp_type_it, "i", HOFFSET(particle_t, i), H5T_NATIVE_INT);

        H5Tinsert(particle_comp_type_it, "ux", HOFFSET(particle_t, ux), H5T_NATIVE_FLOAT);
        H5Tinsert(particle_comp_type_it, "uy", HOFFSET(particle_t, uy), H5T_NATIVE_FLOAT);
        H5Tinsert(particle_comp_type_it, "uz", HOFFSET(particle_t, uz), H5T_NATIVE_FLOAT);
        H5Tinsert(particle_comp_type_it, "w", HOFFSET(particle_t, w), H5T_NATIVE_FLOAT);
#endif

        //MSBprintf("%d: offset %ld %ld %ld \n",mpi_rank,offset,total_particles, numparticles);
            hid_t filespace = H5Screate_simple(1, (hsize_t *)&total_particles, NULL);
            H5Sselect_hyperslab(filespace, H5S_SELECT_SET, (hsize_t *)&offset, NULL, (hsize_t *)&numparticles, NULL);

            //   if(!mpi_rank )
            //    io_log("++Particle H5Sselect_hyperslab) ");

            //plist_id = H5P_DEFAULT;
	        hid_t io_plist_id = H5Pcreate(H5P_DATASET_XFER);

#ifndef N_FILE_N_PROCESS
            #ifdef HAS_INDEPENDENT_IO
                if(!mpi_rank) printf("\n ###\n VPIC Independent I/O! \n ###\n");
	    	    H5Pset_dxpl_mpio(io_plist_id, H5FD_MPIO_INDEPENDENT);
            #else
                H5Pset_dxpl_mpio(io_plist_id, H5FD_MPIO_COLLECTIVE);
            #endif
#endif

            #ifdef H5_ASYNC
            H5Pset_dxpl_async(io_plist_id, true);
            #endif
            //        hsize_t linearspace_count_temp = numparticles;
            //hid_t linearspace = H5Screate_simple(1, &linearspace_count_temp, NULL);

            hsize_t memspace_count_temp;
            hid_t memspace;
#ifdef HAS_PARTICLE_COMP
            memspace_count_temp = numparticles ;
            memspace = H5Screate_simple(1, &memspace_count_temp, NULL);
#else
            memspace_count_temp = numparticles * 8;
            memspace = H5Screate_simple(1, &memspace_count_temp, NULL);
            //MSB   printf("%d: hyperslab %ld\n",mpi_rank, np_local);
            hsize_t memspace_start = 0, memspace_stride = 8, memspace_count = np_local;
            H5Sselect_hyperslab(memspace, H5S_SELECT_SET, &memspace_start, &memspace_stride, &memspace_count, NULL);
#endif
            el1 = uptime() - el1;
            //   if(!mpi_rank || mpi_rank == 2047 )
            //    io_log("Particle TimeHDF5Open): " << el1 << " s"); //Easy to handle results for scripts
            double el2 = uptime();
 int ierr;

#define WRITE_H5_FILE(group_id_p, data_buf_p, type_p, dname_p){\
   hid_t dset_id = H5Dcreate(group_id_p, dname_p, type_p, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); \
   H5Dwrite(dset_id, type_p, memspace, filespace, io_plist_id, data_buf_p);\
   H5Dclose(dset_id);\
}


//MPI_Info_set(info, "romio_cb_write", "disable");
#define WRITE_MPI_FILE(dname_p, offset_p, data_buf_p, count_p, type_p){\
   MPI_File fh;\
   MPI_Status status;\
   sprintf(fname, "%s%s/%s_%ld_%s.h5", daos, subparticle_scratch, sp->name, step_for_viou, dname_p); \
   if(mpi_rank == 0) printf("fname= %s \n", fname);\
   MPI_Info info;\
   MPI_Info_create(&info);\
   MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_WRONLY | MPI_MODE_CREATE, info, &fh);\
   MPI_File_write_at(fh, offset_p, data_buf_p, count_p,type_p, &status);\
   MPI_Info_free(&info);\
   MPI_File_close(&fh);\
}

#ifdef HAS_PARTICLE_COMP
    hid_t dset_id = H5Dcreate(group_id, "ParticleComp", particle_comp_type_it, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dset_id, particle_comp_type_it, memspace, filespace, io_plist_id, sp->p);
    H5Dclose(dset_id);
#else
#ifdef TEST_MPIIO
            //Here we don't use the stripe but just for performance test
         if(!mpi_rank) printf("Test MPI-IO\n");
	        WRITE_MPI_FILE("dX", offset * sizeof(float), Pf, numparticles,  MPI_FLOAT);
            WRITE_MPI_FILE("dY", offset * sizeof(float), Pf, numparticles,  MPI_FLOAT);
            WRITE_MPI_FILE("dZ", offset * sizeof(float), Pf, numparticles,  MPI_FLOAT);
            WRITE_MPI_FILE("i",  offset * sizeof(int),   Pf, numparticles,  MPI_INT);
            WRITE_MPI_FILE("ux", offset * sizeof(float), Pf, numparticles,  MPI_FLOAT);
            WRITE_MPI_FILE("uy", offset * sizeof(float), Pf, numparticles,  MPI_FLOAT);
            WRITE_MPI_FILE("uz", offset * sizeof(float), Pf, numparticles,  MPI_FLOAT);
            WRITE_MPI_FILE("q",  offset * sizeof(float), Pf, numparticles,  MPI_FLOAT);
#else
            #ifndef N_FILE_N_PROCESS
               if(!mpi_rank) printf("Test HDF5-IO Single \n");
            #else
               if(!mpi_rank) printf("Test HDF5-IO N Files N Process\n");
            #endif
               //    if(!mpi_rank )
               //     io_log("++Particle Starting to write ) ");
            WRITE_H5_FILE(group_id,  Pf, H5T_NATIVE_FLOAT, "dX")
            WRITE_H5_FILE(group_id,  Pf+1, H5T_NATIVE_FLOAT, "dY")
            WRITE_H5_FILE(group_id,  Pf+2, H5T_NATIVE_FLOAT, "dZ")
            WRITE_H5_FILE(group_id,  Pi+3, H5T_NATIVE_INT, "i")
            WRITE_H5_FILE(group_id,  Pf+4, H5T_NATIVE_FLOAT, "ux")
            WRITE_H5_FILE(group_id,  Pf+5, H5T_NATIVE_FLOAT, "uy")
            WRITE_H5_FILE(group_id,  Pf+6, H5T_NATIVE_FLOAT, "uz")
            WRITE_H5_FILE(group_id,  Pf+7, H5T_NATIVE_FLOAT, "q")
 #endif
#endif
            el2 = uptime() - el2;
            io_log("Particle TimeHDF5Write: " << el2 << " s");

            double el3 = uptime();
            H5Sclose(memspace);
            H5Sclose(filespace);
            H5Pclose(file_plist_id);
	        H5Pclose(io_plist_id);
            H5Gclose(group_id);



            H5Fclose(file_id);

            #ifdef H5_ASYNC
            H5VLasync_finalize();
            #endif
            el3 = uptime() - el3;
            io_log("Particle TimeHDF5Close: " << el3 << " s");

        }

        void dump_hydro(
            const char *fbase,
            int step,
            hydro_array_t* hydro_array,
            species_t* sp,
            interpolator_array_t* interpolator_array,
            grid_t* grid,
            int ftag
        )
        {
            size_t step_for_viou = step;

            //#define DUMP_INFO_DEBUG 1
            int mpi_size, mpi_rank;
            MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
            MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

            if (!sp)
            {
                ERROR(("Invalid species"));
            }

            clear_hydro_array(hydro_array);
            accumulate_hydro_p(hydro_array, sp, interpolator_array);
            synchronize_hydro_array(hydro_array);

            char hname[256];
            char hydro_scratch[128];
            char subhydro_scratch[128];

            sprintf(hydro_scratch, "./%s", "hydro_hdf5");
            FileUtils::makeDirectory(hydro_scratch);
            sprintf(subhydro_scratch, "%s/T.%zu/", hydro_scratch, step_for_viou);
            FileUtils::makeDirectory(subhydro_scratch);

            sprintf(hname, "%s%s/hydro_%s_%zu.h5", daos, subhydro_scratch, sp->name, step_for_viou);
            double el1 = uptime();
            hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);

            /*
            if(H5Pset_libver_bounds(plist_id, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST) < 0){
                    exit(-1);
            }*/
            //if((fid = H5Fcreate(FILENAME, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id)) < 0)
            //    ERROR_RETURN;

            H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
            #ifdef METADATA_COLL_WRITE
                if(!mpi_rank) printf("Enable collective metadata write !\n");
                H5Pset_coll_metadata_write(plist_id, TRUE);
            #endif
            hid_t file_id = H5Fcreate(hname, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
            H5Pclose(plist_id);

            sprintf(hname, "Timestep_%zu", step_for_viou);
            hid_t group_id = H5Gcreate(file_id, hname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

            el1 = uptime() - el1;
            io_log("TimeHDF5Open: " << el1 << " s"); //Easy to handle results for scripts
            double el2 = uptime();

            // Create a variable list of field values to output.
            //size_t numvars = std::min(global->fdParams.output_vars.bitsum(), total_field_variables);
            //size_t *varlist = new size_t[numvars];

            //for (size_t i(0), c(0); i < total_field_variables; i++)
            //    if (global->fdParams.output_vars.bitset(i))
            //        varlist[c++] = i;

            //printf("\nBEGIN_OUTPUT: numvars = %zu \n", numvars);


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
            UNVOXEL(mpi_rank, mpi_rank_x, mpi_rank_y, mpi_rank_z, grid->gpx, grid->gpy, grid->gpz);

            global_offset[0] = (grid->nx) * mpi_rank_x;
            global_offset[1] = (grid->ny) * mpi_rank_y;
            global_offset[2] = (grid->nz) * mpi_rank_z;

            global_count[0] = (grid->nx);
            global_count[1] = (grid->ny);
            global_count[2] = (grid->nz);

#ifdef DUMP_INFO_DEBUG
            printf("global size   = %llu %llu %llu \n", hydro_global_size[0], hydro_global_size[1], hydro_global_size[2]);
            printf("global_offset = %llu %llu %llu \n", global_offset[0], global_offset[1], global_offset[2]);
            printf("global_count  = %llu %llu %llu \n", global_count[0], global_count[1], global_count[2]);
            printf("mpi-rank = %d, rank index = (%d, %d, %d) \n", mpi_rank, mpi_rank_x, mpi_rank_y, mpi_rank_z);
            fflush(stdout);
#endif

            hid_t filespace = H5Screate_simple(3, hydro_global_size, NULL);
            hid_t memspace = H5Screate_simple(3, hydro_local_size, NULL);
            hid_t dataspace_id;
            /*
            typedef struct hydro {
              float jx, jy, jz, rho; // Current and charge density => <q v_i f>, <q f>
              float px, py, pz, ke;  // Momentum and K.E. density  => <p_i f>, <m c^2 (gamma-1) f>
              float txx, tyy, tzz;   // Stress diagonal            => <p_i v_j f>, i==j
              float tyz, tzx, txy;   // Stress off-diagonal        => <p_i v_j f>, i!=j
              float _pad[2];         // 16-byte align
            } hydro_t;
            */

#if HAS_HYDRO_COMP

            if(!mpi_rank)
              printf("Using Hydro Compund type !\n");

            hid_t  hydro_comp_type_it = H5Tcreate (H5T_COMPOUND, sizeof(hydro_t));

            H5Tinsert(hydro_comp_type_it,"jx", HOFFSET(hydro_t, jx), H5T_NATIVE_FLOAT);
            H5Tinsert(hydro_comp_type_it,"jy", HOFFSET(hydro_t, jy), H5T_NATIVE_FLOAT);
            H5Tinsert(hydro_comp_type_it,"jz", HOFFSET(hydro_t, jz), H5T_NATIVE_FLOAT);
            H5Tinsert(hydro_comp_type_it,"rho", HOFFSET(hydro_t, rho), H5T_NATIVE_FLOAT);

            H5Tinsert(hydro_comp_type_it,"px", HOFFSET(hydro_t, px), H5T_NATIVE_FLOAT);
            H5Tinsert(hydro_comp_type_it,"py", HOFFSET(hydro_t, py), H5T_NATIVE_FLOAT);
            H5Tinsert(hydro_comp_type_it,"pz", HOFFSET(hydro_t, pz), H5T_NATIVE_FLOAT);
            H5Tinsert(hydro_comp_type_it,"ke", HOFFSET(hydro_t, ke), H5T_NATIVE_FLOAT);

            H5Tinsert(hydro_comp_type_it,"txx", HOFFSET(hydro_t, txx), H5T_NATIVE_FLOAT);
            H5Tinsert(hydro_comp_type_it,"tyy", HOFFSET(hydro_t, tyy), H5T_NATIVE_FLOAT);
            H5Tinsert(hydro_comp_type_it,"tzz", HOFFSET(hydro_t, tzz), H5T_NATIVE_FLOAT);

            H5Tinsert(hydro_comp_type_it,"tyz", HOFFSET(hydro_t, tyz), H5T_NATIVE_FLOAT);
            H5Tinsert(hydro_comp_type_it,"tzx", HOFFSET(hydro_t, tzx), H5T_NATIVE_FLOAT);
            H5Tinsert(hydro_comp_type_it,"txy", HOFFSET(hydro_t, txy), H5T_NATIVE_FLOAT);

            //   if(!file_exist_flag){
            dset_id = H5Dcreate(group_id, "HYDRO_DATA_COMP_T", hydro_comp_type_it, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);  // dcpd_id MSB
              //    }else{
              //   dset_id = H5Dopen(group_id, "HYDRO_DATA_COMP_T", H5P_DEFAULT);
              //  }
            dataspace_id = H5Dget_space(dset_id);

            hydro_t *temp_buf = (hydro_t *)malloc(sizeof(hydro_t) * (grid->nx) * (grid->ny) * (grid->nz));
            temp_buf_index = 0;
            for (size_t i(1); i < grid->nx + 1; i++)
                {
                    for (size_t j(1); j < grid->ny + 1; j++)
                    {
                        for (size_t k(1); k < grid->nz + 1; k++)
                        {
                            temp_buf[temp_buf_index].jx = _hydro(i, j, k).jx;
                            temp_buf[temp_buf_index].jy = _hydro(i, j, k).jy;
                            temp_buf[temp_buf_index].jz = _hydro(i, j, k).jz;
                            temp_buf[temp_buf_index].rho = _hydro(i, j, k).rho;
                            temp_buf[temp_buf_index].px = _hydro(i, j, k).px;
                            temp_buf[temp_buf_index].py = _hydro(i, j, k).py;
                            temp_buf[temp_buf_index].pz = _hydro(i, j, k).pz;
                            temp_buf[temp_buf_index].ke = _hydro(i, j, k).ke;
                            temp_buf[temp_buf_index].txx = _hydro(i, j, k).txx;
                            temp_buf[temp_buf_index].tyy = _hydro(i, j, k).tyy;
                            temp_buf[temp_buf_index].tzz = _hydro(i, j, k).tzz;
                            temp_buf[temp_buf_index].tyz = _hydro(i, j, k).tyz;
                            temp_buf[temp_buf_index].tzx = _hydro(i, j, k).tzx;
                            temp_buf[temp_buf_index].txy = _hydro(i, j, k).txy;
                            temp_buf_index = temp_buf_index + 1;
                        }
                    }
                }
            H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, global_offset, NULL, global_count, NULL);
            H5Dwrite(dset_id, hydro_comp_type_it, memspace, dataspace_id, plist_id,  temp_buf);
            H5Sclose(dataspace_id);
            H5Dclose(dset_id);
            H5Tclose(hydro_comp_type_it);
            free(temp_buf);
#else
            float *temp_buf = (float *)malloc(sizeof(float) * (grid->nx) * (grid->ny) * (grid->nz));

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
                            temp_buf[temp_buf_index] = _hydro(i, j, k).ATTRIBUTE_NAME;                                     \
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
#endif

            el2 = uptime() - el2;
            io_log("TimeHDF5Write: " << el2 << " s");

            double el3 = uptime();

            //Write metadata (geo original and geo dx/dy/dz) for ArrayUDF
            /*
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
            H5Aclose(va_geo_attribute_id);*/

            free(temp_buf);
            H5Sclose(filespace);
            H5Sclose(memspace);
            H5Pclose(plist_id);
            H5Gclose(group_id);
            H5Fclose(file_id);

            el3 = uptime() - el3;
            io_log("TimeHDF5Close: " << el3 << " s");

            if (mpi_rank == 0)
            {
                char output_xml_file[128];
                sprintf(output_xml_file, "./%s/%s%s%s", "hydro_hdf5", "hydro-", sp->name, ".xdmf");
                char dimensions_3d[128];
                sprintf(dimensions_3d, "%llu %llu %llu", hydro_global_size[0], hydro_global_size[1], hydro_global_size[2]);
                char dimensions_4d[128];
                sprintf(dimensions_4d, "%llu %llu %llu %d", hydro_global_size[0], hydro_global_size[1], hydro_global_size[2], 3);
                char orignal[128];
                sprintf(orignal, "%f %f %f", grid->x0, grid->y0, grid->z0);
                char dxdydz[128];
                sprintf(dxdydz, "%f %f %f", grid->dx, grid->dy, grid->dz);

                // TODO: remove or let user set
                int hydro_interval = 1;

                // TODO: remove this dependence on number of steps
                int nframes = num_step / hydro_interval + 1;

                const int tframe = tframe_map[sp->id];

#ifdef DUMP_INFO_DEBUG
                printf("         meta file : %s \n", output_xml_file);
                printf(" array dims per var: %s \n", dimensions_3d);
                printf("array dims all vars: %s \n", dimensions_4d);
                printf("            orignal: %s \n", orignal);
                printf("             dxdydz: %s \n", dxdydz);
                printf("            nframes: %d \n", nframes);
                printf("    hydro_fields_interval: %d \n", hydro_interval);
                printf("       current step: %zu \n", step_for_viou);
                printf("    Simulation time: %f \n", grid->t0);
                printf("             tframe: %d \n", tframe);
#endif

                // TODO: why doesnt this just use the cstr?
                char speciesname_new[128];
                sprintf(speciesname_new, "hydro_%s", sp->name);
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
                    create_file_with_header(output_xml_file, dimensions_3d, orignal, dxdydz, nframes, hydro_interval);
                    if (tframe == (nframes - 1))
                    {
                        invert_hydro_xml_item(output_xml_file, speciesname_new, step_for_viou, dimensions_4d, dimensions_3d, 1);
                    }
                    else
                    {
                        invert_hydro_xml_item(output_xml_file, speciesname_new, step_for_viou, dimensions_4d, dimensions_3d, 0);
                    }
                }
                tframe_map[sp->id]++;
            }
        }
};
#endif

#ifdef VPIC_ENABLE_OPENPMD
class OpenPMDDump : public Dump_Strategy {
    public:
        //openPMD::Series* series;
        using Dump_Strategy::Dump_Strategy; // inherit constructor

        //std::string file_type = ".h5";
        std::string file_type = ".bp";

        void dump_fields(
            const char *fbase,
            int step,
            grid_t* grid,
            field_array_t* field_array,
            int ftag
        )
        {
            std::cout << "Writing openPMD data" << std::endl;

            //if (series == nullptr) {
                std::cout << "init series" << std::endl;
                openPMD::Series series = openPMD::Series(
                        fbase,
                        openPMD::AccessType::CREATE,
                        MPI_COMM_WORLD
                );
            //}

            std::cout << "Writing iteration " << step << std::endl;
            auto i = series.iterations[ step ];
            // TODO: it would be nice to set these...
            //series.setAuthor( "Axel Huebl <a.huebl@hzdr.de>");
            //series.setMachine( "Hall Probe 5000, Model 3");
            i.setAttribute( "vacuum", true);

            auto cB = i.meshes["B"];
            auto E = i.meshes["E"];
            auto J = i.meshes["J"];
            auto Tca = i.meshes["Tca"];
            auto Emat = i.meshes["Emat"];
            auto Fmat = i.meshes["Fmat"];
            auto Rho = i.meshes["Rho"];
            auto DivErr = i.meshes["DivErr"];

            // record components
            auto Cbx = cB["x"];
            auto Cby = cB["y"];
            auto Cbz = cB["z"];

            auto Ex = E["x"];
            auto Ey = E["y"];
            auto Ez = E["z"];

            auto Jx = J["x"];
            auto Jy = J["y"];
            auto Jz = J["z"];

            auto Tcax = Tca["x"];
            auto Tcay = Tca["y"];
            auto Tcaz = Tca["z"];

            auto Ematx = Emat["x"];
            auto Ematy = Emat["y"];
            auto Ematz = Emat["z"];

            auto Fmatx = Fmat["x"];
            auto Fmaty = Fmat["y"];
            auto Fmatz = Fmat["z"];

            auto RhoB = Rho["B"];
            auto RhoF = Rho["F"];

            auto DivEErr = DivErr["E"];
            auto DivBErr = DivErr["B"];

            // TODO: set unitDimension so the anaylsis software knows what fields
            // things are
            //
            // // TODO: add timers for the convert and for the write

            size_t gnx = (grid->nx * grid->gpx);
            size_t gny = (grid->ny * grid->gpy);
            size_t gnz = (grid->nz * grid->gpz);
            openPMD::Extent global_extent = {gnx, gny, gnz};

            openPMD::Datatype datatype = openPMD::determineDatatype<float>();
            openPMD::Dataset dataset = openPMD::Dataset(datatype, global_extent);

            Cbx.resetDataset(dataset);
            Cby.resetDataset(dataset);
            Cbz.resetDataset(dataset);

            Ex.resetDataset(dataset);
            Ey.resetDataset(dataset);
            Ez.resetDataset(dataset);

            Jx.resetDataset(dataset);
            Jy.resetDataset(dataset);
            Jz.resetDataset(dataset);

            Tcax.resetDataset(dataset);
            Tcay.resetDataset(dataset);
            Tcaz.resetDataset(dataset);

            Ematx.resetDataset(dataset);
            Ematy.resetDataset(dataset);
            Ematz.resetDataset(dataset);

            Fmatx.resetDataset(dataset);
            Fmaty.resetDataset(dataset);
            Fmatz.resetDataset(dataset);

            RhoB.resetDataset(dataset);
            RhoF.resetDataset(dataset);

            DivEErr.resetDataset(dataset);
            DivBErr.resetDataset(dataset);

            // TODO: hoist this conversion code, as is it used elsewhere
            // Convert rank to local x/y/z
            int rx, ry, rz;
            UNVOXEL(rank, rx, ry, rz, grid->gpx, grid->gpy, grid->gpz);

            size_t nx = grid->nx;
            size_t ny = grid->ny;
            size_t nz = grid->nz;

            // NOTE: this assumes a static mesh decomposition in nx/ny/nz
            size_t global_offset_x = (nx) * rx;
            size_t global_offset_y = (ny) * ry;
            size_t global_offset_z = (nz) * rz;

            openPMD::Offset chunk_offset = {global_offset_x, global_offset_y, global_offset_z};
            openPMD::Extent chunk_extent = {nx, ny, nz};

            std::cout << "Local offset " <<
                " x: " << global_offset_x  <<
                " y: " << global_offset_y  <<
                " z: " << global_offset_z  <<
                std::endl;

            // Store a local copy of the data which we pull out of the AoS
            std::vector<float> cbx_data;
            std::vector<float> cby_data;
            std::vector<float> cbz_data;

            std::vector<float> ex_data;
            std::vector<float> ey_data;
            std::vector<float> ez_data;

            std::vector<float> jx_data;
            std::vector<float> jy_data;
            std::vector<float> jz_data;

            std::vector<float> tcax_data;
            std::vector<float> tcay_data;
            std::vector<float> tcaz_data;

            // TODO: these are material_id (ints not floats)
            std::vector<float> ematx_data;
            std::vector<float> ematy_data;
            std::vector<float> ematz_data;

            std::vector<float> fmatx_data;
            std::vector<float> fmaty_data;
            std::vector<float> fmatz_data;
            // end todo

            std::vector<float> rhob_data;
            std::vector<float> rhof_data;

            std::vector<float> divb_data;
            std::vector<float> dive_data;

            size_t nv = nx * ny * nz;

            // TODO: resize here will zero out the data which we don't need, we
            // could change to a different semantic to avoid this
            cbx_data.resize(nv);
            cby_data.resize(nv);
            cbz_data.resize(nv);

            ex_data.resize(nv);
            ey_data.resize(nv);
            ez_data.resize(nv);

            jx_data.resize(nv);
            jy_data.resize(nv);
            jz_data.resize(nv);

            tcax_data.resize(nv);
            tcay_data.resize(nv);
            tcaz_data.resize(nv);

            ematx_data.resize(nv);
            ematy_data.resize(nv);
            ematz_data.resize(nv);

            fmatx_data.resize(nv);
            fmaty_data.resize(nv);
            fmatz_data.resize(nv);

            rhob_data.resize(nv);
            rhof_data.resize(nv);

            divb_data.resize(nv);
            dive_data.resize(nv);

            // TODO: make this AoS to SoA conversion a function

            // We could do 1D here, but we don't really care about the ghosts, and we
            // can thread over nz/ny (collapsed?)
            // Go over non-ghosts and grab just that data into a dense array
            for (size_t k = 1; k < grid->nz + 1; k++)
            {
                for (size_t j = 1; j < grid->ny + 1; j++)
                {
                    for (size_t i = 1; i < grid->nx + 1; i++)
                    {
                        int local_index  = VOXEL(i-1, j-1, k-1, grid->nx-2, grid->ny-2, grid->nz-2);
                        int global_index = VOXEL(i, j, k, grid->nx, grid->ny, grid->nz);

                        cbx_data[local_index] = field_array->f[global_index].cbx;
                        cby_data[local_index] = field_array->f[global_index].cby;
                        cbz_data[local_index] = field_array->f[global_index].cbz;

                        ex_data[local_index] = field_array->f[global_index].ex;
                        ey_data[local_index] = field_array->f[global_index].ey;
                        ez_data[local_index] = field_array->f[global_index].ez;

                        jx_data[local_index] = field_array->f[global_index].jfx;
                        jy_data[local_index] = field_array->f[global_index].jfy;
                        jz_data[local_index] = field_array->f[global_index].jfz;

                        tcax_data[local_index] = field_array->f[global_index].tcax;
                        tcay_data[local_index] = field_array->f[global_index].tcay;
                        tcaz_data[local_index] = field_array->f[global_index].tcaz;

                        ematx_data[local_index] = field_array->f[global_index].ematx;
                        ematy_data[local_index] = field_array->f[global_index].ematy;
                        ematz_data[local_index] = field_array->f[global_index].ematz;

                        fmatx_data[local_index] = field_array->f[global_index].fmatx;
                        fmaty_data[local_index] = field_array->f[global_index].fmaty;
                        fmatz_data[local_index] = field_array->f[global_index].fmatz;

                        rhob_data[local_index] = field_array->f[global_index].rhob;
                        rhof_data[local_index] = field_array->f[global_index].rhof;

                        dive_data[local_index] = field_array->f[global_index].div_e_err;
                        divb_data[local_index] = field_array->f[global_index].div_b_err;
                    }
                }
            }

            Cbx.storeChunk( cbx_data, chunk_offset, chunk_extent);
            Cby.storeChunk( cby_data, chunk_offset, chunk_extent);
            Cbz.storeChunk( cbz_data, chunk_offset, chunk_extent);

            Ex.storeChunk( ex_data, chunk_offset, chunk_extent);
            Ey.storeChunk( ey_data, chunk_offset, chunk_extent);
            Ez.storeChunk( ez_data, chunk_offset, chunk_extent);

            Jx.storeChunk( jx_data, chunk_offset, chunk_extent);
            Jy.storeChunk( jy_data, chunk_offset, chunk_extent);
            Jz.storeChunk( jz_data, chunk_offset, chunk_extent);

            Tcax.storeChunk( tcax_data, chunk_offset, chunk_extent);
            Tcay.storeChunk( tcay_data, chunk_offset, chunk_extent);
            Tcaz.storeChunk( tcaz_data, chunk_offset, chunk_extent);

            Ematx.storeChunk( ematx_data, chunk_offset, chunk_extent);
            Ematy.storeChunk( ematy_data, chunk_offset, chunk_extent);
            Ematz.storeChunk( ematz_data, chunk_offset, chunk_extent);

            Fmatx.storeChunk( fmatx_data, chunk_offset, chunk_extent);
            Fmaty.storeChunk( fmaty_data, chunk_offset, chunk_extent);
            Fmatz.storeChunk( fmatz_data, chunk_offset, chunk_extent);

            RhoB.storeChunk( rhob_data, chunk_offset, chunk_extent);
            RhoF.storeChunk( rhof_data, chunk_offset, chunk_extent);

            DivEErr.storeChunk( dive_data, chunk_offset, chunk_extent);
            DivBErr.storeChunk( divb_data, chunk_offset, chunk_extent);

            series.flush();
        }

        void dump_particles(
            const char *fbase,
            species_t* sp,
            grid_t* grid,
            int step,
            interpolator_array_t* interpolator_array,
            int ftag
        )
        {
            std::string full_file_name = fbase + file_type;

            std::cout << "writing particles to " << full_file_name << std::endl;

            //if (series == nullptr) {
                openPMD::Series series = openPMD::Series(
                        full_file_name,
                        openPMD::AccessType::CREATE,
                        MPI_COMM_WORLD
                );
            //}

            auto i = series.iterations[ step ];

            // TODO: set these
            i.setTime( (float)step );
            i.setDt(1.0);
            i.setTimeUnitSI(1.0);

            auto& p = i.particles[sp->name];

            const int np = sp->np;

            // TODO: this could be a function call as it's used elsewhere (in hdf5)
            unsigned long long total_particles, offset;
            unsigned long long numparticles = np;
            MPI_Allreduce(&numparticles, &total_particles, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
            MPI_Scan(&numparticles, &offset, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
            offset -= numparticles;

            openPMD::Extent global_extent = {total_particles};
            openPMD::Datatype datatype = openPMD::determineDatatype<float>();
            openPMD::Dataset dataset = openPMD::Dataset(datatype, global_extent);

            auto px = p["position"]["x"];
            auto pxo = p["positionOffset"]["x"];

            auto py = p["position"]["y"];
            auto pyo = p["positionOffset"]["y"];

            auto pz = p["position"]["z"];
            auto pzo = p["positionOffset"]["z"];

            auto ux = p["velocity"]["x"];
            auto uy = p["velocity"]["y"];
            auto uz = p["velocity"]["z"];

            px.resetDataset(dataset);
            pxo.resetDataset(dataset);

            py.resetDataset(dataset);
            pyo.resetDataset(dataset);

            pz.resetDataset(dataset);
            pzo.resetDataset(dataset);

            ux.resetDataset(dataset);
            uy.resetDataset(dataset);
            uz.resetDataset(dataset);
            // convert data to SoA, allowing the user to chunk the operation

            // TODO: Add code the convert to global offsets
#ifndef PMD_MAX_IO_CHUNK // in particles
#define PMD_MAX_IO_CHUNK 16777216; // 512MB total write
#endif
            const int max_chunk = PMD_MAX_IO_CHUNK;

            // Loop over all particles in chunks
            for (int i = 0; i < np; i += max_chunk)
            {
                // We have to be careful as the last chunk may not be full
                // Find how many are left and do that many
                size_t to_write = std::min(np-i, max_chunk);

                // Convert the chunk ready to write
                std::vector<float> x_pos;
                std::vector<float> x_off;
                x_pos.resize(to_write);
                x_off.resize(to_write);

                std::vector<float> y_pos;
                std::vector<float> y_off;
                y_pos.resize(to_write);
                y_off.resize(to_write);

                std::vector<float> z_pos;
                std::vector<float> z_off;
                z_pos.resize(to_write);
                z_off.resize(to_write);

                std::vector<float> ux_pos;
                ux_pos.resize(to_write);

                std::vector<float> uy_pos;
                uy_pos.resize(to_write);

                std::vector<float> uz_pos;
                uz_pos.resize(to_write);

                for (int j = 0; j < to_write; j++)
                {
                    // TODO: do I need to center the particles?
                    auto& particle = sp->p[i+j];

                    x_pos[j] = particle.dx;
                    y_pos[j] = particle.dy;
                    z_pos[j] = particle.dz;

                    ux_pos[j] = particle.ux;
                    uy_pos[j] = particle.uy;
                    uz_pos[j] = particle.uz;

                    std::array<int, 4> gi = global_particle_index(particle.i, grid, rank);
                    x_off[j] = (float)gi[1];
                    y_off[j] = (float)gi[2];
                    z_off[j] = (float)gi[3];
                }

                // Base offset plus i to account for chunks
                auto o = openPMD::Offset{offset + i};
                auto e = openPMD::Extent{to_write};
                px.storeChunk(x_pos, o, e);
                pxo.storeChunk(x_off, o, e);

                py.storeChunk(y_pos, o, e);
                pyo.storeChunk(y_off, o, e);

                pz.storeChunk(z_pos, o, e);
                pzo.storeChunk(z_off, o, e);

                ux.storeChunk(ux_pos, o, e);
                uy.storeChunk(uy_pos, o, e);
                uz.storeChunk(uz_pos, o, e);

                series.flush();
            }
        }
        void dump_hydro(
            const char *fbase,
            int step,
            hydro_array_t* hydro_array,
            species_t* sp,
            interpolator_array_t* interpolator_array,
            grid_t* grid,
            int ftag
        )
        {
            std::string full_file_name = fbase + file_type;

            std::cout << "OpenPMD dumping hydro to " << full_file_name << std::endl;

            //if (series == nullptr) {
                openPMD::Series series = openPMD::Series(
                        full_file_name,
                        openPMD::AccessType::CREATE,
                        MPI_COMM_WORLD
                );
            //}

            auto i = series.iterations[ step ];

            // TODO: set these
            i.setTime( (float)step );
            i.setDt(1.0);
            i.setTimeUnitSI(1.0);

            if( !sp ) ERROR(( "Invalid species \"%s\"", sp->name ));

            // TODO: do we want each backend to have to explicitly call these
            // manually? Or, as it is common, should we hoist it to the VPIC
            // call-site
            clear_hydro_array( hydro_array );
            accumulate_hydro_p( hydro_array, sp, interpolator_array );
            synchronize_hydro_array( hydro_array );

            if( !fbase ) ERROR(( "Invalid filename" ));

            if( rank==0 )
                MESSAGE(("Dumping \"%s\" hydro fields to \"%s\"",sp->name,fbase));

            // Write data
            //float jx, jy, jz, rho; // Current and charge density => <q v_i f>, <q f>
            //float px, py, pz, ke;  // Momentum and K.E. density  => <p_i f>, <m c^2 (gamma-1) f>
            //float txx, tyy, tzz;   // Stress diagonal            => <p_i v_j f>, i==j
            //float tyz, tzx, txy;   // Stress off-diagonal        => <p_i v_j f>, i!=j
            auto J = i.meshes["J"];
            auto P = i.meshes["P"];
            auto T = i.meshes["T"];
            auto _Ke = i.meshes["Ke"];
            auto _Rho = i.meshes["Rho"];

            auto Jx = J["x"];
            auto Jy = J["y"];
            auto Jz = J["z"];

            auto Px = P["x"];
            auto Py = P["y"];
            auto Pz = P["z"];

            auto Txx = T["xx"];
            auto Tyy = T["yy"];
            auto Tzz = T["zz"];
            auto Tyz = T["yz"];
            auto Tzx = T["zx"];
            auto Txy = T["xy"];

            auto Rho = _Rho["rho"]; // TODO: bad name..
            auto Ke = _Ke["ke"]; // TODO: bad name..

            size_t gnx = (grid->nx * grid->gpx);
            size_t gny = (grid->ny * grid->gpy);
            size_t gnz = (grid->nz * grid->gpz);
            openPMD::Extent global_extent = {gnx, gny, gnz};

            openPMD::Datatype datatype = openPMD::determineDatatype<float>();
            openPMD::Dataset dataset = openPMD::Dataset(datatype, global_extent);

            Jx.resetDataset(dataset);
            Jy.resetDataset(dataset);
            Jz.resetDataset(dataset);

            Px.resetDataset(dataset);
            Py.resetDataset(dataset);
            Pz.resetDataset(dataset);

            Txx.resetDataset(dataset);
            Tyy.resetDataset(dataset);
            Tzz.resetDataset(dataset);
            Tyz.resetDataset(dataset);
            Tzx.resetDataset(dataset);
            Txy.resetDataset(dataset);

            Rho.resetDataset(dataset);
            Ke.resetDataset(dataset);

            // TODO: hoist this conversion code, as is it used elsewhere
            // Convert rank to local x/y/z
            int rx, ry, rz;
            UNVOXEL(rank, rx, ry, rz, grid->gpx, grid->gpy, grid->gpz);

            size_t nx = grid->nx;
            size_t ny = grid->ny;
            size_t nz = grid->nz;

            // NOTE: this assumes a static mesh decomposition in nx/ny/nz
            size_t global_offset_x = (nx) * rx;
            size_t global_offset_y = (ny) * ry;
            size_t global_offset_z = (nz) * rz;

            openPMD::Offset chunk_offset = {global_offset_x, global_offset_y, global_offset_z};
            openPMD::Extent chunk_extent = {nx, ny, nz};

            std::cout << "Local offset " <<
                " x: " << global_offset_x  <<
                " y: " << global_offset_y  <<
                " z: " << global_offset_z  <<
                std::endl;

            std::vector<float> jx_data;
            std::vector<float> jy_data;
            std::vector<float> jz_data;

            std::vector<float> px_data;
            std::vector<float> py_data;
            std::vector<float> pz_data;

            std::vector<float> txx_data;
            std::vector<float> tyy_data;
            std::vector<float> tzz_data;
            std::vector<float> tyz_data;
            std::vector<float> tzx_data;
            std::vector<float> txy_data;

            std::vector<float> rho_data;
            std::vector<float> ke_data;

            size_t nv = nx * ny * nz;

            jx_data.resize(nv);
            jy_data.resize(nv);
            jz_data.resize(nv);

            px_data.resize(nv);
            py_data.resize(nv);
            pz_data.resize(nv);

            txx_data.resize(nv);
            tyy_data.resize(nv);
            tzz_data.resize(nv);
            tyz_data.resize(nv);
            tzx_data.resize(nv);
            txy_data.resize(nv);

            rho_data.resize(nv);
            ke_data.resize(nv);

            // Transpose AoS to SoAs
            for (size_t k = 1; k < grid->nz + 1; k++)
            {
                for (size_t j = 1; j < grid->ny + 1; j++)
                {
                    for (size_t i = 1; i < grid->nx + 1; i++)
                    {
                        int local_index  = VOXEL(i-1, j-1, k-1, grid->nx-2, grid->ny-2, grid->nz-2);
                        int global_index = VOXEL(i, j, k, grid->nx, grid->ny, grid->nz);

                        jx_data[local_index] = hydro_array->h[global_index].jx;
                        jy_data[local_index] = hydro_array->h[global_index].jy;
                        jz_data[local_index] = hydro_array->h[global_index].jz;

                        px_data[local_index] = hydro_array->h[global_index].px;
                        py_data[local_index] = hydro_array->h[global_index].py;
                        pz_data[local_index] = hydro_array->h[global_index].pz;

                        txx_data[local_index] = hydro_array->h[global_index].txx;
                        tyy_data[local_index] = hydro_array->h[global_index].tyy;
                        tzz_data[local_index] = hydro_array->h[global_index].tzz;
                        tyz_data[local_index] = hydro_array->h[global_index].tyz;
                        tzx_data[local_index] = hydro_array->h[global_index].tzx;
                        txy_data[local_index] = hydro_array->h[global_index].txy;

                        rho_data[local_index] = hydro_array->h[global_index].rho;
                        ke_data[local_index] = hydro_array->h[global_index].ke;
                    }
                }
            }

            Jx.storeChunk( jx_data, chunk_offset, chunk_extent);
            Jy.storeChunk( jy_data, chunk_offset, chunk_extent);
            Jz.storeChunk( jz_data, chunk_offset, chunk_extent);

            Px.storeChunk( px_data, chunk_offset, chunk_extent);
            Py.storeChunk( py_data, chunk_offset, chunk_extent);
            Pz.storeChunk( pz_data, chunk_offset, chunk_extent);

            Txx.storeChunk( txx_data, chunk_offset, chunk_extent);
            Tyy.storeChunk( tyy_data, chunk_offset, chunk_extent);
            Tzz.storeChunk( tzz_data, chunk_offset, chunk_extent);
            Tyz.storeChunk( tyz_data, chunk_offset, chunk_extent);
            Tzx.storeChunk( tzx_data, chunk_offset, chunk_extent);
            Txy.storeChunk( txy_data, chunk_offset, chunk_extent);

            Rho.storeChunk( rho_data, chunk_offset, chunk_extent);
            Ke.storeChunk( ke_data, chunk_offset, chunk_extent);

            series.flush();
        }
};
#endif

/*
   template <typename Policy = BinaryDump>
   struct IODump : private Policy {
   using Policy::dump_particles;
   using Policy::dump_fields;
   using Policy::dump_hydro;
   };
   */

#endif
