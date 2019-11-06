#ifndef Dump_Strategy_h
#define Dump_Strategy_h

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
#endif

#ifdef VPIC_ENABLE_OPENPMD
#include <openPMD/openPMD.hpp>
#endif

class Dump_Strategy {
    public:
    int rank, nproc;

    Dump_Strategy(int _rank, int _nproc) : rank(_rank), nproc(_nproc) { } // empty
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
        BinaryDump(int _rank, int _nproc) : Dump_Strategy(_rank, _nproc){ } // empty

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
class HDF5Dump : public Dump_Strategy {
    public:
        using Dump_Strategy::Dump_Strategy; // inherit constructor
#define DUMP_DIR_FORMAT "./%s"

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
        void dump_fields(
            const char *fbase,
            int step,
            grid_t* grid,
            field_array_t* field_array,
            int ftag
        )
        {
            size_t step_for_viou = step;

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
            FileUtils::makeDirectory(field_scratch);
            sprintf(subfield_scratch, "%s/T.%zu/", field_scratch, step_for_viou);
            FileUtils::makeDirectory(subfield_scratch);

            sprintf(fname, "%s/%s_%zu.h5", subfield_scratch, "fields", step_for_viou);
            double el1 = uptime();
            hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
            H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
            hid_t file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
            H5Pclose(plist_id);

            sprintf(fname, "Timestep_%zu", step_for_viou);
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

            // Convert rank to local decomposition
            int rx, ry, rz;
            UNVOXEL(mpi_rank, rx, ry, rz, grid->gpx, grid->gpy, grid->gpz);

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
            printf("global size   = %d  %d %d \n", field_global_size[0], field_global_size[1], field_global_size[2]);
            printf("global_offset = %d %d %d \n", global_offset[0], global_offset[1], global_offset[2]);
            printf("global_count  = %d  %d %d \n", global_count[0], global_count[1], global_count[2]);
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
            FileUtils::makeDirector(particle_scratch);
            sprintf(subparticle_scratch, "%s/T.%ld/", particle_scratch, step_for_viou);
            FileUtils::makeDirector(subparticle_scratch);

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

            // Don't need, can just use H5S_ALL
            //hsize_t linearspace_count_temp = numparticles;
            //hid_t linearspace = H5Screate_simple(1, &linearspace_count_temp, NULL);

            plist_id = H5Pcreate(H5P_DATASET_XFER);

            //Comment out for test only
            H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
            H5Sselect_hyperslab(filespace, H5S_SELECT_SET, (hsize_t *)&offset, NULL, (hsize_t *)&numparticles, NULL);

            hsize_t memspace_start = 0, memspace_stride = 8, memspace_count = np_local;
            H5Sselect_hyperslab(memspace, H5S_SELECT_SET, &memspace_start, &memspace_stride, &memspace_count, NULL);

            el1 = uptime() - el1;
            //sim_log("Particle TimeHDF5Open): " << el1 << " s"); //Easy to handle results for scripts

            //double el2 = uptime();

            // This point offset is silly, and loses the type safety (pf+1)
            hid_t dset_id = H5Dcreate(group_id, "dX", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            int ierr = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, Pf);
            H5Dclose(dset_id);

            dset_id = H5Dcreate(group_id, "dY", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            ierr = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, Pf + 1);
            H5Dclose(dset_id);

            dset_id = H5Dcreate(group_id, "dZ", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            ierr = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, Pf + 2);
            H5Dclose(dset_id);

#define OUTPUT_CONVERT_GLOBAL_ID 1
#ifdef OUTPUT_CONVERT_GLOBAL_ID
            // TODO: make a function out of this too, its used in openpmd
            std::vector<int> global_pi;
            global_pi.reserve(numparticles);
            // TODO: this could be parallel
            for (int i = 0; i < numparticles; i++)
            {
                int local_i = sp->p[i].i;

                int ix, iy, iz, rx, ry, rz;

                // Convert rank to local x/y/z
                UNVOXEL(rank, rx, ry, rz, grid->gpx, grid->gpy, grid->gpz);

                // Calculate local ix/iy/iz
                UNVOXEL(local_i, ix, iy, iz, grid->nx+2, grid->ny+2, grid->nz+2);

                // Convert ix/iy/iz to global
                int gix = ix + (grid->nx * (rx));
                int giy = iy + (grid->ny * (ry));
                int giz = iz + (grid->nz * (rz));

                // calculate global grid sizes
                int gnx = grid->nx * grid->gpx;
                int gny = grid->ny * grid->gpy;
                int gnz = grid->nz * grid->gpz;

                // TODO: find a better way to account for the hard coded ghosts in VOXEL
                int global_i = VOXEL(gix, giy, giz, gnx-2, gny-2, gnz-2);

                //std::cout << rank << " local i " << local_i << " becomes " << global_i << std::endl;
                global_pi[i] = global_i;
            }

#undef UNVOXEL
            dset_id = H5Dcreate(group_id, "i", H5T_NATIVE_INT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            ierr = H5Dwrite(dset_id, H5T_NATIVE_INT, H5S_ALL, filespace, plist_id, global_pi.data());
            H5Dclose(dset_id);

#else
            dset_id = H5Dcreate(group_id, "i", H5T_NATIVE_INT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            ierr = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, Pi + 3);
            H5Dclose(dset_id);
#endif

            dset_id = H5Dcreate(group_id, "Ux", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            ierr = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, Pf + 4);
            H5Dclose(dset_id);

            dset_id = H5Dcreate(group_id, "Uy", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            ierr = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, Pf + 5);
            H5Dclose(dset_id);

            dset_id = H5Dcreate(group_id, "Uz", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            ierr = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, Pf + 6);
            H5Dclose(dset_id);

            dset_id = H5Dcreate(group_id, "q", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            ierr = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, Pf + 7);
            H5Dclose(dset_id);

            //el2 = uptime() - el2;
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

            clear_hydro_array(hydro_array);
            accumulate_hydro_p(hydro_array, sp, interpolator_array);
            synchronize_hydro_array(hydro_array);

            char hname[256];
            char hydro_scratch[128];
            char subhydro_scratch[128];

            sprintf(hydro_scratch, "./%s", "hydro_hdf5");
            FileUtils::makeDirector(hydro_scratch);
            sprintf(subhydro_scratch, "%s/T.%zu/", hydro_scratch, step_for_viou);
            FileUtils::makeDirector(subhydro_scratch);

            sprintf(hname, "%s/hydro_%s_%zu.h5", subhydro_scratch, speciesname, step_for_viou);
            double el1 = uptime();
            hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
            H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
            hid_t file_id = H5Fcreate(hname, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
            H5Pclose(plist_id);

            sprintf(hname, "Timestep_%zu", step_for_viou);
            hid_t group_id = H5Gcreate(file_id, hname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

            el1 = uptime() - el1;
            //sim_log("TimeHDF5Open: " << el1 << " s"); //Easy to handle results for scripts
            //double el2 = uptime();

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
            printf("global size   = %d %d %d \n", hydro_global_size[0], hydro_global_size[1], hydro_global_size[2]);
            printf("global_offset = %d %d %d \n", global_offset[0], global_offset[1], global_offset[2]);
            printf("global_count  = %d %d %d \n", global_count[0], global_count[1], global_count[2]);
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

            //el2 = uptime() - el2;
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
        static openPMD::Series* series;
        using Dump_Strategy::Dump_Strategy; // inherit constructor
        void dump_fields(
            const char *fbase,
            int step,
            grid_t* grid,
            field_array_t* field_array,
            int ftag
        )
        {
            std::cout << "Writing openPMD data" << std::endl;

            if (series == nullptr) {
                std::cout << "init series" << std::endl;
                series = new openPMD::Series(
                        fbase,
                        openPMD::AccessType::CREATE,
                        MPI_COMM_WORLD
                        );
            }

            std::cout << "Writing itration " << step << std::endl;
            auto i = series->iterations[ step ];
            // TODO: it would be nice to set these...
            //series.setAuthor( "Axel Huebl <a.huebl@hzdr.de>");
            //series.setMachine( "Hall Probe 5000, Model 3");
            i.setAttribute( "vacuum", true);

            auto cB = i.meshes["B"];
            auto E = i.meshes["E"];
            auto J = i.meshes["J"];

            // record components
            auto cbx = cB["x"];
            auto cby = cB["y"];
            auto cbz = cB["z"];

            auto Ex = E["x"];
            auto Ey = E["y"];
            auto Ez = E["z"];

            auto Jx = J["x"];
            auto Jy = J["y"];
            auto Jz = J["z"];

            // TODO: set unitDimension so the anaylsis software knows what fields
            // things are

            size_t gnx = (grid->nx * grid->gpx);
            size_t gny = (grid->ny * grid->gpy);
            size_t gnz = (grid->nz * grid->gpz);
            openPMD::Extent global_extent = {gny, gny, gnz};

            openPMD::Datatype datatype = openPMD::determineDatatype<float>();
            openPMD::Dataset dataset = openPMD::Dataset(datatype, global_extent);

            cbx.resetDataset(dataset);
            cby.resetDataset(dataset);
            cbz.resetDataset(dataset);

            Ex.resetDataset(dataset);
            Ey.resetDataset(dataset);
            Ez.resetDataset(dataset);

            Jx.resetDataset(dataset);
            Jy.resetDataset(dataset);
            Jz.resetDataset(dataset);

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

            size_t nv = nx * ny * nz;

            cbx_data.reserve(nv);
            cby_data.reserve(nv);
            cbz_data.reserve(nv);

            ex_data.reserve(nv);
            ey_data.reserve(nv);
            ez_data.reserve(nv);

            jx_data.reserve(nv);
            jy_data.reserve(nv);
            jz_data.reserve(nv);

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
                    }
                }
            }

            cbx.storeChunk( cbx_data, chunk_offset, chunk_extent);
            cby.storeChunk( cby_data, chunk_offset, chunk_extent);
            cbz.storeChunk( cbz_data, chunk_offset, chunk_extent);

            Ex.storeChunk( ex_data, chunk_offset, chunk_extent);
            Ey.storeChunk( ey_data, chunk_offset, chunk_extent);
            Ez.storeChunk( ez_data, chunk_offset, chunk_extent);

            Jx.storeChunk( jx_data, chunk_offset, chunk_extent);
            Jy.storeChunk( jy_data, chunk_offset, chunk_extent);
            Jz.storeChunk( jz_data, chunk_offset, chunk_extent);

            series->flush();
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
            if (series == nullptr) {
                std::cout << "init series" << std::endl;
                series = new openPMD::Series(
                        fbase,
                        openPMD::AccessType::CREATE,
                        MPI_COMM_WORLD
                        );
            }

            auto i = series->iterations[ step ];

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

            px.resetDataset(dataset);
            pxo.resetDataset(dataset);

            // convert data to SoA, allowing the user to chunk the operation
            const int max_chunk = 32768*8; // 1MB SoA
            // Loop over all particles in chunks
            for (int i = 0; i < np; i += max_chunk)
            {
                // We have to be careful as the last chunk may not be full
                // Find how many are left and do that many
                size_t to_write = std::min(np-i, max_chunk);

                // Convert the chunk ready to write
                std::vector<float> x_pos;
                std::vector<float> x_off;
                x_pos.reserve(to_write);
                x_off.reserve(to_write);

                for (int j = 0; j < to_write; j++)
                {
                    // TODO: do I need to center the particles?
                    auto& particle = sp->p[i+j];
                    x_pos[j] = particle.dx;
                    std::array<int, 4> gi = global_particle_index(particle.i, grid, rank);
                    x_off[j] = (float)gi[1];
                }

                // Base offset plus i to account for chunks
                auto o = openPMD::Offset{offset + i};
                auto e = openPMD::Extent{to_write};
                px.storeChunk(x_pos, o, e);
                pxo.storeChunk(x_off, o, e);
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
