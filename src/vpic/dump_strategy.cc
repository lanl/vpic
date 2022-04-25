// BinaryDump::BinaryDump(int _rank, int _nproc) : Dump_Strategy(_rank, _nproc)
//{
//// empty
//}
#include "dump_strategy.h"

void checkpt_dump_strategy(const Dump_Strategy *ds)
{
    CHECKPT(ds, 1);
    // CHECKPT_SYM(ds->dump_fields);
    // CHECKPT_SYM(ds->dump_hydro);
    // CHECKPT_SYM(ds->dump_particles);
}

Dump_Strategy *
restore_dump_strategy(void)
{
    Dump_Strategy *ds;
    RESTORE(ds);
    // RESTORE_SYM(ds->dump_fields);
    // RESTORE_SYM(ds->dump_hydro);
    // RESTORE_SYM(ds->dump_particles);
    return ds;
}

Dump_Strategy *
new_dump_strategy(DumpStrategyID dump_strategy_id)
{
    Dump_Strategy *ds;
    int i;
    MALLOC(ds, 1);
    CLEAR(ds, 1);

    // Do any post init/restore simulation modifications
    switch (dump_strategy_id)
    {
    case DUMP_STRATEGY_BINARY:
        //::cout << "DUMP_STRATEGY_BINARY \n";
        ds = new BinaryDump(rank(), nproc());
        break;
    case DUMP_STRATEGY_HDF5:
#ifdef VPIC_ENABLE_HDF5
        // std::cout << "DUMP_STRATEGY_HDF5 \n";
        ds = new HDF5Dump(rank(), nproc());
#else
        std::cout << "HDF5Dump is not enabled \n";

#endif
        break;
    case DUMP_STRATEGY_OPENPMD:
#ifdef VPIC_ENABLE_OPENPMD
        // std::cout << "DUMP_STRATEGY_OPENPMD \n";
        ds = new OpenPMDDump(rank(), nproc());
#else
        std::cout << "OpenPMDDump is not enabled \n";
#endif
        break;
    default:
        break;
    }

    // REGISTER_OBJECT(ds, checkpt_dump_strategy, restore_dump_strategy, NULL);
    return ds;
}

void delete_dump_strategy(Dump_Strategy *ds)
{
    if (!ds)
        return;
    UNREGISTER_OBJECT(ds);
    FREE(ds);
}

void BinaryDump::dump_fields(
    const char *fbase,
    int step,
    grid_t *grid,
    field_array_t *field_array,
    int ftag)
{
    char fname[256];
    FileIO fileIO;
    int dim[3];

    if (!fbase)
        ERROR(("Invalid filename"));

    if (rank == 0)
        MESSAGE(("Dumping fields to \"%s\"", fbase));

    if (ftag)
        sprintf(fname, "%s.%li.%i", fbase, (long)step, rank);
    else
        sprintf(fname, "%s.%i", fbase, rank);

    FileIOStatus status = fileIO.open(fname, io_write);
    if (status == fail)
        ERROR(("Could not open \"%s\".", fname));

    /* IMPORTANT: these values are written in WRITE_HEADER_V0 */
    size_t nxout = grid->nx;
    size_t nyout = grid->ny;
    size_t nzout = grid->nz;
    float dxout = grid->dx;
    float dyout = grid->dy;
    float dzout = grid->dz;

    WRITE_HEADER_V0(dump_type::field_dump, -1, 0, fileIO, step, rank, nproc);

    dim[0] = grid->nx + 2;
    dim[1] = grid->ny + 2;
    dim[2] = grid->nz + 2;
    WRITE_ARRAY_HEADER(field_array->f, 3, dim, fileIO);
    fileIO.write(field_array->f, dim[0] * dim[1] * dim[2]);
    if (fileIO.close())
        ERROR(("File close failed on dump fields!!!"));
}

void BinaryDump::dump_particles(
    const char *fbase,
    species_t *sp,
    grid_t *grid,
    int step,
    interpolator_array_t *interpolator_array,
    int ftag)
{
    char fname[256];
    FileIO fileIO;
    int dim[1], buf_start;
    static particle_t *ALIGNED(128) p_buf = NULL;

    // TODO: reconcile this with MAX_IO_CHUNK, and update Cmake option
    // description to explain what backends use it
#define PBUF_SIZE 32768 // 1MB of particles

    if (!sp)
        ERROR(("Invalid species name \"%s\".", sp->name));

    if (!fbase)
        ERROR(("Invalid filename"));

    if (!p_buf)
        MALLOC_ALIGNED(p_buf, PBUF_SIZE, 128);

    if (rank == 0)
        MESSAGE(("Dumping \"%s\" particles to \"%s\"", sp->name, fbase));

    if (ftag)
        sprintf(fname, "%s.%li.%i", fbase, (long)step, rank);
    else
        sprintf(fname, "%s.%i", fbase, rank);
    FileIOStatus status = fileIO.open(fname, io_write);
    if (status == fail)
        ERROR(("Could not open \"%s\"", fname));

    /* IMPORTANT: these values are written in WRITE_HEADER_V0 */
    size_t nxout = grid->nx;
    size_t nyout = grid->ny;
    size_t nzout = grid->nz;
    float dxout = grid->dx;
    float dyout = grid->dy;
    float dzout = grid->dz;

    WRITE_HEADER_V0(dump_type::particle_dump, sp->id, sp->q / sp->m, fileIO, step, rank, nproc);

    dim[0] = sp->np;
    WRITE_ARRAY_HEADER(p_buf, 1, dim, fileIO);

    // Copy a PBUF_SIZE hunk of the particle list into the particle
    // buffer, timecenter it and write it out. This is done this way to
    // guarantee the particle list unchanged while not requiring too
    // much memory.

    // FIXME: WITH A PIPELINED CENTER_P, PBUF NOMINALLY SHOULD BE QUITE
    // LARGE.

    particle_t *sp_p = sp->p;
    sp->p = p_buf;
    int sp_np = sp->np;
    sp->np = 0;
    int sp_max_np = sp->max_np;
    sp->max_np = PBUF_SIZE;
    for (buf_start = 0; buf_start < sp_np; buf_start += PBUF_SIZE)
    {
        sp->np = sp_np - buf_start;
        if (sp->np > PBUF_SIZE)
            sp->np = PBUF_SIZE;
        COPY(sp->p, &sp_p[buf_start], sp->np);
        center_p(sp, interpolator_array);
        fileIO.write(sp->p, sp->np);
    }
    sp->p = sp_p;
    sp->np = sp_np;
    sp->max_np = sp_max_np;

    if (fileIO.close())
        ERROR(("File close failed on dump particles!!!"));
}
void BinaryDump::dump_hydro(
    const char *fbase,
    int step,
    hydro_array_t *hydro_array,
    species_t *sp,
    interpolator_array_t *interpolator_array,
    grid_t *grid,
    int ftag)
{
    char fname[256];
    FileIO fileIO;
    int dim[3];

    if (!sp)
        ERROR(("Invalid species \"%s\"", sp->name));

    clear_hydro_array(hydro_array);
    accumulate_hydro_p(hydro_array, sp, interpolator_array);
    synchronize_hydro_array(hydro_array);

    if (!fbase)
        ERROR(("Invalid filename"));

    if (rank == 0)
        MESSAGE(("Dumping \"%s\" hydro fields to \"%s\"", sp->name, fbase));

    if (ftag)
        sprintf(fname, "%s.%li.%i", fbase, (long)step, rank);
    else
        sprintf(fname, "%s.%i", fbase, rank);
    FileIOStatus status = fileIO.open(fname, io_write);
    if (status == fail)
        ERROR(("Could not open \"%s\".", fname));

    /* IMPORTANT: these values are written in WRITE_HEADER_V0 */
    size_t nxout = grid->nx;
    size_t nyout = grid->ny;
    size_t nzout = grid->nz;
    float dxout = grid->dx;
    float dyout = grid->dy;
    float dzout = grid->dz;

    WRITE_HEADER_V0(dump_type::hydro_dump, sp->id, sp->q / sp->m, fileIO, step, rank, nproc);

    dim[0] = grid->nx + 2;
    dim[1] = grid->ny + 2;
    dim[2] = grid->nz + 2;
    WRITE_ARRAY_HEADER(hydro_array->h, 3, dim, fileIO);
    fileIO.write(hydro_array->h, dim[0] * dim[1] * dim[2]);
    if (fileIO.close())
        ERROR(("File close failed on dump hydro!!!"));
}
