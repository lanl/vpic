#ifndef dump_h
#define dump_h

#include <array>

// TODO: should this be an enum?
namespace dump_type {
  const int grid_dump = 0;
  const int field_dump = 1;
  const int hydro_dump = 2;
  const int particle_dump = 3;
  const int restart_dump = 4;
  const int history_dump = 5;
} // namespace

// TODO: namesapce?
std::array<int, 4> global_particle_index(int local_i, grid_t* grid, int rank)
{
    int ix, iy, iz, rx, ry, rz;
    // Convert rank to local x/y/z
    UNVOXEL(rank, rx, ry, rz, grid->gpx, grid->gpy, grid->gpz);

    // Calculate local ix/iy/iz
    UNVOXEL(local_i, ix, iy, iz, grid->nx+2, grid->ny+2, grid->nz+2);

    // Account for the "first" ghost cell
    ix = ix - 1;
    iy = iy - 1;
    iz = iz - 1;

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

    return { global_i, gix, giy, giz };
}
#endif
