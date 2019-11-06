#ifndef dump_h
#define dump_h

#include <array>
#include "../grid/grid.h"

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
std::array<int, 4> global_particle_index(int local_i, grid_t* grid, int rank);
#endif
