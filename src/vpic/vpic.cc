/*
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Heavily revised and extended from earlier V4PIC versions
 *
 */

#include "vpic.h"
#include "dump_strategy.h"

/* Note that, when a vpic_simulation is created (and thus registered
   with the checkpt service), it is created empty; none of the simulation
   objects on which it depends have been created yet. (These get created
   as the simulation configuration is assembled in the input deck.) As
   such, vpic_simulation pointers do not point to older objects.  We
   thus need to use the checkpt_fptr functions and write a reanimator.

   FIXME: We could avoid this by giving each one of the objects pointed
   to proper resize semantics so we could then create the objects during
   vpic_simulation construction (as opposed to after it). */

void checkpt_vpic_simulation(const vpic_simulation *vpic)
{
  CHECKPT(vpic, 1);
  // CHECKPT_PTR(vpic->dump_strategy);
  CHECKPT_PTR(vpic->entropy);
  CHECKPT_PTR(vpic->sync_entropy);
  CHECKPT_PTR(vpic->grid);
  CHECKPT_FPTR(vpic->material_list);
  CHECKPT_FPTR(vpic->field_array);
  CHECKPT_FPTR(vpic->interpolator_array);
  CHECKPT_FPTR(vpic->accumulator_array);
  CHECKPT_FPTR(vpic->hydro_array);
  CHECKPT_FPTR(vpic->species_list);
  CHECKPT_FPTR(vpic->particle_bc_list);
  CHECKPT_FPTR(vpic->emitter_list);
  CHECKPT_FPTR(vpic->collision_op_list);
}

vpic_simulation *
restore_vpic_simulation(void)
{
  vpic_simulation *vpic;
  RESTORE(vpic);
  // RESTORE_PTR(vpic->dump_strategy);
  RESTORE_PTR(vpic->entropy);
  RESTORE_PTR(vpic->sync_entropy);
  RESTORE_PTR(vpic->grid);
  RESTORE_FPTR(vpic->material_list);
  RESTORE_FPTR(vpic->field_array);
  RESTORE_FPTR(vpic->interpolator_array);
  RESTORE_FPTR(vpic->accumulator_array);
  RESTORE_FPTR(vpic->hydro_array);
  RESTORE_FPTR(vpic->species_list);
  RESTORE_FPTR(vpic->particle_bc_list);
  RESTORE_FPTR(vpic->emitter_list);
  RESTORE_FPTR(vpic->collision_op_list);

  // printf("\n\n\n vpic->dump_strategy_id is restored ! \n\n\n");
  vpic->dump_strategy = new_dump_strategy(vpic->dump_strategy_id);
  vpic->dump_strategy->num_step = vpic->num_step;

  //   // Do any post init/restore simulation modifications
  //   switch (vpic->dump_strategy_id)
  //   {
  //   case DUMP_STRATEGY_BINARY:
  //     //::cout << "DUMP_STRATEGY_BINARY \n";
  //     vpic->dump_strategy = new BinaryDump(rank(), nproc());
  //     break;
  //   case DUMP_STRATEGY_HDF5:
  // #ifdef VPIC_ENABLE_HDF5
  //     // std::cout << "DUMP_STRATEGY_HDF5 \n";
  //     vpic->dump_strategy = new HDF5Dump(rank(), nproc());
  // #else
  //     std::cout << "HDF5Dump is not enabled \n";

  // #endif
  //     break;
  //   case DUMP_STRATEGY_OPENPMD:
  // #ifdef VPIC_ENABLE_OPENPMD
  //     // std::cout << "DUMP_STRATEGY_OPENPMD \n";
  //     vpic->dump_strategy = new OpenPMDDump(rank(), nproc());
  // #else
  //     std::cout << "OpenPMDDump is not enabled \n";
  // #endif
  //     break;
  //   default:
  //     break;
  //   }

  return vpic;
}

void reanimate_vpic_simulation(vpic_simulation *vpic)
{
  REANIMATE_FPTR(vpic->material_list);
  REANIMATE_FPTR(vpic->field_array);
  REANIMATE_FPTR(vpic->interpolator_array);
  REANIMATE_FPTR(vpic->accumulator_array);
  REANIMATE_FPTR(vpic->hydro_array);
  REANIMATE_FPTR(vpic->species_list);
  REANIMATE_FPTR(vpic->particle_bc_list);
  REANIMATE_FPTR(vpic->emitter_list);
  REANIMATE_FPTR(vpic->collision_op_list);
}

vpic_simulation::vpic_simulation()
{
  // TODO: why is this a good idea?
  // Is this just trying to 0 initialize everything?
  // CLEAR( this, 1 );

  // Now done in the class def / header
  ///* Set non-zero defaults */
  // verbose = 1;
  // num_comm_round = 3;
  // num_div_e_round = 2;
  // num_div_b_round = 2;

#if defined(VPIC_USE_PTHREADS) // Pthreads case.
  int n_rng = serial.n_pipeline;
  if (n_rng < thread.n_pipeline)
    n_rng = thread.n_pipeline;

#elif defined(VPIC_USE_OPENMP) // OpenMP case.
  int n_rng = omp_helper.n_pipeline;

#else // Error case.
#error "VPIC_USE_OPENMP or VPIC_USE_PTHREADS must be specified"

#endif

  // int                           n_rng = serial.n_pipeline;
  // if( n_rng<thread.n_pipeline ) n_rng = thread.n_pipeline;

  // # if defined(CELL_PPU_BUILD) && defined(USE_CELL_SPUS)
  //   if( n_rng<spu.n_pipeline    ) n_rng = spu.n_pipeline;
  // # endif

  n_rng++;

  entropy = new_rng_pool(n_rng, 0, 0);
  sync_entropy = new_rng_pool(n_rng, 0, 1);
  grid = new_grid();
  dump_strategy = new_dump_strategy(dump_strategy_id);

  REGISTER_OBJECT(this, checkpt_vpic_simulation,
                  restore_vpic_simulation, reanimate_vpic_simulation);

  // Initialize the dump strategy to use the binary dumpin, assuming the user
  // may overwrite this later
  // dump_strategy = std::unique_ptr<Dump_Strategy>(new BinaryDump( rank(), nproc() ));
  // enable_binary_dump();

  // TODO: this this still makes sense now we have a dump strategy
  // #ifdef VPIC_ENABLE_HDF5
  // Default init hdf5 dump flags
  // field_interval = 1;
  // hydro_interval = 1;
  // field_dump_flag = field_dump_flag_t();
  // hydro_dump_flag = hydro_dump_flag_t();
  // #endif

  field_interval = 1;
  hydro_interval = 1;
}

vpic_simulation::~vpic_simulation()
{
  UNREGISTER_OBJECT(this);
  delete_collision_op_list(collision_op_list);
  delete_emitter_list(emitter_list);
  delete_particle_bc_list(particle_bc_list);
  delete_species_list(species_list);
  delete_hydro_array(hydro_array);
  delete_accumulator_array(accumulator_array);
  delete_interpolator_array(interpolator_array);
  delete_field_array(field_array);
  delete_material_list(material_list);
  delete_grid(grid);
  delete_rng_pool(sync_entropy);
  delete_rng_pool(entropy);
}
