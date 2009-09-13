/*
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Heavily revised and extended from earlier V4PIC versions
 *
 */
 
#include "vpic.hxx"
 
/* At some point, entropy_pool should become fully fleshed out and
   publicly exposed. */

void
checkpt_entropy_pool( const mt_rng_t ** pool ) {
  int r;
  for( r=0; pool[r]; r++ ) /**/;
  CHECKPT( pool, r+1 );
  for( r=0; pool[r]; r++ ) CHECKPT_PTR( pool[r] );
}

mt_rng_t **
restore_entropy_pool( void ) {
  mt_rng_t ** pool;
  int r;
  RESTORE( pool );
  for( r=0; pool[r]; r++ ) RESTORE_PTR( pool[r] );
  return pool;
}

mt_rng_t **
new_entropy_pool( int n_rng,
                  int seed,
                  int sync ) {
  mt_rng_t ** pool;
  int r;
  if( n_rng<1 ) ERROR(( "Bad args" ));
  seed = (sync ? world_size : world_rank) + (world_size+1)*n_rng*seed;
  MALLOC( pool, n_rng+1 );
  for( r=0; r<n_rng; r++ ) pool[r] = new_mt_rng( seed + (world_size+1)*r );
  pool[n_rng] = NULL; 
  REGISTER_OBJECT( pool, checkpt_entropy_pool, restore_entropy_pool, NULL );
  return pool;
}

void
delete_entropy_pool( mt_rng_t ** pool ) {
  int r;
  if( !pool ) return;
  UNREGISTER_OBJECT( pool );
  for( r=0; pool[r]; r++ ) delete_mt_rng( pool[r] );
  FREE( pool );
}

/* Note that, when a vpic_simulation is created (and thus registered
   with the checkpt service), it is created empty; none of the simulation
   objects on which it depends have been created yet. (These get created
   as the simulation configuration is assembled in the input deck.) As
   such, vpic_simulation pointers do not point to older objects.  We
   thus need to use the checkpt_fptr functions and write a reanimator.

   FIXME: We could avoid this by giving each one of the objects pointed
   to proper resize semantics so we could then create the objects during
   vpic_simulation construction (as opposed to after it). */

void
checkpt_vpic_simulation( const vpic_simulation * vpic ) {
  CHECKPT( vpic, 1 );
  CHECKPT_PTR( vpic->rng );
  CHECKPT_PTR( vpic->sync_rng );
  CHECKPT_PTR( vpic->grid );
  CHECKPT_FPTR( vpic->material_list );
  CHECKPT_FPTR( vpic->field_array );
  CHECKPT_FPTR( vpic->interpolator_array );
  CHECKPT_FPTR( vpic->accumulator_array );
  CHECKPT_FPTR( vpic->hydro_array );
  CHECKPT_FPTR( vpic->species_list );
  CHECKPT_FPTR( vpic->particle_bc_list );
  CHECKPT_FPTR( vpic->emitter_list );
  CHECKPT_FPTR( vpic->collision_op_list );
}

vpic_simulation *
restore_vpic_simulation( void ) {
  vpic_simulation * vpic;
  RESTORE( vpic );
  RESTORE_PTR( vpic->rng );
  RESTORE_PTR( vpic->sync_rng );
  RESTORE_PTR( vpic->grid );
  RESTORE_FPTR( vpic->material_list );
  RESTORE_FPTR( vpic->field_array );
  RESTORE_FPTR( vpic->interpolator_array );
  RESTORE_FPTR( vpic->accumulator_array );
  RESTORE_FPTR( vpic->hydro_array );
  RESTORE_FPTR( vpic->species_list );
  RESTORE_FPTR( vpic->particle_bc_list );
  RESTORE_FPTR( vpic->emitter_list );
  RESTORE_FPTR( vpic->collision_op_list );
  return vpic;
}

void
reanimate_vpic_simulation( vpic_simulation * vpic ) {
  REANIMATE_FPTR( vpic->material_list );
  REANIMATE_FPTR( vpic->field_array );
  REANIMATE_FPTR( vpic->interpolator_array );
  REANIMATE_FPTR( vpic->accumulator_array );
  REANIMATE_FPTR( vpic->hydro_array );
  REANIMATE_FPTR( vpic->species_list );
  REANIMATE_FPTR( vpic->particle_bc_list );
  REANIMATE_FPTR( vpic->emitter_list );
  REANIMATE_FPTR( vpic->collision_op_list );
}


vpic_simulation::vpic_simulation() {
  CLEAR( this, 1 );

  /* Set non-zero defaults */
  verbose = 1;
  num_comm_round = 3;
  num_div_e_round = 2;
  num_div_b_round = 2;

  int                           n_rng = serial.n_pipeline;
  if( n_rng<thread.n_pipeline ) n_rng = thread.n_pipeline;
# if defined(CELL_PPU_BUILD) && defined(USE_CELL_SPUS)
  if( n_rng<spu.n_pipeline    ) n_rng = spu.n_pipeline;
# endif
  n_rng++; 

  /**/ rng = new_entropy_pool( n_rng, 0, 0 );
  sync_rng = new_entropy_pool( n_rng, 0, 1 );
  grid     = new_grid();

  REGISTER_OBJECT( this, checkpt_vpic_simulation,
                   restore_vpic_simulation, reanimate_vpic_simulation );
}
 
vpic_simulation::~vpic_simulation() {
  UNREGISTER_OBJECT( this );
  delete_emitter_list( emitter_list );
  delete_particle_bc_list( particle_bc_list );
  delete_species_list( species_list );
  delete_hydro_array( hydro_array );
  delete_accumulator_array( accumulator_array );
  delete_interpolator_array( interpolator_array );
  delete_field_array( field_array );
  delete_material_list( material_list );
  delete_grid( grid );
  delete_entropy_pool( sync_rng );
  delete_entropy_pool( rng );
}
 
