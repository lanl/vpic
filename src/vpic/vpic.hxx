/* 
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Heavily revised and extended from earlier V4PIC versions
 *
 */

#ifndef _vpic_hxx_
#define _vpic_hxx_

#include <emitter.h>
#include <boundary.h>

#ifndef USER_GLOBAL_SIZE
#define USER_GLOBAL_SIZE 16384
#endif

class vpic_simulation {
public:
  vpic_simulation();
  ~vpic_simulation();
  void initialize( int argc, char **argv );
  void restart( const char *filebase );
  void modify_runparams( const char *fname );
  int advance( void );

  // some helpers that need to be exposed to
  // remove direct MPI calls

  inline double rank(void) {
    return (double)mp_rank(grid->mp);
  }

  inline void barrier(void) {
    mp_barrier(grid->mp);
  }

private:

  // Directly initialized by user; saved in a restart dump

  int step;                 // Number of steps taken so far
  int num_step;             // Number of steps to take
  int status_interval;      // How often to print status messages
  int clean_div_e_interval; // How often to clean div e
  int clean_div_b_interval; // How often to clean div b
  int sync_shared_interval; // How often to synchronize shared faces

  double quota; 
  int restart_interval; 
  int hydro_interval; 
  int field_interval; 
  int particle_interval; 

  // Helper initialized by user; restart-saved

  mt_handle rng;               // Random helpers (seed defaults to rank)
  material_t *material_list;   // Material helpers
  grid_t *grid;                // Grid helpers
  field_t * ALIGNED(16) field; // Grid helpers / field helpers
  species_t *species_list;     // Species helpers / particle helpers
  emitter_t *emitter_list;     // Particle emitters

  // Helper initialized by user; not restart-saved (can be derived from above)

  material_coefficient_t * ALIGNED(16) material_coefficient;
  /**/                                        // Material helpers
  interpolator_t * ALIGNED(128) interpolator; // Grid helpers
  accumulator_t * ALIGNED(128) accumulator;   // Grid helpers
  hydro_t * ALIGNED(16) hydro;                // Grid helpers
  species_t **species_lookup;                 // Species / particle helpers

  // Internal use only variables; restart saved

  double p_time; // Time spent pushing particles since last status update
  double g_time; // Time spent processing guard list since last status update
  double f_time; // Time spent processing fields since last status update
  double u_time; // Time spent in user functions since last status update

  // User defined variables; restart saved

  char user_global[USER_GLOBAL_SIZE];
  // Note: user_global is aliased with user_global_t (see deck_wrapper.cxx)

  ///////////////
  // Dump helpers

  // Text dumps
  void dump_energies( const char *fname, int append = 1 );
  void dump_materials( const char *fname );
  void dump_species( const char *fname );

  // Binary dumps
  void dump_grid( const char *fbase );
  void dump_fields( const char *fbase, int fname_tag = 1 );
  void dump_hydro( const char *sp_name, const char *fbase,
                   int fname_tag = 1 );
  void dump_particles( const char *sp_name, const char *fbase,
                       int fname_tag = 1 );
  void dump_restart( const char *fbase, int fname_tag = 1 );

  ///////////////
  // Grid helpers

  void IUO_allocate_fields(void);

  // The below functions automatically create partition simple grids with
  // simple boundary conditions on the edges.

  inline void define_periodic_grid( double xl,  double yl,  double zl,
                                    double xh,  double yh,  double zh,
				    double gnx, double gny, double gnz,
				    double gpx, double gpy, double gpz ) {
    partition_periodic_box( grid, xh-xl, yh-yl, zh-zl,
                            (int)gnx, (int)gny, (int)gnz,
                            (int)gpx, (int)gpy, (int)gpz );
    grid->x0 += xl;
    grid->y0 += yl;
    grid->z0 += zl;
    IUO_allocate_fields();
  }

  inline void define_absorbing_grid( double xl,  double yl,  double zl,
                                     double xh,  double yh,  double zh,
                                     double gnx, double gny, double gnz,
                                     double gpx, double gpy, double gpz,
                                     int pbc ) {
    partition_absorbing_box( grid, xh-xl, yh-yl, zh-zl,
                             (int)gnx, (int)gny, (int)gnz,
                             (int)gpx, (int)gpy, (int)gpz,
                             pbc );
    grid->x0 += xl;
    grid->y0 += yl;
    grid->z0 += zl;
    IUO_allocate_fields();
    IUO_allocate_fields();
  }

  inline void define_reflecting_grid( double xl,  double yl,  double zl,
                                      double xh,  double yh,  double zh,
                                      double gnx, double gny, double gnz,
                                      double gpx, double gpy, double gpz ) {
    partition_metal_box( grid, xh-xl, yh-yl, zh-zl,
                         (int)gnx, (int)gny, (int)gnz,
                         (int)gpx, (int)gpy, (int)gpz );
    grid->x0 += xl;
    grid->y0 += yl;
    grid->z0 += zl;
    IUO_allocate_fields();
    IUO_allocate_fields();
  }

  // The below macros allow custom domains to be created

  // Creates a particle reflecting metal box in the local domain
  inline void size_domain( double lnx, double lny, double lnz ) {
    size_grid(grid,(int)lnx,(int)lny,(int)lnz);
    IUO_allocate_fields();
  }

  // Attaches a local domain boundary to another domain
  inline void join_domain( int boundary, double rank ) {
    join_grid( grid, boundary, (int)rank );
  }

  // Sets the field boundary condition of a local domain boundary
  inline void set_domain_field_bc( int boundary, int fbc ) {
    set_fbc( grid, boundary, fbc );
  }

  // Sets the particle boundary condition of a local domain boundary
  inline void set_domain_particle_bc( int boundary, int pbc ) {
    set_pbc( grid, boundary, pbc );
  }

  ///////////////////
  // Material helpers

  inline material_id define_material( const char *name,
				      double eps,
				      double mu = 1,
				      double sigma = 0 ) {
    return new_material(name,
			eps,  eps,  eps,
			mu,   mu,   mu,
			sigma,sigma,sigma,
			&material_list);
  }
  
  inline material_id define_material( const char *name,
				      double epsx,   double epsy,  double epsz,
				      double mux,    double muy,   double muz,
				      double sigmax,
                                      double sigmay,
				      double sigmaz ) {
    return new_material(name,
			epsx,  epsy,  epsz,
			mux,   muy,   muz,
			sigmax,sigmay,sigmaz,
			&material_list);
  }

  inline void finalize_materials(void) {
    material_coefficient = new_material_coefficients( grid, material_list );
  }

  inline material_id lookup_material( const char *name ) {
    material_t *m = find_material_name(name,material_list);
    return m==NULL ? invalid_material_id : m->id;
  }

  //////////////////
  // Species helpers

  inline species_id define_species( const char *name,
				    double q_m,
				    double max_local_np,
				    double sort_interval=25.,
				    double max_local_nm=-1. ) {
    // Compute a reasonble number of movers if user did not specify
    // Based on the twice the number of particles expected to hit the boundary
    // of a wpdt=0.2 / dx=lambda species in a 3x3x3 domain
    if( max_local_nm<=-1 ) {
      max_local_nm = 2*max_local_np/25;
      if( max_local_nm<16 ) max_local_nm = 16;
    }
    return new_species( name, (float)q_m, (int)max_local_np, (int)max_local_nm,
			(int)sort_interval, &species_list );
  }

  inline void finalize_species(void) {
    if( species_list==NULL ) return;
    species_lookup = new_species_lookup( species_list );
  }

  inline species_id lookup_species( const char *name ) {
    species_t *sp = find_species_name(name,species_list);
    return sp==NULL ? invalid_species_id : sp->id;
  }

  // BJA - Putting this here in first cut of collisionality addition in
  // order to facilitate implementation of collision models in input deck.  
  // This is a first cut at the problem; eventually, we will redesign this 
  // so that it is more user-transparent and supports user-defined collision 
  // operators. 

  inline species_t * find_species( const char *name ) {
     return find_species_name( name, species_list );  
  }


  ////////////////
  // Field helpers

  // Field helpers are provided by macros in deck_wrapper.cxx

  ///////////////////
  // Particle helpers

  error_code inject_particle( species_id id,
                              double x,  double y,  double z,
                              double ux, double uy, double uz,
                              double q,  double age = 0 );

  // It is more efficient to use species id for repeated particle
  // injection as it does not require the string compare costs for
  // species name lookup
  inline error_code inject_particle( const char *name,
                                     double x,  double y,  double z,
                                     double ux, double uy, double uz,
                                     double q,  double age = 0 ) {
    species_t *sp = find_species_name(name,species_list);
    if( sp==NULL ) return ERROR_CODE("Invalid species name");
    return inject_particle( sp->id, x, y, z, ux, uy, uz, q, age );
  }

  //////////////////////////////////
  // Random number generator helpers

  inline void seed_rand( double seed ) {
    mt_srand(rng,(int)seed);
  }

  // Uniform random number on (low,high) (open interval)
  inline double uniform_rand( double low, double high ) {
    return low + (high-low)*mt_drand(rng);
  }

  // Maxwellian random number with standard deviation dev
  inline double maxwellian_rand( double dev ) {
    return dev*mt_normal_drand(rng);
  }

  ////////////////////////
  // Miscellaneous helpers

  inline double nproc(void) {
    return (double)mp_nproc(grid->mp);
  }

  inline double time00( void ) {
    return mp_time00(grid->mp);
  }

  inline double elapsed( void ) {
    return mp_elapsed(grid->mp);
  }

  inline void abort( double code ) {
    mp_abort((int)code,grid->mp);
  }

  // Truncate "a" to the nearest integer multiple of "b"
  inline double trunc_granular( double a, double b ) {
    return b*int(a/b);
  }

  // Compute the remainder of a/b
  inline double remainder( double a, double b ) {
    return drem(a,b);
  }

  // Compute the Courant length on a regular mesh
  inline double courant_length( double lx, double ly, double lz,
				double nx, double ny, double nz ) {
    double w0, w1 = 0;
    if( nx>1 ) w0 = nx/lx, w1 += w0*w0;
    if( ny>1 ) w0 = ny/ly, w1 += w0*w0;
    if( nz>1 ) w0 = nz/lz, w1 += w0*w0;
    return sqrt(1/w1);
  }

  ////////////////////////////////////////////////////////////
  // User input deck provided functions (see deck_wrapper.cxx)

  void user_initialization( int argc, char **argv );
  void user_particle_injection(void);
  void user_current_injection(void);
  void user_field_injection(void);
  void user_diagnostics(void);
};

#endif // _vpic_hxx_
