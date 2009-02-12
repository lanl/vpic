/*
 * Written by:
 *   Kevin J. Bowers, Ph.D.
 *   Plasma Physics Group (X-1)
 *   Applied Physics Division
 *   Los Alamos National Lab
 * March/April 2004 - Heavily revised and extended from earlier V4PIC versions
 *
 * snell - revised to add new dumps, 20080310
 *
 */
 
#ifndef _vpic_hxx_
#define _vpic_hxx_
 
// FIXME: INCLUDES ONCE ALL IS CLEANED UP
#include "../emitter/emitter.h"
#include "../boundary/boundary.h"
#include "../species_advance/standard/spa.h"
#include <FileIO.hxx>
#include <BitField.hxx>
#include <vector>
 
#include <stdio.h>
 
#ifndef USER_GLOBAL_SIZE
#define USER_GLOBAL_SIZE 16384
#endif
 
#ifndef NVARHISMX
#define NVARHISMX 250
#endif
//  #include "dumpvars.h"
 
typedef FileIO FILETYPE;

const uint32_t all			(0xffffffff);
const uint32_t electric		(1<<0 | 1<<1 | 1<<2);
const uint32_t div_e_err	(1<<3);
const uint32_t magnetic		(1<<4 | 1<<5 | 1<<6);
const uint32_t div_b_err	(1<<7);
const uint32_t tca			(1<<8 | 1<<9 | 1<<10);
const uint32_t rhob			(1<<11);
const uint32_t current		(1<<12 | 1<<13 | 1<<14);
const uint32_t rhof			(1<<15);
const uint32_t emat			(1<<16 | 1<<17 | 1<<18);
const uint32_t nmat			(1<<19);
const uint32_t fmat			(1<<20 | 1<<21 | 1<<22);
const uint32_t cmat			(1<<23);

const size_t total_field_variables(24);
const size_t total_field_groups(12); // this counts vectors, tensors etc...
// These bits will be tested to determine which variables to output
const size_t field_indeces[12] = { 0, 3, 4, 7, 8, 11, 12, 15, 16, 19, 20, 23 };

struct FieldInfo {
	char name[128];
	char degree[128];
	char elements[128];
	char type[128];
	size_t size;
}; // struct FieldInfo

const uint32_t current_density	(1<<0 | 1<<1 | 1<<2);
const uint32_t charge_density	(1<<3);
const uint32_t momentum_density	(1<<4 | 1<<5 | 1<<6);
const uint32_t ke_density		(1<<7);
const uint32_t stress_tensor	(1<<8 | 1<<9 | 1<<10 | 1<<11 | 1<<12 | 1<<13);
/* May want to use these instead
const uint32_t stress_diagonal 		(1<<8 | 1<<9 | 1<<10);
const uint32_t stress_offdiagonal	(1<<11 | 1<<12 | 1<<13);
*/

const size_t total_hydro_variables(14);
const size_t total_hydro_groups(5); // this counts vectors, tensors etc...
// These bits will be tested to determine which variables to output
const size_t hydro_indeces[5] = { 0, 3, 4, 7, 8 };

struct HydroInfo {
	char name[128];
	char degree[128];
	char elements[128];
	char type[128];
	size_t size;
}; // struct FieldInfo

/*----------------------------------------------------------------------------
 * DumpFormat Enumeration
----------------------------------------------------------------------------*/
enum DumpFormat {
  band = 0,
  band_interleave = 1
}; // enum DumpFormat

/*----------------------------------------------------------------------------
 * DumpParameters Struct
----------------------------------------------------------------------------*/
struct DumpParameters {

  void output_variables(uint32_t mask) {
    output_vars.set(mask);
  } // output_variables

  BitField output_vars;

  size_t stride_x;
  size_t stride_y;
  size_t stride_z;

  DumpFormat format;

  char name[128];
  char baseDir[128];
  char baseFileName[128];

}; // struct DumpParameters

class vpic_simulation {
public:
  vpic_simulation();
  ~vpic_simulation();
  void initialize( int argc, char **argv );
  void restart( const char *filebase );
  void modify_runparams( const char *fname );
  int advance( void );
  void finalize( void );
 
  // some helpers that need to be exposed to
  // remove direct MPI calls
 
  inline double rank(void) {
    return (double)mp_rank(grid->mp);
  }
 
  inline void barrier(void) {
    mp_barrier(grid->mp);
  }
 
  inline void * grid_mp(void) {
  	return (void *)grid->mp;
  }
 
private:
 
  // Directly initialized by user; saved in a restart dump
 
  int verbose;              // Should system be verbose
  int step;                 // Number of steps taken so far
  int num_step;             // Number of steps to take
  int num_comm_round;       // Num comm round
  int status_interval;      // How often to print status messages
  int clean_div_e_interval; // How often to clean div e
  int clean_div_b_interval; // How often to clean div b
  int sync_shared_interval; // How often to synchronize shared faces
 
  double quota;
  int restart_interval;
  int hydro_interval;
  int field_interval;
  int particle_interval;
 
  size_t nxout, nyout, nzout;
  float dxout, dyout, dzout;

  int ndfld;
  int ndhyd;
  int ndpar;
  int ndhis;
  int ndgrd;
  int head_option;
  int istride;
  int jstride;
  int kstride;
  int stride_option;
  int pstride;
  int nprobe;
  int ijkprobe[NVARHISMX][4];
  float xyzprobe[NVARHISMX][3];
  int block_dump;
  int stepdigit;
  int rankdigit;
  int ifenergies;
 
  // Helper initialized by user; restart-saved
 
  mt_rng_t *rng;               // Random helpers (seed defaults to rank)
  material_t *material_list;   // Material helpers
  grid_t *grid;                // Grid helpers
  species_t *species_list;     // Species helpers / particle helpers
  emitter_t *emitter_list;     // Particle emitters
  field_advance_t * field_advance;
 
  // FIXME: Backward compatible hacks
  field_t * ALIGNED(128) field;
 
  // Helper initialized by user; not restart-saved (can be derived from above)
 
  /**/
  interpolator_t * ALIGNED(128) interpolator;
  accumulator_t * ALIGNED(128) accumulator;
  hydro_t * ALIGNED(128) hydro;
 
  // Internal use only variables; restart saved
 
  double p_time; // Time spent pushing particles since last status update
  double s_time; // Time spent performance sorting particles
  double g_time; // Time spent processing guard list since last status update
  double f_time; // Time spent processing fields since last status update
  double u_time; // Time spent in user functions since last status update
 
  // User defined variables; restart saved
 
  char user_global[USER_GLOBAL_SIZE];
  // Note: user_global is aliased with user_global_t (see deck_wrapper.cxx)
 
  /*----------------------------------------------------------------------------
   * Check Sums
   ---------------------------------------------------------------------------*/
#if defined(ENABLE_OPENSSL)
  void output_checksum_fields();
  void output_checksum_species(const char * species);
#endif // ENABLE_OPENSSL

  ///////////////
  // Dump helpers
 
  int dump_mkdir(const char * dname);
  int dump_cwd(char * dname, size_t size);

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
 
  // convenience functions for simlog output
  void create_field_list(char * strlist, DumpParameters & dumpParams);
  void create_hydro_list(char * strlist, DumpParameters & dumpParams);

  void print_hashed_comment(FileIO & fileIO, const char * comment);
  void global_header(const char * base,
  	std::vector<DumpParameters *> dumpParams);

  void field_header(const char * fbase, DumpParameters & dumpParams);
  void hydro_header(const char * speciesname, const char * hbase,
    DumpParameters & dumpParams);

  void field_dump(DumpParameters & dumpParams);
  void hydro_dump(const char * speciesname, DumpParameters & dumpParams);

  ///////////////
  // Grid helpers
 
  // The below functions automatically create partition simple grids with
  // simple boundary conditions on the edges.
 
  // FIXME: THE TIMESTEP SHOULD HAVE HELPERS LIKE THIS
 
  inline void
  define_periodic_grid( double xl,  double yl,  double zl,
                        double xh,  double yh,  double zh,
                        double gnx, double gny, double gnz,
                        double gpx, double gpy, double gpz ) {
    partition_periodic_box( grid, xl, yl, zl, xh, yh, zh,
                            (int)gnx, (int)gny, (int)gnz,
                            (int)gpx, (int)gpy, (int)gpz );
  }
 
  inline void
  define_absorbing_grid( double xl,  double yl,  double zl,
                         double xh,  double yh,  double zh,
                         double gnx, double gny, double gnz,
                         double gpx, double gpy, double gpz, int pbc ) {
    partition_absorbing_box( grid, xl, yl, zl, xh, yh, zh,
                             (int)gnx, (int)gny, (int)gnz,
                             (int)gpx, (int)gpy, (int)gpz,
                             pbc );
  }
 
  inline void
  define_reflecting_grid( double xl,  double yl,  double zl,
                          double xh,  double yh,  double zh,
                          double gnx, double gny, double gnz,
                          double gpx, double gpy, double gpz ) {
    partition_metal_box( grid, xl, yl, zl, xh, yh, zh,
                         (int)gnx, (int)gny, (int)gnz,
                         (int)gpx, (int)gpy, (int)gpz );
  }
 
  // The below macros allow custom domains to be created
 
  // Creates a particle reflecting metal box in the local domain
  inline void
  size_domain( double lnx, double lny, double lnz ) {
    size_grid(grid,(int)lnx,(int)lny,(int)lnz);
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
 
  inline material_id
  define_material( const char *name,
                   double eps,
                   double mu = 1,
                   double sigma = 0,
                   double zeta = 0 ) {
    return new_material( name,
			 eps,   eps,   eps,
			 mu,    mu,    mu,
			 sigma, sigma, sigma,
                         zeta,  zeta,  zeta,
			 &material_list );
  }
 
  inline material_id
  define_material( const char *name,
                   double epsx,        double epsy,       double epsz,
                   double mux,         double muy,        double muz,
                   double sigmax,      double sigmay,     double sigmaz,
		   double zetax = 0 ,  double zetay = 0,  double zetaz = 0 ) {
    return new_material( name,
			 epsx,   epsy,   epsz,
			 mux,    muy,    muz,
			 sigmax, sigmay, sigmaz,
                         zetax,  zetay,  zetaz,
			 &material_list );
  }
 
  inline material_id
  lookup_material( const char * name ) {
    material_t *m = find_material_name(name,material_list);
    return m==NULL ? invalid_material_id : m->id;
  }
 
  //////////////////////
  // Setup field advance
 
  inline void
  finalize_field_advance( field_advance_methods_t * fam = standard_field_advance ) {
    int nx1 = grid->nx + 1, ny1 = grid->ny+1, nz1 = grid->nz+1;
 
    field_advance = new_field_advance( grid, material_list, fam );
    field        = field_advance->f; // FIXME: Temporary hack
 
    interpolator = new_interpolator(grid);
    accumulator  = new_accumulators(grid);
    hydro        = new_hydro(grid);
 
    // Pre-size communications buffers
    // This is done to get 99.9% of memory allocation over with before
    // the simulation starts running
    mp_size_recv_buffer(BOUNDARY(-1, 0, 0),ny1*nz1*sizeof(hydro[0]),grid->mp);
    mp_size_recv_buffer(BOUNDARY( 1, 0, 0),ny1*nz1*sizeof(hydro[0]),grid->mp);
    mp_size_recv_buffer(BOUNDARY( 0,-1, 0),nz1*nx1*sizeof(hydro[0]),grid->mp);
    mp_size_recv_buffer(BOUNDARY( 0, 1, 0),nz1*nx1*sizeof(hydro[0]),grid->mp);
    mp_size_recv_buffer(BOUNDARY( 0, 0,-1),nx1*ny1*sizeof(hydro[0]),grid->mp);
    mp_size_recv_buffer(BOUNDARY( 0, 0, 1),nx1*ny1*sizeof(hydro[0]),grid->mp);
 
    mp_size_send_buffer(BOUNDARY(-1, 0, 0),ny1*nz1*sizeof(hydro[0]),grid->mp);
    mp_size_send_buffer(BOUNDARY( 1, 0, 0),ny1*nz1*sizeof(hydro[0]),grid->mp);
    mp_size_send_buffer(BOUNDARY( 0,-1, 0),nz1*nx1*sizeof(hydro[0]),grid->mp);
    mp_size_send_buffer(BOUNDARY( 0, 1, 0),nz1*nx1*sizeof(hydro[0]),grid->mp);
    mp_size_send_buffer(BOUNDARY( 0, 0,-1),nx1*ny1*sizeof(hydro[0]),grid->mp);
    mp_size_send_buffer(BOUNDARY( 0, 0, 1),nx1*ny1*sizeof(hydro[0]),grid->mp);
  }
 
  //////////////////
  // Species helpers
 
  // FIXME: SILLY PROMOTIONS
  inline species_t *
  define_species( const char *name,
                  double q_m,
                  double max_local_np,
                  double max_local_nm,
                  double sort_interval,
                  double sort_out_of_place ) {
    // Compute a reasonble number of movers if user did not specify
    // Based on the twice the number of particles expected to hit the boundary
    // of a wpdt=0.2 / dx=lambda species in a 3x3x3 domain
    if( max_local_nm<=-1 ) {
      max_local_nm = 2*max_local_np/25;
      if( max_local_nm<16*(MAX_PIPELINE+1) )
        max_local_nm = 16*(MAX_PIPELINE+1);
    }
    return new_species( name, (float)q_m, (int)max_local_np, (int)max_local_nm,
			(int)sort_interval, (int)sort_out_of_place, &species_list );
  }
 
  // BJA - Putting this here in first cut of collisionality addition in
  // order to facilitate implementation of collision models in input deck.
  // This is a first cut at the problem; eventually, we will redesign this
  // so that it is more user-transparent and supports user-defined collision
  // operators.
 
  inline species_t *
  find_species( const char *name ) {
     return find_species_name( name, species_list );
  }
 
  ////////////////
  // Field helpers
 
  // Field helpers are provided by macros in deck_wrapper.cxx
 
  ///////////////////
  // Particle helpers
 
  // Note: Don't use injection with aging during initialization
 
  void
  inject_particle( species_t * sp,
                   double x,  double y,  double z,
                   double ux, double uy, double uz,
                   double q,  double age, int update_rhob );
 
  // Inject particle raw is for power users!
  // No nannyism _at_ _all_:
  // - Availability of free stoarge is _not_ checked.
  // - Particle displacements and voxel index are _not_ for validity.
  // - The rhob field is _not_ updated.
  // - Injection with displacment may use up movers (i.e. don't use
  //   injection with displacement during initialization).
  // This injection is _ultra_ _fast_.
 
  inline void
  inject_particle_raw( species_t * sp,
                       float dx, float dy, float dz, int32_t i,
                       float ux, float uy, float uz, float q ) {
    particle_t * p = sp->p + (sp->np++);
    p->dx = dx; p->dy = dy; p->dz = dz; p->i = i;
    p->ux = ux; p->uy = uy; p->uz = uz; p->q = q;
  }
 
  inline void
  inject_particle_raw( species_t * sp,
                       float dx, float dy, float dz, int32_t i,
                       float ux, float uy, float uz, float q,
                       float dispx, float dispy, float dispz,
                       int update_rhob ) {
    particle_t       * p  = sp->p  + (sp->np++);
    particle_mover_t * pm = sp->pm + sp->nm;
    p->dx = dx; p->dy = dy; p->dz = dz; p->i = i;
    p->ux = ux; p->uy = uy; p->uz = uz; p->q = q;
    pm->dispx = dispx; pm->dispy = dispy; pm->dispz = dispz; pm->i = sp->np-1;
    if( update_rhob ) { p->q=-q; accumulate_rhob( field, p, grid ); p->q=q; }
    sp->nm += move_p( sp->p, pm, accumulator, grid );
  }
 
 
  //////////////////////////////////
  // Random number generator helpers
 
  inline void seed_rand( double seed ) {
    seed_mt_rng( rng, (int)seed );
  }
 
  // Uniform random number on (low,high) (open interval)
  // FIXME: THINK CAREFULLY ABOUT THE INTERVAL FINITE PRECISION HERE!
  inline double uniform_rand( double low, double high ) {
    double dx = mt_drand(rng);
    return low*(1-dx) + high*dx;
  }
 
  // Maxwellian random number with standard deviation dev
  inline double maxwellian_rand( double dev ) {
    return dev*mt_drandn(rng);
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
    return drem(a,b);   // remainder(a,b);
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
  void user_particle_collisions(void);
};
 
#endif // _vpic_hxx_
