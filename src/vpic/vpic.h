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

#ifndef vpic_h
#define vpic_h

#include <vector>
#include <cmath>
#include <functional>
#include <cassert>
#include <algorithm>

#include "../boundary/boundary.h"
#include "../collision/collision.h"
#include "../emitter/emitter.h"
// FIXME: INCLUDES ONCE ALL IS CLEANED UP
#include "../util/io/FileIO.h"
#include "../util/bitfield.h"
#include "../util/checksum.h"
#include "../util/system.h"
#include "dumpmacros.h"

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

// TODO: should this be an enum?
namespace dump_type {
  const int grid_dump = 0;
  const int field_dump = 1;
  const int hydro_dump = 2;
  const int particle_dump = 3;
  const int restart_dump = 4;
  const int history_dump = 5;
} // namespace


class vpic_simulation {
public:
  vpic_simulation();
  ~vpic_simulation();
  void initialize( int argc, char **argv );
  void modify( const char *fname );
  int advance( void );
  void finalize( void );

  #ifdef VPIC_GLOBAL_PARTICLE_ID
  // TODO: move these somewhere more sensible
  int predicate_count(species_t* sp, std::function <bool (int)> f);
  int predicate_count(species_t* sp, std::function <bool (particle_t)> f);

  // TOOD: those specialized in together should probably be wrapped in a class
  void predicate_copy(species_t* sp_from, species_t* sp_to, std::function <bool (int)> f);
  void predicate_copy(species_t* sp_from, species_t* sp_to, std::function <bool (particle_t)> f);
  #endif

protected:

  // Directly initialized by user

  int verbose;              // Should system be verbose
  int num_step;             // Number of steps to take
  int num_comm_round;       // Num comm round
  int status_interval;      // How often to print status messages
  int clean_div_e_interval; // How often to clean div e
  int num_div_e_round;      // How many clean div e rounds per div e interval
  int clean_div_b_interval; // How often to clean div b
  int num_div_b_round;      // How many clean div b rounds per div b interval
  int sync_shared_interval; // How often to synchronize shared faces

  // FIXME: THESE INTERVALS SHOULDN'T BE PART OF vpic_simulation
  // THE BIG LIST FOLLOWING IT SHOULD BE CLEANED UP TOO

  double quota;
  int checkpt_interval;
  int hydro_interval;
  int field_interval;
  int particle_interval;

  size_t nxout, nyout, nzout;
  size_t px, py, pz;
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

  // Helper initialized by user

  /* There are enough synchronous and local random number generators
     to permit the host thread plus all the pipeline threads for one
     dispatcher to simultaneously produce both synchronous and local
     random numbers.  Keeping the synchronous generators in sync is
     the generator users responsibility. */

  rng_pool_t           * entropy;            // Local entropy pool
  rng_pool_t           * sync_entropy;       // Synchronous entropy pool
  grid_t               * grid;               // define_*_grid et al
  material_t           * material_list;      // define_material
  field_array_t        * field_array;        // define_field_array
  interpolator_array_t * interpolator_array; // define_interpolator_array
  accumulator_array_t  * accumulator_array;  // define_accumulator_array
  hydro_array_t        * hydro_array;        // define_hydro_array
  species_t            * species_list;       // define_species /
                                             // species helpers
  particle_bc_t        * particle_bc_list;   // define_particle_bc /
                                             // boundary helpers
  emitter_t            * emitter_list;       // define_emitter /
                                             // emitter helpers
  collision_op_t       * collision_op_list;  // collision helpers

  // User defined checkpt preserved variables
  // Note: user_global is aliased with user_global_t (see deck_wrapper.cxx)

  char user_global[USER_GLOBAL_SIZE];

  /*----------------------------------------------------------------------------
   * Diagnostics
   ---------------------------------------------------------------------------*/
  double poynting_flux(double e0);

  /*----------------------------------------------------------------------------
   * Check Sums
   ---------------------------------------------------------------------------*/
#if defined(ENABLE_OPENSSL)
  void output_checksum_fields();
  void checksum_fields(CheckSum & cs);
  void output_checksum_species(const char * species);
  void checksum_species(const char * species, CheckSum & cs);
#endif // ENABLE_OPENSSL

  void print_available_ram() {
    SystemRAM::print_available();
  } // print_available_ram

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


#ifdef VPIC_GLOBAL_PARTICLE_ID
// TODO: merge back down to one function
// TODO: template out the functor type
// TODO: find a way to specify if we want to predicate on particle array, or
// particle index
  template<typename Predicate>
  void dump_particles_predicate(
      const char *sp_name,
      const char *fbase,
      int ftag,
      //const std::function <bool (int)>& f = nullptr
      //const std::function <bool (particle_t)>& f = nullptr
      const Predicate& f
  )
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

      WRITE_HEADER_V0( dump_type::particle_dump, sp->id, sp->q/sp->m, fileIO );

      int count_true = predicate_count(sp, f);
      std::cout << "copying " << count_true << " of " << sp->np << std::endl;

      dim[0] = count_true;
      WRITE_ARRAY_HEADER( p_buf, 1, dim, fileIO );

      // Copy a PBUF_SIZE hunk of the particle list into the particle
      // buffer, timecenter it and write it out. This is done this way to
      // guarantee the particle list unchanged while not requiring too
      // much memory.

      // FIXME: WITH A PIPELINED CENTER_P, PBUF NOMINALLY SHOULD BE QUITE
      // LARGE.

      // Make a second species array, and space to hold the IDs too
      particle_t* ALIGNED(128) _p;
      size_t*     ALIGNED(128) _p_id;
      MALLOC_ALIGNED( _p,    count_true, 128 );
      MALLOC_ALIGNED( _p_id, count_true, 128 );

      species_t _sp = *sp; // FIXME: Is this copy careful/safe enough?
      _sp.p = _p;
      _sp.np = count_true;
      _sp.p_id = _p_id;

      // TODO: why do we need to update max_np?
      _sp.max_np = PBUF_SIZE;

      // Copy the right particles over, that meet the predicate
      predicate_copy(sp, &_sp, f);

      // Use that instead below
      for( buf_start=0; buf_start<count_true; buf_start += PBUF_SIZE ) {
          _sp.np = count_true-buf_start;

          if( _sp.np > PBUF_SIZE ) {
              _sp.np = PBUF_SIZE;
          }

          //COPY( _sp.p, &sp_p[buf_start], _sp.np );

          center_p( &_sp, interpolator_array );

          fileIO.write( _p, _sp.np );
      }

      // append ID array at the end of the file
      if(sp->has_ids) {
          dim[0] = count_true;
          WRITE_ARRAY_HEADER( sp->p_id, 1, dim, fileIO );
          // Maybe do this write in batches of PBUF_SIZE as well
          fileIO.write(_sp.p_id, count_true);
      }

      // Free the particle array and ID array (sp done by scope)
      FREE_ALIGNED( _p );
      FREE_ALIGNED( _p_id );

      if( fileIO.close() ) ERROR(("File close failed on dump particles!!!"));
  }
#endif



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

  ///////////////////
  // Useful accessors

  inline void
  barrier() { mp_barrier(); }

  inline double
  time() {
    return grid->t0 + (double)grid->dt*(double)grid->step;
  }

  inline int64_t &
  step() {
   return grid->step;
  }

  inline field_t &
  field( const int v ) {
    return field_array->f[ v ];
  }

  inline int
  voxel( const int ix, const int iy, const int iz ) {
    return ix + grid->sy*iy + grid->sz*iz;
  }

  inline field_t &
  field( const int ix, const int iy, const int iz ) {
    return field_array->f[ voxel(ix,iy,iz) ];
  }

  inline interpolator_t &
  interpolator( const int v ) {
    return interpolator_array->i[ v ];
  }

  inline interpolator_t &
  interpolator( const int ix, const int iy, const int iz ) {
    return interpolator_array->i[ voxel(ix,iy,iz) ];
  }

  inline hydro_t &
  hydro( const int v ) {
    return hydro_array->h[ v ];
  }

  inline hydro_t &
  hydro( const int ix, const int iy, const int iz ) {
    return hydro_array->h[ voxel(ix,iy,iz) ];
  }

  inline rng_t *
  rng( const int n ) {
    return entropy->rng[n];
  }

  inline rng_t *
  sync_rng( const int n ) {
    return sync_entropy->rng[n];
  }

  ///////////////
  // Grid helpers

  inline void
  define_units( float cvac,
                float eps0 ) {
    grid->cvac = cvac;
    grid->eps0 = eps0;
  }

  inline void
  define_timestep( float dt, double t0 = 0, int64_t step = 0 ) {
    grid->t0   = t0;
    grid->dt   = (float)dt;
    grid->step = step;
  }

  // The below functions automatically create partition simple grids with
  // simple boundary conditions on the edges.

  inline void
  define_periodic_grid( double xl,  double yl,  double zl,
                        double xh,  double yh,  double zh,
                        double gnx, double gny, double gnz,
                        double gpx, double gpy, double gpz ) {
	px = size_t(gpx); py = size_t(gpy); pz = size_t(gpz);
    partition_periodic_box( grid, xl, yl, zl, xh, yh, zh,
                            (int)gnx, (int)gny, (int)gnz,
                            (int)gpx, (int)gpy, (int)gpz );
  }

  inline void
  define_absorbing_grid( double xl,  double yl,  double zl,
                         double xh,  double yh,  double zh,
                         double gnx, double gny, double gnz,
                         double gpx, double gpy, double gpz, int pbc ) {
	px = size_t(gpx); py = size_t(gpy); pz = size_t(gpz);
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
	px = size_t(gpx); py = size_t(gpy); pz = size_t(gpz);
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

  inline material_t *
  define_material( const char * name,
                   double eps,
                   double mu = 1,
                   double sigma = 0,
                   double zeta = 0 ) {
    return append_material( material( name,
                                      eps,   eps,   eps,
                                      mu,    mu,    mu,
                                      sigma, sigma, sigma,
                                      zeta,  zeta,  zeta ), &material_list );
  }

  inline material_t *
  define_material( const char * name,
                   double epsx,        double epsy,       double epsz,
                   double mux,         double muy,        double muz,
                   double sigmax,      double sigmay,     double sigmaz,
		   double zetax = 0 ,  double zetay = 0,  double zetaz = 0 ) {
    return append_material( material( name,
                                      epsx,   epsy,   epsz,
                                      mux,    muy,    muz,
                                      sigmax, sigmay, sigmaz,
                                      zetax,  zetay,  zetaz ), &material_list );
  }

  inline material_t *
  lookup_material( const char * name ) {
    return find_material_name( name, material_list );
  }

  inline material_t *
  lookup_material( material_id id ) {
    return find_material_id( id, material_list );
  }

  //////////////////////
  // Field array helpers

  // If fa is provided, define_field_advance will use it (and take ownership
  // of the it).  Otherwise the standard field array will be used with the
  // optionally provided radition damping parameter.

  inline void
  define_field_array( field_array_t * fa = NULL, double damp = 0 ) {
    int nx1 = grid->nx + 1, ny1 = grid->ny+1, nz1 = grid->nz+1;

    if( grid->nx<1 || grid->ny<1 || grid->nz<1 )
      ERROR(( "Define your grid before defining the field array" ));
    if( !material_list )
      ERROR(( "Define your materials before defining the field array" ));

    field_array        = fa ? fa :
                         new_standard_field_array( grid, material_list, damp );
    interpolator_array = new_interpolator_array( grid );
    accumulator_array  = new_accumulator_array( grid );
    hydro_array        = new_hydro_array( grid );

    // Pre-size communications buffers. This is done to get most memory
    // allocation over with before the simulation starts running

    mp_size_recv_buffer(grid->mp,BOUNDARY(-1, 0, 0),ny1*nz1*sizeof(hydro_t));
    mp_size_recv_buffer(grid->mp,BOUNDARY( 1, 0, 0),ny1*nz1*sizeof(hydro_t));
    mp_size_recv_buffer(grid->mp,BOUNDARY( 0,-1, 0),nz1*nx1*sizeof(hydro_t));
    mp_size_recv_buffer(grid->mp,BOUNDARY( 0, 1, 0),nz1*nx1*sizeof(hydro_t));
    mp_size_recv_buffer(grid->mp,BOUNDARY( 0, 0,-1),nx1*ny1*sizeof(hydro_t));
    mp_size_recv_buffer(grid->mp,BOUNDARY( 0, 0, 1),nx1*ny1*sizeof(hydro_t));

    mp_size_send_buffer(grid->mp,BOUNDARY(-1, 0, 0),ny1*nz1*sizeof(hydro_t));
    mp_size_send_buffer(grid->mp,BOUNDARY( 1, 0, 0),ny1*nz1*sizeof(hydro_t));
    mp_size_send_buffer(grid->mp,BOUNDARY( 0,-1, 0),nz1*nx1*sizeof(hydro_t));
    mp_size_send_buffer(grid->mp,BOUNDARY( 0, 1, 0),nz1*nx1*sizeof(hydro_t));
    mp_size_send_buffer(grid->mp,BOUNDARY( 0, 0,-1),nx1*ny1*sizeof(hydro_t));
    mp_size_send_buffer(grid->mp,BOUNDARY( 0, 0, 1),nx1*ny1*sizeof(hydro_t));
  }

  // Other field helpers are provided by macros in deck_wrapper.cxx

  //////////////////
  // Species helpers

  // FIXME: SILLY PROMOTIONS
  inline species_t *
  define_species( const char *name,
                  double q,
                  double m,
                  double max_local_np,
                  double max_local_nm,
                  double sort_interval,
                  double sort_out_of_place ) {
    // Compute a reasonble number of movers if user did not specify
    // Based on the twice the number of particles expected to hit the boundary
    // of a wpdt=0.2 / dx=lambda species in a 3x3x3 domain
    if( max_local_nm<0 ) {
      max_local_nm = 2*max_local_np/25;
      if( max_local_nm<16*(MAX_PIPELINE+1) )
        max_local_nm = 16*(MAX_PIPELINE+1);
    }
    return append_species( species( name, (float)q, (float)m,
                                    (size_t)max_local_np, (size_t)max_local_nm,
                                    (int)sort_interval, (int)sort_out_of_place,
                                    grid ), &species_list );
  }

  inline species_t *
  find_species( const char *name ) {
     return find_species_name( name, species_list );
  }

  inline species_t *
  find_species( int32_t id ) {
     return find_species_id( id, species_list );
  }

  // enum class Tracertype { Copy, Move };
  #ifdef VPIC_GLOBAL_PARTICLE_ID
  // REVIEW: make the name a std::string? or match define_species?
  inline species_t * make_tracers_by_percentage(species_t* parentspecies, const float percentage, const Tracertype copyormove, const char* tracername) {
    if((percentage <= 0.) || (percentage > 100.)) {
      ERROR(( "%f is a bad percentage to select tracers", percentage ));
    }
    // Implemented in species_advance.cc
    species_t* tracerspecies =  tracerspecies_by_skip(parentspecies, 100./percentage, copyormove, std::string(tracername), species_list, grid);
    return append_species(tracerspecies, &species_list);
  }
  inline species_t * make_tracers_by_nth(species_t* parentspecies, const float nth, const Tracertype copyormove, const char* tracername) {
    if(nth < 1.) {
      ERROR(( "%f is a bad stride to select every nth particle as tracers", nth ));
    }
    species_t* tracerspecies =  tracerspecies_by_skip(parentspecies, nth, copyormove, std::string(tracername), species_list, grid);
    return append_species(tracerspecies, &species_list);
  }
  inline species_t * make_n_tracers(species_t* parentspecies, const float n, const Tracertype copyormove, const char* tracername) {
    if(!parentspecies) ERROR(( "Invalid parent species" ));
    if((n < 1.) || (n > parentspecies->np)) {
      ERROR(( "%f is a bad number of tracers", n ));
    }
    species_t* tracerspecies =  tracerspecies_by_skip(parentspecies, parentspecies->np/n, copyormove, std::string(tracername), species_list, grid);
    return append_species(tracerspecies, &species_list);
  }

  // versions without user supplied name
  inline species_t * make_tracers_by_percentage(species_t* parentspecies, const float percentage, const Tracertype copyormove) {
    if(!parentspecies) ERROR(( "Invalid parent species" ));
    std::string name = make_tracer_name_unique(std::string(parentspecies->name) + std::string("_tracer"), species_list);
    return make_tracers_by_percentage(parentspecies, percentage, copyormove, name.c_str());
  }
  inline species_t * make_tracers_by_nth(species_t* parentspecies, const float nth, const Tracertype copyormove) {
    if(!parentspecies) ERROR(( "Invalid parent species" ));
    std::string name = make_tracer_name_unique(std::string(parentspecies->name) + std::string("_tracer"), species_list);
    return make_tracers_by_nth(parentspecies, nth, copyormove, name.c_str());
  }
  inline species_t * make_n_tracers(species_t* parentspecies, const float n, const Tracertype copyormove) {
    if(!parentspecies) ERROR(( "Invalid parent species" ));
    std::string name = make_tracer_name_unique(std::string(parentspecies->name) + std::string("_tracer"), species_list);
    return make_n_tracers(parentspecies, n, copyormove, name.c_str());
  }
  #else
  make_tracers_by_percentage(species_t* parentspecies, const float percentage, const Tracertype copyormove, const char* tacername) {
    ERROR(( "If you want to use make_tracers_by_percentage you need to compile with GLOBAL_PARTICLE_ID" ));
  }
  make_tracers_by_nth(species_t* parentspecies, const float nth, const Tracertype copyormove, const char* tacername) {
    ERROR(( "If you want to use make_tracers_by_nth you need to compile with GLOBAL_PARTICLE_ID" ));
  }
  make_n_tracers(species_t* parentspecies, const float n, const Tracertype copyormove, const char* tacername) {
    ERROR(( "If you want to use make_n_tracers you need to compile with GLOBAL_PARTICLE_ID" ));
  }
  make_tracers_by_percentage(species_t* parentspecies, const float percentage, const Tracertype copyormove) {
    ERROR(( "If you want to use make_tracers_by_percentage you need to compile with GLOBAL_PARTICLE_ID" ));
  }
  make_tracers_by_nth(species_t* parentspecies, const float nth, const Tracertype copyormove) {
    ERROR(( "If you want to use make_tracers_by_nth you need to compile with GLOBAL_PARTICLE_ID" ));
  }
  make_n_tracers(species_t* parentspecies, const float n, const Tracertype copyormove) {
    ERROR(( "If you want to use make_n_tracers you need to compile with GLOBAL_PARTICLE_ID" ));
  }
  #endif

  // Note: Don't use injection with aging during initialization
  // Defaults in the declaration below enable backwards compatibility.

  void
  inject_particle( species_t * sp,
                   double x,  double y,  double z,
                   double ux, double uy, double uz,
                   double w,  double age = 0, int update_rhob = 1 );

  // Inject particle raw is for power users!
  // No nannyism _at_ _all_:
  // - Availability of free stoarge is _not_ checked.
  // - Particle displacements and voxel index are _not_ for validity.
  // - The rhob field is _not_ updated.
  // - Injection with displacment may use up movers (i.e. don't use
  //   injection with displacement during initialization).
  // This injection is _ultra_ _fast_.

  inline void
  inject_particle_raw( species_t * RESTRICT sp,
                       float dx, float dy, float dz, int32_t i,
                       float ux, float uy, float uz, float w )
  {
    particle_t * RESTRICT p = sp->p + sp->np;
    #ifdef VPIC_GLOBAL_PARTICLE_ID
    if(sp->has_ids) {
      size_t * RESTRICT p_id = sp->p_id + sp->np;
      *p_id = sp->generate_particle_id( sp->np, sp->max_np );
    }
    #endif
    p->dx = dx; p->dy = dy; p->dz = dz; p->i = i;
    p->ux = ux; p->uy = uy; p->uz = uz; p->w = w;
    sp->np++;
  }

  // This variant does a raw inject and moves the particles

  inline void
  inject_particle_raw( species_t * RESTRICT sp,
                       float dx, float dy, float dz, int32_t i,
                       float ux, float uy, float uz, float w,
                       float dispx, float dispy, float dispz,
                       int update_rhob )
  {
    particle_t       * RESTRICT p  = sp->p  + sp->np;
    particle_mover_t * RESTRICT pm = sp->pm + sp->nm;
    #ifdef VPIC_GLOBAL_PARTICLE_ID
    if(sp->has_ids) {
      size_t           * RESTRICT p_id = sp->p_id + sp->np;
      *p_id = sp->generate_particle_id( sp->np, sp->max_np );
    }
    #endif
    p->dx = dx; p->dy = dy; p->dz = dz; p->i = i;
    p->ux = ux; p->uy = uy; p->uz = uz; p->w = w;
    pm->dispx = dispx; pm->dispy = dispy; pm->dispz = dispz; pm->i = sp->np-1;
    if( update_rhob ) accumulate_rhob( field_array->f, p, grid, -sp->q );
    sp->nm += move_p( sp->p, pm, accumulator_array->a, grid, sp->q );
    sp->np++;
  }

  //////////////////////////////////
  // Random number generator helpers

  // seed_rand seed the all the random number generators.  The seed
  // used for the individual generators is based off the user provided
  // seed such each local generator in each process (rng[0:r-1]) gets
  // a unique seed.  Each synchronous generator (sync_rng[0:r-1]) gets a
  // unique seed that does not overlap with the local generators
  // (common across each process).  Lastly, all these seeds are such
  // that, no individual generator seeds are reused across different
  // user seeds.
  // FIXME: MTRAND DESPERATELY NEEDS A LARGER SEED SPACE!

  inline void seed_entropy( int base ) {
    seed_rng_pool( entropy,      base, 0 );
    seed_rng_pool( sync_entropy, base, 1 );
  }

  // Uniform random number on (low,high) (open interval)
  // FIXME: IS THE INTERVAL STILL OPEN IN FINITE PRECISION
  //        AND IS THE OPEN INTERVAL REALLY WHAT USERS WANT??
  inline double uniform( rng_t * rng, double low, double high ) {
    double dx = drand( rng );
    return low*(1-dx) + high*dx;
  }

  // Normal random number with mean mu and standard deviation sigma
  inline double normal( rng_t * rng, double mu, double sigma ) {
    return mu + sigma*drandn( rng );
  }

  /////////////////////////////////
  // Emitter and particle bc helpers

  // Note that append_emitter is hacked to silently returne if the
  // emitter is already in the list.  This allows things like:
  //
  // define_surface_emitter( my_emitter( ... ), rgn )
  // ... or ...
  // my_emit_t * e = my_emit( ... )
  // define_surface_emitter( e, rgn )
  // ... or ...
  // my_emit_t * e = define_emitter( my_emit( ... ) )
  // define_surface_emitter( e, rng )
  // ...
  // All to work.  (Nominally, would like define_surface_emitter
  // to evaluate to the value of e.  But, alas, the way
  // define_surface_emitter works and language limitations of
  // strict C++ prevent this.)

  inline emitter_t *
  define_emitter( emitter_t * e ) {
    return append_emitter( e, &emitter_list );
  }

  inline particle_bc_t *
  define_particle_bc( particle_bc_t * pbc ) {
    return append_particle_bc( pbc, &particle_bc_list );
  }

  inline collision_op_t *
  define_collision_op( collision_op_t * cop ) {
    return append_collision_op( cop, &collision_op_list );
  }

  ////////////////////////
  // Miscellaneous helpers

  inline void abort( double code ) {
    nanodelay(2000000000); mp_abort((((int)code)<<17)+1);
  }

  // Truncate "a" to the nearest integer multiple of "b"
  inline double trunc_granular( double a, double b ) { return b*int(a/b); }

  // Compute the remainder of a/b
  inline double remainder( double a, double b ) { return std::remainder(a,b); }
  // remainder(a,b);

  // Compute the Courant length on a regular mesh
  inline double courant_length( double lx, double ly, double lz,
				double nx, double ny, double nz ) {
    double w0, w1 = 0;
    if( nx>1 ) w0 = nx/lx, w1 += w0*w0;
    if( ny>1 ) w0 = ny/ly, w1 += w0*w0;
    if( nz>1 ) w0 = nz/lz, w1 += w0*w0;
    return sqrt(1/w1);
  }

  //////////////////////////////////////////////////////////
  // These friends are used by the checkpt / restore service

  friend void checkpt_vpic_simulation( const vpic_simulation * vpic );
  friend vpic_simulation * restore_vpic_simulation( void );
  friend void reanimate_vpic_simulation( vpic_simulation * vpic );

  ////////////////////////////////////////////////////////////
  // User input deck provided functions (see deck_wrapper.cxx)

  void user_initialization( int argc, char **argv );
  void user_particle_injection(void);
  void user_current_injection(void);
  void user_field_injection(void);
  void user_diagnostics(void);
  void user_particle_collisions(void);
};


#endif // vpic_h
